#!/usr/bin/env python

"""
Create Position per Sample Mask Matrix from BED Files

This script aggregates per-sample flagged position BED files into a single
mask matrix where rows are positions and columns are samples. 

Each sample-specific BED file (*.flagged-pos.bed) contains positions that should
be masked only for that particular sample.

The output is a TSV with 1/0 values indicating whether a position should be kept (1) 
or masked (0) for each sample.

Command-line Arguments
----------------------
bed-files : str
    Space-separated list of sample-specific BED files

Usage
-----
create_mask_matrix.py \\
    --bed-files sample1.flagged-pos.bed sample2.flagged-pos.bed ... \\
"""

import logging
from pathlib import Path

import click
import pandas as pd

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s",
    level=logging.INFO,
    datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("create_mask_matrix")

def add_bed_positions(bed_df: pd.DataFrame, 
                      sample_name: str, 
                      masked_positions: set, 
                      mask_data: list) -> int:
    """
    Add BED positions to mask data for a specific sample.
    
    Parameters
    ----------
    bed_df : pd.DataFrame
        BED dataframe with CHROM, START, END, FILTER columns
    sample_name : str
        Sample name to apply masking to
    masked_positions : set
        Set tracking already-masked (CHROM, POS, SAMPLE) tuples
    mask_data : list
        List to append mask entries to
    
    Returns
    -------
    int
        Number of new entries added
    """
    entries_added = 0
    
    for row in bed_df.itertuples():
        for pos in range(row.START, row.END + 1):
            key = (row.CHROM, pos, sample_name)
            if key not in masked_positions:
                mask_data.append({
                    'CHROM': row.CHROM,
                    'POS': pos,
                    'SAMPLE': sample_name,
                    'KEEP': 0,
                })
                masked_positions.add(key)
                entries_added += 1
    
    return entries_added

def prepare_matrix(mask_data: list) -> pd.DataFrame:
    """
    Prepare the mask matrix dataframe from collected mask data.
    
    Parameters
    ----------
    mask_data : list
        List of dictionaries with CHROM, POS, SAMPLE, KEEP keys
    
    Returns
    -------
    pd.DataFrame
        DataFrame in long format ready for pivoting to matrix
    """
    # Create dataframe from collected data
    mask_df = pd.DataFrame(mask_data)
    
    # Pivot to matrix format: rows = (CHROM, POS), columns = samples
    mask_matrix = mask_df.pivot_table(
        index=['CHROM', 'POS'],
        columns='SAMPLE',
        values='KEEP',
        fill_value=1,    # Positions not in a sample's BED = we should keep
        aggfunc='first'  # In case of duplicates, take first
    )

    # Reset index to make CHROM and POS regular columns
    mask_matrix = mask_matrix.reset_index()
    
    # Sort by chromosome and position for readability
    mask_matrix = mask_matrix.sort_values(['CHROM', 'POS']).reset_index(drop=True)

    return mask_matrix


def create_mask_matrix(bed_files: list) -> None:
    """
    Create a position per sample mask matrix from per-sample BED files.
    
    Parameters
    ----------
    bed_files : list
        List of paths to sample-specific BED files
    """
    LOG.info(f"Processing {len(bed_files)} BED files...")
    
    # Filter out all_samples aggregated file if present
    sample_bed_files = [f for f in bed_files if 'all_samples' not in Path(f).stem]
    
    LOG.info(f"Found {len(sample_bed_files)} sample-specific BED files")
    
    # Collect all sample's names
    all_samples = set()
    for bed_file in sample_bed_files:
        sample_name = Path(bed_file).stem.replace('.flagged-pos', '')
        all_samples.add(sample_name)
    
    # Track already-masked positions using a set for lookups
    masked_positions = set()  # Set of (CHROM, POS, SAMPLE) tuples to avoid duplicates within sample
    mask_data = []
    
    # Process sample-specific BED files
    sample_count = 0
    for bed_file in sample_bed_files:
        sample_name = Path(bed_file).stem.replace('.flagged-pos', '')
        
        try:
            bed_df = pd.read_csv(bed_file, sep="\t", header=None,
                                names=["CHROM", "START", "END", "FILTER"])
            
            if bed_df.empty:
                LOG.info(f"No flagged positions for {sample_name}")
                continue
            
            # Add sample-specific positions
            sample_specific = add_bed_positions(bed_df, sample_name, masked_positions, mask_data)
            
            LOG.info(f"{sample_name}: {len(bed_df)} positions in BED, {sample_specific} unique entries added")
            sample_count += sample_specific
            
        except Exception as e:
            LOG.warning(f"Could not process {bed_file}: {e}")
            continue
    
    LOG.info(f"Total mask entries: {sample_count}")
    
    # Handle case where no positions need masking
    if not mask_data:
        LOG.info("No positions to mask across all samples, creating empty matrix")
        empty_df = pd.DataFrame(columns=['CHROM', 'POS'])
        empty_df.to_csv("flagged_positions.mask.tsv.gz", sep="\t", index=False, compression='gzip')
        return
    
    # Create dataframe from collected data
    mask_matrix = prepare_matrix(mask_data)
    
    # Save to compressed TSV
    mask_matrix.to_csv("flagged_positions.mask.tsv.gz", sep="\t", index=False, compression='gzip')
    LOG.info(f"Mask matrix saved to: flagged_positions.mask.tsv.gz")


@click.command()
@click.option('--bed-files', multiple=True, required=True, type=click.Path(exists=True),
              help='Sample-specific flagged position BED files')
def main(bed_files: tuple):
    """
    Create mask matrix from sample BED files.
    """
    try:
        create_mask_matrix(list(bed_files))
    except Exception as e:
        LOG.error(f"Error creating mask matrix: {e}")
        raise


if __name__ == '__main__':
    main()