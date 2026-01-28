#!/usr/bin/env python

"""
Create Position per Sample Mask Matrix from BED Files

This script aggregates per-sample flagged position BED files into a single
mask matrix where rows are positions and columns are samples.

The output is a TSV with TRUE/FALSE values indicating whether a position
should be masked (depth set to 0) for a given sample.

Command-line Arguments
----------------------
bed-files : str
    Space-separated list of sample BED files (*.flagged-pos.bed)
output : str
    Output mask matrix file path (TSV, gzipped)

Usage
-----
create_mask_matrix.py \\
    --bed-files sample1.flagged-pos.bed sample2.flagged-pos.bed ... \\
    --output flagged_positions.mask.tsv.gz
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


def create_mask_matrix(bed_files: list) -> None:
    """
    Create a position per sample mask matrix from per-sample BED files.
    
    Parameters
    ----------
    bed_files : list
        List of paths to sample-specific BED files
    output_file : str
        Output path for the mask matrix (TSV, gzipped)
    """
    LOG.info(f"Processing {len(bed_files)} BED files...")
    # Collect all positions and their associated samples
    mask_data = []
    
    for bed_file in bed_files:
        # Extract sample name from filename (e.g., "sample1.flagged-pos.bed" -> "sample1")
        sample_name = Path(bed_file).stem.replace('.flagged-pos', '')
        
        try:
            # Read BED file (CHROM, START, END, FILTER)
            bed_df = pd.read_csv(bed_file, sep="\t", header=None,
                                names=["CHROM", "START", "END", "FILTER"])
            
            if bed_df.empty:
                LOG.info(f"No flagged positions for {sample_name}")
                continue
            
            # Expand each region to individual positions
            for _, row in bed_df.iterrows():
                # BED files are 1-based inclusive (non-standard format)
                # For SNVs: START == END, so we just take that position
                # For ranges: include both START and END
                for pos in range(row['START'], row['END'] + 1):
                    mask_data.append({
                        'CHROM': row['CHROM'],
                        'POS': pos,
                        'SAMPLE': sample_name,
                        'KEEP': 0,  # Indicate position should be masked
                    })
            
            LOG.info(f"Added {len(bed_df)} positions for {sample_name}")
            
        except Exception as e:
            LOG.warning(f"Could not process {bed_file}: {e}")
            continue
    
    # Handle case where no positions need masking
    if not mask_data:
        LOG.info("No positions to mask across all samples, creating empty matrix")
        empty_df = pd.DataFrame(columns=['CHROM', 'POS'])
        empty_df.to_csv("flagged_positions.mask.tsv.gz", sep="\t", index=False, compression='gzip')
        return
    
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
    
    LOG.info(f"Created mask matrix: {len(mask_matrix)} positions × {len(mask_matrix.columns)-2} samples")
    
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