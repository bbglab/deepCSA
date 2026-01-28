#!/usr/bin/env python

"""
Extract Sample-Specific Masked BED - Per-Sample Filter Extraction Script

This script processes a per-sample MAF file and extracts sample-specific filters into a BED format.
Unlike filter_cohort.py which handles cohort-level filters, this focuses on sample-level artifacts.

Command-line Arguments
----------------------
maf-file : str
    Path to the per-sample MAF file (can be gzipped).
json-filters : str
    Path to json file with list of filter criteria to consider for masking.
sample-id : str
    Sample identifier for output naming.

Usage
-----
sample_flagged_positions_2bed.py \\
    --maf-file sample.filtered.tsv.gz \\
    --json-filters filters.txt \\
    --sample-id sample_name

"""
import json
import logging
from pathlib import Path

import click
import pandas as pd
from read_utils import custom_na_values
from utils_filter import expand_filter_column

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s", 
    level=logging.INFO, 
    datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("extract_sample_masked_bed")

def extract_flagged_regions_bed(maf_df: pd.DataFrame, name: str, FILTERS: list[str]) -> pd.DataFrame | None:
    """
    Returns a BED file with the regions discarded, including the list of filters applied to each mutation.
    Creates a properly formatted BED file with 0-based coordinates and half-open intervals.

    Parameters
    ----------
    maf_df : pd.DataFrame
        Input MAF dataframe with filter columns. POS column should contain 1-based coordinates.
    bed_format : bool
        If True, output coordinates are converted to 0-based with half-open intervals [start, end).

    Returns
    -------
    pd.DataFrame
        A BED dataframe with discarded mutations and filters applied to each region.
        Output coordinates are 0-based with half-open intervals [start, end).
    """
    # List of filter columns you want to check for
    filter_columns = [f"FILTER.{f}" for f in FILTERS if f in ','.join(list(maf_df.columns))]

    maf_df_filters = maf_df[maf_df[filter_columns].any(axis=1)]

    if maf_df_filters.empty:
        LOG.warning("No mutations were flagged based on the applied filters.")
        return 

    # Create BED-like dataframe with filter columns
    bed_df = maf_df_filters[["CHROM", "POS"] + filter_columns]

    # Transform to long format
    _bed_melt = (pd.melt(bed_df,
                    id_vars=["CHROM", "POS"],
                    value_vars=filter_columns,
                    var_name="FILTERS")
            .query("value == True")
            )

    LOG.info("Mutations flagged: %s", _bed_melt.shape[0])

    # Aggregate filters per position
    bed_annotated = (
                _bed_melt
                .drop_duplicates()
                .groupby(["CHROM","POS"])["FILTERS"]
                .agg(','.join)
                .reset_index()
                .rename(columns={"POS": "START"})
    )

    # This creates 1-based inclusive coordinates -> not BED format
    # This is because this bed will be used to mask positions in depth files which are 1-based
    bed_annotated["END"] = bed_annotated["START"]

    LOG.info("Unique regions flagged: %s", bed_annotated.shape[0])

    # Write BED file without header or index
    (bed_annotated[["CHROM", "START", "END", "FILTERS"]]
        .sort_values(by=["CHROM", "START"])
        .to_csv(f"{name}.flagged-pos.bed", sep="\t", header=False, index=False)
    )

def extract_flagged_positions(maf_file: str, json_filters: str, json_filters_somatic, sample_name: str) -> None:
    """
    Extract sample-specific masked positions from a per-sample MAF file.

    Parameters
    ----------
    maf_file : str
        Path to input MAF file
    json_filters : str
        Path to json file with list of filter criteria to consider for masking.
    sample_name : str
        Sample identifier for output naming.
    """
    LOG.info(f"Processing sample: {sample_name}")
    
    # Load MAF dataframe
    maf_df = pd.read_csv(maf_file, sep='\t', header=0, na_values=custom_na_values)
    LOG.info(f"Loaded {maf_df.shape[0]} mutations for sample {sample_name}")
    
    # Expand FILTER column to individual boolean columns
    maf_df = expand_filter_column(maf_df)
    
    # Load the filter criteria from the JSON file
    with open(json_filters, 'r') as file:
        filter_data = json.load(file)
    with open(json_filters_somatic, 'r') as file:
        filter_data_somatic = json.load(file)
    
    # Extract the FILTER list from the JSON structure
    filter_criteria = filter_data.get("FILTER", []) + filter_data_somatic.get("FILTER", [])
    LOG.info(f"Loaded filter criteria: {filter_criteria}")
    
    # Get all the filters to consider for this sample from the filter criteria and filter criteria somatic
    FILTERS = [f.replace("notcontains ", "") for f in filter_criteria if "notcontains" in f]

    # Extract sample-specific masked positions to BED file
    LOG.info(f"Extracting sample-specific filters: {FILTERS}")
    
    # Extract flagged regions (will create {sample_name}.flagged-pos.bed)
    extract_flagged_regions_bed(maf_df, sample_name, FILTERS)
    
    # Rename to the requested output name if different
    generated_bed = f"{sample_name}.flagged-pos.bed"
    if Path(generated_bed).exists():
        LOG.info(f"Sample-specific masked BED file created: {generated_bed}")
    else:
        LOG.warning(f"No masked positions found for sample {sample_name}")


@click.command()
@click.option('--maf-file', required=True, type=click.Path(exists=True), 
              help='Input per-sample MAF file (TSV, can be gzipped)')
@click.option('--json-filters', required=True, type=click.Path(), 
              help='Path to json file with list of filter criteria to consider for masking.')
@click.option('--json-filters-somatic', required=True, type=click.Path(), 
              help='Path to json file with list of filter criteria to consider for masking.')
@click.option('--sample-name', required=True, type=str, 
              help='Sample identifier')
def main(maf_file: str, json_filters: str, json_filters_somatic: str, sample_name: str):
    """
    CLI wrapper for extracting sample-specific masked positions.
    """
    try:
        extract_flagged_positions(maf_file, json_filters, json_filters_somatic, sample_name)
    except Exception as e:
        LOG.error(f"Error processing sample {sample_name}: {e}")
        raise


if __name__ == '__main__':
    main()
