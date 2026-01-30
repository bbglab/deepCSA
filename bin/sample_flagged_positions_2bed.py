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
from utils_filter import expand_filter_column, extract_flagged_regions_bed, load_filter_criteria

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s", 
    level=logging.INFO, 
    datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("extract_sample_masked_bed")


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
    
    # Load filter criteria (excluding cohort-level filters)
    FILTERS = load_filter_criteria(json_filters, json_filters_somatic, only_cohort_filters=False)

    # Extract sample-specific masked positions to BED file
    LOG.info(f"Extracting sample-specific filters: {FILTERS}")
    
    # Extract flagged regions (will create {sample_name}.flagged-pos.bed)
    extract_flagged_regions_bed(maf_df, sample_name, FILTERS)
    
    LOG.info(f"Sample-specific masked BED file created: {sample_name}.flagged-pos.bed")


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
