#!/usr/bin/env python

import logging
import click
import pandas as pd
import json
from read_utils import custom_na_values
from utils_filter import expand_filter_column, extract_flagged_regions_bed, load_filter_criteria

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s",
    level=logging.INFO,
    datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("write_mafs")

@click.command()
@click.option('--maf-file', required=True, type=click.Path(exists=True), help='Input gzipped MAF file (TSV)')
@click.option('--groups-json', required=True, type=click.Path(exists=True), help='Optional JSON file with group/sample mapping')
@click.option('--filters', required=False, type=str, default='', help='Comma-separated list of filter criteria')
@click.option('--somatic-filters', required=False, type=str, default='', help='Comma-separated list of somatic filter criteria')
def main(maf_file, groups_json, filters: str, somatic_filters: str):
    maf_df = pd.read_csv(maf_file, compression='gzip', header=0, sep='\t', na_values=custom_na_values)
    maf_df["SAMPLE_ID"] = maf_df["SAMPLE_ID"].astype(str)

    with open(groups_json, 'r') as file:
        groups_info = json.load(file)

    for group_name, samples in groups_info.items():
        samples = [str(x) for x in samples]
        maf_df[maf_df["SAMPLE_ID"].isin(samples)].sort_values(by=["CHROM", "POS"]).to_csv(
            f"{group_name}.filtered.tsv.gz",
            sep="\t",
            header=True,
            index=False
        )
    
    # Expand FILTER column for BED extraction
    maf_df = expand_filter_column(maf_df)
    
    # Determine which cohort-level filters to extract based on configuration
    # Load filter criteria (only cohort-level filters)
    cohort_filters = load_filter_criteria(filters, somatic_filters, only_cohort_filters=True)
    
    # Extract flagged regions to BED file (cohort-wide, applies to all samples)
    LOG.info("Extracting cohort-wide flagged positions to BED file")
    extract_flagged_regions_bed(maf_df, "shared_cohort", cohort_filters)

if __name__ == '__main__':
    main()
