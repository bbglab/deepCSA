#!/usr/bin/env python

"""
Write MAFs - Create MAF file per sample and group. Extract flagged regions to BED file.

This script processes a Mutation Annotation Format (MAF) file to separate it
into individual MAF files for group defined in a provided JSON file.
It also extracts regions flagged by specified filters into a BED file for downstream analysis.
We obtain the flagged regions BED for both individual samples and the entire cohort.

Command-line Arguments
----------------------
maf-file : str
    Path to the gzipped input MAF file.
groups-json : str
    Path to the JSON file containing group/sample mapping.
filters : str
    Comma-separated list of filter criteria to apply.
somatic-filters : str
    Comma-separated list of somatic filter criteria to apply.

Authors
-------
Author  : Ferriol Calvet (@FerriolCalvet)
Email   : ferriol.calvet@irbbarcelona.org

Contributors
------------
- Marta Huertas - @m-huertasp (marta.huertas@irbbarcelona.org)

Usage
-----
python write_mafs.py \\
    --maf-file input.maf.tsv.gz \\
    --groups-json groups.json \\
    --filters "filter1,filter2,filter3" \\
    --somatic-filters "somatic_filter1,somatic_filter2"
"""
import logging
import click
import pandas as pd
import json
from read_utils import custom_na_values
from utils_filter import extract_flagged_regions_bed, load_filter_criteria

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s",
    level=logging.INFO,
    datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("write_mafs")

# Constants
# Define flags that are applied for all samples equally.
COHORT_FLAGS = ["cohort_n_rich", "cohort_n_rich_uni", "repetitive_variant", "gnomAD_SNP", "nanoseq_snp", "nanoseq_noise", "not_covered", "not_in_exons"]

def create_group_mafs(maf_df: pd.DataFrame, group_name: str, samples: list[str]) -> None:
    """
    Create separate MAF files for each group defined in the groups_info.

    Parameters
    ----------
    maf_df : pd.DataFrame
        Input MAF DataFrame
    group_name : str
        Name of the group
    samples : list[str]
        List of sample IDs belonging to the group
    """
    # Create MAF for the group by filtering the samples
    group_maf = maf_df[maf_df["SAMPLE_ID"].isin(samples)].sort_values(by=["CHROM", "POS"])
    # Generate output file
    group_maf.to_csv(f"{group_name}.filtered.tsv.gz", sep="\t", header=True, index=False, compression='gzip')
    LOG.info(f"Written MAF for group '{group_name}' with {len(group_maf)} records.")

def create_sample_flagged_bed(maf_df: pd.DataFrame, samples: list[str], filter_criteria: list[str]) -> None:
    """
    Extract the flagged regions for each sample and write them to BED files.

    Parameters
    ----------
    maf_df : pd.DataFrame
        Input MAF DataFrame
    samples : list[str]
        List of sample IDs to process
    filter_criteria : list[str]
        List of filter criteria to apply when extracting flagged regions
    """
    for sample in samples:
        sample_maf = maf_df[maf_df["SAMPLE_ID"] == sample]
        LOG.info(f"Extracting flagged positions from {sample}.")
        extract_flagged_regions_bed(sample_maf, sample, filter_criteria)  


@click.command()
@click.option('--maf-file', required=True, type=click.Path(exists=True), help='Input gzipped MAF file (TSV)')
@click.option('--groups-json', required=True, type=click.Path(exists=True), help='Optional JSON file with group/sample mapping')
@click.option('--filters', required=False, type=str, default='', help='Comma-separated list of filter criteria')
@click.option('--somatic-filters', required=False, type=str, default='', help='Comma-separated list of somatic filter criteria')
def main(maf_file, groups_json, filters: str, somatic_filters: str):
    maf_df = pd.read_csv(maf_file, compression='gzip', header=0, sep='\t', na_values=custom_na_values)
    maf_df["SAMPLE_ID"] = maf_df["SAMPLE_ID"].astype(str)
    filter_criteria = load_filter_criteria(filters, somatic_filters)
    
    with open(groups_json, 'r') as file:
        groups_info = json.load(file)

    # Obtain the set of all samples to be analyzed
    all_samples = set(groups_info.get("all_samples"))
    # Obtain the set of samples for which we have information in the MAF file
    available_samples = set(maf_df["SAMPLE_ID"].unique())

    requested_n_available_samples = available_samples.intersection(all_samples)

    if len(all_samples) != len(requested_n_available_samples):
        missing_samples = all_samples - available_samples

        error = f"Some SAMPLE_IDs listed in the features table have no matching entries in the MAF file. Missing SAMPLE_IDs: {', '.join(missing_samples)}"
        LOG.error(error)
        raise ValueError(error)
        
    for group_name, samples in groups_info.items():
        samples = [str(x) for x in samples]
        # Extract group-specific MAFs
        create_group_mafs(maf_df, group_name, samples)

    # Extract flagged regions for each sample in the group if not already done
    create_sample_flagged_bed(maf_df, list(all_samples), filter_criteria)

    # Extract flagged regions to BED file (cohort-wide, applies to all samples)
    LOG.info("Extracting all flagged positions to BED file.")
    extract_flagged_regions_bed(maf_df, "all_samples", filter_criteria, "all-")

    # Extract flagged regions applied only for all the cohort, not sample-specific
    LOG.info("Extracting only cohort-wide flagged positions to BED file.")
    cohort_wide_filters = [f for f in COHORT_FLAGS if f in filter_criteria]
    extract_flagged_regions_bed(maf_df, "all_samples", cohort_wide_filters, "cohort-wide-")

if __name__ == '__main__':
    main()
