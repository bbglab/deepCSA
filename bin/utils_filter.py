#!/usr/bin/env python
import json
import logging
import pandas as pd
from pathlib import Path
"""
Utility functions for extracting filters from a MAF DataFrame.
"""

LOG = logging.getLogger(__name__)


def load_filter_criteria(filters: str, somatic_filters: str, only_cohort_filters: bool = False) -> list[str]:
    """
    Parse filter criteria from comma-separated strings.
    
    Parameters
    ----------
    filters : str
        Comma-separated list of filter criteria
    somatic_filters : str
        Comma-separated list of somatic filter criteria
    only_cohort_filters : bool
        If True, only include cohort-level filters from the returned list
    
    Returns
    -------
    list[str]
        List of filter names to apply
    """
    # Parse comma-separated strings into lists
    filter_list = [f.strip() for f in filters.split(',') if f.strip()]
    somatic_filter_list = [f.strip() for f in somatic_filters.split(',') if f.strip()]
    
    # Combine both lists
    all_filters = filter_list + somatic_filter_list
    
    # Define cohort-level filters
    COHORT_FILTERS = {
        'repetitive_variant', 'repetitive_mapping_variant',
        'cohort_n_rich', 'cohort_n_rich_uni',
        'cohort_n_rich_threshold', 'other_sample_SNP', 'gnomAD_SNP'
    }
    
    # Filter based on cohort_filters flag
    if only_cohort_filters:
        result = [f.replace("notcontains ", "") for f in all_filters if f.startswith("notcontains ") and f.replace("notcontains ", "") in COHORT_FILTERS]
    else:
        result = [f.replace("notcontains ", "") for f in all_filters if f.startswith("notcontains ") and f.replace("notcontains ", "") not in COHORT_FILTERS]
        
    LOG.info(f"Loaded {len(result)} filter criteria: {result}")
    return result

def expand_filter_column(maf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Expands the FILTER column by creating new columns for each unique filter.
    Each new column indicates if the corresponding filter is present (True/False).
    """
    # Split FILTER column once per row and convert to set for O(1) lookup
    filter_sets = maf_df["FILTER"].str.split(";").apply(lambda x: set(x) if x != [''] else set())
    
    # Get all unique filter values (excluding empty strings)
    all_filters = set(
        filter_val 
        for filter_val in maf_df["FILTER"].str.split(";").explode().unique() 
        if filter_val and filter_val != ''
    )
    
    # Ensure "not_covered" and "not_in_exons" exist
    required_filters = {"not_covered", "not_in_exons"}
    all_filters.update(required_filters)

    # Create boolean columns efficiently
    for filt in sorted(all_filters):
        maf_df[f"FILTER.{filt}"] = filter_sets.apply(lambda x: filt in x)

    return maf_df

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

    maf_df_filters = maf_df[maf_df[filter_columns].any(axis=1)] if filter_columns else pd.DataFrame()

    if maf_df_filters.empty:
        LOG.warning("No mutations were flagged based on the applied filters. Creating empty BED file.")
        # Create empty BED file to satisfy pipeline requirements
        Path(f"{name}.flagged-pos.bed").touch()
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

    # The idea is to filter depth files at these positions, so make END = START (1-based)
    bed_annotated["END"] = bed_annotated["START"]

    LOG.info("Unique regions flagged: %s", bed_annotated.shape[0])

    # Write BED file without header or index
    (bed_annotated[["CHROM", "START", "END", "FILTERS"]]
        .sort_values(by=["CHROM", "START"])
        .to_csv(f"{name}.flagged-pos.bed", sep="\t", header=False, index=False)
    )