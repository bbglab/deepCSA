#!/usr/local/bin/python

import argparse
import numpy as np
import pandas as pd
from read_utils import custom_na_values
from utils import add_filter
import logging

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("filter_cohort")

# Globals
FILTERS = ["cohort_n_rich", "cohort_n_rich_uni", "other_sample_SNP", "repetitive_variant"]

def correct_vaf(maf):
    """
    Computes ``VAF_CORRECTED`` for the subset of variants satisfying 0 < VAF < 0.2.
    !!!!! THIS FUNCTION IS NOT USED IN THE FINAL PIPELINE !!!!!

    Parameters
    ----------
    maf : pandas.DataFrame
        Input MAF dataframe

    Returns
    -------
    pandas.DataFrame

    Notes
    -----
    The output is a copy of the input MAF dataframe with two new columns:

    - ``VAF_CORRECTED``: Corrected VAF for the subset of variants satisfying 0 < VAF < 0.2
    - ``IS_VAF_CORRECTED``: Boolean indicating whether the VAF has been corrected
    """

    # TODO revise the 0.2 VAF threshold to see if it can be kept across datasets
    df = maf[(0 < maf["VAF"]) & (maf["VAF"] < 0.2)][["SAMPLE_ID", "MUT_ID", "VAF", "DEPTH"]].sort_values("DEPTH").copy()
    N = df.shape[0]

    df["VAF_ROLLING_MEAN"] = df["VAF"].rolling(N // 25).mean()
    df["VAF_ROLLING_STD"] = df["VAF"].rolling(N // 25).std()

    stable_mean = df["VAF_ROLLING_MEAN"].values[-1]
    stable_std = df["VAF_ROLLING_STD"].values[-1]

    df["VAF_CORRECTED"] = df.apply(
        lambda r: (r["VAF"] - r["VAF_ROLLING_MEAN"]) * (stable_std / r["VAF_ROLLING_STD"]) + stable_mean, axis=1
    )
    df = maf.merge(df[["VAF_CORRECTED", "MUT_ID", "SAMPLE_ID"]], on=["MUT_ID", "SAMPLE_ID"], how="outer")
    df["IS_VAF_CORRECTED"] = ~df["VAF_CORRECTED"].isnull()
    df.loc[~df["IS_VAF_CORRECTED"], "VAF_CORRECTED"] = df[~df["IS_VAF_CORRECTED"]]["VAF"].values
    return df

# Functions
def filter_repetitive_variants(maf_df, repetitive_variant_treshold, somatic_vaf_boundary):
    """
    Filters out repetitive variants from the MAF dataframe. A variant is considered repetitive if it appears in at least
    ``repetitive_variant_treshold`` samples

    Parameters
    ----------
    maf_df : pandas.DataFrame
        MAF dataframe
    repetitive_variant_treshold : int
        Minimum number of samples a variant must appear in to be considered repetitive

    Returns
    -------
    pandas.DataFrame
        MAF dataframe with a new column 'repetitive_variant' that flags repetitive variants
    """

    # TODO revise these numbers, the repetitive_variant_treshold is the boundary at which we start considering a mutation as "repetitive"
    max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))

    if max_samples < repetitive_variant_treshold:
        LOG.warning("Not enough samples to identify potential repetitive variants!")

        return maf_df

    # work with already filtered df + somatic only to explore potential artifacts
    # Filter somatic variants with VAF <= somatic_vaf_boundary
    maf_df_somatic = maf_df.loc[maf_df["VAF"] <= somatic_vaf_boundary][["MUT_ID", "SAMPLE_ID"]].reset_index(drop=True)

    # Group by 'MUT_ID' and count occurrences
    maf_df_somatic_pivot = maf_df_somatic.groupby("MUT_ID").size().reset_index(name="count")

    # Store repetitive variants
    repetitive_variants = maf_df_somatic_pivot[maf_df_somatic_pivot["count"] >= repetitive_variant_treshold]["MUT_ID"]
    LOG.info(f"{len(repetitive_variants)} repetitive_variants")

    maf_df["repetitive_variant"] = maf_df["MUT_ID"].isin(repetitive_variants)

    maf_df["FILTER"] = maf_df[["FILTER", "repetitive_variant"]].apply(
        lambda x: add_filter(x["FILTER"], x["repetitive_variant"], "repetitive_variant"), axis=1
    )
    maf_df = maf_df.drop("repetitive_variant", axis=1)

    return maf_df


def filter_cohort_n_rich(maf_df, somatic_vaf_boundary):
    """
    Filters out cohort_n_rich variants from the MAF dataframe

    Parameters
    ----------
    maf_df : pandas.DataFrame
        MAF dataframe
    somatic_vaf_boundary : float
        VAF boundary to consider a variant as somatic
    """

    max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))

    if max_samples < 2:
        # TODO add basic logger instead of print
        LOG.warning("Not enough samples to identify cohort_n_rich mutations!")
        return maf_df

    # work with already filtered df + somatic only to explore potential artifacts
    # Filter somatic variants with VAF <= somatic_vaf_boundary
    maf_df_f_somatic = maf_df[maf_df["VAF"] <= somatic_vaf_boundary][["MUT_ID", "SAMPLE_ID", "FILTER"]].reset_index(
        drop=True
    )

    n_rich_vars_df = maf_df_f_somatic[maf_df_f_somatic["FILTER"].str.contains("n_rich")].groupby("MUT_ID").size()
    n_rich_vars = list(n_rich_vars_df[n_rich_vars_df > 1].index)

    maf_df["cohort_n_rich"] = maf_df["MUT_ID"].isin(n_rich_vars)
    LOG.info(f"{maf_df['cohort_n_rich'].sum()} muts flagged as cohort_n_rich") # output the number of mutations flagged as cohort_n_rich

    maf_df["FILTER"] = maf_df[["FILTER", "cohort_n_rich"]].apply(
        lambda x: add_filter(x["FILTER"], x["cohort_n_rich"], "cohort_n_rich"), axis=1
    )
    maf_df = maf_df.drop("cohort_n_rich", axis=1)

    # if the variant appeared flagged as n_rich in a single sample it is also filtered out from all other samples
    n_rich_vars_uni = list(n_rich_vars_df[n_rich_vars_df > 0].index)

    maf_df["cohort_n_rich_uni"] = maf_df["MUT_ID"].isin(n_rich_vars_uni)
    LOG.info(f"{maf_df['cohort_n_rich_uni'].sum()} muts flagged as cohort_n_rich_uni")

    maf_df["FILTER"] = maf_df[["FILTER", "cohort_n_rich_uni"]].apply(
        lambda x: add_filter(x["FILTER"], x["cohort_n_rich_uni"], "cohort_n_rich_uni"), axis=1
    )

    maf_df = maf_df.drop("cohort_n_rich_uni", axis=1)

    return maf_df


def filter_other_samples_snp(maf_df, somatic_vaf_boundary):
    """
    Filters out SNPs from other samples from the MAF dataframe

    Parameters
    ----------
    maf_df : pandas.DataFrame
        MAF dataframe
    somatic_vaf_boundary : float
        VAF boundary to consider a variant as somatic
    """

    # this is if we were to consider both unique and no-unique variants
    germline_vars_all_samples = maf_df.loc[maf_df["VAF"] > somatic_vaf_boundary, "MUT_ID"].unique()
    LOG.debug(f"{len(germline_vars_all_samples)} using all germline variants of all samples")

    # this is if we were to consider only unique germline
    unique_germline_vars = (
        maf_df.loc[maf_df["VAF"] > somatic_vaf_boundary][["MUT_ID", "SAMPLE_ID"]]
        .groupby("MUT_ID")
        .size()
        .sort_values(ascending=False)
        .to_frame("n_samples")
        .reset_index()
    )
    LOG.debug(f"{len(unique_germline_vars)} using only germline variants unique to a single sample")

    maf_df["other_sample_SNP"] = np.where(
        (maf_df["MUT_ID"].isin(germline_vars_all_samples)) & (maf_df["VAF"] <= somatic_vaf_boundary), True, False
    )

    LOG.info(f"{maf_df['other_sample_SNP'].sum()} other_sample_SNP")

    maf_df["FILTER"] = maf_df[["FILTER", "other_sample_SNP"]].apply(
        lambda x: add_filter(x["FILTER"], x["other_sample_SNP"], "other_sample_SNP"), axis=1
    )
    maf_df = maf_df.drop("other_sample_SNP", axis=1)

    return maf_df


def explode_filter_column(maf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Expands the FILTER column, creating new columns for each unique filter.
    """
    for filt in pd.unique(maf_df["FILTER"].str.split(";").explode()):
        maf_df[f"FILTER.{filt}"] = maf_df["FILTER"].apply(lambda x: filt in x.split(";"))

    return maf_df


def get_discarded_regions(
    maf_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Returns a BED file with the regions discarded, including the list of filters applied to each mutation.

    Parameters
    ----------
    maf_df : pd.DataFrame
        Input MAF dataframe with filter columns.

    Returns
    -------
    pd.DataFrame
        A BED dataframe with discarded mutations and filters applied to each region.
    """

    # List of filter columns you want to check for
    filter_columns = [f"FILTER.{f}" for f in FILTERS if f in ','.join(list(maf_df.columns))]

    LOG.info(f"Filters applied: {filter_columns}")
    
    maf_df_filters = maf_df[maf_df[filter_columns].any(axis=1)]
    
    if maf_df_filters.empty:
        LOG.warning("No mutations were discarded based on the applied filters.")
        return

    # Output only CHROM, START, and END columns to a BED file
    bed_df = maf_df_filters[["CHROM", "POS"] + filter_columns]

    _bed_melt = (pd.melt(bed_df, 
                        id_vars=["CHROM", "POS"], 
                        value_vars=filter_columns, 
                        var_name="FILTERS")
                .query("value == True")
                )
    LOG.info(f"Discarded mutations: {_bed_melt.shape[0]}")

    bed_annotated = (
                _bed_melt
                .drop_duplicates()
                .groupby(["CHROM","POS"])["FILTERS"]
                .apply(lambda x: ','.join(x))
                .reset_index()
                .rename(columns={"POS": "START"})
    )

    bed_annotated["END"] = bed_annotated["START"]

    LOG.info(f"Discarded unique regions: {bed_annotated.shape[0]}")

    # Write the BED file without headers or index
    bed_annotated[["CHROM", "START", "END", "FILTERS"]].to_csv(f"{maf_df.name}_discarded.bed", sep="\t", header=False, index=False)


def main(maf_path: str, samp_name: str, repetitive_variant_thr: int, somatic_vaf_boundary: float):
    """
    Script to process a MAF (Mutation Annotation Format) file.
    It filters out repetitive variants, cohort_n_rich variants, and SNPs from other samples.

    Parameters
    ----------
    maf_df_file : str
        Path to the input MAF file (gzipped)
    samp_name : str
        Name for the output sample
    repetitive_variant_threshold : int
        Minimum number of occurrences to mark a variant as repetitive
    somatic_vaf_boundary : float
        VAF threshold to distinguish somatic mutations
    """
    LOG.info("Starting filtering cohort...")

    # Load data
    maf = pd.read_csv(maf_path, compression="gzip", header=0, sep="\t", na_values=custom_na_values)

    # Correct VAF -- NOT USED
    # maf_df = correct_vaf(maf_df)
    # print("VAF corrected")

    # Filter repetitive variants
    maf = filter_repetitive_variants(maf, repetitive_variant_thr, somatic_vaf_boundary)

    # Filter cohort_n_rich variants
    maf = filter_cohort_n_rich(maf, somatic_vaf_boundary)

    # Filter SNPs from other samples
    maf = filter_other_samples_snp(maf, somatic_vaf_boundary)

    # Expand FILTER column
    maf = explode_filter_column(maf)

    maf.name = samp_name

    # Save the final DataFrame
    maf.to_csv(f"{samp_name}.cohort.filtered.tsv.gz", sep="\t", header=True, index=False)

    # Get discarded mutations
    get_discarded_regions(maf)

    # TODO add basic logger instead of print
    LOG.info("Filtering complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a MAF (Mutation Annotation Format) file.")

    parser.add_argument("maf_path", type=str, help="Path to the input MAF file (gzipped)")
    parser.add_argument("samp_name", type=str, help="Name for the output sample")
    parser.add_argument("repetitive_variant_thr", type=int, help="Threshold for repetitive variants")
    parser.add_argument("somatic_vaf_boundary", type=float, help="VAF threshold to distinguish somatic mutations")

    # Parse the arguments
    args = parser.parse_args()

    # Call main function with parsed arguments
    main(args.maf_path, args.samp_name, args.repetitive_variant_thr, args.somatic_vaf_boundary)
