#!/usr/bin/env python

"""
Filter Cohort - MAF Processing and Flagging Script

This script processes a Mutation Annotation Format (MAF) file to filter variants by specific criteria and 
generates a final filtered MAF.

Command-line Arguments
----------------------
maf_df-file : str
    Path to the gzipped input MAF file.
sample_name : str
    Output sample name.
repetitive-variant-threshold : int
    Minimum occurrences threshold to flag a repetitive variant.
somatic-vaf-boundary : float
    VAF threshold to classify somatic mutations.
n-rich-cohort-proportion : float
    Proportion threshold for n-rich cohort filtering.

Authors
-------
Author  : Ferriol Calvet (@FerriolCalvet)
Email   : ferriol.calvet@irbbarcelona.org

Contributors
------------
- Raquel Blanco - @rblancomi (raquel.blanco@irbbarcelona.org)
- Federica Brando - @FedericaBrando (federica.brando@irbbarcelona.org)
- Marta Huertas - @m-huertasp (marta.huertas@irbbarcelona.org)

Usage
-----
python filter_cohort.py \\
    --maf-df-file input.maf.tsv.gz \\
    --sample-name sample1 \\
    --repetitive-variant-threshold 5 \\
    --somatic-vaf-boundary 0.4 \\
    --n-rich-cohort-proportion 0.1

"""
import logging
from pathlib import Path

import click
import pandas as pd
from utils import add_filter
from read_utils import custom_na_values
from utils_filter import expand_filter_column

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("filter_cohort")

def flag_repetitive_variants(maf_df: pd.DataFrame,
                             repetitive_variant_threshold: int,
                             somatic_vaf_boundary: float) -> pd.DataFrame:
    """
    Flags filter column for repetitive variants from the MAF dataframe. A variant is considered repetitive if it appears in at least
    ``repetitive_variant_threshold`` samples. Additionally, variants that consistently appear at the same position in reads
    are also flagged as ``repetitive_mapping_variant``.

    Parameters
    ----------
    maf_df : pd.DataFrame
        MAF dataframe
    repetitive_variant_threshold : int
        Minimum number of samples a variant must appear in to be considered repetitive
    somatic_vaf_boundary : float
        VAF threshold to classify somatic mutations

    Returns
    -------
    pandas.DataFrame
        MAF dataframe with a new column 'repetitive_variant' that flags repetitive variants
    """
    max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))

    n_samples = list(range(repetitive_variant_threshold, max_samples + 1))
    if n_samples == 0:
        LOG.warning("Not enough samples to identify potential repetitive variants!")

        return maf_df

    # Work with already filtered df + somatic only to explore potential artifacts
    # take only variant and sample info from the df
    maf_df_f_somatic = maf_df.loc[maf_df["VAF"] <= somatic_vaf_boundary][["MUT_ID","SAMPLE_ID", "PMEAN", "PSTD"]].reset_index(drop = True)

    # Group by 'MUT_ID' and count occurrences
    maf_df_f_somatic_pivot = maf_df_f_somatic.groupby("MUT_ID").size().reset_index(name="count")

    # Store repetitive variants
    repetitive_variants = maf_df_f_somatic_pivot[maf_df_f_somatic_pivot["count"] >= repetitive_variant_threshold]["MUT_ID"]
    LOG.info("%s repetitive_variants", len(repetitive_variants))

    # Flag repetitive variants in the original dataframe
    maf_df["repetitive_variant"] = maf_df["MUT_ID"].isin(repetitive_variants)
    maf_df["FILTER"] = maf_df[["FILTER", "repetitive_variant"]].apply(
        lambda x: add_filter(x["FILTER"], x["repetitive_variant"], "repetitive_variant"), axis=1
    )
    maf_df = maf_df.drop("repetitive_variant", axis=1)

    # Use the position in read information to filter repetitive variants with a fixed position (likely artifacts)
    maf_df_f_somatic_pos_info = maf_df_f_somatic[~(maf_df_f_somatic["PMEAN"].isna()) & 
                                                    (maf_df_f_somatic["PMEAN"] != -1) &
                                                    (maf_df_f_somatic["PSTD"] == 0)]
    
    # Check if there are any repetitive variants with a fixed position
    if maf_df_f_somatic_pos_info.shape[0] == 0:
        LOG.info("No repetitive variants with fixed position found.")
        return maf_df

    # Count unique PMEAN values for each MUT_ID
    maf_df_f_somatic_compiled_pos = maf_df_f_somatic_pos_info.groupby("MUT_ID")["PMEAN"].nunique().reset_index()

    # Identify variants always found in the same position -> only one PMEAN value
    variants_with_rep_position = maf_df_f_somatic_compiled_pos[(maf_df_f_somatic_compiled_pos["PMEAN"] == 1)]["MUT_ID"]
    LOG.info("Variants always found in the same position: %d", len(variants_with_rep_position))

    # Intersect with previously identified repetitive variants
    variants_with_rep_position = set(variants_with_rep_position).intersection(set(repetitive_variants))
    LOG.info("Repetitive variants always found in the same position: %d", len(variants_with_rep_position))

    # Flag these variants in the maf dataframe
    maf_df["repetitive_mapping_variant"] = maf_df["MUT_ID"].isin(variants_with_rep_position)
    LOG.info("%s muts flagged as repetitive_mapping_variant", maf_df["repetitive_mapping_variant"].sum())
    
    maf_df["FILTER"] = maf_df[["FILTER","repetitive_mapping_variant"]].apply(lambda x: add_filter(x["FILTER"], x["repetitive_mapping_variant"], "repetitive_mapping_variant"),
                                                                        axis = 1
                                                                    )
    maf_df = maf_df.drop("repetitive_mapping_variant", axis = 1)

    return maf_df

def flag_cohort_n_rich(maf_df: pd.DataFrame,
                       n_rich_cohort_proportion: float,
                       somatic_vaf_boundary: float) -> pd.DataFrame:
    """
    Flags FILTER column for cohort_n_rich variants from the MAF dataframe. 

    Parameters
    ----------
    maf_df : pandas.DataFrame
        MAF dataframe
    n_rich_cohort_proportion : float
        Proportion of samples to consider a variant as cohort_n_rich
    somatic_vaf_boundary : float
        VAF threshold to classify somatic mutations

    Returns
    -------
    maf_df : pandas.DataFrame
        MAF dataframe with cohort_n_rich variants flagged
    """
    LOG.info("Flagging cohort_n_rich...")

    max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))
    if max_samples < 2:
        LOG.warning("Not enough samples to identify cohort_n_rich mutations!")
        return maf_df
    
    number_of_samples = max(2, (max_samples * n_rich_cohort_proportion) // 1)
    LOG.info(f"Flagging mutations that are n_rich in at least: {number_of_samples} samples as cohort_n_rich")

    # Work with already filtered df to explore potential artifacts
    # take only variant and sample info from the df.
    maf_df_f = maf_df[["MUT_ID", "SAMPLE_ID", "VAF_Ns", "FILTER"]].reset_index(drop = True)

    # Aggregate n_rich variants
    n_rich_vars_df = (
        maf_df_f[maf_df_f["FILTER"].str.contains("n_rich")]
        .groupby("MUT_ID")
        .agg(
            N_rich_frequency=('SAMPLE_ID', 'count'),
            VAF_Ns_threshold=('VAF_Ns', 'min')
        )
        )
    
    # Flag variants that are n_rich in at least number_of_samples samples -> cohort_n_rich
    n_rich_vars = set(n_rich_vars_df[n_rich_vars_df['N_rich_frequency'] >= number_of_samples].index)

    maf_df["cohort_n_rich"] = maf_df["MUT_ID"].isin(n_rich_vars)
    LOG.info("%s muts flagged as cohort_n_rich", maf_df["cohort_n_rich"].sum())

    maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich"], "cohort_n_rich"),
                                                                            axis = 1
                                                                        )
    
    # Flag variants that are n_rich in at least 1 sample -> cohort_n_rich_uni
    n_rich_vars_uni = set(n_rich_vars_df[n_rich_vars_df['N_rich_frequency'] > 0].index)

    maf_df["cohort_n_rich_uni"] = maf_df["MUT_ID"].isin(n_rich_vars_uni)
    LOG.info("%s muts flagged as cohort_n_rich_uni", maf_df["cohort_n_rich_uni"].sum())

    maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich_uni"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich_uni"], "cohort_n_rich_uni"),
                                                                        axis = 1
                                                                    )
    
    # Flag variants that exceed the VAF_Ns threshold -> cohort_n_rich_threshold
    maf_df = maf_df.merge(n_rich_vars_df, on = 'MUT_ID', how = 'left')
    maf_df['N_rich_frequency'] = maf_df['N_rich_frequency'].fillna(0)
    maf_df['VAF_Ns_threshold'] = maf_df['VAF_Ns_threshold'].fillna(1.1)

    maf_df["cohort_n_rich_threshold"] = maf_df["VAF_Ns"] >= maf_df['VAF_Ns_threshold']
    LOG.info("%s muts flagged as cohort_n_rich_threshold", maf_df["cohort_n_rich_threshold"].sum())

    maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich_threshold"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich_threshold"], "cohort_n_rich_threshold"),
                                                                        axis = 1
                                                                    )
    # Drop temporary columns
    maf_df = maf_df.drop(["cohort_n_rich", "cohort_n_rich_uni", "cohort_n_rich_threshold"], axis = 1)
 
    return maf_df

def flag_other_samples_snp(maf_df,
                           somatic_vaf_boundary: float) -> pd.DataFrame:
    """
    Filters out SNPs from other samples from the MAF dataframe

    Parameters
    ----------
    maf_df : pandas.DataFrame
        MAF dataframe
    somatic_vaf_boundary : float
        VAF boundary to consider a variant as somatic

    Returns
    -------
    maf_df : pandas.DataFrame
        MAF dataframe with other sample SNPs flagged
    """
    LOG.info("Flagging SNPs from other samples...")
    # Get all germline variants from all samples, consider both unique and non-unique variants
    germline_vars_all_samples = maf_df.loc[(maf_df["VAF"] > somatic_vaf_boundary) &
                                    (maf_df["VAF_AM"] > somatic_vaf_boundary) &
                                    (maf_df["vd_VAF"] > somatic_vaf_boundary),
                                    "MUT_ID"].unique()
    
    LOG.info(f"Using all germline variants of all samples, total: {len(germline_vars_all_samples)} variants.")

    # Identify variants that are germline in other samples but somatic in the current sample
    maf_df["other_sample_SNP"] = False
    maf_df.loc[(maf_df["MUT_ID"].isin(germline_vars_all_samples)) &
            (maf_df["VAF"] <= somatic_vaf_boundary), "other_sample_SNP"] = True
    LOG.info("%s muts flagged as other_sample_SNP", maf_df['other_sample_SNP'].sum())

    # Flag variants that are germline in other samples but somatic in the current sample
    maf_df["FILTER"] = maf_df[["FILTER","other_sample_SNP"]].apply(
                                                    lambda x: add_filter(x["FILTER"], x["other_sample_SNP"], "other_sample_SNP"),
                                                    axis = 1
                                                )
    maf_df = maf_df.drop("other_sample_SNP", axis = 1)

    return maf_df

def flag_gnomad_snp(maf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Flags gnomAD SNPs in the MAF dataframe

    Parameters
    ----------
    maf_df : pandas.DataFrame
        MAF dataframe

    Returns
    -------
    maf_df : pandas.DataFrame
        MAF dataframe with gnomAD SNPs flagged
    """
    LOG.info("Flagging gnomAD SNPs...")

    # Flag gnomAD SNPs
    if "gnomAD_SNP" in maf_df.columns:
        maf_df["gnomAD_SNP"] = maf_df["gnomAD_SNP"].replace({"True": True, "False": False, '-' : False}).fillna(False).astype(bool)
        LOG.info("Out of %d positions, %d are gnomAD SNPs", maf_df["gnomAD_SNP"].shape[0], maf_df["gnomAD_SNP"].sum())
        
        maf_df["FILTER"] = maf_df[["FILTER","gnomAD_SNP"]].apply(
                                                                    lambda x: add_filter(x["FILTER"], x["gnomAD_SNP"], "gnomAD_SNP"),
                                                                    axis = 1
                                                                )
        maf_df = maf_df.drop("gnomAD_SNP", axis = 1)

    return maf_df

def flag_vaf_ns_threshold(maf_df: pd.DataFrame, vaf_ns_threshold: float) -> pd.DataFrame:

    """
    Flag variants that have a proportion of Ns higher than vaf_ns_threshold

    Parameters
    ----------
    maf_df : pandas.DataFrame
        MAF dataframe
    vaf_ns_threshold : float
        VAF of Ns threshold for filtering variants

    Returns
    -------
    maf_df : pandas.DataFrame
        MAF dataframe with variants exceeding VAF_Ns threshold flagged
    """
    LOG.info("Flagging variants with proportion of Ns higher than VAF_Ns threshold...")

    maf_df["high_n_vaf"] = maf_df[["VAF_Ns", "VAF_Ns_AM"]].ge(vaf_ns_threshold).any(axis=1)
    LOG.info("%s muts flagged as high_n_vaf", maf_df["high_n_vaf"].sum())

    maf_df["FILTER"] = maf_df[["FILTER","high_n_vaf"]].apply(
                                                                lambda x: add_filter(x["FILTER"], x["high_n_vaf"], "high_n_vaf"),
                                                                axis = 1
                                                            )
    maf_df = maf_df.drop("high_n_vaf", axis = 1)

    return maf_df

def flag_distorted_expanded(maf_df: pd.DataFrame) -> pd.DataFrame:
    """
    If there is a column named VAF_distorted_expanded_sq, add a filter flag for variants with distorted VAF distribution.

    Parameters
    ----------
    maf_df : pandas.DataFrame
        MAF dataframe

    Returns
    -------
    maf_df : pandas.DataFrame
        MAF dataframe with variants having distorted VAF distribution flagged
    """
    LOG.info("Flagging variants with distorted VAF distribution...")

    if 'VAF_distorted_expanded_sq' in maf_df.columns:
        maf_df["FILTER"] = maf_df[["FILTER","VAF_distorted_expanded_sq"]].apply(
                                                lambda x: add_filter(x["FILTER"], x["VAF_distorted_expanded_sq"], "VAF_distorted_expanded_sq"),
                                                axis = 1
                                                                                )

    return maf_df

def flag_maf(maf_df: pd.DataFrame, sample_name: str,
               repetitive_variant_threshold: int,
               somatic_vaf_boundary: float,
               n_rich_cohort_proportion: float,
               vaf_ns_threshold: float) -> None:
    """
    Script to process a MAF (Mutation Annotation Format) file.
    It filters out repetitive variants, cohort_n_rich variants, and SNPs from other samples.

    Parameters
    ----------
    maf_df : pd.DataFrame
        Input MAF dataframe
    sample_name : str
        Name for the output sample
    repetitive_variant_threshold : int
        Minimum number of occurrences to flag a variant as repetitive
    somatic_vaf_boundary : float
        VAF threshold to distinguish somatic mutations
    n_rich_cohort_proportion : float
        Proportion threshold for n-rich cohort filtering
    vaf_ns_threshold : float
        VAF of Ns threshold for filtering variants
    """
    ## Filter repetitive variants,
    # both based on frequency and including information on position in read
    maf_df = flag_repetitive_variants(maf_df, repetitive_variant_threshold, somatic_vaf_boundary)

    ## Filter cohort_n_rich variants
    maf_df = flag_cohort_n_rich(maf_df, n_rich_cohort_proportion, somatic_vaf_boundary)

    ## Filter variants with high VAF of Ns
    maf_df = flag_vaf_ns_threshold(maf_df, vaf_ns_threshold)

    ## Filter other sample's SNP
    maf_df = flag_other_samples_snp(maf_df, somatic_vaf_boundary)

    ## Filter gnomad SNP
    maf_df = flag_gnomad_snp(maf_df)

    ## Optionally, flag variants with distorted VAF distribution if the corresponding column is present
    maf_df = flag_distorted_expanded(maf_df)

    ## Expand FILTER column
    maf_df = expand_filter_column(maf_df)

    ## Save final filtered MAF
    maf_df.to_csv(f"{sample_name}.cohort.filtered.tsv.gz",
                  sep = "\t",
                  header = True,
                  index = False)
    
    LOG.info("Cohort flagging complete!")

@click.command()
@click.option('--maf-df-file', required=True, type=click.Path(exists=True), help='Input gzipped MAF file (TSV)')
@click.option('--sample-name', required=True, type=str, help='Sample name for output file')
@click.option('--repetitive-variant-threshold', required=True, type=int, help='Threshold for repetitive variants')
@click.option('--somatic-vaf-boundary', required=True, type=float, help='VAF boundary for somatic variants')
@click.option('--n-rich-cohort-proportion', required=True, type=float, help='Proportion for n-rich cohort filtering')
@click.option('--vaf-ns-threshold', required=False, type=float, default=0.1, help='VAF of Ns threshold for filtering variants')
def main(maf_df_file: str, sample_name: str, repetitive_variant_threshold: int,
         somatic_vaf_boundary: float, n_rich_cohort_proportion: float, vaf_ns_threshold: float):
    """
    CLI wrapper for flag_maf function.
    """
    # Load MAF dataframe
    maf_df = pd.read_csv(maf_df_file, compression='gzip', header=0, sep='\t', na_values=custom_na_values)
    LOG.debug(f"{maf_df_file}")
    # Flag MAF file
    flag_maf(maf_df,
        sample_name,
        repetitive_variant_threshold,
        somatic_vaf_boundary, 
        n_rich_cohort_proportion,
        vaf_ns_threshold)
    

if __name__ == '__main__':
    main()