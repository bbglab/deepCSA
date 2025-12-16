#!/usr/bin/env python

"""
Filter Cohort - MAF Processing and Flagging Script

This script processes a Mutation Annotation Format (MAF) file to filter variants by specific criteria and 
generates a final filtered MAF along with an output of flagged regions in BED format.

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
# MISSING!

"""
import logging

import click
import pandas as pd
from utils import add_filter
from read_utils import custom_na_values

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("filter_cohort")

@click.command()
@click.option('--maf-df-file', required=True, type=click.Path(exists=True), help='Input gzipped MAF file (TSV)')
@click.option('--sample-name', required=True, type=str, help='Sample name for output file')
@click.option('--repetitive-variant-threshold', required=True, type=int, help='Threshold for repetitive variants')
@click.option('--somatic-vaf-boundary', required=True, type=float, help='VAF boundary for somatic variants')
@click.option('--n-rich-cohort-proportion', required=True, type=float, help='Proportion for n-rich cohort filtering')

def flag_repetitive_variants(maf_df, repetitive_variant_threshold, somatic_vaf_boundary) -> pd.DataFrame:
    """
    Filters out repetitive variants from the MAF dataframe. A variant is considered repetitive if it appears in at least
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

    # work with already filtered df + somatic only to explore potential artifacts
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
    maf_df["FILTER"] = maf_df[["FILTER","repetitive_mapping_variant"]].apply(lambda x: add_filter(x["FILTER"], x["repetitive_mapping_variant"], "repetitive_mapping_variant"),
                                                                        axis = 1
                                                                    )
    maf_df = maf_df.drop("repetitive_mapping_variant", axis = 1)

    return maf_df


def main(maf_df_file, sample_name, repetitive_variant_threshold, somatic_vaf_boundary, n_rich_cohort_proportion):
    maf_df = pd.read_csv(maf_df_file, compression='gzip', header=0, sep='\t', na_values=custom_na_values)
    sequenced_genes = list(pd.unique(maf_df["SYMBOL"]))
    
    #######
    ###  Filter repetitive variants,
    ###     both based on frequency and including information on position in read
    #######
    maf_df = flag_repetitive_variants(maf_df, repetitive_variant_threshold, somatic_vaf_boundary)

    #######
    ###  Filter cohort_n_rich
    #######

    max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))

    if max_samples < 2:
        print("Not enough samples to identify cohort_n_rich mutations!")

    else:
        number_of_samples = max(2, (max_samples * n_rich_cohort_proportion) // 1)
        print(f"flagging mutations that are n_rich in at least: {number_of_samples} samples as cohort_n_rich")

        # work with already filtered df + somatic only to explore potential artifacts
        # take only variant and sample info from the df
        maf_df_f_somatic = maf_df[["MUT_ID", "SAMPLE_ID", "VAF_Ns", "FILTER"]].reset_index(drop = True)

        n_rich_vars_df = maf_df_f_somatic[maf_df_f_somatic["FILTER"].str.contains("n_rich")].groupby("MUT_ID")[
                                                ['SAMPLE_ID', 'VAF_Ns']
                                            ].agg({'SAMPLE_ID' : len, 'VAF_Ns' : min})
        n_rich_vars_df = n_rich_vars_df.rename({'SAMPLE_ID' : 'N_rich_frequency', 'VAF_Ns' : 'VAF_Ns_threshold'}, axis = 'columns')

        n_rich_vars = list(n_rich_vars_df[n_rich_vars_df['N_rich_frequency'] >= number_of_samples].index)

        maf_df["cohort_n_rich"] = maf_df["MUT_ID"].isin(n_rich_vars)

        maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich"], "cohort_n_rich"),
                                                                            axis = 1
                                                                        )
        maf_df = maf_df.drop("cohort_n_rich", axis = 1)



        # if the variant appeared flagged as n_rich in a single sample it is also filtered out from all other samples
        n_rich_vars_uni = list(n_rich_vars_df[n_rich_vars_df['N_rich_frequency'] > 0].index)

        maf_df["cohort_n_rich_uni"] = maf_df["MUT_ID"].isin(n_rich_vars_uni)

        maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich_uni"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich_uni"], "cohort_n_rich_uni"),
                                                                            axis = 1
                                                                        )
        maf_df = maf_df.drop("cohort_n_rich_uni", axis = 1)


        # if the variant appeared flagged as n_rich in a single sample it is also filtered out from all other samples
        maf_df = maf_df.merge(n_rich_vars_df, on = 'MUT_ID', how = 'left')
        maf_df['N_rich_frequency'] = maf_df['N_rich_frequency'].fillna(0)
        maf_df['VAF_Ns_threshold'] = maf_df['VAF_Ns_threshold'].fillna(1.1)

        maf_df["cohort_n_rich_threshold"] = maf_df["VAF_Ns"] >= maf_df['VAF_Ns_threshold']

        maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich_threshold"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich_threshold"], "cohort_n_rich_threshold"),
                                                                            axis = 1
                                                                        )
        maf_df = maf_df.drop("cohort_n_rich_threshold", axis = 1)






    #######
    ###  Filter other sample's SNP
    #######

    # this is if we were to consider both unique and no-unique variants
    germline_vars_all_samples = maf_df.loc[(maf_df["VAF"] > somatic_vaf_boundary) &
                                        (maf_df["VAF_AM"] > somatic_vaf_boundary) &
                                        (maf_df["vd_VAF"] > somatic_vaf_boundary),
                                        "MUT_ID"].unique()
    print(len(germline_vars_all_samples), "using all germline variants of all samples")


    maf_df["other_sample_SNP"] = False
    maf_df.loc[(maf_df["MUT_ID"].isin(germline_vars_all_samples)) &
                (maf_df["VAF"] <= somatic_vaf_boundary), "other_sample_SNP"] = True

    maf_df["FILTER"] = maf_df[["FILTER","other_sample_SNP"]].apply(
                                                    lambda x: add_filter(x["FILTER"], x["other_sample_SNP"], "other_sample_SNP"),
                                                    axis = 1
                                                )
    maf_df = maf_df.drop("other_sample_SNP", axis = 1)



    #######
    ###  Filter gnomad SNP
    #######

    if "gnomAD_SNP" in maf_df.columns:
        maf_df["gnomAD_SNP"] = maf_df["gnomAD_SNP"].replace({"True": True, "False": False, '-' : False}).fillna(False).astype(bool)
        print("Out of ", maf_df["gnomAD_SNP"].shape[0], "positions", maf_df["gnomAD_SNP"].sum(), "are gnomAD SNPs (>0.1)")
        maf_df["FILTER"] = maf_df[["FILTER","gnomAD_SNP"]].apply(
                                                                    lambda x: add_filter(x["FILTER"], x["gnomAD_SNP"], "gnomAD_SNP"),
                                                                    axis = 1
                                                                )


    for filt in pd.unique(maf_df["FILTER"].str.split(";").explode()):
        maf_df[f"FILTER.{filt}"] = maf_df["FILTER"].apply(lambda x: filt in x.split(";"))

    for filtt in [ "not_covered", "not_in_exons"]:
        if f"FILTER.{filtt}" not in maf_df.columns:
            maf_df[f"FILTER.{filtt}"] = False


    maf_df.to_csv(f"{sample_name}.cohort.filtered.tsv.gz",
                  sep = "\t",
                  header = True,
                  index = False)

if __name__ == '__main__':
    main()
