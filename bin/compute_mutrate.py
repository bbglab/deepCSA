#!/usr/bin/env python3

import pandas as pd
import sys

# -- Auxiliary functions -- #

def mutrate_sample(maf_df, depths_df):
    """
    computes mutation rate per Mb sequenced
    per sample
    """

    # muts per group
    nmuts_df = maf_df.groupby("SAMPLE_ID").size().to_frame("N_MUTS").reset_index()

    # depth per group
    depth_samples_df = depths_df[["CHROM", "POS", "SAMPLE_ID", "DEPTH"]].groupby("SAMPLE_ID")["DEPTH"].sum()

    # merge
    mutrate_df = nmuts_df.merge(depth_samples_df, on = "SAMPLE_ID")

    # calculate Mb mutation rate
    mutrate_df["MUTRATE_MB"] = (mutrate_df["N_MUTS"] / mutrate_df["DEPTH"] * 1000000).astype("float")
    mutrate_df["canonical_SYMBOL"] = "ALL_GENES"

    return mutrate_df

def mutrate_gene(maf_df, depths_df):
    """
    computes mutation rate per Mb sequenced
    per gene
    """

    # muts per group
    nmuts_df = maf_df.groupby("canonical_SYMBOL").size().to_frame("N_MUTS").reset_index()

    # depth per group
    depth_gene_df = depths_df.groupby("canonical_SYMBOL")["DEPTH"].sum()

    # merge
    mutrate_df = nmuts_df.merge(depth_gene_df, on = "canonical_SYMBOL")

    # calculate Mb mutation rate
    mutrate_df["MUTRATE_MB"] = (mutrate_df["N_MUTS"] / mutrate_df["DEPTH"] * 1000000).astype("float")
    mutrate_df["SAMPLE_ID"] = "ALL_SAMPLES"

    return mutrate_df

def mutrate_sample_gene(maf_df, depths_df):
    """
    computes mutation rate per Mb sequenced
    per sample and gene
    """

    # muts per group
    nmuts_df = maf_df.groupby(["canonical_SYMBOL", "SAMPLE_ID"]).size().to_frame("N_MUTS").reset_index()

    # depth per group
    depth_sample_gene_df = depths_df.groupby(["canonical_SYMBOL", "SAMPLE_ID"])["DEPTH"].sum()

    # merge
    mutrate_df = nmuts_df.merge(depth_sample_gene_df, on = ["canonical_SYMBOL", "SAMPLE_ID"])

    # calculate Mb mutation rate
    mutrate_df["MUTRATE_MB"] = (mutrate_df["N_MUTS"] / mutrate_df["DEPTH"] * 1000000).astype("float")

    return mutrate_df


# -- Main function -- #

def compute_mutrate(maf_path, depths_path, annot_panel_path, output_prefix):

    # File loading
    depths_df = pd.read_csv(depths_path, sep = "\t")
    annot_panel_df = pd.read_csv(annot_panel_path, sep = "\t")
    maf_df = pd.read_csv(maf_path, sep = "\t")

    # Subset MAF and depths with panel
    depths_df = depths_df.melt(id_vars = ["CHROM", "POS"], var_name = "SAMPLE_ID", value_name = "DEPTH")
    depths_df["SAMPLE_ID"] = depths_df.apply(lambda row: row["SAMPLE_ID"].split(".")[0], axis = 1)
    depths_subset_df = depths_df.merge(annot_panel_df[["CHROM", "POS", "GENE"]].drop_duplicates(),
                                       on = ["CHROM", "POS"], how = "inner").rename({"GENE": "canonical_SYMBOL"}, axis = 1)

    ## double filtering: panel and consequence (some parts of the panels intersect)
    maf_subset_df = maf_df.merge(annot_panel_df[["CHROM", "POS"]].drop_duplicates(),
                                 on = ["CHROM", "POS"], how = "inner")
    maf_subset_df = maf_subset_df.loc[maf_subset_df["Consequence_broader"].isin(annot_panel_df["IMPACT"].unique())]

    # Compute mutation rates
    mutrate_samples_df = mutrate_sample(maf_subset_df, depths_subset_df)
    mutrate_genes_df = mutrate_gene(maf_subset_df, depths_subset_df)
    mutrate_samples_genes_df = mutrate_sample_gene(maf_subset_df, depths_subset_df)
    mutrate_df = pd.concat([mutrate_samples_df, mutrate_genes_df, mutrate_samples_genes_df])

    # Save
    mutrate_df.to_csv(f"{output_prefix}.mutrates.tsv", sep = "\t", index = False)


if __name__ == '__main__':
    maf_path = sys.argv[1]
    depths_path = sys.argv[2]
    annot_panel_path = sys.argv[3]
    output_prefix = sys.argv[4]

    compute_mutrate(maf_path, depths_path, annot_panel_path, output_prefix)


