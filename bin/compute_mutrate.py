#!/usr/local/bin/python

import sys
import pandas as pd

# -- Auxiliary functions -- #

def mutrate_sample(maf_df, depths_df, sample_name, method = 'unique'):
    """
    computes mutation rate per Mb sequenced
    per sample
    """
    print(maf_df.shape)
    unique_maf = maf_df[["SAMPLE_ID", "MUT_ID", "ALT_DEPTH"]].drop_duplicates()
    print(unique_maf.shape)

    # make sure to count each mutation only once (avoid annotation issues)
    n_muts = unique_maf.shape[0]
    print(n_muts)

    # make sure to count each mutation only once (avoid annotation issues)
    n_muts_per_sample = unique_maf.groupby(by = ["SAMPLE_ID", "MUT_ID"] ).agg({"ALT_DEPTH" : "sum" }).reset_index()
    n_mutated_reads = n_muts_per_sample["ALT_DEPTH"].sum()
    n_mutated_reads2 = unique_maf["ALT_DEPTH"].sum()
    print(n_muts_per_sample, n_mutated_reads, n_mutated_reads2)


    sample_features = {"N_MUTS" : n_muts,
                        "N_MUTATED" : n_mutated_reads,
                        "DEPTH" : depths_df[f"{sample_name}"].sum()
                        }
    sample_features["MUTRATE_MB"] = ( sample_features["N_MUTS"] / sample_features["DEPTH"] * 1000000 ).astype(float)
    sample_features["MUTREADSRATE_MB"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH"] * 1000000 ).astype(float)
    sample_features["GENE"] = "ALL_GENES"

    return pd.DataFrame([sample_features])


def mutrate_gene(maf_df, depths_df, sample_name, method = 'unique'):
    """
    computes mutation rate per Mb sequenced
    per gene
    """
    unique_maf = maf_df[["SAMPLE_ID", "GENE", "MUT_ID", "ALT_DEPTH"]].drop_duplicates()

    # make sure to count each mutation only once (avoid annotation issues)
    n_muts_gene = unique_maf.groupby(by = ["GENE"] ).agg({"ALT_DEPTH" : "count" })
    n_muts_gene.columns = ["N_MUTS"]
    print(n_muts_gene)

    # make sure to count each mutation only once (avoid annotation issues)
    n_mutated_reads = unique_maf.groupby(by = ["GENE"] ).agg({"ALT_DEPTH" : "sum" })
    n_mutated_reads.columns = ["N_MUTATED"]
    print(n_mutated_reads)

    depths_gene_df = depths_df.groupby("GENE").agg({f"{sample_name}" : "sum" })
    depths_gene_df.columns = ["DEPTH"]

    print(n_muts_gene)
    print(depths_gene_df)

    mut_rate_mut_reads_df = n_muts_gene.merge(n_mutated_reads, on = "GENE")
    mut_depths_df = depths_gene_df.merge(mut_rate_mut_reads_df, on = "GENE")

    mut_depths_df["MUTRATE_MB"] = (mut_depths_df["N_MUTS"] / mut_depths_df["DEPTH"] * 1000000).astype(float)
    mut_depths_df["MUTREADSRATE_MB"] = (mut_depths_df["N_MUTATED"] / mut_depths_df["DEPTH"] * 1000000).astype(float)

    print(mut_depths_df)

    return mut_depths_df.reset_index()




# -- Main function -- #
def compute_mutrate(maf_path, depths_path, annot_panel_path, sample_name, panel_v):

    # File loading
    depths_df = pd.read_csv(depths_path, sep = "\t")
    annot_panel_df = pd.read_csv(annot_panel_path, sep = "\t")
    maf_df = pd.read_csv(maf_path, sep = "\t")

    depths_df = depths_df.drop("CONTEXT", axis = 1)

    # Subset MAF and depths with panel
    depths_subset_df = depths_df.merge(annot_panel_df[["CHROM", "POS", "GENE"]].drop_duplicates(),
                                        on = ["CHROM", "POS"], how = "inner")

    del depths_df
    del annot_panel_df

    ## double filtering: panel and consequence (some parts of the panels intersect)
    # maf_subset_df = maf_df.merge(annot_panel_df[["CHROM", "POS"]].drop_duplicates(),
    #                                 on = ["CHROM", "POS"], how = "inner")
    # maf_subset_df = maf_subset_df.loc[maf_subset_df["Consequence_broader"].isin(annot_panel_df["IMPACT"].unique())]

    # Compute mutation rates
    mutrate_samples_df = mutrate_sample(maf_df, depths_subset_df, sample_name.split('.')[0])
    mutrate_genes_df = mutrate_gene(maf_df, depths_subset_df, sample_name.split('.')[0])
    mutrate_df = pd.concat([mutrate_samples_df, mutrate_genes_df])
    mutrate_df["SAMPLE_ID"] = sample_name.split('.')[0]
    mutrate_df["REGIONS"] = panel_v

    # Save
    mutrate_df[["SAMPLE_ID", "GENE", "REGIONS",
                "N_MUTS", "N_MUTATED", "DEPTH",
                "MUTRATE_MB", "MUTREADSRATE_MB"]].to_csv(f"{sample_name}.{panel_v}.mutrates.tsv",
                                                            sep = "\t",
                                                            header = True,
                                                            index = False
                                                            )



if __name__ == '__main__':
    # TODO reimplement with click
    maf_path = sys.argv[1]
    depths_path = sys.argv[2]
    annot_panel_path = sys.argv[3]
    sample_name = sys.argv[4]
    panel_version = sys.argv[5]

    compute_mutrate(maf_path, depths_path, annot_panel_path, sample_name, panel_version)


