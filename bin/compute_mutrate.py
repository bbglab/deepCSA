#!/usr/local/bin/python

import sys
import pandas as pd

# -- Auxiliary functions -- #

def mutrate_sample(maf_df, depths_df, depths_adj_df, sample_name, type_list = False):
    """
    Computes a sample's global mutation rate. Returns the mutation rate
    both per Mb and Kb sequenced, both non-adjusted and adjusted by panel
    composition.
    """

    # filter by mutation type according to type_list
    if not type_list:
        unique_maf = maf_df[["SAMPLE_ID", "MUT_ID", "ALT_DEPTH"]].drop_duplicates()
        types_included = 'all_types'
    else:
        maf_df = maf_df[maf_df['TYPE'].isin(type_list)]
        unique_maf = maf_df[["SAMPLE_ID", "MUT_ID", "ALT_DEPTH"]].drop_duplicates()
        types_included = '-'.join(sorted(type_list))

    # count number of mutations and mutated reads in the sample
    ## make sure to count each mutation only once (avoid annotation issues)
    n_muts = unique_maf.shape[0]
    ## make sure to count each mutation only once (avoid annotation issues)
    n_muts_per_sample = unique_maf.groupby(by = ["SAMPLE_ID", "MUT_ID"] ).agg({"ALT_DEPTH" : "sum" }).reset_index()
    n_mutated_reads = n_muts_per_sample["ALT_DEPTH"].sum()
    n_mutated_reads2 = unique_maf["ALT_DEPTH"].sum()
    print(n_muts_per_sample, n_mutated_reads, n_mutated_reads2)

    # mutation rate metrics
    sample_features = {"N_MUTS" : n_muts,
                        "N_MUTATED" : n_mutated_reads,
                        "DEPTH" : depths_df[f"{sample_name}"].sum(),
                        "DEPTH_ADJUSTED": depths_adj_df[f"{sample_name}"].sum()
                        }
    sample_features["MUTRATE_MB"] = ( sample_features["N_MUTS"] / sample_features["DEPTH"] * 1000000 ).astype(float)
    sample_features["MUTRATE_MB_ADJUSTED"] = ( sample_features["N_MUTS"] / sample_features["DEPTH_ADJUSTED"] * 1000000 ).astype(float)
    sample_features["MUTREADSRATE_MB"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH"] * 1000000 ).astype(float)
    sample_features["MUTREADSRATE_MB_ADJUSTED"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH_ADJUSTED"] * 1000000 ).astype(float)

    sample_features["MUTRATE_KB"] = ( sample_features["N_MUTS"] / sample_features["DEPTH"] * 1000 ).astype(float)
    sample_features["MUTRATE_KB_ADJUSTED"] = ( sample_features["N_MUTS"] / sample_features["DEPTH_ADJUSTED"] * 1000 ).astype(float)
    sample_features["MUTREADSRATE_KB"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH"] * 1000 ).astype(float)
    sample_features["MUTREADSRATE_KB_ADJUSTED"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH_ADJUSTED"] * 1000 ).astype(float)

    sample_features["GENE"] = "ALL_GENES"
    sample_features["MUTTYPES"] = types_included


    return pd.DataFrame([sample_features])


def mutrate_gene(maf_df, depths_df, depths_adj_df, sample_name, type_list = False):
    """
    Computes each gene mutation rate. Returns the mutation rate
    both per Mb and Kb sequenced, both non-adjusted and adjusted by panel
    composition.
    """

    # filter by mutation type according to type_list
    if not type_list:
        unique_maf = maf_df[["SAMPLE_ID", "GENE", "MUT_ID", "ALT_DEPTH"]].drop_duplicates()
        types_included = 'all_types'
    else:
        maf_df = maf_df[maf_df['TYPE'].isin(type_list)]
        unique_maf = maf_df[["SAMPLE_ID", "GENE", "MUT_ID", "ALT_DEPTH"]].drop_duplicates()
        types_included = '-'.join(sorted(type_list))

    # count number of mutations and mutated reads per gene
    # make sure to count each mutation only once (avoid annotation issues)
    n_muts_gene = unique_maf.groupby(by = ["GENE"] ).agg({"ALT_DEPTH" : "count" })
    n_muts_gene.columns = ["N_MUTS"]

    # make sure to count each mutation only once (avoid annotation issues)
    n_mutated_reads = unique_maf.groupby(by = ["GENE"] ).agg({"ALT_DEPTH" : "sum" })
    n_mutated_reads.columns = ["N_MUTATED"]

    depths_gene_df = depths_df.groupby("GENE").agg({f"{sample_name}" : "sum" })
    depths_gene_df.columns = ["DEPTH"]
    depths_adj_gene_df = depths_adj_df.groupby("GENE").agg({f"{sample_name}" : "sum" })
    depths_adj_gene_df.columns = ["DEPTH_ADJUSTED"]

    mut_rate_mut_reads_df = n_muts_gene.merge(n_mutated_reads, on = "GENE")
    depths_depthsadj_gene_df = depths_gene_df.merge(depths_adj_gene_df, on = "GENE")
    ## merge so that mutation rate is computed although the number of mutations is NA (meaning, zero)
    mut_depths_df = depths_depthsadj_gene_df.merge(mut_rate_mut_reads_df, on = "GENE", how = 'left')
    mut_depths_df = mut_depths_df.fillna(0) # I think this is not needed

    # mutation rate metrics
    mut_depths_df["MUTRATE_MB"] = (mut_depths_df["N_MUTS"] / mut_depths_df["DEPTH"] * 1000000).astype(float)
    mut_depths_df["MUTRATE_MB_ADJUSTED"] = (mut_depths_df["N_MUTS"] / mut_depths_df["DEPTH_ADJUSTED"] * 1000000).astype(float)
    mut_depths_df["MUTRATE_KB"] = (mut_depths_df["N_MUTS"] / mut_depths_df["DEPTH"] * 1000).astype(float)
    mut_depths_df["MUTRATE_KB_ADJUSTED"] = (mut_depths_df["N_MUTS"] / mut_depths_df["DEPTH_ADJUSTED"] * 1000).astype(float)

    mut_depths_df["MUTREADSRATE_MB"] = (mut_depths_df["N_MUTATED"] / mut_depths_df["DEPTH"] * 1000000).astype(float)
    mut_depths_df["MUTREADSRATE_MB_ADJUSTED"] = (mut_depths_df["N_MUTATED"] / mut_depths_df["DEPTH_ADJUSTED"] * 1000000).astype(float)
    mut_depths_df["MUTREADSRATE_KB"] = (mut_depths_df["N_MUTATED"] / mut_depths_df["DEPTH"] * 1000).astype(float)
    mut_depths_df["MUTREADSRATE_KB_ADJUSTED"] = (mut_depths_df["N_MUTATED"] / mut_depths_df["DEPTH_ADJUSTED"] * 1000).astype(float)

    mut_depths_df["MUTTYPES"] = types_included

    return mut_depths_df.reset_index()



# -- Main function -- #
def compute_mutrate(maf_path, depths_path, annot_panel_path, sample_name, panel_v):

    # File loading
    maf_df = pd.read_csv(maf_path, sep = "\t")
    depths_df = pd.read_csv(depths_path, sep = "\t")
    depths_df = depths_df.drop("CONTEXT", axis = 1)
    annot_panel_df = pd.read_csv(annot_panel_path, sep = "\t")

    # Subset depths with panel
    ## mode 1: each position counts one
    depths_subset_df = depths_df.merge(annot_panel_df[["CHROM", "POS", "GENE"]].drop_duplicates(),
                                        on = ["CHROM", "POS"], how = "inner")
    ## mode 2 (adjusted): each position counts as many times it contributes to the panel
    depths_df[sample_name.split('.')[0]] = depths_df[sample_name.split('.')[0]] / 3   # the depth per position can contribute to three different mutations
    depths_subset_adj_df = depths_df.merge(annot_panel_df[["CHROM", "POS", "GENE"]], on = ["CHROM", "POS"], how = "inner")

    del depths_df
    del annot_panel_df

    # Compute mutation rates
    ## sample mutation rate
    mutrate_sample_allmuts_df = mutrate_sample(maf_df, depths_subset_df, depths_subset_adj_df, sample_name.split('.')[0])
    mutrate_sample_snvs_df = mutrate_sample(maf_df, depths_subset_df, depths_subset_adj_df, sample_name.split('.')[0], ["SNV"])
    mutrate_sample_nonsnvs_df = mutrate_sample(maf_df, depths_subset_df, depths_subset_adj_df, sample_name.split('.')[0], ["INSERTION", "DELETION", "COMPLEX", "MNV"])
    mutrate_sample_indels_df = mutrate_sample(maf_df, depths_subset_df, depths_subset_adj_df, sample_name.split('.')[0], ["INSERTION", "DELETION"])
    ## per gene mutation rate
    mutrate_genes_allmuts_df = mutrate_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name.split('.')[0])
    mutrate_genes_snvs_df = mutrate_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name.split('.')[0], ["SNV"])
    mutrate_genes_nonsnvs_df = mutrate_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name.split('.')[0], ["INSERTION", "DELETION", "COMPLEX", "MNV"])
    mutrate_genes_indels_df = mutrate_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name.split('.')[0], ["INSERTION", "DELETION"])

    mutrate_df = pd.concat([mutrate_sample_allmuts_df, mutrate_sample_snvs_df, mutrate_sample_nonsnvs_df, mutrate_sample_indels_df,
                            mutrate_genes_allmuts_df, mutrate_genes_snvs_df, mutrate_genes_nonsnvs_df, mutrate_genes_indels_df])
    mutrate_df["SAMPLE_ID"] = sample_name.split('.')[0]
    mutrate_df["REGIONS"] = panel_v

    # Save
    mutrate_df[["SAMPLE_ID", "GENE", "REGIONS", "MUTTYPES",
                "N_MUTS", "N_MUTATED", "DEPTH",
                "MUTRATE_MB", "MUTRATE_MB_ADJUSTED", "MUTRATE_KB", "MUTREADSRATE_KB_ADJUSTED",
                "MUTREADSRATE_MB", "MUTREADSRATE_MB_ADJUSTED", "MUTREADSRATE_KB", "MUTREADSRATE_KB_ADJUSTED"]].to_csv(f"{sample_name.split('.')[0]}.{panel_v}.mutrates.tsv",
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


