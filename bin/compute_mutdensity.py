#!/usr/bin/env python

import click
import pandas as pd
from read_utils import custom_na_values

# TODO: bump pandas to 2.2.3

# -- Auxiliary functions -- #

def mutdensity_sample(maf_df, depths_df, depths_adj_df, sample_name, type_list = False):
    """
    Computes a sample's global mutation density. Returns the mutation density
    per Mb, non-adjusted and adjusted by panel
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
    n_muts_per_sample = unique_maf.groupby(by = ["SAMPLE_ID", "MUT_ID"] ).agg({"ALT_DEPTH" : "sum" }).reset_index()
    n_mutated_reads = n_muts_per_sample["ALT_DEPTH"].sum()
    n_mutated_reads2 = unique_maf["ALT_DEPTH"].sum()
    print(n_muts_per_sample, n_mutated_reads, n_mutated_reads2)

    # mutation density metrics
    sample_features = {"N_MUTS" : n_muts,
                        "N_MUTATED" : n_mutated_reads,
                        "DEPTH" : depths_df.drop_duplicates(subset = ["CHROM", "POS"])[f"{sample_name}"].sum(),
                        "DEPTH_ADJUSTED": depths_adj_df[f"{sample_name}"].sum() # they should be the same for all impacts not for subsets of impacts
                        }
    sample_features["MUTDENSITY_MB"] = ( sample_features["N_MUTS"] / sample_features["DEPTH"] * 1000000 ).astype(float)
    sample_features["MUTDENSITY_MB_ADJUSTED"] = ( sample_features["N_MUTS"] / sample_features["DEPTH_ADJUSTED"] * 1000000 ).astype(float)
    sample_features["MUTREADSRATE_MB"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH"] * 1000000 ).astype(float)
    sample_features["MUTREADSRATE_MB_ADJUSTED"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH_ADJUSTED"] * 1000000 ).astype(float)

    sample_features["GENE"] = "ALL_GENES"
    sample_features["MUTTYPES"] = types_included


    return pd.DataFrame([sample_features])


def mutdensity_gene(maf_df, depths_df, depths_adj_df, sample_name, type_list = False):
    """
    Computes each gene mutation density. Returns the mutation density
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
    ## merge so that mutation density is computed although the number of mutations is NA (meaning, zero)
    mut_depths_df = depths_depthsadj_gene_df.merge(mut_rate_mut_reads_df, on = "GENE", how = 'left')
    mut_depths_df = mut_depths_df.fillna(0) # I think this is not needed

    # mutation density metrics
    mut_depths_df["MUTDENSITY_MB"] = (mut_depths_df["N_MUTS"] / mut_depths_df["DEPTH"] * 1000000).astype(float)
    mut_depths_df["MUTDENSITY_MB_ADJUSTED"] = (mut_depths_df["N_MUTS"] / mut_depths_df["DEPTH_ADJUSTED"] * 1000000).astype(float)

    mut_depths_df["MUTREADSRATE_MB"] = (mut_depths_df["N_MUTATED"] / mut_depths_df["DEPTH"] * 1000000).astype(float)
    mut_depths_df["MUTREADSRATE_MB_ADJUSTED"] = (mut_depths_df["N_MUTATED"] / mut_depths_df["DEPTH_ADJUSTED"] * 1000000).astype(float)

    mut_depths_df["MUTTYPES"] = types_included

    return mut_depths_df.reset_index()



# -- Main function -- #
def compute_mutdensity(maf_path, depths_path, annot_panel_path, sample_name, panel_v):

    # File loading
    maf_df = pd.read_csv(maf_path, sep = "\t", na_values = custom_na_values)
    depths_df = pd.read_csv(depths_path, sep = "\t")
    depths_df = depths_df.drop("CONTEXT", axis = 1)
    annot_panel_df = pd.read_csv(annot_panel_path, sep = "\t", na_values = custom_na_values)

    # Subset depths with panel
    ## mode 1: each position counts one (once per gene, be careful that it might be duplicated in different genes)
    depths_subset_df = depths_df.merge(annot_panel_df[["CHROM", "POS", "GENE"]].drop_duplicates(),
                                        on = ["CHROM", "POS"], how = "inner")
    ## mode 2 (adjusted): each position counts as many times it contributes to the panel
    depths_df[sample_name] = depths_df[sample_name] / 3   # the depth per position can contribute to three different mutations
    depths_subset_adj_df = depths_df.merge(annot_panel_df[["CHROM", "POS", "GENE"]], on = ["CHROM", "POS"], how = "inner")

    ## mode 3 (adjusted): each position counts as many times it contributes to the panel, but ONLY ONCE PER SAMPLE
    depths_subset_adj_sample_df = depths_df.merge(annot_panel_df.drop_duplicates(subset = ["CHROM", "POS", "REF", "ALT"])[["CHROM", "POS"]],
                                                    on = ["CHROM", "POS"], how = "inner")

    del depths_df
    del annot_panel_df

    # Compute mutation densities
    ## sample mutation density
    mutdensity_sample_allmuts_df = mutdensity_sample(maf_df, depths_subset_df, depths_subset_adj_sample_df, sample_name)
    mutdensity_sample_snvs_df = mutdensity_sample(maf_df, depths_subset_df, depths_subset_adj_sample_df, sample_name, ["SNV"])
    mutdensity_sample_indels_df = mutdensity_sample(maf_df, depths_subset_df, depths_subset_adj_sample_df, sample_name, ["INSERTION", "DELETION"])
    mutdensity_sample_snvsnindels_df = mutdensity_sample(maf_df, depths_subset_df, depths_subset_adj_sample_df, sample_name, ["SNV", "INSERTION", "DELETION"])

    ## per gene mutation density
    mutdensity_genes_allmuts_df = mutdensity_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name)
    mutdensity_genes_snvs_df = mutdensity_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name, ["SNV"])
    mutdensity_genes_indels_df = mutdensity_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name, ["INSERTION", "DELETION"])
    mutdensity_genes_snvsnindels_df = mutdensity_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name, ["SNV", "INSERTION", "DELETION"])

    mutdensity_df = pd.concat([mutdensity_sample_allmuts_df, mutdensity_sample_snvs_df, mutdensity_sample_snvsnindels_df, mutdensity_sample_indels_df,
                            mutdensity_genes_allmuts_df, mutdensity_genes_snvs_df, mutdensity_genes_snvsnindels_df, mutdensity_genes_indels_df])
    mutdensity_df["SAMPLE_ID"] = sample_name
    mutdensity_df["REGIONS"] = panel_v

    # Save
    mutdensity_df[["SAMPLE_ID", "GENE", "REGIONS", "MUTTYPES",
                "DEPTH",
                "N_MUTS", "N_MUTATED",
                "MUTDENSITY_MB", "MUTDENSITY_MB_ADJUSTED",
                "MUTREADSRATE_MB", "MUTREADSRATE_MB_ADJUSTED",
                ]].to_csv(f"{sample_name}.{panel_v}.mutdensities.tsv",
                                                            sep = "\t",
                                                            header = True,
                                                            index = False
                                                            )


@click.command()
@click.option('--maf_path', type=click.Path(exists=True), required=True, help='Path to the MAF file.')
@click.option('--depths_path', type=click.Path(exists=True), required=True, help='Path to the depths file.')
@click.option('--annot_panel_path', type=click.Path(exists=True), required=True, help='Path to the annotation panel file.')
@click.option('--sample_name', type=str, required=True, help='Sample name.')
@click.option('--panel_version', type=str, required=True, help='Panel version.')
def main(maf_path, depths_path, annot_panel_path, sample_name, panel_version):
    """
    CLI entry point for computing mutation densities.
    """
    compute_mutdensity(maf_path, depths_path, annot_panel_path, sample_name, panel_version)


if __name__ == '__main__':

    main()