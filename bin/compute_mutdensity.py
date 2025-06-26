#!/usr/bin/env python

"""
Mutation density computation script.
Mutation density is a metric that quantifies the number of mutations per megabase (Mb) of sequenced DNA.
This script computes mutation density for a given sample, all possible genes and a list of consequence type groups.
It calculates mutation density per Mb, both adjusted and non-adjusted by the number of sites available for each consequence type.
The results are saved to a TSV file.
"""


import click
import pandas as pd
from read_utils import custom_na_values

# TODO: bump pandas to 2.2.3

# -- Auxiliary functions -- #

MUTDENSITY_IMPACT_GROUPS = [False, ["SNV"] , ["INSERTION", "DELETION"], ["SNV", "INSERTION", "DELETION"]]

def mutdensity_sample(maf_df, depths_df, depths_adj_df, sample_name):
    """
    Computes a sample's global mutation density. Returns the mutation density
    per Mb, non-adjusted and adjusted by panel
    composition.
    """

    impact_group_results = list()

    # mutation density depth information
    sample_features_depth = {"DEPTH" : depths_df.drop_duplicates(subset = ["CHROM", "POS"])[f"{sample_name}"].sum(),
                                "DEPTH_ADJUSTED": depths_adj_df[f"{sample_name}"].sum()
                                }

    for type_list in MUTDENSITY_IMPACT_GROUPS:
        if not type_list:
            unique_maf = maf_df[["SAMPLE_ID", "MUT_ID", "ALT_DEPTH"]].drop_duplicates()
            types_included = 'all_types'
        else:
            unique_maf = maf_df[maf_df['TYPE'].isin(type_list)][["SAMPLE_ID", "MUT_ID", "ALT_DEPTH"]].copy().drop_duplicates()
            types_included = '-'.join(sorted(type_list))

        # count number of mutations and mutated reads in the sample
        ## make sure to count each mutation only once (avoid annotation issues)
        n_muts = unique_maf.shape[0]
        n_muts_per_sample = unique_maf.groupby(by = ["SAMPLE_ID", "MUT_ID"] ).agg({"ALT_DEPTH" : "sum" }).reset_index()
        n_mutated_reads = n_muts_per_sample["ALT_DEPTH"].sum()
        print(n_muts, n_mutated_reads)

        # mutation density metrics
        sample_features = dict().update(sample_features_depth)
        sample_features["N_MUTS"] = n_muts
        sample_features["N_MUTATED"] = n_mutated_reads
        
        sample_features["MUTDENSITY_MB"] = ( sample_features["N_MUTS"] / sample_features["DEPTH"] * 1000000 ).astype(float)
        sample_features["MUTDENSITY_MB_ADJUSTED"] = ( sample_features["N_MUTS"] / sample_features["DEPTH_ADJUSTED"] * 1000000 ).astype(float)
        sample_features["MUTREADSRATE_MB"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH"] * 1000000 ).astype(float)
        sample_features["MUTREADSRATE_MB_ADJUSTED"] = ( sample_features["N_MUTATED"] / sample_features["DEPTH_ADJUSTED"] * 1000000 ).astype(float)

        sample_features["GENE"] = "ALL_GENES"
        sample_features["MUTTYPES"] = types_included

        impact_group_results.append(pd.DataFrame([sample_features]))

    # concatenate results for all impact groups
    mutdensity_sample = pd.concat(impact_group_results)

    return mutdensity_sample


def mutdensity_gene(maf_df, depths_df, depths_adj_df, sample_name):
    """
    Computes each gene mutation density. Returns the mutation density
    both per Mb and Kb sequenced, both non-adjusted and adjusted by panel
    composition.
    """

    impact_group_results = list()

    for type_list in MUTDENSITY_IMPACT_GROUPS:
        # filter by mutation type according to type_list
        if not type_list:
            unique_maf = maf_df[["SAMPLE_ID", "GENE", "MUT_ID", "ALT_DEPTH"]].drop_duplicates()
            types_included = 'all_types'
        else:
            unique_maf = maf_df[maf_df['TYPE'].isin(type_list)][["SAMPLE_ID", "GENE", "MUT_ID", "ALT_DEPTH"]].copy().drop_duplicates()
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
        impact_group_results.append(mut_depths_df.reset_index())

    # concatenate results for all impact groups
    mutdensity_per_gene = pd.concat(impact_group_results)

    return mutdensity_per_gene


def load_n_process_inputs(maf_path, depths_path, annot_panel_path, sample_name):
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

    return maf_df, depths_subset_df, depths_subset_adj_df, depths_subset_adj_sample_df



# -- Main function -- #
def compute_mutdensity(maf_path, depths_path, annot_panel_path, sample_name, panel_v):
    """
    Computes mutation density for a given sample based on MAF, depths, and annotation panel files.
    The function calculates mutation density per Mb and Kb, both adjusted and non-adjusted by
    the panel composition. It saves the results to a TSV file.
    """

    maf_df, depths_subset_df, depths_subset_adj_df, depths_subset_adj_sample_df = load_n_process_inputs(maf_path, depths_path, annot_panel_path, sample_name)

    # Compute mutation densities
    ## sample mutation density
    mutdensity_sample_df = mutdensity_sample(maf_df, depths_subset_df, depths_subset_adj_sample_df, sample_name)

    ## per gene mutation density
    mutdensity_genes_df = mutdensity_gene(maf_df, depths_subset_df, depths_subset_adj_df, sample_name)

    mutdensity_df = pd.concat([mutdensity_sample_df, mutdensity_genes_df])

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