#!/usr/bin/env python

# TODO: bump pandas to 2.2.3

import sys
import pandas as pd
from read_utils import custom_na_values


# -- Main function -- #
def compute_mutepithelium_exon(sample_n, mutations_file, denominator_file):
    '''
    Function that compute the fraction of mutated epithelium per exon given
    the number of reads and unique genomes in an exon and a dataframe of mutations that is subsetted to count mutations in each individual exon.
    '''

    # File loading
    denominator_data = pd.read_csv(denominator_file, sep = "\t", header = 0, na_values = custom_na_values)
    mutations = pd.read_csv(mutations_file, sep = "\t", header = 0, na_values = custom_na_values)

    mutations_list = []
    mutated_reads_list = []
    for ind, row in denominator_data.iterrows():
        numerator = mutations[ (mutations["CHROM"] == row["CHROM"]) & (mutations["POS"] >= row["START"]) & (mutations["POS"] <= row["END"]) ]
        mutations_list.append(numerator.shape[0])
        mutated_reads_list.append(numerator["ALT_DEPTH"].sum())

    mutated_epithelium_data = denominator_data.copy()
    mutated_epithelium_data["MUTATIONS"] = mutations_list
    mutated_epithelium_data["MUTATED_READS"] = mutated_reads_list
    mutated_epithelium_data["LB_MUTATED_EPITHELIUM"] = mutated_epithelium_data["MUTATED_READS"] / mutated_epithelium_data["UNIQ_READS"]
    mutated_epithelium_data["UB_MUTATED_EPITHELIUM"] = mutated_epithelium_data["MUTATED_READS"] / mutated_epithelium_data["UNIQ_GENOMES"]


    orig_columns = list(mutated_epithelium_data.columns)
    mutated_epithelium_data["SAMPLE_ID"] = sample_n# ame.split('.')[0]
    mutated_epithelium_data = mutated_epithelium_data[['SAMPLE_ID'] + orig_columns]

    # mutated_epithelium_data["REGIONS"] = panel_v

    # CHROM	START	END	GENE	EXON	UNIQ_GENOMES	UNIQ_READS	MUTATIONS	MUTATED_READS	LB_MUTATED_EPITHELIUM	UB_MUTATED_EPITHELIUM
    # # Save
    # mutrate_df[["SAMPLE_ID", "GENE", "REGIONS",
    #             "N_MUTS", "N_MUTATED", "DEPTH",
    #             "MUTRATE_MB", "MUTREADSRATE_MB"]].to_csv(f"{sample_name}.{panel_v}.mutrates.tsv",
    #                                                         sep = "\t",
    #                                                         header = True,
    #                                                         index = False
    #                                                         )

    mutated_epithelium_data.to_csv(f"{sample_n}.exon.mutated_epithelium.tsv",
                                                                sep = "\t",
                                                                header = True,
                                                                index = False
                                                                )

    return mutated_epithelium_data

def inclusion_exclusion_init(plist_init):
    plist_init = list(plist_init)
    filtered_plist = [x for x in plist_init if x != 0]

    return inclusion_exclusion(filtered_plist)



def inclusion_exclusion(plist):
    """A more efficient version of the algorithm"""

    plist = list(plist)

    n = len(plist)
    if n > 2:
        p1 = inclusion_exclusion(plist[: n // 2])
        p2 = inclusion_exclusion(plist[n // 2: ])
        return inclusion_exclusion([p1, p2])
    if n == 2:
        return plist[0] + plist[1] - plist[0] * plist[1]
    if n == 1:
        return plist[0]



def compute_mut_epi_per_gene(sample_n, mut_epi_per_exon):
    gene_grouped = mut_epi_per_exon.groupby(by = "GENE")[["LB_MUTATED_EPITHELIUM", "UB_MUTATED_EPITHELIUM",
                                                            "UNIQ_GENOMES", "UNIQ_READS"]].sum().reset_index()
    gene_grouped.columns = ["GENE", "LB_MUTATED_EPITHELIUM_SUM", "UB_MUTATED_EPITHELIUM_SUM", "UNIQ_GENOMES_SUM", "UNIQ_READS_SUM"]

    gene_grouped_pie = mut_epi_per_exon.groupby(by = "GENE").agg({"LB_MUTATED_EPITHELIUM" : inclusion_exclusion_init , "UB_MUTATED_EPITHELIUM" : inclusion_exclusion_init,
                                                                    "UNIQ_GENOMES" : max, "UNIQ_READS" : max }).reset_index()
    gene_grouped_pie.columns = ["GENE", "LB_MUTATED_EPITHELIUM_PIE", "UB_MUTATED_EPITHELIUM_PIE", "UNIQ_GENOMES_MAX", "UNIQ_READS_MAX"]

    gene_mutated_epithelium = gene_grouped_pie.merge(gene_grouped, on = ["GENE"])

    orig_columns = list(gene_mutated_epithelium.columns)
    gene_mutated_epithelium["SAMPLE_ID"] = sample_n
    gene_mutated_epithelium = gene_mutated_epithelium[['SAMPLE_ID'] + orig_columns]

    gene_mutated_epithelium.to_csv(f"{sample_n}.gene.mutated_epithelium.tsv",
                                                                sep = "\t",
                                                                header = True,
                                                                index = False
                                                                )

    return gene_mutated_epithelium




def compute_mut_epi_per_sample(sample_n, mut_epi_per_gene):

    sample_metrics_dict = {}

    gene_metrics = ["LB_MUTATED_EPITHELIUM_PIE", "UB_MUTATED_EPITHELIUM_PIE",
                        "LB_MUTATED_EPITHELIUM_SUM", "UB_MUTATED_EPITHELIUM_SUM"]
    gene_reads_metrics = ["UNIQ_GENOMES_MAX", "UNIQ_READS_MAX",
                            "UNIQ_GENOMES_SUM", "UNIQ_READS_SUM"]

    method = 'PIE'
    for gene_metric in gene_metrics:
        sample_metrics_dict[f'{gene_metric}_G_{method}'] = inclusion_exclusion_init(mut_epi_per_gene[gene_metric])

    method = 'SUM'
    for gene_metric in gene_metrics:
        sample_metrics_dict[f'{gene_metric}_G_{method}'] = mut_epi_per_gene[gene_metric].sum()

    method = 'MAX'
    for gene_reads in gene_reads_metrics:
        sample_metrics_dict[f'{gene_reads}_G_{method}'] = mut_epi_per_gene[gene_reads].max()

    method = 'SUM'
    for gene_reads in gene_reads_metrics:
        sample_metrics_dict[f'{gene_reads}_G_{method}'] = mut_epi_per_gene[gene_reads].sum()

    output_df = pd.DataFrame.from_dict(sample_metrics_dict, orient = 'index').T

    orig_columns = list(output_df.columns)
    output_df = output_df.reset_index()
    output_df["index"] = sample_n
    output_df.columns = ['SAMPLE_ID'] + orig_columns

    output_df.to_csv(f"{sample_n}.sample.mutated_epithelium.tsv",
                                    sep = "\t",
                                    header = True,
                                    index = False
                                    )

    return output_df


if __name__ == '__main__':
    # TODO reimplement with click
    sample_name = sys.argv[1]
    mutations_file = sys.argv[2]
    read_counts_file = sys.argv[3]


    mut_epi_exon = compute_mutepithelium_exon(sample_name, mutations_file, read_counts_file)# , sample_name, panel_version)
    mut_epi_gene = compute_mut_epi_per_gene(sample_name, mut_epi_exon)
    mut_epi_sample = compute_mut_epi_per_sample(sample_name, mut_epi_gene)



