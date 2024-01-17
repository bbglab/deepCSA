#!/usr/bin/env python3

import sys
import click
import pandas as pd
import numpy as np
import concurrent.futures

from bgreference import hg38

from utils import contexts_no_change

def to_int_if_possible(string):
    try:
        int(string)
        return True
    except ValueError:
        return False

cb = dict(zip('ACGT', 'TGCA'))

def transform_context_no_query(chr_, pos, trinuc, mut):
    ref, alt = tuple(mut.split('/'))
    ref_triplet = trinuc
    # ref_triplet = hg38(chr_, pos-1, size=3)
    if ref_triplet[1] not in ['C', 'T']:
        ref_triplet = ''.join(list(map(lambda x: cb[x], ref_triplet[::-1])))
        alt = cb[alt]
    return ref_triplet + '>' + alt


def get_non_ref(l, letters = {"A", "C", "G", "T"}):
    return letters - set(l)

def get_trinucl_in_row(x, context_size = 3):
    return hg38(x["CHROM"], x["POS"] - 1, size = 3)

def generate_all_sites_with_context(input_bedfile):

    positions_df = pd.read_csv(input_bedfile, sep = "\t", header = None)

    first_coord = positions_df.iloc[0,1]
    if to_int_if_possible(first_coord):
        positions_df = positions_df.iloc[:,:3]

    # it means there is a header, and we don't want it
    else:
        positions_df = positions_df.iloc[1:,:3]

    positions_df.columns = ["CHROM", "START", "END"]
    positions_df["CHROM"] = positions_df["CHROM"].astype(str).str.replace("chr", "")
    positions_df[["START", "END"]] = positions_df[["START", "END"]].astype(int)


    positions_df["POS"] = [ list(range(x, y+1)) for x, y in positions_df[["START", "END"]].values ]
    positions_df = positions_df.explode("POS").reset_index(drop = True)
    positions_df = positions_df[["CHROM", "POS"]]

    positions_df["CONTEXT"] = positions_df.apply(get_trinucl_in_row, axis = 1)

    # TODO revise if this is a big problem, how is it possible that positions have Ns?
    positions_df["CONTEXT"] = positions_df["CONTEXT"].apply(lambda x : '-' if 'N' in x else x)
    positions_df = positions_df[positions_df["CONTEXT"] != '-'].reset_index(drop = True)

    positions_df["REF"] = positions_df["CONTEXT"].apply( lambda x : x[1])
    positions_df["ALT"] = positions_df["REF"].apply(get_non_ref)

    positions_df = positions_df.explode("ALT").reset_index(drop = True)

    positions_df = positions_df[["CHROM", "POS", "REF", "ALT", "CONTEXT"]]

    positions_df["CONTEXT_MUT"] = positions_df.apply(lambda x: transform_context_no_query(x["CHROM"], x["POS"], x['CONTEXT'], f'{x["REF"]}/{x["ALT"]}'),
                                                                                            axis = 1
                                                                                            )

    return positions_df[["CHROM", "POS", "REF", "ALT", "CONTEXT_MUT"]]




def compute_mutabilities(sample_name, depths_file, mut_profile_file, regions_bedfile, out_mutability, adjusted = True):
    """
    Compute mutational profile from the input data
          ***Remember to add some pseudocounts to the computation***

        Required information:
            Annotated mutations observed
        Output:
            Mutational profile per sample, possibility to add pseudocounts to prevent some probabilities from being 0
    """

    # Load mutation profiles file
    mut_probability = pd.read_csv(mut_profile_file, sep="\t", header=0)
    mut_probability.columns = ["CONTEXT_MUT"] + list(mut_probability.columns[1:])
    mut_prob_samples = list(mut_probability.columns[1:])


    # Load depths file
    coverage_info_all_samples = pd.read_csv(depths_file, sep="\t", header=0)
    coverage_info_all_samples.columns = ["CHROM", "POS"] + list(coverage_info_all_samples.columns[2:])
    coverage_info_all_samples["CHROM"] = coverage_info_all_samples["CHROM"].astype(str)
    coverage_info_all_samples = coverage_info_all_samples.drop("CONTEXT", axis = 1)
    coverage_samples = list(coverage_info_all_samples.columns[2:])

    # Either load or compute a dataframe with all sites annotated
    # all_sites_to_be_included_annotated = generate_all_sites_with_context(regions_bedfile)
    all_sites_to_be_included_annotated = pd.read_csv(regions_bedfile, sep="\t", header=0, usecols = ["CHROM", "POS", "REF", "ALT", "GENE", "CONTEXT_MUT"])

    print("All loaded")

    # make sure that all the files have the chr prefix in the files
    if not str(all_sites_to_be_included_annotated["CHROM"].iloc[0]).startswith("chr"):
        all_sites_to_be_included_annotated["CHROM"] = "chr" + all_sites_to_be_included_annotated["CHROM"].astype(str)

    if not str(coverage_info_all_samples["CHROM"].iloc[0]).startswith("chr"):
        coverage_info_all_samples["CHROM"] = "chr" + coverage_info_all_samples["CHROM"].astype(str)


    if sample_name == "all_samples":
        coverage_info_all_samples["all_samples"] = coverage_info_all_samples[coverage_samples].sum(axis = 1)
        coverage_info_all_samples = coverage_info_all_samples.drop(coverage_samples, axis = 1)

    all_sites_context_n_coverage = all_sites_to_be_included_annotated.merge(coverage_info_all_samples,
                                                                            on = ["CHROM", "POS"],
                                                                            how = 'left')
    all_sites_context_n_coverage = all_sites_context_n_coverage.fillna(0)

    coverage_and_probability = all_sites_context_n_coverage.merge(mut_probability,
                                                                    on = "CONTEXT_MUT",
                                                                    suffixes = ("_coverage", "_probability"),
                                                                    how = 'left')
    print("Coverage and probability merged")
    del all_sites_context_n_coverage
    del coverage_info_all_samples
    del all_sites_to_be_included_annotated
    print("Remove unnecessary objects")

    # if sample_name != "all_samples":

    coverage_and_probability[f"{sample_name}_adjusted_probability"] = coverage_and_probability[f"{sample_name}_probability"] * coverage_and_probability[f"{sample_name}_coverage"]
    # coverage_and_probability[f"{samp}_adjusted_probability"] = coverage_and_probability[f"{samp}_adjusted_probability"] * norm_mut_burden_sample

    list_of_samples_columns = [f"{sample_name}_adjusted_probability" ]
    all_samples_mutability_per_site = coverage_and_probability[ ['CHROM', 'POS', 'REF', 'ALT', 'GENE'] + list_of_samples_columns ].fillna(0)
    all_samples_mutability_per_site = all_samples_mutability_per_site.sort_values(
                                                                by = ['CHROM', 'POS', 'REF', 'ALT']).reset_index(drop = True)
    sample_mutability_info = all_samples_mutability_per_site[['CHROM', 'POS', 'REF', 'ALT', 'GENE', f"{sample_name}_adjusted_probability"]]

    sample_mutability_info_store = sample_mutability_info.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'],
                                                        keep='first').drop("GENE", axis = 1)
    print("Mutabilities computed")

    # print(sample_mutability_info_store[sample_mutability_info_store[f"{sample_name}_adjusted_probability"] > 0].head(50))
    sample_mutability_info_store.to_csv(out_mutability,
                                            sep = "\t",
                                            header = False,
                                            index = False)
    print("mutability 1 stored")

    sample_mutability_info.to_csv(f"{out_mutability}.with_gene",
                                            sep = "\t",
                                            header = True,
                                            index = False)
    print("mutability 1 with gene stored")

    return sample_mutability_info


    # else:

    #     samples = sorted(list(set(mut_prob_samples).intersection(coverage_samples)))
    #     for samp in samples:
    #         # # this should be a single value of the relative mutational burden
    #         # # of this sample
    #         # norm_mut_burden_sample = snv_counts_per_sample[snv_counts_per_sample["SAMPLE_ID"] == samp]["NORM_COUNT"].values[0]

    #         coverage_and_probability[f"{samp}_adjusted_probability"] = coverage_and_probability[f"{samp}_probability"] * coverage_and_probability[f"{samp}_coverage"]
    #         # coverage_and_probability[f"{samp}_adjusted_probability"] = coverage_and_probability[f"{samp}_adjusted_probability"] * norm_mut_burden_sample

    #     list_of_samples_columns = [ x for x in coverage_and_probability.columns if x.endswith("_adjusted_probability") ]
    #     all_samples_mutability_per_site = coverage_and_probability[ ['CHROM', 'POS', 'REF', 'MUT'] + list_of_samples_columns ].fillna(0)
    #     all_samples_mutability_per_site_indexed = all_samples_mutability_per_site.set_index(['CHROM', 'POS', 'REF', 'MUT'])

    #     all_samples_mutability_per_site_all_samples = all_samples_mutability_per_site_indexed.sum(axis = 1).reset_index()
    #     all_samples_mutability_per_site_all_samples.columns = ['CHROM', 'POS', 'REF', 'MUT', 'all_samples_adjusted_probability']
    #     all_samples_mutability_per_site_all_samples = all_samples_mutability_per_site_all_samples.sort_values(
    #                                                                 by = ['CHROM', 'POS', 'REF', 'MUT']).reset_index(drop = True)
    #     all_samples_mutability_per_site_all_samples

    # all_samples_mutability_per_site = all_samples_mutability_per_site.sort_values(
    #                                                             by = ['CHROM', 'POS', 'REF', 'MUT']).reset_index(drop = True)

    # for samp in [x.replace("_adjusted_probability", "") for x in all_samples_mutability_per_site.columns]:
    #     sample_mutability_info = all_samples_mutability_per_site[['CHROM', 'POS', 'REF', 'MUT', f"{samp}_adjusted_probability"]]
    #     sample_mutability_info.to_csv(f"{samp}.mutability_per_site.tsv",
    #                                 sep = "\t",
    #                                 header = False,
    #                                 index = False)

def process_gene_chunk(chunk, sample_name, mutation_matrix_per_sample_gene_int):
    processed_data_chunk = pd.DataFrame()

    for gene, gene_specific in chunk.groupby("GENE"):
        if gene in mutation_matrix_per_sample_gene_int.index:
            correction_factor = mutation_matrix_per_sample_gene_int.loc[gene, "COUNT"]
        else:
            correction_factor = 0

        gene_specific[f"{sample_name}_adjusted_probability"] *= correction_factor
        processed_data_chunk = pd.concat((processed_data_chunk,
                                        gene_specific[['CHROM', 'POS', 'REF', 'ALT', f"{sample_name}_adjusted_probability"]]
                                        ),
                                        axis=0)

    return processed_data_chunk


def parallel_process_chunks(sample_name, mutation_matrix_per_sample_gene_int, sample_mutability_info_int, chunk_size=200):

    # Define the list of genes to process in chunks
    genes_mutability = pd.unique(sample_mutability_info_int["GENE"])
    genes_mutated = pd.unique(mutation_matrix_per_sample_gene_int.index.values)
    common_genes = list(set(genes_mutability) & set(genes_mutated))
    gene_chunks = np.array_split(common_genes, len(common_genes) // chunk_size + 1)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        chunks = [sample_mutability_info_int[sample_mutability_info_int["GENE"].isin(gene_chunk)]
                for gene_chunk in gene_chunks]
        print("chunks mutability created")

        chunks_matrix = [mutation_matrix_per_sample_gene_int[mutation_matrix_per_sample_gene_int.index.isin(gene_chunk)]
                for gene_chunk in gene_chunks]
        print("chunks matrix created")

        print("execution ready to start")
        results = list(executor.map(process_gene_chunk, chunks, [sample_name] * len(chunks), chunks_matrix))
        print("execution finished")

    processed_data = pd.DataFrame()
    for result in results:
        processed_data = pd.concat((processed_data, result), axis=0)

    return processed_data


def adjust_mutabilities(sample_name, mutation_matrix_file, mutability_info, out_mutability_adjusted):

    # Load mutation profiles file
    mutation_matrix = pd.read_csv(mutation_matrix_file, sep="\t", header=0)
    print("mutation matrix loaded")

    mutation_matrix_per_sample_gene = mutation_matrix.groupby(by = ["SYMBOL"])["MUT_ID"].count().reset_index()
    print("group by worked")
    mutation_matrix_per_sample_gene.columns = ["SYMBOL", "COUNT"]
    mutation_matrix_per_sample_gene = mutation_matrix_per_sample_gene.fillna(0)
    mutation_matrix_per_sample_gene = mutation_matrix_per_sample_gene.set_index("SYMBOL")
    mutation_matrix_per_sample_gene.index.name = "SYMBOL"
    mutation_matrix_per_sample_gene = mutation_matrix_per_sample_gene / mutation_matrix_per_sample_gene.min()
    print("normalization with minimum worked")


    # Parallel processing
    corrected_data = parallel_process_chunks(sample_name, mutation_matrix_per_sample_gene, mutability_info)


    print("Corrected mutabilities")
    del mutability_info
    print("Removing second table")
    corrected_data = corrected_data.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT', f"{sample_name}_adjusted_probability"],
                                                    ascending = (True, True, True, True, False)
                                                    )
    print("Corrected data sorted")

    corrected_data = corrected_data.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'], keep='first')
    print("Corrected data duplicates removed")

    corrected_data.to_csv(out_mutability_adjusted,
                            sep = "\t",
                            header = False,
                            index = False)



@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--mutation_matrix', type=click.Path(exists=True), help='Input mutations matrix file')
@click.option('--depths', type=click.Path(exists=True), help='Input depths file')
@click.option('--profile', type=click.Path(exists=True), help='Input mutational profile dataframe file.')
@click.option('--bedfile', type=click.Path(exists=True), help='BED file of the regions.')
@click.option('--out_mutability', type=click.Path(), help='Output mutability file.')
@click.option('--adjust_local_rate', is_flag=True, help='Generate an additional file with the mutabilities adjusted by the local mutation rate of each gene.')
# @click.option('--json_filters', type=click.Path(exists=True), help='Input mutation filtering criteria file')
# @click.option('--method', type=click.Choice(['unique', 'multiple']), default='unique')
# @click.option('--pseud', type=float, default=0.5)
# @click.option('--mutation_matrix', type=click.Path(exists=True), help='Mutation matrix file (for profile mode)')
# @click.option('--trinucleotide_counts', type=click.Path(exists=True), help='Trinucleotide counts file (for profile mode)')
# @click.option('--plot', is_flag=True, help='Generate plot and save as PDF')


# def main(mode, sample_name, mut_file, out_matrix, json_filters, method, pseud, mutation_matrix, trinucleotide_counts, out_profile, plot):
def main(sample_name, mutation_matrix, depths, profile, bedfile, out_mutability, adjust_local_rate):
    click.echo(f"Computing the mutabilities...")
    # click.echo(f"Using the pseudocount: {pseud}")
    mutabilities_with_gene = compute_mutabilities(sample_name, depths, profile, bedfile, out_mutability)
    click.echo("Mutabilities computed.")
    if adjust_local_rate:
        adjust_mutabilities(sample_name, mutation_matrix, mutabilities_with_gene, f"{out_mutability}.adjusted")

if __name__ == '__main__':
    main()

