#!/usr/bin/env python3


import sys
import click
import json
import pandas as pd

from utils import contexts_formatted
from utils import filter_maf


def compute_mutation_matrix(sample_name, mutations_file, mutation_matrix, json_filters, method = 'unique', pseudocount = 0):
    """
    Compute mutational profile from the input data
          ***Remember to add some pseudocounts to the computation***

        Required information:
            Annotated mutations observed
        Output:
            Mutational profile per sample, possibility to add pseudocounts to prevent some probabilities from being 0
    """

    if method not in ['unique', 'multiple']:
        print('a')
        exit(1)

    # Load your MAF DataFrame (raw_annotated_maf)
    raw_annotated_maf = pd.read_csv(mutations_file, sep = "\t", header = 0)

    # Load the filter criteria from the JSON file
    with open(json_filters, 'r') as file:
        filter_criteria = json.load(file)

    if len(filter_criteria) > 0:
        # Filter the annotated maf using the described filters
        annotated_maf = filter_maf(raw_annotated_maf, filter_criteria)
    else:
        annotated_maf = raw_annotated_maf



    # create the matrix in the desired order
    empty_matrix = pd.DataFrame(index = contexts_formatted)

    # if method == 'unique':
    #     # make sure to count each mutation only once (avoid annotation issues)
    #     annotated_minimal_maf = annotated_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID"]].drop_duplicates().reset_index(drop = True)
    #     annotated_minimal_maf.columns = ["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID"]

    #     # count the mutations per sample and per context
    #     counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["MUT_ID"].count().reset_index()
    #     counts_x_sample_matrix = counts_x_sample_context_long.pivot(index = "CONTEXT_MUT", columns = "SAMPLE_ID", values = "MUT_ID")


    # elif method == 'multiple':
    #     # make sure to count each mutation only once (avoid annotation issues)
    #     annotated_minimal_maf = annotated_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "ALT_DEPTH"]].drop_duplicates().reset_index(drop = True)
    #     annotated_minimal_maf.columns = ["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "ALT_DEPTH"]

    #     # count the mutations per sample and per context
    #     counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["ALT_DEPTH"].sum().reset_index()
    #     counts_x_sample_matrix = counts_x_sample_context_long.pivot(index = "CONTEXT_MUT", columns = "SAMPLE_ID", values = "ALT_DEPTH")
    if method == 'unique':
        # make sure to count each mutation only once (avoid annotation issues)
        annotated_minimal_maf = annotated_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID"]].drop_duplicates().reset_index(drop = True)
        annotated_minimal_maf.columns = ["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID"]

        # count the mutations per sample and per context
        counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["CONTEXT_MUT"])["MUT_ID"].count().reset_index()
        counts_x_sample_matrix = counts_x_sample_context_long.set_index("CONTEXT_MUT")
        counts_x_sample_matrix.columns = [sample_name]


    elif method == 'multiple':
        # make sure to count each mutation only once (avoid annotation issues)
        annotated_minimal_maf = annotated_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "ALT_DEPTH"]].drop_duplicates().reset_index(drop = True)
        annotated_minimal_maf.columns = ["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "ALT_DEPTH"]

        # count the mutations per sample and per context
        counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["CONTEXT_MUT"])["ALT_DEPTH"].sum().reset_index()
        counts_x_sample_matrix = counts_x_sample_context_long.set_index("CONTEXT_MUT")
        counts_x_sample_matrix.columns = [sample_name]

    counts_x_sample_matrix = pd.concat( (empty_matrix, counts_x_sample_matrix) , axis = 1)
    counts_x_sample_matrix = counts_x_sample_matrix.fillna(0)
    counts_x_sample_matrix.index.name = "CONTEXT_MUT"

    # add a pseudocount if desired
    counts_x_sample_matrix = counts_x_sample_matrix + pseudocount

    counts_x_sample_matrix.to_csv(mutation_matrix,
                                    header = True,
                                    index = True,
                                    sep = "\t")


def compute_mutation_profile(sample_name, mutation_matrix_file, trinucleotide_counts_file, json_output, plot):
    """
    Compute mutational profile from the input data

        Required information:
            Mutation matrix
            Trinucleotide content of the sequenced region (depth-aware or non-depth-aware)
        Output:
            Mutational profile per sample
    """

    # Load your mutation matrix
    mutation_matrix = pd.read_csv(mutation_matrix_file, sep = "\t", header = 0)
    mutation_matrix = mutation_matrix.set_index("CONTEXT_MUT")
    # sample_name = mutation_matrix.columns[0]

    trinucleotide_counts = pd.read_csv(trinucleotide_counts_file, sep = "\t", header = None)
    trinucleotide_counts.columns = ["CONTEXT", sample_name]
    trinucleotide_counts = trinucleotide_counts.set_index("CONTEXT")

    trinuc_per_contextmuts = []
    for x in contexts_formatted:
        trinuc_per_contextmuts.append(trinucleotide_counts.loc[x[:3], sample_name])

    contextmut_depth = pd.DataFrame({sample_name : trinuc_per_contextmuts, "CONTEXT_MUT": contexts_formatted})
    contextmut_depth = contextmut_depth.set_index("CONTEXT_MUT")

    # divide
    mut_probability = mutation_matrix.divide( contextmut_depth )

    # normalize
    mut_probability = mut_probability / mut_probability.sum()

    # reindex to ensure the right order
    mut_probability = mut_probability.reindex(contexts_formatted)
    mut_probability.index.name = "CONTEXT_MUT"
    mut_probability = mut_probability.reset_index()

    # # Load the filter criteria from the JSON file
    # with open("sample.json", "w") as outfile:
    #     json.dump(dictionary, outfile)

    mut_probability.to_csv(json_output,
                            header = True,
                            index = False,
                            sep = "\t")

    if plot:
        from utils_plot import plot_profile

        max_freq = max(mut_probability[sample_name]) * 1.1

        order_mag = 100
        size_step = 10
        ylabs = []
        while len(ylabs) < 2:
            ylabs = [ x/order_mag for x in range(0, round(max_freq * order_mag) + 1 ) if (x % size_step) == 0 ]
            if size_step == 10:
                size_step = 5
            else:
                order_mag *= 10
                size_step = 10

        plot_profile(dict(zip(mut_probability["CONTEXT_MUT"], mut_probability[sample_name])),
                        title=f'{sample_name}',
                        ylabels = ylabs,
                        ymax = max_freq,
                        output_f = f'{sample_name}.profile.pdf')



@click.command()
@click.argument('mode', type=click.Choice(['matrix', 'profile']))
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--mut_file', type=click.Path(exists=True), help='Input mutation file')
@click.option('--out_matrix', type=click.Path(), help='Output mutation matrix file')
@click.option('--json_filters', type=click.Path(exists=True), help='Input mutation filtering criteria file')
@click.option('--method', type=click.Choice(['unique', 'multiple']), default='unique')
@click.option('--pseud', type=float, default=0.5)

@click.option('--mutation_matrix', type=click.Path(exists=True), help='Mutation matrix file (for profile mode)')
@click.option('--trinucleotide_counts', type=click.Path(exists=True), help='Trinucleotide counts file (for profile mode)')
@click.option('--out_profile', type=click.Path(), help='JSON output file (for profile mode)')
@click.option('--plot', is_flag=True, help='Generate plot and save as PDF')

def main(mode, sample_name, mut_file, out_matrix, json_filters, method, pseud, mutation_matrix, trinucleotide_counts, out_profile, plot):

    if mode == 'matrix':
        click.echo(f"Running in matrix mode...")
        click.echo(f"Using the method: {method}")
        click.echo(f"Using the pseudocount: {pseud}")
        compute_mutation_matrix(sample_name, mut_file, out_matrix, json_filters, method, pseud)

    elif mode == 'profile':
        click.echo(f"Running in profile mode...")
        compute_mutation_profile(sample_name, mutation_matrix, trinucleotide_counts, out_profile, plot)
        click.echo("Profile computation completed.")

    else:
        click.echo("Unrecognized mode. Choose 'matrix' or 'profile'.")

if __name__ == '__main__':
    main()

