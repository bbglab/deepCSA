#!/usr/bin/env python3


import sys
import click
import json
import pandas as pd

from utils import contexts_formatted
from utils import filter_maf


def compute_mutation_matrix(mutations_file, mutation_matrix, json_filters, method = 'unique', pseudocount = 0):
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

    if method == 'unique':
        # make sure to count each mutation only once (avoid annotation issues)
        annotated_minimal_maf = annotated_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID"]].drop_duplicates().reset_index(drop = True)
        annotated_minimal_maf.columns = ["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID"]

        # count the mutations per sample and per context
        counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["MUT_ID"].count().reset_index()
        counts_x_sample_matrix = counts_x_sample_context_long.pivot(index = "CONTEXT_MUT", columns = "SAMPLE_ID", values = "MUT_ID")


    elif method == 'multiple':
        # make sure to count each mutation only once (avoid annotation issues)
        annotated_minimal_maf = annotated_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "ALT_DEPTH"]].drop_duplicates().reset_index(drop = True)
        annotated_minimal_maf.columns = ["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "ALT_DEPTH"]

        # count the mutations per sample and per context
        counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["ALT_DEPTH"].sum().reset_index()
        counts_x_sample_matrix = counts_x_sample_context_long.pivot(index = "CONTEXT_MUT", columns = "SAMPLE_ID", values = "ALT_DEPTH")

    counts_x_sample_matrix = pd.concat( (empty_matrix, counts_x_sample_matrix) , axis = 1)
    counts_x_sample_matrix = counts_x_sample_matrix.fillna(0)
    counts_x_sample_matrix.index.name = "CONTEXT_MUT"

    # add a pseudocount if desired
    counts_x_sample_matrix = counts_x_sample_matrix + pseudocount

    counts_x_sample_matrix.to_csv(mutation_matrix,
                                    header = True,
                                    index = True,
                                    sep = "\t")


@click.command()
@click.argument('mode', type=click.Choice(['matrix']), )
@click.option('--mut_file', type=click.Path(exists=True), help='Input mutation file')
@click.option('--out_matrix', type=click.Path(), help='Output mutation matrix file')
@click.option('--json_filters', type=click.Path(exists=True), help='Input mutation filtering criteria file')
@click.option('--method', type=click.Choice(['unique', 'multiple']), default='unique')
@click.option('--pseud', type=float, default=0.5)

def main(mode, mut_file, out_matrix, json_filters, method, pseud):
    if mode == 'matrix':
        click.echo(f"Using the method: {method}")
        click.echo(f"Using the pseudocount: {pseud}")
        compute_mutation_matrix(mut_file, out_matrix, json_filters, method, pseud)
    else:
        click.echo("Unrecognized mode. Choose 'matrix'.")

if __name__ == '__main__':
    main()

