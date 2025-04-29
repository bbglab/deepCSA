#!/usr/bin/env python


import click
import pandas as pd
from read_utils import custom_na_values


def combine_mutations_n_signatures(mutations_file, signature_probabilities_files, output_file):
    """
    INFO
    """

    mutations = pd.read_csv(mutations_file, sep = "\t", header = 0, na_values = custom_na_values)

    sig_probs_matrix = pd.read_csv(signature_probabilities_files, sep = "\t", header = 0)

    muts_with_sbs = mutations.merge(sig_probs_matrix,
                                    left_on = ["CONTEXT_MUT_SIGPRO"],
                                    right_on = ["MutationType"],
                                    how = 'left')

    muts_with_sbs.drop(["MutationType"], axis = 1).to_csv(f"{output_file}",
                                                            header=True,
                                                            index=False,
                                                            sep="\t")




@click.command()
@click.option('--mutations', type=click.Path(exists=True), help='Input muttations file')
@click.option('--signature-probabilities', type=click.Path(exists=True), help='File listing decomposed mutation probability files.')
@click.option('--output', type=click.Path(), help='Output annotated mutations file')
# @click.option('--all-groups', type=click.Path(exists=True), help='JSON groups file')


def main(mutations, signature_probabilities, output):
    click.echo(f"Combining mutations MAF with signature probabilities information...")
    combine_mutations_n_signatures(mutations, signature_probabilities, output)

if __name__ == '__main__':
    main()

