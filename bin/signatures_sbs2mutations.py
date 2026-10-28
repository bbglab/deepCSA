#!/usr/bin/env python

# TODO: bump pandas to 2.2.3

import click
import pandas as pd


def process_signatures(signature_probabilities_files):
    """
    Processes mutation probability signature files and outputs per-sample probability files.
    """

    sig_probs_matrix = pd.concat((pd.read_csv(file.strip(), sep='\t', header=0) for file in open(signature_probabilities_files, 'r')), axis=0)
    ## this is equivalent to the line above
    # sig_probs_matrix = pd.DataFrame()
    # with open(signature_probabilities_files, 'r') as file:
    #     for line in file.readlines():
    #         sig_file = line.strip()
    #         sig_probs_matrix_single = pd.read_csv(sig_file, sep = "\t", header = 0)
    #         sig_probs_matrix = pd.concat((sig_probs_matrix, sig_probs_matrix_single), axis = 0)
    sig_probs_matrix = sig_probs_matrix.fillna(0)

    for sample in sig_probs_matrix["Sample Names"].unique():
        sample_df = sig_probs_matrix[sig_probs_matrix["Sample Names"] == sample].drop(["Sample Names"], axis = 1).copy()

        sample_df.to_csv(f"{sample}.decomposed_probabilities.tsv",
                            header=True,
                            index=False,
                            float_format=None,
                            sep="\t")



@click.command()
@click.option('--signature-probabilities', type=click.Path(exists=True), help='File listing decomposed mutation probability files.')

def main(signature_probabilities):
    click.echo("Combining signature probabilities of all samples and groups...")
    process_signatures(signature_probabilities)

if __name__ == '__main__':
    main()
