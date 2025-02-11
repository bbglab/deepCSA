#!/usr/local/bin/python


import click
import pandas as pd


def process_signatures(signature_probabilities_files):
    """
    INFO
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
        sig_probs_matrix[sig_probs_matrix["Sample Names"] == sample].drop(["Sample Names"], axis = 1).to_csv(f"{sample}.decomposed_probabilities.tsv",
                                                                                                header=True,
                                                                                                index=False,
                                                                                                sep="\t")



@click.command()
@click.option('--signature-probabilities', type=click.Path(exists=True), help='File listing decomposed mutation probability files.')

def main(signature_probabilities):
    click.echo("Combining signature probabilities of all samples and groups...")
    process_signatures(signature_probabilities)

if __name__ == '__main__':
    main()

