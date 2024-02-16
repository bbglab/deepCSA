#!/usr/local/bin/python


import click
import pandas as pd


def combine_mutations_n_signatures(mutations_file, signature_probabilities_files, output_file):
    """
    INFO
    """

    mutations = pd.read_csv(mutations_file, sep = "\t", header = 0)

    sig_probs_matrix = pd.concat((pd.read_csv(file.strip(), sep='\t', header=0) for file in open(signature_probabilities_files, 'r')), axis=0)
    ## this is equivalent to the line above
    # sig_probs_matrix = pd.DataFrame()
    # with open(signature_probabilities_files, 'r') as file:
    #     for line in file.readlines():
    #         sig_file = line.strip()
    #         sig_probs_matrix_single = pd.read_csv(sig_file, sep = "\t", header = 0)
    #         sig_probs_matrix = pd.concat((sig_probs_matrix, sig_probs_matrix_single), axis = 0)
    sig_probs_matrix = sig_probs_matrix.fillna(0)

    sig_probs_matrix.to_csv("test_signature_probs.tsv",
                                header=True,
                                index=False,
                                sep="\t")

    muts_with_sbs = mutations.merge(sig_probs_matrix,
                                    left_on = ["SAMPLE_ID", "CONTEXT_MUT_SIGPRO"],
                                    right_on = ["Sample Names", "MutationType"],
                                    how = 'left')

    muts_with_sbs.drop(["Sample Names", "MutationType"], axis = 1).to_csv(f"{output_file}",
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

