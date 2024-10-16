#!/usr/local/bin/python


import click
import pandas as pd
from read_utils import custom_na_values


def compute_relative_mutrate(mutrate_file, output_file):
    """
    INFO
    """

    mutrate_df = pd.read_csv(mutrate_file, sep = "\t", header = 0, na_values = custom_na_values)

    synonymous_mutrates_all_samples = mutrate_df[(mutrate_df["MUTTYPES"] == "SNV")].reset_index(drop = True)

    relative_synonymous_mutrates_all_samples = synonymous_mutrates_all_samples[['GENE', 'N_MUTS']].set_index(['GENE']) \
                                                    / synonymous_mutrates_all_samples["N_MUTS"].sum()
    relative_synonymous_mutrates_all_samples.columns = ["REL_MUTRATE"]

    relative_synonymous_mutrates_all_samples.reset_index()[["GENE", "REL_MUTRATE"]].to_csv(f"{output_file}",
                                                            header=["GENE", "SYNONYMOUS_MUTS"],
                                                            index=False,
                                                            sep="\t")


@click.command()
@click.option('--mutations', type=click.Path(exists=True), help='Input muttations file')
@click.option('--signature-probabilities', type=click.Path(exists=True), help='File listing decomposed mutation probability files.')
@click.option('--output', type=click.Path(), help='Output annotated mutations file')
# @click.option('--all-groups', type=click.Path(exists=True), help='JSON groups file')


def main(mutrates, output):
    click.echo(f"Computing the relative mutation rate...")
    compute_relative_mutrate(mutrates, output)

if __name__ == '__main__':
    main()

