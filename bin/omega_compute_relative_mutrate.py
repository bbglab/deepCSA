#!/usr/local/bin/python


import click
import pandas as pd
from read_utils import custom_na_values


def compute_relative_mutrate(mutrate_file, output_file):
    """
    INFO
    """

    mutrate_df = pd.read_csv(mutrate_file, sep = "\t", header = 0, na_values = custom_na_values)

    synonymous_mutrates_all_samples = mutrate_df[(mutrate_df["MUTTYPES"] == "SNV") & 
                                                  (mutrate_df["GENE"] != "ALL_GENES")].reset_index(drop = True)

    relative_synonymous_mutrates_all_samples = synonymous_mutrates_all_samples[['GENE', 'MUTRATE_MB_ADJUSTED']].set_index(['GENE']) \
                                                    / synonymous_mutrates_all_samples["MUTRATE_MB_ADJUSTED"].sum()
    relative_synonymous_mutrates_all_samples.columns = ["REL_MUTRATE"]

    relative_synonymous_mutrates_all_samples.reset_index()[["GENE", "REL_MUTRATE"]].to_csv(f"{output_file}",
                                                            header=["GENE", "SYNONYMOUS_MUTS"],
                                                            index=False,
                                                            sep="\t")


@click.command()
@click.option('--mutrates', type=click.Path(exists=True), help='Input mutation rate file')
@click.option('--output', type=click.Path(), help='Output file')


def main(mutrates, output):
    click.echo(f"Computing the relative mutation rate...")
    compute_relative_mutrate(mutrates, output)

if __name__ == '__main__':
    main()

