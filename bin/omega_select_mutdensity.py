#!/usr/bin/env python


import click
import pandas as pd
from read_utils import custom_na_values


def select_syn_mutdensity(mutdensity_file, output_file, mode):
    """
    INFO
    """

    mutdensity_df = pd.read_csv(mutdensity_file, sep = "\t", header = 0, na_values = custom_na_values)

    synonymous_mutdensities_all_samples = mutdensity_df[(mutdensity_df["MUTTYPES"] == "SNV") &
                                                        (mutdensity_df["GENE"] != "ALL_GENES") &
                                                        ~(mutdensity_df["GENE"].str.contains("--"))].reset_index(drop = True)

    if mode == 'mutations':
        synonymous_mutdensities_genes = synonymous_mutdensities_all_samples[['GENE', 'MUTDENSITY_MB_ADJUSTED']]
    elif mode == 'mutated_reads':
        synonymous_mutdensities_genes = synonymous_mutdensities_all_samples[['GENE', 'MUTREADSRATE_MB_ADJUSTED']]

    synonymous_mutdensities_genes.columns = ["GENE", "MUTDENSITY"]
    synonymous_mutdensities_genes.to_csv(f"{output_file}",
                                        header=True,
                                        index=False,
                                        sep="\t")


@click.command()
@click.option('--mutdensities', type=click.Path(exists=True), help='Input mutation density file')
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--mode', type=click.Choice(['mutations', 'mutated_reads']), default='mutations')
def main(mutdensities, output, mode):
    click.echo("Selecting the gene synonymous mutation densities...")
    select_syn_mutdensity(mutdensities, output, mode)

if __name__ == '__main__':
    main()

