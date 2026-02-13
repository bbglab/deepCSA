#!/usr/bin/env python


import click
import pandas as pd
from read_utils import custom_na_values


def select_syn_mutdensity(mutdensity_file, output_file, mode):
    """
    This function selects the synonymous mutation densities for all genes
    from the mutation density file of all samples.
    
    right now the use of mode is not implemented,
    since we only compute one type of synonymous mutation densities.
    """

    mutdensity_df = pd.read_csv(mutdensity_file, sep = "\t", header = 0, na_values = custom_na_values)

    synonymous_mutdensities_all_samples = mutdensity_df[(mutdensity_df["SAMPLE_ID"] == 'all_samples') &
                                                        ~(mutdensity_df["GENE"].str.contains("--"))
                                                        ]["synonymous"].reset_index(drop = True)

    synonymous_mutdensities_genes = synonymous_mutdensities_all_samples[['GENE', 'synonymous']]

    # TODO implement these different modes if appropriate
    # if mode == 'mutations':
    #     synonymous_mutdensities_genes = synonymous_mutdensities_all_samples[['GENE', 'synonymous']]
    # elif mode == 'mutated_reads':
    #     synonymous_mutdensities_genes = synonymous_mutdensities_all_samples[['GENE', 'synonymous']]

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

