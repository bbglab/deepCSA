#!/usr/bin/env python


import click
import pandas as pd
from read_utils import custom_na_values


def select_syn_mutrate(mutrate_file, output_file, mode):
    """
    INFO
    """

    mutrate_df = pd.read_csv(mutrate_file, sep = "\t", header = 0, na_values = custom_na_values)

    synonymous_mutrates_all_samples = mutrate_df[(mutrate_df["MUTTYPES"] == "SNV") &
                                                    (mutrate_df["GENE"] != "ALL_GENES")].reset_index(drop = True)

    if mode == 'mutations':
        synonymous_mutrates_genes = synonymous_mutrates_all_samples[['GENE', 'MUTRATE_MB_ADJUSTED']]
    elif mode == 'mutated_reads':
        synonymous_mutrates_genes = synonymous_mutrates_all_samples[['GENE', 'MUTREADSRATE_MB_ADJUSTED']]

    ## FIXME not sure if this is really needed since when called through main()
    # the input would have already been forced to be either of the two options
    # it might still be useful in case this was not called from main()
    else:
        print('unknown mode, please enter either mutations or mutated_reads')
        exit(1) ## FIXME not sure if this is the right code to exit with

    synonymous_mutrates_genes.columns = ["GENE", "MUTRATE"]
    synonymous_mutrates_genes.to_csv(f"{output_file}",
                                        header=True,
                                        index=False,
                                        sep="\t")


@click.command()
@click.option('--mutrates', type=click.Path(exists=True), help='Input mutation rate file')
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--mode', type=click.Choice(['mutations', 'mutated_reads']), default='mutations')

def main(mutrates, output, mode):
    click.echo("Selecting the gene synonymous mutation rates...")
    select_syn_mutrate(mutrates, output, mode)

if __name__ == '__main__':
    main()

