#!/usr/bin/env python3

import sys
import click
import pandas as pd

from utils import contexts_no_change


def compute_trinucleotides(sample_name, depths_file, pseudocount = 0):
    """
    Compute mutational profile from the input data
          ***Remember to add some pseudocounts to the computation***

        Required information:
            Annotated mutations observed
        Output:
            Mutational profile per sample, possibility to add pseudocounts to prevent some probabilities from being 0
    """

    # Load your MAF DataFrame (raw_annotated_maf)
    depths_annotated = pd.read_csv(depths_file, sep = "\t", header = 0)

    # create the matrix in the desired order
    empty_matrix = pd.DataFrame(index = contexts_no_change)

    samples = [ x for x in depths_annotated.columns if x not in ["CHROM", "POS", "CONTEXT"] ]
    trinucleotides_per_sample = depths_annotated.groupby(by = "CONTEXT")[samples].sum().fillna(0)
    trinucleotides_per_sample = pd.concat( (empty_matrix, trinucleotides_per_sample) , axis = 1)
    if '-' in trinucleotides_per_sample.index:
        trinucleotides_per_sample = trinucleotides_per_sample.drop("-", axis = 0)
    trinucleotides_per_sample.index.name = "CONTEXT"
    trinucleotides_per_sample = trinucleotides_per_sample.reset_index(drop = False)
    trinucleotides_per_sample.columns = ["CONTEXT", sample_name]

    trinucleotides_per_sample[sample_name] = trinucleotides_per_sample[sample_name] + pseudocount


    trinucleotides_per_sample[["CONTEXT", sample_name]].to_csv(f"{sample_name}.trinucleotides.tsv.gz",
                                                                sep = "\t",
                                                                header = True,
                                                                index = False)


@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--depths_file', type=click.Path(exists=True), help='Input depths file')
# @click.option('--out_matrix', type=click.Path(), help='Output mutation matrix file')
# @click.option('--json_filters', type=click.Path(exists=True), help='Input mutation filtering criteria file')
# @click.option('--method', type=click.Choice(['unique', 'multiple']), default='unique')
@click.option('--pseud', type=float, default=0.5)

# @click.option('--mutation_matrix', type=click.Path(exists=True), help='Mutation matrix file (for profile mode)')
# @click.option('--trinucleotide_counts', type=click.Path(exists=True), help='Trinucleotide counts file (for profile mode)')
# @click.option('--out_profile', type=click.Path(), help='JSON output file (for profile mode)')
# @click.option('--plot', is_flag=True, help='Generate plot and save as PDF')

# def main(mode, sample_name, mut_file, out_matrix, json_filters, method, pseud, mutation_matrix, trinucleotide_counts, out_profile, plot):
def main(sample_name, depths_file, pseud):
    click.echo(f"Running the trinucleotide computation...")
    click.echo(f"Using the pseudocount: {pseud}")
    compute_trinucleotides(sample_name, depths_file, pseud )
    click.echo("Trinucleotides computation completed.")

if __name__ == '__main__':
    main()

