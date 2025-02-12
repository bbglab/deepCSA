#!/usr/local/bin/python


import click
import pandas as pd
import numpy as np

@click.group()
def cli():
    pass

@cli.command()
@click.argument('depths')
@click.option('--file', required=True, type=click.Path(exists=True), help='Input file with depths.')
@click.option('--proportion', required=True, type=float, help='Proportion to downsample depth values.')
def downsample(depths, file, proportion):
    """
    Downsample depth values in the provided file based on the given proportion.

    FILE format: tab-delimited with columns: chromosome, position, depth.
    """

    if not (0 < proportion <= 1):
        click.echo("Error: Proportion must be between 0 and 1.")
        return

    # Load the input file
    df = pd.read_csv(file, sep='\t', header=0)
    sample_name = df.columns[-1]

    # Downsample the depth values by multiplying and rounding
    df[sample_name] = ((df[sample_name] * proportion) // 1).astype(int)

    # Output the downsampled data
    output_file = f"{sample_name}.depths.downsampled.tsv.gz"
    df.to_csv(output_file, sep='\t', index=False)

    click.echo(f"Downsampled file saved as: {output_file}")



@cli.command()
@click.argument('mutations')
@click.option('--file', required=True, type=click.Path(exists=True), help='Input file with mutations info.')
@click.option('--proportion', required=True, type=float, help='Proportion to downsample depth values.')
@click.option('--samplename', required=True, type=str, help='Sample name for output file.')
def downsample_mutations(mutations, file, proportion, samplename):
    """
    Downsample mutation depths and filter mutations based on VAF.

    FILE format: tab-delimited with columns:
    """

    if not (0 < proportion <= 1):
        click.echo("Error: Proportion must be between 0 and 1.")
        return

    # Load the input file
    df = pd.read_csv(file, sep='\t', header=0)

    # Downsample the depth values by multiplying and rounding
    df["DEPTH"] = ((df["DEPTH"] * proportion) // 1).astype(int)

    # Use Poisson distribution to decide if mutation remains
    df['upd_ALT_DEPTH'] = df.apply(lambda row: np.random.poisson(row['VAF'], row['DEPTH']).sum(), axis=1)
    df['ALT_DEPTH'] = df['upd_ALT_DEPTH']
    df['retain_mutation'] = df['upd_ALT_DEPTH'] > 0

    # Filter to retain only mutations
    filtered_df = df[df['retain_mutation']].drop(columns=['retain_mutation', 'upd_ALT_DEPTH'])
    filtered_df["VAF"] = df['ALT_DEPTH'] / df["DEPTH"]

    # Output the downsampled data
    output_file = f"{samplename}.downsampled.mutations.tsv"
    filtered_df.to_csv(output_file, sep='\t', index=False)

    click.echo(f"Downsampled mutations file saved as: {output_file}")

if __name__ == '__main__':
    cli()
