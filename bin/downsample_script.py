#!/usr/bin/env python


import click
import pandas as pd
import numpy as np


# Perform uniform sampling without replacement from a pool of 0s and 1s
def uniform_sampling(row, suffix=''):
    ref_depth = row[f'orig_DEPTH{suffix}'] - row[f'ALT_DEPTH{suffix}']
    pool = [1] * row[f'ALT_DEPTH{suffix}'] + [0] * ref_depth
    sampled = np.random.choice(pool, size=row[f'DEPTH{suffix}'], replace=False)
    return sampled.sum()

@click.group()
def cli():
    pass

@cli.command(name='depths')
@click.option('--file', required=True, type=click.Path(exists=True), help='Input file with depths.')
@click.option('--proportion', required=True, type=float, help='Proportion to downsample depth values.')
def downsample_depths(file, proportion):
    """
    Downsample depth values in the provided file based on the given proportion.

    FILE format: tab-delimited with columns: chromosome, position, depth.
    """

    if not (0 < proportion <= 1):
        click.echo("Error: Proportion must be between 0 and 1.")
        return

    # Load the input file
    df = pd.read_csv(file, sep='\t', header=0)
    sample_name = file.split('.')[0]  # Extract sample name from file name
    column_names = df.columns[3:]

    # Downsample the depth values by multiplying and rounding
    df[column_names] = ((df[column_names] * proportion) // 1).astype(int)

    # Output the downsampled data
    output_file = f"{sample_name}.depths.downsampled.tsv.gz"
    df.to_csv(output_file, sep='\t', index=False)

    click.echo(f"Downsampled file saved as: {output_file}")



@cli.command(name='mutations')
@click.option('--file', required=True, type=click.Path(exists=True), help='Input file with mutations info.')
@click.option('--proportion', required=True, type=float, help='Proportion to downsample depth values.')
@click.option('--samplename', required=True, type=str, help='Sample name for output file.')
def downsample_mutations(file, proportion, samplename):
    """
    Downsample mutation depths and filter mutations based on VAF.

    FILE format: tab-delimited with columns:
    """

    if not (0 < proportion <= 1):
        click.echo("Error: Proportion must be between 0 and 1.")
        return

    # Load the input file
    df = pd.read_csv(file, sep='\t', header=0)
    df = df.rename(columns={"DEPTH": "orig_DEPTH"})

    # Downsample the depth values by multiplying and rounding
    df["DEPTH"] = ((df["orig_DEPTH"] * proportion) // 1).astype(int)

    # Use Poisson distribution to decide if mutation remains
    #df['upd_ALT_DEPTH'] = df.apply(lambda row: np.random.poisson(row['VAF'], row['DEPTH']).sum(), axis=1)
    df['ALT_DEPTH'] = df.apply(uniform_sampling, suffix='', axis=1)
    df['retain_mutation'] = df['ALT_DEPTH'] > 0

    # Filter to retain only mutations
    filtered_df = df[df['retain_mutation']].drop(columns=['retain_mutation', 'orig_DEPTH'])
    filtered_df["VAF"] = df['ALT_DEPTH'] / df["DEPTH"]

    if filtered_df.empty:
        click.echo("No mutations retained after downsampling.")
        return

    if "DEPTH_AM" in filtered_df.columns and "ALT_DEPTH_AM" in filtered_df.columns:
        # Downsample the all-molecules counts as well
        df = filtered_df.rename(columns={"DEPTH_AM": "orig_DEPTH_AM"})

        # Downsample the depth values by multiplying and rounding
        df["DEPTH_AM"] = ((df["orig_DEPTH_AM"] * proportion) // 1).astype(int)

        # Use Poisson distribution to decide if mutation remains
        df['ALT_DEPTH_AM'] = df.apply(uniform_sampling, suffix='_AM', axis=1)
        
        # the AM counts cannot be lower than the original counts, so we take the max of the two
        df['ALT_DEPTH_AM'] = df[['ALT_DEPTH', 'ALT_DEPTH_AM']].max(axis=1)

        # Filter to retain only mutations
        filtered_df = df.drop(columns=['orig_DEPTH_AM'])
        filtered_df["VAF_AM"] = df['ALT_DEPTH_AM'] / df["DEPTH_AM"]

    # Output the downsampled data
    output_file = f"{samplename}.downsampled.mutations.tsv"
    filtered_df.to_csv(output_file, sep='\t', index=False)

    click.echo(f"Downsampled mutations file saved as: {output_file}")

if __name__ == '__main__':
    cli()
