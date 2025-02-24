#!/usr/local/bin/python

import click
import pandas as pd

def blacklist_mutations(mutations_file, blacklist_file, output_file):
    """
    Remove mutations listed in the blacklist file from the mutations file.
    """
    # Load the mutations DataFrame
    mutations_df = pd.read_csv(mutations_file, sep="\t", header=0)

    # Load the blacklist mutation IDs
    with open(blacklist_file, 'r') as f:
        blacklist_ids = set(line.strip() for line in f)

    # Filter out the blacklisted mutations
    filtered_mutations_df = mutations_df[~mutations_df['MUT_ID'].isin(blacklist_ids)]

    # Write the filtered mutations to the output file
    filtered_mutations_df.to_csv(output_file, sep="\t", index=False)

@click.command()
@click.option('--mutations_file', type=click.Path(exists=True), help='Input mutations file')
@click.option('--blacklist_file', type=click.Path(exists=True), help='File with list of mutation IDs to remove')
@click.option('--output_file', type=click.Path(), help='Output file for filtered mutations')
def main(mutations_file, blacklist_file, output_file):
    click.echo(f"Blacklisting mutations...")
    blacklist_mutations(mutations_file, blacklist_file, output_file)

if __name__ == '__main__':
    main()