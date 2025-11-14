#!/usr/bin/env python

import os
import click
import pandas as pd
from read_utils import custom_na_values


def merge_cohort_files(output_file):
    """
    Merge all .tsv.gz files in the current directory into a single file.
    """
    directory = '.'
    file_names = [file for file in os.listdir(directory) if file.endswith('.tsv.gz')]

    dfs = []
    for file_name in file_names:
        file_path = os.path.join(directory, file_name)
        df = pd.read_csv(file_path, compression='gzip', header=0, sep='\t', na_values=custom_na_values)
        dfs.append(df)
    concatenated_df = pd.concat(dfs, ignore_index=True)

    concatenated_df.to_csv(output_file, sep="\t", header=True, index=False)


@click.command()
@click.option('--output_file', type=click.Path(), required=True, help='Output file name for the merged cohort.')
def main(output_file):
    """
    CLI entry point for merging cohort files.
    """
    click.echo(f"Merging files from the current directory.")
    merge_cohort_files(output_file)
    click.echo(f"Merged cohort saved to: {output_file}")


if __name__ == '__main__':
    main()


