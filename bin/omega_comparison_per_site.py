#!/usr/local/bin/python

import pandas as pd
import gzip
import click

def load_data(panel_file, mutations_file, mutabilities_file):
    # Load captured panel
    panel = pd.read_csv(panel_file, sep="\t", compression="gzip")
    # Load mutations
    mutations = pd.read_csv(mutations_file, sep="\t")
    if "EFFECTIVE_MUTS" not in mutations.columns:
        mutations["EFFECTIVE_MUTS"] = 1
    # Load mutabilities
    mutabilities = pd.read_csv(mutabilities_file, sep="\t", compression="gzip")
    return panel, mutations, mutabilities

def get_sample_column(mutabilities):
    # Identify the sample column in mutabilities (assumes it's the last column excluding CHROM, POS, etc.)
    non_data_columns = ['CHROM', 'POS', 'CONTEXT_MUT', 'GENE', 'IMPACT']
    sample_columns = [col for col in mutabilities.columns if col not in non_data_columns]
    if len(sample_columns) != 1:
        raise ValueError(f"Expected 1 sample column in mutabilities but found {len(sample_columns)}: {sample_columns}")
    return sample_columns[0]

def compute_by_size(panel, mutations, mutabilities, size, sample_column):
    # Merge mutations and panel to associate protein positions
    mutations = mutations.merge(panel[['CHROM', 'POS', 'REF', 'ALT', 'Protein_position', 'Amino_acids', 'GENE', 'Feature']],
                                on=['CHROM', 'POS', 'REF', 'ALT'], how='inner')

    # Mutabilities: Merge with panel to associate protein positions and transcript info
    mutabilities = mutabilities.merge(panel[['CHROM', 'POS', 'REF', 'ALT', 'CONTEXT_MUT', 'Protein_position', 'Amino_acids', 'GENE', 'Feature']],
                                        on=['CHROM', 'POS', 'CONTEXT_MUT', 'GENE'], how='inner')
    mutabilities = mutabilities.drop('CONTEXT_MUT', axis = 'columns')

    # Define grouping based on size
    if size == "site":
        group_cols = ['CHROM', 'POS', 'REF', 'ALT', 'GENE',]
    elif size == "aminoacid":
        group_cols = ['GENE', 'Feature', 'Protein_position']
    elif size == "aminoacid_change":
        group_cols = ['GENE', 'Feature', 'Protein_position', 'Amino_acids']
    else:
        raise ValueError(f"Invalid size: {size}. Choose 'site', 'aminoacid', or 'aminoacid_change'.")

    # Observed counts
    observed = mutations.groupby(group_cols)['EFFECTIVE_MUTS'].sum().reset_index()
    observed.columns = group_cols + ['OBSERVED_MUTS']

    # Accumulated mutability
    mutability = mutabilities.groupby(group_cols)[sample_column].sum().reset_index()
    mutability.columns = group_cols + ['EXPECTED_MUTS']

    # Merge results
    result = pd.merge(observed, mutability, on=group_cols, how='outer').fillna(0)
    result["OBS-EXP"] = (result["OBSERVED_MUTS"] - result["EXPECTED_MUTS"]).fillna(0)
    result["(OBS-EXP)/EXP"] = ((result["OBSERVED_MUTS"] - result["EXPECTED_MUTS"])/result["EXPECTED_MUTS"]).fillna(0)

    return result




@click.command()
@click.option('--mutations-file', type=click.Path(exists=True), required=True, help="Path to mutations file.")
@click.option('--mutabilities-file', type=click.Path(exists=True), required=True, help="Path to mutabilities file (gzip compressed).")
@click.option('--panel-file', type=click.Path(exists=True), required=True, help="Path to captured panel file (gzip compressed).")
@click.option('--size', type=click.Choice(['site', 'aminoacid', 'aminoacid_change'], case_sensitive=False), default='aminoacid_change', help="Level of aggregation.")
@click.option('--output-file', type=click.Path(writable=True), required=True, help="Path to output file.")
def main(panel_file, mutations_file, mutabilities_file, size, output_file):
    """Compute comparison between observed mutations and accumulated mutabilities."""
    panel, mutations, mutabilities = load_data(panel_file, mutations_file, mutabilities_file)

    # Get sample column dynamically
    sample_column = get_sample_column(mutabilities)
    click.echo(f"Detected sample column: {sample_column}")

    result = compute_by_size(panel, mutations, mutabilities, size, sample_column)
    result.to_csv(output_file, sep="\t", index=False)
    click.echo(f"Results written to {output_file}")

if __name__ == "__main__":
    main()
