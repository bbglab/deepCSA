#!/usr/bin/env python


import click
import pandas as pd
from read_utils import custom_na_values

# -- Main function -- #
def prepare_depths(depths_path, annot_panel_path, filee):

    # File loading
    depths_df = pd.read_csv(depths_path, sep = "\t")
    depths_df = depths_df.drop("CONTEXT", axis = 1)
    annot_panel_df = pd.read_csv(annot_panel_path, sep = "\t", na_values = custom_na_values)
    depth_column_name = depths_df.columns[2]

    # Subset depths with panel
    ## mode 1: each position counts one
    depths_subset_df = depths_df.merge(annot_panel_df[["CHROM", "POS", "GENE"]].drop_duplicates(),
                                        on = ["CHROM", "POS"], how = "inner")
    del depths_df
    del annot_panel_df

    depths_gene_df = depths_subset_df.groupby("GENE").agg({depth_column_name: 'mean'})
    depths_gene_df.columns = ["DEPTH"]
    depths_gene_df.to_csv(f"{filee}", sep = '\t', header = False, index = True)
    # ARID1A  10859.146
    # BAP1    15240.987
    # MTOR    12738.001
    # PBRM1   13411.108




@click.command()
@click.option('--depths-path', required=True, type=click.Path(exists=True), help='Input depths file (TSV)')
@click.option('--annot-panel-path', required=True, type=click.Path(exists=True), help='Input annotation panel file (TSV)')
@click.option('--output', required=True, type=click.Path(), help='Output file for gene depths (TSV)')
def main(depths_path, annot_panel_path, output):
    """
    Prepare gene depths from depths and annotation panel files.
    """
    prepare_depths(depths_path, annot_panel_path, output)

if __name__ == '__main__':
    main()


