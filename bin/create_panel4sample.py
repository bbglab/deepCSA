#!/usr/bin/env python


import click
import pandas as pd

# -- Main function -- #

def create_panel4sample(compact_annot_panel_path, depths_path, panel_name, min_depth):

    # Load captured panel and depths
    compact_annot_panel_df = pd.read_csv(compact_annot_panel_path, sep = "\t")

    depths_df = pd.read_csv(depths_path, sep = "\t")
    depths_df.iloc[:, 1:] = depths_df.iloc[:, 1:].astype(int)

    # Generate panel subset per sample based on min_depth threshold
    min_depth = int(min_depth)
    for sample in depths_df.columns[2:]:

        sample_depth = depths_df[["CHROM", "POS", sample]]
        sample_depth = sample_depth.loc[sample_depth[sample] >= min_depth]
        sample_panel = compact_annot_panel_df.merge(sample_depth, on = ["CHROM", "POS"], how = "inner").drop(sample, axis = 1)

        sample_panel.to_csv(f"{sample.split('.')[0]}.{panel_name}.tsv", index = False, sep = "\t")



@click.command()
@click.option('--compact-annot-panel-path', required=True, type=click.Path(exists=True), help='Input compact annotation panel file (TSV)')
@click.option('--depths-path', required=True, type=click.Path(exists=True), help='Input depths file (TSV)')
@click.option('--panel-name', required=True, type=str, help='Panel name for output files')
@click.option('--min-depth', required=True, type=int, help='Minimum depth threshold')
def main(compact_annot_panel_path, depths_path, panel_name, min_depth):
    create_panel4sample(compact_annot_panel_path, depths_path, panel_name, min_depth)

if __name__ == '__main__':
    main()
