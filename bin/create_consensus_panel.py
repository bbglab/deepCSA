#!/usr/local/bin/python

import pandas as pd
import sys

# -- Main function -- #

def create_consensus_panel(compact_annot_panel_path, depths_path, consensus_min_depth):

    # Load captured panel and depths
    compact_annot_panel_df = pd.read_csv(compact_annot_panel_path, sep = "\t")

    depths_df = pd.read_csv(depths_path, sep = "\t")
    depths_df.iloc[:, 1:] = depths_df.iloc[:, 1:].astype(int)

    # Keep only positions having the minimum consensus depth in all the samples
    consensus_min_depth = int(consensus_min_depth)
    min_depths_df = depths_df.loc[(depths_df.iloc[:, 2:] >= consensus_min_depth).all(axis = 1)]

    # Filter captured panel to only keep minimally covered positions
    consensus_panel = compact_annot_panel_df.merge(min_depths_df[["CHROM", "POS"]], on = ["CHROM", "POS"], how = "inner")
    version = compact_annot_panel_path.split("/")[-1].split(".")[2]
    consensus_panel.to_csv(f"consensus_panel.{version}.tsv", sep = "\t", index = False)

if __name__ == '__main__':
    compact_annot_panel_path = sys.argv[1]
    depths_path = sys.argv[2]
    min_depth = sys.argv[3]

    create_consensus_panel(compact_annot_panel_path, depths_path, min_depth)
