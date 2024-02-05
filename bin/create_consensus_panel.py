#!/usr/local/bin/python

import pandas as pd
import sys

# -- Main function -- #

def create_consensus_panel(compact_annot_panel_path, depths_path, version, consensus_min_depth):

    # Load captured panel and depths
    compact_annot_panel_df = pd.read_csv(compact_annot_panel_path, sep = "\t")

    depths_df = pd.read_csv(depths_path, sep = "\t")
    depths_df.iloc[:, 1:] = depths_df.iloc[:, 1:].astype(int)

    # Keep only positions having the minimum consensus depth in all the samples
    consensus_min_depth = int(consensus_min_depth)
    min_depths_df = depths_df.loc[(depths_df.iloc[:, 2:] >= consensus_min_depth).all(axis = 1)]

    # Filter captured panel to only keep minimally covered positions
    consensus_panel = compact_annot_panel_df.merge(min_depths_df[["CHROM", "POS"]], on = ["CHROM", "POS"], how = "inner")
    consensus_panel.to_csv(f"consensus.{version}.tsv", sep = "\t", index = False)

if __name__ == '__main__':
    compact_annot_panel_path = sys.argv[1]
    depths_path = sys.argv[2]
    version = sys.argv[3]
    min_depth = sys.argv[4]


    create_consensus_panel(compact_annot_panel_path, depths_path, version, min_depth)
