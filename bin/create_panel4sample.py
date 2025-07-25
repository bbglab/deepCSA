#!/usr/bin/env python

import pandas as pd
import sys

# -- Main function -- #

def create_panel4sample(compact_annot_panel_path, depths_path, panel_name, min_depth):

    # Load captured panel and depths
    compact_annot_panel_df = pd.read_csv(compact_annot_panel_path, sep = "\t")
    # version = compact_annot_panel_path.split("/")[-1].split(".")[2]

    depths_df = pd.read_csv(depths_path, sep = "\t")
    depths_df.iloc[:, 1:] = depths_df.iloc[:, 1:].astype(int)

    # Generate panel subset per sample based on min_depth threshold
    min_depth = int(min_depth)
    for sample in depths_df.columns[2:]:

        sample_depth = depths_df[["CHROM", "POS", sample]]
        sample_depth = sample_depth.loc[sample_depth[sample] >= min_depth]
        sample_panel = compact_annot_panel_df.merge(sample_depth, on = ["CHROM", "POS"], how = "inner").drop(sample, axis = 1)

        sample_panel.to_csv(f"{sample.split('.')[0]}.{panel_name}.tsv", index = False, sep = "\t")


if __name__ == '__main__':
    compact_annot_panel_path = sys.argv[1]
    depths_path = sys.argv[2]
    panel_v = sys.argv[3]
    min_depth = sys.argv[4]

    create_panel4sample(compact_annot_panel_path, depths_path, panel_v, min_depth)
