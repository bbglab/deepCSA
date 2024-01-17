#!/usr/bin/env python3

import pandas as pd
import os
import sys

# -- Auxiliary functions -- #

panel_impact_dict = {

    "protein_affecting": ["nonsense", "missense", "protein_altering_variant",
                          "transcript_amplification", "coding_sequence_variant",
                         "essential_splice"
                         ],

    "non_protein_affecting": ["synonymous", "intron_variant",
                              "non_coding_exon_region", "splice_region",
                              "non_genic_variant", "non_coding_transcript_variant"
                              ],

    "exons_splice_sites": ["missense", "essential_splice", "coding_sequence_variant",
                          "synonymous",
                          "nonsense", "protein_altering_variant", # unclear
                           "transcript_amplification", "non_coding_exon_region", # unclear
                           "non_coding_transcript_variant"   # unclear
                          ],

    "introns_intergenic": ["splice_regions", "intron_variant",
                          "non_genic_variant"   # unclear
                          ]

}

# -- Main function -- #

def create_panel_versions(compact_annot_panel_path, output_path):

    # Load VEP annotated panel, already compacted to have one variant per site
    ## requires column named IMPACT with consequence type
    compact_annot_panel_df = pd.read_csv(compact_annot_panel_path, sep = "\t")

    # Create panel versions
    for version in panel_impact_dict:

        panel_version = compact_annot_panel_df.loc[compact_annot_panel_df["IMPACT"].isin(panel_impact_dict[version])]
        panel_version.to_csv(f"{output_path}.compact.{version}.tsv",
                             sep = "\t", index = False)

    # Store complete panel (better change this way of using this version in nextflow)
    version = "all"
    compact_annot_panel_df.to_csv(f"{output_path}.compact.{version}.tsv",
                             sep = "\t", index = False)

if __name__ == '__main__':
    compact_annot_panel_path = sys.argv[1]
    output_path = sys.argv[2]

    create_panel_versions(compact_annot_panel_path, output_path)
