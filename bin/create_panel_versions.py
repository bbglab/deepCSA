#!/usr/bin/env python3

import pandas as pd
import os
import sys

# -- Auxiliary functions -- #

panel_impact_dict = {

    "protein_affecting": ["nonsense", "missense",
                            "essential_splice",
                            "protein_altering_variant",     # probably not appear
                            "transcript_amplification",     # probably not appear
                            "coding_sequence_variant"       # probably not appear AMBIGUOUS TODO
                            ],

    "non_protein_affecting": ["synonymous", "intron_variant",
                                "non_coding_exon_region", "splice_region",
                                "non_genic_variant", "non_coding_transcript_variant"
                                ],

    "exons_splice_sites": ["nonsense", "missense",
                            "essential_splice",
                            "synonymous",
                            "coding_sequence_variant"

                        #    "protein_altering_variant", # unclear   # probably not appear
                        #    "transcript_amplification",             # probably not appear
                        #    "non_coding_exon_region", # unclear
                        #    "non_coding_transcript_variant"   # unclear
                            ],

    "introns_intergenic": ["splice_regions", "intron_variant",
                            "non_genic_variant"   # unclear
                            ]

## CONTENT of the kidney regions
#       1 IMPACT
#    2787 essential_splice
#    7276 nonsense
#  119061 missense

#   36745 synonymous

#  174240 non_coding_exon_region

#  685875 intron_variant

#   74913 non_coding_transcript_variant
#  436551 non_genic_variant

#   14662 splice_region

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
