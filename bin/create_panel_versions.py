#!/usr/bin/env python3

"""
create_panel_versions_polars.py

Generates multiple VEP annotation panel subsets based on the 'IMPACT' column
using the high-performance Polars library.

Usage:
    python create_panel_versions_polars.py <input_tsv> <output_prefix>
"""

import polars as pl
import click
import os
import sys

PANEL_IMPACT_DICT = {

    "protein_affecting": ["nonsense", "missense",
                            "essential_splice",
                            "protein_altering_variant",     # probably not appear
                            "transcript_amplification",     # probably not appear
                            # "coding_sequence_variant",       # probably not appear AMBIGUOUS TODO
                            # "splice_region_variant",
                            # "splice_donor_region_variant",
                            # "splice_polypyrimidine_tract_variant",
                            ],

    "non_protein_affecting": ["synonymous", "intron_variant",
                                "non_coding_exon_region",
                                "non_genic_variant", "non_coding_transcript_variant",
                                # "splice_region_variant",
                                # "splice_donor_region_variant",
                                # "splice_polypyrimidine_tract_variant",
                                ],

    "exons_splice_sites": ["nonsense", "missense",
                            "essential_splice",
                            "synonymous",
                            "coding_sequence_variant",
                            "splice_region_variant",
                            "splice_donor_region_variant",
                            "splice_polypyrimidine_tract_variant",
                        #    "protein_altering_variant", # unclear   # probably not appear
                        #    "transcript_amplification",             # probably not appear
                        #    "non_coding_exon_region", # unclear
                        #    "non_coding_transcript_variant"   # unclear
                            ],

    "introns_intergenic": ["intron_variant",
                            "non_genic_variant",
                            # "splice_donor_region_variant",
                            # "splice_polypyrimidine_tract_variant",
                            ],

    "synonymous": ["synonymous"]

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


@click.command()
@click.argument("input_path", type=click.Path(exists=True))
@click.argument("output_prefix", type=str)
def create_panel_versions(input_path: str, output_prefix: str) -> None:
    """
    Generates panel subsets from a VEP-annotated file using Polars.

    \b
    INPUT_PATH: Path to the annotated TSV file.
    OUTPUT_PREFIX: Prefix for the output files (e.g., 'output/panel').
    """
    try:
        df = pl.read_csv(input_path, separator="\t")
    except Exception as e:
        click.echo(f"Error reading input file: {e}", err=True)
        sys.exit(1)

    if "IMPACT" not in df.columns:
        click.echo("ERROR: 'IMPACT' column not found in input file.", err=True)
        sys.exit(1)

    for version_name, impact_values in PANEL_IMPACT_DICT.items():
        filtered = df.filter(pl.col("IMPACT").is_in(impact_values))
        filtered.write_csv(f"{output_prefix}.{version_name}.tsv", separator="\t")

    # Write the full file as a version
    df.write_csv(f"{output_prefix}.all.tsv", separator="\t")

    click.echo("Panel versions generated successfully.")

if __name__ == "__main__":
    create_panel_versions()