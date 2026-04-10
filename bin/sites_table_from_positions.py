#!/usr/bin/env python


import click
import pandas as pd

from bgreference import hg38, mm39

assembly_name2function = {"hg38": hg38,
                            "mm39": mm39
                            }

# -- Auxiliary functions -- #
def get_sequence_in_row(chrom, start, length, genome = hg38):
    return genome(chrom, start, size = length)

def get_non_ref(l, letters = {"A", "C", "G", "T"}):
    return letters - set(l)


# -- Main function -- #
def generate_all_sites_4VEP(input_positions, genome, output_prefix, chunk_size=0):

    # read CHROM,POS positions file; check dtypes
    positions_df = pd.read_csv(input_positions, sep = "\t", header = 0,
                                names = ["CHROM", "POS"], dtype = {"CHROM" : str, "POS" : int}
                                )

    genome_func = assembly_name2function[genome]

    # assign REF and all possible ALTs to each positions; add MUTATION and STRAND columnS to meet VEP standards
    positions_df["SEQ"] = positions_df.apply(
        lambda x: get_sequence_in_row(x["CHROM"], x["POS"], 1, genome_func), axis = 1)

    positions_df["ALT"] = positions_df["SEQ"].apply(get_non_ref)
    positions_df = positions_df.explode("ALT").reset_index(drop = True)
    positions_df["MUTATION"] = positions_df["SEQ"].astype(str) + "/" + positions_df["ALT"].astype(str)
    positions_df["STRAND"] = "+"

    # add chr prefix to CHROM
    positions_df["CHROM"] = "chr" + positions_df["CHROM"]

    output_df = positions_df[['CHROM', 'POS', 'POS', 'MUTATION', 'STRAND']]

    # save as chunks or single file
    if chunk_size > 0:
        for i, start in enumerate(range(0, len(output_df), chunk_size)):
            chunk = output_df.iloc[start:start + chunk_size]
            chunk.to_csv(f"{output_prefix}.chunk{i}.tsv",
                         header=False, index=False, sep="\t")
    else:
        output_df.to_csv(f"{output_prefix}.chunk1.tsv",
                         header=False, index=False, sep="\t")



@click.command()
@click.option('--input-positions', required=True, type=click.Path(exists=True), help='Input positions file (TSV)')
@click.option('--genome-assembly', required=True, type=click.Choice(['hg38', 'mm39']), help='Genome assembly (hg38, mm39)')
@click.option('--output-prefix', required=True, type=str, help='Output file prefix (produces <prefix>.chunk<N>.tsv)')
@click.option('--chunk-size', default=0, type=int, help='Number of rows per chunk (0 = no chunking)')
def main(input_positions, genome_assembly, output_prefix, chunk_size):
    generate_all_sites_4VEP(input_positions, genome_assembly, output_prefix, chunk_size)

if __name__ == '__main__':
    main()
