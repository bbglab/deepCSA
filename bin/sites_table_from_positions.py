#!/usr/bin/env python


import click
import pandas as pd

from bgreference import hg38, hg19, mm10, mm39

assembly_name2function = {"hg38": hg38,
                            "hg19": hg19,
                            "mm10": mm10,
                            "mm39": mm39
                            }

# -- Auxiliary functions -- #
def get_sequence_in_row(chrom, start, length, genome = hg38):
    return genome(chrom, start, size = length)

def get_non_ref(l, letters = {"A", "C", "G", "T"}):
    return letters - set(l)


# -- Main function -- #
def generate_all_sites_4VEP(input_positions, genome, output_file_with_sites):

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

    # save
    positions_df[['CHROM', 'POS', 'POS', 'MUTATION', 'STRAND']].to_csv(output_file_with_sites,
                                                                        header = False,
                                                                        index = False,
                                                                        sep = "\t")



@click.command()
@click.option('--input-positions', required=True, type=click.Path(exists=True), help='Input positions file (TSV)')
#@click.option('--genome-assembly', required=True, type=click.Choice(['hg38', 'hg19', 'mm10', 'mm39']), help='Genome assembly (hg38, hg19, mm10, mm39)')
@click.option('--genome-assembly', required=True, type=click.Choice(['hg38', 'mm39']), help='Genome assembly (hg38, mm39)')
@click.option('--output-file-with-sites', required=True, type=click.Path(), help='Output file for sites (TSV)')
def main(input_positions, genome_assembly, output_file_with_sites):
    generate_all_sites_4VEP(input_positions, genome_assembly, output_file_with_sites)

if __name__ == '__main__':
    main()
