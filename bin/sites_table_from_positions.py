#!/usr/bin/env python3

import sys
import pandas as pd

from bgreference import hg38


# -- Auxiliary functions -- #
def get_sequence_in_row(chrom, start, length = 1):
    return hg38(chrom, start, size = length)

def get_non_ref(l, letters = {"A", "C", "G", "T"}):
    return letters - set(l)


# -- Main function -- #
def generate_all_sites_4VEP(input_positions, output_file_with_sites):

    # read CHROM,POS positions file; check dtypes
    positions_df = pd.read_csv(input_positions, sep = "\t", header = 0,
                                names = ["CHROM", "POS"], dtype = {"CHROM" : str, "POS" : int}
                                )

    # assign REF and all possible ALTs to each positions; add MUTATION and STRAND columnS to meet VEP standards
    positions_df["SEQ"] = positions_df.apply(
        lambda x: get_sequence_in_row(x["CHROM"], x["POS"], length = 1), axis = 1)
    positions_df["ALT"] = positions_df["SEQ"].apply(get_non_ref)
    positions_df = positions_df.explode("ALT").reset_index(drop = True)
    positions_df["MUTATION"] = positions_df["SEQ"].astype(str) + "/" + positions_df["ALT"].astype(str)
    positions_df["STRAND"] = "+"

    # save
    positions_df[['CHROM', 'POS', 'POS', 'MUTATION', 'STRAND']].to_csv(output_file_with_sites,
                                                                        header = False,
                                                                        index = False,
                                                                        sep = "\t")


if __name__ == '__main__':
    # Input
    # input_positions = "/workspace/datasets/prominent/metadata/regions/data/oncodrivefml/kidneypanel4oncodrivefml.bed5.bed"
    input_positions = sys.argv[1]

    # Output
    # output_file_with_sites = ./test/preprocessing/KidneyPanel.sites4VEP.tsv"
    output_file_with_sites = sys.argv[2]

    generate_all_sites_4VEP(input_positions, output_file_with_sites)
