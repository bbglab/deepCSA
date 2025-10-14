#!/usr/bin/env python

import click
import pandas as pd
import re
from utils import add_filter, to_int_if_possible
from read_utils import custom_na_values


def filter_panel(mutations_df, positions_df, filtername, positive):
    """
    If positive is True, marks mutations inside the panel.
    If positive is False, marks mutations outside the panel.
    """
    col = "in_panel" if positive else "not_in_panel"
    val = True if positive else False
    fill = False if positive else True

    # Column "in_panel" == True if positive True, column "not_in_panel" == False if positive False
    positions_df[col] = val
    mutations_df = mutations_df.merge(positions_df, on=["CHROM", "POS"], how='left')
    # If positive True, fill NaN with False ("in_panel" == False), if positive False, fill NaN with True ("not_in_panel" == True)
    mutations_df[col] = mutations_df[col].fillna(fill)
    mutations_df["FILTER"] = mutations_df[["FILTER", col]].apply(
                                                            lambda x: add_filter(x["FILTER"], x[col], filtername),
                                                            axis=1
                                                            )
    return mutations_df.drop(col, axis=1)

def remove_non_canonical_chromosomes(positions_df):
    """ Keeps only canonical chromosomes"""
    # Regular expression to match standard chromosomes
    # Only works for human and mouse genomes
    canonical_pattern = re.compile(r'^(chr)?([0-9][0-9]?|X|Y|M(T)?)$', re.IGNORECASE)
    positions_df = positions_df[positions_df["CHROM"].apply(lambda x: canonical_pattern.match(x) is not None)].reset_index(drop=True)

    return positions_df

@click.command()
@click.option('--sample-maf-file', required=True, type=click.Path(exists=True), help='Input sample MAF file (TSV)')
@click.option('--bedfile', required=True, type=click.Path(exists=True), help='Input BED file')
@click.option('--filtername', required=True, type=str, help='Filter name to use')
@click.option('--positive', is_flag=True, default=False, help='Use positive filter (default: negative filter)')

def main(sample_maf_file, bedfile, filtername, positive):
    sample_maf = pd.read_csv(sample_maf_file, sep = '\t', header = 0, na_values = custom_na_values)

    current_filters = pd.unique(sample_maf["FILTER"].astype(str).str.split(";").explode())
    if filtername in current_filters:
        print("Not filtering with this BED file since the provided filter name is already present.")
        exit(1)

    # read BED file
    panel_reg = pd.read_csv(bedfile, sep = "\t", header = None)
    # check if there is a header or not
    first_coord = panel_reg.iloc[0,1]
    if to_int_if_possible(first_coord):
        panel_reg = panel_reg.iloc[:,:3]
    # it means there is a header, and we don't want it
    else:
        panel_reg = panel_reg.iloc[1:,:3]
    panel_reg.columns = ["CHROM", "START", "END"]
    panel_reg["CHROM"] = panel_reg["CHROM"].astype(str)
    panel_reg[["START", "END"]] = panel_reg[["START", "END"]].astype(int)
    panel_reg["POS"] = [ list(range(x+1, y+1)) for x, y in panel_reg[["START", "END"]].values ]
    positions_df = panel_reg.explode("POS").reset_index(drop = True)
    positions_df = positions_df[["CHROM", "POS"]].drop_duplicates()
    positions_df = remove_non_canonical_chromosomes(positions_df)

    # adjust the CHROM field to adapt to the way it is being represented in the mutations list
    if sample_maf.iloc[0,0].startswith("chr") and not positions_df.iloc[0,0].startswith("chr"):
        positions_df["CHROM"] = "chr" + positions_df["CHROM"]
    elif not sample_maf.iloc[0,0].startswith("chr") and positions_df.iloc[0,0].startswith("chr"):
        positions_df["CHROM"] = positions_df["CHROM"].str.replace("chr", "")
    if positive:
        filtered_maf = filter_panel(sample_maf, positions_df, filtername, positive = True)
    else:
        filtered_maf = filter_panel(sample_maf, positions_df, filtername, positive = False)

    filtered_maf.to_csv(f"{'.'.join(sample_maf_file.split('.')[:-2])}.filtered.tsv.gz",
                                            sep = "\t",
                                            header = True,
                                            index = False)

if __name__ == '__main__':
    main()