#!/usr/bin/env python



import sys
import pandas as pd
from utils import add_filter, to_int_if_possible
from read_utils import custom_na_values



def negative_filter_panel_regions(mutations_df, positions_df, filtername):
    """
    Negative filter
    """

    positions_df["not_in_panel"] = False

    # filter mutations not inside the panel
    mutations_df = mutations_df.merge(positions_df, on = ["CHROM", "POS"], how = 'left')
    mutations_df["not_in_panel"] = mutations_df["not_in_panel"].fillna(True)
    mutations_df["FILTER"] = mutations_df[["FILTER","not_in_panel"]].apply(
                                                                    lambda x: add_filter(x["FILTER"], x["not_in_panel"], filtername),
                                                                    axis = 1
                                                                    )

    return mutations_df.drop("not_in_panel", axis = 1)




def filter_panel_regions(mutations_df, positions_df, filtername):
    """
    Positive filter
    """

    positions_df["in_panel"] = True

    # filter mutations not inside the panel
    mutations_df = mutations_df.merge(positions_df, on = ["CHROM", "POS"], how = 'left')
    mutations_df["in_panel"] = mutations_df["in_panel"].fillna(False)
    mutations_df["FILTER"] = mutations_df[["FILTER","in_panel"]].apply(
                                                                    lambda x: add_filter(x["FILTER"], x["in_panel"], filtername),
                                                                    axis = 1
                                                                    )

    return mutations_df.drop("in_panel", axis = 1)



## TODO reimplement with click
sample_maf_file = sys.argv[1]
bedfile = sys.argv[2]
filtername = sys.argv[3]
positive = False

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
positions_df = positions_df.reset_index(drop = True)

# adjust the CHROM field to adapt to the way it is being represented in the mutations list
if sample_maf.iloc[0,0].startswith("chr") and not positions_df.iloc[0,0].startswith("chr"):
    positions_df["CHROM"] = "chr" + positions_df["CHROM"]

elif not sample_maf.iloc[0,0].startswith("chr") and positions_df.iloc[0,0].startswith("chr"):
    positions_df["CHROM"] = positions_df["CHROM"].str.replace("chr", "")


if positive:
    filtered_maf = filter_panel_regions(sample_maf, positions_df, filtername)
else:
    filtered_maf = negative_filter_panel_regions(sample_maf, positions_df, filtername)

filtered_maf.to_csv(f"{'.'.join(sample_maf_file.split('.')[:-2])}.filtered.tsv.gz",
                                        sep = "\t",
                                        header = True,
                                        index = False)
