#!/usr/local/bin/python



# TODO
# add plotting modules to bgreference container
import sys
import json
import pandas as pd
import numpy as np
from utils import to_int_if_possible
from read_utils import custom_na_values

# import seaborn as sns
# import matplotlib.pyplot as plt


panel_file = sys.argv[1]
bedfile = sys.argv[2]
expand = int(sys.argv[3]) # 30


panel_data = pd.read_table(panel_file)
hotspots_bed = pd.read_table(bedfile, header = None, sep = '\t')
hotspots_bed.columns = ["CHROM", "START", "END", "NAME"]


new_data = pd.DataFrame()
hotspots_names = []

current_chr = ''
for ind, row in hotspots_bed.iterrows():
    try :
        if row["CHROM"] != current_chr:
            current_chr = row["CHROM"]
            chr_data = panel_data[panel_data["CHROM"] == current_chr]
            print("updating_chr to", current_chr)

        ind_start = np.where(chr_data["POS"] == row["START"])[0][0]
        ind_end = np.where(chr_data.iloc[ind_start:,:]["POS"] == row["END"])[0][0]
        gene = chr_data.iloc[ind_start]["GENE"]

        ind_end += ind_start

        upd_start = ind_start - expand*3

        # shrink the extension to fit within the same gene
        while chr_data.iloc[upd_start]["GENE"] != gene:
            upd_start += 3

        upd_end = ind_end + expand*3
        while chr_data.iloc[upd_end]["GENE"] != gene:
            upd_end -= 3

        print(ind_start, ind_end)
        print(upd_start, upd_end)

        hotspot_data = chr_data.iloc[upd_start: upd_end, :].copy()
        hotspot_data["GENE"] = row["NAME"]

        new_data = pd.concat((new_data, hotspot_data))
        hotspots_names.append(row["NAME"])

    except:
        print(row)


final_data = pd.concat((panel_data, new_data))

final_data.to_csv("exons_consensus_panel.with_hotspots.tsv",
                                        sep = "\t",
                                        header = True,
                                        index = False)

with open("hotspot_names.json", 'w') as f:
    json.dump({ x : [x] for x in hotspots_names }, f, indent=4)
