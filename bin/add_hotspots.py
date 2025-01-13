#!/usr/local/bin/python



# # TODO
# # add plotting modules to bgreference container
# import sys
# import json
# import pandas as pd
# import numpy as np
# from utils import to_int_if_possible
# from read_utils import custom_na_values

# # import seaborn as sns
# # import matplotlib.pyplot as plt


# panel_file = sys.argv[1]
# bedfile = sys.argv[2]
# expand = int(sys.argv[3])


# panel_data = pd.read_table(panel_file)
# hotspots_bed = pd.read_table(bedfile, header = None, sep = '\t')
# hotspots_bed.columns = ["CHROM", "START", "END", "NAME"]


# new_data = pd.DataFrame()
# hotspots_names = []

# current_chr = ''
# for ind, row in hotspots_bed.iterrows():
#     try :
#         if row["CHROM"] != current_chr:
#             current_chr = row["CHROM"]
#             chr_data = panel_data[panel_data["CHROM"] == current_chr]
#             print("updating_chr to", current_chr)

#         ind_start = np.where(chr_data["POS"] == row["START"])[0][0]
#         ind_end = np.where(chr_data.iloc[ind_start:,:]["POS"] == row["END"])[0][0]
#         gene = chr_data.iloc[ind_start]["GENE"]

#         ind_end += ind_start

#         # shrink the extension to fit within the same gene
#         upd_start = ind_start - expand*3
#         while chr_data.iloc[upd_start]["GENE"] != gene:
#             upd_start += 3

#         upd_end = ind_end + expand*3
#         while chr_data.iloc[upd_end]["GENE"] != gene:
#             upd_end -= 3

#         print(ind_start, ind_end)
#         print(upd_start, upd_end)

#         hotspot_data = chr_data.iloc[upd_start: upd_end, :].copy()
#         hotspot_data["GENE"] = row["NAME"]

#         new_data = pd.concat((new_data, hotspot_data))
#         hotspots_names.append(row["NAME"])

#     except:
#         print(row)


# final_data = pd.concat((panel_data, new_data))

# final_data.to_csv("exons_consensus_panel.with_hotspots.tsv",
#                                         sep = "\t",
#                                         header = True,
#                                         index = False)

# with open("hotspot_names.json", 'w') as f:
#     json.dump({ x : [x] for x in hotspots_names }, f, indent=4)

import sys
import json
import pandas as pd
import numpy as np
from utils import to_int_if_possible
from read_utils import custom_na_values

panel_file = sys.argv[1]
bedfile = sys.argv[2]  # BED file (optional, can be None)
expand = int(sys.argv[3])  # Expansion factor for defining regions
use_bed = bool(int(sys.argv[4]))  # Whether to use the BED file for exon definitions

# Read input data
panel_data = pd.read_table(panel_file)

if use_bed:
    # Use BED file for exon definition
    exons_bed = pd.read_table(bedfile, header=None, sep="\t")
    exons_bed.columns = ["CHROM", "START", "END", "GENE"]
else:
    # Use breakpoints to define exons
    exons_bed = []
    for chrom, gene_group in panel_data.groupby(["CHROM", "GENE"]):
        positions = gene_group["POS"].sort_values().to_numpy()
        starts = [positions[0]]
        ends = []

        # Identify breaks
        for i in range(1, len(positions)):
            if positions[i] - positions[i-1] > 1:
                ends.append(positions[i-1])
                starts.append(positions[i])
        ends.append(positions[-1])

        # Append to exons
        exons_bed.extend([{"CHROM": chrom[0], "START": s, "END": e, "GENE": chrom[1]} for s, e in zip(starts, ends)])
    exons_bed = pd.DataFrame(exons_bed)

print(exons_bed)

## Process panel data to create new genes per exon
new_data = pd.DataFrame()
hotspots_names = []

current_chr = ""
exon_counters = {}  # Dictionary to track exon numbers for each gene
for ind, row in exons_bed.iterrows():
    try:
        if row["CHROM"] != current_chr:
            current_chr = row["CHROM"]
            chr_data = panel_data[panel_data["CHROM"] == current_chr]
            print("Updating chromosome to:", current_chr)

        # Get start and end indices
        ind_start = np.searchsorted(chr_data["POS"], row["START"])
        ind_end = np.searchsorted(chr_data["POS"], row["END"], side="right") - 1
        gene = row["GENE"]

        # Initialize exon counter for the gene if not already present
        if gene not in exon_counters:
            exon_counters[gene] = 1
        exon_number = exon_counters[gene]

        if use_bed:
            # Expand within limits of the same gene only if BED file is used
            upd_start = max(ind_start - expand * 3, 0)
            while chr_data.iloc[upd_start]["GENE"] != gene and upd_start < len(chr_data):
                upd_start += 1

            upd_end = min(ind_end + expand * 3, len(chr_data) - 1)
            while chr_data.iloc[upd_end]["GENE"] != gene and upd_end >= 0:
                upd_end -= 1
        else:
            # No expansion for exon boundaries if BED file is not used
            upd_start = ind_start
            upd_end = ind_end

        # Extract hotspot data and modify gene names
        hotspot_data = chr_data.iloc[upd_start: upd_end + 1, :].copy()
        hotspot_data["GENE"] = f"{gene}--exon_{exon_number}"

        new_data = pd.concat((new_data, hotspot_data))
        hotspots_names.append(f"{gene}--exon_{exon_number}")

        # Increment exon counter for the gene
        exon_counters[gene] += 1

    except Exception as e:
        print(f"Error processing row {row}: {e}")

# Combine panel data with new data and save
final_data = pd.concat((panel_data, new_data))

final_data.to_csv("exons_consensus_panel_with_hotspots.tsv",
                  sep="\t", header=True, index=False)

with open("hotspot_names.json", "w") as f:
    json.dump({x: [x] for x in hotspots_names}, f, indent=4)
