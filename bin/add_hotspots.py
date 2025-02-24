#!/usr/local/bin/python


import sys
import json
import pandas as pd
import numpy as np


panel_file = sys.argv[1]
bedfile = sys.argv[2]  # BED file (optional, can be None)
expand = int(sys.argv[3])  # Expansion factor for defining regions
autoexons = bool(int(sys.argv[4]))  # Whether to use the BED file for exon definitions

# Read input data
panel_data = pd.read_table(panel_file)

if not autoexons:
    # Use BED file for exon definition
    exons_bed = pd.read_table(bedfile, header=None, sep="\t")
    # Check if a name column is provided
    if exons_bed.shape[1] > 3:
        exons_bed = exons_bed.iloc[:,:4]
        exons_bed.columns = ["CHROM", "START", "END", "NAME"]
    else:
        exons_bed = exons_bed.iloc[:,:3]
        exons_bed.columns = ["CHROM", "START", "END"]
        exons_bed["NAME"] = None
    print(exons_bed)

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
        exons_bed.extend([{"CHROM": chrom[0], "START": s, "END": e, "NAME": None, "GENE": chrom[1]} for s, e in zip(starts, ends)])
    exons_bed = pd.DataFrame(exons_bed)


## Process panel data to create new genes per exon
new_data = pd.DataFrame()
hotspots_names = []

current_chr = ""
exon_counters = {}  # Dictionary to track exon numbers for each gene
for ind, row in exons_bed.iterrows():
    try:
        if row["CHROM"] != current_chr:
            current_chr = row["CHROM"]
            chr_data = panel_data[panel_data["CHROM"] == current_chr]#.reset_index(drop = True)
            print("Updating chromosome to:", current_chr)

        # Get start and end indices
        ind_start = np.where(chr_data["POS"] == row["START"])[0][0]
        ind_end = np.where(chr_data.iloc[ind_start:,:]["POS"] == row["END"])[0][0]
        gene = chr_data.iloc[ind_start]["GENE"]

        ind_end += ind_start


        # Initialize exon counter for the gene if not already present
        if gene not in exon_counters:
            exon_counters[gene] = 1

        # Determine the region name
        if pd.notna(row.get("NAME", None)):
            region_name = row["NAME"]
        else:
            exon_number = exon_counters[gene]
            region_name = f"{gene}--exon_{exon_number}"

        if not autoexons:
            # Expand within limits of the same gene only if BED file is used
            upd_start = max(ind_start - expand * 3, 0)
            while upd_start < len(chr_data) and chr_data.iloc[upd_start]["GENE"] != gene:
                upd_start += 1

            upd_end = min(ind_end + expand * 3, len(chr_data) - 1)
            while upd_end >= 0 and chr_data.iloc[upd_end]["GENE"] != gene:
                upd_end -= 1

        else:
            # No expansion for exon boundaries if BED file is not used
            upd_start = ind_start
            upd_end = ind_end

        # Extract hotspot data and modify gene names
        hotspot_data = chr_data.iloc[upd_start: upd_end + 1, :].copy()
        hotspot_data["GENE"] = region_name

        new_data = pd.concat((new_data, hotspot_data))
        hotspots_names.append(region_name)
        print("Small region added:", region_name)

        # Increment exon counter for the gene (if name is generated)
        if row.get("NAME") is None:
            exon_counters[gene] += 1

    except Exception as e:
        print(f"Error processing row {row}: {e}")

# Combine panel data with new data and save
final_data = pd.concat((panel_data, new_data))

final_data.to_csv("exons_consensus_panel_with_hotspots.tsv",
                    sep="\t", header=True, index=False)

with open("hotspot_names.json", "w") as f:
    json.dump({x: [x] for x in hotspots_names}, f, indent=4)
