#!/usr/bin/env python

import click
import json
import pandas as pd
import numpy as np

@click.command()
@click.option('--panel_file', type=click.Path(exists=True), required=True, help='Path to the panel file.')
@click.option('--bedfile', type=click.Path(), default=None, help='Optional BED file for exon definitions.')
@click.option('--expand', type=int, default=0, help='Expansion factor for defining regions.')
@click.option('--autoexons', is_flag=True, help='Use automated exon definition from the panel.')
def main(panel_file, bedfile, expand, autoexons):

    # Read input data
    panel_data = pd.read_table(panel_file)

    exons_bed = pd.DataFrame()
    exons_bed_provided = pd.DataFrame()

    if bedfile != 'None':
        print("Using the BED file provided")
        # Use BED file for exon definition
        exons_bed_provided = pd.read_table(bedfile, header=None, sep="\t")
        # Check if a name column is provided
        if exons_bed_provided.shape[1] > 3:
            exons_bed_provided = exons_bed_provided.iloc[:,:4]
            exons_bed_provided.columns = ["CHROM", "START", "END", "NAME"]
        else:
            exons_bed_provided = exons_bed_provided.iloc[:,:3]
            exons_bed_provided.columns = ["CHROM", "START", "END"]
            exons_bed_provided["NAME"] = None
    else:
        print("NOT using any BED file")

    if autoexons:
        print("Automated generation of exons from the panel")
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
    else:
        print("NO automated generation of exons from the panel")


    # Concatenate the previous two BED dataframes
    exons_bed = pd.concat((exons_bed, exons_bed_provided))

    if exons_bed.shape[0] == 0:
        raise ValueError("No regions to extend")

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

            if ("--exon_" not in region_name) and (expand > 0):
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


if __name__ == '__main__':
    main()
