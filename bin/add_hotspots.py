#!/usr/bin/env python

import click
import json
import pandas as pd
import numpy as np

@click.command()
@click.option('--panel_file', type=click.Path(exists=True), required=True, help='Path to the panel file.')
@click.option('--autoexons', type=click.Path(exists=True), default=None, help='BED file for automated exon definitions.')
@click.option('--autodomains', type=click.Path(exists=True), default=None, help='BED file for automated domain definitions.')
@click.option('--custom', type=click.Path(exists=True), default=None, help='BED file for custom region definitions.')
def main(panel_file, autoexons, autodomains, custom):

    # Read input data
    panel_data = pd.read_table(panel_file)

    exons_bed = pd.DataFrame()

    for file in [custom, autoexons, autodomains]:
        if file is not None:
            print(f"Using the BED file: {file}")
            # Read BED file for exon or domain definitions
            bed_data = pd.read_table(file, header=None, sep="\t")
            # Check if a name column is provided
            if bed_data.shape[1] > 3:
                bed_data = bed_data.iloc[:,:4]
                bed_data.columns = ["CHROM", "START", "END", "NAME"]
            else:
                bed_data = bed_data.iloc[:,:3]
                bed_data.columns = ["CHROM", "START", "END"]
                bed_data["NAME"] = None
            exons_bed = pd.concat((exons_bed, bed_data))
        else:
            print("NOT using any BED file")

    if exons_bed.shape[0] == 0:
        raise ValueError("No regions to extend")

    ## Process panel data to create new genes per exon
    new_data = pd.DataFrame()
    hotspots_names = []

    current_chr = ""
    region_counters = {}  # Dictionary to track exon numbers for each gene
    for _ind, row in exons_bed.iterrows():
        try:
            if row["CHROM"] != current_chr:
                current_chr = row["CHROM"]
                chr_data = panel_data[panel_data["CHROM"] == current_chr]#.reset_index(drop = True)
                print("Updating chromosome to:", current_chr)

            # Get start and end indices
            start_matches = np.where(chr_data["POS"] == row["START"])[0]
            start_found = len(start_matches) > 0
            if start_found:
                ind_start = start_matches[0]
            else:
                # Use the first position greater than or equal to START, or the last if none
                start_matches = np.where(chr_data["POS"] >= row["START"])[0]
                if len(start_matches) > 0:
                    ind_start = start_matches[0]
                    start_found = True
                else:
                    print("No start match for row:", row)
                    continue

            # If you reach ths point, ind_start is defined and valid
            # Search for END position starting from ind_start
            search_end = chr_data.iloc[ind_start:,:]
            end_matches = np.where(search_end["POS"] == row["END"])[0]

            end_found = len(end_matches) > 0
            if end_found:
                ind_end = end_matches[-1]
            else:
                # Use the last position less than or equal to END, or the last available
                end_matches = np.where(search_end["POS"] <= row["END"])[-1]
                if len(end_matches) > 0:
                    ind_end = end_matches[-1]
                    end_found = True
                else:
                    print("No end match found for row:", row)
                    continue

            # if you reach this point, start and end are defined and valid
            gene = chr_data.iloc[ind_start]["GENE"]

            ind_end += ind_start


            # Initialize exon counter for the gene if not already present
            if gene not in region_counters:
                region_counters[gene] = 1

            # Determine the region name
            if pd.notna(row.get("NAME", None)):
                region_name = row["NAME"]
            else:
                region_number = region_counters[gene]
                region_name = f"{gene}--region{region_number}"

            # if ("--exon_" not in region_name) and (expand > 0):
            #     # Expand within limits of the same gene only if BED file is used
            #     upd_start = max(ind_start - expand * 3, 0)
            #     while upd_start < len(chr_data) and chr_data.iloc[upd_start]["GENE"] != gene:
            #         upd_start += 1

            #     upd_end = min(ind_end + expand * 3, len(chr_data) - 1)
            #     while upd_end >= 0 and chr_data.iloc[upd_end]["GENE"] != gene:
            #         upd_end -= 1

            # else:
            #     # No expansion for exon boundaries if BED file is not used
            upd_start = ind_start
            upd_end = ind_end

            # Extract hotspot data and modify gene names
            hotspot_data = chr_data.iloc[upd_start: upd_end + 1, :].copy()
            hotspot_data["GENE"] = region_name

            new_data = pd.concat((new_data, hotspot_data))
            hotspots_names.append(region_name)
            print("Small region added:", region_name)

            # Increment region counter for the gene (if name is generated)
            if row.get("NAME") is None:
                region_counters[gene] += 1

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
