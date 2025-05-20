#!/usr/bin/env python


import pandas as pd
import numpy as np
import sys

from read_utils import custom_na_values
muttype_conversion_map = {
                'G/A': 'C/T',
                'G/C': 'C/G',
                'G/T': 'C/A',
                'A/G': 'T/C',
                'A/T': 'T/A',
                'A/C': 'T/G',
            }


def load_chr_data_chunked(filepath, chrom, chunksize=1_000_000):
    """
    Loads data for a specific chromosome from a large VEP output file in chunks.

    Args:
        filepath (str): Path to the VEP output file.
        chrom (str): Chromosome to filter.
        chunksize (int): Number of rows per chunk.

    Returns:
        pd.DataFrame: Filtered DataFrame for the chromosome.
    """
    reader = pd.read_csv(filepath, sep="\t", na_values=custom_na_values, chunksize=chunksize, dtype={'CHROM': str})
    chr_data = []
    for chunk in reader:
        filtered = chunk[chunk["CHROM"] == chrom]
        if not filtered.empty:
            chr_data.append(filtered)
    return pd.concat(chr_data) if chr_data else pd.DataFrame()


def customize_panel_regions(VEP_output_file, custom_regions_file, customized_output_annotation_file,
                            simple = True
                            ):
    """
    Modifies annotations in a VEP output file based on custom genomic regions.

    - For each region in the custom regions file, identifies the corresponding slice
      in the VEP output.
    - Updates gene names and impact values for the region.
    - Saves both the modified annotation file and a record of added regions.

    Args:
        VEP_output_file (str): Path to the full VEP output file (TSV).
        custom_regions_file (str): Custom region definitions (tab-delimited).
        customized_output_annotation_file (str): Output file for updated annotations.
        simple (bool): If True, outputs simplified annotations; else adds more fields.
    """

    # simple = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID'          , 'GENE', 'IMPACT'                                              , 'CONTEXT_MUT', 'CONTEXT']
    # rich   = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'STRAND', 'GENE', 'IMPACT', 'Feature', 'Protein_position', 'Amino_acids', 'CONTEXT_MUT', 'CONTEXT']

    custom_regions_df = pd.read_table(custom_regions_file)
    added_regions_df = pd.DataFrame()
    current_chr = ""
    chr_data = pd.DataFrame()

    for _, row in custom_regions_df.iterrows():
        try:
            if row["CHROM"] != current_chr:
                current_chr = row["CHROM"]
                chr_data = load_chr_data_chunked(VEP_output_file, current_chr)

                print("Updating chromosome to:", current_chr)

            # Get start and end indices
            ind_start = np.where(chr_data["POS"] == row["START"])[0][0]
            ind_end = np.where(chr_data.iloc[ind_start:,:]["POS"] == row["END"])[0][-1]
            # gene = chr_data.iloc[ind_start]["GENE"]

            ind_end += ind_start

            upd_start = ind_start
            upd_end = ind_end

            # Extract hotspot data and modify gene names
            hotspot_data = chr_data.iloc[upd_start: upd_end + 1, :].copy().drop("IMPACT", axis = 'columns')
            hotspot_data["index_col"] = hotspot_data.index
            original_df_start = hotspot_data.index[0]
            original_df_end = hotspot_data.index[-1]
            hotspot_data["GENE"] = row["NAME"]

            # Split the string into individual entries
            entries = row["IMPACTFUL"].split(',')

            # Create a DataFrame
            impactful_df = pd.DataFrame(entries, columns = ["MUT_ID_pyr"])
            impactful_df["IMPACT"] = row["IMPACT"]

            all_pos_len = hotspot_data.shape[0]

            # convert MUT_ID to C>N and T>N only
            hotspot_data["MUT_ID_pyr"] = hotspot_data["MUT_ID"].replace(muttype_conversion_map, regex = True)
            hotspot_data = hotspot_data.merge(impactful_df,
                                                on = ["MUT_ID_pyr"],
                                                how = 'outer'
                                                )
            all_pos_len_after = hotspot_data.shape[0]

            # TODO add an error raise
            if all_pos_len != all_pos_len_after:
                print("Some of the mutations provided are not in the desired region")
                print(hotspot_data[hotspot_data["POS"].isna()]["MUT_ID_pyr"])
                hotspot_data = hotspot_data[~(hotspot_data["POS"].isna())].reset_index(drop = True)
                hotspot_data["POS"] = hotspot_data["POS"].astype(int)

            hotspot_data["IMPACT"] = hotspot_data["IMPACT"].fillna(row["NEUTRAL"])
            hotspot_data = hotspot_data.sort_values(by='index_col').drop("index_col", axis = 'columns')

            ## Insert modified rows back into the df
            if simple:
                chr_data.loc[original_df_start: original_df_end, ["GENE", "IMPACT"]] = hotspot_data[["GENE", "IMPACT"]].values
            else:
                print("Getting Feature to '-'")
                hotspot_data["Feature"] = '-'
                chr_data.loc[original_df_start: original_df_end, ["GENE", "IMPACT", "Feature"]] = hotspot_data[["GENE", "IMPACT", "Feature"]].values


            added_regions_df = pd.concat((added_regions_df, hotspot_data))
            print("Small region added:", row["NAME"])

        except Exception as e:
            print(f"Error processing row {row}: {e}")

    chr_data = chr_data.drop_duplicates(
        subset=['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'CONTEXT_MUT', 'CONTEXT', 'IMPACT'],
        keep='first'
    )
    chr_data.to_csv(customized_output_annotation_file, header=True, index=False, sep="\t")


    added_regions_df = added_regions_df.drop_duplicates(subset = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID',
                                                                    'GENE', 'CONTEXT_MUT', 'CONTEXT', 'IMPACT'],
                                                        keep = 'first')
    added_regions_df.to_csv('added_regions.tsv',
                                header = True,
                                index = False,
                                sep = "\t")


if __name__ == '__main__':
    # Input
    VEP_output_file = sys.argv[1]

    custom_regions_file = sys.argv[2]

    # Output
    customized_output_annotation_file = sys.argv[3]

    simple = eval(sys.argv[4])

    customize_panel_regions(VEP_output_file, custom_regions_file, customized_output_annotation_file,
                                    simple
                                    )

