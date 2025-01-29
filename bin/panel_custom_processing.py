#!/usr/local/bin/python


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


def customize_panel_regions(VEP_output_file, custom_regions_file, customized_output_annotation_file,
                            simple = True
                            ):
    """
    # TODO
    explain what this function does
    """
    # simple = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID'          , 'GENE', 'IMPACT'                                              , 'CONTEXT_MUT', 'CONTEXT']
    # rich   = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'STRAND', 'GENE', 'IMPACT', 'Feature', 'Protein_position', 'Amino_acids', 'CONTEXT_MUT', 'CONTEXT']
    all_possible_sites = pd.read_csv(VEP_output_file, sep = "\t",
                                        na_values = custom_na_values)
    print("all possible sites loaded")

    custom_regions_df = pd.read_table(custom_regions_file)

    added_regions_df = pd.DataFrame()

    current_chr = ""
    for ind, row in custom_regions_df.iterrows():
        try:
            if row["CHROM"] != current_chr:
                current_chr = row["CHROM"]
                chr_data = all_possible_sites[all_possible_sites["CHROM"] == current_chr]
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
                all_possible_sites.loc[original_df_start: original_df_end, ["GENE", "IMPACT"]] = hotspot_data[["GENE", "IMPACT"]].values
            else:
                print("Getting Feature to '-'")
                hotspot_data["Feature"] = '-'
                all_possible_sites.loc[original_df_start: original_df_end, ["GENE", "IMPACT", "Feature"]] = hotspot_data[["GENE", "IMPACT", "Feature"]].values

            added_regions_df = pd.concat((added_regions_df, hotspot_data))
            print("Small region added:", row["NAME"])

        except Exception as e:
            print(f"Error processing row {row}: {e}")

    all_possible_sites = all_possible_sites.drop_duplicates(subset = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID',
                                                                    'GENE', 'CONTEXT_MUT', 'CONTEXT', 'IMPACT'],
                                                            keep = 'first')
    all_possible_sites.to_csv(customized_output_annotation_file,
                                        header = True,
                                        index = False,
                                        sep = "\t")

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

