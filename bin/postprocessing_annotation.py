#!/usr/local/bin/python

import sys
import pandas as pd

from bgreference import hg38, hg19, mm10, mm39

from utils import vartype
from utils_context import transform_context
from utils_impacts import *
from read_utils import custom_na_values

assembly_name2function = {"hg38": hg38,
                            "hg19": hg19,
                            "mm10": mm10,
                            "mm39": mm39}







def VEP_annotation_to_single_row(df_annotation):
    """
    Process Ensembl VEP output to get a single consequence per gene
    Select always the most deleterious
    """
    # update the first column name to ID
    df_annotation.columns = ["MUT_ID"] + list(df_annotation.columns[1:])

    # select a subset of the columns
    if "CANONICAL" in df_annotation.columns:

        df_annotation_small = df_annotation[['MUT_ID',
                                                'Consequence',
                                                'SYMBOL',
                                                'CANONICAL']]
    else:
        df_annotation_small = df_annotation[['MUT_ID',
                                                'Consequence',
                                                'SYMBOL']]


    df_annotation_small = df_annotation_small.drop_duplicates()


    # add a new column containing a single consequence per variant
    df_annotation_small["Consequence_single"] = df_annotation_small["Consequence"].apply(most_deleterious_within_variant)

    # assign a numerical value to each consequence according to its rank in damaging consequence
    df_annotation_small["NUM_Consequence"] = df_annotation_small["Consequence_single"].map(consequence_rank_dict)

    # sort the data by the ID and the IMPACT in numerical format.
    # The smaller the value of NUM_IMPACT the bigger the impact.
    if "CANONICAL" in df_annotation_small.columns:
        print("prioritizing canonical")
        df_annotation_small_sorted = df_annotation_small.sort_values(by = ['MUT_ID', "NUM_Consequence", "CANONICAL"],
                                                                        ascending = (True, True, False)
                                                                    )
    else:
        df_annotation_small_sorted = df_annotation_small.sort_values(by = ['MUT_ID', "NUM_Consequence"],
                                                                        ascending = (True, True)
                                                                    )


    df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID'],
                                                                                    keep='first')
    returned_df = df_annotation.iloc[df_annotation_small_highest_impact.index.values,:].copy()

    # TODO see if we can reduce this code by outputting the variable directly
    returned_df = returned_df.reset_index(drop = True)

    # we return the dataframe with all the original columns of the VEP file
    return returned_df



def VEP_annotation_to_single_row_only_canonical(df_annotation):
    """
    Process Ensembl VEP output to get a single consequence per gene
    Use only the information in the canonical transcript!!
    (select always the most deleterious)
    """
    # update the first column name to ID
    df_annotation.columns = ["MUT_ID"] + list(df_annotation.columns[1:])

    # select a subset of the columns
    if "CANONICAL" not in df_annotation.columns:
        return None

    df_annotation = df_annotation[df_annotation["CANONICAL"] == "YES"].reset_index(drop = True)
    df_annotation_small = df_annotation[['MUT_ID',
                                            'Consequence',
                                            'SYMBOL']]


    df_annotation_small = df_annotation_small.drop_duplicates()

    # add a new column containing a single consequence per variant
    df_annotation_small["Consequence_single"] = df_annotation_small["Consequence"].apply(most_deleterious_within_variant)

    # assign a numerical value to each consequence according to its rank in damaging consequence
    df_annotation_small["NUM_Consequence"] = df_annotation_small["Consequence_single"].map(consequence_rank_dict)

    # sort the data by the ID and the IMPACT in numerical format.
    # The smaller the value of NUM_IMPACT the bigger the impact.
    df_annotation_small_sorted = df_annotation_small.sort_values(by = ['MUT_ID', "NUM_Consequence"],
                                                                    ascending = (True, True)
                                                                )

    df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID'],
                                                                                    keep='first')
    returned_df = df_annotation.iloc[df_annotation_small_highest_impact.index.values,:].copy()
    returned_df = returned_df.reset_index(drop = True)

    returned_df.columns = ["MUT_ID"] + [f"canonical_{x}" for x in returned_df.columns[1:]]
    # we return the dataframe with all the original columns of the VEP file
    return returned_df







def vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file,
                             hotspots_file = None,
                             all_ = False, assembly = 'hg38'):
    """
    # TODO
    explain what this function does
    """

    all_possible_sites = pd.read_csv(VEP_output_file, sep = "\t", header = 0, na_values = custom_na_values)

    if all_ :
        all_possible_sites[["CHROM", "POS", "MUT" ]] = all_possible_sites.iloc[:,0].str.split("_", expand = True)
        all_possible_sites[["REF", "ALT"]] = all_possible_sites["MUT"].str.split("/", expand = True)
        all_possible_sites["POS"] = all_possible_sites["POS"].astype(int)
    else:
        all_possible_sites[["CHROM:POS", "MUT" ]] = all_possible_sites.iloc[:,0].str.split("_", expand = True)
        all_possible_sites[["CHROM", "POS"]] = all_possible_sites["CHROM:POS"].str.split(":", expand = True)
#        mut_split  = all_possible_sites["MUT"].str.split(">", n = 1, expand = True)
        all_possible_sites[["REF", "ALT"]] = all_possible_sites["MUT"].str.split(">", n = 1, expand = True)
#        print(all_possible_sites.head())
        all_possible_sites["POS"] = all_possible_sites["POS"].astype(int)

    # TODO: Is it robust enough to use columns names here?
#    all_possible_sites = all_possible_sites[['#Uploaded_variation', 'Consequence', 'SYMBOL',
#                                                'CHROM', 'POS', 'REF', 'ALT', 'MUT']]

    all_possible_sites["TYPE"] = all_possible_sites[["REF", "ALT"]].apply(vartype, axis = 1)

    # if the annotation already contains a single consequence per gene this function does not do much
    # but if it contains multiple variants per gene it keeps only the most deleterious
    annotated_variants = VEP_annotation_to_single_row(all_possible_sites)
    annotated_variants_only_canonical = VEP_annotation_to_single_row_only_canonical(all_possible_sites)
    if annotated_variants_only_canonical is not None:
        gnomad_repeated_columns = [x for x in annotated_variants_only_canonical.columns if 'gnomAD' in x ]
        annotated_variants_only_canonical_cleaned = annotated_variants_only_canonical.drop(gnomad_repeated_columns, axis = 'columns')
        annotated_variants = annotated_variants.merge(annotated_variants_only_canonical_cleaned, on = "MUT_ID", how = 'left')
        annotated_variants['canonical_Consequence_single'] = annotated_variants['canonical_Consequence'].apply(most_deleterious_within_variant)
        annotated_variants['canonical_Consequence_broader'] = annotated_variants['canonical_Consequence_single'].apply(lambda x: GROUPING_DICT[x])
        annotated_variants['canonical_Protein_affecting'] = annotated_variants['canonical_Consequence_broader'].apply(lambda x: PROTEIN_AFFECTING_DICT[x])



    # TODO: agree on a consensus for these broader consequence types
    # add a new column containing a broader  consequence per variant
    annotated_variants['Consequence_single'] = annotated_variants['Consequence'].apply(most_deleterious_within_variant)
    annotated_variants['Consequence_broader'] = annotated_variants['Consequence_single'].apply(lambda x: GROUPING_DICT[x])
    annotated_variants['Protein_affecting'] = annotated_variants['Consequence_broader'].apply(lambda x: PROTEIN_AFFECTING_DICT[x])


    # add context type to all SNVs
    # remove context from the other substitution types
    genome_func = assembly_name2function[assembly]
    annotated_variants["CONTEXT_MUT"] = annotated_variants.apply(lambda x: transform_context(x["CHROM"], x["POS"], f'{x["REF"]}/{x["ALT"]}', genome_func) if x["TYPE"] == "SNV" else "-", axis = 1)
    annotated_variants["MUTTYPE"] = annotated_variants["CONTEXT_MUT"].apply(lambda x: f"{x[1]}>{x[4]}" if x != '-' else '-')
    annotated_variants["CONTEXT_MUT_SIGPRO"] = annotated_variants["CONTEXT_MUT"].apply(lambda x: f"{x[0]}[{x[1]}>{x[4]}]{x[2]}" if x != '-' else '-')

    annotated_variants_columns = [x for x in annotated_variants.columns if x.replace("canonical_", "") not in ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'CHROM:POS', 'MUT'] ]
    annotated_variants_reduced = annotated_variants.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] ).reset_index(drop = True)


    if hotspots_file is not None:
        hotspots_def_df = pd.read_table(hotspots_file, header = 0, sep = '\t')
        new_hostpot_columns = [x for x in hotspots_def_df.columns if x not in ['CHROM', 'POS', "MUTTYPE"] ]
        annotated_variants_reduced = annotated_variants_reduced.merge(hotspots_def_df, on = ['CHROM', 'POS', "MUTTYPE"], how = 'left')
        annotated_variants_columns += new_hostpot_columns

    annotated_variants_reduced = annotated_variants_reduced[ annotated_variants_columns ]

    annotated_variants_reduced.to_csv(all_possible_sites_annotated_file,
                                        header = True,
                                        index = False,
                                        sep = "\t")


if __name__ == '__main__':
    # Input
    #VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
    VEP_output_file = sys.argv[1]

    # Output
    #all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"
    all_possible_sites_annotated_file = sys.argv[2]

    try:
        assembly_name = sys.argv[3]
        if assembly_name not in ["hg38", "hg19", "mm10", "mm39"]:
            print("invalid assembly name")
            exit(1)
    except:
        print("No assembly name provided, using hg38 as default.")
        assembly_name = 'hg38'


    if len(sys.argv) > 4:
        print("Using the provided value:", end = "\t")
        try:
            all_sep = eval(f"{sys.argv[4]}")
            print(all_sep)
        except:
            print("You should provide either True or False as the fourth argument.")
            exit(1)
    else:
        all_sep = False


    if len(sys.argv) > 5:
        print("Using the provided value:", end = "\t")
        try:
            hotspots_annotation_file = f"{sys.argv[5]}"
            print(hotspots_annotation_file)
            hotspots_annotation_df = pd.read_table(hotspots_annotation_file)
            print(hotspots_annotation_df.head())
        except:
            print("You should provide the path to the hotspots file as the fifth argument.")
            exit(1)
    else:
        hotspots_annotation_file = None


    vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file, hotspots_annotation_file, all_sep, assembly_name)

