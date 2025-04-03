#!/usr/bin/env python


import pandas as pd
import numpy as np
import sys

from itertools import product
from bgreference import hg38, hg19, mm10, mm39

from utils import vartype
from utils_context import canonical_channels, transform_context
from utils_impacts import *
from read_utils import custom_na_values

assembly_name2function = {"hg38": hg38,
                            "hg19": hg19,
                            "mm10": mm10,
                            "mm39": mm39
                            }


def VEP_annotation_to_single_row(df_annotation, keep_genes = False):
    """
    Process Ensembl VEP output to get a single consequence per gene
    Select always the most deleterious
    """

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

    if keep_genes:
        df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID'
                                                                                                ,'SYMBOL']
                                                                                            , keep='first'
                                                                                            )
    else:
        df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID']
                                                                                            , keep='first'
                                                                                            )

    returned_df = df_annotation.iloc[df_annotation_small_highest_impact.index.values,:].copy()

    # TODO see if we can reduce this code by outputting the variable directly
    returned_df = returned_df.reset_index(drop = True)

    # returned_df["Consequence"] = returned_df["Consequence"].apply(get_single_annotation)

    # we return the dataframe with all the original columns of the VEP file
    return returned_df


def safe_transform_context(row):
    if pd.isna(row["POS"]) or pd.isna(row["CHROM"]) or pd.isna(row["REF"]) or pd.isna(row["ALT"]):
        return "UNKNOWN"
    return transform_context(row["CHROM"], row["POS"], f'{row["REF"]}/{row["ALT"]}', chosen_assembly)


def VEP_annotation_to_single_row_only_canonical(df_annotation, keep_genes = False):
    """
    Process Ensembl VEP output to get a single consequence per gene
    Use only the information in the canonical transcript!!
    (select always the most deleterious)
    """

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

    if keep_genes:
        df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID'
                                                                                                ,'SYMBOL']
                                                                                            , keep='first'
                                                                                            )
    else:
        df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID']
                                                                                            , keep='first'
                                                                                            )

    returned_df = df_annotation.iloc[df_annotation_small_highest_impact.index.values,:].copy()
    returned_df = returned_df.reset_index(drop = True)

    returned_df.columns = [f"canonical_{x}" for x in returned_df.columns]
    # we return the dataframe with all the original columns of the VEP file
    return returned_df



def process_chunk(chunk, chosen_assembly, using_canonical):
    print("all possible sites loaded")
    chunk.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'Feature', 'Consequence', 'Protein_position', 'Amino_acids', 'STRAND', 'SYMBOL', 'CANONICAL', 'ENSP']

    if using_canonical:
        annotated_variants = VEP_annotation_to_single_row_only_canonical(chunk, keep_genes= True)
        if annotated_variants is not None:
            annotated_variants.columns = [ x.replace("canonical_", "") for x in annotated_variants.columns]
            print("Using only canonical transcript annotations for the panel")
        else:
            annotated_variants = VEP_annotation_to_single_row(chunk, keep_genes= True)
            print("CANONICAL was not available in the panel annotation.")
            print("Using most deleterious consequence for the panel")
    else:
        annotated_variants = VEP_annotation_to_single_row(chunk, keep_genes= True)
        print("Using most deleterious consequence for the panel")

    del chunk
    gc.collect()
    annotated_variants[annotated_variants.columns[1:]] = annotated_variants[annotated_variants.columns[1:]].fillna('-')
    print("VEP to single row working")


    # TODO: agree on a consensus for these broader consequence types
    # add a new column containing a broader  consequence per variant
    annotated_variants['Consequence_single'] = annotated_variants['Consequence'].apply(most_deleterious_within_variant)
    annotated_variants['Consequence_broader'] = annotated_variants['Consequence_single'].apply(lambda x: GROUPING_DICT[x])
    print("Consequence to IMPACT working")

    # add context type to all SNVs
    # remove context from the other substitution types
    
    annotated_variants["CONTEXT_MUT"] = annotated_variants.apply(lambda row: safe_transform_context(row, chosen_assembly), axis=1)
    print("Context added")

    annotated_variants["CONTEXT"] = annotated_variants["CONTEXT_MUT"].apply(lambda x: x[:3])

    annotated_variants_reduced = annotated_variants[['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'STRAND', 'SYMBOL', 'Consequence_broader', 'Feature', 'Protein_position', 'Amino_acids', 'CONTEXT_MUT', 'CONTEXT']]
    annotated_variants_reduced.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'STRAND', 'GENE', 'IMPACT', 'Feature', 'Protein_position', 'Amino_acids', 'CONTEXT_MUT', 'CONTEXT']
    annotated_variants_reduced = annotated_variants_reduced.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] )
    print("Annotation sorted")

    return annotated_variants_reduced

def vep2summarizedannotation_panel(VEP_output_file, all_possible_sites_annotated_file,
                                    assembly = 'hg38',
                                    using_canonical = True
                                    ):
    """
    # TODO
    explain what this function does
    """
    chosen_assembly = assembly_name2function[assembly]
    chunk_size = 100000

    reader = pd.read_csv(VEP_output_file, sep="\t", header=None, na_values=custom_na_values, chunksize=chunk_size)
    
    with open(f"{all_possible_sites_annotated_file}_rich.tsv", "w") as rich_out_file, \
         open(f"{all_possible_sites_annotated_file}.tsv", "w") as simple_out_file:
        
        for i, chunk in enumerate(reader):
            processed_chunk = process_chunk(chunk, chosen_assembly, using_canonical)
            
            rich_out_file.write(processed_chunk.to_csv(header=(i == 0), index=False, sep="\t"))
            simple_out_file.write(processed_chunk[['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'IMPACT', 'CONTEXT_MUT', 'CONTEXT']]
                                  .to_csv(header=(i == 0), index=False, sep="\t"))
            
            del processed_chunk
            gc.collect()


if __name__ == '__main__':
    # Input
    # VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
    VEP_output_file = sys.argv[1]

    assembly = sys.argv[2]

    # Output
    # all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"
    all_possible_sites_annotated_file = sys.argv[3]

    only_canonical = bool(sys.argv[4])

    vep2summarizedannotation_panel(VEP_output_file, all_possible_sites_annotated_file, assembly, only_canonical)


