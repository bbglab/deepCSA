#!/usr/bin/env python

import sys
import click
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


muttype_conversion_map = {
                'G>A': 'C>T',
                'G>C': 'C>G',
                'G>T': 'C>A',
                'A>G': 'T>C',
                'A>T': 'T>A',
                'A>C': 'T>G',
            }


def get_canonical_mutid(mutid):
    elements__ = mutid.split("_")
    mutation_change = elements__[-1]
    upd_mutation_change = muttype_conversion_map.get(mutation_change, mutation_change)
    return "_".join(elements__[:-1] + [upd_mutation_change])





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
                                all_ = False, assembly = 'hg38',
                                gnomad_af_threshold = 0.001):
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
        all_possible_sites[["REF", "ALT"]] = all_possible_sites["MUT"].str.split(">", n = 1, expand = True)
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
        _annotated_variants_only_canonical_cleaned = annotated_variants_only_canonical.drop(gnomad_repeated_columns, axis = 'columns')

        try :
            annotated_variants_only_canonical_cleaned = _annotated_variants_only_canonical_cleaned.drop(["canonical_Existing_variation", "canonical_PHENO", "canonical_SOMATIC"], axis = 'columns')
        except:
            annotated_variants_only_canonical_cleaned = _annotated_variants_only_canonical_cleaned

        annotated_variants = annotated_variants.merge(annotated_variants_only_canonical_cleaned, on = "MUT_ID", how = 'left')
        annotated_variants['canonical_Consequence_single'] = annotated_variants['canonical_Consequence'].apply(most_deleterious_within_variant)
        annotated_variants['canonical_Consequence_broader'] = annotated_variants['canonical_Consequence_single'].apply(lambda x: GROUPING_DICT[x])
        annotated_variants['canonical_Protein_affecting'] = annotated_variants['canonical_Consequence_broader'].apply(lambda x: PROTEIN_AFFECTING_DICT[x])
        del annotated_variants_only_canonical
        del _annotated_variants_only_canonical_cleaned
        del annotated_variants_only_canonical_cleaned


    # TODO: agree on a consensus for these broader consequence types
    # add a new column containing a broader consequence per variant
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

    if any("gnomAD" in x for x in annotated_variants_reduced.columns):
        gnomad_columns = [x for x in annotated_variants_reduced.columns if 'gnomAD' in x ]
        annotated_variants_reduced[gnomad_columns] = annotated_variants_reduced[gnomad_columns].replace("-", 0).astype(float)

        # add a column to flag the variants considered to be SNPs based on gnomad information
        annotated_variants_reduced["gnomAD_SNP"] = (annotated_variants_reduced['gnomADe_AF'] > gnomad_af_threshold) | (annotated_variants_reduced['gnomADg_AF'] > gnomad_af_threshold)
        annotated_variants_columns += ["gnomAD_SNP"]


    if hotspots_file is not None:
        hotspots_def_df = pd.read_table(hotspots_file, header = 0, sep = '\t').drop_duplicates()
        new_hostpot_columns = [x for x in hotspots_def_df.columns if x not in ['CHROM', 'POS', "MUTTYPE"] ]
        annotated_variants_reduced = annotated_variants_reduced.merge(hotspots_def_df, on = ['CHROM', 'POS', "MUTTYPE"], how = 'left')
        annotated_variants_columns += new_hostpot_columns

    # Check for duplicates
    duplicates = annotated_variants_reduced[annotated_variants_columns][
        annotated_variants_reduced[annotated_variants_columns].duplicated(keep=False)
    ]
    if not duplicates.empty:
        print(f"Found {duplicates.shape[0]} duplicated rows:")
        print(duplicates)
        sys.exit("Execution stopped due to duplicate rows.")

    annotated_variants_reduced = annotated_variants_reduced[annotated_variants_columns]

    annotated_variants_reduced["MUT_ID_pyr"] = annotated_variants_reduced["MUT_ID"].apply(lambda x : get_canonical_mutid(x))

    annotated_variants_reduced.to_csv(all_possible_sites_annotated_file,
                                        header = True,
                                        index = False,
                                        sep = "\t")




@click.command()
@click.argument('vep_output_file', type=click.Path(exists=True))
@click.argument('all_possible_sites_annotated_file', type=click.Path())
@click.option('--assembly-name', default='hg38', show_default=True, type=click.Choice(['hg38', 'hg19', 'mm10', 'mm39']), help='Reference genome assembly name')
@click.option('--all-underscore', is_flag=True, default=False, show_default=True, help='Whether to use _ to separate all parts of MUT_ID (default: False)')
@click.option('--hotspots-annotation-file', default=None, type=click.Path(exists=False), help='Path to hotspots annotation file')
@click.option('--gnomad-af-threshold', default=0.001, show_default=True, type=float, help='gnomAD allele frequency threshold')
def main(vep_output_file, all_possible_sites_annotated_file, assembly_name, all_underscore, hotspots_annotation_file, gnomad_af_threshold):
    """Summarize VEP annotation with optional hotspots and gnomAD AF threshold."""
    vep2summarizedannotation(
        vep_output_file,
        all_possible_sites_annotated_file,
        hotspots_annotation_file,
        all_underscore,
        assembly_name,
        gnomad_af_threshold
    )

if __name__ == '__main__':
    main()

