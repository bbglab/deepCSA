#!/usr/local/bin/python


import pandas as pd
import numpy as np
import sys

from itertools import product
from bgreference import hg38, hg19, mm10

from utils_context import canonical_channels, transform_context

assembly_name2function = {"hg38": hg38,
                            "hg19": hg19,
                            "mm10": mm10}

def vartype(x,
            letters = ['A', 'T', 'C', 'G'],
            len_SV_lim = 100
            ):
    """
    Define the TYPE of a variant
    """
    if ">" in (x["REF"] + x["ALT"]) or "<" in (x["REF"] + x["ALT"]):
        return "SV"

    elif len(x["REF"]) > (len_SV_lim+1) or len(x["ALT"]) > (len_SV_lim+1) :
        return "SV"

    elif x["REF"] in letters and x["ALT"] in letters:
        return "SNV"

    elif len(x["REF"]) == len(x["ALT"]):
        return "MNV"

    elif x["REF"] == "-" or ( len(x["REF"]) == 1 and x["ALT"].startswith(x["REF"]) ):
        return "INSERTION"

    elif x["ALT"] == "-" or ( len(x["ALT"]) == 1 and x["REF"].startswith(x["ALT"]) ):
        return "DELETION"

    return "COMPLEX"





GROUPING_DICT = {

    'transcript_ablation': 'nonsense',
    'splice_acceptor_variant': 'nonsense',
    'splice_donor_variant': 'nonsense',
    'stop_gained': 'nonsense',
    'frameshift_variant': 'nonsense',
    'stop_lost': 'nonsense',
    'start_lost': 'nonsense',

    'missense_variant': 'missense',
    'inframe_insertion': 'missense',
    'inframe_deletion': 'missense',

    'splice_donor_variant': 'essential_splice',
    'splice_acceptor_variant': 'essential_splice',


    'splice_region_variant': 'splice_region',
    'splice_donor_5th_base_variant': 'splice_region',
    'splice_donor_region_variant': 'splice_region',
    'splice_polypyrimidine_tract_variant': 'splice_region',

    'synonymous_variant': 'synonymous',
    'incomplete_terminal_codon_variant': 'synonymous',
    'start_retained_variant': 'synonymous',
    'stop_retained_variant': 'synonymous',

    'protein_altering_variant' : 'protein_altering_variant', ##
    'transcript_amplification' : 'transcript_amplification', ##
    'coding_sequence_variant': 'coding_sequence_variant', ##


    'mature_miRNA_variant': 'non_coding_exon_region', # TODO fix this
    '5_prime_UTR_variant': 'non_coding_exon_region',
    '3_prime_UTR_variant': 'non_coding_exon_region',
    'non_coding_transcript_exon_variant': 'non_coding_exon_region',

    'NMD_transcript_variant': 'non_coding_exon_region',

    'intron_variant': 'intron_variant',

    'non_coding_transcript_variant' : 'non_coding_transcript_variant',
    'upstream_gene_variant': 'non_genic_variant',
    'downstream_gene_variant': 'non_genic_variant',
    'TFBS_ablation': 'non_genic_variant',
    'TFBS_amplification': 'non_genic_variant',
    'TF_binding_site_variant': 'non_genic_variant',
    'regulatory_region_ablation': 'non_genic_variant',
    'regulatory_region_amplification': 'non_genic_variant',
    'feature_elongation': 'non_genic_variant',
    'regulatory_region_variant': 'non_genic_variant',
    'feature_truncation': 'non_genic_variant',
    'intergenic_variant': 'non_genic_variant',
    '-'  : '-'

}

PROTEIN_AFFECTING_DICT = {

    'nonsense' : 'protein_affecting',
    'missense' : 'protein_affecting',
    'essential_splice' : 'protein_affecting',
    'splice_region' : 'ambiguous',
    'synonymous' : 'non_protein_affecting',

    'protein_altering_variant' : 'protein_affecting',
    'transcript_amplification' : 'protein_affecting',
    'coding_sequence_variant' : 'ambiguous',

    'non_coding_exon_region' : 'non_protein_affecting',
    'intron_variant' : 'non_protein_affecting',
    'non_coding_transcript_variant' : 'non_protein_affecting',
    'non_genic_variant' : 'non_protein_affecting',
    '-'  : '-',
}


CONSEQUENCES_LIST = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'splice_donor_5th_base_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant',
    '-'
]




CONSEQUENCES_LIST_WITHIN = [
    'NMD_transcript_variant',
    'non_coding_transcript_exon_variant',
    'non_coding_transcript_variant',

    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'splice_donor_5th_base_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'intron_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant',
    '-'
]



consequence_rank_dict = {consequence : rank for rank, consequence in enumerate(CONSEQUENCES_LIST)}
rank_consequence_dict = {rank : consequence for rank, consequence in enumerate(CONSEQUENCES_LIST)}



consequence_rank_dict_within = {consequence : rank for rank, consequence in enumerate(CONSEQUENCES_LIST_WITHIN)}
rank_consequence_dict_within = {rank : consequence for rank, consequence in enumerate(CONSEQUENCES_LIST_WITHIN)}

def most_deleterious_within_variant(impact_vep_string):
    """
    to be used when summarizing the different consquences assigned to a same variable in the same transcript
    here we change for example the relevance of NMD_transcript_variant, since we do not want it to make it very damaging
    """
    # TODO: revise if we need to have a try and except or it is better to make sure that the consequence
    # dictionary and ranks correspond to the correct ensembl version?
    try :
        all_consequences = impact_vep_string.split(",")
        all_consequences_ranks = map(lambda x: consequence_rank_dict_within[x], all_consequences)
        return rank_consequence_dict_within[min(all_consequences_ranks)]
    except:
        return '-'


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

    returned_df.columns = ["MUT_ID"] + [f"canonical_{x}" for x in returned_df.columns[1:]]
    # we return the dataframe with all the original columns of the VEP file
    return returned_df






def vep2summarizedannotation_panel(VEP_output_file, all_possible_sites_annotated_file, assembly = 'hg38'):
    """
    # TODO
    explain what this function does
    """
    all_possible_sites = pd.read_csv(VEP_output_file, sep = "\t",
                                        header = None)
    print("all possible sites loaded")
    all_possible_sites.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'Consequence', 'SYMBOL']

    annotated_variants = VEP_annotation_to_single_row(all_possible_sites, keep_genes= True)
    del all_possible_sites
    print("VEP to single row working")


    # annotated_variants_only_canonical = VEP_annotation_to_single_row_only_canonical(all_possible_sites, keep_genes= True)
    # if annotated_variants_only_canonical is not None:
    #     annotated_variants = annotated_variants.merge(annotated_variants_only_canonical, on = "MUT_ID", how = 'left')
    #     annotated_variants_only_canonical = annotated_variants_only_canonical[annotated_variants_only_canonical.columns[1:]].fillna('-')
    #     annotated_variants['canonical_Consequence_single'] = annotated_variants['canonical_Consequence'].apply(most_deleterious_within_variant)
    #     annotated_variants['canonical_Consequence_broader'] = annotated_variants['canonical_Consequence_single'].apply(lambda x: GROUPING_DICT[x])
    #     annotated_variants['canonical_Protein_affecting'] = annotated_variants['canonical_Consequence_broader'].apply(lambda x: PROTEIN_AFFECTING_DICT[x])



    # TODO: agree on a consensus for these broader consequence types
    # add a new column containing a broader  consequence per variant
    annotated_variants['Consequence_single'] = annotated_variants['Consequence'].apply(most_deleterious_within_variant)
    annotated_variants['Consequence_broader'] = annotated_variants['Consequence_single'].apply(lambda x: GROUPING_DICT[x])
    print("Consequence to IMPACT working")

    # add context type to all SNVs
    # remove context from the other substitution types
    chosen_assembly = assembly_name2function[assembly]
    annotated_variants["CONTEXT_MUT"] = annotated_variants.apply(lambda x: transform_context(x["CHROM"], x["POS"], f'{x["REF"]}/{x["ALT"]}', chosen_assembly) , axis = 1)
    print("Context added")

    annotated_variants["CONTEXT"] = annotated_variants["CONTEXT_MUT"].apply(lambda x: x[:3])

    annotated_variants_reduced = annotated_variants[['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'SYMBOL', 'Consequence_broader', 'CONTEXT_MUT', 'CONTEXT']]
    annotated_variants_reduced.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'IMPACT', 'CONTEXT_MUT', 'CONTEXT']
    annotated_variants_reduced = annotated_variants_reduced.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] )
    print("Annotation sorted")


    annotated_variants_reduced.to_csv(all_possible_sites_annotated_file,
                                        header = True,
                                        index = False,
                                        sep = "\t")


if __name__ == '__main__':
    # Input
    # VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
    VEP_output_file = sys.argv[1]

    assembly = sys.argv[2]

    # Output
    # all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"
    all_possible_sites_annotated_file = sys.argv[3]




    vep2summarizedannotation_panel(VEP_output_file, all_possible_sites_annotated_file, assembly)


