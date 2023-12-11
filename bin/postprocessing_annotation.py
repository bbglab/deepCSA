#!/usr/bin/env python3


import pandas as pd
import numpy as np
import sys

from itertools import product
from bgreference import hg38



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





cb = dict(zip('ACGT', 'TGCA'))

def canonical_channels():

    subs = [''.join(z) for z in product('CT', 'ACGT') if z[0] != z[1]]
    flanks = [''.join(z) for z in product('ACGT', repeat=2)]
    contexts_tuples = [(a, b) for a, b in product(subs, flanks)]
    sorted_contexts_tuples = sorted(contexts_tuples, key=lambda x: (x[0], x[1]))
    sorted_contexts = [b[0] + a[0] + b[1] + '>' + a[1] for a, b in sorted_contexts_tuples]
    return sorted_contexts


def transform_context(chr_, pos, mut):
    ref, alt = tuple(mut.split('/'))
    ref_triplet = hg38(chr_, pos-1, size=3)
    if ref_triplet[1] not in ['C', 'T']:
        ref_triplet = ''.join(list(map(lambda x: cb[x], ref_triplet[::-1])))
        alt = cb[alt]
    return ref_triplet + '>' + alt





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


    'mature_miRNA_variant': 'non_coding_exon_region',
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
    all_consequences = impact_vep_string.split(",")
    all_consequences_ranks = map(lambda x: consequence_rank_dict_within[x], all_consequences)
    return rank_consequence_dict_within[min(all_consequences_ranks)]





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







def vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file, all_ = False):
    """
    # TODO
    explain what this function does
    """

    all_possible_sites = pd.read_csv(VEP_output_file, sep = "\t", header = 0)

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
        annotated_variants = annotated_variants.merge(annotated_variants_only_canonical, on = "MUT_ID", how = 'left')
        annotated_variants_only_canonical = annotated_variants_only_canonical[annotated_variants_only_canonical.columns[1:]].fillna('-')
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
    annotated_variants["CONTEXT_MUT"] = annotated_variants.apply(lambda x: transform_context(x["CHROM"], x["POS"], f'{x["REF"]}/{x["ALT"]}') if x["TYPE"] == "SNV" else "-", axis = 1)

#    annotated_variants_reduced = annotated_variants[['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'SYMBOL', 'IMPACT', 'CONTEXT']]
#    annotated_variants_reduced.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'IMPACT', 'CONTEXT_MUT']
    annotated_variants_columns = [x for x in annotated_variants.columns if x.replace("canonical_", "") not in ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'CHROM:POS', 'MUT'] ]
    annotated_variants_reduced = annotated_variants.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] ).reset_index(drop = True)
    annotated_variants_reduced = annotated_variants_reduced[ annotated_variants_columns ]

    annotated_variants_reduced.to_csv(all_possible_sites_annotated_file,
                                        header = True,
                                        index = False,
                                        sep = "\t")


if __name__ == '__main__':
    # Input
    #VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
    VEP_output_file = sys.argv[1]
    if len(sys.argv) >= 3:
        print("Using the provided value:", end = "\t")
        try:
            all_sep = eval(f"{sys.argv[3]}")
            print(all_sep)
        except:
            print("You should provide either True or False as the third argument.")
            exit(1)
    else:
        all_sep = False

    # Output
    #all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"
    all_possible_sites_annotated_file = sys.argv[2]



    vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file, all_sep)


