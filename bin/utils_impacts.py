GROUPING_DICT = {

    'transcript_ablation': 'nonsense',
    'stop_gained': 'nonsense',
    'frameshift_variant': 'nonsense',
    'stop_lost': 'nonsense',
    'start_lost': 'nonsense',

    'missense_variant': 'missense',
    'missense': 'missense',                                     ## customly added
    'inframe_insertion': 'missense',
    'inframe_deletion': 'missense',

    'splice_donor_variant': 'essential_splice',
    'splice_acceptor_variant': 'essential_splice',
    'splice_donor_5th_base_variant': 'essential_splice',


    'splice_region_variant': 'splice_region_variant',
    'splice_donor_region_variant': 'splice_region_variant',
    'splice_polypyrimidine_tract_variant': 'splice_region_variant',

    'synonymous_variant': 'synonymous',
    'synonymous': 'synonymous',                                 ## customly added
    'incomplete_terminal_codon_variant': 'synonymous',
    'start_retained_variant': 'synonymous',
    'stop_retained_variant': 'synonymous',

    'protein_altering_variant' : 'protein_altering_variant', ##
    'transcript_amplification' : 'transcript_amplification', ##
    'coding_sequence_variant': 'coding_sequence_variant', ##


    '5_prime_UTR_variant': 'non_coding_exon_region',
    '3_prime_UTR_variant': 'non_coding_exon_region',
    'non_coding_transcript_exon_variant': 'non_coding_exon_region',

    'NMD_transcript_variant': 'non_coding_exon_region',

    'intron_variant': 'intron_variant',

    'non_coding_transcript_variant' : 'non_coding_transcript_variant',
    'mature_miRNA_variant': 'non_coding_transcript_variant',

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

    'synonymous' : 'non_protein_affecting',

    'protein_altering_variant' : 'protein_affecting',
    'transcript_amplification' : 'protein_affecting',
    'coding_sequence_variant' : 'ambiguous',

    'splice_region_variant': 'ambiguous',
    'splice_donor_region_variant': 'ambiguous',
    'splice_polypyrimidine_tract_variant': 'ambiguous',

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
    'missense',                                 ## customly added
    'protein_altering_variant',
    'splice_region_variant',
    'splice_donor_5th_base_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'synonymous',                               ## customly added

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
    'mature_miRNA_variant',

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
    'missense',                                 ## customly added
    'protein_altering_variant',
    'splice_region_variant',
    'splice_donor_5th_base_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'synonymous',                               ## customly added
    'coding_sequence_variant',
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


broadimpact_grouping_dict = {
        "missense": ["missense"],
        "nonsense": ["nonsense"],
        "essential_splice": ["essential_splice"],
        "splice_region_variant": ["splice_region_variant"],
        "truncating": ["nonsense", "essential_splice"],
        "essential_splice_plus": ["essential_splice", "splice_region_variant"],
        "truncating_plus": ["nonsense", "essential_splice", "splice_region_variant"],
        "nonsynonymous_splice": ["missense", "nonsense", "essential_splice"]
    }