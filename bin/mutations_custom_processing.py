#!/usr/local/bin/python


import pandas as pd
import numpy as np
import sys

from utils_impacts import *
from read_utils import custom_na_values



def customize_annotations(mutation_summary_file, custom_regions_file,
                            customized_mutations_output,
                            simple = True
                            ):
    """
    # TODO
    explain what this function does
    """
    # simple = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID'          , 'GENE', 'IMPACT'                                              , 'CONTEXT_MUT', 'CONTEXT']
    # rich   = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'STRAND', 'GENE', 'IMPACT', 'Feature', 'Protein_position', 'Amino_acids', 'CONTEXT_MUT', 'CONTEXT']
    mutation_summary = pd.read_csv(mutation_summary_file, sep = "\t", na_values = custom_na_values)
    print("all possible sites loaded")

    # alternatively we could take only the MUT_ID but
    # the one of the mutations is different from
    # the one coming from the panels
    custom_regions_df = pd.read_table(custom_regions_file)[["MUT_ID", "GENE", "IMPACT"]].drop_duplicates().reset_index(drop = True)
    custom_regions_df["MUT_ID"] = custom_regions_df["MUT_ID"].str.replace("/", ">", regex = 'False').str.replace("_", ":", n = 1)
    custom_regions_df.columns = ["MUT_ID", "SYMBOL_new", "Consequence_new"]

    mutation_summary_custom = mutation_summary.merge(custom_regions_df,
                                                        on = ["MUT_ID"],
                                                        suffixes = ('_old', ''),
                                                        how = 'left'
                                                        )

    mutations_summary_final = mutation_summary_custom[mutation_summary_custom["SYMBOL_new"].isna()].drop(["SYMBOL_new", "Consequence_new"], axis = 'columns')
    if mutations_summary_final.shape[0] == mutation_summary.shape[0]:
        print("No mutations to change")
        mutations_summary_final.to_csv(customized_mutations_output,
                                        header = True,
                                        index = False,
                                        sep = "\t")
        exit(0)
        # terminate the script here since there is nothing to change
    else:
        print(f"We are going to change: {mutations_summary_final.shape[0] - mutation_summary.shape[0]} mutations")


    updated_mutations_old = mutation_summary_custom[~(mutation_summary_custom["SYMBOL_new"].isna())]

    columns_to_null = ['Gene', 'Feature', 'Feature_type',
                        'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons',
                        'DISTANCE', 'STRAND', 'FLAGS',
                        'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'ENSP',

                        'canonical_Gene', 'canonical_Feature', 'canonical_Feature_type',
                        'canonical_cDNA_position', 'canonical_CDS_position', 'canonical_Protein_position',
                        'canonical_Amino_acids', 'canonical_Codons',
                        'canonical_DISTANCE', 'canonical_STRAND', 'canonical_FLAGS',
                        'canonical_SYMBOL_SOURCE', 'canonical_HGNC_ID', 'canonical_CANONICAL', 'canonical_ENSP',

                        ]

    columns_to_remain = ['MUT_ID', 'Location', 'Allele',
                            'Existing_variation', 'IMPACT',
                            'gnomADe_AF', 'gnomADe_AFR_AF', 'gnomADe_AMR_AF',
                            'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_NFE_AF',
                            'gnomADe_OTH_AF', 'gnomADe_SAS_AF', 'gnomADg_AF', 'gnomADg_AFR_AF',
                            'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF', 'gnomADg_EAS_AF',
                            'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF', 'gnomADg_OTH_AF',
                            'gnomADg_SAS_AF', 'CLIN_SIG', 'SOMATIC', 'PHENO',
                            'canonical_Location', 'canonical_Allele',
                            'canonical_Existing_variation', 'canonical_IMPACT',
                            'canonical_CLIN_SIG', 'canonical_SOMATIC', 'canonical_PHENO',
                            'CONTEXT_MUT', 'MUTTYPE', 'CONTEXT_MUT_SIGPRO',
        ]

#     columns_to_update = ['SYMBOL', 'canonical_SYMBOL',
#                             'Consequence', 'canonical_Consequence'
#                             ]

#     columns_to_redefine = ['canonical_Consequence_single', 'canonical_Consequence_broader', 'canonical_Protein_affecting',
#                             'Consequence_single', 'Consequence_broader', 'Protein_affecting'
#                             ]

    updated_mutations = updated_mutations_old[columns_to_remain].copy()
    updated_mutations[columns_to_null] = '-'


    updated_mutations['SYMBOL'] = updated_mutations_old["SYMBOL_new"]
    updated_mutations['canonical_SYMBOL'] = updated_mutations_old["SYMBOL_new"]
    updated_mutations['Consequence'] = updated_mutations_old["Consequence_new"]
    updated_mutations['canonical_Consequence'] = updated_mutations_old["Consequence_new"]

    updated_mutations['Consequence_single'] = updated_mutations['Consequence'].apply(most_deleterious_within_variant)
    updated_mutations['Consequence_broader'] = updated_mutations['Consequence_single'].apply(lambda x: GROUPING_DICT[x])
    updated_mutations['Protein_affecting'] = updated_mutations['Consequence_broader'].apply(lambda x: PROTEIN_AFFECTING_DICT[x])
    updated_mutations['canonical_Consequence_single'] = updated_mutations['canonical_Consequence'].apply(most_deleterious_within_variant)
    updated_mutations['canonical_Consequence_broader'] = updated_mutations['canonical_Consequence_single'].apply(lambda x: GROUPING_DICT[x])
    updated_mutations['canonical_Protein_affecting'] = updated_mutations['canonical_Consequence_broader'].apply(lambda x: PROTEIN_AFFECTING_DICT[x])


    mutations_summary_final = pd.concat((mutations_summary_final, updated_mutations)).reset_index(drop = True)
    mutations_summary_final.to_csv(customized_mutations_output,
                                        header = True,
                                        index = False,
                                        sep = "\t")


if __name__ == '__main__':
    # Input
    # VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
    mutations_file = sys.argv[1]

    custom_regions_file = sys.argv[2]

    # Output
    # all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"
    customized_output_file = sys.argv[3]

    customize_annotations(mutations_file, custom_regions_file, customized_output_file)
