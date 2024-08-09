# from itertools import product
import json
import pandas as pd

custom_na_values_list = [" ", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan",
                    "1.#IND", "1.#QNAN", "<NA>", "N/A",
                    # "NA",
                    "NULL", "NaN", "None",
                    "n/a", "nan", "null "]

custom_na_values_list = ["Allele", "Gene", "Feature", "Feature_type", "Amino_acids", "Codons",
                            "Existing_variation", "FLAGS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID",
                            "canonical_Allele", "canonical_Gene", "canonical_Feature", "canonical_Feature_type",
                            "canonical_Amino_acids", "canonical_Codons", "canonical_Existing_variation",
                            "canonical_FLAGS", "canonical_SYMBOL", "canonical_SYMBOL_SOURCE", "canonical_HGNC_ID",
                            "INVENTED_COLUMN"]

custom_na_values = {col : custom_na_values_list for col in custom_na_values_list}


def read_maf(file, options = dict()):

    custom_na_values = [" ", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan",
                    "1.#IND", "1.#QNAN", "<NA>", "N/A",
                    # "NA",
                    "NULL", "NaN", "None",
                    "n/a", "nan", "null "]

    return pd.read_table(file, header = 0, sep='\t', na_values = custom_na_values, **options)
