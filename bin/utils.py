# from itertools import product
import json
import pandas as pd


def add_filter(old_filt, add_filt, filt_name):
    """
    old filt is the current FILTER field value
    add_filt is a boolean, either True or False
    filt_name is the name that should be added in the FILTER field in case the add_filt value is True
    """
    if add_filt:
        if old_filt == "PASS":
            return filt_name
        old_filt += ";" + filt_name

    return ";".join( sorted(old_filt.split(";")) )

def to_int_if_possible(string):
    try:
        int(string)
        return int(string)
    except ValueError:
        return None


def filter_maf(maf_df, filter_criteria):
    '''
    Filter a MAF dataframe with filtering information coming from a JSON file.
    {
    'FILTER' : 'notcontains n_rich', 'VAF' : 'le 0.35', 'TYPE' : 'SNV'
    }

    '''

    # Define mappings for operators used in criteria
    operators = {
        'eq': lambda x, y: x == y,
        'ne': lambda x, y: x != y,
        'lt': lambda x, y: x < y,
        'le': lambda x, y: x <= y,
        'gt': lambda x, y: x > y,
        'ge': lambda x, y: x >= y,
        'not': lambda x, y: x != y,
        'notcontains': lambda x, y: ~x.str.contains(y), # (~maf_df["FILTER"].str.contains("not_in_panel"))
        'contains': lambda x, y: x.str.contains(y)
    }

    # Apply filters based on criteria from the JSON file
    for col, criterion in filter_criteria.items():
        if col in maf_df.columns:

            if ' ' in criterion:
                operator, value = criterion.split(maxsplit=1)

                if len(operator) == 2 and operator in operators:
                    # 'VAF' : 'le 0.35'
                    pref_len = maf_df.shape[0]
                    maf_df = maf_df[operators[operator](maf_df[col], float(value))]
                    print(f"Applying {col}:{criterion} filter implied going from {pref_len} mutations to {maf_df.shape[0]} mutations.")

                elif operator in operators:
                    # 'FILTER' : 'notcontains n_rich',
                    pref_len = maf_df.shape[0]
                    maf_df = maf_df[operators[operator](maf_df[col], value)]
                    print(f"Applying {col}:{criterion} filter implied going from {pref_len} mutations to {maf_df.shape[0]} mutations.")

                else:
                    print(f"We have no filtering criteria defined for {col}:{criterion} filter.")


            else:
                # 'TYPE' : 'SNV'
                pref_len = maf_df.shape[0]
                maf_df = maf_df[maf_df[col] == criterion]
                print(f"Applying {col}:{criterion} filter implied going from {pref_len} mutations to {maf_df.shape[0]} mutations.")

    return maf_df







contexts_formatted = ['ACA>A', 'ACC>A', 'ACG>A', 'ACT>A', 'CCA>A', 'CCC>A', 'CCG>A', 'CCT>A', 'GCA>A', 'GCC>A', 'GCG>A', 'GCT>A', 'TCA>A', 'TCC>A', 'TCG>A', 'TCT>A',
                        'ACA>G', 'ACC>G', 'ACG>G', 'ACT>G', 'CCA>G', 'CCC>G', 'CCG>G', 'CCT>G', 'GCA>G', 'GCC>G', 'GCG>G', 'GCT>G', 'TCA>G', 'TCC>G', 'TCG>G', 'TCT>G',
                        'ACA>T', 'ACC>T', 'ACG>T', 'ACT>T', 'CCA>T', 'CCC>T', 'CCG>T', 'CCT>T', 'GCA>T', 'GCC>T', 'GCG>T', 'GCT>T', 'TCA>T', 'TCC>T', 'TCG>T', 'TCT>T',
                        'ATA>A', 'ATC>A', 'ATG>A', 'ATT>A', 'CTA>A', 'CTC>A', 'CTG>A', 'CTT>A', 'GTA>A', 'GTC>A', 'GTG>A', 'GTT>A', 'TTA>A', 'TTC>A', 'TTG>A', 'TTT>A',
                        'ATA>C', 'ATC>C', 'ATG>C', 'ATT>C', 'CTA>C', 'CTC>C', 'CTG>C', 'CTT>C', 'GTA>C', 'GTC>C', 'GTG>C', 'GTT>C', 'TTA>C', 'TTC>C', 'TTG>C', 'TTT>C',
                        'ATA>G', 'ATC>G', 'ATG>G', 'ATT>G', 'CTA>G', 'CTC>G', 'CTG>G', 'CTT>G', 'GTA>G', 'GTC>G', 'GTG>G', 'GTT>G', 'TTA>G', 'TTC>G', 'TTG>G', 'TTT>G']
contexts_no_change = ['ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT',
                        'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT', 'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT']

# subs = [''.join(z) for z in itertools.product('CT', 'ACGT') if z[0] != z[1]]
# flanks = [''.join(z) for z in itertools.product('ACGT', repeat=2)]
# contexts_unformatted = sorted([(a, b) for a, b in itertools.product(subs, flanks)], key=lambda x: (x[0], x[1]))
# contexts_no_change = [b[0]+a[0]+b[1] for a, b in contexts_unformatted]
# contexts_formatted = [b[0]+a[0]+b[1]+'>'+a[1] for a, b in contexts_unformatted]
