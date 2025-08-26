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
    Filter a MAF dataframe with filtering information coming from a list of tuples.
    This can be either a dictionary transformed to list with the .items() method or by directly creating a list of tuples.
    [('VAF', 'le 0.3'), ('VAF_AM', 'le 0.3'), ('vd_VAF', 'le 0.3'),
    ('DEPTH', 'ge 40'), ('FILTER', 'notcontains n_rich'),
    ('FILTER', 'notcontains cohort_n_rich_uni'), ('FILTER', 'notcontains NM20'),
    ('FILTER', 'notcontains no_pileup_support'), ('FILTER', 'notcontains other_sample_SNP'),
    ('FILTER', 'notcontains low_mappability')]
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
        'notcontains': lambda x, y: x.apply(lambda z : y not in z.split(";")), # (~maf_df["FILTER"].str.contains("not_in_panel"))
        'contains': lambda x, y: x.apply(lambda z : y in z.split(";"))
    }

    # Apply filters based on criteria from the JSON file
    for col, criterion in filter_criteria:

        if isinstance(criterion, bool):
            pref_len = maf_df.shape[0]
            maf_df = maf_df[maf_df[col] == criterion]
            print(f"Applying {col}:{criterion} filter implied going from {pref_len} mutations to {maf_df.shape[0]} mutations.")

        elif ' ' in criterion:
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



def vartype(x,
            letters = ['A', 'T', 'C', 'G'],
            len_SV_lim = 50
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




contexts_formatted = ['ACA>A', 'ACC>A', 'ACG>A', 'ACT>A', 'CCA>A', 'CCC>A', 'CCG>A', 'CCT>A', 'GCA>A', 'GCC>A', 'GCG>A', 'GCT>A', 'TCA>A', 'TCC>A', 'TCG>A', 'TCT>A',
                        'ACA>G', 'ACC>G', 'ACG>G', 'ACT>G', 'CCA>G', 'CCC>G', 'CCG>G', 'CCT>G', 'GCA>G', 'GCC>G', 'GCG>G', 'GCT>G', 'TCA>G', 'TCC>G', 'TCG>G', 'TCT>G',
                        'ACA>T', 'ACC>T', 'ACG>T', 'ACT>T', 'CCA>T', 'CCC>T', 'CCG>T', 'CCT>T', 'GCA>T', 'GCC>T', 'GCG>T', 'GCT>T', 'TCA>T', 'TCC>T', 'TCG>T', 'TCT>T',
                        'ATA>A', 'ATC>A', 'ATG>A', 'ATT>A', 'CTA>A', 'CTC>A', 'CTG>A', 'CTT>A', 'GTA>A', 'GTC>A', 'GTG>A', 'GTT>A', 'TTA>A', 'TTC>A', 'TTG>A', 'TTT>A',
                        'ATA>C', 'ATC>C', 'ATG>C', 'ATT>C', 'CTA>C', 'CTC>C', 'CTG>C', 'CTT>C', 'GTA>C', 'GTC>C', 'GTG>C', 'GTT>C', 'TTA>C', 'TTC>C', 'TTG>C', 'TTT>C',
                        'ATA>G', 'ATC>G', 'ATG>G', 'ATT>G', 'CTA>G', 'CTC>G', 'CTG>G', 'CTT>G', 'GTA>G', 'GTC>G', 'GTG>G', 'GTT>G', 'TTA>G', 'TTC>G', 'TTG>G', 'TTT>G']
contexts_no_change = ['ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT',
                        'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT', 'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT']

contexts_formatted_sigprofiler = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T',
                                    'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T',
                                    'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T',
                                    'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T',
                                    'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T',
                                    'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']


# subs = [''.join(z) for z in itertools.product('CT', 'ACGT') if z[0] != z[1]]
# flanks = [''.join(z) for z in itertools.product('ACGT', repeat=2)]
# contexts_unformatted = sorted([(a, b) for a, b in itertools.product(subs, flanks)], key=lambda x: (x[0], x[1]))
# contexts_no_change = [b[0]+a[0]+b[1] for a, b in contexts_unformatted]
# contexts_formatted = [b[0]+a[0]+b[1]+'>'+a[1] for a, b in contexts_unformatted]


def inclusion_exclusion(plist):

    """A more efficient version of the algorithm"""

    n = len(plist)
    if n > 2:
        p1 = inclusion_exclusion(plist[: n // 2])
        p2 = inclusion_exclusion(plist[n // 2: ])
        return inclusion_exclusion([p1, p2])
    if n == 2:
        return plist[0] + plist[1] - plist[0] * plist[1]
    if n == 1:
        return plist[0]
