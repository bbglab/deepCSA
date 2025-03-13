#!/opt/conda/bin/python


import os
import click
from functools import reduce

import numpy as np
import pandas as pd
from scipy import stats

from utils import inclusion_exclusion
from utils_impacts import broadimpact_grouping_dict

def compute_prior(alt_read_count, depth):

    mean_p = (1 + alt_read_count.sum()) / depth.sum()
    return {'alpha': 1.0, 'beta': 1 / mean_p}


def compute_excess(omega_file, chosen_impacts_omega=['missense', 'truncating']):
    """
    intent:
    function that computes the excess of dN/dS for each gene and each impact
    """

    omega_df = pd.read_csv(omega_file, sep='\t')
    omega_df = omega_df[omega_df['impact'].isin(chosen_impacts_omega)]
    omega_df['excess_dnds'] = omega_df['dnds'].apply(lambda w: (w - 1) / w if w >= 1 else 0)
    omega_df['excess_lower'] = omega_df['lower'].apply(lambda w: (w - 1) / w if w >= 1 else 0)
    omega_df['excess_upper'] = omega_df['upper'].apply(lambda w: (w - 1) / w if w >= 1 else 0)

    res = {}
    for gene in omega_df['gene']:
        res[gene] = {}
        for csqn in chosen_impacts_omega:
            try:
                res[gene][csqn] = tuple(omega_df[(omega_df['impact'] == csqn) & (omega_df['gene'] == gene)][['excess_lower', 'excess_dnds', 'excess_upper']].values[0])
            except:
                res[gene][csqn] = (0., 0., 0.)

    return res


def get_binomial_low(x, n, alpha=0.05):

    res = stats.beta(x, n - x + 1).ppf(alpha/2)
    res[x == 0.] = 0.
    return res


def get_binomial_high(x, n, alpha=0.05):

    return stats.beta(x + 1, n - x).ppf(1 - alpha/2)


def compute_covered_genomes(alt_read_count, depth):

    # compute Clopper-Pearson bounds
    lower_bound = get_binomial_low(alt_read_count, depth)
    upper_bound = get_binomial_high(alt_read_count, depth)

    # direct VAF
    vaf = alt_read_count / depth

    res = {
        'GENOMES': inclusion_exclusion(vaf),
        'GENOMES_LOW': inclusion_exclusion(lower_bound),
        'GENOMES_HIGH': inclusion_exclusion(upper_bound)
        }

    return res


def snv_am(sample, somatic_mutations_file, omega_file, chosen_impacts=['missense', 'truncating']):

    # parse mutations
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    consequences_of_interest = [ broadimpact_grouping_dict[csqn] for csqn in chosen_impacts]
    consequences_of_interest = set([item for sublist in consequences_of_interest for item in sublist])

    mutations = somatic_mutations[(somatic_mutations['TYPE'] == 'SNV')
                                    & (somatic_mutations['Consequence_broader'].isin(consequences_of_interest))
                                    ].reset_index(drop = True)

    gene_chr_dict = {x: y  for x, y in somatic_mutations[['GENE', 'CHROM']].drop_duplicates().values}

    # parse dN/dS excess
    excess_dict = compute_excess(omega_file, chosen_impacts_omega = chosen_impacts)
    genes_with_omega = list(excess_dict.keys())

    # select mutations
    res_dict = {}
    res_dict['gene'] = []
    res_dict['impact'] = []
    res_dict['chr'] = []
    for suffix in ['TOTAL', 'UPPER', 'MEAN', 'LOWER']:
        res_dict[f'GENOMES_SNV_AM_{suffix}'] = []

    genes_with_mutations_n_omega = sorted(set(genes_with_omega).intersection(set(mutations['GENE'].unique())))

    for gene in genes_with_mutations_n_omega:

        for csqn in chosen_impacts:

            # total

            df_all = mutations[(mutations['GENE'] == gene) & (mutations['Consequence_broader'].isin(broadimpact_grouping_dict[csqn]))]
            df_all = df_all.sort_values(by=['VAF_AM'], ascending=False)
            df_all.reset_index(drop=True, inplace=True)

            # upper

            ind = int(excess_dict[gene][csqn][2] * df_all.shape[0])
            try:
                assert(ind >= 1)
                df_upper = df_all.loc[:ind - 1].copy()
            except:
                df_upper = pd.DataFrame(columns=df_all.columns)

            # mean

            ind = int(excess_dict[gene][csqn][1] * df_all.shape[0])
            try:
                assert(ind >= 1)
                df_mean = df_all.loc[:ind - 1].copy()
            except:
                df_mean = pd.DataFrame(columns=df_all.columns)

            # lower

            ind = int(excess_dict[gene][csqn][0] * df_all.shape[0])
            try:
                assert(ind >= 1)
                df_lower = df_all.loc[:ind - 1].copy()
            except:
                df_lower = pd.DataFrame(columns=df_all.columns)

            res_dict['gene'] = res_dict['gene'] + [gene]
            res_dict['impact'] = res_dict['impact'] + [csqn]
            res_dict['chr'] = res_dict['chr'] + [gene_chr_dict[gene]]

            # summarize in a dictionary

            df_dict = {'TOTAL': df_all, 'UPPER': df_upper, 'MEAN': df_mean, 'LOWER': df_lower}

            for suffix in ['TOTAL', 'UPPER', 'MEAN', 'LOWER']:

                if df_dict[suffix].shape[0] == 0:
                    res_dict[f'GENOMES_SNV_AM_{suffix}'] = res_dict[f'GENOMES_SNV_AM_{suffix}'] + [0.]

                else:
                    depth = df_dict[suffix]['DEPTH_AM'].astype(np.float32).values
                    alt_read_count = df_dict[suffix]['ALT_DEPTH_AM'].astype(np.float32).values
                    res = compute_covered_genomes(alt_read_count, depth)
                    res_dict[f'GENOMES_SNV_AM_{suffix}'] = res_dict[f'GENOMES_SNV_AM_{suffix}'] + [res['GENOMES']]

    df = pd.DataFrame(res_dict)
    all_columnss = list(df.columns)
    df['SAMPLE'] = sample
    return df[['SAMPLE'] + all_columnss]


def indel_am(sample, somatic_mutations_file):

    # parse mutations
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[(somatic_mutations['TYPE'].isin(['INSERTION', 'DELETION']))].reset_index(drop = True)
    gene_chr_dict = {x: y  for x, y in somatic_mutations[['GENE', 'CHROM']].drop_duplicates().values}

    # select mutations

    res_dict = {}
    res_dict['gene'] = []
    res_dict['chr'] = []
    res_dict['impact'] = []
    res_dict['GENOMES_INDEL_AM_TOTAL'] = []

    for gene in mutations['GENE'].unique():

        res_dict['gene'] = res_dict['gene'] + [gene]
        res_dict['chr'] = res_dict['chr'] + [gene_chr_dict[gene]]
        res_dict['impact'] = res_dict['impact'] + ['indel']

        # total

        df_all = mutations[(mutations['GENE'] == gene)]
        df_all = df_all.sort_values(by=['VAF_AM'], ascending=False)
        df_all.reset_index(drop=True, inplace=True)

        # summarize in a dictionary

        df_dict = {'TOTAL': df_all}

        if df_dict['TOTAL'].shape[0] == 0:

            res_dict[f'GENOMES_INDEL_AM_TOTAL'] = res_dict[f'GENOMES_INDEL_AM_TOTAL'] + [0.]

        else:

            depth = df_dict['TOTAL']['DEPTH_AM'].astype(np.float32).values
            alt_read_count = df_dict['TOTAL']['ALT_DEPTH_AM'].astype(np.float32).values
            res = compute_covered_genomes(alt_read_count, depth)
            res_dict[f'GENOMES_INDEL_AM_TOTAL'] = res_dict[f'GENOMES_INDEL_AM_TOTAL'] + [res['GENOMES']]

    df = pd.DataFrame(res_dict)
    all_columnss = list(df.columns)
    df['SAMPLE'] = sample
    return df[['SAMPLE'] + all_columnss]


def snv_nd(sample, somatic_mutations_file, omega_file, chosen_impacts=['missense', 'truncating']):

    # parse mutations
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    consequences_of_interest = [ broadimpact_grouping_dict[csqn] for csqn in chosen_impacts]
    consequences_of_interest = set([item for sublist in consequences_of_interest for item in sublist])

    mutations = somatic_mutations[
        (somatic_mutations['TYPE'] == 'SNV')
        & (somatic_mutations['Consequence_broader'].isin(consequences_of_interest))
        ].reset_index(drop = True)

    gene_chr_dict = { x : y  for x, y in somatic_mutations[['GENE', 'CHROM']].drop_duplicates().values }

    mutations = mutations[mutations['ALT_DEPTH_ND'] > 0]

    # parse dN/dS excess
    excess_dict = compute_excess(omega_file, chosen_impacts_omega = chosen_impacts)
    genes_with_omega = list(excess_dict.keys())

    # select mutations

    res_dict = {}
    res_dict['gene'] = []
    res_dict['impact'] = []
    res_dict['chr'] = []
    for suffix in ['TOTAL', 'UPPER', 'MEAN', 'LOWER']:
        res_dict[f'GENOMES_SNV_ND_{suffix}'] = []

    genes_with_mutations_n_omega = sorted(set(genes_with_omega).intersection(set(mutations['GENE'].unique())))

    for gene in genes_with_mutations_n_omega:

        for csqn in chosen_impacts:

            # total

            df_all = mutations[(mutations['GENE'] == gene) & (mutations['Consequence_broader'].isin(broadimpact_grouping_dict[csqn]))]
            df_all = df_all.sort_values(by=['VAF_ND'], ascending=False)
            df_all.reset_index(drop=True, inplace=True)

            # upper

            ind = int(excess_dict[gene][csqn][2] * df_all.shape[0])
            try:
                assert(ind >= 1)
                df_upper = df_all.loc[:ind - 1].copy()
            except:
                df_upper = pd.DataFrame(columns=df_all.columns)

            # mean

            ind = int(excess_dict[gene][csqn][1] * df_all.shape[0])
            try:
                assert(ind >= 1)
                df_mean = df_all.loc[:ind - 1].copy()
            except:
                df_mean = pd.DataFrame(columns=df_all.columns)

            # lower

            ind = int(excess_dict[gene][csqn][0] * df_all.shape[0])
            try:
                assert(ind >= 1)
                df_lower = df_all.loc[:ind - 1].copy()
            except:
                df_lower = pd.DataFrame(columns=df_all.columns)

            res_dict['gene'] = res_dict['gene'] + [gene]
            res_dict['impact'] = res_dict['impact'] + [csqn]
            res_dict['chr'] = res_dict['chr'] + [gene_chr_dict[gene]]

            # summarize in a dictionary

            df_dict = {'TOTAL': df_all, 'UPPER': df_upper, 'MEAN': df_mean, 'LOWER': df_lower}

            for suffix in ['TOTAL', 'UPPER', 'MEAN', 'LOWER']:

                if df_dict[suffix].shape[0] == 0:

                    res_dict[f'GENOMES_SNV_ND_{suffix}'] = res_dict[f'GENOMES_SNV_ND_{suffix}'] + [0.]

                else:

                    depth = df_dict[suffix]['DEPTH_ND'].astype(np.float32).values
                    alt_read_count = df_dict[suffix]['ALT_DEPTH_ND'].astype(np.float32).values
                    res = compute_covered_genomes(alt_read_count, depth)
                    res_dict[f'GENOMES_SNV_ND_{suffix}'] = res_dict[f'GENOMES_SNV_ND_{suffix}'] + [res['GENOMES']]

    df = pd.DataFrame(res_dict)
    all_columnss = list(df.columns)
    df['SAMPLE'] = sample
    return df[['SAMPLE'] + all_columnss]


@click.command()
@click.option('--sample')
@click.option('--somatic-mutations-file')
@click.option('--omega-file')
def run_all(sample, somatic_mutations_file, omega_file):
    """
    intent:
    script that combines the 3 covered epithelium functions,
    merges the output in a standardized way, so that we just
    get one output in the mutated_genomes_from_vaf/main.nf
    nextflow module

    """

    df_snv_am = snv_am(sample, somatic_mutations_file, omega_file)
    df_indel_am = indel_am(sample, somatic_mutations_file)
    df_snv_nd = snv_nd(sample, somatic_mutations_file, omega_file)

    all_dataframes = [df_snv_am, df_indel_am, df_snv_nd]
    df_merged = reduce(lambda  left, right: pd.merge(left, right, on=['chr', 'gene', 'impact', 'SAMPLE'], how='outer'), all_dataframes)
    df_merged.to_csv(f'./{sample}.covered_genomes_summary.tsv', sep='\t', index=False)


if __name__ == '__main__':

    run_all()

