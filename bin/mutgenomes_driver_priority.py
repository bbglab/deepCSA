#!/opt/conda/bin/python


import os
import click

import numpy as np
import pandas as pd
from scipy import stats

from utils import inclusion_exclusion



def compute_prior(alt_read_count, depth):

    mean_p = (1 + alt_read_count.sum()) / depth.sum()
    return {'alpha': 1.0, 'beta': 1 / mean_p}


def compute_excess(omega_file):

    omega_df = pd.read_csv(omega_file, sep='\t')
    omega_df = omega_df[omega_df['impact'].isin(['missense', 'nonsense'])]
    omega_df['excess_dnds'] = omega_df['dnds'].apply(lambda w: (w - 1) / w if w >= 1 else 0)
    omega_df['excess_lower'] = omega_df['lower'].apply(lambda w: (w - 1) / w if w >= 1 else 0)
    omega_df['excess_upper'] = omega_df['upper'].apply(lambda w: (w - 1) / w if w >= 1 else 0)

    res = {}
    for gene in omega_df['gene']:
        res[gene] = {}
        for csqn in ['missense', 'nonsense']:
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


def compute_covered_genomes(alt_read_count, depth, bootstrap_N=100):

    # compute prior
    prior = compute_prior(alt_read_count, depth)
    alpha = prior['alpha']
    beta = prior['beta']

    # compute Clopper-Pearson bounds
    lower_bound = get_binomial_low(alt_read_count, depth)
    upper_bound = get_binomial_high(alt_read_count, depth)

    # direct VAF
    vaf = alt_read_count / depth

    # moderated VAF
    mod_vaf = (alt_read_count + alpha) / (depth + beta + alpha)

    # proportion of mutated genomes based on average depth
    proportion_average_depth = np.sum(alt_read_count) / np.mean(depth)

    res = {
        'vaf': inclusion_exclusion(vaf),  # sparse
        'mod_vaf': inclusion_exclusion(mod_vaf),  # non-sparse
        'average_depth': proportion_average_depth  # global
        }

    return res


@click.group()
def cli():
    pass


@cli.command()
@click.option('--sample')
@click.option('--outfolder')
@click.option('--somatic-mutations-file')
@click.option('--omega-file')
def snv_am(sample, outfolder, somatic_mutations_file, omega_file):

    # parse mutations
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[(somatic_mutations['TYPE'] == 'SNV')
                                    & (somatic_mutations['canonical_Consequence_broader'].isin(['missense', 'nonsense']))
                                    ].reset_index(drop = True)

    mutations['LOWER_VAF_AM'] = mutations.apply(lambda r: stats.beta(r['ALT_DEPTH_AM'], r['DEPTH_AM'] - r['ALT_DEPTH_AM'] + 1).ppf(0.05/2), axis=1)

    gene_chr_dict = { x : y  for x, y in somatic_mutations[['canonical_SYMBOL', 'CHROM']].drop_duplicates().values }

    # parse dN/dS excess
    excess_dict = compute_excess(omega_file)
    genes_with_omega = list(excess_dict.keys())

    # select mutations
    res_dict = {}
    res_dict['gene'] = []
    res_dict['impact'] = []
    res_dict['chr'] = []
    for suffix in ['total', 'upper', 'mean', 'lower']:
        res_dict[f'vaf_{suffix}'] = []
        res_dict[f'mod_vaf_{suffix}'] = []
        res_dict[f'average_depth_{suffix}'] = []

    genes_with_mutations_n_omega = sorted(set(genes_with_omega).intersection(set(mutations['canonical_SYMBOL'].unique())))
    for gene in genes_with_mutations_n_omega:

        for csqn in ['missense', 'nonsense']:

            # total

            df_all = mutations[(mutations['canonical_SYMBOL'] == gene) & (mutations['canonical_Consequence_broader'] == csqn)]
            df_all = df_all.sort_values(by=['LOWER_VAF_AM'], ascending=False)
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

            df_dict = {'total': df_all, 'upper': df_upper, 'mean': df_mean, 'lower': df_lower}

            for suffix in ['total', 'upper', 'mean', 'lower']:

                if df_dict[suffix].shape[0] == 0:

                    res_dict[f'vaf_{suffix}'] = res_dict[f'vaf_{suffix}'] + [0.]
                    res_dict[f'mod_vaf_{suffix}'] = res_dict[f'mod_vaf_{suffix}'] + [0.]
                    res_dict[f'average_depth_{suffix}'] = res_dict[f'average_depth_{suffix}'] + [0.]

                else:

                    depth = df_dict[suffix]['DEPTH_AM'].astype(np.float32).values
                    alt_read_count = df_dict[suffix]['ALT_DEPTH_AM'].astype(np.float32).values
                    res = compute_covered_genomes(alt_read_count, depth)

                    res_dict[f'vaf_{suffix}'] = res_dict[f'vaf_{suffix}'] + [res['vaf']]
                    res_dict[f'mod_vaf_{suffix}'] = res_dict[f'mod_vaf_{suffix}'] + [res['mod_vaf']]
                    res_dict[f'average_depth_{suffix}'] = res_dict[f'average_depth_{suffix}'] + [res['average_depth']]

    df = pd.DataFrame(res_dict)
    df["sample"] = sample
    df.to_csv(f'{outfolder}/{sample}.snv.covered_genomes_summary.tsv', sep='\t', index=False)

@cli.command()
@click.option('--sample')
@click.option('--outfolder')
@click.option('--somatic-mutations-file')
def indel_am(sample, outfolder, somatic_mutations_file):

    # parse mutations
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[(somatic_mutations['TYPE'].isin(['INSERTION', 'DELETION']))].reset_index(drop = True)

    # mutations = mutations[mutations['canonical_Consequence_broader'].isin(['missense', 'nonsense'])]
    mutations['LOWER_VAF_AM'] = mutations.apply(lambda r: stats.beta(r['ALT_DEPTH_AM'], r['DEPTH_AM'] - r['ALT_DEPTH_AM'] + 1).ppf(0.05/2), axis=1)

    gene_chr_dict = { x : y  for x, y in somatic_mutations[['canonical_SYMBOL', 'CHROM']].drop_duplicates().values }


    # select mutations

    res_dict = {}
    res_dict['gene'] = []
    res_dict['chr'] = []
    for suffix in ['total']:
        res_dict[f'vaf_{suffix}'] = []
        res_dict[f'mod_vaf_{suffix}'] = []
        res_dict[f'average_depth_{suffix}'] = []

    for gene in mutations['canonical_SYMBOL'].unique():

        res_dict['gene'] = res_dict['gene'] + [gene]
        res_dict['chr'] = res_dict['chr'] + [gene_chr_dict[gene]]

        # total

        df_all = mutations[(mutations['canonical_SYMBOL'] == gene)]
        df_all = df_all.sort_values(by=['LOWER_VAF_AM'], ascending=False)
        df_all.reset_index(drop=True, inplace=True)

        # summarize in a dictionary

        df_dict = {'total': df_all}

        for suffix in ['total']:

            if df_dict[suffix].shape[0] == 0:

                res_dict[f'vaf_{suffix}'] = res_dict[f'vaf_{suffix}'] + [0.]
                res_dict[f'mod_vaf_{suffix}'] = res_dict[f'mod_vaf_{suffix}'] + [0.]
                res_dict[f'average_depth_{suffix}'] = res_dict[f'average_depth_{suffix}'] + [0.]

            else:

                depth = df_dict[suffix]['DEPTH_AM'].astype(np.float32).values
                alt_read_count = df_dict[suffix]['ALT_DEPTH_AM'].astype(np.float32).values
                res = compute_covered_genomes(alt_read_count, depth)

                res_dict[f'vaf_{suffix}'] = res_dict[f'vaf_{suffix}'] + [res['vaf']]
                res_dict[f'mod_vaf_{suffix}'] = res_dict[f'mod_vaf_{suffix}'] + [res['mod_vaf']]
                res_dict[f'average_depth_{suffix}'] = res_dict[f'average_depth_{suffix}'] + [res['average_depth']]

    df = pd.DataFrame(res_dict)
    df["sample"] = sample
    df.to_csv(f'{outfolder}/{sample}.indel.covered_genomes_summary.tsv', sep='\t', index=False)



@click.command()
@click.option('--sample')
@click.option('--outfolder')
@click.option('--somatic-mutations-file')
@click.option('--omega-file')
def cli_nonduplex(sample, outfolder, somatic_mutations_file, omega_file):

    # parse mutations
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[
        (somatic_mutations['TYPE'] == 'SNV')
        & (somatic_mutations['canonical_Consequence_broader'].isin(['missense', 'nonsense']))
        ].reset_index(drop = True)

    gene_chr_dict = { x : y  for x, y in somatic_mutations[['canonical_SYMBOL', 'CHROM']].drop_duplicates().values }

    mutations = mutations[mutations['ALT_DEPTH_ND'] > 0]

    mutations['LOWER_VAF_ND'] = mutations.apply(lambda r: stats.beta(
        r['ALT_DEPTH_ND'],
        r['DEPTH_ND'] - r['ALT_DEPTH_ND'] + 1).ppf(0.05/2), axis=1)

    # parse dN/dS excess
    excess_dict = compute_excess(omega_file)

    # select mutations

    res_dict = {}
    res_dict['gene'] = []
    res_dict['impact'] = []
    res_dict['chr'] = []
    for suffix in ['total', 'upper', 'mean', 'lower']:
        res_dict[f'vaf_{suffix}'] = []
        res_dict[f'mod_vaf_{suffix}'] = []
        res_dict[f'average_depth_{suffix}'] = []

    for gene in mutations['canonical_SYMBOL'].unique():

        for csqn in ['missense', 'nonsense']:

            # total

            df_all = mutations[(mutations['canonical_SYMBOL'] == gene) & (mutations['canonical_Consequence_broader'] == csqn)]
            df_all = df_all.sort_values(by=['LOWER_VAF_ND'], ascending=False)
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

            df_dict = {'total': df_all, 'upper': df_upper, 'mean': df_mean, 'lower': df_lower}

            for suffix in ['total', 'upper', 'mean', 'lower']:

                if df_dict[suffix].shape[0] == 0:

                    res_dict[f'vaf_{suffix}'] = res_dict[f'vaf_{suffix}'] + [0.]
                    res_dict[f'mod_vaf_{suffix}'] = res_dict[f'mod_vaf_{suffix}'] + [0.]
                    res_dict[f'average_depth_{suffix}'] = res_dict[f'average_depth_{suffix}'] + [0.]

                else:

                    depth = df_dict[suffix]['DEPTH_ND'].astype(np.float32).values
                    alt_read_count = df_dict[suffix]['ALT_DEPTH_ND'].astype(np.float32).values
                    res = compute_covered_genomes(alt_read_count, depth)

                    res_dict[f'vaf_{suffix}'] = res_dict[f'vaf_{suffix}'] + [res['vaf']]
                    res_dict[f'mod_vaf_{suffix}'] = res_dict[f'mod_vaf_{suffix}'] + [res['mod_vaf']]
                    res_dict[f'average_depth_{suffix}'] = res_dict[f'average_depth_{suffix}'] + [res['average_depth']]

    df = pd.DataFrame(res_dict)
    df["sample"] = sample
    df.to_csv(f'{outfolder}/{sample}.nd.covered_genomes_summary.nonduplex.tsv', sep='\t', index=False)


if __name__ == '__main__':

    cli()

