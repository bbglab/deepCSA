import os
import sys
import tqdm
import math
from collections import namedtuple
from multiprocessing import Pool
import glob
import click

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
from scipy import stats

import sys
sys.path.append('/workspace/projects/bladder_ts/notebooks/manuscript_figures/')
from consensus_variables import *







import os
import sys
import tqdm
import math
from collections import namedtuple
from multiprocessing import Pool
import glob
import click

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
from scipy import stats



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



def gather():

    # parse all samples

    df_total = []
    for fn in glob.glob('./covered_genomes_snv_AM/*.covered_genomes_summary.tsv'):
        sample = os.path.basename(fn).split('.')[0]
        df = pd.read_csv(fn, sep='\t')
        df['sample'] = sample
        df_total.append(df)

    df_total = pd.concat(df_total, axis=0)

    # load clinical information

    clinical_fn = '/workspace/projects/bladder_ts/data/complete_cohort/samples_metadata/' \
                  'complete_cohort_bladder.discarded_histofindings_lowmuts.clinical_variables_extended.no_transposed.tsv'

    clinical_df = pd.read_csv(clinical_fn, sep='\t')

    # create table per sample-gene-csqn_type

    df = df_total[['chr', 'sample', 'gene', 'impact', 'vaf_lower', 'vaf_mean', 'vaf_upper', 'vaf_total']]
    df = df.rename(columns={'vaf_lower': 'mg_lower', 'vaf_mean': 'mg_mean', 'vaf_upper': 'mg_upper', 'vaf_total': 'mg_total'})
    df = df[~(df['sample'] == 'all_samples')]
    df.reset_index(drop=True, inplace=True)
    df.fillna(0, inplace=True)

    # per gene and csqn_type

    df_clinical = pd.merge(df, clinical_df, left_on='sample', right_on='SAMPLE_ID', how='left')

    for k in ['lower', 'mean', 'upper', 'total']:

        df_clinical[f'me_all_double_hit_{k}'] = df_clinical[f'mg_{k}'].values

        df_clinical[f'me_double_or_single_hit_{k}'] = df_clinical.apply(
            lambda r: r[f'mg_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else (2 - r[f'mg_{k}']) * r[f'mg_{k}'],
            axis=1)

        df_clinical[f'me_all_single_hit_{k}'] = df_clinical.apply(
            lambda r: r[f'mg_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else 2 * r[f'mg_{k}'],
            axis=1)

    # aggregate csqn types

    df_grouped = df.groupby(by=['sample', 'gene', 'chr']).agg({
        k: lambda x: inclusion_exclusion(list(x)) for k in df.columns if k.startswith('mg')
    }).reset_index(drop=False)

    df_grouped_clinical = pd.merge(df_grouped, clinical_df, left_on='sample', right_on='SAMPLE_ID', how='left')

    for k in ['lower', 'mean', 'upper', 'total']:

        df_grouped_clinical[f'me_all_double_hit_{k}'] = df_grouped_clinical[f'mg_{k}'].values

        df_grouped_clinical[f'me_double_or_single_hit_{k}'] = df_grouped_clinical.apply(
            lambda r: r[f'mg_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else (2 - r[f'mg_{k}']) * r[f'mg_{k}'],
            axis=1)

        df_grouped_clinical[f'me_all_single_hit_{k}'] = df_grouped_clinical.apply(
            lambda r: r[f'mg_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else 2 * r[f'mg_{k}'],
            axis=1)

    # aggregate genes

    df_sample_total = df_grouped_clinical.groupby(by=['sample']).agg({
        k: lambda x: inclusion_exclusion(list(x)) for k in df_grouped_clinical.columns if (k.startswith('mg') or k.startswith('me'))
    }).reset_index(drop=False)

    # append clinical data

    df_sample_total_clinical = pd.merge(df_sample_total, clinical_df, left_on='sample', right_on='SAMPLE_ID', how='left')

    # sort

    df_clinical.sort_values(by=['sample', 'gene', 'impact'], inplace=True)
    df_grouped_clinical.sort_values(by=['sample', 'gene'], inplace=True)
    df_sample_total_clinical.sort_values(by=['sample'], inplace=True)

    # dump

    df_clinical.to_csv('./mutated_genomes.snv_AM.csqn_type.tsv', sep='\t', index=False)
    df_grouped_clinical.to_csv('./mutated_genomes.snv_AM.grouped.tsv', sep='\t', index=False)
    df_sample_total_clinical.to_csv('./mutated_genomes.snv_AM.sample_total.tsv', sep='\t', index=False)

    return df_clinical, df_grouped_clinical, df_sample_total_clinical




df_clinical, df_grouped_clinical, df_sample_total_clinical = gather()























def gather():

    # parse all samples

    df_total = []
    for fn in glob.glob('./covered_genomes_indel_AM/*.covered_genomes_summary.tsv'):
        sample = os.path.basename(fn).split('.')[0]
        df = pd.read_csv(fn, sep='\t')
        df['sample'] = sample
        df_total.append(df)

    df_total = pd.concat(df_total, axis=0)

    # load clinical information

    clinical_fn = '/workspace/projects/bladder_ts/data/complete_cohort/samples_metadata/' \
                  'complete_cohort_bladder.discarded_histofindings_lowmuts.clinical_variables_extended.no_transposed.tsv'

    clinical_df = pd.read_csv(clinical_fn, sep='\t')

    # create table per sample-gene-csqn_type

    df = df_total[['chr', 'sample', 'gene', 'vaf_total']]
    df = df.rename(columns={'vaf_total': 'mg_total'})
    df = df[~(df['sample'] == 'all_samples')]
    df.reset_index(drop=True, inplace=True)
    df.fillna(0, inplace=True)

    # per gene and csqn_type

    df_clinical = pd.merge(df, clinical_df, left_on='sample', right_on='SAMPLE_ID', how='left')

    for k in ['total']:

        df_clinical[f'me_all_double_hit_{k}'] = df_clinical[f'mg_{k}'].values

        df_clinical[f'me_double_or_single_hit_{k}'] = df_clinical.apply(
            lambda r: r[f'mg_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else (2 - r[f'mg_{k}']) * r[f'mg_{k}'],
            axis=1)

        df_clinical[f'me_all_single_hit_{k}'] = df_clinical.apply(
            lambda r: r[f'mg_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else 2 * r[f'mg_{k}'],
            axis=1)

    # aggregate genes

    df_sample_total = df_clinical.groupby(by=['sample']).agg({
        k: lambda x: inclusion_exclusion(list(x)) for k in df_clinical.columns if (k.startswith('mg') or k.startswith('me'))
    }).reset_index(drop=False)

    # append clinical data

    df_sample_total_clinical = pd.merge(df_sample_total, clinical_df, left_on='sample', right_on='SAMPLE_ID', how='left')

    # sort

    df_clinical.sort_values(by=['sample', 'gene'], inplace=True)
    df_sample_total_clinical.sort_values(by=['sample'], inplace=True)

    # dump

    df_clinical.to_csv('./mutated_genomes.indel_AM.grouped.tsv', sep='\t', index=False)
    df_sample_total_clinical.to_csv('./mutated_genomes.indel_AM.sample_total.tsv', sep='\t', index=False)

    return df_clinical, df_sample_total_clinical






df_clinical, df_sample_total_clinical = gather()




























df_consensus = pd.read_csv(f'{deepcsa_run_dir}/createpanels/consensuspanels/consensus.all.tsv', sep='\t')

# gene chromosome dictionary
gene_chr_df = df_consensus.groupby(by=['GENE']).agg({'CHROM': 'first'}).reset_index()
gene_chr_dict = dict(zip(gene_chr_df['GENE'].values, gene_chr_df['CHROM'].values))


def compute_prior(alt_read_count, depth):

    mean_p = (1 + alt_read_count.sum()) / depth.sum()
    return {'alpha': 1.0, 'beta': 1 / mean_p}


def compute_excess(sample):

    fn = os.path.join(deepcsa_run_dir, 'omegagloballoc', f'output_mle.{sample}.multi.global_loc.tsv')
    omega_df = pd.read_csv(fn, sep='\t')
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


def inclusion_exclusion(plist):

    """Recursive version that exploits mutual independence of events"""

    n = len(plist)
    if n > 2:
        p1 = inclusion_exclusion(plist[: n // 2])
        p2 = inclusion_exclusion(plist[n // 2: ])
        return inclusion_exclusion([p1, p2])
    if n == 2:
        return plist[0] + plist[1] - plist[0] * plist[1]
    if n == 1:
        return plist[0]


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
def snv_am(sample, outfolder):

    # parse mutations

    somatic_mutations_file = f'{deepcsa_run_dir}/somaticmutations/{sample}.somatic.mutations.tsv'
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[
        ~(somatic_mutations['FILTER'].str.contains("not_in_panel"))
        & (somatic_mutations['canonical_Protein_affecting'] == 'protein_affecting')
        & (somatic_mutations['TYPE'] == 'SNV')]

    mutations = mutations[mutations['canonical_Consequence_broader'].isin(['missense', 'nonsense'])]
    mutations['LOWER_VAF_AM'] = mutations.apply(lambda r: stats.beta(r['ALT_DEPTH_AM'], r['DEPTH_AM'] - r['ALT_DEPTH_AM'] + 1).ppf(0.05/2), axis=1)

    # parse dN/dS excess

    excess_dict = compute_excess(sample)

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
    df.to_csv(f'{outfolder}/{sample}.covered_genomes_summary.tsv', sep='\t', index=False)


@cli.command()
@click.option('--sample')
@click.option('--outfolder')
def indel_am(sample, outfolder):

    # parse mutations

    somatic_mutations_file = f'{deepcsa_run_dir}/somaticmutations/{sample}.somatic.mutations.tsv'
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[
        ~(somatic_mutations['FILTER'].str.contains("not_in_panel"))
        & (somatic_mutations['canonical_Protein_affecting'] == 'protein_affecting')
        & (somatic_mutations['TYPE'].isin(['INSERTION', 'DELETION']))]

    # mutations = mutations[mutations['canonical_Consequence_broader'].isin(['missense', 'nonsense'])]
    mutations['LOWER_VAF_AM'] = mutations.apply(lambda r: stats.beta(r['ALT_DEPTH_AM'], r['DEPTH_AM'] - r['ALT_DEPTH_AM'] + 1).ppf(0.05/2), axis=1)

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
    df.to_csv(f'{outfolder}/{sample}.covered_genomes_summary.tsv', sep='\t', index=False)



@click.command()
@click.option('--sample')
@click.option('--outfolder')
def cli_nonduplex(sample, outfolder):

    # parse mutations

    somatic_mutations_file = f'{deepcsa_run_dir}/somaticmutations/{sample}.somatic.mutations.tsv'
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[
        ~(somatic_mutations['FILTER'].str.contains("not_in_panel"))
        & (somatic_mutations['canonical_Protein_affecting'] == 'protein_affecting')
        & (somatic_mutations['TYPE'] == 'SNV')]

    mutations = mutations[mutations['canonical_Consequence_broader'].isin(['missense', 'nonsense'])]

    mutations['ALT_DEPTH_NONDUPLEX'] = mutations.apply(lambda x: x['ALT_DEPTH_AM'] - x['ALT_DEPTH'], axis=1)
    mutations['DEPTH_NONDUPLEX'] = mutations.apply(lambda x: x['DEPTH_AM'] - x['DEPTH'], axis=1)

    mutations = mutations[mutations['ALT_DEPTH_NONDUPLEX'] > 0]

    mutations['LOWER_VAF_NONDUPLEX'] = mutations.apply(lambda r: stats.beta(
        r['ALT_DEPTH_NONDUPLEX'],
        r['DEPTH_NONDUPLEX'] - r['ALT_DEPTH_NONDUPLEX'] + 1).ppf(0.05/2), axis=1)

    # parse dN/dS excess

    excess_dict = compute_excess(sample)

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
            df_all = df_all.sort_values(by=['LOWER_VAF_NONDUPLEX'], ascending=False)
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

                    depth = df_dict[suffix]['DEPTH_NONDUPLEX'].astype(np.float32).values
                    alt_read_count = df_dict[suffix]['ALT_DEPTH_NONDUPLEX'].astype(np.float32).values
                    res = compute_covered_genomes(alt_read_count, depth)

                    res_dict[f'vaf_{suffix}'] = res_dict[f'vaf_{suffix}'] + [res['vaf']]
                    res_dict[f'mod_vaf_{suffix}'] = res_dict[f'mod_vaf_{suffix}'] + [res['mod_vaf']]
                    res_dict[f'average_depth_{suffix}'] = res_dict[f'average_depth_{suffix}'] + [res['average_depth']]

    df = pd.DataFrame(res_dict)
    df.to_csv(f'{outfolder}/{sample}.covered_genomes_summary.nonduplex.tsv', sep='\t', index=False)


if __name__ == '__main__':

    cli()

