#!/opt/conda/bin/python


import glob
import click

import pandas as pd


from utils import inclusion_exclusion


def gather_snvs(metadata_fn):

    # parse all samples

    df_total = []
    for fn in glob.glob('./covered_genomes_snv_AM/*.covered_genomes_summary.tsv'):
        df = pd.read_csv(fn, sep='\t')
        df_total.append(df)

    df_total = pd.concat(df_total, axis=0)

    # load clinical information
    clinical_df = pd.read_csv(metadata_fn, sep='\t')[['SAMPLE_ID', 'SEX']]

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

    df_gene = df.groupby(by=['sample', 'gene', 'chr']).agg({
        k: lambda x: inclusion_exclusion(list(x)) for k in df.columns if k.startswith('mg')
    }).reset_index(drop=False)

    df_gene_clinical = pd.merge(df_gene, clinical_df, left_on='sample', right_on='SAMPLE_ID', how='left')

    for k in ['lower', 'mean', 'upper', 'total']:

        df_gene_clinical[f'me_all_double_hit_{k}'] = df_gene_clinical[f'mg_{k}'].values

        df_gene_clinical[f'me_double_or_single_hit_{k}'] = df_gene_clinical.apply(
            lambda r: r[f'mg_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else (2 - r[f'mg_{k}']) * r[f'mg_{k}'],
            axis=1)

        df_gene_clinical[f'me_all_single_hit_{k}'] = df_gene_clinical.apply(
            lambda r: r[f'mg_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else 2 * r[f'mg_{k}'],
            axis=1)

    # aggregate genes

    df_sample_total = df_gene_clinical.groupby(by=['sample']).agg({
        k: lambda x: inclusion_exclusion(list(x)) for k in df_gene_clinical.columns if (k.startswith('mg') or k.startswith('me'))
    }).reset_index(drop=False)

    # sort

    df_clinical.sort_values(by=['sample', 'gene', 'impact'], inplace=True)
    df_gene_clinical.sort_values(by=['sample', 'gene'], inplace=True)
    df_sample_total.sort_values(by=['sample'], inplace=True)

    # dump

    df_clinical.to_csv('./mutated_genomes.snv_AM.csqn_type.tsv', sep='\t', index=False)
    df_gene_clinical.to_csv('./mutated_genomes.snv_AM.grouped.tsv', sep='\t', index=False)
    df_sample_total.to_csv('./mutated_genomes.snv_AM.sample_total.tsv', sep='\t', index=False)

    return df_clinical, df_gene_clinical, df_sample_total






def gather_indels(metadata_fn):

    # parse all samples
    df_total = []
    for fn in glob.glob('./covered_genomes_indel_AM/*.covered_genomes_summary.tsv'):
        df = pd.read_csv(fn, sep='\t')
        df_total.append(df)

    df_total = pd.concat(df_total, axis=0)

    # load clinical information
    clinical_df = pd.read_csv(metadata_fn, sep='\t')[['SAMPLE_ID', 'SEX']]

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



@click.group()
def cli():
    pass


@cli.command()
@click.option('--metadata-file')
def gather_all(metadata_file):
    df_clinical, df_gene_clinical, df_sample_total_clinical = gather_snvs(metadata_file)
    df_clinical, df_sample_total_clinical = gather_indels(metadata_file)


if __name__ == '__main__':

    cli()

