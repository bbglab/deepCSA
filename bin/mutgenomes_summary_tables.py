#!/usr/bin/env python


import glob
import click

import pandas as pd


from utils import inclusion_exclusion


@click.command()
@click.option('--metadata-file')
def gather(metadata_file):
    remove_sex_chromosome = False

    # parse all samples

    df_all = []
    for fn in glob.glob('./*.covered_genomes_summary.tsv'):
        df = pd.read_csv(fn, sep='\t')
        df_all.append(df)

    df_all = pd.concat(df_all, axis=0)

    # load clinical information
    # Try reading the metadata file with different separators to find the correct one
    for sep in ['\t', ',']:
        try:
            clinical_df_tmp = pd.read_csv(metadata_file, sep=sep)
            if 'SAMPLE_ID' in clinical_df_tmp.columns and 'SEX' in clinical_df_tmp.columns:
                clinical_df = clinical_df_tmp[['SAMPLE_ID', 'SEX']]
                print(f"""Metadata file read with separator: {sep}""")
                break
        except Exception:
            continue
    else:
        print("Metadata input file does not contain 'SAMPLE_ID' and 'SEX' columns")
        print("Without proper metadata no values can be computed in the X chromosome")
        clinical_df = df_all.copy()
        clinical_df["SEX"] = "UNKNOWN"
        clinical_df = clinical_df[['SAMPLE', 'SEX']].copy()
        clinical_df.rename(columns={'SAMPLE': 'SAMPLE_ID'}, inplace=True)
        remove_sex_chromosome = True


    df_all = df_all[~(df_all['SAMPLE'] == 'all_samples')]
    df_all.reset_index(drop=True, inplace=True)
    df_all.fillna(0, inplace=True)

    ## create table per sample-gene-csqn_type

    df_all_annotated = pd.merge(df_all, clinical_df, left_on='SAMPLE', right_on='SAMPLE_ID', how='left')

    # snvs

    for label in ['GENOMES_SNV_AM', 'GENOMES_SNV_ND']:
        for k in ['LOWER', 'MEAN', 'UPPER', 'TOTAL']:

            df_all_annotated[f'CELLS_DOUBLE_HIT_{'_'.join(label.split('_')[1:])}_{k}'] = df_all_annotated[f'{label}_{k}'].values

            df_all_annotated[f'CELLS_SINGLE_HIT_{'_'.join(label.split('_')[1:])}_{k}'] = df_all_annotated.apply(
                lambda r: r[f'{label}_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else 2 * r[f'{label}_{k}'],
                axis=1)

    # indels

    label = 'GENOMES_INDEL_AM'
    k = 'TOTAL'

    df_all_annotated[f'CELLS_DOUBLE_HIT_{'_'.join(label.split('_')[1:])}_{k}'] = df_all_annotated[f'{label}_{k}'].values

    df_all_annotated[f'CELLS_SINGLE_HIT_{'_'.join(label.split('_')[1:])}_{k}'] = df_all_annotated.apply(
        lambda r: r[f'{label}_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else 2 * r[f'{label}_{k}'],
        axis=1)

    ## collapse csqn types

    df_gene = df_all.groupby(by=['SAMPLE', 'gene', 'chr']).agg({
        k: lambda x: inclusion_exclusion(list(x)) for k in df_all.columns if k.startswith('GENOMES')
    }).reset_index(drop=False)

    df_gene_annotated = pd.merge(df_gene, clinical_df, left_on='SAMPLE', right_on='SAMPLE_ID', how='left')

    # snvs
    for label in ['GENOMES_SNV_AM', 'GENOMES_SNV_ND']:
        for k in ['LOWER', 'MEAN', 'UPPER', 'TOTAL']:

            df_gene_annotated[f'CELLS_DOUBLE_HIT_{'_'.join(label.split('_')[1:])}_{k}'] = df_gene_annotated[f'{label}_{k}'].values

            df_gene_annotated[f'CELLS_SINGLE_HIT_{'_'.join(label.split('_')[1:])}_{k}'] = df_gene_annotated.apply(
                lambda r: r[f'{label}_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else 2 * r[f'{label}_{k}'],
                axis=1)

    # indels
    label = 'GENOMES_INDEL_AM'
    k = 'TOTAL'

    df_gene_annotated[f'CELLS_DOUBLE_HIT_{'_'.join(label.split('_')[1:])}_{k}'] = df_gene_annotated[f'{label}_{k}'].values

    df_gene_annotated[f'CELLS_SINGLE_HIT_{'_'.join(label.split('_')[1:])}_{k}'] = df_gene_annotated.apply(
        lambda r: r[f'{label}_{k}'] if ((r['chr'] == 'chrX') and (r['SEX'] == 'M')) else 2 * r[f'{label}_{k}'],
        axis=1)

    # whenever the information of sex is not available, set chromosomes cannot be trusted and we remove them from the tables
    if remove_sex_chromosome:
        df_all_annotated = df_all_annotated[~(df_all_annotated['chr'] == 'chrX')]
        df_gene_annotated = df_gene_annotated[~(df_gene_annotated['chr'] == 'chrX')]
        # df_all_annotated = df_all_annotated[~(df_all_annotated['chr'].isin(['chrX', 'chrY']))]
        # df_gene_annotated = df_gene_annotated[~(df_gene_annotated['chr'].isin(['chrX', 'chrY']))]

    ## collapse genes
    df_sample_annotated = df_gene_annotated.groupby(by=['SAMPLE']).agg({
        k: lambda x: inclusion_exclusion(list(x)) for k in df_gene_annotated.columns if (k.startswith('GENOMES') or k.startswith('CELLS'))
    }).reset_index(drop=False)

    # sort
    df_all_annotated.sort_values(by=['SAMPLE', 'gene', 'impact'], inplace=True)
    df_gene_annotated.sort_values(by=['SAMPLE', 'gene'], inplace=True)
    df_sample_annotated.sort_values(by=['SAMPLE'], inplace=True)

    # dump
    df_all_annotated.to_csv('./covered_genomes_cells.tsv', sep='\t', index=False)
    df_gene_annotated.to_csv('./covered_genomes_cells.genewise.grouped.tsv', sep='\t', index=False)
    df_sample_annotated.to_csv('./covered_genomes_cells.samplewise.tsv', sep='\t', index=False)


if __name__ == '__main__':

    gather()

