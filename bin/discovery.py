#!/usr/bin/env python
# -*- coding: utf-8 -*-


import warnings
warnings.filterwarnings("ignore")

# import tqdm
import click

import numpy as np
import pandas as pd
import os
import tqdm
import functools
import sys

import pandas as pd

import matplotlib.pyplot as plt

import numpy as np

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors




PROTEIN_AFFECTING_SET = {
    'nonsense',
    'missense',
    'essential_splice',
    'protein_altering_variant',
    'transcript_amplification',
}

def prob_min_uniform_sample_below_cut(N, n, cut):
    """
    Probability that the minimum of a sampling of n elements from [N] is lower or equal than cut
    """
    arr = np.array([np.log(N - cut - k) - np.log(N - k) for k in range(n)])
    return 1 - np.exp(np.sum(arr))

def prob_min_uniform_sample_below_cut_vec(N, n, cut):
    """
    Probability that the minimum of a sampling of n elements from [N] is lower or equal than cut
    Vectorized version for speed up.
    Computed for arrays of (N, n, cut) triples.
    """
    N = np.asarray(N, dtype=int)
    n = np.asarray(n, dtype=int)
    cut = np.asarray(cut, dtype=int)

    out = np.zeros_like(N, dtype=float)
    out[cut >= N] = 1.0
    valid = (N > 0) & (n > 0) & (cut > 0) & (cut < N)

    for j in np.where(valid)[0]:
        k = np.arange(n[j])
        log_terms = np.log(N[j] - cut[j] - k) - np.log(N[j] - k)
        out[j] = 1 - np.exp(log_terms.sum())

    return out

def load_mutations(somatic_mutations_file, protein_affecting_type):
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[
        ~(somatic_mutations['FILTER'].str.contains("not_in_exons"))
        & (somatic_mutations['canonical_Protein_affecting'] == protein_affecting_type)
        & (somatic_mutations['TYPE'] == 'SNV')
    ]

    mutations_lite = mutations[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID', 'ALT_DEPTH', 'ALT_DEPTH_AM', 'SYMBOL', 'canonical_Consequence_broader']]
    mutations_lite = mutations_lite.groupby(by=['CHROM', 'POS', 'REF', 'ALT']
                                            ).agg({'ALT_DEPTH_AM': 'sum', 'ALT_DEPTH': 'sum'}).reset_index()
    return mutations_lite


def get_aachange_format(r):

    if r['Protein_position'] == '-':
        return '-'
    elif len(r['Amino_acids']) == 1:
        return r['Amino_acids'] + r['Protein_position'] + r['Amino_acids']
    else:
        return r['Amino_acids'].split('/')[0] + r['Protein_position'] + r['Amino_acids'].split('/')[1]

        
def collect_vep(vep_filename):
    vep_panel = pd.read_csv(vep_filename, sep='\t')
    dg = vep_panel[['CHROM', 'POS', 'REF', 'ALT', 'Protein_position', 'Amino_acids', 'GENE']].copy()
    dg['POS'] = dg['POS'].astype(int)
    dg['AACHANGE'] = dg.apply(get_aachange_format, axis=1)
    return dg


def load_panel(consensus_panel_file, depths_file, samples_group, vep_annotations):

    df_panel = pd.read_csv(consensus_panel_file, sep='\t')
    df_depth = pd.read_csv(depths_file, sep='\t')

    df_panel = pd.merge(df_panel, df_depth[['CHROM', 'POS', samples_group]], on=['CHROM', 'POS'], how='left')
    df_panel.rename(columns={samples_group: 'DEPTH'}, inplace=True)
    df_panel = pd.merge(df_panel, vep_annotations[['CHROM', 'POS', 'REF', 'ALT', 'AACHANGE', 'GENE']], 
                              left_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'], 
                              right_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'],
                              how='left')
    df_panel = df_panel[df_panel['AACHANGE'] != '-']
    df_panel.dropna(inplace=True)
    df_panel['RESIDUE'] = df_panel['AACHANGE'].apply(lambda s: s[:-1])
    return df_panel




def empirical_discovery_index_curve(gene, mutations_dict, df_panel_dict,
                                    subsampling_rates,
                                    replicates=100, sites='genomic'):

    df = mutations_dict[sites]
    df = df[df['GENE'] == gene]

    dg = df_panel_dict[sites]
    dg = dg[dg['GENE'] == gene]

    size = dg.shape[0]
    mean_depth = df['DEPTH'].mean()
    
    x, mean, err_low, err_high = [], [], [], []
    
    unique_dict = {}
    for i, p in tqdm.tqdm(enumerate(subsampling_rates)):
        dist_bernoulli = tfd.Bernoulli(probs=df[f'UNIQUE_RATE_{i}'].values)
        unique_mutations = np.sum(dist_bernoulli.sample(sample_shape=(100,)), axis=1)
        y = list(unique_mutations / size)
        mean += [np.mean(y)]
        err_low += [np.percentile(y, 2.5)]
        err_high += [np.percentile(y, 97.5)]
        x += [mean_depth * p]
    mean += [df.shape[0] / size]
    err_low += [df.shape[0] / size]
    err_high += [df.shape[0] / size]
    x += [mean_depth]

    return x, mean, err_low, err_high


def plot_empirical_discovery(gene, mutations_dict, df_panel_dict, subsampling_rates, sites='genomic'):
    x, mean, err_low, err_high = empirical_discovery_index_curve(gene, mutations_dict, df_panel_dict, subsampling_rates, sites = sites)
    plt.scatter(x, mean, s=50)
    for i, m in enumerate(x):
        plt.vlines(m, err_low[i], err_high[i])
    plt.xscale('log')
    plt.title(gene)
    plt.show()



# theoretical neutral vs empirical discovery curves
# Create a PDF to save the plots
def main_empirical(sample,
                   mutations_dict, df_panel_dict,
                   omega_mutability_file, relative_mutability_file,
                   subsampling_rates,
                   output_folder,
                   sites='genomic', impact = "protein_affecting", logscale=False, genes_list = None):

    mutations_lite = mutations_dict[sites]
    mutations_lite['VAF'] = mutations_lite.apply(lambda r: r['ALT_DEPTH']/r['DEPTH'], axis=1)

    # relative_mutability_file = os.path.join(deepCSA_folder, 'processing_files', 'relativemutability', 'all_samples.all.tsv.gz')
    # omega_mutability_file = os.path.join(deepCSA_folder, 'selection','omega', 'preprocessing', f'mutability_per_sample_gene_context.{sample}.tsv')


    # retrieve relative mutability
    mutability_raw = pd.read_csv(relative_mutability_file, sep='\t')

    mutability_raw = pd.merge(mutability_raw, df_panel, on=['CHROM', 'POS', 'CONTEXT_MUT', 'GENE','IMPACT'], how='left')
    mutability_raw = mutability_raw.rename(columns={sample: "MUTABILITY"})
    # collect VEP annotations

    df_panel = pd.merge(df_panel, vep[['CHROM', 'POS', 'REF', 'ALT', 'AACHANGE', 'GENE']], 
                              left_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'], 
                              right_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'],
                              how='left')

    if genes_list is None:
        genes_list = df_panel['GENE'].unique()

    for gene in tqdm.tqdm(genes_list):

        try:
            synonymous_mutation_rate = pd.read_csv(omega_mutability_file, sep='\t')
            synonymous_mutation_rate = synonymous_mutation_rate[synonymous_mutation_rate['GENE'] == gene]
            mutability_gene = mutability_raw[mutability_raw['GENE'] == gene]
            mutability_gene = pd.merge(mutability_gene, synonymous_mutation_rate[['CONTEXT_MUT', sample]], on=['CONTEXT_MUT'], how='left')
            mutability_gene.rename(columns={sample: 'MUTRATE'}, inplace=True)

            mutability_gene = pd.merge(mutability_gene, df_panel[['CHROM', 'POS', 'REF', 'ALT', 'AACHANGE']], on=['CHROM', 'POS', 'REF', 'ALT'], how='left')

            # discard positions in non-CDS regions, probably splicing and intronic
            mutability_gene = mutability_gene[(mutability_gene['AACHANGE'] != '-') & (~mutability_gene['AACHANGE'].isnull())]

            # keep only protein affecting mutation sites
            if impact == "protein_affecting":
                mutability_gene = mutability_gene[mutability_gene["IMPACT"].isin(PROTEIN_AFFECTING_SET)]
            else:
                mutability_gene = mutability_gene[mutability_gene["IMPACT"] == 'synonymous']

            mutability_gene['RESIDUE'] = mutability_gene['AACHANGE'].apply(lambda s: s[:-1])
            print(mutability_gene.head())
            if sites == 'residue':
                mutability_gene = mutability_gene.groupby(['GENE', 'RESIDUE']).agg({'MUTRATE': 'sum', 'DEPTH': 'mean'}).reset_index()
            elif sites == 'genomic':
                mutability_gene = mutability_gene.groupby(['GENE', 'POS']).agg({'MUTRATE': 'sum', 'DEPTH': 'mean'}).reset_index()
            
            mutability_gene['MUTABILITY'] = mutability_gene.apply(lambda s: (s['MUTRATE'] / s['DEPTH']), axis=1)

            if sites == 'residue':
                mutability_gene = pd.merge(mutability_gene, mutations_lite[['GENE', 'RESIDUE', 'VAF']], on=['GENE', 'RESIDUE'], how='left')
            elif sites == 'genomic':
                mutability_gene = pd.merge(mutability_gene, mutations_lite[['GENE', 'POS', 'VAF']], on=['GENE', 'POS'], how='left')


            # neutral rate
            mutability_gene['RATE_NEUTRAL'] = mutability_gene['MUTABILITY']
            total_neutral_rate = mutability_gene['RATE_NEUTRAL'].sum()

            # compute saturation theoretical
            y_unique_neutral = []
            if sites == 'genomic':
                x_theoretical = np.logspace(3, 8, num=100)
            elif sites == 'residue':
                x_theoretical = np.logspace(3, 7, num=100)

            for depth in x_theoretical:
                unique_neutral = np.sum(1 - np.exp(-mutability_gene['RATE_NEUTRAL'].values * depth)) / mutability_gene.shape[0]
                y_unique_neutral.append(unique_neutral)

            # compute empirical discovery index curve

            x_empirical, mean, err_low, err_high = empirical_discovery_index_curve(gene, mutations_dict, df_panel_dict, subsampling_rates, sites=sites)

            # plot
            fig, ax1 = plt.subplots(figsize=(2,2))
            ax1.set_xscale('log')
            if logscale:
                ax1.set_yscale('log')

            # empirical

            ax1.scatter(x_empirical[:-1], mean[:-1], label='downsampling', color='brown', s=5)
            ax1.scatter(x_empirical[-1], mean[-1], label='observed', color='white', edgecolors='brown', alpha=1, s=100)
            for i, m in enumerate(x_empirical[:-2]):
                plt.vlines(m, err_low[i], err_high[i], color='brown', lw=1)

            # theoretical

            ax1.plot(x_theoretical, y_unique_neutral, color='grey', lw=2, label='neutral theoretical', alpha=0.5)  # neutral
            # ax1.tick_params(axis='y', labelcolor=color)
            if sites == 'residue':
                ax1.set_ylabel('proportion of\nmutated residues')
            elif sites == 'genomic':
                ax1.set_ylabel('proportion of\nmutated nucleotides')
            
            ax1.set_xlabel('depth per residue')
            # ax1.vlines(5e5, 0, 1., linestyles='dashed', color='maroon', label='cohort', alpha=0.3)

            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)

            ax1.set_xlim(x_theoretical[0], x_theoretical[-1])

            # ax1.legend(loc=(1,0))

            plt.title(gene + " (" + impact + ")")

            if logscale:
                plt.savefig(f'{output_folder}/proportion_mutated_sites_{sites}_logscale_{gene}.{impact}.png', bbox_inches='tight', dpi=300)
                print(f"Plot saved to {output_folder}")
            else:
                plt.savefig(f'{output_folder}/proportion_mutated_sites_{sites}_{gene}.{impact}.png', bbox_inches='tight', dpi=300)
                print(f"Plot saved to {output_folder}")

            plt.show()

        except:

            print(gene)
            continue
main_empirical('CohortCha_TimepointT0', sites='genomic', impact = impact, logscale=False, genes_list = ["DNMT3A","TET2","PPM1D","TP53","CHEK2","ASXL1"])
main_empirical('CohortCha_TimepointT0', sites='residue', impact=impact, logscale=False,genes_list = ["DNMT3A","TET2","PPM1D","TP53","CHEK2","ASXL1"])














@click.command()
@click.option("--somatic-mutations-file", required=True, type=click.Path(exists=True),help="Path to the somatic mutations file")
@click.option("--vep-file", required=True, type=click.Path(exists=True),help="Path to the VEP annotation file")
@click.option("--consensus-panel-file", required=True, type=click.Path(exists=True),help="Path to the exons consensus panel file")
@click.option("--depths-file", required=True, type=click.Path(exists=True),help="Path to the depths file")
@click.option("--residue", is_flag=True, show_default=True, default=False, help="either genomic or residue based sites")
@click.option("--group-name", type=str, default="all_samples", show_default=True, help="Name of the group/sample to be used in the code")
@click.option("--impact", type=str, default="protein_affecting", show_default=True, help="either protein_affecting or non_protein_affecting positions")
def cli(somatic_mutations_file, vep_file, consensus_panel_file, depths_file, residue, group_name, impact):

    print(f"Analyzing {group_name}")
    mutations = load_mutations(somatic_mutations_file, impact)
    print("Mutations loaded")
    vep = collect_vep(vep_file)
    df_panel = load_panel(consensus_panel_file, depths_file, group_name, vep)

    print("Panel loaded")
    mutations_lite = pd.merge(mutations, 
                          df_panel[['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'GENE', 'AACHANGE']],
                          on=['CHROM', 'POS', 'REF', 'ALT'], 
                          how='left')
    mutations_lite.dropna(axis=0, inplace=True)  # drop nan sites
    mutations_lite = mutations_lite[mutations_lite['AACHANGE'] != '-']
    mutations_lite['RESIDUE'] = mutations_lite['AACHANGE'].apply(lambda s: s[:-1])

    if residue:
        mutations_lite = mutations_lite.groupby(['GENE', 'RESIDUE']).agg({'ALT_DEPTH': 'sum', 'DEPTH': 'mean'}).reset_index()    
    else:
        mutations_lite = mutations_lite.groupby(['GENE', 'POS']).agg({'ALT_DEPTH': 'sum', 'DEPTH': 'mean'}).reset_index()
    
    subsampling_rates = np.logspace(-2, np.log10(0.9), num=20)
    print("Starting subsampling")

    depth = mutations_lite["DEPTH"].to_numpy()
    alt_depth = mutations_lite["ALT_DEPTH"].to_numpy()
    print(mutations_lite.shape)
    print(depth.shape, alt_depth.shape)

    # for i, p in tqdm.tqdm(enumerate(subsampling_rates), total=len(subsampling_rates)):
    for i, p in enumerate(subsampling_rates):
        n = np.floor(p * depth).astype(int)
        mutations_lite[f"UNIQUE_RATE_{i}"] = prob_min_uniform_sample_below_cut_vec(depth, n, alt_depth)
    
    if residue:
        out_file = f"{group_name}_mutations_residue_rates.{impact}.tsv"
    else:
        out_file = f"{group_name}_mutations_genomic_rates.{impact}.tsv"

    mutations_lite.to_csv(out_file, sep="\t", index=False)        


    output_folder = f'{group_name}.curves.{"residue" if residue else "genomic"}_{impact}'
    os.makedirs(output_folder, exist_ok=True)

    # # load mutations
    mutations_genomic = None #pd.read_csv(f'{base_folder}/analysis/{deepCSA_run}_deepCSA_CH_I_wSP/{deepCSA_run}_CH_I_wSP_custom_analysis/saturation_mutagenesis/saturation_kinetics/{sample_group}_mutations_genomic_rates.{impact}.tsv', sep='\t')
    mutations_residue = mutations_lite


    # df_panel represents the total number of mutable sites,
    # either genomic or residue sites
    if impact == "protein_affecting":
        df_panel = df_panel[df_panel["IMPACT"].isin(PROTEIN_AFFECTING_SET)]
    else:
        df_panel = df_panel[df_panel["IMPACT"] == "synonymous"]    


    df_panel_genomic = df_panel.groupby(['POS', 'GENE']).agg({'DEPTH': 'mean'}).reset_index()
    df_panel_residue = df_panel.groupby(['RESIDUE',  'GENE']).agg({'DEPTH': 'mean'}).reset_index()
    df_panel_dict = {
        'genomic': df_panel_genomic,
        'residue': df_panel_residue
        }

    mutations_dict = {
        'genomic': mutations_genomic,
        'residue': mutations_residue
    }
    subsampling_rates = np.logspace(-2, np.log10(0.9), num=20)

    # plot_empirical_discovery('TP53', mutations_dict, df_panel_dict, subsampling_rates









if __name__ == '__main__':
    cli()