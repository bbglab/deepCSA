#!/usr/bin/python

import warnings
warnings.filterwarnings("ignore")

import os
import tqdm
import functools
import sys
import argparse

import pandas as pd

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from discovery import collect_vep, load_panel, PROTEIN_AFFECTING_SET

#### Define pathes ####
base_folder = "/data/bbg/nobackup/prominent/chip"
deepCSA_run = '2026-04-14'
deepCSA_folder = f'{base_folder}/deepCSA/runs/{deepCSA_run}_deepCSA/{deepCSA_run}_CH_I_wSP_filtered_muts'
#### Define parameters ####
impact = "protein_affecting"
sample1 = 'CohortCha_TimepointT0'
sample2 = 'CohortCha_TimepointT1'
genes_list = ["PPM1D","TP53","CHEK2", "DNMT3A", "ASXL1", "TET2"]
print("Loading mutations")
mutations_genomic1 = pd.read_csv(f'{base_folder}/analysis/{deepCSA_run}_deepCSA_CH_I_wSP/{deepCSA_run}_CH_I_wSP_custom_analysis/saturation_mutagenesis/saturation_kinetics/{sample1}_mutations_genomic_rates.{impact}.tsv', sep='\t')
mutations_residue1 = pd.read_csv(f'{base_folder}/analysis/{deepCSA_run}_deepCSA_CH_I_wSP/{deepCSA_run}_CH_I_wSP_custom_analysis/saturation_mutagenesis/saturation_kinetics/{sample1}_mutations_residue_rates.{impact}.tsv', sep='\t')
mutations_genomic2 = pd.read_csv(f'{base_folder}/analysis/{deepCSA_run}_deepCSA_CH_I_wSP/{deepCSA_run}_CH_I_wSP_custom_analysis/saturation_mutagenesis/saturation_kinetics/{sample2}_mutations_genomic_rates.{impact}.tsv', sep='\t')
mutations_residue2 = pd.read_csv(f'{base_folder}/analysis/{deepCSA_run}_deepCSA_CH_I_wSP/{deepCSA_run}_CH_I_wSP_custom_analysis/saturation_mutagenesis/saturation_kinetics/{sample2}_mutations_residue_rates.{impact}.tsv', sep='\t')
subsampling_rates = np.logspace(-2, np.log10(0.9), num=20)
output_folder = f'{base_folder}/analysis/{deepCSA_run}_deepCSA_CH_I_wSP/{deepCSA_run}_CH_I_wSP_custom_analysis/saturation_mutagenesis/saturation_kinetics/plots/{sample1}_vs_{sample2}/'
os.makedirs(output_folder, exist_ok=True)
# compute empirical discovery index curves per gene
print("Loading panel")

vep = collect_vep(deepCSA_folder)
df_panel1 = load_panel(deepCSA_folder, sample1, vep)
df_panel2 = load_panel(deepCSA_folder, sample1, vep)

if impact == "protein_affecting":
    df_panel1 = df_panel1[df_panel1["IMPACT"].isin(PROTEIN_AFFECTING_SET)]
    df_panel2 = df_panel2[df_panel2["IMPACT"].isin(PROTEIN_AFFECTING_SET)]
else:
    df_panel1 = df_panel1[df_panel1["IMPACT"] == "synonymous"]
    df_panel2 = df_panel2[df_panel2["IMPACT"] == "synonymous"]

df_panel_genomic1 = df_panel1.groupby(['POS', 'GENE']).agg({'DEPTH': 'mean'}).reset_index()
df_panel_residue1 = df_panel1.groupby(['RESIDUE',  'GENE']).agg({'DEPTH': 'mean'}).reset_index()

df_panel_genomic2 = df_panel2.groupby(['POS', 'GENE']).agg({'DEPTH': 'mean'}).reset_index()
df_panel_residue2 = df_panel2.groupby(['RESIDUE',  'GENE']).agg({'DEPTH': 'mean'}).reset_index()

df_panel_dict1 = {
       'genomic': df_panel_genomic1,
        'residue': df_panel_residue1
        }

mutations_dict1 = {
        'genomic': mutations_genomic1,
        'residue': mutations_residue1
        }

df_panel_dict2 = {
       'genomic': df_panel_genomic2,
        'residue': df_panel_residue2
        }

mutations_dict2 = {
        'genomic': mutations_genomic2,
        'residue': mutations_residue2
        }

def empirical_discovery_index_curve(gene, mutations_dict, df_panel_dict, replicates=100, sites='genomic'):

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

def plot_empirical_discovery(gene):
    x, mean, err_low, err_high = empirical_discovery_index_curve(gene, sites='genomic')
    plt.scatter(x, mean, s=50)
    for i, m in enumerate(x):
        plt.vlines(m, err_low[i], err_high[i])
    plt.xscale('log')
    plt.title(gene)
    plt.show()

def main_empirical_2samples(sample1, sample2, sites = 'genomic', impact = "protein_affecting", logscale=False, genes_list = None):

    mutations_lite1 = mutations_dict1[sites]
    mutations_lite1['VAF'] = mutations_lite1.apply(lambda r: r['ALT_DEPTH']/r['DEPTH'], axis=1)
    mutations_lite2 = mutations_dict2[sites]
    mutations_lite2['VAF'] = mutations_lite2.apply(lambda r: r['ALT_DEPTH']/r['DEPTH'], axis=1)

    df_panel = pd.read_csv(os.path.join(deepCSA_folder, 'regions', 'consensuspanels', 'consensus.exons_splice_sites.tsv'), sep='\t')

    # include depth per site

    df_depth1 = pd.read_csv(os.path.join(deepCSA_folder, 'depths', 'individual', f'{sample1}.depths.annotated.tsv.gz'), sep='\t')
    df_depth2 = pd.read_csv(os.path.join(deepCSA_folder, 'depths', 'individual', f'{sample2}.depths.annotated.tsv.gz'), sep='\t')
    df_panel1 = pd.merge(df_panel, df_depth1[['CHROM', 'POS', sample1]], on=['CHROM', 'POS'], how='left')
    df_panel2 = pd.merge(df_panel, df_depth2[['CHROM', 'POS', sample2]], on=['CHROM', 'POS'], how='left')
    df_panel1.rename(columns={sample1: 'DEPTH'}, inplace=True)
    df_panel2.rename(columns={sample2: 'DEPTH'}, inplace=True)
    
    # retrieve relative mutability

    mutability_raw1 = pd.read_csv(os.path.join(deepCSA_folder, 'processing_files', 'absolutemutabilities', f'mutabilities_per_site.{sample1}.tsv.gz'), 
                                 sep='\t')
    mutability_raw1 = pd.merge(mutability_raw1, df_panel1, on=['CHROM', 'POS', 'CONTEXT_MUT', 'GENE','IMPACT'], how='left')
    mutability_raw1 = mutability_raw1.rename(columns={sample1: "MUTABILITY"})

    mutability_raw2 = pd.read_csv(os.path.join(deepCSA_folder, 'processing_files', 'absolutemutabilities', f'mutabilities_per_site.{sample2}.tsv.gz'), 
                                 sep='\t')
    mutability_raw2 = pd.merge(mutability_raw2, df_panel2, on=['CHROM', 'POS', 'CONTEXT_MUT', 'GENE','IMPACT'], how='left')
    mutability_raw2 = mutability_raw2.rename(columns={sample2: "MUTABILITY"})

    # collect VEP annotations

    df_panel = pd.merge(df_panel, vep[['CHROM', 'POS', 'REF', 'ALT', 'AACHANGE', 'GENE']], 
                              left_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'], 
                              right_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'],
                              how='left')

    if genes_list is None:
        genes_list = df_panel['GENE'].unique()

    for gene in tqdm.tqdm(genes_list):
#        try:
            print(gene)
            synonymous_mutation_rate1 = pd.read_csv(os.path.join(deepCSA_folder, 'selection', 'omega', 'preprocessing', f'mutability_per_sample_gene_context.{sample1}.tsv'), sep='\t')
            synonymous_mutation_rate1 = synonymous_mutation_rate1[synonymous_mutation_rate1['GENE'] == gene]
            synonymous_mutation_rate2 = pd.read_csv(os.path.join(deepCSA_folder, 'selection','omega', 'preprocessing', f'mutability_per_sample_gene_context.{sample2}.tsv'), sep='\t')
            synonymous_mutation_rate2 = synonymous_mutation_rate2[synonymous_mutation_rate2['GENE'] == gene]
            print("Syn rate done")

            mutability_gene1 = mutability_raw1[mutability_raw1['GENE'] == gene]
            mutability_gene1 = pd.merge(mutability_gene1, synonymous_mutation_rate1[['CONTEXT_MUT', sample1]], on=['CONTEXT_MUT'], how='left')
            mutability_gene1.rename(columns={f'{sample1}': 'MUTRATE'}, inplace=True)
            mutability_gene1 = pd.merge(mutability_gene1, df_panel[['CHROM', 'POS', 'REF', 'ALT', 'AACHANGE']], on=['CHROM', 'POS', 'REF', 'ALT'], how='left')         
            mutability_gene1 = mutability_gene1[(mutability_gene1['AACHANGE'] != '-') & (~mutability_gene1['AACHANGE'].isnull())]
            mutability_gene1['RESIDUE'] = mutability_gene1['AACHANGE'].apply(lambda s: s[:-1])

            mutability_gene2 = mutability_raw2[mutability_raw2['GENE'] == gene]
            mutability_gene2 = pd.merge(mutability_gene2, synonymous_mutation_rate2[['CONTEXT_MUT', sample2]], on=['CONTEXT_MUT'], how='left')
            mutability_gene2.rename(columns={f'{sample2}': 'MUTRATE'}, inplace=True)
            mutability_gene2 = pd.merge(mutability_gene2, df_panel[['CHROM', 'POS', 'REF', 'ALT', 'AACHANGE']], on=['CHROM', 'POS', 'REF', 'ALT'], how='left')         
            mutability_gene2 = mutability_gene2[(mutability_gene2['AACHANGE'] != '-') & (~mutability_gene2['AACHANGE'].isnull())]
            mutability_gene2['RESIDUE'] = mutability_gene2['AACHANGE'].apply(lambda s: s[:-1])

            # keep only required impact sites
            if impact == "protein_affecting":
                mutability_gene1 = mutability_gene1[mutability_gene1["IMPACT"].isin(PROTEIN_AFFECTING_SET)]
                mutability_gene2 = mutability_gene2[mutability_gene2["IMPACT"].isin(PROTEIN_AFFECTING_SET)]
            else:
                mutability_gene1 = mutability_gene1[mutability_gene1["IMPACT"] == 'synonymous']
                mutability_gene2 = mutability_gene2[mutability_gene2["IMPACT"] == 'synonymous']

            if sites == 'residue':
                mutability_gene1 = mutability_gene1.groupby(['GENE', 'RESIDUE']).agg({'MUTRATE': 'sum', 'DEPTH': 'mean'}).reset_index()
                mutability_gene2 = mutability_gene2.groupby(['GENE', 'RESIDUE']).agg({'MUTRATE': 'sum', 'DEPTH': 'mean'}).reset_index()
            elif sites == 'genomic':
                mutability_gene1 = mutability_gene1.groupby(['GENE', 'POS']).agg({'MUTRATE': 'sum', 'DEPTH': 'mean'}).reset_index()
                mutability_gene2 = mutability_gene2.groupby(['GENE', 'POS']).agg({'MUTRATE': 'sum', 'DEPTH': 'mean'}).reset_index()
            
            mutability_gene1['MUTABILITY'] = mutability_gene1.apply(lambda s: (s['MUTRATE'] / s['DEPTH']), axis=1)
            mutability_gene2['MUTABILITY'] = mutability_gene2.apply(lambda s: (s['MUTRATE'] / s['DEPTH']), axis=1)

            if sites == 'residue':
                mutability_gene1 = pd.merge(mutability_gene1, mutations_lite1[['GENE', 'RESIDUE', 'VAF']], on=['GENE', 'RESIDUE'], how='left')
                mutability_gene2 = pd.merge(mutability_gene2, mutations_lite2[['GENE', 'RESIDUE', 'VAF']], on=['GENE', 'RESIDUE'], how='left')
            elif sites == 'genomic':
                mutability_gene1 = pd.merge(mutability_gene1, mutations_lite1[['GENE', 'POS', 'VAF']], on=['GENE', 'POS'], how='left')
                mutability_gene2 = pd.merge(mutability_gene2, mutations_lite2[['GENE', 'POS', 'VAF']], on=['GENE', 'POS'], how='left')
            print("Mutability done")
                        
            # neutral rate
            
            mutability_gene1['RATE_NEUTRAL'] = mutability_gene1['MUTABILITY']
            total_neutral_rate1 = mutability_gene1['RATE_NEUTRAL'].sum()
            mutability_gene2['RATE_NEUTRAL'] = mutability_gene2['MUTABILITY']
            total_neutral_rate2 = mutability_gene2['RATE_NEUTRAL'].sum()
            
            # compute saturation theoretical
            print("Rates per gene done")

            y_unique_neutral1 = []
            y_unique_neutral2 = []
            if sites == 'genomic':
                x_theoretical = np.logspace(3, 8, num=100)
            elif sites == 'residue':
    #            if gene == "PPM1D":
    #                x_theoretical = np.logspace(3, 6, num=100)
    #            else:
                x_theoretical = np.logspace(3, 7, num=100)
            
            for depth in x_theoretical:
                unique_neutral1 = np.sum(1 - np.exp(-mutability_gene1['RATE_NEUTRAL'].values * depth)) / mutability_gene1.shape[0]
                y_unique_neutral1.append(unique_neutral1)
                unique_neutral2 = np.sum(1 - np.exp(-mutability_gene2['RATE_NEUTRAL'].values * depth)) / mutability_gene2.shape[0]
                y_unique_neutral2.append(unique_neutral2)

            # compute empirical discovery index curve

            x_empirical1, mean1, err_low1, err_high1 = empirical_discovery_index_curve(gene, mutations_dict1, df_panel_dict1, sites=sites)
            x_empirical2, mean2, err_low2, err_high2 = empirical_discovery_index_curve(gene, mutations_dict2, df_panel_dict2, sites=sites)
            print("Empirical done")

            # plot

            fig, ax1 = plt.subplots(figsize=(2,2))
            ax1.set_xscale('log')
            if logscale:
                ax1.set_yscale('log')

            # empirical

            ax1.scatter(x_empirical1[:-1], mean1[:-1], label='downsampling1', color='black', s=5)
            ax1.scatter(x_empirical2[:-1], mean2[:-1], label='downsampling2', color='brown', s=5)
            ax1.scatter(x_empirical1[-1], mean1[-1], label='observed1', color='white', edgecolors='black', alpha=1, s=100)
            ax1.scatter(x_empirical2[-1], mean2[-1], label='observed2', color='white', edgecolors='brown', alpha=1, s=100)
            for i, m in enumerate(x_empirical1[:-2]):
                plt.vlines(m, err_low1[i], err_high1[i], color='black', lw=1)
            for i, m in enumerate(x_empirical2[:-2]):
                plt.vlines(m, err_low2[i], err_high2[i], color='brown', lw=1)
            print("Empirical plotted")

            # theoretical

            ax1.plot(x_theoretical, y_unique_neutral1, color='grey', lw=2, label='neutral theoretical1', alpha=0.5)  # neutral
            ax1.plot(x_theoretical, y_unique_neutral2, color='sandybrown', lw=2, label='neutral theoretical2', alpha=0.5)  # neutral
            # ax1.tick_params(axis='y', labelcolor=color)
            if sites == 'residue':
                ax1.set_ylabel('proportion mutated residues')
            elif sites == 'genomic':
                ax1.set_ylabel('proportion mutated nucleotides')
            print("Theoretical plotted")
            
            ax1.set_xlabel('depth per residue')
            # ax1.vlines(5e5, 0, 1., linestyles='dashed', color='maroon', label='cohort', alpha=0.3)

            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            
            ax1.legend(loc=(1,0))
            ylims = {
                "PPM1D": (0, 0.4),
            }
            if gene in ylims:
                ax1.set_ylim(ylims[gene])
            plt.title(gene)
            ax1.set_xlim(10**3, max(x_theoretical))
            
            if logscale:
                plt.savefig(f'{output_folder}/proportion_mutated_sites_{sites}_logscale_{gene}.{impact}.png', bbox_inches='tight', dpi=300)
                print(f"Plot saved to {output_folder}")
            else:
                plt.savefig(f'{output_folder}/proportion_mutated_sites_{sites}_{gene}.{impact}.png', bbox_inches='tight', dpi=300)
                print(f"Plot saved to {output_folder}")

            plt.show()
        
#        except:
            
#            print(gene)
#            continue




def main(args):
    # path with mutations


    print("Plotting genomic...")
    main_empirical_2samples(sample1, sample2, sites='genomic', logscale=False,genes_list = genes_list)
    print("Plotting residue ...")
    main_empirical_2samples(sample1, sample2, sites='residue', logscale=False,genes_list = genes_list)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)