#!/usr/bin/env python

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
from scipy.stats import zscore
from matplotlib.backends.backend_pdf import PdfPages

import math
import click
import json
from scipy.stats import fisher_exact
from matplotlib.patches import Patch
from scipy.stats import linregress

from utils import load_samples_info


import matplotlib as mpl
from read_utils import custom_na_values
from utils_plot import metrics_colors_dictionary, plots_general_config


mpl.rcParams.update({
    'axes.titlesize'    : plots_general_config["title_fontsize"],       # Title font size
    'axes.labelsize'    : plots_general_config["xylabel_fontsize"],     # X and Y axis labels
    'xtick.labelsize'   : plots_general_config["xyticks_fontsize"],     # X tick labels
    'ytick.labelsize'   : plots_general_config["xyticks_fontsize"],     # Y tick labels
    'legend.fontsize'   : plots_general_config["legend_fontsize"],      # Legend text
    'figure.titlesize'  : plots_general_config["title_fontsize"],       # Figure suptitle (if used)
})

dot_size = plots_general_config["dot_size_scplot"]
# ### Functions ###

# # Function to filter mutdensity dataframes

# def filter_mutdensities(mutdensity_df, region, mode = 'per_gene', samples = None) :
#     # Filter mutdensity_df
#     if mode == 'per_gene' :
#         filt_mutden_df = mutdensity_df[
#                          (mutdensity_df['SAMPLE_ID'] == 'all_samples')
#                        & (mutdensity_df['MUTTYPES'] == 'SNV') 
#                        & (mutdensity_df['REGIONS'] == region)
#                        & (mutdensity_df['GENE'] != 'ALL_GENES')
#         ].reset_index(drop=True)
    
#     elif mode == 'per_sample' and samples is not None :
#         filt_mutden_df = mutdensity_df[
#                          (mutdensity_df['GENE'] == 'ALL_GENES')
#                        & (mutdensity_df['MUTTYPES'] == 'SNV')
#                        & (mutdensity_df['REGIONS'] == region)
#                        & (mutdensity_df['SAMPLE_ID'].isin(samples))
#         ].reset_index(drop=True)
    
#     elif mode == 'per_sample' and samples is None :
#         raise ValueError("For mode 'per_sample', a list of samples must be provided.")
    
#     else : 
#         raise ValueError(
#             f"Invalid mode: {mode} - Expected value: 'per_sample' or 'per_gene'"
#         )
#     return(filt_mutden_df)

# # Function to calculate z-score from log10 

# def z_score_log10(mutden_df, mode = 'per_sample') :
#     # exclude samples with mutdensity value of 0 - to avoid inf in log10 calculation
#     if mode == 'per_sample' :
#         zero_cases = mutden_df.loc[mutden_df['MUTDENSITY_MB'] == 0, 'SAMPLE_ID'].tolist()
#     elif mode == 'per_gene' :
#         zero_cases = mutden_df.loc[mutden_df['MUTDENSITY_MB'] == 0, 'GENE'].tolist()

#     mut_den_nonzero = mutden_df[mutden_df['MUTDENSITY_MB'] > 0].copy().reset_index(drop=True)

#     # Log10 transformation and z-score 
#     mut_den_nonzero['log10_MUTDENSITY_MB'] = np.log10(mut_den_nonzero['MUTDENSITY_MB'])
#     mut_den_nonzero['zscore'] = zscore(mut_den_nonzero['log10_MUTDENSITY_MB'])
    
#     return(zero_cases, mut_den_nonzero)

# # Function to plot results from z-score distribution

# def plt_violin_omega_qc(mutdensity_zscore_df, mode = 'per_gene', zero_cases_flag = None):
    
#     # Setup subplots
#     fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=False)

#     metrics = ['MUTDENSITY_MB', 'zscore']
#     titles = ['MUTDENSITY_MB', 'Z-score of log10(MUTDENSITY_MB)']

#     for ax, metric, title in zip(axes, metrics, titles):
#         # Violin plot
#         sns.violinplot(data=mutdensity_zscore_df, y=metric, inner=None, color='white', ax=ax)

#         # Jittered points
#         jitter_strength = 0.1
#         x_base = 0
#         xs = x_base + np.random.uniform(-jitter_strength, jitter_strength, size=len(mutdensity_zscore_df))
#         ax.scatter(xs, mutdensity_zscore_df[metric], color='skyblue', s=40, zorder=3)

#         # Threshold lines for zscore
#         if metric == 'zscore':
#             ax.axhline(y=2, color='grey', linestyle='--', linewidth=1)
#             ax.axhline(y=-2, color='grey', linestyle='--', linewidth=1)

#         # Label outliers
#         texts = []
#         for i, row in mutdensity_zscore_df.iterrows():
#             if abs(row['zscore']) > 2:
#                 x = xs[i]
#                 y = row[metric]

#                 if mode == 'per_gene':
#                     texts.append(ax.text(x + 0.02, y, row['GENE'],
#                                          color='red', fontsize=8, va='center'))
#                 elif mode == 'per_sample':
#                     texts.append(ax.text(x + 0.02, y, row['SAMPLE_ID'],
#                                          color='red', fontsize=8, va='center'))

#         if texts:
#             adjust_text(texts, ax=ax,
#                         only_move={'points':'y', 'text':'y'},
#                         arrowprops=None)  # arrows optional

#         ax.set_xticks([])

#         if mode == 'per_gene':
#             ax.set_xlabel(mutdensity_zscore_df['SAMPLE_ID'].iloc[0])
#         if mode == 'per_sample':
#             ax.set_xlabel(mutdensity_zscore_df['GENE'].iloc[0]) 

#         ax.set_ylabel(metric)
#         ax.set_title(title)

#     # Add zero-value genes as a side annotation
#     if zero_cases_flag:
#         zero_text = "Samples with MUTDENSITY_MB=0:\n" + ", ".join(zero_cases_flag)
#         fig.text(0.95, 0.5, zero_text, fontsize=9, color='orange', rotation=90,
#                  va='center', ha='left', wrap=True)

#     # Main title
#     plt.suptitle(f'{mode} - {mutdensity_zscore_df['REGIONS'].iloc[0]}', fontsize=14)
#     plt.tight_layout(rect=[0, 0.03, 0.92, 0.95])  # leave space for side text
    
#     return(fig)

# # main function






# # Compare syn muts observed vs estimated by global loc

def loading_data(observed_syn_file, estimated_syn_file):

    omega_dict = {}
    with open(observed_syn_file, "r") as f:
        observed_files = [line.strip() for line in f if line.strip()]
        for file in observed_files:
            data = pd.read_csv(file, sep = "\t")
            data["SAMPLE"] = file.split('.')[1]
            omega_dict[file] = data
    omega_syn_muts_df = pd.concat(omega_dict.values())


    omega_dict = {}
    with open(estimated_syn_file, "r") as f:
        estimated_files = [line.strip() for line in f if line.strip()]
        for file in estimated_files:
            data = pd.read_csv(file, sep = "\t")
            data["SAMPLE"] = file.split('.')[1]
            omega_dict[file] = data

    omega_syn_muts_gloc_df = pd.concat(omega_dict.values())


    syn_muts_df = omega_syn_muts_df.merge(omega_syn_muts_gloc_df,
                                            on = ["GENE", "SAMPLE"], suffixes = ('_loc','_gloc'),
                                            how = 'outer'
                                            )
    
    # TODO: dangerous column renaming here, make sure to check it in case something is updated in omega
    syn_muts_df.columns = ['GENE', 'Observed_synonymous', 'SAMPLE', 'Estimated_synonymous']
    syn_muts_df = syn_muts_df[['GENE', 'SAMPLE', 'Observed_synonymous', 'Estimated_synonymous']]
    syn_muts_df['Observed_synonymous'] = syn_muts_df['Observed_synonymous'].fillna(0)
    
    # remove all subgenic references
    syn_muts_df = syn_muts_df[~(syn_muts_df["GENE"].str.contains("--"))]

    return syn_muts_df

 
def add_diagonal_line(dff, ax):
    # Add x = y line
    min_val = min(dff['Observed_synonymous'].min(), dff['Estimated_synonymous'].min())
    max_val = max(dff['Observed_synonymous'].max(), dff['Estimated_synonymous'].max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--')


def plot_general(df, output_prefix):
    plt.figure(figsize=(4, 4))
    ax = sns.scatterplot(
        data=df,
        x='Observed_synonymous',
        y='Estimated_synonymous',
        hue='GENE',
        palette='tab10',
        legend = False
    )
    add_diagonal_line(df, ax)
    plt.title('Across genes and samples')
    plt.xlabel('Observed_synonymous')
    plt.ylabel('Estimated_synonymous')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f"{output_prefix}.general_scatter.png", bbox_inches='tight', dpi=300)
    plt.close()
    return None


def plot_per_gene_separated(df, output_prefix):
    genes = df['GENE'].unique()
    samples = df['SAMPLE'].unique()
    palette = sns.color_palette('Set2', n_colors=len(samples))
    sample2color = dict(zip(samples, palette))
    pdf_path = f"{output_prefix}.per_gene_separated.pdf"
    with PdfPages(pdf_path) as pdf:
        # Legend page
        fig_legend, ax_legend = plt.subplots(figsize=(6, 2))
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=sample2color[s], label=s, markersize=8) for s in samples]
        ax_legend.legend(handles=handles, title='SAMPLE', loc='center', ncol=4)
        ax_legend.axis('off')
        fig_legend.tight_layout()
        pdf.savefig(fig_legend)
        plt.close(fig_legend)
        # Data pages
        for gene in genes:
            gene_df = df[df['GENE'] == gene]
            fig, ax = plt.subplots(figsize=(3, 3))
            sns.scatterplot(
                data=gene_df,
                x='Observed_synonymous',
                y='Estimated_synonymous',
                hue='SAMPLE',
                palette=sample2color,
                s=dot_size,
                ax=ax,
                legend=False
            )
            add_diagonal_line(gene_df, ax)
            ax.set_title(f'{gene}',)
            ax.set_xlabel('Observed_synonymous')
            ax.set_ylabel('Estimated_synonymous')
            ax.grid(True, linestyle='--', alpha=0.7)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    return pdf_path


def get_gene_color_mapping(df):
    genes = sorted(df['GENE'].unique())
    palette = sns.color_palette('husl', n_colors=len(genes))
    return dict(zip(genes, palette))


def plot_per_sample_separated(df, output_prefix, gene2color=None):
    samples = df['SAMPLE'].unique()
    genes = sorted(df['GENE'].unique())
    if gene2color is None:
        gene2color = get_gene_color_mapping(df)
    pdf_path = f"{output_prefix}.per_sample_separated.pdf"
    with PdfPages(pdf_path) as pdf:
        # Legend page
        fig_legend, ax_legend = plt.subplots(figsize=(6, 2))
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=gene2color[g], label=g, markersize=8) for g in genes]
        ax_legend.legend(handles=handles, title='GENE', loc='center', ncol=3)
        ax_legend.axis('off')
        fig_legend.tight_layout()
        pdf.savefig(fig_legend)
        plt.close(fig_legend)
        # Data pages
        for sample in samples:
            sample_df = df[df['SAMPLE'] == sample]
            fig, ax = plt.subplots(figsize=(3, 3))
            sns.scatterplot(
                data=sample_df,
                x='Observed_synonymous',
                y='Estimated_synonymous',
                hue='GENE',
                palette=gene2color,
                s=dot_size,
                ax=ax,
                legend=False
            )
            add_diagonal_line(sample_df, ax)
            ax.set_title(f'{sample}')
            ax.set_xlabel('Observed_synonymous')
            ax.set_ylabel('Estimated_synonymous')
            ax.grid(True, linestyle='--', alpha=0.7)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    return pdf_path



## Define functions
def plot_samples_grid(df, output_prefix, gene2color=None):
    
    genes = sorted(df['GENE'].unique())
    if gene2color is None:
        gene2color = get_gene_color_mapping(df)

    samples = df['SAMPLE'].unique()
    n_samples = len(samples)
    ncols = min(4, n_samples)
    nrows = math.ceil(n_samples / ncols)
    fig, axs = plt.subplots(nrows, ncols, figsize=(2.5 * ncols, 2.5 * nrows), sharey=False)
    axs = axs.flatten() if type(axs) is np.ndarray else [axs]
    sample_results = []
    for i, sample in enumerate(samples):
        sample_df = df[df['SAMPLE'] == sample]
        subset_df = sample_df.fillna(0)
        sns.scatterplot(data=subset_df,
                        x='Observed_synonymous',
                        y='Estimated_synonymous',
                        hue='GENE',
                        palette=gene2color,
                        s=dot_size,
                        ax=axs[i],
                        legend=False)
        min_val = -0.05
        max_val = max(subset_df["Observed_synonymous"].max()*1.1, subset_df["Estimated_synonymous"].max()*1.1)
        axs[i].plot([min_val, max_val], [min_val, max_val], color='lightgrey', linestyle='--')
        axs[i].set_xlim(min_val, max_val)
        axs[i].set_ylim(min_val, max_val)
        xtickss = axs[i].get_xticks()
        ytickss = axs[i].get_yticks()
        if len(ytickss) > len(xtickss):
            axs[i].set_xticks(xtickss)
            axs[i].set_yticks(xtickss)
        else:
            axs[i].set_xticks(ytickss)
            axs[i].set_yticks(ytickss)
        axs[i].set_xlim(min_val, max_val)
        axs[i].set_ylim(min_val, max_val)
        x1 = subset_df[["Observed_synonymous", "Estimated_synonymous"]].dropna()["Observed_synonymous"].values
        y1 = subset_df[["Observed_synonymous", "Estimated_synonymous"]].dropna()["Estimated_synonymous"].values
        try:
            res = linregress(x=x1, y=y1)
            r_squared = res.rvalue**2
            p_value = res.pvalue
            if p_value < 0.05:
                p_value = " < 0.05"
            else:
                p_value = " ≥ 0.05"
            axs[i].annotate(f'R-squared: {r_squared:.2f} (p{p_value})', xy=(0.25, 0.98), xycoords='axes fraction',
                            fontsize=plots_general_config["annots_fontsize"], color='black')
            sample_results.append({'SAMPLE': sample,
                                   'R_squared': r_squared,
                                   'p_value': res.pvalue,
                                   'number_mutations': subset_df['Observed_synonymous'].sum().item(),
                                  })
        except ValueError:
            r_squared = "NA"
            p_value = "NA"
            axs[i].annotate(f'R-squared: {r_squared} (p-value: {p_value})', xy=(0.25, 0.98), xycoords='axes fraction',
                            fontsize=plots_general_config["annots_fontsize"], color='black')
            sample_results.append({'SAMPLE': sample,
                                   'R_squared': None,
                                   'p_value': 1,
                                   'number_mutations': subset_df['Observed_synonymous'].sum().item()})
        add_diagonal_line(subset_df, axs[i])
        axs[i].grid(True, linestyle='--', alpha=0.7)
        axs[i].set_ylabel(f"Estimated_synonymous")
        axs[i].set_xlabel(f"Observed_synonymous")
        axs[i].set_title(f"{sample}")
        axs[i].spines['right'].set_visible(False)
        axs[i].spines['top'].set_visible(False)
    # Hide unused subplots
    for j in range(i+1, len(axs)):
        fig.delaxes(axs[j])
    sample_results_df = pd.DataFrame(sample_results)
    plt.tight_layout()
    fig.savefig(f"{output_prefix}.samples_grid.png", bbox_inches='tight', dpi=300)
    plt.close(fig)
    # Save results to TSV
    sample_results_df.to_csv(f"{output_prefix}.samples_grid_results.tsv", sep='\t', index=False)
    # Save non-significant correlations
    nonsig = sample_results_df[(sample_results_df['p_value'] >= 0.05) | (sample_results_df['p_value'].isnull())]
    nonsig.to_csv(f"{output_prefix}.samples_grid_nonsignificant.tsv", sep='\t', index=False)
    return sample_results_df


def plot_genes_grid(df, output_prefix):
    genes = df['GENE'].unique()
    n_genes = len(genes)
    ncols = min(4, n_genes)
    nrows = math.ceil(n_genes / ncols)
    fig, axs = plt.subplots(nrows, ncols, figsize=(2.5 * ncols, 2.5 * nrows), sharey=False)
    axs = axs.flatten() if type(axs) is np.ndarray else [axs]
    gene_results = []
    for i, gene in enumerate(genes):
        gene_df = df[df['GENE'] == gene]
        subset_df = gene_df.fillna(0)
        sns.scatterplot(data=subset_df,
                        x='Observed_synonymous',
                        y='Estimated_synonymous',
                        hue='SAMPLE',
                        palette='Set2',
                        s=dot_size,
                        ax=axs[i], legend=False)
        min_val = -0.05
        max_val = max(subset_df["Observed_synonymous"].max()*1.1, subset_df["Estimated_synonymous"].max()*1.1)
        axs[i].set_xlim(min_val, max_val)
        axs[i].set_ylim(min_val, max_val)
        xtickss = axs[i].get_xticks()
        ytickss = axs[i].get_yticks()
        if len(ytickss) > len(xtickss):
            axs[i].set_xticks(xtickss)
            axs[i].set_yticks(xtickss)
        else:
            axs[i].set_xticks(ytickss)
            axs[i].set_yticks(ytickss)
        axs[i].set_xlim(min_val, max_val)
        axs[i].set_ylim(min_val, max_val)
        x1 = subset_df[["Observed_synonymous", "Estimated_synonymous"]].dropna()["Observed_synonymous"].values
        y1 = subset_df[["Observed_synonymous", "Estimated_synonymous"]].dropna()["Estimated_synonymous"].values
        try:
            res = linregress(x=x1, y=y1)
            r_squared = res.rvalue**2
            p_value = res.pvalue
            if p_value < 0.05:
                p_value = " < 0.05"
            else:
                p_value = " ≥ 0.05"
            axs[i].annotate(f'R-squared: {r_squared:.2f} (p{p_value})', xy=(0.25, 1), xycoords='axes fraction',
                            fontsize=plots_general_config["annots_fontsize"], color='black')
            gene_results.append({'GENE': gene,
                                'R_squared': r_squared,
                                'p_value': res.pvalue,
                                'number_mutations': subset_df['Observed_synonymous'].sum().item(),
                                })
        except ValueError:
            r_squared = "NA"
            p_value = "NA"
            axs[i].annotate(f'R-squared: {r_squared} (p-value: {p_value})', xy=(0.25, 1), xycoords='axes fraction',
                            fontsize=plots_general_config["annots_fontsize"], color='black')
            gene_results.append({'GENE': gene,
                                'R_squared': None,
                                'p_value': 1,
                                'number_mutations': subset_df['Observed_synonymous'].sum().item()})
        add_diagonal_line(subset_df, axs[i])
        axs[i].grid(True, linestyle='--', alpha=0.7)
        axs[i].set_ylabel(f"Estimated_synonymous")
        axs[i].set_xlabel(f"Observed_synonymous")
        axs[i].set_title(f"{gene}")
        axs[i].spines['right'].set_visible(False)
        axs[i].spines['top'].set_visible(False)
    # Hide unused subplots
    for j in range(i+1, len(axs)):
        fig.delaxes(axs[j])
    gene_results_df = pd.DataFrame(gene_results)
    plt.tight_layout()
    fig.savefig(f"{output_prefix}.genes_grid.png", bbox_inches='tight', dpi=300)
    plt.close(fig)
    # Save results to TSV
    gene_results_df.to_csv(f"{output_prefix}.genes_grid_results.tsv", sep='\t', index=False)
    # Save non-significant correlations
    nonsig = gene_results_df[(gene_results_df['p_value'] >= 0.05) | (gene_results_df['p_value'].isnull())]
    nonsig.to_csv(f"{output_prefix}.genes_grid_nonsignificant.tsv", sep='\t', index=False)
    return gene_results_df






def plot_progressive_samples(sample_results_df, output_prefix):
    sorted_results_df = sample_results_df.sort_values(by='R_squared', ascending=True)
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(sorted_results_df['SAMPLE'], sorted_results_df['R_squared'], marker='o', color='b', label='R-squared')
    ax1.set_xlabel("Samples (sorted by R-squared)")
    ax1.set_ylabel("R-squared", color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.set_xticks(range(len(sorted_results_df['SAMPLE'])))
    ax1.set_xticklabels(sorted_results_df['SAMPLE'], rotation=90)
    ax2 = ax1.twinx()
    sns.scatterplot(
        x=range(len(sorted_results_df['number_mutations'])),
        y=sorted_results_df['number_mutations'],
        ax=ax2,
        color='red',
        label='Number of mutations',
        s=dot_size
    )
    ax2.set_ylabel("Number of mutations", color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax1.grid(axis='x', linestyle='--', alpha=0.7)
    plt.title("R-squared and P-values Across Samples")
    fig.savefig(f"{output_prefix}.progressive_samples.png", bbox_inches='tight', dpi=300)
    plt.close(fig)
    return None


def plot_progressive_genes(gene_results_df, output_prefix):
    sorted_results_df = gene_results_df.sort_values(by='R_squared', ascending=True)
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(sorted_results_df['GENE'], sorted_results_df['R_squared'], marker='o', color='b', label='R-squared')
    ax1.set_xlabel("Genes (sorted by R-squared)")
    ax1.set_ylabel("R-squared", color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.set_xticks(range(len(sorted_results_df['GENE'])))
    ax1.set_xticklabels(sorted_results_df['GENE'], rotation=90)
    ax2 = ax1.twinx()
    sns.scatterplot(
        x=range(len(sorted_results_df['number_mutations'])),
        y=sorted_results_df['number_mutations'],
        ax=ax2,
        color='red',
        label='Number of mutations',
        s=dot_size
    )
    ax2.set_ylabel("Number of mutations", color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax1.legend(loc='upper left')
    ax2.legend(loc='lower right')
    ax1.grid(axis='x', linestyle='--', alpha=0.7)
    plt.title("R-squared and P-values Across Genes")
    fig.savefig(f"{output_prefix}.progressive_genes.png", bbox_inches='tight', dpi=300)
    plt.close(fig)
    return None







@click.command()
@click.option("--observed-syn", required=True, type=click.Path(exists=True),
              help="File containing an omega preprocessing synonymous mutation file per row (observed)")
@click.option("--estimated-syn", required=True, type=click.Path(exists=True),
              help="File containing an omega preprocessing synonymous mutation file per row (expected)")
@click.option("--output-prefix", required=True, type=str,
              help="Prefix for the outputs")

def main(observed_syn, estimated_syn, output_prefix):
    samples_names, group_names, _all_unit_names, _all_units_groups = load_samples_info()
    syn_muts_df = loading_data(observed_syn, estimated_syn)

    # TODO: decide if some filters of genes or sample groups need to be applied here
    syn_muts_df = syn_muts_df[syn_muts_df["GENE"] != 'ALL_GENES']

    syn_muts_df_all_samples = syn_muts_df[syn_muts_df["SAMPLE"] == 'all_samples']
    plot_per_sample_separated(syn_muts_df_all_samples, output_prefix + ".all_samples")
    syn_muts_df = syn_muts_df[syn_muts_df["SAMPLE"] != 'all_samples']

    # General scatter
    gene2color = get_gene_color_mapping(syn_muts_df)
    

    # Sample level
    syn_muts_df_samples = syn_muts_df[(syn_muts_df["SAMPLE"].isin(samples_names))]

    plot_general(syn_muts_df_samples, output_prefix + ".samples")
    # Per-gene and per-sample separated (PDF)
    plot_per_gene_separated(syn_muts_df_samples, output_prefix + ".samples")
    plot_per_sample_separated(syn_muts_df_samples, output_prefix + ".samples", gene2color=gene2color)
    # Grid plots (PNG) and progressive (PNG)
    sample_results_df = plot_samples_grid(syn_muts_df_samples, output_prefix + ".samples", gene2color=gene2color)
    gene_results_df = plot_genes_grid(syn_muts_df_samples, output_prefix + ".samples")
    plot_progressive_samples(sample_results_df, output_prefix + ".samples")
    plot_progressive_genes(gene_results_df, output_prefix + ".samples")


    # Group level
    syn_muts_df_groups = syn_muts_df[(syn_muts_df["SAMPLE"].isin(group_names))]
    if syn_muts_df_groups.shape[0] == 0:
        print("No group-level data found. Skipping group-level analyses.")
        return
    plot_general(syn_muts_df_groups, output_prefix + ".groups")
    plot_per_gene_separated(syn_muts_df_groups, output_prefix + ".groups")
    plot_per_sample_separated(syn_muts_df_groups, output_prefix + ".groups", gene2color=gene2color)
    sample_results_df_groups = plot_samples_grid(syn_muts_df_groups, output_prefix + ".groups", gene2color=gene2color)
    gene_results_df_groups = plot_genes_grid(syn_muts_df_groups, output_prefix + ".groups")
    plot_progressive_samples(sample_results_df_groups, output_prefix + ".groups")
    plot_progressive_genes(gene_results_df_groups, output_prefix + ".groups")


if __name__ == '__main__':
    main()
