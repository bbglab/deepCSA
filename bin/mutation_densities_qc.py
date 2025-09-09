#!/usr/bin/env python

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
from scipy.stats import zscore
from matplotlib.backends.backend_pdf import PdfPages
import click
import json


### Functions ###

# Function to load sample names

def load_samples_names(samples_file):

    with open(samples_file, 'r') as f:
        samples_definition = json.load(f)
    all_samples_names = list(samples_definition.keys())

    return all_samples_names

# Function to filter mutdensity dataframes

def filter_mutdensities(mutdensity_df, region, mode = 'per_gene', samples = None) :
    # Filter mutdensity_df
    if mode == 'per_gene' :
        filt_mutden_df = mutdensity_df[
                         (mutdensity_df['SAMPLE_ID'] == 'all_samples')
                       & (mutdensity_df['MUTTYPES'] == 'SNV') 
                       & (mutdensity_df['REGIONS'] == region)
                       & (mutdensity_df['GENE'] != 'ALL_GENES')
        ].reset_index(drop=True)
    
    elif mode == 'per_sample' and samples is not None :
        filt_mutden_df = mutdensity_df[
                         (mutdensity_df['GENE'] == 'ALL_GENES')
                       & (mutdensity_df['MUTTYPES'] == 'SNV')
                       & (mutdensity_df['REGIONS'] == region)
                       & (mutdensity_df['SAMPLE_ID'].isin(samples))
        ].reset_index(drop=True)
    
    elif mode == 'per_sample' and samples is None :
        raise ValueError("For mode 'per_sample', a list of samples must be provided.")
    
    else : 
        raise ValueError(
            f"Invalid mode: {mode} - Expected value: 'per_sample' or 'per_gene'"
        )
    return(filt_mutden_df)

# Function to calculate z-score from log10 

def z_score_log10(mutden_df, mode = 'per_sample') :
    # exclude samples with mutdensity value of 0 - to avoid inf in log10 calculation
    if mode == 'per_sample' :
        zero_cases = mutden_df.loc[mutden_df['MUTDENSITY_MB'] == 0, 'SAMPLE_ID'].tolist()
    elif mode == 'per_gene' :
        zero_cases = mutden_df.loc[mutden_df['MUTDENSITY_MB'] == 0, 'GENE'].tolist()

    mut_den_nonzero = mutden_df[mutden_df['MUTDENSITY_MB'] > 0].copy().reset_index(drop=True)

    # Log10 transformation and z-score 
    mut_den_nonzero['log10_MUTDENSITY_MB'] = np.log10(mut_den_nonzero['MUTDENSITY_MB'])
    mut_den_nonzero['zscore'] = zscore(mut_den_nonzero['log10_MUTDENSITY_MB'])
    
    return(zero_cases, mut_den_nonzero)

# Function to plot results from z-score distribution

def plt_violin_omega_qc(mutdensity_zscore_df, mode = 'per_gene', zero_cases_flag = None):
    
    # Setup subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=False)

    metrics = ['MUTDENSITY_MB', 'zscore']
    titles = ['MUTDENSITY_MB', 'Z-score of log10(MUTDENSITY_MB)']

    for ax, metric, title in zip(axes, metrics, titles):
        # Violin plot
        sns.violinplot(data=mutdensity_zscore_df, y=metric, inner=None, color='white', ax=ax)

        # Jittered points
        jitter_strength = 0.1
        x_base = 0
        xs = x_base + np.random.uniform(-jitter_strength, jitter_strength, size=len(mutdensity_zscore_df))
        ax.scatter(xs, mutdensity_zscore_df[metric], color='skyblue', s=40, zorder=3)

        # Threshold lines for zscore
        if metric == 'zscore':
            ax.axhline(y=2, color='grey', linestyle='--', linewidth=1)
            ax.axhline(y=-2, color='grey', linestyle='--', linewidth=1)

        # Label outliers
        texts = []
        for i, row in mutdensity_zscore_df.iterrows():
            if abs(row['zscore']) > 2:
                x = xs[i]
                y = row[metric]

                if mode == 'per_gene':
                    texts.append(ax.text(x + 0.02, y, row['GENE'],
                                         color='red', fontsize=8, va='center'))
                elif mode == 'per_sample':
                    texts.append(ax.text(x + 0.02, y, row['SAMPLE_ID'],
                                         color='red', fontsize=8, va='center'))

        if texts:
            adjust_text(texts, ax=ax,
                        only_move={'points':'y', 'text':'y'},
                        arrowprops=None)  # arrows optional

        ax.set_xticks([])

        if mode == 'per_gene':
            ax.set_xlabel(mutdensity_zscore_df['SAMPLE_ID'].iloc[0])
        if mode == 'per_sample':
            ax.set_xlabel(mutdensity_zscore_df['GENE'].iloc[0]) 

        ax.set_ylabel(metric)
        ax.set_title(title)

    # Add zero-value genes as a side annotation
    if zero_cases_flag:
        zero_text = "Samples with MUTDENSITY_MB=0:\n" + ", ".join(zero_cases_flag)
        fig.text(0.95, 0.5, zero_text, fontsize=9, color='orange', rotation=90,
                 va='center', ha='left', wrap=True)

    # Main title
    plt.suptitle(f'{mode} - {mutdensity_zscore_df['REGIONS'].iloc[0]}', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.92, 0.95])  # leave space for side text
    
    return(fig)

# main function

@click.command()
@click.option("--input-file", required=True, type=click.Path(exists=True),
              help="Directory containing mutdensity/all_mutdensities.tsv")
@click.option("--output-dir", required=True, type=click.Path(writable=True),
              help="Directory where output files will be written")
@click.option("--panel", required=True, type=click.Path(exists=True),
              help="File with a list of panel genes (one gene per line)")
@click.option("--samples", required=True, type=click.Path(exists=True),
              help="json file with the list of samples to be included in the analysis")

def main(input_file, output_dir, panel, samples) :

    mutden_df = pd.read_table(input_file, sep='\t')
    panel = pd.read_table(panel)

    panel_genes_init = panel['GENE'].unique().tolist()
    del panel

    sample_names = load_samples_names(samples)

    ## 1. Initial general filtering by panel genes
    panel_genes = panel_genes_init + ['ALL_GENES'] # is this an issue for the per_gene omega?
    mutden_df_panel = mutden_df[mutden_df['GENE'].isin(panel_genes)].reset_index(drop=True)

    ## 2. Generating csv with z-score and list of flagged samples/genes
    regions = ['synonymous', 'non_protein_affecting']
    mode_list = ['per_sample', 'per_gene']

    for reg in regions : 
        for mod in mode_list: 
            filt_mutden = filter_mutdensities(mutden_df_panel, reg, mod, sample_names)
            zero_cases_flag, mutden_zscore = z_score_log10(filt_mutden, mod)

            ## Generate csv with zscore result and store as csv
            cols_of_interest = ['SAMPLE_ID', 'GENE', 'REGIONS', 'MUTTYPES', 'MUTDENSITY_MB', 'log10_MUTDENSITY_MB', 'zscore']
            mutden_zscore[cols_of_interest].to_csv(f'{output_dir}/{reg}_{mod}_mutdensities_zscore.csv', index=False)

            if mod == 'per_sample':
                flagged_cases = mutden_zscore[
                                        (mutden_zscore['zscore'] > 2) 
                                        | (mutden_zscore['zscore'] < -2)
                ][['SAMPLE_ID', 'zscore']]
            elif mod == 'per_gene':
                flagged_cases = mutden_zscore[
                                        (mutden_zscore['zscore'] > 2) 
                                        | (mutden_zscore['zscore'] < -2)
                ][['GENE', 'zscore']]

            # Store flagged cases in a csv
            
            flagged_cases['reason_exclusion'] = np.where(flagged_cases['zscore'] > 2, 'high_mutdensity - zscore > 2', 'low_mutdensity - zscore < -2')
            flagged_cases = flagged_cases.drop(columns=['zscore'])
            flagged_cases.columns = ['ID', 'reason_exclusion']

            zero_cases = pd.DataFrame(zero_cases_flag, columns=['ID'])
            zero_cases['reason_exclusion'] = 'mutdensity = 0'  

            flagged_cases_all = pd.concat([flagged_cases, zero_cases], ignore_index=True)
            
            flagged_cases_all.to_csv(
                f"{output_dir}/flagged_cases_{reg}_{mod}.csv", index=False, header = False
            )

            # Violin plots z-score
            with PdfPages(f"{output_dir}/{reg}_{mod}_omegaqc_summary_plot.pdf") as pdf:
                fig = plt_violin_omega_qc(mutden_zscore, mode=mod, zero_cases_flag = zero_cases_flag)

                # If your function doesn't return fig, you can get current figure:
                if fig is None:
                    print(f'Error: unable to plot {reg} {mod}')
                    continue
                            
                # Save to PDF
                pdf.savefig(fig)

                # Close to free memory before next loop
                plt.close(fig)

if __name__ == '__main__':
    main()
