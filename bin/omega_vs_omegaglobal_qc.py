#!/usr/bin/env python

# Needed basic packages
import click
import os
import pandas as pd  
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import seaborn as sns
import json
from scipy import stats
from scipy.stats import pearsonr


# Functions
def filter_flagged_genes_per_group(df, flagged_cases_df):
    # Optimized: Vectorized filtering using a multi-index instead of .apply row-by-row
    flagged_index = pd.MultiIndex.from_frame(flagged_cases_df[['cohort', 'gene']].rename(columns={'cohort': 'sample'}))
    df_index = pd.MultiIndex.from_frame(df[['sample', 'gene']])
    
    filtered_df = df[~df_index.isin(flagged_index)].copy()
    return filtered_df


def apply_correlation_and_plotting(df, samples_group, output_dir):

    os.makedirs(output_dir, exist_ok=True) # Ensure output dir exists

    for group in sorted(samples_group):
        subset = df[df['sample'] == group].copy()
        subset[['dnds_omega', 'dnds_omegaglobal']] = subset[['dnds_omega', 'dnds_omegaglobal']].fillna(0)
        print(f'Subset shape for {group}: {subset.shape[0]}')
        
        if subset.empty:
            print(f'No data available for sample: {group}')
            continue

        fig, axes = plt.subplots(1, 2, figsize=(8, 3), sharex=False, sharey=False)
        
        impacts = ['missense', 'truncating']
        
        for i, impact in enumerate(impacts):
            ax = axes[i]
            impact_subset = subset[subset['impact'] == impact]
            impact_subset_reg = impact_subset[(impact_subset['dnds_omega'] > 0) & (impact_subset['dnds_omegaglobal'] > 0)]
            
            if impact_subset.empty:
                ax.set_title(f'No {impact} data for {group}')
                continue
                
            # Scatter plot
            sns.scatterplot(data=impact_subset, x='dnds_omega', y='dnds_omegaglobal', hue='gene', ax=ax)
            
            # Regression line
            sns.regplot(data=impact_subset_reg, x='dnds_omega', y='dnds_omegaglobal', scatter=False, color='red', label='Regression Line', ax=ax)
            
            # X=Y line
            all_vals = pd.concat([impact_subset['dnds_omega'], impact_subset['dnds_omegaglobal']])
            min_val, max_val = all_vals.min(), all_vals.max()
            ax.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', label='x=y')
            
            if impact_subset_reg.shape[0] > 2:
                # Compute correlation
                corr, pval = pearsonr(impact_subset_reg['dnds_omega'], impact_subset_reg['dnds_omegaglobal'])
                
                ax.set_title(f'{impact.capitalize()} mutations\nPearson R: {corr:.3f} (p={pval:.3e})')
            else:
                ax.set_title(f'{impact.capitalize()} mutations\nNot enough data for correlation')

            ax.set_xlabel('dN/dS (Omega)')
            ax.set_ylabel('dN/dS (Omegaglobal)')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.legend().remove()

        fig.suptitle(f'Comparison for {group}')
        fig.tight_layout()

        output_path = os.path.join(output_dir, f"{group}_omega_vs_omegaglobal_qc_summary_plot.pdf")
        with PdfPages(output_path) as pdf:
            pdf.savefig(fig)

        plt.close(fig)

    return


@click.command()
@click.option("--input-omega-file", required=True, type=click.Path(exists=True),
              help="Directory containing selection/omega/all_omegas.tsv")

@click.option("--input-omegaglobal-file", required=True, type=click.Path(exists=True),
              help="Directory containing selection/omegaglobal/all_omegaglobal.tsv")

@click.option("--output-dir", required=True, type=click.Path(writable=True),
              help="Directory where output files will be written")

@click.option("--flagged-genes-omega", required=True, type=click.Path(exists=True),
              help="Directory where flagged_omega genes were stored from /qc/omega_flagged/debug.syn_flagged_gene.tsv")


# Main function
def main(input_omega_file, input_omegaglobal_file, output_dir, flagged_genes_omega): 
    
    # Read tables
    omega_df = pd.read_table(input_omega_file, sep='\t')
    omegaglobal_df = pd.read_table(input_omegaglobal_file, sep='\t')
    flagged_genes_df = pd.read_table(flagged_genes_omega, sep='\t')

    # Extract sample groups
    sample_groups = flagged_genes_df['cohort'].unique().tolist()

    # Apply basic filter for omega and omegaglobal table
    omega_filtered_df = omega_df[~(omega_df['gene'].str.contains("--")) & # discards exons if they are included in the analysis
                    (omega_df['impact'].isin(['missense', 'truncating']))].copy()
    print('Filtered omega table:')
    print(omega_filtered_df.shape)

   # Apply basic filter for omega and omegaglobal table
    omegaglobal_filtered_df = omegaglobal_df[~(omegaglobal_df['gene'].str.contains("--")) & # discards exons if they are included in the analysis
                    (omegaglobal_df['impact'].isin(['missense', 'truncating']))].copy()
    print('Filtered omegaglobal table:')
    print(omegaglobal_filtered_df.shape)

    
    # Merge omega and omegaglobal tables
    omega_vs_omegaglobal_df = omega_filtered_df.merge(omegaglobal_filtered_df, on=['sample', 'gene', 'impact'], how='outer', suffixes=('_omega', '_omegaglobal'))
    print('Merged omega and omegaglobal table:')
    print(omega_vs_omegaglobal_df.shape)

    # Filter flagged omega genes within groups
    filtered_omega_vs_omegaglobal_df = filter_flagged_genes_per_group(omega_vs_omegaglobal_df, flagged_genes_df)

    # Export filtered table
    print('Filtered table from flagged genes, shape:')
    print(filtered_omega_vs_omegaglobal_df.shape)

    filtered_output_path = os.path.join(output_dir, "filtered_omega_omegaglobal_table.tsv")
    filtered_omega_vs_omegaglobal_df.to_csv(filtered_output_path, sep='\t', index=False)

    # Apply function to plot omega values vs omegaglobal per group and compute correlation
    apply_correlation_and_plotting(filtered_omega_vs_omegaglobal_df, sample_groups, output_dir)

if __name__ == '__main__':
    main()