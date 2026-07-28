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
def filter_omega_tables(df):
    # Apply basic filter for omega/omegaglobal tables
    filtered_omega = df[~(df['gene'].str.contains("--")) & # discards exons if they are included in the analysis
                    (df['impact'].isin(['missense', 'truncating']))].copy()
    print('Filtered omega/omegaglobal table:')
    print(filtered_omega.shape)
    return filtered_omega

def process_dndscv_table(dndscv_df):

    # Separate missense and truncating variants for dNdScv table in two separated tables
    missense_df = dndscv_df[['sample', 'gene_name', 'n_mis', 'wmis_cv', 'qmis_cv']].copy()
    missense_df = missense_df.rename(columns={'n_mis': 'mutations', 'wmis_cv': 'dnds', 'qmis_cv': 'pvalue', 'gene_name': 'gene'})
    missense_df['impact'] = 'missense'

    truncating_df = dndscv_df[['sample', 'gene_name', 'n_non', 'n_spl', 'wnon_cv', 'qtrunc_cv']].copy()
    truncating_df['n_trunc'] = truncating_df['n_non'] + truncating_df['n_spl']
    truncating_df = truncating_df.rename(columns={'n_trunc': 'mutations', 'wnon_cv': 'dnds', 'qtrunc_cv': 'pvalue', 'gene_name': 'gene'})
    truncating_df['impact'] = 'truncating'

    # Then concat the two tables
    dndscv_cv_df = pd.concat([missense_df, truncating_df.drop(columns=['n_non', 'n_spl'])], axis=0)
    print('Processed dndscv table: ')
    print(dndscv_cv_df.shape)
    
    return dndscv_cv_df

def process_and_analyze_comparison(
    primary_df, 
    comparison_df, 
    flagged_genes_df, 
    sample_groups, 
    output_dir, 
    comp_label
):
    """
    Unified helper to merge, filter flagged genes, export TSV, and apply plotting/correlation.
    """
    suffix = f'_{comp_label}'
    tsv_filename = f"filtered_omega_{comp_label}_table.tsv"

    # Merge comparison tables (note primary_df is always omega_df)
    merged_df = primary_df.merge(
        comparison_df, 
        on=['sample', 'gene', 'impact'], 
        how='outer', 
        suffixes=('_omega', suffix)
    )
    print(f'Merged omega and {comp_label} table shape:', merged_df.shape)

    # Filter flagged genes
    filtered_df = filter_flagged_genes_per_group(merged_df, flagged_genes_df)
    print(f'Filtered table shape ({comp_label}):', filtered_df.shape)

    # Export filtered table
    filtered_output_path = os.path.join(output_dir, tsv_filename)
    filtered_df.to_csv(filtered_output_path, sep='\t', index=False)

    # Compute correlation and plot (passing comp_label to handle dynamic columns and filenames)
    apply_correlation_and_plotting(
        df=filtered_df, 
        samples_group=sample_groups, 
        output_dir=output_dir, 
        comp_label=comp_label
    )

    return


def filter_flagged_genes_per_group(df, flagged_cases_df):
    # Optimized: Vectorized filtering using a multi-index instead of .apply row-by-row
    flagged_index = pd.MultiIndex.from_frame(flagged_cases_df[['cohort', 'gene']].rename(columns={'cohort': 'sample'}))
    df_index = pd.MultiIndex.from_frame(df[['sample', 'gene']])
    
    filtered_df = df[~df_index.isin(flagged_index)].copy()
    return filtered_df


def apply_correlation_and_plotting(df, samples_group, output_dir, comp_label):
    """
    Plots dNdS comparisons dynamically based on the target dataset (comp_label).
    
    comp_label: e.g. 'omegaglobal' or 'dndscv'
    Column expected: 'dnds_omega' vs f'dnds_{comp_label}'
    """

    os.makedirs(output_dir, exist_ok=True) # Ensure output dir exists

    # Define dynamic column name and display labels
    x_col = 'dnds_omega'
    y_col = f'dnds_{comp_label}'
    target_display_name = "dNdScv" if comp_label == "dndscv" else "Omega Global"

    for group in sorted(samples_group):
        subset = df[df['sample'] == group].copy()

        # Fill missing dNdS values with 0 dynamically
        if x_col in subset.columns and y_col in subset.columns:
            subset[[x_col, y_col]] = subset[[x_col, y_col]].fillna(0)
        print(f'Subset shape for {group} ({target_display_name}): {subset.shape[0]}')
                
        if subset.empty:
            print(f'No data available for sample: {group}')
            continue

        fig, axes = plt.subplots(1, 2, figsize=(8, 3), sharex=False, sharey=False)
        
        impacts = ['missense', 'truncating']
        
        for i, impact in enumerate(impacts):
            ax = axes[i]
            impact_subset = subset[subset['impact'] == impact]
            impact_subset_reg = impact_subset[(impact_subset[x_col] > 0) & (impact_subset[y_col] > 0)]
            
            if impact_subset.empty:
                ax.set_title(f'No {impact} data for {group}')
                continue
                
            # Scatter plot
            sns.scatterplot(data=impact_subset, x=x_col, y=y_col, hue='gene', ax=ax, legend=False)
            
            # Regression line
            sns.regplot(data=impact_subset_reg, x=x_col, y=y_col, scatter=False, color='red', label='Regression Line', ax=ax)
            
            # X=Y line
            all_vals = pd.concat([impact_subset[x_col], impact_subset[y_col]])
            min_val, max_val = all_vals.min(), all_vals.max()
            ax.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', label='x=y')
            
            if impact_subset_reg.shape[0] > 2:
                # Compute correlation
                corr, pval = pearsonr(impact_subset_reg[x_col], impact_subset_reg[y_col])
                
                ax.set_title(f'{impact.capitalize()} mutations\nPearson R: {corr:.3f} (p={pval:.3e})')
            else:
                ax.set_title(f'{impact.capitalize()} mutations\nNot enough data for correlation')

            ax.set_xlabel('dN/dS (Omega)')
            ax.set_ylabel(f'dN/dS ({target_display_name})')

        fig.suptitle(f'Comparison for {group}')
        fig.tight_layout()

        output_path = os.path.join(output_dir, f"{group}_omega_vs_{comp_label}_qc_summary_plot.pdf")
        with PdfPages(output_path) as pdf:
            pdf.savefig(fig)

        plt.close(fig)

    return



@click.command()
@click.option("--input-omega-file", required=True, type=click.Path(exists=True),
              help="Directory containing selection/omega/all_omegas.tsv")

@click.option("--input-omegaglobal-file", required=True, type=click.Path(exists=True),
              help="Directory containing selection/omegaglobal/all_omegaglobal.tsv")

@click.option("--input-dndscv-file", required=False, default=None, type=click.Path(exists=True), #set to default and required false/none to avoid processing when is missing
              help="Directory containing selection/dndscv/cv/all_dNdScv.cv.tsv")

@click.option("--output-dir", required=True, type=click.Path(writable=True),
              help="Directory where output files will be written")

@click.option("--flagged-genes-omega", required=True, type=click.Path(exists=True),
              help="Directory where flagged_omega genes were stored from /qc/omega_flagged/debug.syn_flagged_gene.tsv")

@click.option("--defined_groups", required=True, type=click.Path(exists=True),
              help="User defined sample groups in the analysis")



# Main function
def main(input_omega_file, input_omegaglobal_file, input_dndscv_file, output_dir, flagged_genes_omega, defined_groups): 

    # Read tables; initialize dndscv to avoid errors
    omega_df = pd.read_table(input_omega_file, sep='\t')
    omegaglobal_df = pd.read_table(input_omegaglobal_file, sep='\t')
    flagged_genes_df = pd.read_table(flagged_genes_omega, sep='\t')

    # Read groups from imported json file and load it as list
    with open(defined_groups, 'r') as file:
        json_data = json.load(file)
    sample_groups = list(json_data.keys())
    print(sample_groups)

    # Apply basic filters to omega/omegaglobal tables
    omega_filtered_df = filter_omega_tables(omega_df)
    omegaglobal_filtered_df = filter_omega_tables(omegaglobal_df)

    # Initialize a dictionary for comparisons (label used : table to compare)
    comparisons = {
        'omegaglobal': omegaglobal_filtered_df
    }

    # Add dndscv to comparisons dictionary if provided
    if input_dndscv_file:
        print(f"Processing optional dNdScv file: {input_dndscv_file}")
        dndscv_df = pd.read_table(input_dndscv_file, sep='\t')
        dndscv_cv_df = process_dndscv_table(dndscv_df)
        comparisons['dndscv'] = dndscv_cv_df

    # Execute helper function for each comparison in dictionary:

    for comp_label, comp_df in comparisons.items():
        process_and_analyze_comparison(
            primary_df=omega_filtered_df,
            comparison_df=comp_df,
            flagged_genes_df=flagged_genes_df,
            sample_groups=sample_groups,
            output_dir=output_dir,
            comp_label=comp_label
        )


if __name__ == '__main__':
    main()
