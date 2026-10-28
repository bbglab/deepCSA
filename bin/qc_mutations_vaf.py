#!/usr/bin/env python

"""
Plot VAF and Mutation Density vs Depth with hyperbolic curves.

This script generates scatter plots showing the relationship between sequencing depth
and VAF or mutation density metrics. It overlays hyperbolic curves (N/depth for N=1,2,3...)
and reports counts of data points following these curves.
"""

import click
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from read_utils import custom_na_values
from utils_plot import plots_general_config
import matplotlib as mpl


mpl.rcParams.update({
    'axes.titlesize': plots_general_config["title_fontsize"],
    'axes.labelsize': plots_general_config["xylabel_fontsize"],
    'xtick.labelsize': plots_general_config["xyticks_fontsize"],
    'ytick.labelsize': plots_general_config["xyticks_fontsize"],
    'legend.fontsize': plots_general_config["legend_fontsize"],
    'figure.titlesize': plots_general_config["title_fontsize"],
})


def add_hyperbolic_curves(ax, max_depth, max_n=10, alpha=0.3, linewidth=1):
    """
    Add hyperbolic curves N/depth to a plot.
    
    Parameters:
    -----------
    ax : matplotlib axis
        The axis to plot on
    max_depth : float
        Maximum depth value for plotting
    max_n : int
        Maximum value of N for N/depth curves
    alpha : float
        Transparency of the curves
    linewidth : float
        Width of the curve lines
    
    Returns:
    --------
    dict : Dictionary mapping N values to the curve coordinates
    """
    depth_range = np.linspace(1, max_depth, 1000)
    curves = {}
    
    for n in range(1, max_n + 1):
        vaf_curve = n / depth_range
        # Only plot if VAF is reasonable (< 1.0)
        valid_mask = vaf_curve <= 1.0
        if valid_mask.sum() > 0:
            ax.plot(depth_range[valid_mask], vaf_curve[valid_mask], 
                   '--', alpha=alpha, linewidth=linewidth, 
                   color='gray',
                #    label=f'{n}/depth' if n <= 3 else ''
                   )
            curves[n] = (depth_range[valid_mask], vaf_curve[valid_mask])
    
    return curves


def count_mutations_on_curves(depth, vaf, max_n=10, tolerance=0.05):
    """
    Count how many mutations follow each hyperbolic curve within tolerance.
    
    Parameters:
    -----------
    depth : array-like
        Depth values for mutations
    vaf : array-like
        VAF values for mutations
    max_n : int
        Maximum N value to check
    tolerance : float
        Relative tolerance for considering a point on a curve
    
    Returns:
    --------
    dict : Dictionary with counts for each N value
    """
    counts = {}
    assigned = np.zeros(len(depth), dtype=bool)
    
    for n in range(1, max_n + 1):
        expected_vaf = n / depth
        # Calculate relative difference
        relative_diff = np.abs(vaf - expected_vaf) / (expected_vaf + 1e-10)
        # Count mutations within tolerance that haven't been assigned yet
        on_curve = (relative_diff <= tolerance) & ~assigned & (expected_vaf <= 1.0)
        counts[n] = on_curve.sum()
        assigned |= on_curve
    
    counts['other'] = (~assigned).sum()
    return counts


def plot_vaf_vs_depth_per_site(maf_df, output_pdf, sample_name, max_n=10,
                               depth_col = 'DEPTH',
                               vaf_col='VAF',
                               alt_depth_col='ALT_DEPTH'
                               ):
    """
    Plot VAF vs depth for each mutation site.
    
    Parameters:
    -----------
    maf_df : DataFrame
        MAF dataframe with DEPTH and VAF columns
    output_pdf : PdfPages
        PDF file to save plots
    sample_name : str
        Name of the sample
    max_n : int
        Maximum N for hyperbolic curves
    """
    # Filter valid data
    plot_data = maf_df[[depth_col, vaf_col, alt_depth_col]].dropna()
    plot_data = plot_data[(plot_data[depth_col] > 0) & (plot_data[vaf_col] >= 0) & (plot_data[vaf_col] <= 1.0)]
    
    if plot_data.empty:
        print(f"No valid data for VAF vs depth plot")
        return
    

    # Count mutations on curves
    counts = count_mutations_on_curves(
        plot_data[depth_col].values, 
        plot_data[vaf_col].values, 
        max_n=max_n
    )

    # Add hyperbolic curves
    max_depth = plot_data[depth_col].quantile(0.99)

    for x_max in [max_depth, max_depth/2, max_depth/4, max_depth/8]:
        for y_max in [0.001, 0.005, 0.5, 1.0]:

            # Create figure
            fig, ax = plt.subplots(figsize=(6, 6))
            
            # Scatter plot
            ax.scatter(plot_data[depth_col], plot_data[vaf_col], 
                    alpha=0.5, s=10, edgecolors='none', c='steelblue')

            add_hyperbolic_curves(ax, max_depth, max_n=max_n)
            
            # Add text with counts
            text_str = "Mutations on curves:\n"
            for n in range(1, min(max_n + 1, 6)):
                text_str += f"  {n}/depth: {counts[n]}\n"
            if max_n > 5:
                remaining = sum(counts[n] for n in range(6, max_n + 1))
                text_str += f"  6-{max_n}/depth: {remaining}\n"
            text_str += f"  Other: {counts['other']}"
            
            ax.text(0.98, 0.98, text_str, transform=ax.transAxes,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                fontsize=plots_general_config["annots_fontsize"])
            
            ax.set_xlabel('Depth (reads)', fontsize=plots_general_config["xlabel_fontsize"])
            ax.set_ylabel('VAF', fontsize=plots_general_config["ylabel_fontsize"])
            ax.set_title(f'{sample_name} - VAF vs Depth per Site (N={len(plot_data)})', 
                        fontsize=plots_general_config["title_fontsize"])
            ax.set_xlim(0, x_max)
            ax.set_ylim(0, y_max)
            ax.legend(loc='upper right', fontsize=plots_general_config["legend_fontsize"])
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            output_pdf.savefig()
            plt.close()


def plot_vaf_depth_heatmap(maf_df, output_pdf, sample_name):
    """
    Plot a 2D histogram (heatmap) of VAF vs Depth.
    
    Parameters:
    -----------
    maf_df : DataFrame
        MAF dataframe with DEPTH and VAF columns
    """
    df = maf_df[['DEPTH', 'VAF', 'VAF_AM', 'DEPTH_AM']].dropna()
    df = df[(df['DEPTH'] > 0) & (df['VAF'] >= 0) & (df['VAF'] <= 1.0)]
    
    plt.figure(figsize=(6,6))
    sns.displot(data = df,
                    x = "DEPTH",
                    y = "VAF",
                    #hue = "Protein_affecting",
                    binwidth=(300, .00005),
                    cbar=True,
                    cbar_kws={"label": 'Number of mutations per tile'}
                )
    plt.ylim(-0.0001,0.005)
    plt.xlim(0,15000)
    plt.xlabel('Depth (reads)', fontsize=plots_general_config["xlabel_fontsize"])
    plt.ylabel('VAF', fontsize=plots_general_config["ylabel_fontsize"])
    plt.title(f'{sample_name} - VAF vs Depth per Site (N={len(df)})', 
                fontsize=plots_general_config["title_fontsize"])
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()
    plt.show()


    plt.figure(figsize=(6,6))
    sns.displot(data = df,
                    x = "DEPTH_AM",
                    y = "VAF_AM",
                    #hue = "Protein_affecting",
                    binwidth=(300, .00005),
                    cbar=True,
                    cbar_kws={"label": 'Number of mutations per tile'}
                )
    plt.ylim(-0.0001,0.005)
    plt.xlim(0,15000)
    plt.xlabel('Depth AM (reads)', fontsize=plots_general_config["xlabel_fontsize"])
    plt.ylabel('VAF AM', fontsize=plots_general_config["ylabel_fontsize"])
    plt.title(f'{sample_name} - VAF AM vs Depth AM per mutated site\n(N={len(df)})', 
                fontsize=plots_general_config["title_fontsize"])
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()
    plt.show()


    plt.figure(figsize=(6,6))
    sns.displot(data = df,
                    x = "VAF",
                    y = "VAF_AM",
                    #hue = "Protein_affecting",
                    binwidth=(.0005, .0005),
                    cbar=True,
                    cbar_kws={"label": 'Number of mutations per tile'}
                )
    plt.ylim(-0.005,0.4)
    plt.xlim(-0.005,0.4)
    plt.xlabel('VAF', fontsize=plots_general_config["xlabel_fontsize"])
    plt.ylabel('VAF AM', fontsize=plots_general_config["ylabel_fontsize"])
    plt.title(f'{sample_name} - VAF AM vs VAF per mutated site\n(N={len(df)})', 
                fontsize=plots_general_config["title_fontsize"])
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()
    plt.show()


    plt.figure(figsize=(6,6))
    sns.displot(data = df,
                    x = "VAF",
                    y = "VAF_AM",
                    #hue = "Protein_affecting",
                    binwidth=(.00005, .00005),
                    cbar=True,
                    cbar_kws={"label": 'Number of mutations per tile'}
                )
    plt.ylim(-0.0001,0.005)
    plt.xlim(-0.0001,0.005)
    plt.xlabel('VAF', fontsize=plots_general_config["xlabel_fontsize"])
    plt.ylabel('VAF AM', fontsize=plots_general_config["ylabel_fontsize"])
    plt.title(f'{sample_name} - VAF AM vs VAF per mutated site\n(N={len(df)})', 
                fontsize=plots_general_config["title_fontsize"])
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()
    plt.show()




def get_top_mutations_all_samples(maf_file, output_prefix, top_n=50):
    """
    Get top N mutations by VAF and save to CSV.
    
    Parameters:
    -----------
    maf_file : str
        Path to MAF file
    output_prefix : str
        Prefix for output CSV file
    top_n : int
        Number of top mutations to select
    """
    
    maf_df = pd.read_csv(maf_file, sep='\t', na_values=custom_na_values)
    most_repeated_mutations = maf_df["MUT_ID"].value_counts().reset_index()
    most_repeated_mutations.columns = ["MUT_ID", "COUNT"]
    repeated_maf_df = maf_df.merge(most_repeated_mutations,
                                   on="MUT_ID")[['COUNT',
                                                 'MUT_ID',
                                                 'CONTEXT_MUT', 'canonical_SYMBOL',
                                                 'canonical_Consequence_single',
                                                 'canonical_Protein_position',
                                                 'canonical_Amino_acids'
                                                 ]].drop_duplicates().sort_values(by='COUNT', ascending=False)

    output_csv_path = f"{output_prefix}.cohort_most_repeated_mutations.tsv"
    repeated_maf_df.iloc[:200,:].to_csv(output_csv_path, sep='\t', index=False)
    print(f"Most repeated mutations saved to {output_csv_path}")

    most_repeated_mutations = maf_df.groupby(by =["MUT_ID",
                                                  "canonical_SYMBOL",
                                                  "canonical_Consequence_single",
                                                  'canonical_Protein_position',
                                                  "canonical_Amino_acids"
                                                  ]).agg({"VAF": "mean",
                                                            "ALT_DEPTH": "mean",
                                                            "SAMPLE_ID": "size",
                                                            "FILTER": lambda x: ';'.join(sorted(x.unique()))
                                                            }).reset_index()
    most_repeated_mutations = most_repeated_mutations[most_repeated_mutations["SAMPLE_ID"] > 1]

    if not most_repeated_mutations.empty:
        for criteria in ['VAF', 'ALT_DEPTH']:
            top_mutations = most_repeated_mutations.nlargest(top_n, criteria)
            
            output_csv_path = f"{output_prefix}.cohort_top_{top_n}_mutations_by_{criteria}.tsv"
            top_mutations.to_csv(output_csv_path, sep='\t', index=False)
            print(f"Top {top_n} mutations by {criteria} saved to {output_csv_path}")


def get_top_mutations(maf_file, output_prefix, top_n=50):
    """
    Get top N mutations by VAF and save to CSV.
    
    Parameters:
    -----------
    maf_file : str
        Path to MAF file
    output_prefix : str
        Prefix for output CSV file
    top_n : int
        Number of top mutations to select
    """
    
    maf_df = pd.read_csv(maf_file, sep='\t', na_values=custom_na_values)
    for criteria in ['VAF', 'ALT_DEPTH']:
        maf_df = maf_df.dropna(subset=[criteria])
        top_mutations = maf_df.nlargest(top_n, criteria)
        top_mutations_small = top_mutations[['MUT_ID', #'CHROM', 'POS', 'REF', 'ALT',
                                             'FILTER',
                                             'CONTEXT_MUT', 'canonical_SYMBOL',
                                             'canonical_Consequence_single',
                                             'canonical_Protein_position',
                                             'canonical_Amino_acids',
                                             'VAF', 'DEPTH', 'ALT_DEPTH',
                                             'VAF_AM', 'DEPTH_AM', 'ALT_DEPTH_AM',
                                             ]]
        
        output_csv_path = f"{output_prefix}.top_{top_n}_mutations_by_{criteria}.tsv"
        top_mutations_small.to_csv(output_csv_path, sep='\t', index=False)
        print(f"Top {top_n} mutations by {criteria} saved to {output_csv_path}")
        
def plot_vaf_vs_vafam_histogram (maf_df, output_pdf):
    """
    Plot VAF distribution compared to VAF_AM in a histogram.
    
    Parameters:
    -----------
    maf_df : DataFrame
        MAF dataframe containing VAF and VAF_AM columns
    """
    df_long = maf_df.melt(value_vars=["VAF", "VAF_AM"], 
                       var_name="VAF_type", 
                       value_name="VAF_value")

    fig, ax = plt.subplots(figsize=(5, 4))
    sns.histplot(data=df_long, x="VAF_value", hue="VAF_type", 
                 ax=ax, palette=["dodgerblue", "orange"],
                 log_scale=True, element="step")

    plt.title('VAF vs VAF_AM distribution', 
                fontsize=plots_general_config["title_fontsize"])
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()
    plt.show()


def vaf_pseudocount(alt_depth, depth, weight, prior_vaf=None):    
    return (alt_depth + prior_vaf * weight) / (depth + weight)

def plot_vaf_pseudocount_curve(maf_df, output_pdf, suffix=''):
    """
    Plot VAF distribution compared to VAF_AM in a histogram.
    
    Parameters:
    -----------
    maf_df : DataFrame
        MAF dataframe containing VAF and VAF_AM columns
    """

    fig, axes = plt.subplots(3, 3, figsize=(8, 8))
    axes = axes.flatten()
    average_depth = maf_df[f'DEPTH{suffix}'].mean()
    dg = maf_df
    prior_vaf = dg[f'ALT_DEPTH{suffix}'].sum() / dg[f'DEPTH{suffix}'].sum()
    fig.suptitle(f'VAF{suffix} with pseudocounts, prior vaf={prior_vaf:.2e}, avg. depth={average_depth:.1f}')
    for i, weight_prop in enumerate([0, 0.2, 0.5, 0.75, 1, 1.25, 1.5, 2, 3]):
        weight = weight_prop * average_depth // 1
        dg['VAF_PSEUDO'] = vaf_pseudocount(dg[f'ALT_DEPTH{suffix}'], dg[f'DEPTH{suffix}'], weight, prior_vaf=prior_vaf)
        axes[i].scatter(dg[f'DEPTH{suffix}'], dg['VAF_PSEUDO'], s=3, alpha=0.1)
        rho = np.corrcoef(np.log(dg[f'DEPTH{suffix}']), np.log(dg['VAF_PSEUDO']))[0, 1]
        axes[i].set_xlabel(f'DEPTH{suffix}')
        axes[i].set_ylabel(f'VAF_PSEUDO', fontsize=6)
        axes[i].set_yscale('log')
        axes[i].set_xscale('log')
        axes[i].set_title(f"{weight_prop}|{weight}\nrho={rho:.2f}, \nprop_mut_tissue={dg['VAF_PSEUDO'].sum():.2f}", fontsize=6)
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()
    plt.show()
    


@click.command()
@click.option('--sample_name', type=str, required=True, help='Name of the sample')
@click.option('--maf_file', type=click.Path(exists=True), required=False, help='MAF file with mutations')
@click.option('--output_prefix', type=str, required=True, help='Prefix for output files')
@click.option('--max_n', type=int, default=10, help='Maximum N for hyperbolic curves')
def main(sample_name, maf_file, output_prefix, max_n):
    """
    Generate VAF and mutation density vs depth plots with hyperbolic curves.
    
    Creates scatter plots showing relationships between depth and:
    - VAF (variant allele frequency) per site
    - Mutation density per gene (with different mutation type flavors)
    
    Overlays hyperbolic curves N/depth and reports mutation counts.
    """
    output_pdf_path = f"{output_prefix}.mutations_vaf.pdf"
    
    # Plot VAF vs depth per site
    if maf_file:
        print(f"Generating VAF vs depth plot from {maf_file}")
        maf_df = pd.read_csv(maf_file, sep='\t', na_values=custom_na_values)
        if maf_df.shape[0] < 5:
            print(f"There are less than 5 mutations in MAF file {maf_file}. Skipping VAF vs depth plot.")
        else:
            with PdfPages(output_pdf_path) as pdf:
                plot_vaf_vs_depth_per_site(maf_df, pdf, sample_name, max_n=max_n)
                print(f"VAF vs depth plot complete")
                plot_vaf_depth_heatmap(maf_df, pdf, sample_name)
                print(f"VAF vs depth heatmap complete")
                plot_vaf_vs_vafam_histogram(maf_df, pdf)
                print(f"VAF vs VAF_AM histogram complete")
                plot_vaf_pseudocount_curve(maf_df, pdf, suffix='')
                print(f"VAF pseudocount vs depth plot complete")
                plot_vaf_pseudocount_curve(maf_df, pdf, suffix='_AM')
                print(f"VAF AM pseudocount vs depth plot complete")

            print(f"Plots saved to {output_pdf_path}")

    get_top_mutations(maf_file, output_prefix)
    if len(maf_df["SAMPLE_ID"].unique()) > 1:
        get_top_mutations_all_samples(maf_file, output_prefix)


if __name__ == '__main__':
    try :
        main()
    except Exception as e:
        print(f"An error occurred: {e}")
