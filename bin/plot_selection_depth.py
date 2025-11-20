#!/usr/bin/env python

"""
Plot selection metrics (Omega, OncodriveFML) vs Depth.

This script generates scatter plots showing the relationship between sequencing depth
and selection metrics without hyperbolic curves.
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


def plot_omega_vs_depth(omega_df, depth_df, output_pdf, sample_name):
    """
    Plot omega values vs average depth per gene.
    
    Parameters:
    -----------
    omega_df : DataFrame
        Omega dataframe with gene-level omega values
    depth_df : DataFrame
        Depth dataframe with per-gene average depths
    output_pdf : PdfPages
        PDF file to save plots
    sample_name : str
        Name of the sample
    """
    # Merge omega with depth information
    plot_data = omega_df.merge(depth_df, left_on='gene', right_on='GENE', how='inner')
    
    # Get different impact types
    impacts = plot_data['impact'].unique() if 'impact' in plot_data.columns else []
    
    for impact in impacts:
        impact_data = plot_data[plot_data['impact'] == impact]
        impact_data = impact_data[(impact_data['MEAN_GENE_DEPTH'] > 0) & (impact_data['dnds'] > 0)]
        
        if impact_data.empty:
            continue
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Color by significance
        colors = ['red' if p < 0.05 else 'gray' for p in impact_data['pvalue']]
        
        # Scatter plot
        ax.scatter(impact_data['MEAN_GENE_DEPTH'], impact_data['dnds'], 
                  alpha=0.6, s=40, c=colors, edgecolors='black', linewidths=0.5)
        
        # Add reference line at omega = 1 (neutral)
        ax.axhline(y=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Neutral (ω=1)')
        
        ax.set_xlabel('Average Gene Depth (reads)', fontsize=plots_general_config["xlabel_fontsize"])
        ax.set_ylabel('dN/dS (ω)', fontsize=plots_general_config["ylabel_fontsize"])
        ax.set_title(f'{sample_name} - Omega vs Depth ({impact}, N={len(impact_data)} genes)', 
                    fontsize=plots_general_config["title_fontsize"])
        
        max_depth = impact_data['MEAN_GENE_DEPTH'].quantile(0.99)
        ax.set_xlim(0, max_depth)
        
        # Add legend for colors
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='red', label='p < 0.05'),
                          Patch(facecolor='gray', label='p ≥ 0.05')]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=plots_general_config["legend_fontsize"])
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_pdf.savefig()
        plt.close()


def plot_oncodrivefml_vs_depth(ofml_df, depth_df, output_pdf, sample_name):
    """
    Plot OncodriveFML scores vs average depth per gene.
    
    Parameters:
    -----------
    ofml_df : DataFrame
        OncodriveFML dataframe with gene-level scores
    depth_df : DataFrame
        Depth dataframe with per-gene average depths
    output_pdf : PdfPages
        PDF file to save plots
    sample_name : str
        Name of the sample
    """
    # Merge OncodriveFML with depth information
    plot_data = ofml_df.merge(depth_df, left_on='SYMBOL', right_on='GENE', how='inner')
    plot_data = plot_data[(plot_data['MEAN_GENE_DEPTH'] > 0)]
    
    if plot_data.empty:
        print(f"No valid data for OncodriveFML vs depth plot")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter plot - color by significance if available
    if 'QVALUE' in plot_data.columns:
        colors = ['red' if q < 0.1 else 'gray' for q in plot_data['QVALUE']]
    else:
        colors = 'steelblue'
    
    ax.scatter(plot_data['MEAN_GENE_DEPTH'], plot_data['SCORE'], 
              alpha=0.6, s=40, c=colors, edgecolors='black', linewidths=0.5)
    
    ax.set_xlabel('Average Gene Depth (reads)', fontsize=plots_general_config["xlabel_fontsize"])
    ax.set_ylabel('OncodriveFML Score', fontsize=plots_general_config["ylabel_fontsize"])
    ax.set_title(f'{sample_name} - OncodriveFML vs Depth (N={len(plot_data)} genes)', 
                fontsize=plots_general_config["title_fontsize"])
    
    max_depth = plot_data['MEAN_GENE_DEPTH'].quantile(0.99)
    ax.set_xlim(0, max_depth)
    
    # Add legend for colors if q-values available
    if 'QVALUE' in plot_data.columns:
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='red', label='q < 0.1'),
                          Patch(facecolor='gray', label='q ≥ 0.1')]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=plots_general_config["legend_fontsize"])
    
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()


@click.command()
@click.option('--sample_name', type=str, required=True, help='Name of the sample')
@click.option('--omega_file', type=click.Path(exists=True), required=False, help='Omega (dN/dS) results file')
@click.option('--oncodrivefml_file', type=click.Path(exists=True), required=False, help='OncodriveFML results file')
@click.option('--depth_file', type=click.Path(exists=True), required=True, help='Depth per gene file')
@click.option('--output_prefix', type=str, required=True, help='Prefix for output files')
def main(sample_name, omega_file, oncodrivefml_file, depth_file, output_prefix):
    """
    Generate selection metric vs depth plots.
    
    Creates scatter plots showing relationships between depth and:
    - Omega (dN/dS) per gene
    - OncodriveFML scores per gene
    
    No hyperbolic curves are added to these plots.
    """
    output_pdf_path = f"{output_prefix}.selection_depth.pdf"
    
    # Load depth data
    depth_df = pd.read_csv(depth_file, sep='\t', na_values=custom_na_values)
    
    with PdfPages(output_pdf_path) as pdf:
        # Plot omega vs depth per gene
        if omega_file:
            print(f"Generating omega vs depth plots")
            omega_df = pd.read_csv(omega_file, sep='\t', na_values=custom_na_values)
            plot_omega_vs_depth(omega_df, depth_df, pdf, sample_name)
        
        # Plot OncodriveFML vs depth per gene
        if oncodrivefml_file:
            print(f"Generating OncodriveFML vs depth plots")
            ofml_df = pd.read_csv(oncodrivefml_file, sep='\t', na_values=custom_na_values)
            plot_oncodrivefml_vs_depth(ofml_df, depth_df, pdf, sample_name)
    
    print(f"Plots saved to {output_pdf_path}")


if __name__ == '__main__':
    main()
