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
                   color='gray', label=f'{n}/depth' if n <= 3 else '')
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


def plot_vaf_vs_depth_per_site(maf_df, output_pdf, sample_name, max_n=10):
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
    plot_data = maf_df[['DEPTH', 'VAF', 'ALT_DEPTH']].dropna()
    plot_data = plot_data[(plot_data['DEPTH'] > 0) & (plot_data['VAF'] >= 0) & (plot_data['VAF'] <= 1.0)]
    
    if plot_data.empty:
        print(f"No valid data for VAF vs depth plot")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter plot
    ax.scatter(plot_data['DEPTH'], plot_data['VAF'], 
              alpha=0.5, s=20, edgecolors='none', c='steelblue')
    
    # Add hyperbolic curves
    max_depth = plot_data['DEPTH'].quantile(0.99)
    add_hyperbolic_curves(ax, max_depth, max_n=max_n)
    
    # Count mutations on curves
    counts = count_mutations_on_curves(
        plot_data['DEPTH'].values, 
        plot_data['VAF'].values, 
        max_n=max_n
    )
    
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
    ax.set_xlim(0, max_depth)
    ax.set_ylim(0, 1.0)
    ax.legend(loc='upper right', fontsize=plots_general_config["legend_fontsize"])
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()


def plot_mutdensity_vs_depth(mutdensity_df, depth_df, output_pdf, sample_name, max_n=10):
    """
    Plot mutation density vs average depth per gene.
    
    Parameters:
    -----------
    mutdensity_df : DataFrame
        Mutation density dataframe with per-gene information
    depth_df : DataFrame
        Depth dataframe with per-gene average depths
    output_pdf : PdfPages
        PDF file to save plots
    sample_name : str
        Name of the sample
    max_n : int
        Maximum N for hyperbolic curves
    """
    # Merge mutation density with depth information
    plot_data = mutdensity_df.merge(depth_df, on='GENE', how='inner')
    
    # Filter to exclude ALL_GENES summary
    plot_data = plot_data[plot_data['GENE'] != 'ALL_GENES']
    
    # Get different mutation type flavors
    mut_types = plot_data['MUTTYPES'].unique() if 'MUTTYPES' in plot_data.columns else ['all_types']
    
    for mut_type in mut_types:
        type_data = plot_data[plot_data['MUTTYPES'] == mut_type] if 'MUTTYPES' in plot_data.columns else plot_data
        type_data = type_data[(type_data['MEAN_GENE_DEPTH'] > 0) & (type_data['MUTDENSITY_MB'] >= 0)]
        
        if type_data.empty:
            continue
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Scatter plot
        ax.scatter(type_data['MEAN_GENE_DEPTH'], type_data['MUTDENSITY_MB'], 
                  alpha=0.6, s=40, edgecolors='black', linewidths=0.5, c='coral')
        
        # For mutation density, we need to adapt the hyperbolic curves
        # mutation_density ≈ N_mutations / depth
        # If we normalize: density_normalized = mutations / (depth / 1e6) = mutations * 1e6 / depth
        max_depth = type_data['MEAN_GENE_DEPTH'].quantile(0.99)
        
        # Add reference curves for mutation density
        # These represent densities that would result from N mutations at various depths
        depth_range = np.linspace(1, max_depth, 1000)
        
        # Calculate reasonable N values based on data
        max_density = type_data['MUTDENSITY_MB'].quantile(0.99)
        for n_muts in [1, 2, 5, 10, 20, 50]:
            density_curve = (n_muts / depth_range) * 1e6
            if density_curve[0] <= max_density * 1.5:
                ax.plot(depth_range, density_curve, '--', alpha=0.3, 
                       linewidth=1, color='gray', 
                       label=f'{n_muts} muts' if n_muts in [1, 2, 5, 10] else '')
        
        ax.set_xlabel('Average Gene Depth (reads)', fontsize=plots_general_config["xlabel_fontsize"])
        ax.set_ylabel('Mutation Density (per Mb)', fontsize=plots_general_config["ylabel_fontsize"])
        
        type_label = mut_type.replace('-', ', ')
        ax.set_title(f'{sample_name} - Mutation Density vs Depth\n({type_label}, N={len(type_data)} genes)', 
                    fontsize=plots_general_config["title_fontsize"])
        ax.set_xlim(0, max_depth)
        ax.set_ylim(0, max_density * 1.1)
        ax.legend(loc='upper right', fontsize=plots_general_config["legend_fontsize"])
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_pdf.savefig()
        plt.close()


@click.command()
@click.option('--sample_name', type=str, required=True, help='Name of the sample')
@click.option('--maf_file', type=click.Path(exists=True), required=False, help='MAF file with mutations')
@click.option('--mutdensity_file', type=click.Path(exists=True), required=False, help='Mutation density file')
@click.option('--depth_file', type=click.Path(exists=True), required=False, help='Depth per gene file')
@click.option('--output_prefix', type=str, required=True, help='Prefix for output files')
@click.option('--max_n', type=int, default=10, help='Maximum N for hyperbolic curves')
def main(sample_name, maf_file, mutdensity_file, depth_file, output_prefix, max_n):
    """
    Generate VAF and mutation density vs depth plots with hyperbolic curves.
    
    Creates scatter plots showing relationships between depth and:
    - VAF (variant allele frequency) per site
    - Mutation density per gene (with different mutation type flavors)
    
    Overlays hyperbolic curves N/depth and reports mutation counts.
    """
    output_pdf_path = f"{output_prefix}.vaf_mutdensity_depth.pdf"
    
    with PdfPages(output_pdf_path) as pdf:
        # Plot VAF vs depth per site
        if maf_file:
            print(f"Generating VAF vs depth plot from {maf_file}")
            maf_df = pd.read_csv(maf_file, sep='\t', na_values=custom_na_values)
            plot_vaf_vs_depth_per_site(maf_df, pdf, sample_name, max_n=max_n)
        
        # Plot mutation density vs depth per gene
        if mutdensity_file and depth_file:
            print(f"Generating mutation density vs depth plots")
            mutdensity_df = pd.read_csv(mutdensity_file, sep='\t', na_values=custom_na_values)
            depth_df = pd.read_csv(depth_file, sep='\t', na_values=custom_na_values)
            plot_mutdensity_vs_depth(mutdensity_df, depth_df, pdf, sample_name, max_n=max_n)
    
    print(f"Plots saved to {output_pdf_path}")


if __name__ == '__main__':
    main()
