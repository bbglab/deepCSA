#!/usr/bin/env python

import click
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

from utils_context import triplet_context_iterator
from utils_impacts import broadimpact_grouping_dict_with_synonymous


# Utils
# define canonical order for trinucleotide contexts

triplet_contexts = triplet_context_iterator()


def get_correction_factor(sample_name, trinucleotide_counts_df, mutability_df, flat=False):
    """
    The following function computes the normalization factor so that the mutation density
        that we will compute downstream matches the equivalent to 1 mut/(Mb-genome)
        on average across the genome.
    """

    # vector of triplet count in 96-channel canonical sorting
    triplet_counts_dict = dict(zip(trinucleotide_counts_df['CONTEXT'].values, trinucleotide_counts_df['COUNT'].values))
    l = [triplet_counts_dict[c[:3]] for c in triplet_contexts]
    triplet_counts = np.array(l)

    # genome length in Mb
    genome_length = sum(triplet_counts) / (3 * 1e6)

    # vector of relative mutabilities in 96-channel canonical sorting
    if flat:
        relative_mutability = (1 / 96) * np.ones(96)
    else:
        d = dict(zip(mutability_df['CONTEXT_MUT'], mutability_df[f'{sample_name}.all']))
        relative_mutability = np.array([d.get(k, 0) for k in triplet_contexts])

    # return correction factor a.k.a. "alpha hat"
    correction_factor = genome_length / np.dot(relative_mutability, triplet_counts)
    
    return correction_factor, relative_mutability


def mutation_density(sample_name, depths_file, somatic_mutations_file, mutability_file, panel_file, trinucleotide_counts_file, flat=False):

    depths_df = pd.read_csv(depths_file, sep='\t')
    somatic_mutations_df = pd.read_csv(somatic_mutations_file, sep='\t')
    mutability_df = pd.read_csv(mutability_file, sep='\t')
    panel_df = pd.read_csv(panel_file, sep='\t')
    trinucleotide_counts_df = pd.read_csv(trinucleotide_counts_file, sep='\t')

    genes = panel_df['GENE'].unique()
    consequences = broadimpact_grouping_dict_with_synonymous.keys()

    res = pd.DataFrame(index=genes, columns=consequences)

    for csqn, csqn_set in broadimpact_grouping_dict_with_synonymous.items():
        
        for gene in panel_df['GENE'].unique():
            
            # compute vector of sum of depths per trinucleotide context
            # tailored to the specific gene-impact target
            region_df = panel_df[(panel_df['IMPACT'].isin(csqn_set)) & (panel_df['GENE'] == gene)].copy()

            # counting every position once
            dh = pd.merge(region_df[['CHROM', 'POS']],
                          depths_df[['CHROM', 'POS', 'CONTEXT', sample_name]],
                          on=['CHROM', 'POS'], how='left')
            depth_sum_df = dh.groupby(by='CONTEXT').agg({sample_name: 'sum'}).reset_index()
            depth_region_dict = dict(zip(depth_sum_df['CONTEXT'], depth_sum_df[sample_name]))
            depth_region_vect = np.array([depth_region_dict.get(k[:3], 0) for k in triplet_contexts])


            # compute correction factor "alpha hat"
            correction_factor, relative_mutability = get_correction_factor(sample_name, trinucleotide_counts_df, mutability_df, flat=flat)


            # compute effective length
            effective_length = correction_factor * np.dot(relative_mutability, depth_region_vect)

            try:
                assert(effective_length != 0)
            except AssertionError:
                res.loc[gene, csqn] = None
                continue
            
            # observed somatic mutations

            n = somatic_mutations_df[(somatic_mutations_df['IMPACT'].isin(csqn_set)) & (somatic_mutations_df['GENE'] == gene)].shape[0]

            res.loc[gene, csqn] = n / effective_length
    
    return res


def logfoldchange_plot(sample_name, df, df_flat):
    """
    Heatmap representing log-FC between mutability adjusted and flat mutation density.
    This can be useful as a diagnostic tool of the mutation density calculation.
    """

    dh = (df + 0.01) / (df_flat + 0.01)  # Include a pseudocount
    dh = np.log2(dh.astype(float))
    with PdfPages(f'{sample_name}.logfoldchangeplot.pdf') as pdf:

        plt.figure(figsize=(8, 6)) # Set the size of the plot
        sns.heatmap(
            dh,
            annot=False,       
            cmap='coolwarm',   # Choose a colormap (e.g., 'coolwarm', 'viridis', 'YlGnBu')
            fmt=".2f",         # Format the annotations to two decimal places
            linewidths=.5,     # Add lines between cells for better readability
            cbar=True,         # Show the color bar
            square=True,       # Make the cells square
            vmin=-1, vmax=1, 
            cbar_kws={'label': 'log2 fold-change'}
        )

        # Set the title and labels
        plt.title(f'{sample_name}\nlog2 fold-change mutability adjusted vs flat', fontsize=16)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)

        # Display the plot
        plt.tight_layout() # Adjust plot to ensure everything fits
        pdf.savefig()
        plt.close()


@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--depths_file', type=str, help='Depth per position in panel')
@click.option('--somatic_mutations_file', type=click.Path(exists=True), help='Observed somatic mutations')
@click.option('--mutability_file', type=click.Path(exists=True), help='Relative mutability per trinucleotide context')
@click.option('--panel_file', type=click.Path(exists=True), help='All possible SNVs with csqn type')
@click.option('--trinucleotide_counts_file', type=click.Path(exists=True), help='How many trinucleotides of each type are in the entire genome')
def main(sample_name, depths_file, somatic_mutations_file, mutability_file, panel_file, trinucleotide_counts_file):

    click.echo(f"Running the mutability adjusted mutation density estimation...")
    
    # main calculations
    res = mutation_density(sample_name, depths_file, somatic_mutations_file, mutability_file, panel_file, trinucleotide_counts_file, flat=False)
    res_flat = mutation_density(sample_name, depths_file, somatic_mutations_file, mutability_file, panel_file, trinucleotide_counts_file, flat=True)
    logfoldchange_plot(sample_name, res, res_flat)

    # save results
    res["SAMPLE"] = sample_name
    res_flat["SAMPLE"] = sample_name
    res[['SAMPLE'] + [col for col in res.columns if col != 'SAMPLE']].to_csv(f'{sample_name}.mutdensities.tsv', sep='\t')
    res_flat[['SAMPLE'] + [col for col in res_flat.columns if col != 'SAMPLE']].to_csv(f'{sample_name}.mutdensities_flat.tsv', sep='\t')



if __name__ == '__main__':
    main()