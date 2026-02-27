#!/usr/bin/env python

import sys
import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from utils import contexts_no_change


def compute_trinucleotides(sample_name, depths_file, pseudocount = 0):
    """
    Compute mutational profile from the input data
          ***Remember to add some pseudocounts to the computation***

        Required information:
            Annotated mutations observed
        Output:
            Mutational profile per sample, possibility to add pseudocounts to prevent some probabilities from being 0
    """

    # Load your MAF DataFrame (raw_annotated_maf)
    depths_annotated = pd.read_csv(depths_file, sep = "\t", header = 0)

    # create the matrix in the desired order
    empty_matrix = pd.DataFrame(index = contexts_no_change)

    samples = [ x for x in depths_annotated.columns if x not in ["CHROM", "POS", "CONTEXT"] ]
    trinucleotides_per_sample = depths_annotated.groupby(by = "CONTEXT")[samples].sum().fillna(0)
    trinucleotides_per_sample = pd.concat( (empty_matrix, trinucleotides_per_sample) , axis = 1)
    if '-' in trinucleotides_per_sample.index:
        trinucleotides_per_sample = trinucleotides_per_sample.drop("-", axis = 0)
    trinucleotides_per_sample.index.name = "CONTEXT"
    trinucleotides_per_sample = trinucleotides_per_sample.reset_index(drop = False)
    trinucleotides_per_sample.columns = ["CONTEXT", sample_name]

    trinucleotides_per_sample[sample_name] = trinucleotides_per_sample[sample_name] + pseudocount


    trinucleotides_per_sample[["CONTEXT", sample_name]].to_csv(f"{sample_name}.trinucleotides.tsv.gz",
                                                                sep = "\t",
                                                                header = True,
                                                                index = False)

    return trinucleotides_per_sample[["CONTEXT", sample_name]]


def plot_trinucleotide_proportions(wgs_counts_file, sample_trinucleotides, sample_name):
    wgs_counts = pd.read_table(wgs_counts_file)
    wgs_counts.columns = ['CONTEXT', 'COUNT_WGS']


    counts_all = wgs_counts.copy()
    counts_all = counts_all.merge(sample_trinucleotides, on = 'CONTEXT')
    counts_all = counts_all.set_index("CONTEXT")
    proportions_all = counts_all / counts_all.sum()

    proportions_all_plot = proportions_all.copy()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
       
    rmse = np.sqrt(((proportions_all_plot["COUNT_WGS"] - proportions_all_plot[sample_name])**2).mean())
    
    # Scatter plot
    sns.scatterplot(data=proportions_all_plot,
                    x="COUNT_WGS",
                    y=sample_name,
                    hue="CONTEXT",
                    legend=False,
                    ax=ax)
    
    # Identity line (x = y)
    lims = [
        min(proportions_all_plot["COUNT_WGS"].min(), proportions_all_plot[sample_name].min()),
        max(proportions_all_plot["COUNT_WGS"].max(), proportions_all_plot[sample_name].max()),
    ]
    ax.plot(lims, lims, '--', color='gray')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    # Annotate points with CONTEXT
    for namee, row in proportions_all_plot.iterrows():
        ax.text(row["COUNT_WGS"], row[sample_name], namee,
                fontsize=6, alpha=0.7)
    
    ax.set_xlabel("Proportion WGS")
    ax.set_ylabel(f"Proportion {sample_name}")
    ax.set_title(f"{sample_name} (RMSE: {rmse:.4f})")

    plt.tight_layout()
    plt.savefig(f"{sample_name}.trinucleotide_content_comparisons.proportions.png", dpi=100)
    plt.show()



@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--depths_file', type=click.Path(exists=True), help='Input depths file')
@click.option('--ref_wgs_trinucleotides', default = None, type=click.Path(exists=True), help='File with trinucleotide counts for the whole-genome')
@click.option('--pseud', type=float, default=0.5)
def main(sample_name, depths_file, ref_wgs_trinucleotides, pseud):
    click.echo(f"Running the trinucleotide computation...")
    click.echo(f"Using the pseudocount: {pseud}")
    sample_trinuc = compute_trinucleotides(sample_name, depths_file, pseud )
    click.echo("Trinucleotides computation completed.")
    
    click.echo("Plotting trinucleotides.")
    if ref_wgs_trinucleotides is not None:
        plot_trinucleotide_proportions(ref_wgs_trinucleotides, sample_trinuc, sample_name)

if __name__ == '__main__':
    main()

