#!/usr/bin/env python


import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from read_utils import custom_na_values


def plot_single_needle(gene, snvs_maf_obs,
                        seq_info_df,
                        sample = ''
                        ):

    snvs_maf_obs_gene = snvs_maf_obs[(snvs_maf_obs["canonical_SYMBOL"] == gene)
                                        & (snvs_maf_obs["TYPE"] == "SNV")]

    # Subset gene
    seq_info_df_gene = seq_info_df[seq_info_df["Gene"] == gene]
    gene_len = len(seq_info_df_gene.Seq.values[0])
    pos_gene = pd.DataFrame({"Pos": range(1, gene_len + 1)})

    # Get per-position SNV mutations count
    obs_snv_count_gene = snvs_maf_obs_gene.groupby("canonical_Protein_position").size().reset_index(name='Count').astype(int)
    obs_snv_count_gene.columns = ["Pos", "Count"]
    obs_snv_count_gene = pos_gene.merge(obs_snv_count_gene, how="left", on="Pos")


    # Determine the maximum y-value for this gene
    max_y_gene = obs_snv_count_gene["Count"].max()


    # Create a single figure with two subplots
    fig, axs = plt.subplots(1, 1, figsize=(10, 3), sharex=True)


    # Plot the data for observed and randomized mutations in separate subplots
    plotting_needle_from_counts(obs_snv_count_gene, axs, max_y=max_y_gene)

    # Set titles for each subplot
    axs.set_title(f'{gene} {sample} - {int(obs_snv_count_gene["Count"].sum()):,} observed SNVs')

    # Adjust layout to prevent overlap
    plt.tight_layout()

    return fig


def plotting_needle_from_counts(data_gene,
                                ax = None,
                                max_y=None,
                                col_pos_track='#003366',
                                label_pos_track='observed',
                                col_hv_lines='grey'):

    # Determine max_y if not provided
    if max_y is None:
        max_y = data_gene["Count"].max()

    # Calculate the precise margin to add to the y-axis
    marker_size = 60
    marker_radius = np.sqrt(marker_size / np.pi)
    marker_margin = marker_radius * 0.5  # Convert marker radius to a suitable margin

    # Add the margin to max_y
    max_y += marker_margin


    # If ax is not provided, create a new figure and axis
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 3))
    else:
        fig = None  # No need to create a figure if ax is provided


    ax.vlines(data_gene["Pos"], ymin=0, ymax=data_gene["Count"], color=col_hv_lines, lw=1, zorder=1, alpha=0.5)
    ax.scatter(data_gene["Pos"], data_gene["Count"], color='white', zorder=3, lw=1, ec="white")  # To cover the overlapping needle top part
    ax.scatter(data_gene["Pos"].values, data_gene["Count"].values, color=col_pos_track, zorder=4,
                alpha=0.7, lw=0.1, ec="black", s=30, label=label_pos_track)

    # Remove the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Set the y-axis limit
    ax.set_ylim(0, max_y)

    # Add labels
    ax.set_xlabel('Position')
    ax.set_ylabel('Count')

    # If a new figure was created, show it
    if fig is not None:
        plt.show()

    return fig


def manager(mutations_file, o3d_seq_file, sample_name, sample_name_out):
    # Load your MAF DataFrame (raw_annotated_maf)
    maf = pd.read_csv(mutations_file, sep = "\t", header = 0, na_values = custom_na_values)
    o3d_seq_df = pd.read_csv(o3d_seq_file, sep = "\t", header = 0)

    maf_f = maf[maf["canonical_Protein_position"] != '-'].reset_index(drop = True)
    del maf

    gene_order = sorted(pd.unique(maf_f["canonical_SYMBOL"]))


    os.makedirs(f"{sample_name_out}")

    # Loop over each gene to plot
    for geneeee in gene_order:
        print(geneeee)
        try :
            fig_needles = plot_single_needle(geneeee, maf_f,
                                            o3d_seq_df,
                                            sample=sample_name
                                            )
            fig_needles.savefig(f"{sample_name_out}/{geneeee}.needles.pdf", bbox_inches='tight')
            plt.close()

        except Exception as exe:
            print(geneeee)
            print(exe)





# @click.command()
# @click.option('--sample_name', type=str, help='Name of the sample being processed.')
# @click.option('--mut_file', type=click.Path(exists=True), help='Input mutation file')
# @click.option('--out_maf', type=click.Path(), help='Output MAF file')
# @click.option('--json_filters', type=click.Path(exists=True), help='Input mutation filtering criteria file')
# @click.option('--req_plots', type=click.Path(exists=True), help='Column names to output')
# # @click.option('--plot', is_flag=True, help='Generate plot and save as PDF')

# def main(sample_name, mut_file, out_maf, json_filters, req_plots): # , plot):
#     click.echo(f"Subsetting MAF file...")
#     subset_mutation_dataframe(sample_name, mut_file, out_maf, json_filters, req_plots)

# if __name__ == '__main__':
#     main()


sample_name_  = sys.argv[1]
mut_file     = sys.argv[2]
o3d_seq_file_ = sys.argv[3]
sample_name_out_ = sys.argv[4]
# out_maf      = sys.argv[3]
# json_filters = sys.argv[4]
# req_plots    = sys.argv[5]



if __name__ == '__main__':
    # maf = subset_mutation_dataframe(mut_file, json_filters)
    # plot_manager(sample_name, maf, req_plots)
    manager(mut_file, o3d_seq_file_, sample_name_, sample_name_out_)
