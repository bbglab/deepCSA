#!/usr/bin/env python

import click
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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



def get_counts_per_position_n_consequence(somatic_maf_file, expanded = False):
    somatic_maf = pd.read_table(somatic_maf_file, na_values = custom_na_values)

    somatic_maf_clean = somatic_maf[(somatic_maf["TYPE"] == 'SNV')
                                    & (~somatic_maf["FILTER.not_in_exons"])
                                    & (somatic_maf['canonical_Protein_position'] != '-')
                                    ].reset_index(drop = True)
    somatic_maf_clean['canonical_Protein_position'] = somatic_maf_clean['canonical_Protein_position'].astype(int)
    if expanded:
        somatic_maf_clean.loc[somatic_maf_clean["ALT_DEPTH"] > 5, "ALT_DEPTH"] = 5
        counts_per_position = somatic_maf_clean.groupby(by = ["SAMPLE_ID", "canonical_SYMBOL", 'canonical_Consequence_broader', 'canonical_Protein_position'])['ALT_DEPTH'].sum().to_frame('Count').reset_index()
        click.echo("Using expanded counts")

    else:
        counts_per_position = somatic_maf_clean.groupby(by = ["SAMPLE_ID", "canonical_SYMBOL", 'canonical_Consequence_broader', 'canonical_Protein_position'])['ALT_DEPTH'].size().to_frame('Count').reset_index()

    counts_per_position.columns = ["SAMPLE_ID", 'Gene', 'Consequence', 'Pos', 'Count']

    return counts_per_position


def plot_count_track(count_df,
                     gene_start,
                     gene_end,
                        axes,
                        colors_dict,
                        ax=0,
                        alpha=1,
                        indel=False,
                        n_batches = 10
                    ):
    """
    Plots needle plot of mutation counts by position.

    Parameters:
        count_df: DataFrame with ['Pos', 'Consequence', 'Count'] columns
        gene_start: Start position of the gene
        gene_end: End position of the gene
        axes: matplotlib axes array
        colors_dict: dictionary mapping consequence -> color
        ax: index of subplot
        alpha: transparency
        indel: whether to include 'indel' consequence
        n_batches: number of batches to split the data into for plotting

    """
    assert gene_start < gene_end, "gene_start must be less than gene_end"

    gene_len = gene_end - gene_start + 1
    # Shuffle the data and split into batches
    shuffled_df = count_df.sample(frac=1, random_state=42).reset_index(drop=True)
    batches = np.array_split(shuffled_df, n_batches)

    legend_list = []
    pos_df = pd.DataFrame({"Pos" : range(gene_start, gene_end + 1)})

    for batch_idx, batch in enumerate(batches):
        for cnsq in ['nonsense', 'missense', 'synonymous']:
            if indel == False and cnsq == "indel":
                continue
            count_cnsq_df = batch[batch["Consequence"] == cnsq].reset_index(drop=True)
            count_cnsq_df = pos_df.merge(count_cnsq_df, on="Pos", how="left")
            axes[ax].vlines(count_cnsq_df["Pos"], ymin=0, ymax=count_cnsq_df["Count"], lw=1, zorder=1, alpha=0.5, color=colors_dict["hv_lines_needle"])
            axes[ax].scatter(count_cnsq_df["Pos"], count_cnsq_df["Count"], s=20, color='white', zorder=3, lw=0.1, ec="none") # To cover the overlapping needle top part
            if cnsq not in legend_list:
                axes[ax].scatter(count_cnsq_df["Pos"].values, count_cnsq_df["Count"].values, zorder=4,
                                    alpha=alpha, lw=0.1, ec="none", s=20, label= "Truncating" if cnsq == 'nonsense' else cnsq.capitalize(), color=colors_dict[cnsq])
                legend_list.append(cnsq)
            else:
                axes[ax].scatter(count_cnsq_df["Pos"].values, count_cnsq_df["Count"].values, zorder=4,
                                    alpha=alpha, lw=0.1, ec="none", s=20, color=colors_dict[cnsq])

    axes[ax].spines['right'].set_visible(False)
    axes[ax].spines['top'].set_visible(False)
    axes[ax].set_ylabel("Mutation count")
    axes[ax].set_xlabel("Protein position")

    # Add right Y axis with proportion labels
    total_mutations = count_df["Count"].sum()
    if total_mutations > 0:
        ax_right = axes[ax].twinx()
        left_ticks = axes[ax].get_yticks()
        right_ticks = left_ticks / total_mutations
        ax_right.set_ylim(axes[ax].get_ylim())
        ax_right.set_yticks(left_ticks)
        ax_right.set_yticklabels([f"{rt:.2f}" for rt in right_ticks])
        ax_right.set_ylabel("Proportion of mutations")
    else:
        ax_right = axes[ax].twinx()
        ax_right.set_ylim(axes[ax].get_ylim())
        ax_right.set_yticks(axes[ax].get_yticks())
        ax_right.set_yticklabels(["0.00"] * len(axes[ax].get_yticks()))
        ax_right.set_ylabel("Proportion of mutations")

    axes[ax].set_xlim(gene_start - gene_len * 0.02, gene_end + gene_len * 0.02)

def plot_stacked_bar_track_binned(count_df,
                                    gene_start,
                                    gene_end,
                                    axes,
                                    colors_dict,
                                    ax=0,
                                    alpha=1,
                                    indel=False,
                                    min_bin_size=3,
                                    num_bins = 100):
    """
    Plots stacked barplot of mutation counts binned by position.

    Parameters:
        count_df: DataFrame with ['Pos', 'Consequence', 'Count'] columns
        gene_start: Start position of the gene
        gene_end: End position of the gene
        axes: matplotlib axes array
        colors_dict: dictionary mapping consequence -> color
        ax: index of subplot
        alpha: transparency
        indel: whether to include 'indel' consequence
        min_bin_size: minimum size of non-overlapping bins
        num_bins: number of bins to create
    """
    assert gene_start < gene_end, "gene_start must be less than gene_end"

    gene_len = gene_end - gene_start + 1

    # Compute bin_size or result to default
    candidate_bin_size = max(1, gene_len // num_bins)
    bin_size = max(min_bin_size, candidate_bin_size)

    valid_consequences = ['nonsense', 'missense', 'synonymous']
    if indel:
        valid_consequences.append('indel')

    filtered_df = count_df[count_df["Consequence"].isin(valid_consequences)].copy()

    # Assign bin start
    filtered_df["Bin"] = ((filtered_df["Pos"] - 1) // bin_size) * bin_size + 1

    # Group and pivot
    binned_df = (
        filtered_df
        .groupby(["Bin", "Consequence"])["Count"]
        .sum()
        .unstack(fill_value=0)
        .reindex(columns=valid_consequences, fill_value=0)
    )

    # Ensure all bins are represented
    all_bins = list(range((gene_start-1)//bin_size * bin_size + 1, gene_end + 1, bin_size))
    binned_df = binned_df.reindex(all_bins, fill_value=0)

    # Plot stacked bars
    bottom = np.zeros(len(binned_df))
    for cnsq in valid_consequences:
        axes[ax].bar(
            binned_df.index,
            binned_df[cnsq],
            bottom=bottom,
            width=bin_size * 0.8,
            align="edge",
            color=colors_dict.get(cnsq, 'gray'),
            alpha=alpha,
            label="Truncating" if cnsq == 'nonsense' else cnsq.capitalize(),
            linewidth=0
        )
        bottom += binned_df[cnsq].values


    # Clean up plot
    axes[ax].spines['right'].set_visible(False)
    axes[ax].spines['top'].set_visible(False)
    axes[ax].set_ylabel(f"Mutation count\n({bin_size} AA bin)")
    axes[ax].set_xlabel("Protein position")

    # Add right Y axis with proportion labels
    total_mutations = count_df["Count"].sum()
    if total_mutations > 0:
        ax_right = axes[ax].twinx()
        left_ticks = axes[ax].get_yticks()
        right_ticks = left_ticks / total_mutations
        ax_right.set_ylim(axes[ax].get_ylim())
        ax_right.set_yticks(left_ticks)
        ax_right.set_yticklabels([f"{rt:.2f}" for rt in right_ticks])
        ax_right.set_ylabel("Proportion of mutations")
    else:
        ax_right = axes[ax].twinx()
        ax_right.set_ylim(axes[ax].get_ylim())
        ax_right.set_yticks(axes[ax].get_yticks())
        ax_right.set_yticklabels(["0.00"] * len(axes[ax].get_yticks()))
        ax_right.set_ylabel("Proportion of mutations")

    axes[ax].set_xlim(gene_start - gene_len * 0.02, gene_end + gene_len * 0.02)


def manager(mutations_file, outdir, expand_muts = False):

    counts_per_position = get_counts_per_position_n_consequence(mutations_file, expanded=expand_muts)

    gene_order = sorted(pd.unique(counts_per_position["Gene"]))


    # Loop over each gene to plot
    for gene in gene_order:
        print(gene)
        try:
            mut_count_df = counts_per_position[(counts_per_position["Gene"] == gene)]
            mut_count_df = mut_count_df.groupby(by=["Pos", "Consequence"])["Count"].sum().reset_index()

            fig, ax = plt.subplots(1, 1, figsize=(5, 1.2))
            plot_count_track(
                mut_count_df,
                gene_start=mut_count_df["Pos"].min(),
                gene_end=mut_count_df["Pos"].max(),
                axes=[ax], ax=0,
                colors_dict=metrics_colors_dictionary, indel=False, alpha=0.7
            )
            ax.set_title(f"{gene}")
            plt.savefig(f"{outdir}/{gene}{'_expanded' if expand_muts else ''}.needle.pdf", bbox_inches='tight', dpi=100)
            plt.show()
            plt.close()

        except Exception as exe:
            print(gene)
            print(exe)

        # stacked version
        try:
            mut_count_df = counts_per_position[(counts_per_position["Gene"] == gene)]
            mut_count_df = mut_count_df.groupby(by=["Pos", "Consequence"])["Count"].sum().reset_index()

            fig, ax = plt.subplots(1, 1, figsize=(5, 1.2))
            plot_stacked_bar_track_binned(
                count_df=mut_count_df,
                gene_start=mut_count_df["Pos"].min(),
                gene_end=mut_count_df["Pos"].max(),
                axes=[ax], ax=0,
                colors_dict=metrics_colors_dictionary,
                alpha=1,
                indel=False
            )
            ax.set_title(f"{gene}")
            plt.savefig(f"{outdir}/{gene}{'_expanded' if expand_muts else ''}.stacked.pdf", bbox_inches='tight', dpi=100)
            plt.show()
            plt.close()

        except Exception as exe:
            print(gene)
            print(exe)





@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--mut_file', type=click.Path(exists=True), help='Input mutations file')
@click.option('--o3d_seq_file', type=click.Path(exists=True), help='Input Oncodrive3D sequence df file')
@click.option('--outdir', type=click.Path(), help='Output path for plots')
def main(sample_name, mut_file, o3d_seq_file, outdir):
    click.echo("Generating simple needle plots...")
    manager(mut_file, outdir)

    click.echo("Generating expanded needle plots...")
    manager(mut_file, outdir, expand_muts=True)

if __name__ == '__main__':
    main()
