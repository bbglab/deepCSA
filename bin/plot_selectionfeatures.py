#!/usr/local/bin/python


import pandas as pd
import json
import os, sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

pd.set_option('display.max_columns', None)

# Suppress warnings
import warnings
pd.options.mode.chained_assignment = None
warnings.filterwarnings("ignore", message="FixedFormatter should only be used together with FixedLocator")


import pandas as pd
import os, sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import linregress,norm
import tabix
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
    maf = pd.read_csv(mutations_file, sep = "\t", header = 0)
    o3d_seq_df = pd.read_csv(o3d_seq_file, sep = "\t", header = 0)

    maf_f = maf[maf["canonical_Protein_position"] != '-'].reset_index(drop = True)
    del maf

    gene_order = sorted(pd.unique(maf_f["canonical_SYMBOL"]))


    os.makedirs(f"{sample_name_out}")

    # Loop over each gene to plot
    for geneeee in gene_order:
        print(geneeee)
        fig_needles = plot_single_needle(geneeee, maf_f,
                                        o3d_seq_df,
                                        sample=sample_name
                                        )

        fig_needles.savefig(f"{sample_name_out}/{geneeee}.needles.pdf", bbox_inches='tight')
        plt.close()




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


# Init
run_name = "all_samples"

deepcsa_run_dir = "/workspace/nobackup/bladder_ts/results/2024-06-20_deepCSA"
maf_file = os.path.join(deepcsa_run_dir, f"writemaf/{run_name}.filtered.tsv.gz")
o3d_datasets = "/workspace/nobackup/scratch/oncodrive3d/datasets_240506"
o3d_annotations = "/workspace/nobackup/scratch/oncodrive3d/annotations_240506"
fig3_data = "/workspace/projects/bladder_ts/notebooks/manuscript_figures/Fig3/data"

gene_order = ["KMT2D","EP300","ARID1A","CREBBP","NOTCH2","KMT2C","STAG2","RB1",
              "RBM10","KDM6A","TP53","FGFR3","CDKN1A","FOXQ1",
             # "PIK3CA","TERT"
             ]



# Oncodrive3D

o3d_seq_df = pd.read_table(f"{o3d_datasets}/seq_for_mut_prob.tsv")
o3d_annot_df = pd.read_csv(f"{o3d_annotations}/uniprot_feat.tsv", sep="\t")

o3d_prob = f"{deepcsa_run_dir}/oncodrive3d/run/{run_name}/{run_name}.miss_prob.processed.json"
o3d_prob = json.load(open(o3d_prob, encoding="utf-8"))

o3d_score = f"{deepcsa_run_dir}/oncodrive3d/run/{run_name}/{run_name}.3d_clustering_pos.csv"
o3d_score = pd.read_csv(o3d_score)[["Gene", "Pos", "Score", "Score_obs_sim", "pval", "C", "C_ext"]]


# MAF

maf_df = pd.read_csv(maf_file, sep = "\t")
maf_df_f, snv_df, trunc_df, synon_df, miss_df = preprocess_maf(maf_df)





for gene in gene_order:

    plot_pars = init_plot_pars()
    gene_len = len(o3d_seq_df[o3d_seq_df["Gene"] == gene].Seq.values[0])
    if gene_len > 400:
        # if gene_len > 2500:
        #     gene_len = 2500
        ref_ratio = 500 / gene_len
        fsize_x = round(gene_len / 26)
    else:
        ref_ratio = 1
        fsize_x = 15

    plot_pars["fsize"] = fsize_x, 10
    plot_pars["ofml_cbar_coord"] = (0.04, 0.36, 1.1*ref_ratio, 0.85)


    generate_plot(gene,
                  o3d_seq_df,
                  o3d_annot_df,
                  o3d_prob,
                  o3d_score,
                  maf_df_f,
                  snv_df,
                  trunc_df,
                  synon_df,
                  miss_df,
                  plot_pars,
                  output_path="Fig3_plots/Fig3b")







# Plot
# ====

def init_plot_pars():

    plot_pars = {"fsize"                : (20,10),
                 "hspace"               : 0.1,                    # General space between all tracks
                 "ofml_cbar_coord"      : (0.01, 0.35, 0.9, 1),   # Box to anchor coordinate for OncodriveFML color bar
                 "track_title_x_coord"  : 0.83,                   # x-coordinate (respect to protein len) for track txt title
                 "score_txt_x_coord"    : 1.13,                   # as track title but for track score txt
                 "track_title_fontsize" : 14,
                 "ylabel_fontsize"      : 13.5,
                 "xlabel_fontsize"      : 13.5,
                 "ylabel_pad"           : 38,
                 "ticksize"             : 10.5,
                 "legend_fontsize"      : 12,
                 "legend_frameon"       : True,
                 "dpi"                  : 300
                }
    plot_pars["colors"] = {"ofml"        : "viridis_r",
                           "omega_trunc" : "#FA5E32",
                           "omega_synon" : "#89E4A2",
                           "omega_miss"  : "#FABE4A",
                           "o3d_score"   : "#6DBDCC",
                           "o3d_cluster" : "skyblue",
                           "o3d_prob"    : "darkgray",
                           "frameshift"  : "#E4ACF4",
                           "inframe"     : "C5",
                           "hv_lines"    : "lightgray"      # General horizontal and vertical lines (e.g., needle plot vline)
                           }
    plot_pars["h_ratios"] = {"omega_trunc" : 0.3,
                             "omega_synon" : 0.3,
                             "omega_miss"  : 0.3,
                             "space1"      : 0.015,
                             "o3d"         : 0.5,
                             "space2"      : 0.015,
                             "ofml"        : 0.5,
                             "space3"      : 0.015,
                             "indels"      : 0.5,
                             "space4"      : 0.015,
                             "domain"      : 0.07
                            }

    return plot_pars


def plot_count_track(pos_track,
                     axes,
                     ax=0,
                     neg_track=None,
                     gene_len=None,
                     col_pos_track="gray",
                     col_neg_track="gray",
                     label_pos_track=None,
                     label_neg_track=None,
                     ymargin=None,
                     hv_lines="lightgray"):

    axes[ax].vlines(pos_track["Pos"], ymin=0, ymax=pos_track["Count"], color=hv_lines, lw=1, zorder=1, alpha=0.5)
    axes[ax].scatter(pos_track["Pos"], pos_track["Count"], color='white', zorder=3, lw=1, ec="white") # To cover the overlapping needle top part
    axes[ax].scatter(pos_track["Pos"].values, pos_track["Count"].values, color=col_pos_track, zorder=4,
                     alpha=0.7, lw=0.1, ec="black", s=60, label=label_pos_track)

    if isinstance(neg_track, pd.DataFrame):
        if isinstance(gene_len, int):
            axes[ax].hlines(0, xmin=0, xmax=gene_len, color="gray", lw=0.6, zorder=1)
        axes[ax].vlines(neg_track["Pos"], ymin=-neg_track["Count"], ymax=0, color=hv_lines, lw=1, zorder=1, alpha=0.5)
        axes[ax].scatter(neg_track["Pos"], -neg_track["Count"], color='white', zorder=3, lw=1, ec="white") # To cover the overlapping needle top part
        axes[ax].scatter(pos_track["Pos"].values, -neg_track["Count"].values, color=col_neg_track, zorder=4,
                         alpha=0.7, lw=0.1, ec="black", s=60, label=label_neg_track)
        axes[ax].set_yticklabels(abs(axes[ax].get_yticks()))

    if ymargin is not None:
        axes[ax].set_ylim(-np.max(neg_track["Count"])-ymargin, np.max(pos_track["Count"])+ymargin)


def add_ax_text(score, pvalue, y_text, x_text, ax, y_shift=0.2, y_adjust=0, equal_less=True):

    ax.text(x_text, y_text+(y_text*y_shift)+y_adjust,
            fr'$\mathit{{Score}}$ = {np.round(score, 2)}', ha='center', va='center', fontsize=13.5, color="black")

    equal = "≤" if equal_less else "="
    ax.text(x_text, y_text-(y_text*y_shift)+y_adjust,
            fr'$\mathit{{p}}$-value {equal} {pvalue}', ha='center', va='center', fontsize=13.5, color="black")


def get_transcript_ids(gene, maf_df_f, o3d_seq_df):

    canonical_tr = maf_df_f[maf_df_f["SYMBOL"] == gene].canonical_Feature.unique()[0]
    o3d_tr = o3d_seq_df[o3d_seq_df["Gene"] == gene].Ens_Transcr_ID.values[0]

    return canonical_tr, o3d_tr


def plot_gene_selection_signals(maf_trunc_count_gene,
                                maf_miss_count_gene,
                                maf_synon_count_gene,
                                frameshift_indels_count_gene,
                                inframe_indels_count_gene,
                                ofml_muts_score_gene,
                                o3d_score_gene,
                                o3d_prob_gene,
                                plot_pars,
                                domain_df=None,
                                add_track_title_text=True,
                                add_score_text=False,
                                rm_spines=False,
                                light_spines=False,
                                output_path=None):


    gene_len = len(o3d_prob_gene)
    fig, axes = plt.subplots(len(plot_pars["h_ratios"]), 1,
                             figsize=plot_pars["fsize"],
                             sharex=True,
                             gridspec_kw={'hspace': plot_pars["hspace"],
                                          'height_ratios': plot_pars["h_ratios"].values()})
    colors = plot_pars["colors"]


    # Omega
    # =====

    # Omega trunc
    ax=0
    plot_count_track(maf_trunc_count_gene, axes, ax=ax,  gene_len=gene_len, col_pos_track=colors["omega_trunc"])
    if add_score_text:
        add_ax_text(score=np.round(omega_truncating_score, 2),
                    pvalue=omega_truncating_pvalue,
                    y_text=np.max(maf_synon_count_gene["Count"])/2,
                    x_text=gene_len*plot_pars["score_txt_x_coord"],
                    ax=axes[ax],
                    y_shift=0.9,
                    y_adjust=2.4)
    if add_track_title_text:
        axes[ax].text(gene_len*plot_pars["track_title_x_coord"], 4.7,
                      fr'$\mathbf{{Omega}}$ $\mathbf{{Truncating}}$',
                      ha='center', va='center', fontsize=plot_pars["track_title_fontsize"], color="black")
    axes[ax].set_ylabel('Truncating\ncount', fontsize=plot_pars["ylabel_fontsize"], rotation=45, labelpad=plot_pars["ylabel_pad"], va='center')

    n_max = np.max(maf_trunc_count_gene["Count"])
    j = 12
    i = 24
    axes[ax].set_ylim(-n_max/i, n_max + n_max/j)

    axes[ax].yaxis.set_major_locator(MaxNLocator(integer=True, nbins=3))
    axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])
    if rm_spines:
        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)
    elif light_spines:
        axes[ax].spines['top'].set_color(colors["hv_lines"])
        axes[ax].spines['right'].set_color(colors["hv_lines"])

    # Omega synonym
    ax=1
    plot_count_track(maf_synon_count_gene, axes, ax=ax, gene_len=gene_len, col_pos_track=colors["omega_synon"], hv_lines=colors["hv_lines"])
    axes[ax].set_ylabel('Synonymous\ncount', fontsize=plot_pars["ylabel_fontsize"], rotation=45, labelpad=plot_pars["ylabel_pad"], va='center')

    n_max = np.max(maf_synon_count_gene["Count"])
    j = 10
    i = 85.714
    axes[ax].set_ylim(-n_max/i, n_max + n_max/j)

    axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])
    axes[ax].yaxis.set_major_locator(MaxNLocator(integer=True, nbins=3))
    axes[ax].spines['top'].set_visible(False)
    if rm_spines:
        axes[ax].spines['right'].set_visible(False)
    elif light_spines:
        axes[ax].spines['right'].set_color(colors["hv_lines"])

    # Omega miss
    ax=2
    plot_count_track(maf_miss_count_gene, axes, ax=ax, gene_len=gene_len, col_pos_track=colors["omega_miss"], hv_lines=colors["hv_lines"])
    y_text=np.max(maf_miss_count_gene["Count"])/2
    if add_score_text:
        add_ax_text(score=np.round(omega_misss_score, 2),
                    pvalue=omega_misss_pvalue,
                    y_text=y_text,
                    x_text=gene_len*plot_pars["score_txt_x_coord"],
                    y_shift=0.3,
                    ax=axes[ax],
                    y_adjust=4)
    if add_track_title_text:
        axes[ax].text(gene_len*plot_pars["track_title_x_coord"], y_text+y_text*0.6935,
                      fr'$\mathbf{{Omega}}$ $\mathbf{{Missense}}$',
                      ha='center', va='center', fontsize=plot_pars["track_title_fontsize"], color="black")
    axes[ax].set_ylabel('Missense\ncount', fontsize=plot_pars["ylabel_fontsize"], rotation=45, labelpad=plot_pars["ylabel_pad"], va='center')

    n_max = np.max(maf_miss_count_gene["Count"])
    j = 10.74
    i = 14.5
    axes[ax].set_ylim(-n_max/i, n_max + n_max/j)

    axes[ax].spines['top'].set_visible(False)
    axes[ax].yaxis.set_major_locator(MaxNLocator(integer=True, nbins=3))
    axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])
    if rm_spines:
        axes[ax].spines['right'].set_visible(False)
    elif light_spines:
        axes[ax].spines['right'].set_color(colors["hv_lines"])


    # Oncodrive3D
    # ===========

    ax=4
    axes[ax].plot(range(1, gene_len+1), o3d_score_gene["O3D_score_norm"], zorder=2, color=colors["o3d_score"], lw=1, label="Clustering score")
    axes[ax].plot(range(1, gene_len+1), -o3d_prob_gene, zorder=3, color=colors["o3d_prob"], lw=1, label="Missense mut probability")
    axes[ax].fill_between(o3d_score_gene['Pos'], 0, o3d_score_gene["O3D_score_norm"], where=(o3d_score_gene['Cluster'] == 1),
                          color=colors["o3d_cluster"], alpha=0.3, label='Cluster', zorder=0, lw=2)
    if add_score_text:
        add_ax_text(score=np.round(o3d_score, 2),
                    pvalue=o3d_pvalue,
                    y_text=(np.max(o3d_score_gene["O3D_score_norm"]) + np.max(o3d_prob_gene))/2,
                    x_text=gene_len*plot_pars["score_txt_x_coord"],
                    ax=axes[ax],
                    y_shift=0.15,
                    y_adjust=-0.0075,
                    equal_less=True)
    if add_track_title_text:
        axes[ax].text(gene_len*plot_pars["track_title_x_coord"], 0.015,
                      fr'$\mathbf{{3D}}$-$\mathbf{{Clustering}}$',
                      ha='center', va='center', fontsize=plot_pars["track_title_fontsize"], color="black")

    axes[ax].set_ylabel('Score &\nprobability', fontsize=plot_pars["ylabel_fontsize"], rotation=45, labelpad=plot_pars["ylabel_pad"], va='center')
    axes[ax].yaxis.set_major_locator(MaxNLocator(integer=False, nbins=4))
    tick_labels = [np.round(abs(tick), 3) for tick in axes[ax].get_yticks()]
    axes[ax].set_yticklabels(tick_labels)
    axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])

    axes[ax].legend(loc="upper left", fontsize=plot_pars["legend_fontsize"], frameon=plot_pars["legend_frameon"])
    #axes[ax].legend(bbox_to_anchor=[0, 0.35, 1, 0])
    if rm_spines:
        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)
    elif light_spines:
        axes[ax].spines['top'].set_color(colors["hv_lines"])
        axes[ax].spines['right'].set_color(colors["hv_lines"])


    # OncodriveFML
    # ============

    ax=6
    axes[ax].vlines(ofml_muts_score_gene["Pos"].values, ymin=0, ymax=ofml_muts_score_gene["Count"], color=colors["hv_lines"], lw=1, zorder=1, alpha=0.5)
    fml_scatter = axes[ax].scatter(ofml_muts_score_gene["Pos"].values, ofml_muts_score_gene["Count"].values,
                                   c=ofml_muts_score_gene["CADD_score"].values, cmap=colors["ofml"], zorder=4,
                                   alpha=0.7, lw=0.1, ec="black", s=60, label="SNV")

    y_text=np.round(np.max(ofml_muts_score_gene["Count"])/2)
    if add_score_text:
        add_ax_text(score=np.round(ofml_score, 2),
                    pvalue=ofml_pvalue,
                    y_text=y_text,
                    x_text=gene_len*plot_pars["score_txt_x_coord"],
                    ax=axes[ax],
                    y_adjust=3,
                    equal_less=True)

    if add_track_title_text:
        axes[ax].text(gene_len*plot_pars["track_title_x_coord"], y_text+y_text*0.8,
                      fr'$\mathbf{{Functional}}$ $\mathbf{{Impact}}$',
                      ha='center', va='center', fontsize=plot_pars["track_title_fontsize"], color="black")

    inset_ax = inset_axes(axes[ax], width="10%", height="10%", loc='center left',
                          bbox_to_anchor=plot_pars["ofml_cbar_coord"], bbox_transform=axes[ax].transAxes, borderpad=0)
    cbar = plt.colorbar(fml_scatter, cax=inset_ax, orientation='horizontal')
    cbar.set_label('Impact score', fontsize=plot_pars["legend_fontsize"])
    cbar.ax.tick_params(labelsize=plot_pars["ticksize"])

    n_max = np.max(ofml_muts_score_gene["Count"])
    j = 14.5
    i = 19.33
    axes[ax].set_ylim(-n_max/i, n_max + n_max/j)

    axes[ax].set_ylabel('SNV count', fontsize=plot_pars["ylabel_fontsize"], rotation=45, labelpad=plot_pars["ylabel_pad"], va='center')
    axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])

    if rm_spines:
        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)
    elif light_spines:
        axes[ax].spines['top'].set_color(colors["hv_lines"])
        axes[ax].spines['right'].set_color(colors["hv_lines"])


    # Frameshift enrichment
    # =====================

    ax=8
    plot_count_track(frameshift_indels_count_gene, axes, ax=ax, neg_track=inframe_indels_count_gene, gene_len=gene_len,
                     col_pos_track=colors["frameshift"], col_neg_track=colors["inframe"], label_pos_track="Frameshift", label_neg_track="Inframe",
                     ymargin=1, hv_lines=colors["hv_lines"])
    if add_score_text:
        add_ax_text(score=np.round(indels_score, 2),
                    pvalue=np.round(indels_pvalue, 2),
                    y_text=(np.max(frameshift_indels_count_gene["Count"]) + np.max(inframe_indels_count_gene["Count"]))/2,
                    x_text=gene_len*plot_pars["score_txt_x_coord"],
                    ax=axes[ax],
                    y_shift=0.24,
                    y_adjust=0.4,
                    equal_less=False)
    if add_track_title_text:
        axes[ax].text(gene_len*plot_pars["track_title_x_coord"], 3.9,
                      fr'$\mathbf{{Frameshift}}$ $\mathbf{{enrichment}}$',
                      ha='center', va='center', fontsize=plot_pars["track_title_fontsize"], color="black")
    axes[ax].legend(loc="upper left", fontsize=plot_pars["legend_fontsize"], frameon=plot_pars["legend_frameon"])
    axes[ax].set_ylabel('Indels count', fontsize=plot_pars["ylabel_fontsize"], rotation=45, labelpad=plot_pars["ylabel_pad"], va='center')

    axes[ax].yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    tick_labels = [abs(int(tick)) for tick in axes[ax].get_yticks()]
    axes[ax].set_yticklabels(tick_labels)

    if rm_spines:
        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)
    elif light_spines:
        axes[ax].spines['top'].set_color(colors["hv_lines"])
        axes[ax].spines['right'].set_color(colors["hv_lines"])


    # Domain
    # ======

    if isinstance(domain_df, pd.DataFrame):
        ax=10
        domain_color_dict = {}

        for n, name in enumerate(domain_df["Description"].unique()):
            domain_color_dict[name] = f"C{n}"

        n = 0
        added_domain = []
        for i, row in domain_df.iterrows():
            if pd.Series([row["Description"], row["Begin"], row["End"]]).isnull().any():
                continue

            name = row["Description"]
            start = int(row["Begin"])
            end = int(row["End"])
            axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=domain_color_dict[name])
            if name not in added_domain:
                y = -0.04
                axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                added_domain.append(name)
        axes[ax].set_yticks([])
    else:
        axes[10].remove()


    # Spaces
    # =====
    axes[3].remove()
    axes[5].remove()
    axes[7].remove()
    axes[9].remove()


    # X axes
    # ======
    axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])
    axes[ax].tick_params(axis='x', labelsize=plot_pars["ticksize"])
    axes[ax].set_xlabel('Protein position', fontsize=plot_pars["xlabel_fontsize"])

    if output_path is not None:
        print(f"Saving {output_path}")
        plt.savefig(output_path, dpi=plot_pars["dpi"], bbox_inches='tight')
    plt.show()


def generate_plot(gene,
                  o3d_seq_df,
                  o3d_annot_df,
                  o3d_prob,
                  o3d_score,
                  maf_df_f,
                  snv_df,
                  trunc_df,
                  synon_df,
                  miss_df,
                  plot_pars,
                  domain_df=None,
                  add_track_title_text=False,
                  add_score_text=False,
                  rm_spines=True,
                  light_spines=False,
                  output_path=None):

    # Prepare data
    # ============

    # O3D gene
    o3d_score_df_gene, o3d_score_gene, o3d_prob_gene, o3d_top_pos_score_gene, uni_id, gene_pos = get_o3d_gene_data(gene,
                                                                                                                   o3d_seq_df,
                                                                                                                   o3d_prob,
                                                                                                                   o3d_score)

    # MAF gene
    snv_count_gene = get_maf_pos_count(gene, snv_df, gene_pos)
    trunc_count_gene = get_maf_pos_count(gene, trunc_df, gene_pos)
    synon_count_gene = get_maf_pos_count(gene, synon_df, gene_pos)
    miss_count_gene = get_maf_pos_count(gene, miss_df, gene_pos)

    # OFML gene
    ofml_score_gene = get_ofml_score_gene(gene, fig3_data)

    # Indels
    frameshift_indels_df, inframe_indels_df = get_frameshift_indels_maf(maf_df_f)
    frameshift_indels_count_gene, inframe_indels_count_gene = get_frameshift_indels_gene(gene,
                                                                                         frameshift_indels_df,
                                                                                         inframe_indels_df,
                                                                                         gene_pos)

    # Domain
    domain_gene = o3d_annot_df[(o3d_annot_df["Gene"] == gene) &
                               (o3d_annot_df["Type"] == "DOMAIN") &
                               (o3d_annot_df["Evidence"] == "Pfam")].reset_index(drop=True)

    # Transcripts and Uniprot ID
    canonical_tr, o3d_tr = get_transcript_ids(gene, maf_df_f, o3d_seq_df)
    print(f"> {gene} - {canonical_tr} - {o3d_tr} - {uni_id}")


    # Plot
    # ====

    if output_path is not None:
        output_path = f"{output_path}/{gene}.png"
    plot_gene_selection_signals(trunc_count_gene,
                                miss_count_gene,
                                synon_count_gene,
                                frameshift_indels_count_gene,
                                inframe_indels_count_gene,
                                ofml_score_gene,
                                o3d_score_df_gene,
                                o3d_prob_gene,
                                plot_pars,
                                domain_gene,
                                add_track_title_text,
                                add_score_text,
                                rm_spines,
                                light_spines,
                                output_path)





# Prepare data
# ============

def get_o3d_gene_data(gene,
                      seq_df,
                      prob_dict,
                      score_df):

    # Subset gene
    seq_df_gene = seq_df[seq_df["Gene"] == gene]
    gene_len = len(seq_df_gene.Seq.values[0])
    gene_pos = pd.DataFrame({"Pos" : range(1, gene_len+1)})
    uni_id, af_f = seq_df_gene[["Uniprot_ID", "F"]].values[0]
    prob_gene = np.array(prob_dict[f"{uni_id}-F{af_f}"])
    score_gene_df = score_df[score_df["Gene"] == gene].reset_index(drop=True)
    score_gene, top_pos_score_gene = score_df[["Score", "Score_obs_sim"]].values[0]

    ## O3D score vector
    score_gene_df = gene_pos.merge(score_gene_df[["Pos", "Score_obs_sim", "C", "C_ext"]], how="left", on="Pos")

    # Don't include Extended clusters
    score_gene_df["C"] = (score_gene_df["C"] == 1) & (score_gene_df["C_ext"] == 0)
    score_gene_df["C"] = score_gene_df["C"].astype(int)
    score_gene_df = score_gene_df.drop(columns=["C_ext"])

    score_gene_df.columns = "Pos", "O3D_score", "Cluster"
    score_gene_df["O3D_score"] = score_gene_df["O3D_score"].fillna(0)
    score_gene_df["Cluster"] = score_gene_df["Cluster"].fillna(0)
    score_gene_df["O3D_score_norm"] = score_gene_df["O3D_score"] / sum(score_gene_df["O3D_score"])

    return score_gene_df, score_gene, prob_gene, top_pos_score_gene, uni_id, gene_pos


def preprocess_maf(maf_df):

    maf_df["CLEAN_SAMPLE_ID"] = maf_df["SAMPLE_ID"].apply(lambda x: "_".join(x.split("_")[1:3]))

    # Reduce the number of samples and the number of mutations
    samples_histo_findings = ['P19_0017_BDO_01', 'P19_0017_BTR_01',
                          'P19_0032_BDO_01', 'P19_0032_BTR_01',
                          'P19_0044_BDO_01', 'P19_0044_BTR_01']
    maf_df_f = maf_df.loc[(maf_df["VAF"] <= 0.35) &
                          # (maf_df["FILTER.repetitive_variant"] == False) & # filter not well defined yet; may hide hotspots
                          (~maf_df["FILTER.not_in_panel"]) &
                          (~maf_df["FILTER.no_pileup_support"]) & # avoid variants w/o VAF recomputed
                          (~maf_df["FILTER.n_rich"]) &
                          (~maf_df["FILTER.low_mappability"]) &
                          (~maf_df["FILTER.other_sample_SNP"]) &
                          (~maf_df["SAMPLE_ID"].isin(samples_histo_findings))
                         ].reset_index(drop = True)

    # SNV
    snvs_maf = maf_df_f[(maf_df_f["TYPE"] == 'SNV') &
                        (maf_df_f["canonical_SYMBOL"].isin(gene_order)) &
                        (maf_df_f["canonical_Protein_position"] != '-' )
                       ].reset_index(drop = True)
    snvs_maf["canonical_Protein_position"] = snvs_maf["canonical_Protein_position"].astype(int)

    # Omega
    maf_df_trunc = snvs_maf[snvs_maf["canonical_Consequence_broader"].isin(['nonsense', "essential_splice"])]
    maf_df_synon = snvs_maf[snvs_maf["canonical_Consequence_broader"] == 'synonymous']
    maf_df_miss = snvs_maf[snvs_maf["canonical_Consequence_broader"] == 'missense']

    return maf_df_f, snvs_maf, maf_df_trunc, maf_df_synon, maf_df_miss


def get_maf_pos_count(gene, maf_df, gene_pos):

    # Subset gene and cols
    cols = ["SYMBOL",
            "TYPE",
            "Feature",
            "canonical_Feature",
            "canonical_Consequence",
            "canonical_Consequence_broader",
            "canonical_Protein_position"]

    df_gene = maf_df.loc[maf_df["SYMBOL"] == gene, cols]

    # Get per-position mutations count
    df_gene_count = df_gene.groupby("canonical_Protein_position").apply(lambda x: len(x)).reset_index().astype(int)
    df_gene_count.columns = "Pos", "Count"
    df_gene_count = gene_pos.merge(df_gene_count, how="left", on="Pos")

    return df_gene_count


def get_ofml_score_gene(gene, fig3_data_path):

    mutations_in_gene = pd.read_table(f"{fig3_data_path}/mutations_scored.{gene}.tsv")
    ofml_muts_score_gene = mutations_in_gene.groupby(by = "canonical_Protein_position").agg( { "MUT_ID" : 'count',
                                                                                               "CADDscore" : 'mean'}).reset_index()
    ofml_muts_score_gene.columns = "Pos", "Count", "CADD_score"

    return ofml_muts_score_gene


def get_frameshift_indels_maf(maf_df_f):

    # Filter and somatic only
    indel_maf_df = maf_df_f.loc[
                        (maf_df_f["TYPE"].isin(["INSERTION", "DELETION"]))
                        ].reset_index(drop = True)
    indel_maf_df["INDEL_LENGTH"] = (indel_maf_df["REF"].str.len() - indel_maf_df["ALT"].str.len()).abs()
    indel_maf_df["INDEL_INFRAME"] = [ x % 3 == 0 for x in indel_maf_df["INDEL_LENGTH"] ]
    indel_maf_df.loc[indel_maf_df["INDEL_LENGTH"] >= 15 , "INDEL_INFRAME"] = False
    indel_maf_df.loc[indel_maf_df["canonical_Consequence_broader"] == 'nonsense' , "INDEL_INFRAME"] = False
    indel_maf_df.loc[indel_maf_df["canonical_Consequence_broader"] == 'essential_splice' , "INDEL_INFRAME"] = False

    # Observed frameshift indels + inframes of length >= 5 AA
    frameshift_indels = indel_maf_df[(~indel_maf_df["INDEL_INFRAME"]) &
                                 (indel_maf_df["canonical_Protein_position"] != '-' )].reset_index(drop = True)

    # Inframe indels
    inframe_indels = indel_maf_df[(indel_maf_df["INDEL_INFRAME"]) &
                                  (indel_maf_df["canonical_Protein_position"] != '-' )].reset_index(drop = True)

    return frameshift_indels, inframe_indels


def get_frameshift_indels_gene(gene, frameshift_indels_df, inframe_indels_df, gene_pos):

    ## Frameshift indels
    frameshift_indels_gene = frameshift_indels_df[frameshift_indels_df["canonical_SYMBOL"] == gene]

    # Count by pos considering first pos as protein pos
    frameshift_indels_gene.canonical_Protein_position = frameshift_indels_gene.canonical_Protein_position.apply(lambda x: x.split("-")[0])
    frameshift_indels_gene.canonical_Protein_position = frameshift_indels_gene.canonical_Protein_position[
        frameshift_indels_gene.canonical_Protein_position.apply(lambda x: x.isdigit())]
    frameshift_indels_count_gene = frameshift_indels_gene.groupby("canonical_Protein_position").apply(lambda x: len(x)).reset_index().astype(int)
    frameshift_indels_count_gene.columns = "Pos", "Count"
    frameshift_indels_count_gene = gene_pos.merge(frameshift_indels_count_gene, how="left", on="Pos")

    ## Inframe indels
    inframe_indels_gene = inframe_indels_df[inframe_indels_df["canonical_SYMBOL"] == gene]

    # Count
    inframe_indels_gene.canonical_Protein_position = inframe_indels_gene.canonical_Protein_position.apply(lambda x: x.split("-")[0])
    inframe_indels_gene.canonical_Protein_position = inframe_indels_gene.canonical_Protein_position[
        inframe_indels_gene.canonical_Protein_position.apply(lambda x: x.isdigit())]
    inframe_indels_count_gene = inframe_indels_gene.groupby("canonical_Protein_position").apply(lambda x: len(x)).reset_index().astype(int)
    inframe_indels_count_gene.columns = "Pos", "Count"
    inframe_indels_count_gene = gene_pos.merge(inframe_indels_count_gene, how="left", on="Pos")

    return frameshift_indels_count_gene, inframe_indels_count_gene

seq_regions_df = pd.read_table("/workspace/nobackup/scratch/oncodrive3d/datasets_240506/seq_for_mut_prob.tsv")
tb = tabix.open("/workspace/datasets/CADD/v1.6/hg38/whole_genome_SNVs.tsv.gz")


def get_scores_gene(gene_info):
    chromosome = gene_info["Chr"].values[0]
    reverse_strand = gene_info["Reverse_strand"].values[0] == 1
    exon_coords = eval(gene_info["Exons_coord"].values[0])
    if reverse_strand:
        exon_coords = [(y, x) for x, y in exon_coords]


    position_mutation = dict()
    for start, stop in exon_coords:
        for e in tb.query(str(chromosome), start, stop):
            if int(e[1]) not in position_mutation:
                position_mutation[int(e[1])] = dict()
            position_mutation[int(e[1])][e[3]] = float(e[5])

    start2end_coords = { x:y for x, y in enumerate(sorted(position_mutation.keys(), reverse = not reverse_strand))}

    return position_mutation, start2end_coords




def score_mutations_from_gene(snv_data):
    scores_muts = []
    for ind, row in snv_data.iterrows():
        try:
            scores_muts.append(scores_all_mutations_in_gene[row["POS"]].get(row["ALT"], 0))
            print("found")
        except:
            scores_muts.append(0)
            print(row["Consequence"], "not found")
    snv_data["CADDscore"] = scores_muts
    return snv_data



# for gene in ["TP53"]:
for gene in gene_order:
    gene_info = seq_regions_df[seq_regions_df["Gene"] == gene]
    scores_all_mutations_in_gene, scores_of_protein_position = get_scores_gene(gene_info)

    mutations_in_gene = snvs_maf[snvs_maf["canonical_SYMBOL"] == gene].copy()
    scored_mutations_in_gene = score_mutations_from_gene(mutations_in_gene)

    scored_mutations_in_gene.to_csv(f"{data_dir}/mutations_scored.{gene}.tsv", sep = '\t', header= True, index = False)


