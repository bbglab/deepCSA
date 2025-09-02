#!/usr/bin/env python

import os
import time
import requests
import json
import warnings
import click
import pandas as pd
import seaborn as sns
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as mpatches



#####
# Define functions
#####


def get_normal_maf(path_maf, gene_list = None, only_protein_pos=True, truncating=True):

    # TODO: check for correct filtering
    maf_df = pd.read_table(path_maf)

    if gene_list is None:
        gene_list = maf_df["canonical_SYMBOL"].unique().tolist()

    # Final filter
    maf_df_f = maf_df[
        (maf_df["TYPE"].isin(["SNV", "INSERTION", "DELETION"])
        ) & (maf_df["canonical_SYMBOL"].isin(gene_list))
        ].reset_index(drop = True)

    if only_protein_pos:
        maf_df_f = maf_df_f[maf_df_f["canonical_Protein_position"] != '-' ]

    maf_df_f["plotting_broad_consequence"] = maf_df_f["canonical_Consequence_broader"]
    maf_df_f.loc[(maf_df_f["TYPE"].isin(["INSERTION", "DELETION"])), "plotting_broad_consequence"] = "indel"
    maf_df_f["plotting_broad_consequence"] = maf_df_f["plotting_broad_consequence"].replace("splice_region_variant", "splicing")

    if truncating:
        maf_df_f["plotting_broad_consequence"] = maf_df_f["plotting_broad_consequence"].replace(
            {"nonsense": "truncating", "essential_splice": "truncating"}
            )

    maf_df_f['Alt_amino_acid'] = np.where(
        maf_df_f['plotting_broad_consequence'] == 'missense',
        maf_df_f['canonical_Amino_acids'].str.split('/').str[-1],
        np.nan
    )

    # Parse
    cols = ["canonical_SYMBOL",
            "canonical_Feature",
            "canonical_Protein_position",
            "plotting_broad_consequence",
            "Alt_amino_acid",
            "CHROM",
            "POS",
            "DEPTH",
            "ALT_DEPTH"]
    maf_df_f = maf_df_f[cols].rename(columns={
        "canonical_SYMBOL" : "Gene",
        "canonical_Feature" : "Ens_transcript_ID",
        "canonical_Protein_position" : "Pos",
        "plotting_broad_consequence" : "Consequence",
        "POS" : "DNA_POS"}
        )

    # Extract first pos if multiple ones are defined
    maf_df_f['Pos'] = maf_df_f['Pos'].str.extract(r'^(\d+)')
    maf_df_f['Pos'] = pd.to_numeric(maf_df_f['Pos'], errors='coerce')

    if only_protein_pos:
        maf_df_f = maf_df_f.dropna(subset="Pos").reset_index(drop=True)
        maf_df_f.Pos= maf_df_f.Pos.astype(int)

    return maf_df_f, gene_list



# Load cancer clustering
# ----------------------

def get_o3d_gene_data(
    gene,
    seq_df,
    o3d_pos_df
    ):

    # Subset gene
    seq_df_gene = seq_df[seq_df["Gene"] == gene]
    protein_len = len(seq_df_gene.Seq.values[0])
    gene_pos = pd.DataFrame({"Pos" : range(1, protein_len+1)})
    uni_id, af_f = seq_df_gene[["Uniprot_ID", "F"]].values[0]
    score_gene_df = o3d_pos_df[o3d_pos_df["Gene"] == gene].reset_index(drop=True)

    ## O3D score vector
    score_gene_df = gene_pos.merge(score_gene_df[["Pos", "Score_obs_sim", "C", "C_ext"]], how="left", on="Pos")

    # Don't include Extended clusters
    score_gene_df["C"] = (score_gene_df["C"] == 1) & (score_gene_df["C_ext"] == 0)
    score_gene_df["C"] = score_gene_df["C"].astype(int)
    score_gene_df = score_gene_df.drop(columns=["C_ext"])

    score_gene_df.columns = ["Pos", "O3D_score", "Cluster"]
    score_gene_df["O3D_score"] = score_gene_df["O3D_score"].fillna(0)
    score_gene_df["Cluster"] = score_gene_df["Cluster"].fillna(0)

    return score_gene_df



warnings.filterwarnings("ignore", message="Warning: converting a masked element to nan.", category=UserWarning)
warnings.filterwarnings("ignore", message="FixedFormatter should only be used together with FixedLocator", category=UserWarning)


plot_pars = {
    "fsize"                     : (12.3, 8.3),
    "hspace"                    : 0.1,                    # General space between all tracks
    "track_title_x_coord"       : 0.83,                   # x-coordinate (respect to protein len) for track txt title
    "track_title_fontsize"      : 14,
    "ylabel_fontsize"           : 13.5,
    "xlabel_fontsize"           : 13.5,
    "ylabel_pad"                : 38,
    "ticksize"                  : 10.25,
    "legend_fontsize"           : 12,
    "legend_frameon"            : False,
    "j_margin"                  : 15,
    "i_margin"                  : 45,
    "c_margin"                  : 0,
    "y_labels_coord"            : [-0.05, 0.5],
    "cnsq_bbox_to_anchor"       : [1.1123, 0.84],
    "o3d_bbox_to_anchor"        : [1.27, 1.35],
    "saturation_bbox_to_anchor" : [1.086, 8.5, 1.4, 5.05],
    "depth_bbox_to_anchor"      : [1.086, 6.1, 1.4, 5.05],
    "depth_res_bbox_to_anchor"  : [1.086, -2, 1.4, 5.05],
    "depth_0_bbox_to_anchor"    : [1.246, -2], # 3.9 5.9
    "sse_bbox_to_anchor"        : [1.2645, 2.4],
    "domain_sel_bbox_to_anchor" : [1.086, -3, 1.4, 4.8],
    "domain_y_bbox_to_anchor"   : 5,
    "domain_x_bbox_to_anchor"   : {},
    "legend_depth_fontsize"     : 10.5,
    "txt_fontsize"              : 9,
    "len_txt_thr"               : 2400,
    "dpi"                       : 300,
    "colors"                    : {
                                "truncating"  : "#FB8E6F",  # group nonsense + splicing
                                "nonsense"    : "#FB8E6F",
                                "synonymous"  : "#ACECBD",
                                "missense"    : "#FBD180",
                                "o3d_score"   : "#6DBDCC",
                                "o3d_cluster" : "#DAF0F8", # "skyblue",
                                "indel"       : "#ECC4F7",
                                "splicing"    : "#A1C5DF",
                                "vlines"      : "lightgray",
                                "plddt_pacc"  : "#f7f3ff",
                                "site_selection" : "#b3b3ff",
                                "site_selection_hits" : "#4d4dff",               # #4d4dff #1a1aff
                                "exon_selection" : {
                                                    "missense" : "#FBD180",
                                                    "truncating" : "#FB8E6F"
                                                    },
                                "depth_below_thr" : "#d9d9d9"
                                },
    "sse_colors"                : {
                                'Coil'   : "#D5E8D4",
                                'Helix'  : "#F7CAC9",
                                'Ladder' : "#A7C7E7"
                                },
    "sse_lw"                    : 1
}


def get_exon_depth_saturation(gene_depth, gene_mut, dna=False):

    # Exon average depth
    gene_depth = gene_depth.copy()
    exon_depth = gene_depth.groupby("EXON_RANK").apply(
        lambda x: 0 if sum(x.COVERED) == 0 else sum(x.DEPTH) / sum(x.COVERED)).reset_index().rename(columns = {0 : "DEPTH"})
    exon_depth["START_PROT_POS"] = gene_depth.groupby("EXON_RANK").apply(lambda x: x.PROT_POS.min())

    # Exon saturation
    if dna:
        exon_depth["COVERED"] = gene_depth.groupby("EXON_RANK").apply(lambda x: sum(x.COVERED)).values
        gene_depth["MUTATED"] = gene_depth.DNA_POS.isin(gene_mut.DNA_POS.unique()).astype(int)
        exon_depth["MUTATED"] = gene_depth.groupby("EXON_RANK").apply(lambda x: sum(x.MUTATED)).values
        exon_depth["SATURATION"] = exon_depth["MUTATED"] / exon_depth["COVERED"]
    else:

        # If an amino acids from the same codon are in two exons, consider the position only in the second one
        # (but look at booth to know if the residue is in the panel and if it is mutated)
        exon_depth_prot = gene_depth.groupby("PROT_POS").apply(lambda x: (x.COVERED.max(), x.EXON_RANK.max())).reset_index()
        exon_depth_prot[["COVERED", "EXON_RANK"]] = pd.DataFrame(exon_depth_prot[0].tolist(), index=exon_depth_prot.index)
        exon_depth_prot = exon_depth_prot.drop(columns=[0])

        exon_depth["COVERED"] = exon_depth_prot.groupby("EXON_RANK").apply(lambda x: sum(x.COVERED)).values
        exon_depth_prot["MUTATED"] = exon_depth_prot["PROT_POS"].isin(gene_mut["Pos"].unique()).astype(int)
        exon_depth["MUTATED"] = exon_depth_prot.groupby("EXON_RANK").apply(lambda x: sum(x.MUTATED)).values
        exon_depth["SATURATION"] = exon_depth["MUTATED"] / exon_depth["COVERED"]

    # check_mutated_masked = exon_depth_prot[(exon_depth_prot["MUTATED"] == 1) & (exon_depth_prot["COVERED"] == 0)]
    # check_mutated_masked["GENE"] = gene_depth.GENE.unique()[0]

    return exon_depth


def get_exon_mid_prot_pos(exon_info, prot_len):

    lst_end_pos = []
    for i in range(len(exon_info)):
        exon = exon_info.iloc[i]
        start_pos = int(exon.START_PROT_POS)
        end_pos = int(exon_info.iloc[i+1].START_PROT_POS) if i < len(exon_info) -1 else prot_len
        lst_end_pos.append(end_pos)
    exon_info["END_PROT_POS"] = lst_end_pos
    exon_info["MID_PROT_POS"] = (exon_info["START_PROT_POS"] + exon_info["END_PROT_POS"]) / 2

    return exon_info


def get_res_coverage(dna_prot_df):

    res_depth = dna_prot_df.dropna(subset=["PROT_POS"])[["PROT_POS", "COVERED", "DEPTH"]].reset_index(drop=True)
    res_depth = res_depth.groupby("PROT_POS").mean().reset_index()

    return res_depth


def add_consecutive_numbers(nums, max_n):

    result = []
    for i in range(len(nums)):
        result.append(nums[i])
        # Check if the current number is the start of a consecutive sequence
        if i < len(nums) - 1 and nums[i] + 1 != nums[i + 1]:
            result.append(nums[i] + 1)

    # Add the last consecutive number after the final element
    result.append(nums[-1] + 1)

    return result


def where_plus(condition):
    """
    Util function to extend the color of mpl filling to the next position.
    """

    ix = np.where(condition)[0]

    if len(ix) > 0:
        ix = add_consecutive_numbers(ix, max_n=len(condition))
        if len(condition) in ix:
            ix.remove(len(condition))
        boolean_vector = np.zeros(len(condition), dtype=bool)
        boolean_vector[ix] = True

        return pd.Series(boolean_vector)

    else:
        return condition


def plot_count_track(
    count_df,
    protein_len,
    axes,
    colors_dict,
    ax=0,
    negative=False,
    label_pos_track=None,
    label_neg_track=None,
    ymargin=None,
    alpha=1,
    indel=False,
    n_batches = 10
    ):

    # Shuffle the data and split into batches
    shuffled_df = count_df.sample(frac=1, random_state=42).reset_index(drop=True)
    batches = np.array_split(shuffled_df, n_batches)

    legend_list = []
    pos_df = pd.DataFrame({"Pos" : range(1, protein_len+1)})

    for batch_idx, batch in enumerate(batches):
        for cnsq in ['indel', 'truncating', 'nonsense', 'missense', 'synonymous', 'splicing']:

            if (not indel and cnsq == "indel") or (cnsq not in batch["Consequence"].unique()) :
                continue

            count_cnsq_df = batch[batch["Consequence"] == cnsq].reset_index(drop=True)
            count_cnsq_df = pos_df.merge(count_cnsq_df, on="Pos", how="left")

            if negative:
                axes[ax].vlines(count_cnsq_df["Pos"], ymin=-count_cnsq_df["Count"], ymax=0, lw=1, zorder=1, alpha=0.5, color=colors_dict["vlines"])
                axes[ax].scatter(count_cnsq_df["Pos"], -count_cnsq_df["Count"], color='white', zorder=3, lw=1, ec="white") # To cover the overlapping needle top part
                if cnsq not in legend_list:
                    axes[ax].scatter(
                        count_cnsq_df["Pos"].values, -count_cnsq_df["Count"].values, zorder=4,
                        alpha=alpha, lw=0.1, ec="black", s=60, label=cnsq.capitalize(), color=colors_dict[cnsq]
                        )
                    legend_list.append(cnsq)

                else:
                    axes[ax].scatter(
                        count_cnsq_df["Pos"].values, -count_cnsq_df["Count"].values, zorder=4,
                        alpha=alpha, lw=0.1, ec="black", s=60, color=colors_dict[cnsq]
                        )

            else:
                axes[ax].vlines(count_cnsq_df["Pos"], ymin=0, ymax=count_cnsq_df["Count"], lw=1, zorder=1, alpha=0.5, color=colors_dict["vlines"])
                axes[ax].scatter(count_cnsq_df["Pos"], count_cnsq_df["Count"], color='white', zorder=3, lw=1, ec="white") # To cover the overlapping needle top part
                if cnsq not in legend_list:
                    axes[ax].scatter(
                        count_cnsq_df["Pos"].values, count_cnsq_df["Count"].values, zorder=4,
                        alpha=alpha, lw=0.1, ec="black", s=60, label=cnsq.capitalize(), color=colors_dict[cnsq]
                        )
                    legend_list.append(cnsq)
                else:
                    axes[ax].scatter(
                        count_cnsq_df["Pos"].values, count_cnsq_df["Count"].values, zorder=4,
                        alpha=alpha, lw=0.1, ec="black", s=60, color=colors_dict[cnsq]
                        )


def get_domain_selection_gene(domain_selection_all_genes, domain_gene, gene, sort_by_impact=True):

    domain_selection_gene = domain_selection_all_genes.merge(domain_gene.rename(columns={"NAME": "gene"})[["gene", "DOMAIN_ID", "Begin", "End"]], how="inner")
    domain_selection_gene = domain_selection_gene.sort_values(["Begin", "impact"]).reset_index(drop=True)
    if sort_by_impact:
        average_dnds = domain_selection_gene.groupby("DOMAIN_ID").apply(lambda x: x["dnds"].mean()).to_dict()
        domain_selection_gene["average_dnds"] = domain_selection_gene["DOMAIN_ID"].map(average_dnds)
        domain_selection_gene = domain_selection_gene.sort_values(["average_dnds", "impact"], ascending=[False, True]).reset_index(drop=True)
    domain_selection_gene["x_pos"] = domain_selection_gene["DOMAIN_ID"].astype("category").cat.set_categories(domain_selection_gene["DOMAIN_ID"].unique()).cat.codes
    domain_selection_gene["x_pos"] += np.where(domain_selection_gene["impact"] == "missense", -0.15, 0.15)

    return domain_selection_gene


def color_bar(
            axes,
            ax,
            cmap,
            norm,
            general_ticks=False,
            bbox_to_anchor=[1.086, -2, 1.4, 5.05],
            label='Cbar label',
            fontsize=12,
            tick_size=10.25,
            tick_span=0,
            borderpad=10
            ):

            # Color bar
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            inset_ax = inset_axes(
                axes[ax], width="10%", height="10%", loc='center left',
                bbox_to_anchor=bbox_to_anchor, bbox_transform=axes[ax].transAxes, borderpad=borderpad
                )
            cbar = plt.colorbar(sm, cax=inset_ax, orientation='horizontal')
            cbar.set_label(label, fontsize=fontsize, labelpad=7)
            cbar.ax.xaxis.set_label_position('top')

            cbar.ax.tick_params(labelsize=tick_size)
            span = norm.vmax - norm.vmin
            left_tick = norm.vmin + tick_span * span
            mid_tick = norm.vmin + span / 2
            right_tick = norm.vmax - tick_span * span
            tick_values = [left_tick, mid_tick, right_tick]
            cbar.set_ticks(tick_values)
            if general_ticks:
                cbar.set_ticklabels([tick for tick in ["Low", "High"]])


def get_equivalent_alt_prot_id(gene, seq_df, alt_seq_df):

    seq = seq_df[seq_df["Gene"] == gene].Seq.unique()[0]
    alt_seq = alt_seq_df[alt_seq_df["Gene"] == gene].Seq.unique()[0]

    if seq == alt_seq:
        return alt_seq_df[alt_seq_df["Gene"] == gene].Uniprot_ID.unique()[0]


def avg_per_pos_ddg(mut_df, ddg_prot):
    """
    Compute per-position average stability change upon mutations (DDG).
    """

    mut_df = mut_df.copy()
    mut_df = mut_df[mut_df["Consequence"] == "missense"].reset_index(drop=True)

    ddg_vec = np.repeat(0., len(ddg_prot.keys()))
    for pos, df in mut_df.groupby('Pos'):
        pos = str(pos)
        obs_mut = df.Alt_amino_acid
        if pos in ddg_prot:
            ddg_pos = ddg_prot[pos]
            ddg_pos = np.mean([ddg_pos[mut] for mut in obs_mut])
            ddg_vec[int(pos)-1] = ddg_pos

    return pd.DataFrame({"Pos": ddg_prot.keys(), "Stability_change": ddg_vec}).astype({"Pos": int, "Stability_change": float})


def get_ddg_gene(gene, gene_mut, o3d_annotations, o3d_seq_df,
                    # o3d_alt_seq_df
                    ):

    uni_id = o3d_seq_df[o3d_seq_df["Gene"] == gene].Uniprot_ID.unique()[0]
    ddg_path = f"{o3d_annotations}/stability_change/{uni_id}_ddg.json"
    # alt_uni_id = get_equivalent_alt_prot_id(gene, o3d_seq_df, o3d_alt_seq_df)
    # alt_ddg_path = os.path.join(o3d_annotations, "stability_change", f"{alt_uni_id}_ddg.json")
    if os.path.isfile(ddg_path):
        ddg_gene = json.load(open(ddg_path))
        ddg_gene = avg_per_pos_ddg(gene_mut, ddg_gene)
    # elif isinstance(alt_uni_id, str) and os.path.isfile(alt_ddg_path):
    #     ddg_gene = json.load(open(alt_ddg_path))
    #     ddg_gene = avg_per_pos_ddg(gene_mut, ddg_gene)
    else:
        ddg_gene = None

    return ddg_gene


def plot_gene_selection(mut_count_df,
                        mut_df,
                        o3d_df,
                        pdb_tool_df,
                        domain_df,
                        coverage_df,
                        site_selection_df,
                        exon_selection_df,
                        domain_selection_df,
                        gene,
                        max_depth,
                        protein_len,
                        plot_pars,
                        title,
                        ddg_df=None,
                        thr_selection=1e-5,
                        lst_tracks=["Mut_count", "Site_selection", "Res_depth", "Domain"],
                        default_track_order=False,
                        save=False,
                        filename="gene_selection.png"):

    h_ratios_dict = {
        "Mut_count": 0.20,
        "3d_clustering": 0.15,
        "Site_selection": 0.15,
        "Exon_selection": 0.15,
        "Exon_saturation": 0.03,
        "Exon_depth": 0.03,
        "Res_depth": 0.03,
        "Solvent_accessibility": 0.05,
        "Stability_change": 0.05,
        "Secondary_structure": 0.03,
        "Domain": 0.03,
        "Domain_selection": 0.035,
    }

    lst_tracks = [track.capitalize() for track in lst_tracks]
    if default_track_order:
        lst_tracks = [track for track in h_ratios_dict.keys() if track in lst_tracks]
    h_ratios = [h_ratios_dict[track] for track in lst_tracks]

    fig, axes = plt.subplots(
        len(h_ratios), 1,
        figsize=plot_pars["fsize"],
        sharex=True,
        gridspec_kw={'hspace': plot_pars["hspace"],
                    'height_ratios': h_ratios}
        )

    n_max = np.max(mut_count_df["Count"])
    n_3d_max = np.max(o3d_df["O3D_score"])


    # Track mut
    # ---------
    if "Mut_count" in lst_tracks:
        ax = lst_tracks.index("Mut_count")

        plot_count_track(mut_count_df, protein_len=protein_len, axes=axes, ax=ax, colors_dict=plot_pars["colors"])
        axes[ax].set_ylabel('Mutations', rotation=0, va='center', ha='right', fontsize=plot_pars["ylabel_fontsize"])
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])

        #axes[ax].set_ylim(0 - n_max/plot_pars["i_margin"], n_max + n_max/plot_pars["j_margin"])
        axes[ax].yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
        axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])

        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)


    # Track site selection
    # --------------------
    if "Site_selection" in lst_tracks:
        ax = lst_tracks.index("Site_selection")

        site_selection_df = pd.DataFrame({"Protein_position" : np.arange(protein_len)+1}).merge(site_selection_df, how="left")
        axes[ax].plot(
            site_selection_df.Protein_position, site_selection_df.Selection, zorder=1,
            color=plot_pars["colors"]["site_selection"], lw=1
            )
        axes[ax].fill_between(
            site_selection_df.Protein_position, 0, site_selection_df.Selection,
            color=plot_pars["colors"]["site_selection"], alpha=0.4, zorder=0, lw=1.5
            )

        site_selection_hits_df = pd.DataFrame({"Protein_position" : np.arange(protein_len)+1}).merge(
            site_selection_df[site_selection_df["p_value"] < thr_selection].reset_index(drop=True), how="left")
        axes[ax].fill_between(
            site_selection_hits_df.Protein_position, 0, site_selection_hits_df.Selection,
            color=plot_pars["colors"]["site_selection_hits"], alpha=1, zorder=2, lw=1.5, label="Significant"
            )

        n_max = np.max(site_selection_df.Selection)
        ax_ylim_min = 0 #if plot_pars["c_margin"] == 0 else 0 - n_max/c_margin
        ax_ylim_max = n_max + n_max/plot_pars["j_margin"]
        axes[ax].set_ylim(ax_ylim_min, ax_ylim_max)

        axes[ax].set_ylabel('Site selection', rotation=0, va='center', ha='right', fontsize=plot_pars["ylabel_fontsize"])  # 'Site selection\n(dN/dS)    '
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])
        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)
        ax_legend = axes[ax].legend()
        ax_legend.get_frame().set_linewidth(0.5)


    # Track 3D clustering
    # -------------------
    if "3d_clustering" in lst_tracks:
        ax = lst_tracks.index("3d_clustering")

        axes[ax].plot(np.array(range(protein_len-1))+1, o3d_df["O3D_score"], zorder=2, color=plot_pars["colors"]["o3d_score"], lw=1, label="Clustering score")
        axes[ax].fill_between(o3d_df['Pos'], 0, n_3d_max, where=(o3d_df['Cluster'] == 1),
                            color=plot_pars["colors"]["o3d_cluster"], alpha=1, label='Cluster', zorder=0, lw=1.5)

        axes[ax].set_ylabel('3D-clustering\n(missense mutations)', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])

        axes[ax].set_ylim(0 - n_3d_max/(plot_pars["j_margin"] * 5), n_3d_max + n_3d_max/plot_pars["j_margin"])
        axes[ax].yaxis.set_major_locator(MaxNLocator(integer=True, nbins=3))
        axes[ax].tick_params(axis='y', labelsize=plot_pars["ticksize"])

        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)

        ax_legend = axes[ax].legend()
        ax_legend.get_frame().set_linewidth(0.5)
        # axes[ax].legend(fontsize=plot_pars["legend_fontsize"], frameon=plot_pars["legend_frameon"],
        #                 bbox_to_anchor=plot_pars["o3d_bbox_to_anchor"], title = "3D-clustering",
        #                 title_fontsize=plot_pars["legend_fontsize"])


    # Track exon selection
    # ---------------------
    if "Exon_selection" in lst_tracks:
        ax = lst_tracks.index("Exon_selection")

        exon_coverage = get_exon_depth_saturation(coverage_df, mut_df)
        exon_coverage = get_exon_mid_prot_pos(exon_coverage, protein_len)
        exon_info = exon_coverage[["EXON_RANK", "START_PROT_POS", "END_PROT_POS", "MID_PROT_POS"]]
        exon_selection_df = exon_selection_df.merge(exon_info)
        custom_legend = []

        for impact in ["missense", "truncating"]:

            exon_selection_impact = exon_selection_df[exon_selection_df["impact"] == impact]
            exon_selection_impact_hits = exon_selection_impact[exon_selection_impact["pvalue"] < thr_selection].reset_index(drop=True)
            axes[ax].scatter(
                exon_selection_impact["MID_PROT_POS"], exon_selection_impact["dnds"],
                zorder=3, color=plot_pars["colors"]["exon_selection"][impact], s=60, lw=0.1, ec="black"
                )
            axes[ax].scatter(
                exon_selection_impact_hits["MID_PROT_POS"], exon_selection_impact_hits["dnds"],
                zorder=3, color=plot_pars["colors"]["exon_selection"][impact], s=60, lw=1, ec="black"
                )

            exon_selection_impact = exon_selection_impact.merge(exon_info, how="outer").fillna(0).sort_values("EXON_RANK").reset_index(drop=True)
            exon_selection_impact["impact"] = impact
            custom_legend.append(Line2D([0], [0], linestyle='-', marker='o', color=plot_pars["colors"]["exon_selection"][impact],
                                        label=impact.capitalize(), markersize=7.5, markeredgewidth=0.1, markeredgecolor="black"))

            for i, exon in exon_selection_impact.iterrows():
                start_pos = exon.START_PROT_POS
                end_pos = exon.END_PROT_POS
                start_dnds = 0 if i == 0 else exon_selection_impact.iloc[i-1].dnds
                dnds = exon.dnds
                axes[ax].vlines(start_pos, ymin=start_dnds, ymax=dnds, lw=1, zorder=1, alpha=1, color=plot_pars["colors"]["exon_selection"][impact])
                axes[ax].hlines(dnds, xmin=start_pos, xmax=end_pos, lw=1, zorder=1, alpha=1, color=plot_pars["colors"]["exon_selection"][impact])

        axes[ax].set_ylabel('Exon selection\n(dN/dS)', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)

        custom_legend.append(Line2D([0], [0], linestyle='-', marker='o', color="white",
                                    label="Significant", markersize=7.5, markeredgewidth=1.1, markeredgecolor="black"))
        ax_legend = axes[ax].legend(handles=custom_legend)
        ax_legend.get_frame().set_linewidth(0.5)


    # Track exon saturation
    # ---------------------
    if "Exon_saturation" in lst_tracks:
        ax = lst_tracks.index("Exon_saturation")

        norm = mcolors.Normalize(vmin=0, vmax=1)
        cmap = plt.get_cmap("hot_r")
        colors = cmap(norm(exon_coverage["SATURATION"].values))

        for i in range(len(exon_coverage)):
            exon = exon_coverage.iloc[i]
            start_pos = int(exon.START_PROT_POS)
            end_pos = int(exon.END_PROT_POS)

            if exon.DEPTH == 0:
                axes[ax].fill_between(
                    [start_pos, end_pos],
                    -0.2, 1.2,
                    color = plot_pars["colors"]["depth_below_thr"],
                    label=" ")
            else:
                axes[ax].fill_between(
                    [start_pos, end_pos],
                    -0.2, 1.2, color=colors[i])

            axes[ax].vlines(start_pos, ymin=-0.2, ymax=1.2, lw=0.4, zorder=1, alpha=1, color="black")
            if i == len(exon_coverage)-1:
                axes[ax].vlines(end_pos, ymin=-0.2, ymax=1.2, lw=0.4, zorder=1, alpha=1, color="black")

        axes[ax].set_yticks([])
        axes[ax].set_ylabel('Exon saturation', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].set_ylim(-0.2, 1.2)

        color_bar(
            axes,
            ax,
            cmap,
            norm,
            general_ticks=False,
            bbox_to_anchor=plot_pars["depth_res_bbox_to_anchor"],
            label='Exon saturation',
            fontsize=plot_pars["legend_fontsize"],
            tick_size=plot_pars["ticksize"],
            tick_span=0,
            borderpad=0)


    # Track exon depth
    # ----------------
    if "Exon_depth" in lst_tracks:
        ax = lst_tracks.index("Exon_depth")

        norm = mcolors.Normalize(vmin=200, vmax=max_depth)
        cmap = plt.get_cmap("viridis_r")
        colors = cmap(norm(exon_coverage["DEPTH"].values))
        for i in range(len(exon_coverage)):
            exon = exon_coverage.iloc[i]
            start_pos = int(exon.START_PROT_POS)
            end_pos = int(exon.END_PROT_POS)

            if exon.DEPTH == 0:
                axes[ax].fill_between(
                    [start_pos, end_pos],
                    -0.2, 1.2,
                    color = plot_pars["colors"]["depth_below_thr"],
                    label=" ")
            else:
                axes[ax].fill_between(
                    [start_pos, end_pos],
                    -0.2, 1.2, color=colors[i])

            axes[ax].vlines(start_pos, ymin=-0.2, ymax=1.2, lw=0.4, zorder=1, alpha=1, color="black")
            if i == len(exon_coverage)-1:
                axes[ax].vlines(end_pos, ymin=-0.2, ymax=1.2, lw=0.4, zorder=1, alpha=1, color="black")
            # if "Res_depth" in lst_tracks:
            #     axes[ax+1].vlines(start_pos, ymin=-0.2, ymax=1.2, lw=0.4, zorder=3, alpha=1, color="black")
            #     if i == len(exon_coverage)-1:
            #         axes[ax+1].vlines(end_pos, ymin=-0.2, ymax=1.2, lw=0.4, zorder=3, alpha=1, color="black")

        axes[ax].set_yticks([])
        axes[ax].set_ylabel('Exon depth', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].set_ylim(-0.2, 1.2)

        # Legend
        if (exon_coverage['DEPTH'] == 0).any():
            gray_patch = mpatches.Patch(color=plot_pars["colors"]["depth_below_thr"])

            legend_depth = axes[ax].legend(handles=[gray_patch],
                                        fontsize=plot_pars["legend_fontsize"],
                                        frameon=plot_pars["legend_frameon"],
                                        bbox_to_anchor=plot_pars["depth_0_bbox_to_anchor"],
                                        title="Depth < 200  ", title_fontsize=plot_pars["legend_depth_fontsize"])

            plt.setp(legend_depth.get_title(), fontsize=plot_pars["legend_fontsize"])
            handle = legend_depth.get_patches()[0]  # For patches (like bars in bar plots)
            handle.set_edgecolor('black')
            handle.set_linewidth(0.75)

        if "Res_depth" not in lst_tracks:
            color_bar(
                axes,
                ax,
                cmap,
                norm,
                general_ticks=False,
                bbox_to_anchor=plot_pars["depth_res_bbox_to_anchor"],
                label='Depth',
                fontsize=plot_pars["legend_fontsize"],
                tick_size=plot_pars["ticksize"],
                tick_span=0,
                borderpad=0)


    # Track per-residue depth
    # -----------------------
    if "Res_depth" in lst_tracks:
        ax = lst_tracks.index("Res_depth")

        res_coverage = get_res_coverage(coverage_df)
        res_coverage["DEPTH"] = res_coverage["DEPTH"].apply(lambda x: np.nan if x == 0 else x)
        norm = mcolors.Normalize(vmin=200, vmax=max_depth)
        cmap = plt.get_cmap("viridis_r")
        colors = cmap(norm(res_coverage["DEPTH"].values))

        # Fill for depth = 0
        axes[ax].fill_between(
            res_coverage['PROT_POS'], -0.2, 1.2, where=where_plus(res_coverage['DEPTH'].fillna(0) == 0),
            color=plot_pars["colors"]["depth_below_thr"], alpha=1, label=' '
            )

        # Fill for depth > 0
        for i in range(len(res_coverage) - 1):
            axes[ax].fill_between(
                [res_coverage["PROT_POS"].iloc[i], res_coverage["PROT_POS"].iloc[i+1]],
                -0.2, 1.2,
                color=colors[i]
            )

        axes[ax].set_yticks([])
        axes[ax].set_ylabel('Residue depth', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].set_ylim(-0.2, 1.2)
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])

        color_bar(
            axes,
            ax,
            cmap,
            norm,
            general_ticks=False,
            bbox_to_anchor=plot_pars["depth_res_bbox_to_anchor"],
            label='Depth',
            fontsize=plot_pars["legend_fontsize"],
            tick_size=plot_pars["ticksize"],
            tick_span=0,
            borderpad=0)

        if "Exon_depth" not in lst_tracks and (res_coverage['DEPTH'].fillna(0) == 0).any():
            gray_patch = mpatches.Patch(color=plot_pars["colors"]["depth_below_thr"])

            legend_depth = axes[ax].legend(handles=[gray_patch],
                                        fontsize=plot_pars["legend_fontsize"],
                                        frameon=plot_pars["legend_frameon"],
                                        bbox_to_anchor=plot_pars["depth_0_bbox_to_anchor"],
                                        title="Depth < 200  ", title_fontsize=plot_pars["legend_depth_fontsize"])

            plt.setp(legend_depth.get_title(), fontsize=plot_pars["legend_fontsize"])
            handle = legend_depth.get_patches()[0]  # For patches (like bars in bar plots)
            handle.set_edgecolor('black')
            handle.set_linewidth(0.75)


    # Track pACC
    # ----------
    if "Solvent_accessibility" in lst_tracks and isinstance(pdb_tool_df, pd.DataFrame):
        ax = lst_tracks.index("Solvent_accessibility")

        max_pacc = np.max(pdb_tool_df["pACC"].fillna(0))
        axes[ax].fill_between(
            pdb_tool_df["Pos"], 0, pdb_tool_df["pACC"].fillna(0),
            zorder=2, color=sns.color_palette("pastel")[4], alpha=0.35
            )
        axes[ax].plot(
            pdb_tool_df['Pos'], pdb_tool_df["pACC"].fillna(0),
            label="pACC", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5
            )
        axes[ax].set_ylabel('Solvent\naccessibility', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)
        axes[ax].set_ylim(0 - max_pacc/30, max_pacc + max_pacc/10)


    # Plot stability change
    # ---------------------
    if "Stability_change" in lst_tracks and isinstance(ddg_df, pd.DataFrame):
        ax = lst_tracks.index("Stability_change")

        max_value, min_value = ddg_df["Stability_change"].max(), ddg_df["Stability_change"].min()
        axes[ax].fill_between(ddg_df['Pos'], 0, ddg_df["Stability_change"], zorder=1,
                                color="white")
        axes[ax].fill_between(ddg_df['Pos'], 0, ddg_df["Stability_change"], zorder=1,
                                color=sns.color_palette("pastel")[4], alpha=0.35)
        axes[ax].plot(ddg_df['Pos'], ddg_df["Stability_change"],
                        label="Stability change", zorder=2, color=sns.color_palette("tab10")[4], lw=0.5)
        axes[ax].set_ylabel('Stability change (kcal/mol)', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].spines['top'].set_visible(False)
        axes[ax].spines['right'].set_visible(False)


    # Track SSE
    # ---------
    if "Secondary_structure" in lst_tracks and isinstance(pdb_tool_df, pd.DataFrame):
        ax = lst_tracks.index("Secondary_structure")

        for n, sse in enumerate(('Coil', 'Helix', 'Ladder')):
            axes[ax].fill_between(
                pdb_tool_df["Pos"], 0, 1, where=(pdb_tool_df['SSE'] == sse),
                zorder=2, color=plot_pars["sse_colors"][sse],
                alpha=1, label=sse, lw=plot_pars["sse_lw"]
                )
        axes[ax].set_yticks([])
        axes[ax].legend(fontsize=plot_pars["legend_fontsize"], frameon=plot_pars["legend_frameon"],
                        bbox_to_anchor=plot_pars["sse_bbox_to_anchor"], title = "Secondary structure",
                        title_fontsize=plot_pars["legend_fontsize"],
                        handleheight=0.65, handlelength=2)
        axes[ax].set_ylabel('Secondary structure', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].set_ylim(0, 1)


    # Track Domain
    # ------------
    if "Domain" in lst_tracks and isinstance(domain_df, pd.DataFrame):
        ax = lst_tracks.index("Domain")

        domain_color_dict = {}
        for n, name in enumerate(domain_df["DOMAIN_ID"].unique()):
            domain_color_dict[name] = f"C{n}"

        n = 0
        added_domain = []
        for i, row in domain_df.iterrows():
            if pd.Series([row["DOMAIN_ID"], row["Begin"], row["End"]]).isnull().any():
                continue

            name = row["DOMAIN_ID"]
            start = int(row["Begin"])
            end = int(row["End"])

            if name not in added_domain and (protein_len >= plot_pars["len_txt_thr"]):
                axes[ax].fill_between(range(start, end + 1), -0.5, 0.45, alpha=0.5, color=domain_color_dict[name], label=name, lw=0.5)
            else:
                axes[ax].fill_between(range(start, end + 1), -0.5, 0.45, alpha=0.5, color=domain_color_dict[name], lw=0.5)

            if name not in added_domain:
                if protein_len < plot_pars["len_txt_thr"]:
                    y = -0.04
                    axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                added_domain.append(name)
        axes[ax].set_yticks([])
        if (protein_len >= plot_pars["len_txt_thr"]) and len(added_domain) > 0:
            pass
            # axes[ax].legend(fontsize=plot_pars["legend_fontsize"], frameon=plot_pars["legend_frameon"],
            #                 bbox_to_anchor=(plot_pars["domain_x_bbox_to_anchor"][gene], plot_pars["domain_y_bbox_to_anchor"]), title = "Domain",
            #                 title_fontsize=plot_pars["legend_fontsize"],
            #                 handleheight=0.67, handlelength=2)
        axes[ax].set_ylabel('Domain', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].set_ylim(-0.5, 0.45)


    # Track Domain Selection
    # ----------------------
    if "Domain_selection" in lst_tracks and isinstance(domain_selection_df, pd.DataFrame):
        ax = lst_tracks.index("Domain_selection")

        norm = mcolors.Normalize(vmin=1, vmax=domain_selection_df.dnds.max())
        cmap = plt.get_cmap("Reds")
        colors = cmap(norm(domain_selection_df["dnds"].values))

        for i, row in domain_selection_df.iterrows():
            if pd.Series([row["dnds"], row["Begin"], row["End"]]).isnull().any():
                continue

            dnds = row["dnds"]
            start = int(row["Begin"])
            end = int(row["End"])
            y = (0.02, 0.5) if row["impact"] == "missense" else (-0.5, -0.04)

            axes[ax].fill_between(range(start, end + 1), y[0], y[1], color=colors[i])

        # xlim = axes[ax].get_xlim()
        # axes[ax].hlines(0, 0-(protein_len*0.2), protein_len+(protein_len*0.2),  lw=0.7, zorder=1, alpha=1, color="black")
        # axes[ax].set_xlim(xlim)

        axes[ax].set_yticks([0.25, -0.25])
        axes[ax].set_yticklabels(["M", "T"])
        axes[ax].set_ylabel('Domain selection', fontsize=plot_pars["ylabel_fontsize"], rotation=0, va='center', ha='right')
        axes[ax].yaxis.set_label_coords(plot_pars["y_labels_coord"][0], plot_pars["y_labels_coord"][1])
        axes[ax].set_ylim(-0.5, 0.5)

        color_bar(
            axes,
            ax,
            cmap,
            norm,
            general_ticks=False,
            bbox_to_anchor=plot_pars["domain_sel_bbox_to_anchor"],
            label='Domain selection',
            fontsize=plot_pars["legend_fontsize"],
            tick_size=plot_pars["ticksize"],
            tick_span=0,
            borderpad=0)


    # Legend
    # ======

    handles_0, labels_0 = axes[0].get_legend_handles_labels()
    handles, labels = [], []
    for h, l in zip(handles_0, labels_0):
        if l not in labels:
            handles.append(h)
            labels.append(l)
    fig.legend(
        handles, labels, loc='upper right', fontsize=plot_pars["legend_fontsize"],
        frameon=plot_pars["legend_frameon"], bbox_to_anchor=plot_pars["cnsq_bbox_to_anchor"],
        title="Mutation consequence", title_fontsize=plot_pars["legend_fontsize"]
        )

    if gene is not None:
        fig.suptitle(title, y=0.93)

    if save:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

    if "Domain" in lst_tracks and isinstance(domain_df, pd.DataFrame):
        return domain_color_dict



def plot_domain_selection(
    domain_selection_in_gene,
    gene,
    color_map,
    title="Domain selection",
    show_domain_legend=False,
    legend2_coord=(1.05, 1),
    save=False,
    filename="domain_selection.png",
    pvalue_threshold = 0.05
    ):

    custom_markers = {"missense": "o", "truncating": "D"}  # 'o' = Circle, 'D' = Rotated Square
    custom_size = {"missense": 200, "truncating": 145}  # Larger for missense


    domain_selection_in_gene = domain_selection_in_gene.drop(columns=["Begin", "End"]).drop_duplicates()

    # this dictionary is defined in the gene plot with the gene--domain name
    domain_selection_in_gene["color"] = domain_selection_in_gene["DOMAIN_ID"].map(color_map)

    domain_selection_in_gene["edge_width"] = domain_selection_in_gene["pvalue"].apply(lambda p: 1.5 if p < pvalue_threshold else 0.2)

    fig = plt.figure(figsize=(5, 2.5))

    # H-line and V-line
    plt.xlim(domain_selection_in_gene["x_pos"].min() - 0.5, domain_selection_in_gene["x_pos"].max() + 0.5)
    plt.hlines(1, *plt.xlim(), lw=1, zorder=1, alpha=1, color="gray", linestyles="dashed")
    plt.vlines(domain_selection_in_gene["x_pos"], ymin=domain_selection_in_gene["lower"], ymax=domain_selection_in_gene["upper"], lw=1, zorder=1, alpha=1, color="black")

    # Scatter plot
    for gene_domain_name, color in color_map.items():
        for impact, marker in custom_markers.items():
            subset = domain_selection_in_gene[(domain_selection_in_gene["DOMAIN_ID"] == gene_domain_name)
                                                & (domain_selection_in_gene["impact"] == impact)]

            plt.scatter(
                subset["x_pos"], subset["dnds"], lw=0, edgecolor="black",
                color="white", marker=marker, s=custom_size[impact]
            )
            plt.scatter(
                subset["x_pos"], subset["dnds"],
                color=color, marker=marker, s=custom_size[impact],
                edgecolor="black", lw=subset["edge_width"], alpha=0.7
            )

    ax = plt.gca()
    ax.set_xticks(range(len(domain_selection_in_gene["DOMAIN_ID"].unique())))
    ax.set_xticklabels(domain_selection_in_gene["DOMAIN_ID"].unique(), rotation=45, ha="right")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel("dN/dS")
    plt.title(title)

    # Legends
    selection_legend = [mpatches.Patch(color=color, label=domain_id) for domain_id, color in color_map.items()]
    custom_size_legend = {"missense": 100, "truncating": 70}
    impact_legend = [
        plt.scatter(
            [], [], marker=custom_markers[impact], s=custom_size_legend[impact],
            color="white", edgecolor="black", label=impact
            ) for impact in custom_markers
        ]


    legend2 = plt.legend(handles=impact_legend, title="Impact" if show_domain_legend else None, bbox_to_anchor=legend2_coord, loc="upper left", frameon=False)
    if show_domain_legend:
        legend1 = plt.legend(handles=selection_legend, title="Domain ID", bbox_to_anchor=(1.05, 0.65), loc="upper left", frameon=False)
        ax.add_artist(legend1)

    if save:
        fig.savefig(filename, dpi=300, bbox_inches='tight')

    plt.show()


def mut_count_horizontal_barplot(
    df,
    color_dict,
    figsize=(3.5, 0.25),
    title=None,
    save=False,
    filename="mut_count_horizontal_barplot.png"
    ):
    order = ["truncating", "missense", "synonymous"]
    df = df.set_index("Consequence").loc[order]#.reset_index()

    fig, ax = plt.subplots(figsize=figsize)
    bottom = 0 #np.zeros(len(df))  # Track bottom for stacking

    for consequence in order:
        bars = ax.barh(0, df.loc[consequence], left=bottom, color=color_dict[consequence], label=consequence.capitalize())

        for bar in bars:
            width = 0.75
            if df.loc[consequence].item() > 0:  # Only label nonzero values
                ax.text(bottom + df.loc[consequence].item() / 2 , width, str(int(df.loc[consequence])),
                        ha='center', va='center', fontsize=10)

        bottom += df.loc[consequence].item()  # Update bottom for stacking

    ax.set_xlabel("Number of mutations")
    # ax.set_title(title if title else "Mutation Count per Sample", pad=10)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticks([])  # Remove y-axis ticks
    ax.set_yticklabels([])  # Remove y-axis ticks
    # ax.legend()

    if save:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()


def get_selection_groups(df, thr_selection):

    df = df.copy().rename(columns={"Protein_position": "Pos"})
    g0 = df[df["Selection"] == 0].copy()
    g1 = df[(df["Selection"] != 0) & (df["p_value"] >= thr_selection)].copy()

    g2 = df[df["p_value"] < thr_selection].copy()
    g2 = g2.sort_values("Selection").reset_index(drop=True)

    g0["Group"] = "G0"
    g1["Group"] = "G1"
    g2["Group"] = "G2"

    return pd.concat((g0, g1, g2)).reset_index(drop=True)


def plot_feat_by_selection_group(df, figsize=(10, 4), save=False, filename="selection_groups_feat.png",):

    color_dict = {
    "G0": "gray",
    "G1": "#6baed6",
    "G2": "blue"
}
    rename_labels = {
        "G0": f"Not observed\n(N={df[df['Group'] == 'G0'].shape[0]})",
        "G1": f"Not significant\n(N={df[df['Group'] == 'G1'].shape[0]})",
        "G2": f"Significant\n(N={df[df['Group'] == 'G2'].shape[0]})"
    }

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=False)

    # Violin + Jitter plot for pACC
    sns.violinplot(ax=axes[0], x="Group", y="pACC", data=df, color="#d9d9d9", inner=None)
    sns.stripplot(ax=axes[0], x="Group", y="pACC", data=df,
                palette=color_dict, jitter=True, size=4, alpha=0.7)
    axes[0].set_title("Solvent Accessibility by Selection")
    axes[0].set_xlabel(None)
    axes[0].set_ylabel("Solvent accessibility")
    axes[0].set_ylim(-2,102)

    # Violin + Jitter plot for Stability Change
    sns.violinplot(ax=axes[1], x="Group", y="Stability_change", data=df, color="#d9d9d9", inner=None)
    sns.stripplot(ax=axes[1], x="Group", y="Stability_change", data=df,
                palette=color_dict, jitter=True, size=4, alpha=0.7)
    axes[1].set_title("Stability Change by Selection")
    axes[1].set_xlabel(None)
    axes[1].set_ylabel("Stability change")

    # Overlay median values
    medians = df.groupby("Group")["pACC"].median()
    for i, median in enumerate(medians):
        axes[0].scatter(i, median, zorder=3, color='black', marker='x', s=80, label="Median" if i == 0 else "")
    medians = df.groupby("Group")["Stability_change"].median()
    for i, median in enumerate(medians):
        axes[1].scatter(i, median, zorder=3, color='black', marker='x', s=80, label="Median" if i == 0 else "")

    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    axes[0].set_xticklabels([rename_labels[label] for label in color_dict.keys()])
    axes[1].set_xticklabels([rename_labels[label] for label in color_dict.keys()])

    # # Improve layout
    # plt.tight_layout()

    if save:
        fig.savefig(filename, dpi=300, bbox_inches='tight')

    plt.show()


def plot_all_domain_selection(df, color_map, figsize=(10, 3), show_domain_legend=True, save=False, filename="all_domain_selection.png"):

    custom_markers = {"missense": "o", "truncating": "D"}  # 'o' = Circle, 'D' = Rotated Square
    custom_size = {"missense": 150, "truncating": 100}  # Larger for missense
    df["color"] = df["domain_id"].map(color_map) 
    df["edge_width"] = df["pvalue"].apply(lambda p: 1.5 if p < 0.05 else 0.2)

    fig = plt.figure(figsize=figsize)

    # H-line and V-line
    plt.xlim(df["x_pos"].min() - 0.5, df["x_pos"].max() + 0.5) 
    plt.hlines(1, *plt.xlim(), lw=1, zorder=1, alpha=1, color="gray", linestyles="dashed")
    plt.vlines(df["x_pos"], ymin=df["lower"], ymax=df["upper"], lw=1, zorder=1, alpha=1, color="black")

    # Scatter plot
    for gene, color in color_map.items():
        for impact, marker in custom_markers.items():
            subset = df[(df["gene"] == gene) & (df["impact"] == impact)]
            plt.scatter(
                subset["x_pos"], subset["dnds"], lw=0,
                color="white", marker=marker, s=custom_size[impact]
            )
            plt.scatter(
                subset["x_pos"], subset["dnds"],
                color=color, marker=marker, s=custom_size[impact], 
                edgecolor="black", lw=subset["edge_width"], alpha=0.7
            )

    ax = plt.gca()
    ax.set_xticks(range(len(df["domain_id"].unique()))) 
    ax.set_xticklabels(df["domain_id"].unique(), rotation=45, ha="right") 
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlabel("Domain ID")
    plt.ylabel("dN/dS")
    plt.title(f"Top  domain selection")

    # Legends
    selection_legend = [mpatches.Patch(color=color, label=domain_id) for domain_id, color in color_map.items()]
    custom_size_legend = {"missense": 100, "truncating": 70}
    impact_legend = [plt.scatter(
        [], [], marker=custom_markers[impact], s=custom_size_legend[impact], 
        color="white", edgecolor="black", label=impact) for impact in custom_markers]

    legend1 = plt.legend(handles=selection_legend, title="Gene", bbox_to_anchor=(1.05, 0.65), loc="upper left", frameon=False)
    legend2 = plt.legend(handles=impact_legend, title="Impact", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)
    
    if show_domain_legend:
        plt.gca().add_artist(legend1)
        
    if save:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    
def format_domain_selection(df, xshift=0.15, sort_by="average_dnds", palette=sns.color_palette("Paired")):
    
    df = df.copy()
    df["domain_id"] = df["gene"].str.split("--").apply(lambda x : x[1])
    df["gene"] = df["gene"].str.split("--").apply(lambda x: x[0])

    if sort_by == "average_dnds":
        sort_map = df.groupby("domain_id").apply(lambda x: x["dnds"].mean()).to_dict()
    else:
        sort_map = df.groupby("domain_id").apply(lambda x: x[x["impact"] == sort_by]["dnds"].max()).to_dict()

    df["sort"] = df["domain_id"].map(sort_map)
    df = df.sort_values(["sort", "impact"], ascending=[False, True]).reset_index(drop=True)

    df["x_pos"] = df["domain_id"].astype("category").cat.set_categories(df["domain_id"].unique()).cat.codes
    df["x_pos"] += np.where(df["impact"] == "missense", -xshift, xshift)
    
    gene_color_dict = {gene: mcolors.to_hex(color) for gene, color in zip(df.gene.unique(), palette)}
    
    return df, gene_color_dict


# TODO:
# - Check that max res depth equal max exon depth, if not use the one that's actually included in the plot

def plotting_wrapper(maf, exons_depth, o3d_df, exon_selection, domain_selection, site_selection, list_genes, output_directory = '.', o3d_seq_df=None, o3d_pdb_tool_df=None, domain=None):
    # for gene in ["TP53"]:
    for gene in list_genes:
        try:
            plotting_single_gene(gene, maf, exons_depth, o3d_df, exon_selection, domain_selection, site_selection, output_directory, o3d_seq_df, o3d_pdb_tool_df, domain)
        except Exception as e:
            print(f"Saturation plots for gene: {gene} did not work.")
            print("Error", e)


    try :
        # TODO revise if here some lines can be skipped,
        # and the formatting of domain selection can be done only once
        df, gene_color_dict = format_domain_selection(domain_selection, xshift=0.22)
        df = df.copy()[:60]
        plot_all_domain_selection(df, gene_color_dict, figsize=(16,6), save=True, filename=f"{output_directory}/domain_selection.missense_truncating.pdf")
        for impact in ["missense", "truncating"]:
            df, gene_color_dict = format_domain_selection(domain_selection[domain_selection["impact"] == impact], sort_by=impact, xshift=0)[:60]
            plot_all_domain_selection(df, gene_color_dict, figsize=(14,6), save=True, filename=f"{output_directory}/domain_selection.{impact}.pdf")

    except Exception as e:
        print(f"Plots for domain selection at cohort level did not work.")
        print("Error", e)




indels = False

def plotting_single_gene(gene, maf, exons_depth, o3d_df,
                            exon_selection, domain_selection, site_selection,
                            output_dir, o3d_seq_df, o3d_pdb_tool_df, domain):

    # Data
    # ====

    # Mut
    gene_mut = maf[maf["Gene"] == gene].drop(columns="Gene").reset_index(drop=True)
    gene_mut_count = gene_mut.groupby(['Consequence', 'Pos']).size().reset_index(name='Count')
    gene_mut_cnsq_count = gene_mut.groupby(['Consequence']).size().reset_index(name='Count')

    # Oncodrive3D
    o3d_gene_df = get_o3d_gene_data(gene, o3d_seq_df, o3d_df)
    uni_id = o3d_seq_df[o3d_seq_df["Gene"] == gene]["Uniprot_ID"].values[0]
    pdb_tool_gene = o3d_pdb_tool_df[o3d_pdb_tool_df["Uniprot_ID"] == uni_id].reset_index(drop=True)

    if not indels :
        gene_mut_count = gene_mut_count[gene_mut_count["Consequence"] != "indel"].reset_index(drop=True)
        gene_mut = gene_mut[gene_mut["Consequence"] != "indel"].reset_index(drop=True)

    prot_len = int(exons_depth[exons_depth["GENE"] == gene]["PROT_POS"].max())
    domain_gene = domain[domain["GENE"] == gene].reset_index(drop=True)

    all_exons_depth = (
                        exons_depth
                        .groupby(["GENE", "EXON_RANK"], as_index=False)
                        .agg(total_depth=("DEPTH", "sum"), total_covered=("COVERED", "sum"))
                        .assign(DEPTH=lambda df: np.where(df["total_covered"] == 0, np.nan, df["total_depth"] / df["total_covered"]))
                        .drop(columns=["total_depth", "total_covered"])
                    )

    max_depth = all_exons_depth["DEPTH"].max()

    gene_exons_depth = exons_depth[exons_depth["GENE"] == gene].sort_values("DNA_POS")
    gene_exons_depth = gene_exons_depth.dropna(subset=["EXON_RANK"]).reset_index(drop=True)

    exon_selection_gene = exon_selection[exon_selection["gene"].str.startswith(gene)]
    exon_selection_gene = exon_selection_gene.sort_values("exon_rank").reset_index(drop=True).rename(columns={"exon_rank" : "EXON_RANK"})
    exon_selection_gene = exon_selection_gene[["EXON_RANK", "impact", "dnds", "lower", "upper", "pvalue"]]

    site_selection_gene = site_selection[site_selection["GENE"] == gene].sort_values("Protein_position").reset_index(drop=True)

    domain_selection_gene = get_domain_selection_gene(domain_selection, domain_gene, gene)

    ddg_gene = None
    # ddg_gene = get_ddg_gene(gene, gene_mut, o3d_annotations, o3d_seq_df)


    # Plot
    # ====

    plot_pars["fsize"] = (12.3, 16)
    lst_tracks=[
        "Mut_count",
        "3d_clustering",
        "Site_selection",
        "Exon_selection",
        "Exon_saturation",
        "Exon_depth",
        "Res_depth",
        "Solvent_accessibility",
        "Stability_change",
        "Secondary_structure",
        "Domain",
        "Domain_selection"
        ]

    # FIXME
    # this should point to a list of all the genes for which these tracks are not available,
    # otherwise, add a check that makes sure the gene is/isnot suitable for having these tracks
    if gene in ["KDM6A", "STAG2"]:
        plot_pars["fsize"] = (12.3, 12)
        [lst_tracks.remove(track) for track in [
            "3d_clustering",
            "Cancer_3d_clustering",
            "Solvent_accessibility",
            "Stability_change",
            "Secondary_structure",
            ] if track in lst_tracks
        ]

    domain_color_dict_gene = plot_gene_selection(
        mut_count_df=gene_mut_count,
        mut_df=gene_mut,
        o3d_df=o3d_gene_df,
        pdb_tool_df=pdb_tool_gene,
        domain_df=domain_gene,
        coverage_df=gene_exons_depth,
        site_selection_df=site_selection_gene,
        exon_selection_df=exon_selection_gene,
        domain_selection_df=domain_selection_gene,
        ddg_df=ddg_gene,
        protein_len=prot_len,
        gene=gene,
        max_depth=max_depth,
        plot_pars=plot_pars,
        title=f"{gene}",
        lst_tracks=lst_tracks,
        thr_selection=0.00001,
        default_track_order=False,
        save=True,
        filename=f"{output_dir}/{gene}.saturation_all.png"
        )




    # Side plots
    # ==========

    if domain_selection_gene.shape[0] > 0:
        plot_domain_selection(
            domain_selection_gene,
            gene,
            domain_color_dict_gene,
            title=f"{gene}\nDomain selection",
            legend2_coord=(0.75, 1),
            save=True,
            filename=f"{output_dir}/{gene}.domain_selection.png"
            )
    else:
        print("Not plotting domains, no omega value for domains")


    mut_count_horizontal_barplot(
        gene_mut_cnsq_count, plot_pars["colors"],
        title=f"{gene}",
        save=True,
        filename=f"{output_dir}/{gene}.normal.stacked_horizontal.png"
        )



    # Selection groups feat
    # =====================
    site_selection_gene_grouped = get_selection_groups(site_selection_gene, thr_selection=0.00001)
    site_selection_gene_grouped = site_selection_gene_grouped.merge(pdb_tool_gene[["Pos", "pACC"]])

    if ddg_gene is not None:
        site_selection_gene_grouped = site_selection_gene_grouped.merge(ddg_gene)
        plot_feat_by_selection_group(site_selection_gene_grouped, save=True, filename=f"{output_dir}/{gene}.selection_groups.png")



#####
# Load cohort data
#####
def data_loading(sample_name, domain, exons_depth):

    mutations_file = f"{sample_name}.somatic.mutations.tsv"

    # Count each mutation only ones if it appears in multiple reads
    maf, list_of_genes = get_normal_maf(mutations_file, truncating=True)


    ## Positive selection
    site_selection = f"{sample_name}.aminoacid.comparison.tsv.gz"

    omega_file = f"output_mle.{sample_name}.global_loc.tsv"
    o3d_df_file = f"{sample_name}.3d_clustering_pos.csv"


    # Omega data
    omega_table = pd.read_table(omega_file)

    omega_subgenic_table = omega_table[(omega_table["gene"].str.contains("--")) &
                                        (omega_table["impact"].isin(["missense", "truncating"]))
                                    ].reset_index(drop=True)

    # Exon selection
    omega_exon = omega_subgenic_table[omega_subgenic_table["gene"].str.contains('--exon_')]
    exon_selection = omega_exon.merge(
        exons_depth.rename(columns={"EXON_ID" : "gene", "EXON_RANK" : "exon_rank"}).dropna()[["gene", "exon_rank"]].drop_duplicates().reset_index(drop=True),
        how="left").sort_values(["gene", "exon_rank", "impact"]).reset_index(drop=True)
    print("> exon_selection:", exon_selection.shape )

    # Domain selection
    domain_selection = omega_subgenic_table.copy()
    domain_selection = domain_selection[domain_selection["gene"].isin(domain["NAME"].values)].reset_index(drop=True)
    print("> domain_selection:", domain_selection.shape )


    # Per-site selection
    site_selection = pd.read_table(site_selection).rename(
        columns={"OBS/EXP" : "Selection"}).drop(
        columns=["OBSERVED_MUTS", "EXPECTED_MUTS"])
    site_selection.loc[site_selection["Selection"] < 0, "Selection"] = 0
    site_selection = site_selection[site_selection["Protein_position"] != "-"].reset_index(drop=True)
    site_selection["Protein_position"] = site_selection["Protein_position"].astype(int)


    # Oncodrive3D data
    o3d_df = pd.read_csv(o3d_df_file)[["Gene", "Pos", "Score", "Score_obs_sim", "pval", "C", "C_ext"]]

    return maf, exons_depth, o3d_df, exon_selection, domain_selection, site_selection, list_of_genes



def get_reference_data(domain_file, exons_depths):
    #####
    # Get reference data
    #####

    # Oncodrive3D
    o3d_seq_df = pd.read_table("seq_for_mut_prob.tsv")
    o3d_pdb_tool_df = pd.read_table("pdb_tool_df.tsv")

    # Domain annotations
    domain = pd.read_table(domain_file, sep="\t", low_memory=False)

    exons_depths_df = pd.read_table(exons_depths, sep="\t", low_memory=False)

    return o3d_seq_df, o3d_pdb_tool_df, domain, exons_depths_df



@click.command()
@click.option('--sample_name', type=str, default='all_samples', help='Name of the sample being processed.')
@click.option('--outdir', type=click.Path(), help='Output path for plots')
@click.option('--domain_file', type=click.Path(), help='Path to the domain file')
@click.option('--exons_depths', type=click.Path(), help='Path to the exons depths file')
def main(sample_name, outdir, domain_file, exons_depths):
    click.echo("Starting to run plot saturation results...")

    try :
        o3d_seq_df, o3d_pdb_tool_df, domain, exons_depths_df = get_reference_data(domain_file, exons_depths)
        click.echo("Reference data loaded...")
        maf, exons_depth, o3d_df, exon_selection, domain_selection, site_selection, gene_list = data_loading(sample_name, domain, exons_depths_df)
        click.echo("Sample data loaded...")
        plotting_wrapper(maf, exons_depth, o3d_df, exon_selection, domain_selection, site_selection, gene_list, outdir, o3d_seq_df, o3d_pdb_tool_df, domain)

    except Exception as e:
        print("Error in the process", e)

if __name__ == '__main__':
    main()

