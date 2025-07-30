#!/usr/bin/env python

import click
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import norm

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from read_utils import custom_na_values
from utils_plot import metrics_colors_dictionary
from plot_selection_omega import plot_omega_vertical, build_counts_from_df_complete



def generate_all_side_figures(sample,
                                outdir = '.',
                                gene_list = None,
                                tools = ["oncodrivefml", "omega_trunc", "omega_mis", "excess_indels"]
                                ):

    maf = pd.read_table(f"{sample}.somatic.mutations.tsv", na_values = custom_na_values)
    snvs_maf = maf[maf["TYPE"] == "SNV"].reset_index(drop = True)

    possible_genes = []
    if "oncodrivefml" in tools:
        oncodrivefml_data = pd.read_table(f"{sample}-oncodrivefml.tsv.gz")
        oncodrivefml_data = oncodrivefml_data[["GENE_ID", "Z-SCORE", "Q_VALUE", "AVG_SCORE_OBS", "POPULATION_MEAN", "STD_OF_MEANS"]]
        oncodrivefml_data.columns = ["GENE", "OncodriveFML", "pvalue", "OBSERVED_MEAN", "BACKGROUND_MEAN", "BACKGROUND_STD"]
        oncodrivefml_genes = list(pd.unique(oncodrivefml_data["GENE"]))
        possible_genes += oncodrivefml_genes


    if "omega_trunc" in tools or "omega_mis" in tools:
        omega_data = pd.read_table(f"output_mle.{sample}.tsv")
        omega_data = omega_data[(omega_data["impact"].isin(['missense', 'truncating']))
                                    & ~(omega_data["gene"].str.contains('--'))        # select only genes
                                ]
        if "omega_trunc" in tools :
            omega_truncating = omega_data[omega_data["impact"] == "truncating"].reset_index(drop = True)[["gene", "mutations", "dnds", "pvalue", "lower", "upper"]]
            omega_truncating.columns = ["GENE", "mutations_trunc", "omega_trunc", "pvalue", "lower", "upper"]
            omega_truncating_genes = list(pd.unique(omega_truncating["GENE"]))
            possible_genes += omega_truncating_genes

        if "omega_mis" in tools :
            omega_missense = omega_data[omega_data["impact"] == "missense"].reset_index(drop = True)[["gene", "mutations", "dnds", "pvalue", "lower", "upper"]]
            omega_missense.columns = ["GENE", "mutations_mis", "omega_mis", "pvalue", "lower", "upper"]
            omega_missense_genes = list(pd.unique(omega_truncating["GENE"]))
            possible_genes += omega_missense_genes


    if "excess_indels" in tools:
        indels_data = pd.read_table(f"{sample}.sample.indels.tsv",
                                        sep = '\t',
                                        header = 0)
        indels_genes = list(pd.unique(indels_data["SYMBOL"]))
        possible_genes += indels_genes

    gene_list = list(set(possible_genes).intersection(set(snvs_maf["canonical_SYMBOL"].unique())))

    for genee in gene_list:
        print(genee, end = '\t')
        try :

            if "oncodrivefml" in tools:
                # there is no run of oncodrivefml with ALL_GENES
                if genee in oncodrivefml_genes:
                    oncodrivefml_gene_data = oncodrivefml_data[oncodrivefml_data["GENE"] == genee].to_dict(orient='records')[0]

                    fig_gene_fml = plot_oncodrivefml_side(oncodrivefml_gene_data)
                    fig_gene_fml.savefig(f"{outdir}/{genee}.{sample}.oncodrivefml.pdf", bbox_inches='tight', dpi = 100)
                    plt.show()
                    plt.close()
                    print("ofml done", end = '\t')

            if "omega_trunc" in tools:
                if genee in omega_truncating_genes and genee in omega_missense_genes:
                    omega_df = build_counts_from_df_complete(genee, snvs_maf, omega_truncating, omega_missense)

                    fig_gene_omega = plot_omega_vertical(omega_df)
                    fig_gene_omega.savefig(f"{outdir}/{genee}.{sample}.omega.pdf", bbox_inches='tight', dpi = 100)
                    plt.show()
                    plt.close()
                    print("omega done", end = '\t')

            if "excess_indels" in tools:
                if genee in indels_genes:
                    indel_data_gene = indels_data[indels_data["SYMBOL"] == genee].to_dict(orient='records')[0]

                    fig_gene_indel = plotting_indels_side(indel_data_gene)
                    fig_gene_indel.savefig(f"{outdir}/{genee}.{sample}.indels.pdf", bbox_inches='tight', dpi = 100)
                    plt.show()
                    plt.close()
                    print("indels done")

        except Exception as exe:
            print("failed processing of")
            print(genee)
            print(exe)






def plot_oncodrivefml_side(geneee_data):
    legend_fontsize = 12
    xaxis_fontsize = 12

    # Extract the necessary values
    observed_mean = geneee_data['OBSERVED_MEAN']
    background_mean = geneee_data['BACKGROUND_MEAN']
    background_std = geneee_data['BACKGROUND_STD']
    p_value = geneee_data['pvalue']

    observed_color = metrics_colors_dictionary["ofml_score"]

    deviation = abs(observed_mean - background_mean)

    # Calculate the Z-score
    z_score = (observed_mean - background_mean) / background_std

    # Generate a range of values for the x-axis
    x = np.linspace(background_mean - 1*background_std - deviation, background_mean + 1*background_std + deviation, 1000)
    # Generate the normal distribution based on background mean and std
    y = norm.pdf(x, background_mean, background_std)

    # Compute half of the maximum value
    mid_y = max(y) / 2

    # Plot the normal distribution vertically
    fig, ax = plt.subplots(figsize=(3, 1.75))

    background_color = 'dimgrey'
    background_color_line = 'dimgrey'
    ax.plot(x, y, color = background_color)
    ax.fill_betweenx(y, x, color = background_color, alpha=0.2) #, label = 'Randomized\nmeans' )

    # Arrows with annotations
    arrow_props = dict(facecolor= observed_color, edgecolor = observed_color, arrowstyle='<->')

    if z_score > 0:
        ax.set_xlim(background_mean - (1*background_std + deviation) / 2, background_mean + (1*background_std + deviation))
        # ax.text(background_mean + background_std, mid_y * 1.75, f'Randomized\nmeans', color=background_color, ha='left', va='top')
        #ax.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.2, 1.1), labelcolor = background_color)

        # # Set integer tick labels on the x-axis
        # x_ticks = np.linspace(background_mean - (1*background_std + deviation) / 2, background_mean + (1*background_std + deviation), num=3)
        # x_ticks_int = np.round(x_ticks).astype(int)
        # ax.set_xticks(x_ticks_int)

        # Arrow 1
        ax.annotate('', xy=(background_mean, mid_y), xytext=(observed_mean, mid_y), arrowprops=arrow_props)

    else:
        ax.set_xlim(background_mean - (1*background_std + deviation), background_mean + (1*background_std + deviation) / 2)
        # ax.text(background_mean - background_std, mid_y * 1.75, f'Randomized\nmeans', color=background_color, ha='left', va='top')
        # Add legend without border
        #ax.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.2, 1.1), labelcolor = background_color)

        # # Set integer tick labels on the x-axis
        # x_ticks = np.linspace(background_mean - (1*background_std + deviation), background_mean + (1*background_std + deviation) / 2, num=3)
        # x_ticks_int = np.round(x_ticks).astype(int)
        # ax.set_xticks(x_ticks_int)

        # Arrow
        ax.annotate('', xy=(background_mean, mid_y), xytext=(observed_mean, mid_y), arrowprops=arrow_props)



    # Legend
    legend_handles = [Patch(facecolor=background_color_line, alpha=0.2, edgecolor='none', label='Randomized means'),
                        Line2D([0], [0], color="black", linestyle='--', label='Observed mean')]
    legend = ax.legend(handles=legend_handles, frameon=False, loc='upper right', bbox_to_anchor=(1.6, 1.1), fontsize=legend_fontsize)

    # Adjust the color of the text labels in the legend
    for text, color in zip(legend.get_texts(), ["black", "black"]):
        text.set_color(color)


    ax.set_xticks([])

    # Add a vertical line for the observed mean
    ax.axvline(observed_mean, color="black", linestyle='--', ymin=0, ymax=0.5, label="Observed mean")

    # Add a label for the observed mean
    if p_value == 1e-6:
        #ax.text(observed_score*1.1, max(y)/2, text, ha='left', va='center', fontsize=text_fontsize, color=observed_color)
        #ax.text(observed_score*1.1, max(y)/2 - 0.35*(max(y)/2), fr'$\mathit{{p}}$-value < {pvalue}', ha='left', va='center', fontsize=text_fontsize, color=observed_color)
        ax.text(observed_mean * 1.03, mid_y, f'$Score$ = {z_score:.2f}',
                color=observed_color, ha='left', va='center', fontsize = 13)
        ax.text(observed_mean * 1.03, mid_y - 0.35*mid_y, f'$p$-value < {p_value:.2g}',
                color=observed_color, ha='left', va='center', fontsize = 13)
    else:
        ax.text(observed_mean * 1.03, mid_y, f'$Score$ = {z_score:.2f}',
                color=observed_color, ha='left', va='center', fontsize = 13)
        ax.text(observed_mean * 1.03, mid_y - 0.35*mid_y, f'$p$-value = {p_value:.2g}',
                color=observed_color, ha='left', va='center', fontsize = 13)

    # Set labels and title
    ax.set_xlabel('Impact score', fontsize = xaxis_fontsize)

    # Hide the bottom, right, and top borders (spines) of the plot
    for spine in ['left', 'right', 'top']:
        ax.spines[spine].set_visible(False)

    # Hide the entire y-axis
    ax.yaxis.set_visible(False)

    # Display the plot
    plt.show()

#     # Annotate the Z-score, formula and p-value with LaTeX
#     formula = (f'$Z = \\frac{{{observed_mean:.2f} - {background_mean:.2f}}}{{{background_std:.2f}}} = {z_score:.2f}$\n'
#                f'$p$-value = {p_value:.2g}')
#     plt.text(observed_mean - 2, mid_y * 1.25, formula,
#              color='black', ha='center', va='center')

    return fig





def plotting_indels_side(data_gene):

    # Bar labels and values
    labels = ['MULTIPLE_3.non_protein_affecting', 'NOT_MULTIPLE_3.non_protein_affecting',
                'NON_TRUNCATING.protein_affecting', 'TRUNCATING.protein_affecting']
    values = [ data_gene[x] for x in labels ]

    # Colors:
    truncating_indels_color = metrics_colors_dictionary["frameshift"]
    nontruncating_color = metrics_colors_dictionary["inframe"]
    colors = ['none', 'none', nontruncating_color, truncating_indels_color]
    edgecolors = [nontruncating_color, truncating_indels_color, nontruncating_color, truncating_indels_color]

    # Create the plot with a smaller height
    fig, ax = plt.subplots(figsize=(3, 1.75))  # Adjust the height and width

    # Plot the bars horizontally with the specified colors and hatches
    bars = ax.barh(labels, values, color=colors, edgecolor=edgecolors, linewidth = 2)

    # Remove the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


    max_syn = max(values) * 1.55

    # Vertical lines
    line_separation = max(values) * 0.03
    text_of_line_separation = line_separation # / 2
    line_pos1 = max(values[:2]) + line_separation
    line_pos2 = max(values[2:]) + line_separation
    line_margin = 0.04
    ax.axvline(x=line_pos1, ymin = line_margin, ymax = 0.5 - line_margin, color='black',
                # linestyle='--',
                linewidth=1)
    ax.axvline(x=line_pos2, ymin = 0.5 + line_margin, ymax = 1 - line_margin, color='black',
                # linestyle='--',
                linewidth=1)

    # Annotation for "PA" and "NPA" next to vertical lines
    ax.text(line_pos1 + 2.2*text_of_line_separation, 0.5, 'Non\ncoding', ha='center', va='center', fontsize=9,
            rotation = 270
            # bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.5)
            )
    ax.text(line_pos2 + text_of_line_separation, 2.5, 'Coding', ha='center', va='center', fontsize=9,
            rotation = 270
            # bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.5)
            )


    # # Arrows with annotations
    # arrow_props = dict(facecolor='black', arrowstyle='<-')

    # # Arrow 1
    # ax.annotate('', xy=(line_pos1 + 15, 0.5), xytext=(max_syn/1.25, 1.2),
    #             arrowprops=arrow_props)

    # # Arrow 2
    # ax.annotate('', xy=(line_pos2, 2.5), xytext=(max_syn, 1.75),
    #             arrowprops=arrow_props)

    # formula_ratio = (f'$\\frac{{{data_gene["pa_TRUNC/NOTTRUNC"]:.2f}}}{{{data_gene["Npa_NM3/M3"]:.2f}}} = {data_gene["pa/Npa"]:.2f}$')
    # formula_ratio = (f'$\\frac{{{data_gene["pa_TRUNC/NOTTRUNC"]:.2f}}}{{{data_gene["Npa_NM3/M3"]:.2f}}} = {data_gene["pa/Npa"]:.2f}$')
    ax.text(max_syn, 1.5, f'$Score$ = {data_gene["pa/Npa"]:.2f}',
            color=truncating_indels_color, ha='center', va='center', fontsize = 13)

    pvalue_tag = f'$p$-value = {data_gene["pvalue"]:.2g}'
    ax.text(max_syn, 0.8, pvalue_tag,
            color=truncating_indels_color,
            fontsize=12, ha='center', va='center'
            )

    # Add title and labels
    ax.set_xlabel('Count')

    # Customize the ticks and labels
    ax.xaxis.set_ticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # Update labels
    plt.xlabel('Number of indels')
    plt.ylabel('')

    return fig






# ## __Define functions__

def plot_all_positive_selection(omega_truncating,
                                omega_missense,
                                indels_panel_df,
                                oncodrive3d_data_scores,
                                oncodrivefml_data,
                                gene_order,
                                title = None,
                                pvalue_thres = 0.05,
                                linewidth_def = 0.6
                                ):

    num_genes = len(gene_order)

    # Create the figure and subplots
    fig, (ax1, ax2, ax3, ax4, ax5, ) = plt.subplots(5, 1, figsize=(2.5, 2.9), gridspec_kw={'height_ratios': [5, 5, 5, 5, 5]})
    if title:
        fig.suptitle(title)

    # Separate data based on significance of p-value
    omega_truncating_sig = omega_truncating[omega_truncating["pvalue"] <= pvalue_thres].reset_index(drop = True)
    omega_truncating_notsig = omega_truncating[omega_truncating["pvalue"] > pvalue_thres].reset_index(drop = True)

    # Plot the second bar plot in the middle subplot
    sns.barplot(data=omega_truncating_notsig, x='GENE', y='omega_trunc',
                ax=ax1, alpha=1,
                fill = False,
                legend = False,
                linewidth = linewidth_def,
                order = gene_order,
                color = metrics_colors_dictionary["omega_trunc"])


    # Plot the second bar plot in the middle subplot
    sns.barplot(data=omega_truncating_sig, x='GENE', y='omega_trunc',
                ax=ax1, alpha=1,
                legend = False,
                linewidth = linewidth_def,
                order = gene_order,
                color = metrics_colors_dictionary["omega_trunc"],
                edgecolor = None
                )

    ax1.set_xlabel('')
    ax1.set_ylabel('dN/dS of\ntruncating', rotation = 0, labelpad=17, verticalalignment = 'center')
    ax1.set_xticklabels([])  # Hide x-axis labels on the middle plot
    ax1.axhline(1, color='black', linestyle='--')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)




    # Separate data based on significance of p-value
    omega_missense_sig = omega_missense[omega_missense["pvalue"] <= pvalue_thres].reset_index(drop = True)
    omega_missense_notsig = omega_missense[omega_missense["pvalue"] > pvalue_thres].reset_index(drop = True)

    # Plot the second bar plot in the middle subplot
    sns.barplot(data=omega_missense_notsig, x='GENE', y='omega_mis',
                ax=ax2, alpha=1,
                fill = False,
                legend = False,
                linewidth = linewidth_def,
                order = gene_order,
                color = metrics_colors_dictionary["omega_miss"])


    # Plot the second bar plot in the middle subplot
    sns.barplot(data=omega_missense_sig, x='GENE', y='omega_mis',
                ax=ax2, alpha=1,
                legend = False,
                linewidth = linewidth_def,
                order = gene_order,
                color = metrics_colors_dictionary["omega_miss"],
                edgecolor = None
                )

    ax2.set_xlabel('')
    ax2.set_ylabel('dN/dS of\nmissense', rotation = 0, labelpad=17, verticalalignment = 'center')
    ax2.set_xticklabels([])  # Hide x-axis labels on the upper plot
    ax2.axhline(1, color='black', linestyle='--')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)



    df = oncodrive3d_data_scores
    name_metric = "o3d_score"
    variable_name = df.columns[1]
    max_score = df[variable_name].max()
    for j, gene in enumerate(gene_order):
        try:
            value_original = df.loc[df['GENE'] == gene, variable_name].values[0]
            value =  value_original / max_score * 5
            pvalue = df.loc[df['GENE'] == gene, 'pvalue'].values[0]
            color = metrics_colors_dictionary[name_metric] if pvalue < pvalue_thres else 'none'
            edgecolor = metrics_colors_dictionary[name_metric]
            size = value * 12  # Scale size for better visualization
            ax3.scatter(j, 0, s=size, color=color, edgecolors=edgecolor, linewidth = linewidth_def, alpha=0.9,)
        except Exception as e:
            print("Gene", gene, "failed because of", e)
            continue

    # Set axis labels
    ax3.set_xticks(range(len(gene_order)))
    ax3.set_xticklabels([])  # Hide x-axis labels on the middle plot
    ax3.set_yticks([])
    ax3.set_yticklabels([])
    ax3.set_ylabel('3D\nclustering', rotation = 0, labelpad=17,
                    verticalalignment = 'center',
                    horizontalalignment = 'right'
                    )
    ax3.set_ylim(-0.5, 0.5)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Set title and labels
    ax3.set_xlabel('')




    df = oncodrivefml_data
    name_metric = "ofml_score"
    variable_name = df.columns[1]
    max_score = df[variable_name].max()
    for j, gene in enumerate(gene_order):
        try:
            value_original = df.loc[df['GENE'] == gene, variable_name].values[0]
            value =  value_original / max_score * 5
            pvalue = df.loc[df['GENE'] == gene, 'pvalue'].values[0]
            color = metrics_colors_dictionary[name_metric] if pvalue < pvalue_thres else 'none'
            edgecolor = metrics_colors_dictionary[name_metric]
            size = value * 12  # Scale size for better visualization
            if size > 0:
                ax4.scatter(j, 0, s=size, color=color, edgecolors=edgecolor, linewidth = linewidth_def, alpha=0.9)
            else:
                ax4.scatter(j, 0, s=-size, color=color, edgecolors=edgecolor, linewidth = linewidth_def, linestyle = '--', alpha=0.9)
        except Exception as e:
            print("Gene", gene, "failed because of", e)
            continue

    # Set axis labels
    ax4.set_xticks(range(len(gene_order)))
    ax4.set_xticklabels([])
    ax4.set_yticks([])
    ax4.set_yticklabels([])
    ax4.set_ylabel('Functional\nimpact\nbias', rotation = 0, labelpad=17,
                    verticalalignment = 'center',
                    horizontalalignment = 'right'
                    )
    ax4.set_ylim(-0.5, 0.5)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)




    # Separate data based on significance of p-value
    indels_panel_df_sig = indels_panel_df[indels_panel_df["pvalue"] <= pvalue_thres].reset_index(drop = True)
    indels_panel_df_notsig = indels_panel_df[indels_panel_df["pvalue"] > pvalue_thres].reset_index(drop = True)

    # Plot the second bar plot in the middle subplot
    sns.barplot(data=indels_panel_df_notsig, x='GENE', y='Indels_score',
                ax=ax5, alpha=1, #0.6,
                fill = False,
                legend = False,
                linewidth = linewidth_def,
                order = gene_order,
                color = metrics_colors_dictionary["frameshift"])


    # Plot the second bar plot in the middle subplot
    sns.barplot(data=indels_panel_df_sig, x='GENE', y='Indels_score',
                ax=ax5, alpha=1,
                legend = False,
                linewidth = linewidth_def,
                order = gene_order,
                color = metrics_colors_dictionary["frameshift"],
                edgecolor = None
                )

    ax5.set_xlabel('')
    ax5.set_ylabel('Excess of\nframeshift\nindels', rotation = 0,
                    verticalalignment = 'center',
                    horizontalalignment = 'right'
                    )
    ax5.set_xticks(range(num_genes))
    ax5.set_xticklabels(gene_order, rotation=90)
    ax5.axhline(1, color='black', linestyle='--')
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)


    # Set consistent x-axis limits for all subplots
    separation = 2
    ax1.set_xlim([-separation, num_genes - 1 + separation])
    ax2.set_xlim([-separation, num_genes - 1 + separation])
    ax3.set_xlim([-separation, num_genes - 1 + separation])
    ax4.set_xlim([-separation, num_genes - 1 + separation])
    ax5.set_xlim([-separation, num_genes - 1 + separation])

    return fig


def get_all_data(sample, outdir, pvaluee = 0.05):

    oncodrivefml_data = pd.read_table(f"{sample}-oncodrivefml.tsv.gz")
    oncodrivefml_data = oncodrivefml_data[["GENE_ID", "Z-SCORE", "Q_VALUE", "AVG_SCORE_OBS", "POPULATION_MEAN", "STD_OF_MEANS"]]
    oncodrivefml_data.columns = ["GENE", "OncodriveFML", "pvalue", "OBSERVED_MEAN", "BACKGROUND_MEAN", "BACKGROUND_STD"]


    omega_data = pd.read_table(f"output_mle.{sample}.tsv")
    omega_data = omega_data[(omega_data["impact"].isin(['missense', 'truncating']))
                                    & ~(omega_data["gene"].str.contains('--'))
                                ]
    omega_truncating = omega_data[omega_data["impact"] == "truncating"].reset_index(drop = True)[["gene", "dnds", "pvalue", "lower", "upper"]]
    omega_truncating.columns = ["GENE", "omega_trunc", "pvalue", "lower", "upper"]
    truncating_decreasing = list(omega_truncating.sort_values("omega_trunc", ascending= False)["GENE"].values)
    print("Truncating\n", truncating_decreasing)

    omega_missense = omega_data[omega_data["impact"] == "missense"].reset_index(drop = True)[["gene", "dnds", "pvalue", "lower", "upper"]]
    omega_missense.columns = ["GENE", "omega_mis", "pvalue", "lower", "upper"]
    missense_decreasing = list(omega_missense.sort_values("omega_mis", ascending= False)["GENE"].values)
    print("Missense\n", missense_decreasing)

    # merge omegas to decide sorting
    omega_df = omega_truncating[["GENE", "omega_trunc", "pvalue"]].merge(omega_missense[["GENE", "omega_mis", "pvalue"]],
                                                                            on = ["GENE"],
                                                                            suffixes = ("_trunc", "_mis"))
    omega_df["mean_omega"] = omega_df[["omega_trunc", "omega_mis"]].mean(axis = 1)
    omega_df["any_signif"] = omega_df[["pvalue_trunc", "pvalue_mis"]].apply(lambda x: (x < 0.05).any(), axis = 1)
    global_omega_decreasing = list(omega_df.sort_values(by = ["any_signif", "mean_omega"], ascending = False)["GENE"].values)
    global_omega_decreasing.remove("ALL_GENES")

    print("Global\n", global_omega_decreasing)
    if len(global_omega_decreasing) > 20:
        print("Keeping top 20 genes")
        global_omega_decreasing = global_omega_decreasing[:20]

    positively_selected_trunc = omega_truncating[(omega_truncating["pvalue"] < pvaluee) &
                                                    (omega_truncating["omega_trunc"] > 1)
                                                ]["GENE"].values
    positively_selected_mis = omega_missense[(omega_missense["pvalue"] < pvaluee) &
                                                (omega_missense["omega_mis"] > 1)
                                            ]["GENE"].values

    all_positively_selected = set(positively_selected_trunc).union(set(positively_selected_mis))
    print( "all_positively_selected", sorted(all_positively_selected))
    positively_selected_both = set(positively_selected_trunc).intersection(set(positively_selected_mis))
    print( "positively_selected_both", sorted(positively_selected_both))
    positively_selected_trunc_only = set(positively_selected_trunc) - set(positively_selected_mis)
    print( "positively_selected_trunc_only", sorted(positively_selected_trunc_only))
    positively_selected_mis_only = set(positively_selected_mis) - set(positively_selected_trunc)
    print( "positively_selected_mis_only", sorted(positively_selected_mis_only))


    oncodrive3d_data = pd.read_table(f"{sample}.3d_clustering_genes.csv", sep = ',')
    oncodrive3d_data_scores = oncodrive3d_data[["Gene", "Score_obs_sim_top_vol", "qval"]]
    oncodrive3d_data_scores.columns = ["GENE", "Oncodrive3D", 'pvalue']


    indels_data = pd.read_table(f"{sample}.sample.indels.tsv")
    indels_panel_df = indels_data[["SYMBOL", "pa/Npa", "pvalue"]]
    indels_panel_df.columns = ["GENE", "Indels_score", "pvalue"]


    figuree = plot_all_positive_selection(omega_truncating,
                                            omega_missense,
                                            indels_panel_df,
                                            oncodrive3d_data_scores,
                                            oncodrivefml_data,
                                            global_omega_decreasing,
                                            title = sample,
                                            pvalue_thres = pvaluee)


    figuree.savefig(f"{outdir}/{sample}.positive_selection_summary.pdf", bbox_inches='tight')




@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--outdir', type=click.Path(), help='Output path for plots')
def main(sample_name, outdir):
    click.echo("Plotting omega results...")
    generate_all_side_figures(sample_name, outdir)
    get_all_data(sample_name, outdir)

if __name__ == '__main__':
    main()





