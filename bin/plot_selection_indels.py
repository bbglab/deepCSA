#!/usr/bin/env python


import sys
import os
import pandas as pd
# import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import linregress,norm
# # import tabix
# import matplotlib.cm as cm
# import matplotlib.colors as mcolors
# from matplotlib.patches import Patch
# from matplotlib.lines import Line2D
from read_utils import custom_na_values


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
# sample_name_out_ = sys.argv[4]
# out_maf      = sys.argv[3]
# json_filters = sys.argv[4]
# req_plots    = sys.argv[5]



def generate_all_side_figures(sample,
                                gene_list = None):

    possible_genes = []
    indels_data = pd.read_table(f"{sample}.sample.indels.tsv",
                                    sep = '\t',
                                    header = 0, na_values = custom_na_values)
    indels_genes = list(pd.unique(indels_data["SYMBOL"]))
    possible_genes += indels_genes

    gene_list = list(set(possible_genes))


    os.makedirs(f"{sample}.plots")

    for genee in gene_list:
        print(genee)
        try :

            indel_data_gene = indels_data[indels_data["SYMBOL"] == genee].to_dict(orient='records')[0]
            fig_gene_indel = plotting_indels_side(indel_data_gene)
            fig_gene_indel.savefig(f"{sample}.plots/{genee}.{sample}.indels_side.pdf", bbox_inches='tight')
            plt.show()
            plt.close()

        except Exception as exe:
            print(genee)
            print(exe)







def plotting_indels_side(data_gene):

    # Bar labels and values
    labels = ['MULTIPLE_3.non_protein_affecting', 'NOT_MULTIPLE_3.non_protein_affecting',
                'NON_TRUNCATING.protein_affecting', 'TRUNCATING.protein_affecting']
    values = [ data_gene[x] for x in labels ]

    # Colors:
    truncating_indels_color = colors_dictionary["frameshift"]
    nontruncating_color = colors_dictionary["inframe"]
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






# generate_all_side_figures("all_samples", deepcsa_run_dir,
#                           # store_plots_dir = "Fig3/plots/2024-07-08_side"
#                          )



colors_dictionary = {"ofml"        : "viridis_r",
                     "ofml_score"  : "#6A33E0",
                     "omega_trunc" : "#FA5E32",
                     "omega_synon" : "#89E4A2",
                     "omega_miss"  : "#FABE4A",
                     "o3d_score"   : "#6DBDCC",
                     "o3d_cluster" : "skyblue",
                     "o3d_prob"    : "darkgray",
                     "frameshift"  : "#E4ACF4",
                     "inframe"     : "C5",
                     "hv_lines"    : "lightgray" # horizontal and vertical lines
                    }



if __name__ == '__main__':
    # maf = subset_mutation_dataframe(mut_file, json_filters)
    # plot_manager(sample_name, maf, req_plots)
    # manager(mut_file, o3d_seq_file_, sample_name_, sample_name_out_)
    generate_all_side_figures(sample_name_)
