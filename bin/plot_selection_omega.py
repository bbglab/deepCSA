#!/usr/bin/env python


import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress,norm
# import tabix
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
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




def generate_all_side_figures(sample,
                                mut_file,
                                omega_file,
                                # gene_list = None,
                                tools = ["omega_trunc", "omega_mis"]
                                ):

    snvs_maf = pd.read_table(mut_file, na_values = custom_na_values)
    snvs_maf = snvs_maf[snvs_maf["TYPE"] == "SNV"].reset_index(drop = True)

    possible_genes = []

    omega_data = pd.read_table(omega_file)
    omega_data = omega_data[omega_data["impact"].isin(['missense', 'truncating'])]
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


    gene_list = list(set(possible_genes))

    gene_list = [x for x in gene_list if x in list(pd.unique(snvs_maf["canonical_SYMBOL"]))]

    os.makedirs(f"{sample}.plots")

    for genee in gene_list:
        print(genee)
        try :
            if "omega_trunc" in tools:
                if genee in omega_truncating_genes and genee in omega_missense_genes:
                    omega_df = build_counts_from_df_complete(genee, snvs_maf, omega_truncating, omega_missense)

                    fig_gene_omega = plot_omega_side_complete(omega_df)
                    fig_gene_omega.savefig(f"{sample}.plots/{genee}.{sample}.omega_side.pdf", bbox_inches='tight')
                    plt.show()
                    plt.close()

        except Exception as exe:
            print(genee)
            print(exe)







def build_counts_from_df_complete(genee, snvs_maf, omega_truncating, omega_missense):

    trunc_omega = float(omega_truncating[omega_truncating["GENE"] == genee]["omega_trunc"].values[0])
    trunc_pvalue = float(omega_truncating[omega_truncating["GENE"] == genee]["pvalue"].values[0])

    miss_omega = float(omega_missense[omega_missense["GENE"] == genee]["omega_mis"].values[0])
    miss_pvalue = float(omega_missense[omega_missense["GENE"] == genee]["pvalue"].values[0])
    snvs_gene = snvs_maf[snvs_maf["canonical_SYMBOL"] == genee].reset_index(drop = True)


    # Calculate counts based on canonical consequences
    truncating_count = float(omega_truncating[omega_truncating["GENE"] == genee]["mutations_trunc"].values[0])
    missense_count = float(omega_missense[omega_missense["GENE"] == genee]["mutations_mis"].values[0])
    synonymous_count = snvs_gene[snvs_gene["canonical_Consequence_broader"].isin(["synonymous"])].shape[0]

    # Compute
    expected_missense = (1 - ((miss_omega - 1) / miss_omega)) * missense_count
    expected_truncating = (1 - ((trunc_omega - 1) / trunc_omega)) * truncating_count


    # Create a dataframe from the counts and expected values
    data = {
        'type': ['truncating', 'synonymous', 'missense'],
        'number_obs': [truncating_count, synonymous_count, missense_count],
        'expected': [expected_truncating, None, expected_missense],
        'omega': [trunc_omega, None, miss_omega],
        'pvalue': [trunc_pvalue, None, miss_pvalue]
    }
    df = pd.DataFrame(data)

    # Print the final dataframe
    return df






def plot_omega_side_complete(df):

    consequence_order = ['truncating', "synonymous", "missense"]

    # Define colors for each type
    colors = {
        'truncating': colors_dictionary["omega_trunc"],
        'missense': colors_dictionary["omega_miss"],
        'synonymous': colors_dictionary["omega_synon"]
    }

    t_omega = df[df['type'] == 'truncating']['omega'].item()
    t_pvalue = df[df['type'] == 'truncating']['pvalue'].item()
    m_omega = df[df['type'] == 'missense']['omega'].item()
    m_pvalue = df[df['type'] == 'missense']['pvalue'].item()


    # Plot for truncating, missense, and synonymous
    fig, ax = plt.subplots(figsize=(4, 2))

#     # Bar positions
#     y_positions = np.arange(len(consequence_order)) * 2  # Add separation between bars

#     # Plot the observed values as filled bars
#     for i, cons in enumerate(consequence_order):
#         obs_value = df[df['type'] == cons]['number_obs'].values[0]
#         ax.barh(y_positions[i], obs_value, color=colors[cons], label=cons)

    sns.barplot(data=df, y='type', x='number_obs',
                hue = 'type', legend = False,
                hue_order = consequence_order,
                order=consequence_order, palette=[colors[x] for x in consequence_order],
                ax = ax
                )

    # Overlay the expected values with shaded portions
    leg = True
    for i, row in df.iterrows():
        if row['expected'] is not None:
            # plt.barh(row['type'], row['expected'], color='gray', alpha=0.5, edgecolor='none')
            # plt.barh(row['type'], row['expected'], color='none', edgecolor=colors[row['type']], hatch='//', linewidth=0)
            if leg:
                ax.barh(row['type'], row['expected'], color='none', edgecolor="gray", hatch='////', linewidth=0, label = 'expected')
                leg = False
            else:
                ax.barh(row['type'], row['expected'], color='none', edgecolor="gray", hatch='////', linewidth=0)



    # Remove the right and top spines
    # ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Customize the ticks and labels
    ax.xaxis.set_ticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])


    # Text annotations
    text_separation = max(df['number_obs']) * 0.03
    deviation_for_centering = text_separation * 12

    # text_of_line_separation = line_separation / 3
    text_pos1 = df[df['type'].isin(consequence_order[:2])][['number_obs', 'expected']].max().max() + text_separation + deviation_for_centering
    text_pos2 = df[df['type'].isin(consequence_order[1:])][['number_obs', 'expected']].max().max() + text_separation + deviation_for_centering
    # print(text_pos1, text_pos2)

    ax.text(text_pos1, 0.05, r'$\omega_{trunc}$ = ' + f"{t_omega:.2f}",
            # transform=ax.transAxes,
            fontsize=13, ha='center', va='bottom',
            color=colors['truncating'])
                   #f'$p$-value = {p_value:.2g}')

    ax.text(text_pos1, 0.45, f'$p$-value < {max(t_pvalue, 1e-6):.1e}',
            # transform=ax.transAxes,
            fontsize=12, ha='center', va='bottom',
            color=colors['truncating'])

    ax.text(text_pos2, 1.95, r'$\omega_{mis}$ = ' + f"{m_omega:.2f}",
            # transform=ax.transAxes,
            fontsize=13, ha='center', va='bottom',
            color=colors['missense'])
    ax.text(text_pos2, 2.35, f'$p$-value < {max(m_pvalue, 1e-6):.1e}',
            # transform=ax.transAxes,
            fontsize=12, ha='center', va='bottom',
            color=colors['missense'])

    # Update labels
    plt.xlabel('Number of Mutations')
    # plt.xscale('log')
    plt.ylabel('')
    plt.legend(loc='right', frameon=False)
    # plt.legend(loc='right', frameon=False, fontsize=12, handlelength=3, markerscale=2)  # Increase legend text and marker size
    # plt.show()
    return fig




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

    sample_name_    = sys.argv[1]
    mut_filee        = sys.argv[2]
    omega_filee      = sys.argv[3]
    # o3d_seq_file_ = sys.argv[3]
    # sample_name_out_ = sys.argv[4]
    # out_maf      = sys.argv[3]
    # json_filters = sys.argv[4]
    # req_plots    = sys.argv[5]

    generate_all_side_figures(sample_name_,
                                mut_filee,
                                omega_filee)
