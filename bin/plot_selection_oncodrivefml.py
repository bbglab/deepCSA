#!/opt/conda/bin/python


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


sample_name_  = sys.argv[1]
mut_file     = sys.argv[2]
o3d_seq_file_ = sys.argv[3]
# sample_name_out_ = sys.argv[4]
# out_maf      = sys.argv[3]
# json_filters = sys.argv[4]
# req_plots    = sys.argv[5]



def generate_all_side_figures(sample,
                                gene_list = None,
                                tools = ["oncodrivefml", "oncodrive3d", "omega_trunc", "omega_mis", "excess_indels"]
                                ):

    snvs_maf = pd.read_table(f"{sample}.somatic.mutations.tsv", na_values = custom_na_values)

    possible_genes = []
    if "oncodrivefml" in tools:
        oncodrivefml_data = pd.read_table(f"{sample}-oncodrivefml.tsv.gz")
        oncodrivefml_data = oncodrivefml_data[["GENE_ID", "Z-SCORE", "Q_VALUE", "AVG_SCORE_OBS", "POPULATION_MEAN", "STD_OF_MEANS"]]
        oncodrivefml_data.columns = ["GENE", "OncodriveFML", "pvalue", "OBSERVED_MEAN", "BACKGROUND_MEAN", "BACKGROUND_STD"]
        oncodrivefml_genes = list(pd.unique(oncodrivefml_data["GENE"]))
        possible_genes += oncodrivefml_genes


    if "omega_trunc" in tools or "omega_mis" in tools:
        omega_data = pd.read_table(f"output_mle.{sample}.tsv")
        omega_data = omega_data[omega_data["impact"].isin(['missense', 'truncating'])]
        if "omega_trunc" in tools :
            omega_truncating = omega_data[omega_data["impact"] == "truncating"].reset_index(drop = True)[["gene", "dnds", "pvalue", "lower", "upper"]]
            omega_truncating.columns = ["GENE", "omega_trunc", "pvalue", "lower", "upper"]
            omega_truncating_genes = list(pd.unique(omega_truncating["GENE"]))
            possible_genes += omega_truncating_genes

        if "omega_mis" in tools :
            omega_missense = omega_data[omega_data["impact"] == "missense"].reset_index(drop = True)[["gene", "dnds", "pvalue", "lower", "upper"]]
            omega_missense.columns = ["GENE", "omega_mis", "pvalue", "lower", "upper"]
            omega_missense_genes = list(pd.unique(omega_truncating["GENE"]))
            possible_genes += omega_missense_genes


    if "excess_indels" in tools:
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

            if "oncodrivefml" in tools:
                # there is no run of oncodrivefml with ALL_GENES
                if genee in oncodrivefml_genes:
                    oncodrivefml_gene_data = oncodrivefml_data[oncodrivefml_data["GENE"] == genee].to_dict(orient='records')[0]

                    fig_gene_fml = plot_oncodrivefml_side(oncodrivefml_gene_data)

                    fig_gene_fml.savefig(f"{sample}.plots/{genee}.{sample}.oncodrivefml_side.pdf", bbox_inches='tight')
                    plt.show()
                    plt.close()

            if "omega_trunc" in tools:
                if genee in omega_truncating_genes and genee in omega_missense_genes:
                    omega_df = build_counts_from_df_complete(genee, snvs_maf, omega_truncating, omega_missense)

                    fig_gene_omega = plot_omega_side_complete(omega_df)
                    fig_gene_omega.savefig(f"{sample}.plots/{genee}.{sample}.omega_side.pdf", bbox_inches='tight')
                    plt.show()
                    plt.close()


            if "excess_indels" in tools:
                if genee in indels_genes:
                    indel_data_gene = indels_data[indels_data["SYMBOL"] == genee].to_dict(orient='records')[0]
                    fig_gene_indel = plotting_indels_side(indel_data_gene)
                    fig_gene_indel.savefig(f"{sample}.plots/{genee}.{sample}.indels_side.pdf", bbox_inches='tight')
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
    snvs_gene = snvs_maf[snvs_maf["SYMBOL"] == genee].reset_index(drop = True)


    # Calculate counts based on canonical consequences
    truncating_count = snvs_gene[snvs_gene["canonical_Consequence_broader"].isin(['nonsense', "essential_splice"])].shape[0]
    missense_count = snvs_gene[snvs_gene["canonical_Consequence_broader"].isin(['missense'])].shape[0]
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

    # Plot the observed values as filled bars
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
    text_pos1 = max(df[df['type'].isin(consequence_order[:2])]['number_obs']) + text_separation + deviation_for_centering
    text_pos2 = max(df[df['type'].isin(consequence_order[1:])]['number_obs']) + text_separation + deviation_for_centering
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





def plot_oncodrivefml_side(geneee_data):
    legend_fontsize = 12
    xaxis_fontsize = 12

    # Extract the necessary values
    observed_mean = geneee_data['OBSERVED_MEAN']
    background_mean = geneee_data['BACKGROUND_MEAN']
    background_std = geneee_data['BACKGROUND_STD']
    p_value = geneee_data['pvalue']

    observed_color = colors_dictionary["ofml_score"]

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
