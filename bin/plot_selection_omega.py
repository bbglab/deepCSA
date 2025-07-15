#!/usr/bin/env python


import click
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
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


def generate_all_side_figures(sample,
                                mut_file,
                                omega_file,
                                plots_output_dir
                                ):

    maf = pd.read_table(mut_file, na_values = custom_na_values)
    snvs_maf = maf[maf["TYPE"] == "SNV"].reset_index(drop = True)

    possible_genes = []

    omega_data = pd.read_table(omega_file)
    omega_data = omega_data[omega_data["impact"].isin(['missense', 'truncating'])]

    omega_truncating = omega_data[omega_data["impact"] == "truncating"].reset_index(drop = True)[["gene", "mutations", "dnds", "pvalue", "lower", "upper"]]
    omega_truncating.columns = ["GENE", "mutations_trunc", "omega_trunc", "pvalue", "lower", "upper"]
    omega_truncating_genes = list(pd.unique(omega_truncating["GENE"]))
    possible_genes += omega_truncating_genes

    omega_missense = omega_data[omega_data["impact"] == "missense"].reset_index(drop = True)[["gene", "mutations", "dnds", "pvalue", "lower", "upper"]]
    omega_missense.columns = ["GENE", "mutations_mis", "omega_mis", "pvalue", "lower", "upper"]
    omega_missense_genes = list(pd.unique(omega_truncating["GENE"]))
    possible_genes += omega_missense_genes


    gene_list = list(set(possible_genes).intersection(set(snvs_maf["canonical_SYMBOL"].unique())))


    for genee in gene_list:
        print(genee)
        try :
            if genee in omega_truncating_genes and genee in omega_missense_genes:
                omega_df = build_counts_from_df_complete(genee, snvs_maf, omega_truncating, omega_missense)

                fig_gene_omega = plot_omega_vertical(omega_df)
                fig_gene_omega.savefig(f"{plots_output_dir}/{genee}.{sample}.omega.pdf", bbox_inches='tight', dpi = 100)
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
    print(df)

    # Print the final dataframe
    return df






def plot_omega_vertical(df,
                        ymax = None,
                        bar_width=0.8,
                        figsize=(1.4, 1.17),
                        between_text = 1.5,
                        withinbartext_off = 1.8,
                        text_off = 0.5,
                        min_pvalue = 1e-6,
                        gene = None,
                        legenddd = True
                        ):
    consequence_order = ['truncating', 'missense', 'synonymous',]

    # Define colors
    colors = {
        'truncating': metrics_colors_dictionary["omega_trunc"],
        'missense': metrics_colors_dictionary["omega_miss"],
        'synonymous': metrics_colors_dictionary["omega_synon"]
    }

    # Filter relevant data
    df = df[df['type'].isin(consequence_order)]

    t_obs = df[df['type'] == 'truncating']['number_obs'].item()
    t_omega = df[df['type'] == 'truncating']['omega'].item()
    t_pvalue = df[df['type'] == 'truncating']['pvalue'].item()

    m_obs = df[df['type'] == 'missense']['number_obs'].item()
    m_omega = df[df['type'] == 'missense']['omega'].item()
    m_pvalue = df[df['type'] == 'missense']['pvalue'].item()

    s_obs = df[df['type'] == 'synonymous']['number_obs'].item()  # Added synonymous mutations

    # Compute x positions for bars
    spacing_factor = bar_width * 1.1  # Adjust spacing based on bar width
    x_positions = np.arange(len(consequence_order)) * spacing_factor

    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # **Matplotlib Barplot**
    ax.bar(x_positions,
            [t_obs, m_obs, s_obs],
            color=[colors[x] for x in consequence_order],
            width=bar_width,
            edgecolor='none')

    # Overlay expected values as hatched bars (only for truncating & missense)
    for i, row in df.iterrows():
        if row['type'] != 'synonymous':  # No hatch for synonymous
            if legenddd:
                ax.bar(x_positions[consequence_order.index(row['type'])], row['expected'],
                        color='none', edgecolor="grey", hatch= '//////',
                        linewidth=0,
                        width=bar_width,
                        label = 'expected'
                        )
                legenddd = False
            else:
                ax.bar(x_positions[consequence_order.index(row['type'])], row['expected'],
                        color='none', edgecolor="grey", hatch='//////',
                        linewidth=0,
                        width=bar_width
                        )

    # Remove top/right spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Customize ticks and labels
    ax.set_xticks(x_positions)
    ax.set_xticklabels([])
    # ax.set_yticklabels(ax.get_yticklabels())

    plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))


    # ax.spines['left'].set_visible(False)
    ax.set_ylabel('Number of mutations')

    # Positioning text annotations
    between_text_offset = max(max(df['number_obs'].max(), df['expected'].max()) * 0.1, between_text)
    text_offset = max(max(df['number_obs'].max(), df['expected'].max()) * 0.02, text_off)
    within_bar_text_offset = max(max(df['number_obs'].max(), df['expected'].max()) * 0.1, withinbartext_off)

    for i, row in df.iterrows():
        x_pos = x_positions[consequence_order.index(row['type'])]
        y_pos = max(row['number_obs'], row['expected']) + text_offset
        y_pos_low = max(row['number_obs'], row['expected']) - within_bar_text_offset
        omega_value = t_omega if row['type'] == 'truncating' else (m_omega if row['type'] == 'missense' else None)
        p_value = t_pvalue if row['type'] == 'truncating' else (m_pvalue if row['type'] == 'missense' else None)

        # Omega annotation (above the bar) - Only for truncating/missense
        if omega_value is not None:
            excess_mutss = row["number_obs"]*((omega_value-1)/omega_value)
            ax.text(x_pos, y_pos + between_text_offset,
                    rf'dN/dS={omega_value:.2f}',
                    fontsize=plots_general_config["annots_fontsize"], ha='center', va='bottom',
                    color='black'
                    )

            # P-value annotation (below omega)
            ax.text(x_pos, y_pos,
                    f'$p$<{min_pvalue:.1e}' if p_value < min_pvalue else (f'$p$={p_value:.1e}' if p_value < 0.01 else f'$p$={p_value:.2f}'),
                    fontsize=plots_general_config["annots_fontsize"], ha='center', va='bottom',
                    color='black'
                    )

            # Add excess mutations in bar
            if excess_mutss >= 1:
                ax.text(x_pos, y_pos_low,
                        f'{excess_mutss:,.0f}',
                        fontsize=plots_general_config["annots_fontsize"], ha='center', va='bottom', color= 'black')

        else:
            mutations = row['number_obs']
            ax.text(x_pos,
                    y_pos,
                    rf'{mutations:.0f}',
                    fontsize=plots_general_config["annots_fontsize"], ha='center', va='bottom', color='gray')


    plt.legend(frameon=False, bbox_to_anchor = (1,1))

    if ymax is not None:
        plt.ylim(0,ymax)

    if gene is not None:
        plt.title(gene, pad = 12)

    return fig





@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--mut_file', type=click.Path(exists=True), help='Input mutations file')
@click.option('--omega_file', type=click.Path(exists=True), help='Input omega results file')
@click.option('--outdir', type=click.Path(), help='Output path for plots')
def main(sample_name, mut_file, omega_file, outdir):
    click.echo("Plotting omega results...")
    generate_all_side_figures(sample_name, mut_file, omega_file, outdir)

if __name__ == '__main__':
    main()
