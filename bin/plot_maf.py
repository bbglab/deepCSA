#!/usr/bin/env python


# import click
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from utils import filter_maf
from read_utils import custom_na_values


def subset_mutation_dataframe(mutations_file, json_filters):
    """
    INFO
    """
    # Load your MAF DataFrame (raw_annotated_maf)
    raw_annotated_maf = pd.read_csv(mutations_file, sep = "\t", header = 0, na_values = custom_na_values)

    data_tuples = []

    # Load the filter criteria from the JSON file
    with open(json_filters, 'r') as file:
        filter_criteria = json.load(file)

        # Convert the dictionary into a list of tuples
        for key, value in filter_criteria.items():
            if isinstance(value, list):
                for item in value:
                    data_tuples.append((key, item))
            else:
                data_tuples.append((key, value))
        print(data_tuples)


    if len(data_tuples) > 0:
        # Filter the annotated maf using the described filters
        print("MAF subset")
        return filter_maf(raw_annotated_maf, data_tuples)

    return raw_annotated_maf


def plot_mutations_per_sample(sample_name, maf, sample_column_name = "SAMPLE_ID"):
    muts_per_sample = maf.groupby(sample_column_name).size().to_frame("NUM_MUTS").reset_index()

    # Calculate the length of the longest sample name
    max_label_length = muts_per_sample[sample_column_name].astype(str).str.len().max()

    # Determine the rotation angle and adjust figure size accordingly
    rotation_angle = 90 if max_label_length > 7 else 30  # Adjust threshold as needed
    # fig_height = 4 + (max_label_length / 10) * 2  # Adjust multiplier as needed
    # # Calculate the figure width based on the number of samples
    # fig_width = min(18, max(2, len(muts_per_sample) * 0.5))

    # # Create the figure and axis
    # fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    fig, ax = plt.subplots(figsize=(18, 4))
    sns.barplot(data = muts_per_sample, x = sample_column_name, y = "NUM_MUTS",
                ax = ax, palette = ["salmon"], dodge = False)
    plt.xticks(rotation = rotation_angle)
    plt.xlabel("")
    plt.ylabel("Number of mutations", fontsize = 14)
    ax.set_title(f"Somatic mutations in {sample_name}")
    plt.tight_layout()

    return fig


def plot_mutations_per_gene(sample_name, maf, parameters = {}):
    default_vals = {'gene_column_name' : 'SYMBOL', 'minimum_muts' : 5, 'annotation' : True}
    default_vals.update(parameters)

    muts_per_gene = maf.groupby(default_vals['gene_column_name']).size().to_frame("NUM_MUTS")
    muts_per_gene = muts_per_gene[muts_per_gene["NUM_MUTS"] >= default_vals['minimum_muts']].reset_index()

    # # Calculate the length of the longest sample name
    # max_label_length = muts_per_gene[default_vals['gene_column_name']].str.len().max()

    # # Determine the rotation angle and adjust figure size accordingly
    # rotation_angle = 30 if max_label_length > 10 else 0  # Adjust threshold as needed
    # fig_height = 4 + (max_label_length / 10) * 2  # Adjust multiplier as needed
    # # Calculate the figure width based on the number of samples
    # fig_width = min(18, max(2, len(muts_per_gene) * 0.5))

    # # Create the figure and axis
    # fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig, ax = plt.subplots(figsize=(18, 4))

    muts_per_gene = muts_per_gene.sort_values(by = "NUM_MUTS", ascending = False).reset_index(drop = True)

    # Plot the barplot
    sns.barplot(data=muts_per_gene,
                x = default_vals['gene_column_name'], y="NUM_MUTS",
                ax=ax,
                color="salmon"
                )
    if default_vals['annotation']:
        # Calculate the maximum height of the bars
        max_height = max(bar.get_height() for bar in ax.patches)

        # Set the spacing above the bars for annotations
        annotation_spacing = max_height * 0.05  # Adjust this value as needed

        # Increase y-axis limits to accommodate annotations
        y_min, y_max = ax.get_ylim()
        ax.set_ylim(y_min, y_max + annotation_spacing)

        # Add annotations above each bar
        for bar in ax.patches:
            height = bar.get_height()
            ax.annotate(format(height, '.0f'),
                        (bar.get_x() + bar.get_width() / 2, height),
                        ha='center', va='bottom', fontsize=12, xytext=(0, 3),
                        textcoords='offset points')

    ax.set_xlabel("Genes", fontsize=14)
    ax.set_ylabel("Number mutations", fontsize=14)

    plt.xticks(fontsize=10, rotation=90)
    plt.yticks(fontsize=12)
    plt.tick_params(axis='x', which='major', labelsize=14)
    ax.set_title(f"Mutations per gene in {sample_name}")
    plt.tight_layout()

    return fig


def plot_filter_stats(df, x_axis_group, hue_group, logy=False, stacked=True, fig_width=20, fig_height=5, max_groups_per_plot=30):
    """
    Only admits two variables to combine, ideally
    the x_axis_group variable should be sample or
    gene symbol
    """

    # Get the unique values of x_axis_group
    unique_x_values = sorted(df[x_axis_group].unique())

    # Check if the number of unique values exceeds the maximum
    num_unique_x_values = len(unique_x_values)
    if num_unique_x_values > max_groups_per_plot:
        # Calculate the number of plots needed
        num_plots = int(np.ceil(num_unique_x_values / max_groups_per_plot))

        # Split the unique x values into chunks for each plot
        x_chunks = np.array_split(unique_x_values, num_plots)

        # Iterate over the x chunks and create separate plots
        fig_list = []
        for i, x_chunk in enumerate(x_chunks):
            # Create a new figure for each plot
            fig, ax = plt.subplots(1, 1)
            fig.set_size_inches(fig_width, fig_height)
            plt.xticks(rotation=90)
            plt.title(f"Plotting mutations per {x_axis_group} colored by {hue_group} (Plot {i+1}/{num_plots})")

            # Filter the dataframe to include only the rows with x values in the current chunk
            subset_df = df[df[x_axis_group].isin(x_chunk)]

            # Plot the data for the current chunk
            if stacked:
                gp_df = subset_df.groupby([x_axis_group, hue_group], dropna=False).size().to_frame('number_mutations').reset_index().pivot(
                    columns=hue_group, index=x_axis_group, values="number_mutations")
                gp_df.plot(kind='bar', stacked=True, ax=ax)
            else:
                gp_df = subset_df.groupby([x_axis_group, hue_group], dropna=False).size().to_frame('number_mutations').reset_index()
                gp_df = gp_df.fillna("FALSE")
                sns.barplot(data=gp_df, x=x_axis_group, y="number_mutations", hue=hue_group, ax=ax)
                if logy:
                    ax.set_yscale('log')

            plt.tight_layout()
            fig_list.append(fig)
            plt.close()

        return fig_list

    else:
        # Plot all x values in a single plot
        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(fig_width, fig_height)
        plt.xticks(rotation=90)
        plt.title(f"Plotting mutations per {x_axis_group} colored by {hue_group}")

        if stacked:
            gp_df = df.groupby([x_axis_group, hue_group], dropna=False).size().to_frame('number_mutations').reset_index().pivot(
                columns=hue_group, index=x_axis_group, values="number_mutations")
            gp_df.plot(kind='bar', stacked=True, ax=ax)
        else:
            gp_df = df.groupby([x_axis_group, hue_group], dropna=False).size().to_frame('number_mutations').reset_index()
            gp_df = gp_df.fillna("FALSE")
            sns.barplot(data=gp_df, x=x_axis_group, y="number_mutations", hue=hue_group, ax=ax)
            if logy:
                ax.set_yscale('log')

        plt.tight_layout()

        return [fig]




def filter_wrapper(sample_name, maf, parameters = {}):
    var_to_plot = parameters.get('variable', "SAMPLE_ID")
    filters_to_plot = parameters.get('filter_list', ['no_pileup_support', 'n_rich'])

    fig_list = []

    for filt in filters_to_plot:
        if f"FILTER.{filt}" not in maf.columns:
            continue
        fig = plot_filter_stats(maf, var_to_plot, f"FILTER.{filt}")
        if type(fig) == list:
            for subfig in fig:
                fig_list.append(subfig)
        else:
            fig_list.append(fig)

    return fig_list



def variable_plot_wrapper(sample_name, maf, parameters = {}):
    var_to_plot = parameters.get('variable', "SAMPLE_ID")
    vars_to_plot = parameters.get('columns_list', ['canonical_Consequence_broader', 'TYPE'])

    fig_list = []

    for varp in vars_to_plot:
        fig = plot_filter_stats(maf, var_to_plot, f"{varp}")
        if type(fig) == list:
            for subfig in fig:
                fig_list.append(subfig)
        else:
            fig_list.append(fig)

    return fig_list



dict_plotname2func = {
    "per_gene" : plot_mutations_per_gene,
    "per_sample" : plot_mutations_per_sample,
    "filter_stats" : filter_wrapper,
    "plot_stats" : variable_plot_wrapper
}

def plot_manager(sample_name, maf, plotting_criteria_file):
    # Load the filter criteria from the JSON file
    with open(plotting_criteria_file, 'r') as file:
        plotting_criteria = json.load(file)

    for suffix, criteria_list in plotting_criteria.items():
        suffix = suffix.strip('.')

        with PdfPages(f'{sample_name}.{suffix}.pdf') as pdf:
            for criterion in criteria_list:
                if criterion.startswith("filter_stats"):
                    main_criterion = criterion.split(' ')[0]
                    gene_or_sample = criterion.split(' ')[1]
                    filters = criterion.split(' ')[2].split(",")
                    fig_lisst = dict_plotname2func[main_criterion](sample_name, maf, parameters = {'variable' : gene_or_sample, 'filter_list' : filters})
                    for ff in fig_lisst: pdf.savefig(ff)

                elif criterion.startswith("plot_stats"):
                    main_criterion = criterion.split(' ')[0]
                    gene_or_sample = criterion.split(' ')[1]
                    cols = criterion.split(' ')[2].split(",")
                    fig_lisst = dict_plotname2func[main_criterion](sample_name, maf, parameters = {'variable' : gene_or_sample, 'columns_list' : cols})
                    for ff in fig_lisst: pdf.savefig(ff)

                elif ' ' in criterion:
                    main_criterion = criterion.split(' ')[0]
                    additional_params = { x.split(":")[0] : x.split(":")[1] for x in criterion.split(' ')[1:] }
                    fig1 = dict_plotname2func[main_criterion](sample_name, maf, parameters = additional_params)
                    pdf.savefig()
                    plt.close()


                else :
                    fig1 = dict_plotname2func[criterion](sample_name, maf)
                    pdf.savefig()
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


sample_name  = sys.argv[1]
mut_file     = sys.argv[2]
out_maf      = sys.argv[3]
json_filters = sys.argv[4]
req_plots    = sys.argv[5]



if __name__ == '__main__':
    maf = subset_mutation_dataframe(mut_file, json_filters)
    plot_manager(sample_name, maf, req_plots)
