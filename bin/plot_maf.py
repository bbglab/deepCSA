#!/usr/local/bin/python


# import click
import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from utils import filter_maf


def subset_mutation_dataframe(mutations_file, json_filters):
    """
    INFO
    """
    # Load your MAF DataFrame (raw_annotated_maf)
    raw_annotated_maf = pd.read_csv(mutations_file, sep = "\t", header = 0)


    # Load the filter criteria from the JSON file
    with open(json_filters, 'r') as file:
        filter_criteria = json.load(file)

    if len(filter_criteria) > 0:
        # Filter the annotated maf using the described filters
        print("MAF subset")
        return filter_maf(raw_annotated_maf, filter_criteria)

    return raw_annotated_maf


def plot_mutations_per_sample(sample_name, maf, sample_column_name = "SAMPLE_ID"):
    muts_per_sample = maf.groupby(sample_column_name).size().to_frame("NUM_MUTS").reset_index()

    # Calculate the length of the longest sample name
    max_label_length = muts_per_sample[sample_column_name].str.len().max()

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

    muts_per_gene = muts_per_gene.sort_values(by = "NUM_MUTS").reset_index(drop = True)

    # Plot the barplot
    sns.barplot(data=muts_per_gene,
                x = default_vals['gene_column_name'], y="NUM_MUTS",
                ax=ax,
                color="salmon"
                #, order=gene_order
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

    plt.xticks(fontsize=10, rotation=30)
    plt.yticks(fontsize=12)
    plt.tick_params(axis='x', which='major', labelsize=14)
    ax.set_title(f"Mutations per gene in {sample_name}")
    plt.tight_layout()

    return fig


dict_plotname2func = {
    "per_gene" : plot_mutations_per_gene,
    "per_sample" : plot_mutations_per_sample
}

def plot_manager(sample_name, maf, plotting_criteria_file):
    # Load the filter criteria from the JSON file
    with open(plotting_criteria_file, 'r') as file:
        plotting_criteria = json.load(file)

    for suffix, criteria_list in plotting_criteria.items():
        suffix = suffix.strip('.')

        with PdfPages(f'{sample_name}.{suffix}.pdf') as pdf:
            for criterion in criteria_list:
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
