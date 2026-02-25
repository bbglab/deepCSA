#!/usr/bin/env python

"""
Plot depth of all genes and specific genes per sample groups.

This script plots the average depth at all consensus exons across all genes and for specific genes of interest, either in the panel or in a custom subset of genes, stratified by sample groups defined in the metadata.
The output is stored in a pdf file.
"""

import click
import pandas as pd
from utils_plot import plots_general_config
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

mpl.rcParams.update({
    'axes.titlesize': plots_general_config["title_fontsize"],
    'axes.labelsize': plots_general_config["xylabel_fontsize"],
    'xtick.labelsize': plots_general_config["xyticks_fontsize"],
    'ytick.labelsize': plots_general_config["xyticks_fontsize"],
    'figure.titlesize': plots_general_config["title_fontsize"],
})

separator2character = {
    'tab' : '\t',
    'comma' : ','
}


def plot_depth_per_group(df, group_col, data_type, pdf):
    '''
    Function to plot depth within a group of samples in all genes or a specific subset of genes
    '''
 
    col_name = group_col[0] if isinstance(group_col, list) else group_col

    # Get the number of unique categories to plot
    num_categories = df[col_name].nunique() + 2
    plt.figure(figsize=(num_categories, 4))

    ax = sns.boxplot(data=df, x=col_name, y="MEAN_GENE_DEPTH", hue=col_name, showfliers=False, showmeans=False,legend=False)
    ax = sns.stripplot(data=df, x=col_name, y="MEAN_GENE_DEPTH", color='grey', alpha=0.5, size=4, legend=False)    

    if data_type == 'all_genes':
        plt.title(f"Average Depth for {data_type} in {col_name} group", fontsize=plots_general_config["title_fontsize"])
    elif data_type == 'gene':
        gene = df['GENE'].iloc[0]
        plt.title(f"Average Depth for {gene} in {col_name} group", fontsize=plots_general_config["title_fontsize"])
    else:
        print(f"Unknown data type: {data_type}. Title will not be set.")
    
    plt.xlabel('', fontsize=plots_general_config["xylabel_fontsize"])
    plt.ylabel(f"Average Cons Exons Depth", fontsize=plots_general_config["xylabel_fontsize"])
    plt.yticks( fontsize=plots_general_config["yticks_fontsize"])
    plt.xticks(fontsize=plots_general_config["xticks_fontsize"])
    plt.tick_params(axis='x', rotation=90)
    plt.tight_layout()
    pdf.savefig()
    plt.close()
    plt.show()

    return 


@click.command()
@click.option('--table-filename', required=True, type=click.Path(exists=True), help='Input features table file')
@click.option('--depth-table', required=True, type=click.Path(exists=True), help='Input depth table file')
@click.option('--separator', required=True, type=click.Choice(['tab', 'comma']), help='Separator used in features table: tab or comma')
@click.option('--unique-identifier', default=None, type=str, help='Unique identifier column name')
@click.option('--groups', default=None, type=str, help='List of columns with grouping information')
@click.option('--custom-genes', required=False, type=str, help='Comma separated list of custom genes')
@click.option('--output_prefix', type=str, required=True, help='Prefix for output files')


def main(table_filename, depth_table, unique_identifier, separator, groups, custom_genes, output_prefix):

    sep_char = separator2character[separator]
    
    # Read tables
    features_table = pd.read_table(table_filename, header=0, sep=sep_char)
    depth_table = pd.read_table(depth_table, header=0, sep="\t")

    # Process panel genes 
    panel_genes = sorted(set(depth_table['GENE'].unique()))
    if custom_genes:
        print(f'Custom genes provided, plotting custom genes only: {custom_genes}')
        custom_gene_list = [g.strip() for g in custom_genes.split(",")]
        panel_genes = sorted(set(custom_gene_list) & set(panel_genes))
    
    else:
        print(f'No custom genes provided, plotting all genes in the panel: {panel_genes}')

    output_name = f"{output_prefix}.plot_depth_per_group.pdf"

    # groups may contain lists of lists, but all formatted into a string
    groups_of_interest_init = [group.strip().strip(",").split(",") for group in groups.replace("[", ";;;").replace("]", "").split(";;;")] if groups else []

    groups_of_interest = []
    for comparison in groups_of_interest_init:
        comparison_group_clean = [item.strip() for item in comparison]
        comparison_group = [item for item in comparison_group_clean if item != '']
        if len(comparison_group) > 0:
            groups_of_interest.append(comparison_group)

    uniq_name = unique_identifier if unique_identifier else "sample"

    print(f"Processing data for the groups of interest: {groups_of_interest}")

    with PdfPages(output_name) as pdf:
        for group in groups_of_interest:
            print(f"Processing {group} group, type: {type(group)}")
            metadata_group_df = features_table[[uniq_name, str(group[0])]]
            merged_depth_df = pd.merge(metadata_group_df, depth_table, how='left', left_on=uniq_name, right_on='SAMPLE_ID')
            print(merged_depth_df[[uniq_name, str(group[0]), 'MEAN_GENE_DEPTH']].head())
            print('Length of the processed table', len(merged_depth_df))

            # Plot depth of all samples for each group
            plot_depth_per_group(merged_depth_df, group, 'all_genes', pdf)

            # Plot depths per gene (defined cutom genes or genes in the panel) for each group
            merged_depth_df = merged_depth_df[merged_depth_df['GENE'].isin(panel_genes)]
            for gene in panel_genes:
                gene_data = merged_depth_df[merged_depth_df['GENE'] == gene]
                print('Length of gene data for gene', gene, ':', len(gene_data))
                plot_depth_per_group(gene_data, group, 'gene', pdf)

            print(f"Plots saved as {output_name}")

if __name__ == "__main__":
    main()

'''
Example usage:
python depth_group_comparison.py \
       --table-filename metadata_table_all_with_bacterial_signatures.tsv \
        --depth-table all_samples.exons_cons.depth_per_gene_per_sample.tsv \
        --separator tab \
        --unique-identifier Sample_Name \
        --groups "[ ["Sample_Group"], ["cancer"], ["Age_onset"], ["Cancer_age_group"] , ["Bacterial_Signatures_identified"]]" \
        --custom-genes APC,BRAF,FBXW7,KRAS,PIK3CA,SMAD4,TP53' \
        --output_prefix depth_group_comparison
'''