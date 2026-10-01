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
import ast
import re

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

    if data_type == 'ALL_GENES':
        plot_df = df[df['GENE'] == 'ALL_GENES']
        title = f"Average Depth for ALL_GENES in {col_name} group"

    else: # data_type is the Gene Name
        plot_df = df[df['GENE'] == data_type]
        title = f"Average Depth for {data_type} in {col_name} group"
    
    if plot_df.empty:
        print(f"No data available for {data_type} in {col_name} group. Skipping plot.")
        return
    
    sns.boxplot(data=plot_df, x=col_name, y="MEAN_GENE_DEPTH", hue=col_name, showfliers=False, showmeans=False, legend=False)
    sns.stripplot(data=plot_df, x=col_name, y="MEAN_GENE_DEPTH", color='grey', alpha=0.5, size=4, legend=False)    

    plt.title(title, fontsize=plots_general_config["title_fontsize"])                
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
@click.option('--depth-gene-sample', required=True, type=click.Path(exists=True), help='Input depth file per gene per sample')
@click.option('--depth-sample', required=True, type=click.Path(exists=True), help='Input depth file per sample')
@click.option('--separator', required=True, type=click.Choice(['tab', 'comma']), help='Separator used in features table: tab or comma')
@click.option('--unique-identifier', default=None, type=str, help='Unique identifier column name')
@click.option('--groups', required=True, type=str, help='List of columns with grouping information')
@click.option('--custom-genes', required=False, type=str, help='Comma separated list of custom genes')


def main(table_filename, depth_gene_sample, depth_sample, unique_identifier, separator, groups, custom_genes):

    sep_char = separator2character[separator]
    uniq_name = unique_identifier if unique_identifier else "sample"
    
    # Read tables
    features_table = pd.read_table(table_filename, header=0, sep=sep_char)
    depth_genes_samples = pd.read_table(depth_gene_sample, header=0, sep="\t")
    depth_per_sample = pd.read_table(depth_sample, header=0, sep="\t")

    # Process panel genes 
    panel_genes = sorted(set(depth_genes_samples['GENE'].unique()))
    if custom_genes:
        print(f'Custom genes provided, plotting custom genes only: {custom_genes}')
        custom_gene_list = [g.strip() for g in custom_genes.split(",")]
        panel_genes = sorted(set(custom_gene_list) & set(panel_genes))
    
    else:
        print(f'No custom genes provided, plotting all genes in the panel: {panel_genes}')


    # Process depth per sample to add the 'ALL_GENES' depth value per sample
    print('Processing per sample depth table to add the ALL_GENES depth column: ')
    depth_per_sample['GENE'] = 'ALL_GENES'
    depth_per_sample = depth_per_sample.rename(columns={'avg_depth_sample': 'MEAN_GENE_DEPTH'})
    
    print('Depth per sample table after adding ALL_GENES value in GENES column:')
    print(depth_per_sample.head())
    
    depth_genes_samples = pd.concat([depth_genes_samples, depth_per_sample], ignore_index=True)
    print('Added ALL_GENES depth values to depth genes samples table:')
    print('Length of table:', len(depth_genes_samples))
    print(depth_genes_samples.head())

    # Merge data with metadata table to get the group information for each sample
    merged_depth_df = pd.merge(features_table, depth_genes_samples, how='left', left_on=uniq_name, right_on='SAMPLE_ID')
    print('Merged depth and metadata table:')
    print(merged_depth_df.head())
    print('Length of merged table:', len(merged_depth_df))

    # Process groups
    # First clean the string (adds quotes if the shell stripped them)
    cleaned = re.sub(r'(?<!["\'])\b([a-zA-Z0-9_ ]*)\b(?!["\'])', r'"\1"', groups)

    # Then parse into a list of lists
    raw_groups = ast.literal_eval(cleaned) if groups else []

    # Finally flatten and get unique items
    # This nested loop goes through each sub-list and each item within it
    groups_of_interest = []
    for sublist in raw_groups:
        if type(sublist) == list:
            for item in sublist:
                groups_of_interest.append(item)
        else:
            groups_of_interest.append(sublist)
    groups_of_interest = list(set(groups_of_interest))


    print(f"Processing data for the groups of interest: {groups_of_interest}")

    # Handle groups so each group has its own plot in all and individual genes and stored in the same pdf file per group
    for group in groups_of_interest:
        group_formatted = group.replace(" ", "_")
        output_name = f"{group_formatted}.plot_depth_per_group.pdf"
        
        with PdfPages(output_name) as pdf:    
            print(f"Processing {group} group, type: {type(group)}")
            group_table = merged_depth_df[[uniq_name, group, 'GENE', 'MEAN_GENE_DEPTH']]
            print(group_table.head())
            print('Length of the processed table', len(group_table))

            # Plot depth of all samples for each group
            plot_depth_per_group(group_table, group, 'ALL_GENES', pdf)

            # Plot depths per gene (defined cutom genes or genes in the panel) for each group
            gene_table = group_table[group_table['GENE'].isin(panel_genes)]
            for gene in panel_genes:
                gene_table = group_table[group_table['GENE'] == gene]
                print('Length of gene data for gene', gene, ':', len(gene_table))
                plot_depth_per_group(gene_table, group, gene, pdf)

            print(f"Plots saved as {output_name}")

if __name__ == "__main__":
    main()

'''
Example usage:
python depth_group_comparison.py \
       --table-filename metadata_table_all_with_bacterial_signatures.tsv \
        --depth-table all_samples.exons_cons.depth_per_gene_per_sample.tsv \
        --depth-sample all_samples.exons_cons.depth_per_sample.tsv \
        --separator tab \
        --unique-identifier Sample_Name \
        --groups "[ ["Sample_Group"], ["cancer"], ["Age_onset"], ["Cancer_age_group"] , ["Bacterial_Signatures_identified"]]" \
        --custom-genes APC,BRAF,FBXW7,KRAS,PIK3CA,SMAD4,TP53' \
'''