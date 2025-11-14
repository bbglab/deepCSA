#!/usr/bin/env python

import click
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


from matplotlib.backends.backend_pdf import PdfPages



def filter_data_from_config(dataa, config):
    filtered_data = dataa.copy()
    for key, value in config.items():
        # print("Filtering:", key, value)
        # print(filtered_data.shape)
        filtered_data = filtered_data[filtered_data[key] == value]
        # print(filtered_data.shape)
    return filtered_data

def mut_density_heatmaps(data, genes_list, samples_list, outdir, prefix = '',
                            config_datasets = {
                                "all" : ({"MUTTYPES": 'all_types', "REGIONS": 'all'}, 'MUTDENSITY_MB'),
                                "all protein-affecting" : ({"MUTTYPES": 'all_types', "REGIONS": 'protein_affecting'}, 'MUTDENSITY_MB'),
                                "all non-protein-affecting" : ({"MUTTYPES": 'all_types', "REGIONS": 'non_protein_affecting'}, 'MUTDENSITY_MB'),
                                "SNVs" : ({"MUTTYPES": 'SNV', "REGIONS": 'all'}, 'MUTDENSITY_MB'),
                                "SNVs protein-affecting" : ({"MUTTYPES": 'SNV', "REGIONS": 'protein_affecting'}, 'MUTDENSITY_MB'),
                                "SNVs non-protein-affecting" : ({"MUTTYPES": 'SNV', "REGIONS": 'non_protein_affecting'}, 'MUTDENSITY_MB'),
                                "INDELs" : ({"MUTTYPES": 'DELETION-INSERTION', "REGIONS": 'all'}, 'MUTDENSITY_MB'),
                                "INDELs protein-affecting" : ({"MUTTYPES": 'DELETION-INSERTION', "REGIONS": 'protein_affecting'}, 'MUTDENSITY_MB'),
                                "INDELs non-protein-affecting" : ({"MUTTYPES": 'DELETION-INSERTION', "REGIONS": 'non_protein_affecting'}, 'MUTDENSITY_MB')
                            }
                        ):
    """
    Take a DataFrame and create mutation density heatmaps, storing all plots in a single PDF file.
    Allows control over figure size and sample/gene ordering.
    SAMPLE_ID       GENE    REGIONS MUTTYPES        DEPTH   N_MUTS  N_MUTATED       MUTDENSITY_MB   MUTDENSITY_MB_ADJUSTED MUTREADSRATE_MB  MUTREADSRATE_MB_ADJUSTED
    """
    data = data[data["SAMPLE_ID"].isin(samples_list)].reset_index(drop = True)

    if data.shape[0] == 0:
        print("No data available for the selected samples/groups")
        return
    
    pdf_filename = f"{outdir}/{prefix}mut_density_heatmaps.pdf"
    with PdfPages(pdf_filename) as pdf:

        for title, (config, value) in config_datasets.items():
            print("Creating heatmap for:", title, config, value)
            filtered_data = filter_data_from_config(data, config)
            # print(filtered_data[['GENE', 'SAMPLE_ID', value]].head())
            # Create a pivot table for the heatmap
            heatmap_data = filtered_data.pivot_table(index='GENE', columns='SAMPLE_ID', values=value)
            heatmap_data = heatmap_data.reindex(index=genes_list, columns=samples_list)
            heatmap_data = heatmap_data.replace([np.inf, -np.inf], np.nan).fillna(0).astype(float)

            heatmap_data.to_csv(f"{outdir}/{title.replace(' ', '')}_{value}.tsv", sep='\t', index=True)

            # Standard heatmap
            plt.figure(figsize=(max(10, 0.1*len(samples_list)), max(8, 0.05*len(genes_list))))
            sns.heatmap(heatmap_data, cmap='viridis', cbar_kws={'label': value})
            plt.title(f"{title}")
            plt.xlabel('Sample')
            plt.ylabel('Gene')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            pdf.savefig()
            plt.close()

            # Clustermap
            g = sns.clustermap(heatmap_data, cmap='viridis', figsize=(max(10, 0.1*len(samples_list)), max(8, 0.05*len(genes_list))))
            plt.suptitle(f"{title} - Clustermap")
            pdf.savefig(g.fig)
            plt.close(g.fig)



def adj_mut_density_heatmaps(data, genes_list, samples_list, outdir, prefix = '',
                                config_datasets = {
                                    "all" : ({"MUTTYPES": 'all_types', "REGIONS": 'all'}, 'MUTDENSITY_MB'),
                                    "all protein-affecting" : ({"MUTTYPES": 'all_types', "REGIONS": 'protein_affecting'}, 'MUTDENSITY_MB'),
                                    "all non-protein-affecting" : ({"MUTTYPES": 'all_types', "REGIONS": 'non_protein_affecting'}, 'MUTDENSITY_MB'),
                                    "SNVs" : ({"MUTTYPES": 'SNV', "REGIONS": 'all'}, 'MUTDENSITY_MB'),
                                    "SNVs protein-affecting" : ({"MUTTYPES": 'SNV', "REGIONS": 'protein_affecting'}, 'MUTDENSITY_MB'),
                                    "SNVs non-protein-affecting" : ({"MUTTYPES": 'SNV', "REGIONS": 'non_protein_affecting'}, 'MUTDENSITY_MB'),
                                    "INDELs" : ({"MUTTYPES": 'DELETION-INSERTION', "REGIONS": 'all'}, 'MUTDENSITY_MB'),
                                    "INDELs protein-affecting" : ({"MUTTYPES": 'DELETION-INSERTION', "REGIONS": 'protein_affecting'}, 'MUTDENSITY_MB'),
                                    "INDELs non-protein-affecting" : ({"MUTTYPES": 'DELETION-INSERTION', "REGIONS": 'non_protein_affecting'}, 'MUTDENSITY_MB')
                                }
                            ):
    """
    Take a DataFrame and create mutation density heatmaps, storing all plots in a single PDF file.
    Allows control over figure size and sample/gene ordering.
    """
    data = data[data["SAMPLE_ID"].isin(samples_list)].reset_index(drop = True)

    if data.shape[0] == 0:
        print("No data available for the selected samples/groups")
        return
    
    pdf_filename = f"{outdir}/{prefix}mut_density_heatmaps.pdf"
    with PdfPages(pdf_filename) as pdf:

        for title, (config, value) in config_datasets.items():
            print("Creating heatmap for:", title, config, value)
            filtered_data = filter_data_from_config(data, config)
            # print(filtered_data[['GENE', 'SAMPLE_ID', value]].head())
            # Create a pivot table for the heatmap
            heatmap_data = filtered_data.pivot_table(index='GENE', columns='SAMPLE_ID', values=value)
            heatmap_data = heatmap_data.reindex(index=genes_list, columns=samples_list)
            heatmap_data = heatmap_data.replace([np.inf, -np.inf], np.nan).fillna(0).astype(float)

            heatmap_data.to_csv(f"{outdir}/{title.replace(' ', '')}_{value}.tsv", sep='\t', index=True)

            # Standard heatmap
            plt.figure(figsize=(max(10, 0.1*len(samples_list)), max(8, 0.05*len(genes_list))))
            sns.heatmap(heatmap_data, cmap='viridis', cbar_kws={'label': value})
            plt.title(f"{title}")
            plt.xlabel('Sample')
            plt.ylabel('Gene')
            plt.xticks(rotation=45)
            plt.tight_layout()
            pdf.savefig()
            plt.close()

            # Clustermap
            g = sns.clustermap(heatmap_data, cmap='viridis', figsize=(max(10, 0.1*len(samples_list)), max(8, 0.05*len(genes_list))))
            plt.suptitle(f"{title} - Clustermap")
            pdf.savefig(g.fig)
            plt.close(g.fig)


def plotting_manager(outdir, samples_list, genes_list, out_prefix,
                     requested_plots = [], data = [],
                     ):
    dict_plotname2func = {
        "mutdensities" : mut_density_heatmaps,
        "adjmutdensities" : adj_mut_density_heatmaps,

        # "adjusted_mut_density_heatmaps" : adj_mut_density_heatmaps,
        # "per_sample" : plot_mutations_per_sample,
        # "filter_stats" : filter_wrapper,
        # "plot_stats" : variable_plot_wrapper
    }
    
    for index, plot in enumerate(requested_plots):
        if plot in dict_plotname2func:
            dict_plotname2func[plot](data[index], samples_list, genes_list, outdir, out_prefix)
        else:
            print("requested plot not available:", plot)



def loading_data(panel_regions, mutdensities, adjusted_mutdensities, samples_json, all_groups_json):
    loaded_data = []
    loaded_data_names = []

    panel_regions_df = pd.read_csv(panel_regions, sep="\t")
    panel_genes = sorted(panel_regions_df['GENE'].unique().tolist())
    del panel_regions_df

    # this should be imported from utils.py
    # since it is something that might be of interest in other scenarios
    with open(samples_json, 'r') as f:
        samples_definition = json.load(f)
    all_samples_names = list(samples_definition.keys())
    
    # this should be imported from utils.py
    # since it is something that might be of interest in other scenarios
    with open(all_groups_json, 'r') as f:
        groups_definition = json.load(f)
    all_groups_names = list(groups_definition.keys())

    if mutdensities:
        mut_dens = pd.read_csv(mutdensities, sep="\t")
        loaded_data.append(mut_dens)
        loaded_data_names.append('mutdensities')

    if adjusted_mutdensities:
        adj_mut_dens = pd.read_csv(adjusted_mutdensities, sep="\t")
        loaded_data.append(adj_mut_dens)
        loaded_data_names.append('adjmutdensities')


    return panel_genes, all_samples_names, all_groups_names, loaded_data_names, loaded_data


@click.command()
@click.option('--outdir', type=click.Path(), help='Output path for plots')
@click.option('--panel-regions', type=click.Path(), help='Panel regions file')
@click.option('--mutdensities', default = None, type=click.Path(), help='Mutation densities file')
@click.option('--adjusted-mutdensities', type=click.Path(), help='Mutation densities file')
@click.option('--samples-json', type=click.Path(), help='Samples JSON file')
@click.option('--all-groups-json', type=click.Path(), help='All groups JSON file')
def main(outdir, panel_regions, samples_json, all_groups_json, mutdensities, adjusted_mutdensities):
    click.echo("Exploring variability...")
    genes_list, samples_list, all_groups_list, data_string, data_objects = loading_data(panel_regions, mutdensities, adjusted_mutdensities, samples_json, all_groups_json)
    groups_names = [ x for x in all_groups_list if x not in samples_list ]

    try:
        plotting_manager(outdir, genes_list, samples_list, "samples.", data_string, data_objects)
    except Exception as e:
        print("Error in the process", e)
    
    try:
        plotting_manager(outdir, genes_list, groups_names, "groups.", data_string, data_objects)
    except Exception as e:
        print("Error in the process", e)


if __name__ == '__main__':
    main()

