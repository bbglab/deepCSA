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

def mut_density_heatmaps(data, genes_list, samples_list, outdir,
                         config_datasets = {"SNVs" : ({"MUTTYPES": 'SNV', "REGIONS": 'all'}, 'MUTDENSITY_MB'),
                                            "SNVs protein-affecting" : ({"MUTTYPES": 'SNV', "REGIONS": 'protein_affecting'}, 'MUTDENSITY_MB'),
                                            "SNVs non-protein-affecting" : ({"MUTTYPES": 'SNV', "REGIONS": 'non_protein_affecting'}, 'MUTDENSITY_MB')
                                            }):
    """
    Take a DataFrame and create mutation density heatmaps, storing all plots in a single PDF file.
    Allows control over figure size and sample/gene ordering.
    SAMPLE_ID       GENE    REGIONS MUTTYPES        DEPTH   N_MUTS  N_MUTATED       MUTDENSITY_MB   MUTDENSITY_MB_ADJUSTED MUTREADSRATE_MB  MUTREADSRATE_MB_ADJUSTED
    """

    pdf_filename = f"{outdir}/mut_density_heatmaps.pdf"
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


def plotting_manager(outdir, samples_list, genes_list, data,
                     requested_plots = ["mut_density_heatmaps"]):
    dict_plotname2func = {
        "mut_density_heatmaps" : mut_density_heatmaps,
        # "adjusted_mut_density_heatmaps" : adj_mut_density_heatmaps,
        # "per_sample" : plot_mutations_per_sample,
        # "filter_stats" : filter_wrapper,
        # "plot_stats" : variable_plot_wrapper
    }
    
    for plot in requested_plots:
        if plot in dict_plotname2func:
            dict_plotname2func[plot](data, samples_list, genes_list, outdir)
        else:
            print("requested plot not available")

def loading_data(panel_regions, mutdensities, samples_json):
    panel_regions_df = pd.read_csv(panel_regions, sep="\t")
    panel_genes = sorted(panel_regions_df['GENE'].unique().tolist())
    del panel_regions_df

    # this should be imported from utils.py
    # since it is something that might be of interest in other scenarios
    with open(samples_json, 'r') as f:
        samples_definition = json.load(f)
    all_samples_names = list(samples_definition.keys())
    
    mut_dens = pd.read_csv(mutdensities, sep="\t")
    mut_dens = mut_dens[mut_dens["SAMPLE_ID"].isin(all_samples_names)]

    return panel_genes, all_samples_names, mut_dens

@click.command()
@click.option('--outdir', type=click.Path(), help='Output path for plots')
@click.option('--panel-regions', type=click.Path(), help='Panel regions file')
@click.option('--mutdensities', type=click.Path(), help='Mutation densities file')
@click.option('--samples_json', type=click.Path(), help='Samples JSON file')
def main(outdir, panel_regions, mutdensities, samples_json):
    click.echo("Exploring variability...")
    genes_list, samples_list, mutation_densities = loading_data(panel_regions, mutdensities, samples_json)

    try:
        plotting_manager(outdir, genes_list, samples_list, mutation_densities)
    except Exception as e:
        print("Error in the process", e)

if __name__ == '__main__':
    main()

