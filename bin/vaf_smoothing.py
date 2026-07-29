#!/usr/bin/env python


import click
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from read_utils import custom_na_values
from utils_plot import plots_general_config
import matplotlib as mpl




def compute_priors(mutdensity_file, samples):
    """
    From simple adjusted version
    This function needs to be removed eventually when the use of the
    updated adjusted mutation density is functional in omega
    """

    mutdensity_df = pd.read_csv(mutdensity_file, sep = "\t", header = 0, na_values = custom_na_values)
    npa_mutdensities = mutdensity_df[(mutdensity_df["MUTTYPES"] == "all_types")
                                     & (mutdensity_df["GENE"] == "ALL_GENES")
                                     & (mutdensity_df["SAMPLE_ID"].isin(samples))
                                     ].reset_index(drop = True)
    npa_mutdensities['prior'] = npa_mutdensities['N_MUTATED'] / npa_mutdensities['DEPTH']
    npa_mutdensities_values = npa_mutdensities[['SAMPLE_ID', 'prior']]

    npa_mutdensities_values.to_csv(f"sample_specific_priors.tsv",
                                   header=True,
                                   index=False,
                                   sep="\t")



def apply_correction_per_sample(mutations_table, priors_per_sample, depths_per_sample,
                                samples,
                                weights_list = [0, 0.2, 0.5, 0.75, 1, 1.25, 1.5, 2, 3]
                                ):
    """
    Apply the pseudocount correction to the VAF values in the mutations table.
    """

    mutation_tables_versions = []

    for sample_id in samples:
        prior_vaf = priors_per_sample.loc[priors_per_sample['SAMPLE_ID'] == sample_id, 'prior'].values[0]
        sample_mutations = mutations_table[mutations_table['SAMPLE_ID'] == sample_id][
            ['SAMPLE_ID', 'MUT_ID', 'ALT_DEPTH', 'DEPTH', 'ALT_DEPTH_AM', 'DEPTH_AM']
        ]
        sample_depths = depths_per_sample[depths_per_sample['SAMPLE_ID'] == sample_id]
        average_depth = sample_depths['DEPTH']

        for weight_prop in weights_list:
            weight = weight_prop * average_depth // 1
            sample_mutations[f'VAF_PSEUDO_{weight_prop}'] = vaf_pseudocount(
                                        sample_mutations['ALT_DEPTH'],
                                        sample_mutations['DEPTH'],
                                        weight,
                                        prior_vaf=prior_vaf
                                    )
            sample_mutations[f'VAF_AM_PSEUDO_{weight_prop}'] = vaf_pseudocount(
                                                    sample_mutations['ALT_DEPTH_AM'],
                                                    sample_mutations['DEPTH_AM'],
                                                    weight,
                                                    prior_vaf=prior_vaf
                                                )
            
            mutation_tables_versions.append(sample_mutations[['SAMPLE_ID', 'MUT_ID'] + [x for x in sample_mutations.columns if x.startswith('VAF_PSEUDO') or x.startswith('VAF_AM_PSEUDO')]])
    
    return pd.concat(mutation_tables_versions, ignore_index=True)




def vaf_pseudocount(alt_depth, depth, weight, prior_vaf=None):    
    return (alt_depth + prior_vaf * weight) / (depth + weight)

def plot_vaf_pseudocount_curve(maf_df, output_pdf, suffix=''):
    """
    Plot VAF distribution compared to VAF_AM in a histogram.
    
    Parameters:
    -----------
    maf_df : DataFrame
        MAF dataframe containing VAF and VAF_AM columns
    """

    fig, axes = plt.subplots(3, 3, figsize=(8, 8))
    axes = axes.flatten()
    dg = maf_df
    fig.suptitle(f'VAF{suffix} with pseudocounts')
    for i, weight_prop in enumerate([0, 0.2, 0.5, 0.75, 1, 1.25, 1.5, 2, 3]):
        axes[i].scatter(dg[f'DEPTH{suffix}'], dg[f'VAF{suffix}_PSEUDO'], s=3, alpha=0.1)
        rho = np.corrcoef(np.log(dg[f'DEPTH{suffix}']), np.log(dg[f'VAF{suffix}_PSEUDO']))[0, 1]
        axes[i].set_xlabel(f'DEPTH{suffix}')
        axes[i].set_ylabel(f'VAF{suffix}_PSEUDO', fontsize=6)
        axes[i].set_yscale('log')
        axes[i].set_xscale('log')
        axes[i].set_title(f"{weight_prop}\nrho={rho:.2f}, \nprop_mut_tissue={dg[f'VAF{suffix}_PSEUDO'].sum():.2f}", fontsize=6)
    plt.tight_layout()
    output_pdf.savefig()
    plt.close()
    plt.show()
    



@click.command()
@click.option('--mutdensities', type=click.Path(exists=True), help='Input mutation density file')
@click.option('--mutations', type=click.Path(), help='Mutations file')
@click.option('--depth-sample', type=click.Path(), help='Depth per sample file')

def main(mutdensities, mutations, depth_sample):
    """
    Main function to execute the VAF smoothing pipeline.
    """
    click.echo("Selecting the gene synonymous mutation densities...")
    mutations_table = pd.read_csv(mutations, sep = "\t", header = 0, na_values = custom_na_values)
    depths_per_sample = pd.read_csv(depth_sample, sep = "\t", header = 0, na_values = custom_na_values)

    samples_list = mutations_table['SAMPLE_ID'].unique()

    priors_per_sample = compute_priors(mutdensities, samples_list)
    all_mutations_table = apply_correction_per_sample(mutations_table, priors_per_sample, depths_per_sample, samples_list)

    # with PdfPages("vaf_pseudocounts.pdf") as pdf:
    #     plot_vaf_pseudocount_curve(all_mutations_table, priors_per_sample, pdf, suffix='')

if __name__ == '__main__':
    main()

