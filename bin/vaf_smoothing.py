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
                                     & (mutdensity_df["REGIONS"] == "non_protein_affecting")
                                     & (mutdensity_df["SAMPLE_ID"].isin(samples))
                                     ].reset_index(drop = True)
    npa_mutdensities['prior'] = npa_mutdensities['N_MUTATED'] / npa_mutdensities['DEPTH']
    npa_mutdensities_values = npa_mutdensities[['SAMPLE_ID', 'prior']]

    npa_mutdensities_values.to_csv(f"sample_specific_priors.tsv",
                                   header=True,
                                   index=False,
                                   sep="\t")
    return npa_mutdensities_values



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
            ['SAMPLE_ID', 'MUT_ID', 'ALT_DEPTH', 'DEPTH', 'VAF', 'ALT_DEPTH_AM', 'DEPTH_AM', 'VAF_AM']
        ]
        sample_depths = depths_per_sample[depths_per_sample['SAMPLE_ID'] == sample_id]
        average_depth = sample_depths['avg_depth_sample'].values[0]

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

            mutation_tables_versions.append(sample_mutations)
    
    return pd.concat(mutation_tables_versions, ignore_index=True)


def vaf_pseudocount(alt_depth, depth, weight, prior_vaf=None):    
    return (alt_depth + prior_vaf * weight) / (depth + weight)

def plot_vaf_pseudocount_curve(maf_df, samples, prior_vafs, output_pdf, suffix=''):
    """
    Plot VAF distribution compared to VAF_AM in a histogram.
    
    Parameters:
    -----------
    maf_df : DataFrame
        MAF dataframe containing VAF and VAF_AM columns
    """

    
    for sample in samples:
        sample_prior_vaf = prior_vafs.loc[prior_vafs['SAMPLE_ID'] == sample, 'prior'].values[0]
        fig, axes = plt.subplots(3, 3, figsize=(8, 8), sharex=True, sharey=True)
        axes = axes.flatten()
        fig.suptitle(f'{sample} : VAF{suffix} with pseudocounts')

        samples_maf_df = maf_df[maf_df['SAMPLE_ID'] == sample]
        for i, weight_prop in enumerate([0, 0.2, 0.5, 0.75, 1, 1.25, 1.5, 2, 3]):
            axes[i].scatter(samples_maf_df[f'DEPTH{suffix}'], samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}'], s=3, alpha=0.1)
            rho = np.corrcoef(np.log(samples_maf_df[f'DEPTH{suffix}']), np.log(samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}']))[0, 1]
            # axes[i].axhline(y=sample_prior_vaf, color='r', linestyle='--')
            axes[i].set_xlabel(f'DEPTH{suffix}')
            axes[i].set_ylabel(f'VAF{suffix}_PSEUDO', fontsize=6)
            axes[i].set_yscale('log')
            axes[i].set_xscale('log')
            axes[i].set_title(f"{weight_prop}\nrho={rho:.2f}, \nprop_mut_tissue={samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}'].sum():.2f}", fontsize=6)
        plt.tight_layout()
        output_pdf.savefig()
        plt.close()
        plt.show()




def calc_vaf_distance_summary(
    df,
    weights,
    vaf_col="VAF",
    target_col="VAF_AM",
    pseudo_prefix="VAF_PSEUDO",
    depth_col="DEPTH",
    depth_quantile=0.1,
    large_clone_threshold=1e-2,
):
    """
    Calculate distance between VAF_PSEUDO and target VAF column
    for low-depth clones and large clones.
    """

    depth_cutoff = df[depth_col].quantile(depth_quantile)

    df_low_depth = df[(df[depth_col] < depth_cutoff) & (df['ALT_DEPTH'] < 2)].copy()
    df_large = df[df[target_col] > large_clone_threshold].copy()

    # Reference: original VAF vs target VAF
    ref_dist_log_low_depth = np.abs(np.log10(df_low_depth[vaf_col]) - np.log10(df_low_depth[target_col])).sum()
    ref_dist_log_large = np.abs(np.log10(df_large[vaf_col]) - np.log10(df_large[target_col])).sum()
    ref_dist_low_depth = np.abs(df_low_depth[vaf_col] - df_low_depth[target_col]).sum()
    ref_dist_large = np.abs(df_large[vaf_col] - df_large[target_col]).sum()

    results = []

    for w in weights:

        col = f"{pseudo_prefix}_{w}"

        dist_log_low_depth = np.abs(np.log10(df_low_depth[col]) -  np.log10(df_low_depth[target_col]))
        dist_log_large = np.abs(np.log10(df_large[col]) -  np.log10(df_large[target_col]))
        dist_low_depth = np.abs(df_low_depth[col] -  df_low_depth[target_col])
        dist_large = np.abs(df_large[col] -  df_large[target_col])

        results.append({
            "w": w,
            "sum_dist_log_low_depth": dist_log_low_depth.sum(),
            "sum_dist_low_depth": dist_low_depth.sum(),
            "mean_dist_log_low_depth": dist_log_low_depth.mean(),
            "sum_dist_log_large_clones": dist_log_large.sum(),
            "sum_dist_large_clones": dist_large.sum(),
            "mean_dist_log_large_clones": dist_log_large.mean(),
            "n_low_depth": len(df_low_depth),
            "n_large_clones": len(df_large),
            "depth_cutoff": depth_cutoff,
            "ref_dist_log_low_depth": ref_dist_log_low_depth,
            "ref_dist_log_large_clones": ref_dist_log_large,
            "ref_dist_low_depth": ref_dist_low_depth,
            "ref_dist_large_clones": ref_dist_large,
        })

    return pd.DataFrame(results)



def plot_vaf_distance_summary(
    dist_df,
    use_mean=False,
    figsize=(7, 4),
    log_dist=False,
    log_y=False,
):
    """
    Plot pseudocount distance curves and reference VAF vs VAF_AM distances.
    """

    if use_mean:
        y_low = "mean_dist_log_low_depth"
        y_large = "mean_dist_large_clones_log"
        ylabel = "Mean log-distance"
    else:
        if log_dist:
            ref_low = dist_df["ref_dist_log_low_depth"].iloc[0]
            ref_large = dist_df["ref_dist_log_large_clones"].iloc[0]
            y_low = "sum_dist_log_low_depth"
            y_large = "sum_dist_log_large_clones"
            ylabel = "Sum log distance"
        else:
            ref_low = dist_df["ref_dist_low_depth"].iloc[0]
            ref_large = dist_df["ref_dist_large_clones"].iloc[0]
            y_low = "sum_dist_low_depth"
            y_large = "sum_dist_large_clones"
            ylabel = "Sum distance"
            


    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(
        dist_df["w"],
        dist_df[y_low],
        "-o",
        lw=2,
        ms=5,
        label="Low-depth clones"
    )

    ax.plot(
        dist_df["w"],
        dist_df[y_large],
        "-o",
        lw=2,
        ms=5,
        label="Large clones"
    )

    ax.axhline(
        ref_low,
        color="C0",
        linestyle="--",
        linewidth=2,
        label="VAF vs VAF_AM (low-depth)"
    )

    ax.axhline(
        ref_large,
        color="C1",
        linestyle="--",
        linewidth=2,
        label="VAF vs VAF_AM (large)"
    )

    ax.set_xlabel("Pseudocount weight (w)")
    ax.set_ylabel(ylabel)

    if log_y:
        ax.set_yscale("log")

    ax.legend(frameon=False)
    sns.despine()
    plt.tight_layout()

    return fig, ax



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

    samples_list = sorted(mutations_table['SAMPLE_ID'].unique())

    priors_per_sample = compute_priors(mutdensities, samples_list)

    all_mutations_table = apply_correction_per_sample(mutations_table, priors_per_sample, depths_per_sample, samples_list)
    priors_n_depth = priors_per_sample.merge(depths_per_sample, on='SAMPLE_ID', how='left')

    mutations_plus_info = all_mutations_table.merge(priors_n_depth, on='SAMPLE_ID', how='left')
    mutations_plus_info = mutations_plus_info[['SAMPLE_ID', 'MUT_ID', 'prior', 'avg_depth_sample'] + [x for x in all_mutations_table.columns if x not in ['SAMPLE_ID', 'MUT_ID', 'prior', 'avg_depth_sample']]]
    mutations_plus_info.to_csv(f"mutations_with_smoothed_VAF.tsv",
                                header=True,
                                index=False,
                                sep="\t")

    with PdfPages("vaf_pseudocounts.pdf") as pdf:
        plt.figure(figsize=(8, 6))
        plt.hist(priors_per_sample['prior'], bins=10)
        plt.title('Distribution of sample-specific priors')
        plt.xlabel('Prior VAF')
        plt.ylabel('Frequency')
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        plot_vaf_pseudocount_curve(all_mutations_table, samples_list, priors_per_sample, pdf, suffix='')
        plot_vaf_pseudocount_curve(all_mutations_table, samples_list, priors_per_sample, pdf, suffix='_AM')

        weights = [0, 0.2, 0.5, 0.75, 1, 1.25, 1.5, 2, 3]

        dist_df = calc_vaf_distance_summary(
            all_mutations_table,
            weights=weights,
            vaf_col="VAF",
            target_col="VAF_AM",
            pseudo_prefix="VAF_PSEUDO",
            depth_col="DEPTH",
            depth_quantile=0.2,
            large_clone_threshold=1e-2,
        )

        fig, ax = plot_vaf_distance_summary(
            dist_df,
            log_dist=False,   
            use_mean=False,
            log_y=False
        )
        pdf.savefig()
        plt.close()


if __name__ == '__main__':
    main()

