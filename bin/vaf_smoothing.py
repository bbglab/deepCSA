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
    npa_priors = npa_mutdensities[['SAMPLE_ID', 'prior']]

    return npa_priors

# Reference: original VAF vs target VAF
def compute_column_distances(df_subset, vaf_col, target_col):
    """Compute the distance between two columns in a DataFrame."""
    dist_log = np.abs(np.log10(df_subset[vaf_col]) - np.log10(df_subset[target_col]))
    dist = np.abs(df_subset[vaf_col] - df_subset[target_col])
    reldist = np.abs(df_subset[vaf_col] - df_subset[target_col]) / df_subset[target_col]

    # return dist_log, dist, reldist
    return reldist




def apply_correction_per_sample(mutations_table, priors_per_sample,
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
        ].copy()
        sample_mutations["avg_depth_mutations"] = sample_mutations["DEPTH"].mean()
        sample_mutations["avg_depth_am_mutations"] = sample_mutations["DEPTH_AM"].mean()
        for weight_prop in weights_list:
            sample_mutations[f'VAF_PSEUDO_{weight_prop}'] = vaf_pseudocount(
                                        sample_mutations['ALT_DEPTH'],
                                        sample_mutations['DEPTH'],
                                        weight_prop,
                                        prior_vaf=prior_vaf
                                    )
            sample_mutations[f'VAF_AM_PSEUDO_{weight_prop}'] = vaf_pseudocount(
                                                    sample_mutations['ALT_DEPTH_AM'],
                                                    sample_mutations['DEPTH_AM'],
                                                    weight_prop,
                                                    prior_vaf=prior_vaf
                                                )

        mutation_tables_versions.append(sample_mutations)
    
    return pd.concat(mutation_tables_versions, ignore_index=True)


def vaf_pseudocount(alt_depth, depth, weight, prior_vaf=None):    
    return (alt_depth + prior_vaf * weight) / (depth + weight)

def plot_vaf_pseudocount_curve(maf_df, samples, output_pdf, suffix='',
                               weights_list=[0, 0.2, 0.5, 0.75, 1, 1.25, 1.5, 2, 3]):
    """
    Plot VAF distribution compared to VAF_AM in a histogram.
    
    Parameters:
    -----------
    maf_df : DataFrame
        MAF dataframe containing VAF and VAF_AM columns
    """

    n_cols = 5 if len(weights_list) > 9 else 3
    n_rows = len(weights_list) // n_cols if len(weights_list) % n_cols == 0 else len(weights_list) // n_cols + 1
    for sample in samples:
        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(n_cols * 2.5, n_rows * 2.5), sharex=True, sharey=True)
        axes = axes.flatten()
        fig.suptitle(f'{sample} : VAF{suffix} with pseudocounts')

        samples_maf_df = maf_df[maf_df['SAMPLE_ID'] == sample]
        for i, weight_prop in enumerate(weights_list):
            axes[i].scatter(samples_maf_df[f'DEPTH{suffix}'], samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}'], s=3, alpha=0.1)
            rho = np.corrcoef(np.log(samples_maf_df[f'DEPTH{suffix}']), np.log(samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}']))[0, 1]
            axes[i].set_xlabel(f'DEPTH{suffix}')
            axes[i].set_ylabel(f'VAF{suffix}_PSEUDO', fontsize=6)
            axes[i].set_yscale('log')
            axes[i].set_xscale('log')
            axes[i].set_title(f"{weight_prop}\nrho={rho:.2f}, \nprop_mut_tissue={samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}'].sum():.2f}", fontsize=6)
        plt.tight_layout()
        output_pdf.savefig()
        plt.close()
        plt.show()





def plot_vaf_pseudocount_curve_pair(maf_df, samples, selected_weights, suffix=''):
    """
    Plot VAF distribution compared to VAF_AM in a histogram.
    
    Parameters:
    -----------
    maf_df : DataFrame
        MAF dataframe containing VAF and VAF_AM columns
    """
    with PdfPages("vaf_pre_post_comparison.pdf") as output_pdf:

        n_cols = 2
        n_rows = 1
        for sample in samples:
            fig, axes = plt.subplots(n_rows, n_cols,
                                    figsize=(n_cols * 2.5, n_rows * 2.5), sharex=True, sharey=True)
            axes = axes.flatten()
            fig.suptitle(f'{sample} : VAF{suffix} with pseudocounts', fontsize = 7)

            samples_maf_df = maf_df[maf_df['SAMPLE_ID'] == sample]
            sample_selected_weight = selected_weights.loc[selected_weights['SAMPLE_ID'] == sample, 'selected_weight'].values[0]

            mean_vaf = samples_maf_df[f'VAF{suffix}'].mean()

            for i, weight_prop in enumerate([0, sample_selected_weight.item()]):
                axes[i].scatter(samples_maf_df[f'DEPTH{suffix}'], samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}'], s=3, alpha=0.1)
                rho = np.corrcoef(np.log(samples_maf_df[f'DEPTH{suffix}']), np.log(samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}']))[0, 1]
                axes[i].axhline(y=mean_vaf, color='red', linestyle='--', alpha=0.5, linewidth=0.5)
                axes[i].set_xlabel(f'DEPTH{suffix}', fontsize=6)
                axes[i].set_ylabel(f'VAF{suffix}_PSEUDO', fontsize=6)
                axes[i].set_yscale('log')
                axes[i].set_xscale('log')
                axes[i].tick_params(labelsize=6)
                axes[i].set_title(f"{weight_prop}\nrho={rho:.2f}, \nprop_mut_tissue={samples_maf_df[f'VAF{suffix}_PSEUDO_{weight_prop}'].sum():.2f}", fontsize=6)
            plt.tight_layout()
            output_pdf.savefig()
            plt.close()
            plt.show()



def plot_vaf_pseudocount_curve_groups(sample, df, df_low, df_large, df_high, suffix='', pdf_to_store=None):
    """
    Plot VAF distribution compared to VAF_AM in a histogram.
    
    Parameters:
    -----------
    sample : str
        Sample ID
    df : DataFrame
        MAF dataframe containing VAF and VAF_AM columns
    df_low : DataFrame
        MAF dataframe containing low depth clones
    df_large : DataFrame
        MAF dataframe containing large clones
    df_high : DataFrame
        MAF dataframe containing high depth clones
    suffix : str
        Suffix to append to column names
    pdf_to_store : PdfPages
        PDF pages object to store the plot
    """

    n_cols = 4
    n_rows = 1
    fig, axes = plt.subplots(n_rows, n_cols,
                            figsize=(n_cols * 2.5, n_rows * 2.5), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(f'{sample} groups : VAF{suffix} vs DEPTH{suffix}', fontsize = 7)

    names = ["ALL", "Low Depth", "Large Clones", "High Depth"]
    for i, samples_maf_df in enumerate([df, df_low, df_large, df_high]):
        axes[i].scatter(samples_maf_df[f'DEPTH{suffix}'], samples_maf_df[f'VAF{suffix}'], s=3, alpha=0.1)
        axes[i].set_xlabel(f'DEPTH{suffix}', fontsize=6)
        axes[i].set_ylabel(f'VAF{suffix}', fontsize=6)
        axes[i].set_yscale('log')
        axes[i].set_xscale('log')
        axes[i].tick_params(labelsize=6)
        axes[i].set_title(names[i])
    plt.tight_layout()
    if pdf_to_store is not None:
        pdf_to_store.savefig()
    plt.close()
    plt.show()



def calc_vaf_distance_summary(
    df,
    weights,
    target_col="VAF_AM",
    pseudo_prefix="VAF_PSEUDO",
    depth_col="DEPTH",
    depth_quantile=0.1,
    large_clone_quantile=0.95,
    output_pdf=None
    ):
    """
    Calculate distance between VAF_PSEUDO and target VAF column
    for low-depth clones and large clones.
    """

    depth_cutoff_low = df[depth_col].quantile(depth_quantile)
    depth_cutoff_high = df[depth_col].quantile(1 - depth_quantile)

    df_low_depth = df[(df[depth_col] < depth_cutoff_low)
                      & (df['ALT_DEPTH'] < 2)].copy()

    # define large clones except from low depth clones
    large_clone_quantile_vaf_value = df[(df[depth_col] > depth_cutoff_low)
                                        ][target_col].quantile(large_clone_quantile)

    df_large = df[(df[depth_col] > depth_cutoff_low)
                  & (df[target_col] > large_clone_quantile_vaf_value)].copy()

    df_large_ALTDEPTH = df[(df[depth_col] > depth_cutoff_low)
                            & (df["ALT_DEPTH"] > 1)
                            ].copy()

    df_high_depth = df[(df[depth_col] > depth_cutoff_high)].copy()

    plot_vaf_pseudocount_curve_groups(df_high_depth["SAMPLE_ID"].iloc[0], df, df_low_depth, df_large, df_high_depth, '_AM', output_pdf)

    print(f"Low depth clones: {df_low_depth.shape[0]}")
    print(f"Large clones: {df_large.shape[0]}")
    print(f"Large clones ALTDEPTH: {df_large_ALTDEPTH.shape[0]}")
    print(f"High depth clones: {df_high_depth.shape[0]}")

    results = []

    for w in weights:

        col = f"{pseudo_prefix}_{w}"

        reldist_low_depth       = compute_column_distances(df_low_depth, col, target_col)
        reldist_large           = compute_column_distances(df_large, col, target_col)
        reldist_largeALTDEPTH   = compute_column_distances(df_large_ALTDEPTH, col, target_col)
        reldist_high_depth      = compute_column_distances(df_high_depth, col, target_col)

        results.append({
            "w": w,

            "mean_reldist_low_depth": reldist_low_depth.mean().round(6),
            "mean_reldist_large_clones": reldist_large.mean().round(6),
            "mean_reldist_largeALTDEPTH": reldist_largeALTDEPTH.mean().round(6),
            "mean_reldist_high_depth": reldist_high_depth.mean().round(6),

            "n_total": len(df),
            "n_low_depth": len(df_low_depth),
            "n_large_clones": len(df_large),
            "n_largeALTDEPTH": len(df_large_ALTDEPTH),
            "n_high_depth": len(df_high_depth),
            "depth_cutoff_low": depth_cutoff_low,
            "depth_cutoff_high": depth_cutoff_high,
            "large_clone_definition_vaf": large_clone_quantile_vaf_value.round(6)
        })

    return pd.DataFrame(results)



def plot_vaf_distance_summary(
    dist_df,
    figsize=(7, 4),
    log_y=False,
    sample_name = "sample"
    ):
    """
    Plot pseudocount distance curves and reference VAF vs VAF_AM distances.
    """

    y_low       = "mean_reldist_low_depth"
    y_large     = "mean_reldist_large_clones"
    y_large_alt = "mean_reldist_largeALTDEPTH"
    y_high_depth= "mean_reldist_high_depth"
    ylabel      = "Mean relative distance of updated VAF"


    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(dist_df["w"], dist_df[y_low],
            "-o", lw=2, ms=5, label=f"Low-depth clones (n={dist_df['n_low_depth'].iloc[0]})"
            )

    ax.plot(dist_df["w"], dist_df[y_large],
            "-o", lw=2, ms=5, label=f"Large clones (n={dist_df['n_large_clones'].iloc[0]})"
            )

    if dist_df['n_largeALTDEPTH'].iloc[0] > 0:
        ax.plot( dist_df["w"], dist_df[y_large_alt],
                "-o", lw=2, ms=5, label=f"Large clones ALTDEPTH (n={dist_df['n_largeALTDEPTH'].iloc[0]})"
                )

    if dist_df['n_high_depth'].iloc[0] > 0:
        ax.plot(dist_df["w"], dist_df[y_high_depth],
                "-o", lw=2, ms=5, label=f"Large clones HIGH_DEPTH (n={dist_df['n_high_depth'].iloc[0]})"
                )

    ax.set_xlabel("Pseudocount weight (w)")
    ax.set_ylabel(ylabel)

    if log_y:
        ax.set_yscale("log")

    plt.title(f"{sample_name} (n={dist_df['n_total'].iloc[0]})\nlow_depth cutoff: {dist_df['depth_cutoff_low'].iloc[0]:.2f}, large clone definition: {dist_df['large_clone_definition_vaf'].iloc[0]:.2e}")
    ax.legend(frameon=False)
    sns.despine()
    plt.tight_layout()

    return fig, ax



def select_vafpseudo_per_sample(selected_weights_df, mutations_table):
    """
    Select the best VAF_PSEUDO column for each sample based on the selected weights.
    """

    selected_vafpseudo_columns = []
    columns_to_keep = ['SAMPLE_ID', 'MUT_ID', 'ALT_DEPTH', 'DEPTH', 'VAF', 'ALT_DEPTH_AM', 'DEPTH_AM', 'VAF_AM', 'prior', 'avg_depth_sample']
    for _, row in selected_weights_df.iterrows():
        sample_id = row['SAMPLE_ID']
        selected_weight = row['selected_weight']
        selected_column = f'VAF_PSEUDO_{selected_weight}'
        selected_vafpseudo_columns.append((sample_id, selected_column, selected_weight))

    # Create a new DataFrame with the selected VAF_PSEUDO columns
    selected_vafpseudo_df = pd.DataFrame(selected_vafpseudo_columns, columns=['SAMPLE_ID', 'selected_column', 'selected_weight'])

    # Merge with the original mutations table to get the selected VAF_PSEUDO values
    merged_df = mutations_table.merge(selected_vafpseudo_df, on='SAMPLE_ID', how='left')

    # Create a new column for the final selected VAF_PSEUDO values
    merged_df['selected_VAF_PSEUDO'] = merged_df.apply(lambda x: x[x['selected_column']], axis=1)
    merged_df = merged_df[columns_to_keep + ['selected_weight', 'selected_VAF_PSEUDO', 'avg_depth_mutations', 'avg_depth_am_mutations']]

    # Save the final DataFrame with selected VAF_PSEUDO values
    merged_df.to_csv("mutations_with_final_selected_VAF_PSEUDO.tsv.gz", sep="\t", index=False)

    corrections_summary_df = merged_df[['SAMPLE_ID', 'prior', 'avg_depth_sample', 'selected_weight',
                                        'avg_depth_mutations', 'avg_depth_am_mutations']].drop_duplicates()

    # Save the final DataFrame with a summary of the corrections applied per sample
    corrections_summary_df.to_csv("corrections_per_sample_summary.tsv", sep="\t", index=False)


@click.command()
@click.option('--mutdensities', type=click.Path(exists=True), help='Input mutation density file')
@click.option('--mutations', type=click.Path(), help='Mutations file')
@click.option('--depth-sample', type=click.Path(), help='Depth per sample file')
def main(mutdensities, mutations, depth_sample):
    """
    Main function to execute the VAF smoothing pipeline.
    """
    weights = [x for x in np.arange(0, 20001, 500)]

    click.echo("Loading data...")
    mutations_table = pd.read_csv(mutations, sep = "\t", header = 0, na_values = custom_na_values)
    mutations_0_vaf = mutations_table[~(mutations_table["VAF"] > 0)
                                          | ~(mutations_table["VAF_AM"] > 0)
                                          ]
    if not mutations_0_vaf.empty:
        click.echo(f"Warning: {mutations_0_vaf.shape[0]} mutations have VAF or VAF_AM <= 0 and will be excluded from the analysis.")
        click.echo("These mutations are:")
        click.echo(mutations_0_vaf[['SAMPLE_ID', 'MUT_ID', 'VAF', 'VAF_AM']])
        mutations_table = mutations_table[(mutations_table["VAF"] > 0)
                                          & (mutations_table["VAF_AM"] > 0)
                                          ].reset_index(drop = True)

    depths_per_sample = pd.read_csv(depth_sample, sep = "\t", header = 0, na_values = custom_na_values)

    samples_list = sorted(mutations_table['SAMPLE_ID'].unique())

    click.echo("Computing priors...")
    priors_per_sample = compute_priors(mutdensities, samples_list)

    click.echo("Applying corrections...")
    all_mutations_table = apply_correction_per_sample(mutations_table, priors_per_sample, samples_list, weights_list=weights)
    priors_n_depth = priors_per_sample.merge(depths_per_sample, on='SAMPLE_ID', how='left')

    mutations_plus_info = all_mutations_table.merge(priors_n_depth, on='SAMPLE_ID', how='left')
    mutations_plus_info = mutations_plus_info[['SAMPLE_ID', 'MUT_ID', 'prior', 'avg_depth_sample'] + [x for x in all_mutations_table.columns if x not in ['SAMPLE_ID', 'MUT_ID', 'prior', 'avg_depth_sample']]]
    mutations_plus_info.to_csv(f"mutations_with_smoothed_VAF.tsv.gz",
                                header=True,
                                index=False,
                                sep="\t")

    # click.echo("Plotting VAF pseudocount curves...")
    # with PdfPages("vaf_pseudocounts_curves.pdf") as pdf:
    #     plt.figure(figsize=(8, 6))
    #     plt.hist(priors_per_sample['prior'], bins=10)
    #     plt.title('Distribution of sample-specific priors')
    #     plt.xlabel('Prior VAF')
    #     plt.ylabel('Frequency')
    #     plt.tight_layout()
    #     pdf.savefig()
    #     plt.close()

    #     plot_vaf_pseudocount_curve(all_mutations_table, samples_list, pdf, suffix='', weights_list=weights)
    #     plot_vaf_pseudocount_curve(all_mutations_table, samples_list, pdf, suffix='_AM', weights_list=weights)

    click.echo("Plotting pseudocount curves weight comparison...")
    with PdfPages("vaf_groups_per_sample.pdf") as pdf_separation:
        with PdfPages("vaf_pseudocount_weights_comparison.pdf") as pdf:

            selected_weights_dict = {}
            try :
                for sample in samples_list:
                    sample_mutations = all_mutations_table[all_mutations_table['SAMPLE_ID'] == sample]

                    dist_df = calc_vaf_distance_summary(
                        sample_mutations,
                        weights=weights,
                        target_col="VAF_AM",
                        pseudo_prefix="VAF_AM_PSEUDO",
                        depth_col="DEPTH_AM",
                        depth_quantile=0.2,
                        large_clone_quantile=0.90,
                        output_pdf=pdf_separation
                    )
                    dist_df["mean_reldist_diff_low_over_large_clones"] = (dist_df["mean_reldist_low_depth"] / dist_df["mean_reldist_large_clones"]).round(4)
                    dist_df["mean_reldist_diff_low_over_high_depth"] = (dist_df["mean_reldist_low_depth"] / dist_df["mean_reldist_high_depth"]).round(4)
                    dist_df.to_csv(f"vaf_distance_summary_{sample}.tsv", sep="\t", index=False)

                    #get minimum distance weight for low depth clones
                    try:
                        max_weight_below_threshold = dist_df.loc[np.where(dist_df['mean_reldist_large_clones'] < 0.2)[0][-1], 'w']
                    except IndexError:
                        click.echo(f"Warning: No weight found for sample {sample} where mean_reldist_large_clones < 0.2. Using the smallest non-0 weight.")
                        max_weight_below_threshold = weights[1]  # or some default value

                    selected_weights_dict[sample] = max_weight_below_threshold

                    print(f"Sample {sample}: Maximum weight below threshold for low depth clones: {max_weight_below_threshold}")

                    fig, ax = plot_vaf_distance_summary( dist_df, log_y=False, sample_name = sample)
                    pdf.savefig()
                    plt.close()


                dist_df = calc_vaf_distance_summary(
                    all_mutations_table,
                    weights=weights,
                    target_col="VAF_AM",
                    pseudo_prefix="VAF_AM_PSEUDO",
                    depth_col="DEPTH_AM",
                    depth_quantile=0.2,
                    large_clone_quantile=0.90,
                    output_pdf=pdf_separation
                )
                dist_df["mean_reldist_diff_low_over_large_clones"] = (dist_df["mean_reldist_low_depth"] / dist_df["mean_reldist_large_clones"]).round(4)
                dist_df["mean_reldist_diff_low_over_largeALTDEPTH"] = (dist_df["mean_reldist_low_depth"] / dist_df["mean_reldist_largeALTDEPTH"]).round(4)
                dist_df["mean_reldist_diff_low_over_high_depth"] = (dist_df["mean_reldist_low_depth"] / dist_df["mean_reldist_high_depth"]).round(4)
                dist_df.to_csv(f"vaf_distance_summary_all_samples.tsv", sep="\t", index=False)

                fig, ax = plot_vaf_distance_summary(dist_df, log_y=False, sample_name = 'all_samples')
                pdf.savefig()
                plt.close()

                #get minimum distance weight for low depth clones
                try:
                    max_weight_below_threshold = dist_df.loc[np.where(dist_df['mean_reldist_large_clones'] < 0.2)[0][-1], 'w']
                except IndexError:
                    click.echo(f"Warning: No weight found for sample {sample} where mean_reldist_large_clones < 0.2. Using the smallest non-0 weight.")
                    max_weight_below_threshold = weights[1]  # or some default value

                selected_weights_dict['all_samples'] = max_weight_below_threshold

                # store selected weights in a tsv file
                selected_weights_df = pd.DataFrame(list(selected_weights_dict.items()), columns=['SAMPLE_ID', 'selected_weight'])

            except Exception as e:
                click.echo(f"Error in plotting VAF distance summary: {e}")

    plot_vaf_pseudocount_curve_pair(all_mutations_table, samples_list, selected_weights_df, suffix='_AM')

    click.echo("Storing summary of priors and weights used.")
    select_vafpseudo_per_sample(selected_weights_df, mutations_plus_info)

    click.echo("VAF smoothing pipeline completed successfully.")



if __name__ == '__main__':
    main()

