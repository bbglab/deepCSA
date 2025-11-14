#!/usr/bin/env python

"""
Script for checking contamination between samples.
"""


import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from read_utils import custom_na_values

# TODO: bump pandas to 2.2.3

# -- Auxiliary functions -- #

# # Plot heatmap
# def plot_shared_variants_heatmap(shared_counts):
#     plt.figure(figsize=(23, 20))
#     sns.heatmap(shared_counts,
#                 annot=shared_counts.applymap(lambda x: x if x > 30 else ""),
#                 fmt="",
#                 cmap="Blues",
#                 linewidths=0.5)
#     #sns.clustermap(shared_counts, annot=shared_counts.applymap(lambda x: x if x > 0 else ""), fmt="", cmap="Blues", linewidths=0.5)
#     #sns.heatmap(shared_counts, annot=True, cmap='Blues', fmt='g', linewidths=0.5)
#     plt.xlabel("Germline Samples")
#     plt.ylabel("Somatic Samples")
#     plt.title("Shared Variants Between Somatic and Germline Samples")
#     plt.show()

# # Example usage:
# # plot_shared_variants_heatmap(shared_variants_matrix)



# Assuming somatic_variants and germline_variants are loaded as pandas DataFrames
def compute_shared_variants(somatic_variants, germline_variants):
    """
    # Example usage:
    # shared_variants_matrix = compute_shared_variants(somatic_variants, germline_variants)
    """

    unique_somatic_samples = sorted(somatic_variants['SAMPLE_ID'].unique())
    unique_germline_samples = sorted(germline_variants['SAMPLE_ID'].unique())

    # Create a DataFrame to store counts
    shared_counts = pd.DataFrame(index=unique_somatic_samples, columns=unique_germline_samples).fillna(0)

    # Iterate through each somatic sample
    for somatic_sample in unique_somatic_samples:
        somatic_mutations = set(somatic_variants[somatic_variants['SAMPLE_ID'] == somatic_sample]['MUT_ID'])

        # Compare with germline mutations of all other samples
        for germline_sample in unique_germline_samples:
            germline_mutations = set(germline_variants[germline_variants['SAMPLE_ID'] == germline_sample]['MUT_ID'])

            # Count shared mutations
            shared_counts.loc[somatic_sample, germline_sample] = len(somatic_mutations & germline_mutations)

    return shared_counts




def contamination_detection(maf_file, somatic_maf_file):

    maf_df = pd.read_table(maf_file, na_values=custom_na_values)
    print(maf_df.shape)
    maf_df = maf_df[~(maf_df["FILTER.not_covered"])
                    & (maf_df["TYPE"] == 'SNV')
                    ].reset_index()
    print(maf_df.shape)

    somatic_maf_df = pd.read_table(somatic_maf_file, na_values=custom_na_values)
    print(somatic_maf_df.shape)
    somatic_maf_df = somatic_maf_df[(somatic_maf_df["TYPE"] == 'SNV')]
    print(somatic_maf_df.shape)


    # this is if we were to consider both unique and no-unique variants
    germline_vars_all_samples = maf_df.loc[(maf_df["VAF"] > 0.3) & (maf_df["vd_VAF"] > 0.3) & (maf_df["VAF_AM"] > 0.3),
                                            ["SAMPLE_ID", "MUT_ID"]].drop_duplicates()

    print(germline_vars_all_samples["MUT_ID"].shape)
    print(len(germline_vars_all_samples["MUT_ID"].unique()))


    somatic_variants = somatic_maf_df[["SAMPLE_ID", "MUT_ID"]]
    print(somatic_variants.shape)


    all_variants = maf_df[["SAMPLE_ID", "MUT_ID"]]
    print(all_variants.shape)


    ## Somatic vs Germline

    shared_variants_somatic2germline_matrix = compute_shared_variants(somatic_variants, germline_vars_all_samples)

    plt.figure(figsize=(18, 15))

    # Compute total number of germline mutations per sample
    germline_counts = germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_variants_somatic2germline_matrix.columns)

    # Create custom column labels with germline mutation counts
    col_labels = [f"(n={germline_counts[col]}) {col}" for col in shared_variants_somatic2germline_matrix.columns]

    sns.heatmap(
        shared_variants_somatic2germline_matrix,
        annot=shared_variants_somatic2germline_matrix.applymap(lambda x: x if x > 30 else ""),
        fmt="",
        cmap="Blues",
        cbar_kws={'label': 'Shared Mutations'},
        xticklabels=col_labels,
        yticklabels=shared_variants_somatic2germline_matrix.index,
        linewidths=0.5,
        annot_kws={"color": "black", "fontsize": 10}  # <- Set annotation text color and size here
    )

    plt.xlabel("Germline Samples", fontsize=14)
    plt.ylabel("Somatic Samples", fontsize=14)
    plt.title("Somatic mutations that are germline in other samples", fontsize=16)
    plt.savefig("somatic_vs_germline.pdf", bbox_inches = 'tight', dpi = 100)
    plt.show()





    ## All vs Germline

    shared_all_vs_germline_variants_matrix = compute_shared_variants(all_variants, germline_vars_all_samples)

    # Compute total number of germline mutations per sample
    germline_counts = germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_all_vs_germline_variants_matrix.columns)


    normalized_shared_all_vs_germline_variants_matrix = shared_all_vs_germline_variants_matrix.divide(germline_counts, axis=1)



    # Count shared mutations between somatic and germline samples

    plt.figure(figsize=(18, 15))


    # Create custom column labels with germline mutation counts
    col_labels = [f"(n={germline_counts[col]}) {col}" for col in normalized_shared_all_vs_germline_variants_matrix.columns]


    sns.heatmap(normalized_shared_all_vs_germline_variants_matrix,
                annot=normalized_shared_all_vs_germline_variants_matrix.applymap(lambda x: round(x,2) if (x > 0.8) and (x < 1) else ""),
                fmt="",
                cmap="Blues",
                cbar_kws={'label': 'Shared Mutations'},
                xticklabels=col_labels, yticklabels=normalized_shared_all_vs_germline_variants_matrix.index,
                annot_kws={"color": "white", "fontsize": 10},
                linewidths=0.5)

    plt.xlabel("Germline Samples", fontsize = 14)
    plt.ylabel("All mutations samples", fontsize = 14)
    plt.title("All mutations that are germline in other samples", fontsize = 16)
    plt.savefig("allmutations_vs_germline.pdf", bbox_inches = 'tight', dpi = 100)
    plt.show()







    ## Germline vs Germline

    shared_germline_variants_matrix = compute_shared_variants(germline_vars_all_samples, germline_vars_all_samples)

    plt.figure(figsize=(18, 15))

    # Compute total number of germline mutations per sample
    germline_counts = germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_germline_variants_matrix.columns)

    # Create custom column labels with germline mutation counts
    col_labels = [f"(n={germline_counts[col]}) {col}" for col in shared_germline_variants_matrix.columns]


    sns.heatmap(shared_germline_variants_matrix,
                annot=shared_germline_variants_matrix.applymap(lambda x: x if x < 0 else ""),
                fmt="",
                cmap="Blues",
                cbar_kws={'label': 'Shared Mutations'},
                xticklabels=col_labels, yticklabels=shared_germline_variants_matrix.index,
                linewidths=0.5,
                annot_kws={"fontsize": 8}
            )

    plt.xlabel("Germline Samples", fontsize = 14)
    plt.ylabel("Germline Samples", fontsize = 14)
    plt.title("Germline mutations that are germline in other samples", fontsize = 16)
    plt.savefig("germline_vs_germline.pdf", bbox_inches = 'tight', dpi = 100)
    plt.show()





    # Compute total number of germline mutations per sample
    germline_counts = germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_germline_variants_matrix.columns)

    normalized_share_germline_vs_germline_variants_matrix = shared_germline_variants_matrix.divide(germline_counts, axis=1)


    plt.figure(figsize=(18, 15))


    # Create custom column labels with germline mutation counts
    col_labels = [f"(n={germline_counts[col]}) {col}" for col in normalized_share_germline_vs_germline_variants_matrix.columns]


    sns.heatmap(normalized_share_germline_vs_germline_variants_matrix,
                annot=normalized_share_germline_vs_germline_variants_matrix.applymap(lambda x: round(x,2) if (x > 0.8) and (x < 1) else ""),
                fmt="",
                cmap="Blues",
                cbar_kws={'label': 'Shared Mutations'},
                xticklabels=col_labels, yticklabels=normalized_share_germline_vs_germline_variants_matrix.index,
                annot_kws={"color": "white", "fontsize": 10},
                linewidths=0.5)

    plt.xlabel("Germline Samples", fontsize = 14)
    plt.ylabel("Germline samples", fontsize = 14)
    plt.savefig("normalized.germline_vs_germline.pdf", bbox_inches = 'tight', dpi = 100)
    plt.show()




    ## Somatic vs Remaining Germline
    shared_somatic_to_non_shared_germline = shared_all_vs_germline_variants_matrix - shared_germline_variants_matrix


    # Compute total number of germline mutations per sample
    germline_counts = germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_somatic_to_non_shared_germline.columns)


    total_germline_available_per_sample = (germline_counts - shared_germline_variants_matrix)

    shared_somatic_to_non_shared_germline_proportion = (shared_somatic_to_non_shared_germline / total_germline_available_per_sample).fillna(0)



    plt.figure(figsize=(22, 18))

    sns.heatmap(shared_somatic_to_non_shared_germline_proportion,
                annot=shared_somatic_to_non_shared_germline_proportion.applymap(lambda x: round(x, 2) if x > 0.45 else ""),
                fmt="",
                cmap="Blues",
                cbar_kws={'label': 'Shared Mutations'},
                # xticklabels=col_labels,
                yticklabels=shared_somatic_to_non_shared_germline_proportion.index,
                annot_kws={"color": "black", "fontsize": 10},
                linewidths=0.5)

    plt.xlabel("Non-shared germline", fontsize = 14)
    plt.ylabel("Somatic", fontsize = 14)
    plt.title("Somatic mutations that are germline in other samples", fontsize = 16)
    plt.savefig("contamination.somatic_vs_remaininggermline.pdf", bbox_inches = 'tight', dpi = 100)
    plt.show()





    max_prop_per_sample = shared_somatic_to_non_shared_germline_proportion.max(axis = 'columns')
    max_prop_per_sample[max_prop_per_sample>0.5]

    max_prop_per_sample.to_csv(f"contaminated_samples.tsv",
                            header = True,
                            sep  = '\t',
                            index = True)

    contaminated_samples = list(max_prop_per_sample[max_prop_per_sample>0.5].index.values)
    contaminated_samples




    ## Exploration of contaminated samples
    for sample, max_val in max_prop_per_sample[max_prop_per_sample>0.5].reset_index().values:
        sample_vals = shared_somatic_to_non_shared_germline_proportion.loc[sample,:]
        source_sampleid = sample_vals[sample_vals == max_val].index.values[0]


        print(f'{sample} has {max_val:.2f} proportion of the germline variants of {source_sampleid} as with a VAF not corresponding to germline variants.')
        print()

        subseeeet = maf_df[["SAMPLE_ID", "MUT_ID", "VAF", "ALT_DEPTH"]]
        p_dest = subseeeet[subseeeet["SAMPLE_ID"] == sample].drop("SAMPLE_ID", axis = 1)

        p_source_germ = germline_vars_all_samples[germline_vars_all_samples["SAMPLE_ID"] == source_sampleid]
        p_source = subseeeet[(subseeeet["SAMPLE_ID"] == source_sampleid)
                                & (subseeeet["MUT_ID"].isin(p_source_germ["MUT_ID"].values)) ].drop("SAMPLE_ID", axis = 1)

        merged_samples = p_dest.merge(p_source,
                                on = ["MUT_ID"],
                                suffixes = ("_dest", "_source"),
                                how = 'right'
                                )

        merged_samples.to_csv(f"{source_sampleid}.germline_variants_in.{sample}.tsv",
                                header = True,
                                sep  = '\t',
                                index = False)

        plt.scatter(x = merged_samples["VAF_dest"].fillna(0),
                    y = merged_samples["VAF_source"].fillna(0),
                    # color = ['blue' if x == 0 else 'red' for x in merged_samples["VAF_dest"].fillna(0)]
                )

        plt.xscale('log')
        # plt.yscale('log')
        plt.xlabel("VAF_dest    "   + sample)
        plt.ylabel("VAF_source  " + source_sampleid)
        plt.show()




@click.command()
@click.option('--maf_path', type=click.Path(exists=True), required=True, help='Path to the MAF file.')
@click.option('--somatic_maf', type=click.Path(exists=True), required=True, help='Path to the depths file.')
def main(maf_path, somatic_maf):
    """
    CLI entry point for computing mutation densities.
    """
    contamination_detection(maf_path, somatic_maf)



if __name__ == '__main__':

    main()
