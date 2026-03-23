#!/usr/bin/env python

"""
Script for checking contamination between samples.
"""


import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from read_utils import custom_na_values


# Assuming somatic_variants and germline_variants are loaded as pandas DataFrames
def compute_shared_variants(somatic_variants, germline_variants):
    """
    # Example usage:
    # shared_variants_matrix = compute_shared_variants(somatic_variants, germline_variants)
    """

    unique_somatic_samples = sorted(somatic_variants['SAMPLE_ID'].unique())
    unique_germline_samples = sorted(germline_variants['SAMPLE_ID'].unique())

    # Create a DataFrame to store counts (avoid .fillna downcasting warning)
    shared_counts = pd.DataFrame(0, index=unique_somatic_samples, columns=unique_germline_samples, dtype=int)

    # Iterate through each somatic sample
    for somatic_sample in unique_somatic_samples:
        somatic_mutations = set(somatic_variants[somatic_variants['SAMPLE_ID'] == somatic_sample]['MUT_ID'])

        # Compare with germline mutations of all other samples
        for germline_sample in unique_germline_samples:
            germline_mutations = set(germline_variants[germline_variants['SAMPLE_ID'] == germline_sample]['MUT_ID'])

            # Count shared mutations
            shared_counts.loc[somatic_sample, germline_sample] = len(somatic_mutations & germline_mutations)

    return shared_counts




def contamination_detection_between_samples(maf_df, somatic_maf_df):

    # this is if we were to consider both unique and no-unique variants
    vaf_threshold = 0.2
    germline_vars_all_samples = maf_df.loc[(maf_df["VAF"] > vaf_threshold) & (maf_df["vd_VAF"] > vaf_threshold) & (maf_df["VAF_AM"] > vaf_threshold),
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

    # Build annotation DataFrame without using deprecated DataFrame.applymap
    mask = shared_variants_somatic2germline_matrix > 30
    annot = shared_variants_somatic2germline_matrix.where(mask)
    # convert selected values to nullable int then to string, replace missing with empty string
    annot = annot.round(0).astype('Int64').astype(str).replace('<NA>', '').fillna('')

    sns.heatmap(
        shared_variants_somatic2germline_matrix,
        annot=annot,
        fmt="",
        cmap="Blues",
        cbar_kws={'label': 'Shared Mutations'},
        xticklabels=col_labels,
        yticklabels=shared_variants_somatic2germline_matrix.index,
        linewidths=0.5,
        annot_kws={"color": "black", "fontsize": 10}
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


    # Annotation: show rounded values only when 0.8 < x < 1
    cond = (normalized_shared_all_vs_germline_variants_matrix > 0.8) & (normalized_shared_all_vs_germline_variants_matrix < 1)
    annot = normalized_shared_all_vs_germline_variants_matrix.where(cond)
    annot = annot.round(2).astype('string').fillna('')

    sns.heatmap(normalized_shared_all_vs_germline_variants_matrix,
                annot=annot,
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


    # Annotation: follow original logic (keep values where < 0, else blank)
    mask = shared_germline_variants_matrix < 0
    annot = shared_germline_variants_matrix.where(mask)
    annot = annot.astype('string').fillna('')

    sns.heatmap(shared_germline_variants_matrix,
                annot=annot,
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


    cond = (normalized_share_germline_vs_germline_variants_matrix > 0.8) & (normalized_share_germline_vs_germline_variants_matrix < 1)
    annot = normalized_share_germline_vs_germline_variants_matrix.where(cond)
    annot = annot.round(2).astype('string').fillna('')

    sns.heatmap(normalized_share_germline_vs_germline_variants_matrix,
                annot=annot,
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

    # Those cases where the number of mutations is smaller than 5 are set to 0
    shared_somatic_to_non_shared_germline[shared_somatic_to_non_shared_germline < 5] = 0

    # Compute total number of germline mutations per sample
    germline_counts = germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_somatic_to_non_shared_germline.columns)


    total_germline_available_per_sample = (germline_counts - shared_germline_variants_matrix)

    shared_somatic_to_non_shared_germline_proportion = (shared_somatic_to_non_shared_germline / total_germline_available_per_sample).fillna(0)



    plt.figure(figsize=(22, 18))

    cond = shared_somatic_to_non_shared_germline_proportion > 0.45
    annot = shared_somatic_to_non_shared_germline_proportion.where(cond)
    annot = annot.round(2).astype('string').fillna('')

    sns.heatmap(shared_somatic_to_non_shared_germline_proportion,
                annot=annot,
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
    plt.close()


    plt.figure(figsize=(22, 18))

    cond = shared_somatic_to_non_shared_germline > 0
    annot = shared_somatic_to_non_shared_germline.where(cond)
    # convert to nullable int then string, replace missing with empty string
    annot = annot.round(0).astype('Int64').astype('string').replace('<NA>', '').fillna('')

    sns.heatmap(shared_somatic_to_non_shared_germline,
                annot=annot,
                fmt="",
                cmap="Blues",
                cbar_kws={'label': 'Shared Mutations'},
                # xticklabels=col_labels,
                yticklabels=shared_somatic_to_non_shared_germline.index,
                annot_kws={"color": "black", "fontsize": 10},
                linewidths=0.5)

    plt.xlabel("Non-shared germline", fontsize = 14)
    plt.ylabel("Somatic", fontsize = 14)
    plt.title("Somatic mutations that are germline in other samples (count)", fontsize = 16)
    plt.savefig("contamination.somatic_vs_remaininggermline.numbers.pdf", bbox_inches = 'tight', dpi = 100)
    plt.show()
    plt.close()


    max_prop_per_sample = shared_somatic_to_non_shared_germline_proportion.max(axis = 'columns')

    ## Exploration of contaminated samples
    receiver_source_pairs = []
    for sample, max_val in max_prop_per_sample[max_prop_per_sample>0.5].reset_index().values:
        sample_vals = shared_somatic_to_non_shared_germline_proportion.loc[sample,:]
        sample_vals_count = shared_somatic_to_non_shared_germline.loc[sample,:]

        source_sampleids = sample_vals[sample_vals == max_val].index.values
        source_sampleid = source_sampleids[0]
        receiver_source_pairs.append((sample, round(max_val,3),
                                      list(zip([sample_vals_count[x].item() for x in source_sampleids], source_sampleids))))

        print(f'{sample} has {max_val:.2f} proportion of the germline variants of {source_sampleid} as with a VAF not corresponding to germline variants.')
        print(f'Shared variants count: {sample_vals_count[source_sampleid]}')
        print()


        subseeeet = maf_df[["SAMPLE_ID", "MUT_ID", 'canonical_SYMBOL', "ALT_DEPTH", "DEPTH", "VAF", 'canonical_Consequence_broader', 'FILTER']]
        p_dest = subseeeet[subseeeet["SAMPLE_ID"] == sample].drop("SAMPLE_ID", axis = 1)

        p_source_germ = germline_vars_all_samples[germline_vars_all_samples["SAMPLE_ID"] == source_sampleid]
        p_source = subseeeet[(subseeeet["SAMPLE_ID"] == source_sampleid)
                                & (subseeeet["MUT_ID"].isin(p_source_germ["MUT_ID"].values))
                                ].drop("SAMPLE_ID", axis = 1)

        merged_samples = p_dest.merge(p_source,
                                on = ["MUT_ID", 'canonical_SYMBOL', 'canonical_Consequence_broader'],
                                suffixes = ("_dest", "_source"),
                                how = 'right'
                                )

        merged_samples.sort_values(by =["VAF_dest"], ascending=False
                                   ).to_csv(f"{source_sampleid}.germline_variants_in.{sample}.tsv",
                                            header = True,
                                            sep  = '\t',
                                            index = False)
        
        # plt.figure(figsize=(8, 6))
        # plt.scatter(x = merged_samples["VAF_dest"].fillna(0),
        #             y = merged_samples["VAF_source"].fillna(0),
        #             # color = ['blue' if x == 0 else 'red' for x in merged_samples["VAF_dest"].fillna(0)]
        #         )

        # plt.xscale('log')
        # # plt.yscale('log')
        # plt.xlabel("VAF_dest    "   + sample)
        # plt.ylabel("VAF_source  " + source_sampleid)
        # plt.savefig(f"{source_sampleid}_germline_in_{sample}_VAF_scatter.pdf", bbox_inches = 'tight', dpi = 100)
        # plt.show()

    # Store contamination results
    contamination_detailed_df = pd.DataFrame(receiver_source_pairs,
                                             columns=["SAMPLE_ID", "MAX_PROPORTION_GERMLINE_FROM_SOURCE", "SOURCE_SAMPLEID_COUNTS"])
    contamination_detailed_df.to_csv(f"contaminated_samples.detailed.tsv",
                                               header = True,
                                               sep  = '\t',
                                               index = False)

    if contamination_detailed_df.empty:
        print("No contaminated samples detected.")
        return
    contamination_detailed_df_long = contamination_detailed_df.explode("SOURCE_SAMPLEID_COUNTS")
    expanded_df = pd.DataFrame(contamination_detailed_df_long["SOURCE_SAMPLEID_COUNTS"].tolist())
    expanded_df.columns = ["SHARED_VARIANT_COUNT", "SOURCE_SAMPLEID"]
    contamination_detailed_df_long["SHARED_VARIANT_COUNT"] = expanded_df["SHARED_VARIANT_COUNT"].values
    contamination_detailed_df_long["SOURCE_SAMPLEID"] = expanded_df["SOURCE_SAMPLEID"].values
    contamination_detailed_df_long = contamination_detailed_df_long.drop("SOURCE_SAMPLEID_COUNTS", axis = 1)
    contamination_detailed_df_long.to_csv(f"contaminated_samples.detailed.long.tsv",
                                               header = True,
                                               sep  = '\t',
                                               index = False)



def data_loading(maf_path, somatic_maf_path):
    maf_df = pd.read_table(maf_path, na_values=custom_na_values)
    print(maf_df.shape)
    maf_df = maf_df[~(maf_df["FILTER.not_covered"])
                     & (maf_df["TYPE"] == 'SNV')
                    ].reset_index()
    print(maf_df.shape)

    somatic_maf_df = pd.read_table(somatic_maf_path, na_values=custom_na_values)
    print(somatic_maf_df.shape)
    somatic_maf_df = somatic_maf_df[(somatic_maf_df["TYPE"] == 'SNV')]
    print(somatic_maf_df.shape)
    return maf_df, somatic_maf_df


def contamination_detection_in_snps(maf):

    snp_positions_maf = maf[maf["FILTER.gnomAD_SNP"]][
        ["SAMPLE_ID", "MUT_ID", "VAF"]
        ].reset_index(drop = True)
    
    # being very restrictive in the VAF to count the occurrences of potentially contaminated mutations
    somatic_snp_positions_maf = snp_positions_maf[snp_positions_maf["VAF"] < 0.05].reset_index(drop = True)
    germline_snp_positions_maf = snp_positions_maf[snp_positions_maf["VAF"] >= 0.05].reset_index(drop = True)

    unique_SNP_positions = snp_positions_maf["MUT_ID"].unique()
    number_unique_SNP_positions = len(unique_SNP_positions)

    sample_SNP_mutation_freq = []
    for sample in snp_positions_maf["SAMPLE_ID"].unique():
        germline_count = len(germline_snp_positions_maf[germline_snp_positions_maf["SAMPLE_ID"] == sample])
        somatic_count = len(somatic_snp_positions_maf[somatic_snp_positions_maf["SAMPLE_ID"] == sample])
        remaining_germline = number_unique_SNP_positions-germline_count
        sample_SNP_mutation_freq.append([sample,
                                         germline_count,
                                         remaining_germline,
                                         somatic_count,
                                         somatic_count / remaining_germline if remaining_germline > 0 else 1
                                         ])
    sample_SNP_mutation_freq_df = pd.DataFrame(sample_SNP_mutation_freq)
    sample_SNP_mutation_freq_df.columns = ["SAMPLE_ID", "germline_count", "remaining_germline", "somatic_count", "prop_somatic_SNPs"]

    # identify outliers in the "prop_somatic_SNPs" column
    sample_SNP_mutation_freq_df = sample_SNP_mutation_freq_df.sort_values(by = "prop_somatic_SNPs", ascending = False)
    sample_SNP_mutation_freq_df.to_csv("sample_SNP_mutation_freq.tsv", header = True, sep = '\t', index = False)



@click.command()
@click.option('--maf_path', type=click.Path(exists=True), required=True, help='Path to the MAF file.')
@click.option('--somatic_maf', type=click.Path(exists=True), required=True, help='Path to the filtered somatic mutations file.')
def main(maf_path, somatic_maf):
    """
    CLI entry point for assessing contamination between samples using germline and somatic mutations.
    """
    
    maf_df, somatic_maf_df = data_loading(maf_path, somatic_maf)

    print("Running contamination analysis between samples")
    contamination_detection_between_samples(maf_df, somatic_maf_df)

    print("Running general contamination analysis")
    contamination_detection_in_snps(maf_df)



if __name__ == '__main__':

    main()
