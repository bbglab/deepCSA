#!/usr/bin/env python

"""
Script for checking contamination between samples.
"""


import click
import matplotlib.pyplot as plt
import pandas as pd
import logging
import seaborn as sns
import operator
from read_utils import custom_na_values
from utils_filter import germline_mask, somatic_mask

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s",
    level=logging.INFO,
    datefmt="%m/%d/%Y %I:%M:%S %p"
)

LOG = logging.getLogger("check_contamination")

# Constants
GERMLINE_LABEL = "Germline Samples"
SOMATIC_LABEL = "Somatic Samples"
OPS = {
    ">": operator.gt,
    "<": operator.lt,
    ">=": operator.ge,
    "<=": operator.le
}

# Assuming somatic_variants and germline_variants are loaded as pandas DataFrames
def compute_shared_variants(somatic_variants, germline_variants):
    """Count mutations shared between each somatic sample and each germline sample.

    Parameters
    ----------
    somatic_variants : pd.DataFrame
        Variant table with at least ``SAMPLE_ID`` and ``MUT_ID`` columns,
        providing the somatic side of the comparison.
    germline_variants : pd.DataFrame
        Variant table with at least ``SAMPLE_ID`` and ``MUT_ID`` columns,
        providing the germline side of the comparison.

    Returns
    -------
    pd.DataFrame
        Integer matrix indexed by somatic ``SAMPLE_ID`` (rows) and germline
        ``SAMPLE_ID`` (columns); each cell is the number of ``MUT_ID`` values
        shared between that pair of samples.
    """

    unique_somatic_samples = sorted(somatic_variants["SAMPLE_ID"].unique())
    unique_germline_samples = sorted(germline_variants["SAMPLE_ID"].unique())

    # Create a DataFrame to store counts (avoid .fillna downcasting warning)
    shared_counts = pd.DataFrame(0, index=unique_somatic_samples, columns=unique_germline_samples, dtype=int)

    # Iterate through each somatic sample
    for somatic_sample in unique_somatic_samples:
        somatic_mutations = set(somatic_variants[somatic_variants["SAMPLE_ID"] == somatic_sample]["MUT_ID"])

        # Compare with germline mutations of all other samples
        for germline_sample in unique_germline_samples:
            germline_mutations = set(germline_variants[germline_variants["SAMPLE_ID"] == germline_sample]["MUT_ID"])

            # Count shared mutations
            shared_counts.loc[somatic_sample, germline_sample] = len(somatic_mutations & germline_mutations)

    return shared_counts

def create_heatmap(variants_matrix: pd.DataFrame,
                   annot: pd.DataFrame,
                   col_labels: list,
                   xlabel: str,
                   ylabel: str,
                   title: str,
                   output_file: str,
                   annot_kws_color: str = "black",
                   size: tuple[int, int] = (18, 15)):
    """Create heatmap for specified set of mutations."""
    plt.figure(figsize=size)
    sns.heatmap(
        variants_matrix,
        annot=annot,
        fmt="",
        cmap="Blues",
        cbar_kws={"label": "Shared Mutations"},
        xticklabels=col_labels,
        yticklabels=variants_matrix.index,
        linewidths=0.5,
        annot_kws={"color": annot_kws_color, "fontsize": 10},
    )

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.savefig(output_file, bbox_inches="tight", dpi=100)
    plt.close()


def prepare_datasets(maf_df, somatic_maf_df, somatic_vaf_boundary):
    """Prepare datasets for contamination analysis.

    Parameters
    ----------
    maf_df : pd.DataFrame
        Full mutation table for all samples (used to derive germline variants and the
        all-variants set), with at least ``SAMPLE_ID``, ``MUT_ID``, ``VAF``, ``vd_VAF`` and
        ``VAF_AM`` columns.
    somatic_maf_df : pd.DataFrame
        Filtered somatic mutation table, with at least ``SAMPLE_ID`` and ``MUT_ID`` columns.
    somatic_vaf_boundary : float
        VAF threshold passed to ``germline_mask`` to identify germline variants (a variant is
        germline when all of ``VAF``/``vd_VAF``/``VAF_AM`` exceed it).

    Returns
    -------
    tuple
        Tuple containing:
        - germline_vars_all_samples: DataFrame of germline variants across all samples.
        - somatic_variants: DataFrame of somatic variants.
        - all_variants: DataFrame of all variants.
    """
    # Consider both unique and non-unique variants when collecting germline variants
    germline_vars_all_samples = maf_df.loc[
        germline_mask(maf_df, somatic_vaf_boundary), ["SAMPLE_ID", "MUT_ID"]
    ].drop_duplicates()
    LOG.info(f"Total variants: {germline_vars_all_samples.shape}")
    LOG.info(f"Unique germline variants: {len(germline_vars_all_samples['MUT_ID'].unique())}")
    
    somatic_variants = somatic_maf_df[["SAMPLE_ID", "MUT_ID"]]
    LOG.info(f"Somatic variants: {somatic_variants.shape}")
    all_variants = maf_df[["SAMPLE_ID", "MUT_ID"]]
    LOG.info(f"All variants: {all_variants.shape}")

    return germline_vars_all_samples, somatic_variants, all_variants


def two_way_comparison(df_a: pd.DataFrame, df_b: pd.DataFrame, annotation_threshold: float, operator_str: str) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Detect cross-sample contamination, compare one set of mutations with another set of mutations.

        Parameters
        ----------
        df_a : pd.DataFrame
            First DataFrame of variants (e.g. somatic).
        df_b : pd.DataFrame
            Second DataFrame of variants (e.g. germline).
        annotation_threshold : float
            Threshold for annotating shared variants.
        operator_str : str
            String representing the comparison operator (e.g., ">", "<", ">=", "<=").
    """
    shared_df = compute_shared_variants(df_a, df_b,)

    # Compute total number of counts of mutations per sample -> determined by the first dataframe (df_a)
    counts = (
        df_a["SAMPLE_ID"].value_counts().reindex(shared_df.columns)
    )

    # Create custom column labels with germline mutation counts
    col_labels = [f"(n={counts[col]}) {col}" for col in shared_df.columns]

    # Build annotation DataFrame without using deprecated DataFrame.applymap
    mask = OPS[operator_str](shared_df, annotation_threshold)
    annot = shared_df.where(mask)
    # convert selected values to nullable int then to string, replace missing with empty string
    annot = annot.round(0).astype("Int64").astype(str).replace("<NA>", "").fillna("")

    return shared_df, annot, col_labels

def contamination_detection_between_samples(maf_df, somatic_maf_df, somatic_vaf_boundary):
    """Detect cross-sample contamination by comparing somatic and germline mutations.

    Parameters
    ----------
    maf_df : pd.DataFrame
        Full mutation table for all samples (used to derive germline variants and the
        all-variants set), with at least ``SAMPLE_ID``, ``MUT_ID``, ``VAF``, ``vd_VAF`` and
        ``VAF_AM`` columns.
    somatic_maf_df : pd.DataFrame
        Filtered somatic mutation table, with at least ``SAMPLE_ID`` and ``MUT_ID`` columns.
    somatic_vaf_boundary : float
        VAF threshold passed to ``germline_mask`` to identify germline variants (a variant is
        germline when all of ``VAF``/``vd_VAF``/``VAF_AM`` exceed it).
    """
    # Prepare datasets
    germline_vars_all_samples, somatic_variants, all_variants = prepare_datasets(maf_df, somatic_maf_df, somatic_vaf_boundary)

    ## Somatic vs Germline
    shared_variants_somatic2germline_matrix, somatic_vs_germline_annot, somatic_vs_germline_col_labels = two_way_comparison(germline_vars_all_samples, somatic_variants, annotation_threshold=0.3, operator_str=">")

    create_heatmap(
        variants_matrix=shared_variants_somatic2germline_matrix,
        annot=somatic_vs_germline_annot,
        col_labels=somatic_vs_germline_col_labels,
        xlabel=GERMLINE_LABEL,
        ylabel=SOMATIC_LABEL,
        title="Somatic mutations that are germline in other samples",
        output_file="somatic_vs_germline.pdf")

    ## Germline vs Germline
    shared_germline_variants_matrix, germline_vs_germline_annot, germline_vs_germline_col_labels = two_way_comparison(germline_vars_all_samples, germline_vars_all_samples, annotation_threshold=0, operator_str="<")

    create_heatmap(
        variants_matrix=shared_germline_variants_matrix,
        annot=germline_vs_germline_annot,
        col_labels=germline_vs_germline_col_labels,
        xlabel=GERMLINE_LABEL,
        ylabel=GERMLINE_LABEL,
        title="Germline mutations that are germline in other samples",
        output_file="germline_vs_germline.pdf")

    ## All vs Germline
    shared_all_vs_germline_variants_matrix = compute_shared_variants(all_variants, germline_vars_all_samples)

    # Compute total number of germline mutations per sample
    germline_counts = (
        germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_all_vs_germline_variants_matrix.columns)
    )

    normalized_shared_all_vs_germline_variants_matrix = shared_all_vs_germline_variants_matrix.divide(
        germline_counts, axis=1
    )

    # Count shared mutations between somatic and germline samples

    # Create custom column labels with germline mutation counts
    col_labels = [
        f"(n={germline_counts[col]}) {col}" for col in normalized_shared_all_vs_germline_variants_matrix.columns
    ]

    # Annotation: show rounded values only when 0.8 < x < 1
    cond = (normalized_shared_all_vs_germline_variants_matrix > 0.8) & (
        normalized_shared_all_vs_germline_variants_matrix < 1
    )
    annot = normalized_shared_all_vs_germline_variants_matrix.where(cond)
    annot = annot.round(2).astype("string").fillna("")

    create_heatmap(
        variants_matrix=normalized_shared_all_vs_germline_variants_matrix,
        annot=annot,
        col_labels=col_labels,
        xlabel=GERMLINE_LABEL,
        ylabel="All mutations samples",
        title="All mutations that are germline in other samples",
        output_file="allmutations_vs_germline.pdf",
        annot_kws_color="white")
    

    # Show total number of germline mutations

    germline_counts = (
        germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_germline_variants_matrix.columns)
    )

    normalized_share_germline_vs_germline_variants_matrix = shared_germline_variants_matrix.divide(
        germline_counts, axis=1
    )

    # Create custom column labels with germline mutation counts
    col_labels = [
        f"(n={germline_counts[col]}) {col}" for col in normalized_share_germline_vs_germline_variants_matrix.columns
    ]

    cond = (normalized_share_germline_vs_germline_variants_matrix > 0.8) & (
        normalized_share_germline_vs_germline_variants_matrix < 1
    )
    annot = normalized_share_germline_vs_germline_variants_matrix.where(cond)
    annot = annot.round(2).astype("string").fillna("")

    create_heatmap(
        variants_matrix=normalized_share_germline_vs_germline_variants_matrix,
        annot=annot,
        col_labels=col_labels,
        xlabel=GERMLINE_LABEL,
        ylabel=GERMLINE_LABEL,
        title="Normalized germline mutations that are germline in other samples",
        output_file="normalized.germline_vs_germline.pdf",
        annot_kws_color="white")

    ## Somatic vs Remaining Germline
    shared_somatic_to_non_shared_germline = shared_all_vs_germline_variants_matrix - shared_germline_variants_matrix

    # Those cases where the number of mutations is smaller than 5 are set to 0
    shared_somatic_to_non_shared_germline[shared_somatic_to_non_shared_germline < 5] = 0

    # Compute total number of germline mutations per sample
    germline_counts = (
        germline_vars_all_samples["SAMPLE_ID"].value_counts().reindex(shared_somatic_to_non_shared_germline.columns)
    )

    total_germline_available_per_sample = germline_counts - shared_germline_variants_matrix

    shared_somatic_to_non_shared_germline_proportion = (
        shared_somatic_to_non_shared_germline / total_germline_available_per_sample
    ).fillna(0)

    cond = shared_somatic_to_non_shared_germline_proportion > 0.45
    annot = shared_somatic_to_non_shared_germline_proportion.where(cond)
    annot = annot.round(2).astype("string").fillna("")

    create_heatmap(
        variants_matrix=shared_somatic_to_non_shared_germline_proportion,
        annot=annot,
        col_labels=col_labels,
        xlabel="Non-shared germline",
        ylabel="Somatic",
        title="Somatic mutations that are germline in other samples",
        output_file="contamination.somatic_vs_remaininggermline.pdf",
        annot_kws_color="white",
        size=(22, 18)
    )
    cond = shared_somatic_to_non_shared_germline > 0
    annot = shared_somatic_to_non_shared_germline.where(cond)
    # convert to nullable int then string, replace missing with empty string
    annot = annot.round(0).astype("Int64").astype("string").replace("<NA>", "").fillna("")

    create_heatmap(
        variants_matrix=shared_somatic_to_non_shared_germline,
        annot=annot,
        col_labels=col_labels,
        xlabel="Non-shared germline",
        ylabel="Somatic",
        title="Somatic mutations that are germline in other samples (count)",
        output_file="contamination.somatic_vs_remaininggermline.numbers.pdf",
        size=(22, 18)
    )

    max_prop_per_sample = shared_somatic_to_non_shared_germline_proportion.max(axis="columns")

    ## Exploration of contaminated samples
    receiver_source_pairs = []
    for sample, max_val in max_prop_per_sample[max_prop_per_sample > 0.5].reset_index().values:
        sample_vals = shared_somatic_to_non_shared_germline_proportion.loc[sample, :]
        sample_vals_count = shared_somatic_to_non_shared_germline.loc[sample, :]

        source_sampleids = sample_vals[sample_vals == max_val].index.values
        source_sampleid = source_sampleids[0]
        receiver_source_pairs.append(
            (
                sample,
                round(max_val, 3),
                list(zip([sample_vals_count[x].item() for x in source_sampleids], source_sampleids)),
            )
        )

        print(
            f"{sample} has {max_val:.2f} proportion of the germline variants of {source_sampleid} as with a VAF not corresponding to germline variants."
        )
        print(f"Shared variants count: {sample_vals_count[source_sampleid]}")
        print()

        subseeeet = maf_df[
            [
                "SAMPLE_ID",
                "MUT_ID",
                "canonical_SYMBOL",
                "ALT_DEPTH",
                "DEPTH",
                "VAF",
                "canonical_Consequence_broader",
                "FILTER",
            ]
        ]
        p_dest = subseeeet[subseeeet["SAMPLE_ID"] == sample].drop("SAMPLE_ID", axis=1)

        p_source_germ = germline_vars_all_samples[germline_vars_all_samples["SAMPLE_ID"] == source_sampleid]
        p_source = subseeeet[
            (subseeeet["SAMPLE_ID"] == source_sampleid) & (subseeeet["MUT_ID"].isin(p_source_germ["MUT_ID"].values))
        ].drop("SAMPLE_ID", axis=1)

        merged_samples = p_dest.merge(
            p_source,
            on=["MUT_ID", "canonical_SYMBOL", "canonical_Consequence_broader"],
            suffixes=("_dest", "_source"),
            how="right",
        )

        merged_samples.sort_values(by=["VAF_dest"], ascending=False).to_csv(
            f"{source_sampleid}.germline_variants_in.{sample}.tsv", header=True, sep="\t", index=False
        )

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
    contamination_detailed_df = pd.DataFrame(
        receiver_source_pairs, columns=["SAMPLE_ID", "MAX_PROPORTION_GERMLINE_FROM_SOURCE", "SOURCE_SAMPLEID_COUNTS"]
    )
    contamination_detailed_df.to_csv("contaminated_samples.detailed.tsv", header=True, sep="\t", index=False)

    if contamination_detailed_df.empty:
        print("No contaminated samples detected.")
        return
    contamination_detailed_df_long = contamination_detailed_df.explode("SOURCE_SAMPLEID_COUNTS")
    expanded_df = pd.DataFrame(contamination_detailed_df_long["SOURCE_SAMPLEID_COUNTS"].tolist())
    expanded_df.columns = ["SHARED_VARIANT_COUNT", "SOURCE_SAMPLEID"]
    contamination_detailed_df_long["SHARED_VARIANT_COUNT"] = expanded_df["SHARED_VARIANT_COUNT"].values
    contamination_detailed_df_long["SOURCE_SAMPLEID"] = expanded_df["SOURCE_SAMPLEID"].values
    contamination_detailed_df_long = contamination_detailed_df_long.drop("SOURCE_SAMPLEID_COUNTS", axis=1)
    contamination_detailed_df_long.to_csv("contaminated_samples.detailed.long.tsv", header=True, sep="\t", index=False)


def data_loading(maf_path, somatic_maf_path):
    """Load the full and somatic MAF tables, keeping only covered SNVs.

    Parameters
    ----------
    maf_path : str
        Path to the full MAF file; rows flagged ``FILTER.not_covered`` are dropped and only
        ``TYPE == "SNV"`` rows are kept.
    somatic_maf_path : str
        Path to the filtered somatic MAF file; only ``TYPE == "SNV"`` rows are kept.

    Returns
    -------
    tuple
        ``(maf_df, somatic_maf_df)`` — the filtered full mutation table and the filtered
        somatic mutation table, both as ``pd.DataFrame``.
    """
    maf_df = pd.read_table(maf_path, na_values=custom_na_values)
    print(maf_df.shape)
    maf_df = maf_df[~(maf_df["FILTER.not_covered"]) & (maf_df["TYPE"] == "SNV")].reset_index()
    print(maf_df.shape)

    somatic_maf_df = pd.read_table(somatic_maf_path, na_values=custom_na_values)
    print(somatic_maf_df.shape)
    somatic_maf_df = somatic_maf_df[(somatic_maf_df["TYPE"] == "SNV")]
    print(somatic_maf_df.shape)
    return maf_df, somatic_maf_df


def contamination_detection_in_snps(maf):
    """Estimate per-sample contamination from the VAF distribution at known SNP positions.

    Restricts to gnomAD SNP positions, splits them into somatic-looking and germline-looking
    sets by VAF, computes the per-sample proportion of SNP positions that look somatic, writes
    the resulting table, and plots its distribution across samples.

    Parameters
    ----------
    maf : pd.DataFrame
        Full mutation table with at least ``SAMPLE_ID``, ``MUT_ID``, ``VAF``, ``vd_VAF``, ``VAF_AM``, and the boolean
        ``FILTER.gnomAD_SNP`` column.
    """
    snp_positions_maf = maf[maf["FILTER.gnomAD_SNP"]][["SAMPLE_ID", "MUT_ID", "VAF", "vd_VAF", "VAF_AM"]].reset_index(
        drop=True
    )

    # being very restrictive in the VAF to count the occurrences of potentially contaminated mutations
    contamination_vaf_threshold = 0.05
    somatic_snp_positions_maf = snp_positions_maf.loc[
        somatic_mask(snp_positions_maf, contamination_vaf_threshold)
    ].reset_index(drop=True)
    germline_snp_positions_maf = snp_positions_maf.loc[
        germline_mask(snp_positions_maf, contamination_vaf_threshold)
    ].reset_index(drop=True)

    unique_SNP_positions = snp_positions_maf["MUT_ID"].unique()
    number_unique_SNP_positions = len(unique_SNP_positions)

    sample_SNP_mutation_freq = []
    for sample in snp_positions_maf["SAMPLE_ID"].unique():
        germline_count = len(germline_snp_positions_maf[germline_snp_positions_maf["SAMPLE_ID"] == sample])
        somatic_count = len(somatic_snp_positions_maf[somatic_snp_positions_maf["SAMPLE_ID"] == sample])
        remaining_germline = number_unique_SNP_positions - germline_count
        sample_SNP_mutation_freq.append(
            [
                sample,
                germline_count,
                remaining_germline,
                somatic_count,
                somatic_count / remaining_germline if remaining_germline > 0 else 1,
            ]
        )
    sample_SNP_mutation_freq_df = pd.DataFrame(sample_SNP_mutation_freq)
    sample_SNP_mutation_freq_df.columns = [
        "SAMPLE_ID",
        "germline_count",
        "remaining_germline",
        "somatic_count",
        "prop_somatic_SNPs",
    ]

    # identify outliers in the "prop_somatic_SNPs" column
    sample_SNP_mutation_freq_df = sample_SNP_mutation_freq_df.sort_values(by="prop_somatic_SNPs", ascending=False)
    sample_SNP_mutation_freq_df.to_csv("sample_SNP_mutation_freq.tsv", header=True, sep="\t", index=False)

    plt.figure(figsize=(6, 3))
    sns.violinplot(data=sample_SNP_mutation_freq_df, x="prop_somatic_SNPs", fill=False, color="lightgray", inner=None)
    sns.swarmplot(data=sample_SNP_mutation_freq_df, x="prop_somatic_SNPs", color="black", size=3)

    plt.title("Proportion of all SNPs across samples\ndetected as somatic")
    plt.xlabel("Proportion of somatic SNPs per sample")
    plt.ylabel("Density")
    plt.savefig("sample_SNP_mutation_freq.pdf", dpi=300, bbox_inches="tight")
    plt.close()


@click.command()
@click.option("--maf_path", type=click.Path(exists=True), required=True, help="Path to the MAF file.")
@click.option(
    "--somatic_maf", type=click.Path(exists=True), required=True, help="Path to the filtered somatic mutations file."
)
@click.option(
    "--somatic-vaf-boundary",
    type=float,
    default=0.3,
    show_default=True,
    help="VAF boundary for somatic variants; a variant with all of VAF/vd_VAF/VAF_AM above it is germline.",
)
def main(maf_path, somatic_maf, somatic_vaf_boundary):
    """Assess cross-sample contamination using germline and somatic mutations.

    Loads the input tables and runs both the between-samples and the SNP-based contamination
    analyses, writing their tables and plots to the current working directory.

    Parameters
    ----------
    maf_path : str
        Path to the full MAF file.
    somatic_maf : str
        Path to the filtered somatic mutations file.
    somatic_vaf_boundary : float, optional
        VAF boundary separating somatic from germline variants. Default is 0.3.
    """

    maf_df, somatic_maf_df = data_loading(maf_path, somatic_maf)

    print("Running contamination analysis between samples")
    contamination_detection_between_samples(maf_df, somatic_maf_df, somatic_vaf_boundary)

    print("Running general contamination analysis")
    contamination_detection_in_snps(maf_df)


if __name__ == "__main__":
    main()
