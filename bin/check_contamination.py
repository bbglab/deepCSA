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
CONTAMINATION_PROPORTION_THRESHOLD = 0.5
# Deliberately more restrictive than the --somatic-vaf-boundary used for the between-samples
# analysis: this is the numerator of a QC metric, so only confidently low-VAF calls should count.
SNP_CONTAMINATION_VAF_THRESHOLD = 0.05
VAF_COLUMNS = ["VAF", "vd_VAF", "VAF_AM"]
VARIANT_DETAIL_COLUMNS = [
    "SAMPLE_ID",
    "MUT_ID",
    "canonical_SYMBOL",
    "ALT_DEPTH",
    "DEPTH",
    "VAF",
    "canonical_Consequence_broader",
    "FILTER",
]
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
                   col_labels: list | None,
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
        xticklabels=col_labels if col_labels is not None else False,
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


def two_way_comparison(df_a: pd.DataFrame, df_b: pd.DataFrame, annotation_threshold: float, operator_str: str, normalize: bool) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Detect cross-sample contamination, compare one set of mutations with another set of mutations.

    Parameters
    ----------
        df_a : pd.DataFrame
            First DataFrame of variants (e.g. somatic); its samples become the rows of the matrix.
        df_b : pd.DataFrame
            Second DataFrame of variants (e.g. germline); its samples become the columns of the
            matrix and provide the per-sample counts used for the labels and the normalization.
        annotation_threshold : float
            Threshold for annotating shared variants. Ignored when ``normalize`` is True, where
            cells are annotated when the proportion falls between 0.8 and 1.
        operator_str : str
            String representing the comparison operator (e.g., ">", "<", ">=", "<=").
        normalize : bool
            Whether to divide the shared counts by the number of variants of each df_b sample.
    """
    shared_df = compute_shared_variants(df_a, df_b)

    # Compute total number of counts of mutations per sample -> the columns of shared_df are the samples of df_b
    counts = (
        df_b["SAMPLE_ID"].value_counts().reindex(shared_df.columns)
    )

    if normalize:
        shared_df = shared_df.divide(
        counts, axis=1
    )

    # Create custom column labels with mutation counts
    col_labels = [f"(n={counts[col]}) {col}" for col in shared_df.columns]

    # Build annotation DataFrame without using deprecated DataFrame.applymap
    if normalize:
        mask = (shared_df > 0.8) & (
        shared_df < 1)
    else:
        mask = OPS[operator_str](shared_df, annotation_threshold)

    annot = shared_df.where(mask)
    if normalize:
        annot = annot.round(2).astype("string").fillna("")
    else:
        # convert selected values to nullable int then to string, replace missing with empty string
        annot = annot.round(0).astype("Int64").astype("string").replace("<NA>", "").fillna("")

    return shared_df, annot, col_labels

def find_contaminated_pairs(proportion_matrix: pd.DataFrame, counts_matrix: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Identify receiver samples carrying the germline variants of another sample.

    A sample is flagged as a receiver when its highest proportion of non-shared germline
    variants coming from any single other sample exceeds ``threshold``. All sources tied at
    that maximum are reported.

    Parameters
    ----------
    proportion_matrix : pd.DataFrame
        Matrix of the proportion of each source sample's non-shared germline variants that are
        present as non-germline variants in each receiver sample (receivers as rows, sources as
        columns).
    counts_matrix : pd.DataFrame
        Matrix of the corresponding raw shared variant counts, with the same shape and labels.
    threshold : float
        Minimum proportion above which a receiver sample is considered contaminated.

    Returns
    -------
    pd.DataFrame
        One row per contaminated receiver, with columns ``SAMPLE_ID``,
        ``MAX_PROPORTION_GERMLINE_FROM_SOURCE`` and ``SOURCE_SAMPLEID_COUNTS``, the latter
        holding the list of ``(shared_variant_count, source_sample_id)`` pairs tied at the
        maximum proportion.
    """
    max_prop_per_sample = proportion_matrix.max(axis="columns")

    receiver_source_pairs = []
    for sample, max_val in max_prop_per_sample[max_prop_per_sample > threshold].items():
        sample_vals = proportion_matrix.loc[sample, :]
        sample_vals_count = counts_matrix.loc[sample, :]

        source_sampleids = sample_vals[sample_vals == max_val].index.values
        receiver_source_pairs.append(
            (
                sample,
                round(max_val, 3),
                list(zip([sample_vals_count[x].item() for x in source_sampleids], source_sampleids)),
            )
        )

        LOG.info(
            f"{sample} has {max_val:.2f} proportion of the germline variants of {source_sampleids[0]} "
            f"with a VAF not corresponding to germline variants "
            f"(shared variants count: {sample_vals_count[source_sampleids[0]]})."
        )

    return pd.DataFrame(
        receiver_source_pairs, columns=["SAMPLE_ID", "MAX_PROPORTION_GERMLINE_FROM_SOURCE", "SOURCE_SAMPLEID_COUNTS"]
    )


def contaminated_pairs_to_long(contaminated_pairs: pd.DataFrame) -> pd.DataFrame:
    """Expand the tied-source lists into one row per receiver/source pair.

    Parameters
    ----------
    contaminated_pairs : pd.DataFrame
        Output of ``find_contaminated_pairs``; must be non-empty.

    Returns
    -------
    pd.DataFrame
        Long-format table with columns ``SAMPLE_ID``, ``MAX_PROPORTION_GERMLINE_FROM_SOURCE``,
        ``SHARED_VARIANT_COUNT`` and ``SOURCE_SAMPLEID``.
    """
    contaminated_pairs_long = contaminated_pairs.explode("SOURCE_SAMPLEID_COUNTS")
    contaminated_pairs_long[["SHARED_VARIANT_COUNT", "SOURCE_SAMPLEID"]] = pd.DataFrame(
        contaminated_pairs_long["SOURCE_SAMPLEID_COUNTS"].tolist(), index=contaminated_pairs_long.index
    )

    return contaminated_pairs_long.drop(columns="SOURCE_SAMPLEID_COUNTS")


def export_germline_variants_in_receiver(maf_df: pd.DataFrame,
                                         germline_vars_all_samples: pd.DataFrame,
                                         receiver_sample: str,
                                         source_sample: str):
    """Write the source sample's germline variants alongside their VAF in the receiver sample.

    Parameters
    ----------
    maf_df : pd.DataFrame
        Full mutation table for all samples, with at least the ``VARIANT_DETAIL_COLUMNS``.
    germline_vars_all_samples : pd.DataFrame
        Germline variants across all samples, with at least ``SAMPLE_ID`` and ``MUT_ID`` columns.
    receiver_sample : str
        Sample suspected of being contaminated.
    source_sample : str
        Sample suspected of being the contamination source.
    """
    variant_details = maf_df[VARIANT_DETAIL_COLUMNS]

    receiver_variants = variant_details[variant_details["SAMPLE_ID"] == receiver_sample].drop("SAMPLE_ID", axis=1)

    source_germline = germline_vars_all_samples[germline_vars_all_samples["SAMPLE_ID"] == source_sample]
    source_variants = variant_details[
        (variant_details["SAMPLE_ID"] == source_sample) & (variant_details["MUT_ID"].isin(source_germline["MUT_ID"].values))
    ].drop("SAMPLE_ID", axis=1)

    merged_samples = receiver_variants.merge(
        source_variants,
        on=["MUT_ID", "canonical_SYMBOL", "canonical_Consequence_broader"],
        suffixes=("_dest", "_source"),
        how="right",
    )

    merged_samples.sort_values(by=["VAF_dest"], ascending=False).to_csv(
        f"{source_sample}.germline_variants_in.{receiver_sample}.tsv", header=True, sep="\t", index=False
    )


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
    shared_variants_somatic2germline_matrix, somatic_vs_germline_annot, somatic_vs_germline_col_labels = two_way_comparison(somatic_variants, germline_vars_all_samples, annotation_threshold=30, operator_str=">", normalize=False)

    create_heatmap(
        variants_matrix=shared_variants_somatic2germline_matrix,
        annot=somatic_vs_germline_annot,
        col_labels=somatic_vs_germline_col_labels,
        xlabel=GERMLINE_LABEL,
        ylabel=SOMATIC_LABEL,
        title="Somatic mutations that are germline in other samples",
        output_file="somatic_vs_germline.pdf")

    ## Germline vs Germline
    shared_germline_variants_matrix, germline_vs_germline_annot, germline_vs_germline_col_labels = two_way_comparison(germline_vars_all_samples, germline_vars_all_samples, annotation_threshold=0, operator_str="<", normalize=False)

    create_heatmap(
        variants_matrix=shared_germline_variants_matrix,
        annot=germline_vs_germline_annot,
        col_labels=germline_vs_germline_col_labels,
        xlabel=GERMLINE_LABEL,
        ylabel=GERMLINE_LABEL,
        title="Germline mutations that are germline in other samples",
        output_file="germline_vs_germline.pdf")

    ## All vs Germline
    normalized_shared_all_vs_germline_variants_matrix, all_vs_germline_annot, all_vs_germline_col_labels = two_way_comparison(all_variants, germline_vars_all_samples, annotation_threshold=0.3, operator_str=">", normalize=True)

    create_heatmap(
        variants_matrix=normalized_shared_all_vs_germline_variants_matrix,
        annot=all_vs_germline_annot,
        col_labels=all_vs_germline_col_labels,
        xlabel=GERMLINE_LABEL,
        ylabel="All mutations samples",
        title="All mutations that are germline in other samples",
        output_file="allmutations_vs_germline.pdf",
        annot_kws_color="white")
    

    # Germline mutations
    normalized_germline_vs_germline_variants_matrix, norm_germline_vs_germline_annot, norm_germline_vs_germline_col_labels = two_way_comparison(germline_vars_all_samples, germline_vars_all_samples, annotation_threshold=0, operator_str="<", normalize=True)

    create_heatmap(
        variants_matrix=normalized_germline_vs_germline_variants_matrix,
        annot=norm_germline_vs_germline_annot,
        col_labels=norm_germline_vs_germline_col_labels,
        xlabel=GERMLINE_LABEL,
        ylabel=GERMLINE_LABEL,
        title="Normalized germline mutations that are germline in other samples",
        output_file="normalized.germline_vs_germline.pdf",
        annot_kws_color="white")

    ## Somatic vs Remaining Germline
    shared_all_vs_germline_variants_matrix = compute_shared_variants(all_variants, germline_vars_all_samples)
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
        col_labels=None,
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
        col_labels=None,
        xlabel="Non-shared germline",
        ylabel="Somatic",
        title="Somatic mutations that are germline in other samples (count)",
        output_file="contamination.somatic_vs_remaininggermline.numbers.pdf",
        size=(22, 18)
    )

    ## Exploration of contaminated samples
    contaminated_pairs = find_contaminated_pairs(
        shared_somatic_to_non_shared_germline_proportion,
        shared_somatic_to_non_shared_germline,
        CONTAMINATION_PROPORTION_THRESHOLD,
    )
    contaminated_pairs.to_csv("contaminated_samples.detailed.tsv", header=True, sep="\t", index=False)

    if contaminated_pairs.empty:
        LOG.info("No contaminated samples detected.")
        return

    for receiver in contaminated_pairs.itertuples():
        # Only the first of the sources tied at the maximum proportion is reported in detail
        source_sampleid = receiver.SOURCE_SAMPLEID_COUNTS[0][1]
        export_germline_variants_in_receiver(maf_df, germline_vars_all_samples, receiver.SAMPLE_ID, source_sampleid)

    contaminated_pairs_to_long(contaminated_pairs).to_csv(
        "contaminated_samples.detailed.long.tsv", header=True, sep="\t", index=False
    )


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


def prepare_snp_datasets(maf: pd.DataFrame, vaf_threshold: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Split the known SNP positions into somatic-looking and germline-looking variants.

    Germline is defined as the complement of somatic rather than through ``germline_mask``, so
    that the two sets partition the SNP rows. A variant whose VAF estimates disagree (some at or
    below the threshold, some above) is therefore counted as germline: it is not a confident
    somatic call, so it must not be left out of both sets and silently inflate the denominator of
    ``prop_somatic_SNPs``.

    Rows with an undefined VAF cannot be classified either way and are dropped, since
    ``somatic_mask`` is False for them and the germline complement would otherwise absorb them.

    Parameters
    ----------
    maf : pd.DataFrame
        Full mutation table with at least ``SAMPLE_ID``, ``MUT_ID``, the ``VAF_COLUMNS`` and the
        boolean ``FILTER.gnomAD_SNP`` column.
    vaf_threshold : float
        Upper bound (inclusive) on all VAF estimates for a variant to be called somatic.

    Returns
    -------
    tuple
        Tuple containing:
        - snp_positions_maf: DataFrame of all classifiable variants at known SNP positions.
        - somatic_snp_positions_maf: DataFrame of the somatic ones.
        - germline_snp_positions_maf: DataFrame of the remaining ones.
    """
    snp_positions_maf = maf.loc[maf["FILTER.gnomAD_SNP"], ["SAMPLE_ID", "MUT_ID", *VAF_COLUMNS]]

    undefined_vaf = snp_positions_maf[VAF_COLUMNS].isna().any(axis="columns")
    if undefined_vaf.any():
        LOG.warning(f"Discarding {undefined_vaf.sum()} SNP variants with an undefined {VAF_COLUMNS} VAF.")
    snp_positions_maf = snp_positions_maf.loc[~undefined_vaf]
    LOG.info(f"SNP variants: {snp_positions_maf.shape}")

    is_somatic = somatic_mask(snp_positions_maf, vaf_threshold)
    somatic_snp_positions_maf = snp_positions_maf.loc[is_somatic]
    germline_snp_positions_maf = snp_positions_maf.loc[~is_somatic]
    LOG.info(f"Somatic SNP variants: {somatic_snp_positions_maf.shape}")
    LOG.info(f"Germline SNP variants: {germline_snp_positions_maf.shape}")

    return snp_positions_maf, somatic_snp_positions_maf, germline_snp_positions_maf


def compute_snp_somatic_proportion(snp_positions_maf: pd.DataFrame,
                                   somatic_snp_positions_maf: pd.DataFrame,
                                   germline_snp_positions_maf: pd.DataFrame) -> pd.DataFrame:
    """Compute, per sample, the proportion of non-germline SNP positions that look somatic.

    Parameters
    ----------
    snp_positions_maf : pd.DataFrame
        All classifiable variants at known SNP positions, as returned by ``prepare_snp_datasets``.
    somatic_snp_positions_maf : pd.DataFrame
        The somatic subset of ``snp_positions_maf``.
    germline_snp_positions_maf : pd.DataFrame
        The germline subset of ``snp_positions_maf``.

    Returns
    -------
    pd.DataFrame
        One row per sample with columns ``SAMPLE_ID``, ``germline_count``, ``remaining_germline``,
        ``somatic_count`` and ``prop_somatic_SNPs``, sorted by the latter in descending order.
    """
    samples = snp_positions_maf["SAMPLE_ID"].unique()
    number_unique_snp_positions = snp_positions_maf["MUT_ID"].nunique()

    germline_count = germline_snp_positions_maf["SAMPLE_ID"].value_counts().reindex(samples, fill_value=0)
    somatic_count = somatic_snp_positions_maf["SAMPLE_ID"].value_counts().reindex(samples, fill_value=0)

    # NOTE: the denominator counts every SNP position of the cohort that is not germline in this
    # sample, including positions where the sample has no variant at all. Those can never end up in
    # the numerator, so `prop_somatic_SNPs` is driven as much by how many SNP positions were called
    # in a sample as by how many of them look somatic. It is comparable across samples of a cohort
    # sequenced with the same panel, not an absolute contamination rate.
    remaining_germline = number_unique_snp_positions - germline_count

    sample_snp_mutation_freq_df = pd.DataFrame(
        {
            "germline_count": germline_count,
            "remaining_germline": remaining_germline,
            "somatic_count": somatic_count,
            "prop_somatic_SNPs": (somatic_count / remaining_germline).where(remaining_germline > 0, 1),
        }
    ).rename_axis("SAMPLE_ID").reset_index()

    return sample_snp_mutation_freq_df.sort_values(by="prop_somatic_SNPs", ascending=False)


def create_snp_proportion_plot(sample_snp_mutation_freq_df: pd.DataFrame, output_file: str):
    """Plot the distribution of the per-sample proportion of somatic SNPs across the cohort.

    Parameters
    ----------
    sample_snp_mutation_freq_df : pd.DataFrame
        Per-sample table with a ``prop_somatic_SNPs`` column, as returned by
        ``compute_snp_somatic_proportion``.
    output_file : str
        Path of the plot to write.
    """
    plt.figure(figsize=(6, 3))
    sns.violinplot(data=sample_snp_mutation_freq_df, x="prop_somatic_SNPs", fill=False, color="lightgray", inner=None)
    sns.swarmplot(data=sample_snp_mutation_freq_df, x="prop_somatic_SNPs", color="black", size=3)

    plt.title("Proportion of all SNPs across samples\ndetected as somatic")
    plt.xlabel("Proportion of somatic SNPs per sample")
    plt.ylabel("Density")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


def contamination_detection_in_snps(maf):
    """Estimate per-sample contamination from the VAF distribution at known SNP positions.

    Restricts to gnomAD SNP positions, splits them into somatic-looking and germline-looking
    sets by VAF, computes the per-sample proportion of SNP positions that look somatic, writes
    the resulting table, and plots its distribution across samples.

    Parameters
    ----------
    maf : pd.DataFrame
        Full mutation table with at least ``SAMPLE_ID``, ``MUT_ID``, the ``VAF_COLUMNS`` and the
        boolean ``FILTER.gnomAD_SNP`` column.
    """
    snp_positions_maf, somatic_snp_positions_maf, germline_snp_positions_maf = prepare_snp_datasets(
        maf, SNP_CONTAMINATION_VAF_THRESHOLD
    )

    sample_snp_mutation_freq_df = compute_snp_somatic_proportion(
        snp_positions_maf, somatic_snp_positions_maf, germline_snp_positions_maf
    )
    sample_snp_mutation_freq_df.to_csv("sample_SNP_mutation_freq.tsv", header=True, sep="\t", index=False)

    create_snp_proportion_plot(sample_snp_mutation_freq_df, "sample_SNP_mutation_freq.pdf")


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
