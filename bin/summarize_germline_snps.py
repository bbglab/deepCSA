#!/usr/bin/env python

"""
Summarize germline SNPs from the clean and somatic mutation tables.

The germline mutations of each sample are obtained by subtracting the somatic
calls from the clean calls (an anti-join on SAMPLE_ID and MUT_ID). Each germline
variant is then assigned to a VAF-based genotype bin, and a PCA is run on the
samples x variants genotype matrix restricted to the heterozygous-like bins.

Outputs (written to the current working directory):
    {output_prefix}.germline.mutations.tsv        Germline mutations with their genotype bin and pathogenic flag.
    {output_prefix}.pathogenic_snps.tsv           Candidate pathogenic germline SNPs with their annotation.
    {output_prefix}.pathogenic_snps_summary.tsv   Number of candidate pathogenic germline SNPs per sample.
    {output_prefix}.ancestry_inference.tsv        Mean gnomAD population AF profile per sample and most likely ethnic group.
    {output_prefix}.germline_snps_summary.pdf     VAF bins, pathogenic SNPs per sample, ancestry heatmap, PCA elbow and PC1 vs PC2 scatter.
"""

import click
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from utils_filter import germline_mask

# VAF cuts separating the genotype bins: (0, 0.25], (0.25, 0.75], (0.75, 1]
GENOTYPE_BINS = [0, 0.25, 0.75, 1]
# Genotype bins kept for the PCA (heterozygous-like VAF ranges)
GENOTYPE_BINS_FOR_PCA = [1, 2]
NUMBER_OF_COMPONENTS = 6
# Minimum VAF for a germline variant to be considered a SNP at all
MIN_SNP_VAF = 0.25
# gnomAD allele frequency columns used to assess rarity; all the available ones
# must be below the threshold for a variant to be considered rare
GNOMAD_AF_COLUMNS = ["gnomADg_AF", "gnomADe_AF"]
# Population-specific gnomAD allele frequency columns used to infer the most
# likely ethnic group of each sample; the genome columns are preferred and the
# exome ones are used as fallback
GNOMAD_GENOME_POP_COLUMNS = {
    "AFR": "gnomADg_AFR_AF",
    "AMI": "gnomADg_AMI_AF",
    "AMR": "gnomADg_AMR_AF",
    "ASJ": "gnomADg_ASJ_AF",
    "EAS": "gnomADg_EAS_AF",
    "FIN": "gnomADg_FIN_AF",
    "MID": "gnomADg_MID_AF",
    "NFE": "gnomADg_NFE_AF",
    "OTH": "gnomADg_OTH_AF",
    "SAS": "gnomADg_SAS_AF",
}
GNOMAD_EXOME_POP_COLUMNS = {
    "AFR": "gnomADe_AFR_AF",
    "AMR": "gnomADe_AMR_AF",
    "ASJ": "gnomADe_ASJ_AF",
    "EAS": "gnomADe_EAS_AF",
    "FIN": "gnomADe_FIN_AF",
    "NFE": "gnomADe_NFE_AF",
    "OTH": "gnomADe_OTH_AF",
    "SAS": "gnomADe_SAS_AF",
}
# Value of the canonical_Protein_affecting column for protein-affecting variants
# (nonsense, missense and essential splice consequences)
PROTEIN_AFFECTING_VALUE = "protein_affecting"
# Columns reported in the per-sample pathogenic SNPs table
PATHOGENIC_SNP_COLUMNS = [
    "SAMPLE_ID",
    "MUT_ID",
    "canonical_SYMBOL",
    "canonical_Consequence_broader",
    "canonical_Amino_acids",
    "VAF",
    "gnomADg_AF",
    "gnomADe_AF",
]


def subtract_somatic_mutations(clean_mutations, somatic_mutations):
    """
    Remove the somatic mutations from the clean mutations to keep the germline ones.

    Parameters
    ----------
    clean_mutations : pd.DataFrame
        Clean mutation table (all samples), with at least SAMPLE_ID, MUT_ID and VAF.
    somatic_mutations : pd.DataFrame
        Somatic mutation table (all samples), with at least SAMPLE_ID, MUT_ID and VAF.

    Returns
    -------
    pd.DataFrame
        Rows of the clean table whose (SAMPLE_ID, MUT_ID) pair is not present in
        the somatic table.
    """
    merged = clean_mutations.merge(
        somatic_mutations[["SAMPLE_ID", "MUT_ID", "VAF"]],
        on=["SAMPLE_ID", "MUT_ID"],
        suffixes=("", "_somatic"),
        how="outer",
    )
    germline_mutations = merged[~(merged["VAF_somatic"].notnull())].reset_index(drop=True)
    return germline_mutations.drop(columns=["VAF_somatic"])


def assign_genotype_bins(germline_mutations):
    """
    Bin the VAF of each germline mutation into genotype classes.

    Bin 0: VAF <= 0, bin 1: (0, 0.25], bin 2: (0.25, 0.75], bin 3: (0.75, 1],
    bin 4: VAF > 1.
    """
    germline_mutations = germline_mutations.copy()
    germline_mutations["genotype_bin"] = np.digitize(
        germline_mutations["VAF"], right=True, bins=GENOTYPE_BINS
    )
    return germline_mutations


def build_pca_matrix(germline_mutations):
    """
    Build the samples x variants genotype matrix used for the PCA.

    Only the mutations in the heterozygous-like genotype bins are kept; samples
    without a mutation get a 0.

    Returns
    -------
    pd.DataFrame
        Matrix with one row per sample and one column per germline mutation.
    """
    subset = germline_mutations[germline_mutations["genotype_bin"].isin(GENOTYPE_BINS_FOR_PCA)]
    pca_data = subset.pivot_table(
        index="SAMPLE_ID", columns="MUT_ID", values="genotype_bin", aggfunc="first"
    )
    return pca_data.fillna(0)


def run_pca(pca_data, number_of_components=NUMBER_OF_COMPONENTS):
    """
    Scale the genotype matrix and run the PCA on it.

    Parameters
    ----------
    pca_data : pd.DataFrame
        Samples x variants genotype matrix.
    number_of_components : int, optional
        Requested number of principal components. It is capped by the number of
        samples and variants available.

    Returns
    -------
    pca_df : pd.DataFrame
        PC coordinates per sample.
    pca : sklearn.decomposition.PCA
        Fitted PCA object.
    """
    n_components = min(number_of_components, pca_data.shape[0], pca_data.shape[1])
    X_scaled = StandardScaler().fit_transform(pca_data)
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(X_scaled)
    pca_df = pd.DataFrame(
        pca_result,
        columns=[f"PC{i+1}" for i in range(n_components)],
        index=pca_data.index,
    )
    return pca_df, pca


def plot_vaf_bins(germline_mutations, pdf):
    """Plot the VAF distribution of the germline mutations per genotype bin."""
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(data=germline_mutations, x="genotype_bin", y="VAF", showfliers=False, ax=ax)
    sns.stripplot(data=germline_mutations, x="genotype_bin", y="VAF", ax=ax)
    ax.set_xlabel("Genotype bin")
    ax.set_ylabel("VAF")
    ax.set_title("Germline mutations VAF per genotype bin")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_elbow(pca, pdf):
    """Plot the variance explained by each principal component."""
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(pca.explained_variance_ratio_, marker="o", linestyle="-")
    ax.set_ylim(0, pca.explained_variance_ratio_.max() + 0.05)
    ax.set_xlabel("PC")
    ax.set_ylabel("Percent of variance explained")
    ax.set_title("PCA elbow plot")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_pca_scatter(pca_df, pca, pdf):
    """Scatter plot of PC1 vs PC2 with the samples annotated."""
    component_x, component_y = "PC1", "PC2"
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(
        pca_df[component_x],
        pca_df[component_y],
        s=80,
        alpha=0.85,
        edgecolors="k",
        linewidths=0.3,
    )
    for sample in pca_df.index:
        ax.annotate(
            sample,
            (pca_df.loc[sample, component_x], pca_df.loc[sample, component_y]),
            fontsize=6,
            ha="left",
            va="bottom",
        )
    ax.set_xlabel(f"{component_x} ({pca.explained_variance_ratio_[0] * 100:.1f}% variance explained)")
    ax.set_ylabel(f"{component_y} ({pca.explained_variance_ratio_[1] * 100:.1f}% variance explained)")
    ax.set_title("PCA of SNPs within Panel Genes")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def flag_pathogenic_snps(germline_mutations, gnomad_af_threshold):
    """
    Flag candidate pathogenic germline SNPs based on the available annotation.

    A germline SNP is considered a candidate pathogenic variant when it complies
    with the germline criteria (VAF, vd_VAF and VAF_AM above the germline
    threshold), has a VAF greater than 0.25 (otherwise it is not a SNP in any
    sense), is protein-affecting (nonsense, missense or essential splice
    according to the canonical_Protein_affecting column) and rare in the
    population (all the available gnomAD allele frequencies below the
    threshold).

    Parameters
    ----------
    germline_mutations : pd.DataFrame
        Germline mutation table with the VEP annotation columns.
    gnomad_af_threshold : float
        Maximum gnomAD allele frequency for a variant to be considered rare.

    Returns
    -------
    pd.DataFrame
        The input table with an additional IS_PATHOGENIC boolean column.
    """
    germline_mutations = germline_mutations.copy()
    is_protein_affecting = (
        germline_mutations["canonical_Protein_affecting"] == PROTEIN_AFFECTING_VALUE
    )

    available_af_columns = [
        col for col in GNOMAD_AF_COLUMNS if col in germline_mutations.columns
    ]
    if available_af_columns:
        af_values = germline_mutations[available_af_columns].apply(pd.to_numeric, errors="coerce")
        is_rare = (af_values.fillna(0) <= gnomad_af_threshold).all(axis="columns")
    else:
        print("No gnomAD allele frequency columns found; the rarity criterion is skipped")
        is_rare = pd.Series(True, index=germline_mutations.index)

    # The germline criteria and the minimum SNP VAF are hard requirements
    is_snp = germline_mask(germline_mutations, 0) & (germline_mutations["VAF"] > MIN_SNP_VAF)
    germline_mutations["IS_PATHOGENIC"] = is_protein_affecting & is_rare & is_snp
    return germline_mutations


def summarize_pathogenic_snps(germline_mutations):
    """
    Build the per-sample and per-variant tables of candidate pathogenic germline SNPs.

    Returns
    -------
    pathogenic_snps : pd.DataFrame
        One row per candidate pathogenic germline SNP with its main annotation.
    pathogenic_snps_summary : pd.DataFrame
        One row per sample with the number of candidate pathogenic germline SNPs.
    """
    pathogenic_snps = germline_mutations[germline_mutations["IS_PATHOGENIC"]].reset_index(drop=True)
    reported_columns = [col for col in PATHOGENIC_SNP_COLUMNS if col in pathogenic_snps.columns]
    pathogenic_snps = pathogenic_snps[reported_columns]

    pathogenic_snps_summary = (
        pathogenic_snps.groupby("SAMPLE_ID", as_index=False)
        .size()
        .rename(columns={"size": "N_PATHOGENIC_SNPS"})
    )
    return pathogenic_snps, pathogenic_snps_summary


def plot_pathogenic_snps(pathogenic_snps_summary, pdf):
    """Bar plot with the number of candidate pathogenic germline SNPs per sample."""
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(pathogenic_snps_summary)), 6))
    ax.bar(pathogenic_snps_summary["SAMPLE_ID"], pathogenic_snps_summary["N_PATHOGENIC_SNPS"])
    ax.set_xlabel("Sample")
    ax.set_ylabel("Number of candidate pathogenic germline SNPs")
    ax.set_title("Candidate pathogenic germline SNPs per sample")
    ax.set_xticks(range(len(pathogenic_snps_summary)))
    ax.set_xticklabels(pathogenic_snps_summary["SAMPLE_ID"], rotation=90, fontsize=6)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def infer_sample_ancestry(germline_mutations):
    """
    Infer the most likely ethnic group of each sample from the gnomAD
    population-specific allele frequencies.

    For each sample, the mean population-specific gnomAD allele frequency of its
    confident germline SNPs (VAF greater than 0.25) is computed per population;
    the population with the highest mean allele frequency is reported as the
    most likely ethnic group. The rationale is that a sample carrying variants
    that are common in a given population is more likely to come from that
    population. The OTH (other) population is reported but excluded from the
    most likely group assignment since it is not an actual ethnic group.

    Returns
    -------
    ancestry_profiles : pd.DataFrame
        One row per sample with the mean population-specific gnomAD allele
        frequency per population and the MOST_LIKELY_POPULATION column.
    """
    population_columns = {
        population: column
        for population, column in GNOMAD_GENOME_POP_COLUMNS.items()
        if column in germline_mutations.columns
    }
    if not population_columns:
        population_columns = {
            population: column
            for population, column in GNOMAD_EXOME_POP_COLUMNS.items()
            if column in germline_mutations.columns
        }
    if not population_columns:
        print("No gnomAD population-specific allele frequency columns found; skipping the ancestry inference")
        return pd.DataFrame()

    confident_snps = germline_mutations[germline_mutations["VAF"] > MIN_SNP_VAF]
    if confident_snps.empty:
        print("No confident germline SNPs (VAF > 0.25) found; skipping the ancestry inference")
        return pd.DataFrame()

    af_profiles = confident_snps[list(population_columns.values())].apply(pd.to_numeric, errors="coerce")
    ancestry_profiles = (
        af_profiles.assign(SAMPLE_ID=confident_snps["SAMPLE_ID"])
        .groupby("SAMPLE_ID", as_index=False)
        .mean()
        .rename(columns={column: population for population, column in population_columns.items()})
    )

    assignable_populations = [pop for pop in population_columns if pop != "OTH"]
    ancestry_profiles["MOST_LIKELY_POPULATION"] = ancestry_profiles[assignable_populations].idxmax(axis="columns")
    return ancestry_profiles


def plot_ancestry_heatmap(ancestry_profiles, pdf):
    """Heatmap of the mean population-specific gnomAD AF per sample."""
    populations = [col for col in ancestry_profiles.columns if col not in ("SAMPLE_ID", "MOST_LIKELY_POPULATION")]
    heatmap_data = ancestry_profiles.set_index("SAMPLE_ID")[populations]
    fig, ax = plt.subplots(figsize=(1.0 + 0.5 * len(populations), 1.0 + 0.4 * len(heatmap_data)))
    sns.heatmap(heatmap_data, annot=True, fmt=".4f", cmap="viridis", ax=ax, cbar_kws={"label": "Mean gnomAD population AF"})
    ax.set_xlabel("gnomAD population")
    ax.set_ylabel("Sample")
    ax.set_title("Ancestry inference from gnomAD population AFs")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


@click.command()
@click.option(
    "--clean_maf",
    type=click.Path(exists=True),
    required=True,
    help="Path to the clean mutations file (all samples).",
)
@click.option(
    "--somatic_maf",
    type=click.Path(exists=True),
    required=True,
    help="Path to the somatic mutations file (all samples).",
)
@click.option(
    "--output_prefix",
    default="",
    show_default=True,
    help="Prefix for the output files.",
)
@click.option(
    "--gnomad-af-threshold",
    default=0.001,
    show_default=True,
    type=float,
    help="Maximum gnomAD allele frequency for a germline SNP to be considered rare.",
)
def main(clean_maf, somatic_maf, output_prefix, gnomad_af_threshold):
    """Derive the germline mutations and summarize them with VAF bins and a PCA."""
    clean_mutations = pd.read_table(clean_maf)
    somatic_mutations = pd.read_table(somatic_maf)
    print(f"Clean mutations: {clean_mutations.shape}")
    print(f"Somatic mutations: {somatic_mutations.shape}")

    germline_mutations = assign_genotype_bins(
        subtract_somatic_mutations(clean_mutations, somatic_mutations)
    )
    print(f"Germline mutations after subtracting the somatic calls: {germline_mutations.shape}")

    germline_mutations = flag_pathogenic_snps(germline_mutations, gnomad_af_threshold)
    germline_mutations.to_csv(f"{output_prefix}.germline.mutations.tsv", sep="\t", index=False)

    pathogenic_snps, pathogenic_snps_summary = summarize_pathogenic_snps(germline_mutations)
    pathogenic_snps.to_csv(f"{output_prefix}.pathogenic_snps.tsv", sep="\t", index=False)
    pathogenic_snps_summary.to_csv(f"{output_prefix}.pathogenic_snps_summary.tsv", sep="\t", index=False)
    print(f"Candidate pathogenic germline SNPs: {pathogenic_snps.shape[0]}")

    ancestry_profiles = infer_sample_ancestry(germline_mutations)
    if not ancestry_profiles.empty:
        ancestry_profiles.to_csv(f"{output_prefix}.ancestry_inference.tsv", sep="\t", index=False)
        print(f"Most likely population per sample:\n{ancestry_profiles[['SAMPLE_ID', 'MOST_LIKELY_POPULATION']].to_string(index=False)}")

    with PdfPages(f"{output_prefix}.germline_snps_summary.pdf") as pdf:
        if germline_mutations.empty:
            print("No germline mutations left after subtracting the somatic calls; skipping the plots")
            return
        plot_vaf_bins(germline_mutations, pdf)
        if not pathogenic_snps_summary.empty:
            plot_pathogenic_snps(pathogenic_snps_summary, pdf)
        if not ancestry_profiles.empty:
            plot_ancestry_heatmap(ancestry_profiles, pdf)

        pca_data = build_pca_matrix(germline_mutations)
        if pca_data.empty:
            print("No germline mutations in the genotype bins used for the PCA; skipping the PCA")
            return
        pca_df, pca = run_pca(pca_data)
        plot_elbow(pca, pdf)
        if pca_df.shape[1] >= 2:
            plot_pca_scatter(pca_df, pca, pdf)
        else:
            print("Only one principal component available; skipping the PC1 vs PC2 scatter")


if __name__ == "__main__":
    main()
