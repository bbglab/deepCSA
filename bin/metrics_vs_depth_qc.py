#!/usr/bin/env python

import json
from pathlib import Path

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import linregress


def load_group_samples(group_definition_file, group_name):
    with open(group_definition_file, "r") as f:
        groups = json.load(f)
    if group_name not in groups:
        raise ValueError(f"Group '{group_name}' not present in {group_definition_file}")
    return [str(x) for x in groups[group_name]]


def _find_column(columns, candidates):
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def load_depth_table(depth_gene_sample_file, samples):
    depth_df = pd.read_csv(depth_gene_sample_file, sep="\t", header=0)
    sample_col = _find_column(depth_df.columns, ["SAMPLE_ID", "sample"])
    gene_col = _find_column(depth_df.columns, ["GENE", "gene"])
    depth_col = _find_column(depth_df.columns, ["MEAN_GENE_DEPTH", "mean_gene_depth"])

    if sample_col is None or gene_col is None or depth_col is None:
        raise ValueError("Depth table must contain SAMPLE_ID, GENE and MEAN_GENE_DEPTH columns.")

    depth_df = depth_df.rename(columns={sample_col: "SAMPLE_ID", gene_col: "GENE", depth_col: "MEAN_GENE_DEPTH"})
    depth_df["SAMPLE_ID"] = depth_df["SAMPLE_ID"].astype(str)
    depth_df["GENE"] = depth_df["GENE"].astype(str)
    depth_df = depth_df[depth_df["SAMPLE_ID"].isin(samples)].copy().reset_index(drop=True)
    depth_df["MEAN_GENE_DEPTH"] = pd.to_numeric(depth_df["MEAN_GENE_DEPTH"], errors="coerce")
    return depth_df


def load_mutdensity(mutdensity_file, samples, metric_name):
    mut = pd.read_csv(mutdensity_file, sep="\t", header=0)
    required = {"SAMPLE_ID", "GENE", "MUTDENSITY_MB"}
    if not required.issubset(set(mut.columns)):
        raise ValueError(f"Mutdensity file {mutdensity_file} must contain {sorted(required)}")

    mut = mut[(mut["SAMPLE_ID"].astype(str).isin(samples)) & (mut["GENE"] != "ALL_GENES")].copy()
    mut["SAMPLE_ID"] = mut["SAMPLE_ID"].astype(str)
    mut["GENE"] = mut["GENE"].astype(str)
    mut["metric_value"] = pd.to_numeric(mut["MUTDENSITY_MB"], errors="coerce")
    mut["metric_name"] = metric_name
    if "REGIONS" in mut.columns:
        mut["metric_name"] = mut["metric_name"] + "." + mut["REGIONS"].fillna("unknown").astype(str)
    return mut[["SAMPLE_ID", "GENE", "metric_name", "metric_value"]]


def load_omegas(omegas_file, samples):
    omega = pd.read_csv(omegas_file, sep="\t", header=0)
    sample_col = _find_column(omega.columns, ["sample", "SAMPLE_ID"])
    gene_col = _find_column(omega.columns, ["gene", "GENE"])
    dnds_col = _find_column(omega.columns, ["dnds", "omega"])
    if sample_col is None or gene_col is None or dnds_col is None:
        raise ValueError(f"Omega file {omegas_file} must contain sample, gene and dnds/omega columns")

    omega = omega.rename(columns={sample_col: "SAMPLE_ID", gene_col: "GENE", dnds_col: "metric_value"})
    omega["SAMPLE_ID"] = omega["SAMPLE_ID"].astype(str)
    omega["GENE"] = omega["GENE"].astype(str).str.split("--").str[0]
    omega = omega[omega["SAMPLE_ID"].isin(samples)].copy()
    omega["metric_value"] = pd.to_numeric(omega["metric_value"], errors="coerce")
    omega["metric_name"] = "omega_globalloc"
    return omega[["SAMPLE_ID", "GENE", "metric_name", "metric_value"]]


def summarize_effect(df):
    rows = []
    work = df.dropna(subset=["MEAN_GENE_DEPTH", "metric_value"]).copy()
    all_groups = [("all_samples", work)] + [(s, g) for s, g in work.groupby("SAMPLE_ID")]
    for group_name, gdf in all_groups:
        n = len(gdf)
        if n < 3:
            rows.append(
                {
                    "sample_scope": group_name,
                    "n_points": n,
                    "pearson_r": np.nan,
                    "spearman_r": np.nan,
                    "slope_metric_per_depth": np.nan,
                    "linreg_pvalue": np.nan,
                }
            )
            continue
        lr = linregress(gdf["MEAN_GENE_DEPTH"], gdf["metric_value"])
        rows.append(
            {
                "sample_scope": group_name,
                "n_points": n,
                "pearson_r": gdf["MEAN_GENE_DEPTH"].corr(gdf["metric_value"], method="pearson"),
                "spearman_r": gdf["MEAN_GENE_DEPTH"].corr(gdf["metric_value"], method="spearman"),
                "slope_metric_per_depth": lr.slope,
                "linreg_pvalue": lr.pvalue,
            }
        )
    return pd.DataFrame(rows)


def summarize_missingness(depth_df, merged_df):
    m = depth_df.merge(
        merged_df[["SAMPLE_ID", "GENE", "metric_value"]],
        on=["SAMPLE_ID", "GENE"],
        how="left",
    )
    m["is_missing_metric"] = m["metric_value"].isna()
    by_gene = (
        m.groupby("GENE")
        .apply(
            lambda g: pd.Series(
                {
                    "n_samples_total": int(len(g)),
                    "n_missing_metric": int(g["is_missing_metric"].sum()),
                    "missing_fraction": float(g["is_missing_metric"].mean()),
                    "median_depth_missing": g.loc[g["is_missing_metric"], "MEAN_GENE_DEPTH"].median(),
                    "median_depth_nonmissing": g.loc[~g["is_missing_metric"], "MEAN_GENE_DEPTH"].median(),
                    "min_depth_missing": g.loc[g["is_missing_metric"], "MEAN_GENE_DEPTH"].min(),
                    "max_depth_missing": g.loc[g["is_missing_metric"], "MEAN_GENE_DEPTH"].max(),
                }
            )
        )
        .reset_index()
    )
    return by_gene.sort_values(["missing_fraction", "n_missing_metric"], ascending=[False, False]).reset_index(drop=True)


def plot_scatter_per_sample(df, group_name, metric_name, output_pdf):
    samples = sorted(df["SAMPLE_ID"].dropna().unique().tolist())
    if not samples:
        return

    sns.set_style("whitegrid")
    per_page = 12
    ncols = 4
    nrows = 3
    with PdfPages(output_pdf) as pdf:
        for start in range(0, len(samples), per_page):
            page_samples = samples[start : start + per_page]
            fig, axes = plt.subplots(nrows, ncols, figsize=(18, 12), sharex=False, sharey=False)
            axes = axes.flatten()
            for i, sample in enumerate(page_samples):
                ax = axes[i]
                sdf = df[df["SAMPLE_ID"] == sample].dropna(subset=["MEAN_GENE_DEPTH", "metric_value"])
                if sdf.empty:
                    ax.set_title(f"{sample} (no data)")
                    continue
                sns.scatterplot(
                    data=sdf,
                    x="MEAN_GENE_DEPTH",
                    y="metric_value",
                    alpha=0.5,
                    s=16,
                    linewidth=0,
                    ax=ax,
                )
                if len(sdf) >= 3:
                    sns.regplot(
                        data=sdf,
                        x="MEAN_GENE_DEPTH",
                        y="metric_value",
                        scatter=False,
                        line_kws={"color": "darkred", "linewidth": 1.2},
                        ax=ax,
                    )
                ax.set_title(f"{sample} (n={len(sdf)})", fontsize=9)
                ax.set_xlabel("MEAN_GENE_DEPTH")
                ax.set_ylabel(metric_name)

            for j in range(len(page_samples), len(axes)):
                axes[j].axis("off")

            fig.suptitle(f"{group_name} | {metric_name} vs depth (per sample)", fontsize=14)
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig)
            plt.close(fig)


@click.command()
@click.option("--mutdensity-file", required=True, type=click.Path(exists=True))
@click.option("--depth-gene-sample-file", required=True, type=click.Path(exists=True))
@click.option("--group-definition", required=True, type=click.Path(exists=True))
@click.option("--group-name", required=True, type=str)
@click.option("--output-dir", required=True, type=click.Path())
@click.option("--adjusted-mutdensity-file", required=False, type=click.Path(exists=True))
@click.option("--omegas-file", required=False, type=click.Path(exists=True))
def main(
    mutdensity_file,
    depth_gene_sample_file,
    group_definition,
    group_name,
    output_dir,
    adjusted_mutdensity_file,
    omegas_file,
):
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    samples = load_group_samples(group_definition, group_name)
    depth_df = load_depth_table(depth_gene_sample_file, samples)

    metric_frames = [load_mutdensity(mutdensity_file, samples, "mutdensity")]
    if adjusted_mutdensity_file:
        try:
            metric_frames.append(load_mutdensity(adjusted_mutdensity_file, samples, "adjusted_mutdensity"))
        except Exception as e:
            print(f"Warning: skipping adjusted mutdensity file {adjusted_mutdensity_file}: {e}")
    if omegas_file:
        try:
            metric_frames.append(load_omegas(omegas_file, samples))
        except Exception as e:
            print(f"Warning: skipping omegas file {omegas_file}: {e}")

    metrics_df = pd.concat(metric_frames, ignore_index=True)

    status_rows = []
    for metric_name, mdf in metrics_df.groupby("metric_name"):
        merged = depth_df.merge(mdf, on=["SAMPLE_ID", "GENE"], how="left")
        merged["metric_name"] = metric_name
        summary_effect = summarize_effect(merged)
        summary_effect.insert(0, "metric_name", metric_name)
        summary_effect.to_csv(outdir / f"{group_name}.{metric_name}.depth_effect_summary.tsv", sep="\t", index=False)

        missingness = summarize_missingness(depth_df, merged)
        missingness.insert(0, "metric_name", metric_name)
        missingness.to_csv(outdir / f"{group_name}.{metric_name}.depth_missingness_by_gene.tsv", sep="\t", index=False)

        plot_scatter_per_sample(
            merged,
            group_name=group_name,
            metric_name=metric_name,
            output_pdf=outdir / f"{group_name}.{metric_name}.depth_scatter_per_sample.pdf",
        )

        status_rows.append(
            {
                "group_name": group_name,
                "metric_name": metric_name,
                "n_depth_rows": int(len(depth_df)),
                "n_metric_rows": int(len(mdf)),
                "n_merged_nonmissing": int(merged["metric_value"].notna().sum()),
            }
        )

    pd.DataFrame(status_rows).to_csv(outdir / f"{group_name}.metrics_vs_depth_qc.status.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
