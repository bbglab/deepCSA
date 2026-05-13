#!/usr/bin/env python

import json
from typing import Iterable

import click
import numpy as np
import pandas as pd


def _load_json_keys(path: str) -> set[str]:
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    return set(data.keys())


def _find_column(df: pd.DataFrame, candidates: Iterable[str], label: str) -> str:
    lower_map = {col.lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
        lowered = candidate.lower()
        if lowered in lower_map:
            return lower_map[lowered]
    raise ValueError(f"Missing required {label} column. Tried: {', '.join(candidates)}")


def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    if pvals.size == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = np.arange(1, pvals.size + 1)
    adjusted = np.empty_like(pvals, dtype=float)
    adjusted[order] = pvals[order] * pvals.size / ranked
    adjusted_sorted = np.minimum.accumulate(adjusted[order][::-1])[::-1]
    adjusted[order] = adjusted_sorted
    return np.clip(adjusted, 0.0, 1.0)


def _sample_category(sample: str, group_names: set[str], sample_names: set[str]) -> str:
    if sample == "all_samples":
        return "all_samples"
    if sample in group_names:
        return "group"
    if sample in sample_names:
        return "sample"
    return "unknown"


@click.command()
@click.option("--omegas-file", type=click.Path(exists=True), required=True)
@click.option("--samples-json", type=click.Path(exists=True), required=True)
@click.option("--groups-json", type=click.Path(exists=True), required=True)
@click.option("--output", type=click.Path(), required=True)
def main(omegas_file: str, samples_json: str, groups_json: str, output: str) -> None:
    df = pd.read_table(omegas_file)

    gene_col = _find_column(df, ["gene", "GENE"], "gene")
    sample_col = _find_column(df, ["sample", "SAMPLE"], "sample")
    pvalue_col = _find_column(df, ["pvalue", "p_value", "pval", "p-value"], "p-value")

    df[pvalue_col] = pd.to_numeric(df[pvalue_col], errors="coerce")
    df["pvalue_adj"] = np.nan

    sample_names = _load_json_keys(samples_json)
    group_names = _load_json_keys(groups_json) - {"all_samples"}

    sample_categories = df[sample_col].astype(str).map(
        lambda sample: _sample_category(sample, group_names, sample_names)
    )
    region_levels = df[gene_col].astype(str).str.contains("--", na=False).map(
        lambda has_subgenic: "subgenic" if has_subgenic else "gene"
    )

    for sample_group in ("all_samples", "group", "sample"):
        for region_level in ("gene", "subgenic"):
            mask = (sample_categories == sample_group) & (region_levels == region_level)
            valid_mask = mask & df[pvalue_col].notna()
            if not valid_mask.any():
                continue
            adjusted = _benjamini_hochberg(df.loc[valid_mask, pvalue_col].to_numpy())
            df.loc[valid_mask, "pvalue_adj"] = adjusted

    unknown_samples = sample_categories[sample_categories == "unknown"].unique()
    if len(unknown_samples) > 0:
        click.echo(
            "Warning: omega results contain samples not found in groups/samples JSON: "
            f"{', '.join(sorted(unknown_samples))}"
        )

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()
