#!/usr/bin/env python
"""
Annotate omegas file with flagged gene/sample pairs.

Takes an omegas TSV and one or more compiled flagged-cases TSVs and adds two
columns to the omegas file: `flagged` (True/False) and `flag_reason` (string or
empty).

"""

import sys
import click
from pathlib import Path
import pandas as pd
from typing import List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns



def load_flagged_tables(paths: List[Path]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read flagged files and return a normalized DataFrame with possible columns:
    SAMPLE, GENE, reason_exclusion
    """
    pieces = []
    for p in paths:
        if not p.exists():
            print(f"Warning: flagged file {p} does not exist, skipping.")
            continue        
        df = pd.read_csv(p, sep="\t", header=0, dtype=str)
        pieces.append(df)
    if not pieces:
        return pd.DataFrame(columns=["sample", "gene", "reason_exclusion"]), pd.DataFrame(columns=["sample", "gene", "reason_exclusion"])

    all_df = pd.concat(pieces, ignore_index=True, sort=False)

    # These are the columns that all_df will contain:
    # ID      reason_exclusion        cohort  regions criteria

    all_df = all_df.rename(columns={
        "cohort": "sample",
        "ID": "gene",
    })
    syn_df = all_df[all_df["regions"] == "synonymous"][["sample", "gene", "reason_exclusion"]].reset_index(drop=True)
    npa_df = all_df[all_df["regions"] == "non_protein_affecting"][["sample", "gene", "reason_exclusion"]].reset_index(drop=True)
    return syn_df, npa_df


def parse_flag_pairs_from_id(id_str: str) -> Tuple[Optional[str], Optional[str]]:
    """Try to extract (sample, gene) or (gene, sample) from id string with separators.
    Returns (sample, gene) where either may be None.
    """
    for sep in [":", "|", ",", "-"]:
        if sep in id_str:
            parts = [p.strip() for p in id_str.split(sep) if p.strip()]
            if len(parts) == 2:
                # Try both orders: (sample,gene) and (gene,sample). We'll return both
                return parts[0], parts[1]
    return None, None


def plot_flagged_summary(syn_flagged: pd.DataFrame, npa_flagged: pd.DataFrame, output_prefix: str = 'flagged', top_n: int = 50) -> None:
    """Plot stacked bar summaries of flagged genes for synonymous and non-protein-affecting.

    Uses the gene order from `syn_flagged` (most frequently flagged first) to
    order the bars. Bars are stacked by `reason_exclusion` and colors for a
    small set of known reasons are defined a priori; any additional reasons
    will get colors from the seaborn color palette.
    """
    # Ensure dataframes exist
    syn = syn_flagged.copy() if syn_flagged is not None else pd.DataFrame(columns=['sample', 'gene', 'reason_exclusion'])
    npa = npa_flagged.copy() if npa_flagged is not None else pd.DataFrame(columns=['sample', 'gene', 'reason_exclusion'])

    # Determine gene order from synonymously flagged genes (descending frequency)
    syn_order = syn['gene'].value_counts().index.tolist()

    # Include any genes present in npa but not in syn at the end
    npa_genes = [g for g in npa['gene'].unique() if g not in syn_order]
    order = syn_order + npa_genes

    # Count occurrences per gene x reason
    syn_counts = syn.groupby(['gene', 'reason_exclusion']).size().unstack(fill_value=0)
    npa_counts = npa.groupby(['gene', 'reason_exclusion']).size().unstack(fill_value=0)

    # Reindex to the chosen order; keep only top_n genes by total syn+npa counts
    total_counts = (syn.groupby('gene').size().add(npa.groupby('gene').size(), fill_value=0)).sort_values(ascending=False)
    top_genes = total_counts.head(top_n).index.tolist()
    order = [g for g in order if g in top_genes]

    syn_counts = syn_counts.reindex(order, fill_value=0)
    npa_counts = npa_counts.reindex(order, fill_value=0)

    # Colors: predefined mapping for known reasons
    predefined = {
        'mutdensity = 0': '#8c8c8c',
        'high_mutdensity - zscore > 2': '#d62728',
        'low_mutdensity - zscore < -2': '#1f77b4',
    }

    # Collect all reasons in consistent order (predefined first)
    reasons = []
    for r in predefined.keys():
        if r in set(list(syn_counts.columns) + list(npa_counts.columns)):
            reasons.append(r)
    # add any other reasons found
    other_reasons = [r for r in sorted(set(list(syn_counts.columns) + list(npa_counts.columns))) if r not in reasons]
    reasons.extend(other_reasons)

    # Build color list: use predefined colors then palette for others
    colors = []
    palette = sns.color_palette('tab10', n_colors=max(3, len(other_reasons)))
    other_iter = iter(palette)
    for r in reasons:
        if r in predefined:
            colors.append(predefined[r])
        else:
            colors.append(next(other_iter))

    # Plotting stacked horizontal barplots (two subplots stacked vertically)
    fig, axes = plt.subplots(2, 1, figsize=(10, max(6, 0.25 * len(order) * 2)), sharex=True)

    # Synonymous stacked bar
    if not syn_counts.empty:
        syn_counts[reasons].plot(kind='barh', stacked=True, color=colors, ax=axes[0])
        axes[0].set_title('Synonymous flagged cases (by gene)')
        axes[0].set_xlabel('Count')
    else:
        axes[0].text(0.5, 0.5, 'No synonymous flagged cases', ha='center')

    # NPA stacked bar
    if not npa_counts.empty:
        npa_counts[reasons].plot(kind='barh', stacked=True, color=colors, ax=axes[1])
        axes[1].set_title('Non-protein-affecting flagged cases (by gene)')
        axes[1].set_xlabel('Count')
    else:
        axes[1].text(0.5, 0.5, 'No non-protein-affecting flagged cases', ha='center')

    plt.tight_layout()
    out_png = f"{output_prefix}_cases_summary.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    # Also write TSV summary
    summary = pd.DataFrame({'gene': order})
    for r in reasons:
        summary[r] = syn_counts.get(r, 0).reindex(order, fill_value=0).values + npa_counts.get(r, 0).reindex(order, fill_value=0).values
    summary.to_csv(f"{output_prefix}_cases_counts.flagged.tsv", sep='\t', index=False)
    print(f"Wrote flagged summary plot {out_png} and table {output_prefix}_cases_counts.flagged.tsv")


def annotate(omegas: pd.DataFrame, flagged: pd.DataFrame) -> pd.DataFrame:
    """Annotate omegas DataFrame with flagged info.

    The returned DataFrame includes two new columns: `flagged` (bool) and
    `flag_reason` (string).
    """
    # omegas has these columns:
    # gene	sample	impact	mutations	dnds	pvalue	lower	upper

    # and flagged has:
    # sample, gene, reason_exclusion


    annotated_omegas = omegas.merge(flagged, how="left",
                                   on=["sample", "gene"])
    annotated_omegas["flagged"] = annotated_omegas["reason_exclusion"].notnull()
    annotated_omegas["flag_reason"] = annotated_omegas["reason_exclusion"].fillna("")

    return annotated_omegas.drop(columns=["reason_exclusion"])





@click.command()
@click.option('--omegas-file', required=True, type=click.Path(exists=True), help='TSV file with omegas; must contain SAMPLE and GENE columns')
@click.option('--compiled-flagged-files', required=True, type=click.Path(exists=True),
              help='Path to a text file that lists one flagged-case TSV file per line (each path may be absolute or relative).')
@click.option('--output', default='annotated_omegas.tsv', type=click.Path(), help='Output annotated TSV file')
def main(omegas_file: str, compiled_flagged_files: str, output: str) -> None:
    """Annotate omegas with flagged gene/sample pairs (CLI entrypoint using click)."""

    omegas_path = Path(omegas_file)

    # `compiled_flagged_files` must be a path to a text file listing one flagged
    # file path per line. Read and collect all non-empty lines.
    cf_path = Path(compiled_flagged_files)
    with cf_path.open() as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    flagged_paths = [Path(l) for l in lines]

    # Read omegas
    omegas = pd.read_csv(omegas_path, sep="\t", header=0, dtype=str).fillna("")

    syn_flagged, npa_flagged = load_flagged_tables(flagged_paths)

    syn_flagged.to_csv("debug.syn_flagged.tsv", sep="\t", index=False)
    npa_flagged.to_csv("debug.npa_flagged.tsv", sep="\t", index=False)

    plot_flagged_summary(syn_flagged, npa_flagged)

    annotated = annotate(omegas, syn_flagged)

    annotated.to_csv(output, sep="\t", index=False)
    print(f"Wrote annotated omegas to {output}")


if __name__ == "__main__":
    raise SystemExit(main())
