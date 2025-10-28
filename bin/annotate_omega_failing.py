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


def plot_sample_flagged_summary(syn_flagged: pd.DataFrame, npa_flagged: pd.DataFrame, output_prefix: str = 'flagged_samples', top_n: int = 50) -> None:
    """Plot stacked bar summaries of how many genes per sample fail QC.

    Uses sample order derived from total flagged genes (descending) and stacks
    by reason_exclusion. Output mirror of plot_flagged_summary but across
    samples/cohorts.
    """
    syn = syn_flagged.copy() if syn_flagged is not None else pd.DataFrame(columns=['sample', 'gene', 'reason_exclusion'])
    npa = npa_flagged.copy() if npa_flagged is not None else pd.DataFrame(columns=['sample', 'gene', 'reason_exclusion'])

    # For each sample and reason count unique genes flagged
    syn_uniq = syn.groupby(['sample', 'reason_exclusion'])['gene'].nunique().unstack(fill_value=0)
    npa_uniq = npa.groupby(['sample', 'reason_exclusion'])['gene'].nunique().unstack(fill_value=0)

    total_per_sample = syn.groupby('sample')['gene'].nunique().add(npa.groupby('sample')['gene'].nunique(), fill_value=0)
    top_samples = total_per_sample.sort_values(ascending=False).head(top_n).index.tolist()

    syn_uniq = syn_uniq.reindex(top_samples, fill_value=0)
    npa_uniq = npa_uniq.reindex(top_samples, fill_value=0)

    # reasons and colors reuse logic from gene plot
    predefined = {
        'mutdensity = 0': '#8c8c8c',
        'high_mutdensity - zscore > 2': '#d62728',
        'low_mutdensity - zscore < -2': '#1f77b4',
    }
    reasons = []
    for r in predefined.keys():
        if r in set(list(syn_uniq.columns) + list(npa_uniq.columns)):
            reasons.append(r)
    other_reasons = [r for r in sorted(set(list(syn_uniq.columns) + list(npa_uniq.columns))) if r not in reasons]
    reasons.extend(other_reasons)
    palette = sns.color_palette('tab10', n_colors=max(3, len(other_reasons)))
    other_iter = iter(palette)
    colors = [predefined[r] if r in predefined else next(other_iter) for r in reasons]

    # Build combined dataframe per sample (sum syn+npa reasons)
    combined = pd.DataFrame(index=top_samples)
    for r in reasons:
        combined[r] = syn_uniq.get(r, 0).reindex(top_samples, fill_value=0).values + npa_uniq.get(r, 0).reindex(top_samples, fill_value=0).values

    if combined.empty:
        print('No flagged samples to summarize')
        return

    fig, ax = plt.subplots(figsize=(10, max(6, 0.3 * len(top_samples))))
    combined[reasons].plot(kind='barh', stacked=True, color=colors, ax=ax)
    ax.set_title('Top samples by number of flagged genes')
    ax.set_xlabel('Number of unique genes failing QC')
    plt.tight_layout()
    out_png = f"{output_prefix}_samples_summary.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    combined.reset_index().to_csv(f"{output_prefix}_samples_counts.tsv", sep='\t', index=False)
    print(f"Wrote sample flagged summary plot {out_png} and table {output_prefix}_samples_counts.tsv")


def plot_flagged_heatmap(syn_flagged: pd.DataFrame, npa_flagged: pd.DataFrame, output_prefix: str = 'flagged_heatmap', top_n_genes: int = 200, top_n_samples: int = 200) -> None:
    """Create a genes x samples heatmap showing failure reasons per cell.

    Cells are blank/white when no failure; otherwise colored by reason. If
    multiple reasons exist for the same (gene,sample) we choose the first in
    predefined priority, otherwise pick one deterministically.
    """
    syn = syn_flagged.copy() if syn_flagged is not None else pd.DataFrame(columns=['sample', 'gene', 'reason_exclusion'])
    npa = npa_flagged.copy() if npa_flagged is not None else pd.DataFrame(columns=['sample', 'gene', 'reason_exclusion'])

    # Combine flagged tables: keep all reason_exclusion entries; if duplicates keep first
    all_flagged = pd.concat((syn, npa), ignore_index=True, sort=False)
    if all_flagged.empty:
        print('No flagged entries to create heatmap')
        return

    # choose top genes and samples by total counts
    gene_order = all_flagged['gene'].value_counts().head(top_n_genes).index.tolist()
    sample_order = all_flagged['sample'].value_counts().head(top_n_samples).index.tolist()

    # build pivot with one reason per cell; if multiple reasons, join with ';'
    pivot = all_flagged.groupby(['gene', 'sample'])['reason_exclusion'].apply(lambda s: ';'.join(sorted(set([x for x in s if x])))).unstack(fill_value='')
    pivot = pivot.reindex(index=gene_order, columns=sample_order, fill_value='')

    # determine reason to color: pick first of ';' separated if multiple, else ''
    # ensure we only split strings; non-string or empty -> ''
    cell_reason = pivot.fillna('').astype(str)
    # apply per-column to avoid potential type-checker issues with applymap
    cell_reason = cell_reason.apply(lambda col: col.map(lambda v: v.split(';')[0] if v else ''))

    # collect reasons and colors (reuse predefined)
    predefined = {
        'mutdensity = 0': '#8c8c8c',
        'high_mutdensity - zscore > 2': '#d62728',
        'low_mutdensity - zscore < -2': '#1f77b4',
    }
    unique_reasons = sorted(set([r for r in cell_reason.values.flatten() if r]))
    other_reasons = [r for r in unique_reasons if r not in predefined]
    # convert palette colors to hex strings so we consistently use str color codes
    from matplotlib.colors import to_hex
    palette = [to_hex(c) for c in sns.color_palette('tab20', n_colors=max(1, len(other_reasons)))]
    reason_colors = {r: predefined[r] for r in predefined if r in unique_reasons}
    for i, r in enumerate(other_reasons):
        reason_colors[r] = palette[i]

    # Map reasons to integers for heatmap; 0 = blank
    reason_to_int = {r: idx+1 for idx, r in enumerate(sorted(reason_colors.keys()))}
    int_matrix = cell_reason.replace('', pd.NA).apply(lambda col: col.map(lambda v: reason_to_int.get(v, pd.NA)))

    # Build colormap with white at 0 then listed colors
    from matplotlib.colors import ListedColormap
    cmap_list = ['white'] + [reason_colors[r] for r in sorted(reason_colors.keys())]
    cmap = ListedColormap(cmap_list)

    plt.figure(figsize=(max(8, 0.3 * len(sample_order)), max(6, 0.15 * len(gene_order))))
    ax = sns.heatmap(int_matrix.fillna(0).astype(int), cmap=cmap, cbar=False, linewidths=0.2)
    # set ticks labels
    ax.set_yticklabels(int_matrix.index, rotation=0)
    ax.set_xticklabels(int_matrix.columns, rotation=90)
    plt.title('Flagged genes (rows) x samples (columns) - colored by failure reason')
    plt.tight_layout()
    out_png = f"{output_prefix}.png"
    plt.savefig(out_png, dpi=200)
    plt.close()

    # also write legend mapping
    legend_df = pd.DataFrame([{'reason': r, 'color': reason_colors[r]} for r in sorted(reason_colors.keys())])
    legend_df.to_csv(f"{output_prefix}_legend.tsv", sep='\t', index=False)
    print(f"Wrote heatmap {out_png} and legend {output_prefix}_legend.tsv")


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

    # keep debug outputs for inspection
    syn_flagged.to_csv("debug.syn_flagged.tsv", sep="\t", index=False)
    npa_flagged.to_csv("debug.npa_flagged.tsv", sep="\t", index=False)

    # Combine flagged tables for annotation and plotting
    flagged_combined = pd.concat((syn_flagged, npa_flagged), ignore_index=True, sort=False)

    if flagged_combined.empty:
        print('No flagged entries found; skipping plots and annotating with no flags.')
        
        # Ensure an empty DataFrame with expected columns is passed
        flagged_for_annot = pd.DataFrame(columns=['sample', 'gene', 'reason_exclusion'])
    else:
        # Gene-based summary
        try:
            plot_flagged_summary(syn_flagged, npa_flagged)
        except Exception as e:
            print(f"Warning: plot_flagged_summary failed: {e}")
        # Sample-based summary
        try:
            plot_sample_flagged_summary(syn_flagged, npa_flagged)
        except Exception as e:
            print(f"Warning: plot_sample_flagged_summary failed: {e}")
        # Heatmap
        try:
            plot_flagged_heatmap(syn_flagged, npa_flagged)
        except Exception as e:
            print(f"Warning: plot_flagged_heatmap failed: {e}")

        flagged_for_annot = flagged_combined.drop_duplicates(subset=['sample', 'gene'], keep='first')

    annotated = annotate(omegas, flagged_for_annot)

    annotated.to_csv(output, sep="\t", index=False)
    print(f"Wrote annotated omegas to {output}")


if __name__ == "__main__":
    raise SystemExit(main())
