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



def load_flagged_tables(paths: List[Path]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
        return pd.DataFrame(columns=["sample", "gene", "reason_exclusion"]), pd.DataFrame(columns=["sample", "gene", "reason_exclusion"]), pd.DataFrame(columns=["sample", "gene", "reason_exclusion"]), pd.DataFrame(columns=["sample", "gene", "reason_exclusion"])

    all_df = pd.concat(pieces, ignore_index=True, sort=False)

    # These are the columns that all_df will contain:
    # ID      reason_exclusion        cohort  regions criteria

    all_df = all_df.rename(columns={
        "cohort": "sample",
        "ID": "gene",
    })
    syn_sample_df = all_df[(all_df["regions"] == "synonymous")
                            & (all_df["criteria"] == "per_sample")
                            ][["sample", "gene", "reason_exclusion"]].reset_index(drop=True)
    syn_sample_df.rename(columns={"sample": "cohort", "gene": "sample"}, inplace=True)
    syn_gene_df = all_df[(all_df["regions"] == "synonymous")
                            & (all_df["criteria"] == "per_gene")
                            ][["sample", "gene", "reason_exclusion"]].reset_index(drop=True)
    syn_gene_df.rename(columns={"sample": "cohort"}, inplace=True)

    npa_sample_df = all_df[(all_df["regions"] == "non_protein_affecting")
                            & (all_df["criteria"] == "per_sample")
                            ][["sample", "gene", "reason_exclusion"]].reset_index(drop=True)
    npa_sample_df.rename(columns={"sample": "cohort", "gene": "sample"}, inplace=True)

    npa_gene_df = all_df[(all_df["regions"] == "non_protein_affecting")
                            & (all_df["criteria"] == "per_gene")
                            ][["sample", "gene", "reason_exclusion"]].reset_index(drop=True)
    npa_gene_df.rename(columns={"sample": "cohort"}, inplace=True)

    return syn_sample_df, syn_gene_df, npa_sample_df, npa_gene_df


def plot_flagged_summary(syn_flagged: pd.DataFrame, npa_flagged: pd.DataFrame,
                         output_prefix: str = 'flagged_gene', top_n: int = 50,
                         variable : str = 'gene') -> None:
    """Plot stacked bar summaries of flagged genes for synonymous and non-protein-affecting.

    Uses the gene order from `syn_flagged` (most frequently flagged first) to
    order the bars. Bars are stacked by `reason_exclusion` and colors for a
    small set of known reasons are defined a priori; any additional reasons
    will get colors from the seaborn color palette.
    """
    # Ensure dataframes exist
    syn = syn_flagged.copy()
    npa = npa_flagged.copy()

    # Determine gene order from synonymously flagged genes (descending frequency)
    syn_order = syn[variable].value_counts().index.tolist()

    # Include any genes present in npa but not in syn at the end
    npa_genes = [g for g in npa[variable].unique() if g not in syn_order]
    order = syn_order + npa_genes

    # Count occurrences per gene x reason
    syn_counts = syn.groupby([variable, 'reason_exclusion']).size().unstack(fill_value=0)
    npa_counts = npa.groupby([variable, 'reason_exclusion']).size().unstack(fill_value=0)

    # Reindex to the chosen order; keep only top_n genes by total syn+npa counts
    total_counts = (syn.groupby(variable).size().add(npa.groupby(variable).size(), fill_value=0)).sort_values(ascending=False)
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
    syn_reasons = list(syn_counts.columns)
    npa_reasons = list(npa_counts.columns)

    # Build color list: use predefined colors then palette for others
    syn_colors = []
    for r in syn_reasons:
        syn_colors.append(predefined.get(r, 'grey'))

    npa_colors = []
    for r in npa_reasons:
        npa_colors.append(predefined.get(r, 'grey'))

    # Plotting stacked horizontal barplots (two subplots stacked vertically)
    fig, axes = plt.subplots(2, 1, figsize=(10, max(6, 0.25 * len(order) * 2)), sharex=True)

    # Synonymous stacked bar
    if not syn_counts.empty:
        syn_counts[syn_reasons].plot(kind='barh', stacked=True, color=syn_colors, ax=axes[0])
        axes[0].set_title('Synonymous flagged cases (by gene)')
        axes[0].set_xlabel('Count')
    else:
        axes[0].text(0.5, 0.5, 'No synonymous flagged cases', ha='center')

    # NPA stacked bar
    if not npa_counts.empty:
        npa_counts[npa_reasons].plot(kind='barh', stacked=True, color=npa_colors, ax=axes[1])
        axes[1].set_title('Non-protein-affecting flagged cases (by gene)')
        axes[1].set_xlabel('Count')
    else:
        axes[1].text(0.5, 0.5, 'No non-protein-affecting flagged cases', ha='center')

    plt.tight_layout()
    fig.savefig(f"{output_prefix}.{variable}.cases_summary.png", dpi=300)
    plt.close(fig)

    # # Also write TSV summary
    # summary = pd.DataFrame({variable: order})
    # for r in reasons:
    #     summary[r] = syn_counts.get(r, 0).reindex(order, fill_value=0).values + npa_counts.get(r, 0).reindex(order, fill_value=0).values
    # summary.to_csv(f"{output_prefix}_cases_counts.flagged.tsv", sep='\t', index=False)
    # print(f"Wrote flagged summary plot {out_png} and table {output_prefix}_cases_counts.flagged.tsv")


def annotate(omegas: pd.DataFrame, gene_flagged: pd.DataFrame, sample_flagged: pd.DataFrame ) -> pd.DataFrame:
    """Annotate omegas DataFrame with flagged info.

    The returned DataFrame includes two new columns: `flagged` (bool) and
    `flag_reason` (string).
    """
    # omegas has these columns:
    # gene	sample	impact	mutations	dnds	pvalue	lower	upper

    # and flagged has:
    # sample, gene, reason_exclusion

    omegas["original_gene"] = omegas["gene"].copy()
    omegas["gene"] = omegas["gene"].str.split("--").str[0]
    annotated_omegas_gene = omegas.merge(gene_flagged, how="left",
                                         left_on=["sample", "gene"],
                                         right_on=["cohort", "gene"],
                                         suffixes = ("", "_gene")
                                         ).drop(columns=["cohort"])
    annotated_omegas = annotated_omegas_gene.merge(sample_flagged.drop(columns=["cohort"]),
                                                   how="left",
                                                   on=["sample"],
                                                   suffixes = ("", "_sample")
                                                   )
    
    annotated_omegas["gene"] = annotated_omegas["original_gene"]
    annotated_omegas["flag_reason"] = annotated_omegas["reason_exclusion"].fillna("") + annotated_omegas["reason_exclusion_sample"].fillna("")
    annotated_omegas["flagged"] = annotated_omegas["flag_reason"] != ""    
    annotated_omegas = annotated_omegas.drop(columns=["original_gene", "reason_exclusion", "reason_exclusion_sample"])

    return annotated_omegas[['gene', 'sample', 'impact', 'mutations',
                             'dnds', 'pvalue', 'lower', 'upper',
                             'pvalue_adj',
                             'flagged', 'flag_reason']]


def plot_flagged_heatmap(syn_flagged: pd.DataFrame, npa_flagged: pd.DataFrame,
                            output_prefix: str = 'flagged_heatmap', top_n_genes: int = 200,
                            top_n_samples: int = 200, mode: str = 'combined',
                            variable: str = 'gene'
                            ) -> None:

    """Create a genes x samples heatmap showing failure reasons per cell.

    mode: 'combined' (default) aggregates syn+npa into one heatmap. mode: 'separate'
    will produce two heatmaps (suffixes '_syn' and '_npa').
    """
    def _heatmap_one(df: pd.DataFrame, prefix: str, variable: str) -> None:
        syn = df.copy()
        if syn.empty:
            print(f'No flagged entries to create heatmap for {prefix}')
            return

        gene_order = syn[variable].value_counts().head(top_n_genes).index.tolist()
        sample_order = syn['cohort'].value_counts().head(top_n_samples).index.tolist()

        pivot = syn.groupby([variable, 'cohort'])['reason_exclusion'].apply(lambda s: ';'.join(sorted(set([x for x in s if x])))).unstack(fill_value='')
        pivot = pivot.reindex(index=gene_order, columns=sample_order, fill_value='')
        print(pivot)

        cell_reason = pivot.fillna('').astype(str)
        cell_reason = cell_reason.apply(lambda col: col.map(lambda v: v.split(';')[0] if v else ''))

        predefined = {
            'mutdensity = 0': '#8c8c8c',
            'high_mutdensity - zscore > 2': '#d62728',
            'low_mutdensity - zscore < -2': '#1f77b4',
        }
        unique_reasons = sorted(set([r for r in cell_reason.values.flatten() if r]))
        reason_colors = {r: predefined[r] for r in predefined if r in unique_reasons}

        reason_to_int = {r: idx+1 for idx, r in enumerate(sorted(reason_colors.keys()))}
        int_matrix = cell_reason.replace('', pd.NA).apply(lambda col: col.map(lambda v: reason_to_int.get(v, pd.NA)))

        from matplotlib.colors import ListedColormap
        cmap_list = ['white'] + [reason_colors[r] for r in sorted(reason_colors.keys())]
        cmap = ListedColormap(cmap_list)

        plt.figure(figsize=(max(8, 0.3 * len(sample_order)), max(6, 0.15 * len(gene_order))))
        ax = sns.heatmap(int_matrix.fillna(0).astype(int), cmap=cmap, cbar=False, linewidths=0.2)
        ax.set_yticklabels(int_matrix.index, rotation=0)
        ax.set_xticklabels(int_matrix.columns, rotation=90)
        plt.title('Flagged genes (rows) x cohorts (columns) - colored by failure reason')
        plt.tight_layout()
        out_png = f"{prefix}.{variable}.png"
        plt.savefig(out_png, dpi=200)
        plt.close()

    # combined df
    all_df = pd.concat((syn_flagged.copy() if syn_flagged is not None else pd.DataFrame(columns=['cohort',variable,'reason_exclusion']),
                        npa_flagged.copy() if npa_flagged is not None else pd.DataFrame(columns=['cohort',variable,'reason_exclusion'])),
                       ignore_index=True, sort=False)

    if mode == 'combined':
        if all_df.empty:
            print('No flagged entries to create heatmap')
            return
        _heatmap_one(all_df, output_prefix, variable)
    elif mode == 'separate':
        _heatmap_one(syn_flagged.copy() if syn_flagged is not None else pd.DataFrame(columns=['cohort',variable,'reason_exclusion']), f"{output_prefix}_syn", variable)
        _heatmap_one(npa_flagged.copy() if npa_flagged is not None else pd.DataFrame(columns=['cohort',variable,'reason_exclusion']), f"{output_prefix}_npa", variable)
    else:
        raise ValueError(f"Unknown mode for plot_flagged_heatmap: {mode}")





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
    omegas = pd.read_csv(omegas_path, sep="\t", header=0)

    syn_flagged_sample, syn_flagged_gene, npa_flagged_sample, npa_flagged_gene = load_flagged_tables(flagged_paths)

    # keep debug outputs for inspection
    syn_flagged_sample.to_csv("debug.syn_flagged_sample.tsv", sep="\t", index=False)
    syn_flagged_gene.to_csv("debug.syn_flagged_gene.tsv", sep="\t", index=False)
    npa_flagged_sample.to_csv("debug.npa_flagged_sample.tsv", sep="\t", index=False)
    npa_flagged_gene.to_csv("debug.npa_flagged_gene.tsv", sep="\t", index=False)

    if syn_flagged_sample.empty and npa_flagged_sample.empty:
        print('No flagged entries found; skipping plots and annotating with no flags.')

    else:
        try:
            # Gene-based summary
            plot_flagged_summary(syn_flagged_gene, npa_flagged_gene, variable = 'gene')
        except Exception as e:
            print(f"Warning: plot_flagged_summary failed: {e}")

        try:
            # Sample-based summaries (separate syn/npa)
            plot_flagged_summary(syn_flagged_sample, npa_flagged_sample, variable = 'sample')
        except Exception as e:
            print(f"Warning: plot_flagged_summary failed: {e}")

        try:
            # Heatmaps (separate syn/npa)
            plot_flagged_heatmap(syn_flagged_gene,   npa_flagged_gene   , mode='separate', variable='gene'  )
            plot_flagged_heatmap(syn_flagged_sample, npa_flagged_sample , mode='separate', variable='sample')
        except Exception as e:
            print(f"Warning: plot_flagged_summary failed: {e}")

    annotated = annotate(omegas, syn_flagged_gene, syn_flagged_sample)

    annotated.to_csv(output, sep="\t", index=False)
    print(f"Wrote annotated omegas to {output}")


if __name__ == "__main__":
    raise SystemExit(main())
