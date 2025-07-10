#!/usr/bin/env python


import polars as pl
import click

def create_consensus_panel(compact_annot_panel_path, depths_path, version, consensus_min_depth, compliance_threshold, genes=None):
    gene_list = [g.strip() for g in genes.split(",")] if genes else None

    # Load captured panel and depths
    compact_annot_panel_df = pl.read_csv(compact_annot_panel_path, separator="\t")
    depths_df = pl.read_csv(depths_path, separator="\t")

    # Ensure numeric columns are int
    for col in depths_df.columns[1:]:
        depths_df = depths_df.with_columns(pl.col(col).cast(pl.Int64))

    # Compliance: True if depth >= consensus_min_depth
    compliance_df = depths_df.select(depths_df.columns[2:]) >= consensus_min_depth

    # Calculate row compliance (fraction of samples meeting min depth)
    row_compliance = compliance_df.sum(axis=1) / compliance_df.width

    # Filter rows passing compliance threshold
    passing_rows = row_compliance >= compliance_threshold
    min_depths_df = depths_df.filter(passing_rows)

    # Filter captured panel to only keep minimally covered positions
    consensus_panel = compact_annot_panel_df.join(
        min_depths_df.select(["CHROM", "POS"]),
        on=["CHROM", "POS"],
        how="inner"
    )
    if gene_list and version != 'all':
        consensus_panel = consensus_panel.filter(pl.col("SYMBOL").is_in(gene_list))
    consensus_panel.write_csv(f"consensus.{version}.tsv", separator="\t")

    # Filter failing columns only for rows that pass the compliance threshold
    compliance_df_passing = compliance_df.filter(passing_rows)
    failing_mask = ~compliance_df_passing
    failing_columns_counts = []
    for row_idx, row in enumerate(failing_mask.rows()):
        for col_idx, failed in enumerate(row):
            if failed:
                failing_columns_counts.append({
                    "Row": row_idx,
                    "SAMPLE_ID": compliance_df_passing.columns[col_idx],
                    "Failed": True
                })

    if failing_columns_counts:
        failing_columns_counts_df = pl.DataFrame(failing_columns_counts)
        failure_counts_filtered = (
            failing_columns_counts_df.groupby("SAMPLE_ID")
            .count()
            .rename({"count": "FAILING_COUNT"})
        )
        failure_counts_filtered.write_csv(f"failing_consensus.{version}.tsv", separator="\t")


@click.command()
@click.option('--compact_annot_panel_path', type=click.Path(exists=True), required=True, help='Path to the compact annotation panel file.')
@click.option('--depths_path', type=click.Path(exists=True), required=True, help='Path to the depths file.')
@click.option('--version', type=str, required=True, help='Panel version.')
@click.option('--consensus_min_depth', type=int, required=True, help='Minimum depth for consensus.')
@click.option('--compliance_threshold', type=float, default=0.8, show_default=True, help='Compliance threshold (fraction of samples required to meet min depth).')
@click.option('--genes', default=None, help='Comma-separated list of genes to filter (e.g., TP53,KRAS,BRAF).')
def main(compact_annot_panel_path, depths_path, version, consensus_min_depth, compliance_threshold, genes):
    """
    CLI entry point for creating a consensus panel.
    """
    create_consensus_panel(compact_annot_panel_path, depths_path, version, consensus_min_depth, compliance_threshold, genes)

if __name__ == '__main__':
    main()