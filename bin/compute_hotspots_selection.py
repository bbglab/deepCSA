#!/usr/bin/env python

import click
import numpy as np
import pandas as pd
import scipy.stats as stats
import glob
import os

from utils import MIN_NONZERO_PVALUE
from omega_comparison_per_site import poisson_pvalue, benjamini_hochberg

def load_panel_hotspots(panel_file, hotspots_file):
    panel = pd.read_csv(panel_file, sep="\t", compression="gzip")
    hotspots = pd.read_csv(hotspots_file, sep="\t")
    
    # Merge hotspots with panel to find which sites are hotspots
    # Hotspots have CHROM, POS, MUTTYPE. MUTTYPE is REF>ALT.
    hotspots['REF'] = hotspots['MUTTYPE'].str.split('>').str[0]
    hotspots['ALT'] = hotspots['MUTTYPE'].str.split('>').str[1]
    
    # Some hotspots might have MUTTYPE = '-' for non-SNV, but let's assume standard format for now
    hotspots_panel = panel.merge(hotspots, on=['CHROM', 'POS', 'REF', 'ALT'], how='inner')
    
    return hotspots_panel

def process_comparison(comparison_file, hotspots_panel, size_type, output_prefix):
    comp_df = pd.read_csv(comparison_file, sep="\t", compression="gzip")
    
    if size_type == "site":
        group_cols = ['CHROM', 'POS', 'REF', 'ALT', 'GENE']
    elif size_type == "aminoacid":
        group_cols = ['GENE', 'Feature', 'Protein_position']
    elif size_type == "aminoacid_change":
        group_cols = ['GENE', 'Feature', 'Protein_position', 'Amino_acids']
    else:
        raise ValueError(f"Unknown size type: {size_type}")

    # Determine which entries are hotspots
    # Get the unique hotspots for the current grouping
    hotspots_grouped = hotspots_panel[group_cols].drop_duplicates()
    hotspots_grouped['Hotspot'] = 'Yes'

    # Merge with the comparison data
    comp_df = comp_df.merge(hotspots_grouped, on=group_cols, how='left')
    comp_df['Hotspot'] = comp_df['Hotspot'].fillna('No')

    # Group by GENE and Hotspot
    grouped_size = comp_df.groupby(['GENE', 'Hotspot']).size().reset_index(name='Count')
    grouped_size_mut = comp_df.groupby(['GENE', 'Hotspot']).agg({'OBSERVED_MUTS': lambda x: (x != 0).sum()}).reset_index(name = 'CountDiffMutatedSites')
    grouped_size = grouped_size.merge(grouped_size_mut, on=['GENE', 'Hotspot'])
    grouped = comp_df.groupby(['GENE', 'Hotspot'])[['OBSERVED_MUTS', 'EXPECTED_MUTS']].sum().reset_index()

    grouped["OBS/EXP"] = (grouped["OBSERVED_MUTS"] / grouped["EXPECTED_MUTS"]).fillna(0)
    grouped["OBS/EXP"] = grouped["OBS/EXP"].replace([np.inf, -np.inf], 0)
    grouped["p_value"] = grouped.apply(lambda row: poisson_pvalue(row["OBSERVED_MUTS"], row["EXPECTED_MUTS"]), axis=1)
    
    grouped["p_value"] = grouped["p_value"].replace(0, MIN_NONZERO_PVALUE)
    grouped = grouped.merge(grouped_size, on=['GENE', 'Hotspot'])
    grouped = grouped[["GENE", "Hotspot", "Count", "CountDiffMutatedSites", "OBSERVED_MUTS", "EXPECTED_MUTS", "OBS/EXP", "p_value"]]

    # TO BE FIXED: Adjusted p-values are not being calculated correctly. The following code is commented out for now.
    # grouped["p_value_adj"] = np.nan
    # for gene, gene_df in grouped.groupby("GENE", dropna=False):
    #     valid_mask = gene_df["p_value"].notna()
    #     if not valid_mask.any():
    #         continue
    #     ## something is wrong here with the adjusted p-values, they are not being assigned correctly. Let's fix that.
    #     adjusted = benjamini_hochberg(gene_df.loc[valid_mask, "p_value"].to_numpy())
    #     grouped.loc[gene_df.index[valid_mask], "p_value_adj"] = adjusted

    # Write output
    output_file = f"{output_prefix}.{size_type}.hotspots_selection.tsv.gz"
    grouped.to_csv(output_file, sep="\t", index=False)
    click.echo(f"Results for size '{size_type}' written to {output_file}")


@click.command()
@click.option('--comparisons', multiple=True, type=click.Path(exists=True), required=True, help="Path to comparison files.")
@click.option('--panel-file', type=click.Path(exists=True), required=True, help="Path to captured panel file (gzip compressed).")
@click.option('--hotspots-file', type=click.Path(exists=True), required=True, help="Path to hotspots file.")
@click.option('--output-prefix', type=str, required=True, help="Output file prefix.")
def main(comparisons, panel_file, hotspots_file, output_prefix):
    """Compute selection for known hotspots based on site comparison output."""
    hotspots_panel = load_panel_hotspots(panel_file, hotspots_file)

    for comp_file in comparisons:
        filename = os.path.basename(comp_file)
        if ".site.comparison." in filename:
            size_type = "site"
        elif ".aminoacid_change.comparison." in filename:
            size_type = "aminoacid_change"
        elif ".aminoacid.comparison." in filename:
            size_type = "aminoacid"
        else:
            click.echo(f"Warning: Could not determine size type from filename {filename}. Skipping.")
            continue
        
        click.echo(f"Processing size: {size_type} from {comp_file}")
        process_comparison(comp_file, hotspots_panel, size_type, output_prefix)

if __name__ == "__main__":
    main()
