#!/usr/bin/env python

# Plotting gene saturation metrics


import sys
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import click
import os


pd.set_option('display.max_columns', None)



metrics_colors_dictionary = {"ofml"        : "viridis_r",
                                "ofml_score"  : "#6A33E0",
                                "omega_trunc" : "#FA5E32",
                                "omega_synon" : "#89E4A2",
                                "omega_miss"  : "#FABE4A",
                                "o3d_score"   : "#6DBDCC",
                                "o3d_cluster" : "skyblue",
                                "o3d_prob"    : "darkgray",
                                "frameshift"  : "#E4ACF4",
                                "inframe"     : "C5",
                                "hv_lines"    : "lightgray", # horizontal and vertical lines,
                                "hv_lines_needle" : "gray",
                                "needle_obs"  : "#003366",
                                "omega_miss_tert" : "#f5840c",
                                "omega_synon_tert": "#378c12",
                                "nonsense" : "#FA5E32",
                                "synonymous" : "#89E4A2",
                                "missense"  : "#FABE4A",
                                "indel"       : "#ECC4F7",
                                "splicing"    : "#A1C5DF",
                                }


# Define consistent order and base colors
consequence_order = ['nonsense', 'missense']
base_colors = {  'nonsense': metrics_colors_dictionary["nonsense"],
                 'missense': metrics_colors_dictionary["missense"]}


# Define color gradients for each consequence type, mapped to frequency bins
from matplotlib.colors import to_rgb, to_hex

def make_gradient_dict(color):
    rgb = np.array(to_rgb(color))
    factors = [1.4, 1.0, 0.6]  # lightest to darkest
    freq_bins = ['1', '2', '3+']
    return {freq: to_hex(np.clip(rgb * f, 0, 1)) for freq, f in zip(freq_bins, factors)}

color_gradients = {k: make_gradient_dict(v) for k, v in base_colors.items()}

 

def create_regions_df(rich_panel_file, expanded_panel_file, consensus_panel_file,
                      genes_list = []):
    consensus_df = pd.read_table(consensus_panel_file)
    if genes_list:
        consensus_df = consensus_df[consensus_df["GENE"].isin(genes_list)]
    rich_panel_df = pd.read_table(rich_panel_file)
    consensus_enriched = consensus_df.merge(rich_panel_df, how = 'left')
    expanded_panel_df = pd.read_table(expanded_panel_file)
    
    return consensus_enriched.drop(['GENE'], axis = 1).merge(expanded_panel_df, how = 'inner')




# Grouping logic as a function
def group_mutations(df, mode="aminoacid", count_mutations=True):
    """
    Group the input DataFrame by the specified mode.
    mode: 'aminoacid', 'protein_position', 'nucleotide_change', 'nucleotide_position'
    count_mutations: if True, count mutations (for mutation tables); if False, just group (for panel tables)

    define these as different options within a function:
    this should be coupled with a proper processing of the consensus_enriched_expanded table
    so that it is also grouped by in the same way, there is no need for counting in there

    """
    if mode == "aminoacid":
        group_cols = ["GENE", "IMPACT", "Protein_position", "Amino_acids"]
    elif mode == "protein_position":
        group_cols = ["GENE", "IMPACT", "Protein_position"]
    elif mode == "nucleotide_change":
        group_cols = ["GENE", "IMPACT", "CHROM", "POS", "REF", "ALT"]
    elif mode == "nucleotide_position":
        group_cols = ["GENE", "IMPACT", "CHROM", "POS"]
    else:
        raise ValueError(f"Unknown grouping mode: {mode}")

    if count_mutations:
        grouped = df.groupby(by=group_cols)["SAMPLE_ID"].size().to_frame("Count").reset_index()
    else:
        grouped = df[group_cols].drop_duplicates().reset_index(drop=True)

    return grouped


def compute_proportion_per_consequence_type(mutations_info,
                                            impacts = [['missense'],['nonsense']],
                                            regions = None
                                            ):

    if regions is not None:
        regions_list = regions
    else:
        regions_list = sorted(mutations_info["GENE"].unique())

    regions_results = []
    for region in regions_list:
        region_terms = region.split("--")
        gene_name = region_terms[0]
        segment_name = region_terms[1] if len(region_terms) > 1 else gene_name
        if '_ENSE0' in segment_name:
            region_type = 'exon'
            segment_name = int(segment_name.split("_")[1])
        elif segment_name != gene_name:
            region_type = 'domain'
            segment_name = segment_name
        else:
            region_type = 'gene'

        mutations_info_gene = mutations_info[mutations_info["GENE"] == region]
        for consequence in impacts:
            gene_cnsq_sites = mutations_info_gene[mutations_info_gene["IMPACT"].isin(consequence) ]
            total_sites = gene_cnsq_sites.shape[0]
            if total_sites > 0:
                mutated_sites = gene_cnsq_sites[gene_cnsq_sites["Count"] > 0].shape[0]
                regions_results.append((gene_name, segment_name, region_type, ",".join(consequence), total_sites, mutated_sites, mutated_sites/total_sites))

    results_df = pd.DataFrame(regions_results)
    results_df.columns = ["GENE_NAME", "SEGMENT_NAME", "SEGMENT_TYPE", "IMPACT", "TOTAL", "MUTATED", "PROPORTION_MUTATED"]
    return results_df



# Compute proportion per consequence type, stratified by mutation frequency (1, 2, 3+)
def compute_proportion_per_consequence_type_by_frequency(mutations_info,
                                                        impacts=[['missense'], ['nonsense']],
                                                        regions=None):
    """
    Compute the proportion of mutated positions per consequence type, stratified by mutation frequency (1, 2, 3+).
    Returns a DataFrame with columns:
    [GENE_NAME, SEGMENT_NAME, SEGMENT_TYPE, IMPACT, FREQ_BIN, TOTAL, MUTATED, PROPORTION_MUTATED]
    """
    if regions is not None:
        regions_list = regions
    else:
        regions_list = sorted(mutations_info["GENE"].unique())

    freq_bins = [(1, 1, '1'), (2, 2, '2'), (3, float('inf'), '3+')]
    results = []
    for region in regions_list:
        region_terms = region.split("--")
        gene_name = region_terms[0]
        segment_name = region_terms[1] if len(region_terms) > 1 else gene_name
        if '_ENSE0' in segment_name:
            region_type = 'exon'
            segment_name = int(segment_name.split("_")[1])
        elif segment_name != gene_name:
            region_type = 'domain'
            segment_name = segment_name
        else:
            region_type = 'gene'

        mutations_info_gene = mutations_info[mutations_info["GENE"] == region]
        for consequence in impacts:
            gene_cnsq_sites = mutations_info_gene[mutations_info_gene["IMPACT"].isin(consequence)]
            total_sites = gene_cnsq_sites.shape[0]
            for min_c, max_c, freq_label in freq_bins:
                if total_sites > 0:
                    mutated_sites_in_freq = gene_cnsq_sites[(gene_cnsq_sites["Count"] >= min_c) & (gene_cnsq_sites["Count"] <= max_c)].shape[0]
                    mutated_sites = gene_cnsq_sites[(gene_cnsq_sites["Count"] >= min_c)].shape[0]
                    results.append((gene_name, segment_name, region_type, ",".join(consequence), freq_label, total_sites, mutated_sites, mutated_sites/total_sites, mutated_sites_in_freq, mutated_sites_in_freq/total_sites if total_sites > 0 else 0))
    results_df = pd.DataFrame(results)
    results_df.columns = ["GENE_NAME", "SEGMENT_NAME", "SEGMENT_TYPE", "IMPACT", "FREQ_BIN", "TOTAL", "MUTATED", "PROPORTION_MUTATED", "MUTATED_IN_FREQ", "PROPORTION_MUTATED_IN_FREQ"]
    return results_df




def plot_genes(df, mode=None):
    df_genes = df[df["SEGMENT_TYPE"] == "gene"]

    g = sns.catplot(
        data=df_genes,
        x="GENE_NAME",
        y="PROPORTION_MUTATED",
        row="IMPACT",
        hue  ='IMPACT',
        kind="bar",
        order=sorted(df_genes["GENE_NAME"].unique()),
        palette=base_colors,
        height=1.5,
        aspect=4.5,
        sharey=True
    )

    g.set_axis_labels("Gene", "")
    g.set_titles("")
    g.set(ylim=(0, 1))

    g._legend.remove()
    g.fig.text(0, 0.5, "Proportion mutated", va='center', rotation='vertical')

    for ax in g.axes.flatten():
        ax.tick_params(axis='x', rotation=45)

    # Save plot to file
    suffix = f"_{mode}" if mode else ""
    plot_path = f"{plots_dir}/saturation_genes{suffix}.pdf"
    g.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close(g.fig)

def plot_domains(df, genes = None, mode=None):
    if genes is None:
        genes = df["GENE_NAME"].unique()

    seg_type = 'domain'
    suffix = f"_{mode}" if mode else ""
    pdf_path = f"{plots_dir}/saturation_domains_all{suffix}.pdf"
    with PdfPages(pdf_path) as pdf:
        for gene in genes:
            df_gene = df[(df["GENE_NAME"] == gene) & (df["SEGMENT_TYPE"] == "domain")]
            if df_gene.empty:
                continue
            g = sns.catplot(
                data=df_gene,
                x="SEGMENT_NAME",
                y="PROPORTION_MUTATED",
                row="IMPACT",
                hue = 'IMPACT',
                kind="bar",
                order=sorted(df_gene["SEGMENT_NAME"].unique()),
                palette=base_colors,
                height=1.5,
                aspect=4.5,
                sharey=True
            )
            g.fig.suptitle(f"{gene} - {seg_type.capitalize()}", y=1)
            g.set_axis_labels("Domains", "")
            g.set_titles("")
            g.set(ylim=(0, 1))
            g._legend.remove()
            g.fig.text(0, 0.5, "Proportion mutated", va='center', rotation='vertical')
            for ax in g.axes.flatten():
                ax.tick_params(axis='x', rotation=45)
            pdf.savefig(g.fig, bbox_inches='tight', dpi=300)
            plt.close(g.fig)

def plot_exons(df, genes = None, mode=None):
    if genes is None:
        genes = df["GENE_NAME"].unique()

    seg_type = 'exon'
    suffix = f"_{mode}" if mode else ""
    pdf_path = f"{plots_dir}/saturation_exons_all{suffix}.pdf"
    with PdfPages(pdf_path) as pdf:
        for gene in genes:
            df_gene = df[(df["GENE_NAME"] == gene) & (df["SEGMENT_TYPE"] == "exon")]
            if df_gene.empty:
                continue
            g = sns.catplot(
                data=df_gene,
                x="SEGMENT_NAME",
                y="PROPORTION_MUTATED",
                row="IMPACT",
                hue = 'IMPACT',
                kind="bar",
                order=sorted(df_gene["SEGMENT_NAME"].unique()),
                palette=base_colors,
                height=1.5,
                aspect=4.5,
                sharey=True
            )
            g.fig.suptitle(f"{gene} - {seg_type.capitalize()}", y=1)
            g.set_axis_labels("Exons", "")
            g.set_titles("")
            g.set(ylim=(0, 1))
            g._legend.remove()
            g.fig.text(0.05, 0.5, "Proportion mutated", va='center', rotation='vertical')
            for ax in g.axes.flatten():
                ax.tick_params(axis='x', rotation=0)
            pdf.savefig(g.fig, bbox_inches='tight', dpi=300)
            plt.close(g.fig)


# Frequency-stratified domain plot
def plot_domains_by_freq(df, genes = None, mode=None):
    if genes is None:
        genes = df["GENE_NAME"].unique()

    freq_bin_order = ['3+', '2', '1']
    seg_type = 'domain'
    suffix = f"_{mode}" if mode else ""
    pdf_path = f"{plots_dir}/saturation_domains_byfreq_all{suffix}.pdf"
    with PdfPages(pdf_path) as pdf:
        for gene in genes:
            df_gene = df[(df["GENE_NAME"] == gene) & (df["SEGMENT_TYPE"] == "domain")]
            if df_gene.empty:
                continue
            segs = sorted(df_gene["SEGMENT_NAME"].unique())
            impacts = consequence_order
            x = np.arange(len(segs))
            fig, axes = plt.subplots(len(impacts), 1, figsize=(max(9, len(segs)*0.5), 4*len(impacts)), sharex=True)
            if len(impacts) == 1:
                axes = [axes]
            legend_patches = []
            for j, impact in enumerate(impacts):
                ax = axes[j]
                bottoms = np.zeros(len(segs))
                for freq_label in freq_bin_order:
                    df_impact = df_gene[(df_gene["IMPACT"] == impact) & (df_gene["FREQ_BIN"] == freq_label)]
                    heights = []
                    for seg in segs:
                        row = df_impact[df_impact["SEGMENT_NAME"] == seg]
                        heights.append(row["PROPORTION_MUTATED_IN_FREQ"].iloc[0] if not row.empty else 0)
                    bar = ax.bar(x, heights, width=0.8,
                                 color=color_gradients[impact][freq_label],
                                 label=f"{freq_label}x" if j == 0 else None,
                                 bottom=bottoms)
                    if j == 0:
                        legend_patches.append(bar)
                    bottoms += np.array(heights)
                ax.set_ylabel(f"Proportion mutated\n({impact})")
                ax.set_ylim(0, 1)
            # Only add legend to the top axis, with one entry per freq bin
            axes[0].legend([p[0] for p in legend_patches], [f"{f}" for f in freq_bin_order], title="Frequency", frameon=False)
            axes[-1].set_xticks(x)
            axes[-1].set_xticklabels(segs, rotation=45)
            axes[-1].set_xlabel("Domains")
            fig.suptitle(f"{gene} - Domains (by mutation frequency)")
            plt.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig, bbox_inches='tight', dpi=300)
            plt.close(fig)

# Frequency-stratified exon plot
def plot_exons_by_freq(df, genes = None, mode=None):
    if genes is None:
        genes = df["GENE_NAME"].unique()
    freq_bin_order = ['3+', '2', '1']
    seg_type = 'exon'
    suffix = f"_{mode}" if mode else ""
    pdf_path = f"{plots_dir}/saturation_exons_byfreq_all{suffix}.pdf"
    with PdfPages(pdf_path) as pdf:
        for gene in genes:
            df_gene = df[(df["GENE_NAME"] == gene) & (df["SEGMENT_TYPE"] == "exon")]
            if df_gene.empty:
                continue
            segs = sorted(df_gene["SEGMENT_NAME"].unique())
            impacts = consequence_order
            x = np.arange(len(segs))
            fig, axes = plt.subplots(len(impacts), 1, figsize=(max(9, len(segs)*0.5), 4*len(impacts)), sharex=True)
            if len(impacts) == 1:
                axes = [axes]

            legend_patches = []
            for j, impact in enumerate(impacts):
                ax = axes[j]
                bottoms = np.zeros(len(segs))
                for freq_label in freq_bin_order:
                    df_impact = df_gene[(df_gene["IMPACT"] == impact) & (df_gene["FREQ_BIN"] == freq_label)]
                    heights = []
                    for seg in segs:
                        row = df_impact[df_impact["SEGMENT_NAME"] == seg]
                        heights.append(row["PROPORTION_MUTATED_IN_FREQ"].iloc[0] if not row.empty else 0)
                    bar = ax.bar(x, heights, width=0.8,
                                 color=color_gradients[impact][freq_label],
                                 label=f"{freq_label}x" if j == 0 else None,
                                 bottom=bottoms)
                    if j == 0:
                        legend_patches.append(bar)
                    bottoms += np.array(heights)
                ax.set_ylabel(f"Proportion mutated\n({impact})")
                ax.set_ylim(0, 1)
            axes[0].legend([p[0] for p in legend_patches], [f"{f}" for f in freq_bin_order], title="Frequency", frameon=False)
            axes[-1].set_xticks(x)
            axes[-1].set_xticklabels(segs, rotation=0)
            axes[-1].set_xlabel("Exons")
            fig.suptitle(f"{gene} - Exons (by mutation frequency)")
            plt.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig, bbox_inches='tight', dpi=300)
            plt.close(fig)



# Call all three plotting functions from the same table
def plot_all_saturation_tables(df, mode=None):
    plot_genes(df, mode)

    genes_list = df["GENE_NAME"].unique()
    genes_list = genes_list[:200]
    plot_domains(df, genes=genes_list, mode=mode)
    plot_exons(df, genes=genes_list, mode=mode)


# Call all three frequency-stratified plotting functions from the same table
def plot_all_saturation_tables_by_freq(df, mode=None):
    genes_list = df["GENE_NAME"].unique()
    genes_list = genes_list[:200]
    
    plot_domains_by_freq(df, genes=genes_list, mode=mode)
    plot_exons_by_freq(df, genes=genes_list, mode=mode)


def generate_all_saturation_plots(consensus_enriched_expanded, somatic_maf_clean,
                                    grouping_modes=["aminoacid", "protein_position", "nucleotide_change", "nucleotide_position"]):
    """
    Generate plots for each grouping mode. Accepts pre-loaded consensus regions table and cleaned somatic maf.
    """
    for mode in grouping_modes:
        # merge mutations with panel for this mode
        consensus_with_mutations = consensus_enriched_expanded.merge(somatic_maf_clean[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID']], how='inner')

        # Group mutations
        consensus_with_mutation_counts = group_mutations(consensus_with_mutations, mode=mode, count_mutations=True)
        consensus_enriched_grouped = group_mutations(consensus_enriched_expanded, mode=mode, count_mutations=False)
        consensus_with_mutations_info = consensus_enriched_grouped.merge(consensus_with_mutation_counts, how='left')
        consensus_with_mutations_info["Count"] = consensus_with_mutations_info["Count"].fillna(0).astype(int)

        # Compute proportions
        proportion_results = compute_proportion_per_consequence_type(consensus_with_mutations_info)
        proportion_results_by_freq = compute_proportion_per_consequence_type_by_frequency(consensus_with_mutations_info)

        # Plot original (non-frequency) plots
        plot_all_saturation_tables(proportion_results, mode)

        # Plot frequency-stratified plots (for exons and domains only)
        plot_all_saturation_tables_by_freq(proportion_results_by_freq, mode)


@click.command()
@click.option('--rich-panel', default='captured_panel.tab.compact_rich.tsv', help='Rich panel TSV')
@click.option('--expanded-panel', default='exons_consensus_panel_with_hotspots.tsv', help='Expanded panel TSV')
@click.option('--consensus-panel', default='consensus.exons_splice_sites.tsv', help='Consensus panel TSV')
@click.option('--maf', default='all_samples.somatic.mutations.tsv', help='Somatic MAF TSV')
@click.option('--plots-dir', default='plot', help='Output directory for plots')
@click.option('--genes', default='', help='Comma-separated list of genes to include (optional)')
@click.option('--grouping-modes', default='aminoacid,protein_position,nucleotide_change,nucleotide_position', help='Comma-separated grouping modes to run')
def cli(rich_panel, expanded_panel, consensus_panel, maf, plots_dir, genes, grouping_modes):
    """CLI entry point to generate saturation plots."""
    os.makedirs(plots_dir, exist_ok=True)
    # update module-level plots_dir so plotting functions write to the requested directory
    globals()['plots_dir'] = plots_dir

    genes_list = [g.strip() for g in genes.split(',')] if genes else []

    consensus_enriched_expanded = create_regions_df(rich_panel, expanded_panel, consensus_panel, genes_list=genes_list)

    somatic_maf = pd.read_table(maf)
    somatic_maf_clean = somatic_maf[(somatic_maf["TYPE"] == 'SNV')
                                    & (~somatic_maf["FILTER.not_in_exons"]) 
                                    & (somatic_maf['canonical_Protein_position'] != '-')].reset_index(drop=True)

    grouping_modes_list = [m.strip() for m in grouping_modes.split(',')]

    try :
        # Run generation
        generate_all_saturation_plots(consensus_enriched_expanded, somatic_maf_clean, grouping_modes=grouping_modes_list)
    except Exception as e:
        print(f"Error during plot generation: {e}", file=sys.stderr)


if __name__ == '__main__':
    cli()
