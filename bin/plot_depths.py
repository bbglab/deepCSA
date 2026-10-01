#!/usr/bin/env python


import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Plotting size limits to keep the PDF tractable for large panels / cohorts.
# Any section that exceeds these limits is skipped (a message is printed to the log).
MAX_CATEGORIES_PER_PLOT = 30     # max samples/genes in a single boxplot/stripplot page
MAX_REGIONS_PER_PLOT = 200       # max regions in a single region boxplot/stripplot page
MAX_LINES_PER_PLOT = 100         # max lines (hue categories) in a lineplot
MAX_LINEPLOT_POINTS = 50000      # max total points (lines x categories) in a lineplot
MAX_HEATMAP_CELLS = 25000        # max cells (rows x cols) in a heatmap
MAX_CLUSTERMAP_DIM = 300         # max matrix dimension for a clustermap
MAX_FIG_WIDTH = 40               # cap on figure width (inches)
MAX_FIG_HEIGHT = 30              # cap on figure height (inches)


def _section_allowed(ok, section, reason):
    """Return True if a plot section should be drawn; otherwise log why it is skipped."""
    if not ok:
        print(f"[plot_depths] Skipping plot section '{section}': {reason}")
    return ok


def annotate_depth_region(row, dp_df):

    # look in the DEPTH file for the positions whithin this region
    dp_df_f = dp_df.loc[(dp_df.CHROM == row.CHROM) & (dp_df.POS >= row.START) & (dp_df.POS <= row.END)]

    # annotate total sequenced bps EXON_SEQ per SAMPLE_ID
    samples = dp_df_f.columns[~dp_df_f.columns.isin(["CHROM","POS"])]
    depths = []
    size = []
    for s in samples:
        depths.append(dp_df_f.loc[:,s].sum())
        size.append(dp_df_f.loc[:,s].count())

    row["SAMPLE_ID"] = list(samples)
    row["EXON_SEQ"] = depths
    row["EXON_SIZE"] = size


    return row


def load_depths_and_panel(sample_name, depth_file, panel_bed6_file):
    """
    Loads and prepares depth and panel data, returning all necessary objects for downstream analysis.
    """
    depth_df = pd.read_csv(f"{depth_file}", sep="\t", header=0)
    depth_df = depth_df.drop(["CONTEXT"], axis='columns')

    samples_list = depth_df.columns[~depth_df.columns.isin(["CHROM", "POS"])]

    stats_per_sample = pd.DataFrame(depth_df[samples_list].describe())
    stats_per_sample.to_csv(f"{sample_name}.depth_per_sample.stats.tsv", sep='\t', float_format="%.3f", header=True, index=True)

    avgdepth_per_sample_names = pd.DataFrame(depth_df[samples_list].mean().T).reset_index()
    avgdepth_per_sample_names.columns = ["SAMPLE_ID", "avg_depth_sample"]
    avgdepth_per_sample_names.to_csv(f"{sample_name}.avgdepth_per_sample.tsv", sep='\t', float_format="%.3f", header=True, index=False)

    # Sort the list of samples by depth, ascending
    samples_list = list(avgdepth_per_sample_names.sort_values(by=["avg_depth_sample"], ascending=False)["SAMPLE_ID"].values)

    # Load BED file with probes info
    bed6_probes_df = pd.read_csv(panel_bed6_file, sep="\t", header=0).reset_index()
    bed6_probes_df["EXON"] = bed6_probes_df["ELEMENT"].astype(str) + "_" + bed6_probes_df["index"].astype(str)
    bed6_probes_df = bed6_probes_df[["CHROMOSOME", "START", "END", "ELEMENT", "EXON"]]
    bed6_probes_df.columns = ["CHROM", "START", "END", "GENE", "EXON"]
    bed6_probes_df["CHROM"] = 'chr' + bed6_probes_df["CHROM"].astype(str)

    # Store bed6_probes_df
    bed6_probes_df.to_csv(f"{sample_name}.bed6_probes_df.tsv", sep='\t', float_format="%.3f", header=True, index=False)

    bed6_probes_df = bed6_probes_df.apply(lambda row: annotate_depth_region(row, depth_df), axis = 1)
    bed6_probes_df = bed6_probes_df.explode(["SAMPLE_ID", "EXON_SEQ", "EXON_SIZE"])
    bed6_probes_df = bed6_probes_df.dropna().reset_index(drop = True)

    # compute mean depth per gene
    bed6_probesByGene_df = bed6_probes_df.groupby(["GENE", "SAMPLE_ID"]).agg({"EXON_SEQ" : 'sum',
                                                                              "EXON_SIZE" : 'sum'}).reset_index()
    bed6_probesByGene_df["MEAN_GENE_DEPTH"] = bed6_probesByGene_df["EXON_SEQ"] / bed6_probesByGene_df["EXON_SIZE"]
    bed6_probesByGene_df_to_store = bed6_probesByGene_df.copy()
    bed6_probesByGene_df_to_store.columns = ["GENE", "SAMPLE_ID", "GENE_SEQ", "GENE_SIZE", "MEAN_GENE_DEPTH"]
    bed6_probesByGene_df_to_store.to_csv(f"{sample_name}.depth_per_gene_per_sample.tsv", sep = '\t', float_format="%.3f", header = True, index = False)

    bed6_probesByGene_df = bed6_probesByGene_df.drop(["EXON_SEQ", "EXON_SIZE"], axis = 'columns')

    average_depth_per_gene = bed6_probesByGene_df_to_store.groupby(by = ["GENE"] )[ ["GENE_SIZE", "MEAN_GENE_DEPTH"]].mean().reset_index()
    genes_list = list(average_depth_per_gene.sort_values(by = ["MEAN_GENE_DEPTH"], ascending = False)["GENE"].values)
    average_depth_per_gene.to_csv(f"{sample_name}.avgdepth_per_gene.tsv", sep = '\t', float_format="%.3f", header = True, index = False)

    return depth_df, samples_list, avgdepth_per_sample_names, bed6_probes_df, bed6_probesByGene_df, genes_list


def load_regions_and_compute_depths(sample_name, depth_df, region_file, bed6_probes_df=None):
    """
    Load the chromosomal regions file and compute the mean depth of the covered
    positions within each region, per sample and across all samples.

    The region file is a header-less TSV with columns: CHROM, START, END, REGION.
    Only the positions present in the depth file (i.e. the covered/sequenced
    positions) that fall within each region are used to compute the mean depth.

    If bed6_probes_df (the panel probes with GENE coordinates) is provided, an
    extra TSV is written listing the genes that fall within each region.
    """
    regions_df = pd.read_csv(region_file, sep="\t", header=None, names=["CHROM", "START", "END", "REGION"])
    regions_df["CHROM"] = regions_df["CHROM"].astype(str)
    # Ensure the CHROM column matches the depth file naming (with 'chr' prefix)
    regions_df["CHROM"] = regions_df["CHROM"].apply(lambda x: x if x.startswith("chr") else "chr" + x)

    # For each region, annotate the total sequenced bps (EXON_SEQ) and the number
    # of covered positions (EXON_SIZE) per SAMPLE_ID, reusing the same logic as for genes.
    regions_df = regions_df.apply(lambda row: annotate_depth_region(row, depth_df), axis=1)
    regions_df = regions_df.explode(["SAMPLE_ID", "EXON_SEQ", "EXON_SIZE"])
    regions_df = regions_df.dropna().reset_index(drop=True)

    # Keep only the region-sample pairs that actually have coverage (at least one
    # covered position). This drops regions that are not covered in any sample.
    regions_df = regions_df[regions_df["EXON_SIZE"] > 0].reset_index(drop=True)

    # Mean depth per region per sample
    regions_per_sample = regions_df.groupby(["REGION", "SAMPLE_ID"]).agg({"EXON_SEQ": "sum", "EXON_SIZE": "sum"}).reset_index()
    regions_per_sample["MEAN_REGION_DEPTH"] = regions_per_sample["EXON_SEQ"] / regions_per_sample["EXON_SIZE"]
    regions_per_sample_to_store = regions_per_sample[["REGION", "SAMPLE_ID", "EXON_SEQ", "EXON_SIZE", "MEAN_REGION_DEPTH"]].copy()
    regions_per_sample_to_store.columns = ["REGION", "SAMPLE_ID", "REGION_SEQ", "REGION_SIZE", "MEAN_REGION_DEPTH"]
    regions_per_sample_to_store.to_csv(f"{sample_name}.depth_per_region_per_sample.tsv", sep='\t', float_format="%.3f", header=True, index=False)

    # Mean depth per region across all samples
    regions_by_region = regions_df.groupby(["REGION"]).agg({"EXON_SEQ": "sum", "EXON_SIZE": "sum"}).reset_index()
    regions_by_region["MEAN_REGION_DEPTH"] = regions_by_region["EXON_SEQ"] / regions_by_region["EXON_SIZE"]
    regions_by_region_to_store = regions_by_region[["REGION", "EXON_SEQ", "EXON_SIZE", "MEAN_REGION_DEPTH"]].copy()
    regions_by_region_to_store.columns = ["REGION", "REGION_SEQ", "REGION_SIZE", "MEAN_REGION_DEPTH"]
    regions_by_region_to_store.to_csv(f"{sample_name}.avgdepth_per_region.tsv", sep='\t', float_format="%.3f", header=True, index=False)

    # Map the panel genes to the regions they fall within. A gene is assigned to a
    # region if any of its exons overlaps the region (same CHROM and coordinate overlap).
    if bed6_probes_df is not None:
        # Collapse the panel probes to one row per gene (chromosome + full span).
        genes_coords = bed6_probes_df.groupby(["CHROM", "GENE"]).agg(
            GENE_START=("START", "min"), GENE_END=("END", "max")).reset_index()
        genes_by_region = {}
        for _, reg in regions_df[["CHROM", "START", "END", "REGION"]].drop_duplicates().iterrows():
            genes_in_region = genes_coords[
                (genes_coords["CHROM"] == reg["CHROM"]) &
                (genes_coords["GENE_START"] <= reg["END"]) &
                (genes_coords["GENE_END"] >= reg["START"])
            ]["GENE"].unique()
            genes_by_region[reg["REGION"]] = sorted(genes_in_region)
        genes_per_region = pd.DataFrame(
            [{"REGION": r, "GENES": ",".join(genes)} for r, genes in genes_by_region.items()],
            columns=["REGION", "GENES"]
        )
        if len(genes_per_region) > 0:
            genes_per_region = genes_per_region.sort_values(by=["REGION"]).reset_index(drop=True)
        genes_per_region.to_csv(f"{sample_name}.genes_per_region.tsv", sep='\t', float_format="%.3f", header=True, index=False)

    # Order the regions by increasing (alphabetical) region name for plotting.
    regions_list = list(regions_by_region_to_store.sort_values(by=["REGION"], ascending=True)["REGION"].values)

    return regions_df, regions_per_sample, regions_list


def general_plotting(sample_name, samples_list, bed6_probesByGene_df, genes_list, regions_per_sample=None, regions_list=None, avgdepth_per_sample_names=None):

    with PdfPages(f'{sample_name}.depths.pdf') as pdf:

        ######
        ## Depth per sample
        ######
        if _section_allowed(len(samples_list) <= MAX_CATEGORIES_PER_PLOT,
                            "Depth per sample (boxplot)",
                            f"{len(samples_list)} samples > {MAX_CATEGORIES_PER_PLOT}"):
            fig, ax1 = plt.subplots(1, 1)
            fig.set_size_inches(min(max(0.5*len(samples_list), 10), 20), 5)
            sns.boxplot(data = bed6_probesByGene_df, x = "SAMPLE_ID", y = "MEAN_GENE_DEPTH", ax = ax1, order = samples_list, showfliers = False)
            sns.stripplot(data = bed6_probesByGene_df, x = "SAMPLE_ID", y = "MEAN_GENE_DEPTH", ax = ax1, order = samples_list, jitter = True,
                        alpha = 0.5, size = 4)
            ax1.set_title("Depth per SAMPLE_ID")
            ax1.tick_params(axis = 'x', labelrotation = 90)
            plt.tight_layout()
            pdf.savefig()
            plt.close()


        ######
        ## Depth per gene
        ######
        if _section_allowed(len(genes_list) <= MAX_CATEGORIES_PER_PLOT,
                            "Depth per gene (boxplot)",
                            f"{len(genes_list)} genes > {MAX_CATEGORIES_PER_PLOT}"):
            fig, ax2 = plt.subplots(1, 1)
            fig.set_size_inches(10, 5)
            sns.boxplot(data = bed6_probesByGene_df,
                        x = "GENE", y = "MEAN_GENE_DEPTH",
                        ax = ax2, showfliers = False,
                    )
            sns.stripplot(data = bed6_probesByGene_df,
                        x = "GENE", y = "MEAN_GENE_DEPTH",
                        ax = ax2,
                        jitter = True,
                        alpha = 0.5, size = 4)
            ax2.set_title("Depth per GENE")
            ax2.tick_params(axis = 'x', labelrotation = 90)
            plt.tight_layout()
            pdf.savefig()
            plt.close()


        ####
        ## PLOT DEPTHS per sample
        ####

        # Calculate the number of unique SAMPLE_ID values
        num_samples = len(samples_list)

        # Define the maximum number of SAMPLE_ID values per plot
        max_samples_per_plot = 30

        # Check if the number of SAMPLE_ID values exceeds the maximum
        if num_samples > max_samples_per_plot:
            # Calculate the number of plots needed
            num_plots = int(np.ceil(num_samples / max_samples_per_plot))

            # Split the SAMPLE_ID values into chunks for each plot
            sample_chunks = np.array_split(samples_list, num_plots)

            # Iterate over the sample chunks and create separate plots
            for i, sample_chunk in enumerate(sample_chunks):
                # Create a new figure and axis for each plot
                fig, ax = plt.subplots(figsize=(15, 5))

                # Filter the dataframe to include only the SAMPLE_ID values in the current chunk
                subset_df = bed6_probesByGene_df[bed6_probesByGene_df['SAMPLE_ID'].isin(sample_chunk)]

                # Plot the boxplot and stripplot for the current SAMPLE_ID values
                sns.boxplot(data=subset_df, x="SAMPLE_ID", y="MEAN_GENE_DEPTH", ax=ax, order=sample_chunk, showfliers = False)
                sns.stripplot(data=subset_df, x="SAMPLE_ID", y="MEAN_GENE_DEPTH", ax=ax, order=sample_chunk,
                            jitter=True, alpha=0.5, size=4)

                # Set plot title and rotate x-axis labels
                ax.set_title(f"Depth per SAMPLE_ID (Plot {i+1}/{num_plots})")
                ax.tick_params(axis='x', labelrotation=90)

                # Ensure tight layout and display the plot
                plt.tight_layout()
                pdf.savefig()
                plt.close()
        else:
            # If the number of SAMPLE_ID values does not exceed the maximum, plot all SAMPLE_ID values in a single plot
            fig, ax = plt.subplots(figsize=(15, 5))
            sns.boxplot(data=bed6_probesByGene_df, x="SAMPLE_ID", y="MEAN_GENE_DEPTH", ax=ax, order=samples_list, showfliers = False)
            sns.stripplot(data=bed6_probesByGene_df, x="SAMPLE_ID", y="MEAN_GENE_DEPTH", ax=ax, order=samples_list,
                        jitter=True, alpha=0.5, size=4)
            ax.set_title("Depth per SAMPLE_ID")
            ax.tick_params(axis='x', labelrotation=90)
            plt.tight_layout()
            pdf.savefig()
            plt.close()




        ####
        ## PLOT DEPTHS per gene
        ####

        # Calculate the number of genes
        num_genes = len(genes_list)


        # Define the maximum number of genes per plot
        max_genes_per_plot = 30

        # Check if the number of genes exceeds the maximum
        if num_genes > max_genes_per_plot:
            # Calculate the number of plots needed
            num_plots = int(np.ceil(num_genes / max_genes_per_plot))

            # Split the genes into chunks for each plot
            gene_chunks = np.array_split(genes_list, num_plots)

            # Iterate over the gene chunks and create separate plots
            for i, gene_chunk in enumerate(gene_chunks):
                # Create a new figure and axis for each plot
                fig, ax = plt.subplots(figsize=(15, 6))

                # Filter the dataframe to include only the genes in the current chunk
                subset_df = bed6_probesByGene_df[bed6_probesByGene_df['GENE'].isin(gene_chunk)]

                # Plot the boxplot and stripplot for the current genes
                sns.boxplot(data=subset_df, x="GENE", y="MEAN_GENE_DEPTH", ax=ax, order = gene_chunk, showfliers = False)
                sns.stripplot(data=subset_df, x="GENE", y="MEAN_GENE_DEPTH", ax=ax, order = gene_chunk,
                            jitter=True, alpha=0.5, size=4)

                # Set plot title and rotate x-axis labels
                ax.set_title(f"Depth per GENE (Plot {i+1}/{num_plots})")
                ax.tick_params(axis='x', labelrotation=90)

                # Ensure tight layout and display the plot
                plt.tight_layout()
                pdf.savefig()
                plt.close()

        else:
            # If the number of genes does not exceed the maximum, plot all genes in a single plot
            fig, ax = plt.subplots(figsize=(15, 6))
            sns.boxplot(data=bed6_probesByGene_df, x="GENE", y="MEAN_GENE_DEPTH", ax=ax, order = genes_list, showfliers = False)
            sns.stripplot(data=bed6_probesByGene_df, x="GENE", y="MEAN_GENE_DEPTH", ax=ax, order = genes_list,
                        jitter=True, alpha=0.5, size=4)
            ax.set_title("Depth per GENE")
            ax.tick_params(axis='x', labelrotation=90)
            plt.tight_layout()
            pdf.savefig()
            plt.close()



        ######
        ## Mean depth per gene
        ######
        if _section_allowed(len(samples_list) <= MAX_LINES_PER_PLOT and
                            len(samples_list) * len(genes_list) <= MAX_LINEPLOT_POINTS,
                            "Mean depth per gene (lineplot)",
                            f"{len(samples_list)} lines / {len(samples_list)*len(genes_list)} points too large"):
            plt.figure(figsize = (14,6))
            bed6_probesByGene_df_per_gene = bed6_probesByGene_df.copy()
            bed6_probesByGene_df_per_gene["GENE"] = pd.Categorical(bed6_probesByGene_df_per_gene["GENE"], categories=genes_list, ordered=True)
            ax = sns.lineplot(data = bed6_probesByGene_df_per_gene, x = "GENE", y = "MEAN_GENE_DEPTH", alpha = .7,
                            hue = "SAMPLE_ID",
                            hue_order = samples_list,
                            legend = False
                            )
            ax.set_title("Mean DEPTH per GENE")
            ax.tick_params(axis = 'x', labelrotation = 90)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
            del bed6_probesByGene_df_per_gene


        ######
        ## Mean depth per sample
        ######
        if _section_allowed(len(genes_list) <= MAX_LINES_PER_PLOT and
                            len(genes_list) * len(samples_list) <= MAX_LINEPLOT_POINTS,
                            "Mean depth per sample (lineplot)",
                            f"{len(genes_list)} lines / {len(genes_list)*len(samples_list)} points too large"):
            plt.figure(figsize = (14,6))
            bed6_probesByGene_df_per_sample = bed6_probesByGene_df.copy()
            bed6_probesByGene_df_per_sample["SAMPLE_ID"] = pd.Categorical(bed6_probesByGene_df_per_sample["SAMPLE_ID"], categories=samples_list, ordered=True)
            ax = sns.lineplot(data = bed6_probesByGene_df_per_sample, x = "SAMPLE_ID", y = "MEAN_GENE_DEPTH", alpha = .7,
                            hue = "GENE",
                            hue_order = genes_list,
                            legend = False
                            )
            ax.set_title("Mean DEPTH per SAMPLE_ID")
            ax.tick_params(axis = 'x', labelrotation = 90)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
            del bed6_probesByGene_df_per_sample


        # ######
        # ## Depth per gene per sample boxplots
        # ######

        # sns.set_theme(style='white')
        # g = sns.FacetGrid(data = bed6_probesByGene_df, col = "SAMPLE_ID", col_wrap = 4, height = 4, col_order = samples_list)
        # g.map(sns.boxplot, "GENE", "MEAN_GENE_DEPTH", showfliers = False, order = genes_list)
        # g.map(sns.stripplot, "GENE", "MEAN_GENE_DEPTH", jitter = True, order = genes_list, alpha = 0.5)
        # g.tick_params('x', labelrotation = 90)
        # plt.tight_layout()
        # pdf.savefig()
        # plt.close()

        # Add a heatmap: MEAN_GENE_DEPTH per gene per sample
        if _section_allowed(len(genes_list) * len(samples_list) <= MAX_HEATMAP_CELLS,
                            "Mean gene depth heatmap",
                            f"{len(genes_list)}x{len(samples_list)} cells > {MAX_HEATMAP_CELLS}"):
            heatmap_data = bed6_probesByGene_df.pivot(index="GENE", columns="SAMPLE_ID", values="MEAN_GENE_DEPTH")
            heatmap_data = heatmap_data.reindex(index=genes_list, columns=samples_list)
            heatmap_data = heatmap_data.astype(float)
            plt.figure(figsize=(min(max(10, 0.4*len(samples_list)), MAX_FIG_WIDTH),
                                min(max(8, 0.3*len(genes_list)), MAX_FIG_HEIGHT)))
            sns.heatmap(data = heatmap_data, cmap="viridis", cbar_kws={"label": "Mean Gene Depth"})
            plt.title("Mean Gene Depth per Gene per Sample")
            plt.xlabel("Sample ID")
            plt.ylabel("Gene")
            plt.tight_layout()
            pdf.savefig()
            plt.close()

            if (heatmap_data.shape[0] > 1 and heatmap_data.shape[1] > 1
                    and max(heatmap_data.shape) <= MAX_CLUSTERMAP_DIM):
                # Clustermap
                sns.clustermap(data=heatmap_data, cmap="viridis",
                               figsize=(min(max(10, 0.4*len(samples_list)), MAX_FIG_WIDTH),
                                        min(max(8, 0.3*len(genes_list)), MAX_FIG_HEIGHT)))
                pdf.savefig()
                plt.close()


        ######
        ## Depth per chromosomal region
        ######
        if regions_per_sample is not None and regions_list is not None and len(regions_list) > 0:

            # Depth per REGION (boxplot + stripplot)
            if _section_allowed(len(regions_list) <= MAX_REGIONS_PER_PLOT,
                                "Depth per region (boxplot)",
                                f"{len(regions_list)} regions > {MAX_REGIONS_PER_PLOT}"):
                fig, ax = plt.subplots(1, 1)
                fig.set_size_inches(min(max(0.5*len(regions_list), 10), 20), 5)
                sns.boxplot(data = regions_per_sample, x = "REGION", y = "MEAN_REGION_DEPTH", ax = ax, order = regions_list, showfliers = False)
                sns.stripplot(data = regions_per_sample, x = "REGION", y = "MEAN_REGION_DEPTH", ax = ax, order = regions_list, jitter = True,
                            alpha = 0.5, size = 4)
                ax.set_title("Depth per REGION")
                ax.tick_params(axis = 'x', labelrotation = 90)
                plt.tight_layout()
                pdf.savefig()
                plt.close()


            # Mean depth per REGION (lineplot per sample)
            if _section_allowed(len(samples_list) <= MAX_LINES_PER_PLOT and
                                len(samples_list) * len(regions_list) <= MAX_LINEPLOT_POINTS,
                                "Mean depth per region (lineplot)",
                                f"{len(samples_list)} lines / {len(samples_list)*len(regions_list)} points too large"):
                plt.figure(figsize = (14,6))
                regions_per_sample_plot = regions_per_sample.copy()
                regions_per_sample_plot["REGION"] = pd.Categorical(regions_per_sample_plot["REGION"], categories=regions_list, ordered=True)
                ax = sns.lineplot(data = regions_per_sample_plot, x = "REGION", y = "MEAN_REGION_DEPTH", alpha = .7,
                                hue = "SAMPLE_ID",
                                hue_order = samples_list,
                                legend = False
                                )
                ax.set_title("Mean DEPTH per REGION")
                ax.tick_params(axis = 'x', labelrotation = 90)
                plt.tight_layout()
                pdf.savefig()
                plt.close()
                del regions_per_sample_plot


            # Heatmap: MEAN_REGION_DEPTH per region per sample
            if _section_allowed(len(regions_list) * len(samples_list) <= MAX_HEATMAP_CELLS,
                                "Mean region depth heatmap",
                                f"{len(regions_list)}x{len(samples_list)} cells > {MAX_HEATMAP_CELLS}"):
                region_heatmap_data = regions_per_sample.pivot(index="REGION", columns="SAMPLE_ID", values="MEAN_REGION_DEPTH")
                region_heatmap_data = region_heatmap_data.reindex(index=regions_list, columns=samples_list)
                region_heatmap_data = region_heatmap_data.astype(float)
                plt.figure(figsize=(min(max(10, 0.4*len(samples_list)), MAX_FIG_WIDTH),
                                    min(max(8, 0.3*len(regions_list)), MAX_FIG_HEIGHT)))
                sns.heatmap(data = region_heatmap_data, cmap="viridis", cbar_kws={"label": "Mean Region Depth"})
                plt.title("Mean Region Depth per Region per Sample")
                plt.xlabel("Sample ID")
                plt.ylabel("Region")
                plt.tight_layout()
                pdf.savefig()
                plt.close()


            # Depth per REGION normalized by the sample's average depth
            regions_norm = regions_per_sample.merge(avgdepth_per_sample_names, on="SAMPLE_ID", how="left")
            regions_norm["MEAN_REGION_DEPTH_NORM"] = regions_norm["MEAN_REGION_DEPTH"] / regions_norm["avg_depth_sample"]

            # Boxplot + stripplot of the normalized depth per region
            if _section_allowed(len(regions_list) <= MAX_REGIONS_PER_PLOT,
                                "Depth per region normalized (boxplot)",
                                f"{len(regions_list)} regions > {MAX_REGIONS_PER_PLOT}"):
                fig, ax = plt.subplots(1, 1)
                fig.set_size_inches(min(max(0.5*len(regions_list), 10), 20), 5)
                sns.boxplot(data = regions_norm, x = "REGION", y = "MEAN_REGION_DEPTH_NORM", ax = ax, order = regions_list, showfliers = False)
                sns.stripplot(data = regions_norm, x = "REGION", y = "MEAN_REGION_DEPTH_NORM", ax = ax, order = regions_list, jitter = True,
                            alpha = 0.5, size = 4)
                ax.set_title("Depth per REGION (normalized by sample's average depth)")
                ax.tick_params(axis = 'x', labelrotation = 90)
                plt.tight_layout()
                pdf.savefig()
                plt.close()


            # Mean normalized depth per REGION (lineplot per sample)
            if _section_allowed(len(samples_list) <= MAX_LINES_PER_PLOT and
                                len(samples_list) * len(regions_list) <= MAX_LINEPLOT_POINTS,
                                "Mean depth per region normalized (lineplot)",
                                f"{len(samples_list)} lines / {len(samples_list)*len(regions_list)} points too large"):
                plt.figure(figsize = (14,6))
                regions_norm_plot = regions_norm.copy()
                regions_norm_plot["REGION"] = pd.Categorical(regions_norm_plot["REGION"], categories=regions_list, ordered=True)
                ax = sns.lineplot(data = regions_norm_plot, x = "REGION", y = "MEAN_REGION_DEPTH_NORM", alpha = .7,
                                hue = "SAMPLE_ID",
                                hue_order = samples_list,
                                legend = False
                                )
                ax.set_title("Mean DEPTH per REGION (normalized by sample's average depth)")
                ax.tick_params(axis = 'x', labelrotation = 90)
                plt.tight_layout()
                pdf.savefig()
                plt.close()
                del regions_norm_plot


            # Heatmap: normalized MEAN_REGION_DEPTH per region per sample
            if _section_allowed(len(regions_list) * len(samples_list) <= MAX_HEATMAP_CELLS,
                                "Mean region depth normalized heatmap",
                                f"{len(regions_list)}x{len(samples_list)} cells > {MAX_HEATMAP_CELLS}"):
                region_norm_heatmap_data = regions_norm.pivot(index="REGION", columns="SAMPLE_ID", values="MEAN_REGION_DEPTH_NORM")
                region_norm_heatmap_data = region_norm_heatmap_data.reindex(index=regions_list, columns=samples_list)
                region_norm_heatmap_data = region_norm_heatmap_data.astype(float)
                plt.figure(figsize=(min(max(10, 0.4*len(samples_list)), MAX_FIG_WIDTH),
                                    min(max(8, 0.3*len(regions_list)), MAX_FIG_HEIGHT)))
                sns.heatmap(data = region_norm_heatmap_data, cmap="viridis", cbar_kws={"label": "Mean Region Depth (normalized)"})
                plt.title("Mean Region Depth per Region per Sample (normalized by sample's average depth)")
                plt.xlabel("Sample ID")
                plt.ylabel("Region")
                plt.tight_layout()
                pdf.savefig()
                plt.close()


def process_within_gene_depths(sample_name, depth_df, bed6_probes_df, bed6_probesByGene_df, genes_list, samples_list, avgdepth_per_sample_names):
    """
    Function to process and plot within-gene depths.
    """
    # Annotate MEAN_GENE_DEPTH normalized by average MEAN_GENE_DEPTH per GENE and SAMPLE_ID
    bed6_probes_df["EXON_DEPTH"] = bed6_probes_df["EXON_SEQ"] / bed6_probes_df["EXON_SIZE"]
    bed6_probes_df["exon_cov_norm_gene_depth"] = bed6_probes_df.apply(
        lambda row: row.EXON_DEPTH / bed6_probesByGene_df.loc[
            (bed6_probesByGene_df["SAMPLE_ID"] == row["SAMPLE_ID"]) &
            (bed6_probesByGene_df.GENE == row.GENE), "MEAN_GENE_DEPTH"].values[0], axis=1)

    # Build dataframe of unique positions
    unique_pos_df = depth_df[["CHROM", "POS"]].drop_duplicates()

    exons_list_col = []
    for ind, row in unique_pos_df.iterrows():
        vals = bed6_probes_df.loc[
            (bed6_probes_df["CHROM"] == row.CHROM) &
            (bed6_probes_df["START"] <= row["POS"]) &
            (row["POS"] <= bed6_probes_df["END"]), "EXON"].values
        exons_list_col.append(vals[0] if len(vals) > 0 else None)

    unique_pos_df["EXON"] = exons_list_col
    pos_to_fake_pos = {y: x for x, y in enumerate(pd.unique(unique_pos_df["POS"]))}
    unique_pos_df["fake_pos"] = unique_pos_df["POS"].replace(pos_to_fake_pos)
    unique_pos_df.to_csv(f"{sample_name}.positions_covered.with_fake_pos.tsv", sep='\t', header=True, index=False)



    with PdfPages(f'{sample_name}.depths_within_gene.pdf') as pdf:

        sns.set_theme(style='white')
        g = sns.FacetGrid(data = bed6_probes_df, col = "SAMPLE_ID", col_wrap = 4, height = 4, col_order = samples_list)
        g.map(sns.lineplot, "EXON", "EXON_DEPTH", alpha = .7)
        g.set_xticklabels("", "")
        # g.add_legend()
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        for g in genes_list:
            plt.figure(figsize = (12, 6))
            data = bed6_probes_df.loc[bed6_probes_df.GENE == g]
            ax = sns.lineplot(data = data, x = "EXON", y = "EXON_DEPTH", alpha = .7,
                            #   hue = "PROJECT_NAME",
                                hue = "SAMPLE_ID", hue_order = samples_list,
                            # palette = sample_colors,
                            legend = False
                            )
            ax.axhline(data.EXON_DEPTH.mean(), linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
            ax.axhspan(data.EXON_DEPTH.mean(), data.EXON_DEPTH.mean()+2*data.EXON_DEPTH.std(), facecolor = '0.8')
            ax.axhspan(data.EXON_DEPTH.mean(), max(0, data.EXON_DEPTH.mean()-2*data.EXON_DEPTH.std()), facecolor = '0.8')
            ax.set_title(g)
            ax.set_xticks([])
            # ax.legend(fontsize = 8, bbox_to_anchor= (1.05,1))
            plt.tight_layout()
            

            pdf.savefig()
            plt.close()






        for g in genes_list:
            fig, (ax1, ax2) = plt.subplots(1, 2)
            fig.set_size_inches(17, 5)

            data = bed6_probes_df.loc[bed6_probes_df.GENE == g]

            sns.lineplot(data = data, x = "EXON", y = "EXON_DEPTH", alpha = .7, hue = "SAMPLE_ID", hue_order = samples_list, ax = ax1, legend = False)
            ax1.axhline(data.EXON_DEPTH.mean(), linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
            ax1.axhspan(data.EXON_DEPTH.mean(), data.EXON_DEPTH.mean()+2*data.EXON_DEPTH.std(), facecolor = '0.8')
            ax1.axhspan(data.EXON_DEPTH.mean(), data.EXON_DEPTH.mean()-2*data.EXON_DEPTH.std(), facecolor = '0.8')
            ax1.set_title(g)
            ax1.set_xticks([])

            sns.lineplot(data = data, x = "EXON", y = "exon_cov_norm_gene_depth", alpha = .7, hue = "SAMPLE_ID", hue_order = samples_list, ax = ax2,
                        legend = False)
            ax2.axhline(1, linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
            ax2.set_title(f"{g} (normalized by GENE's mean EXON_DEPTH)")
            ax2.set_xticks([])
            # ax2.legend(fontsize = 8, bbox_to_anchor= (1.25,1))
            plt.tight_layout()
            
            pdf.savefig()
            plt.close()


        ## read again the uniq_pos_df
        unique_pos_df = pd.read_csv(f"{sample_name}.positions_covered.with_fake_pos.tsv", sep = '\t', header = 0)

        depth_df_melted =  pd.melt(depth_df, id_vars=["CHROM", "POS"], value_vars=samples_list,
                                                    var_name='SAMPLE_ID', value_name='READS')
        depth_df_exonlab = depth_df_melted.merge(unique_pos_df, on = ["CHROM", "POS"]).drop("POS", axis = 'columns')
        depth_df_exonlab.columns = ["CHROM", "SAMPLE_ID", "READS", "EXON", "fake_pos"]
        depth_df_exonlab.to_csv(f"{sample_name}.depth_df_exonlab.tsv", sep = '\t', header = True, index = False)


        depth_df_exonlab_exondepth = depth_df_exonlab.merge(bed6_probes_df, on = ["SAMPLE_ID", "EXON", "CHROM"]).drop("CHROM", axis = 'columns')
        depth_df_exonlab_exondepth.to_csv(f"{sample_name}.depth_df_exonlab_exondepth.tsv", sep = '\t', header = True, index = False)



        bed6_probesByGene_df.columns = ["GENE", "SAMPLE_ID", "GENE_DEPTH"]
        depth_df_exonlab_exondepth_genedepth = depth_df_exonlab_exondepth.merge(bed6_probesByGene_df, on = ["GENE", "SAMPLE_ID"])
        print(depth_df_exonlab_exondepth_genedepth.shape)
        print(depth_df_exonlab_exondepth_genedepth.head())
        depth_df_exonlab_exondepth_genedepth.to_csv(f"{sample_name}.depth_df_exonlab_exondepth_genedepth.tsv",
                                                            sep = '\t', header = True, index = False)


        depth_df_exonlab_exondepth_genedepth_sample_depth = depth_df_exonlab_exondepth_genedepth.merge(avgdepth_per_sample_names, on = 'SAMPLE_ID')
        print(depth_df_exonlab_exondepth_genedepth_sample_depth.shape)
        print(depth_df_exonlab_exondepth_genedepth_sample_depth.head())
        depth_df_exonlab_exondepth_genedepth_sample_depth.to_csv(f"{sample_name}.depth_df_exonlab_exondepth_genedepth_sample_depth.tsv",
                                                                    sep = '\t', header = True, index = False)
        # depth_df_exonlab_exondepth_genedepth_sample_depth

        depth_df_exonlab_exondepth_genedepth_sample_depth = depth_df_exonlab_exondepth_genedepth_sample_depth[ (depth_df_exonlab_exondepth_genedepth_sample_depth["EXON_DEPTH"] > 0)
                                                                                                                & (depth_df_exonlab_exondepth_genedepth_sample_depth["GENE_DEPTH"] > 0)
                                                                                                                ].reset_index(drop = True)



        depth_df_exonlab_exondepth_genedepth_sample_depth["cov_norm_exon_depth"] = depth_df_exonlab_exondepth_genedepth_sample_depth["READS"] \
                                                                                        / depth_df_exonlab_exondepth_genedepth_sample_depth["EXON_DEPTH"]
        depth_df_exonlab_exondepth_genedepth_sample_depth["cov_norm_gene_depth"] = depth_df_exonlab_exondepth_genedepth_sample_depth["READS"] \
                                                                                        / depth_df_exonlab_exondepth_genedepth_sample_depth["GENE_DEPTH"]
        depth_df_exonlab_exondepth_genedepth_sample_depth["cov_norm_sample_depth"] = depth_df_exonlab_exondepth_genedepth_sample_depth["READS"] \
                                                                                        / depth_df_exonlab_exondepth_genedepth_sample_depth["avg_depth_sample"]

        depth_df_exonlab_exondepth_genedepth_sample_depth.to_csv(f"{sample_name}.depth_df_exonlab_exondepth_genedepth_sample_depth_averages.tsv",
                                                                    sep = '\t', header = True, index = False)

        print(depth_df_exonlab_exondepth_genedepth_sample_depth.head())


        # plot no normalized and mean EXON_DEPTH-normalized together
        for g in genes_list:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
            fig.set_size_inches(20, 12)  # Adjust the figure size as needed

            data = depth_df_exonlab_exondepth_genedepth_sample_depth.loc[depth_df_exonlab_exondepth_genedepth_sample_depth.GENE == g]

            try :
                sing_sampl_data = data[data["SAMPLE_ID"] == data["SAMPLE_ID"].iloc[0]]

                # Find the break points for the intron-exon boundary
                prev_start = ""
                exon_boundary_pos = []
                for ind, r in sing_sampl_data.iloc[1:,:].iterrows():
                    if prev_start != r["EXON"]:
                        exon_boundary_pos.append(r["fake_pos"])
                    prev_start = r["EXON"]
            except:
                exon_boundary_pos = []


            # Plot READS vs. fake_pos
            sns.lineplot(data=data, x="fake_pos", y="READS", alpha=.7, hue="SAMPLE_ID", hue_order=samples_list, ax=ax1, legend=False)
            ax1.axhline(data.EXON_DEPTH.mean(), linestyle="--", linewidth=2, color="black", alpha=0.5)
            ax1.axhspan(data.EXON_DEPTH.mean() - 2 * data.EXON_DEPTH.std(), data.EXON_DEPTH.mean() + 2 * data.EXON_DEPTH.std(),
                        facecolor='0.8')
            for ver_pos in exon_boundary_pos:
                ax1.axvline(x=ver_pos, linestyle="--", color="black", alpha=0.3)
            ax1.set_title(f"{g} positions - coverage")
            ax1.set_xticks([])

            # Plot cov_norm_exon_depth vs. fake_pos
            sns.lineplot(data=data, x="fake_pos", y="cov_norm_exon_depth", alpha=.7, hue="SAMPLE_ID", hue_order=samples_list, ax=ax2,
                        legend = False)
            ax2.axhline(1, linestyle="--", linewidth=2, color="black", alpha=0.5)
            for ver_pos in exon_boundary_pos:
                ax2.axvline(x=ver_pos, linestyle="--", color="black", alpha=0.3)
            ax2.set_title(f"{g} positions (normalized by exons's mean DEPTH)")
            ax2.set_xticks([])
    #        ax2.legend(fontsize=8, bbox_to_anchor=(1.25, 1))


            # Plot cov_norm_gene_depth vs. fake_pos
            sns.lineplot(data=data, x="fake_pos", y="cov_norm_gene_depth", alpha=.7, hue="SAMPLE_ID", hue_order=samples_list, ax=ax3, legend=False)
            ax3.axhline(1, linestyle="--", linewidth=2, color="black", alpha=0.5)
            for ver_pos in exon_boundary_pos:
                ax3.axvline(x=ver_pos, linestyle="--", color="black", alpha=0.3)
            ax3.set_title(f"{g} positions (normalized by gene's mean DEPTH)")
            ax3.set_xticks([])

            # Plot cov_norm_sample_depth vs. fake_pos
            sns.lineplot(data=data, x="fake_pos", y="cov_norm_sample_depth", alpha=.7, hue="SAMPLE_ID", hue_order=samples_list, ax=ax4, legend=False)
            ax4.axhline(1, linestyle="--", linewidth=2, color="black", alpha=0.5)
            for ver_pos in exon_boundary_pos:
                ax4.axvline(x=ver_pos, linestyle="--", color="black", alpha=0.3)
            ax4.set_title(f"{g} positions (normalized by sample's mean DEPTH)")
            ax4.set_xticks([])

            plt.tight_layout()
            pdf.savefig()
            plt.close()


def process_depths(sample_name, depth_file, panel_bed6_file, panel_name, plot_within_gene, region_file=None):
    """
    Main function to process depths and generate plots.
    """
    depth_df, samples_list, avgdepth_per_sample_names, bed6_probes_df, bed6_probesByGene_df, genes_list = load_depths_and_panel(sample_name, depth_file, panel_bed6_file)

    # Compute the mean depth per chromosomal region (if a region file is provided)
    regions_per_sample = None
    regions_list = None
    if region_file:
        _, regions_per_sample, regions_list = load_regions_and_compute_depths(sample_name, depth_df, region_file, bed6_probes_df)

    general_plotting(sample_name, samples_list, bed6_probesByGene_df, genes_list, regions_per_sample, regions_list, avgdepth_per_sample_names)

    ####
    ## Until here all the plots and all the information was at the gene or sample level.
    # Now it starts to contain within gene information, exons, normalized scores, ...
    ####
    if plot_within_gene:
        process_within_gene_depths(sample_name, depth_df, bed6_probes_df, bed6_probesByGene_df, genes_list, samples_list, avgdepth_per_sample_names)



@click.command()
@click.option('--sample_name', type=str, required=True, help='Name of the sample being processed.')
@click.option('--depth_file', type=click.Path(exists=True), required=True, help='Input depth file.')
@click.option('--panel_bed6_file', type=click.Path(exists=True), required=True, help='Input BED6 file for the panel.')
@click.option('--panel_name', type=str, required=True, help='Name of the panel.')
@click.option('--plot_within_gene', type=bool, default=False, help='Whether to plot within-gene information.')
@click.option('--region_file', type=click.Path(exists=True), default=None, help='Optional TSV file defining chromosomal regions (CHROM, START, END, REGION) to compute the mean depth per region.')
def main(sample_name, depth_file, panel_bed6_file, panel_name, plot_within_gene, region_file):
    """
    CLI entry point for processing depths and generating plots.
    """
    process_depths(sample_name, depth_file, panel_bed6_file, panel_name, plot_within_gene, region_file)


if __name__ == '__main__':
    main()
