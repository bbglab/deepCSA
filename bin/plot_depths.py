#!/usr/local/bin/python


# import click
import sys
# import os
# import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


# from copy import deepcopy

# import math
# from matplotlib.patches import Circle
# from matplotlib.collections import PatchCollection
# import matplotlib.patches as Patch, mpatches


# import re
# import glob
# import itertools


sample_name  = sys.argv[1]
depth_file   = sys.argv[2]
panel_bed6_file  = sys.argv[3]
panel_name   = sys.argv[4]

#json_groups  = sys.argv[5]
#req_plots    = sys.argv[6]



#maf_df = pd.read_csv(f"{mut_file}", sep = "\t", header = 0)


#unique_samples = maf_df[["SAMPLE_ID", "PROJECT_NAME"]].drop_duplicates().reset_index(drop = True)
#print(unique_samples.shape)


# panel = ["ARID1A", "BAP1", "MTOR",
#          "PBRM1", "PIK3CA", "PTEN",
#          "SETD2", "TP53", "VHL"]

# colors = ["#2f4f4f", "#7f0000", "#191970",
#           "#006400", "#bdb76b", "#ff0000",
#           "#ffa500", "#ffff00", "#0000cd"
#           # ,
#           # "#00ff00", "#00fa9a", "#da70d6",
#           # "#ff00ff", "#1e90ff", "#87ceeb",
#           # "#ffb6c1"
#          ]

# gene2color = dict(zip(panel, colors))


depth_df = pd.read_csv(f"{depth_file}", sep = "\t", header = 0)
depth_df = depth_df.drop(["CONTEXT"], axis = 'columns')
print(depth_df.head())

samples_list = depth_df.columns[~depth_df.columns.isin(["CHROM","POS"])]

stats_per_sample = pd.DataFrame(depth_df[samples_list].describe())
stats_per_sample.to_csv(f"{sample_name}.depth_per_sample.stats.tsv", sep = '\t', header = True, index = True)

avgdepth_per_sample = pd.DataFrame(depth_df[samples_list].mean().T)
avgdepth_per_sample.columns = ["avg_depth_sample"]
avgdepth_per_sample.to_csv(f"{sample_name}.avgdepth_per_sample.tsv", sep = '\t', header = True, index = True)
avgdepth_per_sample_names = avgdepth_per_sample.reset_index()
avgdepth_per_sample_names.columns = ["SAMPLE_ID", "avg_depth_sample"]

####
# load BED file with probes info
####
bed6_probes_df = pd.read_csv(panel_bed6_file, sep = "\t", header = 0).reset_index()
bed6_probes_df["EXON"] = bed6_probes_df["ELEMENT"] + "_" + bed6_probes_df["index"].astype(str)
bed6_probes_df = bed6_probes_df[["CHROMOSOME", "START", "END", "ELEMENT", "EXON"]]
bed6_probes_df.columns = ["CHROM", "START", "END", "GENE", "EXON"]
bed6_probes_df["CHROM"] = 'chr' + bed6_probes_df["CHROM"].astype(str)

# store bed6_probes_df
bed6_probes_df.to_csv(f"{sample_name}.bed6_probes_df.tsv", sep = '\t', header = True, index = False)



####
# Build dataframe of unique positions
####
unique_pos_df = depth_df[["CHROM", "POS"]].drop_duplicates()

exons_list_col = []
for ind, row in unique_pos_df.iterrows():
    vals = bed6_probes_df.loc[(bed6_probes_df["CHROM"] == row.CHROM)
                                & (bed6_probes_df["START"] <= row["POS"])
                                & (row["POS"] <= bed6_probes_df["END"]), "EXON"].values
    if len(vals) > 0:
        exons_list_col.append(vals[0])
    else:
        print(row)
        exons_list_col.append(None)

unique_pos_df["EXON"] = exons_list_col

pos_to_fake_pos = {  y : x  for x, y in enumerate(pd.unique(unique_pos_df["POS"])) }

unique_pos_df["fake_pos"] = unique_pos_df["POS"].replace(pos_to_fake_pos)

# fake position means that all positions in the dataframe will be right one
# after the other according to the fake_pos columns
unique_pos_df.to_csv(f"{sample_name}.positions_covered.with_fake_pos.tsv", sep = '\t', header = True, index = False)
del unique_pos_df




def annotate_avgdepth_region(row, dp_df):

    # look in the DEPTH file for the positions whithin this region
    dp_df_f = dp_df.loc[(dp_df.CHROM == row.CHROM) & (dp_df.POS >= row.START) & (dp_df.POS <= row.END)]

    # annotate mean EXON_DEPTH per SAMPLE_ID
    samples = dp_df_f.columns[~dp_df_f.columns.isin(["CHROM","POS"])]
    depths = []
    size = []
    for s in samples:
        depths.append(dp_df_f.loc[:,s].mean())
        size.append(dp_df_f.loc[:,s].count())

    row["EXON_DEPTH"] = depths
    row["SAMPLE_ID"] = list(samples)
    row["EXON_SIZE"] = size


    return row


bed6_probes_df = bed6_probes_df.apply(lambda row: annotate_avgdepth_region(row, depth_df), axis = 1)
bed6_probes_df = bed6_probes_df.drop(["START", "END"], axis = 'columns').explode(["EXON_DEPTH", "SAMPLE_ID", "EXON_SIZE"])
bed6_probes_df = bed6_probes_df.dropna().reset_index(drop = True)



# TODO FIXME
# this EXON_DEPTH is not fully correct since we are doing the mean of exons which might not all have the same length
bed6_probes_df["EXON_SEQ"] = bed6_probes_df["EXON_DEPTH"] * bed6_probes_df["EXON_SIZE"]
bed6_probesByGene_df = bed6_probes_df.groupby(["GENE", "SAMPLE_ID"]).agg({"EXON_SEQ" : 'sum', "EXON_DEPTH" : 'mean',  "EXON_SIZE" : 'sum'}).reset_index()
bed6_probesByGene_df["MEAN_GENE_DEPTH"] = bed6_probesByGene_df["EXON_SEQ"] / bed6_probesByGene_df["EXON_SIZE"]
bed6_probesByGene_df.to_csv(f"{sample_name}.depth_per_gene_per_sample.tsv", sep = '\t', header = True, index = False)

bed6_probesByGene_df = bed6_probesByGene_df.drop(["EXON_DEPTH", "EXON_SEQ", "EXON_SIZE"], axis = 'columns')

# annotate MEAN_GENE_DEPTH normalized by average MEAN_GENE_DEPTH per GENE and SAMPLE_ID
bed6_probes_df["exon_cov_norm_gene_depth"] = bed6_probes_df.apply(
    lambda row: row.EXON_DEPTH / bed6_probesByGene_df.loc[(bed6_probesByGene_df["SAMPLE_ID"] == row["SAMPLE_ID"]) &
                                                        (bed6_probesByGene_df.GENE == row.GENE), "MEAN_GENE_DEPTH"].values[0], axis = 1)






with PdfPages(f'{sample_name}.depths.pdf') as pdf:


    sns.set_theme(style='white')
    g = sns.FacetGrid(data = bed6_probes_df, col = "SAMPLE_ID", col_wrap = 4, height = 4, col_order = samples_list)
    g.map(sns.lineplot, "EXON", "EXON_DEPTH", alpha = .7)
    g.set_xticklabels("", "")
    # g.add_legend()
    plt.tight_layout()
    plt.show()

    pdf.savefig()
    plt.close()

    genes = bed6_probes_df.GENE.unique()
    for g in genes:
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
        plt.show()

        pdf.savefig()
        plt.close()




# plt.figure(figsize = (20, 5))
# ax = sns.boxplot(data = bed6_probes_df, x = "GENE", y = "EXON_DEPTH", hue = "SAMPLE_ID", order = panel, hue_order = samples_list,
#                  # palette = sample_colors
#                 )
# ax.legend(fontsize = 8, bbox_to_anchor= (1,1))


# plt.figure(figsize = (20, 5))
# ax = sns.boxplot(data = bed6_probes_df, x = "SAMPLE_ID", y = "EXON_DEPTH", hue = "GENE", hue_order = panel, palette = colors, order = samples_list)
# ax.legend(fontsize = 8, bbox_to_anchor= (1,1))
# plt.xticks(rotation = 90)
# plt.show()



    sns.set_theme(style='white')
    g = sns.FacetGrid(data = bed6_probes_df, col = "SAMPLE_ID", col_wrap = 4, height = 4, col_order = samples_list)
    g.map(sns.boxplot, "GENE", "EXON_DEPTH", showfliers = False)
    g.map(sns.stripplot, "GENE", "EXON_DEPTH", jitter = True, alpha = 0.5)
    g.tick_params('x', labelrotation = 90)
    # g.add_legend()
    plt.tight_layout()
    plt.show()
    pdf.savefig()
    plt.close()



    ######
    ## Depth per sample
    ######
    fig, ax1 = plt.subplots(1, 1)
    fig.set_size_inches(10, 5)
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
    fig, ax2 = plt.subplots(1, 1)
    fig.set_size_inches(10, 5)
    sns.boxplot(data = bed6_probesByGene_df, x = "GENE", y = "MEAN_GENE_DEPTH", ax = ax2, showfliers = False,
                # order = panel, palette = colors
            )
    sns.stripplot(data = bed6_probesByGene_df, x = "GENE", y = "MEAN_GENE_DEPTH", ax = ax2,
                    # order = panel, palette = colors
                # hue = "PROJECT_NAME",
                # palette = colors,
                jitter = True,
                alpha = 0.5, size = 4)
    ax2.set_title("Depth per GENE")
    ax2.tick_params(axis = 'x', labelrotation = 90)
    plt.tight_layout()
    pdf.savefig()
    plt.close()



    # same plot as above but collapsing by patient-GENE mean instead of plotting with the EXON EXON_DEPTH per patient and GENE
    # bed6_probesByGene_df = bed6_probes_df.groupby(["GENE", "SAMPLE_ID", "PROJECT_NAME"])["EXON_DEPTH"].mean().reset_index()


    # Calculate the number of unique SAMPLE_ID values
    num_samples = bed6_probesByGene_df['SAMPLE_ID'].nunique()

    # Define the maximum number of SAMPLE_ID values per plot
    max_samples_per_plot = 30

    # Check if the number of SAMPLE_ID values exceeds the maximum
    if num_samples > max_samples_per_plot:
        # Calculate the number of plots needed
        num_plots = int(np.ceil(num_samples / max_samples_per_plot))

        # Split the SAMPLE_ID values into chunks for each plot
        sample_chunks = np.array_split(bed6_probesByGene_df['SAMPLE_ID'].unique(), num_plots)

        # Iterate over the sample chunks and create separate plots
        for i, sample_chunk in enumerate(sample_chunks):
            # Create a new figure and axis for each plot
            fig, ax = plt.subplots(figsize=(10, 5))

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
        fig, ax = plt.subplots(figsize=(10, 5))
        sns.boxplot(data=bed6_probesByGene_df, x="SAMPLE_ID", y="MEAN_GENE_DEPTH", ax=ax, order=samples_list)
        sns.stripplot(data=bed6_probesByGene_df, x="SAMPLE_ID", y="MEAN_GENE_DEPTH", ax=ax, order=samples_list,
                    jitter=True, alpha=0.5, size=4)
        ax.set_title("Depth per SAMPLE_ID")
        ax.tick_params(axis='x', labelrotation=90)
        plt.tight_layout()
        pdf.savefig()
        plt.close()






    ####
    # Copied chunk
    ####

    # Calculate the number of genes
    unique_genes = sorted(bed6_probesByGene_df['GENE'].unique())
    num_genes = len(unique_genes)


    # Define the maximum number of genes per plot
    max_genes_per_plot = 30

    # Check if the number of genes exceeds the maximum
    if num_genes > max_genes_per_plot:
        # Calculate the number of plots needed
        num_plots = int(np.ceil(num_genes / max_genes_per_plot))

        # Split the genes into chunks for each plot
        gene_chunks = np.array_split(bed6_probesByGene_df['GENE'].unique(), num_plots)

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
        sns.boxplot(data=bed6_probesByGene_df, x="GENE", y="MEAN_GENE_DEPTH", ax=ax, order = unique_genes, showfliers = False)
        sns.stripplot(data=bed6_probesByGene_df, x="GENE", y="MEAN_GENE_DEPTH", ax=ax, order = unique_genes,
                    jitter=True, alpha=0.5, size=4)
        ax.set_title("Depth per GENE")
        ax.tick_params(axis='x', labelrotation=90)
        plt.tight_layout()
        pdf.savefig()
        plt.close()



    ######
    ## Mean depth per gene
    ######
    plt.figure(figsize = (14,6))
    ax = sns.lineplot(data = bed6_probesByGene_df, x = "GENE", y = "MEAN_GENE_DEPTH", alpha = .7,
                    hue = "SAMPLE_ID",
                    hue_order = samples_list,
                    legend = False)
    ax.set_title("Mean DEPTH per GENE")
    # ax.legend(fontsize = 8, bbox_to_anchor= (1,1))
    ax.tick_params(axis = 'x', labelrotation = 90)
    plt.tight_layout()
    plt.show()
    pdf.savefig()
    plt.close()



    ######
    ## Mean depth per sample
    ######
    plt.figure(figsize = (14,6))
    ax = sns.lineplot(data = bed6_probesByGene_df, x = "SAMPLE_ID", y = "MEAN_GENE_DEPTH", alpha = .7,
                    hue = "GENE",
                    hue_order = unique_genes,
                    # hue_order = panel, palette = colors,
                    legend = False
                    )
    ax.set_title("Mean DEPTH per SAMPLE_ID")
    # ax.legend(fontsize = 8, bbox_to_anchor= (1,1))
    ax.tick_params(axis = 'x', labelrotation = 90)
    plt.tight_layout()
    plt.show()
    pdf.savefig()
    plt.close()



    # plot no normalized and mean EXON_DEPTH-normalized together
    genes = bed6_probes_df.GENE.unique()
    for g in genes:
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
        plt.show()
        pdf.savefig()
        plt.close()



    ####
    ## Until here all the plots and all the information is by exon
    ####


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
    genes = depth_df_exonlab_exondepth_genedepth_sample_depth.GENE.unique()
    for g in genes:
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



#     # In[99]:


#     # plot no normalized and mean DEPTH-normalized together
#     genes = depth_df_exonlab_exondepth_genedepth_sample_depth.GENE.unique()
#     for g in genes:
#         fig, (ax1, ax2) = plt.subplots(1, 2)
#         fig.set_size_inches(17, 5)

#         data = depth_df_exonlab_exondepth_genedepth_sample_depth.loc[depth_df_exonlab_exondepth_genedepth_sample_depth.GENE == g]

#         sns.lineplot(data = data, x = "fake_pos", y = "READS", alpha = .7, hue = "SAMPLE_ID", hue_order = samples_list, ax = ax1, legend = False)
#         ax1.axhline(data.DEPTH.mean(), linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#         ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()+2*data.DEPTH.std(), facecolor = '0.8')
#         ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()-2*data.DEPTH.std(), facecolor = '0.8')
#         ax1.set_title(f"{g} positions")
#         ax1.set_xticks([])

#         sns.lineplot(data = data, x = "fake_pos", y = "cov_norm_gene_depth", alpha = .7, hue = "SAMPLE_ID", hue_order = samples_list, ax = ax2)
#         ax2.axhline(1, linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#         ax2.set_title(f"{g} positions (normalized by GENE's mean DEPTH)")
#         ax2.set_xticks([])
#         ax2.legend(fontsize = 8, bbox_to_anchor= (1.25,1))

#         plt.tight_layout()
#         plt.show()
#         pdf.savefig()
#         plt.close()


#     # In[121]:


#     # plot no normalized and mean DEPTH-normalized together
#     genes = depth_df_exonlab_exondepth_genedepth_sample_depth.GENE.unique()
#     for g in genes:
#         fig, (ax1, ax2) = plt.subplots(1, 2)
#         fig.set_size_inches(17, 5)

#         data = depth_df_exonlab_exondepth_genedepth_sample_depth.loc[depth_df_exonlab_exondepth_genedepth_sample_depth.GENE == g].reset_index(drop = True)

#         sns.lineplot(data = data, x = "fake_pos", y = "READS", alpha = .7, hue = "SAMPLE_ID", hue_order = samples_list, ax = ax1, legend = False)
#         ax1.axhline(data.DEPTH.mean(), linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#         ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()+2*data.DEPTH.std(), facecolor = '0.8')
#         ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()-2*data.DEPTH.std(), facecolor = '0.8')
#         ax1.set_title(f"{g} positions")
#         ax1.set_xticks([])

#         sns.lineplot(data = data, x = "fake_pos", y = "cov_norm_sample_depth", alpha = .7, hue = "SAMPLE_ID", hue_order = samples_list, ax = ax2)
#         ax2.axhline(1, linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#         ax2.set_title(f"{g} positions (normalized by SAMPLE_ID's mean DEPTH)")
#         ax2.set_xticks([])

#         sing_sampl_data = data[data["SAMPLE_ID"] == data["SAMPLE_ID"][0]]
#         prev_start = ""
#         for ind, r in sing_sampl_data.iloc[1:,:].iterrows():
#             if prev_start != r["START"]:
#                 ax2.axvline(x = r["fake_pos"], linestyle = "--", color = "black", alpha = 0.3)
#             prev_start = r["START"]

#         # ax2.legend(fontsize = 8, bbox_to_anchor= (1.25,1))
#         ax2.legend('', frameon=False)

#         plt.tight_layout()
#         plt.show()
#         pdf.savefig()
#         plt.close()


#     # ## __Splitting by project__

#     # In[102]:


#     # plot no normalized and mean DEPTH-normalized together
#     genes = depth_df_exonlab_exondepth_genedepth_sample_depth.GENE.unique()
#     projects = depth_df_exonlab_exondepth_genedepth_sample_depth.PROJECT_NAME.unique()
#     for g in genes:

#         data_gene = depth_df_exonlab_exondepth_genedepth_sample_depth.loc[depth_df_exonlab_exondepth_genedepth_sample_depth.GENE == g]
#         for proj in projects:
#             data = data_gene.loc[data_gene.PROJECT_NAME == proj]

#             fig, (ax1, ax2) = plt.subplots(1, 2)
#             fig.set_size_inches(17, 5)

#             sns.lineplot(data = data, x = "fake_pos", y = "READS", alpha = .7, hue = "SAMPLE_ID",
#                         # hue_order = samples_list,
#                         ax = ax1, legend = False)
#             ax1.axhline(data.DEPTH.mean(), linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#             ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()+2*data.DEPTH.std(), facecolor = '0.8')
#             ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()-2*data.DEPTH.std(), facecolor = '0.8')
#             ax1.set_title(f"{g} positions in {proj}")
#             ax1.set_xticks([])

#             sns.lineplot(data = data, x = "fake_pos", y = "cov_norm_exon_depth", alpha = .7, hue = "SAMPLE_ID",
#                         # hue_order = samples_list,
#                         ax = ax2)
#             ax2.axhline(1, linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#             ax2.set_title(f"{g} positions in {proj} (normalized by exons's mean DEPTH)")
#             ax2.set_xticks([])
#             ax2.legend(fontsize = 8, bbox_to_anchor= (1.25,1))

#             plt.tight_layout()
#             plt.show()
#             pdf.savefig()
#             plt.close()


#     # In[105]:


#     # plot no normalized and mean DEPTH-normalized together
#     genes = depth_df_exonlab_exondepth_genedepth_sample_depth.GENE.unique()
#     projects = depth_df_exonlab_exondepth_genedepth_sample_depth.PROJECT_NAME.unique()
#     for g in genes:

#         data_gene = depth_df_exonlab_exondepth_genedepth_sample_depth.loc[depth_df_exonlab_exondepth_genedepth_sample_depth.GENE == g]
#         for proj in projects:
#             data = data_gene.loc[data_gene.PROJECT_NAME == proj]

#             fig, (ax1, ax2) = plt.subplots(1, 2)
#             fig.set_size_inches(17, 5)

#             sns.lineplot(data = data, x = "fake_pos", y = "READS", alpha = .7, hue = "SAMPLE_ID",
#                         # hue_order = samples_list,
#                         ax = ax1, legend = False)
#             ax1.axhline(data.DEPTH.mean(), linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#             ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()+2*data.DEPTH.std(), facecolor = '0.8')
#             ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()-2*data.DEPTH.std(), facecolor = '0.8')
#             ax1.set_title(f"{g} positions in {proj}")
#             ax1.set_xticks([])

#             sns.lineplot(data = data, x = "fake_pos", y = "cov_norm_gene_depth", alpha = .7, hue = "SAMPLE_ID",
#                         # hue_order = samples_list,
#                         ax = ax2)
#             ax2.axhline(1, linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#             ax2.set_title(f"{g} positions in {proj} (normalized by GENE's mean DEPTH)")
#             ax2.set_xticks([])
#             ax2.legend(fontsize = 8, bbox_to_anchor= (1.25,1))

#             plt.tight_layout()
#             plt.show()
#             pdf.savefig()
#             plt.close()


#     # In[106]:


#     # plot no normalized and mean DEPTH-normalized together
#     genes = depth_df_exonlab_exondepth_genedepth_sample_depth.GENE.unique()
#     projects = depth_df_exonlab_exondepth_genedepth_sample_depth.PROJECT_NAME.unique()
#     for g in genes:

#         data_gene = depth_df_exonlab_exondepth_genedepth_sample_depth.loc[depth_df_exonlab_exondepth_genedepth_sample_depth.GENE == g]
#         for proj in projects:
#             data = data_gene.loc[data_gene.PROJECT_NAME == proj]

#             fig, (ax1, ax2) = plt.subplots(1, 2)
#             fig.set_size_inches(17, 5)

#             sns.lineplot(data = data, x = "fake_pos", y = "READS", alpha = .7, hue = "SAMPLE_ID",
#                         # hue_order = samples_list,
#                         ax = ax1, legend = False)
#             ax1.axhline(data.DEPTH.mean(), linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#             ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()+2*data.DEPTH.std(), facecolor = '0.8')
#             ax1.axhspan(data.DEPTH.mean(), data.DEPTH.mean()-2*data.DEPTH.std(), facecolor = '0.8')
#             ax1.set_title(f"{g} positions in {proj}")
#             ax1.set_xticks([])

#             sns.lineplot(data = data, x = "fake_pos", y = "cov_norm_sample_depth", alpha = .7, hue = "SAMPLE_ID",
#                         # hue_order = samples_list,
#                         ax = ax2)
#             ax2.axhline(1, linestyle = "--", linewidth = 2, color = "black", alpha = 0.5)
#             ax2.set_title(f"{g} positions in {proj} (normalized by SAMPLE_ID's mean DEPTH)")
#             ax2.set_xticks([])
#             ax2.legend(fontsize = 8, bbox_to_anchor= (1.25,1))

#             plt.tight_layout()
#             plt.show()
#             pdf.savefig()
#             plt.close()










# def subset_mutation_dataframe(mutations_file, json_filters):
#     """
#     INFO
#     """
#     # Load your MAF DataFrame (raw_annotated_maf)
#     raw_annotated_maf = pd.read_csv(mutations_file, sep = "\t", header = 0)


#     # Load the filter criteria from the JSON file
#     with open(json_filters, 'r') as file:
#         filter_criteria = json.load(file)

#     if len(filter_criteria) > 0:
#         # Filter the annotated maf using the described filters
#         print("MAF subset")
#         return filter_maf(raw_annotated_maf, filter_criteria)

#     return raw_annotated_maf


# def plot_mutations_per_sample(sample_name, maf, sample_column_name = "SAMPLE_ID"):
#     muts_per_sample = maf.groupby(sample_column_name).size().to_frame("NUM_MUTS").reset_index()

#     # Calculate the length of the longest SAMPLE_ID name
#     max_label_length = muts_per_sample[sample_column_name].str.len().max()

#     # Determine the rotation angle and adjust figure size accordingly
#     rotation_angle = 90 if max_label_length > 7 else 30  # Adjust threshold as needed
#     # fig_height = 4 + (max_label_length / 10) * 2  # Adjust multiplier as needed
#     # # Calculate the figure width based on the number of samples
#     # fig_width = min(18, max(2, len(muts_per_sample) * 0.5))

#     # # Create the figure and axis
#     # fig, ax = plt.subplots(figsize=(fig_width, fig_height))

#     fig, ax = plt.subplots(figsize=(18, 4))
#     sns.barplot(data = muts_per_sample, x = sample_column_name, y = "NUM_MUTS",
#                 ax = ax, palette = ["salmon"], dodge = False)
#     plt.xticks(rotation = rotation_angle)
#     plt.xlabel("")
#     plt.ylabel("Number of mutations", fontsize = 14)
#     ax.set_title(f"Somatic mutations in {sample_name}")
#     plt.tight_layout()

#     return fig


# def plot_mutations_per_gene(sample_name, maf, parameters = {}):
#     default_vals = {'gene_column_name' : 'SYMBOL', 'minimum_muts' : 5, 'annotation' : True}
#     default_vals.update(parameters)

#     muts_per_gene = maf.groupby(default_vals['gene_column_name']).size().to_frame("NUM_MUTS")
#     muts_per_gene = muts_per_gene[muts_per_gene["NUM_MUTS"] >= default_vals['minimum_muts']].reset_index()

#     # # Calculate the length of the longest SAMPLE_ID name
#     # max_label_length = muts_per_gene[default_vals['gene_column_name']].str.len().max()

#     # # Determine the rotation angle and adjust figure size accordingly
#     # rotation_angle = 30 if max_label_length > 10 else 0  # Adjust threshold as needed
#     # fig_height = 4 + (max_label_length / 10) * 2  # Adjust multiplier as needed
#     # # Calculate the figure width based on the number of samples
#     # fig_width = min(18, max(2, len(muts_per_gene) * 0.5))

#     # # Create the figure and axis
#     # fig, ax = plt.subplots(figsize=(fig_width, fig_height))
#     fig, ax = plt.subplots(figsize=(18, 4))

#     muts_per_gene = muts_per_gene.sort_values(by = "NUM_MUTS", ascending = False).reset_index(drop = True)

#     # Plot the barplot
#     sns.barplot(data=muts_per_gene,
#                 x = default_vals['gene_column_name'], y="NUM_MUTS",
#                 ax=ax,
#                 color="salmon"
#                 )
#     if default_vals['annotation']:
#         # Calculate the maximum height of the bars
#         max_height = max(bar.get_height() for bar in ax.patches)

#         # Set the spacing above the bars for annotations
#         annotation_spacing = max_height * 0.05  # Adjust this value as needed

#         # Increase y-axis limits to accommodate annotations
#         y_min, y_max = ax.get_ylim()
#         ax.set_ylim(y_min, y_max + annotation_spacing)

#         # Add annotations above each bar
#         for bar in ax.patches:
#             height = bar.get_height()
#             ax.annotate(format(height, '.0f'),
#                         (bar.get_x() + bar.get_width() / 2, height),
#                         ha='center', va='bottom', fontsize=12, xytext=(0, 3),
#                         textcoords='offset points')

#     ax.set_xlabel("Genes", fontsize=14)
#     ax.set_ylabel("Number mutations", fontsize=14)

#     plt.xticks(fontsize=10, rotation=30)
#     plt.yticks(fontsize=12)
#     plt.tick_params(axis='x', which='major', labelsize=14)
#     ax.set_title(f"Mutations per GENE in {sample_name}")
#     plt.tight_layout()

#     return fig


# def plot_filter_stats(df, x_axis_group, hue_group, logy = False, stacked = True, fig_width = 20, fig_height = 5):
#     """
#     Only admits two variables to combine, ideally
#     the x_axis_group variable should be SAMPLE_ID or
#     GENE symbol
#     """

#     # create groupby object and count instances
#     ## displayed stacked values
#     if stacked:
#         gp_df = df.groupby([x_axis_group, hue_group], dropna = False).size().to_frame('number_mutations').reset_index().pivot(
#         columns=hue_group, index=x_axis_group, values="number_mutations")

#     # displayed in separated bars
#     else:
#         gp_df = df.groupby([x_axis_group, hue_group], dropna = False).size().to_frame('number_mutations').reset_index()
#         gp_df = gp_df.fillna("FALSE")

#     # plot
#     fig, ax = plt.subplots(1, 1)
#     fig.set_size_inches(fig_width, fig_height)
#     plt.xticks(rotation = 90)

#     if stacked:
#         gp_df.plot(kind = 'bar', stacked = True, ax = ax)
#     else:
#         sns.barplot(data = gp_df, x = x_axis_group, y = "number_mutations", hue = hue_group, ax = ax)
#         if logy:
#             ax.set_yscale('log')

#     return None



# def filter_wrapper(sample_name, maf, parameters = {}):
#     var_to_plot = parameters.get('variable', "SAMPLE_ID")
#     filters_to_plot = parameters.get('filter_list', ['not_in_panel', 'no_pileup_support', 'n_rich'])

#     for filt in filters_to_plot:
#         plot_filter_stats(maf, var_to_plot, filt)
#         pdf.savefig()
#         plt.close()






# dict_plotname2func = {
#     "per_gene" : plot_mutations_per_gene,
#     "per_sample" : plot_mutations_per_sample,
#     "filter_stats" : filter_wrapper
# }

# def plot_manager(sample_name, maf, plotting_criteria_file):
#     # Load the filter criteria from the JSON file
#     with open(plotting_criteria_file, 'r') as file:
#         plotting_criteria = json.load(file)

#     for suffix, criteria_list in plotting_criteria.items():
#         suffix = suffix.strip('.')

#         with PdfPages(f'{sample_name}.{suffix}.pdf') as pdf:
#             for criterion in criteria_list:
#                 if criterion.startswith("filter_stats"):
#                     main_criterion = criterion.split(' ')[0]
#                     gene_or_sample = criterion.split(' ')[1]
#                     filters = criterion.split(' ')[2].split(",")
#                     fig1 = dict_plotname2func[main_criterion](sample_name, maf, parameters = {'variable' : gene_or_sample, 'filter_list' : filters})


#                 elif ' ' in criterion:
#                     main_criterion = criterion.split(' ')[0]
#                     additional_params = { x.split(":")[0] : x.split(":")[1] for x in criterion.split(' ')[1:] }
#                     fig1 = dict_plotname2func[main_criterion](sample_name, maf, parameters = additional_params)
#                     pdf.savefig()
#                     plt.close()


#                 else :
#                     fig1 = dict_plotname2func[criterion](sample_name, maf)
#                     pdf.savefig()
#                     plt.close()



# # @click.command()
# # @click.option('--sample_name', type=str, help='Name of the SAMPLE_ID being processed.')
# # @click.option('--mut_file', type=click.Path(exists=True), help='Input mutation file')
# # @click.option('--out_maf', type=click.Path(), help='Output MAF file')
# # @click.option('--json_filters', type=click.Path(exists=True), help='Input mutation filtering criteria file')
# # @click.option('--req_plots', type=click.Path(exists=True), help='Column names to output')
# # # @click.option('--plot', is_flag=True, help='Generate plot and save as PDF')

# # def main(sample_name, mut_file, out_maf, json_filters, req_plots): # , plot):
# #     click.echo(f"Subsetting MAF file...")
# #     subset_mutation_dataframe(sample_name, mut_file, out_maf, json_filters, req_plots)

# # if __name__ == '__main__':
# #     main()


# sample_name  = sys.argv[1]
# mut_file     = sys.argv[2]
# out_maf      = sys.argv[3]
# json_filters = sys.argv[4]
# req_plots    = sys.argv[5]



# if __name__ == '__main__':
#     maf = subset_mutation_dataframe(mut_file, json_filters)
#     plot_manager(sample_name, maf, req_plots)
