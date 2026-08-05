#!/opt/conda/bin/Rscript --vanilla


library(Hmisc)
library(tidyr)
library(stringr)
library(dplyr, warn = FALSE)
library(ggplot2)
library(jsonlite)
library(Biostrings)
library(optparse)
library(R.utils)
library(data.table)

## Read it from a TSV and format it as json

## command line arguments
option_list = list(
  make_option(c("-n", "--samplename"), type="character", default=NULL,
              help="sample name/identifier of the run", metavar="character"),
  make_option(c("-m", "--mutations"), type="character", default=NULL,
              help="mutation dataset file name", metavar="character"),
  make_option(c("-o", "--outputprefix"), type="character", default="dNdScv_output",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--depths"), type="character", default=NULL,
              help="depths dataset file name", metavar="character"),
  make_option(c("-c", "--consensus_bed"), type="character", default=NULL,
              help="consensus bed file name", metavar="character"),
  make_option(c("-w", "--wgs_counts"), type="character", default=NULL,
              help="WGS trinucleotide counts file name", metavar="character"),
  make_option(c("-p", "--panel_version"), type="character", default=NULL,
              help="panel version", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


sample_name <- opt$samplename
mutations <- opt$mutations
depths <- opt$depths
consensus_bed <- opt$consensus_bed
wgs_counts <- opt$wgs_counts
output_name <- opt$outputprefix
panel_version <- opt$panel_version

print(paste("Sample name: ", sample_name))
print(paste("Mutations file: ", mutations))
print(paste("Depths file: ", depths))
print(paste("Consensus bed file: ", consensus_bed))
print(paste("WGS counts file: ", wgs_counts))
print(paste("Output name: ", output_name))
print(paste("Panel version: ", panel_version))

path2out = "wgs_mutrate"

get_genome_content <- function(genome_content){
  #' Get genome trinucleotide content
  #'
  #' This function reads trinucleotide counts from a JSON genome composition file
  #' and computes the 96 pyrimidine-centered contexts. 

  #' As input uses path to the json file with genome trinucleotide content

  counts_df <- read.table(genome_content, header=TRUE, sep="\t")
  colnames(counts_df) <- c("CONTEXT", "N_sites_genome")
  message("Contexts in genome sites = ", nrow(counts_df))
  return(counts_df)
}  


get_consensus_sites_depth <- function(sample, depths_file, consensus_bed){
  #' Get depth per position for positions in consensus panel
  #'
  #' This function intersects consensus panel with file with annotated depth per position

  # Load depth data
  dt_pos <- fread(depths_file)
  colnames(dt_pos) <- c("CHROM", "POS", "CONTEXT", "DEPTH")
  dt_pos[, `:=`(
    start = POS,
    end   = POS
  )]
  setkey(dt_pos, CHROM, start, end)
  # Load consensus bed
  dt_bed <- fread(consensus_bed, col.names = c("CHROM", "start", "end"))
  setkey(dt_bed, CHROM, start, end)
  # Overlap
  hits <- foverlaps(dt_pos, dt_bed, type = "any", nomatch = 0L)
  consensus_depth <- hits[, .(CHROM, POS, CONTEXT, DEPTH)]
  return(consensus_depth)
}

get_mutations_and_sites <- function(sample_name, sites_file, mutations_file, consensus_bed, genome_sites_df){
  #' Get number of mutations per sample and normalize panel content to whole genome content
  #'
  #' This function gets number of mutations per sample in each context and normalizes panel content to whole genome content 
  result <- NULL
  print(sample_name)

  df_sites = get_consensus_sites_depth(sample_name, sites_file, consensus_bed)
  print("1")
  df_sites$depth = as.numeric(df_sites$DEPTH)
  print("2")
  df_sites_panel_agg = df_sites %>%
    group_by(CONTEXT) %>% 
    summarise(N = sum(DEPTH))
  print("3")
  df_sites_panel_agg = as.data.frame(df_sites_panel_agg)
  colnames(df_sites_panel_agg) <- c("CONTEXT", "N_sites_panel")
  message("Contexts in panel sites = ", nrow(df_sites_panel_agg))
  print(head(df_sites_panel_agg))

  df_sites = merge(genome_sites_df, df_sites_panel_agg, by="CONTEXT")
  print(head(df_sites))
  message("Contexts in panel and genome sites = ", nrow(df_sites))
  df_sites <- df_sites %>% mutate(proportion_genome=N_sites_genome/sum(N_sites_genome))
  df_sites <- df_sites %>% mutate(proportion_panel=N_sites_panel/sum(N_sites_panel))
  df_sites$ratio2genome = df_sites$proportion_panel/df_sites$proportion_genome
  print(head(df_sites))
  print(sum(df_sites$N_sites_panel))

  df_mutations = read.table(mutations_file, header=TRUE, sep="\t")
  df_mutations = df_mutations[df_mutations$TYPE=="SNV",]
  df_mutations$CONTEXT = str_split_i(df_mutations$CONTEXT_MUT, ">", 1)
  df_mutations_agg = df_mutations %>%
    group_by(CONTEXT) %>% 
    summarise(N_mut = n())
  print(head(df_mutations_agg))
  df_mutations = as.data.frame(df_mutations)
  result_sample = merge(df_sites, df_mutations_agg, by="CONTEXT", all=TRUE)
  # if some contexts are absent in mutataions - keep them but put mutation number to 0
  if (nrow(result_sample[is.na(result_sample$N_mut),]) > 0){
    result_sample[is.na(result_sample$N_mut),]$N_mut <- 0
  }
  print(head(result_sample))  
  message("Contexts in df with mutations = ", nrow(result_sample))
  result_sample$N_mut_corrected = result_sample$N_mut * result_sample$ratio2genome
  sample_out <- c(sample_name, sum(result_sample$N_mut), sum(result_sample$N_mut_corrected), sum(df_sites$N_sites_panel))
  result <- rbind(result, sample_out) 

  return(result)
}

#Download df with genome trinucleotide contexts
df_sites_genome_agg <- get_genome_content(wgs_counts)

#Get number of mutations and normalize panel contetnt to genome content
result <- get_mutations_and_sites(sample_name, depths, mutations, consensus_bed, df_sites_genome_agg)
print(head(result))


result <- as.data.frame(result)
colnames(result) <- c("sample", "N_mut", "N_mut_corrected", "DEPTH")
result$N_mut_corrected = as.numeric(result$N_mut_corrected)
result$DEPTH = as.numeric(result$DEPTH)
result$N_mut = as.numeric(result$N_mut)
result$mutrate_observed = result$N_mut_corrected/result$DEPTH
result$mutrate_observed_per_MB <- result$mutrate_observed * 10**6

# why the indices of the two boundaries are different? is this desired or a typo?
result$mutrate_CI_high <- apply(result, 1, function(x) binconf(as.numeric(x["N_mut_corrected"]), as.numeric(x["DEPTH"]), alpha=0.05, method=c("wilson","exact","asymptotic","all"), include.x=FALSE, include.n=FALSE, return.df=FALSE)[3])
result$mutrate_CI_low <- apply(result, 1, function(x) binconf(as.numeric(x["N_mut_corrected"]), as.numeric(x["DEPTH"]), alpha=0.05, method=c("wilson","exact","asymptotic","all"), include.x=FALSE, include.n=FALSE, return.df=FALSE)[2])
result$Muts_per_cell <- result$mutrate_observed*2*sum(df_sites_genome_agg$N_sites_genome)
result$panel_version <- panel_version
print((result))

message("Output will be written to ", paste0(output_name, "_mutrates_results.tsv"))
write.table(result, paste0(output_name, "_mutrates_results.tsv"), sep="\t", row.names = FALSE, quote = FALSE)

# setwd(path2out)

# jpeg(filename=paste("mutrate_trint_corrected.with_nanoseq.jpeg", sep=""), width=30, height=15, res=300, units='cm')

# ggplot(result, aes(x=sample, y=mutrate_observed)) +
#   geom_bar(stat="identity", position="dodge", fill="grey") +
#   geom_errorbar(aes(x=sample, ymin=mutrate_CI_low, ymax=mutrate_CI_high), width=0.4,  alpha=0.9, linewidth=1, position=position_dodge(.9)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust=1)) +
#   geom_text(aes(label = N_mut, x = sample, y = mutrate_observed), position = position_dodge(width = 0.9), vjust = -0.5, hjust = -0.1) +
#   xlab("") +
#   ylab("Mutation rate") +
#   theme(legend.position="bottom")
# dev.off()

# jpeg(filename=paste("mutrate_trint_corrected.jpeg", sep=""), width=30, height=15, res=300, units='cm')
# result<-result[result$protocol != "Nanoseq_Sanger",]
# ggplot(result, aes(x=sample, y=mutrate_observed)) +
#   geom_bar(stat="identity", position="dodge", fill="grey") +
#   geom_errorbar(aes(x=sample, ymin=mutrate_CI_low, ymax=mutrate_CI_high), width=0.4,  alpha=0.9, linewidth=1, position=position_dodge(.9)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust=1)) +
#   geom_text(aes(label = N_mut, x = sample, y = mutrate_observed), position = position_dodge(width = 0.9), vjust = -0.5, hjust = -0.1) +
#   xlab("") +
#   ylab("Mutation rate") +
#   theme(legend.position="bottom")
# dev.off()
