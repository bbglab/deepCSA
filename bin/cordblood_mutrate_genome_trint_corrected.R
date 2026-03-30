library(Hmisc)
library(tidyr)
library(stringr)
library(dplyr, warn = FALSE)
library(ggplot2)
library(jsonlite)
library(Biostrings)
library(data.table)

genome_content_json <- "../data/genome_counts_tribases.json"

deepCSA_run <- "2025-12-15"
add_nanoseq_for_comparison = TRUE #if you want to add cord blood sequenced with Nanoseq for comparison

# root_dir <- "../data/cord_blood_run"
root_dir <- "/data/bbg/nobackup2/prominent/duplex_seq_tests/error_rate/cord_blood/bbg/deepCSA/2026-03-26_deepUMIcaller_and_dupcaller_4_paper"
deepCSA_run_dir <- root_dir # paste0(root_dir,"deepCSA/", deepCSA_run)

path2sites = paste0( deepCSA_run_dir, "/depths/individual/")
path2mutations = paste0( deepCSA_run_dir, "/mutations/clean_somatic/")
consensus_bed <- paste0( deepCSA_run_dir, "/regions/consensuspanels/consensus.all.bed")

path2out = "results/cordblood_mutrate"

get_genome_content <- function(genome_content_json){
  #' Get genome trinucleotide content
  #'
  #' This function reads trinucleotide counts from a JSON genome composition file
  #' and computes the 96 pyrimidine-centered contexts. 

  #' As input uses path to the json file with genome trinucleotide content

  genome_json <- read_json(genome_content_json)
  df_sites_genome <- data.frame(Context = names(genome_json), sites = unlist(genome_json), stringsAsFactors = FALSE)
  df_sites_genome <- df_sites_genome %>%
    filter(!grepl("N",Context))
  df_sites_genome$Compl_context <- as.character(reverseComplement(DNAStringSet(df_sites_genome$Context)))
  df_sites_genome$CONTEXT <- ifelse(substr(df_sites_genome$Context,2,2) %in% c("T","C"),
                       df_sites_genome$Context,
                       df_sites_genome$Compl_context)
  df_sites_genome_agg = df_sites_genome %>%
    group_by(CONTEXT) %>% 
    summarise(N_sites_genome= sum(sites))
  df_sites_genome_agg <- as.data.frame(df_sites_genome_agg)  
  message("Contexts in genome sites = ", nrow(df_sites_genome_agg))
  return(df_sites_genome_agg)
}  


get_consensus_sites_depth <- function(sample, depth_path, consensus_bed){
  #' Get depth per position for positions in consensus panel
  #'
  #' This function intersects consensus panel with file with annotated depth per position

  # Load depth data
  dt_pos <- fread(paste0(depth_path, sample, ".depths.annotated.tsv.gz"))
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

get_mutations_and_sites <- function(path2sites, path2mutations, consensus_bed, genome_sites_df){
  #' Get number of mutations per sample and normalize panel content to whole genome content
  #'
  #' This function gets number of mutations per sample in each context and normalizes panel content to whole genome content 
  result <- NULL
  sites_files <- list.files(path=path2sites, pattern=glob2rx("*.depths.annotated.tsv.gz"))
  sites_files <- sites_files[sites_files != 'all_samples.depths.annotated.tsv.gz']
  for(file in sites_files){
    sample_name = str_split_i(file, ".depths", 1)
    print(sample_name)

    df_sites = get_consensus_sites_depth(sample_name, path2sites, consensus_bed)
    df_sites$depth = as.numeric(df_sites$DEPTH)
    df_sites_panel_agg = df_sites %>%
      group_by(CONTEXT) %>% 
      summarise(N = sum(DEPTH))
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

    df_mutations = read.table(paste(path2mutations, sample_name, ".somatic.mutations.tsv", sep=""), header=TRUE, sep="\t")
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
    sample_out <- c(sample_name, sum(result_sample$N_mut), sum(result_sample$N_mut_corrected), sum(df_sites$N_sites_panel), "sample", sample_name)
    result <- rbind(result, sample_out) 

  }
  return(result)
}

#Download df with genome trinucleotide contexts
df_sites_genome_agg <- get_genome_content(genome_content_json)

#Get number of mutations and normalize panel contetnt to genome content
result <- get_mutations_and_sites(path2sites, path2mutations, consensus_bed, df_sites_genome_agg)
print(head(result))

# Add nanoseq rate for comparison if needed
if (add_nanoseq_for_comparison == TRUE){
  result <- rbind(result, c("PD48442_cordblood_nanoseqv2", 41, 38.43131656, 2799554062, "Nanoseq_Sanger", "PD48442"))
  result <- rbind(result, c("PD47269_cordblood_nanoseqv2", 26, 28.29844724, 2003725667, "Nanoseq_Sanger", "PD47269"))
}
result <- as.data.frame(result)
colnames(result) <- c("sample", "N_mut", "N_mut_corrected", "DEPTH", "protocol", "donor_id")
result$N_mut_corrected = as.numeric(result$N_mut_corrected)
result$DEPTH = as.numeric(result$DEPTH)
result$N_mut = as.numeric(result$N_mut)
result$mutrate_observed = result$N_mut_corrected/result$DEPTH
result$mutrate_observed_per_MB <- result$mutrate_observed * 10**6

# why the indices of the two boundaries are different? is this desired or a typo?
result$mutrate_CI_high <- apply(result, 1, function(x) binconf(as.numeric(x["N_mut_corrected"]), as.numeric(x["DEPTH"]), alpha=0.05, method=c("wilson","exact","asymptotic","all"), include.x=FALSE, include.n=FALSE, return.df=FALSE)[3])
result$mutrate_CI_low <- apply(result, 1, function(x) binconf(as.numeric(x["N_mut_corrected"]), as.numeric(x["DEPTH"]), alpha=0.05, method=c("wilson","exact","asymptotic","all"), include.x=FALSE, include.n=FALSE, return.df=FALSE)[2])
result$Muts_per_cell <- result$mutrate_observed*2*sum(df_sites_genome_agg$N_sites_genome)
print((result))

message("Output will be written to ", path2out)
if (!dir.exists(path2out)) {
  dir.create(path2out, recursive = TRUE)
}

write.csv(result, paste0(path2out,"/mutrates_results.tsv"),row.names = FALSE, quote = FALSE)

setwd(path2out)

jpeg(filename=paste("mutrate_trint_corrected.with_nanoseq.jpeg", sep=""), width=30, height=15, res=300, units='cm')

ggplot(result, aes(x=sample, y=mutrate_observed)) +
  geom_bar(stat="identity", position="dodge", fill="grey") +
  geom_errorbar(aes(x=sample, ymin=mutrate_CI_low, ymax=mutrate_CI_high), width=0.4,  alpha=0.9, linewidth=1, position=position_dodge(.9)) +
  theme_bw() +
  facet_grid(~factor(protocol, levels=c("tests","IDT","TWS","Nanoseq_Sanger")), scales="free_x", space="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  geom_text(aes(label = N_mut, x = sample, y = mutrate_observed), position = position_dodge(width = 0.9), vjust = -0.5, hjust = -0.1) +
  xlab("") +
  ylab("Mutation rate") +
  theme(legend.position="bottom")
dev.off()

jpeg(filename=paste("mutrate_trint_corrected.jpeg", sep=""), width=30, height=15, res=300, units='cm')
result<-result[result$protocol != "Nanoseq_Sanger",]
ggplot(result, aes(x=sample, y=mutrate_observed)) +
  geom_bar(stat="identity", position="dodge", fill="grey") +
  geom_errorbar(aes(x=sample, ymin=mutrate_CI_low, ymax=mutrate_CI_high), width=0.4,  alpha=0.9, linewidth=1, position=position_dodge(.9)) +
  theme_bw() +
  facet_grid(~factor(protocol, levels=c("tests","IDT","TWS","Nanoseq_Sanger")), scales="free_x", space="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  geom_text(aes(label = N_mut, x = sample, y = mutrate_observed), position = position_dodge(width = 0.9), vjust = -0.5, hjust = -0.1) +
  xlab("") +
  ylab("Mutation rate") +
  theme(legend.position="bottom")
dev.off()