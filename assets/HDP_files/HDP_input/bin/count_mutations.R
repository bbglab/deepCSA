library(stringi)
library(stringr)
library(jsonlite)  # For reading JSON

args = commandArgs(trailingOnly=TRUE)
profile_WGS_path <- args[1]
mutations_path <- args[2]
output_file <- args[3]
json_file <- args[4]  # Path to the JSON file

# Read and parse the JSON file
json_data <- fromJSON(json_file)
samples <- names(json_data)  # Extract dictionary keys as sample names

# Initialize lists
Nmut_list <- NULL
C2A_list <- NULL
C2T_list <- NULL
C2G_list <- NULL
T2A_list <- NULL
T2C_list <- NULL
T2G_list <- NULL

# Iterate over the sample names
for (sample in samples) {
  # Construct file paths based on the sample name
  file <- file.path(profile_WGS_path, paste0(sample, ".all.profile.tsv.matrix.WGS"))
  
  if (file.exists(file)) {
    print(sample)
    sample_counts <- read.table(file, header=TRUE)
    sample_counts$mutation <- paste(substr(sample_counts$CONTEXT_MUT, 2, 2), substr(sample_counts$CONTEXT_MUT, 5, 5), sep="")
    
    # Calculate counts
    Nmut <- sum(sample_counts[, 2])
    C2A <- sum(sample_counts[sample_counts$mutation == "CA", ][, 2])
    C2T <- sum(sample_counts[sample_counts$mutation == "CT", ][, 2])
    C2G <- sum(sample_counts[sample_counts$mutation == "CG", ][, 2])
    T2A <- sum(sample_counts[sample_counts$mutation == "TA", ][, 2])
    T2C <- sum(sample_counts[sample_counts$mutation == "TC", ][, 2])
    T2G <- sum(sample_counts[sample_counts$mutation == "TG", ][, 2])
    
    # Append to lists
    Nmut_list <- c(Nmut_list, Nmut)
    C2A_list <- c(C2A_list, C2A)
    C2T_list <- c(C2T_list, C2T)
    C2G_list <- c(C2G_list, C2G)
    T2A_list <- c(T2A_list, T2A)
    T2C_list <- c(T2C_list, T2C)
    T2G_list <- c(T2G_list, T2G)
  } else {
    warning(paste("File not found for sample:", sample))
  }
}

# Save results or further processing
# You can write Nmut_list, C2A_list, etc., to the output_file if required.

result = data.frame(SAMPLE_ID = samples, Total_mut_genome = Nmut_list, CtoA_genome = C2A_list, CtoT_genome = C2T_list, CtoG_genome = C2G_list, TtoA_genome = T2A_list, TtoC_genome = T2C_list, TtoG_genome = T2G_list)
write.csv(result,output_file, row.names=FALSE)
print(head(result))


#id_chanels_df <- read.table("/home/mandrianova/Prominent/Bladder/scripts/HDP_sigExtraction/data/ID_channels.txt", header=TRUE)
#print(head(id_chanels_df))
#id_chanels <- id_chanels_df$sigProfiler

#samples <- NULL
#result_df <- NULL
#for (file_name in list.files(path=mutations_path, pattern="*_01.somatic.mutations.tsv", full.names = T)){
    #print(file_name)
    #input_df = read.csv(file_name, sep="\t", header=T)
    #sample <- str_split_i(str_split_i(file_name, "/",8), "[[.]]", 1)
    #input_df$Sample=sample
    #input_df = input_df[input_df$TYPE == "DELETION"|input_df$TYPE == "INSERTION",]
    #input = input_df[,c("Sample","CHROM","POS","REF","ALT")]
    #colnames(input) <- c('Sample', 'chr', 'pos',  'REF', 'ALT')
    #input$pos <- as.numeric(as.character(input$pos))
    #result_df <- rbind(result_df, c(sample, nrow(input)))
#}

#result_df=as.data.frame(result_df)
#colnames(result_df)<-c("SAMPLE_ID", "N_indels")
#result_df$N_indels = as.numeric(result_df$N_indels)
#print(head(result_df))

#final_result = merge(result, result_df, by="SAMPLE_ID")
#print(head(final_result))

#write.csv(final_result,"/workspace/nobackup2/prominent/bladder/signatures/2024-07-29_deepCSA/Mutation_stats_per_genome_per_sample_with_indels.med.all.csv", row.names=FALSE)


