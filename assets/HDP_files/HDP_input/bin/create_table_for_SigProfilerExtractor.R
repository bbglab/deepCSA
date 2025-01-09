library(stringi)
library(stringr)
library(jsonlite)  # For reading JSON
library(utils)

args = commandArgs(trailingOnly=TRUE)
profile_WGS_path <- args[1]
output_file <- args[2]
json_file <- args[3]  # Path to the JSON file

# Read and parse the JSON file
json_data <- fromJSON(json_file)
samples <- names(json_data)  # Extract dictionary keys as sample names

result_counts <- NULL

for (sample in samples) {
  # Construct file path based on the sample name
  file <- file.path(profile_WGS_path, paste0(sample, ".all.profile.tsv.matrix.WGS"))
  
  if (file.exists(file)) {
    print(sample)
    sample_counts <- read.table(file, header=TRUE)
    colnames(sample_counts)[2] <- sample
    sample_counts[, 2] <- round(sample_counts[, 2])
    
    if (is.null(result_counts)) {
      result_counts <- sample_counts
    } else {
      result_counts <- merge(result_counts, sample_counts, by="CONTEXT_MUT")
    }
  } else {
    warning(paste("File not found for sample:", sample))
  }
}

# Add mutation types to the result
result_counts$Mutation.Types <- paste(
  substr(result_counts$CONTEXT_MUT, 1, 1), "[",
  substr(result_counts$CONTEXT_MUT, 2, 2), ">",
  substr(result_counts$CONTEXT_MUT, 5, 5), "]",
  substr(result_counts$CONTEXT_MUT, 3, 3), sep=""
)

print(colnames(result_counts))

# Reorder columns to place Mutation.Types first
result_counts <- result_counts[, c(ncol(result_counts), 2:(ncol(result_counts) - 1))]

# Write the result to the output file
write.table(result_counts, output_file, row.names=FALSE, sep="\t")

