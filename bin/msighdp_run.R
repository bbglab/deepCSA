#!/usr/bin/Rscript --vanilla


# This is designed to run in a singularity/apptainer container.
# Do not call it directly from the command line.
#
# Suggested call for simple example
# singularity exec <my_R_container> Rscript --vanilla test_mSigHdp.R 123
#
# But normally something like:
# nice singularity exec <my_R_container> Rscript --vanilla test_mSigHdp.R 123 >& log.txt &

# We need to tell R to look for packages first inside the container
# but also to be able to look outside the container (the "host"
# system). This works only if the host OS is the same as
# the container OS, or if the libraries in the host system do not
# have compiled code. In the last case, you would have to
# convert the container to a sandbox, add the necessary libraries,
# and then convert back to a .sif.
.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))

# args have the arguments passed to Rscript from bash script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Please provide the random seed as an argument")
}
seed <- as.numeric(args[1])
input_file <- args[2]
output_home <- args[3]
K_guess <- as.numeric(args[4])
num_cpus <-  as.numeric(args[5])

sigprofiler <- TRUE


library(ICAMS)
library(mSigHdp)


if (!dir.exists(output_home)) {
    if (!dir.create(path = output_home, recursive = TRUE)) {
        stop("unable to create ", output_home)
    }
}

checkpoint_dir <- file.path(output_home, "checkpoint")
if (!dir.exists(checkpoint_dir)) {
    if (!dir.create(path = checkpoint_dir, recursive = TRUE)) {
        stop("unable to create ", checkpoint_dir)
    }
}


# Load data
data <- read.table(input_file, header = TRUE, sep = "\t")

# Check if Mutation_type and Trinucleotide columns exist
if (!("Mutation type" %in% colnames(data)) && !("Trinucleotide" %in% colnames(data))) {
    if ("CONTEXT_MUT" %in% colnames(data)) {
        if (sigprofiler) {
            # Define sample columns
            sample_columns <- colnames(data)[colnames(data) != "CONTEXT_MUT"]

            # Create new columns
            data$`Mutation type` <- substr(data$CONTEXT_MUT, 3, 5)
            data$Trinucleotide <- paste(substr(data$CONTEXT_MUT, 1, 1), substr(data$CONTEXT_MUT, 3, 3), substr(data$CONTEXT_MUT, 7, 7), sep = '')

            # Reorder columns
            data <- data[,c("Mutation type", "Trinucleotide", sample_columns)]
            data <- data[order(data[, 1], data[, 2]), ]
            print(head(data))

            # Write to output file
            reformatted_file <- file.path(output_home, "reformatted_input.tsv")
            write.table(data, file = reformatted_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

        } else {
            # Define sample columns
            sample_columns <- colnames(data)[colnames(data) != "CONTEXT_MUT"]

            # Create new columns
            data$`Mutation type` <- paste(substr(data$CONTEXT_MUT, 2, 2), ">", substr(data$CONTEXT_MUT, 4, 4), sep = '')
            data$Trinucleotide <- substr(data$CONTEXT_MUT, 1, 3)

            # Reorder columns
            data <- data[,c("Mutation type", "Trinucleotide", sample_columns)]
            data <- data[order(data[, 1], data[, 2]), ]

            # Write to output file
            reformatted_file <- file.path(output_home, "reformatted_input.tsv")
            write.table(data, file = reformatted_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
    } else {
        stop("Column 'CONTEXT_MUT' not found in the input data.")
    }
} else {
    # If Mutation_type and Trinucleotide columns exist, keep the original data
    reformatted_file <- file.path(output_home, "reformatted_input.tsv")
    write.table(data, file = reformatted_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}




# All inputs and outputs are in the current working directory
input_catalog <- ICAMS::ReadCatalog(file = reformatted_file)

print(input_catalog)

# Set working directory to checkpoint_dir directory so that the checkpoint Rdata
# files will be saved there
setwd(checkpoint_dir)

message("Start running mSigHdp")
message("Using seed ", seed)
message("Using K_guess ", K_guess)
message("mSigHdp version: ")
print(packageVersion("mSigHdp"))

# The parameters used here are for testing purpose only to save time.
# For recommended parameters to use, please refer to this paper
# Mo Liu, Yang Wu, Nanhai Jiang, Arnoud Boot, Steven G. Rozen,
# mSigHdp: hierarchical Dirichlet process mixture modeling
# for mutational signature discovery,
# https://doi.org/10.1093/nargab/lqad005

# TODO
# revise whether we should change any of the parameters here for a better optimization
mSigHdp::RunHdpxParallel(
    input.catalog           = input_catalog,
    seedNumber              = seed,
    K.guess                 = K_guess,
    out.dir                 = output_home,
    multi.types             = FALSE,
    burnin                  = 5000,
    burnin.multiplier       = 20,
    post.n                  = 200,
    post.space              = 100,
    num.child.process       = max(20, num_cpus),
    CPU.cores               = num_cpus,
    high.confidence.prop    = 0.3,
    gamma.alpha             = 1,
    gamma.beta              = 20,
    checkpoint              = TRUE,
    verbose                 = TRUE,
    downsample_threshold    = 3000
)

message("Finished running mSigHdp")
