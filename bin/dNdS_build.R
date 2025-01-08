#!/opt/conda/bin/Rscript --vanilla

# Load required library
library(optparse)
library(dndscv)

# Define command-line options
option_list <- list(
    make_option(c("-b", "--biomart_cds"), type="character", help="Path to Biomart CDS file", metavar="character"),
    make_option(c("-r", "--reference_genome"), type="character", help="Path to reference genome file", metavar="character"),
    make_option(c("-o", "--outfile"), type="character", help="Output file for RefCDS", metavar="character")
)

# Parse the command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Call buildref function with provided arguments
buildref(opt$biomart_cds, opt$reference_genome, outfile = opt$outfile)
