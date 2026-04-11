#!/opt/conda/bin/Rscript --vanilla

## USAGE
# Rscript dNdScv_single_sample.R --inputfile ../../0initial_processing/data/PILOT5/custom_files/tws/all_muts/all_below035_4dNdScv.txt --outputfile ../results/all_below035.tsv --samplename all --genelist genes.txt --genedepth genes_coverage.txt

# Load required libraries
library(optparse)
library(GenomicRanges)
library(dplyr)
library(httr)
library(utils)


option_list = list(
    make_option(c("-b", "--biomartquery"), type="character", default=NULL,
                help="EnsemblBioMart query file", metavar="character"),
    make_option(c("-o", "--outputfile"), type="character", default=NULL,
                help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Define the BioMart endpoint (archived Ensembl January 2024 version)
biomart_url <- "http://jan2024.archive.ensembl.org/biomart/martservice"

# Load and encode the query from 'biomartQuery.txt'
biomart_query <- paste0(readLines(opt$biomartquery), collapse = "")
encoded_query <- URLencode(biomart_query)

# Make the request to BioMart
response <- GET(paste0(biomart_url, "?query=", encoded_query))

# Check for success
if (response$status_code != 200) {
  stop("Error: Failed to retrieve data from BioMart. Status code: ", response$status_code)
}

# Save and filter the response
biomart_output <- content(response, type = "text", encoding = "UTF-8")
biomart_data <- read.delim(textConnection(biomart_output), header = TRUE, sep = "\t")

# Remove rows with empty CDS information
exon_file <- subset(biomart_data, !is.na(Genomic.coding.start) & Genomic.coding.start != "")
print(head(exon_file))

# Save the updated exon data
write.table(
  exon_file,
  file = opt$outputfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)