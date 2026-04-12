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
    make_option(c("-g", "--genelist"), type="character", default=NULL,
                help="File containing comma-separated gene names", metavar="character"),
    make_option(c("-o", "--outputfile"), type="character", default=NULL,
                help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read genes
genes <- readLines(opt$genelist)
genes <- unlist(strsplit(genes, ","))

# Define the BioMart endpoint (archived Ensembl January 2024 version)
biomart_url <- "http://jan2024.archive.ensembl.org/biomart/martservice"

# Function to build query for a batch of genes
build_query <- function(gene_batch) {
  gene_string <- paste0(gene_batch, collapse = ",")
  query <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<!DOCTYPE Query>',
    '<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >',
    '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >',
    '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>',
    '<Filter name = "biotype" value = "protein_coding"/>',
    '<Filter name = "external_gene_name" value = "', gene_string, '"/>',
    '<Filter name = "mane_select" excluded = "0"/>',
    '<Attribute name = "ensembl_gene_id" />',
    '<Attribute name = "external_gene_name" />',
    '<Attribute name = "ensembl_peptide_id" />',
    '<Attribute name = "chromosome_name" />',
    '<Attribute name = "genomic_coding_start" />',
    '<Attribute name = "genomic_coding_end" />',
    '<Attribute name = "cds_start" />',
    '<Attribute name = "cds_end" />',
    '<Attribute name = "cds_length" />',
    '<Attribute name = "strand" />',
    '<Attribute name = "ensembl_transcript_id" />',
    '<Attribute name = "exon_chrom_start" />',
    '<Attribute name = "exon_chrom_end" />',
    '</Dataset>',
    '</Query>'
  )
  return(query)
}

# Initialize output file (remove if it exists to ensure fresh start)
if (file.exists(opt$outputfile)) {
  file.remove(opt$outputfile)
}

# Process genes in batches of 300
batch_size <- 300
batches <- split(genes, ceiling(seq_along(genes)/batch_size))
has_written_data <- FALSE

for (i in seq_along(batches)) {
  message(sprintf("Processing batch %d of %d", i, length(batches)))
  batch_query <- build_query(batches[[i]])
  encoded_query <- URLencode(batch_query)
  
  response <- GET(paste0(biomart_url, "?query=", encoded_query))
  
  if (response$status_code != 200) {
    stop("Error: Failed to retrieve data from BioMart for batch ", i, ". Status code: ", response$status_code)
  }
  
  biomart_output <- content(response, type = "text", encoding = "UTF-8")
  
  # Check if response is empty (BioMart sometimes returns empty responses for small valid queries)
  if (nchar(trimws(biomart_output)) > 0) {
    batch_data <- read.delim(textConnection(biomart_output), header = TRUE, sep = "\t")
    
    # Remove rows with empty CDS information
    exon_file <- subset(batch_data, !is.na(Genomic.coding.start) & Genomic.coding.start != "")
    
    if (nrow(exon_file) > 0) {
      if (!has_written_data) {
        print(head(exon_file))
        has_written_data <- TRUE
      }
      
      # Save the updated exon data to file iteratively using append
      write.table(
        exon_file,
        file = opt$outputfile,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE
      )
    }
  }
}

if (!has_written_data) {
  stop("No data retrieved from BioMart for any batch.")
}