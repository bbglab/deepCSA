#!/opt/conda/bin/Rscript --vanilla

## USAGE
# Rscript dNdScv_single_sample.R --inputfile ../../0initial_processing/data/PILOT5/custom_files/tws/all_muts/all_below035_4dNdScv.txt --outputfile ../results/all_below035.tsv --samplename all --genelist genes.txt --genedepth genes_coverage.txt

# Load required libraries
library(optparse)
library(GenomicRanges)
library(dplyr)
library(biomaRt)
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

# Define the BioMart ensembl connection (archived Ensembl January 2024 version)
ensembl <- useEnsembl(
  biomart = "genes", 
  dataset = "hsapiens_gene_ensembl", 
  host = "https://www.ensembl.org"
)

# Define query attributes
attributes <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "ensembl_peptide_id",
  "chromosome_name",
  "genomic_coding_start",
  "genomic_coding_end",
  "cds_start",
  "cds_end",
  "cds_length",
  "strand",
  "ensembl_transcript_id",
  "exon_chrom_start",
  "exon_chrom_end"
)

# Initialize output file (remove if it exists to ensure fresh start)
if (file.exists(opt$outputfile)) {
  file.remove(opt$outputfile)
}

# Process genes in batches of 250
batch_size <- 250
batches <- split(genes, ceiling(seq_along(genes)/batch_size))
has_written_data <- FALSE

for (i in seq_along(batches)) {
  message(sprintf("Processing batch %d of %d", i, length(batches)))
  
  # Add a 30 second delay between queries (except for the first batch) to avoid overloading the server
  if (i > 1) {
    Sys.sleep(30)
  }
  
  # Fetch data with biomaRt, retrying up to 3 times on failure
  batch_data <- NULL
  success <- FALSE
  attempts <- 0
  max_attempts <- 3
  
  while (!success && attempts < max_attempts) {
    attempts <- attempts + 1
    
    batch_data <- tryCatch({
      getBM(
        attributes = attributes,
        filters = c("chromosome_name", "biotype", "external_gene_name", "transcript_is_canonical"),
        values = list(
          c(as.character(1:22), "X", "Y"),
          "protein_coding",
          batches[[i]],
          TRUE
        ),
        mart = ensembl
      )
    }, error = function(e) {
      message(sprintf("Attempt %d failed: %s", attempts, e$message))
      NULL
    })
    
    if (!is.null(batch_data)) {
      success <- TRUE
    } else {
      if (attempts < max_attempts) {
        message("Retrying in 30 seconds...")
        Sys.sleep(30)
      } else {
        stop("Error: Failed to retrieve data from BioMart for batch ", i, " after ", max_attempts, " attempts.")
      }
    }
  }
  
  # Check if response isn't empty
  if (nrow(batch_data) > 0) {
    # Remove rows with empty CDS information
    exon_file <- subset(batch_data, !is.na(genomic_coding_start) & genomic_coding_start != "")
    
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