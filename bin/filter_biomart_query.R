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
    make_option(c("-r", "--bedfile"), type="character", default=NULL,
                help="BED file to use for subsetting the regions",
                metavar="character"),
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

# # Read the exon file
# exon_file <- read.delim(opt$biomartoutput, header = TRUE, stringsAsFactors = FALSE)

# Convert exon file to a GRanges object
exons_gr <- GRanges(
  seqnames = exon_file$Chromosome.scaffold.name,
  ranges = IRanges(start = as.integer(exon_file$Genomic.coding.start),
                    end = as.integer(exon_file$Genomic.coding.end)),
  strand = exon_file$Strand,
  gene_id = exon_file$Gene.stable.ID,
  gene_name = exon_file$Gene.name,
  protein_id = exon_file$Protein.stable.ID,
  transcript_id = exon_file$Transcript.stable.ID,
  CDS_start = exon_file$CDS.start,
  CDS_end = exon_file$CDS.end,
  CDS_length = exon_file$CDS.Length,
  exon_region_start = exon_file$Exon.region.start..bp.,
  exon_region_end = exon_file$Exon.region.end..bp.
)

# Read the BED file with regions of interest
bed_file <- read.delim(opt$bedfile, header = FALSE, stringsAsFactors = FALSE)
colnames(bed_file) <- c("chrom", "start", "end")
bed_file$chrom <- gsub("^chr", "", bed_file$chrom)
regions_gr <- GRanges(seqnames = bed_file$chrom, ranges = IRanges(start = bed_file$start, end = bed_file$end))

# Find overlaps between exons and regions of interest
overlaps <- findOverlaps(exons_gr, regions_gr)

# Filter the exons that overlap the regions of interest
filtered_exons <- exons_gr[queryHits(overlaps)]

# Adjust CDS coordinates
adjusted_exons <- filtered_exons %>%
  as.data.frame() %>%
  group_by(transcript_id) %>%
  arrange(CDS_start) %>%
  mutate(
    CDS_start = cumsum(end - start + 1) - (end - start + 1) + 1,
    CDS_end = cumsum(end - start + 1),
    CDS_length = sum(end - start + 1),
    Exon.region.start..bp. = start,
    Exon.region.end..bp. = end
  ) %>%
  ungroup()

# Save the updated exon data
write.table(
  adjusted_exons,
  file = opt$outputfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
