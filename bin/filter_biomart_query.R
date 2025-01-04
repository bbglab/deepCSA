#!/opt/conda/bin/Rscript --vanilla

## USAGE
# Rscript dNdScv_single_sample.R --inputfile ../../0initial_processing/data/PILOT5/custom_files/tws/all_muts/all_below035_4dNdScv.txt --outputfile ../results/all_below035.tsv --samplename all --genelist genes.txt --genedepth genes_coverage.txt

# Load required libraries
library(optparse)
library(GenomicRanges)
library(dplyr)


option_list = list(
  make_option(c("-r", "--bedfile"), type="character", default=NULL,
              help="BED file to use for subsetting the regions",
              metavar="character"),
  make_option(c("-b", "--biomartoutput"), type="character", default=NULL,
              help="EnsemblBioMart output file", metavar="character"),
  make_option(c("-o", "--outputfile"), type="character", default=NULL,
              help="output file name [default= %default]", metavar="character"),
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);





# Read the exon file
exon_file <- read.delim(opt$biomartoutput, header = TRUE, stringsAsFactors = FALSE)

# Convert exon file to a GRanges object
exons_gr <- GRanges(
  seqnames = exon_file$Chromosome.scaffold.name,
  ranges = IRanges(start = exon_file$Genomic.coding.start, end = exon_file$Genomic.coding.end),
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
    CDS_start = cumsum(width(ranges)) - width(ranges) + 1,
    CDS_end = cumsum(width(ranges)),
    CDS_length = sum(width(ranges)),
    Exon.region.start..bp. = start(ranges),
    Exon.region.end..bp. = end(ranges)
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





# if a file with coverage per gene  file is provided use that list of genes
# otherwise use all genes
if ( !is.null(opt$genedepth) ){
  # read the file into a character vector
  genes_coverage <- read.table(opt$genedepth, header = FALSE, col.names = c("GENE", "AVG_DEPTH"))
  genes_coverage <- unique.data.frame(genes_coverage)
  genes_with_info <- genes_coverage$GENE

  # create a named vector of the mean coverage values
  mean_coverage <- setNames(genes_coverage$AVG_DEPTH, genes_coverage$GENE)
  print(paste("Running dNdS with duplex coverage information."))

} else {
  mean_coverage = NULL
  genes_with_info = NULL
  print(paste("Running dNdS without information on duplex coverage."))
}


# if a genelist file is provided use that list of genes
# otherwise use all genes
if (!is.null(opt$genelist)){
  # read the file into a character vector
  genes <- readLines(opt$genelist)
  # remove empty strings due to empty lines
  genes <- genes[nzchar(genes)]

  if (!is.null(genes_with_info)){
    genes <- intersect(genes_with_info, genes)
    print("Keeping only the genes with information on duplex coverage")
  }

  print(paste("Running targeted dNdS in", length(genes), "genes."))
} else {
  genes = genes_with_info

  if (!is.null(genes_with_info)){
    print("Only the genes with information on duplex coverage")
  } else {
    print(paste("Running dNdS for all genes."))
  }

}


# Loads the covs object
load(opt$covariates)

# Identify genes that are in 'genes' but not in the row names of 'covs'
missing_genes <- setdiff(genes, rownames(covs))

# Print the missing genes, if any
if (length(missing_genes) > 0) {
  print("These genes are in the 'genes' list but not in 'covs':")
  print(missing_genes)
} else {
  print("All requested genes are present in 'covs'.")
}

# Check that all the "requested" genes are in the covariates file
genes <- intersect(rownames(covs), genes)
print("Keeping only the genes with in the covariates")


transcripts_file = opt$referencetranscripts


genes_sample_overlap = data.frame()
muts = read.table(opt$inputfile, sep = "\t", header = F)
print("##")
print(opt$samplename)
print("##")

colnames(muts) = c("sampleID", "chr", "pos", "ref", "mut")
dim(muts)

# Now filter to keep only the SNV
ourSNVs = muts[sapply(muts$ref, is_SNV), ]
ourSNVs = ourSNVs[sapply(ourSNVs$mut, is_SNV), ]
if (opt$snvsonly) {
  muts = ourSNVs
  print(paste(dim(ourSNVs)[1], "SNVs"))
} else {
  print(paste(dim(ourSNVs)[1], "SNVs"))
  print(paste(dim(muts)[1] - dim(ourSNVs)[1], "indels"))

}



dndsout = dndscv(muts,
                 refdb=transcripts_file,
                 gene_list=genes,
                 cv = covs,
                 max_muts_per_gene_per_sample = Inf,
                 max_coding_muts_per_sample = Inf,
                 dc = mean_coverage
)


# Check if dndsout$sel_cv is non-empty before proceeding
dnds_genes <- dndsout$sel_cv
if (!is.null(dnds_genes) && nrow(dnds_genes) > 0) {
  # Add theta_ind if nbregind exists in dndsout
  if ("nbregind" %in% names(dndsout)) {
    if (!is.null(dndsout$nbregind) && !is.null(dndsout$nbregind$theta)) {
      dnds_genes <- cbind(list("theta_ind" = dndsout$nbregind$theta), dnds_genes)
    }
  }

  # Add theta
  if (!is.null(dndsout$nbreg$theta)) {
    dnds_genes <- cbind(list("theta" = dndsout$nbreg$theta), dnds_genes)
  }

  # Add sample name
  dnds_genes <- cbind(list("sample" = opt$samplename), dnds_genes)

  # Write to file if dnds_genes is still valid
  if (nrow(dnds_genes) > 0) {
    write.table(dnds_genes,
                file = opt$outputfile,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
  }
}

# Handle dndsout$globaldnds
dnds_genes <- dndsout$globaldnds
if (!is.null(dnds_genes) && nrow(dnds_genes) > 0) {
  dnds_genes <- cbind(list("sample" = opt$samplename), dnds_genes)

  write.table(dnds_genes,
              file = paste(opt$outputfile, 'globaldnds', sep = ''),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}

# Handle dndsout$sel_loc
dnds_genes <- dndsout$sel_loc
if (!is.null(dnds_genes) && nrow(dnds_genes) > 0) {
  dnds_genes <- cbind(list("sample" = opt$samplename), dnds_genes)

  write.table(dnds_genes,
              file = paste(opt$outputfile, 'loc', sep = ''),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}
