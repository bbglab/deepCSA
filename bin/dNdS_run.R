#!/opt/conda/bin/Rscript --vanilla

## USAGE
# Rscript dNdScv_single_sample.R --inputfile ../../0initial_processing/data/PILOT5/custom_files/tws/all_muts/all_below035_4dNdScv.txt --outputfile ../results/all_below035.tsv --samplename all --genelist genes.txt --genedepth genes_coverage.txt


library(optparse)
library(dndscv)

is_SNV <- function(x){
  if (x == "A" || x == "C" || x == "T" || x == "G"){
    return(TRUE)
  } else {
    return(FALSE)
  }
}



option_list = list(
  make_option(c("-n", "--samplename"), type="character", default=NULL,
              help="sample name/identifier of the run", metavar="character"),
  make_option(c("-i", "--inputfile"), type="character", default=NULL,
              help="mutation dataset file name", metavar="character"),
  make_option(c("-o", "--outputprefix"), type="character", default="dNdScv_output",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-r", "--referencetranscripts"), type="character",
              default="RefCDS.rda",
              help="Annotation reference file [default= %default]", metavar="character"),
  make_option(c("-c", "--covariates"), type="character",
              default="covariates_hg19_hg38_epigenome_pcawg.rda",
              help="Human GRCh38 covariates file [default= %default]", metavar="character"),
  make_option(c("-g", "--genelist"), type="character",
              default=NULL,
              help="Gene list file [default= %default]", metavar="character"),
  make_option(c("-d", "--genedepth"), type="character",
              default=NULL,
              help="Gene depth file (2 columns: GENE\tAVG_DEPTH) [default= %default]", metavar="character"),
  make_option(c("-s", "--snvsonly"), type="logical",
              default=FALSE,
              help="Only use SNVs for the analysis [default= %default]", metavar="logical")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);






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

# CDKN2A.p16INK4a
# Loads the covs object
load(opt$covariates)
load(opt$referencetranscripts)

# Remove CDKN2A.p14arf row and
#   rename CDKN2A.p16INK4a to CDKN2A
if ("CDKN2A.p14arf" %in% rownames(covs)) {
  covs <- covs[rownames(covs) != "CDKN2A.p14arf", ]
}
if ("CDKN2A.p16INK4a" %in% rownames(covs)) {
  rownames(covs)[rownames(covs) == "CDKN2A.p16INK4a"] <- "CDKN2A"
}

reference_genes <- intersect(
  unique(rownames(covs)),
  unique(gr_genes$names)
  )

# Identify genes that are in 'genes' but not in the row names of 'covs'
missing_genes <- setdiff(genes, reference_genes)

# Print the missing genes, if any
if (length(missing_genes) > 0) {
  print("These genes are in the 'genes' list but not in 'covs':")
  print(missing_genes)
} else {
  print("All requested genes are present in 'covs'.")
}

# Check that all the "requested" genes are in the covariates file
genes <- intersect(reference_genes, genes)
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
                file = paste(opt$outputprefix, '.dNdScv.cv.tsv', sep = ''),
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
              file = paste(opt$outputprefix, '.dNdScv.globaldnds.tsv', sep = ''),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}

# Handle dndsout$sel_loc
dnds_genes <- dndsout$sel_loc
if (!is.null(dnds_genes) && nrow(dnds_genes) > 0) {
  dnds_genes <- cbind(list("sample" = opt$samplename), dnds_genes)

  write.table(dnds_genes,
              file = paste(opt$outputprefix, '.dNdScv.loc.tsv', sep = ''),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}
