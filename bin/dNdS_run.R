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
  make_option(c("-o", "--outputfile"), type="character", default=NULL,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-r", "--referencetranscripts"), type="character",
              default="/workspace/projects/prominent/analysis/dNdScv/data/reference_files/RefCDS_human_latest_intogen.rda",
              help="Annotation reference file [default= %default]", metavar="character"),
  make_option(c("-c", "--covariates"), type="character",
              default="/workspace/projects/prominent/analysis/dNdScv/data/reference_files/covariates_hg19_hg38_epigenome_pcawg.rda",
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






# Loads the covs object
load(opt$covariates)
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


dnds_genes = dndsout$sel_cv

if ("nbregind" %in% names(dndsout)) {
  dnds_genes = cbind(list("theta_ind" = dndsout$nbregind$theta), dnds_genes)
}

dnds_genes = cbind(list("theta" = dndsout$nbreg$theta), dnds_genes)
dnds_genes = cbind(list("sample" = opt$samplename), dnds_genes)


write.table(dnds_genes,
            file = opt$outputfile,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)




dnds_genes = dndsout$globaldnds
dnds_genes = cbind(list("sample" = opt$samplename), dnds_genes)


write.table(dnds_genes,
            file = paste(opt$outputfile, 'globaldnds', sep = ''),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)




dnds_genes = dndsout$sel_loc

dnds_genes = cbind(list("sample" = opt$samplename), dnds_genes)


write.table(dnds_genes,
            file = paste(opt$outputfile, 'loc', sep = ''),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
