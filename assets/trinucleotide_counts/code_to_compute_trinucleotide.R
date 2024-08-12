library(BSgenome.Hsapiens.UCSC.hg38)

chromosomes = paste0("chr", c(1:22, "X", "Y"))
chrom_counts = trinucleotideFrequency(getSeq(Hsapiens, chromosomes))
sums = colSums(chrom_counts)

purine_index = substr(names(sums), 2,2) %in% c("A", "G")

purines = sums[purine_index]
pyrimidines = sums[!purine_index]

names(purines) = names(purines) |> 
  DNAStringSet() |> 
  reverseComplement() |> 
  as.character()

trinuc_counts = purines[names(pyrimidines)] + pyrimidines
print(trinuc_counts)

trinuc_df <- data.frame(CONTEXT = names(trinuc_counts), COUNT = as.numeric(trinuc_counts), row.names = NULL)

# Write the data frame into a TSV file without quotes
write.table(trinuc_df, file = "trinuc_counts.homo_sapiens.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)






library(BSgenome.Mmusculus.UCSC.mm10)
chromosomes = paste0("chr", c(1:19, "X", "Y"))
chrom_counts = trinucleotideFrequency(getSeq(Mmusculus, chromosomes))
sums = colSums(chrom_counts)

purine_index = substr(names(sums), 2,2) %in% c("A", "G")

purines = sums[purine_index]
pyrimidines = sums[!purine_index]

names(purines) = names(purines) |> 
  DNAStringSet() |> 
  reverseComplement() |> 
  as.character()

trinuc_counts = purines[names(pyrimidines)] + pyrimidines
trinuc_df <- data.frame(CONTEXT = names(trinuc_counts), COUNT = as.numeric(trinuc_counts), row.names = NULL)

# Write the data frame into a TSV file without quotes
write.table(trinuc_df, file = "trinuc_counts.mus_musculus.mm10.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)






library(BSgenome.Mmusculus.UCSC.mm39)
chromosomes = paste0("chr", c(1:19, "X", "Y"))
chrom_counts = trinucleotideFrequency(getSeq(Mmusculus, chromosomes))
sums = colSums(chrom_counts)

purine_index = substr(names(sums), 2,2) %in% c("A", "G")

purines = sums[purine_index]
pyrimidines = sums[!purine_index]

names(purines) = names(purines) |>
  DNAStringSet() |>
  reverseComplement() |>
  as.character()

trinuc_counts = purines[names(pyrimidines)] + pyrimidines
trinuc_df <- data.frame(CONTEXT = names(trinuc_counts), COUNT = as.numeric(trinuc_counts), row.names = NULL)

# Write the data frame into a TSV file without quotes
write.table(trinuc_df, file = "trinuc_counts.mus_musculus.mm39.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

