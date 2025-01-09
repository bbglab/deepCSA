library(stringr)
args = commandArgs(trailingOnly=TRUE)
mutation_counts_file <- args[1]
out_csv_name <- args[2]
out_rds_name <- args[3]

mutation_counts <- read.table(mutation_counts_file, header = T)
rownames(mutation_counts) <- mutation_counts[["Mutation.Types"]]
mutation_counts <- mutation_counts[,c(2:ncol(mutation_counts))]
colnames(mutation_counts) <- str_split_i(colnames(mutation_counts), "[.]", 1)
print(mutation_counts[1:5, 1:5])
mutation_counts <- t(mutation_counts)
print(mutation_counts[1:5, 1:5])
write.table(mutation_counts,out_csv_name, row.names=TRUE, sep="\t",quote = FALSE)
saveRDS(mutation_counts, file = out_rds_name)
