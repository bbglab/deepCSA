library("stringr")

args = commandArgs(trailingOnly=TRUE)
input_file <- args[1]
out_csv_name <- args[2]
out_rds_name <- args[3]

df = read.table(input_file, header =TRUE, sep=",")
df$sample = df$SAMPLE_ID
df$individual = df$SAMPLE_ID
df$timepoint = "L"

df_result = df[,c("sample", "individual", "timepoint")]
#df_result = df[,c("sample")]

write.table(df_result,out_csv_name, row.names=FALSE, sep="\t",quote = FALSE)
saveRDS(df_result,file=out_rds_name)
