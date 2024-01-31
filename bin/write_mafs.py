#!/usr/local/bin/python

import sys
import pandas as pd

maf_df_file = sys.argv[1]

maf_df = pd.read_csv(maf_df_file, compression='gzip', header = 0, sep='\t')  # Read gzipped TSV

unique_samples = pd.unique(maf_df["SAMPLE_ID"])

for sample in unique_samples:
    maf_df[maf_df["SAMPLE_ID"] == sample].to_csv(f"{sample}.filtered.tsv.gz",
                                                    sep = "\t",
                                                    header = True,
                                                    index = False)

maf_df.sort_values(by = ["CHROM", "POS"]).to_csv(f"all_samples.filtered.tsv.gz",
                                                    sep = "\t",
                                                    header = True,
                                                    index = False)
