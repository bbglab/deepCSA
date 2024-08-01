#!/usr/local/bin/python

import sys
import pandas as pd
import json

maf_df_file = sys.argv[1]

maf_df = pd.read_csv(maf_df_file, compression='gzip', header = 0, sep='\t', na_filter = False)  # Read gzipped TSV
maf_df["SAMPLE_ID"] = maf_df["SAMPLE_ID"].astype(str)

if len(sys.argv) > 2:
    json_file = sys.argv[2]

    with open(json_file, 'r') as file:
            groups_info = json.load(file)

    for group_name, samples in groups_info.items():
        samples = [ str(x) for x in samples ]
        maf_df[maf_df["SAMPLE_ID"].isin(samples)].sort_values(by = ["CHROM", "POS"]).to_csv(f"{group_name}.filtered.tsv.gz",
                                                                                                sep = "\t",
                                                                                                header = True,
                                                                                                index = False)


else:
    unique_samples = pd.unique(maf_df["SAMPLE_ID"])
    unique_samples = [ str(x) for x in unique_samples ]

    for sample in unique_samples:
        maf_df[ maf_df["SAMPLE_ID"] == sample].to_csv(f"{sample}.filtered.tsv.gz",
                                                        sep = "\t",
                                                        header = True,
                                                        index = False)

    maf_df.sort_values(by = ["CHROM", "POS"]).to_csv(f"all_samples.filtered.tsv.gz",
                                                        sep = "\t",
                                                        header = True,
                                                        index = False)
