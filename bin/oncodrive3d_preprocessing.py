#!/usr/local/bin/python

import sys
import pandas as pd


maf_df_file = sys.argv[1]
vep_output_all = sys.argv[2]
sample = sys.argv[3]


vep_data = pd.read_table(vep_output_all, na_filter = False)
mutationsdata = pd.read_table(maf_df_file, na_filter = False)

if "Tumor_Sample_Barcode" in mutationsdata.columns:
    reduced_mutationsdata = mutationsdata[["Tumor_Sample_Barcode", "MUT_ID"]]
else:
    reduced_mutationsdata = mutationsdata[["SAMPLE_ID", "MUT_ID"]]
    reduced_mutationsdata.columns = ["Tumor_Sample_Barcode", "MUT_ID"]

uniq_mut_ids = reduced_mutationsdata["MUT_ID"].unique()
all_vep_data = vep_data[vep_data["#Uploaded_variation"].isin(uniq_mut_ids)].reset_index(drop = True)

all_vep_data_sample = reduced_mutationsdata.merge(all_vep_data,
                                                    left_on=["MUT_ID"],
                                                    right_on=["#Uploaded_variation"],
                                                    how='left'
                                                    )

print(all_vep_data_sample.shape)
all_vep_data_sample.to_csv(f"{sample}.mutations.raw_vep.tsv", sep = '\t', header = True, index = False )

