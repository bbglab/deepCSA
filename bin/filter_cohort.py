#!/usr/local/bin/python

import sys
import pandas as pd
from utils import add_filter

maf_df_file = sys.argv[1]
samp_name = sys.argv[2]
repetitive_variant_treshold = int(sys.argv[3])
somatic_vaf_boundary = float(sys.argv[4])

maf_df = pd.read_csv(maf_df_file, compression='gzip', header = 0, sep='\t')  # Read gzipped TSV

sequenced_genes = list(pd.unique(maf_df["SYMBOL"]))



# def correct_vaf(maf):

#     """
#     Computes VAF_CORRECTED for the subset of variants satisfying 0 < VAF < 0.2
#     Returns the input MAF with two new columns:
#         VAF_CORRECTED with new corrected VAF else it copies the VAF
#         IS_VAF_CORRECTED with a boolean that indicates whether the VAF has been corrected
#     """

#     # TODO revise the 0.2 VAF threshold to see if it can be kept across datasets
#     df = maf[(0 < maf['VAF']) & (maf['VAF'] < 0.2)][['SAMPLE_ID', 'MUT_ID', 'VAF', 'DEPTH']]
#     df = df.sort_values('DEPTH')
#     N  = df.shape[0]
#     df['VAF_ROLLING_MEAN'] = df['VAF'].rolling(N // 25).mean()
#     df['VAF_ROLLING_STD'] = df['VAF'].rolling(N // 25).std()
#     stable_mean = df['VAF_ROLLING_MEAN'].values[-1]
#     stable_std  = df['VAF_ROLLING_STD'].values[-1]
#     df['VAF_CORRECTED'] = df.apply(lambda r: (r['VAF'] - r['VAF_ROLLING_MEAN']) * (stable_std / r['VAF_ROLLING_STD']) + stable_mean, axis=1)
#     df = maf.merge(df[['VAF_CORRECTED', 'MUT_ID', 'SAMPLE_ID']],
#                                     on=['MUT_ID', 'SAMPLE_ID'],
#                                     how='outer')
#     df['IS_VAF_CORRECTED'] = ~df['VAF_CORRECTED'].isnull()
#     df.loc[~df['IS_VAF_CORRECTED'], 'VAF_CORRECTED'] = df[~df['IS_VAF_CORRECTED']]['VAF'].values
#     return df



# #######
# ###  Add a corrected VAF column
# #######
# maf_df = correct_vaf(maf_df)
# print("VAF corrected")




#######
###  Filter repetitive variants
#######

# TODO revise these numbers, the repetitive_variant_treshold is the boundary at which we start considering a mutation as "repetitive"
max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))

n_samples = list(range(repetitive_variant_treshold, max_samples + 1))
if len(n_samples) == 0:
    print("Not enough samples to identify potential repetitive variants!")

else:

    # work with already filtered df + somatic only to explore potential artifacts
    # take only variant and sample info from the df
    maf_df_f_somatic = maf_df.loc[maf_df["VAF"] <= somatic_vaf_boundary][["MUT_ID","SAMPLE_ID"]].reset_index(drop = True)

    # add counter column
    maf_df_f_somatic["count"] = 1
    maf_df_f_somatic_pivot = maf_df_f_somatic.groupby("MUT_ID")["count"].sum().reset_index()

    repetitive_variants = maf_df_f_somatic_pivot[maf_df_f_somatic_pivot["count"] >= repetitive_variant_treshold]["MUT_ID"]

    maf_df["repetitive_variant"] = maf_df["MUT_ID"].isin(repetitive_variants)

    maf_df["FILTER"] = maf_df[["FILTER","repetitive_variant"]].apply(lambda x: add_filter(x["FILTER"], x["repetitive_variant"], "repetitive_variant"),
                                                                        axis = 1
                                                                    )
    maf_df = maf_df.drop("repetitive_variant", axis = 1)







#######
###  Filter other sample's SNP
#######

# this is if we were to consider both unique and no-unique variants
germline_vars_all_samples = maf_df.loc[maf_df["VAF"] > somatic_vaf_boundary, "MUT_ID"].unique()

# this is if we were to consider only unique germline
germline_vars = maf_df.loc[maf_df["VAF"] > somatic_vaf_boundary][["MUT_ID", "SAMPLE_ID"]].groupby(
                                                            "MUT_ID").size().sort_values(ascending = False).to_frame("n_samples").reset_index()
unique_germline_vars = germline_vars.loc[germline_vars["n_samples"] == 1]["MUT_ID"].unique()



print(len(germline_vars_all_samples), "using all germline variants of all samples")
print(len(unique_germline_vars), "using only germline variants unique to a single sample")

maf_df["other_sample_SNP"] = False
maf_df.loc[(maf_df["MUT_ID"].isin(germline_vars_all_samples)) &
            (maf_df["VAF"] <= somatic_vaf_boundary), "other_sample_SNP"] = True

maf_df["FILTER"] = maf_df[["FILTER","other_sample_SNP"]].apply(
                                                lambda x: add_filter(x["FILTER"], x["other_sample_SNP"], "other_sample_SNP"),
                                                axis = 1
                                            )
maf_df = maf_df.drop("other_sample_SNP", axis = 1)


for filt in pd.unique(maf_df["FILTER"].str.split(";").explode()):
    maf_df[f"FILTER.{filt}"] = maf_df["FILTER"].str.contains(filt)


maf_df.to_csv(f"{samp_name}.cohort.filtered.tsv.gz",
                        sep = "\t",
                        header = True,
                        index = False)
