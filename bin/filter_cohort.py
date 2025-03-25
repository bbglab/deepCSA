#!/usr/local/bin/python

import sys
import pandas as pd
from utils import add_filter
from read_utils import custom_na_values

maf_df_file = sys.argv[1]
samp_name = sys.argv[2]
repetitive_variant_threshold = int(sys.argv[3])
somatic_vaf_boundary = float(sys.argv[4])
n_rich_cohort_proportion = float(sys.argv[5])

maf_df = pd.read_csv(maf_df_file, compression='gzip', header = 0, sep='\t', na_values = custom_na_values)

sequenced_genes = list(pd.unique(maf_df["SYMBOL"]))





#######
###  Filter repetitive variants
#######

# TODO revise these numbers, the repetitive_variant_threshold is the boundary at which we start considering a mutation as "repetitive"
max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))

n_samples = list(range(repetitive_variant_threshold, max_samples + 1))
if len(n_samples) == 0:
    print("Not enough samples to identify potential repetitive variants!")

else:

    # work with already filtered df + somatic only to explore potential artifacts
    # take only variant and sample info from the df
    maf_df_f_somatic = maf_df.loc[maf_df["VAF"] <= somatic_vaf_boundary][["MUT_ID","SAMPLE_ID"]].reset_index(drop = True)

    # add counter column
    maf_df_f_somatic["count"] = 1
    maf_df_f_somatic_pivot = maf_df_f_somatic.groupby("MUT_ID")["count"].sum().reset_index()

    repetitive_variants = maf_df_f_somatic_pivot[maf_df_f_somatic_pivot["count"] >= repetitive_variant_threshold]["MUT_ID"]

    maf_df["repetitive_variant"] = maf_df["MUT_ID"].isin(repetitive_variants)

    maf_df["FILTER"] = maf_df[["FILTER","repetitive_variant"]].apply(lambda x: add_filter(x["FILTER"], x["repetitive_variant"], "repetitive_variant"),
                                                                        axis = 1
                                                                    )
    maf_df = maf_df.drop("repetitive_variant", axis = 1)




#######
###  Filter cohort_n_rich
#######

max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))

if max_samples < 2:
    print("Not enough samples to identify cohort_n_rich mutations!")

else:
    number_of_samples = max(2, (max_samples * n_rich_cohort_proportion) // 1)
    print(f"flagging mutations that are n_rich in at least: {number_of_samples} samples as cohort_n_rich")

    # work with already filtered df + somatic only to explore potential artifacts
    # take only variant and sample info from the df
    maf_df_f_somatic = maf_df[["MUT_ID", "SAMPLE_ID", "VAF_Ns", "FILTER"]].reset_index(drop = True)

    n_rich_vars_df = maf_df_f_somatic[maf_df_f_somatic["FILTER"].str.contains("n_rich")].groupby("MUT_ID")[
                                            ['SAMPLE_ID', 'VAF_Ns']
                                        ].agg({'SAMPLE_ID' : len, 'VAF_Ns' : min})
    n_rich_vars_df = n_rich_vars_df.rename({'SAMPLE_ID' : 'N_rich_frequency', 'VAF_Ns' : 'VAF_Ns_threshold'}, axis = 'columns')

    n_rich_vars = list(n_rich_vars_df[n_rich_vars_df['N_rich_frequency'] >= number_of_samples].index)

    maf_df["cohort_n_rich"] = maf_df["MUT_ID"].isin(n_rich_vars)

    maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich"], "cohort_n_rich"),
                                                                        axis = 1
                                                                    )
    maf_df = maf_df.drop("cohort_n_rich", axis = 1)



    # if the variant appeared flagged as n_rich in a single sample it is also filtered out from all other samples
    n_rich_vars_uni = list(n_rich_vars_df[n_rich_vars_df['N_rich_frequency'] > 0].index)

    maf_df["cohort_n_rich_uni"] = maf_df["MUT_ID"].isin(n_rich_vars_uni)

    maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich_uni"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich_uni"], "cohort_n_rich_uni"),
                                                                        axis = 1
                                                                    )
    maf_df = maf_df.drop("cohort_n_rich_uni", axis = 1)


    # if the variant appeared flagged as n_rich in a single sample it is also filtered out from all other samples
    maf_df = maf_df.merge(n_rich_vars_df, on = 'MUT_ID', how = 'left')
    maf_df['N_rich_frequency'] = maf_df['N_rich_frequency'].fillna(0)
    maf_df['VAF_Ns_threshold'] = maf_df['VAF_Ns_threshold'].fillna(1.1)

    maf_df["cohort_n_rich_threshold"] = maf_df["VAF_Ns"] >= maf_df['VAF_Ns_threshold']

    maf_df["FILTER"] = maf_df[["FILTER","cohort_n_rich_threshold"]].apply(lambda x: add_filter(x["FILTER"], x["cohort_n_rich_threshold"], "cohort_n_rich_threshold"),
                                                                        axis = 1
                                                                    )
    maf_df = maf_df.drop("cohort_n_rich_threshold", axis = 1)






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
# print(len(unique_germline_vars), "using only germline variants unique to a single sample")

maf_df["other_sample_SNP"] = False
maf_df.loc[(maf_df["MUT_ID"].isin(germline_vars_all_samples)) &
            (maf_df["VAF"] <= somatic_vaf_boundary), "other_sample_SNP"] = True

maf_df["FILTER"] = maf_df[["FILTER","other_sample_SNP"]].apply(
                                                lambda x: add_filter(x["FILTER"], x["other_sample_SNP"], "other_sample_SNP"),
                                                axis = 1
                                            )
maf_df = maf_df.drop("other_sample_SNP", axis = 1)



#######
###  Filter gnomad SNP
#######

if "gnomAD_SNP" in maf_df.columns:
    maf_df["FILTER"] = maf_df[["FILTER","gnomAD_SNP"]].apply(
                                                                lambda x: add_filter(x["FILTER"], x["gnomAD_SNP"], "gnomAD_SNP"),
                                                                axis = 1
                                                            )


for filt in pd.unique(maf_df["FILTER"].str.split(";").explode()):
    maf_df[f"FILTER.{filt}"] = maf_df["FILTER"].apply(lambda x: filt in x.split(";"))


maf_df.to_csv(f"{samp_name}.cohort.filtered.tsv.gz",
                        sep = "\t",
                        header = True,
                        index = False)
