#!/usr/bin/env python


import click
import pandas as pd
from utils import add_filter
from read_utils import custom_na_values

@click.command()
@click.option('--maf-df-file', required=True, type=click.Path(exists=True), help='Input gzipped MAF file (TSV)')
@click.option('--sample-name', required=True, type=str, help='Sample name for output file')
@click.option('--repetitive-variant-threshold', required=True, type=int, help='Threshold for repetitive variants')
@click.option('--somatic-vaf-boundary', required=True, type=float, help='VAF boundary for somatic variants')
@click.option('--n-rich-cohort-proportion', required=True, type=float, help='Proportion for n-rich cohort filtering')
def main(maf_df_file, sample_name, repetitive_variant_threshold, somatic_vaf_boundary, n_rich_cohort_proportion):
    maf_df = pd.read_csv(maf_df_file, compression='gzip', header=0, sep='\t', na_values=custom_na_values)
    sequenced_genes = list(pd.unique(maf_df["SYMBOL"]))


    #######
    ###  Filter repetitive variants,
    ###     both based on frequency and including information on position in read
    #######

    max_samples = len(pd.unique(maf_df["SAMPLE_ID"]))

    n_samples = list(range(repetitive_variant_threshold, max_samples + 1))
    if len(n_samples) == 0:
        print("Not enough samples to identify potential repetitive variants!")

    else:

        # work with already filtered df + somatic only to explore potential artifacts
        # take only variant and sample info from the df
        maf_df_f_somatic = maf_df.loc[maf_df["VAF"] <= somatic_vaf_boundary][["MUT_ID","SAMPLE_ID", "PMEAN", "PSTD"]].reset_index(drop = True)

        # add counter column
        maf_df_f_somatic["count"] = 1
        maf_df_f_somatic_pivot = maf_df_f_somatic.groupby("MUT_ID")["count"].sum().reset_index()

        repetitive_variants = maf_df_f_somatic_pivot[maf_df_f_somatic_pivot["count"] >= repetitive_variant_threshold]["MUT_ID"]
        print("Repetitive variants: ", len(repetitive_variants))

        maf_df["repetitive_variant"] = maf_df["MUT_ID"].isin(repetitive_variants)

        maf_df["FILTER"] = maf_df[["FILTER","repetitive_variant"]].apply(lambda x: add_filter(x["FILTER"], x["repetitive_variant"], "repetitive_variant"),
                                                                            axis = 1
                                                                        )
        maf_df = maf_df.drop("repetitive_variant", axis = 1)



        # use the position in read information to filter repetitive variants with a fixed position (likely artifacts)
        maf_df_f_somatic_pos_info = maf_df_f_somatic[~(maf_df_f_somatic["PMEAN"].isna()) & 
                                                        (maf_df_f_somatic["PMEAN"] != -1) &
                                                        (maf_df_f_somatic["PSTD"] == 0)]
        
        if maf_df_f_somatic_pos_info.shape[0] > 0:
            maf_df_f_somatic_compiled_pos = maf_df_f_somatic_pos_info.groupby("MUT_ID")["PMEAN"].nunique().reset_index()

            variants_with_rep_position = maf_df_f_somatic_compiled_pos[(maf_df_f_somatic_compiled_pos["PMEAN"] == 1)]["MUT_ID"]
            print("Variants always found in the same position: ", len(variants_with_rep_position))

            variants_with_rep_position = set(variants_with_rep_position).intersection(set(repetitive_variants))
            print("Repetitive variants always found in the same position: ", len(variants_with_rep_position))

            maf_df["repetitive_mapping_variant"] = maf_df["MUT_ID"].isin(variants_with_rep_position)

            maf_df["FILTER"] = maf_df[["FILTER","repetitive_mapping_variant"]].apply(lambda x: add_filter(x["FILTER"], x["repetitive_mapping_variant"], "repetitive_mapping_variant"),
                                                                                axis = 1
                                                                            )
            maf_df = maf_df.drop("repetitive_mapping_variant", axis = 1)




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
    germline_vars_all_samples = maf_df.loc[(maf_df["VAF"] > somatic_vaf_boundary) &
                                        (maf_df["VAF_AM"] > somatic_vaf_boundary) &
                                        (maf_df["vd_VAF"] > somatic_vaf_boundary),
                                        "MUT_ID"].unique()
    print(len(germline_vars_all_samples), "using all germline variants of all samples")


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
        maf_df["gnomAD_SNP"] = maf_df["gnomAD_SNP"].replace({"True": True, "False": False, '-' : False}).fillna(False).astype(bool)
        print("Out of ", maf_df["gnomAD_SNP"].shape[0], "positions", maf_df["gnomAD_SNP"].sum(), "are gnomAD SNPs (>0.1)")
        maf_df["FILTER"] = maf_df[["FILTER","gnomAD_SNP"]].apply(
                                                                    lambda x: add_filter(x["FILTER"], x["gnomAD_SNP"], "gnomAD_SNP"),
                                                                    axis = 1
                                                                )


    for filt in pd.unique(maf_df["FILTER"].str.split(";").explode()):
        maf_df[f"FILTER.{filt}"] = maf_df["FILTER"].apply(lambda x: filt in x.split(";"))

    for filtt in [ "not_covered", "not_in_exons"]:
        if f"FILTER.{filtt}" not in maf_df.columns:
            maf_df[f"FILTER.{filtt}"] = False


    maf_df.to_csv(f"{sample_name}.cohort.filtered.tsv.gz",
                  sep = "\t",
                  header = True,
                  index = False)

if __name__ == '__main__':
    main()
