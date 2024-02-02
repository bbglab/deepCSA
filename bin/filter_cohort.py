#!/usr/local/bin/python

import sys
import pandas as pd
from utils import add_filter

maf_df_file = sys.argv[1]
samp_name = sys.argv[2]
intogen_mutations = sys.argv[3]

maf_df = pd.read_csv(maf_df_file, compression='gzip', header = 0, sep='\t')  # Read gzipped TSV

sequenced_genes = list(pd.unique(maf_df["SYMBOL"]))
somatic_vaf_boundary = 0.35



#######
###  Filter repetitive variants
#######

# TODO revise these numbers, the repetitive_variant_treshold is the boundary at which we start considering a mutation as "artifactual"
repetitive_variant_treshold = 2
max_samples = len(pd.unique(maf_df["SAMPLE_ID"])) + 1

n_samples = list(range(repetitive_variant_treshold, max_samples))
if len(n_samples) == 0:
    print("Not enough samples to identify potential repetitive variants!")

else:
    # Repetitive variants
    intogen_muts_df = pd.read_csv(intogen_mutations, sep = "\t")

    # only keep panel muts
    intogen_muts_df = intogen_muts_df.loc[intogen_muts_df["SYMBOL"].isin(sequenced_genes)].reset_index(drop = True)
    intogen_muts_df["CHR"] = intogen_muts_df["CHR"].astype(str)

    # edit mutation field so that the format is the same as in the MAF file
    intogen_muts_df["MUTATION_MAFformat"] = intogen_muts_df.apply(
        lambda row: "chr"+row["MUTATION"].split(":")[0]+":"+row["MUTATION"].split(":")[1]+"_"+row["MUTATION"].split(":")[2], axis = 1)





    ####
    #  SNVs
    ####


    ##
    # Revise these functions based on the preprocessed information already available in the MAF file
    ##

    # edit variant IDs in MAF file
    # I am finally using a function based on Ferriol's approach which I think is safer
    # he takes the Ensembl info in Allele, which is the Ensemble REF and Location, which is the Ensembl updated location
    def correct_MAFid_indels(row):

        # is an insertion
        if row["TYPE"] == "INSERTION":

            # change REF to "-" and omit first nucleotide in ALT
            return row["CHROM"]+":"+row["Location"].split(":")[1].split("-")[0]+"_->"+row["Allele"]

        # is a deletion
        elif row["TYPE"] == "DELETION":

            # change ALT to "-" and omit first nucleotide in REF
            return row["CHROM"]+":"+row["Location"].split(":")[1].split("-")[0]+"_"+row["REF"][len(row["ALT"]):]+">"+row["Allele"]

        return row["MUT_ID"]



    maf_df["MUT_ID_IntoGen_format"] = maf_df.apply(lambda row: correct_MAFid_indels(row), axis = 1)


    # save number of samples per gene to use afterwards to calculate percentages
    intogen_muts_samples_panel_df = intogen_muts_df.groupby("SYMBOL")["SAMPLES"].sum().to_frame("n_samples_intogen").reset_index()
    intogen_muts_samples_panel_df = intogen_muts_samples_panel_df.loc[intogen_muts_samples_panel_df["SYMBOL"].isin(sequenced_genes)]





    # work with already filtered df + somatic only to explore potential artifacts
    # take only variant and sample info from the df
    maf_df_f_somatic = maf_df.loc[(maf_df["VAF"] < somatic_vaf_boundary)
                                    & (maf_df["TYPE"] == "SNV")
                                    ][["MUT_ID_IntoGen_format","SAMPLE_ID"]].reset_index(drop = True)

    # add counter column
    maf_df_f_somatic["count"] = 1
    # pivot by variation and sample
    maf_df_f_somatic_pivot = maf_df_f_somatic.loc[maf_df_f_somatic["MUT_ID_IntoGen_format"].notna()].pivot(
        index = 'MUT_ID_IntoGen_format', columns = 'SAMPLE_ID', values = "count")
    # reorder
    # maf_df_f_somatic_pivot = maf_df_f_somatic_pivot[sample_order]
    # sum to get the number of times the same mutation appears in the cohort
    maf_df_f_somatic_pivot["n_samples_sequenced"] = maf_df_f_somatic_pivot.notnull().sum(axis = 1)




    df_dict = {}
    for n in n_samples:

        # take mutations appearing in n samples
        rep_muts_n_df = maf_df_f_somatic_pivot.loc[maf_df_f_somatic_pivot["n_samples_sequenced"] == n]["n_samples_sequenced"].reset_index()
        rep_muts_n_df.columns = ["MUTATION_MAFformat", "n_samples_sequenced"]
        rep_muts_n = rep_muts_n_df["MUTATION_MAFformat"]

        # annotate consequence and symbol from maf_df to identify the gene in the panel
        symb_df = maf_df.loc[maf_df["MUT_ID_IntoGen_format"].isin(rep_muts_n)
                            ][["MUT_ID_IntoGen_format", "SYMBOL", "Consequence_single", "Consequence_broader"]].drop_duplicates()
        symb_df.columns = ["MUTATION_MAFformat", "SYMBOL", "Consequence_single", "Consequence_broader"]
        rep_muts_n_df = rep_muts_n_df.merge(symb_df, on = "MUTATION_MAFformat")

        # look in intogen the number of samples in which this mutation appears
        intogen_muts_df_sset = intogen_muts_df.loc[intogen_muts_df["MUTATION_MAFformat"].isin(rep_muts_n)][["MUTATION_MAFformat", "SYMBOL", "SAMPLES"]]
        intogen_muts_df_sset_gby = intogen_muts_df_sset.groupby(["MUTATION_MAFformat", "SYMBOL"])["SAMPLES"].sum().to_frame("n_samples_intogen").reset_index()
        intogen_muts_df_sset_gby["n_samples_sequenced"] = n
        intogen_muts_df_sset_gby = rep_muts_n_df.merge(intogen_muts_df_sset_gby, how = "left", on = ["MUTATION_MAFformat", "SYMBOL", "n_samples_sequenced"])

        df_dict[n] = intogen_muts_df_sset_gby

    potential_artifact_df_snv = pd.concat(df_dict.values()).sort_values(by = ["n_samples_sequenced"], ascending = False)





    ####
    #  Indels
    ####
    # work with already filtered df + somatic only to explore potential artifacts
    # take only variant and sample info from the df

    maf_df_f_somatic = maf_df.loc[(maf_df["VAF"] < somatic_vaf_boundary)
                                    & (maf_df["TYPE"].isin(["DELETION", "INSERTION"]))
                                    ][["MUT_ID_IntoGen_format", "SAMPLE_ID"]].reset_index(drop = True).drop_duplicates() # drop duplis because there are some variants annot as indels which are not
    maf_df_f_somatic["count"] = 1
    maf_df_f_somatic_pivot = maf_df_f_somatic.pivot(index = 'MUT_ID_IntoGen_format', columns = 'SAMPLE_ID', values = "count")
    # maf_df_f_somatic_pivot = maf_df_f_somatic_pivot[sample_order]
    maf_df_f_somatic_pivot["n_samples_sequenced"] = maf_df_f_somatic_pivot.notnull().sum(axis = 1)



    df_dict = {}

    for n in n_samples:

        rep_muts_n_df = maf_df_f_somatic_pivot.loc[maf_df_f_somatic_pivot["n_samples_sequenced"] == n]["n_samples_sequenced"].reset_index()
        rep_muts_n_df.columns = ["MUTATION_MAFformat", "n_samples_sequenced"]
        rep_muts_n = rep_muts_n_df["MUTATION_MAFformat"]

        symb_df = maf_df.loc[maf_df["MUT_ID_IntoGen_format"].isin(
            rep_muts_n)][["MUT_ID_IntoGen_format", "SYMBOL", "Consequence_single", "Consequence_broader"]].drop_duplicates()
        symb_df.columns = ["MUTATION_MAFformat", "SYMBOL", "Consequence_single", "Consequence_broader"]
        rep_muts_n_df = rep_muts_n_df.merge(symb_df, on = "MUTATION_MAFformat")

        intogen_muts_df_sset = intogen_muts_df.loc[intogen_muts_df["MUTATION_MAFformat"].isin(rep_muts_n)][["MUTATION_MAFformat", "SYMBOL", "SAMPLES"]]
        intogen_muts_df_sset_gby = intogen_muts_df_sset.groupby(["MUTATION_MAFformat", "SYMBOL"])["SAMPLES"].sum().to_frame("n_samples_intogen").reset_index()
        intogen_muts_df_sset_gby["n_samples_sequenced"] = n
        intogen_muts_df_sset_gby = rep_muts_n_df.merge(intogen_muts_df_sset_gby, how = "left", on = ["MUTATION_MAFformat", "SYMBOL", "n_samples_sequenced"])

        df_dict[n] = intogen_muts_df_sset_gby

    potential_artifact_df_indel = pd.concat(df_dict.values()).sort_values(by = ["n_samples_sequenced"], ascending = False)

    potential_artifact_muts = pd.concat([potential_artifact_df_snv, potential_artifact_df_indel]).drop_duplicates()["MUTATION_MAFformat"].values

    maf_df["repetitive_variant"] = maf_df["MUT_ID_IntoGen_format"].isin(potential_artifact_muts)


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
            (maf_df["VAF"] <= 0.35), "other_sample_SNP"] = True

maf_df["FILTER"] = maf_df[["FILTER","other_sample_SNP"]].apply(
                                                lambda x: add_filter(x["FILTER"], x["other_sample_SNP"], "other_sample_SNP"),
                                                axis = 1
                                            )
maf_df = maf_df.drop("other_sample_SNP", axis = 1)



maf_df.to_csv(f"{samp_name}.cohort.filtered.tsv.gz",
                        sep = "\t",
                        header = True,
                        index = False)
