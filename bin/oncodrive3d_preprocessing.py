#!/usr/bin/env python
# TODO: bump pandas to 2.2.3

import click
import pandas as pd
from read_utils import custom_na_values


@click.command()
@click.option('--maf-df-file', required=True, type=click.Path(exists=True), help='Input MAF dataframe file')
@click.option('--vep-output-all', required=True, type=click.Path(exists=True), help='Input VEP output file')
@click.option('--sample', required=True, type=str, help='Sample name for output file')
def main(maf_df_file, vep_output_all, sample):
    vep_data = pd.read_table(vep_output_all, na_values = custom_na_values)
    mutationsdata = pd.read_table(maf_df_file, na_values = custom_na_values)

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


if __name__ == '__main__':
    main()