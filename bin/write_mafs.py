#!/usr/bin/env python


import click
import pandas as pd
import json
from read_utils import custom_na_values

@click.command()
@click.option('--maf-file', required=True, type=click.Path(exists=True), help='Input gzipped MAF file (TSV)')
@click.option('--groups-json', required=True, type=click.Path(exists=True), help='Optional JSON file with group/sample mapping')
def main(maf_file, groups_json):
    maf_df = pd.read_csv(maf_file, compression='gzip', header=0, sep='\t', na_values=custom_na_values)
    maf_df["SAMPLE_ID"] = maf_df["SAMPLE_ID"].astype(str)

    with open(groups_json, 'r') as file:
        groups_info = json.load(file)

    for group_name, samples in groups_info.items():
        samples = [str(x) for x in samples]
        maf_df[maf_df["SAMPLE_ID"].isin(samples)].sort_values(by=["CHROM", "POS"]).to_csv(
            f"{group_name}.filtered.tsv.gz",
            sep="\t",
            header=True,
            index=False
        )

if __name__ == '__main__':
    main()
