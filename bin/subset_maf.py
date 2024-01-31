#!/usr/local/bin/python


import click
import json
import pandas as pd

from utils import filter_maf


def subset_mutation_dataframe(sample_name, mutations_file, mutations_file_out, json_filters, requested_fields):
    """
    INFO
    """
    # Load your MAF DataFrame (raw_annotated_maf)
    raw_annotated_maf = pd.read_csv(mutations_file, sep = "\t", header = 0)

    # Load the filter criteria from the JSON file
    with open(json_filters, 'r') as file:
        filter_criteria = json.load(file)

    if len(filter_criteria) > 0:
        # Filter the annotated maf using the described filters
        annotated_maf = filter_maf(raw_annotated_maf, filter_criteria)
        print("MAF subset")
    else:
        annotated_maf = raw_annotated_maf

    # Load the filter criteria from the JSON file
    with open(requested_fields, 'r') as file:
        output_format = json.load(file)

    if output_format:
        header_ = output_format.get("header", False)
        columns = output_format.get("columns", annotated_maf.columns.values)
        colnames = output_format.get("colnames", columns)

        if header_ :
            annotated_maf[columns].to_csv(mutations_file_out,
                                        header=colnames,
                                        index=False,
                                        sep="\t")
        else:
            annotated_maf[columns].to_csv(mutations_file_out,
                                        header=None,
                                        index=False,
                                        sep="\t")

    else:
        annotated_maf.to_csv(mutations_file_out,
                                header=True,
                                index=False,
                                sep="\t")



@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--mut_file', type=click.Path(exists=True), help='Input mutation file')
@click.option('--out_maf', type=click.Path(), help='Output MAF file')
@click.option('--json_filters', type=click.Path(exists=True), help='Input mutation filtering criteria file')
@click.option('--req_fields', type=click.Path(exists=True), help='Column names to output')
# @click.option('--plot', is_flag=True, help='Generate plot and save as PDF')

def main(sample_name, mut_file, out_maf, json_filters, req_fields): # , plot):
    click.echo(f"Subsetting MAF file...")
    subset_mutation_dataframe(sample_name, mut_file, out_maf, json_filters, req_fields)

if __name__ == '__main__':
    main()

