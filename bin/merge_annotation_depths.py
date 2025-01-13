#!/usr/local/bin/python


import click
import json
import pandas as pd

def annotate_depths(annotation_file, depths_file, json_f, input_file):
    """
    INFO
    """

    raw_depths = pd.read_csv(depths_file, sep = "\t", header = 0)
    annotations = pd.read_csv(annotation_file, sep = "\t", header = 0)

    annotated_depths = raw_depths.merge(annotations, on = ["CHROM", "POS"], how = 'left')

    del raw_depths
    del annotations

    print(annotated_depths.head())

    annotated_depths["CONTEXT"] = annotated_depths["CONTEXT"].fillna('-')
    sample_columns = [x for x in annotated_depths.columns if x not in ["CHROM", "POS", "CONTEXT"] ]
    annotated_depths = annotated_depths[["CHROM", "POS", "CONTEXT"] + sample_columns]

    input_csv = pd.read_table(input_file, sep = ',', header = 0)
    bam2sample_dict = dict(zip(input_csv["bam"].astype(str).apply(lambda x: x.split("/")[-1]),
                               input_csv["sample"].astype(str)
                               )
                            )
    
    sample_columns_correct = [ bam2sample_dict[str(x)] for x in sample_columns ]
    annotated_depths.columns = ["CHROM", "POS", "CONTEXT"] + sample_columns_correct
    annotated_depths.to_csv("all_samples_indv.depths.tsv.gz",
                                header=True,
                                index=False,
                                sep="\t")

    if json_f:
        with open(json_f, 'r') as file:
            groups_info = json.load(file)

        for group_name, samples in groups_info.items():
            annotated_depths[group_name] = annotated_depths.loc[:,samples].sum(axis=1)
            annotated_depths[["CHROM", "POS", "CONTEXT", group_name]].to_csv(f"{group_name}.depths.annotated.tsv.gz",
                                                                                sep = "\t",
                                                                                header = True,
                                                                                index = False)

    else:
        for sample in sample_columns_correct:
            annotated_depths[["CHROM", "POS", "CONTEXT", f"{sample}"]].to_csv(f"{sample}.depths.annotated.tsv.gz",
                                                                                sep = "\t",
                                                                                header = True,
                                                                                index = False)

        annotated_depths["all_samples"] = annotated_depths.iloc[:,3:].sum(axis=1)
        annotated_depths[["CHROM", "POS", "CONTEXT", "all_samples"]].to_csv(f"all_samples.depths.annotated.tsv.gz",
                                                                            sep = "\t",
                                                                            header = True,
                                                                            index = False)



@click.command()
@click.option('--annotation', type=click.Path(exists=True), help='Input annotation file')
@click.option('--depths', type=click.Path(exists=True), help='Input depths file')
@click.option('--json_file', type=click.Path(exists=True), help='JSON groups file')
@click.option('--input_csv', type=click.Path(exists=True), help='Input CSV file')
# @click.option('--output', type=click.Path(), help='Output annotated depths file')

def main(annotation, depths, json_file, input_csv):
    click.echo(f"Annotating depths file...")
    annotate_depths(annotation, depths, json_file, input_csv)

if __name__ == '__main__':
    main()

