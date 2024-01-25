#!/usr/bin/env python3


import click

import pandas as pd




def annotate_depths(annotation_file, depths_file, output):
    """
    INFO
    """
    
    raw_depths = pd.read_csv(depths_file, sep = "\t", header = 0)
    annotations = pd.read_csv(annotation_file, sep = "\t", header = 0)

    annotated_depths = raw_depths.merge(annotations, on = ["CHROM", "POS"], how = 'left')
    
    del raw_depths
    del annotations

    annotated_depths["CONTEXT"] = annotated_depths["CONTEXT"].fillna('-')
    sample_columns = [x for x in annotated_depths.columns if x not in ["CHROM", "POS", "CONTEXT"] ]
    annotated_depths = annotated_depths[["CHROM", "POS", "CONTEXT"] + sample_columns]
    sample_columns_first = [ str(x).split('.')[0] for x in sample_columns ]
    annotated_depths.columns = ["CHROM", "POS", "CONTEXT"] + sample_columns_first
    annotated_depths.to_csv(output,
                                header=True,
                                index=False,
                                sep="\t")


@click.command()
@click.option('--annotation', type=click.Path(exists=True), help='Input annotation file')
@click.option('--depths', type=click.Path(exists=True), help='Input depths file')
@click.option('--output', type=click.Path(), help='Output annotated depths file')

def main(annotation, depths, output):
    click.echo(f"Annotating depths file...")
    annotate_depths(annotation, depths, output)

if __name__ == '__main__':
    main()

