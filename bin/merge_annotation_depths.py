#!/usr/local/bin/python


import click
import json
import pandas as pd
import logging

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("merge_annotation_depths")

# Globals
COLS = ["CHROM", "POS", "CONTEXT"]

# Functions
def preprocess(annots, depths):
    """
    Merge annotation and depths files.
    """
    LOG.info("Preprocessing annotation and depths files...")
    _depths = pd.read_csv(depths, sep = "\t", header = 0)
    _annots = pd.read_csv(annots, sep = "\t", header = 0)

    annot_depth =  _depths.merge(_annots, on = ["CHROM", "POS"], how = 'left').fillna(value={'CONTEXT': '-'})

    sample_columns = annot_depth.columns.difference(COLS).tolist()
    rename_map = {col: col.split('.')[0] for col in sample_columns} # Dict to remove the .*.bam suffix

    LOG.info("Samples: %s", list(rename_map.values())) # List of samples without the .*.bam suffix

    # Place COLS=[CHROM, POS, CONTEXT] columns at the beginning then the rest of the columns
    return annot_depth[COLS + sample_columns].rename(columns=rename_map), list(rename_map.values())


def mask_panel_regions(annotated_depths, regions_to_filter):
    """
    Mask regions flagged in the BED file by assigning depth 0 to them.
    """
    LOG.info("Masking regions --> Assign depth 0 to regions flagged in the BED file")
    regions_to_mask = pd.read_csv(regions_to_filter, usecols=[0,1,3], sep = "\t", names=["CHROM", "POS", "FILTERS"])

    LOG.debug("Regions to mask: %s", regions_to_mask.drop_duplicates().shape[0])
    mask = (annotated_depths.set_index(["CHROM", "POS"]).index.isin(regions_to_mask.set_index(["CHROM", "POS"]).index))
    annotated_depths.loc[mask, annotated_depths.columns.difference(COLS)] = 0
    
    LOG.info("Regions masked in depths file: %s", mask.sum())
    

    return annotated_depths


def output_annotate_dephts(annotated_depths, json_f, samples):
    """
    Output annotated depths file

    Parameters
    ----------
    annotated_depths : pd.DataFrame
        Annotated depths dataframe
    json_f : str
        JSON file with groups information
    samples : list
        List of samples
    """
    LOG.info("Outputting annotated depths file...")

    # Output annotated depths file for all samples
    annotated_depths.to_csv("all_samples_indv.depths.tsv.gz",
                                header=True,
                                index=False,
                                sep="\t")

    try:
        with open(json_f, 'r') as file:
            groups_info = json.load(file)
            LOG.info("JSON file found. Outputting annotated depths file for each group.")
        
        for group_name, samples in groups_info.items():
            annotated_depths[group_name] = annotated_depths.loc[:,samples].sum(axis=1)
            annotated_depths[COLS + [group_name]].to_csv(f"{group_name}.depths.annotated.tsv.gz",
                                                                                sep = "\t",
                                                                                header = True,
                                                                                index = False)   
    except (TypeError, FileNotFoundError):
        LOG.warning("JSON file not found. Outputting annotated depths file for each sample.")
        for sample in samples:
            annotated_depths[COLS + [str(sample)]].to_csv(f"{sample}.depths.annotated.tsv.gz",
                                                                                sep = "\t",
                                                                                header = True,
                                                                                index = False)

        annotated_depths["all_samples"] = annotated_depths.iloc[:,3:].sum(axis=1)
        annotated_depths[COLS + ["all_samples"]].to_csv("all_samples.depths.annotated.tsv.gz",
                                                                            sep = "\t",
                                                                            header = True,
                                                                            index = False)



@click.command()
@click.option('--annotation', type=click.Path(exists=True), help='Input annotation file')
@click.option('--depths', type=click.Path(exists=True), help='Input depths file')
@click.option('--json_file', type=click.Path(exists=True), help='JSON groups file')
@click.option('--regions-to-filter', type=click.Path(), help='BED file with regions to filter')
# @click.option('--output', type=click.Path(), help='Output annotated depths file')
def main(annotation, depths, json_file, regions_to_filter):
    LOG.info("Annotating depths file...")

    # Preprocess annotation and depths files
    annotated_depths, samples  = preprocess(annotation, depths)

    print(annotated_depths.head())

    # Mask regions flagged in the BED file
    annotated_depths = mask_panel_regions(annotated_depths, regions_to_filter)

    output_annotate_dephts(annotated_depths, json_file, samples)

    LOG.info("Done!")

if __name__ == '__main__':
    main()

