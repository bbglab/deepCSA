#!/usr/bin/env python

"""
Merge annotation depths - Depth Annotation and Flagging Script

This script processes a depths annotation file to merge and annotate depth information across samples.

Command-line Arguments
----------------------

Authors
-------
Author  : Ferriol Calvet (@FerriolCalvet)
Email   : ferriol.calvet@irbbarcelona.org

Contributors
------------
- Federica Brando - @FedericaBrando (federica.brando@irbbarcelona.org)
- Marta Huertas - @m-huertasp (marta.huertas@irbbarcelona.org)
"""

import click
import json
import subprocess
import polars as pl
import logging

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p"
)
LOG = logging.getLogger("merge_annotation_depths")

# Globals
COLS = ["CHROM", "POS", "CONTEXT"]

# Functions
def preprocess(annotation_file: str, depths_file: str, input_csv: str) -> tuple[pl.DataFrame, list]:
    """
    Merge annotation and depth files.

    Parameters
    ----------
    annotation_file : str
        Path to the annotation file.
    depths_file : str
        Path to the depths file.
    input_csv : str
        Path to the input CSV file.

    Returns
    -------
    pl.DataFrame
        Merged DataFrame with annotations and depths.
    list
        List of sample column names without suffixes.
    """
    LOG.info("Preprocessing annotation and depths files...")
    _depths = pl.read_csv(depths_file, separator="\t")
    _annots = pl.read_csv(annotation_file, separator="\t")
    _input = pl.read_csv(input_csv, separator=",", infer_schema_length=0)

    if "bam" in _input.columns:
        bam2sample_dict = dict(
            _input.select(
                pl.col("bam").str.split("/").list.last(),
                pl.col("sample")
            ).rows()
        )
    else:
        bam2sample_dict = dict(zip(_input["sample"], _input["sample"]))

    # Drop CONTEXT column if it exists
    _depths = _depths.drop("CONTEXT", strict=False)
    annot_depth = _depths.join(_annots, on=["CHROM", "POS"], how='left').fill_null(pl.lit('-'))

    sample_columns = [col for col in annot_depth.columns if col not in COLS]

    # Ensure all columns but CHROM and CONTEXT are numeric
    annot_depth = annot_depth.with_columns([
        pl.col(col).cast(pl.Int64) if col not in ["CHROM", "CONTEXT"] else pl.col(col)
        for col in annot_depth.columns
    ])

    LOG.info("Samples: %s", list(bam2sample_dict.values()))  # List of samples without the .*.bam suffix

    # Place COLS=[CHROM, POS, CONTEXT] columns at the beginning then the rest of the columns
    return annot_depth.select(COLS + sample_columns).rename(bam2sample_dict), list(bam2sample_dict.values())

def apply_mask_matrix(annotated_depths: pl.DataFrame, mask_matrix_file: str) -> pl.DataFrame:
    """
    Apply position per sample mask matrix to depths.
    
    Sets depth = 0 where mask value is 0 (position should be masked).
    Keeps depth as-is where mask value is 1 (position should be kept).
    
    Parameters
    ----------
    annotated_depths : pl.DataFrame
        Depths dataframe with CHROM, POS, CONTEXT, and sample columns
    mask_matrix_file : str
        Path to mask matrix file (TSV, gzipped)
    
    Returns
    -------
    pl.DataFrame
        Depths with masked positions (where mask=0) set to 0
    """
    LOG.info("Loading mask matrix...")
    mask_df = pl.read_csv(mask_matrix_file, separator="\t")
    
    if mask_df.is_empty():
        LOG.info("Mask matrix is empty, no masking applied")
        return annotated_depths
    
    LOG.info(f"Mask matrix loaded: {len(mask_df)} positions")
    
    
    # Obtain sample columns present in both annotated_depths and mask_df
    sample_cols = [col for col in mask_df.columns if col in annotated_depths.columns and col not in ["CHROM", "POS"]]
    
    LOG.info(f"Applying mask to {len(sample_cols)} samples")
    
    # Join on CHROM and POS. For each sample, multiply depth by mask (0 or 1). If there is no mask, keep original depth (multiply by 1).
    annotated_depths_filtered = (
        annotated_depths.join(
            mask_df.select(["CHROM", "POS", *sample_cols]), 
            on=["CHROM", "POS"], 
            how="left", 
            suffix="_mask"
        )
        .with_columns([
            # If the mask is missing (null), we fill with 1 to keep original depth
            (pl.col(s) * pl.col(f"{s}_mask").fill_null(1)).alias(s)
            for s in sample_cols
        ])
        .drop([f"{col}_mask" for col in sample_cols])
    )
    
    return annotated_depths_filtered

def output_annotate_dephts(annotated_depths: pl.DataFrame, json_f: str, samples: list):
    """
    Output annotated depths file

    Parameters
    ----------
    annotated_depths : pl.DataFrame
        Annotated depths dataframe
    json_f : str
        JSON file with groups information
    samples : list
        List of samples
    """
    LOG.info("Outputting annotated depths file...")

    # Output annotated depths file for all samples
    # Polars functionality to write gzipped files is currently not working, so we write the file and then gzip it with subprocess
    temp_file = "all_samples_indv.depths.tsv"
    annotated_depths.write_csv(temp_file, include_header=True, separator="\t")
    subprocess.run(["gzip", "-f", temp_file], check=True)

    try:
        with open(json_f, 'r') as file:
            groups_info = json.load(file)
            LOG.info("JSON file found. Outputting annotated depths file for each group.")

        # Per group
        for group_name, samples in groups_info.items():
            annotated_depths = annotated_depths.with_columns(
                pl.sum_horizontal([pl.col(sample) for sample in samples]).alias(group_name)
            )
            temp_file = f"{group_name}.depths.annotated.tsv"
            annotated_depths.select(COLS + [group_name]).write_csv(temp_file, separator="\t", include_header=True)
            subprocess.run(["gzip", "-f", temp_file], check=True)
    except (TypeError, FileNotFoundError):
        LOG.warning("JSON file not found. Outputting annotated depths file for each sample.")
        
        # Per sample
        for sample in samples:
            temp_file = f"{sample}.depths.annotated.tsv"
            annotated_depths.select(COLS + [str(sample)]).write_csv(temp_file, separator="\t", include_header=True)
            subprocess.run(["gzip", "-f", temp_file], check=True)

        # All the samples
        annotated_depths = annotated_depths.with_columns(
            pl.sum_horizontal(pl.exclude(COLS)).alias("all_samples"))
        
        temp_file = "all_samples.depths.annotated.tsv"
        annotated_depths.select(COLS + ["all_samples"]).write_csv(temp_file, separator="\t", include_header=True)
        subprocess.run(["gzip", "-f", temp_file], check=True)

@click.command()
@click.option('--annotation', type=click.Path(exists=True), help='Input annotation file')
@click.option('--depths', type=click.Path(exists=True), help='Input depths file')
@click.option('--input-csv', type=click.Path(exists=True), help='Input CSV file')
@click.option('--json_file', type=click.Path(exists=True), help='JSON groups file')
@click.option('--mask-matrix', type=click.Path(exists=True), help='Position per sample mask matrix file (1=keep, 0=mask)')
def main(annotation, depths, input_csv, json_file, mask_matrix):
    LOG.info("Annotating depths file...")

    # Preprocess annotation and depths files
    annotated_depths, samples  = preprocess(annotation, depths, input_csv)

    # Apply position masking
    annotated_depths = apply_mask_matrix(annotated_depths, mask_matrix)

    output_annotate_dephts(annotated_depths, json_file, samples)

    LOG.info("Done!")

if __name__ == '__main__':
    main()

