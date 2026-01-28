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
def preprocess(annotation_file: str, depths_file: str) -> tuple[pd.DataFrame, list]:
    """
    Merge annotation and depth files.

    Parameters
    ----------
    annotation_file : str
        Path to the annotation file.
    depths_file : str
        Path to the depths file.

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with annotations and depths.
    list
        List of sample column names without suffixes.
    """
    LOG.info("Preprocessing annotation and depths files...")
    _depths = pd.read_csv(depths_file, sep = "\t", header = 0)
    _annots = pd.read_csv(annotation_file, sep = "\t", header = 0)

    annot_depth = _depths.merge(_annots, on = ["CHROM", "POS"], how = 'left').fillna(value={'CONTEXT': '-'})

    sample_columns = annot_depth.columns.difference(COLS).tolist()

    rename_map = {col: col.split('.')[0] for col in sample_columns} # Dict to remove the .*.bam suffix

    LOG.info("Samples: %s", list(rename_map.values())) # List of samples without the .*.bam suffix

    # Place COLS=[CHROM, POS, CONTEXT] columns at the beginning then the rest of the columns
    return annot_depth[COLS + sample_columns].rename(columns=rename_map), list(rename_map.values())

def apply_mask_matrix(annotated_depths: pd.DataFrame, mask_matrix_file: str) -> pd.DataFrame:
    """
    Apply position per sample mask matrix to depths.
    
    Sets depth = 0 where mask value is 0 (position should be masked).
    Keeps depth as-is where mask value is 1 (position should be kept).
    
    Parameters
    ----------
    annotated_depths : pd.DataFrame
        Depths dataframe with CHROM, POS, CONTEXT, and sample columns
    mask_matrix_file : str
        Path to mask matrix file (TSV, gzipped)
    
    Returns
    -------
    pd.DataFrame
        Depths with masked positions (where mask=0) set to 0
    """
    LOG.info("Loading mask matrix...")
    mask_df = pd.read_csv(mask_matrix_file, sep="\t", compression='gzip')
    
    if mask_df.empty:
        LOG.info("Mask matrix is empty, no masking applied")
        return annotated_depths
    
    LOG.info(f"Mask matrix loaded: {len(mask_df)} positions")
    
    # Create a copy to avoid modifying the original
    result = annotated_depths.copy()
    
    # Set CHROM and POS as index for alignment
    result.set_index(['CHROM', 'POS'], inplace=True)
    mask_df.set_index(['CHROM', 'POS'], inplace=True)
    
    # Find common samples between depths and mask
    sample_cols = [col for col in mask_df.columns if col in result.columns]
    
    LOG.info(f"Applying mask to {len(sample_cols)} samples")
    
    # Multiply depths by mask (0 or 1) -> only positions with mask=1 will retain their depth
    common_idx = result.index.intersection(mask_df.index)
    result.loc[common_idx, sample_cols] = result.loc[common_idx, sample_cols] * mask_df.loc[common_idx, sample_cols]
    
    # Indicate number of positions masked
    n_positions_masked = len(common_idx)
    LOG.info(f"Masked {n_positions_masked} positions across {len(sample_cols)} samples")
    
    # Reset index to restore original structure
    annotated_depths_filtered = result.reset_index()
    
    return annotated_depths_filtered

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
        # If there are defined groups in the JSON file, output annotated depths file for each group
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
@click.option('--mask-matrix', type=click.Path(exists=True), help='Position per sample mask matrix file (1=keep, 0=mask)')
def main(annotation, depths, json_file, mask_matrix):
    LOG.info("Annotating depths file...")

    # Preprocess annotation and depths files
    annotated_depths, samples  = preprocess(annotation, depths)

    # Apply position masking
    annotated_depths = apply_mask_matrix(annotated_depths, mask_matrix)

    output_annotate_dephts(annotated_depths, json_file, samples)

    LOG.info("Done!")

if __name__ == '__main__':
    main()

