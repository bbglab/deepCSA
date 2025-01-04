#!/usr/local/bin/python

# Standard library imports
import os
import sys
import time

# Third-party imports
import requests
import numpy as np
import pandas as pd
import click



# ---------------------- Data Processing Functions ---------------------- #

def get_normal_maf(path_maf):
    """
    Process MAF file to extract relevant mutations and columns

    Args:
        path_maf (str): Path to MAF file.
        only_protein_pos (bool): Filter for protein-level mutations if True.

    Returns:
        pd.DataFrame: Filtered and renamed MAF DataFrame.
    """
    maf_df = pd.read_table(path_maf)
    print("Initial shape of MAF DataFrame:", maf_df.shape)

    maf_df_filtered = maf_df.loc[
        (~maf_df["FILTER.not_in_panel"]) &
        (maf_df["TYPE"].isin(["SNV", "INSERTION", "DELETION"]))
    ].reset_index(drop=True)

    # Select and rename columns
    cols = [
        "canonical_SYMBOL", "canonical_Feature"
    ]
    maf_df_filtered = maf_df_filtered[cols].rename(columns={
        "canonical_SYMBOL": "GENE",
        "canonical_Feature": "Ens_transcript_ID"
    })
    return maf_df_filtered




def get_cds_len_with_utr(transcript_id):
    """
    Retrieve the length of a transcript (with UTR).

    Args:
        transcript_id (str): Ensembl transcript ID.

    Returns:
        dict: Transcript information with length.
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{transcript_id}?expand=1"

    while True:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if r.ok:
            break
        print("Retrying transcript lookup...")
        time.sleep(5)
    return r.json()


def get_cds_coord(transcript_id):
    """
    Get CDS coordinates for a given transcript.

    Args:
        transcript_id (str): Ensembl transcript ID.

    Returns:
        list: List of CDS mappings.
    """
    length = get_cds_len_with_utr(transcript_id)["length"]
    server = "https://rest.ensembl.org"
    ext = f"/map/cds/{transcript_id}/1..{length}?"

    while True:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if r.ok:
            break
        print("Retrying CDS mapping...")
        time.sleep(5)
    return r.json()["mappings"]


def parse_cds_coord(exon):
    """
    Parse individual CDS coordinates from Ensembl API.

    Args:
        exon (dict): Exon data.

    Returns:
        list: Chromosome, start, end, and strand.
    """
    strand = exon["strand"]
    chrom = exon["seq_region_name"]
    start, end = (exon["start"], exon["end"]) if strand == 1 else (exon["end"], exon["start"])
    return [chrom, start, end, strand]


def get_dna_map_to_protein(coord_df):
    """
    Map DNA positions to protein positions.

    Args:
        coord_df (pd.DataFrame): Coordinate DataFrame.

    Returns:
        pd.DataFrame: Mapping of DNA positions to protein positions.
    """
    strand = coord_df.Strand.unique()[0]
    exons_range = coord_df[["Start", "End"]].values

    exons = np.concatenate([get_dna_exon_pos(exon, strand) for exon in exons_range])
    exons_ix = np.concatenate([get_exon_ix(i, exon, strand) for i, exon in enumerate(exons_range)])

    df = pd.DataFrame({
        "GENE": coord_df["GENE"].unique()[0],
        "CHR": f'chr{coord_df["Chr"].unique()[0]}',
        "DNA_POS": exons,
        "PROT_POS": np.arange(len(exons)) // 3 + 1,
        "REVERSE_STRAND": strand,
        "EXON": exons_ix
    })
    return df


def get_dna_exon_pos(exon_range, strand):
    """
    Get DNA positions for a given exon range.

    Args:
        exon_range (list): Start and end of exon.
        strand (int): Strand orientation (1 or -1).

    Returns:
        np.ndarray: Array of DNA positions.
    """
    positions = np.arange(exon_range[1], exon_range[0] + 1) if strand == -1 else np.arange(exon_range[0], exon_range[1] + 1)
    return positions[::-1] if strand == -1 else positions


def get_exon_ix(i, exon_range, strand):
    """
    Assign exon index for DNA positions.

    Args:
        i (int): Exon index.
        exon_range (list): Start and end of exon.
        strand (int): Strand orientation.

    Returns:
        np.ndarray: Repeated exon index.
    """
    length = len(get_dna_exon_pos(exon_range, strand))
    return np.repeat(i, length)


# ---------------------- DNA and Protein Mapping ---------------------- #

def generate_dna_protein_mapping(maf_file, consensus_file_unique, output_file):
    """
    Main function to process DNA to protein mapping.
    """
    # Load and prepare data
    dna_sites = pd.read_table(consensus_file_unique).rename(columns={"POS": "DNA_POS", "CHROM": "CHR"})
    normal_maf_df = get_normal_maf(maf_file)
    coord_df = normal_maf_df[normal_maf_df["GENE"] != '-'][["GENE", "Ens_transcript_ID"]].drop_duplicates().reset_index(drop=True)

    # Retrieve exon coordinates
    coord_df_lst = []
    for gene, transcript in coord_df.values:
        print("Processing gene:", gene)
        exons = [parse_cds_coord(exon) + [i] for i, exon in enumerate(get_cds_coord(transcript))]
        gene_df = pd.DataFrame(exons, columns=["Chr", "Start", "End", "Strand", "Exon"])
        gene_df["GENE"] = gene
        gene_df["Ens_transcript_ID"] = transcript
        coord_df_lst.append(gene_df)

    coord_df = pd.concat(coord_df_lst)

    # Map DNA to protein positions
    dna_prot_df_lst = []
    for gene in coord_df["GENE"].unique():
        print("Mapping DNA to protein for gene:", gene)
        gene_coord_df = coord_df[coord_df["GENE"] == gene]
        dna_prot_df_lst.append(get_dna_map_to_protein(gene_coord_df))

    dna_prot_df = pd.concat(dna_prot_df_lst)

    # Merge data and save output
    merged_df = dna_sites.merge(dna_prot_df, on=["GENE", "CHR", "DNA_POS"], how="outer")
    merged_df["COVERED"] = merged_df["CONTEXT"].notnull().astype(int)
    merged_df.to_csv(output_file, sep="\t", index=False)



# ---------------------- CLI Interface ---------------------- #

@click.command()
@click.option('--maf', type=click.Path(exists=True), help='Input mutations file')
@click.option('--consensus-file', type=click.Path(exists=True), help='Input consensus panel file')
@click.option('--output', type=click.Path(), help='Output file')

def main(maf, consensus_file, output):
    click.echo("Retrieving the DNA 2 protein mapping...")
    generate_dna_protein_mapping(maf, consensus_file, output)



if __name__ == '__main__':
    main()





"""


# Covered DNA and DNA to protein mapping
# Use real coordinates from gene transcripts

# Get exons coord
# ---------------
def get_cds_len_with_utr(transcript_id):

    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{transcript_id}?expand=1"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    while not r.ok:
        print("Retrying lookup..")
        time.sleep(5)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    return r.json()

def get_cds_coord(transcript_id):

    # The len I need is actually without UTR but it works anyway
    len_cds_with_utr = get_cds_len_with_utr(transcript_id)["length"]
    server = "https://rest.ensembl.org"
    ext = f"/map/cds/{transcript_id}/1..{len_cds_with_utr}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    while not r.ok:
        print("Retrying DCS map..")
        time.sleep(5)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    return r.json()["mappings"]


def parse_cds_coord(exon):

    strand = exon["strand"]
    chrom = exon["seq_region_name"]
    if strand == 1:
        start = exon["start"]
        end = exon["end"]
    else:
        start = exon["end"]
        end = exon["start"]

    return [chrom, start, end, strand]


# Get Exon coord to protein pos
# -----------------------------
def get_dna_exon_pos(exon_range, strand):
    if strand == -1:
         return np.arange(exon_range[1], exon_range[0] + 1)[::-1]
    else:
        return np.arange(exon_range[0], exon_range[1] + 1)

def get_exon_ix(i, exon_range, strand):

    len_exon = len(get_dna_exon_pos(exon_range, strand))

    return np.repeat(i, len_exon)


def get_dna_map_to_protein(coord_df):

    strand = coord_df.Strand.unique()[0]
    exons_range = coord_df[["Start", "End"]].values
    exons = np.concatenate([get_dna_exon_pos(exon, strand) for exon in exons_range])
    exons_ix = np.concatenate([get_exon_ix(i, exon, strand) for i, exon in enumerate(exons_range)])

    prot_pos = np.arange(len(exons)) // 3 + 1

    df = pd.DataFrame({"GENE" : coord_df.Gene.unique()[0],
                       "CHR" : f'chr{coord_df.Chr.unique()[0]}',
                       "DNA_POS" : exons,
                       "PROT_POS" : np.arange(len(exons)) // 3 + 1,
                       "REVERSE_STRAND" : strand,
                       "EXON" : exons_ix})

    return df




# Count each mutation only ones if it appears in multiple reads
normal_maf_df = get_normal_maf(path_normal_maf)


# To generate DNA sites
dna_sites = pd.read_table(consensus_file_unique)
dna_sites = dna_sites.rename(columns={"POS" : "DNA_POS", "CHROM" : "CHR"})

# Get exons coordinates

# Init df for coordinates
coord_df = normal_maf_df[["Gene", "Ens_transcript_ID"]].drop_duplicates().reset_index(drop=True)
# Get coord
coord_df_lst = []
for gene, transcript in coord_df.values:
    print(gene)
    coord_lst = []
    for i, exon in enumerate(get_cds_coord(transcript)):
        coord_lst.append((parse_cds_coord(exon) + [i]))

    gene_coord_df = pd.DataFrame(coord_lst, columns = ["Chr", "Start", "End", "Strand", "Exon"])
    gene_coord_df["Gene"] = gene
    gene_coord_df["Ens_transcript_ID"] = transcript
    coord_df_lst.append(gene_coord_df)

coord_df = pd.concat(coord_df_lst)



# Get a DNA to protein mapping & DNA to GENE annotation

# Map DNA to protein pos, get exons index to protein pos, etc
dna_prot_df_lst = []
for gene in gene_coord_df["Gene"].unique():
    print(gene)
    gene_coord_df = coord_df[coord_df["Gene"] == gene]
    dna_prot_df_lst.append(get_dna_map_to_protein(gene_coord_df))
dna_prot_df = pd.concat(dna_prot_df_lst)

# Merge CDS position with availble sites (not masked)
# and any other site that was included in the panel (splicing sites out of the CDS)
dna_prot_df = dna_sites.merge(dna_prot_df, on=["GENE", "CHR", "DNA_POS"], how="outer")
dna_prot_df["COVERED"] = dna_prot_df["CONTEXT"].notnull().astype(int)

dna_prot_df.to_csv(output_file, sep="\t", index=False)


"""
