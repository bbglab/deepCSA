#!/usr/bin/env python

import time

# Third-party imports
import requests
import numpy as np
import pandas as pd
import click


#####
# Define functions
#####
def get_transcript_gene_from_maf(path_maf, consensus_file):
    maf_df = pd.read_table(path_maf)
    genes_in_consensus = pd.read_table(consensus_file)["GENE"].unique()

    # Final filter
    maf_df_f = maf_df[
        (maf_df["TYPE"].isin(["SNV", "INSERTION", "DELETION"]))
        & (maf_df["canonical_SYMBOL"].isin(genes_in_consensus))
        & (maf_df["canonical_Protein_position"] != '-')
        ]

    gene_transcript_pairs = maf_df_f[["canonical_SYMBOL", "canonical_Feature"]].drop_duplicates().reset_index(drop=True)
    gene_transcript_pairs.columns = ["Gene", "Ens_transcript_ID"]

    return gene_transcript_pairs




# Depth
# =====


def get_tr_lookup(transcript_id, max_iter = 20):
    """
    Get exons coord
    """

    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{transcript_id}?expand=1"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    iter_count = 0
    while not r.ok and iter_count < max_iter:
        print("Retrying lookup... (attempt {}/{})".format(iter_count+1, max_iter))
        time.sleep(5)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        iter_count += 1

    return r.json()


def get_cds_coord(transcript_id, len_cds_with_utr, max_iter = 20):
    """
    Get CDS coordinates
    """
    server = "https://rest.ensembl.org"
    ext = f"/map/cds/{transcript_id}/1..{len_cds_with_utr}?"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    iter_count = 0
    while not r.ok and iter_count < max_iter:
        print("Retrying CDS map... (attempt {}/{})".format(iter_count+1, max_iter))
        time.sleep(5)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        iter_count += 1

    return r.json()["mappings"]


def parse_cds_coord(exon):

    strand = exon["strand"]

    if strand == 1:
        start = exon["start"]
        end = exon["end"]
    else:
        start = exon["end"]
        end = exon["start"]

    if "id" in exon:
        exon_id = exon["id"]
        chrom = f'chr{exon["seq_region_name"]}'

        return exon_id, [chrom, start, end, strand]

    else:
        chrom = exon["seq_region_name"]

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

    strand = coord_df["Strand"].unique()[0]

    exons_range = coord_df[["Start", "End"]].values
    exons = np.concatenate([get_dna_exon_pos(exon, strand) for exon in exons_range])
    exons_ix = np.concatenate([get_exon_ix(i, exon, strand) for i, exon in enumerate(exons_range)])
    prot_pos = np.arange(len(exons)) // 3 + 1

    df = pd.DataFrame({"GENE" : coord_df["Gene"].unique()[0],
                        "CHROM" : f'chr{coord_df["Chr"].unique()[0]}',
                        "DNA_POS" : exons,
                        "PROT_POS" : prot_pos,
                        "REVERSE_STRAND" : strand,
                        "EXON_RANK" : exons_ix,
                        "TRANSCRIPT_ID" : coord_df["Ens_transcript_ID"].unique()[0]})

    return df


def get_prot_coverage(dna_prot_df, gene, filter_masked_depth=True):

    gene_dna_prot_df = dna_prot_df[dna_prot_df["GENE"] == gene]
    gene_dna_prot_df = gene_dna_prot_df.dropna(subset=["PROT_POS"])[["PROT_POS", "COVERED", "DEPTH"]].reset_index(drop=True)
    gene_dna_prot_df = gene_dna_prot_df.groupby("PROT_POS").sum().reset_index()
    gene_dna_prot_df.COVERED = (gene_dna_prot_df.COVERED > 0).astype(int)

    return gene_dna_prot_df


def get_exon_coord_wrapper(gene_n_transcript):

    # Init df for coordinates
    coord_df = gene_n_transcript

    # Get coord
    coord_df_lst = []
    exons_coord_df_lst = []
    for gene, transcript in coord_df.values:
        print("Processing gene:", gene)
        coord_lst = []

        # Get the coord of exons with CDS and UTR as well as the length with UTR and exons ID
        exons_lookup = get_tr_lookup(transcript) # We will use this to get Exons ID
        for i, exon in enumerate(exons_lookup["Exon"]):
            exon_id, exons_coord = parse_cds_coord(exon)
            exons_coord_df_lst.append([f"{gene}--exon_{i+1}_{transcript}_{exon_id}"] + exons_coord)

        # Get the CDS coordinates of the exons removing UTRs to map to protein positions
        for i, exon in enumerate(get_cds_coord(transcript, exons_lookup["length"])):
            coord_lst.append((parse_cds_coord(exon) + [i]))

        gene_coord_df = pd.DataFrame(coord_lst, columns = ["Chr", "Start", "End", "Strand", "Exon_rank"])
        gene_coord_df["Gene"] = gene
        gene_coord_df["Ens_transcript_ID"] = transcript
        coord_df_lst.append(gene_coord_df)

    coord_df = pd.concat(coord_df_lst)
    exons_coord_df = pd.DataFrame(exons_coord_df_lst, columns = ["ID", "Chr", "Start", "End", "Strand"])

    return coord_df, exons_coord_df



def dna2prot_depth(gene_list, coord_df, dna_sites, depth_df):
    """
    Get a DNA to protein mapping of all positions in the provided list of genes
    Add as well coverage info & DNA to GENE annotation
    """

    # Map DNA to protein pos, get exons index to protein pos, etc
    dna_prot_df_lst = []
    for gene in gene_list:
        gene_coord_df = coord_df[coord_df["Gene"] == gene]
        dna_prot_df_lst.append(get_dna_map_to_protein(gene_coord_df))
    dna_prot_df = pd.concat(dna_prot_df_lst)

    # dna_prot_df
    #   contains all the protein positions of the genes in the panel
    #     mapped to its corresponding genomic position

    # Merge CDS position with availble sites (not masked) and depth info
    # and any other site that was included in the panel (splicing sites out of the CDS)
    dna_prot_df = dna_sites.merge(dna_prot_df, on=["GENE", "CHROM", "DNA_POS"], how="outer")
    dna_prot_df["COVERED"] = dna_prot_df["CONTEXT"].notnull().astype(int)

    # fill the depth of the regions outside of the consensus panel
    dna_prot_df = dna_prot_df.merge(depth_df.rename(columns={"POS" : "DNA_POS"}),
                                    how="left", on=["CHROM", "DNA_POS"])

    # effectively, when a position is not part of the consensus panel
    # we put the depth of that position to 0
    dna_prot_df.loc[dna_prot_df["COVERED"] == 0, "DEPTH"] = 0

    return dna_prot_df


def get_dna2prot_depth(gene_n_transcript_info, depth_file, consensus_file):
    """
    This function outputs:
    dna_prot_df:

    exons_coord_df:
        df containing the definition of all exons of the genes/transcripts in the panel
        including UTR regions
    """

    consensus_df = pd.read_table(consensus_file)
    depth_df = pd.read_table(depth_file)

    consensus_df = consensus_df.merge(depth_df[["CHROM", "POS", "CONTEXT"]], on = ["CHROM", "POS"], how = 'left')
    consensus_df = consensus_df.rename(columns={"POS" : "DNA_POS"})

    if "DEPTH" not in depth_df:
        depth_df["DEPTH"] = depth_df.drop(columns=["CHROM", "POS", "CONTEXT"]).mean(1)
    depth_df = depth_df[["CHROM", "POS", "DEPTH"]].rename(columns = {"POS" : "DNA_POS"})

    coord_df, exons_coord_df = get_exon_coord_wrapper(gene_n_transcript_info)
    gene_list = list(coord_df["Gene"].unique())
    dna_prot_df = dna2prot_depth(gene_list, coord_df, consensus_df, depth_df)

    dna_prot_df["EXON_ID"] = dna_prot_df.apply(lambda x: find_exon(x, exons_coord_df), axis=1)

    # fix coordinates order in exons_coord_df for BED
    exons_coord_df_final = exons_coord_df.copy()
    exons_coord_df_final.loc[exons_coord_df_final["Strand"] == -1, "Start"] = exons_coord_df.loc[exons_coord_df["Strand"] == -1, "End"]
    exons_coord_df_final.loc[exons_coord_df_final["Strand"] == -1, "End"] = exons_coord_df.loc[exons_coord_df["Strand"] == -1, "Start"]

    return dna_prot_df, exons_coord_df_final


# Utils function to retrieve exon ID from coordinate

def find_exon(x_coord, exon_coord_df):

    dna_pos, chrom, strand = x_coord["DNA_POS"], x_coord["CHROM"], x_coord["REVERSE_STRAND"]

    if strand == -1:
        matches = exon_coord_df[(exon_coord_df['Chr'] == chrom) & (exon_coord_df['End'] <= dna_pos) & (dna_pos <= exon_coord_df['Start'])]

    else:
        matches = exon_coord_df[(exon_coord_df['Chr'] == chrom) & (exon_coord_df['End'] >= dna_pos) & (dna_pos >= exon_coord_df['Start'])]

    return matches['ID'].values[0] if not matches.empty else np.nan




def get_exon_depth_saturation(gene_depth, gene_mut, dna=False):

    # Exon average depth
    gene_depth = gene_depth.copy()
    exon_depth = gene_depth.groupby("EXON_RANK").apply(
        lambda x: 0 if sum(x.COVERED) == 0 else sum(x.DEPTH) / sum(x.COVERED)).reset_index().rename(columns = {0 : "DEPTH"})
    exon_depth["START_PROT_POS"] = gene_depth.groupby("EXON_RANK").apply(lambda x: x.PROT_POS.min())

    # Exon saturation
    if dna:
        exon_depth["COVERED"] = gene_depth.groupby("EXON_RANK").apply(lambda x: sum(x.COVERED)).values
        gene_depth["MUTATED"] = gene_depth.DNA_POS.isin(gene_mut.DNA_POS.unique()).astype(int)
        exon_depth["MUTATED"] = gene_depth.groupby("EXON_RANK").apply(lambda x: sum(x.MUTATED)).values
        exon_depth["SATURATION"] = exon_depth["MUTATED"] / exon_depth["COVERED"]
    else:

        # If an amino acids from the same codon are in two exons, consider the position only in the second one
        # (but look at booth to know if the residue is in the panel and if it is mutated)
        exon_depth_prot = gene_depth.groupby("PROT_POS").apply(lambda x: (x.COVERED.max(), x.EXON_RANK.max())).reset_index()
        exon_depth_prot[["COVERED", "EXON_RANK"]] = pd.DataFrame(exon_depth_prot[0].tolist(), index=exon_depth_prot.index)
        exon_depth_prot = exon_depth_prot.drop(columns=[0])

        exon_depth["COVERED"] = exon_depth_prot.groupby("EXON_RANK").apply(lambda x: sum(x.COVERED)).values
        exon_depth_prot["MUTATED"] = exon_depth_prot["PROT_POS"].isin(gene_mut["Pos"].unique()).astype(int)
        exon_depth["MUTATED"] = exon_depth_prot.groupby("EXON_RANK").apply(lambda x: sum(x.MUTATED)).values
        exon_depth["SATURATION"] = exon_depth["MUTATED"] / exon_depth["COVERED"]

    # check_mutated_masked = exon_depth_prot[(exon_depth_prot["MUTATED"] == 1) & (exon_depth_prot["COVERED"] == 0)]
    # check_mutated_masked["GENE"] = gene_depth.GENE.unique()[0]

    return exon_depth


def get_exon_mid_prot_pos(exon_info, prot_len):

    lst_end_pos = []
    for i in range(len(exon_info)):
        exon = exon_info.iloc[i]
        start_pos = int(exon.START_PROT_POS)
        end_pos = int(exon_info.iloc[i+1].START_PROT_POS) if i < len(exon_info) -1 else prot_len
        lst_end_pos.append(end_pos)
    exon_info["END_PROT_POS"] = lst_end_pos
    exon_info["MID_PROT_POS"] = (exon_info["START_PROT_POS"] + exon_info["END_PROT_POS"]) / 2

    return exon_info


def get_res_coverage(dna_prot_df):

    res_depth = dna_prot_df.dropna(subset=["PROT_POS"])[["PROT_POS", "COVERED", "DEPTH"]].reset_index(drop=True)
    res_depth = res_depth.groupby("PROT_POS").mean().reset_index()

    return res_depth


def add_consecutive_numbers(nums, max_n):

    result = []
    for i in range(len(nums)):
        result.append(nums[i])
        # Check if the current number is the start of a consecutive sequence
        if i < len(nums) - 1 and nums[i] + 1 != nums[i + 1]:
            result.append(nums[i] + 1)

    # Add the last consecutive number after the final element
    result.append(nums[-1] + 1)

    return result


def where_plus(condition):
    """
    Util function to extend the color of mpl filling to the next position.
    """

    ix = np.where(condition)[0]

    if len(ix) > 0:
        ix = add_consecutive_numbers(ix, max_n=len(condition))
        if len(condition) in ix:
            ix.remove(len(condition))
        boolean_vector = np.zeros(len(condition), dtype=bool)
        boolean_vector[ix] = True

        return pd.Series(boolean_vector)

    else:
        return condition


@click.command()
@click.option('--mutations-file', type=click.Path(exists=True), help='Mutations file')
@click.option('--consensus-file', type=click.Path(exists=True), help='Input consensus panel file')
@click.option('--depths-file', type=click.Path(exists=True), help='Input depths of all samples file')
def main(mutations_file, consensus_file, depths_file):
    click.echo("Starting to run plot saturation results...")

    # Count each mutation only ones if it appears in multiple reads
    gene_n_transcript = get_transcript_gene_from_maf(mutations_file, consensus_file)

    exons_depth, exons_coord_id = get_dna2prot_depth(gene_n_transcript, depths_file, consensus_file)
    print("Exons coordinates and depth computed")
    exons_depth.to_csv("depths_per_position_exon_gene.tsv", header = True, index = False, sep = '\t')

    exons_coordinates_bed_like = exons_coord_id[['Chr', 'Start', 'End', 'ID']]
    exons_coordinates_bed_like.to_csv("panel_exons.bed4.bed", header = False, index = False, sep = '\t')


if __name__ == '__main__':
    main()
