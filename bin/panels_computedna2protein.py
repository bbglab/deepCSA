#!/usr/bin/env python

"""
Panels Compute DNA to Protein - MAF processing to protein positions and 
coverage computation.

This script processes the MAF file to retrieve the transcript-gene pairs,
then retrieves exon, and CDS coordinates and maps DNA positions to protein positions.
It also computes the coverage of each position and generates plots of coverage
per gene.

Arguments:
----------------------
--mutations-file : str 
    Path to the MAF file containing mutations.
--consensus-file : str
    Path to the consensus panel file.
--depths-file : str
    Path to the file containing depth information for all samples.

Authors
----------------------
Author  : Stefano Pellegrini (@St3451)
Email   : stefano.pellegrini@irbbarcelona.org

Contributors
----------------------
- Ferriol Calvet - @FerriolCalvet (ferriol.calvet@irbbarcelona.org)
- Marta Huertas - @m-huertasp (marta.huertas@irbbarcelona.org)
"""

import click
import logging
import requests
import numpy as np
import pandas as pd
import polars as pl
import io
import gzip
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Logging
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s - %(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p"
)
logging.getLogger('matplotlib').setLevel(logging.WARNING) # Avoid too much logging

LOG = logging.getLogger("DNA2protein")


#  Data parsing and processing functions
# ----------------------------------------------------------
def get_transcript_gene_from_maf(path_maf, consensus_file):
    """
    Process MAF file to retrieve gene-transcript pairs for the genes in the consensus panel.

    Parameters
    ----------
    path_maf : str
        Path to the MAF file containing mutations.
    consensus_file : str
        Path to the consensus panel file.

    Returns
    ------------
    gene_transcript_pairs : pandas.DataFrame
        DataFrame with gene-transcript pairs.
    """
    maf_df = pd.read_table(path_maf)
    genes_in_consensus = pd.read_table(consensus_file)["GENE"].unique()

    # Filter MAF (keep only positions with canonical protein position)
    maf_df_f = maf_df[
        (maf_df["TYPE"].isin(["SNV", "INSERTION", "DELETION"]))
        & (maf_df["canonical_SYMBOL"].isin(genes_in_consensus))
        & (maf_df["canonical_Protein_position"] != '-')
        ]

    gene_transcript_pairs = maf_df_f[["canonical_SYMBOL", "canonical_Feature"]].drop_duplicates().reset_index(drop=True)
    gene_transcript_pairs.columns = ["Gene", "Ens_transcript_ID"]

    LOG.info(f"Retrieved {len(gene_transcript_pairs)} gene-transcript pairs from MAF file.")
    return gene_transcript_pairs

# Generator to filter lines before they reach a dataframe
def get_gff_to_generator(release: int = 111, species: str = "homo_sapiens", genome: str = "GRCh38", gff3_file: str | None = None):
    """
    Get GFF file from a local path or Ensembl FTP and filter it on the fly to
    keep only exon and CDS lines.

    Parameters
    ------------
    release : int
        The release number of the Ensembl GFF file to retrieve (default is 111).
    gff3_file : str | None
        Optional path to a local GFF3 or GFF3.GZ file. If provided, this file
        is used instead of downloading from Ensembl.

    Returns
    ------------
    generator
        A generator that yields lines from the GFF file that correspond to exon and CDS features.
    """
    if gff3_file:
        LOG.info(f"Using local GFF3 file: {gff3_file}")

        open_func = gzip.open if gff3_file.endswith(".gz") else open
        mode = "rt"
        try:
            with open_func(gff3_file, mode, encoding="utf-8") as handle:
                for line_str in handle:
                    # Skip comments
                    if line_str.startswith('#'):
                        continue

                    # Only "yield" lines that are exon or CDS
                    parts = line_str.split('\t')
                    if len(parts) > 2 and parts[2] in ["exon", "CDS"]:
                        yield line_str
        except OSError as e:
            error = RuntimeError(f"Failed to read local GFF3 file: {gff3_file}")
            LOG.error(error)
            raise error from e
        return

    url = f"https://ftp.ensembl.org/pub/release-{release}/gff3/{species}/{species.capitalize()}.{genome}.{release}.gff3.gz"

    # Open request
    try:
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()
    except requests.exceptions.Timeout:
        error = RuntimeError(f"Request timed out while fetching Ensembl GFF3 (release {release}) from: {url}")
        LOG.error(error)
        raise error
    except requests.exceptions.HTTPError as e:
        error = RuntimeError(f"HTTP error {response.status_code} while fetching Ensembl GFF3 (release {release}) from: {url}")
        LOG.error(error)
        raise error
    except requests.exceptions.RequestException as e:
        error = RuntimeError(f"Failed to fetch Ensembl GFF3 (release {release}) from: {url}")
        LOG.error(error)
        raise error

    # Decompress the stream on the fly
    with gzip.GzipFile(fileobj=response.raw) as gz:
        for line in gz:
            # Decode bytes to string
            line_str = line.decode('utf-8')

            # Skip comments
            if line_str.startswith('#'):
                continue

            # Only "yield" lines that are exon or CDS
            parts = line_str.split('\t')
            if len(parts) > 2 and parts[2] in ["exon", "CDS"]:
                yield line_str

def gff_to_filtered_df(gene_n_transcript: pd.DataFrame, release: int, species: str, genome: str, gff3_file: str | None = None) -> pd.DataFrame:
    """
    Transforms the yields from get_gff_to_generator into a filtered DataFrame. The reading and filtering is done with polars
    to improve efficiency.

    Parameters
    ------------
    gene_n_transcript : pd.DataFrame
        A DataFrame containing gene and transcript information.
    release : int
        The Ensembl release number to use for GFF file retrieval (default is 111).
    species : str
        The species for which to retrieve GFF data.
    genome : str
        The genome assembly for which to retrieve GFF data.
    gff3_file : str | None
        Optional path to a local GFF3 or GFF3.GZ file. If provided, this file
        is used instead of downloading from Ensembl.

    Returns
    ------------
    pd.DataFrame
        A filtered DataFrame of GFF lines for the specified genes and release.
    """
    # Join the generator into a single buffer for Polars to read
    filtered_buffer = io.StringIO("".join(get_gff_to_generator(release=release, species=species, genome=genome, gff3_file=gff3_file)))

    # Read generator with polars, generate transcript_id column and filter by the genes in the panel
    df = (
        pl.read_csv(
            filtered_buffer,
            separator="\t",
            has_header=False,
            new_columns=["chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"],
            schema_overrides={"start": pl.Int64, "end": pl.Int64, "chr": pl.Utf8, "feature": pl.Utf8, "attributes": pl.Utf8}
        )
        .with_columns([
            pl.col("attributes").str.extract(r"transcript:(ENST\d+)", 1).alias("transcript_id")
        ])
        .filter(
            [pl.col("transcript_id").is_in(gene_n_transcript["Ens_transcript_ID"].to_list())]
        )
    )
    
    LOG.info(f"Filtered GFF data to {df.shape[0]} rows for the specified genes and release.")
    # Transform to pandas for easier manipulation later on
    return df.to_pandas()

def _parse_strand_coords(exon) -> tuple[str, int, int, int]:
    """Extract chromosome, strand-ordered start/end, and strand integer from a GFF row."""
    strand = 1 if exon.strand == "+" else -1
    if strand == 1:
        start, end = exon.start, exon.end
    else:
        start, end = exon.end, exon.start
    return exon.chr, start, end, strand

def parse_exon_coord(exon: pd.Series) -> tuple[str, list]:
    """
    Parses the coordinates of an exon row from the GFF DataFrame and returns the exon ID and its coordinates in a list format.

    Parameters
    ------------
    exon : pandas.Series
        A row from the GFF DataFrame corresponding to an exon feature.

    Returns
    ------------
    exon_id : str
        The ID of the exon extracted from the attributes column.
    exons_coord : list
        A list containing the chromosome, start, end, and strand information of the exon.
    """
    chrom, start, end, strand = _parse_strand_coords(exon)
    
    # Add "chr" prefix to chromosome
    chrom = f'chr{exon.chr}'

    # Extract exon_id using regex
    match = re.search(r"exon_id=([^;]+)", exon.attributes)
    
    # Extract exon_id using regex
    if match:
        exon_id = match.group(1)
        
    else:
        exon_id = ""
        LOG.warning(f"Could not extract exon_id from attributes: {exon.transcript_id} - {exon.attributes}")

    return exon_id, [chrom, start, end, strand]

def get_exon_coord_wrapper(gene_n_transcript: pd.DataFrame, gff_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]: 
    """
    Wrapper function to retrieve exon and CDS coordinates for the genes and transcripts in the panel.

    Parameters
    ------------
    gene_n_transcript : pandas.DataFrame
        A DataFrame containing gene and transcript information for the panel.
    gff_df : pandas.DataFrame
        A DataFrame containing the GFF data filtered for exon and CDS features.

    Returns
    ------------
    final_coord_df : pandas.DataFrame
        A DataFrame containing the CDS coordinates for each gene-transcript pair, along with protein position mapping.
    final_exons_df : pandas.DataFrame
        A DataFrame containing the exon coordinates (including UTRs) for each gene-transcript pair.
    """
    coord_df_lst = []
    exons_coord_df_lst = []

    for gene, transcript in gene_n_transcript.values:
        
        # Get all features for this transcript once to avoid double filtering
        transcript_data = gff_df[gff_df["transcript_id"] == transcript]
        if transcript_data.empty:
            continue

        # Determine strand from the first available entry
        strand = transcript_data["strand"].iloc[0]
        is_positive = (strand == "+")

        # --- Handle exons (with UTRs) ---
        exons_lookup = transcript_data[transcript_data["feature"] == "exon"]
        # Sort based on strand: '+' is low-to-high, '-' is high-to-low
        exons_lookup = exons_lookup.sort_values("start", ascending=is_positive)

        for i, exon in enumerate(exons_lookup.itertuples(index=False)):
            exon_id, exons_coord = parse_exon_coord(exon)
            exons_coord_df_lst.append([f"{gene}--exon_{i+1}_{transcript}_{exon_id}"] + exons_coord)

        # --- Handle CDS ---
        cds_lookup = transcript_data[transcript_data["feature"] == "CDS"]
        cds_lookup = cds_lookup.sort_values("start", ascending=is_positive)

        coord_lst = []
        for i, exon in enumerate(cds_lookup.itertuples(index=False)):
            exons_coord = list(_parse_strand_coords(exon))
            # Include the biological rank (i)
            coord_lst.append(exons_coord + [i])

        if coord_lst:
            gene_coord_df = pd.DataFrame(coord_lst, columns=["Chr", "Start", "End", "Strand", "Exon_rank"])
            gene_coord_df["Gene"] = gene
            gene_coord_df["Ens_transcript_ID"] = transcript
            coord_df_lst.append(gene_coord_df)

    # Combine results
    final_coord_df = pd.concat(coord_df_lst) if coord_df_lst else pd.DataFrame()
    final_exons_df = pd.DataFrame(exons_coord_df_lst, columns=["ID", "Chr", "Start", "End", "Strand"])

    LOG.info(f"Retrieved coordinates for {len(final_coord_df['Gene'].unique())} genes and {len(final_exons_df['ID'].unique())} exons.")
    return final_coord_df, final_exons_df

# Coordinate parsing and DNA-to-protein mapping functions
# ----------------------------------------------------------
def get_dna_exon_pos(exon_range: np.ndarray, strand: int) -> list:
    """
    Get the DNA positions of an exon given its range and strand.
    Parameters
    ------------
    exon_range : list
        A list containing the start and end positions of the exon.
    strand : int
        The strand of the exon (1 for positive strand, -1 for negative strand).

    Returns
    ------------
    list
        A list of DNA positions corresponding to the exon, ordered according to the strand.
    """

    if strand == -1:
        return np.arange(exon_range[1], exon_range[0] + 1)[::-1]
    else:
        return np.arange(exon_range[0], exon_range[1] + 1)


def get_exon_ix(i: int, exon_range: np.ndarray, strand: int):
    """
    Get the exon index for each DNA position in the exon, ordered according to the strand.

    Parameters
    ------------
    i : int
        The exon index (starting from 0).
    exon_range : np.ndarray
        A numpy array containing the start and end positions of the exon.
    strand : int
        The strand of the exon (1 for positive strand, -1 for negative strand).

    Returns
    ------------
    list
        A list of exon indices corresponding to each DNA position in the exon,
        ordered according to the strand.
    """
    len_exon = len(get_dna_exon_pos(exon_range, strand))

    return np.repeat(i, len_exon)


def get_dna_map_to_protein(coord_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get a mapping of DNA positions to protein positions for a given gene-transcript pair.
    
    Parameters
    ------------
    coord_df : pandas.DataFrame
        A DataFrame containing the coordinates of the CDS for a specific gene-transcript pair.

    Returns
    ------------
    pandas.DataFrame
        A DataFrame containing the mapping of DNA positions to protein positions.
    """
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

def find_exon(dna_prot_df: pd.DataFrame, exon_coord_df: pd.DataFrame):
    """
    For each DNA position in *dna_prot_df* the function finds the matching exon
    ID from *exon_coord_df* by merging on chromosome and filtering by positional
    interval.  This avoids a Python-level loop and is orders of magnitude faster
    than the ``apply`` approach for large panels.

    Parameters
    ------------
    dna_prot_df : pd.DataFrame
        DataFrame with at least ``CHROM`` and ``DNA_POS`` columns.
    exon_coord_df : pd.DataFrame
        DataFrame with ``Chr``, ``Start``, ``End``, and ``ID`` columns.

    Returns
    ------------
    np.ndarray
        Array of exon IDs aligned to *dna_prot_df*\'s row order;
        ``np.nan`` where no matching exon is found.
    """
    # Normalise strand-swapped coordinates so lo <= hi in all cases
    exon_df = exon_coord_df[["Chr", "Start", "End", "ID"]].copy()
    exon_df["lo"] = exon_df[["Start", "End"]].min(axis=1)
    exon_df["hi"] = exon_df[["Start", "End"]].max(axis=1)

    # Cross-merge on chromosome; each row of dna_prot_df is paired only with
    # exons on the same chromosome, then filtered to those whose interval
    # contains the DNA position.
    merged = dna_prot_df[["CHROM", "DNA_POS"]].merge(
        exon_df[["Chr", "lo", "hi", "ID"]],
        left_on="CHROM", right_on="Chr",
        how="left"
    )
    in_exon = (merged["DNA_POS"] >= merged["lo"]) & (merged["DNA_POS"] <= merged["hi"])
    # Keep the first matching exon per position
    matched = merged[in_exon].drop_duplicates(subset=["CHROM", "DNA_POS"])

    # Re-align to the original DataFrame row order
    result = dna_prot_df[["CHROM", "DNA_POS"]].merge(
        matched[["CHROM", "DNA_POS", "ID"]], on=["CHROM", "DNA_POS"], how="left"
    )
    return result["ID"].values


# Depth computation and coverage functions
# ----------------------------------------------------------
def dna2prot_depth(gene_list: list, coord_df: pd.DataFrame, dna_sites: pd.DataFrame, depth_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get a DNA to protein mapping of all positions in the provided list of genes
    Add coverage info & DNA to GENE annotation

    Parameters
    ------------
    gene_list : list
        A list of gene names for which to compute the DNA to protein mapping and coverage information.
    coord_df : pandas.DataFrame
        A DataFrame containing the coordinates of the CDS
        for the genes in the panel.
    dna_sites : pandas.DataFrame
        A DataFrame containing the DNA positions that are part of the consensus panel,
        along with their coverage information.
    depth_df : pandas.DataFrame
        A DataFrame containing the depth information for all positions in the genome.

    Returns
    ------------
    pandas.DataFrame
        A DataFrame containing the mapping of DNA positions to protein positions for the specified genes,
        along with coverage and depth information for each position.
    """

    # Map DNA to protein pos, get exons index to protein pos, etc
    dna_prot_df_lst = []
    for gene in gene_list:
        gene_coord_df = coord_df[coord_df["Gene"] == gene]
        dna_prot_df_lst.append(get_dna_map_to_protein(gene_coord_df))
    dna_prot_df = pd.concat(dna_prot_df_lst)

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


def get_dna2prot_depth(gene_n_transcript_info: pd.DataFrame, depth_file: str, consensus_file: str, release: int, species: str, genome: str, gff3_file: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Function to get the DNA to protein mapping for all positions in the provided list of genes,
    along with coverage and depth information for each position, and the definition of all exons of
    the genes/transcripts in the panel including UTR regions.

    Parameters
    ------------
    gene_n_transcript_info : pandas.DataFrame
        A DataFrame containing gene and transcript information for the panel.
    depth_file : str
        Path to the file containing depth information for all positions in the genome.
    consensus_file : str
        Path to the consensus panel file containing DNA positions and coverage information.
    gff3_file : str | None
        Optional path to a local GFF3 or GFF3.GZ file. If provided, this file
        is used instead of downloading from Ensembl.
    
    Returns
    ------------
    dna_prot_df : pandas.DataFrame
        A DataFrame containing the mapping of DNA positions to protein positions for the specified genes,
        along with coverage and depth information for each position.
    exons_coord_df_final : pandas.DataFrame
        A DataFrame containing the coordinates of all exons (including UTRs) for the genes and transcripts
        in the panel, formatted for BED files.
    """

    consensus_df = pd.read_table(consensus_file)
    depth_df = pd.read_table(depth_file)
    gff_df = gff_to_filtered_df(gene_n_transcript_info, release=release, species=species, genome=genome, gff3_file=gff3_file)

    consensus_df = consensus_df.merge(depth_df[["CHROM", "POS", "CONTEXT"]], on = ["CHROM", "POS"], how = 'left')
    consensus_df = consensus_df.rename(columns={"POS" : "DNA_POS"})

    if "DEPTH" not in depth_df.columns:
        depth_df["DEPTH"] = depth_df.drop(columns=["CHROM", "POS", "CONTEXT"]).mean(1)
    depth_df = depth_df[["CHROM", "POS", "DEPTH"]]

    coord_df, exons_coord_df = get_exon_coord_wrapper(gene_n_transcript_info, gff_df)
    gene_list = list(coord_df["Gene"].unique())
    dna_prot_df = dna2prot_depth(gene_list=gene_list, coord_df=coord_df, dna_sites=consensus_df, depth_df=depth_df)

    dna_prot_df["EXON_ID"] = find_exon(dna_prot_df, exons_coord_df)

    # fix coordinates order in exons_coord_df for BED
    exons_coord_df_final = exons_coord_df.copy()
    exons_coord_df_final.loc[exons_coord_df_final["Strand"] == -1, "Start"] = exons_coord_df.loc[exons_coord_df["Strand"] == -1, "End"]
    exons_coord_df_final.loc[exons_coord_df_final["Strand"] == -1, "End"] = exons_coord_df.loc[exons_coord_df["Strand"] == -1, "Start"]

    return dna_prot_df, exons_coord_df_final


# Plots
# ----------------------------------------------------------
def plot_coverage_per_gene(depths_df: pd.DataFrame) -> None:
    """
    Wrapper function to plot coverage per gene for DNA, protein and exon levels.

    Parameters
    ----------
    depths_df : pandas.DataFrame
        DataFrame containing depth and coverage information for DNA, protein and exon positions.
        In this case, it will be the "exons_depth" created in `get_dna2prot_depth` function.
    """
    column_to_subset = {'DNA': 'DNA_POS',
                        'Protein': 'PROT_POS',
                        'Exon': 'EXON_ID'
                        }
    coverage = depths_df.drop_duplicates(subset=['GENE', 'COVERED', 'DNA_POS'])
    for prefix in ["DNA", "Protein", "Exon"]:
        LOG.info(f"Plotting coverage for {prefix}...")
        coverage_element = coverage.drop_duplicates(subset=['GENE', 'COVERED', column_to_subset[prefix]])
        coverage_summary = coverage_element.groupby(['GENE', 'COVERED']).size().reset_index(name='COUNT')
        plot_single_coverage(coverage_summary, prefix)


def plot_single_coverage(coverage_summary: pd.DataFrame, prefix: str, batch_size: int = 5):
    """
    Plot coverage for a single prefix (DNA, protein or exon) and save the results in a PDF file.

    Parameters
    ----------
    coverage_summary : pandas.DataFrame
        DataFrame containing the coverage summary for each gene and coverage status.
    prefix : str
        The prefix to plot (DNA, Protein or Exon).
    batch_size : int, optional
        The number of genes to include in each batch when plotting (default is 5).
    """
    # Generate copy of the coverage summary to avoid modifying the original DataFrame
    coverage_summary_cp = coverage_summary.copy()

    coverage_pivot = coverage_summary_cp.pivot(index='GENE', columns='COVERED', values='COUNT').fillna(0)
    coverage_pivot = coverage_pivot.sort_values(1, ascending=False).reset_index()
    coverage_pivot = coverage_pivot.set_index('GENE')
    genes_list = coverage_pivot.index.tolist()
    coverage_perc = coverage_pivot.div(coverage_pivot.sum(axis=1), axis=0) * 100
    coverage_pivot.to_csv(f"coverage_count_{prefix}.tsv", header=True, index=True, sep='\t')
    coverage_perc.to_csv(f"coverage_perc_{prefix}.tsv", header=True, index=True, sep='\t')

    if len(genes_list) < batch_size:

        fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

        covered_col = True if True in coverage_pivot.columns else (1 if 1 in coverage_pivot.columns else coverage_pivot.columns[-1])

        coverage_pivot.plot(kind='bar', stacked=True, ax=axes[0], color=['#ff9999','#66b3ff'],
                            )

        axes[0].set_title(f'{prefix} : covered vs. non-covered')
        axes[0].set_ylabel(f'Number of {prefix}')
        axes[0].set_xlabel('Gene')
        axes[0].legend(title='Covered', loc='upper right')
        axes[0].tick_params(axis='x', rotation=45)

        coverage_perc[covered_col].plot(kind='bar', ax=axes[1], color='#66b3ff')
        axes[1].set_title(f'{prefix} : percentage of covered')
        axes[1].set_ylabel('Percentage (%)')
        axes[1].set_xlabel('Gene')
        axes[1].set_ylim(0, 100)
        axes[1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        plt.savefig(f"coverage_per_{prefix}.pdf", dpi=300)
    else:
        with PdfPages(f"coverage_per_{prefix}_batches.pdf") as pdf:
            # split into batches N genes
            for i in range(0, len(genes_list), batch_size): 
                batch_genes = genes_list[i:i+batch_size]
                batch_coverage_pivot = coverage_pivot.loc[batch_genes]
                batch_coverage_perc = coverage_perc.loc[batch_genes]

                fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

                covered_col = True if True in batch_coverage_pivot.columns else (1 if 1 in batch_coverage_pivot.columns else batch_coverage_pivot.columns[-1])

                batch_coverage_pivot.plot(kind='bar', stacked=True, ax=axes[0], color=['#ff9999','#66b3ff'])
                axes[0].set_title(f'{prefix} : covered vs. non-covered (genes {i+1}-{min(i+batch_size, len(genes_list))})')
                axes[0].set_ylabel(f'Number of {prefix}')
                axes[0].set_xlabel('Gene')
                axes[0].legend(title='Covered', loc='upper right')
                axes[0].tick_params(axis='x', rotation=45)

                batch_coverage_perc[covered_col].plot(kind='bar', ax=axes[1], color='#66b3ff')
                axes[1].set_title(f'{prefix} : percentage of covered (genes {i+1}-{min(i+batch_size, len(genes_list))})')
                axes[1].set_ylabel('Percentage (%)')
                axes[1].set_xlabel('Gene')
                axes[1].set_ylim(0, 100)
                axes[1].tick_params(axis='x', rotation=45)

                plt.tight_layout()
                pdf.savefig(dpi=300)
                plt.close()

@click.command()
@click.option('--mutations-file', type=click.Path(exists=True), help='Mutations file')
@click.option('--consensus-file', type=click.Path(exists=True), help='Input consensus panel file')
@click.option('--depths-file', type=click.Path(exists=True), help='Input depths of all samples file')
@click.option('--ensembl-species', type=str, default="homo_sapiens", help='Ensembl species name to use for GFF file retrieval (default)')
@click.option('--ensembl-genome', type=str, default="GRCh38", help='Ensembl genome name to use for GFF file retrieval (default)')
@click.option('--ensembl-release', type=int, default=111, help='Ensembl release number to use for GFF file retrieval (default is 111)')
@click.option('--gff3-file', type=click.Path(exists=True), default=None, help='Optional local GFF3(.gz) file. If provided, skips Ensembl download')
def main(mutations_file, consensus_file, depths_file, ensembl_species, ensembl_genome, ensembl_release, gff3_file):
    # Count each mutation only ones if it appears in multiple reads
    gene_n_transcript = get_transcript_gene_from_maf(mutations_file, consensus_file)

    exons_depth, exons_coord_id = get_dna2prot_depth(
        gene_n_transcript,
        depths_file,
        consensus_file,
        ensembl_release,
        ensembl_species,
        ensembl_genome,
        gff3_file,
    )
    LOG.info("Exons coordinates and depth computed!")
    exons_depth.to_csv("depths_per_position_exon_gene.tsv", header = True, index = False, sep = '\t')

    exons_coordinates_bed_like = exons_coord_id[['Chr', 'Start', 'End', 'ID']]
    exons_coordinates_bed_like.to_csv("panel_exons.bed4.bed", header = False, index = False, sep = '\t')

    plot_coverage_per_gene(exons_depth)

    LOG.info("All done!")


if __name__ == '__main__':
    main()
