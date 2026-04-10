#!/usr/bin/env python3
"""
Codon-align panel BED regions for dndscv build_refCDS.R compatibility.

copied from : https://github.com/bbglab/intogen-plus-dsl2/blob/dev/build-refcds/bin/panel_reformat.py

"""

from __future__ import annotations

import bisect
import csv
import gzip
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass, field

import click

FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

OUTPUT_FIELDS = [
    "GENE_ID",
    "SYMBOL",
    "PROTEIN_ID",
    "CHR",
    "START",
    "END",
    "CDS_START",
    "CDS_END",
    "CDS_LEN",
    "STRAND",
    "TRANSCRIPT_ID",
    "EXON_CHR_START",
    "EXON_CHR_END",
]


@dataclass
class Exon:
    """
    A single CDS exon with genomic and transcript-relative coordinates.
    """

    start: int          # CDS region genomic start (1-based, closed)
    end: int            # CDS region genomic end (1-based, closed)
    cds_start: int      # first position within transcript CDS (1-based)
    cds_end: int        # last position within transcript CDS (1-based)
    exon_chr_start: int # full exon genomic start (may extend into UTR)
    exon_chr_end: int   # full exon genomic end (may extend into UTR)


@dataclass
class Gene:
    """A gene with metadata and CDS exons in 5'→3' transcript order.

    Exons are sorted by ``cds_start`` ascending, which gives transcript
    order for both forward and reverse strand genes because ``cds_start=1``
    always marks the 5'-most exon in the transcript.
    """

    gene_id: str
    symbol: str
    protein_id: str
    chrom: str
    strand: int
    cds_len: int
    transcript_id: str
    exons: list[Exon] = field(default_factory=list)

    @classmethod
    def from_rows(cls, gene_id: str, rows: list[dict[str, str]]) -> Gene:
        """Build a Gene from annotation TSV rows sharing the same symbol.

        Gene-level fields are taken from the first row.  One Exon is created
        per row and the list is sorted by ``cds_start`` to enforce transcript
        order.
        """
        if not rows:
            raise ValueError(f"Gene {gene_id!r} has no annotation rows")
        first = rows[0]
        exons: list[Exon] = []
        for r in rows:
            genomic_start = int(r["START"])
            genomic_end = int(r["END"])
            chr_start_raw = r.get("EXON_CHR_START")
            chr_end_raw = r.get("EXON_CHR_END")
            exons.append(
                Exon(
                    start=genomic_start,
                    end=genomic_end,
                    cds_start=int(r["CDS_START"]),
                    cds_end=int(r["CDS_END"]),
                    exon_chr_start=int(chr_start_raw) if chr_start_raw else genomic_start,
                    exon_chr_end=int(chr_end_raw) if chr_end_raw else genomic_end,
                )
            )
        exons.sort(key=lambda e: e.cds_start)
        return cls(
            gene_id=gene_id,
            symbol=first.get("SYMBOL", gene_id),
            protein_id=first.get("PROTEIN_ID", ""),
            chrom=first["CHR"],
            strand=int(first["STRAND"]),
            cds_len=int(first["CDS_LEN"]),
            transcript_id=first.get("TRANSCRIPT_ID", ""),
            exons=exons,
        )


@dataclass
class Fragment:
    """An aligned genomic fragment ready for mock-transcript construction."""

    start: int
    end: int
    exon: Exon

    @property
    def length(self) -> int:
        """Fragment length in bases (1-based closed)."""
        return self.end - self.start + 1

    @property
    def at_exon_start(self) -> bool:
        """True if this fragment reaches the exon's CDS start boundary."""
        return self.start <= self.exon.start

    @property
    def at_exon_end(self) -> bool:
        """True if this fragment reaches the exon's CDS end boundary."""
        return self.end >= self.exon.end


@dataclass
class PanelGene:
    """Panel-covered fragments of a gene in 5'→3' transcript order.

    Built by intersecting a Gene's exons with panel regions and aligning
    each overlap's 5' boundary to a codon boundary.  The ``gene`` reference
    provides access to gene-level metadata for downstream output.
    """

    gene: Gene
    fragments: list[Fragment]

    @classmethod
    def from_gene(
        cls,
        gene: Gene,
        panel_regions: list[tuple[int, int]],
        aligner: "CodonAligner",
    ) -> "PanelGene":
        """Intersect gene exons with panel regions and build a PanelGene.

        Exons are iterated in CDS order.  For each overlap the natural
        junction state with the previously collected fragment is computed
        before calling the aligner, so the 5' decision is made once and
        correctly.

        Within-exon hits are sorted in transcript direction before being
        appended, ensuring the full list stays in 5'→3' order.
        """
        fragments: list[Fragment] = []
        last_fragment: Fragment | None = None

        # Panel regions are sorted by start coordinate.  Extract start
        # coordinates once for bisect lookups.
        panel_starts = [r[0] for r in panel_regions]

        for exon in gene.exons:
            exon_fragments: list[Fragment] = []

            # Use bisect to skip panel regions that end before this exon.
            # Any region with start > exon.end cannot overlap, so we
            # only scan from the first region whose start could reach
            # back to exon.start.
            lo = bisect.bisect_right(panel_starts, exon.end)
            for idx in range(lo):
                p_start, p_end = panel_regions[idx]
                if p_end < exon.start:
                    continue
                ov_start = max(p_start, exon.start)
                ov_end = min(p_end, exon.end)

                is_junction = cls._check_5prime_junction(
                    last_fragment, ov_start, ov_end, exon, gene.strand,
                )

                if gene.strand == 1:
                    result = aligner.align_forward(ov_start, ov_end, exon, is_junction)
                else:
                    result = aligner.align_reverse(ov_start, ov_end, exon, is_junction)

                if result is None:
                    logger.debug(
                        "Discarded fragment: overlap (%d, %d) collapsed after alignment",
                        ov_start, ov_end,
                    )
                    continue
                exon_fragments.append(Fragment(start=result[0], end=result[1], exon=exon))

            # Sort within-exon hits in transcript order before appending.
            if gene.strand == 1:
                exon_fragments.sort(key=lambda f: f.start)
            else:
                exon_fragments.sort(key=lambda f: -f.end)
            if exon_fragments:
                last_fragment = exon_fragments[-1]
            fragments.extend(exon_fragments)

        return cls(gene=gene, fragments=fragments)

    @staticmethod
    def _check_5prime_junction(
        prev: Fragment | None,
        ov_start: int,
        ov_end: int,
        exon: Exon,
        strand: int,
    ) -> bool:
        """Return True if the overlap forms a natural 5' junction with prev.

        All three conditions must hold:
        1. The previous fragment ends at its exon's natural 3' boundary.
        2. The overlap starts at the current exon's natural 5' boundary.
        3. The two exons are consecutive in the transcript CDS.
        """
        if prev is None:
            return False
        if strand == 1:
            prev_at_3prime = prev.end == prev.exon.end
            curr_at_5prime = ov_start == exon.start
        else:
            prev_at_3prime = prev.start == prev.exon.start
            curr_at_5prime = ov_end == exon.end
        cds_adjacent = prev.exon.cds_end + 1 == exon.cds_start
        return prev_at_3prime and curr_at_5prime and cds_adjacent


class CodonAligner:
    """Align the 5' end of panel-exon overlaps to codon boundaries.

    Separate methods for forward and reverse strands.  The general strategy
    is pad-first (include uncovered CDS bases within the exon span) with
    trim as fallback.  At natural exon junctions the 5' coordinate is left
    untouched so that the split codon spanning the junction is preserved.

    Only the 5' boundary is adjusted here.  The 3' boundary is handled
    during the sequential junction pass in TranscriptBuilder.
    """

    @staticmethod
    def _phase(cds_anchor: int, distance: int) -> int:
        """Return the codon phase (0, 1, or 2) of a coordinate."""
        return (cds_anchor + distance - 1) % 3

    # -- Forward strand (5' = low genomic coordinate) -------------------------

    @staticmethod
    def _pad_or_trim_start_fwd(start: int, phase: int, exon_chr_start: int) -> int:
        """Pad start backward into uncovered CDS or trim forward.

        Padding is preferred: extend by 1-2 bases to reach phase 0 within
        the full exon span.  If padding would exceed the exon boundary, trim
        the covered region instead.
        """
        if phase == 1:
            padded = start - 1
            return padded if padded >= exon_chr_start else start + 2
        if phase == 2:
            padded = start - 2
            return padded if padded >= exon_chr_start else start + 1
        return start

    def align_forward(
        self,
        panel_start: int,
        panel_end: int,
        exon: Exon,
        is_5prime_junction: bool,
    ) -> tuple[int, int] | None:
        """Align the 5' (start) of a forward-strand overlap.

        Three cases:
        1. Natural junction: keep start (split codon continuation).
        2. At CDS boundary, no junction: trim leading orphan bases.
        3. Interior cut: pad backward within exon, trim as fallback.

        Returns adjusted (start, end) or None if the region collapses.
        """
        if is_5prime_junction:
            adj_start = panel_start
        elif panel_start == exon.start:
            phase = (exon.cds_start - 1) % 3
            trim = (3 - phase) % 3
            adj_start = panel_start + trim
        else:
            phase = self._phase(exon.cds_start, panel_start - exon.start)
            adj_start = self._pad_or_trim_start_fwd(panel_start, phase, exon.exon_chr_start)

        if adj_start > panel_end:
            return None
        return (adj_start, panel_end)

    # -- Reverse strand (5' = high genomic coordinate) ------------------------

    @staticmethod
    def _pad_or_trim_end_rev(end: int, phase: int, exon_chr_end: int) -> int:
        """Pad end forward into uncovered CDS or trim backward.

        Mirror of ``_pad_or_trim_start_fwd`` for the reverse strand.
        """
        if phase == 1:
            padded = end + 1
            return padded if padded <= exon_chr_end else end - 2
        if phase == 2:
            padded = end + 2
            return padded if padded <= exon_chr_end else end - 1
        return end

    def align_reverse(
        self,
        panel_start: int,
        panel_end: int,
        exon: Exon,
        is_5prime_junction: bool,
    ) -> tuple[int, int] | None:
        """Align the 5' (end) of a reverse-strand overlap.

        Three cases (mirrored from forward):
        1. Natural junction: keep end (split codon continuation).
        2. At CDS boundary, no junction: trim leading orphan bases.
        3. Interior cut: pad forward within exon, trim as fallback.

        Returns adjusted (start, end) or None if the region collapses.
        """
        if is_5prime_junction:
            adj_end = panel_end
        elif panel_end == exon.end:
            phase = (exon.cds_start - 1) % 3
            trim = (3 - phase) % 3
            adj_end = panel_end - trim
        else:
            phase = self._phase(exon.cds_start, exon.end - panel_end)
            adj_end = self._pad_or_trim_end_rev(panel_end, phase, exon.exon_chr_end)

        if panel_start > adj_end:
            return None
        return (panel_start, adj_end)


class PanelLoader:
    """Load panel regions from a BED file."""

    @staticmethod
    def load(path: str) -> dict[str, list[tuple[int, int]]]:
        r"""Read a BED file and return regions indexed by chromosome.

        Coordinates are converted from 0-based half-open (] to 1-based closed [].
        For downstream compatibility with dndscv, chromosome names have "chr" prefix removed if present.

        Notes
        -----
        e.g. a BED line "chr1\t0\t100" becomes (1, 100) for chromosome "1".
        """
        regions: dict[str, list[tuple[int, int]]] = defaultdict(list)
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                if not parts[1].lstrip("-").isdigit() or not parts[2].lstrip("-").isdigit():
                    logger.debug("Skipping header line: %s", line)
                    continue
                regions[parts[0].replace("chr", "")].append((int(parts[1]), int(parts[2])))
        # Sort regions by start coordinate on each chromosome so that
        # downstream code can use bisect to find overlaps efficiently.
        for chrom in regions:
            regions[chrom].sort()
        return dict(regions)


class GeneLoader:
    """Load gene annotations from a TSV file."""

    @staticmethod
    def load(path: str) -> dict[str, Gene]:
        """Read a gene annotation TSV and return Gene objects keyed by GENE_ID.

        Rows are grouped by GENE_ID first, then each group is passed to
        ``Gene.from_rows`` which extracts gene-level metadata from the first
        row and builds Exon objects sorted in transcript order.
        """
        raw_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
        with open(path, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t", fieldnames=OUTPUT_FIELDS):
                raw_rows[row["GENE_ID"]].append(row)
        genes: dict[str, Gene] = {}
        for gene_id, rows in raw_rows.items():
            genes[gene_id] = Gene.from_rows(gene_id, rows)
        return genes


class TranscriptBuilder:
    """Align panel fragments to codon boundaries and build mock transcripts."""

    def __init__(self, aligner: CodonAligner) -> None:
        self._aligner = aligner

    @staticmethod
    def _is_natural_junction(prev: Fragment, frag: Fragment, strand: int) -> bool:
        """Check if two consecutive fragments form a natural CDS continuation.

        A natural junction means:
        1. The previous fragment's 3' end is at the exon's natural boundary.
        2. The current fragment's 5' start is at the exon's natural boundary.
        3. The two exons are adjacent in the original transcript's CDS.
        """
        if strand == 1:
            prev_3_natural = prev.end == prev.exon.end
            curr_5_natural = frag.start == frag.exon.start
        else:
            prev_3_natural = prev.start == prev.exon.start
            curr_5_natural = frag.end == frag.exon.end
        cds_adjacent = prev.exon.cds_end + 1 == frag.exon.cds_start
        return prev_3_natural and curr_5_natural and cds_adjacent

    @staticmethod
    def _ensure_codon_junctions(fragments: list[Fragment], strand: int) -> list[Fragment]:
        """Trim 3' orphan bases at artificial junctions and enforce total CDS length % 3 == 0.

        By the time fragments reach this method, ``align()`` (called from
        ``PanelGene.from_gene``) has already resolved every 5' boundary:
        - Natural junctions: 5' left untouched (split codon preserved).
        - Natural CDS boundary without junction: leading orphans trimmed.
        - Panel cuts into exon interior: pad/trim to phase 0.

        This method only needs to handle the 3' side:
        1. At an artificial junction (panel cut or non-adjacent exons), trim
           trailing orphan bases from the 3' end of the preceding group so
           the accumulated length is a multiple of 3 before the next phase-0
           fragment begins.
        2. At a natural junction (CDS-adjacent exons, both boundaries natural),
           the split codon is biologically real; no trimming is applied.
        3. After the full sequential pass, trim the 3' end of the last
           fragment to make the total CDS length a multiple of 3.
        """
        if not fragments:
            return fragments

        first = fragments[0]
        running_bases = first.length if first.start <= first.end else 0

        for i in range(1, len(fragments)):
            frag = fragments[i]
            prev = fragments[i - 1]

            natural = TranscriptBuilder._is_natural_junction(prev, frag, strand)

            if not natural:
                # Trim trailing orphan bases from the 3' end of the preceding
                # group.  The next fragment already starts at phase 0 (handled
                # by align()), so no 5' adjustment of frag is needed here.
                orphan = running_bases % 3
                if orphan > 0:
                    if strand == 1:
                        prev.end -= orphan
                    else:
                        prev.start += orphan
                    running_bases -= orphan
                    logger.debug("Trimmed %d orphan base(s) at artificial junction", orphan)

            running_bases += frag.length

        # Final: trim total to multiple of 3
        remainder = running_bases % 3
        if remainder > 0 and fragments:
            last = fragments[-1]
            if strand == 1:
                last.end -= remainder
            else:
                last.start += remainder
            logger.debug("Trimmed %d trailing base(s) for total codon alignment", remainder)

        return [f for f in fragments if f.start <= f.end]

    @staticmethod
    def _reindex(gene: Gene, fragments: list[Fragment]) -> list[dict[str, object]]:
        """Assign continuous CDS pointers to ordered fragments."""
        rows: list[dict[str, object]] = []
        cumulative = 0
        for frag in fragments:
            new_cds_start = cumulative + 1
            new_cds_end = cumulative + frag.length
            cumulative = new_cds_end
            rows.append(
                {
                    "GENE_ID": gene.gene_id,
                    "SYMBOL": gene.symbol,
                    "PROTEIN_ID": gene.protein_id,
                    "CHR": gene.chrom,
                    "START": frag.start,
                    "END": frag.end,
                    "CDS_START": new_cds_start,
                    "CDS_END": new_cds_end,
                    "CDS_LEN": 0,
                    "STRAND": gene.strand,
                    "TRANSCRIPT_ID": gene.transcript_id,
                    "EXON_CHR_START": frag.start,
                    "EXON_CHR_END": frag.end,
                }
            )
        for row in rows:
            row["CDS_LEN"] = cumulative
        return rows

    def build(
        self,
        gene: Gene,
        panel_by_chr: dict[str, list[tuple[int, int]]],
    ) -> tuple[list[dict[str, object]], list[Fragment]]:
        """Process a single gene and return mock-transcript rows and fragments.

        The pipeline is three sequential steps, each operating on the
        ``PanelGene`` produced in the first step:

        1. Intersect the gene's exons with the panel and 5'-align each
           overlap to a codon boundary  →  ``PanelGene``
        2. Fix codon junctions between consecutive fragments in transcript
           order, trimming orphan bases at artificial cuts.
        3. Reindex fragments into a continuous mock CDS.

        Returns (output_rows, final_fragments) so callers can inspect which
        fragment boundaries are real vs. artificial for splice site reporting.
        """
        panel_regions = panel_by_chr.get(gene.chrom, [])

        panel_gene = PanelGene.from_gene(gene, panel_regions, self._aligner)
        if not panel_gene.fragments:
            return [], []

        panel_gene.fragments = self._ensure_codon_junctions(panel_gene.fragments, panel_gene.gene.strand)
        if not panel_gene.fragments:
            return [], []

        rows = self._reindex(panel_gene.gene, panel_gene.fragments)
        logger.debug("%s: %d fragment(s), total CDS length = %d", panel_gene.gene.symbol, len(rows), rows[0]["CDS_LEN"])
        return rows, panel_gene.fragments


SPLICE_REPORT_FIELDS = [
    "CHR",
    "POSITION",
    "GENE_SYMBOL",
    "GENE_ID",
    "JUNCTION_TYPE",
    "SPLICE_TYPE",
    "FRAGMENT_BOUNDARY",
]


class SpliceClassifier:
    """Classify splice sites that dndscv will infer from the mock transcript.

    dndscv computes essential splice sites from consecutive CDS exon rows:
      - Donor (5'ss):   prev_END +1, +2, +5   (forward)  or  next_START -1, -2, -5  (reverse)
      - Acceptor (3'ss): next_START -1, -2     (forward)  or  prev_END +1, +2       (reverse)

    Every junction between consecutive output fragments becomes a splice
    site source.  This class determines whether each junction is a REAL
    exon-intron boundary or an ARTIFICIAL panel cut, allowing users to
    post-filter dndscv splice mutation calls.
    """

    @staticmethod
    def classify(
        gene: Gene,
        fragments: list[Fragment],
    ) -> list[dict[str, object]]:
        """Return a list of splice site records for every inter-fragment junction.

        Each record contains the genomic position, gene info, whether the
        junction is real or artificial, and the splice type (donor/acceptor).
        """
        if len(fragments) < 2:
            return []

        records: list[dict[str, object]] = []
        for i in range(len(fragments) - 1):
            prev_frag = fragments[i]
            next_frag = fragments[i + 1]

            # Determine if this junction is a real exon-intron boundary.
            # Real requires: prev fragment reaches its exon's 3' CDS end,
            # next fragment reaches its exon's 5' CDS start, and the two
            # exons are CDS-adjacent in the original transcript.
            if gene.strand == 1:
                prev_at_3prime = prev_frag.at_exon_end
                next_at_5prime = next_frag.at_exon_start
            else:
                prev_at_3prime = prev_frag.at_exon_start
                next_at_5prime = next_frag.at_exon_end
            cds_adjacent = prev_frag.exon.cds_end + 1 == next_frag.exon.cds_start
            is_real = prev_at_3prime and next_at_5prime and cds_adjacent
            junction_type = "real" if is_real else "artificial"

            # Compute the positions dndscv will use, mirroring its
            # get_splicesites() logic.
            if gene.strand == 1:
                donor_positions = [prev_frag.end + 1, prev_frag.end + 2, prev_frag.end + 5]
                acceptor_positions = [next_frag.start - 1, next_frag.start - 2]
            else:
                donor_positions = [next_frag.start - 1, next_frag.start - 2, next_frag.start - 5]
                acceptor_positions = [prev_frag.end + 1, prev_frag.end + 2]

            boundary_desc_donor = f"{prev_frag.end}|{next_frag.start}"
            for pos in donor_positions:
                records.append({
                    "CHR": gene.chrom,
                    "POSITION": pos,
                    "GENE_SYMBOL": gene.symbol,
                    "GENE_ID": gene.gene_id,
                    "JUNCTION_TYPE": junction_type,
                    "SPLICE_TYPE": "donor",
                    "FRAGMENT_BOUNDARY": boundary_desc_donor,
                })
            for pos in acceptor_positions:
                records.append({
                    "CHR": gene.chrom,
                    "POSITION": pos,
                    "GENE_SYMBOL": gene.symbol,
                    "GENE_ID": gene.gene_id,
                    "JUNCTION_TYPE": junction_type,
                    "SPLICE_TYPE": "acceptor",
                    "FRAGMENT_BOUNDARY": boundary_desc_donor,
                })
        return records

@click.command(context_settings={"help_option_names": ["-h", "--help"], "show_default": True})
@click.version_option(version="0.1.0", prog_name="panel_reformat")
@click.option("-b", "--bed", type=click.Path(exists=True), required=True, help="Panel BED file (CHR, START, END).")
@click.option("-g", "--genes", type=click.Path(exists=True), required=True, help="Gene annotation TSV file.")
@click.option("-o", "--output", type=click.Path(), default="-", help="Output TSV path.")
@click.option("-s", "--splice-report", type=click.Path(), default="splice_sites.tsv",
              help="Output TSV with splice site classification (real vs artificial).")
@click.option("-v", "--verbose", is_flag=True, default=False, help="Enable debug logging.")
def cli(bed: str, genes: str, output: str, splice_report: str | None, verbose: bool) -> None:
    """Align panel BED regions to codon boundaries for dndscv build_refCDS.R."""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    panel_by_chr = PanelLoader.load(bed)
    gene_map = GeneLoader.load(genes)
    logger.info("Chromosomes in panel: %s", ", ".join(sorted(panel_by_chr)))
    logger.info("Reference cds region: %d genes.", len(gene_map))

    # Only process genes on chromosomes covered by the panel.
    panel_chroms = set(panel_by_chr)
    candidate_genes = {gid: g for gid, g in gene_map.items() if g.chrom in panel_chroms}
    logger.info("Genes on panel chromosomes: %d (skipping %d).",
                len(candidate_genes), len(gene_map) - len(candidate_genes))

    builder = TranscriptBuilder(CodonAligner())
    splice_classifier = SpliceClassifier()

    out_fh = sys.stdout if output == "-" else open(output, "w", encoding="utf-8")
    try:
        writer = csv.DictWriter(out_fh, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()

        gene_count = 0
        fragment_count = 0
        all_splice_records: list[dict[str, object]] = []
        for gene_id in sorted(candidate_genes):
            gene = candidate_genes[gene_id]
            rows, fragments = builder.build(gene, panel_by_chr)
            for row in rows:
                writer.writerow(row)
            if rows:
                gene_count += 1
                fragment_count += len(rows)
                if splice_report:
                    all_splice_records.extend(
                        splice_classifier.classify(gene, fragments)
                    )

        logger.info("Wrote %d fragment(s) across %d gene(s).", fragment_count, gene_count)
    finally:
        if output != "-":
            out_fh.close()

    if splice_report and all_splice_records:
        real_count = sum(1 for r in all_splice_records if r["JUNCTION_TYPE"] == "real")
        artificial_count = len(all_splice_records) - real_count
        logger.info(
            "Splice report: %d positions (%d real, %d artificial).",
            len(all_splice_records), real_count, artificial_count,
        )
        with open(splice_report, "w", encoding="utf-8") as sfh:
            sw = csv.DictWriter(sfh, fieldnames=SPLICE_REPORT_FIELDS, delimiter="\t")
            sw.writeheader()
            for rec in all_splice_records:
                sw.writerow(rec)


if __name__ == "__main__":
    cli() # pylint: disable=no-value-for-parameter