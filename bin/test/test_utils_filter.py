#!/usr/bin/env python3
"""
Unit tests for utils_filter.py.

Covers:
  - somatic_mask / germline_mask           (each and combined)
  - filter_maf                             (all criterion branches)
  - load_filter_criteria                   (parsing, combining, prefix stripping)
  - expand_filter_column                   (boolean column creation + required columns)
  - extract_flagged_regions_bed            (empty and non-empty BED output)
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

# Add the bin directory to the path to import sibling modules
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils_filter import (
    expand_filter_column,
    extract_flagged_regions_bed,
    filter_maf,
    germline_mask,
    load_filter_criteria,
    somatic_mask,
)

THRESHOLD = 0.3

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_df(rows: list[tuple[float, float, float]]) -> pd.DataFrame:
    """Build a minimal MAF DataFrame with VAF, vd_VAF, and VAF_AM columns."""
    vafs, vd_vafs, vaf_ams = zip(*rows)
    return pd.DataFrame({"VAF": list(vafs), "vd_VAF": list(vd_vafs), "VAF_AM": list(vaf_ams)})


def _make_maf(rows: list[dict]) -> pd.DataFrame:
    """Build a MAF DataFrame from a list of row dicts, preserving column order."""
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# somatic_mask
# ---------------------------------------------------------------------------


class TestSomaticMask(unittest.TestCase):
    """Tests for somatic_mask(maf_df, threshold)."""

    def test_all_below_threshold_is_somatic(self):
        """All three VAF columns strictly below threshold → somatic True."""
        df = _make_df([(0.1, 0.2, 0.05)])
        result = somatic_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [True])

    def test_all_above_threshold_is_not_somatic(self):
        """All three VAF columns strictly above threshold → somatic False."""
        df = _make_df([(0.5, 0.6, 0.4)])
        result = somatic_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [False])

    def test_boundary_equality_is_somatic(self):
        """All three VAF columns exactly equal to threshold → somatic True (≤ is inclusive)."""
        df = _make_df([(THRESHOLD, THRESHOLD, THRESHOLD)])
        result = somatic_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [True])

    def test_asymmetric_is_not_somatic(self):
        """Mixed columns (some ≤ threshold, some > threshold) → somatic False."""
        df = _make_df([(0.1, 0.1, 0.5)])
        result = somatic_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [False])

    def test_multiple_rows(self):
        """Full four-case matrix in a single DataFrame."""
        df = _make_df(
            [
                (0.1, 0.2, 0.05),  # all-below → True
                (0.5, 0.6, 0.4),  # all-above → False
                (THRESHOLD, THRESHOLD, THRESHOLD),  # boundary → True
                (0.1, 0.1, 0.5),  # asymmetric → False
            ]
        )
        result = somatic_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [True, False, True, False])


# ---------------------------------------------------------------------------
# germline_mask
# ---------------------------------------------------------------------------


class TestGermlineMask(unittest.TestCase):
    """Tests for germline_mask(maf_df, threshold)."""

    def test_all_below_threshold_is_not_germline(self):
        """All three VAF columns strictly below threshold → germline False."""
        df = _make_df([(0.1, 0.2, 0.05)])
        result = germline_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [False])

    def test_all_above_threshold_is_germline(self):
        """All three VAF columns strictly above threshold → germline True."""
        df = _make_df([(0.5, 0.6, 0.4)])
        result = germline_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [True])

    def test_boundary_equality_is_not_germline(self):
        """All three VAF columns exactly equal to threshold → germline False (> is exclusive)."""
        df = _make_df([(THRESHOLD, THRESHOLD, THRESHOLD)])
        result = germline_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [False])

    def test_asymmetric_is_not_germline(self):
        """Mixed columns (some ≤ threshold, some > threshold) → germline False."""
        df = _make_df([(0.1, 0.1, 0.5)])
        result = germline_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [False])

    def test_multiple_rows(self):
        """Full four-case matrix in a single DataFrame."""
        df = _make_df(
            [
                (0.1, 0.2, 0.05),  # all-below → False
                (0.5, 0.6, 0.4),  # all-above → True
                (THRESHOLD, THRESHOLD, THRESHOLD),  # boundary → False
                (0.1, 0.1, 0.5),  # asymmetric → False
            ]
        )
        result = germline_mask(df, THRESHOLD)
        self.assertListEqual(result.tolist(), [False, True, False, False])


# ---------------------------------------------------------------------------
# somatic + germline masks are never simultaneously True
# ---------------------------------------------------------------------------


class TestMasksAreNotComplements(unittest.TestCase):
    """Verify that somatic and germline masks are never simultaneously True."""

    def test_no_row_is_true_in_both_masks(self):
        """No row should satisfy both somatic and germline conditions at once."""
        df = _make_df(
            [
                (0.1, 0.2, 0.05),  # all-below
                (0.5, 0.6, 0.4),  # all-above
                (THRESHOLD, THRESHOLD, THRESHOLD),  # boundary
                (0.1, 0.1, 0.5),  # asymmetric
            ]
        )
        somatic = somatic_mask(df, THRESHOLD)
        germline = germline_mask(df, THRESHOLD)
        both_true = (somatic & germline).tolist()
        self.assertListEqual(both_true, [False, False, False, False])

    def test_asymmetric_row_is_false_in_both_masks(self):
        """An asymmetric row must be False in somatic AND False in germline (the 'neither' case)."""
        df = _make_df([(0.1, 0.1, 0.5)])
        self.assertFalse(somatic_mask(df, THRESHOLD).iloc[0])
        self.assertFalse(germline_mask(df, THRESHOLD).iloc[0])


# ---------------------------------------------------------------------------
# filter_maf
# ---------------------------------------------------------------------------


class TestFilterMaf(unittest.TestCase):
    """Tests for filter_maf(maf_df, filter_criteria)."""

    def _base_maf(self) -> pd.DataFrame:
        """Return a small MAF DataFrame exercising all criterion branches."""
        return _make_maf(
            [
                {"MUT_ID": "M1", "VAF": 0.1, "DEPTH": 50, "FILTER": "PASS", "TYPE": "SNV", "FILTER.not_covered": False},
                {
                    "MUT_ID": "M2",
                    "VAF": 0.4,
                    "DEPTH": 30,
                    "FILTER": "n_rich;NM20",
                    "TYPE": "SNV",
                    "FILTER.not_covered": False,
                },
                {
                    "MUT_ID": "M3",
                    "VAF": 0.2,
                    "DEPTH": 60,
                    "FILTER": "low_mappability",
                    "TYPE": "INDEL",
                    "FILTER.not_covered": True,
                },
                {
                    "MUT_ID": "M4",
                    "VAF": 0.05,
                    "DEPTH": 80,
                    "FILTER": "PASS",
                    "TYPE": "SNV",
                    "FILTER.not_covered": False,
                },
            ]
        )

    # --- numeric operator branch (len(operator)==2) ---

    def test_numeric_le_filters_correctly(self):
        """('VAF', 'le 0.3') keeps only rows with VAF ≤ 0.3."""
        df = self._base_maf()
        result = filter_maf(df, [("VAF", "le 0.3")])
        self.assertListEqual(sorted(result["MUT_ID"].tolist()), ["M1", "M3", "M4"])

    def test_numeric_ge_filters_correctly(self):
        """('DEPTH', 'ge 50') keeps only rows with DEPTH ≥ 50."""
        df = self._base_maf()
        result = filter_maf(df, [("DEPTH", "ge 50")])
        self.assertListEqual(sorted(result["MUT_ID"].tolist()), ["M1", "M3", "M4"])

    def test_numeric_lt_filters_correctly(self):
        """('VAF', 'lt 0.2') keeps only rows with VAF < 0.2."""
        df = self._base_maf()
        result = filter_maf(df, [("VAF", "lt 0.2")])
        self.assertListEqual(sorted(result["MUT_ID"].tolist()), ["M1", "M4"])

    def test_multiple_numeric_criteria_are_anded(self):
        """Combining two numeric criteria narrows the result set."""
        df = self._base_maf()
        result = filter_maf(df, [("VAF", "le 0.3"), ("DEPTH", "ge 50")])
        self.assertListEqual(sorted(result["MUT_ID"].tolist()), ["M1", "M3", "M4"])

    # --- notcontains / contains branch ---

    def test_notcontains_excludes_matching_filter_token(self):
        """('FILTER', 'notcontains n_rich') removes rows whose FILTER cell contains 'n_rich'."""
        df = self._base_maf()
        result = filter_maf(df, [("FILTER", "notcontains n_rich")])
        # M2 has 'n_rich' in its FILTER → removed
        self.assertNotIn("M2", result["MUT_ID"].tolist())
        self.assertIn("M1", result["MUT_ID"].tolist())

    def test_notcontains_respects_semicolon_split(self):
        """notcontains splits on ';' so a token that is a prefix of another is handled correctly."""
        df = _make_maf(
            [
                {"MUT_ID": "A", "FILTER": "NM20;PASS"},
                {"MUT_ID": "B", "FILTER": "NM200;PASS"},  # 'NM20' is NOT a token here
                {"MUT_ID": "C", "FILTER": "PASS"},
            ]
        )
        result = filter_maf(df, [("FILTER", "notcontains NM20")])
        self.assertNotIn("A", result["MUT_ID"].tolist())
        self.assertIn("B", result["MUT_ID"].tolist())
        self.assertIn("C", result["MUT_ID"].tolist())

    def test_contains_keeps_only_matching_filter_token(self):
        """('FILTER', 'contains n_rich') keeps only rows whose FILTER cell contains 'n_rich'."""
        df = self._base_maf()
        result = filter_maf(df, [("FILTER", "contains n_rich")])
        self.assertListEqual(result["MUT_ID"].tolist(), ["M2"])

    # --- boolean column branch ---

    def test_boolean_criterion_true_selects_matching_rows(self):
        """('FILTER.not_covered', True) keeps only rows where FILTER.not_covered is True."""
        df = self._base_maf()
        result = filter_maf(df, [("FILTER.not_covered", True)])
        self.assertListEqual(result["MUT_ID"].tolist(), ["M3"])

    def test_boolean_criterion_false_selects_matching_rows(self):
        """('FILTER.not_covered', False) keeps only rows where FILTER.not_covered is False."""
        df = self._base_maf()
        result = filter_maf(df, [("FILTER.not_covered", False)])
        self.assertListEqual(sorted(result["MUT_ID"].tolist()), ["M1", "M2", "M4"])

    # --- plain-value (equality) branch ---

    def test_plain_value_equality_match(self):
        """('TYPE', 'SNV') keeps only rows where TYPE == 'SNV'."""
        df = self._base_maf()
        result = filter_maf(df, [("TYPE", "SNV")])
        self.assertListEqual(sorted(result["MUT_ID"].tolist()), ["M1", "M2", "M4"])

    def test_plain_value_no_match_returns_empty(self):
        """('TYPE', 'NONEXISTENT') returns an empty DataFrame."""
        df = self._base_maf()
        result = filter_maf(df, [("TYPE", "NONEXISTENT")])
        self.assertEqual(len(result), 0)

    # --- no criteria leaves DataFrame unchanged ---

    def test_empty_criteria_returns_all_rows(self):
        """An empty criteria list leaves the DataFrame unchanged."""
        df = self._base_maf()
        result = filter_maf(df, [])
        self.assertEqual(len(result), len(df))


# ---------------------------------------------------------------------------
# load_filter_criteria
# ---------------------------------------------------------------------------


class TestLoadFilterCriteria(unittest.TestCase):
    """Tests for load_filter_criteria(filters, somatic_filters)."""

    def test_extracts_notcontains_entries_from_filters(self):
        """Items starting with 'notcontains ' in filters are returned with prefix stripped."""
        result = load_filter_criteria("notcontains n_rich,notcontains NM20", "")
        self.assertListEqual(sorted(result), ["NM20", "n_rich"])

    def test_extracts_notcontains_entries_from_somatic_filters(self):
        """Items from somatic_filters starting with 'notcontains ' are included."""
        result = load_filter_criteria("", "notcontains low_mappability")
        self.assertListEqual(result, ["low_mappability"])

    def test_combines_both_lists(self):
        """Entries from both arguments are merged before filtering."""
        result = load_filter_criteria("notcontains n_rich", "notcontains NM20")
        self.assertListEqual(sorted(result), ["NM20", "n_rich"])

    def test_non_notcontains_entries_are_excluded(self):
        """Entries that do not start with 'notcontains ' are silently dropped."""
        result = load_filter_criteria("notcontains n_rich,VAF le 0.3,PASS", "")
        self.assertListEqual(result, ["n_rich"])

    def test_empty_strings_return_empty_list(self):
        """Both arguments being empty strings yields an empty list."""
        result = load_filter_criteria("", "")
        self.assertListEqual(result, [])

    def test_whitespace_is_trimmed_around_items(self):
        """Leading/trailing whitespace around comma-separated items is stripped."""
        result = load_filter_criteria("  notcontains n_rich  ,  notcontains NM20  ", "")
        self.assertListEqual(sorted(result), ["NM20", "n_rich"])

    def test_duplicate_entries_are_preserved(self):
        """Duplicates across both arguments are preserved (no deduplication contract)."""
        result = load_filter_criteria("notcontains n_rich", "notcontains n_rich")
        self.assertEqual(result.count("n_rich"), 2)


# ---------------------------------------------------------------------------
# expand_filter_column
# ---------------------------------------------------------------------------


class TestExpandFilterColumn(unittest.TestCase):
    """Tests for expand_filter_column(maf_df)."""

    def _make_filter_df(self, filter_values: list[str]) -> pd.DataFrame:
        """Build a minimal MAF DataFrame with only a FILTER column."""
        return pd.DataFrame({"FILTER": filter_values})

    def test_creates_boolean_column_for_each_token(self):
        """Each unique ';'-delimited token gets its own FILTER.<token> boolean column."""
        df = self._make_filter_df(["n_rich;NM20", "PASS", "NM20"])
        result = expand_filter_column(df)
        self.assertIn("FILTER.n_rich", result.columns)
        self.assertIn("FILTER.NM20", result.columns)
        self.assertIn("FILTER.PASS", result.columns)

    def test_boolean_values_are_correct(self):
        """True only where the token is present in that row's FILTER value."""
        df = self._make_filter_df(["n_rich;NM20", "PASS", "NM20"])
        result = expand_filter_column(df)
        # Row 0: n_rich and NM20 present
        self.assertTrue(result.loc[0, "FILTER.n_rich"])
        self.assertTrue(result.loc[0, "FILTER.NM20"])
        self.assertFalse(result.loc[0, "FILTER.PASS"])
        # Row 1: only PASS present
        self.assertFalse(result.loc[1, "FILTER.n_rich"])
        self.assertTrue(result.loc[1, "FILTER.PASS"])
        # Row 2: only NM20 present
        self.assertFalse(result.loc[2, "FILTER.n_rich"])
        self.assertTrue(result.loc[2, "FILTER.NM20"])

    def test_required_columns_always_exist(self):
        """FILTER.not_covered and FILTER.not_in_exons are always created even if absent in data."""
        df = self._make_filter_df(["PASS", "PASS"])
        result = expand_filter_column(df)
        self.assertIn("FILTER.not_covered", result.columns)
        self.assertIn("FILTER.not_in_exons", result.columns)

    def test_required_columns_are_false_when_token_absent(self):
        """Required columns are all False when neither token appears in the data."""
        df = self._make_filter_df(["PASS", "n_rich"])
        result = expand_filter_column(df)
        self.assertFalse(result["FILTER.not_covered"].any())
        self.assertFalse(result["FILTER.not_in_exons"].any())

    def test_single_token_per_row(self):
        """A FILTER column with no semicolons creates one boolean column per distinct value."""
        df = self._make_filter_df(["alpha", "beta", "alpha"])
        result = expand_filter_column(df)
        self.assertIn("FILTER.alpha", result.columns)
        self.assertIn("FILTER.beta", result.columns)
        self.assertListEqual(result["FILTER.alpha"].tolist(), [True, False, True])
        self.assertListEqual(result["FILTER.beta"].tolist(), [False, True, False])

    def test_all_required_columns_true_when_token_present(self):
        """FILTER.not_covered is True exactly for the rows that contain 'not_covered'."""
        df = self._make_filter_df(["not_covered;n_rich", "PASS", "not_covered"])
        result = expand_filter_column(df)
        self.assertListEqual(result["FILTER.not_covered"].tolist(), [True, False, True])


# ---------------------------------------------------------------------------
# extract_flagged_regions_bed
# ---------------------------------------------------------------------------


class TestExtractFlaggedRegionsBed(unittest.TestCase):
    """Tests for extract_flagged_regions_bed(maf_df, name, filters, specification)."""

    def setUp(self):
        """Switch into a fresh temporary directory for each test; restore on teardown."""
        self._tmpdir = tempfile.mkdtemp()
        self._orig_dir = os.getcwd()
        os.chdir(self._tmpdir)

    def tearDown(self):
        """Restore original working directory."""
        os.chdir(self._orig_dir)

    def _make_expanded_maf(self, rows: list[dict]) -> pd.DataFrame:
        """Build a MAF with CHROM/POS columns then run expand_filter_column."""
        df = pd.DataFrame(rows)
        return expand_filter_column(df)

    # --- empty case ---

    def test_empty_case_creates_empty_bed_file(self):
        """No flagged rows → an empty .bed file is touched and the function returns None."""
        df = self._make_expanded_maf(
            [
                {"CHROM": "chr1", "POS": 100, "FILTER": "PASS"},
                {"CHROM": "chr1", "POS": 200, "FILTER": "PASS"},
            ]
        )
        result = extract_flagged_regions_bed(df, "sample1", ["n_rich"])
        self.assertIsNone(result)
        bed_path = Path("sample1.flagged-pos.bed")
        self.assertTrue(bed_path.exists())
        self.assertEqual(bed_path.stat().st_size, 0)

    def test_empty_case_with_specification_uses_correct_filename(self):
        """specification parameter is included in the BED file name for the empty case."""
        df = self._make_expanded_maf([{"CHROM": "chr1", "POS": 100, "FILTER": "PASS"}])
        extract_flagged_regions_bed(df, "sample1", ["n_rich"], specification="cohort-")
        self.assertTrue(Path("sample1.cohort-flagged-pos.bed").exists())

    def test_empty_case_no_matching_filter_columns(self):
        """When filter names have no corresponding FILTER.* columns, result is empty BED."""
        df = self._make_expanded_maf([{"CHROM": "chr1", "POS": 100, "FILTER": "PASS"}])
        result = extract_flagged_regions_bed(df, "sampleX", ["nonexistent_filter"])
        self.assertIsNone(result)
        self.assertTrue(Path("sampleX.flagged-pos.bed").exists())

    # --- non-empty case ---

    def test_nonempty_case_writes_bed_with_correct_columns(self):
        """BED file has four tab-separated columns: CHROM, START, END, FILTERS."""
        df = self._make_expanded_maf(
            [
                {"CHROM": "chr1", "POS": 500, "FILTER": "n_rich"},
                {"CHROM": "chr2", "POS": 1000, "FILTER": "PASS"},
            ]
        )
        extract_flagged_regions_bed(df, "sample2", ["n_rich"])
        bed_path = Path("sample2.flagged-pos.bed")
        self.assertTrue(bed_path.exists())
        bed = pd.read_csv(bed_path, sep="\t", header=None, names=["CHROM", "START", "END", "FILTERS"])
        self.assertEqual(len(bed), 1)
        row = bed.iloc[0]
        self.assertEqual(row["CHROM"], "chr1")
        self.assertEqual(row["START"], 500)
        self.assertEqual(row["END"], 500)
        self.assertIn("FILTER.n_rich", row["FILTERS"])

    def test_nonempty_case_multiple_filters_joined_with_comma(self):
        """When a position has two active filter flags, FILTERS column is comma-joined."""
        df = self._make_expanded_maf(
            [
                {"CHROM": "chr1", "POS": 300, "FILTER": "n_rich;NM20"},
            ]
        )
        extract_flagged_regions_bed(df, "sample3", ["n_rich", "NM20"])
        bed = pd.read_csv(
            Path("sample3.flagged-pos.bed"), sep="\t", header=None, names=["CHROM", "START", "END", "FILTERS"]
        )
        self.assertEqual(len(bed), 1)
        filters_value = bed.iloc[0]["FILTERS"]
        # Both filter column names should appear, joined by comma
        self.assertIn("FILTER.n_rich", filters_value)
        self.assertIn("FILTER.NM20", filters_value)
        self.assertIn(",", filters_value)

    def test_nonempty_case_multiple_rows_all_written(self):
        """Multiple flagged positions each produce a row in the BED file."""
        df = self._make_expanded_maf(
            [
                {"CHROM": "chr1", "POS": 100, "FILTER": "n_rich"},
                {"CHROM": "chr1", "POS": 200, "FILTER": "NM20"},
                {"CHROM": "chr2", "POS": 50, "FILTER": "PASS"},
            ]
        )
        extract_flagged_regions_bed(df, "sample4", ["n_rich", "NM20"])
        bed = pd.read_csv(
            Path("sample4.flagged-pos.bed"), sep="\t", header=None, names=["CHROM", "START", "END", "FILTERS"]
        )
        self.assertEqual(len(bed), 2)
        self.assertSetEqual(set(bed["START"].tolist()), {100, 200})

    def test_nonempty_case_returns_none(self):
        """Function has no explicit return in the non-empty path, so returns None."""
        df = self._make_expanded_maf([{"CHROM": "chr1", "POS": 100, "FILTER": "n_rich"}])
        result = extract_flagged_regions_bed(df, "sample5", ["n_rich"])
        self.assertIsNone(result)

    def test_nonempty_case_with_specification_uses_correct_filename(self):
        """specification parameter is included in the BED file name for the non-empty case."""
        df = self._make_expanded_maf([{"CHROM": "chr1", "POS": 100, "FILTER": "n_rich"}])
        extract_flagged_regions_bed(df, "sample6", ["n_rich"], specification="cohort-")
        self.assertTrue(Path("sample6.cohort-flagged-pos.bed").exists())

    def test_nonempty_case_bed_is_sorted_by_chrom_and_pos(self):
        """BED rows are sorted by CHROM then POS (ascending)."""
        df = self._make_expanded_maf(
            [
                {"CHROM": "chr2", "POS": 800, "FILTER": "n_rich"},
                {"CHROM": "chr1", "POS": 999, "FILTER": "n_rich"},
                {"CHROM": "chr1", "POS": 100, "FILTER": "n_rich"},
            ]
        )
        extract_flagged_regions_bed(df, "sample7", ["n_rich"])
        bed = pd.read_csv(
            Path("sample7.flagged-pos.bed"), sep="\t", header=None, names=["CHROM", "START", "END", "FILTERS"]
        )
        self.assertEqual(len(bed), 3)
        self.assertEqual(bed.iloc[0]["CHROM"], "chr1")
        self.assertEqual(bed.iloc[0]["START"], 100)
        self.assertEqual(bed.iloc[1]["START"], 999)
        self.assertEqual(bed.iloc[2]["CHROM"], "chr2")


if __name__ == "__main__":
    unittest.main()
