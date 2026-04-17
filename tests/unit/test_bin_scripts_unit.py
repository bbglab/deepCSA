#!/usr/bin/env python3

import csv
import gzip
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
BIN_DIR = REPO_ROOT / "bin"


class TestBinScripts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.has_pandas = importlib.util.find_spec("pandas") is not None

    def _run(self, args, cwd=None):
        return subprocess.run(
            [sys.executable] + args,
            cwd=cwd if cwd is not None else REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

    def test_check_samplesheet_valid_and_duplicate_renaming(self):
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            in_csv = td_path / "samples.csv"
            out_csv = td_path / "samples.valid.csv"

            rows = [
                ["sample", "vcf", "bam"],
                ["S1", "a.vcf.gz", "a.bam"],
                ["S1", "b.vcf.gz", "b.bam"],
            ]
            with in_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerows(rows)

            result = self._run([
                str(BIN_DIR / "check_samplesheet.py"),
                str(in_csv),
                str(out_csv),
                "--bam-required",
            ])
            self.assertEqual(result.returncode, 0, msg=result.stderr)

            with out_csv.open("r", newline="", encoding="utf-8") as handle:
                out_rows = list(csv.DictReader(handle))

            self.assertEqual(out_rows[0]["sample"], "S1")
            self.assertEqual(out_rows[1]["sample"], "S1_T2")

    def test_check_samplesheet_rejects_unsafe_sample_name(self):
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            in_csv = td_path / "samples_bad.csv"
            out_csv = td_path / "samples.valid.csv"

            rows = [
                ["sample", "vcf"],
                ["bad;name", "a.vcf.gz"],
            ]
            with in_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerows(rows)

            result = self._run([
                str(BIN_DIR / "check_samplesheet.py"),
                str(in_csv),
                str(out_csv),
            ])

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("invalid characters", (result.stdout + result.stderr).lower())

    def test_features_table_to_groups_generates_expected_jsons(self):
        if not self.has_pandas:
            self.skipTest("pandas is required for features_1table2groups.py")

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            table_file = td_path / "features.csv"
            table_file.write_text(
                "sample,group,sex\n"
                "S1,A,F\n"
                "S2,A,M\n"
                "S3,B,F\n",
                encoding="utf-8",
            )

            result = self._run(
                [
                    str(BIN_DIR / "features_1table2groups.py"),
                    "--table-filename",
                    str(table_file),
                    "--separator",
                    "comma",
                    "--unique-identifier",
                    "sample",
                    "--groups",
                    "[group],[sex]",
                ],
                cwd=td_path,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)

            samples_json = json.loads((td_path / "samples.json").read_text(encoding="utf-8"))
            groups_json = json.loads((td_path / "groups.json").read_text(encoding="utf-8"))
            all_groups_json = json.loads((td_path / "all_groups.json").read_text(encoding="utf-8"))

            self.assertEqual(samples_json["S1"], ["S1"])
            self.assertIn("all_samples", groups_json)
            self.assertCountEqual(groups_json["all_samples"], ["S1", "S2", "S3"])
            self.assertIn("GroupA", groups_json)
            self.assertCountEqual(groups_json["GroupA"], ["S1", "S2"])
            self.assertIn("SexF", groups_json)
            self.assertIn("S3", all_groups_json)

    def test_merge_cohort_merges_all_tsv_gz_files(self):
        if not self.has_pandas:
            self.skipTest("pandas is required for merge_cohort.py")

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)

            f1 = td_path / "a.tsv.gz"
            f2 = td_path / "b.tsv.gz"
            merged = td_path / "merged.tsv"

            with gzip.open(f1, "wt", encoding="utf-8") as handle:
                handle.write("sample\tvalue\nS1\t1\n")
            with gzip.open(f2, "wt", encoding="utf-8") as handle:
                handle.write("sample\tvalue\nS2\t2\n")

            result = self._run(
                [
                    str(BIN_DIR / "merge_cohort.py"),
                    "--output_file",
                    str(merged),
                ],
                cwd=td_path,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)

            content = merged.read_text(encoding="utf-8")
            self.assertIn("sample\tvalue", content)
            self.assertIn("S1\t1", content)
            self.assertIn("S2\t2", content)


if __name__ == "__main__":
    unittest.main()
