#!/usr/bin/env python3

import os
import re
import unittest
from pathlib import Path


DEFAULT_RUN_DIR = Path("scratchhhh/2026-04-17_add_tests")
TRACE_GLOB = "pipeline_info/execution_trace_*.txt"


def _get_reference_run_dir() -> Path:
    env_value = os.environ.get("DEEPCSA_REFERENCE_RUN_DIR")
    if env_value:
        return Path(env_value)
    return DEFAULT_RUN_DIR


def _extract_families(trace_file: Path) -> set[str]:
    families: set[str] = set()
    lines = trace_file.read_text(encoding="utf-8").splitlines()
    for line in lines[1:]:
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        full_name = cols[3].split(" (")[0]
        families.add(full_name.split(":")[-1])
    return families


def _classify_family(family: str) -> str | None:
    if family in {"SAMPLESHEET_CHECK", "CUSTOM_DUMPSOFTWAREVERSIONS", "MULTIQC"}:
        return "pipeline_info"

    if family in {
        "TABLE2GROUP",
        "FILTERDEPTHS",
        "SITESFROMPOSITIONS",
        "CREATECAPTUREDPANELS",
        "CREATECONSENSUSPANELSALL",
        "CREATECONSENSUSPANELSEXONS",
        "CREATECONSENSUSPANELSINTRONS",
        "CREATECONSENSUSPANELSNONPROTAFFECT",
        "CREATECONSENSUSPANELSPROTAFFECT",
        "CREATECONSENSUSPANELSSYNONYMOUS",
        "SORTPANELCOMPACT",
        "SORTPANELRICH",
        "CUSTOMPROCESSING",
        "CUSTOMPROCESSINGRICH",
        "DOMAINANNOTATION",
        "CUSTOMBEDFILE",
        "COMPAREPANELPROPORTIONS",
        "EXPANDREGIONSALL",
        "EXPANDREGIONSEXONS",
        "EXPANDREGIONSNONPROT",
        "EXPANDREGIONSPROT",
        "EXPANDREGIONSSYNONYMOUS",
        "DNA2PROTEINMAPPING",
        "GROUPGENES",
    }:
        return "regions"

    if family in {
        "ENSEMBLVEP_VEP",
        "POSTPROCESSVEPPANEL",
        "CUSTOMANNOTATION",
        "SUMANNOTATION",
        "VCF2MAF",
        "FILTEREXONS",
        "FILTERPANEL",
        "FILTERNANOSEQSNP",
        "FILTERNANOSEQNOISE",
        "MERGEBATCH",
        "FILTERBATCH",
        "WRITEMAF",
        "CLEANMUTATIONS",
        "SOMATICMUTATIONS",
        "MAF2VCF",
        "CREATEMASKMATRIX",
        "ANNOTATEDEPTHS",
        "QUERYMUTATIONSEXONS",
        "PLOTSOMATICMAF",
        "PLOTMAF",
    }:
        return "mutations"

    if family in {
        "DEPTHSALLCONS",
        "DEPTHSEXONSCONS",
        "DEPTHSNONPROTCONS",
        "DEPTHSPROTCONS",
        "DEPTHSSYNONYMOUSCONS",
        "QUERYDEPTHS",
        "DEPTHSSUMMARY",
        "MUTDENSITY",
        "SUBSETMUTDENSITY",
        "MUTDENSITYADJ",
        "SUBSETMUTDENSITYADJUSTED",
        "SYNMUTDENSITY",
        "SYNMUTREADSDENSITY",
        "PLOTMUTDENSITYQC",
    }:
        return "mutdensity"

    if family in {
        "QUERYMUTATIONS",
        "SUBSETMUTPROFILE",
        "COMPUTETRINUC",
        "COMPUTEMATRIX",
        "COMPUTEPROFILE",
        "CONCATPROFILES",
        "MATRIXCONCATWGS",
        "SIGPROMATRIXGENERATOR",
        "SIGPROFILERASSIGNMENT",
        "SIGPROFILERASSIGNMENTINDELS",
        "SIGPROBS",
        "MUTS2SIGS",
        "COMPARESIGNATURES",
        "PREPARE_INPUT",
        "RUN_HDP_CHAIN_SAMPLING",
        "PROCESS_HDP_RESULTS",
        "REFORMATSIGNATURES",
        "HDPREASSIGNMENT",
    }:
        return "mutational_profile"

    if family in {
        "PREPROCESSDEPTHS",
        "DNDSRUN",
        "SUBSET_DNDS",
        "RELATIVEMUTABILITY",
        "ABSOLUTEMUTABILITIES",
        "ABSOLUTEMUTABILITIESGLOBALLOC",
        "MUTABILITY_BGZIPTABIX",
        "SUBSETMUTABILITY",
        "ONCODRIVEFMLBED",
        "ONCODRIVEFMLSNVS",
        "SUBSETONCODRIVEFML",
        "ONCODRIVE3D_PREPROCESSING",
        "ONCODRIVE3D_RUN",
        "ONCODRIVE3D_PLOT",
        "SUBSETONCODRIVE3D",
        "PREPROCESSING",
        "PREPROCESSINGGLOBALLOC",
        "ESTIMATOR",
        "ESTIMATORGLOBALLOC",
        "SUBSETOMEGA",
        "SUBSETOMEGAMULTI",
        "QUERYPANEL",
        "SITECOMPARISON",
        "SITECOMPARISONMULTI",
        "SITECOMPARISONGLOBALLOC",
        "SITECOMPARISONGLOBALLOCMULTI",
        "EVALOMEGAGLOCESTIMATION",
        "APPLYOMEGAQC",
        "PLOTOMEGA",
        "PLOTOMEGAGLOBALLOC",
        "PLOTSELECTION",
        "INDELS",
        "SUBSETINDELS",
    }:
        return "selection"

    if family in {
        "PLOTMUTATIONSPECIFIC",
        "PLOTNEEDLES",
        "PLOTSATURATION",
        "PLOTSATURATIONPROPORTIONS",
        "PLOTINTERINDIVIDUALVARIABILITY",
    }:
        return "plots"

    return None


class TestExecutedStepCoverage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.run_dir = _get_reference_run_dir()
        if not cls.run_dir.exists():
            raise unittest.SkipTest(
                f"Reference run directory not found: {cls.run_dir}. "
                "Set DEEPCSA_REFERENCE_RUN_DIR to the location of the run outputs."
            )

        trace_files = sorted(cls.run_dir.glob(TRACE_GLOB))
        if not trace_files:
            raise unittest.SkipTest(f"No execution trace found under {cls.run_dir / 'pipeline_info'}")

        cls.trace_file = trace_files[0]
        cls.families = _extract_families(cls.trace_file)

    def test_reference_trace_has_many_families(self):
        # The April 17, 2026 reference run executes a broad workflow footprint.
        self.assertGreaterEqual(len(self.families), 100)

    def test_every_family_is_mapped_to_a_test_category(self):
        unmapped = sorted(f for f in self.families if _classify_family(f) is None)
        self.assertEqual(
            unmapped,
            [],
            msg=(
                "The following executed process families are not yet mapped to a coverage category: "
                + ", ".join(unmapped)
            ),
        )

    def test_each_coverage_category_has_output_evidence(self):
        evidence_patterns = {
            "pipeline_info": [
                "pipeline_info/execution_trace_*.txt",
                "pipeline_info/software_versions.yml",
            ],
            "regions": [
                "regions/**/*.tsv",
                "regions/**/*.bed",
                "group_definition/**/*.json",
                "processing_files/**/create*/*",
            ],
            "mutations": [
                "mutations/**/*.tsv",
                "mutations/**/*.maf",
                "processing_files/**/*maf*",
            ],
            "mutdensity": [
                "mutdensity/**/*.tsv",
                "mutdensity_adjusted/**/*.tsv",
                "depths/**/*.tsv*",
            ],
            "mutational_profile": [
                "mutational_profile/**/*.tsv",
                "signatures/**/*",
            ],
            "selection": [
                "selection/**/*.tsv",
                "selection/**/*.png",
                "selection/**/*.pdf",
            ],
            "plots": [
                "plots/**/*.png",
                "plots/**/*.pdf",
            ],
        }

        observed_categories = {_classify_family(f) for f in self.families if _classify_family(f) is not None}

        for category in sorted(observed_categories):
            patterns = evidence_patterns[category]
            found_any = False
            for pattern in patterns:
                if any(self.run_dir.glob(pattern)):
                    found_any = True
                    break
            self.assertTrue(
                found_any,
                msg=f"No output evidence found for category '{category}' under {self.run_dir}",
            )


if __name__ == "__main__":
    unittest.main()
