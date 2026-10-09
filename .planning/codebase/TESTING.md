---
last_mapped_commit: e7b4ed7374de1197b1c7bcde526d062730c9bd43
last_mapped_at: 2026-10-09
---
# Testing Patterns

**Analysis Date:** 2026-10-09

## Test Framework

This repo has **three test layers**:

1. **Pipeline-level integration tests** — nf-test `nextflow_pipeline` (whole `main.nf` on SLURM)
2. **Module-level process tests** — nf-test `nextflow_process` (single local module; currently only `expand_regions`)
3. **Script-level unit tests** — Python `unittest` (stdlib, no pytest)

**Runner (pipeline):**
- nf-test >= 0.9.2 with plugin `nft-utils@0.0.3` (loaded in `nf-test.config`)
- Config: `nf-test.config` (root) + `tests/nextflow.config` (cluster/executor config)
- Suite: `tests/deepcsa.nf.test`
- Snapshot store: `tests/deepcsa.nf.test.snap`

**Runner (Python):**
- `unittest` via `python -m unittest` or `python -m pytest` (pytest-compatible; no pytest config exists)
- Tests live in `bin/test/`: `test_check_samplesheet.py`, `test_check_contamination.py`, `test_mask_matrix.py`, `test_plot_selectionsideplots.py`, `test_utils_filter.py`

**Run Commands:**

```bash

# Pipeline integration tests (REQUIRES SLURM cluster + Singularity — do not run locally)

nf-test test tests/deepcsa.nf.test                    # whole suite
nf-test test tests/deepcsa.nf.test --tag normal       # single test by tag
nf-test test tests/deepcsa.nf.test --tag omega

# Module-level process tests (single module, fast — NOT auto-discovered)

nf-test test modules/local/expand_regions/tests/main.nf.test

# Python unit tests (run anywhere)

python -m unittest discover -s bin/test -v
python -m unittest bin/test/test_utils_filter.py      # single file
python -m pytest bin/test/ -v                         # pytest alternative
```

**Work directory:** `DEEPCSA_TEST_WORKDIR` env var, default `.nf-test/` in repo root (`nf-test.config`). Test outputs live under `<workDir>/tests/<TEST_ID>/`.

## Test File Organization

**Location:**
- Pipeline tests: `tests/` (separate from code, committed)
- Python unit tests: `bin/test/` (separate `test/` subdir next to the scripts they test)

**Naming:**
- nf-test: one suite file `tests/deepcsa.nf.test` containing all pipeline tests
- Python: `test_<module_under_test>.py` mirroring the module name (`test_utils_filter.py` → `bin/utils_filter.py`)

**Structure:**

```
tests/
├── deepcsa.nf.test        # nf-test suite (5 active tests + commented templates)
├── deepcsa.nf.test.snap   # MD5 snapshots of pipeline outputs
├── nextflow.config        # SLURM executor, Singularity, cluster paths
├── test_data/             # committed local inputs (input.csv, input_maf.csv,
│                          #   input_no_bam.csv, test_mutations.maf)
└── README.md              # how to run, snapshot policy, debugging guide
modules/local/expand_regions/tests/
├── main.nf.test           # module-level nf-test (nextflow_process)
└── main.nf.test.snap      # MD5 snapshots of process outputs
test_data/                 # repo-root fixtures for module tests
├── dummy_file.tsv
└── modules/               # PPM1D BED/TSV fixtures for expand_regions
bin/test/
└── test_*.py              # unittest files
```

## Pipeline Tests (nf-test)

**Suite Organization** (`tests/deepcsa.nf.test`):

```groovy
nextflow_pipeline {
    name "Test DEEPCSA Pipeline"
    script "main.nf"
    config "./nextflow.config"

    test("TEST 1. Basic functionality - MAF-based processing") {
        tag "input_maf"

        when {
            params {
                input = "${projectDir}/tests/test_data/input_maf.csv"
                outdir = "$outputDir"   // nf-test built-in: unique temp dir per test
                input_maf = 'https://raw.githubusercontent.com/bbglab/DeepClone_protocol/main/...'
                use_custom_depths = true
                profileall = true
                signatures = false
            }
        }

        then {
            assert workflow.success : "Pipeline should complete without errors"
            assert path("${outputDir}/mutational_profile").exists()
            // Negative assertions: disabled features must NOT produce output
            assert !path("${outputDir}/mutdensity").exists()
            assert !path("${outputDir}/selection/omega").exists()
            // Snapshot of stable output
            assert snapshot(path("${outputDir}/mutational_profile/all_samples.all.profile.tsv")).match()
        }
    }
}
```

**Active tests (5):**

| Test | Tag | Purpose |
|------|-----|---------|
| TEST 1 | `input_maf` | MAF-based run, minimal features, snapshot of profile TSV |
| TEST 1b | `input_vcf_with_depths` | VCF + custom depths run |
| TEST 2 | `omega` | Omega analysis enabled; schema + structural + snapshot assertions |
| TEST 3 | `input_maf_nodepths_validation` | Pipeline **fails** when `--input_maf` without `--use_custom_depths` |
| TEST 4 | `input_csvnobam_nodepths_validation` | Fails: no BAMs, no depths |
| TEST 5 | `input_csv_no_bam_no_depthsfile_validation` | Fails: depths enabled but no file |

**Assertion patterns:**
- Success: `assert workflow.success : "message"`
- Failure tests: `assert workflow.failed : "message"` (validation tests assert the pipeline exits non-zero)
- Directory existence/non-existence: `path("${outputDir}/x").exists()`
- Snapshots: `snapshot(path(...)).match()` — MD5 of file content stored in `tests/deepcsa.nf.test.snap`
- Content schema checks: read file with `path(...).readLines()`, split on `\t`, assert header columns (`gene`, `sample`, `dnds`, `pvalue_adj`), row column-count consistency, and expected sample membership (TEST 2)
- Row-count assertions for non-deterministic float output: `assert dataLines.size() == 252`

**Handling non-deterministic outputs** (TEST 2 pattern — reuse this for new float-heavy outputs):

```groovy
// Sort by key columns to avoid floating-point order differences
def sortedRounded = filteredRows.sort { line -> /* key columns */ }
    .collect { line ->
        def cols = line.split('\t', -1)
        (4..7).each { i -> snapshotCols[i] = String.format("%.2f", snapshotCols[i] as Double) }
        snapshotCols.join('\t')
    }
assert snapshot(sortedContent.join('\n')).match("omega_results")
```

**Test config** (`tests/nextflow.config`):
- `executor = 'slurm'`, `errorStrategy = 'retry'`, `maxRetries = 2`
- Singularity with `cacheDir`/`libraryDir` on cluster storage
- `validation.ignoreParams = ['input_maf', 'custom_depths_table']` — skips nf-schema file-existence checks for remote HTTP test inputs
- Large remote test data fetched at runtime from `bbglab/DeepClone_protocol` GitHub repo (no local download step)

**Snapshot policy** (from `tests/README.md`):
- Regenerate only from the cluster, never locally: `nf-test test tests/deepcsa.nf.test --update-snapshot` (optionally with `--tag <tag>`)
- **Mandatory after changing default pipeline parameters**
- Review new hashes in `tests/deepcsa.nf.test.snap` before committing

**Debugging a failed pipeline test** (from `tests/README.md`):

```bash
cd <workDir>/tests/<TEST_ID>/work/<HASH>/<HASH>
cat .command.out .command.err .command.sh
bash .command.run          # reproduce exact environment
cat <workDir>/tests/<TEST_ID>/meta/nextflow.log
```

## Module Tests (nf-test `nextflow_process`)

**Current state:** exactly **one** module-level test exists — `modules/local/expand_regions/tests/main.nf.test` for the `EXPAND_REGIONS` process (`modules/local/expand_regions/main.nf`). The other ~45 local modules in `modules/local/` have no module tests.

**Test file convention** (follows the nf-core module test layout):

```
modules/local/<module_name>/
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test        # nextflow_process block
    └── main.nf.test.snap   # snapshots of process output channels
```

**Structure** (`modules/local/expand_regions/tests/main.nf.test`):

```groovy
nextflow_process {
    name "Test EXPAND_REGIONS process"
    script "modules/local/expand_regions/main.nf"   // path relative to repo root
    process "EXPAND_REGIONS"                        // process name inside main.nf

    test("Testing a run without autoexons and autodomains, it should fail") {
        when {
            params {
                autoexons = false
                autodomains = false
                subgenic_bedfile = false
            }
            process {
                """
                // Groovy heredoc: set each input channel of the process
                input[0] = tuple(
                    [ id:'test', single_end:false ],
                    file("${projectDir}/test_data/modules/consensus.exons_splice_sites.PPM1D.tsv")
                )
                input[1] = file("${projectDir}/test_data/modules/PPM1D_domains.bed4.bed")
                input[2] = file("${projectDir}/test_data/modules/PPM1D_exons.bed4.bed")
                input[3] = file("${projectDir}/test_data/dummy_file.tsv")
                """
            }
        }

        then {
            assert !process.success          // negative test: process must fail
        }
    }

    test("Should run with autoexons and autodomains") {
        when {
            params { autoexons = true; autodomains = true; subgenic_bedfile = false }
            process { """ input[0] = ...; input[1] = ... """ }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()   // snapshot ALL output channels
        }
    }
}
```

**Key differences from pipeline tests:**

| Aspect | `nextflow_pipeline` (tests/deepcsa.nf.test) | `nextflow_process` (module tests) |
|---|---|---|
| Scope | whole `main.nf`, all workflows | one process from one module file |
| Inputs | `params {}` (samplesheet, flags) | `process {}` heredoc assigning `input[N]` channels |
| Assertions | `workflow.success/failed`, published dirs | `process.success`, `process.out` channels |
| Snapshots | published output files | `snapshot(process.out).match()` — MD5 per channel element |
| Fixtures | `tests/test_data/` + remote HTTP URLs | repo-root `test_data/` (e.g. `test_data/modules/`) |
| Snapshot file | `tests/deepcsa.nf.test.snap` | `modules/local/<mod>/tests/main.nf.test.snap` |
| Runtime | SLURM + Singularity, minutes–hours | still uses `tests/nextflow.config` (SLURM), but seconds–minutes |

**Snapshot format** (`main.nf.test.snap`): JSON keyed by test name → `content` (per-channel MD5 hashes, e.g. `"panel_increased": [["...tsv:md5,6ad5da..."]]`) → `meta` (nf-test/Nextflow versions) → `timestamp`. The `meta` block records tool versions, so snapshots are regenerated when nf-test/Nextflow versions change.

**⚠️ Discovery caveat:** `nf-test.config` sets `testsDir "tests"` and `ignore 'modules/nf-core/**/*', 'subworkflows/nf-core/**/*'`. This means:
- Running bare `nf-test test` discovers only `tests/deepcsa.nf.test` — module tests under `modules/local/**/tests/` are **not** auto-discovered and must be run by explicit path.
- The `ignore` patterns exclude nf-core modules/subworkflows from testing (they have their own CI upstream), but local modules are *not* ignored — they're just outside `testsDir`.
- Module tests still load `tests/nextflow.config` (via `configFile` in `nf-test.config`), so they submit to SLURM and use Singularity like pipeline tests. They are not runnable on a laptop.

**Stub-run pattern:** the expand_regions test file contains a commented-out `stub = true` test template (`config { stub = true }` inside `then {}`) — the standard nf-core pattern for testing module stub blocks without executing the real script. Revive it when `main.nf` gains a `stub:` block.

**Adding a module test (checklist):**
1. Create `modules/local/<module>/tests/main.nf.test` with a `nextflow_process` block; `script` path is relative to repo root, `process` is the process name in `main.nf`
2. Add small fixtures under repo-root `test_data/modules/` (keep them tiny — they're committed)
3. In `when { process { """...""" } }`, assign every declared input channel (`input[0]`, `input[1]`, ...) using `file("${projectDir}/...")`
4. Assert `process.success` (or `!process.success` for expected-failure cases) and `snapshot(process.out).match()`
5. Generate the snapshot: `nf-test test modules/local/<module>/tests/main.nf.test --update-snapshot` (from the cluster)
6. Commit both `main.nf.test` and `main.nf.test.snap`

## Python Unit Tests (unittest)

**Suite Organization** (`bin/test/test_utils_filter.py`):

```python
#!/usr/bin/env python3
"""Module docstring listing what is covered."""

import sys
import tempfile
import unittest
from pathlib import Path
import pandas as pd

# Add the bin directory to the path to import sibling modules

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils_filter import filter_maf, somatic_mask

THRESHOLD = 0.3

# ---------------------------------------------------------------------------

# Helpers

# ---------------------------------------------------------------------------

def _make_df(rows):
    """Build a minimal MAF DataFrame."""
    return pd.DataFrame(...)

class TestSomaticMask(unittest.TestCase):
    """Tests for somatic_mask(maf_df, threshold)."""

    def test_all_below_threshold_is_somatic(self):
        """All three VAF columns strictly below threshold → somatic True."""
        ...
```

**Patterns:**
- **Path bootstrap (required):** every test file starts with `sys.path.insert(0, str(Path(__file__).parent.parent))` because `bin/` scripts import each other as top-level modules
- **Class per function under test:** `TestSomaticMask`, `TestFilterMaf` — one `TestCase` class per function, named `Test<FunctionName>`
- **Docstring = assertion explanation:** each test method's docstring states the expected behavior in plain language
- **setUp/tearDown with temp dirs:** `tempfile.mkdtemp()` in `setUp`, `shutil.rmtree` + `os.chdir` restore in `tearDown` (`bin/test/test_mask_matrix.py`)
- **Fixture-builder helpers:** module-level `_make_df()` / `_make_maf()` functions and instance methods like `create_mock_bed_file(sample_name, positions)` that write small synthetic files
- **Section comment banners** (`# ---- basic behaviour ----`) to group related tests
- **SystemExit handling:** for CLI scripts that call `sys.exit`, tests wrap calls in `try/except SystemExit` and `self.fail(...)` on unexpected exits (`bin/test/test_check_samplesheet.py`)

**What is tested (113 test methods/classes across 5 files):**
- `check_samplesheet.py` — sample-name validation (security: shell-injection and flag-injection prevention)
- `check_contamination.py` — `compute_shared_variants` on minimal synthetic DataFrames
- `create_mask_matrix.py` + `merge_annotation_depths.apply_mask_matrix` — position × sample mask matrices from synthetic BED files
- `utils_filter.py` — somatic/germline masks, `filter_maf` branches, criteria parsing, BED extraction
- `plot_selectionsideplots.py` — plotting helpers

## Mocking

**Framework:** none — `unittest.mock` is **not used** anywhere in `bin/test/`.

**Patterns:**
- Instead of mocking, tests build **real minimal in-memory fixtures** (small pandas DataFrames) or **real tiny files** in temp dirs (synthetic BED/CSV files)
- The word "mock" appears only in helper names like `create_mock_bed_file` — these create genuine files, not mock objects

**What to Mock:** nothing currently; follow the existing approach — construct minimal real inputs.

**What NOT to Mock:** pandas DataFrames and small text files are cheap to create for real; do not introduce `MagicMock` for them.

## Fixtures and Factories

**Test Data:**

```python

# In-memory DataFrame factory (bin/test/test_utils_filter.py)

def _make_df(rows: list[tuple[float, float, float]]) -> pd.DataFrame:
    vafs, vd_vafs, vaf_ams = zip(*rows)
    return pd.DataFrame({"VAF": list(vafs), "vd_VAF": list(vd_vafs), "VAF_AM": list(vaf_ams)})

# On-disk fixture writer (bin/test/test_mask_matrix.py)

def create_mock_bed_file(self, sample_name, positions):
    bed_file = f"{sample_name}.flagged-pos.bed"
    with open(bed_file, 'w') as f:
        for chrom, start, end, filter_val in positions:
            f.write(f"{chrom}\t{start}\t{end}\t{filter_val}\n")
```

**Location:**
- Python: fixtures are generated inline in test files (no `fixtures/` directory)
- Pipeline: committed small inputs in `tests/test_data/` (`input.csv`, `input_maf.csv`, `input_no_bam.csv`, `test_mutations.maf`); large cohort data fetched at runtime from `https://raw.githubusercontent.com/bbglab/DeepClone_protocol/main/test_datasets/deepCSA/...`

## Coverage

**Requirements:** None enforced. No coverage tooling configured (no `pytest.ini`, `tox.ini`, `codecov.yml`, no CI workflows — `.github/` does not exist).

**View Coverage:**

```bash
python -m pytest bin/test/ --cov=bin --cov-report=term-missing   # if pytest-cov installed
```

## Test Types

**Unit Tests:**
- Python `unittest` in `bin/test/` — pure-function tests on synthetic DataFrames/files, no network, no cluster

**Integration Tests:**
- nf-test pipeline runs in `tests/deepcsa.nf.test` — full `main.nf` execution on SLURM with Singularity containers, validating published output trees, file schemas, row counts, and MD5 snapshots

**E2E Tests:**
- The nf-test suite *is* the E2E layer (entire pipeline end-to-end). No browser/UI tests (not applicable).

## Common Patterns

**Async Testing:** Not applicable (no async code).

**Error Testing:**

```python

# Python: CLI scripts exit non-zero; tests catch SystemExit (bin/test/test_check_samplesheet.py)

try:
    check_samplesheet(input_file, output_file)
    self.fail("Invalid sample name was accepted")   # inverted for invalid-input cases
except SystemExit:
    pass  # expected
```

```groovy
// nf-test: assert the whole pipeline fails on invalid params (tests/deepcsa.nf.test, TEST 3)
then {
    assert workflow.failed : "Pipeline should fail when --input_maf is set without --use_custom_depths"
}
```

**Adding a new pipeline test (checklist):**
1. Add a `test("TEST N. <description>")` block to `tests/deepcsa.nf.test` with a unique `tag`
2. Set params in `when { params { ... } }`; use `$outputDir` for `outdir`
3. If params point to remote HTTP files, add them to `validation.ignoreParams` in `tests/nextflow.config`
4. Assert success/failure, directory presence/absence, then snapshot only deterministic outputs (round floats, sort rows first)
5. Run from the cluster; update snapshots with `--update-snapshot` if outputs changed intentionally

**Adding a new Python unit test (checklist):**
1. Create `bin/test/test_<module>.py` with the `sys.path.insert` bootstrap
2. One `TestCase` class per function under test; docstring each test with the expected behavior
3. Build inputs with `_make_*` helper factories or temp-dir files; clean up in `tearDown`
4. Run: `python -m unittest discover -s bin/test -v`

## Gaps to Be Aware Of

- **No CI:** `.github/` does not exist — neither test layer runs automatically on push/PR
- **~71 of 80 `bin/` scripts have no unit tests** — only 5 modules are covered
- **~45 of 46 local Nextflow modules have no module tests** — only `expand_regions` has a `nextflow_process` test; module tests are also not auto-discovered by `nf-test.config` (`testsDir "tests"`)
- **No pytest config / coverage gate** — coverage is unmeasured
- **nf-test requires SLURM + Singularity** — cannot run the integration suite on a laptop or in plain CI without a cluster runner
- **Commented-out test templates** at the bottom of `tests/deepcsa.nf.test` (MAF + precomputed depths integration test, trace-based process-count assertions) are ready-made patterns to revive when test assets improve

---

*Testing analysis: 2026-10-09*
