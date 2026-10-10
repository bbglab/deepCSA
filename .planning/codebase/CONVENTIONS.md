---
last_mapped_commit: 3c83c5b8aa313077e0ce43239a7b818b281b32ea
last_mapped_at: 2026-10-09
---
# Coding Conventions

**Analysis Date:** 2026-10-09

## Language Mix

The codebase is a Nextflow (DSL2) pipeline with three script languages in `bin/`:

- **Python** (~76 scripts) — dominant language for analysis/plotting logic
- **R** (4 scripts) — `bin/dNdS_run.R`, `bin/mutgenomes_expected_mutrisk.R`, `bin/mutrate_genome_trinuc_corrected.R`, `bin/signatures_msighdp_run.R`
- **Groovy/Nextflow** — `main.nf`, `workflows/`, `subworkflows/local/`, `modules/local/`
- **Shell** (1 script) — `bin/createcustombed.sh`

## Naming Patterns

**Files (Python):**
- `snake_case.py` for all scripts: `check_samplesheet.py`, `mut_profile.py`, `plot_depths.py`
- Shared helper modules prefixed `utils_` or suffixed `_utils`: `bin/read_utils.py`, `bin/utils.py`, `bin/utils_plot.py`, `bin/utils_filter.py`, `bin/utils_context.py`, `bin/utils_impacts.py`
- Test files: `test_<module>.py` in `bin/test/`

**Files (Nextflow):**
- Processes live in `modules/local/<snake_case_process>/main.nf` (e.g. `modules/local/compute_profile/main.nf`)
- Subworkflows live in `subworkflows/local/<snake_case_group>/main.nf` (e.g. `subworkflows/local/omega/main.nf`); multi-file groups use `subworkflows/local/<group>/main` as entry
- One workflow file: `workflows/deepcsa.nf`

**Functions:**
- `snake_case`: `compute_mutation_matrix()`, `validate_and_transform()`, `check_samplesheet()`
- Private methods prefixed `_`: `RowChecker._validate_sample()`, `_validate_vcf_format()` (`bin/check_samplesheet.py`)
- Test helpers prefixed `_`: `_make_df()`, `_make_maf()` (`bin/test/test_utils_filter.py`)

**Classes:**
- `PascalCase`: `RowChecker` (`bin/check_samplesheet.py`), `TestSampleNameValidation` (`bin/test/test_check_samplesheet.py`)

**Variables/Constants:**
- `snake_case` locals and module-level config: `custom_na_values` (`bin/read_utils.py`)
- `UPPER_SNAKE_CASE` for module constants: `MAX_CATEGORIES_PER_PLOT`, `MAX_HEATMAP_CELLS` (`bin/plot_depths.py`), `MIN_NONZERO_PVALUE` (`bin/utils.py`), `VALID_FORMATS_BAM`/`VALID_FORMATS_VCF` (`bin/check_samplesheet.py`)

**Nextflow processes:**
- `UPPERCASE` process names: `COMPUTE_PROFILE`, `INPUT_CHECK`
- Subworkflows `UPPERCASE` with aliases on import: `MUTATIONAL_PROFILE as MUTPROFILEALL` (`workflows/deepcsa.nf`)

## Code Style

**Formatting:**
- Black configured in `pyproject.toml`: `line-length = 120`, `target_version = ["py37", "py38", "py39", "py310"]`
- isort with `profile = "black"` (`pyproject.toml`)
- `.editorconfig`: 4-space indent, LF line endings, UTF-8, final newline, trimmed trailing whitespace; 2-space for `*.md`, `*.yml`, `*.yaml`, `*.html`, `*.css`, `*.scss`, `*.js`
- **Reality check:** only `bin/check_samplesheet.py` is consistently Black-formatted (it is nf-core template code). Most other `bin/` scripts predate the config and use pandas-style `=` spacing (`sep = "\t"`), inconsistent quoting, and occasional 2-space indents. **Do not reformat existing scripts wholesale; match the style of the file you are editing. New nf-core-template-derived code must be Black-formatted.**

**Linting:**
- No repo-level flake8/ruff/pylint config. `.devcontainer/devcontainer.json` enables pylint/flake8 in the dev container only.
- `markdownlint` config in `.markdownlint.json`: `MD013` (line length) off, `MD024` siblings-only.

## Shebangs & Execution Model

- Python: `#!/usr/bin/env python3` (nf-core template scripts) or `#!/usr/bin/env python` (older scripts)
- R: `#!/opt/conda/bin/Rscript --vanilla` (hard-coded conda path — matches the `bbglab/deepcsa-core` container)
- Scripts in `bin/` are invoked by Nextflow processes directly by name (Nextflow auto-stages `bin/` onto `PATH`), e.g. `mut_profile.py profile --sample_name ...` in `modules/local/compute_profile/main.nf`

## CLI Argument Parsing

**Two coexisting styles — prefer `click` for new scripts:**

1. **click (dominant, newer scripts):** decorators with `click.Choice`, `click.Path(exists=True)` for input validation, `is_flag=True` for booleans. Example: `bin/mut_profile.py` (`@click.command()`, `@click.argument('mode', type=click.Choice(['matrix', 'profile']))`), `bin/plot_depths.py`
2. **argparse (nf-core template scripts):** `parse_args(argv=None)` function + `main(argv=None)` + `if __name__ == "__main__: sys.exit(main())`. Example: `bin/check_samplesheet.py`

**Rules for new scripts:**
- Use `click.Path(exists=True)` so missing inputs fail fast
- Use `click.Choice` for enumerated modes
- Keep a `main()` entry point guarded by `if __name__ == '__main__':` (71 of 80 Python scripts do this) so functions stay importable by tests

## Import Organization

**Order (observed in `bin/check_samplesheet.py`, `bin/mut_profile.py`):**
1. Standard library (`sys`, `argparse`, `csv`, `logging`, `re`, `pathlib`)
2. Third-party (`click`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `polars`)
3. Local sibling modules (`from utils import contexts_formatted`, `from read_utils import custom_na_values`)

**Path Aliases:**
- None. Sibling imports rely on Nextflow staging all of `bin/` into the task working directory, so `import utils` works at runtime. Unit tests replicate this with `sys.path.insert(0, str(Path(__file__).parent.parent))` (`bin/test/test_utils_filter.py`).

## Error Handling

**Patterns:**
- **Validation scripts (nf-core template style):** raise `AssertionError` with descriptive f-string messages from validator methods; catch in the caller, log with `logger.critical(...)`, and `sys.exit(1)` (`bin/check_samplesheet.py`). Missing input file → `sys.exit(2)`.
- **click scripts:** rely on `click.Path(exists=True)` for input validation; domain errors often just `print(...)` + `exit(1)` (e.g. `bin/mut_profile.py` `compute_mutation_matrix`).
- **Fail-fast over recovery:** pipeline scripts generally do not catch exceptions; Nextflow's `errorStrategy = 'retry'` with `maxRetries = 2` (`tests/nextflow.config`) handles transient failures at the process level.
- **Defensive data handling:** explicit `fillna(0)` after reindexing, percentile clipping of outliers (99.5th percentile of `ALT_DEPTH` in `bin/mut_profile.py`), and `custom_na_values` lists passed to `pd.read_csv` (`bin/read_utils.py`).

**Do this for new validation code:** raise `AssertionError`/`ValueError` with a message that names the offending value and the allowed set, and let the CLI wrapper convert to a non-zero exit.

## Logging

**Framework:** mixed — no single standard.

- `logging` module with `logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")` and a module-level `logger = logging.getLogger()` — `bin/check_samplesheet.py` (nf-core template pattern)
- `click.echo()` for user-facing progress in click scripts (`bin/mut_profile.py`)
- Plain `print()` in most analysis scripts, sometimes with a `[scriptname]` prefix: `print(f"[plot_depths] Skipping plot section '{section}': {reason}")` (`bin/plot_depths.py`)

**When to log:** skip conditions with reasons, percentile/clip decisions, mode and parameter echo at startup. Keep messages one-line; Nextflow captures stdout/stderr per task.

## Comments

**When to Comment:**
- Section banners in test/Nextflow files using `/* ==== ... ==== */` blocks (`tests/deepcsa.nf.test`, `main.nf`)
- `# TODO` / `# FIXME` markers are used liberally (~20 occurrences across `bin/`), e.g. `bin/concat_sbs_probs.py:3` ("TODO: bump pandas to 2.2.3"), `bin/postprocessing_annotation.py:146`
- Explanatory comments for non-obvious numeric choices (plot size limits in `bin/plot_depths.py`)

**Docstrings:**
- **Google style** in nf-core template code: `Args:`, `Returns:`, `Attributes:`, `Raises:` sections, plus `Example:` blocks (`bin/check_samplesheet.py`)
- **NumPy style** in bbglab-written code: `Parameters\n----------`, `Returns\n-------` (`bin/utils.py`, `bin/test/test_mask_matrix.py`)
- Module-level docstrings summarizing purpose (`"""Provide a command line tool to validate and transform tabular samplesheets."""`)
- **For new code:** match the host file's style; Google style for anything derived from nf-core templates.

## Function Design

**Size:** no enforced limit; legacy scripts contain very long functions and files (`bin/plot_gene_saturation.py` is 1546 lines). New code should keep functions single-purpose.

**Parameters:** plain positional args for 2–3 params; keyword args with defaults for options (`def check_samplesheet(file_in, file_out, bam_required=False)`). R scripts use `optparse` long flags (`--inputfile`, `--outputprefix`).

**Return Values:** DataFrames in/out for data-transform functions; `None` + file writes for plot/report functions. Functions that can fail return `None` and let the caller decide (`compute_mutation_profile` returning `None` in `bin/mut_profile.py`).

## Module Design

**Exports:** flat modules; no `__all__`. Tests import named functions directly: `from utils_filter import filter_maf, somatic_mask` (`bin/test/test_utils_filter.py`).

**Barrel Files:** none. Shared constants live in dedicated modules (`bin/read_utils.py` for NA values / MAF reading, `bin/utils.py` for MAF filters and variant typing, `bin/utils_plot.py` for plotting helpers).

## Nextflow-Specific Conventions

**Process definition** (`modules/local/compute_profile/main.nf`):
- `tag "$meta.id"` for traceability
- Resource `label`s: `cpu_low`/`mem_low`/`process_high_memory` plus domain label `deepcsa_core` (container pinned in `conf/modules.config` via `withLabel: deepcsa_core { container = "docker.io/bbglab/deepcsa-core:0.1.0" }`)
- `input:`/`output:`/`script:`/`stub:` blocks; outputs use `emit:` names and `optional:true` where conditional
- `versions.yml` written via heredoc `cat <<-END_VERSIONS` and emitted `topic: versions` in every process
- `task.ext.args` / `task.ext.prefix` consumed with `?: ""` defaults for configurability from `conf/modules.config`
- `stub:` blocks present in 97 of 48 module dirs' main.nf files (widespread) for fast testing

**Config layering:** `nextflow.config` → `conf/base.config`, `conf/modules.config`, `conf/results_outputs.config`, tool configs (`conf/tools/*.config`), profile configs (`conf/test.config`, `conf/exome.config`, `conf/mice.config`, ...). Site-specific paths isolated in `conf/general_files_IRB.config`.

**Publishing:** output routing centralized in `conf/results_outputs.config` with `saveAs` filters that drop `versions.yml` from published dirs.

## R Conventions

- Shebang `#!/opt/conda/bin/Rscript --vanilla`
- `optparse` for CLI (`make_option(c("-n", "--samplename"), ...)`) — `bin/dNdS_run.R`
- `snake_case` function names (`is_SNV` is a legacy exception)
- Usage example in a header comment block

---

*Convention analysis: 2026-10-09*
