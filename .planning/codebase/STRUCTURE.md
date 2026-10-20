---
last_mapped_commit: 3c83c5b8aa313077e0ce43239a7b818b281b32ea
last_mapped_at: 2026-10-09
---
# Codebase Structure

**Analysis Date:** 2026-10-09

## Directory Layout

```
deepCSA/
├── main.nf                  # Entry point (BBGTOOLS wrapper workflow)
├── nextflow.config          # Params defaults, profiles, registries
├── nextflow_schema.json     # nf-schema param validation
├── modules.json             # Pinned nf-core module/subworkflow versions
├── nf-test.config           # nf-test runner config
├── tower.yml                # Seqera Platform report display config
├── pyproject.toml           # Black/isort config for bin/ scripts
├── workflows/               # Pipeline workflow layer
│   └── deepcsa.nf           # THE pipeline (workflow DEEPCSA, ~800 lines)
├── subworkflows/
│   ├── local/               # 23 local subworkflows (analysis stages)
│   │   ├── depthanalysis/   # Depth computation + filtering
│   │   ├── createpanels/    # Consensus panel creation
│   │   ├── mutationpreprocessing/  # VEP annotation, MAF building
│   │   ├── mutationdensity/ mutationprofile/ mutability/
│   │   ├── omega/ oncodrivefml/ oncodrive3d/ oncodriveclustl/ dnds/
│   │   ├── signatures/ signatures_hdp/ mutatedcells/ indels/
│   │   ├── regressions/ enrichpanels/ adjmutdensity/
│   │   ├── plotdepths/ plotting_qc/ plottingsummary/
│   │   └── utils_nfcore_deepcsa/   # Init: banner, versions, methods text
│   └── nf-core/             # Vendored nf-core subworkflows (5)
├── modules/
│   ├── local/               # 97 local process modules (48 top-level dirs)
│   │   ├── bbgtools/        # bbglab tool wrappers: omega/, oncodrive3d/,
│   │   │                    #   oncodrivefml/, oncodriveclustl/,
│   │   │                    #   bbgregressions/, sitecomparison/
│   │   ├── createpanels/    # captured/, consensus/, compare/, custombedfile/
│   │   ├── dnds/ signatures/ downsample/ plot/ ...
│   │   └── <tool>/main.nf   # One process per main.nf
│   └── nf-core/             # Vendored nf-core modules (multiqc, dumpsoftwareversions)
├── bin/                     # ~90 analysis scripts (auto on PATH in tasks)
│   ├── *.py                 # Python CLIs (click + pandas)
│   ├── *.R                  # R scripts (dNdS_run.R, mutrate_genome_trinuc_corrected.R)
│   ├── utils*.py            # Shared helpers (utils.py, utils_plot.py, read_utils.py...)
│   ├── test/                # pytest unit tests (test_*.py)
│   └── saturation_mutagenesis/  # Saturation kinetics scripts + notebooks
├── conf/                    # Nextflow config includes
│   ├── base.config          # Resource limits, labels, error strategy
│   ├── modules.config       # Per-process ext.args + publishDir (78 withName blocks)
│   ├── results_outputs.config  # Curated final publish paths
│   ├── test.config, test_real.config, exome.config, mice.config, local.config
│   ├── modes/               # Analysis presets: basic, clonal_structure, get_signatures
│   └── tools/               # Per-tool configs: omega, oncodrive3d, oncodrivefml,
│                            #   mutdensity, panels, regressions, hdp_sig_extraction
├── assets/                  # Static reference data & templates
│   ├── schema_input.json    # Samplesheet JSON schema
│   ├── multiqc_config.yml, email_template.*, sendmail_template.txt, slackreport.json
│   ├── omega_consequences_groupings.json, chromosome_bands/, trinucleotide_counts/
│   ├── omega-covariates/, regressions/, build_datasets/, assess_panel/
│   ├── example_inputs/      # Example samplesheets/features tables
│   └── useful_scripts/      # Standalone helper scripts/notebooks (not run by pipeline)
├── docs/                    # usage.md, output.md, metrics.md, tools.md, input_scenarios.md...
├── tests/                   # nf-test pipeline tests
│   ├── deepcsa.nf.test      # Main pipeline test
│   ├── deepcsa.nf.test.snap # Snapshot file
│   ├── nextflow.config      # Test-specific config
│   └── test_data/           # Test fixtures (partially gitignored)
├── test_data/               # Module-level test fixtures
├── scratchhhh/              # Developer scratch (gitignored)
└── .planning/               # GSD planning docs (this directory)
```

## Directory Purposes

**`workflows/`:**
- Purpose: Top-level pipeline workflows
- Contains: `deepcsa.nf` — the single `DEEPCSA` workflow with all channel wiring
- Key files: `workflows/deepcsa.nf`

**`subworkflows/local/`:**
- Purpose: Reusable analysis stages composed of modules
- Contains: One `main.nf` per stage with `include` + `workflow X { take/main/emit }`
- Key files: `subworkflows/local/omega/main.nf`, `subworkflows/local/depthanalysis/main.nf`, `subworkflows/local/mutationpreprocessing/main.nf`

**`subworkflows/nf-core/`:**
- Purpose: Vendored nf-core subworkflows, pinned via `modules.json`
- Contains: `utils_nextflow_pipeline`, `utils_nfcore_pipeline`, `utils_nfvalidation_plugin`, `vcf_annotate_ensemblvep`, `vcf_annotate_ensemblvep_panel`
- Key files: `subworkflows/nf-core/utils_nfcore_pipeline/main.nf`

**`modules/local/`:**
- Purpose: Local process definitions (one `process` per `main.nf`)
- Contains: 97 `main.nf` files; `bbgtools/` sub-tree wraps bbglab tools (omega, oncodrive*, regressions)
- Key files: `modules/local/bbgtools/omega/estimator/main.nf`, `modules/local/createpanels/consensus/main.nf`, `modules/local/table2groups/main.nf`

**`modules/nf-core/`:**
- Purpose: Vendored nf-core modules
- Contains: `multiqc`, `custom/dumpsoftwareversions`
- Key files: `modules/nf-core/multiqc/main.nf`

**`bin/`:**
- Purpose: Executable analysis scripts; Nextflow adds this dir to task PATH automatically
- Contains: Python CLIs (click, pandas, pybedtools, polars), R scripts, shared `utils*.py` helpers, `test/` pytest suite, `saturation_mutagenesis/` extras
- Key files: `bin/utils.py`, `bin/read_utils.py`, `bin/mut_density_simple.py`, `bin/create_consensus_panel.py`, `bin/dNdS_run.R`

**`conf/`:**
- Purpose: Nextflow configuration includes
- Contains: `base.config` (resources/labels/retry), `modules.config` (per-process `ext.args` + publishDir), `results_outputs.config` (final outputs), `tools/*.config`, `modes/*.config`, environment configs (`test.config`, `exome.config`, `mice.config`, `local.config`)
- Key files: `conf/base.config`, `conf/modules.config`

**`assets/`:**
- Purpose: Static data, schemas, templates shipped with the pipeline
- Contains: samplesheet schema, MultiQC config, notification templates, reference datasets (trinucleotide counts, chromosome bands, omega covariates, regressions configs), example inputs
- Key files: `assets/schema_input.json`, `assets/omega_consequences_groupings.json`, `assets/placeholder_no_file.tsv` (used as empty-channel placeholder in `workflows/deepcsa.nf`)

**`tests/` + `test_data/`:**
- Purpose: nf-test pipeline tests and fixtures
- Contains: `tests/deepcsa.nf.test` (+ snapshot), `tests/nextflow.config`, module fixtures in `test_data/modules/`
- Key files: `tests/deepcsa.nf.test`, `nf-test.config`

**`docs/`:**
- Purpose: User/developer documentation
- Contains: `usage.md`, `output.md`, `metrics.md`, `tools.md`, `input_scenarios.md`, `file_formatting.md`, `issue_resolution.md`, `test_data.md`, `images/`

## Key File Locations

**Entry Points:**
- `main.nf`: Pipeline entry; runs `PIPELINE_INITIALISATION` then `BBGTOOLS` → `DEEPCSA`
- `workflows/deepcsa.nf`: All pipeline logic; the file to read to understand data flow

**Configuration:**
- `nextflow.config`: Params defaults, container registries, profiles (docker/singularity/conda/test/debug)
- `nextflow_schema.json`: Param definitions/validation (nf-schema)
- `conf/base.config`: Resource labels (`cpu_low/medium/high`, `mem_low`), error strategy
- `conf/modules.config`: Per-process `ext.args`, `ext.prefix`, publishDir overrides
- `conf/results_outputs.config`: Curated user-facing output paths
- `conf/modes/*.config`: Preset param bundles (`basic`, `clonal_structure`, `get_signatures`)
- `conf/tools/*.config`: Per-tool param bundles
- `assets/schema_input.json`: Samplesheet schema

**Core Logic:**
- `bin/*.py`: All computational logic (mutation density, profiles, panels, omega QC, plotting)
- `bin/utils.py`, `bin/utils_filter.py`, `bin/utils_impacts.py`, `bin/utils_context.py`, `bin/utils_plot.py`, `bin/read_utils.py`: Shared Python helpers — import these, don't duplicate
- `bin/*.R`: dNdScv and mutation-rate R analyses

**Testing:**
- `tests/deepcsa.nf.test`: End-to-end pipeline test (nf-test + snapshots)
- `bin/test/test_*.py`: Python unit tests (pytest)
- `nf-test.config`: Test runner config (profile `test,singularity`, ignores `modules/nf-core/**`)

## Naming Conventions

**Files:**
- Nextflow modules/subworkflows: always `main.nf` inside a lowercase snake_case tool directory: `modules/local/bbgtools/omega/estimator/main.nf`
- Python scripts: `snake_case.py` matching the tool/analysis name: `create_consensus_panel.py`, `mut_density_adjusted.py`
- R scripts: `snake_case.R` with tool prefix: `dNdS_run.R`, `signatures_msighdp_run.R`
- Configs: lowercase with `.config`: `base.config`, `results_outputs.config`

**Directories:**
- Module dirs: lowercase, no separators or underscores preferred: `createpanels`, `filterdepths`, `mutatedcells` (some legacy underscores: `omega_covariates`, `select_mutdensity`)
- bbglab tool wrappers grouped under `modules/local/bbgtools/<tool>/<subcommand>/`

**Code symbols:**
- Nextflow processes: `UPPER_SNAKE_CASE` verbs/phrases: `CREATECONSENSUSPANELS`, `OMEGA_ESTIMATOR`, `COMPUTEDEPTHS`
- Workflows: `UPPER_SNAKE_CASE`: `DEEPCSA`, `OMEGA_ANALYSIS`, `DEPTH_ANALYSIS`
- Aliases when reusing: descriptive suffixes: `MUTDENSITYALL`, `MUTDENSITYPROT`, `OMEGAMULTI`, `PREPROCESSINGGLOBALLOC`

## Where to Add New Code

**New analysis stage (multi-process):**
- Subworkflow: `subworkflows/local/<stage_name>/main.nf` — follow `subworkflows/local/depthanalysis/main.nf` structure (`take:`/`main:`/`emit:`)
- Wire it in `workflows/deepcsa.nf` with an alias include and a `params.<flag>`-gated call

**New single process:**
- Module: `modules/local/<category>/<tool>/main.nf` (use `bbgtools/<tool>/` for bbglab tools)
- Declare `conda "..."` + `container 'docker://...'`, `tag "$meta.id"`, meta-tuple inputs, named `emit:` outputs, `versions.yml` with `topic: versions`, and a `stub:` block
- Add per-process `ext.args`/publishDir in `conf/modules.config`; final outputs in `conf/results_outputs.config`

**New analysis script:**
- Python: `bin/<snake_case_name>.py` — use click for CLI, pandas for tables, import shared helpers from `bin/utils.py` / `bin/read_utils.py` (same dir, so plain `import` works)
- R: `bin/<snake_case_name>.R`
- The script becomes callable by name in any process script block (Nextflow puts `bin/` on PATH)

**New params:**
- Add default in `nextflow.config` params block
- Add schema entry in `nextflow_schema.json`
- If tool-specific, also set in `conf/tools/<tool>.config`

**New tests:**
- Pipeline-level: `tests/deepcsa.nf.test` (nf-test); module fixtures under `test_data/modules/`
- Python unit tests: `bin/test/test_<module>.py` (pytest)

**New static reference data:**
- `assets/<topic>/` (e.g. `assets/trinucleotide_counts/`, `assets/omega-covariates/`)

**New utilities:**
- Shared helpers: `bin/utils*.py` (extend existing files rather than creating new ones when the topic matches)

## Special Directories

**`bin/`:**
- Purpose: Executables auto-mounted on task PATH by Nextflow
- Generated: No
- Committed: Yes (except `bin/__pycache__/`, `bin/saturation_mutagenesis/` notebooks partially gitignored)

**`modules/nf-core/`, `subworkflows/nf-core/`:**
- Purpose: Vendored nf-core components pinned by `modules.json` (nf-core tools install/update)
- Generated: Managed by `nf-core modules install` — do not hand-edit
- Committed: Yes

**`assets/`:**
- Purpose: Static pipeline data/templates
- Generated: No (except `assets/HDP_files*` which is gitignored)
- Committed: Yes (with gitignored exceptions: `assets/useful_scripts/*.ipynb`, `assets/HDP_files*`)

**`scratchhhh/`:**
- Purpose: Developer scratch outputs (SigProfiler runs, plots, test notes)
- Generated: Yes (ad hoc)
- Committed: No (gitignored)

**`work/`, `results/`, `.nf-test/`, `testing/`:**
- Purpose: Nextflow/nf-test runtime outputs
- Generated: Yes
- Committed: No (gitignored)

**`.planning/`:**
- Purpose: GSD planning documents (this analysis)
- Generated: By GSD commands
- Committed: Per team convention

---

*Structure analysis: 2026-10-09*
