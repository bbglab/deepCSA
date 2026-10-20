---
last_mapped_commit: 3c83c5b8aa313077e0ce43239a7b818b281b32ea
last_mapped_at: 2026-10-09
---
# Technology Stack

**Analysis Date:** 2026-10-09

## Languages

**Primary:**
- **Nextflow (DSL2)** — pipeline orchestration. `main.nf` (entry point, `nextflow.enable.dsl = 2`), `workflows/deepcsa.nf` (main `DEEPCSA` workflow, ~800 lines), `subworkflows/local/**` (22 subworkflows), `modules/local/**` (~60 local modules). Manifest requires `nextflowVersion = '!>=25.04.2'` (`nextflow.config`).
- **Python 3** — analysis scripts in `bin/*.py` (~90 scripts). CLI via `click`, data via `pandas`/`polars`, plotting via `matplotlib`/`seaborn`, stats via `scipy`/`statsmodels`. Container images run Python 3.10.x (e.g. conda spec `python=3.10.17` in `modules/local/createpanels/captured/main.nf`); `modules/local/samplesheet_check.nf` pins Python 3.8.3. `pyproject.toml` targets py37–py310 for lint config only.

**Secondary:**
- **R** — 4 scripts: `bin/dNdS_run.R` (dNdScv), `bin/mutgenomes_expected_mutrisk.R`, `bin/mutrate_genome_trinuc_corrected.R`, `bin/signatures_msighdp_run.R` (mSigHdp). Bioconductor-heavy (see R libraries below).
- **Groovy** — implicit in Nextflow; used directly in `subworkflows/nf-core/utils_nfcore_pipeline/main.nf` (email/notification logic, `java.io.File`, Groovy template engine).
- **Bash/Shell** — process `script:` blocks throughout modules; `bin/createcustombed.sh`.
- **YAML** — regression configs (`assets/regressions/configs/*.yml`), parsed by `bin/regressions_editconfig.py` (PyYAML).

## Runtime

**Environment:**
- Nextflow >= 25.04.2 on JVM (Java 17+ implied by Nextflow 25.x)
- Nextflow plugin: `nf-schema@2.3.0` (declared in `plugins` block of `nextflow.config`) — used for params validation (`paramsSummaryMap` in `workflows/deepcsa.nf`) and schema-based input checks
- Execution engines (profiles in `nextflow.config`): Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda, Mamba. Default container registry set to `quay.io` but all module `container` directives point at `docker.io` / `biocontainers` explicitly.

**Package Manager:**
- **No application-level package manager / lockfile.** Dependencies are pinned per-process via `container` directives in module files and inline `conda` directives (e.g. `modules/local/computedepths/main.nf` → `environment.yml` with `bioconda::samtools=1.18`).
- `pyproject.toml` is **lint config only** (Black line-length 120, isort black profile) — not an installable package.
- nf-test plugins: `nft-utils@0.0.3` loaded in `nf-test.config`.

## Frameworks

**Core:**
- **nf-core pipeline template conventions** — boilerplate subworkflows `subworkflows/nf-core/utils_nextflow_pipeline`, `utils_nfcore_pipeline`, `utils_nfvalidation_plugin`; `PIPELINE_INITIALISATION` in `subworkflows/local/utils_nfcore_deepcsa/main.nf` (banner, params summary, completion email/notifications).
- **nf-core modules (pinned via `modules.json`, branch master, git_sha `911696ea0b62df80e900ef244d7867d177971f73`):**
  - `modules/nf-core/custom/dumpsoftwareversions` — software version capture for MultiQC
  - `modules/nf-core/ensemblvep` (`vep`, `veppanel`) — variant annotation
  - `modules/nf-core/multiqc` — aggregated QC report
  - `modules/nf-core/tabix` (`bgziptabix`, `bgziptabixquery`) — indexed TSV querying, reused under many aliases (QUERYDEPTHS, QUERYPANEL, DEPTHS.*CONS)
- **nf-core subworkflows:** `subworkflows/nf-core/vcf_annotate_ensemblvep`, `vcf_annotate_ensemblvep_panel`

**Testing:**
- **nf-test** — `nf-test.config` (testsDir `tests`, workDir `.nf-test` or `$DEEPCSA_TEST_WORKDIR`, profile `test,singularity`, ignores `modules/nf-core/**` and `subworkflows/nf-core/**`). Main test: `tests/deepcsa.nf.test` with snapshot file `tests/deepcsa.nf.test.snap`. Per `docs/test_data.md`, tests must run on a SLURM cluster with Singularity.

**Build/Dev:**
- No build step. Linting: Black + isort via `pyproject.toml` (applies to `bin/check_samplesheet.py` per file comment).
- Process resource/error policy centralized in `conf/base.config` (labels `cpu_low`/`cpu_medium`/`cpu_high`/`mem_low`, `error_ignore`, `error_retry`; exponential-backoff retry on exit codes 104, 130–145).

## Key Dependencies

**Critical (Python, imported across `bin/`):**
- `click` — CLI framework for nearly all `bin/*.py` scripts
- `pandas` — tabular mutation/depth/panel data (`bin/utils.py`, `bin/read_utils.py`, `bin/filter_cohort.py`, …)
- `polars` — high-performance panel processing (`bin/create_consensus_panel.py`, `bin/create_panel_versions.py`; conda spec `conda-forge::polars=1.30.0`)
- `matplotlib`, `seaborn` — all PDF plotting (`bin/utils_plot.py`, `bin/plot_*.py`)
- `scipy` — statistical tests (`bin/compute_hotspots_selection.py`, `bin/signatures_musical.py`)
- `statsmodels` — covariate regression models (`bin/omega_covariates_definitions.py`)
- `PyYAML` — regression config editing (`bin/regressions_editconfig.py`)
- `pybedtools` — BED interval ops (in `createpanels` containers)

**Critical (R):**
- `dndscv` — dN/dS selection model (`bin/dNdS_run.R`)
- `BSgenome.Hsapiens.UCSC.hg38`, `GenomicRanges`, `IRanges`, `Biostrings` — genome sequence/ranges (`bin/mutgenomes_expected_mutrisk.R`)
- `tidyverse`, `data.table`, `ggplot2`, `jsonlite`, `optparse`, `Hmisc`, `R.utils`, `abind`
- `ICAMS`, `mSigHdp` — signature extraction (`bin/signatures_msighdp_run.R`)

**Infrastructure (bioinformatics tools, containerized):**
- Ensembl VEP 111 (default; also supports 102/108 via `params.vep_cache_version` conditional containers in `modules/nf-core/ensemblvep/*/main.nf`)
- samtools 1.18 (`modules/local/computedepths`)
- tabix 1.11 (`modules/nf-core/tabix`)
- MultiQC 1.20 (`modules/nf-core/multiqc`)
- bbglab method containers — see INTEGRATIONS.md

## Configuration

**Environment:**
- All runtime config via Nextflow params — no `.env` files. Schema-validated by `nextflow_schema.json` (nf-schema) and `assets/schema_input.json` (samplesheet columns).
- Hardened env block in `nextflow.config`: `PYTHONNOUSERSITE=1`, `R_PROFILE_USER=/.Rprofile`, `R_ENVIRON_USER=/.Renviron`, `JULIA_DEPOT_PATH=/usr/local/share/julia`, `BGDATA_OFFLINE=TRUE`, `HOME=/tmp` — prevents host Python/R libraries leaking into containers.
- Site-specific reference paths in `conf/general_files_IRB.config` (IRB cluster: `/data/bbg/datasets/...` for VEP cache, FASTA, COSMIC, CADD, dNdScv, Oncodrive3D, NanoSeq masks, GFF3) + `singularity.cacheDir`/`libraryDir`.
- Test env var: `DEEPCSA_TEST_WORKDIR` (overrides nf-test work dir, `nf-test.config`).

**Build/Run config files:**
- `nextflow.config` — params defaults, profiles, plugins, env, timeline/report/trace/dag reporting, manifest (name `bbglab/deepCSA`, version `1.0.1.dev`)
- `conf/base.config` — global resource limits (`max_memory 950.GB`, `max_cpus 196`, `max_time 30.d`), CPU/memory tiers, retry strategy, per-process `withName` overrides
- `conf/modules.config` — `ext.args` per module, `publishDir` routing (`${params.outdir}/processing_files/...`), label→container mapping (`deepcsa_core` → `docker.io/bbglab/deepcsa-core:0.1.0`, `bbgregressions` → `docker.io/rblancomi/bbgregressions:dev`); includes `conf/tools/*.config` and `conf/results_outputs.config`
- `conf/tools/` — tool-specific params: `panels.config`, `omega.config`, `mutdensity.config`, `oncodrive3d.config`, `oncodrivefml.config`, `hdp_sig_extraction.config`, `regressions.config`
- `conf/modes/` — analysis-mode profiles: `basic.config`, `clonal_structure.config`, `get_signatures.config`
- Other profiles: `test`, `test_real`, `test_regressions`, `mice`, `exome`, `irbcluster`, `local`, `debug`, `gitpod`, plus one per container engine
- `nextflow_schema.json` — full parameter schema (input/output, genome, grouping, hotspot, etc.)
- `tower.yml` — Seqera Platform report display config

## Platform Requirements

**Development:**
- Nextflow >= 25.04.2 + JVM
- A container engine (Docker locally; Singularity on cluster) or Conda/Mamba
- For tests: SLURM cluster with Singularity, queues `bbg_cpu_zen4,irb_cpu_zen4` (`tests/nextflow.config`); local execution not supported per `docs/test_data.md`

**Production:**
- HPC via SLURM (IRB cluster profile `irbcluster` / `tests/nextflow.config`); any Nextflow executor works in principle (no executor hard-coded outside test/local configs)
- Shared filesystem for reference data (`/data/bbg/datasets/...`) and Singularity image cache (`/data/bbg/datasets/pipelines/nextflow_containers`)
- Resource ceiling: 196 CPUs, 950 GB RAM, 30 days per process (defaults in `nextflow.config`)

---

*Stack analysis: 2026-10-09*
