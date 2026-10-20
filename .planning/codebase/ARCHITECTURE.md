---
last_mapped_commit: 3c83c5b8aa313077e0ce43239a7b818b281b32ea
last_mapped_at: 2026-10-09
---
<!-- refreshed: 2026-10-09 -->

# Architecture

**Analysis Date:** 2026-10-09

## System Overview

deepCSA (`bbglab/deepCSA`) is a **Nextflow DSL2 pipeline** for analysis of the clonal structure of tissues using duplex-sequencing data. It is a workflow-orchestration architecture: Groovy/Nextflow processes orchestrate containerized Python/R scripts.

```text
┌─────────────────────────────────────────────────────────────────────┐
│                       ENTRY POINT / ORCHESTRATION                    │
│   `main.nf`  →  workflow BBGTOOLS  →  DEEPCSA                       │
│   `workflows/deepcsa.nf` (798 lines — the pipeline "spine")         │
├──────────────────┬──────────────────┬───────────────────────────────┤
│  SUBWORKFLOWS    │  CONFIG LAYER    │  PARAM VALIDATION             │
│  `subworkflows/  │  `nextflow.config`│  nf-schema plugin            │
│   local/*` (23)  │  `conf/*.config` │  `nextflow_schema.json`       │
└────────┬─────────┴──────────────────┴───────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│                        PROCESS MODULES                               │
│  `modules/local/**/main.nf` (97 process definitions, 48 top dirs)   │
│  `modules/nf-core/**` (multiqc, custom/dumpsoftwareversions)        │
│  `subworkflows/nf-core/**` (utils, vep annotation)                  │
└────────┬────────────────────────────────────────────────────────────┘
         │  invokes via PATH (bin/ auto-mounted by Nextflow)
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     ANALYSIS SCRIPTS (bin/)                          │
│  ~90 Python/R scripts + shared utils (`bin/utils*.py`)              │
│  Run inside containers: docker.io/bbglab/*, quay.io registry        │
└────────┬────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  OUTPUT: ${params.outdir}/  (mutdensity/, mutational_profile/,       │
│  plots/, depths/, pipeline_info/ — publish paths in conf/*.config)  │
└─────────────────────────────────────────────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| `main.nf` | Entry point; declares `BBGTOOLS` wrapper workflow and runs `PIPELINE_INITIALISATION` | `main.nf` |
| `DEEPCSA` workflow | The entire pipeline logic: channel wiring, feature flags, fan-out of analyses | `workflows/deepcsa.nf` |
| `PIPELINE_INITIALISATION` | Banner, version validation, params summary log | `subworkflows/local/utils_nfcore_deepcsa/main.nf` |
| Local subworkflows | Compose modules into logical analysis stages (depths, panels, omega, signatures...) | `subworkflows/local/*/main.nf` |
| Local modules | Single `process` definitions with container/conda, script invocation, versions.yml | `modules/local/**/main.nf` |
| nf-core modules/subworkflows | Vendored community components (MultiQC, VEP annotation, utils) | `modules/nf-core/`, `subworkflows/nf-core/` |
| `bin/` scripts | Actual analysis logic (Python with click/pandas, some R) | `bin/*.py`, `bin/*.R` |
| Config layer | Params defaults, resources, publish paths, tool/mode presets | `nextflow.config`, `conf/*.config` |

## Pattern Overview

**Overall:** nf-core-style DSL2 pipeline with local (non-nf-core) module tree.

**Key Characteristics:**
- **Single mega-workflow:** All pipeline logic lives in `workflow DEEPCSA` inside `workflows/deepcsa.nf`. `main.nf` is a thin wrapper (`BBGTOOLS` just calls `DEEPCSA`).
- **Feature-flag fan-out:** Analyses are toggled by boolean params (`params.omega`, `params.oncodrivefml`, `params.dnds`, `params.signatures`, `params.mutationdensity`, `params.profileall`, ...) checked with `if` blocks in the workflow body.
- **Subworkflow reuse via aliasing:** The same subworkflow is included multiple times with aliases, e.g. `MUTATION_DENSITY as MUTDENSITYALL/PROT/NONPROT/SYNONYMOUS` (`workflows/deepcsa.nf:33-37`).
- **Meta-map tuples:** Data flows as `tuple val(meta), path(file)` where `meta.id` is the sample/group key.
- **Channel aggregation with `collectFile`:** Per-sample outputs are flattened and concatenated into cohort-level files (`all_mutdensities.tsv`, `all_profile_stabilities.tsv`) directly in the workflow.
- **Result accumulation channel:** `positive_selection_results` is progressively `.join(..., remainder: true)`-ed with each positive-selection tool's output, then filtered and fed to `PLOTTINGSUMMARY`.

## Layers

**Orchestration layer:**
- Purpose: Entry point and pipeline initialization
- Location: `main.nf`
- Contains: `BBGTOOLS` wrapper workflow, `PIPELINE_INITIALISATION` invocation
- Depends on: `workflows/deepcsa.nf`, `subworkflows/local/utils_nfcore_deepcsa`
- Used by: `nextflow run` CLI

**Workflow layer (the spine):**
- Purpose: All channel wiring, branching, and analysis fan-out
- Location: `workflows/deepcsa.nf`
- Contains: `workflow DEEPCSA` (~700 lines of channel logic), imports of all subworkflows/modules
- Depends on: every subworkflow in `subworkflows/local/`, many modules in `modules/local/` and `modules/nf-core/`
- Used by: `main.nf`

**Subworkflow layer:**
- Purpose: Reusable multi-process analysis stages
- Location: `subworkflows/local/<stage>/main.nf` (23 workflows)
- Contains: `include` of modules + `workflow X { take: main: emit: }` blocks
- Key stages: `depthanalysis`, `createpanels`, `mutationpreprocessing`, `mutationdensity`, `mutationprofile`, `mutability`, `omega`, `oncodrivefml`, `oncodrive3d`, `oncodriveclustl`, `dnds`, `indels`, `signatures`, `signatures_hdp`, `mutatedcells`, `regressions`, `plotdepths`, `plotting_qc`, `plottingsummary`, `enrichpanels`, `adjmutdensity`, `input_check.nf`
- Depends on: `modules/local/**`, `modules/nf-core/**`
- Used by: `workflows/deepcsa.nf`

**Module (process) layer:**
- Purpose: One Nextflow `process` per file, wrapping a script/tool invocation
- Location: `modules/local/<category>/<tool>/main.nf` (97 `main.nf` files)
- Contains: process definition with `conda`, `container`, `input/output` (meta tuples, named `emit:`), `script:` heredoc, `stub:` block, `versions.yml` emission
- Depends on: `bin/` scripts (resolved via Nextflow's automatic `bin/` PATH) and container images
- Used by: subworkflows and directly by `workflows/deepcsa.nf`

**Script layer:**
- Purpose: Actual computation logic
- Location: `bin/*.py`, `bin/*.R`, `bin/saturation_mutagenesis/`
- Contains: Python 3 CLI scripts (click + pandas + pybedtools/polars), R scripts (dNdScv, mutrate), shared helpers (`bin/utils.py`, `bin/utils_context.py`, `bin/utils_filter.py`, `bin/utils_impacts.py`, `bin/utils_plot.py`, `bin/read_utils.py`)
- Depends on: conda envs declared per-process (`python=3.10`, `pybedtools`, `polars`, `click`)
- Used by: module processes

**Config layer:**
- Purpose: Defaults, resources, publishing, tool presets
- Location: `nextflow.config` (params defaults + profiles), `conf/base.config` (resource limits, labels, retry strategy), `conf/modules.config` (78 `withName` overrides: `ext.args`, `publishDir`), `conf/results_outputs.config` (final publish paths), `conf/tools/*.config` (per-tool: omega, oncodrive3d, oncodrivefml, mutdensity, panels, regressions, hdp), `conf/modes/*.config` (presets: `basic`, `clonal_structure`, `get_signatures`)
- Depends on: `nextflow_schema.json` (nf-schema validation), `assets/schema_input.json` (samplesheet schema)
- Used by: Nextflow launcher

## Data Flow

### Primary Request Path

1. **Input validation** — samplesheet CSV validated against `assets/schema_input.json` via nf-schema; `INPUT_CHECK` → `SAMPLESHEET_CHECK` builds `sample_inputs_ch` of `[meta, vcf, bam]` (`subworkflows/local/input_check.nf`). Alternative entry: `--input_maf` + `--use_custom_depths` converts MAF→VCF via `INPUTMAF2VCF` (`workflows/deepcsa.nf:224-238`).
2. **Grouping definition** — `TABLE2GROUP` parses the features table into JSON group definitions (`modules/local/table2groups/main.nf`); group keys are extracted with `JsonSlurper` in the workflow (`workflows/deepcsa.nf:247-260`).
3. **Depth computation** — `DEPTHANALYSIS` computes per-position depths from BAMs (`COMPUTEDEPTHS`) or accepts a custom depths table, then filters by minimum depth (`subworkflows/local/depthanalysis/main.nf`).
4. **Panel creation** — `CREATEPANELS` builds consensus panels (all/exons/protein-coding/non-coding/synonymous BEDs + VEP-annotated panels) from depths + WGS trinucleotide counts (`subworkflows/local/createpanels/main.nf`, `modules/local/createpanels/{captured,consensus,compare,custombedfile}`).
5. **Mutation preprocessing** — `MUT_PREPROCESSING` annotates VCFs with Ensembl VEP (`subworkflows/nf-core/vcf_annotate_ensemblvep*`), blacklists artifacts, builds a mask matrix; emits `somatic_mafs` (`subworkflows/local/mutationpreprocessing/main.nf`).
6. **Depth annotation & enrichment** — `ANNOTATEDEPTHS` merges depths with panels; optional `DOWNSAMPLEDEPTHS`; `ENRICHPANELS` expands panels with subgenic regions and DNA→protein/domain mappings (`subworkflows/local/enrichpanels/main.nf`).
7. **Parallel analyses (feature-flag gated)** — each consumes `somatic_mutations` + relevant consensus BED/panel + depths:
   - Mutational profiles: `MUTPROFILEALL/NONPROT/EXONS/INTRONS`
   - Mutation density: `MUTDENSITYALL/PROT/NONPROT/SYNONYMOUS` + adjusted variant `MUTDENSITYADJUSTED` → `DNDSPROXY`
   - Mutability: `MUTABILITYALL/NONPROT` (feeds oncodrivefml/3d/clustl)
   - Positive selection: `ONCODRIVEFMLALL`, `ONCODRIVE3D`, `ONCODRIVECLUSTL`, `DNDS`, `OMEGA` (+ `OMEGAMULTI`, `OMEGANONPROT` variants), `INDELSSELECTION`
   - Clonal structure: `EXPECTEDMUTATEDCELLS`, `MUTATEDCELLSVAF`, `VAFSMOOTHING`
   - Signatures: `MAF2VCF` → `SIGPROMATRIXGENERATOR` → `SIGNATURESALL/...` → `MUTS2SIGS`
8. **QC & summary** — `PLOTTINGQC` flags failing omegas/QC metrics; results are combined into `positive_selection_results_ready` and passed to `PLOTTINGSUMMARY` (`workflows/deepcsa.nf:640-712`).
9. **Regressions (optional)** — `REGRESSIONSMUTDENSITY/OMEGA/OMEGAGLOB` correlate metrics with metadata (`subworkflows/local/regressions/main.nf`).
10. **Publishing** — `publishDir` paths defined per-process in `conf/modules.config` and `conf/results_outputs.config`; `CUSTOM_DUMPSOFTWAREVERSIONS` collates `versions.yml` topics; MultiQC report at the end.

### Aggregation Pattern (secondary flow)

1. Per-sample process emits `tuple val(meta), path(result.tsv)`.
2. Workflow maps to file only: `.out.mutdensities.map{ it -> it[1] }.flatten()`.
3. `.collectFile(name: "all_*.tsv", storeDir: "${params.outdir}/...", skip: 1, keepHeader: true)` writes the cohort file directly to the outdir (`workflows/deepcsa.nf:316-322, 395-410`).

**State Management:**
- No persistent state; everything is Nextflow channels + task work dirs.
- Placeholder channels (`channel.value(file("${projectDir}/assets/placeholder_no_file.tsv"))`) initialize optional result channels so downstream `.join(..., remainder: true)` never blocks (`workflows/deepcsa.nf:141-148`).
- `scratchhhh/` is a gitignored developer scratch area (not part of the pipeline).

## Key Abstractions

**Meta map (`meta`):**
- Purpose: Sample/group identity carried with every file tuple
- Examples: `workflows/deepcsa.nf` (`meta.id`, `meta.batch`), `subworkflows/local/input_check.nf` (`create_input_channel`)
- Pattern: `tuple val(meta), path(file)`; `meta.id` used for file prefixes via `task.ext.prefix`

**Process module:**
- Purpose: One tool invocation, fully self-contained
- Examples: `modules/local/createpanels/consensus/main.nf`, `modules/local/bbgtools/omega/estimator/main.nf`
- Pattern: `tag "$meta.id"`, `conda "..."`, `container 'docker://bbglab/...'`, `label 'cpu_medium'`, named `emit:` outputs, `versions.yml` with `topic: versions`, `stub:` block for testing

**Feature-flag params:**
- Purpose: Enable/disable analysis branches
- Examples: `nextflow.config` params block (`omega`, `omega_multi`, `omega_globalloc`, `oncodrivefml`, `oncodrive3d`, `dnds`, `signatures`, `mutationdensity`, `profileall`, `regressions`, `downsample`, ...)
- Pattern: derived booleans in the workflow (`run_mutabilities`, `run_mutdensity`, `run_profile_all` at `workflows/deepcsa.nf:206-208`)

**Grouping JSONs:**
- Purpose: Define sample/gene groupings used by omega, regressions, plotting
- Examples: `modules/local/table2groups/main.nf`, `assets/omega_consequences_groupings.json`
- Pattern: `TABLE2GROUP` emits `json_samples`/`json_groups`/`json_allgroups`; keys parsed in workflow to build `samples_keys_ch`/`group_keys_ch`

## Entry Points

**`nextflow run main.nf` (or repo root):**
- Location: `main.nf`
- Triggers: CLI with `--input` samplesheet (or `--input_maf`), profile (`docker`/`singularity`/`conda`), optional `-c conf/modes/<mode>.config`
- Responsibilities: init (banner, version check, params validation via `PIPELINE_INITIALISATION`), then run `DEEPCSA`

**Tests:**
- Location: `tests/deepcsa.nf.test` (nf-test with snapshot `tests/deepcsa.nf.test.snap`), `bin/test/test_*.py` (pytest)
- Triggers: `nf-test test` (config in `nf-test.config`, profile `test,singularity`)

## Architectural Constraints

- **Threading:** Nextflow manages parallelism; per-process resources via labels `cpu_low` (2 cpus), `cpu_medium` (4), `cpu_high` (8), `mem_low` (1 GB) in `conf/base.config`; global caps `params.max_cpus/max_memory/max_time`; per-process `withName` overrides for heavy steps (e.g. `CREATEPANELS:SITESFROMPOSITIONS` 8 GB).
- **Error handling:** `errorStrategy` retries with exponential backoff for exit codes 130–145 and 104 (`conf/base.config`); labels `error_ignore` and `error_retry` opt processes out/in.
- **Containers:** Default registry `quay.io` (`docker.registry`/`singularity.registry` in `nextflow.config`); bbglab tools pinned images (e.g. `docker.io/bbglab/omega:0.2.1`, `bbglab/deepcsa_bed:latest`). Processes declare both `conda` and `container`.
- **bin/ on PATH:** Nextflow auto-prepends `bin/` to the task PATH — scripts in `bin/` are callable by name from any process script block. Do not duplicate scripts elsewhere.
- **Global state:** None beyond `params`; workflow-level Groovy variables (`positive_selection_results`, `all_compiled_omegas`, ...) are channel builders inside `DEEPCSA`.
- **Vendored nf-core code:** Only 2 nf-core modules (`multiqc`, `custom/dumpsoftwareversions`) and 5 nf-core subworkflows are vendored and pinned in `modules.json`; nf-core dirs are excluded from nf-test (`ignore 'modules/nf-core/**/*'` in `nf-test.config`).
- **Circular imports:** None observed; dependency direction is strictly `main.nf → workflows → subworkflows → modules → bin/`.

## Anti-Patterns

### Monolithic workflow file

**What happens:** `workflows/deepcsa.nf` is ~800 lines containing all channel wiring, JSON parsing, aggregation, and branching.
**Why it's wrong:** Hard to test in isolation; any change risks breaking unrelated branches; merge conflicts likely.
**Do this instead:** Follow the existing subworkflow pattern — group related processes into `subworkflows/local/<stage>/main.nf` with `take/main/emit` (see `subworkflows/local/depthanalysis/main.nf` as the model) and keep `workflows/deepcsa.nf` as wiring only.

### Inline Groovy logic in the workflow

**What happens:** `JsonSlurper` parsing of group JSONs and result-channel reshaping happen inline in `workflows/deepcsa.nf:247-260, 713-720`.
**Why it's wrong:** Untestable outside Nextflow; mixes orchestration with data transformation.
**Do this instead:** Push parsing into a `bin/` script invoked by a module process, or a small subworkflow.

## Error Handling

**Strategy:** Process-level retry with exponential backoff; validation fails fast.

**Patterns:**
- `errorStrategy` retry on transient exit codes (130–145, 104) with `sleep(Math.pow(2, task.attempt) * 200)` (`conf/base.config`)
- Labels `error_ignore` / `error_retry` for per-process overrides (`conf/base.config`)
- Explicit `error "..."` guards in the workflow for invalid param combinations (e.g. `--input_maf` without `--use_custom_depths`, `--omega_covariates` without `--omega`) (`workflows/deepcsa.nf:210-218`)
- nf-schema validation of params and samplesheet (`nextflow_schema.json`, `assets/schema_input.json`)
- Python scripts raise/exit non-zero on bad input; shared helpers in `bin/utils.py`, `bin/utils_filter.py`

## Cross-Cutting Concerns

**Logging:** Nextflow task logs; `PIPELINE_INITIALISATION` prints banner + params summary (`paramsSummaryMap` from nf-schema plugin); MultiQC aggregates QC (`assets/multiqc_config.yml`).
**Validation:** nf-schema for params (`params.validate_params`), JSON schema for samplesheet (`assets/schema_input.json`), `SAMPLESHEET_CHECK` module.
**Authentication:** None (HPC/cluster execution; no external API auth). Seqera Platform integration via `tower.yml` (report display only).
**Versioning:** Every process emits `versions.yml` (topic `versions`), collated by `CUSTOM_DUMPSOFTWAREVERSIONS`; pipeline version in `params.version` checked at init.
**Publishing:** Two-tier publish config — default per-process paths in `conf/modules.config`, curated final outputs in `conf/results_outputs.config`; `versions.yml` never published (`saveAs` filter).

---

*Architecture analysis: 2026-10-09*
