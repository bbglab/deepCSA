---
last_mapped_commit: 3c83c5b8aa313077e0ce43239a7b818b281b32ea
last_mapped_at: 2026-10-09
---
# External Integrations

**Analysis Date:** 2026-10-09

## APIs & External Services

**None at runtime.** The pipeline is fully offline: VEP runs with `--cache --offline` (`params.vep_params` in `nextflow.config`), and `BGDATA_OFFLINE=TRUE` is forced in the `env` block. No HTTP calls, no cloud SDKs, no database drivers anywhere in `bin/` or module scripts.

**External data resources (consumed as files, not APIs):**

| Resource | Purpose | Default param (`nextflow.config`) | Site override |
| --- | --- | --- | --- |
| Ensembl VEP cache (v111, GRCh38, homo_sapiens) | Variant annotation | `vep_cache = ".vep"` (staged into process) | `conf/general_files_IRB.config`: `/data/bbg/datasets/vep` |
| COSMIC v3.4 SBS signatures (GRCh38) | Signature fitting | `cosmic_ref_signatures` | IRB: COSMIC v3.5 at `/data/bbg/datasets/COSMIC_signatures/` |
| COSMIC v3.4 ID signatures (GRCh37) | Indel signature fitting | `indel_ref_signatures` | IRB: v3.5 |
| CADD v1.7 whole-genome SNV scores (+ .tbi) | OncodriveFML | `cadd_scores`, `cadd_scores_ind` | IRB: `/data/bbg/datasets/CADD/v1.7/hg38/` |
| Genome FASTA (GRCh38 no-alt masked) | Reference sequence | `fasta` (user-supplied) | IRB: `/data/bbg/datasets/genomes/GRCh38/...` |
| dNdScv Biomart reference (Ensembl v111) | dN/dS CDS mapping | `dnds_biomart_ref` | IRB: MANE version |
| dNdScv covariates (epigenome/PCAWG, .rda) | dN/dS covariates | `dnds_covariates` | IRB path |
| Omega covariates TSV (hg19/hg38 epigenome PCAWG) | Covariate-aware Omega | `omega_covariates_cov_file` (default ships in `assets/omega-covariates/`) | IRB path |
| Oncodrive3D datasets + annotations | 3D selection | `datasets3d`, `annotations3d` | IRB: dated snapshots (`datasets-260603`) |
| NanoSeq SNP/noise masks (BED) | Artifact filtering | `nanoseq_snp`, `nanoseq_noise` | IRB paths |
| GFF3 (Homo_sapiens.GRCh38.111) | DNA→protein/domain mapping | `gff3_file` | IRB path |
| Whole-genome trinucleotide counts | Mutability normalization | `wgs_trinuc_counts` (ships in `assets/trinucleotide_counts/`) | IRB path |
| bbgdomains annotated TSV | Subgenic/domain regions | `domains_file` (IRB only) | IRB path |
| gnomAD allele frequencies | Germline filtering | Embedded in VEP annotation (`--af_gnomadg --af_gnomade`), threshold `gnomad_af_threshold = 1e-3` | — |

Consumption pattern: `channel.fromPath(params.X, checkIfExists: true)` in `workflows/deepcsa.nf` and `subworkflows/local/*/main.nf` — all file-based, no network fetch.

## Data Storage

**Databases:**
- None. All I/O is TSV/MAF/VCF/BED/PDF files on a shared filesystem.

**File Storage:**
- Local/shared filesystem only. Outputs published under `params.outdir` with routing defined in `conf/modules.config` (`${params.outdir}/processing_files/<process-prefix>/`) and `conf/results_outputs.config` (curated outputs: `depths/`, `plots/`, `qc/`, `selection/`, `mutdensity/`, `signatures/`, …). `publish_dir_mode = 'copy'` by default.
- Some cohort-level aggregates written directly via `collectFile(storeDir: "${params.outdir}/...")` (e.g. `workflows/deepcsa.nf:341`, `subworkflows/local/dnds/main.nf:39-45`).

**Caching:**
- Nextflow work dir (default `work/`); nf-test work dir `.nf-test` or `$DEEPCSA_TEST_WORKDIR` (`nf-test.config`).
- Singularity image cache: `singularity.cacheDir`/`libraryDir` = `/data/bbg/datasets/pipelines/nextflow_containers` (`conf/general_files_IRB.config`, `tests/nextflow.config`).
- No application-level cache service (no Redis/Memcached).

## Authentication & Identity

**Auth Provider:** None. No authentication anywhere — pipeline runs on trusted HPC infrastructure. Notification endpoints are the only externally-addressable values (see below).

## Monitoring & Observability

**Error Tracking:**
- None (no Sentry etc.). Failure handling is Nextflow-native: retry with exponential backoff on exit codes 104/130–145, `maxRetries 3` (`conf/base.config`); labels `error_ignore`/`error_retry`.

**Logs / Reports:**
- Nextflow built-ins enabled in `nextflow.config`: `timeline`, `report`, `trace`, `dag` → timestamped files under `${params.outdir}/pipeline_info/`.
- **MultiQC** aggregates per-process `versions.yml` (emitted with `topic: versions` by ~76 local modules) plus custom `workflow_summary_mqc.yaml` and `methods_description_mqc.yaml` (`workflows/deepcsa.nf:762-790`, config `assets/multiqc_config.yml`).
- `CUSTOM_DUMPSOFTWAREVERSIONS` (`modules/nf-core/custom/dumpsoftwareversions`) captures the full software environment.

**Notifications (outgoing):**
- **Email** — completion/failure email via local `sendmail`/`mail` binaries, HTML rendered from `assets/sendmail_template.txt` / `assets/email_template.html` (`subworkflows/nf-core/utils_nfcore_pipeline/main.nf:296-318`). Triggered by `params.email` / `params.email_on_fail`.
- **Slack / MS Teams** — `imNotification` posts to `params.hook_url`; Slack format (`assets/slackreport.json`) when URL contains `hooks.slack.com`, otherwise Adaptive Cards (`assets/adaptivecard.json`) for Teams (`subworkflows/nf-core/utils_nfcore_pipeline/main.nf:403-404`).

## CI/CD & Deployment

**Hosting:**
- Source: GitHub (`https://github.com/bbglab/deepCSA`, manifest in `nextflow.config`). DOI: `dx.doi.org/10.17504/protocols.io.dm6gp1jodgzp/v2`.
- No `.github/` directory — no GitHub Actions workflows detected in this checkout.

**CI Pipeline:**
- **nf-test** is the test harness (`nf-test.config`, `tests/deepcsa.nf.test` + snapshot). Per `docs/test_data.md` and `tests/nextflow.config`, tests execute on the IRB SLURM cluster (queues `bbg_cpu_zen4,irb_cpu_zen4`) with Singularity; results CSVs (`tests/2026-10-03_results.csv`, `tests/2026-10-09_results.csv`) are committed. No hosted CI runner config found in-repo.
- Container Dockerfile recipes live out-of-repo: `https://github.com/bbglab/containers-recipes` (`docs/tools.md`).

**Deployment model:**
- `nextflow run bbglab/deepCSA -profile <engine>,<site>` — no packaging, no container registry push automation in-repo.

## Environment Configuration

**Required env vars:**
- None mandatory. Optional: `DEEPCSA_TEST_WORKDIR` (nf-test work dir), `HOME=/tmp` and the isolation vars (`PYTHONNOUSERSITE`, `R_PROFILE_USER`, `R_ENVIRON_USER`, `JULIA_DEPOT_PATH`, `BGDATA_OFFLINE`) are set by the pipeline itself (`nextflow.config` `env` block).

**Required params (schema-enforced, `nextflow_schema.json`):**
- `input` (samplesheet CSV, validated against `assets/schema_input.json`) and `outdir` are the only `required` entries. `fasta` required unless reference paths come from a site profile. `input_maf` alternative input mode requires `use_custom_depths = true` (`workflows/deepcsa.nf:189-201`).

**Secrets location:**
- No secrets in repo. Notification hook URLs are passed at runtime via `--hook_url`; email addresses via `--email`. No `.env` files present.

## Webhooks & Callbacks

**Incoming:**
- None.

**Outgoing:**
- Slack webhook (`hooks.slack.com`) or MS Teams Adaptive Card POST, only when `params.hook_url` is set (`subworkflows/nf-core/utils_nfcore_pipeline/main.nf`).
- Local `sendmail` invocation (not an HTTP callback).

## Container Images (external registries)

All pulled from Docker Hub / biocontainers at task launch; registry override via `docker.registry`/`singularity.registry` (default `quay.io` in `nextflow.config`, but module directives hard-code `docker.io`/`biocontainers`).

| Image | Used by | Defined in |
| --- | --- | --- |
| `docker.io/bbglab/deepcsa-core:0.1.0` | ~60 `deepcsa_core`-labeled Python processes | `conf/modules.config:714` |
| `docker.io/bbglab/deepcsa_bed:latest` | Panel consensus (pybedtools/polars) | `modules/local/createpanels/consensus/main.nf` |
| `docker.io/bbglab/omega:0.2.1` | Omega preprocess/mutabilities/estimator | `modules/local/bbgtools/omega/*/main.nf` |
| `docker.io/ferriolcalvet/omegacovariates:v0.1.0` | Covariate-aware Omega | `modules/local/omega_covariates/run/main.nf` |
| `docker.io/spellegrini87/oncodrive3d:1.0.9-light` / `-chimerax` | Oncodrive3D run/plots | `modules/local/bbgtools/oncodrive3d/*/main.nf` |
| `docker.io/ferriolcalvet/oncodrivefml:latest` | OncodriveFML | `modules/local/bbgtools/oncodrivefml/main.nf` |
| `docker.io/ferriolcalvet/oncodriveclustl:latest` | OncodriveCLUSTL | `modules/local/bbgtools/oncodriveclustl/main.nf` |
| `docker.io/ferriolcalvet/dnds:latest` | dNdScv buildref/run | `modules/local/dnds/*/main.nf` |
| `docker.io/ferriolcalvet/sigprofiler_assignment:1.1.3` | Signature fitting | `modules/local/signatures/sigprofiler/assignment/*/main.nf` |
| `docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5` | SBS96 matrix generation | `modules/local/signatures/sigprofiler/matrixgenerator/main.nf` |
| `docker.io/ferriolcalvet/sigprofilerassignment` | SigProfilerExtractor | `modules/local/signatures/sigprofiler/extractor/main.nf` |
| `docker.io/ferriolcalvet/msighdp:latest`, `hdp_wrapper`, `musical:latest` | HDP/MUSICAL signature extraction | `modules/local/signatures/*/main.nf` |
| `docker.io/rblancomi/bbgregressions:dev` | bbgregressions (label-based) | `conf/modules.config:719` |
| `docker.io/ferriolcalvet/runningr:v1` | R mutation-density scaling | `modules/local/mut_density/wgscaled/main.nf` |
| `docker.io/ferriolcalvet/saturation:v0.1.0` | Saturation kinetics | `modules/local/saturation_kinetics/compute/main.nf` |
| `docker.io/axelrosendahlhuber/expected_mutrate:latest` | Expected mutated cells | `modules/local/mutated_cells_expected/main.nf` |
| `docker.io/ferranmuinos/test_mutated_genomes` | Mutated genomes from VAF | `modules/local/mutated_genomes_from_vaf/main.nf` |
| `biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0` (also 102/108 variants) | VEP annotation | `modules/nf-core/ensemblvep/*/main.nf` |
| `biocontainers/samtools:1.18--h50ea8bc_1` | Depth computation | `modules/local/computedepths/main.nf` |
| `biocontainers/tabix:1.11--hdfd78af_0` | Indexed TSV queries | `modules/nf-core/tabix/*/main.nf` |
| `biocontainers/multiqc:1.20--pyhdfd78af_0` | QC report | `modules/nf-core/multiqc/main.nf` |
| `biocontainers/pybedtools:0.9.1--py38he0f268d_0` | Custom BED handling | `modules/local/createpanels/custombedfile/main.nf` |
| `biocontainers/python:3.8.3` | Samplesheet check | `modules/local/samplesheet_check.nf` |

Conda fallbacks exist for a handful of processes (`-profile conda`/`mamba`): inline specs in `modules/local/createpanels/*/main.nf` (`python=3.10.17`, `pybedtools=0.12.0`, `polars=1.30.0`, `click=8.2.1`, `gcc_linux-64=15.1.0`) and `modules/local/computedepths/environment.yml` (`bioconda::samtools=1.18`).

---

*Integration audit: 2026-10-09*
