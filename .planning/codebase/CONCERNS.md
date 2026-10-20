---
last_mapped_commit: 3c83c5b8aa313077e0ce43239a7b818b281b32ea
last_mapped_at: 2026-10-09
---
# Codebase Concerns

**Analysis Date:** 2026-10-09

## Tech Debt

**Unpinned container images (`latest` / untagged):**
- Issue: Many local modules reference floating or untagged Docker images, so rebuilds can silently pick up different tool versions and break reproducibility.
- Files: `modules/local/*/main.nf`, `subworkflows/local/*/main.nf`, `conf/modules.config`
  - `docker://bbglab/deepcsa_bed:latest` (`modules/local/createpanels/consensus/main.nf`)
  - `docker.io/ferriolcalvet/dnds:latest`, `docker.io/ferriolcalvet/oncodrivefml:latest`, `docker.io/ferriolcalvet/oncodriveclustl:latest`, `docker.io/ferriolcalvet/musical:latest`, `docker.io/ferriolcalvet/msighdp:latest`, `docker.io/axelrosendahlhuber/expected_mutrate:latest`
  - Untagged: `docker.io/ferriolcalvet/hdp_wrapper`, `docker.io/ferranmuinos/test_mutated_genomes`, `docker.io/ferriolcalvet/sigprofilerassignment`
  - Dev-tagged: `docker.io/rblancomi/bbgregressions:dev` (`conf/modules.config`)
- Impact: Non-reproducible results across runs/nodes; a pushed `latest` image can break the pipeline without any repo change.
- Fix approach: Pin every container to an immutable version tag (pattern already used by `docker.io/bbglab/omega:0.2.1`, `docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5`, `docker.io/spellegrini87/oncodrive3d:1.0.9-light`).

**No Python dependency manifest for `bin/` scripts:**
- Issue: `pyproject.toml` only configures Black/isort; there is no `requirements.txt` or dependency list for the ~86 Python scripts in `bin/` (pandas, numpy, polars, click, matplotlib, etc.).
- Files: `pyproject.toml`, `bin/*.py`
- Impact: Local development and unit tests (`bin/test/`) depend on whatever is in the current conda env; version drift between the `deepcsa-core` container and local envs causes subtle failures.
- Fix approach: Add a pinned `requirements.txt` (or `[project.dependencies]` in `pyproject.toml`) matching the `bbglab/deepcsa-core` container contents.

**Pending pandas 2.2.3 upgrade:**
- Issue: Two scripts carry `# TODO: bump pandas to 2.2.3`, indicating known incompatibilities with newer pandas.
- Files: `bin/concat_sbs_probs.py:3`, `bin/mut_density_simple.py:16`
- Impact: Blocks dependency modernization; mixed pandas/polars codebase (`bin/create_consensus_panel.py`, `bin/create_panel_versions.py`, `bin/merge_annotation_depths.py`, `bin/panel_custom_processing.py`, `bin/panels_computedna2protein.py` use polars) increases maintenance surface.
- Fix approach: Audit pandas-2.x breaking changes in these scripts, upgrade, and remove the TODOs.

**Duplicated annotation post-processing logic:**
- Issue: `bin/postprocessing_annotation.py` (263 lines) and `bin/panel_postprocessing_annotation.py` (215 lines) are near-duplicates with the same TODOs repeated in both (`bin/postprocessing_annotation.py:85,190`, `bin/panel_postprocessing_annotation.py:70,164`).
- Impact: Bug fixes must be applied twice; the two files can drift (one already imports `utils.vartype`, the other does not).
- Fix approach: Extract shared consequence-mapping/muttype-conversion logic into `bin/utils_impacts.py` or a new shared module and parameterize the differences.

**Stub CHANGELOG and incomplete nf-core template remnants:**
- Issue: `CHANGELOG.md` still contains the template placeholder `## v1.0dev - [date]` with empty Added/Fixed/Dependencies/Deprecated sections, while `nextflow.config` manifest declares `version = '1.0.1.dev'`. `subworkflows/local/utils_nfcore_deepcsa/main.nf:199-201` still has the placeholder Zenodo DOI TODO.
- Files: `CHANGELOG.md`, `nextflow.config:392`, `subworkflows/local/utils_nfcore_deepcsa/main.nf`
- Impact: Release history is not tracked; citation text may render incomplete.
- Fix approach: Backfill changelog per release; register Zenodo DOI and fill `toolCitationText`/`toolBibliographyText`.

**Scratch/experiment artifacts in repo root:**
- Issue: A 1.1 GB `scratchhhh/` directory (SigProfiler outputs, VCFs) sits in the repo root. It is gitignored (`.gitignore` line `scratchhhh/`) but its name suggests it was renamed to dodge cleanup; stray files `nf-2eI0igpGfSq5RI-reports.tsv`, `.nextflow.log*` (10 rotated logs, ~2 MB total) also accumulate at root. `bin/explore_saturation.ipynb` (an exploration notebook) is git-tracked inside `bin/`.
- Files: `scratchhhh/`, `nf-2eI0igpGfSq5RI-reports.tsv`, `.nextflow.log*`, `bin/explore_saturation.ipynb`
- Impact: Confuses repo structure; risk of accidentally committing large outputs; notebooks in `bin/` mix exploration code with pipeline scripts.
- Fix approach: Move scratch data outside the repo or into a dedicated `scratch/` (already gitignored); delete rotated logs; move exploration notebooks to `assets/useful_scripts/` (where `.ipynb` files are already the convention and gitignored).

**Unexplained/undocumented functions:**
- Issue: Several functions carry explicit "explain what this does" TODOs.
- Files: `bin/mutations_custom_processing.py:17`, `bin/mut_density_adjusted_dnds.py:12-13` (also requests a log file output)
- Impact: Onboarding cost and higher risk of misuse.
- Fix approach: Write docstrings and add the requested statistics log output.

## Known Bugs

**Five failing nf-test process tests (as of 2026-10-09 run):**
- Symptoms: `tests/2026-10-09_results.csv` records 5 FAILED / 20 PASSED:
  - `modules/local/filtermaf/tests/main.nf.test` — "Should run and emit cohort filtered mutations" (FILTER_BATCH)
  - `modules/local/group_genes/tests/main.nf.test` — GROUP_GENES genes-to-groups JSON
  - `modules/local/mut_density/simple/tests/main.nf.test` — MUTATION_DENSITY all-panel densities
  - `modules/local/sig_matrix_concat/tests/main.nf.test` — MATRIX_CONCAT concatenated WGS matrices
  - `modules/local/sitesfrompositions/tests/main.nf.test` — SITESFROMPOSITIONS panel sites chunk
- Files: `modules/local/filtermaf/`, `modules/local/group_genes/`, `modules/local/mut_density/`, `modules/local/sig_matrix_concat/`, `modules/local/sitesfrompositions/`
- Trigger: Run `nf-test` suite (see `nf-test.config`); results CSVs in `tests/` (gitignored via `tests/2*results*`).
- Workaround: None — these indicate current regressions or stale snapshots (`tests/deepcsa.nf.test.snap` may also need regeneration).
- Note: The results CSV references test files (e.g. `modules/local/blacklistmuts/tests/main.nf.test`) that no longer exist in the working tree — only 2 `*.nf.test` files are git-tracked under `modules/local`. The suite was run against a different tree state; reconcile test files with the results.

**Bedtools coordinate hack:**
- Symptoms: `bin/createcustombed.sh:21` contains `#HACK to handle issue arq5x/bedtools2#359` — an awk workaround mutating BED start/end coordinates.
- Files: `bin/createcustombed.sh`
- Trigger: Custom BED file creation path (`use_custom_bedfile = true`).
- Workaround: In place; will need revisiting if bedtools version changes (container `bbglab/deepcsa_bed:latest` is unpinned, compounding risk).

**Dummy channel value in OMEGA subworkflow:**
- Symptoms: `subworkflows/local/omega/main.nf:68` — `// FIXME here I am using bedfile as a dummy value channel` passed into PREPROCESSING.
- Files: `subworkflows/local/omega/main.nf`
- Impact: Fragile channel wiring; refactoring the OMEGA inputs can silently misalign channels.
- Fix approach: Replace with an explicit value channel or remove the unused input.

## Security Considerations

**No secrets handling in repo (good), but no CI either:**
- Risk: No GitHub workflows exist (`.nf-core.yml` template `skip: [github, ci]`), so no automated linting/tests/secret scanning run on pushes.
- Files: `.nf-core.yml`, absent `.github/workflows/`
- Current mitigation: `.pre-commit-config.yaml` exists locally (prettier/markdownlint only).
- Recommendations: Add at minimum a CI workflow running `nf-test` and the `bin/test/` unittest suite; add nf-core `linting.yml` equivalent.

**External subprocess calls:**
- Risk: `bin/merge_annotation_depths.py:161,175,183,191` shells out to `gzip` via `subprocess.run(..., check=True)` (no `shell=True`, arguments are list-form — low risk).
- Files: `bin/merge_annotation_depths.py`
- Current mitigation: List-argument invocation, `check=True`.
- Recommendations: None urgent; keep avoiding `shell=True` with user-derived paths.

**Hardcoded reference/cache versions:**
- Risk: `vep_cache_version = 111` and `dnds_biomart_ref = "homo_sapiens.v111.canonical.biomart.tsv"` are hardcoded defaults (`nextflow.config:138,159`); annotation results silently depend on cache availability and version.
- Files: `nextflow.config`, `conf/`
- Current mitigation: Parameters are overridable via CLI/config.
- Recommendations: Document the required VEP cache version in `docs/usage.md` and validate cache presence at pipeline init.

## Performance Bottlenecks

**Memory-hungry annotation/mapping processes:**
- Problem: `DNA2PROTEINMAPPING` requests `30.GB * task.attempt` and is set to `errorStrategy = 'ignore'` after retries; `PLOTMUTATIONSPECIFIC` and `ONCODRIVEFMLSNVS` request `24.GB * task.attempt`.
- Files: `conf/base.config:211-216,220-222`, `conf/tmp_quick_fixes.config:9-12`
- Cause: Whole-panel tables loaded into pandas in one shot (e.g. `bin/panels_computedna2protein.py` uses polars, but several annotation scripts still use pandas `read_table` on full files).
- Improvement path: Extend the existing chunking mechanism — `params.panel_sites_chunk_size` (default 1,000,000, `nextflow.config:119`, consumed at `conf/modules.config:723` for `SITESFROMPOSITIONS`) — to other site-level processes; prefer polars/lazy evaluation in the largest scripts.

**Silent output loss on OOM:**
- Problem: Because `DNA2PROTEINMAPPING` and several plot processes end in `errorStrategy = 'ignore'`, an OOM failure produces missing downstream outputs rather than a pipeline failure.
- Files: `conf/base.config:201-222`, `conf/modules.config:227` (`PLOTOMEGA`), `conf/tmp_quick_fixes.config`
- Improvement path: Convert to `retry` with a final `fail` (or emit an explicit empty-output marker) so missing results are detectable.

## Fragile Areas

**Positional column renaming in omega QC:**
- Files: `bin/omega_syn_qc.py:187` — `# TODO: dangerous column renaming here` — assigns `syn_muts_df.columns = [...]` positionally after a merge with suffixes `_loc`/`_gloc`.
- Why fragile: Any column-order change in the OMEGA container output (`bbglab/omega:0.2.1`) silently mislabels observed/estimated synonymous counts.
- Safe modification: Rename by explicit column-name mapping and assert expected columns before assignment.
- Test coverage: None (no unit test for this script).

**Column-name robustness in annotation post-processing:**
- Files: `bin/postprocessing_annotation.py:162` (`# TODO: Is it robust enough to use columns names here?`), `bin/postprocessing_annotation.py:146` (bare `# TODO`), `bin/utils_impacts.py:196` (try/except around consequence handling flagged for revision), `bin/utils_context.py:33` (`# TODO remove this try-except`).
- Why fragile: VEP output format changes break parsing; bare `except:` clauses (`bin/postprocessing_annotation.py:178`, `bin/utils_context.py:42`, `bin/utils_impacts.py:202`, `bin/mutgenomes_driver_priority.py:30,43,134,259`, `bin/plot_depths.py:687`) swallow all errors including `KeyboardInterrupt`.
- Safe modification: Replace bare `except:` with narrow exception types; validate expected columns at read time.
- Test coverage: Only `bin/test/` covers `utils_filter`, `mask_matrix`, `check_samplesheet`, `check_contamination`, `plot_selectionsideplots` — none of the fragile files above.

**Monolithic plotting scripts:**
- Files: `bin/plot_gene_saturation.py` (1546 lines, with `# TODO:`/`# FIXME` at lines 30, 1234, 1250, 1333), `bin/mutgenomes_expected_mutrisk.R` (836 lines, TODO at 49, FIXME at 218), `bin/plot_selectionsideplots.py` (796 lines), `bin/plot_depths.py` (771 lines).
- Why fragile: Single-file scripts mixing parsing, computation, and plotting; hard to test in isolation.
- Safe modification: Extract pure-computation helpers (pattern already proven by `bin/utils_filter.py` + `bin/test/test_utils_filter.py`).
- Test coverage: None for these scripts.

**Giant central workflow:**
- Files: `workflows/deepcsa.nf` (798 lines, 89 `params.` references, ~40 subworkflow includes with repeated aliasing of the same subworkflow, e.g. `MUTATION_DENSITY` included 5×, `MUTATIONAL_PROFILE` 4×, `OMEGA_ANALYSIS` 4×, `SIGNATURES` 4×).
- Why fragile: Any channel rename ripples across all aliased invocations; conditional logic (`if (run_profile_all)` etc.) is hard to test.
- Safe modification: Keep changes local to one aliased block; verify with `tests/deepcsa.nf.test` pipeline-level test after edits.

**Config-level failure masking:**
- Files: `conf/tmp_quick_fixes.config` (entire file is patches: `DNA2PROTEINMAPPING` retry-then-ignore, `COMPARE_SIGNATURES` ignore, `EXPECTEDMUTATEDCELLS` ignore), `conf/base.config:201-222`, `conf/exome.config:136-147`.
- Why fragile: "Quick fixes" have become permanent infrastructure; failures are invisible.
- Safe modification: Before removing an ignore, confirm the upstream bug is fixed; add explicit warning logs when a process is skipped.

## Scaling Limits

**Whole-cohort in-memory operations:**
- Current capacity: Cohort-level MAF/depth tables processed as single DataFrames; `DNA2PROTEINMAPPING` already needs 30 GB.
- Limit: Cohorts with large panels/WGS-scale site counts will exceed typical node memory in annotation and plotting steps.
- Scaling path: Roll out `panel_sites_chunk_size` chunking (implemented for `SITESFROMPOSITIONS`, `conf/modules.config:723`) to other processes; shard per-chromosome where possible.

**Test runtime:**
- Current capacity: Process-level nf-tests take ~20-80 s each (per `tests/2026-10-09_results.csv`); the pipeline-level test (`tests/deepcsa.nf.test`) runs full `main.nf`.
- Limit: Adding tests for all 48 `modules/local/` processes at current per-test cost makes CI slow.
- Scaling path: Use minimal test data (`tests/test_data/` is only 20 KB — good), parallelize nf-test, and consider tagging slow tests.

## Dependencies at Risk

**`latest`-tagged method containers (see Tech Debt above):**
- Risk: `dnds`, `oncodrivefml`, `oncodriveclustl`, `musical`, `msighdp`, `deepcsa_bed`, `expected_mutrate`, `bbgregressions:dev`, plus three untagged images.
- Impact: Silent behavioral change of core selection/signature methods.
- Migration plan: Pin versions; publish versioned images for `hdp_wrapper`, `test_mutated_genomes`, `sigprofilerassignment`.

**nf-core modules nearly absent:**
- Risk: `modules.json` tracks only 2 nf-core modules (`custom/dumpsoftwareversions`, `multiqc`) at a single git SHA (`911696ea0b62df80e900ef244d7867d177971f73`); everything else is local code in `modules/local/` (48 processes) and `subworkflows/local/` (24 subworkflows).
- Impact: No upstream fixes/updates flow in; local modules lack `meta.yml` in many cases (e.g. `modules/local/blacklistmuts/` has only `main.nf`).
- Migration plan: Incrementally port stable local processes to nf-core module standards (container pinning + `meta.yml` + tests).

**Mixed pandas/polars codebase:**
- Risk: pandas pinned mentally to <2.2.3 (TODOs in `bin/concat_sbs_probs.py`, `bin/mut_density_simple.py`) while newer scripts use polars.
- Impact: Upgrade friction; two APIs to maintain.
- Migration plan: Complete pandas 2.2.3 bump, then standardize new development on polars for large tables.

## Missing Critical Features

**No CI pipeline:**
- Problem: `.nf-core.yml` skips all GitHub workflows; nothing runs tests or linting on push/PR.
- Blocks: Reliable regression detection (5 tests are already failing unnoticed in local runs).

**No per-module tests for the vast majority of local modules:**
- Problem: 48 processes in `modules/local/` but only 1 `*.nf.test` file exists in the tree (`modules/local/expand_regions/tests/main.nf.test`) plus 2 tracked elsewhere; `subworkflows/local/` has none (only nf-core vendored subworkflows have tests).
- Blocks: Safe refactoring of `workflows/deepcsa.nf` and local modules.

**Pipeline-level test coverage of feature flags:**
- Problem: `tests/deepcsa.nf.test` covers basic run, omega, and MAF-input validation; the many boolean feature switches in `nextflow.config` (e.g. `oncodrive3d`, `dnds`, `indels`, `signatures`, `omega_covariates`, `downsample`, `regressions`, `contamination`) lack dedicated pipeline tests.
- Blocks: Confidence that optional branches still work after workflow edits.

## Test Coverage Gaps

**`bin/` Python scripts:**
- What's not tested: ~81 of 86 scripts have no unit test; only `bin/test/test_check_samplesheet.py`, `bin/test/test_check_contamination.py`, `bin/test/test_mask_matrix.py`, `bin/test/test_plot_selectionsideplots.py`, `bin/test/test_utils_filter.py` exist (unittest style, `sys.path.insert` sibling imports).
- Files: `bin/*.py`, `bin/test/`
- Risk: Silent numeric errors in selection statistics (omega/dNdS QC scripts) and plotting code.
- Priority: High for `bin/omega_syn_qc.py`, `bin/postprocessing_annotation.py`, `bin/utils_impacts.py`, `bin/mutgenomes_driver_priority.py` (all contain bare excepts or flagged-fragile logic).

**Local Nextflow modules:**
- What's not tested: 47 of 48 `modules/local/` processes; all 24 `subworkflows/local/` subworkflows.
- Files: `modules/local/`, `subworkflows/local/`
- Risk: Regressions like the 5 currently failing tests recur undetected.
- Priority: High for `filtermaf`, `group_genes`, `mut_density`, `sig_matrix_concat`, `sitesfrompositions` (currently failing); Medium for processes on the default execution path (`createpanels`, `mutationpreprocessing`).

**R scripts:**
- What's not tested: `bin/dNdS_run.R`, `bin/mutgenomes_expected_mutrisk.R` (836 lines with FIXME), `bin/mutrate_genome_trinuc_corrected.R`, `bin/signatures_msighdp_run.R` have no test harness.
- Files: `bin/*.R`
- Risk: Statistical errors in expected-mutation-risk and dN/dS computations go unnoticed.
- Priority: Medium.

---

*Concerns audit: 2026-10-09*
