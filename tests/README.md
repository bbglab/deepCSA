## DEEPCSA Tests

This directory contains reproducible tests for the DEEPCSA Nextflow pipeline using nf-test. Tests run on the IRB cluster via SLURM and cover a minimal-profile run, an omega-enabled run with snapshot validation, and a parameter validation failure case.

### Overview
- Framework: nf-test (>= 0.9.2), plugin `nft-utils@0.0.3`
- Orchestrator: Nextflow
- Executor: SLURM (configured in `tests/nextflow.config`)
- Container engine: Singularity
- Snapshot store: `tests/deepcsa.nf.test.snap`
- Test spec: `tests/deepcsa.nf.test`
- nf-test work directory: set via `DEEPCSA_TEST_WORKDIR` env var (default: `.nf-test/` in the repo root)
- Test output directory: managed automatically by nf-test per-test via `$outputDir` (derived from `workDir`; no manual configuration needed)

### Structure
- `nf-test.config`: nf-test root configuration (work dir, profile, ignored paths)
- `tests/deepcsa.nf.test`: nf-test suite with three active test cases
  - "TEST 1. Basic functionality - MAF-based processing" (tag: `normal`)  
  - "TEST 2.Omega analysis test run" (tag: `omega`)
  - "TEST 3. Should fail when --input_maf is provided without --use_custom_depths" (tag: `input_maf_validation`)
- `tests/deepcsa.nf.test.snap`: Snapshot file managed by nf-test
- `tests/nextflow.config`: Test-specific Nextflow config (SLURM, Singularity, cluster paths, param overrides)
- `tests/test_data/`: Local committed input files used by the tests (`input.csv`, `test_mutations.maf`). The MAF and depths files for the main tests are fetched at runtime directly from [bbglab/DeepClone_protocol](https://github.com/bbglab/DeepClone_protocol) — no download step required.

### Running Tests

> **Prerequisite:** Tests require a SLURM cluster with Singularity available. All jobs are submitted to SLURM automatically — do not run from a local machine, as resource limits will not be met. If you are not on the IRB cluster, follow [Configuring for usercase](#configuring-for-usercase) first.

Run the whole suite:
```bash
nf-test test tests/deepcsa.nf.test
```

Run a single test by tag:
```bash
nf-test test tests/deepcsa.nf.test --tag normal
nf-test test tests/deepcsa.nf.test --tag omega
nf-test test tests/deepcsa.nf.test --tag input_maf_validation
```

The `tests/nextflow.config` is loaded automatically via `nf-test.config`. It handles SLURM submission, Singularity, and all cluster-specific paths — no extra flags needed.

### Snapshots
Snapshots are stored in `tests/deepcsa.nf.test.snap` and encode expected outputs for deterministic comparison.

- Update snapshots after intentional output or parameter changes:
```bash
nf-test test tests/deepcsa.nf.test --update-snapshot
```

- Update only one test:
```bash
nf-test test tests/deepcsa.nf.test --tag omega --update-snapshot
```

> ⚠️ Always run snapshot regeneration from the cluster (not locally). Snapshots encode MD5 hashes of pipeline outputs — the pipeline must complete successfully before updating them. After updating, review the new hashes in `tests/deepcsa.nf.test.snap` and commit the file.
> ⚠️ The tests are using **default parameters for the pipeline**. If the default parameters are modified, it is important to run snapshot regeneration.

### Current Assertions per Test
- `normal` (Minimal features):
  - Pipeline succeeds
  - `mutational_profile/` output directory exists
  - `mutdensity/`, `omega/`, `oncodrivefml/`, `oncodrive3d/` directories do **not** exist
  - Snapshot of `mutational_profile/all_samples.all.profile.tsv` (MD5)

- `omega` (Omega analysis):
  - Pipeline succeeds
  - `mutational_profile/` and `omega/` directories exist
  - `omega/all_omegas.tsv` exists
  - `oncodrivefml/`, `oncodrive3d/` directories do **not** exist
- Structural checks on `all_omegas.tsv`: header contains `gene`, `sample`, `dnds`, `pvalue_adj`; all rows contain same columns, all samples are present.
  - Snapshot of `mutational_profile/all_samples.all.profile.tsv` (MD5)

- `input_maf_validation` (Parameter validation):
  - Pipeline **fails** when `--input_maf` is set without `--use_custom_depths = true`

Note on non-deterministic outputs: omega metrics contain floating-point values that can vary slightly across runs/environments. For now, we assert file structure and line count. When needed, we will switch to content checks with numeric tolerance (e.g., rounding selected columns before comparison) to keep tests robust while validating semantics.

### Cleaning Outputs
Test outputs and nf-test working files are both stored under `workDir` (`DEEPCSA_TEST_WORKDIR`, or `.nf-test/` by default). To wipe everything:
```bash
rm -rf "${DEEPCSA_TEST_WORKDIR:-.nf-test}"
```
Individual test output directories live inside that path under `tests/<TEST_ID>/`.

### Debugging Failures
When a process fails under nf-test, the output includes a work directory like:
```
<workDir>/tests/<TEST_ID>/work/<HASH>/<HASH>
```
where `<workDir>` is the path set in `nf-test.config`.
Investigate there:
```bash
cd <workDir>/tests/<TEST_ID>/work/<HASH>/<HASH>
cat .command.out
cat .command.err
cat .command.sh
```
Reproduce the exact command environment:
```bash
bash .command.run
```

Also check the test’s Nextflow log:
```bash
cat <workDir>/tests/<TEST_ID>/meta/nextflow.log
```

### Configuring for usercase

The tests are set up for the IRB cluster but can be adapted to any SLURM environment. There are four things to change:

#### 1. Work directory — `DEEPCSA_TEST_WORKDIR`
Set the `DEEPCSA_TEST_WORKDIR` environment variable to a scratch or fast-storage path on your cluster before running tests:
```bash
export DEEPCSA_TEST_WORKDIR="/scratch/your_site/nf-tests"
```
If the variable is not set, nf-test defaults to `.nf-test/` in the repo root (already covered by `.gitignore`). No file editing required.

#### 2. `tests/nextflow.config` — executor and queue
Update the SLURM queue in `tests/nextflow.config`:
```groovy
process {
    executor = 'slurm'
    queue    = 'your_queue_name'   // replace with your cluster's queue
    ...
}
```
Test output directories are managed automatically by nf-test via `$outputDir` and do not need to be configured manually.

#### 3. `tests/nextflow.config` — reference file paths
The file includes `../conf/general_files_IRB.config`, which supplies all cluster-specific reference paths (genome FASTA, VEP cache, CADD scores, Singularity container cache, etc.). Either:
- Create your own `conf/general_files_<yoursite>.config` with the equivalent paths and point `tests/nextflow.config` to it, or
- Override each path directly in the `params {}` block of `tests/nextflow.config`.

The parameters you will need to provide are:

| Parameter | Description |
|---|---|
| `fasta` | GRCh38 reference genome FASTA |
| `vep_cache` | Ensembl VEP cache directory |
| `cadd_scores` / `cadd_scores_ind` | CADD scores TSV + index |
| `cosmic_ref_signatures` | COSMIC SBS signatures file |
| `nanoseq_snp` / `nanoseq_noise` | NanoSeq masking BED files |
| `dnds_biomart_ref` / `dnds_covariates` | dNdScv reference files |
| `datasets3d` / `annotations3d` | Oncodrive3D datasets |
| `singularity.cacheDir` / `singularity.libraryDir` | Singularity image cache |

#### 4. Test data — no download required
The MAF file and precomputed depths table are fetched **at runtime** directly from [bbglab/DeepClone_protocol](https://github.com/bbglab/DeepClone_protocol). Following the same convention used by nf-core pipelines (e.g. [nf-core/fastquorum](https://github.com/nf-core/fastquorum)), the remote URLs are set directly in the test params blocks inside `tests/deepcsa.nf.test`:

> ⚠️ **Requirement:** The `bbglab/DeepClone_protocol` repository must be **publicly accessible** for Nextflow to fetch the files at runtime. If the repository is private, the pipeline will fail with a "No such file or directory" error for the remote URLs.

```groovy
input_maf           = 'https://raw.githubusercontent.com/bbglab/DeepClone_protocol/main/test_datasets/deepCSA/testdata/maf/all_samples.somatic.mutations.maf'
use_custom_depths   = true
custom_depths_table = 'https://raw.githubusercontent.com/bbglab/DeepClone_protocol/main/test_datasets/deepCSA/testdata/depth/all_samples_indv.depths.tsv.gz'
```

Because nf-schema 2.x validates `file-path` params for local existence, `tests/nextflow.config` includes:
```groovy
validation {
    ignoreParams = ['input_maf', 'custom_depths_table']
}
```
This tells nf-schema to skip the file-existence check for these two params when running tests, allowing HTTP URLs to pass through to Nextflow's native remote file handling. No manual path editing or pre-download is required.

### Adding a New Test
1. Duplicate a `test("...")` block in `tests/deepcsa.nf.test`.
2. Give it a new `tag` and configure `when { params { ... } }`.
3. Add assertions in `then { ... }`.
4. Run the test and create/update snapshots:
   ```bash
   nf-test test tests/deepcsa.nf.test --tag <your-tag> --update-snapshot
   ```

### Troubleshooting
- Snapshot mismatches: Inspect the diff in nf-test output and update snapshots if the change is expected.
- Header-related errors: Ensure the pipeline stage writes headers with `keepHeader: true` when concatenating files.
- Containers/paths: Validate container image availability and local paths used in `nextflow.config`.

### Versions
The snapshot records tool versions under `meta`. Example (as of writing):
- nf-test: 0.9.5
- Nextflow: 25.10.4
