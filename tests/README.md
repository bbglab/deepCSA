## DEEPCSA Tests

This directory contains reproducible tests for the DEEPCSA Nextflow pipeline using nf-test. Tests run on the IRB cluster via SLURM and cover a minimal-profile run, an omega-enabled run with snapshot validation, and a parameter validation failure case.

### Overview
- Framework: nf-test (>= 0.9.2), plugin `nft-utils@0.0.3`
- Orchestrator: Nextflow
- Executor: SLURM (configured in `tests/nextflow.config`)
- Container engine: Singularity
- Snapshot store: `tests/deepcsa.nf.test.snap`
- Test spec: `tests/deepcsa.nf.test`
- nf-test work directory: configured in `nf-test.config` (IRB default: `/scratch/bbg/work/deepCSA/nf-tests` — customize for your environment)  
- Test output directory: configured in `tests/nextflow.config` (IRB default: `/scratch/bbg/work/deepCSA/nf-tests-outputs` — customize for your environment)  

### Structure
- `nf-test.config`: nf-test root configuration (work dir, profile, ignored paths)
- `tests/deepcsa.nf.test`: nf-test suite with three active test cases
  - "TEST 1. Basic functionality - MAF-based processing" (tag: `normal`)  
  - "TEST 2.Omega analysis test run" (tag: `omega`)
  - "TEST 3. Should fail when --input_maf is provided without --use_custom_depths" (tag: `input_maf_validation`)
- `tests/deepcsa.nf.test.snap`: Snapshot file managed by nf-test
- `tests/nextflow.config`: Test-specific Nextflow config (SLURM, Singularity, cluster paths, param overrides)
- `tests/test_data/`: Input dataset used by the tests (points to cluster paths; adapt `tests/test_data/input.csv` for your environment — see [Configuring for usercase](#configuring-for-usercase))

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

- Update snapshots after intentional output changes:
```bash
nf-test test tests/deepcsa.nf.test --update-snapshot
```

- Update only one test:
```bash
nf-test test tests/deepcsa.nf.test --tag omega --update-snapshot
```

> ⚠️ Always run snapshot regeneration from the cluster (not locally). Snapshots encode MD5 hashes of pipeline outputs — the pipeline must complete successfully before updating them. After updating, review the new hashes in `tests/deepcsa.nf.test.snap` and commit the file.

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
  - Structural checks on `all_omegas.tsv`: header contains `gene`, `sample`, `dnds`; all rows contain same columns, all samples are present.
  - Snapshot of `mutational_profile/all_samples.all.profile.tsv` (MD5)

- `input_maf_validation` (Parameter validation):
  - Pipeline **fails** when `--input_maf` is set without `--use_custom_depths = true`

Note on non-deterministic outputs: omega metrics contain floating-point values that can vary slightly across runs/environments. For now, we assert file structure and line count. When needed, we will switch to content checks with numeric tolerance (e.g., rounding selected columns before comparison) to keep tests robust while validating semantics.

### Cleaning Outputs
Test outputs go to the directory set in `params.outdir` inside `tests/nextflow.config`. If you want a clean run:
```bash
rm -rf /path/to/nf-tests-outputs/tests_results_normal
rm -rf /path/to/nf-tests-outputs/tests_results_omega
```

nf-test working files (Nextflow work dirs, logs) are under the path set as `workDir` in `nf-test.config`. These can be deleted safely between runs:
```bash
rm -rf /path/to/nf-tests-workdir
```

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

#### 1. `nf-test.config` — work directory
Change `workDir` to a scratch or fast-storage path on your cluster:
```groovy
workDir "/path/to/nf-tests-workdir"
```

#### 2. `tests/nextflow.config` — executor and output directory
Update the SLURM queue and the test output path:
```groovy
process {
    executor = 'slurm'
    queue    = 'your_queue_name'   // replace with your cluster's queue
    ...
}
params {
    outdir = "/path/to/nf-tests-outputs"
    ...
}
```

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
| `dnds_ref_transcripts` / `dnds_covariates` | dNdScv reference files |
| `datasets3d` / `annotations3d` | Oncodrive3D datasets |
| `singularity.cacheDir` / `singularity.libraryDir` | Singularity image cache |

#### 4. `tests/test_data/input.csv` — input data paths
This file lists per-sample VCF paths. The paths need to be added for the pipelie to work but the file used is the MAF file for all samples, added via `tests/nextflow.config`:
```csv
sample,maf
P19_0001_BDO_01,/path/to/test_datasets/vcfs/P19_0001_BDO_01.maf
...
```
In the same direction, update the depths file path in `tests/nextflow.config`:
```groovy
custom_depths_table         = "path/to/all_samples_indv.depths.tsv.gz"
input_maf                   = "path/to/all_samples.maf"
```

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
