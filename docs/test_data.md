# bbglab/deepCSA: Test data

deepCSA ships a minimal nf-test suite ([tests/deepcsa.nf.test](../tests/deepcsa.nf.test)) that exercises the main input scenarios and validation paths. This document describes where the test data lives, what it contains, and how it is consumed by the tests.

## Where the test data lives

The reference test datasets are hosted in the [bbglab/DeepClone_protocol](https://github.com/bbglab/DeepClone_protocol) repository, under `test_datasets/deepCSA/testdata/`:

```
test_datasets/deepCSA/testdata/
├── maf/
│   └── all_samples.somatic.mutations.maf      # cohort-level MAF (3 samples)
├── depth/
│   └── all_samples_indv.depths.tsv.gz         # precomputed per-position depths table
└── input_vcfs/
    ├── P19_0002_BDO_01.vcf
    ├── P19_0002_BTR_01.vcf
    └── P19_0003_BDO_01.vcf
```

The three test samples (`P19_0002_BDO_01`, `P19_0002_BTR_01`, `P19_0003_BDO_01`) come from a bladder duplex-sequencing experiment and are large enough to exercise the panel, mutational-profile, depth, and omega code paths while keeping runtimes short.

Locally committed inputs under [tests/test_data/](../tests/test_data/) only contain the small CSV samplesheets and one toy MAF used by the validation-failure tests:

| File | Purpose |
|---|---|
| `input.csv` | Samplesheet with `sample,vcf,bam` columns referring to internal IRB paths (not used by the public CI tests). |
| `input_maf.csv` | Samplesheet with `sample,vcf` columns pointing to remote VCFs from `bbglab/DeepClone_protocol`. Used by the MAF-input test. |
| `input_no_bam.csv` | Same as `input_maf.csv`, used by the VCF-without-BAM tests. |
| `test_mutations.maf` | Tiny MAF used only by the parameter-validation failure tests (3, 4, 5). |

## Remote-fetching convention

Following the convention used by nf-core pipelines (e.g. [nf-core/fastquorum](https://github.com/nf-core/fastquorum)), the MAF and depths files are fetched **at runtime** directly from `bbglab/DeepClone_protocol` rather than pre-downloaded:

```groovy
input_maf           = 'https://raw.githubusercontent.com/bbglab/DeepClone_protocol/main/test_datasets/deepCSA/testdata/maf/all_samples.somatic.mutations.maf'
use_custom_depths   = true
custom_depths_table = 'https://raw.githubusercontent.com/bbglab/DeepClone_protocol/main/test_datasets/deepCSA/testdata/depth/all_samples_indv.depths.tsv.gz'
```

> ⚠️ The `bbglab/DeepClone_protocol` repository must remain **publicly accessible** for Nextflow to fetch these files at runtime. If access is restricted the tests fail with a "No such file or directory" error.

Because nf-schema 2.x validates `file-path` parameters for local existence, [tests/nextflow.config](../tests/nextflow.config) excludes the remote-URL parameters from that check:

```groovy
validation {
    ignoreParams = ['input_maf', 'custom_depths_table']
}
```

## How tests map to input scenarios

The five nf-test cases cover all three [input scenarios](input_scenarios.md) plus three validation-failure paths:

| Test | Scenario covered | Inputs |
|---|---|---|
| TEST 1 — basic MAF processing | Scenario 3 (cohort MAF + depths) | `input_maf.csv` + remote MAF + remote depths |
| TEST 1b — VCF + depths | Scenario 2 (VCF + precomputed depths) | `input_no_bam.csv` + remote depths |
| TEST 2 — omega run | Scenario 3 with `omega = true` | same as TEST 1 |
| TEST 3 — `--input_maf` without `--use_custom_depths` | Validation failure | `input_maf.csv` + local toy MAF |
| TEST 4 — VCF samplesheet without BAMs and `use_custom_depths = false` | Validation failure | `input_no_bam.csv` |
| TEST 5 — `use_custom_depths = true` without a depths table | Validation failure | `input_no_bam.csv` |

Snapshots (MD5 of selected outputs) are stored in [tests/deepcsa.nf.test.snap](../tests/deepcsa.nf.test.snap).

## Running and updating the tests

For the full execution-environment notes (SLURM, Singularity, `DEEPCSA_TEST_WORKDIR`, configuring for a non-IRB site) see [tests/README.md](../tests/README.md). The short version:

```bash
nf-test test tests/deepcsa.nf.test                                      # run the whole suite
nf-test test tests/deepcsa.nf.test --tag omega                          # run a single test
nf-test test tests/deepcsa.nf.test --update-snapshot                    # regenerate snapshots
```

The tests must be run on a SLURM cluster with Singularity; local execution is not supported because resource limits won't be met. Snapshots must be regenerated whenever default pipeline parameters change.
