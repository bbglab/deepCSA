## DEEPCSA Tests

This directory contains reproducible tests for the DEEPCSA Nextflow pipeline using nf-test. It covers a minimal-profile run and an omega-enabled run with snapshot validation.

### Overview
- Framework: nf-test (>= 0.9.2)
- Orchestrator: Nextflow
- Snapshot store: `tests/deepcsa.nf.test.snap`
- Test spec: `tests/deepcsa.nf.test`

### Structure
- `tests/deepcsa.nf.test`: nf-test suite with two test cases
  - "Minimal features test run" (tag: `normal`)
  - "Omega analysis test run" (tag: `omega`)
- `tests/deepcsa.nf.test.snap`: Snapshot file managed by nf-test
- `tests/test_data/`: Small input dataset used by the tests
- Output (created by tests):
  - `${projectDir}/tests_results_normal`
  - `${projectDir}/tests_results_omega`

### Running Tests
Run the whole suite:
```bash
nf-test test tests/deepcsa.nf.test
```

Run a single test by tag:
```bash
nf-test test tests/deepcsa.nf.test --tag normal
nf-test test tests/deepcsa.nf.test --tag omega
```

Use the repo’s nextflow config (default in test spec):
```bash
nf-test test tests/deepcsa.nf.test -c nextflow.config
```

### Snapshots
Snapshots are stored in `tests/deepcsa.nf.test.snap` and encode expected outputs for deterministic comparison.

- Update snapshots after intentional output changes:
```bash
nf-test test tests/deepcsa.nf.test --update-snapshot
# or older CLI may require: --updateSnapshot
```

- Update only one test:
```bash
nf-test test tests/deepcsa.nf.test --tag omega --update-snapshot
```

### Current Assertions per Test
- `normal` (Minimal features):
  - Pipeline succeeds
  - Only compute profile artifacts are produced
  - Snapshot of `all_samples.all.profile.tsv` (MD5)

- `omega` (Omega analysis):
  - Pipeline succeeds
  - Omega directory and `all_omegas.tsv` exist
  - Structural checks on `all_omegas.tsv` header and row count (`lines.size() == 59`)
  - Snapshot of the deterministic profile file (MD5)

Note on non-deterministic outputs: omega metrics contain floating-point values that can vary slightly across runs/environments. For now, we assert file structure and line count. When needed, we will switch to content checks with numeric tolerance (e.g., rounding selected columns before comparison) to keep tests robust while validating semantics.

### Cleaning Outputs
If you want a clean run:
```bash
rm -rf tests_results_normal tests_results_omega
```

### Debugging Failures
When a process fails under nf-test, the output includes a work directory like:
```
.nf-test/tests/<TEST_ID>/work/<HASH>/<HASH>
```
Investigate there:
```bash
cd .nf-test/tests/<TEST_ID>/work/<...>
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
cat .nf-test/tests/<TEST_ID>/meta/nextflow.log
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
- nf-test: 0.9.2
- Nextflow: 25.04.6
