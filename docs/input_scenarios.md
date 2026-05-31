# bbglab/deepCSA: Input scenarios

deepCSA supports three input scenarios depending on what you already have available (BAMs, mutations as VCFs, or a cohort-level MAF together with a precomputed depths table). All scenarios still require the standard samplesheet CSV passed via `--input`.

Sample naming rules apply to every scenario: avoid `.` in sample names and prefer text-like names instead of purely numeric ones. See [File formatting](file_formatting.md) for details on each file.

## Scenario summary

| Scenario | `--input` columns | Depth source | Extra flags |
|---|---|---|---|
| 1. VCF + BAM (default) | `sample,vcf,bam` | Computed from BAMs | — |
| 2. VCF + precomputed depths | `sample,vcf` | `--custom_depths_table` | `--use_custom_depths true` |
| 3. Cohort MAF + precomputed depths | `sample,vcf` (metadata only) | `--custom_depths_table` | `--input_maf <maf>` + `--use_custom_depths true` |

The pipeline validates these combinations at start-up and stops with an explicit error if `--input_maf` is set without `--use_custom_depths true` (see [workflows/deepcsa.nf](../workflows/deepcsa.nf)).

## Scenario 1 — VCF + BAM (default)

Use this scenario when you have per-sample variant calls and the BAM files that were used to produce them.

```csv
sample,vcf,bam
sample1,sample1.filtered.vcf,sample1.sorted.bam
sample2,sample2.filtered.vcf,sample2.sorted.bam
```

The pipeline derives per-position sequencing depth directly from the BAMs (subworkflow `depthanalysis`). No extra flag is needed.

## Scenario 2 — VCF + precomputed depths

Use this scenario when you already have a depths table (for example produced by a previous deepCSA run, or by an external tool) and you want to skip BAM-based pileup.

```csv
sample,vcf
sample1,sample1.filtered.vcf
sample2,sample2.filtered.vcf
```

```console
params {
    use_custom_depths   = true
    custom_depths_table = '/path/to/precomputed_depths_table.tsv'
}
```

Notes:

- The depths-table column names must match the sample names declared in the `sample` column of the input CSV.
- `custom_depths_table` may be TSV or CSV but must follow the per-position depth layout that deepCSA expects.
- If the file is missing or unreadable the pipeline fails immediately.

See [Usage — Using a precomputed depths table](usage.md#using-a-precomputed-depths-table) for additional notes on how columns are matched and on preparing the file from a previous deepCSA run.

## Scenario 3 — Cohort MAF + precomputed depths

Use this scenario when all mutations for the cohort are already consolidated in a single MAF/TSV file and you also have the matching precomputed depths table.

```console
params {
    input                = "samplesheet.csv"
    input_maf            = "cohort_mutations.maf"
    use_custom_depths    = true
    custom_depths_table  = "precomputed_depths.tsv"
}
```

```bash
nextflow run bbglab/deepCSA \
    --input        samplesheet.csv \
    --outdir       results/ \
    --input_maf    cohort_mutations.maf \
    --use_custom_depths    true \
    --custom_depths_table  precomputed_depths.tsv \
    -profile <DESIRED_PROFILE>
```

What happens under the hood:

1. The MAF file is split into one VCF per unique `SAMPLE_ID` by `INPUTMAF2VCF` (script [assets/useful_scripts/deepcsa_maf2samplevcfs.py](../assets/useful_scripts/deepcsa_maf2samplevcfs.py)).
2. The per-sample VCFs are published under `<outdir>/processing_files/input_vcfs/`.
3. The rest of the pipeline runs as in Scenario 2.

The standard `--input` samplesheet is still required, because it provides the sample metadata used by other pipeline steps. The `SAMPLE_ID` values in the MAF must match the `sample` column of the samplesheet.

For the expected MAF columns (deepCSA-generated MAF vs external MAF) see [Usage — MAF file format](usage.md#maf-file-format).

## Related parameters

| Parameter | Purpose |
|---|---|
| `input` | Samplesheet CSV with `sample,vcf[,bam]` columns. Always required. |
| `input_maf` | Cohort-level MAF file (Scenario 3). Requires `use_custom_depths = true`. |
| `use_custom_depths` | Skip BAM-based depth computation. Required for Scenarios 2 and 3. |
| `custom_depths_table` | Path to the precomputed per-position depths table. Required when `use_custom_depths = true`. |

Custom-mutation workflows (e.g. forcing your own filter list) layer on top of these scenarios. See [Usage — Custom mutation calls](usage.md#custom-mutation-calls----option-1-building-input-vcfs-and-providing-them-via-normal-input) for the advanced options.
