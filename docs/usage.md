# bbglab/deepCSA: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Table of contents

- [Introduction](#introduction)
- [How to run the pipeline](#how-to-run-the-pipeline)
- [Input scenarios](#input-scenarios)
- [Samplesheet input](#samplesheet-input)
- [Available genomes](#available-genomes)
- [Proposed run modes](#proposed-run-modes)
  - [Initial run. Data exploration](#initial-run-data-exploration)
  - [Clonal structure definition](#clonal-structure-definition-complete-run-with-a-focus-on-positive-selection-at-the-cohort-level)
  - [Mutational processes in alternative genomic regions](#mutational-processes-in-alternative-genomic-regions-partial-run-with-a-focus-on-mutational-processessignatures)
  - [Interindividual variability and sample comparison](#interindividual-variability-and-sample-comparison-complete-run-with-downstream-steps-for-computation-of-linear-regressions-to-compare-different-samplesgroups-based-on-clinical-variables-or-sample-metadata)
- [Definition of structural parameters](#definition-of-structural-parameters)
- [Additional customizable parameters](#additional-customizable-parameters)
- [Custom mutation calls](#custom-mutation-calls)
- [MAF file as input (alternative input mode)](#maf-file-as-input-alternative-input-mode)

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## How to run the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run bbglab/deepCSA --outdir <OUTDIR> -profile <DESIRED PROFILE> --input samplesheet.csv
```

For more information on how to run Nextflow pipelines check a more detailed explanation [below](#running-the-pipeline) in this same document or check the [Nextflow](https://www.nextflow.io/docs/latest/index.html) or [nf-core](https://nf-co.re) community documentations.

## Input scenarios

deepCSA accepts three different input combinations: per-sample VCF + BAM (default), per-sample VCF + a precomputed depths table, or a cohort-level MAF + a precomputed depths table. The sections below describe each piece in detail; for a concise summary of the three modes and when to use each, see [Input scenarios](input_scenarios.md).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

Example:

```csv
sample,vcf,bam
sample1,sample1.high.filtered.vcf,sample1.sorted.bam
sample2,sample2.high.filtered.vcf,sample2.sorted.bam
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Sample names cannot contain dots (`.`). Ideally the sample name should have a _Python string-like_ format, it should not be a single number.  |
| `vcf` | Full path to VCF file containing all the mutations called in your sample. It should be uncompressed and with the VCF format field complying with the expected format. See [custom mutation calling](#custom-mutation-calls) below in case the input is not coming from deepUMIcaller.                                                             |
| `bam` | Full path to BAM file containing the duplex aligned reads that were used for the variant calling.                                                             |

An [example samplesheet](../assets/example_inputs/input_example.csv) has been provided with the pipeline.

### Input files from deepUMIcaller

If your mutations were called with [deepUMIcaller](https://github.com/bbglab/deepUMIcaller), use the **final** per-sample outputs produced by that pipeline:

- **VCF**: the filtered duplex VCF for each sample (commonly named `*.duplex.filtered.vcf` or `*.high.filtered.vcf`), typically found under the `mutations_vcf/` output directory. The VCF should be uncompressed.
- **BAM**: the **duplex-consensus** BAM used for calling those variants (commonly named `*.sorted.bam` in `sortbamduplexcons/`). This BAM must be aligned to the same reference genome as the VCF.

Batch/sample guidance:

- **One row per library/run**: if you have multiple sequencing libraries for the same biological sample and want to **aggregate** them, keep the same `sample` name across rows and list each matching VCF/BAM pair.
- **Separate batches**: if you want to compare batches or runs, keep distinct sample names and (optionally) add a `BATCH` or `RUN_ID` column in the feature groups table to stratify the analysis.

## Available genomes

deepCSA pipeline heavily relies on bgreference and bgdata tools so the use of this pipeline is limited to those genomes available in these packages. In particular, the default containers that are being used already have the hg38 and mm39 genomes cached, if you want to use any other genome, open an issue and we will address it as soon as we can.

## Proposed run modes

These are 4 different ways of running the pipeline, each of them serving for a specific purpose, a list of expected outputs is provided in each run mode section.

### Initial run. Data exploration

It will provide:

- Definition of regions to analyze
- Depth per sample and/or per gene
- Somatic mutations
  - Needle plots
- Mutational profile
- Mutational signatures

```console
params {
    plot_depths   = true
    signatures    = true
    profileall    = true
}
```

### Clonal structure definition. Complete run with a focus on positive selection at the cohort-level

It will provide:

- All the previously described outputs plus...
- Mutation density
- Positive selection per gene
  multiple positive selection metrics
  - Per gene, all samples together
  - Per gene, per group of samples
  - Per gene, per sample

```console
params {
    mutationdensity             = true

    profileall                  = true

    oncodrivefml                = true
    oncodriveclustl             = true

    oncodrive3d                 = true
    o3d_raw_vep                 = true
    o3d_plot                    = true


    omega                       = true
    omega_multi                 = true
    omega_globalloc             = true
    omega_mutabilities          = true
    site_comparison_grouping    = 'all'
    omega_plot                  = true

    create_subgenic_regions     = true
    autodomains                 = true
    autoexons                   = true

    mutated_cells_vaf           = true
    mutepi_genes_to_recode      = null

    indels                      = true

    signatures                  = true
}
```

### Mutational processes in alternative genomic regions. Partial run with a focus on mutational processes/signatures

It will provide:

- Same as initial run (even it can be ignored)
- Mutational profile and mutational signatures based on:
  - All genomic regions
  - Only exonic regions
  - Only non-protein affecting regions (synonymous mutations and intronic, intergenic)
  - Intronic and intergenic regions

```console
params {
    mutationdensity             = true

    profileall                  = true
    profilenonprot              = true
    profileexons                = true
    profileintrons              = true

    signatures                  = true
}
```

### Interindividual variability and sample comparison. Complete run with downstream steps for computation of linear regressions to compare different samples/groups based on clinical variables or sample metadata

It will provide:

- Same as complete clonal structure definition +
- Computation of univariate and multivariate linear regressions between clonal structure metrics and clonal selection

```console
params {
    mutationdensity             = true

    profileall                  = true

    omega                       = true
    omega_multi                 = true
    omega_globalloc             = true
    omega_mutabilities          = true
    site_comparison_grouping    = 'all'
    omega_plot                  = true

    create_subgenic_regions     = true
    autodomains                 = true
    autoexons                   = true

    regressions                 = true
    // additional regression parameters, see nextflow_schema.json for more info
      ...
}
```

### Feature toggles (turning steps on/off)

All pipeline parameters are defined in `nextflow_schema.json` and exposed via `--help`. Most analysis steps can be enabled/disabled with boolean flags such as `mutationdensity`, `omega`, `oncodrive3d`, `signatures`, and `regressions`. If a step is turned off, its output directories are not produced. For a full list of parameters run:

```bash
nextflow run bbglab/deepCSA --help
```

## Definition of structural parameters

Before running, ensure container images can be pulled (or are already available) and that all reference datasets are accessible from your execution environment.

### Data sources and reference downloads

- **Ensembl VEP cache**  
  Download from the [Ensembl VEP cache](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) page and set `vep_cache`, `vep_cache_version`, `vep_species`, and `vep_genome` accordingly.

- **CADD scores**  
  Use the ["All possible SNVs" GRCh38/hg38 file](https://cadd.gs.washington.edu/download). Provide both the compressed TSV (`cadd_scores`) and its tabix index (`cadd_scores_ind`).

- **COSMIC signatures**  
  Download the SBS signatures (context size = 96) for your genome build from the [COSMIC signatures downloads page](https://cancer.sanger.ac.uk/signatures/downloads/). Set `cosmic_ref_signatures`.

- **dNdScv reference data**  
  Download `covariates_hg19_hg38_epigenome_pcawg.rda` from the [dNdScv](https://github.com/im3sanger/dndscv) reference data (also mirrored by [IntOGen](https://intogen.org/download)). Set `dnds_covariates`.
  Run the scripts for the generation of the dNdScv required reference input that you can find here: https://github.com/bbglab/deepCSA/tree/dev/assets/build_datasets/dndscv

- **Oncodrive3D datasets**  
  Build datasets and annotations following the [Oncodrive3D dataset instructions](https://github.com/bbglab/oncodrive3d?tab=readme-ov-file#building-datasets). Provide the resulting `datasets3d` and `annotations3d` directories.

- **Trinucleotide counts**  
  Provide a `wgs_trinuc_counts` file with the total count of each trinucleotide in your reference genome (see `assets/trinucleotide_counts/` for the expected format).

- **DNA2PROTEINMAPPING GFF3 (optional)**  
  By default, the pipeline fetches a GFF3 file from Ensembl FTP at runtime. For local/offline use, download the matching file from  
  `https://ftp.ensembl.org/pub/release-<release>/gff3/<species>/`  
  (e.g., `Homo_sapiens.GRCh38.111.gff3.gz`) and provide it via `gff3_file`.

- **Domain definitions**  
  Supply a Pfam/InterPro domain file (see [file formatting](file_formatting.md#domain-definition-file)).

- **NanoSeq masks (optional)**  
  See [Nanoseq genomic masks](#nanoseq-genomic-masks) below.

### Mandatory parameter configuration

See [File formatting docs](file_formatting.md) for more details on the structure of files that can be provided to deepCSA.

```console
params {

    fasta                      = null

    cosmic_ref_signatures      = "COSMIC_v3.4_SBS_GRCh38.txt"
    wgs_trinuc_counts          = "assets/trinucleotide_counts/trinuc_counts.homo_sapiens.tsv"

    // oncodrivefml (only for human; could be adapted to others)
    cadd_scores                = "CADD/v1.7/hg38/whole_genome_SNVs.tsv.gz"
    cadd_scores_ind            = "CADD/v1.7/hg38/whole_genome_SNVs.tsv.gz.tbi"

    // dnds
    // dnds_biomart_ref is a biomart TSV; deepCSA dynamically builds a per-run RefCDS_custom.rda
    // by intersecting it with the panel BED (replaces the previously required static
    // RefCDS_*.rda transcripts file). See assets/build_datasets/dndscv/instructions.txt
    // for how to regenerate the biomart export.
    dnds_biomart_ref           = "biomart_export.tsv"
    dnds_covariates            = "covariates_hg19_hg38_epigenome_pcawg.rda"

    // GFF3 annotation for the genome assembly, consumed when building exon/domain panels
    gff3_file                  = "Homo_sapiens.GRCh38.111.gff3.gz"

    // oncodrive3d + fancy plots
    datasets3d                 = "oncodrive3d/datasets"
    annotations3d              = "oncodrive3d/annotations"
    domains_file               = "pfam.tsv"


    vep_cache                  = ".vep"

    // Ensembl VEP for homo_sapiens, but should be adjusted accordingly to species and cache version
    vep_genome                 = "GRCh38"
    vep_species                = "homo_sapiens"
    vep_cache_version          = 111
    vep_out_format             = "tab"
    vep_params                 = "--no_stats --cache --offline --symbol --protein --canonical --af_gnomadg --af_gnomade"
    vep_params_panel           = "--no_stats --cache --offline --symbol --protein --canonical"
}
```

<!-- TODO: revise if this could be used in mice http://mendel.stanford.edu/SidowLab/downloads/gerp/ -->
### Optional parameters configuration

See [File formatting docs](file_formatting.md) for more details on the structure of files that can be provided to deepCSA.

```console
params {

    // definition of gene groups
    // could be fixed or dynamic based on the study
    custom_groups               = false
    custom_groups_file          = null
    custom_groups_separator     = 'tab'

    // customize the annotation of certain regions i.e. TERT promoter mutations, other non-coding drivers...
    customize_annotation        = false
    custom_annotation_tsv       = ''


    // define a set of common known hotspots
    hotspots_annotation         = false
    hotspots_definition_file    = ''


    // definition of specific regions within genes with specific interest on computing dN/dS
    subgenic_bedfile      = null

    // define a file of mutations that should not be trusted
    //  and you want to remove from all the analysis
    blacklist_mutations        = null
}
```

#### Nanoseq genomic masks

These files identify sites overlapping common SNPs and noisy or variable genomic regions, as described in [Abascal et al, 2021](https://www.nature.com/articles/s41586-021-03477-4) and used in the [Nanoseq pipeline](https://github.com/cancerit/NanoSeq). Two BED files are available to be used:

- Nanoseq SNP: Common SNP positions that should be excluded from analysis
- Nanoseq Noise: Regions with high noise or variability

Enable them with:

```console
params {
    nanoseq_snp   = "SNP_GRCh38.wgns.bed.gz"
    nanoseq_noise = "NOISE_GRCh38.wgns.bed.gz"
}
```

Both files are available for GRCh38 at the [shared folder](https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL) from Iñigo Martincorena's group, at the Wellcome Sanger Institute.

## Additional customizable parameters

In addition to several files that can be provided as input listed in the [optional files parameters](#optional-parameters-configuration), there are some more parameters that allow for specific tunnings of the analysis.

### Minimum depth thresholds

There are several depth thresholds that can be defined in the pipeline, I will list them below from the most strict to the least strict.

- consensus_panel_min_depth  = 500

For a given genomic position to be included in the so called "consensus panel" this position needs to have a depth of at least `consensus_panel_min_depth` in at least 80% of the samples. This should always be the highest value among all the depth thresholds and it should be big enough to classify a mutation as somatic vs germline. It should be at least 40.

- sample_panel_min_depth     = 40

This value impacts the creation of sample specific panels that capture which genomics positions have been sequenced to at least this depth in each specific sample. This should be big enough to classify a mutation as somatic vs germline. It should be at least 40.

- mutation_depth_threshold   = 40

This value is used for filtering the mutations by depth. Meaning that if a mutation does not reach this minimum sequencing depth it will not be kept for further analysis. This value should be big enough to be able to classify a mutation as somatic vs germline, and reach a trustworthy computation of the mutation frequency. It should be at least 40.

- use_custom_minimum_depth    = 0

This value is the less stringent depth threshold and is used in the first step of computing the positions that may be part of the so called "panels". This value indicates the minimum average depth at a given position for this position to be kept for the posterior depth analysis and definition on panels. The main use of this value should be to reduce the size of the files that are being processed afterwards. This can be set to 20 or more very safely.

### VAF-distortion filter

- vaf_distortion_threshold   = 3

Mutations whose ratio `VAF_AM / VAF` (all-molecules VAF over duplex VAF) exceeds this threshold are flagged as VAF-distorted during mutation filtering. Lower values are more conservative.

### Using a precomputed depths table

If you already have a precomputed table with per-position depths for your cohort (for example produced by a previous run or an external tool), you can instruct the pipeline to use that table instead of re-computing depths from the BAM files. This can save time and compute resources when depth computation has been performed once and re-used.

Set the following parameters in your `nextflow.config` or pass them on the command line:

```console
params {
  use_custom_depths     = true
  custom_depths_table   = '/path/to/precomputed_depths_table.tsv'
}
```

Notes and requirements:

- `use_custom_depths` (boolean): when `true`, the pipeline will use the file pointed by `custom_depths_table` instead of computing depths from BAMs.
- `custom_depths_table` (string / file path): path to the precomputed depths table. It should be an absolute or relative path accessible from the running environment. The file may be TSV or CSV but should follow the same layout expected by deepCSA (per-position depth across samples). If `use_custom_depths=true` and the file is missing or unreadable the pipeline will fail.
  - If your input.csv file contains only `sample` and `vcf` columns, the columns of the depths table have to be the same ones as the sample names indicated in the sample column of the input.csv file.
  - If your input.csv file contains `sample`, vcf and bam columns, the columns of the depths table have to be the same as the name of the BAM files of each sample in the input.csv file.
- Make sure that you remove the column CONTEXT from the table in case you are starting with the all_samples individual depths table that is outputted by deepCSA. Check out the assets/useful_scripts/downsample_depths.ipynb file for an example on how to prepare the input for this parameter.

### Plotting controls

If you are running with sample groups (see [Feature groups](file_formatting.md#feature-groups)), you can control whether plotting steps generate cohort-only outputs or include all group-level plots:

```console
params {
    plot_only_allsamples = true  // only cohort-level plots
}
```

Set `plot_only_allsamples = false` to generate per-group plots alongside the cohort summaries.

### Positive selection with non-protein-affecting profiles

deepCSA can compute positive selection metrics using a **non-protein-affecting** mutational profile (synonymous + intronic + intergenic mutations). This is useful to compare selection metrics against a background that excludes protein-altering events.

```console
params {
    positive_selection_non_protein_affecting = true
}
```

When enabled, outputs will be labeled with the `.non_prot_aff` suffix in the corresponding selection directories (e.g., `omega/`, `omegagloballoc/`).

## Custom mutation calls -- option 1 (building input VCFs and providing them via normal input)

If you want to run deepCSA with your own mutation calls, this is also possible. Reasons behind this would be:

- the variant calling was not done using deepUMIcaller.
- you came up with a set of mutations that you trust and want to force them as the ones to be used for the analysis.

### Step 1: Generate properly formatted input VCFs

For this, you will need to generate a VCF file per sample with the same format as that expected by deepCSA using the following script that you can find in the deepCSA repository in the following relative path:

`assets/useful_scripts/deepcsa_maf2samplevcfs.py`

The script itself contains this brief explanation on the usage and required parameters:

```{console}
#######
# This script converts a mutations file (TSV format) to one or multiple VCF-formatted files.
#######

#######
# Usage:
#######
## If your sample names are NOT in a column called SAMPLE_ID,
## you can use the --sample-name-column option to specify it.

# if the maf is from deepCSA, use this one
# usage: python deepcsa_maf2samplevcfs.py --mutations-file all_samples.somatic.mutations.tsv --output-dir ./test/ --maf-from-deepcsa

# if the maf file is not from deepCSA, use this one
# usage: python deepcsa_maf2samplevcfs.py --mutations-file all_samples.somatic.mutations.tsv --output-dir ./test/



#######
# Mandatory columns in input mutations:
#######

# if the maf is from deepCSA, it must contain the following columns, as they were originally generated
# ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

# if the maf file is not from deepCSA, then it MUST contain the following columns
# ['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'ALT_DEPTH']
# where:
#     DEPTH indicates the total number of duplex reads sequenced at the position where the mutation occurs
#     ALT_DEPTH indicates the total number of duplex reads supporting the variant at the same position
```

### Step 2: Prepare input.csv file

Make sure to prepare the input.csv file with matching the correct VCF-BAM files for each sample.

If you want to run deepCSA as a basic user and ensure that mutations are properly filtered stop here.

### (optional; advanced users) Step 3: Force no filtering of variants

In case you are following these steps to run deepCSA with a set of mutations that you already filtered and trust there is one last thing that you should do.

When running the pipeline you should set the following parameters:

```console
params {
    no_filter               = true
    filter_criteria         = []
    filter_criteria_somatic = []
}
```

## MAF file as input (alternative input mode)

In addition to the standard per-sample VCF input described above, deepCSA supports providing all mutations from an entire cohort in a single **MAF file** via the `--input_maf` parameter.

### When to use this mode

- You have a cohort-level MAF/TSV file with mutations already called for multiple samples and want to run the downstream deepCSA analysis (mutational profiles, signatures, positive selection, etc.) without preparing one VCF per sample manually.
  - Ideally this file should be generated with the same format as it is generated within deepCSA.
- You already have a precomputed depths table for your cohort (e.g. produced by a previous deepCSA run) or have the information to generate it.

### Requirements

`--input_maf` **must** be used together with `--use_custom_depths true` and a valid `--custom_depths_table`, because BAM-based depth computation is not performed in this mode.

The standard `--input` (samplesheet CSV) is still required to supply sample metadata used by other pipeline steps.

### MAF file format

The MAF file must be tab-separated with a `.maf` extension. The mandatory columns depend on the origin of the file:

| Origin | Mandatory columns | Optional columns |
|---|---|---|
| deepCSA output | `CHROM`, `POS`, `REF`, `ALT`, `FILTER`, `INFO`, `FORMAT`, `SAMPLE`, `SAMPLE_ID` | - |
| External / non-deepCSA | `CHROM`, `POS`, `REF`, `ALT`, `DEPTH`, `ALT_DEPTH`, `SAMPLE_ID` | `DEPTH_AM`, `ALT_DEPTH_AM` |

The `SAMPLE_ID` column identifies each individual sample; internally, the pipeline will generate one VCF file per unique value in that column. The names of the samples in this table should be the same as those indicated in the `input.csv` file.

### Running the pipeline with `--input_maf`

```bash
nextflow run bbglab/deepCSA \
    --input        samplesheet.csv \
    --outdir       results/ \
    --input_maf    cohort_mutations.maf \
    --use_custom_depths    true \
    --custom_depths_table  precomputed_depths.tsv \
    -profile <DESIRED_PROFILE>
```

```console
params {
    input                = "samplesheet.csv"
    outdir               = "results/"
    input_maf            = "cohort_mutations.maf"
    use_custom_depths    = true
    custom_depths_table  = "precomputed_depths.tsv"
}
```

### What happens under the hood

1. The MAF file is passed to `deepcsa_maf2samplevcfs.py`, which converts the cohort-level file into individual per-sample VCFs compatible with the deepCSA input format.
2. The per-sample VCFs are published under `<outdir>/processing_files/input_vcfs/`.
3. The rest of the pipeline proceeds identically to a standard run using per-sample VCFs.

> **Note:** If `--input_maf` is provided without `--use_custom_depths true`, the pipeline will stop immediately with an error message rather than silently ignoring the MAF file.
