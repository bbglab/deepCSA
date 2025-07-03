# bbglab/deepCSA: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Table of contents

- [Introduction](#introduction)
- [How to run the pipeline](#how-to-run-the-pipeline)
- [Samplesheet input](#samplesheet-input)
- [Proposed run modes](#proposed-run-modes)
  - [Initial run. Data exploration](#initial-run-data-exploration)
  - [Clonal structure definition](#clonal-structure-definition-complete-run-with-a-focus-on-positive-selection-at-the-cohort-level)
  - [Mutational processes in alternative genomic regions](#mutational-processes-in-alternative-genomic-regions-partial-run-with-a-focus-on-mutational-processessignatures)
  - [Interindividual variability and sample comparison](#interindividual-variability-and-sample-comparison-complete-run-with-downstream-steps-for-computation-of-linear-regressions-to-compare-different-samplesgroups-based-on-clinical-variables-or-sample-metadata)
- [Definition of structural parameters](#definition-of-structural-parameters)
- [Additional customizable parameters](#additional-customizable-parameters)
- [Custom mutation calls](#custom-mutation-calls)
- [Running the pipeline](#running-the-pipeline)
- [Core Nextflow arguments](#core-nextflow-arguments)
- [Custom configuration](#custom-configuration)
- [Running in the background](#running-in-the-background)
- [Nextflow memory requirements](#nextflow-memory-requirements)

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## How to run the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run bbglab/deepCSA --outdir <OUTDIR> -profile <DESIRED PROFILE> --input samplesheet.csv
```

For more information on how to run Nextflow pipelines check a more detailed explanation [below](#running-the-pipeline) in this same document or check the [Nextflow](https://www.nextflow.io/docs/latest/index.html) or [nf-core](https://nf-co.re) community documentations.

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

## Proposed run modes

<!-- TODO decide if this could be added as different profiles/modes for parameter definition -->

- Initial run. Data exploration
- Clonal structure definition. Complete run with a focus on positive selection at the cohort-level.
- Mutational processes in alternative genomic regions. Partial run with a focus on mutational processes/signatures.
- Interindividual variability and sample comparison. Complete run with downstream steps for computation of linear regressions to compare different samples/groups based on clinical variables or sample metadata.

### Initial run. Data exploration

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

- All the previously described outputs plus...
- Mutation density
- Positive selection per gene
  multiple positive selection metrics
  - Per gene, all samples together
  - Per gene, per group of samples
  - Per gene, per sample

```console
params {
    mutationrate                = true

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

    omega_withingene            = true
    omega_autodomains           = true
    omega_autoexons             = true
    
    mutated_cells_vaf           = true
    mutepi_genes_to_recode      = null

    indels                      = true

    signatures                  = true
}
```

### Mutational processes in alternative genomic regions. Partial run with a focus on mutational processes/signatures

- Same as initial run (even it can be ignored)
- Mutational profile and mutational signatures based on:
  - All genomic regions
  - Only exonic regions
  - Only non-protein affecting regions (synonymous mutations and intronic, intergenic)
  - Intronic and intergenic regions

```console
params {
    mutationrate                = true

    profileall                  = true
    profilenonprot              = true
    profileexons                = true
    profileintrons              = true

    signatures                  = true
}
```

### Interindividual variability and sample comparison. Complete run with downstream steps for computation of linear regressions to compare different samples/groups based on clinical variables or sample metadata

- Same as complete clonal structure definition +
- Computation of univariate and multivariate linear regressions between clonal structure metrics and clonal selection

```console
params {
    mutationrate                = true

    profileall                  = true

    omega                       = true
    omega_multi                 = true
    omega_globalloc             = true
    omega_mutabilities          = true
    site_comparison_grouping    = 'all'
    omega_plot                  = true

    omega_withingene            = true
    omega_autodomains           = true
    omega_autoexons             = true
    
    regressions                 = true
    // additional regression parameters, see nextflow_schema.json for more info
      ...
}
```

## Definition of structural parameters

- Container pulling (either prior to running the pipeline or directly as the pipeline runs)
- Generation of Oncodrive3D datasets (see: [Oncodrive3D repo datasets building process](https://github.com/bbglab/oncodrive3d?tab=readme-ov-file#building-datasets))

- Download of additional specific datasets
  - Ensembl VEP (see: [Ensembl VEP docs](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache)). Modify accordingly your `nextflow.config` vep parameters, `vep_cache`, `vep_cache_version`, etc.
  <!-- TODO we should revise if we can provide more specific information on how to download the cache -->
  - CADD scores (see: [CADD downloads page](cadd.gs.washington.edu/download) "All possible SNVs of GRCh38/hg38" file)
  - COSMIC signatures (i.e. [COSMIC signatures downloads page](https://cancer.sanger.ac.uk/signatures/downloads/) (select context size = 96 and your desired species of interest))
  <!-- * dNdScv datasets (see: ) -->
  
### Mandatory parameter configuration

See [File formatting docs](file_formatting.md) for more details on the structure of files that can be provided to deepCSA.

```console
params {

    fasta                       = null

    cosmic_ref_signatures      = "COSMIC_v3.4_SBS_GRCh38.txt"
    wgs_trinuc_counts          = "assets/trinucleotide_counts/trinuc_counts.homo_sapiens.tsv"

    // oncodrivefml (only for human; could be adapted to others)
    cadd_scores                = "CADD/v1.7/hg38/whole_genome_SNVs.tsv.gz"
    cadd_scores_ind            = "CADD/v1.7/hg38/whole_genome_SNVs.tsv.gz.tbi"

    // dnds
    dnds_ref_transcripts       = "RefCDS_human_latest_intogen.rda"
    dnds_covariates            = "covariates_hg19_hg38_epigenome_pcawg.rda"

    // oncodrive3d + fancy plots
    datasets3d                 = "oncodrive3d/datasets"
    annotations3d              = "oncodrive3d/annotations"


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
    omega_hotspots_bedfile      = null

    // define a file of mutations that should not be trusted
    //  and you want to remove from all the analysis
    blacklist_mutations        = null
}
```

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

## Custom mutation calls

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

## Additional Nextflow documentation

### Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run bbglab/deepCSA --input ./samplesheet.csv --outdir ./results --fasta genome.fa -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run bbglab/deepCSA -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
fasta: 'genome.fa'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull bbglab/deepCSA
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [bbglab/deepCSA releases page](https://github.com/bbglab/deepCSA/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

### Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

### Custom configuration

#### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here (example nf-core/rnaseq)](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

#### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

#### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
