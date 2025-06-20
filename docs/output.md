# bbglab/deepCSA: Output

## Introduction

This document describes the output produced by the pipeline.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Directory Structure](#directory-structure)
- [Input and configuration](#input-and-configuration)
- [Depth analysis](#depth-analysis)
- [Mutation preprocessing](#mutation-preprocessing)
- [Basic analysis](#basic-analysis)
- [Intermediate outputs](#intermediate-outputs)
- [Positive selection](#positive-selection)
- [Site selection metrics](#site-selection-metrics)
- [Additional clonal structure metrics](#additional-clonal-structure-metrics)
- [Mutational signatures](#mutational-signatures)
- [Plotting funcionalities](#plotting-funcionalities)
- [Additional outputs](#additional-outputs)


## Directory Structure

The directory structure listed below will be created in the results directory after the pipeline has finished.
The structure captures the maximum diversity of created outputs, but when only certain run options are turned on, not all directories will be generated.
All paths are relative to the top-level results directory.

```{console}
{outdir}
├──absolutemutabilities
├──absolutemutabilitiesgloballoc
├──annotatedepths
├──clean_germline_somatic
├──clean_somatic
├──computematrix
├──computeprofile
├──createpanels
│   ├── consensus
│   │   └── <region_type>.consensus.bed
│   │   └── <region_type>.consensus.tsv
│   ├── captured
│   │   └── <region_type>.captured.bed
│   │   └── <region_type>.captured.tsv
│   └── sample
│       └── <region_type>.<sample>.bed
│       └── <region_type>.<sample>.tsv
├──customannotation
├──customprocessing
├──customprocessingrich
├──depthssummary
├──dna2proteinmapping
├──domainannotation
├──expandregions
├──filterexons
├──germline_somatic
├──groupgenes
├──indels
├──matrixconcatwgs
├──multiqc
├──mutability
├──mutatedcellsfromvafam
├──mutatedgenomesfromvafam
├──mutrate
├──muts2sigs
├──omega
│   ├── preprocessing
│   │   └── syn_muts.<sample>
│   │   └── mutabilities.<sample>
│   └── output_mle.<sample>.tsv
├──omegagloballoc
│   ├── preprocessing
│   │   └── syn_muts.<sample>
│   │   └── mutabilities.<sample>
│   └── output_mle.<sample>.tsv
├──oncodrive3d
│   ├── run
│       └── <sample>
│   └── plot
│       └── <sample>
├──oncodrivefmlsnvs
├──pipeline_info
├──plotmaf
├──plotneedles
│   └── <sample>
│       └── <gene-sample needles>
├──plotselection
├──plotsomaticmaf
├──postprocessveppanel
├──signatures_hdp
│   └── output.<region_type>
│       └── <sample>
├──sigprobs
├──sigprofilerassignment
│   └── output.<region_type>
│       └── <sample>
├──sitecomparison
├──sitecomparisongloballoc
├──sitecomparisongloballocmulti
├──sitecomparisonmulti
├──sitesfrompositions
├──sumannotation
├──synmutrate
├──synmutreadsrate
└──table2group
work/
.nextflow.log
```

## Input and configuration

See Usage docs for extensive explanation on required inputs and format. Including documentation on parameters to run on for 4 different suggested running modes.

## Depth analysis

### Key role

- Computation of depth per sample for each specific position
  Most analysis may be influenced by sequencing depth, it is essential to correct for these values.

- Definition of regions to analyze
  Only genomic areas that have been properly covered across samples will be used for the analysis.

### Outputs

- sitesfrompositions
- postprocessveppanel
- createpanels
- annotatedepths
- depthssummary

Optional:

- dna2proteinmapping
- domainannotation
- customprocessing
- customprocessingrich

## Mutation preprocessing

### Key role

- VCF annotation: Annotate mutations with Ensembl VEP.
- VCF to MAF conversion: Convert VCFs to MAF, define VAF, and merge with annotation.
- Custom region annotation: Allow user to define different consequence types for specific regions.
- Hotspot annotation: Add known hotspots to mutation annotation.
- Filtering:
  - Filter mutations at the sample level (e.g., VAF distortion).
  - Filter at the cohort level (e.g., other_sample_SNP, repetitive_variant, not_covered, not_in_exons).
- Blacklist mutations if activated (see assets for example).
- Downsample mutations if activated.


### Outputs

- sumannotation
- customannotation
- germline_somatic
- clean_somatic
- clean_germline_somatic

## Basic analysis

### Key role

- Mutation density computation
  Correct the number of mutations observed by the number of sequenced nucleotides.

- Mutational profile computation
  Capture the mutation probability of each trinucleotide. Represent it in three different normalization conditions.



### Outputs

- computematrix
- computeprofile
- mutrate

## Intermediate outputs

### Key role

- Matrix concatenation
  Combine WGS-renomralized matrices for mutational signature analysis.

- Mutability calculation
  Compute relative mutabilities using depths and mutational profile.

- Choose synonymous mutation rates for downstream analysis.

### Outputs

- matrixconcatwgs
- mutability
- synmutrate
- synmutreadsrate

## Positive selection

### Key role

- Compute multiple positive selection metrics
  This is done at the cohort-level, but also for each sample or group of samples.

- OncodriveFML: Detects functional impact bias in observed mutations.

- Oncodrive3D: Identifies 3D protein regions with mutation clustering, using relative mutabilities and raw VEP annotation.

- Omega: dN/dS-based, quantifies selection pressure in defined regions (genes, exons, domains, hotspots, etc.).

- Indels: Analysis of indel selection.

### Outputs

- omega
- omegagloballoc
- oncodrive3d
- oncodrivefmlsnvs
- indels

## Site selection metrics

### Key role

- Compute absolute mutabilities for each position.

- Compare the observed number of mutations per site to the expected number of mutations and estimate a site selection value.

### Outputs

- absolutemutabilities
- absolutemutabilitiesgloballoc

- sitecomparison
- sitecomparisongloballoc
- sitecomparisongloballocmulti
- sitecomparisonmulti

## Additional clonal structure metrics

### Key role

- VAF-based definition of the number of mutated genomes.

### Outputs

- mutatedcellsfromvafam
- mutatedgenomesfromvafam


## Mutational signatures

### Key role

- Signature assignment: Use SigProfilerAssignment with optional custom signatures.
- HDP: Hierarchical Dirichlet Process for signature extraction.
- (Pending) Signature extraction: SigProfilerExtractor support.

### Outputs

- signatures_hdp
- sigprofilerassignment
- sigprobs
- muts2sigs

## Plotting funcionalities

### Key role

- Plotting basic statistics of numbers and distribution of mutations in genes.

- Optionally think on adding more plots.

### Outputs

- plotmaf
- plotneedles
- plotselection
- plotsomaticmaf



## Additional outputs

### Key role

- Definition of groups, expanded regions and other metrics related with the full pipeline execution.

### Outputs

- table2group
- groupgenes
- expandregions
- filterexons
- multiqc
- pipeline_info
