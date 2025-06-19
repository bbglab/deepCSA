# bbglab/deepCSA: Output

## Introduction

This document describes the output produced by the pipeline.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Directory Structure](#directory-structure)
- [Input and configuration](#)
- [Depth analysis](#)
- [Mutation preprocessing](#)
- [Basic analysis](#)
- [Positive selection](#)
- [Mutational signatures](#)
- [Plots along the way](#)

## Directory Structure

The directory structure listed below will be created in the results directory after the pipeline has finished.
The structure captures the maximum diversity of created outputs, but when only certain run options are turned on, not all directories will be generated.
All paths are relative to the top-level results directory.

```
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

- Comp
  Most

### Outputs

- sumannotation
- customannotation
- germline_somatic
- clean_somatic
- clean_germline_somatic

## Basic analysis

### Key role

- Comp
  Most

### Outputs

- computematrix
- computeprofile
- mutrate

## Intermediate outputs

### Key role

- Comp
  Most

### Outputs

- matrixconcatwgs
- mutability
- synmutrate
- synmutreadsrate

## Positive selection

### Key role

- Comp
  Most

### Outputs

- omega
- omegagloballoc
- oncodrive3d
- oncodrivefmlsnvs
- indels

## Site selection metrics

### Key role

- Comp
  Most

### Outputs

- absolutemutabilities
- absolutemutabilitiesgloballoc

- sitecomparison
- sitecomparisongloballoc
- sitecomparisongloballocmulti
- sitecomparisonmulti

## Additional clonal structure metrics

### Key role

- Comp
  Most

### Outputs

- mutatedcellsfromvafam
- mutatedgenomesfromvafam


## Mutational signatures

### Key role

- Comp
  Most

### Outputs

- signatures_hdp
- sigprofilerassignment
- sigprobs
- muts2sigs

## Plotting funcionalities

### Key role

- Comp
  Most

### Outputs

- plotmaf
- plotneedles
- plotselection
- plotsomaticmaf



## Additional cohort description outputs

### Key role

- Comp
  Most

### Outputs

- table2group
- groupgenes








Define regions to analyze
Compute depth to define regions to analyze
Predefine a
MAX BED file
Predefine a minimum average depth across samples

Annotate functional impact of all possible positions
Most deleterious or canonical only
Customize the impact of certain positions
Via: custom_annotation
custom_annotation_tsv

Create all panel, and the different impact-based regions.
This contain all the positions sequenced at least once that comply with the impact criteria.

Define well-covered regions (create consensus panels)
Define minimum depth to be included
Define minimum number of samples that must comply with it
Report which samples have included positions below minimum depth.
? Report which positions were excluded from the panel.

Create sample panels
These are based on a threshold that applies only to the depth of the sample.
They could potentially be used for signatures or something like this but the reality is that we don't use them.

Plot depths of the genes/samples
Use the panels defined above to decide which regions to plot.
Report depth per sample, per gene sample, also relative values per exons and more.
Generate a pdf with tons of pages, maybe generating two pdfs would be better. One short and summarized and another big one with all the details

Generate depth file per sample.

Depths can be downsampled uniformly. The downsampled depth is deterministic, repeating the downsampling would give the same depth values.

Subset those depth files to each possible panel regions for downstream analysis.

Preprocess mutations vcfs and filter mutations.
Annotate mutations with Ensembl VEP from VCF
Process the concatenated annotations with python
Both only canonical + most_deleterious
SHOULD ADD the potential SNP filter here
Customize regions if needed
If the user wants to define a different consequence type to a specific region of the genome, this can be done by providing a file with the corresponding details.
Add information on which mutations are known hotspots
This can be done with the information extracted from IntoGen or from any other source and has a specific input format.
Convert the VCF to MAF and define the VAF and all the variables inferred from the VCF also merge with the annotation.
Working with VarDict output after our custom pileup recounting
SHOULD BE ADAPTED to work with any VCF, maybe this is useful:
<https://gist.github.com/FerriolCalvet/a5cb90064272197948d63cbe593433a7>

Add filters at the level of sample
(currently adding VAF distortion here, but it could be moved to deepUMIcaller)

Compute filters at the level of cohort:
other_sample_SNP
We use same germline threshold as in somatic mutations filter to separate between germline and somatic. If any somatic mutation in a sample is detected as germline in another, it is labelled with this flag.
cohort_n_rich
cohort_n_rich_threshold
cohort_n_rich
cohort_n_rich_uni
Repetitive_variant
If variant is present in more than N samples, it is flagged as repetitive_variant.
Define not_covered flag:
This mutation is outside the consensus panel so it is not well-covered across all the samples
Define not_in_exons flag:
This mutation is not in exonic regions
Blacklist mutations (if activated)
Add example in assets.
Downsample mutations (if activated)

Summary plot of mutations:
Currently using all mutations
Would be nice to add distortion plots here
Would be nice with somatic mutations only, to report what is continuing the process.
Generating a needle plot per sample & gene.

DNA2PROTEIN mapping
Revise the code for this and how useful it is for generating the plots.
We can probably remove it to start using the panel rich, but check how easy it is to make the plots then.

If the user provided some clinical variables to generate groups, the samples of each group, are pooled together as if they were a single sample with all those mutations and all the depth.

DATA PREPROCESSING FINALIZES HERE

DATA ANALYSIS

Mutational profile
Chose panel where to focus the analysis:
All
Exons only
Introns and intergenic only (or non protein affecting only)
Denominator (COMPUTETRINUC)
Count every position in the panel once
Count every position in the panel as many times as it was sequenced
Use the sequencing depth of each sample. (default)
Numerator (COMPUTEMATRIX)
Unique number of mutations (count each mutation only once per sample)
(default)
Mutated reads aka expanded = True
(count each mutation as many times as reads we find in each sample)
Outputs:
Plot of the 96 channels mutational profile.
Mutational profile TSV file with the proportions.
Matrix WGS
Matrix of theoretical counts observed if the sequenced region was the WGS.
Using the mutational profile and the trinucleotide counts of the WGS we get an estimation on which would be the distribution of the number of mutations observed in each sample.
Potential outputs:
We could use the profile of all_samples or a specific group to use as a pseudocount for the individual profiles, maybe not here but in the posterior analysis or something like this. When the counts are very small this would help to smooth the mutation probability between trinucleotides, and when the counts are bigger this should have minimal impact, we could do it in a way that the background accounts for 20 mutations in total and then add the real counts on top.

Mutabilities (relative mutabilities)
Take the depths, specific regions, the mutational profile and the observed mutations.
With the depths - mutational profile we obtain a mutability that is a completely relative metric to quantify the differences in mutational probability between two different positions in the sample.
The observed mutations are not uniformly distributed between genes so the information on genes can also be used to weight the “distribution” of mutations between the different genes, this can be used by ofml when not randomizing mutations within a gene, but within a group of genes.
(ofml by groups not fully tested yet, but yes implemented, and it is probably working)
The impact of this weighted distribution of mutabilities between genes has no impact on the usage of mutabilities at the level of each gene independently.

Compute mutation rate density
Numerator
Number of mutations -> mutation rate
Number of mutated reads -> mutreads rate
Denominator
Count every position with the full depth at that position.
Adjusted mutation rate count every site with ⅓ of the depth at a given position, if a site has different consequence types, it is differently balanced between PA and NPA.
Computed for all impacts, PA, NPA and synonymous and for each of them we computed with and without indels and also with only indels.

Mutated genomes
VAF
This takes into account the value of omega to estimate the number of mutations in excess for each gene sample. With this number of mutations in excess, E, the VAF_AM (if available) of each mutation is used for ranking the mutations in descending order within that gene sample and the top E mutations are defined as drivers. This is done for each consequence type, and then we take the VAFs of those mutations and apply the principle of inclusion exclusion to obtain a value of mutated genomes for that consequence type, gene and sample, taking into account the random probability of two mutations occurring in the same genome.
This step requires a SAMPLE_ID and a SEX column defined in the metadata file provided to the pipeline. SEX values need to be encoded as M and F.
Reads
No need to describe it, this has been outdated. I will describe it as soon as possible once everything else has been properly documented.
Expected
MutRisk
Signatures
SigProfilerAssignment
You can provide a custom set of signatures to do the decomposition.
SigProfilerExtractor
Pending implementation, waiting for a SigProfilerExtractor working container
HDP
Implemented, pending group definition to get per group plots.
Positive selection.
We run several methods to quantify clonal selection in the samples:

Oncodrivefml
It is based on measuring the functional impact bias. Using a predefined set of impact scores for all possible mutations, are the mutations we see enriched in high impact mutations? See paper for more details.

Oncodrive3d
It is based on identifying protein regions that are close in the 3D space where there are more mutations than expected by the neutral mutagenesis. See paper for more details.
This process consists in running Oncodrive3D code using the relative mutabilities which assign a given mutation probability to each of the possible mutations in a gene.
In general the option of passing the raw VEP annotation to O3D should be activated since this was O3D chooses the annotation that corresponds to the protein structure used for computing the volumes in the 3D space and all this.
As part of the output genes are scored but also positions within those genes.
In addition to the table output the user can also ask for automated plots summarizing the results in the cohort, and also a plot per gene that shows how the mutability is distributed along the gene together with the missense mutations the clustering score and more information about the protein sequence.
Optionally one can ask for a plot of the 3d structure with the residues that are in clusters colored.

Oncodriveclustl
Postponed

Omega
(logistics: requires profileall, mutrate, omega)
Anything can be defined as a region as long as it has a specific genomic coordinates that correspond to it. It can be a single continuous element (a single exon, or a domain that is continuous in the genomic space) or a discontinuous region, such as a gene composed of multiple exons or multi-exon or discontinuous domains.
There are different ways of defining the genes where to run the analysis
ExpandREGIONS
Options to handle this:
Omega_withingene
Omega_hotspots_bedfile
hotspot_expansion
omega_autodomains
Omega_autoexons
If omega_withingene yes, but no other information on how to do it, then error.
If no mutation of the corresponding impact, then no omega is computed.
Loc
If no synonymous mutation, then omega is not computed
Global loc
The number of expected synonymous mutations is derived from a file of mutation rate per MB per gene. This must be provided as required by omega and within deepCSA this corresponds to the synonymous mutation rate when taking all the samples of the cohort together.
