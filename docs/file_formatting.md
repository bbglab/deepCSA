# File formats of inputs

Many of the mandatory or optional files that need to be provided as input have specific requirements.

For most of them there is an example version in the `assets/example_inputs` directory, but this document is intended to explain the details on each of the required files.

## Required input files

### BAM file

A BAM file should be provided for each sample containing only the duplex reads belonging to that particular sample. Note that this will be used for normalizing the mutation counts, so it is essential that this contains only the duplex reads or those reads that were used for calling mutations.

### VCF file

VCF file needs to be a single sample VCF with the format field following the specific deepCSA requirements.

More details on this can be found in the [Usage docs](usage.md#custom-mutation-calls).

## Optional input files

Check the `assets/example_inputs` directory for full information on which are the expected file formats for the different required and mandatory inputs.

### Regions BED file

You can provide a BED file that restricts the genomic regions that will be used for the analysis.
Be careful when providing this file since a too restrictive file focussed only in the exonic regions will result in excluding some potentially useful data from the analysis.

#### Related parameters

```console
params {
    use_custom_bedfile          = false
    custom_bedfile              = null
}
```

### Feature groups

By default all the metrics are computed for each individual sample and for the entire cohort of samples, but additionally this can be computed for groups of samples that need to be indicated using the parameters listed below.

Long story short:

* a comma or tab separated table containing at least a row with information on the sample names formatted equally as those provided as part of the input.csv
* indicate whether the previous file is comma or tab separated
* a dictionary to define which column contains the name of the samples and then combinations of columns that will be used for the groupings. Note that the columns indicated for groupings will be treated as categorical variables even if they contains numerical values. Note also that there should not be missing values in any of the columns used for the groupings, instead replace them (i.e. replace by Unknown)

#### Related parameters

```console
params {
    features_table              = null
    features_table_separator    = 'comma'
    features_table_dict         = ''
}
```

### Gene groups

We also allow the computation of positive selection for groups of genes, particularly for omega, but this might be expanded to other methods. For this to happen, the `custom_groups` parameter needs to be set to `true` and then the corresponding groups file and its separator should be provided. Find an example below.

This is an example on how this file could look like.

```bash
chr15q  chr15q  IDH2    SIN3A
chr17p  chr17p  MAP2K4  NCOR1   TP53    USP6
chr17q  chr17q  AXIN2   BRCA1   CDK12   ERBB2   NF1     PPM1D   RNF43   SOX9    SPOP    SRSF2   STAT3   VEZF1
chr1p   chr1p   ARID1A  CDKN2C  MTOR    NOTCH2  NRAS    RPL22   RPL5    SPEN
chr19p  chr19p  KEAP1   NOTCH3  SMARCA4 STK11
```

#### Related parameters

```console
params {
    custom_groups               = false
    custom_groups_file          = null
    custom_groups_separator     = 'tab'
}
```

### Hotspots file

This file will be used for the annotation of some specific mutations as known hotspots.

We generally use a file with a list of hotspots retrieved from highly frequent mutations in cancer, but this could be changed by any other file with the desired information tailored to each specific use case.

You should provide a file with the following structure (MUTTYPE should be pyrimidine centric, or '-' in case of non SNV variants).

```bash
CHROM   POS     MUTTYPE HOTSPOT_NAME
chr12   25245351        C>A     KRAS_G12C_Missense
chr12   25245350        C>A     KRAS_G12V_Missense
chr7    55191822        T>G     EGFR_L858R_Missense
chr12   25245350        C>T     KRAS_G12D_Missense
```

#### Related parameters

```console
params {
    hotspots_annotation         = false
    hotspots_definition_file    = ''
}
```

### Omega regions BED file

We can also customize omega by providing a file with a specific definition of genomic regions in which omega will be computed independently. For example this could be useful for a custom definition of protein domains, or for a definition of other genomic elements within genes in which the user wants to have an independent measurement of omega.

The accompaining hotspot_expansion parameter is used for expanding the definition of this region to analyse a given number of aminoacids on each extreme. The goal of this is to be able to capture the signal focused on the area of interest but also including some more neighbouring areas to potenitally increase the stability of the values in case that the targeted area is very small. In summary, except for the cases in this the `omega_hotspots_bedfile` defined regions are very small <10 AA the `hotspot_expansion` parameter should be set to 0.

#### Related parameters

```console
params {
    omega_hotspots_bedfile      = null
    hotspot_expansion           = 0
}
```

### Blacklist mutations

The user can provide a file with mutations that are very likely artifacts and should be removed from any analysis. This could be a list compiled across mutiple projects and coming from past knowledge. The format of this list is a single MUT_ID per line, with MUT_ID being generated in the same format that we use within deepCSA.

```bash
chr1:11107296_C>CA
chr1:11107450_C>A
chr1:11108379_T>A
```

#### Related parameters

```console
params {
    blacklist_mutations         = null
}
```

### Customize annotations at specific positions

Document this custom_regions has to be a TSV file with the following columns:

* chromosome  start   end gene_name    impactful_mutations [neutral_impact] [new_impact]

* chromosome start and end indicate the region that is being customized

* gene_name           : is the name of the region that is being added, make sure that it does not coincide with the name of any other gene.

* impactful_mutations : is a comma-separated list of SNVs that need to be labelled with the value indicated in new_impact, format: chr5_1294991_C/T, with pyrimidine based definition

* neutral_impact      : (optional, default; synonymous)

* new_impact          : (optional, default: missense) is the impact that the mutations listed in impactful_mutations will receive.

```console
CHROM   START   END     NAME    IMPACTFUL       NEUTRAL IMPACT
chr5    1294942 1295289 TERTpromoter    chr5_1295017_C/A,chr5_1295057_C/T,chr5_1295061_C/G,chr5_1295061_C/T,chr5_1295081_C/T,chr5_1295138_C/A,chr5_1294991_C/T,chr5_1295054_C/T,chr5_1295057_C/A,chr5_1295081_C/G,chr5_1295115_C/T,chr5_1295132_C/T,chr5_1295211_C/T,chr5_1295219_C/T,chr5_1295259_C/T,chr5_1295018_C/T,chr5_1295032_C/G,chr5_1295043_C/A,chr5_1295034_C/A,chr5_1295127_C/T,chr5_1295138_C/T,chr5_1295090_C/T,chr5_1295113_C/A,chr5_1295046_T/G,chr5_1295135_C/T,chr5_1295113_C/T   synonymous      missense
```

#### Related parameters

```console
params {
    customize_annotation        = false
    custom_annotation_tsv       = ''
}
```

### cosmic_ref_signatures

### wgs_trinuc_counts

### cadd_scores

### cadd_scores_ind

### dnds_ref_transcripts

### dnds_covariates

### datasets3d

### annotations3d
