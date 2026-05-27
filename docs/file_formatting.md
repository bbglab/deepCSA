# File formats of inputs

Many of the mandatory or optional files that need to be provided as input have specific requirements.

For most of them there is an example version in the `assets/example_inputs` directory, but this document is intended to explain the details on each of the required files.

## Required input files

### BAM file

A BAM file should be provided for each sample containing only the duplex reads belonging to that particular sample. Note that this will be used for normalizing the mutation counts, so it is essential that this contains only the duplex reads or those reads that were used for calling mutations.

### VCF file

VCF file needs to be a single sample VCF with the format field following the specific deepCSA requirements.

More details on this can be found in the [Usage docs](usage.md#custom-mutation-calls).

### Domain definition file

This must be a tab separated file with these 4 columns, each of them is self-explanatory, the only comment is that the coordinates for Begin and End must be protein coordinates.

```bash
Ens_Transcr_ID	Begin	End	NAME
ENST00000240079	30	167	CCDC53
ENST00000447540	2482	2520	BRK
ENST00000447540	2555	2598	BRK
ENST00000447540	878	1150	SNF2_rel_dom
ENST00000447540	1183	1296	Helicase_C
ENST00000447540	774	826	Chromo
ENST00000447540	691	750	Chromo
ENST00000256319	8	118	Erg28
ENST00000278618	125	244	ACPS
```

We suggest getting this information from sources such as Pfam, InterPro or (in a near future) our in-house catalog of domains, [bbgdomains](https://github.com/bbglab/bbgdomains) which is still work in progress.

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

See example below.

#### Related parameters

```console
params {
    features_table              = null
    features_table_separator    = 'comma'
    features_unique_identifier  = null
    features_groups_list        = null
}
```

#### Example features file and parameters

```csv
SAMPLE_ID,BLADDER_LOCATION,SEX,SMOKING_STATUS
donor1_BDO,dome,M,never
donor1_BTR,trigone,M,never
donor2_BDO,dome,F,former
donor2_BTR,trigone,F,former
donor3_BDO,dome,F,former
donor4_BDO,dome,M,former
donor4_BTR,trigone,M,former
donor5_BDO,dome,F,never
donor5_BTR,trigone,F,never
```

```console
params {
    features_table              = "features_example.csv"
    features_table_separator    = 'comma'
    features_unique_identifier  = 'SAMPLE_ID'
    features_groups_list        = "[ ["BLADDER_LOCATION"], ["BLADDER_LOCATION", "SEX"], ["BLADDER_LOCATION", "SEX", "SMOKING_STATUS"] ]"
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

The accompanying hotspot_expansion parameter is used for expanding the definition of this region to analyse a given number of aminoacids on each extreme. The goal of this is to be able to capture the signal focused on the area of interest but also including some more neighbouring areas to potenitally increase the stability of the values in case that the targeted area is very small. In summary, except for the cases in this the `subgenic_bedfile` defined regions are very small <10 AA the `hotspot_expansion` parameter should be set to 0.

#### Related parameters

```console
params {
    subgenic_bedfile      = null
    hotspot_expansion     = 0
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

Path to the COSMIC SBS signature reference file used by SigProfilerAssignment. Use the SBS 96 context file for your genome build (e.g., `COSMIC_v3.4_SBS_GRCh38.txt`). The file is a tab-delimited matrix where the first column encodes mutation context and each additional column corresponds to a signature.

### wgs_trinuc_counts

Tab-delimited file with two columns:

```text
CONTEXT	COUNT
ACA	118979126
ACC	67570313
...
```

The file represents the **total number of occurrences of each trinucleotide** in the reference genome. The pipeline provides a default example in `assets/trinucleotide_counts/`.

### cadd_scores

Path to the CADD "All possible SNVs" file (BGZIP-compressed TSV). This file is used for OncodriveFML scoring.

Recommended download: [CADD downloads](https://cadd.gs.washington.edu/download) → "All possible SNVs of GRCh38/hg38".

### cadd_scores_ind

Tabix index (`.tbi`) for the `cadd_scores` file. If you need to generate it:

```bash
bgzip -c whole_genome_SNVs.tsv > whole_genome_SNVs.tsv.gz
tabix -s 1 -b 2 -e 2 whole_genome_SNVs.tsv.gz
```

### dnds_ref_transcripts

https://github.com/bbglab/deepCSA/tree/dev/assets/build_datasets/dndscv

### dnds_covariates

dNdScv covariates file, usually `covariates_hg19_hg38_epigenome_pcawg.rda`. This provides covariate regression terms for mutation rate modeling.

### datasets3d

Directory containing precomputed Oncodrive3D datasets (structure and mutation mapping information). Build using the [Oncodrive3D dataset builder](https://github.com/bbglab/oncodrive3d?tab=readme-ov-file#building-datasets).

### annotations3d

Directory containing Oncodrive3D annotation datasets (protein annotations, stability data, etc.). Use the same build process as `datasets3d` to ensure compatibility.

### gff3_file

Optional local GFF3 file used by the DNA2PROTEINMAPPING step. If not provided, the pipeline downloads the GFF3 from Ensembl. If provided, it must match the Ensembl release, species, and genome build you are using (compressed `.gff3.gz` files are supported).

## Examples

### Blacklist mutations

```
chr1:11107296_C>CA
chr1:11107450_C>A
chr1:11108379_T>A
```

### Gene grouping

```
chr15q  chr15q  IDH2    SIN3A
chr17p  chr17p  MAP2K4  NCOR1   TP53    USP6
```

### Custom annotation

See `assets/example_inputs/custom_regions.example.tsv` for a full example of a custom region annotation file.

### Omega hotspots / subgenic regions

Provide a BED file with 3 or 4 columns (`CHROM`, `START`, `END`, optional `NAME`):

```
chr7    55191765    55191840    EGFR_L858R_region
chr12   25245300    25245380    KRAS_G12_region
```

You can expand these regions with `hotspot_expansion` and optionally generate complements with `subgenic_regions_complement`.
