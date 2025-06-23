# File formats of inputs

Many of the mandatory or optional files that need to be provided as input have specific requirements.

For most of them there is an example version in the `assets/example_inputs` directory, but this document is intended to explain the details on each of the required files.


## Required input files

### BAM file

### VCF file



## Optional input files

Check the `assets/example_inputs` directory for full information on which are the expected file formats for the different required and mandatory inputs.


### Regions BED file

### Feature groups

### Gene groups

### Hotspots file

### cosmic_ref_signatures

### wgs_trinuc_counts

### cadd_scores

### cadd_scores_ind

### dnds_ref_transcripts

### dnds_covariates

### datasets3d

### annotations3d

### omega_hotspots_bedfile

### Custom TSV regions annotation file should contain this

Document this custom_regions has to be a TSV file with the following columns:

- chromosome  start   end gene_name    impactful_mutations [neutral_impact] [new_impact]

- chromosome start and end indicate the region that is being customized

- gene_name           : is the name of the region that is being added, make sure that it does not coincide with the name of any other gene.

- impactful_mutations : is a comma-separated list of SNVs that need to be labelled with the value indicated in new_impact, format: chr5_1294991_C/T, with pyrimidine based definition

- neutral_impact      : (optional, default; synonymous)

- new_impact          : (optional, default: missense) is the impact that the mutations listed in impactful_mutations will receive.
