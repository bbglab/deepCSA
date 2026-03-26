#!/bin/bash

# from an existing valid mutcatalog_simulator output directory, generate the deepCSA input files
# it generates required deepCSA input files

DEEPCSA_RUN_DIR="/data/bbg/nobackup/prominent/SantPauCH/deepCSA/sp_vhio/custom/2026-03-23_Sp_CHa_H_T0s_only_custom"
OUTPUT_DIR="${DEEPCSA_RUN_DIR}/fake_mutations_poisson_blood"

ls ${OUTPUT_DIR}/vcf/*.vcf | xargs -n 1 basename | rev | cut -d '.' -f2- | rev > ${OUTPUT_DIR}/vcf/sample_names.txt

echo "sample,vcf,bam" > ${OUTPUT_DIR}/vcf/deepCSA_input.csv
ls ${OUTPUT_DIR}/vcf/*.vcf | xargs -n 1 realpath > ${OUTPUT_DIR}/vcf/vcf_paths.txt
paste -d , ${OUTPUT_DIR}/vcf/sample_names.txt ${OUTPUT_DIR}/vcf/vcf_paths.txt >> ${OUTPUT_DIR}/vcf/deepCSA_input.csv

echo "SAMPLE_ID,original_sample,omega,depth_correction,impact,replicate" > ${OUTPUT_DIR}/vcf/deepCSA_input_metadata.csv
paste -d , <(cat ${OUTPUT_DIR}/vcf/sample_names.txt) <(cat ${OUTPUT_DIR}/vcf/sample_names.txt | rev | sed 's/_/,/; s/_/,/; s/_/,/' | rev) >> ${OUTPUT_DIR}/vcf/deepCSA_input_metadata.csv

rm ${OUTPUT_DIR}/vcf/sample_names.txt
rm ${OUTPUT_DIR}/vcf/vcf_paths.txt