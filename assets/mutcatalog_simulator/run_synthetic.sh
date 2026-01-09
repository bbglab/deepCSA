#!/bin/bash

MUTCATALOG_SIMULATOR_HOME="/home/fcalvet/projects/deepCSA/assets/mutcatalog_simulator"

DEEPCSA_RUN_DIR="/data/bbg/nobackup/bladder_ts/results/2025-12-19_deepCSA_updated_run"
RUN_CONFIG="${MUTCATALOG_SIMULATOR_HOME}/test_config.json"
OUTPUT_DIR="${DEEPCSA_RUN_DIR}/fake_mutations_poisson"

mkdir -p ${OUTPUT_DIR}/maf
mkdir -p ${OUTPUT_DIR}/vcf

# generate MAFs

python ${MUTCATALOG_SIMULATOR_HOME}/synthetic_maf.py \
    --deepcsa_run_dir  "${DEEPCSA_RUN_DIR}" \
    --run_config "${RUN_CONFIG}" \
    --output_dir "${OUTPUT_DIR}/maf"

# generate VCFs

CORES=$(($(nproc) - 1))

for file in $(ls ${OUTPUT_DIR}/maf/*.tsv | grep -v '_depths.tsv'); do
    python ${MUTCATALOG_SIMULATOR_HOME}/deepcsa_maf2samplevcfs.py \
        --mutations-file "$file" \
        --sample-name-column "SAMPLE_ID" \
        --output-dir "${OUTPUT_DIR}/vcf" &

    # This block limits the number of background jobs
    if [[ $(jobs -r | wc -l) -ge $CORES ]]; then
        wait -n  # Wait for at least one job to finish
    fi

done

wait

# generate the deepCSA input files
ls ${OUTPUT_DIR}/vcf/*_??.vcf | xargs -n 1 basename | rev | cut -d '.' -f2- | rev > ${OUTPUT_DIR}/vcf/sample_names.txt

echo "sample,vcf,bam" > ${OUTPUT_DIR}/vcf/deepCSA_input.csv
ls ${OUTPUT_DIR}/vcf/*_??.vcf | xargs -n 1 realpath > ${OUTPUT_DIR}/vcf/vcf_paths.txt
paste -d , ${OUTPUT_DIR}/vcf/sample_names.txt ${OUTPUT_DIR}/vcf/vcf_paths.txt >> ${OUTPUT_DIR}/vcf/deepCSA_input.csv

echo "SAMPLE_ID,original_sample,omega,depth_correction,replicate" > ${OUTPUT_DIR}/vcf/deepCSA_input_metadata.csv
paste -d , <(cat ${OUTPUT_DIR}/vcf/sample_names.txt) <(cat ${OUTPUT_DIR}/vcf/sample_names.txt | rev | sed 's/_/,/; s/_/,/; s/_/,/' | rev) >> ${OUTPUT_DIR}/vcf/deepCSA_input_metadata.csv

rm ${OUTPUT_DIR}/vcf/sample_names.txt
rm ${OUTPUT_DIR}/vcf/vcf_paths.txt