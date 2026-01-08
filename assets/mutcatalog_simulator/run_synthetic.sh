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

for file in ${OUTPUT_DIR}/maf/*.tsv; do
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