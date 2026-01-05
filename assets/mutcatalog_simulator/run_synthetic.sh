#!/bin/bash

DEEPCSA_RUN_DIR="/data/bbg/nobackup/lung_duplex/analysis/fullcohortnormal/2025-09-19_PEACE"
RUN_CONFIG="./test_config.json"
OUTPUT_DIR="${DEEPCSA_RUN_DIR}/fake_mutations_poisson"

mkdir -p ${OUTPUT_DIR}/maf
mkdir -p ${OUTPUT_DIR}/vcf

# generate MAFs

python synthetic_maf.py \
    --deepcsa_run_dir  "${DEEPCSA_RUN_DIR}" \
    --run_config "${RUN_CONFIG}" \
    --output_dir "${OUTPUT_DIR}/maf"

# generate VCFs

CORES=$(($(nproc) - 1))

for file in ${OUTPUT_DIR}/maf/*.tsv; do
    python ~/mutcatalog_simulator/deepcsa_maf2samplevcfs.py \
        --mutations-file "$file" \
        --output-dir "${OUTPUT_DIR}/vcf" &

    # This block limits the number of background jobs
    if [[ $(jobs -r | wc -l) -ge $CORES ]]; then
        wait -n  # Wait for at least one job to finish
    fi

done

wait