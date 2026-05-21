#!/bin/bash

MUTCATALOG_SIMULATOR_HOME="/home/fmuinos/deepCSA/assets/mutcatalog_simulator"

# DEEPCSA_RUN_DIR="/data/bbg/nobackup/prominent/SantPauCH/deepCSA/sp_vhio/custom/2026-03-23_Sp_CHa_H_T0s_only_custom"
# DEEPCSA_RUN_DIR="/data/bbg/nobackup/prominent/SantPauCH/deepCSA/sp_vhio/custom/2026-04-13_Sp_vhio_H_T0sonly_custom07"
DEEPCSA_RUN_DIR="/data/bbg/nobackup/prominent/SantPauCH/deepCSA/sp_vhio/custom/2026-05-05_sp_vhio_H_T0sonly_custom"
RUN_CONFIG="${MUTCATALOG_SIMULATOR_HOME}/test_config_blood.json"
OUTPUT_DIR="${DEEPCSA_RUN_DIR}/fake_mutations_poisson_blood"

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