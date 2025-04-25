#!/bin/bash

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <downsample_proportion> <num_runs>"
    exit 1
fi


## Fill all the required variables below,
#  and check if you also need to define a specific profile
#  in order to use your parameters of choice, otherwise you can use alternative ways of defining them
#  (e.g., using a config file or passing them directly in the command line)

DOWNSAMPLE_PROPORTION=$1
NUM_RUNS=$2
CONFIG_FILE="slurm_executor.config"
OUTDIR_BASE="downsampling"

# Format downsample proportion to match desired suffix (e.g., 0.1 -> 01, 0.05 -> 005)
SUFFIX=$(echo "$DOWNSAMPLE_PROPORTION" | sed 's/\.//')  # Remove "0." prefix

WORK_DIR="regressions_${SUFFIX}"
echo "Working directory: $WORK_DIR"

mkdir -p $WORK_DIR;

cd $WORK_DIR;


# Run the first iteration and capture the session ID
OUTDIR="${OUTDIR_BASE}/${SUFFIX}_1"
echo "Running pipeline iteration 1 with downsample_proportion=${DOWNSAMPLE_PROPORTION}"

nextflow run bbglab/deepCSA \
    -profile singularity,irbcluster \
    -work-dir "${WORK_DIR}" \
    -c "${CONFIG_FILE}" \
    --downsample_proportion "${DOWNSAMPLE_PROPORTION}" \
    --outdir "${OUTDIR}" \
    -revision downsampling \
    -ansi-log false \
    -latest


SESSION_ID=$(grep 'Session UUID' .nextflow.log | rev | cut -d ' ' -f 1 | rev)

echo "Captured Session ID: $SESSION_ID"

# Run the pipeline for the remaining runs
for i in $(seq 2 "$NUM_RUNS"); do
    OUTDIR="${OUTDIR_BASE}/${SUFFIX}_${i}"
    echo "Running pipeline iteration $i with downsample_proportion=${DOWNSAMPLE_PROPORTION}"
    
    nextflow run bbglab/deepCSA \
        -profile singularity,irbcluster \
        -work-dir "${WORK_DIR}" \
        -c "${CONFIG_FILE}" \
        --downsample_proportion "${DOWNSAMPLE_PROPORTION}" \
        --outdir "${OUTDIR}" \
        -revision downsampling \
        -resume $SESSION_ID

done

echo "All $NUM_RUNS runs completed successfully!"

# Remove work directories after execution
rm -r ${WORK_DIR}*
echo "Cleaned up work directories."
