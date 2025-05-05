#!/bin/bash
#SBATCH -J SigProfilerExtractor
#SBATCH -o /data/bbg/nobackup/work/tests/sigprofiler_extractor_%j.out
#SBATCH -e /data/bbg/nobackup/work/tests/sigprofiler_extractor_%j.err
#SBATCH --cpus-per-task 20
#SBATCH -t 07:00:00
#SBATCH --mem 40096M
#SBATCH -p bbg_cpu_zen4,irb_cpu_zen4

. "/home/fcalvet/miniforge3/etc/profile.d/conda.sh"
. "/home/fcalvet/miniforge3/etc/profile.d/mamba.sh"

mamba activate sigprofiler_env

SIGPRO_script="/home/fcalvet/projects/deepCSA/assets/useful_scripts"



DEEPCSARUN="/data/bbg/nobackup/prominent/kidney/pancancerpanel/2025-04-30_LCMs"
signatures_output="${DEEPCSARUN}/sigprofilerextractor_auto_subset1"
matrix_name="${DEEPCSARUN}/matrixconcatwgs/subsets/2025-05-02.below03_artif.sp.tsv"

mkdir ${signatures_output};

${SIGPRO_script}/signatures_sigprofilerextractor.py matrix ${signatures_output} ${matrix_name} --ref_genome GRCh38 --max_sig 10 --nmf_replicates 100 --cpu 20

