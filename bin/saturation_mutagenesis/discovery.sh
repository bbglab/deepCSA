#!/bin/bash
#SBATCH -t 100:00:00 
#SBATCH --mem 100gb
#SBATCH -J saturation_discovery
#SBATCH --array=0-2

FILES=("all_samples" CohortCha_TimepointT0 CohortCha_TimepointT1)
GROUP=${FILES[$SLURM_ARRAY_TASK_ID]}

source activate prominent_R4.4_irb

python discovery_with_params.py --group-name $GROUP --deepcsa-dir '/data/bbg/nobackup/prominent/chip/analysis/2026-04-14_CH_I_wSP_custom'
python discovery_with_params.py --residue --group-name $GROUP --deepcsa-dir '/data/bbg/nobackup/prominent/chip/analysis/2026-04-14_CH_I_wSP_custom'
python discovery_with_params.py --impact non_protein_affecting --group-name $GROUP --deepcsa-dir '/data/bbg/nobackup/prominent/chip/analysis/2026-04-14_CH_I_wSP_custom'
python discovery_with_params.py --impact non_protein_affecting --residue --group-name $GROUP --deepcsa-dir '/data/bbg/nobackup/prominent/chip/analysis/2026-04-14_CH_I_wSP_custom'