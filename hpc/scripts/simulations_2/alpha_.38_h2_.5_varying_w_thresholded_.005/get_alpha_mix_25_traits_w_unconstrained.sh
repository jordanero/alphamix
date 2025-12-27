#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-4:00
#SBATCH -p short 
#SBATCH --mem=70000
#SBATCH --account=price
#SBATCH --array=25,37#0-279
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_25_traits_w_unconstrained_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_25_traits_w_unconstrained_%A_%a.err


set -e
set -x
module load conda/miniforge3

FIRST_SEED=$(( (SLURM_ARRAY_TASK_ID % 40) * 25 ))
LAST_SEED=$((FIRST_SEED + 24))
W_ARRAY=(0 0.05 0.25 0.5 0.75 0.95 1)
w=${W_ARRAY[$((SLURM_ARRAY_TASK_ID / 40))]}
echo ${w}
echo ${FIRST_SEED}
echo ${LAST_SEED}


SCRATCH_DIR=/n/scratch/users/j/jor6523/alpha_mix/simulations_2/w_${w}_alpha_.38_h2_.5_thresholded_.005
RESULTS_DIR=/n/groups/price/jordan/h2xancestry/data/simulations/alpha_mix_2/w_${w}_alpha_.38_h2_.5_thresholded_.005
LDSC_DIR=/n/groups/price/jordan/pkgs/ldsc
LDSCORE_PREFIX=/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_grid/maf_grid.
BASELINELD_PREFIX=/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/baselineLD2.2_no_maf_bins/baselineLD_no_maf_bins.
WEIGHTS_PREFIX=/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC.
BED_PREFIX=/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC

mkdir -p ${RESULTS_DIR}


## step 4 estimate alpha mix parameters

LDSCORE_PREFIX="${LDSCORE_PREFIX%?}"
BASELINELD_PREFIX="${BASELINELD_PREFIX%?}"

conda run -p /home/jor6523/rh9_envs/multisusie python -c \
"import sys;
sys.path.append('/n/groups/price/jordan/h2xancestry/scripts');
import ldsc_utils;
import pandas as pd;
ldsc_utils.get_alpha_from_grid_across_traits_pmix_thresholded(
    results_dir = '${RESULTS_DIR}',
    traits = [f'seed_{seed}' for seed in range(${FIRST_SEED}, ${LAST_SEED} + 1)],
    in_sample_annotation_prefix_list = ['${LDSCORE_PREFIX}', '${BASELINELD_PREFIX}'], 
    grid_annotation_prefix = '/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_grid/maf_grid',
    out_path_prefix = '${RESULTS_DIR}/seeds_${FIRST_SEED}_${LAST_SEED}',
    european_maf_prefix = '/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC',
    african_maf_prefix = '/n/groups/price/jordan/h2xancestry/data/af/1000G_EUR_Phase3/1000G.EUR.QC',
    african_maf_suffix = 'MAF_afr_v2.txt',
    pmix_threshold = 0,
    skip_jackknife = True,
    allow_w_out_of_range = True
)"
