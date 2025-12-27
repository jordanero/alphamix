#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-06:00                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=20000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH --array=23,45#0-153
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_genetic_components_chr_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_genetic_components_chr_%A_%a.err

W_ARRAY=(0 0.05 0.25 0.5 0.75 0.95 1)
CHROM=$(($SLURM_ARRAY_TASK_ID % 22))
CHROM=$((CHROM + 1))
w=${W_ARRAY[$SLURM_ARRAY_TASK_ID / 22]}

SCRATCH_DIR=/n/scratch/users/j/jor6523/alpha_mix/simulations_2/w_${w}_alpha_.38_h2_.5_thresholded_.005
EFFECT_SIZE_ROOT_DIR=/n/groups/price/jordan/h2xancestry/data/simulations/effect_sizes_2/w_${w}_alpha_.38_h2_.5_thresholded_.005

mkdir -p ${SCRATCH_DIR}

/n/groups/price/jordan/pkgs/plink2 \
    --bpfile /n/groups/price/UKBiobank/UKB_pgen/highinfo_snps/UKB.highinfo.${CHROM}\
    --score-list ${EFFECT_SIZE_ROOT_DIR}/effect_size_file_list.txt cols=+scoresums zs\
    --memory 19000 \
    --out ${SCRATCH_DIR}/chr${CHROM}_scores 
