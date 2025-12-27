#!/bin/bash
#SBATCH -c 4                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=20000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH --array=1-70
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/gwas_chr19_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/gwas_chr19_%A_%a.err
#SBATCH --dependency=afterok:21465331

CHROM=19
BLOCK=$(( ((SLURM_ARRAY_TASK_ID - 1) % 10)))
SEED_START=$(((BLOCK * 100) + 3))
SEED_END=$(((BLOCK * 100) + 102))
W_ARRAY=(0 0.05 0.25 0.5 0.75 0.95 1)
W=${W_ARRAY[(($SLURM_ARRAY_TASK_ID - 1) / 10)]}

SCRATCH_DIR=/n/scratch/users/j/jor6523/alpha_mix/simulations_2/w_${W}_alpha_.38_h2_.5_thresholded_.05

/n/groups/price/jordan/pkgs/plink2 \
    --bpfile /n/groups/price/UKBiobank/UKB_pgen/highinfo_snps/UKB.highinfo.${CHROM} \
    --keep /n/groups/price/UKBiobank/UKB_pgen/keep_files/British.norel.keep.txt \
    --pheno ${SCRATCH_DIR}/phenos.txt \
    --glm zs omit-ref allow-no-covars \
    --no-input-missing-phenotype \
    --memory 19000 \
    --threads 4 \
    --pheno-col-nums ${SEED_START}-${SEED_END} \
    --out ${SCRATCH_DIR}/chr${CHROM}
