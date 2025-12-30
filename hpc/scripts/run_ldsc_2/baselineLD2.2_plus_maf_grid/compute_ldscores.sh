#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-00:30                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=20000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH --array=1-22
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/run_ldsc/1000G_EUR_Phase3/baselineLD2.2_plus_maf_grid/compute_ldscores_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/run_ldsc/1000G_EUR_Phase3/baselineLD2.2_plus_maf_grid/compute_ldscores_%A_%a.err

LDSC_DIR=/n/groups/price/jordan/pkgs/ldsc
BED_PREFIX=/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC
CHROM=$SLURM_ARRAY_TASK_ID
ANNOT_PREFIX=/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_grid/maf_grid
REGRESSION_SNP_KEEP_FILE=/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/list.txt

python ${LDSC_DIR}/ldsc.py \
    --print-snps ${REGRESSION_SNP_KEEP_FILE} \
    --ld-wind-cm 1 \
    --out ${ANNOT_PREFIX}.${CHROM} \
    --bfile ${BED_PREFIX}.${CHROM} \
    --annot ${ANNOT_PREFIX}.${CHROM}.annot.gz \
    --l2
