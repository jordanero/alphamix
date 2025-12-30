#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=15000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH --array=1-50
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid_yri/run_ldsc_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid_yri/run_ldsc_%A_%a.err
#SBATCH --dependency=afterok:22362736


LDSC_DIR=/n/groups/price/jordan/pkgs/ldsc
SS_DIR=/n/groups/price/ldsc/sumstats_formatted_2024/sumstats/
LDSCORE_PREFIX=/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_grid_yri/maf_grid_yri.
BASELINELD_PREFIX=/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/baselineLD2.2_no_maf_bins/baselineLD_no_maf_bins.
OUT_DIR=/n/groups/price/jordan/h2xancestry/data/real_traits_2/baselineLD2.2_plus_maf_grid_yri
WEIGHTS_PREFIX=/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC.
BED_PREFIX=/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC
TRAIT_LIST=/n/groups/price/jordan/h2xancestry/data/sumstats/jordans_trait_list.txt

LINE_NUM=$((${SLURM_ARRAY_TASK_ID}+1))
TRAIT_LINE=$(head ${TRAIT_LIST} -n ${LINE_NUM} | tail -n 1)
IFS=$'\t' read -ra TRAIT <<< "$TRAIT_LINE"


module load conda/miniforge3
conda run -p /home/jor6523/.conda/envs/ldsc python \
    ${LDSC_DIR}/ldsc.py \
    --set-jk-seps-wo-filtering \
    --h2 ${SS_DIR}/${TRAIT}.sumstats \
    --ref-ld-chr ${LDSCORE_PREFIX},${BASELINELD_PREFIX} \
    --w-ld-chr ${WEIGHTS_PREFIX} \
    --overlap-annot \
    --frqfile-chr ${BED_PREFIX}. \
    --print-coefficients \
    --print-delete-vals \
    --out ${OUT_DIR}/${TRAIT}