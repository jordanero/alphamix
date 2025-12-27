#!/bin/bash

SCRIPT_DIR=/n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/

for CHR in {2..22}; do
    NSNPS_2=$(wc -l /n/groups/price/UKBiobank/UKB_pgen/highinfo_snps/UKB.highinfo.2.bim | awk '{print $1}')
    NSNPS_CHR=$(wc -l /n/groups/price/UKBiobank/UKB_pgen/highinfo_snps/UKB.highinfo.${CHR}.bim | awk '{print $1}')
    hours=$(echo "scale=2; 15 * ${NSNPS_CHR} / ${NSNPS_2}" | bc | awk -F. '{print ($2 && $2!=0) ? $1+1 : $1}')

    cp ${SCRIPT_DIR}/gwas_chr1.sh ${SCRIPT_DIR}/gwas_chr${CHR}.sh
    sed -i "s/chr1/chr${CHR}/g" ${SCRIPT_DIR}/gwas_chr${CHR}.sh
    sed -i "s/CHROM=1/CHROM=${CHR}/g" ${SCRIPT_DIR}/gwas_chr${CHR}.sh
    sed -i "s/0-15:00/0-${hours}:00/g" ${SCRIPT_DIR}/gwas_chr${CHR}.sh
    if [ $hours -lt 13 ]; then
        sed -i "s/medium/short/g" ${SCRIPT_DIR}/gwas_chr${CHR}.sh
    fi
done

for CHR in {1..22}; do
    sbatch ${SCRIPT_DIR}/gwas_chr${CHR}.sh
done
