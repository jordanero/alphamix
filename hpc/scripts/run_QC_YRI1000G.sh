#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-04:00                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=40000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/one_time_scipts/run_QC_YRI1000G_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/one_time_scipts/run_QC_YRI1000G_%A_%a.err

module load plink2/1.90

for CHR in {1..22};
do
  POP="YRI"
  DIR=/n/groups/price/jordan/h2xancestry/data/plink/1000G_${POP}

  awk '{if ($3=="YRI") {print $1, $1} }' /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/TGP2261.txt > ${DIR}/list_id.txt
  wget https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 
  mv ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${DIR}/
  gunzip ${DIR}/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  awk '{print $2}' /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.${CHR}.bim > ${DIR}/1000G.EUR.QC.${CHR}.snplist

  # make files that contain the 1000G eur reference snps
  plink \
    --vcf ${DIR}/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \
    --keep ${DIR}/list_id.txt \
    --keep-allele-order \
    --extract ${DIR}/1000G.EUR.QC.${CHR}.snplist \
    --make-bed \
    --out ${DIR}/1000G.$POP.$CHR

  plink \
    --keep-allele-order \
    --bfile ${DIR}/1000G.$POP.$CHR \
    --freq \
    --out ${DIR}/1000G.$POP.$CHR  
  done
