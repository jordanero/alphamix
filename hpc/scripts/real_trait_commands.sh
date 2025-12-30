#!/bin/bash

BASE_DIR=

# Make MAF annotations
# consider adding AoU jupyter notebooks to get unadmixed allele frequencies from
# AoU v7.1
# consider adding the liftover script to go from AoU hg38 to baselineLD GRCh37
sbatch ${BASE_DIR}/scripts/run_QC_YRI1000G.sh
# run code in jupyter notebook under "Get unadmixed maf for 1000G EUR reference snps" in ${BASE_DIR}/scripts/make_v2_african_maf_annotations.ipynb

# Make figure 3
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/compute_ldscores.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/run_ldsc.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/coarse_grid_enrichment.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/maf_bins_enrichment.sh

# main real trait results
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_trait_specific.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_common_domain.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_common_domain_trait_specific.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_pass_traits.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_ukbiobank_traits.sh

# baselineLD no maf bins
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2/run_ldsc.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2/get_alpha_from_grid_w_unconstrained.sh

# YRI maf bins
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid_yri/compute_ldscores.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid_yri/run_ldsc.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid_yri/get_alpha_from_grid_w_unconstrained.sh

# MAF grid only
sbatch ${BASE_DIR}/scripts/run_ldsc_2/maf_grid/run_ldsc.sh
sbatch ${BASE_DIR}/scripts/run_ldsc_2/maf_grid/get_alpha_from_grid_w_unconstrained.sh

