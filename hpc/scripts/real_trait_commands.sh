#!/bin/bash
BASE_DIR=${BASE_DIR}
# main real trait results
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/run_ldsc.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_trait_specific.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_common_domain.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_common_domain_trait_specific.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_pass_traits.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_ukbiobank_traits.sh

# baselineLD no maf bins
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2/run_ldsc.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2/get_alpha_from_grid_w_unconstrained.sh

# YRI maf bins
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid_yri/compute_ldscores.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid_yri/run_ldsc.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid_yri/get_alpha_from_grid_w_unconstrained.sh

# MAF grid only
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/maf_grid/run_ldsc.sh
#sbatch ${BASE_DIR}/scripts/run_ldsc_2/maf_grid/get_alpha_from_grid_w_unconstrained.sh

# common domain
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_25_traits_w_unconstrained_common_domain.sh
