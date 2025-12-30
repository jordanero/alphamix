#!/bin/bash

## Most simulation experiments proceed in five distinct slurm job submissions.
# The first four correspond to different steps of the generative model, and the
# last corresponds to inference. The generative model is split across four
# job submissions allow for different amounts of parallelization in different
# steps.

BASE_DIR=/n/groups/price/jordan/h2xancestry

#### MAIN SIMS?
python ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/simulate_effects.py
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_genetic_components_chr.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_pheno_from_genetic_components.sh
bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/start_gwas_jobs.sh 
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_25_traits_w_unconstrained.sh

#### Sims varying T
###### T = MAC1
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/simulate_effects.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/get_genetic_components_chr.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/get_pheno_from_genetic_components.sh
bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/start_gwas_jobs.sh 
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/get_alpha_mix_25_traits_w_unconstrained.sh

###### T = 0.05
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/simulate_effects.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/get_genetic_components_chr.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/get_pheno_from_genetic_components.sh
bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/start_gwas_jobs.sh 
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/get_alpha_mix_25_traits_w_unconstrained.sh


#### Sims varying polygenicity
###### n_causal = 5000
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/simulate_effects.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/get_genetic_components_chr.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/get_pheno_from_genetic_components.sh
bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/start_gwas_jobs.sh 
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/get_alpha_mix_25_traits_w_unconstrained.sh

###### n_causal = 20000
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/simulate_effects.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/get_genetic_components_chr.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/get_pheno_from_genetic_components.sh
bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/start_gwas_jobs.sh 
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/get_alpha_mix_25_traits_w_unconstrained.sh

#### Sims doing single trait analyses
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_genetic_components_chr.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_pheno_from_genetic_components.sh
bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/start_gwas_jobs.sh 
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_single_trait_w_unconstrained.sh

#### sims using UKBB british MAF for generative model
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/simulate_effects.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/get_genetic_components_chr.sh
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/get_pheno_from_genetic_components.sh
bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/start_gwas_jobs.sh 
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/get_alpha_mix_25_traits_w_unconstrained.sh

# common domain
sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_25_traits_w_unconstrained_common_domain.sh
