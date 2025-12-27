#!/bin/bash

## Most simulation experiments proceed in five distinct slurm job submissions.
# The first four correspond to different steps of the generative model, and the
# last corresponds to inference. The generative model is split across four
# job submissions allow for different amounts of parallelization in different
# steps.

BASE_DIR=/n/groups/price/jordan/h2xancestry

## simulations varying w without thresholding
## ran code to generate effect sizes in ${BASE_DIR}/analyses/simulate_phenotypes.ipynb
##sbatch ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w/get_pheno.sh
##sbatch ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w/get_pheno_from_genetic_components.sh
##bash ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w/start_gwas_jobs.sh
##sbatch ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w/get_alpha_25_traits.sh
##sbatch ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w/get_alpha_single_trait.sh
#
## simulations with thresholding at .002
## ran code to generate effect sizes in ${BASE_DIR}/analyses/simulate_phenotypes.ipynb
##sbatch ${BASE_DIR}/scripts/simulations/w_1_alpha_.38_h2_.5_thresholded_.002/get_pheno.sh
##sbatch ${BASE_DIR}/scripts/simulations/w_1_alpha_.38_h2_.5_thresholded_.002/get_pheno_from_genetic_components.sh
##bash ${BASE_DIR}/scripts/simulations/w_1_alpha_.38_h2_.5_thresholded_.002/start_gwas_jobs.sh
#sbatch ${BASE_DIR}/scripts/simulations/w_1_alpha_.38_h2_.5_thresholded_.002/get_alpha_25_traits.sh
#
## simulations with thresholding at .05
## ran code to generate effect sizes in ${BASE_DIR}/analyses/simulate_phenotypes.ipynb
##sbatch ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w_thresholded_.05/get_pheno.sh
##sbatch ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w_thresholded_.05/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w_thresholded_.05/start_gwas_jobs.sh
##sbatch ${BASE_DIR}/scripts/simulations/alpha_.38_h2_.5_varying_w_thresholded_.05/get_alpha_25_traits.sh



## SIMULATIONS TAKE 2

#### MAIN SIMS?
#python ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/simulate_effects.py
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_genetic_components_chr.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/start_gwas_jobs.sh 
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_25_traits.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_25_traits_w_unconstrained.sh

#### Sims varying T
###### T = MAC1
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/simulate_effects.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/get_genetic_components_chr.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/start_gwas_jobs.sh 
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/get_alpha_mix_25_traits_w_unconstrained.sh
###### T = 0.05
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/simulate_effects.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/get_genetic_components_chr.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/start_gwas_jobs.sh 
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.05/get_alpha_mix_25_traits_w_unconstrained.sh


#### Sims varying polygenicity
###### n_causal = 5000
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/simulate_effects.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/get_genetic_components_chr.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/start_gwas_jobs.sh 
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/get_alpha_mix_25_traits_w_unconstrained.sh

###### n_causal = 20000
#OLD_DIR=${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/
#NEW_DIR=${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/
#mkdir -p ${NEW_DIR}
#for f in simulate_effects.sh get_genetic_components_chr.sh get_pheno_from_genetic_components.sh start_gwas_jobs.sh gwas_chr1.sh get_alpha_mix_25_traits_w_unconstrained.sh; do
#    cp ${OLD_DIR}/${f} ${NEW_DIR}/${f}
#    sed -i "s/5000/20000/g" ${NEW_DIR}/${f}
#done
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/simulate_effects.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/get_genetic_components_chr.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/start_gwas_jobs.sh 
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_20000/get_alpha_mix_25_traits_w_unconstrained.sh

#### Sims doing single trait analyses
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_genetic_components_chr.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/start_gwas_jobs.sh 
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_.005/get_alpha_mix_single_trait_w_unconstrained.sh

#### Sims using YRI MAF for generative model
#OLD_DIR=${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/
#NEW_DIR=${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_yri_maf/
#mkdir -p ${NEW_DIR}
#for f in simulate_effects.sh get_genetic_components_chr.sh get_pheno_from_genetic_components.sh start_gwas_jobs.sh gwas_chr1.sh get_alpha_mix_25_traits_w_unconstrained.sh; do
#    cp ${OLD_DIR}/${f} ${NEW_DIR}/${f}
#done
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_yri_maf/simulate_effects.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_yri_maf/get_genetic_components_chr.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_yri_maf/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_yri_maf/start_gwas_jobs.sh 
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_yri_maf/get_alpha_mix_25_traits_w_unconstrained.sh


#### sims using UKBB british MAF for generative model
#OLD_DIR=${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/
#NEW_DIR=${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/
#mkdir -p ${NEW_DIR}
#for f in simulate_effects.sh get_genetic_components_chr.sh get_pheno_from_genetic_components.sh start_gwas_jobs.sh gwas_chr1.sh get_alpha_mix_25_traits_w_unconstrained.sh; do
#    cp ${OLD_DIR}/${f} ${NEW_DIR}/${f}
#done
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/simulate_effects.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/get_genetic_components_chr.sh
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/get_pheno_from_genetic_components.sh
#bash ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/start_gwas_jobs.sh 
#sbatch ${BASE_DIR}/scripts/simulations_2/alpha_.38_h2_.5_varying_w_ukbb_british_maf/get_alpha_mix_25_traits_w_unconstrained.sh
