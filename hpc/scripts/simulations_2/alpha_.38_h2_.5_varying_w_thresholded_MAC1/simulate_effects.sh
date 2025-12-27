#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=100000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/simulate_effects_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_thresholded_MAC1/simulate_effects_%A_%a.err

set -e
set -x
module load conda/miniforge3

SIM_SUFFIX='alpha_.38_h2_.5_thresholded_MAC1'

conda run -p /home/jor6523/rh9_envs/multisusie python -c \
"import sys;
sys.path.append('/n/groups/price/jordan/h2xancestry/scripts');
import ldsc_utils;
ldsc_utils.simulate_effects(
    sim_suffix = '${SIM_SUFFIX}',
    pmix_threshold = 'MAC1' 
)"