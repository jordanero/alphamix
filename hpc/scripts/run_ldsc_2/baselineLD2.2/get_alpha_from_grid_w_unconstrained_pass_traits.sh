#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=100000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2/get_alpha_from_grid_w_unconstrained_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2/get_alpha_from_grid_w_unconstrained_%A_%a.err

BASELINELD_PREFIX=/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/baselineLD.
BASELINELD_PREFIX="${BASELINELD_PREFIX%?}"
LDSC_RESULTS_DIR=/n/groups/price/jordan/h2xancestry/data/real_traits_2/baselineLD2.2
OUT_DIR=/n/groups/price/jordan/h2xancestry/data/real_traits_2/baselineLD2.2

python -c "import sys;
sys.path.append('/n/groups/price/jordan/h2xancestry/scripts');
import ldsc_utils;
import pandas as pd;
traits = pd.read_csv('/n/groups/price/jordan/h2xancestry/data/sumstats/jordans_trait_list.txt', sep = '\t').iloc[:,0];
ldsc_utils.get_alpha_from_grid_across_traits_pmix_thresholded(
    results_dir = '${LDSC_RESULTS_DIR}',
    traits = traits,
    in_sample_annotation_prefix_list = ['${LDSCORE_PREFIX}', '${BASELINELD_PREFIX}'], 
    grid_annotation_prefix = '/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_grid/maf_grid',
    out_path_prefix = '${OUT_DIR}/nh2_selected_w_unconstrained',
    european_maf_prefix = '/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC',
    african_maf_prefix = '/n/groups/price/jordan/h2xancestry/data/af/1000G_EUR_Phase3/1000G.EUR.QC',
    african_maf_suffix = 'MAF_afr_v2.txt',
    pmix_threshold = None,
    allow_w_out_of_range = True
)"
