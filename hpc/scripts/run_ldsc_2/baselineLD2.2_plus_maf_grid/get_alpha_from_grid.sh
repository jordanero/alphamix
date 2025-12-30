#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
#SBATCH -p medium # Partition to run in
#SBATCH --mem=100000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_%A_%a.err

LDSCORE_PREFIX=/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_grid/maf_grid.
BASELINELD_PREFIX=/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/baselineLD2.2_no_maf_bins/baselineLD_no_maf_bins.
LDSCORE_PREFIX="${LDSCORE_PREFIX%?}"
BASELINELD_PREFIX="${BASELINELD_PREFIX%?}"
LDSC_RESULTS_DIR=/n/groups/price/jordan/h2xancestry/data/real_traits_2/baselineLD2.2_plus_maf_grid
OUT_DIR=/n/groups/price/jordan/h2xancestry/data/real_traits_2/baselineLD2.2_plus_maf_grid

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
    out_path_prefix = '${OUT_DIR}/nh2_selected',
    european_maf_prefix = '/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC',
    african_maf_prefix = '/n/groups/price/jordan/h2xancestry/data/af/1000G_EUR_Phase3/1000G.EUR.QC',
    african_maf_suffix = 'MAF_afr_v2.txt',
    pmix_threshold = 0
)"
