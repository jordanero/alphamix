#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=100000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH --array=0-49
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_common_domain_trait_specific_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/get_alpha_from_grid_w_unconstrained_common_domain_trait_specific_%A_%a.err


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
trait = traits[${SLURM_ARRAY_TASK_ID}]
grid_filter_func = lambda g: float(g.split('_')[-1].split(',')[0][1:]) >= 0.04;
ldsc_utils.get_alpha_from_grid_across_traits_pmix_thresholded(
    results_dir = '${LDSC_RESULTS_DIR}',
    traits = [trait],
    in_sample_annotation_prefix_list = ['${LDSCORE_PREFIX}', '${BASELINELD_PREFIX}'], 
    grid_annotation_prefix = '/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_grid/maf_grid',
    out_path_prefix = f'${OUT_DIR}/{trait}_w_unconstrained_common_domain',
    european_maf_prefix = '/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC',
    african_maf_prefix = '/n/groups/price/jordan/h2xancestry/data/af/1000G_EUR_Phase3/1000G.EUR.QC',
    african_maf_suffix = 'MAF_afr_v2.txt',
    pmix_threshold = None,
    allow_w_out_of_range = True,
    grid_name_filter_function = grid_filter_func,
    n_chrom = 22
)"