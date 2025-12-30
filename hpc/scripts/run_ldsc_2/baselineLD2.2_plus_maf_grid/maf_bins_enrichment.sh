#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=100000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH --array=0-49
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/maf_bins_enrichment_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/run_ldsc_2/baselineLD2.2_plus_maf_grid/maf_bins_enrichment_%A_%a.err

module load conda/miniforge3

conda run -p /home/jor6523/rh9_envs/multisusie python -c \
"import sys;
sys.path.append('/n/groups/price/jordan/h2xancestry/scripts');
import ldsc_utils;
import pandas as pd;
traits = pd.read_csv('/n/groups/price/jordan/h2xancestry/data/sumstats/jordans_trait_list.txt', sep = '\t').iloc[:,0].tolist()
trait = traits[${SLURM_ARRAY_TASK_ID}]
ldsc_utils.get_enrichment_fast(
    results_prefix_list = [f'/n/groups/price/jordan/h2xancestry/data/real_traits_2/baselineLD2.2_plus_maf_grid/{trait}'],
    in_sample_annotation_prefix_list = ['/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_grid/maf_grid', '/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/baselineLD2.2_no_maf_bins/baselineLD_no_maf_bins'], 
    out_of_sample_annotation_prefix = '/n/groups/price/jordan/h2xancestry/data/ldscores/1000G_EUR_Phase3/maf_bins_stratified_by_afr_maf_corrected/maf_bins_stratified_by_afr_maf_corrected',
    out_path_list = [f'/n/groups/price/jordan/h2xancestry/data/real_traits_2/baselineLD2.2_plus_maf_grid/{trait}.maf_bin_stratified_by_afr_maf_corrected_enrichment.txt'],
    maf_prefix = '/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC',
    pruning_labels = None,
    per_allele = True,
    n_chrom = 22,
    jackknife_M = True
)"