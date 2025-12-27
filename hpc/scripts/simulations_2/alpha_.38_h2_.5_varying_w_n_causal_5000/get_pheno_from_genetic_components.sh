#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node
#SBATCH -t 0-03:00                         # Runtime in D-HH:MM format
#SBATCH -p short # Partition to run in
#SBATCH --mem=70000 # Memory total in MB (for all cores)
#SBATCH --account=price
#SBATCH --array=0-6
#SBATCH -o /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/get_pheno_from_genetic_components_%A_%a.out
#SBATCH -e /n/groups/price/jordan/h2xancestry/scripts/simulations_2/alpha_.38_h2_.5_varying_w_n_causal_5000/get_pheno_from_genetic_components_%A_%a.err
#SBATCH --dependency=afterok:21737081

module load conda/miniforge3
W_ARRAY=(0 0.05 0.25 0.5 0.75 0.95 1)
W=${W_ARRAY[$SLURM_ARRAY_TASK_ID]}

SCRATCH_DIR=/n/scratch/users/j/jor6523/alpha_mix/simulations_2/w_${W}_alpha_.38_h2_.5_n_causal_5000

conda run -p /home/jor6523/rh9_envs/multisusie python -c "import pandas as pd; import tqdm; import numpy as np; w = ${W};
brits = pd.read_csv('/n/groups/price/UKBiobank/UKB_pgen/keep_files/British.norel.keep.txt', sep = '\t', header = None).iloc[:, 0]
for chr in tqdm.tqdm(range(1,23)):
    score_chr = pd.read_csv(
        f'${SCRATCH_DIR}/chr{chr}_scores.sscore.zst',
        sep = '\t'
    ).set_index(
        'IID'
    ).loc[
        brits
    ].reset_index(
    )
    if chr == 1:
        iid = score_chr.IID
        score_sum = np.zeros((score_chr.shape[0], 1000))
    assert score_chr.IID.equals(iid)
    score_sum += score_chr[[f'SCORE{seed + 1}_SUM' for seed in range(1000)]].to_numpy()

genetic_variance = np.var(score_sum, axis = 0)
rng = np.random.default_rng(1001)
env_component = rng.normal(0, 1, score_sum.shape)
env_component_scaled = env_component * np.sqrt(genetic_variance / np.var(env_component, axis = 0))
full_phenotype = score_sum + env_component_scaled

pheno_df = pd.DataFrame(
    full_phenotype, 
    columns = [f'seed_{i}' for i in range(1000)]
).assign(
    FID = brits,
    IID = brits 
)
pheno_df[pheno_df.columns[1000:].tolist() + pheno_df.columns[:1000].tolist()].to_csv(
    f'${SCRATCH_DIR}/phenos.txt',
    sep = '\t', 
    index = False
)
"