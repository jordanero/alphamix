import numpy as np
import pandas as pd
import tqdm
import os
import plotly.express as px
from scipy.optimize import differential_evolution
import sys
sys.path.append('/n/groups/price/jordan/h2xancestry/scripts')
import ldsc_utils

W_ARRAY = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]
EFFECT_SIZE_ROOT_DIR = '/n/groups/price/jordan/h2xancestry/data/simulations/effect_sizes_2'
os.makedirs(EFFECT_SIZE_ROOT_DIR, exist_ok = True)



maf_df_list = []
for chrom in tqdm.tqdm(range(1,23)):
    afr_maf_chrom = pd.read_csv(
        f'/n/groups/price/jordan/h2xancestry/data/af/1000G_EUR_Phase3/1000G.EUR.QC.{chrom}.MAF_afr_v2.txt',
        sep = '\t'
    )
    eur_maf_chrom = pd.read_csv(
        f'/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.{chrom}.frq',
        sep = '\\s+'
    ).rename(
        columns = {'MAF' : 'MAF_eur'}
    ).assign(
        MAF_afr = afr_maf_chrom.MAF_afr,
        MAF_YRI = afr_maf_chrom.MAF_YRI,
    ).assign(
        MAF_afr = lambda df: np.minimum(df.MAF_afr, 1 - df.MAF_afr),
        MAF_yri = lambda df: np.minimum(df.MAF_YRI, 1 - df.MAF_YRI),
    ).assign(
        MAF_eur = lambda df: np.minimum(df.MAF_eur, 1 - df.MAF_eur)
    )

    maf_df_list.append(eur_maf_chrom)
maf_df = pd.concat(
    maf_df_list
)

ukb_maf = pd.concat([pd.read_parquet(f'/n/groups/price/UKBiobank/UKB_pgen/highinfo_snps/MAFs/mafs.British.{chrom}.parquet').assign(CHR = chrom) for chrom in range(1,23)])
ukb_intersection_maf_forward = ukb_maf.rename(
    columns = {'SNP' : 'SNP_ukb'}
).assign(
    SNP = ukb_maf.SNP.str.split('.', expand = True).iloc[:,0],
    A1 = ukb_maf.SNP.str.split('.', expand = True).iloc[:,2],
    A2 = ukb_maf.SNP.str.split('.', expand = True).iloc[:,1]
).rename(
    columns = {'MAF' : 'MAF_ukb'}
).merge(
    maf_df,
    how = 'inner',
).assign(
    match = 'forward'
)
ukb_intersection_maf_flipped = ukb_maf.rename(
    columns = {'SNP' : 'SNP_ukb'}
).assign(
    SNP = ukb_maf.SNP.str.split('.', expand = True).iloc[:,0],
    A1 = ukb_maf.SNP.str.split('.', expand = True).iloc[:,1],
    A2 = ukb_maf.SNP.str.split('.', expand = True).iloc[:,2]
).rename(
    columns = {'MAF' : 'MAF_ukb'}
).merge(
    maf_df,
    how = 'inner',
).assign(
    match = 'flipped'
)
A1 = ukb_intersection_maf_flipped.A1
ukb_intersection_maf_flipped.A1 = ukb_intersection_maf_flipped.A2
ukb_intersection_maf_flipped.A2 = A1
ukb_intersection = pd.concat([ukb_intersection_maf_forward, ukb_intersection_maf_flipped])
ukb_intersection_eur_common = ukb_intersection.query('MAF_eur > 0.05')

alpha = -.38
h2 = 0.5
n_causal = 10000

for w in W_ARRAY:
    os.makedirs(f'{EFFECT_SIZE_ROOT_DIR}/w_{w}_alpha_.38_h2_.5_thresholded_.005', exist_ok = True)

    pd.DataFrame({'f' : [f'{EFFECT_SIZE_ROOT_DIR}/w_{w}_alpha_.38_h2_.5_thresholded_.005/seed_{seed}_beta.tsv' for seed in range(1000)]}).to_csv(
        f'{EFFECT_SIZE_ROOT_DIR}/w_{w}_alpha_.38_h2_.5_thresholded_.005/effect_size_file_list.txt',
        sep = '\t',
        header = False,
        index = False
    )

    for seed in tqdm.tqdm(range(1000)):
        rng = np.random.default_rng(seed)
        causal_effects = ukb_intersection_eur_common.sample(int(n_causal)).assign(
            p_mix = lambda df: df.MAF_afr * w + df.MAF_eur * (1 - w)
        ).assign(
            p_mix = lambda df: np.maximum(df.p_mix, .005),
            variance = lambda df: (df.p_mix * (1 - df.p_mix)) ** alpha,
            effect = lambda df: rng.normal(0, np.sqrt(df.variance), size = len(df))
        ).assign(
            effect = lambda df: [0 if not np.isfinite(effect) else effect for effect in df.effect],
            expected_h2 = lambda df: 2 * df.MAF_eur * (1 - df.MAF_eur) * (df.effect ** 2),
            scaled_effect = lambda df: df.effect * np.sqrt(h2 / np.sum(df.expected_h2)),
            expected_h2_post_scaling = lambda df: 2 * df.MAF_eur * (1 - df.MAF_eur) * (df.scaled_effect ** 2),
        )
        causal_effects[
            ['SNP_ukb', 'A1', 'scaled_effect']
        ].to_csv(
            f'{EFFECT_SIZE_ROOT_DIR}/w_{w}_alpha_.38_h2_.5_thresholded_.005/seed_{seed}_beta.tsv',
            sep = '\t',
            header = False,
            index = False
        )
