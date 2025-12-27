import pandas as pd
import numpy as np
import admix
import subprocess
import argparse

parser = argparse.ArgumentParser("simple_example")
parser.add_argument('--pvar_path', type = str)
parser.add_argument('--psam_path', type = str)
parser.add_argument('--sample_keep_path', type = str)
parser.add_argument('--lanc_path', type = str)
parser.add_argument('--bfile_prefix', type = str)
parser.add_argument('--chromosome', type = str)
parser.add_argument('--out_dir', type = str)
args = parser.parse_args()

sample_keep = pd.read_csv(args.sample_keep_path, header = None, sep = ' ')
keep_id_set = set(sample_keep.iloc[:,1])
indiv_meta = pd.read_csv(args.psam_path, sep = '\t', dtype = int, usecols = [0])
in_keep_idx = np.flatnonzero([i in keep_id_set for i in indiv_meta.iloc[:,0]])

lanc = admix.io.read_lanc(args.lanc_path)
afr_dosage = np.sum(lanc.dask()[:, in_keep_idx] == 0, axis = 2).compute()

snp_meta = pd.read_csv(args.pvar_path, sep = '\t')

bp_start = 0
n_blocks = 1
index_start = 0


for i in range(afr_dosage.shape[0] - 1):
    if (not np.all(afr_dosage[i] == afr_dosage[i+1])):
        bp_stop = snp_meta.iloc[i + 1].POS - 1
        unadmixed_idx = np.nonzero(np.all(afr_dosage[index_start:i+1] == 2, axis = 0))
        indiv_meta.iloc[in_keep_idx].iloc[unadmixed_idx].assign(
            FID = 0
        )[
            ['FID', '#IID']
        ].to_csv(
            'temp.keep',
            sep = ' ',
            header = False,
            index = False
        )


        subprocess.run(["plink",
                        "--keep-allele-order",
                        "--bfile", args.bfile_prefix,
                        "--memory", "10000",
                        "--freq",
                        "--keep", "temp.keep",
                        "--from-bp", str(int(bp_start)),
                        "--to-bp", str(int(bp_stop)),
                        "--chr", args.chromosome,
                        "--out", f'temp{n_blocks}'],
            stdout = subprocess.DEVNULL
        )

        bp_start = snp_meta.iloc[i + 1].POS
        index_start = i+1
        n_blocks += 1

# need to do one more for everything leftover
unadmixed_idx = np.nonzero(np.all(afr_dosage[index_start:] == 2, axis = 0))
indiv_meta.iloc[in_keep_idx].iloc[unadmixed_idx].assign(
    FID = 0
)[
    ['FID', '#IID']
].to_csv(
    'temp.keep',
    sep = ' ',
    header = False,
    index = False
)
subprocess.run(["plink",
                "--keep-allele-order",
                "--bfile", args.bfile_prefix,
                "--memory", "10000",
                "--freq",
                "--keep", "temp.keep",
                "--from-bp", str(int(bp_start)),
                "--chr", args.chromosome,
                "--out", f'temp{n_blocks}'])
n_blocks += 1

af_df_list = []
for i in range(1, n_blocks):
    af_df = pd.read_csv(f'temp{i}.frq', delim_whitespace = True)
    af_df_list.append(af_df)
af_df = pd.concat(af_df_list)
af_df.to_csv(
    f'{args.out_dir}/chr{args.chromosome}_afr_unadmixed.frq',
    sep = '\t',
    index = False
)
