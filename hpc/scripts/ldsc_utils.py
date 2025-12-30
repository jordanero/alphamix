import numpy as np
import pandas as pd
import tqdm

BASE_DIR='/n/groups/price/jordan/h2xancestry'
KG_PLINK_DIR='/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files'
UKB_PGEN_DIR='/n/groups/price/UKBiobank/UKB_pgen/highinfo_snps/MAFs'


def get_alpha_from_grid_across_traits_pmix_thresholded(
    results_dir, traits, in_sample_annotation_prefix_list, grid_annotation_prefix,
    out_path_prefix, european_maf_prefix, african_maf_prefix, african_maf_suffix, n_chrom = 22,
    grid_name_filter_function = None, allow_w_out_of_range = False, pmix_threshold = None,
    alpha_mix_suffix = 'alpha_mix', get_residuals = False, null_alpha_mix = False, workers = 1,
    skip_jackknife = False):

    from scipy.optimize import differential_evolution, minimize_scalar

    pmix_getter = Pmix_getter(
        grid_annotation_prefix = grid_annotation_prefix,
        european_maf_prefix = european_maf_prefix,
        african_maf_prefix = african_maf_prefix,
        african_maf_suffix = african_maf_suffix,
        n_chrom = n_chrom,
        block_path = results_dir + '/' + traits[0] + '.blocks',
        pmix_threshold = pmix_threshold
    )

    expectation_b2_labels_list = []
    expectation_b2_jk_list = []
    expectation_b2_list = []
    print('Calculating grid b2 expectations')
    for i in tqdm.tqdm(range(len(traits))):
        expectation_b2_labels, expectation_b2_jk, expectation_b2 = get_enrichment_fast(
            results_prefix_list = [results_dir + '/' + traits[i]],
            in_sample_annotation_prefix_list = in_sample_annotation_prefix_list,
            out_of_sample_annotation_prefix = grid_annotation_prefix,
            out_path_list = [None],
            maf_prefix = european_maf_prefix,
            pruning_labels = None,
            n_chrom = n_chrom,
            per_allele = True,
            return_b2_jk_estimates = True,
            jackknife_M = True 
        )

        # drop the base annotation label
        assert expectation_b2_labels[0] == 'base'
        expectation_b2_labels = expectation_b2_labels[1:]
        expectation_b2_jk = expectation_b2_jk[1:,:]
        expectation_b2 = expectation_b2[1:]

        if grid_name_filter_function is not None:
            grid_mask = np.flatnonzero([grid_name_filter_function(g) for g in expectation_b2_labels])
            expectation_b2_labels = [expectation_b2_labels[i] for i in grid_mask]
            expectation_b2_jk = expectation_b2_jk[grid_mask,:]
            expectation_b2 = expectation_b2[grid_mask]

        expectation_b2_labels_list.append(expectation_b2_labels)
        expectation_b2_jk_list.append(expectation_b2_jk)
        expectation_b2_list.append(expectation_b2)

        assert all([l == l0 for (l, l0) in zip(expectation_b2_labels, expectation_b2_labels_list[0])])
        


    if grid_name_filter_function is not None:
        pmix_getter.apply_grid_name_filter_function(grid_name_filter_function)
        assert all([l == l0 for (l, l0) in zip(pmix_getter.grid_labels, expectation_b2_labels_list[0])])
    
    print('Estimating alpha-mix parameters')
    if not null_alpha_mix:
        alpha_point_estimate = differential_evolution(
            func = alpha_mix_loss_across_traits,
            bounds = (
                (-1.5, 1.5),
                (-0.5, 1.5)
            ),
            args = (pmix_getter, expectation_b2_list, allow_w_out_of_range),
            workers = workers
        )
    else:
        loss = alpha_mix_loss_across_traits([0, 0], pmix_getter, expectation_b2_list, allow_w_out_of_range)
        alpha_point_estimate = Opt(0, 0, loss)


    if skip_jackknife:
        pd.DataFrame({
            'alpha' : [alpha_point_estimate.x[0]],
            'w' : alpha_point_estimate.x[1],
            'loss' : alpha_point_estimate.fun
        }).to_csv(
            f'{out_path_prefix}.{alpha_mix_suffix}.skip_jackknife.txt', 
            sep = '\t',
            index = False
        )
        return

    n_blocks = expectation_b2_jk_list[0].shape[1]
    jk_result_list = []
    for jk_idx in tqdm.tqdm(range(n_blocks), desc = 'jackknifing alpha-mix parameters'):
        pmix_getter.set_jk_idx(jk_idx)
        if not null_alpha_mix:
            alpha_jk_estimate = differential_evolution(
                func = alpha_mix_loss_across_traits,
                bounds = (
                    (-1.5, 1.5),
                    (-0.5, 1.5)
                ),
                args = (pmix_getter, [e[:,jk_idx] for e in expectation_b2_jk_list]),
                workers = workers
            )
        else:
            loss = alpha_mix_loss_across_traits([0, 0], pmix_getter, [e[:,jk_idx] for e in expectation_b2_jk_list], allow_w_out_of_range)
            alpha_jk_estimate = Opt(0, 0, loss)
        alpha, w = alpha_jk_estimate.x
        loss = alpha_jk_estimate.fun
        jk_result_list.append(pd.DataFrame({
            'jk_idx' : [jk_idx],
            'alpha' : alpha,
            'w' : w,
            'loss' : loss,
        }))

    min_param_df = pd.concat(jk_result_list)
    min_param_df.assign(
        fake = 1
    ).groupby(
        'fake'
    ).apply(lambda df: pd.DataFrame({
        'alpha_mean' : [get_jk_mean_and_se(df.alpha)[0]],
        'alpha_se' : get_jk_mean_and_se(df.alpha)[1],
        'alpha_point_estimate' : alpha_point_estimate.x[0],
        'w_mean' : get_jk_mean_and_se(df.w)[0],
        'w_se' : get_jk_mean_and_se(df.w)[1],
        'w_point_estimate' : alpha_point_estimate.x[1],
        'loss_mean' : [get_jk_mean_and_se(df.loss)[0]],
        'loss_se' : get_jk_mean_and_se(df.loss)[1],
        'loss_point_estimate' : alpha_point_estimate.fun
    })
    ).to_csv(
        f'{out_path_prefix}.{alpha_mix_suffix}.txt', 
        sep = '\t',
        index = False
    )

    if null_alpha_mix:
        return

    varying_w_param_df_list = []
    n_w_points = 101
    w_array = np.linspace(0, 1, n_w_points)
    for w_idx in tqdm.tqdm(range(n_w_points), desc = 'estimating alpha with fixed w'):
        pmix_getter.set_jk_idx(None)
        w = w_array[w_idx]

        alpha_point_estimate = minimize_scalar(
            fun = lambda alpha, *args: alpha_mix_loss_across_traits([alpha, w], *args),
            method = 'brent',
            args = (
                pmix_getter,
                expectation_b2_list,
                allow_w_out_of_range
            ),
            bracket = (-1.5, 1.5)
        )

        jk_result_list = []
        for jk_idx in range(n_blocks):
            pmix_getter.set_jk_idx(jk_idx)
            alpha_jk_estimate = minimize_scalar(
                fun = lambda alpha, *args: alpha_mix_loss_across_traits([alpha, w], *args),
                method = 'brent',
                args = (
                    pmix_getter,
                    [e[:,jk_idx] for e in expectation_b2_jk_list],
                    allow_w_out_of_range
                ),
                bracket = (-1.5, 1.5)
            )
            jk_result_list.append(pd.DataFrame({
                'jk_idx' : [jk_idx],
                'alpha' : alpha_jk_estimate.x,
                'w' : w,
                'loss' : alpha_jk_estimate.fun
            }))

        jk_result_df = pd.concat(jk_result_list)
        varying_w_param_df_list.append(pd.DataFrame({
            'w' : [w],
            'alpha_mean' : get_jk_mean_and_se(jk_result_df.alpha)[0],
            'alpha_se' : get_jk_mean_and_se(jk_result_df.alpha)[1],
            'alpha_point_estimate' : alpha_point_estimate.x,
            'loss_mean' : [get_jk_mean_and_se(jk_result_df.loss)[0]],
            'loss_se' : get_jk_mean_and_se(jk_result_df.loss)[1],
            'loss_point_estimate' : alpha_point_estimate.fun
        }))

    varying_w_param_df = pd.concat(varying_w_param_df_list)
    varying_w_param_df.to_csv(
        f'{out_path_prefix}.{alpha_mix_suffix}_varying_w.txt', 
        sep = '\t',
        index = False
    )


class Pmix_getter:
    def __init__(self, grid_annotation_prefix, european_maf_prefix, african_maf_prefix, african_maf_suffix, n_chrom = 22, jackknife_M = False, block_path = None, moment_order = ['Expectation_P_A', 'Expectation_P_E','Expectation_P_E_P_A', 'Expectation_P_A_sq', 'Expectation_P_E_sq'], pmix_threshold = None):
        
        # get the annotation df 
        grid_labels = pd.read_csv(f'{grid_annotation_prefix}.1.annot.gz', sep = '\t').columns.tolist()[4:]
        grid_df_with_maf = make_annotation_df(
            annotation_prefix_list = [grid_annotation_prefix], 
            n_chrom = n_chrom, 
            european_maf_prefix = european_maf_prefix, 
            african_maf_prefix = african_maf_prefix, 
            african_maf_suffix = african_maf_suffix,
            pruning_labels = grid_labels,
            block_path = block_path
        ).assign(
            grid_bin = lambda df: pd.from_dummies(df[grid_labels], default_category = 'lf')
        )

        # get maf moments without thresholding p_mix
        maf_sums_leave_one_jk_in_df = grid_df_with_maf[
            ['block', 'grid_bin', 'MAF_eur', 'MAF_afr']
        ].groupby(
            ['block', 'grid_bin']
        ).apply(
            lambda df: pd.DataFrame({
                'Expectation_P_A': [np.sum(df.MAF_afr)],
                'Expectation_P_E': np.sum(df.MAF_eur),
                'Expectation_P_E_P_A': np.sum(df.MAF_afr * df.MAF_eur),
                'Expectation_P_A_sq' : np.sum(df.MAF_afr ** 2),
                'Expectation_P_E_sq' : np.sum(df.MAF_eur ** 2),
                'n' : df.shape[0],
                'min_P_A' : np.min(df.MAF_afr),
                'min_P_E' : np.min(df.MAF_eur),
                'max_P_A' : np.max(df.MAF_afr),
                'max_P_E' : np.max(df.MAF_eur)
            })
        ).droplevel(2)

        n_blocks = maf_sums_leave_one_jk_in_df.index.get_level_values('block').max() + 1
        maf_sums_leave_one_jk_in = np.zeros((n_blocks, len(grid_labels), len(moment_order)))
        n_leave_one_jk_in = np.zeros((n_blocks, len(grid_labels)))
        for i in range(n_blocks):
            for j in range(len(grid_labels)):
                if (i, grid_labels[j]) in maf_sums_leave_one_jk_in_df.index:
                    maf_sums_leave_one_jk_in[i, j] = maf_sums_leave_one_jk_in_df.loc[i,grid_labels[j]][moment_order].to_numpy()
                    n_leave_one_jk_in[i, j] = maf_sums_leave_one_jk_in_df.loc[i,grid_labels[j]]['n']

        maf_sums = np.sum(maf_sums_leave_one_jk_in, axis = 0)
        maf_sums_leave_one_jk_out = maf_sums - maf_sums_leave_one_jk_in

        n_sums = np.sum(n_leave_one_jk_in, axis = 0)
        n_sums_leave_one_jk_out = n_sums - n_leave_one_jk_in

        maf_means = maf_sums / n_sums[:,np.newaxis]
        maf_means_leave_one_jk_out = maf_sums_leave_one_jk_out / n_sums_leave_one_jk_out[:,:,np.newaxis]

        maf_maxs_and_mins = maf_sums_leave_one_jk_in_df.groupby(
            'grid_bin'
        ).apply(
            lambda df: pd.DataFrame({
                'max_P_A' : [np.max(df.max_P_A)],
                'min_P_A' : np.min(df.min_P_A),
                'max_P_E' : np.max(df.max_P_E),
                'min_P_E' : np.min(df.min_P_E)
            })
        ).loc[
            grid_labels
        ]

        if pmix_threshold == 0:
            pmix_threshold = None

        self.annotation_df_grouped = grid_df_with_maf.groupby('grid_bin')
        self.maf_means_leave_one_jk_out = maf_means_leave_one_jk_out
        self.maf_means = maf_means
        self.n_sums_leave_one_jk_out = n_sums_leave_one_jk_out
        self.n_sums = n_sums
        self.grid_labels = grid_labels
        self.moment_order = moment_order
        self.P_A_mins = maf_maxs_and_mins.min_P_A.to_numpy()
        self.P_E_mins = maf_maxs_and_mins.min_P_E.to_numpy()
        self.P_A_maxs = maf_maxs_and_mins.max_P_A.to_numpy()
        self.P_E_maxs = maf_maxs_and_mins.max_P_E.to_numpy()
        self.pmix_threshold = pmix_threshold
        self.jk_idx = None


    # how does jackknifing interact with this function?
    def get_pmix_moments(self, w):

        # check which cells need to be recomputed
        if self.pmix_threshold is not None:
            if w <= 1 and w >= 0:
                recomputation_needed = w * self.P_A_mins + (1 - w) * self.P_E_mins < self.pmix_threshold
            elif w > 1:
                recomputation_needed = w * self.P_A_mins + (1 - w) * self.P_E_maxs < self.pmix_threshold
            elif w < 0:
                recomputation_needed = w * self.P_A_maxs + (1 - w) * self.P_E_mins < self.pmix_threshold

    
        # get pmix times one minux pmix from precomputed moments
        if self.jk_idx is None:
            maf_means = self.maf_means
        else:
            maf_means = self.maf_means_leave_one_jk_out[self.jk_idx]

        P_E = maf_means[:,np.flatnonzero([g == 'Expectation_P_E' for g in self.moment_order])[0]]
        P_A = maf_means[:,np.flatnonzero([g == 'Expectation_P_A' for g in self.moment_order])[0]]
        P_E_P_A = maf_means[:,np.flatnonzero([g == 'Expectation_P_E_P_A' for g in self.moment_order])[0]]
        P_A_sq = maf_means[:,np.flatnonzero([g == 'Expectation_P_A_sq' for g in self.moment_order])[0]]
        P_E_sq = maf_means[:,np.flatnonzero([g == 'Expectation_P_E_sq' for g in self.moment_order])[0]]

        pmix_times_one_minus_pmix =  w * P_A + \
            (1 - w) * P_E - \
            (w ** 2) * P_A_sq - \
            (2 * w * (1 - w)) * P_E_P_A - \
            ((1 - w) ** 2) * P_E_sq


        # replace pmix(1-pmix) for moments that may be affected by thresholding
        if self.pmix_threshold is not None:
            for i in np.flatnonzero(recomputation_needed):
                annotation_df = self.annotation_df_grouped.get_group(self.grid_labels[i])
                if self.jk_idx is not None:
                    annotation_df = annotation_df.query(
                        f'block != {self.jk_idx}'
                    )
                pmix_i = w * annotation_df.MAF_afr + (1 - w) * annotation_df.MAF_eur
                pmix_i = np.maximum(pmix_i, self.pmix_threshold)
                pmix_times_one_minus_pmix_i = pmix_i * (1 - pmix_i)
                pmix_times_one_minus_pmix[i] = np.mean(pmix_times_one_minus_pmix_i)

        return pmix_times_one_minus_pmix
            
    def get_n_variants(self):
        if self.jk_idx is None:
            return self.n_sums
        else:
            return self.n_sums_leave_one_jk_out[self.jk_idx]

    def set_jk_idx(self, jk_idx):
        self.jk_idx = jk_idx

    def apply_grid_name_filter_function(self, grid_name_filter_function):
        grid_mask = np.flatnonzero([grid_name_filter_function(g) for g in self.grid_labels])
        self.grid_labels = [self.grid_labels[i] for i in grid_mask]
        self.maf_means_leave_one_jk_out = self.maf_means_leave_one_jk_out[:,grid_mask,:]
        self.maf_means = self.maf_means[grid_mask,:]
        self.n_sums_leave_one_jk_out = self.n_sums_leave_one_jk_out[:,grid_mask]
        self.n_sums = self.n_sums[grid_mask]
        self.annotation_df_grouped = self.annotation_df_grouped.filter(lambda x: x.grid_bin.iloc[0] in self.grid_labels)
        self.P_A_mins = self.P_A_mins[grid_mask]
        self.P_E_mins = self.P_E_mins[grid_mask]
        self.P_A_maxs = self.P_A_maxs[grid_mask]
        self.P_E_maxs = self.P_E_maxs[grid_mask]


def get_enrichment_fast(
    results_prefix_list, in_sample_annotation_prefix_list, out_of_sample_annotation_prefix,
    out_path_list, maf_prefix, pruning_labels = None, n_chrom = 22, per_allele = False, 
    return_b2_jk_estimates = False, return_b2_point_estimates = False, jackknife_M = False):
    
    from functools import reduce, partial


    M_overlap_chrom_list = []
    M_overlap_per_allele_chrom_list = []

    if jackknife_M:
        assert len(results_prefix_list) == 1
        jk_blocks = pd.read_csv(f'{results_prefix_list[0]}.blocks', sep = '\t').assign(block = lambda df: df.block.astype('Int64'))
        if n_chrom == 1: # this could break and is just for debugging
            n_blocks = 17
        else:
            n_blocks = jk_blocks.block.max() + 1
        M_overlap_leave_one_jk_in_list = [0 for _ in range(n_blocks)]
        M_overlap_per_allele_leave_one_jk_in_list = [0 for _ in range(n_blocks)]
    else:
        M_overlap_chrom_list = []
        M_overlap_per_allele_chrom_list = []

    for chrom in range(1,n_chrom + 1):

        # read and merge annotation dfs for the chromosome
        annotation_df = reduce(
            lambda df1, df2: df1.merge(df2, how = 'inner'), 
            [pd.read_csv(f'{prefix}.{chrom}.annot.gz', sep = '\t') for prefix in in_sample_annotation_prefix_list]
        ) 

        # prune annotation dfs for the chromosome to common variants
        if pruning_labels is not None:
            annotation_df = annotation_df[np.sum(annotation_df[pruning_labels], axis = 1) == 1]
        if maf_prefix is not None:
            if isinstance(maf_prefix, str):
                frq_function = lambda c: f'{maf_prefix}.{c}.frq'
            else:
                frq_function = maf_prefix
            maf_df = pd.read_csv(
                frq_function(chrom),
                sep = '\\s+'
            ).assign(
                MAF = lambda df: np.minimum(df.MAF, 1 - df.MAF)
            )
            assert all(maf_df.SNP == annotation_df.SNP)
            annotation_df = annotation_df[maf_df.MAF > .05]
            maf_df = maf_df[maf_df.MAF > .05]
            maf = maf_df.MAF.to_numpy()

            ldsc_annotation_labels = annotation_df.columns.tolist()[4:]
        
        if out_of_sample_annotation_prefix is not None:
            out_of_sample_annotation_df = pd.read_csv(
                f'{out_of_sample_annotation_prefix}.{chrom}.annot.gz', 
                sep = '\t'
            ).drop(
                columns = ['CM']
            )
            enrichment_annotation_labels = out_of_sample_annotation_df.columns.tolist()[3:]
            if 'base' not in enrichment_annotation_labels:
                enrichment_annotation_labels = ['base'] + enrichment_annotation_labels
            annotation_df = annotation_df.merge(out_of_sample_annotation_df, how = 'left')
        else:
            enrichment_annotation_labels = ldsc_annotation_labels

        if jackknife_M:
            annotation_df = annotation_df.merge(jk_blocks, how = 'left')
            annotation_df.iloc[annotation_df.shape[0]-1, np.flatnonzero(annotation_df.columns == 'block')] = n_blocks - 1
            annotation_df.block = annotation_df.block.bfill()
            maf_df = maf_df.merge(annotation_df[['SNP', 'block']], validate = 'one_to_one', indicator = True)
            assert(all(maf_df._merge == 'both'))
            maf_df = maf_df.drop(columns = ['_merge'])

        # count the number of variants in each intersection of annotations
        enrichment_annotations_numpy = annotation_df[
            enrichment_annotation_labels
        ].to_numpy()
        n_enrichment_annotations = enrichment_annotations_numpy.shape[1]
        ldsc_annotations_numpy = annotation_df[
            ldsc_annotation_labels
        ].to_numpy()
        n_ldsc_annotations = ldsc_annotations_numpy.shape[1]
        M_overlap_chrom_list.append(enrichment_annotations_numpy.T.dot(ldsc_annotations_numpy))

        if per_allele:
            variance = 2 * maf * (1 - maf)
            M_overlap_per_allele_chrom_list.append(np.dot(enrichment_annotations_numpy.T / variance, ldsc_annotations_numpy))

        if jackknife_M:
            for block in annotation_df.block.unique():
                enrichment_annotations_numpy_block = annotation_df[annotation_df.block == block][enrichment_annotation_labels].to_numpy()
                ldsc_annotations_numpy_block = annotation_df[annotation_df.block == block][ldsc_annotation_labels].to_numpy()
                M_overlap_jk = enrichment_annotations_numpy_block.T.dot(ldsc_annotations_numpy_block)
                M_overlap_leave_one_jk_in_list[block] += M_overlap_jk
                if per_allele:
                    maf_block = maf_df[annotation_df.block == block].MAF.to_numpy()
                    variance_block = 2 * maf_block * (1 - maf_block)
                    M_overlap_per_allele_jk = np.dot(enrichment_annotations_numpy_block.T / variance_block, ldsc_annotations_numpy_block)
                    M_overlap_per_allele_leave_one_jk_in_list[block] += M_overlap_per_allele_jk

    M_overlap = sum(M_overlap_chrom_list)
    M_overlap_per_allele = sum(M_overlap_per_allele_chrom_list)
    if jackknife_M:
        M_overlap_leave_one_jk_out = sum(M_overlap_leave_one_jk_in_list) - np.concatenate([np.expand_dims(M, 0) for M in M_overlap_leave_one_jk_in_list], axis = 0)
        if per_allele:
            M_overlap_per_allele_leave_one_jk_out = sum(M_overlap_per_allele_leave_one_jk_in_list) - np.concatenate([np.expand_dims(M, 0) for M in M_overlap_per_allele_leave_one_jk_in_list], axis = 0)
    else:
        M_overlap_leave_one_jk_out = np.concatenate([np.expand_dims(M, 0) for M in M_overlap_chrom_list], axis = 0)
        if per_allele:
            M_overlap_per_allele_leave_one_jk_out = np.concatenate([np.expand_dims(M, 0) for M in M_overlap_per_allele_chrom_list], axis = 0)

    if len(results_prefix_list) == 1:
        iter_object = range(1)
    else:
        iter_object = tqdm.tqdm(range(len(results_prefix_list)))

    for trait_index in iter_object:
        results_prefix = results_prefix_list[trait_index]
        out_path = out_path_list[trait_index]
        for line in open(f'{results_prefix}.log'):
            if line.startswith('Categories: '):
                coefficient_labels = line.strip().split('Categories: ')[1].split(' ')
                coefficient_labels = [c[:-4] for c in coefficient_labels]
                break
        jk_coefficients = pd.read_csv(
            f'{results_prefix}.part_delete', 
            sep = ' ', 
            header = None,
            names = coefficient_labels
        )[
            [c.replace('[', '').replace(']', '') for c in ldsc_annotation_labels]
        ].to_numpy().T[:,:n_blocks] # should be n_annot by n_jk 

        point_estimate_coefficients = pd.read_csv(
            f'{results_prefix}.results', 
            sep = '\t'
        ).assign(
            Category = lambda df: df.Category.str[:-4]
        ).set_index(
            'Category'
        ).loc[
            ldsc_annotation_labels, 'Coefficient'
        ].to_numpy()

        annotation_h2_jk = np.einsum('jel,lj->ej', M_overlap_leave_one_jk_out, jk_coefficients)
        #annotation_h2_jk = M_overlap.dot(jk_coefficients) # should be n_annot by n_jk
        annotation_h2_point_estimates = M_overlap.dot(point_estimate_coefficients)

        base_enrichment_idx = np.flatnonzero([c == 'base' for c in enrichment_annotation_labels])[0]
        base_ldsc_idx = np.flatnonzero([c == 'base' for c in ldsc_annotation_labels])[0]
        total_h2 = annotation_h2_point_estimates[base_enrichment_idx]
        total_h2_jk = annotation_h2_jk[base_enrichment_idx] 
        prop_h2 = annotation_h2_point_estimates / total_h2
        prop_h2_jk = annotation_h2_jk / total_h2_jk
        prop_snps = M_overlap[:,base_ldsc_idx] / M_overlap[base_enrichment_idx, base_ldsc_idx]
        prop_snps_jk = M_overlap_leave_one_jk_out[:,:,base_ldsc_idx].T / M_overlap_leave_one_jk_out[:,base_enrichment_idx, base_ldsc_idx]
        enrichment = prop_h2 / prop_snps
        enrichment_jk = prop_h2_jk / prop_snps_jk
        #enrichment_jk = np.transpose(prop_h2_jk.T / prop_snps)

        prop_h2_mean, prop_h2_se = zip(*[get_jk_mean_and_se(prop_h2_jk[i]) for i in range(prop_h2_jk.shape[0])])
        enrichment_mean, enrichment_se = zip(*[get_jk_mean_and_se(enrichment_jk[i]) for i in range(enrichment_jk.shape[0])])
        coefficient_mean, coefficient_se = zip(*[get_jk_mean_and_se(jk_coefficients[i]) for i in range(jk_coefficients.shape[0])])

        result_df = pd.DataFrame({
            'Category' : enrichment_annotation_labels,
            'prop_snps' : prop_snps,
            'prop_h2' : prop_h2_mean,
            'prop_h2_se' : prop_h2_se,
            'enrichment' : enrichment_mean,
            'enrichment_se' : enrichment_se,
            'prop_h2_point_estimate' : prop_h2,
            'enrichment_point_estimate' : enrichment
        })

        if per_allele:
            #expectation_b2_jk = np.transpose(M_overlap_per_allele.dot(jk_coefficients).T / M_overlap[:,base_ldsc_idx])
            #expectation_b2_point_estimates = np.transpose(M_overlap_per_allele.dot(point_estimate_coefficients).T / M_overlap[:,base_ldsc_idx])

            expectation_b2_jk = np.einsum('jel,lj->ej', M_overlap_per_allele_leave_one_jk_out, jk_coefficients) / M_overlap_leave_one_jk_out[:,:,base_ldsc_idx].T
            expectation_b2_point_estimates = np.dot(M_overlap_per_allele, point_estimate_coefficients) / M_overlap[:,base_ldsc_idx]

            if return_b2_jk_estimates:
                return (enrichment_annotation_labels, expectation_b2_jk, expectation_b2_point_estimates)
            if return_b2_point_estimates:
                return (enrichment_annotation_labels, expectation_b2_point_estimates)

            expectation_b2_normalized = expectation_b2_point_estimates / expectation_b2_point_estimates[base_enrichment_idx]
            expectation_b2_normalized_jk = expectation_b2_jk / expectation_b2_jk[base_enrichment_idx]

            expectation_b2_mean, expectation_b2_se = zip(*[get_jk_mean_and_se(expectation_b2_jk[i]) for i in range(n_enrichment_annotations)])
            expectation_b2_normalized_mean, expectation_b2_normalized_se = zip(*[get_jk_mean_and_se(expectation_b2_normalized_jk[i]) for i in range(n_enrichment_annotations)])

            result_df = result_df.assign(
                expectation_b2 = expectation_b2_mean,
                expectation_b2_se = expectation_b2_se,
                expectation_b2_normalized = expectation_b2_normalized_mean,
                expectation_b2_normalized_se = expectation_b2_normalized_se,
                expectation_b2_point_estimate = expectation_b2_point_estimates,
                expectation_b2_normalized_point_estimate = expectation_b2_normalized
            )
        else:
            result_df = result_df.assign(
                coefficient = coefficient_mean,
                coefficient_se = coefficient_se
            )

        result_df.to_csv(
            out_path, sep = '\t', index = False
        )


def get_jk_mean_and_se(jk_stats):
    '''
    jk_stats is a 1 dimensional vector. each element is the statistic 
    computed leaving out each block one at a time. it returns a list where 
    the first element is the jackknife mean and the second element is the
    jackknife standard error
    '''
    n = len(jk_stats)
    tn_bar = np.mean(jk_stats)
    Q = np.sum((jk_stats - tn_bar) ** 2)
    jk_var = Q * (n - 1) / n
    jk_se = np.sqrt(jk_var)
    return [tn_bar, jk_se]


def make_annotation_df(
    annotation_prefix_list,
    n_chrom = 22,
    european_maf_prefix = None,
    african_maf_prefix = None,
    african_maf_suffix = None,
    pruning_labels = None,
    sanitize_columns = False,
    block_path = None):

    from functools import reduce

    annotation_df_list = []

    for chrom in range(1, n_chrom + 1):
        annotation_df = reduce(lambda df1, df2: df1.merge(df2, how = 'inner'), [pd.read_csv(f'{prefix}.{chrom}.annot.gz', sep = '\t') for prefix in annotation_prefix_list]) 
        if sanitize_columns:
            annotation_df.columns = [c.replace('[', '').replace(']', '') for c in annotation_df.columns]
        if pruning_labels is not None:
            annotation_df = annotation_df[np.sum(annotation_df[pruning_labels], axis = 1) == 1]
        if european_maf_prefix is not None:
            eur_maf_chrom = pd.read_csv(f'{european_maf_prefix}.{chrom}.frq', sep = '\\s+').assign(
                MAF_eur = lambda df: np.minimum(df.MAF, 1 - df.MAF)
            )
            annotation_df = annotation_df.merge(eur_maf_chrom[['SNP', 'MAF_eur']], how = 'left', validate = 'one_to_one')
        if african_maf_prefix is not None:
            afr_maf_chrom = pd.read_csv(f'{african_maf_prefix}.{chrom}.{african_maf_suffix}', sep = '\t').assign(
                MAF_afr = lambda df: np.minimum(df.MAF_afr, 1 - df.MAF_afr)
            )
            annotation_df = annotation_df.merge(afr_maf_chrom[['SNP', 'MAF_afr']], how = 'left', validate = 'one_to_one')
        annotation_df_list.append(annotation_df)
    annotation_df = pd.concat(annotation_df_list)

    if block_path is not None:
        jk_blocks = pd.read_csv(block_path, sep = '\t').assign(block = lambda df: df.block.astype('Int64'))
        annotation_df = annotation_df.merge(jk_blocks, how = 'left')
        #n_blocks = np.nanmax(annotation_df.block.to_numpy()) #why did this stop working?
        n_blocks = annotation_df.block.max(skipna = True) + 1
        annotation_df.iloc[annotation_df.shape[0]-1, np.flatnonzero(annotation_df.columns == 'block')] = n_blocks - 1
        annotation_df.block = annotation_df.block.bfill()

    return annotation_df


def get_maf_df():
    maf_df_list = []
    for chrom in tqdm.tqdm(range(1,23)):
        afr_maf_chrom = pd.read_csv(
            f'{BASE_DIR}/data/af/1000G_EUR_Phase3/1000G.EUR.QC.{chrom}.MAF_afr_v2.txt',
            sep = '\t'
        )
        eur_maf_chrom = pd.read_csv(
            f'{KG_PLINK_DIR}/1000G.EUR.QC.{chrom}.frq',
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

    ukb_maf = pd.concat([pd.read_parquet(f'{UKB_PGEN_DIR}/mafs.British.{chrom}.parquet').assign(CHR = chrom) for chrom in range(1,23)])
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

    return ukb_intersection


class Opt:
    def __init__(self, alpha, c, loss):
        self.x = [alpha, c]
        self.fun = loss


def alpha_mix_predictions_2(params, pmix_getter, allow_w_out_of_range = False, pmix_times_one_minus_pmix = None):
    alpha, log_sigma_g_squared, w = params
    sigma_g_squared = np.exp(log_sigma_g_squared)
    if (w < 0 or w > 1) and not allow_w_out_of_range:
        return np.nan
    if pmix_times_one_minus_pmix is None:
        pmix_times_one_minus_pmix = pmix_getter.get_pmix_moments(pmix_threshold = pmix_getter.pmix_threshold, w = w)
    if (pmix_times_one_minus_pmix <= 0).any():
        return np.nan
    return sigma_g_squared * ((2 * pmix_times_one_minus_pmix) ** alpha)


def alpha_mix_loss_2(params, pmix_getter, expect_b2, error_scaling = 1e12, allow_w_out_of_range = False, pmix_times_one_minus_pmix = None, n_variants = None):
    if pmix_times_one_minus_pmix is not None and (pmix_times_one_minus_pmix <= 0).any():
        return error_scaling * 1e20
    predictions = alpha_mix_predictions_2(params, pmix_getter, allow_w_out_of_range = allow_w_out_of_range, pmix_times_one_minus_pmix = pmix_times_one_minus_pmix)
    if np.isnan(predictions).any():
        return error_scaling * 1e20
    if n_variants is None:
        n_variants = pmix_getter.get_n_variants()
    error = (predictions - expect_b2) ** 2 
    error_weighted = error * n_variants / np.sum(n_variants)
    return error_scaling * np.sum(error_weighted)


def alpha_mix_loss_across_traits(alpha_and_w, pmix_getter, expectation_b2_list, allow_w_out_of_range = False):
    alpha = alpha_and_w[0]
    w = alpha_and_w[1]
    loss = 0
    loss_per_variant = 0
    s = []
    pmix_times_one_minus_pmix = pmix_getter.get_pmix_moments(w)
    if any(pmix_times_one_minus_pmix <= 0):
        return 1e20
    n_variants = pmix_getter.get_n_variants()
    for trait_idx in range(len(expectation_b2_list)):
        sigma_g_squared = np.sum(n_variants * expectation_b2_list[trait_idx] * ((2 * pmix_times_one_minus_pmix) ** alpha)) /\
            np.sum(n_variants * ((2 * pmix_times_one_minus_pmix) ** (2 * alpha)))
        if sigma_g_squared <= 0:
            return 1e20
        loss += alpha_mix_loss_2(params = [alpha, np.log(sigma_g_squared), w], pmix_getter = None, expect_b2 = expectation_b2_list[trait_idx], error_scaling = 1e12, allow_w_out_of_range = allow_w_out_of_range, pmix_times_one_minus_pmix = pmix_times_one_minus_pmix, n_variants = n_variants)
    if not np.isfinite(loss):
        return 1e20
    return loss


def simulate_effects(sim_suffix: str, n_causal_variants: int = 10000, pmix_threshold: float = 0.005, use_YRI_MAF = False, use_ukbb_british_maf = False):
    """Simulate effect sizes across a grid of w values.

    Parameters
    ----------
    sim_suffix : str
        Suffix for the simulation directory.
    n_causal_variants : int, default 10000
        Number of causal variants to sample per simulation.
    pmix_threshold : float, default 0.005
        Minimum threshold applied to p_mix before computing variance.
    """
    import numpy as np
    import pandas as pd
    import tqdm
    import os

    W_ARRAY = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]
    EFFECT_SIZE_ROOT_DIR = f'{BASE_DIR}/data/simulations/effect_sizes_2'
    os.makedirs(EFFECT_SIZE_ROOT_DIR, exist_ok = True)

    maf_df_list = []
    for chrom in tqdm.tqdm(range(1,23)):
        afr_maf_chrom = pd.read_csv(
            f'{BASE_DIR}/data/af/1000G_EUR_Phase3/1000G.EUR.QC.{chrom}.MAF_afr_v2.txt',
            sep = '\t'
        )
        eur_maf_chrom = pd.read_csv(
            f'{KG_PLINK_DIR}/1000G.EUR.QC.{chrom}.frq',
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

    ukb_maf = pd.concat([pd.read_parquet(f'{UKB_PGEN_DIR}/mafs.British.{chrom}.parquet').assign(CHR = chrom) for chrom in range(1,23)])
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


    if use_YRI_MAF:
        ukb_intersection_eur_common = ukb_intersection_eur_common.drop(
            columns = ['MAF_afr']
        ).rename(
            columns = {'MAF_YRI' : 'MAF_afr'}
        )
    if use_ukbb_british_maf:
        ukb_intersection_eur_common = ukb_intersection_eur_common.drop(
            columns = ['MAF_eur']
        ).assign(
            MAF_eur = lambda df: df.MAF_ukb
        )

    alpha = -.38
    h2 = 0.5
    n_causal = n_causal_variants

    if pmix_threshold == 'MAC1':
        pmix_threshold = ukb_intersection_eur_common.query('MAF_eur > 0.05').query('MAF_afr > 0').MAF_afr.min()

    for w in W_ARRAY:
        os.makedirs(f'{EFFECT_SIZE_ROOT_DIR}/w_{w}_{sim_suffix}', exist_ok = True)

        pd.DataFrame({'f' : [f'{EFFECT_SIZE_ROOT_DIR}/w_{w}_{sim_suffix}/seed_{seed}_beta.tsv' for seed in range(1000)]}).to_csv(
            f'{EFFECT_SIZE_ROOT_DIR}/w_{w}_{sim_suffix}/effect_size_file_list.txt',
            sep = '\t',
            header = False,
            index = False
        )


        for seed in tqdm.tqdm(range(1000)):
            rng = np.random.default_rng(seed)
            causal_effects = ukb_intersection_eur_common.sample(int(n_causal)).assign(
                p_mix = lambda df: df.MAF_afr * w + df.MAF_eur * (1 - w)
            ).assign(
                p_mix = lambda df: np.maximum(df.p_mix, pmix_threshold),
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
                f'{EFFECT_SIZE_ROOT_DIR}/w_{w}_{sim_suffix}/seed_{seed}_beta.tsv',
                sep = '\t',
                header = False,
                index = False
            )


# Not yet in paper but likely will be
def get_alpha_from_snps_across_traits_true_effects(seed_start, seed_stop, w, out_path):
    from scipy.optimize import differential_evolution

    maf_df = get_maf_df()

    beta_df_list = []
    for seed in range(seed_start, seed_stop + 1):
        beta_df = pd.read_csv(
            f'{BASE_DIR}/data/simulations/effect_sizes_2/w_{w}_alpha_.38_h2_.5_thresholded_.005/seed_{seed}_beta.tsv', 
            sep = '\t',
            header = None,
            names = ['SNP_ukb', 'A1', 'BETA']
        ).assign(
            seed = seed
        )
        beta_df_list.append(beta_df)
    beta_df = pd.concat(
        beta_df_list
    ).merge(
        maf_df,
        how = 'left'
    ).assign(
        beta_sq = lambda df: df.BETA ** 2
    )
    betasq_list = []
    p_a_list = []
    p_e_list = []
    for seed in range(seed_start, seed_stop + 1):
        current_beta_df = beta_df.query('seed == @seed')
        betasq_list.append(current_beta_df.beta_sq.to_numpy())
        p_a_list.append(current_beta_df.MAF_afr.to_numpy())
        p_e_list.append(current_beta_df.MAF_eur.to_numpy())

    alpha_point_estimate = differential_evolution(
        func = alpha_mix_loss_across_traits_snps,
        bounds = (
            (-1.5, 1.5),
            (-0.5, 1.5)
        ),
        args = (betasq_list, p_a_list, p_e_list),
    )

    pd.DataFrame({
        'alpha_hat' : [alpha_point_estimate.x[0]],
        'w_hat' : [alpha_point_estimate.x[1]],
        'loss' : alpha_point_estimate.fun
    }).to_csv(
        out_path, 
        sep = '\t',
        index = False
    )


# Not yet in paper but likely will be
def get_alpha_from_grid_across_traits_true_effects(seed_start, seed_stop, w, out_path):
    from scipy.optimize import differential_evolution, minimize_scalar

    maf_df = get_maf_df()

    grid_annotations = pd.read_csv(f'{BASE_DIR}/data/ldscores/1000G_EUR_Phase3/maf_grid/maf_grid.22.annot.gz', sep = '\t')
    afr_maf_bins = sorted(list(set([float(s.split('_')[-1].split(',')[1][:-1]) for s in grid_annotations.columns[grid_annotations.columns.str.contains('afr_maf_')]])))
    eur_maf_bins = sorted(list(set([float(s.split('_')[-4].split(',')[1][:-1]) for s in grid_annotations.columns[grid_annotations.columns.str.contains('eur_maf_')]])))
    afr_maf_bins = [-1] + afr_maf_bins + [1]
    eur_maf_bins = [-1] + eur_maf_bins + [1]

    beta_df_list = []
    for seed in range(seed_start, seed_stop + 1):
        beta_df = pd.read_csv(
            f'{BASE_DIR}/data/simulations/effect_sizes_2/w_{w}_alpha_.38_h2_.5_thresholded_.005/seed_{seed}_beta.tsv', 
            sep = '\t',
            header = None,
            names = ['SNP_ukb', 'A1', 'BETA']
        ).assign(
            seed = seed
        )
        beta_df_list.append(beta_df)
    beta_df = pd.concat(
        beta_df_list
    ).merge(
        maf_df,
        how = 'left'
    ).assign(
        beta_sq = lambda df: df.BETA ** 2
    ).assign(
        afr_maf_bin = lambda df: pd.cut(df.MAF_afr, bins = afr_maf_bins, labels = range(len(afr_maf_bins) - 1)),
        eur_maf_bin = lambda df: pd.cut(df.MAF_eur, bins = eur_maf_bins, labels = range(len(eur_maf_bins) - 1))
    ).assign(
        combined_bin = lambda df: df.afr_maf_bin.astype(str) + '_' + df.eur_maf_bin.astype(str),
    )

    # each list one list per trait. Each trait specific list contains one list per bin
    betasq_list = []
    p_a_list = []
    p_e_list = []
    n_list = []

    for seed in range(seed_start, seed_stop + 1):
        current_trait_beta_df = beta_df.query('seed == @seed')
        combined_bins = current_trait_beta_df.combined_bin.unique()
        current_trait_beta_df_grouped = current_trait_beta_df.groupby('combined_bin')
        current_trait_betasq_list = []
        current_trait_p_a_list = []
        current_trait_p_e_list = []
        current_trait_n_list = []

        for combined_bin in combined_bins:
            current_trait_beta_df_bin = current_trait_beta_df_grouped.get_group(combined_bin)
            current_trait_betasq_list.append(current_trait_beta_df_bin.beta_sq.to_numpy())
            current_trait_p_a_list.append(current_trait_beta_df_bin.MAF_afr.to_numpy())
            current_trait_p_e_list.append(current_trait_beta_df_bin.MAF_eur.to_numpy())
            current_trait_n_list.append(current_trait_beta_df_bin.shape[0])
        betasq_list.append(current_trait_betasq_list)
        p_a_list.append(current_trait_p_a_list)
        p_e_list.append(current_trait_p_e_list)
        n_list.append(current_trait_n_list)

    alpha_point_estimate = differential_evolution(
        func = alpha_mix_loss_across_traits_without_precomputing_moments,
        bounds = (
            (-1.5, 1.5),
            (-0.5, 1.5)
        ),
        args = (betasq_list, p_a_list, p_e_list, n_list),
    )

    pd.DataFrame({
        'alpha_hat' : [alpha_point_estimate.x[0]],
        'w_hat' : [alpha_point_estimate.x[1]],
        'loss' : alpha_point_estimate.fun
    }).to_csv(
        out_path, 
        sep = '\t',
        index = False
    )


# Not yet in paper but likely will be
def alpha_mix_loss_across_traits_without_precomputing_moments(alpha_and_w, betasq_list, p_a_list, p_e_list, n_list):
    MAF_FOR_MAC_1 = 1.489e-05
    alpha = alpha_and_w[0]
    w = alpha_and_w[1]
    loss = 0
    loss_per_variant = 0
    s = []


    for trait_idx in range(len(betasq_list)):
        trait_pmix_times_one_minus_pmix_bin_mean_list = []
        trait_bin_n_list = []
        trait_betasq_mean_list = []
        for bin_idx in range(len(p_a_list[trait_idx])):
            p_a = p_a_list[trait_idx][bin_idx]
            p_e = p_e_list[trait_idx][bin_idx]
            p_mix = p_a * w + p_e * (1 - w)
            pmix = np.maximum(p_mix, .005)
            pmix_times_one_minus_pmix = pmix * (1 - pmix)
            trait_pmix_times_one_minus_pmix_bin_mean_list.append(np.mean(pmix_times_one_minus_pmix))
            trait_bin_n_list.append(n_list[trait_idx][bin_idx])
            trait_betasq_mean_list.append(np.mean(betasq_list[trait_idx][bin_idx]))

        trait_pmix_times_one_minus_pmix_means = np.array(trait_pmix_times_one_minus_pmix_bin_mean_list)
        trait_bin_ns = np.array(trait_bin_n_list)
        trait_betasq_means = np.array(trait_betasq_mean_list)

        sigma_g_squared = np.sum(trait_bin_ns * trait_betasq_means * ((2 * trait_pmix_times_one_minus_pmix_means) ** alpha)) /\
            np.sum(trait_bin_ns * ((2 * trait_pmix_times_one_minus_pmix_means) ** (2 * alpha)))
        if sigma_g_squared <= 0:
            return 1e20
        preds = sigma_g_squared * ((2 * trait_pmix_times_one_minus_pmix_means) ** alpha)
        error = (preds - trait_betasq_means) ** 2
        error_weighted = error * trait_bin_ns / np.sum(trait_bin_ns)
        loss += np.sum(error_weighted)


    if not np.isfinite(loss):
        assert False
    return loss


# Not yet in paper but likely will be
def alpha_mix_loss_across_traits_snps(alpha_and_w, betasq_list, p_a_list, p_e_list):
    MAF_FOR_MAC_1 = 1.489e-05
    alpha = alpha_and_w[0]
    w = alpha_and_w[1]
    loss = 0
    loss_per_variant = 0
    s = []

    # INSERT PMIX CALCULATIONS HERE

    for trait_idx in range(len(betasq_list)):
        p_mix = p_a_list[trait_idx] * w + p_e_list[trait_idx] * (1 - w)
        pmix = np.maximum(p_mix, .005)
        pmix_times_one_minus_pmix = pmix * (1 - pmix)

        # this never trigers because of the pmix thresholding above
        if any(pmix_times_one_minus_pmix <= 0):
            return 1e20

        sigma_g_squared = np.sum(betasq_list[trait_idx] * ((2 * pmix_times_one_minus_pmix) ** alpha)) /\
            np.sum(((2 * pmix_times_one_minus_pmix) ** (2 * alpha)))
        if sigma_g_squared <= 0:
            return 1e20
        preds = sigma_g_squared * ((2 * pmix_times_one_minus_pmix) ** alpha)
        loss += np.sum((preds - betasq_list[trait_idx]) ** 2)
    if not np.isfinite(loss):
        return 1e20
    return loss
