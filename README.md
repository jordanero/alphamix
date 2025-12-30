# Per-allele disease and complex trait effect sizes are predominantly African MAF-dependent in European populations

This directory contains scripts and notebooks used to perform the analyses described in "Per-allele disease and complex trait effect sizes are predominantly African MAF-dependent in European populations". 

The analyses performed were split across the All of Us Researcher Workbench and an academic high performance computing cluster. To make this directory easy to parse, I've removed lots of code which didn't produce a result included in the final manuscript. It's possible that some code was left out in this process, so please submit an issue in this repository or email Jordan Rossen (jordanerossen@gmail.com) for questions regarding reproducibility or code usage.

The file structure of this directory is set up accordingly:

```
.
├── all_of_us
│   ├── notebooks/
│   └── scripts/
├── hpc
│   ├── ldsc/
│   ├── notebooks/
│   └── scripts/
├── local/
└── README.md
```

``all_of_us`` contains scripts and jupyter notebooks which were run on the All of Us (AoU) Researcher Workbench. These correspond to analyses of non-synoymous variation.

``hpc`` contains scripts and notebooks which were run on the Harvard Medical School O2 computing cluster. These correspond to analyses of GWAS effects. ``hpc`` additionally contains a modified version of the linkage-disequilibrium score regression (LDSC) software which was used in this project.

``local`` contains R scripts which were run locally to make figures for analyses described in both ``all_of_us`` and ``hpc``.

## All of Us Researcher Workbench analyses
``all_us_of`` contains jupyter notebooks and scripts which were run on the AoU Researcher Workbench. These notebooks focus on analyses of non-synonymous variation and correspond to the first two figures of the preprint. The directory looks like this:
```
├── notebooks
│   ├── 01_parse_vat.ipynb
│   ├── 02_make_plinks.ipynb
│   ├── 03_get_unadmixed_allele_frequencies.ipynb
│   └── 04_analyze_nonsynonymous_variation.ipynb
└── scripts
    └── get_unadmixed_allele_frequencies.py
```
Each jupyter notebook in ``all_of_us`` is described below with references to files in ``scripts``. 

1. ``01_parse_vat.ipynb`` - Extracts variant metadata from the All of Us variant annotation table using dsub. This is used to select SNPs based on MAF in ``02_make_plinks.ipynb`` and to create the non-synonymous status annotation used in ``04_analyze_nonsynonymous_variation.ipynb``. Parsing this table is tricky because it's very large. 

2. ``02_make_plinks.ipynb`` - Creates ancestry-specific plink .bed/.fam/.bim/.frq files with reduced variant sets using dsub. Inputs are the AoU preprocessed ACAF v8 plink files (available to all controlled tier users), MAFs from the parsed variant annotation table from ``01_parse_vat.ipynb``, and plink sample keep files.

3. ``03_get_unadmixed_allele_frequencies.ipynb`` - Estimates unadmixed African allele frequencies by restricting allele frequency calculation to African ancestry individuals with two African haplotypes locally using dsub. Inputs are RFMix2 local ancestry calls (we use local ancestry calls from Mandla and Shi et al. 2025 medRxiv), and plink and MAF files from ``02_make_plinks.ipynb``. The dsub call in this notebook runs ``get_unadmixed_allele_frequencies.py`` from ``scripts``. This notebook modifies the merged maf parquet file from ``02_make_plinks.ipynb``

4. ``04_analyze_nonsynonymous_variation.ipynb`` - Performs analyses to create data plotted in manuscript figures. Headers detail which figure each data corresponds to in v1 of preprint. Inputs are the maf file from ``02_make_plinks.ipynb`` and non-synonymous annotations from ``01_parse_vat.ipynb``.

## High-performance computing cluster analyses
``hpc`` contains python code and jupyter notebooks which was run on the Havard Medical School O2 HPC cluster. These analyses focus on real and simulated GWAS data and correspond to figures 3-6 of the preprint. ``hpc`` looks like this:
```
├── ldsc/
├── notebooks/
└── scripts
    ├── real_trait_commands.sh
    ├── simulation_commands.sh
    ├── ldsc_utils.py
    ├── run_ldsc_2/
    ├── run_QC_YRI1000G.sh
    └── simulations_2/
```

In ``hpc``, the bulk of the work is done in ``scripts``, not ``analyses`` (unlike in ``all_of_us``). ``hpc`` also contains ``ldsc`` which is a modified copy of the linkage-disequilbrium score regression (ldsc) software package used in this project.

If you're just interested in seeing how $\alpha_{mix}$ is fit, see the ``get_alpha_from_grid_across_traits_pmix_thresholded`` function in ``ldsc_utils.py``. If you're interested in understanding in greater detail, read on.


``scripts`` contains two master bash scripts, ``real_trait_commands.sh`` and ``simulation_commands.sh``. These scripts contain commands that start slurm jobs to perform most of the analyses for Figure 3-6. ``real_trait_commands.sh`` contains commands to analyze real GWAS summary statistics, and primarily starts slurm jobs to run scripts in ``run_ldsc_2``. ``simulation_commands.sh`` contains commands to analyze simulated GWAS summary statistics and primarily starts slurm jobs to run scripts in ``run_ldsc_2``.

``ldsc_utils.py`` contains functions to fit the $\alpha_{mix}$ model and more.

Note that many scripts contain hard-coded absolute paths to the computing cluster environment. Users adapting this code will need to update these paths to match their local directory structure.

``ldsc`` is a modified copy of the ldsc software package. The primary difference is that it has an additional flag ``--set-jk-seps-wo-filtering`` which splits regression SNPs into jackknife blocks before filtering based on presence in the summary statistics file. This allows SNPs to be split into jackknife blocks which are consistent across GWAS traits, even if the GWAS traits contain different SNPs. 

## Local analyses
``local`` contains R scripts used for plotting (``figure_*.R``) and a script used to liftover All of Us GRCh38 allele frequencies to GRCh37.