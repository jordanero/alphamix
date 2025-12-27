# Per-allele disease and complex trait effect sizes are predominantly African MAF-dependent in European populations

This directory contains scripts and notebooks used to perform the analyses described in "Per-allele disease and complex trait effect sizes are predominantly African MAF-dependent in European populations". 

The analyses performed were split across the All of Us Researcher Workbench and an academic high performance computing cluster. To make this directory easy to parse, I've extracted code used to produce all published results. Please submit an issue in this repository or contact the corresponding author for questions regarding reproducibility or code usage.

The file structure of this directory is set up accordingly:

```
├── all_of_us
│   ├── notebooks
│   │   ├── 01_parse_vat.ipynb
│   │   ├── 02_make_plinks.ipynb
│   │   ├── 03_get_unadmixed_allele_frequencies.ipynb
│   │   └── 04_analyze_nonsynonymous_variation.ipynb
│   └── scripts
│       └── get_unadmixed_allele_frequencies.py
├── hpc
├── plotting
│   ├── figure_1.R
│   ├── figure_2.R
│   ├── figure_3_M_jk.R
│   ├── figure_4.R
│   ├── figure_5.R
│   ├── figure_6_supplementary.R
│   └── figure_6.R
└── README.md
```

``all_of_us`` contains scripts and jupyter notebooks which were run on the All of Us (AoU) Researcher Workbench. These correspond to analyses of non-synoymous variation.

``hpc`` contains scripts which were run on the Harvard Medical School O2 computing cluster. These correspond to analyses of GWAS effects.

``plotting`` contains R scripts which were run locally to make figures for analyses described in both ``all_of_us`` and ``hpc``.

## All of Us Researcher Workbench analyses
``all_us_of`` contains jupyter notebooks and scripts which were run on the AoU Researcher Workbench. These notebooks focus on analyses of non-synonymous variation and correspond to the first two figures of the preprint. We describe each jupyter notebook in ``all_of_us`` in turn, with references, as necessary to files in ``scripts``. Most of the Jupyter notebooks contain code which starts dsub jobs (dsub is the AoU workbench equivalent of a HPC job scheduler).

1. ``01_parse_vat.ipynb`` - Extracts variant metadata from the All of Us variant annotation table using dsub. This is used to select SNPs based on MAF in ``02_make_plinks.ipynb`` and to create the non-synonymous status annotation used in ``04_analyze_nonsynonymous_variation.ipynb``. Parsing this table is tricky because it's very large. 

2. ``02_make_plinks.ipynb`` - Creates ancestry-specific plink .bed/.fam/.bim/.frq files with reduced variant sets using dsub. Inputs are the AoU preprocessed ACAF v8 plink files (available to all controlled tier users), MAFs from the parsed variant annotation table from ``01_parse_vat.ipynb``, and plink sample keep files (not released because they might count as individual-level data).

3. ``03_get_unadmixed_allele_frequencies.ipynb`` - Estimates unadmixed African allele frequencies by restricting allele frequency calculation to African ancestry individuals with two African haplotypes locally using dsub. Inputs are RFMix2 local ancestry calls (we use local ancestry calls from Mandla and Shi et al. 2025 medRxiv), and plink and MAF files from ``02_make_plinks.ipynb``. The dsub call in this notebook runs ``get_unadmixed_allele_frequencies.py`` from ``scripts``. This script modifies the mergerd maf parquet file from ``02_make_plinks.ipynb``

4. ``04_analyze_nonsynonymous_variation.ipynb`` - Performs analyses to create data plotted in manuscript figures. Headers detail which figure each data corresponds to in v1 of preprint. Inputs are the maf file from ``02_make_plinks.ipynb`` and non-synonymous annotations which I HAVE TO FIND THE CODE FOR STILL!!!!

## High-performance computing cluster analyses
``hpc`` contains python code and jupyter notebooks which was run on the Havard Medical School O2 HPC cluster. These analyses focus on real and simulated GWAS data and correspond to figures 3-6 of the preprint. In ``hpc``, the bulk of the work is done in ``scripts``, not ``analyses`` (unlike in ``all_of_us``). 

``scripts`` contains two "master" bash scripts, ``simulation_commands.sh`` and ``real_trait_commands.sh``, which contain commands that start slurm jobs to perform most of the analyses for Figure 3-6. Most of the time, the commands in the master scripts submit scripts contained in ``scripts/simulations`` or ``scripts/run_ldsc_2``.  ``scripts/ldsc_utils.py`` contains functions to fit the $\alpha_{mix}$ model and more.

