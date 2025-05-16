# RhP–SEIRS Weibo Dataset & Reproduction Pipeline  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15378206.svg)](https://doi.org/10.5281/zenodo.15378206)
[![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15378372.svg)](https://doi.org/10.5281/zenodo.15378372)
[![License: MIT](https://img.shields.io/badge/Code-MIT-green.svg)](LICENSE)
[![Data License: CC BY-NC 4.0](https://img.shields.io/badge/Data-CC%20BY--NC%204.0-lightgrey.svg)](DATA_LICENSE)
---

## 1 Project overview
This repository contains **everything you need to reproduce every table, figure and
numerical result in the paper**.  
The core resources are

* **A nearly unbiased crawl of 8.7 million Weibo posts** (2020-2025)  
  *full text removed, user-IDs irreversibly randomised; precise timestamps recovered from
IDs, max error ≤ 1 s*  
* **A fitted weekly rhythm–modulation function** derived from that corpus  
* **A CMA-ES → PSO → LSQ pipeline written in MATLAB R2024b** that calibrates a four-compartment
RhP–SEIRS model to real cascades  
* **Best-fit parameters, diagnostic plots and global-sensitivity analyses**

The complete data–processing chain is documented in the paper; only essential usage
instructions are repeated below.

---

## 2 Directory layout

```text
best_parameter/                     # best-fit parameters for each city (MAT)
    Chengdu.mat
    Kashgar.mat
    Qingdao.mat
cascades/                           # benchmark cascades (timestamps only)
    Chengdu_out.csv
    Kashgar_out.csv
    Qingdao_out.csv
CMA-PSO-LSQ_pipeline/               # end-to-end fitting pipeline
    Chengdu_ode_pipeline_wtNLL.m
    Kashgar_ode_pipeline_wtNLL.m
    Qingdao_ode_pipeline_NLL.m
ode_check/
    ode_final_check_chengdu.m
    ode_final_check_kashgar.m
    ode_final_check_qingdao.m
rhythm/
    preprocess-1-merge.py
    preprocess-2-interpolation.py
    preprocess-3-final.py
    raw_posts_idp/                  # the original data after slicing
    preprocessed_data/              # preprocessed data and spline parameters 
sensitivity_analysis/
    ode_check_final_prcc.m
    ode_check_final_S1_ST.m
    ode_check_final_sobol_dynamic.m
    sample/                         # analysis sample 
```

---

## 3 Quick start — reproduce Qingdao in one line

```bash
# from the repository root
matlab -batch "addpath(genpath('.')); CMA_PSO_LSQ_pipeline/Qingdao_ode_pipeline_wtNLL"
```

This will:
* load cascades/Qingdao_out.csv,
* read the rhythm baseline from rhythm/,
* run CMA-ES → PSO → LSQ,
* save the fitted θ to best_parameter/, and
* plot the predicted vs. observed cumulative cascade.

Run the analogous script for Kashgar or Chengdu to replicate all results in the
paper.


---

## 4 Requirements

software	version
MATLAB	R2024b (earlier 2023a-2024a work but not tested exhaustively)
optional Python	≥ 3.9, only needed if you wish to regenerate rhythm/ from the raw crawl; pip install pandas pyarrow tqdm cryptography

All MATLAB scripts are self-contained (no toolboxes beyond Statistics and Optimization).

---

## 5 Data provenance & privacy
* posts_uidp_part\*.parquet contains timestamps and numerical metadata only.
Original text, screen names, profile URLs etc. have been permanently removed.
* The column "weibo_id_rand" is a cryptographically random, one-way mapping of the original Weibo IDs.
The secret key will never be disclosed; re-identification is infeasible.
* The three cascade files originate from the MIT-licensed Weibo-COV corpus but retain
timestamps only.
Obtain the full corpus from the original authors if you need text or user fields.

This design satisfies the “irreversible anonymisation” requirement of China PIPL and the
“pseudonymisation” guidance (EDPB 01/2025).

---

## 6 How to cite

* Code DOI  : 10.5281/zenodo.15378206   (this repository)
* Data DOI  : 10.5281/zenodo.15378372   (all versions)

---

## 7 Licence
* **Source code** — MIT License (see 'LICENSE')
* **Dataset (rhythm/raw_posts_idp/ & rhythm/preprocessed_data/)** — CC BY-NC 4.0 (see 'DATA_LICENSE')

---

## 8 Contact

Questions, bug reports or suggestions → Yushi Sun
sun-yushi@outlook.com

