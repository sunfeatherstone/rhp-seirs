# RhP–SEIRS Weibo Dataset & Reproduction Pipeline  
*(companion artefact to the manuscript “Rhythm-modulated epidemic-like diffusion on social media”)*
---

## 1 Project overview
This repository contains **everything you need to reproduce every table, figure and
numerical result in the paper**.  
The core resources are

* **A nearly unbiased crawl of 8.7 million Weibo posts** (2019-2025)  
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


├── best_parameter/                   # Best-fit θ for each city (MAT files)
│   ├── Chengdu.mat
│   ├── Kashgar.mat
│   └── Qingdao.mat
├── cascades/                         # Three benchmark cascades extracted from
│   ├── Chengdu_out.csv               #   Weibo-COV; only the “created_at” column
│   ├── Kashgar_out.csv               #   is retained
│   └── Qingdao_out.csv
├── CMA-PSO-LSQ_pipeline/             # End-to-end fitting pipeline
│   ├── Chengdu_ode_pipeline_wtNLL.m
│   ├── Kashgar_ode_pipeline_wtNLL.m
│   └── Qingdao_ode_pipeline_NLL.m
├── ode_check/                        # Quick plotting scripts using the best θ
│   ├── ode_final_check_chengdu.m
│   ├── ode_final_check_kashgar.m
│   └── ode_final_check_qingdao.m
├── rhythm/                           # Raw & processed rhythm data
│   ├── posts_uidp.parquet            # 8.7 M posts (text stripped, IDs randomised)
│   ├── posts_uidp_interpolation_alignment.parquet
│   ├── posts_uidp_interpolation_alignment_bucket_count_2h.csv
│   └── posts_uidp_interpolation_alignment_spline_pp.mat
└── sensitivity_analysis/             # PRCC, Sobol, ST-index
├── ode_check_final_prcc.m
├── ode_check_final_S1_ST.m
└── ode_check_final_sobol_dynamic.m

---

## 3 Quick start — reproduce Chengdu in one line

```bash
# from the repository root
matlab -batch "addpath(genpath('.')); CMA_PSO_LSQ_pipeline/Chengdu_ode_pipeline_wtNLL"

This will
	1.	load cascades/Chengdu_out.csv,
	2.	read the rhythm baseline from rhythm/,
	3.	run CMA-ES → PSO → LSQ,
	4.	save the fitted θ to best_parameter/, and
	5.	plot the predicted vs. observed cumulative cascade.

Run the analogous script for Kashgar or Qingdao to replicate all results in the
paper.

⸻

4 Requirements

software	version
MATLAB	R2024b (earlier 2023a-2024a work but not tested exhaustively)
optional Python	≥ 3.9, only needed if you wish to regenerate rhythm/ from the raw crawl; pip install pandas pyarrow tqdm cryptography

All MATLAB scripts are self-contained (no toolboxes beyond Statistics and
Optimization).

⸻

5 Data provenance & privacy
	•	posts_uidp.parquet contains timestamps and numerical metadata only.
Original text, screen names, profile URLs etc. have been permanently removed.
	•	The column uid_rand is a cryptographically random, one-way mapping of the original
Weibo IDs.
The secret key will never be disclosed; re-identification is infeasible.
	•	The three cascade files originate from the MIT-licensed Weibo-COV corpus but retain
timestamps only.
Obtain the full corpus from the original authors if you need text or user fields.

This design satisfies the “irreversible anonymisation” requirement of China PIPL and the
“pseudonymisation” guidance (EDPB 01/2025).

⸻

6 How to cite

Data DOI: coming soon
Code DOI: coming soon

A full citation string will appear here once the Zenodo (code) and Mendeley Data
(dataset) records are registered.

⸻

7 Licence
	•	Data — Creative Commons CC BY-NC 4.0
(share alike, non-commercial, keep attribution)
	•	Code — MIT License

⸻

8 Contact

Questions, bug reports or suggestions → Yushi Sun
sun-yushi@outlook.com

⸻

Last updated: 2025-05-10

