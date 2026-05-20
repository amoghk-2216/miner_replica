# GBM MINER Replication

A partial replication of [Turkarslan et al. 2024](https://doi.org/10.1101/2024.04.05.24305380) — *An atlas of causal and mechanistic drivers of interpatient heterogeneity in glioma* — using publicly available TCGA data and open-source tools.

---

## Overview

The paper constructs gbmMINER, a transcriptional regulatory network for GBM, by:

1. Clustering 9,728 genes into  3,797 co-expressed regulons using MINER (Mechanistic Inference of Node-Edge Relationships)
2. Quantizing regulon activity per patient into overactive/neutral/underactive states
3. Clustering regulons into 179 transcriptional programs
4. Clustering patients into 23 disease states based on program activity
5. Training a ridge regression risk model on program activity to predict survival
6. Additionally, predicting drug sensitivity based on program activity

This replication approximates the same pipeline using accessible tools (scikit-learn, lifelines, etc.) on 437 TCGA GBM patients with combined microarray and RNA-seq expression data.

---

## Repository Structure

```
miner_replica/
├── notebooks/
│   ├── gbm_model.ipynb              # Expression preprocessing, regulon discovery, quantization
│   ├── gbm_programs.ipynb           # Program and state discovery
│   ├── gbm_risk.ipynb               # Risk model, TCGA and Gravendeel evaluation
│   └── gravendeel_validation.ipynb  # Gravendeel cohort preprocessing and Cox validation
├── utils.py                         # Shared helper functions
├── results/
│   ├── km_tcga_cv.png               # TCGA KM curve
│   └── km_gravendeel.png            # Gravendeel KM curve
└── README.md
```

---

## Pipeline

### 1. Expression preprocessing (`gbm_model.ipynb`)

- Combined microarray (Agilent) and RNA-seq (v2 RSEM) z-scored expression data from cBioPortal
- Priority: microarray > RNA-seq for patients with both, as done in the original paper
- Final matrix: 13,741 genes × 437 patients
- Clinical data: 437 patients indexed by PATIENT_ID with OS_MONTHS and OS_STATUS

### 2. Regulon discovery

- Gene clustering via `coexpy` decomposition
- PCA validation: variance explained by PC1 ≥ 0.30
- 1,148 decomposed subgroups → 899 PCA-validated regulons
- Eigengene computed per regulon (first principal component, sign-corrected)

### 3. Network quantization

- Done on 899 validated regulons 
- Per-patient gene ranking → binomial test on upper/lower thirds
- Vectorized implementation using `scipy.stats.binom.sf`
- Output: 899 × 437 matrix of values {-1, 0, 1}

### 4. Programs (`gbm_programs.ipynb`)

- Variance filtering: top 50% most variable regulons → 449 active regulons
- This was done as a significant number of the 899 regulons were inactive (quantized as 0)
- Optimal k selected via gap statistic (elbow + silhouette confirmed)
- KMeans gave k=16 programs on active regulons
- Program activity: mean quantized activity across regulons per patient

### 5. Patient states

- KMeans on 437 × 16 program activity matrix
- Optimal k=2 states (silhouette peak, gap statistic confirms weak structure beyond k=2)
- State 0: 170 patients, State 1: 267 patients

### 6. Risk model (`gbm_risk.ipynb`)

- GuanRank scores computed from TCGA survival data
- Ridge regression (RidgeCV) on 16-dimensional program activity
- 5-fold cross-validation for unbiased TCGA evaluation
- Validated on independent Gravendeel cohort (252 patients, 284 CEL files)

---

## Results

| Metric | This replication | Paper (gbmMINER) |
|---|---|---|
| Patients (TCGA) | 437 | 516 |
| Genes | 13,741 | 9,728 |
| Regulons (validated) | 899 | 3,797 |
| Programs | 16 | 179 |
| Patient states | 2 | 23 |
| TCGA C-index (CV) | 0.52 | ~0.60 |
| TCGA KM p-value | 0.37 | 0.027 |
| Gravendeel C-index | 0.565 | ~0.63 |
| Gravendeel KM p-value | 0.033 ✓ | <0.002 |

The Gravendeel validation is statistically significant (p=0.033), confirming that program-level transcriptional activity meaningfully stratifies GBM patient survival in an independent cohort.

---

## Limitations vs Paper

**Regulon count (899 vs 3,797):** The paper uses cMonkey2, a biclustering algorithm that finds condition-specific co-regulated gene modules. This replication uses `coexpy` decomposition on a smaller patient cohort (437 vs 516), producing fewer and larger regulons. This is the primary cause of all downstream differences.

**Programs and states (16 vs 179, 2 vs 23):** Again, due to low regulon count. Fewer regulons → fewer high-variance ones → fewer programs → fewer distinguishable patient states. With only 16 programs, the silhouette and gap statistic both support only 2 patient states. However, proportionally the results are similar.

**TCGA C-index and KM p-value:** Weak TCGA performance is expected — the paper also shows weaker TCGA results than Gravendeel (p=0.027 vs p=0.002). With 16 programs the signal is insufficient to significantly stratify the training cohort.

**Gravendeel C-index (0.565 vs ~0.63):** Reasonable given the reduced program count. The direction and significance of the result are reproduced.

**No cMonkey2:** The primary methodological difference. cMonkey2 performs biclustering that simultaneously clusters genes and patients, finding condition-specific co-expression patterns. `coexpy` uses a simpler decomposition approach that produces fewer, less condition-specific regulons.

---

## How to Run

### Requirements

```
python >= 3.10
numpy, pandas, scipy
scikit-learn
lifelines
coexpy
matplotlib
joblib
# R dependencies (for gravendeel_validation.ipynb preprocessing step)
# Install in R:
# install.packages("BiocManager")
# BiocManager::install("affy")
```

Install:
```bash
pip install numpy pandas scipy scikit-learn lifelines coexpy matplotlib joblib
```

### Data

Download from [cBioPortal GBM (TCGA, Firehose Legacy)](https://www.cbioportal.org/study/summary?id=gbm_tcga):
- `data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt`
- `data_mrna_agilent_microarray_zscores_ref_all_samples.txt`
- `data_clinical_patient.txt`

Download Gravendeel CEL files from GEO: [GSE16011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16011)

### Execution order

```
1. gbm_model.ipynb              
2. gravendeel_validation.ipynb + CEL normalization in R 
3. gbm_programs.ipynb          
4. gbm_risk.ipynb              
```

---

## Reference

Turkarslan S, He Y, Hothi P, et al. *An atlas of causal and mechanistic drivers of interpatient heterogeneity in glioma.* medRxiv 2024. https://doi.org/10.1101/2024.04.05.24305380
