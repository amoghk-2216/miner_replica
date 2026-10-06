# gbmMINER replication (GBM)

Replication of the gbmMINER transcriptional regulatory network for glioblastoma (Turkarslan et al., medRxiv 2024) on TCGA-GBM.
Notebooks 02 to 07 are the same as in the pan-cancer pipeline (to be uploaded); notebook 01 is specific to GBM because it combines RNA-seq and
microarray patients, as the paper does.

## Layout

```
config.py                 cancer, data root and per-cancer settings to customise the pipeline
driver_tiers.py           published driver genes per cancer (used by notebooks 02, 06, 07)
notebooks/                01 to 07, run in order
data/gbm/gbm_data/        inputs: gbm_tcga/ (cBioPortal), rnaseq_gdc/, mirna_gdc/, gbm_gravendeel/ and the files written by 01 and 02
data/gbm/gbm_output/      everything the notebooks write from 02 onwards
reference/                copy of miner's network dictionaries
go-basic.obo, goa_human.gaf   Gene Ontology files (also read by the pan-cancer project, so they stay at the top level)
RESULTS.md                numbers from the current run, how they compare with the paper, and what does not replicate
results/                  the smaller result tables of the current run (see results/README.md)
```

The `data/` folder (inputs and all files the notebooks write), `reference/` and the Gene Ontology files are not uploaded to this repository.
The results are summarised in `RESULTS.md`, the main result tables are in `results/`, and the notebooks show their outputs.

## Notebooks

| # | Notebook | What it does |
|---|---|---|
| 01 | `01_gbm_expression` | GDC RNA-seq + Affymetrix microarray, z-scored within platform, then combined; drops genes whose two platform measurements correlate below `MIN_CROSS_PLATFORM_R` (0 in the current run) |
| 01b | `01b_gbm_mirna` | miRNA expression matrix, used by the miRNA section of 02 |
| 02 | `02_regulons_clusters_states` | mutation matrix (MC3) with pathway mutations, clusters, regulons, causal flows, miRNA regulators, states, programs |
| 03 | `03_clinical` | survival table from the cBioPortal clinical file |
| 04 | `04_disease_relevance` | survival- or hallmark-associated regulons, co-expression validated in Gravendeel |
| 05 | `05_risk_model` | gene, regulon and program risk models; held-out TCGA test and Gravendeel |
| 06 | `06_regulators` | TF–TF network and master regulators per state |
| 07 | `07_candidate_regulators` | evidence-based ranking of transcription factors (an extension, not in the paper) |

Start Jupyter from `notebooks/`. The MC3 mutation file is read from the pan-cancer project folder by default (`MAF_PATH` in `config.py`).

## Differences from the paper

- Expression: genes whose microarray and RNA-seq measurements correlate negatively across the 262 patients measured on both are dropped (589 of 10,379); the paper applies no such filter.

- Each patient's expression is mean-centred after z-scoring (notebook 02); the paper describes no such step.

- Cohort: 550 patients; the paper's 516 are all among them. The paper's training data also include four post-mortem control samples, which are not used here.
- Mutations: genes and NCI-Nature (`PID_`) pathways, as in the paper.

- miRNA regulators come from the miRNA target sets shipped with miner (last section of notebook 02), not from the FIRM pipeline, and use only the patients with miRNA data.

- Validation: Gravendeel only; the paper accepts replication in any of three microarray cohorts (Gravendeel, Rembrandt, GSE7696). Gravendeel contains gliomas of grades I to IV: co-expression validation (04) uses all 276 samples, where the paper uses grades II to IV only; the survival test of the risk models (05) uses the 155 GBM (grade IV) patients.

- Disease-relevant regulons: here (survival-associated or hallmark-enriched) and co-expression validated. In the paper's Supplementary Table 1 every disease-relevant regulon is survival-associated in TCGA and again in an independent cohort, so the rule here is looser. Applying the paper's rule gives 1,147 regulons against its 1,083 (last cell of notebook 04); notebooks 05 to 07 use the looser list of 1,292.

- Hallmark association uses Lin similarity (the paper: Jiang-Conrath).

- Risk models use repeated 70/10/20 splits and cross-validation (the paper: one split), and the gene model uses z-scores (the paper: TPM). The paper's third test cohort (XCELSIOR) is not public and is not used.

- Notebook 06 restricts causal flows to published driver genes and drops barely expressed TFs. It also builds the TF–TF network once for all TFs, where the paper builds one per state and prunes it. **The master regulators do not replicate the paper's**: here the same few TFs lead in almost every state (see `RESULTS.md`).


