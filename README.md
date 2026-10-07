# gbmMINER replication (GBM)

Replication of the gbmMINER transcriptional regulatory network for glioblastoma (Turkarslan et al., medRxiv 2024) on TCGA-GBM.
Notebooks 02 to 07 follow the pan-cancer pipeline (to be uploaded), with additions for GBM: pathway mutations, miRNA
regulators and the Gravendeel validation cohort. Notebook 01 is specific to GBM because it combines RNA-seq and microarray
patients, as the paper does.

**Preprocessing (notebooks 01 and 01b)**
RNA-seq data from the GDC are combined with Affymetrix microarray data to give 550 patients (284 with RNA-seq, 266 with
microarray only). Expression is z-scored within each platform to put the two on a common scale, and the 589 genes whose
RNA-seq and microarray measurements correlate negatively are dropped, leaving 9,790. Notebook 01b builds the miRNA
expression matrix.

**Clustering and causal mapping (notebook 02)**
MINER clusters genes on the basis of co-expression across patients. Within a cluster, the genes that carry a binding site
for the same transcription factor, and whose expression correlates with that TF, form a regulon. Somatic mutations in genes
and pathways are then linked to regulons through their regulators (causal flows); miRNA regulators are tested separately.
Each regulon is scored in every patient as overactive, neutral or underactive. Patients with similar activity profiles
across the network are grouped into transcriptional states, which may correspond to disease subtypes, and regulons whose
activity varies together across patients are grouped into transcriptional programs.

**Clinical data (notebook 03)**
Overall survival (time and status) is taken from the TCGA clinical file on cBioPortal; patients without survival
information are dropped (594 patients, 492 deaths).

**Disease relevance (notebook 04)**
Each regulon is tested for survival association (Cox regression on its eigengene), for enrichment of cancer hallmarks
(Gene Ontology terms, Lin similarity >= 0.8) and for replication in an independent cohort (Gravendeel, GSE16011).
The paper's rule for a disease-relevant regulon, as applied in its Supplementary Table 1, has four parts:

1) Co-expression in TCGA: variance explained by PC1 >= 0.3 and permutation p < 0.05
2) Survival association in TCGA: Cox hazard ratio p <= 0.05
3) Replication of co-expression in the independent cohort: PC1 >= 0.3 and permutation p <= 0.05
4) Replication of the survival association in the independent cohort: Cox p <= 0.05, same direction of effect

Parts 2 to 4 are applied by the notebook and give 1,147 disease-relevant regulons (the paper: 1,083). Part 1 is not a
filter in the notebook; checked afterwards, all 1,147 pass it (3,872 of all 3,961 regulons do). Notebooks 05 to 07 use an
earlier, looser list of 1,292 regulons: (survival-associated or hallmark-enriched) and co-expression validated.

**Risk models (notebook 05)**
Three ridge regression models predict a patient's survival risk (GuanRank score) from gene expression, regulon activity or
program activity. Each is trained on 70% of the TCGA patients, tuned on 10% and tested on the remaining 20%, and then
tested without refitting on the 155 GBM patients of Gravendeel. Performance is reported as the C-index and as Kaplan-Meier
curves of the predicted high- and low-risk groups.

**Master regulators (notebook 06)**
Causal flows are restricted to published driver genes, and transcription factors that are barely expressed are dropped.
A TF-TF network is built with LASSO regression: each TF's expression is predicted from the other TFs that have a binding
site for it. For each transcriptional state, the TFs linked to that state are ranked by the ratio of out-degree to
in-degree within the network, and those at the top are called master regulators. The paper fits a separate network for
each state; here one network is fitted for all TFs, and the result does not replicate the paper's (see `RESULTS.md`).

**Candidate regulators (notebook 07)**
An extension that is not in the paper: TFs are ranked by combining several lines of evidence with Robust Rank Aggregation.
Not yet run for this dataset.

**Results**
`RESULTS.md` gives the numbers, compares them with the paper, and states what does not replicate and the known issue with
the transcriptional states, half of which follow the measurement platform. The main result tables are in `results/`.

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


