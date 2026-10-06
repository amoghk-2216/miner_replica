# Results: combined-platform run (6 Oct 2026)

Run settings: RNA-seq and microarray patients z-scored within platform and combined, as in the paper;
`MIN_CROSS_PLATFORM_R = 0` (genes whose two platform measurements correlate negatively are dropped);
pathway mutations and miRNA regulators included; all genes mutated in at least 5% of patients tested (`RESTRICT_TO_DRIVERS = False`).

Notebooks 01 to 06 have been run. **Notebook 07 has not been run yet for this run**; its section is a placeholder.

## Reproducibility note

Notebook 02 was run twice on the same inputs. The outputs displayed in the notebook are from the second run; every
number in this file and everything in notebooks 03 to 06 comes from the first. The two runs have identical co-expression
clusters and identical transcriptional states (the same patients in each of the 24). The second run has one regulon
fewer (3,960 against 3,961; a six-gene regulon is absent), 194 programs against 193 (186 identical) and 26,986 causal
flows against 26,988. Because regulons are numbered in order, the regulon IDs above 3590 differ by one between the runs.
The pipeline is therefore reproducible to within one regulon, not exactly.

## Network size, against the paper

The paper's counts are from its Results section, Table 1 and Supplementary Tables 1, 5 and 9. Rows marked "disease-relevant"
count only what acts on disease-relevant regulons, which is how the paper reports them; they use the 1,147 regulons that are
disease-relevant under the paper's rule.

| | This run | Paper | Note |
|---|---|---|---|
| Patients | 550 (284 RNA-seq, 266 microarray only) | 516 | all 516 of the paper's patients are among the 550; the paper's Methods give 536 tumours + 4 post-mortem controls |
| Genes entering MINER | 9,790 | 10,263 | |
| Regulons | 3,961 | 3,797 | |
| Mean genes per regulon | 11.0 | 11 | |
| Genes in any regulon | 6,701 | 9,728 | not explained; see open questions |
| Transcriptional states | 24 | 23 | |
| Transcriptional programs | 193 | 179 | |
| Disease-relevant regulons, paper's rule | 1,147 | 1,083 | survival-associated in TCGA and Gravendeel, co-expression validated |
| Disease-relevant regulons, looser rule used by notebooks 05 to 07 | 1,292 | n/a | see next section |
| Survival-associated programs (Cox p ≤ 0.05 on program activity) | 35 of 193 (25 high-risk, 10 low-risk) | 58 of 179 (51 high-risk, 7 low-risk) | about 10 of 193 expected by chance |
| TFs regulating a disease-relevant regulon | 333 | 306 | |
| Genes in disease-relevant regulons | 1,954 | 2,049 | |
| Mutated genes with a causal flow to a disease-relevant regulon | 23 of 44 tested | 1,017 in its Supplementary Table 2 | the paper's "30" is its curated list of prognostic markers (Supplementary Table 5); see below |
| Mutated pathways with a causal flow to a disease-relevant regulon | 145 | 49 | cannot be traced: the paper's Supplementary Table 2 lists no pathway flows |
| miRNA regulators | 14 | 73 | different method (miner target sets, not FIRM) |

Causal flows (mutation → regulator → regulon) passing the paper's four criteria: 26,988 over all regulons, of which 9,072
act on a disease-relevant regulon. The mutated genes with the most flows are ATRX, IDH1, TP53 and NF1, as in the paper.
**Causal flows against the paper's Supplementary Table 2.** The filter is the same: the smallest values that pass are identical
in both (mutation → regulator p ≤ 0.05, regulator → regulon Spearman p ≤ 0.05, regulon stratification |t| of about 2, at least
half of the edges aligned). What differs is how many mutations are tested:

- The paper's table has 16,382 flows from 1,017 mutated genes onto 1,074 regulons. Here 44 genes pass the 5% mutation filter
  in the MC3 calls, and 23 of them have a flow onto a disease-relevant regulon. The paper's Methods state the same 5% rule, so its
  mutation calls must differ from MC3 (its list includes genes such as KRTAP4-11 and HSD17B7P2 that do not reach 5% here).
- The genes with the most flows agree: ATRX, NF1, IDH1 and TP53 in the paper; ATRX, IDH1, TP53 and NF1 here.
- Of the 30 prognostic markers the paper highlights (its Supplementary Table 5), 27 pass the mutation filter here and 14 have a
  flow onto a disease-relevant regulon: AHNAK2, APOB, ATRX, EGFR, IDH1, KEL, NF1, PIK3CA, PTEN, RB1, SPTA1, SYNE1, TP53, TTN.
- All flows in the paper's table are from mutated genes. Here pathways give 22,352 of the 26,988 flows, so the gene-level flows
  (4,636) are the part that can be compared with the paper.

## Disease relevance (notebook 04)

| Step | Regulons |
|---|---|
| All | 3,961 |
| Survival-associated (Cox p ≤ 0.05) | 1,354 |
| Survival-associated or hallmark-enriched | 1,536 |
| Co-expression validated in Gravendeel | 2,379 (60%) |
| Final: (survival or hallmark) and validated | 1,292 |

The Gravendeel cohort used for co-expression validation contains gliomas of all grades (276 samples).

**This rule is looser than the paper's.** In the paper's Supplementary Table 1 every one of the 1,083 disease-relevant
regulons has a significant Cox p-value in TCGA and a significant Cox p-value in at least one independent cohort
(1,071 of them in Gravendeel), in addition to validated co-expression. Here:

- survival association is not required to replicate in Gravendeel;
- 138 of the 1,292 regulons enter through hallmark enrichment alone (1,154 are survival-associated and validated);
- the paper validates in any of three microarray cohorts (Gravendeel, Rembrandt, GSE7696) and uses grade II to IV gliomas
  (252 Gravendeel patients); this run uses Gravendeel only and all 276 glioma samples, including 8 of grade I.

**Count under the paper's rule** (last cell of notebook 04; Cox regression of survival on each regulon's eigengene in the
264 Gravendeel patients of grades II to IV, 239 deaths):

| Step | Regulons |
|---|---|
| Survival-associated in TCGA (Cox p ≤ 0.05) | 1,354 |
| and co-expression validated in Gravendeel | 1,154 |
| and survival-associated in Gravendeel (Cox p ≤ 0.05) | 1,148 |
| and with the same direction of effect in both cohorts | **1,147** |

So 1,147 regulons are disease-relevant under the paper's rule, against the paper's 1,083. They are listed in
`final_disease_relevant_regulons_paper_rule.csv`. Notebooks 05 to 07 still use the looser list of 1,292.

Two points on reading this. The survival replication is a weak filter in this cohort: 85% of all 3,961 regulons are
survival-associated in Gravendeel, because it mixes glioma grades with very different survival. The direction of effect is
the more informative check: it agrees with TCGA for 99% of the regulons that are survival-associated in TCGA, against 78%
of those that are not.

## Risk models (notebook 05)

Ridge regression on GuanRank scores; features and penalty chosen on a validation split. The Gravendeel test uses the
155 GBM (grade IV) patients with survival data (147 deaths), as the paper's test set is GBM only.

| Model | Features used | TCGA test C-index, 10 splits (mean ± sd) | TCGA 5-fold C-index | Gravendeel GBM C-index | Gravendeel log-rank p |
|---|---|---|---|---|---|
| Gene | 25 | 0.558 ± 0.035 | 0.575 | 0.613 | 2.1e-07 |
| Regulon | 2,000 | 0.544 ± 0.049 | 0.565 | 0.621 | 1.8e-04 |
| Program | 193 | 0.532 ± 0.042 | 0.527 | 0.562 | 1.6e-04 |

- On TCGA the three models are within one standard deviation of each other.
- On Gravendeel the gene and regulon models are level and the program model is lower.
- The paper's Supplementary Table 8 shows the same ordering on these two cohorts. Rebuilding its Kaplan-Meier data gives
  log-rank p = 0.0014 (gene), 0.0081 (regulon) and 0.027 (program) on the TCGA test set, and < 0.0001, 0.0014 and 0.0021
  on Gravendeel. The paper's program model is the strongest only in its third cohort (XCELSIOR), which is not public.
- A diagnostic with continuous instead of three-level activity gave 0.604 (regulon) and 0.572 (program) on Gravendeel,
  so the three-level quantisation is not what limits the program model.
- Differences from the paper's setup: its gene model uses TPM values (here z-scores), it reports one 70/10/20 split
  (here ten, plus 5-fold cross-validation), and it tests on 113 TCGA and 150 Gravendeel patients.

## Known issue: half of the transcriptional states follow platform

- 12 of the 24 states are at least 80% one platform (274 of 550 patients); three are 100% one platform.
- The two platform groups do not differ in age, survival, Karnofsky score or mutation frequencies (IDH1 6.5% vs 6.7%).
  They differ in year of diagnosis and tissue source site, so the split is technical.
- IDH1-mutant patients fall into two states, one per platform (state 9: 4 of 6 mutated; state 19: 8 of 20 mutated;
  7% overall).
- The mean level of every regulon is the same on both platforms. The signal is in the fine structure: for 90% of
  patients the most similar other patient is on the same platform (about 50% expected).
- A stricter gene filter reduces it only at a large cost in genes (r ≥ 0.7: 1,230 genes, 63%).
- Normalising samples does not remove it: after quantile-normalising every sample to the same distribution and z-scoring the
  genes again, the nearest neighbour is still on the same platform for 94% of patients (the current matrix: 94%). miner's own
  preprocessing function cannot be used for this test; its quantile route fails on input that is already z-scored.
- Expression subtype (the paper's Supplementary Table 9) confirms it. The two platform groups have a similar subtype mix
  (microarray vs RNA-seq: classical 79 vs 61, mesenchymal 73 vs 82, proneural 50 vs 48, neural 35 vs 50, G-CIMP 22 vs 16).
  The states do follow subtype (states 0 to 5 are mostly mesenchymal, 13 proneural, 18 neural), but one subtype is split
  into separate states by platform: classical into state 15 (94% classical, 94% microarray) and state 14 (68% classical,
  81% RNA-seq); G-CIMP into state 19 (16 G-CIMP patients, 79% microarray) and state 9 (4 of 6 G-CIMP, all RNA-seq).

What this affects: the states, the programs (built using the state clustering), the per-state master regulators of
notebook 06 and, through them, one of the five evidence lines of notebook 07. Regulons, causal flows, disease relevance
and the gene and regulon risk models do not use the states.

The paper combines the platforms in the same way and does not report the platform composition of its states.
A microarray-only run (528 patients) is planned as a sensitivity check.

## Open questions

- Regulons here cover 6,701 genes; the paper's cover 9,728, although the number of regulons and their mean size agree.
- Notebook 02 uses miner's default `min_correlation = 0.2` (the value in miner's example analysis). The paper's Methods give a
  minimum regulator-to-regulon Spearman correlation of 0.3 as a later filtering criterion, which is not applied here.
- The paper tests far more mutated genes (1,017 with a causal flow) than pass the 5% filter here (44); its mutation calls are not MC3.
- 145 mutated pathways have a causal flow to a disease-relevant regulon, against 49 in the paper, whose supplementary tables list no pathway flows.
- The 5% mutation filter keeps five genes mutated in 17 of 350 patients (4.9%); the cut-off should be checked.
- The paper's training data include four post-mortem control samples; this run has none.

## Master regulators (notebook 06)

- Driver filter: 26,231 of the 26,988 causal flows are kept. They start from 9 driver genes (ATRX, EGFR, IDH1, NF1, PIK3CA, PIK3R1, PTEN, RB1, TP53) and 172 pathways;
  pathways are not subject to the driver filter.
- TF expression filter: 279 of the 294 causal regulators are expressed in enough tumours and enter the TF-TF network
  (14,023 edges).
- 204 of the 294 causal regulators are among the 306 disease-relevant TFs of the paper (its Supplementary Table 3).
- 74 TFs are a master regulator of at least one state; 51 of them are among the paper's 306 TFs.
  Every state has between 11 and 17 master regulators.

**This result does not replicate the paper's.** The paper reports 67 master regulators covering 20 of its 23 states and
states that most are specific to one state, with only a few shared (it names MAFB, RUNX1, STAT4, SOX2, SOX9 and E2F1).
Here KLF15 and SP1 are master regulators of all 24 states, TCF3 of 23, FOXC2 of 20, NR3C1 of 19, MAZ of 18 and MSX2 of 16;
31 of the 74 appear in one state only. None of the six TFs the paper names is a master regulator here (RUNX1, STAT4, SOX9
and E2F1 are in the network; MAFB and SOX2 are not among the causal regulators).

The likely cause is how notebook 06 differs from the paper's Methods:

- the paper fits the LASSO models separately for each state, on that state's list of TFs; the notebook fits them once
  over all 279 TFs and takes each state's subset of the edges;
- the notebook has no pruning step;
- the notebook's list of TFs for a state is barely selective: regulons are taken at an uncorrected p ≤ 0.05, which gives
  a median of 260 of the 279 TFs per state (range 182 to 274), so every state sees almost the whole network;
- the notebook adds a selection step (miner's `unmix` and `select_optimal_tfs`) before ranking by out-degree / in-degree.

The paper does not define "causally associated with the state" and its script is not linked in the preprint, so an exact
reimplementation is not possible from the text alone.

Platform share of each state, to be read with the known issue above (states at or above 80% are marked):

| State | Patients | Largest platform share | Master regulators | First five |
|---|---|---|---|---|
| 0 | 38 | 68% RNA-seq | 14 | FOXC2, MAZ, SP1, TCF3, KLF15 |
| 1 | 37 | 65% microarray | 15 | MSX2, FOXC2, MAZ, KLF15, SP1 |
| 2 | 20 | 100% microarray **(skewed)** | 17 | MSX2, FOXC2, MAZ, NKX2-5, TCF3 |
| 3 | 20 | 75% microarray | 16 | MSX2, MAZ, TCF3, KLF15, SP1 |
| 4 | 18 | 56% microarray | 13 | MSX2, FOXC2, SP1, TCF3, KLF15 |
| 5 | 8 | 50% microarray | 14 | MSX2, FOXC2, MAZ, TCF3, KLF15 |
| 6 | 53 | 92% RNA-seq **(skewed)** | 17 | MSX2, FOXC2, MAZ, SP1, TCF3 |
| 7 | 21 | 100% RNA-seq **(skewed)** | 16 | FOXC2, MAZ, KLF15, TCF3, SP1 |
| 8 | 10 | 90% microarray **(skewed)** | 17 | MSX2, FOXC2, MAZ, TCF3, SP1 |
| 9 | 6 | 100% RNA-seq **(skewed)** | 16 | MSX2, FOXC2, MAZ, SP1, TCF3 |
| 10 | 33 | 88% RNA-seq **(skewed)** | 14 | MSX2, MAZ, SP1, TCF3, HSF1 |
| 11 | 14 | 86% microarray **(skewed)** | 13 | ZNF35, FOXC2, MAZ, TEAD3, SP1 |
| 12 | 13 | 62% microarray | 14 | MAZ, FOXC2, TEAD3, SP1, TCF3 |
| 13 | 31 | 52% microarray | 16 | ZNF35, MAZ, FOXC2, TEAD3, SP1 |
| 14 | 48 | 81% RNA-seq **(skewed)** | 13 | FOXC2, MAZ, TCF3, SP1, TEAD3 |
| 15 | 18 | 94% microarray **(skewed)** | 16 | FOXC2, NKX2-5, SP1, TCF3, KLF15 |
| 16 | 28 | 93% RNA-seq **(skewed)** | 14 | MSX2, FOXC2, TCF3, SP1, KLF15 |
| 17 | 25 | 72% microarray | 14 | MSX2, FOXC2, MAZ, TEAD3, SP1 |
| 18 | 15 | 73% microarray | 17 | MSX2, FOXC2, TCF3, SP1, KLF15 |
| 19 | 47 | 79% microarray | 14 | FOXC2, MAZ, KLF15, TCF3, SP1 |
| 20 | 11 | 91% microarray **(skewed)** | 11 | MSX2, MAZ, KLF11, NKX2-5, SP1 |
| 21 | 12 | 83% microarray **(skewed)** | 16 | MSX2, FOXC2, TCF3, SP1, NR3C1 |
| 22 | 15 | 53% microarray | 14 | ZNF263, MSX2, FOXC2, MAZ, TCF3 |
| 23 | 9 | 67% microarray | 11 | MSX2, TFAP2A, SP1, EGR1, KLF15 |

`NF1` and `VBP1` appear in the list because miner's binding-site database has entries under those names; NF1 here
is most likely the nuclear factor I motif family, not the neurofibromin gene. Check before reporting either as a TF.

## Candidate regulators (notebook 07)

Not yet run for this run. This notebook is an extension and is not part of the paper.
