# Result tables

Copies of the smaller result files of the combined-platform run (6 Oct 2026). The numbers in `../RESULTS.md` come from these files.
Regulon IDs are those of this run; genes and regulators are Ensembl IDs unless the file uses symbols.

| File | Written by | Contents |
|---|---|---|
| `regulons.json`, `regulonDf.csv` | 02 | the 3,961 regulons: member genes, and one row per (regulon, TF, gene) |
| `transcriptional_states.json` | 02 | patients of each of the 24 states |
| `transcriptional_programs.json` | 02 | regulons of each of the 193 programs |
| `programs_vs_states.csv` | 02 | program activity in each state |
| `state_coherence_check.csv` | 02 | within-state similarity of each state against random patient sets |
| `mutation_filter_log.csv`, `pathway_mutation_log.csv` | 02 | mutated genes and pathways considered, with frequencies |
| `filteredCausalResults.csv` | 02 | causal flows (mutation → regulator → regulon) passing all four criteria |
| `mirna_regulonDf.csv`, `mirna_regulon_regulators.csv`, `mirna_filteredCausalResults.csv` | 02 | miRNA regulons, their regulators and causal flows |
| `regulon_cox_results.csv` | 04 | Cox hazard ratio and p-value per regulon |
| `disease_relevant_regulons_survival_or_hallmark.csv` | 04 | regulon IDs, survival-associated or hallmark-enriched (no header) |
| `gravendeel_regulon_coherence.csv`, `gravendeel_coexp_validated_regulons.csv` | 04 | co-expression validation in Gravendeel |
| `final_disease_relevant_regulons.csv` | 04 | the 1,292 disease-relevant regulon IDs (no header) |
| `gravendeel_regulon_cox_results.csv` | 04 | Cox hazard ratio and p-value per regulon in Gravendeel (grades II to IV) |
| `final_disease_relevant_regulons_paper_rule.csv` | 04 | the 1,147 regulon IDs that are disease-relevant under the paper's rule (no header) |
| `program_cox_results.csv` | 04 or 05 | Cox hazard ratio and p-value per program (program activity -1 / 0 / +1) |
| `risk_model_summary.csv`, `tcga_cv_summary.csv`, `tcga_repeated_splits.csv` | 05 | risk model settings and C-indices |
| `driverCausalResults.csv` | 06 | causal flows kept by the driver filter |
| `tf_tf_network.csv`, `tf_network_degrees.csv` | 06 | TF–TF network edges and each TF's in- and out-degree |
| `master_regulators_per_state.csv` | 06 | master regulators of each state |

Not included because of size: regulon eigengenes and the over/under-activity matrices (8 to 42 MB each).
