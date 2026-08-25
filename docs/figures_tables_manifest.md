# Manuscript tables and figures — approved working set

## Main manuscript

### Table 1
**File:** `results/tables/MANUSCRIPT_Table1_dataset_characteristics.csv`

**Proposed title:** Characteristics of the human TRPA1 inhibitory-activity
dataset derived from ChEMBL 37.

### Table 2
**File:** `results/tables/MANUSCRIPT_Table2_scaffold_benchmark.csv`

**Proposed title:** Performance of molecular representations and regressors
under three scaffold-aware cross-validation partitions.

**Note:** report mean ± SD across the three partitions. Do not also include
the candidate benchmark figure in the final manuscript.

### Figure 1
**File:** `results/figures/MANUSCRIPT_Fig1_workflow.png`

**Proposed caption:** Study workflow. Exact IC50 records for human TRPA1
were obtained from ChEMBL 37, standardized and aggregated by compound.
Four molecular representations were evaluated with two regressors under
scaffold-aware cross-validation. Out-of-fold predictions were then used
to examine the association of prediction error with disagreement among
published potency values, chemical similarity to the training set and the
validation strategy.

### Figure 2
**File:** `results/figures/MANUSCRIPT_Fig2_variability_vs_error.png`

**Proposed caption:** Association between dispersion of published pIC50
values for the same compound and absolute out-of-fold prediction error.
Points show Spearman correlation coefficients; horizontal intervals show
95% scaffold-cluster bootstrap confidence intervals. Analyses were
performed separately for compounds represented in at least two assays and
in at least two source documents.

### Figure 3
**File:** `results/figures/MANUSCRIPT_Fig3_similarity_vs_error.png`

**Proposed caption:** Association between maximum Morgan/Tanimoto similarity
to training compounds and absolute out-of-fold prediction error. Points
show Spearman correlation coefficients and horizontal intervals show 95%
scaffold-cluster bootstrap confidence intervals.

### Figure 4
**File:** `results/figures/MANUSCRIPT_Fig4_random_vs_scaffold_R2.png`

**Proposed caption:** Comparison of random and scaffold-aware validation.
Points show mean R² across three partitions and error bars show the
standard deviation. Random splitting produced systematically higher
apparent predictive performance.

## Candidate only

`results/figures/CANDIDATE_Fig2_scaffold_benchmark_R2.png`

This is useful for the supervisor draft but should not be submitted
together with Table 2 because it duplicates the same benchmark result.

## Internal statistics, not a manuscript table

`results/tables/MANUSCRIPT_internal_exact_statistics.csv`

Use it to write exact values and p-values; do not automatically insert it
as a third table.
