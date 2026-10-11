# Methods fact sheet — попередній виконаний експеримент

**Редакція довідки: ORIGINAL-COMPARISON-2026-10-08.**  
**Поточну задачу й наступні дії визначає тільки [STATUS.md](../STATUS.md).**

Це технічна довідка про вже виконаний блок «чотири готові описи структури
× Random Forest/XGBoost» та його додаткові аналізи. Вона не є новим
протоколом змагання нейромереж і не підтверджує, що це змагання завершене.
Старі рукописи — документи попереднього етапу, не план поточної роботи.

Наведені параметри походять із попередніх перевірок коду та збережених
метаданих. Під час цієї передачі контексту навчання, відтворення результатів
або новий повний аудит не виконувалися. Перед використанням чисел у новому
експерименті звіряти конкретні вхідні файли й версії коду.
GitHub для AI — тільки читання.

## Data source and selection

- Frozen data: **ChEMBL 37**, target **CHEMBL6007, human TRPA1**.
- Endpoint: **IC50**, exact `standard_relation = "="`, pChEMBL required.
- All retained standard units in the frozen record-level table: **nM**.
- Record-level table: `data/raw/trpa1_current_api_raw.csv`.
- Molecule-level table: `data/processed/trpa1_primary_dataset.csv`.
- Extended source snapshot: `data/raw/trpa1_raw_snapshots/ChEMBL_37_20260806T144919Z/`.

Chronology documented in the previous manuscript and root `reproducibility_notes.md`:
the original molecule-level benchmark input was named `trpa1_antagonists.csv`
and attributed to ChEMBL 36; subsequent ChEMBL 37 comparison retained the same
1645 compounds and aggregate values. The later record-level analysis used
2196 measurements corresponding to these compounds. Do not confuse this
chronology with the aggregation mechanism or perform a new live extraction.

## Dataset size

- Activity records: **2196**.
- Standardized compounds: **1645**.
- Assays: **97**; ChEMBL documents: **55**; document years: **2010–2025**.
- Unique Bemis–Murcko scaffolds: **544**.
- One measurement: **1196** compounds; at least two measurements: **449**.
- At least two assays: **393** compounds; at least two documents: **52**.

Source: `results/tables/MANUSCRIPT_Table1_dataset_characteristics.csv`.

## Structure standardization and aggregation

Documented sequence: RDKit parsing; salt removal/main-fragment selection;
canonical isomeric SMILES and InChIKey; grouping by standardized structure;
median aggregation of reported pChEMBL values; Bemis–Murcko scaffolds.
Historical implementation evidence includes `scripts/resolve_delta.py` and
`scripts/chembl37_delta.py`. Do not run their live-API extraction as a rebuild.

The frozen table contains standardized SMILES, InChIKey, scaffold, pChEMBL
median/minimum/maximum/SD, measurement counts and other original columns.
**Still open:** one canonical offline `scripts/build_primary_dataset.py`
reproducing the frozen final table. Its absence does not mean that no historical
standardization code exists. Preserve row order for compatibility with stored
embeddings and predictions; report discrepancies rather than forcing a match.

## Prediction target

`pchembl_median` is the median of all retained IC50 pChEMBL records for a
standardized compound. It is not an assay-balanced median and not the IC50 of
one standardized protocol. The alternative median of compound–assay medians
is a proposed sensitivity analysis, not a completed change to the target.

## Molecular representations

1. **Morgan ECFP4:** RDKit Morgan generator, radius 2, 2048 bits,
   chirality not included.
2. **RDKit-15:** MolWt, MolLogP, MolMR, TPSA, NumHAcceptors, NumHDonors,
   NumRotatableBonds, NumAromaticRings, RingCount, FractionCSP3,
   HeavyAtomCount, NumAliphaticRings, NumSaturatedRings, NumHeteroatoms,
   LabuteASA.
3. **ChemBERTa-CLS:** `DeepChem/ChemBERTa-77M-MTR`, frozen, first-token
   representation, 384 dimensions, **max_length = 128**. The previous value
   of 512 in this sheet belonged to an incorrect provenance attribution;
   do not restore it. The generation notebook does not pin a ChemBERTa revision.
4. **MolFormer-Mean:** `ibm/MoLFormer-XL-both-10pct`, revision
   `7b12d946c181a37f6012b9dc3b002275de070314`, Transformers **4.44.2**,
   `trust_remote_code=True`, `deterministic_eval=True`, attention-mask-aware
   mean pooling, 768 dimensions, **max_length = 202**, frozen parameters.

**Generation notebook:** `scripts/GroupKFold_CV.ipynb`.
**Full comparison notebook:** `scripts/grid_benchmark.ipynb`.
`scripts/Grid_Benchmark.py` is a historical notebook-to-text export, not the
canonical standalone executable. A CLI conversion is not required merely
because the actual executable is a Colab notebook.

`embeddings_all.npz` stores the `rdkit`, `cb_cls`, `mf_mean`, `y` and `scaffold`
arrays. Its historical checksum is in `grid_final_metadata_20260801-152155.json`.
The file is not present in the checked repository; an external Google Drive
copy was mentioned by the author but has not been verified in this update.
It is useful for an exact historical-grid rerun; do not regenerate it by default.
The model names, MolFormer revision and generation code are already located.

## Regressors

**Random Forest:** `RandomForestRegressor`, `n_estimators=500`,
`random_state=42`, `n_jobs=-1`; remaining parameters recorded in run metadata.

**XGBoost:** `XGBRegressor`, `n_estimators=500`, `max_depth=6`,
`learning_rate=0.05`, `objective="reg:squarederror"`, `eval_metric="rmse"`,
`tree_method="hist"`, `random_state=42`, `n_jobs=-1`.

Each representation was evaluated with both regressors: all eight
combinations in this historical block were executed. This is not a claim
that the broader neural-versus-classical comparison is complete. D-MPNN and fine-tuning experiments are exploratory
and are not the basis for a general claim that neural models lost.

## Scaffold-aware validation and metrics

Five folds grouped by Bemis–Murcko scaffold; three partitions. Partition 0
uses deterministic GroupKFold; partitions 1 and 2 reassign scaffold groups
using seeds 1001 and 1002 with fold-size balancing. All eight combinations
share assignments. Each compound has one out-of-fold prediction per partition.

RMSE, R² and Spearman correlation are calculated on pooled predictions within
each partition and then summarized by mean and SD across partitions.
MAE is additionally reported in the random-versus-scaffold comparison.
A training-fold mean prediction is the baseline.

## Additional analyses (H1–H5 are internal file labels only)

### H1 — prediction on unseen scaffolds

The completed analysis compared frozen RF/XGBoost + Morgan scaffold
predictions with the mean baseline. Source: `FINAL_H1_scaffold_performance.csv`.

### H2 — published-value disagreement and prediction error

First collapse records by median within `compound × assay`; calculate dispersion
only for compounds in at least two assays (n=393). Separately collapse within
`compound × document` and use compounds in at least two documents (n=52).
Dispersion measures: range, sample SD, median absolute deviation, median
pairwise absolute difference. These groupings do not by themselves prove
independent technical repeats or independent laboratories.

Relate dispersion to absolute out-of-fold error averaged across scaffold
partitions using Spearman correlation. Use 1500 scaffold-cluster bootstrap
iterations for 95% CIs; partial rank correlation controls median pIC50 and
assay/document count. Holm correction applies to the four dispersion measures
within each model and aggregation level. Repeatedly measured subsets are
small/non-random; correlations are not a causal variance decomposition.

### H3 — auxiliary ligand-versus-random-decoy classification

Random unmatched ChEMBL controls; Morgan ECFP4; stratified 80/20 split,
seed 34; ten-fold stratified ROC AUC evaluation; RF, RBF-SVM and a neural
model labeled FFNN in the frozen summary. This is not an exact replication
of Mihai et al. or prospective validation of quantitative potency prediction.

`analyze_h1_h5.py` reads the frozen H3 summary, not a new classifier run.
The historical FFNN backend remains unverified because the original code
allowed an MLP fallback. Correcting the script does not retroactively establish
the backend that produced saved results. See root `reproducibility_notes.md`.

### H4 — chemical similarity and prediction error

For each test compound and fold, obtain maximum Morgan/Tanimoto similarity
to training compounds; average similarities and absolute errors across the
three partitions. Spearman correlation, 1500 scaffold-cluster bootstrap
iterations, partial rank correlation controlling median pIC50.
A weak association is not a calibrated individual uncertainty estimate.

### H5 — random versus scaffold validation

Same molecule-level data, Morgan features and model settings; shuffled
five-fold KFold with seeds 1000, 1001, 1002 compared with frozen scaffold
partitions. Metrics: RMSE, MAE, R² and Spearman.

The saved comparison includes 5000 scaffold-cluster bootstrap iterations and
10000 permutations; source `scripts/test_h5_validation_difference.py` and
`results/tables/FINAL_H5_significance.csv`. Manuscript v5 emphasizes differences
and CIs; do not automatically restore the removed H5 `p_Holm=0.00030` wording.
Historical statistics files remain unchanged. A separate environment manifest
for the significance script was not established; do not invent one.

## Software versions recorded for completed runs

| Package | Main grid | Additional analysis metadata |
|---|---|---|
| Python | 3.12.13 | 3.13.13 |
| NumPy | 2.0.2 | 2.4.6 |
| pandas | 2.2.2 | 3.0.3 |
| scikit-learn | 1.6.1 | 1.8.0 |
| SciPy | 1.16.3 | 1.17.1 |
| XGBoost | 3.3.0 | 3.2.0 |
| RDKit | 2026.03.4 | 2026.03.2 |

Sources: `grid_final_metadata_20260801-152155.json` and
`FINAL_H1_H5_METADATA.json`. These are historical environments, not a request
to update packages or assert identical environments for every separate run.

## Canonical result files

All paths below are relative to `results/tables/`:

```text
grid_final_results_20260801-152155.csv
grid_final_oof_20260801-152155.csv
grid_final_fold_assignments_20260801-152155.csv
grid_final_metadata_20260801-152155.json
FINAL_H1_scaffold_performance.csv
FINAL_H2_variability_vs_error.csv
FINAL_H3_mihai_replication_summary.csv
FINAL_H4_similarity_vs_error.csv
FINAL_H5_random_vs_scaffold.csv
FINAL_H5_random_oof_morgan.csv
FINAL_H5_significance.csv
FINAL_H1_H5_METADATA.json
```

## Completed corrections versus remaining work

`scripts/analyze_h1_h5.py` and `FINAL_H5_random_oof_morgan.csv` already exist.
The former request to rename/re-run the analysis to create them is obsolete.
The ChemBERTa token limit and exact MolFormer provenance above are resolved.

Offline dataset reconstruction and the two proposed robustness analyses remain
open as recorded in `STATUS.md`. The 27 records from 12 later-excluded assays
were retained in the main model input; v5 explicitly states this limitation.
Do not claim a robustness result that has not been generated.

Do not overwrite frozen files or edit historical paths/checksums merely to
make them portable. New outputs must be separate, with explicit provenance.

## Current neural-run audit: consult STATUS.md first

The preceding sections are descriptions of historical runs, not a to-do list.
Previous inspection identified an incomplete best-state restoration in
`scripts/MolFormer.ipynb` (head-only state while fine-tuning the backbone)
and differences in validation/reporting in `scripts/D_MPNN_CV.ipynb`.
Their fixes and impact on performance have not been established by a new run.
Neither these concerns nor the frozen-descriptor grid establish that neural
networks win or lose. Current model selection and training protocol remain open.
