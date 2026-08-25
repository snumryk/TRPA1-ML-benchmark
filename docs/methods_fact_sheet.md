# Methods fact sheet — TRPA1 manuscript

**Purpose:** source of truth for drafting the “Methods” section.  
**Status:** numbers frozen in GitHub on 2026-08-25.

## Approved research question

> Наскільки точно можна за молекулярною структурою прогнозувати агреговані за літературними даними значення pIC50 інгібування TRPA1 людини для хімічних каркасів, не представлених під час навчання, і чи пов’язана похибка такого прогнозу з розбіжністю опублікованих значень для тієї самої сполуки?

## Data source and selection

- Database: **ChEMBL 37**.
- Release date recorded in repository metadata: **2026-05-01**.
- Target: **CHEMBL6007, human TRPA1**.
- Endpoint: **IC50**.
- Relation: **exact values only (`standard_relation = "="`)**.
- pChEMBL value required.
- All retained standard units in the frozen raw table: **nM**.
- Frozen record-level table:
  `data/raw/trpa1_current_api_raw.csv`.
- Main molecule-level table:
  `data/processed/trpa1_primary_dataset.csv`.

## Dataset size

- Activity records: **2196**.
- Standardized compounds: **1645**.
- Assays: **97**.
- ChEMBL documents: **55**.
- Document years: **2010–2025**.
- Unique Bemis–Murcko scaffolds: **544**.
- Compounds with one measurement: **1196**.
- Compounds with at least two measurements: **449**.
- Compounds represented in at least two assays: **393**.
- Compounds represented in at least two documents: **52**.

## Structure standardization and aggregation

The frozen molecule-level table contains standardized SMILES, InChIKey,
Bemis–Murcko scaffold, median/minimum/maximum/SD of pChEMBL, and the number
of measurements.

Historical project documentation describes the following pipeline:

1. RDKit parsing of molecular structures.
2. Removal of salts / selection of the main molecular fragment.
3. Canonical isomeric SMILES generation.
4. Standardized InChIKey generation.
5. Deduplication by standardized structure.
6. Aggregation of reported pChEMBL values by the **median**.

**Important reproducibility gap:** the repository currently does not contain
one canonical script that rebuilds `trpa1_primary_dataset.csv` from the raw
records. The frozen input file and its benchmark SHA-256 are available, but
before formal submission the build script should be located or recreated.
This does not prevent drafting the manuscript, but it must not be hidden.

## Prediction target

- Target variable: `pchembl_median`.
- Meaning: median of retained reported pChEMBL values for the same
  standardized compound.
- It is an aggregated literature-derived value, not potency measured in one
  standardized protocol.

## Molecular representations

1. **Morgan ECFP4**
   - RDKit Morgan generator.
   - Radius: 2.
   - Vector size: 2048 bits.
   - Chirality: not included.

2. **RDKit-15**
   - MolWt, MolLogP, MolMR, TPSA, NumHAcceptors, NumHDonors,
     NumRotatableBonds, NumAromaticRings, RingCount, FractionCSP3,
     HeavyAtomCount, NumAliphaticRings, NumSaturatedRings,
     NumHeteroatoms, LabuteASA.

3. **ChemBERTa-CLS**
   - Model recorded in the repository code:
     `DeepChem/ChemBERTa-77M-MTR`.
   - Frozen pretrained model; no task-specific fine-tuning in the main grid.
   - CLS-token embedding, 384 dimensions.
   - SMILES tokenization with truncation, maximum length 512.

4. **MolFormer-Mean**
   - Frozen pretrained MoLFormer embedding.
   - Mean pooling, 768 dimensions.
   - Exact model revision must be copied from the original embedding-
     generation environment before submission; the benchmark metadata
     preserves the SHA-256 of `embeddings_all.npz`.

## Regressors

### Random Forest

- `RandomForestRegressor`
- `n_estimators = 500`
- `random_state = 42`
- `n_jobs = -1`
- Other parameters: scikit-learn defaults recorded in run metadata.

### XGBoost

- `XGBRegressor`
- `n_estimators = 500`
- `max_depth = 6`
- `learning_rate = 0.05`
- `objective = "reg:squarederror"`
- `eval_metric = "rmse"`
- `tree_method = "hist"`
- `random_state = 42`
- `n_jobs = -1`

## Scaffold-aware validation

- Five folds.
- Group variable: Bemis–Murcko scaffold.
- Three scaffold partitions.
- Partition 0: deterministic `GroupKFold`.
- Partitions 1–2: independently reassigned scaffold groups using seeds
  **1001** and **1002**, while balancing fold sizes.
- The same fold assignments were used for every representation–regressor
  combination.
- No scaffold was allowed in both training and validation within a fold.
- One out-of-fold prediction was saved for every compound in each
  partition.

## Metrics

Primary and secondary metrics:

- RMSE.
- R².
- Spearman rank correlation.
- MAE was additionally used in H5.

Pooled metrics were calculated from all out-of-fold predictions in each
partition. Mean and standard deviation across the three partitions were
then reported.

## H1 — prediction on unseen scaffolds

H1 was evaluated from the frozen scaffold OOF predictions. The model was
compared with a mean-prediction baseline; both RF and XGBoost were clearly
better than the baseline.

## H2 — disagreement of reported potency values and model error

1. Raw records were first collapsed to the median within each
   `compound × assay`.
2. Dispersion across assays was calculated only for compounds represented
   in at least two assays (**n = 393**).
3. A stricter analysis first collapsed records within
   `compound × document`, then calculated dispersion across documents for
   compounds represented in at least two documents
   (**n = 52**).
4. Dispersion measures:
   - range;
   - sample SD;
   - median absolute deviation;
   - median pairwise absolute difference.
5. Model error: absolute OOF error averaged across the three scaffold
   partitions.
6. Association: Spearman ρ.
7. 95% confidence intervals: scaffold-cluster bootstrap,
   **1500 iterations**.
8. Partial Spearman analysis controlled for median pIC50 and the number of
   assays/documents.
9. The repeated-measurement subset is non-random and enriched for more
   extensively studied compounds; this must be stated as a limitation.

## H3 — ligand-versus-random-decoy classification

- Antagonists versus randomly sampled ChEMBL decoys.
- Morgan ECFP4, 2048 bits.
- Random 80/20 stratified split, `random_state = 34`.
- Ten-fold stratified cross-validation for ROC AUC.
- RF, RBF-SVM and FFNN.
- This is a secondary comparison of task difficulty, not the main potency-
  prediction model.

## H4 — applicability-domain analysis

- For every test compound in every scaffold fold, the maximum Morgan
  Tanimoto similarity to all training compounds was calculated.
- Similarity was averaged across the three scaffold partitions.
- Association between similarity and absolute OOF error:
  Spearman ρ.
- 95% confidence intervals: scaffold-cluster bootstrap,
  1500 iterations.
- Partial analysis controlled for median pIC50.

## H5 — random versus scaffold validation

- Same molecule-level dataset.
- Same Morgan representation and RF/XGBoost settings.
- Random validation: shuffled five-fold `KFold`.
- Partition seeds: 1000, 1001 and 1002.
- Compared with the three frozen scaffold partitions.
- Metrics: RMSE, MAE, R² and Spearman ρ.

## Software versions for the main 4 × 2 benchmark

- Python 3.12.13
- NumPy 2.0.2
- pandas 2.2.2
- scikit-learn 1.6.1
- SciPy 1.16.3
- XGBoost 3.3.0
- RDKit 2026.03.4

## Canonical result files

- `results/tables/grid_final_results_20260801-152155.csv`
- `results/tables/grid_final_oof_20260801-152155.csv`
- `results/tables/grid_final_fold_assignments_20260801-152155.csv`
- `results/tables/grid_final_metadata_20260801-152155.json`
- `results/tables/FINAL_H2_variability_vs_error.csv`
- `results/tables/FINAL_H3_mihai_replication_summary.csv`
- `results/tables/FINAL_H4_similarity_vs_error.csv`
- `results/tables/FINAL_H5_random_vs_scaffold.csv`

## Repository corrections to make

1. Rename `scripts/analyze_hypotesis.py` to
   `scripts/analyze_h1_h5.py`.
2. Replace internal references to the misspelled filename.
3. Run the corrected script once to generate the currently absent audit
   file `results/tables/FINAL_H5_random_oof_morgan.csv`.
4. Preserve the frozen benchmark OOF file; do not regenerate it with newer
   package versions unless explicitly performing a new sensitivity run.
