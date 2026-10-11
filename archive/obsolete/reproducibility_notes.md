# Reproducibility notes — TRPA1 benchmark

**Updated:** 27 August 2026

This file records historical names, exact model provenance, and distinctions
between scientific blockers and optional engineering cleanup.

## 1. Dataset chronology

The main molecule-level benchmark input was originally created from ChEMBL 36
and stored as:

```text
trpa1_antagonists.csv
```

The current repository name is:

```text
data/processed/trpa1_primary_dataset.csv
```

A ChEMBL 37 re-extraction with the same strict filter and simplified RDKit
standardization produced:

```text
1645 compounds in ChEMBL 36
1645 compounds in ChEMBL 37
1645 overlapping compounds
0 new compounds
0 removed compounds
0 changed aggregate values
```

Therefore the frozen molecule-level benchmark dataset was unchanged by the
release transition.

Relevant files:

```text
scripts/chembl37_delta.py
scripts/resolve_delta.py
data/raw/trpa1_current_api_metadata.json
data/raw/trpa1_current_api_raw.csv
data/processed/trpa1_current_api_aggregated.csv
data/processed/trpa1_primary_dataset.csv
```

## 2. Historical metadata must remain historical

`results/tables/grid_final_metadata_20260801-152155.json` records the actual
Google Drive filename used during the benchmark:

```text
trpa1_antagonists.csv
```

Do not rewrite frozen metadata merely because the repository file was later
renamed. Document the mapping instead.

Windows-style paths in frozen metadata are also historical provenance. New
scripts should write POSIX-style repository-relative paths, but old metadata
should not be cosmetically edited.

## 3. Line endings and SHA-256

The repository uses:

```gitattributes
* text=auto eol=lf
```

This is appropriate and should remain unchanged.

A CSV written with CRLF on Google Drive/Windows and the same CSV stored with
LF in Git may have identical parsed content but different byte-level SHA-256.

Rules for future generated CSV files:

```python
frame.to_csv(path, index=False, lineterminator="\n")
```

Future metadata should record:

- byte-level SHA-256 of the exact file used;
- row count;
- column names and order;
- optionally a normalized-content hash based on LF line endings.

Do not modify the historical benchmark hashes.

## 4. Embeddings provenance

Generation notebook:

```text
scripts/GroupKFold_CV.ipynb
```

ChemBERTa:

```text
model = DeepChem/ChemBERTa-77M-MTR
pooling = CLS / first token
max_length = 128
dimensions = 384
```

MolFormer:

```text
model = ibm/MoLFormer-XL-both-10pct
revision = 7b12d946c181a37f6012b9dc3b002275de070314
transformers = 4.44.2
trust_remote_code = true
deterministic_eval = true
pooling = attention-mask-aware mean
max_length = 202
dimensions = 768
```

`embeddings_all.npz` contains:

```text
rdkit
cb_cls
mf_mean
y
scaffold
```

Historical SHA-256:

```text
682d13c98d0565b8e9c79342cb985f27476c43e0846bf2f0af911400685c726e
```

The NPZ file is useful for an exact rerun of the full representation grid.
It is not required to draft the manuscript or rerun Morgan-only H1–H5.

## 5. Canonical benchmark executable

Canonical executable notebook:

```text
scripts/grid_benchmark.ipynb
```

Historical text export:

```text
scripts/Grid_Benchmark.py
```

The latter should be labeled legacy/not executable unless replaced with a
clean standalone script. This is engineering cleanup, not evidence that the
saved OOF predictions or metrics are invalid.

## 6. H3 terminology and backend

The H3 analysis is not an exact replication of Mihai et al. because it uses:

- Morgan fingerprints instead of MNA descriptors;
- a different dataset size;
- randomly sampled, unmatched ChEMBL decoys.

Correct terminology:

```text
auxiliary reproduction of the ligand-versus-random-decoy setup
```

The historical script silently fell back from TensorFlow/Keras FFNN to
`sklearn.neural_network.MLPClassifier`. This is unacceptable because the
reported label `FFNN` then does not uniquely identify the trained model.

The corrected script must fail explicitly when TensorFlow is unavailable.
Historical frozen results should not be relabeled as a confirmed Keras FFNN
unless the original execution environment or stdout proves the backend.

RF and SVM results alone are already sufficient to demonstrate the
task-framing contrast.

## 7. Target aggregation robustness check

Frozen benchmark target:

```text
pchembl_median = median of all retained activity records for a compound
```

Potential alternative:

```text
1. median within compound × assay
2. median across assay-level medians
```

The alternative gives equal weight to assays instead of equal weight to
individual records. It is a sensitivity analysis, not a silent replacement of
the frozen benchmark target.

## 8. Post-audit sensitivity check

The original benchmark used all 2196 retained records.

Later assay audit classified:

```text
primary:     1953 records
sensitivity: 216 records
excluded:      27 records
```

The 27 excluded records belong to 12 assays and involve 22 compounds.

A minimal robustness run should:

```text
remove the 27 records
rebuild molecule-level targets
rerun only Morgan + RF and Morgan + XGBoost
compare scaffold-CV metrics
```

No full 4 × 2 rerun is required.
