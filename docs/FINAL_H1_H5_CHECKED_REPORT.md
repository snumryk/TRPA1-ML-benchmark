# FINAL CHECK — approved H1–H5

- **H1 supported:** scaffold-CV R² is about 0.57, clearly above the mean baseline near zero.
- **H2 supported, but weak across assays:** rho 0.127 (RF) and 0.164 (XGB), n=393. Across different documents the association is stronger — rho 0.347 and 0.449 — but n=52 and the subset is selected.
- **H3 supported by the existing replication:** accuracy 98.64–99.09%, CV AUC 0.9968–0.9983 for ligand versus random decoy classification.
- **H4 supported:** nearest-training-compound similarity correlates negatively with error, rho -0.177 (RF) and -0.177 (XGB).
- **H5 supported:** RF R² 0.572 scaffold vs 0.640 random; XGB 0.567 vs 0.638.

**Interpretation:** random splitting is the strongest source of optimism; chemical distance has a smaller but reproducible effect; disagreement among published potency values explains only a modest part of residual error. Do not claim a universal prediction ceiling or causal assay effect.
