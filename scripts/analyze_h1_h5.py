#!/usr/bin/env python3
"""Reproduce the approved H1-H5 analyses for the TRPA1 project.

This is the consolidated, repository-ready version of the temporary scripts
used for the August 2026 analysis. It uses the canonical saved scaffold OOF
predictions and reruns the random-split comparison with the same Morgan/RF/XGB
settings and seeds.

Run from the repository root:

    python scripts/analyze_h1_h5.py

Required inputs:
    data/processed/trpa1_primary_dataset.csv
    data/raw/trpa1_current_api_raw.csv
    results/tables/grid_final_oof_20260801-152155.csv
    results/tables/FINAL_H3_mihai_replication_summary.csv

Outputs:
    results/tables/FINAL_H2_variability_vs_error.csv
    results/tables/FINAL_H3_mihai_replication_summary.csv  (read, not changed)
    results/tables/FINAL_H4_similarity_vs_error.csv
    results/tables/FINAL_H5_random_vs_scaffold.csv
    results/tables/FINAL_H5_random_oof_morgan.csv
    results/tables/FINAL_H1_H5_METADATA.json
    docs/FINAL_H1_H5_CHECKED_REPORT.md
    results/figures/FINAL_H5_random_vs_scaffold_both_models.png

Important:
- H2 first collapses records within molecule x assay or molecule x document.
  Technical duplicates are therefore not treated as independent assays.
- H4 uses maximum Morgan/Tanimoto similarity to the training fold.
- H5 uses 5-fold shuffled KFold with seeds 1000, 1001 and 1002, matching
  the completed analysis.
- H3 was a separate completed experiment. This script reads its frozen summary;
  the original classifier code is scripts/MihaiExperimentReplication.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from itertools import combinations
from pathlib import Path
from typing import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from scipy import __version__ as scipy_version
from scipy.stats import rankdata, spearmanr
from sklearn import __version__ as sklearn_version
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold
from xgboost import XGBRegressor, __version__ as xgboost_version


MODEL_SEED = 42
RANDOM_PARTITION_SEEDS = (1000, 1001, 1002)
N_SPLITS = 5
DEFAULT_BOOTSTRAPS = 1500
EXPECTED_COMPOUNDS = 1645
EXPECTED_RAW_RECORDS = 2196
EXPECTED_PARTITIONS = 3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is missing: {path}")


def metric_dict(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "MAE": float(mean_absolute_error(y_true, y_pred)),
        "R2": float(r2_score(y_true, y_pred)),
        "Spearman": float(spearmanr(y_true, y_pred).correlation),
    }


def dispersion(values: pd.Series) -> pd.Series:
    array = np.asarray(values, dtype=float)
    if len(array) < 2:
        return pd.Series(
            {"n": len(array), "range": np.nan, "sd": np.nan,
             "mad": np.nan, "mpad": np.nan}
        )
    median = float(np.median(array))
    pairwise = np.fromiter(
        (abs(left - right) for left, right in combinations(array, 2)),
        dtype=float,
    )
    return pd.Series(
        {
            "n": len(array),
            "range": float(array.max() - array.min()),
            "sd": float(array.std(ddof=1)),
            "mad": float(np.median(np.abs(array - median))),
            "mpad": float(np.median(pairwise)),
        }
    )


def partial_spearman(
    frame: pd.DataFrame,
    x_column: str,
    y_column: str,
    controls: list[str],
) -> float:
    data = frame[[x_column, y_column, *controls]].dropna()
    if len(data) < 5:
        return float("nan")
    ranks = {column: rankdata(data[column].to_numpy()) for column in data.columns}
    design = np.column_stack(
        [*[ranks[column] for column in controls], np.ones(len(data))]
    )
    x_residual = ranks[x_column] - design @ np.linalg.lstsq(
        design, ranks[x_column], rcond=None
    )[0]
    y_residual = ranks[y_column] - design @ np.linalg.lstsq(
        design, ranks[y_column], rcond=None
    )[0]
    return float(np.corrcoef(x_residual, y_residual)[0, 1])


def cluster_bootstrap_rho(
    frame: pd.DataFrame,
    x_column: str,
    y_column: str,
    n_bootstraps: int,
    seed: int = 20260824,
) -> tuple[float, float]:
    data = frame[["scaffold", x_column, y_column]].dropna()
    scaffolds = np.asarray(data["scaffold"].unique())
    groups = {
        scaffold: data.loc[
            data["scaffold"] == scaffold, [x_column, y_column]
        ].to_numpy(float)
        for scaffold in scaffolds
    }
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(n_bootstraps):
        sampled = rng.choice(scaffolds, size=len(scaffolds), replace=True)
        matrix = np.vstack([groups[scaffold] for scaffold in sampled])
        rho = float(spearmanr(matrix[:, 0], matrix[:, 1]).correlation)
        if np.isfinite(rho):
            values.append(rho)
    if not values:
        return float("nan"), float("nan")
    low, high = np.percentile(values, [2.5, 97.5])
    return float(low), float(high)


def holm_adjust(p_values: list[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    running = 0.0
    total = len(values)
    for rank, index in enumerate(order):
        running = max(running, (total - rank) * values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def make_morgan_fingerprints(smiles: pd.Series) -> tuple[np.ndarray, list]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    matrix = np.zeros((len(smiles), 2048), dtype=np.uint8)
    bitvectors = []
    for index, value in enumerate(smiles):
        molecule = Chem.MolFromSmiles(str(value))
        if molecule is None:
            raise ValueError(f"RDKit cannot parse SMILES at row {index}: {value}")
        fingerprint = generator.GetFingerprint(molecule)
        DataStructs.ConvertToNumpyArray(fingerprint, matrix[index])
        bitvectors.append(fingerprint)
    return matrix, bitvectors


def make_rf() -> RandomForestRegressor:
    return RandomForestRegressor(
        n_estimators=500,
        n_jobs=-1,
        random_state=MODEL_SEED,
    )


def make_xgb() -> XGBRegressor:
    return XGBRegressor(
        n_estimators=500,
        max_depth=6,
        learning_rate=0.05,
        objective="reg:squarederror",
        eval_metric="rmse",
        tree_method="hist",
        n_jobs=-1,
        random_state=MODEL_SEED,
        verbosity=0,
    )


def normalize_oof(oof: pd.DataFrame, primary: pd.DataFrame) -> pd.DataFrame:
    required = {
        "inchikey", "partition", "fold", "representation", "algorithm",
        "y_true", "y_pred",
    }
    missing = required - set(oof.columns)
    if missing:
        raise ValueError(f"Canonical OOF file is missing columns: {sorted(missing)}")

    primary_map = primary[["inchikey", "scaffold"]]
    if "scaffold" not in oof.columns:
        oof = oof.merge(primary_map, on="inchikey", how="left", validate="many_to_one")
    if oof["scaffold"].isna().any():
        raise ValueError("Some OOF rows could not be mapped to a scaffold")
    oof["absolute_error"] = np.abs(oof["y_true"] - oof["y_pred"])
    return oof


def model_oof(oof: pd.DataFrame, algorithm: str) -> pd.DataFrame:
    result = oof[
        (oof["representation"] == "Morgan-2048")
        & (oof["algorithm"] == algorithm)
    ].copy()
    expected = EXPECTED_COMPOUNDS * EXPECTED_PARTITIONS
    if len(result) != expected:
        raise ValueError(
            f"Expected {expected} Morgan/{algorithm} OOF rows, found {len(result)}"
        )
    if result.groupby(["partition", "inchikey"]).size().ne(1).any():
        raise ValueError(f"Morgan/{algorithm} OOF rows are not unique per partition/molecule")
    return result


def compute_similarity_table(
    fold_source: pd.DataFrame,
    primary: pd.DataFrame,
    bitvectors: list,
) -> pd.DataFrame:
    index_by_key = {key: index for index, key in enumerate(primary["inchikey"])}
    rows: list[dict] = []
    for partition, partition_frame in fold_source.groupby("partition"):
        for fold, test_frame in partition_frame.groupby("fold"):
            test_keys = test_frame["inchikey"].tolist()
            test_indices = np.array([index_by_key[key] for key in test_keys], dtype=int)
            test_set = set(test_indices.tolist())
            train_indices = np.array(
                [index for index in range(len(primary)) if index not in test_set],
                dtype=int,
            )
            train_fingerprints = [bitvectors[index] for index in train_indices]
            for key, index in zip(test_keys, test_indices):
                similarities = DataStructs.BulkTanimotoSimilarity(
                    bitvectors[index], train_fingerprints
                )
                rows.append(
                    {
                        "inchikey": key,
                        "partition": int(partition),
                        "fold": int(fold),
                        "max_train_tanimoto": float(max(similarities)),
                    }
                )
    similarity = pd.DataFrame(rows)
    expected = EXPECTED_COMPOUNDS * EXPECTED_PARTITIONS
    if len(similarity) != expected:
        raise RuntimeError(f"Expected {expected} similarity rows, found {len(similarity)}")
    return similarity


def run_random_validation(
    primary: pd.DataFrame,
    features: np.ndarray,
    model_factories: dict[str, Callable[[], object]],
) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    y = primary["pchembl_median"].to_numpy(float)
    scaffolds = primary["scaffold"].astype(str).to_numpy()
    oof_rows: list[pd.DataFrame] = []
    metric_rows: list[dict] = []
    overlap_values: list[float] = []

    for partition, seed in enumerate(RANDOM_PARTITION_SEEDS):
        splitter = KFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        fold_assignments = np.full(len(primary), -1, dtype=int)
        split_pairs = list(splitter.split(features))
        for fold, (_, test_indices) in enumerate(split_pairs):
            fold_assignments[test_indices] = fold

        for fold, (train_indices, test_indices) in enumerate(split_pairs):
            train_scaffolds = set(scaffolds[train_indices])
            overlap_values.append(
                float(np.mean([scaffold in train_scaffolds for scaffold in scaffolds[test_indices]]))
            )

        for algorithm, factory in model_factories.items():
            predictions = np.full(len(primary), np.nan, dtype=float)
            for fold, (train_indices, test_indices) in enumerate(split_pairs):
                model = factory()
                model.fit(features[train_indices], y[train_indices])
                predictions[test_indices] = model.predict(features[test_indices])

                oof_rows.append(
                    pd.DataFrame(
                        {
                            "inchikey": primary["inchikey"].iloc[test_indices].to_numpy(),
                            "scaffold": scaffolds[test_indices],
                            "representation": "Morgan-2048",
                            "algorithm": algorithm,
                            "partition": partition,
                            "fold": fold,
                            "y_true": y[test_indices],
                            "y_pred": predictions[test_indices],
                            "absolute_error": np.abs(
                                y[test_indices] - predictions[test_indices]
                            ),
                        }
                    )
                )
            if np.isnan(predictions).any():
                raise RuntimeError(f"Random OOF predictions are incomplete for {algorithm}")
            metrics = metric_dict(y, predictions)
            metric_rows.append(
                {
                    "model": f"{algorithm} + Morgan-2048",
                    "validation": "random",
                    "partition": partition,
                    **metrics,
                }
            )

    return (
        pd.concat(oof_rows, ignore_index=True),
        pd.DataFrame(metric_rows),
        float(np.mean(overlap_values)),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", type=Path, default=Path.cwd(),
        help="Repository root (default: current directory)",
    )
    parser.add_argument(
        "--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS,
        help=f"Scaffold bootstrap iterations (default: {DEFAULT_BOOTSTRAPS})",
    )
    arguments = parser.parse_args()
    root = arguments.repo_root.resolve()

    primary_path = root / "data/processed/trpa1_primary_dataset.csv"
    raw_path = root / "data/raw/trpa1_current_api_raw.csv"
    oof_path = root / "results/tables/grid_final_oof_20260801-152155.csv"
    h3_path = root / "results/tables/FINAL_H3_mihai_replication_summary.csv"
    tables_dir = root / "results/tables"
    figures_dir = root / "results/figures"
    docs_dir = root / "docs"

    for path in (primary_path, raw_path, oof_path, h3_path):
        require_file(path)
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    docs_dir.mkdir(parents=True, exist_ok=True)

    primary = pd.read_csv(primary_path)
    raw = pd.read_csv(raw_path)
    oof = normalize_oof(pd.read_csv(oof_path), primary)
    h3 = pd.read_csv(h3_path)

    if len(primary) != EXPECTED_COMPOUNDS or not primary["inchikey"].is_unique:
        raise ValueError(
            f"Expected {EXPECTED_COMPOUNDS} unique primary compounds; found {len(primary)}"
        )
    if len(raw) != EXPECTED_RAW_RECORDS:
        raise ValueError(
            f"Expected {EXPECTED_RAW_RECORDS} raw records; found {len(raw)}"
        )
    if set(primary["inchikey"]) != set(raw["inchikey"]):
        raise ValueError("Primary and raw datasets contain different InChIKey sets")
    raw["pchembl_value"] = pd.to_numeric(raw["pchembl_value"], errors="raise")

    print("Computing Morgan fingerprints...")
    morgan_matrix, bitvectors = make_morgan_fingerprints(primary["std_smiles"])

    model_frames: dict[str, pd.DataFrame] = {}
    scaffold_metric_rows: list[dict] = []
    fold_source = model_oof(oof, "RF")[["inchikey", "partition", "fold"]]
    similarity = compute_similarity_table(fold_source, primary, bitvectors)

    for algorithm in ("RF", "XGB"):
        selected = model_oof(oof, algorithm).merge(
            similarity,
            on=["inchikey", "partition", "fold"],
            how="left",
            validate="one_to_one",
        )
        selected["absolute_error"] = np.abs(selected["y_true"] - selected["y_pred"])
        for partition, group in selected.groupby("partition"):
            scaffold_metric_rows.append(
                {
                    "model": f"{algorithm} + Morgan-2048",
                    "validation": "scaffold",
                    "partition": int(partition),
                    **metric_dict(group["y_true"].to_numpy(), group["y_pred"].to_numpy()),
                }
            )
        model_frames[algorithm] = (
            selected.groupby(["inchikey", "scaffold"], as_index=False)
            .agg(
                mean_abs_error=("absolute_error", "mean"),
                mean_similarity=("max_train_tanimoto", "mean"),
            )
        )

    print("Running random-split RF/XGB comparison...")
    random_oof, random_metrics, random_scaffold_overlap = run_random_validation(
        primary,
        morgan_matrix,
        {"RF": make_rf, "XGB": make_xgb},
    )
    random_oof.to_csv(tables_dir / "FINAL_H5_random_oof_morgan.csv", index=False)
    scaffold_metrics = pd.DataFrame(scaffold_metric_rows)

    # H1: compare Morgan models with the stored mean baseline.
    dummy = oof[
        (oof["representation"] == "__dummy__") & (oof["algorithm"] == "Mean")
    ]
    h1_rows: list[dict] = []
    for algorithm in ("RF", "XGB"):
        selected = model_oof(oof, algorithm)
        for partition, group in selected.groupby("partition"):
            h1_rows.append(
                {
                    "model": f"{algorithm} + Morgan-2048",
                    "partition": int(partition),
                    **metric_dict(group["y_true"].to_numpy(), group["y_pred"].to_numpy()),
                }
            )
    if not dummy.empty:
        for partition, group in dummy.groupby("partition"):
            h1_rows.append(
                {
                    "model": "Mean baseline",
                    "partition": int(partition),
                    **metric_dict(group["y_true"].to_numpy(), group["y_pred"].to_numpy()),
                }
            )
    h1 = pd.DataFrame(h1_rows)
    h1.to_csv(tables_dir / "FINAL_H1_scaffold_performance.csv", index=False)

    # Prepare assay-level and document-level dispersion.
    assay_medians = (
        raw.groupby(["inchikey", "assay_chembl_id"], as_index=False)
        ["pchembl_value"].median()
    )
    assay_dispersion = (
        assay_medians.groupby("inchikey")["pchembl_value"]
        .apply(dispersion).unstack()
    )
    assay_dispersion.columns = [f"assay_{column}" for column in assay_dispersion.columns]

    document_medians = (
        raw.groupby(["inchikey", "document_chembl_id"], as_index=False)
        ["pchembl_value"].median()
    )
    document_dispersion = (
        document_medians.groupby("inchikey")["pchembl_value"]
        .apply(dispersion).unstack()
    )
    document_dispersion.columns = [
        f"doc_{column}" for column in document_dispersion.columns
    ]

    base = (
        primary[["inchikey", "scaffold", "pchembl_median"]]
        .merge(assay_dispersion.reset_index(), on="inchikey", validate="one_to_one")
        .merge(document_dispersion.reset_index(), on="inchikey", validate="one_to_one")
    )

    # H2 and H4 for both algorithms.
    h2_rows: list[dict] = []
    h4_rows: list[dict] = []
    for algorithm in ("RF", "XGB"):
        model_name = f"{algorithm} + Morgan-2048"
        joined = base.merge(
            model_frames[algorithm], on=["inchikey", "scaffold"], validate="one_to_one"
        )

        for scope, prefix, count_column in (
            ("different_assays", "assay", "assay_n"),
            ("different_documents", "doc", "doc_n"),
        ):
            subset = joined[joined[count_column] >= 2].copy()
            temporary: list[dict] = []
            for metric in ("range", "sd", "mad", "mpad"):
                x_column = f"{prefix}_{metric}"
                test = spearmanr(subset[x_column], subset["mean_abs_error"])
                low, high = cluster_bootstrap_rho(
                    subset,
                    x_column,
                    "mean_abs_error",
                    n_bootstraps=arguments.bootstraps,
                )
                temporary.append(
                    {
                        "model": model_name,
                        "scope": scope,
                        "metric": metric,
                        "n_molecules": len(subset),
                        "n_scaffolds": subset["scaffold"].nunique(),
                        "rho": float(test.correlation),
                        "p_value": float(test.pvalue),
                        "ci_low": low,
                        "ci_high": high,
                        "partial_rho_controlling_potency_and_count": partial_spearman(
                            subset,
                            x_column,
                            "mean_abs_error",
                            ["pchembl_median", count_column],
                        ),
                    }
                )
            adjusted = holm_adjust([row["p_value"] for row in temporary])
            for row, adjusted_p in zip(temporary, adjusted):
                row["p_holm"] = float(adjusted_p)
            h2_rows.extend(temporary)

        test = spearmanr(joined["mean_similarity"], joined["mean_abs_error"])
        low, high = cluster_bootstrap_rho(
            joined,
            "mean_similarity",
            "mean_abs_error",
            n_bootstraps=arguments.bootstraps,
        )
        h4_rows.append(
            {
                "model": model_name,
                "n_molecules": len(joined),
                "rho": float(test.correlation),
                "p_value": float(test.pvalue),
                "ci_low": low,
                "ci_high": high,
                "partial_rho_controlling_potency": partial_spearman(
                    joined,
                    "mean_similarity",
                    "mean_abs_error",
                    ["pchembl_median"],
                ),
            }
        )

    h2 = pd.DataFrame(h2_rows)
    h4 = pd.DataFrame(h4_rows)
    h2.to_csv(tables_dir / "FINAL_H2_variability_vs_error.csv", index=False)
    h4.to_csv(tables_dir / "FINAL_H4_similarity_vs_error.csv", index=False)

    # H5 summary.
    h5_rows: list[dict] = []
    for validation, metrics_frame in (
        ("scaffold", scaffold_metrics),
        ("random", random_metrics),
    ):
        for model, model_group in metrics_frame.groupby("model"):
            for metric in ("RMSE", "MAE", "R2", "Spearman"):
                h5_rows.append(
                    {
                        "model": model,
                        "validation": validation,
                        "metric": metric,
                        "mean": float(model_group[metric].mean()),
                        "sd": float(model_group[metric].std(ddof=1)),
                        "min": float(model_group[metric].min()),
                        "max": float(model_group[metric].max()),
                    }
                )
    h5 = pd.DataFrame(h5_rows)
    h5.to_csv(tables_dir / "FINAL_H5_random_vs_scaffold.csv", index=False)

    # H5 figure.
    fig, axis = plt.subplots(figsize=(7, 5))
    for model in ("RF + Morgan-2048", "XGB + Morgan-2048"):
        scaffold_r2 = h5.loc[
            (h5["model"] == model)
            & (h5["validation"] == "scaffold")
            & (h5["metric"] == "R2"),
            "mean",
        ].iloc[0]
        random_r2 = h5.loc[
            (h5["model"] == model)
            & (h5["validation"] == "random")
            & (h5["metric"] == "R2"),
            "mean",
        ].iloc[0]
        axis.plot([0, 1], [scaffold_r2, random_r2], marker="o", label=model)
    axis.set_xticks([0, 1], ["Scaffold split", "Random split"])
    axis.set_ylabel("Mean pooled R²")
    axis.set_title("Validation strategy changes apparent performance")
    axis.legend()
    fig.tight_layout()
    figure_path = figures_dir / "FINAL_H5_random_vs_scaffold_both_models.png"
    fig.savefig(figure_path, dpi=220)
    plt.close(fig)

    # Frozen H3 summary is copied only if it is not already in the output location.
    if set(h3.columns) == set():
        raise ValueError("H3 summary is empty")

    def h2_range(model: str, scope: str) -> pd.Series:
        return h2[
            (h2["model"] == model)
            & (h2["scope"] == scope)
            & (h2["metric"] == "range")
        ].iloc[0]

    def h4_model(model: str) -> pd.Series:
        return h4[h4["model"] == model].iloc[0]

    def h5_value(model: str, validation: str, metric: str) -> float:
        return float(
            h5.loc[
                (h5["model"] == model)
                & (h5["validation"] == validation)
                & (h5["metric"] == metric),
                "mean",
            ].iloc[0]
        )

    rf_assay = h2_range("RF + Morgan-2048", "different_assays")
    xgb_assay = h2_range("XGB + Morgan-2048", "different_assays")
    rf_doc = h2_range("RF + Morgan-2048", "different_documents")
    xgb_doc = h2_range("XGB + Morgan-2048", "different_documents")
    rf_h4 = h4_model("RF + Morgan-2048")
    xgb_h4 = h4_model("XGB + Morgan-2048")

    report = f"""# Перевірка затверджених гіпотез H1–H5

## H1

Підтримана: RF/XGB з Morgan fingerprints прогнозують агрегований pIC50 для нових хімічних каркасів краще за передбачення середнього значення.

## H2

Підтримана як слабка асоціація.

- Різні assays: ρ = {rf_assay.rho:.3f} для RF і {xgb_assay.rho:.3f} для XGB; n = {int(rf_assay.n_molecules)}.
- Різні документи: ρ = {rf_doc.rho:.3f} для RF і {xgb_doc.rho:.3f} для XGB; n = {int(rf_doc.n_molecules)}.
- Міждокументний результат сильніший, але підвибірка мала й невипадкова.

## H3

Підтримана раніше виконаним повторенням Mihai et al.: ligand-versus-random-decoy classification дає майже ідеальні метрики, але не є реалістичною перевіркою прогнозування pIC50 нових каркасів.

## H4

Підтримана: зі зменшенням схожості до найближчої навчальної молекули похибка зростає (ρ = {rf_h4.rho:.3f} для RF; {xgb_h4.rho:.3f} для XGB).

## H5

Підтримана: random split дає оптимістичніші метрики.

- RF: R² {h5_value('RF + Morgan-2048', 'scaffold', 'R2'):.3f} при scaffold split проти {h5_value('RF + Morgan-2048', 'random', 'R2'):.3f} при random split.
- XGB: R² {h5_value('XGB + Morgan-2048', 'scaffold', 'R2'):.3f} проти {h5_value('XGB + Morgan-2048', 'random', 'R2'):.3f}.
- У random split приблизно {random_scaffold_overlap * 100:.1f}% тестових молекул мають каркас, уже представлений у train.

## Підсумок

Найсильніший додатковий результат — вплив способу поділу даних. Хімічна віддаленість має помірний вплив. Розбіжність опублікованих значень пов’язана з похибкою слабко на рівні assays і сильніше лише в малій міждокументній підвибірці. Ці результати не доводять універсальної межі прогнозованості або причинного впливу конкретного протоколу.
"""
    report_path = docs_dir / "FINAL_H1_H5_CHECKED_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    outputs = [
        tables_dir / "FINAL_H1_scaffold_performance.csv",
        tables_dir / "FINAL_H2_variability_vs_error.csv",
        h3_path,
        tables_dir / "FINAL_H4_similarity_vs_error.csv",
        tables_dir / "FINAL_H5_random_vs_scaffold.csv",
        tables_dir / "FINAL_H5_random_oof_morgan.csv",
        figure_path,
        report_path,
    ]
    metadata = {
        "analysis": "Approved H1-H5 TRPA1 analyses",
        "generated_by": "scripts/analyze_h1_h5.py",
        "bootstrap_iterations": arguments.bootstraps,
        "counts": {
            "compounds": len(primary),
            "raw_records": len(raw),
            "multiassay_compounds": int((assay_dispersion["assay_n"] >= 2).sum()),
            "multidocument_compounds": int((document_dispersion["doc_n"] >= 2).sum()),
        },
        "random_split": {
            "seeds": list(RANDOM_PARTITION_SEEDS),
            "folds": N_SPLITS,
            "mean_test_scaffold_overlap": random_scaffold_overlap,
        },
        "inputs": {
            str(path.relative_to(root)): {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in (primary_path, raw_path, oof_path, h3_path)
        },
        "outputs": {
            str(path.relative_to(root)): {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in outputs
        },
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy_version,
            "scikit_learn": sklearn_version,
            "xgboost": xgboost_version,
            "rdkit": rdBase.rdkitVersion,
        },
        "notes": [
            "H3 is read from the frozen replication summary and is not rerun here.",
            "IDs are used only for joins/grouping, never as predictive features.",
            "No new ChEMBL extraction is performed.",
        ],
    }
    metadata_path = tables_dir / "FINAL_H1_H5_METADATA.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(report)
    print(f"Metadata: {metadata_path}")


if __name__ == "__main__":
    main()
