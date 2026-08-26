#!/usr/bin/env python3
"""Cluster-aware significance test for H5: random split vs scaffold split.

Run after:
    python scripts/analyze_h1_h5.py --repo-root . --bootstraps 1500

Required inputs:
    data/processed/trpa1_primary_dataset.csv
    results/tables/grid_final_oof_20260801-152155.csv
    results/tables/FINAL_H5_random_oof_morgan.csv

Output:
    results/tables/FINAL_H5_significance.csv

The unit resampled/permuted is the Bemis–Murcko scaffold, not an individual
molecule. This respects the dependence among compounds sharing a scaffold.

For each validation strategy, metrics are calculated separately for each of
the three saved partitions and then averaged, matching the existing H5 table.

Positive ``random_improvement`` means that random splitting produced a more
favourable result:
- RMSE and MAE: scaffold value minus random value;
- R2: random value minus scaffold value.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ALGORITHMS = ("RF", "XGB")
REPRESENTATION = "Morgan-2048"
EXPECTED_PARTITIONS = (0, 1, 2)


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is missing: {path}")


def holm_adjust(p_values: list[float]) -> np.ndarray:
    """Holm step-down adjustment, preserving original order."""
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted_sorted = np.maximum.accumulate(
        (len(values) - np.arange(len(values))) * values[order]
    )
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return adjusted


def prepare_oof(
    frame: pd.DataFrame,
    primary: pd.DataFrame,
    validation: str,
    algorithm: str,
) -> pd.DataFrame:
    required = {
        "inchikey", "partition", "representation", "algorithm",
        "y_true", "y_pred",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"{validation} OOF is missing columns: {sorted(missing)}"
        )

    selected = frame[
        (frame["representation"] == REPRESENTATION)
        & (frame["algorithm"] == algorithm)
    ].copy()

    if selected.empty:
        raise ValueError(
            f"No {validation} rows found for {algorithm} + {REPRESENTATION}"
        )

    if "scaffold" not in selected.columns:
        selected = selected.merge(
            primary[["inchikey", "scaffold"]],
            on="inchikey",
            how="left",
            validate="many_to_one",
        )

    if selected["scaffold"].isna().any():
        raise ValueError(f"{validation}: some rows have no scaffold")

    selected["partition"] = selected["partition"].astype(int)
    found_partitions = tuple(sorted(selected["partition"].unique()))
    if found_partitions != EXPECTED_PARTITIONS:
        raise ValueError(
            f"{validation}: expected partitions {EXPECTED_PARTITIONS}, "
            f"found {found_partitions}"
        )

    duplicate = selected.duplicated(["inchikey", "partition"])
    if duplicate.any():
        raise ValueError(
            f"{validation}: duplicate molecule × partition rows found"
        )

    expected_rows = len(primary) * len(EXPECTED_PARTITIONS)
    if len(selected) != expected_rows:
        raise ValueError(
            f"{validation}: expected {expected_rows} rows, found {len(selected)}"
        )

    if set(selected["inchikey"]) != set(primary["inchikey"]):
        raise ValueError(
            f"{validation}: molecule set differs from primary dataset"
        )

    selected["squared_error"] = (
        selected["y_true"].astype(float) - selected["y_pred"].astype(float)
    ) ** 2
    selected["absolute_error"] = np.abs(
        selected["y_true"].astype(float) - selected["y_pred"].astype(float)
    )
    selected["y_sum"] = selected["y_true"].astype(float)
    selected["y_sq_sum"] = selected["y_true"].astype(float) ** 2
    selected["validation"] = validation
    return selected


def aggregate_by_scaffold_partition(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.groupby(["scaffold", "partition"], as_index=False)
        .agg(
            sse=("squared_error", "sum"),
            sae=("absolute_error", "sum"),
            n=("inchikey", "size"),
            sum_y=("y_sum", "sum"),
            sum_y2=("y_sq_sum", "sum"),
        )
    )


def align_arrays(
    scaffold_frame: pd.DataFrame,
    random_frame: pd.DataFrame,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray]]:
    keys = ["scaffold", "partition"]
    merged = scaffold_frame.merge(
        random_frame,
        on=keys,
        how="inner",
        suffixes=("_scaffold", "_random"),
        validate="one_to_one",
    )

    expected = scaffold_frame.shape[0]
    if len(merged) != expected or len(random_frame) != expected:
        raise ValueError("Scaffold and random OOF summaries do not align")

    for field in ("n", "sum_y", "sum_y2"):
        left = merged[f"{field}_scaffold"].to_numpy(float)
        right = merged[f"{field}_random"].to_numpy(float)
        if not np.allclose(left, right, rtol=0, atol=1e-12):
            raise ValueError(
                f"Outcome/sample summaries differ between validations: {field}"
            )

    scaffolds = np.asarray(sorted(merged["scaffold"].unique()), dtype=object)
    scaffold_index = {value: index for index, value in enumerate(scaffolds)}
    partition_index = {
        value: index for index, value in enumerate(EXPECTED_PARTITIONS)
    }
    shape = (len(scaffolds), len(EXPECTED_PARTITIONS))

    common = {
        "n": np.zeros(shape, dtype=float),
        "sum_y": np.zeros(shape, dtype=float),
        "sum_y2": np.zeros(shape, dtype=float),
    }
    scaffold_values = {
        "sse": np.zeros(shape, dtype=float),
        "sae": np.zeros(shape, dtype=float),
    }
    random_values = {
        "sse": np.zeros(shape, dtype=float),
        "sae": np.zeros(shape, dtype=float),
    }

    for row in merged.itertuples(index=False):
        i = scaffold_index[row.scaffold]
        j = partition_index[int(row.partition)]
        common["n"][i, j] = float(row.n_scaffold)
        common["sum_y"][i, j] = float(row.sum_y_scaffold)
        common["sum_y2"][i, j] = float(row.sum_y2_scaffold)
        scaffold_values["sse"][i, j] = float(row.sse_scaffold)
        scaffold_values["sae"][i, j] = float(row.sae_scaffold)
        random_values["sse"][i, j] = float(row.sse_random)
        random_values["sae"][i, j] = float(row.sae_random)

    return scaffolds, common, {
        "scaffold_sse": scaffold_values["sse"],
        "scaffold_sae": scaffold_values["sae"],
        "random_sse": random_values["sse"],
        "random_sae": random_values["sae"],
    }


def mean_partition_metrics(
    indices: np.ndarray,
    common: dict[str, np.ndarray],
    values: dict[str, np.ndarray],
) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}

    for validation in ("scaffold", "random"):
        rmse_values: list[float] = []
        mae_values: list[float] = []
        r2_values: list[float] = []

        for partition_position in range(len(EXPECTED_PARTITIONS)):
            n = common["n"][indices, partition_position].sum()
            sum_y = common["sum_y"][indices, partition_position].sum()
            sum_y2 = common["sum_y2"][indices, partition_position].sum()
            sse = values[f"{validation}_sse"][
                indices, partition_position
            ].sum()
            sae = values[f"{validation}_sae"][
                indices, partition_position
            ].sum()

            if n <= 1:
                raise ValueError("Bootstrap/permutation sample is too small")

            sst = sum_y2 - (sum_y ** 2) / n
            if sst <= 0:
                raise ValueError("Outcome variance is zero")

            rmse_values.append(float(np.sqrt(sse / n)))
            mae_values.append(float(sae / n))
            r2_values.append(float(1.0 - sse / sst))

        result[validation] = {
            "RMSE": float(np.mean(rmse_values)),
            "MAE": float(np.mean(mae_values)),
            "R2": float(np.mean(r2_values)),
        }

    return result


def improvement(metrics: dict[str, dict[str, float]]) -> dict[str, float]:
    return {
        "RMSE": metrics["scaffold"]["RMSE"] - metrics["random"]["RMSE"],
        "MAE": metrics["scaffold"]["MAE"] - metrics["random"]["MAE"],
        "R2": metrics["random"]["R2"] - metrics["scaffold"]["R2"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="Repository root (default: current directory)",
    )
    parser.add_argument(
        "--bootstraps",
        type=int,
        default=5000,
        help="Scaffold bootstrap iterations (default: 5000)",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=10000,
        help="Scaffold label-swap permutations (default: 10000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260826,
        help="Random seed",
    )
    args = parser.parse_args()

    if args.bootstraps < 100:
        raise ValueError("--bootstraps must be at least 100")
    if args.permutations < 100:
        raise ValueError("--permutations must be at least 100")

    root = args.repo_root.resolve()
    primary_path = root / "data/processed/trpa1_primary_dataset.csv"
    scaffold_oof_path = (
        root / "results/tables/grid_final_oof_20260801-152155.csv"
    )
    random_oof_path = (
        root / "results/tables/FINAL_H5_random_oof_morgan.csv"
    )
    output_path = root / "results/tables/FINAL_H5_significance.csv"

    for path in (primary_path, scaffold_oof_path, random_oof_path):
        require_file(path)

    primary = pd.read_csv(primary_path)
    required_primary = {"inchikey", "scaffold"}
    missing_primary = required_primary - set(primary.columns)
    if missing_primary:
        raise ValueError(
            f"Primary dataset is missing: {sorted(missing_primary)}"
        )
    if primary["inchikey"].duplicated().any():
        raise ValueError("Primary dataset has duplicate InChIKeys")

    scaffold_oof = pd.read_csv(scaffold_oof_path)
    random_oof = pd.read_csv(random_oof_path)
    rng = np.random.default_rng(args.seed)

    rows: list[dict] = []

    for algorithm in ALGORITHMS:
        scaffold_selected = prepare_oof(
            scaffold_oof, primary, "scaffold", algorithm
        )
        random_selected = prepare_oof(
            random_oof, primary, "random", algorithm
        )

        scaffold_summary = aggregate_by_scaffold_partition(scaffold_selected)
        random_summary = aggregate_by_scaffold_partition(random_selected)
        scaffolds, common, values = align_arrays(
            scaffold_summary, random_summary
        )

        all_indices = np.arange(len(scaffolds), dtype=int)
        observed_metrics = mean_partition_metrics(
            all_indices, common, values
        )
        observed_improvement = improvement(observed_metrics)

        bootstrap_values = {
            metric: np.empty(args.bootstraps, dtype=float)
            for metric in ("RMSE", "MAE", "R2")
        }
        for iteration in range(args.bootstraps):
            sampled = rng.integers(
                0, len(scaffolds), size=len(scaffolds), endpoint=False
            )
            current = improvement(
                mean_partition_metrics(sampled, common, values)
            )
            for metric in bootstrap_values:
                bootstrap_values[metric][iteration] = current[metric]

        permutation_values = {
            metric: np.empty(args.permutations, dtype=float)
            for metric in ("RMSE", "MAE", "R2")
        }
        for iteration in range(args.permutations):
            swap = rng.random(len(scaffolds)) < 0.5
            permuted = {
                key: array.copy() for key, array in values.items()
            }
            for error_type in ("sse", "sae"):
                left_key = f"scaffold_{error_type}"
                right_key = f"random_{error_type}"
                left = values[left_key]
                right = values[right_key]
                permuted[left_key][swap, :] = right[swap, :]
                permuted[right_key][swap, :] = left[swap, :]

            current = improvement(
                mean_partition_metrics(all_indices, common, permuted)
            )
            for metric in permutation_values:
                permutation_values[metric][iteration] = current[metric]

        model_rows: list[dict] = []
        for metric in ("RMSE", "MAE", "R2"):
            observed = observed_improvement[metric]
            low, high = np.percentile(
                bootstrap_values[metric], [2.5, 97.5]
            )
            null = permutation_values[metric]
            p_value = (
                1.0 + np.sum(np.abs(null) >= abs(observed))
            ) / (args.permutations + 1.0)

            model_rows.append(
                {
                    "model": f"{algorithm} + {REPRESENTATION}",
                    "metric": metric,
                    "scaffold_mean": observed_metrics["scaffold"][metric],
                    "random_mean": observed_metrics["random"][metric],
                    "random_improvement": observed,
                    "ci_low": float(low),
                    "ci_high": float(high),
                    "p_scaffold_permutation": float(p_value),
                    "n_compounds": len(primary),
                    "n_scaffolds": len(scaffolds),
                    "n_partitions": len(EXPECTED_PARTITIONS),
                    "bootstrap_iterations": args.bootstraps,
                    "permutation_iterations": args.permutations,
                    "seed": args.seed,
                }
            )

        adjusted = holm_adjust(
            [row["p_scaffold_permutation"] for row in model_rows]
        )
        for row, adjusted_p in zip(model_rows, adjusted):
            row["p_holm_within_model"] = float(adjusted_p)

        rows.extend(model_rows)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)

    print(result.to_string(index=False))
    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    main()
