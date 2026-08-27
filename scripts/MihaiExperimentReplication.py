#!/usr/bin/env python3
"""Auxiliary reproduction of a ligand-versus-random-decoy TRPA1 setup.

This is NOT an exact replication of Mihai et al.:
- Morgan ECFP4 fingerprints replace MNA descriptors;
- the active and decoy sets differ in size;
- random ChEMBL decoys are not property-matched;
- validation is random/stratified rather than scaffold-held-out.

Scientific purpose:
Demonstrate that ligand-versus-random-unmatched-decoy classification can be
nearly trivial and therefore must not be confused with scaffold-aware pIC50
regression or prospective virtual-screening validation.

The script deliberately has NO silent neural-network fallback. TensorFlow is
required for the Keras FFNN. If it is unavailable, the run stops explicitly.
"""

from __future__ import annotations

import argparse
import json
import platform
import random
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from sklearn import __version__ as sklearn_version
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.svm import SVC

SEED = 34
FP_SIZE = 2048
N_SPLITS = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="Repository root; defaults to the current directory.",
    )
    parser.add_argument(
        "--actives",
        type=Path,
        default=Path("data/processed/trpa1_primary_dataset.csv"),
    )
    parser.add_argument(
        "--decoys",
        type=Path,
        default=Path("data/processed/decoys_clean.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/H3_auxiliary_random_decoy_summary.csv"),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("results/tables/H3_auxiliary_random_decoy_metadata.json"),
    )
    return parser.parse_args()


def resolve(root: Path, value: Path) -> Path:
    return value if value.is_absolute() else root / value


def smiles_to_fp(smiles: str, generator) -> np.ndarray:
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    fingerprint = generator.GetFingerprint(molecule)
    array = np.zeros((FP_SIZE,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fingerprint, array)
    return array


def metrics(y_true: np.ndarray, y_pred: np.ndarray, y_prob: np.ndarray) -> dict:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    tpr = tp / (tp + fn)
    tnr = tn / (tn + fp)
    return {
        "TPR": float(tpr),
        "TNR": float(tnr),
        "Accuracy": float((tp + tn) / (tp + tn + fp + fn)),
        "balanced_accuracy": float((tpr + tnr) / 2),
        "FPR": float(fp / (tn + fp)),
        "NPV": float(tn / (tn + fn)),
        "test_AUC": float(roc_auc_score(y_true, y_prob)),
    }


def cv_auc_sklearn(
    factory: Callable[[], object],
    features: np.ndarray,
    labels: np.ndarray,
) -> tuple[float, float]:
    splitter = StratifiedKFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=SEED,
    )
    values = []
    for train_idx, test_idx in splitter.split(features, labels):
        model = factory()
        model.fit(features[train_idx], labels[train_idx])
        probability = model.predict_proba(features[test_idx])[:, 1]
        values.append(roc_auc_score(labels[test_idx], probability))
    return float(np.mean(values)), float(np.std(values, ddof=1))


def make_rf() -> RandomForestClassifier:
    return RandomForestClassifier(
        n_estimators=50,
        max_depth=90,
        max_features="sqrt",
        min_samples_split=2,
        min_samples_leaf=1,
        random_state=SEED,
        n_jobs=-1,
    )


def make_svm() -> SVC:
    return SVC(
        C=8,
        gamma=0.001,
        kernel="rbf",
        probability=True,
        random_state=SEED,
    )


def main() -> None:
    args = parse_args()
    root = args.repo_root.resolve()
    active_path = resolve(root, args.actives)
    decoy_path = resolve(root, args.decoys)
    output_path = resolve(root, args.output)
    metadata_path = resolve(root, args.metadata)

    for path in (active_path, decoy_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    actives = pd.read_csv(active_path)
    decoys = pd.read_csv(decoy_path)

    for name, frame in (("actives", actives), ("decoys", decoys)):
        if "std_smiles" not in frame.columns:
            raise ValueError(f"{name} file has no std_smiles column")
        if frame["std_smiles"].isna().any():
            raise ValueError(f"{name} file contains missing SMILES")

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=FP_SIZE,
    )
    active_fp = np.vstack(
        [smiles_to_fp(value, generator) for value in actives["std_smiles"]]
    )
    decoy_fp = np.vstack(
        [smiles_to_fp(value, generator) for value in decoys["std_smiles"]]
    )

    features = np.vstack([active_fp, decoy_fp])
    labels = np.concatenate(
        [
            np.ones(len(active_fp), dtype=np.uint8),
            np.zeros(len(decoy_fp), dtype=np.uint8),
        ]
    )

    train_x, test_x, train_y, test_y = train_test_split(
        features,
        labels,
        test_size=0.20,
        random_state=SEED,
        stratify=labels,
    )

    rows: list[dict] = []

    for name, factory in (("RF", make_rf), ("SVM", make_svm)):
        model = factory()
        model.fit(train_x, train_y)
        predicted = model.predict(test_x)
        probability = model.predict_proba(test_x)[:, 1]
        current = metrics(test_y, predicted, probability)
        cv_mean, cv_sd = cv_auc_sklearn(factory, features, labels)
        rows.append(
            {
                "model": name,
                **current,
                "CV_AUC_mean": cv_mean,
                "CV_AUC_sd": cv_sd,
                "backend": "scikit-learn",
                "task": "ligand_vs_random_unmatched_decoy",
            }
        )

    # Fail explicitly rather than silently substituting sklearn MLP.
    try:
        import tensorflow as tf
        from tensorflow.keras.layers import Dense, Dropout, Input
        from tensorflow.keras.models import Sequential
    except ImportError as exc:
        raise RuntimeError(
            "TensorFlow is required for the Keras FFNN. "
            "No sklearn MLP fallback is permitted."
        ) from exc

    random.seed(SEED)
    np.random.seed(SEED)
    tf.keras.utils.set_random_seed(SEED)
    try:
        tf.config.experimental.enable_op_determinism()
    except Exception:
        # Some TensorFlow versions/platforms do not expose this function.
        pass

    def make_ffnn():
        model = Sequential(
            [
                Input(shape=(FP_SIZE,)),
                Dense(750, activation="relu"),
                Dropout(0.6),
                Dense(1, activation="sigmoid"),
            ]
        )
        model.compile(optimizer="adam", loss="binary_crossentropy")
        return model

    ffnn = make_ffnn()
    ffnn.fit(train_x, train_y, epochs=30, batch_size=16, verbose=0)
    probability = ffnn.predict(test_x, verbose=0).ravel()
    predicted = (probability >= 0.5).astype(np.uint8)
    current = metrics(test_y, predicted, probability)

    splitter = StratifiedKFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=SEED,
    )
    cv_values = []
    for fold, (train_idx, test_idx) in enumerate(
        splitter.split(features, labels),
        start=1,
    ):
        tf.keras.backend.clear_session()
        tf.keras.utils.set_random_seed(SEED + fold)
        model = make_ffnn()
        model.fit(
            features[train_idx],
            labels[train_idx],
            epochs=30,
            batch_size=16,
            verbose=0,
        )
        fold_probability = model.predict(
            features[test_idx],
            verbose=0,
        ).ravel()
        cv_values.append(
            roc_auc_score(labels[test_idx], fold_probability)
        )

    rows.append(
        {
            "model": "Keras-FFNN",
            **current,
            "CV_AUC_mean": float(np.mean(cv_values)),
            "CV_AUC_sd": float(np.std(cv_values, ddof=1)),
            "backend": f"tensorflow-{tf.__version__}",
            "task": "ligand_vs_random_unmatched_decoy",
        }
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False, lineterminator="\n")

    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "analysis_name": (
            "Auxiliary reproduction of ligand-versus-random-unmatched-decoy "
            "classification"
        ),
        "not_an_exact_replication": True,
        "seed": SEED,
        "n_splits": N_SPLITS,
        "n_actives": int(len(active_fp)),
        "n_decoys": int(len(decoy_fp)),
        "fingerprint": {
            "type": "Morgan ECFP4",
            "radius": 2,
            "size": FP_SIZE,
        },
        "inputs": {
            "actives": str(active_path.relative_to(root)),
            "decoys": str(decoy_path.relative_to(root)),
        },
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "rdkit": rdBase.rdkitVersion,
            "scikit_learn": sklearn_version,
            "tensorflow": tf.__version__,
        },
    }
    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(result.to_string(index=False))
    print(f"\nSaved: {output_path}")
    print(f"Metadata: {metadata_path}")


if __name__ == "__main__":
    main()
