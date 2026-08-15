ЩО ЦЕЙ СКРИПТ РОБИТЬ
============================
Чесно порівнює 4 способи опису молекули, прогоняючи КОЖЕН з них
ОБОМА алгоритмами на ТИХ САМИХ розбиттях даних. Виходить таблиця 4х2.

Додатково робить кілька НЕЗАЛЕЖНИХ розбиттів (partitions), бо одне
розбиття може випадково виявитись зручним або незручним.

Зберігає прогноз для кожної молекули - він потрібен наступним кроком,
щоб перевірити чи різниця між методами справжня чи це шум.

ЩО ПЕРЕВІРЯЄ ПЕРЕД РОБОТОЮ (fail-fast)
======================================
  * що рядки CSV і NPZ описують ТІ САМІ молекули (за y + scaffold);
  * що всі SMILES читаються (жодного мовчазного пропуску);
  * що немає NaN/inf, форми масивів правильні, InChIKey унікальні;
  * що в жодному розбитті каркас не потрапив і в навчання, і в перевірку;
  * що кожна молекула рівно один раз побувала в перевірці.
Якщо щось не так - скрипт зупиняється з поясненням, а не рахує тихо далі.


%pip install -q rdkit

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from xgboost import XGBRegressor

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

# Налаштування

RDLogger.DisableLog("rdApp.*")

from google.colab import drive
drive.mount("/content/drive")

from pathlib import Path

MODEL_SEED = 42
PARTITION_SEED = 1000
N_SPLITS = 5
N_PARTITIONS = 3

DATA_DIR = Path("/content/drive/MyDrive/trpa1_project")

CSV_FILE = DATA_DIR / "trpa1_antagonists.csv"
NPZ_FILE = DATA_DIR / "embeddings_all.npz"

OUT_PREFIX = "grid_final_v1"

# Допоміжне

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(msg: str) -> None:
    print(f"\n!!! ЗУПИНКА: {msg}\n")
    sys.exit(1)

# 1. Аргументи і файли

DATA_DIR = Path("/content/drive/MyDrive/trpa1_project").resolve()

CSV_FILE = DATA_DIR / "trpa1_antagonists.csv"
NPZ_FILE = DATA_DIR / "embeddings_all.npz"

OUT_PREFIX = "grid_final"

# Перевірка папки та вхідних файлів
if not DATA_DIR.is_dir():
    raise FileNotFoundError(f"Папки не існує: {DATA_DIR}")

for file_path in (CSV_FILE, NPZ_FILE):
    if not file_path.is_file():
        raise FileNotFoundError(f"Не знайдено файл: {file_path}")

# Унікальні імена вихідних файлів
STAMP = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")

OUT_GRID = DATA_DIR / f"{OUT_PREFIX}_results_{STAMP}.csv"
OUT_OOF = DATA_DIR / f"{OUT_PREFIX}_oof_{STAMP}.csv"
OUT_FOLDS = DATA_DIR / f"{OUT_PREFIX}_fold_assignments_{STAMP}.csv"
OUT_META = DATA_DIR / f"{OUT_PREFIX}_metadata_{STAMP}.json"

print("=" * 78)
print("СІТКА: 4 ОПИСИ МОЛЕКУЛИ × 2 АЛГОРИТМИ")
print("=" * 78)
print(f"Папка даних: {DATA_DIR}")
print(f"CSV: {CSV_FILE.name}   sha256={sha256(CSV_FILE)[:16]}...")
print(f"NPZ: {NPZ_FILE.name}   sha256={sha256(NPZ_FILE)[:16]}...")
print("\nВихідні файли:")
print(f"  {OUT_GRID.name}")
print(f"  {OUT_OOF.name}")
print(f"  {OUT_FOLDS.name}")
print(f"  {OUT_META.name}")

Завантаження та перевірка

df = pd.read_csv(CSV_FILE)
data = np.load(NPZ_FILE, allow_pickle=True)

REQUIRED_COLS = {"inchikey", "std_smiles", "pchembl_median", "scaffold"}
missing = REQUIRED_COLS - set(df.columns)
if missing:
    fail(f"у CSV бракує колонок: {sorted(missing)}")

REQUIRED_KEYS = {"rdkit", "cb_cls", "mf_mean", "y", "scaffold"}
missing = REQUIRED_KEYS - set(data.files)
if missing:
    fail(f"у NPZ бракує масивів: {sorted(missing)}")

y = np.asarray(data["y"], dtype=float).reshape(-1)
n = len(y)

if len(df) != n:
    fail(f"у CSV {len(df)} рядків, у NPZ {n} значень")

if df["inchikey"].duplicated().any():
    fail("у CSV є повторні InChIKey - набір не дедуплікований")

if not np.isfinite(y).all():
    fail("у цільовій змінній є NaN або inf")

# ПЕРЕВІРКА ВІДПОВІДНОСТІ РЯДКІВ
Найнебезпечніша помилка: embeddings від однієї молекули, структура
від іншої. Перевіряємо ДВІ незалежні колонки. Збіг обох на всіх
рядках практично виключає перемішування.

print("\nПеревірка відповідності рядків CSV <-> NPZ:")

y_csv = df["pchembl_median"].to_numpy(dtype=float)
if not np.allclose(y_csv, y, rtol=0, atol=1e-9):
    bad = int(np.flatnonzero(~np.isclose(y_csv, y, rtol=0, atol=1e-9))[0])
    fail(f"значення pchembl не збігаються, перша розбіжність у рядку {bad}")
print(f"  [OK] pchembl_median збігається для всіх {n} рядків")

scaf_npz = np.asarray(data["scaffold"]).astype(str).reshape(-1)
scaf_csv = df["scaffold"].astype(str).to_numpy()
if len(scaf_npz) != n:
    fail(f"у NPZ {len(scaf_npz)} каркасів на {n} молекул")
if not np.array_equal(scaf_csv, scaf_npz):
    bad = int(np.flatnonzero(scaf_csv != scaf_npz)[0])
    fail(f"каркаси не збігаються, перша розбіжність у рядку {bad}:\n"
         f"       CSV='{scaf_csv[bad]}'\n       NPZ='{scaf_npz[bad]}'")
print(f"  [OK] scaffold збігається для всіх {n} рядків (рядкове порівняння)")

if "inchikey" in data.files:
    ik_npz = np.asarray(data["inchikey"]).astype(str).reshape(-1)
    if not np.array_equal(df["inchikey"].astype(str).to_numpy(), ik_npz):
        fail("InChIKey у NPZ не збігається з CSV")
    print("  [OK] inchikey збігається (найсильніша перевірка)")
else:
    print("  [i]  inchikey у NPZ відсутній - при наступному витягуванні")
    print("       embeddings варто його зберігати")

scaffolds = np.asarray(scaf_csv, dtype="<U").astype(str)
if pd.isna(scaffolds).any():
    fail("є відсутні каркаси")
if any(str(s).strip() in ("", "nan", "None") for s in scaffolds):
    fail("є порожні каркаси")
print(f"\nМолекул: {n}   Каркасів: {len(np.unique(scaffolds))}")


# 3. Представлення

EXPECTED_DIMS = {"RDKit-15": 15, "ChemBERTa-CLS": 384, "MolFormer-Mean": 768}
NPZ_KEY = {"RDKit-15": "rdkit", "ChemBERTa-CLS": "cb_cls", "MolFormer-Mean": "mf_mean"}

REPRESENTATIONS: dict[str, np.ndarray] = {}
for name, key in NPZ_KEY.items():
    X = np.asarray(data[key], dtype=np.float32)
    if X.shape != (n, EXPECTED_DIMS[name]):
        fail(f"{name}: очікувалась форма {(n, EXPECTED_DIMS[name])}, отримано {X.shape}")
    if not np.isfinite(X).all():
        fail(f"{name}: знайдено NaN або inf")
    REPRESENTATIONS[name] = X

# Morgan рахуємо тут; будь-який нечитабельний SMILES = зупинка
print("\nРахую Morgan fingerprints (achiral ECFP4, 2048 біт)...")
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def morgan(smi, idx):
    if not isinstance(smi, str) or not smi.strip():
        fail(f"порожній SMILES у рядку {idx}")
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        fail(f"RDKit не зміг розібрати SMILES у рядку {idx}: {smi}")
    return np.asarray(mfpgen.GetFingerprint(mol), dtype=np.uint8)

X_morgan = np.vstack([morgan(s, i) for i, s in enumerate(df["std_smiles"])])
empty = int((X_morgan.sum(axis=1) == 0).sum())
if empty:
    fail(f"{empty} молекул дали порожній fingerprint")
print(f"  [OK] Morgan: {X_morgan.shape}, порожніх немає")

REPRESENTATIONS = {"Morgan-2048": X_morgan, **REPRESENTATIONS}

# 4. Алгоритми

def make_rf():
    return RandomForestRegressor(n_estimators=500, n_jobs=-1,
                                 random_state=MODEL_SEED)

def make_xgb():
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

ALGORITHMS = {"RF": make_rf, "XGB": make_xgb}

# 5. Розбиття: кілька незалежних partitions
GroupKFold детермінований - він завжди дає одне й те саме розбиття.
Щоб отримати інші, перемішуємо каркаси і розкладаємо їх по частинах жадібно (найбільші спершу), балансуючи розміри.

def make_partition(scaffolds, n_splits, seed=None):
    uniq, counts = np.unique(scaffolds, return_counts=True)

    if len(uniq) < n_splits:
        raise ValueError(
            f"Недостатньо унікальних каркасів: {len(uniq)} "
            f"для {n_splits} folds"
        )

    if seed is None:
        # Partition 0: стандартний детермінований GroupKFold
        gkf = GroupKFold(n_splits=n_splits)

        folds = np.full(len(scaffolds), -1, dtype=int)

        for fold_id, (_, test_idx) in enumerate(
            gkf.split(
                X=np.zeros(len(scaffolds)),
                groups=scaffolds,
            )
        ):
            folds[test_idx] = fold_id

        return folds

    # Додаткові partitions:
    # каркаси випадково впорядковуються всередині груп однакового розміру,
    # після чого жадібно розподіляються між folds для балансування
    rng = np.random.default_rng(seed)

    permutation = rng.permutation(len(uniq))
    uniq_permuted = uniq[permutation]
    counts_permuted = counts[permutation]

    order = np.argsort(-counts_permuted, kind="stable")

    fold_load = np.zeros(n_splits, dtype=int)
    scaffold_to_fold = {}

    for index in order:
        fold_id = int(np.argmin(fold_load))

        scaffold = uniq_permuted[index]
        scaffold_to_fold[scaffold] = fold_id
        fold_load[fold_id] += counts_permuted[index]

    return np.array(
        [scaffold_to_fold[scaffold] for scaffold in scaffolds],
        dtype=int,
    )


if N_PARTITIONS < 1:
    raise ValueError("N_PARTITIONS має бути щонайменше 1")

partitions = []

for partition_id in range(N_PARTITIONS):
    partition_seed = (
        None
        if partition_id == 0
        else PARTITION_SEED + partition_id
    )

    fold_assignments = make_partition(
        scaffolds=scaffolds,
        n_splits=N_SPLITS,
        seed=partition_seed,
    )

    partitions.append(fold_assignments)


print(
    f"\nРозбиттів: {N_PARTITIONS}, "
    f"у кожному {N_SPLITS} folds"
)
print(
    "  partition 0 = стандартний GroupKFold; "
    "partitions 1+ = додаткові перемішані scaffold partitions"
)


# Перевірка кожного partition
for partition_id, fold_assignments in enumerate(partitions):
    test_count = np.zeros(n, dtype=int)

    for fold_id in range(N_SPLITS):
        test_idx = np.flatnonzero(fold_assignments == fold_id)
        train_idx = np.flatnonzero(fold_assignments != fold_id)

        if len(test_idx) == 0:
            raise RuntimeError(
                f"Partition {partition_id}, fold {fold_id}: "
                "порожня test-вибірка"
            )

        train_scaffolds = set(scaffolds[train_idx])
        test_scaffolds = set(scaffolds[test_idx])
        overlap = train_scaffolds & test_scaffolds

        if overlap:
            raise RuntimeError(
                f"Partition {partition_id}, fold {fold_id}: "
                f"{len(overlap)} каркасів є одночасно "
                "у train і test"
            )

        test_count[test_idx] += 1

    if not np.all(test_count == 1):
        bad_count = int(np.sum(test_count != 1))
        raise RuntimeError(
            f"Partition {partition_id}: {bad_count} молекул "
            "не потрапили в test рівно один раз"
        )

    fold_sizes = [
        int(np.sum(fold_assignments == fold_id))
        for fold_id in range(N_SPLITS)
    ]

    print(
        f"  partition {partition_id}: "
        f"розміри test folds {fold_sizes} — витоку немає"
    )

# 6. Прогін

def metrics(y_true, y_pred):
    if np.ptp(y_pred) == 0:
        rho = float("nan")
    else:
        rho = float(spearmanr(y_true, y_pred).correlation)

    return {
        "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "R2": float(r2_score(y_true, y_pred)),
        "Spearman": rho,
    }


# Фіксовані checkpoint-файли.
# Вони навмисно без timestamp, щоб запуск можна було продовжити.
CHECKPOINT_GRID = DATA_DIR / f"{OUT_PREFIX}_checkpoint_results.csv"
CHECKPOINT_OOF = DATA_DIR / f"{OUT_PREFIX}_checkpoint_oof.csv"

COMBOS = [
    (rep_name, alg_name)
    for rep_name in REPRESENTATIONS
    for alg_name in ALGORITHMS
]
COMBOS.append(("__dummy__", "Mean"))


# Завантажуємо вже завершені результати, якщо runtime раніше обірвався.
if CHECKPOINT_GRID.exists():
    checkpoint_grid = pd.read_csv(CHECKPOINT_GRID)
    rows = checkpoint_grid.to_dict("records")

    completed = set(
        zip(
            checkpoint_grid["representation"],
            checkpoint_grid["algorithm"],
            checkpoint_grid["partition"],
        )
    )

    print(
        f"Знайдено checkpoint: "
        f"{len(completed)} завершених combo-partition запусків"
    )
else:
    rows = []
    completed = set()


if CHECKPOINT_OOF.exists():
    checkpoint_oof = pd.read_csv(CHECKPOINT_OOF)
    oof_records = checkpoint_oof.to_dict("records")
else:
    oof_records = []


total_jobs = len(COMBOS) * N_PARTITIONS
finished_jobs = len(completed)

print("=" * 78)
print(
    f"ЗАПУСК: {len(COMBOS)} комбінацій × "
    f"{N_PARTITIONS} partitions × {N_SPLITS} folds"
)
print(f"Вже завершено: {finished_jobs} із {total_jobs} combo-partition jobs")
print("=" * 78)


for combo_index, (rep_name, alg_name) in enumerate(COMBOS, start=1):
    is_dummy = rep_name == "__dummy__"

    label = (
        "Baseline (середнє)"
        if is_dummy
        else f"{alg_name} + {rep_name}"
    )

    X = X_morgan if is_dummy else REPRESENTATIONS[rep_name]

    print(f"\n[{combo_index}/{len(COMBOS)}] {label}")

    for partition_id, fold_assignments in enumerate(partitions):
        job_key = (rep_name, alg_name, partition_id)

        if job_key in completed:
            print(f"  partition {partition_id}: вже готовий, пропускаю")
            continue

        partition_start = time.perf_counter()

        oof = np.full(n, np.nan, dtype=float)
        fold_metrics = []

        for fold_id in range(N_SPLITS):
            fold_start = time.perf_counter()

            test_idx = np.flatnonzero(fold_assignments == fold_id)
            train_idx = np.flatnonzero(fold_assignments != fold_id)

            if is_dummy:
                model = DummyRegressor(strategy="mean")
            else:
                model = ALGORITHMS[alg_name]()

            model.fit(X[train_idx], y[train_idx])
            predictions = model.predict(X[test_idx])

            if not np.isfinite(predictions).all():
                raise RuntimeError(
                    f"{label}, partition {partition_id}, fold {fold_id}: "
                    "модель повернула NaN або inf"
                )

            oof[test_idx] = predictions
            fold_result = metrics(y[test_idx], predictions)
            fold_metrics.append(fold_result)

            fold_seconds = time.perf_counter() - fold_start

            print(
                f"    fold {fold_id + 1}/{N_SPLITS}: "
                f"RMSE={fold_result['RMSE']:.3f}, "
                f"R2={fold_result['R2']:.3f}, "
                f"{fold_seconds:.1f} с"
            )

            del model, predictions
            gc.collect()

        if np.isnan(oof).any():
            raise RuntimeError(
                f"{label}, partition {partition_id}: "
                "не всі OOF-прогнози заповнені"
            )

        pooled = metrics(y, oof)

        result_row = {
            "label": label,
            "representation": rep_name,
            "algorithm": alg_name,
            "partition": partition_id,

            "RMSE_pooled": pooled["RMSE"],
            "R2_pooled": pooled["R2"],
            "Spearman_pooled": pooled["Spearman"],

            "RMSE_foldmean": float(
                np.mean([m["RMSE"] for m in fold_metrics])
            ),
            "RMSE_foldsd": float(
                np.std([m["RMSE"] for m in fold_metrics], ddof=1)
            ),

            "R2_foldmean": float(
                np.mean([m["R2"] for m in fold_metrics])
            ),
            "R2_foldsd": float(
                np.std([m["R2"] for m in fold_metrics], ddof=1)
            ),

            "Spearman_foldmean": float(
                np.nanmean([m["Spearman"] for m in fold_metrics])
            ),
            "Spearman_foldsd": float(
                np.nanstd(
                    [m["Spearman"] for m in fold_metrics],
                    ddof=1,
                )
            ),

            "R2_per_fold": ";".join(
                f"{m['R2']:.4f}" for m in fold_metrics
            ),
        }

        rows.append(result_row)

        for molecule_index in range(n):
            oof_records.append({
                "inchikey": df["inchikey"].iat[molecule_index],
                "partition": partition_id,
                "fold": int(fold_assignments[molecule_index]),
                "representation": rep_name,
                "algorithm": alg_name,
                "label": label,
                "y_true": float(y[molecule_index]),
                "y_pred": float(oof[molecule_index]),
            })

        # Checkpoint одразу після кожного partition.
        pd.DataFrame(rows).to_csv(CHECKPOINT_GRID, index=False)
        pd.DataFrame(oof_records).to_csv(CHECKPOINT_OOF, index=False)

        completed.add(job_key)

        partition_seconds = time.perf_counter() - partition_start

        print(
            f"  partition {partition_id}: ГОТОВО — "
            f"RMSE={pooled['RMSE']:.3f}, "
            f"R2={pooled['R2']:.3f}, "
            f"Spearman={pooled['Spearman']:.3f}; "
            f"{partition_seconds / 60:.1f} хв"
        )


grid = pd.DataFrame(rows)

print("\n" + "=" * 78)
print("УСІ КОМБІНАЦІЇ ЗАВЕРШЕНО")
print(f"Checkpoint results: {CHECKPOINT_GRID}")
print(f"Checkpoint OOF:     {CHECKPOINT_OOF}")
print("=" * 78)

# 7. Таблиця

print("\n" + "=" * 78)
print("ГОЛОВНА ТАБЛИЦЯ")
print("(RMSE - головна метрика, менше = краще; усереднено по partitions)")
print("=" * 78)
print(f"{'Опис молекули':<18} {'RF: RMSE':>12} {'R2':>7} | {'XGB: RMSE':>12} {'R2':>7}")
print("-" * 78)
for rep in REPRESENTATIONS:
    cells = []
    for alg in ALGORITHMS:
        s = grid[(grid.representation == rep) & (grid.algorithm == alg)]
        cells.append((s["RMSE_pooled"].mean(), s["R2_pooled"].mean()))
    print(f"{rep:<18} {cells[0][0]:>12.3f} {cells[0][1]:>7.3f} | "
          f"{cells[1][0]:>12.3f} {cells[1][1]:>7.3f}")
b = grid[grid.representation == "__dummy__"]
print("-" * 78)
print(f"{'Baseline (середнє)':<18} {b['RMSE_pooled'].mean():>12.3f} "
      f"{b['R2_pooled'].mean():>7.3f}")
print("-" * 78)

if N_PARTITIONS > 1:
    print("\nСтабільність між розбиттями (RMSE_pooled: мін - макс):")
    for rep in REPRESENTATIONS:
        for alg in ALGORITHMS:
            s = grid[(grid.representation == rep) & (grid.algorithm == alg)]
            print(f"  {alg} + {rep:<18} {s['RMSE_pooled'].min():.3f} - "
                  f"{s['RMSE_pooled'].max():.3f}")

print("\n[!] Середні по рядках/стовпцях свідомо НЕ друкуються:")
print("    два алгоритми - замало щоб усереднення вважати ефектом опису.")
print("    Порівнювати треба попарно при однаковому алгоритмі -")
print("    це наступний крок (статистична перевірка).")

# 8. Збереження

grid.to_csv(OUT_GRID, index=False)
pd.DataFrame(oof_records).to_csv(OUT_OOF, index=False)

fold_rows = []
for p, folds in enumerate(partitions):
    for i in range(n):
        fold_rows.append({"inchikey": df["inchikey"].iat[i],
                          "scaffold": scaffolds[i],
                          "partition": p, "fold": int(folds[i])})
pd.DataFrame(fold_rows).to_csv(OUT_FOLDS, index=False)

import sklearn, scipy, xgboost, rdkit
meta = {
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "csv_file": str(CSV_FILE), "csv_sha256": sha256(CSV_FILE),
    "npz_file": str(NPZ_FILE), "npz_sha256": sha256(NPZ_FILE),
    "n_compounds": int(n),
    "n_scaffolds": int(len(np.unique(scaffolds))),
    "n_partitions": N_PARTITIONS,
    "n_splits": N_SPLITS,
    "model_seed": MODEL_SEED,
    "partition_seed_base": PARTITION_SEED,
    "note_on_seeds": ("model_seed впливає лише на RF/XGB; partition 0 - "
                      "детермінований GroupKFold, partitions 1+ згенеровані "
                      "з partition_seed_base"),
    "primary_metric": "RMSE_pooled (на всіх out-of-fold прогнозах)",
    "secondary_metrics": ["R2_pooled", "Spearman_pooled", "fold-wise mean/sd (ddof=1)"],
    "representations": {k: int(v.shape[1]) for k, v in REPRESENTATIONS.items()},
    "morgan": {"radius": 2, "fpSize": 2048, "includeChirality": False},
    "rf_params": make_rf().get_params(),
    "xgb_params": make_xgb().get_params(),
    "leakage_checks_passed": True,
    "versions": {
        "python": platform.python_version(), "numpy": np.__version__,
        "pandas": pd.__version__, "sklearn": sklearn.__version__,
        "scipy": scipy.__version__, "xgboost": xgboost.__version__,
        "rdkit": rdkit.__version__,
    },
}
OUT_META.write_text(json.dumps(meta, indent=2, ensure_ascii=False, default=str),
                    encoding="utf-8")

print("\n" + "=" * 78)
print("ГОТОВО")
print("=" * 78)
for f in (OUT_GRID, OUT_OOF, OUT_FOLDS, OUT_META):
    print(f"  {f.name}")