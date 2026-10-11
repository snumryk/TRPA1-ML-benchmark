import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import (mean_squared_error, mean_absolute_error, r2_score,
                             roc_auc_score, matthews_corrcoef, balanced_accuracy_score)
from scipy.stats import spearmanr
from xgboost import XGBRegressor, XGBClassifier
from lightgbm import LGBMRegressor, LGBMClassifier

RANDOM_SEED = 42
THRESHOLD = 7.0
np.random.seed(RANDOM_SEED)

# ── Load + featurize ───────────────────────────────────────────
df = pd.read_csv('trpa1_antagonists.csv')
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def smiles_to_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = mfpgen.GetFingerprint(mol)
    arr = np.zeros((2048,), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

print("Computing fingerprints...")
df['fp'] = df['std_smiles'].apply(smiles_to_fp)
df = df[df['fp'].notna()].copy()

def get_xy(subset):
    return np.vstack(subset['fp'].values), subset['pchembl_median'].values

X_train, y_train = get_xy(df[df['split'] == 'train'])
X_test,  y_test  = get_xy(df[df['split'] == 'test'])

y_train_cls = (y_train >= THRESHOLD).astype(int)
y_test_cls  = (y_test  >= THRESHOLD).astype(int)

rng = np.random.default_rng(RANDOM_SEED)

def bootstrap_ci(metric_fn, y_true, y_pred, n=1000):
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(y_true), len(y_true))
        try:
            if len(np.unique(y_true[idx])) < 2 and metric_fn is roc_auc_score:
                continue
            vals.append(metric_fn(y_true[idx], y_pred[idx]))
        except Exception:
            continue
    return np.percentile(vals, [2.5, 97.5])

# ── Define models ──────────────────────────────────────────────
regressors = {
    'RF':       RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=RANDOM_SEED),
    'XGBoost':  XGBRegressor(n_estimators=500, max_depth=6, learning_rate=0.05,
                             subsample=0.8, colsample_bytree=0.8,
                             n_jobs=-1, random_state=RANDOM_SEED),
    'LightGBM': LGBMRegressor(n_estimators=500, max_depth=6, learning_rate=0.05,
                              subsample=0.8, colsample_bytree=0.8,
                              n_jobs=-1, random_state=RANDOM_SEED, verbose=-1),
}

classifiers = {
    'RF':       RandomForestClassifier(n_estimators=500, n_jobs=-1, random_state=RANDOM_SEED,
                                       class_weight='balanced'),
    'XGBoost':  XGBClassifier(n_estimators=500, max_depth=6, learning_rate=0.05,
                              subsample=0.8, colsample_bytree=0.8,
                              n_jobs=-1, random_state=RANDOM_SEED,
                              eval_metric='logloss'),
    'LightGBM': LGBMClassifier(n_estimators=500, max_depth=6, learning_rate=0.05,
                               subsample=0.8, colsample_bytree=0.8,
                               n_jobs=-1, random_state=RANDOM_SEED, verbose=-1),
}

# ── Regression ─────────────────────────────────────────────────
print("\n" + "="*70)
print("REGRESSION RESULTS")
print("="*70)
print(f"{'Model':<10} {'RMSE':>7} {'MAE':>7} {'R2':>7} {'Spearman':>10} {'R2 95% CI':>20}")
print("-"*70)

reg_results = {}
for name, model in regressors.items():
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, pred))
    mae  = mean_absolute_error(y_test, pred)
    r2   = r2_score(y_test, pred)
    rho  = spearmanr(y_test, pred).correlation
    ci   = bootstrap_ci(r2_score, y_test, pred)
    reg_results[name] = pred
    print(f"{name:<10} {rmse:>7.3f} {mae:>7.3f} {r2:>7.3f} {rho:>10.3f}   [{ci[0]:.3f}, {ci[1]:.3f}]")

# ── Classification ─────────────────────────────────────────────
print("\n" + "="*70)
print(f"CLASSIFICATION RESULTS (active = pchembl >= {THRESHOLD})")
print("="*70)
print(f"{'Model':<10} {'AUC':>7} {'MCC':>7} {'BalAcc':>8} {'AUC 95% CI':>20}")
print("-"*70)

for name, model in classifiers.items():
    model.fit(X_train, y_train_cls)
    prob = model.predict_proba(X_test)[:, 1]
    pred = model.predict(X_test)
    auc  = roc_auc_score(y_test_cls, prob)
    mcc  = matthews_corrcoef(y_test_cls, pred)
    bacc = balanced_accuracy_score(y_test_cls, pred)
    ci   = bootstrap_ci(roc_auc_score, y_test_cls, prob)
    print(f"{name:<10} {auc:>7.3f} {mcc:>7.3f} {bacc:>8.3f}   [{ci[0]:.3f}, {ci[1]:.3f}]")

print("\n" + "="*70)
print("Best classical baseline is the bar GNNs must beat.")
print("="*70)