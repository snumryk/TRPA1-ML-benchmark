import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import (mean_squared_error, mean_absolute_error, r2_score,
                             roc_auc_score, matthews_corrcoef, balanced_accuracy_score)
from scipy.stats import spearmanr

RANDOM_SEED = 42
THRESHOLD = 7.0
np.random.seed(RANDOM_SEED)

# ── 1. Load data ───────────────────────────────────────────────
df = pd.read_csv('trpa1_antagonists.csv')
print(f"Loaded {len(df)} antagonist compounds")

# ── 2. Morgan fingerprints (ECFP4, radius=2, 2048 bits) ───────
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def smiles_to_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = mfpgen.GetFingerprint(mol)
    arr = np.zeros((2048,), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

print("Computing Morgan fingerprints...")
df['fp'] = df['std_smiles'].apply(smiles_to_fp)
df = df[df['fp'].notna()].copy()
print(f"After fingerprinting: {len(df)} compounds")

# ── 3. Train/val/test split (use existing scaffold split) ─────
def get_xy(subset):
    X = np.vstack(subset['fp'].values)
    y = subset['pchembl_median'].values
    return X, y

train = df[df['split'] == 'train']
val   = df[df['split'] == 'val']
test  = df[df['split'] == 'test']

X_train, y_train = get_xy(train)
X_val,   y_val   = get_xy(val)
X_test,  y_test  = get_xy(test)

print(f"\nSplit sizes: train={len(train)}, val={len(val)}, test={len(test)}")

# ── 4. REGRESSION ──────────────────────────────────────────────
print("\n" + "="*50)
print("REGRESSION (predicting pchembl_median)")
print("="*50)

rf_reg = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=RANDOM_SEED)
rf_reg.fit(X_train, y_train)

y_pred = rf_reg.predict(X_test)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
mae  = mean_absolute_error(y_test, y_pred)
r2   = r2_score(y_test, y_pred)
rho  = spearmanr(y_test, y_pred).correlation

print(f"  RMSE:     {rmse:.3f}")
print(f"  MAE:      {mae:.3f}")
print(f"  R2:       {r2:.3f}")
print(f"  Spearman: {rho:.3f}")

# Bootstrap CI for R2
boot_r2 = []
rng = np.random.default_rng(RANDOM_SEED)
for _ in range(1000):
    idx = rng.integers(0, len(y_test), len(y_test))
    boot_r2.append(r2_score(y_test[idx], y_pred[idx]))
ci_low, ci_high = np.percentile(boot_r2, [2.5, 97.5])
print(f"  R2 95% CI: [{ci_low:.3f}, {ci_high:.3f}]")

# ── 5. CLASSIFICATION (threshold 7.0) ──────────────────────────
print("\n" + "="*50)
print(f"CLASSIFICATION (active = pchembl >= {THRESHOLD})")
print("="*50)

y_train_cls = (y_train >= THRESHOLD).astype(int)
y_test_cls  = (y_test  >= THRESHOLD).astype(int)

print(f"  Train balance: {y_train_cls.mean()*100:.1f}% active")
print(f"  Test balance:  {y_test_cls.mean()*100:.1f}% active")

rf_cls = RandomForestClassifier(n_estimators=500, n_jobs=-1, random_state=RANDOM_SEED,
                                class_weight='balanced')
rf_cls.fit(X_train, y_train_cls)

y_prob = rf_cls.predict_proba(X_test)[:, 1]
y_pred_cls = rf_cls.predict(X_test)

auc = roc_auc_score(y_test_cls, y_prob)
mcc = matthews_corrcoef(y_test_cls, y_pred_cls)
bacc = balanced_accuracy_score(y_test_cls, y_pred_cls)

print(f"  AUC-ROC:           {auc:.3f}")
print(f"  MCC:               {mcc:.3f}")
print(f"  Balanced accuracy: {bacc:.3f}")

# Bootstrap CI for AUC
boot_auc = []
for _ in range(1000):
    idx = rng.integers(0, len(y_test_cls), len(y_test_cls))
    if len(np.unique(y_test_cls[idx])) < 2:
        continue
    boot_auc.append(roc_auc_score(y_test_cls[idx], y_prob[idx]))
ci_low, ci_high = np.percentile(boot_auc, [2.5, 97.5])
print(f"  AUC 95% CI:        [{ci_low:.3f}, {ci_high:.3f}]")

print("\n" + "="*50)
print("INTERPRETATION")
print("="*50)
print(f"  Regression R2 > 0.4  → data is learnable")
print(f"  Classification AUC > 0.75 → model finds real signal")
print(f"  This is your BASELINE. All modern models compare to this.")