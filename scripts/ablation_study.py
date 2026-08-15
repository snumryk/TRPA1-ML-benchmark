"""
Ablation study: do RDKit descriptors alone explain D-MPNN+RDKit performance?
If yes -> GNN adds nothing, improvement came from descriptors.
If no  -> GNN genuinely contributes structural signal.
"""
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error, roc_auc_score, matthews_corrcoef
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

SEED = 42
THRESHOLD = 7.0

# ── 1. Same 15 RDKit descriptors as D-MPNN+RDKit ─────────────
RDKIT_DESCS = [
    'MolWt', 'MolLogP', 'MolMR', 'TPSA',
    'NumHAcceptors', 'NumHDonors', 'NumRotatableBonds',
    'NumAromaticRings', 'RingCount', 'FractionCSP3',
    'HeavyAtomCount', 'NumAliphaticRings', 'NumSaturatedRings',
    'NumHeteroatoms', 'LabuteASA',
]

def compute_features(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    feats = []
    for name in RDKIT_DESCS:
        func = getattr(Descriptors, name, None)
        try:
            feats.append(float(func(mol)))
        except Exception:
            feats.append(0.0)
    return np.array(feats, dtype=np.float32)

# ── 2. Load data, same scaffold split ─────────────────────────
df = pd.read_csv('trpa1_antagonists.csv')
df['rdkit_feats'] = df['std_smiles'].apply(compute_features)
df = df[df['rdkit_feats'].notna()].copy()

def get_xy(subset):
    X = np.vstack(subset['rdkit_feats'].values)
    y = subset['pchembl_median'].values
    return X, y

X_train, y_train = get_xy(df[df['split'] == 'train'])
X_val, y_val     = get_xy(df[df['split'] == 'val'])
X_test, y_test   = get_xy(df[df['split'] == 'test'])

# Scale features (important for MLP, not for tree models)
scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_val_s   = scaler.transform(X_val)
X_test_s  = scaler.transform(X_test)

print(f"Features: {len(RDKIT_DESCS)} RDKit descriptors")
print(f"Train: {len(X_train)}, Val: {len(X_val)}, Test: {len(X_test)}")

# ── 3. Models on RDKit descriptors only ───────────────────────
models = {
    'RF (15 RDKit)': RandomForestRegressor(
        n_estimators=500, n_jobs=-1, random_state=SEED),
    'XGBoost (15 RDKit)': XGBRegressor(
        n_estimators=500, max_depth=6, learning_rate=0.05,
        n_jobs=-1, random_state=SEED),
    'MLP (15 RDKit)': MLPRegressor(
        hidden_layer_sizes=(128, 64), activation='relu',
        max_iter=500, early_stopping=True, validation_fraction=0.1,
        random_state=SEED),
}

def evaluate(y_true, y_pred):
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    rho = spearmanr(y_true, y_pred).correlation
    y_cls = (y_true >= THRESHOLD).astype(int)
    p_cls = (y_pred >= THRESHOLD).astype(int)
    auc = roc_auc_score(y_cls, y_pred)
    mcc = matthews_corrcoef(y_cls, p_cls)
    return rmse, r2, rho, auc, mcc

print(f"\n{'='*75}")
print("ABLATION: RDKit descriptors only (15 features) vs full models")
print(f"{'='*75}")
print(f"{'Model':<25} {'RMSE':>7} {'R2':>7} {'Spearman':>10} {'AUC':>7} {'MCC':>7}")
print("-"*75)

for name, model in models.items():
    # MLP uses scaled features
    if 'MLP' in name:
        model.fit(X_train_s, y_train)
        preds = model.predict(X_test_s)
    else:
        model.fit(X_train, y_train)
        preds = model.predict(X_test)
    rmse, r2, rho, auc, mcc = evaluate(y_test, preds)
    print(f"{name:<25} {rmse:>7.3f} {r2:>7.3f} {rho:>10.3f} {auc:>7.3f} {mcc:>7.3f}")

# ── 4. Comparison with previous results ───────────────────────
print(f"\n{'='*75}")
print("FULL COMPARISON (ablation context)")
print(f"{'='*75}")
print(f"{'Model':<25} {'RMSE':>7} {'R2':>7} {'Spearman':>10} {'AUC':>7} {'Features'}")
print("-"*75)
print(f"{'RF (Morgan 2048)':<25} {'0.795':>7} {'0.232':>7} {'0.552':>10} {'0.795':>7}  Morgan ECFP4")
print(f"{'XGBoost (Morgan 2048)':<25} {'0.781':>7} {'0.257':>7} {'0.547':>10} {'0.790':>7}  Morgan ECFP4")
print(f"{'D-MPNN (graph only)':<25} {'0.823':>7} {'0.176':>7} {'0.554':>10} {'0.761':>7}  Learned from graph")
print(f"{'D-MPNN+RDKit (ES)':<25} {'0.763':>7} {'0.292':>7} {'0.621':>10} {'0.810':>7}  Graph + 15 RDKit")

# Re-print ablation results
for name, model in models.items():
    if 'MLP' in name:
        preds = model.predict(X_test_s)
    else:
        preds = model.predict(X_test)
    rmse, r2, rho, auc, mcc = evaluate(y_test, preds)
    print(f"{name:<25} {rmse:>7.3f} {r2:>7.3f} {rho:>10.3f} {auc:>7.3f}  15 RDKit only")

print("-"*75)
print("\nINTERPRETATION:")
print("If 'XGBoost (15 RDKit)' R2 ~ 0.29 -> GNN adds nothing, descriptors explain all")
print("If 'XGBoost (15 RDKit)' R2 ~ 0.15 -> GNN genuinely adds structural signal")