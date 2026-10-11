"""
GroupKFold CV (groups = scaffolds) for Morgan fingerprints AND D-MPNN.
Brings these two model types onto the same statistical footing as the
transformer/RDKit models already done in Colab.

Two parts:
  PART A: Morgan ECFP4 (2048) with RF and XGBoost — fast
  PART B: D-MPNN (graph only) and D-MPNN+RDKit — slower (5 trainings each)
"""
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, Descriptors
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_squared_error, roc_auc_score, matthews_corrcoef
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

SEED = 42
THRESHOLD = 7.0

df = pd.read_csv('trpa1_antagonists.csv')
y = df['pchembl_median'].values
scaffolds = df['scaffold'].values
print(f"Compounds: {len(df)}, scaffolds: {len(np.unique(scaffolds))}")

def metrics(y_true, y_pred):
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    rho = spearmanr(y_true, y_pred).correlation
    y_cls = (y_true >= THRESHOLD).astype(int)
    if len(np.unique(y_cls)) < 2:
        return rmse, r2, rho, np.nan, np.nan
    auc = roc_auc_score(y_cls, y_pred)
    mcc = matthews_corrcoef(y_cls, (y_pred >= THRESHOLD).astype(int))
    return rmse, r2, rho, auc, mcc

def ms(vals):
    vals = np.array(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    return np.mean(vals), np.std(vals)

gkf = GroupKFold(n_splits=5)
all_results = {}

# ══════════════════════════════════════════════════════════════
# PART A: Morgan fingerprints (fast)
# ══════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("PART A: Morgan ECFP4 (2048 bits) — RF and XGBoost")
print("="*70)

mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def smi_to_morgan(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return np.zeros(2048, dtype=np.int8)
    fp = mfpgen.GetFingerprint(mol)
    arr = np.zeros((2048,), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

X_morgan = np.vstack([smi_to_morgan(s) for s in df['std_smiles']])
print(f"Morgan matrix: {X_morgan.shape}")

morgan_models = {
    'RF (Morgan)':  lambda: RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=SEED),
    'XGB (Morgan)': lambda: XGBRegressor(n_estimators=500, max_depth=6, learning_rate=0.05,
                                          n_jobs=-1, random_state=SEED),
}

for name, make_model in morgan_models.items():
    fm = {'rmse': [], 'r2': [], 'rho': [], 'auc': [], 'mcc': []}
    for tr, te in gkf.split(X_morgan, y, groups=scaffolds):
        m = make_model()
        m.fit(X_morgan[tr], y[tr])
        p = m.predict(X_morgan[te])
        for k, v in zip(['rmse','r2','rho','auc','mcc'], metrics(y[te], p)):
            fm[k].append(v)
    all_results[name] = {k: ms(fm[k]) for k in fm}
    r = all_results[name]
    print(f"  {name:<16} R2={r['r2'][0]:.3f}±{r['r2'][1]:.3f}  "
          f"RMSE={r['rmse'][0]:.3f}±{r['rmse'][1]:.3f}  "
          f"AUC={r['auc'][0]:.3f}±{r['auc'][1]:.3f}")

# ══════════════════════════════════════════════════════════════
# PART B: D-MPNN (slower — 5 trainings per config)
# ══════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("PART B: D-MPNN (graph only) and D-MPNN+RDKit")
print("="*70)
print("(5 trainings each with early stopping — ~15 min total on CPU)")

import torch
import lightning as L
from lightning.pytorch.callbacks import EarlyStopping
from chemprop.data import MoleculeDatapoint, MoleculeDataset, build_dataloader
from chemprop.nn import BondMessagePassing, MeanAggregation, RegressionFFN
from chemprop.models import MPNN

RDKIT_DESCS = [
    'MolWt', 'MolLogP', 'MolMR', 'TPSA',
    'NumHAcceptors', 'NumHDonors', 'NumRotatableBonds',
    'NumAromaticRings', 'RingCount', 'FractionCSP3',
    'HeavyAtomCount', 'NumAliphaticRings', 'NumSaturatedRings',
    'NumHeteroatoms', 'LabuteASA',
]

def compute_rdkit(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return np.zeros(len(RDKIT_DESCS), dtype=np.float32)
    return np.array([float(getattr(Descriptors, n)(mol)) for n in RDKIT_DESCS], dtype=np.float32)

smiles = df['std_smiles'].values

def run_dmpnn_cv(use_rdkit):
    """Run 5-fold GroupKFold for D-MPNN, optionally with RDKit extra features."""
    fm = {'rmse': [], 'r2': [], 'rho': [], 'auc': [], 'mcc': []}

    for fold, (tr, te) in enumerate(gkf.split(smiles, y, groups=scaffolds)):
        # Further split train into train/val for early stopping (80/20 of train)
        rng = np.random.default_rng(SEED + fold)
        tr_shuffled = rng.permutation(tr)
        n_val = int(0.15 * len(tr_shuffled))
        val_idx = tr_shuffled[:n_val]
        fit_idx = tr_shuffled[n_val:]

        def make_points(indices):
            pts = []
            for i in indices:
                mol = Chem.MolFromSmiles(smiles[i])
                if mol is None:
                    continue
                if use_rdkit:
                    x_d = compute_rdkit(smiles[i])
                    pts.append(MoleculeDatapoint(mol=mol, y=np.array([y[i]]), x_d=x_d))
                else:
                    pts.append(MoleculeDatapoint(mol=mol, y=np.array([y[i]])))
            return pts

        fit_dset = MoleculeDataset(make_points(fit_idx))
        val_dset = MoleculeDataset(make_points(val_idx))
        test_dset = MoleculeDataset(make_points(te))

        fit_loader = build_dataloader(fit_dset, batch_size=64, shuffle=True, seed=SEED)
        val_loader = build_dataloader(val_dset, batch_size=64, shuffle=False)
        test_loader = build_dataloader(test_dset, batch_size=64, shuffle=False)

        mp = BondMessagePassing(d_h=300, depth=3)
        agg = MeanAggregation()
        if use_rdkit:
            ffn = RegressionFFN(input_dim=300 + len(RDKIT_DESCS))
        else:
            ffn = RegressionFFN()
        model = MPNN(message_passing=mp, agg=agg, predictor=ffn,
                     warmup_epochs=2, init_lr=1e-4, max_lr=1e-3, final_lr=1e-4)

        early = EarlyStopping(monitor='val_loss', patience=10, mode='min')
        trainer = L.Trainer(max_epochs=100, accelerator='cpu',
                            enable_progress_bar=False, enable_model_summary=False,
                            logger=False, callbacks=[early])
        trainer.fit(model, fit_loader, val_loader)

        preds = torch.cat(trainer.predict(model, test_loader)).numpy().ravel()
        actuals = np.array([y[i] for i in te if Chem.MolFromSmiles(smiles[i]) is not None])
        for k, v in zip(['rmse','r2','rho','auc','mcc'], metrics(actuals, preds)):
            fm[k].append(v)
        print(f"    fold {fold+1}/5 done (R2={fm['r2'][-1]:.3f})")

    return {k: ms(fm[k]) for k in fm}

print("\nD-MPNN (graph only)...")
all_results['D-MPNN'] = run_dmpnn_cv(use_rdkit=False)
r = all_results['D-MPNN']
print(f"  D-MPNN: R2={r['r2'][0]:.3f}±{r['r2'][1]:.3f}  AUC={r['auc'][0]:.3f}±{r['auc'][1]:.3f}")

print("\nD-MPNN+RDKit...")
all_results['D-MPNN+RDKit'] = run_dmpnn_cv(use_rdkit=True)
r = all_results['D-MPNN+RDKit']
print(f"  D-MPNN+RDKit: R2={r['r2'][0]:.3f}±{r['r2'][1]:.3f}  AUC={r['auc'][0]:.3f}±{r['auc'][1]:.3f}")

# ══════════════════════════════════════════════════════════════
# FINAL: combined table with Colab results
# ══════════════════════════════════════════════════════════════
print("\n" + "="*95)
print("COMPLETE CV TABLE (all models, 5-fold GroupKFold by scaffold) — mean ± std")
print("="*95)
print(f"{'Model':<22} {'R2':>14} {'RMSE':>14} {'Spearman':>14} {'AUC':>14}")
print("-"*95)

# From Colab (transformer + RDKit models)
colab_results = {
    'XGB CB-CLS only':   {'r2': (0.539, 0.084), 'rmse': (0.633, 0.057), 'rho': (0.725, 0.064), 'auc': (0.866, 0.043)},
    'XGB MF-Mean only':  {'r2': (0.529, 0.085), 'rmse': (0.639, 0.042), 'rho': (0.703, 0.066), 'auc': (0.857, 0.045)},
    'XGB CB-CLS+RDKit':  {'r2': (0.536, 0.097), 'rmse': (0.634, 0.057), 'rho': (0.728, 0.073), 'auc': (0.869, 0.049)},
    'RF (RDKit-15)':     {'r2': (0.478, 0.163), 'rmse': (0.664, 0.053), 'rho': (0.695, 0.089), 'auc': (0.848, 0.051)},
}

# Merge: Colab first (best), then local (Morgan, D-MPNN)
merged = {**colab_results, **all_results}
# Sort by mean R2 descending
for name in sorted(merged, key=lambda k: merged[k]['r2'][0], reverse=True):
    r = merged[name]
    print(f"{name:<22} {r['r2'][0]:>7.3f}±{r['r2'][1]:.3f} "
          f"{r['rmse'][0]:>7.3f}±{r['rmse'][1]:.3f} "
          f"{r['rho'][0]:>7.3f}±{r['rho'][1]:.3f} "
          f"{r['auc'][0]:>7.3f}±{r['auc'][1]:.3f}")

print("-"*95)

import json
save = {k: {m: list(v[m]) for m in ['r2','rmse','rho','auc']} for k, v in all_results.items()}
with open('morgan_dmpnn_cv_results.json', 'w') as f:
    json.dump(save, f, indent=2)
print("\nSaved morgan_dmpnn_cv_results.json")