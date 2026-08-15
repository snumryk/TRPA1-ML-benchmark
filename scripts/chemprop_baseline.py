"""
D-MPNN+RDKit with Early Stopping and proper LR warmup.
Addresses reviewer criticism of unfair comparison.
"""
import numpy as np
import pandas as pd
import torch
import lightning as L
from lightning.pytorch.callbacks import EarlyStopping
from rdkit import Chem
from rdkit.Chem import Descriptors
from chemprop.data import MoleculeDatapoint, MoleculeDataset, build_dataloader
from chemprop.nn import BondMessagePassing, MeanAggregation, RegressionFFN
from chemprop.models import MPNN
from sklearn.metrics import (r2_score, mean_squared_error, mean_absolute_error,
                             roc_auc_score, matthews_corrcoef, balanced_accuracy_score)
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

SEED = 42
THRESHOLD = 7.0
L.seed_everything(SEED)

# ── 1. RDKit descriptors ──────────────────────────────────────
RDKIT_DESCS = [
    'MolWt', 'MolLogP', 'MolMR', 'TPSA',
    'NumHAcceptors', 'NumHDonors', 'NumRotatableBonds',
    'NumAromaticRings', 'RingCount', 'FractionCSP3',
    'HeavyAtomCount', 'NumAliphaticRings', 'NumSaturatedRings',
    'NumHeteroatoms', 'LabuteASA',
]

def compute_rdkit_features(mol):
    feats = []
    for name in RDKIT_DESCS:
        func = getattr(Descriptors, name, None)
        if func is None:
            feats.append(0.0)
        else:
            try:
                feats.append(float(func(mol)))
            except Exception:
                feats.append(0.0)
    return np.array(feats, dtype=np.float32)

# ── 2. Load data ──────────────────────────────────────────────
df = pd.read_csv('trpa1_antagonists.csv')

def make_datapoints(subset):
    points = []
    for _, row in subset.iterrows():
        mol = Chem.MolFromSmiles(row['std_smiles'])
        if mol is None:
            continue
        x_d = compute_rdkit_features(mol)
        points.append(
            MoleculeDatapoint(mol=mol, y=np.array([row['pchembl_median']]), x_d=x_d)
        )
    return points

train_points = make_datapoints(df[df['split'] == 'train'])
val_points   = make_datapoints(df[df['split'] == 'val'])
test_points  = make_datapoints(df[df['split'] == 'test'])

train_dset = MoleculeDataset(train_points)
val_dset   = MoleculeDataset(val_points)
test_dset  = MoleculeDataset(test_points)

train_loader = build_dataloader(train_dset, batch_size=64, shuffle=True, seed=SEED)
val_loader   = build_dataloader(val_dset, batch_size=64, shuffle=False)
test_loader  = build_dataloader(test_dset, batch_size=64, shuffle=False)

print(f"Train: {len(train_dset)}, Val: {len(val_dset)}, Test: {len(test_dset)}")

# ── 3. Build model with Chemprop's built-in LR schedule ───────
mp = BondMessagePassing(d_h=300, depth=3)
agg = MeanAggregation()
ffn = RegressionFFN(input_dim=300 + len(RDKIT_DESCS))

model = MPNN(
    message_passing=mp,
    agg=agg,
    predictor=ffn,
    warmup_epochs=2,       # linear warmup from init_lr to max_lr
    init_lr=0.0001,        # starting LR
    max_lr=0.001,          # peak LR after warmup
    final_lr=0.0001,       # LR decays back to this
)

print(f"Parameters: {sum(p.numel() for p in model.parameters()):,}")

# ── 4. Train with Early Stopping ──────────────────────────────
early_stop = EarlyStopping(
    monitor='val_loss',
    patience=10,           # stop if val_loss doesn't improve for 10 epochs
    mode='min',
    verbose=True,
)

trainer = L.Trainer(
    max_epochs=100,        # upper bound, early stopping will cut shorter
    accelerator='cpu',
    enable_progress_bar=True,
    enable_model_summary=False,
    logger=False,
    callbacks=[early_stop],
)

print("\nTraining D-MPNN+RDKit (max 100 epochs, early stopping patience=10)...")
trainer.fit(model, train_loader, val_loader)

stopped_epoch = trainer.current_epoch + 1
print(f"\nTraining stopped at epoch {stopped_epoch}")

# ── 5. Predict + metrics ──────────────────────────────────────
results = trainer.predict(model, test_loader)
preds = torch.cat(results).numpy().ravel()
actuals = np.array([p.y[0] for p in test_points])

# Regression
rmse = np.sqrt(mean_squared_error(actuals, preds))
mae  = mean_absolute_error(actuals, preds)
r2   = r2_score(actuals, preds)
rho  = spearmanr(actuals, preds).correlation

rng = np.random.default_rng(SEED)
boot_r2 = []
for _ in range(1000):
    idx = rng.integers(0, len(actuals), len(actuals))
    boot_r2.append(r2_score(actuals[idx], preds[idx]))
r2_ci = np.percentile(boot_r2, [2.5, 97.5])

# Classification
y_cls = (actuals >= THRESHOLD).astype(int)
p_cls = (preds >= THRESHOLD).astype(int)
auc = roc_auc_score(y_cls, preds)
mcc = matthews_corrcoef(y_cls, p_cls)
bacc = balanced_accuracy_score(y_cls, p_cls)

boot_auc = []
for _ in range(1000):
    idx = rng.integers(0, len(y_cls), len(y_cls))
    if len(np.unique(y_cls[idx])) < 2:
        continue
    boot_auc.append(roc_auc_score(y_cls[idx], preds[idx]))
auc_ci = np.percentile(boot_auc, [2.5, 97.5])

# ── 6. Full comparison table ──────────────────────────────────
print(f"\n{'='*75}")
print("FINAL HEAD-TO-HEAD: All models on scaffold split")
print(f"{'='*75}")
print(f"{'Model':<20} {'RMSE':>7} {'R2':>7} {'Spearman':>10} {'AUC':>7} {'MCC':>7}")
print("-"*75)
print(f"{'RF':<20} {'0.795':>7} {'0.232':>7} {'0.552':>10} {'0.795':>7} {'0.360':>7}")
print(f"{'XGBoost':<20} {'0.781':>7} {'0.257':>7} {'0.547':>10} {'0.790':>7} {'0.378':>7}")
print(f"{'LightGBM':<20} {'0.782':>7} {'0.255':>7} {'0.547':>10} {'0.772':>7} {'0.352':>7}")
print(f"{'D-MPNN':<20} {'0.823':>7} {'0.176':>7} {'0.554':>10} {'0.761':>7} {'0.324':>7}")
print(f"{'D-MPNN+RDKit(50ep)':<20} {'0.784':>7} {'0.253':>7} {'0.619':>10} {'0.804':>7} {'0.353':>7}")
print(f"{'D-MPNN+RDKit(ES)':<20} {rmse:>7.3f} {r2:>7.3f} {rho:>10.3f} {auc:>7.3f} {mcc:>7.3f}")
print("-"*75)
print(f"D-MPNN+RDKit(ES) stopped at epoch: {stopped_epoch}")
print(f"D-MPNN+RDKit(ES) R2 95% CI:  [{r2_ci[0]:.3f}, {r2_ci[1]:.3f}]")
print(f"D-MPNN+RDKit(ES) AUC 95% CI: [{auc_ci[0]:.3f}, {auc_ci[1]:.3f}]")
print(f"\nvs best classical (XGBoost):")
print(f"  R2 delta:       {r2 - 0.257:+.3f}")
print(f"  AUC delta:      {auc - 0.795:+.3f}")
print(f"  Spearman delta: {rho - 0.552:+.3f}")