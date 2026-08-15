"""
ChemBERTa-2 pretrained embeddings for TRPA1 antagonists.
Two pooling strategies: CLS token vs Mean Pooling (Gemini's suggestion).
Combined with RDKit descriptors and classical baselines.
"""
import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from sklearn.metrics import r2_score, mean_squared_error, roc_auc_score, matthews_corrcoef
from scipy.stats import spearmanr
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

SEED = 42
THRESHOLD = 7.0

# ── 1. Load ChemBERTa ─────────────────────────────────────────
MODEL_NAME = "DeepChem/ChemBERTa-77M-MTR"
print(f"Loading {MODEL_NAME}...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
chemberta = AutoModel.from_pretrained(MODEL_NAME)
chemberta.eval()
print("ChemBERTa loaded.")

# ── 2. Load dataset ───────────────────────────────────────────
df = pd.read_csv('trpa1_antagonists.csv')
print(f"Dataset: {len(df)} compounds")

# ── 3. Extract BOTH CLS and Mean Pooling embeddings ───────────
print("Extracting embeddings (CLS + Mean Pooling)...")

cls_embeddings = []
mean_embeddings = []

for smi in tqdm(df['std_smiles'], desc="Embedding"):
    tokens = tokenizer(smi, return_tensors="pt", padding=True,
                       truncation=True, max_length=512)
    with torch.no_grad():
        output = chemberta(**tokens)
    hidden = output.last_hidden_state  # shape: (1, seq_len, 384)

    # CLS: first token
    cls_emb = hidden[:, 0, :].numpy().ravel()

    # Mean pooling: average over all tokens (excluding padding)
    mask = tokens['attention_mask'].unsqueeze(-1).float()  # (1, seq_len, 1)
    mean_emb = (hidden * mask).sum(dim=1) / mask.sum(dim=1)  # (1, 384)
    mean_emb = mean_emb.numpy().ravel()

    cls_embeddings.append(cls_emb)
    mean_embeddings.append(mean_emb)

print(f"CLS embedding shape:  {cls_embeddings[0].shape}")
print(f"Mean embedding shape: {mean_embeddings[0].shape}")

# ── 4. RDKit descriptors ──────────────────────────────────────
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
        return np.zeros(len(RDKIT_DESCS))
    feats = []
    for name in RDKIT_DESCS:
        try:
            feats.append(float(getattr(Descriptors, name)(mol)))
        except Exception:
            feats.append(0.0)
    return np.array(feats, dtype=np.float32)

df['rdkit'] = df['std_smiles'].apply(compute_rdkit)
df['cls_emb'] = cls_embeddings
df['mean_emb'] = mean_embeddings

# ── 5. Prepare splits ─────────────────────────────────────────
train = df[df['split'] == 'train']
test  = df[df['split'] == 'test']

y_train = train['pchembl_median'].values
y_test  = test['pchembl_median'].values

def stack(series):
    return np.vstack(series.values)

# All feature combinations
features = {
    'CLS(384)':           (stack(train['cls_emb']),  stack(test['cls_emb'])),
    'MeanPool(384)':      (stack(train['mean_emb']), stack(test['mean_emb'])),
    'RDKit(15)':          (stack(train['rdkit']),     stack(test['rdkit'])),
    'CLS+RDKit':          (np.hstack([stack(train['cls_emb']),  stack(train['rdkit'])]),
                           np.hstack([stack(test['cls_emb']),   stack(test['rdkit'])])),
    'MeanPool+RDKit':     (np.hstack([stack(train['mean_emb']), stack(train['rdkit'])]),
                           np.hstack([stack(test['mean_emb']),  stack(test['rdkit'])])),
}

print(f"\nTrain: {len(y_train)}, Test: {len(y_test)}")

# ── 6. Run all combinations ───────────────────────────────────
def evaluate(name, y_true, y_pred):
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    rho = spearmanr(y_true, y_pred).correlation
    y_cls = (y_true >= THRESHOLD).astype(int)
    auc = roc_auc_score(y_cls, y_pred)
    mcc = matthews_corrcoef(y_cls, (y_pred >= THRESHOLD).astype(int))
    return rmse, r2, rho, auc, mcc

print(f"\n{'='*85}")
print("ChemBERTa-2 EMBEDDING BENCHMARK (CLS vs Mean Pooling)")
print(f"{'='*85}")
print(f"{'Model':<35} {'RMSE':>6} {'R2':>6} {'Spearman':>9} {'AUC':>6} {'MCC':>6}")
print("-"*85)

all_results = []

for feat_name, (X_tr, X_te) in features.items():
    # XGBoost for all
    xgb = XGBRegressor(n_estimators=500, max_depth=6, learning_rate=0.05,
                       n_jobs=-1, random_state=SEED)
    xgb.fit(X_tr, y_train)
    preds = xgb.predict(X_te)
    rmse, r2, rho, auc, mcc = evaluate(feat_name, y_test, preds)
    label = f"XGB + {feat_name}"
    print(f"  {label:<35} {rmse:>6.3f} {r2:>6.3f} {rho:>9.3f} {auc:>6.3f} {mcc:>6.3f}")
    all_results.append((label, rmse, r2, rho, auc, mcc))

# Also RF on best embedding combo
for feat_name in ['MeanPool+RDKit', 'CLS+RDKit']:
    X_tr, X_te = features[feat_name]
    rf = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=SEED)
    rf.fit(X_tr, y_train)
    preds = rf.predict(X_te)
    rmse, r2, rho, auc, mcc = evaluate(feat_name, y_test, preds)
    label = f"RF + {feat_name}"
    print(f"  {label:<35} {rmse:>6.3f} {r2:>6.3f} {rho:>9.3f} {auc:>6.3f} {mcc:>6.3f}")
    all_results.append((label, rmse, r2, rho, auc, mcc))

# ── 7. Full comparison ────────────────────────────────────────
print(f"\n{'='*85}")
print("FULL COMPARISON: all experiments to date")
print(f"{'='*85}")
print(f"{'Model':<35} {'RMSE':>6} {'R2':>6} {'Spearman':>9} {'AUC':>6}")
print("-"*85)

prev = [
    ("RF (Morgan 2048)",           0.795, 0.232, 0.552, 0.795),
    ("XGB (Morgan 2048)",          0.781, 0.257, 0.547, 0.790),
    ("D-MPNN (graph only)",        0.823, 0.176, 0.554, 0.761),
    ("D-MPNN+RDKit (ES)",          0.763, 0.292, 0.621, 0.810),
    ("RF (15 RDKit) *prev best*",  0.724, 0.363, 0.624, 0.800),
]
for name, rmse, r2, rho, auc in prev:
    print(f"  {name:<35} {rmse:>6.3f} {r2:>6.3f} {rho:>9.3f} {auc:>6.3f}")

print()
for name, rmse, r2, rho, auc, mcc in all_results:
    marker = " <-- NEW BEST" if r2 > 0.363 else ""
    print(f"  {name:<35} {rmse:>6.3f} {r2:>6.3f} {rho:>9.3f} {auc:>6.3f}{marker}")

print("-"*85)
best = max(all_results, key=lambda x: x[2])
print(f"\nBest new: {best[0]} (R2={best[2]:.3f})")
print(f"vs RF(15 RDKit): R2 delta = {best[2]-0.363:+.3f}")
print(f"vs XGB(Morgan):  R2 delta = {best[2]-0.257:+.3f}")