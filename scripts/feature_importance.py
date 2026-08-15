"""
Feature importance for TRPA1 potency — biological/SAR insight.
Three complementary views:
  1. RF impurity-based importance (fast, but biased to continuous features)
  2. Permutation importance with CV folds (honest, with error bars)
  3. Signed Spearman correlation of each descriptor with pIC50 (direction of effect)
Model: RF on 15 RDKit descriptors (interpretable, not the most accurate).
"""
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GroupKFold
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

SEED = 42

RDKIT_DESCS = [
    'MolWt', 'MolLogP', 'MolMR', 'TPSA',
    'NumHAcceptors', 'NumHDonors', 'NumRotatableBonds',
    'NumAromaticRings', 'RingCount', 'FractionCSP3',
    'HeavyAtomCount', 'NumAliphaticRings', 'NumSaturatedRings',
    'NumHeteroatoms', 'LabuteASA',
]

# Human-readable names for the paper
DESC_LABELS = {
    'MolWt': 'Molecular weight',
    'MolLogP': 'Lipophilicity (LogP)',
    'MolMR': 'Molar refractivity',
    'TPSA': 'Polar surface area',
    'NumHAcceptors': 'H-bond acceptors',
    'NumHDonors': 'H-bond donors',
    'NumRotatableBonds': 'Rotatable bonds',
    'NumAromaticRings': 'Aromatic rings',
    'RingCount': 'Total rings',
    'FractionCSP3': 'Fraction sp3 C',
    'HeavyAtomCount': 'Heavy atoms',
    'NumAliphaticRings': 'Aliphatic rings',
    'NumSaturatedRings': 'Saturated rings',
    'NumHeteroatoms': 'Heteroatoms',
    'LabuteASA': 'Accessible surface area',
}

# ── Load ──────────────────────────────────────────────────────
df = pd.read_csv('trpa1_antagonists.csv')

def compute_rdkit(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return np.full(len(RDKIT_DESCS), np.nan)
    return np.array([float(getattr(Descriptors, n)(mol)) for n in RDKIT_DESCS], dtype=np.float32)

X = np.vstack(df['std_smiles'].apply(compute_rdkit).values)
y = df['pchembl_median'].values
scaffolds = df['scaffold'].values

print(f"Data: {X.shape}, {len(np.unique(scaffolds))} scaffolds")

# ── 1. RF impurity importance (averaged over CV folds) ─────────
gkf = GroupKFold(n_splits=5)
impurity_imp = []
perm_imp = []

for tr, te in gkf.split(X, y, groups=scaffolds):
    rf = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=SEED)
    rf.fit(X[tr], y[tr])
    impurity_imp.append(rf.feature_importances_)

    # Permutation importance on the held-out fold
    perm = permutation_importance(rf, X[te], y[te], n_repeats=10,
                                  random_state=SEED, n_jobs=-1)
    perm_imp.append(perm.importances_mean)

impurity_imp = np.array(impurity_imp)
perm_imp = np.array(perm_imp)

imp_mean = impurity_imp.mean(axis=0)
imp_std  = impurity_imp.std(axis=0)
perm_mean = perm_imp.mean(axis=0)
perm_std  = perm_imp.std(axis=0)

# ── 2. Signed Spearman correlation (direction of effect) ───────
# Compute on full data (direction is a descriptive statistic)
directions = []
for i in range(len(RDKIT_DESCS)):
    valid = ~np.isnan(X[:, i])
    rho = spearmanr(X[valid, i], y[valid]).correlation
    directions.append(rho)
directions = np.array(directions)

# ── 3. Assemble and rank ──────────────────────────────────────
result = pd.DataFrame({
    'descriptor': [DESC_LABELS[d] for d in RDKIT_DESCS],
    'impurity_imp': imp_mean,
    'impurity_std': imp_std,
    'perm_imp': perm_mean,
    'perm_std': perm_std,
    'spearman_with_pIC50': directions,
})
result = result.sort_values('perm_imp', ascending=False).reset_index(drop=True)

pd.set_option('display.width', 140)
pd.set_option('display.max_columns', None)

print("\n" + "="*100)
print("FEATURE IMPORTANCE FOR TRPA1 ANTAGONIST POTENCY (RF on 15 RDKit descriptors, 5-fold CV)")
print("="*100)
print(f"{'Descriptor':<26} {'Permut.Imp':>16} {'Impurity Imp':>16} {'Spearman r':>12} {'Effect'}")
print("-"*100)
for _, row in result.iterrows():
    direction = "↑ potency" if row['spearman_with_pIC50'] > 0.05 else \
                "↓ potency" if row['spearman_with_pIC50'] < -0.05 else "~ neutral"
    print(f"{row['descriptor']:<26} "
          f"{row['perm_imp']:>7.4f}±{row['perm_std']:.4f} "
          f"{row['impurity_imp']:>7.4f}±{row['impurity_std']:.4f} "
          f"{row['spearman_with_pIC50']:>12.3f}  {direction}")

print("-"*100)

# ── 4. SAR narrative for the paper ────────────────────────────
top5 = result.head(5)
print("\nTOP-5 DRIVERS OF TRPA1 ANTAGONIST POTENCY (by permutation importance):")
for i, (_, row) in enumerate(top5.iterrows(), 1):
    direction = "higher" if row['spearman_with_pIC50'] > 0.05 else \
                "lower" if row['spearman_with_pIC50'] < -0.05 else "context-dependent"
    print(f"  {i}. {row['descriptor']}: {direction} values → more potent "
          f"(Spearman r={row['spearman_with_pIC50']:.3f})")

# Save for the paper
result.to_csv('feature_importance_results.csv', index=False)
print(f"\nSaved feature_importance_results.csv")