"""
Replication of Mihai et al. (2020) on ChEMBL 36 data.
Task: TRPA1 antagonist (class 1) vs random ChEMBL decoy (class 0)
Features: Morgan ECFP4 2048 bits (proxy for their MNA level-3 1376)
Split: random 80/20 stratified, random_state=34
Models: RF, SVM, FFNN with their exact hyperparameters
Metrics: TPR, TNR, ACC, bACC, FPR, NPV, ROC AUC + 10-fold CV
"""
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import confusion_matrix, roc_auc_score
import warnings
warnings.filterwarnings('ignore')

SEED = 34  # Mihai's random_state

# ── 1. Load data ───────────────────────────────────────────────
ant = pd.read_csv('trpa1_antagonists.csv')
dec = pd.read_csv('decoys_clean.csv')
print(f"Antagonists (class 1): {len(ant)}")
print(f"Decoys (class 0):      {len(dec)}")
print(f"Ratio: {len(ant)/len(dec):.2f}:1")

# ── 2. Morgan fingerprints ─────────────────────────────────────
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def smi_to_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = mfpgen.GetFingerprint(mol)
    arr = np.zeros((2048,), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

print("\nComputing fingerprints...")
ant_fps = [smi_to_fp(s) for s in ant['std_smiles']]
dec_fps = [smi_to_fp(s) for s in dec['std_smiles']]

# Remove any None (failed parse)
ant_valid = [(fp, 1) for fp in ant_fps if fp is not None]
dec_valid = [(fp, 0) for fp in dec_fps if fp is not None]
print(f"Valid fingerprints: {len(ant_valid)} antagonists, {len(dec_valid)} decoys")

all_data = ant_valid + dec_valid
X = np.vstack([d[0] for d in all_data])
y = np.array([d[1] for d in all_data])
print(f"Combined: X={X.shape}, class 1={y.sum()}, class 0={len(y)-y.sum()}")

# ── 3. Random 80/20 stratified split ──────────────────────────
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.20, random_state=SEED, stratify=y)
print(f"\nTrain: {len(y_train)} (class 1: {y_train.sum()}, class 0: {len(y_train)-y_train.sum()})")
print(f"Test:  {len(y_test)} (class 1: {y_test.sum()}, class 0: {len(y_test)-y_test.sum()})")

# ── 4. Metrics function (exact same as Mihai Table 1) ─────────
def compute_metrics(y_true, y_pred, y_prob):
    cm = confusion_matrix(y_true, y_pred)
    tn, fp, fn, tp = cm.ravel()
    tpr = tp / (tp + fn) if (tp + fn) else 0
    tnr = tn / (tn + fp) if (tn + fp) else 0
    acc = (tp + tn) / (tp + tn + fp + fn)
    bacc = (tpr + tnr) / 2
    fpr = fp / (tn + fp) if (tn + fp) else 0
    npv = tn / (tn + fn) if (tn + fn) else 0
    auc = roc_auc_score(y_true, y_prob)
    return {
        'TPR': tpr, 'TNR': tnr, 'ACC': acc, 'bACC': bacc,
        'FPR': fpr, 'NPV': npv, 'AUC': auc
    }

def print_metrics(name, m):
    print(f"\n  {name}:")
    print(f"    TPR (sensitivity): {m['TPR']*100:6.2f}%")
    print(f"    TNR (specificity): {m['TNR']*100:6.2f}%")
    print(f"    ACC:               {m['ACC']*100:6.2f}%")
    print(f"    bACC:              {m['bACC']*100:6.2f}%")
    print(f"    FPR:               {m['FPR']*100:6.2f}%")
    print(f"    NPV:               {m['NPV']*100:6.2f}%")
    print(f"    ROC AUC:           {m['AUC']:.4f}")

# ── 5. 10-fold CV function ────────────────────────────────────
def cv_auc(model_fn, X, y):
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
    aucs = []
    for tr_idx, te_idx in cv.split(X, y):
        model = model_fn()
        model.fit(X[tr_idx], y[tr_idx])
        prob = model.predict_proba(X[te_idx])[:, 1]
        aucs.append(roc_auc_score(y[te_idx], prob))
    return np.mean(aucs), np.std(aucs)

# ══════════════════════════════════════════════════════════════
# MODEL 1: RANDOM FOREST (Mihai's exact params)
# ══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("MODEL 1: RANDOM FOREST")
print("  n_estimators=50, max_depth=90, max_features=sqrt")
print("  min_samples_split=2, min_samples_leaf=1, random_state=34")
print("="*60)

def make_rf():
    return RandomForestClassifier(
        n_estimators=50, max_depth=90, max_features='sqrt',
        min_samples_split=2, min_samples_leaf=1, random_state=SEED)

rf = make_rf()
rf.fit(X_train, y_train)
rf_pred = rf.predict(X_test)
rf_prob = rf.predict_proba(X_test)[:, 1]
rf_m = compute_metrics(y_test, rf_pred, rf_prob)
print_metrics("Test set", rf_m)

rf_cv_mean, rf_cv_std = cv_auc(make_rf, X, y)
print(f"\n  10-fold CV Mean AUC: {rf_cv_mean:.4f} (+/-{rf_cv_std:.4f})")

# ══════════════════════════════════════════════════════════════
# MODEL 2: SVM (Mihai's exact params)
# ══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("MODEL 2: SVM")
print("  C=8, gamma=0.001, kernel=rbf")
print("="*60)

def make_svm():
    return SVC(C=8, gamma=0.001, kernel='rbf',
               probability=True, random_state=SEED)

svm = make_svm()
svm.fit(X_train, y_train)
svm_pred = svm.predict(X_test)
svm_prob = svm.predict_proba(X_test)[:, 1]
svm_m = compute_metrics(y_test, svm_pred, svm_prob)
print_metrics("Test set", svm_m)

svm_cv_mean, svm_cv_std = cv_auc(make_svm, X, y)
print(f"\n  10-fold CV Mean AUC: {svm_cv_mean:.4f} (+/-{svm_cv_std:.4f})")

# ══════════════════════════════════════════════════════════════
# MODEL 3: FFNN (Mihai's architecture)
# ══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("MODEL 3: FFNN")
print("  Architecture: 2048 -> 750 (ReLU, dropout=0.6) -> 1 (sigmoid)")
print("  Loss: binary_crossentropy, epochs=30, batch=16")
print("="*60)

# Try TensorFlow first, fall back to sklearn MLP
try:
    import tensorflow as tf
    tf.random.set_seed(SEED)
    from tensorflow.keras.models import Sequential
    from tensorflow.keras.layers import Dense, Dropout, Input
    USE_TF = True
    print("  Using TensorFlow/Keras")
except ImportError:
    USE_TF = False
    print("  TensorFlow not available, using sklearn MLPClassifier")

if USE_TF:
    def make_ffnn():
        m = Sequential([
            Input(shape=(2048,)),
            Dense(750, activation='relu'),
            Dropout(0.6),
            Dense(1, activation='sigmoid'),
        ])
        m.compile(optimizer='adam', loss='binary_crossentropy')
        return m

    ffnn = make_ffnn()
    ffnn.fit(X_train, y_train, epochs=30, batch_size=16, verbose=0)
    ffnn_prob = ffnn.predict(X_test, verbose=0).ravel()
    ffnn_pred = (ffnn_prob >= 0.5).astype(int)
    ffnn_m = compute_metrics(y_test, ffnn_pred, ffnn_prob)
    print_metrics("Test set", ffnn_m)

    # 10-fold CV for FFNN
    print("\n  10-fold CV (this takes a few minutes)...")
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
    ffnn_cv_aucs = []
    for fold, (tr_idx, te_idx) in enumerate(cv.split(X, y)):
        tf.random.set_seed(SEED)
        m = make_ffnn()
        m.fit(X[tr_idx], y[tr_idx], epochs=30, batch_size=16, verbose=0)
        prob = m.predict(X[te_idx], verbose=0).ravel()
        auc = roc_auc_score(y[te_idx], prob)
        ffnn_cv_aucs.append(auc)
        print(f"    fold {fold+1}/10: AUC={auc:.4f}")
    ffnn_cv_mean = np.mean(ffnn_cv_aucs)
    ffnn_cv_std = np.std(ffnn_cv_aucs)
    print(f"\n  10-fold CV Mean AUC: {ffnn_cv_mean:.4f} (+/-{ffnn_cv_std:.4f})")

else:
    from sklearn.neural_network import MLPClassifier
    def make_mlp():
        return MLPClassifier(
            hidden_layer_sizes=(750,), activation='relu',
            max_iter=30, batch_size=16, random_state=SEED)
    mlp = make_mlp()
    mlp.fit(X_train, y_train)
    mlp_prob = mlp.predict_proba(X_test)[:, 1]
    mlp_pred = mlp.predict(X_test)
    ffnn_m = compute_metrics(y_test, mlp_pred, mlp_prob)
    print_metrics("Test set", ffnn_m)
    ffnn_cv_mean, ffnn_cv_std = cv_auc(make_mlp, X, y)
    print(f"\n  10-fold CV Mean AUC: {ffnn_cv_mean:.4f} (+/-{ffnn_cv_std:.4f})")

# ══════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ══════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("COMPARISON: Mihai et al. (2020) vs Our Replication (ChEMBL 36)")
print("="*70)
print(f"\nMihai: 371 actives + 127 decoys, MNA level-3, random 80/20")
print(f"Ours:  1645 actives + 560 decoys, Morgan ECFP4, random 80/20")
print(f"\n{'':>8} {'':>8} {'--- Mihai ---':>24} {'--- Ours ---':>24}")
print(f"{'Model':<8} {'Metric':<8} {'Value':>12}         {'Value':>12}")
print("-"*70)
print(f"{'RF':<8} {'ACC':>8} {'99.00%':>12}         {rf_m['ACC']*100:>11.2f}%")
print(f"{'RF':<8} {'bACC':>8} {'98.00%':>12}         {rf_m['bACC']*100:>11.2f}%")
print(f"{'RF':<8} {'TPR':>8} {'100.00%':>12}         {rf_m['TPR']*100:>11.2f}%")
print(f"{'RF':<8} {'TNR':>8} {'96.00%':>12}         {rf_m['TNR']*100:>11.2f}%")
print(f"{'RF':<8} {'AUC CV':>8} {'0.9936':>12}         {rf_cv_mean:>11.4f}")
print(f"{'':>8} {'':>8}")
print(f"{'SVM':<8} {'ACC':>8} {'90.00%':>12}         {svm_m['ACC']*100:>11.2f}%")
print(f"{'SVM':<8} {'bACC':>8} {'88.00%':>12}         {svm_m['bACC']*100:>11.2f}%")
print(f"{'SVM':<8} {'TPR':>8} {'92.00%':>12}         {svm_m['TPR']*100:>11.2f}%")
print(f"{'SVM':<8} {'TNR':>8} {'84.00%':>12}         {svm_m['TNR']*100:>11.2f}%")
print(f"{'SVM':<8} {'AUC CV':>8} {'0.9354':>12}         {svm_cv_mean:>11.4f}")
print(f"{'':>8} {'':>8}")
print(f"{'FFNN':<8} {'ACC':>8} {'88.00%':>12}         {ffnn_m['ACC']*100:>11.2f}%")
print(f"{'FFNN':<8} {'bACC':>8} {'85.33%':>12}         {ffnn_m['bACC']*100:>11.2f}%")
print(f"{'FFNN':<8} {'TPR':>8} {'90.67%':>12}         {ffnn_m['TPR']*100:>11.2f}%")
print(f"{'FFNN':<8} {'TNR':>8} {'80.00%':>12}         {ffnn_m['TNR']*100:>11.2f}%")
print(f"{'FFNN':<8} {'AUC CV':>8} {'0.9354':>12}         {ffnn_cv_mean:>11.4f}")

print("\n" + "="*70)
print("INTERPRETATION")
print("="*70)
print("If our numbers are close to Mihai's -> their result is reproducible")
print("  and the high performance comes from the TASK FRAMING (ligand vs")
print("  random decoy), not from MNA descriptors specifically.")
print("If our numbers are much lower -> descriptor choice matters more")
print("  than task framing, and MNA may genuinely outperform Morgan here.")