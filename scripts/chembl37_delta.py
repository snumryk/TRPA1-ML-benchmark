"""
Step 2: pull TRPA1 IC50 from current API (ChEMBL 37) using the EXACT same
filter + standardization as the original v36 pull. Compare by InChIKey to
the saved v36 file to identify NEW compounds = temporal holdout set.
"""
import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import SaltRemover

TRPA1_HUMAN = "CHEMBL6007"
remover = SaltRemover.SaltRemover()

# ---- 1. Pull activities (SAME filter as v36) ----
print("Pulling TRPA1 IC50 activities from ChEMBL 37...")
activity = new_client.activity
acts = activity.filter(
    target_chembl_id=TRPA1_HUMAN,
    standard_type="IC50",
    standard_relation="=",
    pchembl_value__isnull=False,
).only([
    'molecule_chembl_id', 'canonical_smiles', 'pchembl_value',
    'standard_value', 'standard_units', 'assay_type',
    'assay_chembl_id', 'document_chembl_id', 'document_year',
])
acts = list(acts)
print(f"Raw activities pulled: {len(acts)}")

df_raw = pd.DataFrame(acts)
print(f"Unique compounds (raw): {df_raw['molecule_chembl_id'].nunique()}")
print(f"Columns available: {list(df_raw.columns)}")

# ---- 2. Standardize (SAME pipeline as v36) ----
def standardize(smiles):
    if smiles is None:
        return None, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mol = remover.StripMol(mol, dontRemoveEverything=True)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    can = Chem.MolToSmiles(mol, canonical=True)
    ik = Chem.MolToInchiKey(mol)
    return can, ik

df_raw['pchembl_value'] = pd.to_numeric(df_raw['pchembl_value'], errors='coerce')
df_raw = df_raw.dropna(subset=['pchembl_value', 'canonical_smiles'])

std_smiles, inchikeys = [], []
for smi in df_raw['canonical_smiles']:
    cs, ik = standardize(smi)
    std_smiles.append(cs)
    inchikeys.append(ik)
df_raw['std_smiles'] = std_smiles
df_raw['inchikey'] = inchikeys
df_raw = df_raw.dropna(subset=['inchikey'])

# ---- 3. Deduplicate by (InChIKey) with median pchembl ----
df_v37 = df_raw.groupby('inchikey').agg(
    std_smiles=('std_smiles', 'first'),
    molecule_chembl_id=('molecule_chembl_id', 'first'),
    pchembl_median=('pchembl_value', 'median'),
    n_measurements=('pchembl_value', 'count'),
    year_min=('document_year', 'min'),
    year_max=('document_year', 'max'),
    n_documents=('document_chembl_id', 'nunique'),
    n_assays=('assay_chembl_id', 'nunique'),
).reset_index()
print(f"\nChEMBL 37 unique compounds after standardization: {len(df_v37)}")

# ---- 4. Compare to saved v36 file ----
df_v36 = pd.read_csv('trpa1_antagonists.csv')
v36_keys = set(df_v36['inchikey'])
v37_keys = set(df_v37['inchikey'])

new_keys = v37_keys - v36_keys
overlap = v37_keys & v36_keys
lost = v36_keys - v37_keys  # in v36 but not v37 (shouldn't happen much)

print(f"\n{'='*60}")
print("TEMPORAL DELTA: ChEMBL 36 -> ChEMBL 37")
print(f"{'='*60}")
print(f"  v36 compounds (saved file):     {len(v36_keys)}")
print(f"  v37 compounds (fresh pull):     {len(v37_keys)}")
print(f"  Overlap (in both):              {len(overlap)}")
print(f"  NEW in v37 (temporal holdout):  {len(new_keys)}")
print(f"  In v36 but not v37:             {len(lost)}")

# ---- 5. Characterize the new compounds ----
df_new = df_v37[df_v37['inchikey'].isin(new_keys)].copy()
if len(df_new) > 0:
    print(f"\nNEW compound set characterization:")
    print(f"  pchembl_median range: {df_new['pchembl_median'].min():.2f} - {df_new['pchembl_median'].max():.2f}")
    print(f"  pchembl_median mean:  {df_new['pchembl_median'].mean():.2f}")
    print(f"  year_min range: {df_new['year_min'].min()} - {df_new['year_min'].max()}")
    print(f"  compounds with >=2 measurements: {(df_new['n_measurements']>=2).sum()}")
    df_new.to_csv('trpa1_v37_new.csv', index=False)
    print(f"\n  Saved {len(df_new)} new compounds to trpa1_v37_new.csv")

    # verdict
    if len(new_keys) >= 100:
        print(f"\n  VERDICT: {len(new_keys)} new compounds — SUFFICIENT for temporal validation.")
    elif len(new_keys) >= 40:
        print(f"\n  VERDICT: {len(new_keys)} new compounds — MARGINAL but usable.")
    else:
        print(f"\n  VERDICT: {len(new_keys)} new compounds — TOO FEW; use document_year temporal split instead.")
else:
    print("\n  No new compounds found — API may still serve v36, or no new TRPA1 data.")