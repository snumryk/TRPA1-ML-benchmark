from chembl_webresource_client.new_client import new_client
import pandas as pd
import random

random.seed(34)  # Mihai used random_state=34

# ── 1. Load our existing TRPA1 compounds to exclude them ──────
trpa1 = pd.read_csv('trpa1_human_clean.csv')
trpa1_chembl_ids = set(trpa1['molecule_chembl_id'].dropna())
print(f"TRPA1 compounds to exclude: {len(trpa1_chembl_ids)}")

# ── 2. Get the total count of molecules in ChEMBL ─────────────
molecule = new_client.molecule
# We pull molecules that have a canonical SMILES and are small molecules
print("Fetching a pool of random ChEMBL molecules...")

# Strategy: molecules have IDs like CHEMBL1, CHEMBL2, ... up to ~CHEMBL5M+
# We sample random integer IDs and fetch those molecules
N_DECOYS_TARGET = 560
pool = []
attempts = 0
max_attempts = 3000

# ChEMBL 36 has ~2.4M compounds; IDs are sparse, so we oversample
mol_client = new_client.molecule

checked = 0
while len(pool) < N_DECOYS_TARGET and attempts < max_attempts:
    # Random ChEMBL ID in a plausible range
    rand_id = f"CHEMBL{random.randint(1, 4000000)}"
    attempts += 1
    if rand_id in trpa1_chembl_ids:
        continue
    try:
        mol = mol_client.get(rand_id)
    except Exception:
        continue
    if mol is None:
        continue
    # Need a molecule with a canonical SMILES
    structures = mol.get('molecule_structures')
    if not structures:
        continue
    smi = structures.get('canonical_smiles')
    if not smi:
        continue
    pool.append({
        'molecule_chembl_id': rand_id,
        'canonical_smiles': smi,
        'pref_name': mol.get('pref_name'),
    })
    checked += 1
    if len(pool) % 50 == 0:
        print(f"  collected {len(pool)} decoys (after {attempts} attempts)...")

print(f"\nDone. Collected {len(pool)} decoys after {attempts} attempts.")

df_decoys = pd.DataFrame(pool)
df_decoys.to_csv('decoys_raw.csv', index=False)
print(f"Saved: decoys_raw.csv")
print(f"\nSample decoys:")
print(df_decoys.head(10).to_string())