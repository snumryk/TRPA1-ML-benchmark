"""
Standardize decoys: salt strip, canonical SMILES, InChIKey.
Remove any decoy that collides with TRPA1 set.
Remove internal duplicates.
Output: decoys_clean.csv
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import SaltRemover, Descriptors

# ── 1. Load ────────────────────────────────────────────────────
decoys = pd.read_csv('decoys_raw.csv')
print(f"Raw decoys loaded: {len(decoys)}")

# ── 2. Standardize ─────────────────────────────────────────────
remover = SaltRemover.SaltRemover()

def standardize(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None, None
        mol = remover.StripMol(mol, dontRemoveEverything=True)
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        can = Chem.MolToSmiles(mol, canonical=True)
        ik = Chem.MolToInchiKey(mol)
        mw = Descriptors.MolWt(mol)
        if not can or not ik:
            return None, None, None
        return can, ik, mw
    except Exception:
        return None, None, None

print("Standardizing...")
results = [standardize(s) for s in decoys['canonical_smiles']]
decoys['std_smiles'] = [r[0] for r in results]
decoys['inchikey'] = [r[1] for r in results]
decoys['mw'] = [r[2] for r in results]

failed = decoys['std_smiles'].isna().sum()
print(f"Failed to parse: {failed}")
decoys = decoys[decoys['std_smiles'].notna()].copy()

# ── 3. Remove collisions with TRPA1 set ───────────────────────
trpa1 = pd.read_csv('trpa1_antagonists.csv')
trpa1_keys = set(trpa1['inchikey'].dropna())
before = len(decoys)
decoys = decoys[~decoys['inchikey'].isin(trpa1_keys)].copy()
collisions = before - len(decoys)
print(f"Collisions with TRPA1 antagonists removed: {collisions}")

# Also check against full dataset (includes agonists)
trpa1_full = pd.read_csv('trpa1_human_clean.csv')
trpa1_all_keys = set(trpa1_full['inchikey'].dropna())
before = len(decoys)
decoys = decoys[~decoys['inchikey'].isin(trpa1_all_keys)].copy()
extra_collisions = before - len(decoys)
print(f"Collisions with agonists removed: {extra_collisions}")

# ── 4. Remove internal duplicates ──────────────────────────────
before = len(decoys)
decoys = decoys.drop_duplicates(subset='inchikey').copy()
print(f"Internal duplicates removed: {before - len(decoys)}")

# ── 5. Summary ─────────────────────────────────────────────────
decoys.to_csv('decoys_clean.csv', index=False)
print(f"\nFinal clean decoys: {len(decoys)}")
print(f"MW range: {decoys['mw'].min():.0f} - {decoys['mw'].max():.0f}, "
      f"mean={decoys['mw'].mean():.0f}")
print(f"Saved: decoys_clean.csv")
print(f"\nTarget ratio check:")
print(f"  Antagonists: 1645")
print(f"  Decoys: {len(decoys)}")
print(f"  Ratio: {1645/len(decoys):.2f}:1 (Mihai had 2.92:1)")