"""
Standardize the 560 raw decoys as a SEPARATE, explicit step.
Mirrors exactly what we did for antagonists.
Output: decoys_clean.csv with std_smiles + inchikey.
Also removes any decoy that accidentally collides with a TRPA1 compound by InChIKey.
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import SaltRemover

decoys = pd.read_csv('decoys_raw.csv')
print(f"Raw decoys loaded: {len(decoys)}")

remover = SaltRemover.SaltRemover()

def standardize(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        mol = remover.StripMol(mol, dontRemoveEverything=True)
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        can = Chem.MolToSmiles(mol, canonical=True)
        ik  = Chem.MolToInchiKey(mol)
        if not can or not ik:
            return None, None
        return can, ik
    except Exception:
        return None, None

results = decoys['canonical_smiles'].apply(standardize)
decoys['std_smiles'] = [r[0] for r in results]
decoys['inchikey']   = [r[1] for r in results]

failed = decoys['std_smiles'].isna().sum()
print(f"Failed to standardize (unparseable): {failed}")
decoys = decoys[decoys['std_smiles'].notna()].copy()

# Cross-check: drop any decoy whose InChIKey matches a TRPA1 compound
trpa1 = pd.read_csv('trpa1_human_clean.csv')
trpa1_keys = set(trpa1['inchikey'].dropna())
before = len(decoys)
decoys = decoys[~decoys['inchikey'].isin(trpa1_keys)].copy()
collisions = before - len(decoys)
print(f"Decoys colliding with TRPA1 set (removed): {collisions}")

# Also drop internal duplicate decoys (same InChIKey twice)
before = len(decoys)
decoys = decoys.drop_duplicates(subset='inchikey').copy()
print(f"Internal duplicate decoys removed: {before - len(decoys)}")

decoys.to_csv('decoys_clean.csv', index=False)
print(f"\nFinal clean decoys: {len(decoys)}")
print(f"Saved: decoys_clean.csv")