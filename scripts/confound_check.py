"""
Confound check: do the top potency-driving descriptors correlate with publication year?
If TPSA/heteroatoms correlate strongly with year_min, part of the "SAR signal"
may reflect medicinal-chemistry trends over time, not TRPA1 biology.
"""
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

df = pd.read_csv('trpa1_antagonists.csv')
print(f"Compounds: {len(df)}")
print(f"Year range: {df['year_min'].min():.0f} - {df['year_min'].max():.0f}")
print(f"Compounds with year data: {df['year_min'].notna().sum()}")

# Top descriptors from feature importance
TOP_DESCS = ['TPSA', 'NumHeteroatoms', 'MolMR', 'NumHAcceptors', 'MolLogP', 'NumHDonors']

def compute_desc(smi, name):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return np.nan
    try:
        return float(getattr(Descriptors, name)(mol))
    except Exception:
        return np.nan

# Only rows with year data
sub = df[df['year_min'].notna()].copy()
year = sub['year_min'].values
pchembl = sub['pchembl_median'].values

print(f"\n{'='*85}")
print("CONFOUND CHECK: descriptor vs YEAR vs POTENCY")
print(f"{'='*85}")
print(f"{'Descriptor':<20} {'corr(desc,year)':>16} {'corr(desc,pIC50)':>18} {'Interpretation'}")
print("-"*85)

# First: does potency itself correlate with year?
rho_year_pot = spearmanr(year, pchembl).correlation
print(f"{'pIC50 itself':<20} {'—':>16} {rho_year_pot:>18.3f}  "
      f"{'potency rises over time' if rho_year_pot > 0.1 else 'no strong year trend'}")
print("-"*85)

for name in TOP_DESCS:
    vals = np.array([compute_desc(s, name) for s in sub['std_smiles']])
    valid = ~np.isnan(vals)
    rho_dy = spearmanr(vals[valid], year[valid]).correlation
    rho_dp = spearmanr(vals[valid], pchembl[valid]).correlation

    # Flag if descriptor-year correlation is a meaningful fraction of descriptor-potency
    if abs(rho_dy) > 0.3:
        flag = "⚠️ year confound possible"
    elif abs(rho_dy) > 0.15:
        flag = "~ mild year trend"
    else:
        flag = "✓ year-independent"

    print(f"{name:<20} {rho_dy:>16.3f} {rho_dp:>18.3f}  {flag}")

print("-"*85)

# ── Partial correlation: descriptor vs potency, CONTROLLING for year ──
# If partial corr stays strong after removing year, the signal is real.
print(f"\n{'='*85}")
print("PARTIAL CORRELATION: descriptor vs pIC50, controlling for year")
print("(If partial ≈ raw correlation → signal is NOT explained by year)")
print(f"{'='*85}")

def partial_corr(x, y, z):
    """Spearman partial correlation of x,y controlling for z."""
    from scipy.stats import rankdata
    xr, yr, zr = rankdata(x), rankdata(y), rankdata(z)
    # Residualize x and y on z (linear on ranks)
    def resid(a, b):
        A = np.vstack([b, np.ones_like(b)]).T
        coef, *_ = np.linalg.lstsq(A, a, rcond=None)
        return a - A @ coef
    xres = resid(xr, zr)
    yres = resid(yr, zr)
    return np.corrcoef(xres, yres)[0, 1]

print(f"{'Descriptor':<20} {'raw corr':>12} {'partial (|year)':>16} {'change'}")
print("-"*85)
for name in TOP_DESCS:
    vals = np.array([compute_desc(s, name) for s in sub['std_smiles']])
    valid = ~np.isnan(vals)
    raw = spearmanr(vals[valid], pchembl[valid]).correlation
    partial = partial_corr(vals[valid], pchembl[valid], year[valid])
    change = abs(raw - partial)
    verdict = "signal holds" if change < 0.1 else "partly year-driven"
    print(f"{name:<20} {raw:>12.3f} {partial:>16.3f}  {verdict}")

print("-"*85)
print("\nBOTTOM LINE:")
print("If partial correlations stay close to raw → SAR signal is real TRPA1 biology.")
print("If they drop sharply → part of the signal is publication-era medchem trend.")