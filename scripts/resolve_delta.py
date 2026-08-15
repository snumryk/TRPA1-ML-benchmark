"""
resolve_delta.py

Compare the saved ChEMBL 36 TRPA1 IC50 dataset with the current
ChEMBL REST API release.

Expected input file:
    trpa1_antagonists.csv

The script checks:
1. Current ChEMBL API database version.
2. New and removed standardized compounds by InChIKey.
3. Changes in aggregated pChEMBL values.
4. Changes in measurement, assay and document counts.
5. Presence of recent documents.

Important:
The saved CSV must really originate from ChEMBL 36.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import SaltRemover


# ============================================================
# Configuration
# ============================================================

TRPA1_HUMAN = "CHEMBL6007"

V36_FILE = Path("trpa1_antagonists.csv")
V37_AGGREGATED_FILE = Path("trpa1_current_api_aggregated.csv")
V37_RAW_FILE = Path("trpa1_current_api_raw.csv")
CHANGED_FILE = Path("trpa1_v36_v37_changed.csv")
NEW_FILE = Path("trpa1_v37_new_compounds.csv")
REMOVED_FILE = Path("trpa1_v37_removed_compounds.csv")
METADATA_FILE = Path("trpa1_current_api_metadata.json")

STATUS_URL = "https://www.ebi.ac.uk/chembl/api/data/status.json"

remover = SaltRemover.SaltRemover()


# ============================================================
# Helpers
# ============================================================

def get_chembl_status() -> dict:
    """Return the current ChEMBL API status payload."""
    response = requests.get(STATUS_URL, timeout=30)
    response.raise_for_status()

    payload = response.json()

    if not isinstance(payload, dict):
        raise RuntimeError(
            f"Unexpected ChEMBL status payload type: {type(payload)}"
        )

    return payload


def standardize(smiles: str | None) -> tuple[str | None, str | None]:
    """
    Apply the same simplified standardization used for the original dataset:
    parse -> salt removal -> largest fragment -> canonical SMILES -> InChIKey.
    """
    if smiles is None or pd.isna(smiles):
        return None, None

    mol = Chem.MolFromSmiles(str(smiles))

    if mol is None:
        return None, None

    try:
        mol = remover.StripMol(
            mol,
            dontRemoveEverything=True,
        )

        fragments = Chem.GetMolFrags(
            mol,
            asMols=True,
            sanitizeFrags=True,
        )

        if len(fragments) > 1:
            mol = max(
                fragments,
                key=lambda fragment: fragment.GetNumHeavyAtoms(),
            )

        canonical_smiles = Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=True,
        )

        inchikey = Chem.MolToInchiKey(mol)

        return canonical_smiles, inchikey

    except Exception:
        return None, None


def numeric_series(df: pd.DataFrame, column: str) -> pd.Series:
    """Safely convert a DataFrame column to numeric values."""
    return pd.to_numeric(
        df[column],
        errors="coerce",
    )


def values_changed(
    left: pd.Series,
    right: pd.Series,
    tolerance: float = 1e-10,
) -> pd.Series:
    """Compare numeric series while treating paired NaN values as equal."""
    left_numeric = pd.to_numeric(left, errors="coerce")
    right_numeric = pd.to_numeric(right, errors="coerce")

    return ~np.isclose(
        left_numeric,
        right_numeric,
        rtol=0,
        atol=tolerance,
        equal_nan=True,
    )


# ============================================================
# 1. Check API release
# ============================================================

print("=" * 70)
print("ChEMBL API STATUS")
print("=" * 70)

status = get_chembl_status()

db_version = status.get("chembl_db_version")
release_date = status.get("chembl_release_date")

print(f"Database version: {db_version}")
print(f"Release date:     {release_date}")
print(f"API status:       {status.get('status')}")
print(f"Activities:       {status.get('activities')}")
print(f"Query time UTC:   {datetime.now(timezone.utc).isoformat()}")

if db_version != "ChEMBL_37":
    raise RuntimeError(
        f"Expected ChEMBL_37, but API currently reports {db_version}."
    )


# ============================================================
# 2. Load saved ChEMBL 36 compound dataset
# ============================================================

print("\n" + "=" * 70)
print("LOADING SAVED ChEMBL 36 DATASET")
print("=" * 70)

if not V36_FILE.exists():
    raise FileNotFoundError(
        f"Input file not found: {V36_FILE.resolve()}"
    )

df_v36 = pd.read_csv(V36_FILE)

required_v36_columns = {
    "inchikey",
    "std_smiles",
    "pchembl_median",
}

missing_v36_columns = required_v36_columns - set(df_v36.columns)

if missing_v36_columns:
    raise ValueError(
        "The saved v36 file is missing required columns: "
        f"{sorted(missing_v36_columns)}"
    )

df_v36["inchikey"] = (
    df_v36["inchikey"]
    .astype(str)
    .str.strip()
)

df_v36 = df_v36[
    df_v36["inchikey"].notna()
    & (df_v36["inchikey"] != "")
    & (df_v36["inchikey"].str.lower() != "nan")
].copy()

duplicate_v36_keys = df_v36["inchikey"].duplicated().sum()

if duplicate_v36_keys:
    raise ValueError(
        f"The saved v36 file contains {duplicate_v36_keys} "
        "duplicate InChIKeys. Expected one row per compound."
    )

print(f"File:                   {V36_FILE.resolve()}")
print(f"Rows:                   {len(df_v36)}")
print(f"Unique InChIKeys:       {df_v36['inchikey'].nunique()}")
print(f"Columns:                {list(df_v36.columns)}")


# ============================================================
# 3. Pull current ChEMBL 37 activities
# ============================================================

print("\n" + "=" * 70)
print("PULLING ChEMBL 37 TRPA1 IC50 ACTIVITIES")
print("=" * 70)

activity = new_client.activity

query = activity.filter(
    target_chembl_id=TRPA1_HUMAN,
    standard_type="IC50",
    standard_relation="=",
    pchembl_value__isnull=False,
).only([
    "activity_id",
    "molecule_chembl_id",
    "canonical_smiles",
    "pchembl_value",
    "standard_value",
    "standard_units",
    "standard_relation",
    "standard_type",
    "assay_type",
    "assay_chembl_id",
    "document_chembl_id",
    "document_year",
    "data_validity_comment",
])

activities = list(query)

if not activities:
    raise RuntimeError(
        "The API returned no TRPA1 IC50 activities."
    )

df_raw = pd.DataFrame(activities)

print(f"Raw activities:         {len(df_raw)}")
print(
    "Raw molecule IDs:      "
    f"{df_raw['molecule_chembl_id'].nunique()}"
)
print(f"Returned columns:       {list(df_raw.columns)}")


# ============================================================
# 4. Clean and standardize current records
# ============================================================

df_raw["pchembl_value"] = pd.to_numeric(
    df_raw["pchembl_value"],
    errors="coerce",
)

df_raw["document_year"] = pd.to_numeric(
    df_raw["document_year"],
    errors="coerce",
)

df_raw = df_raw.dropna(
    subset=[
        "pchembl_value",
        "canonical_smiles",
    ]
).copy()

standardized = df_raw["canonical_smiles"].apply(standardize)

df_raw["std_smiles"] = standardized.apply(lambda value: value[0])
df_raw["inchikey"] = standardized.apply(lambda value: value[1])

failed_standardization = df_raw["inchikey"].isna().sum()

df_raw = df_raw.dropna(
    subset=[
        "std_smiles",
        "inchikey",
    ]
).copy()

print(f"Failed standardization: {failed_standardization}")
print(f"Usable raw activities: {len(df_raw)}")

df_raw.to_csv(
    V37_RAW_FILE,
    index=False,
)

print(f"Saved raw v37 records:  {V37_RAW_FILE.resolve()}")


# ============================================================
# 5. Aggregate ChEMBL 37 by standardized InChIKey
# ============================================================

df_v37 = (
    df_raw.groupby(
        "inchikey",
        as_index=False,
    )
    .agg(
        std_smiles=("std_smiles", "first"),
        molecule_chembl_id=("molecule_chembl_id", "first"),
        pchembl_median=("pchembl_value", "median"),
        pchembl_min=("pchembl_value", "min"),
        pchembl_max=("pchembl_value", "max"),
        pchembl_std=("pchembl_value", "std"),
        n_measurements=("pchembl_value", "count"),
        year_min=("document_year", "min"),
        year_max=("document_year", "max"),
        n_documents=("document_chembl_id", "nunique"),
        n_assays=("assay_chembl_id", "nunique"),
    )
)

df_v37.to_csv(
    V37_AGGREGATED_FILE,
    index=False,
)

print("\nChEMBL 37 aggregated dataset:")
print(f"Unique compounds:       {len(df_v37)}")
print(f"Measurements:           {df_v37['n_measurements'].sum()}")
print(f"Unique assays:          {df_raw['assay_chembl_id'].nunique()}")
print(
    f"Unique documents:       "
    f"{df_raw['document_chembl_id'].nunique()}"
)
print(
    f"Document year range:    "
    f"{df_raw['document_year'].min()} – "
    f"{df_raw['document_year'].max()}"
)
print(
    f"Saved aggregated data:  "
    f"{V37_AGGREGATED_FILE.resolve()}"
)


# ============================================================
# 6. Compare compound identities
# ============================================================

v36_keys = set(df_v36["inchikey"])
v37_keys = set(df_v37["inchikey"])

new_keys = v37_keys - v36_keys
removed_keys = v36_keys - v37_keys
overlap_keys = v36_keys & v37_keys

print("\n" + "=" * 70)
print("COMPOUND-LEVEL DELTA: ChEMBL 36 -> ChEMBL 37")
print("=" * 70)

print(f"v36 compounds:          {len(v36_keys)}")
print(f"v37 compounds:          {len(v37_keys)}")
print(f"Overlap:                {len(overlap_keys)}")
print(f"New in v37:             {len(new_keys)}")
print(f"Removed from v37:       {len(removed_keys)}")

df_new = df_v37[
    df_v37["inchikey"].isin(new_keys)
].copy()

df_removed = df_v36[
    df_v36["inchikey"].isin(removed_keys)
].copy()

df_new.to_csv(NEW_FILE, index=False)
df_removed.to_csv(REMOVED_FILE, index=False)

print(f"Saved new compounds:    {NEW_FILE.resolve()}")
print(f"Saved removed compounds:{REMOVED_FILE.resolve()}")


# ============================================================
# 7. Compare aggregate values for overlapping compounds
# ============================================================

candidate_columns = [
    "pchembl_median",
    "pchembl_min",
    "pchembl_max",
    "pchembl_std",
    "n_measurements",
    "year_min",
    "year_max",
    "n_documents",
    "n_assays",
]

common_columns = [
    column
    for column in candidate_columns
    if column in df_v36.columns
    and column in df_v37.columns
]

print("\nComparable aggregate columns:")
print(common_columns)

comparison = df_v36[
    ["inchikey", "std_smiles"] + common_columns
].merge(
    df_v37[
        ["inchikey", "std_smiles"] + common_columns
    ],
    on="inchikey",
    how="outer",
    suffixes=("_v36", "_v37"),
    indicator=True,
)

comparison["any_aggregate_change"] = False

print("\n" + "=" * 70)
print("AGGREGATE-VALUE CHANGES")
print("=" * 70)

for column in common_columns:
    changed_column = f"{column}_changed"

    comparison[changed_column] = values_changed(
        comparison[f"{column}_v36"],
        comparison[f"{column}_v37"],
    )

    comparison["any_aggregate_change"] |= comparison[changed_column]

    changed_count = int(
        comparison.loc[
            comparison["_merge"] == "both",
            changed_column,
        ].sum()
    )

    print(f"{column:20s}: {changed_count} compounds changed")

changed_overlap = comparison[
    (comparison["_merge"] == "both")
    & comparison["any_aggregate_change"]
].copy()

changed_overlap.to_csv(
    CHANGED_FILE,
    index=False,
)

print(
    f"\nOverlapping compounds with any changed aggregate: "
    f"{len(changed_overlap)}"
)
print(f"Saved changed records:  {CHANGED_FILE.resolve()}")


# ============================================================
# 8. Detailed measurement and potency changes
# ============================================================

if "n_measurements" in common_columns:
    old_count = numeric_series(
        comparison,
        "n_measurements_v36",
    )
    new_count = numeric_series(
        comparison,
        "n_measurements_v37",
    )

    overlap_mask = comparison["_merge"] == "both"

    print("\nMeasurement count changes:")
    print(
        "Increased:             "
        f"{int(((new_count > old_count) & overlap_mask).sum())}"
    )
    print(
        "Decreased:             "
        f"{int(((new_count < old_count) & overlap_mask).sum())}"
    )
    print(
        "Unchanged:             "
        f"{int(((new_count == old_count) & overlap_mask).sum())}"
    )

    print("\nTotal measurements reconstructed from aggregates:")
    print(
        "v36:                   "
        f"{numeric_series(df_v36, 'n_measurements').sum()}"
    )
    print(
        "v37:                   "
        f"{numeric_series(df_v37, 'n_measurements').sum()}"
    )


if "pchembl_median" in common_columns:
    old_median = numeric_series(
        comparison,
        "pchembl_median_v36",
    )
    new_median = numeric_series(
        comparison,
        "pchembl_median_v37",
    )

    overlap_mask = comparison["_merge"] == "both"
    delta = new_median - old_median
    overlap_delta = delta[overlap_mask]

    print("\nMedian pChEMBL delta:")
    print(overlap_delta.describe())

    print(
        "|delta| > 0.01:        "
        f"{int((overlap_delta.abs() > 0.01).sum())}"
    )
    print(
        "|delta| > 0.10:        "
        f"{int((overlap_delta.abs() > 0.10).sum())}"
    )
    print(
        "|delta| > 0.50:        "
        f"{int((overlap_delta.abs() > 0.50).sum())}"
    )


# ============================================================
# 9. Recent document diagnostic
# ============================================================

recent = df_raw[
    df_raw["document_year"] >= 2025
].copy()

print("\n" + "=" * 70)
print("RECENT-DOCUMENT DIAGNOSTIC")
print("=" * 70)

print(f"Activities from 2025+:  {len(recent)}")
print(
    "Compounds from 2025+:  "
    f"{recent['inchikey'].nunique()}"
)
print(
    "Documents from 2025+:  "
    f"{recent['document_chembl_id'].nunique()}"
)
print(
    "Maximum document year: "
    f"{df_raw['document_year'].max()}"
)

if not recent.empty:
    recent_documents = (
        recent[
            [
                "document_chembl_id",
                "document_year",
                "assay_chembl_id",
                "inchikey",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            [
                "document_year",
                "document_chembl_id",
            ]
        )
    )

    print("\nRecent records:")
    print(recent_documents.to_string(index=False))


# ============================================================
# 10. Save metadata and print verdict
# ============================================================

metadata = {
    "extraction_time_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "chembl_db_version": db_version,
    "chembl_release_date": release_date,
    "target_chembl_id": TRPA1_HUMAN,
    "standard_type": "IC50",
    "standard_relation": "=",
    "requires_pchembl": True,
    "raw_activity_count_v37": int(len(df_raw)),
    "compound_count_v36": int(len(v36_keys)),
    "compound_count_v37": int(len(v37_keys)),
    "overlap_compounds": int(len(overlap_keys)),
    "new_compounds_v37": int(len(new_keys)),
    "removed_compounds_v37": int(len(removed_keys)),
    "changed_overlapping_compounds": int(len(changed_overlap)),
    "maximum_document_year": (
        None
        if pd.isna(df_raw["document_year"].max())
        else int(df_raw["document_year"].max())
    ),
}

METADATA_FILE.write_text(
    json.dumps(
        metadata,
        indent=2,
        ensure_ascii=False,
    ),
    encoding="utf-8",
)

print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

if len(new_keys) > 0:
    print(
        f"ChEMBL 37 contains {len(new_keys)} new standardized "
        "TRPA1 IC50 compounds."
    )
    print(
        "These compounds may be considered for a release-based "
        "temporal holdout after additional assay-quality checks."
    )

elif len(changed_overlap) > 0:
    print(
        "No new standardized TRPA1 IC50 compounds were found, "
        f"but {len(changed_overlap)} existing compounds have "
        "changed aggregate data."
    )
    print(
        "This is not an independent compound-level temporal "
        "holdout because the structures already existed in v36."
    )

else:
    print(
        "No new standardized compounds and no aggregate-level "
        "changes were detected for this strict human TRPA1 "
        "IC50 subset."
    )
    print(
        "For this filtered subset, ChEMBL 36 and ChEMBL 37 "
        "appear identical."
    )

print(f"\nMetadata saved to:      {METADATA_FILE.resolve()}")