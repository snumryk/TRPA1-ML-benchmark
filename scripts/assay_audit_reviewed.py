"""
assay_audit_reviewed.py

Build a reproducible, review-ready assay audit table for the curated human
TRPA1 IC50 activity dataset.

The script does NOT decide that an assay is mechanistically homogeneous.
Automatic keyword flags only prioritize assays for manual review.

Typical usage (PowerShell):
    python assay_audit_reviewed.py \
        --input trpa1_current_api_raw.csv \
        --output trpa1_assay_audit.csv

Outputs:
    trpa1_assay_audit.csv             one row per assay, ready for manual review
    trpa1_assay_audit_summary.csv     automatic flag coverage summary
    trpa1_assay_audit_metadata.json   provenance and run metadata
    trpa1_assay_cache.json            cached raw ChEMBL assay records

Requirements:
    pandas
    requests
    chembl_webresource_client
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import unicodedata
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd
import requests
from chembl_webresource_client.new_client import new_client


STATUS_URL = "https://www.ebi.ac.uk/chembl/api/data/status.json"
DEFAULT_INPUT = "trpa1_current_api_raw.csv"
DEFAULT_OUTPUT = "trpa1_assay_audit.csv"
DEFAULT_SUMMARY = "trpa1_assay_audit_summary.csv"
DEFAULT_METADATA = "trpa1_assay_audit_metadata.json"
DEFAULT_CACHE = "trpa1_assay_cache.json"

REQUIRED_INPUT_COLUMNS = {"assay_chembl_id", "inchikey"}

# Each rule is (human-readable label, regular expression).
# Matching is case-insensitive because text is normalized to lowercase first.
KEYWORD_RULES: dict[str, list[tuple[str, str]]] = {
    # Readout / modality
    "calcium": [
        ("calcium", r"\bcalcium\b"),
        ("calcium_flux", r"\bcalcium[ -]flux\b"),
        ("calcium_mobilization", r"\bcalcium[ -]mobili[sz]ation\b"),
        ("intracellular_ca", r"\bintracellular\s+ca(?:2\+|\+2)?\b"),
        ("ca2+", r"\bca\s*(?:2\+|\+2)\b"),
        ("fluo-4", r"\bfluo[ -]?4\b"),
        ("fura-2", r"\bfura[ -]?2\b"),
        ("flipr", r"\bflipr\b"),
    ],
    "electrophysiology": [
        ("patch_clamp", r"\bpatch[ -]?clamp\b"),
        ("electrophysiology", r"\belectrophysiolog\w*\b"),
        ("whole_cell", r"\bwhole[ -]?cell\b"),
        ("voltage_clamp", r"\bvoltage[ -]?clamp\b"),
        ("membrane_current", r"\bmembrane\s+current\b"),
        ("ionic_current", r"\bionic\s+current\b"),
        ("current_amplitude", r"\bcurrent\s+amplitude\b"),
        ("holding_potential", r"\bholding\s+potential\b"),
    ],
    "fluorescence": [
        ("fluorescence", r"\bfluorescen\w*\b"),
        ("fluorometric", r"\bfluorometric\b"),
        ("fluo-4", r"\bfluo[ -]?4\b"),
        ("fura-2", r"\bfura[ -]?2\b"),
        ("flipr", r"\bflipr\b"),
    ],
    "binding": [
        ("binding", r"\bbinding\b"),
        ("radioligand", r"\bradioligand\b"),
        ("displacement", r"\bdisplacement\b"),
    ],

    # Challenge agonists / activators
    "aitc": [
        ("AITC", r"\baitc\b"),
        ("allyl_isothiocyanate", r"\ballyl\s+isothiocyanate\b"),
        ("mustard_oil", r"\bmustard\s+oil\b"),
    ],
    "cinnamaldehyde": [
        ("cinnamaldehyde", r"\bcinnamaldehyde\b"),
        ("cinnamald", r"\bcinnamald\w*\b"),
    ],
    "jt010": [
        ("JT010", r"\bjt[ -]?010\b"),
    ],
    "other_named_agonist": [
        ("acrolein", r"\bacrolein\b"),
        ("allicin", r"\ballicin\b"),
        ("methylglyoxal", r"\bmethylglyoxal\b"),
        ("4-HNE", r"\b(?:4[ -]?)?hydroxynonenal\b|\b4[ -]?hne\b"),
        ("carvacrol", r"\bcarvacrol\b"),
        ("thymol", r"\bthymol\b"),
        ("menthol", r"\bmenthol\b"),
    ],
    "agonist_context": [
        ("agonist", r"\bagonist\w*\b"),
        ("activation", r"\bactivat(?:e|ed|es|ing|ion)\b"),
        ("evoked", r"\bevoked\b"),
        ("induced", r"\binduced\b"),
        ("stimulated", r"\bstimulat(?:e|ed|es|ing|ion)\b"),
    ],

    # Protocol order / pharmacological direction
    "preincubation": [
        ("preincubation", r"\bpre[ -]?incubat\w*\b"),
        ("pretreatment", r"\bpre[ -]?treat\w*\b"),
        ("before_agonist", r"\bbefore\s+(?:agonist|aitc|allyl\s+isothiocyanate)\b"),
    ],
    "coapplication": [
        ("coapplication", r"\bco[ -]?applicat\w*\b"),
        ("simultaneous", r"\bsimultaneous(?:ly)?\b"),
        ("concomitant", r"\bconcomitant(?:ly)?\b"),
    ],
    "postactivation": [
        ("after_activation", r"\bafter\s+(?:channel\s+)?activat\w*\b"),
        ("following_activation", r"\bfollowing\s+(?:channel\s+)?activat\w*\b"),
        ("post_activation", r"\bpost[ -]?activat\w*\b"),
        ("sustained_current", r"\bsustained\s+current\b"),
    ],
    "inhibition": [
        ("antagonist", r"\bantagonis(?:t|ts|m|tic)\b"),
        ("inhibition", r"\binhibit\w*\b"),
        ("block", r"\bblock(?:ed|er|ers|ing|ade|ades)?\b"),
        ("IC50", r"\bic\s*50\b"),
    ],
    "direct_agonism": [
        ("agonist_activity", r"\bagonist\s+activity\b"),
        ("as_an_agonist", r"\bas\s+an?\s+agonist\b"),
        ("compound_induced_activation", r"\b(?:compound|test\s+compound)[ -]induced\s+activat\w*\b"),
        ("activation_by_compound", r"\bactivat\w*\s+by\s+(?:the\s+)?(?:compound|test\s+compound)\b"),
    ],

    # Construct / species / system
    "mutant_or_chimera": [
        ("mutant", r"\bmutant\w*\b"),
        ("mutation", r"\bmutation\w*\b"),
        ("chimera", r"\bchimer\w*\b"),
        ("variant", r"\bvariant\b"),
        ("C621", r"\bc621\b"),
        ("C641", r"\bc641\b"),
        ("C665", r"\bc665\b"),
        ("N855", r"\bn855\b"),
    ],
    "wild_type": [
        ("wild_type", r"\bwild[ -]?type\b"),
        ("WT", r"\bwt\b"),
    ],
    "human": [
        ("human", r"\bhuman\b"),
        ("hTRPA1", r"\bhtrpa1\b"),
    ],
    "hek293": [
        ("HEK293", r"\bhek[ -]?293\w*\b"),
    ],
    "cho": [
        ("CHO", r"\bcho(?:[ -]?k1)?\b"),
    ],
    "xenopus": [
        ("Xenopus", r"\bxenopus\b"),
        ("oocyte", r"\boocyte\w*\b"),
    ],
}

COMPILED_RULES: dict[str, list[tuple[str, re.Pattern[str]]]] = {
    category: [(label, re.compile(pattern)) for label, pattern in rules]
    for category, rules in KEYWORD_RULES.items()
}

# Fields copied from the ChEMBL assay record when available. Fetching the full
# record rather than using .only([...]) makes the script more tolerant of API
# schema changes; absent fields remain blank.
ASSAY_FIELDS = [
    "assay_chembl_id",
    "description",
    "assay_type",
    "assay_category",
    "assay_test_type",
    "assay_organism",
    "assay_tax_id",
    "assay_strain",
    "assay_tissue",
    "assay_cell_type",
    "assay_subcellular_fraction",
    "bao_format",
    "bao_label",
    "confidence_score",
    "relationship_type",
    "target_chembl_id",
    "document_chembl_id",
    "cell_chembl_id",
    "tissue_chembl_id",
    "variant_id",
    "assay_variant_accession",
    "assay_variant_mutation",
]

MANUAL_COLUMNS_DEFAULTS: dict[str, str] = {
    "review_status": "unreviewed",
    "manual_assay_modality": "",
    "manual_functional_direction": "",
    "manual_challenge_agonist": "",
    "manual_application_order": "",
    "manual_readout": "",
    "manual_expression_system": "",
    "manual_channel_construct": "",
    "manual_mechanistic_interpretability": "",
    "manual_desensitization_risk": "",
    "manual_include_strict_subset": "",
    "manual_exclusion_reason": "",
    "manual_notes": "",
}


class AuditError(RuntimeError):
    """Raised when the audit cannot be completed reliably."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a review-ready assay audit table for a raw TRPA1 activity CSV. "
            "Automatic flags are screening aids, not final biological labels."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Raw activity CSV")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Assay audit CSV")
    parser.add_argument("--summary", default=DEFAULT_SUMMARY, help="Flag summary CSV")
    parser.add_argument("--metadata", default=DEFAULT_METADATA, help="Run metadata JSON")
    parser.add_argument("--cache", default=DEFAULT_CACHE, help="ChEMBL assay JSON cache")
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Ignore cached assay records and fetch all assays again",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="API attempts per assay (default: 3)",
    )
    parser.add_argument(
        "--retry-delay",
        type=float,
        default=1.0,
        help="Initial retry delay in seconds; exponential backoff is used",
    )
    parser.add_argument(
        "--allow-api-errors",
        action="store_true",
        help="Write partial output and exit successfully even if some assays failed",
    )
    parser.add_argument(
        "--manual-input",
        default=None,
        help=(
            "Optional previous audit CSV whose manual_* columns should be preserved. "
            "If omitted and --output already exists, that file is used automatically."
        ),
    )
    parser.add_argument(
        "--reset-manual",
        action="store_true",
        help="Do not preserve manual annotations from an existing audit CSV",
    )
    return parser.parse_args()


def normalize_text(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = unicodedata.normalize("NFKC", str(value)).lower()
    text = text.replace("−", "-").replace("–", "-").replace("—", "-")
    return re.sub(r"\s+", " ", text).strip()


def match_rule(text: str, category: str) -> list[str]:
    matches: list[str] = []
    for label, pattern in COMPILED_RULES[category]:
        if pattern.search(text):
            matches.append(label)
    return sorted(set(matches))


def join_unique(values: Iterable[Any]) -> str:
    clean = {
        str(value).strip()
        for value in values
        if pd.notna(value) and str(value).strip() and str(value).strip().lower() != "nan"
    }
    return ";".join(sorted(clean))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_raw_activities(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path.resolve()}")

    raw = pd.read_csv(path, low_memory=False)
    missing = REQUIRED_INPUT_COLUMNS - set(raw.columns)
    if missing:
        raise AuditError(
            "Input CSV is missing required columns: " + ", ".join(sorted(missing))
        )

    raw = raw.copy()
    raw["assay_chembl_id"] = raw["assay_chembl_id"].astype("string").str.strip()
    raw["inchikey"] = raw["inchikey"].astype("string").str.strip()

    invalid_assay = raw["assay_chembl_id"].isna() | raw["assay_chembl_id"].eq("")
    if invalid_assay.any():
        print(f"WARNING: dropping {int(invalid_assay.sum())} rows without assay_chembl_id")
        raw = raw.loc[~invalid_assay].copy()

    invalid_key = raw["inchikey"].isna() | raw["inchikey"].eq("")
    if invalid_key.any():
        print(
            "WARNING: "
            f"{int(invalid_key.sum())} rows have no InChIKey; they remain in measurement "
            "counts but cannot contribute to unique-compound counts."
        )

    for numeric_column in ("pchembl_value", "document_year"):
        if numeric_column in raw.columns:
            raw[numeric_column] = pd.to_numeric(raw[numeric_column], errors="coerce")

    return raw


def aggregate_raw_by_assay(raw: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for assay_id, group in raw.groupby("assay_chembl_id", sort=True, dropna=False):
        row: dict[str, Any] = {
            "assay_chembl_id": str(assay_id),
            "n_measurements": int(len(group)),
            "n_compounds": int(group["inchikey"].nunique(dropna=True)),
        }

        optional_nunique = {
            "n_activity_ids": "activity_id",
            "n_molecule_chembl_ids": "molecule_chembl_id",
            "n_documents_raw": "document_chembl_id",
        }
        for output_column, source_column in optional_nunique.items():
            if source_column in group.columns:
                row[output_column] = int(group[source_column].nunique(dropna=True))

        optional_join = {
            "document_chembl_ids_raw": "document_chembl_id",
            "assay_types_raw": "assay_type",
            "data_validity_comments_raw": "data_validity_comment",
        }
        for output_column, source_column in optional_join.items():
            if source_column in group.columns:
                row[output_column] = join_unique(group[source_column])

        if "pchembl_value" in group.columns:
            values = group["pchembl_value"].dropna()
            row.update(
                {
                    "pchembl_n": int(values.count()),
                    "pchembl_min": float(values.min()) if not values.empty else None,
                    "pchembl_median": float(values.median()) if not values.empty else None,
                    "pchembl_max": float(values.max()) if not values.empty else None,
                    "pchembl_mean": float(values.mean()) if not values.empty else None,
                    "pchembl_std": float(values.std(ddof=1)) if len(values) >= 2 else None,
                }
            )

        if "document_year" in group.columns:
            years = group["document_year"].dropna()
            row["year_min_raw"] = int(years.min()) if not years.empty else None
            row["year_max_raw"] = int(years.max()) if not years.empty else None
            row["document_years_raw"] = ";".join(
                str(int(value)) for value in sorted(years.unique())
            )

        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["n_compounds", "n_measurements", "assay_chembl_id"],
        ascending=[False, False, True],
    )


def get_chembl_status(timeout: float = 30.0) -> dict[str, Any]:
    response = requests.get(STATUS_URL, timeout=timeout)
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, dict):
        raise AuditError(f"Unexpected ChEMBL status payload: {type(payload)}")
    return payload


def load_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuditError(f"Cannot read cache {path}: {exc}") from exc

    if not isinstance(payload, dict):
        raise AuditError(f"Cache must contain a JSON object: {path}")
    return {
        str(key): value
        for key, value in payload.items()
        if isinstance(value, dict)
    }


def save_cache(path: Path, cache: dict[str, dict[str, Any]]) -> None:
    path.write_text(
        json.dumps(cache, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )


def fetch_one_assay(
    assay_id: str,
    retries: int,
    retry_delay: float,
) -> dict[str, Any]:
    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            # Fetch the complete record. This avoids making the script fail when an
            # optional field name is absent from a future API version.
            result = list(new_client.assay.filter(assay_chembl_id=assay_id))
            if len(result) != 1:
                raise AuditError(
                    f"Expected exactly one assay record for {assay_id}, got {len(result)}"
                )
            record = dict(result[0])
            record["assay_chembl_id"] = assay_id
            return record
        except Exception as exc:  # API/network boundary: retry and preserve details.
            last_error = exc
            if attempt < retries:
                time.sleep(retry_delay * (2 ** (attempt - 1)))

    raise AuditError(f"Failed to retrieve {assay_id}: {last_error}")


def fetch_assays(
    assay_ids: Sequence[str],
    cache_path: Path,
    refresh_cache: bool,
    retries: int,
    retry_delay: float,
) -> tuple[pd.DataFrame, list[str]]:
    cache = {} if refresh_cache else load_cache(cache_path)
    records: list[dict[str, Any]] = []
    errors: list[str] = []

    total = len(assay_ids)
    for index, assay_id in enumerate(assay_ids, start=1):
        if assay_id in cache and not refresh_cache:
            record = dict(cache[assay_id])
            source = "cache"
        else:
            try:
                record = fetch_one_assay(assay_id, retries, retry_delay)
                cache[assay_id] = record
                save_cache(cache_path, cache)  # checkpoint after every successful fetch
                source = "api"
            except AuditError as exc:
                message = str(exc)
                errors.append(message)
                record = {
                    "assay_chembl_id": assay_id,
                    "description": "",
                    "api_error": message,
                }
                source = "error"

        selected = {field: record.get(field) for field in ASSAY_FIELDS}
        selected["assay_chembl_id"] = assay_id
        selected["api_fetch_source"] = source
        selected["api_error"] = record.get("api_error", "")
        records.append(selected)

        if index == 1 or index % 10 == 0 or index == total:
            print(f"  assays processed: {index}/{total}")

    return pd.DataFrame(records), errors


def annotate_keywords(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()

    # Search biological phrases in a combined field, while preserving the original
    # description separately for manual review.
    text_columns = [
        "description",
        "assay_organism",
        "assay_cell_type",
        "assay_tissue",
        "bao_label",
        "assay_variant_mutation",
    ]
    available = [column for column in text_columns if column in result.columns]
    combined = result[available].fillna("").astype(str).agg(" | ".join, axis=1)
    result["normalized_search_text"] = combined.map(normalize_text)

    for category in COMPILED_RULES:
        matches = result["normalized_search_text"].map(
            lambda text, current=category: match_rule(text, current)
        )
        result[f"flag_{category}"] = matches.map(bool)
        result[f"matches_{category}"] = matches.map(lambda values: ";".join(values))

    description = result.get("description", pd.Series("", index=result.index)).fillna("")
    normalized_description = description.map(normalize_text)
    result["description_length"] = normalized_description.str.len()
    result["flag_missing_or_sparse_description"] = result["description_length"] < 25

    assay_type = result.get("assay_type", pd.Series("", index=result.index)).fillna("")
    bao_label = result.get("bao_label", pd.Series("", index=result.index)).fillna("")
    functional_metadata = (
        assay_type.astype(str).str.upper().eq("F")
        | bao_label.astype(str).str.contains("functional", case=False, na=False)
    )
    result["flag_functional_metadata"] = functional_metadata

    # Candidate flags are intentionally conservative in naming. They indicate where
    # manual review should start; they do not establish mechanistic homogeneity.
    result["candidate_calcium_AITC_review"] = (
        result["flag_calcium"] & result["flag_aitc"]
    )
    result["candidate_functional_inhibition_review"] = (
        result["flag_functional_metadata"] & result["flag_inhibition"]
    )
    result["candidate_strict_review"] = (
        result["candidate_calcium_AITC_review"]
        & result["candidate_functional_inhibition_review"]
        & ~result["flag_mutant_or_chimera"]
    )

    # Priority is a triage tool. Missing API data and ambiguous descriptions come first,
    # followed by large candidate assays.
    result["auto_review_priority"] = "normal"
    result.loc[
        result["api_fetch_source"].eq("error")
        | result["flag_missing_or_sparse_description"],
        "auto_review_priority",
    ] = "high"
    result.loc[
        result["candidate_strict_review"]
        & ~result["flag_missing_or_sparse_description"]
        & ~result["api_fetch_source"].eq("error"),
        "auto_review_priority",
    ] = "candidate_strict"

    for column, default in MANUAL_COLUMNS_DEFAULTS.items():
        if column not in result.columns:
            result[column] = default

    return result



def preserve_manual_annotations(
    audit: pd.DataFrame,
    previous_path: Path | None,
) -> pd.DataFrame:
    """Merge prior manual review columns by assay_chembl_id."""
    if previous_path is None or not previous_path.exists():
        return audit

    previous = pd.read_csv(previous_path, low_memory=False)
    if "assay_chembl_id" not in previous.columns:
        raise AuditError(
            f"Previous audit lacks assay_chembl_id: {previous_path.resolve()}"
        )

    previous = previous.copy()
    previous["assay_chembl_id"] = (
        previous["assay_chembl_id"].astype("string").str.strip()
    )
    if previous["assay_chembl_id"].duplicated().any():
        raise AuditError(
            f"Previous audit contains duplicate assay_chembl_id values: "
            f"{previous_path.resolve()}"
        )

    available_manual = [
        column for column in MANUAL_COLUMNS_DEFAULTS if column in previous.columns
    ]
    if not available_manual:
        return audit

    prior_manual = previous[["assay_chembl_id", *available_manual]].copy()
    result = audit.drop(columns=available_manual, errors="ignore").merge(
        prior_manual,
        on="assay_chembl_id",
        how="left",
        validate="one_to_one",
    )

    for column, default in MANUAL_COLUMNS_DEFAULTS.items():
        if column not in result.columns:
            result[column] = default
        else:
            result[column] = result[column].fillna(default)

    print(
        f"Preserved manual annotations from: {previous_path.resolve()} "
        f"({len(available_manual)} columns)"
    )
    return result

def unique_compound_coverage(raw: pd.DataFrame, assay_ids: Iterable[str]) -> int:
    assay_set = {str(value) for value in assay_ids}
    return int(
        raw.loc[raw["assay_chembl_id"].isin(assay_set), "inchikey"].nunique(dropna=True)
    )


def build_summary(audit: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    categories = [
        "calcium",
        "electrophysiology",
        "fluorescence",
        "binding",
        "aitc",
        "cinnamaldehyde",
        "jt010",
        "other_named_agonist",
        "agonist_context",
        "preincubation",
        "coapplication",
        "postactivation",
        "inhibition",
        "direct_agonism",
        "mutant_or_chimera",
        "wild_type",
        "human",
        "hek293",
        "cho",
        "xenopus",
        "functional_metadata",
        "missing_or_sparse_description",
        "candidate_calcium_AITC_review",
        "candidate_functional_inhibition_review",
        "candidate_strict_review",
    ]

    rows: list[dict[str, Any]] = []
    for category in categories:
        column = category if category.startswith("candidate_") else f"flag_{category}"
        if column not in audit.columns:
            continue
        mask = audit[column].fillna(False).astype(bool)
        selected = audit.loc[mask]
        rows.append(
            {
                "category": category,
                "n_assays": int(mask.sum()),
                "sum_assay_compound_counts_nonunique": int(
                    selected["n_compounds"].fillna(0).sum()
                ),
                "unique_compounds_union": unique_compound_coverage(
                    raw, selected["assay_chembl_id"]
                ),
                "n_measurements": int(selected["n_measurements"].fillna(0).sum()),
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["n_assays", "unique_compounds_union", "category"],
        ascending=[False, False, True],
    )


def package_version(package_name: str) -> str | None:
    try:
        return importlib_metadata.version(package_name)
    except importlib_metadata.PackageNotFoundError:
        return None


def ordered_columns(df: pd.DataFrame) -> list[str]:
    preferred = [
        "assay_chembl_id",
        "n_compounds",
        "n_measurements",
        "n_activity_ids",
        "n_molecule_chembl_ids",
        "pchembl_n",
        "pchembl_min",
        "pchembl_median",
        "pchembl_max",
        "pchembl_mean",
        "pchembl_std",
        "description",
        "description_length",
        "assay_type",
        "bao_label",
        "assay_category",
        "assay_test_type",
        "assay_organism",
        "assay_cell_type",
        "assay_tissue",
        "confidence_score",
        "target_chembl_id",
        "document_chembl_id",
        "document_chembl_ids_raw",
        "document_years_raw",
        "year_min_raw",
        "year_max_raw",
        "assay_variant_accession",
        "assay_variant_mutation",
        "api_fetch_source",
        "api_error",
        "auto_review_priority",
        "candidate_strict_review",
        "candidate_calcium_AITC_review",
        "candidate_functional_inhibition_review",
        "flag_functional_metadata",
        "flag_missing_or_sparse_description",
    ]

    flag_columns = sorted(
        column
        for column in df.columns
        if column.startswith("flag_") and column not in preferred
    )
    match_columns = sorted(column for column in df.columns if column.startswith("matches_"))
    manual_columns = list(MANUAL_COLUMNS_DEFAULTS)
    link_columns = ["assay_api_url", "document_api_url"]

    chosen: list[str] = []
    for column in preferred + flag_columns + match_columns + manual_columns + link_columns:
        if column in df.columns and column not in chosen:
            chosen.append(column)

    remaining = [column for column in df.columns if column not in chosen]
    return chosen + remaining


def main() -> int:
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    summary_path = Path(args.summary)
    metadata_path = Path(args.metadata)
    cache_path = Path(args.cache)

    if args.retries < 1:
        raise AuditError("--retries must be at least 1")
    if args.retry_delay < 0:
        raise AuditError("--retry-delay cannot be negative")

    print("=" * 76)
    print("TRPA1 ASSAY AUDIT")
    print("=" * 76)

    raw = load_raw_activities(input_path)
    assay_counts = aggregate_raw_by_assay(raw)
    assay_ids = assay_counts["assay_chembl_id"].astype(str).tolist()

    print(f"Input:                 {input_path.resolve()}")
    print(f"Activity rows:         {len(raw)}")
    print(f"Unique compounds:      {raw['inchikey'].nunique(dropna=True)}")
    print(f"Unique assays:         {len(assay_ids)}")

    status = get_chembl_status()
    print(f"ChEMBL DB version:     {status.get('chembl_db_version')}")
    print(f"ChEMBL release date:   {status.get('chembl_release_date')}")

    print("\nRetrieving assay records (cache is used unless --refresh-cache is set)...")
    assay_records, api_errors = fetch_assays(
        assay_ids=assay_ids,
        cache_path=cache_path,
        refresh_cache=args.refresh_cache,
        retries=args.retries,
        retry_delay=args.retry_delay,
    )

    audit = assay_counts.merge(
        assay_records,
        on="assay_chembl_id",
        how="left",
        validate="one_to_one",
    )
    audit = annotate_keywords(audit)

    previous_manual_path: Path | None = None
    if not args.reset_manual:
        if args.manual_input:
            previous_manual_path = Path(args.manual_input)
        elif output_path.exists():
            previous_manual_path = output_path
    audit = preserve_manual_annotations(audit, previous_manual_path)

    audit["assay_api_url"] = audit["assay_chembl_id"].map(
        lambda value: f"https://www.ebi.ac.uk/chembl/api/data/assay/{value}.json"
    )
    if "document_chembl_id" in audit.columns:
        audit["document_api_url"] = audit["document_chembl_id"].fillna("").map(
            lambda value: (
                f"https://www.ebi.ac.uk/chembl/api/data/document/{value}.json"
                if str(value).strip()
                else ""
            )
        )

    priority_order = {"high": 0, "candidate_strict": 1, "normal": 2}
    audit["_priority_order"] = audit["auto_review_priority"].map(priority_order).fillna(9)
    audit = audit.sort_values(
        ["_priority_order", "n_compounds", "n_measurements", "assay_chembl_id"],
        ascending=[True, False, False, True],
    ).drop(columns="_priority_order")

    audit = audit[ordered_columns(audit)]
    summary = build_summary(audit, raw)

    # utf-8-sig opens cleanly in Microsoft Excel on Windows.
    audit.to_csv(output_path, index=False, encoding="utf-8-sig")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    metadata = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(input_path.resolve()),
        "input_sha256": sha256_file(input_path),
        "output_file": str(output_path.resolve()),
        "summary_file": str(summary_path.resolve()),
        "cache_file": str(cache_path.resolve()),
        "chembl_status": status,
        "python_version": sys.version,
        "package_versions": {
            "pandas": package_version("pandas"),
            "requests": package_version("requests"),
            "chembl_webresource_client": package_version("chembl-webresource-client"),
        },
        "n_activity_rows": int(len(raw)),
        "n_unique_compounds": int(raw["inchikey"].nunique(dropna=True)),
        "n_assays_expected": int(len(assay_ids)),
        "n_assays_written": int(len(audit)),
        "n_api_errors": int(len(api_errors)),
        "api_errors": api_errors,
        "keyword_rules": KEYWORD_RULES,
        "methodological_warning": (
            "Automatic flags are keyword-based triage only. In particular, "
            "calcium plus AITC does not establish a homogeneous antagonist assay; "
            "application order, readout, expression system, construct, direct agonism, "
            "and desensitization risk require manual review."
        ),
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("\n" + "=" * 76)
    print("OUTPUT")
    print("=" * 76)
    print(f"Audit table:           {output_path.resolve()}")
    print(f"Flag summary:          {summary_path.resolve()}")
    print(f"Metadata:              {metadata_path.resolve()}")
    print(f"Assay cache:           {cache_path.resolve()}")
    print(f"API errors:            {len(api_errors)}")

    candidate = audit["candidate_strict_review"].fillna(False)
    candidate_assays = int(candidate.sum())
    candidate_compounds = unique_compound_coverage(
        raw, audit.loc[candidate, "assay_chembl_id"]
    )
    print("\nAUTOMATIC TRIAGE ONLY")
    print(f"Candidate strict-review assays:    {candidate_assays}")
    print(f"Unique compounds in their union:   {candidate_compounds}")
    print(
        "This is not yet the strict or homogeneous dataset. Complete the manual_* "
        "columns before inclusion/exclusion decisions."
    )

    if api_errors and not args.allow_api_errors:
        print(
            "\nERROR: Some assay records could not be retrieved. Files were saved for "
            "diagnosis, but the audit is incomplete. Rerun the script; successful "
            "records are cached. Use --allow-api-errors only intentionally.",
            file=sys.stderr,
        )
        return 2

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AuditError, FileNotFoundError, pd.errors.ParserError, requests.RequestException) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc