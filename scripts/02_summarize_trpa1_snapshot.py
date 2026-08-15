from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


SNAPSHOT_ROOT = Path("trpa1_raw_snapshots")
TARGET_ID = "CHEMBL6007"
TARGET_ORGANISM = "Homo sapiens"
EXPECTED_ASSAYS = 97
SCRIPT_VERSION = "2.1.0"


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    if isinstance(value, str):
        return not value.strip()
    if isinstance(value, (list, dict, tuple, set)):
        return len(value) == 0
    return False


def text(value: Any) -> str:
    return "" if is_missing(value) else str(value).strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path.name}: очікувався JSON-об'єкт")
    return data


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"{path.name}, рядок {line_number}: некоректний JSON"
                ) from error
            if not isinstance(record, dict):
                raise ValueError(
                    f"{path.name}, рядок {line_number}: очікувався JSON-об'єкт"
                )
            records.append(record)
    return records


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def first_present(frame: pd.DataFrame, primary: str, fallback: str) -> pd.Series:
    first = frame[primary].map(text) if primary in frame else pd.Series("", index=frame.index)
    second = frame[fallback].map(text) if fallback in frame else pd.Series("", index=frame.index)
    return first.mask(first.eq(""), second).astype("string")


def parse_bool(value: Any) -> bool | None:
    if is_missing(value):
        return None
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "1.0", "yes", "y"}:
        return True
    if normalized in {"false", "0", "0.0", "no", "n"}:
        return False
    return None


def find_snapshot(explicit: Path | None) -> Path:
    if explicit is not None:
        snapshot = explicit.resolve()
        if not snapshot.is_dir():
            raise FileNotFoundError(snapshot)
        return snapshot

    candidates = []
    if SNAPSHOT_ROOT.exists():
        candidates = sorted(
            path
            for path in SNAPSHOT_ROOT.iterdir()
            if path.is_dir() and (path / "SNAPSHOT_COMPLETE.txt").exists()
        )
    if not candidates:
        raise FileNotFoundError("Немає завершеного snapshot. Спочатку запусти Скрипт 1.")
    return candidates[-1].resolve()


def validate_snapshot(snapshot: Path) -> dict[str, Any]:
    required = [
        "SNAPSHOT_COMPLETE.txt",
        "snapshot_manifest.json",
        "activities_raw.jsonl",
        "assays_raw.jsonl",
    ]
    absent = [name for name in required if not (snapshot / name).exists()]
    if absent:
        raise FileNotFoundError(f"У snapshot немає файлів: {absent}")
    if (snapshot / "SNAPSHOT_FAILED.txt").exists():
        raise RuntimeError("Snapshot має markers COMPLETE і FAILED одночасно")

    manifest = read_json(snapshot / "snapshot_manifest.json")
    if text(manifest.get("target_chembl_id")) != TARGET_ID:
        raise ValueError("Manifest містить не CHEMBL6007")
    if int(manifest.get("n_assays", -1)) != EXPECTED_ASSAYS:
        raise ValueError("Manifest не містить рівно 97 assays")

    hashes = manifest.get("file_sha256", {})
    if not isinstance(hashes, dict):
        raise ValueError("Manifest не містить file_sha256")
    for name, expected in hashes.items():
        path = snapshot / str(name)
        if not path.exists() or sha256(path) != str(expected):
            raise ValueError(f"Файл відсутній або пошкоджений: {name}")
    return manifest


def build_records(raw: list[dict[str, Any]]) -> pd.DataFrame:
    frame = pd.DataFrame(raw)
    required = {"activity_id", "assay_chembl_id", "target_chembl_id"}
    if frame.empty or not required.issubset(frame.columns):
        raise ValueError("activities_raw.jsonl порожній або має неправильну схему")

    for column in required:
        frame[column] = frame[column].map(text).astype("string")
    if frame["activity_id"].eq("").any() or frame["activity_id"].duplicated().any():
        raise ValueError("Є порожні або повторні activity_id")
    if not frame["target_chembl_id"].eq(TARGET_ID).all():
        raise ValueError("Є records не для CHEMBL6007")

    organism = frame.get("target_organism", pd.Series("", index=frame.index)).map(text)
    if (organism.ne("") & organism.ne(TARGET_ORGANISM)).any():
        raise ValueError("Є records не для Homo sapiens")

    frame["effective_type"] = first_present(frame, "standard_type", "type").str.upper()
    frame["effective_relation"] = first_present(
        frame, "standard_relation", "relation"
    )
    for column in ["standard_value", "standard_upper_value", "pchembl_value"]:
        if column not in frame:
            frame[column] = pd.NA
        frame[f"_{column}"] = pd.to_numeric(frame[column], errors="coerce")

    comments = pd.Series("", index=frame.index, dtype="string")
    for column in ["activity_comment", "standard_text_value", "text_value"]:
        if column in frame:
            comments += " " + frame[column].map(text).astype("string")

    frame["has_standard_text_value"] = frame.get(
        "standard_text_value",
        pd.Series("", index=frame.index),
    ).map(text).ne("")

    frame["exact_ic50"] = (
        frame["effective_type"].eq("IC50")
        & frame["effective_relation"].eq("=")
        & frame["_standard_value"].notna()
        & frame["_standard_upper_value"].isna()
    )
    frame["non_exact_ic50"] = (
    frame["effective_type"].eq("IC50")
    & ~frame["exact_ic50"]
    )
    frame["interval_ic50"] = (
        frame["effective_type"].eq("IC50")
        & frame["_standard_upper_value"].notna()
    )
    frame["right_censored_ic50"] = (
        frame["effective_type"].eq("IC50")
        & frame["effective_relation"].isin([">", ">=", ">>"])
    )
    frame["left_censored_ic50"] = (
        frame["effective_type"].eq("IC50")
        & frame["effective_relation"].isin(["<", "<=", "<<"])
    )
    frame["text_contains_inactive"] = comments.str.contains(
        r"\binactive\b|\bnot\s+active\b", case=False, regex=True, na=False
    )
    frame["single_point_inhibition"] = (
        frame["effective_type"].str.contains("INHIB", case=False, na=False)
        & ~frame["effective_type"].eq("IC50")
    )
    frame["missing_structure"] = frame.get(
        "canonical_smiles", pd.Series("", index=frame.index)
    ).map(text).eq("")
    frame["potential_duplicate_true"] = frame.get(
        "potential_duplicate", pd.Series(None, index=frame.index)
    ).map(parse_bool).eq(True)
    frame["has_validity_comment"] = frame.get(
        "data_validity_comment", pd.Series("", index=frame.index)
    ).map(text).ne("")
    frame["missing_target_organism"] = organism.eq("")

    if "activity_properties" not in frame:
        frame["activity_properties"] = None
    frame["activity_properties_json"] = frame["activity_properties"].map(
        lambda value: ""
        if is_missing(value)
        else json.dumps(value, ensure_ascii=False, sort_keys=True)
    )

    flags = [
        "exact_ic50",
        "non_exact_ic50",
        "interval_ic50",
        "right_censored_ic50",
        "left_censored_ic50",
        "text_contains_inactive",
        "single_point_inhibition",
        "missing_structure",
        "potential_duplicate_true",
        "has_validity_comment",
        "missing_target_organism",
        "has_standard_text_value",
    ]
    frame["descriptive_flags"] = frame.apply(
        lambda row: ";".join(flag for flag in flags if bool(row[flag])), axis=1
    )

    columns = [
        "activity_id", "assay_chembl_id", "document_chembl_id",
        "molecule_chembl_id", "canonical_smiles", "target_chembl_id",
        "target_organism", "standard_type", "type", "effective_type",
        "standard_relation", "relation", "effective_relation", "standard_value",
        "standard_upper_value", "standard_units", "pchembl_value", "value",
        "upper_value", "units", "activity_comment", "standard_text_value",
        "text_value", "data_validity_comment", "potential_duplicate",
        "activity_properties_json", *flags, "descriptive_flags",
    ]
    for column in columns:
        if column not in frame:
            frame[column] = pd.NA
    return frame[columns].copy()


def build_assay_overview(records: pd.DataFrame, raw_assays: list[dict[str, Any]]) -> pd.DataFrame:
    assays = pd.DataFrame(raw_assays)
    if "assay_chembl_id" not in assays or len(assays) != EXPECTED_ASSAYS:
        raise ValueError("assays_raw.jsonl не містить рівно 97 assays")
    assays["assay_chembl_id"] = assays["assay_chembl_id"].map(text)
    if assays["assay_chembl_id"].duplicated().any():
        raise ValueError("У assay metadata є дублікати")

    summary = records.groupby("assay_chembl_id").agg(
        n_records=("activity_id", "size"),
        n_unique_molecules=("molecule_chembl_id", "nunique"),
        n_exact_ic50=("exact_ic50", "sum"),
        n_interval_ic50=("interval_ic50", "sum"),
        n_right_censored_ic50=("right_censored_ic50", "sum"),
        n_left_censored_ic50=("left_censored_ic50", "sum"),
        n_text_inactive=("text_contains_inactive", "sum"),
        n_single_point_inhibition=("single_point_inhibition", "sum"),
        n_missing_structure=("missing_structure", "sum"),
        n_potential_duplicate=("potential_duplicate_true", "sum"),
    ).reset_index()

    keep = [
        "assay_chembl_id", "description", "assay_type", "assay_organism",
        "document_chembl_id", "confidence_score", "bao_format", "cell_chembl_id",
    ]
    for column in keep:
        if column not in assays:
            assays[column] = pd.NA
    result = assays[keep].merge(summary, on="assay_chembl_id", how="left")
    count_columns = [column for column in result if column.startswith("n_")]
    result[count_columns] = result[count_columns].fillna(0).astype(int)
    return result.sort_values("assay_chembl_id")


def build_counts(records: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, column in [
        ("activity_type", "effective_type"),
        ("relation", "effective_relation"),
        ("standard_text_value", "standard_text_value"),
    ]:
        counts = records[column].map(text).replace("", "<MISSING>").value_counts()
        rows.extend(
            {"group": group, "value": value, "n_records": int(count)}
            for value, count in counts.items()
        )
    for flag in [
        "exact_ic50", "interval_ic50", "right_censored_ic50",
        "left_censored_ic50", "text_contains_inactive",
        "single_point_inhibition", "missing_structure",
        "potential_duplicate_true", "has_validity_comment",
        "missing_target_organism",
    ]:
        rows.append(
            {"group": "descriptive_flag", "value": flag, "n_records": int(records[flag].sum())}
        )
    return pd.DataFrame(rows)


def build_field_coverage(raw: list[dict[str, Any]]) -> pd.DataFrame:
    total = len(raw)
    rows = []
    for field in sorted({key for record in raw for key in record}):
        present = sum(not is_missing(record.get(field)) for record in raw)
        rows.append(
            {"field": field, "n_present": present, "n_missing": total - present,
             "coverage_percent": round(present / total * 100, 3)}
        )
    return pd.DataFrame(rows).sort_values(
        ["coverage_percent", "field"], ascending=[False, True]
    )

def build_activity_properties_coverage(
    raw: list[dict[str, Any]],
) -> pd.DataFrame:
    property_rows: list[dict[str, Any]] = []

    for record in raw:
        properties = record.get("activity_properties")

        if is_missing(properties):
            continue

        if isinstance(properties, dict):
            properties = [properties]

        if not isinstance(properties, list):
            continue

        for property_item in properties:
            if isinstance(property_item, dict):
                property_rows.append(property_item)

    if not property_rows:
        return pd.DataFrame(
            columns=["field", "n_present", "n_property_rows", "coverage_percent"]
        )

    fields = sorted({
        key
        for property_item in property_rows
        for key in property_item
    })

    rows = []
    total = len(property_rows)

    for field in fields:
        present = sum(
            not is_missing(property_item.get(field))
            for property_item in property_rows
        )
        rows.append({
            "field": field,
            "n_present": present,
            "n_property_rows": total,
            "coverage_percent": round(present / total * 100, 3),
        })

    return pd.DataFrame(rows).sort_values(
        ["coverage_percent", "field"],
        ascending=[False, True],
    )

def main() -> int:
    output_dir: Path | None = None
    try:
        parser = argparse.ArgumentParser(
            description="Огляд raw TRPA1 snapshot без labels та агрегації."
        )
        parser.add_argument("snapshot_dir", nargs="?", type=Path)
        args = parser.parse_args()

        snapshot = find_snapshot(args.snapshot_dir)
        manifest = validate_snapshot(snapshot)
        raw_records = read_jsonl(snapshot / "activities_raw.jsonl")
        raw_assays = read_jsonl(snapshot / "assays_raw.jsonl")
        records = build_records(raw_records)
        assays = build_assay_overview(records, raw_assays)

        if len(records) != int(manifest.get("n_activity_records", -1)):
            raise ValueError("Кількість records не збігається з manifest")
        grouped_counts = records.groupby("assay_chembl_id").size().astype(int).to_dict()
        expected = {str(k): int(v) for k, v in manifest["activity_counts_by_assay"].items()}
        actual = {assay_id: int(grouped_counts.get(assay_id, 0)) for assay_id in expected}
        extra_assays = sorted(set(grouped_counts) - set(expected))
        if actual != expected or extra_assays:
            raise ValueError("Per-assay counts не збігаються з manifest")

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        output_dir = snapshot / f"step2_overview_{stamp}"
        output_dir.mkdir(exist_ok=False)

        write_csv(records, output_dir / "records_overview.csv")
        write_csv(assays, output_dir / "assay_overview.csv")
        write_csv(build_counts(records), output_dir / "counts_overview.csv")
        write_csv(build_field_coverage(raw_records), output_dir / "field_coverage.csv")
        write_csv(build_activity_properties_coverage(raw_records), output_dir / "activity_properties_field_coverage.csv",)

        review_mask = (
            records["interval_ic50"] | records["right_censored_ic50"]
            | records["non_exact_ic50"]
            | records["left_censored_ic50"] | records["text_contains_inactive"]
            | records["single_point_inhibition"]
            | records["has_standard_text_value"]
        )
        write_csv(records.loc[review_mask], output_dir / "non_exact_and_textual_records.csv")

        molecules = records["molecule_chembl_id"].map(text)
        summary = {
            "script_version": SCRIPT_VERSION,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_snapshot": str(snapshot),
            "chembl_release": text(manifest.get("chembl_release")),
            "n_activity_records": len(records),
            "n_assays": int(records["assay_chembl_id"].nunique()),
            "n_unique_molecules": int(molecules[molecules.ne("")].nunique()),
            "n_exact_ic50": int(records["exact_ic50"].sum()),
            "n_non_exact_ic50": int(records["non_exact_ic50"].sum()),
            "n_interval_ic50": int(records["interval_ic50"].sum()),
            "n_right_censored_ic50": int(records["right_censored_ic50"].sum()),
            "n_left_censored_ic50": int(records["left_censored_ic50"].sum()),
            "n_text_contains_inactive": int(records["text_contains_inactive"].sum()),
            "n_single_point_inhibition": int(records["single_point_inhibition"].sum()),
            "n_missing_structure": int(records["missing_structure"].sum()),
            "n_potential_duplicate": int(records["potential_duplicate_true"].sum()),
            "n_validity_comment": int(records["has_validity_comment"].sum()),
            "n_non_exact_or_textual": int(review_mask.sum()),
            "labels_created": False,
            "aggregation_created": False,
        }
        write_json(output_dir / "step2_summary.json", summary)
        (output_dir / "step2_report.txt").write_text(
            "\n".join(f"{key}: {value}" for key, value in summary.items()) + "\n",
            encoding="utf-8",
        )
        (output_dir / "STEP2_COMPLETE.txt").write_text(
            "Step 2 completed. No labels were created.\n", encoding="utf-8"
        )

        print("\nГОТОВО")
        print(f"Snapshot: {snapshot}")
        print(f"Результати: {output_dir.resolve()}")
        print(f"Records: {len(records)}")
        print(f"Exact IC50: {summary['n_exact_ic50']}")
        print(f"Non-exact/textual: {summary['n_non_exact_or_textual']}")
        print("Labels та агрегація не створювалися.")
        return 0

    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr)
        if output_dir is not None and output_dir.exists():
            (output_dir / "STEP2_FAILED.txt").write_text(error, encoding="utf-8")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
