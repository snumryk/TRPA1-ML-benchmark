from __future__ import annotations

import csv
import hashlib
import json
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


INPUT_ASSAY_AUDIT = Path("trpa1_assay_audit.csv")
OUTPUT_ROOT = Path("trpa1_raw_snapshots")
EXPECTED_ASSAY_COUNT = 97
EXPECTED_CHEMBL_RELEASE = "ChEMBL_37"
TARGET_CHEMBL_ID = "CHEMBL6007"
TARGET_ORGANISM = "Homo sapiens"
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
PAGE_SIZE = 1000
TIMEOUT = (15, 120)
REQUEST_PAUSE = 0.10
SCRIPT_VERSION = "1.0.0"


def text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_json(path: Path, data: Any) -> None:
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def save_jsonl_record(handle: Any, record: dict[str, Any]) -> None:
    handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def read_assay_ids(path: Path) -> tuple[list[str], int]:
    if not path.exists():
        raise FileNotFoundError(f"Не знайдено {path.resolve()}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"assay_chembl_id", "target_chembl_id", "assay_organism"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Відсутні колонки: {sorted(missing)}")

        assay_ids: list[str] = []
        blank_organism = 0

        for line_number, row in enumerate(reader, start=2):
            assay_id = text(row["assay_chembl_id"])
            target_id = text(row["target_chembl_id"])
            organism = text(row["assay_organism"])

            if not assay_id:
                raise ValueError(f"Рядок {line_number}: порожній assay ID")
            if target_id != TARGET_CHEMBL_ID:
                raise ValueError(
                    f"Рядок {line_number}: {assay_id} має target {target_id!r}"
                )
            if organism and organism != TARGET_ORGANISM:
                raise ValueError(
                    f"Рядок {line_number}: {assay_id} має organism {organism!r}"
                )
            if not organism:
                blank_organism += 1

            assay_ids.append(assay_id)

    if len(set(assay_ids)) != len(assay_ids):
        raise ValueError("У assay audit є дублікати assay ID")
    if len(assay_ids) != EXPECTED_ASSAY_COUNT:
        raise ValueError(
            f"Знайдено {len(assay_ids)} assays, очікувалося {EXPECTED_ASSAY_COUNT}"
        )

    return sorted(assay_ids), blank_organism


def make_session() -> requests.Session:
    retries = Retry(
        total=5,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retries)
    session = requests.Session()
    session.mount("https://", adapter)
    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": f"TRPA1-raw-snapshot/{SCRIPT_VERSION}",
        }
    )
    return session


def get_json(
    session: requests.Session,
    url: str,
    params: dict[str, Any] | None = None,
) -> dict[str, Any]:
    response = session.get(url, params=params, timeout=TIMEOUT)
    response.raise_for_status()
    try:
        payload = response.json()
    except ValueError as error:
        raise RuntimeError(f"API повернув не JSON: {response.url}") from error
    if not isinstance(payload, dict):
        raise RuntimeError(f"API повернув не JSON-об'єкт: {response.url}")
    time.sleep(REQUEST_PAUSE)
    return payload


def fetch_status(session: requests.Session) -> dict[str, Any]:
    return get_json(session, f"{BASE_URL}/status.json")


def check_release(status: dict[str, Any]) -> str:
    release = text(status.get("chembl_db_version"))
    if release != EXPECTED_CHEMBL_RELEASE:
        raise RuntimeError(
            f"ChEMBL release = {release!r}, очікувався "
            f"{EXPECTED_CHEMBL_RELEASE!r}"
        )
    return release


def fetch_target(session: requests.Session) -> dict[str, Any]:
    target = get_json(session, f"{BASE_URL}/target/{TARGET_CHEMBL_ID}.json")
    if text(target.get("target_chembl_id")) != TARGET_CHEMBL_ID:
        raise RuntimeError("Target endpoint повернув інший target ID")
    if text(target.get("organism")) != TARGET_ORGANISM:
        raise RuntimeError("CHEMBL6007 не підтверджено як Homo sapiens target")
    return target


def fetch_assays(
    session: requests.Session,
    assay_ids: list[str],
    output_path: Path,
) -> int:
    blank_organism = 0
    temp_path = output_path.with_suffix(".jsonl.part")

    with temp_path.open("w", encoding="utf-8", newline="\n") as output:
        for index, assay_id in enumerate(assay_ids, start=1):
            record = get_json(session, f"{BASE_URL}/assay/{assay_id}.json")

            if text(record.get("assay_chembl_id")) != assay_id:
                raise RuntimeError(f"API повернув не той assay для {assay_id}")
            if text(record.get("target_chembl_id")) != TARGET_CHEMBL_ID:
                raise RuntimeError(f"Assay {assay_id} має інший target")

            organism = text(record.get("assay_organism"))
            if organism and organism != TARGET_ORGANISM:
                raise RuntimeError(f"Assay {assay_id} має organism {organism!r}")
            if not organism:
                blank_organism += 1

            save_jsonl_record(output, record)
            print(f"Assays: {index}/{len(assay_ids)} — {assay_id}")

    temp_path.replace(output_path)
    return blank_organism


def fetch_one_assay_activities(
    session: requests.Session,
    assay_id: str,
    output: Any,
    seen_activity_ids: set[str],
) -> tuple[int, int]:
    offset = 0
    expected_total: int | None = None
    written = 0
    blank_organism = 0

    while True:
        payload = get_json(
            session,
            f"{BASE_URL}/activity.json",
            params={
                "assay_chembl_id": assay_id,
                "limit": PAGE_SIZE,
                "offset": offset,
                "order_by": "activity_id",
            },
        )
        records = payload.get("activities")
        page_meta = payload.get("page_meta")
        if not isinstance(records, list) or not isinstance(page_meta, dict):
            raise RuntimeError(f"Некоректна activity page для {assay_id}")

        total = int(page_meta.get("total_count", len(records)))
        if expected_total is None:
            expected_total = total
        elif total != expected_total:
            raise RuntimeError(f"total_count змінився під час {assay_id}")

        if not records:
            break

        for record in records:
            activity_id = text(record.get("activity_id"))
            if text(record.get("assay_chembl_id")) != assay_id:
                raise RuntimeError(f"У {assay_id} потрапив record іншого assay")
            if text(record.get("target_chembl_id")) != TARGET_CHEMBL_ID:
                raise RuntimeError(f"Activity {activity_id} має інший target")

            organism = text(record.get("target_organism"))
            if organism and organism != TARGET_ORGANISM:
                raise RuntimeError(f"Activity {activity_id} має {organism!r}")
            if not organism:
                blank_organism += 1

            if not activity_id:
                raise RuntimeError(f"У {assay_id} є record без activity_id")
            if activity_id in seen_activity_ids:
                raise RuntimeError(f"Повторився activity_id={activity_id}")

            seen_activity_ids.add(activity_id)
            save_jsonl_record(output, record)
            written += 1

        offset += len(records)
        if expected_total is not None and written >= expected_total:
            break

    expected_total = expected_total or 0
    if written != expected_total:
        raise RuntimeError(
            f"{assay_id}: записано {written}, API повідомив {expected_total}"
        )
    return written, blank_organism


def fetch_all_activities(
    session: requests.Session,
    assay_ids: list[str],
    output_path: Path,
) -> tuple[dict[str, int], int, int]:
    counts: dict[str, int] = {}
    seen_activity_ids: set[str] = set()
    blank_organism_total = 0
    temp_path = output_path.with_suffix(".jsonl.part")

    with temp_path.open("w", encoding="utf-8", newline="\n") as output:
        for index, assay_id in enumerate(assay_ids, start=1):
            count, blank_organism = fetch_one_assay_activities(
                session, assay_id, output, seen_activity_ids
            )
            counts[assay_id] = count
            blank_organism_total += blank_organism
            print(
                f"Activities: {index}/{len(assay_ids)} — "
                f"{assay_id}: {count} records"
            )

    temp_path.replace(output_path)
    return counts, len(seen_activity_ids), blank_organism_total


def main() -> int:
    output_dir: Path | None = None
    try:
        assay_ids, audit_blank_organism = read_assay_ids(INPUT_ASSAY_AUDIT)
        session = make_session()

        start_status = fetch_status(session)
        release = check_release(start_status)

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        output_dir = OUTPUT_ROOT / f"{release}_{stamp}"
        output_dir.mkdir(parents=True, exist_ok=False)

        copied_audit = output_dir / INPUT_ASSAY_AUDIT.name
        shutil.copy2(INPUT_ASSAY_AUDIT, copied_audit)
        (output_dir / "assay_ids.txt").write_text(
            "\n".join(assay_ids) + "\n", encoding="utf-8"
        )
        save_json(output_dir / "chembl_status_start.json", start_status)
        save_json(output_dir / "target_CHEMBL6007.json", fetch_target(session))

        assays_path = output_dir / "assays_raw.jsonl"
        live_assay_blank_organism = fetch_assays(
            session, assay_ids, assays_path
        )

        activities_path = output_dir / "activities_raw.jsonl"
        counts, total_activities, activity_blank_organism = fetch_all_activities(
            session, assay_ids, activities_path
        )

        end_status = fetch_status(session)
        if check_release(end_status) != release:
            raise RuntimeError("ChEMBL release змінився під час запуску")
        save_json(output_dir / "chembl_status_end.json", end_status)

        manifest = {
            "script_version": SCRIPT_VERSION,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "chembl_release": release,
            "target_chembl_id": TARGET_CHEMBL_ID,
            "target_organism": TARGET_ORGANISM,
            "n_assays": len(assay_ids),
            "n_activity_records": total_activities,
            "activity_counts_by_assay": counts,
            "audit_blank_assay_organism": audit_blank_organism,
            "live_blank_assay_organism": live_assay_blank_organism,
            "activity_blank_target_organism": activity_blank_organism,
            "input_audit_sha256": sha256(INPUT_ASSAY_AUDIT),
            "query": {
                "filter": "assay_chembl_id only",
                "page_size": PAGE_SIZE,
                "order_by": "activity_id",
                "cache_used": False,
                "no_filter_on_type_relation_or_pchembl": True,
            },
            "file_sha256": {
                "assays_raw.jsonl": sha256(assays_path),
                "activities_raw.jsonl": sha256(activities_path),
                INPUT_ASSAY_AUDIT.name: sha256(copied_audit),
            },
        }
        save_json(output_dir / "snapshot_manifest.json", manifest)
        (output_dir / "SNAPSHOT_COMPLETE.txt").write_text(
            "Snapshot completed successfully.\n", encoding="utf-8"
        )

        print("\nГОТОВО")
        print(f"Папка: {output_dir.resolve()}")
        print(f"Assays: {len(assay_ids)}")
        print(f"Activity records: {total_activities}")
        print("Labels та агрегація не створювалися.")
        return 0

    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr)
        if output_dir is not None and output_dir.exists():
            (output_dir / "SNAPSHOT_FAILED.txt").write_text(
                error, encoding="utf-8"
            )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
