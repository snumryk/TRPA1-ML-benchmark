#!/usr/bin/env python3
"""Create exactly 50,000 UNLABELED ChEMBL 37 molecules AFTER P/A exclusions.

Windows PowerShell (from experiments/02_classification/code):
    python build_u_v1.py --self-test
    python build_u_v1.py --download

If chembl_37_chemreps.txt.gz was downloaded already:
    python build_u_v1.py --chemreps "C:\\path\\chembl_37_chemreps.txt.gz"

Requires Python >=3.10 and rdkit (already installed in the trpa1 environment).
No pandas, no database, no ML training. P v1 and all source files are READ-ONLY.

Researcher-approved scope:
  * Universe: the ChEMBL 37 chemical-representations file (ALL compounds).
  * Historic molecular identity: RDKit parse -> SaltRemover -> largest
    remaining fragment -> canonical isomeric SMILES and full InChIKey.
  * Deduplicate on that full historical InChIKey.
  * Exclude the 2383 approved P v1 keys, 9 explicitly declined BindingDB
    candidates, and the separate researcher-selected agonist A v1 set.
  * EXCLUDE first, THEN uniformly sample exactly 50,000 unique structures,
    fixed seed 20261108 (without reducing the requested sample size).
  * Preserve the original ChEMBL SMILES and ChEMBL ID; U is NOT label=0.
  * NO further scientific/physicochemical filters, no re-audit of P.

Scientific limits:
  - This pool is not experimentally validated as inactive.
  - It may contain unrecognized human-TRPA1 inhibitors and known TRPA1
    compounds outside P/A. Excluding all other TRPA1-tested records was only
    a proposal, NOT approved and NOT applied.
  - A v1 is an evidence-based subset, not an exhaustive list of all agonists.
  - The planned 30k/8k/12k U_train/U_val/U_eval allocation is NOT frozen here:
    the eventual split must be coordinated with the P scaffold/source split.
  - Sampling is uniform over unique surviving standardized structures in
    ChEMBL 37, NOT over all imaginable chemical space and NOT cluster-balanced.

On success, writes a NEW directory and ZIP in:
  experiments/02_classification/data/unlabeled_U_v1_50000_excluding_agonists/
    trpa1_unlabeled_U_v1.csv
    manifest.json
No existing output is overwritten and no Git operation is performed.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import random
import sys
import tempfile
import urllib.request
import zipfile
from datetime import datetime, timezone

from rdkit import Chem, rdBase
from rdkit.Chem import SaltRemover

CHEMBL_RELEASE = 37
SOURCE_URL = ("https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/"
              "releases/chembl_37/chembl_37_chemreps.txt.gz")
P_REL = Path("experiments/02_classification/data/positive_P_v1_2383")
U_REL = Path("experiments/02_classification/data/unlabeled_U_v1_50000_excluding_agonists")
A_REL = Path("experiments/02_classification/data/agonist_A_v1")
A_NAME = "trpa1_known_agonists_v1.csv"
SOURCE_REL = Path("experiments/02_classification/data/chembl_37_chemreps.txt.gz")
APPROVED_P_SHA = "be9974cb3b49694d823b8ea0a261e610974dcdc257701fcf48289befa2eafb88"
TARGET_SIZE = 50_000
SEED = 20261108
EXPECTED_DECLINED = 9
REMOVER = SaltRemover.SaltRemover()
OUTPUT_COLUMNS = ["inchikey", "smiles", "chembl_id", "original_smiles", "source_inchi_key"]


def digest_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def digest_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def key_digest(keys) -> str:
    return digest_bytes(("\n".join(sorted(set(keys))) + "\n").encode("utf-8"))


def require(ok: bool, message: str):
    if not ok:
        raise ValueError(message)


def repo_root() -> Path:
    for start in (Path(__file__).resolve().parent, Path.cwd()):
        for candidate in (start, *start.parents):
            if (candidate / P_REL / "trpa1_positive_v1.csv").is_file():
                return candidate
    raise FileNotFoundError("Cannot find P v1. Use --root C:\\trpa1_data")


def load_approved_p(folder: Path) -> tuple[set[str], set[str], dict[str, str]]:
    p_csv = folder / "trpa1_positive_v1.csv"
    p_manifest = folder / "manifest.json"
    raw = p_csv.read_bytes()
    manifest_bytes = p_manifest.read_bytes()
    with io.StringIO(raw.decode("utf-8-sig"), newline="") as f:
        reader = csv.DictReader(f)
        require(reader.fieldnames is not None and {"inchikey", "label", "smiles"} <= set(reader.fieldnames),
                "P CSV does not have expected columns.")
        rows = list(reader)
    keys = [row["inchikey"].strip() for row in rows]
    require(len(rows) == 2383 and len(set(keys)) == 2383 and all(keys),
            "P v1 must contain exactly 2383 unique structures.")
    require(set(row["label"].strip() for row in rows) == {"1"},
            "P v1 labels differ from the approved positive set.")
    require(key_digest(keys) == APPROVED_P_SHA,
            "P v1 membership differs from the RESEARCHER-APPROVED key set. Stop.")
    manifest = json.loads(manifest_bytes)
    require(manifest.get("positive_key_set_sha256") == APPROVED_P_SHA,
            "P v1 manifest identity does not match the approved set.")
    require(manifest.get("counts", {}).get("positive_molecules") == 2383,
            "P v1 manifest molecule count is not 2383.")
    declined = manifest.get("declined_BindingDB_candidates")
    require(isinstance(declined, dict) and len(declined) == EXPECTED_DECLINED,
            "P manifest must list the nine researcher-declined BindingDB candidates.")
    exclusions = set(declined.keys())
    require(not exclusions & set(keys), "Excluded candidate appears in the approved P set.")
    return set(keys), exclusions, {"p_csv_sha256": digest_bytes(raw),
                                    "p_manifest_sha256": digest_bytes(manifest_bytes)}


def load_agonists(folder: Path) -> tuple[set[str], dict[str, str]]:
    """Read the researcher-selected A v1. No reclassification or new criteria."""
    csv_path = folder / A_NAME
    manifest_path = folder / "manifest.json"
    raw, meta_raw = csv_path.read_bytes(), manifest_path.read_bytes()
    with io.StringIO(raw.decode("utf-8-sig"), newline="") as buf:
        reader = csv.DictReader(buf)
        require(reader.fieldnames is not None and "inchikey" in reader.fieldnames,
                "Agonist CSV missing inchikey column.")
        rows = list(reader)
    keys = [r["inchikey"].strip() for r in rows]
    require(bool(keys) and all(keys) and len(keys) == len(set(keys)),
            "Agonist A v1 must have nonempty unique full InChIKeys.")
    meta = json.loads(meta_raw)
    require(meta.get("compounds") == len(keys), "Agonist A v1 count disagrees with manifest.")
    require(meta.get("output_file_sha256", {}).get(A_NAME) == digest_bytes(raw),
            "Agonist A v1 file hash disagrees with manifest.")
    return set(keys), {"a_csv_sha256": digest_bytes(raw),
                       "a_manifest_sha256": digest_bytes(meta_raw),
                       "a_keys_sha256": key_digest(keys)}


def historical_standardize(smiles: str) -> tuple[str, str]:
    """Exact historic TRPA1 processing, no added structural exclusion criteria."""
    if not smiles:
        return "", ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "", ""
        mol = REMOVER.StripMol(mol, dontRemoveEverything=True)
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if not fragments:
            return "", ""
        if len(fragments) > 1:
            mol = max(fragments, key=lambda f: f.GetNumHeavyAtoms())
        return (Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
                Chem.MolToInchiKey(mol))
    except (ValueError, RuntimeError, IndexError):
        return "", ""


def source_rows(path: Path):
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = set(reader.fieldnames or [])
        require({"chembl_id", "canonical_smiles", "standard_inchi_key"} <= fields,
                "Not a ChEMBL chemreps TSV (required: chembl_id, canonical_smiles, standard_inchi_key).")
        yield from reader


def select_pool(rows, p_keys: set[str], declined_keys: set[str], agonist_keys: set[str],
                size: int, seed: int, progress: bool = False):
    rng = random.Random(seed)
    seen = set()
    reservoir = []
    stats = dict(rows_read=0, missing_smiles=0, unparseable_smiles=0,
                 duplicate_standardized_keys=0, excluded_p=0, excluded_declined=0,
                 excluded_agonists=0, eligible_unique=0)
    for r in rows:
        stats["rows_read"] += 1
        original = (r.get("canonical_smiles") or "").strip()
        if not original:
            stats["missing_smiles"] += 1
            continue
        standard, key = historical_standardize(original)
        if not key:
            stats["unparseable_smiles"] += 1
            continue
        if key in seen:
            stats["duplicate_standardized_keys"] += 1
            continue
        seen.add(key)
        if key in p_keys:
            stats["excluded_p"] += 1
            continue
        if key in declined_keys:
            stats["excluded_declined"] += 1
            continue
        if key in agonist_keys:
            stats["excluded_agonists"] += 1
            continue
        stats["eligible_unique"] += 1
        candidate = {"inchikey": key, "smiles": standard,
                     "chembl_id": (r.get("chembl_id") or "").strip(),
                     "original_smiles": original,
                     "source_inchi_key": (r.get("standard_inchi_key") or "").strip()}
        if len(reservoir) < size:
            reservoir.append(candidate)
        else:
            j = rng.randrange(stats["eligible_unique"])
            if j < size:
                reservoir[j] = candidate
        if progress and stats["rows_read"] % 100000 == 0:
            print(f"  processed {stats['rows_read']:,} source rows; "
                  f"{stats['eligible_unique']:,} eligible unique structures", flush=True)
    require(stats["eligible_unique"] >= size,
            f"Only {stats['eligible_unique']} eligible unique structures; need {size}.")
    reservoir.sort(key=lambda r: r["inchikey"])
    require(len({r["inchikey"] for r in reservoir}) == size,
            "Internal error: selected pool has duplicate structures.")
    require(not {r["inchikey"] for r in reservoir} & (p_keys | declined_keys | agonist_keys),
            "Internal error: excluded identity appeared in U.")
    return reservoir, stats


def csv_bytes(rows: list[dict]) -> bytes:
    with io.StringIO(newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        return f.getvalue().encode("utf-8-sig")


def download_chemreps(target: Path):
    """Download one pinned release file if absent. Never replace existing input."""
    require(not target.exists(), f"Source already exists: {target}; refusing overwrite.")
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".partial")
    require(not part.exists(), f"Incomplete download exists: {part}; review/delete manually.")
    print(f"Downloading ChEMBL {CHEMBL_RELEASE} chemreps (approx. 293 MB)...\n{SOURCE_URL}", flush=True)
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "TRPA1-ML-benchmark/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=90) as response, part.open("xb") as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
        # Read the gzip header, but do not yet decompress the full file.
        with gzip.open(part, "rb") as check:
            require(bool(check.read(128)), "Downloaded source is empty/corrupt gzip.")
        part.replace(target)
    except Exception:
        if part.exists():
            part.unlink()
        raise
    print(f"Source cached: {target}", flush=True)


def self_test():
    a = historical_standardize("CCO")
    b = historical_standardize("CCO.[Na+]")
    assert a == b, "Historical salt handling changed."
    left = historical_standardize("C[C@H](O)C(=O)O")
    right = historical_standardize("C[C@@H](O)C(=O)O")
    assert left[1] and right[1] and left[1] != right[1], "Stereo identities merged."
    assert historical_standardize("invalid_bad_smiles") == ("", "")
    test_rows = [
        {"chembl_id": "CHEMBL1", "canonical_smiles": "CCO", "standard_inchi_key": ""},
        {"chembl_id": "CHEMBL2", "canonical_smiles": "CCO.[Na+]", "standard_inchi_key": ""},
        {"chembl_id": "CHEMBL3", "canonical_smiles": "CCN", "standard_inchi_key": ""},
        {"chembl_id": "CHEMBL4", "canonical_smiles": "CC(=O)O", "standard_inchi_key": ""},
        {"chembl_id": "CHEMBL5", "canonical_smiles": "", "standard_inchi_key": ""},
    ]
    p = {historical_standardize("CCN")[1]}
    akeys = {historical_standardize("CCO")[1]}
    # With CCO excluded, just one eligible compound (acetic acid) remains.
    selected1, stats = select_pool(test_rows, p, set(), akeys, 1, SEED)
    selected2, _ = select_pool(test_rows, p, set(), akeys, 1, SEED)
    assert selected1 == selected2
    assert stats["excluded_p"] == 1 and stats["duplicate_standardized_keys"] == 1
    assert stats["excluded_agonists"] == 1
    assert stats["missing_smiles"] == 1 and len(selected1) == 1
    assert all(k not in (p | akeys) for k in [r["inchikey"] for r in selected1])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "test.tsv.gz"
        with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["chembl_id", "canonical_smiles", "standard_inchi_key"], delimiter="\t")
            w.writeheader()
            w.writerows(test_rows)
        assert select_pool(source_rows(path), p, set(), akeys, 1, SEED)[0] == selected1
    # Exclusions are processed BEFORE random sampling, so sample size stays exact.
    more = test_rows + [{"chembl_id": "CHEMBL6", "canonical_smiles": "CCC",
                         "standard_inchi_key": ""}]
    picked, st = select_pool(more, p, set(), akeys, 2, SEED)
    assert len(picked) == 2 and st["excluded_agonists"] == 1
    assert not set(x["inchikey"] for x in picked) & (p | akeys)
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        # Match an independently built agonist A v1 manifest, without running chemistry.
        import csv as _csv
        with (folder / A_NAME).open("w", encoding="utf-8-sig", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=["inchikey"])
            w.writeheader(); w.writerow({"inchikey": next(iter(akeys))})
        raw = (folder / A_NAME).read_bytes()
        (folder / "manifest.json").write_text(json.dumps({"compounds": 1,
            "output_file_sha256": {A_NAME: digest_bytes(raw)}}), encoding="utf-8")
        actual, _ = load_agonists(folder)
        assert actual == akeys
    print("Self-tests: PASS (historical RDKit, exclusions BEFORE sampling, exact size, reproducibility, gzip TSV, agonist manifest).")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, help="Repository root, e.g. C:\\trpa1_data")
    parser.add_argument("--chemreps", type=Path, help="Local chembl_37_chemreps.txt.gz or uncompressed TSV")
    parser.add_argument("--agonists", type=Path, help="Directory containing approved agonist A v1 CSV and manifest")
    parser.add_argument("--download", action="store_true", help="Download the pinned ChEMBL 37 chemreps if absent")
    parser.add_argument("--self-test", action="store_true", help="Small offline tests; no project files read/changed")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0

    root = args.root.resolve() if args.root else repo_root()
    p_folder = root / P_REL
    p_keys, declined, p_hashes = load_approved_p(p_folder)
    a_folder = args.agonists.resolve() if args.agonists else root / A_REL
    a_keys, a_hashes = load_agonists(a_folder)
    source = args.chemreps.resolve() if args.chemreps else root / SOURCE_REL
    out = root / U_REL
    zip_path = out.with_suffix(".zip")
    require(not out.exists() and not zip_path.exists(),
            "Output already exists. Refusing to overwrite U (original U is preserved).")
    if not source.exists():
        if args.download:
            download_chemreps(source)
        else:
            raise FileNotFoundError(
                f"Missing {source}. Run with --download or provide --chemreps PATH.")
    print(f"Source: {source}\nP: {len(p_keys)}; declined P: {len(declined)}; A: {len(a_keys)}\n"
          f"Sampling exactly {TARGET_SIZE:,} U AFTER exclusions, seed={SEED}.", flush=True)
    source_hash_before = digest_file(source)
    sample, stats = select_pool(source_rows(source), p_keys, declined, a_keys,
                                TARGET_SIZE, SEED, progress=True)
    require(source_hash_before == digest_file(source), "ChEMBL input changed during processing.")
    require(p_hashes == {"p_csv_sha256": digest_file(p_folder / 'trpa1_positive_v1.csv'),
                         "p_manifest_sha256": digest_file(p_folder / 'manifest.json')},
            "P input changed during processing.")
    require(a_hashes["a_csv_sha256"] == digest_file(a_folder / A_NAME) and
            a_hashes["a_manifest_sha256"] == digest_file(a_folder / 'manifest.json'),
            "A input changed during processing.")

    csv_content = csv_bytes(sample)
    manifest = {
        "version": "U-v1-ChEMBL37-50000-excluding-agonist-A-v1-2026-11-08",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": f"ChEMBL {CHEMBL_RELEASE} chemical representations",
        "source_url": SOURCE_URL,
        "source_path": str(source),
        "source_sha256": source_hash_before,
        "seed": SEED,
        "sampling": "uniform reservoir sampling over all eligible unique standardized full InChIKeys",
        "standardization": "historical RDKit parse, SaltRemover, largest fragment, isomeric canonical SMILES, InChIKey",
        "p_keys_excluded": len(p_keys), "declined_bindingdb_keys_excluded": len(declined),
        "agonist_a_keys_excluded": len(a_keys), "agonist_a_membership_sha256": a_hashes["a_keys_sha256"],
        "agonist_a_input_hashes": a_hashes,
        "p_membership_sha256": key_digest(p_keys),
        "p_input_hashes": p_hashes,
        "observed_source_counts": stats,
        "u_molecules": len(sample), "u_keys_sha256": key_digest(r["inchikey"] for r in sample),
        "u_csv_sha256": digest_bytes(csv_content),
        "python": sys.version,
        "rdkit": rdBase.rdkitVersion,
        "script_sha256": digest_file(Path(__file__)),
        "label": "U/unlabeled; no negative or zero class label assigned",
        "split_status": "NOT_ASSIGNED; coordinate any future U_train/U_val/U_eval split with P scaffold/source-aware protocol",
        "proposed_split_not_applied": {"U_train": 30000, "U_val": 8000, "U_eval": 12000},
        "known_limitation": [
            "Approved P v1, nine declined BindingDB candidates and agonist A v1 keys were excluded BEFORE sampling.",
            "Agonist A v1 is only a subset of documented agonists, not a comprehensive catalog.",
            "Other TRPA1-tested molecules were NOT excluded because that additional rule has not been approved.",
            "ChEMBL is an annotated bioactivity library, not a representative sample of all conceivable chemicals.",
            "U may contain experimentally unrecognized or other known TRPA1 inhibitors.",
            "No scaffold-aware train/validation/test partition was assigned.",
            "Original ChEMBL SMILES and identifiers are retained, and historical standardization was applied only for U molecular identity.",
        ],
    }
    out.mkdir(parents=True, exist_ok=False)
    try:
        (out / "trpa1_unlabeled_U_v1.csv").write_bytes(csv_content)
        (out / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        with zipfile.ZipFile(zip_path, "x", compression=zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()):
                z.write(f, f.name)
    except BaseException:
        (out / "INCOMPLETE.txt").write_text("Interrupted: do not use this partial dataset.\n", encoding="utf-8")
        raise
    print("\nU v1 created successfully")
    for key, val in stats.items():
        print(f"  {key}: {val}")
    print(f"  U_molecules: {len(sample)}\n  overlap_with_approved_P: 0\n"
          f"  overlap_with_9_declined: 0\n  overlap_with_agonists_A_v1: 0\n"
          f"  model_negative_labels_assigned: 0\n"
          f"  split_assigned: False\n  source_bytes_unchanged: True\n  P_bytes_unchanged: True")
    print(f"Output: {out}\nZIP: {zip_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, csv.Error, EOFError, zipfile.BadZipFile) as exc:
        print(f"ERROR: {exc}\nExisting P and source files were not modified.", file=sys.stderr)
        raise SystemExit(2)
