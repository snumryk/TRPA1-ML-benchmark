
#!/usr/bin/env python3
"""Propose an auditable P/U train/val/test split for TRPA1 experiment 02.

Research design (PROPOSAL, not a pre-approved frozen split):
  * Read the APPROVED P v1 (2,383) and U v1 (50,000) AS IS.
  * Use stored SMILES only to calculate non-chiral Bemis-Murcko scaffolds;
    never change SMILES, identity, activity labels or dataset membership.
  * Allocate each scaffold present anywhere in P or U to ONE partition.
    For acyclic molecules, whose Murcko scaffold is empty, use separate
    molecular-identity groups; explicitly report this limitation.
  * Audit shared ChEMBL/BindingDB evidence documents across P partitions.
    --source-mode strict additionally keeps every cited P document in one
    partition (which can make the split infeasible).
  * U targets 30,000/8,000/12,000 are HARD; P targets ~1,903/240/240 are
    approximate due to grouping. Do not break groups to force counts.
  * Default is a READ-ONLY PREVIEW. --write creates a NEW output folder/ZIP;
    it NEVER changes input CSVs/manifests or GitHub.

Usage (PowerShell from experiments/02_classification/code):
    python trpa1_split_scaffold_source_v1.py --self-test
    python trpa1_split_scaffold_source_v1.py
    python trpa1_split_scaffold_source_v1.py --write
    python trpa1_split_scaffold_source_v1.py --source-mode strict

Dependencies: Python 3.10+; rdkit already used in the project's trpa1 env.
No internet, pandas, sklearn, training, new filters or class-0 labels.

A split is not a result of an experiment. The researcher must review the
printed P deviations, source overlaps and acyclic-scaffold diagnostics before
using a generated split for model fitting or reporting it in a manuscript.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import random
import sys
import tempfile
import zipfile

from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold

BASE = Path("experiments/02_classification/data")
P_DIR = BASE / "positive_P_v1_2383"
U_DIR = BASE / "unlabeled_U_v1_50000_excluding_agonists"
OUT_DIR = BASE / "split_P2383_U50000_scaffold_source_proposal_v1"
P_APPROVED_SHA = "be9974cb3b49694d823b8ea0a261e610974dcdc257701fcf48289befa2eafb88"
U_APPROVED_SHA = "ff7fa29dcb24d9f6c470f247b2b6689a04be32635de08545389619e136739842"
SPLITS = ("train", "val", "test")
DEFAULT_U = (30000, 8000, 12000)
DEFAULT_P = (1903, 240, 240)  # proposal only; not approved as final counts
DEFAULT_SEED = 20261109


def require(condition: bool, explanation: str):
    if not condition:
        raise ValueError(explanation)


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def key_hash(keys) -> str:
    return sha256(("\n".join(sorted(set(keys))) + "\n").encode("utf-8"))


def read_csv(path: Path, required: set[str]):
    raw = path.read_bytes()
    with io.StringIO(raw.decode("utf-8-sig"), newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        require(fields and len(fields) == len(set(fields)), f"{path}: invalid/duplicate headers")
        require(required <= set(fields), f"{path}: missing columns {required-set(fields)}")
        rows = list(reader)
    require(rows and all(None not in r and all(v is not None for v in r.values()) for r in rows),
            f"{path}: malformed or empty CSV")
    return rows, fields, sha256(raw)


def find_root() -> Path:
    for start in (Path(__file__).resolve().parent, Path.cwd()):
        for candidate in (start, *start.parents):
            if (candidate / P_DIR / "trpa1_positive_v1.csv").is_file():
                return candidate
    raise FileNotFoundError("Cannot find approved P v1; pass --root C:\\trpa1_data")


def load_data(root: Path):
    p_path, u_path = root / P_DIR, root / U_DIR
    pm = json.loads((p_path / "manifest.json").read_text(encoding="utf-8-sig"))
    um = json.loads((u_path / "manifest.json").read_text(encoding="utf-8-sig"))
    p, pfields, phash = read_csv(p_path / "trpa1_positive_v1.csv", {"inchikey", "smiles", "label", "membership_origin"})
    u, ufields, uhash = read_csv(u_path / "trpa1_unlabeled_U_v1.csv", {"inchikey", "smiles", "chembl_id"})
    ev, _, evhash = read_csv(p_path / "measurements.csv", {"inchikey", "source", "document_id_original", "supporting_patent_from_previous_review"})

    pk, uk = [r["inchikey"].strip() for r in p], [r["inchikey"].strip() for r in u]
    require(len(p) == 2383 and len(set(pk)) == 2383 and all(pk), "P membership/count invalid")
    require(len(u) == 50000 and len(set(uk)) == 50000 and all(uk), "U membership/count invalid")
    require(not set(pk) & set(uk), "P and U identity overlap; no split generated")
    require(set(r["label"].strip() for r in p) == {"1"}, "P contains nonpositive label")
    require("label" not in ufields or all(r["label"].strip() not in {"0", "1"} for r in u),
            "U has unexpected binary labels")
    require(key_hash(pk) == P_APPROVED_SHA == pm.get("positive_key_set_sha256"),
            "P membership does not match frozen approved P v1")
    require(key_hash(uk) == U_APPROVED_SHA == um.get("u_keys_sha256"),
            "U membership does not match frozen U v1")
    require(phash == pm.get("output_sha256", {}).get("trpa1_positive_v1.csv"), "P CSV SHA256 differs from manifest")
    require(evhash == pm.get("output_sha256", {}).get("measurements.csv"), "P measurements SHA256 differs from manifest")
    require(uhash == um.get("u_csv_sha256"), "U CSV SHA256 differs from manifest")
    require(len(ev) == 2943 and pm.get("counts", {}).get("measurement_rows") == 2943,
            "P measurement count differs")
    require(all(r["inchikey"].strip() in set(pk) for r in ev), "Measurements reference unknown P keys")

    # One bibliographic document can support several molecules; count each
    # P molecule once per document, NOT once per database activity row.
    refs = defaultdict(set)
    for r in ev:
        doc = (r.get("supporting_patent_from_previous_review") or "").strip()
        if not doc:
            doc = (r.get("document_id_original") or "").strip()
        if doc:
            refs[r["inchikey"].strip()].add(r["source"].strip().upper() + ":" + doc.upper())
    summary = {
        "p_csv_sha256": phash, "u_csv_sha256": uhash,
        "p_measurements_sha256": evhash,
        "p_keys_sha256": key_hash(pk), "u_keys_sha256": key_hash(uk),
        "p_missing_document_refs": sum(not refs[k] for k in pk),
    }
    return p, pfields, u, ufields, refs, summary



def scaffold_key(smiles: str, inchikey: str) -> tuple[str, str, bool]:
    try:
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)

        require(
            mol is not None,
            f"RDKit cannot parse stored SMILES for {inchikey}"
        )

        # Work on a temporary molecule, never modify input data.
        # Murcko scaffolds are intentionally non-stereochemical.
        # Clearing stereo BEFORE scaffold extraction avoids the
        # RDKit Canon.cpp 'bad bond stereo' failure.
        scaffold_mol = Chem.Mol(mol)
        Chem.RemoveStereochemistry(scaffold_mol)

        scaffold = MurckoScaffold.MurckoScaffoldSmiles(
            mol=scaffold_mol,
            includeChirality=False
        )

    except (RuntimeError, ValueError) as exc:
        raise ValueError(
            f"Scaffold calculation failed for {inchikey}: {exc}"
        ) from exc

    if scaffold:
        return "S:" + scaffold, scaffold, False

    return "ACYCLIC:" + inchikey, "", True


def prepare_groups(p, u, refs, source_mode):
    groups = {}
    details = {"P": {}, "U": {}}
    for kind, rows in (("P", p), ("U", u)):
        for row in rows:
            key = row["inchikey"].strip()
            gid, scaffold, acyclic = scaffold_key(row["smiles"].strip(), key)
            item = groups.setdefault(gid, {"p": [], "u": [], "docs": set()})
            item[kind.lower()].append(key)
            if kind == "P":
                item["docs"].update(refs.get(key, ()))
            details[kind][key] = (gid, scaffold, acyclic)

    if source_mode == "strict":
        parent = {g: g for g in groups}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def join(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                if ra > rb:
                    ra, rb = rb, ra
                parent[rb] = ra

        doc_first = {}
        for key, (gid, _, _) in details["P"].items():
            for doc in refs.get(key, ()):
                if doc in doc_first:
                    join(gid, doc_first[doc])
                else:
                    doc_first[doc] = gid
        merged = {}
        for gid, g in groups.items():
            root = find(gid)
            dest = merged.setdefault(root, {"p": [], "u": [], "docs": set()})
            dest["p"].extend(g["p"])
            dest["u"].extend(g["u"])
            dest["docs"].update(g["docs"])
        groups = merged
        # Map original scaffold groups to their source-connected component.
        remap = {gid: find(gid) for gid in parent}
        details = {kind: {key: (remap[gid], scaffold, acyclic)
                          for key, (gid, scaffold, acyclic) in info.items()}
                   for kind, info in details.items()}
    return groups, details


def source_audit(p, assignment, refs):
    by_doc = defaultdict(lambda: [0, 0, 0])
    for r in p:
        key = r["inchikey"].strip()
        for doc in refs.get(key, ()):
            by_doc[doc][SPLITS.index(assignment[("P", key)])] += 1
    shared = {doc: counts for doc, counts in by_doc.items() if sum(n > 0 for n in counts) > 1}
    return by_doc, shared


def allocate(groups, ptarget, utarget, seed: int, trials: int, refs, p, source_mode):
    require(sum(ptarget) == sum(len(g["p"]) for g in groups.values()), "P target sum must equal P size")
    require(sum(utarget) == sum(len(g["u"]) for g in groups.values()), "U target sum must equal U size")
    p_groups = [(gid, g) for gid, g in groups.items() if g["p"]]
    u_groups = [(gid, g) for gid, g in groups.items() if not g["p"]]
    best = None
    best_quality = None

    for attempt in range(trials):
        rng = random.Random(seed + attempt * 7919)
        p_sorted = sorted(p_groups, key=lambda item: (-len(item[1]["p"]), rng.random()))
        u_sorted = sorted(u_groups, key=lambda item: (-len(item[1]["u"]), rng.random()))
        pc = [0, 0, 0]
        uc = [0, 0, 0]
        destinations = {}
        source_placed = defaultdict(set)
        for gid, g in p_sorted:
            pn, un = len(g["p"]), len(g["u"])
            order = list(range(3))
            rng.shuffle(order)
            def price(j):
                next_p = pc[:]
                next_p[j] += pn
                overflow = max(0, uc[j] + un - utarget[j])
                doc_cross = sum(bool(source_placed[d] - {j}) for d in g["docs"])
                return (sum(((next_p[k]-ptarget[k])**2)/max(1, ptarget[k]) for k in range(3))
                        + 1000 * (overflow**2) / max(1, utarget[j])
                        + (0.8 * doc_cross if source_mode == "audit" else 0))
            j = min(order, key=price)
            destinations[gid] = j
            pc[j] += pn
            uc[j] += un
            for d in g["docs"]:
                source_placed[d].add(j)

        # Fill remaining U capacity with whole U-only scaffold groups.
        # Acyclic molecules are singleton groups; exact totals are often
        # attainable but are never manufactured by breaking a ring scaffold.
        for gid, g in u_sorted:
            un = len(g["u"])
            choices = list(range(3))
            rng.shuffle(choices)
            fitting = [j for j in choices if uc[j] + un <= utarget[j]]
            if fitting:
                j = max(fitting, key=lambda k: (utarget[k]-uc[k], rng.random()))
            else:
                j = min(choices, key=lambda k: (uc[k]+un-utarget[k], rng.random()))
            destinations[gid] = j
            uc[j] += un

        assignment = {}
        for gid, group in groups.items():
            part = SPLITS[destinations[gid]]
            for key in group["p"]:
                assignment[("P", key)] = part
            for key in group["u"]:
                assignment[("U", key)] = part
        by_doc, shared_docs = source_audit(p, assignment, refs)
        deviation_u = sum(abs(x-y) for x, y in zip(uc, utarget))
        deviation_p = sum(abs(x-y) for x, y in zip(pc, ptarget))
        quality = (deviation_u, deviation_p, len(shared_docs))
        if best_quality is None or quality < best_quality:
            best_quality = quality
            best = (assignment, pc[:], uc[:], by_doc, shared_docs, attempt + 1)
        if quality == (0, 0, 0):
            break
    return best, best_quality


def csv_data(rows, fields):
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue().encode("utf-8-sig")


def save_output(root, out, p, pfields, u, ufields, info, refs, summary,
                target_p, target_u, source_mode, seed, trials, results):
    assignment, pc, uc, by_doc, shared_docs, winning_trial = results
    zip_path = out.with_name(out.name + ".zip")
    require(not out.exists() and not zip_path.exists(), f"Output already exists: {out} / {zip_path}; refusing overwrite")
    outputs = {}
    table = []
    for kind, rows in (("P", p), ("U", u)):
        for r in rows:
            key = r["inchikey"].strip()
            gid, scaffold, acyclic = info[kind][key]
            part = assignment[(kind, key)]
            table.append(dict(collection=kind, inchikey=key, split=part,
                              murcko_scaffold=scaffold, acyclic_no_scaffold=int(acyclic),
                              component_id=hashlib.sha256(gid.encode()).hexdigest()[:20]))
        for partition in SPLITS:
            subset = [dict(r, split=partition, murcko_scaffold=info[kind][r["inchikey"].strip()][1])
                      for r in rows if assignment[(kind, r["inchikey"].strip())] == partition]
            name = f"{kind}_{'eval' if (kind == 'U' and partition == 'test') else partition}.csv"
            outputs[name] = csv_data(subset, list((pfields if kind == "P" else ufields)) + ["split", "murcko_scaffold"])
    outputs["split_assignments.csv"] = csv_data(table, list(table[0]))
    source_rows = [dict(source_document=doc, P_train=cnt[0], P_val=cnt[1], P_test=cnt[2],
                        present_in_n_splits=sum(n > 0 for n in cnt))
                   for doc, cnt in sorted(by_doc.items())]
    outputs["p_source_distribution.csv"] = csv_data(source_rows, list(source_rows[0])) if source_rows else b""
    manifest = {
        "status": "SPLIT_PROPOSAL_NOT_ML_TRAINED",
        "protocol": "Nonchiral Bemis-Murcko scaffold-group assignment across P+U; acyclic singleton by original InChIKey",
        "source_mode": source_mode,
        "source_note": "Source information only available for P evidence; U ChEMBL chemreps have no matching assay document provenance.",
        "acyclic_note": "Scaffold-free structures cannot be protected as a chemical family by the Murcko key; individually grouped, not filtered.",
        "seed": seed, "trials_requested": trials, "winning_trial": winning_trial,
        "p_target": dict(zip(SPLITS, target_p)), "u_target": dict(zip(SPLITS, target_u)),
        "p_actual": dict(zip(SPLITS, pc)), "u_actual": dict(zip(SPLITS, uc)),
        "p_source_documents": len(by_doc), "p_documents_present_in_multiple_splits": len(shared_docs),
        "p_molecules_without_document_reference": summary["p_missing_document_refs"],
        "acyclic_without_murcko_scaffold": {
            kind: sum(row[2] for row in info[kind].values()) for kind in ("P", "U")},
        "input_sha256": {k: v for k, v in summary.items() if k.endswith("sha256")},
        "rdkit_version": rdBase.rdkitVersion, "python_version": sys.version,
        "output_sha256": {k: sha256(v) for k, v in outputs.items()},
        "no_new_molecular_filter": True,
        "no_negative_class_assigned": True,
        "caution": "Do NOT claim true specificity/precision or independent negative labels from U; source overlap in audit mode must be interpreted.",
    }
    outputs["split_manifest.json"] = (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    out.mkdir(parents=True, exist_ok=False)
    try:
        for name, raw in outputs.items():
            (out / name).write_bytes(raw)
        with zipfile.ZipFile(zip_path, "x", zipfile.ZIP_DEFLATED) as z:
            for path in sorted(out.iterdir()):
                z.write(path, path.name)
    except BaseException:
        (out / "INCOMPLETE.txt").write_text("Output interrupted; do not use this split.\n", encoding="utf-8")
        raise
    return out, zip_path


def self_test():
    a, b, c = scaffold_key("Cc1ccccc1", "a"), scaffold_key("Oc1ccccc1", "b"), scaffold_key("CCO", "c")
    require(a[0] == b[0] and a[1] == "c1ccccc1", "Murcko ring grouping failed")
    require(c[2] and c[0] == "ACYCLIC:c", "Acyclic handling failed")
    p = [dict(inchikey=f"p{i}", smiles=smiles) for i, smiles in enumerate([
        "Cc1ccccc1", "c1ccncc1", "C1CCCCC1", "c1ccoc1", "c1ccsc1", "c1ccc2ccccc2c1"])]
    u = [dict(inchikey=f"u{i}", smiles=smiles) for i, smiles in enumerate([
        "Oc1ccccc1", "Fc1ccncc1", "CO", "CCO", "CCCC", "CCN", "CCCl", "CCBr", "C#CC", "CCF"])]
    refs = defaultdict(set, {"p0": {"ChEMBL:DOCA"}, "p1": {"ChEMBL:DOCA"}, "p2": {"BINDINGDB:US1"}})
    for mode in ("audit", "strict"):
        groups, info = prepare_groups(p, u, refs, mode)
        result1, qual1 = allocate(groups, (4, 1, 1), (6, 2, 2), 17, 10, refs, p, mode)
        result2, qual2 = allocate(groups, (4, 1, 1), (6, 2, 2), 17, 10, refs, p, mode)
        require(result1[:3] == result2[:3] and qual1 == qual2, "Determinism failed")
        assn = result1[0]
        require(assn[("P", "p0")] == assn[("U", "u0")], "Cross P/U scaffold leaked")
        require(assn[("P", "p1")] == assn[("U", "u1")], "Cross P/U second scaffold leaked")
        if mode == "strict":
            require(assn[("P", "p0")] == assn[("P", "p1")], "Strict source grouping failed")
            require(not result1[4], "Strict source documents crossed partitions")
    with tempfile.TemporaryDirectory() as td:
        groups, info = prepare_groups(p, u, refs, "audit")
        result, quality = allocate(groups, (4, 1, 1), (6, 2, 2), 17, 10, refs, p, "audit")
        pfields, ufields = ["inchikey", "smiles"], ["inchikey", "smiles"]
        folder = Path(td) / "proposed_split"
        written, archive = save_output(Path(td), folder, p, pfields, u, ufields,
                                       info, refs, {"p_missing_document_refs": 0, "p_csv_sha256": "test", "u_csv_sha256": "test"},
                                       (4, 1, 1), (6, 2, 2), "audit", 17, 10, result)
        with zipfile.ZipFile(archive) as z:
            require(set(z.namelist()) >= {"P_train.csv", "U_eval.csv", "split_assignments.csv", "split_manifest.json"},
                    "Output ZIP incomplete")
            split_map = list(csv.DictReader(io.StringIO(z.read("split_assignments.csv").decode("utf-8-sig"))))
            require(len(split_map) == len(p) + len(u), "Output membership loss")
        try:
            save_output(Path(td), folder, p, pfields, u, ufields,
                        info, refs, {"p_missing_document_refs": 0},
                        (4, 1, 1), (6, 2, 2), "audit", 17, 10, result)
        except ValueError:
            pass
        else:
            raise AssertionError("Existing output was overwritten")
    print("Self-tests: PASS (Murcko, acyclic, P/U isolation, source modes, reproducibility, ZIP, no overwrite)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, help="Repository root (e.g. C:\\trpa1_data)")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--trials", type=int, default=24)
    ap.add_argument("--p-val", type=int, default=240, help="Approximate P_val target (proposal)")
    ap.add_argument("--p-test", type=int, default=240, help="Approximate P_test target (proposal)")
    ap.add_argument("--u-train", type=int, default=DEFAULT_U[0])
    ap.add_argument("--u-val", type=int, default=DEFAULT_U[1])
    ap.add_argument("--u-eval", type=int, default=DEFAULT_U[2])
    ap.add_argument("--source-mode", choices=("audit", "strict"), default="audit",
                    help="audit: report/softly discourage P source crossings; strict: group each P document")
    ap.add_argument("--out", type=Path, help="NEW output folder; existing path refused")
    ap.add_argument("--write", action="store_true", help="Create split files after diagnostics")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    require(args.trials >= 1, "--trials must be positive")
    ptarget = (2383 - args.p_val - args.p_test, args.p_val, args.p_test)
    utarget = (args.u_train, args.u_val, args.u_eval)
    require(all(n > 0 for n in ptarget) and all(n > 0 for n in utarget), "Targets must all be positive")
    require(sum(utarget) == 50000, "U targets must add to exactly 50000")
    root = args.root.resolve() if args.root else find_root()
    p, pfields, u, ufields, refs, summary = load_data(root)
    print("Input membership and CSV SHA256: MATCH approved P v1 / frozen U v1")
    print("Scaffold calculation from stored SMILES (no dataset modifications)...")
    groups, detail = prepare_groups(p, u, refs, args.source_mode)
    results, quality = allocate(groups, ptarget, utarget, args.seed, args.trials, refs, p, args.source_mode)
    assn, pc, uc, sources, shared, winning_trial = results
    require(len(assn) == 52383, "Missing split assignments")
    # Never violate P/U scaffold isolation across partitions.
    assigned_groups = defaultdict(set)
    for kind, d in detail.items():
        for key, (group_id, _, _) in d.items():
            assigned_groups[group_id].add(assn[(kind, key)])
    require(all(len(v) == 1 for v in assigned_groups.values()), "Scaffold/source group leakage")
    if args.source_mode == "strict":
        require(not shared, "Strict mode produced shared P source documents")
    print(f"Groups: {len(groups)}; trial: {winning_trial}/{args.trials}; mode: {args.source_mode}")
    print(f"P train/val/test: {pc} (proposed targets {ptarget})")
    print(f"U train/val/eval: {uc} (proposed targets {utarget})")
    print("P source documents shared between partitions:", len(shared), "/", len(sources))
    print("P without identifiable evidence documents:", summary["p_missing_document_refs"])
    print("Acyclic, no Murcko scaffold (P, U):",
          sum(v[2] for v in detail["P"].values()), sum(v[2] for v in detail["U"].values()))
    print("P/U identical Murcko scaffold crossing partitions: 0 (by construction)")
    if quality[0] != 0:
        print("WARNING: Exact U allocation is infeasible with this heuristic/these groups; no scaffold is split.")
    if shared:
        print("WARNING: Some P evidence sources occur in several partitions; see source audit before ML.")
    if not args.write:
        print("READ-ONLY PREVIEW. No files written. Use --write only after reviewing this split proposal.")
        return 0
    require(all(pc) and all(uc), "Cannot write: empty partition")
    require(quality[0] == 0, "Cannot write: requested U targets not met without scaffold leakage")
    out = args.out.resolve() if args.out else root / OUT_DIR
    require(not out.is_relative_to(root / P_DIR) and not out.is_relative_to(root / U_DIR),
            "Refusing output inside P/U input directories")
    output, zipped = save_output(root, out, p, pfields, u, ufields, detail, refs, summary,
                                 ptarget, utarget, args.source_mode, args.seed, args.trials, results)
    require(sha256((root / P_DIR / "trpa1_positive_v1.csv").read_bytes()) == summary["p_csv_sha256"], "P input changed during run")
    require(sha256((root / U_DIR / "trpa1_unlabeled_U_v1.csv").read_bytes()) == summary["u_csv_sha256"], "U input changed during run")
    print("SPLIT PROPOSAL WRITTEN:", output)
    print("ZIP:", zipped)
    print("Original P/U memberships, labels, measurements and files were not edited.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted; input files were not changed.", file=sys.stderr)
        raise SystemExit(130)
    except (OSError, ValueError, KeyError, RuntimeError, json.JSONDecodeError) as exc:
        print("ERROR:", exc, file=sys.stderr)
        raise SystemExit(2)
