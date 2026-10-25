#!/usr/bin/env python3
"""Compare cached BindingDB human-TRPA1 records with the original 1645 dataset.

Place in experiments/02_classification/code/ and run:
    python compare_bindingdb_surechembl.py
Optional, separate, network-dependent SureChEMBL check:
    python compare_bindingdb_surechembl.py --surechembl

Offline mode reads the completed BindingDB extraction made by
trpa1_pubchem_bindingdb_check.py (folder OR ZIP); it does NOT repeat that download.
The overall external run may be incomplete because PubChem was skipped. A valid
BindingDB-specific completion marker and its data hash are required.

This is a versioned source-comparison audit, NOT automatic construction of P/U.
No Git operations. No existing datasets are modified. All new results have a
unique directory. SureChEMBL patent mentions never become positive labels.
Dependencies: pandas, rdkit; requests only with --surechembl.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import html
import io
import json
import math
import re
import sys
import time
import zipfile
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any

try:
    import pandas as pd
    from rdkit import Chem, rdBase
    from rdkit.Chem import SaltRemover
except ImportError as exc:
    raise SystemExit("Missing dependency. Run: python -m pip install pandas rdkit") from exc

VERSION = "1.0.0"
REVIEW_DATE = "2026-10-17"
RAW_MEMBER = "raw/bindingdb_target_records.json"
COMPLETE_MEMBER = "raw/bindingdb_complete.json"
# This is a manually reviewed, versioned evidence lookup. It is not a parser
# pretending to independently validate every patent structure or assay.
SOURCE_EVIDENCE = {
    "US10626112": ("US10626112B2", "Experimental Examples 1-1 and 1-2; Tables 40 and 41", "Human TRPA1 expressed in HEK293T; AITC-induced calcium response; inhibitory IC50 in micromolar."),
    "US12391660": ("US12391660B2", "FLIPR assay/Table 3.1; patch-clamp assay/Table 3.2", "Explicit human TRPA1 inhibitory assays; initial compound response also monitored. Calcium signal and patch-clamp endpoints remain separate records."),
    "US20240101561": ("US20240101561A1", "Assay A; Table 1", "Human TRPA1-expressing HEK293; AITC FLIPR; table heading hTRPA1 IC50 (nM)."),
    "US20250154139": ("US20250154139A1", "Example 25: TRPA1 inhibitor activity assay; Table A", "Inducible human TRPA1 CHO cells; SyncroPatch; carvacrol-evoked current. Most values are reported as potency bounds, not exact IC50."),
    "US20250270198": ("US20250270198A1", "Example 6: TRPA1 inhibitor activity assay; activity table", "Inducible human TRPA1 CHO cells; SyncroPatch; carvacrol-evoked current. BindingDB assay description is missing; primary patent supplies the assay context."),
    "US11701353": ("US11891403B2", "Assay A; Table 1, examples 1-24", "BindingDB Patent Number disagrees with its ligand names. Names identify US11891403; all 24 example/value pairs match Table 1 of US11891403B2. Keep raw identifier unchanged and record the resolved supporting patent separately."),
}
CORRECTED_TABLE_NM = [53,15,23,32,39,45,60,62,73,80,81,115,145,191,227,237,248,153,214,140,155,74,52,100]
REF_RAW_SHA = "ddc59d451b645cadd01deceb8557b6ae07867f9969f12da649942711538959dc"
REF_BASE_KEY_SHA = "d5afeb066d87f12ec18509c68fe9cd75a7568910ed84005358de483d55e61d9e"
REF_NEW_KEY_SHA = "efe791b46151b77392d4be7e6683632658d6ba1d437253d55f76d073631eaf7c"
FIELDS = {
    "BindingDB Reactant_set_id":"rsid", "BindingDB MonomerID":"monomer_id",
    "BindingDB Ligand Name":"ligand_name", "Ligand SMILES":"original_smiles",
    "Ligand InChI Key":"bindingdb_inchikey",
    "Target Name":"target_name",
    "Target Source Organism According to Curator or DataSource":"organism",
    "UniProt (SwissProt) Primary ID of Target Chain 1":"accession",
    "Number of Protein Chains in Target (>1 implies a multichain complex)":"target_chains",
    "IC50 (nM)":"ic50_raw_nM", "Ki (nM)":"ki_raw_nM", "Kd (nM)":"kd_raw_nM", "EC50 (nM)":"ec50_raw_nM",
    "Patent Number":"patent_raw", "Article DOI":"article_doi", "PMID":"pmid",
    "ChEMBL ID of Ligand":"chembl_id", "Curation/DataSource":"curation_source",
    "BindingDB Entry DOI":"bindingdb_entry_doi", "PubChem AID":"pubchem_aid",
}
REMOVER = SaltRemover.SaltRemover()
rdBase.DisableLog("rdApp.info")
rdBase.DisableLog("rdApp.warning")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def key_digest(keys) -> str:
    return sha256(("\n".join(sorted(set(keys))) + "\n").encode("utf-8"))


def text(v: Any) -> str:
    return "" if v is None else str(v).strip()


@lru_cache(maxsize=30000)
def standardize(smiles: str) -> tuple[str, str]:
    """Historical procedure used for both sides, never trusting supplied keys."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "", ""
        mol = REMOVER.StripMol(mol, dontRemoveEverything=True)
        fs = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if not fs:
            return "", ""
        mol = max(fs, key=lambda x: x.GetNumHeavyAtoms())
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True), Chem.MolToInchiKey(mol)
    except Exception:
        return "", ""


def parse_ic50(value: str) -> tuple[str, float | None]:
    match = re.fullmatch(r"\s*(<=|>=|<|>|=|≤|≥)?\s*(\d+(?:\.\d*)?|\.\d+)(?:([eE][+-]?\d+))?\s*", value)
    if not match:
        return "unparseable_or_missing", None
    n = float(match.group(2) + (match.group(3) or ""))
    if not math.isfinite(n) or n <= 0:
        return "invalid", None
    return {"≤":"<=", "≥":">="}.get(match.group(1), match.group(1)) or "=", n


def member_bytes(source: Path, member: str) -> bytes:
    if source.is_dir():
        return (source / member).read_bytes()
    with zipfile.ZipFile(source) as z:
        # Support both flat ZIPs and ZIPs with one enclosing run directory.
        names = [n for n in z.namelist() if n == member or n.endswith("/" + member)]
        if len(names) != 1:
            raise ValueError(f"Expected exactly one {member} in {source.name}; got {len(names)}")
        return z.read(names[0])


def find_repo() -> Path | None:
    for anchor in [Path.cwd(), Path(__file__).resolve().parent]:
        for parent in [anchor, *anchor.parents]:
            if (parent / "experiments/01_regression/data/trpa1_primary_dataset.csv").is_file():
                return parent
    return None


def discover_input(results: Path, kind: str) -> Path | None:
    if not results.is_dir():
        return None
    matches = []
    for p in results.iterdir():
        if kind == "bindingdb":
            if not p.name.startswith("TRPA1_external_check_"):
                continue
            try:
                marker = json.loads(member_bytes(p, COMPLETE_MEMBER))
                if marker.get("info", {}).get("retrieval_complete_for_export") is True:
                    matches.append(p)
            except (OSError, ValueError, zipfile.BadZipFile, KeyError):
                continue
        else:
            if not ("chembl" in p.name.lower() and "check_" in p.name.lower()):
                continue
            try:
                marker = json.loads(member_bytes(p, "summary.json"))
                if marker.get("complete") is True and marker.get("human_target") == "CHEMBL6007":
                    matches.append(p)
            except (OSError, ValueError, zipfile.BadZipFile, KeyError):
                continue
    # Deterministic selection by dated run name; prefer a folder over its ZIP.
    return sorted(matches, key=lambda p: (p.stem if p.suffix==".zip" else p.name, p.is_dir()))[-1] if matches else None


def load_bindingdb(source: Path) -> tuple[list[dict], dict, str]:
    meta = json.loads(member_bytes(source, COMPLETE_MEMBER))
    if meta.get("info", {}).get("retrieval_complete_for_export") is not True:
        raise ValueError("BindingDB extraction is not marked complete. Refusing partial counts.")
    raw = member_bytes(source, RAW_MEMBER)
    digest = sha256(raw)
    if not meta.get("selected_sha256") or digest != meta["selected_sha256"]:
        raise ValueError("BindingDB raw-data hash does not match completion marker.")
    data = json.loads(raw)
    if not isinstance(data, list) or not all(isinstance(x, dict) for x in data):
        raise ValueError("Unexpected BindingDB extraction schema.")
    expected = meta.get("info", {}).get("TRPA1_related_records")
    if expected is not None and len(data) != int(expected):
        raise ValueError("BindingDB row count does not match completion marker.")
    return data, meta, digest


def split_ids(s: str) -> set[str]:
    return set(re.findall(r"CHEMBL\d+", s))


def evidence(r: dict) -> dict:
    pat = r["patent_raw"]
    res = {"supporting_patent":"", "source_locator":"", "source_url":"",
           "source_review":"bindingdb_record_only", "source_note":"Primary-source protocol not rechecked in this audit."}
    if pat not in SOURCE_EVIDENCE:
        return res
    support, locator, note = SOURCE_EVIDENCE[pat]
    if pat == "US11701353":
        m = re.search(r"US11891403,\s*Example\s+(\d+)", r["ligand_name"])
        valid = bool(m and 1 <= int(m.group(1)) <= 24 and r["ic50_relation"] == "=" and
                     math.isclose(r["ic50_numeric_nM"], CORRECTED_TABLE_NM[int(m.group(1))-1], rel_tol=0, abs_tol=1e-9))
        if not valid:
            res.update(source_review="unresolved_patent_mismatch", source_note=note)
            return res
        status = "human_inhibition_protocol_reviewed_patent_id_resolved_by_24_value_table"
    else:
        status = "human_inhibition_protocol_reviewed"
    res.update(supporting_patent=support, source_locator=locator,
               source_url=f"https://patents.google.com/patent/{support}/en",
               source_review=status, source_note=note)
    return res


def export_csv(path: Path, rows: list[dict], empty_cols: list[str] | None = None):
    # UTF-8 with BOM for Windows Excel, LF consistently for hashes.
    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=empty_cols or ["status"])
    df.to_csv(path, index=False, encoding="utf-8-sig", lineterminator="\n")


def join_values(rows: list[dict], field: str) -> str:
    return " | ".join(sorted({text(r.get(field)) for r in rows if text(r.get(field))}))


def compare(raw: list[dict], baseline: pd.DataFrame, chembl_activities: list[dict] | None):
    required = {"std_smiles", "inchikey", "molecule_chembl_id"}
    if not required.issubset(baseline.columns):
        raise ValueError(f"Baseline lacks columns: {sorted(required-set(baseline.columns))}")
    if baseline.inchikey.duplicated().any():
        raise ValueError("Baseline contains duplicate InChIKeys.")
    calc = [standardize(s) for s in baseline.std_smiles]
    changed = [(sup, got[1]) for sup, got in zip(baseline.inchikey, calc) if sup != got[1]]
    if changed:
        raise ValueError(f"Historical standardization does not reproduce {len(changed)} baseline identities under RDKit {rdBase.rdkitVersion}. No comparison performed.")
    bkeys = set(baseline.inchikey)
    bconn = {k[:14] for k in bkeys}
    bids = set().union(*(split_ids(s) for s in baseline.molecule_chembl_id))
    ckeys = None
    if chembl_activities is not None:
        ckeys = {standardize(text(r.get("canonical_smiles")))[1] for r in chembl_activities
                 if r.get("target_chembl_id") == "CHEMBL6007" and r.get("target_organism") == "Homo sapiens"}
        ckeys.discard("")
    all_rows = []
    for i, rawrow in enumerate(raw, 1):
        r = {dest: text(rawrow.get(src)) for src, dest in FIELDS.items()}
        r["assay_ids"] = " | ".join(rawrow.get("_entry_assay_ids") or [])
        r["assay_description"] = html.unescape(" | ".join(rawrow.get("_assay_descriptions") or []))
        smi, key = standardize(r["original_smiles"])
        rel, val = parse_ic50(r["ic50_raw_nM"])
        human = (r["organism"].casefold() in {"human", "homo sapiens"} and
                 r["accession"] == "O75762" and r["target_chains"] == "1")
        r.update(standard_smiles=smi, standard_inchikey=key, ic50_relation=rel, ic50_numeric_nM=val,
                 strict_human_single_chain=human,
                 cxsmiles_annotation_present="|" in r["original_smiles"],
                 supplied_key_differs=bool(key and r["bindingdb_inchikey"] and key != r["bindingdb_inchikey"]),
                 novelty_status="",
                 same_chembl_identifier=bool(split_ids(r["chembl_id"]) & bids),
                 same_connectivity_in_baseline=bool(key and key[:14] in bconn),
                 also_in_full_chembl_by_historical_key=(key in ckeys if ckeys is not None and key else None))
        if not human:
            status = "excluded_not_explicit_human_single_chain"
        elif rel not in {"=", "<", "<="}:
            status = "excluded_no_exact_or_upper_bound_ic50"
        elif not key:
            status = "excluded_invalid_structure"
        elif key in bkeys:
            status = "already_in_primary_by_full_key"
        else:
            if r["same_chembl_identifier"]:
                status = "review_same_chembl_id_different_representation"
            elif r["same_connectivity_in_baseline"]:
                status = "review_same_connectivity_stereo_or_isotope_variant"
            else:
                status = "new_to_primary_conservative"
        r["novelty_status"] = status
        r.update(evidence(r))
        r["structure_review"] = "original_CXSMILES_preserved_check_relative_stereo" if r["cxsmiles_annotation_present"] else "not_manually_checked_against_patent_drawing"
        all_rows.append(r)
        if i % 500 == 0:
            print(f"  Compared {i:,}/{len(raw):,} BindingDB records.", flush=True)
    new = [r for r in all_rows if r["novelty_status"] == "new_to_primary_conservative"]
    groups: dict[str, list[dict]] = {}
    for r in new:
        groups.setdefault(r["standard_inchikey"], []).append(r)
    molecules = []
    for key, rows in sorted(groups.items()):
        exact = [r["ic50_numeric_nM"] for r in rows if r["ic50_relation"] == "="]
        reviewed = any(r["source_review"].startswith("human_inhibition_protocol_reviewed") for r in rows)
        molecules.append({
            "inchikey":key, "standard_smiles":rows[0]["standard_smiles"],
            "bindingdb_monomer_ids":join_values(rows,"monomer_id"),
            "bindingdb_reactant_set_ids":join_values(rows,"rsid"),
            "n_records":len(rows), "evidence_type":"exact_ic50_available" if exact else "upper_bound_only",
            "exact_ic50_min_nM":min(exact) if exact else None,
            "exact_ic50_max_nM":max(exact) if exact else None,
            "original_ic50_expressions_nM":join_values(rows,"ic50_raw_nM"),
            "source_level_human_inhibition_reviewed":reviewed,
            "patents_raw":join_values(rows,"patent_raw"), "supporting_patents":join_values(rows,"supporting_patent"),
            "article_dois":join_values(rows,"article_doi"),
            "source_urls":join_values(rows,"source_url"),
            "cxsmiles_annotation_present":any(r["cxsmiles_annotation_present"] for r in rows),
            "original_smiles":join_values(rows,"original_smiles"),
            "also_in_full_chembl_by_historical_key":rows[0]["also_in_full_chembl_by_historical_key"],
            "automatic_P_label_assigned":False,
        })
    srcgroups: dict[tuple[str,str], list[dict]] = {}
    for r in new:
        srcgroups.setdefault((r["patent_raw"],r["article_doi"]),[]).append(r)
    sources = []
    for (pat,doi),rows in sorted(srcgroups.items()):
        keys={r["standard_inchikey"] for r in rows}
        exkeys={r["standard_inchikey"] for r in rows if r["ic50_relation"]=="="}
        sources.append({"patent_raw":pat, "article_doi":doi, "records":len(rows),
                        "unique_structures":len(keys), "with_exact_ic50":len(exkeys),
                        "upper_bound_only":len(keys-exkeys),
                        "supporting_patent":join_values(rows,"supporting_patent"),
                        "source_review":join_values(rows,"source_review"),
                        "source_locator":join_values(rows,"source_locator"),
                        "source_url":join_values(rows,"source_url"),
                        "review_note":join_values(rows,"source_note"), "review_date":REVIEW_DATE})
    counts = {
        "baseline_molecules":len(baseline), "baseline_identity_mismatches":len(changed),
        "bindingdb_TRPA1_related_records":len(raw),
        "strict_human_single_chain_records":sum(r["strict_human_single_chain"] for r in all_rows),
        "human_exact_or_upper_bound_ic50_records":sum(r["strict_human_single_chain"] and r["ic50_relation"] in {"=","<","<="} for r in all_rows),
        "new_records_after_conservative_comparison":len(new), "new_unique_structures":len(molecules),
        "new_with_exact_ic50":sum(r["evidence_type"]=="exact_ic50_available" for r in molecules),
        "new_upper_bound_only":sum(r["evidence_type"]=="upper_bound_only" for r in molecules),
        "new_source_level_human_inhibition_reviewed":sum(r["source_level_human_inhibition_reviewed"] for r in molecules),
        "new_primary_protocol_not_rechecked":sum(not r["source_level_human_inhibition_reviewed"] for r in molecules),
        "new_with_CXSMILES_annotation":sum(r["cxsmiles_annotation_present"] for r in molecules),
        "new_also_in_full_chembl_by_historical_key":sum(r["also_in_full_chembl_by_historical_key"] is True for r in molecules) if ckeys is not None else None,
        "record_disposition":dict(sorted(__import__('collections').Counter(r["novelty_status"] for r in all_rows).items())),
    }
    return all_rows, new, molecules, sources, counts


def surechembl_check(output: Path, cache: Path, baseline: pd.DataFrame, new_molecules: list[dict]) -> dict:
    """Bounded local API check, not an exhaustive target/ligand evidence extraction.

    Search patent titles for TRPA1 (up to 200 hits, status records truncation),
    then export chemistry ONLY for the six source-reviewed patents above.
    No inferred inhibitor labels. Cache every successful response, stop service
    on first failure, no hundreds of retries. The live API is untested here.
    API contract: https://www.surechembl.org/api/v3/api-docs
    """
    try:
        import requests
    except ImportError:
        return {"status":"INCOMPLETE", "error":"Install requests: python -m pip install requests"}
    cache.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "TRPA1-academic-source-audit/1.0"
    last = [0.0]
    def fetch(path: str, params: dict, suffix: str) -> bytes:
        request_key = sha256(json.dumps({"path":path,"params":params},sort_keys=True).encode())
        dest = cache/(request_key+suffix)
        if dest.exists():
            return dest.read_bytes()
        time.sleep(max(0, 1.5-(time.monotonic()-last[0])))
        last[0] = time.monotonic()
        response = session.post("https://www.surechembl.org/api"+path, params=params, timeout=(10,60))
        response.raise_for_status()  # fail fast on 429 instead of hammering server
        payload=response.content
        if suffix == ".zip":
            with zipfile.ZipFile(io.BytesIO(payload)) as z:
                if not any(n.lower().endswith('.csv') for n in z.namelist()):
                    raise ValueError("SureChEMBL export ZIP has no CSV.")
        else:
            env=json.loads(payload)
            if env.get("status")!="OK" or not isinstance(env.get("data"),dict):
                raise ValueError(f"Unexpected SureChEMBL envelope: {str(env)[:160]}")
        tmp=dest.with_suffix(dest.suffix+'.part'); tmp.write_bytes(payload); tmp.replace(dest)
        return payload
    result={"status":"INCOMPLETE", "scope":"title query TRPA1 (maximum 200 hits), chemistry exports of six reviewed patents only; NOT an exhaustive SureChEMBL scan", "query":"ttl:TRPA1", "live_network_tested_by_assistant":False, "validated_new_inhibitors":None}
    leads=[]; pages=[]; documents=[]; exports=[]
    known_keys=set(baseline.inchikey)|{r['inchikey'] for r in new_molecules}
    try:
        for page in [1,2]:
            raw=fetch('/search/content',{'query':'ttl:TRPA1','page':page,'itemsPerPage':100},'.json')
            payload=json.loads(raw)['data']; results=payload.get('results')
            if not isinstance(results,dict) or not isinstance(results.get('documents'),list):
                raise ValueError('Unexpected SureChEMBL search result schema; raw response kept in cache.')
            docs=results['documents']; documents.extend(docs); pages.append(payload)
            total=int(results['total_hits'])
            result['reported_title_hits']=total
            if not docs or len(documents)>=total: break
        result['fetched_title_documents']=len(documents)
        result['title_search_truncated']=len(documents)<result['reported_title_hits']
        (output/'surechembl_patent_search.json').write_text(json.dumps(pages,ensure_ascii=False,indent=2),encoding='utf-8')
        for support in sorted({v[0] for v in SOURCE_EVIDENCE.values()}):
            m=re.fullmatch(r'([A-Z]{2})(\d+)([A-Z]\d)',support)
            if not m: raise ValueError(f'Unsupported patent identifier: {support}')
            doc_id='-'.join(m.groups())
            print(f'  SureChEMBL chemistry export: {doc_id}',flush=True)
            raw=fetch('/export/document-chemistry',{'docID':doc_id,'output_type':'csv'},'.zip')
            count=0
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                for member in z.namelist():
                    if not member.lower().endswith('.csv'):continue
                    reader=csv.DictReader(io.StringIO(z.read(member).decode('utf-8-sig')))
                    names={str(k).strip().lower():k for k in reader.fieldnames or []}
                    scol=names.get('smiles') or names.get('canonical_smiles')
                    if not scol: raise ValueError(f'No recognized SMILES column in {member}; raw export cached.')
                    for row in reader:
                        count+=1
                        original=text(row.get(scol)); smi,key=standardize(original)
                        leads.append({'patent':support,'surechembl_id':row.get(names.get('id','')) or row.get(names.get('surechembl_id','')) or '',
                                      'original_smiles':original,'standard_smiles':smi,'inchikey':key,
                                      'absent_from_baseline_and_bindingdb_candidates':bool(key and key not in known_keys),
                                      'activity_evidence':'PATENT_MENTION_ONLY_NO_COMPOUND_ACTIVITY_LINK',
                                      'automatic_P_label_assigned':False})
            exports.append({'patent':support,'export_rows':count,'payload_sha256':sha256(raw)})
        result['status']='TARGETED_EXPORTS_COMPLETE_NOT_EXHAUSTIVE'
    except (requests.RequestException,OSError,ValueError,KeyError,zipfile.BadZipFile) as exc:
        result['error']=f'{type(exc).__name__}: {str(exc)[:350]}'
        print('  SureChEMBL stopped; cached successes kept. '+result['error'],flush=True)
    result['exports']=exports
    result['patent_mentions_rows']=len(leads)
    result['new_mentioned_structure_keys']=len({r['inchikey'] for r in leads if r['absent_from_baseline_and_bindingdb_candidates']})
    if leads: export_csv(output/'surechembl_structure_leads.csv',leads)
    return result


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--baseline',type=Path,help='Original molecule-level CSV; normally discovered in repo.')
    ap.add_argument('--bindingdb',type=Path,help='Completed prior external-check folder or ZIP.')
    ap.add_argument('--chembl',type=Path,help='Completed live ChEMBL check folder or ZIP; optional full-snapshot overlap check.')
    ap.add_argument('--output-dir',type=Path,help='NEW output directory; must not exist.')
    ap.add_argument('--surechembl',action='store_true',help='Optional bounded, network-dependent SureChEMBL search/exports.')
    args=ap.parse_args()
    repo=find_repo()
    resultroot=repo/'experiments/02_classification/results' if repo else Path.cwd()/'results'
    bpath=args.baseline or (repo/'experiments/01_regression/data/trpa1_primary_dataset.csv' if repo else None)
    external=args.bindingdb or discover_input(resultroot,'bindingdb')
    chsource=args.chembl or discover_input(resultroot,'chembl')
    if bpath is None or not bpath.is_file():
        ap.error('Cannot find original dataset; use --baseline PATH or place script in the repo code directory.')
    if external is None:
        ap.error('Cannot find completed BindingDB cache; use --bindingdb PATH_TO_PRIOR_FOLDER_OR_ZIP. Do not redownload the export.')
    print(f'BindingDB input: {external}\nBaseline:        {bpath}\nRDKit:           {rdBase.rdkitVersion}\nDefault mode is OFFLINE; no PubChem access.',flush=True)
    if args.surechembl:
        print('SureChEMBL mode requested: bounded network check AFTER the offline comparison.',flush=True)
    raw,meta,rawsha=load_bindingdb(external)
    bbytes=bpath.read_bytes(); bframe=pd.read_csv(io.BytesIO(bbytes),dtype=str,keep_default_na=False)
    activities=None; chmeta=None; chsha=None
    if chsource:
        chmeta=json.loads(member_bytes(chsource,'summary.json'))
        if chmeta.get('complete') is not True or chmeta.get('human_target')!='CHEMBL6007':
            raise ValueError('ChEMBL optional input is incomplete or the wrong target.')
        cb=member_bytes(chsource,'human_trpa1_activities.json'); chsha=sha256(cb); activities=json.loads(cb)
    allrows,new,molecules,sources,counts=compare(raw,bframe,activities)
    output=args.output_dir or resultroot/('bindingdb_surechembl_comparison_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ'))
    output.mkdir(parents=True,exist_ok=False)
    try:
        export_csv(output/'new_bindingdb_molecules.csv',molecules,['inchikey','standard_smiles'])
        export_csv(output/'new_bindingdb_records.csv',new,['rsid','standard_inchikey'])
        export_csv(output/'comparison_audit.csv',allrows,['rsid','novelty_status'])
        export_csv(output/'source_summary.csv',sources,['patent_raw','unique_structures'])
        digest=key_digest([r['inchikey'] for r in molecules])
        reference = ('MATCH' if digest==REF_NEW_KEY_SHA else 'DIFFERENT_CHECK_TABLES_AND_RDKIT_VERSION') if rawsha==REF_RAW_SHA and key_digest(bframe.inchikey)==REF_BASE_KEY_SHA else 'DIFFERENT_INPUT_NOT_REFERENCE_RUN'
        sc={'status':'NOT_RUN','validated_new_inhibitors':None,'reason':'Offline BindingDB comparison only. SureChEMBL data access did not succeed in the assistant environment; not evidence of zero compounds.'}
        if args.surechembl:
            cacheroot=(repo/'experiments/02_classification/data/surechembl_cache') if repo else output.parent/'surechembl_cache'
            sc=surechembl_check(output,cacheroot,bframe,molecules)
        summary={
            'script_version':VERSION,'created_utc':datetime.now(timezone.utc).isoformat(),
            'comparison_complete':True, 'bindingdb_release':meta.get('info',{}).get('release'),
            'bindingdb_cache_path':str(external.resolve()),'bindingdb_raw_sha256':rawsha,
            'baseline_path':str(bpath.resolve()),'baseline_file_sha256':sha256(bbytes),
            'baseline_key_set_sha256':key_digest(bframe.inchikey),
            'chembl_full_snapshot_path':str(chsource.resolve()) if chsource else None,
            'chembl_full_snapshot_sha256':chsha,
            'python':sys.version,'rdkit':rdBase.rdkitVersion,'pandas':pd.__version__,
            'script_sha256':sha256(Path(__file__).read_bytes()),
            'counts':counts,'new_structure_key_set_sha256':digest,'reference_key_set_comparison':reference,
            'surechembl':sc,
            'methods':[
                'Both sides: RDKit parse; SaltRemover; largest fragment; canonical isomeric SMILES; recompute InChIKey.',
                'Human filter: organism Human/Homo sapiens, UniProt O75762, one target chain. Multichain, unknown and nonhuman entries excluded from conservative main count.',
                'Candidate IC50 relations =, <, <= with strictly positive finite values in nM. No potency cutoff. > bounds do not constitute positive evidence by this rule.',
                'Exclude baseline full-key matches. Different-key matches by existing ChEMBL ID or InChI connectivity block go to review, not to the new count. This conservative rule also sets aside possible new stereoisomers; it is not exhaustive chemical-identity curation.',
                'Novelty is relative to the historical 1645 primary dataset, not novelty to chemistry or to all ChEMBL.',
                'Manual source-level patent review is encoded with URL/locator and dated. This is not automated validation of every structure drawing, assay outcome, stereoisomer, selectivity or direct antagonism.',
                'Original CXSMILES is preserved; its stereochemical annotations must be considered before final ML input. Current comparison uses historical simple isomeric SMILES for compatibility.',
                'These outputs propose additions for curation; no final P/U labels or model results are generated.'
            ]}
        (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        zp=output.with_suffix('.zip')
        with zipfile.ZipFile(zp,'x',zipfile.ZIP_DEFLATED) as z:
            for p in sorted(output.iterdir()):
                if p.is_file(): z.write(p,p.name)
        print('\nRESULTS')
        for k,v in counts.items():
            if not isinstance(v,dict): print(f'  {k}: {v}')
        print(f'  Reference key-set check: {reference}\n  SureChEMBL: {sc["status"]}\nResults: {output}\nZIP: {zp}')
        if sc['status']=='INCOMPLETE':
            print('BindingDB comparison succeeded. SureChEMBL failed/partial; do not interpret its count as complete.')
        return 0
    except BaseException as exc:
        (output/'INCOMPLETE.json').write_text(json.dumps({'error':str(exc),'complete':False}),encoding='utf-8')
        raise

if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (OSError,ValueError,KeyError,zipfile.BadZipFile) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        raise SystemExit(1)
