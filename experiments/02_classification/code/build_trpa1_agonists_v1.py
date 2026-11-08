#!/usr/bin/env python3
"""Collect evidence-backed human TRPA1 agonists from frozen ChEMBL 37 and BindingDB 202610.

Run from experiments/02_classification/code/:
    python build_trpa1_agonists_v1.py --self-test
    python build_trpa1_agonists_v1.py

Inputs (READ ONLY):
    experiments/02_classification/results/chembl_live_check_20261011T181249_632407Z.zip
    experiments/02_classification/results/bindingdb_surechembl_comparison_20261025T145527_934656Z.zip
    experiments/02_classification/data/positive_P_v1_2383/trpa1_positive_v1.csv

ChEMBL assay descriptions are fetched from the ChEMBL API (READ ONLY) only for
EC50 candidate assays, with a strict ChEMBL_37 release check. To repeat offline,
pass --assays-cache PATH pointing to chembl_assay_metadata.json from this output.

This script identifies a documented subset of experimentally reported agonists,
NOT all possible activators. To avoid mislabeling antagonists, it accepts only
EC50 (=, <, <=; nM; positive number) in an assay that explicitly describes
TRPA1 agonism / activation, with human target metadata. It does NOT use EC50
alone to imply agonism. It retains overlapping P molecules as dual-action records.

No potency cutoff, new structural standardization, class N=0, or edits to P/U.
The agonist evidence is derived from database assay annotations; this is NOT
an independent primary-paper/compound-drawing audit. Outputs are PROPOSED for
researcher review and not an automatic change to the existing U dataset.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile

from rdkit import Chem, rdBase
from rdkit.Chem import SaltRemover

P_REL = Path('experiments/02_classification/data/positive_P_v1_2383/trpa1_positive_v1.csv')
CHEMBL_REL = Path('experiments/02_classification/results/chembl_live_check_20261011T181249_632407Z.zip')
BDB_REL = Path('experiments/02_classification/results/bindingdb_surechembl_comparison_20261025T145527_934656Z.zip')
OUT_REL = Path('experiments/02_classification/data/agonist_A_v1')
REMOVER = SaltRemover.SaltRemover()
API = 'https://www.ebi.ac.uk/chembl/api/data/'

# Deliberately restrictive assay interpretation. EC50 can describe inhibition.
AG = re.compile(r'\bagonist\s+activity\b|\bactivation\s+of\s+(?:human\s+)?TRPA1\b|\bTRPA1\s+activation\b', re.I)
NOT_AG = re.compile(r'\bantagonis\w*|\binhibit\w*|\bdesensiti\w*|\bdecrease\b|\breduction\b|\bblock(?:ade|er|ing)?\b', re.I)
NOT_HUMAN = re.compile(r'\bunknown origin\b|\bmutant\b|\bchimer\w*|\brat\b|\bmouse\b|\bmurine\b', re.I)
# Potentiation of another agonist's induced signal is not clear direct agonism.
CHALLENGE = re.compile(r'\b(?:AITC|menthol|cinnamaldehyde|mustard oil|allyl isothiocyanate)[ -]induced\b', re.I)
MEASUREMENT = re.compile(r'^\s*(<=|>=|<|>|=|≤|≥)?\s*((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$')


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(cond, msg):
    if not cond:
        raise ValueError(msg)


def root_path():
    for start in (Path(__file__).resolve().parent, Path.cwd()):
        for p in (start, *start.parents):
            if (p / P_REL).is_file():
                return p
    raise FileNotFoundError('Repository not found. Supply --root C:\\trpa1_data')


def read_zip(path: Path, filename: str) -> bytes:
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n == filename or n.endswith('/' + filename)]
        require(len(names) == 1, f'{path}: expected one {filename}; found {names}')
        return z.read(names[0])


def csv_rows(raw: bytes) -> list[dict]:
    with io.StringIO(raw.decode('utf-8-sig'), newline='') as buf:
        reader = csv.DictReader(buf)
        require(reader.fieldnames is not None and len(reader.fieldnames) == len(set(reader.fieldnames)), 'Invalid CSV headers')
        records = list(reader)
    require(all(None not in r and all(v is not None for v in r.values()) for r in records), 'Malformed CSV row')
    return records


def csv_bytes(rows: list[dict]) -> bytes:
    require(bool(rows), 'Cannot write empty result')
    fields = list(rows[0])
    require(all(set(r) == set(fields) for r in rows), 'Inconsistent CSV rows')
    with io.StringIO(newline='') as buf:
        writer = csv.DictWriter(buf, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
        return buf.getvalue().encode('utf-8-sig')


def standardize(smiles: str) -> tuple[str, str]:
    """EXACT historical identity: parse, salt removal, largest fragment, isomeric SMILES."""
    if not smiles:
        return '', ''
    try:
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return '', ''
            mol = REMOVER.StripMol(mol, dontRemoveEverything=True)
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
            if not frags:
                return '', ''
            if len(frags) > 1:
                mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True), Chem.MolToInchiKey(mol)
    except (ValueError, RuntimeError, IndexError):
        return '', ''


def parse_value(value: str) -> tuple[str, float | None]:
    m = MEASUREMENT.fullmatch(str(value or ''))
    if not m:
        return '', None
    f = float(m[2])
    if not math.isfinite(f) or f <= 0:
        return '', None
    return {'≤': '<=', '≥': '>='}.get(m[1], m[1]) or '=', f


def context_status(description: str) -> str:
    """Only direction explicitly labelled as TRPA1 agonist/activation is accepted."""
    d = str(description or '').split('\t')[-1].strip()
    if not d:
        return 'missing_assay_description'
    if NOT_HUMAN.search(d):
        return 'ambiguous_species_or_variant'
    if NOT_AG.search(d):
        return 'inhibitory_or_desensitization_context'
    if CHALLENGE.search(d):
        return 'other_agonist_induced_signal'
    if not AG.search(d):
        return 'not_explicit_agonist_assay'
    return 'explicit_agonist_assay'


def http_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={'User-Agent': 'TRPA1-agonists-frozen-ChEMBL37/1.0', 'Accept': 'application/json'})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                data = json.load(resp)
            require(isinstance(data, dict), f'Unexpected API payload at {url}')
            return data
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError('Unreachable')


def get_assays(ids: set[str], cache: Path | None) -> tuple[dict, bytes]:
    """Use fully saved assay metadata or fetch ONLY when API reports ChEMBL_37."""
    if cache:
        raw = cache.read_bytes()
        payload = json.loads(raw)
        require(payload.get('chembl_db_version') == 'ChEMBL_37', 'Assay cache is not ChEMBL_37.')
        assays = payload.get('assays', {})
    else:
        status = http_json(API + 'status.json')
        require(status.get('chembl_db_version') == 'ChEMBL_37',
                'ChEMBL API release is not 37. STOP; do not mix releases. Use --assays-cache from ChEMBL_37.')
        assays = {}
        for i, start in enumerate(range(0, len(ids), 25), 1):
            batch = sorted(ids)[start:start + 25]
            url = API + 'assay.json?' + urllib.parse.urlencode({
                'assay_chembl_id__in': ','.join(batch), 'limit': 1000})
            result = http_json(url)
            for a in result.get('assays', []):
                aid = a.get('assay_chembl_id')
                if aid in ids:
                    assays[aid] = a
            print(f'  ChEMBL assay metadata batches: {i}, retrieved {len(assays)}/{len(ids)}', flush=True)
        payload = {'chembl_db_version': 'ChEMBL_37', 'assays': assays,
                   'metadata_source': 'ChEMBL API 37; each assay ID was requested explicitly'}
        raw = (json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2) + '\n').encode('utf-8')
    missing = sorted(ids - assays.keys())
    require(not missing, f'Assay descriptions missing for {len(missing)} ChEMBL assays: {missing[:12]}')
    return assays, raw


def candidate_chembl_activities(records: list[dict]) -> tuple[list[dict], Counter]:
    out, reasons = [], Counter()
    for r in records:
        # Activity records have no reliable activation/inhibition direction alone.
        if str(r.get('standard_type', '')).upper() != 'EC50':
            continue
        reasons['chembl_ec50_records'] += 1
        if r.get('target_chembl_id') != 'CHEMBL6007' or r.get('target_organism') not in ('Homo sapiens', None):
            reasons['wrong_target_or_species'] += 1
            continue
        rel = str(r.get('standard_relation') or '').strip()
        val = str(r.get('standard_value') or '').strip()
        parsed_rel, number = parse_value(rel + val)
        if parsed_rel not in ('=', '<', '<=') or number is None:
            reasons['not_exact_or_upper_ec50'] += 1
            continue
        if r.get('standard_units') != 'nM':
            reasons['not_nM'] += 1
            continue
        if str(r.get('data_validity_comment') or '') not in ('', 'Manually validated'):
            reasons['quality_flagged'] += 1
            continue
        if not r.get('canonical_smiles') or not r.get('assay_chembl_id'):
            reasons['missing_smiles_or_assay'] += 1
            continue
        out.append(r)
    return out, reasons


def collect(chembl: list[dict], assays: dict, bdb: list[dict], p_keys: set[str]):
    accepted, source_counts, flags = [], Counter(), Counter()
    # Preserve every accepted source row (no averaging of censored measurements).
    for r in chembl:
        aid = str(r['assay_chembl_id'])
        a = assays[aid]
        direction = context_status(a.get('description', ''))
        if a.get('target_chembl_id') not in (None, 'CHEMBL6007'):
            direction = 'assay_target_mismatch'
        if a.get('assay_variant_mutation') or a.get('assay_variant_accession'):
            direction = 'assay_variant'
        source_counts['chembl_' + direction] += 1
        if direction != 'explicit_agonist_assay':
            continue
        smiles, key = standardize(str(r.get('canonical_smiles') or ''))
        if not key:
            flags['chembl_failed_structure'] += 1
            continue
        rel, value = parse_value(str(r.get('standard_relation') or '') + str(r.get('standard_value') or ''))
        accepted.append(dict(inchikey=key, smiles=smiles, source='ChEMBL',
            record_id=str(r['activity_id']), molecule_id=str(r.get('molecule_chembl_id') or ''),
            original_smiles=str(r.get('canonical_smiles') or ''), assay_id=aid,
            assay_description=str(a.get('description') or ''), relation=rel, ec50_nM=value,
            original_ec50_expression_nM=str(r.get('standard_relation') or '') + str(r.get('standard_value') or ''),
            supporting_document=str(r.get('document_chembl_id') or ''),
            source_url='https://www.ebi.ac.uk/chembl/explore/activity/' + str(r['activity_id']),
            in_inhibitor_P_v1=(key in p_keys),
            raw_record_json=json.dumps(r, ensure_ascii=False, sort_keys=True, separators=(',', ':'))))
    for r in bdb:
        source_counts['bindingdb_rows_seen'] += 1
        human = (str(r.get('organism') or '').casefold() in ('human', 'homo sapiens') and
                 r.get('accession') == 'O75762' and r.get('target_chains') == '1')
        if not human:
            source_counts['bindingdb_nonhuman_or_not_single_chain'] += 1
            continue
        relation, val = parse_value(str(r.get('ec50_raw_nM') or ''))
        if relation not in ('=', '<', '<=') or val is None:
            source_counts['bindingdb_no_exact_or_upper_ec50'] += 1
            continue
        direction = context_status(r.get('assay_description', ''))
        source_counts['bindingdb_' + direction] += 1
        if direction != 'explicit_agonist_assay':
            continue
        smiles, key = standardize(str(r.get('original_smiles') or ''))
        if not key:
            flags['bindingdb_failed_structure'] += 1
            continue
        if r.get('standard_inchikey') and r['standard_inchikey'] != key:
            flags['bindingdb_recomputed_key_differs'] += 1
            # Do not silently attach a different identity to a measured value.
            continue
        accepted.append(dict(inchikey=key, smiles=smiles, source='BindingDB',
            record_id=str(r.get('rsid') or ''), molecule_id=str(r.get('monomer_id') or ''),
            original_smiles=str(r.get('original_smiles') or ''), assay_id=str(r.get('assay_ids') or ''),
            assay_description=str(r.get('assay_description') or ''), relation=relation, ec50_nM=val,
            original_ec50_expression_nM=str(r.get('ec50_raw_nM') or ''),
            supporting_document=str(r.get('article_doi') or r.get('patent_raw') or ''),
            source_url=('https://doi.org/' + r['article_doi']) if r.get('article_doi') else
                       ('https://patents.google.com/patent/' + r['patent_raw'] + '/en' if r.get('patent_raw') else ''),
            in_inhibitor_P_v1=(key in p_keys),
            raw_record_json=json.dumps(r, ensure_ascii=False, sort_keys=True, separators=(',', ':'))))
    # Same measurement imported into both databases is not treated as an independent replicate.
    by_key = defaultdict(list)
    for r in accepted:
        by_key[r['inchikey']].append(r)
    molecules = []
    for key, rows in sorted(by_key.items()):
        chosen = min(rows, key=lambda r: (r['source'] != 'ChEMBL', r['record_id']))
        molecules.append(dict(inchikey=key, smiles=chosen['smiles'],
            in_inhibitor_P_v1=(key in p_keys), sources=';'.join(sorted(set(r['source'] for r in rows))),
            source_record_ids=';'.join(sorted(set(r['source'] + ':' + r['record_id'] for r in rows))),
            number_of_database_rows=len(rows)))
    accepted.sort(key=lambda r: (r['inchikey'], r['source'], r['record_id']))
    return molecules, accepted, dict(source_counts), dict(flags)


def self_test():
    assert context_status('Agonist activity at human TRPA1, increased Ca2+ influx') == 'explicit_agonist_assay'
    assert context_status('Antagonist activity at human TRPA1 assessed as EC50') == 'inhibitory_or_desensitization_context'
    assert context_status('Agonist activity at TRPA1 (unknown origin)') == 'ambiguous_species_or_variant'
    assert context_status('Agonist activity at human TRPA1 assessed as inhibition of AITC current') == 'inhibitory_or_desensitization_context'
    assert context_status('Agonist activity at human TRPA1 assessed as increase in menthol-induced calcium signal') == 'other_agonist_induced_signal'
    assert parse_value('<=1000') == ('<=', 1000.0)
    assert parse_value('>1000')[0] == '>'
    assert standardize('CCO.[Na+]') == standardize('CCO')
    assert standardize('C[C@H](O)C(=O)O')[1] != standardize('C[C@@H](O)C(=O)O')[1]
    chembl = [dict(activity_id='12', assay_chembl_id='CHEMBL1', molecule_chembl_id='CHEMBL2',
                   canonical_smiles='CCO', standard_relation='=', standard_value='1000',
                   standard_units='nM', standard_type='EC50', target_chembl_id='CHEMBL6007',
                   target_organism='Homo sapiens', data_validity_comment='', document_chembl_id='CHEMBL3')]
    assays = {'CHEMBL1': {'target_chembl_id': 'CHEMBL6007', 'description': 'Agonist activity at human TRPA1 by calcium influx'}}
    b = [dict(rsid='777', monomer_id='5', original_smiles='CCO', standard_inchikey=standardize('CCO')[1],
              organism='Homo sapiens', accession='O75762', target_chains='1', ec50_raw_nM='1000',
              assay_description='Agonist activity at human TRPA1', assay_ids='1', article_doi='')]
    c, stats = candidate_chembl_activities(chembl)
    assert len(c) == 1 and stats['chembl_ec50_records'] == 1
    m, e, _, _ = collect(c, assays, b, set())
    assert len(m) == 1 and len(e) == 2 and m[0]['sources'] == 'BindingDB;ChEMBL'
    b[0]['assay_description'] = 'Antagonist activity at human TRPA1'
    m, e, _, _ = collect(c, assays, b, set())
    assert len(m) == 1 and len(e) == 1
    print('Self-tests: PASS (agonist-vs-antagonist EC50, species, stereo, salts, duplicates, provenance).')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--root', type=Path, help='Repository root, e.g. C:\\trpa1_data')
    p.add_argument('--chembl-snapshot', type=Path, help='Frozen chembl_live_check ZIP')
    p.add_argument('--bindingdb-comparison', type=Path, help='Frozen BindingDB comparison ZIP')
    p.add_argument('--assays-cache', type=Path, help='Existing chembl_assay_metadata.json (offline repeat)')
    p.add_argument('--self-test', action='store_true')
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    root = args.root.resolve() if args.root else root_path()
    chembl_path = (args.chembl_snapshot.resolve() if args.chembl_snapshot else root / CHEMBL_REL)
    binding_path = (args.bindingdb_comparison.resolve() if args.bindingdb_comparison else root / BDB_REL)
    p_path = root / P_REL
    out = root / OUT_REL
    zip_out = out.with_suffix('.zip')
    require(not out.exists() and not zip_out.exists(), 'Output already exists. Do not overwrite earlier datasets.')
    print('Reading frozen ChEMBL and BindingDB inputs...', flush=True)
    chembl_json = json.loads(read_zip(chembl_path, 'human_trpa1_activities.json'))
    chembl_status = json.loads(read_zip(chembl_path, 'status_start.json'))
    bdb_summary = json.loads(read_zip(binding_path, 'summary.json'))
    require(chembl_status.get('chembl_db_version') == 'ChEMBL_37', 'Input ChEMBL snapshot is not release 37.')
    require(bdb_summary.get('comparison_complete') is True and bdb_summary.get('bindingdb_release') == '202610',
            'Input BindingDB comparison is not the approved complete 202610 release.')
    require(len(chembl_json) == 6745, 'Unexpected ChEMBL snapshot size; do not silently substitute another input.')
    bind_records = csv_rows(read_zip(binding_path, 'comparison_audit.csv'))
    require(len(bind_records) == 4230, 'Unexpected BindingDB cached record count.')
    p_rows = csv_rows(p_path.read_bytes())
    p_keys = {r['inchikey'] for r in p_rows}
    require(len(p_rows) == len(p_keys) == 2383, 'P v1 changed; stop instead of modifying it.')
    chembl_candidates, chembl_counts = candidate_chembl_activities(chembl_json)
    ids = {str(r['assay_chembl_id']) for r in chembl_candidates}
    print(f'ChEMBL EC50 candidate rows: {len(chembl_candidates)}; assay descriptions to load: {len(ids)}', flush=True)
    assays, assay_bytes = get_assays(ids, args.assays_cache.resolve() if args.assays_cache else None)
    molecules, evidence, directions, flags = collect(chembl_candidates, assays, bind_records, p_keys)
    require(bool(molecules), 'Zero documented agonists found. Do not create an empty A dataset.')
    require(all(r['inchikey'] for r in molecules), 'Missing molecular identity.')
    tables = {'trpa1_known_agonists_v1.csv': csv_bytes(molecules),
              'agonist_evidence.csv': csv_bytes(evidence),
              'chembl_assay_metadata.json': assay_bytes}
    manifest = dict(version='A-v1-proposed-2026-11-08', status='SOURCE_ANNOTATION_EVIDENCE_NOT_MANUALLY_VALIDATED',
                    compounds=len(molecules), source_measurements=len(evidence),
                    overlap_inhibitor_P_v1=sum(m['in_inhibitor_P_v1'] for m in molecules),
                    chimbl_ec50_input_counts=dict(chembl_counts),
                    direction_decisions=directions, technical_flags=flags,
                    selection='EC50 = or upper bound in explicitly activating human TRPA1 assay; no potency cutoff',
                    assay_direction='based on ChEMBL/BindingDB descriptions, not independent primary-publication verification',
                    reference_sources=['ChEMBL 37 snapshot', 'BindingDB 202610 comparison'],
                    input_file_sha256={'chembl_zip': sha(chembl_path.read_bytes()),
                                       'bindingdb_zip': sha(binding_path.read_bytes()),
                                       'p_csv': sha(p_path.read_bytes())},
                    output_file_sha256={name: sha(data) for name,data in tables.items()},
                    rdkit_version=rdBase.rdkitVersion, python_version=sys.version,
                    limitations=['EC50 can represent inhibition; only explicit agonist assays are included.',
                                 'Other agonists measured without EC50, or unclear assays, may be absent.',
                                 'Evidence can be duplicated across ChEMBL and BindingDB.',
                                 'Molecules also in P v1 are preserved in both sets; neither role cancels the other.',
                                 'This script does not relabel P, change U, or classify any agonist as a negative.'])
    out.mkdir(parents=True, exist_ok=False)
    try:
        for name,data in tables.items():
            (out / name).write_bytes(data)
        (out / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        with zipfile.ZipFile(zip_out, 'x', compression=zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()):
                z.write(f, f.name)
    except BaseException:
        (out / 'INCOMPLETE.txt').write_text('Interrupted output. Do not use this incomplete A dataset.\n', encoding='utf-8')
        raise
    print('Agonist A v1 candidate set created:', len(molecules), 'unique structures; source rows:', len(evidence))
    print('Overlap with inhibitor P v1:', manifest['overlap_inhibitor_P_v1'])
    print('Output:', out, '\nZIP:', zip_out)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, FileNotFoundError, KeyError, zipfile.BadZipFile, urllib.error.URLError) as exc:
        raise SystemExit(f'ERROR: {exc}\nInputs were not modified. No P/U labels were changed.')
