#!/usr/bin/env python3
"""Build a conservative, evidence-linked human-TRPA1 positive set (P), version 1.

Run from experiments/02_classification/code/ on Windows:
    python build_trpa1_positive.py
    python build_trpa1_positive.py --self-test
Dependencies: pandas, rdkit. No network, GPU, Git operations, or new API extraction.

SCOPE: the historical 1645 ChEMBL structures + the reproduced 747 BindingDB
additions. This is curation of that frozen pool, NOT a new exhaustive search.
Inputs are the completed ChEMBL check and the reproduced BindingDB comparison
(folder or ZIP), plus the historical baseline CSV. Original files stay unchanged.

RULES (not a potency threshold):
- P requires a valid =, <, or <= IC50 > 0 in nM, linked to human WT TRPA1
  inhibition. > bounds, missing fits, agonism and desensitization are not P evidence.
- Use the saved assay audit, explicit exclusions and source-protocol review.
  An assay explicitly saying 'unknown origin' is not human evidence.
- RDKit Cleanup + known salt removal; retain stereo; reject mixtures left after
  salt removal. Do not invent absolute stereo from CXSMILES or compound names.
- For this strict first version, named stereomixtures/conflicting identities/CX annotations anywhere
  in human records for a candidate put that candidate on HOLD. This is not a
  statement that the compound is inactive or that all annotations are erroneous.
- Explicit human inhibitory inactive outcomes and measured human agonism are
  retained and put the candidate on HOLD. Inactivity after washout or in a
  mutant is not a contradiction. An IC50 > bound alone is not inactivity.
- At least one acceptable measurement is enough if there is no global HOLD.
  Never average bounds, apply the discarded BayeshERG cutoff/SD rule, or count
  database rows as independent experimental replicates.
- HOLD candidates are NOT negatives and must not be quietly put into U.

Outputs: trpa1_positive_v1.csv; evidence.csv; candidate_audit.csv; summary.json;
plus a ZIP of these. Evidence includes the original strings, inclusion decisions,
source IDs and descriptions. The 'final' set means a versioned conservative core,
not a guarantee that every underlying publication or chemical drawing is correct.

Manual review basis:
- shared/assays/assay_table.csv, Git blob 297129a53a3596eb99913fb956c1753fbe870c8e,
  commit fb0fb9a42c086f235a31c3c20c5dbc978a486f6e (97 assays; 12 excluded).
- compare_bindingdb_surechembl.py source-protocol audit (2026-10-17), six patents.
- Additional review here reads frozen assay descriptions; it does not purport to
  re-read every patent drawing or to prove direct molecular binding/selectivity.
- CXSMILES format: https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html
  The 'r' extension denotes absence of the absolute stereoconfiguration flag.
- A previously reported possible lost inequality at ChEMBL activity 14750391
  is conservatively held with imported copies, not silently corrected.
- Unspecified stereo in plain SMILES is retained as unspecified; that is NOT
  a claim that an isolated stereoisomer was experimentally characterized.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import zipfile

try:
    import pandas as pd
    from rdkit import Chem, rdBase
    from rdkit.Chem import SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError as exc:
    raise SystemExit('Install dependencies: python -m pip install pandas rdkit') from exc

VERSION = '1.0.0'
BASELINE_REL = 'experiments/01_regression/data/trpa1_primary_dataset.csv'
EXPERIMENT_REL = 'experiments/02_classification'
REF_BASE_KEYS = 'd5afeb066d87f12ec18509c68fe9cd75a7568910ed84005358de483d55e61d9e'
REF_ADDITION_KEYS = 'efe791b46151b77392d4be7e6683632658d6ba1d437253d55f76d073631eaf7c'
REF_CHEMBL_RAW = '75c1659f597cbeea1406fde7425ccab129176540824a63588923872022cc6384'
REF_BINDINGDB_RAW = 'ddc59d451b645cadd01deceb8557b6ae07867f9969f12da649942711538959dc'
# Output references are checks, NEVER used to assign labels or force a count.
REF_OUTPUT_KEYS = 'c033c5cc99bfbdfde11a75e5835024f83c87c076ead34dbfc1aad059fcf87dc3'
REF_OUTPUT_LABELS = '43f4473373e0a27994f96441afddea49cf61da9d1baeae9ab800d64c43f419cd'
REF_AUDIT_CONTENT = '7b88291822fce9ee3375da5d433a3d4a969dd309422ef4a470b57025adb77c1f'  # parsed fields, not CSV newlines
REVIEWED_ASSAYS = set('''CHEMBL1049520 CHEMBL1118226 CHEMBL1118227 CHEMBL1118233 CHEMBL1118234 CHEMBL1118499 CHEMBL1118500 CHEMBL1118501 CHEMBL1118734 CHEMBL1118735 CHEMBL1118736 CHEMBL1118737 CHEMBL1118738 CHEMBL1118739 CHEMBL1118740 CHEMBL1944172 CHEMBL1949457 CHEMBL1949463 CHEMBL2071975 CHEMBL2071976 CHEMBL2328091 CHEMBL2422134 CHEMBL2433975 CHEMBL3222891 CHEMBL3299550 CHEMBL3299551 CHEMBL3299560 CHEMBL3299561 CHEMBL3300585 CHEMBL3363074 CHEMBL3364861 CHEMBL3706106 CHEMBL3757221 CHEMBL3757748 CHEMBL3789590 CHEMBL3791636 CHEMBL3858029 CHEMBL3887573 CHEMBL3888572 CHEMBL3888736 CHEMBL4022611 CHEMBL4022625 CHEMBL4022626 CHEMBL4036813 CHEMBL4036817 CHEMBL4036824 CHEMBL4049011 CHEMBL4190621 CHEMBL4190624 CHEMBL4190628 CHEMBL4190659 CHEMBL4347272 CHEMBL4347276 CHEMBL4347295 CHEMBL4347296 CHEMBL4828020 CHEMBL4837096 CHEMBL4880071 CHEMBL4880072 CHEMBL4885651 CHEMBL4885657 CHEMBL5096638 CHEMBL5096639 CHEMBL5126028 CHEMBL5142855 CHEMBL5210282 CHEMBL5210580 CHEMBL5318377 CHEMBL5318380 CHEMBL5360087 CHEMBL5577360 CHEMBL5577372 CHEMBL5577374 CHEMBL5577382 CHEMBL5621665 CHEMBL5680872 CHEMBL5680880 CHEMBL5680929 CHEMBL5680930 CHEMBL5708945 CHEMBL5731467 CHEMBL5732578 CHEMBL5733104 CHEMBL5734964 CHEMBL5735474 CHEMBL5735475 CHEMBL5735476 CHEMBL5735477 CHEMBL5735742 CHEMBL5736402 CHEMBL5737447 CHEMBL5737489 CHEMBL5737786 CHEMBL5739043 CHEMBL5739380 CHEMBL6104158 CHEMBL6104159'''.split())
EXCLUDED_ASSAYS = {
    'CHEMBL5096638': 'agonist_desensitization',
    'CHEMBL5096639': 'agonist_desensitization',
    'CHEMBL4190624': 'chimeric_channel',
    'CHEMBL4190621': 'chimeric_channel',
    'CHEMBL4036817': 'chimeric_channel',
    'CHEMBL5732578': 'previous_audit_duplicate_evidence',
    **dict.fromkeys('CHEMBL5680930 CHEMBL1949463 CHEMBL2328091 CHEMBL4347295 '
                    'CHEMBL4347296 CHEMBL5577382'.split(),
                    'excluded_by_previous_assay_audit'),
}
# These are source-level protocol reviews, NOT per-compound structure validation.
PATENTS = {
    'US10626112': ('US10626112B2', 'Experimental Examples 1-1/1-2; Tables 40/41'),
    'US12391660': ('US12391660B2', 'FLIPR Table 3.1; patch clamp Table 3.2'),
    'US20240101561': ('US20240101561A1', 'Assay A; Table 1'),
    'US20250154139': ('US20250154139A1', 'Example 25; Table A'),
    'US20250270198': ('US20250270198A1', 'Example 6; activity table'),
    'US11701353': ('US11891403B2', 'Assay A; Table 1 examples 1-24; corrected source ID'),
}
CORRECTED_NM = [53,15,23,32,39,45,60,62,73,80,81,115,145,191,227,237,248,153,214,140,155,74,52,100]
NEGATIVE_WORDS = {'not active', 'inactive', 'not_active', 'not-active', 'no activity'}
BAD_CONTEXT = re.compile(r'\bmutant\b|\bmutation\b|\bchimera\b|\bchimeric\b|\boS5\b|desensiti[sz]', re.I)
INHIBITION = re.compile(r'\b(?:antagonist|antagonism|inhibition|inhibitory|blocking|blockade)\b', re.I)
AG_START = re.compile(r'^\s*(?:agonist activity|activation of|activation at)', re.I)
NAMED_MIXTURE = re.compile(r'\bracem(?:ic|ate|ates)|\brac\s*[-(]|\(\s*\+\s*/\s*-\s*\)|\(±\)|\bmixture\b', re.I)
# Previous project audit reported a possible lost '>' in Rooney's compound 22.
# This build does NOT silently change '=' to '>' or claim a fresh PDF check.
# Quarantine this recorded value and its imported copies pending primary review.
PENDING_VALUE_KEYS = {
    ('HHCCDKWJHHNZND-UHFFFAOYSA-N', 'CHEMBL3299551'):
        'previous_audit_possible_lost_relation_14750391_Rooney2014',
}

REMOVER = SaltRemover.SaltRemover()


def txt(v) -> str:
    return '' if v is None else str(v).strip()


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parsed_audit_digest(rows) -> str:
    # Excludes derived 'standard_smiles'; all source strings/decisions used here
    # remain protected. Different CSV newlines do not alter this digest.
    fields = ['rsid','monomer_id','ligand_name','original_smiles','target_name',
        'organism','accession','target_chains','ic50_raw_nM','ec50_raw_nM',
        'patent_raw','article_doi','pmid','chembl_id','curation_source',
        'bindingdb_entry_doi','assay_ids','assay_description','standard_inchikey','novelty_status']
    canonical = [{k: r[k] for k in fields} for r in sorted(rows,key=lambda r:r['rsid'])]
    return digest(json.dumps(canonical, sort_keys=True, ensure_ascii=False,
                             separators=(',',':')).encode('utf-8'))


def key_digest(keys) -> str:
    return digest(('\n'.join(sorted(set(keys))) + '\n').encode('utf-8'))


def normalized(s: str) -> str:
    return ' '.join(s.lower().split())


def read_member(source: Path, name: str) -> bytes:
    if source.is_dir():
        return (source / name).read_bytes()
    with zipfile.ZipFile(source) as z:
        hits = [p for p in z.namelist() if p == name or p.endswith('/' + name)]
        if len(hits) != 1:
            raise ValueError(f'{source}: expected one {name}, got {len(hits)}')
        return z.read(hits[0])


def csv_rows(raw: bytes) -> list[dict]:
    return pd.read_csv(io.BytesIO(raw), dtype=str, keep_default_na=False).to_dict('records')


def numeric(s) -> float | None:
    try:
        n = float(s)
        return n if math.isfinite(n) and n > 0 else None
    except (TypeError, ValueError):
        return None


def parse_ic50(s: str) -> tuple[str, float | None]:
    m = re.fullmatch(r'\s*(<=|>=|<|>|=|≤|≥)?\s*(\d+(?:\.\d*)?|\.\d+)([eE][+-]?\d+)?\s*', s)
    if not m:
        return '', None
    return {'≤':'<=', '≥':'>='}.get(m[1], m[1]) or '=', numeric(m[2] + (m[3] or ''))


@lru_cache(maxsize=30000)
def structure(smiles: str) -> dict:
    """Historical key for tracing; normalized structure only when unambiguous."""
    out = dict(historical_key='', smiles='', inchikey='', connectivity_key='',
               structure_issue='', unspecified_tetra_centers=0)
    if not smiles:
        out['structure_issue'] = 'missing_structure'
        return out
    try:
        params = Chem.SmilesParserParams()
        params.parseName = False
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles, params)
        if mol is None:
            out['structure_issue'] = 'invalid_structure'
            return out
        old = REMOVER.StripMol(Chem.Mol(mol), dontRemoveEverything=True)
        fragments = Chem.GetMolFrags(old, asMols=True, sanitizeFrags=True)
        old = max(fragments, key=lambda m: m.GetNumHeavyAtoms())
        with rdBase.BlockLogs():
            out['historical_key'] = Chem.MolToInchiKey(old)
        # Do NOT silently erase extensions, including |r|. Keep raw text in ledger.
        if '|' in smiles:
            out['structure_issue'] = ('relative_stereo_annotation' if re.search(r'\|r(?:[,:|])', smiles)
                                      else 'other_CXSMILES_annotation')
            return out
        mol = rdMolStandardize.Cleanup(mol)
        mol = REMOVER.StripMol(mol, dontRemoveEverything=True)
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(fragments) != 1:
            out['structure_issue'] = 'unresolved_multicomponent_structure'
            return out
        mol = fragments[0]
        if any(a.GetAtomicNum() == 0 or a.HasQuery() for a in mol.GetAtoms()):
            out['structure_issue'] = 'query_or_dummy_atoms'
            return out
        out['smiles'] = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        with rdBase.BlockLogs():
            out['inchikey'] = Chem.MolToInchiKey(mol)
        if not out['inchikey']:
            out['structure_issue'] = 'missing_inchikey'
            return out
        out['connectivity_key'] = out['inchikey'].split('-')[0]
        out['unspecified_tetra_centers'] = sum(
            label == '?' for _, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                                                  useLegacyImplementation=False))
        return out
    except (ValueError, RuntimeError, IndexError) as exc:
        out['structure_issue'] = 'structure_error:' + type(exc).__name__
        return out


def human_context(description: str, mutation: str = '') -> str:
    if mutation:
        return 'variant_channel'
    if BAD_CONTEXT.search(description):
        return 'mutant_chimera_or_desensitization'
    if 'unknown origin' in description.lower():
        return 'channel_species_explicitly_unknown'
    return ''


def chembl_assay_status(r: dict) -> str:
    aid, desc = txt(r.get('assay_chembl_id')), txt(r.get('assay_description'))
    if (txt(r.get('target_chembl_id')) != 'CHEMBL6007' or
        txt(r.get('target_organism')).lower() != 'homo sapiens' or
        txt(r.get('target_tax_id')) not in {'9606', '9606.0'}):
        return 'not_human_TRPA1'
    issue = human_context(desc, txt(r.get('assay_variant_mutation')))
    if issue:
        return issue
    if aid in EXCLUDED_ASSAYS:
        return EXCLUDED_ASSAYS[aid]
    if aid not in REVIEWED_ASSAYS:
        return 'assay_not_in_reviewed_97'
    return ''


def clear_human_inhibition(description: str) -> bool:
    return (not human_context(description) and 'human' in description.lower()
            and 'trpa1' in description.lower() and bool(INHIBITION.search(description))
            and not AG_START.search(description))


def biological_flags(r: dict) -> set[str]:
    """Only use actual human WT assay meaning, not generic target labels."""
    flags = set()
    if (txt(r.get('target_organism')).lower() != 'homo sapiens' or
        txt(r.get('target_chembl_id')) != 'CHEMBL6007' or
        human_context(txt(r.get('assay_description')), txt(r.get('assay_variant_mutation')))):
        return flags
    desc = txt(r.get('assay_description'))
    words = {normalized(txt(r.get(k))) for k in ('standard_text_value','text_value','activity_comment')}
    if (words & NEGATIVE_WORDS and clear_human_inhibition(desc) and
        not re.search(r'wash.?out|wash.?off', desc, re.I) and
        txt(r.get('assay_chembl_id')) not in EXCLUDED_ASSAYS):
        flags.add('explicit_inhibitory_inactive_outcome')
    if (AG_START.search(desc) and 'human' in desc.lower() and
        txt(r.get('standard_type')).upper() in {'EC50','PEC50'} and
        txt(r.get('standard_relation')) == '=' and numeric(r.get('standard_value'))):
        flags.add('also_measured_human_agonism')
    return flags


def binding_protocol(r: dict, by_aid: dict, by_description: dict) -> tuple[str,str,str,str]:
    """Return failure reason, support URL, locator, linked ChEMBL assay."""
    desc = txt(r.get('assay_description'))
    issue = human_context(desc + ' ' + txt(r.get('target_name')))
    if issue:
        return issue, '', '', ''
    aids = sorted(set(re.findall(r'CHEMBL\d+', desc)))
    if not aids and '\t' in desc:
        aids = by_description.get(normalized(desc.split('\t', 1)[1]), [])
    known = [aid for aid in aids if aid in by_aid]
    if known:
        accepted = [aid for aid in known if not chembl_assay_status(by_aid[aid])]
        if not accepted:
            return 'linked_chembl_assay_not_eligible', '', ';'.join(known), ';'.join(known)
        aid = accepted[0]
        return '', f'https://www.ebi.ac.uk/chembl/explore/assay/{aid}', 'saved ChEMBL assay audit', aid
    patent = txt(r.get('patent_raw'))
    if patent in PATENTS:
        support, locator = PATENTS[patent]
        if patent == 'US11701353':
            m = re.search(r'US11891403,\s*Example\s+(\d+)', txt(r.get('ligand_name')))
            relation, value = parse_ic50(txt(r.get('ic50_raw_nM')))
            if not (m and 1 <= int(m[1]) <= 24 and relation == '=' and value is not None and
                    math.isclose(value, CORRECTED_NM[int(m[1])-1], rel_tol=0, abs_tol=1e-9)):
                return 'unresolved_patent_identifier_conflict', '', locator, ''
        return '', f'https://patents.google.com/patent/{support}/en', locator, ''
    # US9186360: human CHO/AITC inhibitory protocol is fully present in the cache.
    # This is explicitly a DATABASE-PROTOCOL review, not a new primary PDF review.
    if patent == 'US9186360' and clear_human_inhibition(desc):
        return '', 'https://patents.google.com/patent/US9186360B2/en', 'BindingDB saved assay 7368_1; database protocol reviewed', ''
    return 'unreviewed_inhibition_protocol', '', '', ''


def source_ids(r: dict) -> set[str]:
    return set(re.findall(r'CHEMBL\d+', txt(r.get('chembl_id'))))


def joined(values) -> str:
    return ';'.join(sorted({txt(v) for v in values if txt(v)}))


def find_root() -> Path:
    for start in (Path(__file__).resolve().parent, Path.cwd()):
        for parent in (start, *start.parents):
            if (parent / BASELINE_REL).is_file():
                return parent
    raise ValueError('Cannot find repository baseline. Use --repo-root C:\\trpa1_data')


def discover(results: Path, kind: str) -> Path:
    candidates = []
    for p in results.iterdir() if results.is_dir() else []:
        if not (p.is_dir() or p.suffix.lower() == '.zip'):
            continue
        if kind == 'chembl' and not ('chembl_live_check_' in p.name.lower() or 'trpa1_live_check_' in p.name.lower()):
            continue
        if kind == 'comparison' and not p.name.startswith('bindingdb_surechembl_comparison_'):
            continue
        try:
            meta = json.loads(read_member(p, 'summary.json'))
            good = (meta.get('complete') is True and meta.get('human_target') == 'CHEMBL6007') if kind == 'chembl' else meta.get('comparison_complete') is True
            if good:
                candidates.append(p)
        except (OSError, ValueError, KeyError, zipfile.BadZipFile):
            pass
    if not candidates:
        option = '--chembl-run' if kind == 'chembl' else '--comparison-run'
        raise ValueError(f'No completed {kind} run in {results}. Specify {option} path_to_folder_or_zip')
    def timestamp_key(p):
        m=re.search(r'(20\d{6}T\d{6}(?:_\d+)?Z)',p.name)
        return (m[1] if m else p.stem, p.is_dir())
    return max(candidates, key=timestamp_key)


def load_inputs(root: Path, chembl: Path, comparison: Path, baseline: Path) -> tuple:
    base_raw = baseline.read_bytes()
    base = csv_rows(base_raw)
    base_keys = [r['inchikey'] for r in base]
    if len(base_keys) != len(set(base_keys)) or key_digest(base_keys) != REF_BASE_KEYS:
        raise ValueError('Baseline identities differ from the reviewed 1645 pool. No silent relabelling.')
    for r in base:
        if structure(r['std_smiles'])['historical_key'] != r['inchikey']:
            raise ValueError('RDKit baseline identity mismatch: ' + r['inchikey'])
    cj = json.loads(read_member(chembl, 'summary.json'))
    bj = json.loads(read_member(comparison, 'summary.json'))
    cr = read_member(chembl, 'human_trpa1_activities.json')
    br = read_member(comparison, 'comparison_audit.csv')
    nr = read_member(comparison, 'new_bindingdb_molecules.csv')
    if cj.get('complete') is not True or cj.get('human_target') != 'CHEMBL6007':
        raise ValueError('ChEMBL download is incomplete or is not human TRPA1.')
    if bj.get('comparison_complete') is not True or bj.get('reference_key_set_comparison') != 'MATCH':
        raise ValueError('BindingDB comparison is not the completed reproduced result.')
    if digest(cr) != REF_CHEMBL_RAW or bj.get('bindingdb_raw_sha256') != REF_BINDINGDB_RAW:
        raise ValueError('Different source snapshots. This v1 recipe is frozen; review new inputs separately.')
    if bj.get('baseline_key_set_sha256') != REF_BASE_KEYS:
        raise ValueError('BindingDB comparison used a different baseline.')
    c, b, additions = json.loads(cr), csv_rows(br), csv_rows(nr)
    if REF_AUDIT_CONTENT and parsed_audit_digest(b) != REF_AUDIT_CONTENT:
        raise ValueError('BindingDB source record content changed from reviewed input; do not force a match.')
    new_keys = [r['inchikey'] for r in additions]
    if (len(new_keys) != len(set(new_keys)) or key_digest(new_keys) != REF_ADDITION_KEYS or
        set(base_keys) & set(new_keys)):
        raise ValueError('Addition identity list differs from the reproduced 747 or overlaps baseline.')
    if len(c) != cj['total_human_activity_records'] or len(b) != bj['counts']['bindingdb_TRPA1_related_records']:
        raise ValueError('Input row count mismatch.')
    seen_ids = set()
    for r in b:
        if r['rsid'] in seen_ids:
            raise ValueError('Duplicate BindingDB record ID; review input: ' + r['rsid'])
        seen_ids.add(r['rsid'])
        got = structure(r['original_smiles'])['historical_key']
        if got != r['standard_inchikey']:
            raise ValueError('BindingDB identity differs from reproduced audit: ' + r['rsid'])
    if {r['standard_inchikey'] for r in b if r['novelty_status']=='new_to_primary_conservative'} != set(new_keys):
        raise ValueError('BindingDB molecule list and record audit disagree.')
    metadata = dict(baseline_path=str(baseline), baseline_file_sha256=digest(base_raw),
                    baseline_key_sha256=key_digest(base_keys), chembl_run=str(chembl),
                    chembl_raw_sha256=digest(cr), comparison_run=str(comparison),
                    comparison_audit_sha256=digest(br), additions_sha256=digest(nr),
                    bindingdb_original_raw_sha256=REF_BINDINGDB_RAW,
                    bindingdb_parsed_content_sha256=parsed_audit_digest(b),
                    chembl_version=cj.get('chembl_version'), bindingdb_release=bj.get('bindingdb_release'))
    return base, additions, c, b, metadata


def build(base: list[dict], additions: list[dict], c: list[dict], b: list[dict]) -> tuple:
    base_keys = {r['inchikey'] for r in base}
    additions_keys = {r['inchikey'] for r in additions}
    candidates = base_keys | additions_keys
    evidence, flags, ids_to_keys = [], defaultdict(set), defaultdict(set)
    monomer_keys = defaultdict(set)
    by_aid, by_description = {}, defaultdict(list)
    assay_fields = ('assay_description','assay_variant_mutation','target_organism',
                    'target_chembl_id','target_tax_id')
    for r in c:
        aid = txt(r.get('assay_chembl_id'))
        if aid in by_aid and any(txt(r.get(f)) != txt(by_aid[aid].get(f)) for f in assay_fields):
            raise ValueError('Inconsistent metadata within assay ' + aid)
        by_aid[aid] = r
        s = normalized(txt(r.get('assay_description')))
        if aid not in by_description[s]:
            by_description[s].append(aid)
        key = structure(txt(r.get('canonical_smiles')))['historical_key']
        if key:
            ids_to_keys[txt(r.get('molecule_chembl_id'))].add(key)
    for r in b:
        key = r['standard_inchikey']
        if key and r['organism'].lower() in {'human','homo sapiens'} and r['accession']=='O75762':
            for cid in source_ids(r):
                ids_to_keys[cid].add(key)
            if txt(r.get('monomer_id')):
                monomer_keys[r['monomer_id']].add(key)
    for cid, keys in ids_to_keys.items():
        if cid and len(keys) > 1:
            for key in keys & candidates:
                flags[key].add('cross_database_structure_identifier_conflict')

    for monomer, keys in monomer_keys.items():
        if len(keys) > 1:
            for key in keys & candidates:
                flags[key].add('BindingDB_monomer_identity_conflict')

    # Retain ALL source records of candidate identities, even ones rejected below.
    # This prevents a bad/contradictory row disappearing behind an accepted duplicate.
    for r in c:
        st = structure(txt(r.get('canonical_smiles')))
        key = st['historical_key']
        if key not in candidates:
            continue
        flags[key].update(biological_flags(r))
        relation, value = txt(r.get('standard_relation')), numeric(r.get('standard_value'))
        reason = chembl_assay_status(r)
        reason = reason or ('not_numeric_inhibitory_IC50' if txt(r.get('standard_type')).upper() != 'IC50' else '')
        reason = reason or ('not_positive_IC50_relation' if relation not in {'=','<','<='} or value is None else '')
        reason = reason or ('non_nM_units' if txt(r.get('standard_units')) != 'nM' else '')
        reason = reason or ('database_validity_flag' if txt(r.get('data_validity_comment')) not in {'','Manually validated'} else '')
        reason = reason or PENDING_VALUE_KEYS.get((key, txt(r.get('assay_chembl_id'))), '')
        reason = reason or st['structure_issue']
        if st['structure_issue']:
            flags[key].add(st['structure_issue'])
        evidence.append(dict(source='ChEMBL', record_id=txt(r['activity_id']), candidate_key=key,
            original_smiles=txt(r.get('canonical_smiles')), standard_smiles=st['smiles'],
            standardized_key=st['inchikey'], molecule_id=txt(r.get('molecule_chembl_id')),
            assay_id=txt(r.get('assay_chembl_id')), document_id=txt(r.get('document_chembl_id')),
            source_url='https://www.ebi.ac.uk/chembl/explore/activity/' + txt(r['activity_id']),
            source_locator='frozen ChEMBL record + saved 97-assay review',
            supporting_document=txt(r.get('document_chembl_id')),
            assay_description=txt(r.get('assay_description')), endpoint=txt(r.get('standard_type')),
            relation=relation, value_original=txt(r.get('standard_value')),
            units=txt(r.get('standard_units')), text_result=joined(r.get(k) for k in ('standard_text_value','activity_comment','text_value')),
            validity_flag=txt(r.get('data_validity_comment')), potential_duplicate=txt(r.get('potential_duplicate')),
            ligand_name=txt(r.get('molecule_pref_name')), source_organism=txt(r.get('target_organism')),
            source_target=txt(r.get('target_chembl_id')), source_accession='',
            source_chembl_ids=txt(r.get('molecule_chembl_id')), source_article_doi='',
            source_patent='', source_bindingdb_inchikey='',
            positive_evidence=not bool(reason), record_exclusion_reason=reason,
            structure_issue=st['structure_issue'], unspecified_tetra_centers=st['unspecified_tetra_centers']))
    # Index ChEMBL flags by exact identity + assay, so an imported copy cannot
    # bypass the original validity warning just because BindingDB omitted it.
    quality = defaultdict(set)
    for r in c:
        key = structure(txt(r.get('canonical_smiles')))['historical_key']
        if key in candidates and txt(r.get('standard_type')).upper() == 'IC50':
            q = txt(r.get('data_validity_comment'))
            if q and q != 'Manually validated':
                quality[(key, txt(r.get('assay_chembl_id')))].add(q)
    for r in b:
        key = r['standard_inchikey']
        if key not in candidates:
            continue
        st = structure(r['original_smiles'])
        human = (r['organism'].lower() in {'human','homo sapiens'} and r['accession']=='O75762'
                 and r['target_chains']=='1')
        if human and st['structure_issue']:
            flags[key].add(st['structure_issue'])
        name_mixture = human and bool(NAMED_MIXTURE.search(txt(r.get('ligand_name'))))
        if name_mixture:
            flags[key].add('named_racemate_or_mixture')
        erel, eval_ = parse_ic50(txt(r.get('ec50_raw_nM')))
        edesc = txt(r.get('assay_description')).split('\t')[-1]
        if (human and not human_context(edesc) and 'human' in edesc.lower()
            and AG_START.search(edesc) and erel=='=' and eval_ is not None):
            flags[key].add('also_measured_human_agonism')
        relation, value = parse_ic50(r['ic50_raw_nM'])
        issue, url, locator, linked = binding_protocol(r, by_aid, by_description)
        reason = ('not_explicit_human_single_chain' if not human else '')
        reason = reason or ('not_positive_IC50_relation' if relation not in {'=','<','<='} or value is None else '')
        reason = reason or issue
        reason = reason or ('original_ChEMBL_validity_flag' if quality.get((key, linked)) else '')
        reason = reason or PENDING_VALUE_KEYS.get((key, linked), '')
        reason = reason or st['structure_issue']
        reason = reason or ('named_racemate_or_mixture' if name_mixture else '')
        if r['novelty_status'].startswith('review_'):
            reason = reason or r['novelty_status']
            flags[key].add(r['novelty_status'])
        evidence.append(dict(source='BindingDB', record_id=r['rsid'], candidate_key=key,
            original_smiles=r['original_smiles'], standard_smiles=st['smiles'], standardized_key=st['inchikey'],
            molecule_id=r['monomer_id'], assay_id=r['assay_ids'], document_id=r['patent_raw'] or r['article_doi'],
            source_url=url or (f'https://doi.org/{r["bindingdb_entry_doi"]}' if r['bindingdb_entry_doi'] else ''),
            supporting_document=(url.split('/patent/')[1].split('/')[0] if '/patent/' in url
                                 else (r['patent_raw'] or r['article_doi'])),
            source_locator=locator, assay_description=r['assay_description'], endpoint='IC50', relation=relation,
            value_original=r['ic50_raw_nM'], units='nM', text_result='',
            validity_flag=joined(quality.get((key,linked),set())), potential_duplicate=r['curation_source'],
            ligand_name=r['ligand_name'], source_organism=r['organism'], source_target=r['target_name'],
            source_accession=r['accession'], source_chembl_ids=r['chembl_id'],
            source_article_doi=r['article_doi'], source_patent=r['patent_raw'],
            source_bindingdb_inchikey=r['bindingdb_inchikey'],
            positive_evidence=not bool(reason), record_exclusion_reason=reason,
            structure_issue=st['structure_issue'], unspecified_tetra_centers=st['unspecified_tetra_centers']))
    by_key = defaultdict(list)
    for r in evidence:
        by_key[r['candidate_key']].append(r)
    positives, audit = [], []
    for key in sorted(candidates):
        rr = by_key[key]
        good = [r for r in rr if r['positive_evidence']]
        reasons = set(flags[key])
        if not good:
            reasons.add('no_eligible_positive_measurement')
        mapped = {r['standardized_key'] for r in good}
        if len(mapped) > 1:
            reasons.add('inconsistent_normalized_identity')
        status = 'HOLD' if reasons else 'P'
        origin = 'historical_ChEMBL_1645' if key in base_keys else 'BindingDB_additions_747'
        audit.append(dict(candidate_key=key, origin_pool=origin, status=status,
            reasons=joined(reasons), record_exclusion_reasons=joined(r['record_exclusion_reason'] for r in rr),
            supporting_database_rows=len(good), all_database_rows=len(rr),
            chembl_ids=joined(r['molecule_id'] for r in rr if r['source']=='ChEMBL'),
            bindingdb_ids=joined(r['molecule_id'] for r in rr if r['source']=='BindingDB'),
            identity_conflict_details=joined(
                'ChEMBL:'+cid+' -> '+','.join(sorted(keys))
                for cid,keys in ids_to_keys.items() if key in keys and len(keys)>1),
            exclude_from_U=True))
        for r in rr:
            r['molecule_status'] = status
            r['used_for_P'] = status == 'P' and r['positive_evidence']
        if status != 'P':
            continue
        exemplar = min(good, key=lambda r: (r['source']!='ChEMBL', r['record_id']))
        exact = [parse_ic50(r['value_original'])[1] for r in good if r['relation']=='=']
        exact = [x for x in exact if x is not None]
        positives.append(dict(inchikey=exemplar['standardized_key'], smiles=exemplar['standard_smiles'],
            connectivity_key=exemplar['standardized_key'].split('-')[0], label=1, historical_key=key,
            origin_pool=origin, supporting_sources=joined(r['source'] for r in good),
            chembl_ids=joined(r['molecule_id'] for r in good if r['source']=='ChEMBL'),
            bindingdb_ids=joined(r['molecule_id'] for r in good if r['source']=='BindingDB'),
            evidence_ids=joined(r['source']+':'+r['record_id'] for r in good),
            documents=joined(r['supporting_document'] for r in good), n_supporting_database_rows=len(good),
            has_exact_ic50=bool(exact), has_upper_bound=any(r['relation'] in {'<','<='} for r in good),
            min_exact_ic50_nM=min(exact) if exact else '', max_exact_ic50_nM=max(exact) if exact else '',
            unspecified_tetra_centers=exemplar['unspecified_tetra_centers']))
    positives.sort(key=lambda r:r['inchikey'])
    evidence.sort(key=lambda r:(r['candidate_key'],r['source'],r['record_id']))
    if len({r['inchikey'] for r in positives}) != len(positives):
        raise ValueError('Cleanup merged candidates. Review normalization before finalizing P.')
    if len(audit) != len(candidates) or not positives:
        raise ValueError('Incomplete candidate audit or empty P.')
    if not all(r['label']==1 and r['n_supporting_database_rows']>0 and
               not flags[r['historical_key']] for r in positives):
        raise ValueError('Internal positive-set validation failed.')
    if not all(r['relation'] in {'=','<','<='} and r['units']=='nM' and
               not r['record_exclusion_reason'] for r in evidence if r['used_for_P']):
        raise ValueError('An unsupported measurement entered the positive set.')
    return positives, evidence, audit


def hold_category(reasons: str) -> str:
    rs=set(reasons.split(';'))
    if rs & {'relative_stereo_annotation','other_CXSMILES_annotation',
             'named_racemate_or_mixture','unresolved_multicomponent_structure'}:
        return 'stereochemistry_or_mixture'
    if rs & {'cross_database_structure_identifier_conflict','BindingDB_monomer_identity_conflict',
             'inconsistent_normalized_identity'}:
        return 'structure_identifier_conflict'
    if rs & {'also_measured_human_agonism','explicit_inhibitory_inactive_outcome'}:
        return 'biological_conflict_or_mixed_action'
    return 'no_eligible_evidence_or_pending_record_check'


def export_csv(path: Path, rows: list[dict]):
    if not rows:
        raise ValueError('Refusing an unstructured empty output: ' + path.name)
    pd.DataFrame(rows).to_csv(path, index=False, encoding='utf-8', lineterminator='\n', float_format='%.12g')


def self_test():
    assert parse_ic50('>30') == ('>',30.0)
    assert parse_ic50('< 1000') == ('<',1000.0)
    assert parse_ic50('1.2e3') == ('=',1200.0)
    assert parse_ic50('0')[1] is None and parse_ic50('NaN')[1] is None
    assert parse_ic50('Not Determined')[1] is None
    assert structure('C[C@H](O)F |r|')['structure_issue']=='relative_stereo_annotation'
    assert not structure('C[C@H](O)F')['structure_issue']
    assert structure('C[C@H](O)F')['inchikey'] != structure('C[C@@H](O)F')['inchikey']
    assert structure('CCN.Cl')['inchikey']==structure('CCN')['inchikey']
    assert structure('CC(O)F')['unspecified_tetra_centers'] == 1
    assert NAMED_MIXTURE.search('rac-(1R,2R)-compound')
    assert NAMED_MIXTURE.search('Racemic compound')
    assert not NAMED_MIXTURE.search('tracer compound')
    assert structure('CCN.CCCO')['structure_issue']=='unresolved_multicomponent_structure'
    assert human_context('human TRPA1 mutant')
    assert human_context('TRPA1 unknown origin')
    assert not human_context('Human TRPA1, AITC inhibition')
    t=dict(target_chembl_id='CHEMBL6007',target_organism='Homo sapiens',target_tax_id='9606',
           assay_description='Antagonist activity at human TRPA1',standard_text_value='Not active')
    assert 'explicit_inhibitory_inactive_outcome' in biological_flags(t)
    assert not biological_flags({**t,'assay_description':'Agonist activity at human TRPA1'})
    assert not biological_flags({**t,'assay_description':'Antagonist activity at human TRPA1 after compound washout'})
    assert not biological_flags({**t,'assay_description':'Antagonist activity at human TRPA1 mutant'})
    assert not biological_flags({**t,'target_organism':'Rattus norvegicus'})
    assert not biological_flags({**t,'standard_text_value':'Not Determined'})
    assert 'also_measured_human_agonism' in biological_flags({**t,'standard_text_value':'',
        'assay_description':'Agonist activity at human TRPA1','standard_type':'EC50',
        'standard_relation':'=','standard_value':'1000'})
    # Synthetic end-to-end checks use no remote data.
    def cr(record_id, smiles='CCN', **kw):
        return dict(activity_id=record_id, canonical_smiles=smiles,
            molecule_chembl_id='CHEMBL999000001', assay_chembl_id='CHEMBL1118499',
            target_chembl_id='CHEMBL6007', target_organism='Homo sapiens', target_tax_id='9606',
            assay_description='Antagonist activity at human TRPA1',
            standard_type='IC50',standard_relation='=',standard_value='100',
            standard_units='nM',**kw)
    key=structure('CCN')['historical_key']
    base=[{'inchikey':key}]
    pos,ev,au=build(base,[],[cr('1')],[])
    assert len(pos)==1 and au[0]['status']=='P'
    # A > bound never creates a P label, nor does it contradict an exact positive.
    r=cr('2');r['standard_relation']='>'
    # Empty positive sets deliberately fail the builder's invariant. Test alongside a control.
    control=cr('3','CCCO');control['molecule_chembl_id']='CHEMBL999000002'
    control_key=structure('CCCO')['historical_key']
    pos,_,au=build(base+[{'inchikey':control_key}],[],[r,control],[])
    assert next(x['status'] for x in au if x['candidate_key']==key)=='HOLD'
    pos,_,au=build(base,[],[cr('1'),r],[]);assert len(pos)==1
    low=cr('1');low['standard_value']='1000000'
    pos,_,_=build(base,[],[low],[]);assert len(pos)==1  # No hidden potency threshold.
    bound=cr('1');bound['standard_relation']='<'
    pos,_,_=build(base,[],[bound],[]);assert pos[0]['has_upper_bound'] and not pos[0]['has_exact_ic50']
    neg=cr('2');neg.update(standard_type='Activity',standard_value='',standard_text_value='Not active')
    pos,_,au=build(base+[{'inchikey':control_key}],[],[cr('1'),neg,control],[])
    assert next(x['status'] for x in au if x['candidate_key']==key)=='HOLD'
    mutant=cr('2');mutant.update(assay_chembl_id='CHEMBL_TEST_MUTANT',
                               assay_description='Inhibition of human TRPA1 mutant')
    pos,_,_=build(base,[],[cr('1'),mutant],[]);assert len(pos)==1
    # One bad BindingDB row must not be hidden by a good ChEMBL row.
    def br(smiles='CCN', **kw):
        r={k:'' for k in ('rsid','monomer_id','ligand_name','original_smiles','bindingdb_inchikey',
            'target_name','organism','accession','target_chains','ic50_raw_nM','ec50_raw_nM',
            'patent_raw','article_doi','pmid','chembl_id','curation_source','bindingdb_entry_doi',
            'assay_ids','assay_description','standard_inchikey','novelty_status')}
        r.update(rsid='B1',monomer_id='M1',original_smiles=smiles,
                 target_name='TRPA1',organism='Human',accession='O75762',target_chains='1',
                 ic50_raw_nM='100',patent_raw='US9186360',
                 assay_description='Human TRPA1 inhibition in CHO cells',
                 standard_inchikey=structure(smiles)['historical_key'],novelty_status='already_in_primary_by_full_key')
        r.update(kw);return r
    pos,_,au=build(base+[{'inchikey':control_key}],[],[cr('1'),control],
                   [br(ligand_name='Racemic compound')])
    assert next(x['status'] for x in au if x['candidate_key']==key)=='HOLD'
    pos,_,au=build(base+[{'inchikey':control_key}],[],[cr('1'),control],
                   [br('CCO',chembl_id='CHEMBL999000001')])
    assert next(x['status'] for x in au if x['candidate_key']==key)=='HOLD'
    # Imported copy cannot bypass its source's quality flag.
    flagged=cr('1',data_validity_comment='Outside typical range')
    copyrow=br(assay_description='ChEMBL_x (CHEMBL1118499)\tAntagonist activity at human TRPA1')
    pos,ev,au=build(base+[{'inchikey':control_key}],[],[flagged,control],[copyrow])
    assert next(x['status'] for x in au if x['candidate_key']==key)=='HOLD'
    assert any(x['record_exclusion_reason']=='original_ChEMBL_validity_flag' for x in ev)
    print('SELF-TEST: PASS')


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--repo-root',type=Path)
    p.add_argument('--chembl-run',type=Path)
    p.add_argument('--comparison-run',type=Path)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--output-dir',type=Path,help='New, not-yet-existing directory')
    p.add_argument('--self-test',action='store_true')
    args=p.parse_args()
    if args.self_test:
        self_test(); return
    root=args.repo_root.resolve() if args.repo_root else find_root()
    results=root/EXPERIMENT_REL/'results'
    chembl=args.chembl_run or discover(results,'chembl')
    comparison=args.comparison_run or discover(results,'comparison')
    baseline=args.baseline or root/BASELINE_REL
    print(f'Build {VERSION}; RDKit {rdBase.rdkitVersion}; OFFLINE')
    print(f'ChEMBL:    {chembl}\nBindingDB: {comparison}\nBaseline:  {baseline}')
    base, additions, c, b, inputs=load_inputs(root,chembl,comparison,baseline)
    print('Inputs verified. Reviewing 2392 candidate identities and their source records...')
    positives,evidence,audit=build(base,additions,c,b)
    out=args.output_dir or results/('positive_set_v1_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ'))
    if out.exists():
        raise ValueError(f'Refusing to overwrite {out}')
    out.mkdir(parents=True)
    try:
        export_csv(out/'trpa1_positive_v1.csv',positives)
        export_csv(out/'evidence.csv',evidence)
        export_csv(out/'candidate_audit.csv',audit)
        keys=key_digest(r['inchikey'] for r in positives)
        labels=digest(('\n'.join(r['candidate_key']+'\t'+r['status']+'\t'+r['reasons'] for r in audit)+'\n').encode())
        match='NOT_SET' if not REF_OUTPUT_KEYS else ('MATCH' if keys==REF_OUTPUT_KEYS and labels==REF_OUTPUT_LABELS else 'MISMATCH')
        counts=dict(candidate_molecules=len(audit),positive_molecules=len(positives),
                    from_historical_baseline=sum(r['origin_pool']=='historical_ChEMBL_1645' for r in positives),
                    from_bindingdb_additions=sum(r['origin_pool']=='BindingDB_additions_747' for r in positives),
                    held_molecules=sum(r['status']=='HOLD' for r in audit),
                    positive_with_exact_ic50=sum(r['has_exact_ic50'] for r in positives),
                    positive_upper_bound_only=sum(not r['has_exact_ic50'] for r in positives),
                    supporting_database_rows=sum(r['used_for_P'] for r in evidence),
                    audit_database_rows=len(evidence))
        reason_counts=Counter(x for r in audit if r['status']=='HOLD' for x in r['reasons'].split(';') if x)
        summary=dict(version=VERSION,created_utc=datetime.now(timezone.utc).isoformat(),
            complete=match!='MISMATCH',reference_result=match,rdkit=rdBase.rdkitVersion,pandas=pd.__version__,
            python=sys.version,script_sha256=digest(Path(__file__).read_bytes()),inputs=inputs,
            counts=counts,hold_reasons_nonexclusive=dict(sorted(reason_counts.items())),
            hold_groups_exclusive=dict(sorted(Counter(hold_category(r['reasons'])
                for r in audit if r['status']=='HOLD').items())),
            rules_are_conservative_v1_not_a_universal_standard=True,
            pending_record_checks={'ChEMBL:14750391':
                'Possible lost inequality reported previously; primary PDF not rechecked in this build. '
                'Record and matched BindingDB copy are quarantined, not corrected.'},
            positive_key_sha256=keys,decision_sha256=labels,
            output_sha256={n:digest((out/n).read_bytes()) for n in ('trpa1_positive_v1.csv','evidence.csv','candidate_audit.csv')},
            policy=['Scope is frozen historical 1645 + reproduced 747, not an exhaustive fresh search.',
                'No potency cutoff. Positive IC50 relation =, <, <= only; bounds are not averaged.',
                'Review uses saved human assay descriptions and previous source-level audits.',
                'Strict core holds CXSMILES, explicitly named stereomixtures, and identity conflicts.',
                'HOLD is not inactive and not automatically U. Reasons are not mutually exclusive.',
                'No median activity or independent-experiment count is inferred from imported source rows.',
                'Use connectivity_key for grouping; do not split same molecule/series across train/test.',
                'No SureChEMBL mentions, random decoys, or new model results are included.'])
        (out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        zip_path=out.with_suffix('.zip')
        with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()): z.write(f,f.name)
        print('\nRESULTS')
        for k,v in counts.items():print(f'  {k}: {v}')
        print(f'  Reference result: {match}\n  P identity SHA256: {keys}\nResults: {out}\nZIP: {zip_path}')
        if match=='MISMATCH':
            raise SystemExit('Reference mismatch: keep the outputs for review; do not use them as reproduced P.')
    except (Exception,KeyboardInterrupt) as exc:
        (out/'INCOMPLETE.json').write_text(json.dumps({'complete':False,'error':str(exc)},ensure_ascii=False,indent=2),encoding='utf-8')
        raise


if __name__=='__main__':
    try:
        main()
    except (OSError,ValueError,zipfile.BadZipFile) as exc:
        raise SystemExit('ERROR: '+str(exc))
