#!/usr/bin/env python3
"""OFFLINE DIAGNOSTICS ONLY: the previously reproduced 747 BindingDB additions.

Place in experiments/02_classification/code/ and run:
    python audit_bindingdb_747.py --self-test
    python audit_bindingdb_747.py
Optional: --comparison PATH_TO_FOLDER_OR_ZIP --baseline PATH_TO_CSV

DEPENDENCY: rdkit (standard library for everything else); Python >= 3.10.

AUTHORITY / SCOPE
- The researcher, not this script, decides scientific inclusion criteria.
- No P/N/U/HOLD labels, potency threshold, new dataset, automatic exclusions,
  network access, Git commands, input edits, or application of the withdrawn
  build_trpa1_positive.py rules.
- Exactly one diagnostic row is retained for EVERY input candidate, including
  malformed structures, annotations, conflicting IDs and incomplete evidence.
- The historical parse -> SaltRemover -> largest fragment -> InChIKey procedure
  is repeated ONLY for identity comparison, on in-memory copies. No Cleanup,
  tautomer canonicalization, stereo removal or enumeration is introduced.
- All original SMILES/CXSMILES, inequalities and source metadata are preserved.
  Parsed keys do NOT fully represent CXSMILES/enhanced stereochemical meaning.
- Source-review assertions are copied from the prior comparison, NOT reverified
  here. Keyword matches are navigational hints, not biological classifications.
- No primary patent chemical drawing/table is independently checked by this
  offline program. A database row linking a structure to an IC50 is not proof
  that the source publication has been independently verified.

OUTPUTS (one new dated folder + ZIP, never overwrite an existing path)
  molecule_audit.csv: all 747 candidate rows + diagnostic observations.
  record_audit.csv: ALL cached BindingDB rows with the same saved candidate key,
                    including rows not selected by the previous comparison.
  summary.json: counts, input hashes, versions and declared limitations.
  summary.txt: short Ukrainian explanation and observation-code dictionary.

API documentation used only to implement structural diagnostics:
  https://www.rdkit.org/docs/RDKit_Book.html#cxsmiles-cxsmarts-extensions
  https://www.rdkit.org/docs/source/rdkit.Chem.SaltRemover.html
The optional original BindingDB cache is read locally to compare raw fields.
Its absence is reported, not replaced by a claim of verification.
No documentation lookup or remote request occurs at runtime.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import html
import io
import json
import math
from pathlib import Path
import re
import sys
import tempfile
import zipfile

try:
    from rdkit import Chem, rdBase
    from rdkit.Chem import SaltRemover
except ImportError as exc:
    raise SystemExit('Install in the active environment: python -m pip install rdkit') from exc

VERSION = '1.0.0'
BASELINE_REL = Path('experiments/01_regression/data/trpa1_primary_dataset.csv')
RESULTS_REL = Path('experiments/02_classification/results')
EXPECTED_CANDIDATES = 747
EXPECTED_KEYS = 'efe791b46151b77392d4be7e6683632658d6ba1d437253d55f76d073631eaf7c'
EXPECTED_BASE_KEYS = 'd5afeb066d87f12ec18509c68fe9cd75a7568910ed84005358de483d55e61d9e'
MEMBERS = ('new_bindingdb_molecules.csv', 'new_bindingdb_records.csv',
           'comparison_audit.csv', 'source_summary.csv', 'summary.json')
REMOVER = SaltRemover.SaltRemover()

# These are observations, NEVER reasons to delete/accept a molecule.
NOTES = {
    'raw_bindingdb_record_missing': 'Запис не знайдений у незалежному сирому кеші BindingDB; рядок аудиту збережено.',
    'raw_bindingdb_fields_differ': 'Поля попередньої таблиці відрізняються від відповідного сирого запису BindingDB; відмінності наведено.',
    'smiles_parse_failed': 'RDKit не прочитав структуру; початковий запис збережено.',
    'historical_key_recheck_failed': 'Не вдалося повторно обчислити історичний ключ.',
    'saved_key_not_reproduced': 'Повторно обчислений ключ не збігається зі збереженим.',
    'supplied_bindingdb_key_differs': 'Ключ BindingDB відрізняється від обчисленого; це може бути різниця представлення, не доведена помилка.',
    'cx_extension_present': 'Є CXSMILES-доповнення; сам факт його наявності не означає проблеми.',
    'cx_r_token_present': 'У вихідному доповненні є r; точний текст збережено, абсолютну конфігурацію не домислено.',
    'enhanced_stereo_groups_present': 'RDKit прочитав додаткові стереогрупи; їх збережено в діагностиці.',
    'cx_changes_plain_smiles_key': 'Повний CXSMILES і запис без доповнення дають різні історичні ключі; доповнення не можна мовчки відкинути.',
    'multiple_fragments_after_salt_removal': 'Після історичного видалення солей залишаються кілька фрагментів; набір не змінено.',
    'mixture_word_in_name': 'У назві/синонімах є позначка суміші чи рацемату; це текстова підказка, не встановлений склад зразка.',
    'human_target_fields_not_all_match': 'Не всі поля цього запису відповідають Human/Homo sapiens, O75762, 1 chain. Інші записи молекули не скасовуються.',
    'inhibition_terms_not_found_in_description': 'В описі не знайдено одночасно TRPA1 та слів інгібування; це не доказ відсутності інгібування.',
    'assay_description_missing': 'У кеші немає опису тесту; попередні посилання на патент збережено.',
    'primary_protocol_not_previously_reviewed': 'У вхідному файлі немає попереднього підтвердження протоколу на рівні першоджерела.',
    'variant_or_nonhuman_terms_in_description': 'Опис згадує інший вид/варіант; слово може стосуватися контрольного досліду, потрібен контекст.',
    'desensitization_term_in_description': 'Опис згадує десенситизацію; це лише підказка для читання.',
    'ic50_missing_or_unparseable': 'Числовий IC50 не розібрано; вихідний текст залишено.',
    'ic50_lower_bound_record': 'Є нижня межа IC50 (> або >=); її не названо неактивністю.',
    'ic50_cached_parse_disagrees': 'Нове технічне читання IC50 не збігається з розбором у попередній таблиці.',
    'ec50_also_present': 'Є також поле EC50; без контексту це не доказ агонізму і не причина виключення.',
    'raw_and_supporting_patent_differ': 'Первинний номер патенту і раніше вказане першоджерело різні; старе пояснення збережено, нової корекції немає.',
    'identifier_maps_to_multiple_saved_keys': 'Один ID пов’язаний із кількома ключами в кеші; усі пов’язані ключі та записи наведені.',
    'candidate_key_in_baseline': 'Збережений ключ кандидата вже є в основному наборі.',
    'recomputed_key_in_baseline': 'Повторно обчислений історичний ключ є в основному наборі.',
    'connectivity_block_in_baseline': 'Збігається перший блок InChIKey; це не доведена тотожність стереоізомерів.',
    'chembl_identifier_in_baseline': 'ChEMBL ID цього запису є в основному наборі.',
    'candidate_smiles_key_not_reproduced': 'SMILES у таблиці кандидатів не відтворює її ключ.',
    'candidate_original_cx_missing_from_saved_plain_smiles': 'У сирому записі є CXSMILES, у попередньому стандартному записі його немає. Зміст не вважається автоматично збереженим.',
}



# Original BindingDB field names. These are not eligibility rules.
RAW_FIELDS = {
    'BindingDB Reactant_set_id':'rsid', 'BindingDB MonomerID':'monomer_id',
    'BindingDB Ligand Name':'ligand_name', 'Ligand SMILES':'original_smiles',
    'Ligand InChI Key':'bindingdb_inchikey', 'Target Name':'target_name',
    'Target Source Organism According to Curator or DataSource':'organism',
    'UniProt (SwissProt) Primary ID of Target Chain 1':'accession',
    'Number of Protein Chains in Target (>1 implies a multichain complex)':'target_chains',
    'IC50 (nM)':'ic50_raw_nM', 'Ki (nM)':'ki_raw_nM', 'Kd (nM)':'kd_raw_nM',
    'EC50 (nM)':'ec50_raw_nM', 'Patent Number':'patent_raw', 'Article DOI':'article_doi',
    'PMID':'pmid', 'ChEMBL ID of Ligand':'chembl_id', 'Curation/DataSource':'curation_source',
    'BindingDB Entry DOI':'bindingdb_entry_doi', 'PubChem AID':'pubchem_aid',
}
RAW_MEMBER = 'raw/bindingdb_target_records.json'


def raw_fields(record: dict) -> dict:
    row = {dest: '' if record.get(src) is None else str(record[src]).strip()
           for src, dest in RAW_FIELDS.items()}
    row['assay_ids'] = ' | '.join(record.get('_entry_assay_ids') or [])
    row['assay_description'] = html.unescape(' | '.join(record.get('_assay_descriptions') or []))
    return row


def load_raw_cache(explicit: Path | None, results: Path, meta: dict):
    """Optional independent raw-row comparison. Never download missing caches."""
    expected = meta.get('bindingdb_raw_sha256')
    candidates = [explicit.resolve()] if explicit else []
    if not explicit:
        saved = str(meta.get('bindingdb_cache_path',''))
        name = saved.replace('\\','/').rsplit('/',1)[-1]
        if saved and Path(saved).exists():
            candidates.append(Path(saved))
        if name:
            candidates.extend([results/name, results/Path(name).stem])
        if results.is_dir():
            candidates.extend(sorted(results.glob('TRPA1_external_check_*'), reverse=True))
    seen = set()
    for path in candidates:
        if path in seen or not path.exists():
            continue
        seen.add(path)
        try:
            payload = read_member(path, RAW_MEMBER)
            if not expected or sha(payload) != expected:
                if explicit:
                    raise ValueError('Explicit raw BindingDB cache does not match the comparison hash.')
                continue
            rows = json.loads(payload)
            if not isinstance(rows,list) or not all(isinstance(r,dict) for r in rows):
                raise ValueError('Malformed raw BindingDB cache.')
            index = {}
            for r in rows:
                ident = raw_fields(r)['rsid']
                if not ident or ident in index:
                    raise ValueError('Missing/duplicate rsid in raw BindingDB cache.')
                index[ident] = raw_fields(r)
            return index, dict(status='MATCHED_FROZEN_CACHE', path=str(path), sha256=sha(payload), rows=len(rows)), (path,payload)
        except (OSError,KeyError,zipfile.BadZipFile,ValueError):
            if explicit:
                raise
    if explicit:
        raise ValueError('Specified raw BindingDB cache was not usable.')
    return None, dict(status='NOT_AVAILABLE_NOT_CHECKED',
        explanation='Comparison CSV diagnostics still run; original BindingDB raw rows were not independently rechecked.'), None


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def key_hash(keys) -> str:
    return sha(('\n'.join(sorted(set(keys))) + '\n').encode('utf-8'))


def js(obj) -> str:
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def ids(text: str) -> set[str]:
    return set(re.findall(r'CHEMBL\d+', text or ''))


def parse_csv(data: bytes, required: set[str], name: str) -> list[dict]:
    reader = csv.DictReader(io.StringIO(data.decode('utf-8-sig')))
    if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError(f'{name}: missing or repeated CSV columns.')
    if required - set(reader.fieldnames):
        raise ValueError(f'{name}: missing columns {sorted(required-set(reader.fieldnames))}')
    rows = list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError(f'{name}: malformed CSV row. No report created.')
    return rows


def read_member(path: Path, name: str) -> bytes:
    if path.is_dir():
        return (path / name).read_bytes()
    with zipfile.ZipFile(path) as z:
        found = [n for n in z.namelist() if n == name or n.endswith('/' + name)]
        if len(found) != 1:
            raise ValueError(f'{path.name}: expected one {name}, found {len(found)}.')
        return z.read(found[0])


def bundle(path: Path) -> dict[str, bytes]:
    return {n: read_member(path, n) for n in MEMBERS}


def bundle_signature(data: dict[str, bytes]) -> str:
    # Ignore run paths/timestamps in summary.json, but never ignore evidence edits.
    return sha(js({n: sha(data[n]) for n in MEMBERS if n != 'summary.json'}).encode())


def discover(results: Path) -> Path:
    choices, invalid = [], []
    if not results.is_dir():
        raise ValueError(f'No results directory: {results}. Use --comparison PATH.')
    for p in sorted(results.iterdir()):
        if not p.name.startswith(('bindingdb_surechembl_comparison_', 'TRPA1_external_expansion_')):
            continue
        if not (p.is_dir() or p.suffix.lower() == '.zip'):
            continue
        try:
            data = bundle(p)
            meta = json.loads(data['summary.json'])
            if meta.get('comparison_complete') is not True:
                raise ValueError('comparison_complete is not true')
            choices.append((p, bundle_signature(data)))
        except (OSError, ValueError, KeyError, zipfile.BadZipFile) as exc:
            invalid.append(f'{p.name}: {exc}')
    if invalid:
        print('Some incomplete/unreadable comparisons were not selected:\n  ' + '\n  '.join(invalid))
    if not choices:
        raise ValueError('No complete comparison found. Use --comparison PATH_TO_FOLDER_OR_ZIP.')
    if len({s for _, s in choices}) != 1:
        raise ValueError('Different comparison contents found; refusing to choose silently.\n'
                         + '\n'.join(str(p) for p, _ in choices) + '\nUse --comparison PATH.')
    # A directory and ZIP with identical evidence are interchangeable.
    return max((p for p, _ in choices), key=lambda p: (p.stem if p.suffix == '.zip' else p.name, p.is_dir()))


def find_root() -> Path | None:
    for start in (Path.cwd(), Path(__file__).resolve().parent):
        for p in (start, *start.parents):
            if (p / BASELINE_REL).is_file():
                return p
    return None


def parse_ic50(text: str) -> tuple[str, float | None]:
    m = re.fullmatch(r'\s*(<=|>=|<|>|=|≤|≥)?\s*((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*', text or '')
    if not m:
        return 'unparseable_or_missing', None
    value = float(m[2])
    if not math.isfinite(value) or value <= 0:
        return 'invalid', None
    return {'≤': '<=', '≥': '>='}.get(m[1], m[1]) or '=', value


def number(text: str) -> float | None:
    try:
        n = float(text)
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def historical_key(mol) -> tuple[str, int]:
    """An identity CHECK on a copy; NOT a new stored/model input structure."""
    stripped = REMOVER.StripMol(Chem.Mol(mol), dontRemoveEverything=True)
    frags = Chem.GetMolFrags(stripped, asMols=True, sanitizeFrags=True)
    if not frags:
        return '', 0
    chosen = max(frags, key=lambda f: f.GetNumHeavyAtoms())
    return Chem.MolToInchiKey(chosen), len(frags)


@lru_cache(maxsize=15000)
def inspect_structure(smiles: str) -> dict:
    """Return observations only. Preserve CX text, including unsupported extensions."""
    out = dict(parse_ok=False, parse_error='', historical_key='',
               fragments_before='', fragments_after_salt_removal='', cx_extension='',
               cx_r_token_present=False, stereo_groups_json='[]',
               assigned_tetra_centers='', unspecified_tetra_centers='',
               plain_core_historical_key='', cx_changes_plain_key=False)
    if '|' in smiles:
        out['cx_extension'] = smiles[smiles.find('|'):smiles.rfind('|')+1]
        out['cx_r_token_present'] = bool(re.search(r'(?:\||,)r(?=[:;,|]|$)', out['cx_extension']))
    try:
        with rdBase.BlockLogs():
            # Same parser as the reproduced comparison; never manually strip CX.
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                out['parse_error'] = 'RDKit returned no molecule'
                return out
            out['parse_ok'] = True
            out['fragments_before'] = len(Chem.GetMolFrags(mol))
            out['historical_key'], out['fragments_after_salt_removal'] = historical_key(mol)
            groups = [{'type': str(g.GetGroupType()), 'atom_indices': [a.GetIdx() for a in g.GetAtoms()]}
                      for g in mol.GetStereoGroups()]
            out['stereo_groups_json'] = js(groups)
            centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
            out['assigned_tetra_centers'] = sum(label != '?' for _, label in centers)
            out['unspecified_tetra_centers'] = sum(label == '?' for _, label in centers)
            if out['cx_extension']:
                plain = Chem.MolFromSmiles(smiles.split('|', 1)[0].strip())
                if plain is not None:
                    out['plain_core_historical_key'] = historical_key(plain)[0]
                    out['cx_changes_plain_key'] = bool(out['historical_key'] and
                        out['plain_core_historical_key'] and
                        out['historical_key'] != out['plain_core_historical_key'])
    except (ValueError, RuntimeError, IndexError) as exc:
        out['parse_error'] = f'{type(exc).__name__}: {exc}'
    return out


# Text matches are clues only: "agonist" may describe the stimulating reagent.
MIXTURE = re.compile(r'\bracem(?:ate|ic)\b|\brac[- ]|\bmixture\b|\(\+\s*/\s*-\)|\(±\)', re.I)
INHIBITION = re.compile(r'\binhibit\w*|\bantagonis\w*|\bblock(?:ade|er|ing)?\b', re.I)
VARIANT = re.compile(r'\bmutant\w*|\bchimer\w*|\bmouse\b|\bmurine\b|\brat\b|\bunknown origin\b', re.I)
PATENT = re.compile(r'\b(?:US|WO|EP)\d{6,}(?:[AB]\d)?\b')
EXAMPLE = re.compile(r'\b((?:US|WO|EP)\d{6,}(?:[AB]\d)?)[, ]+\s*(?:(?:Example|Compound|Ex\.)\s*)?([0-9]+[A-Za-z]?(?:[.-][0-9A-Za-z]+)?)', re.I)


def patent_base(s: str) -> str:
    return re.sub(r'[AB]\d$', '', (s or '').replace('-', '').upper())


def compact(values) -> str:
    return ' | '.join(sorted({str(v) for v in values if str(v)}))


def excerpt(text: str) -> str:
    m = re.search(r'TRPA1|inhibit|antagonis', text, re.I)
    if not m:
        return text[:320]
    return text[max(0, m.start()-80):m.start()+350]


def validate_inputs(data: dict[str, bytes], baseline_bytes: bytes, fixed_scope: bool = True):
    meta = json.loads(data['summary.json'])
    if meta.get('comparison_complete') is not True:
        raise ValueError('Input comparison was not completed.')
    m = parse_csv(data['new_bindingdb_molecules.csv'], {'inchikey','standard_smiles','bindingdb_reactant_set_ids'}, 'candidates')
    a = parse_csv(data['comparison_audit.csv'], {'rsid','monomer_id','standard_inchikey','original_smiles',
        'organism','accession','target_chains','ic50_raw_nM','assay_description','source_review'}, 'all records')
    selected = parse_csv(data['new_bindingdb_records.csv'], set(a[0]) if a else {'rsid'}, 'selected records')
    src = parse_csv(data['source_summary.csv'], {'patent_raw','article_doi','source_review'}, 'sources')
    base = parse_csv(baseline_bytes, {'inchikey','std_smiles','molecule_chembl_id'}, 'baseline')
    for label, rows, field in [('candidates',m,'inchikey'), ('baseline',base,'inchikey'),
                                ('all records',a,'rsid'), ('selected records',selected,'rsid')]:
        if not rows or any(not r[field] for r in rows) or len({r[field] for r in rows}) != len(rows):
            raise ValueError(f'{label}: missing or repeated {field}; refusing to silently deduplicate.')
    keys = {r['inchikey'] for r in m}
    basekeys = {r['inchikey'] for r in base}
    if fixed_scope and (len(keys) != EXPECTED_CANDIDATES or key_hash(keys) != EXPECTED_KEYS):
        raise ValueError('This is not the reproduced set of 747 additions. No selection rules were applied.')
    if fixed_scope and (len(base) != 1645 or key_hash(basekeys) != EXPECTED_BASE_KEYS):
        raise ValueError('Baseline differs from the historical 1645 identities. Nothing was changed.')
    if meta.get('new_structure_key_set_sha256') != key_hash(keys):
        raise ValueError('Candidate key-set hash disagrees with comparison metadata.')
    if meta.get('baseline_key_set_sha256') != key_hash(basekeys):
        raise ValueError('Baseline key-set hash disagrees with comparison metadata.')
    counts = meta.get('counts', {})
    for field, actual in [('new_unique_structures',len(m)), ('baseline_molecules',len(base)),
                           ('bindingdb_TRPA1_related_records',len(a)),
                           ('new_records_after_conservative_comparison',len(selected))]:
        if field in counts and int(counts[field]) != actual:
            raise ValueError(f'Input count mismatch: {field}.')
    by_id = {r['rsid']:r for r in a}
    for r in selected:
        if r['rsid'] not in by_id or r != by_id[r['rsid']]:
            raise ValueError(f'Selected row {r["rsid"]} disagrees with comparison_audit.csv.')
        if r['standard_inchikey'] not in keys:
            raise ValueError('Selected record points outside the candidate set.')
    expected_links = defaultdict(set)
    for r in selected:
        expected_links[r['standard_inchikey']].add(r['rsid'])
    for r in m:
        actual = {v.strip() for v in r['bindingdb_reactant_set_ids'].split('|') if v.strip()}
        if actual != expected_links[r['inchikey']]:
            raise ValueError(f'Candidate-to-record links disagree for {r["inchikey"]}.')
    return m, a, selected, src, base, meta


def analyse(molecules: list[dict], all_records: list[dict], selected: list[dict], base: list[dict], raw_index: dict | None = None):
    original_snapshot = js([molecules, all_records, selected, base])
    keys = {r['inchikey'] for r in molecules}
    basekeys = {r['inchikey'] for r in base}
    baseconn = {k[:14] for k in basekeys}
    baseids = set().union(*(ids(r['molecule_chembl_id']) for r in base))
    selected_ids = {r['rsid'] for r in selected}
    baseline_recheck = []
    # Technical identity check ONLY; no biological re-audit of the baseline.
    for r in base:
        s = inspect_structure(r['std_smiles'])
        if s['historical_key'] != r['inchikey']:
            baseline_recheck.append({'inchikey':r['inchikey'], 'recomputed_key':s['historical_key'],
                                     'error':s['parse_error']})
    monomers, chembl = defaultdict(lambda: defaultdict(set)), defaultdict(lambda: defaultdict(set))
    for r in all_records:
        if r['monomer_id'] and r['standard_inchikey']:
            monomers[r['monomer_id']][r['standard_inchikey']].add(r['rsid'])
        for cid in ids(r.get('chembl_id','')):
            if r['standard_inchikey']:
                chembl[cid][r['standard_inchikey']].add(r['rsid'])
    record_out, by_key = [], defaultdict(list)
    for r in sorted(all_records, key=lambda x:(x['standard_inchikey'], x['rsid'])):
        key = r['standard_inchikey']
        if key not in keys:
            continue  # scope, not a scientific filter: only the fixed 747 are audited
        s = inspect_structure(r['original_smiles'])
        notes = set()
        if not s['parse_ok']:
            notes.add('smiles_parse_failed')
        if not s['historical_key']:
            notes.add('historical_key_recheck_failed')
        elif s['historical_key'] != key:
            notes.add('saved_key_not_reproduced')
        supplied = r.get('bindingdb_inchikey','')
        supplied_relation = ('not_supplied' if not supplied else 'not_comparable' if not s['historical_key']
            else 'same' if supplied == s['historical_key'] else 'same_connectivity_block' if
            supplied[:14] == s['historical_key'][:14] else 'different_connectivity_block')
        if supplied_relation in {'same_connectivity_block','different_connectivity_block'}:
            notes.add('supplied_bindingdb_key_differs')
        if s['cx_extension']:
            notes.add('cx_extension_present')
            if '|' not in r.get('standard_smiles',''):
                notes.add('candidate_original_cx_missing_from_saved_plain_smiles')
        if s['cx_r_token_present']:
            notes.add('cx_r_token_present')
        if s['stereo_groups_json'] != '[]':
            notes.add('enhanced_stereo_groups_present')
        if s['cx_changes_plain_key']:
            notes.add('cx_changes_plain_smiles_key')
        if s['fragments_after_salt_removal'] and s['fragments_after_salt_removal'] > 1:
            notes.add('multiple_fragments_after_salt_removal')
        mixture_matches = sorted({m.group(0) for m in MIXTURE.finditer(r.get('ligand_name',''))})
        if mixture_matches:
            notes.add('mixture_word_in_name')
        human_fields = (r['organism'].strip().casefold() in {'human','homo sapiens'}
                        and r['accession'].strip() == 'O75762' and r['target_chains'].strip() == '1')
        if not human_fields:
            notes.add('human_target_fields_not_all_match')
        desc = r['assay_description']
        inhibition_terms = sorted({m.group(0) for m in INHIBITION.finditer(desc)})
        if not desc.strip():
            notes.add('assay_description_missing')
        elif not ('trpa1' in desc.lower() and inhibition_terms):
            notes.add('inhibition_terms_not_found_in_description')
        variant_terms = sorted({m.group(0) for m in VARIANT.finditer(desc)})
        if variant_terms:
            notes.add('variant_or_nonhuman_terms_in_description')
        if re.search(r'\bdesensiti\w*',desc,re.I):
            notes.add('desensitization_term_in_description')
        prior_review = r['source_review'].startswith('human_inhibition_protocol_reviewed')
        if not prior_review:
            notes.add('primary_protocol_not_previously_reviewed')
        relation, value = parse_ic50(r['ic50_raw_nM'])
        if value is None:
            notes.add('ic50_missing_or_unparseable')
        if relation in {'>','>='}:
            notes.add('ic50_lower_bound_record')
        saved_value = number(r.get('ic50_numeric_nM',''))
        if (relation != r.get('ic50_relation','') or
            (value is None) != (saved_value is None) or
            (value is not None and saved_value is not None and
             not math.isclose(value,saved_value,rel_tol=1e-12,abs_tol=0))):
            notes.add('ic50_cached_parse_disagrees')
        if r.get('ec50_raw_nM','').strip():
            notes.add('ec50_also_present')
        if r.get('patent_raw') and r.get('supporting_patent') and patent_base(r['patent_raw']) != patent_base(r['supporting_patent']):
            notes.add('raw_and_supporting_patent_differ')
        conflicts = {}
        for idname, mapping in [('BindingDB:'+r['monomer_id'],monomers.get(r['monomer_id'],{})),
                *[('ChEMBL:'+cid,chembl[cid]) for cid in sorted(ids(r.get('chembl_id','')))]]:
            if len(mapping)>1:
                conflicts[idname] = {k:sorted(v) for k,v in sorted(mapping.items())}
        if conflicts:
            notes.add('identifier_maps_to_multiple_saved_keys')
        overlap_ids = sorted(ids(r.get('chembl_id','')) & baseids)
        if key in basekeys:
            notes.add('candidate_key_in_baseline')
        if s['historical_key'] in basekeys:
            notes.add('recomputed_key_in_baseline')
        if key[:14] in baseconn or (s['historical_key'] and s['historical_key'][:14] in baseconn):
            notes.add('connectivity_block_in_baseline')
        if overlap_ids:
            notes.add('chembl_identifier_in_baseline')
        raw_status, raw_differences = 'NOT_CHECKED_CACHE_UNAVAILABLE', {}
        if raw_index is not None:
            raw = raw_index.get(r['rsid'])
            if raw is None:
                raw_status = 'RECORD_NOT_FOUND'
                notes.add('raw_bindingdb_record_missing')
            else:
                raw_differences = {k: {'comparison_value':r.get(k,''), 'raw_value':v}
                                   for k,v in raw.items() if r.get(k,'') != v}
                raw_status = 'DIFFERENCES_FOUND' if raw_differences else 'FIELDS_MATCH'
                if raw_differences:
                    notes.add('raw_bindingdb_fields_differ')
        d = dict(r)  # preserve every input column; add diagnostics under audit_* only
        d.update(audit_raw_cache_check=raw_status,audit_raw_differences_json=js(raw_differences),
            audit_candidate_key=key, audit_was_in_reproduced_757_records=r['rsid'] in selected_ids,
            audit_database_human_target_fields_match=human_fields,
            audit_description_has_human_term=bool(re.search(r'\bhuman\b|Homo sapiens',desc,re.I)),
            audit_description_inhibition_terms=compact(inhibition_terms),
            audit_description_variant_terms=compact(variant_terms), audit_description_excerpt=excerpt(desc),
            audit_prior_protocol_review_present=prior_review,
            audit_protocol_review_repeated_in_this_run=False,
            audit_primary_structure_measurement_link='NOT_INDEPENDENTLY_CHECKED',
            audit_ligand_name_example_mentions=js([{'patent':x[0], 'compound_label':x[1]}
                for x in EXAMPLE.findall(r.get('ligand_name',''))]),
            audit_ligand_name_patent_mentions=compact(PATENT.findall(r.get('ligand_name',''))),
            audit_ic50_relation=relation,audit_ic50_value_nM='' if value is None else value,
            audit_supplied_key_relation=supplied_relation, audit_mixture_name_matches=compact(mixture_matches),
            audit_identifier_conflicts_json=js(conflicts),audit_baseline_chembl_id_overlap=compact(overlap_ids),
            audit_observation_codes=';'.join(sorted(notes)))
        d.update({'audit_structure_'+k:v for k,v in s.items()})
        record_out.append(d)
        by_key[key].append(d)
    molecule_out = []
    for m in sorted(molecules, key=lambda r:r['inchikey']):
        rr = by_key[m['inchikey']]
        if not rr:
            raise ValueError(f'Candidate {m["inchikey"]} has no cached records.')
        original = [r for r in rr if r['audit_was_in_reproduced_757_records']]
        observations = set().union(*(set(r['audit_observation_codes'].split(';'))-{''} for r in rr))
        if inspect_structure(m['standard_smiles'])['historical_key'] != m['inchikey']:
            observations.add('candidate_smiles_key_not_reproduced')
        d = dict(m)
        d.update(audit_reproduced_records=len(original), audit_all_linked_records=len(rr),
            audit_additional_cached_record_ids=compact(r['rsid'] for r in rr if not r['audit_was_in_reproduced_757_records']),
            audit_all_linked_record_ids=compact(r['rsid'] for r in rr),
            audit_all_structures_parseable=all(r['audit_structure_parse_ok'] for r in rr),
            audit_all_historical_keys_match=all(r['audit_structure_historical_key']==m['inchikey'] for r in rr),
            audit_has_exact_ic50_in_reproduced_records=any(r['audit_ic50_relation']=='=' and r['audit_ic50_value_nM']!='' for r in original),
            audit_has_upper_bound_in_reproduced_records=any(r['audit_ic50_relation'] in {'<','<='} and r['audit_ic50_value_nM']!='' for r in original),
            audit_has_prior_source_protocol_review=any(r['audit_prior_protocol_review_present'] for r in original),
            audit_primary_structure_measurement_link='NOT_INDEPENDENTLY_CHECKED',
            audit_observation_codes=';'.join(sorted(observations)),
            audit_observations_uk=' | '.join(NOTES[k] for k in sorted(observations)),
            researcher_decision='')
        molecule_out.append(d)
    if js([molecules,all_records,selected,base]) != original_snapshot:
        raise AssertionError('Internal error: an input object was modified.')
    if len(molecule_out)!=len(molecules) or {r['inchikey'] for r in molecule_out}!=keys:
        raise AssertionError('Internal error: candidate membership changed.')
    counts = {
        'candidate_molecules_input':len(molecules), 'candidate_molecules_in_report':len(molecule_out),
        'molecules_removed':0, 'molecules_added':0, 'scientific_labels_assigned':0,
        'reproduced_addition_records':len(selected), 'all_cached_records_for_these_candidates':len(record_out),
        'baseline_molecules_read_only':len(base), 'baseline_identity_recheck_differences':len(baseline_recheck),
        'candidate_keys_overlapping_baseline':len(keys & basekeys),
        'with_exact_ic50_in_reproduced_records':sum(r['audit_has_exact_ic50_in_reproduced_records'] for r in molecule_out),
        'upper_bound_only_in_reproduced_records':sum(r['audit_has_upper_bound_in_reproduced_records'] and not r['audit_has_exact_ic50_in_reproduced_records'] for r in molecule_out),
        'with_prior_source_level_protocol_review':sum(r['audit_has_prior_source_protocol_review'] for r in molecule_out),
        'without_prior_source_level_protocol_review':sum(not r['audit_has_prior_source_protocol_review'] for r in molecule_out),
        'primary_compound_structure_links_independently_verified_in_this_run':0,
        'raw_bindingdb_rows_independently_compared':sum(r['audit_raw_cache_check'] in {'FIELDS_MATCH','DIFFERENCES_FOUND'} for r in record_out),
        'raw_bindingdb_rows_with_differences':sum(r['audit_raw_cache_check']=='DIFFERENCES_FOUND' for r in record_out),
    }
    observations = Counter(c for r in molecule_out for c in r['audit_observation_codes'].split(';') if c)
    return molecule_out,record_out,counts,dict(sorted(observations.items())),baseline_recheck


def csv_bytes(rows: list[dict]) -> bytes:
    if not rows:
        raise ValueError('Unexpected empty audit table.')
    fieldnames = list(rows[0])
    if any(set(r)!=set(fieldnames) for r in rows):
        raise ValueError('Inconsistent output columns.')
    stream=io.StringIO(newline='')
    writer=csv.DictWriter(stream,fieldnames=fieldnames,lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return stream.getvalue().encode('utf-8-sig')


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--root',type=Path,help='Local repository root (normally found automatically).')
    p.add_argument('--comparison',type=Path,help='Reproduced comparison folder or ZIP, not the positive_set_v1 output.')
    p.add_argument('--bindingdb-cache',type=Path,help='Optional original external-check folder/ZIP; normally detected by its hash.')
    p.add_argument('--baseline',type=Path,help='Historical 1645-molecule CSV; read-only.')
    p.add_argument('--output-dir',type=Path,help='Exact NEW audit folder; must not already exist.')
    p.add_argument('--self-test',action='store_true',help='Run embedded offline tests only; do not read project data.')
    args=p.parse_args()
    if args.self_test:
        self_test();return 0
    root=args.root.resolve() if args.root else find_root()
    if args.baseline:
        baseline=args.baseline.resolve()
    elif root:
        baseline=root/BASELINE_REL
    else:
        p.error('Repository not found. Use --root PATH or --baseline PATH and --comparison PATH.')
    comparison=args.comparison.resolve() if args.comparison else discover(root/RESULTS_REL if root else baseline.parent)
    output_parent=(root/RESULTS_REL) if root else comparison.parent
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out=args.output_dir.resolve() if args.output_dir else output_parent/('bindingdb_747_audit_'+stamp)
    zipout=out.with_name(out.name+'.zip')
    if out.exists() or zipout.exists():
        raise ValueError('Output path already exists. Refusing to overwrite.')
    print('DIAGNOSTICS ONLY: no exclusions, no P/N/U labels, no network.')
    print('Comparison:',comparison);print('Baseline:  ',baseline);print('RDKit:     ',rdBase.rdkitVersion)
    data=bundle(comparison);basebytes=baseline.read_bytes()
    molecules,records,selected,sources,base,meta=validate_inputs(data,basebytes)
    print(f'Auditing {len(molecules)} candidates; baseline is read-only.',flush=True)
    raw_index,raw_info,raw_input=load_raw_cache(args.bindingdb_cache,output_parent,meta)
    print('Original BindingDB cache:',raw_info['status'])
    mo,ro,counts,observations,base_recheck=analyse(molecules,records,selected,base,raw_index)
    # Integrity guard: no caller input may have changed during the run.
    if baseline.read_bytes()!=basebytes or bundle(comparison)!=data:
        raise ValueError('An input changed during the audit. No completed output written.')
    if raw_input and read_member(raw_input[0],RAW_MEMBER)!=raw_input[1]:
        raise ValueError('Raw BindingDB cache changed during audit. No completed output written.')
    csvs={'molecule_audit.csv':csv_bytes(mo),'record_audit.csv':csv_bytes(ro)}
    summary={
        'script_version':VERSION,'created_utc':datetime.now(timezone.utc).isoformat(),
        'diagnostics_complete':True,'scientific_acceptance_performed':False,
        'primary_sources_reverified':False,'network_used':False,
        'python':sys.version,'rdkit':rdBase.rdkitVersion,'script_sha256':sha(Path(__file__).read_bytes()),
        'comparison_path':str(comparison),'baseline_path':str(baseline),
        'input_member_sha256':{k:sha(v) for k,v in data.items()},'baseline_file_sha256':sha(basebytes),
        'baseline_key_set_sha256':key_hash(r['inchikey'] for r in base),
        'candidate_keys_input_sha256':key_hash(r['inchikey'] for r in molecules),
        'candidate_keys_output_sha256':key_hash(r['inchikey'] for r in mo),
        'candidate_membership_unchanged':True,'input_bytes_unchanged':True,
        'counts':counts,'observations_per_molecule_nonexclusive':observations,
        'original_bindingdb_cache_check':raw_info,
        'baseline_identity_recheck_details':base_recheck,
        'source_review_table_from_previous_comparison':sources,
        'input_bindingdb_release':meta.get('bindingdb_release'),
        'scope': 'Fixed previously reproduced 747 structures; diagnostic joins use saved full InChIKeys. Other-key identifier conflicts are reported, not merged.',
        'limits':[
            'A successful SMILES parse/identity match does not validate experimental stereochemistry.',
            'Prior source-level review is inherited, not repeated; primary compound drawings and value tables are not verified here.',
            'Keyword matches, CXSMILES, racemate names and missing fields do not delete or relabel a candidate.',
            'No claim of comprehensive conflict coverage beyond the supplied BindingDB comparison and baseline.',
            'The withdrawn positive_set_v1 audit and its P/HOLD decisions are not inputs.',
        ],
        'output_sha256':{k:sha(v) for k,v in csvs.items()},
    }
    lines=['АУДИТ 747 ДОДАТКОВИХ СТРУКТУР BINDINGDB — БЕЗ ВІДСІВУ', '',
        f"Вхідних кандидатів: {counts['candidate_molecules_input']}",
        f"Рядків у molecule_audit.csv: {counts['candidate_molecules_in_report']}",
        'Виключено: 0. Додано: 0. Наукових міток призначено: 0.',
        f"Вимірювань у відтвореному додаванні: {len(selected)}; усіх пов’язаних кешованих рядків: {len(ro)}.",
        f"Попередня перевірка протоколу першоджерела є для {counts['with_prior_source_level_protocol_review']} кандидатів; відсутня для {counts['without_prior_source_level_protocol_review']}.",
        'Це успадкований висновок попереднього аудиту, а не повторне читання статті/патенту.',
        'Перевірка полів за оригінальним сирим кешем BindingDB: '+raw_info['status']+'.',
        'Незалежне зіставлення кожного хімічного рисунка з номером сполуки та IC50 НЕ виконувалося.',
        'NOT_INDEPENDENTLY_CHECKED означає межу цього аудиту, а не 747 помилкових молекул.',
        'Порожня колонка researcher_decision залишена для рішення дослідника.',
        '', 'СПОСТЕРЕЖЕННЯ (категорії перетинаються; не є правилами виключення):']
    lines += [f'{n}: {NOTES[k]} [{k}]' for k,n in observations.items()]
    lines += ['', 'ФАЙЛИ:', 'molecule_audit.csv — один рядок на кожного кандидата; первинні поля збережено.',
        'record_audit.csv — усі пов’язані записи: вихідний SMILES/CXSMILES, IC50, опис, ID, посилання й діагностика.',
        'summary.json — контрольні суми, версії, підсумки та таблиця попереднього аудиту джерел.',
        '', 'Вхідні файли залишилися побайтово незмінними. Набір ChEMBL не переоцінено біологічно.',
        'Цей запуск не формує остаточного позитивного датасету. Жодне спостереження не скасовує інші вимірювання молекули.']
    out.mkdir(parents=True,exist_ok=False)
    try:
        for name,content in csvs.items():
            (out/name).write_bytes(content)
        (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        (out/'summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8-sig')
        with zipfile.ZipFile(zipout,'x',zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()):
                z.write(f,f.name)
    except BaseException:
        (out/'INCOMPLETE.txt').write_text('Output writing was interrupted; do not use partial files.\n',encoding='utf-8')
        raise
    print('\nRESULTS')
    for name,value in counts.items():
        print(f'  {name}: {value}')
    print('  candidate_membership_unchanged: True')
    print('  input_bytes_unchanged: True')
    print('Results:',out);print('ZIP:    ',zipout)
    return 0


def self_test():
    """Small deterministic tests; no project writes or scientific curation."""
    assert parse_ic50('<1000')==('<',1000.0)
    assert parse_ic50('>3e3')==('>',3000.0)
    assert parse_ic50('≤ 100')==('<=',100.0)
    assert parse_ic50('0')[1] is None and parse_ic50('Not determined')[1] is None
    assert inspect_structure('not_a_smiles')['parse_ok'] is False
    assert inspect_structure('C[C@H](O)C(=O)O |r|')['cx_r_token_present'] is True
    assert inspect_structure('C[C@H](O)C(=O)O |r|')['parse_ok'] is True
    assert inspect_structure('CCO |(0,0,;1,0,;2,0,)|')['cx_r_token_present'] is False
    assert inspect_structure('CCO.CCC')['fragments_after_salt_removal']==2
    basekey=inspect_structure('c1ccccc1')['historical_key']
    key=inspect_structure('C[C@H](O)C(=O)O |r|')['historical_key']
    base=[dict(inchikey=basekey,std_smiles='c1ccccc1',molecule_chembl_id='CHEMBL1')]
    r=dict(rsid='1',monomer_id='2',ligand_name='Example 7, racemate',
        original_smiles='C[C@H](O)C(=O)O |r|',bindingdb_inchikey=key,
        target_name='TRPA1',organism='Homo sapiens',accession='O75762',target_chains='1',
        ic50_raw_nM='<1000',ec50_raw_nM='',ic50_relation='<',ic50_numeric_nM='1000',
        standard_inchikey=key,standard_smiles='C[C@H](O)C(=O)O',chembl_id='',
        assay_description='Human TRPA1 inhibition; activated by agonist AITC',
        source_review='bindingdb_record_only',patent_raw='',supporting_patent='')
    molecule=dict(inchikey=key,standard_smiles='C[C@H](O)C(=O)O',bindingdb_reactant_set_ids='1')
    mo,ro,c,obs,_=analyse([molecule],[r],[r],base)
    assert len(mo)==1 and c['molecules_removed']==0 and c['scientific_labels_assigned']==0
    assert mo[0]['researcher_decision']=='' and 'cx_extension_present' in obs and 'mixture_word_in_name' in obs
    assert ro[0]['original_smiles']==r['original_smiles'] and ro[0]['ic50_raw_nM']=='<1000'
    assert 'label' not in mo[0] and 'status' not in mo[0]
    mo2,ro2,_,_,_=analyse([molecule],[r],[r],base,{'1':{'rsid':'1','ic50_raw_nM':'123'}})
    assert ro2[0]['audit_raw_cache_check']=='DIFFERENCES_FOUND' and len(mo2)==1
    assert 'raw_bindingdb_fields_differ' in mo2[0]['audit_observation_codes']
    invalid=dict(r,rsid='2',original_smiles='not_a_smiles',ic50_raw_nM='>10',ic50_relation='>',ic50_numeric_nM='10')
    mo,ro,c,obs,_=analyse([molecule],[r,invalid],[r],base)
    assert len(mo)==1 and len(ro)==2 and c['molecules_removed']==0
    assert 'smiles_parse_failed' in obs and 'ic50_lower_bound_record' in obs
    conflict=dict(r,rsid='3',standard_inchikey=basekey,original_smiles='c1ccccc1')
    mo,ro,c,obs,_=analyse([molecule],[r,conflict],[r],base)
    assert 'identifier_maps_to_multiple_saved_keys' in obs and len(mo)==1
    assert 'status' not in mo[0] and mo[0]['researcher_decision']==''
    payload={'molecules.csv':csv_bytes([molecule])}
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);folder=root/'folder';folder.mkdir()
        (folder/'molecules.csv').write_bytes(payload['molecules.csv'])
        archive=root/'test.zip'
        with zipfile.ZipFile(archive,'w') as z:
            z.writestr('enclosing/molecules.csv',payload['molecules.csv'])
        assert read_member(folder,'molecules.csv')==read_member(archive,'molecules.csv')
    print('Self-tests: PASS (bounds, CXSMILES, invalid structures, conflicts, no exclusions, folder/ZIP).')


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print('\nInterrupted. Inputs were not modified.',file=sys.stderr)
        raise SystemExit(130)
    except (ValueError,OSError,KeyError,zipfile.BadZipFile) as exc:
        print(f'ERROR: {exc}\nInputs were not modified; no scientific decisions were applied.',file=sys.stderr)
        raise SystemExit(2)