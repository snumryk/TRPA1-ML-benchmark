#!/usr/bin/env python3
"""Freeze ONLY the researcher-approved P v1: 1645 ChEMBL + 738 BindingDB.

Run from experiments/02_classification/code/:
    python freeze_trpa1_positive_v1.py
    python freeze_trpa1_positive_v1.py --self-test

Python >=3.10; standard library ONLY. Offline. No RDKit, model training or Git.

This is assembly, NOT curation. It keeps all 1645 historical identities and
all 747 reproduced additions EXCEPT the nine explicitly declined by the
researcher. The nine keys are an immutable decision list below, not a newly
applied source-review/stereo/quality filter. No assay, potency, conflict,
CXSMILES or salt criteria are introduced. Existing SMILES are copied exactly.
Original structure strings, units, inequalities and source columns are retained.
The withdrawn 1664-molecule dataset is NEVER an input.

Inputs already in the repository:
  experiments/01_regression/data/trpa1_primary_dataset.csv
  shared/data/trpa1_current_api_raw.csv
  experiments/02_classification/results/
    bindingdb_surechembl_comparison_20261025T145527_934656Z.zip (or folder)

Outputs: one NEW data/positive_P_v1_2383/ folder and ZIP:
  trpa1_positive_v1.csv  -- exactly 2383 rows, all label=1
  measurements.csv      -- original supporting database rows, not averages
  manifest.json         -- approved decision, input/output hashes and counts

No scientific interpretation is implied by an integrity check. If the frozen
input identities/counts differ, STOP rather than silently repairing/filtering.
The two JSON columns preserve original input fields, including historical
annotations. Historical split/automatic_P_label_assigned fields inside those
JSON values are provenance only; they do not define this experiment's labels
or train/test split. Only the top-level label=1 is the present label.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import zipfile

VERSION = 'P-v1-approved-2383-2026-11-07'
BASELINE = Path('experiments/01_regression/data/trpa1_primary_dataset.csv')
RAW = Path('shared/data/trpa1_current_api_raw.csv')
EXPERIMENT = Path('experiments/02_classification')
COMPARISON_NAME = 'bindingdb_surechembl_comparison_20261025T145527_934656Z'
BASE_KEYS_HASH = 'd5afeb066d87f12ec18509c68fe9cd75a7568910ed84005358de483d55e61d9e'
ADDITION_KEYS_HASH = 'efe791b46151b77392d4be7e6683632658d6ba1d437253d55f76d073631eaf7c'
# The researcher explicitly declined these nine candidates on 2026-11-07.
# DO NOT replace this list with a new general exclusion rule.
DECLINED = {
    'DTGKSKDOIYIVQL-UHFFFAOYSA-N': 'BindingDB 50645704',
    'GZLJEJQXXLODND-UHFFFAOYSA-N': 'BindingDB 191698',
    'KSIOKTDMWIMTCN-UHFFFAOYSA-N': 'BindingDB 191696',
    'LFYXNXGVLGKVCJ-FBIMIBRVSA-N': 'BindingDB 50645717',
    'NZHGWWWHIYHZNX-CSKARUKUSA-N': 'BindingDB 21613 (tranilast)',
    'OMPIXGOSBCVTGS-UHFFFAOYSA-N': 'BindingDB 191699',
    'OYZUXDMMZRXREK-UHFFFAOYSA-N': 'BindingDB 191701',
    'PJTIUFKDAKIUEV-UHFFFAOYSA-N': 'BindingDB 191659',
    'RYYVLZVUVIJVGH-UHFFFAOYSA-N': 'BindingDB 10849 (caffeine)',
}
APPROVED_KEYS_HASH = 'be9974cb3b49694d823b8ea0a261e610974dcdc257701fcf48289befa2eafb88'


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def key_hash(keys) -> str:
    return digest(('\n'.join(sorted(set(keys))) + '\n').encode('utf-8'))


def json_text(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def parse_csv(data: bytes, columns: set[str], name: str) -> list[dict]:
    reader = csv.DictReader(io.StringIO(data.decode('utf-8-sig'), newline=''))
    fields = reader.fieldnames or []
    require(fields and len(fields) == len(set(fields)), f'{name}: invalid header.')
    require(columns <= set(fields), f'{name}: missing columns {columns-set(fields)}.')
    rows = list(reader)
    require(rows and all(None not in r and all(v is not None for v in r.values())
                         for r in rows), f'{name}: empty or malformed CSV.')
    return rows


def csv_bytes(rows: list[dict]) -> bytes:
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode('utf-8-sig')


def read_member(path: Path, name: str) -> bytes:
    if path.is_dir():
        return (path/name).read_bytes()
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n == name or n.endswith('/'+name)]
        require(len(names) == 1, f'{path}: expected one {name}, got {len(names)}.')
        return z.read(names[0])


def unique(rows: list[dict], field: str, name: str) -> set[str]:
    values = [r[field] for r in rows]
    require(all(values) and len(set(values)) == len(values),
            f'{name}: missing/duplicate {field}; nothing was silently deduplicated.')
    return set(values)


def assemble(base: list[dict], chembl: list[dict], additions: list[dict],
             binding: list[dict], declined: set[str]):
    """Copy approved members and their records. Never interpret chemistry."""
    bk = unique(base, 'inchikey', 'Baseline')
    ak = unique(additions, 'inchikey', 'Additions')
    unique(chembl, 'activity_id', 'ChEMBL records')
    unique(binding, 'rsid', 'BindingDB records')
    require(not bk & ak, 'Baseline/additions overlap; no automatic merging allowed.')
    require(declined <= ak and not declined & bk, 'Declined keys not in additions.')
    require({r['inchikey'] for r in chembl} == bk,
            'ChEMBL record identities differ from baseline; no automatic filtering.')
    require({r['standard_inchikey'] for r in binding} == ak,
            'BindingDB record identities differ from the 747-candidate input.')
    by_c, by_b = defaultdict(list), defaultdict(list)
    for r in chembl:
        by_c[r['inchikey']].append(r)
    for r in binding:
        by_b[r['standard_inchikey']].append(r)
    # Verify the frozen links, not biological eligibility.
    for r in base:
        require(len(by_c[r['inchikey']]) == int(r['n_measurements']),
                f'Historical measurement count differs for {r["inchikey"]}.')
    for r in additions:
        linked = {v.strip() for v in r['bindingdb_reactant_set_ids'].split('|') if v.strip()}
        require(linked == {v['rsid'] for v in by_b[r['inchikey']]},
                f'Addition-to-record links differ for {r["inchikey"]}.')

    molecules, measurements = [], []
    for source, rows, groups in [('ChEMBL', base, by_c), ('BindingDB', additions, by_b)]:
        for r in sorted(rows, key=lambda x: x['inchikey']):
            key = r['inchikey']
            if source == 'BindingDB' and key in declined:
                continue  # The ONLY membership exclusion: explicit researcher list.
            records = sorted(groups[key], key=lambda x: x['activity_id'] if source == 'ChEMBL' else x['rsid'])
            originals, evidence_ids, refs = set(), [], set()
            for a in records:
                is_c = source == 'ChEMBL'
                ident = a['activity_id'] if is_c else a['rsid']
                evidence_id = source + ':' + ident
                original = a['canonical_smiles'] if is_c else a['original_smiles']
                originals.add(original)
                evidence_ids.append(evidence_id)
                doc = a['document_chembl_id'] if is_c else a.get('patent_raw', '') or a.get('article_doi', '')
                refs.add(doc)
                measurements.append(dict(
                    inchikey=key, evidence_id=evidence_id, source=source,
                    source_record_id=ident,
                    source_molecule_id=a['molecule_chembl_id'] if is_c else a['monomer_id'],
                    original_smiles=original,
                    stored_standard_smiles=a['std_smiles'] if is_c else a['standard_smiles'],
                    standard_endpoint=a['standard_type'] if is_c else 'IC50',
                    standard_relation=a['standard_relation'] if is_c else a['ic50_relation'],
                    standard_value_nM=a['standard_value'] if is_c else a['ic50_numeric_nM'],
                    standard_units=a['standard_units'] if is_c else 'nM',
                    original_value=a.get('value', a['standard_value']) if is_c else a['ic50_raw_nM'],
                    original_units=a.get('units', a['standard_units']) if is_c else 'nM',
                    original_relation=a.get('relation', a['standard_relation']) if is_c else a['ic50_relation'],
                    assay_id=a['assay_chembl_id'] if is_c else a.get('assay_ids', ''),
                    document_id_original=doc,
                    supporting_patent_from_previous_review='' if is_c else a.get('supporting_patent',''),
                    source_url=('https://www.ebi.ac.uk/chembl/explore/activity/'+ident) if is_c else a.get('source_url',''),
                    original_record_json=json_text(a)))
            molecules.append(dict(
                inchikey=key, smiles=r['std_smiles'] if source=='ChEMBL' else r['standard_smiles'],
                label=1, membership_origin=source,
                original_smiles_json=json_text(sorted(originals)),
                evidence_ids_json=json_text(evidence_ids),
                source_references_json=json_text(sorted(refs-{''})),
                n_database_records=len(records),
                input_molecule_record_json=json_text(r)))
    molecules.sort(key=lambda r: r['inchikey'])
    measurements.sort(key=lambda r: (r['inchikey'], r['source'], r['source_record_id']))
    require({r['inchikey'] for r in molecules} == bk | (ak-declined), 'Assembly membership error.')
    return molecules, measurements


def find_root() -> Path:
    for start in (Path.cwd(), Path(__file__).resolve().parent):
        for p in (start, *start.parents):
            if (p/BASELINE).is_file():
                return p
    raise ValueError('Repository not found. Supply --root C:\\trpa1_data .')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--root', type=Path)
    p.add_argument('--comparison', type=Path, help='The reproduced comparison ZIP or folder.')
    p.add_argument('--output-dir', type=Path, help='Exact NEW output folder; existing paths are refused.')
    p.add_argument('--self-test', action='store_true')
    args = p.parse_args()
    if args.self_test:
        self_test()
        return 0
    root = args.root.resolve() if args.root else find_root()
    comparison = args.comparison.resolve() if args.comparison else root/EXPERIMENT/'results'/(COMPARISON_NAME+'.zip')
    if not args.comparison and not comparison.exists():
        comparison = root/EXPERIMENT/'results'/COMPARISON_NAME
    payloads = {'baseline.csv': (root/BASELINE).read_bytes(),
                'chembl_records.csv': (root/RAW).read_bytes()}
    for name in ['new_bindingdb_molecules.csv','new_bindingdb_records.csv','source_summary.csv','summary.json']:
        payloads[name] = read_member(comparison, name)
    base = parse_csv(payloads['baseline.csv'], {'inchikey','std_smiles','n_measurements'}, 'baseline')
    cr = parse_csv(payloads['chembl_records.csv'],
        {'inchikey','std_smiles','canonical_smiles','activity_id','molecule_chembl_id',
         'document_chembl_id','assay_chembl_id','standard_type','standard_relation','standard_value','standard_units'}, 'ChEMBL records')
    add = parse_csv(payloads['new_bindingdb_molecules.csv'],
        {'inchikey','standard_smiles','bindingdb_reactant_set_ids'}, 'additions')
    br = parse_csv(payloads['new_bindingdb_records.csv'],
        {'standard_inchikey','standard_smiles','original_smiles','rsid','monomer_id',
         'ic50_relation','ic50_numeric_nM','ic50_raw_nM'}, 'BindingDB records')
    info = json.loads(payloads['summary.json'])
    require(info.get('comparison_complete') is True, 'Input comparison was incomplete.')
    require(len(base)==1645 and key_hash(r['inchikey'] for r in base)==BASE_KEYS_HASH, 'Baseline identities differ.')
    require(len(add)==747 and key_hash(r['inchikey'] for r in add)==ADDITION_KEYS_HASH, 'Addition identities differ.')
    require(len(cr)==2196 and len(br)==757, 'Frozen measurement counts differ; nothing was changed.')
    before = json_text([base,cr,add,br])
    molecules, measurements = assemble(base, cr, add, br, set(DECLINED))
    require(json_text([base,cr,add,br]) == before, 'Internal error: input data changed.')
    require(len(molecules)==2383 and key_hash(r['inchikey'] for r in molecules)==APPROVED_KEYS_HASH,
            'Approved P v1 identity check failed.')
    require(len(measurements)==2943, 'Unexpected measurement count.')
    # The historical baseline's exact stored structures are not normalized anew.
    output_map = {r['inchikey']: r for r in molecules}
    require(all(output_map[r['inchikey']]['smiles']==r['std_smiles'] for r in base), 'Baseline structures changed.')
    require(all(output_map[r['inchikey']]['smiles']==r['standard_smiles'] for r in add if r['inchikey'] not in DECLINED),
            'BindingDB stored structures changed.')
    require((root/BASELINE).read_bytes()==payloads['baseline.csv'] and (root/RAW).read_bytes()==payloads['chembl_records.csv'],
            'An input changed during this run.')
    require(all(read_member(comparison,n)==payloads[n] for n in payloads if n not in {'baseline.csv','chembl_records.csv'}),
            'Comparison input changed during this run.')
    out = args.output_dir.resolve() if args.output_dir else root/EXPERIMENT/'data'/'positive_P_v1_2383'
    zipout = out.with_name(out.name+'.zip')
    require(not out.exists() and not zipout.exists(), 'Output already exists; refusing overwrite.')
    tables = {'trpa1_positive_v1.csv':csv_bytes(molecules), 'measurements.csv':csv_bytes(measurements)}
    counts = dict(positive_molecules=2383, from_ChEMBL=1645, from_BindingDB=738,
        approved_BindingDB_exclusions=9, additional_exclusions=0,
        chembl_measurement_rows=2196, bindingdb_measurement_rows=747, measurement_rows=2943)
    manifest = dict(version=VERSION, created_utc=datetime.now(timezone.utc).isoformat(),
        decision='Researcher explicitly approved all 1645 historical ChEMBL plus the 747 reproduced BindingDB additions minus exactly nine listed candidates.',
        declined_BindingDB_candidates=DECLINED,
        exclusion_reason='Researcher declined these nine source-unreviewed candidates; NOT evidence of inactivity.',
        counts=counts, positive_key_set_sha256=key_hash(output_map),
        baseline_membership_and_stored_smiles_unchanged=True, bindingdb_stored_smiles_unchanged=True,
        input_bytes_unchanged=True, new_chemical_processing=False, new_scientific_filters=False,
        full_primary_source_reverification_performed=False, input_paths=dict(baseline=str(root/BASELINE),
        chembl_records=str(root/RAW), bindingdb_comparison=str(comparison)),
        input_sha256={k:digest(v) for k,v in payloads.items()},
        input_bindingdb_release=info.get('bindingdb_release'),
        inherited_source_review_table=parse_csv(payloads['source_summary.csv'], {'source_review'}, 'source summary'),
        script_sha256=digest(Path(__file__).read_bytes()), python=sys.version,
        output_sha256={k:digest(v) for k,v in tables.items()},
        notes=[
            'Database records are not necessarily independent experiments. No cross-database independence is assumed.',
            'Original fields/CXSMILES are preserved; stored model SMILES are copied from the previously reproduced inputs, not reprocessed.',
            'All 1645 baseline members remain, including any earlier audit flags. No claim that those flags have been resolved.',
            'No averaging of BindingDB inequalities, no new potency cutoff, no stereochemical exclusion.',
            'Historical metadata inside input_molecule_record_json are provenance, not a new train/test split or label.',
            'Nine excluded candidates are not negative examples and must not automatically become U.',
            'This replaces the unapproved 1664-molecule proposal; that proposal is not an input.'])
    out.mkdir(parents=True, exist_ok=False)
    try:
        for name, content in tables.items():
            (out/name).write_bytes(content)
        (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        with zipfile.ZipFile(zipout,'x',zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()):
                z.write(f,f.name)
    except BaseException:
        (out/'INCOMPLETE.txt').write_text('Writing interrupted; do not use this partial output.\n',encoding='utf-8')
        raise
    print('APPROVED P v1')
    for key, value in counts.items():
        print(f'  {key}: {value}')
    print('  approved_membership: MATCH\n  original_strings_preserved: True\n  input_bytes_unchanged: True')
    print('Results:',out,'\nZIP:',zipout)
    return 0


def self_test():
    base=[dict(inchikey='A',std_smiles='C[C@H](O)C(=O)O',n_measurements='1')]
    cr=[dict(inchikey='A',std_smiles=base[0]['std_smiles'],canonical_smiles='C[C@H](O)C(=O)O',
        activity_id='1',molecule_chembl_id='CHEMBL1',document_chembl_id='CHEMBL2',
        assay_chembl_id='CHEMBL3',standard_type='IC50',standard_relation='=',
        standard_value='100000',standard_units='nM',value='100',units='uM',relation='=')]
    add=[dict(inchikey='B',standard_smiles='CCO',bindingdb_reactant_set_ids='2'),
         dict(inchikey='C',standard_smiles='CC',bindingdb_reactant_set_ids='3')]
    br=[dict(standard_inchikey=k,standard_smiles=s,original_smiles=s+' |r|',rsid=i,monomer_id=i,
         ic50_relation='<',ic50_numeric_nM='1000',ic50_raw_nM='<1000',
         source_review='flag_does_not_change_membership',patent_raw='US1',supporting_patent='US2')
        for k,s,i in [('B','CCO','2'),('C','CC','3')]]
    snapshot=json_text([base,cr,add,br])
    m,e=assemble(base,cr,add,br,{'C'})
    require([r['inchikey'] for r in m]==['A','B'], 'Wrong test membership.')
    require(all(r['label']==1 for r in m), 'Wrong labels.')
    require(e[1]['original_value']=='<1000' and e[1]['original_smiles']=='CCO |r|', 'Original text lost.')
    require(json.loads(e[1]['original_record_json'])==br[0], 'Original fields lost.')
    require(json_text([base,cr,add,br])==snapshot, 'Inputs mutated.')
    with tempfile.TemporaryDirectory() as tmp:
        folder=Path(tmp)/'folder';folder.mkdir();data=csv_bytes(m)
        (folder/'m.csv').write_bytes(data)
        archive=Path(tmp)/'test.zip'
        with zipfile.ZipFile(archive,'w') as z:z.writestr('nested/m.csv',data)
        require(read_member(folder,'m.csv')==read_member(archive,'m.csv'), 'ZIP/folder mismatch.')
    print('Self-tests: PASS (explicit membership; original fields/CXSMILES/bounds; no input changes; ZIP/folder).')


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print('Interrupted. Inputs unchanged.',file=sys.stderr)
        raise SystemExit(130)
    except (OSError,ValueError,KeyError,zipfile.BadZipFile) as exc:
        print(f'ERROR: {exc}\nNo new inclusion/exclusion rules were applied.',file=sys.stderr)
        raise SystemExit(2)
