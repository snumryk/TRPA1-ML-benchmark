"""
Витягує метадані двох документів з яких походять 4 critical assays.
Треба знати які саме патенти шукати.
"""
from chembl_webresource_client.new_client import new_client

DOCS = ["CHEMBL5727828", "CHEMBL5727829"]

doc_client = new_client.document

for did in DOCS:
    d = list(doc_client.filter(document_chembl_id=did))
    if not d:
        print(f"{did}: не знайдено")
        continue
    d = d[0]
    print("=" * 70)
    print(f"{did}")
    print("=" * 70)
    for k in ["title", "doc_type", "patent_id", "year", "journal",
              "authors", "abstract", "doi", "pubmed_id"]:
        v = d.get(k)
        if v:
            print(f"  {k:12}: {str(v)[:200]}")