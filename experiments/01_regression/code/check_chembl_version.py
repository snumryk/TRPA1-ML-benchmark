from chembl_webresource_client.new_client import new_client

# chembl_release endpoint holds version metadata
rel = new_client.chembl_release
records = list(rel.all())
print(f"Total release records: {len(records)}")
print("\nMost recent releases:")
for r in records[-3:]:
    print(f"  {r.get('chembl_release')}: {r.get('creation_date')}")