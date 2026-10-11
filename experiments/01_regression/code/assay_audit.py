"""
assay_audit.py

Мета простими словами:
Наші 1645 молекул виміряні всього в 97 різних експериментах.
Цей скрипт витягує опис кожного з 97 експериментів і автоматично
позначає ключові слова: яким методом міряли (кальцій / струм),
який активатор використовували (AITC та ін.), в якому порядку
додавали речовину. Результат — таблиця trpa1_assay_audit.csv,
яку далі переглядаєш оком.

Вхідний файл: trpa1_current_api_raw.csv (сирі дані з попереднього кроку,
                                          містять assay_chembl_id для кожного виміру)
"""

import pandas as pd
from chembl_webresource_client.new_client import new_client

RAW_FILE = "trpa1_current_api_raw.csv"
OUT_FILE = "trpa1_assay_audit.csv"

# ============================================================
# 1. Дізнаємось які 97 експериментів у нас є і скільки молекул у кожному
# ============================================================
raw = pd.read_csv(RAW_FILE)
print(f"Загалом вимірів: {len(raw)}")
print(f"Унікальних молекул: {raw['inchikey'].nunique()}")

# Скільки молекул (і вимірів) припадає на кожен експеримент
assay_counts = raw.groupby("assay_chembl_id").agg(
    n_measurements=("inchikey", "count"),
    n_compounds=("inchikey", "nunique"),
).reset_index().sort_values("n_compounds", ascending=False)

assay_ids = assay_counts["assay_chembl_id"].tolist()
print(f"Унікальних експериментів (assays): {len(assay_ids)}\n")

# ============================================================
# 2. Витягуємо опис кожного експерименту з ChEMBL
# ============================================================
print("Витягую описи експериментів з ChEMBL...")
assay_client = new_client.assay

records = []
for i, aid in enumerate(assay_ids, 1):
    try:
        a = assay_client.filter(assay_chembl_id=aid).only([
            "assay_chembl_id",
            "description",
            "assay_type",
            "assay_organism",
            "assay_tissue",
            "assay_cell_type",
            "bao_label",           # стандартна назва типу експерименту
            "confidence_score",    # наскільки надійно прив'язано до мішені
            "document_chembl_id",
        ])
        a = list(a)
        if a:
            records.append(a[0])
        else:
            records.append({"assay_chembl_id": aid, "description": "(not found)"})
    except Exception as e:
        records.append({"assay_chembl_id": aid, "description": f"(error: {e})"})
    if i % 10 == 0:
        print(f"  {i}/{len(assay_ids)}")

df = pd.DataFrame(records)

# Додаємо кількість молекул до кожного опису
df = df.merge(assay_counts, on="assay_chembl_id", how="left")

# ============================================================
# 3. Автоматичне позначення ключових слів у описі
#    (проста перевірка чи є слово в тексті опису)
# ============================================================
def has_any(text, words):
    """Повертає True якщо в тексті є хоч одне зі слів зі списку."""
    if not isinstance(text, str):
        return False
    low = text.lower()
    return any(w in low for w in words)

desc = df["description"].fillna("")

# Метод вимірювання
df["flag_calcium"]     = desc.apply(lambda t: has_any(t, ["calcium", "ca2+", "ca(2+)", "fluo", "fura"]))
df["flag_electrophys"] = desc.apply(lambda t: has_any(t, ["patch", "electrophysiolog", "current", "whole-cell", "voltage clamp"]))
df["flag_fluorescence"]= desc.apply(lambda t: has_any(t, ["fluorescen", "fluo-4", "fura-2"]))

# Який активатор (агоніст) використовували
df["flag_AITC"]        = desc.apply(lambda t: has_any(t, ["aitc", "allyl isothiocyanate", "mustard"]))
df["flag_cinnamald"]   = desc.apply(lambda t: has_any(t, ["cinnamaldehyde", "cinnamald"]))
df["flag_other_agonist"]= desc.apply(lambda t: has_any(t, ["agonist", "activation", "activated", "evoked", "induced"]))

# Порядок додавання / тип дії
df["flag_preincub"]    = desc.apply(lambda t: has_any(t, ["preincub", "pre-incub", "pretreat", "pre-treat"]))
df["flag_antagonist"]  = desc.apply(lambda t: has_any(t, ["antagonist", "inhibit", "block"]))
df["flag_agonist_self"]= desc.apply(lambda t: has_any(t, ["agonist activity", "as an agonist", "channel activation by"]))

# Мутантний чи дикий тип
df["flag_mutant"]      = desc.apply(lambda t: has_any(t, ["mutant", "mutation", "c621", "n855", "chimera"]))

# ============================================================
# 4. Зведення — що маємо загалом
# ============================================================
print("\n" + "="*70)
print("ЗВЕДЕННЯ ПО ЕКСПЕРИМЕНТАХ")
print("="*70)
print(f"Усього експериментів: {len(df)}")
print(f"\nЗа методом (кількість експериментів, що згадують):")
print(f"  кальцій:          {df['flag_calcium'].sum()}")
print(f"  електрофізіологія:{df['flag_electrophys'].sum()}")
print(f"  флуоресценція:    {df['flag_fluorescence'].sum()}")
print(f"\nЗа активатором:")
print(f"  AITC/гірчична олія:{df['flag_AITC'].sum()}")
print(f"  cinnamaldehyde:   {df['flag_cinnamald'].sum()}")
print(f"  інший/загальний:  {df['flag_other_agonist'].sum()}")
print(f"\nЗа типом дії:")
print(f"  preincubation:    {df['flag_preincub'].sum()}")
print(f"  антагоніст/інгіб.:{df['flag_antagonist'].sum()}")
print(f"  мутант/химера:    {df['flag_mutant'].sum()}")

# Скільки МОЛЕКУЛ припадає на найбільші експерименти
print(f"\n{'='*70}")
print("НАЙБІЛЬШІ ЕКСПЕРИМЕНТИ (за кількістю молекул):")
print("="*70)
top = df.sort_values("n_compounds", ascending=False).head(15)
for _, r in top.iterrows():
    d = str(r["description"])[:90] if isinstance(r["description"], str) else "(no description)"
    print(f"  {r['n_compounds']:>4} молекул | {r['assay_chembl_id']} | {d}")

# ============================================================
# 5. Зберігаємо таблицю для ручного перегляду
# ============================================================
# Впорядковуємо колонки зручно: спершу кількість молекул, опис, потім прапорці
col_order = [
    "assay_chembl_id", "n_compounds", "n_measurements", "description",
    "bao_label", "assay_organism", "assay_cell_type", "confidence_score",
    "flag_calcium", "flag_electrophys", "flag_fluorescence",
    "flag_AITC", "flag_cinnamald", "flag_other_agonist",
    "flag_preincub", "flag_antagonist", "flag_agonist_self", "flag_mutant",
    "document_chembl_id",
]
col_order = [c for c in col_order if c in df.columns]
df = df[col_order].sort_values("n_compounds", ascending=False)

df.to_csv(OUT_FILE, index=False)
print(f"\n{'='*70}")
print(f"Збережено: {OUT_FILE}")
print(f"Тепер цей файл можна відкрити і переглянути оком усі 97 експериментів.")
print("="*70)

# ============================================================
# 6. Оцінка: скільки молекул у найбільшій ОДНОРІДНІЙ групі
#    (кальцій + AITC — найпоширеніший чистий сценарій)
# ============================================================
homogeneous = df[df["flag_calcium"] & df["flag_AITC"]]
if len(homogeneous) > 0:
    n_homog_compounds = raw[raw["assay_chembl_id"].isin(homogeneous["assay_chembl_id"])]["inchikey"].nunique()
    print(f"\nОРІЄНТОВНО: експериментів 'кальцій + AITC': {len(homogeneous)}")
    print(f"           молекул у них: ~{n_homog_compounds}")
    print("(Це попередня оцінка розміру чистої однорідної групи —")
    print(" точну визначиш після ручного перегляду таблиці.)")
else:
    print("\nЕкспериментів з явним поєднанням 'кальцій + AITC' в описі не знайдено —")
    print("треба дивитись описи вручну, ключові слова можуть бути іншими.")