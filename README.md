# TRPA1 ML Benchmark

> **Головна точка входу до проєкту.**  
> Перед будь-якою роботою новий дослідник або AI-асистент має спочатку прочитати цей README, а вже потім відкривати окремі CSV, скрипти та старі звіти.
>
> **Остання звірка стану:** 20 серпня 2026 року.

---

## 1. Проєкт одним реченням

Мета проєкту — передбачати **reported inhibitory potency human TRPA1** за структурою молекули та перевірити, чи покращується прогноз, якщо разом зі структурою врахувати **умови біологічного assay**.

Простими словами:

```text
Модель 1: структура молекули → pIC50

Модель 2: структура молекули + умови експерименту → pIC50
```

Головне нове питання:

> Чи пояснюють тип assay, агоніст, клітинна система, час преінкубації та інші параметри протоколу частину розкиду опублікованих IC50/pIC50 для human TRPA1?

---

## 2. Чому assay context важливий

TRPA1 — складний іонний канал. Виміряне IC50 може залежати не лише від структури інгібітора, а й від:

- активатора каналу: AITC, cinnamaldehyde, Ca²⁺, Zn²⁺ та ін.;
- концентрації або рівня активатора, наприклад EC50 чи EC80;
- клітинної системи: CHO, HEK293 та їхніх похідних;
- способу реєстрації: calcium fluorescence, radiometric calcium, electrophysiology;
- порядку додавання речовин;
- часу контакту тестової сполуки з клітинами;
- конструкта каналу;
- fit-процедури та інших деталей протоколу.

Тому цільова змінна цього проєкту — не «універсальна афінність молекули до TRPA1», а:

> **reported inhibitory potency measured in heterogeneous human TRPA1 assays.**

Концептуально:

```text
reported pIC50
≈ внесок структури молекули
+ внесок умов assay
+ експериментальний шум
```

---

## 3. Джерело даних і зафіксований snapshot

Основне джерело — **ChEMBL 37**, release date **2026-05-01**.

Target:

```text
CHEMBL6007 = human TRPA1
```

Строгий activity filter:

```text
standard_type = IC50
standard_relation = "="
pchembl_value is not null
target_chembl_id = CHEMBL6007
```

Зафіксований набір:

| Показник | Значення |
|---|---:|
| Raw activities після API filter | 2 203 |
| Придатні activity records після стандартизації | 2 196 |
| Унікальні стандартизовані молекули | 1 645 |
| Assays | 97 |
| ChEMBL documents | 55 |
| Роки документів | 2010–2025 |
| Assays із `confidence_score = 9` | 87 |
| Assays із `confidence_score = 8` | 10 |

Ролі 97 assays у `data/assays/assay_table.csv`:

| Роль | Assays | Значення |
|---|---:|---|
| `primary` | 76 | консервативне ядро для assay-aware аналізу |
| `sensitivity` | 9 | homologous-target / score-8 assays, потенційно корисні для додаткової перевірки |
| `excluded` | 12 | механістично непридатні, mutant/chimera, direct agonism/desensitization, дублікати або інші причини виключення |

### Критичне правило

**Не робити нового ChEMBL extraction і не перебудовувати dataset з поточного API без окремого рішення автора проєкту.**

Джерело істини:

```text
data/raw/trpa1_raw_snapshots/ChEMBL_37_20260806T144919Z/
```

---

## 4. Що вже зроблено

### Завершено

- зафіксовано ChEMBL 37 snapshot і provenance;
- стандартизовано структури;
- сформовано molecule-level dataset із 1 645 молекул;
- виконано аудит усіх 97 assays;
- створено бібліографічний registry для 55 ChEMBL documents;
- автоматично розібрано assay descriptions на окремі параметри протоколу;
- виконано structure-only benchmark:
  - Morgan ECFP4;
  - RDKit-15 descriptors;
  - ChemBERTa embeddings;
  - MolFormer embeddings;
  - RF і XGBoost;
  - 3 scaffold partitions × 5 folds;
- збережено fold assignments, out-of-fold predictions, metadata та результати;
- знайдено дві великі paired-assay пари для однакових compounds.

### Попередньо виконано, але не готово для сильних висновків

- D-MPNN та D-MPNN + descriptors;
- fine-tuning ChemBERTa/MolFormer;
- feature importance;
- paired-assay causal interpretation;
- protocol families, побудовані лише за автоматично витягнутими полями.

### Ще не виконано

- фінальний QC protocol columns у record-level assay-aware dataset;
- чесне порівняння:
  - structure-only model;
  - structure + assay conditions model;
- matched homogeneous-subset analysis;
- document-held-out stress test;
- остаточне статистичне парне порівняння assay-aware моделей.

---

## 5. Що показав уже виконаний structure-only benchmark

Повна сітка `4 molecular representations × 2 regressors` виконана й лежить у `results/tables/`.

Коректний попередній висновок:

> Morgan fingerprints і ChemBERTa embeddings показали близьку якість. Стабільної переваги складніших learned representations над Morgan не встановлено.

Некоректні твердження:

```text
"Старі алгоритми однозначно виграли."
"GNN остаточно програла."
"Нейромережі не працюють для TRPA1."
```

D-MPNN runs мають технічні обмеження: early stopping, seed variability і нерівний tuning budget.

---

## 6. Поточний головний експеримент

### Наукове питання

Чи покращується прогноз окремого assay-level pIC50, якщо моделі дати інформацію про assay protocol?

### Вхідний файл

```text
data/processed/trpa1_assay_aware_primary_conservative.csv
```

Поточний стан файла:

| Показник | Значення |
|---|---:|
| Activity rows | 1 946 |
| Унікальні молекули | 1 444 |
| Assays | 76 |
| Documents | 44 |
| Колонки | 35 |
| Роки | 2010–2025 |
| Confidence score | усі рядки = 9 |
| Units | усі рядки = nM |

Один рядок означає:

```text
конкретна молекула
+ конкретний assay
+ конкретне вимірювання IC50/pIC50
+ доступні параметри протоколу
```

### Порівняння

#### Model 1 — structure-only control

```text
canonical_smiles → molecular representation → predicted pic50
```

#### Model 2 — assay-aware

```text
canonical_smiles
+ method
+ agonist category
+ agonist level/concentration
+ cell line
+ preincubation
+ application order
+ construct
→ predicted pic50
```

Обидві моделі повинні:

- навчатися на тих самих records;
- використовувати ті самі folds;
- не мати exact molecule/scaffold leakage;
- оцінюватися парно на тих самих out-of-fold predictions.

Primary metrics:

```text
RMSE
MAE
R²
Spearman correlation
```

---

## 7. Два різні ML-файли — не плутати

### `data/processed/trpa1_primary_dataset.csv`

```text
одна молекула = один рядок
```

- 1 645 молекул;
- pChEMBL агрегований медіаною між різними вимірюваннями та assays;
- використаний у вже виконаному structure-only benchmark;
- не містить окремого assay context для кожного вимірювання.

### `data/processed/trpa1_assay_aware_primary_conservative.csv`

```text
одне вимірювання в одному assay = один рядок
```

- одна молекула може повторюватися;
- зберігає `assay_chembl_id`, `document_chembl_id` і protocol fields;
- використовується для поточного assay-aware експерименту.

**Не замінювати один файл іншим.**

---

## 8. Як таблиці пов’язані між собою

### Exact activity record

```text
trpa1_assay_aware_primary_conservative.activity_id
=
trpa1_current_api_raw.activity_id
```

### Assay

```text
trpa1_assay_aware_primary_conservative.assay_chembl_id
=
assay_table.assay_id
=
trpa1_assay_audit.assay_chembl_id
=
ключ у trpa1_assay_cache.json
```

### ChEMBL document / source

Прямий join:

```text
trpa1_assay_aware_primary_conservative.document_chembl_id
=
source_registry.document_id
```

Альтернативний маршрут:

```text
assay_chembl_id
→ assay_table.assay_id
→ assay_table.source_id
→ source_registry.source_id
```

### Молекула

```text
molecule_chembl_id
canonical_smiles
```

Для standardized structure identity в molecule-level dataset використовується `inchikey`.

---

## 9. Структура репозиторію

```text
TRPA1-ML-benchmark/
│
├── README.md
├── TRPA1_ML_benchmark.ipynb
├── env.yml
│
├── data/
│   ├── raw/
│   │   ├── raw_trpa1_human.csv
│   │   ├── trpa1_current_api_raw.csv
│   │   ├── trpa1_current_api_metadata.json
│   │   └── trpa1_raw_snapshots/
│   │       └── ChEMBL_37_20260806T144919Z/
│   │
│   ├── processed/
│   │   ├── trpa1_primary_dataset.csv
│   │   ├── trpa1_assay_aware_primary_conservative.csv
│   │   ├── trpa1_current_api_aggregated.csv
│   │   ├── trpa1_human_clean.csv
│   │   ├── trpa1_agonists.csv
│   │   ├── decoys_raw.csv
│   │   ├── decoys_clean.csv
│   │   └── trpa1_v36_v37_*.csv
│   │
│   └── assays/
│       ├── assay_table.csv
│       ├── trpa1_assay_audit.csv
│       ├── trpa1_assay_audit_metadata.json
│       ├── trpa1_assay_audit_summary.csv
│       ├── trpa1_assay_cache.json
│       ├── protocol_decomposition.csv
│       ├── protocol_families.csv
│       ├── protocol_metadata.json
│       └── paired_contrasts.csv
│
├── sources/
│   ├── source_registry.csv
│   └── papers/
│
├── results/
│   ├── tables/
│   ├── figures/
│   └── checkpoints/
│
├── scripts/
│   ├── 01_fetch_trpa1_raw_snapshot.py
│   ├── 02_summarize_trpa1_snapshot.py
│   ├── assay_audit_reviewed.py
│   ├── protocol_decomposition.py
│   ├── Grid_Benchmark.py
│   ├── chemberta_benchmark.py
│   ├── chemprop_baseline.py
│   ├── morgan_dmpnn_cv.py
│   ├── MihaiExperimentReplication.py
│   ├── feature_importance.py
│   ├── confound_check.py
│   ├── ablation_study.py
│   ├── prepare_decoys.py
│   ├── resolve_delta.py
│   └── інші exploratory/legacy scripts
│
└── docs/
    ├── paper_plan.md
    ├── technical_report_v2.md
    ├── Assey database state.docx
    ├── Assey database state.odt
    ├── Повний аудит 97 assays інгібування human TRPA1.pdf
    └── archive/
```

---

## 10. Роль основних файлів

### `data/raw/`

Незмінені або максимально близькі до сирих дані.

#### `trpa1_raw_snapshots/ChEMBL_37_20260806T144919Z/`

Зафіксований source-of-truth snapshot:

- `activities_raw.jsonl` — raw activity responses;
- `assays_raw.jsonl` — raw assay responses;
- `assay_ids.txt` — список 97 assays;
- `target_CHEMBL6007.json` — target metadata;
- `snapshot_manifest.json` — склад і hashes snapshot;
- `chembl_status_start.json`, `chembl_status_end.json` — версія ChEMBL;
- `SNAPSHOT_COMPLETE.txt` — marker успішного завершення.

#### `trpa1_current_api_raw.csv`

Очищена record-level таблиця з 2 196 activity records. Головний join для `activity_id`.

#### `raw_trpa1_human.csv`

Розширений старіший raw export. Використовувати лише коли потрібного поля немає в `trpa1_current_api_raw.csv`, і перевіряти provenance.

### `data/processed/`

Готові або проміжні datasets.

#### `trpa1_primary_dataset.csv`

Molecule-level aggregated dataset для завершеного structure-only benchmark.

#### `trpa1_assay_aware_primary_conservative.csv`

Record-level dataset для поточного assay-aware експерименту.

#### `trpa1_current_api_aggregated.csv`

Агрегація поточного strict ChEMBL pull за standardized InChIKey. Технічний intermediate.

#### `trpa1_human_clean.csv`

Історичний clean molecule-level dataset. Не вважати автоматично актуальнішим за `trpa1_primary_dataset.csv`.

#### `trpa1_agonists.csv`

Окрема історична/додаткова вибірка agonist records. Не є входом поточного inhibitory-potency experiment.

#### `decoys_raw.csv`, `decoys_clean.csv`

Decoys для replication ligand-vs-random-decoy classification. Не використовуються в regression pIC50 experiment.

#### `trpa1_v36_v37_*.csv`

Результати release-to-release comparison. Це version-stability check, не external validation.

### `data/assays/`

Assay metadata, курація і protocol decomposition.

#### `assay_table.csv`

Людська карта 97 assays:

```text
assay_id
n_compounds
dataset_role
document_id
year
assay_summary
source_id
```

#### `trpa1_assay_audit.csv`

Технічна таблиця — один рядок на assay. Містить ChEMBL metadata, automatic flags та manual-review columns.

#### `trpa1_assay_cache.json`

Повні кешовані ChEMBL assay records. Використовується для перевірки descriptions, cell type, confidence score, target та document ID.

#### `protocol_decomposition.csv`

Один рядок на assay. Автоматично витягнуті protocol fields, evidence fragments і confidence labels.

#### `protocol_families.csv`

Групування assays за однаковими **витягнутими** полями. Це не доводить біологічну однорідність.

#### `paired_contrasts.csv`

Кандидатні assay pairs із shared compounds. Результати потребують source-level verification.

### `sources/`

#### `source_registry.csv`

Бібліографічний/provenance registry для 55 ChEMBL documents:

```text
source_id
document_id
year
assay_ids
source_type
title
identifier
status
local_source
notes
```

Статуси:

- `verified` — source mapping підтверджений;
- `partial` — source відомий, але protocol або exact mapping неповний;
- `recheck` — суперечність або непідтверджена патентна сім’я;
- `unresolved` — точне джерело ще не встановлено.

**`verified source` не означає автоматично `verified assay protocol`.**

#### `papers/`

Локальні статті та Supporting Information. Наявність PDF не означає, що assay mapping уже перевірений.

### `results/`

- `tables/` — metrics, OOF predictions, folds, metadata;
- `checkpoints/` — neural-network checkpoints;
- `figures/` — рисунки.

Canonical record конкретного benchmark run:

```text
grid_final_results_<timestamp>.csv
grid_final_oof_<timestamp>.csv
grid_final_fold_assignments_<timestamp>.csv
grid_final_metadata_<timestamp>.json
```

Файли з `checkpoint` у назві можуть дублювати timestamped outputs; звіряти SHA.

### `docs/`

Робочі тексти та історичні звіти.

- `technical_report_v2.md` — найповніший історичний технічний звіт;
- `paper_plan.md` — план structure-only representation paper;
- DOCX/ODT/PDF — assay audit у людському форматі.

Якщо старий звіт конфліктує з README або фактичними CSV/metadata, пріоритет має:

```text
CSV/JSON data
→ run metadata
→ code
→ README
→ technical reports
→ chat history
```

### `scripts/`

Скрипти різного рівня готовності. Частина є exploratory або Colab-oriented.

Ключові:

- `01_fetch_trpa1_raw_snapshot.py` — створення snapshot; **не запускати для поточного проєкту без окремого рішення**;
- `02_summarize_trpa1_snapshot.py` — опис snapshot;
- `assay_audit_reviewed.py` — assay-level audit;
- `protocol_decomposition.py` — parsing protocol fields;
- `Grid_Benchmark.py` — завершений structure-only 4×2 benchmark;
- `morgan_dmpnn_cv.py`, `chemprop_baseline.py` — graph-model experiments;
- `MihaiExperimentReplication.py` — ligand-vs-decoy replication;
- `feature_importance.py`, `confound_check.py`, `ablation_study.py` — supplementary analyses.

---

## 11. `trpa1_primary_dataset.csv`: опис колонок

| Колонка | Значення |
|---|---|
| `inchikey` | standardized molecular identity |
| `standard_type` | activity endpoint; у strict subset — `IC50` |
| `std_smiles` | standardized canonical isomeric SMILES |
| `molecule_chembl_id` | representative ChEMBL molecule ID |
| `pchembl_median` | median pChEMBL між усіма retained measurements молекули |
| `pchembl_min` | minimum reported pChEMBL |
| `pchembl_max` | maximum reported pChEMBL |
| `pchembl_std` | SD між measurements; порожнє при одному measurement |
| `n_measurements` | кількість measurements для молекули |
| `assay_types` | ChEMBL assay-type codes; **не надійна biological modality** |
| `year_min` | найраніший document year |
| `year_max` | найпізніший document year |
| `scaffold` | Bemis–Murcko scaffold |
| `split` | історичний single train/test split; не заміна current CV folds |

`pchembl_median` — агрегована літературна мітка, а не protocol-specific potency.

---

## 12. `trpa1_assay_aware_primary_conservative.csv`: опис усіх колонок

### Ідентифікатори та provenance

| Колонка | Значення |
|---|---|
| `activity_id` | primary key конкретного ChEMBL activity record; найкращий join із raw CSV |
| `record_id` | ChEMBL compound-record identifier; provenance |
| `molecule_chembl_id` | ChEMBL molecule ID |
| `parent_molecule_chembl_id` | parent molecule ID; у поточному файлі фактично дублює molecule ID |
| `canonical_smiles` | структура для побудови molecular representation |
| `assay_chembl_id` | ID конкретного assay |
| `document_chembl_id` | ID ChEMBL document/source record |
| `document_year` | рік документа |

### Activity value

| Колонка | Значення |
|---|---|
| `standard_value_num` | IC50 у standard units; **не подавати моделі, бо це target leakage** |
| `standard_units` | units; у поточному файлі `nM` |
| `pic50` | обчислений target: `-log10(IC50 in molar)` |
| `reported_pchembl` | ChEMBL-reported pChEMBL; **не predictive feature** |
| `pchembl_delta` | `reported_pchembl - pic50`; QC різниці округлення |

### Protocol fields

| Колонка | Значення |
|---|---|
| `method` | normalized modality: calcium fluorescence, electrophysiology, radiometric calcium, unknown |
| `agonist` | точний normalized challenge activator |
| `agonist_model_category` | укрупнена категорія агоніста для моделювання |
| `agonist_level` | EC50, EC80 тощо |
| `agonist_concentration_uM` | числова концентрація активатора в µM |
| `cell_line` | нормалізована cell-line family, зараз переважно CHO/HEK293 |
| `compound_preincubation_min` | час контакту test compound до agonist addition |
| `activity_time_min` | змішане time field; **не використовувати до semantic QC** |
| `application_order` | порядок додавання, наприклад `compound_before_agonist` |
| `protocol_fingerprint` | composite diagnostic label; **не model feature у поточному стані** |
| `construct` | wild type, mutant/chimera або not reported |

### QC та курація

| Колонка | Значення |
|---|---|
| `confidence_score` | ChEMBL target confidence |
| `data_validity_comment` | ChEMBL validity note |
| `potential_duplicate` | duplicate flag |
| `source_review_status` | рівень protocol/source review |
| `source_note` | людська нотатка про protocol verification |
| `compound_assay_measurement_count` | число measurements для тієї самої molecule × assay |
| `is_compound_assay_replicate` | чи має molecule × assay більше одного record |
| `include_primary_conservative` | inclusion flag для conservative subset |
| `include_broad_sensitivity` | inclusion flag для broad/sensitivity logic |
| `primary_exclusion_reasons` | причини виключення з conservative subset |
| `broad_exclusion_reasons` | причини виключення з broad subset |

---

## 13. Що таке `protocol_fingerprint`

Формат:

```text
record_species
| construct
| assay_modality
| challenge_agonist
| agonist_level
| agonist_concentration
| cell_line
| application_order
| compound_preincubation_min
```

Приклад:

```text
human | not_reported | calcium_fluorescence | NA | NA | NA | HEK293 | NA | NA
```

Це лише компактна мітка для групування assays.

Окремі protocol columns уже існують. Розбирати fingerprint назад не потрібно.

### Поточне обмеження

Fingerprint був сформований раніше, ніж частину окремих protocol fields виправили. У поточному CSV він не узгоджується з окремими колонками для частини assays.

Правило:

```text
Не використовувати protocol_fingerprint як ML feature.
```

Після завершення QC його можна:

- видалити з model-input table;
- або перегенерувати з фінальних окремих колонок.

---

## 14. Поточне покриття protocol fields

Покриття на рівні 76 assays:

| Поле | Assays із корисним значенням |
|---|---:|
| `method` | 65 / 76 |
| `agonist` / `agonist_model_category` | 52 / 76 |
| `cell_line` | 44 / 76 |
| `agonist_concentration_uM` | 12 / 76 |
| `application_order` | 8 / 76 |
| `agonist_level` | 7 / 76 |
| `compound_preincubation_min` | 7 / 76 |
| `activity_time_min` | 6 / 76 |
| `construct` | 1 / 76 |

Не оцінювати coverage лише за кількістю рядків: великі patent assays повторюють однакові protocol values для сотень молекул.

---

## 15. Відомі проблеми assay-aware table

Перед моделюванням потрібен assay-level QC.

### 1. `protocol_fingerprint` застарів

Не використовувати до регенерації.

### 2. `activity_time_min` семантично неоднорідний

Колонка може змішувати:

- preincubation;
- measurement duration;
- time after treatment;
- cell-expression/incubation duration.

Також є floating artifacts на кшталт `10.002` і `4.9998`.

Не використовувати як feature до розділення на конкретні типи часу.

### 3. Частина `compound_preincubation_min` пропущена

У ChEMBL descriptions є явні значення, які parser не завжди переніс у фінальний CSV.

Відомі кандидати для ручної перевірки:

```text
CHEMBL3789590
CHEMBL3791636
CHEMBL4022611
CHEMBL4049011
CHEMBL5360087
```

### 4. `method` потребує точкового QC

Наприклад, `[45]Ca²⁺ influx by microbeta plate count` не слід називати fluorescence assay.

### 5. `cell_line` надто грубо нормалізована

Поточні сімейства CHO/HEK293 можуть приховувати:

```text
CHO-K1
CHO-TREX
HEK293F
HEK293T
HEK293-TREx
T-REx-293
IMR-90
```

Бажано мати дві колонки:

```text
cell_line_family
cell_line_exact
```

### 6. Replicates

Не трактувати кілька measurements однієї molecule × assay як незалежні біологічні observations без наперед визначеного правила агрегації.

### 7. Provenance семи вилучених records

У 76 primary assays у сирому наборі було 1 953 records, а у conservative file — 1 946. Причина вилучення семи records має бути відтворювано зафіксована в коді або audit output.

---

## 16. Які колонки дозволено подавати моделі

### Structure input

```text
canonical_smiles
```

Із нього генеруються:

- Morgan fingerprint;
- RDKit descriptors;
- frozen embedding;
- graph representation.

### Target

```text
pic50
```

### Candidate assay-aware features після QC

```text
method
agonist_model_category
agonist_level
agonist_concentration_uM
cell_line
compound_preincubation_min
application_order
construct
```

### Не використовувати як predictive features

```text
standard_value_num
reported_pchembl
pchembl_delta
activity_id
record_id
assay_chembl_id
document_chembl_id
source_review_status
source_note
protocol_fingerprint
```

Причини:

- target leakage;
- запам’ятовування assay/document identity;
- текстова або ручна інформація, недоступна для справжнього нового prediction;
- дублювання protocol fields.

IDs можна використовувати лише для joins, grouping і leakage-safe splitting.

---

## 17. Правила валідації assay-aware experiment

1. Structure-only і assay-aware models використовують ті самі records.
2. Exact molecule не може потрапити одночасно в train і test.
3. Основне групування — molecular scaffold; exact molecules автоматично залишаються в одному fold.
4. Assay ID, document ID і source ID не є features.
5. Unknown protocol values не вгадуються; вони залишаються explicit `unknown`.
6. Feature preprocessing fit лише на train fold.
7. Усі out-of-fold predictions зберігаються.
8. Різниця між моделями оцінюється парно на тих самих records/scaffolds.
9. Покращення predictive performance не доводить причинний ефект конкретного protocol variable.
10. Document-held-out evaluation — окремий stress test, не заміна scaffold validation.

---

## 18. Найближчий план роботи

### Крок 1 — QC protocol fields

Перевірити 76 assays, а не 1 946 рядків:

```text
method
agonist_model_category
agonist_level
agonist_concentration_uM
cell_line
compound_preincubation_min
application_order
construct
```

### Крок 2 — виправити record-level table

- внести manual corrections;
- розділити `activity_time_min` на семантично окремі поля або не використовувати;
- додати `cell_line_exact`;
- перегенерувати або видалити `protocol_fingerprint`;
- зафіксувати причину вилучення семи records;
- зберегти metadata і SHA-256.

### Крок 3 — заморозити experimental protocol

До перегляду результатів зафіксувати:

- molecular representation;
- estimator;
- folds;
- protocol features;
- missing-value handling;
- primary metrics;
- statistical comparison.

### Крок 4 — Model 1 проти Model 2

```text
Model 1: structure only
Model 2: structure + verified assay features
```

### Крок 5 — secondary analyses

- matched homogeneous subset;
- document-held-out stress test;
- paired-assay analysis після exact patent mapping.

---

## 19. Інструкції для AI-асистента в новому чаті

1. Спочатку прочитай цей README.
2. Не реконструюй стан проєкту з пам’яті або старого чату.
3. Не роби нового ChEMBL extraction.
4. Для чисел і стану роботи використовуй:
   ```text
   data files
   → metadata
   → code
   → README
   → reports
   → chat
   ```
5. Чітко розрізняй:
   ```text
   DONE
   PRELIMINARY
   PLANNED
   BROKEN
   ```
6. Не називай source mapping verified, якщо збігається лише назва.
7. Не заповнюй unknown protocol fields припущеннями.
8. Не використовуй `protocol_fingerprint` як готову ML feature.
9. Не повертайся до structure-only benchmark як до «наступного експерименту» — він уже виконаний.
10. Поточна робота:
   ```text
   QC trpa1_assay_aware_primary_conservative.csv
   → structure-only vs assay-aware prediction
   ```

---

## 20. Відтворюваність і відомі engineering gaps

### `env.yml`

Поточний файл не містить повного ML environment. Перед фінальним запуском потрібно зафіксувати щонайменше:

```text
python
numpy
pandas
scipy
rdkit
scikit-learn
xgboost
transformers
torch
chemprop
```

### `embeddings_all.npz`

Файл використовувався structure-only benchmark, але не зберігається в репозиторії. Run metadata містить його SHA-256.

### Colab-oriented scripts

Деякі скрипти містять:

- Google Drive paths;
- notebook magic;
- exploratory code;
- historical filenames.

Перед фінальним відтворюваним запуском їх слід перетворити на звичайні CLI scripts із repo-relative paths.

---

## 21. Короткий checkpoint для нового чату

```text
PROJECT:
Human TRPA1 inhibitory-potency ML with assay heterogeneity.

DATA SOURCE:
Frozen ChEMBL 37 snapshot, target CHEMBL6007, exact IC50 records.

COMPLETED:
Structure-only molecular-representation benchmark and 97-assay audit.

CURRENT DATASET:
data/processed/trpa1_assay_aware_primary_conservative.csv
1946 records, 1444 molecules, 76 primary assays.

CURRENT TASK:
QC protocol columns, then compare:
structure-only vs structure + assay conditions.

DO NOT:
re-extract ChEMBL,
use protocol_fingerprint as feature,
use assay/document IDs as features,
confuse aggregated molecule-level and record-level datasets.
```
