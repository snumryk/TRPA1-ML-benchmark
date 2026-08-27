# TRPA1 ML Benchmark

> **Канонічна точка входу до проєкту.**  
> Цей README відображає затверджений стан статті за `STATUS.md` і `docs/paper_plan.md`.  
> **Оновлено:** 27 серпня 2026 року.

---

## 1. Проєкт одним реченням

Проєкт оцінює, наскільки точно за молекулярною структурою можна прогнозувати **агреговане за літературними даними pIC50 інгібування TRPA1 людини для нових хімічних каркасів**, і перевіряє, як із похибкою прогнозу пов’язані розбіжність опублікованих значень, хімічна віддаленість сполуки та спосіб поділу даних.

Основна задача:

```text
стандартизована молекулярна структура
→ модель
→ агреговане літературне pIC50 human TRPA1
```

Основна перевірка виконується для Bemis–Murcko-каркасів, відсутніх у навчальній частині відповідного fold.

---

## 2. Затверджене наукове питання

> **Наскільки точно можна за молекулярною структурою прогнозувати агреговані за літературними даними значення pIC50 інгібування TRPA1 людини для хімічних каркасів, не представлених під час навчання, і чи пов’язана похибка такого прогнозу з розбіжністю опублікованих значень для тієї самої сполуки?**

Безпосередня мета проєкту — завершити повний рукопис для подання до *Fiziologichnyi Zhurnal*.

---

## 3. Джерело даних і заморожений набір

Основне джерело — **ChEMBL 37**, target `CHEMBL6007` — human TRPA1.

Строгий activity filter:

```text
target_chembl_id = CHEMBL6007
standard_type = IC50
standard_relation = "="
pchembl_value is not null
```

Зафіксований стан даних:

| Показник | Значення |
|---|---:|
| Придатні activity records | 2 196 |
| Стандартизовані сполуки | 1 645 |
| Bemis–Murcko scaffolds | 544 |
| Assays | 97 |
| ChEMBL documents | 55 |
| Роки документів | 2010–2025 |

Джерело істини для сирих даних:

```text
data/raw/trpa1_raw_snapshots/ChEMBL_37_20260806T144919Z/
```

Заморожена record-level таблиця:

```text
data/raw/trpa1_current_api_raw.csv
```

Основний molecule-level dataset:

```text
data/processed/trpa1_primary_dataset.csv
```

Перехід від 2 196 записів до 1 645 сполук є **агрегацією**, а не відбором 1 645 окремих біологічних експериментів. Записи об’єднано за стандартизованою молекулярною структурою, а цільову змінну `pchembl_median` визначено як медіану збережених значень pChEMBL для сполуки.

У рукописі цей показник позначається як агреговане pIC50, оскільки всі включені activity records мають endpoint `IC50`. Це **літературно агрегована potency**, а не значення, отримане в одному стандартизованому assay-протоколі.

### Критичне правило

**Не виконувати нового ChEMBL extraction і не перебудовувати набір із поточного API без окремого рішення автора проєкту.**

---

## 4. Затверджені гіпотези H1–H5

### H1 — прогноз для нових каркасів

Молекулярна структура дозволяє прогнозувати агреговане pIC50 краще за передбачення середнього значення навіть для каркасів, відсутніх у train.

**Стан:** підтримана. Для RF/XGBoost + Morgan середнє scaffold-CV `R² ≈ 0,57`.

### H2 — розбіжність опублікованих значень

Більша розбіжність assay-level pIC50 пов’язана з більшою похибкою прогнозу.

**Стан:** виявлено слабку позитивну асоціацію на широкій multi-assay підвибірці:

```text
393 multi-assay compounds
ρ = 0,127–0,164
```

Для 52 сполук, представлених у кількох документах, асоціація сильніша (`ρ = 0,347–0,449`), але ця підвибірка мала й невипадкова.

### H3 — ligand-versus-random-decoy classification

Проста класифікація «TRPA1-ліганд проти випадкового decoy» створює значно оптимістичніше враження, ніж scaffold-регресія pIC50.

**Стан:** підтримана як допоміжний контраст. Отримано `AUC = 0,9968–0,9983`, але ця задача не доводить практичної якості prospective virtual screening.

### H4 — хімічна віддаленість

Похибка зростає зі зменшенням подібності тестової сполуки до найближчої навчальної молекули.

**Стан:** підтримана; для RF і XGBoost `ρ ≈ −0,177`.

### H5 — random split проти scaffold split

Random split завищує видиму якість моделі порівняно з поділом за хімічними каркасами.

**Стан:** підтримана:

```text
RF:       R² 0,572 scaffold → 0,640 random
XGBoost:  R² 0,567 scaffold → 0,638 random
```

При random split приблизно 79% тестових молекул мають scaffold, уже присутній у train.

---

## 5. Що завершено

- зафіксовано ChEMBL 37 snapshot і provenance;
- сформовано заморожений molecule-level dataset із 1 645 сполук;
- виконано structure-only benchmark:
  - Morgan ECFP4;
  - RDKit-15 descriptors;
  - frozen ChemBERTa embeddings;
  - frozen MolFormer embeddings;
  - Random Forest;
  - XGBoost;
- виконано три незалежні 5-кратні scaffold partitions;
- використано однакові fold assignments для всіх моделей;
- збережено out-of-fold predictions і run metadata;
- виконано аудит 97 assays і source registry для 55 ChEMBL documents;
- виконано затверджені аналізи H1–H5;
- для H2 і H4 застосовано scaffold-cluster bootstrap;
- для H5 виконано cluster bootstrap і permutation comparison random-versus-scaffold validation;
- підготовлено робочу версію розділу «Матеріали та методи».

---

## 6. Поточний головний технічний розрив

У репозиторії є заморожені сирі дані та готовий файл:

```text
data/processed/trpa1_primary_dataset.csv
```

але немає одного канонічного executable-скрипту, який відтворює весь перехід:

```text
record-level ChEMBL data
→ стандартизація структур
→ molecule-level aggregation
→ Bemis–Murcko scaffolds
→ trpa1_primary_dataset.csv
```

Історична документація описує такий pipeline:

1. RDKit parsing молекулярних структур.
2. Видалення солей або вибір основного фрагмента.
3. Генерація canonical isomeric SMILES.
4. Генерація стандартизованого InChIKey.
5. Deduplication за стандартизованою структурою.
6. Агрегація pChEMBL медіаною.
7. Побудова Bemis–Murcko scaffold.

### Наступний крок

Створити канонічний скрипт побудови основного датасету, робоча назва:

```text
scripts/03_build_primary_dataset.py
```

Скрипт має:

- працювати лише із замороженими файлами репозиторію;
- не звертатися до поточного ChEMBL API;
- мати явні критерії включення та fail-fast перевірки;
- детерміновано стандартизувати структури;
- документувати всі вилучені або непридатні записи;
- формувати `trpa1_primary_dataset.csv` із зафіксованою схемою колонок;
- перевіряти очікувані 2 196 records, 1 645 compounds і 544 scaffolds;
- записувати metadata, версії бібліотек, параметри та SHA-256 входів/виходів;
- або точно відтворити заморожений dataset, або явно показати й пояснити кожну розбіжність.

**Не змінювати мовчки визначення target, правила стандартизації чи склад набору лише для отримання збігу.**

---

## 7. Основний benchmark

### Молекулярні представлення

```text
Morgan ECFP4, radius 2, 2048 bits, без chirality
RDKit-15 physicochemical descriptors
ChemBERTa-CLS, 384 dimensions
MolFormer-Mean, 768 dimensions
```

### Регресійні алгоритми

```text
RandomForestRegressor
XGBRegressor
```

### Основна валідація

```text
3 scaffold partitions
× 5 folds
Bemis–Murcko scaffold grouping
однакові folds для всіх моделей
збережені OOF predictions
```

Основні метрики:

```text
RMSE
R²
Spearman correlation
```

MAE додатково використано для H5.

Коректний узагальнений висновок:

> Morgan fingerprints і ChemBERTa embeddings показали близьку якість. Стабільної переваги складніших learned representations над Morgan не встановлено.

Некоректні висновки:

```text
«Нейронні мережі програли»
«Random Forest є найкращою можливою моделлю»
«Модель уже знаходить нові антагоністи»
```

---

## 8. Що не є головним напрямом цієї статті

Наступні задачі можуть бути окремою майбутньою роботою, але не повинні затримувати поточний рукопис:

- повна assay-aware модель із protocol features;
- прогноз pIC50 для конкретного assay;
- повна ручна реконструкція всіх assay-протоколів;
- новий D-MPNN run лише для спроби перевершити RF;
- fine-tuning ChemBERTa або MolFormer;
- prospective screening;
- wet-lab validation;
- остаточне закриття всіх bibliographic `recheck/unresolved` записів, які не використовуються як джерела тверджень у рукописі.

Assay audit залишається важливим для опису неоднорідності й обмежень набору, але **assay-aware prediction більше не є центральним експериментом цієї статті**.

---

## 9. Що не можна стверджувати

- що assay heterogeneity створює універсальну «стелю» R²;
- що встановлено причинний вплив конкретного assay-протоколу;
- що neural representations принципово непридатні для TRPA1;
- що random forest остаточно переміг усі інші алгоритми;
- що ligand-versus-decoy AUC характеризує реальну якість прогнозу potency;
- що модель prospectively підтвердила нові антагоністи TRPA1.

---

## 10. Канонічні файли

### Контекст і план статті

```text
README.md
STATUS.md
docs/paper_plan.md
docs/methods_fact_sheet.md
docs/TRPA1_METHODS_AUDIT_v2.md
docs/TRPA1_METHODS_AUDIT_v2.docx
```

### Дані

```text
data/raw/trpa1_raw_snapshots/ChEMBL_37_20260806T144919Z/
data/raw/trpa1_current_api_raw.csv
data/processed/trpa1_primary_dataset.csv
data/assays/assay_table.csv
sources/source_registry.csv
```

### Основний benchmark

```text
results/tables/grid_final_results_20260801-152155.csv
results/tables/grid_final_oof_20260801-152155.csv
results/tables/grid_final_fold_assignments_20260801-152155.csv
results/tables/grid_final_metadata_20260801-152155.json
```

### Аналізи H1–H5

```text
results/tables/FINAL_H1_scaffold_performance.csv
results/tables/FINAL_H2_variability_vs_error.csv
results/tables/FINAL_H3_mihai_replication_summary.csv
results/tables/FINAL_H4_similarity_vs_error.csv
results/tables/FINAL_H5_random_vs_scaffold.csv
results/tables/FINAL_H5_random_oof_morgan.csv
results/tables/FINAL_H5_significance.csv
results/tables/FINAL_H1_H5_METADATA.json
docs/FINAL_H1_H5_CHECKED_REPORT.md
scripts/analyze_h1_h5.py
scripts/test_h5_validation_difference.py
```

---

## 11. Пріоритет джерел істини

Якщо файли суперечать один одному, використовувати такий порядок:

```text
заморожені CSV/JSON і snapshot
→ run metadata та hashes
→ executable code
→ STATUS.md
→ docs/paper_plan.md
→ README.md
→ робочі звіти
→ історія чатів
```

README є навігатором і стислим відображенням затвердженого стану, але конкретні числові твердження слід звіряти з канонічними таблицями та metadata.

---

## 12. Інструкції для нового чату або AI-асистента

Перед роботою прочитати:

```text
README.md
STATUS.md
docs/paper_plan.md
docs/methods_fact_sheet.md
results/tables/FINAL_H1_H5_METADATA.json
```

Після цього:

1. Не змінювати затверджене центральне питання без прямого погодження автора.
2. Не робити нового ChEMBL extraction.
3. Не плутати 2 196 record-level measurements із 1 645 molecule-level rows.
4. Не перегенеровувати заморожені benchmark outputs новими версіями бібліотек без окремого sensitivity run.
5. Чітко розрізняти `DONE`, `PRELIMINARY`, `PLANNED` і `BROKEN`.
6. Не повертати assay-aware prediction у центр поточної статті.
7. Не подавати H3 як точну prospective virtual-screening validation.
8. Усі нові результати супроводжувати input/output hashes, environment metadata та відтворюваним кодом.

---

## 13. Поточна виробнича послідовність

```text
1. Канонічний raw → primary dataset build script
2. Перевірка відтворення counts, columns і hashes
3. Остаточні рисунки та зведені таблиці
4. Завершення «Матеріалів та методів»
5. Написання «Результатів та обговорення» за H1–H5
6. Вступ, висновки, українське й англійське резюме
7. Повний рукопис для наукового керівника
```

Поточне безпосереднє завдання:

```text
створити scripts/03_build_primary_dataset.py
```
