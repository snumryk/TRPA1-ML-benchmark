# TRPA1 ML Benchmark

> **Канонічна точка входу до проєкту.**  
> Цей README відображає затверджений стан статті за `STATUS.md` і
> `docs/paper_plan.md`.
>
> **Оновлено:** 27 серпня 2026 року.

---

## 1. Проєкт одним реченням

Проєкт оцінює, наскільки точно за молекулярною структурою можна прогнозувати
**агреговане за літературними даними pIC50 інгібування TRPA1 людини для нових
хімічних каркасів**, і перевіряє, як із похибкою прогнозу пов’язані
розбіжність опублікованих значень, хімічна віддаленість сполуки та спосіб
поділу даних.

Основна задача:

```text
стандартизована молекулярна структура
→ модель
→ агреговане літературне pIC50 human TRPA1
```

Основна оцінка виконується для Bemis–Murcko-каркасів, відсутніх у навчальній
частині відповідного fold.

---

## 2. Затверджене наукове питання

> **Наскільки точно можна за молекулярною структурою прогнозувати агреговані
> за літературними даними значення pIC50 інгібування TRPA1 людини для
> хімічних каркасів, не представлених під час навчання, і чи пов’язана
> похибка такого прогнозу з розбіжністю опублікованих значень для тієї самої
> сполуки?**

Безпосередня мета — завершити повний рукопис для подання до
*Fiziologichnyi Zhurnal*.

---

## 3. Походження і версія даних

### Історична послідовність

1. Початковий molecule-level dataset для benchmark був сформований із
   ChEMBL 36 і зберігався під назвою:

   ```text
   trpa1_antagonists.csv
   ```

2. Для ChEMBL 37 було повторено той самий строгий IC50-filter і ту саму
   спрощену RDKit-стандартизацію.

3. Перевірка ChEMBL 36 → ChEMBL 37 показала:

   ```text
   compounds in ChEMBL 36: 1645
   compounds in ChEMBL 37: 1645
   overlap:                 1645
   new:                     0
   removed:                 0
   changed aggregates:      0
   ```

4. Тому заморожений benchmark dataset і його результати залишилися
   незмінними, а для provenance та подальших аналізів у репозиторії
   зафіксовано ChEMBL 37 snapshot.

### Поточні канонічні назви

```text
історична назва:
trpa1_antagonists.csv

поточна назва того самого основного molecule-level набору:
data/processed/trpa1_primary_dataset.csv
```

Історичні run metadata не переписувати: вони мають зберігати фактичну назву
файла, використану під час запуску benchmark.

### Строгий activity filter

```text
target_chembl_id = CHEMBL6007
standard_type = IC50
standard_relation = "="
pchembl_value is not null
```

### Заморожений стан

| Показник | Значення |
|---|---:|
| Придатні activity records | 2 196 |
| Стандартизовані сполуки | 1 645 |
| Bemis–Murcko scaffolds | 544 |
| Assays | 97 |
| ChEMBL documents | 55 |
| Роки документів | 2010–2025 |

Record-level таблиця:

```text
data/raw/trpa1_current_api_raw.csv
```

Основний molecule-level dataset:

```text
data/processed/trpa1_primary_dataset.csv
```

Заморожений ChEMBL 37 snapshot:

```text
data/raw/trpa1_raw_snapshots/ChEMBL_37_20260806T144919Z/
```

### Що означає target

Один рядок у `trpa1_primary_dataset.csv` відповідає одній стандартизованій
молекулярній структурі. Записи об’єднано за InChIKey, а `pchembl_median`
визначено як медіану всіх retained pChEMBL-вимірювань цієї сполуки.

Оскільки всі включені записи мають endpoint `IC50`, у рукописі
`pchembl_median` позначається як агреговане pIC50.

Це **літературно агрегована potency**, а не значення, отримане в одному
стандартизованому assay-протоколі.

### Критичне правило

**Не виконувати нового ChEMBL extraction і не перебудовувати основний набір
із поточного API без окремого рішення автора проєкту.**

---

## 4. Що завершено

- зафіксовано ChEMBL 37 snapshot і provenance;
- підтверджено, що molecule-level dataset ChEMBL 36 і ChEMBL 37 ідентичний
  для 1 645 стандартизованих сполук;
- виконано structure-only benchmark:
  - Morgan ECFP4;
  - RDKit-15 descriptors;
  - frozen ChemBERTa embeddings;
  - frozen MolFormer embeddings;
  - Random Forest;
  - XGBoost;
- виконано три незалежні 5-кратні scaffold partitions;
- використано однакові fold assignments для всіх pipelines;
- збережено out-of-fold predictions, fold assignments і run metadata;
- виконано аудит 97 assays;
- створено source registry для 55 ChEMBL documents;
- виконано затверджені аналізи H1–H5;
- для H2 і H4 застосовано scaffold-cluster bootstrap;
- для H5 виконано cluster bootstrap і permutation comparison;
- підготовлено робочу версію розділу «Матеріали та методи».

---

## 5. Затверджені гіпотези H1–H5

### H1 — прогноз для нових каркасів

Молекулярна структура дозволяє прогнозувати агреговане pIC50 краще за
передбачення середнього значення навіть для каркасів, відсутніх у train.

**Результат:** підтримана. Для RF/XGBoost + Morgan середнє scaffold-CV
`R² ≈ 0,57`.

### H2 — розбіжність опублікованих значень

Серед сполук, виміряних у кількох assays, більша розбіжність assay-level
pIC50 пов’язана з більшою похибкою прогнозу.

**Результат:** виявлено слабку позитивну асоціацію:

```text
393 multi-assay compounds
ρ = 0,127–0,164
```

Для 52 сполук, представлених у кількох документах, асоціація сильніша
(`ρ = 0,347–0,449`), але ця підвибірка мала й невипадкова.

Технічні записи спочатку агрегуються медіаною всередині
`compound × assay` або `compound × document`.

### H3 — ligand-versus-random-decoy classification

Допоміжне відтворення постановки
«TRPA1-ліганд проти випадкової ChEMBL decoy-сполуки» дає майже ідеальні
метрики:

```text
CV AUC = 0,9968–0,9983
```

Це підтримує висновок, що випадкові непідібрані decoys роблять задачу значно
простішою за scaffold-регресію pIC50. Цей аналіз **не є точною реплікацією**
Mihai et al. і не доводить практичної якості prospective virtual screening.

### H4 — хімічна віддаленість

Чим менша Tanimoto similarity тестової молекули до найближчої навчальної
молекули, тим більша похибка.

**Результат:** виявлено невеликий, але стабільний зв’язок:
`ρ ≈ −0,177` для RF і XGBoost.

### H5 — random split проти scaffold split

Random split завищує видиму якість моделі:

```text
RF:       R² 0,572 scaffold → 0,640 random
XGBoost:  R² 0,567 scaffold → 0,638 random
```

При random split приблизно 79% тестових молекул мають scaffold, уже
присутній у train.

---

## 6. Молекулярні представлення і provenance embeddings

### Morgan ECFP4

```text
RDKit Morgan generator
radius = 2
fpSize = 2048
includeChirality = false
```

### RDKit-15

П’ятнадцять фізико-хімічних дескрипторів RDKit.

### ChemBERTa-CLS

```text
model = DeepChem/ChemBERTa-77M-MTR
pooling = first-token / CLS representation
dimensions = 384
max_length = 128
fine-tuning = no
```

### MolFormer-Mean

```text
model = ibm/MoLFormer-XL-both-10pct
revision = 7b12d946c181a37f6012b9dc3b002275de070314
transformers = 4.44.2
trust_remote_code = true
deterministic_eval = true
pooling = attention-mask-aware mean pooling
dimensions = 768
max_length = 202
fine-tuning = no
```

### Що таке `embeddings_all.npz`

Це NumPy-архів із заздалегідь обчисленими матрицями:

```text
rdkit   → 1645 × 15
cb_cls  → 1645 × 384
mf_mean → 1645 × 768
y
scaffold
```

Історичний benchmark input SHA-256:

```text
682d13c98d0565b8e9c79342cb985f27476c43e0846bf2f0af911400685c726e
```

Файл потрібен для **точного повторного запуску повної 4 × 2 сітки без
повторної генерації embeddings**. Він не потрібен для написання рукопису,
аналізів H1–H5 на Morgan або перевірки основного датасету.

Генерація цього файла відтворюється в:

```text
scripts/GroupKFold_CV.ipynb
```

---

## 7. Алгоритми і валідація

### Регресори

```text
RandomForestRegressor
XGBRegressor
```

### Основна scaffold-aware validation

```text
3 scaffold partitions
× 5 folds
groups = Bemis–Murcko scaffold
однакові folds для всіх pipelines
по одному OOF-прогнозу на молекулу в кожному partition
```

Основні метрики:

```text
RMSE
R²
Spearman correlation
```

MAE додатково використано для H5.

Коректний загальний висновок:

> Morgan fingerprints і ChemBERTa embeddings показали близьку якість.
> Стабільної переваги складніших learned representations над Morgan не
> встановлено.

Некоректні твердження:

```text
«Нейронні мережі програли»
«Random Forest є найкращою можливою моделлю»
«Модель уже знаходить нові антагоністи»
```

---

## 8. Канонічний benchmark code

Виконуваний Colab-notebook:

```text
scripts/grid_benchmark.ipynb
```

Він містить:

- завантаження CSV і NPZ;
- сувору перевірку відповідності рядків;
- побудову Morgan fingerprints;
- три scaffold partitions;
- RF/XGBoost;
- mean baseline;
- OOF predictions;
- fold assignments;
- metadata і SHA-256.

Файл:

```text
scripts/Grid_Benchmark.py
```

є історичним notebook-to-text export і не повинен вважатися канонічним
standalone executable.

Це engineering-cleanup, а не підстава сумніватися у вже збережених
benchmark results.

---

## 9. Поточний головний технічний крок

Початковий extraction/standardization pipeline існував до assay-aware
етапу. Його історичні частини збереглися в notebook/scripts, а точний
ChEMBL 36 → ChEMBL 37 результат уже перевірено.

Однак у `main` досі немає одного канонічного repo-relative скрипту, який
із замороженої record-level таблиці відтворює фінальний файл:

```text
data/processed/trpa1_primary_dataset.csv
```

Наступний файл:

```text
scripts/build_primary_dataset.py
```

має:

- читати лише `data/raw/trpa1_current_api_raw.csv`;
- не звертатися до ChEMBL API;
- відтворювати strict IC50 subset;
- виконувати зафіксовану RDKit-стандартизацію;
- агрегувати записи за InChIKey;
- додавати Bemis–Murcko scaffold;
- відтворювати схему і порядок рядків канонічного CSV;
- зберігати metadata, версії пакетів і hashes;
- показувати всі розбіжності, а не приховувати їх.

Це поточний пріоритет перед фіналізацією розділу «Методи».

---

## 10. Дві короткі robustness-перевірки після build-script

Ці перевірки не вимагають повторення повної 4 × 2 сітки.

### A. Пізніше виключені assays

Після assay-аудиту 27 із 2 196 records були віднесені до 12 assays зі
статусом `excluded`.

Мінімальна перевірка:

```text
прибрати 27 records
→ перебудувати molecule-level target
→ повторити лише Morgan + RF/XGBoost
→ порівняти scaffold-CV metrics
```

Мета — перевірити, що головний висновок не залежить від цих 27 records.

### B. Однакова вага для assays

Поточний target:

```text
median усіх retained activity records сполуки
```

Додатковий assay-balanced target:

```text
median усередині compound × assay
→ median між assay-level медіанами
```

Це robustness analysis. Поточний заморожений target не замінювати мовчки і
старі benchmark results не перейменовувати.

---

## 11. Що не є головним напрямом цієї статті

Наступні задачі можуть бути окремою майбутньою роботою, але не повинні
затримувати поточний рукопис:

- повна assay-aware модель із protocol features;
- прогноз pIC50 для конкретного assay;
- повна ручна реконструкція всіх assay-протоколів;
- новий D-MPNN run лише для спроби перевершити RF;
- fine-tuning ChemBERTa або MolFormer;
- prospective screening;
- wet-lab validation;
- остаточне закриття всіх bibliographic `recheck/unresolved` записів, які не
  використовуються як джерела тверджень у рукописі.

Assay audit залишається важливим для опису неоднорідності й обмежень
набору, але assay-aware prediction **не є центральним експериментом цієї
статті**.

---

## 12. Що не можна стверджувати

- що assay heterogeneity створює універсальну «стелю» R²;
- що встановлено причинний вплив конкретного assay-протоколу;
- що neural representations принципово непридатні для TRPA1;
- що random forest остаточно переміг усі інші алгоритми;
- що ligand-versus-random-decoy AUC характеризує реальну якість прогнозу
  potency;
- що модель prospectively підтвердила нові антагоністи TRPA1.

---

## 13. Канонічні файли

### Контекст і план

```text
README.md
STATUS.md
docs/paper_plan.md
docs/methods_fact_sheet.md
docs/TRPA1_METHODS_AUDIT_v2.md
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
scripts/GroupKFold_CV.ipynb
scripts/grid_benchmark.ipynb
results/tables/grid_final_results_20260801-152155.csv
results/tables/grid_final_oof_20260801-152155.csv
results/tables/grid_final_fold_assignments_20260801-152155.csv
results/tables/grid_final_metadata_20260801-152155.json
```

### Аналізи H1–H5

```text
scripts/analyze_h1_h5.py
scripts/test_h5_validation_difference.py
results/tables/FINAL_H1_scaffold_performance.csv
results/tables/FINAL_H2_variability_vs_error.csv
results/tables/FINAL_H3_mihai_replication_summary.csv
results/tables/FINAL_H4_similarity_vs_error.csv
results/tables/FINAL_H5_random_vs_scaffold.csv
results/tables/FINAL_H5_random_oof_morgan.csv
results/tables/FINAL_H5_significance.csv
results/tables/FINAL_H1_H5_METADATA.json
docs/FINAL_H1_H5_CHECKED_REPORT.md
```

---

## 14. Пріоритет джерел істини

Якщо файли суперечать один одному:

```text
заморожені CSV/JSON і snapshot
→ run metadata та hashes
→ виконуваний code/notebooks
→ STATUS.md
→ docs/paper_plan.md
→ README.md
→ робочі звіти
→ історія чатів
```

Історичні metadata не редагувати для «красивішої» назви файла або path:
вони мають відображати фактичний запуск.

---

## 15. Інструкції для нового чату або AI-асистента

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
4. Не перегенеровувати заморожені benchmark outputs новими версіями
   бібліотек без окремого sensitivity run.
5. Не повертати assay-aware prediction у центр поточної статті.
6. Не подавати H3 як точну replication або prospective validation.
7. Усі нові результати супроводжувати input/output hashes, environment
   metadata та відтворюваним кодом.

---

## 16. Поточна виробнича послідовність

```text
1. scripts/build_primary_dataset.py
2. перевірка точного відтворення dataset
3. дві компактні robustness-перевірки
4. остаточна правка «Матеріалів та методів»
5. остаточні рисунки й таблиці
6. «Результати та обговорення»
7. вступ, висновки, українське й англійське резюме
8. повний рукопис для наукового керівника
```
