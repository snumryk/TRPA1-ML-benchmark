# TRPA1 ML Benchmark

> **Канонічна точка входу до проєкту.**  
> Актуальний стан статті та найближчі дії визначаються разом із `STATUS.md`
> і `docs/paper_plan.md`.
>
> **Оновлено:** 28 серпня 2026 року.

## Політика доступу до репозиторію

**AI-асистенти використовують цей репозиторій тільки для читання.**

Заборонено безпосередньо створювати, оновлювати, видаляти або комітити файли,
створювати pull requests чи виконувати будь-які інші write-операції. Усі
підготовлені зміни передаються автору проєкту як окремі файли; публікує їх
лише власник репозиторію.

---

## 1. Проєкт одним реченням

Проєкт оцінює, наскільки точно за молекулярною структурою можна прогнозувати
агреговане за літературними даними pIC50 інгібування TRPA1 людини для
хімічних каркасів, не представлених під час навчання, та як із похибкою
прогнозу пов’язані неоднорідність опублікованих значень, хімічна
віддаленість і спосіб валідації.

## 2. Назва статті

**Затверджена українська назва:**

> **Машинне навчання для прогнозування інгібування TRPA1 людини: валідація
> за хімічними каркасами, хімічна віддаленість і неоднорідність
> літературних даних**

Робочий англійський переклад зберігається у
`docs/TRPA1_STRUCTURED_ABSTRACTS_AND_TITLES_v2.docx` і остаточно
перевіряється під час складання рукопису.

## 3. Мета роботи

Оцінити точність прогнозування агрегованого pIC50 інгібування TRPA1 людини
для нових хімічних каркасів і перевірити, чи пов’язана похибка прогнозу з
розбіжністю опублікованих значень для тієї самої сполуки, її хімічною
віддаленістю від навчальної вибірки та способом поділу даних.

---

## 4. Дані

Основне джерело — **ChEMBL 37**, target `CHEMBL6007` — human TRPA1.

Строгий activity filter:

```text
target_chembl_id = CHEMBL6007
standard_type = IC50
standard_relation = "="
pchembl_value is not null
```

Заморожений стан:

| Показник | Значення |
|---|---:|
| Activity records | 2 196 |
| Стандартизовані сполуки | 1 645 |
| Bemis–Murcko scaffolds | 544 |
| Assays | 97 |
| ChEMBL documents | 55 |
| Роки документів | 2010–2025 |

Основні файли:

```text
data/raw/trpa1_raw_snapshots/ChEMBL_37_20260806T144919Z/
data/raw/trpa1_current_api_raw.csv
data/processed/trpa1_primary_dataset.csv
```

Один рядок `trpa1_primary_dataset.csv` відповідає одній стандартизованій
молекулярній структурі. Цільова змінна `pchembl_median` є медіаною retained
pChEMBL-вимірювань сполуки. Оскільки всі включені записи мають endpoint
`IC50`, у рукописі цей показник позначено як агреговане pIC50.

**Не виконувати нового ChEMBL extraction і не замінювати заморожений набір
даними з поточного API без окремого рішення автора проєкту.**

---

## 5. Що завершено

### Дані та аудит

- зафіксовано ChEMBL 37 snapshot і provenance;
- сформовано molecule-level dataset із 1 645 сполук;
- виконано аудит 97 assays;
- створено source registry для 55 ChEMBL documents;
- підготовлено характеристики набору для рукопису.

### Основний benchmark

Виконано порівняння:

```text
Morgan ECFP4
RDKit-15 descriptors
ChemBERTa-CLS
MolFormer-Mean
×
Random Forest
XGBoost
```

Основна валідація:

```text
3 незалежні scaffold partitions
× 5 folds
groups = Bemis–Murcko scaffold
однакові folds для всіх pipelines
збережені out-of-fold predictions
```

Канонічні notebooks:

```text
scripts/GroupKFold_CV.ipynb
scripts/grid_benchmark.ipynb
```

### Аналізи H1–H5

- H1 — прогноз агрегованого pIC50 для нових каркасів;
- H2 — розбіжність опублікованих значень і похибка прогнозу;
- H3 — допоміжне відтворення постановки ligand-versus-random-decoy;
- H4 — хімічна віддаленість і похибка;
- H5 — random split проти scaffold split.

Збережено канонічні таблиці, out-of-fold predictions, статистичні
перевірки, рисунки та metadata.

### Матеріали рукопису

У репозиторії наявні:

```text
docs/TRPA1_METHODS_AUDIT_v2.md
docs/TRPA1_METHODS_AUDIT_v2.docx
docs/TRPA1_INTRODUCTION_v2.docx
docs/TRPA1_CONCLUSIONS_v1.docx
docs/TRPA1_STRUCTURED_ABSTRACTS_AND_TITLES_v2.docx
docs/TRPA1_REFERENCE_LIST_v1.docx
```

У `sources/papers/` збережено локальні PDF джерел, використаних для
підготовки рукопису.

Розділ «Результати та обговорення»:

```text
репозиторій: docs/TRPA1_RESULTS_AND_DISCUSSION_v1.docx
актуальна робоча версія: TRPA1_RESULTS_AND_DISCUSSION_v2.docx
```

Перед фіналізацією репозиторну `v1` слід замінити актуальною `v2`.

Підготовлено основний набір ілюстрацій і таблиць:

```text
results/figures/MANUSCRIPT_Fig1_workflow.png
results/figures/MANUSCRIPT_Fig2_variability_vs_error.png
results/figures/MANUSCRIPT_Fig3_similarity_vs_error.png
results/figures/MANUSCRIPT_Fig4_random_vs_scaffold_R2.png
results/tables/MANUSCRIPT_Table1_dataset_characteristics.csv
results/tables/MANUSCRIPT_Table2_scaffold_benchmark.csv
```

---

## 6. Основні результати

### H1 — прогноз для нових каркасів

RF/XGBoost + Morgan прогнозують агреговане pIC50 краще за mean baseline.
Середнє R² при валідації за хімічними каркасами становить приблизно 0,57.

### H2 — неоднорідність опублікованих значень

Для 393 сполук, виміряних щонайменше у двох assays, зв’язок між розкидом
pIC50 і похибкою є слабким:

```text
ρ = 0,127–0,164
```

Для 52 сполук із кількох документів зв’язок сильніший:

```text
ρ = 0,347–0,449
```

Ця міждокументна підвибірка мала й невипадкова.

### H3 — ligand-versus-random-decoy

Допоміжна класифікаційна постановка дала:

```text
CV AUC = 0,9968–0,9983
```

Це демонструє легкість задачі з випадковими непідібраними decoys, але не є
точною реплікацією Mihai et al. і не доводить якість prospective virtual
screening.

### H4 — хімічна віддаленість

Менша Tanimoto similarity до найближчої навчальної молекули пов’язана з
більшою похибкою:

```text
ρ ≈ −0,177
```

Зв’язок невеликий, але стабільний для RF і XGBoost.

### H5 — спосіб валідації

Random split дав оптимістичніші оцінки:

```text
RF:       R² 0,572 scaffold → 0,640 random
XGBoost:  R² 0,567 scaffold → 0,638 random
```

При random split приблизно 79% тестових молекул мають scaffold, уже
представлений у train.

---

## 7. Поточний стан рукопису

Наукові аналізи H1–H5 завершені. Окремі розділи, резюме, список літератури,
таблиці та рисунки підготовлено.

**Поточна робота:**

```text
складання всіх частин в один рукопис
→ наскрізна редактура
→ звірка чисел, посилань і термінології
→ форматування за вимогами Fiziologichnyi Zhurnal
→ передача науковому керівнику
```

Не починати новий великий ML-експеримент під час складання рукопису без
окремого рішення автора.

---

## 8. Невиконані технічні задачі

Ці задачі не слід називати завершеними:

1. Канонічний repo-relative build-script:

   ```text
   scripts/build_primary_dataset.py
   ```

   Він має відтворити `trpa1_primary_dataset.csv` із замороженої
   record-level таблиці без нового API extraction.

2. Коротка robustness-перевірка після вилучення 27 records із 12 assays,
   пізніше позначених як `excluded`.

3. Коротка robustness-перевірка з assay-balanced target:

   ```text
   median усередині compound × assay
   → median між assay-level медіанами
   ```

Ці перевірки можуть бути виконані як окремий QC до формального подання, але
не повинні мовчки змінювати заморожений benchmark або вже зафіксовані
результати H1–H5.

---

## 9. Що не можна стверджувати

- що assay heterogeneity визначає універсальну «стелю» R²;
- що встановлено причинний вплив конкретного assay-протоколу;
- що neural representations принципово непридатні для TRPA1;
- що Random Forest є найкращою можливою моделлю;
- що H3 є точною реплікацією Mihai et al.;
- що ligand-versus-random-decoy AUC характеризує реальну якість прогнозу
  potency;
- що модель prospectively підтвердила нові антагоністи TRPA1.

---

## 10. Канонічні файли результатів

```text
results/tables/grid_final_results_20260801-152155.csv
results/tables/grid_final_oof_20260801-152155.csv
results/tables/grid_final_fold_assignments_20260801-152155.csv
results/tables/grid_final_metadata_20260801-152155.json

results/tables/FINAL_H1_scaffold_performance.csv
results/tables/FINAL_H2_variability_vs_error.csv
results/tables/FINAL_H3_mihai_replication_summary.csv
results/tables/FINAL_H4_similarity_vs_error.csv
results/tables/FINAL_H5_random_vs_scaffold.csv
results/tables/FINAL_H5_random_oof_morgan.csv
results/tables/FINAL_H5_significance.csv
results/tables/FINAL_H1_H5_METADATA.json
```

---

## 11. Пріоритет джерел істини

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

1. Використовувати GitHub лише для читання.
2. Не змінювати затверджену назву, мету або план без прямого погодження.
3. Не робити нового ChEMBL extraction.
4. Не плутати 2 196 record-level measurements із 1 645 molecule-level rows.
5. Не повертати assay-aware prediction у центр поточної статті.
6. Не подавати H3 як точну replication або prospective validation.
7. Чітко розрізняти завершені аналізи, робочі рукописні файли та
   невиконаний технічний backlog.
