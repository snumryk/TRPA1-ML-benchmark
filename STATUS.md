# STATUS — TRPA1 ML Project

**Оновлено:** 27 серпня 2026 року  
**Безпосередня мета:** передати науковому керівнику повний рукопис для
подання до *Fiziologichnyi Zhurnal*.  
**Назва статті:** ще не затверджена.

## Затверджене центральне питання

> **Наскільки точно можна за молекулярною структурою прогнозувати агреговані
> за літературними даними значення pIC50 інгібування TRPA1 людини для
> хімічних каркасів, не представлених під час навчання, і чи пов’язана
> похибка такого прогнозу з розбіжністю опублікованих значень для тієї самої
> сполуки?**

## Походження даних

- Початковий molecule-level benchmark dataset був сформований із ChEMBL 36
  під історичною назвою `trpa1_antagonists.csv`.
- Той самий strict IC50-filter і standardization pipeline було повторено для
  ChEMBL 37.
- Перевірка показала 1 645 / 1 645 спільних стандартизованих сполук,
  0 нових, 0 вилучених і 0 змінених aggregate values.
- Поточний канонічний файл:
  `data/processed/trpa1_primary_dataset.csv`.
- Заморожена ChEMBL 37 record-level таблиця містить 2 196 придатних
  IC50-записів із 97 assays і 55 ChEMBL documents.
- Основний benchmark отриманий на агрегованому molecule-level dataset.

## Що завершено

- зафіксовано ChEMBL 37 snapshot;
- підтверджено ідентичність molecule-level наборів ChEMBL 36 і ChEMBL 37;
- виконано structure-only benchmark:
  Morgan, RDKit-15, ChemBERTa, MolFormer × RF/XGBoost;
- виконано три 5-кратні розділення за Bemis–Murcko scaffolds;
- збережено fold assignments і out-of-fold predictions;
- виконано аудит 97 assays і source registry;
- перевірено затверджені гіпотези H1–H5;
- виконано статистичні перевірки H2, H4 і H5;
- підготовлено версію розділу «Матеріали та методи» для source audit;
- у репозиторії є виконувані Colab-notebooks для generation embeddings і
  основного benchmark.

## Затверджені гіпотези та результати

### H1 — прогноз для нових каркасів

**Підтримана.** RF/XGBoost + Morgan прогнозують агреговане pIC50 краще за
mean baseline; середнє scaffold-CV `R² ≈ 0,57`.

### H2 — розбіжність опублікованих значень

**Виявлено слабку позитивну асоціацію на рівні assays.**

- 393 сполуки мають значення щонайменше з двох assays.
- Зв’язок між розкидом pIC50 і похибкою: `ρ = 0,127–0,164`.
- 52 сполуки представлені щонайменше у двох документах.
- Для них зв’язок сильніший: `ρ = 0,347–0,449`, але ця підвибірка мала й
  невипадкова.

Технічні записи всередині одного assay спочатку агрегуються медіаною.

### H3 — допоміжна ligand-versus-random-decoy постановка

**Підтримує methodological contrast.** Випадкова класифікація дає
`AUC = 0,9968–0,9983`, але:

- використовує Morgan замість MNA;
- використовує інший набір і кількість decoys;
- decoys не підібрані за фізико-хімічними властивостями;
- не є точною реплікацією Mihai et al.;
- не доводить practical prospective virtual screening.

Коректна назва: **допоміжне відтворення постановки
ligand-versus-random-decoy**.

### H4 — хімічна віддаленість

**Підтримана як невеликий, але стабільний зв’язок.** Чим менша подібність
до найближчої навчальної молекули, тим більша похибка:
`ρ ≈ −0,177` для RF і XGBoost.

### H5 — random split проти scaffold split

**Підтримана.**

- RF: `R² 0,572 → 0,640`;
- XGBoost: `R² 0,567 → 0,638`;
- приблизно 79% random-test molecules мають scaffold у train.

## Головний висновок

Найсильніше видиму якість моделі змінює спосіб поділу даних. Хімічна
віддаленість демонструє невеликий, але стабільний зв’язок із похибкою.
Розбіжність опублікованих значень пояснює лише невелику частину залишкової
похибки на рівні assays і показує сильніший зв’язок лише в малій
міждокументній підвибірці.

## Що не можна стверджувати

- що assay heterogeneity створює універсальну «стелю» R²;
- що random forest остаточно переміг нейромережі;
- що встановлено причинний вплив конкретного assay-протоколу;
- що H3 є точною реплікацією Mihai et al.;
- що модель уже довела здатність знаходити нові антагоністи prospectively.

## Що більше не є головним напрямом цієї статті

- повна assay-aware модель із неповних protocol fields;
- прогноз pIC50 для конкретного assay;
- пошук нової архітектури лише для перемоги над RF;
- повна реконструкція всіх assay-протоколів перед першим рукописом.

## Підтверджений provenance embeddings

### ChemBERTa

```text
DeepChem/ChemBERTa-77M-MTR
CLS representation
384 dimensions
max_length = 128
```

### MolFormer

```text
ibm/MoLFormer-XL-both-10pct
revision = 7b12d946c181a37f6012b9dc3b002275de070314
transformers = 4.44.2
mean pooling
768 dimensions
max_length = 202
```

Generation notebook:

```text
scripts/GroupKFold_CV.ipynb
```

Canonical benchmark notebook:

```text
scripts/grid_benchmark.ipynb
```

`scripts/Grid_Benchmark.py` є історичним notebook-to-text export, а не
канонічним standalone executable. Це не ставить під сумнів заморожені
результати.

## Поточний головний технічний крок

Підготувати й додати до репозиторію канонічний:

```text
scripts/build_primary_dataset.py
```

Він має відтворити `data/processed/trpa1_primary_dataset.csv` із
`data/raw/trpa1_current_api_raw.csv` без нового API extraction.

Після цього виконати дві компактні robustness-перевірки:

1. Morgan + RF/XGBoost після вилучення 27 records із 12 пізніше виключених
   assays.
2. Morgan + RF/XGBoost із assay-balanced target:
   median усередині `compound × assay`, потім median між assays.

Це не повторення повної 4 × 2 сітки.

## Канонічні файли результатів

```text
results/tables/grid_final_oof_20260801-152155.csv
results/tables/FINAL_H1_scaffold_performance.csv
results/tables/FINAL_H2_variability_vs_error.csv
results/tables/FINAL_H3_mihai_replication_summary.csv
results/tables/FINAL_H4_similarity_vs_error.csv
results/tables/FINAL_H5_random_vs_scaffold.csv
results/tables/FINAL_H5_significance.csv
results/tables/FINAL_H1_H5_METADATA.json
results/figures/FINAL_H5_random_vs_scaffold_both_models.png
docs/FINAL_H1_H5_CHECKED_REPORT.md
scripts/analyze_h1_h5.py
scripts/test_h5_validation_difference.py
```

## Найближчі дії

1. Додати й перевірити `scripts/build_primary_dataset.py`.
2. Виконати дві компактні robustness-перевірки.
3. Внести перевірені corrections у «Матеріали та методи».
4. Побудувати остаточні рисунки H2 і H4 та зведену таблицю.
5. Написати «Результати та обговорення».
6. Завершити повний рукопис.

## Правило перенесення контексту

Новий чат спочатку читає:

```text
README.md
STATUS.md
docs/paper_plan.md
docs/methods_fact_sheet.md
results/tables/FINAL_H1_H5_METADATA.json
```

Після цього він не змінює затверджене питання або план без прямого
погодження автора.
