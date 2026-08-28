# STATUS — TRPA1 ML Project

**Оновлено:** 28 серпня 2026 року  
**Безпосередня мета:** завершити складання повного рукопису, виконати
наскрізну редактуру й передати його науковому керівнику для підготовки
подання до *Fiziologichnyi Zhurnal*.

## Політика GitHub

AI-асистенти мають лише read-only доступ у межах робочого процесу проєкту.
Будь-які зміни готуються як окремі файли; комітить і пушить лише власник
репозиторію.

## Назва статті

> **Машинне навчання для прогнозування інгібування TRPA1 людини: валідація
> за хімічними каркасами, хімічна віддаленість і неоднорідність
> літературних даних**

Українська назва затверджена автором. Англійський варіант перевіряється під
час фінального складання двомовного титульного й резюме-блоку.

## Мета роботи

Оцінити точність прогнозування агрегованого pIC50 інгібування TRPA1 людини
для нових хімічних каркасів і перевірити, чи пов’язана похибка прогнозу з
розбіжністю опублікованих значень для тієї самої сполуки, хімічною
віддаленістю та способом поділу даних.

## Дані

```text
ChEMBL 37
CHEMBL6007 — human TRPA1
exact IC50 records with pChEMBL
2 196 activity records
1 645 standardized compounds
544 Bemis–Murcko scaffolds
97 assays
55 ChEMBL documents
2010–2025
```

Канонічний molecule-level файл:

```text
data/processed/trpa1_primary_dataset.csv
```

## Наукова частина — завершено

- ChEMBL 37 snapshot і provenance;
- стандартизація та molecule-level aggregation;
- structure-only benchmark:
  Morgan, RDKit-15, ChemBERTa, MolFormer × RF/XGBoost;
- три незалежні 5-кратні розділення за хімічними каркасами;
- out-of-fold predictions і fold assignments;
- аудит 97 assays;
- source registry;
- аналізи H1–H5;
- cluster-bootstrap для H2 і H4;
- random-versus-scaffold comparison для H5;
- manuscript tables і figures.

## Основні результати

### H1

RF/XGBoost + Morgan прогнозують агреговане pIC50 для нових каркасів краще
за mean baseline; середнє R² при валідації за каркасами становить приблизно
0,57.

### H2

- 393 multi-assay compounds;
- `ρ = 0,127–0,164`;
- 52 multi-document compounds;
- `ρ = 0,347–0,449`.

Широкий multi-assay зв’язок слабкий. Міждокументний результат сильніший,
але підвибірка мала й невипадкова.

### H3

Допоміжна постановка ligand-versus-random-decoy дала:

```text
CV AUC = 0,9968–0,9983
```

Це methodological contrast, а не точна реплікація Mihai et al. і не
prospective validation.

### H4

Хімічна віддаленість демонструє невеликий, але стабільний зв’язок із
похибкою:

```text
ρ ≈ −0,177
```

### H5

```text
RF:       R² 0,572 scaffold → 0,640 random
XGBoost:  R² 0,567 scaffold → 0,638 random
```

При random split приблизно 79% тестових молекул мають scaffold у train.

## Рукопис — поточний стан

| Компонент | Стан | Файл |
|---|---|---|
| Назва статті | затверджено | `docs/TRPA1_STRUCTURED_ABSTRACTS_AND_TITLES_v2.docx` |
| Українське й англійське резюме | готово | `docs/TRPA1_STRUCTURED_ABSTRACTS_AND_TITLES_v2.docx` |
| Вступ | готово | `docs/TRPA1_INTRODUCTION_v2.docx` |
| Матеріали та методи | готова робоча версія | `docs/TRPA1_METHODS_AUDIT_v2.docx` / `.md` |
| Результати та обговорення | готова актуальна v2; у repo ще v1 | `TRPA1_RESULTS_AND_DISCUSSION_v2.docx` |
| Висновки | готово | `docs/TRPA1_CONCLUSIONS_v1.docx` |
| Список літератури | 20 джерел, готово | `docs/TRPA1_REFERENCE_LIST_v1.docx` |
| Локальні PDF джерел | наявні | `sources/papers/` |
| Таблиці | готово | `results/tables/MANUSCRIPT_*` |
| Рисунки | готово | `results/figures/MANUSCRIPT_*` |
| Єдиний файл рукопису | складається й редагується | локальна робоча версія автора |

### Важлива версійна примітка

У репозиторії зараз лежить:

```text
docs/TRPA1_RESULTS_AND_DISCUSSION_v1.docx
```

Актуальна виправлена версія:

```text
TRPA1_RESULTS_AND_DISCUSSION_v2.docx
```

Перед фінальною фіксацією стану репозиторію `v1` слід замінити `v2` або
чітко перенести `v1` до archive.

## Поточна робота

```text
об’єднати розділи
→ узгодити нумерацію посилань
→ звірити всі числа з канонічними таблицями
→ прибрати внутрішні audit-маркери
→ перевірити українську й англійську термінологію
→ вставити таблиці та рисунки після першої згадки
→ перевірити обсяг і формат журналу
→ передати рукопис науковому керівнику
```

## Невиконаний технічний backlog

Не позначати як завершене:

1. `scripts/build_primary_dataset.py`;
2. robustness run без 27 records із 12 пізніше виключених assays;
3. robustness run з assay-balanced target.

Ці задачі є відтворювальним/QC backlog. Вони не повинні мовчки змінювати
заморожений benchmark. Рішення про їх виконання до формального подання
приймається окремо після складання повного рукопису.

## Що не можна стверджувати

- що assay heterogeneity створює універсальну «стелю» R²;
- що Random Forest остаточно переміг нейромережі;
- що встановлено причинний вплив конкретного assay-протоколу;
- що H3 є точною реплікацією Mihai et al.;
- що модель уже довела здатність prospectively знаходити нові антагоністи.

## Канонічні файли результатів

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

## Наступні дії

1. Завершити єдиний файл рукопису.
2. Замінити репозиторну версію «Результатів та обговорення» на v2.
3. Виконати наскрізну перевірку чисел, посилань, термінів і скорочень.
4. Перевірити відповідність вимогам *Fiziologichnyi Zhurnal*.
5. Передати рукопис науковому керівнику.
6. Окремо вирішити долю трьох технічних QC-задач.

## Правило перенесення контексту

Новий чат спочатку читає:

```text
README.md
STATUS.md
docs/paper_plan.md
docs/methods_fact_sheet.md
results/tables/FINAL_H1_H5_METADATA.json
```

Новий AI-асистент не змінює репозиторій, центральну мету, назву або
зафіксовані результати без прямого погодження автора.
