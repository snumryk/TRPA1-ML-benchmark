# Технічний звіт v2: assay-aware ML-дослідження інгібувальної активності human TRPA1

**Дата зрізу:** 27 липня 2026 року  
**Статус:** робочий технічний звіт для незалежного рецензування; це ще не рукопис статті.  
**Автор проєкту:** аспірант 2-го року, біофізика іонних каналів, Інститут фізіології ім. О.О. Богомольця НАН України.  
**Обчислювальне середовище:** Windows 11 / Google Colab; Python 3.13 локально; RDKit, scikit-learn, XGBoost, Chemprop, Transformers.

## Прохання до незалежного рецензента

Цей документ має бути достатнім для оцінки без доступу до попереднього листування. Просимо критично оцінити:

1. Чи коректно сформульоване центральне наукове питання.
2. Чи не переоцінено роль assay heterogeneity.
3. Чи коректний ML-дизайн і статистика.
4. Чи достатня новизна після врахування робіт 2020, 2023 і 2025 років.
5. Які аналізи є обов’язковими до submission, а які можна винести в майбутню роботу.
6. Чи є фатальні проблеми з курацією ChEMBL, агрегацією IC50, scaffold split або paired-assay analysis.
7. Який найкращий рівень журналу для роботи після запропонованих виправлень.

Просимо не пом’якшувати критику. Особливо цінні конкретні falsification tests: які результати могли б спростувати нашу поточну інтерпретацію.

---

# 0. Резюме проєкту

## 0.1. З чого починалося дослідження

Початкова мета була відносно простою:

> Порівняти класичні молекулярні представлення, pretrained chemical language models і graph neural networks для прогнозування pIC50 human TRPA1 inhibitors.

Було сформовано великий ChEMBL-набір, відтворено ligand-vs-random-decoy класифікацію Mihai et al. (2020), а потім виконано scaffold-aware regression benchmark для Morgan fingerprints, RDKit descriptors, frozen ChemBERTa/MolFormer embeddings і D-MPNN.

Попередній результат: найкращі classical fingerprints і frozen transformer embeddings показали близьку якість, тоді як D-MPNN не дав стабільної переваги.

## 0.2. Що змінило напрямок роботи

Під час біологічної перевірки стало очевидно, що один стовпець `pIC50` об’єднує вимірювання, отримані різними протоколами:

- різні активатори TRPA1;
- різні концентрації та рівні активатора (наприклад EC50 або EC80);
- calcium-flux assays і patch-clamp;
- CHO, HEK293 та інші системи;
- різний порядок додавання речовин;
- різні часи контакту зі сполукою;
- різні способи fit dose-response curve;
- можливий вплив desensitization і state dependence.

Тому центральне питання було переглянуте:

> Чи обмежується прогнозування human TRPA1 inhibitory potency головним чином молекулярним представленням, чи значну частину помилки задає неоднорідність експериментальних протоколів?

## 0.3. Поточна робоча теза статті

Робота має поєднати дві лінії доказу:

1. **Representation benchmark:** складні learned representations поки не демонструють стабільної переваги над Morgan fingerprints під scaffold-aware evaluation.
2. **Assay-aware analysis:** одна й та сама молекула може отримувати систематично різні apparent pIC50 у різних варіантах протоколу.

Наразі друга лінія має сильний, але ще не остаточно інтерпретований результат: знайдено дві великі пари assays із тими самими compounds, де short і long conditions дають систематичний зсув pIC50 приблизно на 0.2–0.3 log unit після агрегації за хімічними каркасами.

## 0.4. Що вже можна сказати, а що ще не можна

### Уже підтримується даними

- Strict human TRPA1 IC50 subset містить 2196 придатних activity records, 1645 стандартизованих compounds, 97 assays і 55 documents.
- Усі 97 assay records прив’язані до target `CHEMBL6007` (human TRPA1); rat target `CHEMBL5160` у model dataset не входить.
- ChEMBL `assay_type` ненадійний: багато очевидних FLIPR і electrophysiology assays позначені як `B`.
- Ligand-vs-random-decoy класифікація є майже тривіальною і дає AUC близько 0.99.
- Під scaffold GroupKFold fingerprints і frozen transformer embeddings мають близьку попередню якість.
- Дві paired-assay comparisons показують узгоджений short-versus-long shift для тих самих compounds у двох окремих патентних документах.

### Ще не доведено

- Що саме assay heterogeneity є причиною «стелі» R² близько 0.55.
- Що paired shift однозначно спричинений саме 10- проти 90-хвилинної compound preincubation.
- Що fine-tuning програв саме через catastrophic forgetting.
- Що feature importance описує справжній загальний SAR TRPA1.
- Що один homogeneous subset обов’язково дасть кращу generalization.
- Що наша модельна робота є першою GNN або першою pIC50 regression для TRPA1.

---

# 1. Переглянуте наукове питання і гіпотези

## 1.1. Головне питання

> Наскільки передбачуваною є reported inhibitory potency human TRPA1 для нових chemical scaffolds, і яка частина невизначеності пов’язана з молекулярною структурою, а яка — з assay protocol?

## 1.2. Первинні гіпотези

### H1. Representation hypothesis

Morgan fingerprints, frozen pretrained SMILES embeddings і D-MPNN не матимуть великої та стабільної різниці після чесної scaffold-aware evaluation, якщо training set залишається відносно малим і неоднорідним.

### H2. Protocol hypothesis

Зміна assay condition може давати systematic shift у reported pIC50 для тієї самої молекули, порівнянний з істотною частиною model RMSE.

### H3. Assay-aware modelling hypothesis

Модель, що враховує protocol covariates або навчається на верифікованому протокольно однорідному subset, може перевершити structure-only model більше, ніж проста заміна RF/XGBoost на складнішу molecular architecture.

### H4. Validation hypothesis

Random split буде істотно оптимістичнішим за scaffold-, document- або time-aware evaluation.

## 1.3. Робоча назва статті

Основний варіант:

**Assay-Aware Benchmarking of Molecular Representations for Human TRPA1 Inhibitory Potency Prediction**

Альтернативний, сильніший, але допустимий лише після завершення assay-aware tests:

**Experimental Protocol Heterogeneity Limits Molecular Machine Learning for Human TRPA1 Inhibitory Potency Prediction**

---

# 2. Контекст літератури і реальна новизна

Початковий звіт переоцінював новизну. Після перевірки літератури необхідно врахувати щонайменше три близькі роботи.

## 2.1. Mihai et al., 2020

Mihai et al. побудували RF, SVM і FFNN для пошуку TRPA1 antagonists у постановці ligand-vs-decoy. Це не регресія potency і не scaffold-aware benchmark, але це вже ML на TRPA1.

**Посилання:** Mihai et al. *Artificial Intelligence Algorithms for Discovering New Active Compounds Targeting TRPA1 Pain Receptors*. AI, 2020, 1, 276–292. DOI: `10.3390/ai1020018`.

## 2.2. Gawalska et al., 2023

Gawalska et al. створили AutoML regression models для TRPA1, PDE4B і PDE8A та використали їх для virtual screening multi-target candidates. Отже, не можна заявляти «перша pIC50 regression для TRPA1».

**Посилання:** Gawalska et al. *Application of automated machine learning in the identification of multi-target-directed ligands blocking PDE4B, PDE8A, and TRPA1...* Molecular Informatics, 2023, 42, e202200214. DOI: `10.1002/minf.202200214`.

## 2.3. Castañeda-Leautaud and Amaro, 2025

У 2025 році опубліковано GNN/message-passing benchmark, який включав BindingDB TRPA1 dataset приблизно з 3020 IC50 records. У цій роботі TRPA1 був поставлений переважно як classification task із threshold 100 nM; repeated SMILES із різними IC50 були залишені. Отже, не можна заявляти «перше застосування GNN до TRPA1».

**Посилання:** Castañeda-Leautaud and Amaro. *Optimal message passing for molecular prediction is simple, attentive and spatial*. Digital Discovery, 2025, 4, 3320–3338. DOI: `10.1039/D5DD00193E`.

## 2.4. Потенційно захищена новизна нашої роботи

Після врахування цих робіт найбільш захищене формулювання таке:

> To our knowledge, this is the first target-focused, assay-aware comparison of classical fingerprints, physicochemical descriptors, frozen pretrained SMILES representations and D-MPNN for human TRPA1 IC50 regression under scaffold-aware validation, supplemented by paired analysis of protocol-associated shifts in measurements of the same compounds.

Це твердження все одно має бути повторно перевірене повним literature search перед submission.

Ключова новизна має бути не «ми використали трансформер», а:

- assay-level curation;
- paired comparison of protocol variants for the same compounds;
- evaluation of whether protocol information explains more variance than model complexity;
- strict human-only filtering;
- scaffold/document-aware validation.

---

# 3. Формування набору даних

## 3.1. Поточне джерело

Поточний API status під час assay audit:

- database: **ChEMBL 37**;
- release date: **2026-05-01**;
- extraction/audit run: липень 2026 року.

ChEMBL 37 глобально містить значні оновлення порівняно з ChEMBL 36, але це не означає, що strict TRPA1 subset обов’язково змінився.

## 3.2. Target і строгий фільтр

Target:

```text
target_chembl_id = CHEMBL6007
target = human TRPA1
```

Strict activity filter:

```text
standard_type = IC50
standard_relation = "="
pchembl_value is not null
target_chembl_id = CHEMBL6007
```

Поточні числа:

- 2203 raw activities після API filter;
- 2196 usable activities після structure standardization і видалення непридатних structure rows;
- 1652 raw molecule IDs;
- 1645 unique standardized compounds;
- 97 assays;
- 55 documents;
- document years 2010–2025.

## 3.3. Номенклатура

У старій версії всі IC50 records називалися «антагоністами». Це надто сильне припущення.

`IC50` у heterogeneous functional assays може відображати:

- справжній allosteric antagonist;
- pore або state-dependent blocker;
- indirect inhibition;
- agonist-induced desensitization;
- assay interference;
- compound toxicity;
- інші protocol-dependent effects.

Тому в основному тексті слід використовувати:

> compounds with reported inhibitory IC50 values against human TRPA1

або коротше:

> human TRPA1 inhibitors in ChEMBL

Термін **functional antagonist** допустимий лише для вручну перевірених assays/compounds із відповідним протоколом.

## 3.4. Structure standardization

Пайплайн:

1. RDKit parsing через `Chem.MolFromSmiles`.
2. Видалення солей через `SaltRemover`.
3. Якщо є кілька фрагментів — вибір найбільшого за кількістю heavy atoms.
4. Canonical isomeric SMILES.
5. InChIKey.
6. Дедуплікація за InChIKey.
7. Основна агрегація target: median pChEMBL.

Для кожного compound зберігалися:

- `pchembl_median`;
- `pchembl_min`;
- `pchembl_max`;
- `pchembl_std`;
- `n_measurements`;
- year range;
- molecule and assay identifiers.

### Важлива незавершена перевірка

Потрібно переконатися, що stereoisomers не колапсують помилково на якомусь етапі та що всі model representations отримують однаковий stereochemical input. Це особливо важливо для TRPA1 patent series, де окремі stereoisomers можуть мати різну активність.

## 3.5. Human і rat не змішуються

TRPA1 має критично важливу species-dependent pharmacology. Для окремих chemotypes одна речовина може бути antagonist у human TRPA1 і agonist/partial agonist або слабкоактивною в rat TRPA1.

Перевірка `trpa1_assay_audit.csv` показала:

- усі 97 records: `target_chembl_id = CHEMBL6007`;
- 87 records: `assay_organism = Homo sapiens`;
- у 10 `assay_organism` не заповнене, але target залишається human;
- rat target `CHEMBL5160` у dataset відсутній.

Деякі generic assay descriptions згадують паралельне тестування human і rat, але це не означає, що rat activity values увійшли в human dataset. У protocol decomposition такі згадки зберігаються лише як warning, а record species жорстко задається target ID.

### Обов’язковий engineering fix

Наступна версія pipeline має не просто припускати human target, а падати з помилкою, якщо у вході присутній будь-який target ID, відмінний від `CHEMBL6007`.

## 3.6. Порівняння ChEMBL 36 і ChEMBL 37

Було проведено повторне API extraction та зіставлення strict subset із раніше збереженим файлом, який вважався ChEMBL 36 dataset.

Результат:

- 1645 unique InChIKeys у saved dataset;
- 1645 unique InChIKeys у ChEMBL 37 pull;
- overlap: 1645;
- new compounds: 0;
- lost compounds: 0;
- total usable measurements: 2196 у кожному;
- per-compound median/min/max/std/count/year aggregates не змінилися, окрім floating-point noise близько `1e-15`.

Коректний висновок:

> No new compounds or activity-level aggregate changes were identified between the saved dataset and ChEMBL 37 for the strict human TRPA1 IC50 subset with exact relation and available pChEMBL values.

Некоректний висновок:

> У ChEMBL 37 взагалі не додано нових TRPA1 data.

Інші activity types, relations, species або records без pChEMBL могли змінитися.

### Provenance caveat

Старий CSV не містив embedded database-version metadata. Тому абсолютне твердження, що він гарантовано був сформований саме до релізу ChEMBL 37, залежить від зовнішнього timestamp/notebook history. У майбутньому status JSON, extraction timestamp і input hash зберігаються автоматично.

---

# 4. Початковий ML-бенчмарк

## 4.1. Molecular representations

| Representation | Dimension | Downstream model / details |
|---|---:|---|
| Morgan ECFP4 | 2048 bits | RF або XGBoost |
| RDKit physicochemical descriptors | 15 | RF |
| ChemBERTa-77M-MTR | 384 | CLS або mean pooling; XGBoost |
| MolFormer-XL-both-10pct | 768 | CLS або mean pooling; XGBoost |
| Molecular graph | learned | D-MPNN |
| Graph + RDKit descriptors | learned + 15 | D-MPNN + descriptors |

15 RDKit descriptors:

- MolWt;
- MolLogP;
- MolMR;
- TPSA;
- HBA/HBD;
- rotatable bonds;
- aromatic/aliphatic/saturated ring counts;
- FractionCSP3;
- heavy atom count;
- heteroatom count;
- LabuteASA.

### Заплановане виправлення

Morgan fingerprint має бути перерахований з явним `includeChirality=True`, а стара і нова версії слід порівняти як sensitivity analysis.

## 4.2. Validation design

Використано:

- deterministic Bemis–Murcko scaffold split: train 1324 / validation 160 / test 161;
- 5-fold `GroupKFold` із scaffold як group;
- 544 unique scaffolds;
- regression metrics: R², RMSE, Spearman;
- secondary threshold analysis: AUC і MCC при `pChEMBL >= 7`.

Scaffold CV був обраний як більш реалістична перевірка переносу на нові chemotypes, ніж random split.

## 4.3. Попередні результати scaffold GroupKFold

| Model | R² mean±SD | RMSE mean±SD | Spearman | AUC |
|---|---:|---:|---:|---:|
| RF + Morgan | **0.547±0.112** | **0.618±0.050** | не збережено через bug | 0.856±0.060 |
| XGB + ChemBERTa CLS | 0.539±0.084 | 0.633±0.057 | 0.725±0.064 | **0.866±0.043** |
| XGB + Morgan | 0.539±0.106 | 0.625±0.046 | не збережено через bug | 0.860±0.052 |
| XGB + MolFormer mean | 0.529±0.085 | 0.639±0.042 | 0.703±0.066 | 0.857±0.045 |
| D-MPNN + RDKit | 0.491±0.092 | 0.663±0.040 | 0.677±0.086 | 0.843±0.049 |
| RF + RDKit-15 | 0.478±0.163 | 0.664±0.053 | 0.695±0.089 | 0.848±0.051 |
| D-MPNN graph-only | 0.410±0.107 | 0.714±0.031 | 0.597±0.089 | 0.794±0.050 |

Попередня інтерпретація:

- Morgan fingerprints і frozen language-model embeddings практично близькі.
- D-MPNN не показав переваги.
- Fold-to-fold variability велика.
- Без paired statistics не можна стверджувати, що перші чотири моделі реально відрізняються.

## 4.4. Single-split results

Single scaffold test set із 161 compounds був складнішим за середній CV fold.

R²:

- RF-Morgan: 0.232;
- XGB-Morgan: 0.257;
- LightGBM-Morgan: 0.255;
- RF-RDKit15: 0.363;
- D-MPNN graph: 0.176;
- D-MPNN+RDKit: 0.292;
- XGB-ChemBERTa CLS: 0.344;
- XGB-ChemBERTa CLS + RDKit: 0.395;
- XGB-MolFormer mean: 0.418;
- ChemBERTa fine-tuning: 0.247;
- MolFormer fine-tuning: 0.193.

Ці single-split numbers не мають бути головною таблицею статті.

## 4.5. Fine-tuning

End-to-end fine-tuning обох pretrained models погіршив single-split performance порівняно з frozen embeddings.

Стара інтерпретація називала це «catastrophic forgetting». Це не доведено. Альтернативні пояснення:

- недостатній або несправедливий tuning;
- learning-rate schedule;
- model-selection variance;
- scaffold-specific overfitting;
- нестабільність half precision;
- невдалий validation split;
- small effective sample size.

Коректне формулювання:

> Under the tested training protocols, end-to-end fine-tuning did not outperform frozen feature extraction.

## 4.6. Відтворення Mihai et al.

Dataset:

- 1645 TRPA1 IC50 compounds як positive class;
- 560 random ChEMBL decoys;
- random 80/20 split;
- Morgan ECFP4 замість proprietary MNA descriptors.

Результати:

| Model | Accuracy | CV AUC |
|---|---:|---:|
| RF | 98.64% | 0.9982 |
| SVM | 98.64% | 0.9968 |
| FFNN | 99.09% | 0.9983 |

Це підтримує висновок, що ligand-vs-random-decoy classification є надто легкою постановкою. Decoys не були property-matched, тому replication не повинна використовуватися як доказ практичної virtual-screening якості.

## 4.7. Відомі проблеми ML-бенчмарку

### Критичні

1. **Morgan Spearman bug:** значення не були коректно перенесені до фінальної таблиці.
2. **Нерівний tuning budget:** різні model families отримали нерівноцінну оптимізацію.
3. **D-MPNN early stopping:** не спрацював у CV; усі folds дійшли до max epochs.
4. **Немає єдиного OOF prediction table** для всіх моделей і paired bootstrap comparison.
5. **Немає document-held-out або temporal stress test** у завершеному вигляді.
6. **Assay context ігнорувався** в target construction.

### Суттєві

7. Decoys не property-matched.
8. Threshold classification є secondary і залежить від довільного cutoff.
9. Fine-tuning conclusions зроблені з одного split.
10. D-MPNN і classical models не були однаково оптимізовані.
11. Acyclic molecules можуть потребувати спеціальної обробки в Murcko grouping.
12. Не проведена systematic uncertainty analysis для multiple measurements.

---

# 5. Чому assay context є біологічно важливим саме для TRPA1

TRPA1 не є простою мішенню, для якої будь-яке IC50 означає одну й ту саму фізичну величину.

## 5.1. Різні способи активації

Електрофільні agonists, такі як AITC, можуть активувати TRPA1 через reversible covalent modification intracellular cysteines. Неелектрофільні ligands можуть взаємодіяти з іншими sites і conformational states.

## 5.2. Calcium dependence і desensitization

TRPA1 є Ca²⁺-permeable channel. Calcium може спочатку potentiating, а потім desensitizing/inactivating response. У calcium imaging observed inhibition може бути сумою:

- direct channel antagonism;
- state-dependent block;
- agonist-induced desensitization;
- depletion/signalling effects;
- cell toxicity або fluorescence interference.

## 5.3. Різні antagonist mechanisms

Відомі antagonists можуть взаємодіяти з різними binding regions або channel states. Тому IC50 проти AITC не обов’язково еквівалентне IC50 проти cinnamaldehyde або іншого agonist.

## 5.4. Species dependence

Human і rat TRPA1 pharmacology може бути якісно різною. Через це model dataset має залишатися strict human-only; species не можна змішувати як звичайні technical replicates.

## 5.5. Наслідок для ML target

Наш target слід визначати так:

> reported inhibitory potency in heterogeneous human TRPA1 assays

а не:

> universal binding affinity of TRPA1 antagonists.

Це не робить dataset непридатним. Але assay protocol стає частиною label definition.

---

# 6. Assay audit

## 6.1. Масштаб

Strict dataset містить:

- 97 assays;
- 55 documents;
- 2196 activities;
- 1645 compounds.

Це зробило assay-level manual audit реалістичним: потрібно розібрати 97 assay descriptions, а не вручну читати 1645 compounds.

## 6.2. Базовий assay audit

Було створено standalone pipeline, який:

- витягує ChEMBL assay metadata;
- кешує API responses;
- зберігає status/version і hashes;
- підраховує compounds та measurements за assay;
- ставить keyword flags;
- створює порожні manual annotation columns;
- зберігає вже внесені manual annotations при повторному запуску.

Виявлена важлива помилка ChEMBL metadata:

- 77 із 97 assays мали `assay_type = B`;
- серед них були очевидні FLIPR calcium та whole-cell patch-clamp experiments.

Тому `assay_type` не використовується як inclusion gate.

## 6.3. Protocol decomposition

Другий pipeline розбирає assay descriptions на параметри:

- modality;
- challenge agonist;
- agonist level;
- agonist concentration;
- cell line;
- record species;
- construct;
- dye loading;
- plate equilibration;
- compound preincubation candidates;
- application order;
- Hill coefficient;
- evidence fragment;
- confidence/ambiguity;
- review priority.

Поточне coverage:

| Field | Recognized assays |
|---|---:|
| assay modality | 74/97 |
| challenge agonist | 50/97 |
| agonist level | 7/97 |
| cell line | 49/97 |
| single unambiguous compound preincubation time | 4/97 |
| application order | 8/97 |

Низьке coverage не означає, що параметр відсутній у першоджерелі. Часто ChEMBL description просто неповний.

## 6.4. Review priority

Поточна автоматична triage:

- critical: 4 assays;
- high: 4 assays;
- medium: 89 assays;
- low: 0.

`critical` assays — дві великі paired comparisons, які формують центральний result.

## 6.5. Protocol signatures

Сформовано 22 **protocol signatures**.

Це означає:

> groups equal by currently extracted fields

а не:

> 22 підтверджені однорідні протоколи.

Кілька найбільших signatures:

| Unique compounds | Assays | Documents | Signature quality | Extracted signature |
|---:|---:|---:|---|---|
| 275 | 2 | 1 | well specified | calcium fluorescence; cinnamaldehyde; EC50; CHO; human |
| 291 | 7 | 4 | low information | electrophysiology; human; інші поля NA |
| 150 | 2 | 1 | well specified | calcium fluorescence; cinnamaldehyde; EC80; CHO; human |
| 214 | 2 | 2 | well specified | calcium fluorescence; cinnamaldehyde; EC80; CHO; human |
| 144 | 8 | 8 | partially specified | calcium fluorescence; AITC; human |
| 106 | 9 | 8 | partially specified | calcium fluorescence; HEK293; human |
| 108 | 1 | 1 | well specified | calcium fluorescence; cinnamaldehyde; EC50; CHO; human |

`NA == NA` не є доказом однакового protocol. Тому `verified_protocol_family` залишається `False` до manual/source review.

---

# 7. Paired-assay analysis

## 7.1. Ідея

Найсильніший спосіб показати protocol-associated variation — порівнювати **ті самі compounds** у двох assays з одного документа.

Це контролює molecular identity: структура не змінюється, змінюється experimental condition.

Search criteria:

- same ChEMBL document;
- at least 20 shared InChIKeys;
- assays мають одну extracted difference або time-related difference;
- pChEMBL агрегується всередині assay за InChIKey;
- delta = assay B − assay A;
- додатково compounds агрегуються за Bemis–Murcko scaffold.

## 7.2. Чотири critical assays

| Document | Patent | Assay | Compounds in assay | Agonist level | Current ChEMBL title time |
|---|---|---|---:|---|---:|
| CHEMBL5727828 | US-10710994-B2 | CHEMBL5735474 | 131 | EC80 | 15 min |
| CHEMBL5727828 | US-10710994-B2 | CHEMBL5735475 | 150 | EC80 | 90 min |
| CHEMBL5727829 | US-10711004-B2 | CHEMBL5735476 | 275 | EC50 | 15 min |
| CHEMBL5727829 | US-10711004-B2 | CHEMBL5735477 | 169 | EC50 | 90 min |

Патентні документи мають однакову назву *Oxadiazole transient receptor potential channel inhibitors*, але це два окремі патентні документи та споріднені medicinal-chemistry series, а не дві незалежні лабораторії.

Patent family information:

- US-10710994-B2 пов’язаний із WO2019182925A1.
- US-10711004-B2 пов’язаний із ранньою oxadiazole family, включно з WO2018162607A1 / EP3592739A1.

## 7.3. Що встановлено з патентного методу

Для family WO2019182925A1 method описує:

1. CHO cells expressing human TRPA1.
2. BD calcium indicator dye: 1 hour at 37 °C.
3. 15 minutes at room temperature після dye loading.
4. Додавання test compounds.
5. Incubation with compounds for **10 minutes or 90 minutes**.
6. Додавання приблизно EC80 cinnamaldehyde.
7. Measurement of blocked cinnamaldehyde-induced calcium influx.
8. Hill fit із fixed coefficient `n = 1.5`.

Отже:

- 15 min у загальному методі — plate/dye equilibration;
- 10/90 min — compound incubation before agonist.

## 7.4. Невирішена проблема mapping

ChEMBL assay descriptions для чотирьох records внутрішньо суперечливі:

- title і phrase після dye loading змінені на 15 або 90 min;
- далі в обох descriptions залишається `compounds for 10 minutes or 90 minutes`.

Тому ще не доведено mapping:

```text
CHEMBL5735474 = 10 min compound preincubation?
CHEMBL5735475 = 90 min compound preincubation?
CHEMBL5735476 = 10 min compound preincubation?
CHEMBL5735477 = 90 min compound preincubation?
```

Найімовірніше, ChEMBL curator розділив short і long compound conditions, але невдало сформував description. Проте до прямої перевірки patent columns/PubChem mapping у статті слід писати **short versus long assay condition**.

## 7.5. Поточні paired results

Після приведення напряму до:

```text
long condition − short condition
```

отримано:

| Document | Shared compounds | Murcko scaffolds | Compound median Δ | Scaffold-median Δ | Direction consistency |
|---|---:|---:|---:|---:|---:|
| CHEMBL5727829 | 169 | 70 | +0.290 | +0.200 | 92% |
| CHEMBL5727828 | 131 | 51 | +0.330 | +0.290 | 95% |

На compound level mean shifts:

- CHEMBL5727829: +0.317, t-based CI +0.278 to +0.357;
- CHEMBL5727828: +0.404 у normalized long-minus-short direction; original console pair was printed in reverse order as −0.404.

Приблизна concentration ratio:

- `10^0.20 ≈ 1.58`;
- `10^0.29 ≈ 1.95`;
- `10^0.33 ≈ 2.14`.

Отже, long condition робить compounds apparent approximately 1.6–2.1-fold more potent.

## 7.6. Статистична caveat у поточному script

Поточний script показує:

- point estimate: median of per-scaffold median deltas;
- bootstrap CI: resampling scaffold clusters, потім concatenation усіх compound deltas і median.

Ці дві statistics не є повністю ідентичними: великі scaffold groups все ще мають більшу вагу в bootstrap sample.

До submission потрібно перерахувати CI саме для **median of scaffold medians**:

1. обчислити median delta для кожного scaffold;
2. bootstrap-resample ці scaffold medians із поверненням;
3. брати median кожного bootstrap sample;
4. percentile 95% CI.

Тому поточні scaffold bootstrap CIs і p-values слід вважати provisional, хоча напрямок і magnitude result навряд чи зникнуть.

## 7.7. Інші caveats paired analysis

- Два документи походять із спорідненої pharmaceutical program; це replication across documents/series, але не independent laboratory replication.
- Murcko scaffold не повністю моделює залежність усередині medicinal-chemistry series.
- Acyclic molecules не слід об’єднувати в один empty scaffold.
- Потрібно перевірити exact values/columns для 5–10 compounds у кожній pair.
- Необхідно підтвердити, що інші assay parameters справді однакові.
- Fixed Hill coefficient 1.5 є додатковою protocol variable і має бути врахований при comparison з іншими documents.

## 7.8. Коректне поточне формулювання

> Two paired sets of human TRPA1 FLIPR measurements showed a reproducible systematic shift in apparent pIC50 between short and long assay conditions for the same compounds.

Після підтвердження mapping можна буде написати:

> Increasing compound preincubation from 10 to 90 minutes was associated with a systematic increase in apparent pIC50.

До mapping слово `associated` безпечніше за `caused`.

---

# 8. Переглянутий дизайн статті

## 8.1. Experiment A: replication of easy ligand-vs-decoy classification

Мета:

- відтворити Mihai et al.;
- показати, що near-perfect metrics не означають accurate potency prediction;
- винести в короткий secondary section або supplement.

## 8.2. Experiment B: representation benchmark

Models:

1. RF + Morgan;
2. XGB + Morgan;
3. RF/XGB + RDKit descriptors;
4. XGB + frozen ChemBERTa;
5. XGB + frozen MolFormer;
6. D-MPNN;
7. D-MPNN + descriptors, опційно.

Primary validation:

- 5-fold scaffold GroupKFold;
- identical folds;
- identical OOF table;
- paired scaffold bootstrap between models;
- modest, declared and reasonably balanced tuning budget.

Primary metrics:

- RMSE;
- MAE;
- R²;
- Spearman.

AUC/MCC — secondary only.

## 8.3. Experiment C: protocol decomposition

Для 97 assays:

- modality;
- agonist;
- agonist level/concentration;
- cell line;
- construct;
- application order;
- preincubation;
- readout;
- fit method;
- evidence and confidence.

Output should distinguish:

- extracted signature;
- manually verified protocol family;
- unknown/ambiguous fields.

## 8.4. Experiment D: paired protocol contrasts

Primary paired contrasts:

- CHEMBL5735474 vs CHEMBL5735475;
- CHEMBL5735476 vs CHEMBL5735477.

Primary estimate:

- median of per-scaffold paired deltas;
- scaffold-level bootstrap CI;
- direction consistency;
- Spearman between conditions;
- Bland–Altman style visualization;
- per-document result, not only pooled result.

## 8.5. Experiment E: assay-aware prediction

Порівняти:

### Model 1

```text
molecular structure only
```

### Model 2

```text
molecular structure + protocol covariates
```

Можливі protocol covariates:

- modality;
- agonist;
- EC level;
- cell line;
- compound incubation;
- document/patent;
- fit method.

Критична вимога: protocol covariates мають бути available for test data; не можна створити leakage через assay/document identity.

## 8.6. Experiment F: verified homogeneous subset

Вибрати найбільший protocol subset, верифікований за першоджерелами.

Не просто `cinnamaldehyde + calcium`, а точна комбінація:

- human target;
- construct;
- cell line;
- FLIPR/calcium readout;
- agonist;
- agonist level/concentration;
- application order;
- compound preincubation;
- fit method.

Порівнювати homogeneous subset не напряму з усіма 1645 compounds, а з repeated matched broad subsamples:

- same sample size;
- approximately matched scaffold count/diversity;
- matched pIC50 distribution;
- same model/hyperparameters;
- same validation principle.

Інакше higher R² може пояснюватися narrow chemical series, а не cleaner labels.

## 8.7. Experiment G: document-aware stress test

Якщо protocol зустрічається в кількох documents:

- train on some documents;
- test on held-out document;
- no exact structures across split;
- report chemical similarity to training set.

Це важливо, бо patent-specific medicinal-chemistry series можуть робити random/scaffold evaluation оптимістичною.

---

# 9. Що завершено і що залишилося

## 9.1. Завершено

- Human TRPA1 target validation.
- Strict IC50 extraction і standardization.
- ChEMBL status/provenance capture.
- ChEMBL 36-like saved dataset vs ChEMBL 37 strict-subset comparison.
- Mihai replication.
- Preliminary scaffold CV benchmark.
- Frozen transformer embeddings.
- D-MPNN and D-MPNN+descriptor runs.
- Preliminary fine-tuning attempts.
- Assay metadata audit for all 97 assays.
- Protocol decomposition with evidence/confidence fields.
- Human/rat parsing correction.
- Identification of 22 extracted protocol signatures.
- Identification of two large paired assay contrasts.
- Preliminary scaffold-aware paired effect estimates.

## 9.2. Обов’язково до submission

1. Перерахувати Morgan Spearman.
2. Зберегти OOF predictions для всіх models.
3. Вирівняти tuning protocol або чітко обмежити claims.
4. Полагодити D-MPNN early stopping.
5. Перерахувати Morgan with chirality.
6. Зробити paired model comparison with scaffold bootstrap.
7. Виправити bootstrap CI для median of scaffold medians.
8. Обробити acyclic compounds без artificial empty-scaffold cluster.
9. Додати hard assertion `target_chembl_id == CHEMBL6007`.
10. Включити construct і agonist concentration у verified protocol comparison.
11. Перевірити exact mapping чотирьох critical assays до 10/90-minute patent columns.
12. Ручно перевірити 4 critical і 4 high-priority assays.
13. Створити хоча б один assay-aware modelling або matched homogeneous-subset analysis.
14. Провести document-aware stress test, якщо структура даних дозволяє.

## 9.3. Бажано, але не обов’язково для першої статті

- Повна ручна курація всіх 97 assays.
- Temporal split за first appearance.
- Additional pretrained model.
- Docking/MD.
- Virtual screening.
- Prospective experimental validation.

Ці пункти не повинні затримувати першу methodology-focused publication, якщо core analyses завершені.

---

# 10. Feature importance і SAR: переглянута позиція

Стара версія звіту робила сильний SAR-висновок із RF feature importance:

- TPSA і heteroatoms positive;
- LogP negative;
- HBD negative;
- aromatic rings positive.

Це не можна називати загальним TRPA1 SAR, тому що correlations можуть бути зумовлені:

- assay/document composition;
- patent-series composition;
- scaffold family;
- target-range restriction;
- medicinal-chemistry optimization history.

Partial correlation controlling only for publication year не усуває assay, document і scaffold confounding.

Коректне використання:

> exploratory dataset-level descriptor associations

Вони можуть увійти в supplement, але не мають бути центральним mechanistic result без within-series або multilevel analysis.

---

# 11. Поточний очікуваний висновок статті

Найбільш захищений варіант, якщо наступні аналізи підтвердять поточну картину:

> Classical fingerprints and frozen pretrained molecular representations showed similar scaffold-level performance for human TRPA1 inhibitory potency prediction, while increased architectural complexity did not provide a consistent advantage. Assay-level decomposition revealed substantial protocol heterogeneity, and paired measurements of the same compounds under short and long assay conditions showed reproducible systematic shifts in apparent pIC50. These findings indicate that experimental context should be modelled explicitly and may be at least as important as molecular representation choice for TRPA1 bioactivity prediction.

Чого не слід писати без додаткового доказу:

> Assay heterogeneity definitively sets an R² ceiling of 0.55.

---

# 12. Пропонована структура рукопису

## Introduction

1. TRPA1 як therapeutic target.
2. Складність ligand pharmacology, calcium modulation, desensitization і species dependence.
3. Попередні ML studies: Mihai 2020, Gawalska 2023, Castañeda-Leautaud 2025.
4. Проблема heterogeneous bioactivity labels.
5. Мета: representation benchmark plus assay-aware analysis.

## Methods

1. ChEMBL 37 extraction.
2. Strict human filtering.
3. Structure standardization and deduplication.
4. Assay metadata audit.
5. Protocol decomposition and manual verification.
6. Molecular representations and models.
7. Scaffold/document validation.
8. Paired assay statistics.
9. Protocol-aware modelling.
10. Reproducibility and provenance.

## Results

1. Dataset and assay heterogeneity.
2. Replication of ligand-vs-decoy classification.
3. Scaffold-aware regression benchmark.
4. Protocol signatures.
5. Paired short-long condition shifts.
6. Assay-aware/homogeneous-subset modelling.
7. Applicability domain.

## Discussion

1. Model complexity versus data-generating process.
2. Biological meaning of protocol-associated shifts.
3. Limitations of ChEMBL assay descriptions.
4. Species strictness.
5. Patent-series dependence.
6. Implications for future TRPA1 screening.

---

# 13. Мінімальний набір figures і tables

## Figure 1

Data flow:

```text
ChEMBL 37
→ strict human TRPA1 IC50
→ structure standardization
→ assay audit
→ protocol decomposition
→ molecular ML + paired assay analysis
```

## Figure 2

Assay landscape:

- compounds per assay;
- documents;
- modalities;
- agonists;
- protocol-field missingness.

## Figure 3

Scaffold CV model comparison with paired confidence intervals.

## Figure 4

Paired short/long conditions:

- identity plot;
- delta distribution;
- scaffold median deltas;
- separate colors for two documents.

## Figure 5

Prediction error or performance versus:

- protocol information;
- homogeneous vs matched broad subsets;
- train similarity.

## Table 1

Dataset and protocol audit characteristics.

## Table 2

Model performance.

## Table 3

Four critical assays and source-verified conditions.

## Table 4

Paired effects by document.

---

# 14. Files and reproducibility assets

Основні наявні файли:

- `trpa1_antagonists.csv`
- `trpa1_current_api_raw.csv`
- `trpa1_assay_audit.csv`
- `trpa1_assay_audit_metadata.json`
- `protocol_decomposition.csv`
- `protocol_families.csv`  
  Назва залишилася legacy; логічно це `protocol_signatures.csv`.
- `paired_contrasts.csv`
- `protocol_metadata.json`
- scripts for extraction, assay audit and protocol decomposition
- model notebooks/scripts for Morgan, descriptors, transformers and Chemprop

До reproducibility package потрібно додати:

- exact environment files;
- model hyperparameter tables;
- fold assignment file;
- all OOF predictions;
- raw ChEMBL status JSON;
- hashes усіх input files;
- manual protocol override table;
- source evidence for four critical assays.

---

# 15. Питання до рецензента

1. Чи є revised framing достатньо новим після Gawalska 2023 і Castañeda-Leautaud 2025?
2. Чи paired short-long result є достатньо сильним центральним result, якщо exact 10/90 mapping буде підтверджено?
3. Чи потрібна повна курація всіх 97 assays, чи достатньо 8 priority assays плюс найбільші protocol signatures?
4. Який найкращий design для protocol-aware model без leakage?
5. Чи доречно використовувати median of scaffold medians як primary paired estimate?
6. Чи краще cluster bootstrap за patent/scaffold hierarchy або mixed-effects model?
7. Як правильно оцінити effective sample size у patent medicinal-chemistry series?
8. Чи достатній matched-subsample design для порівняння homogeneous і broad datasets?
9. Чи варто залишити fine-tuning section, чи перенести його повністю в supplement?
10. Чи варто вилучити Mihai replication з main text?
11. Чи є document-held-out evaluation обов’язковою?
12. Чи має сенс terminology `reported inhibitory potency`, чи потрібне ще обережніше формулювання?
13. Який мінімальний набір виправлень потрібен для українського фахового журналу?
14. Що додатково потрібно для міжнародного журналу рівня *Molecular Informatics*, *AI* або подібного?

---

# 16. Ключові references

1. Mihai DP et al. *Artificial Intelligence Algorithms for Discovering New Active Compounds Targeting TRPA1 Pain Receptors*. AI. 2020;1:276–292. DOI: `10.3390/ai1020018`.
2. Gawalska A et al. *Application of automated machine learning in the identification of multi-target-directed ligands blocking PDE4B, PDE8A, and TRPA1 with potential use in the treatment of asthma and COPD*. Molecular Informatics. 2023;42:e202200214. DOI: `10.1002/minf.202200214`.
3. Castañeda-Leautaud AC, Amaro RE. *Optimal message passing for molecular prediction is simple, attentive and spatial*. Digital Discovery. 2025;4:3320–3338. DOI: `10.1039/D5DD00193E`.
4. Talavera K et al. *Mammalian Transient Receptor Potential TRPA1 Channels: From Structure to Disease*. Physiological Reviews. 2020;100:725–803. DOI: `10.1152/physrev.00005.2019`.
5. Meents JE, Ciotu CI, Fischer MJM. *TRPA1: a molecular view*. Journal of Neurophysiology. 2019;121:427–443. DOI: `10.1152/jn.00524.2018`.
6. Zhao J et al. *Irritant-evoked activation and calcium modulation of the TRPA1 receptor*. Nature. 2020;585:141–145. DOI: `10.1038/s41586-020-2480-9`.
7. Patent US-10710994-B2. *Oxadiazole transient receptor potential channel inhibitors*. 2020. Related PCT publication: WO2019182925A1.
8. Patent US-10711004-B2. *Oxadiazole transient receptor potential channel inhibitors*. 2020. Related family includes WO2018162607A1 / EP3592739A1.
9. ChEMBL 37 release information, EMBL-EBI ChEMBL team, May 2026.

---

# 17. One-paragraph summary for reviewer

This project began as a molecular representation benchmark for human TRPA1 IC50 prediction. A strict ChEMBL 37 subset produced 2196 usable activity records for 1645 standardized compounds across 97 assays. Preliminary scaffold GroupKFold results showed similar performance for Morgan fingerprints and frozen ChemBERTa/MolFormer representations (R² approximately 0.53–0.55), with lower performance for the tested D-MPNN protocols. Subsequent biological audit revealed substantial assay heterogeneity. Protocol decomposition identified two large paired assay contrasts from two separate oxadiazole patent documents, involving the same compounds under short and long FLIPR conditions. Long conditions yielded apparent pIC50 increases of approximately 0.2–0.3 log unit at the scaffold-aggregated level, with 92–95% directional consistency. However, exact mapping of the four ChEMBL assay IDs to 10- versus 90-minute compound preincubation requires final source verification, and the current scaffold bootstrap CI implementation needs correction. The revised proposed article therefore combines a scaffold-aware representation benchmark with assay-level decomposition and paired protocol analysis, testing whether experimental context explains prediction error more strongly than molecular model complexity.
