# План статті: TRPA1 representation benchmark (ревізія 2)

**Статус:** робочий документ. Усі числа з фактично виконаних обчислень. Позначки `[Є]`, `[ТРЕБА]`, `[БАГ]` вказують реальний стан.

**Що змінено проти ревізії 1:** виправлено змішування представлення з алгоритмом; знижено впевненість формулювань до рівня наявних доказів; побічні експерименти винесено в Supplementary; переоцінено терміни.

---

## 0. Позиціонування

### Що стаття стверджує

> Порівняння класичних та навчених молекулярних представлень при однаковому регресорі для прогнозування **агрегованої літературної** інгібувальної активності людського TRPA1, за scaffold-aware валідації, на відтворювано сформованому відкритому наборі.

### Цільова змінна — називати чесно

Після дедуплікації мітка = **медіана pChEMBL по всіх вимірюваннях сполуки з різних assays**. Це не афінність, не потужність у конкретному протоколі, не механізм.

**Формулювання:** *prediction of aggregated reported human TRPA1 inhibitory potency from molecular structure.*

Слово **aggregated** обов'язкове. Без нього створюється враження що pIC50 — однозначна фізична властивість молекули. У нашому наборі 24% сполук виміряні в ≥2 assays, 18% мають мітку усереднену між двома протокольними умовами.

### Заборонені формулювання

| Заборонено | Ким спростовано / чому |
|---|---|
| «перше ML для TRPA1» | Mihai et al. 2020 |
| «перша регресія pIC50 для TRPA1» | Gawalska et al. 2023 |
| «перша GNN для TRPA1» | Castañeda-Leautaud & Amaro 2025 |
| «перший pretrained CLM для TRP-каналу» | Lee et al. 2023 (TRPV1) |
| «representation важливіше за algorithm» як нове | Adamczyk & Czech; Sadeghi et al. 2024 |
| «ми виявили що змішані IC50 шумні» | Landrum & Riniker 2024 |
| «складність моделі не окуповується» як нове | Schiebroek et al. 2026 |
| «**statistically** unresolved tier» | парний тест ще не зроблено |
| «**curated** dataset» | assay-контекст вручну не анотовано |
| «fine-tuning систематично програє» | протоколи не зрівняні |
| «catastrophic forgetting» | механізм не вимірювався |
| «TPSA підвищує потужність» | лише асоціація, мультиколінеарність |
| «виявлено помилку ChEMBL» | mapping patent→assay не підтверджено таблицями |

### Новизна: 2-3/10

Target-specific comparative study. Внесок — відкритий набір з provenance, відтворюваний бенчмарк, характеристика гетерогенності. Не методологічний прорив.

### Журнал

Український фаховий або скромне міжнародне видання.

---

## 1. Назва

**Benchmarking Classical and Learned Molecular Representations for Human TRPA1 Inhibitory Potency Prediction**

Не використовувати «A Curated Dataset» — набір відтворювано зібраний і стандартизований, але не вручну курований за протоколами.

---

## 2. Головний дизайн — сітка представлення × алгоритм

**Це виправлення найбільшої діри ревізії 1.** Порівнювати «Morgan+RF проти MolFormer+XGB» некоректно: одночасно змінюються представлення і алгоритм, ефекти нерозділені.

Чесний дизайн — **повна сітка на однакових фолдах**:

| Представлення | RF | XGBoost |
|---|---|---|
| Morgan ECFP4 (2048) | 0.547±0.112 `[Є]` | 0.539±0.106 `[Є]` |
| RDKit-15 дескриптори | 0.478±0.163 `[Є]` | **`[ТРЕБА]`** |
| ChemBERTa-CLS (384) | **`[ТРЕБА]`** | 0.539±0.084 `[Є]` |
| MolFormer-Mean (768) | **`[ТРЕБА]`** | 0.529±0.085 `[Є]` |

Бракує **3 клітинок з 8**. Усі дешеві — embeddings збережені в `embeddings_all.npz`, треба лише прогнати ті самі фолди.

**Окремо — end-to-end пайплайни** (не представлення, а цілі методи):
- D-MPNN (граф) 0.410±0.107 `[Є, потребує мультисіду]`
- D-MPNN + RDKit-15 0.491±0.092 `[Є, потребує мультисіду]`

Тоді можна казати: *"With the regressor held fixed, pretrained embeddings did not outperform Morgan fingerprints."* Зараз так казати не можна.

---

## 3. Abstract — скелет

1. TRPA1 як мішень; складна фармакологія.
2. Наявні ML-роботи (Mihai 2020, Gawalska 2023, Castañeda-Leautaud 2025) не порівнювали класи представлень систематично.
3. 1645 сполук / 2196 вимірювань / 97 assays з ChEMBL 37; чотири представлення × два регресори на однакових 5-кратних GroupKFold розбиттях за каркасами.
4. Головний результат benchmark.
5. Набір і код відкриті.

**Не включати в abstract:** fine-tuning, pooling, descriptor importance, протокольний case study, реплікацію Mihai.

---

## 4. Methods

### 4.1. Формування набору `[Є]`

- **ChEMBL 37**, реліз **2026-05-01**, через `chembl_webresource_client`. Версію підтверджено прямим запитом `status.json` (24 527 044 активностей).
- Мішень **CHEMBL6007** (human TRPA1).
- Фільтр на етапі запиту: `standard_type='IC50'`, `standard_relation='='`, `pchembl_value__isnull=False`. Цензуровані значення виключені за побудовою.
- **2203** сирі активності → **2196** придатних після стандартизації → **1645** унікальних сполук, **97** assays, **55** документів, роки **2010-2025**.
- **Уточнити в тексті:** з 97 assays **87** мають `assay_organism = Homo sapiens`, **10** не мають анотації організму; усі відібрані за `target_chembl_id = CHEMBL6007`.

**Стандартизація:** RDKit — парсинг, `SaltRemover.StripMol(dontRemoveEverything=True)`, найбільший фрагмент, canonical isomeric SMILES, InChIKey, дедуплікація з **медіаною** pChEMBL.

**Sanity check (не валідація):** A-967079 pchembl_median **7.17** (n=13); HC-030031 **5.21** (n=15) — узгоджуються з літературою.

**Provenance:** зафіксовано час витягу UTC, версію бази, дату релізу, запит. Порівняння з попереднім знімком ChEMBL 36 дало повну ідентичність підмножини (0 нових сполук, 0 змінених агрегатів). Це **version stability**, а не зовнішня валідація — release-to-release валідація для цієї мішені неможлива.

### 4.2. Представлення `[Є]`

| Представлення | Розмірність | Деталі |
|---|---|---|
| Morgan ECFP4 | 2048 | `GetMorganGenerator(radius=2, fpSize=2048)`, **achiral — стандартний baseline** |
| RDKit дескриптори | 15 | перелік нижче |
| ChemBERTa-77M-MTR | 384 | CLS pooling |
| MolFormer-XL-both-10pct | 768 | Mean pooling |

**15 дескрипторів:** MolWt, MolLogP, MolMR, TPSA, NumHAcceptors, NumHDonors, NumRotatableBonds, NumAromaticRings, RingCount, FractionCSP3, HeavyAtomCount, NumAliphaticRings, NumSaturatedRings, NumHeteroatoms, LabuteASA.

**Хіральність:** achiral ECFP4 лишається основним baseline. Chiral ECFP4 — **опційний sensitivity analysis у Supplementary**, не заміна. Міняти представлення після перегляду результатів не можна.

**Відтворюваність MolFormer:** потрібні `transformers==4.44.2` і пін ревізії `7b12d946c181a37f6012b9dc3b002275de070314`. → Supplementary.

### 4.3. Валідація

- Каркаси Bemis-Murcko: **544 унікальних** на 1645 сполук.
- **Перевірено `[Є]`:** ациклічних молекул у наборі **0 з 1645** — порожній каркас не утворює штучної супергрупи, розбиття коректні.
- `GroupKFold(n_splits=5)`, `groups=scaffold`.
- **`[ТРЕБА]`** Однакові фолди для всіх моделей + збереження out-of-fold прогнозів.
- **`[ТРЕБА]`** Повторні scaffold-розбиття (кілька partitions), бо один детермінований GroupKFold дає лише один варіант розподілу 544 каркасів по 5 фолдах.
- **`[ТРЕБА]`** Applicability domain: max Tanimoto кожної тестової молекули до train.

**Метрики в основному тексті:** RMSE, R², Spearman ρ.
**У Supplementary:** AUC і MCC за порогом 7.0 (штучна бінаризація регресійної задачі відкриває зайві питання: чому 7.0, чи збалансовані класи, чи оптимізувався поріг).

### 4.4. Моделі

```
RandomForestRegressor(n_estimators=500, random_state=42)
XGBRegressor(n_estimators=500, max_depth=6, learning_rate=0.05, random_state=42)

D-MPNN (Chemprop 2.2.3):
  BondMessagePassing(d_h=300, depth=3), MeanAggregation()
  RegressionFFN(input_dim=300 [+15 при x_d])
  warmup_epochs=2, init_lr=1e-4, max_lr=1e-3, final_lr=1e-4
```

**`[ТРЕБА]` Три виправлення D-MPNN:**
1. Масштабування `x_d` **всередині кожного train-фолду** (fit scaler лише на train, передавати `X_d_transform`).
2. Early stopping не спрацював — у **кожному** фолді досягнуто `max_epochs=100`. Протокол зламаний.
3. **3-5 сідів.** Schiebroek et al.: медіанна сідова варіативність GNN **0.41** проти **0.07** у RF. Наш результат з одним сідом непридатний до публікації.

**Тюнінг:** зараз нерівноцінний між класами моделей. Або однаковий помірний бюджет через inner CV, або чесне зазначення як обмеження.

### 4.5. Sensitivity analysis агрегації `[ТРЕБА]`

Перевірити що бенчмарк не тримається на способі агрегації міток. Три варіанти набору:
1. основний — медіана всіх вимірювань;
2. одне значення на молекулу за наперед визначеним правилом;
3. без 300 молекул з відомими парними протокольними умовами.

Питання: чи не перевертається ранжування моделей.

---

## 5. Results

### 5.1. Характеристика набору `[Є]`

Склад, розподіл pIC50, кількість вимірювань на сполуку, каркаси, роки.

**Підрозділ «Гетерогенність джерел»** (короткий, у Dataset characterization — **не** окремий центр статті):
- Топ-5 assays дають 49% вимірювань, топ-20 — 81%.
- Домінантний агоніст — cinnamaldehyde (850 сполук), не AITC (265).
- Серед 22 assays з ≥20 вимірювань медіани pIC50 від **4.96 до 8.15**. *Застереження в тексті: цей розкид змішує протокол з різними хімічними серіями і не є оцінкою ефекту протоколу.*
- 393 сполуки (24%) виміряні в ≥2 assays, медіанний внутрішньосполучний розкид **0.310 log**.
- 300 сполук (18%) мають мітку усереднену між двома протокольними умовами.

**Парний контраст — 3-4 речення, обережно:**

Два записи одного FLIPR-протоколу з кожного з двох патентів, ті самі сполуки:

| Документ | Спільних | Каркасів | Δ (медіана scaffold-медіан) | 95% CI | Узгодж. напрямку |
|---|---|---|---|---|---|
| CHEMBL5727829 | 169 | 70 | +0.200 | [+0.150, +0.267] | 92% |
| CHEMBL5727828 | 131 | 51 | +0.290 | [+0.220, +0.350] | 95% |

CI — bootstrap ресемплінгом scaffold-медіан (10 000 повторів).

**Формулювання:** *"Two ChEMBL records of the same assay protocol within a single patent document showed a systematic pIC50 offset of 0.20-0.29 for identical compounds. Per the source patent, compounds were incubated for either 10 or 90 minutes prior to agonist addition, but the mapping of individual ChEMBL assay records to these conditions could not be confirmed from the available metadata."*

**Не писати:** «виявлено помилку ChEMBL», «параметр X спричиняє зсув», «пояснює залишковий шум Landrum». Порівняння з їхнім MAE 0.27 — лише **порівняння масштабу**, і лише в Discussion, одним реченням.

**Мета підрозділу:** показати чому агрегований endpoint треба інтерпретувати обережно. Не окремий внесок статті.

### 5.2. Головна таблиця — сітка представлення × алгоритм

Заповнена сітка з п.2 + два end-to-end пайплайни. Mean±std по фолдах.

**Формулювання поки парний тест не зроблено:**
> The top pipelines showed similar mean cross-validation performance, with differences small relative to fold-to-fold variability.

**Після парного тесту** можна додати слово *statistically*.

**`[БАГ]`** Spearman для Morgan-моделей захардкоджено 0.000 при ручному складанні таблиці — перерахувати.

**`[ТРЕБА]`** Paired bootstrap на OOF-прогнозах, ресемплінг за каркасами, CI для Δmetric між парами.

### 5.3. Applicability domain `[ТРЕБА]`

Розподіл max Tanimoto тестових молекул до train; залежність похибки від схожості.

---

## 6. Discussion

**6.1.** Узгодженість з ширшою літературою: fingerprints конкурентні з навченими представленнями (Adamczyk & Czech; Sadeghi et al.), GNN не дають переваги на малих наборах (Schiebroek et al. 2026). **Ми не претендуємо на новизну цього висновку** — показуємо що він відтворюється на TRPA1.

**6.2.** Специфіка TRPA1: ковалентна активація електрофілами, кілька сайтів зв'язування антагоністів, Ca²⁺-залежна десенситизація, різкі міжвидові відмінності. IC50 може відображати різні механізми — competitive antagonism, pore block, agonist-induced desensitization. Обґрунтування терміну «inhibitory potency».

**6.3.** Обмеження (розгорнуто, п.8).

**6.4.** Напрямок далі: ручна протокольна анотація 97 assays — окрема триваліша робота.

---

## 7. Supplementary

Усе що не витримує головного тексту:

- Реплікація Mihai (RF 98.64%, SVM 98.64%, FFNN 99.09%) — з формулюванням «recreated», не «replicated», і з обмеженнями decoys.
- Fine-tuning (ChemBERTa 0.247, MolFormer 0.193) — «indicative rather than conclusive», не в abstract і не у висновках.
- Pooling (ChemBERTa CLS>Mean, MolFormer Mean>CLS) — **без причинного заголовка**, з визнанням test-set peeking.
- Descriptor associations (TPSA permutation 0.276, Spearman +0.520; контроль на рік 0.520→0.505) — не називати SAR, не будувати Discussion навколо.
- AUC / MCC.
- Chiral ECFP4 sensitivity analysis.
- Технічні нотатки MolFormer.

---

## 8. Limitations

1. Одна мішень, висновки не узагальнюються.
2. Цільова змінна — **агрегована** літературна величина, не властивість молекули.
3. Гіперпараметри тюнились нерівноцінно між класами моделей.
4. D-MPNN: early stopping не спрацював у жодному фолді.
5. Два трансформери однієї архітектурної родини; GROVER, Uni-Mol, graph transformers не тестувались.
6. Pooling обирався з огляду на test performance (test-set peeking) — результат у Supplementary.
7. IC50 не гарантує механізм антагонізму; assay-контекст вручну не анотований для 97 assays.
8. Набір змішує протокольні умови; 18% сполук мають усереднену мітку.
9. Decoys у реплікації Mihai не property-matched.
10. Publication bias.
11. Release-to-release зовнішня валідація неможлива.
12. Один тип scaffold-розбиття (пом'якшується повторними partitions).

---

## 9. Роботи перед написанням Results

| # | Задача | Оцінка |
|---|---|---|
| 1 | Перерахувати Spearman для Morgan | 0.5 год |
| 2 | Заповнити 3 клітинки сітки (XGB+RDKit15, RF+CB-CLS, RF+MF-Mean) | 1 год |
| 3 | Однакові фолди + збереження OOF для всіх моделей | 3 год |
| 4 | Повторні scaffold-partitions | 3 год |
| 5 | D-MPNN: x_d scaling + early stopping + 3-5 сідів × 2 конфігурації | 1-2 дні |
| 6 | Paired bootstrap CI для Δmetric | 4 год |
| 7 | Applicability domain | 3 год |
| 8 | Sensitivity analysis агрегації | 4 год |
| 9 | Рисунки | 1 день |

**Реалістично 5-10 робочих днів**, не два. D-MPNN сам по собі — 30-50 тренувань плюс налагодження.

Intro і Methods (4.1-4.2) можна писати паралельно — вони від цих чисел не залежать.

---

## 10. Карта цитувань

**TRPA1 ML — прямі попередники:**
- Mihai D.P. et al. *AI* 2020, 1(2):276-285. DOI 10.3390/ai1020018
- **Gawalska A. et al. *Molecular Informatics* 2023. DOI 10.1002/minf.202200214** — головний прямий попередник для регресії
- Castañeda-Leautaud & Amaro 2025, arXiv:2509.10871

**Споріднені TRP:**
- Wei X. et al. *Molecules* 2024, 29(2):295
- Lee et al. *Applied Sciences* 2023, 13(9):5617

**Якість даних ChEMBL:**
- Landrum G.A., Riniker S. *JCIM* 2024, 64:1560-1567
- Schiebroek C.C.G., Landrum G.A., Riniker S. *JCIM* 2026, 66:7446-7452
- Kalliokoski T. et al. *PLoS One* 2013, 8:e61007
- Kramer C. et al. *J Med Chem* 2012, 55:5165-5173

**Представлення / бенчмарки:**
- Adamczyk & Czech, arXiv:2508.06199
- Sadeghi S. et al. *BMC Bioinformatics* 2024
- Hasan et al. Rep3Net, arXiv:2512.00521

**Моделі:**
- ChemBERTa-2: Ahmad et al. arXiv:2209.01712
- MolFormer: Ross et al. *Nat Mach Intell* 2022
- Chemprop: Heid E. et al. *JCIM* 2024, 64:9-17
- RDKit, scikit-learn, XGBoost

**TRPA1 фармакологія:**
- Структурні роботи по сайтах зв'язування (cryo-EM)
- Міжвидові відмінності TRPA1
- Gawalska A. et al. *Molecules* 2022, 27:3077

---

## 11. Відкриті дані

- `trpa1_ic50_chembl37.csv` — InChIKey, SMILES, pchembl_median, n_measurements, scaffold, assay/document IDs
- `assay_audit.csv` — 97 assays з метаданими і автоматичними прапорцями
- Скрипти: витяг, стандартизація, бенчмарк, протокольна декомпозиція
- `metadata.json` — версія ChEMBL, дата витягу, точний запит

Відкритий набір з provenance — реальна частина внеску, можливо цінніша за самі метрики.
