# Технічний звіт: ML-бенчмарк для прогнозування потужності антагоністів TRPA1

**Призначення документа:** повний технічний опис виконаного дослідження для незалежного рецензування. Містить усі параметри, результати, фрагменти коду, виявлені проблеми та методологічні застереження. Проза — українською, код і технічні терміни — англійською.

**Прохання до рецензента:** оцінити (1) методологічну коректність, (2) чи достатньо новизни для публікації, (3) чи є невиявлені помилки, (4) що варто зробити далі. Не пом'якшувати критику.

---

## 0. Контекст і мета

- **Автор:** аспірант 2-го року, біофізика іонних каналів (Інститут фізіології ім. Богомольця, Київ). Працює віддалено. Windows 11, Ryzen 7 260, 16 GB RAM, RTX 5060 Laptop (Blackwell, sm_120, 8 GB). Оточення: `mamba env "trpa1"`, Python 3.13.
- **Мета:** підготувати публікабельну статтю (український фаховий журнал або рівень MDPI) для допуску до захисту. НЕ дисертація — саме стаття.
- **Відправна точка:** Mihai et al. (2020), *AI* (MDPI), DOI 10.3390/ai1020018 — приклад роботи мінімальними ресурсами.
- **Наукове питання:** чи дають сучасні архітектури (GNN, pretrained chemical language models) перевагу над класичними методами для прогнозування потужності антагоністів TRPA1 — мішені, для якої такого порівняння ніколи не робили.
- **Методологічний принцип роботи:** строга покрокова діагностика — писати мінімальний скрипт, запускати, дивитись реальний вивід, потім рухатись далі. Ніколи не екстраполювати наперед. Усі результати перехресно перевірялись через Gemini і другий екземпляр Claude як незалежних рецензентів.

---

## 1. Формування набору даних

### 1.1. Джерело і витягування

- База: **ChEMBL 36** (реліз жовтень 2025), через `chembl_webresource_client` (REST API, не SQLite dump — обрано заради простоти).
- Верифіковані ChEMBL target IDs (кілька спочатку були галюцинованими, потім виправлені):
  - TRPA1 human = **CHEMBL6007**
  - TRPA1 rat = CHEMBL5160
  - TRPV1 = CHEMBL4794, TRPM8 = CHEMBL1075319, TRPV3 = CHEMBL5522, TRPV4 = CHEMBL3119, hERG = CHEMBL240
- Сирий pull: **2481 activities, 1822 унікальних сполук**.

**Важливі нюанси API, виявлені емпірично:**
- `pchembl_value` і `standard_value` повертаються як **рядки** (str), не числа — потребують конверсії.
- `action_type` завжди `null`.
- `assay_type` ненадійний ('B' навіть для activation assays) — не використовувати для розділення агоніст/антагоніст.
- Розділення антагоніст/агоніст робили за `standard_type`: IC50 → антагоністи, EC50 → агоністи.

### 1.2. Стандартизація

Пайплайн через RDKit:
1. `Chem.MolFromSmiles` (парсинг);
2. `SaltRemover.StripMol(mol, dontRemoveEverything=True)` (видалення солей);
3. якщо лишається кілька фрагментів — беремо найбільший за `GetNumHeavyAtoms`;
4. `Chem.MolToSmiles(mol, canonical=True)` (canonical SMILES);
5. `Chem.MolToInchiKey` (InChIKey для дедуплікації);
6. дедуплікація за `(InChIKey, standard_type)`, агрегація через **median pchembl**.

Результат: **1827 рядків, 1814 унікальних сполук**. Розподіл: антагоністи (IC50) = **1645**, агоністи (EC50) = **182**.

```python
from rdkit import Chem
from rdkit.Chem import SaltRemover

remover = SaltRemover.SaltRemover()

def standardize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mol = remover.StripMol(mol, dontRemoveEverything=True)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    can = Chem.MolToSmiles(mol, canonical=True)
    ik  = Chem.MolToInchiKey(mol)
    return can, ik
```

### 1.3. Валідація якості даних

Перевірка за еталонними антагоністами (літературні значення збіглися точно):
- **A-967079** (CHEMBL3697701): pchembl_median = **7.17**, n=13 ✓
- **HC-030031** (CHEMBL1086310): pchembl_median = **5.21**, n=15 ✓ (SMILES: `CC(C)c1ccc(NC(=O)Cn2cnc3c2c(=O)n(C)c(=O)n3C)cc1`)
- JT010 (CHEMBL5219660): відсутній у ChEMBL для TRPA1 (дані Takaya 2015 не депоновані — нормально).

Розподіл pchembl_median антагоністів: min 4.00, max 9.31, median 7.32.

### 1.4. Decoys (для відтворення Mihai)

Проблема: витягування випадкових молекул з ChEMBL. Три підходи провалились через зависання API:
- Guessing випадкових ID по одному — тисячі неіснуючих ID, кожен на таймауті.
- `itertools.islice` з великим offset — клієнт гортає всі записи до offset сторінками по 1000 (1100+ запитів).

**Робоче рішення:** forward-only streaming з фільтром MW, потім random subsample.

```python
qs = molecule.filter(
    molecule_properties__mw_freebase__gte=150,
    molecule_properties__mw_freebase__lte=900,
    molecule_type='Small molecule',
).only(['molecule_chembl_id', 'molecule_structures'])
# stream 5000, then random.shuffle + take 560
```

Результат: **560 decoys**, стандартизовані, 0 колізій з TRPA1-набором (перевірено за InChIKey), 0 внутрішніх дублів. MW 152–898 (mean 382; Mihai мали 113–1089, mean 418). Співвідношення антагоністи:decoys = **2.94:1** (Mihai мали 2.92:1).

### 1.5. Фінальні файли

`trpa1_antagonists.csv` (1645, колонки: inchikey, standard_type, std_smiles, molecule_chembl_id, pchembl_median/min/max/std, n_measurements, assay_types, year_min/max, **scaffold**, **split**), `trpa1_agonists.csv` (182), `trpa1_human_clean.csv` (1827), `decoys_clean.csv` (560).

---

## 2. Молекулярні представлення

| Representation | Розмірність | Деталі |
|---|---|---|
| Morgan ECFP4 | 2048 біт | `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)` |
| RDKit descriptors | 15 | фізико-хімічні (список нижче) |
| ChemBERTa-77M-MTR | 384 | CLS та Mean Pooling |
| MolFormer-XL-both-10pct | 768 | CLS та Mean Pooling |

**15 RDKit дескрипторів:** MolWt, MolLogP, MolMR, TPSA, NumHAcceptors, NumHDonors, NumRotatableBonds, NumAromaticRings, RingCount, FractionCSP3, HeavyAtomCount, NumAliphaticRings, NumSaturatedRings, NumHeteroatoms, LabuteASA.

**Витягування ембедингів** (обидва pooling для порівняння):

```python
def extract(smiles, model, tokenizer, pooling, max_length):
    tokens = tokenizer(smiles, return_tensors="pt", truncation=True,
                       padding=True, max_length=max_length).to('cuda')
    with torch.no_grad():
        output = model(**tokens)
    hidden = output.last_hidden_state           # (1, seq_len, dim)
    if pooling == 'cls':
        return hidden[:, 0, :].cpu().numpy().ravel()
    else:  # mean pooling over non-pad tokens
        mask = tokens['attention_mask'].unsqueeze(-1).float().to('cuda')
        return ((hidden * mask).sum(1) / mask.sum(1)).cpu().numpy().ravel()
```

max_length: ChemBERTa 128, MolFormer 202.

---

## 3. Експериментальний дизайн

- **Scaffold split** (Bemis-Murcko generic scaffold): train 1324 / val 160 / test 161. Детерміністичний.
- **544 унікальних scaffolds** на 1645 сполук.
- **Перехресна валідація:** `GroupKFold(n_splits=5)` з `groups=scaffold` — жоден scaffold не потрапляє одночасно в train і test.
- **Метрики:** RMSE, R² (регресія pIC50); Spearman ρ (ранжування); AUC-ROC та MCC для класифікації за порогом pchembl ≥ **7.0** (при цьому порозі ~64% активних).
- **Довірчі інтервали:** bootstrap (1000 resamples) для single-split; mean±std по 5 фолдах для CV.

---

## 4. Відтворення Mihai et al. (2020)

### 4.1. Постановка

Задача Mihai — **бінарна класифікація ligand-vs-decoy** (не потужний-vs-слабкий). Усі 371 антагоніст = клас 1, 127 випадкових = клас 0. Пряма цитата: *"predict TRPA1 ligands with antagonistic activity within any range of potencies, rather than only antagonists with high potency."*

Ми відтворили постановку на ChEMBL 36 (1645 антагоністів + 560 decoys), random 80/20 stratified split (seed=34), з **точними гіперпараметрами** авторів, але Morgan ECFP4 замість MNA level-3 (MNA — пропрієтарні дескриптори PASS, неможливо точно відтворити). Логіка: якщо результат збігається — високий показник йде від легкості задачі, не від конкретних дескрипторів.

### 4.2. Гіперпараметри (з оригіналу)

```python
RandomForestClassifier(n_estimators=50, max_depth=90, max_features='sqrt',
                       min_samples_split=2, min_samples_leaf=1, random_state=34)
SVC(C=8, gamma=0.001, kernel='rbf', probability=True)
# FFNN: 2048 -> 750 (ReLU, dropout=0.6) -> 1 (sigmoid), BCE, 30 epochs, batch=16
```

### 4.3. Результати

| Model | Mihai ACC | Наш ACC | Mihai AUC (CV) | Наш AUC (CV) |
|---|---|---|---|---|
| RF | 99.0% | **98.64%** | 0.9936 | **0.9982** |
| SVM | 90.0% | **98.64%** | 0.9354 | **0.9968** |
| FFNN | 88.0% | **99.09%** | 0.9354 | **0.9983** |

**Висновки:** (1) результат відтворюється; (2) FFNN тепер виграє (у Mihai програвав) — підтверджує їхню гіпотезу що слабкість FFNN була через малий датасет (у нас 2205 vs їхніх 498 зразків); (3) усі три моделі впираються в стелю ~99%, бо задача ligand-vs-random-decoy штучно легка. Це обґрунтувало перехід до реалістичної постановки (регресія потужності, scaffold split).

---

## 5. Основний бенчмарк (регресія потужності)

### 5.1. Гіперпараметри

```python
# Classical
RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=42)
XGBRegressor(n_estimators=500, max_depth=6, learning_rate=0.05, random_state=42)

# D-MPNN (Chemprop 2.2.3)
mp = BondMessagePassing(d_h=300, depth=3)
agg = MeanAggregation()
ffn = RegressionFFN(input_dim=300 + 15)  # +15 when adding RDKit x_d features
model = MPNN(message_passing=mp, agg=agg, predictor=ffn,
             warmup_epochs=2, init_lr=1e-4, max_lr=1e-3, final_lr=1e-4)
# EarlyStopping(monitor='val_loss', patience=10), max_epochs=100
```

Chemprop v2.2.3 API-нюанс: `MoleculeDatapoint(mol=Chem_obj, y=np.array([...]), x_d=features)` — приймає RDKit Mol об'єкт (не SMILES) і `y` як np.ndarray. При додаванні `x_d` треба `RegressionFFN(input_dim=300+n_features)`, інакше shape mismatch (64x315 vs 300x300).

### 5.2. Результати CV (5-fold GroupKFold, scaffold) — головна таблиця

| Model | R² | RMSE | Spearman | AUC | Тип |
|---|---|---|---|---|---|
| RF (Morgan) | **0.547±0.112** | 0.618±0.050 | —¹ | 0.856±0.060 | Fingerprint |
| XGB CB-CLS | 0.539±0.084 | 0.633±0.057 | 0.725±0.064 | **0.866±0.043** | Transformer frozen |
| XGB (Morgan) | 0.539±0.106 | 0.625±0.046 | —¹ | 0.860±0.052 | Fingerprint |
| XGB MF-Mean | 0.529±0.085 | 0.639±0.042 | 0.703±0.066 | 0.857±0.045 | Transformer frozen |
| D-MPNN+RDKit | 0.491±0.092 | 0.663±0.040 | 0.677±0.086 | 0.843±0.049 | GNN+Desc |
| RF (RDKit-15) | 0.478±0.163 | 0.664±0.053 | 0.695±0.089 | 0.848±0.051 | Descriptors |
| D-MPNN (graph) | 0.410±0.107 | 0.714±0.031 | 0.597±0.089 | 0.794±0.050 | GNN |

¹ **BUG:** Spearman для Morgan-моделей у фінальній таблиці захардкоджено 0.000 — Morgan CV рахувався локально, склейка таблиці була ручна, значення втрачено. Потребує перерахунку.

Per-fold R² для D-MPNN graph: [0.254, 0.532, 0.334, 0.405, 0.522]. Для D-MPNN+RDKit: [0.360, 0.582, 0.410, 0.517, 0.588].

### 5.3. Пулінг: CLS vs Mean залежить від pretraining

Ключове спостереження. На single split:
- **ChemBERTa: CLS (0.344) > Mean (0.298).** Причина: модель MTR (Multi-Task Regression) тренована з regression head на CLS-токені.
- **MolFormer: Mean (0.418) > CLS (0.262).** Причина: тренована на masked language modeling без спеціального CLS-завдання, інформація розмазана по токенах.

Gemini рекомендував Mean Pooling для обох — був правий для MolFormer, помилявся для ChemBERTa. Урок: перевіряти емпірично, а не за загальним правилом для RoBERTa.

### 5.4. Ablation (критичний експеримент)

Питання: чи GNN додає щось до RDKit, чи RDKit тягне все? Запустили MLP та RF на самих 15 RDKit дескрипторах (single split):
- RF (15 RDKit): R²=0.363 — **краще** за D-MPNN+RDKit (0.292) на single split.
- MLP (15 RDKit): R²=-0.116 (провал — MLP погано працює на 15 фічах без масштабування/тюнінгу).

Висновок: GNN-компонент не додавав корисного сигналу понад дескриптори; на single split навіть погіршував. Це врятувало від хибного твердження «гібрид оптимальний». (На CV картина вирівнялась — див. 5.2.)

---

## 6. Fine-tuning мовних моделей (frozen vs end-to-end)

### 6.1. ChemBERTa fine-tuning

```python
model = AutoModelForSequenceClassification.from_pretrained(
    "DeepChem/ChemBERTa-77M-MTR", num_labels=1,
    problem_type="regression", ignore_mismatched_sizes=True)
# Замінює 199-task regression head на single output.
TrainingArguments(num_train_epochs=50, per_device_train_batch_size=32,
    learning_rate=2e-5, warmup_ratio=0.1, weight_decay=0.01,
    load_best_model_at_end=True, metric_for_best_model='eval_loss', fp16=True)
# EarlyStoppingCallback(early_stopping_patience=5)
```

Результат: **R²=0.247, Spearman=0.497, AUC=0.739** — гірше за frozen embeddings (0.539).

### 6.2. MolFormer fine-tuning

MolFormer з `fp16=True` давав **NaN** на першій evaluation (linear attention + random feature maps нестабільні в half precision). Рішення від зовнішнього рецензента (перевірено, коректне): `deterministic_eval=True` + backbone завжди в eval mode + двоетапне тренування (Stage 1: head only 3 epochs lr=1e-3; Stage 2: end-to-end backbone lr=5e-6, head lr=1e-4, CosineAnnealingLR).

Результат: **R²=0.193, Spearman=0.567, AUC=0.759** — теж гірше за frozen (0.529).

### 6.3. Висновок (стійкий на двох моделях)

На датасеті ~1324 train молекул fine-tuning трансформерів **систематично програє** frozen extraction. Причина — catastrophic forgetting: ~3.4M (ChemBERTa) / ~44M (MolFormer) параметрів проти малого датасету. Frozen embeddings зберігають знання з 77M/1.1B молекул незайманим, XGBoost вчиться поверх стабільного представлення. Це самостійний публікабельний висновок для low-data regime.

Деталь: у MolFormer val_r2 доходив до 0.54, val_spearman до 0.74, але test R²=0.19 — класичний scaffold overfitting (добре на validation scaffolds, погано на нових test scaffolds).

---

## 7. Feature importance та SAR

RF на 15 RDKit дескрипторах, 5-fold CV. Impurity-based + permutation importance + signed Spearman (напрямок ефекту).

| Дескриптор | Perm. Imp | Impurity | Spearman | Ефект |
|---|---|---|---|---|
| TPSA (полярна поверхня) | **0.276±0.045** | 0.333 | +0.520 | ↑ потужність |
| Heteroatoms | 0.073±0.041 | 0.095 | +0.507 | ↑ |
| FractionCSP3 | 0.058±0.005 | 0.074 | +0.098 | ↑ |
| MolMR | 0.055±0.017 | 0.067 | +0.457 | ↑ |
| HDonors | 0.047±0.020 | 0.028 | **-0.101** | ↓ |
| LogP | 0.046±0.004 | 0.078 | **-0.252** | ↓ |
| Aromatic rings | 0.046±0.024 | 0.014 | +0.501 | ↑ |
| HAcceptors | 0.042±0.023 | 0.077 | +0.541 | ↑ |

**SAR-наратив:** потужність антагоністів TRPA1 визначається насамперед **полярною топологією** — велика TPSA, багато акцепторів, багато гетероатомів підвищують; надмірна ліпофільність і надлишок донорів знижують. Профіль потужного антагоніста: полярні акцепторні групи (карбоніли, N у гетероциклах, сульфони), помірна ліпофільність, мало донорів. Узгоджується з еталонами (A-967079: оксим + фторфеніл; HC-030031: пуриновий каркас з акцепторами).

TPSA домінує в ~4 рази над наступною ознакою; impurity і permutation узгоджуються.

---

## 8. Confound analysis (рік публікації)

Ризик: чи «біологічний» TPSA-сигнал не є артефактом того, що новіші сполуки полярніші за дизайном?

- pIC50 vs year_min: Spearman **0.149** (слабкий тренд).
- TPSA vs year: 0.217; Heteroatoms vs year: 0.206 (мілкі тренди).
- **Partial correlation** (TPSA vs pIC50, controlling for year): 0.520 → **0.505** (падіння 0.015, «signal holds»).

Для всіх топ-дескрипторів partial ≈ raw. **Висновок: SAR-сигнал не пояснюється роком публікації** — реальна структурна закономірність. Це закриває найочевиднішу атаку рецензента.

```python
def partial_corr(x, y, z):  # Spearman partial, x,y controlling for z
    from scipy.stats import rankdata
    xr, yr, zr = rankdata(x), rankdata(y), rankdata(z)
    def resid(a, b):
        A = np.vstack([b, np.ones_like(b)]).T
        coef, *_ = np.linalg.lstsq(A, a, rcond=None)
        return a - A @ coef
    return np.corrcoef(resid(xr, zr), resid(yr, zr))[0, 1]
```

---

## 9. Технічні проблеми та рішення (для відтворюваності)

### 9.1. Blackwell GPU (sm_120)
RTX 5060 Laptop — Blackwell sm_120. PyTorch stable (вкл. 2.9.0) не підтримує. Nightly cu129 нестабільні. **Рішення:** Chemprop на CPU (3.8 s/epoch для 1324 молекул — цілком прийнятно, GPU не потрібен). Локальний torch згодом зламався повністю (`shm.dll` procedure not found — bug Windows-білдів) — трансформери й D-MPNN CV перенесли в Colab.

### 9.2. MolFormer version hell (найдовша боротьба)
Каскад несумісностей кастомного коду MolFormer з різними версіями transformers:
- `ModuleNotFoundError: transformers.onnx` (прибрано в transformers v5)
- `cannot import find_pruneable_heads_and_indices from transformers.pytorch_utils` (видалено в 4.57)
- `cannot import create_bidirectional_mask from transformers.masking_utils` (нова версія коду MolFormer вимагає новішого transformers)
- `MolformerModel object has no attribute get_head_mask`

Корінь: HF автоматично тягне **найновішу** ревізію кастомного коду (hash `a14249e5`), яка конфліктує зі старим transformers; але старий код (hash `7b12d946`) вимагає старих API яких немає в новому transformers.

**Робоче рішення:** `transformers==4.44.2` + **пін старої ревізії коду**:

```python
mf_model = AutoModel.from_pretrained(
    "ibm/MoLFormer-XL-both-10pct",
    deterministic_eval=True, trust_remote_code=True,
    revision="7b12d946c181a37f6012b9dc3b002275de070314")
```

### 9.3. Colab session pitfalls
- rdkit-pypi більше не існує на Python 3.12 → `pip install rdkit`.
- Runtime restart скидає всі pip-пакети і GPU-runtime (треба перевибирати T4).
- Одна нова сесія стартувала на CPU → D-MPNN CV ішов годинами замість хвилин.

---

## 10. Відомі обмеження та методологічні застереження

**Критичні (потребують виправлення перед submission):**
1. **Spearman bug** для Morgan-моделей (0.000 у таблиці) — технічна помилка склейки, перерахувати.
2. **Гіперпараметри не тюнились однаково.** Класика (розумні дефолти) і D-MPNN (стандарт) без GridSearch/Optuna. Рецензент справедливо спитає, чи D-MPNN програв не через відсутність тюнінгу.
3. **Early stopping не спрацював у D-MPNN CV** — `max_epochs=100 reached` у кожному фолді. D-MPNN отримав повні 100 епох. Працює на користь висновку (навіть з форою програв), але training protocol технічно зламаний — треба полагодити моніторинг val_loss.

**Суттєві:**
4. Лише **два трансформери, обидва RoBERTa-подібні** (ChemBERTa, MolFormer) — одна архітектурна родина. Не тестовано GROVER (GNN-transformer), Uni-Mol (3D), graph transformers.
5. **Decoys не property-matched** (та сама проблема що в Mihai) — модель могла вчитись тривіальним властивостям (MW, LogP). Для регресії потужності це менш критично (там decoys не використовуються), але для реплікації — так.
6. **Publication bias:** усі 1645 антагоністів — «справжні» (немає true inactives для регресії).
7. Малий тест (161) на single split — але це закрито через CV.

**Не помилки (для балансу):** курація чиста, scaffold split і GroupKFold коректні, confound check зроблено (рідкість), reference compounds валідовані, реплікація коректна, bootstrap/mean±std — правильна статистика, frozen-vs-fine-tuning порівняння чесне.

---

## 11. Оцінка новизни (deep search по літературі)

**Захищені твердження:**
- Перший systematic benchmark pretrained transformer embeddings + GNN vs класичних методів для TRPA1 antagonist pIC50 regression.
- Перше застосування будь-якої GNN до TRPA1.
- Перше застосування pretrained chemical language models до будь-якого TRP-каналу.
- Найбільший курований набір антагоністів TRPA1 (1645, ChEMBL 36).

**Найближча прилегла робота:**
- **Mihai et al. 2020** — TRPA1, класика, класифікація (наша відправна точка).
- **Wei et al. 2024** (*Molecules*) — TRPV1, лише класичні методи (SVM/Bagging/GBDT/XGBoost), R²~0.78-0.81, інші дані — не порівнянно.
- **Lee et al. 2023** — MT-DTI transformer на TRPV1 для repurposing (трансформер-на-TRP існує, але не TRPA1, не benchmark).
- **Rep3Net** (arXiv:2512.00521, 2025) — **той самий рецепт** (RDKit + ChemBERTa + graph) для PARP1, R²=0.43. Наш pipeline методологічно не новий; новизна — target (TRPA1) + масштаб порівняння.

**Уникати тверджень:** «перший ML на TRPA1» (Mihai зробив), «перший трансформер на TRP» (Lee зробив TRPV1), «representation > algorithm» як загальний висновок (добре відомо: Jiang 2021, Adamczyk & Czech arXiv:2508.06199, Sadeghi 2024).

**Чесна оцінка сили:** головний результат («складність не окуповується на малих даних») — **відомий факт**, вперше показаний на TRPA1. Це інкрементальна новизна. Для українського журналу — достатньо. Для MDPI — на межі, потрібне посилення. Framing вирішує: не «трансформери революціонізують TRPA1», а «систематичний строгий benchmark representation learning для розуміння TRPA1 SAR».

---

## 12. Заплановані наступні кроки

**Виправлення методології:**
- Перерахувати Spearman для всіх моделей на єдиній основі.
- Однакова оптимізація гіперпараметрів (Optuna) для всіх підходів.
- Полагодити early stopping у D-MPNN CV.
- Property comparison actives vs decoys (розподіли MW/LogP/TPSA) + caveat про DUD-E в Discussion.

**Прикладний результат (скринінг):**
- Застосувати найкращу модель до бібліотеки комерційно доступних сполук (Enamine REAL / ZINC) для прогнозу нових антагоністів TRPA1.
- Відібрати 20–50 кандидатів з обґрунтуванням структурної привабливості.
- Це прямо закриває формулювання теми «прогнозування нових лігандів».

**(Опційно, для сильнішого журналу):** додати принципово інший тип моделі (Uni-Mol 3D або GROVER); multi-task на TRP-родину (TRPV1/TRPM8); експериментальна валідація кількох кандидатів у лабораторії (calcium imaging / patch-clamp).

---

## Додаток: зведення всіх числових результатів

**Реплікація Mihai (ligand-vs-decoy, random split):** RF 98.64% / SVM 98.64% / FFNN 99.09% ACC; AUC CV 0.998/0.997/0.998.

**Single-split scaffold (test=161), R²:** RF-Morgan 0.232, XGB-Morgan 0.257, LightGBM-Morgan 0.255, D-MPNN 0.176, D-MPNN+RDKit(ES) 0.292, RF-RDKit15 0.363, XGB-CB-CLS 0.344, XGB-CB-CLS+RDKit 0.395, XGB-MF-Mean 0.418(single), ChemBERTa-FT 0.247, MolFormer-FT 0.193.

**CV (5-fold GroupKFold), R²:** RF-Morgan 0.547±0.112, XGB-CB-CLS 0.539±0.084, XGB-Morgan 0.539±0.106, XGB-MF-Mean 0.529±0.085, D-MPNN+RDKit 0.491±0.092, RF-RDKit15 0.478±0.163, D-MPNN 0.410±0.107.

**Розбіжність single-split vs CV** для трансформерів (напр. CB-CLS 0.344 → 0.539) пояснюється тим, що конкретний single test=161 був одним із важчих; CV усереднює по 5 нарізках і дає надійнішу оцінку.
