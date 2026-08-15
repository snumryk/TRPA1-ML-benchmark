"""
protocol_decomposition.py

ЩО ЦЕЙ СКРИПТ РОБИТЬ (простими словами)
---------------------------------------
Кожен з 97 експериментів має текстовий опис. Скрипт розбирає цей текст на
окремі параметри протоколу: які клітини, який активатор, яка концентрація,
скільки часу чекали, яким приладом міряли.

ВАЖЛИВО — чого скрипт НЕ робить:
  * не оголошує жодну групу "однорідною" (це вирішує людина);
  * не робить висновків про причину різниць;
  * не вигадує числа коли текст неоднозначний — ставить "ambiguous".

Для кожного витягнутого параметра зберігається:
  * саме значення;
  * фрагмент тексту з якого воно взяте (evidence);
  * рівень впевненості.

ВИХІДНІ ФАЙЛИ
-------------
  protocol_decomposition.csv  — по одному рядку на експеримент, усі параметри
  protocol_families.csv       — експерименти згруповані за однаковим протоколом
  paired_contrasts.csv        — пари експериментів що відрізняються одним параметром
  protocol_metadata.json      — provenance

ВХІДНІ ФАЙЛИ (мають лежати поруч)
---------------------------------
  trpa1_assay_audit.csv       — описи 97 експериментів
  trpa1_current_api_raw.csv   — сирі виміри з assay_chembl_id
"""

from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    HAVE_RDKIT = True
except ImportError:
    HAVE_RDKIT = False

try:
    from scipy import stats as sps
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


AUDIT_FILE = Path("trpa1_assay_audit.csv")
RAW_FILE = Path("trpa1_current_api_raw.csv")

OUT_DECOMP = Path("protocol_decomposition.csv")
OUT_FAMILIES = Path("protocol_families.csv")
OUT_PAIRED = Path("paired_contrasts.csv")
OUT_META = Path("protocol_metadata.json")

MIN_SHARED_COMPOUNDS = 20
EXPECTED_TARGET = "CHEMBL6007"   # human TRPA1   # мінімум спільних молекул щоб рахувати пару


# =====================================================================
# ЧАСТИНА 1. Розбір опису на параметри протоколу
# =====================================================================

def _find(text: str, pattern: str, flags=re.I) -> tuple[str | None, str | None]:
    """Шукає патерн. Повертає (значення_групи_1, фрагмент_тексту_доказ)."""
    if not isinstance(text, str):
        return None, None
    m = re.search(pattern, text, flags)
    if not m:
        return None, None
    value = m.group(1) if m.groups() else m.group(0)
    start = max(0, m.start() - 60)
    end = min(len(text), m.end() + 60)
    evidence = text[start:end].replace("\n", " ").strip()
    return value, evidence


def _find_all(text: str, pattern: str, flags=re.I) -> list[str]:
    if not isinstance(text, str):
        return []
    return [m.group(1) if m.groups() else m.group(0)
            for m in re.finditer(pattern, text, flags)]


def parse_protocol(description: str, title_hint: str = "") -> dict:
    """
    Розбирає опис експерименту на параметри.
    Кожен параметр отримує три поля: значення, доказ, впевненість.
    """
    d = description if isinstance(description, str) else ""
    out: dict = {}

    def put(field, value, evidence, confidence):
        out[field] = value
        out[f"{field}__evidence"] = evidence
        out[f"{field}__confidence"] = confidence

    # ---------- метод вимірювання ----------
    modality, ev = None, None
    if re.search(r"\bFLIPR\b|calcium indicator dye|calcium influx|Ca2\+ influx|"
                 r"calcium flux|fluo-?4|fura-?2|calcium mobilization", d, re.I):
        modality, (_, ev) = "calcium_fluorescence", _find(
            d, r"(FLIPR|calcium indicator dye|calcium influx|fluo-?4|fura-?2)")
    if re.search(r"patch[- ]clamp|whole[- ]cell|voltage[- ]clamp|electrophysiolog", d, re.I):
        m2, ev2 = _find(d, r"(patch[- ]clamp|whole[- ]cell|electrophysiolog\w*)")
        if modality is None:
            modality, ev = "electrophysiology", ev2
        else:
            modality, ev = "mixed_calcium_and_ephys", f"{ev} || {ev2}"
    put("assay_modality", modality, ev,
        "explicit_chembl_description" if modality else "unknown")

    # ---------- активатор (агоніст) ----------
    agonists = []
    if re.search(r"\bAITC\b|\ballyl[\s-]*isothiocyanate\b|\bmustard[\s-]*oil\b", d, re.I):
        agonists.append("AITC")
    if re.search(r"cinnamald", d, re.I):
        agonists.append("cinnamaldehyde")
    if re.search(r"\bJT010\b", d, re.I):
        agonists.append("JT010")
    if re.search(r"carvacrol", d, re.I):
        agonists.append("carvacrol")
    if re.search(r"\bNMM\b|N-methylmaleimide", d, re.I):
        agonists.append("NMM")
    if re.search(r"\bH2O2\b|hydrogen peroxide", d, re.I):
        agonists.append("H2O2")
    _, ag_ev = _find(d, r"(AITC|allyl[\s-]*isothiocyanate|cinnamald\w*|JT010|carvacrol)")
    put("challenge_agonist", ";".join(agonists) if agonists else None, ag_ev,
        "explicit_chembl_description" if len(agonists) == 1
        else ("ambiguous_multiple" if len(agonists) > 1 else "unknown"))

    # ---------- рівень агоніста (EC50 / EC80) ----------
    lv = _find_all(d, r"\b(EC\s?(?:50|80|90|20))\b")
    lv = sorted({x.replace(" ", "").upper() for x in lv})
    _, lv_ev = _find(d, r"((?:about\s+an?\s+)?EC\s?(?:50|80|90|20)\s+concentration[^.]{0,60})")
    put("agonist_level", ";".join(lv) if lv else None, lv_ev,
        "explicit_chembl_description" if len(lv) == 1
        else ("ambiguous_multiple" if len(lv) > 1 else "unknown"))

    # ---------- концентрація агоніста ----------
    # Число БЕЗ одиниць не вважаємо концентрацією: "cinnamaldehyde (75)" може бути
    # і 75 uM, і номером сполуки в патенті.
    conc, conc_ev = _find(d, r"(?:cinnamald\w*|AITC|agonist)[^.]{0,80}?"
                             r"(\d+(?:\.\d+)?)\s*(?:uM|µM|μM|nM|mM|micromolar)")
    if conc:
        put("agonist_concentration", conc, conc_ev, "explicit_with_units")
        out["agonist_concentration_unresolved_token"] = None
    else:
        tok, tok_ev = _find(d, r"(?:cinnamald\w*|AITC)\s*\((\d+(?:\.\d+)?)\)")
        put("agonist_concentration", None, tok_ev, "unknown")
        out["agonist_concentration_unresolved_token"] = tok

    # ---------- клітинна лінія ----------
    cells = []
    for pat, name in [(r"\bCHO\b", "CHO"), (r"\bHEK[- ]?293\w*", "HEK293"),
                      (r"\bF-?11\b", "F11"), (r"\bU-?2\s?OS\b", "U2OS"),
                      (r"\bDRG\b|dorsal root ganglion", "DRG_native")]:
        if re.search(pat, d, re.I):
            cells.append(name)
    _, cell_ev = _find(d, r"(CHO|HEK[- ]?293\w*|F-?11|DRG)[^.]{0,50}")
    put("cell_line", ";".join(cells) if cells else None, cell_ev,
        "explicit_chembl_description" if len(cells) == 1
        else ("ambiguous_multiple" if len(cells) > 1 else "unknown"))

    # ---------- вид / конструкт ----------
    species = []
    if re.search(r"\bhuman\b|\bhTRPA1\b", d, re.I):
        species.append("human")
    if re.search(r"\brat\b|\brTRPA1\b", d, re.I):
        species.append("rat")
    if re.search(r"\bmouse\b|\bmTRPA1\b", d, re.I):
        species.append("mouse")
    _, sp_ev = _find(d, r"((?:human|rat|mouse)\s+TRPA1)")
    # УВАГА: дані витягнуто за target_chembl_id=CHEMBL6007 (людський TRPA1),
    # тому вид ЗАПИСУ завжди human. Опис може згадувати й інші види, бо в патенті
    # один текст описує паралельні тести. Це лише попередження, НЕ ознака протоколу.
    put("species_mentioned_in_text", ";".join(species) if species else None, sp_ev,
        "text_mentions_only__not_record_species")
    out["record_species"] = "human"          # за конструкцією вибірки
    out["text_mentions_non_human"] = bool(set(species) - {"human"})

    is_mut = bool(re.search(r"\bmutant\b|\bmutation\b|\bC621\b|\bN855\b|chimera|"
                            r"\b[A-Z]\d{2,4}[A-Z]\b", d))
    is_wt = bool(re.search(r"\bwild[\s-]?type\b|\bWT\b", d))
    if is_mut:
        construct, conf = "mutant_or_chimera", "explicit_chembl_description"
    elif is_wt:
        construct, conf = "wild_type", "explicit_chembl_description"
    else:
        construct, conf = "not_reported", "unknown"   # відсутність згадки != WT
    _, mut_ev = _find(d, r"(mutant|mutation|C621|N855|chimera|wild[\s-]?type)")
    put("construct", construct, mut_ev, conf)

    # ---------- ТРИ РІЗНІ ЧАСИ (критично не плутати!) ----------
    # 1) завантаження барвника
    dye, dye_ev = _find(d, r"(?:loaded with|dye load\w*)[^.]{0,80}?"
                           r"for\s+(\d+(?:\.\d+)?\s*(?:hr|hour|h|min|minute)s?)")
    put("dye_loading_time", dye, dye_ev,
        "explicit_chembl_description" if dye else "unknown")

    # 2) стабілізація планшета при кімнатній температурі
    equil, eq_ev = _find(d, r"followed by\s+(\d+)\s*minutes?\s+at room temperature")
    put("plate_equilibration_min", equil, eq_ev,
        "explicit_chembl_description" if equil else "unknown")

    # 3) інкубація зі сполукою ДО додавання агоніста — головний параметр
    # ВАЖЛИВО: текст може мати форму "10 minutes or 90 minutes" — треба зловити ОБА числа.
    # Спершу вирізаємо всю фразу, потім шукаємо в ній усі числа біля слова minute.
    pre_clause, pre_ev = _find(
        d, r"incubated with (?:the\s+)?compounds?\s+for\s+"
           r"([^.]{0,140}?)(?=\s+at room temperature|\s+prior to|\s+before|\.)")
    if pre_clause is None:
        pre_clause, pre_ev = _find(
            d, r"(?:pre-?incubat\w*|pre-?treat\w*)\s*(?:with[^.]{0,40})?"
               r"(?:for\s+)?([^.]{0,100}?)(?=\s+prior to|\s+before|\s+then|\.)")

    nums = re.findall(r"(\d+(?:\.\d+)?)\s*(?:min|minute)", pre_clause, re.I) if pre_clause else []
    if not nums and pre_clause:
        nums = re.findall(r"\d+", pre_clause)

    if len(nums) == 1:
        put("compound_preincubation_min", nums[0], pre_ev, "explicit_chembl_description")
        out["compound_preincubation_candidates"] = nums[0]
    elif len(nums) > 1:
        # кілька значень у тексті ("10 minutes or 90 minutes") -> НЕ вгадуємо
        put("compound_preincubation_min", None, pre_ev, "ambiguous")
        out["compound_preincubation_candidates"] = ";".join(nums)
    else:
        put("compound_preincubation_min", None, pre_ev, "unknown")
        out["compound_preincubation_candidates"] = None

    # ---------- число в назві assay (окремо! може означати що завгодно) ----------
    title_num, title_ev = _find(d[:80], r"\((\d+)\s*minutes?\)")
    put("assay_title_time_min", title_num, title_ev,
        "inferred_from_assay_title" if title_num else "unknown")

    # ---------- порядок додавання ----------
    order = None
    if re.search(r"prior to (?:adding|addition of)\s+(?:the\s+)?agonist|"
                 r"before (?:adding|addition of)|"
                 r"compounds were added[^.]{0,120}?(?:then|followed by)[^.]{0,80}"
                 r"(?:cinnamald|AITC|agonist)", d, re.I):
        order = "compound_before_agonist"
    elif re.search(r"after (?:activation|adding agonist)|agonist[^.]{0,60}then[^.]{0,40}compound", d, re.I):
        order = "compound_after_agonist"
    _, ord_ev = _find(d, r"([^.]{0,60}prior to adding agonist[^.]{0,40})")
    put("application_order", order, ord_ev,
        "explicit_chembl_description" if order else "unknown")

    # ---------- обробка кривої ----------
    hill, hill_ev = _find(d, r"Hill coefficient[^.]{0,40}?fixed to\s*(\d+(?:\.\d+)?)")
    put("hill_coefficient_fixed", hill, hill_ev,
        "explicit_chembl_description" if hill else "unknown")

    return out


# =====================================================================
# ЧАСТИНА 2. Відбиток протоколу (protocol fingerprint)
# =====================================================================

FINGERPRINT_FIELDS = [
    "record_species", "construct",
    "assay_modality", "challenge_agonist",
    "agonist_level", "agonist_concentration",
    "cell_line", "application_order", "compound_preincubation_min",
]


def build_fingerprint(row: pd.Series) -> str:
    parts = []
    for f in FINGERPRINT_FIELDS:
        v = row.get(f)
        parts.append("NA" if (v is None or (isinstance(v, float) and np.isnan(v))
                              or str(v) == "nan") else str(v))
    return " | ".join(parts)


# =====================================================================
# ЧАСТИНА 3. Пошук парних контрастів
# =====================================================================

def scaffold_of(smiles: str) -> str | None:
    if not HAVE_RDKIT or not isinstance(smiles, str):
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return None
    # RDKit повертає "" для молекул без кілець. Якщо це залишити, УСІ ациклічні
    # молекули стануть однією штучною "хімічною серією" і статистика поїде.
    if not scaf:
        return "ACYCLIC::" + Chem.MolToSmiles(mol, canonical=True)
    return scaf



def bootstrap_scaffold_median_ci(delta, scaffolds, n_boot=10000, seed=42):
    """
    Bootstrap CI для ТІЄЇ САМОЇ величини що й точкова оцінка:
    median of scaffold medians.

    Раніше була помилка: ресемплили каркаси, але потім зливали ВСІ їхні
    молекули разом і брали медіану молекул. Великий каркас з 20 молекул
    важив більше за каркас з однією -> це інша статистика, ніж
    median(per_scaffold_medians). Тепер ресемплимо самі scaffold-медіани.
    """
    rng = np.random.default_rng(seed)
    sdf = pd.DataFrame({"scaffold": scaffolds, "delta": np.asarray(delta, dtype=float)})
    sdf = sdf[sdf["scaffold"].notna() & sdf["scaffold"].ne("")]
    per_scaffold = sdf.groupby("scaffold")["delta"].median().to_numpy()
    k = len(per_scaffold)
    if k < 3:
        return np.nan, np.nan, "insufficient_scaffolds"
    boot = np.empty(n_boot)
    for i in range(n_boot):
        boot[i] = np.median(rng.choice(per_scaffold, size=k, replace=True))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return float(lo), float(hi), "bootstrap_of_scaffold_medians"


def bootstrap_compound_median_ci(delta, n_boot=10000, seed=42):
    """Простий bootstrap CI для медіани на рівні молекул (для порівняння)."""
    delta = np.asarray(delta, dtype=float)
    n = len(delta)
    if n < 3:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        boot[i] = np.median(delta[rng.integers(0, n, size=n)])
    return tuple(float(x) for x in np.percentile(boot, [2.5, 97.5]))


def paired_stats(delta: np.ndarray, scaffolds: list | None = None) -> dict:
    """
    Статистика парної різниці.
    ВАЖЛИВО:
      * mean має свій CI (t-інтервал), median має свій (bootstrap) — вони РІЗНІ;
      * основна оцінка — median з cluster-bootstrap CI за каркасами,
        бо молекули однієї хімічної серії не є незалежними спостереженнями.
    """
    res = {}
    n = len(delta)
    res["n_compounds"] = n
    res["delta_mean"] = float(np.mean(delta))
    res["delta_median"] = float(np.median(delta))
    res["delta_sd"] = float(np.std(delta, ddof=1)) if n > 1 else np.nan

    frac_pos = float(np.mean(delta > 0))
    res["sign_consistency"] = max(frac_pos, 1.0 - frac_pos)
    res["dominant_direction"] = "B_higher" if frac_pos >= 0.5 else "A_higher"
    res["frac_abs_gt_0.5"] = float(np.mean(np.abs(delta) > 0.5))

    # --- CI для СЕРЕДНЬОГО (t-інтервал, рівень молекул) ---
    if HAVE_SCIPY and n > 2:
        se = res["delta_sd"] / np.sqrt(n)
        lo, hi = sps.t.interval(0.95, n - 1, loc=res["delta_mean"], scale=se)
        res["mean_ci95_low_compound_level"] = float(lo)
        res["mean_ci95_high_compound_level"] = float(hi)
        res["wilcoxon_p_compound_level"] = float(sps.wilcoxon(delta).pvalue)
    else:
        res["mean_ci95_low_compound_level"] = np.nan
        res["mean_ci95_high_compound_level"] = np.nan
        res["wilcoxon_p_compound_level"] = np.nan

    # --- CI для медіани на рівні МОЛЕКУЛ ---
    lo_c, hi_c = bootstrap_compound_median_ci(delta)
    res["median_ci95_low_compound_level"] = lo_c
    res["median_ci95_high_compound_level"] = hi_c

    # --- рівень каркасів ---
    if scaffolds is not None and len(scaffolds) == n:
        sdf = pd.DataFrame({"scaf": scaffolds, "d": delta}).dropna(subset=["scaf"])
        per_scaf = sdf.groupby("scaf")["d"].median()
        k = len(per_scaf)
        res["n_scaffolds"] = k
        res["delta_median_scaffold_level"] = float(per_scaf.median()) if k else np.nan
        res["delta_mean_scaffold_level"] = float(per_scaf.mean()) if k else np.nan
        lo_s, hi_s, meth = bootstrap_scaffold_median_ci(delta, scaffolds)
        res["median_ci95_low_scaffold_level"] = lo_s
        res["median_ci95_high_scaffold_level"] = hi_s
        res["median_ci_method"] = meth
        if HAVE_SCIPY and k > 2:
            sd = float(per_scaf.std(ddof=1))
            se = sd / np.sqrt(k)
            lo2, hi2 = sps.t.interval(0.95, k - 1, loc=per_scaf.mean(), scale=se)
            res["mean_ci95_low_scaffold_level"] = float(lo2)
            res["mean_ci95_high_scaffold_level"] = float(hi2)
            res["wilcoxon_p_scaffold_level"] = float(sps.wilcoxon(per_scaf.values).pvalue)
        else:
            res["mean_ci95_low_scaffold_level"] = np.nan
            res["mean_ci95_high_scaffold_level"] = np.nan
            res["wilcoxon_p_scaffold_level"] = np.nan
    else:
        for kk in ["n_scaffolds", "delta_median_scaffold_level", "delta_mean_scaffold_level",
                   "mean_ci95_low_scaffold_level", "mean_ci95_high_scaffold_level",
                   "median_ci95_low_scaffold_level", "median_ci95_high_scaffold_level",
                   "wilcoxon_p_scaffold_level"]:
            res[kk] = np.nan
        res["median_ci_method"] = "no_scaffolds"

    return res


# =====================================================================
# MAIN
# =====================================================================

def main() -> int:
    if not AUDIT_FILE.exists() or not RAW_FILE.exists():
        print(f"ПОМИЛКА: потрібні файли {AUDIT_FILE} і {RAW_FILE}")
        return 1

    audit = pd.read_csv(AUDIT_FILE)
    raw = pd.read_csv(RAW_FILE)

    # --- ЖОРСТКА ПЕРЕВІРКА МІШЕНІ ---
    # TRPA1 людини і щура фармакологічно РІЗНІ: одна речовина може активувати
    # канал одного виду і блокувати іншого. Змішувати не можна.
    tgt_col = next((c for c in ("target_chembl_id", "target_chembl_id_x")
                    if c in raw.columns), None)
    if tgt_col is not None:
        bad = sorted(set(raw[tgt_col].dropna().unique()) - {EXPECTED_TARGET})
        if bad:
            raise RuntimeError(f"У даних є не-людські мішені TRPA1: {bad}. "
                               f"Очікувався лише {EXPECTED_TARGET}.")
        print(f"  Перевірка мішені: усі записи = {EXPECTED_TARGET} (людський TRPA1) OK")
    else:
        print(f"  УВАГА: колонки target_chembl_id немає у {RAW_FILE}.")
        print(f"  Вид записів приймається як human за конструкцією вибірки")
        print(f"  (дані тягнулись фільтром target_chembl_id={EXPECTED_TARGET}),")
        print(f"  але автоматично це НЕ перевірено.")
    if "assay_organism" in audit.columns:
        orgs = set(audit["assay_organism"].dropna().unique()) - {"Homo sapiens"}
        if orgs:
            raise RuntimeError(f"assay_organism містить не-людські записи: {sorted(orgs)}")

    print("=" * 78)
    print("ПРОТОКОЛЬНА ДЕКОМПОЗИЦІЯ")
    print("=" * 78)
    print(f"Експериментів: {len(audit)}")
    print(f"Вимірів:       {len(raw)}")
    print(f"Молекул:       {raw['inchikey'].nunique()}")
    print(f"RDKit: {'є' if HAVE_RDKIT else 'НЕМАЄ (каркаси не рахуватимуться)'}")
    print(f"SciPy: {'є' if HAVE_SCIPY else 'НЕМАЄ (статистика обмежена)'}")

    # ---- 1. розбір ----
    parsed = [parse_protocol(r.get("description", "")) for _, r in audit.iterrows()]
    pdf = pd.DataFrame(parsed)
    keep = ["assay_chembl_id", "n_compounds", "document_chembl_id", "description"]
    for extra in ("target_chembl_id", "assay_organism", "confidence_score"):
        if extra in audit.columns:
            keep.append(extra)
    dec = pd.concat([audit[keep].reset_index(drop=True),
                     pdf.reset_index(drop=True)], axis=1)

    dec["protocol_fingerprint"] = dec.apply(build_fingerprint, axis=1)

    conf_cols = [c for c in dec.columns if c.endswith("__confidence")]
    dec["n_ambiguous_fields"] = dec[conf_cols].apply(
        lambda r: sum(str(v).startswith("ambiguous") for v in r), axis=1)
    dec["n_unknown_fields"] = dec[conf_cols].apply(
        lambda r: sum(str(v) == "unknown" for v in r), axis=1)
    # Пріоритет ручної перевірки. Відсутність розпізнаного значення (unknown)
    # теж потребує перевірки, а не лише явна суперечність (ambiguous).
    def _priority(r):
        big = r["n_compounds"] >= 50
        amb = r["n_ambiguous_fields"] > 0
        unk = r["n_unknown_fields"]
        if amb and big:
            return "critical"
        if amb or (big and unk >= 3):
            return "high"
        if big or unk >= 4:
            return "medium"
        return "low"
    dec["review_priority"] = dec.apply(_priority, axis=1)
    dec["requires_source_review"] = dec["review_priority"].isin(["critical", "high"])

    print("\n--- РОЗБІР ЗАВЕРШЕНО ---")
    for f in ["assay_modality", "challenge_agonist", "agonist_level",
              "cell_line", "compound_preincubation_min", "application_order"]:
        known = dec[f].notna().sum()
        print(f"  {f:<32} розпізнано у {known}/{len(dec)} експериментах")

    n_amb = int(dec["requires_source_review"].sum())
    print("\n  Пріоритет ручної перевірки:")
    for lvl in ["critical", "high", "medium", "low"]:
        sub = dec[dec.review_priority == lvl]
        if len(sub):
            print(f"    {lvl:<9}: {len(sub):>3} експ. ({int(sub.n_compounds.sum())} молекул-вимірів)")

    dec.to_csv(OUT_DECOMP, index=False)
    print(f"  Збережено: {OUT_DECOMP}")

    # ---- 2. родини протоколів ----
    # УВАГА: це СИГНАТУРИ (однакові за РОЗПІЗНАНИМИ полями), а не підтверджені
    # однакові протоколи. Два assays де все NA потрапляють в одну групу лише тому,
    # що ми однаково мало про них знаємо.
    dec["n_known_signature_fields"] = dec[FINGERPRINT_FIELDS].notna().sum(axis=1)
    fam = dec.groupby("protocol_fingerprint").agg(
        n_assays=("assay_chembl_id", "count"),
        assays=("assay_chembl_id", lambda s: ";".join(s)),
        sum_compounds=("n_compounds", "sum"),
        n_documents=("document_chembl_id", "nunique"),
        n_known_fields=("n_known_signature_fields", "max"),
    ).reset_index().sort_values("sum_compounds", ascending=False)
    fam["signature_information"] = np.where(
        fam["n_known_fields"] >= 5, "well_specified",
        np.where(fam["n_known_fields"] >= 3, "partially_specified", "low_information"))
    fam["verified_protocol_family"] = False   # ставиться ЛИШЕ після ручної перевірки

    # справжня кількість унікальних молекул на родину
    uniq = []
    for _, r in fam.iterrows():
        aids = r["assays"].split(";")
        uniq.append(raw[raw.assay_chembl_id.isin(aids)]["inchikey"].nunique())
    fam["unique_compounds"] = uniq
    fam.to_csv(OUT_FAMILIES, index=False)

    print(f"\n--- СИГНАТУР ПРОТОКОЛУ: {len(fam)} ---")
    print("(сигнатура = збіг за РОЗПІЗНАНИМИ полями; НЕ доказ однакового протоколу)")
    for _, r in fam.head(10).iterrows():
        print(f"  {r['unique_compounds']:>5} мол | {r['n_assays']} експ | {r['n_documents']} док | "
              f"{r['signature_information']:<20} | {r['protocol_fingerprint'][:60]}")
    print(f"  Збережено: {OUT_FAMILIES}")

    # ---- 3. каркаси ----
    scaf_map = {}
    if HAVE_RDKIT:
        print("\nОбчислюю каркаси Мурко...")
        sm = raw.drop_duplicates("inchikey")[["inchikey", "std_smiles"]]
        for _, r in sm.iterrows():
            scaf_map[r["inchikey"]] = scaffold_of(r["std_smiles"])

    # ---- 4. парні контрасти ----
    print("\n--- ПОШУК ПАРНИХ КОНТРАСТІВ ---")
    print("(пари експериментів: той самий документ, спільні молекули,")
    print(" відрізняються рівно одним параметром протоколу)")

    contrasts = []
    for doc, grp in dec.groupby("document_chembl_id"):
        if len(grp) < 2:
            continue
        rows = grp.to_dict("records")
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                A, B = rows[i], rows[j]
                # які поля відрізняються
                diff = []
                for f in FINGERPRINT_FIELDS:
                    va, vb = A.get(f), B.get(f)
                    va = None if (va is None or str(va) == "nan") else str(va)
                    vb = None if (vb is None or str(vb) == "nan") else str(vb)
                    if va != vb:
                        diff.append(f)
                # додатково перевіряємо час у назві (він може бути єдиною відмінністю)
                ta, tb = A.get("assay_title_time_min"), B.get("assay_title_time_min")
                ta = None if (ta is None or str(ta) == "nan") else str(ta)
                tb = None if (tb is None or str(tb) == "nan") else str(tb)
                title_differs = (ta != tb)
                eqa, eqb = A.get("plate_equilibration_min"), B.get("plate_equilibration_min")
                eqa = None if (eqa is None or str(eqa) == "nan") else str(eqa)
                eqb = None if (eqb is None or str(eqb) == "nan") else str(eqb)
                equil_differs = (eqa != eqb)

                if len(diff) > 1:
                    continue  # відрізняються більш ніж одним ключовим полем
                if len(diff) == 0 and not (title_differs or equil_differs):
                    continue  # ідентичні — нема контрасту

                dA = raw[raw.assay_chembl_id == A["assay_chembl_id"]] \
                    .groupby("inchikey")["pchembl_value"].median()
                dB = raw[raw.assay_chembl_id == B["assay_chembl_id"]] \
                    .groupby("inchikey")["pchembl_value"].median()
                common = dA.index.intersection(dB.index)
                if len(common) < MIN_SHARED_COMPOUNDS:
                    continue

                delta = (dB[common] - dA[common]).values
                scafs = [scaf_map.get(k) for k in common] if scaf_map else None
                st = paired_stats(delta, scafs)

                if diff:
                    field = diff[0]
                    va, vb = str(A.get(field)), str(B.get(field))
                    conf = "one_difference_among_extracted_fields__others_unverified"
                elif equil_differs:
                    field, va, vb = "plate_equilibration_min", eqa, eqb
                    conf = "plate_equilibration_differs__other_fields_unverified"
                else:
                    field, va, vb = "assay_title_time_min", ta, tb
                    conf = "assay_title_differs__meaning_unverified"

                rec = {
                    "document_chembl_id": doc,
                    "assay_A": A["assay_chembl_id"],
                    "assay_B": B["assay_chembl_id"],
                    "differing_field": field,
                    "value_A": va,
                    "value_B": vb,
                    "contrast_confidence": conf,
                    "A_preincub_candidates": A.get("compound_preincubation_candidates"),
                    "B_preincub_candidates": B.get("compound_preincubation_candidates"),
                    "A_title_time": ta, "B_title_time": tb,
                    "A_equilibration": eqa, "B_equilibration": eqb,
                    "interpretation_note": (
                        "ЩО САМЕ ЗМІНЮВАЛОСЬ — НЕ ПІДТВЕРДЖЕНО. "
                        "Перевірити першоджерело (патент/стаття) перш ніж називати параметр."
                    ),
                }
                rec.update(st)
                if HAVE_SCIPY and len(common) > 2:
                    rec["spearman_between_conditions"] = float(
                        sps.spearmanr(dA[common], dB[common]).correlation)
                contrasts.append(rec)

    if contrasts:
        cdf = pd.DataFrame(contrasts).sort_values("n_compounds", ascending=False)
        cdf.to_csv(OUT_PAIRED, index=False)
        print(f"\nЗнайдено контрастів: {len(cdf)}")
        print()
        for _, r in cdf.iterrows():
            print(f"  {r['assay_A']} vs {r['assay_B']}  ({r['document_chembl_id']})")
            print(f"    відрізняється: {r['differing_field']} = {r['value_A']} -> {r['value_B']}")
            print(f"    статус:        {r['contrast_confidence']}")
            print(f"    спільних молекул: {r['n_compounds']}, каркасів: {r.get('n_scaffolds')}")
            print(f"    -- рівень МОЛЕКУЛ (аналоги не незалежні -> точність завищена) --")
            print(f"       медіана {r['delta_median']:+.3f}  "
                  f"bootstrap CI[{r['median_ci95_low_compound_level']:+.3f},"
                  f" {r['median_ci95_high_compound_level']:+.3f}]")
            print(f"       середнє {r['delta_mean']:+.3f}  "
                  f"t-CI[{r['mean_ci95_low_compound_level']:+.3f},"
                  f" {r['mean_ci95_high_compound_level']:+.3f}]")
            print(f"    -- рівень КАРКАСІВ (ОСНОВНА ОЦІНКА) --")
            print(f"       медіана scaffold-медіан {r['delta_median_scaffold_level']:+.3f}  "
                  f"CI[{r['median_ci95_low_scaffold_level']:+.3f},"
                  f" {r['median_ci95_high_scaffold_level']:+.3f}]")
            print(f"       метод CI: {r['median_ci_method']}")
            if not np.isnan(r.get('wilcoxon_p_scaffold_level', np.nan)):
                print(f"       Wilcoxon за каркасами p = {r['wilcoxon_p_scaffold_level']:.2e}")
            print(f"    напрямок збігається у {r['sign_consistency']*100:.0f}% молекул "
                  f"({r['dominant_direction']})")
            print()
        print(f"  Збережено: {OUT_PAIRED}")
    else:
        print("\nПарних контрастів за заданими умовами не знайдено.")
        cdf = pd.DataFrame()

    # ---- 5. описова статистика між експериментами (НЕ причинна) ----
    print("\n" + "=" * 78)
    print("ОПИСОВА СТАТИСТИКА МІЖ ЕКСПЕРИМЕНТАМИ")
    print("=" * 78)
    print("УВАГА: різниця медіан між експериментами змішує ефект протоколу")
    print("з ефектом різних хімічних серій. Це НЕ доказ впливу протоколу.")
    per_assay = raw.groupby("assay_chembl_id")["pchembl_value"].agg(["median", "count"])
    big = per_assay[per_assay["count"] >= 20].sort_values("median")
    print(f"\nЕкспериментів з >=20 вимірами: {len(big)}")
    print(f"  діапазон медіан: {big['median'].min():.2f} - {big['median'].max():.2f}"
          f"  (розкид {big['median'].max()-big['median'].min():.2f} log)")

    # ---- 6. metadata ----
    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_assays": int(len(audit)),
        "n_activities": int(len(raw)),
        "n_compounds": int(raw["inchikey"].nunique()),
        "n_protocol_families": int(len(fam)),
        "n_paired_contrasts": int(len(cdf)),
        "n_assays_requiring_source_review": n_amb,
        "min_shared_compounds_threshold": MIN_SHARED_COMPOUNDS,
        "fingerprint_fields": FINGERPRINT_FIELDS,
        "rdkit_available": HAVE_RDKIT,
        "scipy_available": HAVE_SCIPY,
        "caveats": [
            "Скрипт НЕ визначає причину різниць між експериментами.",
            "Значення 'compound_preincubation_min' позначене ambiguous коли "
            "в тексті кілька чисел — воно НЕ вгадується.",
            "Число в назві assay ('15 minutes') збережене окремо і його зміст "
            "НЕ підтверджений — може стосуватись барвника, планшета або сполуки.",
            "Різниця медіан між експериментами не є ефектом протоколу.",
            "Статистика на рівні каркасів чесніша за статистику на рівні молекул.",
        ],
    }
    OUT_META.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n" + "=" * 78)
    print("ГОТОВО")
    print("=" * 78)
    print(f"  {OUT_DECOMP}   — параметри кожного експерименту + докази")
    print(f"  {OUT_FAMILIES} — родини однакових протоколів")
    print(f"  {OUT_PAIRED}   — парні контрасти зі статистикою")
    print(f"  {OUT_META}     — provenance")
    print()
    print("НАСТУПНИЙ КРОК (людина, не скрипт):")
    print(f"  1. Відкрити {OUT_DECOMP}, переглянути колонки *__evidence")
    print(f"     для експериментів де requires_source_review = True")
    print("  2. Для парних контрастів — знайти першоджерело (патент) і")
    print("     підтвердити ЩО САМЕ відрізнялось між умовами")
    print("  3. Тільки після цього називати параметр у статті")
    return 0


if __name__ == "__main__":
    sys.exit(main())