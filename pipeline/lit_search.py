"""
체계적 선행연구 검색 (PubMed E-utilities) — 방법론 신규성 확인용.

목적: '약식동원(food-medicine homology) 화합물 × 공간ABM × 적응요법/공존 ×
myCAF/CAF × PDAC'의 통합 조합을 다룬 선행 논문이 있는지 차원별 불린 쿼리로 확인.

출력:
  docs/literature_search/results.json   — 원자료(쿼리별 count·PMID·요약·초록)
  docs/literature_search/report.md      — 사람이 읽는 신규성 판정 리포트

주: NCBI E-utilities는 무키 3req/s 제한 → sleep 0.4s. 남용 금지.
"""
import sys, os, json, time, urllib.parse, urllib.request
from xml.etree import ElementTree as ET

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "docs", "literature_search")
TOOL = "pdac_novelty_check"
EMAIL = "kusi81kim@gmail.com"   # E-utilities 예의상 식별자(연락처)

# ── 교집합 차원별 쿼리 (각 축이 얼마나 이미 다뤄졌는지 측정) ──────────────
# 핵심: 개별 축은 히트가 많을 것(=기존 분야), '통합' 쿼리는 0에 가까울 것(=신규성).
QUERIES = [
    ("D1_adaptive_natural",
     '("adaptive therapy"[tiab] OR "adaptive dosing"[tiab]) AND '
     '("natural product"[tiab] OR "natural compound"[tiab] OR herbal[tiab] OR '
     'phytochemical[tiab] OR "dietary compound"[tiab]) AND '
     '(cancer[tiab] OR tumor[tiab] OR tumour[tiab])'),

    ("D2_abm_natural_tme",
     '("agent-based model"[tiab] OR "agent based model"[tiab] OR '
     '"individual-based model"[tiab] OR "computational model"[tiab]) AND '
     '(herbal[tiab] OR "natural product"[tiab] OR "natural compound"[tiab] OR '
     'phytochemical[tiab]) AND (tumor[tiab] OR tumour[tiab] OR '
     '"tumor microenvironment"[tiab] OR fibroblast[tiab])'),

    ("D3_netpharm_pdac_caf",
     '("network pharmacology"[tiab]) AND '
     '(pancreatic[tiab] OR PDAC[tiab]) AND '
     '(fibroblast[tiab] OR CAF[tiab] OR stroma[tiab])'),

    ("D4_mycaf_natural_pdac",
     '(myCAF[tiab] OR "myofibroblastic CAF"[tiab] OR '
     '"cancer-associated fibroblast"[tiab] OR "cancer associated fibroblast"[tiab]) AND '
     '(herbal[tiab] OR "natural product"[tiab] OR "natural compound"[tiab] OR '
     'phytochemical[tiab] OR ginsenoside[tiab] OR curcumin[tiab]) AND '
     '(pancreatic[tiab] OR PDAC[tiab])'),

    ("D5_foodmed_cancer_model",
     '("food-medicine homology"[tiab] OR "medicine food homology"[tiab] OR '
     '"medicinal food"[tiab] OR "medicine-food"[tiab]) AND '
     '(cancer[tiab] OR tumor[tiab] OR tumour[tiab]) AND '
     '(model[tiab] OR simulation[tiab] OR "network pharmacology"[tiab] OR '
     'computational[tiab])'),

    ("D6_adaptive_pdac_resistance",
     '("adaptive therapy"[tiab] OR "evolutionary therapy"[tiab]) AND '
     '(pancreatic[tiab] OR PDAC[tiab]) AND '
     '(resistance[tiab] OR "mathematical model"[tiab] OR coexistence[tiab])'),

    ("D7_spatial_abm_pdac",
     '("spatial transcriptomics"[tiab] OR Xenium[tiab] OR "spatial omics"[tiab]) AND '
     '("agent-based"[tiab] OR "agent based"[tiab] OR "computational model"[tiab] OR '
     'simulation[tiab]) AND (pancreatic[tiab] OR PDAC[tiab] OR tumor[tiab])'),

    ("D8_coexistence_natural_control",
     '(coexistence[tiab] OR "tumor control"[tiab] OR "containment"[tiab] OR '
     '"controllability"[tiab]) AND '
     '(herbal[tiab] OR "natural product"[tiab] OR "traditional medicine"[tiab] OR '
     '"dietary"[tiab]) AND (cancer[tiab] OR tumor[tiab] OR tumour[tiab])'),

    # ── 통합 쿼리(신규성의 핵심 증거): 세 축을 동시에 요구 ──
    ("INTEGRATED_A_model_natural_caf",
     '("agent-based"[tiab] OR "computational model"[tiab] OR simulation[tiab] OR '
     '"network pharmacology"[tiab]) AND '
     '(herbal[tiab] OR "natural product"[tiab] OR "traditional medicine"[tiab] OR '
     '"medicinal food"[tiab] OR phytochemical[tiab]) AND '
     '("cancer-associated fibroblast"[tiab] OR CAF[tiab] OR myCAF[tiab] OR stroma[tiab])'),

    ("INTEGRATED_B_adaptive_natural_caf",
     '("adaptive therapy"[tiab] OR coexistence[tiab] OR "tumor control"[tiab]) AND '
     '(herbal[tiab] OR "natural product"[tiab] OR "traditional medicine"[tiab] OR '
     'phytochemical[tiab]) AND '
     '(fibroblast[tiab] OR CAF[tiab] OR stroma[tiab] OR microenvironment[tiab])'),
]

RETMAX = 60


def _get(url):
    req = urllib.request.Request(url, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read()


def esearch(term):
    q = urllib.parse.urlencode({
        "db": "pubmed", "term": term, "retmax": RETMAX, "retmode": "json",
        "sort": "relevance", "tool": TOOL, "email": EMAIL})
    data = json.loads(_get(f"{BASE}/esearch.fcgi?{q}"))
    res = data.get("esearchresult", {})
    return int(res.get("count", 0)), res.get("idlist", [])


def efetch_abstracts(pmids):
    """PMID 목록 → {pmid: {title, journal, year, abstract}}."""
    if not pmids:
        return {}
    q = urllib.parse.urlencode({
        "db": "pubmed", "id": ",".join(pmids), "retmode": "xml",
        "tool": TOOL, "email": EMAIL})
    root = ET.fromstring(_get(f"{BASE}/efetch.fcgi?{q}"))
    out = {}
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//PMID") or "?"
        title = "".join(art.find(".//ArticleTitle").itertext()) \
            if art.find(".//ArticleTitle") is not None else ""
        journal = art.findtext(".//Journal/Title") or ""
        year = art.findtext(".//JournalIssue/PubDate/Year") or \
            art.findtext(".//JournalIssue/PubDate/MedlineDate") or "?"
        abst = " ".join("".join(n.itertext())
                        for n in art.findall(".//Abstract/AbstractText"))
        out[pmid] = {"title": title.strip(), "journal": journal.strip(),
                     "year": year, "abstract": abst.strip()}
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = {}
    all_pmids = {}   # pmid -> set(dimension)
    for name, term in QUERIES:
        try:
            count, ids = esearch(term)
        except Exception as e:
            print(f"[ERR] {name}: {e}")
            results[name] = {"term": term, "error": str(e)}
            time.sleep(0.5)
            continue
        print(f"{name:34s} count={count:5d}  fetched={len(ids)}")
        results[name] = {"term": term, "count": count, "pmids": ids}
        for p in ids:
            all_pmids.setdefault(p, set()).add(name)
        time.sleep(0.45)

    # 초록 수집 (union, 배치 200개씩)
    uniq = list(all_pmids)
    print(f"\n[INFO] unique PMIDs across all queries: {len(uniq)}")
    meta = {}
    for i in range(0, len(uniq), 150):
        batch = uniq[i:i + 150]
        try:
            meta.update(efetch_abstracts(batch))
        except Exception as e:
            print(f"[ERR] efetch batch {i}: {e}")
        time.sleep(0.45)

    for p, m in meta.items():
        m["dimensions"] = sorted(all_pmids.get(p, []))

    with open(os.path.join(OUT_DIR, "results.json"), "w", encoding="utf-8") as f:
        json.dump({"queries": results, "articles": meta},
                  f, ensure_ascii=False, indent=2)
    print(f"[OK] wrote results.json ({len(meta)} articles)")

    # ── 리포트 ──
    lines = ["# PDAC 방법론 선행연구 체계검색 (PubMed)\n",
             f"검색 도구: E-utilities · 쿼리 {len(QUERIES)}개 · "
             f"고유 논문 {len(meta)}편 수집\n",
             "\n## 1. 차원별 히트 수 (개별 축이 얼마나 이미 다뤄졌나)\n",
             "\n| 차원 | 쿼리 성격 | PubMed count |",
             "|---|---|---:|"]
    labels = {
        "D1_adaptive_natural": "적응요법 × 천연물",
        "D2_abm_natural_tme": "ABM/계산모델 × 천연물 × TME",
        "D3_netpharm_pdac_caf": "네트워크약리 × PDAC × CAF",
        "D4_mycaf_natural_pdac": "myCAF × 천연물 × PDAC",
        "D5_foodmed_cancer_model": "약식동원 × 암 × 모델",
        "D6_adaptive_pdac_resistance": "적응요법 × PDAC × 내성",
        "D7_spatial_abm_pdac": "공간전사체 × ABM × PDAC",
        "D8_coexistence_natural_control": "공존/통제 × 천연물",
        "INTEGRATED_A_model_natural_caf": "**통합A: 모델×천연물×CAF**",
        "INTEGRATED_B_adaptive_natural_caf": "**통합B: 적응요법×천연물×CAF**",
    }
    for name, _ in QUERIES:
        r = results.get(name, {})
        c = r.get("count", "ERR")
        lines.append(f"| {name} | {labels.get(name, '')} | {c} |")

    lines.append("\n> 해석: 개별 축(D1~D8)은 히트가 많을수록 '이미 확립된 분야'이며 "
                 "신규성 주장 불가. **통합 쿼리(INTEGRATED_A/B) count가 극히 낮으면** "
                 "= 세 축을 묶은 선행연구가 희소 = **당신 방법론의 신규성 근거**.\n")

    # 통합 쿼리 히트 상세
    for intk in ("INTEGRATED_A_model_natural_caf", "INTEGRATED_B_adaptive_natural_caf"):
        r = results.get(intk, {})
        lines.append(f"\n## 2. {labels.get(intk, intk)} — 히트 상세 "
                     f"(count={r.get('count', '?')})\n")
        ids = r.get("pmids", [])[:20]
        if not ids:
            lines.append("- (히트 없음)")
        for p in ids:
            m = meta.get(p, {})
            lines.append(f"- **{m.get('title', '?')}** "
                         f"({m.get('journal', '?')}, {m.get('year', '?')}) "
                         f"[PMID {p}](https://pubmed.ncbi.nlm.nih.gov/{p}/)")

    # 가장 가까운 이웃(다차원 중복 PMID = 여러 축을 동시에 건드림)
    multi = sorted(meta.items(),
                   key=lambda kv: len(kv[1].get("dimensions", [])), reverse=True)
    lines.append("\n## 3. 가장 가까운 이웃 (여러 축에 동시 히트 = 차별화 대상)\n")
    for p, m in multi[:15]:
        dims = ", ".join(m.get("dimensions", []))
        lines.append(f"\n### {m.get('title', '?')}")
        lines.append(f"*{m.get('journal', '?')} ({m.get('year', '?')})* · "
                     f"[PMID {p}](https://pubmed.ncbi.nlm.nih.gov/{p}/) · 축: {dims}")
        ab = m.get("abstract", "")
        lines.append(f"\n{ab[:600]}{'…' if len(ab) > 600 else ''}\n")

    with open(os.path.join(OUT_DIR, "report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print("[OK] wrote report.md")


if __name__ == "__main__":
    main()
