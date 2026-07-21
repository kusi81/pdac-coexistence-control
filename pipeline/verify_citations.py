"""
Introduction 파운데이션 인용의 PMID·서지 확정 (PubMed E-utilities).
각 논문을 정밀 쿼리로 검색 → 상위 후보의 저자·저널·연도·권·페이지·PMID·DOI 출력.
사람이 후보 중 정답을 골라 introduction_draft.md에 반영.
"""
import sys, os, json, time, urllib.parse, urllib.request
from xml.etree import ElementTree as ET

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "pdac_cite_verify"
EMAIL = "kusi81kim@gmail.com"
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "docs", "manuscript")

# ref키 : (설명, 정밀쿼리, 후보수)
TARGETS = [
    ("1_siegel_stats",
     "Siegel 최신 암통계 (췌장암 생존율 인용원)",
     'Siegel RL[au] AND "Cancer statistics"[ti] AND "CA Cancer J Clin"[ta]', 3),
    ("2_pdac_stroma_review",
     "PDAC desmoplasia/stroma 개관 리뷰",
     '(desmoplasia[ti] OR stroma[ti]) AND (pancreatic[ti] OR PDAC[ti]) AND '
     'review[pt] AND microenvironment[tiab]', 4),
    ("3_ohlund_2017",
     "Öhlund iCAF/myCAF J Exp Med 2017",
     'Ohlund D[au] AND fibroblast[tiab] AND pancreatic[tiab] AND '
     '("J Exp Med"[ta]) AND 2017[dp]', 3),
    ("4_biffi_2019",
     "Biffi IL-1 JAK/STAT CAF heterogeneity Cancer Discov 2019",
     'Biffi G[au] AND (fibroblast[tiab] OR CAF[tiab]) AND pancreatic[tiab] AND '
     '2019[dp]', 3),
    ("5_rhim_2014",
     "Rhim stromal elements restrain PDAC Cancer Cell 2014",
     'Rhim AD[au] AND stroma*[tiab] AND (pancreatic[tiab] OR PDAC[tiab]) AND '
     '2014[dp]', 3),
    ("6_ozdemir_2014",
     "Özdemir CAF depletion PDAC reduced survival Cancer Cell 2014",
     'Ozdemir BC[au] AND fibroblast[tiab] AND (pancreatic[tiab] OR PDAC[tiab]) AND '
     '2014[dp]', 3),
    ("9_adaptive_evol",
     "적응요법 진화동역학 임상/모델 (Zhang 등)",
     '"adaptive therapy"[tiab] AND (evolutionary[tiab] OR prostate[tiab] OR '
     'clinical[tiab]) AND (Zhang[au] OR Gatenby[au])', 5),
    ("10_intermittent_npj_2024",
     "intermittent vs continuous adaptive dosing npj Syst Biol Appl 2024",
     '(intermittent[tiab] AND continuous[tiab] AND adaptive[tiab]) AND '
     '(dosing[tiab] OR chemotherapy[tiab]) AND cancer[tiab]', 4),
    ("12_healthstrengthening_2018",
     "health-strengthening herbal network pharmacology cancer Cancers 2018",
     '"health-strengthening"[tiab] AND "network pharmacology"[tiab] AND cancer[tiab]', 3),
    ("13_medfood_anticancer_review",
     "medicine-food homology 항암 리뷰 (Hawthorn/Zanthoxylum)",
     '("medicine food homology"[tiab] OR "homology of medicine and food"[tiab]) AND '
     '(anticancer[tiab] OR antitumor[tiab] OR cancer[tiab]) AND review[pt]', 4),
    ("17_multiscale_abm_2019",
     "Multiscale agent-based hybrid modeling tumor immune microenvironment",
     '"agent-based"[tiab] AND (multiscale[tiab] OR "multi-scale"[tiab]) AND '
     '"tumor immune"[tiab]', 4),
    ("18_spatial_pdac_model",
     "공간전사체 접지 PDAC 계산/ABM 모델 대표",
     '("spatial transcriptomics"[tiab] OR "spatial omics"[tiab]) AND '
     '(pancreatic[tiab] OR PDAC[tiab]) AND (model[tiab] OR "agent-based"[tiab] OR '
     'simulation[tiab] OR computational[tiab])', 5),
]


def _get(url):
    req = urllib.request.Request(url, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read()


def esearch(term, retmax):
    q = urllib.parse.urlencode({
        "db": "pubmed", "term": term, "retmax": retmax, "retmode": "json",
        "sort": "relevance", "tool": TOOL, "email": EMAIL})
    data = json.loads(_get(f"{BASE}/esearch.fcgi?{q}"))
    return data.get("esearchresult", {}).get("idlist", [])


def efetch(pmids):
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
        journal = art.findtext(".//Journal/ISOAbbreviation") or \
            art.findtext(".//Journal/Title") or ""
        year = art.findtext(".//JournalIssue/PubDate/Year") or \
            (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "?")[:4]
        vol = art.findtext(".//JournalIssue/Volume") or ""
        issue = art.findtext(".//JournalIssue/Issue") or ""
        pages = art.findtext(".//Pagination/MedlinePgn") or ""
        authors = []
        for a in art.findall(".//AuthorList/Author")[:3]:
            ln = a.findtext("LastName")
            ini = a.findtext("Initials") or ""
            if ln:
                authors.append(f"{ln} {ini}")
        if art.find(".//AuthorList/Author") is not None and \
                len(art.findall(".//AuthorList/Author")) > 3:
            authors.append("et al.")
        doi = ""
        for aid in art.findall(".//ArticleIdList/ArticleId"):
            if aid.get("IdType") == "doi":
                doi = aid.text or ""
        out[pmid] = {"pmid": pmid, "title": title.strip(), "journal": journal,
                     "year": year, "vol": vol, "issue": issue, "pages": pages,
                     "authors": ", ".join(authors), "doi": doi}
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    report = ["# 파운데이션 인용 검증 결과\n",
              "각 항목의 후보 중 제목·저자·연도가 맞는 것을 골라 "
              "introduction_draft.md References에 반영.\n"]
    allcand = {}
    for key, desc, term, n in TARGETS:
        try:
            ids = esearch(term, n)
            time.sleep(0.4)
            meta = efetch(ids)
            time.sleep(0.4)
        except Exception as e:
            report.append(f"\n## [{key}] {desc}\n- ERROR: {e}")
            continue
        allcand[key] = meta
        report.append(f"\n## [{key}] {desc}")
        report.append(f"`{term}`\n")
        if not meta:
            report.append("- (후보 없음)")
        # esearch 순서 유지
        for pmid in ids:
            m = meta.get(pmid)
            if not m:
                continue
            cite = (f"{m['authors']}. {m['title']} *{m['journal']}* "
                    f"{m['year']};{m['vol']}"
                    f"{'('+m['issue']+')' if m['issue'] else ''}"
                    f"{':'+m['pages'] if m['pages'] else ''}.")
            report.append(f"- **PMID {pmid}** — {cite}"
                          + (f" doi:{m['doi']}" if m['doi'] else ""))
        print(f"{key:32s} candidates={len(meta)}")

    with open(os.path.join(OUT_DIR, "citation_candidates.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(report))
    with open(os.path.join(OUT_DIR, "citation_candidates.json"), "w",
              encoding="utf-8") as f:
        json.dump(allcand, f, ensure_ascii=False, indent=2)
    print("[OK] wrote citation_candidates.md / .json")


if __name__ == "__main__":
    main()
