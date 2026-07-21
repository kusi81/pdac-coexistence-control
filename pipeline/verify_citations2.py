"""라운드2: 1차에서 빈 항목 제목쿼리 재검색 + DOI 정확 추출(참조목록 제외)."""
import sys, os, json, time, urllib.parse, urllib.request
from xml.etree import ElementTree as ET

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "pdac_cite_verify2"
EMAIL = "kusi81kim@gmail.com"

TARGETS = [
    ("3_ohlund_2017", "Öhlund iCAF/myCAF J Exp Med 2017",
     '"Distinct populations of inflammatory fibroblasts and myofibroblastic"[ti]', 2),
    ("6_ozdemir_2014", "Özdemir CAF depletion PDAC 2014",
     '"Depletion of Carcinoma-Associated Fibroblasts and Fibrosis"[ti]', 2),
    ("9b_zhang_prostate_2017", "Zhang 적응요법 전립선 임상 Nat Commun 2017",
     '"Integrating evolutionary dynamics into treatment"[ti]', 2),
    ("2_stroma_review", "PDAC 기질 치료표적 리뷰 (Hosein/Maitra 계열)",
     '(pancreatic[ti]) AND (stroma[ti] OR desmoplas*[ti] OR microenvironment[ti]) '
     'AND review[pt] AND ("Nat Rev"[ta] OR "Nature reviews"[ta])', 5),
]


def _get(url):
    req = urllib.request.Request(url, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read()


def esearch(term, retmax):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax,
                                "retmode": "json", "sort": "relevance",
                                "tool": TOOL, "email": EMAIL})
    return json.loads(_get(f"{BASE}/esearch.fcgi?{q}")
                      ).get("esearchresult", {}).get("idlist", [])


def _doi(art):
    # 논문 자체 DOI만: ELocationID 우선, 없으면 PubmedData/ArticleIdList 첫 doi
    for el in art.findall(".//Article/ELocationID"):
        if el.get("EIdType") == "doi":
            return el.text or ""
    pd = art.find("PubmedData/ArticleIdList")
    if pd is not None:
        for aid in pd.findall("ArticleId"):
            if aid.get("IdType") == "doi":
                return aid.text or ""
    return ""


def efetch(pmids):
    if not pmids:
        return {}
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids),
                                "retmode": "xml", "tool": TOOL, "email": EMAIL})
    root = ET.fromstring(_get(f"{BASE}/efetch.fcgi?{q}"))
    out = {}
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//PMID") or "?"
        t = art.find(".//ArticleTitle")
        title = "".join(t.itertext()) if t is not None else ""
        journal = art.findtext(".//Journal/ISOAbbreviation") or ""
        year = art.findtext(".//JournalIssue/PubDate/Year") or \
            (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "?")[:4]
        vol = art.findtext(".//JournalIssue/Volume") or ""
        issue = art.findtext(".//JournalIssue/Issue") or ""
        pages = art.findtext(".//Pagination/MedlinePgn") or ""
        auth = []
        alist = art.findall(".//AuthorList/Author")
        for a in alist[:3]:
            ln = a.findtext("LastName")
            if ln:
                auth.append(f"{ln} {a.findtext('Initials') or ''}".strip())
        if len(alist) > 3:
            auth.append("et al.")
        out[pmid] = {"pmid": pmid, "title": title.strip(), "journal": journal,
                     "year": year, "vol": vol, "issue": issue, "pages": pages,
                     "authors": ", ".join(auth), "doi": _doi(art)}
    return out


def main():
    for key, desc, term, n in TARGETS:
        print(f"\n## [{key}] {desc}")
        try:
            ids = esearch(term, n)
            time.sleep(0.4)
            meta = efetch(ids)
            time.sleep(0.4)
        except Exception as e:
            print(f"  ERROR: {e}")
            continue
        for pmid in ids:
            m = meta.get(pmid)
            if not m:
                continue
            print(f"  PMID {pmid} | {m['authors']}. {m['title']} "
                  f"{m['journal']} {m['year']};{m['vol']}"
                  f"{'('+m['issue']+')' if m['issue'] else ''}"
                  f"{':'+m['pages'] if m['pages'] else ''}. "
                  f"doi:{m['doi'] or '-'}")


if __name__ == "__main__":
    main()
