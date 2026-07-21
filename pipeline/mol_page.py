"""
🧬 분자 결합 뷰어 — 화학물질이 분자 수준에서 어떤 표적에 어떻게 결합하고,
그 결과 종양 유전자/신호에 어떤 변화를 주는지 시각화.

구성:
  1) 2D 화학구조 (SmilesDrawer, SMILES -> 그림)
  2) 3D 표적 단백질 + 결합부위 (3Dmol.js, RCSB PDB 다운로드)
  3) 결합 기전(잔기·근거) + 하위 유전자 효과 + 시뮬레이션 파라미터 연결

⚠️ 정직성: 천연물-표적 결합 자세는 대부분 실험 공결정이 아니라 도킹(계산) 예측이다.
각 표적에 evidence(experimental/docking/mechanistic)를 표기한다. 3D로 보이는 것은
'표적 단백질 구조'이며, 강조된 잔기는 문헌상 결합부위다(리간드가 그 화합물이 아닐 수
있음). 인터넷 연결이 있어야 PDB/CDN을 불러온다.
"""

import streamlit as st
import streamlit.components.v1 as components

EVIDENCE_BADGE = {
    "experimental": ("🟢 실험 공결정", "#27AE60"),
    "docking":      ("🟠 도킹(계산) 예측", "#E67E22"),
    "mechanistic":  ("🔵 신호전달 기전(구조 미해결)", "#2E86AB"),
}

# --------------------------------------------------------------------------
# 표준화/상용 형태 (참고 정보) — 원물 대신 유효성분을 섭취할 수 있는 상용 제품·의약품.
# ⚠️ 정보 제공일 뿐 복용 권고가 아니며, 항암 치료 중 상호작용 위험이 있어 반드시
#    담당 종양내과의·약사와 상의해야 한다. '의약품'과 '건강기능식품(보충제)'을 구분 표기.
# --------------------------------------------------------------------------
COMMERCIAL_FORMS = {
    "garlic": {
        "type": "보충제(표준화)",
        "form": "숙성마늘추출물(AGE, 예: Kyolic) — S-allylcysteine(SAC) 표준화",
        "why": "생마늘 allicin은 불안정·흡수 변동↑. 숙성 과정이 allicin을 안정·수용성 "
               "SAC로 전환(경구 흡수율 98%+, 혈중 검출 PK 확인). 장용정 allicin 표준화 "
               "제품도 있으나 실제 전신노출은 SAC 쪽이 근거 우세.",
    },
    "mugwort": {
        "type": "의약품",
        "form": "아르테수네이트(Artesunate) — 아르테미시닌의 수용성 유도체(경구·정맥)",
        "why": "항말라리아 승인 의약품이자 종양 임상시험 진행 중(정맥 1상 MTD 18mg/kg 내약 "
               "양호; 대장암 NeoART 2상 200mg/일; HCC 용량증량). 원물 쑥보다 용량·순도 관리 "
               "가능. (아르테미시닌·DHA 보충제도 존재)",
    },
    "ginsenoside_rg3": {
        "type": "의약품(중국 승인)",
        "form": "선이교낭(Shenyi Capsule) — Rg3 단일성분 제제",
        "why": "중국 승인 1류 항암 신약(Z20030044). 폐·유방·소화기암에서 화학요법 **보조제**로 "
               "사용(면역기능 개선·혈관신생 억제 보고). 정확히 라이브러리의 Rg3를 전달.",
    },
    "wild_ginseng": {
        "type": "의약품(중국 승인)",
        "form": "선이교낭(Shenyi Capsule, Rg3) / 표준화 진세노사이드 추출물",
        "why": "산삼 원물 대신 Rg3 단일성분 제제(선이교낭)로 표준화 섭취 가능. 산삼은 "
               "진세노사이드 함량·조성 편차가 크다.",
    },
    "milk_thistle": {
        "type": "의약품/보충제",
        "form": "실리마린 표준화 제제(예: Legalon) — silibinin 함량 표준화",
        "why": "여러 국가에서 간보호 의약품으로 승인·표준화. 원물보다 silibinin 용량 일관.",
    },
    "curcumin": {
        "type": "보충제(고생체이용률)",
        "form": "phytosome형(Meriva)·Theracurmin·BCM-95 등 흡수 개선 표준화 커큐민",
        "why": "일반 커큐민/강황은 경구 흡수가 매우 낮다. 생체이용률 개선 제형이 표준화·"
               "정량 섭취에 유리.",
    },
    "danshen": {
        "type": "의약품(중국 승인)",
        "form": "복방단삼적환(Compound Danshen Dripping Pill) 등 표준화 단삼 제제",
        "why": "단삼 성분(Salvianolic acid B, Tanshinone)을 표준화한 심혈관 의약품. 원물 "
               "대신 정량 섭취 가능(단, 항섬유화/항암은 별개 적응증).",
    },
}

# --------------------------------------------------------------------------
# 분자 라이브러리 (문헌 접지)
# --------------------------------------------------------------------------
MOLECULES = {
    "curcumin": dict(
        label="커큐민 (Curcumin)", evidence="strong",
        pubchem_cid=969516,
        smiles="COC1=CC(=CC=C1O)/C=C/C(=O)CC(=O)/C=C/C2=CC(=C(C=C2)O)OC",
        smiles_ok=True,
        targets=[
            dict(name="STAT3 (SH2 도메인)", pdb="1BG1", evidence="docking",
                 residues=[474, 325, 324, 252, 258, 250],
                 site="SH2 인산펩타이드 결합 포켓",
                 detail="TRP474 (pTyr 인식 핵심 잔기)와 결합, ARG325/GLU324 염다리, "
                        "ILE252/258·ALA250 소수성 접촉 → STAT3 이량체화·인산화 차단"),
            dict(name="COX-2 (PTGS2)", pdb="5IKR", evidence="docking",
                 residues=[], site="사이클로옥시게나제 활성부위",
                 detail="COX-2 활성 억제 → PGE2 생성↓ (CAF 염증성 표현형 완화)"),
        ],
        gene_effects=[("STAT3 표적유전자 (BCL2, MYC, MMP9)", "↓ 하향"),
                      ("COX-2 / PTGS2", "↓ 하향"),
                      ("CCND1 (Cyclin D1)", "↓ 하향"),
                      ("α-SMA (ACTA2, CAF 활성)", "↓ 하향")],
        abm_link="k_caf_activate↓ (CAF 장벽 형성 억제), k_kill↑ (면역 회복)",
    ),
    "gemcitabine": dict(
        label="젬시타빈 (Gemcitabine)", evidence="reference",
        pubchem_cid=60750,
        smiles="C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)(F)F",
        smiles_ok=True,
        targets=[
            dict(name="dCK 활성화효소 (실측)", pdb="1P62", evidence="experimental",
                 residues=[], site="dCK 촉매부위 (젬시타빈 GEO 공결정)",
                 detail="젬시타빈은 전구약물 — dCK(deoxycytidine kinase)가 인산화해 "
                        "활성화한다. 1P62에 젬시타빈이 실제로 결합(GEO, 불소 함유)되어 "
                        "있어 3D에서 약물을 눈으로 확인할 수 있다."),
            dict(name="리보뉴클레오타이드 환원효소 (RRM1)", pdb=None,
                 evidence="mechanistic",
                 residues=[], site="RNR 촉매부위 (인간 공결정 미확보)",
                 detail="활성대사체 dFdCDP가 RNR을 억제 → dNTP 고갈. (효모/대장균 RNR "
                        "공결정은 있으나 인간 RRM1 공결정은 여기 미표시.)"),
            dict(name="DNA (사슬 종결)", pdb=None, evidence="mechanistic",
                 residues=[], site="신생 DNA 가닥",
                 detail="dFdCTP가 DNA에 삽입되어 신장 차단(masked chain termination) "
                        "→ S기 정지·세포사멸"),
        ],
        gene_effects=[("dNTP 풀 (기질 고갈)", "↓ 고갈"),
                      ("DNA 복제 (S기 진행)", "↓ 차단"),
                      ("세포사멸 (apoptosis)", "↑ 유도")],
        abm_link="k_prolif↓ (증식 직접 억제), k_tumor_apoptosis↑ (사멸 유도)",
    ),
    "erlotinib": dict(
        label="에를로티닙 (EGFR 표적 TKI)", evidence="strong",
        pubchem_cid=176870,
        smiles="COCCOC1=CC2=C(C=C1OCCOC)C(=NC=N2)NC1=CC=CC(=C1)C#C",
        smiles_ok=True,
        targets=[
            dict(name="EGFR 인산화효소 도메인", pdb="1M17", evidence="experimental",
                 residues=[769, 766, 721],
                 site="ATP 결합 포켓 (경첩 Met769, 게이트키퍼 Thr766)",
                 detail="quinazoline이 ATP 포켓에 결합, 경첩부 Met769와 수소결합 → EGFR "
                        "자가인산화 차단 → Ras/MAPK·PI3K/Akt↓. **1M17에 erlotinib이 공결정** "
                        "되어 있어 3D에 실제 약물(초록 스틱)이 보임"),
        ],
        gene_effects=[("EGFR 자가인산화", "↓ 차단"), ("Ras/MAPK (ERK)", "↓ 하향"),
                      ("Cyclin D1 (CCND1)", "↓ 하향"),
                      ("젬시타빈-유도 ERBB2 생존신호", "↓ 차단(병용 시너지)")],
        abm_link="k_prolif↓ (증식 억제). 젬시타빈과 병용 시 MAPK 구제신호 차단 → 시너지",
    ),
    "entecavir": dict(
        label="엔테카비르 (B형간염약, 리퍼포징 가설)", evidence="hypothesis",
        pubchem_cid=135398508,
        smiles="C=C1C(CC(C1CO)O)N1C=NC2=C1N=C(N)NC2=O",
        smiles_ok=True,
        targets=[
            dict(name="HBV 중합효소/역전사효소 (실제 표적)", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="HBV pol 활성부위 (항바이러스 — 이것이 실제 표적)",
                 detail="구아노신 유사체가 HBV 중합효소를 억제해 바이러스 복제 차단. "
                        "간섬유화 되돌림은 이 **항바이러스 효과로 원인(HBV)이 제거된 "
                        "간접 결과**이지 성상세포/myCAF 직접 작용이 아님"),
            dict(name="KDM5B 탈메틸화효소 (항암 가설)", pdb="5A1F",
                 evidence="docking", residues=[],
                 site="KDM5B 촉매부위 (in silico 도킹)",
                 detail="in silico 연구에서 ETV가 KDM5B(여러 종양 과발현)를 억제해 증식↓·"
                        "세포사멸↑ 가설 (Pharmaceuticals 2023). **PDAC 근거 없음**"),
        ],
        gene_effects=[("HBV DNA (항바이러스)", "↓ 억제"),
                      ("KDM5B (가설 표적)", "↓ 가설"),
                      ("myCAF/PDAC 직접 효과", "확립된 근거 없음")],
        abm_link="k_prolif 약간↓ (KDM5B 가설). 섬유화 효과는 항바이러스 간접이라 "
                 "PDAC(비바이러스성)엔 해당 없음",
    ),
    "ginsenoside_rg3": dict(
        label="진세노사이드 Rg3", evidence="strong",
        pubchem_cid=9918693,
        smiles="",             # 대형 트리테르페노이드 배당체 — 2D 생략(PubChem 참조)
        smiles_ok=False,
        targets=[
            dict(name="PD-L1 (당화 조절)", pdb="5C3T", evidence="mechanistic",
                 residues=[], site="PD-L1 (직접 공결정 없음)",
                 detail="EGFR 인산화↓ → p-GSK3β↓ → PD-L1 당화↓·분해 촉진 "
                        "→ 표면 PD-L1↓ → T세포 세포독성 회복"),
            dict(name="EGFR 인산화효소", pdb="2ITY", evidence="mechanistic",
                 residues=[], site="EGFR 키나제 도메인",
                 detail="EGFR 신호 하향 → 하위 PD-L1 당화 경로 차단"),
        ],
        gene_effects=[("CD274 (PD-L1) 표면발현", "↓ 하향"),
                      ("PI3K/Akt/mTOR", "↓ 하향"),
                      ("Granzyme B / Perforin (CD8)", "↑ 상향")],
        abm_link="k_kill↑ (PD-L1↓로 면역살상 회복), k_prolif↓ (Akt/mTOR↓)",
    ),
    # ==================== 전통 약재 (근거 수준 상이) ====================
    "garlic": dict(
        label="마늘 (Allicin / DATS)", evidence="strong",
        pubchem_cid=65036,
        smiles="C=CCSS(=O)CC=C",   # 알리신
        smiles_ok=True,
        targets=[
            dict(name="NF-κB / PI3K-Akt (다중 표적)", pdb="1VKX",
                 evidence="mechanistic", residues=[],
                 site="티올 반응성 다중 표적 (특정 공결정 없음)",
                 detail="Allicin·DATS의 반응성 황이 단백질 티올과 반응 → NF-κB·"
                        "PI3K/Akt·MAPK 신호 억제, Bcl-2↓·caspase↑ 세포사멸. "
                        "3D는 NF-κB(p50/p65) 참고 구조일 뿐 결합부위 아님"),
        ],
        gene_effects=[("Bcl-2", "↓ 하향"), ("Caspase-3/9", "↑ 상향"),
                      ("NF-κB 활성", "↓ 하향")],
        abm_link="k_prolif↓ (증식 억제), k_tumor_apoptosis↑ (세포사멸)",
    ),
    "mugwort": dict(
        label="쑥 (Artesunate / DHA)", evidence="strong",
        pubchem_cid=68827,        # 아르테미시닌
        smiles="", smiles_ok=False,   # 엔도퍼옥사이드 — 구조는 PubChem 참조
        targets=[
            dict(name="TGF-β / Smad (CAF 비활성화)", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="TGF-β 신호 경로 (구조 미해결)",
                 detail="아르테수네이트/DHA가 TGF-β 신호 억제 → CAF를 활성→휴지 "
                        "비활성화, 면역억제(TGF-β1·IL-10)↓ (JECCR 2018)"),
        ],
        gene_effects=[("TGF-β1", "↓ 하향"), ("α-SMA (CAF 활성)", "↓ 하향"),
                      ("IL-10 (면역억제)", "↓ 하향")],
        abm_link="k_caf_activate↓ (장벽 억제), k_kill↑, k_prolif↓",
    ),
    "wild_ginseng": dict(
        label="산삼 (Ginsenosides)", evidence="strong",
        pubchem_cid=9918693,
        smiles="", smiles_ok=False,   # 대형 배당체
        targets=[
            dict(name="TGF-β1 / Smad (성상세포 HSC)", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="HSC 활성화 경로 (구조 미해결)",
                 detail="진세노사이드가 성상세포 활성 억제 → 콜라겐·TGF-β1↓ "
                        "(Smad2/3/4 회복), MAPK/NF-κB↓ 세포사멸, Rg3 PD-L1↓"),
        ],
        gene_effects=[("TGF-β1 / 콜라겐", "↓ 하향"),
                      ("PD-L1 (CD274)", "↓ 하향"), ("MAPK / NF-κB", "↓ 하향")],
        abm_link="k_caf_activate↓ (항섬유화), k_kill↑, k_prolif↓",
    ),
    "platycodon": dict(
        label="도라지 (Platycodin D)", evidence="moderate",
        pubchem_cid=162859,
        smiles="", smiles_ok=False,   # 트리테르페노이드 사포닌
        targets=[
            dict(name="PD-1 / STAT3-VEGF-A (CD8)", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="CD8 T세포 PD-1 조절 (구조 미해결)",
                 detail="종양 STAT3-VEGF-A↓ → CD8 표면 PD-1↓ → 침윤·살상↑; "
                        "PUMA↑ 세포사멸 (Front Pharmacol 2022)"),
        ],
        gene_effects=[("PD-1 (CD8 표면)", "↓ 하향"), ("PUMA", "↑ 상향"),
                      ("Bcl-2", "↓ 하향")],
        abm_link="k_kill↑ (PD-1↓), cd8_recruit↑ (침윤↑), k_prolif↓",
    ),
    "sea_cucumber": dict(
        label="해삼 (Frondoside A)", evidence="moderate",
        pubchem_cid=11969532,
        smiles="", smiles_ok=False,   # 트리테르펜 배당체
        targets=[
            dict(name="p53 / caspase 세포사멸", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="세포내 세포사멸 경로 (구조 미해결)",
                 detail="p53↑ → caspase 활성 세포사멸, S/G2M 정지, "
                        "EGFR/Akt/ERK/FAK/MMP-9↓ 전이 억제 (Mar Drugs 2015)"),
        ],
        gene_effects=[("p53", "↑ 상향"), ("Caspase-3/9", "↑ 상향"),
                      ("MMP-9 (전이)", "↓ 하향")],
        abm_link="k_prolif↓ (증식 억제), k_tumor_apoptosis↑ (세포사멸)",
    ),
    "deer_antler": dict(
        label="녹용 (Velvet peptides)", evidence="weak",
        pubchem_cid=None,             # 펩타이드 혼합물 — 단일 구조 없음
        smiles="", smiles_ok=False,
        targets=[
            dict(name="TGF-β / NF-κB (항섬유화) · NK 활성", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="펩타이드 혼합물 — 단일 표적 아님",
                 detail="NK 활성·조혈 촉진(면역지원), TGF-β/NF-κB 항섬유화. "
                        "**직접 항암 근거는 약함 — 보조적** (PMC11125008)"),
        ],
        gene_effects=[("NK 세포 활성", "↑ 상향"), ("TGF-β (섬유화)", "↓ 하향")],
        abm_link="cd8_recruit↑ (면역지원), k_caf_activate↓ (보조적)",
    ),
    # ==================== 간섬유화 항섬유화 약재 (HCC 맥락) ====================
    "danshen": dict(
        label="단삼 (Salvianolic acid B / Tanshinone IIA)", evidence="strong",
        pubchem_cid=11629084,     # Salvianolic acid B
        smiles="", smiles_ok=False,   # 대형 폴리페놀
        targets=[
            dict(name="TGF-β1 / p38·ERK MAPK (성상세포)", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="HSC 활성화 경로 (구조 미해결)",
                 detail="Salvianolic acid B가 TGF-β1·p38/ERK MAPK 억제로 HSC 활성 차단; "
                        "Tanshinone IIA는 활성 HSC 세포사멸 유도 (PMC9385955)"),
        ],
        gene_effects=[("α-SMA / 콜라겐 I", "↓ 하향"), ("TGF-β1", "↓ 하향"),
                      ("활성 HSC 생존", "↓ 세포사멸")],
        abm_link="k_caf_activate↓ (섬유화 억제) — HCC 맥락에서 전암 토양 약화",
    ),
    "milk_thistle": dict(
        label="밀크씨슬 (Silymarin / Silibinin)", evidence="strong",
        pubchem_cid=31553,        # Silibinin
        smiles="", smiles_ok=False,   # 플라보노리그난 (복잡)
        targets=[
            dict(name="HSC 세포주기 / Akt (성상세포)", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="HSC 증식·활성 경로 (구조 미해결)",
                 detail="Silibinin이 α-SMA↓, p27/p53↑·Akt↓로 HSC 세포주기 정지·증식 억제; "
                        "HCC 세포 증식도 억제 (PMC5052367)"),
        ],
        gene_effects=[("α-SMA (HSC 활성)", "↓ 하향"), ("p27 / p53", "↑ 상향"),
                      ("Akt", "↓ 하향")],
        abm_link="k_caf_activate↓, k_prolif↓ — 항섬유화 + 항증식",
    ),
    "astragaloside": dict(
        label="황기 (Astragaloside IV)", evidence="strong",
        pubchem_cid=13943297,     # Astragaloside IV
        smiles="", smiles_ok=False,   # 대형 배당체
        targets=[
            dict(name="TGF-β1 / Smad (성상세포)", pdb=None,
                 evidence="mechanistic", residues=[],
                 site="TGF-β/Smad 경로 (구조 미해결)",
                 detail="Astragaloside IV가 TGF-β1·p-Smad2/3↓, Smad7↑로 HSC 활성·"
                        "콜라겐 형성 억제 (PMC5952439)"),
        ],
        gene_effects=[("TGF-β1 / p-Smad2/3", "↓ 하향"), ("Smad7 (억제자)", "↑ 상향"),
                      ("콜라겐 / MMP-2", "↓ 조절")],
        abm_link="k_caf_activate↓ (TGF-β/Smad 차단) — HCC 항섬유화",
    ),
}

# 물질(성분)-수준 근거 배지 (표적-수준 evidence와 별개)
SUBSTANCE_EVIDENCE = {
    "strong":     ("🟢 강함 (다수 연구)", "#27AE60"),
    "moderate":   ("🟡 보통 (일부 연구)", "#F1C40F"),
    "weak":       ("🟠 약함 (예비·보조)", "#E67E22"),
    "reference":  ("⚪ 기준/벤치마크", "#95A5A6"),
    "hypothesis": ("🔴 가설 (미검증)", "#C0392B"),
}


def _mol2d_html(smiles):
    return f"""
<div style="background:#fff;border-radius:8px;padding:6px;">
  <canvas id="mol2d" style="width:100%;max-width:420px;height:300px;"></canvas>
</div>
<script src="https://unpkg.com/smiles-drawer@2.0.3/dist/smiles-drawer.min.js"></script>
<script>
  var opts = {{width:420, height:300, bondThickness:1.1}};
  var sd = new SmilesDrawer.Drawer(opts);
  SmilesDrawer.parse("{smiles}", function(tree) {{
      sd.draw(tree, "mol2d", "light", false);
  }}, function(err) {{
      document.getElementById("mol2d").insertAdjacentHTML(
        "afterend", "<p style='color:#C0392B'>2D 구조 렌더 실패</p>");
  }});
</script>
"""


def _mol3d_html(pdb, residues):
    resi = ",".join(str(r) for r in residues) if residues else ""
    resi_block = f"""
      v.setStyle({{resi:[{resi}]}},
                 {{stick:{{colorscheme:'redCarbon'}}, cartoon:{{color:'spectrum'}}}});
      v.addResLabels({{resi:[{resi}]}}, {{fontSize:11, showBackground:true}});
    """ if residues else ""
    zoom_block = f"v.zoomTo({{resi:[{resi}]}}); v.zoom(0.7);" if residues else "v.zoomTo();"
    return f"""
<div id="mol3d" style="width:100%;height:460px;position:relative;
     border:1px solid #ddd;border-radius:8px;"></div>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
  var v = $3Dmol.createViewer("mol3d", {{backgroundColor:"white"}});
  $3Dmol.download("pdb:{pdb}", v, {{}}, function() {{
      v.setStyle({{}}, {{cartoon:{{color:"spectrum"}}}});
      v.setStyle({{resn:"HOH"}}, {{}});                    // 물 숨김
      // 결합 리간드(hetero, 물 제외)를 초록 스틱으로
      v.setStyle({{hetflag:true, resn:"HOH", invert:true}},
                 {{stick:{{colorscheme:"greenCarbon"}}}});
      {resi_block}
      {zoom_block}
      v.render();
  }});
</script>
"""


def render():
    st.sidebar.markdown("---")
    st.sidebar.markdown("**분자 선택**")
    mkey = st.sidebar.selectbox(
        "molecule", list(MOLECULES.keys()),
        format_func=lambda k: MOLECULES[k]["label"])
    mol = MOLECULES[mkey]
    tnames = [t["name"] for t in mol["targets"]]
    tsel = st.sidebar.selectbox("표적 단백질", tnames)
    target = next(t for t in mol["targets"] if t["name"] == tsel)

    # ---- header ----
    st.title("🧬 분자 결합 뷰어")
    ev_txt, ev_col = SUBSTANCE_EVIDENCE.get(
        mol.get("evidence", "reference"), ("", "#999"))
    pubchem = (f" · [PubChem CID {mol['pubchem_cid']}]"
               f"(https://pubchem.ncbi.nlm.nih.gov/compound/{mol['pubchem_cid']})"
               if mol.get("pubchem_cid") else " · (단일 화합물 아님)")
    st.markdown(
        f"**{mol['label']}** 이(가) 분자 수준에서 표적에 어떻게 결합하고, "
        "그 결과 종양 유전자·신호에 어떤 변화를 주는지 봅니다.")
    st.markdown(
        f"<span style='background:{ev_col}22;border-left:4px solid {ev_col};"
        f"padding:.2rem .55rem;border-radius:4px'>항암 근거 수준: "
        f"<b>{ev_txt}</b></span>{pubchem}", unsafe_allow_html=True)
    if mol.get("evidence") == "weak":
        st.warning("⚠️ 이 물질은 **직접 항암 근거가 약합니다**(주로 면역지원·항섬유화 보조). "
                   "종양 억제 효과로 해석하지 마세요.")
    elif mol.get("evidence") == "hypothesis":
        st.error("🔴 **미검증 가설입니다.** 이 약물의 췌장암/myCAF 직접 작용은 확립된 근거가 "
                 "없습니다. 아래 표적은 '실제 표적(다른 적응증)'과 '가설 표적'을 구분해 "
                 "표기하며, myCAF와의 직접 결합구조는 존재하지 않습니다.")

    st.info("ℹ️ 3D는 **표적 단백질 구조**(RCSB PDB)이며 강조 잔기는 문헌상 결합부위입니다. "
            "천연물 결합 자세는 대부분 도킹(계산) 예측이라, 각 표적에 근거 수준을 표기합니다. "
            "PDB·구조 라이브러리를 불러오려면 인터넷 연결이 필요합니다.")

    left, right = st.columns([2, 3])

    # ---- 2D structure ----
    with left:
        st.subheader("① 2D 화학구조")
        if mol["smiles_ok"] and mol["smiles"]:
            components.html(_mol2d_html(mol["smiles"]), height=340)
            st.caption(f"SMILES: `{mol['smiles']}`")
        else:
            st.warning("대형 배당체 구조 — 2D 렌더 생략. PubChem에서 확인하세요.")

    # ---- 3D binding site ----
    with right:
        st.subheader("② 3D 표적 + 결합부위")
        badge, color = EVIDENCE_BADGE[target["evidence"]]
        pdb_txt = f"PDB `{target['pdb']}`" if target["pdb"] else "구조 미해결"
        st.markdown(
            f"<span style='background:{color}22;border-left:4px solid {color};"
            f"padding:.2rem .5rem;border-radius:4px'><b>{tsel}</b> · {pdb_txt} · "
            f"{badge}</span>", unsafe_allow_html=True)
        if target["pdb"]:
            components.html(_mol3d_html(target["pdb"], target["residues"]), height=480)
        else:
            st.info("이 표적은 실험 구조가 없어 3D를 표시하지 않습니다 (신호전달 기전).")

    # ---- 구조 해석 가이드 ----
    with st.expander("🧭 3D 구조 해석 가이드 — 무엇을 읽고, 무엇은 결론 내리면 안 되는지"):
        st.markdown("""
**화면 읽는 법**

| 요소 | 의미 |
|---|---|
| 리본(무지개색) | 단백질 뼈대 접힘 — N말단(파랑)→C말단(빨강) |
| 🟥 빨간 스틱 | 강조한 **결합부위 잔기**(문헌/도킹상 약물이 닿는 곳) |
| 🟩 초록 스틱 | 그 PDB에 **공결정된 리간드** 분자 |

⚠️ **초록 스틱 주의**: 🟢**실험** 배지일 때만 초록 = *우리 화합물*(에를로티닙·젬시타빈).
🟠도킹·🔵기전 배지에선 초록은 그 구조에 **우연히 들어있는 참고 분자**입니다 — 예: COX-2(5IKR)엔
메페남산, EGFR(2ITY)엔 게피티닙이 붙어 있고 **커큐민·Rg3가 아닙니다**. (참고 분자는 결합
주머니 위치를 알려주는 용도로만 보세요.)

마우스로 회전·확대해 강조 잔기가 **움푹 파인 주머니(pocket)** 안인지 보세요 = 실제 결합 자리.

**✅ 내릴 수 있는 결론**
- 약물이 단백질의 **어디**에 작용하나(ATP 주머니·SH2·촉매부위 등) → 기전의 *종류*
- druggable **주머니**가 실제 있나
- (실험 공결정만) **진짜 결합 자세**·접촉 잔기 — 실측 증거
- 결합부위가 기능부위와 겹치면 억제가 **구조적으로 타당**

**❌ 내리면 안 되는 결론**
- **결합 세기·필요 용량**(Kd/IC50) — 구조는 안 알려줌. '붙는 그림'≠'세게 붙음'
- **세포·환자 효과** — 좋은 자세 ≠ 실제 작동(생체이용률·off-target 별개)
- **선택성**(다른 단백질에도 붙는지)
- 🟠**도킹**이면 강조 잔기는 *계산 예측*일 뿐, 🔵**기전**이면 애초에 '여기 붙는다' 주장 불가

**핵심: 위 ② 배지가 신뢰도를 가릅니다** — 🟢 실험이면 자세를 믿고, 🟠 도킹이면 가설,
🔵 기전이면 결합 주장이 아니라 경로 작용까지만.
""")

    # ---- binding mechanism ----
    st.subheader("③ 결합 기전")
    st.markdown(f"**결합 부위:** {target['site']}")
    st.markdown(f"{target['detail']}")
    if target["residues"]:
        st.caption("강조 잔기(빨강 스틱): " +
                   ", ".join(str(r) for r in target["residues"]) +
                   " — 문헌상 결합/상호작용 잔기 (도킹 예측인 경우 계산 결과).")

    # ---- downstream gene effects ----
    st.subheader("④ 하위 종양 유전자·신호 변화")
    import pandas as pd
    df = pd.DataFrame(mol["gene_effects"], columns=["유전자 / 경로", "변화"])
    st.dataframe(df, hide_index=True, use_container_width=True)

    # ---- loop back to simulation ----
    st.subheader("⑤ → 교란 시뮬레이션 파라미터 연결")
    st.success(f"이 분자 효과는 ABM에서 **{mol['abm_link']}** 로 인코딩됩니다. "
               "'💊 교란 시뮬레이션' 모드에서 조직 수준 반응(종양 궤적)으로 확인하세요.")
    st.caption("분자(결합) → 신호/유전자 → 세포·조직(ABM) → 종양 성장으로 스케일이 "
               "연결됩니다. 단, 분자 결합 자세와 물질→파라미터 매핑은 모두 문헌 기반 "
               "가설이며 실측 검증이 필요합니다.")

    # ---- 표준화/상용 형태 (참고) ----
    cf = COMMERCIAL_FORMS.get(mkey)
    if cf:
        st.subheader("⑥ 표준화·상용 형태 (참고)")
        st.markdown(f"**[{cf['type']}]** {cf['form']}")
        st.caption(cf["why"])
        st.error("⚠️ **정보 제공일 뿐 복용 권고가 아닙니다.** 원물(생마늘·생쑥 등)을 그대로 "
                 "먹기 어려울 때 유효성분을 표준화해 섭취하는 상용 형태를 '참고'로 정리한 "
                 "것입니다. **항암 치료 중에는 이런 성분이 항암제와 상호작용**(예: 마늘↔"
                 "항응고·CYP효소, 아르테수네이트·화학요법의 골수억제 중첩)할 수 있어, "
                 "복용 여부·용량은 **반드시 담당 종양내과의·약사와 상의**해야 합니다. "
                 "임상시험 중이라는 것은 효과가 입증됐다는 뜻이 아닙니다.")
