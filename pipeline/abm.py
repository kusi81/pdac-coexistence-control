"""
공간 agent-based model (ABM) — 종양-CAF-면역 동역학 + 물질 교란 시뮬레이션.

경로 B: 정적 합성 생성기(synthetic.py)를 시간에 따라 변하는 규칙 기반 모델로 확장.
매 스텝마다 세포가 증식/이동/사멸하고, 그 결과 조직 지도가 진화한다. 최종/중간
스냅샷은 기존 대시보드 지표(장벽 점수 등)로 그대로 판독한다.

핵심 기전 (문헌 접지 — 아래 REFERENCES)
------------------------------------------------
1. 종양 증식: 국소 밀도가 수용력 미만이면 k_prolif 확률로 분열.
2. myCAF 장벽: 종양 병소 가장자리에서 기질세포가 k_caf 확률로 myCAF로 활성화되어
   링을 형성/두껍게 함. (TGF-β 구동 myCAF; Grauel 2020, Dominguez 2020)
3. CD8 이동: 가장 가까운 종양으로 이동하되, 경로 상 myCAF 회랑 밀도에 비례해
   **감쇠**된다 -> 두꺼운 장벽은 CD8를 물리적으로 가둔다.
   (myCAF/LRRC15가 CD8를 주변부에 corral: Cancer Res 2022; SCOTIA Nat Genet 2024)
4. CD8 살상: kill_radius 안의 종양을 k_kill 확률로 사멸.
   (종양 내 CTL ~1 kill / 6h; Boissonnas 2007. 생체 살상력 제한: Halle 2016)
5. CD8 유입/소멸: 주변부에서 신규 유입, 일부 turnover.

물질(한약 성분 등) = 위 파라미터에 대한 **배수 벡터**. 조합은 배수를 곱해 합성.
용량(dose 0~1)은 효과 강도를 선형 보간한다.

정직한 경고
-----------
이 모델은 **가설 탐색 샌드박스**이지 실제 반응 예측이 아니다. 파라미터는 문헌에서
방향과 자릿수만 접지했을 뿐 특정 환자/약물에 피팅되지 않았다. 물질->배수 매핑은
약리 문헌 기반 '가설'이며, 실측 용량-반응으로 검증(경로 A)해야 한다.

REFERENCES (기전 방향의 근거)
- Boissonnas et al. 2007 JCI — 종양 내 CTL 살상률 ~1 cell/6h.
- Halle et al. 2016 Immunity — 생체 CTL 살상력은 제한적·협동적.
- Cancer Res 2022 82(16):2904 — CAF가 CD8 침윤 억제, ICB 저항 부여.
- Dominguez 2020 / Grauel 2020 — TGF-β가 myCAF(LRRC15) 유도.
- Sci Rep 2023 (s41598-023-48073-w) — 커큐민이 CAF α-SMA/COX-2↓ 재프로그래밍.
- ScienceDirect S2949866X23000126 — 진세노사이드 Rg3 항암 기전(PD-L1↓, Akt/mTOR↓).
"""

from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree

CELL_TYPES = ["Tumor", "myCAF", "iCAF", "CD8_T", "Macrophage"]

# ==========================================================================
# 기저 파라미터 (문헌 접지 기본값 — 모두 튜닝 가능)
# ==========================================================================
DEFAULT_PARAMS = dict(
    field_um=1500.0,
    dt_days=0.5,               # 시간 스텝
    # --- 종양 ---
    k_prolif=0.11,             # 분열 확률/day (배가 ~6day; 면역 살상이 경쟁 가능하도록)
    tumor_density_cap=10,      # kill_radius 내 최대 종양 수 (접촉 억제)
    k_tumor_apoptosis=0.01,    # 기저 자연 사멸/day
    # --- myCAF 장벽 ---
    k_caf_activate=0.18,       # 종양 주변 기질->myCAF 활성화 확률/day
    caf_ring_um=150.0,         # 활성화가 일어나는 종양으로부터의 거리대
    caf_cap_per_tumor=0.9,     # 종양당 myCAF 상한 (포화)
    # myCAF의 '포섭된 촉진자' 역할: 국소 myCAF가 종양 증식을 얼마나 높이나.
    # 0=순수 장벽(억제)만, >0=전암 토양·증식촉진(포섭). PDAC 기본 0, HCC는 크게.
    caf_protumor=0.0,          # 증식 부스트 강도 (0~1+)
    caf_boost_ref=18,          # 완전 부스트 기준 국소 myCAF 수 (높을수록 선형·비포화)
    k_caf_death=0.07,          # myCAF 비활성화(→iCAF) turnover/day. 항섬유화가
                               # k_caf_activate를 낮추면 이 turnover로 장벽 질량이 실제 감소
    # --- myCAF 물리적 통제 장벽 (L1: containment 구현) ---
    # 논문 중심가설("제거가 아닌 장벽 활용")을 모델에 직접 구현: 국소 myCAF/ECM 밀도가
    # 종양의 공간확장을 물리적으로 가둔다(딸세포가 기질벽 안쪽으로 못 나감). 0이면 구모델
    # (myCAF가 종양에 물리적 영향 없음 = 리뷰어가 지적한 상태).
    caf_confine=0.8,           # 종양 확장(딸세포 배치) 차단 강도 (0=없음, ~1=강한 가둠)
    caf_confine_ref=3,         # 완전 confinement 기준 국소(20µm) myCAF 수 (실측 밀도 접지)
    caf_pressure=1.2,          # 기계적 배제: 국소 myCAF가 종양 수용력을 지수적으로 낮춤
                               # (기질벽에선 종양이 못 쌓임 → 공간확산 자체를 제한). 0=없음
    # 트레이드오프: 조밀 기질은 약물 침투도 저해 → myCAF 과다는 오히려 통제 악화
    # (면역배제 cd8_barrier_alpha와 함께). 이 세 힘의 균형이 '최적 기질 상태'를 만든다.
    caf_drug_block=0.6,        # 조밀 myCAF의 약물 침투 저해 강도 (0=없음)
    # --- CD8 ---
    cd8_speed_um=120.0,        # 종양 방향 이동 속도 um/day (장벽 없을 때)
    cd8_barrier_alpha=0.9,     # 회랑 myCAF 밀도에 대한 이동 감쇠 계수 (핵심 게이트)
    corridor_um=30.0,          # 장벽 회랑 폭
    kill_radius_um=20.0,       # CD8가 종양을 죽일 수 있는 거리
    k_kill=1.0,                # 접촉 시 종양 사멸 확률/day (~1 kill/day)
    cd8_recruit=30.0,          # 주변부 신규 CD8 유입/day
    k_cd8_death=0.03,          # CD8 turnover/day
    cd8_entry_margin_um=60.0,  # 신규 CD8가 들어오는 조직 가장자리 폭
    # --- 저항성(내성) 아형: adaptive therapy의 핵심 (문헌 접지) ---
    # 젬시타빈 내성은 주로 비유전적·적응적(후성유전·미세환경; PMC5704175) → 소수 사전풀 +
    # 표현형 전환으로 근사. 점돌연변이율(~1e-6)은 ABM 규모·기간에 무의미하므로 미사용.
    init_resistant_frac=0.01,  # 초기 내성 분율 (표현형 내성 소수 사전풀)
    mutation_rate=0.001,       # 감수성->내성 표현형 전환(분열당). 내성은 약물 효과 회피
    resistant_immune_evasion=0.45,  # 내성 클론 면역회피(0~1): PD-L1↑·항원(B2M/MHC-I)
                               # 소실로 CD8 살상 회피 (JEM 2017 28302645; 완전소실은 ~1.0).
    resistance_cost=0.24,      # 내성 적합도 비용: 내성은 약물 없을 때 감수성의 ~76% 속도로
                               # 증식(NSCLC 실측 75.6%, Science Adv abm7212). 이 비용 덕에
                               # 휴약 중 감수성이 내성을 경쟁 억제 → adaptive therapy 성립.
                               # (단 회전율 높으면 비용 없이도 성립: Cancer Res 2021 33172930)
    seed=0,
)

# 교란 가능한 파라미터 (물질 배수가 곱해지는 대상)
PERTURBABLE = ["k_prolif", "k_caf_activate", "k_kill",
               "cd8_recruit", "cd8_speed_um", "k_tumor_apoptosis"]


# ==========================================================================
# 장기 맥락 (CONTEXTS) — 성상세포·TGF-β 축은 공유, 파라미터·의미만 변경.
# PDAC: 결합조직증식이 '장벽(억제)'으로 작동 (caf_protumor=0).
# HCC : 간경화(HSC 활성)가 '전암 토양(촉진)'으로 작동 (caf_protumor 큼) —
#       사용자 통찰 '포섭된 촉진자'를 반영. 문헌: 활성 HSC가 HCC 기질에 침윤해
#       성장인자·사이토카인으로 종양세포 증식을 자극 (PMC7419619).
# ==========================================================================
CONTEXTS = {
    "pdac": dict(
        label="췌장암 (PDAC)",
        params=dict(caf_protumor=0.0, k_prolif=0.11),
        display=dict(tumor="종양 (PDAC)", barrier="myCAF",
                     stellate="췌장성상세포 (PSC)"),
        note="결합조직증식(desmoplasia)이 장벽으로 작동 — 면역배제 우세. "
             "CAF 제거가 오히려 악화된 맥락(Özdemir/Rhim 2014).",
    ),
    "hcc": dict(
        label="간세포암 (HCC) / 간경화",
        params=dict(caf_protumor=0.9, k_prolif=0.11, k_caf_activate=0.22),
        display=dict(tumor="간세포암 (HCC)",
                     barrier="활성 성상세포 (myofibroblast)",
                     stellate="간성상세포 (HSC)"),
        note="간경화(HSC 활성)가 전암 토양 — 섬유화가 증식을 '촉진'(포섭된 촉진자 우세). "
             "HCC의 ~80–90%가 경화된 간에서 발생 (PMC7419619).",
    ),
}


def context_params(context="pdac"):
    """맥락 프리셋을 DEFAULT_PARAMS 위에 병합한 파라미터 dict 반환."""
    p = dict(DEFAULT_PARAMS)
    p.update(CONTEXTS.get(context, {}).get("params", {}))
    return p


# ==========================================================================
# 물질 라이브러리 (문헌 기반 기전 -> 파라미터 배수, dose=1.0 기준)
# 각 항목: 배수 dict + 근거 한 줄. dose로 강도를 선형 보간.
# ==========================================================================
# evidence: strong=다수 in vitro/in vivo, moderate=일부 연구, weak=예비/보조적,
#           reference=표준요법 벤치마크. 근거 수준은 정직하게 표기해 UI에 노출한다.
SUBSTANCES = {
    "curcumin": dict(
        label="커큐민 (Curcumin)", evidence="strong",
        effects={"k_caf_activate": 0.45, "k_kill": 1.20},
        rationale="CAF의 α-SMA·COX-2↓ 재프로그래밍, Treg↓·IFN-γ↑ "
                  "(Sci Rep 2023 s41598-023-48073-w) → 장벽 형성↓, 면역살상↑",
    ),
    "ginsenoside_rg3": dict(
        label="진세노사이드 Rg3", evidence="strong",
        effects={"k_kill": 1.60, "k_prolif": 0.75},
        rationale="PD-L1 글리코실화↓·Granzyme B/perforin↑ → 살상↑; "
                  "PI3K/Akt/mTOR↓ → 증식↓ (ScienceDirect S2949866X23000126)",
    ),
    "gemcitabine": dict(
        label="젬시타빈 (표준 화학요법, 대조 기준)", evidence="reference",
        effects={"k_prolif": 0.55, "k_tumor_apoptosis": 3.0},
        rationale="뉴클레오사이드 유사체 — 증식 억제·세포사멸 유도 (표준요법 벤치마크)",
    ),
    "erlotinib": dict(
        label="에를로티닙 (EGFR 표적 TKI)", evidence="strong",
        effects={"k_prolif": 0.65},
        rationale="EGFR 인산화효소 억제(quinazoline TKI) → Ras/MAPK·PI3K/Akt↓ 증식 억제. "
                  "젬시타빈이 유도하는 MAPK/ERBB2 생존신호를 차단 → 병용 시 3상 생존이득 "
                  "(Cancer Res 2013 73:2221; 췌장암 젬시타빈 병용 승인)",
    ),
    "entecavir": dict(
        label="엔테카비르 (B형간염 항바이러스, 리퍼포징 가설)", evidence="hypothesis",
        effects={"k_prolif": 0.96},   # 매우 미미 — KDM5B 가설뿐, PDAC 근거 없음
        rationale="구아노신 유사체 HBV 중합효소 억제제. 간섬유화 되돌림은 **바이러스 억제로 "
                  "인한 간접효과**(원인 제거)이며 myCAF 직접작용 아님. 항암은 KDM5B 억제 "
                  "in silico 가설(Pharmaceuticals 2023). **췌장암 적용은 미검증 가설**",
    ),
    "gv1001": dict(
        label="GV1001 (텔로머라제 hTERT 펩타이드 백신)", evidence="moderate",
        effects={"cd8_recruit": 1.6, "k_kill": 1.3},
        rationale="hTERT 유래 16-mer 펩타이드 백신 → 텔로머라제⁺ 종양 대상 CD4/CD8 T세포 "
                  "프라이밍(면역 부스트). 비선택 진행성 PDAC은 음성(TeloVac 3상, Lancet "
                  "Oncol 2014)이나 **eotaxin-high 하위군에서 생존 개선**(KG4/2015 3상, BJC "
                  "2023; OS 11.3 vs 7.5m). 면역 프라이밍 제제 — 효과는 T세포가 종양에 "
                  "도달할 수 있어야 발현(장벽 열림 필요)",
    ),
    # ---- 전통 약재 (문헌 기반, 근거 수준 상이 — 반드시 evidence 확인) ----
    "garlic": dict(
        label="마늘 (Allicin/DATS)", evidence="strong",
        effects={"k_prolif": 0.70, "k_tumor_apoptosis": 2.0},
        rationale="Allicin·DATS가 Bcl-2↓·caspase↑ 세포사멸, NF-κB/PI3K-Akt↓; "
                  "췌장암세포(Capan-2) 직접 억제 (Discover Oncol 2025)",
    ),
    "mugwort": dict(
        label="쑥 (Artesunate/DHA, Eupatilin)", evidence="strong",
        effects={"k_caf_activate": 0.40, "k_kill": 1.20, "k_prolif": 0.85},
        rationale="아르테수네이트/DHA가 TGF-β 억제로 CAF를 활성→휴지 비활성화, "
                  "면역억제(TGF-β1/IL-10)↓; eupatilin 항섬유화 (JECCR 2018)",
    ),
    "wild_ginseng": dict(
        label="산삼 (Ginsenosides Rb1/Rg1/Rg3)", evidence="strong",
        effects={"k_caf_activate": 0.50, "k_kill": 1.30, "k_prolif": 0.85},
        rationale="진세노사이드가 성상세포(HSC) TGF-β1·콜라겐↓ 항섬유화(Smad 회복), "
                  "MAPK/NF-κB↓ 세포사멸, Rg3 PD-L1↓ (BMC Complement Med 2016)",
    ),
    "platycodon": dict(
        label="도라지 (Platycodin D)", evidence="moderate",
        effects={"k_kill": 1.40, "cd8_recruit": 1.30, "k_prolif": 0.85},
        rationale="Platycodin D가 PD-1⁺CD8↓로 CD8 침윤·살상↑(VEGF-A/STAT3 경유), "
                  "PUMA↑ 세포사멸 (Front Pharmacol 2022)",
    ),
    "sea_cucumber": dict(
        label="해삼 (Frondoside A)", evidence="moderate",
        effects={"k_prolif": 0.75, "k_tumor_apoptosis": 1.8},
        rationale="Frondoside A가 p53↑·caspase↑ 세포사멸, S/G2M 정지, "
                  "EGFR/Akt/MMP-9↓ 전이 억제 (Mar Drugs 2015; PMC9861588)",
    ),
    "deer_antler": dict(
        label="녹용 (Velvet peptides/Pantocrine)", evidence="weak",
        effects={"cd8_recruit": 1.40, "k_caf_activate": 0.80},
        rationale="녹용 펩타이드가 NK활성·조혈 촉진(면역지원), TGF-β/NF-κB 항섬유화. "
                  "**직접 항암 근거는 약함 — 보조적** (PMC11125008, PMC8814364)",
    ),
    # ---- 간섬유화 항섬유화 약재 (HCC 맥락에서 특히 유효) ----
    "danshen": dict(
        label="단삼 (Salvianolic acid B / Tanshinone IIA)", evidence="strong",
        effects={"k_caf_activate": 0.40},
        rationale="Salvianolic acid B가 성상세포 활성을 TGF-β1/p38·ERK MAPK 억제로 차단, "
                  "Tanshinone IIA는 활성 HSC 세포사멸 유도 (PMC9385955)",
    ),
    "milk_thistle": dict(
        label="밀크씨슬 (Silymarin / Silibinin)", evidence="strong",
        effects={"k_caf_activate": 0.45, "k_prolif": 0.85},
        rationale="Silibinin이 HSC 활성·증식 억제(α-SMA↓, p27/p53↑·Akt↓), "
                  "HCC 세포 증식도 억제 (PMC5052367)",
    ),
    "astragaloside": dict(
        label="황기 (Astragaloside IV)", evidence="strong",
        effects={"k_caf_activate": 0.40},
        rationale="Astragaloside IV가 TGF-β1/Smad 경로 조절(TGF-β1·p-Smad2/3↓, Smad7↑)로 "
                  "성상세포 활성·콜라겐 형성 억제 (PMC5952439)",
    ),
    # ---- 일반 기전 클래스 ----
    "generic_antifibrotic": dict(
        label="일반 항섬유화 (예: TGF-β 억제 계열)", evidence="reference",
        effects={"k_caf_activate": 0.35},
        rationale="myCAF 활성화 억제 → 장벽 붕괴 (기전 클래스)",
    ),
    "generic_immunostim": dict(
        label="일반 면역증강", evidence="reference",
        effects={"cd8_recruit": 1.8, "k_kill": 1.4},
        rationale="CD8 유입·살상능 증강 (기전 클래스)",
    ),
    "generic_cytotoxic": dict(
        label="일반 세포독성", evidence="reference",
        effects={"k_prolif": 0.6, "k_tumor_apoptosis": 2.5},
        rationale="직접 증식 억제·사멸 유도 (기전 클래스)",
    ),
}

# 근거 수준 배지 (UI 표기용)
EVIDENCE_BADGE = {
    "strong":     "🟢 강함(다수 연구)",
    "moderate":   "🟡 보통(일부 연구)",
    "weak":       "🟠 약함(예비·보조)",
    "reference":  "⚪ 기준/벤치마크",
    "hypothesis": "🔴 가설(미검증)",
}

# ==========================================================================
# 독성(환자 부담) 가중치 0~1 — 화학요법·표적치료↑, 천연물↓.
# ⚠️ 예시 값일 뿐 임상 독성 수치가 아니다. 순차/사이클 투여가 누적 부담을 줄이는지
# 비교하기 위한 상대적 지표. 종양 억제만 보던 모델의 맹점(고통·독성)을 보완한다.
# ==========================================================================
TOXICITY = {
    "gemcitabine": 0.85, "erlotinib": 0.60, "entecavir": 0.25, "gv1001": 0.20,
    "generic_cytotoxic": 0.80, "generic_immunostim": 0.30,
    "generic_antifibrotic": 0.20,
    "curcumin": 0.10, "ginsenoside_rg3": 0.15, "garlic": 0.10, "mugwort": 0.20,
    "wild_ginseng": 0.15, "platycodon": 0.10, "sea_cucumber": 0.10,
    "deer_antler": 0.10, "danshen": 0.15, "milk_thistle": 0.10,
    "astragaloside": 0.10,
}


def regimen_toxicity(subs):
    """활성 물질 목록 [(name,dose)]의 순간 독성률 합."""
    return float(sum(TOXICITY.get(n, 0.3) * float(d) for n, d in subs))


def compose_regimen(selected, synergy=0.0):
    """선택된 물질들(name, dose)을 배수 dict로 합성.

    selected: [(name, dose_0_to_1), ...]
    각 물질의 배수를 dose로 1.0(무효과)에서 선형 보간한 뒤 곱한다.
    synergy>0이면 항암 방향(증식·CAF 배수를 추가로 낮춤) 시너지를 약하게 준다.
    반환: {param: multiplier}
    """
    mult = {p: 1.0 for p in PERTURBABLE}
    for name, dose in selected:
        if name not in SUBSTANCES or dose <= 0:
            continue
        for p, m in SUBSTANCES[name]["effects"].items():
            eff = 1.0 + (m - 1.0) * float(np.clip(dose, 0, 1))  # dose 보간
            mult[p] = mult.get(p, 1.0) * eff
    if synergy > 0 and len(selected) >= 2:
        for p in ("k_prolif", "k_caf_activate"):
            mult[p] *= (1.0 - 0.15 * synergy)   # 조합 시 억제 방향 약한 시너지
    return mult


# ==========================================================================
# ABM
# ==========================================================================
class TumorABM:
    def __init__(self, coords, labels, params=None, regimen=None,
                 regimen_subs=None, schedule=None, synergy=0.0,
                 adaptive=False, adapt_on=1.2, adapt_off=0.6):
        """regimen: 정적 배수 dict(레거시). regimen_subs: [(name,dose)] (독성 추적).
        schedule: 시차 순차 사이클 = [{'label','days','substances':[(name,dose)]}, ...]

        adaptive=True: 적응형 on/off 제어(Gatenby). 종양이 초기의 adapt_on 배를 넘으면
        투여, adapt_off 배 밑으로 떨어지면 휴약 → 감수성 세포를 남겨 내성을 억누르고
        (경쟁적 억제) 최소 부담으로 오래 통제(공존)한다.
        """
        self.coords = np.asarray(coords, float).copy()
        self.labels = np.asarray(labels, dtype=object).copy()
        self.p = dict(DEFAULT_PARAMS)
        if params:
            self.p.update(params)
        self.rng = np.random.default_rng(self.p["seed"])
        self.schedule = schedule
        self.regimen_subs = regimen_subs
        self.static_regimen = regimen or {}
        self.synergy = synergy
        self.adaptive = adaptive
        self.adapt_on = adapt_on
        self.adapt_off = adapt_off
        self.drug_on = True
        # 저항성(내성) 라벨: 종양 세포별 내성 여부(coords와 정렬). 비종양=False.
        self.res = np.zeros(len(self.coords), dtype=bool)
        tmask = self._mask("Tumor")
        self.res[tmask] = self.rng.random(int(tmask.sum())) < self.p["init_resistant_frac"]
        self.n0 = int(tmask.sum())
        self.crit_time = None       # 초기의 adapt_on 배 도달 시각(time-to-burden)
        self.t = 0.0
        self.cum_tox = 0.0
        self.history = []
        self._apply_regimen()   # self.eff, self._active_subs, self.phase_label 설정
        self._record()

    def _apply_regimen(self):
        """현재 시점 self.t에 맞는 배수를 self.eff에 반영. 순차 스케줄이면 phase 선택.
        adaptive면 종양 부담 band로 투여 on/off를 먼저 결정(off면 무투여)."""
        subs, label = None, None
        if self.adaptive:
            n_now = int(self._mask("Tumor").sum())
            if n_now >= self.adapt_on * self.n0:
                self.drug_on = True
            elif n_now <= self.adapt_off * self.n0:
                self.drug_on = False
            # band 사이면 직전 상태 유지(hysteresis)
            if not self.drug_on:
                self.eff = dict(self.p); self._active_subs = None
                self.phase_label = "휴약"
                return
        if self.schedule:
            cyc = sum(ph["days"] for ph in self.schedule) or 1
            tt = self.t % cyc
            acc, ph = 0.0, self.schedule[-1]
            for cand in self.schedule:
                acc += cand["days"]
                if tt < acc:
                    ph = cand
                    break
            subs, label = ph["substances"], ph.get("label")
            reg = compose_regimen(subs, synergy=self.synergy)
        elif self.regimen_subs is not None:
            subs = self.regimen_subs
            reg = compose_regimen(subs, synergy=self.synergy)
        else:
            reg = self.static_regimen
        self.eff = dict(self.p)
        for k, m in reg.items():
            if k in self.eff:
                self.eff[k] = self.p[k] * m
        self._active_subs = subs
        self.phase_label = label

    # --- helpers ---
    def _mask(self, ct):
        return self.labels == ct

    def _corridor_density(self, cd8_pts, tumor_pts, caf_tree):
        """각 CD8->최근접종양 경로의 myCAF 회랑 밀도 (장벽점수와 동일 논리)."""
        if len(cd8_pts) == 0 or caf_tree is None:
            return np.zeros(len(cd8_pts))
        cu = self.p["corridor_um"]
        out = np.zeros(len(cd8_pts))
        for i, (p0, p1) in enumerate(zip(cd8_pts, tumor_pts)):
            seg = p1 - p0
            L = np.linalg.norm(seg)
            if L < 1e-9:
                continue
            n_steps = max(2, int(L / (cu / 2)))
            ts = np.linspace(0.15, 0.85, n_steps)
            pts = p0 + ts[:, None] * seg
            n_near = caf_tree.query_ball_point(pts, r=cu, return_length=True)
            out[i] = float(np.mean(n_near))
        return out

    # --- one time step ---
    def step(self):
        dt = self.p["dt_days"]
        field = self.p["field_um"]
        # 현재 시점의 phase에 맞춰 약물 배수 갱신 + 누적 독성 적산
        self._apply_regimen()
        if self._active_subs:
            self.cum_tox += regimen_toxicity(self._active_subs) * dt
        tmask0 = self._mask("Tumor")
        tum = self.coords[tmask0]
        tum_res = self.res[tmask0]
        if len(tum) == 0:
            self.t += dt
            self._record()
            return
        tumor_tree = cKDTree(tum)
        # myCAF 트리를 증식 전에 계산 (포섭된 촉진자 부스트에 사용)
        caf = self.coords[self._mask("myCAF")]
        caf_tree = cKDTree(caf) if len(caf) else None

        # (1) 종양 증식 — 국소 밀도 수용력 + myCAF 물리적 통제(confinement) + 약물침투 저해
        # 내성: 약물 회피(기저 k_prolif·적합도비용). 감수성: 약물 적용하되 조밀 myCAF가 침투 저해.
        new_coords, new_labels, new_res = [], [], []
        kr = self.p["kill_radius_um"]
        cap = self.p["tumor_density_cap"]
        protumor = self.p.get("caf_protumor", 0.0)
        boost_ref = max(1, self.p.get("caf_boost_ref", 8))
        mut = self.p.get("mutation_rate", 0.0)
        base_kp = self.p["k_prolif"]
        drug_reduction = max(0.0, base_kp - self.eff["k_prolif"])  # 약물의 감수성 증식억제분
        kp_res = base_kp * (1.0 - self.p.get("resistance_cost", 0.0))
        confine = self.p.get("caf_confine", 0.0)
        confine_ref = max(1, self.p.get("caf_confine_ref", 12))
        drug_block = self.p.get("caf_drug_block", 0.0)
        pressure = self.p.get("caf_pressure", 0.0)
        use_caf = caf_tree is not None and (protumor > 0 or confine > 0
                                            or drug_block > 0 or pressure > 0)
        for c, is_res in zip(tum, tum_res):
            local = tumor_tree.query_ball_point(c, r=kr, return_length=True)
            caf_here = (caf_tree.query_ball_point(c, r=kr, return_length=True)
                        if use_caf else 0)          # 국소(20µm) myCAF: 약물차단·가둠·압력용
            boost, eff_cap = 1.0, float(cap)
            if protumor > 0 and caf_tree is not None:
                caf2 = caf_tree.query_ball_point(c, r=kr * 2, return_length=True)
                soil = protumor * min(caf2 / boost_ref, 1.0)
                boost = 1.0 + soil
                eff_cap = cap * (1.0 + soil)
            # 기계적 배제: 국소 myCAF가 종양 수용력을 지수적으로 낮춤 → 기질벽에선 못 쌓임
            if pressure > 0 and caf_here:
                eff_cap *= np.exp(-pressure * caf_here / confine_ref)
            if local >= eff_cap:                 # 국소 수용력 → 감수성·내성 경쟁
                continue
            if is_res:
                kp = kp_res
            else:                                # 조밀 myCAF → 약물 침투 저해 → 약효 감소
                pen = np.exp(-drug_block * caf_here / confine_ref) if drug_block > 0 else 1.0
                kp = base_kp - drug_reduction * pen
            if self.rng.random() < kp * boost * dt:
                ang = self.rng.uniform(0, 2 * np.pi)
                off = self.rng.uniform(4, kr) * np.array([np.cos(ang), np.sin(ang)])
                dpos = c + off
                # myCAF 물리적 통제: 딸세포가 기질-밀집 위치로 확장하려 하면 차단(가둠)
                if confine > 0 and caf_tree is not None:
                    caf_d = caf_tree.query_ball_point(dpos, r=kr, return_length=True)
                    if self.rng.random() < confine * min(caf_d / confine_ref, 1.0):
                        continue                 # 기질벽에 막힘 → 확장 못함 (containment)
                new_coords.append(dpos)
                new_labels.append("Tumor")
                daughter_res = bool(is_res) or (self.rng.random() < mut)  # 유전+변이
                new_res.append(daughter_res)

        # (2) myCAF 장벽 활성화: 종양 링 거리대의 기질세포 일부를 myCAF로 전환/증식
        ring = self.p["caf_ring_um"]
        stromal_mask = np.isin(self.labels, ["iCAF", "Macrophage"])
        my_now = int(self._mask("myCAF").sum())
        my_cap = int(self.p["caf_cap_per_tumor"] * len(tum))
        if my_now < my_cap:
            idxs = np.where(stromal_mask)[0]
            if len(idxs):
                d, _ = tumor_tree.query(self.coords[idxs], k=1)
                near = idxs[(d > 40) & (d < ring)]
                prob = self.eff["k_caf_activate"] * dt
                conv = near[self.rng.random(len(near)) < prob]
                # 상한까지만
                conv = conv[: max(0, my_cap - my_now)]
                self.labels[conv] = "myCAF"

        # (2b) myCAF 비활성화 turnover (→iCAF). 활성화(k_caf_activate)와 균형을 이뤄
        # 장벽 질량의 평형을 만든다. 항섬유화제로 활성화가 낮아지면 평형 myCAF가 감소.
        my_idx = np.where(self._mask("myCAF"))[0]
        if len(my_idx):
            dead = my_idx[self.rng.random(len(my_idx)) < self.p["k_caf_death"] * dt]
            self.labels[dead] = "iCAF"

        # (3)+(4) CD8 이동(장벽 감쇠) 및 살상
        cd8_idx = np.where(self._mask("CD8_T"))[0]
        caf = self.coords[self._mask("myCAF")]
        caf_tree = cKDTree(caf) if len(caf) else None
        killed = set()
        if len(cd8_idx):
            cd8_pts = self.coords[cd8_idx]
            d_ct, ti = tumor_tree.query(cd8_pts, k=1)
            partners = tum[ti]
            corr = self._corridor_density(cd8_pts, partners, caf_tree)
            mob = np.exp(-self.eff.get("cd8_barrier_alpha",
                         self.p["cd8_barrier_alpha"]) * corr)  # 0~1 감쇠
            # 이동: 종양 방향 단위벡터 * 속도 * 이동성
            vecs = partners - cd8_pts
            norms = np.linalg.norm(vecs, axis=1, keepdims=True)
            norms[norms < 1e-9] = 1e-9
            step_len = (self.eff["cd8_speed_um"] * dt) * mob[:, None]
            self.coords[cd8_idx] = cd8_pts + (vecs / norms) * np.minimum(step_len, norms)
            # 살상: 이동 후 kill_radius 안의 종양. 내성 클론은 면역회피로 덜 죽음.
            d_after, ti2 = tumor_tree.query(self.coords[cd8_idx], k=1)
            engaged = np.where(d_after <= kr)[0]
            evade = self.p.get("resistant_immune_evasion", 0.0)
            for e in engaged:
                ki = int(ti2[e])
                kk = self.eff["k_kill"] * (1.0 - evade) if tum_res[ki] else self.eff["k_kill"]
                if self.rng.random() < kk * dt:
                    killed.add(ki)

        # (5) 종양 자연 사멸 — 약물 유도 사멸(k_tumor_apoptosis 증가분)은 감수성만.
        tum_all_idx = np.where(self._mask("Tumor"))[0]
        res_all = self.res[tum_all_idx]
        apop_rate = np.where(res_all, self.p["k_tumor_apoptosis"],
                             self.eff["k_tumor_apoptosis"]) * dt
        apop = self.rng.random(len(tum_all_idx)) < apop_rate
        # killed는 tum(로컬 좌표) 인덱스 -> 전역 인덱스로 매핑 (면역살상은 내성 무관)
        tumor_global = tum_all_idx  # 동일 순서 (mask 순서 보존)
        remove = set(int(tumor_global[k]) for k in killed if k < len(tumor_global))
        remove |= set(int(i) for i in tum_all_idx[apop])

        # (6) CD8 turnover + 신규 유입
        cd8_now = np.where(self._mask("CD8_T"))[0]
        die = cd8_now[self.rng.random(len(cd8_now)) < self.p["k_cd8_death"] * dt]
        remove |= set(int(i) for i in die)

        # 적용: 제거 (res 배열도 동일하게)
        if remove:
            keep = np.array([i for i in range(len(self.labels)) if i not in remove])
            self.coords = self.coords[keep]
            self.labels = self.labels[keep]
            self.res = self.res[keep]

        # 신규 종양(증식) 추가 — 딸세포 내성 라벨 포함
        if new_coords:
            self.coords = np.vstack([self.coords, np.array(new_coords)])
            self.labels = np.concatenate([self.labels, np.array(new_labels, dtype=object)])
            self.res = np.concatenate([self.res, np.array(new_res, dtype=bool)])

        # 신규 CD8 유입 (가장자리) — 비종양이므로 res=False
        n_new = self.rng.poisson(self.eff["cd8_recruit"] * dt)
        if n_new > 0:
            m = self.p["cd8_entry_margin_um"]
            pts = self.rng.uniform(0, field, (n_new, 2))
            edge = (pts[:, 0] < m) | (pts[:, 0] > field - m) | \
                   (pts[:, 1] < m) | (pts[:, 1] > field - m)
            pts = pts[edge] if edge.any() else pts
            if len(pts):
                self.coords = np.vstack([self.coords, pts])
                self.labels = np.concatenate(
                    [self.labels, np.array(["CD8_T"] * len(pts), dtype=object)])
                self.res = np.concatenate([self.res, np.zeros(len(pts), dtype=bool)])

        # 경계 클리핑
        self.coords = np.clip(self.coords, 0, field)
        self.t += dt
        self._record()

    def _record(self):
        lab = self.labels
        tmask = lab == "Tumor"
        tum = self.coords[tmask]
        cd8 = self.coords[lab == "CD8_T"]
        n_tumor = int(tmask.sum())
        n_res = int(self.res[tmask].sum())
        # CD8 침윤(진짜 관통): 종양 kill 사거리(25um) 내 CD8 비율.
        infil = 0.0
        if len(tum) and len(cd8):
            d, _ = cKDTree(tum).query(cd8, k=1)
            infil = float(np.mean(d <= 25))
        # 임계부담(초기 adapt_on 배) 첫 도달 시각 = time-to-progression 대리지표
        if self.crit_time is None and self.n0 > 0 and n_tumor >= self.adapt_on * self.n0:
            self.crit_time = round(self.t, 2)
        self.history.append(dict(
            t=round(self.t, 2),
            n_tumor=n_tumor,
            n_tumor_resistant=n_res,
            n_tumor_sensitive=n_tumor - n_res,
            resistant_frac=round(n_res / n_tumor, 3) if n_tumor else 0.0,
            n_myCAF=int((lab == "myCAF").sum()),
            n_CD8=int((lab == "CD8_T").sum()),
            cd8_infiltration=infil,
            cum_toxicity=round(self.cum_tox, 3),
            drug_on=bool(self.drug_on),
            phase=self.phase_label,
        ))

    def run(self, days=20.0, snapshots=(0.0, 0.5, 1.0)):
        """days까지 진행. snapshots(진행률 0~1 목록)에서 조직 스냅샷 저장."""
        n_steps = int(round(days / self.p["dt_days"]))
        snap_steps = {int(round(s * n_steps)): s for s in snapshots}
        snaps = {}
        if 0 in snap_steps:
            snaps[0.0] = (self.coords.copy(), self.labels.copy())
        for i in range(1, n_steps + 1):
            self.step()
            if i in snap_steps:
                snaps[round(self.t, 2)] = (self.coords.copy(), self.labels.copy())
        return self.history, snaps


def simulate(coords, labels, days=20.0, params=None, regimen=None,
             regimen_subs=None, schedule=None, synergy=0.0,
             adaptive=False, adapt_on=1.2, adapt_off=0.6,
             snapshots=(0.0, 0.5, 1.0)):
    """편의 함수: ABM 한 번 실행하고 (history, snapshots) 반환.

    schedule 지정 시 시차 순차 사이클, regimen_subs 지정 시 정적 투여(독성 추적 포함).
    adaptive=True면 종양 부담 band로 투여 on/off (Gatenby 적응요법).
    """
    m = TumorABM(coords, labels, params=params, regimen=regimen,
                 regimen_subs=regimen_subs, schedule=schedule, synergy=synergy,
                 adaptive=adaptive, adapt_on=adapt_on, adapt_off=adapt_off)
    return m.run(days=days, snapshots=snapshots)


def control_metrics(history, n0=None, crit_mult=1.5, dt=0.5):
    """통제(공존) 목적함수 지표. '박멸' 대신 '최소 부담으로 오래 통제'를 평가.

    - ttp: time-to-progression 대리 = 종양이 초기의 crit_mult 배에 처음 도달한 시각
           (도달 안 하면 관측기간 = censored, 더 좋음).
    - peak_frac: 관측기간 중 최대 종양/초기 비.
    - final_resistant_frac: 종료 시 내성 분율(경쟁적 방출 지표, 낮을수록 통제 양호).
    - cum_toxicity: 누적 부담.
    - auc_burden: 종양부담 곡선하 면적/초기(낮을수록 부담 적게 통제).
    - control_score: 통제기간 / (독성+1) — 높을수록 '적은 부담으로 오래 통제'.
    """
    import numpy as np
    n0 = n0 or history[0]["n_tumor"]
    ts = np.array([h["t"] for h in history])
    nt = np.array([h["n_tumor"] for h in history], float)
    ttp = None
    for h in history:
        if h["n_tumor"] >= crit_mult * n0:
            ttp = h["t"]; break
    censored = ttp is None
    ttp_val = ts[-1] if censored else ttp
    _trap = getattr(np, "trapezoid", getattr(np, "trapz", None))
    auc = float(_trap(nt / n0, ts))
    return {
        "ttp_days": float(ttp_val), "progression_censored": censored,
        "peak_frac": float(nt.max() / n0),
        "final_frac": float(nt[-1] / n0),
        "final_resistant_frac": float(history[-1]["resistant_frac"]),
        "cum_toxicity": float(history[-1]["cum_toxicity"]),
        "auc_burden": auc,
        "control_score": float(ttp_val / (history[-1]["cum_toxicity"] + 1.0)),
    }


def build_cycle_schedule(antifibrotic_subs, cytotoxic_subs, phase_days=5):
    """사용자 예시의 3-phase 사이클 구성:
    ① myCAF 약화(항섬유화/면역) → ② 증식 억제(세포독성/표적) → ③ 기질 리셋(항섬유화).
    반복 사이클로 돌아간다."""
    phases = []
    if antifibrotic_subs:
        phases.append(dict(label="① myCAF 약화", days=phase_days,
                           substances=antifibrotic_subs))
    if cytotoxic_subs:
        phases.append(dict(label="② 증식 억제", days=phase_days,
                           substances=cytotoxic_subs))
    if antifibrotic_subs:
        phases.append(dict(label="③ 기질 리셋", days=phase_days,
                           substances=antifibrotic_subs))
    return phases or None


def classify_substances(subs):
    """[(name,dose)]를 ①③ 약화(항섬유화/면역)와 ② 증식억제(세포독성/표적)로 분류.
    면역·항섬유화 효과가 있으면 '약화'로, 순수 증식차단이면 '증식억제'로 보낸다.
    반환: (antifibrotic_immune, cytotoxic_targeted)."""
    anti, cyto = [], []
    for name, dose in subs:
        eff = SUBSTANCES.get(name, {}).get("effects", {})
        immune_or_antifib = (eff.get("k_caf_activate", 1.0) < 1.0 or
                             eff.get("k_kill", 1.0) > 1.0 or
                             eff.get("cd8_recruit", 1.0) > 1.0)
        pure_cyto = ((eff.get("k_prolif", 1.0) < 0.9 or
                      eff.get("k_tumor_apoptosis", 1.0) > 1.2)
                     and not immune_or_antifib)
        if immune_or_antifib:
            anti.append((name, dose))
        else:                       # 순수 세포독성/표적
            cyto.append((name, dose))
    return anti, cyto
