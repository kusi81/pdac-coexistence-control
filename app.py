"""
PDAC 공간 분석 대시보드 — CAF 장벽 / 면역배제 (containment) 검정

Layer 4 spatial pipeline을 브라우저 대시보드로 감싼 것.
데이터 소스(합성 또는 .h5ad 업로드) 선택 → 4개 지표 실행 → 결과·그림·해석 표시.

핵심 설계:
  - 가장 중요한 지표(②방사형 ③장벽점수 근접성)는 numpy/scipy만으로 항상 동작.
  - ①이웃농축 ④동시출현은 squidpy가 있을 때만 추가로 표시(우아한 성능저하).
"""

import io
import json
import sys
import os

# make the pipeline package importable
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline"))

import numpy as np
import pandas as pd
import streamlit as st

from synthetic import make_tissue
from analysis import (
    run_core_metrics, run_squidpy_metrics, squidpy_available,
    platform_ceiling, barrier_verdict, run_rim_panel,
)
from spatial_core import PLATFORM_SPECS
import plots

st.set_page_config(page_title="PDAC 공간 분석 대시보드",
                   page_icon="🔬", layout="wide")

# ── 모드 선택: 기존 공간 분석 vs 신규 교란 시뮬레이션(ABM) ─────────────────
st.sidebar.title("모드")
MODE = st.sidebar.radio(
    "mode", ["🔬 공간 분석", "💊 교란 시뮬레이션", "🎯 통제 최적화", "🧬 분자 결합 뷰어"],
    label_visibility="collapsed")
if MODE == "💊 교란 시뮬레이션":
    import sim_page
    sim_page.render()
    st.stop()
if MODE == "🎯 통제 최적화":
    import opt_page
    opt_page.render()
    st.stop()
if MODE == "🧬 분자 결합 뷰어":
    import mol_page
    mol_page.render()
    st.stop()

LEVEL_BADGE = {
    "weak":      ("🔴 약함(Weak)", "#C0392B"),
    "moderate":  ("🟠 보통(Moderate)", "#E67E22"),
    "supported": ("🟢 지지됨(Supported)", "#27AE60"),
}
VERDICT_BADGE = {
    "interposed":   ("✅ 장벽 개재 (containment 지지)", "#27AE60"),
    "none":         ("➖ 개재 근거 없음", "#7F8C8D"),
    "anti":         ("🔵 경로에서 배제됨(장벽 아님)", "#2E86AB"),
    "insufficient": ("⚠️ 세포 수 부족", "#95A5A6"),
}


# --------------------------------------------------------------------------
# Sidebar — configuration
# --------------------------------------------------------------------------
st.sidebar.title("🔬 설정")

st.sidebar.markdown("**데이터 소스**")
source = st.sidebar.radio(
    "source", ["합성 조직 생성", "실데이터: PDAC Xenium (GSE274673)",
               "실데이터: MIBI-TOF (squidpy)", ".h5ad 업로드"],
    label_visibility="collapsed")

uploaded = None
syn_mode = "contained"
seed = 42
mibitof_lib = None
xenium_bundle = None
# 소스별 세포 유형 매핑 기본값
def_tumor, def_immune, def_barriers = "Tumor", "CD8_T", "myCAF,iCAF"

if source == "합성 조직 생성":
    syn_mode = st.sidebar.selectbox(
        "합성 조직 유형", ["contained", "diffuse"],
        format_func=lambda m: {"contained": "contained — myCAF 링 + CD8 배제",
                               "diffuse": "diffuse — 동일 세포수, 구조 없음"}[m])
    seed = st.sidebar.number_input("seed", value=42, step=1)
    st.sidebar.caption("⚠️ 합성 조직은 지표 검증용입니다 — 수치를 실제 소견으로 보고하지 마세요.")
elif source.startswith("실데이터: PDAC Xenium"):
    if not squidpy_available():
        st.sidebar.error("scanpy/anndata가 필요합니다.")
    import data_loader as _dl
    _samples = _dl.list_xenium_samples()
    if not _samples:
        st.sidebar.warning("받아진 PDAC Xenium 시료가 아직 없습니다. 다운로드/추출 완료 후 "
                           "나타납니다. (진행 중이면 잠시 후 재시도)")
    else:
        _labels = [s[1] for s in _samples]
        _pick = st.sidebar.selectbox(
            f"시료 ({len(_samples)}개 사용가능)", _labels,
            help="GSE274673 췌장암 Xenium 단일세포. 치료 전(naive)/CRT 후 표기.")
        xenium_bundle = _samples[_labels.index(_pick)][0]
    def_tumor, def_immune, def_barriers = "Tumor", "CD8_T", "myCAF,iCAF"
    st.sidebar.caption("실 췌장암 단일세포(µm 좌표). 세포타입은 저자 module-score 마커로 "
                       "자동 주석(annotate_pdac). 플랫폼=xenium 권장.")
elif source.startswith("실데이터: MIBI-TOF"):
    if not squidpy_available():
        st.sidebar.error("squidpy/anndata가 필요합니다.")
    mibitof_lib = st.sidebar.selectbox(
        "시료 (library)", ["point8", "point16", "point23"],
        help="MIBI-TOF 대장암 3개 시료 중 선택 (서로 다른 좌표계라 분리 분석)")
    def_tumor, def_immune, def_barriers = "Epithelial", "Tcell_CD8", "Fibroblast"
    st.sidebar.caption("⚠️ 좌표는 픽셀 단위(µm 아님). 대장암(비-PDAC)이며 Fibroblast는 "
                       "myCAF/iCAF로 세분되지 않은 단일 클러스터입니다.")
else:
    uploaded = st.sidebar.file_uploader(
        "adata.obsm['spatial'](µm) + adata.obs['cell_type'] 필요", type=["h5ad"])
    if not squidpy_available():
        st.sidebar.warning("`.h5ad` 읽기에는 anndata가 필요합니다. 미설치 시 합성 조직을 사용하세요.")

st.sidebar.markdown("---")
platform = st.sidebar.selectbox(
    "플랫폼", list(PLATFORM_SPECS.keys()), index=2,
    help="플랫폼 해상도가 주장 가능한 최대 수위를 결정합니다.")

st.sidebar.markdown("**세포 유형 매핑**")
tumor = st.sidebar.text_input("종양(tumor)", def_tumor)
immune = st.sidebar.text_input("면역(immune)", def_immune)
barriers_str = st.sidebar.text_input("장벽 후보(barriers, 쉼표)", def_barriers)
barriers = tuple(b.strip() for b in barriers_str.split(",") if b.strip())

st.sidebar.markdown("**파라미터**")
radius = st.sidebar.slider("이웃 그래프 반경 (µm)", 20, 150, 50, 5)
corridor_um = st.sidebar.slider("장벽 회랑 폭 (µm)", 10, 60, 30, 5)
n_perms = st.sidebar.select_slider("순열 횟수", [200, 500, 1000, 2000], value=1000)

run_squidpy = st.sidebar.checkbox(
    "squidpy 지표(①이웃농축 ④동시출현) 포함", value=squidpy_available(),
    disabled=not squidpy_available(),
    help="squidpy 미설치 시 비활성화됩니다." if not squidpy_available() else None)

go = st.sidebar.button("▶ 분석 실행", type="primary", use_container_width=True)


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def _load_mibitof_cached(library):
    import data_loader as dl
    adata = dl.load_mibitof()
    ck = dl.find_celltype_key(adata)
    coords, labels = dl.to_coords_labels(adata, celltype_key=ck, library=library)
    return coords, labels


@st.cache_data(show_spinner="PDAC Xenium 로드·주석 중…")
def _load_xenium_cached(bundle):
    import data_loader as dl
    adata = dl.load_xenium_bundle(bundle)
    labels, cov = dl.annotate_pdac(adata)   # 저자 module-score 주석
    return np.asarray(adata.obsm["spatial"], float), labels, cov


def load_data():
    if source == "합성 조직 생성":
        coords, labels, _ = make_tissue(mode=syn_mode, seed=int(seed))
        return coords, labels, f"합성 조직 ({syn_mode}, seed={seed})"
    if source.startswith("실데이터: PDAC Xenium"):
        if xenium_bundle is None:
            return None, None, None
        coords, labels, cov = _load_xenium_cached(xenium_bundle)
        import os as _os
        tag = _os.path.basename(xenium_bundle)
        return coords, labels, f"PDAC Xenium (실 췌장암 단일세포, µm · 마커 {cov:.0%})"
    if source.startswith("실데이터: MIBI-TOF"):
        coords, labels = _load_mibitof_cached(mibitof_lib)
        return coords, labels, f"MIBI-TOF 대장암 · {mibitof_lib} (실데이터, 픽셀좌표)"
    if uploaded is None:
        return None, None, None
    import anndata as ad
    tmp = os.path.join(os.path.dirname(__file__), "data", "_uploaded.h5ad")
    with open(tmp, "wb") as f:
        f.write(uploaded.getbuffer())
    adata = ad.read_h5ad(tmp)
    if "spatial" not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] 좌표가 없습니다.")
    if "cell_type" not in adata.obs:
        raise ValueError("adata.obs['cell_type'] 라벨이 없습니다.")
    return (np.asarray(adata.obsm["spatial"]),
            np.asarray(adata.obs["cell_type"]), f"업로드: {uploaded.name}")


# --------------------------------------------------------------------------
# Header
# --------------------------------------------------------------------------
st.title("PDAC 공간 분석 대시보드")
st.markdown(
    "췌장암 조직에서 **CAF(암연관섬유아세포)가 종양과 T세포 사이를 물리적으로 "
    "가로막는지**(containment / 면역 배제)를 검정합니다. "
    "구성비(abundance)가 아니라 **위치(positioning)** 를 보는 것이 핵심입니다.")

if not go:
    lvl, cav = platform_ceiling(platform)
    badge, color = LEVEL_BADGE[lvl]
    st.info(
        f"왼쪽에서 데이터와 파라미터를 정한 뒤 **▶ 분석 실행**을 누르세요.\n\n"
        f"현재 플랫폼 **{platform}** → containment 주장 수위: {badge}")
    with st.expander("이 대시보드가 답하는 질문과 4개 지표", expanded=True):
        st.markdown("""
| # | 지표 | 질문 | 방향성 | 의존성 |
|---|------|------|:---:|---|
| 1 | 이웃 농축 | 누가 누구 옆에 있나 | ✗ | squidpy |
| 2 | 방사형 프로파일 + 근접성 | 거리에 따라 밀도가 어떻게 변하나 | 부분 | 기본 |
| 3 | **장벽 점수** | **CAF가 종양↔T세포 사이에 있나** | **✓** | 기본 |
| 4 | 동시출현 vs 반경 | 어떤 길이 척도에서 연관되나 | ✗ | squidpy |

**오직 지표 3(장벽 점수)만 containment를 검정합니다.** 나머지는 "서로 가깝다"까지만
말할 수 있고 *개재(interposition)* 개념이 없습니다. 귀무가설이 장벽 세포 **수를 맞추기**
때문에, 이 점수는 *양(abundance)이 아니라 위치*를 측정합니다.
""")
    st.stop()


# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
try:
    coords, labels, src_desc = load_data()
except Exception as e:  # noqa: BLE001
    st.error(f"데이터 로딩 실패: {e}")
    st.stop()

if coords is None:
    if source.startswith("실데이터: PDAC Xenium"):
        st.warning("사용 가능한 PDAC Xenium 시료를 선택하세요 (아직 없으면 다운로드 완료 대기).")
    else:
        st.warning("`.h5ad` 파일을 업로드하세요.")
    st.stop()

with st.spinner("분석 실행 중… (근접성/장벽 점수 순열 검정 포함)"):
    core = run_core_metrics(coords, labels, tumor=tumor, immune=immune,
                            barriers=barriers, corridor_um=corridor_um,
                            n_perms=n_perms, seed=0)
    sq = {}
    if run_squidpy:
        sq = run_squidpy_metrics(coords, labels, platform=platform,
                                 radius=radius, n_perms=n_perms, seed=0)

lvl, cav = platform_ceiling(platform)
badge, color = LEVEL_BADGE[lvl]

st.markdown(f"**데이터:** {src_desc}  |  **플랫폼:** `{platform}`")
if lvl == "weak":
    st.error(f"플랫폼 해상도 상한 → containment 주장: {badge}\n\n{cav}")
elif lvl == "moderate":
    st.warning(f"플랫폼 해상도 상한 → containment 주장: {badge}\n\n{cav}")
else:
    st.success(f"플랫폼 해상도 상한 → containment 주장: {badge}\n\n{cav}")

# --- 2D 단면 한계 배너 (실데이터) -----------------------------------------
if source != "합성 조직 생성":
    st.warning(
        "📐 **2D 단면 한계** — 조직은 3D이나 이 데이터는 **얇은 2D 단면 하나**입니다"
        "(Xenium은 시료당 1절편, z 정보 없음). 따라서: ① **거리 과대추정**(진짜 최근접 "
        "세포가 평면 위/아래일 수 있음) ② **얇은 장벽은 자른 위치에 의존**(myCAF 껍질이 "
        "링/덮개로 다르게 보임) ③ 종양은 3D 관(duct) 구조의 단면.")
    with st.expander("그래도 유효한 이유 & 무엇은 못 보는지"):
        st.markdown(
            "- ✅ **상대·통계 신호는 유효**: 장벽점수·rim 농축은 *같은 2D 단면 안의 무작위 "
            "재라벨 귀무*와 비교(matched-null)하므로, 2D 왜곡이 관측·귀무에 동일 적용 → "
            "**'무작위 기질보다 몰렸나'는 성립**. 대형 단면(수만 세포)은 3D의 좋은 평면 표본.\n"
            "- ❌ **절대 3D 기하는 불가**: 진짜 containment 껍질·절대 거리·경로 연결성은 "
            "한 단면으로 복원 못함 → **'myCAF가 종양을 3D로 완전히 감쌌다'류 절대 주장은 과신 금물**.\n"
            "- 🔬 **3D로 가려면**: 같은 블록의 **연속 절편(serial sections)을 z-정합**하거나 "
            "3D 공간전사체 플랫폼 필요(이 데이터셋엔 없음). 확보 시 3D barrier score로 확장 예정.")


# --- top summary cards: barrier verdicts ----------------------------------
st.subheader("① 결론: 장벽 점수 (containment 검정)")
bcols = st.columns(max(len(barriers), 1))
for col, b in zip(bcols, barriers):
    bs = core["barrier"].get(b)
    verdict, vtext = barrier_verdict(bs)
    vbadge, vcolor = VERDICT_BADGE[verdict]
    with col:
        if bs and bs.get("z_score") is not None:
            col.metric(f"{b} — z", f"{bs['z_score']:+.1f}",
                       help="matched-count 랜덤 기질 대비 z. z>2면 개재 근거.")
            col.markdown(
                f"<div style='padding:.4rem .6rem;border-radius:.4rem;"
                f"background:{vcolor}22;border-left:4px solid {vcolor}'>"
                f"<b>{vbadge}</b><br><small>회랑 밀도 {bs['observed_corridor_density']:.2f} "
                f"vs 귀무 {bs['null_mean']:.2f} · p={bs['p_value']:.3g}</small></div>",
                unsafe_allow_html=True)
        else:
            col.metric(f"{b} — z", "N/A")
            col.markdown(f"<small>{vbadge}</small>", unsafe_allow_html=True)

fig_b = plots.fig_barrier(core["barrier"])
if fig_b:
    st.pyplot(fig_b)
st.caption("귀무가설은 동일한 **수**의 무작위 기질 세포를 가짜 장벽으로 재라벨합니다. "
           "따라서 이 검정은 양이 아니라 **위치**를 봅니다 — 80~90%가 기질인 조직에서 "
           "수를 맞추지 않으면 모든 것이 '장벽'으로 나옵니다.")


# --- 누가 종양에 인접한가 (rim 농축 패널) ----------------------------------
st.subheader("①-b 누가 종양에 인접한가 (rim 농축)")
st.markdown(
    f"장벽 점수가 *특정 후보*(myCAF/iCAF)를 검정한다면, 이 패널은 **모든 세포타입**의 "
    f"종양 주변({corridor_um:g}) 껍질 농축을 재서 *실제로 누가 종양에 인접한지* 데이터로 "
    "보여줍니다. 얇은 rim에 강건한 개선 지표(rim_enrichment).")
with st.spinner("rim 농축 계산 중…"):
    rim_rows = run_rim_panel(coords, labels, tumor=tumor, shell_um=corridor_um,
                             n_perms=n_perms, seed=0)
if rim_rows:
    rcol_l, rcol_r = st.columns([3, 2])
    with rcol_l:
        fig_r = plots.fig_rim(rim_rows, shell_um=corridor_um, tumor=tumor)
        if fig_r:
            st.pyplot(fig_r)
    with rcol_r:
        import pandas as _pd
        df_rim = _pd.DataFrame([{
            "세포타입": r["cell_type"],
            "rim z": round(r["z"], 1),
            "농축배수": round(r["enrichment_ratio"], 2) if r["enrichment_ratio"] else None,
            "판정": r["verdict"],
        } for r in rim_rows])
        st.dataframe(df_rim, hide_index=True, use_container_width=True)
        top = rim_rows[0]
        if top["z"] > 2:
            st.success(f"가장 종양 인접: **{top['cell_type']}** "
                       f"(z={top['z']:+.0f}, {top['enrichment_ratio']:.2f}배)")
        st.caption("z>2 종양 인접(농축) · z<−2 종양에서 배제. matched-null 대비 위치.")
else:
    st.info("rim 농축을 계산할 세포가 부족합니다.")


# --- proximity positive control -------------------------------------------
st.subheader("② 양성 대조 — 근접성 (myCAF가 iCAF보다 종양에 가까운가)")
pt = core.get("proximity_test")
if pt:
    c1, c2, c3 = st.columns(3)
    c1.metric(f"{barriers[0]} → {tumor} 중앙거리", f"{pt['median_dist_a']:.1f} µm")
    c2.metric(f"{barriers[1]} → {tumor} 중앙거리", f"{pt['median_dist_b']:.1f} µm")
    c3.metric("순열 p", f"{pt['p_permutation']:.3g}")
    closer = pt["interpretation"]
    expected = f"{barriers[0]} closer to {tumor}"
    if closer == expected:
        st.success(f"✅ {closer} — 문헌 기대(myCAF가 종양에 더 근접, Nat Genet 2024)와 일치")
    else:
        st.warning(f"⚠️ {closer} — 잘 재현된 문헌과 **불일치**. 세포 타이핑을 먼저 의심하세요.")
else:
    st.info("근접성 검정에는 장벽 후보 2종 이상이 필요합니다.")


# --- radial profile --------------------------------------------------------
st.subheader("③ 방사형 밀도 프로파일 — CD8 면역 배제")
col_l, col_r = st.columns([3, 2])
with col_l:
    st.pyplot(plots.fig_radial(coords, labels, tumor=tumor,
                               targets=tuple(list(barriers) + [immune])))
with col_r:
    rows = []
    for ct, mn in core["median_nn"].items():
        rows.append({"cell type": ct,
                     f"median dist → {tumor} (µm)": None if np.isnan(mn["median_um"]) else round(mn["median_um"], 1),
                     "n": mn["n"]})
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    st.caption("CD8 밀도가 거리에 따라 **상승**하면 종양 병소에서 T세포가 배제된 것.")


# --- tissue map ------------------------------------------------------------
st.subheader("④ 조직 지도 — 구조 확인 & 세포 타이핑 점검")
mcol, tcol = st.columns([3, 2])
with mcol:
    st.pyplot(plots.fig_tissue(coords, labels, title=src_desc))
with tcol:
    counts = pd.Series(core["cell_counts"]).sort_values(ascending=False)
    st.dataframe(counts.rename("cells").reset_index().rename(columns={"index": "cell type"}),
                 hide_index=True, use_container_width=True)
    st.bar_chart(counts)


# --- optional squidpy metrics ---------------------------------------------
if run_squidpy:
    st.subheader("⑤ 이웃 농축 & 동시출현 (squidpy)")
    if sq.get("error"):
        st.warning(f"squidpy 지표 실패: {sq['error']}")
    z = sq.get("nhood_zscore")
    cooc = sq.get("cooccurrence")
    scol1, scol2 = st.columns(2)
    with scol1:
        fig_n = plots.fig_nhood(z)
        if fig_n:
            st.pyplot(fig_n)
        elif sq.get("nhood_error"):
            st.caption(f"이웃 농축 실패: {sq['nhood_error']}")
    with scol2:
        fig_c = plots.fig_cooccurrence(cooc, focus=(tumor, immune))
        if fig_c:
            st.pyplot(fig_c)
        elif sq.get("cooccurrence_error"):
            st.caption(f"동시출현 실패: {sq['cooccurrence_error']}")
    st.caption("이웃 농축·동시출현은 **대칭적**이라 방향성이 없습니다 — "
               "'가깝다'는 말하지만 '사이에 있다(개재)'는 말하지 못합니다.")


# --- interpretation guide + downloads -------------------------------------
with st.expander("실데이터 결과 해석 가이드"):
    st.markdown("""
| 관찰 | 의미 |
|------|------|
| myCAF z > 2, iCAF z ≈ 0 | containment 지지 — 문헌과 일치 |
| 둘 다 z ≈ 0 | 기질 밀도 이상의 구조 없음 |
| myCAF가 iCAF보다 가깝지 **않음** | 믿기 전에 **세포 타이핑을 점검** — 잘 재현된 소견과 모순 |
| 전부 z > 2 | 귀무가설 오설정 가능 — 수 맞추기 확인 |

**정직한 한계:** 직선 경로 가정(실제 T세포 이동은 ECM 섬유·혈관을 따름), 정지 스냅샷,
상관≠인과, 그리고 모든 지표는 세포 라벨 하위에 있음 — 잘못된 주석은 확신에 찬 틀린 기하를 만듭니다.
""")

# downloadable JSON summary (drop non-serializable frames)
summary = {
    "source": src_desc, "platform": platform, "platform_level": lvl,
    "cell_counts": core["cell_counts"],
    "median_nn": {k: {"median_um": (None if np.isnan(v["median_um"]) else v["median_um"]),
                      "n": v["n"]} for k, v in core["median_nn"].items()},
    "proximity_test": pt,
    "barrier": core["barrier"],
}
st.download_button("⬇ 결과 요약 JSON 다운로드",
                   data=json.dumps(summary, indent=2, ensure_ascii=False,
                                   default=lambda o: None),
                   file_name="spatial_results.json", mime="application/json")
