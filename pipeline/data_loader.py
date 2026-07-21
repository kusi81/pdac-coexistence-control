"""
실데이터 로더 — squidpy 내장 공간 데이터셋을 파이프라인 형식(coords + labels)으로.

지금까지 합성 조직으로만 검증하던 장벽 점수를 **실제 공간 단일세포 데이터**에 적용한다.
MIBI-TOF(대장암) 데이터셋은 단일세포 해상도 + 세포타입 주석을 갖춰 barrier score에 적합.

주의
----
- 좌표 단위가 µm가 아닐 수 있다(MIBI-TOF는 픽셀). corridor_um은 좌표 단위 기준이므로
  해상도에 맞춰 해석해야 한다(아래 pixel_size로 스케일 정보 제공).
- 여러 시료(library_id)가 한 AnnData에 합쳐져 있으면 시료별로 분리해 분석해야 한다
  (서로 다른 좌표계가 섞이면 거리 계산이 무의미).
"""

import os
import glob
import numpy as np

# PDAC 마커 패널 (synthetic 생성기와 동일한 문헌 기반 마커) — 실데이터에 세포타입
# 주석이 없을 때 마커 기반으로 타이핑한다(Visium은 spot-level argmax = 근사).
PDAC_MARKERS = {
    "Tumor":      ["KRT19", "EPCAM", "MSLN", "KRT8", "KRT18"],
    "myCAF":      ["ACTA2", "TAGLN", "POSTN", "LRRC15", "COL1A1", "COL1A2"],
    "iCAF":       ["IL6", "CXCL12", "CFD", "PDPN", "DCN", "CLU"],
    "CD8_T":      ["CD8A", "CD3E", "CD3D", "GZMB", "GZMK"],
    "Macrophage": ["CD68", "CD14", "HLA-DRA", "LYZ", "AIF1"],
}

# obs에서 세포타입 라벨로 쓸 후보 키 (데이터셋마다 이름이 다름)
CELLTYPE_KEYS = ["cell_type", "Cluster", "cluster", "celltype",
                 "cell_types", "CellType", "annotation"]
LIBRARY_KEYS = ["library_id", "batch", "sample", "point", "fov"]


def _find_key(adata, candidates):
    for k in candidates:
        if k in adata.obs.columns:
            return k
    return None


def find_celltype_key(adata):
    k = _find_key(adata, CELLTYPE_KEYS)
    if k is None:
        # 폴백: 범주형 컬럼 중 카테고리 수가 적당한 것
        for col in adata.obs.columns:
            s = adata.obs[col]
            if str(s.dtype) == "category" and 2 <= len(s.cat.categories) <= 40:
                return col
        raise ValueError("세포타입 컬럼을 찾지 못했습니다. celltype_key를 직접 지정하세요.")
    return k


def find_library_key(adata):
    return _find_key(adata, LIBRARY_KEYS)


def load_mibitof():
    """squidpy MIBI-TOF 대장암 데이터셋 (인터넷 다운로드 필요)."""
    import squidpy as sq
    return sq.datasets.mibitof()


def load_imc():
    """squidpy IMC(imaging mass cytometry) 데이터셋 (대체용)."""
    import squidpy as sq
    return sq.datasets.imc()


def libraries(adata, library_key=None):
    """시료(library) 목록."""
    library_key = library_key or find_library_key(adata)
    if library_key is None:
        return [None]
    s = adata.obs[library_key]
    return list(s.cat.categories) if str(s.dtype) == "category" else sorted(set(s))


def to_coords_labels(adata, celltype_key=None, library_key=None, library=None,
                     spatial_key="spatial"):
    """AnnData → (coords[N,2], labels[N] str). library 지정 시 해당 시료만."""
    if celltype_key is None:
        celltype_key = find_celltype_key(adata)
    sub = adata
    if library is not None:
        library_key = library_key or find_library_key(adata)
        if library_key is not None:
            sub = adata[adata.obs[library_key].astype(str) == str(library)]
    if spatial_key not in sub.obsm:
        raise ValueError(f"obsm['{spatial_key}'] 좌표가 없습니다.")
    coords = np.asarray(sub.obsm[spatial_key], dtype=float)[:, :2]
    labels = np.asarray(sub.obs[celltype_key]).astype(str)
    return coords, labels


# 저자(FertigLab Xenium PDAC) 검증 module-score 마커 — 2단계 CAF 서브타이핑용.
# ① 대분류 모듈로 세포타입 argmax → ② CAF로 분류된 세포만 서브타입 argmax.
PDAC_COARSE_MODULES = {
    "Tumor":       ["EPCAM", "KRT19", "KRT8", "KRT18", "KRT17", "MUC1",
                    "CEACAM6", "TFF1", "TFF3", "MSLN"],
    "CD8_T":       ["CD8A", "CD3D", "CD3E", "GZMB", "GZMK", "NKG7", "CCL5"],
    "Macrophage":  ["CD68", "CD14", "LYZ", "ITGAM", "AIF1", "C1QA", "C1QB", "CD163"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "PLVAP"],
    "CAF":         ["FAP", "LUM", "DCN", "COL1A1", "COL1A2", "PDGFRB"],  # panCAF
}
CAF_SUBTYPE_MODULES = {
    "myCAF": ["TAGLN", "MYL9", "TPM2", "MMP11", "HOPX", "TWIST1", "SOX4", "ACTA2"],
    "iCAF":  ["CXCL1", "CXCL2", "CCL2", "LMNA", "HAS1", "HAS2", "IL6", "CXCL12"],
    "apCAF": ["CD74", "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "SLPI"],
}


def _module_scores(adata, modules, normalize=True):
    """각 모듈의 z-score 평균을 세포별 점수로. 반환 (scores[N,M], type_list, coverage)."""
    import scipy.sparse as sp
    X = adata.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = np.asarray(X, dtype=float)
    var = {g.upper(): i for i, g in enumerate(map(str, adata.var_names))}
    if normalize:
        lib = X.sum(1, keepdims=True); lib[lib == 0] = 1.0
        X = np.log1p(X / lib * 1e4)
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    types = list(modules)
    scores = np.full((X.shape[0], len(types)), -np.inf)
    npres = ntot = 0
    for j, t in enumerate(types):
        cols = []
        for g in modules[t]:
            ntot += 1
            gi = var.get(g.upper()) or var.get(g.replace("-", "_").upper()) \
                or var.get(g.replace("_", "-").upper())
            if gi is not None:
                cols.append(gi); npres += 1
        if cols:
            scores[:, j] = Z[:, cols].mean(1)
    return scores, types, npres / max(1, ntot)


def annotate_pdac(adata, key="cell_type"):
    """저자 방법(2단계 module-score) 재현: 대분류 후 CAF만 서브타입 배정.

    반환 (labels, coverage). myCAF/iCAF/apCAF를 원리적으로 구분 → argmax보다 신뢰.
    """
    cs, ctypes, cov1 = _module_scores(adata, PDAC_COARSE_MODULES)
    coarse = np.array([ctypes[k] for k in np.argmax(cs, axis=1)])
    labels = coarse.copy()
    caf_mask = coarse == "CAF"
    if caf_mask.any():
        ss, stypes, cov2 = _module_scores(adata, CAF_SUBTYPE_MODULES)
        sub = np.array([stypes[k] for k in np.argmax(ss[caf_mask], axis=1)])
        labels[caf_mask] = sub
    else:
        cov2 = 0.0
    adata.obs[key] = labels
    return labels, (cov1 + cov2) / 2


def annotate_by_markers(adata, markers=None, key="cell_type", normalize=True):
    """세포타입 주석이 없는 실데이터를 마커 기반으로 타이핑(argmax).

    각 타입의 마커 유전자 z-score 평균을 점수로, 최고 점수 타입을 배정.
    Visium(spot)에서는 spot당 우세 타입 근사(= deconvolution argmax의 단순판).
    adata.obs[key]에 라벨을 쓰고 labels 배열을 반환. 반환 (labels, coverage) —
    coverage는 패널 유전자 중 데이터에 존재하는 비율.
    """
    import scipy.sparse as sp
    markers = markers or PDAC_MARKERS
    X = adata.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = np.asarray(X, dtype=float)
    var = {g.upper(): i for i, g in enumerate(map(str, adata.var_names))}

    if normalize:                      # 원시 counts면 total-count + log1p
        lib = X.sum(1, keepdims=True); lib[lib == 0] = 1.0
        X = np.log1p(X / lib * 1e4)
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    Z = (X - mu) / sd

    types = list(markers)
    scores = np.full((X.shape[0], len(types)), -np.inf)
    n_present = n_total = 0
    for j, t in enumerate(types):
        cols = []
        for g in markers[t]:
            n_total += 1
            gi = var.get(g.upper()) or var.get(g.replace("-", "_").upper()) \
                or var.get(g.replace("_", "-").upper())
            if gi is not None:
                cols.append(gi); n_present += 1
        if cols:
            scores[:, j] = Z[:, cols].mean(1)
    labels = np.array([types[k] for k in np.argmax(scores, axis=1)])
    adata.obs[key] = labels
    coverage = n_present / max(1, n_total)
    return labels, coverage


# GSE274673 Xenium PDAC 시료: 부위코드 → (치료상태, GSM, 환자)
XENIUM_TREATMENT = {
    "31076": ("naive", "GSM8454446", "P1"), "37080": ("CRT", "GSM8454447", "P2"),
    "38245": ("CRT", "GSM8454448", "P3"),   "39928": ("naive", "GSM8454449", "P4"),
    "28429": ("naive", "GSM8454450", "P5"), "35406": ("CRT", "GSM8454451", "P6"),
}


def list_xenium_samples(data_dir=None):
    """받아진(추출된) PDAC Xenium 번들 목록 → [(bundle_path, label, treat)].
    label 예: 'P3 CRT치료후 (GSM8454448)'. 치료상태 미상이면 코드만."""
    if data_dir is None:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        data_dir = os.path.join(here, "data", "xenium")
    out = []
    for b in sorted(glob.glob(os.path.join(data_dir, "output-*"))):
        if not os.path.isdir(b):
            continue
        # cell_feature_matrix.h5 가 있어야 로드 가능(추출 완료)
        if not glob.glob(os.path.join(b, "**", "cell_feature_matrix.h5"),
                         recursive=True):
            continue
        code = None
        for cc in XENIUM_TREATMENT:
            if cc in os.path.basename(b):
                code = cc; break
        if code:
            treat, gsm, pt = XENIUM_TREATMENT[code]
            tk = "치료 전(naive)" if treat == "naive" else "CRT 치료후"
            out.append((b, f"{pt} · {tk} ({gsm})", treat))
        else:
            out.append((b, os.path.basename(b)[:30], "unknown"))
    return out


def load_xenium_bundle(bundle_dir):
    """10x Xenium 출력 번들 → AnnData(cells×genes, obsm['spatial']=µm 좌표).

    번들(tar.gz 해제) 안의 cell_feature_matrix.h5(발현) + cells.parquet/csv(좌표)만
    사용한다(transcripts.zarr·이미지는 불필요). Xenium 좌표는 µm 단위라 corridor_um이
    실제 마이크론으로 해석된다(Visium 픽셀과 달리).
    """
    import scanpy as sc
    import pandas as pd
    h5 = glob.glob(os.path.join(bundle_dir, "**", "cell_feature_matrix.h5"),
                   recursive=True)
    if not h5:
        raise FileNotFoundError(f"{bundle_dir}에 cell_feature_matrix.h5 없음")
    adata = sc.read_10x_h5(h5[0])          # gex_only 기본 True (유전자만)
    adata.var_names_make_unique()

    pq = glob.glob(os.path.join(bundle_dir, "**", "cells.parquet"), recursive=True)
    cz = glob.glob(os.path.join(bundle_dir, "**", "cells.csv.gz"), recursive=True)
    if pq:
        cdf = pd.read_parquet(pq[0])
    elif cz:
        cdf = pd.read_csv(cz[0])
    else:
        raise FileNotFoundError("cells.parquet / cells.csv.gz 없음")
    cdf["cell_id"] = cdf["cell_id"].astype(str)
    cdf = cdf.set_index("cell_id")
    obs_ids = adata.obs_names.astype(str)
    if set(obs_ids).issubset(set(cdf.index)):
        cdf = cdf.loc[obs_ids]             # cell_id로 정렬
    # (Xenium은 두 파일의 세포 순서가 동일 → 그래도 안 맞으면 위치 정렬)
    if len(cdf) != adata.n_obs:
        cdf = cdf.iloc[:adata.n_obs]
    adata.obsm["spatial"] = cdf[["x_centroid", "y_centroid"]].to_numpy(float)
    return adata


def load_zhou_pdac(data_dir=None, index=0, extract=True):
    """Figshare Zhou et al. PDAC Visium 데이터셋 로드.

    data/zhou_processed.zip(다운로드됨)을 풀어 index번째 h5ad를 읽는다.
    세포타입 주석이 없으면 annotate_by_markers로 타이핑해야 한다.
    """
    import anndata as ad
    if data_dir is None:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        data_dir = os.path.join(here, "data")
    zpath = os.path.join(data_dir, "zhou_processed.zip")
    outdir = os.path.join(data_dir, "zhou")
    if extract and not os.path.isdir(outdir):
        import zipfile
        with zipfile.ZipFile(zpath) as z:
            z.extractall(outdir)
    files = sorted(glob.glob(os.path.join(outdir, "**", "*.h5ad"), recursive=True))
    if not files:
        raise FileNotFoundError(f"{outdir}에 .h5ad가 없습니다. zip 압축을 확인하세요.")
    return ad.read_h5ad(files[index]), files


def guess_roles(labels, tumor_hint=("epithel", "tumor", "carcinoma", "malig"),
                immune_hint=("cd8", "tcell_cd8", "cytotoxic"),
                barrier_hint=("fibro", "caf", "stroma", "myofibro")):
    """라벨 문자열에서 종양/면역/장벽 역할을 휴리스틱으로 추정."""
    cats = sorted(set(labels))
    low = {c: c.lower() for c in cats}

    def _match(hints):
        for c in cats:
            if any(h in low[c] for h in hints):
                return c
        return None

    return dict(tumor=_match(tumor_hint), immune=_match(immune_hint),
                barrier=_match(barrier_hint), categories=cats)
