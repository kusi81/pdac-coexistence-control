"""
한글 폰트 등록 — matplotlib 그림에서 한글이 네모(□)로 깨지는 것을 방지.

import 만 하면 부작용으로 matplotlib rcParams가 설정된다:
    import fonts   # 한글 폰트가 잡혀 있으면 자동 적용

Windows의 맑은 고딕(Malgun Gothic)을 우선 사용하고, 없으면 다른 흔한 한글
폰트로 폴백한다. 하나도 없으면 조용히 넘어가(그림은 여전히 그려지고, 한글만
네모로 남음) 앱이 죽지 않도록 한다.
"""

import os
import matplotlib
import matplotlib.font_manager as fm

# 우선순위: 파일 경로(Windows) → 폰트 패밀리명(설치돼 있으면)
_FONT_FILES = [
    "C:/Windows/Fonts/malgun.ttf",       # 맑은 고딕 (Windows 기본)
    "C:/Windows/Fonts/NanumGothic.ttf",  # 나눔고딕
    "C:/Windows/Fonts/gulim.ttc",        # 굴림
    "C:/Windows/Fonts/batang.ttc",       # 바탕
    "/usr/share/fonts/truetype/nanum/NanumGothic.ttf",  # Linux
    "/System/Library/Fonts/AppleSDGothicNeo.ttc",       # macOS
]
_FONT_NAMES = ["Malgun Gothic", "NanumGothic", "AppleGothic", "Gulim", "Batang"]

KOREAN_FONT = None


def _register():
    global KOREAN_FONT
    # 1) 파일이 존재하면 font_manager에 등록하고 패밀리명 확보
    for path in _FONT_FILES:
        if os.path.exists(path):
            try:
                fm.fontManager.addfont(path)
                KOREAN_FONT = fm.FontProperties(fname=path).get_name()
                break
            except Exception:
                continue
    # 2) 파일을 못 찾았으면 이미 설치된 패밀리명 중에서 탐색
    if KOREAN_FONT is None:
        installed = {f.name for f in fm.fontManager.ttflist}
        for name in _FONT_NAMES:
            if name in installed:
                KOREAN_FONT = name
                break

    if KOREAN_FONT:
        matplotlib.rcParams["font.family"] = KOREAN_FONT
    # 한글 폰트는 유니코드 마이너스(−)가 없어 축 눈금이 깨지므로 ASCII '-' 사용
    matplotlib.rcParams["axes.unicode_minus"] = False
    return KOREAN_FONT


_register()
