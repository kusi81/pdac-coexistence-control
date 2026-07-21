# PDAC 공간 분석 대시보드 실행 (Windows PowerShell)
# 사용: 이 폴더에서  ./run.ps1
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$py = Join-Path $here ".venv\Scripts\python.exe"

if (-not (Test-Path $py)) {
    Write-Host "가상환경이 없습니다. 생성 중..." -ForegroundColor Yellow
    & py -3.13 -m venv (Join-Path $here ".venv")
    & $py -m pip install --upgrade pip
    & $py -m pip install -r (Join-Path $here "requirements.txt")
}

Write-Host "대시보드를 시작합니다 → http://localhost:8501" -ForegroundColor Green
& $py -m streamlit run (Join-Path $here "app.py")
