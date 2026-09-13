#!/usr/bin/env sh
# Verify the files of code/OCRT_C against SHA256SUMS.txt.
# Installed Mie tables (inputs/, GitHub Release data-v1) are checked as well; tables that are
# not installed are skipped. Regenerate the manifests with: python3 scripts/make_manifests.py
#
# 한국어: code/OCRT_C 의 파일을 SHA256SUMS.txt 와 대조한다. 설치된 Mie 표도 검사하며, 없는 표는 건너뛴다.
# 매니페스트 재생성은 python3 scripts/make_manifests.py 로 한다.
set -eu
cd "$(dirname "$0")"
sha256sum -c --quiet --ignore-missing SHA256SUMS.txt
echo "SHA256SUMS.txt: all present files verified"
