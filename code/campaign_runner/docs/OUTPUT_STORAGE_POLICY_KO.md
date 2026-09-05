# 출력 저장 정책

- 각 full-grid 결과는 검증 후 `CSV.GZ`로 저장한다. 압축은 무손실이며 이 파일이 authoritative source이다.
- 실행 중에는 임시 `.csv`를 만들지만, 구조·finite·격자 검사를 통과한 뒤 `.csv.gz.tmp`로 압축하고 atomic rename한다.
- 원시 CSV는 기본적으로 삭제한다. `--keep-raw`를 주면 보존한다.
- 각 run의 command, 환경변수, binary SHA-256, output SHA-256, 소요시간, 행 수를 `run_manifest.csv`에 저장한다.
- 자료를 재처리하여 그림을 갱신할 때에는 source CSV.GZ, metrics, scatter/residual PNG, plotting script, PLOT_DATA_MAP, manifest를 같은 revision에서 함께 교체한다.
- 전체 저장용량은 첫 16-run pilot 후 실제 평균 압축 크기로 산정한다. 실행 전 추정치를 임의로 고정하지 않는다.
