# 과거 8,045-run 행렬의 현재 상태

현재 확보된 자료는 과거 campaign의 전체 수량과 공통 격자를 입증하지만, 다음 원표는 이 실행 패키지에 포함되어 있지 않다.

- 24개 aerosol model의 정확한 목록과 중앙 aerosol ID
- AOD865 7개 값
- water state 64개의 정확한 Chl/TSM/species/aDOM 값과 중앙 water ID
- 170개 wind-sensitivity run의 정확한 조건
- coupled 5,525개 run의 실제 command line/advanced numerical options
- thin 48-direction exact selection list

따라서 `HISTORICAL_RUN_MATRIX_REQUIRED.csv`는 의도적으로 빈 상태다. 이를 추정하여 채우지 않는다.

원본 `OCRT_paper_migration_2026-08-04.tar.gz`를 확보하면 다음 순서로 진행한다.

```text
python runner/inspect_historical_package.py --archive OCRT_paper_migration_2026-08-04.tar.gz
python runner/build_matrix_from_run_commands.py --commands <extracted run_commands.txt> --output matrices/HISTORICAL_RUN_MATRIX.csv
python runner/run_campaign.py --matrix matrices/HISTORICAL_RUN_MATRIX.csv ...
```

`build_matrix_from_run_commands.py`는 모든 row에 `--pssa`가 있는지 확인하고, water row에는 `--n-mu-water`가 명시되어 있지 않으면 중단한다.
