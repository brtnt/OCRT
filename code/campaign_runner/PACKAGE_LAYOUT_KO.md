# 배포 파일 구성

- `base/`: 최신 migration base와 승인된 policy/code overlay
- `runner/`: setup, 16-worker campaign engine, resume/monitor, historical package inspector
- `config/`: PSSA/parallel/truncation/output policy
- `matrices/`: 즉시 실행 가능한 15-run Rayleigh matrix와 full historical matrix placeholder
- `docs/`: exact matrix blocker, output storage, scientific policy
- `manifests/`: core 및 PART00–PART10 SHA-256/size contract

이 core package 자체에는 1.4 GiB의 PART data를 중복 포함하지 않는다. `REQUIRED_EXTERNAL_DATA_PARTS.csv`에 기록된 exact archive를 함께 둔다.
