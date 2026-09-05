# 실행 패키지의 과학·수치 정책

- PSSA: ON (`--pssa`를 모든 실행에 명시)
- PSSA 범위: 대기 direct-beam TOA→0+; 수중은 기존 plane-parallel geometry
- water particle phase: latest exact FR631, raw direct kernel
- broad water-particle truncation: OFF
- `0–0.005°` local cap: 사용하지 않음
- process-level parallelism: 16 independent OCRT processes
- per-process OpenMP: 1 thread
- `--n-mu-water`: advanced option이며 water run matrix row에서 값을 명시해야 함. runner가 임의 기본값을 넣지 않음.
- 신규 test/debug option: 없음
- `OCRT_DEBUG`: 사용하지 않음
