# Platform-neutral rollback validation

이 디렉터리는 2026-08-26 native-Windows 전용 작업 롤백 후 수행한 재검증 자료이다.

- clean CMake build PASS
- CTest 2/2 PASS
- persistent cache off/build/required byte-identical
- 이전 persistent-cache 정본 대 platform-neutral 최종본: pure water + Red clay 624 geometry byte-identical
- Rrs I/Q/U 산포도 exact 1:1
- cache binary format 5/5 byte-identical
- OpenMP/process row independence PASS
- EAP policy 변경 없음
- persistent-cache core 신규 OS-specific API 0건

Raw full-grid outputs는 외부 validation package에 보존하고, 이 source tree에는 paired source CSV, metrics, figures, logs, source diff를 포함한다.
