# Platform-neutral rollback validation

This directory holds the re-validation done after the rollback of the native-Windows-only work
on 2026-08-26.

- clean CMake build PASS
- CTest 2/2 PASS
- persistent cache off / build / required: byte-identical
- previous persistent-cache tree vs platform-neutral final tree: pure water + red clay, 624
  geometries, byte-identical
- Rrs I/Q/U scatter: exact 1:1
- cache binary format: 5/5 byte-identical
- OpenMP / process row independence: PASS
- EAP policy: unchanged
- new OS-specific API calls in the persistent-cache core: 0

Raw full-grid outputs are kept in the external validation package; this source tree holds the
paired source CSV, metrics, figures, logs and the source diff.

---

# 플랫폼 중립 롤백 검증

이 디렉터리는 2026-08-26 native-Windows 전용 작업을 롤백한 뒤 수행한 재검증 자료다.

- clean CMake 빌드 PASS
- CTest 2/2 PASS
- 영속 캐시 off/build/required 바이트 동일
- 이전 영속 캐시 정본 대 플랫폼 중립 최종본: 순수해수 + red clay 624 기하 바이트 동일
- Rrs I/Q/U 산포도 정확히 1:1
- 캐시 이진 형식 5/5 바이트 동일
- OpenMP/process 행 독립성 PASS
- EAP 정책 변경 없음
- 영속 캐시 코어의 새 OS 전용 API 0건

원시 전체격자 출력은 외부 검증 패키지에 보존하고, 이 소스 트리에는 쌍 원본 CSV, 지표, 그림, 로그, 소스
diff 를 포함한다.
