# OCRT 성능개선 최종 동결 지점

## 최종 포함 범위

1. T_wa incoming-U Q-fold 결함 수정
2. air-side Cox–Munk all-mode FKC
3. 실제 Fourier `m_max` build cap
4. air–water interpolation sort/dedup map memoization
5. versioned persistent Cox–Munk/interface operator cache
6. batch `p2b_row` cache worker-local화 및 row 간 lock 제거
7. persistent-cache core의 ISO C 플랫폼 중립화

## 폐기한 항목

- native-Windows 전용 speed gate
- Windows 전용 runtime 조립·build·PowerShell workflow
- 특정 OS에서만 공식 성능 판정을 인정하는 규칙
- OS별 물리·성능 code path

## 별도 재승인 전까지 금지

- output-only exact-view projection
- contracted `T_aw/T_wa` operator cache 재설계
- atmosphere/ocean subsystem R/T operator cache
- Jacobian partial-rebuild engine
- solver 구조의 추가 성능변경

## 병렬 계약

```text
independent processes
OMP_NUM_THREADS=1 per process
single prebuilder before simulation
persistent cache = immutable/read-only during rows
no first-row builder
no row-to-row wait/poll/lock
required-mode miss = immediate failure
```

EAP species scattering은 production runtime OFF이며 명시 선택은 fail-loud이다.
