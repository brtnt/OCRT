# M50C — 공식 검증 제외 (diagnostic only), 영구 폐기 아님
근거 문서: M50C_POLICY_v0_9.md (사용자 기록; 패키지 외부). 등재: 2026-07-04 (v1.09 커밋 #3).

## 제외 사유 — lineage 오염 (모델 물리 문제 아님)
1. 83-angle 구판 M50C.mie 와 361-angle 재생성판이 동일 계보가 아님 —
   단순 각도 세분화가 아니라 물리적으로 다른 모델 (5/22 실측: bulk Ext_Co
   0.0253 vs 0.0366, P11 최대 30–40% 상이).
2. M50C.inp 수치 성분 fraction 이 주석 fraction 과 불일치 — 참 조성 판별 불가.
3. 통합 workbook 이 83-계보와 361-계보를 혼용.
4. legacy .mie 가 P11-only — vector 검증 부적합 (5/29 OSOAA 세션).
→ 어느 파일이 진짜 M50C 인지 참조 자체가 오염되어 검증 대상 자격 상실.

## 공식 aerosol 검증 subset
T50 + C50 + M80C (ocrt_ac_lut.py MODELS 기본값과 일치 — 확인됨 2026-07-04)

## 복권 조건 (5단계)
① .inp fraction 규약(number vs volume) 해소 → ② 361-angle M50C.mie 재생성
(vrt_solver Mie 생성기 가용 — C50.inp 재현 0.0000% 검증 완료) → ③ M50C.ext 재생성
→ ④ 6SV/OCRT/OSOAA 참조 재실행 → ⑤ lineage 청정 확인 후 재도입.

## 적용된 게이트
- 파일 위치: inputs/deprecated/ (production inputs/ 에서 제거)
- solver: --mie 경로에 "M50C" 포함 시 기본 거부, OCRT_DEBUG=1 로만 우회
  (diagnostic-only 상태의 코드 구현; option-hygiene 패턴)
- 검증 가이드라인 Tier 2 에 게이트 동작 검사 항목 추가
- 한계: 파일명 기반 매칭 — 파일을 개명해 투입하면 우회 가능. lineage 검증의
  최종 방어선은 이 게이트가 아니라 복권 5단계 절차임.
