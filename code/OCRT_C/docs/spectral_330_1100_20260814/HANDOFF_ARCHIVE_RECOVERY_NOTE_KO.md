# 첨부 handoff ZIP 복구·무결성 기록

## 1. 원본

```text
OCRT_SPECTRAL_DATA_330_1100_HANDOFF_2026-08-14(2).zip
SHA-256 a38679873150bbc94f5340f1154ee763d2f563efaf1cd2b526466ecca3376562
```

## 2. 발견된 문제

일반 ZIP reader에서 중앙 디렉터리를 읽을 수 없었다. 압축파일 후반부가 잘려 central-directory와 trailing directory group이 유실된 상태였다.

## 3. 복구 방식

local file header를 순차 검색하여 각 entry의 compressed payload를 복구했다.

복구 결과:

- 파일 866개
- runtime replacement payload: C 197 + Python 197
- 핵심 통합문서와 기술문서
- AHN/aerosol source/provenance의 복구 가능한 부분

## 4. 유실된 것으로 확인된 항목

원 README가 기술한 다음 후반부 directory group은 완전 복구되지 않았다.

```text
04_REPRODUCIBILITY_ARCHIVES
05_VALIDATION_AND_MANIFESTS
06_TOOLS
07_REFERENCES
```

또한 non-runtime AHN plot 1개가 불완전했다.

## 5. 보완조치

- 197 paired payload를 실제 C/Python 파일구조에서 다시 mapping
- 모든 파일의 SHA-256을 재계산
- 196 accepted와 1 conditional/rejected를 decision map으로 재구성
- data coverage/physicality validators를 현행 코드에서 재작성
- acceptance report와 methods/provenance 문서를 복구문서 및 실제자료에서 재작성

## 6. 신뢰범위

복구된 runtime payload는 C/Python paired SHA equality를 만족한다. Accepted 196개는 최종 설치본과 handoff payload가 일치한다. Conditional detritus candidate도 archive copy가 handoff와 일치한다.

원 archive 후반부에 있었을 것으로 예상되는 기존 manifest/tool의 원문 자체를 완전 복구했다고 주장하지 않는다.
