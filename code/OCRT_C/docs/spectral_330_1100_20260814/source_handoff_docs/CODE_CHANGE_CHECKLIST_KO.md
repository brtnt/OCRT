# OCRT 330-1100 nm 코드 변경 체크리스트

## C implementation

- [ ] `rt_iop_organic.h`: Detritus 채택 시 `ORGANIC_WAVELENGTH_MIN_NM` 350 -> 330
- [ ] `main.c`, `rt_water_rt.c`: organic phase 오류문구와 범위검사 330-1100
- [ ] `rt_water_iop.c`: exact endpoint의 extrapolated flag에서 `<=/>=`를 `</>`로 수정
- [ ] `rt_iop_organic.c`: table/phase 범위 밖 silent endpoint hold를 상위 fail-loud로 차단
- [ ] `rt_iop_ahn_mineral.c`: 동일
- [ ] `rt_aerosol.c`: 동일
- [ ] `rt_absorption.c`: 40 x 771 동적 파장축, 330/1100 exact endpoint
- [ ] Mie reader: 16/22-node 파일과 361-angle matrix 정상
- [ ] EAP production path 비활성 확인
- [ ] version/release docs 갱신

## Python implementation

- [ ] `data/` paired files가 C와 byte-identical
- [ ] constituent evaluator가 330-1100 범위를 동적으로 읽음
- [ ] organic phase lower bound를 detritus 채택상태와 일치시킴
- [ ] phase interpolation은 source node 내부에서만 수행
- [ ] out-of-range silent endpoint hold 제거/상위 차단
- [ ] gas xsec 771-point axis 정상
- [ ] EAP 비활성 확인
- [ ] version/release docs 갱신

## Search commands

```bash
grep -RIn "ORGANIC_WAVELENGTH_MIN_NM\|outside 350-1100\|350.*1100" src

grep -RIn "xq <= x\[0\]\|xq >= x\[n - 1\]" src

grep -RIn "751\|350.0" src data tests
```

검색결과를 기계적으로 모두 변경하지 말고, production runtime·test·generator·historical document를 구분한다.
