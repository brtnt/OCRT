#!/usr/bin/env python3
"""EAP 17종 카탈로그 동결표.

이 표의 모든 값은 논문·부속 코드·문헌에서 확정된 것이며 임의 결정이 없다.
출처는 항목마다 주석으로 적는다.
"""
# 공통 상수 -----------------------------------------------------------------
N_MEDIA      = 1.334      # 논문 식 1, 노트북 셀 6
A_SOL_675    = 0.027      # Johnsen et al. 1994, 논문 본문 / 노트북 셀 10  [m^2/mg]
V_EFF        = 0.6        # 논문 본문 / 노트북 셀 6
CI_KG_M3     = 2.0        # 껍질 실수 굴절률이 물보다 커지는 유일한 값 (본 세션 판정)
NORM_NM      = 675.0      # 논문 식 1. 공개 자료의 655 nm는 색인 오류로 규명됨
HILBERT_NM   = (300.0, 1000.0, 1.0)   # 가장자리 오염 회피 (본 세션 판정)
PSD_STEP_UM  = 0.05       # 수렴 확인: 0.10 um 이하에서 결과 동일
PSD_MAX_MULT = 5.5        # 유효직경 배수

# 진핵 15종 + 홍조류 = 16종 --------------------------------------------------
EUK = dict(Vs=0.2, nshell0=1.10, ncore0=1.02)   # 논문 본문 / 역산으로 재확인

# 원핵 1종 (Matthews & Bernard 2013, Biogeosciences 10:8139) -----------------
PROK = dict(Vs=0.5, nshell0=1.12,
            ncore_pow=(333.0, -1.94, 0.82),      # n  = 333*L^-1.94 + 0.82
            kcore_pow=(2.28e7, -4.66, 1.08e-5))  # n' = 2.28e7*L^-4.66 + 1.08e-5

# id, 굴절률 형상 열 이름, 유효직경(um), 종류
CATALOG = [
    ( 0, 'Diatoms (pennate)',             6.0,  'euk'),
    ( 1, 'Chlorophytes',                  8.0,  'euk'),
    ( 2, 'Diatoms (centric)',             6.0,  'euk'),
    ( 3, 'Cryptophytes',                  6.0,  'euk'),
    ( 4, 'Cyanobacteria (blue)',          6.0,  'euk'),
    ( 5, 'Cyanobacteria (red)',           6.0,  'euk'),
    ( 6, 'Dinoflagellates',              24.0,  'euk'),
    ( 7, 'Eustigmatophytes',              6.0,  'euk'),
    ( 8, 'Haptophytes: Pavlovaceae',      6.0,  'euk'),
    ( 9, 'Pelagophytes',                  3.0,  'euk'),
    (10, 'Prasinophytes',                 3.0,  'euk'),
    (11, 'Prochlorococcus',               0.5,  'euk'),
    (12, 'Haptophytes: Prymnesiaceae',    4.0,  'euk'),
    (13, 'Raphidophytes',                24.0,  'euk'),
    (14, 'Rhodophytes',                   6.0,  'euk'),
    (15, 'Synechococcus spp',             1.2,  'euk'),
    (16, 'Microcystis-like (vacuolate)',  5.0,  'prok'),
]

# 종별 전체 유효직경 목록 (논문 Table 1). mie 파일은 전부 생산한다.
DEFF_ALL = {
    'Diatoms (pennate)':[6,12,24,48], 'Chlorophytes':[2,4,6,8],
    'Diatoms (centric)':[6,12,24,48], 'Cryptophytes':[2,6,12,24,48],
    'Cyanobacteria (blue)':[2,6,12,24], 'Cyanobacteria (red)':[2,6,12,24],
    'Dinoflagellates':[2,6,12,24], 'Eustigmatophytes':[2,6,12,24],
    'Haptophytes: Pavlovaceae':[2,6,12,24], 'Pelagophytes':[1,2,3,4],
    'Prasinophytes':[2,3,4,5], 'Prochlorococcus':[0.4,0.5,0.7,0.9],
    'Haptophytes: Prymnesiaceae':[1,2,3,4], 'Raphidophytes':[12,24,48,60],
    'Rhodophytes':[2,6,12,24,48], 'Synechococcus spp':[0.4,0.8,1.2,1.8],
    'Microcystis-like (vacuolate)':[3,4,5,6],
}
if __name__ == '__main__':
    n = sum(len(v) for v in DEFF_ALL.values())
    print('종 %d개, 유효직경 조합 총 %d개' % (len(CATALOG), n))
    for i, name, de, kind in CATALOG:
        print('  id=%2d  %-30s Deff=%-5g %s' % (i, name, de, kind))
