#!/usr/bin/env python3
"""full grid에서 N개 랜덤 부분집합을 뽑아 draft grid를 만든다.

full_grid_design_v1.csv의 행 중 N개를 시드 고정으로 무작위 추출한다. draft는
full의 부분집합이므로, draft 결과는 full 결과의 일부와 정확히 일치한다(같은
case_id). case_id는 원본 그대로 유지한다.

사용법:
    python subsample_grid.py --in full_grid_design_v1.csv --n 100 --out draft_grid_100.csv
    python subsample_grid.py --n 100          # 기본 in/out 사용
"""
import argparse
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', default='full_grid_design_v1.csv',
                    help='원본 full grid CSV')
    ap.add_argument('--n', type=int, default=100, help='추출할 행 수')
    ap.add_argument('--seed', type=int, default=20260720, help='난수 시드')
    ap.add_argument('--out', default='', help='출력 CSV (기본 draft_grid_<n>.csv)')
    args = ap.parse_args()

    import numpy as np
    rng = np.random.default_rng(args.seed)

    with open(args.inp, newline='') as fh:
        rd = csv.reader(fh)
        header = next(rd)
        rows = list(rd)

    n = min(args.n, len(rows))
    idx = rng.choice(len(rows), size=n, replace=False)
    idx.sort()   # 원본 순서 유지
    pick = [rows[i] for i in idx]

    out = args.out or f'draft_grid_{n}.csv'
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(pick)

    print(f'[subsample] {len(rows)}행 중 {n}행 추출 (seed={args.seed}) -> {out}')
    print(f'  case_id 예시: {pick[0][0]}, {pick[1][0]}, ..., {pick[-1][0]}')


if __name__ == '__main__':
    main()
