import numpy as np
SRC="/mnt/user-data/uploads/zo0_iop.txt"
# 자료 읽기 (\r\n, 머리말 건너뜀)
wl=[];aw=[];bw=[]
for line in open(SRC, encoding='latin-1'):
    line=line.strip()
    if not line or line[0] in '/!':
        continue
    p=line.split()
    if len(p)>=3:
        try:
            w=float(p[0]); a=float(p[1]); b=float(p[2])
        except ValueError:
            continue
        wl.append(w); aw.append(a); bw.append(b)
wl=np.array(wl); aw=np.array(aw); bw=np.array(bw)
print("읽은 점 수:", len(wl), " 범위:", wl.min(), "~", wl.max(), " 간격 중앙값:", np.median(np.diff(wl)))

# 1 nm 띠 평균: 각 정수 λ 에 대해 [λ-0.5, λ+0.5] 구간 사다리꼴 적분 / 1nm
def band_avg(lam_int):
    lo, hi = lam_int-0.5, lam_int+0.5
    m = (wl>=lo-1e-9) & (wl<=hi+1e-9)
    x=wl[m]
    if len(x)<2: return None, None
    # 사다리꼴 적분 후 폭으로 나눔
    ia = np.trapezoid(aw[m], x)/(x[-1]-x[0])
    ib = np.trapezoid(bw[m], x)/(x[-1]-x[0])
    return ia, ib

LAM = np.arange(350, 1201)  # 350~1200 nm, 1 nm (1100 요구에 여유)
out=[]
for L in LAM:
    a,b = band_avg(L)
    out.append((L,a,b))

# 검증: 앵커 밴드 띠평균 vs 단일점
print("\n앵커 밴드 띠평균(적분) 확인:")
for L in (443,555,660,850,900,1000,1100):
    a,b=band_avg(L)
    print("  %d nm: aw=%.6e  bw=%.6e"%(L,a,b))

# 새 표 쓰기
with open("/tmp/water_coef_z09_1nm_new.txt","w") as f:
    f.write("# pure-water IOP LUT (lambda_nm  a_w[m-1]  b_w[m-1])\n")
    f.write("# a_w = Pope&Fry 1997 (380-700nm) + Kou et al. 1993 (NIR); b_w = Zhang et al. 2009 (S=38.4 g/kg, T=20C)\n")
    f.write("# Source file: zo0_iop.txt (0.1nm, 200-2450nm). Backscatter uses 0.5*b_w.\n")
    f.write("# REBUILT 2026-07-27: full-range replacement, 0.1nm->1nm trapezoidal band-average over [lam-0.5, lam+0.5].\n")
    for L,a,b in out:
        f.write("%.1f %.6e %.6e\n"%(L,a,b))
print("\n새 표 기록: /tmp/water_coef_z09_1nm_new.txt  (%d 줄, %d~%d nm)"%(len(out),LAM[0],LAM[-1]))
