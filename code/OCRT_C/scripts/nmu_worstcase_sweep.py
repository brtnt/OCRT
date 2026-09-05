#!/usr/bin/env python3
# Pending 4: worst-case in-water n_mu convergence check before production LUTs.
# Tests whether n_mu=48 is converged for the highest-omega / most forward-peaked
# cases. Compares n_mu in {48,64,80}. Isolates n_mu (same truncation) so the
# red/NIR OSOAA-truncation offset cancels. Resumable via /tmp/nmu_sweep.csv.
# n_mu=48 rows must reproduce the golden CSV (cross-check of this setup).
import subprocess, os, csv, time, re

BIN   = "/home/claude/work/build/v2_solver_vk"
LUTD  = "/tmp/harness"                       # stored blend P11 LUTs from the harness
OUT   = "/tmp/nmu_sweep.csv"
BUDGET_S = 240                               # per-call wall budget (resumable)

# from harness build_all log (ExtData ext/sca -> a*,b*; backscatter fraction bb/b)
astar = {412:0.0802,443:0.0613,490:0.0400,555:0.0233,660:0.0124,865:0.0087}
bstar = {412:0.5586,443:0.5548,490:0.5459,555:0.5233,660:0.4737,865:0.3981}
bb_b  = {412:0.0103,443:0.0111,490:0.0125,555:0.0146,660:0.0177,865:0.0211}
A_W   = {412:4.5506e-3,443:7.069e-3,490:1.500e-2,555:5.960e-2,660:4.100e-1,865:4.605}
B_W   = {412:6.650e-3,443:4.872e-3,490:3.164e-3,555:1.859e-3,660:8.875e-4,865:2.763e-4}

# worst-case candidates: highest-omega (555,660 @ Csed50) + highest-g blue (412) + flagged 490
CASES = [(50.0,660),(50.0,555),(50.0,490),(50.0,412)]
NMUS  = [48, 64, 80]

def omega(cs,b):
    A=cs*astar[b]+A_W[b]; B=cs*bstar[b]+B_W[b]; return B/(A+B)

def run(cs,b,nmu):
    A=cs*astar[b]+A_W[b]; B=cs*bstar[b]+B_W[b]
    BB=cs*bb_b[b]*bstar[b]+0.5*B_W[b]        # exact harness value (inert; loader-only)
    lut=f"{LUTD}/ocrt_p11_cs{cs}_{b}.csv"
    cmd=[BIN,"--surface","ocean","--wind-speed","3","--sza","30","--vza","0","--raa","90",
         "--wavelength",str(b),"--pressure","0","--aod","0",
         "--water-model","iop","--iop-a",f"{A:.6f}","--iop-b",f"{B:.6f}",
         "--iop-bb",f"{BB:.6f}","--iop-phase-lut",lut,"--n-mu-water",str(nmu)]
    t0=time.time()
    r=subprocess.run(cmd,capture_output=True,text=True,
                     env={**os.environ,"OMP_NUM_THREADS":"1","OCRT_ADVANCED":"1"},timeout=600)
    m=re.search(r"rrs0minus=([0-9.eE+-]+)",r.stdout)
    return (float(m.group(1)) if m else None), time.time()-t0

done=set()
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT)): done.add((r['Csed'],r['band'],r['n_mu']))
else:
    open(OUT,'w').write("Csed,band,omega,n_mu,rrs,sec\n")

t_start=time.time()
for cs,b in CASES:
    for nmu in NMUS:
        key=(str(cs),str(b),str(nmu))
        if key in done: continue
        if time.time()-t_start>BUDGET_S:
            print(f"[budget {BUDGET_S}s reached -- resumable, re-run]"); raise SystemExit(0)
        rrs,sec=run(cs,b,nmu)
        with open(OUT,'a') as f: f.write(f"{cs},{b},{omega(cs,b):.4f},{nmu},{rrs:.6e},{sec:.0f}\n")
        print(f"Csed={cs} {b}nm w={omega(cs,b):.3f} n_mu={nmu}: rrs={rrs:.6e} [{sec:.0f}s]")

tot=len(CASES)*len(NMUS)
print(f"progress {len(done)+0}/{tot}" if len(done)<tot else "ALL DONE")
