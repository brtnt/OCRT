#!/usr/bin/env python3
"""n_mu self-convergence check for a worst-case (high-omega) fixed-bulk case.
Usage: nmu_converge.py <band>   (Csed fixed = 50, the highest-scattering matrix column)
Reuses the harness-generated blended phase LUT and reproduces A,B (BB is inert).
"""
import subprocess, os, re, sys, glob, csv

band = sys.argv[1]
CSED = 50.0
NMUS = [48, 64, 80]
A_W = {'412':4.5506e-3,'443':7.069e-3,'490':1.500e-2,'555':5.960e-2,'660':4.100e-1,'865':4.605}
B_W = {'412':6.650e-3,'443':4.872e-3,'490':3.164e-3,'555':1.859e-3,'660':8.875e-4,'865':2.763e-4}

f = glob.glob(f'/tmp/harness/extdata/*{band}*')[0]
ext = sca = None
for line in open(f):
    U = line.upper()
    if 'EXTINCTION' in U: ext = float(line.split()[-1])
    elif 'SCATTERING' in U: sca = float(line.split()[-1])
astar, bstar = ext - sca, sca
A = CSED*astar + A_W[band]
B = CSED*bstar + B_W[band]
BB = 0.05*B                                   # inert (loader needs BB<0.5B); result unaffected
lut = f'/tmp/harness/ocrt_p11_cs{CSED}_{band}.csv'
assert os.path.exists(lut), f'phase LUT missing: {lut}'
print(f'band={band} Csed={CSED}  a*={astar:.5f} b*={bstar:.5f}  A={A:.5f} B={B:.5f} BB={BB:.5f}  omega={B/(A+B):.4f}')

res = {}
for nmu in NMUS:
    cmd = ['./build/v2_solver_vk','--surface','ocean','--wind-speed','3','--sza','30',
           '--vza','0','--raa','90','--wavelength',band,'--pressure','0','--aod','0',
           '--water-model','iop','--iop-a',f'{A:.6f}','--iop-b',f'{B:.6f}',
           '--iop-bb',f'{BB:.6f}','--iop-phase-lut',lut,'--n-mu-water',str(nmu)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           env={**os.environ,'OMP_NUM_THREADS':'1','OCRT_ADVANCED':'1'}, timeout=200)
        m = re.search(r'rrs0minus=([0-9.eE+-]+)', r.stdout)
        res[nmu] = float(m.group(1)) if m else None
    except subprocess.TimeoutExpired:
        res[nmu] = None
    print(f'  n_mu={nmu:>3}: rrs={res[nmu]}')

with open(f'/tmp/nmu_conv_{band}.csv','w',newline='') as fp:
    w = csv.writer(fp); w.writerow(['band','Csed','n_mu','rrs'])
    for nmu in NMUS: w.writerow([band, CSED, nmu, res[nmu]])

if all(v is not None for v in res.values()):
    r80 = res[80]
    print(f'  n_mu=48 vs 80: {(res[48]/r80-1)*100:+.3f}%')
    print(f'  n_mu=64 vs 80: {(res[64]/r80-1)*100:+.3f}%')
