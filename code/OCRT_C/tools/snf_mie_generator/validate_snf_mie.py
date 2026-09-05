#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json, math, re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict
import numpy as np
from scipy.interpolate import PchipInterpolator

NUM_RE = re.compile(r'^[\s]*[+-]?(?:\d|\.\d)')
WLS = np.array([0.350,0.400,0.412,0.443,0.470,0.488,0.515,0.550,0.590,0.633,
                0.670,0.694,0.760,0.860,1.240,1.536,1.650,1.950,2.250,3.750])

@dataclass
class MieData:
    path: Path
    n_ang: int
    wavelengths: np.ndarray
    spectral: np.ndarray # nwl,6
    angles: np.ndarray
    blocks: Dict[str,np.ndarray]


def nums(line: str):
    out=[]
    for t in line.replace('D','E').replace('d','e').split():
        try: out.append(float(t))
        except ValueError: break
    return out


def read_mie(path: Path) -> MieData:
    lines=path.read_text(errors='strict').splitlines()
    n_ang=None; first_i=None
    for i,line in enumerate(lines):
        v=nums(line)
        if len(v)==1 and float(v[0]).is_integer():
            n_ang=int(v[0]); first_i=i; break
    if not n_ang: raise ValueError(f'{path}: n_ang missing')
    hi=next((i for i in range(first_i+1,len(lines)) if 'Wlgth' in lines[i]),None)
    if hi is None: raise ValueError(f'{path}: spectral header missing')
    wl=[]; sp=[]; i=hi+1
    while i<len(lines):
        v=nums(lines[i])
        if len(v)>=7:
            wl.append(v[0]); sp.append(v[1:7]); i+=1
        elif wl:
            break
        else: i+=1
    if not wl: raise ValueError(f'{path}: spectral rows missing')
    nwl=len(wl)
    rows=[]
    for line in lines[i:]:
        v=nums(line)
        if len(v)>=nwl+1:
            rows.append(v[:nwl+1])
    if len(rows)!=3*n_ang:
        raise ValueError(f'{path}: phase rows {len(rows)} != {3*n_ang}')
    a=np.array([r[0] for r in rows[:n_ang]],float)
    blocks={}
    for bi,name in enumerate(('P11','P12','P33')):
        rr=rows[bi*n_ang:(bi+1)*n_ang]
        aa=np.array([r[0] for r in rr],float)
        if not np.allclose(aa,a,rtol=0,atol=1e-12):
            raise ValueError(f'{path}: angle mismatch in {name}')
        blocks[name]=np.array([r[1:] for r in rr],float)
    return MieData(path,n_ang,np.array(wl),np.array(sp),a,blocks)


def sha(path: Path):
    h=hashlib.sha256();
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def check_one(m: MieData):
    mu=np.cos(np.deg2rad(m.angles))
    order=np.argsort(mu); mu=mu[order]
    p11=m.blocks['P11'][order,:]
    p12=m.blocks['P12'][order,:]
    p33=m.blocks['P33'][order,:]
    # Native-grid trapezoid reproduces the generator's tabulated g column.
    norm_native=0.5*np.trapezoid(p11,x=mu,axis=0)
    g_native=np.trapezoid(p11*mu[:,None],x=mu,axis=0)/np.trapezoid(p11,x=mu,axis=0)
    # Dense monotone-cubic integration is the scientific resolution check.  It
    # does not treat the 0.5-degree phase samples as piecewise linear in mu.
    aord=np.argsort(m.angles)
    theta=np.linspace(0.0,180.0,7201)
    pd=PchipInterpolator(m.angles[aord],m.blocks['P11'][aord,:],axis=0)(theta)
    tr=np.deg2rad(theta); st=np.sin(tr); ct=np.cos(tr)
    norm=0.5*np.trapezoid(pd*st[:,None],x=tr,axis=0)
    gcalc=np.trapezoid(pd*ct[:,None]*st[:,None],x=tr,axis=0)/np.trapezoid(pd*st[:,None],x=tr,axis=0)
    ext=m.spectral[:,4]; sca=m.spectral[:,5]
    i550=int(np.argmin(abs(m.wavelengths-.550)))
    derived=np.column_stack([ext/ext[i550],sca/sca[i550],sca/ext,g_native])
    sc_diff=np.abs(m.spectral[:,:4]-derived)
    phys12=np.max(np.maximum(np.abs(m.blocks['P12'])-m.blocks['P11'],0.0))
    phys33=np.max(np.maximum(np.abs(m.blocks['P33'])-m.blocks['P11'],0.0))
    idx0=int(np.argmin(np.abs(m.angles-0.0))); idx180=int(np.argmin(np.abs(m.angles-180.0)))
    finite=all(np.isfinite(x).all() for x in (m.wavelengths,m.spectral,m.angles,*m.blocks.values()))
    return {
        'model':m.path.stem,'n_ang':m.n_ang,'n_wl':len(m.wavelengths),
        'angle_min_deg':float(m.angles.min()),'angle_max_deg':float(m.angles.max()),
        'angle_step_max_deg':float(np.max(np.abs(np.diff(np.sort(m.angles))))),
        'wavelength_min_um':float(m.wavelengths.min()),'wavelength_max_um':float(m.wavelengths.max()),
        'finite':bool(finite),'min_ssa':float(m.spectral[:,2].min()),'max_ssa':float(m.spectral[:,2].max()),
        'min_p11':float(m.blocks['P11'].min()),
        'max_norm_abs_error':float(np.max(np.abs(norm-1))),
        'max_norm_abs_error_pct':float(100*np.max(np.abs(norm-1))),
        'max_g_table_vs_native_abs':float(np.max(np.abs(g_native-m.spectral[:,3]))),
        'max_g_table_vs_dense_phase_abs':float(np.max(np.abs(gcalc-m.spectral[:,3]))),
        'max_native_norm_abs_error_pct':float(100*np.max(np.abs(norm_native-1))),
        'max_spectral_internal_abs':float(np.max(sc_diff)),
        'max_mueller_violation_p12':float(phys12),'max_mueller_violation_p33':float(phys33),
        'max_abs_p12_endpoints':float(max(np.abs(m.blocks['P12'][idx0,:]).max(),np.abs(m.blocks['P12'][idx180,:]).max())),
        'sha256':sha(m.path),'bytes':m.path.stat().st_size,
    }


def compare_baseline(ref: MieData, gen: MieData):
    assert ref.n_ang==gen.n_ang and np.array_equal(ref.wavelengths,gen.wavelengths) and np.array_equal(ref.angles,gen.angles)
    p12mask=np.ones_like(ref.blocks['P12'],dtype=bool)
    p12mask[np.isclose(ref.angles,0.0)|np.isclose(ref.angles,180.0),:]=False
    all_ref=np.concatenate([ref.spectral.ravel(),ref.blocks['P11'].ravel(),ref.blocks['P12'].ravel(),ref.blocks['P33'].ravel()])
    all_gen=np.concatenate([gen.spectral.ravel(),gen.blocks['P11'].ravel(),gen.blocks['P12'].ravel(),gen.blocks['P33'].ravel()])
    exact=(all_ref==all_gen)
    return {
      'model':ref.path.stem,
      'spectral_max_abs':float(np.max(np.abs(ref.spectral-gen.spectral))),
      'P11_max_abs':float(np.max(np.abs(ref.blocks['P11']-gen.blocks['P11']))),
      'P12_interior_max_abs':float(np.max(np.abs(ref.blocks['P12'][p12mask]-gen.blocks['P12'][p12mask]))),
      'P12_endpoint_max_abs':float(np.max(np.abs(ref.blocks['P12'][~p12mask]-gen.blocks['P12'][~p12mask]))),
      'P33_max_abs':float(np.max(np.abs(ref.blocks['P33']-gen.blocks['P33']))),
      'exact_numeric_fraction':float(exact.mean()),
      'different_numeric_values':int((~exact).sum()),
      'total_numeric_values':int(exact.size),
      'operational_values_exact':bool(np.max(np.abs(ref.spectral-gen.spectral))==0 and np.max(np.abs(ref.blocks['P11']-gen.blocks['P11']))==0 and np.max(np.abs(ref.blocks['P12'][p12mask]-gen.blocks['P12'][p12mask]))==0 and np.max(np.abs(ref.blocks['P33']-gen.blocks['P33']))==0),
    }


def compare_legacy(old:MieData,new:MieData):
    # Interpolate high-resolution new data in angle onto the legacy grid.
    asc=np.argsort(new.angles)
    rows={'model':old.path.stem,'old_n_ang':old.n_ang,'new_n_ang':new.n_ang}
    for name in ('P11','P12','P33'):
        interp=PchipInterpolator(new.angles[asc],new.blocks[name][asc,:],axis=0,extrapolate=False)
        pred=interp(old.angles)
        d=pred-old.blocks[name]
        rms=np.sqrt(np.mean(d*d)); scale=np.sqrt(np.mean(old.blocks[name]**2))
        rows[name+'_nrmse']=float(rms/max(scale,1e-300))
        rows[name+'_max_abs_over_peak']=float(np.max(np.abs(d))/max(np.max(np.abs(old.blocks[name])),1e-300))
        # Exclude algebraic-zero endpoints for the absolute max report of P12
        if name=='P12':
            mask=~(np.isclose(old.angles,0)|np.isclose(old.angles,180))
            rows['P12_interior_nrmse']=float(np.sqrt(np.mean(d[mask,:]**2))/max(np.sqrt(np.mean(old.blocks[name][mask,:]**2)),1e-300))
    common=min(len(old.wavelengths),len(new.wavelengths))
    assert np.allclose(old.wavelengths[:common],new.wavelengths[:common])
    names=['NorExt','NorSca','SSA','g','Ext','Sca']
    for j,n in enumerate(names): rows[n+'_max_abs']=float(np.max(np.abs(old.spectral[:common,j]-new.spectral[:common,j])))
    return rows


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--generated',type=Path,required=True);ap.add_argument('--current',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    gen={p.stem:read_mie(p) for p in sorted(a.generated.glob('*.mie'))}
    checks=[check_one(gen[k]) for k in sorted(gen)]
    bas=[]
    for k in ('T50','C50','M80C','O99'):
        bas.append(compare_baseline(read_mie(a.current/f'{k}.mie'),gen[k]))
    leg=[]
    for k in ('M95C','M98C'):
        leg.append(compare_legacy(read_mie(a.current/f'{k}.mie'),gen[k]))
    dep=a.current/'deprecated'/'M50C.mie'
    if dep.exists(): leg.append(compare_legacy(read_mie(dep),gen['M50C']))
    result={'checks':checks,'baseline_reproduction':bas,'legacy_replacement':leg}
    (a.out/'validation.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    for name,rows in [('generated_quality.csv',checks),('baseline_reproduction.csv',bas),('legacy_replacement.csv',leg)]:
        with (a.out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    print(json.dumps(result,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
