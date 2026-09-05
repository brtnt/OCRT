"""pytest: sample detritus vs FR631 phase-grid spec (angle verbatim + normalization)."""
import os, math, re
ROOT=os.path.join(os.path.dirname(__file__),'..')
MIE=os.path.join(ROOT,'03_PYTHON_AND_NEW_DATA','Detritus_Stramski2001.mie')
SPEC=os.path.join(ROOT,'01_OCRT_C','docs','PHASE_ANGLE_GRID_FR631_SPEC_KO.md')
def _angles_from_spec():
    txt=open(SPEC,encoding='utf-8').read()
    block=txt.split('```')[1]
    return [l.strip() for l in block.splitlines() if l.strip()]
def _angles_from_mie():
    lines=open(MIE).read().splitlines()
    h=[i for i,l in enumerate(lines) if 'Phase Function (P11)' in l][0]
    out=[]; j=h+2
    while j<len(lines) and 'Phase' not in lines[j] and lines[j].split():
        out.append(lines[j].split()[0]); j+=1
    return out
def test_angle_grid_verbatim_631():
    a=_angles_from_mie(); b=_angles_from_spec()
    assert len(a)==631==len(b)
    assert all(float(x)==float(y) for x,y in zip(a,b))
def test_normalization_sampled():
    import numpy as np
    lines=open(MIE).read().splitlines()
    h=[i for i,l in enumerate(lines) if 'Phase Function (P11)' in l][0]
    wls=[float(x) for x in lines[h+1].split()[1:]]
    M=np.loadtxt(lines[h+2:h+2+631])
    ang=np.radians(M[:,0]); o=np.argsort(ang); ang=ang[o]
    for wl in (0.330,0.443,0.555,0.865,1.100):
        k=min(range(len(wls)),key=lambda i:abs(wls[i]-wl))
        P=M[o,1+k]
        n=np.trapezoid(P*np.sin(ang),ang)
        assert abs(n-2.0)<1e-3, (wl,n)
