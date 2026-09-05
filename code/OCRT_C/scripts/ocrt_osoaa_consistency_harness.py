#!/usr/bin/env python3
# =============================================================================
#  OCRT vs OSOAA  in-water rrs(0-) consistency harness   (post-migration test)
# =============================================================================
#  PURPOSE
#    Session-start sanity check that an OCRT migration/build reproduces the
#    frozen in-water RT behaviour vs OSOAA (HYD.Model 3) across a mineral-TSM
#    matrix.  Regenerates EVERYTHING from a single mineral .mie -- nothing in
#    /tmp is assumed to survive a sandbox reset.
#
#  WORKFLOW
#    mineral.mie --converter--> OSOAA ExtData (vector phase) + bulk ext/sca
#       |-> OCRT side : fixed-bulk total IOP + BLENDED scalar P11 LUT
#       |               (mineral P11 + water Rayleigh P11, scattering-weighted)
#       |-> OSOAA side: HYD.Model 3 = UserProfile(a,b mineral) + ExtData(phase)
#       |               + OSOAA-internal pure-water Z09
#       run both (1 core, sequential, resumable) -> compare -> MAPE + scatter
#
#  HARD-WON TRAPS (do not "simplify" these away)
#    * OSOAA atm OFF is IMPOSSIBLE: AP.MOT 1e-4 -> entire radiance field NaN
#      (degenerate atmosphere poisons the in-water boundary condition).
#      AP.MOT=0.01 (thin) is the FLOOR.  Skylight effect is band-independent
#      and negligible (anchor Csed5/555 agrees 0.06% vs OCRT atm-OFF).
#    * 0- level: Flux.txt uses NUMERIC level 27 (Z=-0.0, Direct_Up=0).
#      Adv_UP "0-" is a TEXT line only -- data rows are numeric 27.
#      Adv_UP cols: level Z vza sca_angle I Q U  -> vza=s[2], I=s[4].
#      rrs(0-,nadir) = I(27, vza=0) / Total_Down(27).
#    * OSOAA_ROOT env var is MANDATORY (else ERROR_4000).
#    * n_mu_water >= 48 : floor for forward-peaked mineral phase. n_mu=32 is ~1.5%
#      under-converged vs n_mu=64 (Brown_earth Csed5/555nm); the floor rises with
#      single-scattering albedo so high-omega bands may need >48. (Was 32; raised
#      2026-06-30 after a convergence study.  Applies to LUT mode too.)
#    * VZA=0 nadir only needs Fourier m=0 (Ed = azimuth integral; nadir Lu at
#      mu=1 kills m>0).  Off-nadir would need water-side Snell + m>0.
#    * RAA=90 : glint-free, convention-invariant.
#    * Blend phase is the OCRT<->OSOAA parity: OSOAA blends Rayleigh INTERNALLY,
#      so OCRT must inject (b_min*P_min + b_w*P_Ray)/(b_min+b_w).
#    * --iop-bb is INERT in this phase-LUT validation path (SOS uses phase-embedded bb/b); only the
#      loader constraint bb < 0.5*b matters.
#    * HYD.Model 3 prints a warning that Chl/Csed are ignored -- EXPECTED.
#    * OMP threads are MOOT on 1 core (OMP=4 == OMP=1).  No process
#      parallelism helps either (nproc=1).
#
#  KNOWN RESIDUAL (not a migration failure)
#    MAPE ~1.6% with red/NIR (660,865) running -2..-5% (OCRT low).  Cause =
#    model3-cap(OCRT) vs OSOAA-truncation in the mineral phase, exposed at
#    red/NIR where b_w is tiny so the blend is mineral-dominated.  This is the
#    GOLDEN baseline; FAIL = materially worse than this.
# =============================================================================
import subprocess, os, sys, re, time, csv
import numpy as np
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))   # NumPy 2.x renamed trapz

# ============================== CONFIG =======================================
OCRT_BIN     = "/home/claude/ws/MIGRATION_FULL_2026-07-04_1917KST/ocrt/build/v2_solver_vk"
OSOAA_ROOT   = "/home/claude/ws/MIGRATION_FULL_2026-07-04_1917KST/osoaa"
CONVERTER    = "/home/claude/ws/MIGRATION_FULL_2026-07-04_1917KST/ocrt/scripts/mie_to_osoaa_extdata.py"   # mie -> ExtData
MINERAL_MIE  = "/home/claude/ws/MIGRATION_FULL_2026-07-04_1917KST/inputs/tsm_ahn/Brown_earth_AHN.mie"                # any of the 4 minerals
BANDS        = [412, 443, 490, 555, 660, 865]               # nm
CSED         = [0.5, 5.0, 50.0]                             # suspended mineral g/m3
# =============================================================================
# WIND-SPEED VALIDITY RULE (Jae, 2026-07-10):
#   OSOAA's rough-surface model has known accuracy issues at LOW WIND.
#   Cross-comparisons against OSOAA are EXCLUDED for wind < 3 m/s.
#   WIND = 3 m/s is the validated minimum (this harness's golden baseline
#   MAPE=1.90% was established at exactly 3 m/s).  Do NOT run OSOAA
#   comparisons at lower winds; low-wind anomalies (e.g. flat-surface
#   water-leaving deficits) are OSOAA-side artifacts, not OCRT defects.
# =============================================================================
SZA, VZA, RAA, WIND = 30.0, 0.0, 90.0, 3                    # nadir, glint-free
if WIND < 3.0:
    sys.exit("ERROR: wind %.1f m/s < 3 m/s - OSOAA comparison excluded by the "
             "wind-speed validity rule (see header). Use WIND >= 3." % WIND)
N_MU         = 48     # in-water Gauss nodes. >=48 recommended: n_mu=32 under-converged
                      # ~1.5% for forward-peaked mineral; high-omega bands may need >48.
AP_MOT       = 0.01                                         # OSOAA thin atm (FLOOR)
DELTA_RAY    = 0.039                                        # water Rayleigh depol (OCRT parity)
WORKDIR      = "/tmp/harness"
OUT_PNG      = "/mnt/user-data/outputs/consistency_ocrt_vs_osoaa.png"
OSOAA_TIMEOUT, OCRT_TIMEOUT = 240, 200
# Frozen Z09 pure water, 6-band (read inputs/water_iop/water_coef_z09_1nm.txt for others)
A_W = {412:4.5506e-3, 443:7.069e-3, 490:1.500e-2, 555:5.960e-2, 660:4.100e-1, 865:4.605}
B_W = {412:6.650e-3, 443:4.872e-3, 490:3.164e-3, 555:1.859e-3, 660:8.875e-4, 865:2.763e-4}
# Optional golden MAPE for regression-guard messaging (None = just report)
# GOLDEN_MAPE: must be re-established at n_mu=48.  The prior 1.60% was at the under-
# converged n_mu=32, where the Csed5/555 anchor +0.06% vs OSOAA was a convergence
# coincidence; at converged n_mu the OCRT-OSOAA difference is ~ -1.3% (true model3-cap
# vs OSOAA-truncation).  Run this harness at n_mu=48 to set the new golden, then restore
# a numeric value here.  None => regression-guard reports MAPE without a stale comparison.
GOLDEN_MAPE  = 1.90   # ESTABLISHED 2026-06-30 @ n_mu=48 (Brown_earth, 18 cases = 6 bands x Csed{0.5,5,50}).
                      # OCRT < OSOAA at ALL 18 cases (systematic negative bias). Baseline offset ~ -1 to -2%
                      # is the known OCRT model-3 phase-cap vs OSOAA phase-truncation difference (prior
                      # analysis). NOTE: high-absorption bands 660/865nm deviate MORE (worst -5.08% @
                      # Csed=50/865nm) -- NOT yet explained; flagged for investigation. Per rule 3 do NOT
                      # tune to close it. anchor Csed=5/555nm = -1.15%. Per-case golden:
                      # golden_ocrt_rrs_nmu48_2026-06-30.csv (bit-level reference file).
# =============================================================================


# ---------------------------------------------------------------- phase utils
def rayleigh_p11(theta_deg, delta=DELTA_RAY):
    """Depolarized Rayleigh scalar phase, normalized (1/2)*int P sin = 1.
       c2=(1-d)/(1+d); P=(1+c2 cos^2)/(1+c2/3).  P(0)=1.4713 for d=0.039."""
    c2 = (1.0 - delta) / (1.0 + delta)
    mu = np.cos(np.deg2rad(theta_deg))
    return (1.0 + c2 * mu * mu) / (1.0 + c2 / 3.0)

def backscatter_ratio(theta_deg, p11):
    """bb/b = (1/2) int_{90}^{180} P11 sin th dth, with P11 normalized."""
    th = np.deg2rad(theta_deg); s = np.sin(th)
    norm = 0.5 * _trapz(p11 * s, th)              # ~1 if normalized
    p = p11 / norm
    back = theta_deg >= 90.0
    return 0.5 * _trapz((p * s)[back], th[back])

# ---------------------------------------------------------------- ExtData I/O
def build_extdata(band, d):
    """Run converter mie->ExtData (cached). Returns path."""
    out = f"{d}/extdata/{band}.extdata"
    if not os.path.exists(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        r = subprocess.run(["python3", CONVERTER, MINERAL_MIE, str(band), out],
                           capture_output=True, text=True)
        if not os.path.exists(out):
            sys.stderr.write(f"ExtData FAIL {band}: {r.stderr[:200]}\n"); return None
    return out

def parse_extdata(path):
    """Return (ext, sca, theta[deg], F11[]) from an OSOAA IMOD=4 ExtData."""
    ext = sca = None; rows = []
    for ln in open(path):
        s = ln.split()
        if "EXTINCTION_COEF" in ln: ext = float(s[-1])
        elif "SCATTERING_COEF" in ln: sca = float(s[-1])
        elif len(s) >= 2:
            try: rows.append((float(s[0]), float(s[1])))
            except ValueError: pass
    a = np.array(rows); return ext, sca, a[:, 0], a[:, 1]

# ---------------------------------------------------------------- OCRT side
def build_ocrt_phase_lut(cs, band, theta, p_min, p_ray, b_min, b_w, d):
    """Blended scalar P11 LUT (mineral + water Rayleigh, scattering-weighted)."""
    out = f"{d}/ocrt_p11_cs{cs}_{band}.csv"
    blend = (b_min * p_min + b_w * p_ray) / (b_min + b_w)
    with open(out, "w") as f:
        f.write("theta_deg,P11\n")
        for t, p in zip(theta, blend): f.write(f"{t:.4f},{p:.6e}\n")
    return out

def run_ocrt(cs, band, A, B, BB, lut):
    """OCRT fixed-bulk single solve (atm OFF). Returns rrs0minus or None."""
    cmd = [OCRT_BIN, "--surface", "ocean", "--wind-speed", str(WIND),
           "--sza", str(SZA), "--vza", str(VZA), "--raa", str(RAA),
           "--wavelength", str(band), "--pressure", "0", "--aod", "0",
           "--water-model", "iop", "--iop-a", f"{A:.6f}", "--iop-b", f"{B:.6f}",
           "--iop-bb", f"{BB:.6f}", "--iop-phase-lut", lut, "--n-mu-water", str(N_MU)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           env={**os.environ, "OMP_NUM_THREADS": "1", "OCRT_ADVANCED": "1"}, timeout=OCRT_TIMEOUT)
    except subprocess.TimeoutExpired:
        return None
    m = re.search(r"rrs0minus=([0-9.eE+-]+)", r.stdout)
    return float(m.group(1)) if m else None

# ---------------------------------------------------------------- OSOAA side
def build_osoaa_profile(cs, band, a_min, b_min, d):
    """HYD.UserProfile: 5 header lines + (depth, a_abs, b_sca) mineral-only rows."""
    out = f"{d}/osoaa_prof_cs{cs}_{band}.txt"
    with open(out, "w") as f:
        f.write(f"# Mineral hydrosol Csed={cs} {band}nm (OSOAA adds Z09 water)\n#\n#\n# depth_m  a_abs  b_sca\n#\n")
        f.write(f"0.0 {a_min:.6f} {b_min:.6f}\n200.0 {a_min:.6f} {b_min:.6f}\n")
    return out

def run_osoaa(cs, band, ext_path, prof_path, d):
    """OSOAA HYD.Model 3 (thin atm AP.MOT). Returns rrs(0-,nadir) or None."""
    res = f"{d}/osoaa/{band}_{cs}"; os.makedirs(res, exist_ok=True)
    env = {**os.environ, "OSOAA_ROOT": OSOAA_ROOT, "HOME": "/home/claude",
           "OSOAA_NO_DIRECT_GLINT": "1"}
    surf = f"{OSOAA_ROOT}/DATABASE/SURF_MATR"
    cmd = [f"{OSOAA_ROOT}/exe/OSOAA_MAIN.exe", "-OSOAA.ResRoot", res, "-OSOAA.Log", "M.Log",
           "-OSOAA.Wa", f"{band/1000:.4f}", "-ANG.Thetas", str(SZA),
           "-ANG.Rad.NbGauss", "48", "-ANG.Mie.NbGauss", "100",
           "-AP.MOT", str(AP_MOT), "-AP.HR", "8.0", "-AER.AOTref", "0.0",
           "-HYD.Model", "3", "-HYD.ExtData", ext_path, "-HYD.UserProfile", prof_path,
           "-PHYTO.Chl", "0.0", "-SED.Csed", "0.0", "-YS.Abs440", "0.0", "-DET.Abs440", "0.0",
           "-SEA.Depth", "200.0", "-SEA.Ind", "1.34", "-SEA.Wind", str(WIND), "-SEA.Dir", surf,
           "-SEA.SurfAlb", "0.0", "-SEA.BotType", "1", "-SEA.BotAlb", "0.0",
           "-OSOAA.View.Phi", str(RAA), "-OSOAA.View.Level", "1",
           "-OSOAA.ResFile.vsVZA", "RESLUM_vsVZA.txt", "-OSOAA.ResFile.Adv.Up", "RESLUM_Adv_UP.txt"]
    try:
        subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=OSOAA_TIMEOUT)
    except subprocess.TimeoutExpired:
        return None
    fadv = f"{res}/Advanced_outputs/RESLUM_Adv_UP.txt"; fflux = f"{res}/Advanced_outputs/Flux.txt"
    if not (os.path.exists(fadv) and os.path.exists(fflux)): return None
    f2 = lambda x: (float(x.replace("D", "E")) if x not in ("NaN","nan") else float("nan"))
    V, I = [], []
    for ln in open(fadv):
        s = ln.split()
        if len(s) >= 7 and s[0] == "27":
            v, ii = f2(s[2]), f2(s[4])
            if v == v and ii == ii and v >= 0: V.append(v); I.append(ii)
    if not V: return None
    V, I = np.array(V), np.array(I); o = np.argsort(V); V, I = V[o], I[o]
    Lu = np.interp(0.0, V, I)
    Ed = None
    for ln in open(fflux):
        s = ln.split()
        if s and s[0] == "27" and len(s) > 4:
            tot, dr = f2(s[4]), f2(s[2]); Ed = tot if tot == tot else dr; break
    return Lu / Ed if Ed else None

# ---------------------------------------------------------------- driver
def load_done(path):
    done = {}
    if os.path.exists(path):
        for l in open(path):
            p = l.strip().split(",")
            if len(p) >= 3 and p[0] != "Csed":
                try: done[(float(p[0]), int(p[1]))] = float(p[2])
                except ValueError: pass
    return done

def main():
    os.makedirs(WORKDIR, exist_ok=True)
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    j = {b: i for i, b in enumerate(BANDS)}
    theta_ray = None

    # --- per-band mineral bulk + phase (from ExtData) ---
    band_data = {}
    for b in BANDS:
        ext_path = build_extdata(b, WORKDIR)
        if not ext_path: sys.exit(f"ExtData build failed for {b}")
        ext, sca, theta, F11 = parse_extdata(ext_path)
        if theta_ray is None:
            theta_ray = theta; p_ray = rayleigh_p11(theta)
        astar, bstar = ext - sca, sca                # mass-specific m2/g
        bb_b = backscatter_ratio(theta, F11)
        band_data[b] = dict(ext_path=ext_path, theta=theta, F11=F11,
                            astar=astar, bstar=bstar, bb_b=bb_b)
        print(f"[{b}nm] a*={astar:.4f} b*={bstar:.4f} bb/b={bb_b:.4f}")

    # --- sweep (resumable) ---
    ocrt_csv, osoaa_csv = f"{WORKDIR}/ocrt_results.csv", f"{WORKDIR}/osoaa_results.csv"
    ocrt_done, osoaa_done = load_done(ocrt_csv), load_done(osoaa_csv)
    if not os.path.exists(ocrt_csv):  open(ocrt_csv, "w").write("Csed,band,rrs0minus\n")
    if not os.path.exists(osoaa_csv): open(osoaa_csv, "w").write("Csed,band,rrs_osoaa\n")
    fo, fs = open(ocrt_csv, "a"), open(osoaa_csv, "a")

    t0 = time.time()
    for cs in CSED:
        for b in BANDS:
            bd = band_data[b]
            a_min = cs * bd["astar"]; b_min = cs * bd["bstar"]
            A = a_min + A_W[b]; B = b_min + B_W[b]
            BB = cs * bd["bb_b"] * bd["bstar"] + 0.5 * B_W[b]   # inert; loader only
            # OCRT
            if (cs, b) not in ocrt_done:
                lut = build_ocrt_phase_lut(cs, b, bd["theta"], bd["F11"], p_ray, b_min, B_W[b], WORKDIR)
                r = run_ocrt(cs, b, A, B, BB, lut)
                fo.write(f"{cs},{b},{r:.6e}\n" if r else f"{cs},{b},NA\n"); fo.flush()
                if r: ocrt_done[(cs, b)] = r
                print(f"  OCRT  Csed={cs:<4} {b}nm  rrs={r:.5e}" if r else f"  OCRT  Csed={cs} {b} FAIL  [{time.time()-t0:.0f}s]")
            # OSOAA
            if (cs, b) not in osoaa_done:
                prof = build_osoaa_profile(cs, b, a_min, b_min, WORKDIR)
                r = run_osoaa(cs, b, bd["ext_path"], prof, WORKDIR)
                fs.write(f"{cs},{b},{r:.6e}\n" if r else f"{cs},{b},NA\n"); fs.flush()
                if r: osoaa_done[(cs, b)] = r
                print(f"  OSOAA Csed={cs:<4} {b}nm  rrs={r:.5e}" if r else f"  OSOAA Csed={cs} {b} FAIL  [{time.time()-t0:.0f}s]")
    fo.close(); fs.close()
    print(f"sweep wall time: {time.time()-t0:.0f}s")

    # --- compare + MAPE + scatter ---
    rows = []
    for cs in CSED:
        for b in BANDS:
            if (cs, b) in ocrt_done and (cs, b) in osoaa_done:
                o, s = ocrt_done[(cs, b)], osoaa_done[(cs, b)]
                rows.append((cs, b, o, s, 100 * (o / s - 1)))
    if not rows: sys.exit("no comparable rows")
    pds = np.array([r[4] for r in rows])
    mape, mx, rmse = np.mean(np.abs(pds)), np.max(np.abs(pds)), np.sqrt(np.mean(pds**2))
    print(f"\nMAPE={mape:.2f}%  max={mx:.2f}%  RMSE={rmse:.2f}%  (n={len(rows)})")
    if GOLDEN_MAPE:
        tag = "PASS" if mape <= GOLDEN_MAPE * 1.05 else ("IMPROVEMENT" if mape < GOLDEN_MAPE * 0.95 else "REGRESSION")
        print(f"regression-guard vs golden {GOLDEN_MAPE:.2f}%: {tag}")

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    col = {CSED[0]:"#2563eb", CSED[1]:"#16a34a", CSED[2]:"#b45309"}; mk = {CSED[0]:"o", CSED[1]:"s", CSED[2]:"^"}
    allx = [r[2] for r in rows]; ally = [r[3] for r in rows]
    lo, hi = min(allx + ally) * 0.8, max(allx + ally) * 1.15
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))
    for cs in CSED:
        xs = [r[2] for r in rows if r[0] == cs]; ys = [r[3] for r in rows if r[0] == cs]
        ax[0].scatter(xs, ys, c=col[cs], marker=mk[cs], s=55, edgecolor="k", lw=0.5, label=f"Csed={cs} g/m3", zorder=3)
        ax[1].scatter(xs, ys, c=col[cs], marker=mk[cs], s=55, edgecolor="k", lw=0.5, zorder=3)
        xb = [r[1] for r in rows if r[0] == cs]; yb = [r[4] for r in rows if r[0] == cs]
        ax[2].plot(xb, yb, marker=mk[cs], color=col[cs], ms=7, lw=1.3, label=f"Csed={cs}")
    ax[0].plot([0, hi], [0, hi], "k--", lw=1, label="1:1"); ax[0].set_xlim(0, max(allx)*1.08); ax[0].set_ylim(0, max(ally)*1.08)
    ax[0].set_xlabel("OCRT rrs(0-) [1/sr]"); ax[0].set_ylabel("OSOAA rrs(0-) [1/sr]"); ax[0].set_title("(a) Linear scatter"); ax[0].legend(fontsize=8, loc="upper left"); ax[0].grid(alpha=.25)
    ax[1].plot([lo, hi], [lo, hi], "k--", lw=1); ax[1].plot([lo, hi], [lo*1.03, hi*1.03], "gray", lw=.7, ls=":"); ax[1].plot([lo, hi], [lo*0.97, hi*0.97], "gray", lw=.7, ls=":", label="+/-3%")
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlim(lo, hi); ax[1].set_ylim(lo, hi)
    ax[1].set_xlabel("OCRT rrs(0-) [1/sr]"); ax[1].set_ylabel("OSOAA rrs(0-) [1/sr]"); ax[1].set_title("(b) Log-log (multi-decade)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.25, which="both")
    ax[2].axhline(0, color="k", lw=.8)
    for y in (1, -1): ax[2].axhline(y, color="gray", lw=.6, ls=":")
    for y in (3, -3): ax[2].axhline(y, color="lightgray", lw=.6, ls="--")
    ax[2].set_xlabel("Wavelength [nm]"); ax[2].set_ylabel("100*(OCRT/OSOAA - 1) [%]"); ax[2].set_title("(c) Relative difference"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.25); ax[2].set_xticks(BANDS)
    fig.suptitle(f"OCRT vs OSOAA  in-water rrs(0-)  |  MAPE={mape:.2f}% max={mx:.2f}% (n={len(rows)}, nadir sza{SZA:.0f} raa{RAA:.0f} wind{WIND})", fontsize=11, y=1.0)
    plt.tight_layout(); plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    print(f"figure: {OUT_PNG}")

if __name__ == "__main__":
    main()
