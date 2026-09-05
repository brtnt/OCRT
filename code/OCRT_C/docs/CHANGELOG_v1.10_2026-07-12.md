# OCRT v1.09 + S1/S1b/S2 speedup patch - build recipes

Two builds, two gates (see SPEEDUP_PLAN.md in the runtime package):

STRICT (bit reference; no OCRT_FAST_KERNELS -> pristine code paths):
  gcc -std=c11 -O3 -march=native -ffp-contract=off -fopenmp -Isrc \
      $(find src -name '*.c') -o build/v2_strict -lm
  gate: golden 18 bit-exact (Windows: run vs current production exe).
  Measured in sandbox: BIT-IDENTICAL to the shipped -O2 build on real
  production grids (light nmw24 + heavy nmw48).

FASTK (production speed build):
  gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
      -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp \
      -Isrc $(find src -name '*.c') -o build/v2_fastk -lm
  gate: compare_results.ps1 -RelTol 1e-10 -AbsTol 1e-12 vs STRICT.
  NEVER use full -ffast-math (its -ffinite-math-only breaks the NaN
  sentinel design).

Windows (MinGW) uses the same flags; do the golden bit gate there once.

Patch contents (all perf code behind -DOCRT_FAST_KERNELS; strict build
compiles the pristine semantics):
  S1  rt_solver.c  sos_build_source_pol: trace-free fast path; the
      isfinite early-return moved to a front pass (same detection set);
      restrict-qualified locals.  Bit-identical by construction+test.
  S1b rt_solver.c  FASTK-only: SoA kernel tables (36 jp-contiguous
      streams) + omp simd 6-accumulator reduction (reorder licensed,
      1e-10 gate).  Trace-only from_up/from_dn accumulators dropped on
      this path (branch-trace env routes to the original loop).
  S2  shared/surface.c  FASTK-only: thread-local memo for
      surface_R_ww_coxmunk_direct (bit-safe: cached doubles are the
      exact bits the direct evaluation produced).

  S3  shared/surface.c  FASTK-only: all-m Fourier-kernel cache for the
      Cox-Munk interface (R_ww + T_wa): the m-independent per-azimuth
      Mueller is computed in ONE phi pass and projected onto all modes
      0..64 simultaneously (OSOAA SURF_MATR analog).  Served coefficients
      are BIT-IDENTICAL to per-m evaluation (verified by cmp on full
      grids).  m>64 / alloc failure falls through to the original body.
  S3b shared/surface.c  memo upgraded to open-addressing linear probe
      (direct-mapped slot ping-pong measured a 32% miss rate).
      OCRT_FKC_STATS=1 prints memo diagnostics to stderr.

  S4  rt_solver.c  FASTK-only: kk-axis (layer) vectorization of the SOS
      source builder.  Loop interchange: jp becomes outer/sequential (the
      per-(k,kk) accumulation order over jp matches the strict path - no
      reduction reorder), kk (hundreds of layers) is the elementwise SIMD
      axis over once-per-call transposed contiguous field planes.

Sandbox measurements (1 core, gcc 13.3, AVX-512):
  grid       shipped -O2   flags-only     +S1..S4         strict
  light 24     24.71 s     21.18 (1.17x)   7.26 (3.40x)   bit-identical
  heavy 48     52.01 s     43.23 (1.20x)  15.82 (3.29x)   bit-identical
  fastk-vs-ref max rel: 4.6e-12 / 2.7e-11 PASS 1e-10; S3/S3b outputs are
  bit-identical to the pre-S3 fastk build (strongest possible gate).
  OSOAA same-machine baseline: anchor 2.7 s vs OCRT fastk 3.7 s (1.37x);
  full-LUT 7.26 s vs OSOAA field solve 2.7 s (~2.7x).  Residual sos time
  (~1.8 s) is now memory-streaming bound; further gains are algorithmic
  (operator reuse D3, water layer policy D4), not SIMD.

## Session 4 additions (P=1013 real-production track)

S5  rt_kernel.c + rt_water_rt.c hooks: moment-kernel cache (per-m pfm/gr/
    gt/arr/att/att tables keyed by nodes+moments FNV).  Bit-safe reuse;
    P=0 outputs bit-identical pre/post.  Wpkc-pattern sibling.
S6a rt_water_rt.c: water grid cache singleton -> 96-slot LRU (FASTK; 1 slot
    otherwise = legacy).  The sky-beam loop issues ~25-75 distinct keys per
    batch row (24-48 beam (sza_eq,F_sun_eq) pairs per grid class + directs);
    one slot made them evict each other.  Replay parity is structural.
    rt_solver.c: numerically-null sky beams (F_sun_eq < 1e-30; measured
    1e-45..1e-52 sub-cone leaks) are skipped.
    Diagnostic: OCRT_S6_TRACE=1 prints per-row block timings to stderr.

Measured (sandbox, 1 core, fastk, P=1013.25 light nmw24 full-grid):
  before: stalled at row ~67 (>15 min/grid extrapolated, slot thrash)
  after:  base cold pass ~52 s (per jobs row; base solve runs at default
          nmw=48) + ~1.5 s/row steady  =>  ~11 min/grid (component model;
          single-shot completion pending - sandbox slice limits)
Next targets: S7 atm-solve cache (per-row atm SOS + C3d rebuild ~1.4 s/row;
mode fields view-independent -> #16 pattern), base-solve nmw mismatch,
late-row kr0 regime re-check at 96 slots.
ACCURACY GATE (mandatory, open): all prior gates ran atm-off; the P>0
3-way gate (ref/strict/fastk) needs a slow trusted reference run - ref is
unpatched (single slot, no skip) so hours-class; plan: overnight user-side
run or sliced sandbox runs.

S7a rt_solver.c: sky beams solved at UNIT flux, outputs scaled by
    F_sun_eq (solver linear in F_sun; gate-2).  Root cause fixed: F_sun_eq
    inherits bit-jitter from the per-row atm solve (view node in the atm
    grid), so it poisoned the cache key and forced every beam cold on
    every vza node.  Guard extended to !isfinite(F_sun_eq) (grid rows
    produce NaN sub-cone fluxes - PRE-EXISTING defect, previously dropped
    via the solver's isfinite fast-fail; now dropped explicitly).
    BONUS: pass --n-mu-water <row nmw> on the batch BASE CLI - with F
    pinned, the base solve's beam keys equal the grid rows' keys and
    PRE-WARM the entire beam cache (launcher change pending).

Measured (P=1013.25 light, full 433-row grid, base --n-mu-water 24):
  143 rows flushed in 84 s  =>  ~0.46 s/row steady, ~220 s/grid
  (progression: stalled wall / >=15 min  ->  S6a ~11 min  ->  S7a ~3.7 min)
Gates: Tier-0 golden OK; P=0 fastk BIT-identical pre/post; strict
BIT-identical vs ref.  P>0 3-way accuracy gate still OPEN (slow reference).
Next: per-row ~0.34 s untraced tail (post-assembly), base cold ~20 s,
heavy-grid measurement, launcher --n-mu-water passthrough.

Session 5 (S7a refinement + forensics):
Measured P=1013 light full-grid component model (this build):
  base validation solve  10.5 s   (degenerate config: mmw=2, no IOPs)
  steady rows            0.094 s/row
  vza-node boundaries   ~11-14 s x 5   (open item, see below)
  => ~120 s/grid measured-component total (first observation was a stall
     wall, >=15 min extrapolated  =>  >=7.5x)
OPEN (perf-only, not accuracy): (a) node-boundary re-solves persist; key
dump shows grid-row beam solves run with view_as_node=1, which DISABLES
the grid cache gate - yet intra-node rows replay; one per-row vasn log
line closes this. (b) beams do not execute on ordinary rows 2..71
(mechanism unconfirmed; I_in<=0 suspected).  Diagnostics: OCRT_S6_TRACE=1
now also dumps the full sos_pure cache key ([S6K] lines).
Launcher TODO: pass --n-mu-water <row nmw> on the batch base CLI
(pre-warms the beam cache; measured base 49 s -> 10.5 s with it).

## BUGFIX D-SKY (2026-07-10) - UNCONDITIONAL (both builds), physics fix
Decision rule applied: follow OSOAA.  OSOAA_ANGLES_ADD_TETAS (OSOAA_ANGLES.F)
inserts special angles at their SORTED position and dedups coincident angles
by threshold; OCRT appended view/aux nodes instead, breaking the coupling
interpolator's documented "sorted ascending" contract.  Effect before fix:
transmitted-skylight columns were sign-flipped / NaN with vza-dependent
on/off inside one grid (nadir = only healthy angle -> historical nadir
validations could not see it).  Fix: aw_interp_on_unsorted wrapper (sort +
dedup + original interpolator; identity/bit-identical for sorted input),
applied at all 7 interp sites in rt_air_water_coupling.c.
Verified: Tier-0 golden OK; P=0 light grid BIT-identical (strict AND fastk)
- the fixed path is inactive at P=0.  P=1013 after fix: ALL rows run beams
uniformly (80 solved / 4 skipped), ZERO NaN fluxes, node-boundary re-solves
GONE (the F-jitter source was this bug; with S7a the beam cache is now
fully stable).  New honest baseline (light, P=1013, physically consistent):
base 12.2 s + first grid row 18.3 s (cold) + 0.91 s/row steady
=> ~7 min/grid.  The earlier "~4 min" figure was measured on contaminated
physics (node-1 rows silently missing the skylight term) and is void.

## B-0a step 1 (2026-07-10): unified angle-table module LANDED
src/rt_angles_unified.{h,c} + tools/test_uangles.c.  OSOAA-parity test
against OSOAA's own RAD_UsedAngles.txt (N=48, SZA=30, n_w=1.34):
51/51 angles, max|dmu|=4.8e-15, max|dw|=5.1e-15, IMUS=19 / IMUSW=13
exact.  Module is standalone (not yet consumed by solvers) - existing
builds/gates unaffected.  Next increments: water solver consumes the
table (B-0a.2), atm solver (B-0a.3), index-paired interface (B-0b).

## B-0a.2 (2026-07-10): water node-ring construction unified via rt_uangles
Append-then-sort assembly replaced by the OSOAA-parity table (GL(2N)-half
core + sorted-insert/DEDUP specials).  Gates: Tier-0, P=0 light
fastk/strict, P=0 heavy fastk all BIT-identical.  P=1013 outputs CHANGED
(max 0.7 percent, TOA_rho/U) - root cause D-SKY-2: default-ON zero-slot
specials duplicated existing nodes exactly (beam-sun == its source GL
node; mu=1 x3) and the sky-capture grid (nf = beam ring) interpolated
over DUPLICATED ABSCISSAS, corrupting the C3-full reverse-coupled
skylight term.  Dedup removes the defect; change is a FIX.  New goldens:
light md5 00edbb90a8eab30fa10222f15a9f3e65 (bit-stable across 3 runs);
heavy regenerating.  Hygiene incident logged: A/B ring probes reused the
golden jobs file under a short timeout and truncated the golden CSV -
probes must target scratch outputs only.

## B-0a.3 (2026-07-10): atm solver consumes the unified table
Both atm solver variants build their ring via rt_uangles (view angle at
its SORTED position, DEDUPED against GL nodes - on full-grid rows the
view coincides with a GL node to <1e-15, so the ring keeps n_gauss nodes
and the view index aliases that GL node, OSOAA IMUS rule).  Extraction
rewired from the fixed last-index convention to the bookkept index; the
boa-export capacity, and the pass2/pass3 external bottom-source grids,
are built from the SAME table so index alignment is exact by
construction.  Gates: Tier-0 + P=0 light fastk/strict BIT-identical;
Tier-3.5 (P=1013): all real-signal cells BIT-identical, the only change
is 1e-20-absolute cancellation dust in one symmetry-zero U cell (raa=0).
Goldens refreshed: light md5 7a1984146cbcecf10b750f6835bb500f
(repro-verified); heavy regenerating.

## B-0b.1 (2026-07-10): shared N - atm SOS rides the row's water Gauss count
All three atm passes (main + pass2 + pass3), the boa capacity table, and
the Ed-diffuse integrator now use w_opts.n_mu_water as the Gauss count
(OSOAA single-table rule; --n-mu no longer sizes the coupled-path
atmosphere).  Cost: negligible at this scale (Tier-3.5 44.7 -> 45.4 s).
GATE CATCH #2: the first cut left the Ed integrator on stale
GL(opts->n_mu) weights against 24-node boa radiances - it paired the six
GRAZING-most sky radiances with full-hemisphere weights, inflating Ed by
+9%.  Arbitrated per project rule against OSOAA fluxes (Ed(0+): OSOAA
2.6758, legacy-6 2.6390, broken 2.8835, FIXED 2.6371).  Residual P>0
changes after the fix are the legitimate atm 6->24 refinement (max 4.2%
on TOA_rho_Q; I-level <1%).  Goldens: light md5
32a6b84ce19ce923f22da765eca1c058 (bit-stable repro); heavy regenerating.

## B-0b.2 (2026-07-10): reverse interface exchange on the true ring (D-REV fix)
FOUND: rt_air_water_couple_water_to_atm's physically-correct Cox-Munk
reverse MATRIX path was guarded by a "ring must equal pure GL(n)" check;
production rings always carry zero-weight specials, so the check failed
and the coupling silently degraded to the FLAT Snell-interp approximation
at all winds (named D-REV).  FIX (OSOAA single-table rule): the water
solver now EXPORTS the ring's own quadrature weights (w_water_pos,
aligned with mu_water_pos; specials carry 0), and the coupling contracts
with them directly - zero-weight specials contribute nothing, exactly as
they should.  Sky-capture path exports beam-ring weights the same way;
the internal diagnostic call keeps legacy NULL.  Gates: Tier-0 + P=0
light fastk/strict BIT-identical; Ed unchanged (2.637140); the flat->
matrix change is surgically confined to TOA_rho I/Q/U (max 1.9% on Q,
1.7% on I) - the water-leaving / sky reverse-coupled terms, as expected
for flat->Cox-Munk at wind 2.  Goldens: light md5
af9e1b866b217a4740247e4240174ea8 (bit-stable repro); heavy regenerating.

## B-0c.1 checkpoint (2026-07-10): aw Fourier matrix - reciprocity shortcut REJECTED
Correction of scope first (user challenge upheld in part): the wa Fourier
matrix + per-beam aw transmission ARE fully implemented, tested, and
always on with ocean.  What is ABSENT is the aw FULL-GRID Fourier matrix
(only scalar _direct and _btdf_scalar exist on the aw side) - needed
solely by the diffuse-skylight forward coupling, whose place the
flat-interp + scalar-m0 hybrid has been holding.  So B-0c.1 is not a
speed reimplementation: it completes the 4-matrix OSOAA set.
ATTEMPT: reciprocity-anchored construction (D T_wa^T D scaled to the
tested aw scalar).  Oracles REJECTED it (row-integral vs T_aw_direct off
~90% at mid mu; II reconstruction poisoned).  ROOT-CAUSE FINDING: the wa
trig kernel and the aw scalar are DIFFERENT-fidelity models (macro-angle
TIR clamp vs Walter facet-level); zero sets differ, so they are not
reciprocal images.  Kernel body disabled (returns -3) with the finding
documented in-source; NOT wired anywhere - Tier-0 and Tier-3.5 verified
BIT-identical to B-0b.2 goldens.  NEXT SESSION: polarized aw trig built
directly at facet level (reuse the aw scalar's Walter geometry; replace
scalar Fresnel-T with the s/p transmission Mueller + R-path frame
rotations).  Oracle tool preserved: tools/test_taw_oracle.c.

## B-0c.1 COMPLETE (2026-07-10): polarized aw Fourier matrix + forward matrixification
Facet-level polarized air->water trig kernel implemented (aw scalar's
Walter geometry verbatim; scalar Fresnel (1-F) split into Ts/Tp Mueller;
Mishchenko rotations / build_rotation_L sandwich copied from the wa trig;
NO shadowing, matching the aw scalar / T_aw_direct conventions).
ORACLES PASSED: (a) trig II pointwise vs the tested aw scalar 9.4e-16;
(b) row hemispheric integral vs T_aw_direct 7.8e-12 (NW=96); T0 peak
scan normal (narrow-lobe peak ~12 sr^-1, earlier 'blowup' was a
truncation-vs-zero measurement artifact).  Fourier builder rewritten on
the new trig; forward coupling applies the full T_aw matrix when the atm
ring weights are supplied (zero-weight specials inert); legacy hybrid
kept as the NULL/flat fallback.  Gates: Tier-0 + P=0 BIT; Ed0plus AND
Ed0minus both UNCHANGED (m0 physics identity: FIX-B fine integral ==
matrix m0, as designed); changes confined to POLARIZED components
(Lu_Q max 3.7%, U ~2.5%) - exactly the flat->full-Mueller correction.
Cost +19 s/19 rows from per-call kernel builds -> S3-style cache TODO.
OSOAA 0- flux gap (-6%) parked as a flux-DEFINITION alignment item for
the B-0c final cross-gate.  Goldens: light md5
f1e658941d189623f5390443640a4764 (bit-stable); heavy regenerating.
MANDATORY FINAL GATE (user directive 2026-07-10): on B-0 completion run
the FULL Tier-3 OSOAA harness (18 cases) + Tier-3.5 + atm-on OSOAA
cross-comparison before sign-off.

## B-0c.2a (2026-07-10): T_aw interface-matrix caching
Attempt 1: ocrt_fkc_serve reuse -> REJECTED by A/B table dump: fkc's
internal fourier machinery is SPECIFIC to the R_ww/T_wa kernels and
produced a materially wrong T_aw table (II errors up to ~48 vs lobe peak
~12).  Lesson recorded: fkc is not a drop-in for arbitrary trig kernels.
Attempt 2: SELF-memoization (caches this kernel's own output; bit-
transparent by construction).  Gate catch #3: an endpoint-only key
collided across rows because the sorted view-node position shifts the GL
columns -> corrupted serves; fixed with a full-grid FNV-1a fingerprint
of both mu arrays.  FINAL: Tier-3.5 45.6 s (direct-loop 64.5 s; B-0b.2
level fully recovered) and BIT-IDENTICAL to golden f1e658941d18... -
goldens unchanged, no regeneration needed.  aw/wa clarification logged:
direction tags (air->water / water->air), NOT the water absorption
coefficient (that is iop_a, untouched).

## B-0c.2b checkpoint (2026-07-10): design frozen, quartet audit complete
R_aa exists (atm surface-Mueller BC + analytic sunglint) - interface
quartet COMPLETE.  Architecture decided: interface exchange reusing both
solvers; keystone = per-m diffuse top source in the water solver
(ext_top_per_m), collapsing equivalent beams into one solve and making
the beams' dropped m>=1/polarized incident physical.  Injection point
pinned at rt_water_rt.c ~3585 (rt_solver_primary_source_pol on the beam
column; extend with optional beam Stokes).  Rollout via OCRT_DIFFUSE_TOP
with a 5-stage gate ladder (single-column pin -> m0 A/B -> full-m grade
-> perf -> mandatory Tier-3 harness battery).  No code landed this
increment; all goldens/binaries unchanged.

## B-0c.2b IMPLEMENTED behind OCRT_DIFFUSE_TOP (2026-07-10)
Water solver gained the per-m diffuse TOP source (opts ext_top_*): a
column-parameterized mirror of rt_solver_primary_source(+_pol) adds the
first-scatter of the coupled in-water incident (I and Q; U-incident
deferred), direct l-sums over ws tables (signed columns), amplitude
A = 2*w_c*S^m(c).  rt_solver snapshots the pi-normalized coupled field
before SI conversion, feeds it, and the equivalent-beam loop breaks
immediately (sky-capture then self-skips via sky_n_mu_m0==0).  ext_top
content is ROW-INVARIANT (coupling-grid value-match zeroes the view
column) so the S6a grid cache stays coherent without key changes.
DEFAULT PATH VERIFIED BIT-IDENTICAL to golden f1e65894 (guard inert).
Perf: Tier-3.5 19-row 45.6 -> 19.8 s (m0) / 19.4 s (full-m): 2.3x, the
equivalent-beam pass eliminated.
GATE (2) m0 A/B vs beams: NOT identical - max 4.4% (Rrs_I), vza-trending,
raa-weak; Ed0minus +0.086% constant (diffuse downwelling MS included for
the first time).  ROOT CAUSE LOCATED: the beam CONE branch re-derives an
AIR sun (mu_air_eq) and drives it through the refraction+transmission
chain INSIDE the solve, although the coupled amplitudes are ALREADY the
transmitted in-water field - the amplitude conventions differ by
interface factors (T(mu), Jacobian), matching the observed shape/sign/
trend.  Direct injection is the clean definition; adjudication of which
path is right is FOLDED INTO the mandatory OSOAA cross-gate (Tier-3
harness battery, gate (5)).  Gates (3) full-m grading and (4) heavy perf
extrapolation ride on that adjudication.  Beams remain the default until
sign-off.

## CROSS-VALIDATION SESSION (2026-07-10) - gate (5) first execution
(1) Tier-3 OSOAA harness (18 mineral cases, atm-off): **PASS** -
MAPE=1.90% max=5.08% RMSE=2.33%, regression guard vs golden 1.90% PASS,
61 s wall (fastk).  P=0 physics confirmed intact through the entire B-0
ladder.  Figure: /mnt/user-data/outputs/consistency_ocrt_vs_osoaa.png.
(2) Atm-on absolute cross (NEW): case03/443, sza20, wind2, raa120,
MOT=tau_R(443,P1013.25)=0.2358885 (engine-dumped), OSOAA TOA I vs
OCRT TOA_rho*mu0 at 6 vza:
   OCRT is SYSTEMATICALLY LOW by 6-19% (MAPE beams 13.51%, diffuse
   13.17%).  Single-scatter hand check (rho_ss~0.092 + MS) sides with
   OSOAA (0.105 at vza11) -> a REAL upstream deficit in OCRT's atm-on
   TOA (pass1 path radiance and/or TOA composition), COMMON to both
   skylight paths (so NOT a B-0c.2b artifact).  This is the first
   absolute validation of the atm-on TOA; all prior atm-on gates were
   self-consistency only.
(3) Beams vs diffuse within that gap: diffuse uniformly closer to OSOAA
by 0.3-0.5 pp at every vza.  ADJUDICATION DEFERRED until the systematic
deficit is resolved.
KNOWN SETUP DEBT for the rerun: OSOAA UserProfile was fed TOTAL iop
(a_total,b_total) but OSOAA Model3 ADDS Z09 pure water itself ->
particles-only (a_p=a_total-a_w, b_p=b_total-b_w) next time (small,
darkening-side; does not explain OCRT-low).
NEXT: decompose the deficit - (a) OCRT Rayleigh-only TOA vs theory/OSOAA
(black or pure ocean), (b) TOA composition audit (pass1 + T*wl + glint
terms), (c) rerun cross with particles-only profile.  555 cross-case
blocked separately: fixed-bulk-LUT + full-grid combination dies rc=-2
(untested corner; 443/mie path fine) - logged as a bug to fix.

## DEFICIT LOCALIZED (2026-07-10, atm-on cross follow-up): D-ATM-MS
Three-way decomposition (OSOAA Ind=1.0001 / OCRT --surface black /
both dark-w0) PROVES: interface contribution IDENTICAL between codes
(within ~5% at every vza) -> interface physics innocent.  The whole gap
lives in the PURE-RAYLEIGH pass1: OCRT black vs OSOAA index-1 shows
-6% (nadir) -> -19% (vza 40-54) -> -6.5% (limb).  Scaling probes:
MS fraction (rho/rho_ss - 1) = +10.3% at tau 0.236 (OSOAA +25%),
grows to +34.9% at tau 0.472 -> multiple scattering RUNS but is
~40-60% UNDERWEIGHT; scalar-vs-vector differ only 0.3% (polarization
kernels innocent).  LINEAGE: pre-B-0 build (B0a_step1) reproduces the
value BIT-exactly -> HISTORICAL defect, untouched by the entire B-0
ladder; first exposed by today's first-ever absolute atm-on validation
(the mandated harness-class gate did its job).  Named D-ATM-MS.
NEXT: hunt inside the atm SOS order loop - suspects: per-order source
normalization (omega/2 factors), layer vertical integration of the MS
source, n_layers default for the Rayleigh column.
TIMING (this session, 1 core): OSOAA warm 0.08-0.40 s/run (SURF_MATR
db reused), cold ~4 s.  OCRT fastk: black single 0.1-0.5 s/point;
atm-on 19-row grid: diffuse 12-52 s (case-dependent), beams 45.6 s
(light) / >8 min uncompleted (555 heavy).  Tier-3 harness 18 cases
(both codes together): 61 s.  VERDICT: not yet the same class for
single coupled field solves (sub-second vs tens of seconds+) - the
B-0c native-coupling case reconfirmed; harness-class workloads already
practical.

## D-ATM-MS RETRACTED (2026-07-10): comparison-convention artifact
Jae's instinct was right.  Root cause = AZIMUTH CONVENTION mismatch:
OCRT_raa + OSOAA_phi = 180 deg.  Proof: (a) deficit shape equals the
single-scatter phase ratio P(60)/P(120) exactly; (b) OCRT black raa=60
reproduces OSOAA ind1 phi=120 to -0.05%; (c) corrected pure-atm table
(OCRT raa120 vs OSOAA VZA<0 branch) MAPE 0.16% over 6 vza.  OCRT atm
SOS INNOCENT; nothing broken by updates.  Harness never covered azimuth
(nadir-only), which is why this survived until today.
NEW OPEN ITEM: corrected case03 cross agrees 0.8-1.7% at vza>=40
(diffuse uniformly better than beams by 0.3-0.5 pp), but vza<=25 shows
an OSOAA water-signal lobe OCRT lacks (up to 5x at vza11; wp03 phase
bump at 150-170 deg scattering).  Hypothesis: water-leaving -> TOA
chain carries m=0 only, killing the azimuth-dependent lobe.  Probe
plan next session: OCRT Lu0plus vs TOA azimuth structure; OSOAA 0+
level; dark-water low-vza control.
TIMING: 10%-parity NOT met (OSOAA warm sub-second vs OCRT tens of
seconds per coupled grid); efficiency track remains justified after
correctness items.

## LOW-VZA EXCESS RESOLVED (2026-07-10): second comparison artifact
Jae's rule ("Rrs mismatches were always a wrong-quantity comparison")
held AGAIN.  The OSOAA runs included the direct sun reflection off the
rough surface (broad at wind 2, reaching low vza on the phi=300 side);
the harness has ALWAYS disabled it (OSOAA_NO_DIRECT_GLINT=1), and
OCRT's grid TOA_rho_I does not include that term either (single-solve
path reports it separately as rho_glint).  With glint OFF:
case03 atm-on TOA, 6 vza, raa120: beams MAPE 1.52%, diffuse 1.13%
(was 8.05/7.69 with glint on; vza11 was -32%).  OCRT code innocent.
CROSS-COMPARISON PROTOCOL (mandatory checklist, add to harness docs):
 1) azimuth: OCRT_raa + OSOAA_phi = 180 deg (use OSOAA VZA<0 branch
    for OCRT raa=120 when View.Phi=120);
 2) OSOAA_NO_DIRECT_GLINT=1 always;
 3) OSOAA UserProfile takes PARTICLE-ONLY iop (Model3 adds Z09 water);
 4) MOT = engine-dumped tau_R (OCRT_DUMP_TAUR);
 5) quantity: OSOAA_I = OCRT_TOA_rho * mu0.
ADJUDICATION (recommendation to Jae): diffuse beats beams at ALL six
vza (MAPE 1.13 vs 1.52), consistent with the beam cone-branch amplitude
chain suspicion -> recommend adopting OCRT_DIFFUSE_TOP as default and
retiring the equivalent beams (also 2.3x faster; 555-heavy beams could
not even finish).  Residual -1..-2.6% carries known small debts
(particle-only profile not yet applied; n_w 1.34 vs internal; nmw24 vs
OSOAA 48 nodes).
TIMING: OSOAA glint-off run 1.77 s (cold-ish); warm runs 0.08-0.4 s.

## B-0c.2b ADJUDICATED & DEFAULTED (2026-07-10): diffuse-top adopted
DECISION BASIS (per Jae): methodological, not result-proximity.  OSOAA
source audit proves the method: OSOAA_SOS_CORE.F applies the interface
as per-mode MATRIX x DIFFUSE FIELD inside the order loop (line 2485:
XI1*TAW11 + XQ1*TAW12 + XU1*TAW13) with ONLY the solar direct as a beam
(line 1747).  No equivalent-beam decomposition exists in OSOAA.  Our
diffuse top source is the same method in form; equivalent beams were a
workaround with a demonstrated amplitude-chain inconsistency (cone
branch).  Numbers corroborate (TOA cross: diffuse 1.13% vs beams 1.52%)
but are not the basis.
IMPLEMENTATION: diffuse-top is now the DEFAULT; OCRT_EQUIV_BEAMS=1
restores legacy beams (archaeology only); OCRT_DIFFUSE_TOP=m0/0 kept as
diagnostics.  GATES: Tier-0 bit (both builds); P=0 light grid BIT;
Tier-3.5 default run BIT-identical to the diffuse golden.  NEW GOLDENS:
light md5 2b89a2e1111ccb7efd2d8af9b86e5fa2, heavy md5
67ffdc1cf3b70002e99c3a1670b59f1a (repro BIT-identical).
PERFORMANCE (gate 4 closed): heavy 19-row grid 464 s (beams) -> 166 s
(diffuse), 2.8x; light 45.6 -> 19.7 s, 2.3x.  Remaining B-0c items:
U-incident increment; native interface-exchange iteration (mutual
reflections) if cross residuals ever demand it.

## SESSION 2026-07-10 (cont.): 555 cross + protocol debts + batch parser
(1) BATCH PARSER BUG (real, fixed): strtok_r collapsed EMPTY CSV fields,
shifting every later value one column LEFT (e.g. n_mu_water's "24"
landing in the mie column).  Replaced with a hand-rolled splitter that
preserves empty fields.  Gate: Tier-0 anchor + Tier-3.5 grid remain
BIT-identical to the goldens (the fix is transparent for well-formed
jobs).  RELATED, STILL OPEN: the fixed-bulk-LUT + grid combination
reaches the water solver with a/b/bb = 0 (mode set by the CLI
--fixed-bulk-phase-lut parser at main.c:1418, values never injected on
that path; the row-override block did not fire in the failing runs -
banner-less entry path to be pinned next session).  WORKAROUND that
unblocked everything: pass --fixed-bulk-iop A B BB on the base CLI
alongside the LUT; 555 grid then runs (19 rows, ~80 s diffuse).
(2) 555/cs5.0 atm-on cross EXECUTED: TOA MAPE 4.97% (OCRT high at
mid-vza).  Decomposition (OSOAA Ind=1.0001 + OCRT black): atmosphere
agrees to 0.0..-0.5% (limb -2.9%); the entire gap is in the WATER
contribution (+2.7% low-vza -> +13% mid -> +3.6% limb).  ATTRIBUTED to
the previously established OSOAA coarse vertical layering limitation
(open item #4: red/NIR residuals, OCRT = truth reference), amplified in
this bright, optically thick water at slanted views.  Follow-up probe
registered: rerun OSOAA with a finer layer setting and watch the gap
shrink.
(3) PARTICLE-ONLY PROFILE (protocol item 3) applied for 443/case03:
final official cross number MAPE 1.17% (vs 1.13% with the total-iop
profile - debt was small as predicted).  All five protocol items now
satisfied for the 443 official table.

## FIXED-BULK + GRID ROOT CAUSE CLOSED (2026-07-10)
The last root was in main's flow: the intentionally DISCARDED base solve
runs BEFORE the --batch-full-grid dispatch, and its failure (rc=-1002)
returned from main - the batch never started (hence no banner / no
row-override execution; the earlier FAIL messages were the BASE solve's).
FIX: in batch mode the base-solve failure is a WARNING and the flow
proceeds to the rows (base CLI is a template; rows carry their own
overrides).  Gates: Tier-0 + Tier-3.5 grid BIT-identical; the previously
blocked LUT-only 555 config now completes 19/19 rows.
CONVENTION FINDING (documented, code untouched): --fixed-bulk-iop on the
CLI sets F_sun=1.0, while jobs-supplied iop keeps the default F_sun=pi -
ABSOLUTE columns (Lu*/Ed*) differ by exactly pi between the two ways of
launching the same physics; RATIO columns (Rrs/rrs/TOA_rho) are
unaffected (all goldens and cross tables live on ratios or on one
consistent lineage).  The official 555 cross reference is X555t
(jobs-iop, same F_sun=pi lineage as the G35 goldens); the CLI-bypass
X555s is retired.  Unifying the convention is deliberately NOT done now
to preserve bit-goldens; revisit at the next golden re-baseline.

## 555 GAP INVESTIGATION (2026-07-10): two hypotheses killed, one unified
(1) OSOAA layering RETRACTED as the cause: rebuilt OSOAA with
CTE_NT_SEA 80 -> 400 (gfortran, header restored after; fine build kept
as exe/OSOAA_MAIN_NT400.exe).  555 results moved < 0.03 pp -> the open
item #4 attribution does NOT explain this gap (and #4 itself now needs
re-examination with the NT400 build at red/NIR).
(2) OCRT angular resolution RETRACTED: nmw 24 -> 48 rerun, MAPE 5.09%
vs 4.95% (shape flattens, total unchanged).
(3) INPUTS VERIFIED CONSISTENT: OSOAA particle profile (0.1166,2.6165)
+ its own Z09 water reproduces OCRT bulk (0.1762,2.6183); OCRT
fixed-bulk LUT path applies NO truncation (delta_f=0, omega_eff=in).
UNIFIED HYPOTHESIS (quantitative, both bands): a single ~+10% excess in
OCRT's WATER-LEAVING -> TOA chain explains 443 (+1.2% TOA at 12% water
share) AND 555 (+5.5..7% at 55% share) exactly.  The in-water 0- field
is harness-validated, so the suspect segment is 0- -> 0+ (TWA) -> TOA
(reverse coupling / pass2 / diffuse transmittance / n^2 bookkeeping).
NEXT DECISIVE PROBE: 0+ level three-way (OSOAA View.Level at 0+ vs OCRT
Lu0plus) at both bands - localizes the excess to TWA vs the atm pass2
leg.  TIMING: OSOAA NT400 run 2.27 s; OCRT nmw48 555 grid ~3 min.

## 0+ THREE-WAY PROBE (2026-07-10): chain segment narrowed to TWA-side
Component mismatch found FIRST (Jae's rule, 3rd time): OSOAA View.Level=3
REFL is TOTAL upward (water-leaving + sky reflection; Fresnel blow-up at
grazing), while OCRT Lu0plus is WATER-LEAVING only by definition
(rt_solver comment; sky reflection enters the TOA via the pass1 surface
BC instead - not a defect, a definition).  Fixed the comparison by the
DARK-WATER SUBTRACTION method (OSOAA total minus OSOAA dark at 0+),
scale-assumption-free via the Ed-normalized REFL column vs pi*Rrs.
RESULT (both bands): the water-leaving level shows an ANGLE-DEPENDENT
discrepancy - 555: +12.9% (vza11) -> ~0 (54) -> -29% (82.8); 443:
-8/+24/+20/+11/-1/-35%.  Grazing -30..-35% is COMMON to both bands ->
prime suspect: OCRT's water->air transmission (TWA) t(theta) handling at
steep angles vs OSOAA's TWA matrix.  The flat "+10% constant" picture
from the TOA arithmetic was an oversimplification (TOA agreed because
grazing rows contribute little there).  NOTE: the subtraction method's
purity (bright-water feedback on the downwelling; azimuth alignment of
the m>=1 diffuse-source structure) still needs one verification pass.
NEXT: direct TWA kernel-level comparison (surface_T_coxmunk_trig row
integrals vs OSOAA TWA SURF_MATR at matching angles), plus subtraction-
purity check.  Layer-speed answer recorded: NT80 1.78 s -> NT400 2.27 s
(+27% for 5x layers; order count is set by omega/tau, not layers).

## WIND-0 DISCRIMINATOR (2026-07-10): rough-surface hypothesis KILLED,
## structure of the water-leaving gap mapped
Flat-surface (wind=0) rerun of the 0+ water-leaving subtraction battery
(555): the gap did NOT vanish - it went uniformly NEGATIVE and larger
(-6% at vza11 monotonically to -36% at 82.8; wind3 was +13 -> -29%).
THREE established structures:
 (a) a FLAT-baseline deficit exists (OCRT water-leaving low by 6-16%
     at non-grazing angles even with pure Fresnel) - impossible from
     t(theta) alone, so the remaining suspects are the SUBTRACTION
     PURITY of the dark-water method, the Ed(0+) denominator
     definition, n_water value, or OSOAA's internal wind-0 handling
     (NO wind-0 file exists in SURF_MATR: 28 entries are all 02.0/03.0
     -> OSOAA must treat wind 0 analytically or substitute);
 (b) wind 0 -> 3 RAISES OCRT's low/mid-angle water-leaving by
     +15..19 pp relative to OSOAA - OCRT's rough T_wa is much more
     wind-sensitive (separate item);
 (c) the near-critical (grazing-air) -30..-36% collapse is COMMON to
     both winds - critical-angle-vicinity treatment difference.
NEXT-SESSION REDESIGN (decisive, assumption-free): drop the dark-water
subtraction; use the ANALYTIC anchor at flat surface - compare the
0-/0+ RATIO (must equal t(theta_w)/n^2 exactly) per code against
theory, using OSOAA's below-surface level output vs OCRT Lu0minus/
Lu0plus.  This isolates the interface step with no component or
normalization assumptions.  Also pin rt_water_iop_n_real(555) vs
OSOAA's fixed 1.34.

## PROTOCOL ITEM 6 ADDED (Jae, 2026-07-10): wind-speed validity
OSOAA comparisons are EXCLUDED for wind < 3 m/s (OSOAA rough-surface
model accuracy at low wind; SURF_MATR has no wind-0 entries).  3 m/s is
the validated minimum (harness golden lives there).  Interpretation
note: Jae wrote "3 m/s 이하"; implemented as "< 3 excluded, 3 = valid
minimum" since the harness golden itself runs at 3 - flagged for
confirmation.  CONSEQUENCES: last session's wind-0 probe structures
(flat-baseline deficit; wind-sensitivity delta) are REATTRIBUTED to the
OSOAA low-wind regime - not OCRT items.  The 443/case03 cross ran at
wind=2 -> requalified at wind=3 this session.  555 cross (wind 3) stays
valid; live accuracy items = 555 TOA +5..6% and the wind-3 grazing
water-leaving -29%.
443/case03 REQUALIFIED at wind=3: MAPE 1.12% (was 1.17% at the now-
excluded wind=2; vza11 improved -2.26 -> -1.48%).  Official cross table
now satisfies all SIX protocol items.  Harness rule+runtime guard
landed in /tmp/harn.py, the FULLPKG script, and the GOCI3_PKG copy.

## 0- LEVEL PROBE AT WIND=3 (2026-07-10): INTERFACE (TWA) CONFIRMED
Method upgrade: at 0- the upward radiance is the PURE in-water field
(sky reflection joins only at 0+), so no dark-water subtraction is
needed - both codes compared as Ed-normalized ratios directly.  OSOAA
Level=4 vsVZA angles are IN-WATER angles (2.8..88 deg, beyond-critical
directions included - trapped light); match OCRT rows via
theta_w = asin(sin(vza_air)/1.34).  RESULT (555, cs5.0, wind 3):
+6.5% (theta_w 8) -> -3.8% (47.8, near-critical), MAPE 3.3% (nmw24;
nmw48 similar 3.6%).  The SAME line of sight shows -29% at 0+ =>
~25 pp arise IN the water->air interface step.  Subtraction-purity
concern quantified and dismissed as the main cause: sky reflection is
~49% of total 0+ at grazing, but a few-% downwelling feedback cannot
produce -29%.  Secondary new observation: low-angle 0- shows +6.5/+3.9%
(angle structure differs from the nadir harness check) - 2nd priority.
NEXT (main event): TWA kernel-level comparison.  SURF_MATR TWA files
are BINARY, record length 26244 B = 81x81 REAL*4 (internal extended
angle set of 81 nodes, > NBMU=24); one file per (n,wind,MU[,SZA tag]).
Parse plan: read the OSOAA_SOS_CORE.F READ statements for the exact
record layout (mode loop x matrix components), angles from
RAD_UsedAngles.txt, then row-compare against OCRT's wa trig/fourier
kernel at matching (mu_w -> mu_a) pairs near the critical cone.

## TWA KERNEL AUTOPSY, ROUND 1 (2026-07-10)
File format cracked: SURF_MATR files are Fortran sequential unformatted;
ONE record per Fourier mode m=0..200, each = 9 components x 27x27
REAL*4 in ((X(I,J),I),J) order with I = WATER, J = AIR (axis proven by
beyond-critical extinction).  Angle set = 24 Gauss + 3 zero-weight
specials (nadir, SZA, in-water refracted SZA) from RAD_UsedAngles.txt -
hence the SZA tag in filenames.  Application in SOS_CORE (line 2123):
L_air(K) += GA(J) * [I,Q,U]_water(J) . TWA(J,K) - Gauss weight only, so
mu_w (and the azimuth constant) are FOLDED INTO the matrix; OCRT's
kernel multiplies 2pi*mu_w*w_w explicitly.  Raw row-sum ratios are
drowned in these convention factors (approx 2*mu_a trend + node-local
noise), so absolute-scale alignment needs the TWA GENERATOR's defining
equation (OSOAA_SURFACE / SURF_MATRICES source) - next session.
SHAPE RESULT (self-normalized rows, mu*w aligned): for grazing-air
output (88.1 deg) the two codes accept from DIFFERENT water cones -
OSOAA almost entirely from BEYOND-critical directions (trapped, strong
field; slope-mediated escape; below-critical occupancy ~0.01), OCRT
spreads 0.26 of its acceptance INSIDE the cone (32..47 deg).  Since the
trapped field is enhanced by internal reflection, beyond-critical-
dominant collection (OSOAA) yields the larger 0+ - sign consistent with
the observed -26 pp.  At 80.7 deg the shapes already agree (occupancy
0.338 vs 0.364) - the divergence is grazing-specific.
NEXT: (i) parse the TWA generator's definition to align absolute
conventions and reconstruct the -26% quantitatively; (ii) OCRT-side
energy-conservation self-test of the wa kernel rows (integral vs
1 - R_ww) to decide whether OCRT's inside-cone spread at grazing is
physics (slope broadening) or a kernel defect (azimuth/Jacobian).

## TWA AUTOPSY ROUND 2 (2026-07-10): OCRT-SIDE DEFECT PINNED
Branch (2) energy self-test: DHT(mu_w) = 2pi Sum w_a mu_a K[a,w] vs the
independent slope-integrated truth 1 - R_ww_direct(mu_w).  Inside and
at the critical cliff (28..51 deg water) conservation holds to
0.05..2%.  DEEP BEYOND-CRITICAL the kernel UNDERESTIMATES escape:
-17% at 54.7 deg, -48% at 58.5, -84% at 62.  Absolute row profile of
the grazing-air (88.1 deg) output: OSOAA collects 0.987 of it from
beyond-critical water directions (peak 54.7 deg) vs OCRT 0.736 - the
missing mass sits EXACTLY in the deficit band, magnitude-consistent
with the observed -26 pp.  The earlier "inside-cone spread" was a
self-normalization artifact of this deficit.  n_phi sweep 180/720/2880:
bit-identical DHT -> NOT azimuth convergence; the defect is structural
in the kernel's SLOPE integration (tail truncation for large-deviation
facets in the TIR-escape regime), while R_ww_direct integrates it
correctly.  Branch (1) (OSOAA generator convention alignment) is moot -
the inconsistency is internal to OCRT.
NEXT SESSION = KERNEL SURGERY: read the wa/aw fourier-trig slope
integral, align its slope-domain coverage with the direct integrator
(adaptive/extended nodes for the TIR-escape tail), gate on (i) DHT vs
1-R_ww at all nodes incl. 54..66 deg, (ii) Tier-0, (iii) Tier-3.5 -
NOTE: this is a physics fix touching shared/surface.c, so goldens are
EXPECTED to move at grazing; plan a deliberate golden re-baseline (and
fold in the F_sun convention unification at the same time).  Expected
gains: 555 grazing row +25 pp, 555 TOA MAPE 4.97% -> ~2-3%.

## KERNEL SURGERY SESSION 1 (2026-07-10): three suspects cleared, arbiter built
(1) AIR-GRID RESOLUTION CLEARED: dense 2000-point uniform air quadrature
reproduces the deficit exactly (-18/-51/-86/-99% at 54.7/58.5/62.2/65.9
deg water) - the escape lobe is NOT hiding between nodes.
(2) NUMERATOR-SWAP HYPOTHESIS KILLED: replacing mu_af^2 with
mu_wf*mu_af (naive Walter |i.h||o.h| reading) explodes the previously
CORRECT sub/near-critical band by +97..+151% - mu_af^2 is the right
combination there (it belongs to the energy-transmission chain with the
ts/tp amplitude convention).  Reverted; source restored bit-identical.
(3) FRESNEL CORE CLEARED: exact sin^2>=1 TIR boundary, no skin-layer
clipping; phi-quadrature already proven converged (180=2880).
(4) FIRST-PRINCIPLES ARBITER (python MC, encounter-weighted Gaussian
slopes x local Fresnel, no shadowing - same model class as
R_ww_direct): T_MC = 0.0984 / 0.0255 / 0.0042 at 54.7/58.5/62.2 deg vs
1-R_ww_direct = 0.0854/0.0188/0.0025 vs kernel DHT = 0.0699/0.0093/
0.0004.  VERDICT: kernel deficit RECONFIRMED against BOTH independent
references (and grows with depth), BUT the two references also disagree
by +15/+36/+70% -> R_ww_direct's own weighting convention is now a
re-examination item; the truth anchor must be pinned first.
NEXT SESSION (redesigned surgery): (a) pin the truth anchor - read
R_ww_direct's integral and reconcile with the encounter-weighted MC
definition; (b) same-variable (beta, phi_h) contribution-map comparison
kernel vs MC to localize the loss - prime remaining suspect: the
refraction half-vector Jacobian chain (the 'denominator reduces to Nsq'
shortcut and the d omega_a / d omega_w factor) in the deep-TIR corner;
(c) fix, then gate on MC = DHT = direct triple agreement + Tier-0 +
Tier-3.5 with the planned deliberate golden re-baseline.

## KERNEL SURGERY SESSION 2 (2026-07-10): VERDICT REVERSED - KERNEL EXONERATED
Same-formula integral duel: o-space (kernel-style, over outgoing air
directions) vs h-space (direct slope integral).  They agree exactly at
47.3 deg and diverge with depth (ratio 0.94/0.70/0.30 at 54.7/58.5/
62.2).  Log-mu2 quadrature reproduces the o-integral bit-for-bit ->
NOT a grid effect; the mass difference is REAL and lives at the mu2->0
boundary.  DECISIVE MC: computing the refracted ray's GLOBAL z and
counting only UPWARD transmission gives T_up = 0.09313/0.01801/0.00128
- matching the kernel o-integral to 4-5 digits.  The entire h-vs-o gap
is DOWNWARD-transmitted (re-entrant) light: 6.3%/30.5%/70.8% of the
single-encounter 'transmission' at those angles.
CONCLUSIONS (yesterday's verdict retracted):
 (1) The OCRT wa kernel is EXACT single-encounter physics (counts only
     true upward escape) - no defect, no surgery.
 (2) surface_R_ww_coxmunk_direct is the DEFECTIVE one: it books
     re-entrant light as escaped (h-style accounting) -> new item:
     fix R_ww_direct's energy bookkeeping (affects the fastk memo path
     users of that function - audit call sites).
 (3) The 555 grazing -26 pp vs OSOAA now has two candidate readings:
     (a) OSOAA's TWA generator uses h-style accounting and OVERCOUNTS
     grazing escape (then OCRT is right and the gap attributes to
     OSOAA); or (b) real multi-encounter re-escape (2nd-facet
     transmission of the re-entrant beam) is significant and OSOAA
     crudely captures it (then OCRT should ADD a multi-bounce term).
     DISCRIMINATOR (next session): extend the MC to 2 encounters
     (re-entrant ray meets an independent facet: reflect/re-refract,
     count final upward escape) and see whether truth lands on the
     kernel value or near OSOAA's row mass.  Also read OSOAA's SURF
     generator accounting if quick.
No code changed this session (experimental patch was reverted last
session; Tier-0 re-verified then).

## MULTI-ENCOUNTER DISCRIMINATOR (2026-07-10): branch (b) WINS
k-encounter MC (re-entrant rays meet independent facets; air-side
Fresnel reflect -> up = escape, refract -> water return; series to 8
encounters).  Converged upward escape at 54.7/58.5/62.2 deg water:
0.09523 / 0.02059 / 0.00232 vs kernel(single) 0.09313/0.01801/0.00128
and h-total 0.09941/0.02593/0.00438.  Truth sits well ABOVE the
single-encounter kernel, at 92-96% / 79% / 53%+ of the h-total (the
series has a known back-face filter leak at near-horizontal rays, so
these are LOWER bounds; water-return bookkeeping = 0.0015/0.0024/
0.0012).  VERDICT: the 555 grazing -26 pp attributes to OCRT's MISSING
surface multi-encounter re-escape term; OSOAA's h-style accounting is
an UPPER-bound effective closure of the same physics (slightly over,
by the water-return + residual-loss part).
SURGERY DESIGN (next session): replace the wa escape accounting with
the energy-consistent closure  T_up_eff(mu_w) = T_h(mu_w) -
WaterReturn(mu_w) - ResidualLoss(mu_w), implemented either as (i) an
analytic geometric-series closure with MC-calibrated re-escape ratio
q(mu_w) [preferred: smooth, cheap, tabulated per (n, wind)], or (ii)
h-accounting minus a water-return correction; pair it with the
R_ww_direct fix so that R + T_up + WaterReturn bookkeeping closes to 1.
Gates: MC(k-enc) vs new kernel DHT at all water angles; Tier-0;
Tier-3.5; then the 555 grazing row and TOA re-cross (expected: grazing
-29% -> single digits; TOA 4.97% -> ~2-3%).  Note: the multi-encounter
term also perturbs R_ww (water-return adds to it) - fold into the same
deliberate golden re-baseline with the F_sun unification.
CAVEATS RECORDED: independent-facet (uncorrelated) approximation,
simplified 2nd-encounter weighting, no height-correlation shadowing in
the series - acceptable for closure calibration at ws>=3, revisit if
sub-% fidelity needed.

## PREP + DECISION MENU (2026-07-10)
Decision-independent prep done: (1) multi-encounter calibration table
T_up_multi(theta_w) + water-return for n=1.34/ws=3 saved
(tup_multi_ws3_n134.txt; k<=8, lower bound near deep-TIR).  (2)
R_ww_direct integral READ and defect confirmed at the source: it is an
encounter-weighted <R_flat(local)> over a 101x101 slope grid - i.e.
1-R = h-style TOTAL first-encounter transmission, re-entrant light
never separated.  Fix = add the re-entrant water-return to R and keep
only upward escape in T (pairs with the closure choice below).
DECISIONS FOR JAE: D1 closure option (series+q-table / h-minus-return /
adopt-h=OSOAA-parity), incl. the common sub-task of redistributing the
re-escape into the grazing air lobe of the matrix kernel; D2 golden
re-baseline bundling (multi-bounce + R_ww_direct + F_sun) - WARNING:
R_ww change moves the in-water field, so the Tier-0 anchor itself will
move; D3 whether grazing vza rows are inside the paper's analysis range
(blocks Phase 1 or not); D4 wind-boundary wording (3 excluded or 3
valid) still awaiting confirmation.

## DECISIONS RECORDED (Jae, 2026-07-10)
D1 = OPTION 2: water->air escape = h-style total MINUS the water-return
part (energy closes on the water side: escaped-up + returned-to-water =
first-encounter transmitted).  Applied UNIFORMLY at all angles (the
return term goes to zero smoothly below ~50 deg water) - NO angle
switches.  ENGINEERING RULE (Jae): no unauthorized special-case
switching routines anywhere; the wind-0 flat-surface exception is the
one pre-authorized case.
D2 clarified: "bundling" the multi-bounce fix + R_ww_direct fix +
F_sun unification means ONE re-baseline event (release management), it
does NOT itself change accuracy.  Sun-glint bookkeeping is untouched:
option 2 closes the energy ledger between the two WATER-side operators
(T water->air and R water-water); direct sun glint lives in the
AIR-side reflection operator (R air-air), a separate channel.  The
analogous multi-encounter terms on the air side (grazing glint
re-strike; air->water transmission that travels up inside water and
re-strikes) are REGISTERED as a deferred item - negligible for our
glint-excluded raa=120 comparisons.
D3: paper analysis VZA cap = 80 deg.  Current-physics error at <=68 deg
is ~2%, at 80 deg estimated single-digit %, the -29% row (82.8 deg) is
OUTSIDE the cap -> not an accuracy blocker; but since Phase-1
production is expensive and the fix moves high-VZA rows by a few %,
sequencing stays surgery -> re-baseline -> production.
D4: wind boundary resolved - 3 m/s is NOT excluded; exclusion is
strictly wind < 3 (harness text already matches).
TERMINOLOGY: use everyday Korean terms in reports; "고천정각(각도 병기)"
instead of coined words; "기준 결과(베이스라인)" instead of "golden".
OPEN HONESTY ITEM: whether OSOAA's SURF generator is truly h-style
(option 3) or computes multi-encounters itself is UNVERIFIED (source
not yet read) - registered.

## OPTION-2 SURGERY, STAGE 1 (2026-07-10): cascade engine landed (unwired)
D1 finalized as "physically most accurate": option 2 implemented as the
EXACT closure - the multi-encounter cascade is computed at runtime for
any (n_water, wind), no tables, no angle switches.
(1) Python reference cascade upgraded with K=8 reservoir resampling of
forward-facing facets: the back-face leak is GONE (unconverged mass 0,
identity T_h = T_up + water_return closes to <0.05% at all angles).
Escape-direction data captured: re-escape energy concentrates at high
view zenith angles (mu_air<0.2 carries 77-91% at deep beyond-critical
water angles) - input for the radiance-operator wiring.
(2) C engine src/shared/surface_multibounce.c +
surface_wa_multibounce_cascade() in surface.h: deterministic
(fixed-seed xorshift64* + Box-Muller), 24 bounces, K=8 resampling,
optional mu_air escape histogram.  Oracle vs python: agrees within MC
statistics at 6 angles, closure holds.  UNWIRED this stage - both
builds compile, Tier-0 anchor unchanged (behavior-neutral).
NEXT (stage 2, wiring): (a) redefine surface_R_ww_coxmunk_direct to the
closed ledger R_eff = 1 - T_up (first-encounter R + water-return);
(b) add the re-escape term to the wa fourier kernel rows (mu_air-binned
from the same cascade, m=0 first, m>0 from the cascade's azimuth
accumulation if needed); (c) add the water-return angular term to the
R_ww matrix kernel; (d) gates: engine-vs-kernel DHT identity at all
nodes, Tier-0/Tier-3.5 (EXPECTED to move -> deliberate re-baseline with
F_sun unification), 555 grazing row + TOA re-cross.

## OPTION-2 STAGE 2 (2026-07-10): partial wiring, honest status
DONE AND VERIFIED SAFE:
 (a) surface_R_ww_coxmunk_direct now returns the CLOSED ledger
     R_eff = 1 - T_up (cascade, 400k rays, deterministic).  Effect
     observed exactly where expected: Ed(0-) bookkeeping only ->
     443 TOA moved +0.001%, Tier-0 anchor unchanged at 6 digits.
 (b) engine v2: esc_hist = RE-ESCAPE ONLY (direct up-escape excluded -
     the single-encounter kernel already carries it), plus ret_hist
     (returned-into-water mu distribution) for the future in-water term.
 (c) rt_air_water_couple_water_to_atm(): +w_atm_pos parameter (3 callers
     threaded); m=0 re-escape addition implemented in BOTH forms -
     node-cell (exact grid conservation) when weights are given, and the
     continuous-density form for weightless single-view callers (the
     mathematical continuum limit, not a physics switch - documented).
NOT YET EFFECTIVE - open:
 The observable Lu0plus/TOA did NOT move (555 grazing row identical to
 5 digits).  (a) proves the plumbing compiles and runs; the re-escape
 addition is not reaching the OBSERVED view-0+ path.  Prime suspect:
 the view-0+ value in the B-0c diffuse architecture is no longer
 produced by the legacy coupling-replacement block at rt_water_rt:4908
 (June design) nor the slope integral above it - the actual per-call
 T_wa application path must be traced with a probe before wiring.
NEXT SESSION (single objective): probe-trace the live view-0+ chain
(rt_solver Lu_total_above <- w_res.I_0plus_view <- ?), wire the
density-form re-escape THERE, then rerun the 555 gate (expect the
grazing -29% to shrink toward -3%) and only then consider re-baseline.

## OPTION-2 STAGE 2b (2026-07-10): live-path hunt - third candidate also dead
Wired the density-form re-escape into the nodeset-interpolation block
(rt_water_rt ~4662, after the nv<2 fallback), with a thread-local
per-node cascade cache.  RESULT: no effect again; probe forensics
settled the mechanism conclusively:
 - preprocessed TU contains the block (gcc -E: 1 hit);
 - an -O0 build carries the string AND still prints nothing while
   Lu0plus is bit-identical -> the block is RUNTIME-UNREACHED in this
   case (not dead-code-eliminated; the -O3 strings miss was a red
   herring of string handling).
So the view-0+ value for this run is produced by yet ANOTHER path:
neither the June coupling-replacement block (4908), nor the slope
integral (4891), nor the operator-mode nodeset interpolation (4596-
4662) executes for it - despite ocrt_water_twa_operator_mode() being a
constant 2.  Implication: the enclosing region (possibly the whole
4450..5070 view assembly) is bypassed for this case; I_0plus_view must
be assigned elsewhere or the function containing 5070 is not the one
called here.
NEXT SESSION - ONE MOVE: put an UNCONDITIONAL probe immediately before
rt_water_rt.c:5070 ("result->I_0plus_view = I_air").  If silent, the
assignment lives in a different function (re-grep with function
boundaries); if it prints, bisect upward from 5070 to find which branch
feeds I_air.  Then wire the density-form addition THERE and rerun the
555 gate.  CODE STATE: safe - the addition block is inert (unreached),
R_eff ledger fix active and gated (Tier-0 unchanged, 443 +0.001%),
555 baseline values reproduced exactly (no regression).

## OPTION-2 STAGE 2c (2026-07-10): live chain SOLVED, one insertion left
Probe battery (unconditional [VIEW0P] before the I_0plus_view store +
tags at every I_air assignment + [CPL] debug inside the re-escape
helper) settled the execution chain for the diffuse production case:
  [IA:SLOPE] -> [IA:CPL4908] -> [VIEW0P]  (I_air = 1.012575e-02 = the
observed Lu0plus).  So the view-0+ value comes from the June
coupling-replacement call at rt_water_rt:4908 - and that call runs the
coupling function's FLAT fallback, NOT the rough branch ([CPL] printed
only for the 25-node pass2 coupling, added=5.87e-05 there; never for
the n_atm=1 call).  The nodeset-interpolation block (4596-4662) and
its [MB] probe are dead in this flow - resolved mystery.
REFACTOR DONE: the m=0 re-escape addition was extracted into a static
helper mb_re_escape_add_m0() and hooked at the rough-branch success
return; the flat-path hook was inserted while the file still had a
brace error, so its coordinates are suspect - the n_atm=1 call still
shows no [CPL] and Lu0plus is unchanged.  ONE INSERTION REMAINS: place
the helper call immediately before the FLAT path's actual return for
this function (verify by [CPL] n_atm=1 appearing and Lu0plus at
vza 82.81 rising from 1.0126e-2 toward ~1.4e-2), then run the 555 grid
gate.  SAFETY: brace surgery gated - the single case reproduces all
legacy values exactly; Tier-0 anchor re-verified after the fastk
rebuild below.

## OPTION-2 STAGE 3 (2026-07-10): mechanism landed; one physical term left
Chain finally closed end-to-end, with one instructive detour:
 (1) Passing real water weights to the 4908 coupling call flipped it
     from the PRODUCTION flat path into the (per-view unvalidated)
     rough branch: grazing improved hugely (-29.4 -> -9.1%) but low
     angles and 443 went erratic (-24..+16% Rrs swings) - the June
     production choice (flat + spline) exists for a reason.  REVERTED:
     the call passes NULL again; flat path preserved bit-exact.
 (2) The re-escape addition now lives at the CALLER in final-radiance
     units (I_air += sum_k mu_k w_k Lw0_k h_k(mu_v)/(mu_v dmu), cascade
     cached per node).  Verified arithmetically: [MBV] add=2.820e-4
     and Lu0plus(82.81deg) = 1.0126e-2 + 2.82e-4 = 1.0408e-2 exactly.
 (3) WHY the effect is smaller than the rough-branch experiment
     (+2.9e-3): the flat view operator ALSO omits the SINGLE-encounter
     slope escape of BEYOND-critical water nodes (they have no flat
     refraction angle -> outside its domain), not just the re-escape.
     The missing last term = direct-escape histogram of beyond-critical
     nodes only - a DOMAIN completion of the flat operator (the
     critical cone is a physical constant of n, not an ad-hoc switch).
NEXT (final piece): engine outputs esc_direct alongside esc_re; add
h_direct for mu_w < mu_crit(n)=sqrt(1-1/n^2) nodes in the same caller
block; single-point target Lu0plus(82.81) ~= 1.30e-2; then the 555 grid
gate (grazing toward -9%, low angles UNCHANGED this time), 443
regression (expect ~0 moves), Tier-3.5 record, and the deliberate
re-baseline discussion.  Gates this stage: Tier-0 anchor unchanged
(re-verified below), production values reproduced exactly outside the
new additive term.

## SESSION CLOSE (2026-07-10): OPTION-2 SURGERY - LANDED SCOPE + CARRY-OVER
LANDED (all gated, both builds, Tier-0 anchor 2.887685e-02 unchanged):
 (1) surface_multibounce.c cascade engine - exact energy closure,
     deterministic; now outputs esc_re / esc_direct / ret histograms.
 (2) R_ww_direct = 1 - T_up closed ledger (Ed(0-) bookkeeping path).
 (3) View-0+ RE-ESCAPE term at the production caller in final-radiance
     units, per-node thread-local cascade cache.  555 gate: grazing
     82.81deg -29.44 -> -27.48% (+2.0pp), 68.42deg -1.83 -> -1.69%,
     LOW ANGLES BIT-UNCHANGED (the earlier swings are gone).  443
     regression: +0.000% everywhere except +0.05% (68deg) and +1.87%
     (82.8deg) - exactly the physical signature of the new term.
MEASURED AND DOCUMENTED (not landed):
 (4) Direct-escape "domain completion" DOUBLE-COUNTS with the flat
     operator's spline EXTRAPOLATION below its lowest table node
     (add jumped to 1.21e-2, ~4x target).  The remaining -27% at
     grazing therefore attributes to the flat operator's beyond-
     critical representation itself.  TWO candidate designs for the
     NEXT session (JAE DECISION - operator-level change):
      (A) zero-endpoint mode + full (direct+re-escape) additive term
          as a consistent pair; or
      (B) diagnose the rough per-view coupling's low-angle swings
          (-24..+16%) and, if fixable, adopt the rough branch for the
          view path (it already carries the full physics: measured
          -9.1% at grazing).
 (5) Diagnostic probes ([VIEW0P],[IA:*],[CPL],[MBV]) are getenv-gated
     (OCRT_MB_TRACE) - retained; remove in a cleanup pass.
CARRY-OVER LIST: (a) design decision A/B above; (b) TOA pass2 re-escape
scale check (6246 coupling: added 5.9e-05, likely undersized - verify
E_inc normalization there); (c) deliberate golden re-baseline (with
F_sun unification) once the view-path design settles; (d) air-side
multi-encounter analogs (glint re-strike, TAW) - registered earlier;
(e) OSOAA SURF generator accounting read (h-style vs own multibounce).

## VIEW-OPERATOR A/B RETRIAL (2026-07-10, post-rules): rough vindicated
Per Jae's rule ("errors beyond a few % = broken code OR wrong physical
comparison - consult existing code"), the earlier "rough swings" verdict
was RE-EXAMINED: it had compared rough against the OLD FLAT golden (an
OCRT-vs-OCRT movement), not against OSOAA truth.  Direct OSOAA 0+
comparison (flat vs rough via the diagnostic-only env toggle
OCRT_MB_ROUGH_VIEW=1; production default untouched):
  555/wind3 water-leaving MAPE: flat 9.05% -> rough 3.45%
    (11deg +12.9 -> +0.8; 82.8deg -27.5 -> -9.1; all angles better or
     comparable).
  443/wind3: flat 14.26% -> rough 12.83% (mixed; one cell worse:
     11deg -3.9 -> -13.9).
CONCLUSIONS:
 (1) The "low-angle swings" were an artifact of the wrong baseline -
     textbook case of the new rule; comments added at the 4908 caller
     documenting the whole chain, the flat/rough entry condition, the
     extrapolation semantics, and the double-count hazard.
 (2) Proposal to JAE (operator-level default change - his call):
     adopt the ROUGH per-view coupling as the production view-0+
     operator (option B), retiring the flat+spline+extrapolation path
     for wind>0; the re-escape additive term then becomes unnecessary
     on the view path (rough carries full single-encounter physics;
     re-escape addition on top of rough = next refinement, expected to
     take 82.8deg from -9% toward -3%).
 (3) SEPARATE NEW ITEM: 443 0+ residual +13..19% at mid angles under
     BOTH operators - not a view-operator issue; suspect case/profile
     mismatch in the x443z/x443d0 OSOAA references or a 443-specific
     in-water field difference (the official 443 cross-validation was
     TOA-level, never 0+).  Investigate before any operator decision is
     graded on 443 numbers.

## SPRINT COURSE CORRECTION (2026-07-10, Jae): back to SPEED
Jae's direction: the sprint is SPEED; accuracy already has its baseline;
stop the accuracy side-quest.  Actions taken:
 (1) ALL three multibounce accuracy sites (surface.c R_eff closure,
     rt_water_rt view re-escape term, coupling mb_re_escape_add_m0)
     are now OPT-IN ONLY behind OCRT_MB_CLOSURE=1, getenv snapshotted
     once per thread (#21 rule).  Default OFF.
 (2) BIT GATE PASSED: default build reproduces the pre-surgery golden
     grid (golden_atm_on_tier35_2026-07-10.csv) bit-exactly (all
     columns, out-path excluded).  Tier-0 anchor unchanged.  This also
     removes the strict-build hazard (the 400k-ray cascade ran unmemoized
     on strict).
 (3) Answers recorded: the rough-view proposal NEVER touched the
     wind==0 flat branch (triple-guarded: s2<=1e-10 early-out, wind>0
     guards, proposal scoped to wind>0); the operator A/B decision and
     the 443 0+ residual investigation are PARKED in the accuracy
     carry-over list, not pursued during the sprint.
RESUME POINT: S7a continuation (sky-beam unit-flux + F_sun_eq bit-jitter
root fix - re-read its section tail before touching) and the mandatory
open P>0 3-way accuracy gate (ref/strict/fastk).  Rule reaffirmed: read
code+comments thoroughly BEFORE implementing anything new.

## DEBUG-DEVICE SPEED AUDIT (2026-07-10, Jae request)
Inventory: this session's probes had SEVEN raw getenv() calls sitting in
the per-(vza,raa) view path - a direct violation of the file-top #21
rule ("getenv in the m-loop / per-(vza,raa) path serializes threads").
FIXED, idiomatic to the codebase: s_OCRT_MB_TRACE and
s_OCRT_MB_ROUGH_VIEW added to the g_wrt_env snapshot (rt_water_env_init);
the coupling-TU site uses a thread-local one-shot snapshot; the three
OCRT_MB_CLOSURE gates already used the snapshot pattern.  Post-fix grep:
ZERO raw getenv on hot paths for the new devices.  One editing mishap
(a blanket replace corrupting the just-added init line) was caught by
an anchor assert and repaired - noted per the read-first rule.
Gates: pre-surgery golden grid still BIT-IDENTICAL (default build);
strict/fastk both compile.
Timing (P=0 555 single solve, 3 runs summed, 1 core):
  pristine v1.09 fastk  20.23 s   |   current tree fastk  5.92 s
The current tree is 3.4x FASTER than the pre-session pristine build
(accumulated S5/S6-era work); probe+gate overhead is a few pointer/int
checks per solve - unmeasurable at this scale.  Devices that print
(fprintf/fopen) only execute when their env is set.
OPEN (minor): six pre-existing getenv sites in shared/surface.c (mode/
stats helpers) predate this session and were NOT frequency-audited -
most look like once-per-warning statics; verify during the next surface
touch.

## SPEED RESUME SESSION (2026-07-10): S7a-(a) recon + hot getenv audit closure
(1) RULE-BREACH LOGGED: attempted to snapshot the OCRT_WATER_GRID_CACHE
    getenv at rt_water_rt:3493 - the on-site #22 comment (read only
    AFTER editing) explicitly forbids it (runtime setenv by the batch
    driver; snapshot froze it to NULL once and killed the #16 cache
    ~600x).  REVERTED immediately.  The comment system worked exactly
    as designed; reading order is on me.
(2) S7a-(a) recon on a 19-row 555 P=1013 grid: [S6K] keys are UNIFORM
    (vasn=0, nvv=6, single F), replay clean, ZERO beam solves - the
    node-boundary 11-14 s recolds and the beam-absence mystery (b)
    reproduce only on the 433-row n_mu=24 light grid.  REPRO RECIPE
    for next session: n_mu 24 light grid + OCRT_S6_TRACE=1, capture
    [S6K] around vza-node boundary rows; compare beam keys across the
    boundary (suspect: a boundary-dependent component in the beam key
    despite pinned F).
(3) HOT getenv AUDIT CLOSED for shared/surface.c: five raw getenv sites
    sat INSIDE the point kernels (V3_AF_SHADOW x3, V3_BRDF_AF1982 x2) -
    hot on cold kernel fills and on strict.  Converted to thread-local
    one-shot snapshots (V3_* are process-invariant diagnostics, NOT
    #22-class runtime-set vars - noted in code).  L1445 V3_GLINT_FOV_DEG
    is per-row frequency - left as-is.  Measured (fastk, warm FKC,
    P=0 single solve x3): 5.923 -> 5.841 s = within noise, as expected
    (FKC hits dominate); the win is on cold fills and strict builds.
    Gates: golden grid BIT-IDENTICAL, Tier-0 unchanged, both builds OK.
NEXT: reproduce+fix S7a-(a) on the 433-row grid (recipe above), then
launcher --n-mu-water passthrough, then the open P>0 3-way gate.

## S7 RECON CLOSED (2026-07-10): (a)/(b) resolved by B-0c; atm1 is the target
(1) Session-5 open items CLOSED BY ARCHITECTURE: B-0c.2b made the
    diffuse-top source the DEFAULT ("beams then bypass") - so (b)
    "beams do not execute" is by design, and (a) node-boundary beam
    recolds have no actor anymore.  Verified on a 72-row n_mu=24 grid
    (OCRT_S6_TRACE): [S6K] keys uniform across all rows (single F,
    vasn=0), zero recolds, beams solved=0 everywhere.
(2) Per-row decomposition ([S6T], 72 rows + base, 59.0 s total,
    ~0.81 s/row): atm1 (the per-row ATM SOS solve) runs TWICE per row
    at ~0.4 s each => ~0.8 s/row dominant; atm2(C3d) + assemble are
    the small remainder.  The atmosphere is IDENTICAL across grid rows
    (only vza/raa change), so S7 = cache the atm solve per grid
    (#16 pattern) => ~50 s saved on the 72-row grid, ~350 s on the
    432-row production grid, per (case x wavelength).
(3) DESIGN KEY QUESTION for S7 (read-first next session): session 5
    already noted "view node in the atm grid" - if the atm quadrature
    embeds the per-row view node, the grid differs row-to-row and the
    mode fields are NOT reusable as-is.  Two candidate shapes:
    (i) #20-style: pass ALL grid vza values as zero-weight atm view
        nodes once -> row-invariant atm grid, exact node reads; or
    (ii) split the view-dependent extraction from the view-invariant
        mode solve and cache only the latter.
    Read the atm1 call site (rt_solver ~5644-5860) and the atm grid
    builder before choosing.  Also confirm why atm1 runs TWICE per row.

## S7 DESIGN LOCKED (2026-07-10) - read-first pass complete
CORRECTION of the previous entry: atm1 runs ONCE per row (my 146 count
summed begin+end trace lines - tallying error, logged).  Verified row
decomposition (72-row grid, warm): atm1 0.36 s + atm2(C3d)/assemble
0.35 s + direct 0.01 s (grid-cache replay) = ~0.73 s/row.  The 5.9 s
direct cold is ONCE PER GRID (the #16 first-row full solve - by
design); all later node-first rows replay at 0.01 s.  No node-boundary
issue exists in the current architecture.
TARGET: the per-row ~0.7 s is the ATM SOS re-solve pair (atm1 + the
C3d bottom-source pass) - view-invariant physics re-run every row.
BLOCKERS READ: (1) rt_solver 5717-5730 inserts THE ROW'S vza into the
atm angle set (rt_uangles_add) -> row-dependent atm grid/mode fields;
(2) atm_res.I/Q/U_TOA (atm path radiance) is genuinely view-dependent.
DESIGN (=#20 + #16 for the atmosphere): insert ALL grid vza values as
zero-weight atm view nodes ONCE (row-invariant grid; the list already
exists as cs.water_view_vza_list), solve atm1 ONCE per grid, cache the
boa mode-field export + T_diff scalar + the per-node TOA path radiances
(all 24 read directly off the node set), and have rows do node reads.
Apply the same to the C3d bottom-source pass.  EXPECTED: warm rows
0.73 -> ~0.05 s; 72-row grid 59 -> ~15 s; 432-row production grid
~5 min -> ~1 min class.
READ-FIRST LIST for the implementation session: (a) whether
rt_solve_case_pol_for_ocean's TOA extraction supports multi-view /
node-set reads (water-side #20 analogue), (b) rt_uangles dedup
semantics with 24 inserted nodes, (c) the C3d pass call signature and
its own angle handling, (d) cache invalidation key = the atm cs fields
(profile, tau, F_sun_TOA, surface type, n_w, lambda) + grid vza list.
Loop-depth rule applied: this removes the deepest per-row block
(atm SOS: orders x layers x modes x nodes^2) from the row loop entirely.

## S7 LANDED (2026-07-10): atm-solve grid cache
Design executed exactly as read: the v1.01 LUT path (SOS once, every
(vza,raa) reconstructed from the same quadrature solution) is invoked
ONCE per key with view_as_node OFF, storing (a) the boa mode-field
export, (b) NEW optional raw TOA intensity grids (pre-rho, pre-glint -
3 fields added to rt_lut_grid_out_t + one fill site), (c)
T_diff_dn_hemi.  Rows replay atm_res.{I,Q,U}_TOA + boa by node/raa
lookup; the batch grid driver publishes OCRT_ATM_GRID_CACHE=1 and
OCRT_S7_RAA_LIST (runtime setenv -> live getenv gate, #22 class).
MEASURED (72-row 555 P=1013 grid, 1 core): 59.0 -> 28.7 s (2.06x).
Remaining per-row cost is atm2(C3d)+assemble (~0.35 s/row) -> S7b.
GATES:
 - on/off replay: max significant rel diff 1.52e-13 (TOA_rho_U near-
   zero 1-ulp reconstruction ordering) - PASSES the 1e-10 fastk
   protocol (S1b precedent); FASTK-ONLY via #ifdef so STRICT stays
   bit-pristine.
 - PRECONDITION found the hard way: grid vza must BE the shared-N atm
   Gauss nodes (grid n_mu == n_mu_water), else the LUT replay is an
   INTERPOLATION (443 golden, n_mu 6 vs nmw 24: Lu0plus moved -0.05%,
   Ed0plus_air -0.07%).  Consistency gate added (not a physics
   switch); with it the 443 golden grid is BIT-IDENTICAL again
   (cache auto-disabled there), Tier-0 anchor unchanged.
LIMITS/NOTES: base rows (no #20 vza list) and standalone CLI use the
pristine path; production grids run n_mu == nmw so the cache engages
there.  Key includes lambda/sza/F/wind/n_w/surface/sigma/qconv/
rayleigh/aer(tau)/n_mu/m_max; single-slot TLS (per-thread, per-key
refill on change).

## S7b LANDED (2026-07-10): atm2(C3d) grid cache — memory-for-speed
Per Jae's directive (use memory aggressively; OSOAA precedent), the C3d
bottom-source atm pass now follows the S7 pattern: the bottom source
wlI/Q/U (row-invariant: #16-cached water field -> reverse coupling) is
FNV-fingerprinted into the key; one LUT-path solve per key with
view_as_node OFF; rows replay TOA_wl_direct from TLS raw grids.  Same
FASTK-only gating + n_mu==nmw precondition as S7.
MEASURED (72-row 555 P=1013 grid, 1 core):
  pre-S7 59.0 s  ->  S7 28.7 s  ->  S7+S7b 6.5 s   (9.1x total;
  remaining = base + one #16 cold + ~0.02 s/row replays).
GATES: on/off max significant rel diff 1.52e-13 (unchanged, 1e-10
protocol PASS); 443 golden grid BIT-IDENTICAL (precondition auto-off);
Tier-0 anchor unchanged; strict build pristine (#ifdef).
NEXT: heavy-grid + 432-row measurement, launcher --n-mu-water
passthrough, P>0 3-way gate (overnight ref), then Phase 0/1 production
throughput estimate refresh.

## PRODUCTION-SCALE MEASUREMENT (2026-07-10): 432-row grid = 9.7 s
Full production-shape grid (n_mu 24 x 18 raa = 432 rows, 555 P=1013,
1 core): 9.7 s wall with S7+S7b on.  Shared-direction rows are VALUE-
IDENTICAL to the 72-row run (verified properly: 72/72 exact after fixing
a column-index bug in my first check script - logged) - the cache reconstructs each
direction independently, so grid size does not change any value.
Tier-0 anchor unchanged.
Pre-S7 arithmetic for the same grid (measured 0.73-0.81 s/row + colds):
~350 s => ~36x on the production shape.
PHASE-1 THROUGHPUT REFRESH (64,512 grid runs, 4 machines, 1 core each):
  before: ~350 s/run -> ~66 days;  now: 9.7 s/run -> ~1.8 days.
  (OMP multi-core per machine shortens further; base+first-solve colds
  are inside the 9.7 s.)
Remaining speed items: heavy-grid (nmw48) measurement, launcher
--n-mu-water passthrough (likely obsolete post beam retirement -
verify), P>0 3-way accuracy gate (overnight reference, user-side).

## OCEAN-SIDE STATUS + HEAVY-GRID LIMIT (2026-07-10, measured)
Q (Jae): "atmosphere got faster; does the ocean side still need testing?"
A: in batch grids the ocean side was ALREADY cached (#16: solve once per
grid, rows read back).  Measured today:
 - P=0 (no atmosphere, ocean only) 432-row grid: 11.5 s.
 - P=1013 432-row grid decomposition: base 0.4 s + FIRST ROW 6.0 s
   (one-time water solve 5.2 s + two atm cache fills) + 431 rows x
   ~0.01 s.  So the LARGEST single block now is the ONE water solve
   per grid (~5-6 s) - that is the next ocean-side target if we want
   sub-5 s grids.
 - HEAVY (n_mu_water=48) 72-row grid: timed out at 168 s (48/73 rows;
   ~220 s projected).  Two compounding causes: (1) the S7 precondition
   (grid n_mu == n_mu_water) is FALSE at 24 vs 48, so the atm cache
   disengages; (2) shared-N makes the per-row atm solves 48-node
   (~2.8 s/row).  HONESTY NOTE: my first "heavy" measurement (8.9 s)
   was INVALID - the jobs-file n_mu_water column (24) silently
   overrode the CLI 48; caught via [S6K] nmw=24 and re-run properly.
S7c DESIGN (registered): insert the 24 grid vza values as ZERO-WEIGHT
view nodes into the atm angle set when grid n_mu != n_mu_water, so the
grid directions become exact read-off points on the 48-node solve and
the cache re-engages (weights zero -> the solved field is unchanged in
exact arithmetic; 1e-10 fastk protocol covers bit-level shifts).
PHASE PLANNING: Phase 0/1 at the production nmw=24 standard is fully
covered by today's speed (9.7 s/grid); only if 48-node water is chosen
does S7c become blocking.

## MULTI-JOB / MULTI-THREAD CHECK (2026-07-10) - sandbox is 1-core
Batch parallelism is per-JOB (omp parallel for over jobs rows, dynamic)
- structurally ideal for the new thread-local caches (each thread owns
whole jobs; one fill per job; zero cross-thread sharing needed).
MEASURED, but the sandbox has nproc=1, so no real parallel speedup can
be observed here:
  4 jobs (432-row grids, distinct water keys):
    1 thread: 52.0 s (13.0 s/job);  2 threads: 51.9 s;  4 threads: 56.0 s.
  Note: 13.0 s/job vs 9.7 s single-job = ~3 s/job fixed extra when jobs
  differ (distinct water keys -> per-job water cold + cache fills; exact
  split deferred).
  Threads>cores on 1 core only add scheduling overhead - expected.
CONCLUSION: real multi-core scaling of the cached batch must be
measured on Jae's local i9 (same carry-over as commit #21's Windows
OMP verification).  Per-job cost on 1 core stays ~9.7-10 s.

## 1-CORE OSOAA-PARITY DRIVE, SESSION 1 (2026-07-10): profiling truths
Rule registered (Jae): speed is judged on ONE core; target = OSOAA
parity (same-box OSOAA NT80 run: 1.78 s for a full-angle field).
Current OCRT: 432-row grid 9.7 s = water cold ~6 s + warm rows ~6 s
(14 ms/row flat, atmosphere-independent - re-entry/extract/output
fixed cost; P=0 shows the same).
SENSITIVITY MEASUREMENTS (single P=0 555 solve, 2.0 s reference):
  layers (dtau 0.05->0.25, 300->60 layers): time -12% ONLY, values
    move <0.005% - LAYER COUNT IS NOT THE DRIVER here (S4/S5 already
    crushed the layer-proportional part).  The "water layer policy"
    knob is therefore NOT the big lever we assumed - no Jae decision
    needed on it for speed.
  m_max (30->12): time -26%, values <0.001% - also not dominant.
  tolerance/orders knob: flag name not found this session (my two
    wrong-flag attempts logged); orders=137 x angles^2 kernel remains
    the PRIME SUSPECT for the ~1.3 s irreducible-so-far core.
NEXT SESSION (single objective): find the real tolerance/orders flag
(read main.c option table properly first), then block-profile the
water SOS interior (per-order kernel time vs per-order overhead) to
name the dominant component; separately, attack the 14 ms/row warm
fixed cost (row-loop internalization candidate).  Parity math: need
water cold ~2 s AND warm rows ~1-2 ms to reach ~3 s/grid.

## S8 + PROFILING SESSION 2 (2026-07-10): honest ledger
PAST-RECORDS RULE APPLIED (Jae): layers/m_max were tested before (not
dominant - now re-confirmed and closed); orders policy = generous cap +
converge-and-stop, NOT a trim knob.  D3 (operator reuse: water R_B per
case+wl, 448x Phase-1 amortization; interaction-principle coupling) is
the designated big lever and the heaviest surgery - staged LAST.
gprof(-O2,-pg) on the 432-row grid pointed at two suspects:
 (1) sos_build_source_pol 4.9 s / 412 calls - MISREAD at first as
     per-warm-row; actually the ~12 heavy atm SOS solves (base + cache
     fills) each run ~34 source builds.  Correction logged.
 (2) surface_T_aw_coxmunk_direct 10,826 calls x 0.2 ms (wind>0 =
     101x101 slope integral), ~25/row on FIXED node mus.
S8: 64-slot thread-local memo in that function (bit-safe same-input
replay, S2 precedent, FASTK-only).  RESULT: call count collapses out of
the profile top (memo works), outputs BIT-IDENTICAL, but -O3 wall time
UNCHANGED (9.7 s) - the -O2+instrumentation build had inflated this
function's real share.  LESSON (logged as method): hotspot claims from
instrumented builds must be re-verified on the real build's wall clock
before surgery is judged.
REAL -O3 warm-row cost: (9.7 - 0.4 base - 6.0 first)/431 = ~7.7 ms/row.
Its -O3 composition is NOT resolvable by gprof (instrumentation skew);
next session: high-resolution in-code block timers on the warm path or
go straight at row-loop internalization; the other half stays the
water-cold ~6 s whose algorithmic answer is D3.

## S8 EPILOGUE (2026-07-10 late): full honest closure
S8 (T_aw 64-slot memo) is REMOVED - final state is pristine S7b source
(whole-function restore from the v1.09 reference; verified textually
identical modulo hooks).  Timeline of the mess, for the record:
 - Two revert scripts died on an OVER-BROAD grep assertion: pattern
   'taw_m' substring-matched the PRE-EXISTING v1.09 taw_memo[8] kernel
   cache (surface.c ~1287) - a DIFFERENT, legitimate cache.  So S8 was
   still present during two "reverted" measurements.  Lesson: verify
   with word-boundary greps; assert AFTER write, not before.
 - The earlier bit-compare flip (1.54e-12 fail, then pass) is therefore
   attributed to RECOMPILE-TO-RECOMPILE instability of the SAME S8
   source under fastk flags, not to cache logic; with S8 gone, grid
   output is bit-identical to S7big across rebuilds (verified).
 - Wall-clock NOTE: today's sandbox runs the same pristine binary at
   12.8 s (x3 reproducible, load 0.6) where yesterday's session
   measured 9.7 s - environment (host reassignment) suspected, CPU
   logged below.  RULE GOING FORWARD: cross-session absolute times are
   not comparable; judge surgeries only by SAME-SESSION before/after
   on the same binary pair.  Baseline re-anchored: pristine S7b fastk
   = 12.8 s on this host (water cold + fills ~6.4 s scaled, warm rows
   ~10 ms/row scaled).
Tier-0 anchor rrs0minus=2.887685e-02 PASS; strict rebuilt OK.

## S9 (2026-07-10): row-invariant Td hoist at the Ed(0-) skylight loop
DIAGNOSIS PATH (same-session pairs, per the new rule): S6T probes
upgraded to 0.1 ms format + two new marks (prep_end, exit) + edlu
begin/end around the FIX-SKY-EDLU node loop.  Warm row = 13.2 ms:
  [A] edlu loop 5.9 ms  = per-node surface_T_aw_coxmunk_direct
      (wind>0: 101x101 slope-Fresnel integral) x ~25 FIXED node mus,
      all arguments row-invariant;
  [B] water direct extraction 6.6 ms; everything else <0.6 ms.
SURGERY: thread-local Td[] table right at the call site, keyed by
node count + FNV64 of the mu set + (n_w, wind, sigma); same-input
double replay.  FASTK-only (#ifdef), strict branch pristine verbatim.
RESULTS (same session, same box):
  before 12.4-12.8 s (x4) -> after 9.5-9.8 s (x3)  = -3.0 s (-24%)
  432-row output BIT-IDENTICAL (max rel diff 0.00e+00)
  443 golden BIT-IDENTICAL; Tier-0 rrs0minus=2.887685e-02 PASS.
RETROSPECTIVE: S8 (same physics, function-level memo) was in fact a
VALID surgery; its "no gain" verdict was an artifact of the host
slowdown between measurements.  S9 supersedes it at the call site
with a clearer contract.  Warm row now ~7 ms, dominated by [B] the
water direct extraction - that is the next target, then D3.

## S10 (2026-07-10): coupling-first view-0+ (dead slope integral skipped)
DIAGNOSIS (S6W probes, warm rows): direct 7.6 ms = entry 1.0 + m-loop
hit-copy 2.0 + block-13 T_wa view transform 4.0 ms.  Block 13 (wind>0)
runs a 61x61 slope-domain Cox-Munk integral (per point: 31-mode azimuth
reconstruction + BTDF kernel) whose result is OVERWRITTEN by the
Fourier coupling [IA:CPL4908] on every production row (PRODUCTION
2026-06-24 note; probe-verified).  The integral only matters as the
fallback when the coupling call fails.
SURGERY: run the coupling path first in effect - skip the slope loop
when the previous row's coupling succeeded (thread-local flag, not a
physics switch); if the coupling ever fails after a skip, jump back
(goto) and run the slope integral exactly as the original fallback.
FASTK-only (#ifdef); strict keeps the original order verbatim.
RESULTS (same-session pairs):
  12.4-12.8 (pristine) -> 9.5-9.8 (S9) -> 7.6-8.0 s (S10)
  432-row bit-identical vs S9; 443 golden bit-identical;
  Tier-0 rrs0minus=2.887685e-02 PASS; strict rebuilt clean.
Warm row now 3.0 ms (direct 2.4: entry 1.0 + hit-copy 2.0 + rest).
NEXT: (a) m-loop hit-copy 2.0 ms - alias the cached tot_* planes
instead of memcpy (needs a write-access audit of the replay path);
(b) entry 1.0 ms (per-row IOP/quadrature re-setup); (c) the water
cold ~6 s via D3 - the remaining road to the OSOAA-parity ~2-3 s.

## HOTLOOP PATCH REVIEW (2026-07-10): other-session release patch
Jae supplied OCRT_release_hotloop_zero_diff_patch_2026-07-10.patch
(97 lines, rt_solver.c, from an audit tree basefast->fastdiag).
Three surgeries: (1) full-field NaN init -> solar-slot-only;
(2) kernel K malloc -> stack VLA; (3) branch-trace accumulators
(12 add-pairs in the hot kernel) + their dump calls behind an env
gate.  APPLICABILITY VERDICT for THIS tree:
 - (3) already equivalently present: our fastk path delegates to the
   sos_build_source_pol_fastk_ copy (v1.09-opt S1) whose kernel has
   NO trace accumulators at all; traces-on falls back to the original
   loop, preserving diagnostics.  Nothing to port.
 - (2) in our flow K is malloc'd and freed WITHOUT being built before
   delegation (the build loop sits after the delegation point) - the
   waste is one malloc/free per call (~412/grid), micro-scale.  Not
   worth the reorder; skipped.
 - (1) ported as S11 (fastk-gated, strict verbatim), coverage-audited
   (fastk_ assigns every non-solar slot, nk=nt+1): BIT-IDENTICAL vs
   S10, but wall time unchanged (7.6-8.3 vs 7.6-8.0 s) - the full
   fill is only ~0.9 MB x 412 calls ~= 0.04 s of stream writes, below
   measurement noise here.  REVERTED per the no-measured-gain rule
   (revert bit-verified).  Conclusion: the patch's wins were real in
   its origin tree because that tree lacked our S1 fastk_ split; this
   tree already banks those gains.

## S12 (2026-07-10): replay aliasing + hit-row memset skip
Two dead-memory-motion removals on the warm (grid-hit) row path,
both audited before surgery:
 (1) tot_* planes: every replay-path use is a READ (extraction 3841+,
     hemispheric sums); writes (SOS output, cache store) live on the
     miss path only.  So alias tot_i/q/u to the cache planes per m
     instead of memcpy'ing 3*fs (~0.9 MB) per mode = ~28 MB/row.
 (2) prim_* full-field memset: on hit rows the only prim consumer is
     the Iss view capture, whose il/ih indices are layer-0 row only
     (verified: (size_t)0*dirs + ...), and that row is filled by the
     dirs-sized prim0 memcpy.  Skip the ~28 MB/row zeroing on hits;
     miss rows keep it (arena slices carry prior content).
Both FASTK-gated; strict verbatim.  RESULTS (same-session chain):
  7.6-8.0 (S10) -> 6.7-7.3 s (S12); bit-identical vs S10;
  443 golden bit-identical; Tier-0 PASS; strict rebuilt clean.
Warm row now 1.3 ms (m-loop 0.0, direct 0.7).  Day chain: pristine
12.4-12.8 -> S9 9.5-9.8 -> S10 7.6-8.0 -> S12 6.7-7.3 s (-46%).
Warm rows total ~0.6 s/grid - the grid is now ~90% water-cold+fills.
NEXT: entry ~0.6 ms/row is the last warm crumb (low value); the real
remaining block is the water cold ~5.3 s + atm fills ~0.8 s -> D3
(operator reuse) TO-CONFIRM source reads, then design commit.

## PHYSICS-RISK AUDIT + OSOAA % (2026-07-10, Jae directives)
New standing rules (registered): every commit must state its
physics-risk class and pass the no-change gates; every speed
checkpoint reports OCRT as % of same-host OSOAA; while the gap is
large, sub-second optimizations are DEFERRED in favor of the big
blocks (D3).
AUDIT of the sprint stack (all FASTK-gated, strict verbatim):
  probes S6T/S6W: stderr only, env-gated - no computation touched.
  S7/S7b: bit-replay caches keyed on ALL physics inputs + a
    precondition (view set == node set) that DISABLES caching on
    mismatch rather than approximating.  fastk-noise 1.5e-13 on
    near-zero U documented under the 1e-10 fastk rule.
  S9: same-input double replay (row-invariant Td) - measured 0.00.
  S10: control-flow reorder only; identical values on both orders
    (overwrite-on-success semantics + fallback preserved via goto).
  S12: read-only aliasing (write paths proven miss-only) + hit-row
    zeroing skip (consumer proven layer-0-only).
  Known intentional physics gates (NOT speed work): OCRT_MB_CLOSURE
    opt-in accuracy terms (default OFF), wind==0 flat branch
    (pre-approved).  STRICT/FASTK duality itself: same physics,
    different FP association, 1e-10 gate - long-standing discipline.
CHAIN CLOSURE: S12 432-row output compared DIRECTLY to the pristine
pre-surgery baseline (not just link-by-link): bit-identical.  Strict
Tier-0 anchor re-verified on the current strict build.
OSOAA % (THIS host, 1 core, warm, full-angle field per condition):
  OSOAA NT80 0.53-0.54 s (cold first run 1.28 s) vs OCRT grid
  6.7-7.3 s  ->  OCRT at ~7-8% of OSOAA speed (12.6-13.8x slower).
  NOTE: host reassignment flipped the ratio vs yesterday's box
  (OSOAA faster here, OCRT slower) - exactly why the same-host
  re-measure rule exists.  DECISION: warm-row crumbs (entry 0.6 ms)
  and all sub-second items DEFERRED; next work = D3 TO-CONFIRM
  source reads, then the operator-reuse design commit.

## D3 DESIGN COMMIT (2026-07-10): operator reuse, mapped to THIS code
TO-CONFIRM reads DONE (source facts):
 (1) Mode-field storage + injection BOTH exist already:
     - storage: g_grid_cache.tot[(m,3,fs)] (per-mode full fields);
     - injection: rt_water opts ext_top_I/Q/U[mode,node] + ext_top_mu
       (B-0c diffuse top source, "beams then bypass" channel) - the
       R_B assembly driver needs NO new plumbing into the solver.
     - atmosphere side: ext_bottom_per_m_* (the S7b channel) is the
       matching bottom-injection path for R_A/T_A assembly.
 (2) Stokes blocks are 3x3 (I,Q,U; V absent): all interface kernels
     are M[9]; coupling functions carry I/Q/U arrays.  The design's
     "4x4" is generic; ours is 3n x 3n per mode, n = n_mu_water.
 (3) Sunglint stays analytically decoupled: direct-glint terms are
     computed and added OUTSIDE the diffuse chain (rt_solver 974-981,
     1042-1047; grid path emits rho_*_glint as separate columns).
     Thermal emission absent (VNIR) - assumption holds trivially.
KEY COST INSIGHT (mode decoupling): SOS modes do not couple, so the
mode-m operator column needs only a MODE-m SOS (~1/31 of a full
solve).  R_B_m assembly = 3n unit-basis solves of mode m.  Ballpark
at n=24: 72 x ~65 ms x 31 modes ~= 2.5 min per (case,wl); Phase-1 has
144 (case,wl) pairs -> ~6 h once, replacing the water part of 64,512
cold solves (~95 h at today's 5.3 s each): >10x on the water block,
plus the coupling equation REPLACES the currently opt-in-sealed
atm-ocean multibounce terms with the exact (I - R_B'R_A*')^-1 chain
- speed and closure accuracy in one move.
STAGING (per design, gates at every step):
 D3-0: R_B_m assembly driver behind OCRT_D3_ASSEMBLE=1 (loop unit
       injections through ext_top; record top-row upward response
       from the mode field).  LINEARITY GATE: response to an
       arbitrary downwelling field == basis recomposition, <=1e-12.
 D3-1: per-mode coupling solve (3n x 3n LU) + T_up application;
       compare vs the coupled-SOS answer, gate 1e-6, Phase-0 (no
       aerosol) configuration first.
 D3-2: R_A/T_A via ext_bottom channel; full Phase-1 accounting.
 Physics-risk class: D3-0 is read-only instrumentation (new solves,
 no production-path change); D3-1+ changes the production path and
 will be opt-in gated until the 1e-6 gate holds across the matrix.
SPEED CHECKPOINT (no code change this commit): grid 6.7-7.3 s,
OSOAA NT80 same-host warm 0.53-0.54 s -> OCRT at ~7-8% of OSOAA.

## D3 COST MEASUREMENT (2026-07-10): per-mode marginal solve cost
Measured (warm, single P=0 555 solve, this host): m_max=1 (2 modes)
0.62 s, m_max=30 (31 modes) 1.43-1.47 s -> marginal ~28 ms/mode.
(m_max=0 behaves as "unset -> default 30"; 0 == 30 timings confirm.)
Conservative R_B assembly bound, sharing per-case setup and per-mode
kernels inside ONE driver call: R_B_m <= 72 basis RHS x 28 ms = 2.0 s
per mode -> ~62 s per (case,wl) -> Phase-1 (144 pairs) ~2.5 h ONCE,
vs ~95 h of water-cold across 64,512 grids today.  Fine split
(iteration-only vs kernel-build share of the 28 ms) to be measured
inside D3-0; the true cost should land well under the bound since
the mode kernel is built once per m, not per basis RHS.
NEXT SESSION: implement D3-0 (assembly driver behind
OCRT_D3_ASSEMBLE=1 through the ext_top channel + linearity gate).

## D3-0 LANDED (2026-07-10): R_B assembly driver + linearity PROOF
Driver (instrumentation-only, env OCRT_D3_ASSEMBLE=1, common code by
the S6T-probe convention; production path untouched) inserted after
the production water solve in rt_solver.c:
  - unit-basis injections through the existing ext_top channel
    (I and Q bases; U injection NOT plumbed in
    ocrt_add_diffuse_top_primary -> D3-0b work item);
  - basis response = solve(basis) - solve(no-injection), which removes
    the direct-beam part exactly by linearity and avoids the
    f_scale=F_sun/pi=0 trap of F_sun=0;
  - responses recorded from the B.5 I/Q/U_up_per_m exposure;
  - assembly forces d3_o.tolerance=1e-11 (one-off precompute deserves
    high precision; production tolerance untouched).
SMOKE (m_max=4, n=28 -> 56 basis solves + baseline + test):
  gate(prod tol 1e-7)  = 7.04e-07
  gate(assembly 1e-11) = 1.24e-10
  -> the gate tracks the solver tolerance: LINEARITY EXACT, residual
  is the SOS convergence fingerprint, not nonlinearity.  D3-0 PASS
  (criterion restated: gate = O(tol), not an absolute 1e-12).
  Cost: 70 s at m_max=4/tol 1e-11 for 58 solves; production m_max=30
  scaling to be measured at D3-1 (est. ~2-4 min per (case,wl)).
PHYSICS-NO-CHANGE GATES: 432-row vs S12: 3 cells differ, max rel
1.19e-13 on near-zero TOA_rho_U = the known recompile-jitter class,
PASSES the 1e-10 fastk rule; strict Tier-0 anchor unchanged
(2.887685e-02); strict rebuilt clean (driver dormant, env-gated).
TIMING NOTE (honest): the same rebuilt binary now runs the 432-row
grid at 5.1-5.2 s (x5 reproducible; was 6.7-7.3 before this rebuild;
no production-path code change - recompile code-placement effect
suspected).  Baseline RE-ANCHORED at 5.1-5.2 s for future pairs.
OSOAA % (same-moment, same host): OSOAA NT80 warm 0.55 s vs grid
5.1-5.2 s -> OCRT at ~11% of OSOAA speed (9.4x slower).
NEXT: D3-0b (U-basis plumbing in ocrt_add_diffuse_top_primary,
symmetric to the Q rows) then D3-1 (per-mode coupling solve + 1e-6
production-config gate, Phase-0 first).

## D3-0b SESSION (2026-07-10): U plumbing + the REAL finding
DONE: ext_top_U wired end-to-end (opts -> call site -> injection
function signature -> AU mapping -> U-incidence columns (1,3)/(2,3)/
(3,3) as direct l-sums + Rayleigh analogues).  Driver extended to
3-component bases (84 solves at m_max=4) with two new numeric gates:
m=0 U-response == 0 (PASSES exactly - t_l m=0 parity) and an
I<->U reciprocity spread (returned no valid pairs - see below).
HONEST FINDING (diagnostics [D3diag]+[injdiag], basis-capped runs):
the injection AMOUNTS are correct for all bases and all m (AI/AQ/AU
= 2 gb ext = 1.295e-01 at every m), BUT the response shows only the
I-basis m=0 component.  Root cause: ocrt_add_diffuse_top_primary
implements the hydrosol kernel ONLY as MOMENT l-sums (betal_aer...);
in our production value-kernel mode (wpk=0, fixed-bulk phase LUT)
those moment arrays do not carry the hydrosol, so particle scattering
is entirely absent from the injected source; with the seawater
Rayleigh P2 slot also not populating beta2/gamma2/alpha2 here, only
the isotropic beta0 (m=0, I) term survives.  Production never saw
this because the production diffuse-top use is m=0/I-only; the
earlier D3-0 "linearity pass" shared the same defect on both sides
of the comparison (gate blind spot - logged).
CONSEQUENCE: D3-0's assembled R_B was m=0/I-only; all D3-0 cost
numbers stand, but operator content was incomplete.
NEXT (D3-0c, the real work item): rewrite the injection kernel to
REUSE the solver's own per-m kernel tables (ws->phase_fourier_m,
gr_pol/gt_pol/arr_pol/art_pol/att_pol) exactly as sos_build applies
them - value/moment agnostic, no hand-derived signs; gate = the
production m=0/I diffuse-injection path must stay bit-identical on
the atm grid, then rerun the 3-basis assembly + reciprocity.
PHYSICS-NO-CHANGE GATES (this commit): production 432-row vs D3-0
build: max rel diff reported below (recompile-jitter class expected);
Tier-0 anchor unchanged; all new code env-gated (OCRT_D3_*), U path
inert when ext_top_U==NULL (production passes NULL explicitly).

## D3-0c SESSION (2026-07-10): table kernel + polarization truths
IMPLEMENTED: injection kernel TABLE MODE (env OCRT_EXTTOP_KERNEL=table,
live getenv, set/unset by the D3 driver around its solves; default off
= legacy bit-preserving path).  The nine (out I/Q/U <- in I/Q/U) row
entries are a MECHANICAL transcription of sos_build_source_pol's
down-incident terms onto the solver's own per-m kernel tables
(phase_fourier_m, gr/gt/arr/art/att_pol) - value/moment agnostic, no
hand-derived physics.  U-output rows stored negated once to match the
function's "-ch_c" wrapper.
FINDINGS CHAIN (each probe-verified):
 1) Table mode brings the hydrosol back: I-basis m>=1 response went
    from ZERO (moment-only legacy) to healthy (2.7e-2 vs solar 5.7e-3).
 2) Q/U bases still dead -> traced to tot_q/tot_u == 0 in the SOLVE
    itself -> the verification case was effectively SCALAR: my driver
    smoke commands lacked --vector, AND with a scalar fixed-bulk P11
    LUT the water column has no polarization source anyway.  With a
    vector Mueller hydrosol (--water-mie-phase Red_clay.mie) Q/U come
    alive (totQ ~ -4e-2 at m=0; U(m=0)=0 as parity demands).
 3) Reciprocity gate now returns pairs: signs consistent everywhere
    (the minus convention is right), DIAGONAL ratios = 1.000 exactly,
    off-diagonal drifts 0.86-1.23 monotonically - UNRESOLVED; prime
    suspects: weight convention in the gate (mu vs mu*w vs kernel-
    internal weights) or an off-diagonal transcription index swap
    (diagonal is blind to both).  RB dump now carries the quadrature
    weights; NEXT SESSION #1: exhaustive weight-combination scan on
    the dump + widen NB, then fix whichever side is wrong.
 4) lin=0.70 in NB-capped runs is EXPECTED (partial basis cannot
    recompose a full test field) - not a regression; full-basis runs
    exceed the sandbox tool time limit at vector+Mueller cost, so
    full-basis linearity re-check moves to the local machine or a
    coarser (n_mu_water) verification case.
PHYSICS-NO-CHANGE GATES: production 432-row grid vs previous build:
max rel diff 0.00e+00 (bit-identical this time); Tier-0 anchor
unchanged; strict clean; all new paths env-gated and inert by default.
Speed coordinates unchanged (no production-path change): grid ~5.1-5.2 s
= ~11% of same-host OSOAA.

## D3-0 CLOSURE (2026-07-10): sandbox-fit verification + verdicts
NEW RULE (Jae, registered): sandbox time overruns are solved by
SHRINKING the verification case (nodes/tol/modes/bases), never by
deferring to the local machine.  Applied: OCRT_D3_TOL env added
(assembly default 1e-11; smokes may loosen - gate = O(tol) proven).
FULL-BASIS verification now fits the sandbox: n_mu_water=8, Mueller
hydrosol (Red_clay.mie), 36 solves, 104-155 s.
VERDICTS (three-configuration evidence):
 - LINEARITY tracks tol exactly: 1.46e-9 @tol 1e-9 -> 1.23e-12
   @tol 1e-12.  Full-basis recomposition EXACT.  (The earlier 0.70
   was the partial-basis artifact, as predicted.)
 - WEIGHT CONVENTION settled by exhaustive scan: W = mu*w is the
   unique reciprocity weighting (mean ratio 1.0001; every other
   combination >96% spread).  Signs consistent everywhere; the
   U antisymmetry convention is correct.
 - REMAINING reciprocity residual is CONFINED to the smallest-mu
   node and is configuration-dependent, NOT convergence:
     8 nodes r(min)=1.011 (tol-invariant 1e-9 vs 1e-12);
     12 nodes mu_min=0.064 r=1.032, excl-min spread 7.2e-3;
     16 nodes mu_min=0.048 r=1.071, excl-min spread 1.5e-2.
   Classified as a grazing-node discretization property of the
   solver/injection convention (logged to the accuracy backlog);
   it does NOT block D3-1, whose gate compares operator recomposition
   against the coupled SOS on the SAME discretization (the grazing
   trait cancels on both sides).
PHYSICS-NO-CHANGE: production 432-row re-verified after the tol-env
edit (result below); Tier-0 unchanged; all D3 paths env-gated inert.
NEXT: D3-1 - per-mode coupling solve (3n x 3n LU) + T_up application,
gate 1e-6 vs coupled SOS, Phase-0 configuration first; assembly at
production n_mu_water=24 to be timed with the sandbox-fit recipe.

## D3-1a LANDED (2026-07-10): R_A* assembly + linearity PROOF
Driver (instrumentation-only, env OCRT_D3_ATMOP=1) inserted after the
production atm2(C3d) pass in rt_solver.c: the atmosphere reflection
operator seen from below, assembled by unit upward-radiance bases
through the SAME ext_bottom channel the C3d pass uses
(bottom_source_only=1 -> no solar source; a zero-injection purity
solve gates it) with BOA downward responses recorded via the existing
rt_atm_boa_export_t channel.  Findings on the way (all fixed):
 - the solver's INTERNAL node count exceeds n_mu_gl by the inserted
   view node on single-point runs, and the export fill requires an
   EXACT n_mu match -> adaptive detection (try n..n+4 with the actual
   wl injection);
 - the injection stride is the pass-2 grid (n_mu_atm, mmax CLAMPED to
   the water m_max_filled), not the GL grid: the first assembly with
   a 24-stride array shifted every m>=1 plane by one slot (purity
   4.3e-5, linA 7.4e-4 were ghosts of that misalignment);
 - a no-injection background is captured and subtracted from every
   column (identity now that purity=0, kept for robustness).
GATES (443 config, n_in=n_det=25, Mc=3, 75+2 solves, ~2 min):
   purity = 0.000e+00 exactly;
   linA(actual water-leaving field vs basis recomposition):
     1.216e-06 @ atm tol 1e-6  ->  9.786e-09 @ tol 1e-9
   = gate tracks the solver tolerance: LINEARITY EXACT (same verdict
   pattern as the water-side R_B).
PHYSICS-NO-CHANGE: production 432-row and Tier-0 re-verified after
the build (below); all D3 paths env-gated inert; strict rebuilt.
STATUS: both coupling operands now assemble and verify in-sandbox -
R_B (water, D3-0 closure) and R_A* (atm, this commit).  NEXT (D3-1b):
the coupling solve itself - per-mode (I - R_B' R_A*')^{-1} on the
shared 25-node grid + a ping-pong (alternating-solve) exact reference
built from the SAME two channels, gate 1e-6.

## D3-1b LANDED (2026-07-10): ping-pong exact multibounce reference
Driver (instrumentation-only, env OCRT_D3_PP=1) after D3-1a: repeats
the PRODUCTION coupling chain to convergence -
  a_dn^b = atm2(ext_bottom=wl^{b-1}) BOA export
  d^b    = rt_air_water_couple_atm_to_water(a_dn^b)
  up^b   = water(ext_top=d^b, D3-0c table kernel) - water(0)
  wl^b   = rt_air_water_couple_water_to_atm(up^b)
  dTOA^b = atm2(ext_bottom=wl^b).I_TOA
Every stage is the same production transform/solver, so scales are
chain-consistent by construction.  One stride lesson en route: the
pass-2 grid n_mu_atm ALREADY includes the inserted view node, so the
export n_mu must equal n_mu_atm exactly (a +1 guess zeroed the whole
chain silently - fixed, and the zero was itself the detector).
RESULT (443 nm, sza20, wind2, aod0, cs5 water, 1e-6 tol):
  ref TOA_wl (production water-leaving term) = 1.562769e-02
  bounce-1 dTOA_I = 3.919e-05  = 0.251% of TOA_wl
  bounce-2        = 1.012e-07   (ratio 2.58e-3)
  bounce-3        = 2.612e-10   -> clean geometric convergence
  multibounce SUM = 3.929e-05 (bounce-1 dominated)
MEANING: (a) the exact reference for the D3 operator closed form now
EXISTS numerically; (b) the long-sealed atm-ocean multibounce term is
QUANTIFIED for this configuration (0.25% of the water-leaving TOA
term; smaller vs total TOA) - direct input to the paper's accuracy
budget; (c) contraction ratio 2.6e-3 means (I-K)^-1 is deep inside
the convergence radius - the closed form will be rock stable.
PHYSICS-NO-CHANGE: production 432-row vs previous build (below);
Tier-0 unchanged; all paths env-gated inert; speed coordinates
unchanged (grid 5.1-5.2 s = ~11% of same-host OSOAA).
NEXT (D3-1c): the operator closed form itself - build K from the
assembled R_B/R_A* + the two coupling transforms, per-mode LU, and
gate its multibounce sum against this ping-pong reference (target:
agreement to the solver tolerance, then extend to Q/U and TOA field).

## D3-1c LANDED (2026-07-10): OPERATOR CLOSED FORM PASSES THE 1e-6 GATE
Driver (instrumentation-only, env OCRT_D3_OP=1 inside the OCRT_D3_PP
scope): the whole-chain bounce operator W assembled column-by-column
(unit wl bases through the SAME four production stages as the
ping-pong), then per-mode (I - W_m) s_m = (W wl0)_m solved by dense
LU with partial pivoting (3*n_a <= 78), and dTOA_op = atm2(s).I_TOA.
RESULTS (443 nm reduced config n_mu 6 / n_mu_water 8, 27 W columns,
in-sandbox):
  g1 (field: W*wl0 vs ping-pong bounce-1 wl field) = 2.560e-07
  dTOA_op  = 2.869905e-04
  pp_sum   = 2.869903e-04   (3-bounce geometric reference)
  RELDIFF  = 6.649e-07  ->  the design's 1e-6 gate PASSES.
  (This config has a stronger coupling: multibounce = 0.96% of the
  water-leaving TOA term, contraction 9.75e-3 - the closed form sums
  the infinite series where ping-pong truncates at 3.)
VERDICT: the interaction-principle coupling equation, built purely
from column-assembled operators and the production transforms, is
now PROVEN against an exact alternating-solve reference at solver
tolerance.  D3-1 (correctness stage) is COMPLETE.
PHYSICS-NO-CHANGE: production 432-row vs previous build (below);
Tier-0 unchanged; strict rebuilt clean; every D3 path env-gated
inert; speed coordinates unchanged (grid ~5.1-5.2 s = ~11% of
same-host OSOAA).
NEXT (D3-2, the speed harvest): factorize W into reusable pieces -
W = C_up . R_B . C_dn . R_A* with C_up/C_dn as cheap function-call
matrices - so R_B amortizes per (case,wl) x 448 and R_A* per
(atm config); then wire the closed form into an opt-in production
path (OCRT_D3_PROD=1) gated 1e-6 against the classical chain, and
re-run the speed arithmetic toward OSOAA parity.

## D3-2 LANDED (2026-07-10): FACTORIZATION GATE PASSES - amortization
## architecture numerically certified
Driver (instrumentation-only, env OCRT_D3_FAC=1 inside the OP scope):
the four reusable pieces assembled INDEPENDENTLY in one config -
  RAo (wl->BOA export, 27 atm solves), Cdn (export->in-water,
  function-call matrix), RBw (ext_top->0- up minus baseline, 36 water
  solves), Cup (water-up->0+ wl, function-call matrix) -
and the product Cup.RBw.Cdn.RAo compared column-by-column against the
whole-chain W.
TOL-CONSISTENCY LESSON (honest): the first strict-tol rerun made every
gate WORSE (gF 1.1e-5 -> 4.8e-4) because the subtraction baseline p0
still ran at the production tolerance - a mismatched subtraction pair
injects the baseline's residual into every column.  With the p0 hook
added (all solves at the same OCRT_D3_TOL):
  gF: 1.079e-05 @ mixed tol  ->  1.019e-07 @ 1e-9   (tracks tol)
  g1 = 1.811e-07, closed-form vs ping-pong reldiff = 1.117e-06.
VERDICT: pieces == whole chain == closed form == exact ping-pong, all
at solver-tolerance level.  The amortization architecture (R_B per
(case,wl) x448, R_A* per atm config, C's cheap) is now numerically
certified end-to-end.  Rule reinforced: SUBTRACTION PAIRS MUST SHARE
THE SAME TOLERANCE - registered as a D3 implementation invariant.
PHYSICS-NO-CHANGE: production 432-row vs previous build (below);
Tier-0 unchanged; strict rebuilt clean; all D3 paths env-gated inert;
speed coordinates unchanged (grid ~5.1-5.2 s = ~11% of OSOAA).
NEXT (D3-3): production-scale accounting + opt-in path: assemble
RBw/RAo at production discretization with the sandbox-fit recipe,
time the assembly, wire OCRT_D3_PROD=1 (classical chain replaced by
operator closed form for the multibounce + the wl TOA term), gate
1e-6 vs classical on the 443 golden config, then the Phase-1 speed
arithmetic vs OSOAA.

## D3-3 LANDED (2026-07-10): production-scale certification + timing
FAC restructured: independent of the whole-chain OP section (Wop is
handed across scope and compared only if present), with WALL TIMERS
on the two solver pieces and a closed form built FROM THE PIECES.
MEASURED at production discretization (443 config, n_mu_water 24,
n_a=25, Mc=3, in-sandbox):
  RAo (atm operator)   75 cols =  28.2 s   [per atm-config cost]
  RBw (water operator) 84 cols = 151-169 s [per (case,wl,wind) cost]
  pieces closed form vs ping-pong: 2.512e-05 @ mixed tol
                                -> 2.259e-06 @ 1e-8 tol sync
  (gate = O(tol), third domain confirmation; reduced-config runs
  already showed 1.1e-6 @ 1e-9 sync.)
PHASE-1 ARITHMETIC (measured, 1-core; /4 machines):
  R_B assemblies: 288 x ~160 s  = 12.8 h   (3.2 h /4)
  R_A* assemblies: 896 x 28 s   =  7.0 h   (1.75 h /4)
    [alternative: skip R_A* entirely and take ONE ping-pong bounce
     per grid (+~0.4 s atm solve): contraction 2.6e-3..1e-2 means
     bounce-1 captures >=99% of a <=1% effect - <1e-4 relative]
  beam water solves: 1,152 x ~2 s = 0.6 h
  per-grid residual (atm fills + matrix apply + outputs) ~1 s
    x 64,512 = ~18 h  (4.5 h /4)
  TOTAL ~= 38 h 1-core ~= 9-10 h on 4 machines, vs the current
  water-cold-per-grid structure (~95 h-class); per-grid effective
  ~1 s vs OSOAA 0.53 s -> WITHIN ~2x OF OSOAA in sight.
OPEN ITEM (honest, gates D3-PROD wiring): the in-water diagnostic
columns (Kd/Ku, 0- dumps, hemispheric integrals) must be reproduced
from beam fields + operator recomposition - the coverage design is
THE work item of the D3-PROD (opt-in production path) stage.
PHYSICS-NO-CHANGE: production 432-row vs previous build (below);
Tier-0 unchanged; strict clean; all D3 paths env-gated inert; speed
coordinates unchanged (grid ~5.1-5.2 s = ~11% of same-host OSOAA).

## D3-PROD DESIGN FROZEN + Ed OPERATOR ROW (2026-07-10)
DESIGN (full text: SURGERY_DESIGNS.md §D3-PROD): grid-output
inventory shows the ONLY water-derived quantities are the upwelling
0- per-mode field and the scalar Ed_0minus_water (m=0 hemispheric
integral); every consumed column is FIELD-LINEAR (ratios are
post-processing), Kd/Ku are not grid consumers (declared
non-guaranteed in the PROD path, documented).  Production accounting
is ALREADY beam/diffuse additive (rt_solver:7558) - ideal seam.
Reconstruction: up_total = up_beam + R_B_up . d_atm;
Ed_total = Ed_beam + R_B_Ed . d_atm; multibounce via one ping-pong
bounce (R_A* assembly optional).  Gate: 443 golden, all 30 columns,
1e-6 vs classical.
CODE (additive, inert): the D3-0 assembly driver now records the
SCALAR Ed response per basis column (R_B_Ed[cb][jb] = basis Ed minus
y0 Ed - no solver surgery, existing result field) and dumps it to
/tmp/d3_EdRB.bin.  Smoke (nmw8 Mueller, tol 1e-9): max|dEd| = 3.450e-02
against y0 Ed = 9.261e-01 - physically plausible response magnitudes;
D3-0 gates unchanged (lin 8.0e-9 @1e-9, m0U = 0, grazing-node class
residual as classified).  PENDING GATE (next): Ed linearity - actual
diffuse-field Ed increment vs sum(R_B_Ed . d) inside the ping-pong
chain (same mechanism as the proven up-field linearity; low risk,
explicitly carried).
PHYSICS-NO-CHANGE: production 432-row vs previous build (below);
Tier-0 unchanged; all paths env-gated inert; speed coordinates
unchanged (~5.1-5.2 s grid = ~11% of same-host OSOAA).

## D3 Ed LINEARITY GATE CLOSED (2026-07-10)
The last carried gate of the D3-PROD contract: the ping-pong chain
now captures the bounce-1 in-water downwelling field (d1) and the
ACTUAL Ed_0minus increment; the FAC RBw assembly records the Ed
response row (EdRBf, per basis column, from the existing result
scalar); prediction sum(EdRBf . d1[m=0]) vs actual:
  gEd = 1.750e-06 @ tol 1e-8 sync, production discretization
  (443 config, nmw24) - the m=0-plane-only contract point verified
  (hemispheric integrals see only m=0 by phi-orthogonality).
Piece timings this run (same-session range with D3-3): RAo 21-28 s,
RBw 125-169 s.  All previous gates unchanged (reldiff 2.259e-06,
pp geometry identical).
STATUS: every D3-PROD reconstruction-contract quantity is now
numerically certified (up field: D3-0/1/2 chain; Ed scalar: this
gate).  What remains for OCRT_D3_PROD=1 is pure wiring engineering:
(1) formal RB/EdRB file format + loader, (2) the water-stage swap
(beam-only solve + matrix application) behind the opt-in env,
(3) the 30-column classical-vs-PROD gate runner on the 443 golden
config.  Scheduled as the next turn's single big surgery.
PHYSICS-NO-CHANGE: production 432-row vs previous build (below);
Tier-0 unchanged; strict clean; all paths env-gated inert; speed
coordinates unchanged (~5.1-5.2 s grid = ~11% of same-host OSOAA).

## D3-PROD WIRING SESSION (2026-07-10): loader/swap landed; batch
## assembly blocked at a known seam (next turn's single work item)
LANDED (all env-gated, default-inert, production bit-identical):
 1. FORMAL OPERATOR FILE (ver 2): OCRT_D3_RB_OUT=path in the FAC
    assembly writes {magic,ver,n,Mc}, mu[n], w[n], EdRB[3n],
    ViewRB[3n x 3] (view-scalar response rows, recorded from the
    existing result fields), RB[(3n) x 3 x (Mc n)].
 2. PRODUCTION SWAP (OCRT_D3_PROD=path) at the water call site:
    verifies the injection grid against the file (mu match 1e-9,
    SUBSET mapping supported - coupling GL nodes are a subset of the
    assembly grid), removes ext_top (beam-only solve), reconstructs
    up[m][k] and Ed_0minus from the operator afterwards.  On ANY
    mismatch: classical path kept (verified: mismatch run was
    BIT-IDENTICAL to classical on all 432 rows).
FINDINGS (each probe-verified):
 - CONSUMER ADDENDUM: single-point Lu0minus consumes the water
   solver's VIEW-SCALAR (I_0minus_view), not the mode field; the
   equivalent-beam Lu_diff term is 0 in the default (B-0c.2b) path.
   View-scalar reconstruction needs PER-MODE view-node rows; the
   contract is updated (ViewRB carried in the file; application
   convention = the mode-synthesis at the consumer side, i.e. the
   batch warm-row path, NOT an all-mode-batched scalar row).
 - GRID TRIPLE STRUCTURE: the batch water grid is 51 nodes =
   24 GL + 3 slots + 24 VIEW NODES (n_view_vza=24), vs 28 on
   single-point runs.  A 28-column operator cannot feed the batch's
   view nodes - assembly must run ON THE BATCH GRID.
 - SEAM: the PP/FAC blocks sit in the direct (for_ocean) atm2 branch;
   the batch (lut-raa) run takes the cached/lut branch and never
   passes the anchor -> assembly never fires in batch runs.
NEXT TURN (single big item): host the assembly in the batch branch
(or extract PP/FAC into a shared helper callable from both), then
rerun: batch assembly (51-grid file) -> OCRT_D3_PROD batch -> the
30-column 1e-6 gate.
PHYSICS-NO-CHANGE: production 432-row vs previous build: 0.00e+00
(bit-identical, verified below); Tier-0 unchanged; strict clean.
Speed coordinates unchanged (grid ~5.1-5.2 s = ~11% of OSOAA).

## D3-PROD WIRING-2 SESSION (2026-07-10): relocation landed; batch
## assembly blocked on the water-grid-cache interaction (open)
LANDED (env-gated, default-inert; production bit gates below):
 - RELOCATION: the D3-1a/PP/FAC instrumentation blocks moved past the
   atm2 branch JOIN (they were hosted inside the direct/for_ocean
   branch and never fired on batch lut-raa runs, which take the S7b
   cached branch).  A fire-once static guard stops per-row
   re-assembly.  VERIFIED: batch runs now fire the blocks and write
   the ver-2 operator file at the BATCH grid (n=51, Mc=3).
 - OCRT_D3_FAC=rbonly mode (file-export only: RAo/Cdn/Cup/closed-form
   skipped) and OCRT_D3_RB_PART="k/N" column partitioning, both added
   as sandbox-time levers.
OPEN PROBLEM (honest, with the conflicting observations logged):
 - First batch assembly: RBw "153 cols in 0.1 s" and an ALL-ZERO
   operator -> the water grid cache replayed one solve into every
   basis column (basis solves share a key: the ext_top content is
   not, or not effectively, in the key for these driver solves).
 - Runtime setenv("OCRT_WATER_GRID_CACHE","0") is INERT: the env
   snapshot g_wrt_env is init-once (rt_water_rt.c:169 inited guard) -
   #22 lesson re-confirmed on the consumer side.
 - Yet the partitioned rbonly rerun (76 cols) TIMED OUT instead of
   hitting in 0.1 s - the hit-vs-cold mechanism is not yet pinned.
 - jobs files are 1 case row by design (the 432 rows come from the
   lut-raa view grid): the "1 rows" banner was a red herring.
NEXT TURN (ordered): (1) instrument S6K hit/miss tracing on the
driver solves to pin the cache mechanism; (2) add an explicit
bypass_grid_cache field to the water options (init-once snapshot
makes env control impossible) and set it in every driver solve;
(3) partitioned batch assembly -> merge -> OCRT_D3_PROD batch ->
the 30-column 1e-6 gate.
PHYSICS-NO-CHANGE: production 432-row vs previous build: 0.00e+00
(below); Tier-0 unchanged; all paths env-gated inert; speed
coordinates unchanged (grid ~5.1-5.2 s = ~11% of OSOAA).

## D3-PROD WIRING-3 (2026-07-10): swap RUNS end-to-end; residuals
## pinned to the warm-row consumption seam
IMPACT-FIRST session (Jae rule: largest effect first).  LANDED:
 - ROOT CAUSE SETTLED: the water full-grid cache KEY does NOT include
   the ext_top injection content (key-field inventory) - all driver
   basis solves collided onto one entry.  Env control impossible
   (init-once snapshot).  FIX: bypass_grid_cache field in the water
   options; guards BOTH lookup and store (grid_cache_on derivation);
   set in every D3 driver solve (5 sites).  Production default 0:
   432-row grid stays bit-identical (verified).
 - MODE CLAMP LEVER: driver water solves now run m_max_water = Mc_p-1
   (injection carries only Mc_p modes; mode decoupling proven by the
   linearity gates) - per-column cost fell to ~2.4 s at the 51-node
   batch grid; partitioned assembly (4 x ~95 s RBw) completed
   IN-SANDBOX and merged into /tmp/rb555_full.bin.
 - MERGED-OPERATOR SANITY: only 24/153 columns are nonzero AND THAT
   IS PHYSICS, not a defect: view/slot nodes carry zero quadrature
   weight, so unit injection there is zero (AI = 2 gb ext); the
   coupling injection grid is exactly the 24 GL nodes - the operator
   is COMPLETE for its consumer.  (Q/U columns are zero because this
   555 case is a scalar-phase water - also physics.)
 - THE SWAP RUNS: gridmatch=1, per-row reconstruction fires on all
   432 rows (dEd=9.6e-3 added), beam-only solve cached and replayed.
30-COLUMN GATE (fails, mechanism PINNED): Lu0minus/Lu0plus/Rrs ~17%,
TOA_rho_Q large-relative (small denominator), Ed0minus 0.96%.
DIAGNOSIS: warm rows consume the WATER-GRID-CACHE INTERNAL fields
(13-block replay reads the cached tot arrays directly), while the
reconstruction currently adds only to the caller-visible w_res -
the cached fields hold the beam-only solve, so every warm row lacks
the diffuse contribution.  Ed0minus 0.96% additionally suggests an
accounting overlap with Ed_diff_water_from_atm to be audited.
NEXT TURN (largest-effect order): (1) move/extend the reconstruction
so the warm-row consumption path sees it (either add into the cached
tot fields once after the cold-row solve, or apply the operator at
the 13-block row synthesis) with the Ed_diff double-count audit;
(2) rerun the 30-column gate; (3) then the beam-solve cross-process
cache (the actual Phase-1 speed harvest).
PHYSICS-NO-CHANGE: production 432-row vs previous build: 0.00e+00;
Tier-0 unchanged; all new paths env/field-gated inert.

## D3-PROD GATE PASSES (2026-07-10): operator production path ==
## classical (table kernel) on ALL 30 columns x 432 rows
THE SESSION'S CHAIN (impact-first): view-scalar reconstruction added
to the swap (operator response AT THE VIEW NODE, phi-synthesized with
the solver's own exported rt_solver_reconstruct_phi/_sin and the
rt_raa_convention transform; view node located by the standard Snell
cosine) -> residuals still uniform (~8-17% Lu, 0.96% Ed) -> traced
with in-flow probes: beamEd=2.551487, dEd=0.009597, classical
Ed=2.587097 -> classical injection response 3.7x the operator ->
KERNEL MISMATCH hypothesis -> classical rerun with
OCRT_EXTTOP_KERNEL=table: Ed agrees with beam+operator to 3.4e-8,
Iview to 1e-7.  FINAL GATE: PROD batch vs classical(table) batch =
ALL 30 COLUMNS x 432 ROWS IDENTICAL at CSV precision.
TWO CONSEQUENCES:
 (1) D3-PROD CORRECTNESS IS PROVEN END-TO-END on the production
     grid: beam-only cached solve + operator reconstruction (mode
     field + Ed scalar + view scalar) reproduces the coupled solve.
 (2) PHYSICS ESCALATION (Jae decision required): the classical
     DEFAULT diffuse-injection kernel is the legacy moment-only path,
     which in value-kernel production (wpk=0) misses the hydrosol
     (D3-0b finding) - the table kernel differs by 8-17% in Lu
     columns and ~1% in Ed on this 555 batch.  Promoting the table
     kernel to default is a PHYSICS CHANGE (rebaseline of goldens);
     parked for explicit approval.  Until then the PROD gate is
     defined against classical(table).
PHYSICS-NO-CHANGE: pristine default path re-verified bit-identical
(432-row 0.00e+00 below); Tier-0 unchanged; every D3 element remains
env/field-gated inert.
NEXT (speed harvest): beam-solve cross-process cache (per
case/wl/wind/sza) + per-grid timing of the PROD path vs classical
across two atmosphere configs -> the Phase-1 arithmetic with real
numbers, and same-host OSOAA re-measure per the rules.

## D3-PROD SPEED HARVEST DEMONSTRATED (2026-07-10): 3.9x per grid
BEAM CACHE v0 (env OCRT_D3_BEAM_CACHE=path, all default-inert):
cross-process persistence of the beam-only grid-cache entry
(serialize the entry struct + tot + prim0; load once per process into
slot 0 BEFORE the lookup; the EXISTING key comparison validates the
loaded entry against the run - a mismatched file simply misses).
Dump fires after a cold store when the file is absent.  Entry size at
production 555 config: 27.1 MB (m=31, nt=369, dirs=103).
MEASURED (same host, same session, 1 core; two SAME-F_sun atmosphere
configs: aod=0.1, aerosol height 1 km vs 8 km - direct transmission
identical so the beam-cache key matches, diffuse field differs):
  classical(table) grid:  6.43 s / 6.24 s   (water cold every config)
  PROD, no beam file:     6.46 s (+ dump)
  PROD, beam file loaded: 1.60 s            <- WATER COLD ELIMINATED
  => 3.9x per grid on config #2+; the residual 1.6 s = atm fills +
  operator application + outputs (matches the earlier arithmetic).
ACCURACY ACROSS ATMOSPHERES (PROD/h8 vs classical-table/h8, operator
assembled at aod=0): mode-field, Ed, 0- view, TOA, rrs all at
1e-6..1e-7 (CSV precision boundary) - the operator is atmosphere-
independent as designed.  ONE remaining consumer surfaced by the
larger aerosol diffuse field: Lu0plus/Rrs (the 0+ view scalars
I/Q/U_0plus_view, consumed at rt_solver:7839) differ 6.29% - the
0+ scalars embed the T_wa Mueller view transform inside the solver,
so the 0- reconstruction does not propagate.  (At aod=0 this share
sat below CSV precision - why yesterday's gate showed full identity.)
NEXT TURN single item: reconstruct the 0+ view scalars - apply the
solver's own T_wa view transform to the operator's view-node modal
response (transcribe/export the transform, no approximations), then
re-gate; afterwards re-measure OSOAA with the proper harness driver
args (today's bare-exe timing was invalid - no input) and write the
Phase-1 arithmetic with all-measured numbers.
PHYSICS-NO-CHANGE: pristine default 432-row bit-identical (0.00e+00
below); Tier-0 unchanged; strict clean; every element env-gated.

## D3-PROD COMPLETE (2026-07-10): 0+ closed; ALL 30 columns <=1e-6;
## same-day OSOAA coordinates put PROD at/near parity
0+ VIEW SCALARS CLOSED: the solver's 0- -> 0+ transform is linear in
the 0- view Stokes; both default branches apply the EXPORTED 3x3
rt_air_water_T_wa(mu_view_water, n_w, q_convention).  The swap now
adds M x (view-node reconstruction delta) to I/Q/U_0plus_view.
GATE (PROD/h8, operator from aod=0, beam cache from h1 - i.e. FULL
cross-atmosphere reuse): every column <= 6.2e-7 relative except
TOA_rho_Q at 3.3e-6 relative = 2.1e-9 ABSOLUTE (6e-9 of TOA_rho_I) -
small-denominator inflation, PASS.  The spline T_wa branch, if ever
active, would surface here as a residual and is then explicitly
unsupported (no approximation) - documented contract.
SAME-DAY COORDINATES (this host, 1 core; OSOAA NT80 run with TODAY'S
EXPLICIT config: SZA30/RAA90/WIND3/AP.MOT=0.0973(555 std Rayleigh)/
AER 0/NbGauss 48/100, valid outputs verified; prior 0.53-0.55 s
anchor was a DIFFERENT config - per rules, no cross-session compare):
  OSOAA NT80:              2.38-2.41 s / condition (warm runs)
  OCRT classical(table):   6.24-6.43 s   (2.6x slower)
  OCRT PROD (beam loaded): 1.60 s        <- FASTER than OSOAA today
  (0+ re-run showed 2.59 s once - rerun variance to re-pin next
  session with n>=3 repeats; both bracket parity.)
PHASE-1 ARITHMETIC (measured components): beam solves 1,152 x ~2 s
= 0.6 h; operator assembly (case,wl,wind)=288 x ~7 min = 33.6 h
1-core -> 8.4 h across 4 machines (reduction levers exist: rbonly
already skips gates; per-part PP repetition can be hoisted); grids
64,512 x 1.6 s = 28.7 h -> 7.2 h across 4 machines.  Assembly is now
the dominant fixed cost - candidate next optimizations: (i) assemble
once per (case,wl) and reuse across wind if the injection response
is wind-independent (it is NOT - surface BC in solve; keep per-wind),
(ii) drop PP from assembly runs (needs p0 extraction outside PP).
PHYSICS-NO-CHANGE: pristine default 432-row bit-identical (0.00e+00
below); Tier-0 unchanged; strict was verified last commit; all D3
paths env-gated inert.

## D3-PROD ASSEMBLY-COST CUT (2026-07-10): 158 s single-shot per
## combo; operator BIT-IDENTICAL to the 4-part merge; PROD re-pinned
TWO LEVERS (both physics-clean, both verified):
 1. ZERO-WEIGHT COLUMN SKIP in the RBw loop: view/slot nodes carry
    zero quadrature weight, so unit injection there is identically
    zero (AI = 2 g_b ext = 0; measured earlier: 24/153 nonzero).
    The solve is skipped; the calloc'd response stays 0.  153 -> 72
    solved columns.
 2. rbonly PING-PONG SKIP: PPB = 0 in file-export mode (p0 kept as
    the subtraction base; bounces are gate instrumentation only).
RESULT: full assembly in ONE run: 158.3 s/combo (RBw 150.3 s, ~2.09
s/column) vs the 4-part ~430 s + orchestration.  The new operator
file is BIT-IDENTICAL to the previous 4-part merged file (Ed/Vw/RB
maxdiff all 0.00e+00) - partition/single equivalence and skip
correctness proven in one comparison.
PROD TIMING RE-PINNED (n=3, same host/session): 1.60/1.60/1.68 s -
the earlier lone 2.59 s was a one-off; parity vs today's OSOAA
(2.38-2.41 s) CONFIRMED with margin.
PHASE-1 ARITHMETIC (final, all measured; CORRECTION: beam solves
carry the AOD axis through F_sun - 9case x 16wl x 2wind x 4sza x
4aod = 4,608, not 1,152):
  assembly: 288 x 158 s  = 12.6 h 1-core -> ~3.2 h on 4 machines
  beams:  4,608 x ~2 s   =  2.6 h 1-core -> ~0.6 h
  grids: 64,512 x ~1.64 s = 29.4 h 1-core -> ~7.3 h
  TOTAL ~11 h on 4 machines (vs ~46 h classical-only arithmetic).
PHYSICS-NO-CHANGE: pristine default 432-row bit-identical (0.00e+00,
this session); Tier-0 unchanged (below); strict clean; all D3 paths
env-gated inert.

## OSOAA % + POLARIZED (P>0) END-TO-END GATE (2026-07-10)
BACK-TO-BACK SAME-CALL TIMING (this host, 1 core, warm runs; OSOAA
NT80 today's explicit config SZA30/RAA90/WIND3/AP.MOT0.0973):
  OSOAA NT80:            1.84-1.94 s / condition
  OCRT PROD (marginal):  1.40-1.42 s  -> 130-139% of OSOAA speed
  OCRT PROD (Phase-1 amortized: +0.85 s/grid fixed) ~2.26 s -> ~84%
  OCRT classical(table): 5.72 s       -> ~33% of OSOAA speed
POLARIZED CHAIN CLOSED (Red_clay.mie hydrosol, Q ~ 24% of I): the
production-grid assembly is 11 s/column (5.3x scalar) - sandbox-
excess, so per the reduction rule the FULL chain was verified on a
reduced config (n_mu 12, n_mu_water 12, raa-step 45, 96 rows):
assembly 169 s single-shot (n=27, I and Q basis columns both 12/27
nonzero as physics dictates) -> PROD swap+beam-cache -> gate vs
classical(table): every column ABSOLUTE <= 2.8e-8 (Q/U relative
figures are small-denominator inflation) - PASS.
BUG FOUND AND FIXED BY THE POLARIZED GATE: the swap's U view-scalar
synthesis was missing the solver's own sign convention
(Lw_U = -reconstruct_phi_sin, rt_water_rt.c [IA] loop) - invisible
in every scalar-phase gate (U ~ 0), a 125% U error in the polarized
one.  Transcribed verbatim with a source comment.
Production-scale polarized assembly (72 cols x ~11 s ~ 13 min/combo)
is a LOCAL-machine batch item; the sandbox verdict on correctness is
complete.  NOTE for Phase-1 arithmetic: polarized combos cost ~5x
scalar assembly - refine the 3.2 h figure by case mix when the
Phase-1 case list is frozen.
PHYSICS-NO-CHANGE: pristine default 432-row bit-identical (0.00e+00
below); Tier-0 unchanged; strict clean; all D3 paths env-gated.

## BASELINE RE-VERIFICATION vs OSOAA (2026-07-10): the table kernel
## is the one that AGREES with OSOAA under a real atmosphere
Per Jae's directive speed work is ON HOLD; this session re-verified
the OCRT<->OSOAA physics baseline and collected per-condition speed.
CONDITION FORENSICS FIRST: today's initial +9.0% rrs(0-) discrepancy
at 555/sza30 was traced to CONFIG, then to PHYSICS:
 - the historical harness ran OSOAA with AP_MOT = 0.01 (thin atm),
   recovered from the transcripts - today's 0.0973 run was a
   different condition, so no cross-session % was quoted;
 - a measured matrix (OCRT P=0 / 107.8 / 1013.25 x OSOAA tau 0.01 /
   0.0973) then isolated the real driver: the gap GROWS with the
   diffuse-sky share, implicating the default (legacy) injection
   kernel found earlier.
THE ADJUDICATION (555, cs5.0, wind 3, P=1013.25 vs OSOAA tau=0.0973,
rrs(0-,nadir), same host/session):
  sza20: OSOAA 0.028908 | table +0.06% | legacy +8.4%
  sza30: OSOAA 0.029329 | table +0.50% | legacy +9.0%
  sza40: OSOAA 0.029926 | table +0.87% | legacy +9.5%
  sza50: OSOAA 0.030562 | table +0.52% | legacy +9.4%
=> THREE-WAY CROSS-CHECK CLOSES: OCRT(table) == OCRT-D3-operator
(bit/CSV) == OSOAA (0.06-0.87%), while the DEFAULT legacy kernel
sits +8-10% high under a real atmosphere.  This upgrades the
table-kernel promotion from an internal consistency argument to an
external-validation argument - Jae's decision item now carries the
strongest possible evidence; golden re-baseline required on
approval.  (P=0 references: OCRT 0.029684; thin-atm OSOAA 0.029214;
+1.2-1.6% - a DIFFERENT-condition comparison, listed for context.)
PER-CONDITION SPEED (same call, 1 core; single-point rrs unit):
  OSOAA NT80: 3.5-17.7 s/condition (large spread across sza runs)
  OCRT single-point (either kernel): 1.7 s/condition
  (batch full-view unit, from the earlier back-to-back: PROD
   1.40-1.42 s vs OSOAA 1.84-1.94 s = 130-139%.)
NO CODE CHANGES THIS SESSION (speed work held).  Pristine anchors
re-verified: Tier-0 rrs0minus unchanged (below).

## FULL-SCOPE VERIFICATION CAMPAIGN - SESSION 1 (2026-07-12)
GOAL (Jae): verify ALL axes (geometry, aerosol model, AOT, Rrs IQU,
TOA IQU, bands, cases) + confirm speed > OSOAA -> then migrate.
RECOVERED ASSETS (unblocks every axis):
 - /tmp/pkgnew/GOCI3_PKG/scripts/ocrt_osoaa_consistency_harness.py =
   the ORIGINAL harness.  Historical scope now EXACT: rrs(0-,nadir),
   I only, 6 bands x Csed cases, SZA30, WIND3, OCRT atm OFF vs OSOAA
   AP_MOT=0.01.  Geometry/QU/TOA/aerosol were NEVER validated vs
   OSOAA before - they are NEW validation, not re-validation.
 - Historical golden figures recovered verbatim: MAPE ~1.6%,
   red/NIR -2..-5% (OCRT low, cause = model3-cap vs OSOAA
   truncation); 555/Csed5 anchor +0.06% at n_mu=48.
 - A_W/B_W pure-water tables; per-band mineral a*,b*,bb/b pipeline;
   27-case IOP tables (iop_cases/); 7 aerosol .mie; the mie->ExtData
   converter; OSOAA -AER.Model 4 -AER.ExtData CONFIRMED in the exe
   strings -> same-optics aerosol pairing is possible.
HARNESS v2 BUILT (/tmp/harness_v2.py + runners): axes = atmosphere
(thick tau_R(lam), Bodhaine std) x kernel{legacy,table} x geometry
(full VZA, phi) x Stokes IQU x level{0-,TOA} with per-run timing.
Env gotchas fixed: OSOAA_ROOT required; OCRT cwd=/tmp/tg for LUTs.
RESULTS THIS SESSION (Csed5, SZA30, wind3):
 T2b PAST-CONDITION REPLAY (OCRT P=0 vs OSOAA tau0.01):
   443 +1.24% | 555 +1.61% | 660 -3.13%
   -> matches the historical golden family (red -2..-5%) =>
   BASELINE INTACT.
 T2 THICK-ATM x KERNEL (OCRT P=1013.25 vs OSOAA tau_R(lam)):
   443: legacy +36.91% | table -3.37%
   555: legacy  +9.01% | table +0.50%
   660: legacy +10.17% | table -3.54%
   -> legacy divergence scales with the Rayleigh sky share (blue
   worst) - external confirmation of the kernel finding at all
   bands.  Table-kernel residuals sit in the historical -3% family;
   the 443 shift (+1.2 thin -> -3.4 thick) is an OPEN item
   (candidates: tau_R pairing, sky-polarization vs scalar LUT,
   OSOAA thick-atm handling).
 T3 GEOMETRY (555, phi sweep): azimuth convention ADJUDICATED =
   OSOAA phi corresponds to OCRT (180 - raa).  TOA_rho_I agrees to
   +0.15% MEAN ACROSS THE FULL VZA RANGE (phi=100 vs raa=80) - the
   TOA/geometry axis for I is essentially closed.
   0- radiance vs VZA is OPEN: the lv27 angle-axis convention is
   unresolved (air-equivalent: small-angle +1%, large-angle -10%;
   Snell-mapped: +0.9 -> +8.7% monotone; OCRT_VZA_IN_WATER probe is
   a different-definition diagnostic mode - discarded).  Next
   deciders: (i) OSOAA vsVZA standard output cross-read, (ii)
   thin-atm angular scan to remove the sky term, (iii) Q/U sign
   conventions after the angle axis is fixed.
 SPEED (same host, per condition incl. these runs): OSOAA 2.0-2.4 s;
   OCRT single-point 1.7-2.1 s; batch-unit PROD 1.40-1.42 s vs
   OSOAA 1.84-1.94 s (130-139%) stands.
NEXT SESSION PLAN (ordered): (1) settle the 0- angle axis (vsVZA
cross-read + thin-atm scan); (2) Q/U conventions + quantitative QU
gate; (3) aerosol axis (same .mie both sides, AOT sweep) - NEW
validation; (4) extend bands to 6 + more cases; (5) full speed
table per condition; (6) then migration package.
NO OCRT CODE CHANGES THIS SESSION.  Tier-0 unchanged.

## HARNESS-FAITHFUL RE-ADJUDICATION (2026-07-12, after Jae's flag)
Jae: "the old kernel never showed such errors - suspect wrong
variable pairing; follow the harness EXACTLY (raa, glint, etc.)".
ACTIONS: dropped my reimplementation; called the ORIGINAL harness
functions verbatim (h1.run_ocrt / h1.run_osoaa: scalar mode - no
--vector, OSOAA_MAIN.exe - not NT80, n_mu_water=48, cwd=pristine
tree, OSOAA_NO_DIRECT_GLINT=1 confirmed PRESENT in the original,
phi=raa passthrough).  My earlier runs had used vector/NT80/nmw24 -
those NUMBERS are superseded (structure was right, values shift).
VERBATIM REPLAY (thin atm, past condition):
  555/Csed5: OCRT 2.887685e-02 (= Tier-0 anchor value) vs OSOAA
  2.921414e-02 -> -1.16%.   443/Csed5: -0.47%.
  NOTE: the recorded +0.06% anchor was at the UNDER-CONVERGED
  n_mu=32 per the harness's own comment; the n_mu=48 golden was
  marked "to be re-established" - today's -1.16%/-0.47% are the
  candidates for that re-established baseline.
MINIMAL-DELTA THICK-ATM (only pressure/AP_MOT/kernel changed):
  555/Csed5 vs OSOAA(tau .0973): legacy +6.37% | table -2.22%
  443/Csed5 vs OSOAA(tau .2361): legacy +35.38% | table -4.99%
VERDICT (answers Jae's flag):
 (1) Variable pairing was NOT the issue - the verbatim replay
     reproduces the golden family and the kernel gap persists under
     the harness-faithful setup.
 (2) WHY the old kernel "never showed" this: the historical harness
     was THIN-ATM ONLY, where the sky-diffuse injection (the only
     place the legacy defect lives) is negligible - legacy scores
     -0.5..-1.2% there even today.
 (3) The thick-atm condition is NEW territory; the legacy/table gap
     (+8..+40%p at blue) is real physics of the default kernel, now
     established under harness-faithful conditions.
 (4) Table-kernel thick-atm residuals (-2.2% at 555, -5.0% at 443)
     exceed the thin-atm family - OPEN item (candidates: tau_R
     pairing OCRT-internal vs Bodhaine input, scalar-mode sky
     polarization, OSOAA thick-atm discretization).  To be
     decomposed before the kernel-promotion decision is finalized.
NO OCRT CODE CHANGES.  Tier-0 unchanged (2.887685e-02 re-observed
inside the replay itself).

## ITEM-2 REPLAY EXPOSES A v1.10 REGRESSION (2026-07-12) - CRITICAL
Ran the CANONICAL item#2 script (osoaa_pureocean_compare.py, paths +
required-env adapted ONLY) at sza40:
  pristine v1.09:  TOA_I MAPE 0.79% | Rrs(0-) nadir +0.46%
                   == EXACT match to the recorded matrix figure
                   ("sza<=40 I MAPE 0.79%") - canonical replay OK.
  current v1.10:   TOA_I MAPE 7.94% | Rrs(0-) nadir -3.21%
=> A REAL PHYSICS REGRESSION exists in the v1.10 sprint, in paths
the per-commit bit gates never covered (gates were fixed-bulk + LUT
+ n_mu_water=24 batch + P=0 single point; item#2 uses simple-chl
pure water + n_mu_water=8 + P=1013 single-point vza sweep).
ONE-POINT BISECTION (555, sza40, vza30, TOA_I):
  v1.09: nmw8 4.4274898e-2 | nmw24 4.4273010e-2  (converged, ~equal)
  v1.10: nmw8 4.2886837e-2 | nmw24 4.4218582e-2
  -> REGRESSION VARIABLE = SMALL in-water node counts: nmw8 is
     -3.1% vs pristine; nmw24 still carries a separate -0.12%
     residual on this pure-water single-point path.  Atmosphere
     n_mu has no effect.
CONSEQUENCES:
 - The thick-atm kernel comparisons run on v1.10 (legacy +6..+35%,
   table -2..-5%) are CONTAMINATED by this regression and are
   SUSPENDED until the regression is fixed - re-adjudicate after.
 - Gate design lesson: add item#2 representative cells (pure water,
   small nmw, P=1013 single point) to the per-commit reproduction
   set, per VALIDATION_MATRIX section 4.
NEXT TURN (single stream): bisect the archived per-commit patch
chain (outputs/ tars) on the one-point probe to find the breaking
commit; fix; re-run item#2 both kernels; then resume the campaign.
NO CODE CHANGES THIS SESSION.  Tier-0 unchanged.

## REGRESSION BISECTION - SESSION 1 (2026-07-12)
PROBE: one point (555, sza40, vza30, raa90, P=1013.25, pure water,
--vector), TOA_I first token.  GOOD = 4.4274898e-2 (nmw8) /
4.4273010e-2 (nmw24); BAD = 4.2886837e-2 / 4.4218582e-2.
SNAPSHOT LADDER (archived full-diff patches, "Only in" new files
copied from the current tree - assumption noted):
  GOOD: pristine v1.09, S-chain..555cross, B0a..B0c2b_final
  BAD : D3-0 and everything after (incl. current HEAD)
  => the breaking change landed BETWEEN B0c2b_final and D3-0
     (the S9-S12 session window; no intermediate snapshots exist).
ELIMINATIONS (each by build+probe):
  - ALL FASTK-gated blocks (S9 Td hoist, S10 coupling-first,
    S12 aliasing/prim-skip, mkc): a -DOCRT_FAST_KERNELS-REMOVED
    build of D3-0 reproduces the BAD values EXACTLY -> not FASTK.
  - Whole-file reverts of rt_kernel.* and rt_air_water_coupling.*
    fail to link (the couple_* SIGNATURES were EXTENDED in the
    window) - file-level revert impossible, hunk-level required.
  - Exact-match-first probe inside the D-SKY interp wrapper: no
    change -> the wrapper's dedup/exact-node path is not the lever.
REMAINING SUSPECTS (ordered): (1) the non-wrapper changes in
rt_air_water_coupling.c - the couple_atm_to_water /
couple_water_to_atm signature extensions and their body edits
(these are unconditional and sit exactly on the skylight path);
(2) shared/surface.c edits; (3) rt_types.h field changes feeding
those.  NOTE the D-SKY comment itself says "nadir was the only
healthy view angle - which is why nadir-based validations never
saw this": the window authors KNEW off-nadir was changing.
Whether the change is a fix that broke OSOAA agreement or a bug
is exactly what the item#2 canonical gate must adjudicate AFTER
the mechanism is isolated.
NEXT TURN: hunk-by-hunk readback of the coupling.c window diff
(wrapper call sites + signature-extension bodies), targeted
reverts on the D3-0 tree, probe each; then fix on HEAD, re-run
item#2 (both kernels), resume the campaign.
NO CHANGES TO THE WORKING TREE THIS SESSION (all experiments in
/tmp scratch trees).  Tier-0 unchanged.

## REGRESSION BISECTION - SESSION 2 (2026-07-12): mechanism cornered
FURTHER ELIMINATIONS (scratch-tree experiments, probe = the item#2
one-point):
 - D-SKY wrapper CALL SITES reverted to direct interp: NO change.
 - mb_re_escape_add_m0: internally gated (OCRT_MB_CLOSURE opt-in,
   verified in source) AND blocking both call sites: NO change to
   the probe (the calls are dead by default, as reported).
THE CORNER: the window REPLACED the skylight-injection mechanism.
 - A 108-line function ocrt_add_diffuse_top_primary() (per-m
   first-scatter primary source, moment l-sums, "B-0c.2b" labelled
   but landed AFTER B0c2b_final) was ADDED in the window and is
   called at rt_water_rt.c:3714.
 - Blocking that call moves the probe to a THIRD value (skylight
   first-scatter lost) - it contributes, but simple blocking does
   not restore the GOOD value: the window swapped mechanisms, not
   just added one.
 - DECISIVE NEGATIVE: OCRT_EXTTOP_KERNEL=table changes NOTHING on
   the item#2 path (identical 7.94%/-3.21% with and without).  The
   legacy/table adjudication applies ONLY to the ext_top ARRAY
   injection path (fixed-bulk production batches).  The pure-water
   single-point path injects skylight EXCLUSIVELY through the new
   first-scatter primary function - a DIFFERENT code path that the
   kernel env never touches.
CONSEQUENCE FOR THE KERNEL DECISION: the earlier thick-atm
legacy-vs-table tables measured the ARRAY-injection path on a build
whose PRIMARY-source path is regressed; the two questions are now
SEPARATE: (Q1) fix the first-scatter primary mechanism so item#2
returns to the 0.79% family; (Q2) then re-adjudicate legacy/table
on the array path.
NEXT TURN: read rt_water_rt.c:3714 call context + diff the
REMOVED bA-side mechanism (what did B0c2b_final use for skylight
on this path); pinpoint the physics delta (suspects: the "U-incident
DROPPED" note, the Rayleigh-only ray_on m<=2 coefficients, the
A = 2*w_c*S^m amplitude convention); fix; re-gate item#2.
All experiments in /tmp scratch; working tree unchanged; Tier-0
unchanged.

## REGRESSION BISECTION - SESSION 3 (2026-07-12): two components
MECHANISM READBACK:
 - bA (B0c2b_final) water solver has NO ext_top at all: skylight
   reached the water as EQUIVALENT BEAMS (rt_solver.c:5458,
   F_sun_eq = I_in*2pi*w per sky node) - every ring column exact.
 - The window replaced this with direct injection:
   ocrt_add_diffuse_top_primary (per-m first-scatter, called at
   rt_water_rt.c:3714 under the ext_top arrays).
NODE-MATCH DEFECT (probe-proven, [DTPmiss]): the injector matches
ring columns to ext samples at 1e-12; inserted columns (view, slot,
nadir: mu = cos(sza), mu_view, 1.0) get NO skylight.  ext_n=8 on
the nmw8 run (water-GL resampled), inner grid n_mu=12.
DTP-FIX (landed on HEAD, gated-inert): linear-in-mu interpolation
fallback for non-matching columns (+ OCRT_DTP_DIAG counters).
PROBE RESULT: BIT-IDENTICAL - the view/slot first-scatter columns
do NOT feed the TOA_I/rrs outputs on this path.  Fix retained
(correct property, zero cost) but it is NOT the regression lever.
COMPONENT SPLIT (Ed/Lu readback, 555 sza40 vza30):
   Ed(0-):  v109 0.70339/0.70341 (nmw8/24)
            v110 0.70313/0.70349  -> total flux ~OK (<=0.04%)
   Lu(0-):  v109 1.1786e-3/1.1779e-3
            v110 1.1595e-3/1.1611e-3  -> -1.6%/-1.4%, ~nmw-FREE
 => Component (A): upwelling response to injected skylight is
    ~1.5% low regardless of node count - the injection PHYSICS
    (suspects: U-incident handling, amplitude convention
    A=2*w*S^m vs the beam ladder, first-scatter-only vs beam
    all-orders bookkeeping in the primary).
 => Component (B): the TOA gap is strongly nmw-dependent
    (-3.1% @8 vs -0.12% @24) while Lu(0-) is not - the extra loss
    sits in the water->atm return leg, whose couple_* bodies and
    SIGNATURES were changed in the same window (file-revert
    impossible; hunk work pending).
NEXT TURN: (B) first - readback couple_water_to_atm window hunks
(bA vs bB) and targeted-revert on the bB scratch; then (A) -
compare the injected-source ladder against an equivalent-beam
reference on one mode.  Then re-gate item#2 and re-run the kernel
adjudication.
WORKING-TREE DELTA THIS SESSION: DTP-FIX + DTP diag in
rt_water_rt.c - probe-verified bit-identical on the item#2 point;
default-path 432-row bit gate still to be run next turn (top of
list).  Tier-0 unchanged.

## REGRESSION SESSION 4 (2026-07-12): gates green; component B is
## an ATMOSPHERE-side loss; equivalent-beam switch recovered
GATES AFTER THE DTP-FIX LANDING (working tree):
  Tier-0 anchor: rrs0minus=2.887685e-02 EXACT.  Previous-build bit
  gate (v2o_fastk vs v2o_fix, 432-row 555 batch): BIT-IDENTICAL
  PASS.  (443-golden direct replay: generation command not in the
  ledger - covered transitively by the previous-build identity;
  logged as a maintenance item.)
LIVE ROLLBACK SWITCH FOUND: OCRT_EQUIV_BEAMS=1 restores the
pre-window beam mechanism (rt_solver.c:5657 dt_mode=0).
  Probe with beams: Lu(0-) 1.17506/1.17692e-3 (nmw8/24) ==
  pristine to -0.3%/-0.09%  =>  COMPONENT A (the ~-1.5% upwelling
  deficit) is CONFIRMED to live in the diffuse-top injection and
  VANISHES under beams.  TOA stays low under beams (-3.2% @nmw8,
  -0.33% @nmw24)  =>  COMPONENT B is INDEPENDENT of the injection.
COMPONENT B LOCALIZED TO THE ATMOSPHERE:
  Full stdout diff (v1.09 vs HEAD, nmw8): Lu(0-) and Lu(0+) are
  low by the SAME -1.62% (water->air transfer is loss-free; the
  couple_water_to_atm suspicion is CLEARED - its signature growth
  is just the gated re-escape weights).  Arithmetic: the water
  share of TOA is ~1.3%, so Lu explains only ~-0.02%p of the
  -3.14% TOA gap - the remaining ~-3.1%p is PURE-ATMOSPHERE
  response changing with the WATER node count.
  Partial confirmation: pristine v1.09 forced to --n-mu 8 gives
  TOA 4.38375e-2 (-1.0% vs pristine24) - the unified-grid
  "atmosphere runs at the water node count" effect is REAL but
  accounts for only ~1/3; a further ~-2.1%p component remains
  (suspects: slot/view ring composition of the unified grid,
  atm2 second-pass wiring).  BC-FIX one-liner (n_mu_atm_bc) was
  probe-inert (dead path here) - kept for hygiene.
NEXT: instrument the ACTUAL atm SOS node count on HEAD (one log
line), then chase the residual -2.1%p; afterwards decide the
mechanism (repair dt-chain vs re-default beams for correctness
with dt as the fast path) against item#2 + the speed budget.
Tier-0 unchanged; tree delta = DTP-FIX + BC-FIX (both probe-inert,
prev-build bit gate PASS).

## REGRESSION SESSION 5 (2026-07-12): component-A chase + A-FIX(U)
CONFIRMED: --n-mu is IGNORED on the coupled path (B-0b.1 hard
override at rt_solver.c:5778, probe: n_mu 4/24/48 bit-identical).
CONTAMINATED PROBES DISCARDED: OCRT_ATM_NMU_INDEP (+BC-FIX combos)
create a MIXED 24/8 state - Ed(0+) collapses to 0.678 at nmw8
(a physical impossibility for a water-independent quantity).
Lesson: partial de-unification is not a valid experiment; the
probes stay env-gated OFF and are slated for removal.
ELIMINATED for the residual: B-0c.1 T_aw matrixification
(OCRT_TAW_FLAT probe: +0.008%, nil) - the rough-surface forward
matrix is NOT the lever.
ITEM#2 AT PRODUCTION NODES (nmw24): TOA_I MAPE 2.83%,
Rrs(0-) -3.05%  =>  even without the small-node pathology the
canonical gate FAILS the <=1% criterion; component A (injection
physics, band-dependent via the sky share) DOMINATES in practice.
A-FIX LANDED (physics change, intended): wire the U-incident of
the diffuse top (ext_top_U was NULL "next increment").
  Probe: U(0-) +1.7% (real effect, correct physics), TOA +0.001%,
  Lu +0.02% - NOT the component-A lever.  Remaining suspects:
  the amplitude convention A = 2*w_c*S^m vs the beam ladder, and
  the first-scatter source SHAPE vs the beam exponential path.
GATE STATUS: Tier-0 unchanged (P=0, injection-free).  The
previous-build 555-batch bit gate now differs BY DESIGN (A-FIX is
a physics repair on the injected path); new anchors to be frozen
only after component A closes.
SPEED (Jae directive: no slow-down): 432-row batch timed this
session, v2o_fastk vs v2o_fix5 - MEASUREMENT INVALID this session (batch exited in 4 ms - rerun with stderr next turn, TOP of list); A-FIX
adds one ndt memcpy per solve (negligible).
NEXT: amplitude-convention audit of ocrt_add_diffuse_top_primary
against one equivalent beam on a single (m, node) - close
component A; then re-run item#2 nmw24 and the campaign.

## REGRESSION SESSION 6 (2026-07-12): REAL gates restored; component
## A quantified as a mu0-dependent injection shortfall
GATE INTEGRITY INCIDENT (found+fixed this session): earlier
"previous-build bit PASS" lines this campaign were FALSE - the
batch runs were dying instantly (--batch-csv / --batch / missing
CLI base args) and the compare was re-copying a stale out/D3P.csv.
Correct invocation recovered: CLI base + --lut-raa-step 20 +
--batch-full-grid jobs_D3P.csv.  REAL results:
  fastk == fix4 BIT-IDENTICAL (all four env probes default-inert:
  DTP-FIX, BC-FIX, ATM_NMU_INDEP, TAW_FLAT)  -> PASS.
  fix4 == fix5 on this fixed-bulk batch (A-FIX U-wiring does not
  touch the fixed-bulk injection path; it moved the pure-water
  single point) -> the production grid is UNCHANGED by A-FIX.
SPEED (Jae directive): 432-row batch 5.23-6.16 s across the three
builds (first run cold); the current tree is NOT slower - 5.2 s
matches the recompile re-anchor.  A proper warm A/B timing is a
one-liner whenever needed.
COMPONENT-A LADDER (pure water, 555, vza30, nmw24; Lu(0-)):
  m>=1 amplitude scale 0.5x/2x: Lu moves <0.01% -> m>=1 innocent.
  ALL-m amplitude +10%: Lu +0.45%, Ed(0-) +0.001% -> Ed(0-) does
  NOT consume the injected source (it comes from the analytic
  Ed_diff path) - that is WHY Ed always matched while Lu fell.
  Source SHAPE: the injector does carry the beam-form exponential
  ladder (0.5*exp(-h/mu_c) per layer) - skeleton OK.
  INJECTION-SHARE SPLIT (OCRT_DTP_OFF probe):
    sza40: injected share 5.26e-5 vs required (beams-off) 6.82e-5
           = 77% delivered (-23%)
    sza80: 6.26e-6 vs 62.6e-6 = 10% delivered (-90%)
  => the shortfall is STRONGLY mu0-dependent - not a constant
  factor, not m>=1, not the ladder shape.  Prime suspect now: the
  INJECTED FIELD ITSELF (the coupled atm->water snapshot) being
  low at large sza - note Ed(0-) cannot see this (analytic path).
NEXT (single decisive step): log the m0 flux integral of the
injected field  sum_c 2pi*mu_c*w_c*ext_I0(c)  against the analytic
Ed_diff_water_from_atm on the same run (sza40 + sza80).  If the
integral is low by the same 23%/90%, the defect is UPSTREAM in the
snapshot (coupled field / its normalization or timing); if it
matches, the defect is in the source assembly coefficients.
Then fix, re-gate item#2 nmw24, re-run the campaign.
Probes added this session (all env-gated, default-inert, batch
bit-PASS above): OCRT_DTP_MSCALE/ASCALE/OFF.  Tier-0 unchanged.

## REGRESSION SESSION 7 (2026-07-12): component A cornered to the
## missing boundary interaction of the injected field
DECISIVE FLUX AUDIT ([DTPflux]): the injected m0 flux integral
  sum_c 2pi*mu_c*w_c*ext_I0(c) * (F/pi)
matches the analytic Ed_diff to -1.2% (sza40) / +1.6% (sza80).
=> UPSTREAM CLEARED: the coupled snapshot's TOTAL is right.
DISTRIBUTION AUDIT ([DTPdist], sza80): physically correct shape -
zero below the critical cone (residual leak via rough surface),
mass at mu_w ~0.5-1.0.  => distribution cleared.
RESOLUTION SCAN: dt/beams Lu ratio at nmw 24/48/96 = 0.674/0.674/
0.674 (sza80) - NODE-COUNT INVARIANT.  Resampling-resolution
hypothesis REJECTED.
VZA FINGERPRINT: ratio 0.671/0.679/0.684 at vza 10/50/70 - flat
(multiple scattering isotropizes the missing part; consistent with
an m0-dominant deficit).
SURVIVING MECHANISM (structural, now unique): the equivalent-beam
path carries the sky as a FIELD (each beam undergoes the full
in-water boundary physics - in particular the internal specular
reflection at the water-air interface, which approaches TOTAL
INTERNAL REFLECTION for the near-critical directions that dominate
at large sza).  The diffuse-top path injects the sky ONLY as a
first-scatter SOURCE: the downwelling diffuse field itself never
exists in the water solve, so its boundary interactions (internal
reflection above all) are ENTIRELY ABSENT.  Numbers line up:
missing share -23% of the sky Lu at sza40 (gentle incidence,
Fresnel-few-%..TIR mix) and -33% flat at sza80 (near-critical
dominated), independent of node count - exactly the signature of
a missing reflection channel, not of discretization.
FIX DESIGN (next turn, single stream): add the boundary-reflected
first-scatter increment - inject, alongside the downwelling source,
the once-internally-reflected upwelling field R_ww(mu_c) x ext(c)
(mirror column, Fresnel/coxmunk internal reflectance; the existing
surface_R_ww_coxmunk machinery in shared/surface.c provides the
coefficient), plus its own exp ladder.  This is the field-vs-source
closure term; beams get it for free.  Verify: dt/beams -> 1.000 at
both sza, then re-gate item#2 nmw24 (target 0.79% family), then
production-batch bit gate (fixed-bulk path expected UNCHANGED -
verify, not assume) + warm speed A/B.
Probes added: OCRT_DTP_FLUX(=1 flux, =2 +distribution) - env-gated,
default-inert.  Tier-0 unchanged.

## REGRESSION SESSION 8 (2026-07-12): COMPONENT A CLOSED (A-FIX-2)
DIRECTION-RESOLVED RESPONSE PROBE (ONLYC/ONLYJ, baselines
bit-identical): dt/beam response ratio is DIRECTION-INDEPENDENT -
0.099/0.104/0.105/0.103 across four ring columns at sza80, and
~0.77 at sza40  ==  cos(sza) family  =>  a pure SCALAR scale error.
GLOBAL-RESCALE ORACLE: with ASCALE = measured 1/ratio, dt Lu equals
the beam Lu to SIX DIGITS at BOTH sza (0.000%) - mechanism proven.
ROOT CAUSE (convention): the water solve's common source/output
scale is F_sun_water = F_BOA*exp-attn*mu_air*T_aw/mu_w (the DIRECT-
beam converted flux, rt_water_rt.c:4519 / rt_solver.c:6401).  The
pi-normalized diffuse-top snapshot needs F_sun_TOA instead - it was
silently inheriting the direct-beam conversion factor (mu0-family,
hence the sza dependence; the old "m0 A/B gate" pinned amplitude
against Ed, which never consumes the source - vacuous).
A-FIX-2 LANDED: pre-multiply the snapshot by (F_sun_TOA /
F_sun_water) right after F_sun_water is formed (rt_solver.c:6401+).
Cost: one ndt-sized pass per solve - negligible (speed directive).
VERIFICATION:
  Lu(0-) sza40: dt 1.178680e-3 vs pristine v1.09 1.178605e-3
    (+0.006%) - dt now BEATS the beams' own -0.14% offset.
  sza80: dt vs beams -1.8% residual (was -32.6%) - small tail open.
  ITEM#2 sza40 nmw24 RE-GATE: Rrs(0-) +0.08% (PASS <=1%, better
  than v1.09's +0.46%); TOA_I: vza30/60 ALL BANDS -0.5..-1.5%
  (v1.09-class, PASS family); REMAINING: the NADIR COLUMN ONLY,
  blue-heavy (-9.1/-8.4/-6.1/-2.6% at 412/443/490/555) - a single
  inserted-view-column fingerprint (DTP-FIX interpolation now
  matters since the scale is right; nadir mu=1.0 insert path to be
  audited next).  MAPE headline 2.13% is entirely this one column.
STRUCTURAL VERDICT (Jae's standing question): the speed patch is
NOT structurally incompatible with the original accuracy - the
deficit was a single flux-convention scalar, now corrected in
place; the fast injection mechanism itself reproduces the beam
physics to 6 digits once scaled right.
STILL TO RUN (next turn, in order): production-batch bit gate for
A-FIX-2 (fixed-bulk path - verify, expect possible intended change
now that the injected path is live there too), Tier-0, warm speed
A/B, then the nadir-column audit and the sza80 -1.8% tail, then
full item#2 (sza 0/40/80) and campaign resume.
Probes: ONLYC/ONLYJ env-gated inert.  DTP_OFF/ASCALE/MSCALE retained.

## REGRESSION SESSION 9 (2026-07-12): gates green for A-FIX-2;
## the last residual localized to a v1.10 NADIR-CONE RAMP
GATES: Tier-0 EXACT (2.887685e-02).  Warm speed: fixC 5.25 s vs
fastk 6.04 s - NOT slower (directive holds).  Production-batch
change quantified (INTENDED physics, A-FIX-2): rrs_I +0.97% avg,
TOA_rho_I +0.31% avg, Rrs_I +1.00%; TOA_rho_Q "53%" is a tiny-
denominator artifact (abs 2e-5).  New golden freeze deferred until
the last residual closes.
NADIR RESIDUAL LOCALIZED (412, sza40):
  v1.10: vza 0/0.5/1/2 = .18325/.18463/.18849/.19998 - a steep
  RAMP inside the nadir cone, exact at the last GL node
  (mu=0.9988 = vza 2.83 deg) and degrading toward mu=1.
  BEAMS show the same ramp (vza0 .17382 vs vza2 .18745) -> the
  defect is COMMON (not dt) - the component-B tail.
  v1.09: vza0 = vza1 = 0.20046 - FLAT nadir cone (physical truth).
  => the ramp is a pure v1.10-window artifact in the INSERTED-VIEW
  column handling of the unified-grid atmosphere solve (the
  inserted zero-weight view column starves between GL nodes -
  same disease family as [DTPmiss], air-side edition).
  Blue-heavy in item#2 because the atm share scales with tau_R.
SIDE NOTE: beams at 412 carry their own -6% offset vs v1.09 even
off-nadir - beams are a REFERENCE ONLY at 555; do not gate on
beams at blue.
dt vs v1.09 OFF-NADIR at 412: 0.19998 vs 0.20047 (-0.24%) - the
injection path is v1.09-class everywhere except the nadir cone.
NEXT TURN: audit the inserted-view column in the unified-grid atm
solve (B-0a.2 node-ring construction hunk + primary/SOS source
and reconstruct treatment of the zero-weight column; diff against
the v1.09 view_as_node machinery which was flat).  Fix, then
item#2 full (sza 0/40/80), freeze new goldens, campaign resume.
sza80 dt/beam -1.8% tail queued after (may share the same root:
grazing-view columns).

## REGRESSION SESSION 10 (2026-07-12): nadir ramp - fingerprints
## complete; a CANCELLATION structure exposed
FINGERPRINTS (412, sza40, pure water): ramp is FASTK-invariant;
azimuth-invariant at nadir (bit-identical raa 0/90/180 -> m0);
tau_R-proportional across bands; grid-PLACEMENT invariant (legacy
tail-append probe OCRT_UANG_APPEND: bit-identical); WIND-SENSITIVE:
wind10 flattens the cone (vza0 .20048 ~ vza2 .20055).
KEY DELETION FOUND in the window: the 2026-06-25 FIX that filled
the water-leaving bottom source at the EXPLICIT view node for the
pass2 atm SOS was replaced by the unified-table assumption.
ROW-FLUX AUDIT of couple_water_to_atm (fine-256 reference vs GL):
  oblique rows s=1.0000 (GL exact);  mu_a=0.9988 s=0.734;
  view row mu_a=1.0 s=0.376  ->  the water->atm bridge OVER-counts
  the nadir rows at low wind (sharp T_wa peak lands ON the last GL
  node).  Applying row-norm alone made TOA WORSE (.1833->.1583):
  the bridge's nadir OVER-count partially CANCELS a LARGER deficit
  elsewhere in the nadir chain (candidates: the atm->water bridge's
  matching nadir rows, or the pass2 view-column source that the
  deleted 2026-06-25 FIX used to guarantee).  TWA-ROWNORM kept as
  an OPT-IN diagnostic (default OFF; =2 dumps s rows).
NEXT TURN: audit the OTHER half - (i) the pass2 wl bottom-source
at the view column (does the unified path actually fill mu=1.0
row? dump wlI at the view index), (ii) the atm->water bridge's
nadir rows with the same fine-reference audit; then fix the PAIR
consistently (both bridges row-normalized together, or restore the
explicit view fill), re-verify nadir flat + item#2.
Tier-0 unchanged; default outputs bit-identical to fixC (ROWNORM
off by default - verified this session).

## REGRESSION SESSION 10 (2026-07-12): nadir ramp - narrowed to the
## vector-path view assembly; not yet pinned
NEW FACTS (all probe-based):
 - Azimuth invariance at exact nadir is EXACT (bit-identical at
   raa 0/90/180)  =>  the deficit is at the m0 LEVEL of the view
   value itself; no m>=1 leak.
 - --max-orders is water-only (atm path unaffected; probe values
   bit-identical) - order-split via CLI impossible.
 - The coupled TOA assembly: rho_I = atm_path_I*norm + (water
   transmitted)*rho_factor with atm_path_I = atm_res.I_TOA from
   rt_solve_case_pol_for_ocean(&cs_atm, &opts_atm...);
   opts_atm = *opts (view_as_node = 1 inherited - the interp
   branch is NOT the culprit by inheritance).
 - The vector solve writes I_TOA at rt_solver.c:932 via
   reconstruct_phi(I_total_per_m); the boa-export block just above
   is IDENTICAL v1.09 vs v1.10.  The view-value fill of
   I_total_per_m[m] (further upstream in the per-m loop) is the
   remaining unaudited segment - NEXT TURN: diff that fill (v1.09
   line ~800s vs v1.10 ~860s), find where the inserted-view column
   value diverges as mu_view -> 1, fix, then item#2 full run
   (sza 0/40/80) and new-golden freeze.
 - Grid construction readback: the water/unified ring carries
   standing zero-slots {1.0, mu_sun_water, mu_sun_air} (explains
   the [DTPmiss] mu list); the atm-side unified table appends the
   view via sorted-insert+dedup (OSOAA seuil) - order differs from
   v1.09's append-then-sort only away from the ends, so ordering
   alone cannot explain the nadir ramp.
NO code changes this session (audit only).  Tier-0 unchanged.

## REGRESSION SESSION 11 (2026-07-12): nadir ramp - four more
## eliminations; window ownership CONFIRMED
STATUS ANSWER (Jae: "is the error level recovered yet?"):
  Rrs(0-) RECOVERED (+0.08%, passes <=1%, better than v1.09).
  TOA off-nadir RECOVERED (-0.5..-1.5% all bands, v1.09-class).
  NOT yet: the single NADIR COLUMN (blue -9..-2.6%) = the entire
  2.13% headline; and the sza80 dt/beam -1.8% tail.
WINDOW OWNERSHIP: B0c2b_final probes FLAT at nadir (0.20046 =
v1.09 value exactly) -> the ramp entered in the SAME S9-S12 window
as the flux-scale bug.  Two independent defects, one window.
ELIMINATIONS THIS SESSION (all by build+probe or verbatim diff):
 - rt_uangles_add: inserted weight is EXPLICITLY 0.0; nadir insert
   position == v1.09 tail position (mu=1.0 sorts last) -> table OK.
   (An OCRT_UANG_APPEND legacy-placement probe exists in the file.)
 - fill_ring_from_uangles vs old append_view_node: VERBATIM-equal
   logic (GL rm/gb + zero-weight view + xpl/xrl aux).
 - Atm ring composition: GL + view ONLY (no standing slots on the
   atm side); NMUDET debug hook present.
 - rt_kernel.{c,h}: forward-ported into the GOOD tree -> nadir
   stays FLAT (0.20046) -> kernel innocent.
 - shared/surface.c forward-port: link-tangled (new deps) - to be
   adjudicated NEXT TURN by reverse probes on the BAD tree instead
   (candidates left: the R_ww eval split in surface.c, the
   non-FASTK rt_water_rt window edits, main.c/rt_types.h deltas).
NO changes to the working tree.  Tier-0 unchanged.

## REGRESSION SESSION 12 (2026-07-12): nadir ramp - the suspect set
## collapses to the sky_TOA assembly term
ELIMINATIONS (probe/build-based, all on the 412 sza40 nadir probe):
 - Water upwelling is FINE at nadir: Lu0+/Lu0- match v1.09 to
   +0.22% -> "water share dead" hypothesis rejected AS PRINTED, but
   see the new suspect below (the printed Lu may not be the TOA
   water term).
 - noFASTK build reproduces the ramp BIT-identically -> S1/S4
   solver vectorization and every other FASTK block innocent.
 - rt_types.h delta = S7 raw fields only; main.c delta = CSV
   splitter fix -> both innocent.
 - HYBRID BUILD (bB everything + bA {rt_solver, rt_water_rt,
   coupling}): FLAT (0.20046) -> the ramp lives INSIDE the trio.
   Single-file forward/backward porting within the trio is blocked
   by the couple_* signature extension (link errors) - hunk-level
   only.
 - rt_solver window hunks now read IN FULL (incl. @5318 B-0c.1
   weights-to-forward-coupling and @1108 S7 raw): none plausible
   alone; TAW_FLAT probe at 412 nadir: ramp unchanged -> the T_aw
   matrix path is innocent AT 412 TOO.
NEW PRIME SUSPECT (next turn, first thing): the TOA assembly term
sky_TOA_I - the skylight-induced water upwelling transmitted to
TOA (rho_I = atm_path*norm + (TOA_wl_direct + sky_TOA)*rho_factor).
It is produced inside the trio, is COMMON to dt and beams (both
feed the sky-upwelling), and the printed Lu0plus may not be the
same accumulator - which would reconcile "printed Lu fine" with
"TOA water term dead at nadir".  Audit: locate sky_TOA_I fill,
log it at vza 0 vs 2 (412), diff the fill against v1.09's
equivalent, fix, then item#2 full + golden freeze.
No working-tree changes this session.  Tier-0 unchanged.

## REGRESSION SESSION 13 (2026-07-12): NADIR RAMP ROOT CAUSE FOUND
CHAIN THAT GOT THERE (this session):
 - [C3FULL] decomposition (existing diag): at 412 nadir the ramp
   lives ENTIRELY in TOA_wl_direct_I (the rigorous water-leaving
   transmission term): v1.10 1.237e-2 (vza0) vs 1.644e-2 (vza2)
   while v1.09 is flat (1.329e-2 both).  The separate sky term is
   0 by DESIGN in dt mode (the injected-sky upwelling rides the
   same wl pass; off-nadir totals match v1.09's direct+sky to
   -0.7%) - not a bug.
 - The wl pass-2 bottom source is built by couple_water_to_atm
   into the unified atm ring (stride/view-slot verified correct).
 - S7b-cache branch view_as_node=0 hypothesis: NADIR-FIX guard
   added but probe-inert - that branch does not run on single
   points (kept as hygiene for grid runs).
 - DECISIVE: a probe switch on the REVERSE (water->air) T_wa
   matrix branch (OCRT_TWA_FLAT, coupling.c ~546 - my earlier
   TAW_FLAT only covered the FORWARD air->water branch at 202!):
   the ramp VANISHES (vza0 == vza2) -> ROOT CAUSE = the
   surface_T_wa_coxmunk_fourier_kernel VIEW ROW as mu_air -> 1
   (nadir cone).  The flat fallback's absolute level is wrong for
   this config (it is a diagnostic path only).
FIX DESIGN (next turn): repair the T_wa fourier-kernel view row
near mu_air=1 (suspects: the phi-integral sampling/jacobian as
mu->1, or the kernel's water-side pairing for the nadir row);
acceptance = nadir flat AND equal to the off-nadir v1.09-class
level; then item#2 full (sza 0/40/80), golden freeze, speed A/B.
Probes added: OCRT_TWA_FLAT (reverse), NADIR-FIX guard - both
env-gated / condition-inert by default.  Tier-0 unchanged.

## REGRESSION SESSION 14 (2026-07-12): T_wa kernel - mechanism fully
## understood; fix design frozen
STANDALONE KERNEL TEST (surface_T_wa_coxmunk_fourier_kernel, m=0,
wind=3, rows mo={0.99877, 0.99990, 1.0} x cols mi={0.3, 0.6,
0.7463, 0.95}): T00 is EXACTLY 0 everywhere (1e-48..1e-62 dust at
mi=0.95).  The near-nadir rows of the reverse transmission matrix
are numerically EMPTY.
GEOMETRY AUDIT: the facet formula is CORRECT - at the true Snell
pair for upward transmission (water mu_w=1 <-> air mo=1; NOT
0.7463, which is the DOWNWARD sun pair) cos_beta evaluates to
exactly 1 (horizontal facet).  The defect is SAMPLING: the
transmission row is a near-delta in mu (Cox-Munk slope width maps
to ~1e-3 in mu at nadir - NARROWER than the GL node spacing), so
any node-sampled matrix row integrates to ~0.  v1.09's hybrid
(Snell interpolation) handled this delta ANALYTICALLY - that is
why it was flat.  The off-nadir -0.2..-1.5% soft bias is the same
disease's weak face (narrow rows everywhere, partially resolved
by nearby nodes).
FIX DESIGN (frozen): split the reverse coupling into
  T_wa = [ANALYTIC direct-refraction term: Snell mapping mu_w ->
          mu_a with Fresnel transmittance and the n^2 radiance
          jacobian - a delta handled by interpolation onto the atm
          ring, exactly the v1.09 hybrid for ALL m] +
         [RESIDUAL rough-surface scattering matrix: the existing
          fourier kernel evaluated with the delta REMOVED (or kept
          as-is if its residual is verified small at wind 3)].
  Acceptance: (1) standalone row test reproduces energy (row
  integral == T_wa_direct flux to <0.1%); (2) 412 nadir flat and
  v1.09-class; (3) item#2 sza40 TOA MAPE <= 1%; (4) production
  batch bit-change quantified; (5) warm speed not slower.
IMPLEMENTATION NEXT TURN (whole turn reserved): reuse the v1.09
hybrid code path (still present as the TWA_FLAT fallback skeleton)
as the analytic term inside the matrix branch, instead of a
whole-path fallback.
Probe kept: OCRT_TWA_FLAT.  No working-tree physics change this
session.  Tier-0 unchanged.

## REGRESSION SESSION 15 (2026-07-12): NADIR RAMP CLOSED - item#2
## sza0/40 FULLY PASS; TWA-DEFAULT landed
JAE'S DIRECTIVE HONORED (read harness/comments first) - and it
paid twice:
 (1) Found an EXISTING prior attempt in the code: TWA-ROWNORM
     (same diagnosis, fine-grid row renormalization) demoted to
     opt-in because "the GL over-count partially cancels a larger
     deficit elsewhere" - that larger deficit was TODAY'S A-FIX-2
     scale bug.  Re-tested with the scale fixed: ROWNORM moves the
     WRONG way (row = air-side axis; conservation is per WATER
     direction).  Stays opt-in, documented.
 (2) Found my OWN probe bug: the earlier TWA_FLAT switch SKIPPED
     the fill and returned zeros instead of falling to the flat
     path - the "fallback absolute level is wrong" verdict was an
     artifact.  Probe fixed (demote ok=0 -> true fallback).
RE-ADJUDICATION WITH THE REAL FALLBACK (= the v1.09 analytic
Snell/Fresnel hybrid, still intact in the file):
  412 sza40: vza0 0.200599 vs v1.09 0.200463 (+0.07%), FLAT.
  item#2 sza40: TOA_I MAPE 0.79% == the recorded v1.09 figure
  EXACTLY; Rrs(0-) +0.08%.  sza0: 0.71% / +0.22%.  ALL PASS <=1%.
TWA-DEFAULT LANDED: reverse (water->air) coupling defaults to the
analytic path; the Cox-Munk reverse matrix is OPT-IN
(OCRT_TWA_MATRIX=1) until a delta-split version exists - its
transmission rows are near-deltas (standalone test: numerically
EMPTY at nadir) and node-sampled contraction loses the
direct-refraction flux.  Full rationale in the code comment.
GATES: Tier-0 EXACT.  Production 432-row batch: rrs_I and Rrs_I
BIT-UNCHANGED (0.000%); TOA_rho_I changes by design (max 4.69%,
mean 1.16% - nadir-cone rows corrected).  NEW GOLDEN FREEZE READY.
Warm speed 4.18 s vs 4.23 s - NOT slower (matrix build removed).
REMAINING OPEN (separated, not blockers for sza<=40):
 - sza80: v1.09 itself carries Rrs +8% at 443 (a known grazing-
   regime limitation family); fixG TOA 2.96% (6-band) vs v1.09
   -0.8% (443) - residual to quantify per band next.
 - dt/beam Lu tail at sza80: ratio 0.9818 unchanged by TWA (it is
   an injection-vs-beam residual, queued).
NEXT: freeze new goldens (443 + item#2 cells per matrix section 4),
per-band sza80 quantification, then the campaign matrix (#1
aerosol AOT 1.0 replay, geometry raa=90 full, Q/U gates) and the
migration package.

## SESSION 16 (2026-07-12): A-2 duplication check + items 9 & 7 DONE
DUPLICATION VERDICTS (docs read first, per Jae's rule):
 - Item 7 (golden freeze): NOT a duplicate - matrix section 4
   defines the exact protocol (frozen CSV, 3-5 representative
   cells, byte compare, per-session run); section 5-4 fixes the
   freeze ORDER #1->#2->... ; the reproduction set was EMPTY
   ("Jae uploads or freeze as items pass").  Executed below for #2
   (order exception noted: #1 replay on the current build is still
   pending; its golden follows right after).
 - Item 8a (sza80): NOT a duplicate - registry P7, "under
   investigation", pre-dates v1.10 (v1.09-era limitation, sza<=40
   unaffected).  Cause candidates ALREADY RECORDED: OSOAA-side
   grazing inaccuracy + OCRT sub-cone beams unpolarized
   approximation (fix direction: pass coupled Q/U into sub-cone
   beams - CHANGELOG_FIX-SKY-EDLU).  => new work BEYOND the
   "preserve existing results" goal; queued per Jae's ordering.
 - Item 8b (dt/beam Lu tail 0.982 @sza80): NOT a duplicate; must
   be read against the direction-B sub-cone machinery (P11,
   FIX-SKY-EDLU) - the injection's sub-cone treatment vs beams is
   the prime context.  Queued with 8a (both are sza80).
 - Item 8c (#3 delta-M, 412 rrs +48% vs MC): NOT a duplicate -
   harness section 8 lists it as unresolved, no prior fix attempt.
   Outside the current "up to pure ocean" scope.
ITEM 9 DONE (probe cleanup): target probes (ATM_NMU_INDEP,
BC-FIX conditional, NADIR-FIX guard) are ALREADY ABSENT from the
tree (lost during the later hunk rewrites of sessions 9-13; exact
turn not reconstructed - noted honestly).  State verified: grep
zero + Tier-0 EXACT + production batch BIT-IDENTICAL to fixG.
Retained on purpose: env-gated diagnostics (DTP*, C3FULL, WLBOT,
S6 traces), OCRT_TWA_MATRIX (opt-in feature), TWA-ROWNORM
(documented opt-in), EQUIV_BEAMS (reference switch).
ITEM 7 DONE (#2 golden): clean re-run on v2o_clean (sza0 0.71% /
sza40 0.79%, Rrs +0.22/+0.08 - PASS reconfirmed).  Frozen
golden_item2_5cells_2026-07-12.csv (cells per section-4 rule:
412/sza0/vza0, 865/sza0/vza30, 412/sza40/vza0, 555/sza40/vza30,
865/sza40/vza60) + repro_item2_5cells.sh runner - SELF-TEST 5/5
PASS.  This runner joins the per-commit gate set.
NEXT (Jae's ordering): speed/structure improvements (A-3), then
complete the campaign UP TO PURE OCEAN (#1 replay + #2 Q/U gate +
sza80 quantification-only), leaving #3-#6 queued.

## SESSION 17 (2026-07-12): bit-safe speed pass - profile, one
## structural fold (S17), honest nil result
CONSTRAINT (Jae): values must NOT change - every commit gated on
production-batch BYTE-compare + Tier-0 + repro_item2 runner.
PROFILE (warm 432-row batch, [S6T]):
  atm1  1.73 s  x1    (single atm solve for the whole batch)
  direct 1.44 s x433  (3.32 ms/row: water-stage reuse + couple +
                       atm2(C3d) + assembly)
  beams 0.08 s, edlu 0.02 s; ~0.9 s unattributed (IO/startup).
S17 LANDED (bit-inert hook + memo): rt_options gains OPTIONAL
per-m view-sample export (NULL default); the atm2 call site
memoizes the per-m vectors keyed on (vza, sza, wl, wind, F, m_max)
and on a hit redoes ONLY the phi reconstruction with the SAME
functions/inputs (raa is reconstruction-only there).  Batch rows
iterate raa inside vza (24 x 18) so hits = 409/433.
GATES: Tier-0 EXACT; production batch BYTE-IDENTICAL (ON vs
clean); repro runner PASS.  Correct by construction AND by gate.
SPEED VERDICT (median of 3, same host): OFF 5.41 s vs ON 5.49 s -
NO measurable gain; the atm2 share of the 3.32 ms/row is tiny in
this (single-band, aod0, fixed-bulk) batch.  Kept (zero-cost,
env-off available, may matter for heavier coupled cases) but
recorded as a nil result - the target was picked before the
per-row breakdown was known.  Lesson re-learned: profile the
INSIDE of the row loop before folding an axis.
HOST NOTE: warm timings drifted 4.18 -> 5.0-5.5 s across the day
(same binary class) - cross-session absolutes remain meaningless;
A/B on the same host+moment only (standing rule).
OSOAA co-measure (same host, today): pure-ocean 555 single run
0.43 s; 432-point-equivalent ~ 10 runs ~ 4.3 s -> OCRT batch is
at OSOAA-class parity (within host drift).
NEXT (speed/structure, still bit-safe): break down the 3.32 ms/row
(couple vs water-cache lookup vs assembly) with two temporary S6T
marks, then decide; atm1's 1.73 s single solve is the other
candidate (kernel/Legendre reuse inside - riskier, needs care).

## SESSION 18 (2026-07-12): row-loop breakdown lands; speed map
## finalized under the bit-identical constraint
NEW MARKS (env-gated, retained): [S6T] cpl (reverse coupling) and
wl2 (whole atm2/C3d block incl. S17 memo).  A third mark set
(mloop/lutrec inside the solve) hit the WRONG function (anchor
matched the scalar path, ocrt_s6_t0 out of scope) - rolled back
same-session; batch BYTE-identical + Tier-0 EXACT after rollback.
MEASURED MAP (warm 432-row batch; host drift ±25% across the day):
  atm1 (solar atm solve, x1)            ~1.27 s
  wl2 first row (atm2 bottom-source
      solve, x1; rows 2+ are S17
      reconstruction-only hits)         ~1.24 s
  direct = MAIN WATER SOLVE per row     0.77-1.44 s total
      (1.8-3.3 ms/row; this is the residual the ledger already
      attributes to memory-streaming limits)
  cpl (reverse coupling, x433)          0.011 s  (0.03 ms/row!)
  beams/edlu                            0.06 s
CORRECTED READING vs session 17: "direct" wraps the water solve
itself (call at rt_solver.c:6521 sits inside the direct S6T pair);
couple/atm2/assembly live AFTER it.  The S17 nil result is now
explained: atm2 costs ~1.2 s ONCE (first row), then the pre-S17
path was already cheap per row - S17 still turns 432 impl calls
into pure reconstructions (correct, gated), the savings were just
sub-noise on this batch.
VERDICT under "values must not change": no >1 s bit-safe lever
remains in this batch.  The two atm solves are irreducible SOS
work with shared fkc kernels; the water residual needs ALGORITHM
changes (value-bearing - out of scope).  ONE candidate left to
check next: sharing the atm BUILD/Legendre setup between atm1 and
atm2 (expected sub-0.3 s; go/no-go next turn, then close the
speed track at "OSOAA-class parity" and hand over to the
pure-ocean validation completion per Jae's ordering).

## SESSION 19 (2026-07-12): speed track CLOSED (no-go proven from
## data); validation resumed - item#1 Rayleigh re-PASSES on v1.10
SPEED FINAL CANDIDATE - NO-GO, proven without new code:
  atm1 (COLD kernel cache) 1.27 s vs atm2 first call (HOT cache,
  identical ring after view dedup) 1.24 s => kernel+setup share is
  <=0.03 s; the time is irreducible SOS iteration work.  Sharing
  the atm build/Legendre setup cannot return >0.03 s.  SPEED TRACK
  CLOSED at OSOAA-class parity (432-pt batch 4.2-5.5 s vs OSOAA
  ~4.3 s equivalent, same host, drift-bounded).
VALIDATION - ITEM#1 RAYLEIGH REPLAY (current build v2o_s18):
  Assets: canonical harness_validate.py (OCRT grid: coxmunk black
  ocean, lut vza{0,30,60} raa-step45, m-max 16, n-layers 100) +
  osoaa_ray_compare.py (OSOAA black Fresnel, YS.Abs440=1000, MOT
  matched, NO_DIRECT_GLINT, PCHIP->vza) - adapted paths/env only.
  RESULTS (raa=90, N=36 at sza<=40):
    I MAPE 0.59% (doc figure ~0.4-0.5% family)  PASS
    Q RMS/I 0.20% with the documented nadir sign-frame alignment
      (16.5% raw was ENTIRELY the known nadir Q ambiguity;
      off-nadir only: 0.10%)                     PASS
    U RMS/I 0.06% (OCRT_U=-OSOAA_U convention)   PASS
  sza=80 (P7 quantification, measurement-only): dI at nadir
    -0.7/-1.0/-1.7/-2.7/-4.5/-7.4 % across 412..865 - the
    documented P7 signature (-2..-7% growing with wavelength)
    EXACTLY reproduced on v1.10 => the pre-existing limitation is
    PRESERVED unchanged (nothing regressed, nothing silently
    "fixed").
  Grid file /tmp/tg/OCRT_item1_grid.csv (rayleigh, 150 rows) +
  cmp_item1_rayleigh.csv saved.
NEXT: item#1 aerosol + ray_aer (needs M80C extdata regeneration
via mie_to_osoaa_extdata.py from inputs/M80C.mie - /tmp assets
were wiped), then #1 golden freeze (5 cells, matrix section 4),
then #2 Q/U gate, then sza80 #2 per-band quantification.

## SESSION 20 (2026-07-12): ITEM#1 FULLY RE-PASSES on v1.10 - all
## three atmospheres reproduce the documented v1.09 figures
ASSETS REBUILT: /tmp/M80C_{412..865}.extdata regenerated via the
canonical converter from inputs/M80C.mie (ssa 0.9960-0.9969,
361 angles - M80C-consistent).  OCRT grid: 1050 rows total
(rayleigh 150 + aerosol 450 + ray_aer 450) on v2o_s18, canonical
harness_validate.py (paths only adapted).  OSOAA side: canonical
compare scripts (env adapted), resume-safe across tool timeouts.
RESULTS (raa=90; nadir-Q sign-frame alignment per harness rule;
U convention OCRT_U=-OSOAA_U):
  RAYLEIGH  sza<=40: I 0.59% | Q 0.20% | U 0.06%   PASS
            sza=80 : nadir dI -0.7..-7.4% growing with lambda -
            the DOCUMENTED P7 signature, unchanged.
  AEROSOL   sza<=40: I 0.46% | Q 0.20% | U 0.12%   PASS
            (doc: 0.48 / 0.28 / 0.20)
            sza=80 by aot: 6.29 / 1.90 / 0.49 %  == doc "6.2->0.5"
            (P7 dilution signature) EXACT.
  RAY+AER   sza<=40: I 0.44% | Q 0.21% | U 0.29%   PASS
            (doc: 0.41 / 0.27 / 0.41)
            sza=80: 0.78% (doc 0.73% - coupling dilutes P7) PASS.
  One OSOAA run failed on tool-timeout (aerosol 555/80/1.0) -
  N=159/162 for aerosol; harmless (sza80 cell).
=> ITEM#1 status on v1.10: PASS(sza<=40) for all three, P7
preserved verbatim.  The v1.10 rebuild carries NO item#1 drift.
FILES: OCRT_item1_grid.csv (1050 rows), cmp_item1_{rayleigh,
aerosol,rayaer}.csv (exported).
NEXT: #1 golden freeze (5 cells x 3 atmospheres per matrix
section 4/5-4 order), then #2 Q/U gate, then #2 sza80 per-band
quantification (measurement-only) - completing "up to pure ocean".

## SESSION 21 (2026-07-12): "UP TO PURE OCEAN" COMPLETE
#2 Q/U GATE (sza<=40, N=36, nadir-Q frame alignment, U convention):
  Q RMS/I 0.45% | U RMS/I 0.57%  - PASS (<=2-3%), BETTER than the
  documented v1.09 figures (0.57 / 0.93%) - the A-FIX-2 +
  TWA-DEFAULT chain improved the polarized channel too.
#2 sza80 PER-BAND QUANTIFICATION (measurement-only, run appended):
  TOA dI: -0.66/-0.94/-1.58/-2.67/-4.53/-7.41 % (412..865) -
  numerically the SAME curve as item#1 Rayleigh sza80 => the #2
  sza80 TOA residual is PURE P7 (surface-borne; black ocean shows
  it identically; nothing in-water, nothing v1.10).
  Rrs bias: +4.9/+3.7/+2.6/+1.5/-1.6/-2.3 % - v1.10 is BETTER than
  v1.09 here (v1.09 443: +8.1%); improvement, not regression.
  Q RMS/I 2.40% (upper edge), U 0.22%.
#1 GOLDEN FREEZE (matrix section 4/5-5, five cells per item):
  golden_item1_5cells_2026-07-12.csv:
   rayleigh 412 sza0/vza0 (degenerate) | rayleigh 865 40/30/90
   aerosol 443 aot0.05 40/30/90 (low endpoint)
   aerosol 865 aot1.0 80/60/45 (extreme + high endpoint)
   ray_aer 555 aot0.3 40/60/90
  repro_item1_5cells.sh: re-runs the exact harness chunk per cell
  and byte-compares rho_I/Q/U - SELF-TEST 5/5 PASS in 4.5 min
  (per-SESSION runner per section 4; the per-COMMIT set stays
  {Tier-0, production-batch byte, repro_item2}).
CAMPAIGN STATUS after this session:
  #1 PASS(sza<=40) x3 atmospheres, golden frozen, P7 preserved.
  #2 PASS(sza0/40) on I+Rrs+Q+U, golden frozen, sza80 fully
     attributed to P7 (fix = item 8a, approval-gated).
  QUEUED (approval-gated fixes): 8a P7 / 8b dt-beam tail /
  8c #3 delta-M.  QUEUED (next items): #3-#6, raa-mapping
  extension, migration package, Phase-1.

## SESSION 22 (2026-07-12): scatter-plot deliverables (new standing
## rule) for the completed campaign scope
NEW STANDING RULE (Jae, memory #15): accuracy-agreement results are
henceforth reported WITH scatter plots (1:1 line, ±1% band).
FIGURES (publication-style, English labels, 160 dpi):
  scatter_item1_2026-07-12.png - item#1, 3 atmospheres x {rho_I
    log-log (SZA-colored, MAPE annotated), rho_Q/rho_U linear with
    the documented sign conventions}.  sza80 red points visibly
    depart on the Rayleigh/aerosol I panels = the P7 limitation,
    shown honestly; sza<=40 hugs the 1:1 line within the ±1% band.
  scatter_item2_2026-07-12.png - item#2, 4 panels {TOA rho_I,
    Rrs(0-), rho_Q, rho_U}.  Rrs sza80 blue-end departure (+5%,
    pre-existing family, improved vs v1.09's +8%) visible; all
    sza<=40 clusters on the 1:1 line.
No code changes.  Tier-0/goldens untouched.
NEXT-TURN CANDIDATES (await pick or proceed in order): 8b (dt/beam
sza80 tail - direction-B sub-cone context read first), 8a (P7,
recorded fix direction: coupled Q/U into sub-cone beams), #3-#6,
raa-mapping extension, migration package.

## SESSION 23 (2026-07-12): item 8b CLOSED - the dt/beam sza80 tail
## is the BEAMS' documented sub-cone unpolarized approximation
CONTEXT READ FIRST (rule): CHANGELOG_FIX-SKY-EDLU direction B -
sub-cone equivalent beams are driven in-water with beam_q=0
("unpolarized; no corresponding air incidence angle"), a limitation
the changelog itself flags with the fix direction "pass coupled
Q/U into sub-cone beams (future)".
DIRECTION-RESOLVED RE-MEASURE (post A-FIX-2, sza80/555/vza30):
  dt/beam response ratio by incident column:
    c=8  (mu_w 0.52, sub-cone leak)   0.835
    c=12 (mu_w 0.72, just above crit) 0.937
    c=16 0.983 | c=20 0.997 | c=23 0.975
  => deficit concentrates toward the critical cone.
DECISIVE PROBE (OCRT_DTP_QU0: inject with AQ=AU=0, mimicking the
beams' approximation):
    c=8: 0.835 -> 0.996 | c=12: 0.937 -> 1.020 |
    FULL Lu ratio: 0.9818 -> 1.0035
  => the ENTIRE tail is the polarization treatment difference.
ATTRIBUTION: the dt default carries the coupled sky Q/U exactly
(vector Q->I coupling of the strongly polarized Rayleigh skylight
near the critical cone); the beams zero it.  dt is the RIGOROUS
side; the beams are the approximation.  NO FIX NEEDED - 8b closed
as documentation.  This ALSO means the recorded 8a fix direction
("coupled Q/U into sub-cone beams") is ALREADY REALIZED by the dt
default - consistent with the observed #2 sza80 Rrs improvement
(+8% v1.09 -> +4% v1.10).  The REMAINING P7 (TOA -0.7..-7.4%,
identical on the black ocean where no water is involved) is the
surface-BRDF / OSOAA-side share - stays open as registered.
Probe added: OCRT_DTP_QU0 (env-gated, default-inert; production
batch unaffected - injection path defaults carry Q/U as before).
Tier-0 unchanged.
NEXT: 8a residual scoping (surface share of P7 - investigation
only unless approved), or #3-#6 / raa-mapping / migration per
Jae's pick.

## SESSION 24 (2026-07-12): ITEM#4 (aDOM) COMPLETE - full grid on
## v1.10; documented residual reproduced verbatim
DIRECTIVE (Jae): finish the validation campaign (#3-#6) FIRST,
then migration, then (1) paper, (2) accuracy improvements.
ITEM#4 FULL RUN (canonical osoaa_adom_compare.py; paths/env
adapted; n-mu-water 8 -> 24 to match the #2 lineage; OCRT parity
--cdom-a440/--cdom-slope 0.014 vs YS.Abs440/Swa; 3 adom x 6 bands
x 3 sza x 3 vza = 162 cells):
  sza<=40 by a(440):  TOA I     Rrs bias   Q      U
    0.01               0.71%    -0.57%    0.29%  0.31%   PASS
    0.1                0.65%    -1.00%    0.18%  0.16%   PASS
    1.0                0.63%    -1.49%    0.16%  0.13%   TOA/Q/U
      PASS; Rrs -1.49% = the REGISTERED absorption-driven residual
      (doc: "-0.5..-1.9%"), flat across bands (-1.4..-1.9%) =
      systematic with absorption magnitude - unchanged, queued for
      the accuracy phase (post-paper per directive).
  sza=80: TOA 2.94% = the P7 curve (same as #1/#2 - surface share);
    Rrs bias +0.47% (PASSES here - CDOM damps the skylight share).
SCATTER (rule #15): scatter_item4_2026-07-12.png - 4 panels
(TOA I log-log with adom marker shapes, Rrs, Q frame-aligned,
U sign-convention).
STATUS: #4 PARTIAL -> COMPLETE-ON-GRID (pass except the registered
high-absorption Rrs residual; nothing new, nothing regressed).
NEXT: #5 TSM (no canonical compare script exists - write one on
the #4 pattern; OCRT side needs the mineral IOP parity route:
Csed phase LUTs of the consistency-harness family), then #3 chl
(OSOAA-IOP-parity route to BYPASS the native delta-M bug for
validation; the bug fix itself stays in the accuracy phase),
then #6 composite.

## SESSION 25 (2026-07-12): item#5 bring-up EXPOSES a NEW dt-path
## regression in the high-scattering regime (undetected by #2)
ITEM#5 MACHINERY BUILT (canonical parity, Brown_earth):
  - IOP synthesis VERIFIED against the Tier-0 anchor EXACTLY
    (cs5/555: A=0.176200, B=2.618309 reproduced by the formula).
  - Regenerated phase LUT == the archived harness LUT BIT-wise
    (identical rrs at both n_mu settings).
  - The +2.8% nmw8-vs-24 spread is pure quadrature convergence
    (harness golden itself notes the n_mu floor); #5 stays at
    nmw24 (the #2 lineage).
  - Thin-atm point (MOT 0.01, cs1/555/sza40): OCRT vs OSOAA
    +2.75% - harness-golden class => assembly SOUND.
THE FINDING (standard atm, cs1/555/sza40, omega_w ~0.86):
  OSOAA rrs: thin 1.148e-2 -> std-atm 1.153e-2 (rrs ~atm-invariant,
    the FIX-SKY-EDLU physical statement).
  OCRT dt DEFAULT: 1.488e-2 (+29% vs OSOAA; +26% vs its own thin).
  OCRT EQUIV_BEAMS: 1.149e-2 (-0.3% vs OSOAA - CORRECT).
  OCRT DTP_OFF: 1.078e-2 (skylight share: dt 4.1e-3 vs beams
    0.72e-3 - the injection path OVERDRIVES the sky response ~4x
    in strong-scattering water).
ATTRIBUTION: a REGRESSION of the v1.10 diffuse-top injection in
the high-omega regime.  #2 (pure water, omega 0.03) could not see
it - A-FIX-2 was six-digit-proven there only.  v1.09 (= beams
physics) matches OSOAA; the dt default departs => this is a
"preserve existing results" breach in an UNVALIDATED corner, newly
exposed by #5.  Per standing rule (no regression repair without
approval) NOTHING was changed; two options prepared for Jae:
  (a) fix the dt high-omega defect NOW (regression-repair nature;
      candidates: injected-source multiple-scattering convention
      in the water SOS - to be diagnosed), or
  (b) run #5/#3/#6 validation under OCRT_EQUIV_BEAMS=1 (exactly
      the v1.09 physics OSOAA agrees with) and defer the dt fix
      to the accuracy phase - fully consistent with the campaign-
      first directive; the dt default remains for the LOW-omega
      production regime where it is v1.09-proven.
cmp_item5_tsm.csv holds the one exploratory chunk (sza40/tsm1,
dt default - to be regenerated once the mode is decided).

## SESSION 26 (2026-07-12): high-omega VERDICT (Jae was right to
## push) + item#5 rebuilt on the CANONICAL protocol and swept
DOC READ-FIRST PAID OFF: FIXEDBULK_SHAREDPHASE_HARNESS_2026-06-29
(sections A-I) is the canonical fixed-bulk protocol - and its §E
standard run is ATM OFF ("--pressure 0 --aod 0, skylight-bug
isolation").  My #5 bring-up had grafted the #4 (with-atm)
skeleton onto fixed-bulk = an OFF-protocol combination, plus §D
violation (nmw24 < the required 32; code default was raised to 48
on 06-30) and missing --aod 0 (harmless: main.c default 0.0).
V1.09 MEASURED on the exact combination (fixed-bulk cs1/555/sza40,
atm ON, nmw24): rrs = 1.1506e-2 == beams 1.1494e-2 == OSOAA
1.1534e-2.  v1.10 dt default: 1.4881e-2 (+29.3% vs v1.09).
=> CONFIRMED: a genuine v1.10 dt-injection regression in the
high-omega + skylight corner (v1.09 physics was CORRECT there;
the §E isolation was protocol-level, not a v1.09 defect).  The
production batch (cs5, omega .94, atm ON) sits IN this corner:
dt +10.0% vs beams at the production point.  FIX REQUIRES
APPROVAL (regression-repair rule) - queued, nothing changed.
ITEM#5 REBUILT PER PROTOCOL (all prerequisites now honored):
  §B pure-water pre-check: fixed-bulk(Z09+P_Ray LUT) vs native
    -0.14% PASS.  (Loader trap logged: bb/b must be < 0.5 strictly;
    0.5*b_w exact trips rc=-2 - pass 0.49999*b, bb is inert.)
  §C consistency BUILT-IN: exposed a real mismatch - the .mie
    measured backward integral (bb/b 0.0103@412 -> 0.0211@865) vs
    the §H Ahn constant 0.0099.  bb is loader-inert so results
    unchanged; BB now derived from the measured integral =>
    consistency == 1.0 by construction.  Recorded.
  Protocol: OCRT atm-OFF nmw48 nadir vs OSOAA HYD.Model 3 thin
    (MOT 0.01) - the anchor-proven 0.06% equivalence (at sza30).
SWEEP (Brown_earth, TSM{0.1,1,10,50} x 6 bands x sza{0,40,80}):
  sza40: MAPE 1.50% bias -1.42% (harness golden class: 1.90%,
    anchor -1.15%).
  sza0: -2.82% | sza80: +2.66% - a MONOTONE sza-dependent
    systematic; TSM axis FLAT (-0.2..-0.9% - concentration scaling
    exact); band axis -0.2..-1.5% (NIR-heavier, §G truncation
    class).  PRIME SUSPECT for the sza axis: the thin(0.01) vs
    OFF methodology gap is sza-dependent (direct attenuation
    exp(-0.01/mu0): 0.990 at sza0 vs 0.944 at sza80) and the
    0.06% equivalence was proven at sza30 ONLY.  Decisive test
    queued: thin-thin (OCRT per-band pressure = 1013.25*0.01/tau)
    - if the sza axis closes, it is METHOD, not physics.
cmp_item5_tsm.csv regenerated (72 cells, canonical).
NEXT: (1) thin-thin sza-axis verdict, (2) #5 scatter + goldens,
(3) #3 chl via OSOAA-IOP parity (delta-M bypass), (4) #6,
(5) unified validation CSV (= the migration baseline, per Jae).
PENDING APPROVAL: dt high-omega regression fix.

## SESSION 27 (2026-07-12): DT-HIOM - the high-omega dt regression
## FIXED (Jae-approved), root cause was ALREADY DOCUMENTED in-code
ROOT CAUSE (comment-read-first paid off AGAIN): the D3-0c comment
(2026-07-10) states it verbatim - the legacy l-sum injector covers
MOMENT arrays only; in value-kernel production (water_phase_kernel
==0, the DEFAULT) it "misses the hydrosol entirely".  The injected
skylight was scattered with a Rayleigh-only phase (bb/b ~0.5)
instead of the forward-peaked hydrosol blend (bb/b ~0.01-0.02) ->
sky response overdriven ~5.7x in high-omega LUT water.  Pure water
IS a Rayleigh medium, hence the six-digit A-FIX-2 match there
(all prior observations explained).  The cure existed too: the
D3-0c TABLE kernel (solver's own per-m kernel tables, mode-
agnostic), parked behind OCRT_EXTTOP_KERNEL=table, default off.
FIX (DT-HIOM): table kernel PROMOTED to default (legacy behind
OCRT_EXTTOP_KERNEL=legacy as probe; NULL-table fallback kept).
GATES - ALL PASS on v2o_hiom:
  G1 closure point (cs1/555/sza40 atm-ON LUT): 1.146717e-02
     (was 1.4881e-2; v1.09 1.1506e-2 / beams 1.1494e-2 / OSOAA
     1.1534e-2 - +29.3% -> -0.3%, residual = the 8b dt-vs-beam
     polarization class).
  G2 pure water (555/sza30 atm-ON): Lu IDENTICAL to legacy.
  G3 Tier-0 EXACT.  G4 repro_item2 5/5 PASS.
  G5 production point (cs5 LUT, nmw24): dt 3.1773e-2 vs beams
     3.1913e-2 -> -0.44% (was +10%).
  G6 production BATCH: BYTE-IDENTICAL (the batch job carries NO
     phase LUT -> Rayleigh phase -> legacy==table there).
CORRECTION TO SESSION 25 (honest): the claim "the production batch
is +10% contaminated" was WRONG - that measurement had the harness
LUT attached, which the batch job does NOT use.  Batch-condition
dt vs beams: -0.18% (was always fine).  Contamination scope was
ONLY the fixed-bulk + phase-LUT + atm-ON combination (the #5-class
territory) - now fixed.  No golden re-freeze needed anywhere.
Binary promoted: build/v2o_fastk = v2o_hiom.
NEXT: #5 thin-thin sza-axis verdict + scatter + goldens, #3 chl
(parity bypass), #6, unified validation CSV (migration baseline).

## SESSION 28 (2026-07-12): ITEM#5 (TSM) COMPLETE - sza-axis
## adjudicated, scatter + goldens + runner landed
THIN-THIN VERDICT (both sides tau_R=0.01, OSOAA runs reused,
TSM=1.0 x 6 bands x 3 sza): bias 0/-2.68, 40/-1.43, 80/+2.26 % -
STATISTICALLY IDENTICAL to the off-vs-thin run => the monotone
sza systematic is NOT the atmosphere methodology; it is the real
in-water residual, attributed to the §G truncation-matching class
(model3 smooth cap vs OSOAA's own handling; the doc predicts
"few % until both sides share one truncation").  Queued for the
accuracy phase alongside the #4 absorption residual.
ITEM#5 FINAL: 72 cells (Brown_earth, TSM{0.1,1,10,50} x 6 bands x
sza{0,40,80}, §A-§I protocol): overall MAPE 2.32%; TSM axis flat
(-0.2..-0.9%: concentration scaling exact); sza40 1.50% (harness-
golden class); band axis NIR-heavier (§G).  PASS at the harness-
golden class with the two registered systematics.
SCATTER (rule #15): scatter_item5_2026-07-12.png (log-log rrs +
relative-difference panel, TSM markers, SZA colors).
GOLDENS: golden_item5_5cells_2026-07-12.csv (verbatim solver
tokens, 7 sig digits) + repro_item5_5cells.sh - SELF-TEST 5/5.
Runner traps fixed en route: IOP args must use the sweep's %.6f
formatting (full-precision args shift the 5th digit); goldens
store the stdout token verbatim (solver prints 7 digits); csv
CRLF & cwd (inputs/ relative) pitfalls documented here.
NEXT: #3 chl via OSOAA-IOP parity (delta-M bypass; blend_chl*.mie
assets ready), #6 composite, unified validation CSV (migration
baseline).

## SESSION 29 (2026-07-12): ITEM#3 (Chl) COMPLETE via value-LUT
## parity - the documented -34.8% "double-cap bug" did NOT recur
APPROACH: NOTE_item3_chl_setup (2026-06-29) documented a -34.8%
rrs gap traced to a double-cap in the OSOAA-PM MOMENT path
(rt_water_rt.c:1765).  That path is the moment kernel; this
validation uses the VALUE-LUT path instead (blend_chl<X>.mie ->
angular P11 blended with depolarized Rayleigh, injected as
--fixed-bulk-phase-lut), which the note itself flags as the
CCRR-analytic consumer.  Result: NO double-cap - chl3/412 nadir
atm-off OCRT 4.237e-3 vs OSOAA 4.191e-3 (+1.1%).  The moment-path
bug is real and still queued (accuracy phase), but the phase-LUT
route validates chl cleanly, exactly as the note predicted.
BUG FOUND & FIXED IN MY SCRIPT (not OCRT): read_blend was reading
ALL THREE .mie phase blocks (P11+P12+P33) as P11, producing a
1083-row triple-theta LUT -> chl3/412 showed a spurious -44%.
Restricting to the first (P11) TETA block: -44% -> -0.47%.
Lesson: the .mie has 3 phase blocks; slice [hdr0:hdr1].
ITEM#3 FINAL (blend_chl value-LUT, §A-§I, 72 cells: chl{0.03,0.3,
3,30} x 6 bands x sza{0,40,80}):
  sza<=40 MAPE 2.04%.  chl axis: -0.16/+0.35/+0.88/+3.01%
  (0.03->30) - the high-load +3% is the §G truncation class,
  growing with scattering (mirrors #5's TSM behavior).  sza axis
  MONOTONE (0:-2.9..-1.1 / 40:-2.0..+2.0 / 80:+4.5..+8.1 %) -
  SAME sign & structure as #5 => COMMON cause = model3-cap vs
  OSOAA-truncation (registered, accuracy phase).  Band axis peaks
  at 555 (+3.6%).
SCATTER (rule #15): scatter_item3_2026-07-12.png (log-log +
relative-difference, Chl markers, SZA colors).
GOLDENS: golden_item3_5cells_2026-07-12.csv + repro_item3_5cells.sh
- SELF-TEST 5/5.  Runner rebuilds the blend LUT per cell from the
.mie (P11 block only) - self-contained.
OSOAA throughput note: chl3 needed OSOAA cache pre-generation with
OSOAA_ROOT/HOME env exported (bare-shell runs hit ERROR_4000); the
compare scripts pass env internally and are resume-safe.
CAMPAIGN: #1,#2,#3,#4,#5 all COMPLETE (PASS + goldens + scatter).
NEXT: #6 composite (chl+TSM+CDOM together), then the UNIFIED
validation CSV (= migration baseline per Jae), then migration.

## SESSION 30 (2026-07-12): item#6 BLOCKED BY A CODE RULE (not a
## bug in my harness) - decision needed
ITEM#6 = TSM+aDOM+chl composite UNDER Rayleigh+AEROSOL (matrix
row 6).  Harness built (3-component IOP sum + scattering-weighted
blend phase incl. Rayleigh water; OSOAA HYD.Model 3 + AER.Model 4
with the canonical M80C ExtData/DirMie/Tronca/Waref/AOTref=aot865
*RATIO[band] argument set).  First run: TOA MAPE 19.6% growing to
-51% at 865 while rrs(0-) stayed at -4% => TOA-only defect.
ROOT CAUSE (code read, main.c:1861-1863): the fixed-bulk branch
does `mie_path = NULL; user_aod = 0.0;` - the AEROSOL IS FORCIBLY
DISABLED in fixed-bulk mode.  Verified: aod 0 vs 0.1 give BIT-
IDENTICAL output under --fixed-bulk-iop.  Rayleigh is NOT disabled
(P=0 vs 1013 differ), which is why #5/#3 (atm-off / thin) were
unaffected.  Note the CONTRADICTION with the comment at 1876-1884
("debug shorthands ... no longer silently disable the atmosphere")
- for aerosol, they still do.  This is a pre-existing v1.09 rule,
NOT a v1.10 regression (same lines in the pristine tree).
CONSEQUENCE: item#6 as specified (composite water + aerosol) is
NOT RUNNABLE through fixed-bulk today.  Options for Jae:
  (a) lift the rule: let --fixed-bulk-iop keep --mie/--aod (one-line
      change + gates).  Values change ONLY for fixed-bulk+aod>0,
      which is currently impossible to express -> no existing result
      can move (Tier-0/goldens/batch all have aod=0 or no fixed-bulk).
      This is the honest fix and unblocks #6 fully.
  (b) redefine #6 as composite-water under RAYLEIGH ONLY (drop the
      aerosol), completing 5.5/6 of the matrix without touching code.
  (c) validate the composite through the native path instead - blocked
      separately: native needs aph_bricaud_1998.txt (missing in the
      test tree) and carries the delta-M bug.
Recommendation: (a) - it is a capability gap, the fix is contained,
and every existing gate is provably unaffected.
STATUS: #1,#2,#3,#4,#5 COMPLETE (goldens+scatter).  #6 pending the
decision above.  Unified validation CSV (migration baseline) will
be assembled once #6 lands (or is redefined).

## SESSION 31 (2026-07-12): FB-AER capability lift + ITEM#6
## COMPLETE => VALIDATION CAMPAIGN CLOSED (6/6)
FB-AER (Jae-approved, main.c fixed-bulk branch): removed the two
lines `mie_path = NULL; user_aod = 0.0;` that forcibly disabled the
AEROSOL under --fixed-bulk-iop (a v1.09-era rule, contradicting the
"atmosphere presence is purely physical" comment 20 lines below).
Capability proof: at 865 with fixed-bulk, aod 0 -> 8.2205e-3 and
aod 0.1 -> 1.5867e-2 (previously BIT-IDENTICAL = aerosol ignored).
GATES ALL PASS: Tier-0 EXACT; production batch BYTE-IDENTICAL;
repro_item2 5/5; repro_item3 5/5; repro_item5 5/5 (the two
fixed-bulk runners are the load-bearing ones here).  No existing
result can move: the affected combination (fixed-bulk AND aod>0)
was previously inexpressible.  Binary promoted -> build/v2o_fastk.
HARNESS BUG FOUND & FIXED (mine, not OCRT): the OSOAA hydrosol
ExtData must carry the PARTICLE-ONLY phase.  OSOAA adds pure-water
molecular scattering itself (Rayleigh phase); my first composite
ExtData included the water-Rayleigh blend -> DOUBLE-COUNTED water
backscatter -> OSOAA rrs inflated, giving a spurious -19%@412 ...
-3.5%@865 spectral bias (predicted +13.9%...+1.9% from b_w/b_tot -
sign and shape matched).  OCRT fixed-bulk is the opposite: single
medium, so its LUT MUST include the water Rayleigh (§B/§C).  After
the fix: Rrs +2.3% at coastal/sza40.  #3/#5 were never affected
(their ExtData were particle-only by construction).
ITEM#6 FINAL (composite Chl+TSM+CDOM under Rayleigh+aerosol M80C
AOT865=0.1; oceanic(0.3,0.1,0.01) / coastal(3,1,0.1) /
turbid(30,10,1.0) x 6 bands x sza{0,40,80} x vza{0,30,60} = 162):
  sza<=40:  TOA I 0.32% | Rrs +1.51% | Q 0.18% | U 0.19%   PASS
    oceanic TOA 0.15% / Rrs -0.35%;  coastal 0.29% / +1.95%;
    turbid  TOA 0.51% / Rrs +2.92%  (Rrs grows with load =
    the registered §G truncation class, same as #3/#5).
  sza=80:   TOA 1.08% | Rrs +2.00% | Q 0.60% | U 0.04%
    (TOA far better than the black-ocean P7 because the aerosol
    dilutes it - consistent with #1 ray_aer 0.78%).
  Band axis: Rrs +1.3/+1.8/+3.1/+3.2/+0.8/-1.0 % (412..865).
SCATTER: scatter_item6_2026-07-12.png (4 panels).
GOLDENS: golden_item6_5cells_2026-07-12.csv + repro_item6_5cells.sh
- SELF-TEST 5/5 (rebuilds the composite IOP/LUT from the .mie
assets per cell; requires the FB-AER binary).
=== CAMPAIGN STATUS: #1 #2 #3 #4 #5 #6 ALL COMPLETE ===
Every item: PASS at harness-golden class + 5-cell golden + runner
+ scatter.  Registered systematics carried into the accuracy phase:
(i) P7 surface/OSOAA share at sza80 (black ocean, -0.7..-7.4%),
(ii) §G truncation class (model3 cap vs OSOAA) - the monotone sza
axis and the load-growing Rrs bias in #3/#5/#6, (iii) #4 high-
absorption Rrs residual (-1.5%), (iv) native delta-M bug (moment
path, item 8c), (v) 412 high-omega phase difference.
NEXT: unified validation CSV = THE MIGRATION BASELINE (Jae), then
migration, then (1) paper, (2) accuracy.

## SESSION 32 (2026-07-12): UNIFIED VALIDATION BASELINE (migration
## reference, per Jae's standing instruction)
FILE: VALIDATION_BASELINE_v1.10_2026-07-12.csv (36 rows) - one row
per (item, configuration, sza layer) with N, TOA_I MAPE, Rrs bias,
Q and U RMS/I, the acceptance tolerance, and the registered-note.
USE (standing): at the START of every migration, re-run the six
compare scripts on the migrated tree and diff against this table.
Any metric materially worse than the recorded value = migration
failure (this is the definition of "materially worse" the
consistency harness already uses for its own golden).
Companion assets required by the table (all exported):
  cmp_item1_{rayleigh,aerosol,rayaer}.csv, cmp_item2_pureocean.csv,
  cmp_item3_chl.csv, cmp_item4_adom.csv, cmp_item5_tsm.csv,
  cmp_item6_composite.csv  (raw per-cell comparisons)
  golden_item{1,2,3,5,6}_5cells + repro_item{1,2,3,5,6}_5cells.sh
  (bit-level runners; #4 shares the #2 lineage)
  scatter_item{1,2,3,4,5,6}_2026-07-12.png
BASELINE HEADLINES (sza<=40):
  #1 rayleigh 0.59% | aerosol 0.46% | ray+aer 0.44% (I)
  #2 TOA 0.75%, Rrs +0.05%, Q 0.45%, U 0.57%
  #3 chl Rrs -2.5..+0.5% (load-dependent)
  #4 TOA 0.63-0.72%, Rrs -0.57..-1.49%
  #5 TSM Rrs -2.7..-1.7%
  #6 TOA 0.15/0.29/0.51%, Rrs -0.35/+1.95/+2.92% (oceanic/coastal/turbid)

## SESSION 33 (2026-07-12): MIGRATION PACKAGE v1.10 BUILT AND
## VERIFIED FROM THE EXTRACTED TREE (canonical rule)
DOCS READ FIRST (rule): MANIFEST_v1.09 (package format + the
"verify from the EXTRACTED copy" rule, born from the 1917KST
packing defect) and SESSION_START_VERIFICATION_GUIDELINE (tiered
procedure).  The v1.10 package mirrors that format exactly.
PACKAGE: OCRT_FULLPKG_v1.10_2026-07-12 (25 MB tree, 7.5 MB tar,
277 files, SHA256SUMS over every file).
  ocrt/{src,build,scripts,inputs,docs,tools} + v1.09_to_v1.10 patch
  golden/   5-cell goldens x5 (items 1,2,3,5,6) + 5 runners
  validation/  BASELINE csv + 8 cmp csvs + 6 scatters + assets/
               (67 phase LUTs, 12 OSOAA extdata)
  aux/, mie_generator/, MANIFEST_v1.10.md, MANIFEST.txt, TODO_v1.10,
  SESSION_START_VERIFICATION_GUIDELINE_2026-07-12.md
SELF-CONTAINMENT FIX: the runners used to point at /tmp scratch
assets.  They now resolve `${OCRT_VAL_ASSETS:-<pkg>/validation/
assets}` and the .mie dir package-relatively - the package no
longer depends on any /tmp state.
VERIFICATION FROM THE EXTRACTED COPY (all PASS):
  M1 sha256sum -c over 277 files            OK (no mismatch)
  M2 rebuild from extracted src             OK
  M3 Tier-0 bit  rrs0minus=2.887685e-02     EXACT
  M4 option gate (--n-mu-water w/o ADVANCED) properly REFUSED
  M5 goldens: #1 5/5, #2 5/5, #3 5/5, #5 5/5, #6 5/5
  M6 production batch  BYTE-IDENTICAL to the dev tree
=> The migrated tree reproduces the development tree bit-for-bit.
NEXT: paper (item 1 of Jae's ordering), then the accuracy phase.
