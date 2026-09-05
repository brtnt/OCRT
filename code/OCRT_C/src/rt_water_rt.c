/* ============================================================================
 * rt_water_rt.c — In-water vector SOS solver (Phase B.3+)
 *
 * STRATEGY
 *   REUSE the atmospheric SOS engine by populating rt_atm_t with in-water
 *   values.  Two scatterer "slots" of the atm engine are repurposed:
 *     Rayleigh slot (ydel, beta2/gamma2/alpha2): pure-seawater depolarized
 *       Rayleigh phase (delta_w = 0.039), carrying the water polarization.
 *       ydel[k] = omega_w (water scattering albedo folded into the layer
 *       mixing fraction, exactly how the atm engine encodes ssa < 1).
 *     Aerosol slot (xdel, betal/gammal/alphal/zetal_aer): the PARTICLE
 *       (hydrosol) phase — scalar moments always, vector Mueller couplings
 *       when available.  xdel[k] = omega_particle.
 *   (History: an early B.3.2 attempt encoded WATER in the aerosol slot and
 *   produced |Q|^2+|U|^2 > I^2 unphysicality; the present split is the
 *   validated arrangement.)
 *
 * PHASE-FUNCTION MACHINERY IN THIS FILE
 *   - HG proxy: P_HG and its exact moments beta_l = (2l+1) g^l.
 *   - Fournier-Forand (FF) analytic particle phase + closed-form bb/b and
 *     its inversion (used for pigment/mineral scalar phases).
 *   - fixed-bulk phase LUT: CSV-tabulated P11 with log-log interpolation.
 *   - Truncation operators: "formal delta-M" HARD angular cut (mechanism
 *     concern documented at the function; the once-recorded "+48%" defect
 *     figure was RETRACTED — no measurement record exists, S-004) and the
 *     OSOAA-style log-linear forward cap (positivity-safe, production).
 *   - Angle-space VALUE kernel (fixed_bulk_direct_phase_fourier): the
 *     phi-quadrature Fourier projection of P11; mathematically equivalent
 *     to the moment kernel at exact moments (proof at the function) and the
 *     production in-water path.  This is the model for the PLANNED aerosol
 *     value kernel (HANDOFF doc).
 * ============================================================================ */

#include <time.h>
#include "rt_windows_compat.h"
#include "rt_water_rt.h"
#include "rt_value_phase.h"
#include "rt_phase_fr631.h"
#include "rt_water_iop.h"
#include "rt_iop_ahn_mineral.h"
#include "rt_iop_organic.h"
#include "rt_spectral_contract.h"
#include "rt_quadrature.h"  /* rt_quadrature_gauss_legendre_pos */
#include "rt_angle_grid.h" /* v1.10 B-0a.2: OSOAA-parity unified angle table */
#include "rt_air_water.h"
#include "shared/surface.h"
#include "rt_types.h"
#include "rt_atm.h"
#include "rt_kernel.h"
#include "rt_solver.h"
#include "rt_sos_operator.h"
#include "rt_raa_convention.h"
#include "rt_air_water_coupling.h"
#include "rt_aerosol.h"        /* rt_aerosol_compute_vector_legendre (#0 .mie vector phase) */
#include "shared/mie_io.h"     /* read_mie_file, build_aerosol_interpolators, eval_aerosol_phase */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------------------
 * v1.09 commit #21: one-shot environment snapshot for the solve path.
 * getenv() in the m-loop / per-(vza,raa) path serializes threads on platforms
 * with a locking getenv (MinGW/Windows CRT), and `static int x=-1; if(x<0)...`
 * lazy init is a data race under OpenMP.  All solve-path reads go through this
 * snapshot, initialized ONCE before any parallel region (batch driver calls it
 * right after its setenv; single-run flows hit the lazy fallback while still
 * single-threaded).  Until inited, every access falls back to plain getenv, so
 * behavior is bit-identical either way. */
typedef struct {
    int inited;
    int   f_OCRT_DEBUG;
    const char *s_OCRT_DEBUG;
    const char *s_OCRT_MB_TRACE;   /* multibounce probe gate (#21 pattern) */
    const char *s_OCRT_MB_ROUGH_VIEW; /* diagnostic-only A/B toggle */
    int   f_OCRT_DUMP_0MINUS;
    const char *s_OCRT_DUMP_0MINUS;
    int   f_OCRT_DUMP_BEAMQ;
    const char *s_OCRT_DUMP_BEAMQ;
    int   f_OCRT_DUMP_FIXEDPHASE;
    const char *s_OCRT_DUMP_FIXEDPHASE;
    int   f_OCRT_DUMP_IOP;
    const char *s_OCRT_DUMP_IOP;
    int   f_OCRT_DUMP_NORMBB;
    const char *s_OCRT_DUMP_NORMBB;
    int   f_OCRT_DUMP_PHASE;
    const char *s_OCRT_DUMP_PHASE;
    int   f_OCRT_DUMP_SSPOL;
    const char *s_OCRT_DUMP_SSPOL;
    int   f_OCRT_DUMP_TWA;
    const char *s_OCRT_DUMP_TWA;
    int   f_OCRT_DUMP_TWA_NODESET;
    const char *s_OCRT_DUMP_TWA_NODESET;
    int   f_OCRT_DUMP_TWA_PROVENANCE;
    const char *s_OCRT_DUMP_TWA_PROVENANCE;
    int   f_OCRT_DUMP_TWA_SPLINE_NODES;
    const char *s_OCRT_DUMP_TWA_SPLINE_NODES;
    int   f_OCRT_DUMP_UPWELLING;
    const char *s_OCRT_DUMP_UPWELLING;
    int   f_OCRT_DUMP_UP_RECON;
    const char *s_OCRT_DUMP_UP_RECON;
    int   f_OCRT_DUMP_VIEW;
    const char *s_OCRT_DUMP_VIEW;
    int   f_OCRT_DUMP_WATER_BOTTOM_FLUX;
    const char *s_OCRT_DUMP_WATER_BOTTOM_FLUX;
    int   f_OCRT_DUMP_WATER_IMS;
    const char *s_OCRT_DUMP_WATER_IMS;
    int   f_OCRT_DUMP_WATER_MIE_MOMENTS;
    const char *s_OCRT_DUMP_WATER_MIE_MOMENTS;
    int   f_OCRT_DUMP_WATER_SOURCE_PREF_TRACE;
    const char *s_OCRT_DUMP_WATER_SOURCE_PREF_TRACE;
    int   f_OCRT_DUMP_WATER_SOURCE_TRACE;
    const char *s_OCRT_DUMP_WATER_SOURCE_TRACE;
    int   f_OCRT_DUMP_WATER_TWA_MODE_CONTRIB;
    const char *s_OCRT_DUMP_WATER_TWA_MODE_CONTRIB;
    int   f_OCRT_DUMP_WATER_TWA_MODE_NODE_CONTRIB;
    const char *s_OCRT_DUMP_WATER_TWA_MODE_NODE_CONTRIB;
    int   f_OCRT_FINAL_MODE_QUEUE_GRID_TABLE;
    const char *s_OCRT_FINAL_MODE_QUEUE_GRID_TABLE;
    int   f_OCRT_FINAL_MODE_QUEUE_MODE;
    const char *s_OCRT_FINAL_MODE_QUEUE_MODE;
    int   f_OCRT_FINAL_MODE_QUEUE_TABLE;
    const char *s_OCRT_FINAL_MODE_QUEUE_TABLE;
    int   f_OCRT_FIXEDBULK_MSFULLOMEGA;
    const char *s_OCRT_FIXEDBULK_MSFULLOMEGA;
    int   f_OCRT_FIXEDBULK_PHASENORM;
    const char *s_OCRT_FIXEDBULK_PHASENORM;
    int   f_OCRT_NO_PHASE_CACHE;
    const char *s_OCRT_NO_PHASE_CACHE;
    int   f_OCRT_ROUTE_NAME;
    const char *s_OCRT_ROUTE_NAME;
    int   f_OCRT_TRACE_MSPEC;
    const char *s_OCRT_TRACE_MSPEC;
    int   f_OCRT_TRACE_WAVELENGTH_NM;
    const char *s_OCRT_TRACE_WAVELENGTH_NM;
    int   f_OCRT_TWA_CRITICAL_ENDPOINT_MODE;
    const char *s_OCRT_TWA_CRITICAL_ENDPOINT_MODE;
    int   f_OCRT_TWA_NODESET_MODE;
    const char *s_OCRT_TWA_NODESET_MODE;
    int   f_OCRT_TWA_OPERATOR_MODE;
    const char *s_OCRT_TWA_OPERATOR_MODE;
    int   f_OCRT_TWA_SPLINE_NODE_MODE;
    const char *s_OCRT_TWA_SPLINE_NODE_MODE;
    int   f_OCRT_TWA_SPLINT_THRESHOLD_DEG;
    const char *s_OCRT_TWA_SPLINT_THRESHOLD_DEG;
    int   f_OCRT_VZA_IN_WATER;
    const char *s_OCRT_VZA_IN_WATER;
    int   f_OCRT_WATER_GRID_CACHE;
    const char *s_OCRT_WATER_GRID_CACHE;
    int   f_OCRT_WATER_MIE_TRUNC_FULLB;
    const char *s_OCRT_WATER_MIE_TRUNC_FULLB;
    int   f_OCRT_WATER_MIE_TRUNC_TRANSPORT_TAU;
    const char *s_OCRT_WATER_MIE_TRUNC_TRANSPORT_TAU;
    int   f_OCRT_WATER_M_NO_EARLY_EXIT;
    const char *s_OCRT_WATER_M_NO_EARLY_EXIT;
    int   f_OCRT_WATER_NADIR_EPS_DEG;
    const char *s_OCRT_WATER_NADIR_EPS_DEG;
    int   f_OCRT_WATER_NADIR_RAA_DEG;
    const char *s_OCRT_WATER_NADIR_RAA_DEG;
    int   f_OCRT_WATER_NADIR_RAA_MODE;
    const char *s_OCRT_WATER_NADIR_RAA_MODE;
    int   f_OCRT_WATER_ZERO_SLOT_MODE;
    const char *s_OCRT_WATER_ZERO_SLOT_MODE;
    int   f_OCRT_WATER_ZERO_SLOT_SOLVE;
    const char *s_OCRT_WATER_ZERO_SLOT_SOLVE;
    int   f_OCRT_VALUE_PHASE_SPLINE;
    const char *s_OCRT_VALUE_PHASE_SPLINE;
    int   f_OCRT_WATER_VALUE_KERNEL_POL;
    const char *s_OCRT_WATER_VALUE_KERNEL_POL;
    const char *s_OCRT_WATER_PARTICLE_KERNEL;
} ocrt_wrt_env_t;
static ocrt_wrt_env_t g_wrt_env;

void rt_water_env_init(void) {
    if (g_wrt_env.inited) return;
    g_wrt_env.s_OCRT_DEBUG = getenv("OCRT_DEBUG");
    g_wrt_env.s_OCRT_MB_TRACE = getenv("OCRT_MB_TRACE");
    g_wrt_env.s_OCRT_MB_ROUGH_VIEW = getenv("OCRT_MB_ROUGH_VIEW");
    g_wrt_env.f_OCRT_DEBUG = (g_wrt_env.s_OCRT_DEBUG != NULL);
    g_wrt_env.s_OCRT_DUMP_0MINUS = getenv("OCRT_DUMP_0MINUS");
    g_wrt_env.f_OCRT_DUMP_0MINUS = (g_wrt_env.s_OCRT_DUMP_0MINUS != NULL);
    g_wrt_env.s_OCRT_DUMP_BEAMQ = getenv("OCRT_DUMP_BEAMQ");
    g_wrt_env.f_OCRT_DUMP_BEAMQ = (g_wrt_env.s_OCRT_DUMP_BEAMQ != NULL);
    g_wrt_env.s_OCRT_DUMP_FIXEDPHASE = getenv("OCRT_DUMP_FIXEDPHASE");
    g_wrt_env.f_OCRT_DUMP_FIXEDPHASE = (g_wrt_env.s_OCRT_DUMP_FIXEDPHASE != NULL);
    g_wrt_env.s_OCRT_DUMP_IOP = getenv("OCRT_DUMP_IOP");
    g_wrt_env.f_OCRT_DUMP_IOP = (g_wrt_env.s_OCRT_DUMP_IOP != NULL);
    g_wrt_env.s_OCRT_DUMP_NORMBB = getenv("OCRT_DUMP_NORMBB");
    g_wrt_env.f_OCRT_DUMP_NORMBB = (g_wrt_env.s_OCRT_DUMP_NORMBB != NULL);
    g_wrt_env.s_OCRT_DUMP_PHASE = getenv("OCRT_DUMP_PHASE");
    g_wrt_env.f_OCRT_DUMP_PHASE = (g_wrt_env.s_OCRT_DUMP_PHASE != NULL);
    g_wrt_env.s_OCRT_DUMP_SSPOL = getenv("OCRT_DUMP_SSPOL");
    g_wrt_env.f_OCRT_DUMP_SSPOL = (g_wrt_env.s_OCRT_DUMP_SSPOL != NULL);
    g_wrt_env.s_OCRT_DUMP_TWA = getenv("OCRT_DUMP_TWA");
    g_wrt_env.f_OCRT_DUMP_TWA = (g_wrt_env.s_OCRT_DUMP_TWA != NULL);
    g_wrt_env.s_OCRT_DUMP_TWA_NODESET = getenv("OCRT_DUMP_TWA_NODESET");
    g_wrt_env.f_OCRT_DUMP_TWA_NODESET = (g_wrt_env.s_OCRT_DUMP_TWA_NODESET != NULL);
    g_wrt_env.s_OCRT_DUMP_TWA_PROVENANCE = getenv("OCRT_DUMP_TWA_PROVENANCE");
    g_wrt_env.f_OCRT_DUMP_TWA_PROVENANCE = (g_wrt_env.s_OCRT_DUMP_TWA_PROVENANCE != NULL);
    g_wrt_env.s_OCRT_DUMP_TWA_SPLINE_NODES = getenv("OCRT_DUMP_TWA_SPLINE_NODES");
    g_wrt_env.f_OCRT_DUMP_TWA_SPLINE_NODES = (g_wrt_env.s_OCRT_DUMP_TWA_SPLINE_NODES != NULL);
    g_wrt_env.s_OCRT_DUMP_UPWELLING = getenv("OCRT_DUMP_UPWELLING");
    g_wrt_env.f_OCRT_DUMP_UPWELLING = (g_wrt_env.s_OCRT_DUMP_UPWELLING != NULL);
    g_wrt_env.s_OCRT_DUMP_UP_RECON = getenv("OCRT_DUMP_UP_RECON");
    g_wrt_env.f_OCRT_DUMP_UP_RECON = (g_wrt_env.s_OCRT_DUMP_UP_RECON != NULL);
    g_wrt_env.s_OCRT_DUMP_VIEW = getenv("OCRT_DUMP_VIEW");
    g_wrt_env.f_OCRT_DUMP_VIEW = (g_wrt_env.s_OCRT_DUMP_VIEW != NULL);
    g_wrt_env.s_OCRT_DUMP_WATER_BOTTOM_FLUX = getenv("OCRT_DUMP_WATER_BOTTOM_FLUX");
    g_wrt_env.f_OCRT_DUMP_WATER_BOTTOM_FLUX = (g_wrt_env.s_OCRT_DUMP_WATER_BOTTOM_FLUX != NULL);
    g_wrt_env.s_OCRT_DUMP_WATER_IMS = getenv("OCRT_DUMP_WATER_IMS");
    g_wrt_env.f_OCRT_DUMP_WATER_IMS = (g_wrt_env.s_OCRT_DUMP_WATER_IMS != NULL);
    g_wrt_env.s_OCRT_DUMP_WATER_MIE_MOMENTS = getenv("OCRT_DUMP_WATER_MIE_MOMENTS");
    g_wrt_env.f_OCRT_DUMP_WATER_MIE_MOMENTS = (g_wrt_env.s_OCRT_DUMP_WATER_MIE_MOMENTS != NULL);
    g_wrt_env.s_OCRT_DUMP_WATER_SOURCE_PREF_TRACE = getenv("OCRT_DUMP_WATER_SOURCE_PREF_TRACE");
    g_wrt_env.f_OCRT_DUMP_WATER_SOURCE_PREF_TRACE = (g_wrt_env.s_OCRT_DUMP_WATER_SOURCE_PREF_TRACE != NULL);
    g_wrt_env.s_OCRT_DUMP_WATER_SOURCE_TRACE = getenv("OCRT_DUMP_WATER_SOURCE_TRACE");
    g_wrt_env.f_OCRT_DUMP_WATER_SOURCE_TRACE = (g_wrt_env.s_OCRT_DUMP_WATER_SOURCE_TRACE != NULL);
    g_wrt_env.s_OCRT_DUMP_WATER_TWA_MODE_CONTRIB = getenv("OCRT_DUMP_WATER_TWA_MODE_CONTRIB");
    g_wrt_env.f_OCRT_DUMP_WATER_TWA_MODE_CONTRIB = (g_wrt_env.s_OCRT_DUMP_WATER_TWA_MODE_CONTRIB != NULL);
    g_wrt_env.s_OCRT_DUMP_WATER_TWA_MODE_NODE_CONTRIB = getenv("OCRT_DUMP_WATER_TWA_MODE_NODE_CONTRIB");
    g_wrt_env.f_OCRT_DUMP_WATER_TWA_MODE_NODE_CONTRIB = (g_wrt_env.s_OCRT_DUMP_WATER_TWA_MODE_NODE_CONTRIB != NULL);
    g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_GRID_TABLE = getenv("OCRT_FINAL_MODE_QUEUE_GRID_TABLE");
    g_wrt_env.f_OCRT_FINAL_MODE_QUEUE_GRID_TABLE = (g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_GRID_TABLE != NULL);
    g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_MODE = getenv("OCRT_FINAL_MODE_QUEUE_MODE");
    g_wrt_env.f_OCRT_FINAL_MODE_QUEUE_MODE = (g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_MODE != NULL);
    g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_TABLE = getenv("OCRT_FINAL_MODE_QUEUE_TABLE");
    g_wrt_env.f_OCRT_FINAL_MODE_QUEUE_TABLE = (g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_TABLE != NULL);
    g_wrt_env.s_OCRT_FIXEDBULK_MSFULLOMEGA = getenv("OCRT_FIXEDBULK_MSFULLOMEGA");
    g_wrt_env.f_OCRT_FIXEDBULK_MSFULLOMEGA = (g_wrt_env.s_OCRT_FIXEDBULK_MSFULLOMEGA != NULL);
    g_wrt_env.s_OCRT_FIXEDBULK_PHASENORM = getenv("OCRT_FIXEDBULK_PHASENORM");
    g_wrt_env.f_OCRT_FIXEDBULK_PHASENORM = (g_wrt_env.s_OCRT_FIXEDBULK_PHASENORM != NULL);
    g_wrt_env.s_OCRT_NO_PHASE_CACHE = getenv("OCRT_NO_PHASE_CACHE");
    g_wrt_env.f_OCRT_NO_PHASE_CACHE = (g_wrt_env.s_OCRT_NO_PHASE_CACHE != NULL);
    g_wrt_env.s_OCRT_ROUTE_NAME = getenv("OCRT_ROUTE_NAME");
    g_wrt_env.f_OCRT_ROUTE_NAME = (g_wrt_env.s_OCRT_ROUTE_NAME != NULL);
    g_wrt_env.s_OCRT_TRACE_MSPEC = getenv("OCRT_TRACE_MSPEC");
    g_wrt_env.f_OCRT_TRACE_MSPEC = (g_wrt_env.s_OCRT_TRACE_MSPEC != NULL);
    g_wrt_env.s_OCRT_TRACE_WAVELENGTH_NM = getenv("OCRT_TRACE_WAVELENGTH_NM");
    g_wrt_env.f_OCRT_TRACE_WAVELENGTH_NM = (g_wrt_env.s_OCRT_TRACE_WAVELENGTH_NM != NULL);
    g_wrt_env.s_OCRT_TWA_CRITICAL_ENDPOINT_MODE = getenv("OCRT_TWA_CRITICAL_ENDPOINT_MODE");
    g_wrt_env.f_OCRT_TWA_CRITICAL_ENDPOINT_MODE = (g_wrt_env.s_OCRT_TWA_CRITICAL_ENDPOINT_MODE != NULL);
    g_wrt_env.s_OCRT_TWA_NODESET_MODE = getenv("OCRT_TWA_NODESET_MODE");
    g_wrt_env.f_OCRT_TWA_NODESET_MODE = (g_wrt_env.s_OCRT_TWA_NODESET_MODE != NULL);
    g_wrt_env.s_OCRT_TWA_OPERATOR_MODE = getenv("OCRT_TWA_OPERATOR_MODE");
    g_wrt_env.f_OCRT_TWA_OPERATOR_MODE = (g_wrt_env.s_OCRT_TWA_OPERATOR_MODE != NULL);
    g_wrt_env.s_OCRT_TWA_SPLINE_NODE_MODE = getenv("OCRT_TWA_SPLINE_NODE_MODE");
    g_wrt_env.f_OCRT_TWA_SPLINE_NODE_MODE = (g_wrt_env.s_OCRT_TWA_SPLINE_NODE_MODE != NULL);
    g_wrt_env.s_OCRT_TWA_SPLINT_THRESHOLD_DEG = getenv("OCRT_TWA_SPLINT_THRESHOLD_DEG");
    g_wrt_env.f_OCRT_TWA_SPLINT_THRESHOLD_DEG = (g_wrt_env.s_OCRT_TWA_SPLINT_THRESHOLD_DEG != NULL);
    g_wrt_env.s_OCRT_VZA_IN_WATER = getenv("OCRT_VZA_IN_WATER");
    g_wrt_env.f_OCRT_VZA_IN_WATER = (g_wrt_env.s_OCRT_VZA_IN_WATER != NULL);
    g_wrt_env.s_OCRT_WATER_GRID_CACHE = getenv("OCRT_WATER_GRID_CACHE");
    g_wrt_env.f_OCRT_WATER_GRID_CACHE = (g_wrt_env.s_OCRT_WATER_GRID_CACHE != NULL);
    g_wrt_env.s_OCRT_WATER_MIE_TRUNC_FULLB = getenv("OCRT_WATER_MIE_TRUNC_FULLB");
    g_wrt_env.f_OCRT_WATER_MIE_TRUNC_FULLB = (g_wrt_env.s_OCRT_WATER_MIE_TRUNC_FULLB != NULL);
    g_wrt_env.s_OCRT_WATER_MIE_TRUNC_TRANSPORT_TAU = getenv("OCRT_WATER_MIE_TRUNC_TRANSPORT_TAU");
    g_wrt_env.f_OCRT_WATER_MIE_TRUNC_TRANSPORT_TAU = (g_wrt_env.s_OCRT_WATER_MIE_TRUNC_TRANSPORT_TAU != NULL);
    g_wrt_env.s_OCRT_WATER_M_NO_EARLY_EXIT = getenv("OCRT_WATER_M_NO_EARLY_EXIT");
    g_wrt_env.f_OCRT_WATER_M_NO_EARLY_EXIT = (g_wrt_env.s_OCRT_WATER_M_NO_EARLY_EXIT != NULL);
    g_wrt_env.s_OCRT_WATER_NADIR_EPS_DEG = getenv("OCRT_WATER_NADIR_EPS_DEG");
    g_wrt_env.f_OCRT_WATER_NADIR_EPS_DEG = (g_wrt_env.s_OCRT_WATER_NADIR_EPS_DEG != NULL);
    g_wrt_env.s_OCRT_WATER_NADIR_RAA_DEG = getenv("OCRT_WATER_NADIR_RAA_DEG");
    g_wrt_env.f_OCRT_WATER_NADIR_RAA_DEG = (g_wrt_env.s_OCRT_WATER_NADIR_RAA_DEG != NULL);
    g_wrt_env.s_OCRT_WATER_NADIR_RAA_MODE = getenv("OCRT_WATER_NADIR_RAA_MODE");
    g_wrt_env.f_OCRT_WATER_NADIR_RAA_MODE = (g_wrt_env.s_OCRT_WATER_NADIR_RAA_MODE != NULL);
    g_wrt_env.s_OCRT_WATER_ZERO_SLOT_MODE = getenv("OCRT_WATER_ZERO_SLOT_MODE");
    g_wrt_env.f_OCRT_WATER_ZERO_SLOT_MODE = (g_wrt_env.s_OCRT_WATER_ZERO_SLOT_MODE != NULL);
    g_wrt_env.s_OCRT_WATER_ZERO_SLOT_SOLVE = getenv("OCRT_WATER_ZERO_SLOT_SOLVE");
    g_wrt_env.f_OCRT_WATER_ZERO_SLOT_SOLVE = (g_wrt_env.s_OCRT_WATER_ZERO_SLOT_SOLVE != NULL);
    g_wrt_env.s_OCRT_VALUE_PHASE_SPLINE = getenv("OCRT_VALUE_PHASE_SPLINE");
    g_wrt_env.f_OCRT_VALUE_PHASE_SPLINE = (g_wrt_env.s_OCRT_VALUE_PHASE_SPLINE && strcmp(g_wrt_env.s_OCRT_VALUE_PHASE_SPLINE, "0") != 0);
    g_wrt_env.s_OCRT_WATER_VALUE_KERNEL_POL = getenv("OCRT_WATER_VALUE_KERNEL_POL");
    g_wrt_env.s_OCRT_WATER_PARTICLE_KERNEL = getenv("OCRT_WATER_PARTICLE_KERNEL");
    /* Stage 3B production contract: untruncated particle Mie phases use the
     * direct theta-linear vector value kernel by default.  The legacy
     * coefficient path remains available for controlled A/B regression via
     * OCRT_WATER_PARTICLE_KERNEL=moment or OCRT_WATER_VALUE_KERNEL_POL=0. */
    if (g_wrt_env.s_OCRT_WATER_PARTICLE_KERNEL &&
        strcmp(g_wrt_env.s_OCRT_WATER_PARTICLE_KERNEL, "moment") == 0) {
        g_wrt_env.f_OCRT_WATER_VALUE_KERNEL_POL = 0;
    } else if (g_wrt_env.s_OCRT_WATER_PARTICLE_KERNEL &&
               strcmp(g_wrt_env.s_OCRT_WATER_PARTICLE_KERNEL, "direct") == 0) {
        g_wrt_env.f_OCRT_WATER_VALUE_KERNEL_POL = 1;
    } else if (g_wrt_env.s_OCRT_WATER_VALUE_KERNEL_POL) {
        g_wrt_env.f_OCRT_WATER_VALUE_KERNEL_POL =
            strcmp(g_wrt_env.s_OCRT_WATER_VALUE_KERNEL_POL, "0") != 0;
    } else {
        g_wrt_env.f_OCRT_WATER_VALUE_KERNEL_POL = 1;
    }
    g_wrt_env.inited = 1;
}

static int water_particle_direct_enabled(int truncation_mode, int ss_mode)
{
    /* The direct theta-linear path supports the untruncated phase and the
     * OSOAA-compatible real-angle forward-cap residual.  Single-scatter
     * correction modes still belong to the legacy coefficient diagnostic
     * path and must never be silently mixed with the direct representation. */
    if ((truncation_mode != 0 && truncation_mode != 1) || ss_mode != 0)
        return 0;
    return g_wrt_env.f_OCRT_WATER_VALUE_KERNEL_POL != 0;
}


/* ---------------------------------------------------------------------------
 * Debug-only environment flags (developer reference: DEBUG_FLAGS.md).
 * A debug flag that can change the physics result is honored ONLY when the
 * master gate OCRT_DEBUG is also set, so a debug toggle can never silently
 * alter a production run. These are NOT user-facing options (absent from CLI
 * --help and user docs); they exist for development / regression isolation.
 * Pure diagnostic dumps (OCRT_DUMP_*) and physics-invariant toggles are read
 * with plain getenv() since they cannot corrupt a result.
 * --------------------------------------------------------------------------- */


/* Step92 diagnostic-only final assembled m=0/m=1 queue-equivalent correction.
 * Default-off.  It never modifies the recursive water-SOS fields; it only adds
 * an externally supplied final assembled-mode correction to the final Lu0+
 * output.  Step91 table coefficients are in the pi-scaled ADV/LUM-equivalent
 * mode-contribution units, so we divide by pi before adding to Lu0+.
 *
 * Env:
 *   OCRT_FINAL_MODE_QUEUE_MODE=table
 *   OCRT_FINAL_MODE_QUEUE_TABLE=/path/to/step91_final_mode_queue_correction_table.csv
 */
typedef struct {
    double wl;
    double vza;
    int comp;  /* 0 I, 1 Q */
    int m;     /* 0 common, 1 odd */
    double coeff;
} ocrt_step92_qrow_t;

static int ocrt_step92_queue_mode(void) {
    const char *s = g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_MODE;
    if (!s || !s[0] || !strcmp(s,"off") || !strcmp(s,"0") || !strcmp(s,"baseline")) return 0;
    if (!strcmp(s,"table") || !strcmp(s,"final_table") || !strcmp(s,"queue_table")) return 1;
    return 0;
}

static int ocrt_step92_load_table(ocrt_step92_qrow_t **rows_out, int *n_out) {
    static ocrt_step92_qrow_t *rows = NULL;
    static int n = -1;
    if (n >= 0) { *rows_out = rows; *n_out = n; return 0; }
    n = 0;
    const char *path = g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_TABLE;
    if (!path || !path[0]) { *rows_out = NULL; *n_out = 0; return -1; }
    FILE *fp = fopen(path, "r");
    if (!fp) { *rows_out = NULL; *n_out = 0; return -2; }
    int cap = 512;
    rows = (ocrt_step92_qrow_t*)calloc((size_t)cap, sizeof(*rows));
    if (!rows) { fclose(fp); return -3; }
    char line[2048];
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); *rows_out = rows; *n_out = n; return 0; }
    while (fgets(line, sizeof(line), fp)) {
        double wl=0.0, vza=0.0, coeff=0.0, norm=0.0;
        int m=0;
        char compstr[16]; compstr[0]='\0';
        int ok = sscanf(line, " %lf,%lf,%15[^,],%d,%lf,%lf,", &wl, &vza, compstr, &m, &coeff, &norm);
        if (ok < 5) continue;
        int comp = -1;
        if (compstr[0] == 'I') comp = 0;
        else if (compstr[0] == 'Q') comp = 1;
        else continue;
        if (!(m == 0 || m == 1)) continue;
        if (n >= cap) {
            cap *= 2;
            ocrt_step92_qrow_t *tmp = (ocrt_step92_qrow_t*)realloc(rows, (size_t)cap*sizeof(*rows));
            if (!tmp) break;
            rows = tmp;
        }
        rows[n].wl = wl; rows[n].vza = vza; rows[n].comp = comp; rows[n].m = m; rows[n].coeff = coeff;
        n++;
    }
    fclose(fp);
    *rows_out = rows; *n_out = n; return 0;
}

static double ocrt_step92_interp_coeff(double wl_nm, double vza_deg, int comp, int m) {
    ocrt_step92_qrow_t *rows = NULL; int n = 0;
    if (ocrt_step92_load_table(&rows, &n) != 0 || !rows || n <= 0) return 0.0;
    double best_wld = 1e99;
    for (int i=0;i<n;i++) if (rows[i].comp==comp && rows[i].m==m) {
        double d = fabs(rows[i].wl - wl_nm);
        if (d < best_wld) best_wld = d;
    }
    if (!(best_wld < 5.0)) return 0.0;
    double xlo=-1e99, xhi=1e99, ylo=0.0, yhi=0.0;
    int have_lo=0, have_hi=0;
    for (int i=0;i<n;i++) {
        if (rows[i].comp!=comp || rows[i].m!=m) continue;
        if (fabs(rows[i].wl - wl_nm) > best_wld + 1e-9) continue;
        double x=rows[i].vza, y=rows[i].coeff;
        if (x <= vza_deg && x > xlo) { xlo=x; ylo=y; have_lo=1; }
        if (x >= vza_deg && x < xhi) { xhi=x; yhi=y; have_hi=1; }
    }
    if (have_lo && have_hi) {
        if (fabs(xhi-xlo) < 1e-12) return ylo;
        double t=(vza_deg-xlo)/(xhi-xlo);
        return ylo + t*(yhi-ylo);
    }
    if (have_lo) return ylo;
    if (have_hi) return yhi;
    return 0.0;
}

static void ocrt_step92_apply_final_mode_queue(double lambda_nm, double vza_deg, double raa_deg,
                                                double *I_air, double *Q_air, double *U_air) {
    (void)U_air;
    if (ocrt_step92_queue_mode() != 1) return;
    const double phi_arg = raa_deg * M_PI / 180.0 + RT_F_SOLAR_PI;
    const double c1 = cos(phi_arg);
    const double dI0 = ocrt_step92_interp_coeff(lambda_nm, vza_deg, 0, 0);
    const double dI1 = ocrt_step92_interp_coeff(lambda_nm, vza_deg, 0, 1);
    const double dQ0 = ocrt_step92_interp_coeff(lambda_nm, vza_deg, 1, 0);
    const double dQ1 = ocrt_step92_interp_coeff(lambda_nm, vza_deg, 1, 1);
    if (I_air) *I_air += (dI0 + dI1 * c1) / M_PI;
    if (Q_air) *Q_air += (dQ0 + dQ1 * c1) / M_PI;
}



/* Step114 diagnostic/provisional final grid/LUT queue table.
 * Default-off. Unlike Step92 mode-contribution tables, this table is already
 * in final Rrs units and is keyed by wavelength/SZA/VZA/RAA. It is intended
 * only for discrete validated grid/LUT products. It does not interpolate; if
 * the exact grid key is missing, no correction is applied.  This prevents an
 * accidental continuous-VZA production claim in the upper-guard region.
 *
 * Env:
 *   OCRT_FINAL_MODE_QUEUE_MODE=grid_table
 *   OCRT_FINAL_MODE_QUEUE_GRID_TABLE=/path/to/step114_runtime_grid_queue_table_Rrs_units.csv
 */
typedef struct {
    double wl, sza, vza, raa;
    double dI, dQ, dU; /* final Rrs-unit deltas */
} ocrt_step114_gridrow_t;

static int ocrt_step114_grid_mode(void) {
    const char *s = g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_MODE;
    if (!s || !s[0]) return 0;
    return (!strcmp(s,"grid_table") || !strcmp(s,"grid") || !strcmp(s,"lut_grid"));
}

static int ocrt_step114_load_grid(ocrt_step114_gridrow_t **rows_out, int *n_out) {
    static ocrt_step114_gridrow_t *rows = NULL;
    static int n = -1;
    if (n >= 0) { *rows_out = rows; *n_out = n; return 0; }
    n = 0;
    const char *path = g_wrt_env.s_OCRT_FINAL_MODE_QUEUE_GRID_TABLE;
    if (!path || !path[0]) { *rows_out = NULL; *n_out = 0; return -1; }
    FILE *fp = fopen(path, "r");
    if (!fp) { *rows_out = NULL; *n_out = 0; return -2; }
    int cap = 2048;
    rows = (ocrt_step114_gridrow_t*)calloc((size_t)cap, sizeof(*rows));
    if (!rows) { fclose(fp); return -3; }
    char line[4096];
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); *rows_out = rows; *n_out = n; return 0; }
    while (fgets(line, sizeof(line), fp)) {
        double wl=0.0, sza=0.0, vza=0.0, raa=0.0, dI=0.0, dQ=0.0, dU=0.0;
        /* Table starts: wavelength_nm,sza_deg,vza_deg,raa_deg,delta_I_Rrs,delta_Q_Rrs,delta_U_Rrs,... */
        int ok = sscanf(line, " %lf,%lf,%lf,%lf,%lf,%lf,%lf,", &wl, &sza, &vza, &raa, &dI, &dQ, &dU);
        if (ok < 7) continue;
        if (n >= cap) {
            cap *= 2;
            ocrt_step114_gridrow_t *tmp = (ocrt_step114_gridrow_t*)realloc(rows, (size_t)cap*sizeof(*rows));
            if (!tmp) break;
            rows = tmp;
        }
        rows[n].wl=wl; rows[n].sza=sza; rows[n].vza=vza; rows[n].raa=raa; rows[n].dI=dI; rows[n].dQ=dQ; rows[n].dU=dU;
        n++;
    }
    fclose(fp);
    *rows_out = rows; *n_out = n; return 0;
}

static double ocrt_step114_angle_diff_deg(double a, double b) {
    double d = fmod(a - b + 540.0, 360.0) - 180.0;
    return fabs(d);
}

static int ocrt_step114_lookup_grid_delta(double wl, double sza, double vza, double raa,
                                           double *dI, double *dQ, double *dU) {
    ocrt_step114_gridrow_t *rows = NULL; int n = 0;
    if (ocrt_step114_load_grid(&rows, &n) != 0 || !rows || n <= 0) return 0;
    const double twl=0.51, tang=0.015;
    for (int i=0;i<n;i++) {
        if (fabs(rows[i].wl - wl) > twl) continue;
        if (fabs(rows[i].sza - sza) > tang) continue;
        if (fabs(rows[i].vza - vza) > tang) continue;
        if (ocrt_step114_angle_diff_deg(rows[i].raa, raa) > tang) continue;
        if (dI) *dI = rows[i].dI;
        if (dQ) *dQ = rows[i].dQ;
        if (dU) *dU = rows[i].dU;
        return 1;
    }
    return 0;
}

static const char *ocrt_debug_env(const char *name) {
    static int dbg = -1;
    if (dbg < 0) dbg = (g_wrt_env.s_OCRT_DEBUG != NULL);
    return dbg ? getenv(name) : NULL;
}



/* Step162 production default: physical-strict water->air TWA route.
 *
 * Historical pointwise/legacy/natural/hybrid interface routes were diagnostic
 * branches that allowed a requested air-side direction to be sampled before
 * the water->air transmission operator.  That ordering is not the physical
 * air-water interface operator used for OSOAA-compatible validation: the
 * water-side radiance field must first be transmitted through T_wa at the
 * water-side nodes, mapped by Snell to air-side nodes, and then sampled at
 * the requested air-side view direction.
 *
 * Therefore the production default is now the OSOAA-compatible clamped SPLINT
 * operator.  Non-physical legacy route switches are intentionally ignored; if
 * a caller tries to request them, the solver remains on the strict path and
 * emits a one-time diagnostic warning to stderr.
 *
 * This is a route-contract change only: no empirical gain, normalization
 * tuning, residual fitting, or queue interpolation is introduced here.
 */
static int ocrt_water_twa_operator_mode(void) {
    const char *m = g_wrt_env.s_OCRT_TWA_OPERATOR_MODE;
    static int warned = 0;
    if (m && m[0] && strcmp(m, "osoaa_spline") != 0 && strcmp(m, "osoaa") != 0 &&
        strcmp(m, "spline") != 0 && strcmp(m, "osoaa_clamped") != 0 &&
        strcmp(m, "clamped") != 0 && strcmp(m, "physical_strict") != 0 &&
        strcmp(m, "strict") != 0) {
        if (!warned) {
            fprintf(stderr, "warning: deprecated OCRT_TWA_OPERATOR_MODE='%s' ignored; using physical_strict/osoaa_clamped\n", m);
            warned = 1;
        }
    }
    return 2;
}

static double ocrt_water_twa_splint_threshold_deg(void) {
    const char *e = g_wrt_env.s_OCRT_TWA_SPLINT_THRESHOLD_DEG;
    if (!e || !e[0]) return 60.0;
    char *endp = NULL;
    double v = strtod(e, &endp);
    if (endp == e || !isfinite(v)) return 60.0;
    if (v < 0.0) v = 0.0;
    if (v > 90.0) v = 90.0;
    return v;
}

static int ocrt_water_twa_operator_mode_for_vza(double vza_deg) {
    (void)vza_deg;
    return ocrt_water_twa_operator_mode();
}


/* Step 55 critical-edge endpoint/extrapolation diagnostic.
 * Applies only to the operator-first TWA path when the requested air-side
 * cosine lies outside the transmitted-node spline interval, which happens near
 * the water-air critical angle.  Default keeps Step 51/52 behaviour.
 *
 *   OCRT_TWA_CRITICAL_ENDPOINT_MODE=splint_extrapolate  default
 *   OCRT_TWA_CRITICAL_ENDPOINT_MODE=linear_endpoint     linear endpoint extrapolation
 *   OCRT_TWA_CRITICAL_ENDPOINT_MODE=clamp_first         clamp below first node to first value
 *   OCRT_TWA_CRITICAL_ENDPOINT_MODE=add_muair0_zero     add (mu_air=0, L=0) to spline grid
 */
static int ocrt_water_twa_critical_endpoint_mode(void) {
    const char *m = g_wrt_env.s_OCRT_TWA_CRITICAL_ENDPOINT_MODE;
    if (!m || !m[0] || strcmp(m, "splint_extrapolate") == 0 ||
        strcmp(m, "spline") == 0 || strcmp(m, "current") == 0) return 0;
    if (strcmp(m, "linear_endpoint") == 0 || strcmp(m, "linear") == 0) return 1;
    if (strcmp(m, "clamp_first") == 0 || strcmp(m, "clamp") == 0) return 2;
    if (strcmp(m, "add_muair0_zero") == 0 || strcmp(m, "add_zero") == 0 ||
        strcmp(m, "mu0_zero") == 0) return 3;
    return 0;
}

static int ocrt_water_twa_spline_node_mode(void) {
    const char *m = g_wrt_env.s_OCRT_TWA_SPLINE_NODE_MODE;
    static int warned = 0;
    if (m && m[0] && strcmp(m, "osoaa_slots") != 0 && strcmp(m, "osoaa") != 0 &&
        strcmp(m, "rmu19") != 0 && strcmp(m, "zero_slots") != 0 &&
        strcmp(m, "physical_strict") != 0 && strcmp(m, "strict") != 0) {
        if (!warned) {
            fprintf(stderr, "warning: deprecated OCRT_TWA_SPLINE_NODE_MODE='%s' ignored; using osoaa_slots\n", m);
            warned = 1;
        }
    }
    return 1;
}

static double ocrt_cubic_spline_eval_sorted(int n_in, const double *x_in, const double *y_in,
                                            double xp, int clamped_end_slopes);


/* Step 56 diagnostic/provisional node-set parity switch.
 * OSOAA's flat TWA operator receives the full OSOAA positive RMU list, which
 * can include zero-weight user/solar directions in addition to weighted Gauss
 * nodes.  OCRT's Step51 operator used weighted water nodes only.  This switch
 * lets us append OSOAA-like user directions to the operator-first node set for
 * critical-edge parity diagnostics.
 *
 *   OCRT_TWA_NODESET_MODE=weighted       default, Step51/52 behaviour
 *   OCRT_TWA_NODESET_MODE=osoaa_user     append mu=1, mu_sun_air, mu_sun_water
 */
static int ocrt_water_twa_nodeset_mode(void) {
    const char *m = g_wrt_env.s_OCRT_TWA_NODESET_MODE;
    static int warned = 0;
    if (m && m[0] && strcmp(m, "osoaa_user") != 0 && strcmp(m, "osoaa_nodes") != 0 &&
        strcmp(m, "osoaa_rmu") != 0 && strcmp(m, "augmented") != 0 &&
        strcmp(m, "physical_strict") != 0 && strcmp(m, "strict") != 0) {
        if (!warned) {
            fprintf(stderr, "warning: deprecated OCRT_TWA_NODESET_MODE='%s' ignored; using osoaa_user\n", m);
            warned = 1;
        }
    }
    return 1;
}


/* Step 58 zero-weight/user-direction solve diagnostic.
 * OSOAA keeps several zero-quadrature-weight RMU slots (mu=1, mu_sun_water,
 * mu_sun_air) in the discrete ordinate list.  Those slots do not contribute as
 * incoming directions to the angular quadrature, but the SOS transport/source
 * sweep still produces a real outgoing field at those directions, which then
 * enters the flat TWA operator and its SPLINT input.  Step56 showed that
 * appending the slots only at the TWA stage and filling them by interpolation
 * is not equivalent.  This guarded route inserts the slots into the water SOS
 * ordinate set itself with gb=0.
 *
 *   OCRT_WATER_ZERO_SLOT_SOLVE=0|off        default
 *   OCRT_WATER_ZERO_SLOT_SOLVE=osoaa|1|on   append mu=1, mu_sun_water, mu_sun_air
 *   OCRT_WATER_ZERO_SLOT_MODE=solve_zero_weight is accepted as an alias.
 */
static int ocrt_water_zero_slot_solve_mode(void) {
    const char *m = g_wrt_env.s_OCRT_WATER_ZERO_SLOT_SOLVE;
    if (!m || !m[0]) m = g_wrt_env.s_OCRT_WATER_ZERO_SLOT_MODE;
    static int warned = 0;
    if (m && m[0] && strcmp(m,"1") != 0 && strcmp(m,"on") != 0 && strcmp(m,"osoaa") != 0 &&
        strcmp(m,"osoaa_slots") != 0 && strcmp(m,"zero_slots") != 0 &&
        strcmp(m,"solve_zero_weight") != 0 && strcmp(m,"solve") != 0 &&
        strcmp(m,"physical_strict") != 0 && strcmp(m,"strict") != 0) {
        if (!warned) {
            fprintf(stderr, "warning: deprecated OCRT_WATER_ZERO_SLOT_SOLVE='%s' ignored; solving zero-weight OSOAA slots\n", m);
            warned = 1;
        }
    }
    return 1;
}

static void ocrt_water_sort_positive_nodes(rt_atm_t *atm) {
    if (!atm || atm->n_mu <= 1) return;
    const int n = atm->n_mu;
    for (int i = 1; i <= n; ++i) {
        for (int j = i + 1; j <= n; ++j) {
            if (atm->rm[+i] > atm->rm[+j]) {
                double tm = atm->rm[+i]; atm->rm[+i] = atm->rm[+j]; atm->rm[+j] = tm;
                double tw = atm->gb[+i]; atm->gb[+i] = atm->gb[+j]; atm->gb[+j] = tw;
            }
        }
    }
    for (int i = 1; i <= n; ++i) {
        atm->rm[-i] = -atm->rm[+i];
        atm->gb[-i] =  atm->gb[+i];
    }
}

/* v1.10 B-0c.2b: add the FIRST-SCATTER source of a per-m diffuse top
 * incident field, column-parameterized mirror of
 * rt_solver_primary_source(+_pol) (direct l-sums over ws tables so
 * signed ring columns are covered, exactly like the v1.059 bq blocks).
 * Amplitude convention: A = 2*w_c*S^m(c) (== F_eq/pi of the equivalent
 * beams; pinned by the m0 A/B gate).  U-incident is DROPPED in this
 * increment (beam machinery had no U either); Q-incident uses the
 * (1,2)/(2,2)/(3,2) mirrors. */
static void ocrt_add_diffuse_top_primary(const rt_atm_t *atm,
                                         const rt_legendre_workspace_t *ws,
                                         int m,
                                         const double *ext_I,
                                         const double *ext_Q,
                                         const double *ext_U,
                                         const double *ext_mu, int ext_n,
                                         double *src_i, double *src_q,
                                         double *src_u) {
    const int nt   = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const int l_max = ws->l_max;
    const int ray_on = (m <= 2) ? 1 : 0;
    const double beta0_m = (m == 0) ? atm->beta0 : 0.0;
    const double beta2 = atm->beta2, gamma2 = atm->gamma2;
    const double alpha2 = atm->alpha2;
    /* v1.10 D3-0b (2026-07-10): U-incident column added (the call-site
     * comment reserved it as "next increment").  Pattern mirrors the
     * verified I/Q columns; negative-index Legendre evaluation at cc
     * absorbs up/down parity.  Correctness enforced numerically by the
     * D3 reciprocity gate and the m=0 U-response==0 check. */
    double *rowI = (double*)malloc((size_t)dirs * 9 * sizeof(double));
    if (!rowI) return;
    double *rowQ = rowI + dirs, *rowU = rowQ + dirs;
    double *rowIq = rowU + dirs, *rowQq = rowIq + dirs, *rowUq = rowQq + dirs;
    double *rowIu = rowUq + dirs, *rowQu = rowIu + dirs, *rowUu = rowQu + dirs;

    for (int c = 1; c <= n_mu; ++c) {
        const double mu_c = atm->rm[+c];
        /* map ring column -> incident sample */
        double AI = 0.0, AQ = 0.0, AU = 0.0;
        for (int e = 0; e < ext_n; ++e) {
            if (fabs(ext_mu[e] - mu_c) < 1e-12) {
                AI = 2.0 * atm->gb[+c] * ext_I[e];
                AQ = ext_Q ? 2.0 * atm->gb[+c] * ext_Q[e] : 0.0;
                AU = ext_U ? 2.0 * atm->gb[+c] * ext_U[e] : 0.0;
                break;
            }
        }
        if (getenv("OCRT_D3_DIAG") && c <= 1)
            fprintf(stderr, "[injdiag] m=%d c=%d AI=%.3e AQ=%.3e AU=%.3e aer=%d b0=%.3e b2=%.3e g2=%.3e a2=%.3e grT=%.3e artT=%.3e\n",
                    m, c, AI, AQ, AU, atm->aerosol_active,
                    beta0_m, beta2, gamma2, alpha2,
                    atm->xrl ? atm->xrl[+1] : -999.0,
                    (ws && ws->gr_pol) ? ws->gr_pol[+2][-3] : -999.0);
        if (AI == 0.0 && AQ == 0.0 && AU == 0.0) {
            /* v1.10 DTP-FIX (2026-07-12): the incident diffuse field is
             * sampled on the water GL nodes only; ring columns that are
             * NOT GL nodes (view / slot / nadir inserts) found no exact
             * match above and previously received NO skylight at all --
             * the item#2 regression (nmw8 -3.1%, nmw24 -0.12%; probe
             * [DTPmiss] showed the missing columns are exactly the
             * inserted ones, e.g. mu=cos(sza), mu_view, 1.0).  The
             * equivalent-beam mechanism this injector replaced fed EVERY
             * column exactly.  Restore that property by linear
             * interpolation of the ext field in mu (clamped at the ends);
             * the diffuse field is smooth in mu, and GL columns still take
             * the exact-match branch above (bit-preserving for them). */
            if (ext_n >= 2) {
                int lo = 0;
                while (lo + 1 < ext_n - 1 && ext_mu[lo + 1] <= mu_c) ++lo;
                const double m0_ = ext_mu[lo], m1_ = ext_mu[lo + 1];
                double t_ = (m1_ != m0_) ? (mu_c - m0_) / (m1_ - m0_) : 0.0;
                if (t_ < 0.0) t_ = 0.0; else if (t_ > 1.0) t_ = 1.0;
                AI = 2.0 * atm->gb[+c] * ((1.0 - t_) * ext_I[lo] + t_ * ext_I[lo + 1]);
                if (ext_Q) AQ = 2.0 * atm->gb[+c] * ((1.0 - t_) * ext_Q[lo] + t_ * ext_Q[lo + 1]);
                if (ext_U) AU = 2.0 * atm->gb[+c] * ((1.0 - t_) * ext_U[lo] + t_ * ext_U[lo + 1]);
            }
            if (getenv("OCRT_DTP_DIAG")) {
                static int dtp_fix = 0;
                fprintf(stderr, "[DTPfix] m=%d c=%d mu_c=%.9f interp AI=%.6e (fix#%d)\n",
                        m, c, mu_c, AI, ++dtp_fix);
            }
            if (AI == 0.0 && AQ == 0.0 && AU == 0.0) continue;
        }
        { const char *ms_ = getenv("OCRT_DTP_MSCALE");   /* v1.10 probe: m>=1 amplitude scale */
          if (ms_ && m > 0) { double f_ = atof(ms_); AI *= f_; AQ *= f_; AU *= f_; } }
        { const char *q0_ = getenv("OCRT_DTP_QU0");   /* v1.10 probe: mimic the beams' sub-cone unpolarized approx (AQ=AU=0) */
          if (q0_) { AQ = 0.0; AU = 0.0; } }
        { const char *as_ = getenv("OCRT_DTP_ASCALE");   /* v1.10 probe: ALL-m amplitude scale */
          if (as_) { double f_ = atof(as_); AI *= f_; AQ *= f_; AU *= f_; } }
        const int cc = -c;                     /* downward incident column */
        /* per-out-direction kernel rows via direct l-sums */
        for (int j = -n_mu; j <= n_mu; ++j) {
            double sI = 0.0, sQ = 0.0, sU = 0.0;
            double qI = 0.0, qQ = 0.0, qU = 0.0;
            double uI = 0.0, uQ = 0.0, uU = 0.0;
            if (atm->aerosol_active) {
                for (int l = m; l <= l_max; ++l) {
                    const double pj = ws->plm[l][j], pc = ws->plm[l][cc];
                    const double rj = ws->rrl[l][j], rc = ws->rrl[l][cc];
                    const double tj = ws->rtl[l][j], tc = ws->rtl[l][cc];
                    const double b = atm->betal_aer[l], g = atm->gammal_aer[l];
                    const double aa = atm->alphal_aer[l], z = atm->zetal_aer[l];
                    sI += pj * pc * b;
                    sQ += rj * pc * g;         /* (2,1) */
                    sU += tj * pc * g;         /* (3,1) */
                    qI += pj * rc * g;         /* (1,2) */
                    qQ += tc * tj * z + rc * rj * aa;   /* (2,2) */
                    qU += tj * rc * aa + rj * tc * z;   /* (3,2) */
                    uI += pj * tc * g;                  /* (1,3) D3-0b */
                    uQ += rj * tc * aa + tj * rc * z;   /* (2,3) D3-0b */
                    uU += tj * tc * aa + rj * rc * z;   /* (3,3) D3-0b */
                }
            }
            const double xplc = atm->xpl[cc], xrlc = atm->xrl ? atm->xrl[cc] : 0.0;
            const double raI = beta0_m + beta2 * atm->xpl[j] * xplc;
            const double raQ = ray_on ? gamma2 * atm->xrl[j] * xplc : 0.0;
            const double raU = ray_on ? gamma2 * atm->xtl[j] * xplc : 0.0;
            const double rqI = ray_on ? gamma2 * atm->xpl[j] * xrlc : 0.0;
            const double rqQ = ray_on ? alpha2 * atm->xrl[j] * xrlc : 0.0;
            const double rqU = ray_on ? alpha2 * atm->xtl[j] * xrlc : 0.0;
            /* v1.10 D3-0c (2026-07-10): TABLE kernel mode (env
             * OCRT_EXTTOP_KERNEL=table, live getenv - runtime-set by the
             * D3 driver; #22 rule respected).  The legacy l-sums above
             * cover MOMENT arrays only; in value-kernel production
             * (wpk=0) they miss the hydrosol entirely (D3-0b finding).
             * The solver's own per-m kernel tables (ws->phase_fourier_m,
             * gr/gt/arr/art/att_pol) are mode-agnostic; the nine entries
             * below are a MECHANICAL transcription of the DOWN-INCIDENT
             * terms of rt_sos_operator_apply_vector (in=c at -c, out=j signed;
             * *_dn lines), with the U-output rows negated once so the
             * existing "ru += -ch_c*(...)" wrapper restores the sign.
             * Default off = legacy bit-preserving path. */
            /* v1.10 DT-HIOM (2026-07-12, Jae-approved regression fix): the
             * TABLE kernel is now the DEFAULT.  The legacy l-sums above
             * cover moment arrays only; in value-kernel production
             * (water_phase_kernel==0, the default) they MISS THE HYDROSOL
             * ENTIRELY (D3-0b/D3-0c finding) - the injected skylight is
             * scattered with a Rayleigh-only phase (bb/b ~ 0.5) instead of
             * the forward-peaked hydrosol blend (bb/b ~ 0.01-0.02),
             * overdriving the sky response ~5.7x in high-omega water
             * (rrs +23..+95% vs OSOAA; production point +10%).  Measured
             * closure with the table kernel: cs1/555/sza40 atm-ON rrs
             * 1.4881e-2 - 1.1467e-2, vs v1.09 1.1506e-2 / beams 1.1494e-2
             * / OSOAA 1.1534e-2.  Pure water (true Rayleigh medium) is
             * numerically IDENTICAL under both kernels (legacy was correct
             * there - the six-digit A-FIX-2 verification stands).
             * OCRT_EXTTOP_KERNEL=legacy restores the old path (probe). */
            { const char *tk = getenv("OCRT_EXTTOP_KERNEL");
              const int tk_table = ws && ws->phase_fourier_m &&
                                   !(tk && !strcmp(tk, "legacy"));
              if (tk_table) {
                if (j > 0) {
                    sI =  ws->phase_fourier_m[+c][-j];
                    sQ =  ws->gr_pol[+c][-j];
                    sU = -ws->gt_pol[+c][-j];              /* -S_U(+j)<-I */
                    qI =  ws->gr_pol[+j][-c];
                    qQ =  ws->arr_pol[+c][-j];
                    qU = +ws->art_pol[+j][-c];             /* -S_U(+j)<-Q = -(-art) */
                    uI = -ws->gt_pol[+j][-c];
                    uQ = +ws->art_pol[+c][-j];
                    uU = -ws->att_pol[+c][-j];             /* -S_U(+j)<-U */
                } else if (j < 0) {
                    sI =  ws->phase_fourier_m[+c][-j];
                    sQ =  ws->gr_pol[+c][-j];
                    sU = -ws->gt_pol[+c][-j];              /* -S_U(-|j|)<-I */
                    qI =  ws->gr_pol[-j][+c];
                    qQ =  ws->arr_pol[+c][-j];
                    qU = -ws->art_pol[-j][+c];             /* -S_U(-|j|)<-Q */
                    uI = +ws->gt_pol[-j][+c];
                    uQ = +ws->art_pol[+c][-j];
                    uU = -ws->att_pol[+c][-j];             /* -S_U(-|j|)<-U */
                } else { sI=sQ=sU=qI=qQ=qU=uI=uQ=uU=0.0; }
                /* row convention: rowX(out I/Q) = S coeff; rowU*(out U)
                 * = -S coeff (wrapper applies the global minus).  The
                 * assignments above already fold that:
                 *   sU/qU/uU hold -(S_U coeff);  sQ/qQ/uQ hold S_Q;
                 *   sI/qI/uI hold S_I -- but careful: legacy names are
                 *   (row of OUTPUT per INCIDENT): sX = X-out <- I-in,
                 *   qX = X-out <- Q-in, uX = X-out <- U-in. */
              }
            }
            rowI [j+n_mu] = sI; rowQ [j+n_mu] = sQ; rowU [j+n_mu] = sU;
            rowIq[j+n_mu] = qI; rowQq[j+n_mu] = qQ; rowUq[j+n_mu] = qU;
            rowIu[j+n_mu] = uI; rowQu[j+n_mu] = uQ; rowUu[j+n_mu] = uU;
            /* stash Rayleigh in-place additions via closure of layer loop:
             * combine below with xdel/ydel per layer. store ray parts: */
            rowI [j+n_mu] = sI;  /* aer only; ray handled per-layer */
            (void)raI;(void)raQ;(void)raU;(void)rqI;(void)rqQ;(void)rqU;
        }
        for (int k = 0; k <= nt; ++k) {
            const double ch_c = 0.5 * exp(-atm->h[k] / mu_c);
            const double xd = atm->xdel[k], yd = atm->ydel[k];
            double *ri = src_i + (size_t)k * (size_t)dirs;
            double *rq = src_q + (size_t)k * (size_t)dirs;
            double *ru = src_u + (size_t)k * (size_t)dirs;
            const double xplc = atm->xpl[cc], xrlc = atm->xrl ? atm->xrl[cc] : 0.0;
            for (int j = -n_mu; j <= n_mu; ++j) {
                if (j == 0) continue;          /* solar slot stays poisoned */
                const double raI = beta0_m + beta2 * atm->xpl[j] * xplc;
                const double raQ = ray_on ? gamma2 * atm->xrl[j] * xplc : 0.0;
                const double raU = ray_on ? gamma2 * atm->xtl[j] * xplc : 0.0;
                const double rqI = ray_on ? gamma2 * atm->xpl[j] * xrlc : 0.0;
                const double rqQ = ray_on ? alpha2 * atm->xrl[j] * xrlc : 0.0;
                const double rqU = ray_on ? alpha2 * atm->xtl[j] * xrlc : 0.0;
                const double xtlc_u = atm->xtl ? atm->xtl[cc] : 0.0;   /* D3-0b */
                /* Stage-2 diffuse-top U-column closure:
                 * ext_top_U uses the same sine-Fourier Stokes-U convention as
                 * the water SOS field.  The existing outer minus on src_u
                 * therefore requires the three molecular U-input coefficients
                 * to be negated here, exactly matching
                 * rt_sos_operator_apply_vector for a downward incident field. */
                const double ruI = ray_on ? -gamma2 * atm->xpl[j] * xtlc_u : 0.0;
                const double ruQ = ray_on ? -alpha2 * atm->xrl[j] * xtlc_u : 0.0;
                const double ruU = ray_on ? -alpha2 * atm->xtl[j] * xtlc_u : 0.0;
                ri[j+n_mu] += ch_c * (AI * (xd*rowI [j+n_mu] + yd*raI)
                                     + AQ * (xd*rowIq[j+n_mu] + yd*rqI)
                                     + AU * (xd*rowIu[j+n_mu] + yd*ruI));
                rq[j+n_mu] += +ch_c * (AI * (xd*rowQ [j+n_mu] + yd*raQ)
                                      + AQ * (xd*rowQq[j+n_mu] + yd*rqQ)
                                      + AU * (xd*rowQu[j+n_mu] + yd*ruQ));
                ru[j+n_mu] += -ch_c * (AI * (xd*rowU [j+n_mu] + yd*raU)
                                      + AQ * (xd*rowUq[j+n_mu] + yd*rqU)
                                      + AU * (xd*rowUu[j+n_mu] + yd*ruU));
            }
        }
    }
    free(rowI);
}

/* 2026-07-14 W-PCHIP (Jae 승인): 수중 모드장(per-m)의 off-node μ 샘플링을
 * 대기 LUT-PCHIP(rt_solver.c pchip_view_at_toa)과 동일 정책으로 교체한다 —
 * monotone cubic (PCHIP, Fritsch-Carlson derivatives, Butland/scipy
 * weighted-harmonic form).  노드 일치점은 t=0/t=1 Hermite 평가로 bit 동일,
 * 노드 범위 밖(mu<=rm[1], mu>=rm[n_mu])과 n_mu<4는 기존 선형 클램프
 * 폴백(bit-parity; 기본 TWA 추가 노드 mu=1.0은 항상 이 폴백을 탄다).
 * 사전 검증(443nm/aCDOM 0.1/SZA40, GL-24 노드 -> off-node 8각,
 * view-as-node 참값 대비): 선형 최대오차 Rrs_I 0.59% / rrs0- 0.006%
 * -> PCHIP 0.009% / 0.0002% (최종량 기준 프록시 테스트).
 * 도함수 헬퍼 2종은 rt_solver.c의 static 원본과 동형 복제다(TU 분리로
 * 직접 재사용 불가; 수식 변경 시 두 곳 동시 갱신 필수). */
static double ocrt_wpchip_deriv_interior(double h0, double h1,
                                         double s0, double s1) {
    if (s0 == 0.0 || s1 == 0.0 || (s0 > 0.0) != (s1 > 0.0)) return 0.0;
    const double w1 = 2.0 * h1 + h0;
    const double w2 = h1 + 2.0 * h0;
    return (w1 + w2) / (w1 / s0 + w2 / s1);
}
static double ocrt_wpchip_deriv_edge(double h0, double h1,
                                     double s0, double s1) {
    /* scipy _edge_case: one-sided three-point estimate + monotonicity clamp */
    double d = ((2.0 * h0 + h1) * s0 - h0 * s1) / (h0 + h1);
    if (d == 0.0 || s0 == 0.0 || (d > 0.0) != (s0 > 0.0)) {
        if (s0 == 0.0 || (d > 0.0) != (s0 > 0.0)) d = 0.0;
    } else if (((s0 > 0.0) != (s1 > 0.0)) && fabs(d) > 3.0 * fabs(s0)) {
        d = 3.0 * s0;
    }
    return d;
}
/* 2026-07-14 W-PCHIP: z=0- 뷰 모드 추출(S-007 지점, 7e 블록)의 내부 구간
 * PCHIP 평가.  y(j) = tot[j + n_mu] (현재 m의 상향 양(+)노드장),
 * x(j) = atm->rm[+j] (오름차순 양의 GL 노드).  호출측이 내부 구간
 * (rm[1] < mu < rm[n_mu]) 및 n_mu >= 4를 보장한다.  도함수·Hermite 식은
 * 대기 pchip_view_at_toa와 동일하다. */
static double ocrt_wpchip_eval_up_nodes(const rt_atm_t *atm, const double *tot,
                                        int n_mu, double mu) {
    int jl = 1, jh = 2;
    for (int jp = 1; jp < n_mu; ++jp) {
        if (atm->rm[+jp] <= mu && mu < atm->rm[+jp+1]) { jl = jp; jh = jp + 1; break; }
    }
    #define OCRT_WPCHIP_UY(j) tot[(size_t)((j) + n_mu)]
    const double x1 = atm->rm[+jl];
    const double x2 = atm->rm[+jh];
    const double y1 = OCRT_WPCHIP_UY(jl);
    const double y2 = OCRT_WPCHIP_UY(jh);
    const double h  = x2 - x1;
    const double s  = (y2 - y1) / h;

    double d1, d2;
    if (jl == 1) {
        const double x3 = atm->rm[+jh + 1];
        const double s2 = (OCRT_WPCHIP_UY(jh + 1) - y2) / (x3 - x2);
        d1 = ocrt_wpchip_deriv_edge(h, x3 - x2, s, s2);
        d2 = ocrt_wpchip_deriv_interior(h, x3 - x2, s, s2);
    } else if (jh == n_mu) {
        const double x0 = atm->rm[+jl - 1];
        const double s0 = (y1 - OCRT_WPCHIP_UY(jl - 1)) / (x1 - x0);
        d1 = ocrt_wpchip_deriv_interior(x1 - x0, h, s0, s);
        d2 = ocrt_wpchip_deriv_edge(h, x1 - x0, s, s0);
    } else {
        const double x0 = atm->rm[+jl - 1];
        const double x3 = atm->rm[+jh + 1];
        const double s0 = (y1 - OCRT_WPCHIP_UY(jl - 1)) / (x1 - x0);
        const double s2 = (OCRT_WPCHIP_UY(jh + 1) - y2) / (x3 - x2);
        d1 = ocrt_wpchip_deriv_interior(x1 - x0, h, s0, s);
        d2 = ocrt_wpchip_deriv_interior(h, x3 - x2, s, s2);
    }
    #undef OCRT_WPCHIP_UY

    const double t  = (mu - x1) / h;
    const double t2 = t * t, t3 = t2 * t;
    return (2.0 * t3 - 3.0 * t2 + 1.0) * y1 + (t3 - 2.0 * t2 + t) * h * d1
         + (-2.0 * t3 + 3.0 * t2) * y2 + (t3 - t2) * h * d2;
}
/* 교체 전 원본 선형 보간(무수정 보존): 범위 밖 클램프 외삽·n_mu<4 폴백 전용. */
static double ocrt_interp_water_mode_at_mu_linear(const rt_atm_t *atm, const double *mode_nodes,
                                                  int n_mu, int mm, double mu) {
    if (!atm || !mode_nodes || n_mu <= 0) return 0.0;
    if (n_mu == 1) return mode_nodes[(size_t)mm*(size_t)n_mu];
    int jl = 1, jh = 2;
    if (mu <= atm->rm[+1]) {
        jl = 1; jh = 2;
    } else if (mu >= atm->rm[+n_mu]) {
        jl = n_mu - 1; jh = n_mu;
    } else {
        for (int jp = 1; jp < n_mu; ++jp) {
            if (atm->rm[+jp] <= mu && mu < atm->rm[+jp+1]) { jl = jp; jh = jp + 1; break; }
        }
    }
    const double ml = atm->rm[+jl];
    const double mh = atm->rm[+jh];
    const double den = mh - ml;
    const double wh = (fabs(den) > 0.0) ? (mu - ml) / den : 0.0;
    const double wl = 1.0 - wh;
    const size_t ia = (size_t)mm*(size_t)n_mu + (size_t)(jl - 1);
    const size_t ib = (size_t)mm*(size_t)n_mu + (size_t)(jh - 1);
    return wl * mode_nodes[ia] + wh * mode_nodes[ib];
}
static double ocrt_interp_water_mode_at_mu(const rt_atm_t *atm, const double *mode_nodes,
                                           int n_mu, int mm, double mu) {
    if (!atm || !mode_nodes || n_mu <= 0) return 0.0;
    /* 범위 밖·저노드는 기존 선형 경로 그대로(bit-parity). */
    if (n_mu < 4 || mu <= atm->rm[+1] || mu >= atm->rm[+n_mu])
        return ocrt_interp_water_mode_at_mu_linear(atm, mode_nodes, n_mu, mm, mu);
    int jl = 1, jh = 2;
    for (int jp = 1; jp < n_mu; ++jp) {
        if (atm->rm[+jp] <= mu && mu < atm->rm[+jp+1]) { jl = jp; jh = jp + 1; break; }
    }
    #define OCRT_WPCHIP_Y(j) mode_nodes[(size_t)mm*(size_t)n_mu + (size_t)((j) - 1)]
    const double x1 = atm->rm[+jl];
    const double x2 = atm->rm[+jh];
    const double y1 = OCRT_WPCHIP_Y(jl);
    const double y2 = OCRT_WPCHIP_Y(jh);
    const double h  = x2 - x1;
    const double s  = (y2 - y1) / h;

    double d1, d2;
    if (jl == 1) {
        const double x3 = atm->rm[+jh + 1];
        const double s2 = (OCRT_WPCHIP_Y(jh + 1) - y2) / (x3 - x2);
        d1 = ocrt_wpchip_deriv_edge(h, x3 - x2, s, s2);
        d2 = ocrt_wpchip_deriv_interior(h, x3 - x2, s, s2);
    } else if (jh == n_mu) {
        const double x0 = atm->rm[+jl - 1];
        const double s0 = (y1 - OCRT_WPCHIP_Y(jl - 1)) / (x1 - x0);
        d1 = ocrt_wpchip_deriv_interior(x1 - x0, h, s0, s);
        d2 = ocrt_wpchip_deriv_edge(h, x1 - x0, s, s0);
    } else {
        const double x0 = atm->rm[+jl - 1];
        const double x3 = atm->rm[+jh + 1];
        const double s0 = (y1 - OCRT_WPCHIP_Y(jl - 1)) / (x1 - x0);
        const double s2 = (OCRT_WPCHIP_Y(jh + 1) - y2) / (x3 - x2);
        d1 = ocrt_wpchip_deriv_interior(x1 - x0, h, s0, s);
        d2 = ocrt_wpchip_deriv_interior(h, x3 - x2, s, s2);
    }
    #undef OCRT_WPCHIP_Y

    /* Cubic Hermite on [x1, x2] — rt_solver.c pchip_view_at_toa와 동일식. */
    const double t  = (mu - x1) / h;
    const double t2 = t * t, t3 = t2 * t;
    return (2.0 * t3 - 3.0 * t2 + 1.0) * y1 + (t3 - 2.0 * t2 + t) * h * d1
         + (-2.0 * t3 + 3.0 * t2) * y2 + (t3 - t2) * h * d2;
}

static double ocrt_linear_endpoint_eval_sorted(int n_in, const double *x_in, const double *y_in, double xp) {
    if (!x_in || !y_in || n_in <= 0) return 0.0;
    if (n_in == 1) return y_in[0];
    int n = n_in;
    double *x = (double*)calloc((size_t)n, sizeof(double));
    double *y = (double*)calloc((size_t)n, sizeof(double));
    if (!x || !y) { free(x); free(y); return 0.0; }
    for (int i = 0; i < n; ++i) { x[i] = x_in[i]; y[i] = y_in[i]; }
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) if (x[i] > x[j]) {
        double tx=x[i]; x[i]=x[j]; x[j]=tx; double ty=y[i]; y[i]=y[j]; y[j]=ty;
    }
    int m = 0;
    for (int i = 0; i < n; ++i) {
        if (m == 0 || fabs(x[i] - x[m-1]) > 1e-14) { x[m]=x[i]; y[m]=y[i]; m++; }
        else y[m-1] = y[i];
    }
    n = m;
    if (n == 1) { double v=y[0]; free(x); free(y); return v; }
    int klo, khi;
    if (xp <= x[0]) { klo = 0; khi = 1; }
    else if (xp >= x[n-1]) { klo = n-2; khi = n-1; }
    else {
        int lo=0, hi=n-1;
        while (hi-lo > 1) { int mid=(hi+lo)/2; if (x[mid] > xp) hi=mid; else lo=mid; }
        klo=lo; khi=hi;
    }
    const double h=x[khi]-x[klo];
    const double t=(fabs(h)>0.0)?(xp-x[klo])/h:0.0;
    const double v=(1.0-t)*y[klo]+t*y[khi];
    free(x); free(y); return v;
}

static double ocrt_clamp_first_eval_sorted(int n_in, const double *x_in, const double *y_in,
                                           double xp, int clamped_end_slopes) {
    if (!x_in || !y_in || n_in <= 0) return 0.0;
    double xmin = x_in[0], ymin = y_in[0], xmax = x_in[0], ymax = y_in[0];
    for (int i = 1; i < n_in; ++i) {
        if (x_in[i] < xmin) { xmin = x_in[i]; ymin = y_in[i]; }
        if (x_in[i] > xmax) { xmax = x_in[i]; ymax = y_in[i]; }
    }
    if (xp <= xmin) return ymin;
    if (xp >= xmax) return ymax;
    return ocrt_cubic_spline_eval_sorted(n_in, x_in, y_in, xp, clamped_end_slopes);
}

static double ocrt_cubic_spline_eval_sorted(int n_in, const double *x_in, const double *y_in,
                                            double xp, int clamped_end_slopes) {
    if (!x_in || !y_in || n_in <= 0) return 0.0;
    if (n_in == 1) return y_in[0];

    /* OSOAA's SOS_INTERPO_SPLINT explicitly sorts X/Y before building the
     * spline.  Keep the same behaviour here. */
    int n = n_in;
    double *x = (double*)calloc((size_t)n, sizeof(double));
    double *y = (double*)calloc((size_t)n, sizeof(double));
    if (!x || !y) { free(x); free(y); return 0.0; }
    for (int i = 0; i < n; ++i) { x[i] = x_in[i]; y[i] = y_in[i]; }
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (x[i] > x[j]) {
                double tx = x[i]; x[i] = x[j]; x[j] = tx;
                double ty = y[i]; y[i] = y[j]; y[j] = ty;
            }
        }
    }
    /* Drop duplicate or non-increasing abscissae defensively. */
    int m = 0;
    for (int i = 0; i < n; ++i) {
        if (m == 0 || fabs(x[i] - x[m-1]) > 1e-14) {
            x[m] = x[i]; y[m] = y[i]; m++;
        } else {
            y[m-1] = y[i];
        }
    }
    n = m;

    const int endpoint_mode = ocrt_water_twa_critical_endpoint_mode();
    if (endpoint_mode == 4 && n >= 1 && x[0] > 1e-14) {
        /* Diagnostic critical-edge rule: include the physical horizon /
         * critical endpoint (mu_air=0, transmitted radiance=0) before
         * applying the same OSOAA-compatible clamped cubic SPLINT. */
        double *xz = (double*)calloc((size_t)(n + 1), sizeof(double));
        double *yz = (double*)calloc((size_t)(n + 1), sizeof(double));
        if (xz && yz) {
            xz[0] = 0.0; yz[0] = 0.0;
            for (int ii = 0; ii < n; ++ii) { xz[ii + 1] = x[ii]; yz[ii + 1] = y[ii]; }
            free(x); free(y);
            x = xz; y = yz; n = n + 1;
        } else {
            free(xz); free(yz);
        }
    }

    if (n == 1) { double v = y[0]; free(x); free(y); return v; }
    if (endpoint_mode != 0) {
        if (xp <= x[0]) {
            double v = y[0];
            if (endpoint_mode == 1 && n >= 2) {
                const double h = x[1] - x[0];
                const double t = (fabs(h) > 0.0) ? (xp - x[0]) / h : 0.0;
                v = (1.0 - t) * y[0] + t * y[1];
            } else if (endpoint_mode == 2) {
                v = y[0];
            } else if (endpoint_mode == 3) {
                v = (x[0] > 1e-14) ? y[0] * (xp / x[0]) : y[0];
            }
            if (endpoint_mode != 4) { free(x); free(y); return v; }
        } else if (xp >= x[n-1]) {
            double v = y[n-1];
            if (endpoint_mode == 1 && n >= 2) {
                const double h = x[n-1] - x[n-2];
                const double t = (fabs(h) > 0.0) ? (xp - x[n-2]) / h : 1.0;
                v = (1.0 - t) * y[n-2] + t * y[n-1];
            } else if (endpoint_mode == 2) {
                v = y[n-1];
            }
            if (endpoint_mode != 4) { free(x); free(y); return v; }
        }
    }

    if (n == 2) {
        const double h = x[1] - x[0];
        const double t = (fabs(h) > 0.0) ? (xp - x[0]) / h : 0.0;
        double v = (1.0 - t) * y[0] + t * y[1];
        free(x); free(y); return v;
    }

    double *y2 = (double*)calloc((size_t)n, sizeof(double));
    double *u  = (double*)calloc((size_t)(n - 1), sizeof(double));
    if (!y2 || !u) {
        int klo = 0, khi = n - 1;
        if (xp <= x[0]) { klo = 0; khi = 1; }
        else if (xp >= x[n-1]) { klo = n-2; khi = n-1; }
        else {
            int lo = 0, hi = n - 1;
            while (hi - lo > 1) {
                int mid = (hi + lo) / 2;
                if (x[mid] > xp) hi = mid; else lo = mid;
            }
            klo = lo; khi = hi;
        }
        const double h = x[khi] - x[klo];
        const double t = (fabs(h) > 0.0) ? (xp - x[klo]) / h : 0.0;
        double v = (1.0 - t) * y[klo] + t * y[khi];
        free(y2); free(u); free(x); free(y); return v;
    }

    if (clamped_end_slopes) {
        const double h0 = x[1] - x[0];
        const double yp1 = (fabs(h0) > 0.0) ? (y[1] - y[0]) / h0 : 0.0;
        y2[0] = -0.5;
        u[0] = (fabs(h0) > 0.0) ? (3.0 / h0) * ((y[1] - y[0]) / h0 - yp1) : 0.0;
    } else {
        y2[0] = 0.0;
        u[0] = 0.0;
    }

    for (int i = 1; i < n - 1; ++i) {
        const double him1 = x[i] - x[i-1];
        const double hi   = x[i+1] - x[i];
        const double den  = x[i+1] - x[i-1];
        if (fabs(him1) <= 0.0 || fabs(hi) <= 0.0 || fabs(den) <= 0.0) {
            y2[i] = 0.0;
            u[i]  = 0.0;
            continue;
        }
        const double sig = him1 / den;
        const double p   = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        const double dydx_hi   = (y[i+1] - y[i]) / hi;
        const double dydx_him1 = (y[i] - y[i-1]) / him1;
        u[i] = (6.0 * (dydx_hi - dydx_him1) / den - sig * u[i-1]) / p;
    }

    if (clamped_end_slopes) {
        const double hn = x[n-1] - x[n-2];
        const double ypn = (fabs(hn) > 0.0) ? (y[n-1] - y[n-2]) / hn : 0.0;
        const double qn = 0.5;
        const double un = (fabs(hn) > 0.0) ? (3.0 / hn) * (ypn - (y[n-1] - y[n-2]) / hn) : 0.0;
        y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);
    } else {
        y2[n-1] = 0.0;
    }

    for (int k = n - 2; k >= 0; --k)
        y2[k] = y2[k] * y2[k+1] + u[k];

    int klo = 0, khi = n - 1;
    if (xp <= x[0]) { klo = 0; khi = 1; }
    else if (xp >= x[n-1]) { klo = n-2; khi = n-1; }
    else {
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            if (x[mid] > xp) hi = mid; else lo = mid;
        }
        klo = lo; khi = hi;
    }
    const double h = x[khi] - x[klo];
    double val;
    if (fabs(h) <= 0.0) {
        val = y[klo];
    } else {
        const double a = (x[khi] - xp) / h;
        const double b = (xp - x[klo]) / h;
        val = a*y[klo] + b*y[khi] + ((a*a*a - a)*y2[klo] + (b*b*b - b)*y2[khi]) * (h*h) / 6.0;
    }
    free(y2); free(u); free(x); free(y);
    return val;
}

/* Diagnostic-only primary source prefactor dump.  Unlike the plain source trace,
 * this reports phase-like source coefficients by dividing the source by ch*xdel
 * for hydrosol-only water.  This is meant to compare against OSOAA source-pref
 * BP trends and is not a production output path. */
static void dump_water_primary_pref_trace_if_requested(int m,
                                                       const rt_atm_t *atm,
                                                       const double *J_I,
                                                       const double *J_Q,
                                                       const double *J_U) {
    const char *path = g_wrt_env.s_OCRT_DUMP_WATER_SOURCE_PREF_TRACE;
    if (!path || !path[0] || !atm || !J_I || !J_Q || !J_U) return;
    FILE *fp = fopen(path, "a+");
    if (!fp) return;
    fseek(fp, 0, SEEK_END);
    long file_len = ftell(fp);
    if (file_len == 0L) {
        fprintf(fp, "solver,stage,m,order,k,j_signed,mu,ch,xdel,ydel,beam_q,phase_I,phase_Q,phase_U,J_I,J_Q,J_U,wavelength_nm\n");
    }
    double wl = NAN;
    const char *wls = g_wrt_env.s_OCRT_TRACE_WAVELENGTH_NM;
    if (wls && wls[0]) wl = atof(wls);
    const int nt = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    for (int ksel_i = 0; ksel_i < 3; ++ksel_i) {
        int k = (ksel_i == 0) ? 0 : ((ksel_i == 1) ? (nt / 2) : nt);
        if (k < 0) k = 0;
        if (k > nt) k = nt;
        const double ch = atm->ch ? atm->ch[k] : NAN;
        const double x  = atm->xdel ? atm->xdel[k] : NAN;
        const double y  = atm->ydel ? atm->ydel[k] : NAN;
        const double denom = ch * x;
        for (int j = -n_mu; j <= n_mu; ++j) {
            if (j == 0) continue;
            size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
            const double pI = (fabs(denom) > 1e-300) ? (J_I[idx] / denom) : NAN;
            const double pQ = (fabs(denom) > 1e-300) ? (J_Q[idx] / denom) : NAN;
            const double pU = (fabs(denom) > 1e-300) ? (J_U[idx] / denom) : NAN;
            fprintf(fp, "water,primary_pref,%d,1,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17g\n",
                    m, k, j, atm->rm[j], ch, x, y, atm->beam_q,
                    pI, pQ, pU, J_I[idx], J_Q[idx], J_U[idx], wl);
        }
    }
    fclose(fp);
}

/* Diagnostic-only primary source dump for the in-water SOS entry point. */
static void dump_water_primary_source_trace_if_requested(int m,
                                                        const rt_atm_t *atm,
                                                        const double *J_I,
                                                        const double *J_Q,
                                                        const double *J_U) {
    const char *path = g_wrt_env.s_OCRT_DUMP_WATER_SOURCE_TRACE;
    if (!path || !path[0] || !atm || !J_I || !J_Q || !J_U) return;
    FILE *fp = fopen(path, "a+");
    if (!fp) return;
    fseek(fp, 0, SEEK_END);
    long file_len = ftell(fp);
    const int nt = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    if (file_len == 0L) {
        fprintf(fp, "solver,stage,m,order,k,j_signed,mu,J_I,J_Q,J_U\n");
    }
    for (int ksel_i = 0; ksel_i < 3; ++ksel_i) {
        int k = (ksel_i == 0) ? 0 : ((ksel_i == 1) ? (nt / 2) : nt);
        if (k < 0) k = 0;
        if (k > nt) k = nt;
        for (int j = -n_mu; j <= n_mu; ++j) {
            if (j == 0) continue;
            size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
            fprintf(fp, "water,primary_source,%d,1,%d,%d,%.17g,%.17e,%.17e,%.17e\n",
                    m, k, j, atm->rm[j], J_I[idx], J_Q[idx], J_U[idx]);
        }
    }
    fclose(fp);
}

/* ---------------------------------------------------------------------------
 * Defaults
 * --------------------------------------------------------------------------- */

rt_water_rt_options_t rt_water_rt_options_default(void) {
    rt_water_rt_options_t o = {
        .max_iterations  = 20,
        .tolerance       = 1.0e-7,
        .n_layers_water  = 0,       /* auto production grid from final τ and Δτ target; debug-only override */
        .tau_max_target  = -1.0,    /* auto base optical depth; debug-only override */
        .max_tau_max_target = 0.0,  /* no production cap unless particle/deep policy sets one */
        .max_z_max_m     = 0.0,     /* no physical-depth ceiling unless particle/deep policy sets one */
        .depth_bottom_tol = 1.0e-8, /* target exp(-transport_ext*z) bottom attenuation */
        .layer_dtau_target = 0.05,  /* production vertical grid target Δτ per layer */
        .n_mu_water      = 48,  /* in-water angular quadrature (positive Gauss nodes).
                                 * RECOMMENDED >= 48 for forward-peaked / high-omega particle
                                 * phases (mineral, phytoplankton) and for LUT generation
                                 * (--output-full-grid).  Convergence study, Brown_earth
                                 * Csed=5 g/m3 @555nm vs n_mu=64 converged value:
                                 *   n_mu=24 +3.0%, n_mu=32 +1.5%, n_mu=48 +0.2%.
                                 * The n_mu floor is case-dependent and rises with single-
                                 * scattering albedo, so high-omega bands may need >48 --
                                 * verify per-case convergence for production LUTs.
                                 * Was 24 (atm-n_mu-matched) through v1.08; raised 2026-06-30
                                 * for in-water accuracy.  NOTE: frozen baselines/regression
                                 * guards that relied on the old default 24 must be re-frozen
                                 * at 48.  Override with --n-mu-water. */
        .view_as_node    = 0,
        .view_vza_deg_list = NULL,   /* #20 */
        .n_view_vza      = 0,
        .m_max_water     = 2,
        .q_convention    = 1,
        /* CDOM 기본값: 비활성 (a_CDOM=0). Slope/ref λ 는 활성화 시 사용. */
        .a_cdom_440_m_inv   = 0.0,
        .S_cdom_nm_inv      = 0.014,
        .cdom_ref_lambda_nm = 440.0,
        .water_input_mode   = RT_WATER_INPUT_UNSET,
        .ccrr_mode          = 0,
        .water_constituent_model = RT_WATER_CONSTITUENT_OCRT,
        .ccrr_chl_mg_m3     = 0.0,
        .organic_phyto_scattering = 0,
        .ccrr_min_g_m3      = 0.0,
        .tsm_species         = AHN_RED_CLAY,
        .organic_phyto_group = ORGANIC_PHYTO_MICRO,
        .detritus_a440_m_inv = 0.0,
        .detritus_slope_nm_inv = ORGANIC_DETRITUS_SLOPE_DEFAULT,
        .ccrr_phase_moments_path = NULL,
        .ccrr_particle_phase_lut_path = NULL,
        .ccrr_particle_phase_case_id = NULL,
        .ccrr_particle_phase_wavelength_nm = 0.0,
        .ccrr_particle_phase_lmax = 10,
        .ccrr_particle_phase_nphi = 720,
        .fixed_bulk_iop_mode = 0,
        .fixed_a_total_m_inv = 0.0,
        .fixed_b_total_m_inv = 0.0,
        .fixed_bb_total_m_inv = 0.0,
        .fixed_bulk_lmax = 0,     /* 0 = per-path default (#23): water-mie 200, LUT 30 */
        .fixed_bulk_phase_model = 0,
        .fixed_bulk_phase_nphi = 720,
        .fixed_bulk_ff_n = 1.18,
        .fixed_bulk_ff_mu = -1.0,
        .fixed_bulk_phase_lut_path = NULL,
        .water_mie_phase_path = NULL,
        .water_mie_moment_mode = 1,
        .water_mie_truncation_mode = 0,
        .water_mie_ss_mode = 0,
        .water_mie_moment_n_mu = 60,
        .fixed_bulk_phase_case_id = NULL,
        .fixed_bulk_phase_wavelength_nm = 0.0,
        .wind_speed = 0.0,
        .cox_munk_sigma_type = 1,
        .water_phase_kernel = 0   /* path B (value kernel + OSOAA cap) by default */
    };
    return o;
}

/* ---------------------------------------------------------------------------
 * CCRR/IOCCG21 scalar particle phase coefficient capacity.
 *
 * Particle coefficients are supplied by rt_iop_ccrr_particle_betal().  That
 * function may use built-in regression arrays or a CCRR-generated
 * phase_moments CSV loaded via --ccrr-phase-moments.  Pure seawater keeps the
 * OCRT vector Rayleigh-like phase path and is not overridden by scalar CCRR
 * water moments.
 * --------------------------------------------------------------------------- */
#define CCRR_PARTICLE_BETAL_LMAX 24
#define FIXED_BULK_BETAL_LMAX 200

/* Henyey-Greenstein Legendre moments for a prescribed backscatter ratio.
 * Exact HG moment identity: g_l = g^l, hence in the (2l+1)-scaled
 * convention used by the kernels:
 *   beta_l = (2l+1) * g^l
 * with g obtained by inverting bb/b(g) (closed form in rt_water_iop.c). */
static int hg_betal_for_target_bb(double bb_over_b, int lmax, double *betal_out)
{
    if (!betal_out || lmax < 0 || lmax > FIXED_BULK_BETAL_LMAX) return -1;
    if (!(bb_over_b > 0.0 && bb_over_b < 0.5)) return -2;
    double g = rt_iop_hg_g_for_backscatter_fraction(bb_over_b);
    double gp = 1.0;
    for (int l = 0; l <= lmax; ++l) {
        betal_out[l] = (2.0 * (double)l + 1.0) * gp;
        gp *= g;
    }
    return 0;
}


/* Henyey-Greenstein phase function, mu = cos(scattering angle):
 *   P_HG(mu; g) = (1 - g^2) / (1 + g^2 - 2 g mu)^{3/2}
 * normalized to (1/2) int_{-1}^{1} P dmu = 1. */
static double hg_phase_p11(double mu, double g)
{
    double den = 1.0 + g*g - 2.0*g*mu;
    if (den <= 1.0e-300) den = 1.0e-300;
    return (1.0 - g*g) / (den * sqrt(den));
}

typedef struct { int n; double *theta_deg; double *p11; } fixed_bulk_phase_table_t;
typedef struct { double t, p; } phase_pair_t;
static int phase_pair_cmp(const void *a, const void *b) {
    const phase_pair_t *pa=(const phase_pair_t*)a, *pb=(const phase_pair_t*)b;
    return (pa->t < pb->t) ? -1 : (pa->t > pb->t);
}
static int split_csv_local(char *line, char **fields, int maxf) {
    int n=0; char *p=line;
    while (p && *p && n<maxf) { fields[n++]=p; char *c=strchr(p, ','); if(!c) break; *c='\0'; p=c+1; }
    for (int i=0;i<n;i++){ char *q=fields[i]; while(*q==' '||*q=='\t'||*q=='\r'||*q=='\n') q++; fields[i]=q; size_t L=strlen(q); while(L>0&&(q[L-1]==' '||q[L-1]=='\t'||q[L-1]=='\r'||q[L-1]=='\n')) q[--L]='\0'; }
    return n;
}
static int csv_find_col(char **names, int n, const char *target) {
    for(int i=0;i<n;i++) if(names[i] && !strcmp(names[i], target)) return i; return -1;
}
static void fixed_bulk_phase_table_free(fixed_bulk_phase_table_t *tab){ if(!tab)return; free(tab->theta_deg); free(tab->p11); tab->theta_deg=NULL; tab->p11=NULL; tab->n=0; }
static int fixed_bulk_phase_lut_load(const char *path, const char *case_id, double wl_nm, fixed_bulk_phase_table_t *tab)
{
    if(!path||!*path||!tab) return -1; memset(tab,0,sizeof(*tab)); FILE *fp=fopen(path,"r"); if(!fp) return -2;
    char line[65536]; if(!fgets(line,sizeof line,fp)){fclose(fp);return -3;} char *cols[64]; int nc=split_csv_local(line,cols,64);
    int itheta=csv_find_col(cols,nc,"theta_deg"); if(itheta<0) itheta=csv_find_col(cols,nc,"theta"); int ip11=csv_find_col(cols,nc,"P11"); if(ip11<0) ip11=csv_find_col(cols,nc,"p11"); int icase=csv_find_col(cols,nc,"case_id"); int iwl=csv_find_col(cols,nc,"wavelength_nm"); if(itheta<0||ip11<0){fclose(fp);return -4;}
    int cap=1024,n=0; phase_pair_t *pairs=(phase_pair_t*)calloc((size_t)cap,sizeof(*pairs)); if(!pairs){fclose(fp);return -5;}
    while(fgets(line,sizeof line,fp)){ if(!line[0]||line[0]=='\n'||line[0]=='\r') continue; char *f[64]; int nf=split_csv_local(line,f,64); if(nf<=itheta||nf<=ip11) continue; if(icase>=0&&case_id&&*case_id){ if(nf<=icase||strcmp(f[icase],case_id)) continue; } if(iwl>=0&&wl_nm>0.0){ if(nf<=iwl) continue; double wl=atof(f[iwl]); if(fabs(wl-wl_nm)>1e-6) continue; } double t=atof(f[itheta]), p=atof(f[ip11]); if(!(t>=0.0&&t<=180.0)||!(p>=0.0)) continue; if(n>=cap){cap*=2; phase_pair_t *np=(phase_pair_t*)realloc(pairs,(size_t)cap*sizeof(*pairs)); if(!np){free(pairs);fclose(fp);return -5;} pairs=np;} pairs[n].t=t; pairs[n].p=p; n++; }
    fclose(fp); if(n<2){free(pairs);return -6;} qsort(pairs,(size_t)n,sizeof(*pairs),phase_pair_cmp); tab->theta_deg=(double*)calloc((size_t)n,sizeof(double)); tab->p11=(double*)calloc((size_t)n,sizeof(double)); if(!tab->theta_deg||!tab->p11){fixed_bulk_phase_table_free(tab);free(pairs);return -5;} tab->n=n; for(int i=0;i<n;i++){tab->theta_deg[i]=pairs[i].t; tab->p11[i]=pairs[i].p;} free(pairs); return 0;
}
static double fixed_bulk_phase_lut_eval(const fixed_bulk_phase_table_t *tab, double mu)
{
    if (!tab || tab->n < 1) return 1.0;
    if (mu >  1.0) mu =  1.0;
    if (mu < -1.0) mu = -1.0;

    const double theta = acos(mu) * 180.0 / M_PI;
    if (theta <= tab->theta_deg[0]) return tab->p11[0];
    if (theta >= tab->theta_deg[tab->n - 1]) return tab->p11[tab->n - 1];

    int lo = 0, hi = tab->n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (tab->theta_deg[mid] <= theta) lo = mid;
        else hi = mid;
    }

    const double t0 = tab->theta_deg[lo];
    const double t1 = tab->theta_deg[hi];
    const double p0 = tab->p11[lo];
    const double p1 = tab->p11[hi];

    if (p0 > 0.0 && p1 > 0.0) {
        double f;
        /* Strongly forward-peaked FF / hydrosol phase functions are much
         * smoother in log(theta)-log(P11) space.  Keep theta=0 as a true
         * endpoint, and use log-theta interpolation only for strictly
         * positive brackets. */
        if (theta > 0.0 && t0 > 0.0 && t1 > 0.0) {
            f = (log(theta) - log(t0)) / (log(t1) - log(t0));
        } else {
            f = (t1 > t0) ? ((theta - t0) / (t1 - t0)) : 0.0;
        }
        if (f < 0.0) f = 0.0;
        if (f > 1.0) f = 1.0;
        return exp(log(p0) * (1.0 - f) + log(p1) * f);
    }

    const double f = (t1 > t0) ? ((theta - t0) / (t1 - t0)) : 0.0;
    return p0 * (1.0 - f) + p1 * f;
}

/* (1/2) * Integral_{lo}^{hi} P11(theta) sin(theta) dtheta on the table grid.
 * Segment endpoints are evaluated with the SAME log-log rule as
 * fixed_bulk_phase_lut_eval (consistency requirement), partial segments at
 * the range boundaries included; trapezoid per segment (the 0.25 factor is
 * 0.5 trapezoid x 0.5 normalization).  Full range [0,180] of a normalized
 * phase gives 1. */
static double fixed_bulk_phase_table_integral_theta(const fixed_bulk_phase_table_t *tab,
                                                   double theta_lo_deg,
                                                   double theta_hi_deg)
{
    if (!tab || tab->n < 2 || theta_hi_deg <= theta_lo_deg) return 0.0;
    double acc = 0.0;
    for (int i = 0; i < tab->n - 1; ++i) {
        double t0 = tab->theta_deg[i];
        double t1 = tab->theta_deg[i + 1];
        if (t1 <= theta_lo_deg || t0 >= theta_hi_deg) continue;
        double a = fmax(t0, theta_lo_deg);
        double b = fmin(t1, theta_hi_deg);
        if (b <= a) continue;
        double f0 = (t1 > t0) ? ((a - t0) / (t1 - t0)) : 0.0;
        double f1 = (t1 > t0) ? ((b - t0) / (t1 - t0)) : 0.0;
        double p0 = tab->p11[i];
        double p1 = tab->p11[i + 1];
        double pa, pb;
        if (p0 > 0.0 && p1 > 0.0) {
            pa = exp(log(p0) * (1.0 - f0) + log(p1) * f0);
            pb = exp(log(p0) * (1.0 - f1) + log(p1) * f1);
        } else {
            pa = p0 * (1.0 - f0) + p1 * f0;
            pb = p0 * (1.0 - f1) + p1 * f1;
        }
        double ar = a * M_PI / 180.0;
        double br = b * M_PI / 180.0;
        acc += 0.25 * (pa * sin(ar) + pb * sin(br)) * (br - ar);
    }
    return acc;
}


/* Fournier-Forand backscatter fraction, closed form (Mobley, Sundman &
 * Boss 2002, App. Opt. 41; original phase: Fournier & Forand 1994):
 *
 *   nu     = (3 - mu_J) / 2
 *   d90    = 4 / (3 (n-1)^2) * sin^2(45 deg) = 2 / (3 (n-1)^2)
 *   bb/b   = 1 - [ 1 - d90^{nu+1} - 0.5 (1 - d90^{nu}) ] / [ (1 - d90) d90^{nu} ]
 *
 * Symbols:
 *   nrel = n = particle refractive index RELATIVE to water (mineral ~1.18,
 *          pigment ~1.045 via pigment_ff_nrel below)
 *   mu_J = Junge (hyperbolic) particle-size-distribution slope, valid
 *          3 < mu_J < 5 (nu in (-1, 0)) */
static double ff_bb_over_b_from_mu(double nrel, double mu_junge)
{
    if (!(nrel > 1.0) || !(mu_junge > 3.0)) return NAN;
    double nu = (3.0 - mu_junge) / 2.0;
    double d90 = 2.0 / (3.0 * (nrel - 1.0) * (nrel - 1.0));
    double d90nu = pow(d90, nu);
    double num = 1.0 - pow(d90, nu + 1.0) - 0.5 * (1.0 - d90nu);
    double den = (1.0 - d90) * d90nu;
    return 1.0 - num / den;
}

static double ff_solve_mu_for_bb(double nrel, double target)
{
    if (!(target > 0.0 && target < 0.5) || !(nrel > 1.0)) return NAN;
    double lo = 3.0001, hi = 4.9990;
    double flo = ff_bb_over_b_from_mu(nrel, lo) - target;
    double fhi = ff_bb_over_b_from_mu(nrel, hi) - target;
    if (!isfinite(flo) || !isfinite(fhi) || flo * fhi > 0.0) return NAN;
    for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (lo + hi);
        double fm = ff_bb_over_b_from_mu(nrel, mid) - target;
        if (!isfinite(fm)) return NAN;
        if (fabs(fm) < 1e-12) return mid;
        if (flo * fm <= 0.0) { hi = mid; fhi = fm; }
        else { lo = mid; flo = fm; }
    }
    return 0.5 * (lo + hi);
}

/* ============================================================================
 * v1.09 commit #9 (PERFORMANCE, bit-identical): cross-beam memo for
 * ff_solve_mu_for_bb.
 *
 * WHY: rt_water_rt_sos_pure() is re-entered once per incident beam (solar +
 * diffuse-sky quadrature), and its pre-m-loop IOP block re-runs
 * ff_solve_mu_for_bb(nrel, bb/b) every time.  Both arguments are fixed by the
 * case IOPs, so all re-entries solve the SAME 100-step bisection, each step
 * re-integrating an FF phase table (~4e6 fixed_bulk_phase_lut_eval calls per
 * solve; 1.24e8 per case measured, ~38% of runtime post-commit-#8).
 *
 * BIT-EXACTNESS: pure function of (nrel, target); the memo returns the very
 * double produced by the first evaluation, so every consumer sees a value
 * bit-identical to an uncached call.  Thread safety follows the g_wpkc
 * pattern (threadprivate; batch OMP threads each keep their own memo).
 * ========================================================================== */
#define FF_MU_MEMO_SLOTS 8
static struct { double nrel, target, mu; int valid; } g_ff_mu_memo[FF_MU_MEMO_SLOTS];
static int g_ff_mu_memo_n = 0;
#pragma omp threadprivate(g_ff_mu_memo, g_ff_mu_memo_n)

static double ff_solve_mu_for_bb_cached(double nrel, double target)
{
    for (int i = 0; i < g_ff_mu_memo_n; ++i)
        if (g_ff_mu_memo[i].valid &&
            g_ff_mu_memo[i].nrel == nrel && g_ff_mu_memo[i].target == target)
            return g_ff_mu_memo[i].mu;
    double mu = ff_solve_mu_for_bb(nrel, target);
    int slot = (g_ff_mu_memo_n < FF_MU_MEMO_SLOTS) ? g_ff_mu_memo_n++ : (FF_MU_MEMO_SLOTS - 1);
    g_ff_mu_memo[slot].nrel = nrel; g_ff_mu_memo[slot].target = target;
    g_ff_mu_memo[slot].mu = mu;     g_ff_mu_memo[slot].valid = 1;
    return mu;
}

/* Fournier-Forand phase function (Fournier & Forand 1994 as given in
 * Mobley-Sundman-Boss 2002), scattering angle psi:
 *
 *   delta      = 4 / (3 (n-1)^2) * sin^2(psi/2),   delta_180 = 4/(3(n-1)^2)
 *   nu         = (3 - mu_J) / 2
 *   beta(psi)  = 1 / (4 pi (1-delta)^2 delta^nu) *
 *                [ nu(1-delta) - (1-delta^nu)
 *                  + ( delta(1-delta^nu) - nu(1-delta) ) / sin^2(psi/2) ]
 *              + (1 - delta_180^nu) / (16 pi (delta_180 - 1) delta_180^nu)
 *                * (3 cos^2 psi - 1)
 *
 * Returned as P11 = 4 pi * beta, i.e. normalized so (1/2)IntP sin = 1
 * analytically; the caller renormalizes NUMERICALLY on the actual grid
 * anyway (generate_ff), so small quadrature error does not leak into the
 * scattering budget.  psi is floored at 1e-9 deg (removable singularity
 * handling of the 1/sin^2 term). */
static double ff_phase_p11_dimensionless(double theta_deg, double nrel, double mu_junge)
{
    if (theta_deg < 1.0e-9) theta_deg = 1.0e-9;
    if (theta_deg > 180.0) theta_deg = 180.0;
    double psi = theta_deg * M_PI / 180.0;
    double nu = (3.0 - mu_junge) / 2.0;
    double s2 = sin(0.5 * psi);
    s2 *= s2;
    double K = 4.0 / (3.0 * (nrel - 1.0) * (nrel - 1.0));
    double delta = K * s2;
    double delta180 = K;
    const double eps = 1.0e-300;
    double s2s = (s2 > eps) ? s2 : eps;
    double ds = (delta > eps) ? delta : eps;
    double dnu = pow(ds, nu);
    double delta_nu = pow(ds, nu);
    double pref = 1.0 / (4.0 * M_PI * (1.0 - ds) * (1.0 - ds) * dnu);
    double bracket = (nu * (1.0 - delta) - (1.0 - delta_nu)
        + (delta * (1.0 - delta_nu) - nu * (1.0 - delta)) / s2s);
    double term1 = pref * bracket;
    double d180_nu = pow(delta180, nu);
    double term2 = ((1.0 - d180_nu) /
        (16.0 * M_PI * (delta180 - 1.0) * d180_nu) *
        (3.0 * cos(psi) * cos(psi) - 1.0));
    double beta_sr = term1 + term2;
    double p11 = 4.0 * M_PI * beta_sr;
    return (isfinite(p11) && p11 > 0.0) ? p11 : 0.0;
}

/* v1.053: pigment (phytoplankton) relative refractive index n_rel(lambda) for the
 * pigment Fournier-Forand phase, replacing the shared mineral-like value 1.18.
 * 3-term Cauchy fit to Phytoplankton.inp visible data (n_abs/n_water, n_water=1.34),
 * size-class-independent: n_rel = A + B/lam^2 + C/lam^4 (lam in nm). Region rule:
 *   lambda <= 694 nm : evaluate Cauchy (fit-interpolation 350-694, fit-extrapolation <350)
 *   lambda  > 694 nm : constant = Cauchy(694) (NIR hold; pigment index featureless there)
 * Pigment n_rel ~ 1.045 (soft particle), distinct from mineral n_rel = 1.18. */
static double pigment_ff_nrel(double lambda_nm)
{
    const double A = 1.04741108, B = -2.64269889e3, C = 5.14760568e8;
    double lam = (lambda_nm <= 694.0) ? lambda_nm : 694.0;
    double x = 1.0 / (lam * lam);
    return A + B * x + C * x * x;
}

static int fixed_bulk_phase_table_generate_ff(double nrel, double mu_junge, fixed_bulk_phase_table_t *tab)
{
    if (!tab || !(nrel > 1.0) || !(mu_junge > 3.0)) return -1;
    memset(tab, 0, sizeof(*tab));
    const int n1 = 500, n2 = 1900, n3 = 3200;
    const int n = 1 + n1 + n2 + n3 + 1;
    tab->theta_deg = (double*)calloc((size_t)n, sizeof(double));
    tab->p11 = (double*)calloc((size_t)n, sizeof(double));
    if (!tab->theta_deg || !tab->p11) { fixed_bulk_phase_table_free(tab); return -2; }
    int k = 0;
    tab->theta_deg[k++] = 1.0e-6;
    for (int i = 0; i < n1; ++i) {
        double f = (double)i / (double)(n1 - 1);
        double t = pow(10.0, -5.0 + f * 5.0);
        if (k > 0 && fabs(t - tab->theta_deg[k-1]) < 1e-12) continue;
        tab->theta_deg[k++] = t;
    }
    for (int i = 1; i <= n2; ++i) tab->theta_deg[k++] = 1.0 + (20.0 - 1.0) * (double)i / (double)n2;
    for (int i = 1; i <= n3; ++i) tab->theta_deg[k++] = 20.0 + (180.0 - 20.0) * (double)i / (double)n3;
    if (tab->theta_deg[k-1] < 180.0) tab->theta_deg[k++] = 180.0;
    tab->n = k;
    for (int i = 0; i < tab->n; ++i) tab->p11[i] = ff_phase_p11_dimensionless(tab->theta_deg[i], nrel, mu_junge);
    double norm = fixed_bulk_phase_table_integral_theta(tab, 0.0, 180.0);
    if (!(norm > 0.0)) { fixed_bulk_phase_table_free(tab); return -3; }
    for (int i = 0; i < tab->n; ++i) tab->p11[i] /= norm;
    return 0;
}

static double fixed_bulk_legendre_P(int l, double x)
{
    if (l <= 0) return 1.0;
    if (l == 1) return x;
    double p0 = 1.0, p1 = x, p = x;
    for (int n = 2; n <= l; ++n) {
        p = ((2.0*n - 1.0)*x*p1 - (n - 1.0)*p0) / (double)n;
        p0 = p1; p1 = p;
    }
    return p;
}

/* NORMALIZED Legendre moment g_l = (1/2) Int P11(mu) P_l(mu) dmu of a
 * (unit-normalized) phase table.  Midpoint rule per theta segment with the
 * GEOMETRIC mean of the endpoint P11 (accurate for the log-space-smooth FF
 * forward peak on the log-refined grid).  Note this returns g_l, NOT the
 * (2l+1)-scaled chi_l; the delta-M f below consumes it directly (f = g_L,
 * Wiscombe convention). */
static double fixed_bulk_phase_table_legendre_moment(const fixed_bulk_phase_table_t *tab, int l)
{
    if (!tab || tab->n < 2 || l < 0) return 0.0;
    double acc = 0.0;
    for (int i = 0; i < tab->n - 1; ++i) {
        double t0 = tab->theta_deg[i];
        double t1 = tab->theta_deg[i+1];
        double p0 = tab->p11[i];
        double p1 = tab->p11[i+1];
        double tm = 0.5*(t0+t1);
        double pm = (p0 > 0.0 && p1 > 0.0) ? sqrt(p0*p1) : 0.5*(p0+p1);
        double mu = cos(tm*M_PI/180.0);
        double Pl = fixed_bulk_legendre_P(l, mu);
        double a = t0*M_PI/180.0, b = t1*M_PI/180.0;
        acc += 0.5 * pm * Pl * sin(tm*M_PI/180.0) * (b-a);
    }
    return acc;
}

static double fixed_bulk_phase_table_cumulative_forward(const fixed_bulk_phase_table_t *tab,
                                                        double theta_hi_deg)
{
    return fixed_bulk_phase_table_integral_theta(tab, 0.0, theta_hi_deg);
}

/* v1.054: accurate normalized Legendre moments chi_l (l=0..L) of an analytic
 * phase table.  chi_l = (1/2) integral_0^pi P11(theta) P_l(cos theta) sin theta dtheta,
 * trapezoidal on the (log-refined-forward) theta grid, normalized so chi_0 = 1.
 * These feed a moment-based (Gibbs-free) phase_fourier kernel under delta-M,
 * replacing the unstable hard-cut + direct phase-value Fourier path. */
static int fixed_bulk_phase_table_chi_moments(const fixed_bulk_phase_table_t *tab, int L, double *chi_out)
{
    if (!tab || tab->n < 2 || L < 0 || !chi_out) return -1;
    for (int l = 0; l <= L; ++l) {
        double acc = 0.0;
        for (int i = 0; i < tab->n - 1; ++i) {
            double t0 = tab->theta_deg[i]   * M_PI / 180.0;
            double t1 = tab->theta_deg[i+1] * M_PI / 180.0;
            double f0 = tab->p11[i]   * fixed_bulk_legendre_P(l, cos(t0)) * sin(t0);
            double f1 = tab->p11[i+1] * fixed_bulk_legendre_P(l, cos(t1)) * sin(t1);
            acc += 0.5 * (f0 + f1) * (t1 - t0);
        }
        chi_out[l] = 0.5 * acc;
    }
    if (!(chi_out[0] > 0.0)) return -2;
    double inv = 1.0 / chi_out[0];
    for (int l = 0; l <= L; ++l) chi_out[l] *= inv;
    return 0;
}

/* "Formal delta-M" — HARD ANGULAR CUT sized by the Wiscombe fraction.
 *
 * Construction:
 *   f        = g_{l_trunc}                     (Wiscombe forward fraction)
 *              clamped to <= 1 - 2*bb/b - eps  (backscatter survival guard)
 *   cut_deg  = smallest tabulated theta whose cumulative forward mass
 *              reaches f * total
 *   P11(theta < cut) := 0                      (hard removal, NOT the
 *                                               Wiscombe delta-subtraction)
 *   renormalize to unit total; report residual bb fraction.
 *
 * MECHANISM CONCERN (why this operator is suspect): zeroing the forward
 * cone removes DIFFERENT angular mass than the Wiscombe delta + (1-f)
 * rescale it stands in for: the residual phase keeps its original
 * backscatter VALUES but the renormalization by the residual total
 * inflates the effective bb/b relative to the similarity-transformed
 * target (out_resid_bbfrac exposes exactly this).
 *
 * HISTORY OF THE "+48%" FIGURE — RETRACTED (S-004 closed 2026-07-05):
 * a "+48% at high-omega 412 nm (native path)" defect entry circulated in
 * the v1.08 handoff docs, but provenance audit found NO underlying
 * measurement anywhere (the citation chain was circular: validation
 * matrix -> harness doc §8 -> a bare one-line assertion), the recorded
 * code pointer ("rt_water_iop.c truncation logic") does not exist in the
 * source, and the native default path (OSOAA cap + value kernel, CAP FIX
 * 2026-06-17) does not even call this function.  User attested
 * (2026-07-05) the error was never observed.  The MECHANISM note above
 * remains valid as an analytic property of the hard cut; any future
 * defect claim against this operator must come with a reproducing
 * command.  DO NOT modify silently either way. */
static int fixed_bulk_phase_table_formal_delta_m(fixed_bulk_phase_table_t *tab,
                                                 int l_trunc,
                                                 double input_bb_over_b,
                                                 double *out_f,
                                                 double *out_cut_deg,
                                                 double *out_resid_bbfrac)
{
    if (!tab || tab->n < 2 || l_trunc < 1) return -1;
    if (!(input_bb_over_b > 0.0 && input_bb_over_b < 0.5)) return -2;
    double f = fixed_bulk_phase_table_legendre_moment(tab, l_trunc);
    if (!(f > 0.0)) return -3;
    const double max_f = 1.0 - 2.0*input_bb_over_b - 1.0e-10;
    if (f > max_f) f = max_f;
    if (f < 0.0) f = 0.0;
    if (!(f < 1.0)) return -4;

    const double total = fixed_bulk_phase_table_integral_theta(tab, 0.0, 180.0);
    const double remove_target = f * total;
    double cut_deg = 0.0;
    for (int i = 0; i < tab->n; ++i) {
        double c = fixed_bulk_phase_table_cumulative_forward(tab, tab->theta_deg[i]);
        if (c >= remove_target) { cut_deg = tab->theta_deg[i]; break; }
    }
    for (int i = 0; i < tab->n; ++i) if (tab->theta_deg[i] < cut_deg) tab->p11[i] = 0.0;
    double resid = fixed_bulk_phase_table_integral_theta(tab, 0.0, 180.0);
    if (!(resid > 0.0)) return -5;
    for (int i = 0; i < tab->n; ++i) tab->p11[i] /= resid;
    double back = fixed_bulk_phase_table_integral_theta(tab, 90.0, 180.0);
    if (out_f) *out_f = f;
    if (out_cut_deg) *out_cut_deg = cut_deg;
    if (out_resid_bbfrac) *out_resid_bbfrac = back;
    return 0;
}

/* Positivity-preserving angular forward-peak cap.  The method uses
 * two physical-angle anchors and a log-linear replacement of the unresolved lobe.  The extreme FF forward peak
 * (theta < T2) is REPLACED by the log-linear line fitted through the two
 * forward anchors (T1, P11(T1)) and (T2, P11(T2)); the backscatter side
 * (theta >= T2) keeps the analytic FF value UNCHANGED and therefore stays
 * positive.  The capped phase is renormalized to (1/2)INT P sin = 1 and the
 * removed forward fraction A is returned (delta-M f equivalent: the physical
 * scattering is preserved by b* = b(1-A), exactly as the existing
 * ccrr_particle_delta_f path does for Wiscombe delta-M).
 *
 * Why this is positivity-safe where Wiscombe delta-M is not: delta-M removes
 * only a forward delta (weight chi_L) and leaves a still-sharp residual whose
 * finite Legendre sum rings NEGATIVE at backscatter.  The angular cap instead
 * replaces the whole forward lobe with a smooth segment and never touches the
 * (analytic, positive) backscatter, so the phase VALUES are positive at every
 * angle.  Consumed by the value-based kernel below it gives positive radiance
 * at any nphi; if later expanded into Legendre moments it converges to a
 * positive function (no fundamental negativity). */
static int fixed_bulk_phase_table_loglinear_cap(fixed_bulk_phase_table_t *tab,
                                            double T1_deg, double T2_deg,
                                            double *out_A)
{
    if (!tab || tab->n < 2 || !(T1_deg > T2_deg) || !(T2_deg > 0.0)) return -1;
    const double total0 = fixed_bulk_phase_table_integral_theta(tab, 0.0, 180.0);
    if (!(total0 > 0.0)) return -2;
    const double T1r = T1_deg * M_PI / 180.0;
    const double T2r = T2_deg * M_PI / 180.0;
    const double P1 = fixed_bulk_phase_lut_eval(tab, cos(T1r));
    const double P2 = fixed_bulk_phase_lut_eval(tab, cos(T2r));
    if (!(P1 > 0.0) || !(P2 > 0.0)) return -3;
    const double AA = (log10(P2) - log10(P1)) / (T2r - T1r);  /* log10 slope per radian */
    for (int i = 0; i < tab->n; ++i) {
        if (tab->theta_deg[i] < T2_deg) {
            double thr = tab->theta_deg[i] * M_PI / 180.0;
            tab->p11[i] = pow(10.0, log10(P2) + AA * (thr - T2r));
        }
        /* theta >= T2: leave analytic backscatter value (positive) */
    }
    const double total1 = fixed_bulk_phase_table_integral_theta(tab, 0.0, 180.0);
    if (!(total1 > 0.0)) return -4;
    double A = 1.0 - total1 / total0;
    if (A < 0.0) A = 0.0;
    if (A > 1.0 - 1.0e-6) A = 1.0 - 1.0e-6;
    for (int i = 0; i < tab->n; ++i) tab->p11[i] /= total1;  /* renormalize total=1 */
    if (out_A) *out_A = A;
    return 0;
}


/* Angle-space VALUE kernel: m-th azimuthal Fourier component of the phase
 * function evaluated directly from P11 values (no Legendre expansion):
 *
 *   pfm[j][k] = (1/nphi) sum_q P11( cosTheta(mu_j, mu_k, phi_q) ) cos(m phi_q)
 *   cosTheta  = mu_j mu_k + sqrt(1-mu_j^2) sqrt(1-mu_k^2) cos(phi)
 *   phi_q     = 2 pi (q + 1/2)/nphi        (midpoint rule; exact for band-
 *                                           limited integrands up to nphi-1)
 *
 * EQUIVALENCE to the moment kernel (why the two paths agree at convergence):
 * spherical-harmonic addition theorem,
 *   P_l(cosTheta) = sum_m eps_m Ptil_l^m(mu_j) Ptil_l^m(mu_k) cos(m phi),
 *   eps_0 = 1, eps_{m>0} = 2,
 * gives (1/2pi) Int P11 cos(m phi) dphi
 *   = sum_l beta_l Ptil_l^m(mu_j) Ptil_l^m(mu_k)   for every m >= 0,
 * i.e. EXACTLY rt_kernel_phase_fourier() with untruncated moments.  The
 * value kernel therefore avoids Legendre-truncation ringing (Gibbs
 * negativity) at the cost of the phi quadrature; it is the production
 * in-water path and the template for the planned atmospheric aerosol value
 * kernel (HANDOFF_aerosol_value_kernel doc — the -20% nadir open defect is
 * the ABSENCE of this construction on the atm side). */

/* ============================================================================
 * v1.09 commit #8 (PERFORMANCE, bit-identical): all-m single-sweep builder for
 * the fixed-bulk direct phase-value Fourier kernel.
 *
 * WHY: fixed_bulk_direct_phase_fourier() was invoked once PER Fourier mode m,
 * but its inner scattering angle cth(j,k,q) — and therefore the phase LUT
 * interpolation P11(cth) — does not depend on m at all; only the cos(m*phi)
 * weight does.  Profiling (Tier-0 golden case, n_mu=48, nphi=720, 31 modes)
 * showed 1.24e8 calls to fixed_bulk_phase_lut_eval, ~37% of runtime, of which
 * 30/31 were exact repetitions.  This builder performs ONE (j,k,q) sweep and
 * accumulates every mode simultaneously.
 *
 * BIT-EXACTNESS ARGUMENT: for each mode m the accumulation
 *     acc_m += p11 * cos((double)m*phi)      (q ascending, then acc_m/nphi)
 * uses the same operands in the same order as the per-m original, so the
 * floating-point result is identical to the letter.  Verified: Tier-0 golden
 * rrs0minus bit-identical pre/post (2.887685e-02), full harness 18/18 bit.
 * The m=0 PHASENORM / NORMBB env branches of the original are NOT replicated
 * here; the caller uses the original function for m=0 whenever those env
 * toggles are active (both default-off diagnostics).
 * ========================================================================== */
static double fixed_bulk_value_eval(const fixed_bulk_phase_table_t *tab,
                                    const rt_value_phase_interp_t *spline,
                                    double mu)
{
    if (spline) {
        double p11 = 0.0;
        rt_value_phase_interp_eval(spline, mu, &p11, NULL, NULL);
        return p11;
    }
    return fixed_bulk_phase_lut_eval(tab, mu);
}

static int fixed_bulk_direct_phase_fourier_allm(int n_mu, const rt_atm_t *atm,
        int m_count, double g, const fixed_bulk_phase_table_t *tab,
        const rt_value_phase_interp_t *spline, int nphi,
        double *out /* [m][j][k] with k offset +n_mu, row stride (2n_mu+1) */)
{
    if(!atm||!out||m_count<1) return -1; if(nphi<16)nphi=16; if(nphi>20000)nphi=20000;
    const double two_pi=2.0*M_PI; const int kw=2*n_mu+1;
    double *acc=(double*)malloc((size_t)m_count*sizeof(double));
    if(!acc) return -1;
    for(int j=0;j<=n_mu;j++){ double mu_j=atm->rm[j]; double sj=sqrt(fmax(0.0,1.0-mu_j*mu_j));
      for(int k=-n_mu;k<=n_mu;k++){ double mu_k=atm->rm[k]; double sk=sqrt(fmax(0.0,1.0-mu_k*mu_k));
        for(int m=0;m<m_count;m++) acc[m]=0.0;
        for(int q=0;q<nphi;q++){ double phi=two_pi*((double)q+0.5)/(double)nphi;
          double cth=mu_j*mu_k+sj*sk*cos(phi); if(cth>1.0)cth=1.0; if(cth<-1.0)cth=-1.0;
          double p11 = tab ? fixed_bulk_value_eval(tab, spline, cth) : hg_phase_p11(cth, g);
          for(int m=0;m<m_count;m++) acc[m] += p11*cos((double)m*phi); }
        for(int m=0;m<m_count;m++)
          out[((size_t)m*(n_mu+1)+(size_t)j)*kw + (k+n_mu)] = acc[m]/(double)nphi;
      } }
    free(acc);
    return 0;
}

static int fixed_bulk_direct_phase_fourier(rt_legendre_workspace_t *ws,const rt_atm_t *atm,int m,double g,const fixed_bulk_phase_table_t *tab,const rt_value_phase_interp_t *spline,int nphi)
{
    if(!ws||!atm||m<0||m>ws->l_max) return -1; if(nphi<16)nphi=16; if(nphi>20000)nphi=20000; const int n_mu=ws->n_mu; const double two_pi=2.0*M_PI;
    for(int j=0;j<=n_mu;j++){ double mu_j=atm->rm[j]; double sj=sqrt(fmax(0.0,1.0-mu_j*mu_j)); for(int k=-n_mu;k<=n_mu;k++){ double mu_k=atm->rm[k]; double sk=sqrt(fmax(0.0,1.0-mu_k*mu_k)); double acc=0.0; for(int q=0;q<nphi;q++){ double phi=two_pi*((double)q+0.5)/(double)nphi; double cth=mu_j*mu_k+sj*sk*cos(phi); if(cth>1.0)cth=1.0; if(cth<-1.0)cth=-1.0; double p11 = tab ? fixed_bulk_value_eval(tab, spline, cth) : hg_phase_p11(cth, g); acc += p11*cos((double)m*phi); } ws->phase_fourier_m[j][k]=acc/(double)nphi; } }
    /* session6 OSOAA-style per-incidence phase normalization: OSOAA renormalizes
     * the hydrosol phase so that (1/2)integral P dmu' = 1 for every incidence on
     * its internal (dense) quadrature.  On the coarse SOS grid the discrete row
     * integral Sum_k w_k phase_fourier_m0[j][k] deviates from 1 (over-counts the
     * forward node, under-counts side/back), damping multiple scattering.  When
     * OCRT_FIXEDBULK_PHASENORM=1 we rescale each m=0 incidence row to unit norm,
     * which is the discrete analogue of OSOAA's normalization.  m>0 rows are left
     * untouched (they carry no net-scattering normalization). default: off. */
    /* session7 DIAGNOSTIC (default-off, production-invariant): per-incidence
     * effective scattering norm and effective backscatter fraction as seen by
     * the discrete SOS grid.  norm_eff(j) = Sum_k gb[k] phase_fourier_m0[j][k]
     * (target 2.0).  back(j) = Sum_{k: mu_k<0} gb[k] phase_fourier_m0[j][k]
     * = 2 * (effective bb/b for incidence j).  This is the quantity that the
     * principle "a,bb total fixes Rrs" hinges on: if the grid-integrated
     * effective bb/b differs from the post-integrated 0.4%-accurate bb/b, then
     * the RT-internal bb is NOT actually correct.  Activated by
     * OCRT_DUMP_NORMBB=1 (independent of PHASENORM). */
    if (m==0) { int dn = g_wrt_env.f_OCRT_DUMP_NORMBB;
      if(dn){ for(int j=0;j<=n_mu;j++){ double mu_j=atm->rm[j]; double sj=sqrt(fmax(0.0,1.0-mu_j*mu_j));
        double s=0.0,sb=0.0;
        /* recompute per-(j,k) integral splitting by TRUE scattering angle theta>90
         * (cos<0), not by node hemisphere mu_k<0, since for side incidence a
         * back-hemisphere node can still be a forward (small-theta) scatter. */
        const double two_pi2=2.0*M_PI; int nphi2=(nphi<256?256:nphi);
        for(int k=-n_mu;k<=n_mu;k++){ double mu_k=atm->rm[k]; double sk=sqrt(fmax(0.0,1.0-mu_k*mu_k));
          double accAll=0.0, accBack=0.0;
          for(int q=0;q<nphi2;q++){ double phi=two_pi2*((double)q+0.5)/(double)nphi2; double cth=mu_j*mu_k+sj*sk*cos(phi); if(cth>1.0)cth=1.0; if(cth<-1.0)cth=-1.0; double p=tab?fixed_bulk_value_eval(tab,spline,cth):hg_phase_p11(cth,g); accAll+=p; if(cth<0.0) accBack+=p; }
          accAll/=(double)nphi2; accBack/=(double)nphi2;
          s += atm->gb[k]*accAll; sb += atm->gb[k]*accBack; }
        fprintf(stderr,"NORMBB n_mu=%d j=%d mu=%.5f norm_eff=%.6f back2=%.6f bb_over_b_eff=%.6f\n",n_mu,j,mu_j,s,sb,0.5*sb); } } }
    if (m==0) { int pn = g_wrt_env.f_OCRT_FIXEDBULK_PHASENORM;
      if(pn){ int dbg = g_wrt_env.f_OCRT_DUMP_VIEW;
        for(int j=0;j<=n_mu;j++){ double s=0.0; for(int k=-n_mu;k<=n_mu;k++){ s += atm->gb[k]*ws->phase_fourier_m[j][k]; }
        if(dbg && (j==0||j==1||j==n_mu)) fprintf(stderr,"PHASENORM j=%d mu=%.4f norm_eff=%.5f\n",j,atm->rm[j],s);
        double ne = s; if(ne>1e-12){ double inv=2.0/ne; for(int k=-n_mu;k<=n_mu;k++) ws->phase_fourier_m[j][k]*=inv; } } } }
    return 0;
}

/* Recompute ONLY the solar-node-coupling entries of phase_fourier_m, i.e. row
 * j=0 and column k=0, where atm->rm[0] = -mu_sun_water is the (beam-dependent)
 * solar direction.  All other nodes (rm[+-1..+-n_mu]) are the fixed GL/view
 * quadrature, so the j,k>=1 block of the kernel is beam-INDEPENDENT and may be
 * cached across beams; only these solar entries change with the incident
 * beam's SZA.  Called after loading the cached GL block on a cache HIT so the
 * result is bit-equivalent to a full per-beam rebuild. */
static int fixed_bulk_direct_phase_fourier_solar(rt_legendre_workspace_t *ws,const rt_atm_t *atm,int m,const fixed_bulk_phase_table_t *tab,const rt_value_phase_interp_t *spline,int nphi)
{
    if(!ws||!atm||!tab||m<0||m>ws->l_max) return -1; if(nphi<16)nphi=16; if(nphi>20000)nphi=20000;
    const int n_mu=ws->n_mu; const double two_pi=2.0*M_PI;
    /* row j=0 (mu_j = rm[0] = -mu_sun) over all k */
    { const int j=0; double mu_j=atm->rm[0]; double sj=sqrt(fmax(0.0,1.0-mu_j*mu_j));
      for(int k=-n_mu;k<=n_mu;k++){ double mu_k=atm->rm[k]; double sk=sqrt(fmax(0.0,1.0-mu_k*mu_k)); double acc=0.0;
        for(int q=0;q<nphi;q++){ double phi=two_pi*((double)q+0.5)/(double)nphi; double cth=mu_j*mu_k+sj*sk*cos(phi); if(cth>1.0)cth=1.0; if(cth<-1.0)cth=-1.0; acc += fixed_bulk_value_eval(tab,spline,cth)*cos((double)m*phi); }
        ws->phase_fourier_m[j][k]=acc/(double)nphi; } }
    /* column k=0 (mu_k = rm[0] = -mu_sun) over all j */
    { const int k=0; double mu_k=atm->rm[0]; double sk=sqrt(fmax(0.0,1.0-mu_k*mu_k));
      for(int j=0;j<=n_mu;j++){ double mu_j=atm->rm[j]; double sj=sqrt(fmax(0.0,1.0-mu_j*mu_j)); double acc=0.0;
        for(int q=0;q<nphi;q++){ double phi=two_pi*((double)q+0.5)/(double)nphi; double cth=mu_j*mu_k+sj*sk*cos(phi); if(cth>1.0)cth=1.0; if(cth<-1.0)cth=-1.0; acc += fixed_bulk_value_eval(tab,spline,cth)*cos((double)m*phi); }
        ws->phase_fourier_m[j][k]=acc/(double)nphi; } }
    return 0;
}

/* ---------------------------------------------------------------------------
 * Snell helpers (local; avoid circular dep with rt_air_water)
 * --------------------------------------------------------------------------- */

static double snell_down(double mu_air, double n_water) {
    double sin_a_sq = fmax(0.0, 1.0 - mu_air * mu_air);
    double sin_w_sq = sin_a_sq / (n_water * n_water);
    if (sin_w_sq >= 1.0) return 0.0;
    return sqrt(1.0 - sin_w_sq);
}

/* ---------------------------------------------------------------------------
 * Single-scattering analytic Lu at view direction (in-water).
 * 사용 식: integrated single-scatter source over [0, τ_max]:
 *   L_SS(μ_v) = (ω · F_sun / (4π)) · P11(θ_scat)
 *              × (1/μ_v + 1/μ_s)^{-1} · (1 - exp(-τ_max·(1/μ_v + 1/μ_s)))
 * where P11 is the pure-water phase function (Rayleigh-like with δ_w=0.039),
 * scattering angle θ_scat between solar and view directions in water.
 * 이는 *infinite-deep + single scattering* 한계에서 정확한 analytic form.
 * --------------------------------------------------------------------------- */
double rt_water_rt_single_scatter_analytic(double mu_sun_w, double mu_view_w,
                                            double theta_scat_deg,
                                            double tau_max,
                                            double omega_w,
                                            double F_sun) {
    if (mu_sun_w <= 0.0 || mu_view_w <= 0.0) return 0.0;
    /* Pure-water Rayleigh-like P11 normalized so ∫P11 dΩ/4π = 1.
     * P11(θ) = (3/(2(2+δ))) [(1-δ)(1+cos²θ) + 2δ]  (Zhang 2009, δ=0.039) */
    const double delta = 0.039;
    double K = 3.0 / (2.0 * (2.0 + delta));
    double cos_th = cos(theta_scat_deg * M_PI / 180.0);
    double P11 = K * ((1.0 - delta) * (1.0 + cos_th * cos_th) + 2.0 * delta);

    double mu_inv_sum = 1.0 / mu_view_w + 1.0 / mu_sun_w;   /* = A */
    /* [v1.04+xsec_norm_cli+1, 2026-05-24] BUG FIX: path_attn 분모에 μ_v factor
     * 누락이었음. 표준 RTE solution at τ=0 (upwelling):
     *   L_u(0-, μ_v) = ∫₀^τ_max J(τ', μ_v) × exp(-τ'/μ_v) / μ_v dτ'
     * where J = (ω/4π) F_sun P11 × exp(-τ/μ_s) (single-scatter source).
     * → L_u(0-, μ_v) = (ω F_sun / 4π) P11 / μ_v × (1 - exp(-τ_max·A)) / A
     *                = (ω F_sun / 4π) P11 × (1 - exp(-τ_max·A)) / (μ_v × A)
     *                = (ω F_sun / 4π) P11 × (1 - exp(-τ_max·A)) / (1 + μ_v/μ_s)
     * Header docstring (rt_water_rt.h line 202-205) 의 식과 일치. 이전 .c 구현
     * (path_attn 분모가 단순 mu_inv_sum) 은 μ_v ↔ μ_s 대칭이라 비물리적이었음.
     * 본 함수는 production RT 미사용 (reference / verification only) 이므로
     * RT 결과 변경 없음 — diagnostic 정확성만 회복. */
    double path_attn  = (1.0 - exp(-tau_max * mu_inv_sum)) / (mu_view_w * mu_inv_sum);

    return (omega_w * F_sun / (4.0 * M_PI)) * P11 * path_attn;
}

/* ---------------------------------------------------------------------------
 * Helper: populate rt_atm_t for in-water Rayleigh-like medium.
 * Caller must rt_atm_alloc() before.
 *
 * Strategy: use the *aerosol slot* (xdel, betal_aer with L_max=2) to encode
 * water phase function with SSA built into xdel layer-mixing. ydel set to 0
 * (Rayleigh OFF) since water is not atmospheric Rayleigh.
 * --------------------------------------------------------------------------- */
static int build_inwater_atm(rt_atm_t *atm, int nt, int n_mu,
                              double tau_total, double omega_w,
                              double omega_particle,
                              double mu_sun_water,
                              int n_mu_quadrature,
                              int particle_L_max,
                              const double *particle_betal,
                              const double *particle_gammal,
                              const double *particle_alphal,
                              const double *particle_zetal)
{
    if (!atm || nt < 1 || n_mu < 1) return -1;
    if (n_mu_quadrature <= 0) n_mu_quadrature = n_mu;
    if (n_mu_quadrature > n_mu) return -1;
    if (!(mu_sun_water > 0.0 && mu_sun_water <= 1.0)) return -1;

    /* Common scalars */
    const double delta_w = 0.039;
    double ron_w = 2.0 * (1.0 - delta_w) / (2.0 + delta_w);

    atm->tau_total      = tau_total;
    atm->depol          = delta_w;
    atm->ssa            = omega_w + omega_particle;  /* metadata only; SOS doesn't read this */
    atm->mu_sun         = mu_sun_water;
    atm->rayleigh_model = RT_RAYLEIGH_MODEL_BODHAINE_1999;  /* unused; placeholder */

    /* Rayleigh slot ON — encode water phase function via atmospheric Rayleigh
     * fields. SSA ω_w is implicit through ydel scaling (layer-mixing trick).
     * 짚어둘 점 (B.3.2 디버깅 결과): aerosol-slot encoding이 |Q|² + |U|² > I²
     * unphysicality 유발. Rayleigh-slot encoding은 atmospheric baseline에서
     * 검증된 path이므로 numerically safe. */
    atm->beta0  = 1.0;
    atm->beta2  = 0.5 * ron_w;
    atm->gamma2 = -ron_w * sqrt(1.5);
    atm->alpha2 =  3.0 * ron_w;

    /* Layer arrays: uniform optical-depth discretization (or debug-gated
     * exponential stretch, OCRT_WATER_STRETCHED_GRID).
     * h[k]  = k * tau / nt  (k=0 -> just-below-surface, k=nt -> bottom)
     * ch[k] = exp(-h[k]/mu_sun)/2 (atm convention, see rt_atm.c header)
     * xdel[k] = omega_particle  (particle/hydrosol slot mixing fraction;
     *           0 when no particles — the old "aerosol slot OFF" comment
     *           predated the particle extension and contradicted the code)
     * ydel[k] = omega_w         (water Rayleigh slot; SOS source picks up
     *           the ssa through these fractions) */
    const double hr_dummy = 1.0;
    const char *stretch_env = ocrt_debug_env("OCRT_WATER_STRETCHED_GRID");
    const int use_stretched_grid = (stretch_env && stretch_env[0] && strcmp(stretch_env, "0") != 0);
    double stretch_s = 0.0;
    if (use_stretched_grid && tau_total > 0.0 && nt > 1) {
        const char *senv = ocrt_debug_env("OCRT_WATER_STRETCH_S");
        stretch_s = (senv && senv[0]) ? atof(senv) : 4.0;
        if (!(stretch_s > 0.0)) stretch_s = 4.0;
        if (stretch_s > 12.0) stretch_s = 12.0;
    }
    for (int k = 0; k <= nt; ++k) {
        double h_k;
        if (stretch_s > 0.0) {
            double x = (double)k / (double)nt;
            h_k = tau_total * expm1(stretch_s * x) / expm1(stretch_s);
        } else {
            h_k = (double)k * tau_total / (double)nt;
        }
        atm->h[k]    = h_k;
        atm->ch[k]   = 0.5 * exp(-h_k / mu_sun_water);
        atm->xdel[k] = omega_particle;
        atm->ydel[k] = omega_w;
        atm->z_km_level[k] = hr_dummy;
        if (k < nt) atm->tau_abs_layer[k] = 0.0;
    }
    atm->tau_abs_total = 0.0;

    /* Direction quadrature — mirror rt_atm_build_rayleigh logic.
     * Gauss-Legendre on (0,1] for n_mu positive nodes; mirror to negative.
     * Solar slot at j=0. */
    double mu_quad[256], w_quad[256];
    if (n_mu > 256 || n_mu_quadrature > 256) return -1;
    int rc = rt_quadrature_gauss_legendre_pos(n_mu_quadrature, mu_quad, w_quad);
    if (rc != 0) return -2;

    atm->rm[0] = -mu_sun_water;
    atm->gb[0] = 0.0;
    for (int j = 1; j <= n_mu; ++j) {
        atm->rm[+j] = 0.0;
        atm->rm[-j] = 0.0;
        atm->gb[+j] = 0.0;
        atm->gb[-j] = 0.0;
    }
    for (int i = 0; i < n_mu_quadrature; ++i) {
        int j_pos = i + 1;
        atm->rm[+j_pos] = mu_quad[i];
        atm->rm[-j_pos] = -mu_quad[i];
        atm->gb[+j_pos] = w_quad[i];
        atm->gb[-j_pos] = w_quad[i];
    }

    /* CCRR/particle slot.  Scalar P11 moments (betal) always; vector Mueller
     * couplings (gammal=P12, alphal/zetal=P22/P33) when provided (#0, the
     * .mie vector phase path).  When the polarized arrays are NULL the slots
     * stay zero (scalar particle, backward-compatible — FF scalar LUT etc.).
     * Pure-water Rayleigh polarization remains in the Rayleigh slot (ydel). */
    if (atm->betal_aer)  { free(atm->betal_aer);  atm->betal_aer  = NULL; }
    if (atm->gammal_aer) { free(atm->gammal_aer); atm->gammal_aer = NULL; }
    if (atm->alphal_aer) { free(atm->alphal_aer); atm->alphal_aer = NULL; }
    if (atm->zetal_aer)  { free(atm->zetal_aer);  atm->zetal_aer  = NULL; }
    atm->aerosol_active = 0;
    atm->L_max          = 0;
    if (omega_particle > 0.0 && particle_L_max >= 0 && particle_betal) {
        int L = particle_L_max;
        size_t bytes = (size_t)(L + 1) * sizeof(double);
        atm->betal_aer  = (double*)calloc((size_t)(L + 1), sizeof(double));
        atm->gammal_aer = (double*)calloc((size_t)(L + 1), sizeof(double));
        atm->alphal_aer = (double*)calloc((size_t)(L + 1), sizeof(double));
        atm->zetal_aer  = (double*)calloc((size_t)(L + 1), sizeof(double));
        if (!atm->betal_aer || !atm->gammal_aer || !atm->alphal_aer || !atm->zetal_aer) return -3;
        memcpy(atm->betal_aer, particle_betal, bytes);
        if (particle_gammal) memcpy(atm->gammal_aer, particle_gammal, bytes);
        if (particle_alphal) memcpy(atm->alphal_aer, particle_alphal, bytes);
        if (particle_zetal)  memcpy(atm->zetal_aer,  particle_zetal,  bytes);
        atm->aerosol_active = 1;
        atm->L_max = L;
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * Determine τ_max for infinite-deep BC.
 * Heuristic: τ_max = K_factor / min(a_w + b_w) over spectrum.
 * K_factor = 5 gives ~99.3% absorption depth (5 e-foldings).
 *
 * For pure water 350-900 nm: min(a_w + b_w) occurs in blue (~450 nm) where
 * a_w ~ 0.0145 m⁻¹, b_w ~ 0.00405 → ext ~ 0.0186 → z_max(τ=5) ~ 269 m.
 * --------------------------------------------------------------------------- */
static double determine_tau_max(double a_w, double b_w, double tau_max_target) {
    if (tau_max_target > 0.0) return tau_max_target;
    /* Production base depth: at least total optical depth 20.  A physical-depth
     * floor and matching vertical-grid policy are applied after the active
     * extinction is known so high-TSM water is not truncated to a shallow slab. */
    (void)a_w; (void)b_w;
    return 20.0;
}

/* ---------------------------------------------------------------------------
 * Main entry: in-water SOS for pure water case
 * --------------------------------------------------------------------------- */
/* ===========================================================================
 * Wavelength-keyed value-kernel cache (path B speed optimization)
 *
 * The value-based azimuth Fourier kernel (fixed_bulk_direct_phase_fourier) is
 * IDENTICAL for every rt_water_rt_sos_pure call at a given wavelength: the
 * direct-beam call and all ~N_mu atm-sky-light beam calls share the same
 * in-water phase function (phase depends on wavelength + constituents, not on
 * the incident beam direction).  Without caching it is rebuilt per call
 * (nphi x (n_mu+1) x (2n_mu+1) x (m_max+1) LUT evals) ~ the dominant cost.
 * This single-entry cache stores the per-m kernel keyed by the full phase
 * signature, so the first beam computes it and every subsequent beam at the
 * same wavelength reuses it.  Memory ~ (m_max+1)*(n_mu+1)*(2n_mu+1) doubles
 * (a few hundred KB) — trading memory for a large speedup, as intended.
 * The water beam loop in rt_solver.c is serial (no OpenMP), so a module-static
 * cache is race-free.  (TODO: multi-entry LUT cache for --batch multi-λ runs;
 * single entry already removes the ~N_mu within-wavelength redundancy.)
 * ======================================================================== */
typedef struct {
    int    valid;
    double lambda_nm, chl, min_c, adom440, adomS, n_water, T_water, nphi_d, mu_view;
    int    n_mu, m_count, method, lut_mode;
    char   lut_path[1024];
    double **kernel_m;   /* [m_count]; each flat (n_mu+1)*(2n_mu+1), idx [j*dirs+(k+n_mu)] */
} water_phase_kernel_cache_t;
static water_phase_kernel_cache_t g_wpkc = {0};

/* Dedicated polarized direct-value kernel cache.
 *
 * Each worker retains multiple complete all-direction kernels so a spectral
 * batch can keep every requested wavelength/active-mixture kernel in RAM.
 * Entries are immutable after construction.  Beam-dependent row/column strips
 * are cached per entry; one worker-private scratch buffer materializes a
 * non-base beam without duplicating a full work buffer for every wavelength.
 */
#define WATER_POL_VALUE_BEAM_CACHE_CAP 128
#define WATER_POL_VALUE_KERNEL_CACHE_CAP 64

typedef struct {
    int valid;
    unsigned long long mu_bits;
    double *strip; /* six x m_count x (row[dirs] + col[n_mu]) */
} water_pol_value_beam_entry_t;

typedef struct {
    int valid;
    char path[512];
    double lambda_nm;
    int trunc, ss, moment_mode, nmg, lmax_req;
    int n_mu, m_count, nphi;
    unsigned long long grid_hash;
    unsigned long long stamp;
    unsigned long long base_mu_bits;
    size_t plane, all;
    size_t strip_plane, strip_all;
    double *base_buf; /* six contiguous full [m][j][k] tables */
    water_pol_value_beam_entry_t beam[WATER_POL_VALUE_BEAM_CACHE_CAP];
    unsigned long long full_build_count;
    unsigned long long phase_hit_count;
    unsigned long long beam_build_count;
    unsigned long long beam_hit_count;
} water_pol_value_kernel_cache_t;

static water_pol_value_kernel_cache_t
    g_wpvkc[WATER_POL_VALUE_KERNEL_CACHE_CAP];
static unsigned long long g_wpvkc_tick = 0;
static double *g_wpvkc_work_buf = NULL;
static size_t g_wpvkc_work_cap = 0; /* doubles */
#ifdef _OPENMP
#pragma omp threadprivate(g_wpvkc, g_wpvkc_tick, g_wpvkc_work_buf, g_wpvkc_work_cap)
#endif

static unsigned long long wpvkc_grid_hash(const rt_atm_t *atm)
{
    /* Hash only the beam-independent direction ring.  atm->rm is an offset
     * pointer indexable as [-n_mu..+n_mu]; hashing from rm[0] as a flat block
     * reads past rm_storage and also folds the solar slot into the key.  The
     * latter defeats the dedicated beam-strip cache. */
    unsigned long long h=1469598103934665603ULL;
    const int n_mu=atm->n_mu;
    const unsigned char *pn=(const unsigned char*)&n_mu;
    for(size_t i=0;i<sizeof n_mu;++i){ h^=(unsigned long long)pn[i]; h*=1099511628211ULL; }
    for(int j=-n_mu;j<=n_mu;++j){
        if(j==0) continue;
        const double mu=atm->rm[j];
        const unsigned char *p=(const unsigned char*)&mu;
        for(size_t i=0;i<sizeof mu;++i){ h^=(unsigned long long)p[i]; h*=1099511628211ULL; }
    }
    return h;
}

static unsigned long long wpvkc_double_bits(double x)
{
    unsigned long long u = 0ULL;
    memcpy(&u, &x, sizeof u);
    return u;
}

static size_t wpvkc_idx(int m, int j, int k, int n_mu)
{
    return ((size_t)m * (size_t)(n_mu + 1) + (size_t)j) *
           (size_t)(2 * n_mu + 1) + (size_t)(k + n_mu);
}

static void wpvkc_entry_free(water_pol_value_kernel_cache_t *c)
{
    if (!c) return;
    for (int i=0; i<WATER_POL_VALUE_BEAM_CACHE_CAP; ++i)
        free(c->beam[i].strip);
    free(c->base_buf);
    memset(c, 0, sizeof *c);
}

static int wpvkc_key_match(const water_pol_value_kernel_cache_t *c,
                           const char *path, double lambda_nm,
                           int trunc, int ss, int moment_mode, int nmg,
                           int lmax_req, const rt_atm_t *atm,
                           int m_count, int nphi)
{
    return c && c->valid && c->lambda_nm==lambda_nm &&
           c->trunc==trunc && c->ss==ss &&
           c->moment_mode==moment_mode && c->nmg==nmg &&
           c->lmax_req==lmax_req && c->n_mu==atm->n_mu &&
           c->m_count==m_count && c->nphi==nphi &&
           c->grid_hash==wpvkc_grid_hash(atm) && c->base_buf &&
           strcmp(c->path,path?path:"")==0;
}

static int wpvkc_entry_reset(water_pol_value_kernel_cache_t *c,
                             const char *path, double lambda_nm,
                             int trunc, int ss, int moment_mode, int nmg,
                             int lmax_req, const rt_atm_t *atm,
                             int m_count, int nphi)
{
    if (!c) return -1;
    wpvkc_entry_free(c);
    const size_t plane=(size_t)(atm->n_mu+1)*(size_t)(2*atm->n_mu+1);
    const size_t all=(size_t)m_count*plane;
    const size_t strip_plane=(size_t)m_count*
        (size_t)((2*atm->n_mu+1)+atm->n_mu);
    c->base_buf=(double*)calloc(6u*all,sizeof(double));
    if(!c->base_buf) { wpvkc_entry_free(c); return -1; }
    if (snprintf(c->path,sizeof c->path,"%s",path?path:"") >=
        (int)sizeof c->path) {
        wpvkc_entry_free(c); return -1;
    }
    c->lambda_nm=lambda_nm; c->trunc=trunc; c->ss=ss;
    c->moment_mode=moment_mode; c->nmg=nmg; c->lmax_req=lmax_req;
    c->n_mu=atm->n_mu; c->m_count=m_count; c->nphi=nphi;
    c->grid_hash=wpvkc_grid_hash(atm); c->plane=plane; c->all=all;
    c->strip_plane=strip_plane; c->strip_all=6u*strip_plane;
    c->stamp=++g_wpvkc_tick;
    c->valid=1;
    return 0;
}

static water_pol_value_kernel_cache_t *wpvkc_find(
        const char *path, double lambda_nm,
        int trunc, int ss, int moment_mode, int nmg, int lmax_req,
        const rt_atm_t *atm, int m_count, int nphi)
{
    for (int i=0; i<WATER_POL_VALUE_KERNEL_CACHE_CAP; ++i) {
        water_pol_value_kernel_cache_t *c=&g_wpvkc[i];
        if (!wpvkc_key_match(c,path,lambda_nm,trunc,ss,moment_mode,nmg,
                             lmax_req,atm,m_count,nphi)) continue;
        c->stamp=++g_wpvkc_tick;
        return c;
    }
    return NULL;
}

static water_pol_value_kernel_cache_t *wpvkc_acquire(
        const char *path, double lambda_nm,
        int trunc, int ss, int moment_mode, int nmg, int lmax_req,
        const rt_atm_t *atm, int m_count, int nphi)
{
    int victim=0;
    for (int i=0; i<WATER_POL_VALUE_KERNEL_CACHE_CAP; ++i) {
        if (!g_wpvkc[i].valid) { victim=i; break; }
        if (g_wpvkc[i].stamp < g_wpvkc[victim].stamp) victim=i;
    }
    if (wpvkc_entry_reset(&g_wpvkc[victim],path,lambda_nm,trunc,ss,
                           moment_mode,nmg,lmax_req,atm,m_count,nphi)!=0)
        return NULL;
    return &g_wpvkc[victim];
}

static double *wpvkc_work_buffer(size_t doubles)
{
    if (g_wpvkc_work_cap >= doubles && g_wpvkc_work_buf)
        return g_wpvkc_work_buf;
    double *q=(double*)realloc(g_wpvkc_work_buf,doubles*sizeof(double));
    if (!q) return NULL;
    g_wpvkc_work_buf=q;
    g_wpvkc_work_cap=doubles;
    return q;
}

static void wpvkc_capture_strip(const water_pol_value_kernel_cache_t *c,
                                const double *full, double *strip)
{
    const int n_mu=c->n_mu;
    const int dirs=2*n_mu+1;
    const int m_count=c->m_count;
    const size_t all=c->all;
    const size_t stride=(size_t)(dirs+n_mu);
    for (int comp=0; comp<6; ++comp) {
        const double *src=full+(size_t)comp*all;
        double *dst=strip+(size_t)comp*c->strip_plane;
        for (int m=0; m<m_count; ++m) {
            double *d=dst+(size_t)m*stride;
            for (int k=-n_mu; k<=n_mu; ++k)
                *d++=src[wpvkc_idx(m,0,k,n_mu)];
            for (int j=1; j<=n_mu; ++j)
                *d++=src[wpvkc_idx(m,j,0,n_mu)];
        }
    }
}

static void wpvkc_apply_strip(const water_pol_value_kernel_cache_t *c,
                              double *full, const double *strip)
{
    const int n_mu=c->n_mu;
    const int dirs=2*n_mu+1;
    const int m_count=c->m_count;
    const size_t all=c->all;
    const size_t stride=(size_t)(dirs+n_mu);
    for (int comp=0; comp<6; ++comp) {
        double *dst=full+(size_t)comp*all;
        const double *src=strip+(size_t)comp*c->strip_plane;
        for (int m=0; m<m_count; ++m) {
            const double *q=src+(size_t)m*stride;
            for (int k=-n_mu; k<=n_mu; ++k)
                dst[wpvkc_idx(m,0,k,n_mu)]=*q++;
            for (int j=1; j<=n_mu; ++j)
                dst[wpvkc_idx(m,j,0,n_mu)]=*q++;
        }
    }
}

static water_pol_value_beam_entry_t *wpvkc_find_beam(
        water_pol_value_kernel_cache_t *c, double mu0)
{
    const unsigned long long bits=wpvkc_double_bits(mu0);
    for (int i=0; i<WATER_POL_VALUE_BEAM_CACHE_CAP; ++i) {
        water_pol_value_beam_entry_t *e=&c->beam[i];
        if (e->valid && e->mu_bits==bits && e->strip) return e;
    }
    return NULL;
}

static int wpvkc_store_beam(water_pol_value_kernel_cache_t *c,
                            double mu0, const double *full)
{
    water_pol_value_beam_entry_t *existing=wpvkc_find_beam(c,mu0);
    if (existing) {
        wpvkc_capture_strip(c,full,existing->strip);
        return 0;
    }
    for (int i=0; i<WATER_POL_VALUE_BEAM_CACHE_CAP; ++i) {
        water_pol_value_beam_entry_t *e=&c->beam[i];
        if (e->valid) continue;
        e->strip=(double*)malloc(c->strip_all*sizeof(double));
        if (!e->strip) return -1;
        e->valid=1;
        e->mu_bits=wpvkc_double_bits(mu0);
        wpvkc_capture_strip(c,full,e->strip);
        return 0;
    }
    return 1; /* correctness unaffected; this beam is simply not retained */
}

static size_t wpvkc_entry_resident_bytes(
        const water_pol_value_kernel_cache_t *c)
{
    if (!c || !c->valid) return 0;
    size_t n=6u*c->all*sizeof(double);
    for (int i=0; i<WATER_POL_VALUE_BEAM_CACHE_CAP; ++i)
        if (c->beam[i].valid && c->beam[i].strip)
            n += c->strip_all*sizeof(double);
    return n;
}

static size_t wpvkc_resident_bytes_all(void)
{
    size_t n=g_wpvkc_work_cap*sizeof(double);
    for (int i=0; i<WATER_POL_VALUE_KERNEL_CACHE_CAP; ++i)
        n += wpvkc_entry_resident_bytes(&g_wpvkc[i]);
    return n;
}

static int wpvkc_entry_count(void)
{
    int n=0;
    for (int i=0; i<WATER_POL_VALUE_KERNEL_CACHE_CAP; ++i)
        if (g_wpvkc[i].valid) ++n;
    return n;
}


/* ---------------------------------------------------------------------------
 * v1.09 commit #24 (full-grid perf): per-(wavelength, .mie path) cache of the
 * water-particle phase decomposition.  read_mie_file() + build_aerosol_
 * interpolators() + rt_aerosol_compute_vector_legendre*() were being re-run on
 * EVERY grid cell (every (vza,raa)), each doing a disk read of the .mie file
 * and an L<=200 Legendre expansion.  Measured cost: ~2.3 s/cell at L200 (half
 * of that the expansion), which both dominated runtime and serialized batch
 * threads on the disk read.  These outputs depend ONLY on (lambda, .mie file,
 * truncation/ss/moment settings), so they are identical across all cells of a
 * grid.  Cache them; the grid loop's per-cell SOS (which uses the cached
 * moments) is the only thing that must actually re-run.
 * threadprivate: each batch thread owns its cache (same pattern as g_wpkc). */
typedef struct {
    int    valid;
    char   path[512];
    double lambda_nm;
    int    trunc, ss, moment_mode, nmg, lmax_req;
    /* cached mie table + interpolator (owned here; freed on key change) */
    mie_data_t              mie;
    aerosol_phase_interp_t  aip;
    int    have_mie, have_aip;
    /* cached decomposition outputs */
    double betal[FIXED_BULK_BETAL_LMAX + 1];
    double gammal[FIXED_BULK_BETAL_LMAX + 1];
    double alphal[FIXED_BULK_BETAL_LMAX + 1];
    double zetal[FIXED_BULK_BETAL_LMAX + 1];
    int    particle_Lmax;
    double delta_f;      /* forward-peak fraction from truncation (0 if none) */
    double A_tr;         /* OSOAA A_TRONCA (2f) if truncation, else 0 */
} water_mie_cache_t;
static water_mie_cache_t g_wmc = {0};
#pragma omp threadprivate(g_wmc)

static int wmc_key_match(const char *path, double lambda_nm, int trunc, int ss,
                         int moment_mode, int nmg, int lmax_req) {
    return g_wmc.valid && path && g_wmc.path[0] &&
           strcmp(g_wmc.path, path) == 0 &&
           g_wmc.lambda_nm == lambda_nm && g_wmc.trunc == trunc &&
           g_wmc.ss == ss && g_wmc.moment_mode == moment_mode &&
           g_wmc.nmg == nmg && g_wmc.lmax_req == lmax_req;
}
static void wmc_free(void) {
    if (g_wmc.have_aip) { aerosol_phase_interp_free(&g_wmc.aip); g_wmc.have_aip = 0; }
    if (g_wmc.have_mie) { mie_data_free(&g_wmc.mie); g_wmc.have_mie = 0; }
    g_wmc.valid = 0;
}

/* Parsed .mie models are cached worker-locally by shared/mie_io.c.  The full
 * model (all bulk and phase wavelengths) remains resident and is reused by
 * both atmospheric-aerosol and in-water particle paths. */

/* ---------------------------------------------------------------------------
 * Stage 3B direct theta-linear particle-phase caches.
 *
 * The production particle path no longer reconstructs a strongly forward-
 * peaked Mie phase from finite generalized moments.  A component table is
 * wavelength-linearly interpolated with one common bracket at its native angle nodes, then represented by
 * rt_value_phase_interp_t (FR631 O(1) lookup or generic theta-linear fallback).
 * The cache owns every table and is thread-private, matching the existing
 * per-worker moment caches.
 * --------------------------------------------------------------------------- */
#define WATER_COMPONENT_VALUE_CACHE_CAP 24
#define WATER_HYD_TRUNC_MU1       0.85
#define WATER_HYD_TRUNC_MU2       0.92
#define WATER_HYD_TRUNC_THRESHOLD 0.10

typedef struct {
    int valid;
    char path[1024];
    double lambda_nm;
    int truncation_mode;
    unsigned long long stamp;
    rt_value_phase_interp_t phase;
    double ssa;
    double g_bulk;
    double A_trunc;
} water_component_value_cache_entry_t;

static water_component_value_cache_entry_t
    g_wcvc[WATER_COMPONENT_VALUE_CACHE_CAP];
static unsigned long long g_wcvc_tick = 0;
static unsigned long long g_wcvc_build_count = 0;
static unsigned long long g_wcvc_hit_count = 0;
#ifdef _OPENMP
#pragma omp threadprivate(g_wcvc, g_wcvc_tick, g_wcvc_build_count, g_wcvc_hit_count)
#endif

static void water_component_value_cache_entry_clear(
        water_component_value_cache_entry_t *e)
{
    if (!e) return;
    if (e->valid || e->phase.n > 0) rt_value_phase_interp_free(&e->phase);
    memset(e, 0, sizeof *e);
}

static int water_component_value_phase_get(
        const char *path, double lambda_nm, int truncation_mode,
        const rt_value_phase_interp_t **phase_out,
        double *ssa_out, double *g_bulk_out, double *A_trunc_out,
        int *cache_hit)
{
    if (phase_out) *phase_out = NULL;
    if (A_trunc_out) *A_trunc_out = 0.0;
    if (cache_hit) *cache_hit = 0;
    if (!path || !path[0] || !phase_out || !isfinite(lambda_nm) ||
        (truncation_mode != 0 && truncation_mode != 1)) return -1;

    const char *disable_cache_env =
        getenv("OCRT_DISABLE_WATER_COMPONENT_PHASE_CACHE");
    const int disable_cache = disable_cache_env && disable_cache_env[0] &&
                              strcmp(disable_cache_env, "0") != 0;

    for (int i = 0; !disable_cache && i < WATER_COMPONENT_VALUE_CACHE_CAP; ++i) {
        water_component_value_cache_entry_t *e = &g_wcvc[i];
        if (!e->valid || e->lambda_nm != lambda_nm ||
            e->truncation_mode != truncation_mode ||
            strcmp(e->path, path) != 0) continue;
        e->stamp = ++g_wcvc_tick;
        ++g_wcvc_hit_count;
        *phase_out = &e->phase;
        if (ssa_out) *ssa_out = e->ssa;
        if (g_bulk_out) *g_bulk_out = e->g_bulk;
        if (A_trunc_out) *A_trunc_out = e->A_trunc;
        if (cache_hit) *cache_hit = 1;
        const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
        if (dump && dump[0] && strcmp(dump, "0") != 0)
            fprintf(stderr,
                    "WCVALUE hit path=%s lambda=%.10g trunc=%d A=%.12g "
                    "n=%d fr631=%d norm0=%.12g g=%.12g bb_b=%.12g "
                    "builds=%llu hits=%llu\n",
                    path, lambda_nm, truncation_mode, e->A_trunc,
                    e->phase.n, e->phase.is_fr631,
                    e->phase.norm_before, e->phase.g_asym,
                    e->phase.bb_b_ratio,
                    g_wcvc_build_count, g_wcvc_hit_count);
        return 0;
    }

    const mie_data_t *mie = NULL;
    mie_data_t mie_tmp = {0};
    int model_cache_hit = 0;
    const int mrc = mie_model_cache_get(path, &mie, &model_cache_hit);
    (void)model_cache_hit;
    if (mrc == 1) {
        if (read_mie_file(path, &mie_tmp) != 0) return -2;
        mie = &mie_tmp;
    } else if (mrc != 0 || !mie) {
        return -2;
    }
    const int n = mie->n_ang;
    double *p11 = (double *)malloc((size_t)n * sizeof(double));
    double *p12 = (double *)malloc((size_t)n * sizeof(double));
    double *p33 = (double *)malloc((size_t)n * sizeof(double));
    if (!p11 || !p12 || !p33) {
        free(p11); free(p12); free(p33);
        if (mrc == 1) mie_data_free(&mie_tmp);
        return -3;
    }
    if (mie_phase_nodes_at_wavelength_linear(
            mie, lambda_nm * 1.0e-3, p11, p12, p33) != 0) {
        free(p11); free(p12); free(p33);
        if (mrc == 1) mie_data_free(&mie_tmp);
        return -4;
    }

    rt_value_phase_interp_t built_raw = {0};
    const int brc = rt_value_phase_interp_build(
        &built_raw, mie->angles, p11, p12, p33, n, 0);
    free(p11); free(p12); free(p33);
    if (brc != 0) {
        if (mrc == 1) mie_data_free(&mie_tmp);
        return -5;
    }

    double ssa = 0.0, g_bulk = 0.0;
    if (mie_bulk_diagnostics_at_wavelength(
            mie, lambda_nm * 1.0e-3, &ssa, &g_bulk) != 0) {
        rt_value_phase_interp_free(&built_raw);
        if (mrc == 1) mie_data_free(&mie_tmp);
        return -6;
    }
    if (mrc == 1) mie_data_free(&mie_tmp);

    rt_value_phase_interp_t built = {0};
    double A_trunc = 0.0;
    if (truncation_mode == 1) {
        const int trc = rt_value_phase_interp_loglinear_truncate(
            &built_raw, WATER_HYD_TRUNC_MU1, WATER_HYD_TRUNC_MU2,
            WATER_HYD_TRUNC_THRESHOLD, &built, &A_trunc);
        rt_value_phase_interp_free(&built_raw);
        if (trc != 0) return -8;
    } else {
        built = built_raw;
        memset(&built_raw, 0, sizeof built_raw);
    }

    if (disable_cache) {
        /* The caller needs stable ownership for the whole solve.  Even under
         * the diagnostic no-cache flag, retain one LRU slot rather than
         * returning a pointer to a temporary allocation. */
    }
    int victim = 0;
    for (int i = 0; i < WATER_COMPONENT_VALUE_CACHE_CAP; ++i) {
        if (!g_wcvc[i].valid) { victim = i; break; }
        if (g_wcvc[i].stamp < g_wcvc[victim].stamp) victim = i;
    }
    water_component_value_cache_entry_t *e = &g_wcvc[victim];
    water_component_value_cache_entry_clear(e);
    if (snprintf(e->path, sizeof e->path, "%s", path) >= (int)sizeof e->path) {
        rt_value_phase_interp_free(&built);
        return -7;
    }
    e->valid = 1;
    e->lambda_nm = lambda_nm;
    e->truncation_mode = truncation_mode;
    e->stamp = ++g_wcvc_tick;
    e->phase = built;
    e->ssa = ssa;
    e->g_bulk = g_bulk;
    e->A_trunc = A_trunc;
    ++g_wcvc_build_count;

    *phase_out = &e->phase;
    if (ssa_out) *ssa_out = e->ssa;
    if (g_bulk_out) *g_bulk_out = e->g_bulk;
    if (A_trunc_out) *A_trunc_out = e->A_trunc;
    const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
    if (dump && dump[0] && strcmp(dump, "0") != 0)
        fprintf(stderr,
                "WCVALUE build path=%s lambda=%.10g trunc=%d A=%.12g "
                "n=%d fr631=%d norm0=%.12g g=%.12g bb_b=%.12g "
                "builds=%llu hits=%llu\n",
                path, lambda_nm, truncation_mode, e->A_trunc,
                e->phase.n, e->phase.is_fr631,
                e->phase.norm_before, e->phase.g_asym,
                e->phase.bb_b_ratio,
                g_wcvc_build_count, g_wcvc_hit_count);
    return 0;
}

#define WATER_MIXED_VALUE_CACHE_CAP 24

typedef struct {
    int valid;
    char path[3][1024];
    double lambda_nm;
    double b[3];
    int truncation_mode;
    int n_component;
    unsigned long long stamp;
    unsigned long long signature;
    rt_value_phase_interp_t phase;
} water_mixed_value_cache_entry_t;

static water_mixed_value_cache_entry_t
    g_wmvc[WATER_MIXED_VALUE_CACHE_CAP];
static unsigned long long g_wmvc_tick = 0;
static unsigned long long g_wmvc_build_count = 0;
static unsigned long long g_wmvc_hit_count = 0;
#ifdef _OPENMP
#pragma omp threadprivate(g_wmvc, g_wmvc_tick, g_wmvc_build_count, g_wmvc_hit_count)
#endif

static void water_mixed_value_cache_entry_clear(water_mixed_value_cache_entry_t *e)
{
    if (!e) return;
    if (e->valid || e->phase.n > 0) rt_value_phase_interp_free(&e->phase);
    memset(e, 0, sizeof *e);
}

static unsigned long long fnv1a_bytes(unsigned long long h,
                                      const void *data, size_t n)
{
    const unsigned char *p = (const unsigned char *)data;
    for (size_t i = 0; i < n; ++i) {
        h ^= (unsigned long long)p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

static unsigned long long water_mixed_value_signature(
        const char *const path[3], const double b[3], double lambda_nm,
        int truncation_mode)
{
    unsigned long long h = 1469598103934665603ULL;
    h = fnv1a_bytes(h, &lambda_nm, sizeof lambda_nm);
    h = fnv1a_bytes(h, &truncation_mode, sizeof truncation_mode);
    for (int i = 0; i < 3; ++i) {
        const char *p = path[i] ? path[i] : "";
        h = fnv1a_bytes(h, p, strlen(p) + 1u);
        h = fnv1a_bytes(h, &b[i], sizeof b[i]);
    }
    return h;
}

static int water_mixed_value_key_match(
        const water_mixed_value_cache_entry_t *e,
        const char *const path[3], const double b[3], double lambda_nm,
        int truncation_mode)
{
    if (!e || !e->valid || e->lambda_nm != lambda_nm ||
        e->truncation_mode != truncation_mode) return 0;
    for (int i = 0; i < 3; ++i) {
        const char *p = path[i] ? path[i] : "";
        if (e->b[i] != b[i] || strcmp(e->path[i], p) != 0) return 0;
    }
    return 1;
}

static int water_mixed_value_phase_get(
        const char *const path[3], const double b[3], double lambda_nm,
        int truncation_mode,
        const rt_value_phase_interp_t **phase_out,
        int *n_component_out, unsigned long long *signature_out,
        int *cache_hit)
{
    if (phase_out) *phase_out = NULL;
    if (n_component_out) *n_component_out = 0;
    if (signature_out) *signature_out = 0ULL;
    if (cache_hit) *cache_hit = 0;
    if (!path || !b || !phase_out || !isfinite(lambda_nm) ||
        (truncation_mode != 0 && truncation_mode != 1)) return -1;

    for (int i = 0; i < WATER_MIXED_VALUE_CACHE_CAP; ++i) {
        water_mixed_value_cache_entry_t *e = &g_wmvc[i];
        if (!water_mixed_value_key_match(
                e, path, b, lambda_nm, truncation_mode)) continue;
        e->stamp = ++g_wmvc_tick;
        ++g_wmvc_hit_count;
        *phase_out = &e->phase;
        if (n_component_out) *n_component_out = e->n_component;
        if (signature_out) *signature_out = e->signature;
        if (cache_hit) *cache_hit = 1;
        return 0;
    }

    double bsum = 0.0;
    int n_component = 0;
    for (int i = 0; i < 3; ++i) {
        if (b[i] > 0.0) {
            if (!path[i] || !path[i][0]) return -2;
            bsum += b[i];
            ++n_component;
        }
    }
    if (!(bsum > 0.0) || n_component < 1) return -3;

    double theta[RT_FR631_N_ANGLE];
    double p11[RT_FR631_N_ANGLE];
    double p12[RT_FR631_N_ANGLE];
    double p33[RT_FR631_N_ANGLE];
    rt_fr631_fill_theta(theta);
    for (int k = 0; k < RT_FR631_N_ANGLE; ++k) {
        p11[k] = 0.0; p12[k] = 0.0; p33[k] = 0.0;
    }

    for (int i = 0; i < 3; ++i) {
        if (!(b[i] > 0.0)) continue;
        const rt_value_phase_interp_t *component = NULL;
        const int rc = water_component_value_phase_get(
            path[i], lambda_nm, truncation_mode,
            &component, NULL, NULL, NULL, NULL);
        if (rc != 0 || !component) return -5;
        const double w = b[i] / bsum;
        for (int k = 0; k < RT_FR631_N_ANGLE; ++k) {
            double q11 = 0.0, q12 = 0.0, q33 = 0.0;
            rt_value_phase_interp_eval_theta(
                component, theta[k], &q11, &q12, &q33);
            p11[k] += w * q11;
            p12[k] += w * q12;
            p33[k] += w * q33;
        }
    }

    rt_value_phase_interp_t built = {0};
    const int brc = rt_value_phase_interp_build(
        &built, theta, p11, p12, p33, RT_FR631_N_ANGLE, 0);
    if (brc != 0) return -6;

    int victim = 0;
    for (int i = 0; i < WATER_MIXED_VALUE_CACHE_CAP; ++i) {
        if (!g_wmvc[i].valid) { victim = i; break; }
        if (g_wmvc[i].stamp < g_wmvc[victim].stamp) victim = i;
    }
    water_mixed_value_cache_entry_t *e = &g_wmvc[victim];
    water_mixed_value_cache_entry_clear(e);
    for (int i = 0; i < 3; ++i) {
        const char *p = path[i] ? path[i] : "";
        if (snprintf(e->path[i], sizeof e->path[i], "%s", p) >=
            (int)sizeof e->path[i]) {
            rt_value_phase_interp_free(&built);
            return -7;
        }
        e->b[i] = b[i];
    }
    e->valid = 1;
    e->lambda_nm = lambda_nm;
    e->truncation_mode = truncation_mode;
    e->n_component = n_component;
    e->stamp = ++g_wmvc_tick;
    e->signature = water_mixed_value_signature(
        path, b, lambda_nm, truncation_mode);
    e->phase = built;
    ++g_wmvc_build_count;

    *phase_out = &e->phase;
    if (n_component_out) *n_component_out = n_component;
    if (signature_out) *signature_out = e->signature;
    const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
    if (dump && dump[0] && strcmp(dump, "0") != 0)
        fprintf(stderr,
                "WCMIXVALUE build lambda=%.10g trunc=%d components=%d n=%d fr631=%d "
                "g=%.12g bb_b=%.12g sig=%016llx builds=%llu hits=%llu\n",
                lambda_nm, truncation_mode, n_component,
                e->phase.n, e->phase.is_fr631,
                e->phase.g_asym, e->phase.bb_b_ratio,
                (unsigned long long)e->signature,
                g_wmvc_build_count, g_wmvc_hit_count);
    return 0;
}

/* ---------------------------------------------------------------------------
 * OCRT organic/TSM component-phase cache.
 *
 * A constituent mixture may use three independent vector phase files
 * (EAP phytoplankton, Stramski detritus and Ahn mineral).  Their normalized
 * P11/P12/P33 -> beta/gamma/alpha/zeta transforms depend only on the phase
 * file, wavelength and moment settings, never on Chl/TSM concentration or on
 * SOS order.  Cache several immutable decompositions per worker so:
 *   - no .mie I/O occurs inside the SOS loop;
 *   - no vector-moment transform is repeated for a different incident beam;
 *   - alternating species/wavelength batch rows do not thrash a single slot.
 * Concentration enters only in the cheap b-weighted linear mixture below. */
#define WATER_COMPONENT_PHASE_CACHE_CAP 24
typedef struct {
    int valid;
    char path[1024];
    double lambda_nm;
    int moment_mode, nmg, lmax;
    unsigned long long stamp;
    double betal[FIXED_BULK_BETAL_LMAX + 1];
    double gammal[FIXED_BULK_BETAL_LMAX + 1];
    double alphal[FIXED_BULK_BETAL_LMAX + 1];
    double zetal[FIXED_BULK_BETAL_LMAX + 1];
    double p11_90, p12_90, p33_90, bb_b_ratio;
} water_component_phase_cache_entry_t;

static water_component_phase_cache_entry_t
    g_wcpc[WATER_COMPONENT_PHASE_CACHE_CAP];
static unsigned long long g_wcpc_tick = 0;
static unsigned long long g_wcpc_build_count = 0;
static unsigned long long g_wcpc_hit_count = 0;
#ifdef _OPENMP
#pragma omp threadprivate(g_wcpc, g_wcpc_tick, g_wcpc_build_count, g_wcpc_hit_count)
#endif

static int water_component_phase_cache_match(
        const water_component_phase_cache_entry_t *e,
        const char *path, double lambda_nm, int moment_mode, int nmg, int lmax) {
    return e && e->valid && path &&
           strcmp(e->path, path) == 0 && e->lambda_nm == lambda_nm &&
           e->moment_mode == moment_mode && e->nmg == nmg && e->lmax == lmax;
}

static int water_component_phase_moments(
        const char *path, double lambda_nm, int moment_mode, int nmg, int lmax,
        double *betal, double *gammal, double *alphal, double *zetal,
        double *p11_90, double *p12_90, double *p33_90, double *bb_b_ratio,
        int *cache_hit) {
    if (cache_hit) *cache_hit = 0;
    if (!path || !path[0] || !betal || !gammal || !alphal || !zetal ||
        lmax < 0 || lmax > FIXED_BULK_BETAL_LMAX) return -1;
    if (nmg <= 1) nmg = 60;

    const char *disable_cache_env =
        getenv("OCRT_DISABLE_WATER_COMPONENT_PHASE_CACHE");
    const int disable_cache = disable_cache_env && disable_cache_env[0] &&
                              strcmp(disable_cache_env, "0") != 0;

    for (int i = 0; !disable_cache && i < WATER_COMPONENT_PHASE_CACHE_CAP; ++i) {
        water_component_phase_cache_entry_t *e = &g_wcpc[i];
        if (!water_component_phase_cache_match(e, path, lambda_nm,
                                               moment_mode, nmg, lmax)) continue;
        e->stamp = ++g_wcpc_tick;
        ++g_wcpc_hit_count;
        memcpy(betal,  e->betal,  (size_t)(lmax + 1) * sizeof(double));
        memcpy(gammal, e->gammal, (size_t)(lmax + 1) * sizeof(double));
        memcpy(alphal, e->alphal, (size_t)(lmax + 1) * sizeof(double));
        memcpy(zetal,  e->zetal,  (size_t)(lmax + 1) * sizeof(double));
        if (p11_90) *p11_90 = e->p11_90;
        if (p12_90) *p12_90 = e->p12_90;
        if (p33_90) *p33_90 = e->p33_90;
        if (bb_b_ratio) *bb_b_ratio = e->bb_b_ratio;
        if (cache_hit) *cache_hit = 1;
        const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
        if (dump && dump[0] && strcmp(dump, "0") != 0)
            fprintf(stderr,
                    "WCPHASE hit path=%s lambda=%.10g L=%d mode=%d nmg=%d builds=%llu hits=%llu\n",
                    path, lambda_nm, lmax, moment_mode, nmg,
                    g_wcpc_build_count, g_wcpc_hit_count);
        return 0;
    }

    mie_data_t mie = {0};
    aerosol_phase_interp_t aip;
    memset(&aip, 0, sizeof aip);
    if (read_mie_file(path, &mie) != 0) return -2;
    if (build_aerosol_interpolators(&mie, lambda_nm * 1.0e-3, &aip) != 0) {
        mie_data_free(&mie);
        return -3;
    }
    const int n_ang = mie.n_ang;
    double *p11 = (double *)malloc((size_t)n_ang * sizeof(double));
    double *p12 = (double *)malloc((size_t)n_ang * sizeof(double));
    double *p33 = (double *)malloc((size_t)n_ang * sizeof(double));
    double *theta = (double *)malloc((size_t)n_ang * sizeof(double));
    if (!p11 || !p12 || !p33 || !theta) {
        free(p11); free(p12); free(p33); free(theta);
        aerosol_phase_interp_free(&aip); mie_data_free(&mie);
        return -4;
    }
    for (int i = 0; i < n_ang; ++i) {
        theta[i] = mie.angles[i];
        eval_aerosol_phase(&aip, theta[i], &p11[i], &p12[i], &p33[i]);
    }

    int rc;
    if (moment_mode == 1) {
        rc = rt_aerosol_compute_vector_legendre_gauss(
            p11, p12, p33, theta, n_ang, lmax, nmg,
            betal, gammal, alphal, zetal);
    } else {
        rc = rt_aerosol_compute_vector_legendre(
            p11, p12, p33, theta, n_ang, lmax, 4096,
            betal, gammal, alphal, zetal);
    }
    double q11 = 0.0, q12 = 0.0, q33 = 0.0;
    eval_aerosol_phase(&aip, 90.0, &q11, &q12, &q33);
    const double qbb = aip.bb_b_ratio;
    free(p11); free(p12); free(p33); free(theta);
    aerosol_phase_interp_free(&aip);
    mie_data_free(&mie);
    if (rc != 0) return -5;

    if (!disable_cache) {
        int victim = 0;
        for (int i = 0; i < WATER_COMPONENT_PHASE_CACHE_CAP; ++i) {
            if (!g_wcpc[i].valid) { victim = i; break; }
            if (g_wcpc[i].stamp < g_wcpc[victim].stamp) victim = i;
        }
        water_component_phase_cache_entry_t *e = &g_wcpc[victim];
        memset(e, 0, sizeof *e);
        if (snprintf(e->path, sizeof e->path, "%s", path) >= (int)sizeof e->path)
            return -6;
        e->valid = 1;
        e->lambda_nm = lambda_nm;
        e->moment_mode = moment_mode;
        e->nmg = nmg;
        e->lmax = lmax;
        e->stamp = ++g_wcpc_tick;
        memcpy(e->betal,  betal,  (size_t)(lmax + 1) * sizeof(double));
        memcpy(e->gammal, gammal, (size_t)(lmax + 1) * sizeof(double));
        memcpy(e->alphal, alphal, (size_t)(lmax + 1) * sizeof(double));
        memcpy(e->zetal,  zetal,  (size_t)(lmax + 1) * sizeof(double));
        e->p11_90 = q11; e->p12_90 = q12; e->p33_90 = q33;
        e->bb_b_ratio = qbb;
    }
    ++g_wcpc_build_count;
    if (p11_90) *p11_90 = q11;
    if (p12_90) *p12_90 = q12;
    if (p33_90) *p33_90 = q33;
    if (bb_b_ratio) *bb_b_ratio = qbb;
    const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
    if (dump && dump[0] && strcmp(dump, "0") != 0)
        fprintf(stderr,
                "WCPHASE build path=%s lambda=%.10g L=%d mode=%d nmg=%d builds=%llu hits=%llu\n",
                path, lambda_nm, lmax, moment_mode, nmg,
                g_wcpc_build_count, g_wcpc_hit_count);
    return 0;
}

/* The Huot spectral allocation needs the phase backscatter ratio at 550 nm
 * as its normalization anchor.  Computing a full beta/gamma/alpha/zeta
 * decomposition only to obtain that scalar would waste work, so cache a
 * lightweight P11/PCHIP integration separately.  This function is called
 * during medium preparation, never from an SOS order loop. */
#define WATER_COMPONENT_RATIO_CACHE_CAP 24
typedef struct {
    int valid;
    char path[1024];
    double lambda_nm;
    double bb_b_ratio;
    unsigned long long stamp;
} water_component_ratio_cache_entry_t;

static water_component_ratio_cache_entry_t
    g_wcrc[WATER_COMPONENT_RATIO_CACHE_CAP];
static unsigned long long g_wcrc_tick = 0;
static unsigned long long g_wcrc_build_count = 0;
static unsigned long long g_wcrc_hit_count = 0;
#ifdef _OPENMP
#pragma omp threadprivate(g_wcrc, g_wcrc_tick, g_wcrc_build_count, g_wcrc_hit_count)
#endif

static int water_component_phase_ratio(const char *path, double lambda_nm,
                                       double *bb_b_ratio, int *cache_hit) {
    if (cache_hit) *cache_hit = 0;
    if (!path || !path[0] || !bb_b_ratio || !isfinite(lambda_nm)) return -1;
    const char *disable_cache_env =
        getenv("OCRT_DISABLE_WATER_COMPONENT_PHASE_CACHE");
    const int disable_cache = disable_cache_env && disable_cache_env[0] &&
                              strcmp(disable_cache_env, "0") != 0;
    for (int i = 0; !disable_cache && i < WATER_COMPONENT_RATIO_CACHE_CAP; ++i) {
        water_component_ratio_cache_entry_t *e = &g_wcrc[i];
        if (!e->valid || e->lambda_nm != lambda_nm || strcmp(e->path, path) != 0)
            continue;
        e->stamp = ++g_wcrc_tick;
        ++g_wcrc_hit_count;
        *bb_b_ratio = e->bb_b_ratio;
        if (cache_hit) *cache_hit = 1;
        const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
        if (dump && dump[0] && strcmp(dump, "0") != 0)
            fprintf(stderr,
                    "WCRATIO hit path=%s lambda=%.10g ratio=%.12g builds=%llu hits=%llu\n",
                    path, lambda_nm, *bb_b_ratio,
                    g_wcrc_build_count, g_wcrc_hit_count);
        return 0;
    }

    const rt_value_phase_interp_t *phase = NULL;
    int phase_hit = 0;
    const int prc = water_component_value_phase_get(
        path, lambda_nm, 0, &phase, NULL, NULL, NULL, &phase_hit);
    if (prc != 0 || !phase) return -2;
    const double ratio = phase->bb_b_ratio;
    if (!(ratio > 0.0 && ratio < 0.5) || !isfinite(ratio)) return -3;
    if (!disable_cache) {
        int victim = 0;
        for (int i = 0; i < WATER_COMPONENT_RATIO_CACHE_CAP; ++i) {
            if (!g_wcrc[i].valid) { victim = i; break; }
            if (g_wcrc[i].stamp < g_wcrc[victim].stamp) victim = i;
        }
        water_component_ratio_cache_entry_t *e = &g_wcrc[victim];
        memset(e, 0, sizeof *e);
        if (snprintf(e->path, sizeof e->path, "%s", path) >= (int)sizeof e->path)
            return -4;
        e->valid = 1;
        e->lambda_nm = lambda_nm;
        e->bb_b_ratio = ratio;
        e->stamp = ++g_wcrc_tick;
    }
    ++g_wcrc_build_count;
    *bb_b_ratio = ratio;
    if (cache_hit) *cache_hit = phase_hit;
    const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
    if (dump && dump[0] && strcmp(dump, "0") != 0)
        fprintf(stderr,
                "WCRATIO build path=%s lambda=%.10g ratio=%.12g source=%s builds=%llu hits=%llu\n",
                path, lambda_nm, ratio, phase_hit ? "value-cache" : "value-build",
                g_wcrc_build_count, g_wcrc_hit_count);
    return 0;
}

static int water_component_phase_ratio_prime(
        const char *path, double lambda_nm, int moment_mode, int nmg, int lmax,
        double *bb_b_ratio) {
    (void)moment_mode; (void)nmg; (void)lmax;
    return water_component_phase_ratio(path, lambda_nm, bb_b_ratio, NULL);
}


/* Cache the final b-weighted particulate vector phase as well as the
 * per-component decompositions.  This removes even the small O(L) mixture
 * rebuild from repeated sky-beam / angular-grid solves at the same water
 * state, and makes the S-009 diagnostic one-per-unique-state. */
#define WATER_MIXED_PHASE_CACHE_CAP 24
typedef struct {
    int valid;
    char path[3][1024];
    double lambda_nm;
    double b[3];
    int moment_mode, nmg, lmax, n_component;
    double recon_min;
    unsigned long long stamp;
    double betal[FIXED_BULK_BETAL_LMAX + 1];
    double gammal[FIXED_BULK_BETAL_LMAX + 1];
    double alphal[FIXED_BULK_BETAL_LMAX + 1];
    double zetal[FIXED_BULK_BETAL_LMAX + 1];
} water_mixed_phase_cache_entry_t;

static water_mixed_phase_cache_entry_t
    g_wmixc[WATER_MIXED_PHASE_CACHE_CAP];
static unsigned long long g_wmixc_tick = 0;
static unsigned long long g_wmixc_build_count = 0;
static unsigned long long g_wmixc_hit_count = 0;
#ifdef _OPENMP
#pragma omp threadprivate(g_wmixc, g_wmixc_tick, g_wmixc_build_count, g_wmixc_hit_count)
#endif

static int water_mixed_phase_key_match(
        const water_mixed_phase_cache_entry_t *e,
        const char *const path[3], const double b[3], double lambda_nm,
        int moment_mode, int nmg, int lmax) {
    if (!e || !e->valid || e->lambda_nm != lambda_nm ||
        e->moment_mode != moment_mode || e->nmg != nmg || e->lmax != lmax)
        return 0;
    for (int i = 0; i < 3; ++i) {
        const char *p = path[i] ? path[i] : "";
        if (e->b[i] != b[i] || strcmp(e->path[i], p) != 0) return 0;
    }
    return 1;
}

static int water_mixed_phase_cache_load(
        const char *const path[3], const double b[3], double lambda_nm,
        int moment_mode, int nmg, int lmax,
        double *betal, double *gammal, double *alphal, double *zetal,
        int *n_component, double *recon_min) {
    const char *disable_cache_env =
        getenv("OCRT_DISABLE_WATER_COMPONENT_PHASE_CACHE");
    const int disable_cache = disable_cache_env && disable_cache_env[0] &&
                              strcmp(disable_cache_env, "0") != 0;
    if (disable_cache) return 0;
    for (int i = 0; i < WATER_MIXED_PHASE_CACHE_CAP; ++i) {
        water_mixed_phase_cache_entry_t *e = &g_wmixc[i];
        if (!water_mixed_phase_key_match(e, path, b, lambda_nm,
                                         moment_mode, nmg, lmax)) continue;
        e->stamp = ++g_wmixc_tick;
        ++g_wmixc_hit_count;
        memcpy(betal, e->betal, (size_t)(lmax + 1) * sizeof(double));
        memcpy(gammal, e->gammal, (size_t)(lmax + 1) * sizeof(double));
        memcpy(alphal, e->alphal, (size_t)(lmax + 1) * sizeof(double));
        memcpy(zetal, e->zetal, (size_t)(lmax + 1) * sizeof(double));
        if (n_component) *n_component = e->n_component;
        if (recon_min) *recon_min = e->recon_min;
        const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
        if (dump && dump[0] && strcmp(dump, "0") != 0)
            fprintf(stderr,
                    "WCMIX hit lambda=%.10g L=%d mode=%d nmg=%d components=%d builds=%llu hits=%llu\n",
                    lambda_nm, lmax, moment_mode, nmg, e->n_component,
                    g_wmixc_build_count, g_wmixc_hit_count);
        return 1;
    }
    return 0;
}

static int water_mixed_phase_cache_store(
        const char *const path[3], const double b[3], double lambda_nm,
        int moment_mode, int nmg, int lmax,
        const double *betal, const double *gammal,
        const double *alphal, const double *zetal,
        int n_component, double recon_min) {
    const char *disable_cache_env =
        getenv("OCRT_DISABLE_WATER_COMPONENT_PHASE_CACHE");
    if (disable_cache_env && disable_cache_env[0] &&
        strcmp(disable_cache_env, "0") != 0) return 0;
    int victim = 0;
    for (int i = 0; i < WATER_MIXED_PHASE_CACHE_CAP; ++i) {
        if (!g_wmixc[i].valid) { victim = i; break; }
        if (g_wmixc[i].stamp < g_wmixc[victim].stamp) victim = i;
    }
    water_mixed_phase_cache_entry_t *e = &g_wmixc[victim];
    memset(e, 0, sizeof *e);
    for (int i = 0; i < 3; ++i) {
        const char *p = path[i] ? path[i] : "";
        if (snprintf(e->path[i], sizeof e->path[i], "%s", p) >=
            (int)sizeof e->path[i]) return -1;
        e->b[i] = b[i];
    }
    e->valid = 1;
    e->lambda_nm = lambda_nm;
    e->moment_mode = moment_mode;
    e->nmg = nmg;
    e->lmax = lmax;
    e->n_component = n_component;
    e->recon_min = recon_min;
    e->stamp = ++g_wmixc_tick;
    memcpy(e->betal, betal, (size_t)(lmax + 1) * sizeof(double));
    memcpy(e->gammal, gammal, (size_t)(lmax + 1) * sizeof(double));
    memcpy(e->alphal, alphal, (size_t)(lmax + 1) * sizeof(double));
    memcpy(e->zetal, zetal, (size_t)(lmax + 1) * sizeof(double));
    ++g_wmixc_build_count;
    const char *dump = getenv("OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE");
    if (dump && dump[0] && strcmp(dump, "0") != 0)
        fprintf(stderr,
                "WCMIX build lambda=%.10g L=%d mode=%d nmg=%d components=%d builds=%llu hits=%llu\n",
                lambda_nm, lmax, moment_mode, nmg, n_component,
                g_wmixc_build_count, g_wmixc_hit_count);
    return 0;
}

static double water_phase_reconstruction_min_back(const double *betal, int lmax) {
    static const double mu_back[] = {
        -1.0, -0.93969262078591, -0.86602540378444,
        -0.70710678118655, -0.5
    };
    double min_value = 1.0e300;
    if (!betal || lmax < 0) return min_value;
    for (size_t i = 0; i < sizeof mu_back / sizeof mu_back[0]; ++i) {
        const double x = mu_back[i];
        double sum = betal[0];
        if (lmax >= 1) sum += betal[1] * x;
        double pm1 = 1.0, p0 = x;
        for (int l = 2; l <= lmax; ++l) {
            const double pl = ((2.0 * l - 1.0) * x * p0 -
                               (l - 1.0) * pm1) / (double)l;
            sum += betal[l] * pl;
            pm1 = p0;
            p0 = pl;
        }
        if (sum < min_value) min_value = sum;
    }
    return min_value;
}

#ifdef _OPENMP
/* The water beam loop within a single case is serial, so a single g_wpkc is
 * race-free per case.  When the DRIVER parallelizes across geometries/cases
 * with OpenMP threads (multi-geometry / multi-case LUT runs), each thread must
 * own its cache: threadprivate gives every thread an independent g_wpkc.  The
 * signature additionally keys on the view direction (mu_view), so a thread
 * reusing its cache across geometries rebuilds when the view node changes
 * (only the within-case ~N_mu-beam redundancy is meant to be cached). */
#pragma omp threadprivate(g_wpkc)
#endif

static void wpkc_free(void) {
    if (g_wpkc.kernel_m) {
        for (int m = 0; m < g_wpkc.m_count; ++m) free(g_wpkc.kernel_m[m]);
        free(g_wpkc.kernel_m);
        g_wpkc.kernel_m = NULL;
    }
    g_wpkc.valid = 0;
    g_wpkc.m_count = 0;
}
static int wpkc_match(double lambda_nm, double chl, double min_c, double adom440, double adomS,
                      double n_water, double T_water, double nphi_d, double mu_view, int n_mu, int m_count,
                      int method, int lut_mode, const char *lut_path) {
    if (!g_wpkc.valid) return 0;
    if (g_wpkc.n_mu != n_mu || g_wpkc.m_count != m_count || g_wpkc.method != method ||
        g_wpkc.lut_mode != lut_mode) return 0;
    if (g_wpkc.lambda_nm != lambda_nm || g_wpkc.chl != chl || g_wpkc.min_c != min_c) return 0;
    if (g_wpkc.adom440 != adom440 || g_wpkc.adomS != adomS) return 0;
    if (g_wpkc.n_water != n_water || g_wpkc.T_water != T_water || g_wpkc.nphi_d != nphi_d) return 0;
    if (g_wpkc.mu_view != mu_view) return 0;
    if (lut_mode && strncmp(g_wpkc.lut_path, lut_path ? lut_path : "", sizeof g_wpkc.lut_path) != 0) return 0;
    return 1;
}
static int wpkc_reset(double lambda_nm, double chl, double min_c, double adom440, double adomS,
                      double n_water, double T_water, double nphi_d, double mu_view, int n_mu, int m_count,
                      int method, int lut_mode, const char *lut_path) {
    wpkc_free();
    const int dirs = 2 * n_mu + 1;
    g_wpkc.kernel_m = (double**)calloc((size_t)m_count, sizeof(double*));
    if (!g_wpkc.kernel_m) return -1;
    for (int m = 0; m < m_count; ++m) {
        g_wpkc.kernel_m[m] = (double*)calloc((size_t)(n_mu + 1) * (size_t)dirs, sizeof(double));
        if (!g_wpkc.kernel_m[m]) { wpkc_free(); return -1; }
    }
    g_wpkc.lambda_nm = lambda_nm; g_wpkc.chl = chl; g_wpkc.min_c = min_c;
    g_wpkc.adom440 = adom440; g_wpkc.adomS = adomS;
    g_wpkc.n_water = n_water; g_wpkc.T_water = T_water; g_wpkc.nphi_d = nphi_d; g_wpkc.mu_view = mu_view;
    g_wpkc.n_mu = n_mu; g_wpkc.m_count = m_count; g_wpkc.method = method; g_wpkc.lut_mode = lut_mode;
    g_wpkc.lut_path[0] = '\0';
    if (lut_mode && lut_path) { strncpy(g_wpkc.lut_path, lut_path, sizeof g_wpkc.lut_path - 1); }
    g_wpkc.valid = 1;
    return 0;
}
static void wpkc_store_m(int m, const rt_legendre_workspace_t *ws, int n_mu) {
    if (!g_wpkc.valid || m < 0 || m >= g_wpkc.m_count) return;
    const int dirs = 2 * n_mu + 1;
    double *dst = g_wpkc.kernel_m[m];
    for (int j = 0; j <= n_mu; ++j)
        for (int k = -n_mu; k <= n_mu; ++k)
            dst[j * dirs + (k + n_mu)] = ws->phase_fourier_m[j][k];
}
static void wpkc_load_m(int m, rt_legendre_workspace_t *ws, int n_mu) {
    if (!g_wpkc.valid || m < 0 || m >= g_wpkc.m_count) return;
    const int dirs = 2 * n_mu + 1;
    const double *src = g_wpkc.kernel_m[m];
    for (int j = 0; j <= n_mu; ++j)
        for (int k = -n_mu; k <= n_mu; ++k)
            ws->phase_fourier_m[j][k] = src[j * dirs + (k + n_mu)];
}


/* ============================================================================
 * v1.09 commit #16: FULL-GRID single-solve cache ("solve once, extract many").
 * Field solution is view-independent (raa: reconstruction-only; vza with
 * view_as_node=0: post-solve interpolation only) -> ONE solve per
 * (sza, band, IOP) serves every grid row.  On a key hit the per-m SOLVE
 * section is skipped (fields restored from cache) and control jumps to the
 * UNMODIFIED extraction tail of the same loop, so replay parity with a cold
 * solve is structural (no duplicated extraction code).  Gated by
 * OCRT_WATER_GRID_CACHE (full-grid driver sets it; "0" disables); requires
 * view_as_node==0 (driver enforces; grid vza ARE quadrature nodes).
 * ========================================================================== */
typedef struct {
    int    valid, m_count, nt, n_mu, dirs;
    int    max_orders_seen_c, all_conv_c;
    double worst_resid_c;
    double sza, lam, Twc, Sgk, nw, Fs, wind, tol;
    int    nmw, mmw, nlw, maxit, vasn, wpk, nvv;  /* nvv: #20 view-node count */
    double a_t, b_t, bb_t, taumax;
    /* Effective IOP and constituent phase signature.  The legacy key only
     * carried fixed-bulk inputs, so two full-grid runs with different Chl/TSM
     * could alias when both used the ordinary constituent branch. */
    double a_iop, b_iop, bb_iop;
    double b_phyto, b_detritus, b_mineral;
    int constituent_mode, constituent_model, organic_group, tsm_species;
    int moment_mode, moment_nmg, phase_lmax, mie_trunc, mie_ss;
    int value_phase_spline, value_kernel_pol;
    char   wmp[512], fbl[512];
    double *tot;   /* [m][3][(nt+1)*dirs] */
    double *prim0; /* [m][3][dirs] (k=0 row) */
} ocrt_grid_cache_t;
/* v1.09-opt S6a: the beam loop issues ~25 distinct cache keys per case
 * (24 sky-beam (sza_eq, F_sun_eq) pairs + the direct solve); a single slot
 * makes them evict each other every row.  FASTK: 28-slot LRU (never
 * evicting the entry just used); non-FASTK: 1 slot = legacy behavior.
 * Replay parity is structural (unmodified extraction tail), so slot count
 * cannot change any output value - only which calls replay vs re-solve. */
#ifdef OCRT_FAST_KERNELS
/* distinct keys per batch row observed: base solve beams (default nmw=48,
 * 48 keys) + grid beams (row nmw, e.g. 24) + direct solves -> ~75.  Slots
 * are lazily allocated (an unused slot costs a struct, no field memory),
 * so 96 is safe; override via OCRT_WGC_SLOTS if memory-tight. */
#define OCRT_WGC_SLOTS 96
#else
#define OCRT_WGC_SLOTS 1
#endif
static ocrt_grid_cache_t g_wgc_slots[OCRT_WGC_SLOTS];
static unsigned long     g_wgc_stamp[OCRT_WGC_SLOTS];
static unsigned long     g_wgc_tick = 0;
static int               ocrt_wgc_cur_i = 0;
#pragma omp threadprivate(g_wgc_slots, g_wgc_stamp, g_wgc_tick, ocrt_wgc_cur_i)
#define g_grid_cache (g_wgc_slots[ocrt_wgc_cur_i])

/* v1.09 commit #19 (OpenMP full-grid batch): invalidate the calling thread's
 * grid cache.  Called at the start of every full-grid run so a thread that
 * processes several batch rows can never replay a previous row's fields,
 * independent of cache-key completeness. */
void rt_water_grid_cache_reset(void) {
    for (int ocrt_i = 0; ocrt_i < OCRT_WGC_SLOTS; ++ocrt_i)
        g_wgc_slots[ocrt_i].valid = 0;
    /* v1.09 #19: also drop the wavelength-keyed phase-kernel cache.  Its key
     * does NOT include the water phase FILE identity, so two batch rows with
     * the same wavelength but different OCRT/IOP .mie phase files would falsely hit.
     * Buffers are freed by the next wpkc_reset(); clearing valid is enough. */
    g_wpkc.valid = 0;
}
static void grid_cache_free(void){ free(g_grid_cache.tot); free(g_grid_cache.prim0); memset(&g_grid_cache,0,sizeof g_grid_cache); }

/* v1.10 D3-PROD beam cache (2026-07-10): cross-process persistence of ONE
 * grid-cache entry (the beam-only water solve of the PROD path).  The beam
 * solve depends only on (case IOPs, wl, wind, sza, F_sun, solver config) -
 * NOT on the atmosphere profile - so across the 56 aerosol/AOD/height combos
 * sharing those, the cold solve amortizes.  Env OCRT_D3_BEAM_CACHE=path:
 * on entry (before lookup) the file is loaded once per process into slot 0;
 * after a cold store, if the file does not exist, the entry is dumped.  The
 * EXISTING key comparison then validates the loaded entry against the run
 * parameters - a mismatched file simply misses (classical cold, no risk). */
static int ocrt_beam_cache_load(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    ocrt_grid_cache_t e; memset(&e, 0, sizeof e);
    size_t nsc = fread(&e, sizeof(ocrt_grid_cache_t), 1, f);
    if (nsc != 1) { fclose(f); return -2; }
    size_t ntot = (size_t)e.m_count * 3u * (size_t)(e.nt + 1) * (size_t)e.dirs;
    size_t npr  = (size_t)e.m_count * 3u * (size_t)e.dirs;
    e.tot   = (double*)malloc(ntot * sizeof(double));
    e.prim0 = (double*)malloc(npr  * sizeof(double));
    if (!e.tot || !e.prim0 ||
        fread(e.tot,   sizeof(double), ntot, f) != ntot ||
        fread(e.prim0, sizeof(double), npr,  f) != npr) {
        free(e.tot); free(e.prim0); fclose(f); return -3;
    }
    fclose(f);
    ocrt_wgc_cur_i = 0;
    grid_cache_free();
    g_grid_cache = e;
    g_grid_cache.valid = 1;
    return 0;
}
static void ocrt_beam_cache_dump(const char *path) {
    if (!g_grid_cache.valid || !g_grid_cache.tot || !g_grid_cache.prim0) return;
    FILE *f = fopen(path, "wb");
    if (!f) return;
    size_t ntot = (size_t)g_grid_cache.m_count * 3u * (size_t)(g_grid_cache.nt + 1) * (size_t)g_grid_cache.dirs;
    size_t npr  = (size_t)g_grid_cache.m_count * 3u * (size_t)g_grid_cache.dirs;
    fwrite(&g_grid_cache, sizeof(ocrt_grid_cache_t), 1, f);
    fwrite(g_grid_cache.tot,   sizeof(double), ntot, f);
    fwrite(g_grid_cache.prim0, sizeof(double), npr,  f);
    fclose(f);
    fprintf(stderr, "[D3beam] cache dumped: %s (m=%d nt=%d dirs=%d, %.1f MB)\n",
            path, g_grid_cache.m_count, g_grid_cache.nt, g_grid_cache.dirs,
            (ntot + npr) * 8.0 / 1048576.0);
}

int rt_water_rt_sos_pure(double sza_deg_air, double vza_deg_air, double raa_deg,
                          double lambda_nm,
                          double T_water_C, double S_water_gkg,
                          double n_water,
                          double F_sun,
                          const rt_water_iop_lut_t* aw_lut,
                          const rt_water_iop_psi_T_lut_t* psi_T_lut,
                          const rt_water_rt_options_t* opts,
                          rt_water_rt_result_t* result)
{
    if (!rt_spectral_wavelength_supported(lambda_nm)) {
        fprintf(stderr,
                "rt_water_rt: wavelength %.12g nm is outside supported %.0f-%.0f nm\n",
                lambda_nm, RT_SPECTRAL_MIN_NM, RT_SPECTRAL_MAX_NM);
        return -1;
    }

    rt_water_env_init();   /* #21: idempotent; snapshot ready before any solve-path read */
    (void)S_water_gkg;  /* salinity not used in B.3 (psi_S deferred) */

    if (!result || !aw_lut) return -1;
    if (!(n_water > 1.0 && n_water < 2.0)) return -1;
    if (!(F_sun > 0.0)) return -1;

    rt_water_rt_options_t o = opts ? *opts : rt_water_rt_options_default();

    /* 1. IOP at requested λ, T — component evaluators
     *
     * 통일된 인터페이스 (rt_iop_t per component → rt_iop_total) 로 합산.
     * CDOM 은 산란이 없으므로 b=bb=0; 향후 phyto/NAP 추가 시 같은 패턴.
     * 본 합산 결과 (a_tot, b_tot, bb_tot) 가 RT 솔버 입력이며,
     * Mueller matrix (phase function) 는 b_tot 안에서 pure-water 만 차지하므로
     * 현재 가중치 처리 불필요. */
    int extrap_dummy = 0;
    rt_iop_components_t comp;
    memset(&comp, 0, sizeof comp);
    int rc_pw = rt_iop_pure_water_eval(aw_lut, psi_T_lut,
                                        lambda_nm, T_water_C,
                                        &comp.pure_water, &extrap_dummy);
    if (rc_pw != 0) return -2;

    rt_iop_cdom_params_t cdom_p = {
        .a440          = o.a_cdom_440_m_inv,
        .S_nm_inv      = o.S_cdom_nm_inv,
        .lambda_ref_nm = o.cdom_ref_lambda_nm
    };
    rt_iop_cdom_eval(&cdom_p, lambda_nm, &comp.cdom);

    const char *organic_phyto_mie_path = NULL;
    const char *organic_detritus_mie_path = NULL;
    const char *ahn_mie_path_used = NULL;
    const int constituent_model_ocrt =
        o.ccrr_mode && o.water_constituent_model == RT_WATER_CONSTITUENT_OCRT;
    const int constituent_model_ccrr =
        o.ccrr_mode && o.water_constituent_model == RT_WATER_CONSTITUENT_CCRR;

    if (constituent_model_ocrt) {
        if (o.organic_phyto_scattering) {
            fprintf(stderr,
                    "rt_water_rt: EAP phytoplankton scattering is disabled by current project policy; "
                    "use common absorption plus detritus phase only\n");
            return -2;
        }
        if (o.ccrr_chl_mg_m3 > 0.0) {
            organic_phyto_group_t group = (organic_phyto_group_t)o.organic_phyto_group;
            if (!rt_iop_organic_ready() || !rt_iop_organic_wavelength_supported(lambda_nm)) {
                fprintf(stderr,
                        "rt_water_rt: OCRT organic model unavailable at lambda=%.10g nm (ready=%d, range %.0f-%.0f nm)\n",
                        lambda_nm, rt_iop_organic_ready(),
                        ORGANIC_WAVELENGTH_MIN_NM, ORGANIC_WAVELENGTH_MAX_NM);
                return -2;
            }
            organic_detritus_mie_path = rt_iop_organic_detritus_phase_path();
            if (!organic_detritus_mie_path || !organic_detritus_mie_path[0]) {
                fprintf(stderr, "rt_water_rt: OCRT detritus phase path resolution failed\n");
                return -2;
            }
            if (!rt_iop_organic_detritus_wavelength_supported(lambda_nm)) {
                fprintf(stderr,
                        "rt_water_rt: active Chl-linked detritus phase supports %.0f-%.0f nm; "
                        "lambda=%.10g nm is unavailable. The 330-1100 recipe candidate "
                        "failed the L=200/SOS acceptance gate and is not active.\n",
                        rt_iop_organic_detritus_wavelength_min_nm(),
                        rt_iop_organic_detritus_wavelength_max_nm(), lambda_nm);
                return -2;
            }

            const int organic_L = (o.fixed_bulk_lmax > 0 &&
                                   o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX)
                                      ? o.fixed_bulk_lmax : 200;
            const int organic_mode = o.water_mie_moment_mode;
            const int organic_nmg = (o.water_mie_moment_n_mu > 1)
                                        ? o.water_mie_moment_n_mu : 60;
            double det_ratio = 0.0, det_ratio550 = 0.0;
            int drp = water_component_phase_ratio_prime(
                organic_detritus_mie_path, lambda_nm,
                organic_mode, organic_nmg, organic_L, &det_ratio);
            if (drp != 0) {
                fprintf(stderr,
                        "rt_water_rt: OCRT detritus phase precomputation failed lambda=%.10g rc=%d\n",
                        lambda_nm, drp);
                return -2;
            }
            if (fabs(lambda_nm - 550.0) <= 1.0e-12) {
                det_ratio550 = det_ratio;
            } else if (water_component_phase_ratio(organic_detritus_mie_path, 550.0,
                                                   &det_ratio550, NULL) != 0) {
                fprintf(stderr,
                        "rt_water_rt: OCRT detritus 550-nm phase-ratio anchor evaluation failed\n");
                return -2;
            }

            int prc = 0;
            if (o.organic_phyto_scattering) {
                organic_phyto_mie_path = rt_iop_organic_phyto_phase_path(group);
                if (!organic_phyto_mie_path || !organic_phyto_mie_path[0]) {
                    fprintf(stderr, "rt_water_rt: OCRT phytoplankton phase path resolution failed\n");
                    return -2;
                }
                double phyto_ratio = 0.0, phyto_ratio550 = 0.0;
                int prp = water_component_phase_ratio_prime(
                    organic_phyto_mie_path, lambda_nm,
                    organic_mode, organic_nmg, organic_L, &phyto_ratio);
                if (prp != 0) {
                    fprintf(stderr,
                            "rt_water_rt: OCRT phytoplankton phase precomputation failed group=%s lambda=%.10g rc=%d\n",
                            rt_iop_organic_group_name(group), lambda_nm, prp);
                    return -2;
                }
                if (fabs(lambda_nm - 550.0) <= 1.0e-12) {
                    phyto_ratio550 = phyto_ratio;
                } else if (water_component_phase_ratio(organic_phyto_mie_path, 550.0,
                                                       &phyto_ratio550, NULL) != 0) {
                    fprintf(stderr,
                            "rt_water_rt: OCRT phytoplankton 550-nm phase-ratio anchor evaluation failed\n");
                    return -2;
                }
                prc = rt_iop_eap_phyto_eval_with_phase_ratios(
                    lambda_nm, o.ccrr_chl_mg_m3, group,
                    phyto_ratio, phyto_ratio550, &comp.eap_phyto);
            } else {
                prc = rt_iop_eap_phyto_absorption_eval(
                    lambda_nm, o.ccrr_chl_mg_m3, group, &comp.eap_phyto.a);
                comp.eap_phyto.b = 0.0;
                comp.eap_phyto.bb = 0.0;
            }
            int drc = rt_iop_detritus_eval_with_phase_ratios(
                lambda_nm, o.ccrr_chl_mg_m3,
                o.detritus_a440_m_inv, o.detritus_slope_nm_inv,
                det_ratio, det_ratio550, &comp.detritus);
            if (prc != 0 || drc != 0) {
                fprintf(stderr,
                        "rt_water_rt: OCRT organic IOP evaluation failed group=%s Chl=%.10g lambda=%.10g prc=%d drc=%d\n",
                        rt_iop_organic_group_name(group), o.ccrr_chl_mg_m3,
                        lambda_nm, prc, drc);
                return -2;
            }
        }
        if (o.ccrr_min_g_m3 > 0.0) {
            int arc = rt_iop_ahn_mineral_eval(lambda_nm, o.ccrr_min_g_m3,
                                               (ahn_species_t)o.tsm_species,
                                               &comp.mineral, &ahn_mie_path_used);
            if (arc != 0) {
                fprintf(stderr, "rt_water_rt: Ahn TSM IOP evaluation failed species=%s C=%.10g lambda=%.10g rc=%d\n",
                        rt_iop_ahn_species_name((ahn_species_t)o.tsm_species),
                        o.ccrr_min_g_m3, lambda_nm, arc);
                return -2;
            }
            if (!ahn_mie_path_used || !ahn_mie_path_used[0]) {
                fprintf(stderr, "rt_water_rt: Ahn TSM phase path resolution failed species=%s\n",
                        rt_iop_ahn_species_name((ahn_species_t)o.tsm_species));
                return -2;
            }
            if (o.ccrr_chl_mg_m3 > 0.0) {
                /* TSM-only remains byte-compatible with v1.16.  In the newly
                 * supported Chl+TSM mixed phase, however, force bb/b to match
                 * the exact P11 representation used in the vector mixture. */
                const int mix_L = (o.fixed_bulk_lmax > 0 &&
                                   o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX)
                                      ? o.fixed_bulk_lmax : 200;
                const int mix_mode = o.water_mie_moment_mode;
                const int mix_nmg = (o.water_mie_moment_n_mu > 1)
                                        ? o.water_mie_moment_n_mu : 60;
                double tsm_ratio = 0.0;
                int trc = water_component_phase_ratio_prime(
                    ahn_mie_path_used, lambda_nm,
                    mix_mode, mix_nmg, mix_L, &tsm_ratio);
                if (trc != 0 || !(tsm_ratio > 0.0 && tsm_ratio < 0.5)) {
                    fprintf(stderr,
                            "rt_water_rt: Ahn TSM phase-ratio precomputation failed species=%s lambda=%.10g rc=%d\n",
                            rt_iop_ahn_species_name((ahn_species_t)o.tsm_species),
                            lambda_nm, trc);
                    return -2;
                }
                comp.mineral.bb = comp.mineral.b * tsm_ratio;
            }
            /* Preserve the v1.16 TSM-only route exactly.  When Chl is also
             * active, the three particulate phases are mixed explicitly below. */
            if (!(o.ccrr_chl_mg_m3 > 0.0)) {
                if (o.water_mie_phase_path && o.water_mie_phase_path[0] &&
                    strcmp(o.water_mie_phase_path, ahn_mie_path_used) != 0) {
                    fprintf(stderr, "rt_water_rt: TSM species phase conflict: '%s' vs '%s'\n",
                            o.water_mie_phase_path, ahn_mie_path_used);
                    return -2;
                }
                o.water_mie_phase_path = ahn_mie_path_used;
            }
        }
    } else if (constituent_model_ccrr) {
        if (rt_iop_ccrr_pigment_eval(lambda_nm, o.ccrr_chl_mg_m3,
                                     &comp.pigment) != 0) return -2;
        if (rt_iop_ccrr_mineral_eval(lambda_nm, o.ccrr_min_g_m3,
                                     &comp.mineral) != 0) return -2;
    }

    rt_iop_t iop = rt_iop_total(&comp);
    double a_w   = comp.pure_water.a;   /* diagnostic */
    double b_w   = comp.pure_water.b;
    double bb_w  = comp.pure_water.bb;
    double a_cdom = comp.cdom.a;
    double a_phyto  = comp.eap_phyto.a;
    double b_phyto  = comp.eap_phyto.b;
    double bb_phyto = comp.eap_phyto.bb;
    double a_det    = comp.detritus.a;
    double b_det    = comp.detritus.b;
    double bb_det   = comp.detritus.bb;
    /* Existing downstream names represent the complete Chl-linked organic
     * particulate budget in OCRT mode, and the legacy pigment budget in CCRR. */
    double a_pig  = comp.pigment.a + a_phyto + a_det;
    double b_pig  = comp.pigment.b + b_phyto + b_det;
    double bb_pig = comp.pigment.bb + bb_phyto + bb_det;
    double a_min  = comp.mineral.a;
    double b_min  = comp.mineral.b;
    double bb_min = comp.mineral.bb;
    double a_tot = iop.a;
    double b_tot = iop.b;
    double bb_tot = iop.bb;

    /* EAP phytoplankton scattering is disabled in the constituent model.
     * Keep a_phyto in the absorption budget, but remove the phytoplankton
     * scattering and backscattering from both the bulk IOP and the organic
     * particulate phase mixture.  Detritus remains the Chl-linked particle
     * scatterer.  This gate is independent of the offline 17-species phase
     * generator integrated in rt_eap_mie_phase.c. */
    if (constituent_model_ocrt && !o.organic_phyto_scattering &&
        (b_phyto > 0.0 || bb_phyto > 0.0)) {
        b_tot  -= b_phyto;
        bb_tot -= bb_phyto;
        b_pig  -= b_phyto;
        bb_pig -= bb_phyto;
        b_phyto = 0.0;
        bb_phyto = 0.0;
        if (b_tot < 0.0 && b_tot > -1.0e-14) b_tot = 0.0;
        if (bb_tot < 0.0 && bb_tot > -1.0e-14) bb_tot = 0.0;
        if (b_pig < 0.0 && b_pig > -1.0e-14) b_pig = 0.0;
        if (bb_pig < 0.0 && bb_pig > -1.0e-14) bb_pig = 0.0;
    }

    if (o.fixed_bulk_iop_mode) {
        if (!(o.fixed_a_total_m_inv >= 0.0) || !(o.fixed_b_total_m_inv > 0.0) ||
            !(o.fixed_bb_total_m_inv > 0.0) ||
            o.fixed_bb_total_m_inv >= 0.5 * o.fixed_b_total_m_inv) return -2;
        a_w = b_w = bb_w = 0.0;
        a_cdom = a_pig = b_pig = bb_pig = a_min = b_min = bb_min = 0.0;
        a_phyto = b_phyto = bb_phyto = a_det = b_det = bb_det = 0.0;
        a_tot = o.fixed_a_total_m_inv;
        b_tot = o.fixed_b_total_m_inv;
        bb_tot = o.fixed_bb_total_m_inv;
    }

    /* Fixed-bulk forward-peak delta-scaling for HydroLight validation.
     *
     * The HydroLight fixed-bulk table supplies only a_total, b_total and
     * bb_total.  Strongly forward-peaked hydrosol phases cannot be represented
     * safely by low-order Legendre HG moments: they ring and can make P11
     * negative at side/back angles.  In fixed-bulk validation mode, remove the
     * unresolved forward component and keep a backscatter-equivalent residual:
     *
     *   b_eff = 2 * bb_total;  phase_eff = isotropic;  c_eff = a_total + b_eff
     *
     * This preserves bb_total, avoids negative P11, and is a validation closure
     * until a HydroLight/DPF-equivalent phase kernel is available.  Diagnostics
     * below still report the input a_total,b_total,bb_total. */
    int fixed_bulk_phase_model = o.fixed_bulk_iop_mode ? o.fixed_bulk_phase_model : 0;
    if (fixed_bulk_phase_model < 0 || fixed_bulk_phase_model > 4) fixed_bulk_phase_model = 0;

    fixed_bulk_phase_table_t fixed_phase_table;
    memset(&fixed_phase_table, 0, sizeof fixed_phase_table);
    rt_value_phase_interp_t fixed_phase_spline = {0};
    int have_fixed_phase_spline = 0;
    double fixed_bulk_delta_f = 0.0;
    double fixed_bulk_ff_betal[FIXED_BULK_BETAL_LMAX + 1];
    int    fixed_bulk_ff_betal_L = -1;   /* >=0 => use moment kernel for fixed-bulk (OCRT_FIXEDBULK_MOMENT) */
    for (int l = 0; l <= FIXED_BULK_BETAL_LMAX; ++l) fixed_bulk_ff_betal[l] = 0.0;

    fixed_bulk_phase_table_t ccrr_particle_phase_table;
    memset(&ccrr_particle_phase_table, 0, sizeof ccrr_particle_phase_table);
    rt_value_phase_interp_t ccrr_phase_spline = {0};
    int have_ccrr_phase_spline = 0;
    double ccrr_particle_delta_f = 0.0;
    int ccrr_particle_lut_mode = (constituent_model_ccrr &&
                                  o.ccrr_particle_phase_lut_path &&
                                  o.ccrr_particle_phase_lut_path[0]);
    /* P1a (v1.051): CCRR particle phase via internally-generated positive FF P11
     * with b-weighted component merge (pigment + mineral), replacing the
     * deprecated L=24 Legendre-moment path (Gibbs -> negative P11). */
    int ccrr_particle_ff_mode = (constituent_model_ccrr && !ccrr_particle_lut_mode &&
                                 !(o.water_mie_phase_path && o.water_mie_phase_path[0]) &&
                                 (b_pig + b_min) > 0.0);

    if (o.fixed_bulk_iop_mode && (fixed_bulk_phase_model == 2 || fixed_bulk_phase_model == 3)) {
        int lrc = fixed_bulk_phase_lut_load(o.fixed_bulk_phase_lut_path,
                                            o.fixed_bulk_phase_case_id,
                                            o.fixed_bulk_phase_wavelength_nm,
                                            &fixed_phase_table);
        if (lrc != 0) {
            fprintf(stderr, "[B3] failed to load fixed-bulk phase LUT '%s' case='%s' wl=%.10g rc=%d\n",
                    o.fixed_bulk_phase_lut_path ? o.fixed_bulk_phase_lut_path : "(null)",
                    o.fixed_bulk_phase_case_id ? o.fixed_bulk_phase_case_id : "",
                    o.fixed_bulk_phase_wavelength_nm, lrc);
            return -4;
        }
        if (fixed_bulk_phase_model == 3) {
            /* v1.06 PRODUCTION DEFAULT: OSOAA smooth forward cap + moment-GSF kernel.
             * Validated vs OSOAA/HL (3-way FF bb030 sweep, SZA=0): all omega0 within
             * the +/-3.2% OSOAA-HL band. Resolves the FF value-kernel non-convergence
             * (the old -17%) and the high-omega0 deficit. Physical bb preserved via
             * b* = b(1-A); no TMS (the fixed-bulk TMS-equivalent over-corrects this
             * bb-preserving path -- falsified, see OCRT_FIXEDBULK_TMS in DEBUG_FLAGS).
             * Debug overrides (all OCRT_DEBUG-gated, default-off):
             *   OCRT_FIXEDBULK_FULLPHASE : full phase, no delta-M (Gibbs diagnostic)
             *   OCRT_FIXEDBULK_MOMENT    : Wiscombe delta-M moment kernel (diagnostic)
             *   OCRT_FIXEDBULK_VALUE     : legacy hard-cut formal delta-M value kernel (pre-v1.06 default)
             *   OCRT_FIXEDBULK_CAPVALUE  : cap + value kernel (no moment; isolates kernel) */
            static int use_moment = -1, use_value = -1, use_full = -1, use_capvalue = -1;
            if (use_moment < 0) use_moment = (ocrt_debug_env("OCRT_FIXEDBULK_MOMENT") != NULL) ? 1 : 0;
            if (use_full < 0) use_full = (ocrt_debug_env("OCRT_FIXEDBULK_FULLPHASE") != NULL) ? 1 : 0;
            if (use_value < 0) use_value = (ocrt_debug_env("OCRT_FIXEDBULK_VALUE") != NULL) ? 1 : 0;
            if (use_capvalue < 0) use_capvalue = (ocrt_debug_env("OCRT_FIXEDBULK_CAPVALUE") != NULL) ? 1 : 0;
            if (use_full) {
                /* Fundamental-fix test: FULL phase via Legendre, NO delta-M removal
                 * (b unchanged, full omega0), like OSOAA's SOS_DECOMPO_LEGENDRE.
                 * betal[l] = (2l+1)*chi[l]. lmax up to FIXED_BULK_BETAL_LMAX=200.
                 * Needs n_mu >~ lmax. FF theta=0 already capped in the LUT. Env-gated. */
                int L = (o.fixed_bulk_lmax > 0 && o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX) ?
                        o.fixed_bulk_lmax : 30;
                double chi[FIXED_BULK_BETAL_LMAX + 1];
                if (fixed_bulk_phase_table_chi_moments(&fixed_phase_table, L, chi) != 0) {
                    fprintf(stderr, "[B3] fixed-bulk chi-moment (fullphase) failed L=%d\n", L);
                    fixed_bulk_phase_table_free(&fixed_phase_table);
                    return -4;
                }
                fixed_bulk_delta_f = 0.0;   /* NO forward removal => full omega0 */
                for (int l = 0; l <= L; ++l)
                    fixed_bulk_ff_betal[l] = (2.0 * (double)l + 1.0) * chi[l];
                fixed_bulk_ff_betal_L = L;
                { int di = g_wrt_env.f_OCRT_DUMP_IOP;
                  if (di) fprintf(stderr, "FULLPHASE(lut-deltam) L=%d f=0 (full omega0, no delta-M; high-L Legendre like OSOAA)\n", L); }
            } else if (use_moment) {
                /* Path A test: Wiscombe delta-M scaled Legendre moments consumed by
                 * the moment kernel (rt_kernel_phase_fourier), instead of the value
                 * kernel. lmax capped at FIXED_BULK_BETAL_LMAX=60. Backscatter may
                 * ring negative (Gibbs) for sharp FF. Env-gated; default unchanged. */
                int L = (o.fixed_bulk_lmax > 0 && o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX) ?
                        o.fixed_bulk_lmax : 30;
                double chi[FIXED_BULK_BETAL_LMAX + 1];
                if (fixed_bulk_phase_table_chi_moments(&fixed_phase_table, L, chi) != 0) {
                    fprintf(stderr, "[B3] fixed-bulk chi-moment computation failed L=%d\n", L);
                    fixed_bulk_phase_table_free(&fixed_phase_table);
                    return -4;
                }
                double f_dm = chi[L];
                if (f_dm < 0.0) f_dm = 0.0;
                if (f_dm > 1.0 - 1.0e-6) f_dm = 1.0 - 1.0e-6;
                fixed_bulk_delta_f = f_dm;
                double one_m_f = 1.0 - f_dm;
                for (int l = 0; l <= L; ++l)
                    fixed_bulk_ff_betal[l] = (2.0 * (double)l + 1.0) * (chi[l] - f_dm) / one_m_f;
                fixed_bulk_ff_betal_L = L;
                { int di = g_wrt_env.f_OCRT_DUMP_IOP;
                  if (di) fprintf(stderr, "MOMENT(lut-deltam) L=%d f_dm=%.6g (delta-M moment kernel; backscatter-Gibbs risk)\n",
                                  L, f_dm); }
            } else if (!use_value) {
                /* PRODUCTION DEFAULT: smooth log-linear forward cap
                 * (T1=31.79 T2=23.07 from inc/OSOAA.h cos 0.85/0.92, OCRT_CAP_T1/T2
                 * overridable) feeding the moment-GSF kernel below. Physical bb
                 * preserved via b* = b(1-A). Matches the CCRR production path. */
                double T1_cap = 6.0, T2_cap = 3.0;  /* CAP FIX 2026-06-17: weak near-delta (was 31.79/23.07) */
                { const char *e1=ocrt_debug_env("OCRT_CAP_T1"), *e2=ocrt_debug_env("OCRT_CAP_T2");
                  if (e1) T1_cap=atof(e1); if (e2) T2_cap=atof(e2); }
                double A_cap = 0.0;
                int crc = fixed_bulk_phase_table_loglinear_cap(&fixed_phase_table, T1_cap, T2_cap, &A_cap);
                if (crc != 0) {
                    fprintf(stderr, "[B3] fixed-bulk log-linear forward cap failed rc=%d\n", crc);
                    fixed_bulk_phase_table_free(&fixed_phase_table);
                    return -4;
                }
                fixed_bulk_delta_f = A_cap;
                /* Route the loglinear-capped, renormalized phase through the moment-GSF
                 * kernel (OSOAA's actual SOS). Capped phase is moment-representable
                 * (no Gibbs), so L up to FIXED_BULK_BETAL_LMAX=200 for backscatter
                 * convergence. OCRT_FIXEDBULK_CAPVALUE instead keeps the value kernel. */
                /* Phase-kernel selection (internal legacy selector, default
                 * 0=value).  Build delta-M moment betal (=> moment kernel) ONLY when
                 * the user explicitly selects moment (water_phase_kernel==1) or sets
                 * the legacy diagnostic env override.  With value selected (default)
                 * leave fixed_bulk_ff_betal_L=-1 so the Gibbs-free value kernel runs.
                 * OCRT_FIXEDBULK_CAPVALUE remains a diagnostic env that also forces value. */
                { if (!use_capvalue && o.water_phase_kernel == 1) {
                      int Lc = (o.fixed_bulk_lmax > 0 && o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX) ?
                               o.fixed_bulk_lmax : 200;
                      double chic[FIXED_BULK_BETAL_LMAX + 1];
                      if (fixed_bulk_phase_table_chi_moments(&fixed_phase_table, Lc, chic) != 0) {
                          fprintf(stderr, "[B3] cap+moment chi-moment failed L=%d\n", Lc);
                          fixed_bulk_phase_table_free(&fixed_phase_table);
                          return -4;
                      }
                      for (int l = 0; l <= Lc; ++l)
                          fixed_bulk_ff_betal[l] = (2.0 * (double)l + 1.0) * chic[l];
                      fixed_bulk_ff_betal_L = Lc;
                  } }
                { int di = g_wrt_env.f_OCRT_DUMP_IOP;
                  if (di) fprintf(stderr, "LOGCAP(lut-deltam) T1=%.4g T2=%.4g A=%.6g (smooth forward cap)\n",
                                  T1_cap, T2_cap, A_cap); }
            } else {
            int L = (o.fixed_bulk_lmax > 0 && o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX) ?
                    o.fixed_bulk_lmax : 6;
            double cut_deg = 0.0, resid_bb = 0.0;
            int frc = fixed_bulk_phase_table_formal_delta_m(&fixed_phase_table, L,
                                                            (b_tot > 0.0) ? (bb_tot / b_tot) : 0.0,
                                                            &fixed_bulk_delta_f,
                                                            &cut_deg,
                                                            &resid_bb);
            if (frc != 0) {
                fprintf(stderr, "[B3] fixed-bulk formal delta-M failed case='%s' wl=%.10g L=%d rc=%d\n",
                        o.fixed_bulk_phase_case_id ? o.fixed_bulk_phase_case_id : "",
                        o.fixed_bulk_phase_wavelength_nm, L, frc);
                fixed_bulk_phase_table_free(&fixed_phase_table);
                return -4;
            }
            { int di = g_wrt_env.f_OCRT_DUMP_IOP;
              if (di) fprintf(stderr, "DELTAM(lut-deltam) L=%d f=%.6g cut_deg=%.4g resid_bb(=bb/b of cut phase)=%.6g\n",
                              L, fixed_bulk_delta_f, cut_deg, resid_bb); }
            { /* session6 diag: dump the post-delta-M fixed_phase_table that the value
               * kernel actually integrates, to inspect P11 positivity/smoothness. */
              static int dumped_fb = 0; const char *dpf = g_wrt_env.s_OCRT_DUMP_FIXEDPHASE;
              if (dpf && !dumped_fb) { dumped_fb = 1; FILE *fp = fopen(dpf, "w");
                if (fp) { fprintf(fp, "theta_deg,P11\n");
                  for (int i = 0; i < fixed_phase_table.n; ++i)
                    fprintf(fp, "%.4f,%.8e\n", fixed_phase_table.theta_deg[i], fixed_phase_table.p11[i]);
                  fclose(fp); } } }
            }
        }
    }

    if (o.fixed_bulk_iop_mode && fixed_bulk_phase_model == 4) {
        double bbob = (b_tot > 0.0) ? (bb_tot / b_tot) : 0.0;
        if (!(bbob > 0.0 && bbob < 0.5)) return -4;
        double nrel = (o.fixed_bulk_ff_n > 1.0) ? o.fixed_bulk_ff_n : 1.18;
        double mu_j = (o.fixed_bulk_ff_mu > 0.0) ? o.fixed_bulk_ff_mu : ff_solve_mu_for_bb_cached(nrel, bbob);
        if (!(mu_j > 0.0) || !isfinite(mu_j)) {
            fprintf(stderr, "[B3] scalar FF: cannot solve mu for n=%.8g bb/b=%.8g\n", nrel, bbob);
            fixed_bulk_phase_table_free(&fixed_phase_table);
            return -4;
        }
        int grc = fixed_bulk_phase_table_generate_ff(nrel, mu_j, &fixed_phase_table);
        if (grc != 0) {
            fprintf(stderr, "[B3] scalar FF phase generation failed n=%.8g mu=%.8g rc=%d\n", nrel, mu_j, grc);
            fixed_bulk_phase_table_free(&fixed_phase_table);
            return -4;
        }
        int L = (o.fixed_bulk_lmax > 0 && o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX) ? o.fixed_bulk_lmax : 10;
        double cut_deg = 0.0, resid_bb = 0.0;
        int frc = fixed_bulk_phase_table_formal_delta_m(&fixed_phase_table, L, bbob,
                                                        &fixed_bulk_delta_f, &cut_deg, &resid_bb);
        if (frc != 0) {
            fprintf(stderr, "[B3] scalar FF formal delta-M failed L=%d bb/b=%.8g rc=%d\n", L, bbob, frc);
            fixed_bulk_phase_table_free(&fixed_phase_table);
            return -4;
        }
    }

    double b_particle_input = b_pig + b_min;
    const double b_tot_untruncated = b_tot;
    const double b_particle_input_untruncated = b_particle_input;
    double ocrt_particle_delta_f = 0.0;
    double ocrt_component_A_trunc[3] = {0.0, 0.0, 0.0};

    /* Stage 4: advanced OSOAA-style broad transform for OCRT constituent
     * Mie phases.
     *
     * CANONICAL POLICY: exact FR631 uses the raw phase with this branch OFF.
     * The reason is not that broad truncation is negligible: matched tests
     * showed that it can transform a large fraction of particle scattering and
     * move directional Rrs/rrs by several percent. A local 0--0.005-degree cap
     * was negligible, but that is a different operation from this broad lobe
     * transform.
     *
     * Each physical component is transformed independently, because a mixture
     * of raw phases does not in general have the same removed forward mass as
     * the b-weighted mixture of individually transformed residuals. The
     * returned residual phase is normalized and the matching transport
     * scattering coefficient is
     *
     *     b_eff = b_raw * (1 - A/2).
     *
     * Backscatter is unchanged: the transform modifies the forward hemisphere,
     * and normalization raises the residual bb/b by 1/(1-A/2). Absorption and
     * every bb component remain untouched. Finite-column comparisons must keep
     * the same pre-transform physical depth in OCRT and the reference model.
     *
     * DOC-REF:
     *   docs/OCRT_MIE_GRID_TRUNCATION_AND_VALIDATION_ARTIFACT_POLICY_2026-08-21.md */
    if (constituent_model_ocrt && !o.fixed_bulk_iop_mode &&
        b_particle_input > 0.0 &&
        (o.water_mie_truncation_mode != 0 || o.water_mie_ss_mode != 0)) {
        if (o.water_mie_ss_mode != 0) {
            fprintf(stderr,
                    "rt_water_rt: nonzero --ocrt-mie-ss-mode is not supported "
                    "by the direct constituent kernel; use ss-mode=0.\n");
            return -2;
        }
        if (o.water_mie_truncation_mode != 1 ||
            !water_particle_direct_enabled(1, 0)) {
            fprintf(stderr,
                    "rt_water_rt: --ocrt-mie-truncation requires the direct "
                    "theta-linear particle kernel.\n");
            return -2;
        }

        struct constituent_trunc_component {
            const char *label;
            const char *path;
            double *b_value;
        } tc[3] = {
            { "phyto",    organic_phyto_mie_path,    &b_phyto },
            { "detritus", organic_detritus_mie_path, &b_det },
            { "tsm",      ahn_mie_path_used,         &b_min }
        };

        double removed_b = 0.0;
        for (int ic = 0; ic < 3; ++ic) {
            const double b_raw = *tc[ic].b_value;
            if (!(b_raw > 0.0)) continue;
            if (!tc[ic].path || !tc[ic].path[0]) {
                fprintf(stderr,
                        "rt_water_rt: missing %s phase path for constituent "
                        "truncation at lambda=%.10g nm\n",
                        tc[ic].label, lambda_nm);
                return -2;
            }
            const rt_value_phase_interp_t *tr_phase = NULL;
            double A = 0.0;
            int phase_hit = 0;
            const int trc = water_component_value_phase_get(
                tc[ic].path, lambda_nm, 1, &tr_phase,
                NULL, NULL, &A, &phase_hit);
            if (trc != 0 || !tr_phase || !(A >= 0.0 && A < 2.0)) {
                fprintf(stderr,
                        "rt_water_rt: %s direct phase truncation failed "
                        "lambda=%.10g rc=%d A=%.12g\n",
                        tc[ic].label, lambda_nm, trc, A);
                return -2;
            }
            const double f = 0.5 * A;
            const double b_eff = b_raw * (1.0 - f);
            removed_b += b_raw - b_eff;
            *tc[ic].b_value = b_eff;
            ocrt_component_A_trunc[ic] = A;
            if (g_wrt_env.f_OCRT_DUMP_IOP) {
                fprintf(stderr,
                        "OCRT_TRUNC_COMPONENT name=%s A=%.12g f=%.12g "
                        "b_raw=%.12g b_eff=%.12g bb_b_residual=%.12g "
                        "cache=%s path=%s\n",
                        tc[ic].label, A, f, b_raw, b_eff,
                        tr_phase->bb_b_ratio,
                        phase_hit ? "hit" : "build", tc[ic].path);
            }
        }

        b_pig = comp.pigment.b + b_phyto + b_det;
        b_particle_input = b_pig + b_min;
        b_tot -= removed_b;
        if (b_tot < 0.0 && b_tot > -1.0e-12) b_tot = 0.0;
        if (!(b_tot >= 0.0) || !(b_particle_input >= 0.0)) {
            fprintf(stderr,
                    "rt_water_rt: invalid post-truncation scattering budget "
                    "b_tot=%.12g b_particle=%.12g\n",
                    b_tot, b_particle_input);
            return -2;
        }
        ocrt_particle_delta_f =
            (b_particle_input_untruncated > 0.0)
                ? removed_b / b_particle_input_untruncated : 0.0;
    }
    if (ccrr_particle_lut_mode && b_particle_input > 0.0) {
        int lrc = fixed_bulk_phase_lut_load(o.ccrr_particle_phase_lut_path,
                                            o.ccrr_particle_phase_case_id,
                                            o.ccrr_particle_phase_wavelength_nm,
                                            &ccrr_particle_phase_table);
        if (lrc != 0) {
            fprintf(stderr, "rt_water_rt: failed to load CCRR positive particle phase LUT '%s' case='%s' wl=%.10g rc=%d\n",
                    o.ccrr_particle_phase_lut_path ? o.ccrr_particle_phase_lut_path : "(null)",
                    o.ccrr_particle_phase_case_id ? o.ccrr_particle_phase_case_id : "",
                    o.ccrr_particle_phase_wavelength_nm, lrc);
            fixed_bulk_phase_table_free(&fixed_phase_table);
            return -4;
        }
        const double particle_bb_over_b = (b_particle_input > 0.0) ? ((bb_pig + bb_min) / b_particle_input) : 0.0;
        int L = (o.ccrr_particle_phase_lmax > 0 && o.ccrr_particle_phase_lmax <= FIXED_BULK_BETAL_LMAX) ?
                o.ccrr_particle_phase_lmax : 10;
        double cut_deg = 0.0, resid_bb = 0.0;
        int frc = fixed_bulk_phase_table_formal_delta_m(&ccrr_particle_phase_table, L, particle_bb_over_b,
                                                        &ccrr_particle_delta_f, &cut_deg, &resid_bb);
        if (frc != 0) {
            fprintf(stderr, "rt_water_rt: CCRR particle LUT formal delta-M failed case='%s' wl=%.10g L=%d bb/b=%.8g rc=%d\n",
                    o.ccrr_particle_phase_case_id ? o.ccrr_particle_phase_case_id : "",
                    o.ccrr_particle_phase_wavelength_nm, L, particle_bb_over_b, frc);
            fixed_bulk_phase_table_free(&fixed_phase_table);
            fixed_bulk_phase_table_free(&ccrr_particle_phase_table);
            return -4;
        }
    }

    /* v1.054: delta-M scaled Legendre moments of the FF particle phase, built in
     * the P1a block below and consumed by the betal/moment kernel (rt_kernel_phase_fourier),
     * replacing the hard-cut + direct phase-value Fourier path that caused the
     * L-dependent multiple-scattering oscillation. */
    double ccrr_ff_betal[FIXED_BULK_BETAL_LMAX + 1];
    int    ccrr_ff_betal_L = -1;
    int    ccrr_ff_value_kernel = 0;   /* 1 = FF consumed via value-based azimuth Fourier (path B); 0 = moment kernel */
    for (int l = 0; l <= FIXED_BULK_BETAL_LMAX; ++l) ccrr_ff_betal[l] = 0.0;

    /* P1a (v1.051): generate CCRR particle phase internally as b-weighted FF P11.
     * Components: pigment FF (mu solved from bb_pig/b_pig) and mineral FF
     * (mu solved from bb_min/b_min); merged on a common theta grid by O1 rule (a):
     *   P_mix(theta) = (b_pig*P_pig + b_min*P_min) / (b_pig + b_min).
     * Forward peak then removed by formal delta-M; phase fed to the kernel via the
     * same direct-Fourier path as the external-LUT mode. */
    if (ccrr_particle_ff_mode && b_particle_input > 0.0) {
        /* v1.053: pigment and mineral FF now use SEPARATE relative refractive
         * indices. Pigment uses the wavelength-dependent Cauchy model
         * (n_rel ~ 1.045, soft particle); mineral keeps n_rel = 1.18. */
        const double nrel_pig = pigment_ff_nrel(lambda_nm);
        const double nrel_min = (o.fixed_bulk_ff_n > 1.0) ? o.fixed_bulk_ff_n : 1.18;
        fixed_bulk_phase_table_t tab_pig, tab_min;
        memset(&tab_pig, 0, sizeof tab_pig);
        memset(&tab_min, 0, sizeof tab_min);
        int have_pig = 0, have_min = 0;
        if (b_pig > 0.0 && bb_pig > 0.0) {
            double bbob_pig = bb_pig / b_pig;
            double mu_pig = ff_solve_mu_for_bb_cached(nrel_pig, bbob_pig);
            if (mu_pig > 3.0 && isfinite(mu_pig) &&
                fixed_bulk_phase_table_generate_ff(nrel_pig, mu_pig, &tab_pig) == 0) have_pig = 1;
        }
        if (b_min > 0.0 && bb_min > 0.0) {
            double bbob_min = bb_min / b_min;
            double mu_min = ff_solve_mu_for_bb_cached(nrel_min, bbob_min);
            if (mu_min > 3.0 && isfinite(mu_min) &&
                fixed_bulk_phase_table_generate_ff(nrel_min, mu_min, &tab_min) == 0) have_min = 1;
        }
        const fixed_bulk_phase_table_t *ref = have_pig ? &tab_pig : (have_min ? &tab_min : NULL);
        if (!ref) {
            fprintf(stderr, "rt_water_rt: CCRR FF particle phase generation failed (b_pig=%.6g b_min=%.6g)\n", b_pig, b_min);
            fixed_bulk_phase_table_free(&tab_pig);
            fixed_bulk_phase_table_free(&tab_min);
            fixed_bulk_phase_table_free(&fixed_phase_table);
            return -4;
        }
        ccrr_particle_phase_table.n = ref->n;
        ccrr_particle_phase_table.theta_deg = (double*)calloc((size_t)ref->n, sizeof(double));
        ccrr_particle_phase_table.p11       = (double*)calloc((size_t)ref->n, sizeof(double));
        if (!ccrr_particle_phase_table.theta_deg || !ccrr_particle_phase_table.p11) {
            fixed_bulk_phase_table_free(&tab_pig);
            fixed_bulk_phase_table_free(&tab_min);
            fixed_bulk_phase_table_free(&ccrr_particle_phase_table);
            fixed_bulk_phase_table_free(&fixed_phase_table);
            return -4;
        }
        double wpig = have_pig ? b_pig : 0.0;
        double wmin = have_min ? b_min : 0.0;
        double wsum = wpig + wmin;
        for (int i = 0; i < ref->n; ++i) {
            ccrr_particle_phase_table.theta_deg[i] = ref->theta_deg[i];
            double ppig = have_pig ? tab_pig.p11[i] : 0.0;
            double pmin = have_min ? tab_min.p11[i] : 0.0;
            ccrr_particle_phase_table.p11[i] = (wpig * ppig + wmin * pmin) / wsum;
        }
        fixed_bulk_phase_table_free(&tab_pig);
        fixed_bulk_phase_table_free(&tab_min);
        double particle_bb_over_b = (bb_pig + bb_min) / b_particle_input;
        (void)particle_bb_over_b;
        int L = (o.ccrr_particle_phase_lmax > 0 && o.ccrr_particle_phase_lmax <= FIXED_BULK_BETAL_LMAX) ?
                o.ccrr_particle_phase_lmax : 10;
        /* Truncation method dispatch (o.water_phase_kernel):
         *   0 = VALUE kernel (path B, default): OSOAA-style angular forward-peak
         *       cap on the analytic merged FF VALUES, then the value-based
         *       azimuth Fourier kernel consumes P11 directly (no Legendre moment
         *       truncation -> no Gibbs -> positive radiance at any nphi).
         *   1 = MOMENT kernel (path A, SHELL/experimental): Wiscombe delta-M
         *       scaled moments -> rt_kernel_phase_fourier.  KNOWN-LIMITED: the
         *       delta-M residual rings negative at backscatter for low OS_NB
         *       (see OSOAA consultation).  TODO: replace with cap -> moments of
         *       the capped phase (OSOAA's actual moment-based recipe, positive
         *       at OS_NB>=64). Kept here only as a method placeholder. */
        if (o.water_phase_kernel == 1) {
            /* --- path A shell: Wiscombe delta-M moment kernel (legacy, negative-prone) --- */
            double chi[FIXED_BULK_BETAL_LMAX + 1];
            if (fixed_bulk_phase_table_chi_moments(&ccrr_particle_phase_table, L, chi) != 0) {
                fprintf(stderr, "rt_water_rt: CCRR FF chi-moment computation failed L=%d\n", L);
                fixed_bulk_phase_table_free(&ccrr_particle_phase_table);
                fixed_bulk_phase_table_free(&fixed_phase_table);
                return -4;
            }
            double f_dm = chi[L];
            if (f_dm < 0.0) f_dm = 0.0;
            if (f_dm > 1.0 - 1.0e-6) f_dm = 1.0 - 1.0e-6;
            ccrr_particle_delta_f = f_dm;
            double one_m_f = 1.0 - f_dm;
            for (int l = 0; l <= L; ++l)
                ccrr_ff_betal[l] = (2.0 * (double)l + 1.0) * (chi[l] - f_dm) / one_m_f;
            ccrr_ff_betal_L = L;
            ccrr_ff_value_kernel = 0;
            {
                static int warned_moment_shell = 0;
                if (!warned_moment_shell) {
                    warned_moment_shell = 1;
                    fprintf(stderr, "rt_water_rt: internal water phase kernel=moment (SHELL: Wiscombe delta-M, "
                                    "backscatter-negative for low OS_NB; TODO cap->moments)\n");
                }
            }
        } else {
            /* --- path B (default): OSOAA angular forward-peak cap + value kernel --- */
            (void)L;
            double A_cap = 0.0;
            /* OSOAA hydrosol forward-cap angles (inc/OSOAA.h: cos 0.85/0.92).
             * Overridable via env for diagnostics (OCRT_CAP_T1/T2). */
            double T1_cap = 6.0, T2_cap = 3.0;  /* CAP FIX 2026-06-17: weak near-delta (was 31.79/23.07) */
            { const char *e1=ocrt_debug_env("OCRT_CAP_T1"), *e2=ocrt_debug_env("OCRT_CAP_T2");
              if (e1) T1_cap=atof(e1); if (e2) T2_cap=atof(e2); }
            if (fixed_bulk_phase_table_loglinear_cap(&ccrr_particle_phase_table, T1_cap, T2_cap, &A_cap) != 0) {
                fprintf(stderr, "rt_water_rt: CCRR FF log-linear forward-peak cap failed\n");
                fixed_bulk_phase_table_free(&ccrr_particle_phase_table);
                fixed_bulk_phase_table_free(&fixed_phase_table);
                return -4;
            }
            ccrr_particle_delta_f = A_cap;   /* physical bb preserved via b* = b(1-A) */
            ccrr_ff_betal_L = -1;            /* signal: isotropic placeholder betal, value override below */
            ccrr_ff_value_kernel = 1;
            /* Diagnostic (OCRT_DUMP_PHASE=path): write the capped merged FF P11(theta)
             * once, to inspect the backscatter shape that drives the BRDF. */
            {
                static int dumped = 0;
                const char *dp = g_wrt_env.s_OCRT_DUMP_PHASE;
                if (dp && !dumped) {
                    dumped = 1;
                    FILE *fp = fopen(dp, "w");
                    if (fp) {
                        fprintf(fp, "# capped merged FF P11(theta) [sr-1]; A_removed=%.6f lambda=%.1f chl=%.4g min=%.4g\n",
                                A_cap, lambda_nm, o.ccrr_chl_mg_m3, o.ccrr_min_g_m3);
                        fprintf(fp, "theta_deg,P11\n");
                        for (int i = 0; i < ccrr_particle_phase_table.n; ++i)
                            fprintf(fp, "%.4f,%.8e\n", ccrr_particle_phase_table.theta_deg[i],
                                    ccrr_particle_phase_table.p11[i]);
                        fclose(fp);
                    }
                }
            }
        }
    }


    /* Pre-scan the hydrosol truncation coefficient for a native .mie path so
     * optical-depth bookkeeping b* = b(1-A/2) is known before ext/omega are
     * built.  In the direct path use the exact same theta-linear residual table
     * that will feed the Fourier kernel; the legacy moment diagnostic keeps its
     * historical Gauss/moment pre-scan. */
    if (o.fixed_bulk_iop_mode && o.water_mie_phase_path &&
        o.water_mie_phase_path[0] && o.water_mie_truncation_mode == 1) {
        if (water_particle_direct_enabled(1, o.water_mie_ss_mode)) {
            const rt_value_phase_interp_t *tr_phase = NULL;
            double Apre = 0.0;
            int phase_hit = 0;
            const int prc = water_component_value_phase_get(
                o.water_mie_phase_path, lambda_nm, 1, &tr_phase,
                NULL, NULL, &Apre, &phase_hit);
            if (prc != 0 || !tr_phase || !(Apre >= 0.0 && Apre < 2.0)) {
                fprintf(stderr,
                        "rt_water_rt: direct .mie truncation pre-scan failed "
                        "path='%s' lambda=%.10g rc=%d A=%.12g\n",
                        o.water_mie_phase_path, lambda_nm, prc, Apre);
                return -2;
            }
            fixed_bulk_delta_f = 0.5 * Apre;
            if (g_wrt_env.f_OCRT_DUMP_IOP)
                fprintf(stderr,
                        "MIE_TRUNC_PRE_DIRECT A=%.12g f=%.12g "
                        "bb_b_residual=%.12g cache=%s\n",
                        Apre, fixed_bulk_delta_f, tr_phase->bb_b_ratio,
                        phase_hit ? "hit" : "build");
        } else {
            mie_data_t mie_pre;
            aerosol_phase_interp_t aip_pre;
            memset(&aip_pre, 0, sizeof aip_pre);
            if (read_mie_file(o.water_mie_phase_path, &mie_pre) != 0 ||
                build_aerosol_interpolators(
                    &mie_pre, lambda_nm * 1.0e-3, &aip_pre) != 0) {
                fprintf(stderr,
                        "rt_water_rt: failed to pre-scan .mie truncation '%s'\n",
                        o.water_mie_phase_path);
                return -2;
            }
            const int Na_pre = mie_pre.n_ang;
            double *P11p = (double *)malloc((size_t)Na_pre * sizeof(double));
            double *P12p = (double *)malloc((size_t)Na_pre * sizeof(double));
            double *P33p = (double *)malloc((size_t)Na_pre * sizeof(double));
            double *thp  = (double *)malloc((size_t)Na_pre * sizeof(double));
            double tmpb[3], tmpg[3], tmpa[3], tmpz[3];
            if (!P11p || !P12p || !P33p || !thp) {
                free(P11p); free(P12p); free(P33p); free(thp);
                aerosol_phase_interp_free(&aip_pre);
                mie_data_free(&mie_pre);
                return -4;
            }
            for (int ii = 0; ii < Na_pre; ++ii) {
                thp[ii] = mie_pre.angles[ii];
                eval_aerosol_phase(
                    &aip_pre, thp[ii], &P11p[ii], &P12p[ii], &P33p[ii]);
            }
            double Apre = 0.0;
            const int nmgpre = (o.water_mie_moment_n_mu > 1)
                                   ? o.water_mie_moment_n_mu : 60;
            const int prc = rt_aerosol_compute_vector_legendre_gauss_truncated(
                P11p, P12p, P33p, thp, Na_pre, 2, nmgpre,
                WATER_HYD_TRUNC_MU1, WATER_HYD_TRUNC_MU2,
                WATER_HYD_TRUNC_THRESHOLD,
                tmpb, tmpg, tmpa, tmpz, &Apre);
            if (prc != 0) {
                fprintf(stderr,
                        "rt_water_rt: .mie truncation pre-scan failed rc=%d\n",
                        prc);
                free(P11p); free(P12p); free(P33p); free(thp);
                aerosol_phase_interp_free(&aip_pre);
                mie_data_free(&mie_pre);
                return -4;
            }
            fixed_bulk_delta_f = 0.5 * Apre;
            if (g_wrt_env.f_OCRT_DUMP_IOP)
                fprintf(stderr,
                        "MIE_TRUNC_PRE_MOMENT A=%.12g f=%.12g nmg=%d\n",
                        Apre, fixed_bulk_delta_f, nmgpre);
            free(P11p); free(P12p); free(P33p); free(thp);
            aerosol_phase_interp_free(&aip_pre);
            mie_data_free(&mie_pre);
        }
    }

    double b_rt = b_tot;
    if (o.fixed_bulk_iop_mode) {
        if (o.water_mie_phase_path && o.water_mie_phase_path[0]) {
            /* #0 vector .mie path.  Default keeps the full phase with full b.
             * The truncation mode pairs the truncated phase
             * moments with b* = b(1-A/2), unless OCRT_WATER_MIE_TRUNC_FULLB=1
             * is set to isolate moment-shape from optical-depth bookkeeping. */
            const char *tr_fullb = g_wrt_env.s_OCRT_WATER_MIE_TRUNC_FULLB;
            if (o.water_mie_truncation_mode == 1 && !(tr_fullb && tr_fullb[0] && strcmp(tr_fullb,"0") != 0))
                b_rt = b_tot * (1.0 - fixed_bulk_delta_f);
            else
                b_rt = b_tot;
        } else if (fixed_bulk_phase_model == 0) {
            b_rt = 2.0 * bb_tot;
        } else if (fixed_bulk_phase_model == 3 || fixed_bulk_phase_model == 4) {
            /* session6 MS-method test (Nakajima-Tanaka 1988 sec 2.3): use the
             * truncated phase for the diffuse field but keep the FULL b (and hence
             * full omega/tau) so that high-order multiple scattering is not
             * over-damped.  OCRT_FIXEDBULK_MSFULLOMEGA=1 enables; default keeps
             * the standard delta-M scaling b*=b(1-f). */
            { const char *msfo = g_wrt_env.s_OCRT_FIXEDBULK_MSFULLOMEGA;
              b_rt = (msfo && msfo[0] && strcmp(msfo,"0") != 0)
                         ? b_tot : b_tot * (1.0 - fixed_bulk_delta_f); }
        } else {
            b_rt = b_tot;
        }
    } else if ((ccrr_particle_lut_mode || ccrr_particle_ff_mode) && b_particle_input > 0.0) {
        b_rt = b_w + b_particle_input * (1.0 - ccrr_particle_delta_f);
    }
    double ext   = a_tot + b_rt;
    double omega = (ext > 0.0) ? (b_rt / ext) : 0.0;
    double omega_w = (ext > 0.0) ? (b_w / ext) : 0.0;
    double b_particle = o.fixed_bulk_iop_mode ? b_rt : ((ccrr_particle_lut_mode || ccrr_particle_ff_mode) ? b_particle_input * (1.0 - ccrr_particle_delta_f) : b_particle_input);
    double omega_particle = (ext > 0.0) ? (b_particle / ext) : 0.0;
    const int ocrt_organic_vector_mix_mode =
        constituent_model_ocrt && o.ccrr_chl_mg_m3 > 0.0 &&
        b_particle_input > 0.0;
    { int di = g_wrt_env.f_OCRT_DUMP_IOP;
      if (di) {
          const double delta_report = o.fixed_bulk_iop_mode
                                          ? fixed_bulk_delta_f
                                          : (constituent_model_ocrt
                                                 ? ocrt_particle_delta_f
                                                 : ccrr_particle_delta_f);
          fprintf(stderr,
                  "IOP a=%.6g b_raw=%.6g b=%.6g bb=%.6g delta_f=%.6g "
                  "b_rt=%.6g omega_in=%.6g omega_eff=%.6g\n",
                  a_tot, b_tot_untruncated, b_tot, bb_tot, delta_report, b_rt,
                  (a_tot + b_tot_untruncated > 0.0
                       ? b_tot_untruncated / (a_tot + b_tot_untruncated) : 0.0),
                  omega);
          if (constituent_model_ocrt && o.water_mie_truncation_mode == 1) {
              fprintf(stderr,
                      "OCRT_TRUNC_SUMMARY particle_b_raw=%.12g "
                      "particle_b_eff=%.12g delta_f=%.12g "
                      "A_phyto=%.12g A_detritus=%.12g A_tsm=%.12g\n",
                      b_particle_input_untruncated, b_particle_input,
                      ocrt_particle_delta_f,
                      ocrt_component_A_trunc[0],
                      ocrt_component_A_trunc[1],
                      ocrt_component_A_trunc[2]);
          }
          if (constituent_model_ocrt) {
              fprintf(stderr,
                      "OCRT_IOP water[a=%.8e b=%.8e bb=%.8e] cdom[a=%.8e] "
                      "phyto[a=%.8e b=%.8e bb=%.8e] detritus[a=%.8e b=%.8e bb=%.8e] "
                      "tsm[a=%.8e b=%.8e bb=%.8e]\n",
                      a_w, b_w, bb_w, a_cdom,
                      a_phyto, b_phyto, bb_phyto,
                      a_det, b_det, bb_det,
                      a_min, b_min, bb_min);
          }
      } }

    if (constituent_model_ccrr && o.ccrr_phase_moments_path && o.ccrr_phase_moments_path[0]) {
        if (rt_iop_ccrr_phase_moments_load(o.ccrr_phase_moments_path) != 0) {
            fprintf(stderr, "rt_water_rt: failed to load CCRR phase moments '%s'\n",
                    o.ccrr_phase_moments_path);
            return -2;
        }
    }

    enum { CCRR_PARTICLE_LMAX = CCRR_PARTICLE_BETAL_LMAX };
    double particle_betal[FIXED_BULK_BETAL_LMAX + 1];
    double pig_betal[CCRR_PARTICLE_LMAX + 1];
    double min_betal[CCRR_PARTICLE_LMAX + 1];
    for (int l = 0; l <= FIXED_BULK_BETAL_LMAX; ++l) {
        particle_betal[l] = 0.0;
    }
    for (int l = 0; l <= CCRR_PARTICLE_LMAX; ++l) {
        pig_betal[l] = 0.0;
        min_betal[l] = 0.0;
    }
    int particle_Lmax = -1;
    /* #0: particle vector-phase Legendre moments (P12→gammal, P22/P33→alphal/zetal).
     * NULL = scalar particle (current FF/CCRR scalar paths).  Set by the .mie
     * vector-phase path (step 2) via rt_aerosol_compute_vector_legendre. */
    const double *particle_gammal = NULL;
    const double *particle_alphal = NULL;
    const double *particle_zetal  = NULL;
    double mie_gammal_buf[FIXED_BULK_BETAL_LMAX + 1];
    double mie_alphal_buf[FIXED_BULK_BETAL_LMAX + 1];
    double mie_zetal_buf [FIXED_BULK_BETAL_LMAX + 1];
    /* Step13 diagnostic NT-TMS parity: OSOAA HYD_TRUNCATION=ON behaves like
     * raw PM for the direct/order-1 source and truncated PM for the MS source.
     * Keep raw moments alongside the truncated moments so primary_source can be
     * built with raw phase while the SOS transport/source loop uses truncated
     * phase + b_eff/tau_tr bookkeeping. Diagnostic only; production default OFF. */
    double mie_raw_betal_buf[FIXED_BULK_BETAL_LMAX + 1];
    double mie_raw_gammal_buf[FIXED_BULK_BETAL_LMAX + 1];
    double mie_raw_alphal_buf[FIXED_BULK_BETAL_LMAX + 1];
    double mie_raw_zetal_buf [FIXED_BULK_BETAL_LMAX + 1];
    int have_mie_raw_moments = 0;
    for (int l = 0; l <= FIXED_BULK_BETAL_LMAX; ++l) {
        mie_raw_betal_buf[l] = 0.0;
        mie_raw_gammal_buf[l] = 0.0;
        mie_raw_alphal_buf[l] = 0.0;
        mie_raw_zetal_buf[l] = 0.0;
    }
    const int particle_direct_mode = water_particle_direct_enabled(
        o.water_mie_truncation_mode, o.water_mie_ss_mode);
    const rt_value_phase_interp_t *particle_value_phase = NULL;
    char particle_value_phase_key[512];
    particle_value_phase_key[0] = '\0';
    int have_particle_value_phase = 0;
    if (ocrt_organic_vector_mix_mode) {
        if (o.water_mie_ss_mode != 0) {
            fprintf(stderr,
                    "rt_water_rt: nonzero --ocrt-mie-ss-mode is not supported "
                    "by the direct constituent kernel.\n");
            return -2;
        }
        if (o.water_mie_truncation_mode == 1 && !particle_direct_mode) {
            fprintf(stderr,
                    "rt_water_rt: constituent real-angle truncation requires "
                    "the direct theta-linear particle kernel.\n");
            return -2;
        }

        if (particle_direct_mode) {
            struct phase_component_direct {
                const char *label;
                const char *path;
                double b;
                double bb;
            } components[3] = {
                { "phyto", organic_phyto_mie_path, b_phyto, bb_phyto },
                { "detritus", organic_detritus_mie_path, b_det, bb_det },
                { "tsm", ahn_mie_path_used, b_min, bb_min }
            };
            const char *mix_paths[3] = {
                components[0].path, components[1].path, components[2].path
            };
            const double mix_b[3] = {
                components[0].b, components[1].b, components[2].b
            };
            int n_component = 0, mix_cache_hit = 0;
            unsigned long long signature = 0ULL;
            const int vrc = water_mixed_value_phase_get(
                mix_paths, mix_b, lambda_nm,
                o.water_mie_truncation_mode, &particle_value_phase,
                &n_component, &signature, &mix_cache_hit);
            if (vrc != 0 || !particle_value_phase) {
                fprintf(stderr,
                        "rt_water_rt: direct OCRT constituent phase build failed "
                        "lambda=%.10g rc=%d\n", lambda_nm, vrc);
                return -3;
            }
            have_particle_value_phase = 1;
            snprintf(particle_value_phase_key, sizeof particle_value_phase_key,
                     "MIX:%016llx", (unsigned long long)signature);

            /* These coefficients are only an allocation contract for the
             * legacy workspace.  Every particle Fourier table is overwritten
             * by the direct vector kernel before primary-source/SOS use. */
            int Ldirect = o.m_max_water;
            if (Ldirect < 2) Ldirect = 2;
            if (Ldirect > FIXED_BULK_BETAL_LMAX)
                Ldirect = FIXED_BULK_BETAL_LMAX;
            for (int l = 0; l <= Ldirect; ++l) {
                particle_betal[l] = 0.0;
                mie_gammal_buf[l] = 0.0;
                mie_alphal_buf[l] = 0.0;
                mie_zetal_buf[l] = 0.0;
            }
            particle_betal[0] = 1.0;
            if (Ldirect >= 1)
                particle_betal[1] = 3.0 * particle_value_phase->g_asym;
            particle_Lmax = Ldirect;
            particle_gammal = mie_gammal_buf;
            particle_alphal = mie_alphal_buf;
            particle_zetal = mie_zetal_buf;

            if (g_wrt_env.f_OCRT_DUMP_IOP) {
                for (int ic = 0; ic < 3; ++ic) {
                    if (!(components[ic].b > 0.0)) continue;
                    const rt_value_phase_interp_t *cp = NULL;
                    double A_component = 0.0;
                    int ch = 0;
                    if (water_component_value_phase_get(
                            components[ic].path, lambda_nm,
                            o.water_mie_truncation_mode, &cp,
                            NULL, NULL, &A_component, &ch) == 0 && cp) {
                        double p11_90 = 0.0, p12_90 = 0.0, p33_90 = 0.0;
                        rt_value_phase_interp_eval_theta(
                            cp, 90.0, &p11_90, &p12_90, &p33_90);
                        fprintf(stderr,
                                "OCRT_PHASE_COMPONENT name=%s weight_b=%.12g "
                                "b=%.12g bb_b_iop=%.12g bb_b_phase=%.12g "
                                "A=%.12g P90=[%.8e %.8e %.8e] "
                                "cache=%s path=%s\n",
                                components[ic].label,
                                components[ic].b / b_particle_input,
                                components[ic].b,
                                components[ic].bb / components[ic].b,
                                cp->bb_b_ratio,
                                A_component,
                                p11_90, p12_90, p33_90,
                                ch ? "hit" : "build", components[ic].path);
                    }
                }
                fprintf(stderr,
                        "OCRT_PHASE_MIX_DIRECT components=%d n=%d fr631=%d "
                        "norm0=%.12g g=%.12g bb_b=%.12g cache=%s key=%s\n",
                        n_component, particle_value_phase->n,
                        particle_value_phase->is_fr631,
                        particle_value_phase->norm_before,
                        particle_value_phase->g_asym,
                        particle_value_phase->bb_b_ratio,
                        mix_cache_hit ? "hit" : "build",
                        particle_value_phase_key);
            }
        } else {
        const int Lmix = (o.fixed_bulk_lmax > 0 &&
                          o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX)
                             ? o.fixed_bulk_lmax : 200;
        const int mix_mode = o.water_mie_moment_mode;
        const int mix_nmg = (o.water_mie_moment_n_mu > 1)
                                ? o.water_mie_moment_n_mu : 60;
        double tmp_b[FIXED_BULK_BETAL_LMAX + 1];
        double tmp_g[FIXED_BULK_BETAL_LMAX + 1];
        double tmp_a[FIXED_BULK_BETAL_LMAX + 1];
        double tmp_z[FIXED_BULK_BETAL_LMAX + 1];

        struct phase_component {
            const char *label;
            const char *path;
            double b;
        } components[3] = {
            { "phyto", organic_phyto_mie_path, b_phyto },
            { "detritus", organic_detritus_mie_path, b_det },
            { "tsm", ahn_mie_path_used, b_min }
        };
        const char *mix_paths[3] = {
            components[0].path, components[1].path, components[2].path
        };
        const double mix_b[3] = {
            components[0].b, components[1].b, components[2].b
        };
        int n_component = 0;
        double recon_min = 0.0;
        const int mixed_cache_hit = water_mixed_phase_cache_load(
            mix_paths, mix_b, lambda_nm, mix_mode, mix_nmg, Lmix,
            particle_betal, mie_gammal_buf, mie_alphal_buf, mie_zetal_buf,
            &n_component, &recon_min);

        if (!mixed_cache_hit) {
            for (int l = 0; l <= Lmix; ++l) {
                particle_betal[l] = 0.0;
                mie_gammal_buf[l] = 0.0;
                mie_alphal_buf[l] = 0.0;
                mie_zetal_buf[l] = 0.0;
            }
            for (int ic = 0; ic < 3; ++ic) {
                if (!(components[ic].b > 0.0)) continue;
                if (!components[ic].path || !components[ic].path[0]) {
                    fprintf(stderr,
                            "rt_water_rt: missing %s vector phase path for positive b=%.10g\n",
                            components[ic].label, components[ic].b);
                    return -2;
                }
                double p11_90 = 0.0, p12_90 = 0.0, p33_90 = 0.0, bbob = 0.0;
                int cache_hit = 0;
                int mrc = water_component_phase_moments(
                    components[ic].path, lambda_nm, mix_mode, mix_nmg, Lmix,
                    tmp_b, tmp_g, tmp_a, tmp_z,
                    &p11_90, &p12_90, &p33_90, &bbob, &cache_hit);
                if (mrc != 0) {
                    fprintf(stderr,
                            "rt_water_rt: %s vector phase decomposition failed path='%s' lambda=%.10g rc=%d\n",
                            components[ic].label, components[ic].path, lambda_nm, mrc);
                    return -3;
                }
                const double weight = components[ic].b / b_particle_input;
                for (int l = 0; l <= Lmix; ++l) {
                    particle_betal[l] += weight * tmp_b[l];
                    mie_gammal_buf[l] += weight * tmp_g[l];
                    mie_alphal_buf[l] += weight * tmp_a[l];
                    mie_zetal_buf[l] += weight * tmp_z[l];
                }
                ++n_component;
                { int di = g_wrt_env.f_OCRT_DUMP_IOP;
                  if (di) fprintf(stderr,
                                  "OCRT_PHASE_COMPONENT name=%s weight_b=%.12g b=%.12g bb_b_iop=%.12g bb_b_phase=%.12g P90=[%.8e %.8e %.8e] cache=%s path=%s\n",
                                  components[ic].label, weight, components[ic].b,
                                  (components[ic].b > 0.0) ?
                                      ((ic == 0 ? bb_phyto : (ic == 1 ? bb_det : bb_min)) /
                                       components[ic].b) : 0.0,
                                  bbob, p11_90, p12_90, p33_90,
                                  cache_hit ? "hit" : "build", components[ic].path); }
            }
            if (n_component < 1) {
                /* With phytoplankton scattering disabled, detritus alone is the
                 * expected Chl-linked particulate phase. */
                fprintf(stderr,
                        "rt_water_rt: OCRT organic vector mixture has no scattering component\n");
                return -3;
            }
            recon_min = water_phase_reconstruction_min_back(particle_betal, Lmix);
            if (water_mixed_phase_cache_store(
                    mix_paths, mix_b, lambda_nm, mix_mode, mix_nmg, Lmix,
                    particle_betal, mie_gammal_buf, mie_alphal_buf, mie_zetal_buf,
                    n_component, recon_min) != 0) {
                fprintf(stderr, "rt_water_rt: failed to cache mixed OCRT organic phase\n");
                return -3;
            }
            if (o.mu_sun_water_override <= 0.0 && recon_min < 0.0) {
                fprintf(stderr,
                        "rt_water_rt: *** S-009 WARNING *** mixed OCRT organic phase: L=%d Legendre reconstruction is negative at sampled back angles (min=%.3e).\n"
                        "  P11/P12/P33 are all included and b-weighted, but a validated multi-component forward-truncation policy is still required for this phase.\n",
                        Lmix, recon_min);
            }
        }
        particle_Lmax = Lmix;
        particle_gammal = mie_gammal_buf;
        particle_alphal = mie_alphal_buf;
        particle_zetal = mie_zetal_buf;

        { int di = g_wrt_env.f_OCRT_DUMP_IOP;
          if (di) fprintf(stderr,
                          "OCRT_PHASE_MIX components=%d L=%d beta0=%.12g beta2=%.12g gamma2=%.12g alpha2=%.12g zeta2=%.12g recon_back_min=%.12g cache=%s\n",
                          n_component, Lmix, particle_betal[0],
                          (Lmix >= 2) ? particle_betal[2] : 0.0,
                          (Lmix >= 2) ? mie_gammal_buf[2] : 0.0,
                          (Lmix >= 2) ? mie_alphal_buf[2] : 0.0,
                          (Lmix >= 2) ? mie_zetal_buf[2] : 0.0,
                          recon_min, mixed_cache_hit ? "hit" : "build"); }
        }
    } else if (o.fixed_bulk_iop_mode) {
        if (fixed_bulk_phase_model == 1 || fixed_bulk_phase_model == 2 || fixed_bulk_phase_model == 3 || fixed_bulk_phase_model == 4) {
            int L = (o.fixed_bulk_lmax > 0 && o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX) ? o.fixed_bulk_lmax : 30;
            if (fixed_bulk_ff_betal_L >= 0) {
                /* moment kernel (OCRT_FIXEDBULK_MOMENT): real delta-M scaled FF moments */
                for (int l = 0; l <= fixed_bulk_ff_betal_L; ++l) particle_betal[l] = fixed_bulk_ff_betal[l];
                particle_Lmax = fixed_bulk_ff_betal_L;
            } else {
                particle_betal[0] = 1.0;
                for (int l = 1; l <= L; ++l) particle_betal[l] = 0.0;
                particle_Lmax = L;
            }
        } else {
            particle_betal[0] = 1.0;
            particle_betal[1] = 0.0;
            particle_betal[2] = 0.0;
            particle_Lmax = 2;
        }
        omega_w = 0.0;
        omega_particle = (ext > 0.0) ? (b_particle / ext) : 0.0;
    } else if (b_particle > 0.0) {
        if (ccrr_particle_lut_mode || ccrr_particle_ff_mode) {
            int L = (o.ccrr_particle_phase_lmax > 0 && o.ccrr_particle_phase_lmax <= FIXED_BULK_BETAL_LMAX) ?
                    o.ccrr_particle_phase_lmax : 10;
            if (ccrr_particle_ff_mode && ccrr_ff_betal_L >= 0) {
                /* v1.054: real delta-M scaled FF moments -> moment-based (Gibbs-free)
                 * phase_fourier via rt_kernel_phase_fourier (no direct phase-value override). */
                for (int l = 0; l <= ccrr_ff_betal_L; ++l) particle_betal[l] = ccrr_ff_betal[l];
                particle_Lmax = ccrr_ff_betal_L;
            } else {
                particle_betal[0] = 1.0;
                for (int l = 1; l <= L; ++l) particle_betal[l] = 0.0;
                particle_Lmax = L;
            }
        } else {
            int Lpig = rt_iop_ccrr_particle_betal(0, CCRR_PARTICLE_LMAX + 1, pig_betal);
            int Lmin = rt_iop_ccrr_particle_betal(1, CCRR_PARTICLE_LMAX + 1, min_betal);
            if (Lpig < 0 || Lmin < 0) return -4;
            particle_Lmax = (Lpig > Lmin) ? Lpig : Lmin;
            if (particle_Lmax > CCRR_PARTICLE_LMAX) particle_Lmax = CCRR_PARTICLE_LMAX;
            for (int l = 0; l <= particle_Lmax; ++l) {
                particle_betal[l] = (b_pig * pig_betal[l] + b_min * min_betal[l]) / b_particle_input;
            }
        }
    }
    /* #0 step2: .mie VECTOR particle phase override.  When an OCRT/IOP .mie
     * phase source is selected, replace the scalar particle phase (betal only)
     * with the full Mueller
     * Legendre moments (betal/gammal/alphal/zetal) computed from the .mie
     * P11/P12/P33 at this wavelength.  Bulk IOP (a,b,bb,omega_particle) stays as
     * set above (from the direct IOP branch).  Reuses the atmospheric aerosol
     * machinery (mie_io PCHIP + rt_aerosol_compute_vector_legendre) for
     * convention consistency with the in-water SOS aerosol-slot kernels. */
    /* v1.06 pol-SS: keep the particle phase interpolator alive past this block
     * so the polarized single-scatter exact-phase correction can evaluate
     * P11/P12/P33 at the view scattering angle. */
    aerosol_phase_interp_t aip; int have_aip = 0; memset(&aip, 0, sizeof aip);
    double *mie_value_kernel_allm = NULL;
    double *mie_value_beamq_allm = NULL;
    int mie_value_kernel_owned = 0;
    int mie_value_kernel_mcount = 0;
    /* v1.09 #24: key for the per-(wavelength,.mie) phase-decomposition cache. */
    int wmc_moment_mode = o.water_mie_moment_mode;
    int wmc_nmg = (o.water_mie_moment_n_mu > 1) ? o.water_mie_moment_n_mu : 60;
    int wmc_lmax_req = (o.fixed_bulk_lmax > 0 && o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX)
                       ? o.fixed_bulk_lmax : 200;
    int wmc_use = (o.water_mie_phase_path && o.water_mie_phase_path[0] && omega_particle > 0.0);
    const int wmc_direct = wmc_use && particle_direct_mode;
    int wmc_hit = (!wmc_direct) && wmc_use &&
                  wmc_key_match(o.water_mie_phase_path, lambda_nm,
                    o.water_mie_truncation_mode, o.water_mie_ss_mode,
                    wmc_moment_mode, wmc_nmg, wmc_lmax_req);
    if (getenv("OCRT_WMC_DEBUG")) fprintf(stderr,
            "[WMC] use=%d direct=%d hit=%d path=%s lam=%.1f valid=%d cachedlam=%.1f\n",
            wmc_use, wmc_direct, wmc_hit,
            o.water_mie_phase_path?o.water_mie_phase_path:"(null)",
            lambda_nm, g_wmc.valid, g_wmc.lambda_nm);
    if (wmc_direct) {
        double phase_ssa = 0.0, phase_g_bulk = 0.0;
        double phase_A_trunc = 0.0;
        int phase_cache_hit = 0;
        const int vrc = water_component_value_phase_get(
            o.water_mie_phase_path, lambda_nm,
            o.water_mie_truncation_mode, &particle_value_phase,
            &phase_ssa, &phase_g_bulk, &phase_A_trunc, &phase_cache_hit);
        if (vrc != 0 || !particle_value_phase) {
            fprintf(stderr,
                    "rt_water_rt: direct water .mie phase build failed path='%s' "
                    "lambda=%.10g rc=%d\n",
                    o.water_mie_phase_path, lambda_nm, vrc);
            return -2;
        }
        have_particle_value_phase = 1;
        unsigned long long h = 1469598103934665603ULL;
        h = fnv1a_bytes(h, o.water_mie_phase_path,
                        strlen(o.water_mie_phase_path) + 1u);
        h = fnv1a_bytes(h, &lambda_nm, sizeof lambda_nm);
        h = fnv1a_bytes(h, &o.water_mie_truncation_mode,
                        sizeof o.water_mie_truncation_mode);
        snprintf(particle_value_phase_key, sizeof particle_value_phase_key,
                 "FILE:%016llx", (unsigned long long)h);

        int Ldirect = o.m_max_water;
        if (Ldirect < 2) Ldirect = 2;
        if (Ldirect > FIXED_BULK_BETAL_LMAX)
            Ldirect = FIXED_BULK_BETAL_LMAX;
        for (int l = 0; l <= Ldirect; ++l) {
            particle_betal[l] = 0.0;
            mie_gammal_buf[l] = 0.0;
            mie_alphal_buf[l] = 0.0;
            mie_zetal_buf[l] = 0.0;
        }
        particle_betal[0] = 1.0;
        if (Ldirect >= 1)
            particle_betal[1] = 3.0 * particle_value_phase->g_asym;
        particle_Lmax = Ldirect;
        particle_gammal = mie_gammal_buf;
        particle_alphal = mie_alphal_buf;
        particle_zetal = mie_zetal_buf;
        if (!constituent_model_ocrt)
            fixed_bulk_delta_f = 0.5 * phase_A_trunc;

        if (g_wrt_env.f_OCRT_DUMP_IOP) {
            double p11_90 = 0.0, p12_90 = 0.0, p33_90 = 0.0;
            rt_value_phase_interp_eval_theta(
                particle_value_phase, 90.0, &p11_90, &p12_90, &p33_90);
            fprintf(stderr,
                    "MIEPHASE_DIRECT '%s' n=%d fr631=%d lambda_nm=%.3f "
                    "ssa=%.6g g_bulk=%.6g A=%.12g norm0=%.12g g_phase=%.12g "
                    "bb_b=%.12g P90=[%.8e %.8e %.8e] cache=%s\n",
                    o.water_mie_phase_path, particle_value_phase->n,
                    particle_value_phase->is_fr631, lambda_nm, phase_ssa,
                    phase_g_bulk, phase_A_trunc,
                    particle_value_phase->norm_before,
                    particle_value_phase->g_asym,
                    particle_value_phase->bb_b_ratio,
                    p11_90, p12_90, p33_90,
                    phase_cache_hit ? "hit" : "build");
        }
    } else if (wmc_hit) {
        /* Restore decomposition + interpolator from cache; skip disk read and
         * the L<=200 Legendre expansion entirely (the whole point of #24). */
        int Lc = g_wmc.particle_Lmax;
        memcpy(particle_betal,  g_wmc.betal,  (size_t)(Lc + 1) * sizeof(double));
        memcpy(mie_gammal_buf,  g_wmc.gammal, (size_t)(Lc + 1) * sizeof(double));
        memcpy(mie_alphal_buf,  g_wmc.alphal, (size_t)(Lc + 1) * sizeof(double));
        memcpy(mie_zetal_buf,   g_wmc.zetal,  (size_t)(Lc + 1) * sizeof(double));
        particle_Lmax   = Lc;
        particle_gammal = mie_gammal_buf;
        particle_alphal = mie_alphal_buf;
        particle_zetal  = mie_zetal_buf;
        fixed_bulk_delta_f = g_wmc.delta_f;
        aip = g_wmc.aip;       /* shallow copy: cache owns the buffers */
        have_aip = 1;          /* freed by wmc, NOT by the cleanup path */
    } else
    if (wmc_use) {
        mie_data_t mie;
        if (read_mie_file(o.water_mie_phase_path, &mie) != 0) {
            fprintf(stderr, "rt_water_rt: failed to read .mie '%s'\n", o.water_mie_phase_path);
            return -2;
        }
        if (build_aerosol_interpolators(&mie, lambda_nm * 1.0e-3, &aip) != 0) {
            fprintf(stderr, "rt_water_rt: build_aerosol_interpolators FAIL (lambda=%.1f nm)\n", lambda_nm);
            mie_data_free(&mie); return -2;
        }
        int Na = mie.n_ang;
        double *P11t = (double*)malloc((size_t)Na * sizeof(double));
        double *P12t = (double*)malloc((size_t)Na * sizeof(double));
        double *P33t = (double*)malloc((size_t)Na * sizeof(double));
        double *tht  = (double*)malloc((size_t)Na * sizeof(double));
        if (!P11t || !P12t || !P33t || !tht) {
            free(P11t); free(P12t); free(P33t); free(tht);
            aerosol_phase_interp_free(&aip); mie_data_free(&mie); return -4;
        }
        for (int i = 0; i < Na; ++i) {
            tht[i] = mie.angles[i];
            eval_aerosol_phase(&aip, tht[i], &P11t[i], &P12t[i], &P33t[i]);
        }
        int Lm = (o.fixed_bulk_lmax > 0 && o.fixed_bulk_lmax <= FIXED_BULK_BETAL_LMAX)
                 ? o.fixed_bulk_lmax : 200;   /* #23 (S-009): water-mie default was 30 —
                                              * sediment phases need L~120-200 (measured). */
        int vrc = 0;
        if (o.water_mie_moment_mode == 1) {
            int nmg = (o.water_mie_moment_n_mu > 1) ? o.water_mie_moment_n_mu : 60;
            if (o.water_mie_truncation_mode == 1 && o.water_mie_ss_mode == 2) {
                int rrc = rt_aerosol_compute_vector_legendre_gauss(
                        P11t, P12t, P33t, tht, Na, Lm, nmg,
                        mie_raw_betal_buf, mie_raw_gammal_buf, mie_raw_alphal_buf, mie_raw_zetal_buf);
                if (rrc != 0) {
                    fprintf(stderr, "rt_water_rt: water .mie raw moment generation for NT-TMS FAIL rc=%d\n", rrc);
                    free(P11t); free(P12t); free(P33t); free(tht);
                    aerosol_phase_interp_free(&aip); mie_data_free(&mie); return -3;
                }
                have_mie_raw_moments = 1;
            }
            if (o.water_mie_truncation_mode == 1) {
                double A_tr = 0.0;
                vrc = rt_aerosol_compute_vector_legendre_gauss_truncated(
                        P11t, P12t, P33t, tht, Na, Lm, nmg,
                        0.85, 0.92, 0.1,
                        particle_betal, mie_gammal_buf, mie_alphal_buf, mie_zetal_buf, &A_tr);
                fixed_bulk_delta_f = 0.5 * A_tr;
            } else {
                vrc = rt_aerosol_compute_vector_legendre_gauss(
                        P11t, P12t, P33t, tht, Na, Lm, nmg,
                        particle_betal, mie_gammal_buf, mie_alphal_buf, mie_zetal_buf);
            }
        } else {
            vrc = rt_aerosol_compute_vector_legendre(P11t, P12t, P33t, tht, Na,
                                                     Lm, 4096,
                                                     particle_betal, mie_gammal_buf,
                                                     mie_alphal_buf, mie_zetal_buf);
        }
        if (vrc != 0) {
            fprintf(stderr, "rt_water_rt: water .mie vector moment generation FAIL rc=%d mode=%d\n",
                    vrc, o.water_mie_moment_mode);
            free(P11t); free(P12t); free(P33t); free(tht);
            aerosol_phase_interp_free(&aip); mie_data_free(&mie); return -3;
        }
        /* v1.09 #23 (S-009 self-diagnosis): rrs is backscatter-dominated, and a
         * strongly forward-peaked phase (g >~ 0.97) is NOT representable by an
         * L<=200 Legendre expansion — the reconstructed phase goes NEGATIVE at
         * back angles (measured: blend g=0.981 min ~ -1.15 at L=30, still
         * negative at L=400; sediment g~0.94 converges by L~120-200).  Check
         * the reconstruction sign at five back angles and warn loudly. */
        {
            const double mus_bk[5] = { -1.0, -0.93969262078591, -0.86602540378444,
                                       -0.70710678118655, -0.5 };
            double recon_min = 1e300;
            for (int t_i = 0; t_i < 5; ++t_i) {
                double x = mus_bk[t_i], pm1 = 1.0, p0 = x;
                double sum = particle_betal[0] + ((Lm >= 1) ? particle_betal[1] * x : 0.0);
                for (int l = 2; l <= Lm; ++l) {
                    double pl = ((2.0 * l - 1.0) * x * p0 - (l - 1.0) * pm1) / (double)l;
                    sum += particle_betal[l] * pl;
                    pm1 = p0; p0 = pl;
                }
                if (sum < recon_min) recon_min = sum;
            }
            if (recon_min < 0.0) {
                fprintf(stderr,
                    "rt_water_rt: *** S-009 WARNING *** water .mie phase '%s': L=%d Legendre\n"
                    "  expansion reconstructs NEGATIVE at back angles (min=%.3e).  rrs from the\n"
                    "  water-mie expansion path is unreliable for this phase.  Use a less\n"
                    "  forward-peaked phase or branch-matched delta truncation\n"
                    "  (--ocrt-mie-truncation/--ocrt-mie-ss-mode or\n"
                    "   --iop-mie-truncation/--iop-mie-ss-mode), and validate\n"
                    "  against an expansion-free path.\n",
                    o.water_mie_phase_path, Lm, recon_min);
            }
        }

        const char *dump_mom = g_wrt_env.s_OCRT_DUMP_WATER_MIE_MOMENTS;
        if (dump_mom && dump_mom[0]) {
            FILE *df = fopen(dump_mom, "w");
            if (df) {
                fprintf(df, "k,alpha,beta11,gamma12,zeta\n");
                for (int kk = 0; kk <= Lm; ++kk) {
                    fprintf(df, "%d,%.17g,%.17g,%.17g,%.17g\n",
                            kk, mie_alphal_buf[kk], particle_betal[kk],
                            mie_gammal_buf[kk], mie_zetal_buf[kk]);
                }
                fclose(df);
            }
        }
        particle_Lmax   = Lm;
        particle_gammal = mie_gammal_buf;
        particle_alphal = mie_alphal_buf;
        particle_zetal  = mie_zetal_buf;
        { int dv = g_wrt_env.f_OCRT_DUMP_IOP;
          if (dv) {
            double p11_90, p12_90, p33_90; eval_aerosol_phase(&aip, 90.0, &p11_90, &p12_90, &p33_90);
            fprintf(stderr,
            "MIEPHASE '%s' Na=%d n_phase_wl=%d lambda_nm=%.3f lambda_um=%.5f ssa=%.4f Lm=%d moment_mode=%d nmg=%d\n"
            "  ang[0]=%.2f ang[Na-1]=%.2f  eval@90: P11=%.4e P12=%.4e P33=%.4e\n"
            "  betal0=%.4f betal2=%.4f gammal2=%.5f alphal2=%.5f zetal2=%.5f bb_b(mie)=%.5f omega_p=%.4f\n",
            o.water_mie_phase_path, Na, mie.n_phase_wl, lambda_nm, lambda_nm*1.0e-3, aip.ssa, Lm,
            o.water_mie_moment_mode, o.water_mie_moment_n_mu,
            mie.angles[0], mie.angles[Na-1], p11_90, p12_90, p33_90,
            particle_betal[0], particle_betal[2], mie_gammal_buf[2],
            mie_alphal_buf[2], mie_zetal_buf[2], aip.bb_b_ratio, omega_particle); } }
        free(P11t); free(P12t); free(P33t); free(tht);
        /* v1.09 #24: hand mie+aip+moments to the threadprivate cache so the
         * remaining grid cells at this wavelength skip the read+expansion.
         * The cache owns aip (and mie) from here; do NOT free them locally. */
        wmc_free();
        g_wmc.mie = mie;  g_wmc.have_mie = 1;
        g_wmc.aip = aip;  g_wmc.have_aip = 1;
        {
            int Lc = Lm;
            memcpy(g_wmc.betal,  particle_betal,  (size_t)(Lc + 1) * sizeof(double));
            memcpy(g_wmc.gammal, mie_gammal_buf,  (size_t)(Lc + 1) * sizeof(double));
            memcpy(g_wmc.alphal, mie_alphal_buf,  (size_t)(Lc + 1) * sizeof(double));
            memcpy(g_wmc.zetal,  mie_zetal_buf,   (size_t)(Lc + 1) * sizeof(double));
            g_wmc.particle_Lmax = Lc;
            g_wmc.delta_f = fixed_bulk_delta_f;
            snprintf(g_wmc.path, sizeof g_wmc.path, "%s", o.water_mie_phase_path);
            g_wmc.lambda_nm = lambda_nm;
            g_wmc.trunc = o.water_mie_truncation_mode; g_wmc.ss = o.water_mie_ss_mode;
            g_wmc.moment_mode = wmc_moment_mode; g_wmc.nmg = wmc_nmg;
            g_wmc.lmax_req = wmc_lmax_req; g_wmc.valid = 1;
        }
        have_aip = 1;   /* aip lives in the cache now (freed via wmc_free) */
    }
    /* A single constituent (most commonly TSM-only) uses the native .mie
     * path rather than the organic mixture block.  Emit the same structured
     * diagnostic so external phase reconstruction does not require adding a
     * dummy second constituent.  Debug-only: no production work or output. */
    if (g_wrt_env.f_OCRT_DUMP_IOP && constituent_model_ocrt &&
        !ocrt_organic_vector_mix_mode && wmc_use && have_aip && b_particle_input > 0.0) {
        double p11_90=0.0,p12_90=0.0,p33_90=0.0;
        eval_aerosol_phase(&aip,90.0,&p11_90,&p12_90,&p33_90);
        const char *label=(b_min>0.0 && b_pig<=0.0)?"tsm":"particle";
        const double bb_iop=bb_pig+bb_min;
        const double rmin=water_phase_reconstruction_min_back(particle_betal,particle_Lmax);
        fprintf(stderr,
            "OCRT_PHASE_COMPONENT name=%s weight_b=1 b=%.12g bb_b_iop=%.12g bb_b_phase=%.12g P90=[%.8e %.8e %.8e] cache=%s path=%s\n",
            label,b_particle_input,bb_iop/b_particle_input,aip.bb_b_ratio,
            p11_90,p12_90,p33_90,wmc_hit?"hit":"build",o.water_mie_phase_path);
        fprintf(stderr,
            "OCRT_PHASE_MIX components=1 L=%d beta0=%.12g beta2=%.12g gamma2=%.12g alpha2=%.12g zeta2=%.12g recon_back_min=%.12g cache=%s\n",
            particle_Lmax,particle_betal[0],particle_Lmax>=2?particle_betal[2]:0.0,
            particle_Lmax>=2?mie_gammal_buf[2]:0.0,
            particle_Lmax>=2?mie_alphal_buf[2]:0.0,
            particle_Lmax>=2?mie_zetal_buf[2]:0.0,rmin,wmc_hit?"hit":"build");
    }
    /* 2. Snell refraction */
    double mu_sun_air  = cos(sza_deg_air * M_PI / 180.0);
    double mu_view_air = cos(vza_deg_air * M_PI / 180.0);
    if (mu_sun_air <= 0.0) return -1;
    /* [FIX-SKY-EDLU 2026-06-28] sub-cone equivalent-beam: take the in-water solar
     * cosine directly when the caller provides it (in-water dir below mu_crit has
     * no real air angle).  Else standard air->water Snell refraction. */
    double mu_sun_water  = (o.mu_sun_water_override > 0.0)
                             ? o.mu_sun_water_override
                             : snell_down(mu_sun_air,  n_water);
    /* 1B (OCRT_VZA_IN_WATER): treat --vza as the IN-WATER view zenith angle
     * directly (mu_view_water = cos(vza)), instead of an in-air angle refracted
     * down. This lets the TMS-corrected view reconstruction sample any in-water
     * direction including super-critical (theta_w>theta_c, TIR region) for the
     * rrs(0-) BRDF. For the above-water (0+) transform we refract back up; a
     * super-critical mu_view_water has NO above-water counterpart (TIR), flagged
     * by mu_view_air=0 and zeroed in the T_wa block below. Env-gated, no physics
     * change to the default path; for HL BRDF comparison. */
    double mu_view_water;
    { int vza_in_water = g_wrt_env.f_OCRT_VZA_IN_WATER;
      if (vza_in_water) {
          mu_view_water = mu_view_air;   /* vza_deg_air reinterpreted as in-water */
          double n_crit = sqrt(1.0 - 1.0/(n_water*n_water));
          double sin_w  = sqrt(fmax(0.0, 1.0 - mu_view_water*mu_view_water));
          double sin_a  = n_water * sin_w;
          mu_view_air = (mu_view_water > n_crit && sin_a < 1.0) ? sqrt(fmax(0.0,1.0 - sin_a*sin_a)) : 0.0;
      } else {
          if (mu_view_air <= 0.0) return -1;
          mu_view_water = snell_down(mu_view_air, n_water);
      }
    }

    /* 3. τ_max
     *   determine_tau_max 의 K_factor 는 약 5 (5 e-foldings ≈ 99.3% 흡수) 로
     *   "충분히 깊은" 의미만 가지며, a_w 와 a_tot (=a_w+a_CDOM) 어느 쪽으로
     *   계산해도 z_max 가 추가 깊어질 뿐 결과 누락은 없다. 그러나 a_CDOM 이
     *   클 경우 ext 변화로 z_max 단축이 가능하므로 a_tot 를 인자로 넘긴다.
     *   (현재 determine_tau_max 본문은 인자를 무시하고 상수 처리하므로
     *    실질 영향은 없으나, 명시성을 위해 a_tot 전달.) */
    double tau_base = determine_tau_max(a_tot, b_w, o.tau_max_target);
    if (!(tau_base > 0.0)) tau_base = 20.0;

    /* Adaptive deep-water depth policy.  Large τ/depth settings are safety
     * caps, not mandatory work.  Start from a base total optical depth, extend
     * only as far as needed for the estimated transport attenuation to fall
     * below depth_bottom_tol, then cap by max_tau_max_target and max_z_max_m.
     *
     * This keeps high-TSM cases from being truncated at τ=15 while avoiding the
     * previous brute-force rule that always forced the full 200 m physical
     * depth and therefore tens of thousands of layers in red/NIR high-TSM
     * cases.  SOS order loops still terminate early on their own convergence
     * criterion. */
    double g_particle = 0.0;
    if (particle_Lmax >= 1) {
        g_particle = particle_betal[1] / 3.0;
        if (g_particle > 0.999) g_particle = 0.999;
        if (g_particle < -0.999) g_particle = -0.999;
    }
    double phase_b_sum = b_particle + b_w;
    double g_eff = (phase_b_sum > 0.0) ? (b_particle * g_particle) / phase_b_sum : 0.0;
    if (g_eff > 0.999) g_eff = 0.999;
    if (g_eff < -0.999) g_eff = -0.999;

    double z_req = (ext > 0.0) ? (tau_base / ext) : 0.0;
    double tol_depth = (o.depth_bottom_tol > 0.0) ? o.depth_bottom_tol : 1.0e-8;
    if (tol_depth < 1.0e-30) tol_depth = 1.0e-30;
    if (tol_depth > 1.0e-2) tol_depth = 1.0e-2;
    double efolds = -log(tol_depth);
    double transport_ext = a_tot + b_rt * fmax(0.0, 1.0 - g_eff);
    if (transport_ext > 0.0) {
        double z_tr = efolds / transport_ext;
        if (z_tr > z_req) z_req = z_tr;
    }
    double tau_max = (ext > 0.0) ? ext * z_req : tau_base;
    if (o.max_tau_max_target > 0.0 && tau_max > o.max_tau_max_target) {
        tau_max = o.max_tau_max_target;
        z_req = (ext > 0.0) ? tau_max / ext : z_req;
    }
    if (o.max_z_max_m > 0.0 && z_req > o.max_z_max_m) {
        z_req = o.max_z_max_m;
        tau_max = (ext > 0.0) ? ext * z_req : tau_max;
    }
    double z_max = (ext > 0.0) ? (tau_max / ext) : z_req;

    const double tau_max_full_for_trunc_diag = tau_max;
    const double ext_full_for_trunc_diag = a_tot + b_tot;
    const double z_full_for_trunc_diag = (ext_full_for_trunc_diag > 0.0) ? (tau_max_full_for_trunc_diag / ext_full_for_trunc_diag) : 0.0;
    (void)z_full_for_trunc_diag;

    if (g_wrt_env.s_OCRT_DUMP_IOP != NULL) {
        fprintf(stderr, "WATER_DEPTH_AUTO tau_base=%.9g tau_max=%.9g z_max=%.9g ext=%.9g transport_ext=%.9g g_eff=%.9g max_tau=%.9g max_z=%.9g tol=%.3g\n",
                tau_base, tau_max, z_max, ext, transport_ext, g_eff,
                o.max_tau_max_target, o.max_z_max_m, tol_depth);
    }
    int n_layers_water_eff = o.n_layers_water;
    if (o.layer_dtau_target > 0.0 && tau_max > 0.0) {
        int n_need = (int)ceil(tau_max / o.layer_dtau_target);
        if (n_need > n_layers_water_eff) n_layers_water_eff = n_need;
    }
    if (n_layers_water_eff < 1) n_layers_water_eff = 1;
    if (g_wrt_env.s_OCRT_DUMP_IOP != NULL) {
        fprintf(stderr, "WATER_LAYER_AUTO n_layers=%d tau_max=%.9g dtau=%.9g\n",
                n_layers_water_eff, tau_max, tau_max / (double)n_layers_water_eff);
    }

    /* OSOAA-compatible hydrosol truncation bookkeeping for fixed-bulk diagnostics.
     * OSOAA builds the physical sea profile from the original bulk extinction
     * (a + b), then replaces the transport/scattering budget by the truncated
     * coefficient b* = b(1-f).  Therefore the effective transport optical
     * depth is tau_tr = z_original * (a + b*), not the original tau_max
     * imposed again after reducing b.  The previous OCRT diagnostic path kept
     * tau_max fixed after b->b*, which changed the physical water depth and
     * did not match OSOAA's _TR bookkeeping.
     */
    if (o.fixed_bulk_iop_mode && o.water_mie_phase_path && o.water_mie_phase_path[0] &&
        o.water_mie_truncation_mode == 1) {
        const char *tr_fullb = g_wrt_env.s_OCRT_WATER_MIE_TRUNC_FULLB;
        const int keep_full_b = (tr_fullb && tr_fullb[0] && strcmp(tr_fullb,"0") != 0);
        if (!keep_full_b && b_tot > 0.0) {
            const double ext_orig = a_tot + b_tot;
            if (ext_orig > 0.0 && ext > 0.0) {
                const double z_orig = tau_max_full_for_trunc_diag / ext_orig;
                const double tau_tr = z_orig * ext;
                double tau_selected = tau_tr;
                const char *tau_mode = g_wrt_env.s_OCRT_WATER_MIE_TRUNC_TRANSPORT_TAU;
                if (tau_mode && tau_mode[0]) {
                    if (strcmp(tau_mode, "full") == 0) {
                        tau_selected = tau_max_full_for_trunc_diag;
                    } else if (strncmp(tau_mode, "mix:", 4) == 0) {
                        double alpha = atof(tau_mode + 4);
                        if (alpha < 0.0) alpha = 0.0;
                        if (alpha > 1.0) alpha = 1.0;
                        tau_selected = tau_tr + alpha * (tau_max_full_for_trunc_diag - tau_tr);
                    } else if (strcmp(tau_mode, "tr") == 0 || strcmp(tau_mode, "trunc") == 0) {
                        tau_selected = tau_tr;
                    }
                }
                tau_max = tau_selected;
                z_max = tau_max / ext;
                if (g_wrt_env.s_OCRT_DUMP_IOP != NULL) {
                    fprintf(stderr, "OSOAA_TRUNC_TAU_MODE mode=%s ext_orig=%.9g ext_tr=%.9g z_orig=%.9g tau_full=%.9g tau_tr=%.9g tau_used=%.9g z_used=%.9g\n",
                            (tau_mode && tau_mode[0]) ? tau_mode : "tr", ext_orig, ext, z_orig,
                            tau_max_full_for_trunc_diag, tau_tr, tau_max, z_max);
                }
            }
        }
    }

    /* 4. Build rt_atm_t for in-water medium.
     *
     * view_as_node is an output/exact-view diagnostic policy: add the requested
     * water-side view direction as a zero-weight node without removing any
     * quadrature node.  Weight zero means it never participates in angular
     * integration, but the source term into that exact direction is evaluated
     * by the same production kernels.  This is diagnostic only; CCRR/IOCCG21
     * comparison defaults to continuous reconstruction. */
    int n_mu_quad = o.n_mu_water;
    int add_target_node = 0;
    if (o.view_as_node && o.n_mu_water >= 2) {
        double mu_probe[256], w_probe[256];
        if (o.n_mu_water > 256) return -4;
        if (rt_quadrature_gauss_legendre_pos(o.n_mu_water, mu_probe, w_probe) != 0) return -4;
        (void)mu_probe; (void)w_probe;
        add_target_node = 1;
    }
    const int add_zero_slot_nodes = ocrt_water_zero_slot_solve_mode() ? 3 : 0;
    /* v1.09 #20: multi view-node direct extraction (S-007 root fix). */
    const int add_view_list = (o.n_view_vza > 0 && o.view_vza_deg_list &&
                               o.n_mu_water >= 2) ? o.n_view_vza : 0;
    int n_mu_alloc = o.n_mu_water + add_target_node + add_view_list + add_zero_slot_nodes;

    rt_atm_t atm;
    memset(&atm, 0, sizeof atm);
    if (rt_atm_alloc(&atm, n_layers_water_eff, n_mu_alloc) != 0) return -4;
    int brc = build_inwater_atm(&atm, n_layers_water_eff, n_mu_alloc,
                                 tau_max, omega_w, omega_particle,
                                 mu_sun_water, n_mu_quad,
                                 particle_Lmax, (particle_Lmax >= 0) ? particle_betal : NULL,
                                 (particle_Lmax >= 0) ? particle_gammal : NULL,
                                 (particle_Lmax >= 0) ? particle_alphal : NULL,
                                 (particle_Lmax >= 0) ? particle_zetal  : NULL);
    if (brc != 0) {
        rt_atm_free(&atm);
        return -4;
    }
    if (add_target_node || add_view_list || add_zero_slot_nodes) {
        /* v1.10 B-0a.2: node-ring construction unified through the
         * OSOAA-parity angle table (rt_angles_unified).  GL(2N)-half
         * core + sorted-insert/dedup specials - the same multiset the
         * previous append-then-sort produced, so outputs stay
         * BIT-identical absent angle coincidences; on a coincidence the
         * table dedups (OSOAA behaviour), which the old path did not. */
        rt_uangles_t uang;
        if (rt_uangles_init(&uang, n_mu_quad) != 0) {
            rt_atm_free(&atm);
            return -4;
        }
        int uidx;
        if (add_target_node)
            (void)rt_uangles_add(&uang, mu_view_water, &uidx);
        if (add_view_list) {
            for (int ivn = 0; ivn < add_view_list; ++ivn) {
                double mu_wv = snell_down(cos(o.view_vza_deg_list[ivn] * M_PI / 180.0),
                                          n_water);
                (void)rt_uangles_add(&uang, mu_wv, &uidx);
            }
        }
        if (add_zero_slot_nodes) {
            (void)rt_uangles_add(&uang, 1.0,          &uidx);
            (void)rt_uangles_add(&uang, mu_sun_water, &uidx);
            (void)rt_uangles_add(&uang, mu_sun_air,   &uidx);
        }
        atm.n_mu = uang.n_total;   /* may shrink if the table deduped */
        for (int jU = 1; jU <= uang.n_total; ++jU) {
            atm.rm[+jU] = +uang.mu[jU-1];
            atm.rm[-jU] = -uang.mu[jU-1];
            atm.gb[+jU] =  uang.w [jU-1];
            atm.gb[-jU] =  uang.w [jU-1];
        }
    }
        if (getenv("OCRT_S6_TRACE")) {
            fprintf(stderr, "[S6R] nmw=%d quad=%d n=%d :", o.n_mu_water, n_mu_quad, atm.n_mu);
            for (int jR = 1; jR <= atm.n_mu; ++jR) fprintf(stderr, " %.12f", atm.rm[+jR]);
            fprintf(stderr, "\n");
        }

    /* 5. Workspace allocation
     * 짚어둘 점: atm.L_max=0 (aerosol slot OFF) 이지만 Rayleigh-like phase는
     * l=2 까지 사용하므로 workspace의 l_max는 *최소 2*로 강제. atmospheric
     * pure-Rayleigh path도 동일 (Rayleigh has β_2, γ_2, α_2 at l=2). */
    rt_legendre_workspace_t ws;
    memset(&ws, 0, sizeof ws);
    int ws_l_max = (atm.L_max >= 2) ? atm.L_max : 2;
    if (rt_legendre_workspace_alloc(&ws, atm.n_mu, ws_l_max) != 0) {
        rt_atm_free(&atm);
        return -4;
    }
    const int water_nt_tms = (o.fixed_bulk_iop_mode && o.water_mie_phase_path && o.water_mie_phase_path[0]
                              && o.water_mie_truncation_mode == 1 && o.water_mie_ss_mode == 2
                              && have_mie_raw_moments && particle_Lmax >= 0);
    rt_legendre_workspace_t ws_primary;
    memset(&ws_primary, 0, sizeof ws_primary);
    if (water_nt_tms && rt_legendre_workspace_alloc(&ws_primary, atm.n_mu, ws_l_max) != 0) {
        rt_legendre_workspace_free(&ws);
        rt_atm_free(&atm);
        return -4;
    }

    /* 6. Allocate Stokes field buffers */
    const int nt   = atm.n_layers;
    const int n_mu = atm.n_mu;
    const int dirs = 2 * n_mu + 1;

    if (g_wrt_env.f_OCRT_VALUE_PHASE_SPLINE) {
        const int nmg = (o.water_mie_moment_n_mu > 1) ? o.water_mie_moment_n_mu : 400;
        if (fixed_phase_table.n >= 3 && rt_value_phase_interp_build(
                &fixed_phase_spline, fixed_phase_table.theta_deg,
                fixed_phase_table.p11, NULL, NULL, fixed_phase_table.n, nmg) == 0)
            have_fixed_phase_spline = 1;
        if (ccrr_particle_phase_table.n >= 3 && rt_value_phase_interp_build(
                &ccrr_phase_spline, ccrr_particle_phase_table.theta_deg,
                ccrr_particle_phase_table.p11, NULL, NULL,
                ccrr_particle_phase_table.n, nmg) == 0)
            have_ccrr_phase_spline = 1;
    }


    const size_t field_sz = (size_t)(nt + 1) * (size_t)dirs;
    /* v1.09 commit #21 (parallel perf): ONE arena allocation replaces the 22
     * per-entry calloc()s of this function plus the 3 per-mode prim buffers.
     * Under OpenMP batch (--batch-full-grid) this function runs per (vza,raa)
     * cell; on Windows/MinGW the CRT heap lock serialized threads and even
     * added convoy overhead (user report: OMP=4 slower than OMP=1).  Slice
     * pointers keep every downstream reference unchanged; calloc semantics for
     * the per-mode buffers are preserved with explicit memset at m-loop entry. */
    int M = o.m_max_water + 1;
    size_t _w_off = 0, _w_total =
        9*(size_t)field_sz + 6*(size_t)M + 6*(size_t)M*(size_t)n_mu + (size_t)dirs;
    double *wbuf = calloc(_w_total, sizeof(double));
    #define WSLICE(n) (wbuf ? (wbuf + (_w_off += (n)) - (n)) : NULL)
    double *src_i = WSLICE((size_t)field_sz);
    double *src_q = WSLICE((size_t)field_sz);
    double *src_u = WSLICE((size_t)field_sz);
    double *tot_i = WSLICE((size_t)field_sz);
    double *tot_q = WSLICE((size_t)field_sz);
    double *tot_u = WSLICE((size_t)field_sz);
    double *prim_i = WSLICE((size_t)field_sz);   /* #21: hoisted from m-loop */
    double *prim_q = WSLICE((size_t)field_sz);
    double *prim_u = WSLICE((size_t)field_sz);
    /* Per-m at view direction storage */
    double *I_m_view = WSLICE((size_t)M);
    double *Iss_m_view = WSLICE((size_t)M); /* v1.054 TMS: single-scatter (primary order) at view direction */
    double *Qss_m_view = WSLICE((size_t)M); /* v1.06 pol-SS: single-scatter Q at view direction */
    double *Uss_m_view = WSLICE((size_t)M); /* v1.06 pol-SS: single-scatter U at view direction */
    double *Q_m_view = WSLICE((size_t)M);
    double *U_m_view = WSLICE((size_t)M);
    /* Per-m at all positive μ for hemispheric integral */
    double *I_m_pos  = WSLICE((size_t)M * (size_t)n_mu);
    double *I_m_neg  = WSLICE((size_t)M * (size_t)n_mu);
    /* m=0 Q at upward positive nodes is needed for water-side internal
     * Fresnel reflection contribution to Ed(0-). */
    double *Q_m_pos  = WSLICE((size_t)M * (size_t)n_mu);
    /* For Kd, Ku: store m=0 mode at level k=1 too */
    double *I_m0_lvl1 = WSLICE((size_t)dirs);
    /* Diagnostic bottom-flux audit (m=0 hemispheric quantities at k=nt).
     * These are π-normalized until the common f_scale conversion below. */
    double bottom_Ed_diffuse_norm = 0.0;
    double bottom_Eu_diffuse_norm = 0.0;
    /* Cox-Munk T_wa (option A): per-m upwelling-node radiance field at z=0-,
     * for ALL m (not just m=0), so the water->air BTDF can be integrated over
     * the full in-water angular field. Indexed [m*n_mu + (jp-1)]. */
    double *I_m_node = WSLICE((size_t)M * (size_t)n_mu);
    double *Q_m_node = WSLICE((size_t)M * (size_t)n_mu);
    double *U_m_node = WSLICE((size_t)M * (size_t)n_mu);
    #undef WSLICE

    int rc = 0;
    int max_orders_seen = 0;
    int all_conv = 1;
    double worst_resid = 0.0;
    double *fb_allm_kernel = NULL;   /* v1.09 #8: per-beam all-m fixed-bulk kernel cache */
    int     fb_allm_mcount = 0; (void)fb_allm_mcount;

    if (!wbuf) {
        rc = -4;
        goto cleanup;
    }

    /* 7. m loop
     * If CCRR mode has no particle scattering, L_max remains 2 even though
     * the validation preset may request m_max=24.  Modes with m>L_max are
     * exactly zero for the current phase expansion; leave their buffers at
     * zero instead of asking the Legendre workspace to compute impossible
     * associated polynomials. */
    const int m_loop_max = (o.m_max_water < ws_l_max) ? o.m_max_water : ws_l_max;
    if (have_particle_value_phase && particle_value_phase) {
        const char *np = getenv("OCRT_WATER_VALUE_NPHI");
        int nphi = np ? atoi(np) : 0;
        if (nphi < 16) {
            const double g = fabs(particle_value_phase->g_asym);
            nphi = (g < 0.90) ? 360 : ((g < 0.95) ? 720 : ((g < 0.98) ? 2880 : 5760));
        }
        mie_value_kernel_mcount = m_loop_max + 1;
        const size_t plane = (size_t)(n_mu + 1) * (size_t)dirs;
        const size_t all = (size_t)mie_value_kernel_mcount * plane;

        const int cache_off = g_wrt_env.f_OCRT_NO_PHASE_CACHE;
        water_pol_value_kernel_cache_t *kc = NULL;
        int phase_cache_hit = 0;
        int beam_cache_hit = 0;
        const char *beam_state = "build";
        if (!cache_off) {
            kc = wpvkc_find(particle_value_phase_key, lambda_nm,
                            o.water_mie_truncation_mode, o.water_mie_ss_mode,
                            wmc_moment_mode, wmc_nmg, wmc_lmax_req,
                            &atm, mie_value_kernel_mcount, nphi);
            phase_cache_hit = (kc != NULL);
        }
        if (cache_off) {
            mie_value_kernel_allm = (double *)malloc(6u * all * sizeof(double));
            if (!mie_value_kernel_allm) { rc = -4; goto cleanup; }
            mie_value_kernel_owned = 1;
            if (rt_aerosol_value_phase_fourier_pol_allm(
                    &atm, mie_value_kernel_mcount, particle_value_phase, nphi,
                    mie_value_kernel_allm, mie_value_kernel_allm + all,
                    mie_value_kernel_allm + 2u*all, mie_value_kernel_allm + 3u*all,
                    mie_value_kernel_allm + 4u*all, mie_value_kernel_allm + 5u*all) != 0) {
                rc = -4; goto cleanup;
            }
        } else if (!kc) {
            kc = wpvkc_acquire(particle_value_phase_key, lambda_nm,
                               o.water_mie_truncation_mode, o.water_mie_ss_mode,
                               wmc_moment_mode, wmc_nmg, wmc_lmax_req,
                               &atm, mie_value_kernel_mcount, nphi);
            if (!kc) { rc = -4; goto cleanup; }
            if (rt_aerosol_value_phase_fourier_pol_allm(
                    &atm, mie_value_kernel_mcount, particle_value_phase, nphi,
                    kc->base_buf, kc->base_buf + all,
                    kc->base_buf + 2u*all, kc->base_buf + 3u*all,
                    kc->base_buf + 4u*all, kc->base_buf + 5u*all) != 0) {
                rc = -4; goto cleanup;
            }
            kc->base_mu_bits = wpvkc_double_bits(atm.rm[0]);
            ++kc->full_build_count;
            if (wpvkc_store_beam(kc, atm.rm[0], kc->base_buf) < 0) {
                rc = -4; goto cleanup;
            }
            ++kc->beam_build_count;
            mie_value_kernel_allm = kc->base_buf;
        } else {
            ++kc->phase_hit_count;
            const unsigned long long mu_bits=wpvkc_double_bits(atm.rm[0]);
            if (mu_bits == kc->base_mu_bits) {
                mie_value_kernel_allm = kc->base_buf;
                ++kc->beam_hit_count;
                beam_cache_hit = 1;
                beam_state = "base";
            } else {
                double *work=wpvkc_work_buffer(6u*all);
                if (!work) { rc=-4; goto cleanup; }
                memcpy(work, kc->base_buf, 6u*all*sizeof(double));
                water_pol_value_beam_entry_t *be = wpvkc_find_beam(kc, atm.rm[0]);
                if (be) {
                    wpvkc_apply_strip(kc, work, be->strip);
                    ++kc->beam_hit_count;
                    beam_cache_hit = 1;
                    beam_state = "hit";
                } else {
                    if (rt_aerosol_value_phase_fourier_pol_solar_allm(
                            &atm, mie_value_kernel_mcount, particle_value_phase, nphi,
                            work, work + all, work + 2u*all, work + 3u*all,
                            work + 4u*all, work + 5u*all) != 0) {
                        rc = -4; goto cleanup;
                    }
                    if (wpvkc_store_beam(kc, atm.rm[0], work) < 0) {
                        rc = -4; goto cleanup;
                    }
                    ++kc->beam_build_count;
                }
                mie_value_kernel_allm = work;
            }
        }
        /* Build the incoming-Q solar column unconditionally.  atm.beam_q is
         * assigned later, after the air-water Fresnel transmission matrix is
         * evaluated; the column itself depends only on the phase and solar
         * direction rm[0], which are already available here. */
        {
            const size_t beam_count = (size_t)mie_value_kernel_mcount *
                                      (size_t)(2 * n_mu + 1);
            mie_value_beamq_allm = (double *)malloc(3u * beam_count * sizeof(double));
            if (!mie_value_beamq_allm) { rc = -4; goto cleanup; }
            if (rt_aerosol_value_phase_fourier_pol_beamq_allm(
                    &atm, mie_value_kernel_mcount, particle_value_phase, nphi,
                    mie_value_beamq_allm,
                    mie_value_beamq_allm + beam_count,
                    mie_value_beamq_allm + 2u * beam_count) != 0) {
                rc = -4; goto cleanup;
            }
        }
        const char *vtrace = getenv("OCRT_VALUE_KERNEL_CACHE_TRACE");
        if (vtrace && vtrace[0] && strcmp(vtrace, "0") != 0) {
            const unsigned long long full_builds=kc?kc->full_build_count:0ULL;
            const unsigned long long phase_hits=kc?kc->phase_hit_count:0ULL;
            const unsigned long long beam_builds=kc?kc->beam_build_count:0ULL;
            const unsigned long long beam_hits=kc?kc->beam_hit_count:0ULL;
            fprintf(stderr,
                    "[VALUE-POL-CACHE] phase=%s beam=%s nmu=%d m=%d nphi=%d "
                    "entry_full_builds=%llu entry_phase_hits=%llu "
                    "entry_beam_builds=%llu entry_beam_hits=%llu "
                    "entries=%d resident=%zu\n",
                    phase_cache_hit ? "hit" : "miss",
                    cache_off ? "disabled" : (beam_cache_hit ? beam_state : "build"),
                    n_mu, mie_value_kernel_mcount, nphi,
                    full_builds, phase_hits, beam_builds, beam_hits,
                    wpvkc_entry_count(), wpvkc_resident_bytes_all());
        }
    }

    /* Path B value-kernel cache: compute the phase signature once and decide
     * HIT (reuse cached per-m kernel) vs MISS (build + store).  Applies only to
     * the value-kernel paths (external LUT or internal FF with value kernel);
     * other phase models build their kernels cheaply per m. */
    int use_value_kernel = (ccrr_particle_lut_mode ||
                            (ccrr_particle_ff_mode && ccrr_ff_value_kernel));
    int wpkc_hit = 0;
    int wpkc_disable = g_wrt_env.f_OCRT_NO_PHASE_CACHE;
    if (use_value_kernel) {
        int wpkc_nphi = (o.ccrr_particle_phase_nphi > 0) ? o.ccrr_particle_phase_nphi : 720;
        double sig_nphi = (double)wpkc_nphi;
        int m_count = m_loop_max + 1;
        if (!wpkc_disable && wpkc_match(lambda_nm, o.ccrr_chl_mg_m3, o.ccrr_min_g_m3,
                       o.a_cdom_440_m_inv, o.S_cdom_nm_inv, n_water, T_water_C,
                       sig_nphi, mu_view_water, n_mu, m_count, o.water_phase_kernel,
                       ccrr_particle_lut_mode, o.ccrr_particle_phase_lut_path)) {
            wpkc_hit = 1;
        } else {
            wpkc_reset(lambda_nm, o.ccrr_chl_mg_m3, o.ccrr_min_g_m3,
                       o.a_cdom_440_m_inv, o.S_cdom_nm_inv, n_water, T_water_C,
                       sig_nphi, mu_view_water, n_mu, m_count, o.water_phase_kernel,
                       ccrr_particle_lut_mode, o.ccrr_particle_phase_lut_path);
        }
    }
    /* Water-side surface internal-reflection feedback: feed the upwelling
     * reflected by Fresnel/TIR at the air-water interface back into the SOS as a
     * top downward BC, so it re-enters and scatters (self-consistent reflecting
     * boundary). Precompute per-node 3x3 internal-reflection Mueller M_Rww at each
     * in-water node; TIR (theta_w>theta_c) -> Rww=1 via rt_air_water_R_ww. Flat
     * surface only (wind=0); wind>0 Cox-Munk feedback not yet implemented.
     *
     * This feedback is physically required and is ALWAYS ON in production — it is
     * NOT a user option. v1.06 fixed a bug where it was off, leaving the internal
     * reflection in the Ed(0-) bookkeeping (Ed_internal_reflect, below) but NOT in
     * the upwelling source — Lu(0-) was under-driven while Ed(0-) was inflated,
     * deflating rrs(0-) (isotropic omega0=0.95: -21% vs the exact Chandrasekhar
     * H-function and HydroLight; +0.1% with feedback on).
     *
     * DEBUG-ONLY (NOT in user CLI/--help/user docs): the disable switch
     * OCRT_DEBUG_AW_FEEDBACK_OFF reproduces the pre-fix (incorrect) behavior for
     * regression/isolation. It is gated behind the master debug flag OCRT_DEBUG,
     * so it is inert in normal use and cannot silently corrupt a production run;
     * disabling requires OCRT_DEBUG=1 AND OCRT_DEBUG_AW_FEEDBACK_OFF=1. */
    static int p0a_fb = -1;
    if (p0a_fb < 0) p0a_fb = 1; /* air-water feedback: always on (physical) */
    double *Rww_M = NULL;
    if (p0a_fb && o.wind_speed <= 0.0) {
        Rww_M = calloc((size_t)n_mu * 9, sizeof(double));
        if (Rww_M) {
            double num = 0.0, den = 0.0;
            for (int jp = 1; jp <= n_mu; ++jp) {
                double M_Rww[9];
                rt_air_water_R_ww(atm.rm[+jp], n_water, o.q_convention, M_Rww);
                memcpy(Rww_M + (size_t)(jp - 1) * 9, M_Rww, 9 * sizeof(double)); /* full 3x3 Mueller */
                num += M_Rww[0] * atm.rm[+jp] * atm.gb[+jp];   /* I->I for r_bar diag */
                den += atm.rm[+jp] * atm.gb[+jp];
            }
            { int di = g_wrt_env.f_OCRT_DUMP_IOP;
              if (di) fprintf(stderr, "P0A_FEEDBACK on: flux-weighted internal reflectance r_bar=%.4f\n", (den>0)?num/den:0.0); }
        }
    }
    /* B2 (2026-06-03): rough-surface (wind>0) water-side internal-reflection
     * feedback. The Cox-Munk angle-coupling kernel depends on the Fourier mode m,
     * so the m-mode kernel Rww_K is (re)built inside the m-loop below; here we
     * only allocate the reusable buffer + the positive-mu node array. wind<=0
     * keeps the per-node flat path (Rww_M) above. Same p0a_fb debug gate. */
    double *Rww_K = NULL, *mu_pos_ww = NULL;
    int n_phi_ww = 128;   /* rough Rww Fourier integration; debug-overridable for parity audits */
    {
        const char *np = ocrt_debug_env("OCRT_RWW_NPHI");
        if (np && np[0]) {
            int v = atoi(np);
            if (v >= 32 && v <= 4096) n_phi_ww = v;
        } else {
            int min_phi = 2 * o.m_max_water + 16;
            if (min_phi > n_phi_ww) n_phi_ww = min_phi;
            if (n_phi_ww < 176) n_phi_ww = 176;
        }
    }
    if (p0a_fb && o.wind_speed > 0.0) {
        Rww_K     = calloc((size_t)n_mu * (size_t)n_mu * 9, sizeof(double));
        mu_pos_ww = calloc((size_t)n_mu, sizeof(double));
        if (Rww_K && mu_pos_ww) {
            for (int j = 1; j <= n_mu; ++j) mu_pos_ww[j - 1] = atm.rm[+j];
        } else { free(Rww_K); free(mu_pos_ww); Rww_K = NULL; mu_pos_ww = NULL; }
    }
    /* v1.059: incident-beam polarization. The solar beam refracted into the
     * water is partially polarized by air-water Fresnel transmission
     * (q_beam = M_T_aw[1,0]/M_T_aw[0,0] at the air-side solar zenith). The
     * pre-v1.059 primary source treated the in-water beam as unpolarized; set
     * atm.beam_q so rt_solver_primary_source(_pol) add the 2nd-column phase
     * coupling of the beam Q. Vector path only (this function) — scalar
     * rt_solve_case leaves beam_q=0, so the §6 scalar regression is invariant.
     * Same q_convention as the surface reflection/transmission Mueller. */
    {
        double M_T_aw[9];
        rt_air_water_T_aw(mu_sun_air, n_water, o.q_convention, M_T_aw);
        atm.beam_q = (M_T_aw[0] != 0.0) ? (M_T_aw[3] / M_T_aw[0]) : 0.0;
        /* [FIX-SKY-EDLU] sub-cone equivalent-beam has no air incidence angle, so
         * the air->water transmission Q-induction is undefined; treat the sub-cone
         * skylight beam as unpolarized (beam_q=0).  Its Q/U contribution to the
         * nadir rrs is second-order (small sub-cone fraction); the intensity Lu
         * (which drives rrs) is carried exactly by mu_sun_water_override above. */
        if (o.mu_sun_water_override > 0.0) atm.beam_q = 0.0;
        int dbq = g_wrt_env.f_OCRT_DUMP_BEAMQ;
        if (dbq) fprintf(stderr, "BEAMQ mu_sun_air=%.4f n_water=%.3f q_beam=%.6e\n",
                         mu_sun_air, n_water, atm.beam_q);
    }
    /* v1.09 commit #16: grid-cache key check (env-gated). */
    const char *gcv_ = getenv("OCRT_WATER_GRID_CACHE");
    /* #22 fix of a #21 regression: this variable is SET at runtime (setenv in
     * the full-grid/batch drivers) AFTER the base solve's lazy snapshot init,
     * so a snapshot read froze it to NULL and silently disabled the #16 grid
     * cache (~600x).  Runtime-controlled variables are snapshot-ineligible;
     * this is a once-per-sos_pure read, so a live getenv costs nothing. */
    int grid_cache_on = (gcv_ && strcmp(gcv_, "0") != 0) && (o.view_as_node == 0)
                        && !o.bypass_grid_cache;   /* D3 driver solves: no lookup, no store */
    if (getenv("OCRT_S6_TRACE"))
        fprintf(stderr,"[S6K] sza=%.9g F=%.6g nmw=%d mmw=%d nlw=%d tol=%g nt=%d nmu=%d dirs=%d tmax=%g a=%.9g b=%.9g bb=%.9g wpk=%d nvv=%d vasn=%d Twc=%g Sgk=%g nw=%.9g wind=%g maxit=%d\n",
                sza_deg_air, F_sun, o.n_mu_water, o.m_max_water, o.n_layers_water,
                o.tolerance, nt, n_mu, dirs, o.tau_max_target,
                o.fixed_a_total_m_inv, o.fixed_b_total_m_inv, o.fixed_bb_total_m_inv,
                o.water_phase_kernel, o.n_view_vza, o.view_as_node,
                T_water_C, S_water_gkg, n_water, o.wind_speed, o.max_iterations);
    { static int beam_loaded = 0;   /* once per process (thread 0 usage pattern) */
      const char *bc = getenv("OCRT_D3_BEAM_CACHE");
      if (bc && bc[0] && !o.bypass_grid_cache && grid_cache_on && !beam_loaded) {
          beam_loaded = 1;
          if (ocrt_beam_cache_load(bc) == 0)
              fprintf(stderr, "[D3beam] cache loaded: %s\n", bc);
      } }
    int grid_hit = 0;
#ifdef OCRT_FAST_KERNELS
    for (int ocrt_s = 0; ocrt_s < OCRT_WGC_SLOTS && !grid_hit; ++ocrt_s) {
        ocrt_wgc_cur_i = ocrt_s;
#endif
    if (grid_cache_on && !o.bypass_grid_cache && g_grid_cache.valid &&
        g_grid_cache.sza==sza_deg_air && g_grid_cache.lam==lambda_nm &&
        g_grid_cache.Twc==T_water_C && g_grid_cache.Sgk==S_water_gkg &&
        g_grid_cache.nw==n_water && g_grid_cache.Fs==F_sun &&
        g_grid_cache.wind==o.wind_speed && g_grid_cache.tol==o.tolerance &&
        g_grid_cache.nmw==o.n_mu_water && g_grid_cache.mmw==o.m_max_water &&
        g_grid_cache.nlw==o.n_layers_water && g_grid_cache.maxit==o.max_iterations &&
        g_grid_cache.vasn==o.view_as_node && g_grid_cache.wpk==o.water_phase_kernel &&
        g_grid_cache.nvv==o.n_view_vza &&
        g_grid_cache.a_t==o.fixed_a_total_m_inv && g_grid_cache.b_t==o.fixed_b_total_m_inv &&
        g_grid_cache.bb_t==o.fixed_bb_total_m_inv && g_grid_cache.taumax==o.tau_max_target &&
        g_grid_cache.a_iop==a_tot && g_grid_cache.b_iop==b_tot &&
        g_grid_cache.bb_iop==bb_tot &&
        g_grid_cache.b_phyto==b_phyto && g_grid_cache.b_detritus==b_det &&
        g_grid_cache.b_mineral==b_min &&
        g_grid_cache.constituent_mode==o.ccrr_mode &&
        g_grid_cache.constituent_model==o.water_constituent_model &&
        g_grid_cache.organic_group==o.organic_phyto_group &&
        g_grid_cache.tsm_species==o.tsm_species &&
        g_grid_cache.moment_mode==o.water_mie_moment_mode &&
        g_grid_cache.moment_nmg==o.water_mie_moment_n_mu &&
        g_grid_cache.phase_lmax==o.fixed_bulk_lmax &&
        g_grid_cache.mie_trunc==o.water_mie_truncation_mode &&
        g_grid_cache.mie_ss==o.water_mie_ss_mode &&
        g_grid_cache.value_phase_spline==g_wrt_env.f_OCRT_VALUE_PHASE_SPLINE &&
        g_grid_cache.value_kernel_pol==g_wrt_env.f_OCRT_WATER_VALUE_KERNEL_POL &&
        g_grid_cache.nt==nt && g_grid_cache.n_mu==n_mu && g_grid_cache.dirs==dirs &&
        !strncmp(g_grid_cache.wmp, o.water_mie_phase_path?o.water_mie_phase_path:"", 511) &&
        !strncmp(g_grid_cache.fbl, o.fixed_bulk_phase_lut_path?o.fixed_bulk_phase_lut_path:"", 511))
        grid_hit = 1;
#ifdef OCRT_FAST_KERNELS
        if (grid_hit) break;
    }
    if (grid_hit) {
        g_wgc_stamp[ocrt_wgc_cur_i] = ++g_wgc_tick;
    } else {
        /* victim = LRU (invalid slots first) */
        int ocrt_v = 0;
        for (int ocrt_s = 0; ocrt_s < OCRT_WGC_SLOTS; ++ocrt_s) {
            if (!g_wgc_slots[ocrt_s].valid) { ocrt_v = ocrt_s; break; }
            if (g_wgc_stamp[ocrt_s] < g_wgc_stamp[ocrt_v]) ocrt_v = ocrt_s;
        }
        ocrt_wgc_cur_i = ocrt_v;
        g_wgc_stamp[ocrt_v] = ++g_wgc_tick;
    }
#endif
    if (grid_cache_on && !grid_hit) {
        grid_cache_free();
        size_t fs=(size_t)(nt+1)*(size_t)dirs;
        g_grid_cache.tot  = (double*)malloc((size_t)(o.m_max_water+1)*3*fs*sizeof(double));
        g_grid_cache.prim0= (double*)malloc((size_t)(o.m_max_water+1)*3*(size_t)dirs*sizeof(double));
        if (!g_grid_cache.tot || !g_grid_cache.prim0) { grid_cache_free(); grid_cache_on = 0; }
    }
    double m_amp_max = 0.0;      /* commit #14: running max mode amplitude */
    int    m_small_streak = 0;   /* commit #14: consecutive sub-tolerance modes */
    int    m_early_exit_at = -1; /* commit #14: -1 = ran full m_loop_max */
    (void)m_early_exit_at;
#ifdef OCRT_FAST_KERNELS
    /* The direct Stage-3B arrays already contain all six particle kernels.
     * Do not initialize the finite-moment particle-kernel cache when they are
     * active; only the molecular Legendre basis is still required per mode. */
    const int ocrt_mkc_ctx = mie_value_kernel_allm ? -1 :
        rt_mkc_begin(&atm, &ws, m_loop_max + 1);
#endif
    if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] wmloop begin abs=%.4f\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
    surface_fkc_set_build_mmax(m_loop_max);   /* v1.11-speed S1: cap R_ww all-m build */
    for (int m = 0; m <= m_loop_max; ++m) {
        if (grid_hit) {
            if (m >= g_grid_cache.m_count) break;   /* cold run early-exited here */
            size_t fs=(size_t)(nt+1)*(size_t)dirs;
#ifdef OCRT_FAST_KERNELS
            /* v1.10 S12 (2026-07-10): REPLAY ALIASING, no copy.
             * Audit (this commit): every tot_* use on the replay path is a
             * READ (view/node extraction, hemispheric sums); writes happen
             * only on the miss path (SOS output 3773+, cache store 3802+),
             * unreachable when grid_hit.  So point tot_* straight at the
             * cache planes for this m instead of streaming 3*fs doubles
             * (~0.9 MB) per mode, ~28 MB per warm row.  The arena slices
             * from WSLICE stay allocated (reused by the next miss row).
             * Strict branch keeps the memcpy verbatim. */
            tot_i = g_grid_cache.tot + ((size_t)m*3+0)*fs;
            tot_q = g_grid_cache.tot + ((size_t)m*3+1)*fs;
            tot_u = g_grid_cache.tot + ((size_t)m*3+2)*fs;
#else
            memcpy(tot_i, g_grid_cache.tot + ((size_t)m*3+0)*fs, fs*sizeof(double));
            memcpy(tot_q, g_grid_cache.tot + ((size_t)m*3+1)*fs, fs*sizeof(double));
            memcpy(tot_u, g_grid_cache.tot + ((size_t)m*3+2)*fs, fs*sizeof(double));
#endif
            max_orders_seen = g_grid_cache.max_orders_seen_c;
            all_conv        = g_grid_cache.all_conv_c;
            worst_resid     = g_grid_cache.worst_resid_c;
            goto grid_replay_extract;
        }
        /* 7a. Compute Legendre + phase Fourier kernels for this m */
        if (rt_legendre_compute_pol(&ws, &atm, m) != 0) {
            fprintf(stderr, "[B3] rt_legendre_compute_pol FAIL m=%d\n", m);
            rc = -3; goto cleanup;
        }
        const int direct_particle_for_m =
            (mie_value_kernel_allm && m < mie_value_kernel_mcount);
        ws.beam_q_direct_valid = 0;
        if (!direct_particle_for_m) {
        /* CCRR scalar particle phase: the intensity kernel used by both
         * rt_solver_primary_source() and rt_sos_operator_apply_vector() lives in
         * ws->phase_fourier_m and is built from betal_aer.  Earlier CCRR
         * ports only populated the polarization kernels (gr/gt/arr/art/att)
         * and left phase_fourier_m at zero, suppressing pigment/mineral
         * scattering by orders of magnitude. */
#ifdef OCRT_FAST_KERNELS
        const int ocrt_mkc_served = (ocrt_mkc_ctx == 1) && rt_mkc_load(m, &ws);
        if (!ocrt_mkc_served) {
#endif
        if (rt_kernel_phase_fourier(&ws, m, atm.betal_aer) != 0) {
            fprintf(stderr, "[B3] rt_kernel_phase_fourier FAIL m=%d\n", m);
            rc = -3; goto cleanup;
        }
#ifdef OCRT_FAST_KERNELS
        if (!ocrt_mkc_served && ocrt_mkc_ctx >= 0)
            rt_mkc_store_pfm(m, &ws);
        }  /* !ocrt_mkc_served (step2) */
#endif
        if (o.fixed_bulk_iop_mode && fixed_bulk_ff_betal_L < 0 && !(o.water_mie_phase_path && o.water_mie_phase_path[0]) && (fixed_bulk_phase_model == 1 || fixed_bulk_phase_model == 2 || fixed_bulk_phase_model == 3 || fixed_bulk_phase_model == 4)) {
            double bbob = (b_tot > 0.0) ? (bb_tot / b_tot) : 0.0;
            double g_hg = rt_iop_hg_g_for_backscatter_fraction(bbob);
            int nphi = (o.fixed_bulk_phase_nphi > 0) ? o.fixed_bulk_phase_nphi : 720;
            const fixed_bulk_phase_table_t *tabp = (fixed_bulk_phase_model == 2 || fixed_bulk_phase_model == 3 || fixed_bulk_phase_model == 4) ? &fixed_phase_table : NULL;
            /* v1.09 commit #8: build ALL m planes once per beam (bit-identical
             * single-sweep, see fixed_bulk_direct_phase_fourier_allm), then
             * serve per-m copies.  Diagnostics envs keep the original path. */
            if (g_wrt_env.s_OCRT_FIXEDBULK_PHASENORM || g_wrt_env.s_OCRT_DUMP_NORMBB) {
                if (fixed_bulk_direct_phase_fourier(&ws, &atm, m, g_hg, tabp, have_fixed_phase_spline ? &fixed_phase_spline : NULL, nphi) != 0) {
                    fprintf(stderr, "[B3] fixed_bulk_direct_phase_fourier FAIL m=%d\n", m);
                    rc = -3; goto cleanup;
                }
            } else {
                const int kw_fb = 2*n_mu+1;
                if (!fb_allm_kernel) {
                    fb_allm_mcount = m_loop_max + 1;
                    fb_allm_kernel = (double*)malloc((size_t)fb_allm_mcount*(n_mu+1)*kw_fb*sizeof(double));
                    if (!fb_allm_kernel ||
                        fixed_bulk_direct_phase_fourier_allm(n_mu, &atm, fb_allm_mcount,
                                                             g_hg, tabp, have_fixed_phase_spline ? &fixed_phase_spline : NULL, nphi, fb_allm_kernel) != 0) {
                        fprintf(stderr, "[B3] fixed_bulk_direct_phase_fourier_allm FAIL\n");
                        rc = -3; goto cleanup;
                    }
                }
                for (int jj=0; jj<=n_mu; jj++)
                    memcpy(&ws.phase_fourier_m[jj][-n_mu],
                           &fb_allm_kernel[((size_t)m*(n_mu+1)+(size_t)jj)*kw_fb],
                           (size_t)kw_fb*sizeof(double));
            }
        } else if ((ccrr_particle_lut_mode || (ccrr_particle_ff_mode && ccrr_ff_value_kernel)) && b_particle > 0.0) {
            /* Value-based direct phase-value Fourier kernel (path B).  Cached
             * across beams at the same wavelength (see wpkc above): build on the
             * first beam (MISS), reuse on the rest (HIT). */
            if (wpkc_hit) {
                int nphi = (o.ccrr_particle_phase_nphi > 0) ? o.ccrr_particle_phase_nphi : 720;
                wpkc_load_m(m, &ws, n_mu);                                    /* beam-independent GL block from cache */
                fixed_bulk_direct_phase_fourier_solar(&ws, &atm, m,
                                                      &ccrr_particle_phase_table, have_ccrr_phase_spline ? &ccrr_phase_spline : NULL, nphi); /* solar coupling for THIS beam */
            } else {
                int nphi = (o.ccrr_particle_phase_nphi > 0) ? o.ccrr_particle_phase_nphi : 720;
                if (fixed_bulk_direct_phase_fourier(&ws, &atm, m, 0.0, &ccrr_particle_phase_table, have_ccrr_phase_spline ? &ccrr_phase_spline : NULL, nphi) != 0) {
                    fprintf(stderr, "[B3] ccrr_particle_phase_lut_fourier FAIL m=%d\n", m);
                    rc = -3; goto cleanup;
                }
                wpkc_store_m(m, &ws, n_mu);
            }
        }
#ifdef OCRT_FAST_KERNELS
        if (!ocrt_mkc_served) {
#endif
        if (rt_kernel_phase_fourier_pol(&ws, m, atm.gammal_aer) != 0) {
            fprintf(stderr, "[B3] rt_kernel_phase_fourier_pol FAIL m=%d\n", m);
            rc = -3; goto cleanup;
        }
        /* Also compute full vector aerosol kernel (arr/art/att) for SOS_pol */
        if (rt_kernel_phase_fourier_aerosol_full(&ws, m,
                                                   atm.alphal_aer,
                                                   atm.zetal_aer) != 0) {
            fprintf(stderr, "[B3] rt_kernel_phase_fourier_aerosol_full FAIL m=%d\n", m);
            rc = -3; goto cleanup;
        }
#ifdef OCRT_FAST_KERNELS
        if (ocrt_mkc_ctx >= 0) rt_mkc_store_pol(m, &ws);
        }  /* !ocrt_mkc_served (steps 4-5) */
#endif
        } /* legacy finite-moment / older value-kernel particle paths */
        if (direct_particle_for_m) {
            const size_t plane=(size_t)(n_mu+1)*(size_t)dirs;
            const size_t all=(size_t)mie_value_kernel_mcount*plane;
            const size_t off=(size_t)m*plane;
            memcpy(ws.pfm_storage,mie_value_kernel_allm+off,plane*sizeof(double));
            memcpy(ws.gr_storage,mie_value_kernel_allm+all+off,plane*sizeof(double));
            memcpy(ws.gt_storage,mie_value_kernel_allm+2u*all+off,plane*sizeof(double));
            memcpy(ws.arr_storage,mie_value_kernel_allm+3u*all+off,plane*sizeof(double));
            memcpy(ws.art_storage,mie_value_kernel_allm+4u*all+off,plane*sizeof(double));
            memcpy(ws.att_storage,mie_value_kernel_allm+5u*all+off,plane*sizeof(double));
            if (mie_value_beamq_allm) {
                const size_t beam_plane = (size_t)dirs;
                const size_t beam_all = (size_t)mie_value_kernel_mcount * beam_plane;
                const size_t beam_off = (size_t)m * beam_plane;
                memcpy(ws.beam_q_to_i_direct,
                       mie_value_beamq_allm + beam_off,
                       beam_plane * sizeof(double));
                memcpy(ws.beam_q_to_q_direct,
                       mie_value_beamq_allm + beam_all + beam_off,
                       beam_plane * sizeof(double));
                memcpy(ws.beam_q_to_u_direct,
                       mie_value_beamq_allm + 2u * beam_all + beam_off,
                       beam_plane * sizeof(double));
                ws.beam_q_direct_valid = 1;
            }
            ws.sos_kernel_pack_valid=0; ws.sos_kernel_pack_m=-1;
        }
        if (rt_sos_operator_prepare(&atm, m, &ws) != 0) {
            fprintf(stderr, "[B3] rt_sos_operator_prepare FAIL m=%d\n", m);
            rc = -3; goto cleanup;
        }

        /* 7b. Primary source
         * Diagnostic NT-TMS parity mode: OSOAA HYD_TRUNCATION=ON source traces
         * show order-1 uses the raw hydrosol phase trend, while IG>=2 source
         * kernels use truncated PM + b_eff/tau_tr.  Build primary source with
         * a raw-moment workspace but keep atm.xdel/ch/h at the effective
         * truncated optical budget.  The subsequent SOS call below still uses
         * the default workspace ws, which is truncated in this route. */
        const rt_atm_t *atm_primary = &atm;
        const rt_legendre_workspace_t *ws_for_primary = &ws;
        rt_atm_t atm_primary_raw;
        if (water_nt_tms) {
            if (rt_legendre_compute_pol(&ws_primary, &atm, m) != 0 ||
                rt_kernel_phase_fourier(&ws_primary, m, mie_raw_betal_buf) != 0 ||
                rt_kernel_phase_fourier_pol(&ws_primary, m, mie_raw_gammal_buf) != 0 ||
                rt_kernel_phase_fourier_aerosol_full(&ws_primary, m, mie_raw_alphal_buf, mie_raw_zetal_buf) != 0) {
                fprintf(stderr, "[B3] NT-TMS raw primary kernel build FAIL m=%d\n", m);
                rc = -3; goto cleanup;
            }
            atm_primary_raw = atm;
            atm_primary_raw.betal_aer  = mie_raw_betal_buf;
            atm_primary_raw.gammal_aer = mie_raw_gammal_buf;
            atm_primary_raw.alphal_aer = mie_raw_alphal_buf;
            atm_primary_raw.zetal_aer  = mie_raw_zetal_buf;
            atm_primary = &atm_primary_raw;
            ws_for_primary = &ws_primary;
        }
        if (rt_solver_primary_source(atm_primary, m, ws_for_primary, src_i) != 0) {
            fprintf(stderr, "[B3] rt_solver_primary_source FAIL m=%d\n", m);
            rc = -3; goto cleanup;
        }
        if (rt_solver_primary_source_pol(atm_primary, m, ws_for_primary, src_q, src_u) != 0) {
            fprintf(stderr, "[B3] rt_solver_primary_source_pol FAIL m=%d\n", m);
            rc = -3; goto cleanup;
        }
        if (o.ext_top_I && o.ext_top_mu && o.ext_top_n > 0 &&
            m <= o.ext_top_m_max) {
            if (!getenv("OCRT_DTP_OFF")) ocrt_add_diffuse_top_primary(atm_primary, ws_for_primary, m,
                &o.ext_top_I[(size_t)m * (size_t)o.ext_top_n],
                o.ext_top_Q ? &o.ext_top_Q[(size_t)m * (size_t)o.ext_top_n] : NULL,
                o.ext_top_U ? &o.ext_top_U[(size_t)m * (size_t)o.ext_top_n] : NULL,  /* D3-0b */
                o.ext_top_mu, o.ext_top_n, src_i, src_q, src_u);
        }
        dump_water_primary_source_trace_if_requested(m, atm_primary, src_i, src_q, src_u);
        dump_water_primary_pref_trace_if_requested(m, atm_primary, src_i, src_q, src_u);

        /* NOTE: atmospheric SOS engine assumes F_sun=π normalization.
         * Output I is *F_sun=π normalized radiance* (i.e. radiance L when F_sun=π).
         * Do NOT multiply src by F_sun here — that creates double normalization.
         * Caller is responsible for scaling output by actual F_TOA at end.
         * For this entry, F_sun input is informational only (matches the
         * normalization of the SS analytic for cross-check). */
        (void)F_sun;  /* unused in source build; applied externally */

        /* 7c. Vertical integration of primary */
grid_replay_extract: ;
#ifdef OCRT_FAST_KERNELS
        /* v1.10 S12: on replay the only prim_* consumer is the Iss view
         * capture (layer-0 row, filled by the dirs-sized memcpy below);
         * k>0 layers are never read, so the full-field zeroing is dead
         * work on hit rows (~28 MB/row).  Miss rows keep it (the primary
         * integration may accumulate; arena slices carry prior content). */
        if (!grid_hit) {
            memset(prim_i, 0, (size_t)field_sz * sizeof(double));   /* #21: arena slice */
            memset(prim_q, 0, (size_t)field_sz * sizeof(double));
            memset(prim_u, 0, (size_t)field_sz * sizeof(double));
        }
#else
        memset(prim_i, 0, (size_t)field_sz * sizeof(double));   /* #21: arena slice */
        memset(prim_q, 0, (size_t)field_sz * sizeof(double));
        memset(prim_u, 0, (size_t)field_sz * sizeof(double));
#endif
        if (!grid_hit &&
            (rt_solver_integrate(&atm, m, src_i, RT_INTEGRATION_METHOD_LINEAR, prim_i) != 0 ||
             rt_solver_integrate(&atm, m, src_q, RT_INTEGRATION_METHOD_LINEAR, prim_q) != 0 ||
             rt_solver_integrate(&atm, m, src_u, RT_INTEGRATION_METHOD_LINEAR, prim_u) != 0)) {
            fprintf(stderr, "[B3] rt_solver_integrate FAIL m=%d\n", m);
            rc = -3; goto cleanup;
        }

        /* 7d. SOS iteration */
        rt_solver_sos_options_t sopts = {
            .max_iterations = o.max_iterations,
            .tolerance      = o.tolerance,
            .acceleration   = RT_SOS_ACCELERATION_PLAIN,
            .save_orders    = 0
        };
        rt_solver_sos_result_t sres = {0};
        int src_rc;
        if (grid_hit) {
            /* replay: fields already restored at loop entry; provide the k=0
             * primary row for the (unmodified) Iss capture below. */
            memcpy(&prim_i[0], g_grid_cache.prim0 + ((size_t)m*3+0)*(size_t)dirs, (size_t)dirs*sizeof(double));
            memcpy(&prim_q[0], g_grid_cache.prim0 + ((size_t)m*3+1)*(size_t)dirs, (size_t)dirs*sizeof(double));
            memcpy(&prim_u[0], g_grid_cache.prim0 + ((size_t)m*3+2)*(size_t)dirs, (size_t)dirs*sizeof(double));
            sres.converged = 1; sres.n_orders_used = 0; sres.final_residual = 0.0;
            src_rc = 0;
        } else {
            if (Rww_K) {
                /* B2 rough: build the Cox-Munk water-side internal-reflection m-mode
                 * kernel for this m, then feed it back as the SOS top downward BC. */
                static int b2flat = -1;
                if (b2flat < 0) b2flat = 0; /* default kernel (physical) */
                if (b2flat) {
                    /* DEBUG (OCRT_DEBUG gated): diagonal kernel carrying the flat
                     * per-node Mueller, R_m[j,j]=M_Rww[j]/(2pi mu_j w_j). Fed through
                     * the rough contraction it MUST reproduce the flat per-node limit
                     * (Lu0plus=1.015187e-02). Isolates solver/contraction from the
                     * Cox-Munk kernel. m-independent (M_Rww has no m dependence). */
                    memset(Rww_K, 0, (size_t)n_mu * (size_t)n_mu * 9 * sizeof(double));
                    for (int j = 1; j <= n_mu; ++j) {
                        double Mf[9];
                        rt_air_water_R_ww(atm.rm[+j], n_water, o.q_convention, Mf);
                        const double az_factor = (m == 0) ? (2.0 * RT_F_SOLAR_PI) : RT_F_SOLAR_PI;
                        double inv = 1.0 / (az_factor * atm.rm[+j] * atm.gb[+j]);
                        double *Kd = Rww_K + ((size_t)(j - 1) * (size_t)n_mu + (size_t)(j - 1)) * 9;
                        for (int e = 0; e < 9; ++e) Kd[e] = Mf[e] * inv;
                    }
                } else {
                    const int rww_rc = surface_R_ww_coxmunk_fourier_kernel(
                                                        mu_pos_ww, n_mu, mu_pos_ww, n_mu,
                                                        m, n_phi_ww, o.wind_speed,
                                                        o.cox_munk_sigma_type, n_water,
                                                        o.q_convention, Rww_K);
                    if (rww_rc != 0) {
                        src_rc = rww_rc;
                    } else {
                        src_rc = rt_solver_sos_pol_intrefl_rough(
                                                         &atm, m, &ws, prim_i, prim_q, prim_u,
                                                         Rww_K, &sopts, tot_i, tot_q, tot_u, &sres);
                    }
                }
                if (b2flat) {
                    src_rc = rt_solver_sos_pol_intrefl_rough(
                                                         &atm, m, &ws, prim_i, prim_q, prim_u,
                                                         Rww_K, &sopts, tot_i, tot_q, tot_u, &sres);
                }
            } else if (Rww_M) {
                src_rc = rt_solver_sos_pol_intrefl(&atm, m, &ws, prim_i, prim_q, prim_u,
                                                   Rww_M, &sopts, tot_i, tot_q, tot_u, &sres);
            } else {
                src_rc = rt_solver_sos_pol(&atm, m, &ws, prim_i, prim_q, prim_u,
                                            &sopts, tot_i, tot_q, tot_u, &sres);
            }
        }
        /* v1.054 TMS: capture single-scatter (primary-order) intensity at the view
         * direction BEFORE freeing prim_i (same linear node interpolation as the
         * total). Reconstructed across m below to form Lu_ss,deltaM for the
         * single-scatter exact-phase (Nakajima-Tanaka) correction. */
        {
            const double mu_t = mu_view_water;
            int jl = 1, jh = 1;
            if (mu_t <= atm.rm[+1]) { jl = 1; jh = 2; }
            else if (mu_t >= atm.rm[+n_mu]) { jl = n_mu - 1; jh = n_mu; }
            else { for (int jp = 1; jp < n_mu; ++jp) if (atm.rm[+jp] <= mu_t && mu_t < atm.rm[+jp+1]) { jl = jp; jh = jp+1; break; } }
            double ml = atm.rm[+jl], mh = atm.rm[+jh], dn = mh - ml;
            double whi = (dn != 0.0) ? (mu_t - ml) / dn : 0.0, wlo = 1.0 - whi;
            size_t il = (size_t)0 * (size_t)dirs + (size_t)(+jl + n_mu);
            size_t ih = (size_t)0 * (size_t)dirs + (size_t)(+jh + n_mu);
            Iss_m_view[m] = wlo * prim_i[il] + whi * prim_i[ih];
            Qss_m_view[m] = wlo * prim_q[il] + whi * prim_q[ih];
            Uss_m_view[m] = wlo * prim_u[il] + whi * prim_u[ih];
        }
        if (grid_cache_on && !grid_hit && g_grid_cache.tot) {
            size_t fs=(size_t)(nt+1)*(size_t)dirs;
            memcpy(g_grid_cache.tot  + ((size_t)m*3+0)*fs, tot_i, fs*sizeof(double));
            memcpy(g_grid_cache.tot  + ((size_t)m*3+1)*fs, tot_q, fs*sizeof(double));
            memcpy(g_grid_cache.tot  + ((size_t)m*3+2)*fs, tot_u, fs*sizeof(double));
            memcpy(g_grid_cache.prim0+ ((size_t)m*3+0)*(size_t)dirs, prim_i, (size_t)dirs*sizeof(double));
            memcpy(g_grid_cache.prim0+ ((size_t)m*3+1)*(size_t)dirs, prim_q, (size_t)dirs*sizeof(double));
            memcpy(g_grid_cache.prim0+ ((size_t)m*3+2)*(size_t)dirs, prim_u, (size_t)dirs*sizeof(double));
            g_grid_cache.m_count = m + 1;
        }
        if (src_rc != 0) {
            fprintf(stderr, "[B3] rt_solver_sos_pol FAIL m=%d, rc=%d\n", m, src_rc);
            rc = -3; goto cleanup;
        }
        if (!grid_hit) {
            if (sres.n_orders_used > max_orders_seen) max_orders_seen = sres.n_orders_used;
            if (!sres.converged) all_conv = 0;
            if (sres.final_residual > worst_resid) worst_resid = sres.final_residual;
        }

        /* 7e. Extract upwelling at z=0- (k=0 level) at view-mu via
         *      LINEAR interpolation from positive Gauss nodes.
         *      사용자 작업 규칙 11: nearest 사용 금지. 양 끝단 외삽 시에도
         *      linear (with edge clamping to neighbor pair) 사용. */
        {
            /* atm.rm[+1..+n_mu] holds positive Gauss nodes; rt_atm_build_rayleigh
             * (and our build_inwater_atm) places them in ascending order:
             *   rm[+1] = smallest μ ≈ 0.04
             *   rm[+n_mu] = largest μ ≈ 0.99 */
            const double mu_target = mu_view_water;
            /* v1.09 #20: if the target direction exists as an (appended
             * zero-weight or coincident Gauss) node, read it directly —
             * no interpolation.  snell_down here and in the node-append
             * block use identical inputs, so the match is exact; the
             * tolerance only guards benign last-ulp effects of sorting. */
            int j_exact = -1;
            for (int jp = 1; jp <= n_mu; ++jp) {
                if (fabs(atm.rm[+jp] - mu_target) < 1e-13) { j_exact = jp; break; }
            }
            if (j_exact > 0) {
                size_t idx_e = (size_t)0 * (size_t)dirs + (size_t)(+j_exact + n_mu);
                I_m_view[m] = tot_i[idx_e];
                Q_m_view[m] = tot_q[idx_e];
                U_m_view[m] = tot_u[idx_e];
            } else if (n_mu >= 4 &&
                       mu_target > atm.rm[+1] && mu_target < atm.rm[+n_mu]) {
                /* 2026-07-14 W-PCHIP (Jae 승인): 내부 구간의 per-m 뷰 추출을
                 * 선형 -> PCHIP으로 교체(대기 LUT-PCHIP과 동일 정책).  범위
                 * 밖 외삽과 n_mu<4는 아래 기존 선형 경로 유지(bit-parity).
                 * 위 j_exact 직접 읽기(#20)가 노드 일치를 이미 처리하므로
                 * 이 분기는 진짜 off-node에서만 발효된다. */
                I_m_view[m] = ocrt_wpchip_eval_up_nodes(&atm, tot_i, n_mu, mu_target);
                Q_m_view[m] = ocrt_wpchip_eval_up_nodes(&atm, tot_q, n_mu, mu_target);
                U_m_view[m] = ocrt_wpchip_eval_up_nodes(&atm, tot_u, n_mu, mu_target);
            } else {
            int j_lo = 1, j_hi = 1;
            if (mu_target <= atm.rm[+1]) {
                /* Below smallest node — linear extrapolation using two lowest */
                j_lo = 1; j_hi = 2;
            } else if (mu_target >= atm.rm[+n_mu]) {
                /* Above largest node — linear extrapolation using two highest */
                j_lo = n_mu - 1; j_hi = n_mu;
            } else {
                /* Bracket: find j such that rm[+j] <= mu_target < rm[+j+1] */
                for (int jp = 1; jp < n_mu; ++jp) {
                    if (atm.rm[+jp] <= mu_target && mu_target < atm.rm[+jp+1]) {
                        j_lo = jp; j_hi = jp + 1; break;
                    }
                }
            }
            double mu_lo = atm.rm[+j_lo];
            double mu_hi = atm.rm[+j_hi];
            double denom = mu_hi - mu_lo;
            double w_hi  = (denom != 0.0) ? (mu_target - mu_lo) / denom : 0.0;
            double w_lo  = 1.0 - w_hi;
            size_t idx_lo = (size_t)0 * (size_t)dirs + (size_t)(+j_lo + n_mu);
            size_t idx_hi = (size_t)0 * (size_t)dirs + (size_t)(+j_hi + n_mu);
            I_m_view[m] = w_lo * tot_i[idx_lo] + w_hi * tot_i[idx_hi];
            Q_m_view[m] = w_lo * tot_q[idx_lo] + w_hi * tot_q[idx_hi];
            U_m_view[m] = w_lo * tot_u[idx_lo] + w_hi * tot_u[idx_hi];
            }
            /* v1.09 commit #14: Fourier-mode early exit.  The azimuthal mode
             * amplitudes of the upwelling field decay rapidly (multiple
             * scattering flattens azimuthal structure); measured on the
             * GOCI-III full-grid repro case (high-omega vector water phase,
             * m_max_water=30) the spectrum falls from |I_0|=4e-2 to 1.5e-8 at
             * m=8 and underflows (1e-63) by m=30 — 23 of 31 modes were pure
             * waste.  Exit when TWO consecutive modes (m>=3, so the Rayleigh
             * band m<=2 is always computed and legacy m_max_water=2 runs are
             * structurally untouched) contribute below opts.tolerance relative
             * to the running mode maximum.  Reuses the SOS tolerance — no new
             * knob.  Escape hatch: OCRT_WATER_M_NO_EARLY_EXIT=1. */
            {
                int no_exit = g_wrt_env.f_OCRT_WATER_M_NO_EARLY_EXIT;
                double amp = fabs(I_m_view[m]);
                if (fabs(Q_m_view[m]) > amp) amp = fabs(Q_m_view[m]);
                if (fabs(U_m_view[m]) > amp) amp = fabs(U_m_view[m]);
                if (amp > m_amp_max) m_amp_max = amp;
                if (!no_exit && m >= 3 && m_amp_max > 0.0) {
                    if (amp < o.tolerance * m_amp_max) m_small_streak++;
                    else                               m_small_streak = 0;
                    if (m_small_streak >= 2) {
                        if (g_wrt_env.s_OCRT_TRACE_MSPEC)
                            fprintf(stderr, "[mspec] early exit after m=%d (streak=2, tol=%.1e)\n", m, o.tolerance);
                        m_early_exit_at = m;
                        break;
                    }
                }
            }
        }

        /* 7e-bis. Cox-Munk T_wa (option A): store this Fourier mode's upwelling
         * radiance at every positive mu-node (z=0-), for ALL m, so the
         * water->air BTDF can be integrated over the full in-water field after
         * the m-loop. (Flat T_wa needs only the single view node; Cox-Munk mixes
         * many in-water directions, hence the full per-node field.) */
        for (int jp = 1; jp <= n_mu; ++jp) {
            size_t i_up = (size_t)0 * (size_t)dirs + (size_t)(+jp + n_mu);
            size_t idx  = (size_t)m * (size_t)n_mu + (size_t)(jp - 1);
            I_m_node[idx] = tot_i[i_up];
            Q_m_node[idx] = tot_q[i_up];
            U_m_node[idx] = tot_u[i_up];
            if (getenv("OCRT_D3_DIAG") && jp <= 2)
                fprintf(stderr, "[extdiag] m=%d jp=%d totI=%.3e totQ=%.3e totU=%.3e\n",
                        m, jp, tot_i[i_up], tot_q[i_up], tot_u[i_up]);
        }

        /* 7f. For hemispheric integrals (m=0 only) capture all μ at k=0
         *      and k=1 (for Kd, Ku log-derivative). */
        if (m == 0) {
            for (int jp = 1; jp <= n_mu; ++jp) {
                size_t i_up = (size_t)0 * (size_t)dirs + (size_t)(+jp + n_mu);
                size_t i_dn = (size_t)0 * (size_t)dirs + (size_t)(-jp + n_mu);
                I_m_pos[(size_t)0 * n_mu + (size_t)(jp - 1)] = tot_i[i_up];
                I_m_neg[(size_t)0 * n_mu + (size_t)(jp - 1)] = tot_i[i_dn];
                Q_m_pos[(size_t)0 * n_mu + (size_t)(jp - 1)] = tot_q[i_up];
            }
            /* k=1 level for Kd, Ku */
            for (int s = -n_mu; s <= n_mu; ++s) {
                I_m0_lvl1[s + n_mu] = tot_i[(size_t)1 * dirs + (size_t)(s + n_mu)];
            }
            /* Bottom level k=nt hemispheric diffuse fluxes.  Ed_bottom also
             * receives the direct beam later; Eu_bottom is diffuse only because
             * the current bottom boundary is black. */
            bottom_Ed_diffuse_norm = 0.0;
            bottom_Eu_diffuse_norm = 0.0;
            for (int jp = 1; jp <= n_mu; ++jp) {
                const double mu_j = atm.rm[+jp];
                const double w_j  = atm.gb[+jp];
                const size_t i_bu = (size_t)nt * (size_t)dirs + (size_t)(+jp + n_mu);
                const size_t i_bd = (size_t)nt * (size_t)dirs + (size_t)(-jp + n_mu);
                bottom_Eu_diffuse_norm += tot_i[i_bu] * mu_j * w_j;
                bottom_Ed_diffuse_norm += tot_i[i_bd] * mu_j * w_j;
            }
            bottom_Eu_diffuse_norm *= 2.0 * M_PI;
            bottom_Ed_diffuse_norm *= 2.0 * M_PI;
        }
    }

    /* 8. Fourier reconstruction in the canonical OCRT public RAA.
     *
     * The atmospheric and reported water Stokes fields share one public
     * reconstruction convention: rt_fourier_reconstruct_* receives RAA_OCRT
     * and applies its internal +pi modal shift.  A legacy pi-RAA conversion
     * mirrored the reported water branches (0 <-> 180) while the atmospheric
     * output retained the public RAA.  Keep local propagation-vector geometry
     * (rt_raa_to_water_scatter_phi) separate from this output reconstruction. */
    /* Step54 diagnostic: exact-nadir RAA is geometrically undefined.
     * OSOAA's user-angle machinery effectively fixes a canonical azimuth at
     * nadir; retaining the public RAA in the Fourier reconstruction injects
     * non-physical RAA dependence into I.  Guarded default-off handler:
     *   OCRT_WATER_NADIR_RAA_MODE=canonical
     *   OCRT_WATER_NADIR_RAA_DEG=90   (default when mode=canonical)
     * Only the water-side reconstruction azimuth is canonicalized; the public
     * metadata/provenance raa_deg is left unchanged. */
    if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] wmloop end abs=%.4f\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
    if (grid_cache_on && !grid_hit && g_grid_cache.tot) {
        g_grid_cache.valid=1; g_grid_cache.nt=nt; g_grid_cache.n_mu=n_mu; g_grid_cache.dirs=dirs;
        g_grid_cache.max_orders_seen_c=max_orders_seen; g_grid_cache.all_conv_c=all_conv; g_grid_cache.worst_resid_c=worst_resid;
        g_grid_cache.sza=sza_deg_air; g_grid_cache.lam=lambda_nm; g_grid_cache.Twc=T_water_C; g_grid_cache.Sgk=S_water_gkg;
        g_grid_cache.nw=n_water; g_grid_cache.Fs=F_sun; g_grid_cache.wind=o.wind_speed; g_grid_cache.tol=o.tolerance;
        g_grid_cache.nmw=o.n_mu_water; g_grid_cache.mmw=o.m_max_water; g_grid_cache.nlw=o.n_layers_water;
        g_grid_cache.maxit=o.max_iterations; g_grid_cache.vasn=o.view_as_node; g_grid_cache.wpk=o.water_phase_kernel; g_grid_cache.nvv=o.n_view_vza;
        g_grid_cache.a_t=o.fixed_a_total_m_inv; g_grid_cache.b_t=o.fixed_b_total_m_inv; g_grid_cache.bb_t=o.fixed_bb_total_m_inv;
        g_grid_cache.taumax=o.tau_max_target;
        g_grid_cache.a_iop=a_tot; g_grid_cache.b_iop=b_tot; g_grid_cache.bb_iop=bb_tot;
        g_grid_cache.b_phyto=b_phyto; g_grid_cache.b_detritus=b_det; g_grid_cache.b_mineral=b_min;
        g_grid_cache.constituent_mode=o.ccrr_mode;
        g_grid_cache.constituent_model=o.water_constituent_model;
        g_grid_cache.organic_group=o.organic_phyto_group;
        g_grid_cache.tsm_species=o.tsm_species;
        g_grid_cache.moment_mode=o.water_mie_moment_mode;
        g_grid_cache.moment_nmg=o.water_mie_moment_n_mu;
        g_grid_cache.phase_lmax=o.fixed_bulk_lmax;
        g_grid_cache.mie_trunc=o.water_mie_truncation_mode;
        g_grid_cache.mie_ss=o.water_mie_ss_mode;
        g_grid_cache.value_phase_spline=g_wrt_env.f_OCRT_VALUE_PHASE_SPLINE;
        g_grid_cache.value_kernel_pol=g_wrt_env.f_OCRT_WATER_VALUE_KERNEL_POL;
        snprintf(g_grid_cache.wmp,512,"%s",o.water_mie_phase_path?o.water_mie_phase_path:"");
        snprintf(g_grid_cache.fbl,512,"%s",o.fixed_bulk_phase_lut_path?o.fixed_bulk_phase_lut_path:"");
        { const char *bc = getenv("OCRT_D3_BEAM_CACHE");
          if (bc && bc[0] && !o.bypass_grid_cache) {
              FILE *tf = fopen(bc, "rb");
              if (tf) fclose(tf);
              else ocrt_beam_cache_dump(bc);
          } }
    }
    double raa_recon_deg = raa_deg;
    {
        const char *nm = g_wrt_env.s_OCRT_WATER_NADIR_RAA_MODE;
        if (nm && (!strcmp(nm, "canonical") || !strcmp(nm, "fixed") || !strcmp(nm, "1"))) {
            double nadir_eps_deg = 1.0e-7;
            const char *ne = g_wrt_env.s_OCRT_WATER_NADIR_EPS_DEG;
            if (ne && *ne) {
                double v = atof(ne);
                if (v >= 0.0 && v < 1.0) nadir_eps_deg = v;
            }
            if (fabs(vza_deg_air) <= nadir_eps_deg) {
                const char *nr = g_wrt_env.s_OCRT_WATER_NADIR_RAA_DEG;
                raa_recon_deg = (nr && *nr) ? atof(nr) : 90.0;
            }
        }
    }
    double dphi_rad = rt_raa_to_atm_fourier_phi(raa_recon_deg);
    double I_diffuse_view = rt_solver_reconstruct_phi    (I_m_view, o.m_max_water, dphi_rad);
    double Q_diffuse_view = rt_solver_reconstruct_phi    (Q_m_view, o.m_max_water, dphi_rad);
    /* U uses the same positive sine reconstruction as atmospheric output. */
    double U_diffuse_view =  rt_solver_reconstruct_phi_sin(U_m_view, o.m_max_water, dphi_rad);
    { int dv = g_wrt_env.f_OCRT_DUMP_VIEW;
      if (dv) { double iss = rt_solver_reconstruct_phi(Iss_m_view, o.m_max_water, dphi_rad);
        fprintf(stderr, "VIEW mu_w=%.8f I_diff_preTMS=%.8e Iss=%.8e ff_mode=%d ff_betalL=%d view_as_node=%d\n",
                mu_view_water, I_diffuse_view, iss, ccrr_particle_ff_mode, ccrr_ff_betal_L, o.view_as_node); } }

    /* v1.054 TMS single-scatter exact-phase correction (Nakajima-Tanaka) for the
     * internal FF particle phase.  I_diffuse_view contains the delta-M single
     * scatter built from the truncated (moment) phase, whose backscatter
     * reconstruction is Gibbs-negative for the very sharp FF phase.  Replace its
     * single-scatter part by the exact-analytic-phase single scatter:
     *   I_diffuse_view' = I_diffuse_view + Iss * (S_exact/S_deltaM - 1)
     * S = omega_w*Pw(Theta) + omega_p*Pp(Theta) is the single-scatter source factor
     * (water part identical in both, cancels; only the particle phase/SSA differ). */
    if (ccrr_particle_ff_mode && b_particle_input > 0.0 && ccrr_ff_betal_L >= 0) {
        double Iss_diffuse_view = rt_solver_reconstruct_phi(Iss_m_view, o.m_max_water, dphi_rad);
        double s0 = sqrt(fmax(0.0, 1.0 - mu_sun_water * mu_sun_water));
        double sv = sqrt(fmax(0.0, 1.0 - mu_view_water * mu_view_water));
        /* v1.06 fix: scattering angle must match the SOS field azimuth
         * convention (verified vs OSOAA/independent: --raa180 -> same-side
         * Θ≈123°, --raa0 -> back-side Θ≈167°):
         *   cosΘ = -mu_sun*mu_view - s0*sv*cos(raa)
         * (prior +s0*sv*cos(raa)-mu_sun*mu_view had the azimuth term sign
         * flipped, placing the TMS phase ratio at the wrong scattering angle;
         * latent because this block is ff_mode-only and the ratio≈1 in the
         * backscatter hemisphere where it was usually sampled). */
        double cos_th = rt_water_scatter_cos_from_public_raa(mu_sun_water, mu_view_water, raa_deg);
        if (cos_th > 1.0) cos_th = 1.0;
        if (cos_th < -1.0) cos_th = -1.0;
        const double dw = 0.039; double Kw = 3.0 / (2.0 * (2.0 + dw));
        double Pw = Kw * ((1.0 - dw) * (1.0 + cos_th * cos_th) + 2.0 * dw);
        double Pp_dm = 0.0;
        for (int l = 0; l <= ccrr_ff_betal_L; ++l)
            Pp_dm += ccrr_ff_betal[l] * fixed_bulk_legendre_P(l, cos_th);
        double Pp_ex = fixed_bulk_phase_lut_eval(&ccrr_particle_phase_table, cos_th);
        double om_w  = (ext > 0.0) ? (b_w / ext) : 0.0;
        double om_pd = (ext > 0.0) ? (b_particle_input * (1.0 - ccrr_particle_delta_f) / ext) : 0.0;
        double ext_ex = a_tot + b_w + b_particle_input;
        double om_pe = (ext_ex > 0.0) ? (b_particle_input / ext_ex) : 0.0;
        double S_dm = om_w * Pw + om_pd * Pp_dm;
        double S_ex = om_w * Pw + om_pe * Pp_ex;
        if (fabs(S_dm) > 1.0e-30) {
            double ratio = S_ex / S_dm;
            I_diffuse_view += Iss_diffuse_view * (ratio - 1.0);
        }
    }

    /* v1.06 fixed-bulk TMS-equivalent (Nakajima-Tanaka single-scatter exact-phase
     * correction) for the lut-deltam OSOAA-cap path. The CCRR TMS above is gated
     * on ccrr_particle_ff_mode (OFF for fixed-bulk), so the fixed-bulk cap+moment
     * single scatter uses the truncated/renormalized phase and over-predicts the
     * backscatter SS by ~1/(1-omega0*A) (visible at low omega0 where SS dominates).
     * Replace the SS part by the exact-phase single scatter:
     *   I += Iss * (S_ex/S_dm - 1),  S = omega0 * P(Theta)
     *   S_dm: delta-M/cap source (omega0_eff, moment-reconstructed capped phase)
     *   S_ex: exact source   (omega0_full, un-capped LUT phase, re-loaded since the
     *                         cap modifies fixed_phase_table in place)
     * Env-gated diagnostic (OCRT_FIXEDBULK_TMS); physics-changing, default-off. */
    {
        static int use_fbtms = -1;
        if (use_fbtms < 0) use_fbtms = (ocrt_debug_env("OCRT_FIXEDBULK_TMS") != NULL) ? 1 : 0;
        /* TMS applies to both the moment-kernel delta-M path (fixed_bulk_ff_betal_L>=0)
         * and the value-kernel formal delta-M path (ff_betal_L<0, truncated phase held
         * in fixed_phase_table).  In the value path the delta-M single-scatter source
         * uses the truncated/renormalized phase evaluated directly from fixed_phase_table
         * (no moment reconstruction), so S_dm = omega0_eff * P_trunc(Theta). */
        const int fbtms_moment = (fixed_bulk_ff_betal_L >= 0);
        const int fbtms_value  = (!fbtms_moment && fixed_phase_table.n >= 2);
        if (use_fbtms && o.fixed_bulk_iop_mode && (fbtms_moment || fbtms_value) && b_tot > 0.0) {
            fixed_bulk_phase_table_t uncap = { 0, NULL, NULL };
            if (fixed_bulk_phase_lut_load(o.fixed_bulk_phase_lut_path,
                                          o.fixed_bulk_phase_case_id,
                                          o.fixed_bulk_phase_wavelength_nm, &uncap) == 0) {
                double Iss_fb = rt_solver_reconstruct_phi(Iss_m_view, o.m_max_water, dphi_rad);
                double s0 = sqrt(fmax(0.0, 1.0 - mu_sun_water * mu_sun_water));
                double sv = sqrt(fmax(0.0, 1.0 - mu_view_water * mu_view_water));
                (void)s0; (void)sv;
                double cth = rt_water_scatter_cos_from_public_raa(mu_sun_water, mu_view_water, raa_deg);
                if (cth > 1.0) cth = 1.0; if (cth < -1.0) cth = -1.0;
                double P_dm = 0.0;
                if (fbtms_moment) {
                    for (int l = 0; l <= fixed_bulk_ff_betal_L; ++l)
                        P_dm += fixed_bulk_ff_betal[l] * fixed_bulk_legendre_P(l, cth);
                } else {
                    /* value-kernel: SOS consumed the truncated phase values directly */
                    P_dm = fixed_bulk_phase_lut_eval(&fixed_phase_table, cth);
                }
                double P_ex = fixed_bulk_phase_lut_eval(&uncap, cth);
                double b_rt  = b_tot * (1.0 - fixed_bulk_delta_f);
                double om_dm = (a_tot + b_rt  > 0.0) ? (b_rt  / (a_tot + b_rt )) : 0.0;
                double om_ex = (a_tot + b_tot > 0.0) ? (b_tot / (a_tot + b_tot)) : 0.0;
                double S_dm = om_dm * P_dm;
                double S_ex = om_ex * P_ex;
                if (fabs(S_dm) > 1.0e-30) {
                    double ratio = S_ex / S_dm;
                    I_diffuse_view += Iss_fb * (ratio - 1.0);
                    { int dft = g_wrt_env.f_OCRT_DUMP_VIEW;
                      if (dft) fprintf(stderr, "FBTMS[%s] cth=%.5f P_ex=%.5e P_dm=%.5e om_ex=%.6f om_dm=%.6f ratio=%.6f Iss=%.6e\n",
                                       fbtms_moment ? "moment" : "value", cth, P_ex, P_dm, om_ex, om_dm, ratio, Iss_fb); }
                }
                fixed_bulk_phase_table_free(&uncap);
            }
        }
    }
    if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] recon_done abs=%.4f\\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
    /* Step 6 diagnostic: water-side IMS/TMS for vector .mie hydrosol truncation.
     *
     * OSOAA HYD_TRUNCATION=ON is not just “truncated PM + b_eff”; it also
     * keeps separate raw/truncated source bookkeeping.  This diagnostic tests
     * the IMS hypothesis on the water side by replacing only the view-direction
     * scalar single-scatter contribution:
     *
     *     I += I_ss,tr * ( S_raw / S_tr - 1 )
     *
     * with S evaluated for a homogeneous fixed-depth column:
     *
     *     S = b * P11(theta) * G(ext,z)
     *     G = (1-exp[-ext*z*(1/mu0+1/muv)])/(ext*(1/mu0+1/muv)).
     *
     * This is diagnostic/reference-parity only.  It is restricted to fixed-bulk
     * hydrosol-only .mie diagnostic runs. NT-TMS is handled earlier by building the primary source from raw moments. */
    if (o.fixed_bulk_iop_mode && o.water_mie_ss_mode == 1 &&
        o.water_mie_truncation_mode == 1 && have_aip &&
        o.water_mie_phase_path && o.water_mie_phase_path[0] &&
        b_tot > 0.0 && b_rt > 0.0 && b_w < 1.0e-12 && particle_Lmax >= 0) {
        double Iss_tr = rt_solver_reconstruct_phi(Iss_m_view, o.m_max_water, dphi_rad);
        double cth = rt_water_scatter_cos_from_public_raa(mu_sun_water, mu_view_water, raa_deg);
        if (cth > 1.0) cth = 1.0;
        if (cth < -1.0) cth = -1.0;
        double theta_deg = acos(cth) * 180.0 / M_PI;
        double P11_raw = 0.0, P12_raw = 0.0, P33_raw = 0.0;
        eval_aerosol_phase(&aip, theta_deg, &P11_raw, &P12_raw, &P33_raw);
        (void)P12_raw; (void)P33_raw;
        double P11_tr = 0.0;
        for (int l = 0; l <= particle_Lmax; ++l)
            P11_tr += particle_betal[l] * fixed_bulk_legendre_P(l, cth);
        double ext_tr  = a_tot + b_rt;
        double ext_raw = a_tot + b_tot;
        double invmu = 1.0/fmax(mu_sun_water, 1e-12) + 1.0/fmax(mu_view_water, 1e-12);
        double G_tr = (ext_tr > 0.0 && z_max > 0.0)
                      ? (1.0 - exp(-ext_tr*z_max*invmu)) / fmax(ext_tr*invmu, 1e-300)
                      : 0.0;
        double G_raw = (ext_raw > 0.0 && z_max > 0.0)
                       ? (1.0 - exp(-ext_raw*z_max*invmu)) / fmax(ext_raw*invmu, 1e-300)
                       : 0.0;
        double S_tr  = b_rt  * P11_tr  * G_tr;
        double S_raw = b_tot * P11_raw * G_raw;
        if (isfinite(S_tr) && isfinite(S_raw) && isfinite(Iss_tr) && fabs(S_tr) > 1.0e-300) {
            double ratio = S_raw / S_tr;
            I_diffuse_view += Iss_tr * (ratio - 1.0);
            { int dwims = g_wrt_env.f_OCRT_DUMP_WATER_IMS;
              if (dwims) fprintf(stderr,
                  "WATER_IMS theta=%.4f cth=%.6f Praw=%.6e Ptr=%.6e braw=%.6e btr=%.6e ext_raw=%.6e ext_tr=%.6e z=%.6e Graw=%.6e Gtr=%.6e ratio=%.6f Iss=%.8e dI=%.8e\n",
                  theta_deg, cth, P11_raw, P11_tr, b_tot, b_rt, ext_raw, ext_tr, z_max, G_raw, G_tr, ratio, Iss_tr, Iss_tr*(ratio-1.0)); }
        }
    }

    { int dv2 = g_wrt_env.f_OCRT_DUMP_VIEW;
      if (dv2) fprintf(stderr, "VIEW I_diff_postTMS=%.8e\n", I_diffuse_view); }

    /* v1.06 polarized single-scatter exact-phase calculation (.mie path).
     * The CCRR scalar single-scatter (TMS) above is gated on
     * ccrr_particle_ff_mode, which is OFF for the fixed-bulk + .mie path. There
     * the view-direction single-scatter Stokes is the linear-interpolated
     * n_mu-quadrature value; I is smooth (interp robust) but the polarized Q/U
     * under-resolve the P12/P33 angular structure. Here the single-scatter Q/U
     * are instead obtained analytically (closed form): the single-scatter Stokes
     * factor as L_ss = [scalar geometry G(mu0,muv,tau)] x [Z(Theta)·F_beam], with
     * G common to I/Q/U. G is read off the n_mu-robust single-scatter I (so I is
     * left unchanged), and Z·F_beam is evaluated exactly (eval_aerosol_phase
     * PCHIP at the exact scattering angle + meridian-plane rotations). No free
     * parameter and no tuning to any reference. Stage 2: general geometry
     * including off-plane via the validated Hovenier rotation (Lrot, sigma1/2);
     * principal plane reduces to sigma=0/pi automatically. Env
     * OCRT_SSPOL_CORR_OFF disables for A/B checks. */
    { static int sspol_off = -1; if (sspol_off < 0) sspol_off = (ocrt_debug_env("OCRT_SSPOL_CORR_OFF") != NULL) ? 1 : 0;
      int sspol_allowed = 1;
      if (o.water_mie_moment_mode == 1) {
          /* Gauss-projected .mie moment mode is a PM-moment parity path.  The
           * legacy exact single-scatter Q/U correction evaluates the raw .mie
           * phase directly at the target angle, bypassing the Gauss moment
           * representation and producing a mixed convention.  Keep it off in
           * this mode unless explicitly re-enabled for diagnostics. */
          sspol_allowed = (ocrt_debug_env("OCRT_SSPOL_CORR_ON") != NULL);
      }
      if (!sspol_off && sspol_allowed && !have_particle_value_phase &&
          have_aip && !ccrr_particle_ff_mode && b_particle > 0.0 && b_w < 1.0e-12) {
        double s_sun  = sqrt(fmax(0.0, 1.0 - mu_sun_water  * mu_sun_water));
        double s_view = sqrt(fmax(0.0, 1.0 - mu_view_water * mu_view_water));
        /* propagation directions: sun O0=(s_sun,0,-mu_sun) downward at phi=0;
         * view Ov upward at physical azimuth phi_v = (180 - raa_deg) deg.
         * (verified: --raa180 -> phi_v=0 same-side, --raa90 -> phi_v=90,
         *  --raa0 -> phi_v=180 back-side.) */
        double phv  = rt_raa_to_water_scatter_phi(raa_deg);
        double O0[3] = { s_sun, 0.0, -mu_sun_water };
        double Ov[3] = { s_view * cos(phv), s_view * sin(phv), mu_view_water };
        double cosT = O0[0]*Ov[0] + O0[1]*Ov[1] + O0[2]*Ov[2];
        if (cosT >  1.0) cosT =  1.0;
        if (cosT < -1.0) cosT = -1.0;
        double Th = acos(cosT) * 180.0 / M_PI;
        double Pp11 = 0.0, Pp12 = 0.0, Pp33 = 0.0;
        eval_aerosol_phase(&aip, Th, &Pp11, &Pp12, &Pp33);
        double bq = atm.beam_q;                  /* Fresnel-transmitted beam Q (sun meridian) */
        /* meridian-plane rotation angles (Hovenier vector convention), ported
         * from the validated independent single-scatter:
         *   eperp(O) = norm(zhat x O) = norm(-Oy, Ox, 0)
         *   n = norm(O0 x Ov);  sig = atan2( (eperp x n)·O, eperp·n )
         *   Z = Lrot(-sig2)·F(Theta)·Lrot(sig1),  Lrot(a)=[[1,0,0],[0,c,s],[0,-s,c]]. */
        double c2s1 = 1.0, s2s1 = 0.0, c2s2 = 1.0, s2s2 = 0.0;  /* default: no rotation */
        double nvec[3] = { O0[1]*Ov[2]-O0[2]*Ov[1], O0[2]*Ov[0]-O0[0]*Ov[2], O0[0]*Ov[1]-O0[1]*Ov[0] };
        double nn = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);
        if (nn > 1.0e-9) {                       /* well-defined scattering plane */
            nvec[0]/=nn; nvec[1]/=nn; nvec[2]/=nn;
            double e0[3] = { -O0[1], O0[0], 0.0 };
            double ne0 = sqrt(e0[0]*e0[0]+e0[1]*e0[1]+e0[2]*e0[2]);
            if (ne0 > 1.0e-12) { e0[0]/=ne0; e0[1]/=ne0; e0[2]/=ne0; } else { e0[0]=0.0; e0[1]=1.0; e0[2]=0.0; }
            double ev[3] = { -Ov[1], Ov[0], 0.0 };
            double nev = sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
            if (nev > 1.0e-12) { ev[0]/=nev; ev[1]/=nev; ev[2]/=nev; } else { ev[0]=0.0; ev[1]=1.0; ev[2]=0.0; }
            double cr0[3] = { e0[1]*nvec[2]-e0[2]*nvec[1], e0[2]*nvec[0]-e0[0]*nvec[2], e0[0]*nvec[1]-e0[1]*nvec[0] };
            double sig1 = atan2(cr0[0]*O0[0]+cr0[1]*O0[1]+cr0[2]*O0[2], e0[0]*nvec[0]+e0[1]*nvec[1]+e0[2]*nvec[2]);
            double crv[3] = { ev[1]*nvec[2]-ev[2]*nvec[1], ev[2]*nvec[0]-ev[0]*nvec[2], ev[0]*nvec[1]-ev[1]*nvec[0] };
            double sig2 = atan2(crv[0]*Ov[0]+crv[1]*Ov[1]+crv[2]*Ov[2], ev[0]*nvec[0]+ev[1]*nvec[1]+ev[2]*nvec[2]);
            c2s1 = cos(2.0*sig1); s2s1 = sin(2.0*sig1);
            c2s2 = cos(2.0*sig2); s2s2 = sin(2.0*sig2);
        }
        /* beam [1,bq,0] (sun meridian) -> scattering frame via Lrot(sig1) */
        double Qb_s = bq * c2s1;
        double Ub_s = -bq * s2s1;
        /* scattering-frame Stokes after F(Theta) (P22=P11) */
        double If = Pp11 + Pp12 * Qb_s;
        double Qf = Pp12 + Pp11 * Qb_s;
        double Uf = Pp33 * Ub_s;
        double Iss_int = rt_solver_reconstruct_phi(Iss_m_view, o.m_max_water, dphi_rad);
        if (fabs(If) > 1.0e-30) {
            double C = Iss_int / If;             /* geometry amplitude (I unchanged) */
            /* rotate scattering frame -> view meridian via Lrot(-sig2) */
            double Qss_ex = C * (c2s2 * Qf - s2s2 * Uf);
            double Uss_ex = C * (s2s2 * Qf + c2s2 * Uf);
            double Qss_int =  rt_solver_reconstruct_phi    (Qss_m_view, o.m_max_water, dphi_rad);
            double Uss_int =  rt_solver_reconstruct_phi_sin(Uss_m_view, o.m_max_water, dphi_rad);
            Q_diffuse_view += (Qss_ex - Qss_int);
            U_diffuse_view += (Uss_ex - Uss_int);
            { int dsp = g_wrt_env.f_OCRT_DUMP_SSPOL;
              if (dsp) fprintf(stderr, "SSPOL Theta=%.3f Pp11=%.4e Pp12=%.4e Pp33=%.4e bq=%.5f c2s1=%.4f c2s2=%.4f Iss=%.4e dQ=%.4e dU=%.4e\n",
                               Th, Pp11, Pp12, Pp33, bq, c2s1, c2s2, Iss_int, Qss_ex - Qss_int, Uss_ex - Uss_int); }
        }
      }
    }

    if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] diag_done abs=%.4f\\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
    /* 9. Hemispheric integrals from m=0 mode (DIFFUSE component only).
     *     Ed_diffuse = ∫_{μ<0} I^{m=0}(μ) · |μ| · 2π dμ
     *     Eu_diffuse = ∫_{μ>0} I^{m=0}(μ) ·  μ  · 2π dμ */
    double Ed_diffuse = 0.0, Eu_diffuse = 0.0;
    double Ed_internal_reflect = 0.0;
    for (int jp = 1; jp <= n_mu; ++jp) {
        double mu_j = atm.rm[+jp];     /* positive Gauss node */
        double w_j  = atm.gb[+jp];
        double I_up = I_m_pos[(size_t)(jp - 1)];
        double Q_up = Q_m_pos[(size_t)(jp - 1)];
        /* Diagnostic (OCRT_DUMP_UPWELLING): m=0 in-water upwelling radiance at
         * native GL node mu_j (= in-water cos(theta_w)). Azimuthally-averaged,
         * so it equals the full BRDF only at SZA=0 (m>=1 vanish). NOTE: this is
         * the delta-M field WITHOUT the view-path TMS exact-phase single-scatter
         * correction, so near-backscatter (near-nadir) may be Gibbs-affected.
         * Env-gated, no physics change; remove before release like OCRT_DUMP_SKY. */
        { int dmp_up = g_wrt_env.f_OCRT_DUMP_UPWELLING;
          if (dmp_up) fprintf(stderr, "UPWELL %.10f %.10e\n", mu_j, I_up); }
        Eu_diffuse += I_up * mu_j * w_j;
        Ed_diffuse += I_m_neg[(size_t)(jp - 1)] * mu_j * w_j;

        /* Water-side internal Fresnel reflection at the interface.
         * The SOS volume solution currently has no top-boundary reflection
         * operator, so I_m_neg(k=0) misses the downward diffuse irradiance
         * produced by reflecting the upwelling water-side radiance back into
         * the water.  For the diagnostic Ed(0-) and rrs denominator, include
         * the m=0 reflected Stokes-I contribution explicitly.
         * [2026-05-31] FIX-C: for wind>0 use the Cox-Munk slope-integrated
         * internal reflectance (rule: no flat-Fresnel at wind>0).  Scalar (I)
         * only — the I->Q internal-reflection polarization coupling is a tiny
         * second-order effect on this m=0 Ed bookkeeping.  Flat Mueller (incl.
         * the I->Q term) is kept only at wind<=0.  Still does not feed back
         * into the SOS source field. */
        double I_ref;
        if (o.wind_speed > 0.0) {
            double Rww = surface_R_ww_coxmunk_direct(mu_j, n_water,
                                                     o.wind_speed,
                                                     o.cox_munk_sigma_type);
            I_ref = Rww * I_up;
        } else {
            double M_Rww[9];
            rt_air_water_R_ww(mu_j, n_water, o.q_convention, M_Rww);
            I_ref = M_Rww[0*3+0] * I_up + M_Rww[0*3+1] * Q_up;
        }
        Ed_internal_reflect += I_ref * mu_j * w_j;
    }
    Eu_diffuse *= 2.0 * M_PI;
    Ed_diffuse *= 2.0 * M_PI;
    Ed_internal_reflect *= 2.0 * M_PI;

    /* 9b. DIRECT BEAM contribution
     *
     *   OCRT/6SV uses F_sun=π convention internally; SOS output I_diffuse is
     *   already in *F_sun=π normalized radiance* form. To return SI-scaled
     *   values (W·m⁻²·nm⁻¹), we apply:
     *     L_SI       = I_diffuse_output × F_sun / π
     *     Ed_diffuse = (hemispheric integral) × F_sun / π
     *
     *   For E_d(0-, direct):  표준 (Mobley 1994 Eq. 4.41)
     *     E_d(0-, direct) = T_F_irrad · E_d(0+, direct) = T_F_irrad · F_sun · μ_sun_air
     *   여기서 T_F_irrad 는 *irradiance* Fresnel transmittance (unpolarized 평균).
     *
     *   이전 (pre-2026-05-23 KST2330) 식: Ed_direct = F_sun × μ_sun_water — 비물리적.
     *     ratio E_d(0-)/E_d(0+) = μ_sun_water/μ_sun_air > 1 (잘못)
     *     예: SZA=30° → 1.071 vs 정확 0.978
     *
     *   For radiance Lu_direct: a δ-function in the -μ_sun_w direction;
     *   contributes to a finite-acceptance view only if numerically aligned
     *   — assumed zero for arbitrary view geometry.
     */
    /* [v1.04+xsec_norm_cli+2, 2026-05-24] A fix 적용 layer 변경:
     * 이전 (v1.04+qaa_debug, 2026-05-23 KST2330): 본 위치에서 Ed_direct 식만 수정
     *   Ed_direct = T_F · F_sun · μ_sun_air
     * 문제: F_sun 입력은 air-side BOA value 그대로라, Lu source normalization 에는
     *       Fresnel transmission 이 반영 안 됨 → Lu(0-) ~9.55 % over 발현.
     *
     * 현재 (옵션 A, 첨부 OCRT_water_RT_QAA_bughunt_..._2026-05-24 §2.1 권장):
     *   rt_solver.c::rt_solve_case_ocean() 에서 water solver 에 들어가는
     *   F_sun 자체를 F_sun_water_direct = F_sun_BOA·μ_air·T_aw/μ_water 로 변환.
     *   → 본 위치는 F_sun_water 가 이미 변환된 값이므로 원래의 표준 식
     *     Ed_direct = F_sun · μ_sun_water 로 복원하면, Ed_direct(0-) 가
     *     자동으로 정확한 T_F·F_sun_BOA·μ_air = 2.660304 W/m²/nm 가 된다.
     *   Lu source 도 변환된 F_sun_water 로 normalize 되어 일관성 확보.
     *
     * 따라서 본 코드는 첨부 §2.1 fix 와 함께 atomic 하게 동작. 둘 중 하나만
     * 적용하면 잘못된 결과 (double Fresnel 또는 Lu over).
     */
    double Ed_direct      = F_sun * mu_sun_water;
    double Eu_direct      = 0.0;
    double Lu_direct_view = 0.0;

    /* 10. Convert diffuse outputs from π-normalized to SI */
    /* P0-A: when the surface internal reflection is fed back into the SOS
     * (Rww_M set), the reflected downward field is already present in
     * Ed_diffuse (I_m_neg at k=0). Adding the post-hoc Ed_internal_reflect
     * again would double-count it in the rrs denominator. */
    if (Rww_M || Rww_K) Ed_internal_reflect = 0.0;   /* B2: rough feedback also in source */
    const double f_scale = F_sun / M_PI;
    double Ed_diffuse_SI = (Ed_diffuse + Ed_internal_reflect) * f_scale;
    double Eu_diffuse_SI = Eu_diffuse * f_scale;

    /* 11. Total Ed, Eu, Lu (direct + diffuse, SI scale) */
    double Ed_total = Ed_direct + Ed_diffuse_SI;
    double Eu_total = Eu_direct + Eu_diffuse_SI;
    { int di = g_wrt_env.f_OCRT_DUMP_IOP;
      if (di) { double IR = Ed_internal_reflect * f_scale, ED = Ed_diffuse * f_scale, EU = Eu_diffuse * f_scale;
        fprintf(stderr, "EDCOMP Ed_direct=%.6e Ed_diffuse=%.6e Ed_intrefl=%.6e Eu_diffuse=%.6e | intrefl/Ed_total=%.4f Eu/Ed=%.4f\n",
                Ed_direct, ED, IR, EU, IR/Ed_total, EU/Ed_total); } }
    double I_view   = Lu_direct_view + I_diffuse_view * f_scale;
    double Q_view   = Q_diffuse_view * f_scale;
    double U_view   = U_diffuse_view * f_scale;
    /* 0- (in-water, just below interface) view Stokes dump for the surface-
     * transmission diagnostic. Env-gated, stderr only; no default-output or
     * physics change. */
    { int d0m = g_wrt_env.f_OCRT_DUMP_0MINUS;
      if (d0m) fprintf(stderr, "STOKES0MINUS Lu0minus=%.10e Qu0minus=%.10e Uu0minus=%.10e\n",
                       I_view, Q_view, U_view); }

    /* Diagnostic-only bottom flux audit.  A physically deep black-bottom proxy
     * should make these ratios small.  If they remain large, the finite slab is
     * participating in the solution and tau/depth is not a valid convergence
     * setting. */
    {
        const char *bfp = g_wrt_env.s_OCRT_DUMP_WATER_BOTTOM_FLUX;
        if (bfp && bfp[0]) {
            double Ed_bottom_diffuse_SI = bottom_Ed_diffuse_norm * f_scale;
            double Eu_bottom_diffuse_SI = bottom_Eu_diffuse_norm * f_scale;
            double Ed_bottom_direct_SI = F_sun * mu_sun_water * exp(-tau_max / fmax(mu_sun_water, 1e-12));
            double Ed_bottom_total_SI = Ed_bottom_direct_SI + Ed_bottom_diffuse_SI;
            FILE *bf = fopen(bfp, "a");
            if (bf) {
                fseek(bf, 0, SEEK_END);
                if (ftell(bf) == 0L) {
                    fprintf(bf, "wavelength_nm,sza_deg,vza_deg,raa_deg,tau_max,z_max,n_layers,dtau_layer,a_total,b_total,omega,orders,conv,Ed_surface,Eu_surface,Lu_view_surface,Ed_bottom_total,Ed_bottom_direct,Ed_bottom_diffuse,Eu_bottom,Ed_bottom_over_Ed_surface,Eu_bottom_over_Eu_surface,Eu_bottom_over_Ed_surface\n");
                }
                fprintf(bf, "%.10g,%.10g,%.10g,%.10g,%.12e,%.12e,%d,%.12e,%.12e,%.12e,%.12e,%d,%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                        lambda_nm, sza_deg_air, vza_deg_air, raa_deg,
                        tau_max, z_max, n_layers_water_eff, tau_max/(double)n_layers_water_eff,
                        a_tot, b_tot, omega, max_orders_seen, all_conv,
                        Ed_total, Eu_total, I_view, Ed_bottom_total_SI,
                        Ed_bottom_direct_SI, Ed_bottom_diffuse_SI, Eu_bottom_diffuse_SI,
                        (Ed_total != 0.0 ? Ed_bottom_total_SI/Ed_total : 0.0),
                        (Eu_total != 0.0 ? Eu_bottom_diffuse_SI/Eu_total : 0.0),
                        (Ed_total != 0.0 ? Eu_bottom_diffuse_SI/Ed_total : 0.0));
                fclose(bf);
            }
        }
    }

    if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] hemi_done abs=%.4f\\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
    /* 12. Kd, Ku from log-derivative across [k=0, k=1]  (SI scale)
     *
     * Definitions (diffuse attenuation coefficients, m^-1):
     *   Kd = -d ln Ed / dz |_{z=0}  ~=  -ln( Ed(z1)/Ed(0) ) / z1
     *   Ku = -d ln Eu / dz |_{z=0}  ~=  -ln( Eu(z1)/Eu(0) ) / z1
     * evaluated by finite difference between the surface level (k=0) and the
     * first optical-depth level k=1 (z1 = dtau_layer / ext).
     *
     * Plane irradiances from the m=0 azimuthal mode only (the m>0 modes
     * integrate to zero over azimuth):
     *   Ed(level) = Ed_direct + 2*pi * sum_{j>0} I^0(level, -mu_j) mu_j w_j
     *   Eu(level) =             2*pi * sum_{j>0} I^0(level, +mu_j) mu_j w_j
     *   Ed_direct(level) = F_sun * mu_s,w * exp(-tau_level / mu_s,w)
     * with F_sun already the BELOW-surface beam (air->water Fresnel + n^2
     * applied upstream), mu_s,w the in-water solar cosine, and the GL sum
     * the discrete 2*pi Int I mu dmu over the respective hemisphere.
     * Sign convention of the Fourier storage: index -j = downwelling,
     * +j = upwelling.  Ku from a one-layer difference is noisy when Eu is
     * near-constant (can even go negative near the surface for upwelling
     * light growing with depth toward the asymptotic regime) — it is a
     * DIAGNOSTIC, not a production product. */
    double Ed_diff_lvl1 = 0.0, Eu_diff_lvl1 = 0.0;
    for (int jp = 1; jp <= n_mu; ++jp) {
        double mu_j = atm.rm[+jp];
        double w_j  = atm.gb[+jp];
        Ed_diff_lvl1 += I_m0_lvl1[-jp + n_mu] * mu_j * w_j;
        Eu_diff_lvl1 += I_m0_lvl1[+jp + n_mu] * mu_j * w_j;
    }
    Ed_diff_lvl1 *= 2.0 * M_PI * f_scale;
    Eu_diff_lvl1 *= 2.0 * M_PI * f_scale;
    double tau_lvl1 = tau_max / (double)nt;
    /* 옵션 A: F_sun 은 이미 air-water Fresnel 변환된 값. 표준 식 사용. */
    double Ed_direct_lvl1 = F_sun * mu_sun_water * exp(-tau_lvl1 / mu_sun_water);
    double Ed_total_lvl1 = Ed_direct_lvl1 + Ed_diff_lvl1;
    double Eu_total_lvl1 = Eu_diff_lvl1;
    double z_lvl1 = tau_lvl1 / ext;
    double Kd = (Ed_total_lvl1 > 0.0 && Ed_total > 0.0) ?
                -log(Ed_total_lvl1 / Ed_total) / fmax(z_lvl1, 1e-12) : 0.0;
    double Ku = (Eu_total_lvl1 > 0.0 && Eu_total > 0.0) ?
                -log(Eu_total_lvl1 / Eu_total) / fmax(z_lvl1, 1e-12) : 0.0;

    if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] m13 abs=%.4f\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
    /* 13. Above-water Stokes via T_wa (water → air) radiance Mueller transform.
     *
     *   L_air(μ_view_air) = T_wa_radiance · L_water(μ_view_water)
     *
     *   T_wa_radiance is computed at the in-water view direction. It already
     *   includes the (1/n_w)² radiance reduction (n² law inverse for upward
     *   transmission).
     *
     *   Note (B.3): T_wa applied to *diffuse* in-water Stokes only.
     *   Direct beam in-water → air is the *original solar beam* (unaffected
     *   by upward transmission since direct is downward-going at solar
     *   direction); does NOT contribute to upward L_u(0+, view).
     *
     *   짚어둘 점: TIR check 불필요 — μ_view_water from Snell-down of μ_view_air
     *   is always above mu_critical (downward refraction has no TIR), so the
     *   reverse direction (upward water→air at same μ_water) is also above
     *   critical and transmits cleanly.
     */
    double I_air, Q_air, U_air;
    /* 1B TIR guard: a super-critical in-water view direction (mu_view_air==0,
     * set above when OCRT_VZA_IN_WATER and theta_w>theta_c) has no above-water
     * radiance — total internal reflection. Zero the 0+ Stokes and skip T_wa. */
    if (mu_view_air <= 0.0) {
        I_air = 0.0; Q_air = 0.0; U_air = 0.0;
        { int dtwap = g_wrt_env.f_OCRT_DUMP_TWA_PROVENANCE;
          if (dtwap) fprintf(stderr,
              "TWA_PROV wl=%.3f sza=%.6f vza=%.6f raa=%.6f branch=TIR_or_invalid mu_air=%.12g mu_w=%.12g n=%.12g I0m=%.12e Q0m=%.12e U0m=%.12e Iair=0 Qair=0 Uair=0 Ed0m=%.12e Ed0p=%.12e rrs0m=%.12e Rrs0p=0\n",
              lambda_nm, sza_deg_air, vza_deg_air, raa_deg, mu_view_air, mu_view_water, n_water,
              I_view, Q_view, U_view, Ed_total, F_sun * mu_sun_air,
              (Ed_total > 0.0) ? I_view / Ed_total : 0.0); }
    } else if (o.wind_speed <= 0.0) {
        int twa_operator_mode = ocrt_water_twa_operator_mode_for_vza(vza_deg_air);
        if (twa_operator_mode != 0) {
            /* OSOAA-compatible flat interface discretisation:
             *   (1) reconstruct L_w at every water-side quadrature node for the
             *       requested azimuth,
             *   (2) apply T_wa(mu_w_node) there,
             *   (3) map each node to air-side mu_a_node,
             *   (4) OSOAA-compatible cubic-spline interpolate the transmitted field to
             *       the requested air mu.
             *
             * This deliberately changes only the output operator order; the
             * underwater SOS field itself is unchanged. */
            double *mu_a_nodes = (double*)calloc((size_t)(n_mu + 8), sizeof(double));
            double *mu_w_nodes = (double*)calloc((size_t)(n_mu + 8), sizeof(double));
            double *ya_I = (double*)calloc((size_t)(n_mu + 8), sizeof(double));
            double *ya_Q = (double*)calloc((size_t)(n_mu + 8), sizeof(double));
            double *ya_U = (double*)calloc((size_t)(n_mu + 8), sizeof(double));
            double *lw_I_nodes = (double*)calloc((size_t)(n_mu + 8), sizeof(double));
            double *lw_Q_nodes = (double*)calloc((size_t)(n_mu + 8), sizeof(double));
            double *lw_U_nodes = (double*)calloc((size_t)(n_mu + 8), sizeof(double));
            int *node_src = (int*)calloc((size_t)(n_mu + 8), sizeof(int));
            int nv = 0;
            const int Mc = m_loop_max + 1;
            for (int jp = 1; jp <= n_mu; ++jp) {
                const double mu_wj = atm.rm[+jp];
                const double mu_aj = rt_air_water_mu_refracted_up(mu_wj, n_water);
                if (mu_aj < 0.0) continue; /* TIR node; not part of transmitted air grid */
                double pmI[64], pmQ[64], pmU[64];
                if (Mc > 64) {
                    /* Current water Fourier caps are far below this.  If a future
                     * run exceeds it, skip the operator mode rather than risk an
                     * overflow. */
                    nv = 0;
                    break;
                }
                for (int mm = 0; mm < Mc; ++mm) {
                    size_t idx = (size_t)mm*(size_t)n_mu + (size_t)(jp - 1);
                    pmI[mm] = I_m_node[idx];
                    pmQ[mm] = Q_m_node[idx];
                    pmU[mm] = U_m_node[idx];
                }
                double Lw_I = rt_solver_reconstruct_phi(pmI, Mc - 1, dphi_rad) * f_scale;
                double Lw_Q = rt_solver_reconstruct_phi(pmQ, Mc - 1, dphi_rad) * f_scale;
                double Lw_U =  rt_solver_reconstruct_phi_sin(pmU, Mc - 1, dphi_rad) * f_scale;
                double M_T_wa_node[9];
                rt_air_water_T_wa(mu_wj, n_water, o.q_convention, M_T_wa_node);
                mu_w_nodes[nv] = mu_wj;
                mu_a_nodes[nv] = mu_aj;
                lw_I_nodes[nv] = Lw_I;
                lw_Q_nodes[nv] = Lw_Q;
                lw_U_nodes[nv] = Lw_U;
                node_src[nv] = 0;
                ya_I[nv] = M_T_wa_node[0*3+0]*Lw_I + M_T_wa_node[0*3+1]*Lw_Q + M_T_wa_node[0*3+2]*Lw_U;
                ya_Q[nv] = M_T_wa_node[1*3+0]*Lw_I + M_T_wa_node[1*3+1]*Lw_Q + M_T_wa_node[1*3+2]*Lw_U;
                ya_U[nv] = M_T_wa_node[2*3+0]*Lw_I + M_T_wa_node[2*3+1]*Lw_Q + M_T_wa_node[2*3+2]*Lw_U;
                nv++;
            }
            if (ocrt_water_twa_spline_node_mode() == 1 && nv < n_mu + 8) {
                /* Step56 diagnostic: OSOAA's RMU array contains three zero-weight
                 * non-quadrature slots for the flat-interface operator: mu=1,
                 * the refracted solar water cosine, and the air-side solar cosine.
                 * They do not carry quadrature weight, but OSOAA_INTERF_MERPLATE
                 * applies TWA at every RMU slot and then calls SOS_INTERPO_SPLINT
                 * on the resulting RMUT/I_DEV arrays.  Insert the same slots in
                 * the OCRT operator-first spline input, reconstructing the water
                 * field by the same linear-in-mu interpolation used for target
                 * extraction. */
                double extra_mu_w[3];
                extra_mu_w[0] = 1.0;
                extra_mu_w[1] = mu_sun_water;
                extra_mu_w[2] = mu_sun_air;
                for (int ex = 0; ex < 3 && nv < n_mu + 8; ++ex) {
                    double mu_wx = extra_mu_w[ex];
                    if (!(mu_wx > 0.0 && mu_wx <= 1.0)) continue;
                    double mu_ax = rt_air_water_mu_refracted_up(mu_wx, n_water);
                    if (mu_ax < 0.0) continue;
                    int dup = 0;
                    for (int ii = 0; ii < nv; ++ii) {
                        if (fabs(mu_a_nodes[ii] - mu_ax) < 1e-12) { dup = 1; break; }
                    }
                    if (dup) continue;
                    int jlo = 1, jhi = 2;
                    if (mu_wx <= atm.rm[+1]) {
                        jlo = 1; jhi = 2;
                    } else if (mu_wx >= atm.rm[+n_mu]) {
                        jlo = n_mu - 1; jhi = n_mu;
                    } else {
                        for (int jp2 = 1; jp2 < n_mu; ++jp2) {
                            if (atm.rm[+jp2] <= mu_wx && mu_wx < atm.rm[+jp2+1]) {
                                jlo = jp2; jhi = jp2 + 1; break;
                            }
                        }
                    }
                    double mlo = atm.rm[+jlo], mhi = atm.rm[+jhi];
                    double den = mhi - mlo;
                    double whi = (den != 0.0) ? (mu_wx - mlo) / den : 0.0;
                    double wlo = 1.0 - whi;
                    double pmI[64], pmQ[64], pmU[64];
                    if (Mc > 64) break;
                    for (int mm = 0; mm < Mc; ++mm) {
                        size_t il = (size_t)mm*(size_t)n_mu + (size_t)(jlo - 1);
                        size_t ih = (size_t)mm*(size_t)n_mu + (size_t)(jhi - 1);
                        pmI[mm] = wlo * I_m_node[il] + whi * I_m_node[ih];
                        pmQ[mm] = wlo * Q_m_node[il] + whi * Q_m_node[ih];
                        pmU[mm] = wlo * U_m_node[il] + whi * U_m_node[ih];
                    }
                    double Lw_I = rt_solver_reconstruct_phi(pmI, Mc - 1, dphi_rad) * f_scale;
                    double Lw_Q = rt_solver_reconstruct_phi(pmQ, Mc - 1, dphi_rad) * f_scale;
                    double Lw_U =  rt_solver_reconstruct_phi_sin(pmU, Mc - 1, dphi_rad) * f_scale;
                    double M_T_wa_node[9];
                    rt_air_water_T_wa(mu_wx, n_water, o.q_convention, M_T_wa_node);
                    mu_a_nodes[nv] = mu_ax;
                    ya_I[nv] = M_T_wa_node[0*3+0]*Lw_I + M_T_wa_node[0*3+1]*Lw_Q + M_T_wa_node[0*3+2]*Lw_U;
                    ya_Q[nv] = M_T_wa_node[1*3+0]*Lw_I + M_T_wa_node[1*3+1]*Lw_Q + M_T_wa_node[1*3+2]*Lw_U;
                    ya_U[nv] = M_T_wa_node[2*3+0]*Lw_I + M_T_wa_node[2*3+1]*Lw_Q + M_T_wa_node[2*3+2]*Lw_U;
                    nv++;
                }
            }
            if (g_wrt_env.s_OCRT_DUMP_TWA_SPLINE_NODES != NULL && nv > 0) {
                const char *path = g_wrt_env.s_OCRT_DUMP_TWA_SPLINE_NODES;
                FILE *fp = fopen((path && path[0]) ? path : "ocrt_twa_spline_nodes.csv", "a");
                if (fp) {
                    static int hdr = 0;
                    if (!hdr) {
                        fprintf(fp, "lambda_nm,sza_deg,vza_deg,raa_deg,node_mode,op_mode,ep_mode,node_rank,mu_air_node,I_air_node,Q_air_node,U_air_node,target_mu_air,target_mu_w,Iair_final,Qair_final,Uair_final\n");
                        hdr = 1;
                    }
                    for (int ii = 0; ii < nv; ++ii) {
                        fprintf(fp, "%.10g,%.10g,%.10g,%.10g,%d,%d,%d,%d,%.17g,%.17e,%.17e,%.17e,%.17g,%.17g,%.17e,%.17e,%.17e\n",
                                lambda_nm, sza_deg_air, vza_deg_air, raa_deg,
                                ocrt_water_twa_spline_node_mode(), twa_operator_mode,
                                ocrt_water_twa_critical_endpoint_mode(), ii,
                                mu_a_nodes[ii], ya_I[ii], ya_Q[ii], ya_U[ii],
                                mu_view_air, mu_view_water, 0.0, 0.0, 0.0);
                    }
                    fclose(fp);
                }
            }
            if (ocrt_water_twa_nodeset_mode() == 1 && nv > 0) {
                const double extra_mus[3] = {1.0, mu_sun_air, mu_sun_water};
                for (int iex = 0; iex < 3; ++iex) {
                    const double mu_wx = extra_mus[iex];
                    if (!(mu_wx > 0.0 && mu_wx <= 1.0)) continue;
                    int dup = 0;
                    for (int iq = 0; iq < nv; ++iq) if (fabs(mu_w_nodes[iq] - mu_wx) < 1e-10) { dup = 1; break; }
                    if (dup) continue;
                    const double mu_ax = rt_air_water_mu_refracted_up(mu_wx, n_water);
                    if (mu_ax < 0.0) continue;
                    if (nv >= n_mu + 8) break;
                    double pmI[64], pmQ[64], pmU[64];
                    if (Mc > 64) break;
                    for (int mm = 0; mm < Mc; ++mm) {
                        pmI[mm] = ocrt_interp_water_mode_at_mu(&atm, I_m_node, n_mu, mm, mu_wx);
                        pmQ[mm] = ocrt_interp_water_mode_at_mu(&atm, Q_m_node, n_mu, mm, mu_wx);
                        pmU[mm] = ocrt_interp_water_mode_at_mu(&atm, U_m_node, n_mu, mm, mu_wx);
                    }
                    double Lw_I = rt_solver_reconstruct_phi(pmI, Mc - 1, dphi_rad) * f_scale;
                    double Lw_Q = rt_solver_reconstruct_phi(pmQ, Mc - 1, dphi_rad) * f_scale;
                    double Lw_U =  rt_solver_reconstruct_phi_sin(pmU, Mc - 1, dphi_rad) * f_scale;
                    double M_T_wa_node[9];
                    rt_air_water_T_wa(mu_wx, n_water, o.q_convention, M_T_wa_node);
                    mu_w_nodes[nv] = mu_wx;
                    mu_a_nodes[nv] = mu_ax;
                    lw_I_nodes[nv] = Lw_I;
                    lw_Q_nodes[nv] = Lw_Q;
                    lw_U_nodes[nv] = Lw_U;
                    node_src[nv] = iex + 1;
                    ya_I[nv] = M_T_wa_node[0*3+0]*Lw_I + M_T_wa_node[0*3+1]*Lw_Q + M_T_wa_node[0*3+2]*Lw_U;
                    ya_Q[nv] = M_T_wa_node[1*3+0]*Lw_I + M_T_wa_node[1*3+1]*Lw_Q + M_T_wa_node[1*3+2]*Lw_U;
                    ya_U[nv] = M_T_wa_node[2*3+0]*Lw_I + M_T_wa_node[2*3+1]*Lw_Q + M_T_wa_node[2*3+2]*Lw_U;
                    nv++;
                }
            }
            if (nv >= 2) {
                const int clamped = (twa_operator_mode == 2);
                const int epmode = ocrt_water_twa_critical_endpoint_mode();
                if (epmode == 3 && nv < n_mu + 8) {
                    mu_a_nodes[nv] = 0.0;
                    ya_I[nv] = 0.0;
                    ya_Q[nv] = 0.0;
                    ya_U[nv] = 0.0;
                    nv++;
                }
                if (epmode == 1) {
                    I_air = ocrt_linear_endpoint_eval_sorted(nv, mu_a_nodes, ya_I, mu_view_air);
                    Q_air = ocrt_linear_endpoint_eval_sorted(nv, mu_a_nodes, ya_Q, mu_view_air);
                    U_air = ocrt_linear_endpoint_eval_sorted(nv, mu_a_nodes, ya_U, mu_view_air);
                } else if (epmode == 2) {
                    I_air = ocrt_clamp_first_eval_sorted(nv, mu_a_nodes, ya_I, mu_view_air, clamped);
                    Q_air = ocrt_clamp_first_eval_sorted(nv, mu_a_nodes, ya_Q, mu_view_air, clamped);
                    U_air = ocrt_clamp_first_eval_sorted(nv, mu_a_nodes, ya_U, mu_view_air, clamped);
                } else {
                    I_air = ocrt_cubic_spline_eval_sorted(nv, mu_a_nodes, ya_I, mu_view_air, clamped);
                    Q_air = ocrt_cubic_spline_eval_sorted(nv, mu_a_nodes, ya_Q, mu_view_air, clamped);
                    U_air = ocrt_cubic_spline_eval_sorted(nv, mu_a_nodes, ya_U, mu_view_air, clamped);
                }
            } else {
                /* Fallback to pointwise if no valid transmitted nodes exist. */
                double M_T_wa[9];
                rt_air_water_T_wa(mu_view_water, n_water, o.q_convention, M_T_wa);
                I_air = M_T_wa[0*3+0]*I_view + M_T_wa[0*3+1]*Q_view + M_T_wa[0*3+2]*U_view;
                Q_air = M_T_wa[1*3+0]*I_view + M_T_wa[1*3+1]*Q_view + M_T_wa[1*3+2]*U_view;
                U_air = M_T_wa[2*3+0]*I_view + M_T_wa[2*3+1]*Q_view + M_T_wa[2*3+2]*U_view;
            }
            {
                const char *mdump = g_wrt_env.s_OCRT_DUMP_WATER_TWA_MODE_CONTRIB;
                if (mdump && mdump[0] && nv >= 2) {
                    FILE *mf = fopen(mdump, "a");
                    if (mf) {
                        static int mhdr = 0;
                        if (!mhdr) {
                            fprintf(mf,
                                "wavelength_nm,sza_deg,vza_deg,raa_deg,m,parity,operator_mode,n_mu,nv,target_mu_air,target_mu_w,n_water,"
                                "coef_IQ,coef_U,Iair_m,Qair_m,Uair_m,Iair_final,Qair_final,Uair_final\n");
                            mhdr = 1;
                        }
                        const double phi_arg_base = dphi_rad + RT_F_SOLAR_PI;
                        const int clamped_mdump = (twa_operator_mode == 2);
                        for (int mm = 0; mm < Mc; ++mm) {
                            double *ymI = (double*)calloc((size_t)nv, sizeof(double));
                            double *ymQ = (double*)calloc((size_t)nv, sizeof(double));
                            double *ymU = (double*)calloc((size_t)nv, sizeof(double));
                            if (!ymI || !ymQ || !ymU) { free(ymI); free(ymQ); free(ymU); break; }
                            double coef_iq = (mm == 0) ? 1.0 : 2.0 * cos((double)mm * phi_arg_base);
                            double coef_u  = (mm == 0) ? 0.0 : -2.0 * sin((double)mm * phi_arg_base);
                            for (int ii = 0; ii < nv; ++ii) {
                                double mu_wi = mu_w_nodes[ii];
                                double LwI_m = ocrt_interp_water_mode_at_mu(&atm, I_m_node, n_mu, mm, mu_wi) * f_scale * coef_iq;
                                double LwQ_m = ocrt_interp_water_mode_at_mu(&atm, Q_m_node, n_mu, mm, mu_wi) * f_scale * coef_iq;
                                double LwU_m = ocrt_interp_water_mode_at_mu(&atm, U_m_node, n_mu, mm, mu_wi) * f_scale * coef_u;
                                double M_T_wa_m[9];
                                rt_air_water_T_wa(mu_wi, n_water, o.q_convention, M_T_wa_m);
                                ymI[ii] = M_T_wa_m[0*3+0]*LwI_m + M_T_wa_m[0*3+1]*LwQ_m + M_T_wa_m[0*3+2]*LwU_m;
                                ymQ[ii] = M_T_wa_m[1*3+0]*LwI_m + M_T_wa_m[1*3+1]*LwQ_m + M_T_wa_m[1*3+2]*LwU_m;
                                ymU[ii] = M_T_wa_m[2*3+0]*LwI_m + M_T_wa_m[2*3+1]*LwQ_m + M_T_wa_m[2*3+2]*LwU_m;
                            }
                            double Im = ocrt_cubic_spline_eval_sorted(nv, mu_a_nodes, ymI, mu_view_air, clamped_mdump);
                            double Qm = ocrt_cubic_spline_eval_sorted(nv, mu_a_nodes, ymQ, mu_view_air, clamped_mdump);
                            double Um = ocrt_cubic_spline_eval_sorted(nv, mu_a_nodes, ymU, mu_view_air, clamped_mdump);
                            {
                                const char *nmdump = g_wrt_env.s_OCRT_DUMP_WATER_TWA_MODE_NODE_CONTRIB;
                                if (nmdump && nmdump[0] && (mm == 0 || mm == 1)) {
                                    FILE *nmf = fopen(nmdump, "a");
                                    if (nmf) {
                                        if (ftell(nmf) == 0L) {
                                            fprintf(nmf,
                                                "wavelength_nm,sza_deg,vza_deg,raa_deg,m,parity,operator_mode,n_mu,nv,target_mu_air,target_mu_w,n_water,coef_IQ,coef_U,node_index,node_src,mu_w,mu_air,LwI_m,LwQ_m,LwU_m,TwaI_m,TwaQ_m,TwaU_m,Iair_m,Qair_m,Uair_m,Iair_final,Qair_final,Uair_final\n");
                                        }
                                        for (int ii = 0; ii < nv; ++ii) {
                                            const double mu_wi2 = mu_w_nodes[ii];
                                            const double LwI_m2 = ocrt_interp_water_mode_at_mu(&atm, I_m_node, n_mu, mm, mu_wi2) * f_scale * coef_iq;
                                            const double LwQ_m2 = ocrt_interp_water_mode_at_mu(&atm, Q_m_node, n_mu, mm, mu_wi2) * f_scale * coef_iq;
                                            const double LwU_m2 = ocrt_interp_water_mode_at_mu(&atm, U_m_node, n_mu, mm, mu_wi2) * f_scale * coef_u;
                                            fprintf(nmf,
                                                "%.10g,%.10g,%.10g,%.10g,%d,%s,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%.17g,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                                                lambda_nm, sza_deg_air, vza_deg_air, raa_deg,
                                                mm, (mm % 2) ? "odd" : "even", twa_operator_mode, n_mu, nv,
                                                mu_view_air, mu_view_water, n_water, coef_iq, coef_u,
                                                ii, node_src[ii], mu_w_nodes[ii], mu_a_nodes[ii],
                                                LwI_m2, LwQ_m2, LwU_m2, ymI[ii], ymQ[ii], ymU[ii],
                                                Im, Qm, Um, I_air, Q_air, U_air);
                                        }
                                        fclose(nmf);
                                    }
                                }
                            }
                            fprintf(mf,
                                "%.10g,%.10g,%.10g,%.10g,%d,%s,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                                lambda_nm, sza_deg_air, vza_deg_air, raa_deg,
                                mm, (mm % 2) ? "odd" : "even", twa_operator_mode, n_mu, nv,
                                mu_view_air, mu_view_water, n_water, coef_iq, coef_u,
                                Im, Qm, Um, I_air, Q_air, U_air);
                            free(ymI); free(ymQ); free(ymU);
                        }
                        fclose(mf);
                    }
                }
            }

            {
                const char *ndump = g_wrt_env.s_OCRT_DUMP_TWA_NODESET;
                if (ndump && ndump[0]) {
                    FILE *nf = fopen(ndump, "a");
                    if (nf) {
                        if (ftell(nf) == 0L) {
                            fprintf(nf,
                                "code,wavelength_nm,sza_deg,vza_deg,raa_deg,operator_mode,nodeset_mode,endpoint_mode,target_mu_air,target_mu_w,n_water,node_index,node_src,mu_w,mu_air,Lw_I,Lw_Q,Lw_U,Twa_I,Twa_Q,Twa_U,final_Iair,final_Qair,final_Uair\n");
                        }
                        for (int ii = 0; ii < nv; ++ii) {
                            fprintf(nf,
                                "OCRT,%.10g,%.10g,%.10g,%.10g,%d,%d,%d,%.12e,%.12e,%.12e,%d,%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                                lambda_nm, sza_deg_air, vza_deg_air, raa_deg,
                                twa_operator_mode, ocrt_water_twa_nodeset_mode(),
                                ocrt_water_twa_critical_endpoint_mode(),
                                mu_view_air, mu_view_water, n_water, ii, node_src[ii],
                                mu_w_nodes[ii], mu_a_nodes[ii],
                                lw_I_nodes[ii], lw_Q_nodes[ii], lw_U_nodes[ii],
                                ya_I[ii], ya_Q[ii], ya_U[ii], I_air, Q_air, U_air);
                        }
                        fclose(nf);
                    }
                }
            }
            free(mu_w_nodes); free(lw_I_nodes); free(lw_Q_nodes); free(lw_U_nodes); free(node_src);
            free(mu_a_nodes); free(ya_I); free(ya_Q); free(ya_U);
            { int dtwap = g_wrt_env.f_OCRT_DUMP_TWA_PROVENANCE;
              if (dtwap) {
                  double M_T_wa[9]; rt_air_water_T_wa(mu_view_water, n_water, o.q_convention, M_T_wa);
                  double Ed0p_tmp = F_sun * mu_sun_air;
                  double rrs_tmp  = (Ed_total > 0.0) ? I_view / Ed_total : 0.0;
                  double Rrs_tmp  = (Ed0p_tmp > 0.0) ? I_air / Ed0p_tmp : 0.0;
                  fprintf(stderr,
                      "TWA_PROV wl=%.3f sza=%.6f vza=%.6f raa=%.6f branch=%s mu_air=%.12g mu_w=%.12g n=%.12g f_scale=%.12e "
                      "I0m=%.12e Q0m=%.12e U0m=%.12e Ed0m=%.12e rrs0m=%.12e "
                      "M00=%.12e M01=%.12e M02=%.12e M10=%.12e M11=%.12e M12=%.12e M20=%.12e M21=%.12e M22=%.12e "
                      "Iair=%.12e Qair=%.12e Uair=%.12e Ed0p=%.12e Rrs0p=%.12e\n",
                      lambda_nm, sza_deg_air, vza_deg_air, raa_deg,
                      (twa_operator_mode == 2) ? "flat_osoaa_clamped" : "flat_osoaa_natural",
                      mu_view_air, mu_view_water, n_water, f_scale,
                      I_view, Q_view, U_view, Ed_total, rrs_tmp,
                      M_T_wa[0], M_T_wa[1], M_T_wa[2], M_T_wa[3], M_T_wa[4], M_T_wa[5], M_T_wa[6], M_T_wa[7], M_T_wa[8],
                      I_air, Q_air, U_air, Ed0p_tmp, Rrs_tmp);
              } }
        } else {
            /* Flat Fresnel T_wa: single refracted view direction. */
            double M_T_wa[9];
            rt_air_water_T_wa(mu_view_water, n_water, o.q_convention, M_T_wa);
            if (g_wrt_env.s_OCRT_MB_TRACE) fprintf(stderr, "[IA:M2]\n");
            I_air = M_T_wa[0*3+0]*I_view + M_T_wa[0*3+1]*Q_view + M_T_wa[0*3+2]*U_view;
            Q_air = M_T_wa[1*3+0]*I_view + M_T_wa[1*3+1]*Q_view + M_T_wa[1*3+2]*U_view;
            U_air = M_T_wa[2*3+0]*I_view + M_T_wa[2*3+1]*Q_view + M_T_wa[2*3+2]*U_view;
            { int dtwap = g_wrt_env.f_OCRT_DUMP_TWA_PROVENANCE;
              if (dtwap) {
                  double Ed0p_tmp = F_sun * mu_sun_air;
                  double rrs_tmp  = (Ed_total > 0.0) ? I_view / Ed_total : 0.0;
                  double Rrs_tmp  = (Ed0p_tmp > 0.0) ? I_air / Ed0p_tmp : 0.0;
                  fprintf(stderr,
                      "TWA_PROV wl=%.3f sza=%.6f vza=%.6f raa=%.6f branch=flat mu_air=%.12g mu_w=%.12g n=%.12g f_scale=%.12e "
                      "I0m=%.12e Q0m=%.12e U0m=%.12e Ed0m=%.12e rrs0m=%.12e "
                      "M00=%.12e M01=%.12e M02=%.12e M10=%.12e M11=%.12e M12=%.12e M20=%.12e M21=%.12e M22=%.12e "
                      "Iair=%.12e Qair=%.12e Uair=%.12e Ed0p=%.12e Rrs0p=%.12e\n",
                      lambda_nm, sza_deg_air, vza_deg_air, raa_deg, mu_view_air, mu_view_water, n_water, f_scale,
                      I_view, Q_view, U_view, Ed_total, rrs_tmp,
                      M_T_wa[0], M_T_wa[1], M_T_wa[2], M_T_wa[3], M_T_wa[4], M_T_wa[5], M_T_wa[6], M_T_wa[7], M_T_wa[8],
                      I_air, Q_air, U_air, Ed0p_tmp, Rrs_tmp);
              } }
        }
    } else {
        /* Cox-Munk water->air transmission (option A), SLOPE-domain integration.
         * The forward direction integral ∫ f_t(w->v) L_w(w) (w·n) dω_w has a
         * near-specular peak far narrower than any practical mu-grid at low wind
         * (convergence needs ~1e5 mu-samples). Change variables to the facet
         * slope (z_x,z_y): the Gaussian slope pdf is smooth, so a modest grid
         * resolves it at every wind. For a fixed view v, each slope fixes the
         * facet normal m and — by refraction — the single in-water direction w
         * that emerges along v (no TIR branch: sin ω_w = sin ω_a / n <= 1/n).
         * Measure conversion (Walter 2007 ∂ω_m/∂ω_w = n² μ_wf / Nsq, and
         * dω_m = cos³β dz_x dz_y):
         *     dω_w = cos³β · Nsq/(n² μ_wf) dz_x dz_y =: J dz_x dz_y.
         * The validated Cox-Munk BTDF kernel (which already carries P(slope)/cos⁴β,
         * the Fresnel transmission and the Stokes rotations) is reused; L_w(w) is
         * reconstructed from the per-node Fourier modes (linear-in-μ interp,
         * rule 11). Modes are pi-normalized; SI f_scale applied afterwards. */
        const double phi_v = rt_raa_to_atm_fourier_phi(raa_deg);
        const double sin_v = sqrt(fmax(0.0, 1.0 - mu_view_air * mu_view_air));
        const double vx = sin_v * cos(phi_v), vy = sin_v * sin(phi_v), vz = mu_view_air;
        double sigma_sq = surface_slope_variance(o.wind_speed, o.cox_munk_sigma_type);
        double sigma_s  = sqrt(0.5 * sigma_sq);          /* per-component slope std */
        const int    NS  = 61;                            /* slope samples per axis */
        const double Lc  = 5.0 * sigma_s;                 /* half-extent (±5 σ_s) */
        const double dz  = (NS > 1) ? (2.0 * Lc / (double)(NS - 1)) : 0.0;
        const int    Mc  = m_loop_max + 1;
        double acc_I = 0.0, acc_Q = 0.0, acc_U = 0.0;
#ifdef OCRT_FAST_KERNELS
        /* v1.10 S10 (2026-07-10): coupling-first view-0+.
         * The slope-domain 61x61 integral below is OVERWRITTEN by the
         * Fourier coupling [IA:CPL4908] whenever that call succeeds (see
         * the PRODUCTION 2026-06-24 note further down; probe-verified
         * 2026-07-10: production rows always take CPL).  S6W probes put
         * this dead integral at ~4.0 ms/row on warm batch rows.  Skip it
         * when the previous row's coupling succeeded; if the coupling
         * ever fails afterwards, jump back and run the slope integral as
         * the original fallback (identical values on both paths).  The
         * success flag is thread-local state, NOT a physics switch; no
         * angle/condition branching is introduced (Jae rule respected).
         * FASTK-only; strict branch keeps the original order verbatim. */
        static _Thread_local int s10_prev_ok = 0;
        int s10_ran_slope = !s10_prev_ok;
        int s10_cpl_hit = 0;
s10_slope_entry:
        if (s10_ran_slope) {
#endif
        for (int iz = 0; iz < NS; ++iz) {
            double zx = -Lc + iz * dz;
            for (int jz = 0; jz < NS; ++jz) {
                double zy = -Lc + jz * dz;
                double r2 = 1.0 + zx*zx + zy*zy, r = sqrt(r2);
                double mx = -zx / r, my = -zy / r, mz = 1.0 / r;   /* facet normal */
                double vdotm = vx*mx + vy*my + vz*mz;              /* cos ω_a (>0) */
                if (vdotm <= 1e-9) continue;
                /* refract v (air) through facet into water -> in-water dir w */
                double tvx = vx - vdotm*mx, tvy = vy - vdotm*my, tvz = vz - vdotm*mz;
                double twx = tvx / n_water, twy = tvy / n_water, twz = tvz / n_water;
                double tw2 = twx*twx + twy*twy + twz*twz;
                if (tw2 >= 1.0) continue;
                double cw_facet = sqrt(1.0 - tw2);                 /* cos ω_w = μ_wf */
                double wx = twx + cw_facet*mx, wy = twy + cw_facet*my, wz = twz + cw_facet*mz;
                if (wz <= 1e-9) continue;                          /* w must go upward */
                double mu_w  = wz;
                double phi_w = atan2(wy, wx);
                double wdotv = wx*vx + wy*vy + wz*vz;              /* cosΨ */
                double Nsq   = n_water*n_water + 1.0 - 2.0*n_water*wdotv;
                double cosb  = 1.0 / r;
                double J = (cosb*cosb*cosb) * Nsq / (n_water*n_water * cw_facet);
                /* L_w(mu_w, phi_w): linear-in-μ mode interp + azimuth sum */
                int jl, jh;
                if (mu_w <= atm.rm[+1])         { jl = 1;        jh = 2; }
                else if (mu_w >= atm.rm[+n_mu]) { jl = n_mu - 1; jh = n_mu; }
                else { jl = 1; jh = 2;
                    for (int jp = 1; jp < n_mu; ++jp)
                        if (atm.rm[+jp] <= mu_w && mu_w < atm.rm[+jp+1]) { jl = jp; jh = jp+1; break; } }
                double ml = atm.rm[+jl], mh = atm.rm[+jh], dn = mh - ml;
                double wh = (dn != 0.0) ? (mu_w - ml)/dn : 0.0, wl = 1.0 - wh;
                double Lw_I = 0.0, Lw_Q = 0.0, Lw_U = 0.0;
                for (int mm = 0; mm < Mc; ++mm) {
                    size_t a = (size_t)mm*(size_t)n_mu + (size_t)(jl-1);
                    size_t b = (size_t)mm*(size_t)n_mu + (size_t)(jh-1);
                    double Imm = wl*I_m_node[a] + wh*I_m_node[b];
                    double Qmm = wl*Q_m_node[a] + wh*Q_m_node[b];
                    double Umm = wl*U_m_node[a] + wh*U_m_node[b];
                    double f = (mm==0) ? 1.0 : 2.0;
                    /* +pi azimuth shift to match rt_solver_reconstruct_phi (OS.f L627
                     * convention: Dphi=0 = forward scatter). The per-m modes carry
                     * this convention, so L_w(mu_w,phi_w) must reconstruct with the
                     * same shifted argument; without it the m>=1 modes pick up the
                     * wrong sign and L_w(view) != I_diffuse_view. (B2 fix 2026-06-03.) */
                    double phi_arg = phi_w + M_PI;
                    double cmp = cos((double)mm*phi_arg), smp = sin((double)mm*phi_arg);
                    Lw_I +=  f*Imm*cmp;
                    Lw_Q +=  f*Qmm*cmp;
                    Lw_U += -f*Umm*smp;
                }
                double dphi = phi_v - phi_w;
                double M_T[9];   /* BTDF kernel already includes P(slope)/cos⁴β */
                surface_T_coxmunk_trig(mu_view_air, mu_w, cos(dphi), sin(dphi),
                                       o.wind_speed, o.cox_munk_sigma_type,
                                       n_water, o.q_convention, M_T);
                double wgt = mu_w * J * dz * dz;
                acc_I += wgt * (M_T[0]*Lw_I + M_T[1]*Lw_Q + M_T[2]*Lw_U);
                acc_Q += wgt * (M_T[3]*Lw_I + M_T[4]*Lw_Q + M_T[5]*Lw_U);
                acc_U += wgt * (M_T[6]*Lw_I + M_T[7]*Lw_Q + M_T[8]*Lw_U);
            }
        }
        if (g_wrt_env.s_OCRT_MB_TRACE) fprintf(stderr, "[IA:SLOPE]\n");
#ifdef OCRT_FAST_KERNELS
        }  /* s10_ran_slope */
#endif
        I_air = acc_I * f_scale;
        Q_air = acc_Q * f_scale;
        U_air = acc_U * f_scale;
        /* PRODUCTION (2026-06-24): the slope-domain view-0+ integral above is
         * sza-biased at near-nadir sun (validated vs OSOAA: -14% at sza=0,
         * non-monotonic for sza<3). Replace its view 0+ result with the
         * Fourier-mode T_wa coupling rt_air_water_couple_water_to_atm, which
         * matches OSOAA to +0.3~1.1% across all view angles and is sza-invariant.
         * Same Cox-Munk BTDF kernel; only the integration differs (real-space
         * slope-grid -> azimuth Fourier).  The final reported field uses the
         * canonical reconstruction: phi_v+pi, I/Q cosine and U positive sine. */
        { int nm=m_loop_max+1; double mu_atm_t[1]={mu_view_air};
          double *muw=(double*)calloc((size_t)n_mu,sizeof(double));
          double *Ia=(double*)calloc((size_t)nm,sizeof(double));
          double *Qa=(double*)calloc((size_t)nm,sizeof(double));
          double *Ua=(double*)calloc((size_t)nm,sizeof(double));
          double *ww_mb=(double*)calloc((size_t)n_mu,sizeof(double));
          if(muw&&ww_mb&&Ia&&Qa&&Ua){ for(int k=0;k<n_mu;k++){ muw[k]=atm.rm[+(k+1)]; ww_mb[k]=atm.gb[+(k+1)]; }
            /* PRODUCTION VIEW-0+ CHAIN (verified by probes 2026-07-10):
             *   [IA:SLOPE] slope integral (overwritten) -> THIS coupling call
             *   [IA:CPL4908] -> stored at I_0plus_view [VIEW0P].
             * w_water_pos=NULL makes the coupling take its FLAT path (the
             * rough branch requires water weights); June production chose
             * flat+spline.  Beyond-critical directions are handled INSIDE the
             * flat path by spline EXTRAPOLATION below the lowest table node -
             * see ocrt_water_twa_critical_endpoint_mode() (Step 55 comment):
             * default splint_extrapolate.  Any additive beyond-critical term
             * must NOT be combined with that extrapolation (double-counts;
             * measured 2026-07-10).  Diagnostic-only toggle (default off):
             * OCRT_MB_ROUGH_VIEW=1 passes real weights -> rough branch, for
             * OSOAA A/B comparison; NOT a production switch. */
            int rcf=rt_air_water_couple_water_to_atm(I_m_node,Q_m_node,U_m_node,
                muw,(g_wrt_env.s_OCRT_MB_ROUGH_VIEW?ww_mb:NULL),n_mu,m_loop_max,n_water,o.q_convention,o.wind_speed,
                o.cox_munk_sigma_type,mu_atm_t,1,NULL,Ia,Qa,Ua); /* w_atm=NULL -> density form */
            if(rcf==0){
              double LI=0.0,LQ=0.0,LU=0.0;
              for(int mm=0;mm<nm;mm++){ double f=(mm==0)?1.0:2.0;
                double pa=phi_v+M_PI; double c=cos((double)mm*pa), s=sin((double)mm*pa);
                LI+=f*Ia[mm]*c; LQ+=f*Qa[mm]*c; LU+= f*Ua[mm]*s; }
              if (g_wrt_env.s_OCRT_MB_TRACE) fprintf(stderr, "[IA:CPL4908]\n");
              I_air=LI*f_scale; Q_air=LQ*f_scale; U_air=LU*f_scale;
#ifdef OCRT_FAST_KERNELS
              s10_cpl_hit = 1; s10_prev_ok = 1;
#endif
              /* v1.10 OPTION-2 RE-ESCAPE (final-radiance units, m=0). */
              /* SPRINT GUARD (2026-07-10): re-escape view term OPT-IN ONLY
               * (OCRT_MB_CLOSURE=1); default off = pre-surgery bit state.
               * wind==0 never enters regardless (flat branch preserved). */
              static _Thread_local int mb_on = -1;
              if (mb_on < 0) mb_on = (getenv("OCRT_MB_CLOSURE") != NULL);
              if (mb_on && o.wind_speed > 0.0 && mu_view_air > 1e-6) {
                  enum { MBV_NB = 64 };
                  static _Thread_local double v_mu[128];
                  static _Thread_local double v_h[128][MBV_NB];
                  static _Thread_local double v_hd[128][MBV_NB];
                  static _Thread_local int    v_n = 0;
                  static _Thread_local double v_nw = 0.0, v_ws = 0.0;
                  if (v_nw != n_water || v_ws != o.wind_speed) {
                      v_n = 0; v_nw = n_water; v_ws = o.wind_speed;
                  }
                  const double dmu = 1.0 / (double)MBV_NB;
                  int bv = (int)(mu_view_air * MBV_NB);
                  if (bv >= MBV_NB) bv = MBV_NB - 1;
                  double addI = 0.0;
                  for (int k = 0; k < n_mu; ++k) {
                      double mu_k = muw[k];
                      double w_k  = ww_mb ? ww_mb[k] : 0.0;
                      if (!(mu_k > 0.0) || !(w_k > 0.0)) continue;
                      int ci = -1;
                      for (int q = 0; q < v_n; ++q)
                          if (v_mu[q] == mu_k) { ci = q; break; }
                      if (ci < 0 && v_n < 128) {
                          double Th_c, Tu_c, Wr_c;
                          if (surface_wa_multibounce_cascade(n_water,
                                  o.wind_speed, mu_k, 300000L,
                                  &Th_c, &Tu_c, &Wr_c,
                                  v_h[v_n], v_hd[v_n], NULL, MBV_NB) == 0) {
                              v_mu[v_n] = mu_k; ci = v_n; v_n++;
                          }
                      }
                      if (ci < 0) continue;
                      /* re-escape (all nodes) + domain-completion direct
                       * escape for BEYOND-critical nodes only (they have no
                       * flat refraction angle; mu_crit is a physical constant
                       * of n_water, not an ad-hoc switch). */
                      /* NOTE(2026-07-10): direct-escape domain completion is
                       * DISABLED - the flat operator already covers beyond-
                       * critical directions via spline EXTRAPOLATION below the
                       * lowest table node, so adding h_direct double-counts
                       * (measured: add jumped to 1.21e-2, 4x target).  Proper
                       * pairing (zero-endpoint mode + full h) is a next-session
                       * design decision.  Re-escape term only, for now. */
                      double hh = v_h[ci][bv];
                      if (hh <= 0.0) continue;
                      double Lw0 = I_m_node[(size_t)k] * f_scale;
                      addI += mu_k * w_k * Lw0 * hh / (mu_view_air * dmu);
                  }
                  if (g_wrt_env.s_OCRT_MB_TRACE)
                      fprintf(stderr, "[MBV] add=%.4e\n", addI);
                  I_air += addI;
              } } }
          free(muw);free(ww_mb);free(Ia);free(Qa);free(Ua); }
#ifdef OCRT_FAST_KERNELS
        if (!s10_cpl_hit) { s10_prev_ok = 0;
            if (!s10_ran_slope) { s10_ran_slope = 1; goto s10_slope_entry; } }
#endif
        { int dtwa = g_wrt_env.f_OCRT_DUMP_TWA;
          if (dtwa) { double Mfw[9]; rt_air_water_T_wa(mu_view_water, n_water, o.q_convention, Mfw);
            double flatI = Mfw[0]*I_view + Mfw[1]*Q_view + Mfw[2]*U_view;
            fprintf(stderr, "TWA wind=%.3f rough_Iair=%.6e flat_Iair=%.6e ratio=%.4f | I_view=%.6e mu_vw=%.5f\n",
                    o.wind_speed, I_air, flatI, (flatI!=0.0)?I_air/flatI:0.0, I_view, mu_view_water); } }
    }

    /* Step92 diagnostic-only: final assembled m=0/m=1 queue-equivalent correction.
     * Applied after the flat/rough TWA path and before Rrs output. */
    ocrt_step92_apply_final_mode_queue(lambda_nm, vza_deg_air, raa_deg, &I_air, &Q_air, &U_air);

    if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] m14 abs=%.4f\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
    /* 14. Above-water downwelling irradiance.
     *
     *   B.3 simple no-atmosphere case: Ed_0plus_air = F_sun × μ_sun_air
     *   (just-above-surface horizontal solar flux, no surface loss applied here).
     *
     *   B.4 (future): atmospheric BOA irradiance from atm SOS output, including
     *   atm diffuse Ed component.
     */
    double Ed_0plus_air = F_sun * mu_sun_air;

    /* Step115A: Step114 grid-table correction is intentionally NOT applied
     * inside the water-only solver.  rt_solve_case_ocean() later resets the
     * air-side Ed0+ denominator to F_sun_BOA*mu0; applying an Rrs-unit delta
     * here uses the pre-reset local denominator and attenuates the requested
     * correction by Ed_internal/Ed_final.  The exact-grid correction is now
     * applied in rt_solver.c after Ed_above_air is finalized. */

    /* 15. r_rs (underwater) and R_rs (above-water) */
    double r_rs_0minus   = (Ed_total > 0.0) ? I_view / Ed_total : 0.0;
    double r_rs_0minus_Q = (Ed_total != 0.0) ? Q_view / Ed_total : 0.0;
    double r_rs_0minus_U = (Ed_total != 0.0) ? U_view / Ed_total : 0.0;
    double R_rs_0plus    = (Ed_0plus_air > 0.0) ? I_air / Ed_0plus_air : 0.0;
    double R_rs_0plus_Q  = (Ed_0plus_air != 0.0) ? Q_air / Ed_0plus_air : 0.0;
    double R_rs_0plus_U  = (Ed_0plus_air != 0.0) ? U_air / Ed_0plus_air : 0.0;

    /* Step-16 diagnostic: provenance of the above-water Rrs path.
     * This dump ties together the exact 0- total Stokes used for rrs0minus,
     * the Stokes vector handed to the flat/rough TWA operator, the TWA output,
     * and the final printed Rrs0plus denominator.  It is diagnostic-only and
     * has no effect unless OCRT_DUMP_TWA_PROVENANCE points to a CSV path. */
    {
        const char *prov_path = g_wrt_env.s_OCRT_DUMP_TWA_PROVENANCE;
        if (prov_path && prov_path[0]) {
            FILE *pf = fopen(prov_path, "a");
            if (pf) {
                long pos = ftell(pf);
                if (pos == 0L) {
                    fprintf(pf,
                        "route,wavelength_nm,wind_speed,sza_deg_air,vza_deg_air,raa_deg,"
                        "mu_sun_air,mu_sun_water,mu_view_air,mu_view_water,"
                        "I_0minus,Q_0minus,U_0minus,Ed_0minus,rrs0minus,"
                        "I_twa_input,Q_twa_input,U_twa_input,"
                        "I_0plus,Q_0plus,U_0plus,Ed_0plus,Rrs0plus,"
                        "twa_I_over_0minus,Ed0plus_over_Ed0minus,orders,conv,tau_max,z_max,a_total,b_total,omega\n");
                }
                const char *route = g_wrt_env.s_OCRT_ROUTE_NAME;
                if (!route || !route[0]) route = "unknown";
                fprintf(pf,
                    "%s,%.10g,%.10g,%.10g,%.10g,%.10g,"
                    "%.12e,%.12e,%.12e,%.12e,"
                    "%.12e,%.12e,%.12e,%.12e,%.12e,"
                    "%.12e,%.12e,%.12e,"
                    "%.12e,%.12e,%.12e,%.12e,%.12e,"
                    "%.12e,%.12e,%d,%d,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                    route, lambda_nm, o.wind_speed, sza_deg_air, vza_deg_air, raa_deg,
                    mu_sun_air, mu_sun_water, mu_view_air, mu_view_water,
                    I_view, Q_view, U_view, Ed_total, r_rs_0minus,
                    I_view, Q_view, U_view,
                    I_air, Q_air, U_air, Ed_0plus_air, R_rs_0plus,
                    (fabs(I_view) > 0.0 ? I_air / I_view : 0.0),
                    (fabs(Ed_total) > 0.0 ? Ed_0plus_air / Ed_total : 0.0),
                    max_orders_seen, all_conv, tau_max, z_max, a_tot, b_tot, omega);
                fclose(pf);
            }
        }
    }

    if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] m16 abs=%.4f\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
    /* 16. Fill result */
    result->max_orders_used  = max_orders_seen;
    result->all_converged    = all_conv;
    result->max_residual     = worst_resid;
    result->lambda_nm_used   = lambda_nm;
    result->T_water_C_used   = T_water_C;
    result->a_w_used         = a_w;
    result->b_w_used         = b_w;
    result->bb_w_used        = bb_w;
    result->a_cdom_used      = a_cdom;
    result->a_pig_used       = a_pig;
    result->a_chl_used       = a_phyto;
    result->b_pig_used       = b_pig;
    result->bb_pig_used      = bb_pig;
    result->a_min_used       = a_min;
    result->b_min_used       = b_min;
    result->bb_min_used      = bb_min;
    result->a_total_used     = a_tot;
    result->b_total_used     = b_tot;
    result->bb_total_used    = bb_tot;
    result->omega_used       = omega;
    result->tau_max_used     = tau_max;
    result->z_max_used       = z_max;

    /* B.5: expose the z=0⁻ upwelling Fourier field (I_m_node, already built for
     * the Cox-Munk T_wa BTDF) to the caller for water→air up-coupling. No-op
     * unless the caller pre-allocated all four buffers (backward compatible). */
    if (result->I_up_per_m && result->Q_up_per_m &&
        result->U_up_per_m && result->mu_water_pos) {
        for (int m2 = 0; m2 < M; ++m2) {
            for (int k = 1; k <= n_mu; ++k) {
                size_t idx = (size_t)m2 * (size_t)n_mu + (size_t)(k - 1);
                /* Store the PHYSICAL upwelling field (= internal SOS field ×
                 * f_scale, f_scale=F_sun/π — same SI scaling the solver applies
                 * to the view in I_view). The caller then applies only T_wa, no
                 * normalization. (SOS single-scatter is included here; at the
                 * view the solver swaps it for the Gibbs-free analytic Iss — that
                 * difference is ~0.01% at lmax≥100.) */
                result->I_up_per_m[idx] = I_m_node[idx] * f_scale;
                result->Q_up_per_m[idx] = Q_m_node[idx] * f_scale;
                result->U_up_per_m[idx] = U_m_node[idx] * f_scale;
            }
        }
        for (int k = 1; k <= n_mu; ++k) {
            result->mu_water_pos[k - 1] = atm.rm[+k];
            if (result->w_water_pos) result->w_water_pos[k - 1] = atm.gb[+k];   /* D3-0c: quadrature weights for operator gates */
            if (result->w_water_pos) result->w_water_pos[k - 1] = atm.gb[+k];
        }
        result->n_mu_water_filled = n_mu;
        result->m_max_filled      = M - 1;   /* modes 0..M-1 (= o.m_max_water) */

        /* Native coupled-LUT export: copy the already-computed view Fourier
         * samples and apply the same linear water->air coupling once per mode.
         * No SOS source, order loop, or phase build is repeated here. */
        result->view_m_max_filled = -1;
        result->view_modes_exact = 0;
        if (result->I_view_per_m && result->Q_view_per_m && result->U_view_per_m &&
            result->I_air_view_per_m && result->Q_air_view_per_m &&
            result->U_air_view_per_m && result->view_m_max_capacity >= M - 1) {
            for (int m2 = 0; m2 < M; ++m2) {
                result->I_view_per_m[m2] = I_m_view[m2] * f_scale;
                result->Q_view_per_m[m2] = Q_m_view[m2] * f_scale;
                result->U_view_per_m[m2] = U_m_view[m2] * f_scale;
            }
            const double mu_air_one[1] = { mu_view_air };
            double *tmpI = (double*)calloc((size_t)M, sizeof(double));
            double *tmpQ = (double*)calloc((size_t)M, sizeof(double));
            double *tmpU = (double*)calloc((size_t)M, sizeof(double));
            if (tmpI && tmpQ && tmpU &&
                rt_air_water_couple_water_to_atm(
                    I_m_node, Q_m_node, U_m_node,
                    atm.rm + 1, atm.gb + 1, n_mu, M - 1,
                    n_water, o.q_convention, o.wind_speed, o.cox_munk_sigma_type,
                    mu_air_one, 1, NULL, tmpI, tmpQ, tmpU) == 0) {
                for (int m2 = 0; m2 < M; ++m2) {
                    result->I_air_view_per_m[m2] = tmpI[m2] * f_scale;
                    result->Q_air_view_per_m[m2] = tmpQ[m2] * f_scale;
                    result->U_air_view_per_m[m2] = tmpU[m2] * f_scale;
                }
                result->view_m_max_filled = M - 1;
                /* The exported samples are the solved Fourier field before any
                 * target-angle-only non-Fourier correction.  Mark them exact
                 * only when the production result is a pure Fourier
                 * reconstruction.  Native angular LUT orchestration falls back
                 * to the authoritative scalar target path otherwise. */
                int target_nonfourier = 0;
                if (ccrr_particle_ff_mode && b_particle_input > 0.0 && ccrr_ff_betal_L >= 0)
                    target_nonfourier = 1;                    /* CCRR NT-TMS */
                if (o.fixed_bulk_iop_mode && o.water_mie_ss_mode == 1 &&
                    o.water_mie_truncation_mode == 1 && have_aip)
                    target_nonfourier = 1;                    /* water IMS */
                if (o.fixed_bulk_iop_mode && ocrt_debug_env("OCRT_FIXEDBULK_TMS"))
                    target_nonfourier = 1;                    /* diagnostic fixed-bulk TMS */
                {
                    const int sspol_off = (ocrt_debug_env("OCRT_SSPOL_CORR_OFF") != NULL);
                    const int sspol_allowed = (o.water_mie_moment_mode != 1) ||
                                              (ocrt_debug_env("OCRT_SSPOL_CORR_ON") != NULL);
                    if (!sspol_off && sspol_allowed && !have_particle_value_phase &&
                        have_aip && !ccrr_particle_ff_mode &&
                        b_particle > 0.0 && b_w < 1.0e-12)
                        target_nonfourier = 1;                /* exact-phase Q/U */
                }
                if (g_wrt_env.s_OCRT_WATER_NADIR_RAA_MODE)
                    target_nonfourier = 1;                    /* target RAA canonicalization */
                result->view_modes_exact =
                    (!target_nonfourier && ocrt_step92_queue_mode() != 1);
            }
            free(tmpI); free(tmpQ); free(tmpU);
        }

        /* B.5 self-check (OCRT_DUMP_UP_RECON): interpolate the EXPOSED nodal
         * field at μ_view_water per mode and compare to the solver's own
         * per-mode view value I_m_view[m]. Same data + same linear interp ⇒
         * max diff must be ~0; a non-zero diff flags a copy/grid/index bug. */
        { int dup = g_wrt_env.f_OCRT_DUMP_UP_RECON;
          if (dup) {
              double mu_t = mu_view_water; int jl = 1, jh = 2;
              if (mu_t <= atm.rm[+1])         { jl = 1;        jh = 2; }
              else if (mu_t >= atm.rm[+n_mu]) { jl = n_mu - 1; jh = n_mu; }
              else { for (int jp = 1; jp < n_mu; ++jp)
                         if (atm.rm[+jp] <= mu_t && mu_t < atm.rm[+jp+1]) { jl = jp; jh = jp+1; break; } }
              double ml = atm.rm[+jl], mh = atm.rm[+jh], dn = mh - ml;
              double wh = (dn != 0.0) ? (mu_t - ml)/dn : 0.0, wl = 1.0 - wh;
              double maxd = 0.0;
              for (int m2 = 0; m2 < M; ++m2) {
                  double a = wl*result->I_up_per_m[(size_t)m2*n_mu+(jl-1)]
                           + wh*result->I_up_per_m[(size_t)m2*n_mu+(jh-1)];
                  double d = fabs(a - I_m_view[m2]*f_scale); if (d > maxd) maxd = d;
              }
              fprintf(stderr, "[UP-RECON] mu_view_water=%.4f n_mu=%d M=%d | max|interp(exposed)_m - I_m_view_m*f_scale|=%.3e (should be ~0)\n",
                      mu_view_water, n_mu, M, maxd);
          }
        }
    }

    result->mu_sun_air       = mu_sun_air;
    result->mu_sun_water     = mu_sun_water;
    result->I_0minus_view    = I_view;
    result->Q_0minus_view    = Q_view;
    result->U_0minus_view    = U_view;
    if (g_wrt_env.s_OCRT_MB_TRACE)
        fprintf(stderr, "[VIEW0P] I_air=%.6e vza=%.3f\n", I_air, vza_deg_air);
    result->I_0plus_view     = I_air;
    result->Q_0plus_view     = Q_air;
    result->U_0plus_view     = U_air;
    result->Ed_0minus_water  = Ed_total;
    result->Eu_0minus_water  = Eu_total;
    result->Ed_0plus_air     = Ed_0plus_air;
    result->r_rs_0minus      = r_rs_0minus;
    result->r_rs_0minus_Q    = r_rs_0minus_Q;
    result->r_rs_0minus_U    = r_rs_0minus_U;
    result->R_rs_0plus       = R_rs_0plus;
    result->R_rs_0plus_Q     = R_rs_0plus_Q;
    result->R_rs_0plus_U     = R_rs_0plus_U;
    result->Kd_0minus        = Kd;
    result->Ku_0minus        = Ku;
    result->Ed_level1_water  = Ed_total_lvl1;
    result->tau_level1_used  = tau_lvl1;
    result->z_level1_used    = z_lvl1;

cleanup:
    free(mie_value_beamq_allm);
    if (mie_value_kernel_owned) free(mie_value_kernel_allm);
    if (have_fixed_phase_spline) rt_value_phase_interp_free(&fixed_phase_spline);
    if (have_ccrr_phase_spline) rt_value_phase_interp_free(&ccrr_phase_spline);
    free(fb_allm_kernel);
    fixed_bulk_phase_table_free(&fixed_phase_table);
    fixed_bulk_phase_table_free(&ccrr_particle_phase_table);
    /* v1.09 #24: aip is owned by the threadprivate cache g_wmc (both on hit and
     * on the store path), so it is NOT freed here — wmc_free() owns its lifetime
     * and releases it on the next key change or at process teardown. */
    /* (was: if (have_aip) aerosol_phase_interp_free(&aip);) */
    free(wbuf);   /* #21: single arena (src/tot/prim/view/pos/node/lvl1) */
    free(Rww_M);
    free(Rww_K); free(mu_pos_ww);
    if (water_nt_tms) rt_legendre_workspace_free(&ws_primary);
    rt_legendre_workspace_free(&ws);
    rt_atm_free(&atm);
    return rc;
}

#ifdef OCRT_TEST_DTP_CLOSURE
void ocrt_test_add_diffuse_top_primary(const rt_atm_t *atm,
                                       const rt_legendre_workspace_t *ws,
                                       int m,
                                       const double *ext_I,
                                       const double *ext_Q,
                                       const double *ext_U,
                                       const double *ext_mu, int ext_n,
                                       double *src_i, double *src_q,
                                       double *src_u) {
    ocrt_add_diffuse_top_primary(atm, ws, m, ext_I, ext_Q, ext_U,
                                 ext_mu, ext_n, src_i, src_q, src_u);
}
#endif
