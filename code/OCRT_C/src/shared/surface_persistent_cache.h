/* ============================================================================
 * shared/surface_persistent_cache.h
 *
 * DOC-REF: docs/PERSISTENT_COXMUNK_INTERFACE_CACHE_KO_2026-08-25.md
 * Versioned, content-addressed persistent cache for immutable Cox-Munk and
 * air-water interface Fourier operators.
 *
 * Concurrency contract:
 *   - cache files are immutable after publication;
 *   - production batch rows only read them;
 *   - no lock file, polling loop, barrier, or row-to-row build dependency;
 *   - cache BUILD mode is rejected by the in-process --batch-full-grid driver.
 *
 * A cache miss in REQUIRED mode is a fail-loud error, never a wait for another
 * row to finish.  BUILD is an explicit, single-prebuilder pre-run step
 * performed before the independent simulation rows are launched.  The cache
 * implementation uses ISO C file I/O only and contains no OS-specific path,
 * lock, process-ID, or filesystem API.
 * ========================================================================== */
#ifndef SHARED_SURFACE_PERSISTENT_CACHE_H
#define SHARED_SURFACE_PERSISTENT_CACHE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define OCRT_SURFACE_PCACHE_FORMAT_VERSION   1u
#define OCRT_SURFACE_PCACHE_ABI_VERSION      2026082501u
#define OCRT_SURFACE_PCACHE_PHYSICS_VERSION  2026082501u
#define OCRT_SURFACE_PCACHE_QFOLD_VERSION    1u

typedef enum {
    OCRT_SURFACE_PCACHE_OFF = 0,
    OCRT_SURFACE_PCACHE_READONLY = 1,
    OCRT_SURFACE_PCACHE_REQUIRED = 2,
    OCRT_SURFACE_PCACHE_BUILD = 3
} ocrt_surface_pcache_mode_t;

typedef enum {
    OCRT_SURFACE_PCACHE_OP_R_AA_RAW_ALLM = 1,
    OCRT_SURFACE_PCACHE_OP_R_WW_RAW_ALLM = 2,
    OCRT_SURFACE_PCACHE_OP_T_AW_FINAL_M  = 3,
    OCRT_SURFACE_PCACHE_OP_T_WA_RAW_ALLM = 4
} ocrt_surface_pcache_operator_t;

typedef struct {
    uint32_t operator_kind;
    int32_t mode_first;
    int32_t mode_count;
    int32_t n_o;
    int32_t n_i;
    int32_t n_phi_quad;
    int32_t sigma_type;
    int32_t q_convention;
    double wind_speed;
    double n_water;
    const double *mu_o;
    const double *mu_i;
} ocrt_surface_pcache_key_t;

/* Read OCRT_SURFACE_PERSIST_CACHE_MODE / _DIR / _TRACE.
 * Idempotent. OCRT main calls this before its OpenMP row region. */
int ocrt_surface_pcache_init_from_env(void);

ocrt_surface_pcache_mode_t ocrt_surface_pcache_mode(void);
const char *ocrt_surface_pcache_mode_name(void);

/* BUILD is a pre-run operation and is deliberately forbidden inside
 * --batch-full-grid. Returns 0 when batch execution is safe, -1 otherwise. */
int ocrt_surface_pcache_batch_guard(void);

/* Load an exact operator payload.
 * Returns: 1 exact hit; 0 disabled/optional miss; negative fail-loud error.
 * REQUIRED mode returns a negative value for a missing, stale, or corrupt file. */
int ocrt_surface_pcache_load(const ocrt_surface_pcache_key_t *key,
                             double *payload, size_t payload_count);

/* Publish an immutable payload in BUILD mode. No-op in other modes.
 * BUILD is a single-prebuilder operation.  The writer uses ISO C file I/O
 * and rename, never waits for another process, and never creates a lock file. */
int ocrt_surface_pcache_store(const ocrt_surface_pcache_key_t *key,
                              const double *payload, size_t payload_count);

#ifdef __cplusplus
}
#endif

#endif /* SHARED_SURFACE_PERSISTENT_CACHE_H */
