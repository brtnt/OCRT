#ifndef OCRT_RT_KERNEL_H
#define OCRT_RT_KERNEL_H

#include <stddef.h>

#include "rt_types.h"

/* Workspace for the angular Fourier representation of a three-Stokes phase
 * matrix.  Direction-indexed rows are offset so that [-n_mu,+n_mu] is valid.
 */
typedef struct {
    int n_mu;
    int l_max;

    double  *plm_storage;
    double **plm;
    double  *rrl_storage;
    double **rrl;
    double  *rtl_storage;
    double **rtl;

    double  *pfm_storage;
    double **phase_fourier_m;
    double  *gr_storage;
    double **gr_pol;
    double  *gt_storage;
    double **gt_pol;
    double  *arr_storage;
    double **arr_pol;
    double  *art_storage;
    double **art_pol;
    double  *att_storage;
    double **att_pol;

    /* Exact particle response to the Fresnel-transmitted direct-beam Q
     * component for the current Fourier mode.  The six standard direct
     * kernels contain the incoming-I and diffuse-vector couplings, but not
     * the complete incoming-Q solar column (I<-Q and U<-Q). */
    double *beam_q_to_i_direct;
    double *beam_q_to_q_direct;
    double *beam_q_to_u_direct;
    int     beam_q_direct_valid;

    /* Reusable storage for the SOS source operator.  These buffers are built
     * outside the scattering-order loop and then treated as immutable. */
    double *sos_source_arena;
    size_t  sos_source_arena_capacity;
    double *sos_kernel_pack;
    size_t  sos_kernel_pack_capacity;
    int     sos_kernel_pack_m;
    int     sos_kernel_pack_valid;
} rt_legendre_workspace_t;

int  rt_legendre_workspace_alloc(rt_legendre_workspace_t *ws,
                                 int n_mu,
                                 int l_max);
void rt_legendre_workspace_free(rt_legendre_workspace_t *ws);

/* Build the scalar and spin-2 angular bases for Fourier order m. */
int rt_legendre_compute(rt_legendre_workspace_t *ws,
                        rt_atm_t *atm,
                        int m);
int rt_legendre_compute_pol(rt_legendre_workspace_t *ws,
                            rt_atm_t *atm,
                            int m);

/* Contract the current angular basis with the phase-matrix moment arrays. */
int rt_kernel_phase_fourier(rt_legendre_workspace_t *ws,
                            int m,
                            const double *betal);
int rt_kernel_phase_fourier_pol(rt_legendre_workspace_t *ws,
                                int m,
                                const double *gammal);
int rt_kernel_phase_fourier_aerosol_full(rt_legendre_workspace_t *ws,
                                         int m,
                                         const double *alphal,
                                         const double *zetal);

#ifdef OCRT_FAST_KERNELS
/* Thread-local immutable cache for the six Fourier kernel tables. */
int  rt_mkc_begin(const rt_atm_t *atm,
                  const rt_legendre_workspace_t *ws,
                  int m_count_hint);
int  rt_mkc_load(int m, rt_legendre_workspace_t *ws);
void rt_mkc_store_pfm(int m, const rt_legendre_workspace_t *ws);
void rt_mkc_store_pol(int m, const rt_legendre_workspace_t *ws);
#endif

#endif /* OCRT_RT_KERNEL_H */
