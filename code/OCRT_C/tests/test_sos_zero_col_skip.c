#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rt_atm.h"
#include "rt_kernel.h"
#include "rt_sos_operator.h"

static uint64_t bits(double x) {
    uint64_t u;
    memcpy(&u, &x, sizeof u);
    return u;
}

static int same_bits(const double *a, const double *b, size_t n, const char *name) {
    for (size_t i = 0; i < n; ++i) {
        if (bits(a[i]) != bits(b[i])) {
            fprintf(stderr, "%s mismatch at %zu: %.17g %.17g\n", name, i, a[i], b[i]);
            return 0;
        }
    }
    return 1;
}

int main(void) {
    enum { N_MU=8, N_LAYERS=9, LMAX=12, M=3 };
    rt_atm_t atm={0};
    rt_legendre_workspace_t ws={0};
    double beta[LMAX+1], gamma[LMAX+1], alpha[LMAX+1], zeta[LMAX+1];
    for (int l=0; l<=LMAX; ++l) {
        beta[l]=(l==0)?1.0:(2.0*l+1.0)*pow(0.71,l);
        gamma[l]=(l<2)?0.0:0.07*pow(-0.53,l-2);
        alpha[l]=(l<2)?0.0:0.11*pow(0.49,l-2);
        zeta[l]=(l<2)?0.0:0.05*pow(-0.43,l-2);
    }
    if (rt_atm_alloc(&atm,N_LAYERS,N_MU)!=0) return 2;
    if (rt_atm_build_aerosol_rayleigh(&atm,0.16,0.0279,0.31,0.94,LMAX,
            beta,gamma,alpha,zeta,cos(35.0*3.14159265358979323846/180.0),
            RT_RAYLEIGH_MODEL_BODHAINE_1999,2.0)!=0) return 3;
    if (rt_legendre_workspace_alloc(&ws,N_MU,LMAX)!=0) return 4;
    if (rt_legendre_compute_pol(&ws,&atm,M)!=0 ||
        rt_kernel_phase_fourier(&ws,M,atm.betal_aer)!=0 ||
        rt_kernel_phase_fourier_pol(&ws,M,atm.gammal_aer)!=0 ||
        rt_kernel_phase_fourier_aerosol_full(&ws,M,atm.alphal_aer,atm.zetal_aer)!=0)
        return 5;

    /* Emulate one inserted output-only ordinate. */
    atm.gb[+3]=0.0;
    atm.gb[-3]=0.0;

    const int dirs=2*N_MU+1;
    const size_t n=(size_t)(N_LAYERS+1)*(size_t)dirs;
    double *mem=(double*)calloc(11*n,sizeof(double));
    if (!mem) return 6;
    double *I=mem+0*n,*Q=mem+1*n,*U=mem+2*n;
    double *v_on_I=mem+3*n,*v_on_Q=mem+4*n,*v_on_U=mem+5*n;
    double *v_off_I=mem+6*n,*v_off_Q=mem+7*n,*v_off_U=mem+8*n;
    double *s_on=mem+9*n,*s_off=mem+10*n;
    for (size_t i=0;i<n;++i) {
        I[i]=0.2+sin(0.0013*(double)(i+1));
        Q[i]=0.03*cos(0.0021*(double)(i+3));
        U[i]=0.02*sin(0.0017*(double)(i+5));
    }

    unsetenv("OCRT_SOS_ZERO_COL_SKIP_OFF");
    if (rt_sos_operator_apply_vector(&atm,M,&ws,I,Q,U,v_on_I,v_on_Q,v_on_U)!=0) return 7;
    if (rt_sos_operator_apply_scalar(&atm,M,&ws,I,s_on)!=0) return 8;
    setenv("OCRT_SOS_ZERO_COL_SKIP_OFF","1",1);
    if (rt_sos_operator_apply_vector(&atm,M,&ws,I,Q,U,v_off_I,v_off_Q,v_off_U)!=0) return 9;
    if (rt_sos_operator_apply_scalar(&atm,M,&ws,I,s_off)!=0) return 10;

    if (!same_bits(v_on_I,v_off_I,n,"vector I") ||
        !same_bits(v_on_Q,v_off_Q,n,"vector Q") ||
        !same_bits(v_on_U,v_off_U,n,"vector U") ||
        !same_bits(s_on,s_off,n,"scalar I")) return 11;

    printf("PASS zero-weight source-column skip: scalar/vector bit-identical values=%zu\n",n);
    free(mem);
    rt_legendre_workspace_free(&ws);
    rt_atm_free(&atm);
    return 0;
}
