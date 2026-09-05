
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rt_atm.h"
#include "rt_kernel.h"
#include "rt_sos_operator.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void ocrt_test_add_diffuse_top_primary(const rt_atm_t *atm,
                                       const rt_legendre_workspace_t *ws,
                                       int m,
                                       const double *ext_I,
                                       const double *ext_Q,
                                       const double *ext_U,
                                       const double *ext_mu, int ext_n,
                                       double *src_i, double *src_q,
                                       double *src_u);

static double max_physical_diff(const double *a, const double *b,
                                int n_layers, int n_mu) {
    const int dirs = 2*n_mu+1;
    double mx = 0.0;
    for (int k=0; k<=n_layers; ++k) {
        for (int j=-n_mu; j<=n_mu; ++j) {
            if (j == 0) continue;
            const int idx = k*dirs + (j+n_mu);
            const double d = fabs(a[idx]-b[idx]);
            if (d > mx) mx = d;
        }
    }
    return mx;
}

static int run_regime(double tau_a) {
    enum { N_MU=8, N_LAYERS=5, LMAX=12 };
    rt_atm_t atm={0};
    rt_legendre_workspace_t ws={0};
    double beta[LMAX+1], gamma[LMAX+1], alpha[LMAX+1], zeta[LMAX+1];

    for (int l=0; l<=LMAX; ++l) {
        const double d=(double)l;
        beta[l]=(l==0)?1.0:(2.0*d+1.0)*pow(0.73,d);
        gamma[l]=(l<2)?0.0:0.13*pow(-0.61,d-2.0);
        alpha[l]=(l<2)?0.0:0.17*pow(0.57,d-2.0);
        zeta[l]=(l<2)?0.0:0.11*pow(-0.49,d-2.0);
    }

    if (rt_atm_alloc(&atm,N_LAYERS,N_MU)!=0) return 2;
    if (rt_atm_build_aerosol_rayleigh(&atm,0.23,0.039,
            tau_a,0.92,LMAX,beta,gamma,alpha,zeta,
            cos(40.0*M_PI/180.0),RT_RAYLEIGH_MODEL_BODHAINE_1999,2.0)!=0) return 3;
    if (rt_legendre_workspace_alloc(&ws,N_MU,LMAX)!=0) return 4;

    const int dirs=2*N_MU+1;
    const size_t nf=(size_t)(N_LAYERS+1)*dirs;
    double *buf=calloc(9*nf,sizeof(double));
    if(!buf)return 5;
    double *prevI=buf+0*nf,*prevQ=buf+1*nf,*prevU=buf+2*nf;
    double *opI=buf+3*nf,*opQ=buf+4*nf,*opU=buf+5*nf;
    double *dtI=buf+6*nf,*dtQ=buf+7*nf,*dtU=buf+8*nf;
    double extI[N_MU],extQ[N_MU],extU[N_MU],mu[N_MU];

    double global_max=0.0;
    for(int m=0;m<=4;m++){
        if(rt_legendre_compute_pol(&ws,&atm,m)!=0 ||
           rt_kernel_phase_fourier(&ws,m,atm.betal_aer)!=0 ||
           rt_kernel_phase_fourier_pol(&ws,m,atm.gammal_aer)!=0 ||
           rt_kernel_phase_fourier_aerosol_full(&ws,m,atm.alphal_aer,atm.zetal_aer)!=0 ||
           rt_sos_operator_prepare(&atm,m,&ws)!=0) return 6;

        for(int c0=0;c0<N_MU;c0++){
          for(int basis=0;basis<3;basis++){
            memset(extI,0,sizeof extI); memset(extQ,0,sizeof extQ); memset(extU,0,sizeof extU);
            for(int c=0;c<N_MU;c++)mu[c]=atm.rm[+(c+1)];
            if(basis==0)extI[c0]=0.73;
            if(basis==1)extQ[c0]=-0.21;
            if(basis==2)extU[c0]=0.17;
            memset(prevI,0,nf*sizeof(double)); memset(prevQ,0,nf*sizeof(double));
            memset(prevU,0,nf*sizeof(double)); memset(dtI,0,nf*sizeof(double));
            memset(dtQ,0,nf*sizeof(double)); memset(dtU,0,nf*sizeof(double));

            const int c=c0+1;
            for(int k=0;k<=N_LAYERS;k++){
                const int idx=k*dirs+(N_MU-c); /* incoming downward -c */
                const double att=exp(-atm.h[k]/atm.rm[+c]);
                prevI[idx]=att*extI[c0];
                prevQ[idx]=att*extQ[c0];
                prevU[idx]=att*extU[c0];
            }
            if(rt_sos_operator_apply_vector(&atm,m,&ws,prevI,prevQ,prevU,
                                            opI,opQ,opU)!=0)return 7;
            ocrt_test_add_diffuse_top_primary(&atm,&ws,m,extI,extQ,extU,mu,N_MU,
                                              dtI,dtQ,dtU);
            const double dI=max_physical_diff(opI,dtI,N_LAYERS,N_MU);
            const double dQ=max_physical_diff(opQ,dtQ,N_LAYERS,N_MU);
            const double dU=max_physical_diff(opU,dtU,N_LAYERS,N_MU);
            if(dI>global_max)global_max=dI;
            if(dQ>global_max)global_max=dQ;
            if(dU>global_max)global_max=dU;
          }
        }
    }

    free(buf); rt_legendre_workspace_free(&ws); rt_atm_free(&atm);
    if(global_max > 1.0e-13) {
        fprintf(stderr,"closure failure tau_a=%.3g max_abs=%.17g\n",tau_a,global_max);
        return 10;
    }
    printf("PASS diffuse-top/source-operator closure tau_a=%.3g max_abs=%.3e\n",
           tau_a,global_max);
    return 0;
}

int main(void) {
    int rc=run_regime(0.0);
    if(rc!=0)return rc;
    rc=run_regime(0.17);
    if(rc!=0)return rc;
    return 0;
}
