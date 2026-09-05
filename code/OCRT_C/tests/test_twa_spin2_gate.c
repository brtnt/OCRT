#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define _USE_MATH_DEFINES
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "rt_air_water_coupling.h"
#include "rt_quadrature.h"
static void run(double wind, const double* mua, int na, double* outI, double* outQ, double* outU){
    const int nmu=48, mmax=2;
    double muw[48], ww[48];
    rt_quadrature_gauss_legendre_pos(nmu, muw, ww);
    size_t sz=(size_t)(mmax+1)*nmu;
    double *I=calloc(sz,sizeof(double)), *Q=calloc(sz,sizeof(double)), *U=calloc(sz,sizeof(double));
    for(int j=0;j<nmu;++j){ I[0*nmu+j]=1.0; Q[2*nmu+j]=-0.5; U[2*nmu+j]=+0.5; }
    size_t sa=(size_t)(mmax+1)*na;
    double *Ia=calloc(sa,sizeof(double)), *Qa=calloc(sa,sizeof(double)), *Ua=calloc(sa,sizeof(double));
    int rc=rt_air_water_couple_water_to_atm(I,Q,U,muw,NULL,nmu,mmax,1.34,1,wind,0,mua,na,NULL,Ia,Qa,Ua);
    if(rc!=0){ fprintf(stderr,"couple rc=%d wind=%g\n",rc,wind); exit(1);}    
    for(int a=0;a<na;++a){ outI[a]=Ia[0*na+a]; outQ[a]=Qa[2*na+a]; outU[a]=Ua[2*na+a]; }
    free(I);free(Q);free(U);free(Ia);free(Qa);free(Ua);
}
int main(void){
    double th[]={2,5,10,15,20,30,40}; int na=7; double mua[7];
    for(int i=0;i<na;++i) mua[i]=cos(th[i]*M_PI/180.0);
    double I3[7],Q3[7],U3[7],I0[7],Q0[7],U0[7];
    run(3.0, mua,na, I3,Q3,U3);
    run(0.0, mua,na, I0,Q0,U0);
    printf("%6s %10s %10s %10s | ratio rough/flat: %8s %8s %8s\n","VZAair","I0_rough","Q2_rough","U2_rough","T_I","T_Q2","T_U2");
    for(int i=0;i<na;++i){
        printf("%6.1f %10.3e %10.3e %10.3e |                  %8.3f %8.3f %8.3f\n",
            th[i],I3[i],Q3[i],U3[i], I3[i]/I0[i], Q3[i]/Q0[i], U3[i]/U0[i]);
    }
    return 0;
}
