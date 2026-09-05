#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rt_angular_basis.h"

int main(void) {
    const int L = 80, n = 1, dirs = 3;
    const double xs[] = {-1.0, -0.9961946980917455, -0.5, 0.0, 0.5, 0.9961946980917455, 1.0};
    double mu_storage[3];
    double *mu = &mu_storage[1];
    double *sstore = calloc((L+1)*dirs, sizeof(double));
    double *rstore = calloc((L+1)*dirs, sizeof(double));
    double *tstore = calloc((L+1)*dirs, sizeof(double));
    double **s = calloc(L+1, sizeof(double*));
    double **r = calloc(L+1, sizeof(double*));
    double **t = calloc(L+1, sizeof(double*));
    double *se = calloc(L+1, sizeof(double));
    double *re = calloc(L+1, sizeof(double));
    double *te = calloc(L+1, sizeof(double));
    if (!sstore || !rstore || !tstore || !s || !r || !t || !se || !re || !te) return 2;
    for (int l=0;l<=L;l++){s[l]=sstore+l*dirs+1;r[l]=rstore+l*dirs+1;t[l]=tstore+l*dirs+1;}
    int fail=0;
    for (size_t ix=0; ix<sizeof(xs)/sizeof(xs[0]); ix++) {
        for (int j=-1;j<=1;j++) mu[j]=xs[ix];
        for (int m=0;m<=16;m++) {
            memset(sstore,0,(L+1)*dirs*sizeof(double));
            memset(rstore,0,(L+1)*dirs*sizeof(double));
            memset(tstore,0,(L+1)*dirs*sizeof(double));
            if (rt_angular_basis_build_scalar(n,L,m,mu,s) ||
                rt_angular_basis_build_spin2(n,L,m,mu,r,t) ||
                rt_angular_basis_eval_scalar(L,m,xs[ix],se) ||
                rt_angular_basis_eval_spin2(L,m,xs[ix],re,te)) return 3;
            for (int l=0;l<=L;l++) {
                if (memcmp(&s[l][0],&se[l],sizeof(double)) ||
                    memcmp(&r[l][0],&re[l],sizeof(double)) ||
                    memcmp(&t[l][0],&te[l],sizeof(double))) {
                    fprintf(stderr,"mismatch x=%.17g m=%d l=%d %.17g %.17g | %.17g %.17g | %.17g %.17g\n",
                            xs[ix],m,l,s[l][0],se[l],r[l][0],re[l],t[l][0],te[l]);
                    fail=1; goto done;
                }
            }
        }
    }
done:
    free(sstore);free(rstore);free(tstore);free(s);free(r);free(t);free(se);free(re);free(te);
    if (!fail) puts("PASS: single-direction basis is bit-identical to table builder");
    return fail;
}
