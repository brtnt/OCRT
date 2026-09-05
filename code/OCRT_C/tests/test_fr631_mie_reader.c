#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "shared/mie_io.h"

static int find_wl(const double *x, int n, double q) {
    for (int i=0;i<n;i++) if (fabs(x[i]-q)<5e-10) return i;
    return -1;
}
static int find_angle(const double *x, int n, double q) {
    for (int i=0;i<n;i++) if (fabs(x[i]-q)<5e-10) return i;
    return -1;
}
int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr,"usage: %s file.mie\n",argv[0]); return 2; }
    mie_data_t m={0}; aerosol_phase_interp_t a={0};
    if (read_mie_file(argv[1],&m)!=0) { fprintf(stderr,"read failed\n"); return 3; }
    if (m.n_ang!=631 || m.n_wl!=771 || m.n_phase_wl!=771) {
        fprintf(stderr,"shape %d %d %d\n",m.n_ang,m.n_wl,m.n_phase_wl); return 4;
    }
    if (fabs(m.wavelengths[0]-0.3300)>1e-12 || fabs(m.wavelengths[770]-1.1000)>1e-12 ||
        fabs(m.phase_wavelengths[0]-0.3300)>1e-12 || fabs(m.phase_wavelengths[770]-1.1000)>1e-12 ||
        fabs(m.angles[0]-180.0)>1e-12 || fabs(m.angles[630]-0.0)>1e-12) return 5;
    if (build_aerosol_interpolators(&m,0.443,&a)!=0) return 6;
    int iw=find_wl(m.phase_wavelengths,m.n_phase_wl,0.443);
    const double tests[]={0.0,0.005,0.2,1.0,5.0,20.0,90.0,180.0};
    for (size_t j=0;j<sizeof(tests)/sizeof(tests[0]);j++) {
        int ia=find_angle(m.angles,m.n_ang,tests[j]);
        if (ia<0 || iw<0) return 7;
        double p11,p12,p33; eval_aerosol_phase(&a,tests[j],&p11,&p12,&p33);
        double e11=BLOCK_AT(m.P11,m.n_phase_wl,ia,iw);
        double e12=BLOCK_AT(m.P12,m.n_phase_wl,ia,iw);
        double e33=BLOCK_AT(m.P33,m.n_phase_wl,ia,iw);
        if (fabs(p11-e11)>1e-12*fmax(1.0,fabs(e11)) ||
            fabs(p12-e12)>1e-12*fmax(1.0,fabs(e12)) ||
            fabs(p33-e33)>1e-12*fmax(1.0,fabs(e33))) {
            fprintf(stderr,"exact-node mismatch theta=%.6f\n",tests[j]); return 8;
        }
    }
    printf("PASS %s n_ang=%d n_wl=%d n_phase_wl=%d bb/b=%.12g g=%.12g\n",
           argv[1],m.n_ang,m.n_wl,m.n_phase_wl,a.bb_b_ratio,a.g_asym);
    aerosol_phase_interp_free(&a); mie_data_free(&m); return 0;
}
