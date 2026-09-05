#include "rt_water_iop.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int close_rel(double got, double ref, double rtol, double atol) {
    return fabs(got-ref) <= atol + rtol*fabs(ref);
}

int main(int argc, char **argv) {
    const char *path = argc > 1 ? argv[1] : "inputs/water_iop/water_coef_z09_1nm.txt";
    rt_water_iop_lut_t lut = {0};
    if (rt_water_iop_lut_load(path, &lut) != 0) return 2;
    if (lut.n != 2250 || lut.lambda_nm[0] != 200.0 || lut.lambda_nm[lut.n-1] != 2449.0) {
        fprintf(stderr, "bad table shape/range: n=%d range=%.1f..%.1f\n",
                lut.n, lut.lambda_nm[0], lut.lambda_nm[lut.n-1]);
        rt_water_iop_lut_free(&lut); return 3;
    }
    struct row { double wl, aw, bw; } refs[] = {
        {380.0, 1.155917e-02, 8.383890e-03},
        {412.0, 4.553622e-03, 5.909969e-03},
        {443.0, 7.067186e-03, 4.331984e-03},
        {490.0, 1.502097e-02, 2.824033e-03},
        {510.0, 3.250344e-02, 2.385585e-03},
        {555.0, 5.963466e-02, 1.672631e-03},
        {620.0, 2.755385e-01, 1.053406e-03},
        {660.0, 4.100068e-01, 8.122936e-04},
        {680.0, 4.653764e-01, 7.176552e-04},
        {709.0, 7.966703e-01, 6.036142e-04},
        {745.0, 2.835900e+00, 4.917794e-04},
        {865.0, 4.605200e+00, 2.656188e-04},
        {970.0, 4.799918e+01, 1.658642e-04},
        {1000.0,4.070878e+01, 1.463765e-04},
        {1100.0,1.908413e+01, 9.904724e-05},
    };
    int fail = 0;
    for (size_t i=0; i<sizeof(refs)/sizeof(refs[0]); ++i) {
        int ea=0, eb=0;
        double aw=rt_water_iop_aw(&lut, refs[i].wl, &ea);
        double bw=rt_water_iop_bw(&lut, refs[i].wl, &eb);
        double bbw=rt_water_iop_bbw(&lut, refs[i].wl, NULL);
        if (ea || eb || !close_rel(aw, refs[i].aw, 1e-13, 1e-15) ||
            !close_rel(bw, refs[i].bw, 1e-13, 1e-15) ||
            !close_rel(bbw, 0.5*refs[i].bw, 1e-13, 1e-15)) {
            fprintf(stderr, "FAIL %.1f aw=%.12e ref=%.12e bw=%.12e ref=%.12e ext=%d/%d\n",
                    refs[i].wl, aw, refs[i].aw, bw, refs[i].bw, ea, eb);
            fail++;
        }
    }
    rt_water_iop_lut_free(&lut);
    if (fail) return 4;
    puts("PASS: pure-water Z09 table 2250 rows, key bands exact, bb=0.5*b");
    return 0;
}
