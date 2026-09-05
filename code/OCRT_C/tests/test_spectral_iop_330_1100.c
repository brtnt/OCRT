#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rt_spectral_contract.h"
#include "rt_water_iop.h"
#include "rt_iop_organic.h"
#include "rt_iop_ahn_mineral.h"
#include "shared/mie_io.h"

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "usage: %s WATER_IOP_DIR TSM_AHN_DIR\n", argv[0]);
        return 2;
    }
    char water_path[2048], psi_path[2048];
    snprintf(water_path, sizeof water_path, "%s/water_coef_z09_1nm.txt", argv[1]);
    snprintf(psi_path, sizeof psi_path, "%s/psi_T_rottgers2014_OCRT.txt", argv[1]);
    rt_water_iop_lut_t w = {0};
    rt_water_iop_psi_T_lut_t p = {0};
    if (rt_water_iop_lut_load(water_path, &w) != 0 ||
        rt_water_iop_psi_T_load(psi_path, &p) != 0 ||
        rt_iop_organic_init_selected(argv[1], ORGANIC_PHYTO_MICRO, 0) != 0 ||
        rt_iop_ahn_mineral_init(argv[2]) != 0) return 3;

    char det_path[2048], min_path[2048];
    snprintf(det_path, sizeof det_path, "%s/Detritus_Stramski2001.mie", argv[1]);
    snprintf(min_path, sizeof min_path, "%s/Red_clay_AHN.mie", argv[2]);
    mie_data_t det_mie = {0}, min_mie = {0};
    if (read_mie_file(det_path, &det_mie) != 0 ||
        read_mie_file(min_path, &min_mie) != 0) return 5;

    const double wavelengths[] = {330,340,349,350,443,555,750,800,865,940,1100};
    puts("wavelength_nm,a_w_T15,b_w,bb_w,a_phyto,b_det,bb_det,a_min,b_min,bb_min");
    for (size_t i = 0; i < sizeof wavelengths / sizeof wavelengths[0]; ++i) {
        double wl = wavelengths[i];
        int ex = 0;
        double aw = rt_water_iop_aw_T(&w, &p, wl, 15.0, &ex);
        double bw = rt_water_iop_bw(&w, wl, &ex);
        double aph = 0.0;
        rt_iop_t det = {0}, min = {0};
        const char *phase = NULL;
        aerosol_phase_interp_t det_w = {0}, det_550 = {0}, min_w = {0};
        const int det_supported = rt_iop_organic_detritus_wavelength_supported(wl);
        if (ex || !rt_spectral_wavelength_supported(wl) ||
            build_aerosol_interpolators(&min_mie, wl * 1.0e-3, &min_w) != 0 ||
            rt_iop_eap_phyto_absorption_eval(wl, 1.0, ORGANIC_PHYTO_MICRO, &aph) != 0 ||
            rt_iop_ahn_mineral_eval(wl, 1.0, AHN_RED_CLAY, &min, &phase) != 0) return 4;
        if (det_supported) {
            if (build_aerosol_interpolators(&det_mie, wl * 1.0e-3, &det_w) != 0 ||
                build_aerosol_interpolators(&det_mie, 0.550, &det_550) != 0 ||
                rt_iop_detritus_eval_with_phase_ratios(
                    wl, 1.0, 0.0, ORGANIC_DETRITUS_SLOPE_DEFAULT,
                    det_w.bb_b_ratio, det_550.bb_b_ratio, &det) != 0) return 4;
        } else {
            rt_iop_t probe = {0};
            if (rt_iop_detritus_eval(wl, 1.0, 0.0,
                                     ORGANIC_DETRITUS_SLOPE_DEFAULT, &probe) != -7)
                return 6;
            det.b = det.bb = NAN;
        }
        min.bb = min.b * min_w.bb_b_ratio;
        aerosol_phase_interp_free(&det_w);
        aerosol_phase_interp_free(&det_550);
        aerosol_phase_interp_free(&min_w);
        printf("%.0f,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
               wl, aw, bw, 0.5*bw, aph, det.b, det.bb, min.a, min.b, min.bb);
    }
    mie_data_free(&det_mie);
    mie_data_free(&min_mie);
    rt_iop_ahn_mineral_free();
    rt_iop_organic_free();
    rt_water_iop_psi_T_free(&p);
    rt_water_iop_lut_free(&w);
    return 0;
}
