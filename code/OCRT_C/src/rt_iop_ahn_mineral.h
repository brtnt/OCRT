/* rt_iop_ahn_mineral.h
 *
 * Ahn (1992/1999) four-species mineral/TSM optical-property adapter for OCRT.
 *
 * Given dry-weight TSM concentration C [g m^-3], species and wavelength:
 *   a(lambda)  = a*(lambda) C
 *   b(lambda)  = b*(lambda) C
 *   bb(lambda) = b(lambda) [bb/b](lambda)
 *
 * a* and b* are read from astarmin_*.txt / bstarmin_*.txt when present.
 * The same values can be recovered from the matching AHN .mie bulk table as
 * a*=Extinct_Co-Scatter_Co and b*=Scatter_Co. bb/b is integrated from P11.
 */
#ifndef RT_IOP_AHN_MINERAL_H
#define RT_IOP_AHN_MINERAL_H

#include "rt_water_iop.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    AHN_RED_CLAY        = 0,
    AHN_BROWN_EARTH     = 1,
    AHN_YELLOW_CLAY     = 2,
    AHN_CALCAREOUS_SAND = 3,
    AHN_SPECIES_COUNT   = 4
} ahn_species_t;

int rt_iop_ahn_mineral_init(const char *data_dir);
void rt_iop_ahn_mineral_free(void);
int rt_iop_ahn_mineral_ready(void);

ahn_species_t rt_iop_ahn_species_from_name(const char *name);
const char *rt_iop_ahn_species_name(ahn_species_t species);
const char *rt_iop_ahn_mineral_phase_path(ahn_species_t species);

int rt_iop_ahn_mineral_eval(double lambda_nm, double tsm_g_m3,
                            ahn_species_t species, rt_iop_t *out,
                            const char **mie_phase_path);

#ifdef __cplusplus
}
#endif
#endif
