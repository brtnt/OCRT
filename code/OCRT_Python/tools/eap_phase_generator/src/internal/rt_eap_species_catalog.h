/* rt_eap_species_catalog.h -- internal only.  Frozen 17-species EAP catalog. */
#ifndef OCRT_RT_EAP_SPECIES_CATALOG_H
#define OCRT_RT_EAP_SPECIES_CATALOG_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    const char   *name;
    const double *deff_um;      /* frozen representative effective diameter */
    const double *vs;           /* shell (chloroplast / chromatoplasm) volume fraction */
    const double *nshell;       /* [nwl] real refractive index, relative to water */
    const double *kshell;       /* [nwl] imaginary part, positive magnitude */
    const double *ncore;
    const double *kcore;
} rt_eap_species_entry_t;

int    rt_eap_catalog_count(void);
double rt_eap_catalog_wl_lo(void);
double rt_eap_catalog_wl_hi(void);
const rt_eap_species_entry_t *rt_eap_catalog_entry(unsigned id);
/* Linear interpolation of the four refractive-index spectra at one wavelength.
 * Returns 0 on success, -1 if the wavelength is outside the catalog range. */
int rt_eap_catalog_ri(unsigned id, double wl_nm,
                      double *nsh, double *ksh, double *nco, double *kco);
#ifdef __cplusplus
}
#endif
#endif
