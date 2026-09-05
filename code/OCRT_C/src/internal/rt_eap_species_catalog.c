#include "rt_eap_species_catalog.h"
#include <math.h>
#include "../generated/rt_eap_species_data.inc"

int rt_eap_catalog_count(void)
{
    return (int)(sizeof k_eap_catalog / sizeof k_eap_catalog[0]);
}
double rt_eap_catalog_wl_lo(void) { return RT_EAP_WL_LO; }
double rt_eap_catalog_wl_hi(void) { return RT_EAP_WL_HI; }

const rt_eap_species_entry_t *rt_eap_catalog_entry(unsigned id)
{
    if ((int)id >= rt_eap_catalog_count()) return NULL;
    return &k_eap_catalog[id];
}

int rt_eap_catalog_ri(unsigned id, double wl_nm,
                      double *nsh, double *ksh, double *nco, double *kco)
{
    const rt_eap_species_entry_t *e = rt_eap_catalog_entry(id);
    if (!e || !isfinite(wl_nm)) return -1;
    if (wl_nm < RT_EAP_WL_LO || wl_nm > RT_EAP_WL_HI) return -1;
    double f = (wl_nm - RT_EAP_WL_LO) / RT_EAP_WL_STEP;
    int i = (int)floor(f);
    if (i < 0) i = 0;
    if (i > RT_EAP_NWL - 2) i = RT_EAP_NWL - 2;
    double t = f - (double)i;
    *nsh = e->nshell[i] + t * (e->nshell[i + 1] - e->nshell[i]);
    *ksh = e->kshell[i] + t * (e->kshell[i + 1] - e->kshell[i]);
    *nco = e->ncore[i]  + t * (e->ncore[i + 1]  - e->ncore[i]);
    *kco = e->kcore[i]  + t * (e->kcore[i + 1]  - e->kcore[i]);
    return 0;
}
