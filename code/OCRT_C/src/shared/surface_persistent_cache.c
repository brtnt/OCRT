/* ============================================================================
 * Versioned persistent Cox-Munk/interface operator cache.
 * DOC-REF: docs/PERSISTENT_COXMUNK_INTERFACE_CACHE_KO_2026-08-25.md
 *
 * File format is explicitly little-endian and includes the full immutable key,
 * two independent 64-bit key hashes, and two payload hashes.  A hash collision
 * cannot create a false hit because every scalar and both direction arrays are
 * compared bit-for-bit after opening the candidate file.
 *
 * The cache is intentionally coordination-free.  Production rows use
 * READONLY/REQUIRED and only open immutable files. BUILD is a separate priming
 * step; no lock, sleep, poll, or dependency on another simulation row exists.
 * ========================================================================== */
#include "surface_persistent_cache.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The persistent-cache implementation deliberately uses only ISO C file I/O.
 * The caller supplies an existing cache directory; the cache layer does not
 * create directories, acquire OS locks, query process IDs, or call a
 * platform-specific filesystem API.  Forward slashes are accepted by the
 * supported C runtimes and keep one source path on all platforms. */
#define OCRT_PCACHE_PATH_CAP 4096
#define OCRT_PCACHE_SEP "/"

#define OCRT_PCACHE_MAGIC "OCRTSC01"
#define OCRT_PCACHE_MAGIC_N 8u

_Static_assert(sizeof(double) == 8, "persistent surface cache requires IEEE-754 binary64 doubles");

typedef struct {
    uint64_t a;
    uint64_t b;
} ocrt_hash128_t;

typedef struct {
    int initialized;
    ocrt_surface_pcache_mode_t mode;
    int trace;
    char dir[OCRT_PCACHE_PATH_CAP];
} ocrt_surface_pcache_config_t;

/* Immutable after the explicit pre-parallel initialization in OCRT main. */
static ocrt_surface_pcache_config_t g_pcache_cfg;
static _Thread_local unsigned long g_pcache_tmp_counter;

static int ocrt_host_little_endian_(void) {
    const uint16_t x = 1u;
    return *((const unsigned char *)&x) == 1u;
}

static void ocrt_hash_init_(ocrt_hash128_t *h) {
    h->a = UINT64_C(1469598103934665603);
    h->b = UINT64_C(7809847782465536322);
}

static void ocrt_hash_bytes_(ocrt_hash128_t *h, const void *data, size_t n) {
    const unsigned char *p = (const unsigned char *)data;
    for (size_t i = 0; i < n; ++i) {
        h->a ^= (uint64_t)p[i];
        h->a *= UINT64_C(1099511628211);
        h->b ^= (uint64_t)(p[i] + (unsigned char)(i * 17u));
        h->b *= UINT64_C(14029467366897019727);
        h->b ^= h->b >> 29;
    }
}

static void ocrt_hash_u32_(ocrt_hash128_t *h, uint32_t v) {
    unsigned char b[4];
    b[0] = (unsigned char)(v);
    b[1] = (unsigned char)(v >> 8);
    b[2] = (unsigned char)(v >> 16);
    b[3] = (unsigned char)(v >> 24);
    ocrt_hash_bytes_(h, b, sizeof b);
}

static void ocrt_hash_u64_(ocrt_hash128_t *h, uint64_t v) {
    unsigned char b[8];
    for (int i = 0; i < 8; ++i) b[i] = (unsigned char)(v >> (8 * i));
    ocrt_hash_bytes_(h, b, sizeof b);
}

static uint64_t ocrt_double_bits_(double v) {
    uint64_t u = 0;
    memcpy(&u, &v, sizeof u);
    return u;
}

static void ocrt_hash_double_(ocrt_hash128_t *h, double v) {
    ocrt_hash_u64_(h, ocrt_double_bits_(v));
}

static void ocrt_pcache_hash_key_(const ocrt_surface_pcache_key_t *k,
                                  ocrt_hash128_t *h) {
    static const char contract[] =
        "OCRT_SURFACE_OPERATOR_CACHE|CM54_SANCER|IQU|TWA_QFOLD_V1";
    ocrt_hash_init_(h);
    ocrt_hash_bytes_(h, contract, sizeof contract - 1u);
    ocrt_hash_u32_(h, OCRT_SURFACE_PCACHE_FORMAT_VERSION);
    ocrt_hash_u32_(h, OCRT_SURFACE_PCACHE_ABI_VERSION);
    ocrt_hash_u32_(h, OCRT_SURFACE_PCACHE_PHYSICS_VERSION);
    ocrt_hash_u32_(h, OCRT_SURFACE_PCACHE_QFOLD_VERSION);
    ocrt_hash_u32_(h, k->operator_kind);
    ocrt_hash_u32_(h, (uint32_t)k->mode_first);
    ocrt_hash_u32_(h, (uint32_t)k->mode_count);
    ocrt_hash_u32_(h, (uint32_t)k->n_o);
    ocrt_hash_u32_(h, (uint32_t)k->n_i);
    ocrt_hash_u32_(h, (uint32_t)k->n_phi_quad);
    ocrt_hash_u32_(h, (uint32_t)k->sigma_type);
    ocrt_hash_u32_(h, (uint32_t)k->q_convention);
    ocrt_hash_double_(h, k->wind_speed);
    ocrt_hash_double_(h, k->n_water);
    for (int i = 0; i < k->n_o; ++i) ocrt_hash_double_(h, k->mu_o[i]);
    for (int i = 0; i < k->n_i; ++i) ocrt_hash_double_(h, k->mu_i[i]);
}

static void ocrt_pcache_hash_payload_(const double *v, size_t n,
                                      ocrt_hash128_t *h) {
    ocrt_hash_init_(h);
    for (size_t i = 0; i < n; ++i) ocrt_hash_double_(h, v[i]);
}

static const char *ocrt_pcache_operator_name_(uint32_t op) {
    switch ((ocrt_surface_pcache_operator_t)op) {
        case OCRT_SURFACE_PCACHE_OP_R_AA_RAW_ALLM: return "raa";
        case OCRT_SURFACE_PCACHE_OP_R_WW_RAW_ALLM: return "rww";
        case OCRT_SURFACE_PCACHE_OP_T_AW_FINAL_M:  return "taw";
        case OCRT_SURFACE_PCACHE_OP_T_WA_RAW_ALLM: return "twa";
        default: return "unknown";
    }
}

static int ocrt_pcache_key_valid_(const ocrt_surface_pcache_key_t *k,
                                  size_t payload_count) {
    if (!k || !k->mu_o || !k->mu_i || payload_count == 0u) return 0;
    if (k->operator_kind < OCRT_SURFACE_PCACHE_OP_R_AA_RAW_ALLM ||
        k->operator_kind > OCRT_SURFACE_PCACHE_OP_T_WA_RAW_ALLM) return 0;
    if (k->mode_first < 0 || k->mode_count < 1 ||
        k->n_o < 1 || k->n_i < 1 || k->n_phi_quad < 4) return 0;
    const size_t pair9 = (size_t)k->n_o * (size_t)k->n_i * 9u;
    return payload_count == pair9 * (size_t)k->mode_count;
}

static int ocrt_parse_mode_(const char *s, ocrt_surface_pcache_mode_t *mode) {
    if (!s || !s[0] || strcmp(s, "off") == 0 || strcmp(s, "0") == 0) {
        *mode = OCRT_SURFACE_PCACHE_OFF; return 0;
    }
    if (strcmp(s, "readonly") == 0 || strcmp(s, "read-only") == 0 ||
        strcmp(s, "ro") == 0 || strcmp(s, "optional") == 0) {
        *mode = OCRT_SURFACE_PCACHE_READONLY; return 0;
    }
    if (strcmp(s, "required") == 0 || strcmp(s, "require") == 0 ||
        strcmp(s, "strict") == 0) {
        *mode = OCRT_SURFACE_PCACHE_REQUIRED; return 0;
    }
    if (strcmp(s, "build") == 0 || strcmp(s, "prime") == 0 ||
        strcmp(s, "prebuild") == 0) {
        *mode = OCRT_SURFACE_PCACHE_BUILD; return 0;
    }
    return -1;
}

int ocrt_surface_pcache_init_from_env(void) {
    if (g_pcache_cfg.initialized) return 0;
    memset(&g_pcache_cfg, 0, sizeof g_pcache_cfg);
    const char *m = getenv("OCRT_SURFACE_PERSIST_CACHE_MODE");
    if (ocrt_parse_mode_(m, &g_pcache_cfg.mode) != 0) {
        fprintf(stderr,
                "error: OCRT_SURFACE_PERSIST_CACHE_MODE must be "
                "off, readonly, required, or build (got '%s')\n",
                m ? m : "");
        return -1;
    }
    const char *t = getenv("OCRT_SURFACE_PERSIST_CACHE_TRACE");
    g_pcache_cfg.trace = (t && t[0] && strcmp(t, "0") != 0 &&
                          strcmp(t, "off") != 0 && strcmp(t, "false") != 0);
    if (g_pcache_cfg.mode != OCRT_SURFACE_PCACHE_OFF) {
        const char *d = getenv("OCRT_SURFACE_PERSIST_CACHE_DIR");
        if (!d || !d[0] || strlen(d) >= sizeof g_pcache_cfg.dir) {
            fprintf(stderr,
                    "error: persistent surface cache mode '%s' requires "
                    "OCRT_SURFACE_PERSIST_CACHE_DIR\n",
                    ocrt_surface_pcache_mode_name());
            return -1;
        }
        /* Directory creation/existence probing is intentionally outside the
         * C solver.  BUILD will fail loudly if the path is not writable;
         * READONLY/REQUIRED will report a cache miss. */
        memcpy(g_pcache_cfg.dir, d, strlen(d) + 1u);
    }
    g_pcache_cfg.initialized = 1;
    if (g_pcache_cfg.trace) {
        fprintf(stderr,
                "[surface-pcache] mode=%s dir=%s format=%u abi=%u physics=%u qfold=%u "
                "coordination=none\n",
                ocrt_surface_pcache_mode_name(),
                g_pcache_cfg.mode == OCRT_SURFACE_PCACHE_OFF ? "-" : g_pcache_cfg.dir,
                OCRT_SURFACE_PCACHE_FORMAT_VERSION,
                OCRT_SURFACE_PCACHE_ABI_VERSION,
                OCRT_SURFACE_PCACHE_PHYSICS_VERSION,
                OCRT_SURFACE_PCACHE_QFOLD_VERSION);
    }
    return 0;
}

ocrt_surface_pcache_mode_t ocrt_surface_pcache_mode(void) {
    if (!g_pcache_cfg.initialized && ocrt_surface_pcache_init_from_env() != 0)
        return OCRT_SURFACE_PCACHE_OFF;
    return g_pcache_cfg.mode;
}

const char *ocrt_surface_pcache_mode_name(void) {
    switch (g_pcache_cfg.mode) {
        case OCRT_SURFACE_PCACHE_OFF: return "off";
        case OCRT_SURFACE_PCACHE_READONLY: return "readonly";
        case OCRT_SURFACE_PCACHE_REQUIRED: return "required";
        case OCRT_SURFACE_PCACHE_BUILD: return "build";
        default: return "invalid";
    }
}

int ocrt_surface_pcache_batch_guard(void) {
    if (ocrt_surface_pcache_init_from_env() != 0) return -1;
    if (g_pcache_cfg.mode == OCRT_SURFACE_PCACHE_BUILD) {
        fprintf(stderr,
                "error: persistent surface cache BUILD mode is forbidden inside "
                "--batch-full-grid; prebuild caches before launching independent "
                "rows, then use mode=required or readonly\n");
        return -1;
    }
    return 0;
}

static int ocrt_pcache_make_path_(const ocrt_surface_pcache_key_t *key,
                                  char *path, size_t cap,
                                  ocrt_hash128_t *key_hash_out) {
    ocrt_hash128_t h;
    ocrt_pcache_hash_key_(key, &h);
    if (key_hash_out) *key_hash_out = h;
    size_t nd = strlen(g_pcache_cfg.dir);
    const int add_sep = nd > 0 && g_pcache_cfg.dir[nd - 1] != '/' &&
                        g_pcache_cfg.dir[nd - 1] != '\\';
    int n = snprintf(path, cap, "%s%socrt_surface_v%u_%s_%016llx%016llx.bin",
                     g_pcache_cfg.dir, add_sep ? OCRT_PCACHE_SEP : "",
                     OCRT_SURFACE_PCACHE_FORMAT_VERSION,
                     ocrt_pcache_operator_name_(key->operator_kind),
                     (unsigned long long)h.a, (unsigned long long)h.b);
    return (n >= 0 && (size_t)n < cap) ? 0 : -1;
}

static int ocrt_write_bytes_(FILE *fp, const void *p, size_t n) {
    return fwrite(p, 1u, n, fp) == n ? 0 : -1;
}
static int ocrt_read_bytes_(FILE *fp, void *p, size_t n) {
    return fread(p, 1u, n, fp) == n ? 0 : -1;
}
static int ocrt_write_u32_(FILE *fp, uint32_t v) {
    unsigned char b[4];
    b[0]=(unsigned char)v; b[1]=(unsigned char)(v>>8);
    b[2]=(unsigned char)(v>>16); b[3]=(unsigned char)(v>>24);
    return ocrt_write_bytes_(fp,b,4u);
}
static int ocrt_read_u32_(FILE *fp, uint32_t *v) {
    unsigned char b[4]; if (ocrt_read_bytes_(fp,b,4u)!=0) return -1;
    *v=(uint32_t)b[0]|((uint32_t)b[1]<<8)|((uint32_t)b[2]<<16)|((uint32_t)b[3]<<24);
    return 0;
}
static int ocrt_write_u64_(FILE *fp, uint64_t v) {
    unsigned char b[8]; for(int i=0;i<8;++i)b[i]=(unsigned char)(v>>(8*i));
    return ocrt_write_bytes_(fp,b,8u);
}
static int ocrt_read_u64_(FILE *fp, uint64_t *v) {
    unsigned char b[8]; if (ocrt_read_bytes_(fp,b,8u)!=0) return -1;
    uint64_t x=0; for(int i=0;i<8;++i)x|=(uint64_t)b[i]<<(8*i); *v=x; return 0;
}
static int ocrt_write_f64_(FILE *fp, double v) {
    return ocrt_write_u64_(fp, ocrt_double_bits_(v));
}
static int ocrt_read_f64_(FILE *fp, double *v) {
    uint64_t u; if (ocrt_read_u64_(fp,&u)!=0) return -1; memcpy(v,&u,sizeof u); return 0;
}

static int ocrt_write_f64_array_(FILE *fp, const double *v, size_t n) {
    if (ocrt_host_little_endian_())
        return fwrite(v, sizeof(double), n, fp) == n ? 0 : -1;
    for (size_t i=0;i<n;++i) if (ocrt_write_f64_(fp,v[i])!=0) return -1;
    return 0;
}
static int ocrt_read_f64_array_(FILE *fp, double *v, size_t n) {
    if (ocrt_host_little_endian_())
        return fread(v, sizeof(double), n, fp) == n ? 0 : -1;
    for (size_t i=0;i<n;++i) if (ocrt_read_f64_(fp,&v[i])!=0) return -1;
    return 0;
}

static int ocrt_required_or_miss_(const char *reason, const char *path) {
    if (g_pcache_cfg.trace || g_pcache_cfg.mode == OCRT_SURFACE_PCACHE_REQUIRED)
        fprintf(stderr, "[surface-pcache] %s: %s\n", reason, path ? path : "-");
    return g_pcache_cfg.mode == OCRT_SURFACE_PCACHE_REQUIRED ? -1 : 0;
}

int ocrt_surface_pcache_load(const ocrt_surface_pcache_key_t *key,
                             double *payload, size_t payload_count) {
    if (ocrt_surface_pcache_init_from_env() != 0) return -1;
    if (g_pcache_cfg.mode == OCRT_SURFACE_PCACHE_OFF) return 0;
    if (!ocrt_pcache_key_valid_(key, payload_count) || !payload) return -1;

    char path[OCRT_PCACHE_PATH_CAP]; ocrt_hash128_t expect_key_hash;
    if (ocrt_pcache_make_path_(key,path,sizeof path,&expect_key_hash)!=0) return -1;
    FILE *fp=fopen(path,"rb");
    if(!fp) return ocrt_required_or_miss_("miss",path);

    int ok=1; char magic[OCRT_PCACHE_MAGIC_N];
    uint32_t fmt=0,abi=0,phys=0,qfold=0,op=0;
    uint32_t mf=0,mc=0,no=0,ni=0,nphi=0,st=0,qc=0;
    uint64_t h1=0,h2=0,count=0,p1=0,p2=0;
    double ws=0.0,nw=0.0;
    if(ocrt_read_bytes_(fp,magic,sizeof magic)!=0 || memcmp(magic,OCRT_PCACHE_MAGIC,sizeof magic)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&fmt)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&abi)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&phys)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&qfold)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&op)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&mf)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&mc)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&no)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&ni)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&nphi)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&st)!=0)ok=0;
    if(ok&&ocrt_read_u32_(fp,&qc)!=0)ok=0;
    if(ok&&ocrt_read_f64_(fp,&ws)!=0)ok=0;
    if(ok&&ocrt_read_f64_(fp,&nw)!=0)ok=0;
    if(ok&&ocrt_read_u64_(fp,&h1)!=0)ok=0;
    if(ok&&ocrt_read_u64_(fp,&h2)!=0)ok=0;
    if(ok&&ocrt_read_u64_(fp,&count)!=0)ok=0;
    if(ok&&ocrt_read_u64_(fp,&p1)!=0)ok=0;
    if(ok&&ocrt_read_u64_(fp,&p2)!=0)ok=0;

    if(ok && (fmt!=OCRT_SURFACE_PCACHE_FORMAT_VERSION ||
              abi!=OCRT_SURFACE_PCACHE_ABI_VERSION ||
              phys!=OCRT_SURFACE_PCACHE_PHYSICS_VERSION ||
              qfold!=OCRT_SURFACE_PCACHE_QFOLD_VERSION ||
              op!=key->operator_kind || mf!=(uint32_t)key->mode_first ||
              mc!=(uint32_t)key->mode_count || no!=(uint32_t)key->n_o ||
              ni!=(uint32_t)key->n_i || nphi!=(uint32_t)key->n_phi_quad ||
              st!=(uint32_t)key->sigma_type || qc!=(uint32_t)key->q_convention ||
              ocrt_double_bits_(ws)!=ocrt_double_bits_(key->wind_speed) ||
              ocrt_double_bits_(nw)!=ocrt_double_bits_(key->n_water) ||
              h1!=expect_key_hash.a || h2!=expect_key_hash.b ||
              count!=(uint64_t)payload_count)) ok=0;

    for(int i=0;ok&&i<key->n_o;++i){double x; if(ocrt_read_f64_(fp,&x)!=0||ocrt_double_bits_(x)!=ocrt_double_bits_(key->mu_o[i]))ok=0;}
    for(int i=0;ok&&i<key->n_i;++i){double x; if(ocrt_read_f64_(fp,&x)!=0||ocrt_double_bits_(x)!=ocrt_double_bits_(key->mu_i[i]))ok=0;}
    if(ok&&ocrt_read_f64_array_(fp,payload,payload_count)!=0)ok=0;
    if(ok){
        ocrt_hash128_t ph; ocrt_pcache_hash_payload_(payload,payload_count,&ph);
        if(ph.a!=p1||ph.b!=p2)ok=0;
    }
    if(ok){int c=fgetc(fp); if(c!=EOF)ok=0;}
    fclose(fp);
    if(!ok){
        if(g_pcache_cfg.mode==OCRT_SURFACE_PCACHE_BUILD) remove(path);
        return ocrt_required_or_miss_("stale-or-corrupt",path);
    }
    if(g_pcache_cfg.trace)
        fprintf(stderr,"[surface-pcache] hit op=%s modes=%d+%d no=%d ni=%d file=%s\n",
                ocrt_pcache_operator_name_(key->operator_kind),key->mode_first,key->mode_count,
                key->n_o,key->n_i,path);
    return 1;
}

int ocrt_surface_pcache_store(const ocrt_surface_pcache_key_t *key,
                              const double *payload, size_t payload_count) {
    if(ocrt_surface_pcache_init_from_env()!=0)return -1;
    if(g_pcache_cfg.mode!=OCRT_SURFACE_PCACHE_BUILD)return 0;
    if(!ocrt_pcache_key_valid_(key,payload_count)||!payload)return -1;

    char path[OCRT_PCACHE_PATH_CAP]; ocrt_hash128_t kh,ph;
    if(ocrt_pcache_make_path_(key,path,sizeof path,&kh)!=0)return -1;
    ocrt_pcache_hash_payload_(payload,payload_count,&ph);
    char tmp[OCRT_PCACHE_PATH_CAP];
    /* BUILD is an explicit single-prebuilder phase.  A payload-hash suffix
     * plus a thread-local counter avoids process-ID and OS-runtime APIs while
     * keeping incomplete files separate from the content-addressed target. */
    int n=snprintf(tmp,sizeof tmp,"%s.tmp.%016llx.%lu",path,
                   (unsigned long long)ph.a,++g_pcache_tmp_counter);
    if(n<0||(size_t)n>=sizeof tmp)return -1;
    FILE *fp=fopen(tmp,"wb");
    if(!fp){fprintf(stderr,"error: cannot create persistent surface cache temp file %s: %s\n",tmp,strerror(errno));return -1;}
    int ok=1;
    if(ocrt_write_bytes_(fp,OCRT_PCACHE_MAGIC,OCRT_PCACHE_MAGIC_N)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,OCRT_SURFACE_PCACHE_FORMAT_VERSION)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,OCRT_SURFACE_PCACHE_ABI_VERSION)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,OCRT_SURFACE_PCACHE_PHYSICS_VERSION)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,OCRT_SURFACE_PCACHE_QFOLD_VERSION)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,key->operator_kind)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,(uint32_t)key->mode_first)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,(uint32_t)key->mode_count)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,(uint32_t)key->n_o)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,(uint32_t)key->n_i)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,(uint32_t)key->n_phi_quad)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,(uint32_t)key->sigma_type)!=0)ok=0;
    if(ok&&ocrt_write_u32_(fp,(uint32_t)key->q_convention)!=0)ok=0;
    if(ok&&ocrt_write_f64_(fp,key->wind_speed)!=0)ok=0;
    if(ok&&ocrt_write_f64_(fp,key->n_water)!=0)ok=0;
    if(ok&&ocrt_write_u64_(fp,kh.a)!=0)ok=0;
    if(ok&&ocrt_write_u64_(fp,kh.b)!=0)ok=0;
    if(ok&&ocrt_write_u64_(fp,(uint64_t)payload_count)!=0)ok=0;
    if(ok&&ocrt_write_u64_(fp,ph.a)!=0)ok=0;
    if(ok&&ocrt_write_u64_(fp,ph.b)!=0)ok=0;
    if(ok&&ocrt_write_f64_array_(fp,key->mu_o,(size_t)key->n_o)!=0)ok=0;
    if(ok&&ocrt_write_f64_array_(fp,key->mu_i,(size_t)key->n_i)!=0)ok=0;
    if(ok&&ocrt_write_f64_array_(fp,payload,payload_count)!=0)ok=0;
    if(fflush(fp)!=0)ok=0;
    if(fclose(fp)!=0)ok=0;
    if(!ok){remove(tmp);fprintf(stderr,"error: failed to write persistent surface cache %s\n",tmp);return -1;}

    /* BUILD mode is an explicit single-prebuilder step.  Publication uses
     * the ISO C rename operation; no lock, poll, wait, or OS-specific API is
     * used.  If an identical target already exists, keep it.  REQUIRED mode
     * validates the full header/key/payload, so an interrupted publication is
     * detected rather than silently consumed. */
    FILE *existing=fopen(path,"rb");
    if(existing){fclose(existing);remove(tmp);}
    else if(rename(tmp,path)!=0){
        existing=fopen(path,"rb");
        if(existing){fclose(existing);remove(tmp);}
        else{fprintf(stderr,"error: cannot publish persistent surface cache %s: %s\n",path,strerror(errno));remove(tmp);return -1;}
    }
    if(g_pcache_cfg.trace)
        fprintf(stderr,"[surface-pcache] stored op=%s modes=%d+%d no=%d ni=%d file=%s\n",
                ocrt_pcache_operator_name_(key->operator_kind),key->mode_first,key->mode_count,
                key->n_o,key->n_i,path);
    return 0;
}
