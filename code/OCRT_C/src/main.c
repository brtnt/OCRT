#ifndef _WIN32
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>

#ifdef _WIN32
/* Native Windows/MinGW mkdir compatibility. Environment/time shims are in
 * rt_windows_compat.h so every translation unit uses the same contract. */
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
#endif
#include "rt_windows_compat.h"

#include "rt_types.h"
#include "rt_solver.h"
#include "rt_io.h"
#include "rt_atm.h"
#include "rt_rayleigh.h"
#include "rt_aerosol.h"
#include "rt_aerosol_runtime.h"
#include "rt_absorption.h"
#include "rt_water_iop.h"
#include "rt_iop_ahn_mineral.h"
#include "rt_iop_organic.h"
#include "rt_water_rt.h"
#include "rt_spectral_contract.h"
#include "rt_quadrature.h"
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "shared/mie_io.h"
#include "shared/surface_persistent_cache.h"

/* Scalar FF validation mode IOP table: ASCII/CSV columns wavelength_nm,a,b.
 * Header aliases accepted for a,b to support hand-written input files. */
typedef struct { double wl, a, b; } scalar_ff_iop_row_t;
static int scalar_ff_iop_cmp_main(const void *pa, const void *pb) {
    const scalar_ff_iop_row_t *a = (const scalar_ff_iop_row_t*)pa;
    const scalar_ff_iop_row_t *b = (const scalar_ff_iop_row_t*)pb;
    return (a->wl > b->wl) - (a->wl < b->wl);
}
static int scalar_ff_col_match(const char *name, const char *want) {
    if (!name || !want) return 0;
    char a[128], b[128];
    size_t i;
    for (i = 0; i + 1 < sizeof(a) && name[i]; ++i) {
        char c = name[i];
        if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
        if (c == '_' || c == '-' || c == '/' || c == '(' || c == ')' || c == '[' || c == ']') c = ' ';
        a[i] = c;
    }
    a[i] = 0;
    for (i = 0; i + 1 < sizeof(b) && want[i]; ++i) {
        char c = want[i];
        if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
        if (c == '_' || c == '-' || c == '/' || c == '(' || c == ')' || c == '[' || c == ']') c = ' ';
        b[i] = c;
    }
    b[i] = 0;
    return strstr(a, b) != NULL;
}
static int scalar_ff_iop_lookup_main(const char *path, double wavelength_nm, double *a_out, double *b_out) {
    if (!path || !*path || !a_out || !b_out || !(wavelength_nm > 0.0)) return -1;
    FILE *fp = fopen(path, "r");
    if (!fp) return -2;
    scalar_ff_iop_row_t *rows = NULL;
    int n = 0, cap = 0;
    int iw = 0, ia = 1, ib = 2;
    char line[8192];
    int header_seen = 0;
    while (fgets(line, sizeof line, fp)) {
        char *p0 = line;
        while (*p0 == ' ' || *p0 == '\t') p0++;
        if (*p0 == '#' || *p0 == '\n' || *p0 == '\r' || !*p0) continue;
        char tmp[8192];
        strncpy(tmp, p0, sizeof(tmp)-1); tmp[sizeof(tmp)-1] = 0;
        char *tok[64]; int nt = 0;
        for (char *t = strtok(tmp, ", \t\r\n"); t && nt < 64; t = strtok(NULL, ", \t\r\n")) tok[nt++] = t;
        if (nt < 3) continue;
        if (!header_seen) {
            char *end = NULL;
            (void)strtod(tok[0], &end);
            if (end == tok[0] || *end != 0) {
                for (int k = 0; k < nt; ++k) {
                    if (scalar_ff_col_match(tok[k], "wavelength") || scalar_ff_col_match(tok[k], "lambda")) iw = k;
                    if (!strcmp(tok[k], "a") || scalar_ff_col_match(tok[k], "a m") || scalar_ff_col_match(tok[k], "a total")) ia = k;
                    if (!strcmp(tok[k], "b") || scalar_ff_col_match(tok[k], "b m") || scalar_ff_col_match(tok[k], "b total")) ib = k;
                }
                header_seen = 1;
                continue;
            }
            header_seen = 1;
        }
        if (iw >= nt || ia >= nt || ib >= nt) continue;
        double wl = atof(tok[iw]);
        double aa = atof(tok[ia]);
        double bb = atof(tok[ib]);
        if (!(wl > 0.0) || !(aa >= 0.0) || !(bb > 0.0)) continue;
        if (n >= cap) {
            cap = cap ? cap * 2 : 64;
            scalar_ff_iop_row_t *nr = (scalar_ff_iop_row_t*)realloc(rows, (size_t)cap * sizeof(*rows));
            if (!nr) { fclose(fp); free(rows); return -3; }
            rows = nr;
        }
        rows[n++] = (scalar_ff_iop_row_t){ wl, aa, bb };
    }
    fclose(fp);
    if (n < 1) { free(rows); return -4; }
    qsort(rows, (size_t)n, sizeof(*rows), scalar_ff_iop_cmp_main);
    if (wavelength_nm < rows[0].wl || wavelength_nm > rows[n-1].wl) { free(rows); return -5; }
    int hi = 0;
    while (hi < n && rows[hi].wl < wavelength_nm) hi++;
    if (hi < n && fabs(rows[hi].wl - wavelength_nm) < 1e-12) {
        *a_out = rows[hi].a; *b_out = rows[hi].b; free(rows); return 0;
    }
    if (hi <= 0 || hi >= n) { free(rows); return -6; }
    int lo = hi - 1;
    double f = (wavelength_nm - rows[lo].wl) / (rows[hi].wl - rows[lo].wl);
    *a_out = rows[lo].a * (1.0 - f) + rows[hi].a * f;
    *b_out = rows[lo].b * (1.0 - f) + rows[hi].b * f;
    free(rows);
    return 0;
}

#define V2_VERSION "OCRT-v1.2-2026-08-16-KST-mie-fr631-direct-truncation"

enum {
    WATER_OPTION_FAMILY_OCRT = 1u << 0,
    WATER_OPTION_FAMILY_CCRR = 1u << 1,
    WATER_OPTION_FAMILY_IOP  = 1u << 2
};

static const char *water_input_mode_name(rt_water_input_mode_t mode) {
    switch (mode) {
        case RT_WATER_INPUT_OCRT: return "ocrt";
        case RT_WATER_INPUT_CCRR: return "ccrr";
        case RT_WATER_INPUT_IOP:  return "iop";
        default:                  return "unset";
    }
}

static unsigned water_option_family_for_mode(rt_water_input_mode_t mode) {
    switch (mode) {
        case RT_WATER_INPUT_OCRT: return WATER_OPTION_FAMILY_OCRT;
        case RT_WATER_INPUT_CCRR: return WATER_OPTION_FAMILY_CCRR;
        case RT_WATER_INPUT_IOP:  return WATER_OPTION_FAMILY_IOP;
        default:                  return 0u;
    }
}

static int water_input_mode_parse(const char *value, rt_water_input_mode_t *out) {
    if (!value || !out) return -1;
    if (!strcmp(value, "ocrt")) *out = RT_WATER_INPUT_OCRT;
    else if (!strcmp(value, "ccrr")) *out = RT_WATER_INPUT_CCRR;
    else if (!strcmp(value, "iop")) *out = RT_WATER_INPUT_IOP;
    else return -1;
    return 0;
}

static int ahn_species_parse_strict(const char *value, ahn_species_t *out) {
    if (!value || !out) return -1;
    if (!strcmp(value, "red_clay")) *out = AHN_RED_CLAY;
    else if (!strcmp(value, "brown_earth")) *out = AHN_BROWN_EARTH;
    else if (!strcmp(value, "yellow_clay")) *out = AHN_YELLOW_CLAY;
    else if (!strcmp(value, "calcareous_sand")) *out = AHN_CALCAREOUS_SAND;
    else return -1;
    return 0;
}

static void usage(FILE *fp, const char *prog) {
    /* Split into multiple fprintf calls to avoid ISO C99 4095-char
     * string-literal length limit (previous single call was 4204 chars
     * and produced [-Woverlength-strings] warning). */
    fprintf(fp,
        "OCRT solver " V2_VERSION " (OCRT/CCRR/IOP water interface)\n"
        "\n"
        "Usage (2026-07-14 M2 실행 모드):\n"
        "  [LUT 격자 모드 - 기본]  %s --sza S --wavelength W [opts]\n"
        "     --vza/--raa 를 지정하지 않으면 균일 vza/raa 격자(기본 2.5도 간격,\n"
        "     vza 0..85도)를 한 번의 solve 로 산출한다.  간격은 --lut-vza-step /\n"
        "     --lut-raa-step 로 조정.  출력은 --output-full-grid F.csv (미지정 시\n"
        "     result<N>.csv 자동, csv 외 확장자는 오류).\n"
        "  [단일 기하 모드]  ... --vza V --raa R  (둘 다 필수; 하나만 주면 오류)\n"
        "     관측 방향을 0-가중 노드로 추가(view-as-node 기법, 자동)해 정확\n"
        "     추출한다.\n"
        "  [배치]  %s --batch IN.csv --output OUT.csv [opts]\n"
        "  모드 정책이 충돌하는 플래그 조합(단일 기하 + 격자/배치 등)은 오류 종료.\n"
        "\n"
        "Geometry / spectral:\n"
        "  --sza DEG          solar zenith angle (deg) [필수]\n"
        "  --wavelength NM    wavelength in nm [필수, 330..1100]\n"
        "  --vza DEG          view zenith angle (deg, 단일 기하 모드)\n"
        "  --raa DEG          OCRT relative azimuth (deg): 180=direct-glint branch,\n"
        "                     0=opposite principal-plane branch (단일 기하 모드)\n"
        "  --surface STR      black | flat | black_fresnel_ocean | ocean  [필수, 2026-07-15]\n"
        "  --wind-speed MS    풍속 [m/s], 거친 Fresnel 해면 경사분산 (0 = 평면 프레넬)\n"
        "                       [ocean/black_fresnel_ocean 에서 필수 -- 0 도 명시, 2026-07-15]\n"
        "  --n-water X        해수 굴절률 (기본 1.34 고정 -- 2026-07-15;\n"
        "                       OSOAA/AF1982 정합.  구 Quan-Fry 자동 기본 폐지)\n"
        "  --decouple-sunglint  직달 선글린트 단일반사 항을 표적 출력에서\n"
        "                       분리·제외 (기본은 포함 -- 2026-07-15 반전;\n"
        "                       AF1982 대조 시 지정)\n"
        "\n"
        "구면 보정 (IPSS):\n"
        "  --pssa              구면 보정 on/off. 알고리즘은 IPSS 하나뿐이다\n"
        "                       (Zhai & Hu 2022, JQSRT 282, 108132): 시선 위 각\n"
        "                       점의 실제 태양빔 기하로 단일산란을 구면 계산하고\n"
        "                       다중/단일 비율을 평면평행에서 가져와 곱한다.\n"
        "                       --sza/--vza/--raa 는 모두 지상 화소 기준으로\n"
        "                       해석한다(위성 L1B solz/senz 규약). 평면평행 극한\n"
        "                       에서 보정량이 0 으로 수렴한다.\n"
        "                       적용 범위(Phase I): 대기 단독/흑해양 경로. 결합\n"
        "                       --water-model 실행과 함께 쓰면 fail-loud.\n"
        "                       구 평균할선 Chapman PSSA 는 v1.11 에서 삭제됐다.\n"
        "\n"
        "Aerosol (--aod-555/--aod-865 > 0 이면 활성; 검증 프로토콜 기본 자동):\n"
        "  --mie FILE         .mie 에어로졸 모델 (aod>0 이면 필수)\n"
        "  --aod-555 X        AOD at 555 nm; internally scaled to each wavelength\n  --aod-865 X        AOD at 865 nm; internally scaled to each wavelength\n  --aod X            compatibility alias: AOD at --aod-ref-wavelength (default 555 nm)\n  --aod-ref-wavelength NM  reference wavelength for --aod\n"
        "  --aer-h-km H       OCRT exponential aerosol scale height [km]\n"
        "                       (default 2.0; positive advanced override only)\n"
        "  규칙: --aod-555 0(기본)에서는 에어로졸 옵션 지정이 오류다.\n"
        "        --aod-555/--aod-865 > 0 에는 --mie 가 필수다.\n"
        "  --pressure P       surface pressure [hPa] (default 1013.25). Rayleigh tau scales as\n"
        "                       P/1013.25; --pressure 0 => NO Rayleigh. Atmosphere has no on/off\n"
        "                       toggle now: Rayleigh is set by pressure, aerosol by AOD (--aod,\n"
        "                       default 0 = no aerosol). For an in-water-only test use --pressure 0.\n"
        "\n"
        "Solver options:\n"
        "  --verbose            extra diagnostic output\n",
        prog, prog);

    fprintf(fp,
        "\n"
        "Atmospheric absorption (기본 활성 -- 2026-07-15):\n"
        "  기체흡수(HITRAN 단면적 + AFGL 프로파일 층별 적분)는 항상 켜져 있다.\n"
        "  완전히 끄려면 6기체 총량을 전부 0 으로 지정한다(레일리 전용 결과와\n"
        "  bit 동일):  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0\n"
        "              --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0\n"
        "  --atm-profile STR    AFGL atmosphere profile (default: usstd76).\n"
        "                         usstd76 | tropical | mlsumm | mlwint |\n"
        "                         sasumm  | sawint   | userdef\n"
        "                       userdef → inputs/afgl_atm/afgl_userdef.dat 를\n"
        "                                  읽는다(사용자가 직접 작성; 경로 고정).\n"
        "\n"
        "  Gas column override — 농도만 넣으면 AFGL 기본 vertical profile\n"
        "  shape 유지하면서 total column 만 scaling 됨. 미지정 시 AFGL default.\n"
        "  권장 (변동성 큰 가스):\n"
        "    --gas-column-h2o G   H2O column [g/cm²]   (alias: --pwv-g-cm)\n"
        "    --gas-column-o3 DU   O3  column [DU]      (1 DU = 2.6868e16 mol/cm²)\n"
        "    --gas-column-no2 DU  NO2 column [DU]\n"
        "  Advanced (well-mixed gases, 일반적으로 default 사용):\n"
        "    --gas-column-o2 N    O2  column [mol/cm²] (default ~4.5e24)\n"
        "    --gas-column-co2 N   CO2 column [mol/cm²] (default ~7.1e21 ≈ 400 ppmv)\n"
        "    --gas-column-ch4 N   CH4 column [mol/cm²] (default ~3.6e19)\n");

    fprintf(fp,
        "\n"
        "In-water RT (used with --surface ocean):\n"
        "  --water-temperature C   pure-water temperature [°C] (default 20)\n"
        "  --water-salinity G      pure-water salinity [g/kg] (default 38.4)\n"
        "\n"
        "  Water input branch (required exactly once with --surface ocean):\n"
        "  --water-model M         ocrt | ccrr | iop\n"
        "                          Omitting it, repeating it, or mixing branch prefixes is an error.\n"
        "\n"
        "  OCRT constituent model (all three values required; explicit zero allowed):\n"
        "  --ocrt-chl X            chlorophyll-a [mg m^-3]\n"
        "  --ocrt-tsm X            inorganic TSM dry mass [g m^-3]\n"
        "  --ocrt-adom440 X        aDOM/CDOM absorption at 440 nm [m^-1]\n"
        "  --ocrt-adom-slope X     aDOM exponential slope [nm^-1] (optional; default 0.014)\n"
        "  --ocrt-phyto-group G    EAP species selector [CURRENTLY DISABLED when Chl>0]\n"
        "                          Legacy names: pico | nano | micro.  The 17-species\n"
        "                          coated-sphere catalog is integrated for offline/API\n"
        "                          generation, but constituent-model EAP scattering is\n"
        "                          disabled until a validated L=200 forward-peak treatment\n"
        "                          is implemented.  Omit this option: Chl then affects\n"
        "                          absorption only and detritus supplies particle scattering.\n"
        "  --ocrt-tsm-species N    red_clay(default)|brown_earth|yellow_clay|calcareous_sand\n"
        "  --ocrt-detritus-a440 X  organic-detritus a(440) [m^-1] (default 0)\n"
        "  --ocrt-detritus-slope X detritus slope [nm^-1] (default 0.0109)\n"
        "\n"
        "  CCRR constituent adapter (all three values required; explicit zero allowed):\n"
        "  --ccrr-chl X            chlorophyll-a [mg m^-3]\n"
        "  --ccrr-tsm X            TSM concentration [g m^-3]\n"
        "  --ccrr-adom440 X        aDOM absorption at 440 nm [m^-1]\n"
        "  --ccrr-adom-slope X     aDOM exponential slope [nm^-1] (optional; default 0.014)\n"
        "                          CCRR converts these inputs to OCRT RT inputs; the RT solver\n"
        "                          itself remains OCRT. Positive Chl uses the bundled Morel\n"
        "                          (1988)/MM01 classic Case-1 absorption closure.\n"
        "\n"
        "  Pure water rule:\n"
        "    In OCRT or CCRR mode, Chl=0, TSM=0 and aDOM440=0 lower to the same\n"
        "    native pure-water path. A model still must be selected explicitly.\n"
        "\n"
        "  Direct total-IOP model (all three values required):\n"
        "  --iop-a A               total absorption a [m^-1], A >= 0\n"
        "  --iop-b B               total scattering b [m^-1], B > 0\n"
        "  --iop-bb BB             total backscattering bb [m^-1], 0 < BB < 0.5 B\n"
        "  --iop-phase-lut F       optional scalar P11 phase LUT\n"
        "  --iop-mie-phase F       optional vector P11/P12/P33 .mie phase\n"
        "\n"
        "  Old --simple-*, unprefixed --chl/--tsm/--adom*, native --cdom-*, and\n"
        "  --fixed-bulk-iop names are rejected to prevent branch ambiguity.\n"
        "  Advanced model-specific phase controls use --ocrt-mie-* or --iop-mie-* prefixes.\n"
        "  --ccrr-phase-moments F  optional CCRR particle-moment override.\n");

    fprintf(fp,
        "  IOP advanced/validation controls:\n"
        "                          HydroLight fixed-total-IOP validation mode: inject\n"
        "                          a_total,b_total,bb_total directly and use scalar\n"
        "                          bulk phase matched to bb_total/b_total.\n"
        "  --iop-mie-moment-mode M  dense|gauss for injected .mie phases.\n"
        "  --ocrt-mie-moment-mode M dense|gauss for packaged OCRT phases.\n"
        "  --iop-mie-moment-nmu N / --ocrt-mie-moment-nmu N (default 400).\n"
        "  --iop-phase-lut F          CSV with theta_deg,P11; optional case_id,wavelength_nm;\n"
        "                              production path is weak-cap lut-deltam.\n"
        "  --iop-mie-truncation / --ocrt-mie-truncation\n"
        "                          broad-transform the selected particle .mie phase.\n"
        "                          --ocrt-* is ADVANCED (OCRT_ADVANCED=1). Canonical\n"
        "                          FR631 scientific reference is truncation OFF; the\n"
        "                          option is for explicitly matched compatibility/\n"
        "                          acceleration tests and changes b, omega, tau and\n"
        "                          the source function. It is never auto-enabled.\n"
        "  --iop-mie-ss-mode N / --ocrt-mie-ss-mode N\n"
        "                          with truncation: 0=none(default), 1=IMS, 2=NT-TMS.\n"
        "  --water-shared-grid     water quadrature shares the atmosphere node set\n"
        "                          (OSOAA-parity single mu grid; excludes --n-mu-water)\n"
        "  --batch-full-grid F     run one ocean full-grid per CSV row, in parallel\n"
        "                          (OpenMP; OMP_NUM_THREADS). Columns: out(req),\n"
        "                          wavelength,sza,wind,aod,aer_h_km,iop_a,iop_b,iop_bb,\n"
        "                          iop_mie_phase,mie,n_mu_water,raa_step,water_model,\n"
        "                          ocrt_chl,ocrt_tsm,ocrt_adom440,ocrt_adom_slope,\n"
        "                          ocrt_phyto_group,ocrt_tsm_species,\n"
        "                          ocrt_detritus_a440,ocrt_detritus_slope,\n"
        "                          ccrr_chl,ccrr_tsm,ccrr_adom440,ccrr_adom_slope;\n"
        "                          empty=CLI base. Branch prefixes may not be mixed.\n"
        "  --water-view-as-node    add requested water-view direction as a zero-weight\n"
        "                          output node when outside the GL range. Default OFF;\n"
        "                          CCRR v0.5.0 uses continuous reconstruction.\n"
        "  --output-full-grid F    with --surface ocean: write water Rrs(0+) I/Q/U\n"
        "                          and rrs(0-) on positive GL view nodes (--n-mu)\n"
        "                          and uniform physical RAA grid (--lut-raa-step).\n"
        "                          In-water quadrature defaults to n_mu_water=96 in\n"
        "                          this mode (interactive default 48 is +0.44%%\n"
        "                          under-converged at high omega; converged ~80,\n"
        "                          96 = margin for forward-peaked phyto phases).\n"
        "\n"
        "Other:\n"
        "  --examples         검증된 주요 실행 예시 6종을 설명과 함께 출력\n"
        "  --output-advanced  --batch 결과 CSV 에 상세 열(기준 대비 오차, 수렴\n"
        "                     차수, 계산 시간)을 추가한다.  기본은 간단 열\n"
        "                     (I/Q/U 반사도만).  구 --output-mode 는 삭제됨.\n"
        "  --help             this help\n"
        "  --version          print version and exit\n"
        "\n"
        "Gates (이 도움말은 일반 옵션만 싣는다 -- Jae 지시 2026-07-15):\n"
        "  --ocrt-mie-truncation과 --n-mu-water는 명시적으로 사용할 때 항상\n"
        "  OCRT_ADVANCED=1이 필요하다. 기타 수치 수렴 노브(--n-mu --n-layers\n"
        "  --max-orders --sos-* --conv-tol --l-max --m-max --water-m-max\n"
        "  --iop-phase-nphi --sigma-model --sigma-type --q-convention)는 기본값이\n"
        "  검증 구성이며, 기본값과 다른 값 지정에 OCRT_ADVANCED=1이 필요하다.\n"
        "  진단 옵션(--debug-* --scalar-ff-* --ccrr-particle-phase-*\n"
        "  --integration-method constant 등)은\n"
        "  OCRT_DEBUG=1 필요이며 역시 제외했다.\n"
        "  게이트 옵션 전체 목록과 설명은 docs/OCRT_USAGE 및 main.c 상단\n"
        "  게이트 일람 주석(4절)에 있다.\n"
        "  실행 예시는  %s --examples  로 표시된다.\n", prog);
}

/* --examples: 검증된 주요 실행 예시 출력.
 * 아래 6개 명령은 OCRT v1.1 정식 빌드에서 실행 검증한다(오류 없이 실행,
 * 소요 시간은 1코어 실측치).  예시를 바꾸면 반드시 재실행으로 검증하고
 * docs/OCRT_USAGE 를 같은 커밋에서 갱신한다(유지보수 규칙, 상단 주석 참조).
 * 명령은 ocrt/ 디렉터리에서 실행해야 한다(inputs/ 상대경로 로드). */
static void usage_examples(FILE *fp, const char *prog) {
    fprintf(fp,
        "OCRT 실행 예시 6종 (v1.1 water-branch contract + CCRR Chl table)\n"
        "공통: ocrt/ 디렉터리에서 실행한다. 재현성 기준은 OMP_NUM_THREADS=1.\n"
        "\n"
        "[1] 해양 단일 기하 -- 한 관측 방향의 원격반사도 (약 4초)\n"
        "  %s --surface ocean --water-model ocrt \\\n"
        "    --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0.1 \\\n"
        "    --ocrt-adom-slope 0.014 --wind-speed 3 \\\n"
        "    --sza 40 --vza 30 --raa 90 --wavelength 443 --pressure 1013.25\n"
        "  CDOM(용존유기물)만 있는 해수 + Rayleigh 대기 + 기체흡수(기본 활성).\n"
        "  OCRT 구성성분 입력은 Chl/TSM/aDOM440 세 값을 모두 명시한다.\n"
        "  aDOM slope 는 선택값이며 생략 시 기본 0.014 nm^-1 이다.\n"
        "  해수 굴절률 1.34 고정, 선글린트 포함이 기본이다(2026-07-15;\n"
        "  분리는 --decouple-sunglint).\n"
        "  관측 방향 (vza=30,\n"
        "  raa=90)을 0-가중 노드로 추가해(view-as-node 자동) 보간 없이 추출한다.\n"
        "  stdout 1행에서 rrs0minus(수면 직하 원격반사도 [1/sr]) 와\n"
        "  Rrs0plus(수면 직상) 를 읽는다.\n"
        "\n"
        "[2] 해양 격자(LUT) 모드 -- 전 각도 장을 한 번에 (15도/30도 격자 72행 약 2분)\n"
        "  %s --surface ocean --water-model ocrt \\\n"
        "    --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0.1 \\\n"
        "    --ocrt-adom-slope 0.014 --wind-speed 3 \\\n"
        "    --sza 40 --wavelength 443 --pressure 1013.25 \\\n"
        "    --lut-vza-step 15 --lut-raa-step 30 --output-full-grid rrs_grid.csv\n"
        "  --vza/--raa 를 빼면 격자 모드다.  vza 0..85도 x raa 균일 격자를\n"
        "  단일 solve 로 산출한다(S7 재생).  출력 파일 미지정 시 result<N>.csv\n"
        "  자동 명명.  수중 절점은 자동으로 96 승격(정확도 마진).\n"
        "  기본 간격 2.5도의 전체 격자는 이보다 훨씬 무겁다.\n",
        prog, prog);
    fprintf(fp,
        "\n"
        "[3] 대기 전용 TOA 반사도 -- Rayleigh + 거친 해면 (약 1초)\n"
        "  %s --surface black_fresnel_ocean --wind-speed 3 \\\n"
        "    --sza 40 --vza 30 --raa 90 --wavelength 555 --pressure 1013.25\n"
        "  대기보정 항 산출의 기본형: Rayleigh + 기체흡수(usstd76 기본 컬럼,\n"
        "  2026-07-15 부터 기본 활성).  stdout 첫 3값 = rho_I rho_Q rho_U\n"
        "  (반사도, F_sun=pi 규약).  --pressure 0 = Rayleigh 없음(기체흡수는\n"
        "  독립 -- 끄는 법은 예시 [5]).\n"
        "\n"
        "[4] 에어로졸 추가 (약 4초)\n"
        "  %s --surface black_fresnel_ocean --wind-speed 3 \\\n"
        "    --sza 40 --vza 30 --raa 90 --wavelength 555 --pressure 1013.25 \\\n"
        "    --mie inputs/M80C.mie --aod-555 0.3\n"
        "  --aod-555/--aod-865 > 0 이면 검증 프로토콜 구성(loglin 절단 + value 커널 + L=80 +\n"
        "  m=16 + 400층 + scale height 2.0 km)이 자동 기본이다 -- 최소 커맨드가\n"
        "  곧 프로토콜(2026-07-15 부터 결합 케이스 포함).  --aod-555/--aod-865 > 0 에는 --mie\n"
        "  가 필수, --aod-555 0(기본)에서 에어로졸 옵션 지정은 오류다.\n"
        "  공식 에어로졸 모델: inputs/ 의 T50, C50, M80C 등 .mie 파일.\n",
        prog, prog);
    fprintf(fp,
        "\n"
        "[5] 기체 흡수 제어 -- 총량 오버라이드와 완전 off (각 약 1초)\n"
        "  기체흡수는 기본 활성이다(HITRAN 단면적 + AFGL usstd76 층별 적분).\n"
        "  (a) 총량 오버라이드(형상 유지, 총량만 스케일; O3/NO2 는 DU,\n"
        "      H2O 는 g/cm2):\n"
        "  %s ... --pressure 1013.25 --gas-column-o3 200\n"
        "  (b) 완전 off = 6기체 전부 0 (레일리 전용 결과와 bit 동일):\n"
        "  %s ... --pressure 1013.25 \\\n"
        "    --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \\\n"
        "    --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0\n"
        "\n"
        "[6] 배치 격자 -- CSV 한 행마다 해양 전체 격자 1개 (2행 예시 약 30초)\n"
        "  printf 'out,wavelength,sza\\n"
        "b443.csv,443,40\\n"
        "b555.csv,555,40\\n"
        "' > batch.csv\n"
        "  %s --surface ocean --water-model ocrt \\\n"
        "    --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0.1 \\\n"
        "    --ocrt-adom-slope 0.014 --wind-speed 3 \\\n"
        "    --sza 40 --wavelength 443 --pressure 1013.25 \\\n"
        "    --lut-vza-step 30 --lut-raa-step 90 --batch-full-grid batch.csv\n"
        "  행 병렬 실행(OpenMP, OMP_NUM_THREADS 로 조절).  CSV 의 빈 칸은\n"
        "  커맨드라인 값을 기본으로 쓴다.  --sza/--wavelength 는 배치에서도\n"
        "  커맨드라인에 필수다(행 미지정 시의 기본값 역할).\n"
        "\n"
        "상세 옵션 설명은 --help 와 docs/OCRT_USAGE 를 본다.\n",
        prog, prog, prog);
}

/* parse a required double argument; exit on missing value */
static double next_double(int argc, char **argv, int *i, const char *flag) {
    if (*i + 1 >= argc) {
        fprintf(stderr, "error: %s requires a value\n", flag);
        exit(2);
    }
    *i += 1;
    return strtod(argv[*i], NULL);
}

static int next_int(int argc, char **argv, int *i, const char *flag) {
    if (*i + 1 >= argc) {
        fprintf(stderr, "error: %s requires a value\n", flag);
        exit(2);
    }
    *i += 1;
    return (int)strtol(argv[*i], NULL, 10);
}

static const char *next_str(int argc, char **argv, int *i, const char *flag) {
    if (*i + 1 >= argc) {
        fprintf(stderr, "error: %s requires a value\n", flag);
        exit(2);
    }
    *i += 1;
    return argv[*i];
}


#define CCRR_CSV_MAX_FIELDS 128
#define CCRR_CSV_MAX_LINE   8192

static int split_csv_simple(char *line, char **fields, int max_fields) {
    int n = 0;
    char *p = line;
    while (p && n < max_fields) {
        fields[n++] = p;
        char *comma = strchr(p, ',');
        if (!comma) break;
        *comma = '\0';
        p = comma + 1;
    }
    for (int i = 0; i < n; ++i) {
        char *q = fields[i];
        while (*q == ' ' || *q == '\t' || *q == '\r' || *q == '\n') ++q;
        fields[i] = q;
        size_t L = strlen(q);
        while (L > 0 && (q[L-1] == ' ' || q[L-1] == '\t' || q[L-1] == '\r' || q[L-1] == '\n')) {
            q[--L] = '\0';
        }
    }
    return n;
}

static int header_index(char **fields, int n, const char *name) {
    for (int i = 0; i < n; ++i) {
        if (!strcmp(fields[i], name)) return i;
    }
    return -1;
}

static double field_double(char **fields, int n, int idx, double fallback) {
    if (idx < 0 || idx >= n || !fields[idx] || !*fields[idx]) return fallback;
    return strtod(fields[idx], NULL);
}

static void warn_value_kernel_exact_nadir_once(void)
{
    static int warned = 0;
    const char *e = getenv("OCRT_WATER_VALUE_KERNEL_POL");
    if (warned || !e || !*e || !strcmp(e, "0")) return;
    warned = 1;
    fprintf(stderr,
        "warning: OCRT_WATER_VALUE_KERNEL_POL at exact VZA=0 deg has a known "
        "polarized output-extraction singularity; use VZA=0.001 deg "
        "(validated safe range 0.001--0.01 deg; avoid <1e-6 deg). "
        "The solver input is not modified automatically.\n");
}

static double normalize_raa_0_360(double raa_deg) {
    while (raa_deg < 0.0) raa_deg += 360.0;
    while (raa_deg >= 360.0) raa_deg -= 360.0;
    return raa_deg;
}

static double ccrr_raa_to_ocrt_raa(double raa_ccrr_deg) {
    /* CCRR v0.5.3 uses a propagation-azimuth convention.
     * OCRT public RAA convention is not changed.
     * Therefore the mapping is applied only while reading CCRR reference rows. */
    return normalize_raa_0_360(180.0 - raa_ccrr_deg);
}


/* -------------------------------------------------------------------------
 * Above-water Rrs angular-grid output for coupled ocean cases.
 *
 * The native coupled LUT path computes the shared atmosphere pass-1, water SOS
 * field, and atmosphere pass-2 once per case/band.  It extracts the authoritative
 * water target Fourier samples once per requested VZA and reconstructs all RAAs
 * without re-entering the coupled solver.  The legacy per-cell replay remains
 * available only for regression comparison via OCRT_NATIVE_COUPLED_LUT_OFF=1.
 * ------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------
 * v1.09 commit #19: OpenMP-parallel batch of independent ocean full-grid runs.
 * Each CSV row describes one full-grid run (one output file per row); rows are
 * embarrassingly parallel and executed with `#pragma omp parallel for`.
 * Schema (header required; unknown columns rejected; empty cell = inherit CLI):
 *   out (required), wavelength, sza, wind, aod, aer_h_km,
 *   iop_a, iop_b, iop_bb (all three or none), iop_mie_phase, mie,
 *   n_mu_water (requires OCRT_ADVANCED=1), raa_step,
 *   water_model,
 *   ocrt_chl, ocrt_tsm, ocrt_adom440, ocrt_adom_slope,
 *   ocrt_phyto_group, ocrt_tsm_species, ocrt_detritus_a440,
 *   ocrt_detritus_slope,
 *   ccrr_chl, ccrr_tsm, ccrr_adom440, ccrr_adom_slope
 * Thread safety: per-row rt_case_t/rt_options_t copies; threadprivate grid
 * cache + per-run rt_water_grid_cache_reset(); OCRT_WATER_GRID_CACHE set once
 * before the parallel region (no setenv inside threads); one output file per
 * row (no contention). */
static int run_ocean_rrs_full_grid_csv(const char *out_csv,
                                       const rt_case_t *base_case,
                                       const rt_options_t *base_opts,
                                       const char *water_aw_lut_path,
                                       const char *water_psi_T_lut_path,
                                       const char *mie_path,
                                       double user_aod,
                                       double user_aod_ref_nm,
                                       const rt_aerosol_runtime_options_t *aer_ropts,
                                       double vza_step_deg,
                                       double vza_max_deg,
                                       double raa_step_deg);

typedef struct {
    char out[512], wphase[512], mie[512];
    char water_model[16], phyto_group[32], tsm_species[64];
    double wl, sza, wind, aod, aer_h, a, b, bb, raa_step;
    double chl, tsm, adom440, adom_slope, det_a440, det_slope;
    int has_wl, has_sza, has_wind, has_aod, has_aerh, has_iop, has_wp, has_mie,
        has_nmw, has_rs, nmw;
    int has_water_model, has_chl, has_tsm, has_adom440, has_adom_slope,
        has_phyto_group, has_tsm_species, has_det_a440, has_det_slope;
    unsigned water_family_mask;
} fg_row_t;

static int fg_col(const char *hdr, const char *name) { return strcmp(hdr, name) == 0; }

static int fg_copy_field(char *dst, size_t cap, const char *src,
                         const char *column, int row_no) {
    if (!dst || cap == 0 || !src) return -1;
    const size_t n = strlen(src);
    if (n >= cap) {
        fprintf(stderr,
                "error: --batch-full-grid row %d column '%s' is too long (max %zu characters)\n",
                row_no, column ? column : "?", cap - 1);
        return -1;
    }
    memcpy(dst, src, n + 1);
    return 0;
}

static int run_ocean_full_grid_batch(const char *list_csv,
                                     const rt_case_t *base_case,
                                     const rt_options_t *base_opts,
                                     const char *water_aw_lut_path,
                                     const char *water_psi_T_lut_path,
                                     const char *base_mie_path,
                                     double base_user_aod,
                                     double base_user_aod_ref_nm,
                                     const rt_aerosol_runtime_options_t *base_aer_ropts,
                                     double vza_step,
                                     double vza_max,
                                     double base_raa_step,
                                     int n_water_explicit,
                                     int wind_explicit /* 2026-07-15 G4 */) {
    /* DOC-REF: docs/PERSISTENT_COXMUNK_INTERFACE_CACHE_KO_2026-08-25.md
     * Persistent surface-cache BUILD is a separate priming step.  Batch
     * simulation rows may only read immutable cache files (required/readonly)
     * or run with persistence off.  No row is allowed to build shared state or
     * wait for another row to publish it. */
    if (ocrt_surface_pcache_init_from_env() != 0 ||
        ocrt_surface_pcache_batch_guard() != 0) {
        return 2;
    }
    FILE *fp = fopen(list_csv, "r");
    if (!fp) { fprintf(stderr, "error: --batch-full-grid cannot open '%s'\n", list_csv); return 1; }
    char line[4096];
    if (!fgets(line, sizeof line, fp)) { fclose(fp); return 1; }
    /* header */
    enum { MAXC = 32 };
    char cols[MAXC][64]; int ncol = 0;
    for (char *tok = strtok(line, ",\r\n"); tok && ncol < MAXC; tok = strtok(NULL, ",\r\n"))
        { snprintf(cols[ncol], sizeof cols[ncol], "%s", tok); ++ncol; }
    int cap = 256, nrow = 0;
    fg_row_t *rows = calloc((size_t)cap, sizeof *rows);
    if (!rows) { fclose(fp); return 1; }
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == '\n' || line[0] == '#') continue;
        if (nrow == cap) { cap *= 2; rows = realloc(rows, (size_t)cap * sizeof *rows);
                           memset(rows + nrow, 0, (size_t)(cap - nrow) * sizeof *rows); }
        fg_row_t *r = &rows[nrow];
        /* v1.10 FIX: hand-rolled splitter that PRESERVES empty fields.
         * strtok_r collapses consecutive delimiters, which shifted every
         * value after an empty column one header to the LEFT (e.g. the
         * n_mu_water value landing in the mie column).  Exposed by the
         * fixed-bulk + full-grid combination (empty phase columns). */
        int c = 0; char *p = line;
        for (; c < ncol; ++c) {
            char *v = p;
            while (*p && *p != ',' && *p != '\r' && *p != '\n') ++p;
            if (*p) *p++ = '\0';
            while (*p == '\r' || *p == '\n') { *p = '\0'; }
            if (!*v) continue;
            if      (fg_col(cols[c], "out"))            { if (fg_copy_field(r->out, sizeof r->out, v, "out", nrow + 1) != 0) { fclose(fp); return 2; } }
            else if (fg_col(cols[c], "wavelength"))     { r->wl = atof(v); r->has_wl = 1; }
            else if (fg_col(cols[c], "sza"))            { r->sza = atof(v); r->has_sza = 1; }
            else if (fg_col(cols[c], "wind"))           { r->wind = atof(v); r->has_wind = 1; }
            else if (fg_col(cols[c], "aod"))            { r->aod = atof(v); r->has_aod = 1; }
            else if (fg_col(cols[c], "aer_h_km"))       { r->aer_h = atof(v); r->has_aerh = 1; }
            else if (fg_col(cols[c], "iop_a"))          { r->a = atof(v); r->has_iop |= 1; r->water_family_mask |= WATER_OPTION_FAMILY_IOP; }
            else if (fg_col(cols[c], "iop_b"))          { r->b = atof(v); r->has_iop |= 2; r->water_family_mask |= WATER_OPTION_FAMILY_IOP; }
            else if (fg_col(cols[c], "iop_bb"))         { r->bb = atof(v); r->has_iop |= 4; r->water_family_mask |= WATER_OPTION_FAMILY_IOP; }
            else if (fg_col(cols[c], "iop_mie_phase"))  { if (fg_copy_field(r->wphase, sizeof r->wphase, v, "iop_mie_phase", nrow + 1) != 0) { fclose(fp); return 2; } r->has_wp = 1; r->water_family_mask |= WATER_OPTION_FAMILY_IOP; }
            else if (fg_col(cols[c], "mie"))            { if (fg_copy_field(r->mie, sizeof r->mie, v, "mie", nrow + 1) != 0) { fclose(fp); return 2; } r->has_mie = 1; }
            else if (fg_col(cols[c], "n_mu_water"))     { r->nmw = atoi(v); r->has_nmw = 1; }
            else if (fg_col(cols[c], "raa_step"))       { r->raa_step = atof(v); r->has_rs = 1; }
            else if (fg_col(cols[c], "water_model"))    { if (fg_copy_field(r->water_model, sizeof r->water_model, v, "water_model", nrow + 1) != 0) { fclose(fp); return 2; } r->has_water_model = 1; }
            else if (fg_col(cols[c], "ocrt_chl"))       { r->chl = atof(v); r->has_chl = 1; r->water_family_mask |= WATER_OPTION_FAMILY_OCRT; }
            else if (fg_col(cols[c], "ocrt_tsm"))       { r->tsm = atof(v); r->has_tsm = 1; r->water_family_mask |= WATER_OPTION_FAMILY_OCRT; }
            else if (fg_col(cols[c], "ocrt_adom440"))   { r->adom440 = atof(v); r->has_adom440 = 1; r->water_family_mask |= WATER_OPTION_FAMILY_OCRT; }
            else if (fg_col(cols[c], "ocrt_adom_slope")){ r->adom_slope = atof(v); r->has_adom_slope = 1; r->water_family_mask |= WATER_OPTION_FAMILY_OCRT; }
            else if (fg_col(cols[c], "ocrt_phyto_group")){ if (fg_copy_field(r->phyto_group, sizeof r->phyto_group, v, "ocrt_phyto_group", nrow + 1) != 0) { fclose(fp); return 2; } r->has_phyto_group = 1; r->water_family_mask |= WATER_OPTION_FAMILY_OCRT; }
            else if (fg_col(cols[c], "ocrt_tsm_species")){ if (fg_copy_field(r->tsm_species, sizeof r->tsm_species, v, "ocrt_tsm_species", nrow + 1) != 0) { fclose(fp); return 2; } r->has_tsm_species = 1; r->water_family_mask |= WATER_OPTION_FAMILY_OCRT; }
            else if (fg_col(cols[c], "ocrt_detritus_a440")){ r->det_a440 = atof(v); r->has_det_a440 = 1; r->water_family_mask |= WATER_OPTION_FAMILY_OCRT; }
            else if (fg_col(cols[c], "ocrt_detritus_slope")){ r->det_slope = atof(v); r->has_det_slope = 1; r->water_family_mask |= WATER_OPTION_FAMILY_OCRT; }
            else if (fg_col(cols[c], "ccrr_chl"))       { r->chl = atof(v); r->has_chl = 1; r->water_family_mask |= WATER_OPTION_FAMILY_CCRR; }
            else if (fg_col(cols[c], "ccrr_tsm"))       { r->tsm = atof(v); r->has_tsm = 1; r->water_family_mask |= WATER_OPTION_FAMILY_CCRR; }
            else if (fg_col(cols[c], "ccrr_adom440"))   { r->adom440 = atof(v); r->has_adom440 = 1; r->water_family_mask |= WATER_OPTION_FAMILY_CCRR; }
            else if (fg_col(cols[c], "ccrr_adom_slope")){ r->adom_slope = atof(v); r->has_adom_slope = 1; r->water_family_mask |= WATER_OPTION_FAMILY_CCRR; }
            else { fprintf(stderr, "error: --batch-full-grid unknown column '%s'\n", cols[c]);
                   fclose(fp); return 2; }
        }
        if (!r->out[0]) { fprintf(stderr, "error: --batch-full-grid row %d missing 'out'\n", nrow + 1);
                          fclose(fp); return 2; }
        if (r->has_iop && r->has_iop != 7) {
            fprintf(stderr, "error: --batch-full-grid row %d needs all of iop_a,iop_b,iop_bb\n", nrow + 1);
            fclose(fp); return 2;
        }
        {
            rt_water_input_mode_t row_mode = (rt_water_input_mode_t)base_case->water_input_mode;
            if (r->has_water_model && water_input_mode_parse(r->water_model, &row_mode) != 0) {
                fprintf(stderr,
                        "error: --batch-full-grid row %d water_model must be ocrt, ccrr or iop\n",
                        nrow + 1);
                fclose(fp); return 2;
            }
            const unsigned expected = water_option_family_for_mode(row_mode);
            if ((r->water_family_mask & ~expected) != 0u) {
                fprintf(stderr,
                        "error: --batch-full-grid row %d mixes option prefixes incompatible with water_model=%s\n",
                        nrow + 1, water_input_mode_name(row_mode));
                fclose(fp); return 2;
            }
            if (r->has_water_model && row_mode != (rt_water_input_mode_t)base_case->water_input_mode) {
                if ((row_mode == RT_WATER_INPUT_OCRT || row_mode == RT_WATER_INPUT_CCRR) &&
                    (!r->has_chl || !r->has_tsm || !r->has_adom440)) {
                    fprintf(stderr,
                            "error: --batch-full-grid row %d changes water_model to %s and must supply all three matching Chl/TSM/aDOM440 columns\n",
                            nrow + 1, water_input_mode_name(row_mode));
                    fclose(fp); return 2;
                }
                if (row_mode == RT_WATER_INPUT_IOP && r->has_iop != 7) {
                    fprintf(stderr,
                            "error: --batch-full-grid row %d changes water_model to iop and must supply iop_a,iop_b,iop_bb\n",
                            nrow + 1);
                    fclose(fp); return 2;
                }
            }
        }
        if (r->has_nmw) {
            const char *adv = getenv("OCRT_ADVANCED");
            if (!adv || !*adv || !strcmp(adv, "0")) {
                fprintf(stderr, "error: --batch-full-grid n_mu_water column requires OCRT_ADVANCED=1\n");
                fclose(fp); return 2;
            }
        }
        ++nrow;
    }
    fclose(fp);
    if (nrow == 0) { fprintf(stderr, "error: --batch-full-grid '%s' has no rows\n", list_csv); return 2; }

    /* Hoist the grid-cache gate out of the threads (setenv is not thread-safe). */
    if (!getenv("OCRT_WATER_GRID_CACHE")) setenv("OCRT_WATER_GRID_CACHE", "1", 1);
    rt_water_env_init();   /* #21: snapshot env once, before threads spawn */

    fprintf(stderr, "# batch-full-grid: %d rows, OMP threads=%d\n", nrow,
#ifdef _OPENMP
            omp_get_max_threads()
#else
            1
#endif
            );
    int fails = 0;
    const int _bdiag = (getenv("OCRT_BATCH_DIAG") && getenv("OCRT_BATCH_DIAG")[0]) ? 1 : 0;
    int    *_row_tid = _bdiag ? (int*)calloc((size_t)nrow, sizeof(int)) : NULL;
    double *_row_sec = _bdiag ? (double*)calloc((size_t)nrow, sizeof(double)) : NULL;
    double *_row_start = _bdiag ? (double*)calloc((size_t)nrow, sizeof(double)) : NULL;
    double *_row_end   = _bdiag ? (double*)calloc((size_t)nrow, sizeof(double)) : NULL;
    double _wall_t0 = 0.0;
    if (_bdiag) {
        struct timespec _ts; clock_gettime(CLOCK_MONOTONIC, &_ts);
        _wall_t0 = _ts.tv_sec + _ts.tv_nsec*1e-9;
        int _nth = 1;
#ifdef _OPENMP
        #pragma omp parallel
        { 
            #pragma omp single
            _nth = omp_get_num_threads();
        }
#endif
        fprintf(stderr, "# BATCH_DIAG: nrow=%d, OMP threads=%d\n", nrow, _nth);
    }
    #pragma omp parallel for schedule(dynamic, 1) reduction(+: fails)
    for (int i = 0; i < nrow; ++i) {
        double _r_t0 = 0.0;
        if (_bdiag) { struct timespec _ts; clock_gettime(CLOCK_MONOTONIC,&_ts); _r_t0=_ts.tv_sec+_ts.tv_nsec*1e-9; _row_start[i]=_r_t0-_wall_t0; }
        rt_case_t    cs   = *base_case;
        rt_options_t opts = *base_opts;
        fg_row_t    *r    = &rows[i];
        if (r->has_wl) {
            cs.wavelength_nm = r->wl;
            /* v1.09 #19 fix: re-derive wavelength-dependent fields exactly as
             * the single-case preparation does (main ~2055).  n_water is the
             * only cs field derived from wavelength; an explicit --n-water on
             * the base CLI stays fixed across rows. */
            if (!n_water_explicit) cs.n_water = rt_water_iop_n_real(cs.wavelength_nm);
        }
        if (r->has_sza)  cs.sza_deg       = r->sza;
        if (r->has_wind) cs.wind_speed    = r->wind;
        /* 2026-07-15 Jae(G4): 행 단위 wind 필수(CLI 에도 행에도 없으면 SKIP). */
        if (!wind_explicit && !r->has_wind) {
            fprintf(stderr, "# batch-full-grid row %d/%d SKIP out=%s "
                    "(error: wind required -- row or CLI --wind-speed)\n",
                    i + 1, nrow, r->out);
            fails += 1;
            continue;
        }
        rt_water_input_mode_t row_mode = (rt_water_input_mode_t)cs.water_input_mode;
        if (r->has_water_model) {
            /* Validated while reading the CSV. */
            (void)water_input_mode_parse(r->water_model, &row_mode);
        }

        if (row_mode != (rt_water_input_mode_t)base_case->water_input_mode) {
            /* A row that changes branch starts from neutral water-property
             * state; it may not inherit values or phase sources from another
             * branch.  The CSV reader already required a complete replacement
             * field set for this case. */
            cs.ccrr_mode = 0;
            cs.fixed_bulk_iop_mode = 0;
            cs.ccrr_chl_mg_m3 = 0.0;
            cs.ccrr_min_g_m3 = 0.0;
            cs.a_cdom_440_m_inv = 0.0;
            cs.S_cdom_nm_inv = 0.014;
            cs.cdom_ref_lambda_nm = 440.0;
            cs.detritus_a440_m_inv = 0.0;
            cs.detritus_slope_nm_inv = ORGANIC_DETRITUS_SLOPE_DEFAULT;
            cs.ccrr_phase_moments_path = NULL;
            cs.ccrr_particle_phase_lut_path = NULL;
            cs.ccrr_particle_phase_case_id = NULL;
            cs.fixed_a_total_m_inv = 0.0;
            cs.fixed_b_total_m_inv = 0.0;
            cs.fixed_bb_total_m_inv = 0.0;
            cs.fixed_bulk_phase_lut_path = NULL;
            cs.water_mie_phase_path = NULL;
        }
        cs.water_input_mode = (int)row_mode;

        if (row_mode == RT_WATER_INPUT_OCRT || row_mode == RT_WATER_INPUT_CCRR) {
            cs.fixed_bulk_iop_mode = 0;
            cs.water_constituent_model = (row_mode == RT_WATER_INPUT_OCRT)
                                           ? RT_WATER_CONSTITUENT_OCRT
                                           : RT_WATER_CONSTITUENT_CCRR;
            if (r->has_chl) cs.ccrr_chl_mg_m3 = r->chl;
            if (r->has_tsm) cs.ccrr_min_g_m3 = r->tsm;
            if (r->has_adom440) cs.a_cdom_440_m_inv = r->adom440;
            if (r->has_adom_slope) cs.S_cdom_nm_inv = r->adom_slope;
            if (r->has_det_a440) cs.detritus_a440_m_inv = r->det_a440;
            if (r->has_det_slope) cs.detritus_slope_nm_inv = r->det_slope;

            if (r->has_phyto_group) {
                organic_phyto_group_t group;
                if (rt_iop_organic_group_parse(r->phyto_group, &group) != 0) {
                    fails += 1;
                    #pragma omp critical
                    fprintf(stderr,
                            "# batch-full-grid row %d/%d rc=2 out=%s (error: ocrt_phyto_group must be pico, nano or micro)\n",
                            i + 1, nrow, r->out);
                    continue;
                }
                cs.organic_phyto_group = (int)group;
            }
            if (r->has_tsm_species) {
                ahn_species_t species = AHN_RED_CLAY;
                if (ahn_species_parse_strict(r->tsm_species, &species) != 0) {
                    fails += 1;
                    #pragma omp critical
                    fprintf(stderr,
                            "# batch-full-grid row %d/%d rc=2 out=%s (error: invalid ocrt_tsm_species)\n",
                            i + 1, nrow, r->out);
                    continue;
                }
                cs.tsm_species = (int)species;
            }

            if (cs.ccrr_chl_mg_m3 < 0.0 || cs.ccrr_min_g_m3 < 0.0 ||
                cs.a_cdom_440_m_inv < 0.0 || cs.detritus_a440_m_inv < 0.0 ||
                !(cs.S_cdom_nm_inv > 0.0) || !(cs.detritus_slope_nm_inv > 0.0)) {
                fails += 1;
                #pragma omp critical
                fprintf(stderr,
                        "# batch-full-grid row %d/%d rc=2 out=%s (error: invalid constituent value)\n",
                        i + 1, nrow, r->out);
                continue;
            }

            const int row_pure =
                (cs.ccrr_chl_mg_m3 == 0.0 &&
                 cs.ccrr_min_g_m3 == 0.0 &&
                 cs.a_cdom_440_m_inv == 0.0);
            cs.ccrr_mode = row_pure ? 0 : 1;
            cs.F_sun = row_pure ? M_PI : 1.0;

            if (row_mode == RT_WATER_INPUT_OCRT) {
                if (r->has_tsm_species && !(cs.ccrr_min_g_m3 > 0.0)) {
                    fails += 1;
                    #pragma omp critical
                    fprintf(stderr,
                            "# batch-full-grid row %d/%d rc=2 out=%s (error: ocrt_tsm_species requires ocrt_tsm>0)\n",
                            i + 1, nrow, r->out);
                    continue;
                }
                if ((r->has_phyto_group || r->has_det_a440 || r->has_det_slope) &&
                    !(cs.ccrr_chl_mg_m3 > 0.0)) {
                    fails += 1;
                    #pragma omp critical
                    fprintf(stderr,
                            "# batch-full-grid row %d/%d rc=2 out=%s (error: OCRT organic options require ocrt_chl>0)\n",
                            i + 1, nrow, r->out);
                    continue;
                }
                if (cs.ccrr_chl_mg_m3 > 0.0 &&
                    !rt_iop_organic_wavelength_supported(cs.wavelength_nm)) {
                    fails += 1;
                    #pragma omp critical
                    fprintf(stderr,
                            "# batch-full-grid row %d/%d rc=2 out=%s (error: OCRT Chl wavelength outside 350-850 nm)\n",
                            i + 1, nrow, r->out);
                    continue;
                }
                if (row_pure) {
                    cs.water_mie_phase_path = NULL;
                } else if (cs.ccrr_min_g_m3 > 0.0 && cs.ccrr_chl_mg_m3 <= 0.0) {
                    cs.water_mie_phase_path = rt_iop_ahn_mineral_phase_path(
                        (ahn_species_t)cs.tsm_species);
                } else {
                    cs.water_mie_phase_path = NULL; /* organic or mixed phase is automatic */
                }
            } else {
                cs.water_mie_phase_path = NULL;
            }
        } else if (row_mode == RT_WATER_INPUT_IOP) {
            cs.ccrr_mode = 0;
            cs.fixed_bulk_iop_mode = 1;
            if (r->has_iop == 7) {
                cs.fixed_a_total_m_inv  = r->a;
                cs.fixed_b_total_m_inv  = r->b;
                cs.fixed_bb_total_m_inv = r->bb;
            }
            if (!(cs.fixed_a_total_m_inv >= 0.0) ||
                !(cs.fixed_b_total_m_inv > 0.0) ||
                !(cs.fixed_bb_total_m_inv > 0.0) ||
                cs.fixed_bb_total_m_inv >= 0.5 * cs.fixed_b_total_m_inv) {
                fails += 1;
                #pragma omp critical
                fprintf(stderr,
                        "# batch-full-grid row %d/%d rc=2 out=%s (error: invalid IOP values)\n",
                        i + 1, nrow, r->out);
                continue;
            }
            cs.F_sun = 1.0;
            if (r->has_wp) cs.water_mie_phase_path = r->wphase;
        } else {
            fails += 1;
            #pragma omp critical
            fprintf(stderr,
                    "# batch-full-grid row %d/%d rc=2 out=%s (error: water_model unset)\n",
                    i + 1, nrow, r->out);
            continue;
        }
        if (r->has_nmw && cs.water_shared_grid) {
            fails += 1;
            #pragma omp critical
            fprintf(stderr, "# batch-full-grid row %d/%d rc=2 out=%s "
                            "(error: n_mu_water column conflicts with --water-shared-grid)\n",
                    i + 1, nrow, r->out);
            continue;
        }
        if (r->has_nmw)  cs.n_mu_water_override  = r->nmw;
        if (r->has_mie)  opts.aerosol_mie_path   = r->mie;
        if (r->has_aod)  { opts.aerosol_aod = r->aod; cs.aerosol_on = (r->aod > 0.0) ? 1 : 0; }
        /* 2026-07-15 Jae(G3): 행 단위 가드 -- 병합 후 aod>0 인데 mie 가
         * 행에도 CLI 에도 없으면 그 행은 실패 처리한다. */
        if (opts.aerosol_aod > 0.0 && !opts.aerosol_mie_path) {
            fprintf(stderr, "# batch-full-grid row %d/%d SKIP out=%s "
                    "(error: aod>0 requires mie -- row or CLI)\n",
                    i + 1, nrow, r->out);
            fails += 1;
            continue;
        }
        if (r->has_aerh) {
            if (!(r->aer_h > 0.0)) {
                fails += 1;
                #pragma omp critical
                fprintf(stderr, "# batch-full-grid row %d/%d rc=2 out=%s (error: aer_h_km must be positive)\n",
                        i + 1, nrow, r->out);
                continue;
            }
            opts.aer_h_km = r->aer_h;
        }
        double rs = r->has_rs ? r->raa_step : base_raa_step;
        const char *row_mie = r->has_mie ? r->mie : base_mie_path;
        const double row_aod = r->has_aod ? r->aod : base_user_aod;
        int rc = run_ocean_rrs_full_grid_csv(r->out, &cs, &opts,
                                             water_aw_lut_path, water_psi_T_lut_path,
                                             row_mie, row_aod, base_user_aod_ref_nm,
                                             base_aer_ropts,
                                             vza_step, vza_max, rs);
        if (rc != 0) fails += 1;
        if (_bdiag) {
            struct timespec _ts; clock_gettime(CLOCK_MONOTONIC,&_ts);
            double _now=_ts.tv_sec+_ts.tv_nsec*1e-9;
            _row_sec[i]=_now-_r_t0;
            _row_end[i]=_now-_wall_t0;
#ifdef _OPENMP
            _row_tid[i]=omp_get_thread_num();
#endif
        }
        #pragma omp critical
        fprintf(stderr, "# batch-full-grid row %d/%d rc=%d out=%s\n", i + 1, nrow, rc, r->out);
    }
    if (_bdiag) {
        struct timespec _ts; clock_gettime(CLOCK_MONOTONIC,&_ts);
        double _wall = (_ts.tv_sec+_ts.tv_nsec*1e-9) - _wall_t0;
        int _maxtid=0; for(int i=0;i<nrow;i++) if(_row_tid[i]>_maxtid)_maxtid=_row_tid[i];
        double _busy_sum=0.0;
        for(int t=0;t<=_maxtid;t++){
            int _cnt=0; double _sec=0.0;
            for(int i=0;i<nrow;i++) if(_row_tid[i]==t){_cnt++;_sec+=_row_sec[i];}
            _busy_sum+=_sec;
            fprintf(stderr,"# BATCH_DIAG thread %d: rows=%d busy=%.3fs\n",t,_cnt,_sec);
        }
        int _NS=1000; double _conc_sum=0.0; int _conc_max=0;
        for(int sIdx=0;sIdx<_NS;sIdx++){
            double _tt=_wall*((double)sIdx+0.5)/(double)_NS;
            int _inflight=0;
            for(int i=0;i<nrow;i++) if(_row_start[i]<=_tt && _tt<_row_end[i]) _inflight++;
            _conc_sum+=_inflight; if(_inflight>_conc_max)_conc_max=_inflight;
        }
        double _mean_conc=_conc_sum/(double)_NS;
        fprintf(stderr,"# BATCH_DIAG wall=%.3fs sum_busy=%.3fs parallel_efficiency=%.2f\n",
                _wall,_busy_sum,_busy_sum/(_wall*(_maxtid+1)+1e-9));
        fprintf(stderr,"# BATCH_DIAG mean_concurrency=%.2f threads_spawned=%d peak_concurrency=%d\n",
                _mean_conc,_maxtid+1,_conc_max);
        fprintf(stderr,"# BATCH_DIAG VERDICT: %s\n",
                (_mean_conc >= (_maxtid+1)*0.7) ? "PARALLEL_OK threads overlap as expected" :
                (_maxtid+1 <= 1) ? "ONLY_1_THREAD check build -fopenmp / OMP_NUM_THREADS / affinity" :
                "SERIALIZED threads spawned but run one-at-a-time -> internal lock/contention");
        fprintf(stderr,"# BATCH_DIAG timeline (row: tid start..end s):\n");
        for(int i=0;i<nrow && i<24;i++)
            fprintf(stderr,"#   row %2d: tid %d  %.3f..%.3f  (%.3fs)\n",
                    i+1,_row_tid[i],_row_start[i],_row_end[i],_row_sec[i]);
        free(_row_tid); free(_row_sec); free(_row_start); free(_row_end);

    }
    free(rows);
    fprintf(stderr, "# batch-full-grid done: %d rows, %d failed\n", nrow, fails);
    return fails ? 1 : 0;
}

static int run_ocean_rrs_full_grid_csv(const char *out_csv,
                                       const rt_case_t *base_case,
                                       const rt_options_t *base_opts,
                                       const char *water_aw_lut_path,
                                       const char *water_psi_T_lut_path,
                                       const char *mie_path,
                                       double user_aod,
                                       double user_aod_ref_nm,
                                       const rt_aerosol_runtime_options_t *aer_ropts,
                                       double vza_step_deg,
                                       double vza_max_deg,
                                       double raa_step_deg) {
    if (!base_case || !base_opts || !water_aw_lut_path) return -1;
    if (base_case->surface != RT_SURFACE_OCEAN) {
        fprintf(stderr, "error: --output-full-grid currently requires --surface ocean\n");
        return 2;
    }
    rt_water_grid_cache_reset();   /* worker-private water field: one independent case */
    rt_solver_all_view_cache_reset(); /* worker-private S7/S7b: no preceding-row state */
    if (!(raa_step_deg > 0.0 && raa_step_deg <= 180.0)) {
        fprintf(stderr, "error: --lut-raa-step/RAA grid step must be in (0,180]\n");
        return 2;
    }

    rt_water_iop_lut_t aw_lut = {0};
    rt_water_iop_psi_T_lut_t psi_T_lut = {0};
    if (rt_water_iop_lut_load(water_aw_lut_path, &aw_lut) != 0) {
        fprintf(stderr, "error: failed to load water a/b LUT '%s'\n", water_aw_lut_path);
        return 1;
    }
    int have_psi_T = (rt_water_iop_psi_T_load(water_psi_T_lut_path, &psi_T_lut) == 0);

    /* Full-grid must use the same aerosol runtime object as the single-view
     * coupled-ocean path. Before this fix the caller passed NULL here, so
     * --mie/--aod were silently ignored in ocean full-grid output. */
    rt_aerosol_input_t aer_fg = {0};
    rt_aerosol_input_t *aer_fg_ptr = NULL;
    const mie_data_t *mie_fg = NULL;
    mie_data_t mie_fg_tmp = {0};
    int mie_fg_tmp_loaded = 0;
    if (base_case->aerosol_on && mie_path && *mie_path && user_aod > 0.0) {
        int model_cache_hit = 0;
        const int mrc = mie_model_cache_get(mie_path, &mie_fg, &model_cache_hit);
        (void)model_cache_hit;
        if (mrc == 1) {
            if (read_mie_file(mie_path, &mie_fg_tmp) != 0) mie_fg = NULL;
            else { mie_fg = &mie_fg_tmp; mie_fg_tmp_loaded = 1; }
        }
        if (mrc < 0 || !mie_fg) {
            fprintf(stderr, "error: full-grid Mie model load('%s') failed\n", mie_path);
            rt_water_iop_lut_free(&aw_lut);
            if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
            return 1;
        }
        rt_aerosol_runtime_options_t ro = aer_ropts ? *aer_ropts : (rt_aerosol_runtime_options_t){
            .aer_L_max = base_opts->aerosol_l_max,
            .theta_cut_deg = base_opts->aerosol_theta_cut_deg,
            .delta_m_N = base_opts->aerosol_delta_m_N,
            .apply_nt_tau = base_opts->aerosol_apply_nt_tau,
            .use_loglin_trunc = 0,
            .loglin_mu1 = 0.8, .loglin_mu2 = 0.94, .loglin_threshold = 0.1
        };
        rt_aerosol_runtime_diag_t ad = {0};
        if (rt_aerosol_runtime_prepare(mie_fg, base_case->wavelength_nm,
                                       user_aod, user_aod_ref_nm, &ro,
                                       &aer_fg, &ad) != 0) {
            fprintf(stderr, "error: full-grid aerosol runtime prepare failed\n");
            if (mie_fg_tmp_loaded) mie_data_free(&mie_fg_tmp);
            rt_water_iop_lut_free(&aw_lut);
            if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
            return 1;
        }
        aer_fg_ptr = &aer_fg;
        fprintf(stderr, "[ocean+aer full-grid] AOD_ref=%.4f@%.1fnm AOD_band=%.4f tau_eff=%.4f ssa=%.4f L_max=%d\n",
                ad.aod_ref, ad.aod_ref_nm, ad.aod_target, aer_fg.tau_a, aer_fg.ssa_a, aer_fg.L_max);
    }

    /* 2026-07-14 M2 (Jae 지시): 해양-대기 통일 균일 격자.  출력 vza 는
     * 0..vza_max 를 vza_step 간격(기본 2.5도)으로 균일 배치한다.  기존
     * GL-노드 vza 규약은 폐기.  각 vza 는 아래에서 0-가중 view 노드(#20)로
     * 수중 solve 에 전달되어 노드 정확 추출되며, 대기부는 S7 캐시가 LUT
     * 재구성(노드 일치 bit-exact, 비노드 PCHIP)으로 셀 간 재사용한다. */
    if (!(vza_step_deg > 0.0 && vza_step_deg <= 90.0)) {
        fprintf(stderr, "error: --lut-vza-step must be in (0,90]\n");
        rt_water_iop_lut_free(&aw_lut);
        if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
        return 2;
    }
    if (!(vza_max_deg > 0.0 && vza_max_deg < 90.0)) vza_max_deg = 85.0;
    int n_vza = (int)floor(vza_max_deg / vza_step_deg + 1e-9) + 1;
    if (n_vza < 1) n_vza = 1;
    if (n_vza > 512) {
        fprintf(stderr, "error: vza grid too fine (%d values > 512); increase --lut-vza-step\n", n_vza);
        rt_water_iop_lut_free(&aw_lut);
        if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
        return 2;
    }
    if (getenv("OCRT_WATER_VALUE_KERNEL_POL"))
        warn_value_kernel_exact_nadir_once();  /* LUT grid includes VZA=0 */
    int n_mu = base_opts->n_mu > 0 ? base_opts->n_mu : 24;   /* 수중/대기 구적 크기(출력 격자와 무관) */

    int n_raa = (int)floor(360.0 / raa_step_deg + 0.5);
    if (n_raa < 1) n_raa = 1;
    FILE *fp = out_csv ? fopen(out_csv, "w") : stdout;
    if (!fp) {
        fprintf(stderr, "error: cannot open full-grid output '%s'\n", out_csv);
        rt_water_iop_lut_free(&aw_lut);
        if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
        return 1;
    }

    fprintf(fp,
        "case_kind,sza_deg,wavelength_nm,wind_speed_ms,n_mu,mu_index,mu_view,vza_deg,raa_deg,"
        "Ed0plus_air,Ed0plus_direct,Ed0plus_diffuse,"
        "Lu0plus_I,Lu0plus_Q,Lu0plus_U,Rrs_I,Rrs_Q,Rrs_U,"
        "Ed0minus_water,Ed0minus_direct,Ed0minus_diffuse,"
        "Eu0minus_water,Eu0minus_direct,Eu0minus_diffuse,"
        "Lu0minus_I,Lu0minus_Q,Lu0minus_U,rrs_I,rrs_Q,rrs_U,"
        "TOA_rho_I,TOA_rho_Q,TOA_rho_U,"
        "T_dir_dn,T_diff_dn_hemi,T_total_dn_hemi,T_diff_dn_dir,"
        "T_dir_up_view,T_diff_up_view,T_total_up_view,TOA_water_signal_I,T_up_rt_valid,"
        "a_total,b_total,bb_total,omega_water,water_orders,water_converged,"
        "a_w,b_w,bb_w,a_chl,"
        "a_phyto_detritus,b_phyto_detritus,bb_phyto_detritus,"
        "a_dom,a_min,b_min,bb_min,Kd0minus\n");

    /* Native coupled angular LUT mode.
     * One atmosphere/water cache-fill solve builds the shared Fourier fields.
     * The authoritative water target projection is evaluated once per VZA;
     * all RAAs in that row are reconstructed without re-entering the coupled
     * solver.  S7/S7b provide the all-angle atmosphere fields and S17a/S17
     * retain exact near-nadir samples.  Opt-out for regression comparison:
     * OCRT_NATIVE_COUPLED_LUT_OFF=1. */
    if (!getenv("OCRT_WATER_GRID_CACHE")) setenv("OCRT_WATER_GRID_CACHE", "1", 1);
    if (!getenv("OCRT_ATM_GRID_CACHE")) setenv("OCRT_ATM_GRID_CACHE", "1", 1);
    { /* v1.10 S7: publish the grid raa list for the atm cache */
        char raabuf[4096]; int off = 0;   /* M2: raa 2.5도 격자(144개) 수용 */
        for (int ir2 = 0; ir2 < n_raa && off < (int)sizeof raabuf - 16; ++ir2)
            off += snprintf(raabuf + off, sizeof raabuf - (size_t)off,
                            "%s%.9g", ir2 ? "," : "", (double)ir2 * raa_step_deg);
        setenv("OCRT_S7_RAA_LIST", raabuf, 1);
    }
    rt_water_env_init();   /* #21 */
    /* v1.09 #20 (S-007 root fix): pass ALL grid vza values as zero-weight
     * in-water view nodes, so every (vza,raa) cell reads its radiance
     * directly from the solved field (OSOAA UserAngFile-equivalent) instead
     * of interpolating between water Gauss nodes.
     * Opt-out (legacy interpolation, regression comparison): OCRT_GRID_VIEW_INTERP=1 */
    double vza_list_fg[512];
    int use_view_nodes = !getenv("OCRT_GRID_VIEW_INTERP") && n_vza <= 512;
    for (int iv = 0; iv < n_vza; ++iv)
        vza_list_fg[iv] = (double)iv * vza_step_deg;   /* M2: 균일 vza 격자 */
    int n_fail = 0;
    rt_result_t *native_results = (rt_result_t*)calloc((size_t)n_vza * (size_t)n_raa,
                                                        sizeof(rt_result_t));
    rt_ocean_lut_grid_out_t native_grid = {0};
    int native_ok = 0;
    const char *allow_replay_env = getenv("OCRT_ALLOW_LEGACY_CELL_REPLAY");
    const int allow_legacy_replay = allow_replay_env && allow_replay_env[0] &&
                                    strcmp(allow_replay_env, "0") != 0;
    int native_rc = -999;
    if (native_results && !getenv("OCRT_NATIVE_COUPLED_LUT_OFF")) {
        double raa_list_fg[512];
        for (int ir=0; ir<n_raa; ++ir) raa_list_fg[ir]=(double)ir*raa_step_deg;
        native_grid.n_vza=n_vza; native_grid.vza_deg=vza_list_fg;
        native_grid.n_raa=n_raa; native_grid.raa_deg=raa_list_fg;
        native_grid.results=native_results;
        native_rc=rt_solve_case_ocean_lut(base_case, base_opts, &aw_lut,
                                          have_psi_T ? &psi_T_lut : NULL,
                                          aer_fg_ptr, &native_grid);
        if (native_rc==0) {
            native_ok=1;
            fprintf(stderr,"[native-coupled-lut] calls=%d water_cold=%d water_views=%d exact_near_nadir=%d cells=%d replay=0\n",
                    native_grid.coupled_calls,native_grid.water_cold_solves,
                    native_grid.water_view_calls,native_grid.near_nadir_exact_rows,n_vza*n_raa);
        }
    }
    if (!native_ok && !allow_legacy_replay) {
        fprintf(stderr,
                "error: all-view coupled LUT is required (rc=%d).  Per-cell VZA/RAA "
                "solver replay is prohibited by the OCRT execution contract.\n"
                "  Diagnostic-only override: OCRT_ALLOW_LEGACY_CELL_REPLAY=1\n",
                native_rc);
        if (fp != stdout) fclose(fp);
        free(native_results);
        if (aer_fg_ptr) rt_aerosol_runtime_free(&aer_fg);
        if (mie_fg_tmp_loaded) mie_data_free(&mie_fg_tmp);
        rt_water_iop_lut_free(&aw_lut);
        if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
        return 3;
    }
    if (!native_ok) {
        fprintf(stderr,
                "warning: diagnostic legacy cell replay enabled rc=%d cells=%d\n",
                native_rc, n_vza*n_raa);
    }
    for (int iv = 0; iv < n_vza; ++iv) {   /* near-nadir -> grazing */
        double vza = vza_list_fg[iv];
        double m = cos(vza * M_PI / 180.0);
        for (int ir = 0; ir < n_raa; ++ir) {
            double raa = (double)ir * raa_step_deg;
            rt_case_t cs = *base_case;
            cs.vza_deg = vza;
            if (use_view_nodes) {                 /* #20: 균일 격자 전체를 0-가중 노드로 */
                cs.water_view_vza_list = vza_list_fg;
                cs.n_water_view_vza    = n_vza;
            }
            cs.raa_deg = raa;
            cs.water_view_as_node = 0;   /* commit #16: grid mode (see above) */
            rt_result_t r = {0};
            int rc = 0;
            if (native_ok) r = native_results[(size_t)iv*(size_t)n_raa + (size_t)ir];
            else rc = rt_solve_case_ocean(&cs, base_opts, &aw_lut,
                                           have_psi_T ? &psi_T_lut : NULL,
                                           aer_fg_ptr, &r);
            if (rc != 0) {
                ++n_fail;
                fprintf(stderr, "warning: full-grid ocean solve failed rc=%d vza=%.6f raa=%.6f\n", rc, vza, raa);
                continue;
            }
            const double Rq = r.R_rs_0plus_Q;
            const double Ru = r.R_rs_0plus_U;
            fprintf(fp,
                "ocean_rrs_grid,%.10g,%.10g,%.10g,%d,%d,%.12g,%.10g,%.10g,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,%.12e,%d,"
                "%.12e,%.12e,%.12e,%.12e,%d,%d,"
                "%.12e,%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,%.12e,%.12e\n",
                cs.sza_deg, cs.wavelength_nm, cs.wind_speed, n_mu, iv, m, vza, raa,
                r.Ed_0plus_air, r.Ed_0plus_direct, r.Ed_0plus_diffuse,
                r.Lu_0plus_view, r.Qu_0plus_view, r.Uu_0plus_view,
                r.R_rs_0plus, Rq, Ru,
                r.Ed_0minus_water, r.Ed_0minus_direct, r.Ed_0minus_diffuse,
                r.Eu_0minus_water, r.Eu_0minus_direct, r.Eu_0minus_diffuse,
                r.Lu_0minus_view, r.Qu_0minus_view, r.Uu_0minus_view,
                r.r_rs_0minus, r.r_rs_0minus_Q, r.r_rs_0minus_U,
                r.rho_I, r.rho_Q, r.rho_U,
                r.T_dir_dn, r.T_diff_dn_hemi, r.T_total_dn_hemi, r.T_diff_dn_dir,
                r.T_dir_up_view, r.T_diff_up_view, r.T_total_up_view,
                r.I_TOA_water_signal, r.T_up_rt_valid,
                r.a_total_used, r.b_total_used, r.bb_total_used, r.omega_water,
                r.n_orders_used, r.converged,
                r.a_w_used, r.b_w_used, r.bb_w_used, r.a_chl_used,
                r.a_pig_used, r.b_pig_used, r.bb_pig_used,
                r.a_cdom_used, r.a_min_used, r.b_min_used, r.bb_min_used,
                r.Kd_0minus);
        }
    }
    if (fp != stdout) fclose(fp);
    free(native_results);

    if (aer_fg_ptr) rt_aerosol_runtime_free(&aer_fg);
    if (mie_fg_tmp_loaded) mie_data_free(&mie_fg_tmp);
    rt_water_iop_lut_free(&aw_lut);
    if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
    fprintf(stderr, "ocean_rrs_full_grid: n_mu=%d n_raa=%d rows=%d fail=%d replay=%d output=%s\n",
            n_vza, n_raa, n_vza*n_raa - n_fail, n_fail, native_ok ? 0 : n_vza*n_raa,
            out_csv ? out_csv : "stdout");
    return n_fail ? 1 : 0;
}

static int run_ccrr_compare_csv(const char *in_csv, const char *out_csv,
                                const char *water_aw_lut_path,
                                const char *water_psi_T_lut_path,
                                const rt_options_t *base_opts,
                                const rt_case_t *base_case) {
    if (!in_csv || !out_csv || !base_opts || !base_case) return -1;

    FILE *fin = fopen(in_csv, "r");
    if (!fin) {
        fprintf(stderr, "error: cannot open CCRR reference CSV '%s'\n", in_csv);
        return 1;
    }
    FILE *fout = fopen(out_csv, "w");
    if (!fout) {
        fprintf(stderr, "error: cannot create CCRR output CSV '%s'\n", out_csv);
        fclose(fin);
        return 1;
    }

    rt_water_iop_lut_t aw_lut = {0};
    rt_water_iop_psi_T_lut_t psi_T_lut = {0};
    if (rt_water_iop_lut_load(water_aw_lut_path, &aw_lut) != 0) {
        fprintf(stderr, "error: failed to load water LUT '%s'\n", water_aw_lut_path);
        fclose(fin); fclose(fout);
        return 1;
    }
    int have_psi_T = (rt_water_iop_psi_T_load(water_psi_T_lut_path, &psi_T_lut) == 0);

    char line[CCRR_CSV_MAX_LINE];
    if (!fgets(line, sizeof line, fin)) {
        fprintf(stderr, "error: CCRR CSV is empty\n");
        rt_water_iop_lut_free(&aw_lut);
        if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
        fclose(fin); fclose(fout);
        return 1;
    }
    char *hfields[CCRR_CSV_MAX_FIELDS];
    int hn = split_csv_simple(line, hfields, CCRR_CSV_MAX_FIELDS);
    int icase = header_index(hfields, hn, "case_id");
    int iwl   = header_index(hfields, hn, "wavelength_nm");
    int ichl  = header_index(hfields, hn, "chl_mg_m3");
    int imin  = header_index(hfields, hn, "min_g_m3");
    int iadom = header_index(hfields, hn, "adom_440_m-1");
    int isdom = header_index(hfields, hn, "adom_s_nm-1");
    int isza  = header_index(hfields, hn, "sza_air_deg");
    int ivza  = header_index(hfields, hn, "vza_air_deg");
    int iraa  = header_index(hfields, hn, "raa_deg");
    int iref_Edm = header_index(hfields, hn, "Ed0minus");
    int iref_Lum = header_index(hfields, hn, "Lu0minus");
    int iref_rrs = header_index(hfields, hn, "rrs");
    int iref_Edp = header_index(hfields, hn, "Ed0plus");
    int iref_Lup = header_index(hfields, hn, "Lu0plus_no_sky");
    int iref_Rrs = header_index(hfields, hn, "Rrs");
    int iref_a   = header_index(hfields, hn, "a_total_m-1");
    int iref_b   = header_index(hfields, hn, "b_total_m-1");
    int iref_bb  = header_index(hfields, hn, "bb_total_m-1");
    if (icase < 0 || iwl < 0 || ichl < 0 || imin < 0 || iadom < 0 || isdom < 0 ||
        isza < 0 || ivza < 0 || iraa < 0) {
        fprintf(stderr, "error: CCRR CSV missing required columns\n");
        rt_water_iop_lut_free(&aw_lut);
        if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
        fclose(fin); fclose(fout);
        return 1;
    }

    fprintf(fout,
        "case_id,wavelength_nm,chl_mg_m3,min_g_m3,adom_440_m-1,adom_s_nm-1,"
        "sza_air_deg,vza_air_deg,raa_ccrr_deg,raa_ocrt_deg,raa_mapping,wind_speed_ms,sunglint_decoupled,"
        "OCRT_Ed0minus,OCRT_Lu0minus,OCRT_rrs,OCRT_Ed0plus,OCRT_Lu0plus_no_sky,OCRT_Rrs,"
        "OCRT_a_total,OCRT_b_total,OCRT_bb_total,"
        "OCRT_a_w,OCRT_b_w,OCRT_bb_w,OCRT_a_dom,"
        "OCRT_a_pig,OCRT_b_pig,OCRT_bb_pig,OCRT_a_min,OCRT_b_min,OCRT_bb_min,"
        "OCRT_orders,OCRT_conv,"
        "ref_Ed0minus,ref_Lu0minus,ref_rrs,ref_Ed0plus,ref_Lu0plus_no_sky,ref_Rrs,"
        "ref_a_total,ref_b_total,ref_bb_total,"
        "rel_Ed0minus_pct,rel_Lu0minus_pct,rel_rrs_pct,rel_Ed0plus_pct,rel_Lu0plus_no_sky_pct,rel_Rrs_pct,"
        "rel_a_total_pct,rel_b_total_pct,rel_bb_total_pct\n");

    int n_rows = 0, n_fail = 0;
    while (fgets(line, sizeof line, fin)) {
        if (!line[0] || line[0] == '\n' || line[0] == '\r') continue;
        char *fields[CCRR_CSV_MAX_FIELDS];
        int nf = split_csv_simple(line, fields, CCRR_CSV_MAX_FIELDS);
        if (nf <= 1) continue;
        rt_case_t cs = *base_case;
        rt_options_t opts = *base_opts;
        cs.surface = RT_SURFACE_OCEAN;
        cs.ccrr_mode = 1;
        cs.F_sun = 1.0;
        cs.cdom_ref_lambda_nm = 440.0;
        cs.wavelength_nm = field_double(fields, nf, iwl, 0.0);
        cs.ccrr_chl_mg_m3 = field_double(fields, nf, ichl, 0.0);
        cs.ccrr_min_g_m3 = field_double(fields, nf, imin, 0.0);
        cs.water_mie_phase_path = (cs.ccrr_min_g_m3 > 0.0)
            ? rt_iop_ahn_mineral_phase_path((ahn_species_t)cs.tsm_species) : NULL;
        cs.a_cdom_440_m_inv = field_double(fields, nf, iadom, 0.0);
        cs.S_cdom_nm_inv = field_double(fields, nf, isdom, 0.014);
        cs.sza_deg = field_double(fields, nf, isza, 0.0);
        cs.vza_deg = field_double(fields, nf, ivza, 0.0);
        double raa_ccrr_deg = field_double(fields, nf, iraa, 0.0);
        double raa_ocrt_deg = ccrr_raa_to_ocrt_raa(raa_ccrr_deg);
        cs.raa_deg = raa_ocrt_deg;
        cs.water_view_as_node = 1;
        cs.wind_speed = 0.0;
        cs.decouple_sunglint = 1;

        rt_result_t res = {0};
        int rc = rt_solve_case_ocean(&cs, &opts, &aw_lut,
                                      (have_psi_T ? &psi_T_lut : NULL),
                                      NULL, &res);
        if (rc != 0) n_fail++;
        double ref_Edm = field_double(fields, nf, iref_Edm, NAN);
        double ref_Lum = field_double(fields, nf, iref_Lum, NAN);
        double ref_rrs = field_double(fields, nf, iref_rrs, NAN);
        double ref_Edp = field_double(fields, nf, iref_Edp, NAN);
        double ref_Lup = field_double(fields, nf, iref_Lup, NAN);
        double ref_Rrs = field_double(fields, nf, iref_Rrs, NAN);
        double ref_a   = field_double(fields, nf, iref_a, NAN);
        double ref_b   = field_double(fields, nf, iref_b, NAN);
        double ref_bb  = field_double(fields, nf, iref_bb, NAN);
        #define RPCT(x,y) ((isfinite(y) && fabs(y) > 0.0) ? (100.0*((x)/(y)-1.0)) : NAN)
        fprintf(fout,
            "%s,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%s,%.10g,%d,"
            "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,"
            "%.10e,%.10e,%.10e,"
            "%.10e,%.10e,%.10e,%.10e,"
            "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,"
            "%d,%d,"
            "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,"
            "%.10e,%.10e,%.10e,"
            "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,"
            "%.10e,%.10e,%.10e\n",
            (icase >= 0 && icase < nf) ? fields[icase] : "",
            cs.wavelength_nm, cs.ccrr_chl_mg_m3, cs.ccrr_min_g_m3,
            cs.a_cdom_440_m_inv, cs.S_cdom_nm_inv, cs.sza_deg, cs.vza_deg,
            raa_ccrr_deg, raa_ocrt_deg, "180-minus-ccrr", cs.wind_speed, cs.decouple_sunglint,
            res.Ed_0minus_water, res.Lu_0minus_view, res.r_rs_0minus,
            res.Ed_0plus_air, res.Lu_0plus_view, res.R_rs_0plus,
            res.a_total_used, res.b_total_used, res.bb_total_used,
            res.a_w_used, res.b_w_used, res.bb_w_used, res.a_cdom_used,
            res.a_pig_used, res.b_pig_used, res.bb_pig_used,
            res.a_min_used, res.b_min_used, res.bb_min_used,
            res.n_orders_used, res.converged,
            ref_Edm, ref_Lum, ref_rrs, ref_Edp, ref_Lup, ref_Rrs,
            ref_a, ref_b, ref_bb,
            RPCT(res.Ed_0minus_water, ref_Edm), RPCT(res.Lu_0minus_view, ref_Lum),
            RPCT(res.r_rs_0minus, ref_rrs), RPCT(res.Ed_0plus_air, ref_Edp),
            RPCT(res.Lu_0plus_view, ref_Lup), RPCT(res.R_rs_0plus, ref_Rrs),
            RPCT(res.a_total_used, ref_a), RPCT(res.b_total_used, ref_b), RPCT(res.bb_total_used, ref_bb));
        #undef RPCT
        n_rows++;
    }
    rt_water_iop_lut_free(&aw_lut);
    if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
    fclose(fin); fclose(fout);
    fprintf(stderr, "ccrr-compare: rows=%d fail=%d output=%s\n", n_rows, n_fail, out_csv);
    return n_fail ? 1 : 0;
}

/* =============================================================================
 * EMBEDDED USAGE MANUAL (v1.09, 2026-07-05) — 재발 방지용 코드 내장판
 * =============================================================================
 * 마이그레이션/패키징에서 문서가 유실되어도 이 파일 하나로 사용법이 전달되도록
 * main() 위에 내장한다.  전체판: docs/OCRT_USAGE_v1.09.md,  옵션 1차 소스: --help.
 * 결함/이력: BUG_SUSPECTS.md + SESSION_START_VERIFICATION_GUIDELINE.
 *
 * ---------------------------------------------------------------------------
 * 0. 빌드 & 실행 위치 (가장 흔한 실패 원인)
 * ---------------------------------------------------------------------------
 *   gcc -std=c11 -O2 -fopenmp -Isrc $(find src -name '*.c') -o build/v2_solver_vk -lm
 *   cd ocrt  <-- 필수. solver 는 inputs/ (water_iop, xsec, afgl_atm)를
 *                cwd 상대경로로 연다. 다른 디렉토리에서 실행하면
 *                "cannot open inputs/..." 로 전면 실패한다.
 *   OMP_NUM_THREADS=1  <-- 재현성 기준(모든 golden 은 단일스레드 실측값).
 *
 * ---------------------------------------------------------------------------
 * 1. 검증된 예제 커맨드 (전부 v1.09 세션 실측; 값까지 재현되어야 정상)
 * ---------------------------------------------------------------------------
 * [A] Tier-0 골든 (해수 fixed-bulk, 기대: rrs0minus=2.887685e-02 bit 일치)
 *   OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface ocean --wind-speed 3 \
 *     --sza 30 --vza 0 --raa 90 --wavelength 555 --pressure 0 --aod-555 0 \
 *     --water-model iop --iop-a 0.176200 --iop-b 2.618309 \
 *     --iop-bb 0.038578 \
 *     --iop-phase-lut <blend_P11.csv>
 *   (blend LUT 는 harness 가 Brown_earth.mie 로부터 재생성; 가이드라인 Tier 0 참조)
 *
 * [B] 에어로졸 OSOAA 대조 프로토콜 (실측 -0.3 ~ -1.3 % @555nm aot0.3 sza40)
 *   OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface black_fresnel_ocean --wind-speed 3 \
 *     --sza 40 --vza 30 --raa 90 --wavelength 555 --pressure 0 \
 *     --mie inputs/M80C.mie --aod-555 0.3 --aer-l-max 80 --m-max 16 \
 *     --aer-h-km 2.0 --n-layers 400 --trunc-aer-loglin \
 *     --aer-phase-kernel value --vector
 *   핵심: 절단은 loglin, 커널은 value 가 production. (delta-fit, delta-M,
 *   nt-tau 는 비검증 연구 경로 -> OCRT_DEBUG=1 게이트, 아래 4절)
 *
 *   [옵션 등급 안내 - v1.09 커밋 #13(2026-07-07): 검증 구성 기본값 승격]
 *   에어로졸 실행(--aod>0)은 검증 프로토콜 구성이 자동 기본이다:
 *     loglin 절단 + value 커널 + aer-l-max 80 + m-max 16 + n-layers 400
 *   따라서 최소 커맨드가 곧 프로토콜이다:
 *     ... --mie inputs/M80C.mie --aod-555 0.3 --aer-h-km 2.0 --vector
 *   (위 [B]의 명시 플래그 버전과 bit 동일 - 검증 완료)
 *   - 수렴 노브 오버라이드(--aer-l-max/--n-layers/--m-max 를 기본값과 다른
 *     값으로): OCRT_ADVANCED=1 필요. 기본값과 같은 값의 재명시는 게이트 없이
 *     통과(기존 문서 커맨드 호환).
 *   - (구판 기술) --vector 는 자동 아님 -- M1(2026-07-14)로 폐기, 벡터가 기본.
 *     2026-07-15 부로 --vector 토큰 자체가 삭제됨(아래 4절 참고).
 *   - Rayleigh-only 실행은 종전 기본(m=2 등) 그대로 - 무영향.
 *
 * [C] 기체 흡수 (Tier-1 골든: 555nm O3 column 9.321e18 = 346.9 DU,
 *     tau_O3 = 3.2292e-02;  602nm 4.8015e-02;  440nm NO2 1.7402e-03)
 *   ... [A]의 커맨드에 --use-absorption 추가.  xsec 는 continuum-merged 판 필수.
 *   컬럼 오버라이드: --gas-column-o3 <DU> 등 (AFGL 프로파일 형상 유지, 총량만 스케일)
 *
 * [D] 해수 3성분 OCRT 모델 (Chl/TSM/aDOM 필수 입력; 0 명시 허용)
 *   ... --surface ocean --water-model ocrt \
 *       --ocrt-chl 0.3 --ocrt-tsm 1.0 --ocrt-adom440 0.05 \
 *       --ocrt-phyto-group micro --ocrt-tsm-species red_clay
 *
 * [E] 대기 LUT (black 표면; 구판의 "--lut 는 --vector 필수" 규칙은 폐기,
 *     두 플래그 모두 2026-07-15 삭제 -- 격자가 기본 모드)
 *   ... --surface black --pressure 1013.25 --vector --lut \
 *       --lut-vza-step 10 --lut-raa-step 30 --lut-vza-max 60
 *   해수 LUT 는 --output-full-grid <파일> (ocean 한정 n_mu_water 자동 96 승격,
 *   단일코어에서 무겁다: 수 분).
 *
 * [F] 검증 하네스 (18케이스, MAPE 골든 1.90 %, 18/18 bit)
 *   cd ocrt && python3 scripts/ocrt_osoaa_consistency_harness.py
 *   주의 3가지: (1) 반드시 ocrt/ 에서 실행, (2) 스크립트 상단 config 경로
 *   (OCRT_BIN, OSOAA_ROOT, CONVERTER, MINERAL_MIE) 를 새 환경에 맞게 수정,
 *   (3) Tier-3 는 OSOAA 설치본이 필요 -- OCRT 단독 패키지에는 OSOAA 가
 *   포함되지 않을 수 있다(마이그레이션 풀패키지에는 포함).
 *
 * ---------------------------------------------------------------------------
 * 2. 출력 규약
 * ---------------------------------------------------------------------------
 *   단일 실행 stdout 1행:  rho_I rho_Q rho_U ... 이름=값 ...
 *   rho = I / mu_s (F_sun = pi 규약).  ocean 이면 rrs0minus / Rrs0plus [1/sr],
 *   Ed, Lu, Kd, Ku(진단용 -- production 산출물 아님) 등이 뒤따른다.
 *
 * ---------------------------------------------------------------------------
 * 3. 환경변수
 * ---------------------------------------------------------------------------
 *   OCRT_DEBUG=1     게이트 옵션 해제(아래 4절) + 진단 출력
 *   OCRT_ADVANCED=1  수렴 파라미터 오버라이드 허용(--n-mu-water 등)
 *   [2026-07-15 Jae 지시] --help 는 일반(무게이트) 옵션만 싣는다.  게이트
 *   옵션의 조건부 표시안은 폐기.  대신 --examples 가 검증된 실행 예시
 *   6종을 출력한다(usage_examples(), 예시 변경 시 실측 재검증 필수).
 *   [2026-07-15 Jae 지시, 옵션 정리 1차] --lut-output(--output-full-grid 와
 *   완전 중복) 삭제.  --output-mode {simple|debug} 삭제 -> --output-advanced
 *   플래그로 대체(기본=간단 열).  --view-as-node 삭제(단일 기하 모드가 강제
 *   활성하므로 순수 선언용이었음).
 *   [2026-07-15 Jae 결정, 옵션 정리 2차 = 1안] --vector/--lut 완전 삭제 +
 *   패키지 내부 스크립트(골든 러너 5, OSOAA 대조 5, 하니스 2, LUT 생산 1,
 *   ps1 도구 2) 동커밋 일괄 수정.  무동작 플래그 제거이므로 골든 값 불변을
 *   러너 비트 재통과로 실증함.  아카이브 사본(*_ORIGINAL_* 등)은 이력
 *   보존 목적으로 미수정 -- 실행 시 미지 인자 오류가 정상이다.
 *   [2026-07-15 Jae 확정, G2(대기) 정리] --tau-r --tau-r-from-input --dump-atm
 *   --rayleigh-model --xsec-dir --use-absorption 삭제.  레일리 = Bodhaine 고정
 *   (hansen-travis 구현 제거; 이로써 τ_R 기준 미결 사안은 Bodhaine 으로 종결,
 *   레거시 H-T 캠페인 값은 CLI 재현 불가 -- 수치는 v1.09 changelog 에 보존).
 *   기체흡수 기본 활성(off = 6기체 0, bit 등가 실측).  대조 스크립트·러너
 *   전 호출에 6기체 0 을 주입해 골든 bit 를 보존함.  --afgl-dir 도 동일
 *   원칙으로 삭제(inputs/afgl_atm 고정; userdef 는 그 안에 배치).
 *   [2026-07-15 Jae 확정, G3(에어로졸) 정리] --no-aerosol(무동작)
 *   --aer-phase-kernel --trunc-aer-loglin 삭제(자동 dispatch 일원화; loglin
 *   dispatch 를 결합 케이스까지 확장 -- 종전 결합 최소 명령은 loglin 미적용,
 *   I 0.51% @555 실측).  --aer-h-km 기본 2.0 지수형; legacy an23 분기 삭제.
 *   가드 신설: aod>0 은 --mie 필수(배치는 행 단위 SKIP), aod=0 에서 에어로졸
 *   옵션 지정은 오류.  .inp 직접 입력 인터페이스는 없음 -- mie_generator
 *   2단계(--mie-gen X.inp out.mie) 유지.
 *   [2026-07-15 Jae 확정, G4(해면·계면) 정리] --surface 필수화(조용한 black
 *   기본 폐지; ccrr-compare 예외).  --wind-speed 는 ocean/black_fresnel_ocean 에서 필수
 *   (0 명시 허용; 배치는 행 단위 SKIP).  --n-water 기본 = 1.34 고정(Quan-Fry
 *   자동 폐지 -- 종전 OSOAA 대조가 굴절률 불일치였음을 발견, 새 기본이 정합).
 *   --sigma-type/--q-convention -> ADVANCED(기본 1 재명시 통과, help 미수록).
 *   선글린트 반전: 기본 포함, --decouple-sunglint 지정 시 분리(구
 *   --no-decouple-sunglint 삭제).  골든 러너는 동결 물리를 명시 핀
 *   (per-λ Quan-Fry --n-water + --decouple-sunglint)하여 bit 보존.
 *   [2026-07-15 Jae 확정, G5 1차(물성·CDOM)] --f-sun 삭제(TOA 조도는 내부
 *   π 규약; 결합 모델에서 대기가 항상 선행 정의).  --water-aw-lut /
 *   --water-psi-T-lut 삭제(경로 고정).  CDOM 3종 세트 규칙 신설(전량 지정
 *   또는 전무; 부분 지정 오류).
 *   [2026-07-15 migration G5-2] SIMPLE 은 명시 --surface ocean 전용이며
 *   native CDOM/fixed-bulk/water-mie 와 배타.  surface·대기 aerosol·wind·
 *   sunglint 를 더 이상 강제 변경하지 않고 F_sun=1 동결 정규화만 유지.
 *   ccrr-particle-phase 부속 3종도 DEBUG 게이트와 일치시킴.
 *   OCRT_DUMP_SKY=1  skylight Ed/Lu 진단 덤프(결함 조사용)
 *
 * ---------------------------------------------------------------------------
 * 4. 게이트 일람 (기본 거부 -> 사유 출력; 이 파일 내 파싱부에 구현)
 *    [2026-07-14 M2/M3 갱신 — Jae 지시]
 * ---------------------------------------------------------------------------
 *   실행 모드(M2): --vza+--raa 동시 지정 = 단일 기하(view-as-node 자동).
 *     하나만 지정 = 오류.  미지정 = LUT 격자 모드(균일 간격, 기본 2.5도,
 *     --lut-vza-step/--lut-raa-step 공통, 해양-대기 통일).  출력 미지정 시
 *     result<N>.csv 자동, csv 외 확장자 오류, 모드 정책 충돌 시 오류 종료.
 *   벡터(M1): 항상 벡터.  스칼라 솔버 rt_solve_case 삭제(2026-07-14).
 *     --vector 는 호환용 no-op 이었음 -> 2026-07-15 완전 삭제.
 *   ADVANCED(기본값 재명시는 통과): --n-mu --n-layers --max-orders
 *     --sos-max-orders --sos-tolerance --conv-tol --l-max --m-max
 *     --water-m-max --n-mu-water --iop-phase-nphi + 에어로졸 노브.
 *   DEBUG: --scalar-ff-* --ccrr-particle-phase-* --simple-compare,
 *     --sos-save-orders --debug-* 계열, --integration-method constant,
 *     --trunc-aer-fit/-m, --nt-tau, "M50C" mie 경로.
 *   항상 거부: --sos-acceleration geometric (S-003), 위상 각도표 [0,180]
 *     미커버 (S-002).
 *
 * ---------------------------------------------------------------------------
 * 5. 입력 파일
 * ---------------------------------------------------------------------------
 *   .mie      : 6SV 형 다중밴드(헤더 Ext/Sca/ssa/asym 표 + P11/P12/P33 각도블록,
 *               각도 0~180 완전 커버 필수).  공식 에어로졸 subset: T50, C50, M80C.
 *   위상 LUT  : CSV "theta_deg,P11" (fixed-bulk 용)
 *   생성기    : mie_generator/vrt_solver --mie-gen X.inp out.mie --dtheta 0.5
 *               (.inp 원본 세트는 v1.09 패키지에 부재 -- MANIFEST 참조)
 *
 * 유지보수 규칙: 옵션/게이트/골든 값을 바꾸면 이 블록과 --help,
 * docs/OCRT_USAGE 를 같은 커밋에서 갱신한다.  주석 안에 별표+슬래시 연쇄를
 * 쓰지 말 것(블록이 조기 종료된다 -- v1.09 세션 실사고).
 * ========================================================================== */

/* 2026-07-14 M2 (Jae 지시) 헬퍼: 출력 파일명 규칙.
 * 확장자 없음 -> ".csv" 자동 부여, ".csv" -> 그대로, 그 외 확장자 -> 오류.
 * 파일명 미지정 -> 실행 디렉터리의 result<N>.csv 중 최대 N+1 자동 선택. */
static int ocrt_resolve_csv_name(const char *given, char *out, size_t cap) {
    if (given) {
        const char *base = strrchr(given, '/');
        base = base ? base + 1 : given;
        const char *dot = strrchr(base, '.');
        if (!dot) { snprintf(out, cap, "%s.csv", given); return 0; }
        if (strcmp(dot, ".csv") == 0) { snprintf(out, cap, "%s", given); return 0; }
        fprintf(stderr, "error: output extension '%s' is not supported (csv only)\n", dot);
        return 2;
    }
    long maxn = 0;
    DIR *d = opendir(".");
    if (d) {
        struct dirent *e;
        while ((e = readdir(d)) != NULL) {
            long n; char tail[8] = {0};
            if (sscanf(e->d_name, "result%ld%7s", &n, tail) == 2 &&
                strcmp(tail, ".csv") == 0 && n > maxn) maxn = n;
        }
        closedir(d);
    }
    snprintf(out, cap, "result%ld.csv", maxn + 1);
    return 0;
}

static int ocrt_advanced_on(void) {
    const char *e = getenv("OCRT_ADVANCED");
    return (e && *e && strcmp(e, "0") != 0);
}

static int ocrt_debug_on(void) {
    const char *d = getenv("OCRT_DEBUG");
    return (d && *d && strcmp(d, "0") != 0);
}

static int ocrt_tsm_data_dir_complete(const char *dir) {
    static const char *const files[] = {
        "Red_clay_AHN.mie", "Brown_earth_AHN.mie",
        "Yellow_clay_AHN.mie", "Calcareous_sand_AHN.mie"
    };
    if (!dir || !dir[0]) return 0;
    for (size_t i = 0; i < sizeof files / sizeof files[0]; ++i) {
        char path[2048];
        struct stat st;
        int n = snprintf(path, sizeof path, "%s/%s", dir, files[i]);
        if (n < 0 || (size_t)n >= sizeof path || stat(path, &st) != 0 || !S_ISREG(st.st_mode))
            return 0;
    }
    return 1;
}

static int ocrt_tsm_try_init(const char *dir, char *resolved, size_t cap) {
    if (!ocrt_tsm_data_dir_complete(dir)) return -1;
    if (rt_iop_ahn_mineral_init(dir) != 0) return -2;
    if (resolved && cap > 0) snprintf(resolved, cap, "%s", dir);
    return 0;
}

/* Resolve the packaged TSM dataset from an explicit environment override,
 * the usual package working directories, or a build/ executable path. */
static int ocrt_tsm_init_for_run(const char *argv0, char *resolved, size_t cap) {
    const char *env = getenv("OCRT_TSM_DIR");
    if (env && env[0]) {
        int rc = ocrt_tsm_try_init(env, resolved, cap);
        if (rc != 0)
            fprintf(stderr, "error: OCRT_TSM_DIR '%s' is missing or malformed\n", env);
        return rc;
    }

    const char *const candidates[] = {
        "inputs/tsm_ahn", "./inputs/tsm_ahn", "../inputs/tsm_ahn",
        "ocrt/inputs/tsm_ahn", "../ocrt/inputs/tsm_ahn"
    };
    for (size_t i = 0; i < sizeof candidates / sizeof candidates[0]; ++i)
        if (ocrt_tsm_try_init(candidates[i], resolved, cap) == 0) return 0;

    if (argv0 && argv0[0]) {
        const char *slash1 = strrchr(argv0, '/');
        const char *slash2 = strrchr(argv0, '\\');
        const char *slash = slash1;
        if (!slash || (slash2 && slash2 > slash)) slash = slash2;
        if (slash) {
            char exe_dir[1536];
            size_t len = (size_t)(slash - argv0);
            if (len < sizeof exe_dir) {
                memcpy(exe_dir, argv0, len);
                exe_dir[len] = 0;
                char candidate[2048];
                snprintf(candidate, sizeof candidate, "%s/../inputs/tsm_ahn", exe_dir);
                if (ocrt_tsm_try_init(candidate, resolved, cap) == 0) return 0;
                snprintf(candidate, sizeof candidate, "%s/inputs/tsm_ahn", exe_dir);
                if (ocrt_tsm_try_init(candidate, resolved, cap) == 0) return 0;
            }
        }
    }
    fprintf(stderr, "error: cannot locate Ahn TSM data; run from the OCRT directory or set OCRT_TSM_DIR\n");
    return -1;
}

static int ocrt_regular_file_in_dir(const char *dir, const char *name) {
    char path[2048];
    struct stat st;
    int n = snprintf(path, sizeof path, "%s/%s", dir, name);
    return n >= 0 && (size_t)n < sizeof path && stat(path, &st) == 0 && S_ISREG(st.st_mode);
}

static int ocrt_organic_data_dir_complete(const char *dir,
                                          organic_phyto_group_t selected_group,
                                          int load_all_phyto) {
    /* EAP_RUNTIME_DISABLED: Chl production requires only the common absorption
     * table plus the Stramski detritus phase.  EAP files are not completeness
     * criteria and are never loaded. */
    (void)selected_group;
    (void)load_all_phyto;
    return dir && dir[0] &&
           ocrt_regular_file_in_dir(dir, "phyto_absorption_default.csv") &&
           ocrt_regular_file_in_dir(dir, "Detritus_Stramski2001.mie");
}

static int ocrt_organic_try_init(const char *dir,
                                 organic_phyto_group_t selected_group,
                                 int load_all_phyto,
                                 char *resolved, size_t cap) {
    if (!ocrt_organic_data_dir_complete(dir, selected_group, load_all_phyto)) return -1;
    if (rt_iop_organic_init_selected(dir, selected_group, load_all_phyto) != 0) return -2;
    if (resolved && cap > 0) snprintf(resolved, cap, "%s", dir);
    return 0;
}

/* Resolve the packaged common-absorption/Stramski dataset. OCRT_ORGANIC_DIR is an
 * explicit override; the normal package location is inputs/water_iop. */
static int ocrt_organic_init_for_run(const char *argv0,
                                     organic_phyto_group_t selected_group,
                                     int load_all_phyto,
                                     char *resolved, size_t cap) {
    const char *env = getenv("OCRT_ORGANIC_DIR");
    if (env && env[0]) {
        int rc = ocrt_organic_try_init(env, selected_group, load_all_phyto, resolved, cap);
        if (rc != 0)
            fprintf(stderr, "error: OCRT_ORGANIC_DIR '%s' is missing or malformed\n", env);
        return rc;
    }
    const char *const candidates[] = {
        "inputs/water_iop", "./inputs/water_iop", "../inputs/water_iop",
        "ocrt/inputs/water_iop", "../ocrt/inputs/water_iop"
    };
    for (size_t i = 0; i < sizeof candidates / sizeof candidates[0]; ++i)
        if (ocrt_organic_try_init(candidates[i], selected_group, load_all_phyto, resolved, cap) == 0) return 0;

    if (argv0 && argv0[0]) {
        const char *slash1 = strrchr(argv0, '/');
        const char *slash2 = strrchr(argv0, '\\');
        const char *slash = slash1;
        if (!slash || (slash2 && slash2 > slash)) slash = slash2;
        if (slash) {
            char exe_dir[1536];
            size_t len = (size_t)(slash - argv0);
            if (len < sizeof exe_dir) {
                memcpy(exe_dir, argv0, len);
                exe_dir[len] = 0;
                char candidate[2048];
                snprintf(candidate, sizeof candidate, "%s/../inputs/water_iop", exe_dir);
                if (ocrt_organic_try_init(candidate, selected_group, load_all_phyto, resolved, cap) == 0) return 0;
                snprintf(candidate, sizeof candidate, "%s/inputs/water_iop", exe_dir);
                if (ocrt_organic_try_init(candidate, selected_group, load_all_phyto, resolved, cap) == 0) return 0;
            }
        }
    }
    fprintf(stderr,
            "error: cannot locate OCRT organic data; run from the OCRT directory or set OCRT_ORGANIC_DIR\n");
    return -1;
}

int main(int argc, char **argv) {
    rt_options_t opts = rt_options_default();
    rt_case_t    cs   = {0};
    cs.surface     = RT_SURFACE_BLACK;
    cs.pressure_hpa = 1013.25;  /* standard atmosphere => Rayleigh present by default. */
    cs.rayleigh_on = 1;         /* placeholder; recomputed from pressure_hpa after arg parse. */
    cs.aerosol_on  = 0;         /* placeholder; recomputed from AOD after arg parse. */
    /* Phase 4 surface defaults — used only when surface != BLACK */
    cs.n_water         = 1.34;
    cs.wind_speed      = 0.0;
    cs.sigma_type      = 1;   /* Cox-Munk 1954 (matches AF1982 reference) */
    cs.q_convention    = 1;   /* Hansen form: Q_kernel=0.5(Rp-Rs), V3 internal RT (gamma2=-ron*sqrt(1.5)) 와 일관 */
    cs.decouple_sunglint = 0; /* 2026-07-15 Jae(G4) 반전: 기본은 선글린트 포함.
                               * 분리는 --decouple-sunglint 명시(AF1982 정합용). */

    /* Phase B (in-water RT) — used only when surface == RT_SURFACE_OCEAN */
    cs.T_water_C    = 20.0;
    cs.S_water_gkg  = 38.4;
    cs.F_sun        = M_PI;   /* OCRT internal convention; user can override */

    /* aDOM defaults shared by the OCRT/CCRR constituent adapters. */
    cs.a_cdom_440_m_inv   = 0.0;
    cs.S_cdom_nm_inv      = 0.014;     /* Bricaud 1981 classical value, nm⁻¹ */
    cs.cdom_ref_lambda_nm = 440.0;
    cs.n_layers_water_override = 0;    /* debug-only: 0 → production automatic */
    cs.tau_max_target_override = 0.0;  /* debug-only base τ: <=0 → production automatic */
    cs.water_max_tau_override = 0.0;   /* debug-only max τ cap: <=0 → production automatic */
    cs.water_max_z_max_m_override = 0.0; /* debug-only max depth cap: <=0 → production automatic */
    cs.water_depth_bottom_tol_override = 0.0; /* debug-only bottom attenuation target */
    cs.water_layer_dtau_target_override = 0.0; /* debug-only: <=0 → production automatic */
    cs.n_mu_water_override = 0;
    cs.water_shared_grid   = 0;        /* #22: --water-shared-grid */        /* 0 → use default 48 (>=48 recommended for in-water,
                                        * incl. LUT mode; n_mu=32 under-converged ~1.5%); CLI --n-mu-water */
    cs.n_mu_water_sky = 0;             /* 0 → same as n_mu_water; CLI --n-mu-water-sky */
    cs.water_m_max_override = 0;       /* 0 → auto; CLI --water-m-max */
    cs.water_max_orders_override = 0;  /* debug-only: 0 → solver automatic cap */
    cs.debug_high_particle_water_max_orders = 0; /* 0 → production default 100 for particle-rich media; debug-only override */
    cs.water_view_as_node = 1;        /* default: continuous reconstruction (CCRR-compatible) */
    cs.water_input_mode = RT_WATER_INPUT_UNSET;
    cs.ccrr_mode = 0;
    cs.water_constituent_model = RT_WATER_CONSTITUENT_OCRT;
    cs.ccrr_chl_mg_m3 = 0.0;
    cs.ccrr_min_g_m3 = 0.0;
    cs.tsm_species = AHN_RED_CLAY;
    cs.organic_phyto_group = ORGANIC_PHYTO_MICRO;
    cs.detritus_a440_m_inv = 0.0;
    cs.detritus_slope_nm_inv = ORGANIC_DETRITUS_SLOPE_DEFAULT;
    cs.ccrr_phase_moments_path = NULL;
    cs.ccrr_particle_phase_lut_path = NULL;
    cs.ccrr_particle_phase_case_id = NULL;
    cs.ccrr_particle_phase_wavelength_nm = 0.0;
    cs.ccrr_particle_phase_lmax = 10;
    cs.ccrr_particle_phase_nphi = 720;
    cs.water_phase_kernel = 0;   /* internal legacy kernel selector; no public CLI override */
    cs.fixed_bulk_iop_mode = 0;
    cs.fixed_a_total_m_inv = 0.0;
    cs.fixed_b_total_m_inv = 0.0;
    cs.fixed_bb_total_m_inv = 0.0;
    cs.fixed_bulk_lmax = 0;   /* #23: 0 = per-path default (water-mie 200, LUT paths 30) */
    cs.fixed_bulk_phase_model = 0;
    cs.fixed_bulk_phase_nphi = 720;
    cs.fixed_bulk_phase_lut_path = NULL;
    cs.water_mie_phase_path = NULL;
    cs.water_mie_moment_mode = 1;
    cs.water_mie_truncation_mode = 0;
    cs.water_mie_ss_mode = 0;
    cs.organic_phyto_scattering = 0;  /* EAP scattering disabled; absorption retained */
    cs.water_mie_moment_n_mu = 400;   /* #23: was 60 — undersampled moments above L~60
                                       * (measured contamination); 400 supports L200. */
    cs.scalar_ff_iop_path = NULL;
    cs.scalar_ff_bb_over_b = 0.0;
    cs.scalar_ff_refractive_index = 1.18;
    cs.scalar_ff_mu = 0.0;
    cs.fixed_bulk_phase_case_id = NULL;
    cs.fixed_bulk_phase_wavelength_nm = 0.0;

    /* Phase B ocean-mode LUT paths (defaults from inputs/water_iop/) */
    const char *water_aw_lut_path    = "inputs/water_iop/water_coef_z09_1nm.txt";
    const char *water_psi_T_lut_path = "inputs/water_iop/psi_T_rottgers2014_OCRT.txt";

    /* Aerosol path inputs (Phase 3, optional) */
    const char *batch_full_grid_in = NULL;   /* v1.09 #19 */
    const char *mie_path = NULL;
    double      user_aod = 0.0;
    double      user_aod_ref_nm = 555.0;
    double      aod_target_used = 0.0;
    double      aod_ratio_used = 0.0;
    int         aer_L_max = 32;
    int         have_aer_L_max = 0;  /* T-V7-2026-05-10: dispatch flag */
    int         have_aer_h = 0;      /* 2026-07-15 G3: aod=0 게이트 판정용 */
    int         have_surface = 0;    /* 2026-07-15 G4: --surface 필수화 */
    int         have_wind = 0;       /* 2026-07-15 G4: ocean/black_fresnel_ocean 풍속 필수화 */
    int         water_model_seen = 0;      /* --water-model supplied exactly once */
    unsigned    water_option_family_mask = 0u;
    int         water_pure_from_zero = 0;  /* OCRT/CCRR selected with Chl=TSM=aDOM440=0 */
    int         simple_adom_s_explicit = 0; /* validate an explicitly supplied slope */
    int         tsm_species_cli_seen = 0;  /* --ocrt-tsm-species supplied */
    int         organic_group_cli_seen = 0;
    int         detritus_cli_seen = 0;
    int         chl_input_seen = 0, tsm_input_seen = 0, adom_input_seen = 0;
    int         iop_a_input_seen = 0, iop_b_input_seen = 0, iop_bb_input_seen = 0;
    int         ccrr_particle_lut_seen = 0;
    int         ccrr_particle_selector_seen = 0;
    int         fixed_bulk_cli_seen = 0;    /* direct bulk-IOP/phase injection branch */
    int         ocrt_mie_control_seen = 0;
    int         iop_mie_control_seen = 0;
    int         ocrt_mie_trunc_seen = 0, ocrt_mie_ss_seen = 0;
    int         iop_mie_trunc_seen = 0, iop_mie_ss_seen = 0;
    int         iop_phase_lut_seen = 0, iop_mie_phase_seen = 0;
    int         have_aer_phase_kernel = 0; /* v1.09 #13: dispatch flag */
    int         have_n_layers = 0;         /* v1.09 #13: dispatch flag */
    double      theta_cut_deg = 0.0;  /* 0 = no delta-fit truncation */
    int         delta_m_N = 0;        /* 0 = no delta-M truncation */
    int         have_trunc_aer_m = 0;  /* user explicitly requested --trunc-aer-m */
    int         have_trunc_aer_fit = 0;/* user explicitly requested --trunc-aer-fit */
    int         apply_nt_tau  = 0;    /* 1 = Wiscombe NT τ transform if f>0 */
    int         use_loglin_trunc = 0; /* 1 = real-angle log10-linear forward-peak truncation */
    int         aer_use_value_kernel = 0; /* --aer-phase-kernel value => Gibbs-free angle-space aerosol phase (default: moment) */

    /* v1.01 LUT mode options */
    int         lut_mode      = 0;     /* 1 → produce angular grid output */
    const char *lut_output    = NULL;  /* file path for LUT CSV; NULL = stdout */
    double      lut_vza_step  = 2.5;
    double      lut_raa_step  = 2.5;   /* M2 (2026-07-14): 기본 격자 간격 2.5도 통일 */
    double      lut_vza_max   = 85.0;  /* 87.5+ degrades; cap at 85 */

    /* v1.01 atmospheric absorption options */
    /* 2026-07-15 Jae 지시(G2 확정): 기체흡수 기본 활성.  구 --use-absorption
     * 플래그 삭제.  완전 off = 6기체 총량 전부 0 지정(레일리 전용 결과와
     * bit 동일 -- 2026-07-15 실측).  압력 0 이어도 기체흡수는 독립 적용되나
     * 비율량(rrs/Rrs/Kd)은 불변이므로 Tier-0 앵커는 유지된다(실측). */
    int                 use_absorption = 1;
    rt_afgl_profile_t   atm_profile    = RT_AFGL_USSTD76;
    const char         *afgl_dir       = "inputs/afgl_atm";
    const char         *xsec_dir       = "inputs/xsec";
    double              gas_col_override[RT_N_GAS] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};

    const char *batch_in = NULL;
    const char *batch_out = NULL;
    const char *ccrr_compare_in = NULL;
    const char *ccrr_compare_out = NULL;

    int have_sza = 0, have_vza = 0, have_raa = 0, have_wl = 0;
    int have_vza_step = 0, have_raa_step = 0;   /* M2: 격자 간격 명시 여부 */
    int single_mode_flag = 0;                   /* M2 잔재.  2026-07-15 --view-as-node 삭제로
                                                 * 항상 0 이나, 모드 판정식의 형태 보존을 위해
                                                 * 변수는 존치(상수 0). */
    int have_n_water = 0;  /* Step57: distinguish explicit --n-water from spectral auto Quan-Fry */
    /* T-V6-2026-05-10: 사용자 가 명시 적으로 --m-max 지정 했는지 추적.
     * argument loop 종료 후 dispatch: 미명시 시 aerosol_on 에 따라
     * Rayleigh-only=2 / aerosol-on=16 자동 설정. */
    int have_m_max = 0;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--debug")) {
            /* master debug gate: must precede any --debug-* option or debug env-var.
             * Sets OCRT_DEBUG so every ocrt_debug_env()/getenv("OCRT_DEBUG") gate
             * (CLI debug options, diagnostic dumps, internal overrides) is honored.
             * Production/validation runs must omit --debug. */
            setenv("OCRT_DEBUG", "1", 1);
            continue;
        }
        if (!strcmp(a, "--help") || !strcmp(a, "-h")) {
            usage(stdout, argv[0]); return 0;
        } else if (!strcmp(a, "--examples")) {
            usage_examples(stdout, argv[0]); return 0;
        } else if (!strcmp(a, "--version")) {
            printf("ocrt " V2_VERSION "\n"); return 0;
        } else if (!strcmp(a, "--sza")) {
            cs.sza_deg = next_double(argc, argv, &i, a); have_sza = 1;
        } else if (!strcmp(a, "--vza")) {
            cs.vza_deg = next_double(argc, argv, &i, a); have_vza = 1;
        } else if (!strcmp(a, "--raa")) {
            cs.raa_deg = next_double(argc, argv, &i, a); have_raa = 1;
        } else if (!strcmp(a, "--wavelength")) {
            cs.wavelength_nm = next_double(argc, argv, &i, a); have_wl = 1;
        } else if (!strcmp(a, "--surface")) {
            const char *s = next_str(argc, argv, &i, a);
            if      (!strcmp(s, "black"))   cs.surface = RT_SURFACE_BLACK;
            else if (!strcmp(s, "lambert")) {
                /* T0-5 (Phase A, 2026-05-09): Lambert surface CLI 차단.
                 * enum 정의되어 있으나 surface-aware solver routing 안 됨.
                 * silent black fallback 으로 sample BC 누락 위험 — 명시적 error. */
                fprintf(stderr,
                        "error: --surface lambert not implemented "
                        "(boundary condition not coded). "
                        "Use --surface black for atm-only or --surface flat/black_fresnel_ocean for a Fresnel ocean boundary.\n");
                return 2;
            }
            else if (!strcmp(s, "black_fresnel_ocean") ||
                     !strcmp(s, "black-fresnel-ocean"))
                cs.surface = RT_SURFACE_BLACK_FRESNEL_OCEAN;
            else if (!strcmp(s, "coxmunk")) {
                fprintf(stderr,
                        "warning: --surface coxmunk is deprecated; use "
                        "--surface black_fresnel_ocean. The surface mode denotes "
                        "a rough Fresnel interface over a black ocean; the slope "
                        "variance model is selected separately.\n");
                cs.surface = RT_SURFACE_BLACK_FRESNEL_OCEAN;
            }
            else if (!strcmp(s, "flat"))    cs.surface = RT_SURFACE_FLAT;
            else if (!strcmp(s, "ocean")) {
                /* Phase B.4 Stage 1 (2026-05-22): in-water RT solver integrated.
                 * Atmosphere coupling is Stage 2 (next milestone) — currently
                 * the in-water solver runs standalone with F_sun as just-
                 * above-surface BOA irradiance (T_aw=1 assumption). */
                cs.surface = RT_SURFACE_OCEAN;
            }
            else { fprintf(stderr, "error: unknown surface '%s'\n", s); return 2; }
            have_surface = 1;
        } else if (!strcmp(a, "--wind-speed")) {
            cs.wind_speed = next_double(argc, argv, &i, a);
            have_wind = 1;
        } else if (!strcmp(a, "--n-water")) {
            cs.n_water = next_double(argc, argv, &i, a);
            have_n_water = 1;
        } else if (!strcmp(a, "--water-temperature")) {
            cs.T_water_C = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--water-salinity")) {
            cs.S_water_gkg = next_double(argc, argv, &i, a);
        /* 2026-07-15 Jae(G5): --water-aw-lut / --water-psi-T-lut 삭제 --
         * 사용자 불변 입력 경로(inputs/water_iop/ 고정, G2 의 xsec/afgl 과
         * 동일 원칙). */
        /* 2026-07-15 Jae(G5): --f-sun 삭제.  결합 모델에서 TOA 조도는 내부
         * 규약(π; 반사도 정규화)이며 사용자가 만질 값이 아니다.  대기는
         * ocean 실행에 항상 선행 정의된다(--pressure 기본 1013.25, 기체흡수
         * 기본 활성).  SIMPLE/고정IOP 경로의 내부 F_sun=1.0 강제는 동결
         * 물리 일부로 존치(Tier-0 앵커 의존). */
        } else if (!strcmp(a, "--cdom-a440") ||
                   !strcmp(a, "--cdom-slope") ||
                   !strcmp(a, "--cdom-ref-lambda")) {
            fprintf(stderr,
                    "error: legacy native CDOM option '%s' was removed; select --water-model ocrt|ccrr and use --ocrt-adom* or --ccrr-adom*\n",
                    a);
            return 2;
        } else if (!strcmp(a, "--water-model")) {
            rt_water_input_mode_t mode = RT_WATER_INPUT_UNSET;
            const char *m = next_str(argc, argv, &i, a);
            if (water_model_seen) {
                fprintf(stderr,
                        "error: --water-model must be specified exactly once (duplicate value '%s')\n",
                        m);
                return 2;
            }
            if (water_input_mode_parse(m, &mode) != 0) {
                fprintf(stderr, "error: --water-model must be ocrt, ccrr or iop\n");
                return 2;
            }
            water_model_seen = 1;
            cs.water_input_mode = (int)mode;
            if (mode == RT_WATER_INPUT_OCRT)
                cs.water_constituent_model = RT_WATER_CONSTITUENT_OCRT;
            else if (mode == RT_WATER_INPUT_CCRR)
                cs.water_constituent_model = RT_WATER_CONSTITUENT_CCRR;
        } else if (!strcmp(a, "--ocrt-chl") || !strcmp(a, "--ccrr-chl")) {
            if (chl_input_seen) {
                fprintf(stderr, "error: Chl was specified more than once; use exactly one model-prefixed Chl option\n");
                return 2;
            }
            chl_input_seen = 1;
            water_option_family_mask |= !strcmp(a, "--ocrt-chl")
                                          ? WATER_OPTION_FAMILY_OCRT
                                          : WATER_OPTION_FAMILY_CCRR;
            cs.ccrr_chl_mg_m3 = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--ocrt-tsm") || !strcmp(a, "--ccrr-tsm")) {
            if (tsm_input_seen) {
                fprintf(stderr, "error: TSM was specified more than once; use exactly one model-prefixed TSM option\n");
                return 2;
            }
            tsm_input_seen = 1;
            water_option_family_mask |= !strcmp(a, "--ocrt-tsm")
                                          ? WATER_OPTION_FAMILY_OCRT
                                          : WATER_OPTION_FAMILY_CCRR;
            cs.ccrr_min_g_m3 = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--ocrt-adom440") || !strcmp(a, "--ccrr-adom440")) {
            if (adom_input_seen) {
                fprintf(stderr, "error: aDOM440 was specified more than once; use exactly one model-prefixed aDOM440 option\n");
                return 2;
            }
            adom_input_seen = 1;
            water_option_family_mask |= !strcmp(a, "--ocrt-adom440")
                                          ? WATER_OPTION_FAMILY_OCRT
                                          : WATER_OPTION_FAMILY_CCRR;
            cs.a_cdom_440_m_inv = next_double(argc, argv, &i, a);
            cs.cdom_ref_lambda_nm = 440.0;
        } else if (!strcmp(a, "--ocrt-adom-slope") || !strcmp(a, "--ccrr-adom-slope")) {
            if (simple_adom_s_explicit) {
                fprintf(stderr, "error: aDOM slope was specified more than once\n");
                return 2;
            }
            simple_adom_s_explicit = 1;
            water_option_family_mask |= !strcmp(a, "--ocrt-adom-slope")
                                          ? WATER_OPTION_FAMILY_OCRT
                                          : WATER_OPTION_FAMILY_CCRR;
            cs.S_cdom_nm_inv = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--ocrt-phyto-group")) {
            organic_phyto_group_t group;
            organic_group_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_OCRT;
            if (rt_iop_organic_group_parse(next_str(argc, argv, &i, a), &group) != 0) {
                fprintf(stderr,
                        "error: unrecognized --ocrt-phyto-group value; use a legacy "
                        "pico/nano/micro name or a canonical eap_* species identifier\n");
                return 2;
            }
            cs.organic_phyto_group = (int)group;
        } else if (!strcmp(a, "--ocrt-tsm-species")) {
            ahn_species_t species = AHN_RED_CLAY;
            const char *name = next_str(argc, argv, &i, a);
            tsm_species_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_OCRT;
            if (ahn_species_parse_strict(name, &species) != 0) {
                fprintf(stderr,
                        "error: --ocrt-tsm-species must be red_clay, brown_earth, yellow_clay or calcareous_sand\n");
                return 2;
            }
            cs.tsm_species = (int)species;
        } else if (!strcmp(a, "--ocrt-detritus-a440")) {
            detritus_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_OCRT;
            cs.detritus_a440_m_inv = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--ocrt-detritus-slope")) {
            detritus_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_OCRT;
            cs.detritus_slope_nm_inv = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--ccrr-phase-moments")) {
            water_option_family_mask |= WATER_OPTION_FAMILY_CCRR;
            cs.ccrr_phase_moments_path = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--ccrr-particle-phase-lut")) {
            if (!ocrt_debug_on()) {
                fprintf(stderr, "error: %s is a diagnostic option and requires OCRT_DEBUG=1\n", a);
                return 2;
            }
            water_option_family_mask |= WATER_OPTION_FAMILY_CCRR;
            ccrr_particle_lut_seen = 1;
            cs.ccrr_particle_phase_lut_path = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--ccrr-particle-phase-case-id")) {
            if (!ocrt_debug_on()) {
                fprintf(stderr, "error: %s is a diagnostic option and requires OCRT_DEBUG=1\n", a);
                return 2;
            }
            water_option_family_mask |= WATER_OPTION_FAMILY_CCRR;
            ccrr_particle_selector_seen = 1;
            cs.ccrr_particle_phase_case_id = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--ccrr-particle-phase-lmax")) {
            if (!ocrt_debug_on()) {
                fprintf(stderr, "error: %s is a diagnostic option and requires OCRT_DEBUG=1\n", a);
                return 2;
            }
            water_option_family_mask |= WATER_OPTION_FAMILY_CCRR;
            ccrr_particle_selector_seen = 1;
            cs.ccrr_particle_phase_lmax = next_int(argc, argv, &i, a);
        } else if (!strcmp(a, "--ccrr-particle-phase-nphi")) {
            if (!ocrt_debug_on()) {
                fprintf(stderr, "error: %s is a diagnostic option and requires OCRT_DEBUG=1\n", a);
                return 2;
            }
            water_option_family_mask |= WATER_OPTION_FAMILY_CCRR;
            ccrr_particle_selector_seen = 1;
            cs.ccrr_particle_phase_nphi = next_int(argc, argv, &i, a);
        } else if (!strcmp(a, "--simple-mode") || !strcmp(a, "--ccrr-mode") ||
                   !strcmp(a, "--chl") || !strcmp(a, "--tsm") ||
                   !strcmp(a, "--adom440") || !strcmp(a, "--adom-slope") ||
                   !strcmp(a, "--simple-chl") || !strcmp(a, "--simple-min") ||
                   !strcmp(a, "--simple-adom440") || !strcmp(a, "--simple-adom-s") ||
                   !strcmp(a, "--simple-adom-slope") ||
                   !strcmp(a, "--ccrr-chl-mg-m3") || !strcmp(a, "--ccrr-min") ||
                   !strcmp(a, "--ccrr-min-g-m3") || !strcmp(a, "--ccrr-adom-440") ||
                   !strcmp(a, "--ccrr-adom-s") || !strcmp(a, "--phyto-group") ||
                   !strcmp(a, "--chl-group") || !strcmp(a, "--tsm-species") ||
                   !strcmp(a, "--detritus-a440") || !strcmp(a, "--detritus-slope")) {
            fprintf(stderr,
                    "error: legacy/unprefixed water option '%s' was removed; select --water-model ocrt|ccrr|iop and use --ocrt-* or --ccrr-* options\n",
                    a);
            return 2;
        } else if (!strcmp(a, "--iop-a")) {
            if (iop_a_input_seen) { fprintf(stderr, "error: --iop-a was specified more than once\n"); return 2; }
            iop_a_input_seen = 1;
            fixed_bulk_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            cs.fixed_a_total_m_inv = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--iop-b")) {
            if (iop_b_input_seen) { fprintf(stderr, "error: --iop-b was specified more than once\n"); return 2; }
            iop_b_input_seen = 1;
            fixed_bulk_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            cs.fixed_b_total_m_inv = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--iop-bb")) {
            if (iop_bb_input_seen) { fprintf(stderr, "error: --iop-bb was specified more than once\n"); return 2; }
            iop_bb_input_seen = 1;
            fixed_bulk_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            cs.fixed_bb_total_m_inv = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--fixed-bulk-iop") || !strcmp(a, "--water-fixed-bulk-iop")) {
            fprintf(stderr,
                    "error: legacy option '%s' was removed; use --water-model iop with --iop-a --iop-b --iop-bb\n",
                    a);
            return 2;
        } else if (!strcmp(a, "--debug-fixed-bulk-lmax")) {
            fixed_bulk_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-fixed-bulk-lmax requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.fixed_bulk_lmax = next_int(argc, argv, &i, a);
        } else if (!strcmp(a, "--debug-high-particle-water-max-orders")) {
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-high-particle-water-max-orders requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.debug_high_particle_water_max_orders = next_int(argc, argv, &i, a);
            if (cs.debug_high_particle_water_max_orders <= 0) {
                fprintf(stderr, "error: --debug-high-particle-water-max-orders must be positive\n");
                return 2;
            }
        } else if (!strcmp(a, "--debug-water-max-orders")) {
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-water-max-orders requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.water_max_orders_override = next_int(argc, argv, &i, a);
            if (cs.water_max_orders_override <= 0) {
                fprintf(stderr, "error: --debug-water-max-orders must be positive\n");
                return 2;
            }
        } else if (!strcmp(a, "--debug-water-tau-max-target")) {
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-water-tau-max-target requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.tau_max_target_override = next_double(argc, argv, &i, a);
            if (cs.tau_max_target_override <= 0.0) {
                fprintf(stderr, "error: --debug-water-tau-max-target must be positive\n");
                return 2;
            }
        } else if (!strcmp(a, "--debug-water-max-tau")) {
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-water-max-tau requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.water_max_tau_override = next_double(argc, argv, &i, a);
            if (cs.water_max_tau_override <= 0.0) {
                fprintf(stderr, "error: --debug-water-max-tau must be positive\n");
                return 2;
            }
        } else if (!strcmp(a, "--debug-water-max-depth") || !strcmp(a, "--debug-water-min-depth")) {
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-water-max-depth requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.water_max_z_max_m_override = next_double(argc, argv, &i, a);
            if (cs.water_max_z_max_m_override <= 0.0) {
                fprintf(stderr, "error: --debug-water-max-depth must be positive\n");
                return 2;
            }
        } else if (!strcmp(a, "--debug-water-bottom-tol")) {
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-water-bottom-tol requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.water_depth_bottom_tol_override = next_double(argc, argv, &i, a);
            if (cs.water_depth_bottom_tol_override <= 0.0 || cs.water_depth_bottom_tol_override >= 1.0) {
                fprintf(stderr, "error: --debug-water-bottom-tol must be in (0,1)\n");
                return 2;
            }
        } else if (!strcmp(a, "--debug-water-layer-dtau")) {
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-water-layer-dtau requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.water_layer_dtau_target_override = next_double(argc, argv, &i, a);
            if (cs.water_layer_dtau_target_override <= 0.0) {
                fprintf(stderr, "error: --debug-water-layer-dtau must be positive\n");
                return 2;
            }
        } else if (!strcmp(a, "--debug-water-n-layers")) {
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-water-n-layers requires OCRT_DEBUG=1\n");
                return 2;
            }
            cs.n_layers_water_override = next_int(argc, argv, &i, a);
            if (cs.n_layers_water_override <= 0) {
                fprintf(stderr, "error: --debug-water-n-layers must be positive\n");
                return 2;
            }
        } else if (!strcmp(a, "--debug-fixed-bulk-phase-model")) {
            fixed_bulk_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            const char *m = next_str(argc, argv, &i, a);
            const char *dbg = getenv("OCRT_DEBUG");
            if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                fprintf(stderr, "error: --debug-fixed-bulk-phase-model requires OCRT_DEBUG=1\n");
                return 2;
            }
            if (!strcmp(m, "delta2bb") || !strcmp(m, "delta") || !strcmp(m, "bb2")) {
                cs.fixed_bulk_phase_model = 0;
            } else if (!strcmp(m, "hg-lut") || !strcmp(m, "hg_direct") || !strcmp(m, "hg")) {
                cs.fixed_bulk_phase_model = 1;
            } else if (!strcmp(m, "lut") || !strcmp(m, "table") || !strcmp(m, "p11-lut")) {
                cs.fixed_bulk_phase_model = 2;
            } else if (!strcmp(m, "lut-deltam") || !strcmp(m, "deltam") || !strcmp(m, "delta-m")) {
                cs.fixed_bulk_phase_model = 3;
            } else {
                fprintf(stderr, "error: unknown --debug-fixed-bulk-phase-model '%s' (expected delta2bb, hg-lut, lut, or lut-deltam)\n", m);
                return 2;
            }
        } else if (!strcmp(a, "--iop-phase-lut")) {
            fixed_bulk_cli_seen = 1;
            iop_phase_lut_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            if (cs.fixed_bulk_phase_model != 2) cs.fixed_bulk_phase_model = 3;  /* CAP FIX 2026-06-17 */
            cs.fixed_bulk_phase_lut_path = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--fixed-bulk-phase-lut")) {
            fprintf(stderr, "error: legacy option --fixed-bulk-phase-lut was removed; use --iop-phase-lut\n");
            return 2;
        } else if (!strcmp(a, "--iop-mie-phase")) {
            iop_mie_phase_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            /* #0: bulk-water particle VECTOR phase (P11/P12/P33) from .mie file.
             * Bulk IOP (a,b,bb) still comes from the IOP branch; only the phase
             * SHAPE (betal/gammal/alphal/zetal) comes from the .mie. */
            cs.water_mie_phase_path = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--water-mie-phase")) {
            fprintf(stderr, "error: ambiguous legacy option --water-mie-phase was removed; use --iop-mie-phase\n");
            return 2;
        } else if (!strcmp(a, "--ocrt-mie-moment-mode") || !strcmp(a, "--iop-mie-moment-mode")) {
            const int is_ocrt_control = !strcmp(a, "--ocrt-mie-moment-mode");
            if (is_ocrt_control) ocrt_mie_control_seen = 1;
            else                 iop_mie_control_seen = 1;
            water_option_family_mask |= is_ocrt_control
                                          ? WATER_OPTION_FAMILY_OCRT
                                          : WATER_OPTION_FAMILY_IOP;
            const char *m = next_str(argc, argv, &i, a);
            if (!strcmp(m, "dense") || !strcmp(m, "pchip-theta")) cs.water_mie_moment_mode = 0;
            else if (!strcmp(m, "gauss") || !strcmp(m, "quadrature") || !strcmp(m, "gauss-mu")) cs.water_mie_moment_mode = 1;
            else { fprintf(stderr, "error: %s must be dense|gauss (got '%s')\n", a, m); return 2; }
            if (cs.water_mie_moment_mode != 1 && !ocrt_advanced_on()) {
                fprintf(stderr,
                        "error: %s dense differs from the validated default gauss and requires OCRT_ADVANCED=1\n", a);
                return 2;
            }
        } else if (!strcmp(a, "--water-mie-moment-mode")) {
            fprintf(stderr, "error: ambiguous legacy option --water-mie-moment-mode was removed; use --ocrt-mie-moment-mode or --iop-mie-moment-mode\n");
            return 2;
        } else if (!strcmp(a, "--ocrt-mie-truncation") || !strcmp(a, "--iop-mie-truncation")) {
            const int is_ocrt_control = !strcmp(a, "--ocrt-mie-truncation");
            if (is_ocrt_control) {
                ocrt_mie_control_seen = 1;
                ocrt_mie_trunc_seen = 1;
            } else {
                iop_mie_control_seen = 1;
                iop_mie_trunc_seen = 1;
            }
            water_option_family_mask |= is_ocrt_control
                                          ? WATER_OPTION_FAMILY_OCRT
                                          : WATER_OPTION_FAMILY_IOP;
            /* ADVANCED BROAD-TRUNCATION POLICY
             *
             * Canonical FR631 production/reference physics uses the raw phase
             * with truncation OFF. The 0--0.005-degree local cap was negligible
             * in validation, whereas OSOAA-style broad truncation can transform
             * a large fraction of particle scattering and alter directional
             * Rrs/rrs by several percent. This option is therefore an explicit
             * compatibility/acceleration mode, never an automatic repair for a
             * coarse grid. Matched comparisons must also match physical depth,
             * b/omega/tau transforms, vector moments and source correction.
             *
             * DOC-REF:
             *   docs/OCRT_MIE_GRID_TRUNCATION_AND_VALIDATION_ARTIFACT_POLICY_2026-08-21.md */
            if (is_ocrt_control && !ocrt_advanced_on()) {
                fprintf(stderr,
                        "error: --ocrt-mie-truncation is an advanced physics option "
                        "and requires OCRT_ADVANCED=1; canonical FR631 default is OFF\n");
                return 2;
            }
            if (is_ocrt_control) {
                fprintf(stderr,
                        "WARNING [OCRT-TRUNCATION]: broad truncation enabled. "
                        "Canonical FR631 scientific reference is raw/OFF. Use this "
                        "advanced mode only with explicitly matched reference-model "
                        "truncation and physical-depth policies.\n");
            }
            cs.water_mie_truncation_mode = 1;
        } else if (!strcmp(a, "--water-mie-truncation")) {
            fprintf(stderr, "error: ambiguous legacy option --water-mie-truncation was removed; use --ocrt-mie-truncation or --iop-mie-truncation\n");
            return 2;
        } else if (!strcmp(a, "--ocrt-mie-ss-mode") || !strcmp(a, "--iop-mie-ss-mode")) {
            const int is_ocrt_control = !strcmp(a, "--ocrt-mie-ss-mode");
            if (is_ocrt_control) {
                ocrt_mie_control_seen = 1;
                ocrt_mie_ss_seen = 1;
            } else {
                iop_mie_control_seen = 1;
                iop_mie_ss_seen = 1;
            }
            water_option_family_mask |= is_ocrt_control
                                          ? WATER_OPTION_FAMILY_OCRT
                                          : WATER_OPTION_FAMILY_IOP;
            cs.water_mie_ss_mode = next_int(argc, argv, &i, a);
            if (cs.water_mie_ss_mode < 0 || cs.water_mie_ss_mode > 2) {
                fprintf(stderr, "error: %s must be 0|1|2 (2 = NT-TMS)\n", a);
                return 2;
            }
        } else if (!strcmp(a, "--water-mie-ss-mode")) {
            fprintf(stderr, "error: ambiguous legacy option --water-mie-ss-mode was removed; use --ocrt-mie-ss-mode or --iop-mie-ss-mode\n");
            return 2;
        } else if (!strcmp(a, "--ocrt-mie-moment-nmu") || !strcmp(a, "--iop-mie-moment-nmu")) {
            const int is_ocrt_control = !strcmp(a, "--ocrt-mie-moment-nmu");
            if (is_ocrt_control) ocrt_mie_control_seen = 1;
            else                 iop_mie_control_seen = 1;
            water_option_family_mask |= is_ocrt_control
                                          ? WATER_OPTION_FAMILY_OCRT
                                          : WATER_OPTION_FAMILY_IOP;
            cs.water_mie_moment_n_mu = next_int(argc, argv, &i, a);
            if (cs.water_mie_moment_n_mu <= 1) {
                fprintf(stderr, "error: %s must be greater than 1\n", a);
                return 2;
            }
            if (cs.water_mie_moment_n_mu != 400 && !ocrt_advanced_on()) {
                fprintf(stderr,
                        "error: %s %d differs from the validated default 400 and requires OCRT_ADVANCED=1\n",
                        a, cs.water_mie_moment_n_mu);
                return 2;
            }
        } else if (!strcmp(a, "--water-mie-moment-nmu")) {
            fprintf(stderr, "error: ambiguous legacy option --water-mie-moment-nmu was removed; use --ocrt-mie-moment-nmu or --iop-mie-moment-nmu\n");
            return 2;
        } else if (!strcmp(a, "--iop-phase-nphi")) {
            fixed_bulk_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            /* ADVANCED (gated behind OCRT_ADVANCED=1): in-water phase-kernel azimuth
             * quadrature.  Default 720 is near-converged (~0.06% vs 2880). */
            const char *adv = getenv("OCRT_ADVANCED");
            if (!adv || !*adv || !strcmp(adv, "0")) {
                fprintf(stderr, "error: --iop-phase-nphi requires OCRT_ADVANCED=1 "
                                "(default 720 is validated)\n");
                return 2;
            }
            cs.fixed_bulk_phase_nphi = next_int(argc, argv, &i, a);
        } else if (!strcmp(a, "--fixed-bulk-phase-nphi")) {
            fprintf(stderr, "error: legacy option --fixed-bulk-phase-nphi was removed; use --iop-phase-nphi\n");
            return 2;
        } else if (!strcmp(a, "--scalar-ff-iop") || !strcmp(a, "--ff-iop") || !strcmp(a, "--ff-scalar-iop")) {
            fixed_bulk_cli_seen = 1;
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            if (!ocrt_debug_on()) {
                fprintf(stderr, "error: %s is a diagnostic option and requires OCRT_DEBUG=1\n", "--scalar-ff-iop");
                return 2;
            }
            cs.fixed_bulk_phase_model = 4;
            cs.scalar_ff_iop_path = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--scalar-ff-bbfrac") || !strcmp(a, "--scalar-ff-bb-over-b") || !strcmp(a, "--ff-bbfrac")) {
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            cs.fixed_bulk_phase_model = 4;
            cs.scalar_ff_bb_over_b = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--scalar-ff-n") || !strcmp(a, "--ff-n")) {
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            cs.fixed_bulk_phase_model = 4;
            cs.scalar_ff_refractive_index = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--scalar-ff-mu") || !strcmp(a, "--ff-mu")) {
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            cs.fixed_bulk_phase_model = 4;
            cs.scalar_ff_mu = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--scalar-ff-lmax") || !strcmp(a, "--ff-lmax")) {
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            cs.fixed_bulk_phase_model = 4;
            cs.fixed_bulk_lmax = next_int(argc, argv, &i, a);
        } else if (!strcmp(a, "--scalar-ff-phase-nphi") || !strcmp(a, "--scalar-ff-nphi") || !strcmp(a, "--ff-phase-nphi") || !strcmp(a, "--ff-nphi")) {
            /* ADVANCED: diagnostic alias of --iop-phase-nphi. */
            const char *adv = getenv("OCRT_ADVANCED");
            if (!adv || !*adv || !strcmp(adv, "0")) {
                fprintf(stderr, "error: %s requires OCRT_ADVANCED=1 (default 720 is validated)\n", a);
                return 2;
            }
            water_option_family_mask |= WATER_OPTION_FAMILY_IOP;
            cs.fixed_bulk_phase_model = 4;
            cs.fixed_bulk_phase_nphi = next_int(argc, argv, &i, a);
        } else if (!strcmp(a, "--n-mu-water")) {
            /* ADVANCED (gated behind OCRT_ADVANCED=1): in-water angular quadrature
             * (positive Gauss nodes).  Default 48 is validated; lowering it under-
             * converges rrs (n_mu=24 +3%, n_mu=32 +1.5%, n_mu=48 +0.2% vs n_mu=64,
             * Brown_earth Csed5/555nm).  >= 48 recommended for forward-peaked / high-
             * omega phases and LUT mode; floor rises with single-scattering albedo. */
            const char *adv = getenv("OCRT_ADVANCED");
            if (!adv || !*adv || !strcmp(adv, "0")) {
                fprintf(stderr, "error: --n-mu-water requires OCRT_ADVANCED=1 "
                                "(default 48 is validated; changing in-water quadrature can degrade accuracy)\n");
                return 2;
            }
            cs.n_mu_water_override = next_int(argc, argv, &i, a);
        } else if (!strcmp(a, "--n-mu-water-sky")) {
            /* coarse in-water grid for atmosphere sky-light coupling beams only
             * (direct beam keeps full n_mu_water). Speeds up coupled phyto runs. */
            cs.n_mu_water_sky = next_int(argc, argv, &i, a);
        } else if (!strcmp(a, "--water-m-max")) {
            cs.water_m_max_override = next_int(argc, argv, &i, a);
            if (!ocrt_advanced_on()) {
                fprintf(stderr, "error: --water-m-max overrides the auto water Fourier cap "
                        "(convergence knob; requires OCRT_ADVANCED=1)\n");
                return 2;
            }
        } else if (!strcmp(a, "--water-shared-grid")) {
            cs.water_shared_grid = 1;
        } else if (!strcmp(a, "--batch-full-grid")) {
            batch_full_grid_in = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--water-view-as-node")) {
            cs.water_view_as_node = 1;
        } else if (!strcmp(a, "--water-view-as-node-off")) {
            /* v1.09 commit #16 diagnostics: force continuous-interpolation view
             * extraction (the full-grid mode) in a standalone run, so grid rows
             * can be verified against independent single-case processes. */
            cs.water_view_as_node = 0;
        } else if (!strcmp(a, "--sigma-model")) {
            const char *v = next_str(argc, argv, &i, a);
            if (!strcmp(v, "ocrt-floor") || !strcmp(v, "ocrt_floor"))
                cs.sigma_type = 1;
            else if (!strcmp(v, "nakajima-tanaka") || !strcmp(v, "nakajima_tanaka"))
                cs.sigma_type = 0;
            else {
                fprintf(stderr, "error: --sigma-model must be ocrt-floor|nakajima-tanaka\n");
                return 2;
            }
            if (cs.sigma_type != 1 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: non-default --sigma-model requires OCRT_ADVANCED=1\n");
                return 2;
            }
        } else if (!strcmp(a, "--sigma-type")) {
            /* Legacy numeric selector. Prefer --sigma-model.  ADVANCED 게이트.
             * 경사분산 모델: 0=Nakajima-Tanaka, 1=OCRT floor law. */
            cs.sigma_type = next_int(argc, argv, &i, a);
            if (cs.sigma_type != 1 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: --sigma-type %d differs from the default 1 "
                        "(OCRT floor slope law; selector requires OCRT_ADVANCED=1)\n",
                        cs.sigma_type);
                return 2;
            }
        } else if (!strcmp(a, "--q-convention")) {
            /* 2026-07-15 Jae(G4): ADVANCED 게이트.  Stokes Q 부호 규약
             * (0=legacy (Rs-Rp)/2, 1=Hansen/Mishchenko=기본). */
            cs.q_convention = next_int(argc, argv, &i, a);
            if (cs.q_convention != 1 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: --q-convention %d differs from the default 1 "
                        "(Hansen/Mishchenko; convention selector requires OCRT_ADVANCED=1)\n",
                        cs.q_convention);
                return 2;
            }
        } else if (!strcmp(a, "--pssa")) {
            /* v1.11 (2026-09-05): --pssa is a plain on/off switch.  IPSS
             * (Zhai & Hu 2022) is the only spherical algorithm in the tree;
             * the legacy average-secant Chapman PSSA and its --pssa-mode
             * selector were REMOVED, so no run can land on the old path. */
            opts.pssa = 1;
        } else if (!strcmp(a, "--pssa-mode")) {
            fprintf(stderr,
                "error: --pssa-mode was removed in v1.11.  IPSS "
                "(Zhai & Hu 2022) is now the only spherical algorithm;\n"
                "       the legacy average-secant Chapman PSSA is gone from "
                "the source tree.  Use plain --pssa.\n");
            return 2;
        } else if (!strcmp(a, "--decouple-sunglint")) {
            /* 2026-07-15 Jae(G4): 반전.  기본은 선글린트 포함(decouple=0)이며,
             * 이 플래그 지정 시에만 직달 선글린트 단일반사 항을 표적 출력에서
             * 사후 분리·제외한다(AF1982 기준 정합값).  구 --no-decouple-sunglint
             * 삭제.  중요·빈용 옵션이므로 게이트 없음(Jae). */
            cs.decouple_sunglint = 1;
        /* 2026-07-15 Jae(G3): --no-aerosol 삭제 -- 빈 분기(무동작)였다.
         * off 는 --aod-555 0(기본)이 담당한다. */
        } else if (!strcmp(a, "--mie")) {
            mie_path = next_str(argc, argv, &i, a);
            /* v1.09 commit #3: M50C is EXCLUDED from official validation
             * (diagnostic-only) per M50C_POLICY_v0_9 — lineage contamination
             * (83- vs 361-angle files are physically different models;
             * .inp fraction ambiguity; P11-only legacy table).  Official
             * aerosol subset: T50 + C50 + M80C.  Reinstatement requires the
             * 5-step lineage-clean procedure (see inputs/deprecated/
             * DEPRECATED_M50C.md).  Case-insensitive filename match;
             * renaming the file bypasses this gate — the gate is a guard,
             * not the lineage proof. */
            for (const char *p = mie_path; *p; ++p) {
                if ((p[0]=='M'||p[0]=='m') && p[1]=='5' && p[2]=='0' &&
                    (p[3]=='C'||p[3]=='c')) {
                    const char *dbg = getenv("OCRT_DEBUG");
                    if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                        fprintf(stderr, "error: M50C aerosol model is diagnostic-only "
                                "(lineage contamination, M50C_POLICY_v0_9; official subset "
                                "T50/C50/M80C). Requires OCRT_DEBUG=1; see "
                                "inputs/deprecated/DEPRECATED_M50C.md\n");
                        return 2;
                    }
                    fprintf(stderr, "# WARNING: M50C is diagnostic-only (lineage "
                            "contamination) — results are NOT valid for official "
                            "validation or production\n");
                    break;
                }
            }
        } else if (!strcmp(a, "--aod-555")) {
            user_aod = next_double(argc, argv, &i, a);
            user_aod_ref_nm = 555.0;
        } else if (!strcmp(a, "--aod-865")) {
            user_aod = next_double(argc, argv, &i, a);
            user_aod_ref_nm = 865.0;
        } else if (!strcmp(a, "--aod-ref-wavelength")) {
            user_aod_ref_nm = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--aod")) {
            /* Compatibility alias: AOD at the configured reference wavelength
             * (555 nm by default), never at the current calculation band. */
            user_aod = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--aod-band")) {
            fprintf(stderr, "error: --aod-band was removed because band-specific AOD breaks spectral consistency; use --aod-555 or --aod-865\n");
            return 2;
        } else if (!strcmp(a, "--trunc-aer-fit")) {
            /* v1.09 commit #4: DEBUG-gated.  The delta-fit/delta-M aerosol
             * truncation family is a NON-VALIDATED research path (production
             * and the validated OSOAA protocol use --trunc-aer-loglin +
             * --aer-phase-kernel value, MAPE 0.48%).  Even with the
             * corrected --nt-tau similarity transform (see rt_aerosol_
             * runtime.c) this family retains ~+-10% truncation-shape
             * residuals vs OSOAA (BUG_SUSPECTS S-001 experiment, items 6-9).
             * Kept for literature-standard delta-M A/B research only. */
            if (!ocrt_debug_on()) {
                fprintf(stderr, "error: --trunc-aer-fit requires OCRT_DEBUG=1 "
                        "(non-validated research truncation; production path is "
                        "--trunc-aer-loglin; see BUG_SUSPECTS S-001)\n");
                return 2;
            }
            theta_cut_deg = next_double(argc, argv, &i, a);
            have_trunc_aer_fit = 1;
        } else if (!strcmp(a, "--trunc-aer-m")) {
            if (!ocrt_debug_on()) {   /* v1.09 commit #4, same rationale */
                fprintf(stderr, "error: --trunc-aer-m requires OCRT_DEBUG=1 "
                        "(non-validated research truncation; production path is "
                        "--trunc-aer-loglin; see BUG_SUSPECTS S-001)\n");
                return 2;
            }
            delta_m_N = next_int(argc, argv, &i, a);
            have_trunc_aer_m = 1;
        /* 2026-07-15 Jae(G3): --trunc-aer-loglin 삭제.  loglin 은 에어로졸 on
         * 시 자동 기본(아래 dispatch 확장 참조; δ 계열 DEBUG 플래그가 오버라이드). */
        /* 2026-07-15 Jae(G3): --aer-phase-kernel 삭제.  에어로졸 on 시 value
         * 커널이 자동 기본(아래 dispatch)이며, moment 는 CLI 도달 불가로 전환. */
        } else if (!strcmp(a, "--aer-h-km")) {
            /* OCRT aerosol vertical profile: n(z) proportional to exp(-z/H).
             * H=2 km is the production model.  A positive advanced override
             * is retained for controlled sensitivity studies; legacy H=0
             * reference-code dispatch has been removed. */
            opts.aer_h_km = next_double(argc, argv, &i, a);
            if (!(opts.aer_h_km > 0.0)) {
                fprintf(stderr, "error: --aer-h-km must be positive (OCRT default: 2.0 km exponential profile)\n");
                return 2;
            }
            have_aer_h = 1;
        } else if (!strcmp(a, "--output-advanced")) {
            /* 2026-07-15 Jae 지시: 구 --output-mode {simple|debug} 를 대체·삭제.
             * 기본(플래그 없음) = 간단 열(I/Q/U 반사도만).
             * --output-advanced = --batch 결과 CSV 에 상세 열(기준 대비 오차,
             * 수렴 차수, 계산 시간) 추가.  "debug" 명칭이 부적절하다는 판단으로
             * 교체.  내부 필드명 output_mode_debug 는 rt_io.c 작성기와의 호환을
             * 위해 유지(동작 동일). */
            opts.output_mode_debug = 1;
        /* 2026-07-15 Jae 지시(1안): --lut 삭제.  M2 에서 격자가 기본 모드이고
         * lut_mode 는 --output-full-grid/기본 모드 판정에서 설정된다. */
        } else if (!strcmp(a, "--output-full-grid")) {
            /* 2026-07-15 Jae 지시: 완전 동일 동작이던 --lut-output 별칭 삭제. */
            lut_output = next_str(argc, argv, &i, a);
            lut_mode = 1;
        } else if (!strcmp(a, "--lut-vza-step") || !strcmp(a, "--vza-gap")) {
            lut_vza_step = next_double(argc, argv, &i, a); have_vza_step = 1;
        } else if (!strcmp(a, "--lut-raa-step") || !strcmp(a, "--raa-gap")) {
            lut_raa_step = next_double(argc, argv, &i, a); have_raa_step = 1;
        } else if (!strcmp(a, "--lut-vza-max")) {
            lut_vza_max = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--atm-profile") || !strcmp(a, "--gas-profile")) {
            const char *s = next_str(argc, argv, &i, a);
            if      (!strcmp(s, "usstd76"))  atm_profile = RT_AFGL_USSTD76;
            else if (!strcmp(s, "tropical")) atm_profile = RT_AFGL_TROPICAL;
            else if (!strcmp(s, "mlsumm"))   atm_profile = RT_AFGL_MLSUMM;
            else if (!strcmp(s, "mlwint"))   atm_profile = RT_AFGL_MLWINT;
            else if (!strcmp(s, "sasumm"))   atm_profile = RT_AFGL_SASUMM;
            else if (!strcmp(s, "sawint"))   atm_profile = RT_AFGL_SAWINT;
            else if (!strcmp(s, "userdef"))  atm_profile = RT_AFGL_USERDEF;
            else { fprintf(stderr, "error: --atm-profile must be one of usstd76|tropical|mlsumm|mlwint|sasumm|sawint|userdef\n"); return 2; }
        /* 2026-07-15 Jae 지시(G2): --afgl-dir 삭제 -- --xsec-dir 와 동일
         * 원칙(사용자 불변 경로).  AFGL 디렉터리는 inputs/afgl_atm 고정.
         * userdef 프로파일은 inputs/afgl_atm/afgl_userdef.dat 에 두어야 한다. */
        /* 2026-07-15 Jae 지시(G2): --xsec-dir 삭제 -- 사용자 불변 경로
         * (inputs/xsec 고정). */
        } else if (!strcmp(a, "--gas-column-h2o") || !strcmp(a, "--pwv-g-cm")) {
            /* H2O column in g/cm² → convert to molecules/cm² (M_h2o=18.015 g/mol)
             * --pwv-g-cm is the alias matching user-facing PWV terminology. */
            double v = next_double(argc, argv, &i, a);
            gas_col_override[RT_GAS_H2O] = v * 6.02214076e23 / 18.01528;
        } else if (!strcmp(a, "--gas-column-o3")) {
            /* O3 column in DU (Dobson Units). 1 DU = 2.6868e16 mol/cm² */
            double v = next_double(argc, argv, &i, a);
            gas_col_override[RT_GAS_O3] = v * 2.6868e16;
        } else if (!strcmp(a, "--gas-column-o2")) {
            /* O2 column in mol/cm² (well-mixed ~4.5e24 default) */
            gas_col_override[RT_GAS_O2] = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--gas-column-co2")) {
            /* CO2 column in mol/cm².  usstd76 적분 기본 7.116e21 (~400 ppmv;
             * 2026-07-15 실측) -- 구 주석의 8.5e21 은 오기였다. */
            gas_col_override[RT_GAS_CO2] = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--gas-column-no2")) {
            /* NO2 column in DU (1 DU = 2.6868e16 mol/cm²) */
            double v = next_double(argc, argv, &i, a);
            gas_col_override[RT_GAS_NO2] = v * 2.6868e16;
        } else if (!strcmp(a, "--gas-column-ch4")) {
            /* CH4 column in mol/cm² */
            gas_col_override[RT_GAS_CH4] = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--nt-tau")) {
            if (!ocrt_debug_on()) {   /* v1.09 commit #4: pairs with the
                 * DEBUG-gated delta-fit/delta-M family above. */
                fprintf(stderr, "error: --nt-tau requires OCRT_DEBUG=1 "
                        "(companion of the research truncation family; "
                        "formula corrected in v1.09 commit #4)\n");
                return 2;
            }
            apply_nt_tau = 1;
        } else if (!strcmp(a, "--aer-l-max")) {
            aer_L_max = next_int(argc, argv, &i, a);
            have_aer_L_max = 1;
        } else if (!strcmp(a, "--no-rayleigh")) {
            cs.pressure_hpa = 0.0;
        } else if (!strcmp(a, "--rayleigh")) {
            const char *v = next_str(argc, argv, &i, a);
            if (!strcmp(v, "off")) cs.pressure_hpa = 0.0;
            else if (!strcmp(v, "on")) { if (cs.pressure_hpa <= 0.0) cs.pressure_hpa = 1013.25; }
            else { fprintf(stderr, "error: --rayleigh must be on|off\n"); return 2; }
        } else if (!strcmp(a, "--pressure") || !strcmp(a, "--surface-pressure")) {
            /* Surface pressure (hPa).  Scales Rayleigh optical depth (P/1013.25);
             * pressure 0 => no Rayleigh.  Replaces the removed --no-rayleigh /
             * --with-rayleigh toggles: atmosphere is now controlled purely by
             * physical quantities (pressure for Rayleigh, AOD for aerosol). */
            cs.pressure_hpa = next_double(argc, argv, &i, a);
        } else if (!strcmp(a, "--n-mu")) {
            /* M2 (2026-07-14, Jae 지시): 수치 노브는 ADVANCED 게이트. 기본값 재명시는 통과. */
            opts.n_mu = next_int(argc, argv, &i, a);
            if (opts.n_mu != 24 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: --n-mu %d differs from the validated default 24 "
                        "(convergence knob; requires OCRT_ADVANCED=1)\n", opts.n_mu);
                return 2;
            }
        } else if (!strcmp(a, "--n-layers")) {
            /* 비에어로졸 기본 40 게이트는 여기서, 에어로졸(400)은 기존
             * v1.09 #13 게이트(파싱 후 블록)에서 검사한다. */
            opts.n_layers = next_int(argc, argv, &i, a);
            have_n_layers = 1;
        } else if (!strcmp(a, "--max-orders")) {
            opts.max_orders = next_int(argc, argv, &i, a);
            if (opts.max_orders != 20 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: --max-orders %d differs from the validated default 20 "
                        "(convergence knob; requires OCRT_ADVANCED=1)\n", opts.max_orders);
                return 2;
            }
        } else if (!strcmp(a, "--l-max")) {
            const char *s = next_str(argc, argv, &i, a);
            opts.l_max = (!strcmp(s, "auto")) ? -1 : (int)strtol(s, NULL, 10);
            if (opts.l_max != -1 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: --l-max %s overrides the auto dispatch "
                        "(convergence knob; requires OCRT_ADVANCED=1)\n", s);
                return 2;
            }
        /* Rayleigh optical depth is fixed to the production Bodhaine model. */
        } else if (!strcmp(a, "--integration-method")) {
            const char *s = next_str(argc, argv, &i, a);
            if      (!strcmp(s, "linear"))   opts.integration_method = RT_INTEGRATION_METHOD_LINEAR;
            else if (!strcmp(s, "constant")) {
                /* DEBUG: 'constant' reproduces v1's grid-divergence pathology (paper/demo only). */
                const char *dbg = getenv("OCRT_DEBUG");
                if (!dbg || !*dbg || !strcmp(dbg, "0")) {
                    fprintf(stderr, "error: --integration-method constant requires OCRT_DEBUG=1 "
                                    "(reproduces v1 grid-divergence pathology; production is linear)\n");
                    return 2;
                }
                opts.integration_method = RT_INTEGRATION_METHOD_CONSTANT;
            }
            else { fprintf(stderr, "error: unknown integration-method '%s' (use linear|constant)\n", s); return 2; }
        } else if (!strcmp(a, "--sos-max-orders")) {
            opts.sos.max_iterations = next_int(argc, argv, &i, a);
            if (opts.sos.max_iterations != 20 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: --sos-max-orders %d differs from the validated default 20 "
                        "(convergence knob; requires OCRT_ADVANCED=1)\n", opts.sos.max_iterations);
                return 2;
            }
        } else if (!strcmp(a, "--sos-tolerance")) {
            opts.sos.tolerance = next_double(argc, argv, &i, a);
            if (opts.sos.tolerance != 1e-7 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: --sos-tolerance %g differs from the validated default 1e-7 "
                        "(convergence knob; requires OCRT_ADVANCED=1)\n", opts.sos.tolerance);
                return 2;
            }
        } else if (!strcmp(a, "--sos-acceleration")) {
            const char *s = next_str(argc, argv, &i, a);
            if      (!strcmp(s, "plain"))     opts.sos.acceleration = RT_SOS_ACCELERATION_PLAIN;
            else if (!strcmp(s, "geometric")) {
                /* v1.09 commit #4 (BUG_SUSPECTS S-003): the geometric value
                 * was ACCEPTED but never implemented anywhere — it silently
                 * behaved as plain.  Rejected now instead of lying. */
                fprintf(stderr, "error: --sos-acceleration geometric is NOT implemented "
                        "(was a silent no-op; S-003). Use plain.\n");
                return 2;
            }
            else { fprintf(stderr, "error: unknown sos-acceleration '%s' (use plain)\n", s); return 2; }
        } else if (!strcmp(a, "--sos-save-orders")) {
            if (!ocrt_debug_on()) {
                fprintf(stderr, "error: %s is a diagnostic option and requires OCRT_DEBUG=1\n", "--sos-save-orders");
                return 2;
            }
            opts.sos.save_orders = 1;
        } else if (!strcmp(a, "--m-max")) {
            opts.fourier_m_max = next_int(argc, argv, &i, a);
            have_m_max = 1;
        /* 2026-07-15 Jae 지시(G2): --tau-r / --tau-r-from-input 삭제.
         * τ_R 은 --pressure 로부터 내부 산출(Bodhaine 1999)만 지원한다.
         * 배치 CSV 의 tau_R_ref 열은 참고 출력용으로 잔존(rt_io.c). */
        /* 2026-07-15 Jae 지시: --view-as-node 삭제.  단일 기하 모드가
         * view_as_node 를 강제 활성하므로(아래 single_geo_req 처리) 순수
         * 선언용이었다.  사용자 혼란 방지 목적. */
        /* 2026-07-15 Jae 지시(1안): --vector 삭제.  vector_mode 는 기본 1
         * (rt_types.h)이며 스칼라 경로는 M1 에서 이미 제거되었다.  패키지 내부
         * 스크립트 14+곳의 --vector 토큰을 같은 커밋에서 일괄 제거함. */
        } else if (!strcmp(a, "--conv-tol")) {
            opts.conv_tol = next_double(argc, argv, &i, a);
            if (opts.conv_tol != 1e-6 && !ocrt_advanced_on()) {
                fprintf(stderr, "error: --conv-tol %g differs from the validated default 1e-6 "
                        "(convergence knob; requires OCRT_ADVANCED=1)\n", opts.conv_tol);
                return 2;
            }
        } else if (!strcmp(a, "--verbose")) {
            opts.verbose = 1;
        } else if (!strcmp(a, "--simple-compare") || !strcmp(a, "--ccrr-compare")) {
            if (!ocrt_debug_on()) {
                fprintf(stderr, "error: %s is a diagnostic option and requires OCRT_DEBUG=1\n", "--simple-compare");
                return 2;
            }
            ccrr_compare_in = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--simple-output") || !strcmp(a, "--ccrr-output")) {
            ccrr_compare_out = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--batch")) {
            batch_in = next_str(argc, argv, &i, a);
        } else if (!strcmp(a, "--output")) {
            batch_out = next_str(argc, argv, &i, a);
        } else {
            fprintf(stderr, "error: unknown argument '%s' (try --help)\n", a);
            return 2;
        }
    }

    /* Public water-input contract (OCRT v1.1 release).
     *
     * For --surface ocean, --water-model ocrt|ccrr|iop is mandatory and may
     * appear exactly once.  Every branch-specific option carries the matching
     * prefix.  The parser records provenance first; this block lowers that
     * provenance to the legacy execution flags consumed by the existing OCRT
     * underwater solver. */
    if (ccrr_compare_in) {
        if (water_model_seen || water_option_family_mask != 0u) {
            fprintf(stderr,
                    "error: --ccrr-compare cannot be combined with a single-run --water-model branch\n");
            return 2;
        }
    } else {
        if (cs.surface == RT_SURFACE_OCEAN) {
            if (!water_model_seen) {
                fprintf(stderr,
                        "error: --surface ocean requires exactly one --water-model selection (ocrt | ccrr | iop)\n");
                return 2;
            }
        } else if (water_model_seen || water_option_family_mask != 0u) {
            fprintf(stderr,
                    "error: --water-model and its branch-specific options require --surface ocean\n");
            return 2;
        }
    }

    if (water_model_seen) {
        const rt_water_input_mode_t mode = (rt_water_input_mode_t)cs.water_input_mode;
        const unsigned expected_family = water_option_family_for_mode(mode);
        const unsigned wrong_family = water_option_family_mask & ~expected_family;
        if (wrong_family != 0u) {
            fprintf(stderr,
                    "error: --water-model %s cannot be combined with options from another water branch (seen-family mask=0x%x)\n",
                    water_input_mode_name(mode), water_option_family_mask);
            return 2;
        }

        /* Reset derived execution flags before lowering the selected branch. */
        cs.ccrr_mode = 0;
        cs.fixed_bulk_iop_mode = 0;

        if (mode == RT_WATER_INPUT_OCRT || mode == RT_WATER_INPUT_CCRR) {
            const char *prefix = (mode == RT_WATER_INPUT_OCRT) ? "--ocrt-" : "--ccrr-";
            if (!chl_input_seen || !tsm_input_seen || !adom_input_seen) {
                fprintf(stderr,
                        "error: --water-model %s requires explicit %schl, %stsm and %sadom440 values (zero is allowed); %sadom-slope is optional\n",
                        water_input_mode_name(mode), prefix, prefix, prefix, prefix);
                return 2;
            }
            if (fixed_bulk_cli_seen || iop_a_input_seen || iop_b_input_seen ||
                iop_bb_input_seen || cs.scalar_ff_iop_path ||
                cs.fixed_bulk_phase_lut_path ||
                (cs.water_mie_phase_path && cs.water_mie_phase_path[0])) {
                fprintf(stderr,
                        "error: --water-model %s cannot use direct-IOP inputs or IOP phase sources\n",
                        water_input_mode_name(mode));
                return 2;
            }
            if (cs.ccrr_chl_mg_m3 < 0.0 || cs.ccrr_min_g_m3 < 0.0 ||
                cs.a_cdom_440_m_inv < 0.0 || cs.detritus_a440_m_inv < 0.0) {
                fprintf(stderr,
                        "error: Chl, TSM, aDOM440 and detritus absorption must be non-negative\n");
                return 2;
            }
            if (simple_adom_s_explicit && !(cs.S_cdom_nm_inv > 0.0)) {
                fprintf(stderr, "error: %sadom-slope must be positive\n", prefix);
                return 2;
            }

            cs.cdom_ref_lambda_nm = 440.0;
            if (!(cs.S_cdom_nm_inv > 0.0)) cs.S_cdom_nm_inv = 0.014;
            cs.water_constituent_model = (mode == RT_WATER_INPUT_OCRT)
                                           ? RT_WATER_CONSTITUENT_OCRT
                                           : RT_WATER_CONSTITUENT_CCRR;

            water_pure_from_zero =
                (cs.ccrr_chl_mg_m3 == 0.0 &&
                 cs.ccrr_min_g_m3 == 0.0 &&
                 cs.a_cdom_440_m_inv == 0.0);

            if (water_pure_from_zero) {
                /* Execute the native pure-water path, not a zero-valued
                 * constituent adapter.  This guarantees OCRT-zero and
                 * CCRR-zero requests share the exact same solver path. */
                cs.ccrr_mode = 0;
                cs.F_sun = M_PI;
                cs.water_mie_phase_path = NULL;
                cs.ccrr_phase_moments_path = NULL;
                cs.ccrr_particle_phase_lut_path = NULL;
                fprintf(stderr,
                        "# water branch: %s selected with Chl=TSM=aDOM440=0; using native pure-water path\n",
                        water_input_mode_name(mode));
            } else {
                cs.ccrr_mode = 1;
                cs.F_sun = 1.0;  /* frozen constituent-interface normalization */
            }

            if (mode == RT_WATER_INPUT_OCRT) {
                if (tsm_species_cli_seen && !(cs.ccrr_min_g_m3 > 0.0)) {
                    fprintf(stderr, "error: --ocrt-tsm-species requires --ocrt-tsm > 0\n");
                    return 2;
                }
                if ((organic_group_cli_seen || detritus_cli_seen) &&
                    !(cs.ccrr_chl_mg_m3 > 0.0)) {
                    fprintf(stderr,
                            "error: --ocrt-phyto-group/--ocrt-detritus-* require --ocrt-chl > 0\n");
                    return 2;
                }
                if (!(cs.detritus_slope_nm_inv > 0.0)) {
                    fprintf(stderr, "error: --ocrt-detritus-slope must be positive\n");
                    return 2;
                }
                if (cs.detritus_slope_nm_inv < ORGANIC_DETRITUS_SLOPE_MIN ||
                    cs.detritus_slope_nm_inv > ORGANIC_DETRITUS_SLOPE_MAX) {
                    fprintf(stderr,
                            "warning: OCRT detritus slope %.8g is outside the source range [%.4g,%.4g] nm^-1\n",
                            cs.detritus_slope_nm_inv,
                            ORGANIC_DETRITUS_SLOPE_MIN, ORGANIC_DETRITUS_SLOPE_MAX);
                }
                if (cs.ccrr_chl_mg_m3 > 0.0 && have_wl &&
                    !rt_iop_organic_wavelength_supported(cs.wavelength_nm)) {
                    fprintf(stderr,
                            "error: OCRT Chl model is validated only over %.0f-%.0f nm (requested %.10g nm)\n",
                            ORGANIC_WAVELENGTH_MIN_NM, ORGANIC_WAVELENGTH_MAX_NM,
                            cs.wavelength_nm);
                    return 2;
                }
                if (cs.ccrr_phase_moments_path || cs.ccrr_particle_phase_lut_path) {
                    fprintf(stderr,
                            "error: CCRR phase options cannot be used with --water-model ocrt\n");
                    return 2;
                }
                if (ocrt_mie_control_seen &&
                    !(cs.ccrr_chl_mg_m3 > 0.0 || cs.ccrr_min_g_m3 > 0.0)) {
                    fprintf(stderr,
                            "error: --ocrt-mie-* controls require --ocrt-chl > 0 or --ocrt-tsm > 0\n");
                    return 2;
                }
                if (ocrt_mie_ss_seen && cs.water_mie_ss_mode != 0 &&
                    !ocrt_mie_trunc_seen) {
                    fprintf(stderr,
                            "error: nonzero --ocrt-mie-ss-mode requires --ocrt-mie-truncation\n");
                    return 2;
                }
            } else {
                if (organic_group_cli_seen || detritus_cli_seen || tsm_species_cli_seen) {
                    fprintf(stderr,
                            "error: OCRT-prefixed particle options cannot be used with --water-model ccrr\n");
                    return 2;
                }
                if (ccrr_particle_selector_seen && !ccrr_particle_lut_seen) {
                    fprintf(stderr,
                            "error: --ccrr-particle-phase-case-id/-lmax/-nphi require --ccrr-particle-phase-lut\n");
                    return 2;
                }
                if (cs.ccrr_phase_moments_path && cs.ccrr_particle_phase_lut_path) {
                    fprintf(stderr,
                            "error: --ccrr-phase-moments and --ccrr-particle-phase-lut are mutually exclusive\n");
                    return 2;
                }
                if ((cs.ccrr_phase_moments_path || cs.ccrr_particle_phase_lut_path) &&
                    !(cs.ccrr_chl_mg_m3 > 0.0 || cs.ccrr_min_g_m3 > 0.0)) {
                    fprintf(stderr,
                            "error: CCRR particle-phase options require --ccrr-chl > 0 or --ccrr-tsm > 0\n");
                    return 2;
                }
                if (cs.ccrr_particle_phase_lmax <= 0 || cs.ccrr_particle_phase_nphi <= 0) {
                    fprintf(stderr, "error: CCRR particle phase lmax/nphi must be positive\n");
                    return 2;
                }
            }
        } else if (mode == RT_WATER_INPUT_IOP) {
            if (chl_input_seen || tsm_input_seen || adom_input_seen ||
                organic_group_cli_seen || detritus_cli_seen ||
                tsm_species_cli_seen || cs.ccrr_phase_moments_path ||
                cs.ccrr_particle_phase_lut_path) {
                fprintf(stderr,
                        "error: --water-model iop cannot use OCRT/CCRR constituent options\n");
                return 2;
            }
            if (!cs.scalar_ff_iop_path &&
                (!iop_a_input_seen || !iop_b_input_seen || !iop_bb_input_seen)) {
                fprintf(stderr,
                        "error: --water-model iop requires explicit --iop-a, --iop-b and --iop-bb values\n");
                return 2;
            }
            if (!cs.scalar_ff_iop_path &&
                (!(cs.fixed_a_total_m_inv >= 0.0) ||
                 !(cs.fixed_b_total_m_inv > 0.0) ||
                 !(cs.fixed_bb_total_m_inv > 0.0) ||
                 cs.fixed_bb_total_m_inv >= 0.5 * cs.fixed_b_total_m_inv)) {
                fprintf(stderr,
                        "error: IOP input requires a>=0, b>0 and 0<bb<0.5*b\n");
                return 2;
            }
            {
                const int n_phase_sources = iop_phase_lut_seen +
                                            iop_mie_phase_seen;
                if (n_phase_sources > 1) {
                    fprintf(stderr,
                            "error: select at most one IOP phase source: --iop-phase-lut or --iop-mie-phase\n");
                    return 2;
                }
            }
            if (iop_mie_control_seen && !iop_mie_phase_seen) {
                fprintf(stderr,
                        "error: --iop-mie-* controls require --iop-mie-phase\n");
                return 2;
            }
            if (iop_mie_ss_seen && cs.water_mie_ss_mode != 0 &&
                !iop_mie_trunc_seen) {
                fprintf(stderr,
                        "error: nonzero --iop-mie-ss-mode requires --iop-mie-truncation\n");
                return 2;
            }
            cs.fixed_bulk_iop_mode = 1;
            cs.F_sun = 1.0;
        } else {
            fprintf(stderr, "error: internal invalid water input mode %d\n", cs.water_input_mode);
            return 2;
        }
    }

    /* EAP_RUNTIME_DISABLED (project policy, 2026-08-22): the 17-species
     * generator remains available for offline/API generation, but EAP species
     * scattering is excluded from the production sensitivity campaign. Never
     * accept a selector and silently substitute another phase. Chl production
     * uses common phytoplankton absorption plus Stramski detritus scattering.
     * DOC-REF: OCRT_EAP_RUNTIME_DISABLED_AND_WINDOWS_RUN_POLICY_2026-08-22. */
    if (organic_group_cli_seen && cs.ccrr_mode &&
        cs.water_constituent_model == RT_WATER_CONSTITUENT_OCRT &&
        cs.ccrr_chl_mg_m3 > 0.0) {
        fprintf(stderr,
            "error: --ocrt-phyto-group is not supported yet (requested '%s').\n"
            "  The 17-species EAP generator is integrated, but EAP scattering is\n"
            "  disabled in the constituent RT path because the current L=200 moment\n"
            "  representation can reconstruct a negative mid-angle P11 for strongly\n"
            "  forward-peaked species.  Drop the option and rerun: phytoplankton then\n"
            "  contributes absorption only, while detritus supplies particle scattering.\n"
            "  Species scattering will be restored only with a validated component-level\n"
            "  truncation or a higher-order RT representation.\n",
            rt_iop_organic_group_name((organic_phyto_group_t)cs.organic_phyto_group));
        return 2;
    }

    /* Positive OCRT Chl requires common phytoplankton absorption plus the
     * Stramski detritus phase. EAP phase tables are deliberately not loaded.
     * Load this read-only pair once before any parallel water solve. */
    if ((cs.ccrr_mode &&
         cs.water_constituent_model == RT_WATER_CONSTITUENT_OCRT &&
         cs.ccrr_chl_mg_m3 > 0.0) || batch_full_grid_in) {
        char organic_dir[1024];
        const organic_phyto_group_t selected_group =
            (organic_phyto_group_t)cs.organic_phyto_group;
        /* EAP_RUNTIME_DISABLED: never request catalog-wide EAP loading. */
        const int load_all_phyto = 0;
        if (ocrt_organic_init_for_run(argv[0], selected_group, load_all_phyto,
                                      organic_dir, sizeof organic_dir) != 0) return 2;
        atexit(rt_iop_organic_free);
        if (cs.ccrr_mode &&
            cs.water_constituent_model == RT_WATER_CONSTITUENT_OCRT &&
            cs.ccrr_chl_mg_m3 > 0.0) {
            fprintf(stderr,
                    "# OCRT organic: Chl=%.10g mg m^-3 group=%s detritus_a440=%.10g m^-1 S_d=%.10g nm^-1 data=%s\n",
                    cs.ccrr_chl_mg_m3,
                    rt_iop_organic_group_name((organic_phyto_group_t)cs.organic_phyto_group),
                    cs.detritus_a440_m_inv, cs.detritus_slope_nm_inv, organic_dir);
        }
    }

    /* Positive TSM requires the packaged Ahn coefficients and matching vector
     * phase.  The diagnostic CCRR CSV reader may encounter positive MIN on any
     * row, so initialize its shared read-only tables up front as well. */
    if ((cs.ccrr_mode &&
         cs.water_constituent_model == RT_WATER_CONSTITUENT_OCRT &&
         cs.ccrr_min_g_m3 > 0.0) || ccrr_compare_in || batch_full_grid_in) {
        char tsm_dir[1024];
        if (ocrt_tsm_init_for_run(argv[0], tsm_dir, sizeof tsm_dir) != 0) return 2;
        atexit(rt_iop_ahn_mineral_free);
        if (cs.ccrr_mode &&
            cs.water_constituent_model == RT_WATER_CONSTITUENT_OCRT &&
            cs.ccrr_min_g_m3 > 0.0 && cs.ccrr_chl_mg_m3 <= 0.0) {
            cs.water_mie_phase_path =
                rt_iop_ahn_mineral_phase_path((ahn_species_t)cs.tsm_species);
            if (!cs.water_mie_phase_path) {
                fprintf(stderr, "error: no Ahn .mie phase for TSM species %s\n",
                        rt_iop_ahn_species_name((ahn_species_t)cs.tsm_species));
                return 2;
            }
            fprintf(stderr, "# TSM Ahn: C=%.10g g m^-3 species=%s data=%s phase=%s\n",
                    cs.ccrr_min_g_m3,
                    rt_iop_ahn_species_name((ahn_species_t)cs.tsm_species),
                    tsm_dir, cs.water_mie_phase_path);
        } else if (cs.ccrr_mode &&
                   cs.water_constituent_model == RT_WATER_CONSTITUENT_OCRT &&
                   cs.ccrr_min_g_m3 > 0.0) {
            fprintf(stderr, "# TSM Ahn: C=%.10g g m^-3 species=%s data=%s (vector phase mixed with organic particles)\n",
                    cs.ccrr_min_g_m3,
                    rt_iop_ahn_species_name((ahn_species_t)cs.tsm_species), tsm_dir);
        }
    }

    if (cs.scalar_ff_iop_path && cs.scalar_ff_iop_path[0]) {
        if (!have_wl || !(cs.wavelength_nm > 0.0)) {
            fprintf(stderr, "error: --scalar-ff-iop requires --wavelength in single-case mode\n");
            return 2;
        }
        if (!(cs.scalar_ff_bb_over_b > 0.0 && cs.scalar_ff_bb_over_b < 0.5)) {
            fprintf(stderr, "error: --scalar-ff-bbfrac must be in (0,0.5)\n");
            return 2;
        }
        double a_ff = 0.0, b_ff = 0.0;
        int lrc = scalar_ff_iop_lookup_main(cs.scalar_ff_iop_path, cs.wavelength_nm, &a_ff, &b_ff);
        if (lrc != 0) {
            fprintf(stderr, "error: failed to read --scalar-ff-iop '%s' at wavelength %.10g nm (rc=%d)\n",
                    cs.scalar_ff_iop_path, cs.wavelength_nm, lrc);
            return 2;
        }
        cs.fixed_bulk_iop_mode = 1;
        cs.fixed_bulk_phase_model = 4;
        cs.fixed_a_total_m_inv = a_ff;
        cs.fixed_b_total_m_inv = b_ff;
        cs.fixed_bb_total_m_inv = b_ff * cs.scalar_ff_bb_over_b;
        if (cs.fixed_bulk_lmax <= 0 || cs.fixed_bulk_lmax == 30) cs.fixed_bulk_lmax = 10;
        if (cs.fixed_bulk_phase_nphi <= 0) cs.fixed_bulk_phase_nphi = 720;
    }

    /* Constituent branches change only the in-water constituent-to-IOP/phase
     * construction.  Surface, atmospheric aerosol, wind and sunglint follow
     * the ordinary OCRT command line.  F_sun=1 is retained in one place solely
     * to preserve the frozen constituent-interface validation normalization. */
    if (cs.ccrr_mode) {
        cs.F_sun = 1.0;
        cs.cdom_ref_lambda_nm = 440.0;
        if (cs.S_cdom_nm_inv <= 0.0) cs.S_cdom_nm_inv = 0.014;
    }

    if (cs.fixed_bulk_iop_mode) {
        /* v1.10 FB-AER (2026-07-12, Jae-approved capability lift): the two
         * lines that used to sit here --
         *     mie_path = NULL;  user_aod = 0.0;
         * -- forcibly DISABLED the aerosol whenever the direct-IOP branch was used,
         * contradicting the "atmosphere presence is purely physical" rule
         * stated ~20 lines below (Rayleigh was never disabled here, only the
         * aerosol).  That made validation-matrix item #6 (composite water under
         * Rayleigh+aerosol) inexpressible: --aod was silently ignored (aod 0 vs
         * 0.1 produced BIT-IDENTICAL output).  Removing the forcing changes
         * results ONLY for the combination "fixed-bulk AND --aod-555/--aod-865 > 0 (or
         * --mie)", which could not be expressed before, so no existing result
         * can move: mie_path/user_aod keep their CLI defaults (NULL / 0.0,
         * declared ~line 1117/1133) when the flags are absent.  Proven by the
         * full gate set at this commit (Tier-0, production batch byte-compare,
         * repro_item2/#3/#5 runners - all unchanged). */
        /* B2 (2026-06-03): preserve CLI --wind-speed for the fixed-bulk path
         * (rough-surface internal reflection, Cox-Munk). Previously hard-forced
         * to 0 here, which made the in-water radiance wind-independent and left
         * the rough internal-reflection code (rt_water_rt.c) dead. Default
         * cs.wind_speed=0 (set ~line 627) keeps flat / regression-safe behavior
         * unless --wind-speed is explicitly passed. */
        cs.decouple_sunglint = 1;
        if (cs.fixed_bulk_lmax <= 0) cs.fixed_bulk_lmax = 0;   /* #23: 0 = per-path default (water-mie 200, LUT paths 30) */
    }

    /* ====================================================================
     * Atmosphere presence is now PURELY PHYSICAL — there is no on/off toggle.
     * Rayleigh is present iff the surface pressure is positive (tau_R scales
     * with P/1013.25, computed in the solver); aerosol is present iff the AOD
     * is positive.  This single derivation replaces every former mode-flag
     * forcing of rayleigh_on/aerosol_on, so water-property branches no longer
     * silently disable the atmosphere.  To run an
     * in-water-only test, pass --pressure 0 (and no --aod).
     * ==================================================================== */
    cs.rayleigh_on = (cs.pressure_hpa > 0.0) ? 1 : 0;
    cs.aerosol_on  = (user_aod      > 0.0) ? 1 : 0;

    /* 2026-07-15 Jae(G3) 가드.  (a) --aod-555/--aod-865 > 0 이면 --mie 필수.  (b) --aod 가
     * 0(기본)이면 에어로졸 옵션 일체 지정 불가.  배치 모드는 행이 aod 를
     * 켤 수 있으므로 (b) 를 CLI 단계에서 걸지 않고, (a) 는 행 병합 지점에서
     * 행 단위로 검사한다(run 루프 참조). */
    if (!batch_in && !batch_full_grid_in) {
        if (user_aod > 0.0 && !mie_path) {
            fprintf(stderr, "error: --aod-555/--aod-865 > 0 requires --mie <file> "
                    "(aerosol phase model; official subset T50/C50/M80C)\n");
            return 2;
        }
        if (user_aod <= 0.0 &&
            (mie_path || have_aer_h || have_aer_L_max ||
             have_trunc_aer_fit || have_trunc_aer_m || apply_nt_tau)) {
            fprintf(stderr, "error: aerosol options require --aod-555/--aod-865 > 0 "
                    "(aerosol is disabled at the default --aod-555 0)\n");
            return 2;
        }
    }

    /* T-V6-2026-05-10: m_max default dispatch — Sweep 2 결과 (V3 자체
     * m_max=16 부터 ρ_I deviation max 0.39%) 기반.
     * 사용자 가 --m-max 미명시 시 aerosol presence 에 따라 자동 dispatch:
     *   Rayleigh-only (aerosol_on=0): m_max=2 (V5 default 유지)
     *   Aerosol-on:                   m_max=16 (production 정확도)
     * 사용자 명시 시 (have_m_max=1) 그대로 따름. */
    if (!have_m_max) {
        if (cs.aerosol_on) {
            opts.fourier_m_max = 16;
            if (opts.verbose) {
                fprintf(stderr, "# m_max default dispatch: aerosol-on → 16\n");
            }
        }
        /* aerosol_on==0 시 default 2 그대로 (rt_types.h RT_OPTS_DEFAULT). */
    }

    /* T-V7-2026-05-10: aer_L_max default dispatch — V7 결정 fact 기반.
     * 사용자 명시 시 (have_aer_L_max=1) 그대로 따름. */
    if (!have_aer_L_max && cs.aerosol_on) {
        aer_L_max = 80;
        if (opts.verbose) {
            fprintf(stderr, "# aer_L_max default dispatch: aerosol-on → 80\n");
        }
    }

    /* Checkpoint-3 candidate: aerosol-only default forward-peak dispatch.
     * Optical domain: pure aerosol path (Rayleigh off), any surface BC.
     * δ-M N=64 is the best current black-surface validation option; NT τ is
     * deliberately NOT enabled by default because the no-NT variant gave lower
     * IQU MAPE in the Stage-4 dispatch audit.  Explicit truncation flags override. */
    /* Aerosol-only default forward-peak dispatch.  Real-angle log10-linear
     * truncation (Potter 1970-style, mass-conserving) is selected as the
     * production default because it gives the lowest aerosol-only IQU MAPE
     * across the FresnelOcean and BlackSurface 180-case validation sets
     * (FresnelOcean MAPE_I 0.432%, BlackSurface MAPE_I 0.300%; vs Wiscombe
     * δ-M N=64 at 1.295% / 0.659% respectively).  Highly absorbing aerosols
     * automatically bypass the operator via the A<0.1 threshold.
     * Explicit truncation flags override. */
    /* 2026-07-15 Jae(G3): 결합(레일리+에어로졸) 케이스까지 확장.  종전에는
     * !rayleigh_on 조건으로 에어로졸 단독에만 자동 적용되어, 결합 최소 명령이
     * loglin 없이 돌았다(명시 플래그 대비 I 0.51% @555 실측) -- #13 주석의
     * "완전 자동 프로토콜" 주장과 불일치.  플래그 삭제(자동설정 전제)에 맞춰
     * 조건을 aerosol_on 전체로 넓힌다.  δ 계열(DEBUG)이 오버라이드. */
    if (cs.aerosol_on &&
        !have_trunc_aer_m && !have_trunc_aer_fit) {
        use_loglin_trunc = 1;
        if (opts.verbose) {
            fprintf(stderr, "# aerosol forward-peak default dispatch: real-angle log10-linear truncation\n");
        }
    }

    /* v1.09 commit #13: remaining validated-protocol defaults + knob gates.
     * With S-005 fixed, an aerosol run now auto-selects the full validated
     * configuration (loglin + value kernel + L=80 + m=16 + 400 layers), so
     * the minimal command "--mie X --aod-555 Y --vector" reproduces the OSOAA
     * protocol bit-for-bit.  Convergence knobs may be re-specified at their
     * default values without any gate (documented commands keep working);
     * CHANGING them requires OCRT_ADVANCED=1 because every accuracy claim
     * (MAPE 0.48%) is tied to these values. */
    if (!have_aer_phase_kernel && cs.aerosol_on) {
        aer_use_value_kernel = 1;
        if (opts.verbose) fprintf(stderr, "# aer-phase-kernel default dispatch: aerosol-on -> value\n");
    }
    if (!have_n_layers && cs.aerosol_on) {
        opts.n_layers = 400;
        if (opts.verbose) fprintf(stderr, "# n_layers default dispatch: aerosol-on -> 400\n");
    }
    /* M2 (2026-07-14): 비에어로졸 실행의 수렴 노브 게이트 (기본값 재명시 통과). */
    if (!cs.aerosol_on && !ocrt_advanced_on()) {
        if (have_n_layers && opts.n_layers != 40) {
            fprintf(stderr, "error: --n-layers %d differs from the validated default 40 "
                    "(convergence knob; requires OCRT_ADVANCED=1)\n", opts.n_layers);
            return 2;
        }
        if (have_m_max && opts.fourier_m_max != 2) {
            fprintf(stderr, "error: --m-max %d differs from the Rayleigh-only default 2 "
                    "(convergence knob; requires OCRT_ADVANCED=1)\n", opts.fourier_m_max);
            return 2;
        }
    }
    if (cs.aerosol_on) {
        const char *adv_ = getenv("OCRT_ADVANCED");
        const int advanced_ = (adv_ && *adv_ && strcmp(adv_, "0") != 0);
        if (!advanced_) {
            if (have_aer_L_max && aer_L_max != 80) {
                fprintf(stderr, "error: --aer-l-max %d differs from the validated default 80 "
                        "(convergence knob; requires OCRT_ADVANCED=1)\n", aer_L_max);
                return 2;
            }
            if (have_n_layers && opts.n_layers != 400) {
                fprintf(stderr, "error: --n-layers %d differs from the validated default 400 "
                        "for aerosol runs (convergence knob; requires OCRT_ADVANCED=1)\n", opts.n_layers);
                return 2;
            }
            if (have_m_max && opts.fourier_m_max != 16) {
                fprintf(stderr, "error: --m-max %d differs from the validated default 16 "
                        "for aerosol runs (convergence knob; requires OCRT_ADVANCED=1)\n", opts.fourier_m_max);
                return 2;
            }
        }
    }

    opts.aerosol_mie_path = mie_path;
    opts.aerosol_aod = user_aod;
    opts.aerosol_aod_ref_nm = user_aod_ref_nm;

    /* (v1.09 commit #13 / BUG_SUSPECTS S-005: the rayleigh_on/aerosol_on
     * derivation was MOVED ABOVE the default-dispatch blocks — it used to sit
     * here, i.e. AFTER the three dispatch blocks that test cs.aerosol_on,
     * which left every auto-default (m_max 16, aer_L_max 80, loglin) dead
     * against the parse-time placeholder aerosol_on=0.  Measured symptom:
     * a plain "--mie X --aod-555 Y --vector" run silently used m_max=2, L=32,
     * no truncation, moment kernel, 40 layers -> -18.7% vs the validated
     * protocol AND conv=0.  See the relocated block before the m_max
     * dispatch. */

    /* ====================================================================
     * LUT mode default in-water quadrature: n_mu_water = 96 (2026-06-30).
     * Convergence (Brown_earth, highest-omega cases Csed=50 @555/660nm):
     * n_mu=48 is +0.44~0.47% above converged; converged by n_mu~80.  The
     * n_mu floor rises with single-scattering albedo, so a production LUT
     * sweep (which includes high-omega entries and even more forward-peaked
     * phyto phases, g~0.97) must not inherit the interactive default 48.
     * 96 = converged 80 + margin.  An explicit OCRT_ADVANCED=1 --n-mu-water
     * still overrides.  Physical basis, not tuning: pure quadrature
     * convergence verified against n_mu=96 self-reference.
     * ==================================================================== */
    if (cs.water_shared_grid) {
        if (cs.n_mu_water_override != 0) {
            fprintf(stderr, "error: --water-shared-grid and --n-mu-water are mutually exclusive\n");
            return 2;
        }
        cs.n_mu_water_override = (opts.n_mu > 0) ? opts.n_mu : 24;
        fprintf(stderr, "# water grid: SHARED with atmosphere (OSOAA-parity single mu set), "
                        "n_mu_water=%d (= n_mu)\n", cs.n_mu_water_override);
    }
    if (lut_mode && cs.surface == RT_SURFACE_OCEAN && cs.n_mu_water_override == 0) {
        /* v1.09 commit #15: PATH-dependent auto n_mu_water (사용자 방향 A).
         * The 96 blanket (above) was measured on the SCALAR fixed-bulk LUT
         * path (sharp forward-peak phase values integrated directly).  The
         * VECTOR IOP .mie phase path expands the phase to L<=~30 Legendre
         * moments, which smooths the angular structure; measured convergence
         * (blend g=0.981, 443nm, near-nadir, m_max=4 fixed):
         *   omega 0.487 / 0.891 / 0.955  ->  n_mu=24 within 0.002% of 32-40.
         * Auto value 32 = measured-converged 24 + one coverage-margin step
         * (single phase family / single band measured so far; lowering to 24
         * is a documented follow-up once coverage widens).  Scalar-LUT keeps
         * the worst-case-verified 96.  Explicit --n-mu-water still overrides
         * (OCRT_ADVANCED). */
        if (cs.water_mie_phase_path) {
            /* v1.09 #18 HOLD + S-007 (2026-07-08): the planned 24 default is
             * ON HOLD.  Grid-mode (view_as_node=0) radiances are extracted by
             * LINEAR INTERPOLATION between in-water quadrature nodes, and that
             * interpolation error — not field convergence — dominates:
             * measured at the vza=58deg grid node vs a direct-solve
             * (view_as_node=1) truth of 3.9454e-03, grid rrs_I is +40% at
             * n_mu=24, -40% at 32, -15% at 48 (Ed, an integral, is converged
             * to 2e-5 at 24; the FIELD itself converges at 24 — standalone
             * 24 vs 32 agree to 0.007%).  Until multi-view-node direct
             * extraction lands (add all refracted grid angles as zero-weight
             * nodes, one solve -> direct values; see S-007), keep 96 so the
             * legacy interpolation error stays at the historical level.  The
             * user-approved 24 default applies to the water QUADRATURE once
             * extraction no longer rides on it. */
            cs.n_mu_water_override = 96;
            fprintf(stderr, "# LUT mode (vector water phase): in-water n_mu_water = 96 "
                            "(S-007 hold: grid view interpolation; see BUG_SUSPECTS)\n");
        } else {
            cs.n_mu_water_override = 96;
            fprintf(stderr, "# LUT mode: in-water n_mu_water default raised to 96 "
                            "(interactive default 48 is +0.4%% under-converged at high omega)\n");
        }
    }

    opts.aerosol_l_max = aer_L_max;
    opts.aerosol_theta_cut_deg = theta_cut_deg;
    opts.aerosol_delta_m_N = delta_m_N;
    opts.aerosol_apply_nt_tau = apply_nt_tau;

    if (ccrr_compare_in) {
        const char *out_csv = ccrr_compare_out ? ccrr_compare_out : batch_out;
        if (!out_csv) {
            fprintf(stderr, "error: --ccrr-compare requires --ccrr-output or --output\n");
            return 2;
        }
        return run_ccrr_compare_csv(ccrr_compare_in, out_csv,
                                    water_aw_lut_path, water_psi_T_lut_path,
                                    &opts, &cs);
    }

    /* ==================================================================
     * 2026-07-14 M2 (Jae 지시): 실행 모드 단일화 및 정책 충돌 검사.
     *   단일 기하 모드: --vza 와 --raa 가 함께 명시(또는 --view-as-node 선언
     *     + 둘 다 명시)되면 view-as-node 정확 추출로 한 방향을 계산한다.
     *   LUT 격자 모드(기본): 기하 미지정 시 균일 vza/raa 격자(기본 2.5도,
     *     --lut-vza-step/--lut-raa-step 로 조정)를 한 번의 solve 로 산출한다.
     *   충돌 규칙: 단일 기하 지시와 격자/배치 지시가 동시에 오면 오류 종료.
     * ================================================================== */
    {
    /* 2026-07-15 Jae(G4): --surface 필수화(분기 결정 옵션이라 조용한 기본 금지).
     * ccrr-compare 모드는 내부에서 ocean 을 강제하므로 예외. */
    if (!have_surface && !ccrr_compare_in) {
        fprintf(stderr, "error: --surface is required "
                "(black | flat | black_fresnel_ocean | ocean)\n");
        return 2;
    }

        const int single_geo_req = (have_vza && have_raa) || single_mode_flag;
        const int grid_req = lut_mode || have_vza_step || have_raa_step;
        if (have_vza != have_raa) {
            fprintf(stderr, "error: single-geometry mode requires BOTH --vza and --raa "
                            "(got only one); omit both for LUT grid mode\n");
            return 2;
        }
        if (single_geo_req && grid_req) {
            fprintf(stderr, "error: conflicting mode flags: single geometry (--vza/--raa)"
                            " vs LUT grid (--output-full-grid/"
                            "--lut-vza-step/--lut-raa-step). Pick one.\n");
            return 2;
        }
        if (single_geo_req && (batch_in || batch_full_grid_in)) {
            fprintf(stderr, "error: conflicting mode flags: single geometry vs batch "
                            "(--batch/--batch-full-grid). Pick one.\n");
            return 2;
        }
        if (batch_in && batch_full_grid_in) {
            fprintf(stderr, "error: conflicting mode flags: --batch vs --batch-full-grid. "
                            "Pick one.\n");
            return 2;
        }
        if (!single_geo_req && !batch_in && !batch_full_grid_in && !ccrr_compare_in) {
            lut_mode = 1;   /* 기본 모드 = LUT 격자 */
            if (single_mode_flag == 0 && have_vza == 0) {
                cs.vza_deg = 0.0; cs.raa_deg = 0.0;   /* 폐기되는 기저 solve 용 기본 기하 */
            }
        }
        /* 단일 기하 모드에서는 노드 추가 추출을 항상 켠다 (M2-4 통합). */
        if (single_geo_req) {
            opts.view_as_node = 1;
            if (fabs(cs.vza_deg) <= 1.0e-15)
                warn_value_kernel_exact_nadir_once();
        }
    }

    /* M2 출력 파일명 규칙 (LUT 격자 모드, 비배치): 미지정 -> result<N>.csv,
     * 확장자 검사(.csv 만 지원). --batch 의 --lut-output 은 디렉터리라 제외. */
    static char ocrt_lut_out_buf[512];
    if (lut_mode && !batch_in) {
        if (ocrt_resolve_csv_name(lut_output, ocrt_lut_out_buf,
                                  sizeof ocrt_lut_out_buf) != 0) return 2;
        lut_output = ocrt_lut_out_buf;
        fprintf(stderr, "# LUT grid output: %s (vza step %.4g deg, raa step %.4g deg, vza max %.4g deg)\n",
                lut_output, lut_vza_step, lut_raa_step, lut_vza_max);
    }

    if (batch_in) {
        if (!batch_out) {
            fprintf(stderr, "error: --batch requires --output\n");
            return 2;
        }
        rt_batch_aerosol_config_t aer_cfg = {
            .mie_path = mie_path,
            .aer_L_max = aer_L_max,
            .theta_cut_deg = theta_cut_deg,
            .delta_m_N = delta_m_N,
            .apply_nt_tau = apply_nt_tau,
            .use_loglin_trunc = use_loglin_trunc
        };
        /* v1.01: propagate LUT options to batch via opts. */
        if (lut_mode) {
            opts.lut_enable = 1;
            opts.lut_output_dir = lut_output;  /* re-used as directory */
            opts.lut_vza_step = lut_vza_step;
            opts.lut_raa_step = lut_raa_step;
            opts.lut_vza_max  = lut_vza_max;
            if (!lut_output) {
                fprintf(stderr, "error: --lut with --batch requires --lut-output <dir>\n");
                return 2;
            }
            /* Create the directory (mkdir, ignore EEXIST). */
            if (mkdir(lut_output, 0755) != 0 && errno != EEXIST) {
                fprintf(stderr, "error: cannot create LUT output dir '%s'\n", lut_output);
                return 1;
            }
        }
        return rt_io_run_batch(batch_in, batch_out, &opts, &cs,
                               (cs.aerosol_on && mie_path) ? &aer_cfg : NULL);
    }

    if (!(have_sza && have_wl)) {
        fprintf(stderr, "error: --sza and --wavelength are required in every mode\n");
        usage(stderr, argv[0]);
        return 2;
    }
    if (!rt_spectral_wavelength_supported(cs.wavelength_nm)) {
        fprintf(stderr,
                "error: wavelength %.12g nm is outside OCRT elastic support %.0f-%.0f nm\n",
                cs.wavelength_nm, RT_SPECTRAL_MIN_NM, RT_SPECTRAL_MAX_NM);
        return 2;
    }
    /* 2026-07-15 Jae(G4): 해수면(ocean/black_fresnel_ocean)은 --wind-speed 필수(0 허용 --
     * 0 = 평면 프레넬 명시).  배치 모드는 행이 wind 를 줄 수 있으므로 행
     * 병합 지점에서 검사한다. */
    if (!batch_in && !batch_full_grid_in &&
        (cs.surface == RT_SURFACE_OCEAN || cs.surface == RT_SURFACE_BLACK_FRESNEL_OCEAN) &&
        !have_wind) {
        fprintf(stderr, "error: --wind-speed is required for --surface ocean/black_fresnel_ocean "
                "(explicit 0 = flat Fresnel)\n");
        return 2;
    }
    rt_result_t res = {0};
    int rc;

    /* v1.02: SOS-integrated gas absorption setup (BEFORE solver call,
     * replaces v1.01 post-RT Beer-Lambert). */
    rt_absorption_t abs_state_local = {0};
    if (use_absorption) {
        int abs_arc = rt_absorption_init(afgl_dir, xsec_dir, atm_profile,
                                          gas_col_override, &abs_state_local);
        if (abs_arc < 0) {
            fprintf(stderr, "error: rt_absorption_init failed rc=%d\n", abs_arc);
            return 1;
        }
        opts.abs_state = &abs_state_local;
        double tau_abs_col = rt_absorption_tau_total(&abs_state_local, cs.wavelength_nm);
        fprintf(stderr, "rt_absorption (SOS-integrated): λ=%.1fnm τ_abs_col=%.4e\n",
                cs.wavelength_nm, tau_abs_col);
        for (int g = 0; g < RT_N_GAS; ++g) {
            double tg = rt_absorption_tau_per_gas(&abs_state_local, g, cs.wavelength_nm);
            fprintf(stderr, "  %s: column=%.3e mol/cm²%s  τ_abs_col=%.4e\n",
                    RT_GAS_NAMES[g], abs_state_local.column_eff[g],
                    abs_state_local.column_overridden[g] ? " (overridden)" : "",
                    tg);
        }
    }

    /* v1.2 LUT-PERF (2026-07-21): dispatch ocean full-grid/batch before
     * the ordinary single-case solve.  The previous control flow executed one
     * complete coupled ocean solve, discarded it, and then entered the grid
     * driver.  That added a full water+atmosphere solution to every LUT job and
     * also duplicated LUT/Mie/runtime initialization. */
    if (batch_full_grid_in || (lut_mode && cs.surface == RT_SURFACE_OCEAN)) {
        if (cs.surface != RT_SURFACE_OCEAN) {
            fprintf(stderr, "error: ocean full-grid mode requires --surface ocean\n");
            return 2;
        }
        rt_aerosol_runtime_options_t fg_ropts = {
            .aer_L_max = aer_L_max, .theta_cut_deg = theta_cut_deg,
            .delta_m_N = delta_m_N, .apply_nt_tau = apply_nt_tau,
            .use_loglin_trunc = use_loglin_trunc,
            .loglin_mu1 = 0.8, .loglin_mu2 = 0.94, .loglin_threshold = 0.1
        };
        if (batch_full_grid_in) {
            return run_ocean_full_grid_batch(batch_full_grid_in, &cs, &opts,
                                             water_aw_lut_path, water_psi_T_lut_path,
                                             mie_path, user_aod, user_aod_ref_nm, &fg_ropts,
                                             lut_vza_step, lut_vza_max,
                                             lut_raa_step, have_n_water,
                                             have_wind /* 2026-07-15 G4 */);
        }
        return run_ocean_rrs_full_grid_csv(lut_output, &cs, &opts,
                                           water_aw_lut_path, water_psi_T_lut_path,
                                           mie_path, user_aod, user_aod_ref_nm, &fg_ropts,
                                           lut_vza_step, lut_vza_max,
                                           lut_raa_step);
    }

    if (cs.surface != RT_SURFACE_OCEAN &&
        cs.aerosol_on && mie_path && user_aod > 0.0) {
        mie_data_t mie = {0};
        if (read_mie_file(mie_path, &mie) != 0) {
            fprintf(stderr, "error: read_mie_file('%s') failed\n", mie_path);
            return 1;
        }

        rt_aerosol_runtime_options_t ropts = {
            .aer_L_max = aer_L_max,
            .theta_cut_deg = theta_cut_deg,
            .delta_m_N = delta_m_N,
            .apply_nt_tau = apply_nt_tau,
            .use_loglin_trunc = use_loglin_trunc,
            .loglin_mu1 = 0.8,
            .loglin_mu2 = 0.94,
            .loglin_threshold = 0.1,
            .use_value_kernel = aer_use_value_kernel
        };
        rt_aerosol_input_t aer = {0};
        rt_aerosol_runtime_diag_t adiag = {0};
        if (rt_aerosol_runtime_prepare(&mie, cs.wavelength_nm, user_aod, user_aod_ref_nm,
                                       &ropts, &aer, &adiag) != 0) {
            fprintf(stderr, "error: aerosol runtime prepare failed\n");
            mie_data_free(&mie);
            return 1;
        }

        aod_target_used = adiag.aod_target;
        aod_ratio_used = adiag.extinction_ratio;
        if (opts.verbose) {
            fprintf(stderr, "# aerosol spectral AOD: ref=%.8g at %.1f nm, extinction_ratio=%.8g, AOD(%.1f nm)=%.8g\n",
                    adiag.aod_ref, adiag.aod_ref_nm, adiag.extinction_ratio,
                    cs.wavelength_nm, adiag.aod_target);
        }
        if (apply_nt_tau) {
            fprintf(stderr, "# NT scaling: f=%.6f τ_a %.6f→%.6f, ω_a %.6f→%.6f\n",
                    adiag.f_fwd, adiag.aod_target, adiag.tau_a_eff,
                    adiag.ssa_raw, adiag.ssa_a_eff);
        }

        rc = rt_solve_case_pol_aerosol(&cs, &opts, &aer, &res);
        rt_aerosol_runtime_free(&aer);
        mie_data_free(&mie);
    } else if (cs.surface == RT_SURFACE_OCEAN) {
        /* 2026-07-15 Jae(G4): Step57 의 Quan-Fry(λ) 자동 기본을 폐지하고
         * 기본을 스칼라 1.34 로 고정한다(--n-water 미지정 시).  근거:
         * OSOAA/AF1982 기준이 1.34 고정이라 대조 정합이 우선이며, 실제로
         * 종전 자동 기본은 OSOAA 대조에서 굴절률 불일치를 만들고 있었다.
         * 골든 러너(07-12 동결)는 동결 당시 물리(Quan-Fry per-λ 값)를
         * --n-water 로 명시 핀하여 bit 를 보존한다.  Quan-Fry 함수
         * rt_water_iop_n_real() 자체는 수중 IOP 내부 용도로 존치. */
        const char *dbg_nw = getenv("OCRT_DUMP_N_WATER");
        if (dbg_nw && *dbg_nw) {
            fprintf(stderr, "OCRT_N_WATER_USED wl=%.6f n=%.12f source=%s\n",
                    cs.wavelength_nm, cs.n_water, have_n_water ? "explicit" : "quan_fry_auto");
        }
        /* ====================================================================
         * Phase B.4 Stage 1-2c: in-water RT with full atm-ocean coupling.
         * cs.F_sun is treated as TOA solar irradiance (Stage 2a+).
         * aerosol is included via aer (Stage 2c).
         * ==================================================================== */
        rt_water_iop_lut_t aw_lut;
        rt_water_iop_psi_T_lut_t psi_T_lut = {0};
        if (rt_water_iop_lut_load(water_aw_lut_path, &aw_lut) != 0) {
            fprintf(stderr, "error: failed to load a_w LUT '%s'\n", water_aw_lut_path);
            return 1;
        }
        int have_psi_T = (rt_water_iop_psi_T_load(water_psi_T_lut_path, &psi_T_lut) == 0);
        if (!have_psi_T) {
            fprintf(stderr,
                    "warn: psi_T LUT not loaded ('%s'); T-correction disabled\n",
                    water_psi_T_lut_path);
        }

        /* Stage 2c: prepare aerosol input (if cs.aerosol_on AND mie file provided) */
        rt_aerosol_input_t aer = {0};
        rt_aerosol_input_t *aer_ptr = NULL;
        mie_data_t mie_ocean = {0};
        int mie_loaded = 0;
        if (cs.aerosol_on && mie_path && user_aod > 0.0) {
            if (read_mie_file(mie_path, &mie_ocean) != 0) {
                fprintf(stderr, "error: read_mie_file('%s') failed\n", mie_path);
                rt_water_iop_lut_free(&aw_lut);
                if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
                return 1;
            }
            mie_loaded = 1;
            rt_aerosol_runtime_options_t ropts = {
                .aer_L_max = aer_L_max,
                .theta_cut_deg = theta_cut_deg,
                .delta_m_N = delta_m_N,
                .apply_nt_tau = apply_nt_tau,
                .use_loglin_trunc = use_loglin_trunc,
                .loglin_mu1 = 0.8,
                .loglin_mu2 = 0.94,
                .loglin_threshold = 0.1
            };
            rt_aerosol_runtime_diag_t adiag_ocean = {0};
            if (rt_aerosol_runtime_prepare(&mie_ocean, cs.wavelength_nm, user_aod, user_aod_ref_nm,
                                            &ropts, &aer, &adiag_ocean) != 0) {
                fprintf(stderr, "error: aerosol runtime prepare failed (ocean mode)\n");
                mie_data_free(&mie_ocean);
                rt_water_iop_lut_free(&aw_lut);
                if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
                return 1;
            }
            aod_target_used = adiag_ocean.aod_target;
            aod_ratio_used = adiag_ocean.extinction_ratio;
            aer_ptr = &aer;
            fprintf(stderr, "[ocean+aer] AOD_ref=%.4f@%.1fnm AOD_band=%.4f τ_eff=%.4f ω_a=%.4f L_max=%d\n",
                    adiag_ocean.aod_ref, adiag_ocean.aod_ref_nm, adiag_ocean.aod_target,
                    aer.tau_a, aer.ssa_a, aer.L_max);
        }
        rc = rt_solve_case_ocean(&cs, &opts, &aw_lut,
                                  (have_psi_T ? &psi_T_lut : NULL),
                                  aer_ptr, &res);

        if (aer_ptr) rt_aerosol_runtime_free(&aer);
        if (mie_loaded) mie_data_free(&mie_ocean);
        rt_water_iop_lut_free(&aw_lut);
        if (have_psi_T) rt_water_iop_psi_T_free(&psi_T_lut);
    } else {
        /* 2026-07-14 M1 (Jae 지시): 스칼라 경로 완전 삭제. 항상 벡터. */
        rc = rt_solve_case_pol(&cs, &opts, &res);
    }
    if (rc != 0) {
        if (batch_full_grid_in) {
            /* v1.10 FIX: the discarded base solve must not kill the batch.
             * The base CLI is only a TEMPLATE for --batch-full-grid rows
             * (e.g. --iop-phase-lut alone, with IOP values coming
             * from the jobs columns, makes the base solve fail -2 by
             * design).  Rows carry their own overrides; warn and proceed. */
            fprintf(stderr, "warning: base solve failed rc=%d; "
                            "proceeding to --batch-full-grid rows\n", rc);
        } else {
            fprintf(stderr, "error: rt_solve_case_pol rc=%d\n", rc);
            return 1;
        }
    }

    /* Ocean full-grid/batch modes are dispatched before the discarded
     * single-case path above (v1.2 LUT-PERF). */

    /* v1.01 LUT(격자) 모드: M2 기본 모드.  (구 "--lut requires --vector"
     * 가드는 두 플래그 삭제 + vector_mode 기본 1 로 사어가 되어 제거함,
     * 2026-07-15.) */
    if (lut_mode) {
        if (lut_vza_step <= 0 || lut_raa_step <= 0 ||
            lut_vza_max <= 0 || lut_vza_max >= 90.0) {
            fprintf(stderr, "error: invalid LUT grid step/max\n");
            return 2;
        }
        int n_vza = (int)(lut_vza_max / lut_vza_step) + 1;
        int n_raa = (int)(360.0 / lut_raa_step);
        double *vza_g = (double*)malloc(sizeof(double) * n_vza);
        double *raa_g = (double*)malloc(sizeof(double) * n_raa);
        for (int i = 0; i < n_vza; ++i) vza_g[i] = (double)i * lut_vza_step;
        for (int i = 0; i < n_raa; ++i) raa_g[i] = (double)i * lut_raa_step;

        rt_lut_grid_out_t lut = {
            .n_vza = n_vza, .vza_deg = vza_g,
            .n_raa = n_raa, .raa_deg = raa_g,
            .rho_I = (double*)calloc((size_t)n_vza * n_raa, sizeof(double)),
            .rho_Q = (double*)calloc((size_t)n_vza * n_raa, sizeof(double)),
            .rho_U = (double*)calloc((size_t)n_vza * n_raa, sizeof(double)),
            .T_diff_dn_dir   = (double*)calloc((size_t)n_vza * n_raa, sizeof(double)),
            .T_sg_up_dir     = (double*)calloc((size_t)n_vza * n_raa, sizeof(double)),
            .T_total_up_dir  = (double*)calloc((size_t)n_vza * n_raa, sizeof(double)),
        };
        if (!lut.rho_I || !lut.rho_Q || !lut.rho_U ||
            !lut.T_diff_dn_dir || !lut.T_sg_up_dir || !lut.T_total_up_dir) {
            fprintf(stderr, "error: LUT output allocation failed\n");
            return 1;
        }

        /* Re-solve with LUT. SOS runs once. */
        if (mie_path) {
            /* Re-initialize aerosol input for LUT solve.
             * Reuse: most users call --lut as a single-case query, separate
             * from a regression batch. Re-running the runtime expansion is
             * acceptable cost (single Mie expansion). */
            mie_data_t mie2 = {0};
            rt_aerosol_input_t aer2 = {0};
            rt_aerosol_runtime_diag_t adiag2 = {0};
            rt_aerosol_runtime_options_t ropts2 = {
                .aer_L_max = aer_L_max,
                .theta_cut_deg = theta_cut_deg,
                .delta_m_N = delta_m_N,
                .apply_nt_tau = apply_nt_tau,
                .use_loglin_trunc = use_loglin_trunc,
                .loglin_mu1 = 0.8, .loglin_mu2 = 0.94, .loglin_threshold = 0.1,
                .use_value_kernel = aer_use_value_kernel
            };
            if (read_mie_file(mie_path, &mie2) != 0 ||
                rt_aerosol_runtime_prepare(&mie2, cs.wavelength_nm,
                                           user_aod, user_aod_ref_nm, &ropts2, &aer2, &adiag2) != 0) {
                fprintf(stderr, "error: LUT aerosol re-init failed\n");
                return 1;
            }
            rc = rt_solve_case_pol_lut(&cs, &opts, &aer2, &res, &lut);
            rt_aerosol_runtime_free(&aer2);
            mie_data_free(&mie2);
        } else {
            rc = rt_solve_case_pol_lut(&cs, &opts, NULL, &res, &lut);
        }
        if (rc != 0) {
            fprintf(stderr, "error: rt_solve_case_pol_lut rc=%d\n", rc);
            return 1;
        }

        FILE *fp = lut_output ? fopen(lut_output, "w") : stdout;
        if (!fp) {
            fprintf(stderr, "error: cannot open LUT output '%s'\n", lut_output);
            return 1;
        }
        /* Single-case LUT CSV — column naming convention identical to batch path
         * (rt_io.c). All transmittances reported at 0+ side (atm-side of surface). */
        fprintf(fp, "sza_deg,wavelength_nm,vza_deg,raa_deg,"
                    "rho_I,rho_Q,rho_U,"
                    "direct_transmittance,irradiance_transmittance,"
                    "diffuse_dn_hemi,diffuse_dn_dir,"
                    "sunglint_up_flux_ratio_diag,toa_up_flux_ratio_diag\n");
        const double T_irrad_case = lut.T_dir_dn + lut.T_diff_dn_hemi;
        for (int iv = 0; iv < n_vza; ++iv) {
            for (int ir = 0; ir < n_raa; ++ir) {
                const int idx = iv * n_raa + ir;
                fprintf(fp, "%.6f,%.3f,%.6f,%.6f,"
                            "%.10e,%.10e,%.10e,"
                            "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                        cs.sza_deg, cs.wavelength_nm,
                        vza_g[iv], raa_g[ir],
                        lut.rho_I[idx], lut.rho_Q[idx], lut.rho_U[idx],
                        lut.T_dir_dn, T_irrad_case,
                        lut.T_diff_dn_hemi,
                        lut.T_diff_dn_dir[idx], lut.T_sg_up_dir[idx],
                        lut.T_total_up_dir[idx]);
            }
        }
        if (lut_output) fclose(fp);
        fprintf(stderr, "rt_lut: emitted %d×%d = %d directions to %s\n",
                n_vza, n_raa, n_vza * n_raa, lut_output ? lut_output : "stdout");
        fprintf(stderr, "  direct_transmittance=%.6e  irradiance_transmittance=%.6e (case-level scalars, 0+ side)\n",
                lut.T_dir_dn, T_irrad_case);
        free(lut.rho_I); free(lut.rho_Q); free(lut.rho_U);
        free(lut.T_diff_dn_dir); free(lut.T_sg_up_dir); free(lut.T_total_up_dir);
        free(vza_g); free(raa_g);
        return 0;
    }
    /* v1.02: SOS-integrated absorption applied INSIDE solver via opts.abs_state.
     * No post-RT correction here (deprecated v1.01 path removed).
     * Cleanup of abs_state happens at end. */

    res.aerosol_aod_ref = user_aod;
    res.aerosol_aod_ref_nm = user_aod_ref_nm;
    res.aerosol_aod_band = aod_target_used;
    res.aerosol_extinction_ratio = aod_ratio_used;
    rt_io_print_single(&cs, &res);
    if (use_absorption) rt_absorption_free(&abs_state_local);
    return 0;
}
