#ifndef OCRT_RT_WINDOWS_COMPAT_H
#define OCRT_RT_WINDOWS_COMPAT_H

/* Native Windows/MinGW portability. Linux/POSIX builds are unchanged.
 * This header deliberately adds no CLI or runtime physics option.
 * DOC-REF: OCRT_EAP_RUNTIME_DISABLED_AND_WINDOWS_RUN_POLICY_2026-08-22. */
#ifdef _WIN32
#include <stdlib.h>
#include <time.h>

static inline int ocrt_setenv_win(const char *name, const char *value, int overwrite)
{
    (void)overwrite;
    return _putenv_s(name, value ? value : "");
}
static inline int ocrt_unsetenv_win(const char *name)
{
    return _putenv_s(name, "");
}
#define setenv ocrt_setenv_win
#define unsetenv ocrt_unsetenv_win

/* Modern MinGW-w64 exposes clock_gettime. Keep a C11 fallback for toolchains
 * that do not expose CLOCK_MONOTONIC. It is used for diagnostics/timing only. */
#ifndef CLOCK_MONOTONIC
#define CLOCK_MONOTONIC 1
static inline int ocrt_clock_gettime_win(int clock_id, struct timespec *ts)
{
    (void)clock_id;
    return timespec_get(ts, TIME_UTC) == TIME_UTC ? 0 : -1;
}
#define clock_gettime ocrt_clock_gettime_win
#endif
#endif

#endif
