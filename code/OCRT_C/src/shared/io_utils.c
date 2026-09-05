/* shared/io_utils.c
 *
 * Implementation of line-based parsing helpers. Sourced from vrt_solver.c v4.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "io_utils.h"

char* read_line(FILE* f, char* buf, int bufsize) {
    if (!fgets(buf, bufsize, f)) return NULL;
    size_t L = strlen(buf);
    while (L > 0 && (buf[L-1] == '\n' || buf[L-1] == '\r')) buf[--L] = '\0';
    return buf;
}

char* my_strdup(const char* s) {
    size_t L = strlen(s) + 1;
    char* p = (char*)malloc(L);
    if (p) memcpy(p, s, L);
    return p;
}

int parse_floats(const char* line, double* out, int max_count) {
    int count = 0;
    const char* s = line;
    while (*s && count < max_count) {
        while (*s && isspace((unsigned char)*s)) s++;
        if (!*s) break;
        /* Numeric-starting? Accept digits, +/- sign, or decimal point. */
        if (!(isdigit((unsigned char)*s) || *s == '-' || *s == '+' || *s == '.')) break;
        char* end = NULL;
        double v = strtod(s, &end);
        if (end == s) break;
        out[count++] = v;
        s = end;
    }
    return count;
}

int line_starts_numeric(const char* line) {
    while (*line && isspace((unsigned char)*line)) line++;
    if (!*line) return 0;
    return isdigit((unsigned char)*line) || *line == '-' || *line == '+' || *line == '.';
}
