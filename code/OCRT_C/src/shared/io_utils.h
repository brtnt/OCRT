/* shared/io_utils.h
 *
 * Line-based file parsing helpers used by both vrt and ocean .mie/.inp readers.
 *
 * Implementation in shared/io_utils.c.
 */
#ifndef OCRT_SHARED_IO_UTILS_H
#define OCRT_SHARED_IO_UTILS_H

#include <stdio.h>

/* Read a line from f into buf[bufsize], stripping trailing \n and \r.
 * Returns buf on success, NULL on EOF. */
char* read_line(FILE* f, char* buf, int bufsize);

/* C11-portable strdup (POSIX strdup not in strict C11, MSVC uses _strdup). */
char* my_strdup(const char* s);

/* Parse up to max_count whitespace-separated floats from line.
 * Stops at first non-numeric token. Returns count parsed. */
int parse_floats(const char* line, double* out, int max_count);

/* Returns 1 if first non-whitespace character of line is digit, +, -, or '.'. */
int line_starts_numeric(const char* line);

#endif /* OCRT_SHARED_IO_UTILS_H */
