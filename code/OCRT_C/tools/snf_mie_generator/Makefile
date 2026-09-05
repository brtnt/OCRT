CC ?= gcc
CFLAGS ?= -std=c11 -O2 -fopenmp -Wall -Wextra
LDLIBS ?= -lm

.PHONY: all clean
all: snf_mie_gen
snf_mie_gen: snf_mie_gen.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@
clean:
	rm -f snf_mie_gen
