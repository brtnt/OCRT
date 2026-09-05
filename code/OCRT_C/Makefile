# OCRT canonical build entry (auto-detectable). Targets: all(=ocrt), test-reader, clean
CC=gcc
CFLAGS=-std=c11 -O3 -march=cascadelake -ffp-contract=fast -fassociative-math -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc
SRC:=$(shell find src -name '*.c')
all: build/ocrt
build/ocrt: $(SRC)
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRC) -o $@ -lm
test-reader: build/test_fr631_mie_reader
build/test_fr631_mie_reader: tests/test_fr631_mie_reader.c src/shared/mie_io.c src/shared/io_utils.c src/shared/numerics.c
	@mkdir -p build
	$(CC) -std=c11 -O2 -Isrc $^ -o $@ -lm
clean:
	rm -f build/ocrt build/test_fr631_mie_reader
.PHONY: all test-reader clean
