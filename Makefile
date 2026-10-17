CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra
LDFLAGS ?=
PREFIX  ?= /usr/local

# Kept out of CFLAGS deliberately. dpkg-buildflags sets CFLAGS, and CFLAGS ?=
# means the environment wins, so anything essential put there is dropped in a
# package build - which would have shipped a single threaded binary quietly.
OPENMP  ?= -fopenmp

# One surface on every machine. A compiler may fuse a*b + c into one
# instruction that rounds once rather than twice, and GCC and clang both do
# wherever the target has one: arm64 always, x86-64 under -march=native. The
# fill's decisions turn on its arithmetic to the last bit, and fused, Apple
# Silicon's surface of the golden square differed from x86-64's in 4,829 of
# 1,442,401 cells, by up to 29 m (danu #194). Baseline x86-64 has no fused
# instruction, so there this changes nothing. := and not ?=, for the reason
# OPENMP is kept out of CFLAGS: nothing the environment sets drops it.
FPFLAGS := -ffp-contract=off

GDAL_CFLAGS := $(shell gdal-config --cflags 2>/dev/null)
GDAL_LIBS   := $(shell gdal-config --libs 2>/dev/null)

ifeq ($(strip $(GDAL_LIBS)),)
$(error gdal-config not found - install libgdal-dev)
endif

# The shared library is the same source compiled without main: one
# implementation, two ways to call it. .dll under MSYS2, .so elsewhere.
ifneq (,$(findstring Msys,$(shell uname -o 2>/dev/null)))
SOEXT := dll
else
SOEXT := so
endif
LIB := libisofill.$(SOEXT)

all: isofill $(LIB)

isofill: src/isofill.c src/isofill.h
	$(CC) $(CFLAGS) $(FPFLAGS) $(OPENMP) $(GDAL_CFLAGS) -o $@ src/isofill.c $(GDAL_LIBS) $(LDFLAGS) $(OPENMP) -lm

$(LIB): src/isofill.c src/isofill.h
	$(CC) $(CFLAGS) $(FPFLAGS) $(OPENMP) -fPIC -shared -DISOFILL_NO_MAIN $(GDAL_CFLAGS) -o $@ src/isofill.c $(GDAL_LIBS) $(LDFLAGS) $(OPENMP) -lm

# A program that links the library and fills a small raster, asserting the
# surface it gets is the one a caller is entitled to: constraints untouched,
# the fill between two contours between them and monotonic, the version the
# header promised. Run against the built library from this directory, which
# is where the loader finds it on both platforms.
# $(OPENMP) here too, not because the test threads anything but so that it can
# ask omp_get_num_procs() what the machine has and check the default against
# it. Built without it, the test's own _OPENMP is undefined while the
# library's is defined, and a test that branched on it asserted the
# single-threaded contract against a threaded library.
# The first pass against the walk it replaced, cell for cell; see the test.
# It includes src/isofill.c to reach the static functions, so it is built from
# the source rather than linked to the library.
tests/pass1_ref: tests/pass1_ref.c src/isofill.c src/isofill.h
	$(CC) $(CFLAGS) $(FPFLAGS) $(OPENMP) $(GDAL_CFLAGS) -o $@ tests/pass1_ref.c $(GDAL_LIBS) $(LDFLAGS) $(OPENMP) -lm

# The second pass against what Laplace's equation must give: zero between a
# coast and the drawn edge, a flat floor in a ring, and random ground against a
# solve in double precision.
tests/pass2_solved: tests/pass2_solved.c $(LIB) src/isofill.h
	$(CC) $(CFLAGS) $(FPFLAGS) $(OPENMP) -Isrc -o $@ tests/pass2_solved.c -L. -lisofill $(LDFLAGS) $(OPENMP) -lm

# The second pass out of core against the same raster solved whole, to the
# bit: the binary run twice, the second time with no memory to speak of.
tests/pass2_ooc: tests/pass2_ooc.c
	$(CC) $(CFLAGS) $(FPFLAGS) $(GDAL_CFLAGS) -o $@ tests/pass2_ooc.c $(GDAL_LIBS) $(LDFLAGS) -lm

tests/run_lib: tests/run_lib.c $(LIB) src/isofill.h
	$(CC) $(CFLAGS) $(FPFLAGS) $(OPENMP) -Isrc -o $@ tests/run_lib.c -L. -lisofill $(LDFLAGS) $(OPENMP)

install: isofill $(LIB)
	install -d $(DESTDIR)$(PREFIX)/bin $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/include
	install -m 755 isofill $(DESTDIR)$(PREFIX)/bin/isofill
	install -m 644 $(LIB) $(DESTDIR)$(PREFIX)/lib/$(LIB)
	install -m 644 src/isofill.h $(DESTDIR)$(PREFIX)/include/isofill.h

# A real check, not a cheerful one: the old target ended in || true and could
# not fail. The binary is run and two things are required of it together -
# the exit status isofill gives after printing usage, which is 2, and the
# first line of that usage text exactly. A loader error can manage one of
# those; not both. Then the library, through the test program.
check: isofill tests/run_lib tests/pass1_ref tests/pass2_solved tests/pass2_ooc
	@out=$$(./isofill 2>&1); st=$$?; \
	if [ $$st -ne 2 ]; then echo "expected exit 2 after usage, got $$st"; exit 1; fi; \
	printf '%s\n' "$$out" | head -1 | grep -qx 'usage: isofill \[options\] <constraints.tif> <out.tif>' \
		|| { echo "first line of output is not the usage line"; printf '%s\n' "$$out" | head -3; exit 1; }; \
	echo "isofill runs: exit 2 and its usage text"
	@./isofill --version | grep -qx "isofill $$(sed -n 's/^#define ISOFILL_VERSION "\(.*\)"/\1/p' src/isofill.h)" \
		&& echo "isofill --version agrees with isofill.h" \
		|| { echo "isofill --version disagrees with isofill.h"; exit 1; }
	@LD_LIBRARY_PATH=. PATH=".:$$PATH" ./tests/run_lib
	@./tests/pass1_ref
	@LD_LIBRARY_PATH=. PATH=".:$$PATH" ./tests/pass2_solved
	@./tests/pass2_ooc

clean:
	rm -f isofill $(LIB) tests/run_lib tests/pass1_ref tests/pass2_solved tests/pass2_ooc

.PHONY: all install check clean
