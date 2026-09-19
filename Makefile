CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra
LDFLAGS ?=
PREFIX  ?= /usr/local

# Kept out of CFLAGS deliberately. dpkg-buildflags sets CFLAGS, and CFLAGS ?=
# means the environment wins, so anything essential put there is dropped in a
# package build - which would have shipped a single threaded binary quietly.
OPENMP  ?= -fopenmp

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
	$(CC) $(CFLAGS) $(OPENMP) $(GDAL_CFLAGS) -o $@ src/isofill.c $(GDAL_LIBS) $(LDFLAGS) $(OPENMP) -lm

$(LIB): src/isofill.c src/isofill.h
	$(CC) $(CFLAGS) $(OPENMP) -fPIC -shared -DISOFILL_NO_MAIN $(GDAL_CFLAGS) -o $@ src/isofill.c $(GDAL_LIBS) $(LDFLAGS) $(OPENMP) -lm

# A program that links the library and fills a small raster, asserting the
# surface it gets is the one a caller is entitled to: constraints untouched,
# the fill between two contours between them and monotonic, the version the
# header promised. Run against the built library from this directory, which
# is where the loader finds it on both platforms.
tests/run_lib: tests/run_lib.c $(LIB) src/isofill.h
	$(CC) $(CFLAGS) -Isrc -o $@ tests/run_lib.c -L. -lisofill $(LDFLAGS)

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
check: isofill tests/run_lib
	@out=$$(./isofill 2>&1); st=$$?; \
	if [ $$st -ne 2 ]; then echo "expected exit 2 after usage, got $$st"; exit 1; fi; \
	printf '%s\n' "$$out" | head -1 | grep -qx 'usage: isofill \[options\] <constraints.tif> <out.tif>' \
		|| { echo "first line of output is not the usage line"; printf '%s\n' "$$out" | head -3; exit 1; }; \
	echo "isofill runs: exit 2 and its usage text"
	@./isofill --version | grep -qx "isofill $$(sed -n 's/^#define ISOFILL_VERSION "\(.*\)"/\1/p' src/isofill.h)" \
		&& echo "isofill --version agrees with isofill.h" \
		|| { echo "isofill --version disagrees with isofill.h"; exit 1; }
	@LD_LIBRARY_PATH=. PATH=".:$$PATH" ./tests/run_lib

clean:
	rm -f isofill $(LIB) tests/run_lib

.PHONY: all install check clean
