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

all: isofill

isofill: src/isofill.c
	$(CC) $(CFLAGS) $(OPENMP) $(GDAL_CFLAGS) -o $@ $< $(GDAL_LIBS) $(LDFLAGS) $(OPENMP) -lm

install: isofill
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 isofill $(DESTDIR)$(PREFIX)/bin/isofill

# A real check, not a cheerful one: the old target ended in || true and could
# not fail. The binary is run and two things are required of it together -
# the exit status isofill gives after printing usage, which is 2, and the
# first line of that usage text exactly. A loader error can manage one of
# those; not both. The same check danu's CI makes of the binary it builds.
check: isofill
	@out=$$(./isofill 2>&1); st=$$?; \
	if [ $$st -ne 2 ]; then echo "expected exit 2 after usage, got $$st"; exit 1; fi; \
	printf '%s\n' "$$out" | head -1 | grep -qx 'usage: isofill \[options\] <constraints.tif> <out.tif>' \
		|| { echo "first line of output is not the usage line"; printf '%s\n' "$$out" | head -3; exit 1; }; \
	echo "isofill runs: exit 2 and its usage text"

clean:
	rm -f isofill

.PHONY: all install check clean
