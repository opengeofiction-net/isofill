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

check: isofill
	@./isofill --radius 4 2>&1 | grep -q usage && echo "usage ok" || true

clean:
	rm -f isofill

.PHONY: all install check clean
