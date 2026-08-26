CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -fopenmp
LDFLAGS ?=
PREFIX  ?= /usr/local

GDAL_CFLAGS := $(shell gdal-config --cflags 2>/dev/null)
GDAL_LIBS   := $(shell gdal-config --libs 2>/dev/null)

ifeq ($(strip $(GDAL_LIBS)),)
$(error gdal-config not found - install libgdal-dev)
endif

all: isofill

isofill: src/isofill.c
	$(CC) $(CFLAGS) $(GDAL_CFLAGS) -o $@ $< $(GDAL_LIBS) $(LDFLAGS) -lm

install: isofill
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 isofill $(DESTDIR)$(PREFIX)/bin/isofill

clean:
	rm -f isofill

.PHONY: all install clean
