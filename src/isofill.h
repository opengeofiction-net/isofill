/*
 * isofill as a library: the in-core fill, arrays in and an array out.
 *
 * The binary's own in-core path calls isofill_run, so a caller of this header
 * gets the surface the binary writes - the same code, not the same idea of it.
 * The out-of-core banding for rasters above --max-mem is the binary's alone;
 * a caller with a raster that size wants the binary.
 *
 * SPDX-License-Identifier: Artistic-1.0 OR GPL-1.0-or-later
 */
#ifndef ISOFILL_H
#define ISOFILL_H

#define ISOFILL_VERSION "0.8.0"

/* The value a constraints cell holds when it is not a constraint and no
 * nodata value is given; also what --no-pass2 leaves in cells it declined.
 * The implementation's NO_ELEV is this, so the two cannot drift. */
#define ISOFILL_NO_ELEV (-32768)

typedef struct {
    int radius;        /* search radius in cells; --radius */
    int barrier;       /* widen constraints for the sight test; --barrier */
    double grad_min;   /* least gradient worth interpolating across; --grad-min */
    int pass2;         /* 1: diffuse the declined cells; 0: --no-pass2 */
    int threads;       /* OpenMP threads; <= 0 leaves the runtime's setting */
} isofill_params;

/* A version string, so a caller can refuse a library that is not the one it
 * was built against; there is no stable ABI. */
const char *isofill_version(void);

/* The defaults the command line uses when a flag is not given - radius 20,
 * barrier 1, grad_min 0.02, pass 2 on, threads left to the runtime. A caller
 * that sets only what it means to set gets the binary's behaviour for the
 * rest, by construction rather than by copying numbers. The binary's own
 * main() starts from this. */
void isofill_params_default(isofill_params *p);

/* What the in-core fill holds for a raster of this size, in megabytes, as the
 * binary reckons it when deciding whether to band. A caller with less than
 * this to spare wants the binary, which bands; the library does not. */
double isofill_whole_mb(int cols, int rows);

/*
 * Fill. constraints is cols*rows floats; a cell is a constraint unless it
 * equals nodata (when has_nodata) or ISOFILL_NO_ELEV (otherwise). mask and water are
 * cols*rows bytes or NULL, with the meaning of --mask and --water. out is
 * cols*rows floats and receives the surface; with pass2 set, cells outside the
 * mask come back as 0, as the binary writes them. Returns the number of cells
 * the first pass set, or a negative value: -1 for bad arguments, -2 for memory.
 */
long long isofill_run(const float *constraints, int has_nodata, double nodata,
                      const unsigned char *mask, const unsigned char *water,
                      int cols, int rows, const isofill_params *params,
                      float *out);

/*
 * The same fill, with the first pass kept.
 *
 * pass1_out is cols*rows floats or NULL. When it is given it receives the
 * surface as the first pass left it - which is what --no-pass2 writes, sentinels
 * and all - before the second pass overwrites the cells the first declined.
 * isofill_run is this with NULL.
 *
 * It exists because the first pass is nearly all of the work: on a 3601x2401
 * raster it is 0.56 s of a 0.67 s run. A caller that wants both the finished
 * surface and a reading of what the first pass could not answer - which is how
 * the editor draws the ground its contours do not describe - would otherwise
 * run the whole thing twice and throw one surface away.
 */
long long isofill_run_ex(const float *constraints, int has_nodata, double nodata,
                         const unsigned char *mask, const unsigned char *water,
                         int cols, int rows, const isofill_params *params,
                         float *out, float *pass1_out);

#endif
