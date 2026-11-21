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

#define ISOFILL_VERSION "0.16.1"

/* The value a constraints cell holds when it is not a constraint and no
 * nodata value is given; also what --no-pass2 leaves in cells it declined.
 * The implementation's NO_ELEV is this, so the two cannot drift. */
#define ISOFILL_NO_ELEV (-32768)

typedef struct {
    int radius;        /* search radius in cells; --radius */
    int barrier;       /* widen constraints for the sight test; --barrier */
    double grad_min;   /* least gradient worth interpolating across; --grad-min */
    int pass2;         /* 1: diffuse the declined cells; 0: --no-pass2 */
    int threads;       /* OpenMP threads; <= 0 is all but two, as --threads */
                       /* 0 used to mean "leave the runtime's setting alone",
                        * and no longer does: the fill sets the count on every
                        * call. A host that has chosen its own count and wants
                        * it kept should pass that number, not a zero. */
    int edge_void;     /* 1: the raster's edge seeds the void, as ground nobody
                        * drew; 0: --no-edge-void, only outside the mask does.
                        * A raster that is a piece of a larger drawn area - an
                        * editor's working set inside its zone - has ground
                        * beyond its edge, and says where it has none by
                        * leaving those edge cells out of the mask */
} isofill_params;

/* A version string, so a caller can refuse a library that is not the one it
 * was built against; there is no stable ABI. */
const char *isofill_version(void);

/* How many threads a given `threads` means: the number itself when positive,
 * and otherwise all but two, which is the command line's default and now the
 * library's. Exposed so a caller can say what it is about to do - the binary
 * prints it - without a second copy of the arithmetic. */
int isofill_resolve_threads(int threads);

/* The defaults the command line uses when a flag is not given - radius 20,
 * barrier 1, grad_min 0.02, pass 2 on, threads all but two. A caller
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
 *
 * -2 is not the whole of what running out of memory does here. The second pass
 * reports it; the first prints and calls exit, which a library has no business
 * doing to the program that loaded it, and converting it is the next thing to
 * do. Until then a caller of this cannot treat a return as the only way out.
 * isofill_diffuse, which is the second pass alone, has no such hole.
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

/*
 * The same fill, with cells held at the caller's values by the second pass.
 *
 * held is cols*rows floats or NULL. A cell is held unless it is
 * ISOFILL_NO_ELEV or NaN; a held cell is written over whatever the first pass
 * left there and the second pass keeps it, and the ground round it settles to
 * it. It counts as ground the first pass answered, so the void flood neither
 * takes it nor crosses it - and so does the ground within the radius of it
 * (the square of side 2 radius + 1) that no constraint reaches: that is
 * declined, and solved, rather than void. Outside the mask and in the water
 * are still held at zero: those win. With pass2 off, held is not used.
 *
 * The first pass does not see held cells - not as constraints, not as
 * barriers - which is what they are for: a graded river's heights, a whisker
 * below the contour beside them, make the steepest pair in reach a nearly flat
 * one, and the ground beside the river fills flat to the edge of the reach.
 * Held instead, the river is a boundary the solve runs down to. It is the
 * caller's to leave those cells out of the constraints; this does not touch
 * them. pass1_out is the first pass as it was, without the held values.
 *
 * There is no nodata argument for held: a caller whose held raster has its
 * own nodata value writes ISOFILL_NO_ELEV or NaN in those cells first, as the
 * binary does for --held.
 *
 * isofill_run_ex is this with NULL. The binary's --held is this.
 */
long long isofill_run_held(const float *constraints, int has_nodata, double nodata,
                           const unsigned char *mask, const unsigned char *water,
                           const float *held, int cols, int rows,
                           const isofill_params *params, float *out, float *pass1_out);

/*
 * The fill's inputs beside the constraints, any of them NULL: what
 * isofill_run_layers takes, so that the next input is a field here rather
 * than another function.
 *
 *   mask, water   cols*rows bytes, as for isofill_run
 *   held          cols*rows floats, as for isofill_run_held
 *   wall          cols*rows bytes: where non-zero, a wall (--wall, 0.16.0) -
 *                 a cliff. The first pass does not see past it, as it does
 *                 not see past a constraint, and never answers it; the second
 *                 does not diffuse across it, a wall's neighbours leaving it
 *                 out of their mean. So the ground either side of a cliff is
 *                 filled from its own side. A wall cell takes the highest of
 *                 its neighbours once they are solved - the cliff's top, where
 *                 OSM draws one. A wall stops the void's flood, which goes
 *                 round a cliff in ground nobody described rather than along
 *                 it; a cell beside a wall that sees nothing past it, with
 *                 constraints in its radius, is declined rather than out of
 *                 reach (0.16.1). A wall cell that is a constraint is a
 *                 constraint.
 */
typedef struct {
    const unsigned char *mask;
    const unsigned char *water;
    const float *held;
    const unsigned char *wall;
} isofill_layers;

/*
 * The fill with all of its inputs. isofill_run_held is this with no wall.
 */
long long isofill_run_layers(const float *constraints, int has_nodata, double nodata,
                             const isofill_layers *layers, int cols, int rows,
                             const isofill_params *params, float *out, float *pass1_out);

/*
 * The held cells written into a surface the first pass left, as
 * isofill_run_held does between its passes: each cell within radius of a held
 * cell (the square of side 2 radius + 1) that the first pass called out of
 * reach becomes declined, and each held cell takes its value. held is
 * cols*rows floats, a cell held unless it is ISOFILL_NO_ELEV or NaN; surface
 * is written in place. isofill_diffuse after it is the rest of
 * isofill_run_held, cell for cell - for a caller that runs the two passes
 * apart, as an editor solving a box does.
 *
 * Returns 0, or -1 for bad arguments and -2 for out of memory.
 */
int isofill_hold(float *surface, const float *held, int cols, int rows, int radius);

/*
 * The second pass on its own, over a surface the first pass has already
 * written. surface is cols*rows floats holding what pass 1 left - sentinels
 * for the cells it declined - and is solved in place. mask and water are
 * cols*rows bytes or NULL and mean what they mean to isofill_run.
 *
 * There is no separate way to say "hold this cell": a cell carrying anything
 * that is not a sentinel is already fixed, and the solve never moves a fixed
 * cell. So a caller that cuts a box out of a larger surface writes the last
 * whole-raster answer into the box's rim, and the second pass runs up to it
 * instead of to the raster's edge.
 *
 * That is an approximation - diffusion is global, and a rim held at yesterday's
 * answer is a boundary the whole-raster solve does not have. How good one is
 * depends on how far a change travels, which is the caller's to measure.
 *
 * isofill_run's own pass 2 is this function's behaviour with the surface it
 * had just filled, which the library test pins.
 *
 * Returns 0, or -1 for bad arguments and -2 for out of memory - and means it:
 * the whole of the second pass reports a failed allocation rather than
 * printing and calling exit. The first pass is not converted yet, which is
 * said on isofill_run where a caller of that will see it; nothing here goes
 * near the first pass.
 */
int isofill_diffuse(float *surface, const unsigned char *mask,
                    const unsigned char *water, int cols, int rows);

/*
 * isofill_diffuse with walls (0.16.0): the second pass not diffusing across
 * them, as isofill_run_layers's does. isofill_diffuse is this with NULL.
 */
int isofill_diffuse_wall(float *surface, const unsigned char *mask,
                         const unsigned char *water, const unsigned char *wall,
                         int cols, int rows);

#endif
