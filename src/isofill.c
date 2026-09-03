/*
 * isofill - fill a raster from rasterised contours
 *
 * Two passes. The first is convert_tile_radius from
 * original/tile_util.c:
 *
 *   1. for every unset cell, gather the constraints within radius which are in
 *      line of sight - the straight line to them crossing no other constraint -
 *      take the pair with the steepest gradient across the cell, and interpolate
 *      between it. If no two differing elevations are in sight, leave the cell
 *      unset. That refusal is the point: it is why water enclosed by a coastline
 *      is never invented, and why a headland casts no plume.
 *
 *   2. interpolate every row and every column between its known values, anchored
 *      at zero one step beyond each end, and blend the two by whichever gradient
 *      is steeper. Cells pass 1 declined then take the interpolation between the
 *      constraints around them - between two coastlines, exactly zero.
 *
 * What is ours rather than the original's:
 *
 *   - only the nearest constraint at each distinct elevation is kept. For a
 *     fixed pair of elevations the steepest gradient uses the nearest of each,
 *     since a closer sample only ever increases (v2-v1)/(d2+d1), so this is
 *     exact and turns an O(k^2) search over every sample into one over the
 *     handful of distinct levels in view.
 *
 *   - a summed area table over the constraint mask, so a cell with nothing
 *     within radius is skipped in constant time. On a sparse zone - 22 drawn
 *     squares in 153 degrees, for one of ours - that is most of the raster.
 *
 *   - threads across the rows of a pass, and, above --max-mem, horizontal bands
 *     with a radius of margin so that memory follows the band rather than the
 *     raster. A cell only ever looks radius away, so the banded result is
 *     identical to the whole raster one, not an approximation of it. The
 *     summed area table is what forces this: eight bytes a cell, twenty
 *     gigabytes on a zone of 2.5 gigapixels. Pass 2 cannot band the same way,
 *     since its columns run the full height, so it goes to a temporary in
 *     column strips and is combined back in row chunks.
 *
 * SPDX-License-Identifier: Artistic-1.0 OR GPL-1.0-or-later
 * Copyright (c) 2017-2020 Thilo Stapff
 * Copyright (c) 2026 Lee Kindness and OpenGeofiction administrators
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gdal.h>
#include <cpl_conv.h>
#include <cpl_string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * An elevation, and the working type of both passes. Float since September
 * 2026: the DEM was whole metres, and the first pass rounded every value it
 * interpolated between two contours while the second rounded everything it
 * filled. A hillshade is a derivative, so that quantisation drew a hard ring
 * around every contour - the terracing - and no amount of smoothing downstream
 * put back what had already been thrown away.
 *
 * The sentinels keep their values rather than becoming NaN. They are exactly
 * representable, they compare exactly, and every test that read them still
 * reads them; NaN would have meant rewriting each one as an ordering question.
 */
typedef float elev;

#define NO_ELEV (-32768)

/*
 * A cell with no constraint at all within the radius. Distinct from NO_ELEV,
 * which is a cell the first pass looked at and declined: this one it never had
 * anything to look at.
 *
 * The second pass fills both alike, which is what the original does - TileUtil's
 * radius_value has a single failure return, and finding nothing falls through to
 * it exactly as finding too little does. Until August 2026 this one was reset to
 * zero afterwards, copying gdal_fillnodata -md rather than the implementation
 * this reproduces, and that put 114,309 cells of zone-tapira at 1 m.
 *
 * The distinction is still recorded, because it is worth being able to see:
 * --no-pass2 writes these sentinels through, and telling "saw nothing" from
 * "saw too little" is how that fault was found.
 */
#define OUT_OF_REACH (-32767)

/*
 * A cell that saw exactly one elevation - a contour, but nothing to interpolate
 * against. Distinct from NO_ELEV, which saw a pair and judged it too flat to
 * trust. The second pass treats all three as unset, so this changes no output;
 * it separates two populations that want different answers. Of the cells the
 * first pass declines inside a described area, which of these they are decides
 * whether a relaxed rule has anything to work with.
 */
#define ONE_LEVEL (-32766)

typedef struct {
    short dx, dy;
    float dist;
    int   first;   /* index into the shared point pool */
    int   n;       /* how many intermediate points */
} Offset;

typedef struct {
    short *px, *py;      /* the point pool, all offsets' lines end to end */
    Offset *off;
    int noff;
} Rays;

/* GDALRasterIO is warn_unused_result; where a failure is not actionable - a
 * write to a temporary we are about to close - say so rather than test it. */
#define IO_(call) do { if ((call) != CE_None) \
        fprintf(stderr, "isofill: raster io failed\n"); } while (0)

/* ------------------------------------------------------------------ rays */

/*
 * make_line_points from the original: the intermediate points between the
 * origin and (x0,y0), both ends excluded. Two cells per step, not one, so the
 * ray is a supercover path - a thin ray slips diagonally between two contour
 * cells and sight leaks through it.
 */
static int line_points(int x0, int y0, short *px, short *py, int cap)
{
    int n = 0, y, x, i;
    double dd;

    if (x0 == 0 && y0 == 0)
        return 0;

    if (abs(x0) > abs(y0)) {                 /* swap the axes */
        n = line_points(y0, x0, px, py, cap);
        for (i = 0; i < n; i++) {
            short t = px[i]; px[i] = py[i]; py[i] = t;
        }
        return n;
    }
    if (y0 < 0) {                            /* mirror north-south */
        n = line_points(x0, -y0, px, py, cap);
        for (i = 0; i < n; i++)
            py[i] = (short) -py[i];
        return n;
    }

    dd = (double) x0 / (double) y0;
    for (y = 0; y < y0; y++) {
        x = (int) floor(dd * (y + 0.5) + 0.5);
        if (!(n > 0 && px[n - 1] == x && py[n - 1] == y) && !(y == 0 && x == 0)) {
            if (n >= cap) return n;
            px[n] = (short) x; py[n] = (short) y; n++;
        }
        if (!(y + 1 == y0 && x == x0)) {
            if (n >= cap) return n;
            px[n] = (short) x; py[n] = (short) (y + 1); n++;
        }
    }
    return n;
}

static int cmp_off(const void *a, const void *b)
{
    const Offset *p = a, *q = b;
    return (p->dist > q->dist) - (p->dist < q->dist);
}

static Rays *rays_build(int radius)
{
    Rays *r = calloc(1, sizeof *r);
    int cap_pts = 4 * radius + 8, dx, dy, i, used = 0;
    short *tmpx = malloc(cap_pts * sizeof *tmpx);
    short *tmpy = malloc(cap_pts * sizeof *tmpy);
    size_t pool = 0;

    r->noff = 0;
    for (dy = -radius; dy <= radius; dy++)
        for (dx = -radius; dx <= radius; dx++)
            if (dx * dx + dy * dy > 0 && dx * dx + dy * dy <= radius * radius)
                r->noff++;
    r->off = malloc(r->noff * sizeof *r->off);

    i = 0;
    for (dy = -radius; dy <= radius; dy++)
        for (dx = -radius; dx <= radius; dx++) {
            if (dx * dx + dy * dy == 0 || dx * dx + dy * dy > radius * radius)
                continue;
            r->off[i].dx = (short) dx;
            r->off[i].dy = (short) dy;
            r->off[i].dist = (float) sqrt((double) (dx * dx + dy * dy));
            i++;
        }
    /* nearest first: the first sample seen at an elevation is then its nearest */
    qsort(r->off, r->noff, sizeof *r->off, cmp_off);

    for (i = 0; i < r->noff; i++)
        pool += 4 * radius + 8;
    r->px = malloc(pool * sizeof *r->px);
    r->py = malloc(pool * sizeof *r->py);

    for (i = 0; i < r->noff; i++) {
        int n = line_points(r->off[i].dx, r->off[i].dy, tmpx, tmpy, cap_pts);
        r->off[i].first = used;
        r->off[i].n = n;
        memcpy(r->px + used, tmpx, n * sizeof *tmpx);
        memcpy(r->py + used, tmpy, n * sizeof *tmpy);
        used += n;
    }
    free(tmpx); free(tmpy);
    return r;
}

static void rays_free(Rays *r)
{
    free(r->px); free(r->py); free(r->off); free(r);
}

/* ------------------------------------------------------------------ pass 1 */

#define MAX_LEVELS 64

typedef struct {
    int rows, cols;
    elev *v;
    unsigned char *is;    /* a constraint: carries a value */
    unsigned char *blk;   /* blocks sight: is, widened by --barrier */
    long long *sat;
} Band;

/*
 * Sight is blocked by the widened mask, not the constraint mask. At 3
 * arcseconds a one cell coastline is a 93 m wall; at 1 arcsecond the same one
 * cell is 31 m, so scaling the radius to hold the search distance left the
 * barrier three times thinner in ground terms and rays began threading gaps
 * which could not exist on the coarser grid. Widening the rasterised line
 * instead would write its elevation into the cells it gained, which is a
 * different and worse thing: what a cell is worth and what it hides are
 * separate questions.
 */
static inline int blocks(const Band *b, int x, int y)
{
    if (x < 0 || y < 0 || x >= b->cols || y >= b->rows)
        return 0;
    return b->blk[(size_t) y * b->cols + x];
}

static void dilate(Band *b, int n)
{
    unsigned char *tmp;
    int x, y, k;

    b->blk = malloc((size_t) b->cols * b->rows);
    if (!b->blk) { fprintf(stderr, "isofill: out of memory\n"); exit(1); }
    memcpy(b->blk, b->is, (size_t) b->cols * b->rows);
    if (n <= 0) return;

    tmp = malloc((size_t) b->cols * b->rows);
    if (!tmp) { fprintf(stderr, "isofill: out of memory\n"); exit(1); }
    for (k = 0; k < n; k++) {           /* n square dilations, separable */
        memcpy(tmp, b->blk, (size_t) b->cols * b->rows);
        for (y = 0; y < b->rows; y++)
            for (x = 0; x < b->cols; x++) {
                size_t i = (size_t) y * b->cols + x;
                if (tmp[i]) continue;
                if ((x > 0 && tmp[i - 1]) || (x + 1 < b->cols && tmp[i + 1])
                    || (y > 0 && tmp[i - b->cols])
                    || (y + 1 < b->rows && tmp[i + b->cols]))
                    b->blk[i] = 1;
            }
    }
    free(tmp);
}

static inline int have(const Band *b, int x, int y)
{
    return x >= 0 && x < b->cols && y >= 0 && y < b->rows && b->is[(size_t) y * b->cols + x];
}

/* constraints within the radius box, in constant time */
static inline long long sat_count(const Band *b, int x0, int y0, int x1, int y1)
{
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 >= b->cols) x1 = b->cols - 1;
    if (y1 >= b->rows) y1 = b->rows - 1;
    if (x1 < x0 || y1 < y0) return 0;
    {
        size_t w = b->cols + 1;
        return b->sat[(size_t)(y1 + 1) * w + (x1 + 1)]
             - b->sat[(size_t) y0      * w + (x1 + 1)]
             - b->sat[(size_t)(y1 + 1) * w +  x0]
             + b->sat[(size_t) y0      * w +  x0];
    }
}

static void build_sat(Band *b)
{
    size_t w = b->cols + 1;
    int x, y;
    b->sat = calloc(w * (b->rows + 1), sizeof *b->sat);
    for (y = 0; y < b->rows; y++)
        for (x = 0; x < b->cols; x++)
            b->sat[(size_t)(y + 1) * w + (x + 1)] =
                b->is[(size_t) y * b->cols + x]
                + b->sat[(size_t) y * w + (x + 1)]
                + b->sat[(size_t)(y + 1) * w + x]
                - b->sat[(size_t) y * w + x];
}

static int g_explain_x = -1, g_explain_y = -1;

static elev radius_value(const Band *b, const Rays *r, int radius,
                          double grad_min, int barrier, int x0, int y0,
                          unsigned char *blocked)
{
    int explain = (x0 == g_explain_x && y0 == g_explain_y);
    elev lev[MAX_LEVELS];
    float dst[MAX_LEVELS];
    int nlev = 0, i, j, k;
    double best = -1e30, out = 0;

    if (have(b, x0, y0))
        return b->v[(size_t) y0 * b->cols + x0];
    if (sat_count(b, x0 - radius, y0 - radius, x0 + radius, y0 + radius) == 0)
        return OUT_OF_REACH;

    for (i = 0; i < r->noff; i++) {
        const Offset *o = &r->off[i];
        int x = x0 + o->dx, y = y0 + o->dy, p, blk = 0;

        if (!have(b, x, y))
            continue;
        for (p = 0; p < o->n; p++) {
            int px = r->px[o->first + p], py = r->py[o->first + p];
            int adx = abs(px - o->dx), ady = abs(py - o->dy);
            /*
             * Near the target use the exact mask, so a constraint is never
             * hidden behind its own widening; everywhere else use the widened
             * one. Chebyshev, because the dilation is square.
             */
            int hit = ((adx > barrier || ady > barrier)
                       ? blocks(b, x0 + px, y0 + py)
                       : have(b, x0 + px, y0 + py));
            if (hit) {
                blk = 1;
                break;
            }
        }
        if (blk)
            continue;
        {
            elev val = b->v[(size_t) y * b->cols + x];
            /* nearest-first, so an elevation already seen is already nearer */
            for (k = 0; k < nlev; k++)
                if (lev[k] == val) break;
            if (k == nlev && nlev < MAX_LEVELS) {
                lev[nlev] = val;
                dst[nlev] = o->dist;
                nlev++;
                if (explain)
                    fprintf(stderr, "    sees %8.2f m at %6.2f cells, offset %+d%+d\n",
                            val, o->dist, o->dx, o->dy);
            }
        }
    }
    (void) blocked;

    for (i = 0; i < nlev; i++)
        for (j = 0; j < nlev; j++) {
            double g;
            if (lev[j] <= lev[i]) continue;
            g = (double)(lev[j] - lev[i]) / (dst[j] + dst[i]);
            if (g > best) {
                best = g;
                out = ((double) lev[i] * dst[j] + (double) lev[j] * dst[i])
                    / (dst[i] + dst[j]);
            }
        }
    if (explain)
        fprintf(stderr, "    %d distinct elevations in sight, steepest %.4f "
                "(floor %.4f) -> %s %.1f\n", nlev, best, grad_min,
                best <= grad_min ? "declined" : "value", out);
    /*
     * Nothing in sight at all, rather than nothing worth interpolating between.
     * The box test above is a square, so it lets through cells whose only
     * constraint is out at radius * sqrt(2); this is the exact answer, and it
     * costs nothing because the rays have already been walked.
     */
    if (nlev == 0)
        return OUT_OF_REACH;
    if (nlev == 1)
        return ONE_LEVEL;
    /*
     * The floor is metres per cell, so what it demands of a pair grows with how
     * far apart they are: at radius 60 it asks for 6 m of fall, at 120 for 12,
     * and on gentle ground that rejects pairs which are perfectly good evidence.
     * It was 0.1 to August 2026, from the original. Measured on the box inside
     * N32E067_Ellarca, dropping it to 0.02 takes the first pass from 32.71% of
     * the described area to 36.42% at radius 60, and from 55.86% to 74.82% at
     * radius 120 - where it was declining 18.95% of the area on this test alone.
     * Below 0.02 nothing further changes, at either radius.
     *
     * The floor wants to be a total drop rather than a rate, or scaled by the
     * pair's separation. That is the better fix and is not this one.
     */
    if (best <= grad_min)
        return NO_ELEV;
    return (elev) out;
}

/* ------------------------------------------------------------------ pass 2 */

static void diffuse(elev *v, const unsigned char *water,
                    const unsigned char *mask, int cols, int rows);

/*
 * Pass 2 the way the original ran it. It never saw more than one tile: 512
 * square with a margin of 128, both passes over the 768 window, and only the
 * central 512 kept - surroundTile and extractSubtile in ElevationTile.pm. That
 * bound is not incidental. Every row and column is anchored at zero beyond its
 * ends, so the length of the line decides how far a value carries, and a whole
 * zone is a row of 86401 cells where the original's was 768. Run unbounded, a
 * single filled cell in the water spreads down its whole row and column: on
 * alved, pass 1 leaves 2.30% of the water carrying elevation and pass 2 takes
 * that to 9.67%.
 *
 * The window is expressed in radii rather than cells so it holds at any
 * resolution - the original's 512 and 128 are 25.6 and 6.4 times its radius of
 * 20.
 */

/*
 * Pass 2 over a raster too large to hold. The column interpolation needs whole
 * columns and the row interpolation needs whole rows, so the columns go to two
 * temporaries in vertical strips, and the rows are then streamed against them.
 * Only the column results are stored - a row is interpolated in memory as it is
 * combined.
 */


/* ------------------------------------------------- pass 2, by diffusion */

/*
 * The linear pass above is continuous in height and not in slope. It
 * interpolates along each row and each column and blends the two, so the
 * surface is piecewise linear and its derivative jumps at every anchor. A
 * hillshade is a derivative, so each of those joints draws an edge, and the
 * edges run along rows and columns because the segments do - which is why
 * sparse ground comes out looking like a circuit diagram rather than a
 * landscape. Measured on a box inside N32E067_Ellarca, where the first pass
 * answers a third of the ground, the artefact covers most of the rest.
 *
 * This solves Laplace's equation instead: every cell the first pass answered is
 * held fixed, and the rest relax to the mean of their four neighbours. The
 * result is smooth in the derivative and has no preferred direction. It is also
 * the same surface to within a metre - on that box the median difference is 0 m
 * and the 95th percentile 3 m - because it is not inventing different terrain,
 * only removing the creases.
 *
 * A cell beyond the edge of the raster counts as zero, which is what the linear
 * pass does by anchoring every row and column one step past its ends, and is
 * what lets open water come out at zero where the sea runs off the side of a
 * zone. Water enclosed by its own coastline needs no such help: bounded by
 * cells at zero on all sides, Laplace gives zero throughout without being told.
 *
 * Solved coarse to fine. Relaxation moves information one cell per sweep, so
 * over a gap five hundred cells wide a fine-grid solve would take thousands of
 * sweeps; each level here starts from the level above and only has to add the
 * detail that level could not carry.
 */


#define DIFFUSE_SWEEPS 24        /* red and black pairs, at each level */
#define DIFFUSE_COARSE 200       /* at the coarsest, which is small and cheap */
#define DIFFUSE_MIN    24        /* stop coarsening below this on either side */

static int is_unset(elev v)
{
    return v == NO_ELEV || v == OUT_OF_REACH || v == ONE_LEVEL;
}

/*
 * One half of a red-black Gauss-Seidel sweep. Red and black cells never
 * neighbour each other, so each half is safely parallel and each reads values
 * the other half has already updated - which converges about twice as fast as
 * Jacobi and needs no second array.
 */
static void gs_half(float *z, const unsigned char *fx, int cols, int rows,
                    int parity, float relax_w)
{
    int y;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (y = 0; y < rows; y++) {
        size_t r = (size_t) y * cols;
        int x;
        for (x = (y + parity) & 1; x < cols; x += 2) {
            size_t k = r + (size_t) x;
            float s;
            if (fx[k]) continue;
            s  = (y > 0)        ? z[k - (size_t) cols] : 0.0f;
            s += (y < rows - 1) ? z[k + (size_t) cols] : 0.0f;
            s += (x > 0)        ? z[k - 1]             : 0.0f;
            s += (x < cols - 1) ? z[k + 1]             : 0.0f;
            z[k] = s * relax_w;
        }
    }
}

static void gs_sweeps(float *z, const unsigned char *fx, int cols, int rows,
                      int n, int level)
{
    const float w = 0.25f;
    int i;
    (void) level;
    for (i = 0; i < n; i++) {
        gs_half(z, fx, cols, rows, 0, w);
        gs_half(z, fx, cols, rows, 1, w);
    }
}

/*
 * Coarsen by two. A coarse cell is fixed only if one of the four under it was,
 * and takes the mean of those rather than of all four - so a contour a single
 * cell wide survives being coarsened instead of being averaged away by the
 * unknown ground around it.
 */
static void diffuse_restrict(const float *z, const unsigned char *fx, int cols, int rows,
                             float *cz, unsigned char *cfx, int ccols, int crows)
{
    int cy;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (cy = 0; cy < crows; cy++) {
        int cx;
        for (cx = 0; cx < ccols; cx++) {
            int y0 = cy * 2, x0 = cx * 2, dy, dx, nf = 0, na = 0;
            double sf = 0, sa = 0;
            for (dy = 0; dy < 2; dy++)
                for (dx = 0; dx < 2; dx++) {
                    int y = y0 + dy, x = x0 + dx;
                    size_t k;
                    if (y >= rows || x >= cols) continue;
                    k = (size_t) y * cols + (size_t) x;
                    sa += z[k]; na++;
                    if (fx[k]) { sf += z[k]; nf++; }
                }
            {
                size_t ck = (size_t) cy * ccols + (size_t) cx;
                cfx[ck] = nf ? 1 : 0;
                cz[ck] = (float) (nf ? sf / nf : (na ? sa / na : 0.0));
            }
        }
    }
}

/* Nearest-neighbour, which the sweeps that follow smooth out. */
static void diffuse_prolong(const float *cz, int ccols, int crows,
                            float *z, int cols, int rows)
{
    int y;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (y = 0; y < rows; y++) {
        int cy = y / 2, x;
        if (cy >= crows) cy = crows - 1;
        for (x = 0; x < cols; x++) {
            int cx = x / 2;
            if (cx >= ccols) cx = ccols - 1;
            z[(size_t) y * cols + (size_t) x] = cz[(size_t) cy * ccols + (size_t) cx];
        }
    }
}

#define DIFFUSE_LEVELS 24

/*
 * Which of the cells the first pass could not reach are void, and which are
 * pockets.
 *
 * Holding every OUT_OF_REACH cell at zero is too blunt. It suppresses the large
 * empty regions, which is wanted, and also punches a pit into every gap between
 * contours that happens to be wider than the radius, which is not: the fill
 * then ramps from the surrounding contour down to zero across sixty cells, a
 * gradient of a fifth, and hillshades as a dark blob with no relation to the
 * terrain. On zone-ellarca that was worse than the artefact it was fixing.
 *
 * What separates them is not the verdict but where the region reaches. The void
 * runs to the edge of the drawn area - ground the contours never described at
 * all; a pocket is enclosed by ground the first pass answered. That is the same
 * distinction demLandClamp.py draws between the sea and an inland lake, by
 * flooding from the edge and keeping what it reaches.
 *
 * Flooded by scanline sweeps rather than a queue: a queue over a region of a
 * hundred million cells wants more memory than the solve does, while forward
 * and backward sweeps carry a flag the length of a run at a time and converge
 * in a handful of passes on shapes this simple.
 */
static unsigned char *mark_void(const unsigned char *oor, const unsigned char *mask,
                                int cols, int rows)
{
    unsigned char *vd = calloc((size_t) cols * rows, 1);
    int y, x, pass, changed;
    if (!vd) { fprintf(stderr, "isofill: out of memory for pass 2\n"); exit(1); }

    /* seeds: outside the drawn area, and the raster edge, which the linear pass
     * anchors at zero one step beyond */
    for (y = 0; y < rows; y++)
        for (x = 0; x < cols; x++) {
            size_t k = (size_t) y * cols + (size_t) x;
            if (mask && !mask[k]) vd[k] = 1;
            else if (oor[k] &&
                     (y == 0 || x == 0 || y == rows - 1 || x == cols - 1)) vd[k] = 1;
        }

    for (pass = 0; pass < 64; pass++) {
        changed = 0;
        for (y = 0; y < rows; y++)
            for (x = 0; x < cols; x++) {
                size_t k = (size_t) y * cols + (size_t) x;
                if (vd[k] || !oor[k]) continue;
                if ((x && vd[k - 1]) || (y && vd[k - (size_t) cols]))
                    { vd[k] = 1; changed = 1; }
            }
        for (y = rows - 1; y >= 0; y--)
            for (x = cols - 1; x >= 0; x--) {
                size_t k = (size_t) y * cols + (size_t) x;
                if (vd[k] || !oor[k]) continue;
                if ((x < cols - 1 && vd[k + 1]) ||
                    (y < rows - 1 && vd[k + (size_t) cols]))
                    { vd[k] = 1; changed = 1; }
            }
        if (!changed) break;
    }
    return vd;
}

/*
 * Solve, given the finest level already classified: z holds the fixed values
 * and zero elsewhere, fx says which is which. Coarse to fine over a pyramid,
 * because relaxation moves information one cell per sweep.
 */
static void solve_pyramid(float *z0, unsigned char *fx0, int cols, int rows)
{
    float *z[DIFFUSE_LEVELS];
    unsigned char *fx[DIFFUSE_LEVELS];
    int lc[DIFFUSE_LEVELS], lr[DIFFUSE_LEVELS];
    int n = 1, i;

    z[0] = z0; fx[0] = fx0; lc[0] = cols; lr[0] = rows;
    while (n < DIFFUSE_LEVELS &&
           lc[n - 1] > DIFFUSE_MIN && lr[n - 1] > DIFFUSE_MIN) {
        int cc = (lc[n - 1] + 1) / 2, cr = (lr[n - 1] + 1) / 2;
        z[n] = malloc((size_t) cc * cr * sizeof *z[n]);
        fx[n] = malloc((size_t) cc * cr);
        if (!z[n] || !fx[n]) { fprintf(stderr, "isofill: out of memory for pass 2\n"); exit(1); }
        diffuse_restrict(z[n - 1], fx[n - 1], lc[n - 1], lr[n - 1],
                         z[n], fx[n], cc, cr);
        lc[n] = cc; lr[n] = cr;
        n++;
    }

    gs_sweeps(z[n - 1], fx[n - 1], lc[n - 1], lr[n - 1], DIFFUSE_COARSE, n - 1);
    for (i = n - 2; i >= 0; i--) {
        size_t m = (size_t) lc[i] * lr[i], j;
        float *fine = (i == 0) ? NULL : malloc(m * sizeof *fine);
        if (i && !fine) { fprintf(stderr, "isofill: out of memory for pass 2\n"); exit(1); }
        if (!fine) {
            /* the finest level is the caller's array. Prolong into the unknowns
             * in place: a whole extra copy of it, at four bytes a cell, is the
             * difference between holding a zone and not */
            int cy, cx;
            for (cy = 0; cy < lr[i]; cy++) {
                int py = cy / 2; if (py >= lr[i + 1]) py = lr[i + 1] - 1;
                for (cx = 0; cx < lc[i]; cx++) {
                    size_t k = (size_t) cy * lc[i] + (size_t) cx;
                    int px = cx / 2; if (px >= lc[i + 1]) px = lc[i + 1] - 1;
                    if (!fx[i][k])
                        z[i][k] = z[i + 1][(size_t) py * lc[i + 1] + (size_t) px];
                }
            }
            (void) j;
        } else {
            diffuse_prolong(z[i + 1], lc[i + 1], lr[i + 1], fine, lc[i], lr[i]);
            for (j = 0; j < m; j++) if (fx[i][j]) fine[j] = z[i][j];
            free(z[i]);
            z[i] = fine;
        }
        gs_sweeps(z[i], fx[i], lc[i], lr[i], DIFFUSE_SWEEPS, i);
    }
    for (i = 1; i < n; i++) { free(z[i]); free(fx[i]); }
}

/* Which cells are held, and at what. */
static void classify(const elev *v, const unsigned char *water,
                     const unsigned char *mask, const unsigned char *vd,
                     float *z, unsigned char *fx, size_t cells)
{
    size_t k;
    for (k = 0; k < cells; k++) {
        if ((mask && !mask[k]) || (vd && vd[k]) || (water && water[k])) {
            fx[k] = 1;
            z[k] = 0.0f;
        } else {
            int unset = is_unset(v[k]);
            fx[k] = (unsigned char) !unset;
            z[k] = unset ? 0.0f : (float) v[k];
        }
    }
}

/*
 * What diffuse needs, in megabytes, to hold a raster of this size: the caller's
 * Int16 buffer, a float and a byte a cell for the solve, and a third again for
 * the pyramid above it. The out-of-reach and void flags are gone by the time
 * the pyramid is built, so they do not set the peak.
 */
static double diffuse_mb(int cols, int rows)
{
    double cells = (double) cols * rows;
    return (cells * (2.0 + 4.0 + 1.0) * (1.0 + 1.0 / 3.0)) / (1024 * 1024);
}

static void diffuse(elev *v, const unsigned char *water,
                    const unsigned char *mask, int cols, int rows)
{
    size_t k, cells = (size_t) cols * rows;
    unsigned char *oor = calloc(cells, 1), *vd;
    float *z = malloc(cells * sizeof *z);
    unsigned char *fx = malloc(cells);

    if (!oor || !z || !fx) { fprintf(stderr, "isofill: out of memory for pass 2\n"); exit(1); }
    for (k = 0; k < cells; k++) oor[k] = (v[k] == OUT_OF_REACH);
    vd = mark_void(oor, mask, cols, rows);
    free(oor);

    classify(v, water, mask, vd, z, fx, cells);
    free(vd);                    /* wanted only while classifying */
    solve_pyramid(z, fx, cols, rows);
    for (k = 0; k < cells; k++) v[k] = (elev) z[k];

    free(z); free(fx);
}

/*
 * The same solve, for a raster too large to hold.
 *
 * Only the finest level is a problem: coarsening by two divides the cells by
 * four, so zone-axian's 2.1 gigapixels are 131 million three levels down, which
 * fits several times over. So the pyramid is entered part way up - the fine
 * rasters are read once, in row groups, and restricted straight into a coarse
 * grid held in memory - and the answer is carried back down band by band.
 *
 * The band pass is not a tiled solve. Each band starts from the coarse answer,
 * which already carries the shape of the whole zone, and only has to add the
 * detail the coarse grid could not hold; a margin of a few hundred rows either
 * side is far more than that correction travels, so the bands agree where they
 * meet. Solving each band from cold would seam, which is the mistake the linear
 * pass's own tiling makes.
 */

static void diffuse_ooc(GDALRasterBandH p1b, GDALRasterBandH mb,
                        GDALRasterBandH wb, GDALRasterBandH ob,
                        int cols, int rows, double max_mem, int margin)
{
    int shift = 0, cc, cr, band, y;
    size_t ccells;
    float *cz;
    unsigned char *cfx, *cvd, *coor, *cout_of_mask;
    elev *v;
    unsigned char *mbuf = NULL, *wbuf = NULL;
    double budget = max_mem * 1024.0 * 1024.0 * 0.45;

    /* how far up the pyramid we have to start for the grid to fit */
    for (shift = 1; shift < 8; shift++) {
        double n = ((double) cols / (1 << shift)) * ((double) rows / (1 << shift));
        if (n * 7.0 <= budget) break;
    }
    cc = (cols + (1 << shift) - 1) >> shift;
    cr = (rows + (1 << shift) - 1) >> shift;
    ccells = (size_t) cc * cr;
    fprintf(stderr, "  pass 2: coarse grid %dx%d, 1 in %d, then %d row bands\n",
            cc, cr, 1 << shift, 1 << shift);
    fprintf(stderr, "isofill: warning - --pass2 diffuse out of core is an "
            "approximation. On zone-ellarca, forced out of core and compared\n"
            "  against the same zone solved whole: 95%% of cells identical, "
            "95th percentile 0 m, but a tenth\n"
            "  of a percent differ by hundreds of metres. Give --max-mem the "
            "room to hold the raster if you can.\n");

    cz = calloc(ccells, sizeof *cz);
    cfx = calloc(ccells, 1);
    coor = malloc(ccells);
    cout_of_mask = malloc(ccells);
    if (!cz || !cfx || !coor || !cout_of_mask) {
        fprintf(stderr, "isofill: out of memory for pass 2\n"); exit(1);
    }
    memset(coor, 1, ccells);
    memset(cout_of_mask, 1, ccells);
    {
        float *acc = calloc(ccells, sizeof *acc);
        int *nfix = calloc(ccells, sizeof *nfix);
        int step = 1 << shift;
        if (!acc || !nfix) { fprintf(stderr, "isofill: out of memory\n"); exit(1); }
        v = malloc((size_t) cols * step * sizeof *v);
        if (mb) mbuf = malloc((size_t) cols * step);
        if (wb) wbuf = malloc((size_t) cols * step);
        if (!v || (mb && !mbuf) || (wb && !wbuf)) {
            fprintf(stderr, "isofill: out of memory\n"); exit(1);
        }
        for (y = 0; y < rows; y += step) {
            int h = (y + step <= rows) ? step : rows - y, yy;
            IO_(GDALRasterIO(p1b, GF_Read, 0, y, cols, h, v, cols, h, GDT_Float32, 0, 0));
            if (mb) IO_(GDALRasterIO(mb, GF_Read, 0, y, cols, h, mbuf, cols, h, GDT_Byte, 0, 0));
            if (wb) IO_(GDALRasterIO(wb, GF_Read, 0, y, cols, h, wbuf, cols, h, GDT_Byte, 0, 0));
            for (yy = 0; yy < h; yy++) {
                size_t r = (size_t) yy * cols, cr_off = (size_t) ((y + yy) >> shift) * cc;
                int x;
                for (x = 0; x < cols; x++) {
                    size_t k = r + (size_t) x, ck = cr_off + (size_t) (x >> shift);
                    if (!(mb && !mbuf[k])) cout_of_mask[ck] = 0;
                    if (v[k] != OUT_OF_REACH) coor[ck] = 0;
                    if (!is_unset(v[k]) && !(wb && wbuf[k]) && !(mb && !mbuf[k])) {
                        acc[ck] += (float) v[k];
                        nfix[ck]++;
                    }
                }
            }
        }
        for (ccells = (size_t) cc * cr, y = 0; (size_t) y < ccells; y++) {
            if (nfix[y]) { cfx[y] = 1; cz[y] = acc[y] / nfix[y]; }
        }
        free(acc); free(nfix);
    }
    for (y = 0; (size_t) y < ccells; y++) if (cout_of_mask[y]) coor[y] = 1;
    cvd = mark_void(coor, mb ? cout_of_mask : NULL, cc, cr);
    for (y = 0; (size_t) y < ccells; y++)
        if (cvd[y] || cout_of_mask[y]) { cfx[y] = 1; cz[y] = 0.0f; }
    free(coor); free(cout_of_mask);

    solve_pyramid(cz, cfx, cc, cr);
    free(cfx);

    /* carry it back down, a band of rows at a time */
    {
        int band_rows = (int) (budget / ((double) cols * 7.0));
        int done = 0;
        if (band_rows < 2 * margin) band_rows = 2 * margin;
        free(v); free(mbuf); free(wbuf);
        for (band = 0; band < rows; band += band_rows) {
            int h = (band + band_rows <= rows) ? band_rows : rows - band;
            int y0 = band - margin, y1 = band + h + margin, bh, oy, x;
            size_t n;
            float *z; unsigned char *fx, *vd;
            if (y0 < 0) y0 = 0;
            if (y1 > rows) y1 = rows;
            bh = y1 - y0; oy = band - y0; n = (size_t) cols * bh;

            v = malloc(n * sizeof *v);
            z = malloc(n * sizeof *z);
            fx = malloc(n);
            vd = malloc(n);
            mbuf = mb ? malloc(n) : NULL;
            wbuf = wb ? malloc(n) : NULL;
            if (!v || !z || !fx || !vd || (mb && !mbuf) || (wb && !wbuf)) {
                fprintf(stderr, "isofill: out of memory for a pass 2 band\n"); exit(1);
            }
            IO_(GDALRasterIO(p1b, GF_Read, 0, y0, cols, bh, v, cols, bh, GDT_Float32, 0, 0));
            if (mb) IO_(GDALRasterIO(mb, GF_Read, 0, y0, cols, bh, mbuf, cols, bh, GDT_Byte, 0, 0));
            if (wb) IO_(GDALRasterIO(wb, GF_Read, 0, y0, cols, bh, wbuf, cols, bh, GDT_Byte, 0, 0));

            /* the void is decided on the coarse grid, where the whole zone is
             * visible; a band on its own cannot see where a region reaches */
            for (y = 0; y < bh; y++)
                for (x = 0; x < cols; x++) {
                    size_t k = (size_t) y * cols + (size_t) x;
                    size_t ck = (size_t) ((y0 + y) >> shift) * cc + (size_t) (x >> shift);
                    vd[k] = (unsigned char) (v[k] == OUT_OF_REACH && cvd[ck]);
                }
            classify(v, wbuf, mbuf, vd, z, fx, n);
            /*
             * The band is solved properly, not merely relaxed. Sweeping alone
             * carries information one cell at a time, and in sparse ground the
             * contours are hundreds of cells apart, so the band would keep
             * whatever the coarse grid gave it and the joins would show.
             *
             * What makes a band solvable on its own is a boundary, and the
             * coarse answer is one: the rows in the margin are held at it, top
             * and bottom, so each band is a complete problem whose edges agree
             * with its neighbours by construction. The margin is then thrown
             * away and only the interior kept.
             */
            for (y = 0; y < bh; y++) {
                int in_margin = (y0 + y < band && y0 > 0) ||
                                (y0 + y >= band + h && y1 < rows);
                for (x = 0; x < cols; x++) {
                    size_t k = (size_t) y * cols + (size_t) x;
                    size_t ck = (size_t) ((y0 + y) >> shift) * cc + (size_t) (x >> shift);
                    if (in_margin) { fx[k] = 1; z[k] = cz[ck]; }
                    else if (!fx[k]) z[k] = cz[ck];
                }
            }
            solve_pyramid(z, fx, cols, bh);

            for (y = 0; y < h; y++)
                for (x = 0; x < cols; x++) {
                    size_t k = (size_t) (oy + y) * cols + (size_t) x;
                    v[k] = (elev) z[k];
                }
            IO_(GDALRasterIO(ob, GF_Write, 0, band, cols, h,
                             v + (size_t) oy * cols, cols, h, GDT_Float32, 0, 0));
            free(v); free(z); free(fx); free(vd); free(mbuf); free(wbuf);
            done += h;
            fprintf(stderr, "\r  pass 2 %d%%", (int) (100.0 * done / rows));
        }
        fprintf(stderr, "\r  pass 2 complete, %d bands with %d margin\n",
                (rows + band_rows - 1) / band_rows, margin);
    }
    free(cz); free(cvd);
}

/* ------------------------------------------------------------------ main */

static void usage(void)
{
    fprintf(stderr,
        "usage: isofill [options] <constraints.tif> <out.tif>\n"
        "  --radius N     search radius in cells (default 20)\n"
        "  --barrier N    widen constraints by N cells for the sight test only,\n"
        "                 not for their values (default 1). At 1 arcsecond a one\n"
        "                 cell line is a third of the wall it was at 3\n"
        "  --grad-min F   least gradient, metres per cell, which counts as a\n"
        "                 slope worth interpolating across (default 0.02)\n"
        "  --no-pass2     leave cells the first pass declined unset\n"
        "  --water FILE   a raster the size of the constraints: where it is not\n"
        "                 zero the cell is water, held at zero by the second\n"
        "                 pass so the fill cannot cross it\n"
        "  --no-reach     accepted and ignored. It used to switch on what is now\n"
        "                 the only behaviour; see the note on reach in README.md\n"
        "  --mask FILE    a raster the size of the constraints: where it is zero\n"
        "                 no cell is filled, though contours there still count as\n"
        "                 evidence for cells outside it\n"
        "  --pass2-margin N   rows of overlap between the second pass's bands\n"
        "                 when the raster will not fit in --max-mem, each\n"
        "                 discarded after use (default 6 * radius)\n"
        "  --max-mem MB   above this, work band by band through a temporary\n"
        "                 beside the output (default 4096)\n"
        "  --threads N    (default: all but two, to leave the box usable)\n");
    exit(2);
}

/*
 * Pass 1 over a band of rows. The band carries a margin of radius above and
 * below, because a cell only ever looks that far, so rows [margin, margin+h)
 * of the band see exactly what they would in the whole raster. The summed area
 * table is the reason this matters: it is eight bytes a cell, and on a zone of
 * 2.5 gigapixels that alone is twenty gigabytes.
 */
static long long pass1_band(GDALRasterBandH sb, GDALRasterBandH mb,
                            const Rays *rays, int radius,
                            double grad_min, int barrier, int has_nd, double nd,
                            int cols, int rows, int y0, int h,
                            elev *out)
{
    unsigned char *mask = NULL;
    int top = y0 - radius, bot = y0 + h + radius;
    int bh, margin, y;
    long long filled = 0;
    Band b;

    if (top < 0) top = 0;
    if (bot > rows) bot = rows;
    bh = bot - top;
    margin = y0 - top;

    b.cols = cols; b.rows = bh;
    b.v = malloc((size_t) cols * bh * sizeof *b.v);
    b.is = malloc((size_t) cols * bh);
    if (!b.v || !b.is) { fprintf(stderr, "isofill: out of memory\n"); exit(1); }
    if (GDALRasterIO(sb, GF_Read, 0, top, cols, bh, b.v, cols, bh,
                     GDT_Float32, 0, 0) != CE_None) {
        fprintf(stderr, "isofill: read failed\n"); exit(1);
    }
    {
        size_t k, n = (size_t) cols * bh;
        for (k = 0; k < n; k++) {
            int isc = has_nd ? (b.v[k] != (elev) nd) : (b.v[k] != NO_ELEV);
            b.is[k] = (unsigned char) isc;
            if (!isc) b.v[k] = NO_ELEV;
        }
    }
    build_sat(&b);
    dilate(&b, barrier);

    /*
     * The mask says where a cell may be filled, not which constraints count -
     * a contour just outside it is still evidence for a cell just inside, so
     * the band is read whole and only the writing is restricted.
     */
    if (mb) {
        mask = malloc((size_t) cols * h);
        if (!mask) { fprintf(stderr, "isofill: out of memory\n"); exit(1); }
        IO_(GDALRasterIO(mb, GF_Read, 0, y0, cols, h, mask, cols, h,
                         GDT_Byte, 0, 0));
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) reduction(+:filled)
#endif
    for (y = 0; y < h; y++) {
        int xx;
        for (xx = 0; xx < cols; xx++) {
            size_t k = (size_t) y * cols + xx;
            elev v;
            if (mask && !mask[k]) {
                /*
                 * Outside the mask, invent nothing - but a constraint drawn
                 * there is still a constraint, and pass 2 interpolates along
                 * whole rows and columns. Drop it and a coastline just beyond
                 * the edge stops anchoring the cells inside.
                 */
                out[k] = have(&b, xx, y + margin)
                    ? b.v[(size_t) (y + margin) * cols + xx] : OUT_OF_REACH;
                continue;
            }
            v = radius_value(&b, rays, radius, grad_min, barrier, xx, y + margin, NULL);
            out[k] = v;
            if (v != NO_ELEV && v != OUT_OF_REACH && v != ONE_LEVEL) filled++;
        }
    }

    free(mask);
    free(b.v); free(b.is); free(b.blk); free(b.sat);
    return filled;
}

int main(int argc, char **argv)
{
    const char *in_path = NULL, *out_path = NULL, *mask_path = NULL;
    const char *water_path = NULL;
    int radius = 20, do_pass2 = 1, threads = 0, i;
    double grad_min = 0.02, max_mem = 4096;
    int p2_margin = -1, barrier = 1;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--radius") && i + 1 < argc) radius = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--grad-min") && i + 1 < argc) grad_min = atof(argv[++i]);
        else if (!strcmp(argv[i], "--max-mem") && i + 1 < argc) max_mem = atof(argv[++i]);
        else if (!strcmp(argv[i], "--mask") && i + 1 < argc) mask_path = argv[++i];
        else if (!strcmp(argv[i], "--water") && i + 1 < argc) water_path = argv[++i];
        else if (!strcmp(argv[i], "--barrier") && i + 1 < argc) barrier = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no-reach")) ;   /* was a switch, now the only behaviour */
        else if (!strcmp(argv[i], "--pass2-margin") && i + 1 < argc) p2_margin = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no-pass2")) do_pass2 = 0;
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) threads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--explain") && i + 2 < argc) {
            g_explain_x = atoi(argv[++i]); g_explain_y = atoi(argv[++i]);
            threads = 1;
        }
        else if (argv[i][0] == '-') usage();
        else if (!in_path) in_path = argv[i];
        else if (!out_path) out_path = argv[i];
        else usage();
    }
    if (!in_path || !out_path) usage();
#ifdef _OPENMP
    /*
     * All but two cores by default. util renders tiles and serves the wiki
     * while this runs, and a fill which takes every core starves them - the
     * headroom matters more than the last 30% of the speedup.
     */
    if (threads <= 0) {
        threads = omp_get_num_procs() - 2;
        if (threads < 1) threads = 1;
    }
    omp_set_num_threads(threads);
    fprintf(stderr, "  %d of %d cores\n", threads, omp_get_num_procs());
#else
    (void) threads;
#endif

    GDALAllRegister();
    {
        GDALDatasetH src = GDALOpen(in_path, GA_ReadOnly);
        GDALDatasetH msk = NULL, wat = NULL;
        GDALRasterBandH sb, mb = NULL, wb = NULL;
        int cols, rows, has_nd = 0;
        double nd, gt[6], whole_mb;
        Rays *rays;
        GDALDatasetH dst;
        char **opts = NULL;
        long long filled = 0;

        if (!src) { fprintf(stderr, "isofill: cannot open %s\n", in_path); return 1; }
        sb = GDALGetRasterBand(src, 1);
        cols = GDALGetRasterXSize(src);
        rows = GDALGetRasterYSize(src);
        nd = GDALGetRasterNoDataValue(sb, &has_nd);
        GDALGetGeoTransform(src, gt);
        if (mask_path) {
            msk = GDALOpen(mask_path, GA_ReadOnly);
            if (!msk) {
                fprintf(stderr, "isofill: cannot open %s\n", mask_path); return 1;
            }
            if (GDALGetRasterXSize(msk) != cols || GDALGetRasterYSize(msk) != rows) {
                fprintf(stderr, "isofill: mask is %dx%d, constraints are %dx%d\n",
                        GDALGetRasterXSize(msk), GDALGetRasterYSize(msk), cols, rows);
                return 1;
            }
            mb = GDALGetRasterBand(msk, 1);
        }
        if (water_path) {
            wat = GDALOpen(water_path, GA_ReadOnly);
            if (!wat) {
                fprintf(stderr, "isofill: cannot open %s\n", water_path); return 1;
            }
            if (GDALGetRasterXSize(wat) != cols || GDALGetRasterYSize(wat) != rows) {
                fprintf(stderr, "isofill: water is %dx%d, constraints are %dx%d\n",
                        GDALGetRasterXSize(wat), GDALGetRasterYSize(wat), cols, rows);
                return 1;
            }
            wb = GDALGetRasterBand(wat, 1);
        }
        rays = rays_build(radius);

        /*
         * Whole raster by default. The original bounded this to a 512 tile with
         * a 128 margin, and reproducing that was worth measuring, but it buys
         * nothing here: on alved it moves the water carrying 10 m or more not at
         * all, and shifts the land by a median of 1 m for the trouble. Left as an
         * option because it is the only thing which bounds how far a value
         * carries when there is no mask to bound it.
         */
        if (p2_margin < 0) p2_margin = 6 * radius;

        /* values, mask, output and the summed area table, in megabytes */
        whole_mb = ((double) cols * rows * 5.0
                    + (double) (cols + 1) * (rows + 1) * 8.0) / (1024 * 1024);

        opts = CSLSetNameValue(opts, "TILED", "YES");
        opts = CSLSetNameValue(opts, "COMPRESS", "DEFLATE");
        /*
         * 3, the floating point predictor, not 2. DEFLATE on raw float bytes
         * compresses badly - the exponent and mantissa interleave into
         * something close to noise - and zone-ellarca's DEM came out at 116 MB
         * against 11 for the Int16 it replaced, ten times rather than the two
         * the extra bytes account for. The float predictor separates the bytes
         * of each value before deflating, which is what the type needs.
         */
        opts = CSLSetNameValue(opts, "PREDICTOR", "3");
        opts = CSLSetNameValue(opts, "BIGTIFF", "IF_SAFER");

        fprintf(stderr, "  %dx%d, radius %d, %d rays, %.0f MB whole (limit %.0f)\n",
                cols, rows, radius, rays->noff, whole_mb, max_mem);

        if (whole_mb <= max_mem) {
            /* small enough to hold: one band, one write, pass 2 in memory */
            elev *out = malloc((size_t) cols * rows * sizeof *out);
            if (!out) { fprintf(stderr, "isofill: out of memory\n"); return 1; }
            filled = pass1_band(sb, mb, rays, radius, grad_min, barrier, has_nd, nd,
                                cols, rows, 0, rows, out);
            fprintf(stderr, "  pass 1 set %lld of %lld cells\n", filled,
                    (long long) cols * rows);
            if (do_pass2) {
                {
                    unsigned char *wbuf = NULL, *mbuf = NULL;
                    if (wb) {
                        wbuf = malloc((size_t) cols * rows);
                        if (!wbuf) { fprintf(stderr, "isofill: out of memory\n"); return 1; }
                        IO_(GDALRasterIO(wb, GF_Read, 0, 0, cols, rows, wbuf, cols, rows,
                                         GDT_Byte, 0, 0));
                    }
                    if (mb) {
                        mbuf = malloc((size_t) cols * rows);
                        if (!mbuf) { fprintf(stderr, "isofill: out of memory\n"); return 1; }
                        IO_(GDALRasterIO(mb, GF_Read, 0, 0, cols, rows, mbuf, cols, rows,
                                         GDT_Byte, 0, 0));
                    }
                    diffuse(out, wbuf, mbuf, cols, rows);
                    free(wbuf); free(mbuf);
                }
                if (mb) {
                    unsigned char *m = malloc((size_t) cols * rows);
                    size_t k, n = (size_t) cols * rows;
                    if (!m) { fprintf(stderr, "isofill: out of memory\n"); return 1; }
                    IO_(GDALRasterIO(mb, GF_Read, 0, 0, cols, rows, m, cols, rows,
                                     GDT_Byte, 0, 0));
                    for (k = 0; k < n; k++) if (!m[k]) out[k] = 0;
                    free(m);
                }
                fprintf(stderr, "  pass 2 complete\n");
            }
            dst = GDALCreate(GDALGetDriverByName("GTiff"), out_path, cols, rows, 1,
                             GDT_Float32, opts);
            GDALSetGeoTransform(dst, gt);
            GDALSetProjection(dst, GDALGetProjectionRef(src));
            if (!do_pass2)
                GDALSetRasterNoDataValue(GDALGetRasterBand(dst, 1), NO_ELEV);
            if (GDALRasterIO(GDALGetRasterBand(dst, 1), GF_Write, 0, 0, cols, rows,
                             out, cols, rows, GDT_Float32, 0, 0) != CE_None)
                fprintf(stderr, "isofill: write failed\n");
            GDALClose(dst);
            free(out);
        } else {
            /*
             * Too large. Pass 1 runs band by band into a temporary, and pass 2
             * reads that back a strip and a row at a time. The temporaries sit
             * beside the output, so the output directory needs the room.
             */
            char t1[4096], t2[4096], t3[4096];
            GDALDatasetH p1;
            int h, y0, band_rows;
            elev *out;
            double per_row;

            snprintf(t1, sizeof t1, "%s.pass1.tif", out_path);
            snprintf(t2, sizeof t2, "%s.pv.tif", out_path);
            snprintf(t3, sizeof t3, "%s.gv.tif", out_path);

            /* rows that fit: values, mask and output at 5 bytes, the table at 8 */
            per_row = ((double) cols * 5.0 + (double) (cols + 1) * 8.0) / (1024 * 1024);
            band_rows = (int) (max_mem / per_row) - 2 * radius;
            if (band_rows < 64) band_rows = 64;
            if (band_rows > rows) band_rows = rows;
            fprintf(stderr, "  banding: %d rows at a time, %d bands\n",
                    band_rows, (rows + band_rows - 1) / band_rows);

            p1 = GDALCreate(GDALGetDriverByName("GTiff"), t1, cols, rows, 1,
                            GDT_Float32, opts);
            if (!p1) { fprintf(stderr, "isofill: cannot create %s\n", t1); return 1; }
            out = malloc((size_t) cols * band_rows * sizeof *out);
            if (!out) { fprintf(stderr, "isofill: out of memory\n"); return 1; }

            for (y0 = 0; y0 < rows; y0 += band_rows) {
                h = (y0 + band_rows <= rows) ? band_rows : rows - y0;
                filled += pass1_band(sb, mb, rays, radius, grad_min, barrier, has_nd, nd,
                                     cols, rows, y0, h, out);
                IO_(GDALRasterIO(GDALGetRasterBand(p1, 1), GF_Write, 0, y0, cols, h,
                                 out, cols, h, GDT_Float32, 0, 0));
                fprintf(stderr, "\r  pass 1 %d%%", (int) (100.0 * (y0 + h) / rows));
            }
            free(out);
            fprintf(stderr, "\r  pass 1 set %lld of %lld cells\n", filled,
                    (long long) cols * rows);

            dst = GDALCreate(GDALGetDriverByName("GTiff"), out_path, cols, rows, 1,
                             GDT_Float32, opts);
            GDALSetGeoTransform(dst, gt);
            GDALSetProjection(dst, GDALGetProjectionRef(src));
            if (do_pass2) {
                /* a multiple of the 256 tile, so a strip write lands on
                 * whole tiles rather than straddling them */
                int strip = (int) ((max_mem / 3) * 1024 * 1024 / (rows * 3.0));
                int chunk = (int) ((max_mem / 5) * 1024 * 1024 / (cols * 9.0));
                strip = (strip / 256) * 256;
                if (strip < 256) strip = 256;
                if (strip > cols) strip = cols;
                chunk = (chunk / 256) * 256;
                if (chunk < 256) chunk = 256;
                if (chunk > rows) chunk = rows;
                fprintf(stderr, "  pass 2: %d column strips, %d row chunks\n",
                        (cols + strip - 1) / strip, (rows + chunk - 1) / chunk);
                if (diffuse_mb(cols, rows) <= max_mem) {
                    /*
                     * The first pass banded and the second does not have to.
                     * They are separate decisions because the two paths are not
                     * of equal standing: the first pass reads each band with a
                     * radius margin, so every interior cell sees its whole
                     * search circle and a banded run is exact. The second pass
                     * banded is an approximation. So where there is room for
                     * the solve but not for the first pass's summed-area table
                     * - which is most of the gap, at eight bytes a cell - take
                     * it, and keep the exact answer.
                     */
                    elev *whole = malloc((size_t) cols * rows * sizeof *whole);
                    unsigned char *wm = NULL, *mm = NULL;
                    if (!whole) { fprintf(stderr, "isofill: out of memory\n"); return 1; }
                    fprintf(stderr, "  pass 2: whole raster, %.0f MB of %.0f\n",
                            diffuse_mb(cols, rows), max_mem);
                    IO_(GDALRasterIO(GDALGetRasterBand(p1, 1), GF_Read, 0, 0,
                                     cols, rows, whole, cols, rows, GDT_Float32, 0, 0));
                    if (wb) {
                        wm = malloc((size_t) cols * rows);
                        if (!wm) { fprintf(stderr, "isofill: out of memory\n"); return 1; }
                        IO_(GDALRasterIO(wb, GF_Read, 0, 0, cols, rows, wm, cols, rows,
                                         GDT_Byte, 0, 0));
                    }
                    if (mb) {
                        mm = malloc((size_t) cols * rows);
                        if (!mm) { fprintf(stderr, "isofill: out of memory\n"); return 1; }
                        IO_(GDALRasterIO(mb, GF_Read, 0, 0, cols, rows, mm, cols, rows,
                                         GDT_Byte, 0, 0));
                    }
                    diffuse(whole, wm, mm, cols, rows);
                    if (mm) {
                        size_t k, nn = (size_t) cols * rows;
                        for (k = 0; k < nn; k++) if (!mm[k]) whole[k] = 0;
                    }
                    IO_(GDALRasterIO(GDALGetRasterBand(dst, 1), GF_Write, 0, 0,
                                     cols, rows, whole, cols, rows, GDT_Float32, 0, 0));
                    free(whole); free(wm); free(mm);
                    fprintf(stderr, "  pass 2 complete\n");
                } else {
                    diffuse_ooc(GDALGetRasterBand(p1, 1), mb, wb,
                                GDALGetRasterBand(dst, 1), cols, rows, max_mem,
                                p2_margin);
                }
            } else {
                elev *row = malloc((size_t) cols * sizeof *row);
                int y;
                GDALSetRasterNoDataValue(GDALGetRasterBand(dst, 1), NO_ELEV);
                for (y = 0; y < rows; y++) {
                    IO_(GDALRasterIO(GDALGetRasterBand(p1, 1), GF_Read, 0, y, cols, 1,
                                     row, cols, 1, GDT_Float32, 0, 0));
                    IO_(GDALRasterIO(GDALGetRasterBand(dst, 1), GF_Write, 0, y, cols, 1,
                                     row, cols, 1, GDT_Float32, 0, 0));
                }
                free(row);
            }
            GDALClose(dst);
            GDALClose(p1);
            VSIUnlink(t1); VSIUnlink(t2); VSIUnlink(t3);
        }

        if (msk) GDALClose(msk);
        if (wat) GDALClose(wat);
        GDALClose(src);
        CSLDestroy(opts);
        rays_free(rays);
    }
    return 0;
}
