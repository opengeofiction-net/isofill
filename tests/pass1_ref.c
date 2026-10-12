/*
 * The first pass against the walk it replaced. radius_value used to walk every
 * offset of the search circle nearest first, a line of sight to each
 * constraint; it now meets the constraints a row at a time and walks to only
 * those that can change the answer. That has to be the same answer to the bit,
 * on every cell: this keeps the old walk, as radius_value_ref, and compares
 * the two over random rasters - contours as random walks, spot heights, enough
 * distinct levels that a cell sees more than MAX_LEVELS of them, radius 20 and
 * 60, barrier 0 to 2.
 *
 * It includes isofill.c rather than linking the library, because what it
 * compares are the static functions inside it.
 */
#define ISOFILL_NO_MAIN
#include "../src/isofill.c"

static int ref_capped = 0;   /* set from every thread, only ever to 1 */

static elev radius_value_ref(const Band *b, const Rays *r, int radius,
                              double grad_min, int barrier, int x0, int y0)
{
    int explain = 0;
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
            for (k = 0; k < nlev; k++)
                if (lev[k] == val) break;
            if (k == nlev && nlev == MAX_LEVELS) ref_capped = 1;
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
    if (nlev == 0)
        return OUT_OF_REACH;
    if (nlev == 1)
        return ONE_LEVEL;
    if (best <= grad_min)
        return NO_ELEV;
    return (elev) out;
}


static unsigned long long seed = 88172645463325252ULL;
static unsigned rnd(void)
{
    seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
    return (unsigned) (seed >> 32);
}

/* a raster of constraints, the rest NO_ELEV */
static void make(elev *v, int cols, int rows, int lines, int levels, int spots)
{
    int i, s, k;
    for (k = 0; k < cols * rows; k++) v[k] = NO_ELEV;
    for (i = 0; i < lines; i++) {
        int x = rnd() % cols, y = rnd() % rows, len = 20 + rnd() % (2 * cols);
        elev lv = (elev) (10 * (rnd() % levels));
        for (s = 0; s < len; s++) {
            v[y * cols + x] = lv;
            switch (rnd() % 8) {          /* a wandering line, mostly onward */
            case 0: case 1: case 2: x++; break;
            case 3: y++; break;
            case 4: y--; break;
            case 5: x++; y++; break;
            case 6: x++; y--; break;
            default: x--; break;
            }
            if (x < 0 || x >= cols || y < 0 || y >= rows) break;
        }
    }
    for (i = 0; i < spots; i++)          /* spot heights: each its own level */
        v[(rnd() % rows) * cols + rnd() % cols] = (elev) (rnd() % 100000) / 7.0f;
    /* with no lines, spot heights alone: gentle ones, each its own level, and
     * one in fifty tall. Which levels a cell keeps once it has MAX_LEVELS in
     * sight then decides its answer - a tall one kept that the nearest-first
     * walk would have dropped is the steepest pair */
    if (!lines && levels == 0)
        for (i = 0; i < cols * rows; i++)
            if (v[i] != NO_ELEV)
                v[i] = (rnd() % 50) ? 100.0f + (elev) i / 1024.0f : 5000.0f + (elev) i;
}

static int compare(int cols, int rows, int radius, int barrier, int lines, int levels, int spots,
                   long long *cells, long long *valued)
{
    elev *v = malloc((size_t) cols * rows * sizeof *v);
    Rays *r = rays_build(radius);
    Band b;
    int x, y;
    make(v, cols, rows, lines, levels, spots);
    b.cols = cols; b.rows = rows; b.v = v;
    b.run = malloc((size_t) cols * rows);
    for (y = 0; y < rows; y++) {
        int next = cols;
        for (x = cols - 1; x >= 0; x--) {
            size_t k = (size_t) y * cols + x;
            if (v[k] != NO_ELEV) next = x;
            b.run[k] = (unsigned char) (next - x > 255 ? 255 : next - x);
        }
    }
    build_sat(&b);
    dilate(&b, barrier);
    {
        int bad = 0;
        long long nc = 0, nv = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) reduction(+:nc, nv) reduction(|:bad)
#endif
        for (y = 0; y < rows; y++) {
            int xx;
            for (xx = 0; xx < cols; xx++) {
                elev a = radius_value(&b, r, radius, 0.02, barrier, xx, y);
                elev e = radius_value_ref(&b, r, radius, 0.02, barrier, xx, y);
                if (memcmp(&a, &e, sizeof a) != 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
                    fprintf(stderr, "pass1_ref: radius %d barrier %d, cell %d,%d: %.9g, the walk gives %.9g\n",
                            radius, barrier, xx, y, a, e);
                    bad = 1;
                }
                nc++;
                if (a != NO_ELEV && a != OUT_OF_REACH && a != ONE_LEVEL) nv++;
            }
        }
        *cells += nc; *valued += nv;
        if (bad) return 1;
    }
    free(v); free(b.run); free(b.blk); free(b.sat); rays_free(r);
    return 0;
}

int main(void)
{
    long long cells = 0, valued = 0;
    int barrier;
    for (barrier = 0; barrier <= 2; barrier++) {
        /* sparse, dense, spot heights among contours, and spot heights alone
         * at one in a hundred - few enough to block little, so a cell sees
         * more distinct levels than MAX_LEVELS holds */
        if (compare(240, 160, 20, barrier, 30, 30, 0, &cells, &valued)) return 1;
        if (compare(240, 160, 20, barrier, 300, 30, 0, &cells, &valued)) return 1;
        if (compare(160, 120, 20, barrier, 40, 300, 2000, &cells, &valued)) return 1;
        if (compare(180, 150, 60, barrier, 60, 50, 0, &cells, &valued)) return 1;
        if (compare(160, 140, 60, barrier, 0, 1, 220, &cells, &valued)) return 1;
        if (compare(160, 140, 60, barrier, 0, 0, 260, &cells, &valued)) return 1;
    }
    if (!ref_capped) {
        fprintf(stderr, "pass1_ref: no cell saw more than %d levels, so the cap went untested\n", MAX_LEVELS);
        return 1;
    }
    printf("first pass: %lld cells, %lld valued, the same to the bit as the walk it replaced, "
           "the level cap among them\n", cells, valued);
    return 0;
}
