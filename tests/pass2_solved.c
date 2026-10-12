/*
 * The second pass solves its equation, not merely relaxes towards it.
 *
 * Laplace's answer is a weighted average of a region's edge, so it can never
 * leave the range of that edge: ground bounded by zero all round is zero, and
 * a ring of contours at one height holds a flat floor. Three cases, each
 * asking that of isofill_diffuse:
 *
 *   - the strip between a coast drawn as a 0 contour and the edge of the drawn
 *     area, with high ground answered just inland. On zone-ellarca the coarse
 *     levels saw the coast and the land as one, and the 24 sweeps at the
 *     finest level left the strip at up to 329 m.
 *   - a hilltop: a closed ring at 2350 m with a wide unanswered floor.
 *   - random ground, against a solve in double precision taken to 1e-9.
 *
 * And the same answer to the bit on one thread as on all of them, which the
 * golden surface in danu relies on: the solve is red-black, and its line
 * search sums tile by tile and adds the tiles in order.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "isofill.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define NO_ELEV ((float) ISOFILL_NO_ELEV)

static int failed = 0;
static double w_strip, w_ring, w_random;
static void fail(const char *what, double got)
{
    fprintf(stderr, "pass2_solved: %s (%.3f)\n", what, got);
    failed = 1;
}

int main(void)
{
    /* the coast strip */
    {
        enum { C = 400, R = 120, COAST = 70, EDGE = 100 };
        static float z[C * R];
        static unsigned char mask[C * R];
        int x, y;
        double worst = 0;
        for (y = 0; y < R; y++)
            for (x = 0; x < C; x++) {
                size_t k = (size_t) y * C + x;
                mask[k] = y < EDGE;
                z[k] = y < COAST - 2 ? 300.0f : y == COAST ? 0.0f : NO_ELEV;
            }
        if (isofill_diffuse(z, mask, NULL, C, R) != 0) { fail("the strip would not solve", 0); return 1; }
        for (y = COAST + 1; y < EDGE; y++)
            for (x = 0; x < C; x++)
                if (fabs(z[(size_t) y * C + x]) > worst) worst = fabs(z[(size_t) y * C + x]);
        w_strip = worst;
        if (worst > 0.02) fail("ground between a 0 coast and the drawn edge is not 0", worst);
    }

    /* the hilltop: a ring two cells wide at 2350, ground at 2000 answered outside it */
    {
        enum { C = 300, R = 300 };
        static float z[C * R];
        double lo = 1e9, hi = -1e9;
        int x, y;
        for (y = 0; y < R; y++)
            for (x = 0; x < C; x++) {
                double d = hypot(x - 150.0, y - 150.0);
                z[(size_t) y * C + x] = d >= 100 && d < 102 ? 2350.0f : d >= 110 ? 2000.0f : NO_ELEV;
            }
        if (isofill_diffuse(z, NULL, NULL, C, R) != 0) { fail("the hilltop would not solve", 0); return 1; }
        for (y = 0; y < R; y++)
            for (x = 0; x < C; x++)
                if (hypot(x - 150.0, y - 150.0) < 100) {
                    double v = z[(size_t) y * C + x];
                    if (v < lo) lo = v;
                    if (v > hi) hi = v;
                }
        w_ring = fmax(2350.0 - lo, hi - 2350.0);
        if (lo < 2349.98 || hi > 2350.02) fail("the floor inside a ring at 2350 is not flat at 2350", lo < 2349.98 ? lo : hi);
    }

    /* random ground against a solve in double precision */
    {
        enum { C = 96, R = 80 };
        static float z[C * R];
        static double ref[C * R];
        static unsigned char fixed[C * R];
        unsigned s = 12345;
        int x, y, it;
        double worst = 0;
        for (y = 0; y < R; y++)
            for (x = 0; x < C; x++) {
                size_t k = (size_t) y * C + x;
                s = s * 1103515245u + 12345u;
                fixed[k] = ((s >> 16) % 23) == 0;
                s = s * 1103515245u + 12345u;
                ref[k] = fixed[k] ? (double) ((s >> 16) % 1000) : 0.0;
                z[k] = fixed[k] ? (float) ref[k] : NO_ELEV;
            }
        for (it = 0; it < 200000; it++) {
            double big = 0;
            for (y = 0; y < R; y++)
                for (x = 0; x < C; x++) {
                    size_t k = (size_t) y * C + x;
                    double v;
                    if (fixed[k]) continue;
                    v = ((y > 0 ? ref[k - C] : 0) + (y < R - 1 ? ref[k + C] : 0)
                         + (x > 0 ? ref[k - 1] : 0) + (x < C - 1 ? ref[k + 1] : 0)) / 4;
                    if (fabs(v - ref[k]) > big) big = fabs(v - ref[k]);
                    ref[k] = v;
                }
            if (big < 1e-9) break;
        }
        if (isofill_diffuse(z, NULL, NULL, C, R) != 0) { fail("the random ground would not solve", 0); return 1; }
        for (y = 0; y < R * C; y++)
            if (fabs(z[y] - ref[y]) > worst) worst = fabs(z[y] - ref[y]);
        w_random = worst;
        if (worst > 0.02) fail("random ground is not the solve", worst);
    }

    /* one thread against all: the coast strip again, which takes cycles and settling both */
#ifdef _OPENMP
    {
        enum { C = 400, R = 120, COAST = 70, EDGE = 100 };
        static float z1[C * R], zn[C * R];
        static unsigned char mask[C * R];
        int x, y, all = omp_get_max_threads();
        for (y = 0; y < R; y++)
            for (x = 0; x < C; x++) {
                size_t k = (size_t) y * C + x;
                mask[k] = y < EDGE;
                z1[k] = zn[k] = y < COAST - 2 ? 300.0f + x : y == COAST ? 0.0f : NO_ELEV;
            }
        omp_set_num_threads(1);
        if (isofill_diffuse(z1, mask, NULL, C, R) != 0) { fail("the strip would not solve on one thread", 0); return 1; }
        omp_set_num_threads(all);
        if (isofill_diffuse(zn, mask, NULL, C, R) != 0) { fail("the strip would not solve", 0); return 1; }
        if (memcmp(z1, zn, sizeof z1) != 0) fail("one thread and all of them disagree", all);
    }
#endif
    if (failed) return 1;
    printf("libisofill %s: the second pass solved - the strip between a 0 coast and the drawn "
           "edge within %.4f m of 0, a ring's floor within %.4f m of the ring, random ground "
           "within %.4f m of the exact answer\n", isofill_version(), w_strip, w_ring, w_random);
    return 0;
}
