/*
 * The library, exercised: fill a small raster between two contours and check
 * what a caller is entitled to. Modest, as the spec said isofill's C tests
 * would be; the golden reference in the danu repository is the one that
 * pins the surface itself.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "isofill.h"

#define COLS 32
#define ROWS 16
#define NO_ELEV ((float) ISOFILL_NO_ELEV)

static int fail(const char *what)
{
    fprintf(stderr, "run_lib: %s\n", what);
    return 1;
}

int main(void)
{
    static float cons[COLS * ROWS], out[COLS * ROWS];
    static unsigned char mask[COLS * ROWS];
    isofill_params p;
    long long filled;
    int x, y;

    if (strcmp(isofill_version(), ISOFILL_VERSION) != 0)
        return fail("isofill_version() is not the header's ISOFILL_VERSION");

    /* start from the binary's defaults and set only what this test means to */
    isofill_params_default(&p);
    if (p.radius != 20 || p.barrier != 1 || p.grad_min != 0.02 || p.pass2 != 1)
        return fail("isofill_params_default() is not the command line's defaults");
    p.threads = 1;
    if (isofill_whole_mb(1000, 1000) < 12.0 || isofill_whole_mb(1000, 1000) > 13.0)
        return fail("isofill_whole_mb(1000, 1000) is not about 12.4 MB");

    /* two contours: 100 m down column 4, 200 m down column 27; everything
     * else unset, and the whole raster inside the mask */
    for (y = 0; y < ROWS; y++)
        for (x = 0; x < COLS; x++) {
            cons[y * COLS + x] = x == 4 ? 100.0f : x == 27 ? 200.0f : NO_ELEV;
            mask[y * COLS + x] = 1;
        }

    filled = isofill_run(cons, 0, 0.0, mask, NULL, COLS, ROWS, &p, out);
    if (filled < 0) return fail("isofill_run returned an error");
    if (filled == 0) return fail("the first pass set nothing between two contours in plain sight");

    /*
     * The interior only. The two columns beside each contour are inside the
     * --barrier dilation, so the first pass sees one level there and the
     * second fills them; and in the raster's top and bottom rows the unset
     * cells anchor toward zero beyond the drawn area, as the second pass is
     * meant to - the constraints themselves are fixed cells and do not move,
     * which the loop above checks on every row. Both are isofill's behaviour,
     * not the library's, and neither is what this asks.
     */
    for (y = 0; y < ROWS; y++)
        if (out[y * COLS + 4] != 100.0f || out[y * COLS + 27] != 200.0f)
            return fail("a constraint did not come back as itself");
    for (y = 3; y < ROWS - 3; y++) {
        float prev = 100.0f;
        for (x = 7; x <= 24; x++) {
            float v = out[y * COLS + x];
            if (v < 100.0f || v > 200.0f) return fail("a cell between the contours is outside them");
            if (v < prev) return fail("the surface between the contours is not monotonic");
            prev = v;
        }
        /* and the barrier-adjacent cells the first pass declined were filled by the second */
        for (x = 5; x <= 6; x++)
            if (out[y * COLS + x] <= -32000.0f) return fail("a declined cell was left as a sentinel");
    }
    /* the caller's arrays were not written */
    if (cons[7 * COLS + 10] != NO_ELEV) return fail("isofill_run wrote to the constraints");

    /*
     * isofill_run_ex hands back the first pass from the same run. What it
     * gives has to be what a second run with pass2 off would write, cell for
     * cell, or a caller reading it is reading something else - and the finished
     * surface has to be the one isofill_run gives, so keeping the first pass
     * costs nothing but the copy.
     */
    {
        static float both[COLS * ROWS], p1[COLS * ROWS], only1[COLS * ROWS];
        long long f_ex, f_1;
        size_t i;
        isofill_params_default(&p);
        p.threads = 1;
        f_ex = isofill_run_ex(cons, 0, 0.0, mask, NULL, COLS, ROWS, &p, both, p1);
        if (f_ex != filled) return fail("isofill_run_ex filled a different number of cells");
        for (i = 0; i < (size_t) COLS * ROWS; i++)
            if (both[i] != out[i]) return fail("isofill_run_ex's surface is not isofill_run's");
        p.pass2 = 0;
        f_1 = isofill_run(cons, 0, 0.0, mask, NULL, COLS, ROWS, &p, only1);
        if (f_1 != filled) return fail("--no-pass2 filled a different number of cells");
        for (i = 0; i < (size_t) COLS * ROWS; i++)
            if (p1[i] != only1[i]) return fail("the kept first pass is not what --no-pass2 writes");
        /* and it really is the first pass, not the finished surface: the cells
         * beside the barrier that the second pass filled are still sentinels */
        if (p1[7 * COLS + 5] > -32000.0f) return fail("the kept first pass has no declined cells in it");
        /* NULL is allowed and is what isofill_run passes */
        isofill_params_default(&p);
        p.threads = 1;
        if (isofill_run_ex(cons, 0, 0.0, mask, NULL, COLS, ROWS, &p, both, NULL) != filled)
            return fail("isofill_run_ex with no pass1_out is not isofill_run");
    }

    /* and bad arguments are refused, not acted on */
    if (isofill_run(NULL, 0, 0.0, NULL, NULL, COLS, ROWS, &p, out) != -1) return fail("NULL constraints accepted");
    p.radius = 0;
    if (isofill_run(cons, 0, 0.0, NULL, NULL, COLS, ROWS, &p, out) != -1) return fail("radius 0 accepted");

    /*
     * isofill_diffuse: the second pass on its own. Given the surface pass 1
     * left, it has to produce what isofill_run produces, cell for cell, or a
     * caller composing the two passes by hand is not composing the same fill.
     */
    {
        static float whole[COLS * ROWS], p1[COLS * ROWS], piece[COLS * ROWS];
        size_t i;
        isofill_params_default(&p);
        p.threads = 1;
        if (isofill_run_ex(cons, 0, 0.0, mask, NULL, COLS, ROWS, &p, whole, p1) < 0)
            return fail("isofill_run_ex returned an error");
        memcpy(piece, p1, sizeof piece);
        if (isofill_diffuse(piece, mask, NULL, COLS, ROWS) != 0)
            return fail("isofill_diffuse returned an error");
        for (i = 0; i < (size_t) COLS * ROWS; i++)
            if (piece[i] != whole[i])
                return fail("isofill_diffuse is not the pass 2 isofill_run does");
        if (isofill_diffuse(NULL, NULL, NULL, COLS, ROWS) != -1)
            return fail("isofill_diffuse accepted a null surface");
    }

    /*
     * And what it is for: a box cut out of a larger surface, with its rim
     * written from the answer the whole raster gave. Both contours here are on
     * the left, so the right-hand two thirds is out of reach of either - the
     * first pass declines it, and a whole-raster second pass floods it from the
     * raster edge and holds it at zero. A rim that says the ground goes on is
     * a different boundary, and the solve has to run up to it.
     */
    {
#define WCOLS 64
        static float far_cons[WCOLS * ROWS], loose[WCOLS * ROWS];
        static float far_p1[WCOLS * ROWS], bounded[WCOLS * ROWS];
        static unsigned char wide_mask[WCOLS * ROWS];
        int yy, xx, differs = 0;
        size_t i, mid = (size_t) (ROWS / 2) * WCOLS;
        for (yy = 0; yy < ROWS; yy++)
            for (xx = 0; xx < WCOLS; xx++) {
                size_t k = (size_t) yy * WCOLS + xx;
                far_cons[k] = xx == 4 ? 100.0f : xx == 10 ? 200.0f : NO_ELEV;
                /* a strip nobody drew, so the mask has something to do. It sits
                 * between the two contours, where the ground either side is in
                 * reach of both: out in the far third it would instead seed the
                 * void flood and hold everything beyond it at zero, which is
                 * the algorithm working and not this scenario */
                wide_mask[k] = (xx == 6 || xx == 7) ? 0 : 1;
            }
        isofill_params_default(&p);
        p.threads = 1;
        if (isofill_run_ex(far_cons, 0, 0.0, wide_mask, NULL, WCOLS, ROWS, &p, loose, far_p1) < 0)
            return fail("the wide raster would not fill");
        if (loose[mid + WCOLS - 2] != 0.0f)
            return fail("the far ground is not zeroed without a rim, so this tests nothing");

        memcpy(bounded, far_p1, sizeof bounded);
        for (yy = 0; yy < ROWS; yy++)
            for (xx = 0; xx < WCOLS; xx++)
                if (xx == 0 || xx == WCOLS - 1 || yy == 0 || yy == ROWS - 1)
                    bounded[(size_t) yy * WCOLS + xx] = 300.0f;
        if (isofill_diffuse(bounded, wide_mask, NULL, WCOLS, ROWS) != 0)
            return fail("the bounded diffuse returned an error");

        for (yy = 0; yy < ROWS; yy++)
            if (bounded[(size_t) yy * WCOLS + WCOLS - 1] != 300.0f)
                return fail("a rim cell moved - the rim was flooded or solved over");
        if (bounded[mid + WCOLS - 2] <= 0.0f)
            return fail("the ground beside the rim was not solved towards it");
        for (i = 0; i < (size_t) WCOLS * ROWS; i++)
            if (bounded[i] != loose[i]) differs++;
        if (differs < ROWS)
            return fail("the rim changed nothing, so it is not bounding the solve");
        /* the undrawn strip is zero, which is what passing the mask through is for */
        for (yy = 1; yy < ROWS - 1; yy++)
            for (xx = 6; xx <= 7; xx++)
                if (bounded[(size_t) yy * WCOLS + xx] != 0.0f)
                    return fail("a cell outside the mask was filled");
    }

    printf("libisofill %s: fills between contours, keeps the first pass, "
           "diffuses a box against a written rim, leaves inputs alone, "
           "refuses nonsense\n",
           isofill_version());
    return 0;
}
