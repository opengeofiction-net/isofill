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
#define NO_ELEV (-32768.0f)

static int fail(const char *what)
{
    fprintf(stderr, "run_lib: %s\n", what);
    return 1;
}

int main(void)
{
    static float cons[COLS * ROWS], out[COLS * ROWS];
    static unsigned char mask[COLS * ROWS];
    isofill_params p = { 20, 1, 0.02, 1, 1 };
    long long filled;
    int x, y;

    if (strcmp(isofill_version(), ISOFILL_VERSION) != 0)
        return fail("isofill_version() is not the header's ISOFILL_VERSION");

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
     * second fills them; and the raster's top and bottom rows anchor toward
     * zero beyond the drawn area, as the second pass is meant to. Both are
     * isofill's behaviour, not the library's, and neither is what this asks.
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

    /* and bad arguments are refused, not acted on */
    if (isofill_run(NULL, 0, 0.0, NULL, NULL, COLS, ROWS, &p, out) != -1) return fail("NULL constraints accepted");
    p.radius = 0;
    if (isofill_run(cons, 0, 0.0, NULL, NULL, COLS, ROWS, &p, out) != -1) return fail("radius 0 accepted");

    printf("libisofill %s: fills between contours, leaves inputs alone, refuses nonsense\n",
           isofill_version());
    return 0;
}
