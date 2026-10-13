/*
 * The second pass out of core against the same raster solved whole: the same,
 * to the bit. The binary is run twice on one raster - with --max-mem room for
 * it, and with a megabyte, so that every array the second pass makes goes to
 * a file - and the two outputs compared cell for cell.
 *
 * The raster is drawn to exercise what the banded approximation got wrong:
 * contours as random walks at many levels, far enough apart that the first
 * pass leaves wide ground unreached; a drawn area that does not cover the
 * raster, so there is void reaching its edge; rings that enclose unreached
 * ground the void must not take; and a coast with sea held at zero.
 *
 * Run from the repository root, with ./isofill built. It works in a directory
 * under tests/, or under $TMPDIR if that is set - not on Windows, where the
 * shell's TMPDIR may be one only MSYS2 can read.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>
#include <gdal.h>
#include <cpl_conv.h>
#ifdef _WIN32
#include <direct.h>
#include <process.h>
#define ISOFILL ".\\isofill.exe"   /* system() is cmd.exe's, which cannot run ./ */
#else
#define ISOFILL "./isofill"
#endif

#define COLS 900
#define ROWS 700
#define NODATA (-32768.0f)

static float cons[ROWS][COLS];
static unsigned char mask[ROWS][COLS], water[ROWS][COLS];

static int fail(const char *what)
{
    fprintf(stderr, "pass2_ooc: %s\n", what);
    return 1;
}

static unsigned long long state = 88172645463325252ULL;
static double rnd(void)
{
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return (double) (state >> 11) / 9007199254740992.0;
}

static void walk(float level, int steps)
{
    double x = 40 + rnd() * (COLS - 80), y = 40 + rnd() * (ROWS - 80), a = rnd() * 6.283;
    int i;
    for (i = 0; i < steps; i++) {
        a += (rnd() - 0.5) * 0.5;
        x += cos(a); y += sin(a);
        if (x < 1 || y < 1 || x >= COLS - 1 || y >= ROWS - 1) break;
        cons[(int) y][(int) x] = level;
    }
}

static void ring(int cx, int cy, int r, float level)
{
    int t;
    for (t = 0; t < 4 * 360; t++) {
        double a = t * 3.14159265358979 / 720;
        int x = cx + (int) lround(r * cos(a)), y = cy + (int) lround(r * sin(a));
        cons[y][x] = level;
    }
}

static int write_tif(const char *path, GDALDataType type, const void *data)
{
    GDALDatasetH d = GDALCreate(GDALGetDriverByName("GTiff"), path, COLS, ROWS, 1, type, NULL);
    double gt[6] = { 0, 1, 0, 0, 0, -1 };
    if (!d) return 1;
    GDALSetGeoTransform(d, gt);
    if (type == GDT_Float32) GDALSetRasterNoDataValue(GDALGetRasterBand(d, 1), NODATA);
    if (GDALRasterIO(GDALGetRasterBand(d, 1), GF_Write, 0, 0, COLS, ROWS, (void *) data,
                     COLS, ROWS, type, 0, 0) != CE_None) { GDALClose(d); return 1; }
    GDALClose(d);
    return 0;
}

static float *read_tif(const char *path)
{
    GDALDatasetH d = GDALOpen(path, GA_ReadOnly);
    float *v = malloc((size_t) COLS * ROWS * sizeof *v);
    if (!d || !v) { if (d) GDALClose(d); free(v); return NULL; }
    if (GDALRasterIO(GDALGetRasterBand(d, 1), GF_Read, 0, 0, COLS, ROWS, v, COLS, ROWS,
                     GDT_Float32, 0, 0) != CE_None) { free(v); v = NULL; }
    GDALClose(d);
    return v;
}

static char *make_dir(char *dir)
{
#ifdef _WIN32
    int i;
    for (i = 0; i < 100; i++) {
        snprintf(dir, 1024, "tests/isofill-pass2-ooc-%d-%d", _getpid(), i);
        if (_mkdir(dir) == 0) return dir;
    }
    return NULL;
#else
    const char *here = getenv("TMPDIR");
    if (here && *here) snprintf(dir, 1024, "%s/isofill-pass2-ooc-XXXXXX", here);
    else snprintf(dir, 1024, "tests/isofill-pass2-ooc-XXXXXX");
    return mkdtemp(dir);
#endif
}

/* everything in it, then it: the runs make nothing but files */
static int clear_dir(const char *dir)
{
    DIR *d = opendir(dir);
    struct dirent *e;
    char path[2048];
    int bad = 0;
    if (!d) return 1;
    while ((e = readdir(d))) {
        if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, "..")) continue;
        snprintf(path, sizeof path, "%s/%s", dir, e->d_name);
        if (remove(path)) bad = 1;
    }
    closedir(d);
    return bad || rmdir(dir);
}

static int leftovers(const char *dir)
{
    DIR *d = opendir(dir);
    struct dirent *e;
    int n = 0;
    if (!d) return 0;
    while ((e = readdir(d)))
        if (!strncmp(e->d_name, "isofill-pass2-", 14)) n++;
    closedir(d);
    return n;
}

int main(void)
{
    char dir[1024], cmd[4096], path[1024];
    float *whole, *ooc;
    int x, y, i, n_diff = 0, ooc_used_files;
    FILE *f;

    if (!make_dir(dir)) return fail("cannot make a working directory");
    GDALAllRegister();

    for (y = 0; y < ROWS; y++)
        for (x = 0; x < COLS; x++) {
            cons[y][x] = NODATA;
            /* the drawn area a rounded box short of every edge but the west,
             * so the void reaches in from three sides */
            mask[y][x] = (x < 820 && y > 60 && y < 640) ? 1 : 0;
            /* sea in the south west, behind a coast drawn at zero */
            water[y][x] = (y > 560 && x < 300) ? 1 : 0;
        }
    for (i = 0; i < 40; i++) walk((float) (25 * (1 + i % 16)), 700);
    ring(450, 330, 90, 500.0f);          /* encloses ground the radius cannot reach */
    ring(650, 450, 60, 200.0f);
    for (x = 0; x < 300; x++) cons[558][x] = 0.0f;

    snprintf(path, sizeof path, "%s/cons.tif", dir);
    if (write_tif(path, GDT_Float32, cons)) return fail("cannot write cons.tif");
    snprintf(path, sizeof path, "%s/mask.tif", dir);
    if (write_tif(path, GDT_Byte, mask)) return fail("cannot write mask.tif");
    snprintf(path, sizeof path, "%s/water.tif", dir);
    if (write_tif(path, GDT_Byte, water)) return fail("cannot write water.tif");

    snprintf(cmd, sizeof cmd, ISOFILL " --radius 20 --threads 2 --mask %s/mask.tif --water %s/water.tif "
             "--max-mem 4096 %s/cons.tif %s/whole.tif 2>%s/whole.log", dir, dir, dir, dir, dir);
    if (system(cmd) != 0) return fail("the whole run failed");
    snprintf(cmd, sizeof cmd, ISOFILL " --radius 20 --threads 2 --mask %s/mask.tif --water %s/water.tif "
             "--max-mem 1 %s/cons.tif %s/ooc.tif 2>%s/ooc.log", dir, dir, dir, dir, dir);
    if (system(cmd) != 0) return fail("the out-of-core run failed");

    /* that it really went out of core: the log says how much of the second
     * pass's arrays were in files, and the surface alone is four bytes a cell */
    {
        double filed = -1;
        char line[512];
        snprintf(path, sizeof path, "%s/ooc.log", dir);
        f = fopen(path, "r");
        if (!f) return fail("cannot read the out-of-core run's log");
        while (fgets(line, sizeof line, f))
            if (strstr(line, "of its arrays were in files")) sscanf(line, " pass 2: %lf", &filed);
        fclose(f);
        ooc_used_files = filed >= 4.0 * COLS * ROWS;
        if (!ooc_used_files) return fail("--max-mem 1 did not put the second pass in files");
    }
    if (leftovers(dir)) return fail("the second pass left its files behind");

    snprintf(path, sizeof path, "%s/whole.tif", dir);
    whole = read_tif(path);
    snprintf(path, sizeof path, "%s/ooc.tif", dir);
    ooc = read_tif(path);
    if (!whole || !ooc) return fail("cannot read the outputs");
    for (i = 0; i < COLS * ROWS; i++)
        if (memcmp(&whole[i], &ooc[i], sizeof *whole)) n_diff++;
    {
        int unset = 0, zero = 0;
        for (i = 0; i < COLS * ROWS; i++) {
            if (whole[i] == NODATA) unset++;
            if (whole[i] == 0.0f) zero++;
        }
        /* the fixture is worth having: the second pass filled ground, and
         * held the void and the sea at zero */
        if (unset) return fail("a cell the second pass left unset");
        if (zero < COLS * ROWS / 10) return fail("the fixture has too little void or sea to test");
    }
    if (n_diff) {
        fprintf(stderr, "pass2_ooc: %d of %d cells differ between the whole and the out-of-core run\n",
                n_diff, COLS * ROWS);
        return 1;
    }
    printf("pass 2 out of core: the same as solved whole, all %d cells to the bit\n", COLS * ROWS);
    free(whole);
    free(ooc);
    if (clear_dir(dir)) return fail("cannot clear the working directory");
    return 0;
}
