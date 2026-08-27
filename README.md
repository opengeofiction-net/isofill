# isofill

Fill a raster from rasterised contour lines: a continuous surface from lines that
only say what the elevation is where they are drawn.

For each unset cell isofill gathers the contours within a search radius which are
in **line of sight** - the straight line to them crossing no other contour - takes
the pair with the steepest gradient across the cell, and interpolates between them.
Where no two differing elevations are in sight it declines to fill the cell at all
rather than extrapolating from one. A second pass then interpolates every row and
column between its known values, anchored at zero beyond both ends, so the cells
the first pass declined settle at the level of whatever surrounds them.

    isofill contours.tif dem.tif
    isofill --radius 60 --barrier 2 contours.tif dem.tif

The input is a raster of burned contour lines - each cell either carries an
elevation or is nodata. `gdal_rasterize -a ele -at` from a vector layer produces
one. The output is `Int16`, on the same grid.

## Why not gdal_fillnodata

`gdal_fillnodata` takes the nearest valid cell in each of eight fixed compass
directions and weights by distance, and it will interpolate from a **single**
sample. That has two consequences on contour data.

It invents ground where there is no evidence for any. A headland casts a plume of
elevation downwind of itself, out into open water, because one coastline cell is
all the fill needs. Bounding it by distance limits how far the plume reaches
without stopping it being drawn.

And it terraces. Whole regions take their value from one nearest sample, so the
surface comes out as plateaus with straight boundaries where the chosen sample
switches from one contour to another. Compared against the process isofill
reimplements, on a common grid, the two agree on the terrain - median slope 1.08%
either way, p90 5.40%, p99 14.04% - and differ on its texture: 42.1% of
neighbouring cells share a value against 33.4%.

That difference barely moves a slope histogram and dominates a hillshade, which is
a derivative and renders every plateau edge as a line. Judge an interpolator for
this job by shading its output, not by its error statistics.

`gdal_grid` is the other obvious tool. It triangulates, taking no account of
contours as barriers, and on the smallest zone tested it took eight hours where
`gdal_fillnodata` took seconds.

## Building

Needs GDAL and a C compiler. OpenMP is used if available.

    make
    sudo make install          # PREFIX=/usr/local

Or build a Debian package:

    dpkg-buildpackage -us -uc -b
    sudo dpkg -i ../isofill_*.deb

Check the built package depends on `libgomp1`. If it does not, `-fopenmp` was
dropped somewhere and the binary is quietly single threaded.

## Options

| | |
| --- | --- |
| `--radius N` | how far a cell may look, in cells. Default 20 |
| `--barrier N` | widen contours for the sight test only, not for their values. Default 1 |
| `--grad-min F` | the least gradient, in metres per cell, worth interpolating across. Default 0.1 |
| `--mask FILE` | fill only where this raster is non-zero |
| `--no-reach` | let the second pass carry values past the radius |
| `--no-pass2` | leave the cells the first pass declined unset |
| `--pass2-tile N` | bound the second pass to a window. Default 0, the whole raster |
| `--max-mem MB` | above this, work band by band. Default 4096 |
| `--threads N` | default: all but two |
| `--explain X Y` | say what one cell could see and what it did with it |

`--radius` and `--grad-min` carry more weight than their size suggests, and both
are about a **distance** rather than a cell count. The defaults are the original's,
set on a 3 arcsecond grid: radius 20 reaches 1,852 m, and 0.1 m per cell is 1.1 m
per km. On a finer grid scale them, or the fill silently searches a shorter
distance and refuses gentler slopes than it was tuned to.

`--barrier` exists for the same reason from the other direction. A one cell contour
is a 93 m wall at 3 arcseconds and a 31 m wall at 1, so holding the search distance
while the grid gets finer leaves rays threading gaps that could not exist on the
coarser grid. Widening the rasterised line instead would write its elevation into
every cell it gained and drag the surface towards each line: what a cell is worth
and what it hides are separate questions.

## Memory, and large rasters

Above `--max-mem` the first pass runs a band of rows at a time, each read with a
margin of `radius`, and the second pass streams. A cell only ever looks `radius`
away, so a banded result is identical to a whole-raster one rather than an
approximation of it - checked at 0 differing cells of 51,854,401.

The bound matters more than it sounds. A summed area table over the constraint
mask, which is what lets a cell with nothing within reach be skipped in constant
time, costs eight bytes a cell: 42 GB on a 86401 x 39601 raster. That same table
is why a sparse raster costs little - a zone 92% of which nobody has drawn fills in
much the time a small dense crop does.

Threading is close to linear: on a 7201 x 7201 raster, 24m37 on one core, 12m20 on
two, 6m08 on four, 4m09 on six.

## Provenance

isofill reimplements the interpolation in the terrain code of
[OGF-terrain-tools](https://github.com/opengeofiction-net/OGF-terrain-tools),
written by Thilo Stapff. The first commit on this repository is that C, imported
unchanged, so what was inherited and what was written here can be told apart from
the history alone. It is kept in `original/` for reference and is not built.

What is ours rather than the original's is in the header of `src/isofill.c`: only
the nearest constraint at each distinct elevation is kept, which is exact and turns
an O(k²) search into one over the handful of levels in view; the summed area table;
banding; threading; and the reach and barrier rules above.

Licensed Artistic 1.0 or GPL 1 or later, inherited from the code it derives from.
See [LICENSE](LICENSE).

## Known

Near a shore, isofill leaves more stray elevation in open water than the process it
reimplements, in streaks a few metres high. Rasterising all touched and
`--barrier 2` together remove a little over a quarter of it. The cause of the rest
is **not known**: three explanations were tested and each refuted by measurement -
the fill walking around the end of an open coastline way, which closed coastlines
still show; an unbounded second pass window, which changes it by less than a fifth
of a percent; and ray density, where the original's own coarser geometry is worse
still. Finding it means differencing the first pass cell by cell against the
original on identical input, rather than more reasoning from its source.
