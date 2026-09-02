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
| `--grad-min F` | the least gradient, in metres per cell, worth interpolating across. Default 0.02 |
| `--mask FILE` | fill only where this raster is non-zero |
| `--no-reach` | let the second pass carry values past the radius |
| `--no-pass2` | leave the cells the first pass declined unset |
| `--water FILE` | a raster the size of the constraints: where it is not zero the cell is water, held at zero by `--pass2 diffuse` |
| `--pass2 WHICH` | `linear`, the default, or `diffuse`. See below |
| `--pass2-tile N` | bound the second pass to a window. Default 0, the whole raster |
| `--max-mem MB` | above this, work band by band. Default 4096. The two passes decide separately; see below |
| `--threads N` | default: all but two |
| `--explain X Y` | say what one cell could see and what it did with it |

`--radius` and `--grad-min` carry more weight than their size suggests, and both
are about a **distance** rather than a cell count. The original's were set on a 3
arcsecond grid: radius 20 reaches 1,852 m, and 0.1 m per cell is 1.1 m per km. On
a finer grid they have to be scaled, or the fill silently searches a shorter
distance and refuses gentler slopes than it was tuned to.

`--radius` was scaled when this was written and `--grad-min` was not, which went
unnoticed until September 2026 because it costs coverage rather than correctness.
At 1 arcsecond a cell is 31 m, so 1.1 m per km is 0.034 m per cell, and the 0.1 it
kept was asking three times the slope the original did. Measured on a box inside
N32E067_Ellarca, dropping it to **0.02** takes the first pass from 32.71% of the
described area to 36.42% at radius 60, and from 55.86% to 74.82% at radius 120 -
where the floor alone was declining 18.95% of the ground. Below 0.02 nothing
changes at either radius, so that is the default rather than the 0.034 the
arithmetic gives.

The floor would be better expressed as a total drop, or scaled by the separation
of the pair it judges: as a rate it asks more of a distant pair than a near one,
which is backwards, since the distant pair is the one that had to survive more
ground to be seen at all.

## The second pass

`--pass2 linear` is the original's: interpolate along every row, along every
column, and blend the two by gradient. It is continuous in height and **not in
slope** - the surface is piecewise linear, so its derivative jumps at every
anchor. A hillshade is a derivative, so each joint draws an edge, and the edges
run along rows and columns because the segments do. On ground the first pass
cannot answer, that reads as rectangular blocks rather than as landform.

`--pass2 diffuse` solves Laplace's equation instead. Every cell the first pass
answered is held fixed and the rest relax to the mean of their four neighbours,
which is smooth in the derivative and has no preferred direction. A cell beyond
the edge of the raster counts as zero, reproducing the way the linear pass
anchors each row and column one step past its ends - that is what lets open
water come out at zero where the sea runs off the side of a zone. Water enclosed
by its own coastline needs no such help: bounded by cells at zero, Laplace gives
zero throughout.

It is solved coarse to fine, because relaxation moves information one cell per
sweep and the gaps are hundreds of cells wide.

Three sets of cells are held at zero rather than solved for, and each was needed:

* **Outside `--mask`.** The linear pass applies the mask afterwards, by zeroing
  what falls outside it. Doing the same here lets the solve run on past the edge
  of the drawn area and then chops the result off at it - on zone-ellarca, a
  1255 m cliff along the envelope.
* **`--water`.** The sea is a boundary, not something to solve for.
* **Out of reach, and connected to ground the contours never described.**
  Laplace has no notion of running out of information: its interior is a
  weighted average of its whole boundary, so a large region with no constraints
  in it takes the mean of everything around. On zone-ellarca that filled 56 km
  of undescribed ground with a 300 m dome, climbing away from the low contour
  edge towards the mountains on the far side.

That last one is why it is a flood rather than a rule about the verdict.
Zeroing *every* `OUT_OF_REACH` cell also punches a pit into each gap between
contours wider than the radius, and the fill then ramps from the surrounding
contour to zero across sixty cells - a gradient of a fifth, which hillshades as
a dark blob unrelated to the terrain, and is worse than the artefact it fixes.
What separates the two is where the region reaches, not what the first pass
said about it. See `mark_void`.

Decaying everything towards zero instead - solving (laplacian - lambda) z - was
tried and abandoned. It cannot tell no data from gentle ground between two
contours, and at any decay length short enough to suppress the first it flattens
the second. Measured on a box inside
N32E067_Ellarca, where the first pass answers a third of the ground: the same
13.5 s as the linear pass, the same surface to within a metre - median
difference 0 m, 95th percentile 3 m - and the longest run of axis-aligned steps
falls from 487 cells to 76.

**The two passes decide banding separately**, because the two banded paths are
not of equal standing. The first pass reads each band with a `--radius` margin,
so every interior cell still sees its whole search circle and a banded run is
**exact** - forcing it to band on zone-ellarca while the second pass holds the
raster gives output identical to the last bit, 100.0000% of cells and a maximum
difference of 0 m. The second pass banded is an approximation. So where there is
room for the solve but not for the first pass's summed area table - most of the
gap, at eight bytes a cell against the solve's seven and a third - it is taken,
and the exact answer kept.

For OGF that means `--max-mem 20480` puts every zone's second pass in memory:
the largest, zone-yuethon at 46801 square, wants 19.0 GB for the solve against
26.5 GB to hold its first pass whole.

Below that it is solved out of core, and **that path is an approximation**. It
enters the pyramid part way up - the fine rasters are read once and restricted
straight into a coarse grid that fits - and carries the answer back down in
bands, each solved with its margin rows held at the coarse answer so
neighbouring bands agree by construction. `--pass2-margin` sets that overlap.

Measured by forcing zone-ellarca out of core and comparing against itself solved
whole: 95% of cells identical, 95th percentile 0 m, and a tenth of a percent
differing by hundreds of metres. Of those, 28% are cells the banded path calls
void and the whole solve does not - the coarse grid deciding connectivity at 1
in 4 rather than cell by cell - and the rest is each band solving a different
problem from the whole raster.

Two things it is **not**, both measured before concluding them. Not the band
joins: the error away from one averages 5.30 m against 6.45 m within a hundred
rows of one, and the worst rows are nowhere near a join. And not the overlap:
widening the margin from 192 rows to 512 changes the output without moving any
of those figures. A first attempt that only relaxed each band from the coarse
answer, rather than solving it, was much worse - p99.9 of 1769 m against 586 -
because relaxation carries information one cell per sweep and the contours here
are hundreds of cells apart.

So give `--max-mem` the room where you can, and treat the banded second pass as
a fallback that gets you a map rather than the map.

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
