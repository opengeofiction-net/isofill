# The original

This is the `TileUtil` directory of OGF-terrain-tools, imported unchanged from
the `thilo-dem-process` tag, which is the last state of the Perl elevation
process before it was replaced in August 2026.

It is here for provenance. isofill reimplements the interpolation this code
performs, so keeping the source it was derived from - as the first commit, before
any work of ours - means the history itself records what was inherited. Nothing
in this directory is built or installed.

The parts that matter are in `tile_util.c`:

| | |
| --- | --- |
| `convert_tile_radius` | the fill: for every unset cell, gather the constraints within `radius` |
| `radius_value` | the rules: only constraints in line of sight, then the steepest visible pair, and no fill at all unless two differing elevations are among them |
| `reachable` | line of sight - the straight line to a constraint must cross no other |
| `make_line_points` | the ray that walk follows, two cells per step rather than one |
| `convert_tile_weighted_linear` | the second pass: each row and column interpolated between its known values, anchored at zero beyond both ends, the two blended towards whichever gradient is steeper |
| `interpolate_linear` | that row and column interpolation |

Two constants carry more weight than their size suggests. `radius` defaults to
20 cells, at the 3 arcsecond grid this ran on, so its reach is 1,852 m. And
`gradMax <= 0.1` - a tenth of a metre of rise per cell, about 1.1 m per km - is
what stops a cell being filled from a single elevation, which is why the sea
outside a coastline stays empty and comes out of the second pass as zero.

See *Admin:Elevation process* in the OGF admin wiki for why this is being
reimplemented: `gdal_fillnodata`, which replaced it, extrapolates from one
sample and terraces the surface into plateaus.
