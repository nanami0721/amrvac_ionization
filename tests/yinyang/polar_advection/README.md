# Matched-resolution polar advection experiment

This test compares standard full-sphere latitude–longitude and fixed Yin–Yang
MHD on the same branch, using the same executable, TVDLF/minmod/three-step
integration, CFL 0.3, GLM and ideal-gas EOS. It does not force either solver to
crash. Four MPI ranks, sequential jobs, no AMR, no polar cell merging.

Exact solution:

rho = 1 + 0.5 exp(-|x - R(-0.3+t,0,1.5)|^2 / 0.18^2),
p = 0.1, v = R(1,0,0), B = R(0.2,0.1,0.3).

The rotation R is identity for the polar trajectory (crossing the north axis at
t=0.3), a 90-degree rotation about y for the equatorial control, and an additional
orientation at longitude 146.25 degrees for the overset-interface control. All
vectors and the density profile are rotated together. Radial ghosts use the exact
solution at the current BC time. True spherical poles use AMRVAC pi-periodic
pole connectivity; longitude is periodic. There is no polar excision or imposed
exact angular boundary in the spherical baseline.

Matched grids (r in [1,2]):

| scale | spherical | Yin–Yang per patch | total cell ratio YY / spherical |
|---|---|---|---|
| 1 | 16x32x64 | 16x20x52 | 1.015625 |
| 2 | 32x64x128 | 32x40x104 | 1.015625 |

Yin–Yang angular extents are theta=[33.75,146.25] degrees and
phi=[-146.25,146.25] degrees. Both angular increments are **identical** to the
spherical baseline: 5.625 degrees (coarse), 2.8125 degrees (fine). Block size is
8x4x4 for both grids. Radial increments likewise match. This matches equatorial
physical resolution; spherical polar cells are intrinsically much smaller.

Snapshots at t=0,0.3,0.6. Report initial dt, total steps, evolution wall time,
full-volume and pulse-region relative L1/L2 error, peak retention, centroid,
velocity/B errors and positivity. Pulse region is a ball of radius 3a=0.54 about
the exact pulse centre; the same physical region is used for all orientations.
All normalized density errors divide by the exact **density excess**, not the
large uniform background. Full-volume relative errors can therefore exceed 100%
when a small density error extends over the entire shell; retain this metric
alongside the pulse-region metric. Peak fraction refers to the analytic amplitude
0.5, and the maximum sampled exact peak is also reported to expose sampling bias.
Centroids use positive density excess within the pulse region.

Global integrals use exact spherical cell volumes; Yin–Yang overlap partition
weights (1/2 inside common coverage, 1 elsewhere) use 16x16 angular subcell
quadrature. Report the residual volume bias. These are diagnostics, not a
conservative exchange scheme. Comparing equal angular resolution does not itself
establish a speedup at equal error.

Build in an external directory with an external source copy. `run_experiment.py`
creates case directories and runs sequentially (`--fine` for scale 2, `--seam`
for orientation 2). `analyze.py SOURCE RUN_ROOT` writes results.json one level
above RUN_ROOT. `plot_results.py SOURCE RUN_ROOT FIGURE_DIR` renders plots after
all timing jobs have finished. Logs and data must remain outside the repository.

Magnetic divergence is additionally estimated by centred differences in local
spherical components, omitting one cell at each patch/grid boundary. It is not
the solver's discrete GLM divergence and does not include the omitted pole or
interface cells; use it together with global Cartesian B error, not as proof
that either solver satisfies a divergence constraint everywhere.

The angular patches deliberately include a substantial overlap so every ghost
stencil has interior donors. Raw angular area counted across both patches exceeds
4π; the total **cell** ratio above includes all duplicate overlap cells. This is
not a claim of equal physical cell sizes near spherical poles. The timings are
local single-run observations; step counts are the more reproducible CFL metric.

After the two-level pilot, an adaptive follow-up (`--accuracy`, named s3) uses
Yin–Yang **40x48x124 per patch**, polar trajectory only. Its aim is to match or
improve the measured fine spherical polar pulse-region L1 error, then compare
actual cost. This is not a third uniform doubling level and was selected after
seeing the first two levels. Do not include it in a uniform-refinement scaling
fit. One follow-up resolution can establish a shared error bound, not the optimal
cost at that bound or equality of every error metric.

Example external execution after building the user module:

```
python run_experiment.py /private/tmp/amrvac-tests/yy-advection-UNIQUE/runs \
  --executable /private/tmp/amrvac-tests/yy-advection-UNIQUE/build-case/amrvac
python run_experiment.py /private/tmp/amrvac-tests/yy-advection-UNIQUE/runs \
  --executable /private/tmp/amrvac-tests/yy-advection-UNIQUE/build-case/amrvac --fine
python run_experiment.py /private/tmp/amrvac-tests/yy-advection-UNIQUE/runs \
  --executable /private/tmp/amrvac-tests/yy-advection-UNIQUE/build-case/amrvac --seam
python run_experiment.py /private/tmp/amrvac-tests/yy-advection-UNIQUE/runs \
  --executable /private/tmp/amrvac-tests/yy-advection-UNIQUE/build-case/amrvac --accuracy
```

Use the same external build-source with AMRVAC_DIR pointing there. Do not run
these independent cases simultaneously when measuring time. Case directories
must be new, to avoid mixing old snapshots and new runs.
