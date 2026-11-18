# Spherical convergence control

This is the same finite-amplitude force-free magnetic transport solution as the
parent suite, using ordinary `spherical` geometry. The pure `yy_basis` and
`yy_position` helpers are used with patch 0 only to evaluate spherical components;
Yin-Yang initialization, communication and interpolation are disabled.

The spherical domain is r=1..2, theta=33.75..146.25 degrees and phi=0..292.5
degrees. The original Yin upper longitude boundary, phi=146.25 degrees, is now
inside this domain. The physical solution is unchanged (`placement=1`). Shifting
the longitude bounds by 146.25 degrees shifts them by exactly 26s cells, preserving
the underlying longitude cell locations. Each resolution has 16s x 20s x 52s
cells. All six faces impose the time-dependent exact solution; optional radial
divB correction applies only on the two radial faces, as in the Yin-Yang test.

Build only in an external case directory with `AMRVAC_DIR` pointing to an external
source copy. Copy this `mod_usr.t` into the case, using a 3D3V makefile. Then run:

```sh
python run_suite.py EXTERNAL_RUNS EXTERNAL_EXECUTABLE --case 1 --direction 1 --scale 1 --limiter woodward --radial-divbfix
```

Repeat with scales 2 and 3 in the same run root. A separate root with a debug
executable and `--debug --limiter woodward --radial-divbfix` checks both solutions
in both directions for three steps. Existing run directories are rejected.

```sh
python analyze.py SOURCE EXTERNAL_SPHERICAL_RUNS EXTERNAL_MATCHING_YINYANG_RUNS
```

The analysis writes `comparison.json` beside the spherical run root. It recomputes
both grids' Cartesian magnetic errors against the same cell-centre exact state
in the same moving sphere of radius 0.54, intersected with the radial shell.
L1 is normalized by the integrated magnetic perturbation, L2 by its squared
integral. Yin-Yang volumes retain the parent suite's overlap weights; ordinary
spherical volumes have no overlap factor. Global-domain errors are not compared
because the domains differ. The angular edges are outside the observation region.

Match numerical schemes, radial correction and CFL between the compared runs.
Equal CFL does not imply equal time-step histories. This is a local spatial
convergence control, not a full-sphere performance or pole-handling benchmark.
