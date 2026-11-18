# Tilted solar dipole demonstration

A potential dipole tilted 20 degrees from global z, in a shell 1 <= r <= 3.
The normalized field is B = [3(m.n)n - m]/(2 r^3), with |m|=1.
The displayed physical scale is R_sun = 1 and polar field = 10 G; the solver
itself uses normalized values. Uniform rho=1, p=0.05 and zero velocity form
an analytic static equilibrium without gravity, heating or wind. Discrete
operators/interpolation do not preserve this equilibrium exactly.

Both radial boundaries prescribe the analytic state. Angular ghosts use the
Yin–Yang exchange. No B0 splitting is used; total B is evolved by MHD+GLM.
Each patch has 32x32x96 cells, with block size 8x8x8. Run on an even number of
MPI ranks, with compilation and all outputs outside the repository.

`plot_dipole.py SOURCE RUN_DIR OUTPUT_DIR` uses the native snapshots at step 0
and step 5, and requires NumPy, SciPy and Matplotlib. Field lines are integrated
by RK4 through interpolated **snapshot fields**, selecting the patch furthest
from its angular edge. At each patch centre the vectors are rotated to global
Cartesian components before interpolation. The surface colour is the analytic
photospheric radial field. Lines start at r=1.04 and stop near the inner shell or
at r=3; reaching the outer plotting boundary does NOT imply truly open flux.
The meridional map uses the numerical snapshot; the final half radial cell is
linearly extrapolated to the displayed outer limit. This is a vacuum dipole,
not a source-surface or solar-wind model.

2026-11-16: four-rank debug run completed 5 steps with empty stderr. Initial
field error <= 1.34e-15 in normalized units. Relative L2 field change after
five steps was < 5e-4; this is not a long-term equilibrium validation.
64 snapshot field lines completed tracing; maximum drift of the analytic
invariant r/sin(theta_m)^2 along a traced line was 8.92e-4.

Artifacts: /private/tmp/amrvac-tests/yinyang-dipole-20261116/.
Source/build basis: codex/yinyang-grid-prep, base 6f5e6f10 plus uncommitted v1.
