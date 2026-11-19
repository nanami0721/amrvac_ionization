# Analytic Yin–Yang CT transfer attribution

Test-only diagnostics on the actual fixed YY spherical mesh, with no evolution
(`it_max=0`) and no modification of production CT routines. Run outside the
repository via `run.py`; the runner copies the working source (including local
changes), records commit/diff/hashes, reuses only compatible cached build
products, and guards serial builds at 4096 MiB and MPI runs at 1024 MiB.
System thresholds remain 65% start / 68% stop; sampling is not a hard quota.
Use the project Python environment and a previous external source cache:

```sh
python3 tests/yinyang/ct_accuracy/run.py /private/tmp/amrvac-tests/UNIQUE \
  --cache-source /private/tmp/amrvac-tests/PREVIOUS/source
python3 tests/yinyang/interface_transport/convergence/memory_guard.py \
  --out /private/tmp/amrvac-tests/UNIQUE/analysis-memory --limit-mib 1024 -- \
  python3 tests/yinyang/ct_accuracy/analyze.py /private/tmp/amrvac-tests/UNIQUE
```

The full matrix has debug n=8,12,14 (np=2), n=14 with quadrature16,
n=8 with np=4, and optimized n=8,12,14. Each patch has n×2n×6n cells
in 12 blocks, `max_blocks=128`. Refinement forecasts use the nearest measured
coarser run and reject estimates over 950 MiB. `--pilot` runs only debug n=8.

Fields: constant Cartesian B; affine divergence-free B=.2(y,z,x); the weak
R=.65 translated magnetic loops at forward t=.404 and reverse t=.522 from
`ct_transport`. Fields are freshly prescribed, not evolved. E=-v×B. Face B
is constructed as the discrete Stokes curl of quadrature-integrated vector
potential, divided by native face area, so the reference is discretely
solenoidal. This differs from a point-B reference or midpoint-A initialization.

Edge modes:

- 0: native midpoint E×length donors, legacy two-Gauss-node receiver integrals.
- 1: integrated donor E, testing line-average-as-midpoint interpretation.
- 2: mode 0 plus an **analytic oracle** correction on receiver edges,
  `midpoint_exact - integral_exact`, to isolate mixed quadrature conventions.
  This is not an implementable reconstruction upgrade or a transport result.

- 3: native midpoint donors and receiver midpoint quadrature, the current runtime convention.
  Constant Cartesian E must recover the native midpoint curl to roundoff.

Ghost stages: 0 reference+closure; 1 interpolation only; 2 interpolation+closure.
Active faces must remain bitwise unchanged. Reference+closure must preserve
B within 1e-11 amplitude units; constrained native divB within 1e-10.

`EDGE_ACCURACY` and `GHOST_ACCURACY` each contain field, mode/stage, 14 sums
and 6 maxima. There are 28 summary rows per case. Detailed E/B records carry
field, mode/stage, patch, global indices, component, reference and actual values
(E additionally includes midpoint). `analyze.py` compares canonical geometric
records for MPI, quadrature and optimized/debug controls. Face samples at
shared block boundaries can be duplicated in RMS sums; comparisons remove
those duplicates. RMS measures are unweighted sample norms, not physical
volume integrals or conservation errors. Curl errors are per unit time and
normalized by B amplitude; they are not full-evolution convergence orders.
All-angular-interface RMS includes field-free samples for compact loops.

Interpret n=8→14 orders as finite-range slopes. The n=14 loop results can be
nonmonotonic as donor stencils and support intersections change. Quadrature
8→16 bounds reference sensitivity, not exactness at compact-support crossings.
Constant-field relative curl errors are undefined; use absolute normalized RMS.

All run logs/data stay outside the repository. Failed pilots and memory-
forecast stops should be retained. The test helper uses a pure output-argument
quadrature subroutine: the local debug compiler reported a recursive
function entry with the original pure function form, despite no recursive
call path. Runtime checks remain enabled; successful debug stderr is empty.

The earlier three-mode campaign used n=16. The four-mode campaign retains both
communication plans; its n=16 forecast exceeded the unchanged 950 MiB ceiling.
The final matrix therefore uses three safe levels, not a raised memory limit.

The face tests now explicitly request `face_moments=.true.`. The runtime default
remains legacy interpolation. Their reference remains independently constructed
from integrated vector potential; constant and affine solenoidal fields must be
reproduced by the new fit to roundoff. The compact loops test truncation and
closure amplification rather than polynomial reproduction. See the progress
report for retained legacy/weighted/unweighted controls and the failed initial
optimized/debug comparison before canonical lattice-node selection was added.
