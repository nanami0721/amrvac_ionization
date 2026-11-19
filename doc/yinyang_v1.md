# Experimental Yin–Yang grid v1

> For the current scope freeze and a fresh-machine handoff, read
> [Mac Studio handoff](../YINYANG_MAC_STUDIO_HANDOFF.md) first. Historical status
> notes below are dated evidence, not the current task queue.

> Recovery note (2026-11-16): the historical temporary artifacts cited below were lost during a reboot. The implementation was reconstructed from session edit records. Fresh accuracy evidence and its validation scope are in `yinyang_accuracy_diagnosis.md`; the historical six-configuration and timing statements below were not all rerun during recovery.

> Development status (2026-11-17): central-gravity support and the hydrostatic
> benchmark passed six-configuration builds, debug/regression/restart checks,
> and 10-crossing spherical/Yin–Yang equilibrium comparisons.
> See `yinyang_gravity_validation.md` for measured scope; open-corona physics remains unvalidated.

Select `call set_coordinate_system('yinyang')` in `usr_init`, before physics
activation. This is a distinct coordinate identifier; each patch reuses the
spherical metric and geometric-source implementation. The Yang frame is the
proper rotation `(x,y,z) -> (-x,z,y)` of the Yin frame.

> CT development update (2026-11-18): actual fixed-grid UCT-HLL integration is
> available only with `yy_ct_experimental=T` and `ct_stage_split=T`. See
> `yinyang_ct_stage_api.md` and `tests/yinyang/ct_static` for static-field tests,
> restrictions and measured errors. Default CT rejection remains intentional.

## Supported scope

- 3D3V, positive-radius shell; two identical fixed level-one meshes.
- Even MPI world size >= 2; contiguous equal groups own Yin and Yang. `mype`,
  `npe`, `icomm`, forest and block identifiers are **patch local** after startup.
  `yy_world_rank`, `yy_world_size`, `yy_patch` identify the global arrangement.
- Ideal total-energy MHD with GLM and FI ideal-gas EOS. Ideal HD and scalar
  transport are available for infrastructure tests.
- Uniform angular spacing and matching radial discretization. The supplied
  tests use TVDLF, minmod and three-step integration.
- Existing radial physical boundary routines; angular external ghosts are
  supplied exclusively by the other patch, including edge/corner ghosts.

AMR, automatic migration/load balancing, local timesteps, non-experimental CT, B0 splitting,
particles, multigrid, angular periodicity, partial-state exchange, extra/auxiliary
state and additional HD/MHD source modules other than optional MHD user gravity are rejected.
Central gravity has passed the equilibrium validation described below; coronal conduction/radiation remain unsupported. User source and boundary callbacks
must use the appropriate patch frame; no automatic transformation of arbitrary
user callbacks is possible. Scalar prescribed velocities are likewise local
components, not automatically a common Cartesian velocity on both patches.

## Exchange and costs

`mod_yinyang` gathers ownership once and caches donor indices, bilinear weights,
basis matrices and peer schedules. Runtime messages carry primitive states;
velocity and magnetic field are interpolated in the global Cartesian frame and
rotated into the receiver basis. Scalar primitives, including GLM psi, use the
same convex weights. Existing physics conversions restore conserved variables.

`getbc` operates on its actual Runge–Kutta state argument. Radial physical BCs
precede overset exchange. Donors must be interior in both angular coordinates;
missing stencils abort with an overlap diagnostic. Radial donors can include
already-filled physical ghosts. There is no recursive ghost-to-ghost angular
interpolation. Both overlapping interiors evolve; this is interpolation coupling,
**not flux-conservative overset coupling**. No global mass/energy repair is applied.

There are no runtime donor searches or per-cell allocations. Messages are batched
by peer using nonblocking point-to-point MPI; block request lists avoid scanning
all requests per block. Bilinear receptors use four samples; donor requests are deduplicated per peer
at schedule construction, and the deduplicated schedule is reused at runtime. Timers report the maximum rank's accumulated primitive
packing time and communication/wait time separately; these maxima may come from
different ranks, and waiting includes load imbalance. Receiver interpolation and
conversion are included in the existing ghost-cell timer, not those two subtimers.

## Output and restart

Output basenames gain `_yin` or `_yang`. Native snapshots retain local spherical
coordinates/components; the Cartesian conversion path rotates Yang positions and
vectors. Existing generic readers must not interpret Yang coordinates as Yin.
The supplied diagnostic reader reconstructs both frames explicitly. Slices,
collapsed output and arbitrary downstream analysis are not validated for a
combined global sphere.

New snapshots have a version-2 `.dat.yy` sidecar (patch, MPI size, mesh,
staggered/background-field flags, step, time and data-file size), plus a
`.dat.yy.complete` record. Before overwriting a snapshot, both patches invalidate
their old completion records. Completion records are published only after both
patch data files and metadata have closed. Restart requires both snapshots,
sidecars and completion records, unchanged MPI size, and the literal filename
placeholder:

```sh
mpirun -np 4 ./amrvac -if 'yy_{patch}0001.dat'
```

The standard snapshot reader checks geometry, extents and block shape; additional
checks reject inconsistent patch times, missing sidecars, incorrect patch/rank
metadata and stale timestamps. Version-2 reads reject missing/stale completion records, wrong field modes and
file-size mismatches. Legacy version-1 sidecars remain readable with a single
world-root notice; they retain their weaker historical guarantees. These checks
are not content checksums or crash-durable atomic filesystem transactions;
same-size corruption is not detected. Timestep, GLM wave speed, wall-time stop,
crash and `savenow` coordination use the world communicator.

## Running validation

All builds and runs must be outside the source repository. The runner copies
source (including untracked files and current edits), stores commit/diff evidence,
and points `AMRVAC_DIR` at the copied source. It requires NumPy, MPI and the
project compiler toolchain. Choose a new work directory:

```sh
python tests/yinyang/run_validation.py \
  --workdir /private/tmp/amrvac-tests/yy-UNIQUE \
  --extended --compile-matrix
```

The MHD `yy_test` namelist selects `initial_case=0` (Cartesian uniform flow/field),
`1` (advected density pulse in the overlap), or `2` (overlap pressure blast).
Separate `tests/yinyang/hd` and `rho` definitions exercise HD uniform flow and a
stationary scalar. Copy these cases to an external directory before using setup
and make. The scalar case checks exchange preservation, not advection accuracy.

`check_mhd.py SOURCE RUN_DIR SNAPSHOT` reports finite-state/pressure checks,
uniform-flow errors, near-edge error, a central-difference divergence estimate,
and weighted volume/mass/energy. The overlap partition is 1 outside the common
region and 1/2 within it, integrated with 16x16 angular subcells. Reported volume
bias makes quadrature error visible. Mass/energy changes also include radial
boundary fluxes; they are not a measured interface-only conservation defect.
The divergence estimate is not the solver's exact discrete divergence operator.
For the entropy case, the Python `diagnose(..., initial_case=1)` API uses the
translated analytic density pulse. Blast cases require positivity/stability
checks, not uniform-flow error assertions.

## Validation recorded on 2026-11-16

Source base: `6f5e6f10e5da5719dfc160da4591ae369450b603` plus this uncommitted
implementation, branch `codex/yinyang-grid-prep`. GNU Fortran 14.2 / local OpenMPI.
Detailed results/build copies are outside the repository at
`/private/tmp/amrvac-tests/yinyang-prep-20261115/`.

- Full modules and dependencies compile for 1D1V, **1D2V**, 1D3V, 2D2V, 2D3V,
  3D3V with debug bounds/FPE/temporary-array checks.
- 2/4/8-rank uniform-flow snapshots agree exactly; restart from step 5 to step 10
  agrees exactly with uninterrupted integration.
- Debug uniform flow runs 200 steps; scalar and HD smoke tests pass. Runtime stderr
  is empty for successful debug tests, including repeated exchange calls.
- Entropy pulse and pressure blast each run 50 steps on four ranks. Both retain
  positive density/pressure. These are smoke/stability tests, not shock accuracy
  or long-term coronal validation.
- At common time 0.026, grids per patch `(16,16,48)`, `(32,32,96)`, `(64,64,192)`
  give uniform-state combined L1 errors approximately `7.6e-5`, `2.2e-5`, `6.1e-6`.
  Maximum near-edge errors are `4.58e-4`, `1.90e-4`, `8.42e-5` respectively.
- Odd MPI size, AMR and CT are rejected. Changed restart MPI size, missing sidecars
  and inconsistent sidecar timestamps are tested failures.
- An optimized four-rank `(32,32,96)`-per-patch, 100-step run completes. This is a
  local timing observation, not a scaling result or a speedup versus spherical.

The v1 evidence supports a working fixed-grid experimental MHD implementation.
Production use still needs application-specific boundaries, resolution studies
for advected magnetic structure/waves and shocks, long-time divergence and
conservation budgets including radial fluxes, and dedicated scaling/memory
measurements. Strict conservative overset reconstruction is intentionally deferred.

## Guarded validation additions

`run_validation.py` uses serial GNU compilation under a 4 GiB process-tree cap;
MPI simulations and separate Python analysis jobs use 1 GiB. System accounting
uses anonymous+wired+physical-compressor pages with 65% start/68% stop thresholds.
Use a Python environment containing NumPy. Compilation and tests are sequential.
The extended runner refuses a refinement whose measured memory estimate exceeds
75% of the simulation cap; this is a resource failure, not a successful skip.
The guard also has a finite `--timeout-seconds` (default 3600); timeout and memory
termination preserve logs and return nonzero.

`checkpoint_faults.py VALIDATED_RUN EXTERNAL_WORKDIR` checks six metadata failure
paths after positive restart validation. Its MPI jobs have a 30-second timeout.
`overlap_diagnostics.py` reconstructs Cartesian magnetic fields from both patches
at identical sphere quadrature points. It reports the common interior sample
count and area, relative field differences, signed flux differences and absolute
normal-flux differences. These are overlap-consistency diagnostics, not a closed
sphere net flux or proof of globally conservative overset coupling.

`ct_preflight` is an ordinary spherical-grid prerequisite test using native CT
and a discrete-curl initialization. Passing it does not enable Yin-Yang CT.
