# Experimental Yin–Yang CT: static magnetic-field evolution

This is an actual fixed-grid MHD time-evolution test, not a probe-face exchange.
It explicitly enables `yy_ct_experimental=T` and `ct_stage_split=T` in methodlist,
with classical total-energy UCT-HLL, unsplit SSPRK3 and continuous radial boundaries.
The flag defaults to false. Non-CT physics, missing staging, incompatible CT
methods, background splitting, AMR, load balancing and internal boundaries are
not admitted by this experimental path. Other solar boundary models are not
implemented by this switch.

## Integration order

- Each RK substage computes all local face fluxes and edge integrals.
- The existing user CT-stage callback runs, then the native Yin–Yang adapter
  synchronizes shared edges, transfers owner overset edges (including radial
  endpoints), and resynchronizes shared copies before applying the curl.
- Ordinary same-patch ghost exchange runs first. In experimental mode the native
  radial boundary callback is replaced by the tested continuous radial policy:
  populate radial cell donors, interpolate angular cell states, complete ordered
  angular faces, close radial ghost faces, then synchronize center B/energy.
- The existing EoS refresh after getbc remains in place. No modifications to
  normal GLM or spherical CT behavior are made when the flag is disabled.

Plans are cached by mesh epoch, built from the actual block directory, and rebuilt
if the epoch changes. Active total energy is preserved; ghost thermal energy is
preserved when center B is replaced. The radial normal face is an anchor for ghost
closure, but its value evolves through the CT curl like other active faces.
Radial electric fields are not artificially zeroed. No velocity reset, field
reset, or equilibrium-restoring source is used in these tests.

## Static cases

The initial density and pressure are one and velocity is zero. Case iprob=0 has
zero B. Case iprob=1 has a uniform Cartesian target B=(0.2,-0.1,0.3), initialized
from the discrete curl of A=B cross r/2. Each patch rotates A into its local basis.
This gives initially divergence-free face fluxes; discrete center values need not
be exactly the analytic uniform vector. Both magnetic drift from the discrete
initial state and error relative to the analytic vector are therefore recorded.

A native spherical zero-field control runs with the same local mesh and time
step. Uniform pressure is not exactly well-balanced by the existing pointwise
geometrical source and finite-volume flux; velocity is therefore measured against
this control, not forced to zero or used as a roundoff-level acceptance criterion.

The short integration is 50 full steps (150 RK substages), dt=0.001, t=0.05.
An additional optimized coarse case advances 500 steps to t=0.5, with load-balance
timing diagnostics enabled. Native CT transfer time is excluded from the
compute-only counter alongside the user stage callback.
Radial mesh is eight cells over r=1..2; angular grids are 16x48 and 32x96 per patch.
Radial resolution is held fixed: this is an angular refinement check and does not
establish a general three-dimensional asymptotic convergence order. Uniform and
stretched radial coordinates are exercised. The runtime records active/seam/ghost
and radial-ghost native divB, maximum and RMS magnetic drift, seam and bulk drift,
spurious speed, peak unmodified EMF, face-field drift, analytic B error, minimum
rho and pressure, and exact RK callback counts. The first two angular cells
bordering a patch side define the seam region; the bulk additionally excludes the
first/last two radial cells. RMS uses local cell volumes and counts both patches'
overlap; it is not a unique-coverage global magnetic-energy integral.

History columns (whitespace-separated, one row per full step):

```
it time div_active div_seam div_ghost div_radial_ghost
B_drift_max speed_max B_drift_seam B_drift_bulk emf_max face_B_drift_max
B_analytic_max B_drift_rms B_analytic_rms rho_min pressure_min
```

Hard checks require finite states, positive rho/p, active divB <=1e-10, correct
stage counts, and zero-field magnetic drift below 1e-13. Magnetic
static-state drift is measured, not hidden by resets or an arbitrary loose pass
threshold. Canonical final states must match bitwise across tested MPI sizes.
Successful debug/optimized stderr must be empty. A successful short static run
does not establish long-time stability, magnetic transport, low-beta accuracy,
restart consistency or realistic solar boundary behavior.

## Reproduction and memory

```
python tests/yinyang/ct_static/run.py /absolute/new/external/directory --cache-source /absolute/external/source
```

All compilation and simulations use the copied external source. The suite retains
source hashes, tracked diff, parent commit, complete logs, parameters and memory
guard evidence. Builds are serial with a 4 GiB process-tree limit; runs use 1 GiB.
System thresholds are 65% start / 68% stop with 0.2 s sampling. max_blocks=128 and
a one-second initial hold expose allocated state. Refinement is gated using the
measured process-tree memory and the ratio of allocated block sizes (including
ghosts), with 25% margin and a 950 MiB predictive ceiling. For the four-rank
fine case, the estimate uses the measured two-rank fine allocation plus twice
the observed coarse-grid rank overhead, again with 25% margin. Guards are sampled,
not hard OS quotas. Debug libraries are checked in all six dimension/component
configurations. Explicit negative cases verify that CT does not run without its
opt-in switch, staging, or compatible edge/curl procedures.
