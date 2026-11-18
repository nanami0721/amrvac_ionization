# Central-gravity equilibrium milestone

Status: validated on 2026-11-17. Six complete debug library configurations,
three repeated debug cases, regressions, diagnostics/partition/restart invariance,
all separated 1/10-crossing equilibria on both grids, and radial refinement passed.
See `doc/yinyang_gravity_validation.md` for metrics and limits. Compilation uses
a user-authorized 4 GiB cap; simulations and analysis retain 1 GiB.

## Configuration

Use `equilibrium_case=0` for nonmagnetic hydrostatic gas (MHD with B=0), 1 for a
potential dipole in uniform gas without gravity, 2 for the combined case, or
3 for a uniform Cartesian field in uniform gas without gravity.
`beta_reference` (default 4) controls the field strength; `output_interval`
(default 1) is measured in sound-crossing times. See the adjacent low_beta suite.
`mhd_gravity` must agree with the case. The existing gravity callback is bound
before MHD activation. Its components are local to each patch: (-1/r^2,0,0).
Only central gravity is targeted by this milestone; arbitrary user gravity
callbacks must still provide correctly transformed local components.

Normalized parameters are R=GM=rho_base=1, c_T^2=0.5, gamma=5/3 and polar field
0.5 tilted by 20 degrees. The hydrostatic reference is rho=exp(2*(1/r-1)),
p=0.5*rho. Evolution retains the total-energy equation; no temperature forcing
or internal state resetting is applied. Initial polar beta is 4.

The shell is 1..2. Both radial ghost boundaries impose the analytic equilibrium
and zero velocity/psi; radial divB correction is disabled. This is an equilibrium
benchmark boundary, not strict face-flux locking or a validated open-corona BC.
The spherical control shifts longitude bounds to place the old Yin seam inside
the patch. CT, B0 splitting, AMR and other additional physics remain unsupported.

## Diagnostics

The optional `usr_fv_stage_diagnostic` hook receives read-only finite-volume
stage data before and after source addition. With no callback it changes no
state and allocates no diagnostic workspace. Flux values retain solver scaling
(-dt/dx for uniform Cartesian geometry, -dt*face area for curvilinear geometry).
The hook is for finite_volume, not every solver/predictor in AMRVAC.

This test restricts budgets to unsplit SSPRK(3,3), with stage multipliers
1/6,2/3,1 on the already time-integrated flux/source increments. It caches 16x16
angular overlap quadrature once, stores independent per-block accumulators,
and performs no MPI collective in the callback. OpenMP blocks have disjoint
accumulator indices. Budget rows record weighted mass/energy, cumulative net
radial inflows, gravity work, all-source energy and all-direction weighted
numerical flux contributions. The weighted flux includes overlap-weight
variations and angular edges: it must not be called physical radial flux.
Non-gravity source energy, including GLM effects, is reported separately.

At initialization, a nonzero-velocity source check verifies gravity's radial
momentum increment and energy work, unchanged transverse momenta, and the global
Cartesian acceleration on both frames. Native divergence/current and their
analytic discrete references are written per rank at each sound crossing (or
every step for debug). The diagnostic does not refresh or change the solution.
The on/off regression must establish that it leaves trajectories unchanged.

Error acceptance uses the common interior sector r=1.2..1.8, theta=60..120 deg,
phi=100..190 deg, with seam/off-seam subsets. The initial targets are Mach RMS
<0.01, Mach max<0.05, and relative density/pressure/B L1 errors <0.02. B=0 cases
use absolute magnetic RMS <1e-12. Full patch/global errors are also reported but
must not be directly compared between different physical domains.

## Execution and memory

Build only in an external source copy. `run.py` defaults to two Yin-Yang ranks
or one spherical rank. It automatically uses the shared memory guard: task-tree
RSS/physical-footprint limit 1 GiB; system memory must be below 65% to start and
below 68% while running. System occupancy sums anonymous pages, wired pages, and physical compressor
pages. File-backed cache is reported separately, along with swap usage.
Compilation alone is authorized up to 4 GiB; runs and analysis remain at 1 GiB. The 0.2 s watchdog is not a
kernel hard quota. Do not raise thresholds automatically or start concurrent
builds/runs/analysis. Guard failures and partial runs are retained.

`suite.py EXTERNAL_WORKDIR` requires successful external `prepare_build.py`
artifacts and executes debug, uniform/transport, diagnostics-invariance, MPI
partition, all three short equilibrium controls and (only after short acceptance)
long runs to 10 sound crossings. It compares restart to the continuous solution
at three sound crossings. A radial-only refinement is attempted only with large
measured memory headroom. Each operation has a separate guarded invocation.
Short-acceptance failures stop long evolution and require diagnosis, not relaxed
thresholds. The complete six-configuration build matrix must pass before this
suite can establish a validated milestone.
