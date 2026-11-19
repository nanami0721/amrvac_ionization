# Experimental Yin–Yang CT: compact magnetic-loop transport

This test advances full ideal MHD using native staged UCT-HLL, SSPRK3,
TVDLF, and van Leer reconstruction. It does not reset velocity, inject
analytic EMFs, or replace evolved magnetic fields. Only test files are added.

## Reference and path

Set rho=1, p=0.1, |v|=1, R=0.35, and peak analytic |B|=epsilon=1e-3.
For q=x-origin-v*t and s=|q|²/R², define A=a0*(1-s)^4*n for s<1,
and A=0 outside. Then B=-8*a0/R²*(1-s)^3*(q cross n).
Here a0=epsilon*R/[(8/sqrt(7))*(6/7)^3]. The actual initial face field
is the native discrete curl of A; centered B comes from native face averaging.
The initial discrepancy from pointwise analytic B is recorded, not subtracted.

The magnetic loop has finite Lorentz force. Rigid translation is an exact
kinematic induction reference, not an exact finite-amplitude full-MHD solution.
An epsilon/2 run checks magnetic back reaction; background discretization
errors in uniform Cartesian flow remain. This is not a low-beta test.

Let phi=146.25 degrees for the seam test and phi=0 for the interior control.
n=(cos(phi),sin(phi),0), v=direction*(-sin(phi),cos(phi),0), and
origin=1.5*n-0.3*v. Evolve t=0..0.6 with dt=0.002 (300 steps).
The center crosses the Yin phi=140.625-degree boundary at t=0.152263
(forward) or t=0.447737 (reverse); the Yang patch covers both trajectories.
The compact analytic support stays inside 1<r<2 throughout. Radial boundaries
use the currently supported continuous experimental CT policy.

## Matrix and memory

Coarse per patch: 16x24x72, blocks 8x12x12.
Fine per patch: 16x32x96, blocks 16x16x16 (fewer radial ghost cells).
Both directions and placements use both angular resolutions. Radial resolution
is fixed: these two levels do not establish a full 3D convergence order.
Rotating the path also changes its alignment to the local mesh; the interior
control can reveal a shared large diffusion error, but cannot uniquely measure
the interpolation contribution by subtracting the two error norms.
Additionally run amplitude/2, 4-rank coarse seam, and optimized coarse seam.
All cases use max_blocks=128. Final sorted text dumps compare both conserved
cell states and block-identified active staggered faces across 2 and 4 ranks.
This is not a snapshot/restart test.

Example (external directory must not exist):

```sh
/Users/nanami/miniconda3-arm64/envs/yt/bin/python3.12 tests/yinyang/ct_transport/run.py /private/tmp/amrvac-tests/ct-transport-UNIQUE --cache-source /path/to/validated/external/source
```

The runner copies current source, records hashes and Git diff/commit, invalidates
changed cached Fortran objects, and builds serially. Build tree limit: 4096 MiB;
runtime: 1024 MiB. System start/stop: 65/68%; 0.2-s sampled monitoring, not an OS
hard quota. Refinement forecast includes 25% margin and must be below 950 MiB.
`--pilot` runs just one debug coarse seam case. `--resume` checks the recorded
source manifest and resumes completed cases; runner revisions are saved separately.
Use a fresh directory after a failed launched case (guard output directories
are intentionally not overwritten). Histories and stderr are retained externally.

## Diagnostics

Every step checks finite states/EMFs, positive rho/p, three CT stages per step,
and max active native face-flux divB <1e-10. Passing this gate means stable,
divergence-preserving evolution; it is not a transport-accuracy acceptance gate.
Ghost and two-angular-cell seam-layer divB are also recorded. Successful stderr
must be empty under the debug runtime checks.

Error norms use Cartesian pointwise analytic B and volume quadrature, with
weight 1/2 wherever both patches cover a point and weight 1 elsewhere. This
partition of unity avoids simply counting overlap twice, but independently
sampled patch quadratures still have quadrature error. It is not evidence of
conservative overset flux exchange. The seam mask width changes with resolution.
Record L1/L2 relative to analytic-field norms, peak fraction, B² centroid,
velocity error, rho/p minima, and the reference energy fraction represented
by the Yin patch. `analyze.py WORKDIR --complete --plot` produces
`analysis.json`, verifies the completed matrix, and renders a summary plot.
`check_reference.py` independently checks the analytic curl by finite
differences, peak normalization, center crossing and radial support.

History columns (zero-based): 0 iteration, 1 time; 2 active divB max,
3 seam divB max, 4 all ghost divB max, 5 radial ghost divB max,
6 B error max, 7 velocity error max, 8 seam B error max, 9 interior B error max,
10 accumulated edge-integral magnitude max, 11 |B| max, 12 |rho-1| max,
13 rho min, 14 p min; 15 integral |error|², 16 integral |B_exact|²,
17 integral |error|, 18 integral |B_exact|; 19..21 integral x*|B|²,
22 integral |B|²; 23..24 seam squared error/reference integrals;
25 Yin weighted reference energy, 26 total weighted reference energy.


## Wider-structure follow-up

`--radius` selects the vector-potential support radius (default 0.35), while
`--radial-min` / `--radial-max` select the shell (defaults 1 and 2).
The runner checks that the translated analytic support remains strictly inside
that shell, and writes `geometry.json`. Resume requires matching geometry.
For a wider loop use `--radius .65 --radial-min .75 --radial-max 2.25`.
The cell counts, angular domain, speed, peak B, integrator and travel time stay
unchanged; widening the shell changes radial spacing. Therefore a matched
narrow-loop control must use the SAME widened shell.

`--seam-pair` runs only the debug forward-seam coarse/fine pair, useful for that
matched narrow control. Analyze it without `--complete` (which requires the
full 11-case suite). Then run `compare_widths.py NARROW_DIR WIDE_DIR`; this
checks that the paired parameter files differ only in radius and records the
relative error reduction in `width-comparison.json` in WIDE_DIR. The independent
reference checker accepts the same three geometry flags.

At r=1.5, the wider support diameter spans about 12 coarse / 16 fine theta
cells (versus 6 / 8 for the narrow loop), and 14 radial cells in the widened
shell (versus 7). These count the support, not the smallest gradient length.
Fixed travel time also means fewer loop widths traversed by the wider loop;
width comparisons diagnose sensitivity to spatial scale, not convergence order.
Only the center and parts of the magnetic structure cross the patch boundary
within the time interval; the entire broad support need not leave the Yin patch.


## Same-orientation native controls and error attribution

`run_controls.py WORKDIR --cache-source /path/to/ct-wide/source` runs a separate
17-case debug matrix, using the wide loop and shell. The cache parent must
retain `debug-d1-p1-s2-np2-a1-memory/summary.json` for the guarded YY forecast.
`--pilot` selects just the coarse native Yin orientation. All builds/runs are
external and serial, with the same 4096/1024 MiB and 65/68% protections.

Native control 0 uses the Yin basis with phi=46.875..234.375 degrees; its phi
spacing is exactly the original YY spacing. Control 1 uses the Yang basis as
a fixed proper rotation of spherical coordinates, theta=14.0625..115.3125 and
phi=-46.875..46.875 degrees. These are independent native spherical wedges,
not a full-sphere native run. Their boundaries are extended past the original
patch seams. Grid origins shift by integer numbers of cells, enabling direct
cell-to-cell comparison without interpolation of the comparison data.

For these controls, the fluid at angular physical ghost boundaries is set to
the uniform Cartesian reference (rho, p and velocity only). Magnetic boundaries
retain native `cont` face handling; radial boundaries remain `cont`. The analytic
magnetic support never reaches the physical boundaries. This avoids incorrectly
extrapolating constant spherical velocity components for a Cartesian flow.
No active velocities/fields are reset, and no analytic electric fields are
injected. The native physical boundaries differ from YY overset boundaries;
remote boundary and residual fluid truncation errors are still possible.

`save_samples=T` saves active states at t=0, .404 and .522, in addition to the
final state. The native reference frame is selected by `reference_patch`; it
is independent of the core YY patch identifier. The new native/sample flags
default off, preserving the existing YY test behavior.

Cases include native Yin/Yang coarse/fine, native reverse fine, YY forward and
reverse samples, and separate changes to dt/2 (native and YY), MC (`woodward`) reconstruction,
radial resolution 16->32, HLL flux choice, and HLL+MC. The setting changes are
experimental tests, not changes to production defaults. Reduced dissipation
must be assessed for local seam artifacts as well as global error.

Run `analyze_controls.py WORKDIR --complete --plot` under the memory guard.
It independently reconstructs native volume weights and exact Cartesian B from
snapshots and verifies the in-solver L2 histories. It also forms an uncoupled
native two-patch diagnostic on the exact same cells/overlap weights as YY,
reports the magnetic energy outside the common mask, and locates both native
and YY values at each YY peak. The combination is a diagnostic quadrature of
two independently evolved native wedges, not a coupled native physical solve.
Its gain/shape decomposition uses C=<Bnum,Bexact>/<Bexact,Bexact>, with
relative L2 squared = (1-C)^2 + ||Bnum-C*Bexact||^2/||Bexact||^2. Shape error
includes phase and distortion; norm differences are not an additive error budget.
The reader supports the E-less three-digit exponent form emitted by ES25.16.


A peak recovery is not automatically a spurious overshoot. The control analysis
also reports full vector error at the YY peak and the native value at exactly
that point, plus positive/negative contributions to the change in squared
error over the common mask. This separates local improvement from global net
error increase. It does not isolate edge exchange from ghost reconstruction.


## Native midpoint representation regression

`run_representation.py EXTERNAL_DIRECTORY --legacy-source EXTERNAL_OLD_SOURCE
--candidate-source EXTERNAL_NEW_SOURCE --reference-results EXTERNAL_CONTROLS`
uses two previously isolated source snapshots with their build caches. It records
source hashes, builds separate debug executables, and uses identical parameters
for TVDLF+van Leer and HLL+MC in both directions on the wide-loop fine grid
(16×32×96 per patch, dt=.002, t=.6). The reference-results argument points to
the same-orientation controls suite containing `yy-p0-s2-d1-hllmc` and its
memory record. It also repeats candidate HLL+MC forward transport at four ranks.
Builds/runs retain the established 4 GiB/1 GiB and 65%/68% memory protections.
The two source directories must be outside the development checkout.

Run `analyze_representation.py EXTERNAL_DIRECTORY --reference-results
EXTERNAL_CONTROLS` under `memory_guard.py`. It verifies unchanged initial samples,
reproduction of retained baseline runs, finite state, positivity, native active
divergence and bitwise MPI cell/face outputs. It reports both L2 and maximum
point error without asserting that either must improve. No performance speedup
is inferred from these single validation runs.

For face reconstruction, `run_representation.py --face-probe` runs only the
forward/reverse HLL+MC two-rank pairs. The legacy source should be the validated
midpoint-EMF version; candidate parameters additionally set
`yy_ct_face_moments=T`. Four-rank fine transport is omitted to retain memory
headroom; the explicit face component and analytic MPI tests cover 2/4/6 and
2/4 ranks respectively. The historical output label `midpoint-*` denotes the
candidate executable in this driver. Analyze with `--face-probe
--legacy-paired-results EXTERNAL_MIDPOINT_TRANSPORT` to compare with the correct
retained baseline. A face-probe result has `mpi_bitwise=null`, not an assertion
of an unrun evolved MPI comparison.

After introducing a configuration switch, `--face-probe --plumbing-only` selects
only default-forward and enabled-reverse transport. Validate those two outputs
with `analyze_representation.py --face-probe --plumbing-reference
EXTERNAL_FULL_FACE_PROBE --reference-results EXTERNAL_CONTROLS`; it requires
bitwise cell/face hashes and identical initial/final diagnostic rows to the
retained four-run campaign. This checks switch wiring, not a new accuracy study.
