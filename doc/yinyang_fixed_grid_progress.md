# Fixed-grid Yin–Yang development progress (2026-11-18)

The full fixed-grid CT/open-corona plan is **not complete**. This increment
hardens validation and checkpoints, adds a common-position overlap diagnostic,
and establishes the ordinary spherical CT prerequisite and a two-pass RK-stage
CT API. Yin–Yang CT, static B0,
open-corona boundaries and Parker/dynamic-event validation remain disabled or
unimplemented. See the latest section for commit/push status.

## Source and evidence

Development worktree: `/Users/nanami/workspace/2608_uawsom/yinyang-recovery-20261116/source`,
branch `codex/yinyang-grid-prep`, base `6f5e6f10`, with preserved pre-existing edits.
External artifacts: `/private/tmp/amrvac-tests/yy-fixed-20261117`.
`baseline-manifest.json`, `baseline.diff` and `baseline-status.txt` preserve the
starting state. `verified/provenance.json` records current core-template hashes
and the content-checked cache invalidation: only `mod_yinyang.t` and
`mod_input_output.t` differ from the gravity milestone's core source.
`verify_current.py` and `regressions.py` in the external directory record exact
build/run commands. Temporary artifacts should be archived externally before a
reboot if long-term retention is required.

## Implemented

- Serial builds guarded at 4 GiB; MPI and separate analysis subprocesses at
  1 GiB; system 65% start and 68% stop. A finite timeout stops verified descendants
  as well as the task process group. NumPy availability is checked before builds.
- Direct transport-runner invocation now activates the memory watchdog.
  Extended validation estimates refinement memory from an observed coarse run
  and refuses an unsafe refinement before launching it.
- Version-2 snapshot sidecars record field representation and file size.
  Completion records are invalidated before overwrite and published after both
  patches close their snapshot/metadata files. Reads reject missing/stale records
  and metadata mismatches. Legacy v1 remains readable with an explicit notice.
  Completion is not fsync-backed transactional durability or a content checksum.
- Independent magnetic reconstructions are compared at the same physical sphere
  sample points. Reports include stencil coverage, area, relative field mismatch,
  signed normal flux difference and absolute normal flux difference. They do not
  measure a closed-sphere net flux or establish global overset conservation.
- Added a native spherical UCT-HLL prerequisite test using discrete-curl initial
  fields and the solver's staggered get_divb operator; no Yin–Yang CT enablement.

## Historical verification (2026-11-17)

The temporary evidence directory below is no longer present. These entries are
retained as historical records; the fresh validation section records checks
repeated against the current source in persistent external storage.

- Six full dependency-library debug configurations passed, including 1D2V.
  Optimized and debug 3D uniform and native-CT cases compiled and ran.
- Uniform 2/4-rank and continuous/restart snapshots matched exactly. A final
  executable comparison also verified legacy-v1 restart, with zero snapshot
  difference. Runtime stderr was empty on successful simulation cases.
- Six checkpoint failure cases passed: missing completion, stale completion,
  wrong data size, wrong field mode in sidecar, unknown metadata version and
  missing metadata. Each terminated in under 2 seconds without watchdog timeout;
  failure logs were bounded. See `checkpoint-faults/results.json`.
- Six three-step debug equilibrium cases (three fields x two geometries) passed.
  This does not repeat the historical ten-crossing equilibrium acceptance.
- Forward/reverse magnetic transport completed. Initial comparison used different
  limiter/diagnostic settings and was not accepted as a regression comparison.
  With the exact historical parameter files, both final snapshot differences
  were zero; see `matched-baseline-comparison.json`.
- Native spherical CT, 1/2 ranks, five three-stage steps: maximum raw divergence
  7.138447277875828e-15 (debug) and 5.13696499869064e-15 (optimized).
  Each run produced 11 diagnostic records. This is a short single-patch invariant
  test, not proof of overset magnetic consistency or long-time accuracy.
- Analytic Python diagnostic tests passed: global uniform-vector preservation,
  invalid-stencil masking and smooth-field reconstruction refinement.
- Uniform GLM overlap relative field mismatch at r=1.25,1.5,1.75 was
  3.278889150475603e-4, 2.758660446397004e-4, 2.380939750065716e-4.
  Each sphere used 1120 common samples from a 64x128 quadrature. These are baseline
  values at one resolution, not an interface-convergence result.
- Startup rejection, process-tree memory termination and timeout termination
  watchdog checks passed. Maximum sampled system occupancy was 51.072%; maximum
  build RSS 1656.125 MiB, simulation/analysis RSS 537.25 MiB; no unexpected
  watchdog stop. Sampling is not a kernel-enforced hard memory quota.

The first baseline analysis attempt used a Python without NumPy and failed;
subsequent checks used the existing scientific Python environment, and the runner
now rejects that environment error before building. Failed-attempt logs remain.

## Fresh staged-CT validation (2026-11-18)

Persistent external evidence:
`/Users/nanami/workspace/2608_uawsom/yinyang-ct-stage-20261118`.
`manifest.json`, `commit.txt` and `source.diff` identify the tested source;
`cache-invalidation.json` records changed-template invalidation in the copied
external build tree. All 173 recorded source/test hashes still match the working
source after validation. Existing source changes were preserved.

The opt-in two-pass API now prepares all blocks' UCT-HLL edge integrals, invokes
one `usr_ct_stage` callback per rank and RK stage, then applies the curl and
original finite-volume/source sequence. The default update path is unchanged.
See `yinyang_ct_stage_api.md` for restrictions and storage lifetime.

- Debug and optimized builds each passed 19 cases, including one intended
  unsupported-contact-CT rejection, and 13 complete-snapshot equality comparisons.
  Coverage includes static/nonuniform moving fields, 1/2/3 MPI ranks, uneven
  block ownership, TVDLF/HLL and a nontrivial stage-edge scaling operation.
- Maximum native raw divergence was 6.37993379755453e-15 (debug) and
  6.713897591972269e-15 (optimized), in the test's units. This is a short native
  spherical CT invariant check, not a Yin–Yang interface result.
- Full dependency libraries compiled in all six debug configurations:
  1D1V, 1D2V, 1D3V, 2D2V, 2D3V and 3D3V.
- Staged dynamic CT at three ranks: ten continuous steps and restart from step
  five produced byte-identical final snapshots. Stage workspace was rebuilt.
- Existing Yin–Yang GLM uniform snapshots were byte-identical at 2/4 ranks and
  matched the retained 2026-11-17 gravity-milestone baseline for both patches.
- All successful simulation stderr files were empty. `git diff --check` passed.
- Sampled system occupancy peaked at 48.674% of 64 GiB, using
  anonymous+wired+physical-compressor accounting. Build process-tree RSS peaked
  at 1627.203 MiB; simulation/analysis trees at 379.938 MiB. No watchdog stopped
  a task. Limits remained 4096/1024 MiB with 65% start and 68% stop thresholds;
  0.2-second sampling is not a kernel-enforced hard quota.

Evidence files: `status.json`, both `stage-suite-*/results.json`,
`followup_checks.py`, `followup-results.json`, per-job memory summaries and
`final-audit.json`. The retained external follow-up script contains exact restart
and GLM commands. No performance improvement or overset CT acceptance is claimed.

## Remaining dependency order

1. M1 face/edge schedules and analytic tests are implemented (see the staged
   sampling section below). Remaining M1 work: cross-resolution common masks
   for physical diagnostics and physical whole-sphere budgets.
2. Complete M2 integration: the explicit edge-integral exchange and angular
   ghost-face constraint APIs are now implemented and tested below. Remaining
   work is integrating the explicit magnetic ghost fill below with duplicate-edge
   synchronization, radial boundary intersections and native face updates.
   The new opt-in `ct_stage_split` route now supplies an all-block boundary
   between edge generation and curl application. The original route remains the
   default. Exchange must use that stage boundary, never a collective inside a
   per-block hook. Do not enable CT by removing yy_validate's rejection or by
   overwriting active face fields after interpolation.
3. Verify overset native divergence, shared-block edges, overlap convergence,
   energy changes, MPI and CT restart before extending physical models.
4. Add radial magnetic/characteristic boundaries and analytic static B0 support.
5. Verify Parker wind, split-monopole wind and a controlled shock-producing event.

The current GLM rejection of CT/B0 remains intentional. M0's small regression
baseline and parts of M1/M4 are available; M2 has the staged update API and tested explicit edge-exchange/ghost-constraint
components. Full MHD CT evolution remains disabled.

## Face/edge schedules and analytic tests (2026-11-18)

`mod_yinyang_staggered` now caches component-specific donor stencils, canonical
block ownership, MPI routing, Cartesian reconstruction coefficients and reusable
buffers for face-normal and edge-tangent point samples. The AMRVAC adapter
`yy_build_staggered_table` uses the real block directory and `ps%x/ps%dx` radial
mesh, including radial stretching. Receiver points may also be supplied directly
for later edge quadrature. Invalid coverage is explicit; no extrapolation or
active-face updates are performed. This API is not connected to CT evolution.

See `tests/yinyang/staggered_schedule/README.md` for indexing conventions,
reconstruction equations, rebuilding/lifetime rules and exact test commands.
The normal-matrix normalization preserves uniform Cartesian vectors; it does
not impose a divergence constraint or global flux conservation. Donor requests
are currently not deduplicated, so production communication scaling is untested.

Fresh analytic evidence is retained in
`/Users/nanami/workspace/2608_uawsom/yinyang-staggered-20261118-c`.
Debug and optimized suites each completed 36 MPI runs: face/edge x uniform/
geometric radial spacing x three resolutions x 2/4/6 ranks. Both patches and all
three component orientations were covered. Canonically ordered reconstructed
outputs were bitwise identical across MPI sizes. Repeated exchange, table
release/rebuild, empty-receiver ranks and radial/angular no-coverage masks passed.
The uniform Cartesian field maximum error was 6.11e-16.

For the smooth divergence-free analytic field `(sin(y),sin(z),sin(x))`, RMS errors
on the staggered boundary rings/corners at `(8,16,48)`, `(16,32,96)` and
`(32,64,192)` cells per patch were:

| Sampling lattice | Radial mesh | Coarse | Medium | Fine | Orders |
|---|---|---:|---:|---:|---|
| Face | uniform | 2.261e-3 | 5.657e-4 | 1.385e-4 | 1.999, 2.031 |
| Face | geometric | 2.227e-3 | 5.504e-4 | 1.339e-4 | 2.016, 2.039 |
| Edge | uniform | 2.307e-3 | 5.702e-4 | 1.384e-4 | 2.016, 2.042 |
| Edge | geometric | 2.265e-3 | 5.540e-4 | 1.338e-4 | 2.032, 2.050 |

These are spatial reconstruction errors on mesh-associated receiver rings, not
time-evolution errors, native discrete divB or proof of edge-integral consistency.
One intermediate diagnostic string exceeded VACPP's line-breaking capability;
that expansion failed before compilation. The string was shortened and the
complete final analytic suite rerun; failed evidence remains in the `-b` folder.

Final AMRVAC adapter evidence:
`/Users/nanami/workspace/2608_uawsom/yinyang-staggered-mesh-20261118-b`.
Twelve runs (debug/optimized x uniform/stretched x 2/4/6 ranks) each checked both
sampling families on 192 actual mesh receiver positions. Maximum discrepancy
between cached and native geometry/uniform-vector reconstruction was 1.56e-15.
All six full debug dependency-library configurations passed. All successful
simulation stderr files were empty, including repeated exchange calls.

Final source/test code hashes match both retained manifests (documentation was
updated after testing). `final-audit.json` records acceptance and memory peaks:
system 49.589% of 64 GiB, build tree 378.797 MiB, runtime tree 379.641 MiB; no
watchdog stop. All compilation, logs and simulation artifacts are external.
No commit or push was performed for this increment.

The point-sampling table was the prerequisite for the next increment below.
It alone did not justify removing the CT rejection.

## Edge-integral exchange and ghost constraints (2026-11-18)

Implemented `mod_yinyang_ct_coupling`: cached two-point Gauss integration along
spherical receiver edges, conversion of native donor `dt*length*E` to `dt*E`,
and theta-then-phi outward angular ghost-face flux constraints including corners.
The latter preserves every active face and its boundary anchors, validates all
used faces/metrics before writing and does not modify radial physical ghosts.

The actual-stage-workspace adapter consists of `yy_ct_stage_plan`,
`yy_prepare_ct_edges` and `yy_exchange_ct_edges`. It checks computed UCT stencil
bounds, missing workspace, stale mesh epoch and incomplete donor coverage. Donors
are all packed before receiver overwrites. The common mesh-directory helper was
refactored and the previous face/edge table integration suite rerun successfully.
See `tests/yinyang/ct_coupling/README.md` for contracts, equations and commands.

Final evidence directories:

- Analytic: `/Users/nanami/workspace/2608_uawsom/yinyang-ct-coupling-20261118-a`.
- Native adapter: `/Users/nanami/workspace/2608_uawsom/yinyang-ct-coupling-mesh-20261118-b`.
- Previous table regression: `/Users/nanami/workspace/2608_uawsom/yinyang-ct-coupling-table-regression-20261118`.

Acceptance:

- 72 analytic MPI runs (debug/optimized, midpoint/exact-line-average donors,
  uniform/geometric radii, three refinements, 2/4/6 ranks) passed. Canonical
  edge-integral outputs were bitwise identical across MPI sizes. Repeated calls,
  zero dt, linear dt scaling, incomplete coverage and empty receivers passed.
- Smooth-field line-integral errors normalized by dt and receiving arc length
  converged at 1.880--2.109 order; fine-grid pairs were 2.080--2.109. The reference
  uses independent 128-panel Simpson integration. Uniform Cartesian E with native
  midpoint donors had maximum normalized error 2.08e-8 from curved-edge
  quadrature. Exact-line-average donors have second-order midpoint-conversion
  error; they are not claimed to preserve constant Cartesian E to roundoff.
- All 16 angular-side flag combinations passed ghost divergence, corner,
  idempotence, exact protected-face preservation and invalid-area atomic-failure
  tests. Analytic spherical flux-divergence maximum was 6.69e-15.
- 12 real AMRVAC mesh runs passed. Maximum donor-length discrepancy versus native
  dsC: 3.47e-16; receiving-integral error versus the uniform-field two-point
  quadrature: 8.56e-17; constrained ghost-cell native get_divb: 7.99e-15. The
  ghost operator was applied to a separate probe state with native surfaceC and
  volumes; it did not mutate a production MHD solution.
- 12 previous-table regression runs passed. All six full debug dependency-library
  configurations compiled. Successful runtime stderr was empty, including
  repeated hot-path calls. Final source/test code hashes match retained manifests;
  `git diff --check` and Python syntax checks passed.
- Final-suite sampled system occupancy peaked at 54.875% of 64 GiB; build/runtime
  process-tree RSS peaks were 123.625/219.907 MiB. No watchdog stop occurred.
  Limits stayed at 4096/1024 MiB and system 65% start/68% stop, sampled at 0.2 s.
  The native directory contains `final-audit.json` and per-job memory evidence.

These tests qualify explicit transfer/constraint components, **not complete
Yin–Yang MHD CT time evolution**. No callback is installed automatically. Remaining
integration must populate magnetic ghost faces, order same-patch ghost/duplicate
edge synchronization, handle radial-wall intersections, synchronize cell-centered
B/energy and validate evolving fields, restart and interface convergence. The
CT/B0 rejection remains; no commit/push was made for this increment.

## Magnetic ghost-face fill (2026-11-18)

The preceding edge-integral/ghost-constraint increment was committed locally as
`0b828120` at the user's request. It was not pushed. The magnetic-fill changes
below remain uncommitted; no push was performed in this task.

`yy_prepare_ct_ghosts` now constructs a cached union of all required angular
magnetic ghost faces, excluding active faces and active-boundary anchors.
`yy_fill_ct_ghosts` packs both donor families before any writes, fills the face
values and optionally applies the outward divergence constraint. Physical radial
ghost cells and unrelated same-patch ghost regions are untouched. The lower
outer face is located from its adjacent center, avoiding undefined x index zero.

Radial wall support uses Br as a patch-invariant scalar on its own face lattice;
this avoids requesting out-of-domain radially centered tangent components.
Tangential faces retain the Cartesian-vector reconstruction. Both radial walls
are exercised with interior/boundary face donors only. A roundoff tolerance was
added to radial endpoint bracketing; missing physical stencils still fail.
See `tests/yinyang/ct_ghost_fill/README.md` for the API and indexing contracts.

Acceptance evidence:

- 36 native-grid MPI runs passed: debug/optimized x uniform/stretched radial mesh
  x three angular resolutions x 2/4/6 ranks. Every run tested constant Cartesian
  B and divergence-free B=(y,z,x), with/without the constraint, and repeated each
  fill three times. All non-receiver faces were bitwise preserved, no ghost donor
  was requested, repeated fills were bitwise reproducible and canonical results
  were bitwise identical across MPI sizes. Unfilled halos were poisoned by a
  large sentinel to expose unintended donor reads.
- Maximum constrained ghost-cell native get_divb was 1.141e-14. This says nothing
  about repairing pre-existing active-cell divergence; those faces are unchanged.
- Angular RMS refinement orders ranged from 1.932 to 2.266 over 16x48, 32x96 and
  48x144 angular grids. Radial resolution was held at eight cells, with two
  radial blocks. The refinement ratios 2 and 1.5 were used explicitly.
- The constraint has a measurable accuracy cost. For B=(y,z,x) on the uniform
  radial mesh, RMS errors without/with the constraint were:

| Angular grid | Without constraint | With constraint |
|---|---:|---:|
| 16x48 | 2.920e-3 | 4.994e-3 |
| 32x96 | 7.651e-4 | 1.038e-3 |
| 48x144 | 3.449e-4 | 4.472e-4 |

  Both sequences decrease at approximately second order; the constrained fine
  error is about 30% larger. These are point-component analytic tests on a
  staggered mesh, not a general 3D asymptotic-order result for area-averaged
  magnetic fields or a full evolving MHD problem.
- The original face/edge table suite (72 MPI runs) and edge-integral/constraint
  suite (72 runs) also passed against the final source. All six full debug
  dependency libraries compiled. Successful runtime stderr remained empty.
- Final source/test code hashes match retained manifests; Python syntax and
  `git diff --check` passed. Final validation system peak was 52.308% of 64 GiB;
  build/runtime process-tree RSS peaks were 364.547/375.563 MiB. Including earlier
  attempts with the default block reserve, runtime peak was 770.016 MiB. No
  watchdog stop occurred and the 4096/1024 MiB, 65%/68% limits were unchanged.

Evidence directories:

- `/Users/nanami/workspace/2608_uawsom/yinyang-ct-ghost-20261118-d`
- `/Users/nanami/workspace/2608_uawsom/yinyang-ct-ghost-table-regression-20261118`
- `/Users/nanami/workspace/2608_uawsom/yinyang-ct-ghost-edge-regression-20261118`

The first attempt's stretched case was rejected because its radial block count
was odd. The next attempts' refinements were blocked by the predictive memory
gate before launch. Final tests use two radial blocks, max_blocks=128 (maximum
54 owned blocks) and a one-second hold to expose short jobs' allocated state to
the 0.2-second watchdog sampling. The upper angular resolution was reduced from
64x192 to 48x144 while retaining three refinement levels; failed/blocked attempt
records are retained separately in the `-a`, `-b`, `-c` directories.

This completes the explicit magnetic ghost-face population plus constraint
component. Full Yin–Yang MHD CT remains disabled until same-patch/duplicate-edge
synchronization, radial boundary ordering and cell-centered B/energy conversion
are integrated and time-evolution/restart accuracy is verified.


## Same-patch CT synchronization and mixed-corner ordering (2026-11-18)

The explicit `yy_ct_patch_plan` now caches shared-edge copies, physical-domain
face halo copies and extended angular-halo copies. Geometrical edge ownership
uses the lower block at each shared transverse boundary, including four-block
corners; it does not depend on rank numbering. `yy_ordered_ct_edges` synchronizes
owners before overset integration and distributes their updated values afterwards.
It rejects nonowner cross-patch receivers before exchange. No EMF averaging is
introduced. `yy_ordered_ct_ghosts` runs ordinary face halo copies, cross-patch
population/closure, then mixed-corner completion. Active magnetic faces and
physical radial ghost cells remain untouched. This is a fixed level-one API;
there is no automatic production-hook registration.

The new test definition is `tests/yinyang/ct_patch_sync`. Retained evidence:
`/Users/nanami/workspace/2608_uawsom/yinyang-ct-patch-20261118-e`.

- All 24 debug/optimized cases passed: two block layouts, uniform/stretched
  radial grids and 2/4/6 MPI ranks. Canonical magnetic face outputs were bitwise
  identical across MPI sizes within each geometry/layout/build.
- Deliberately inconsistent duplicate EMFs became bitwise identical to their
  geometrical owner, including four-block edges. Repeated synchronization was
  idempotent. Invalid edge slots and face flux arrays were unchanged.
- Constant Cartesian E integration after corrupting nonowners had maximum error
  6.592e-17. Ordered exchange left every shared edge copy exactly equal to its
  source. Two invalid-owner cases failed with the intended diagnostic, without
  a watchdog stop or timeout.
- Poisoned mixed-corner face holes were removed by the third phase: 14080 holes
  for 4x8x8 blocks and 9088 for 4x16x16 blocks before that phase. Constant and
  linear Cartesian magnetic fields were tested. Repeated ordered calls matched
  the manually sequenced reference bitwise and preserved all active faces and
  radial physical ghosts. Native angular-ghost divergence, including mixed
  corners, was at most 1.071e-14.
- All successful runtime stderr files were empty. All six full debug
  dimension/component dependency libraries compiled, including 1D2V. Source/test
  hashes matched the retained manifest; Python syntax and git diff checks passed.
  Final build/runtime process-tree RSS peaks were 357.266/243.547 MiB; system peak
  was 51.730% of 64 GiB. Across all attempts the system peak was 51.912% and
  process-tree RSS stayed below 392.188 MiB. No watchdog stop occurred. The
  unchanged limits were 4/1 GiB, 65% start / 68% stop, with 0.2 s sampling. These
  are sampled measurements, not hard OS quotas. See `final-audit.json`.

Earlier attempt records are retained in `-a` through `-d`: a test declaration
needed whitespace for VACPP expansion; the default iprob selected a negative
case before explicit normal-case configuration; iprob then needed its meshlist
placement; and the negative-case harness used the wrong watchdog summary key.
These were test-definition/harness corrections, not relaxed numerical tolerances.

This completes the explicit same-patch synchronization/ghost ordering component.
Radial-wall edge donor closure, physical radial boundary ordering, cell-centered
B/energy consistency and production integration remain before enabling full
Yin–Yang CT. The present tests use native geometry and probe face states, not
an evolving MHD solution, and do not establish CT evolution or shock accuracy.


## CT radial boundaries and B/energy synchronization (2026-11-18)

The explicit boundary/variable completion component is implemented:

- `yy_prepare_ct_edges(...,radial_extrapolate=.true.)` provides one-sided linear
  E_r reconstruction from the nearest two physical radial samples. Radial-wall
  tangent edges and first/last-cell radial edge quadrature are now covered without
  ghost donors. Sampling outside the physical shell is still rejected, angular
  coverage is unchanged, and default calls retain strict non-extrapolating behavior.
- `yy_complete_ct_faces` finishes ordered angular face communication, then applies
  constant radial extrapolation of tangential B and zero-divergence radial ghost
  normal-flux closure. Wall-normal Br anchors remain unchanged. The radial sweep
  runs last over the full angular halo, including all physical mixed corners.
- `mhd_ct_sync_center` reuses native area-weighted face-to-center interpolation.
  Active total energy is unchanged; ghost energy receives the magnetic-energy
  difference to preserve thermal energy. Background splitting, semirelativistic
  MHD and non-total-energy configurations are rejected by this helper.
- `yy_sync_ct_state` first synchronizes active cells and angular ghosts, then
  copies consistent boundary cell variables to physical radial ghosts before
  synchronizing their B/energy. Its non-overlapping slabs visit each ghost once.
  The test varies density, momentum and the additional GLM variable spatially to
  verify the radial source location and protection of other cells.

This is an explicit continuous radial/zero-divergence ghost policy. It differs
from native bc_cont's propagation of neighbouring active divergence. It is not a
characteristic wind outflow, prescribed wall electric field or line-tied boundary.
No active energy correction, pressure floor or low-beta improvement is claimed.
No production callback registration or CT startup restriction has been changed.

The reproducible definition is `tests/yinyang/ct_boundary_sync`; final boundary
results are retained at
`/Users/nanami/workspace/2608_uawsom/yinyang-ct-boundary-20261118-e`.
The full six-dimensional debug dependency rebuild is retained in the preceding
`-d` directory (including complete MHD modules and 1D2V); final validation reuses
those exact-source libraries and verifies all six library targets again.

The suite covers two block layouts, uniform/stretched radial grids, 2/4/6 ranks,
and debug/optimized builds. It tests constant and radially linear Cartesian E,
constant and linear Cartesian B, all physical face halos, nonuniform conserved
cell data, native divergence, active energy preservation, ghost thermal energy,
protected normal boundary anchors, and repeated-call idempotence. Canonical B/E
outputs are required to match bitwise across MPI decompositions. Nonowner edge
receivers and unsupported B0 energy synchronization have explicit negative tests.

An optimized comparison of native face-to-center applied to the full allocation
versus separate active/ghost boxes differed by 1.1–4.4e-16 due to floating-point
loop optimization. Bounded diagnostics in the `-c` directory isolate this from
state drift: repeated ordered synchronization stayed bitwise identical. The
native-operator equivalence test now permits 8 machine epsilons times the field
scale; all MPI, repeat, active-energy and protected-face bitwise checks remain.
Other early attempts corrected a VACPP array-constructor expansion and a missing
physics-procedure import in the test. Failed attempts and diagnostic logs remain
external; no failure was hidden by suppressing stderr.

Production stage/getbc integration is the next step, followed by actual CT static
field, transport, seam-crossing convergence and restart tests. Boundary-model
accuracy for solar-wind or driven-corona runs requires separately selected
physical boundary conditions. Present success verifies explicit operations on
native GLM-MHD geometry with probe staggered states, not evolving Yin–Yang CT.


Final boundary validation and regression results:

- 24 boundary/variable cases and 4 intended rejection cases passed; the previous
  same-patch suite passed 24 cases plus 2 intended rejections, and the original
  integral/constraint suite passed 72 cases with its refinement checks. Total:
  120 successful numerical runs and 6 correctly rejected invalid configurations.
- Maximum edge-integral error was 1.041e-16; physical ghost native divB was
  1.071e-14, including radial/angular corners. Maximum ghost thermal-energy
  difference was 8.882e-16; native center-B comparison differed by at most
  4.441e-16. Active energy and face anchors remained bitwise unchanged, repeated
  operations were bitwise idempotent, and MPI-decomposition hashes matched.
- Successful runtime stderr was empty. All source/test manifests matched the
  final workspace. Full original VACPP procedures and six compiled dimension
  libraries were checked; Python syntax and git diff checks passed. Compiler:
  GNU Fortran 14.2.0, serial debug/optimized builds, one OpenMP thread per rank.
- Across all attempts, builds and regressions, the system memory peak was 56.881%
  of 64 GiB; build/runtime process-tree RSS peaks were 1447.766/284.750 MiB.
  No watchdog stop occurred. Limits remained 4096/1024 MiB, system 65% start and
  68% stop, sampled every 0.2 s. These are sampled process-tree measurements.

Additional retained regression directories:

- `/Users/nanami/workspace/2608_uawsom/yinyang-ct-boundary-patch-regression-20261118`
- `/Users/nanami/workspace/2608_uawsom/yinyang-ct-boundary-edge-regression-20261118`

The final boundary directory contains `final-audit.json` with the aggregated
checks and memory peaks. No commit or push was made for this increment.


## Actual CT static-field evolution (2026-11-18)

The fixed-grid stage and ghost adapters are now connected to the real MHD
integrator behind `&methodlist yy_ct_experimental=T, ct_stage_split=T`. Default
GLM behavior and native spherical CT keep their original paths. Validation
requires classical total energy, UCT-HLL, unsplit SSPRK3, fixed 3D mesh and
continuous radial conditions. Existing restrictions on other Yin–Yang physics
are retained. The new physics synchronization callback is registered only for
supported classical total-energy CT. Plans are cached by mesh epoch.

At each stage, the user hook runs first, then owner/overset/shared-copy EMF
synchronization precedes the native curl. At getbc, ordinary same-patch exchange
precedes radial cell donor copying and angular primitive interpolation; complete
face halos and center B/energy synchronization follow. The old physical radial
callback is bypassed only in the experimental mode. The existing EoS refresh
after getbc remains. No active velocity, magnetic field or EMF is reset to keep
the test static.

`tests/yinyang/ct_static` now evolves uniform rho=p=1, v=0 and either zero B or
Cartesian B=(0.2,-0.1,0.3), initialized by discrete vector-potential curl. Tests
include native spherical zero-field controls, 2/4 ranks, debug/optimized builds,
uniform/stretched radii, two angular meshes and an additional 500-step coarse
run. Short runs have dt=0.001, 50 steps, t=0.05; the longer run reaches t=0.5.
The exact expected RK callback counts are checked: 150 and 1500 respectively.

Uniform pressure is not exactly stationary even with B=0. The original attempted
roundoff-level zero-velocity assertion caught a velocity of 2.46666e-6 after one
step. The added native spherical control gives 2.46811e-6 at the same step; after
50 steps their maxima are 1.16604e-4 (Yin–Yang) and 1.19537e-4 (spherical).
This is consistent with the existing pointwise spherical geometrical sources
not exactly cancelling finite-volume pressure fluxes. The zero-field magnetic
state remains exactly zero. The velocity assertion was replaced by this explicit
control and quantitative reporting, not a velocity reset or an equilibrium fix.
This control does not isolate every source of nonzero-field magnetic drift.

Optimized uniform-radius results, B drift relative to the discrete initial state:

| Cells per patch | Steps / time | B RMS drift / target field strength | Maximum B drift / field strength | Maximum active divB at final time |
|---|---:|---:|---:|---:|
| 8x16x48 | 50 / 0.05 | 0.04180% | 0.30558% | 1.645e-14 |
| 8x32x96 | 50 / 0.05 | 0.01407% | 0.14390% | 3.338e-14 |
| 8x16x48 | 500 / 0.5 | 0.34722% | 2.154% | 1.033e-13 |

The target field strength is sqrt(0.14). The RMS drift decreases by a factor
2.971 on angular refinement, but radial resolution is fixed and there are only
two angular resolutions: no general convergence order is claimed. The volume
weights count both patches' overlap. For the short coarse run the seam maximum
B drift is 8.703e-4, versus 9.390e-4 in the selected bulk; on the finer angular
mesh these are 3.951e-4 and 4.461e-4. The seam region is the first/last two angular
cells, while bulk also excludes two cells at each radial wall.

The 500-step run reaches active/ghost maximum divB 1.034e-13 over its entire
history. Final rho_min=0.999156, p_min=0.999050, and max speed=0.001640. This
verifies the CT divergence invariant through actual evolution with positive
states, but the measured equilibrium drift remains a real accuracy limitation.
It is not evidence of long-time solar-corona, low-beta, transport or restart
qualification. Analytic-field error is also recorded separately from drift to
avoid confusing initial discrete-curl/center reconstruction error with evolution.

Retained final static suite:
`/Users/nanami/workspace/2608_uawsom/yinyang-ct-static-20261118-f`.
The preceding `-e` directory retains the full dependency rebuild and original
static results before the transfer-timer correction.
Earlier attempts are retained in `-a` through `-d`: a VACPP-unbreakable diagnostic
string, the old MHD GLM-only validator, the zero-velocity assumption diagnosed
with a native control, and a conservative four-rank refinement prediction gate.
The four-rank fine task was not launched at that gate. Its revised estimate used
measured two-rank fine memory plus twice the measured coarse rank overhead and
25% margin; the 950 MiB predictive ceiling and all watchdog limits were unchanged.


Final static-test audit:

- 15 actual static-evolution cases and 6 expected configuration rejections passed.
  All successful runtime stderr files were empty. Two/four-rank canonical final
  **cell-centered conserved states** matched bitwise. Native divB uses the actual
  evolving face fields; a Yin–Yang face-snapshot/restart comparison is still pending.
- All six full debug dimension/component libraries compiled, including 1D2V.
  The final code adds transfer time to the communication exclusion in the
  compute-only timer; the 500-step run exercises lb_diagnose=T. All 15 final
  conserved-state hashes match the pre-instrumentation results exactly.
- Native spherical CT regression passed 19 runs (18 normal plus 1 expected
  rejection) and 13 binary snapshot comparisons in each of debug/optimized builds.
  The retained directory is
  `/Users/nanami/workspace/2608_uawsom/yinyang-ct-static-native-regression-20261118`.
  This covers native/staged equivalence, 1/2/3 ranks, nontrivial hook scaling and
  TVDLF/HLL without enabling the Yin–Yang experimental mode.
- Final source/test manifests match the workspace; Python syntax and git diff
  checks passed. The final parent commit is c5578cd64bea3d874dc29658e5aec4a0fecedda8
  on codex/yinyang-grid-prep, plus the retained uncommitted source snapshot.
- Including earlier attempts, all builds and regressions, system memory peaked
  at 55.510% of 64 GiB. Build/runtime process-tree RSS peaks were
  1643.156/623.922 MiB. No watchdog stop occurred. Limits stayed at 4096/1024 MiB,
  system 65% start / 68% stop, sampled every 0.2 s; they are sampled protections,
  not hard OS quotas. `final-audit.json` in the final static directory records
  aggregate checks and resource peaks.

The integration now passes actual static CT divergence preservation checks, with
measured rather than suppressed equilibrium drift. Magnetic transport, full 3D
convergence, low-beta conditions, alternate radial boundary models and Yin–Yang
CT restart remain separate validation tasks. No commit or push was made.


## 2026-11-19: actual CT magnetic-loop transport

Added `tests/yinyang/ct_transport` without changing the CT/physics implementation.
A compact vector potential initializes native staggered face fields, then actual
full MHD advances a weak magnetic loop at unit Cartesian speed through t=0.6
(300 SSPRK3 steps, dt=0.002; TVDLF, van Leer, UCT-HLL). The analytic translated
loop is a kinematic reference, not an exact finite-amplitude MHD solution.
The peak field is 1e-3 at rho=1 and p=0.1. Halving that amplitude changes the
final normalized L2 error by only -1.9003e-7 relative, so magnetic feedback is
not the main cause of the measured error. Uniform-flow discretization error
remains present; this is not a low-beta validation.

The loop center crosses the Yin phi=140.625-degree boundary at t=0.152263
forward and t=0.447737 backward. The entire analytic support stays inside the
radial shell 1<r<2. Native continuous CT radial boundaries are used.
An independent finite-difference curl and geometry check verifies the reference.
The control rotates the trajectory to phi=0, away from the angular seam.
Error quadrature uses weight 1/2 in overlap and 1 elsewhere; it does not measure
conservative overset flux exchange. The rotated control also changes grid
alignment, so subtracting its error does not uniquely isolate interpolation.

| Case | Initial relative L2, coarse / fine | Final relative L2, coarse / fine | Final peak / analytic peak, coarse / fine |
|---|---|---|---|
| Seam forward | 0.13095 / 0.07634 | 0.75981 / 0.65102 | 0.27961 / 0.36681 |
| Seam reverse | 0.13351 / 0.07839 | 0.76798 / 0.65692 | 0.26817 / 0.37800 |
| Interior forward | 0.14345 / 0.08377 | 0.77964 / 0.64919 | 0.23345 / 0.36330 |
| Interior reverse | 0.14345 / 0.08377 | 0.77964 / 0.64919 | 0.23345 / 0.36330 |

Coarse/fine per-patch meshes are 16x24x72 and 16x32x96. Radial resolution is
fixed, so this demonstrates an angular-refinement trend, not full 3D order.
The coarse initial pointwise B errors are already 13–14%; native discrete-curl
and center reconstruction errors are included, not subtracted. Final errors
remain 65–78%, with fine peak retention only 36–38%. Thus the implementation
passes stable transport/divergence checks but is not yet qualified for accurate
transport of this narrow loop. The equally diffusive interior controls show
that the large error is not unique to the overset seam. They do not exclude a
smaller seam contribution.

For the fine seam cases, the B-squared centroid moves 0.59162/0.59163 along the
flow, versus the exact 0.600. Final centroid phi=157.177/135.185 degrees confirms
actual traversal of the Yin boundary in both directions. Max active/ghost
native divB over all cases is 4.256e-17, with no discontinuous growth at crossing.
All density/pressure histories remain positive; fine seam final velocity errors
relative to the uniform reference are about 2.5e-3.

Validation and resource audit:

- Eleven production test cases completed: both directions/placements at both
  resolutions, amplitude halving, four ranks, and an optimized representative.
  A separate debug coarse pilot also completed. Each runtime stderr is empty.
- Final cell-centered conserved states AND block-identified active staggered
  face values agree bitwise between 2/4 ranks on the coarse seam case.
  Debug/optimized normalized errors agree within 1e-8. This is not a restart
  or checkpoint face-layout test.
- No common physics/macros were changed. Builds reused the preceding validated
  external debug/optimized libraries, with source-hash checks and invalidation.
  The preceding six-dimension compilation is not claimed as a new test here.
- The initial fine block layout was rejected before launch by a 972.72 MiB
  memory forecast (>950 MiB). Combining the two radial blocks reduces ghost
  allocation while retaining the same physical mesh; the revised estimate is
  810.60 MiB. No protection thresholds were raised.
- Peak runtime RSS was 740.969 MiB; cached test-build RSS was 97.422 MiB;
  guarded plotting/analysis used 86.719 MiB. System memory peaked at 56.620%
  of 64 GiB. Guards remained 4096 MiB build / 1024 MiB runtime, system 65%
  start / 68% stop, sampled every 0.2 s. No running task hit a guard limit.

External evidence: `/private/tmp/amrvac-tests/ct-transport-20261119-suite`,
including `analysis.json`, `results.json`, `transport-summary.png`, per-step
histories, face/state dumps, memory records, and `final-audit.json`.
Pilot evidence: `/private/tmp/amrvac-tests/ct-transport-20261119-pilot`.
The source snapshot is HEAD c5578cd plus the retained prior uncommitted CT work.
No commit/push was made.

Next accuracy work should separate under-resolution and native scheme damping
from seam interpolation using a more broadly resolved smooth magnetic structure
or a linear wave, before drawing convergence-order or low-beta conclusions.
Actual Yin–Yang CT snapshot/restart and alternative radial boundaries remain
separate validation tasks.


## 2026-11-19: wider magnetic structure with a matched-shell control

The follow-up changes only the test definition/tools: support radius is now a
namelist parameter and the runner records radius/shell geometry. No CT/physics
implementation was modified. The wide case has R=0.65 and the code-unit shell
0.75<r<2.25, which contains the analytic support throughout t=0..0.6.
Mesh sizes stay 16x24x72 / 16x32x96 per patch; the radial spacing is 0.09375
on both meshes. The wider diameter spans about 12/16 theta cells and 14 radial
cells. The smaller gradient lengths remain less resolved than those counts.
The vector potential, peak B=1e-3, gas state, flow speed, dt, reconstruction,
flux scheme and CT/RK stage path are unchanged.

To avoid attributing a domain change to improved resolution, the narrow R=0.35
loop was rerun on exactly the same widened shell. `compare_widths.py` verifies
that each paired parameter file differs ONLY in radius. These matched forward
seam results are:

| Mesh | Narrow initial L2 | Wide initial L2 | Narrow final L2 | Wide final L2 |
|---|---|---|---|---|
| Coarse | 0.13189 | 0.04121 | 0.78015 | 0.35501 |
| Angular fine | 0.07763 | 0.02391 | 0.68198 | 0.23689 |

The wide final errors are 45.5% / 34.7% of the matched narrow errors. This is
strong sensitivity to spatial scale, consistent with appreciable diffusion of
under-resolved structure. At fixed travel time and physical displacement, a
wide structure also travels fewer of its own widths; this width comparison is
not a formal convergence experiment or an attribution of all error to resolution.
The two meshes refine only angles, so they do not establish full 3D order.

| Wide case | Final L2 coarse / fine | Fine peak / analytic peak |
|---|---|---|
| Seam forward | 0.35501 / 0.23689 | 0.78763 |
| Seam reverse | 0.35918 / 0.24064 | 0.86733 |
| Interior forward | 0.34067 / 0.21866 | 0.78310 |
| Interior reverse | 0.34067 / 0.21866 | 0.78310 |

The fine seam/interior L2 gap is 1.82 / 2.19 percentage points, compared with
1.43 / 1.85 on the coarse mesh. That gap does not shrink in this two-level
comparison, despite improvement in all total errors. Rotating the interior
control also rotates the field/flow relative to the local discretization, so
this difference cannot be uniquely attributed to overset interpolation.

A second unresolved feature is nonmonotonic magnetic peaks in the seam runs.
After t=0.1, the fine forward/reverse peaks recover by 0.0401 / 0.0805 of the
analytic peak above their preceding running minima, near t=0.404 / 0.522.
The interior curves are smoother. These maxima do not exceed the analytic peak,
but higher peak retention alone is not an accuracy claim: spatial sampling and
local reconstruction must be checked before assigning a cause. The histories
and plot retain this behavior rather than smoothing it out.

Fine forward/reverse B-squared centroid displacement is 0.60292 / 0.60152,
versus the exact 0.600. Max active divB over the wide matrix is 1.0685e-16;
max ghost divB is 8.5654e-17. All pressure/density histories remain positive.
The amplitude-halving coarse check changes normalized final L2 by -1.5742e-6
relative. The reference remains kinematic weak-field translation, not an exact
finite-amplitude MHD solution, and geometric flow errors remain present.

Validation: 11 wide-loop cases plus 2 matched narrow cases completed. Successful
runtime stderr is empty. Coarse 2/4-rank conserved cell states AND active face
dumps match bitwise; the optimized representative passes the existing error
comparison with debug. Analytic curl/normalization/support checks pass for the
wide and narrow geometries. `--complete` and matched-parameter checks passed.
The caches were reused with unchanged common source; six-dimension builds were
not repeated for these test-only changes.

All simulations/builds stayed outside the repository. Runtime process-tree RSS
peaked at 714.078 MiB, cached test compilation at 106.484 MiB, and guarded analysis
at 84.578 MiB. System memory peaked at 56.187% of 64 GiB. No memory guard stopped
a job; limits stayed at 4096/1024 MiB and system 65% start / 68% stop (0.2-s
sampled protections, not hard quotas).

Evidence:
- `/private/tmp/amrvac-tests/ct-wide-20261119`: full wide matrix, plot, histories,
  face/cell hashes, `geometry.json`, `analysis.json`, `width-comparison.json`,
  `final-audit.json`, source manifests and final test-tool snapshot.
- `/private/tmp/amrvac-tests/ct-narrow-matched-20261119`: paired narrow baseline.

The result supports stable cross-seam transport and substantially reduced error
for a broader structure, but the remaining 22–24% fine-grid errors, nonshrinking
seam/control gap and nonmonotonic peaks prevent a high-accuracy qualification.
A same-orientation native spherical control would help separate the native
scheme's direction dependence from seam coupling. Inspecting the two patch
reconstructions at common physical points around the peak recoveries is another
specific next diagnostic. No commit/push was made.


## 2026-11-19: same-orientation native controls and transport-error attribution

This follow-up adds test-only native controls, sparse snapshots and analysis;
common CT/physics code is unchanged. Seventeen debug cases and two coarse
pilots completed outside the repository. The YY baseline cell AND face dumps
match the preceding wide-loop suite bitwise in both directions, confirming that
the added diagnostics did not alter the scientific result.

### Comparison geometry and boundary qualification

Two independently evolved native spherical wedges retain the respective Yin
and Yang local frames, spacings and cell locations. The Yin wedge extends phi
to 46.875..234.375 degrees. The Yang wedge uses theta=14.0625..115.3125 and
phi=-46.875..46.875 degrees, with the same proper rotation as the actual patch.
Both fully contain the translated analytic magnetic support. Cell origins shift
by integer numbers of original patch cells; comparison data require no spatial
interpolation. This is a pair of native wedge controls, not a full-sphere native
simulation.

The first pilot exposed a boundary trap: constant extrapolation of spherical
velocity components does not preserve uniform Cartesian flow at new angular
physical boundaries. Its max velocity error reached 0.430. The corrected
native-only test BC specifies rho=1, p=0.1 and Cartesian velocity in angular
fluid ghosts, leaving magnetic faces to native cont handling and radial BCs
unchanged. Max velocity error then fell to 0.00444. In that Yin coarse pilot,
final magnetic L2 changed only from 0.34124055 to 0.34124067, so the large
boundary fluid error was remote from the loop. This is evidence for that pilot,
not proof of boundary independence for every configuration.

The strict comparison samples both independent native solutions at exactly
the YY cell positions and uses the same overlap partition weights. Initial
normalized YY/native differences are 3.05e-15 / 2.55e-15. Magnetic energy outside
the common cell mask is at most 2.1e-12 of the YY energy at the recorded times.
The reference remains weak-field kinematic translation in full MHD; active
velocities and fields are not reset and no analytic EMFs are injected.

| Fine-grid baseline | Uncoupled native combination L2 | YY L2 | Direct YY/native field difference L2 |
|---|---|---|---|
| Forward | 0.226432 | 0.236885 | 0.046646 |
| Reverse | 0.222526 | 0.240640 | 0.055490 |

Thus most baseline error already exists without overset coupling, but coupling
introduces measurable changes. Scalar norm differences are not an additive
physical error budget; native physical boundaries also differ from overset
boundaries. The combined native diagnostic is a quadrature over two independent
solutions, not a new coupled native evolution.

### Controlled changes to the spatial method

All following rows keep the fine grid, dt=0.002, peak B, flow and end time fixed.
The native column is the Yang-oriented single wedge, not the combined reference
above; the comparisons within each column isolate method-setting changes.

| Method | Native Yang L2 | YY L2 |
|---|---|---|
| TVDLF + van Leer | 0.226786 | 0.236885 |
| TVDLF + MC (woodward) | 0.173705 | 0.187141 |
| HLL + van Leer | 0.184550 | 0.197280 |
| HLL + MC | 0.142847 | 0.158539 |

Time-step controls, dt=0.002 -> 0.001 at fixed end time:

- Native Yang: 0.2267864841 -> 0.2267862783.
- YY: 0.2368853352 -> 0.2368852496.

The YY half-step run doubles the number of overset updates per physical time;
it also changes the result negligibly. Temporal integration is therefore not
the dominant error in this case. Doubling native radial resolution from 16 to
32 at fixed angles lowers L2 from 0.226786 to 0.210702, indicating a smaller but
nonzero radial spatial contribution. These are sensitivity tests, not a full
three-dimensional convergence-order measurement.

For the YY baseline, the best scalar fit to exact B is C=0.816317 and the
remaining normalized shape error is 0.149584. With the volume inner product,

\[
 C=\frac{\langle B_h,B_*\rangle}{\langle B_*,B_*\rangle},\qquad
 \epsilon_{\rm shape}=\frac{\|B_h-CB_*\|_2}{\|B_*\|_2},\qquad
 \epsilon_B^2=(1-C)^2+\epsilon_{\rm shape}^2.
\]

Thus attenuation contributes about 60% of the squared baseline error; shape
and phase/distortion contribute the remainder. This is a mathematical fit,
not a decomposition into independent physical causes. HLL+MC gives C=0.892193
and shape error 0.116243: both improve. The best tested YY L2 remains 0.158539,
so this is an improvement, not high-accuracy qualification.

### What the current code explains

`src/mod_finite_volume.t:244-251` supplies only cmax to CT for TVDLF, but both
cmax/cmin for HLL. `mhd_get_ct_velocity_hll` in `src/mhd/mod_mhd_phys.t:3989`
then takes symmetric positive/negative bounds for the former and distinct bounds
for the latter. The transverse magnetic transport terms used to build edge
EMFs (`mhd_get_ct_edges_hll`, around lines 8740-8758) have the form

\[
 F_B=\frac{a^+v_LB_L+a^-v_RB_R}{a^++a^-}
       -\underbrace{\frac{a^+a^-}{a^++a^-}}_{D}(B_R-B_L).
\]

Here L/R are reconstructed edge-side quantities and a+/a- are nonnegative
signal bounds. Under the symmetric TVDLF choice, D=cmax/2. HLL changes those
bounds and their weights; it also changes fluid FV fluxes, so the HLL test is
not a pure one-coefficient ablation. In the first-order limit, this jump term
corresponds to a diffusivity of order D*dx. With limited higher-order
reconstruction, the jump and error depend on local smoothness, resolution and
limiter: a single constant numerical resistivity should not be assigned to
these runs.

`src/mod_constrained_transport.t:178-236` reconstructs face B and transverse
velocities using the configured limiter. `src/limiter/mod_limiter.t:173-202`
shows the MC and van Leer slope rules. Their finite-resolution smoothing of
curved profiles/extrema is consistent with the measured attenuation and the
benefit of changing reconstruction. Source/flow errors remain nonzero, but the
reference-field-weighted final velocity error is only about 3e-4..6e-4 in the
baseline controls, and the preceding amplitude-halving tests bound magnetic
feedback. These tests do not exhaustively separate every spatial source term.

CT enforces a discrete divergence invariant through the curl update. It does
not enforce the correct amplitude or shape: scaling or diffusively smoothing
a solenoidal field can leave its divergence zero. Small divB therefore does
not establish small transport error.

### Peak recovery: localized coupling effect, not automatically a false peak

Sparse states at t=0, .404, .522 and .6 allow full-vector comparisons at identical
physical points, rather than comparing only moving global maxima.

- Forward t=.404: YY peak lies in Yin cell [8,16,95], the first active layer
  next to the phi interface. |B|/Bmax is 0.90327 versus 0.75050 in the native
  solution at the SAME point. Their largest component difference is tangential
  Btheta, 0.15397 Bmax.
- Reverse t=.522: YY peak lies in Yang cell [8,0,49], the first active layer
  next to the theta interface. The corresponding amplitudes are 0.93383 versus
  0.76626; the largest change is tangential Bphi, -0.20039 Bmax.

Crucially, full-vector errors at these points are BETTER in YY: 0.03630 versus
0.19028 Bmax forward, and 0.05733 versus 0.21028 Bmax reverse. The nonmonotonic
peak recovery must not be called an unphysical amplification solely from its
shape. Coupling can locally recover less-diffused information from the other
orientation. That mechanism is plausible, but the contributions of mapped
EMFs, magnetic ghosts and ghost divergence closure have not been individually
ablated.

Other cells become worse. At final time, normalized squared-error increases /
decreases over the common mask are 0.005222 / 0.000379 forward and
0.009377 / 0.000987 reverse: the net global effect is unfavorable despite local
improvements. The staggered transfer uses trilinear weights and a weighted
Cartesian constant fit (`src/mod_yinyang_staggered.t:176-223`), followed by
receiver-edge quadrature (`src/mod_yinyang_ct_coupling.t:94-176`). Constant-vector
exactness and divB closure do not guarantee accurate transport of a varying
field or commuting interpolation/curl operators. A mismatched mapped EMF
enters the tangential B curl increment divided by the interface-normal spacing;
this provides a mechanism for localized first-layer changes, not a proven
isolation of the responsible subroutine.

HLL+MC reduces global YY L2 by about 33%, but does not eliminate local interface
structure: maximum point error changes only from 0.30792 to 0.29918 Bmax, and
post-t=.1 peak recovery grows from 0.0401 to 0.0662 Bmax. Those peak changes are
not by themselves a regression, for the reason above; decisions must use the
full analytic error, including local and global metrics.

### Conclusion and evidence status

Nonzero truncation error on a finite grid is generic, but the present 24% error
is not an unavoidable floor. Spatial settings demonstrably lower it at the same
resolution. The remaining priorities are better spatial resolution/reconstruction
and an accuracy audit of coupled edge/ghost updates using full-vector errors,
not simply suppressing peaks or further shrinking dt. No production defaults
were changed, and these smooth weak-field results do not qualify shocks or
low-beta corona runs.

All 17 formal cases completed with positive states and empty runtime stderr.
Native divB stays near roundoff. The independent snapshot-volume norms reproduce
the solver histories; the same-grid mapper verifies integer lattice shifts.
Cached debug builds were used; common physics/macros were not changed, so a new
six-dimension compilation was not claimed. The two new YY baseline dumps agree
bitwise with the previous wide suite. Python syntax and whitespace checks pass.

One initial analysis failed because NumPy rejected valid ES25.16 output such as
3.29-100 (Fortran omits E for three-digit exponents). The reader was corrected;
the failed log remains under `analysis-final-memory`, and subsequent complete
analyses passed. This was a reader failure, not a simulation failure.

External evidence is `/private/tmp/amrvac-tests/ct-controls-20261119-suite`:
`controls-analysis.json`, `controls-summary.png`, sampled/final states, parameters,
all memory/exit records, source hashes and `final-audit.json`. Pilots remain at
`ct-controls-20261119-pilot` and `ct-controls-20261119-nativebc`. All simulations
were serial jobs outside the source repository, using 4096 MiB build / 1024 MiB
runtime limits and system 65% start / 68% stop polling. No guard stopped a run.
No commit/push was made.

Final resource audit for this control campaign (including both pilots): runtime
process-tree RSS peaked at 743.563 MiB, cached compilation at 121.750 MiB, and
system memory at 60.280% of 64 GiB. All 17 formal runs and both pilots exited
successfully; the retained analysis-reader failure is identified separately.

## 2026-11-19: analytic edge/ghost attribution and quadrature control

Added `tests/yinyang/ct_accuracy` only; production source remains unchanged
relative to the preceding transport controls. Ten final cases run outside the
repository: debug n=8,12,14,16, quadrature8/16 at n=16, MPI2/4 at n=8, and
optimized n=8,12,14,16. Each patch has n×2n×6n cells. This is instantaneous
transfer/curl attribution, not another evolved HLL+MC transport experiment.

The face reference is the discrete curl of high-order line-integrated A divided
by native face area. It is discretely solenoidal, unlike merely prescribing
analytic point B. Four fields cover constant Cartesian B, affine solenoidal
B=.2(y,z,x), and the forward/reverse wide-loop snapshots at .404/.522. Native
`get_divb` measures angular ghost layers. Active faces are protected bitwise.

### Edge error and circulation compatibility

The ordinary test feeds exact midpoint E×length donors through the actual
same-patch/overset/same-patch mapping. Error norms below are unweighted sample
RMS, normalized by field amplitude; curl errors are rates (inverse code time),
not magnetic errors after a time step. A sample may repeat on shared block
faces, and loop norms include field-free interface samples.

| Field | Edge RMS n8 -> n16 | Interface curl RMS n8 -> n16 | Edge slope | Curl slope |
|---|---|---|---|---|
| Affine | 1.8433e-4 -> 4.8273e-5 | 1.5956e-3 -> 9.0967e-4 | 1.93 | 0.81 |
| Forward loop | .0041105 -> .0016328 | .031842 -> .020798 | 1.33 | .61 |
| Reverse loop | .0046308 -> .0017044 | .033901 -> .020726 | 1.44 | .71 |

The n14 loop samples are sometimes nonmonotonic. These are finite-range
operator-error slopes, not a claim that the evolved scheme converges at .6–.8
order. Native midpoint curl errors fall close to second order in these tests.
Mapping-minus-midpoint curl change is exactly zero outside the interface band.

A concrete representation mismatch is isolated by constant Cartesian E. Native
edges store midpoint E×length, whereas receiver edges use two-point line
quadrature. The receiver integral itself has fourth-order error in this special
case, but combining it with native midpoint edges gives interface curl error
7.2784e-5 -> 3.5538e-5, approximately first order. The test-only oracle adds
`midpoint_exact - integral_exact` on receiver edges. Curl error then becomes
1.1664e-5 -> 2.9345e-6, matching the native midpoint result to <0.005% and
restoring its near-second-order decrease. This is evidence of mixed quadrature
conventions, not a proposal to use analytic fields in production.

For fine forward/reverse loops the same oracle reduces interface curl RMS by
about 11%/13% (.020798 -> .018485 and .020726 -> .018112). It therefore explains
only part of their error. Feeding integrated donor EMFs alone does not fix the
problem: current mapping interprets their line averages as midpoint components.
A compatible reconstruction must account for both donor and receiver integral
representations; simply increasing receiver quadrature order is insufficient.

### Ghost B and divergence closure

| Field | n16 interpolation RMS | n16 interpolation+closure RMS | Closure change |
|---|---|---|---|
| Affine | 1.7714e-4 | 2.3272e-4 | +31% |
| Forward loop | .0012694 | .0012922 | +1.8% |
| Reverse loop | .0014231 | .0018348 | +29% |

Closure preserves the compatible reference to roundoff and constrained ghost
divB is below 3e-15 in amplitude-normalized code units. It is therefore not
corrupting an already-compatible field. It can amplify an interpolated field's
error while removing divergence. Reverse-loop closed RMS decreases from
.0067145 to .0018348 (finite-range slope1.87), with peak error .19485 -> .04990
of Bpeak; no resolution-independent floor is established. Even constant B has
O(h²) ghost error when exact face averages are interpreted as point components.
The existing constant-point tests do not cover this representation difference.

### Validation, resources and next implementation target

- Canonical geometric MPI2/4 records match exactly.
- Debug/optimized differences are at most 6.0e-15.
- Quadrature8/16 maximum difference is 3.84e-10 in absolute field units,
  much smaller than transfer errors; compact-support crossings are not assumed
  to integrate exactly.
- All ten final simulations return zero with empty stderr. Active-face,
  reference-preservation, constrained-divergence and bulk-curl checks pass.
- Final runtime peak RSS 696.813 MiB; analysis 751.406 MiB; sampled system peak
  including analysis 58.90% of 64 GiB. Build limit4 GiB, runtime/analysis1 GiB,
  system65%/68% protections are retained. These are sampled process-tree peaks.
- Earlier helper compile/runtime failures and a conservative forecast rejection
  are retained in external pilot directories. An intermediate n14 measurement
  allowed n16 within the unchanged forecast ceiling; no limits were raised.

Evidence: `/private/tmp/amrvac-tests/ct-accuracy-20261119-final/`, with
`analysis.json`, `final-audit.json`, source provenance, guarded logs and canonical
field records. Test details and limitations are in `ct_accuracy/README.md`.

The first implementation target should be circulation-compatible edge transfer,
including donor/receiver representation, followed by face-average-aware magnetic
reconstruction and reassessment of closure amplification. Accept changes only
if these operator tests improve AND evolved forward/reverse HLL+MC loops improve
without degrading active-face protection, divergence or MPI consistency. No
production reconstruction upgrade, evolved-error improvement or speedup has
been established by this diagnostic campaign. Earlier peak recovery can still
be locally beneficial; suppressing it is not an acceptance criterion.

## 2026-11-19: native midpoint EMF representation implemented

The runtime mesh adapter now constructs one sample at each receiver edge
midpoint and weights its projected Cartesian reconstruction by native edge
length. Both donor and receiver stage values therefore represent
`dt * length * E_tangent(midpoint)`. This removes the previous mismatch between
native midpoint donor/interior values and two-node receiver quadrature. It does
not turn these values into exact line integrals, raise native CT order, or make
overset circulation globally conservative.

Production changes in this step are confined to `mod_yinyang.t` and
`mod_yinyang_ct_coupling.t`. The low-level builder retains historical Gauss
quadrature by default; the mesh adapter explicitly selects midpoint by default.
Its optional `native_midpoint=.false.` retains the legacy diagnostic path.
One/two-node selection is resolved during plan construction; the shared exchange
kernel has no added allocation or string dispatch. Native spherical CT, face
mapping, ghost divergence closure, stage sequencing and energy sync are unchanged.
No commit/push was requested or performed.

### Operator and evolved accuracy

`ct_accuracy` adds mode3 for the actual runtime midpoint transfer alongside the
legacy modes0–2. At n14 (14×28×84 per patch), seam curl RMS changes are:

| Analytic field | Legacy | Midpoint | Reduction |
|---|---|---|---|
| Constant Cartesian B | 4.07391e-5 | 3.82981e-6 | 90.6% |
| Affine solenoidal B | 1.10498e-3 | 1.07236e-3 | 3.0% |
| Forward loop | .0203277 | .0137923 | 32.2% |
| Reverse loop | .0292191 | .0218047 | 25.4% |

These remain amplitude-normalized, sample RMS curl rates. Constant-field
midpoint transfer recovers the native midpoint curl to roundoff: the extra
interface error is removed, while native quadrature truncation remains. The
legacy edge and unchanged ghost metrics exactly reproduce the retained campaign
at matching cases. Three resolutions n8/12/14, MPI2/4 and optimized/debug agree;
quadrature8/16 sensitivity is 4.17e-10 absolute, below the mapping errors. The
extended two-plan diagnostic forecast rejected n16 under the unchanged 950 MiB
ceiling. Rebuilding a single plan repeatedly increased sampled RSS through
allocator retention; the final test retains two plans and uses the safe n14 cap.
Neither failed forecast launched the n16 simulation or raised memory limits.

`ct_transport/run_representation.py` builds isolated old/new executables and
runs identical wide-loop fine meshes (16×32×96 per patch), dt=.002, t=.6, with
both directions and both spatial configurations. The following are final
volume-weighted relative magnetic L2 errors, not the operator norms above:

| Configuration | Direction | Legacy | Midpoint | Reduction |
|---|---|---|---|---|
| TVDLF + van Leer | Forward | .236885335 | .231467245 | 2.29% |
| TVDLF + van Leer | Reverse | .240640492 | .234869429 | 2.40% |
| HLL + MC | Forward | .158538900 | .149844810 | 5.48% |
| HLL + MC | Reverse | .166044455 | .155309843 | 6.46% |

Final maximum point error / Bpeak also decreases: TVDLF forward .307916 ->
.301856 and reverse .314112 -> .306903; HLL+MC forward .299179 -> .207672 and
reverse .240935 -> .214223. These maxima refer to the final snapshot, not the
maximum over time. This is a useful, bounded improvement and supports retaining
the representation change with HLL+MC as the candidate configuration. It does
not establish global second-order transport or solve magnetic ghost accuracy.
No speedup is claimed from these unreplicated validation timings.

### Regression evidence and resources

All 77 normal simulations passed with zero stderr bytes: eight analytic cases,
nine paired transport/MPI runs, 12 mesh-adapter cases, 24 boundary/variable-sync
cases and 24 same-patch ownership cases. Six invalid-configuration tests produced
the expected bounded rejection. Full debug libraries compiled in all six
configurations, including 1D2V. Tests cover uniform/stretched radial meshes,
split/wide block layouts, 2/4/6 ranks for adapter regressions, zero stage factors,
active-face protection, ghost divergence, thermal-energy and centered-B sync.

Paired initial active-cell samples match bitwise. Retained TVDLF forward/reverse
and HLL+MC forward baselines reproduce their old final histories within 1e-12.
New HLL+MC forward cell and face output hashes match bitwise at 2/4 ranks.
Native active divergence stays at roundoff and all evolved states retain
positive density/pressure. `git diff --check` passes; pre-existing source edits
are preserved.

Final-suite sampled peaks: runtime977.813 MiB, build987.266 MiB,
analysis590.391 MiB; system60.70% of64 GiB. The four-rank transport run is included
in the runtime peak. Runtime/analysis1 GiB, build4 GiB and system65% start/68%
stop guards were retained. These are sampled peaks, not hard memory quotas.

Evidence directories share `/private/tmp/amrvac-tests/ct-midpoint-20261119-`:
`analytic-final`, `transport`, `integration`, `boundary`, `patch`.
`baseline/final-audit.json` aggregates checks/resources and `transport/comparison.json`
contains paired accuracy results. The earlier `analytic` and `analytic2` directories
retain the two forecast rejections. The development checkout contains only source,
test definitions and formal documentation; all build/run outputs remain external.

The next accuracy target remains face-average-aware magnetic ghost reconstruction
and its interaction with divergence closure. The current change resolves the
identified EMF representation inconsistency, but is not a replacement for that
separate work or a proof of global circulation compatibility.

## 2026-11-19: magnetic ghost face-average reconstruction investigation

Implemented an **opt-in** face-average reconstruction, without replacing the
validated runtime default. `yy_ct_face_moments=T` in `&methodlist` requires
experimental Yin–Yang CT; otherwise it defaults to false. The option is owned
by `mod_yinyang` and imported by input parsing. The explicit diagnostic API is
`yy_prepare_ct_ghosts(plan,s,face_moments=.true.)`. Midpoint EMF transfer from the
previous step is unchanged.

### Model and representation

The fitted Cartesian magnetic field is

\[
 \mathbf B(\mathbf x)=\mathbf b+G(\mathbf x-\mathbf x_t)/h,
 \qquad \operatorname{tr}G=0,
\]

with eleven independent coefficients. For each of 24 physical donor faces,

\[
 y_f=\frac1{A_f}\int_f\mathbf B\cdot\mathbf n_f\,dS,
 \qquad \min_{\mathbf a}\sum_f w_f(M_f\mathbf a-y_f)^2.
\]

The moments M_f and target-face moment use four-point Gauss integration in both
surface coordinates and native spherical face areas. Two-pass reorthogonalized
QR computes the target's linear combination of donor values during initialization.
Weights retain trilinear locality with a fixed positive 0.01 floor, keeping
rank at exact lattice nodes. There is no field-dependent fitting/solve or new
allocation inside the time-step exchange. Radial receivers now need all three
component families and 24 samples rather than the old eight scalar-Br samples;
packing therefore uses request component row5. Physical-only one-sided radial
stencils cover radial walls. Fit failure is explicit, not a silent fallback.

An initial unweighted fit reproduced constant/affine fields but made compact-loop
errors worse. The weighted version substantially reduced that regression but did
not eliminate it. Positive weights also exposed roundoff-dependent stencil ties:
optimized/debug builds could select opposite neighbors at an exactly aligned
node. Canonical radial/angular node selection fixed that problem; final canonical
records agree to about 1e-15. The rejected unweighted and pre-tie runs remain
external for audit.

### What improved, and what did not

At the coarse analytic mesh, constant-field closed-ghost RMS changes from
3.0575e-4 to 6.23e-16; affine solenoidal-field RMS from 1.1225e-3 to 5.13e-15.
The reference is independently constructed from line-integrated vector potential
and the discrete Stokes curl. The fit therefore fixes the identified mismatch
between face averages and point components for these polynomial fields.

For the nonlinear loops, normalized closed-ghost RMS is mixed:

| Field / n | Previous mapping | Weighted face fit |
|---|---|---|
| Forward / 8 | .00315977 | .00340852 |
| Reverse / 8 | .00671445 | .00707504 |
| Forward / 14 | .000990685 | .000977258 |
| Reverse / 14 | .002277892 | .002439256 |

At n14, reverse-loop interpolation alone is .00156597, increasing to .00243926
through outward divergence closure (about56%). The previous mapping also has
closure amplification (.00151889 -> .00227789). Polynomial reproduction and
small divergence therefore do not establish nonlinear accuracy. Independently
fitting solenoidal polynomials on different faces does not ensure their shared
cell fluxes sum to zero; the closure is unchanged and still required.

Paired evolved HLL+MC tests keep the preceding midpoint EMFs, 16×32×96 cells
per patch, dt=.002 and t=.6:

| Direction | Previous L2 | Face-fit L2 | Relative change |
|---|---|---|---|
| Forward | .1498448099 | .1498576892 | +0.0086% |
| Reverse | .1553098426 | .1553507525 | +0.0263% |

Final max point errors / Bpeak are .20767217 -> .20768115 and
.21422288 -> .21422287. Initial active-cell samples match bitwise; native active
divergence stays at roundoff and density/pressure remain positive. The baseline
reproduces the previous midpoint campaign. These results provide no evidence of
an evolved accuracy improvement, so the new reconstruction remains explicitly
experimental and **off by default**. No fine-grid four-rank evolved comparison
is claimed for this option; component/analytic MPI checks are separate.

The next promising target is a joint face-flux reconstruction with cell-divergence
constraints, rather than independent face fits followed by outward-only normal
flux repair. That is not implemented here. This step provides a consistent
face-average reproduction option, deterministic geometry and independent tests,
not a resolution of the nonlinear ghost/closure accuracy problem.

### Final validation and retained evidence

The final module-owned option passes eight analytic cases and 12 direct ghost
identity/layout cases. Both debug and optimized explicit fills use independent
Stokes references, physical donors and repeated fills with protected faces.
Earlier same-scheme boundary and ownership regressions contribute 24 cases
each, plus four paired transport runs. Two final configuration-switch transport
probes reproduce the corresponding default and experimental cell/face hashes
and complete initial/final diagnostic rows **bitwise**. This is74 successful
normal simulations across these validation suites, all with empty runtime stderr,
plus seven expected rejection cases. All six full debug dimensional libraries
compile with the final option, including1D2V. Edge EMF analytic metrics exactly
match the retained midpoint campaign; common MHD physics and global-parameter
source files are unchanged in this step.

Final-suite sampled peaks: runtime721.766 MiB, build996.734 MiB,
analysis601.750 MiB. Across all exploratory attempts, system use peaked at60.05%
of64 GiB and compilation at1629.906 MiB, below the retained4 GiB build cap.
Runtime/analysis1 GiB, system65% start/68% stop guards remain unchanged.
An intermediate optimized build was deliberately stopped through its monitor's
cleanup handler when the new option was moved from shared globals into its
owning Yin–Yang module; this avoided unnecessary global-interface rebuilds.
Its recorded SIGTERM/monitor interruption is not a memory-limit failure.

Evidence prefix: `/private/tmp/amrvac-tests/ct-face-20261119-`.
`final2` holds final analytical results; `ghost-final` the final six-dimension
builds and ghost regression; `transport` the four evolved comparisons; `switch`
the final default/opt-in reproduction; `boundary` and `patch` the geometry and
variable-sync regressions. `baseline/final-audit.json` aggregates provenance,
checks and resources. Earlier `pilot`, `weighted`, `affine-suite`, `affine-ties`
and `final` attempts are retained. All artifacts remain outside the checkout.
No commit or push was performed, and existing worktree changes were preserved.

## V1 scope freeze and handoff

The agreed endpoint is a reproducible fixed-grid three-dimensional MHD/CT
prototype with stable bidirectional seam transport, controlled native divergence,
correct variable synchronization, MPI consistency and bounded memory. It does
not require zero seam error, strict second-order behavior in every test, AMR,
or qualification for every solar eruption or low-beta configuration.

Retain midpoint EMF transfer and HLL+MC as the recommended tested configuration.
Keep `yy_ct_face_moments` off by default. Further ghost reconstruction research
is not a prerequisite for handing off this V1. Reopen numerical-method work for
a demonstrated stability/divergence/MPI defect, or a specified scientific case
that fails a predeclared accuracy tolerance with evidence implicating the seam.
The current loop errors are documented limitations, not universal accuracy bounds.
