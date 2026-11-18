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
