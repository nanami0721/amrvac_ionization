# Fixed-grid Yin–Yang development progress (2026-11-18)

The full fixed-grid CT/open-corona plan is **not complete**. This increment
hardens validation and checkpoints, adds a common-position overlap diagnostic,
and establishes the ordinary spherical CT prerequisite and a two-pass RK-stage
CT API. Yin–Yang CT, static B0,
open-corona boundaries and Parker/dynamic-event validation remain disabled or
unimplemented. No commit or push was performed.

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
2. Implement M2: stage-coherent donor electric-field reconstruction, edge-integral
   exchange, ghost-face divergence constraints and native face updates.
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
baseline and parts of M1/M4 are available; M2 has its prerequisite test and staged update API, but no overset edge exchange.

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

The next M2 task is stage-coherent electric-field reconstruction and receiving
edge line integration, followed by divergence-compatible ghost-face closure.
The new point-sampling table alone does not justify removing the CT rejection.
