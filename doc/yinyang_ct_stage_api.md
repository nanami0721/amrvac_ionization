# Experimental staged CT API

This is the RK-stage foundation for overset electric-field exchange. It does not
remove the Yin–Yang CT/B0 restrictions or implement face-field interpolation.

## Configuration and lifetime

`ct_stage_split=.true.` in `&methodlist` selects a two-pass update. The default is
false, preserving the original block-at-a-time route. The experimental route
requires 3D, fixed level-one meshes, global unsplit SSPRK3 stepping, TVDLF or HLL
fluxes and UCT-HLL. AMR, load balancing and local timesteps are rejected. Only
these combinations are covered by the stage test; this is not a blanket
qualification of every MHD source/thermodynamic option.

Each active local block owns cached numerical face fluxes and time-integrated
edge electromotive forces. Arrays are allocated outside threaded kernels and
reused across stages. Workspace is not snapshot state and is regenerated after
restart. It adds storage proportional to local block count; the memory watchdog
must cover the complete MPI process tree.

## Sequence

1. Convert the RK input state, reconstruct faces and compute Riemann fluxes and
   UCT-HLL edge line integrals for every local active block. Do not update the
   destination state or invoke geometric/physical source updates in this pass.
2. Call `usr_ct_stage(qt,qdt,sCT,s,fluxes)` once per rank, outside OpenMP loops.
   Different ranks may own different numbers of blocks. This is the location for
   a future world-communicator Yin–Yang exchange, not inside a block hook.
3. Apply the discrete curl to face fields, then the original finite-volume update,
   geometric source, face-to-center synchronization, small-value handling and
   physical sources in their original order. Store final edge integrals for the
   existing same-patch CT correction path.
4. Apply existing flux/edge corrections and fill ghost cells.

`phys_get_ct_edges` computes edge integrals with the existing face-update signature;
`phys_apply_ct_edges(ixI,ixO,fE,s)` applies their discrete curl. UCT-HLL registers
both. The original `phys_update_faces` entry remains as a wrapper calling both,
so original call sites and other CT choices keep their existing behavior.

`ct_stage_flux` holds `face(spatial,nwflux,ndim)` and `edge(spatial,sdim:3)` arrays.
Only active local blocks are populated. Entries outside the generated stencil
are initialized to zero; zero initialization does not make them valid donors.
The exchange callback may change valid edge line integrals. States and face
fluxes are read-only by contract. The integrals already include substage dt and
physical edge length: do not multiply by those factors again. The old per-block
`usr_set_electric_field` hook runs during edge generation, before this callback.

Register `usr_ct_stage` unconditionally in `usr_init`; `&methodlist` has not yet
been read at that point. Check configuration-dependent requirements in
`usr_set_parameters`, after the main parameter reader. The callback is invoked
only when staging is enabled.

The existing compute-only load-balance timer excludes the new callback duration.
No speedup is claimed: staging trades persistent workspace and a second primitive
conversion for an all-block communication boundary.

## Verification definitions

`tests/yinyang/ct_preflight/stage_suite.py EXECUTABLE EXTERNAL_DIRECTORY` runs under
its process-tree watchdog. It checks native versus staged snapshots, including
staggered data, by complete-file SHA-256 equality at 1/2/3 ranks. Twenty blocks
produce an uneven partition at three ranks. Cases include uniform static fields
and nonuniform discrete-curl fields with finite velocity. A half-strength EMF
applied through the stage hook must match the old per-block hook and must differ
from the unmodified solution. Contact CT with staging must be rejected.

The test does not establish Yin–Yang interface consistency, ghost-face divergence
constraints, overlap convergence or global magnetic-flux conservation.

## Validation evidence

The 2026-11-17 temporary test directories are no longer available. Do not use
those paths as retained evidence. Fresh debug/optimized suites and dimension
builds passed and are recorded in the persistent external directory
`/Users/nanami/workspace/2608_uawsom/yinyang-ct-stage-20261118`.

The reproducible entry point is
`tests/yinyang/ct_preflight/build_stage_suite.py EXTERNAL_DIRECTORY`, with optional
`--cache-source EXTERNAL_SOURCE` and `--resume`. Builds use a copied source tree,
serial compilation and per-job process-tree memory guards. Resume requires exact
source/test manifest agreement. `status.json` records completed tasks, while
`manifest.json`, `source.diff` and `commit.txt` identify the source snapshot.
Full accepted results are listed in `yinyang_fixed_grid_progress.md`.

On 2026-11-18 both suites passed 19 runs and 13 binary comparisons each; all six
full debug dimension libraries compiled. Three-rank staged CT restart was binary
identical to uninterrupted evolution. Existing Yin–Yang GLM matched its retained
baseline at 2/4 ranks. See the progress document and external `final-audit.json`
for divergence and memory measurements. Actual Yin–Yang CT remains disabled.
