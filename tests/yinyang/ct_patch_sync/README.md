# Same-patch CT synchronization and ordered ghost fill

This explicit fixed-level API joins the existing edge/face transfer components.
It does not register production hooks or enable Yin–Yang MHD CT.

`yy_prepare_ct_patch` builds three reusable copy schedules from the actual block
directory. Edge copies use a geometrical owner: the lower block at a shared
transverse boundary, independently in each transverse direction. Four-block
edge corners therefore have one owner independent of MPI rank. Values are copied,
not averaged. Every donor is packed before any receiver is overwritten, including
same-rank peers. Plans must be rebuilt after layout changes/restart.

`yy_ordered_ct_edges` performs same-patch edge synchronization, cross-patch edge
integration, then same-patch synchronization again. The caller must select only
geometrical owners as cross-patch receivers (use `yy_ct_edge_is_owner` on valid
computed edge indices). A nonowner receiver is rejected before communication;
otherwise the final owner broadcast could silently erase the overset update.
The wrapper modifies edge integrals only, which already include dt and length.

`yy_ordered_ct_ghosts` performs three steps:

1. Populate same-patch magnetic face halos from physical-domain faces.
2. Populate cross-patch angular halos and apply outward normal-flux closure.
3. Copy the resulting angular halos into mixed block/physical-boundary corners.

Both face copy schedules protect all active faces, including boundary anchors.
The last phase permits donors in angular physical halos but excludes radial
physical ghost cells. It fills mixed corners from a boundary block that has
already completed its cross-patch closure. Physical radial ghost cells remain
caller-owned. These routines operate on explicit face states; cell-centered B,
energy consistency and the production physical-boundary sequence remain pending.
The cache has persistent packing buffers and performs no per-exchange allocation.
It currently uses all-to-all communication, including local peers; no performance
scaling claim or conservative overset flux claim is made.

## Validation

Run outside the source repository:

```
python tests/yinyang/ct_patch_sync/run.py /absolute/new/external/directory --cache-source /absolute/external/source-cache
```

The runner preserves source hashes, parent commit, tracked diff, parameters and
bounded-console/full-file logs. It invalidates changed template objects in an
external cache copy. Builds run serially with a 4 GiB process-tree limit; MPI jobs
use 1 GiB. The system guard uses 65% start / 68% stop thresholds and 0.2 s sampling.
The test holds allocated state for one second and sets max_blocks=128.

The mesh is 8x16x48 per patch. Block layouts 4x8x8 and 4x16x16 exercise both
four-block shared edges and mixed radial/angular corners. Debug and optimized
builds run with uniform/stretched radii at 2/4/6 ranks (24 successful cases).
Canonical face records must match bitwise across MPI decompositions within each
layout. All six full debug dimension/component libraries are compiled.

The probe tests deliberately inconsistent duplicated EMFs, untouched invalid
edge slots, idempotence, constant Cartesian E integration after corrupting
nonowners, and bitwise equality of all shared copies after overset exchange.
Face tests poison halos, require mixed-corner holes before phase 3, then require
complete population and native get_divb <= 2e-11 on external angular ghost cells,
including mixed corners. Constant and linear Cartesian B=(y,z,x) are used.
Active faces and radial physical ghost cells must remain bitwise unchanged.
Repeated ordered fills must match the manually sequenced reference exactly.
Runtime stderr must be empty in successful cases. Two additional negative cases
must reject a nonowner overset receiver without timing out.

Edge interpolation tests exclude unresolved radial-wall edge stencils. This
suite verifies transfer ordering on supported rho geometry with probe face states,
not a time-dependent MHD run or full CT divergence preservation during evolution.
