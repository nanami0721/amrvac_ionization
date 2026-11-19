# Magnetic ghost-face fill

`yy_prepare_ct_ghosts` builds a caller-owned `yy_ct_ghost_plan` from the actual
fixed Yin–Yang mesh. `yy_fill_ct_ghosts(plan,s,constrain)` fills angular magnetic
face halos and can immediately apply the previously verified divergence closure.
This is an explicit face-state API; it does not enable full MHD CT or update
cell-centered B/energy. CT/B0 startup rejection remains in place.

## Receiver and donor rules

The receiver list is the union of faces surrounding selected angular ghost
cells: theta slabs over interior phi, then phi slabs including physical
angular corners. Every active face, including active-boundary anchors, is
excluded. Radial boundary faces adjacent to angular ghost cells are included;
physical radial ghost **cells** and unrelated same-patch halos are not filled.
The outer lower face is located from its adjacent center when the native x
array has no index zero. Receiver records are cached and deduplicated per block.

Both donor families are packed before either family is exchanged or any face is
written. Donors are physical interior/boundary faces of the other patch, never
unfilled angular ghost faces. The actual-grid test verifies these donor bounds
and poisons all unfilled halos with a large sentinel to expose accidental reads.
Incomplete coverage is fatal at plan construction. It must be addressed by
sufficient angular overlap, not by extrapolation or a partially filled state.

The two component families deliberately differ:

- Radial Br is a scalar under patch rotation at the same physical position.
  Its own face lattice reaches both radial walls. `radial_scalar=.true.` in
  `yy_staggered_build` supplies an eight-sample scalar stencil, without requesting
  radially centered tangential components outside the domain. On the common
  radial-face grid this reduces to angular bilinear interpolation. It is second
  order; it does not preserve an arbitrary constant Cartesian B to roundoff.
- Tangential receiving components use the existing Cartesian-vector
  reconstruction from three independently staggered component lattices, then
  projection onto the receiving face normal. Their radii are cell-centered and
  have complete donor support, including the first and last physical cells.

The shared radial bracketing now tolerates/clamps only coordinate roundoff
within 1e-12 times the radial scale at its endpoints, matching the pre-existing
angular tolerance. This does not admit a physical ghost-layer extrapolation.
Default three-component table behavior remains unchanged away from that tolerance.

`yy_ct_face_plan/build/exchange` in `mod_yinyang_ct_coupling` provide the underlying
component-family interface. `yy_prepare_ct_ghosts` additionally chooses the
protected write set and supplies native indices and geometry. Plans must be
rebuilt after mesh changes or restart; stale epochs are rejected. Buffers are
reused during fill. No hook is registered automatically.

## Constraints and integration limits

With `constrain=.true.`, the fill is followed by the existing theta/phi outward
normal-flux closure. The interpolated tangential ghost faces and active-boundary
normal faces are its inputs. All faces outside the fill write set remain
bitwise unchanged. The constrained ghost divergence is measured with native
`get_divb`, face areas and volumes.

This local constraint can increase the magnetic interpolation error coefficient;
zero ghost divergence alone is not a magnetic-accuracy test. The suite records
both unconstrained/constrained errors and divergence. It does not establish
strict cross-patch flux conservation or accuracy of a full MHD time integrator.
Same-patch ghost/edge synchronization and mixed-corner ordering are now provided
by the explicit wrappers tested in `../ct_patch_sync`. Physical radial boundary
ordering, radial-wall edge closure, cell-centered B/energy synchronization and
production hook integration still precede evolving Yin–Yang CT activation.

## Reproducible tests

```
python tests/yinyang/ct_ghost_fill/run.py /absolute/external/new-directory --cache-source /absolute/external/source-cache
```

The script copies the current source and its uncommitted changes, invalidates
changed cached templates, and builds serially outside the repository. Source
hashes, commands, parameters and full logs are retained externally. Debug and
optimized runs require empty runtime stderr; full debug dependency libraries
are compiled in all six dimension/component configurations.

The tests use separate probe face states on the actual supported rho Yin–Yang
mesh. They exercise all angular ghost faces (not just isolated sample points),
both radial walls, uniform and geometric radial spacing, both patches and
2/4/6 MPI ranks. Fields are a constant Cartesian vector and the divergence-free
linear Cartesian field B=(y,z,x), sampled as normal point components at faces.
Three repeated fills must be bitwise reproducible; all non-receiver faces must
be exactly preserved. Canonically ordered outputs must match across MPI sizes.

Angular grids are 16x48, 32x96 and 48x144, with eight radial cells held fixed.
The two refinement ratios are 2 and 1.5 and are used explicitly in the measured
orders. The chosen analytic fields remove nonlinear radial-field interpolation
as a confounding source of an angular refinement floor. These tests do not claim
a general three-dimensional asymptotic order for arbitrary magnetic fields or
area-averaged face data. Fixed radial finite-volume truncation is not removed
by angular refinement alone. Radial block size is four, giving the even block count
required by AMRVAC's stretched-grid reader; angular block sizes are 16.

Memory safeguards remain 4 GiB for builds, 1 GiB for MPI, system 65% start and
68% stop with 0.2 s sampling. A one-second test hold exposes the allocated state
to multiple samples. `max_blocks=128` accommodates the maximum 54 owned blocks
rather than carrying the default 4000-block reserve on every rank. The runner
also rejects a refinement if its measured-memory prediction exceeds 950 MiB.
These sampled limits are not a hard OS memory quota.

The new experimental face-average path is requested explicitly in this test.
References now use independent Stokes integration of constant/affine-field
vector potentials, rather than point components. These fields are reproduction
identities, so the runner checks absolute error <1e-10 instead of fitting an
order to roundoff. `--identity-only` runs one mesh with uniform/stretched radii
and 2/4/6 ranks in debug/optimized builds; nonlinear refinement is covered in
`ct_accuracy`. No ghost donors, protected-face changes, nonreproducible repeated
fills or nonempty normal runtime stderr are permitted.
