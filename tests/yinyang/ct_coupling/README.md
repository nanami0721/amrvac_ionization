# Edge-integral exchange and angular ghost-face constraints

These are explicit CT coupling APIs, not automatic Yin–Yang CT enablement.
The existing CT/B0 rejection remains. They supply the edge exchange and ghost
constraint operations needed before the full magnetic ghost-fill/update ordering
can be enabled and tested in a time-evolving MHD problem.

## Edge integrals

`mod_yinyang_ct_coupling` adds `yy_ct_edge_plan`, `yy_ct_edge_build`,
`yy_ct_edge_exchange` and `yy_ct_edge_clear`. A plan accepts oriented receiver
edge centers in local spherical coordinates, direction 1/2/3 and coordinate
width dr/dtheta/dphi. Both two-point Gauss nodes must have complete staggered
component donor stencils. Generic plans expose `covered`; uncovered outputs are
zero placeholders, not valid interpolated integrals.

The native AMRVAC quantity is

    F_d = dt * L_d * E_d.

The exchange divides by donor length L_d only, reconstructs the Cartesian vector
`dt*E` at each Gauss node using the preceding staggered-table method, and returns

    F_r = sum(q=1,2) dt*E(x_q) dot e_c(x_q) * J_c(x_q) * Delta_x_c/2,

where `J=(1,r,r*sin(theta))`. No second multiplication or division by dt occurs;
a zero-time increment is well defined. Direction follows the positive local
coordinate axis. Radial edges are straight; theta/phi edges follow their actual
spherical coordinate curves. Vector projection is evaluated independently at
both nodes. Shared nonuniform radial faces determine donor lengths.

Donor inputs are native midpoint approximations to integrals. Treating an exact
line-average input as a point value after division by length adds a second-order
error. The tests cover both input conventions; no fourth-order overall transfer
claim follows from the two-point quadrature.

`yy_prepare_ct_edges` in `mod_yinyang` adapts the actual AMRVAC mesh to a
`yy_ct_stage_plan`. Receiver records are `(local_block_id,i,j,k,direction)` with
native array indices. It validates both receiving and donor indices against the
computed UCT edge stencil. It rejects incomplete coverage (including unresolved
radial-wall stencils), rather than silently mixing exchanged and native edges.

`yy_exchange_ct_edges(plan,fluxes)` reads/writes the actual `ct_stage_flux%edge`
arrays. All donor values are packed before any receiver is overwritten. Invoke
this **once per rank at the all-block RK-stage callback**, not inside a block
loop. Ranks must exchange the same substage; the API does not infer stage time
from data. Receiver records must select the intended interface edges. Changes
to face fluxes and states are outside this API. Repeated exchange does not
allocate buffers. Plans are caller-owned, are not checkpoint data and must be
rebuilt after restart or mesh changes; stale mesh epochs are rejected.

This increment does not install the callback automatically, choose a production
edge-ownership policy, synchronize duplicate same-patch boundary edges after
transfer, or complete electric/magnetic conditions at radial-wall corners.

## Ghost face constraints

`yy_ct_close_ghosts` takes face B values, positive native face areas, array/cell
bounds and the four angular patch-side flags. Its face convention is the native
upper-face index: `face(i,j,k,c)` bounds cell `(i,j,k)` on its positive c side.
Tangential ghost faces must already be populated consistently. The active-domain
boundary normal face is kept as the anchor. For a lower theta ghost cell:

    A_theta^- B_theta^- = A_theta^+ B_theta^+
                           + (A_r^+ B_r^+ - A_r^- B_r^-)
                           + (A_phi^+ B_phi^+ - A_phi^- B_phi^-).

For an upper ghost cell the outward upper normal flux is solved instead. The
routine sweeps from the active boundary outward. Theta slabs with interior phi
are handled first; phi slabs then include the theta/phi corners. Each ghost cell
therefore gets one outward normal face correction. Tangential faces are fixed
inputs during each sweep. The first array index is the innermost long loop.

Only outward angular ghost faces are modified. Active faces, active-boundary
faces, radial physical ghosts and unrelated same-patch ghosts are preserved.
All needed input faces/areas are validated in a dry pass before the first write;
nonfinite input, invalid bounds or nonpositive area returns nonzero `status`
without a partial update. Unused native metric-array corners are not inspected.
Repeated application is idempotent in the tested layouts.

The constraint enforces the finite-volume flux-divergence numerator locally.
It does not make the cross-patch interpolation globally conservative, fix
pre-existing active-cell divergence, or by itself establish magnetic ghost-field
accuracy. Same-patch ghost synchronization, radial boundaries, face-to-center
conversion and energy handling must still be ordered correctly by the caller.

## Validation and reproduction

Run in new external directories, using the established scientific Python:

```
python tests/yinyang/ct_coupling/run.py /absolute/external/new-analytic-directory
python tests/yinyang/ct_coupling/integration.py /absolute/external/new-mesh-directory --cache-source /absolute/external/source-cache
```

The runners preserve copied source, hashes, commands and per-job watchdog
records. Builds are serial and limited to 4 GiB, MPI jobs to 1 GiB; system start
and stop thresholds remain 65% and 68%. Sampling is not a hard memory quota.

The independent suite compares all edge directions and both patches to a
128-panel composite-Simpson line-integral reference. It tests midpoint and
line-average donors, uniform/geometric radial meshes, three resolutions,
2/4/6 MPI ranks, repeated exchange, zero dt, linear dt scaling, incomplete
coverage and empty receivers. Canonical result hashes must match across MPI
sizes. Smooth-field errors are divided by dt and receiver arc length so the
observed order does not gain a spurious factor of grid spacing.

Ghost tests cover all 16 combinations of angular-side flags on spherical areas,
including corners, with arbitrary nonzero input flux imbalance. They check
corrected divergence, exact preservation of every face outside the permitted
write set, idempotence and failure-before-write for invalid metrics.

The AMRVAC integration suite packs native `dt*dsC*E` into actual stage-workspace
types and verifies the receiving integrals. It evaluates ghost constraints with
native `surfaceC` and `get_divb` on a separate probe state. The supported rho
Yin–Yang geometry path supplies the mesh; no production MHD state is changed and
no MHD CT evolution is claimed. Debug/optimized builds cover 2/4/6 ranks and
uniform/stretched radial meshes, and all six full debug dependency libraries are
compiled. Successful runtime stderr must be empty on repeated calls.
