# Face/edge donor schedules and analytic verification

This prerequisite does **not** enable Yin–Yang CT. It builds reusable schedules
for scalar normal face components or tangent edge components, reconstructing a
Cartesian vector at specified receiver points. The existing GLM schedule is
unchanged; no state, active face or ghost face is overwritten.

## Geometry and API

`mod_yinyang_staggered` has no solver-state dependency. Call
`yy_staggered_build` collectively with the complete fixed-grid block directory,
shared radial face coordinates, uniform angular extent, receiver patch and
receiver `(r,theta,phi)` coordinates. `yy_build_staggered_table` in `mod_yinyang`
provides the AMRVAC adapter: it gathers actual block ownership and obtains radial
faces from `ps%x + ps%dx/2`, checking consistency across patches and blocks.

Local cell indices in `request(2:4,:)` are zero-based. Add AMRVAC's `ixMlo`
indices to address actual arrays. Face component c has offset 1 along c and 1/2
along the other axes; edge component c has offset 1/2 along c and 1 along the
other axes, measured from the lower cell corner. Thus array index i represents
the upper face/edge in directions with offset 1. A lower-domain boundary sample
may have index -1. Shared block interfaces use one canonical donor owner.

Each receiver requests eight samples per component (24 scalar samples total).
The component stencils are independently bracketed, including radial offsets;
radial faces may be nonuniform. Full stencils must lie inside the other patch's
physical sampling lattice, including its boundary faces/edges. There is no
extrapolation, clamping to a different physical target or dependency on angular
ghost donors. `covered=.false.` marks receivers lacking any component stencil;
returned zeros for these receivers are placeholders and must be masked. In
particular, a full vector at a radial wall generally needs additional boundary
closure for components centered radially. This module does not supply it.

The owner packs `samples` in `request` order from its scalar face/edge arrays.
`yy_staggered_exchange` executes one world-communicator exchange and returns
Cartesian vectors. Receiver counts may differ or be zero; all ranks must call.
The tables retain request lists, peer counts/displacements, reconstruction
coefficients and communication buffers. Repeated calls do not allocate storage.
`yy_staggered_clear` releases them. The caller must rebuild after restart or any
mesh/ownership change; these tables are not checkpoint state. Requests are not
yet deduplicated; do not infer production communication efficiency from these
small analytic tests.

## Reconstruction used for validation

For donor sample j, let b_j be its scalar spherical component, e_j its unit basis
vector in Cartesian coordinates and w_j its component's trilinear weight. Use

    M = sum_j w_j e_j e_j^T
    v = M^{-1} sum_j w_j e_j b_j.

The inverse and coefficients are cached at construction. This is a weighted
constant Cartesian fit, not independent interpolation of three collocated
components. It reproduces a uniform Cartesian vector to roundoff despite the
staggered positions. Smooth-field accuracy is measured below. A singular or
ill-conditioned normal matrix is rejected. No constrained divergence or flux
conservation property follows from this fit.

Inputs are point components: **not face fluxes/area integrals, nor CT edge
integrals containing dt and edge length**. The later CT implementation must
account for this distinction and evaluate/integrate reconstructed electric
fields along receiving edges. This test does not cover that integration.

## Running (outside the source repository)

Use a Python environment with NumPy (required by the shared memory guard).

```
python tests/yinyang/staggered_schedule/run.py /absolute/external/fresh-analytic-directory
python tests/yinyang/staggered_schedule/integration.py /absolute/external/fresh-mesh-directory --cache-source /absolute/external/previous-source
```

Both runners preserve source copies/hashes, commands, stdout/stderr and watchdog
records. They use serial builds, 4 GiB build and 1 GiB runtime limits, 65% system
start and 68% stop thresholds. The guard samples at 0.2 s; it is not a hard quota.
The analytic runner expands the original VACPP template before compiling it.
The integration runner invalidates changed cached generated sources/objects,
builds the real AMRVAC executable and all six debug dependency libraries.

Analytic cases cover both patches, all six face/edge orientations, angular
boundary rings/corners, three shell radii, uniform/geometric radial spacing,
three mesh resolutions and 2/4/6 ranks. Uniform B=(0.2,-0.3,0.5) checks rotation
and normalization. The smooth divergence-free field (sin y,sin z,sin x) checks
accuracy (not discrete divB). Repeated exchange, explicit uncovered radial/angular
points, empty-receiver ranks, cache destruction/rebuild and MPI-independent
canonical output hashes are checked. Each successful debug run must have empty
stderr, including runtime temporary-array warnings.

The independent driver decodes donor block IDs and indices to construct scalar
samples; the integration test additionally compares cached geometry to actual
AMRVAC arrays and existing `yy_basis`/`yy_position`. The test only samples the
staggered geometry: the solver still runs the supported rho/GLM geometry path.
