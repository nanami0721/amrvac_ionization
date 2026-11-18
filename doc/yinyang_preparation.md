# Yin-Yang coordinate: implementation preparation

Status: design and independent geometry preflight only, 2026-11-15.
No AMRVAC Yin-Yang solver has been implemented, compiled or evolved.
Base commit: 6f5e6f10e5da5719dfc160da4591ae369450b603.
Local branch: codex/yinyang-grid-prep. No commit or push made.

## Decision

Expose `call set_coordinate_system("yinyang")` as a distinct coordinate alongside
`spherical`. Restrict the initial implementation to 3D3V spherical shells with
positive inner radius. Retain existing spherical metric kernels, but make
topology, patch orientation and physical-coordinate conversion explicit.
Do not temporarily change the module-global coordinate inside a kernel: that
would introduce hidden state and OpenMP hazards. Do not duplicate all spherical
operators. Group `case(spherical,yinyang)` only where the local formulas really
are identical; audit every other branch individually. An initialized metric-kind
integer or spherical-metric flag is an alternative for repeated shared dispatch.

## Evidence and references

1. Kageyama & Sato 2004, DOI 10.1029/2004GC000734, Zotero PUAS2AA3,
   attachment BS475Z44, full indexed text (15 pages); PDF page 3 visually checked.
   Equations (2),(4): (x_e,y_e,z_e)=(-x_n,z_n,y_n). Rotation is orthogonal,
   determinant +1 and its own inverse. Section 2 supports reuse of spherical
   metrics. Section 5.2 evolves vector potential A, not an AMRVAC-style face CT
   state. The basic rectangular patches retain approximately 6% overlap even
   as the buffer tends to zero; overlap is not solely a few ghost zones.
2. Luo & Fan 2025, DOI 10.3847/1538-4365/adf8d1, Zotero QLM7XMY5,
   full indexed text read in the preceding review. Section 3.6 describes primitive
   bilinear interpolation and precomputed nonblocking MPI schedules. It does not
   specify an AMRVAC-compatible staggered interface algorithm. Public source was
   not located in the targeted search or HAO highlight; this is not proof that
   no public source exists.
   https://www2.hao.ucar.edu/news/publication-highlight/yin-yang-magnetic-flux-eruption-yin-yang-mfe-code-global-corona
3. Wongwathanarat, Hammer & Mueller 2010, DOI 10.1051/0004-6361/200913435,
   author-hosted preprint, sections 3 and 5 inspected:
   https://wwwmpa.mpa-garching.mpg.de/mpa/publications/preprints/pp2010/MPA2868.pdf
   Relevant warnings: primitive versus conservative interpolation and EOS closure;
   overlap-weighted diagnostics; scalar conservation repairs inconsistent with
   momentum produced negative internal energy in some tests. This is HD evidence,
   not a demonstrated MHD interface solution.
4. Peng, Xiao & Takahashi 2006, DOI 10.1256/qj.05.18: publisher abstract verified,
   full algorithm not yet read. Conservative matching was demonstrated for passive
   transport and shallow water, not for AMRVAC MHD/CT.
   https://rmets.onlinelibrary.wiley.com/doi/abs/10.1256/qj.05.18

## Source audit and required changes

| Area | Owning files / routines | Preparation required |
| --- | --- | --- |
| Configuration | src/mod_geometry.t:set_coordinate_system; src/io/mod_input_output.t:read_par_files | New enum, 3D3V validation, angular units (currently input turns are multiplied by 2*pi), unsplit stepping, overlap/bounds validation |
| Local geometry | src/mod_geometry.t; src/amr/mod_amr_solution_node.t | Reuse spherical face areas, volumes, ds, gradient/curl and ghost geometry; verify complete cases rather than renaming |
| MHD sources | src/mhd/mod_mhd_phys.t:mhd_add_source_geom and semirelativistic/split variants | Shared spherical formulas; unsupported combinations must fail explicitly |
| Patch topology | src/amr/mod_forest.t; mod_amr_neighbors.t; src/mod_connectivity.t | No polar wrap and no angular periodicity for Yin-Yang; separate overset boundary classification and donor mapping |
| Communication | src/mod_comm_lib.t:comm_start; src/mod_ghostcells_update.t:getbc | Explicit local/world ranks, schedule, stage-state input; do not put collectives inside per-block/OpenMP boundary callbacks |
| Global stepping | src/mod_dt.t:setdt; src/mod_advance.t:advect1; src/mod_supertimestepping.t | World-wide dt, cmax for GLM, stop/crash synchronization, RK stage consistency; STS cycles deferred initially |
| Magnetic interfaces | src/mhd/mod_mhd_phys.t:mhd_update_faces_*; src/mod_fix_conserve.t:fix_edges | Define evolved and interpolated faces, edge EMF and divergence budget; existing aligned-AMR correction is not an overset algorithm |
| Coordinates and output | src/io/mod_calculate_xw.t:to_cartesian; src/io/mod_input_output.t | Patch-aware global positions/vectors, disjoint filenames and restart metadata; same local coordinates in both patches represent different physical points |
| Additional physics | src/physics/mod_b0.t, mod_rotating_frame.t, thermal conduction, tracing, topology, emission; particle modules | Global-axis and patch-aware operations; reject unsupported options in the prototype |

## Prototype architecture

Recommended initial tradeoff: one MPI rank group per patch, each with its own
existing tree, using a patch communicator for ordinary block communication.
Keep a world communicator for overset exchange and global decisions. This avoids
an immediate redesign of all root-tree and Morton structures, at the cost of
fixed patch assignment, at least two ranks, and no cross-patch load migration.
It is an experimental architecture, not a final AMR design commitment.

The current startup order is comm_start -> usr_init -> initialize_amrvac/read
parameters. Since usr_init selects coordinates, communicator setup needs a
deliberate insertion point after coordinate selection but before tree allocation.
Any rank-sensitive usr_init behavior must be reviewed. Keep world identity
separate from patch-local mype/npe. Both patch roots otherwise become rank zero,
which affects logging, files and restart. Audit all MPI_COMM_WORLD uses, multigrid
defaults, and MPI_ABORT semantics; assigning icomm alone is insufficient.

For an eventual unified AMR forest, add patch IDs to tree/state ownership,
root indexing, SFC metadata, migration, donor search and I/O. Do not encode a
second patch as an extra spatial dimension or artificially extended longitude.

## Numerical contracts to settle before solver implementation

- Freeze one rotation convention and local-to-global APIs for positions and
  vectors. Use atan2-based angle recovery and an explicit phi interval.
- Freeze whether ghost exchange interpolates primitive or conservative states.
  For the initial ideal-EOS prototype, compare positivity-safe primitive
  reconstruction followed by a single EOS/conservative conversion with a
  conservative alternative. Do not independently interpolate redundant p/e/T.
  Neither point interpolation guarantees global finite-volume conservation.
- Reference vector transfer: rotate each collocated donor vector to global
  Cartesian components, interpolate, rotate into the receptor basis. Interpolating
  donor spherical components then transforming at the receptor is another
  approximation, not numerically identical; compare consistency/order explicitly.
- Reuse radial indexing only when the donor and receptor radial locations match.
  This holds for centre ghosts on identical radial grids, not automatically for
  staggered B components or differently refined AMR blocks.
- Donors must be available at the same stage time. Never use the permanent ps
  state when getbc is operating on another psb RK stage. Fill ordinary donor
  block ghosts as needed; avoid cyclic overset-ghost dependencies. Specify corner
  ordering where radial physical boundaries meet angular overset boundaries.
- Store receiver location, donor patch/rank/block, stencil, weights, vector basis
  transforms, counts/displacements and cache epoch. Use buffers allocated once
  for a fixed mesh, grouped by peer. Rebuild after topology/ownership changes.
- CT needs separate maps for face/edge locations and a decision on which
  interface faces are evolved rather than overwritten. Check divergence in
  interior cells, interface-adjacent cells and the physical composite domain.
- Count each physical volume once using overlap fractions or a defined composite
  ownership integration. Validate shell volume and analytic integrals before
  interpreting conservation residuals. Weighting is not a conservation repair.
- Couple any conservative interface treatment for mass, momentum, energy and
  magnetic fields consistently; include EOS positivity. GLM alone solves neither
  overlap conservation nor interface continuity.

## Initial supported scope and explicit exclusions

First scalar/HD stage: fixed blocks, identical radial grid, ideal EOS, same time
step and integrator, no refinement or dynamic load migration. Then cell-centred
ideal MHD with GLM for coupling tests. CT is a separate milestone, not enabled
by default before face/edge acceptance. Initially reject B0 splitting, particles,
self-gravity/multigrid, rotation, TRAC, radiative transfer, conduction/STS and
other unreviewed couplings. These are staged exclusions, not inherent limitations
of Yin-Yang geometry. Support radial analytic physical boundaries first.

## Acceptance ladder

1. Geometry: point/vector round trips, proper rotation, full sphere coverage,
   all centre/face/edge receiver stencils, angular seam and corner cases,
   positive metric quantities; validate volume integration with overlap removed.
2. Single spherical-patch baseline on the chosen commit, both optimized/debug.
   Identical local metric/source results for spherical versus Yin-Yang patch mode.
3. Scalar solid-body rotation on a shell around an oblique global axis, multiple
   crossings and return periods, at least three resolutions; report norms,
   extrema and mass error. Rotate the patch orientation to expose seam artifacts.
4. HD: uniform global Cartesian flow, oblique shock and off-centre blast, positive
   density/pressure, integral balances including physical boundary fluxes.
5. MHD: globally uniform B, oblique waves and field-loop crossing; quantify
   interface/interior divergence, spurious currents, pressure and reflection.
6. CT acceptance separately: initial discrete divergence and its evolution,
   face consistency and EMF handling, global boundary magnetic flux.
7. Parallel/restart: 2,4,8 ranks (when available), different partitions, restarted
   versus uninterrupted runs, collective ordering and bounded debug stderr.
8. AMR and realistic corona only after these pass. Include refinement, coarsening,
   migration invalidation and mismatched-level donors. Add thermal conduction
   substeps before using coronal energy budgets.

Choose physical tolerances and convergence requirements before runs; do not
retrofit them after observing an error. Shared MHD macro changes additionally
require full module/dependency compilation in all six project configurations,
including 1D2V, even though Yin-Yang itself is restricted to 3D3V.

## Preparation checks performed

`tests/yinyang/geometry_preflight.py` is a standalone stdlib geometry survey,
not an AMRVAC implementation. It checks angular cell-centre ghost stencils at
32x96 and 64x192 interior resolution, 2/3/4 ghost layers and buffer widths
0.51/1/2/3/4 times the nominal theta step: 30 configurations. All surveyed points
have a four-point interior donor stencil; maximum position round-trip error is
3.3306690738754696e-16. Radial ghosts, face/edge maps, interpolated-field accuracy,
AMR, MPI and time evolution have not been tested. This finite survey is not a
proof for arbitrary resolutions or buffers.

Local mpif90 resolves to Homebrew GNU Fortran 14.2.0. Compiler presence is checked;
no baseline rebuild or MPI runtime test has been performed. Project debug.defs
uses -fcheck=all and invalid/zero floating-point traps, among other flags.

## Workspace and next actionable preparation

Source worktree: /private/tmp/amrvac-tests/yinyang-prep-20261115/source.
Evidence and results are sibling directories, outside the source worktree.
The original worktree remains on codex/radsyn-streaming. Its existing unrelated
uncommitted radiation/I-O edits were not imported; a read-only diff/status snapshot
is stored in evidence. This branch therefore assesses the clean recorded commit,
not runtime behavior of those concurrent changes.

Before building, make a separate build-source copy containing this branch's
intended edits and point AMRVAC_DIR at that copy, not at the main repository or
the source worktree. Copy cases into a separate run directory. Record commit,
patch, compiler, flags, rank count and run parameters. Keep logs bounded and
outside source; do not commit/push unless requested.

Next: settle numerical/communicator contracts above, expand the geometry probe
to face/edge and volume integration, and establish a reproducible spherical
baseline. Obtain full Peng et al. paper for conservative-interface design and
seek MFE source/interface documentation if publicly available. These references
do not block the scalar/HD geometry and communication prototype.

## Implementation follow-up (2026-11-16)

The experimental implementation and actual acceptance evidence now live in
`doc/yinyang_v1.md`. Earlier sections above record research/preparation and are
not statements that conservative reconstruction or the deferred physics has been
implemented. The v1 implementation uses primitive interpolation and GLM.
