# CT radial boundaries and cell-variable synchronization

These explicit interfaces extend the fixed-grid Yin–Yang CT components. They do
not register production hooks or remove the Yin–Yang CT startup rejection.

## Radial edges

`yy_prepare_ct_edges(...,radial_extrapolate=.true.)` opts into a one-sided linear
radial closure for the radially centered E_r samples. At either end it uses the
first/last two **physical** radial edge samples. Tangential edge components retain
their existing radial-face stencils. The table admits points only inside the
physical radial shell; angular coverage checks remain unchanged. Default behavior
still rejects incomplete radial stencils. Negative linear extrapolation weights
are intentional; this is not a positivity-preserving or characteristic boundary.
Legacy two-node quadrature can therefore sample first/last radial edges, while
the current native midpoint adapter samples their centers. Tangential edges on radial walls
can also be exchanged without borrowing unfilled radial ghost fields.
Inputs remain dt*length*E, with no division by dt. This supplies donor coverage;
it does not prescribe a conducting, driven or line-tied wall electric field.

## Magnetic radial halos

After cell thermodynamic ghost data have been supplied by ordinary/overset
exchange, call `yy_complete_ct_faces(patch_plan,ghost_plan,s)`. It first invokes
the previously tested ordered angular face fill, then `yy_ct_radial_faces`:

- Extend tangential face B constantly in radius from the adjacent physical cell.
- Preserve the normal face on each radial wall.
- Sweep outward through radial ghost cells and solve the normal face flux from
  zero discrete divergence, using native face areas. The full angular halo is
  included, so radial/angular mixed corners are closed last.

The six-face flux sum is zero to roundoff in those ghosts. Active-cell divergence
is untouched. This is a specified zero-divergence ghost policy, not the native
`bc_cont` policy that propagates a neighbouring active divergence. It is not a
characteristic solar-wind outflow. Tangential constant extrapolation is a boundary
model choice; no second-order general physical boundary accuracy is claimed.

## Center B and energy

Call `yy_sync_ct_state(s,mhd_ct_sync_center)` after the face fill. The callback
uses the existing native face-area-weighted `mhd_face_to_center` operator.

1. Active cells: replace center B, retain total energy exactly, matching the
   existing conservative CT convention.
2. Angular and same-patch ghost cells: replace B and add
   `Delta E = sum((Bnew-Bold)*(Bnew+Bold))/2`, retaining thermal energy.
3. Physical radial ghosts: first copy conserved variables from the now-consistent
   boundary cell, including completed angular halos, then replace B and compensate
   magnetic energy. This ordering avoids propagating the pre-sync active pressure.

Ghost slabs partition the allocation without duplicate visits. Density, momentum
and other variables are preserved except for the explicitly specified constant
radial copying. The helper does not floor pressure, repair pre-existing negative
pressure or change momentum. Preserving ghost pressure is not a conservative
energy correction to active cells and is not a low-beta solver cure. The helper
rejects background splitting, semirelativistic MHD and non-total-energy modes.
It is an explicit API; existing face-to-center call sites are unchanged.

## Tests and limits

`run.py EXTERNAL_DIRECTORY --cache-source EXTERNAL_SOURCE` copies all source and
uncommitted changes into the external directory. Serial compilation uses a 4 GiB
process-tree guard; MPI uses 1 GiB, with system 65% start / 68% stop thresholds,
0.2 s polling, max_blocks=128 and a one-second allocated-state hold. Full logs,
source hashes, provenance, parameters and guard reports stay outside the repo.

The suite uses supported GLM MHD geometry with independent probe staggered states,
not a CT-enabled time evolution. It includes radial-wall/first/last-cell edge
receivers (constant and radially linear Cartesian E), uniform and stretched radial grids, two block layouts, 2/4/6 MPI ranks,
debug/optimized builds, and full debug libraries for all six dimension/component
configurations. Constant and linear Cartesian B test the full face boundary fill.
Checks include native ghost divB at all physical corners, protected active faces,
unchanged active energy, preserved ghost thermal energy (2e-13 absolute tolerance),
native center-B agreement within 8 machine epsilons, and repeated-fill bitwise idempotence. Canonical B/E
records must match across MPI sizes; successful runtime stderr must be empty.

Before production activation the complete stage/ghost callbacks still need wiring,
then actual CT time-evolution, restart and seam-crossing convergence tests. Specific
solar boundary models (line tying/driving/outflow) require separate physical
choices and tests beyond this explicit continuous radial ghost policy.
