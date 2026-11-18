# Native spherical CT prerequisite

This is a 3D3V ordinary spherical patch, not an enabled Yin-Yang CT calculation.
It initializes a uniform Cartesian field through the native discrete curl of
A = B cross r / 2, evolves five three-stage steps with UCT-HLL and TVDLF, and
checks native get_divb over all active cells before/after each full step.

The dimensional test normalization is fixed: shell thickness 1 and B magnitude
sqrt(0.14). The raw-divergence assertion is 1e-10 in those units. It is not an
arbitrary-amplitude normalized divergence test or a long-time validation.
Run both debug/default and 1/2 ranks using an external source/build copy and the
process-tree memory guard. Logs include bounded CT_PREFLIGHT records; runtime
stderr must be empty. All runs in the initial implementation produced 11 records.

This is a prerequisite for later electric-field exchange. It provides no evidence
for overset face-field exchange, patch-interface convergence or global magnetic
flux conservation. Those must pass separate tests before enabling Yin-Yang CT.

## Staged advance tests

After building this test against the current source, run
`python stage_suite.py /absolute/path/to/amrvac /external/new/output-directory`
with the scientific Python environment. Run once per debug/optimized executable.
The suite uses the memory watchdog, includes unequal block counts on three MPI
ranks, and checks whole snapshot files including staggered fields. It also tests
HLL fluxes and stage-level versus block-level EMF scaling. These are ordering and
regression tests; the scaled-EMF experiment is not a physical coronal model.
