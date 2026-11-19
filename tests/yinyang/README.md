# Yin–Yang tests

New Mac Studio / new developer: start with the repository-root
[handoff guide](../../YINYANG_MAC_STUDIO_HANDOFF.md). It distinguishes current
capabilities, pending CT restart/application work, and a cache-free bootstrap.

See `doc/yinyang_v1.md` for scope, geometry, restart and validation details.
Never compile or run inside this source tree. `run_validation.py` copies the
source and MHD case to a new external directory, checks MPI/restart equivalence,
and optionally runs convergence/perturbation tests and all six full-library builds.

`geometry_preflight.py` is an independent geometry-only probe; it does not validate
VACPP expansion, the Fortran exchange, or MHD evolution.

`mhd/` supplies uniform, entropy-pulse and blast initial conditions. `hd/` supplies
uniform-flow/perturbation initial conditions. `rho/` is stationary scalar
preservation. Diagnostics in `check_mhd.py` require NumPy and native snapshots.

For communication-only changes, `compare_exchange.py BASELINE_RUN CANDIDATE_RUN`
compares all stored fields exactly and checks donor reduction diagnostics and
runtime stderr. See `doc/yinyang_exchange_validation.md` for the deduplication
implementation and measured costs.

`interface_transport/` adds exact forward/reverse entropy and force-free magnetic
transport, three-resolution seam convergence, interior controls and seam-inclusive
divergence diagnostics. Results: `doc/yinyang_interface_validation.md`.
