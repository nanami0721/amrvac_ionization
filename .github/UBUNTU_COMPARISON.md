# Ubuntu 24.04 / 26.04 regression comparison

This experiment branch starts directly at upstream AMRVAC 3.3 commit
`91f5091258cd8948f3a189a53c1c19ee15eb1974`. It does not descend from or merge
the fork's modified 3.2 branch. No simulation source, reference logs, input
parameters, tolerances, or existing workflows are changed.

The new workflow runs only in `nanami0721/amrvac_ionization`. A push changing
the harness or workflow on `codex/ubuntu-autotest-comparison` starts the two
jobs. The existing default branch is left unchanged, so GitHub may not expose
the manual dispatch button; use **Re-run all jobs** on the experiment run.
No pull request or upstream push is required.

Both jobs check out the same pinned, clean source commit into `source/` and
run the twelve reported regression cases with `ARCH=debug`, two MPI ranks,
one OpenMP thread, and the original comparison tolerances. Each numerical
comparison is checked directly: successful `make` alone is insufficient.
Cases run serially and preserve independent runtime logs. Numerical failures
do not prevent the remaining cases from running or artifacts from uploading.

Like upstream GNU CI, the environment sets `FC=gfortran-14`, `CC=gcc-14`
and installs Ubuntu's Open MPI packages. AMRVAC actually invokes `mpif90`;
we deliberately retain its default backend rather than setting `OMPI_FC`.
Its actual compiler, package versions, CPU, runner image, and linked libraries
are recorded. Thus this initial experiment compares complete runner
environments, not isolated OS kernels or compiler versions.

The artifacts contain environment/install/build records, per-case output and
reference logs, comparator output, SHA-256 hashes, JSON results, and a summary.
They are retained for 30 days. Red jobs may represent the numerical mismatches
being investigated; consult `summary.md` to distinguish them from build or
runtime errors. This is a targeted twelve-case diagnostic, not the full suite.
