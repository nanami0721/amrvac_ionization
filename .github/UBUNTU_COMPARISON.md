# Ubuntu 24.04 / 26.04 regression comparison

This experiment branch starts directly at upstream AMRVAC 3.3 commit
`91f5091258cd8948f3a189a53c1c19ee15eb1974`. It does not descend from or merge
the fork's modified 3.2 branch. No simulation source, reference logs, input
parameters, tolerances, or existing workflows are changed.

The new workflow runs only in `nanami0721/amrvac_ionization`. A push changing
the harness or workflow on `codex/ubuntu-autotest-comparison` starts the matrix
jobs. The existing default branch is left unchanged, so GitHub may not expose
the manual dispatch button; use **Re-run all jobs** on the experiment run.
No pull request or upstream push is required.

All jobs check out the same pinned, clean source commit into `source/` and
run the twelve reported regression cases with `ARCH=debug`, two MPI ranks,
one OpenMP thread, and the original comparison tolerances. Each numerical
comparison is checked directly: successful `make` alone is insufficient.
Cases run serially and preserve independent runtime logs. Numerical failures
do not prevent the remaining cases from running or artifacts from uploading.

The initial baseline run (39713238936) left `OMPI_FC` unset, matching upstream:
Ubuntu 24.04 used gfortran 13.3 and passed all twelve cases; Ubuntu 26.04 used
gfortran 15.2 and failed all twelve numerical comparisons. Its configuration
is preserved in commit `69171b29d768417bc2913454ae44957302ba55eb`.

The current follow-up explicitly sets `OMPI_FC` and verifies `mpif90`'s actual
compiler major version. It adds Ubuntu 24.04 / GCC 14 and Ubuntu 26.04 / GCC
13 and 14. A small `use mpi` program checks that the selected compiler can
read the prebuilt MPI module and run a two-rank reduction. Compatibility
failures stop before the numerical tests and must not be counted as numerical
mismatches. Ubuntu's MPI packages are retained, and compiler/package/CPU/image
and linked-library versions are recorded. CPU models can vary between hosted
jobs, and compiler minor versions can differ across Ubuntu versions; these
remain limitations when interpreting the comparison.

The artifacts contain environment/install/build records, per-case output and
reference logs, comparator output, SHA-256 hashes, JSON results, and a summary.
They are retained for 30 days. Red jobs may represent the numerical mismatches
being investigated; consult `summary.md` to distinguish them from build or
runtime errors. This is a targeted twelve-case diagnostic, not the full suite.
