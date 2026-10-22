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

The first explicit-backend run (39716141887) ran no numerical cases:
Ubuntu 24/GCC 14 encountered a probe name collision with MPI_Probe, and
Ubuntu 26/GCC 13 and 14 could not read the system-provided mpi.mod format.

The current follow-up fixes the probe name and builds the same OpenMPI
5.0.10 tarball with each selected C/C++/Fortran compiler. Its SHA-256 is
pinned from the official HTTPS download; bundled PMIx, PRRTE, hwloc and
libevent are used. Optional UCX/libfabric/verbs/CUDA are disabled. Build
configuration and logs are retained. No system MPI module is reused.
The matrix includes Ubuntu 24/GCC 13,14 and Ubuntu 26/GCC 13,14,15.
The 13/24 and 15/26 jobs bridge the earlier system-MPI baselines and check
whether changing the MPI build changes the observed regression outcomes.

Each job verifies mpif90's actual backend major version and compiles/runs
a two-rank use-mpi reduction before the numerical cases. Compatibility
failures are not numerical mismatches. Within each OS, compiler versions
are varied; across OS versions, CPU models, compiler minor versions,
glibc, and compiler runtime packages can still differ. This matrix alone
cannot uniquely attribute a difference to the OS or compiler. The complete
source checkout remains pinned and unchanged.

The artifacts contain environment/install/build records, per-case output and
reference logs, comparator output, SHA-256 hashes, JSON results, and a summary.
They are retained for 30 days. Red jobs may represent the numerical mismatches
being investigated; consult `summary.md` to distinguish them from build or
runtime errors. This is a targeted twelve-case diagnostic, not the full suite.

Run 39716632978 completed the Ubuntu 26 cases: GCC 13 passed 12/12;
GCC 14 and 15 each had 12 numerical mismatches. Both Ubuntu 24 jobs
compiled, but all cases aborted during MPI initialization with SIGFPE in
hwloc_libxml_backend_init at topology-xml-libxml.c:363. They produced no
numerical result. The original probe lacked AMRVAC's floating-point traps.

The follow-up applies HWLOC_LIBXML=0 uniformly to use hwloc's built-in XML
parser, preserving AMRVAC's floating-point traps. The preflight probe now
also uses -ffpe-trap=invalid,zero and -fcheck=all. All five combinations are
rerun to validate the workaround and maintain a common runtime setting.
