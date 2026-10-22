# GNU Fortran regression reproducibility

The Ubuntu 26 GCC 15.2 build reproduces the historical passing outputs for the original twelve investigated regression cases with:

```
-ffp-contract=off -fno-tree-loop-vectorize
```

`arch/debug-repro.defs` is an opt-in profile that includes the existing `debug` profile and appends these two flags. It preserves `-O2`, floating-point traps and runtime checks. SLP vectorization remains enabled. The normal `debug` and production profiles are unchanged. Disabling loop vectorization may reduce performance; performance has not been benchmarked.

## Evidence

- Source: `91f5091258cd8948f3a189a53c1c19ee15eb1974`.
- Experiment: https://github.com/nanami0721/amrvac_ionization/actions/runs/39775757459
- Harness: `5e7d6318065fb71bf7a38aa8f06a93748ae5d8f9`.
- Ubuntu 26, packaged GCC 14.3 / 15.2, one shared OpenMPI 5.0.10 build, identical runtime library paths and hashes, two MPI ranks and one OpenMP thread.
- GCC 15.2 with contraction off alone reproduced the earlier Coleman LTE failure exactly.
- Adding `-fno-tree-loop-vectorize` restored all 12/12 cases, with all twelve output logs byte-identical to the historical passing Linux group. Solver sources, inputs, reference logs and tolerances were unchanged.
- The flags were injected by a compiler wrapper during the experiment. This new profile has additionally been checked by make variable expansion to preserve the original flags and append exactly those two options, including in `LINK`.
- This is the twelve-case investigation, not a rerun of the full 535-case suite.

For Coleman LTE, replacing only the vector `log10` function with elementwise scalar libm calls in the same GCC 15 binary also restored the historical output. The probe intercepted 11,696,864 vector calls across the two ranks. Disabling SLP alone did not fix it. Disabling loop vectorization only in `mod_eos*.f` also did not fix it; that scope excludes physics-specific modules such as `mod_hd_eos.f`.

This establishes a floating-point reproducibility configuration for the tested cases, not a claim that scalar arithmetic is physically more accurate. Compiler-generated FMA was absent from the inspected AMRVAC executables; this does not describe internal instructions in precompiled shared libraries.

## Use

Use a fresh build with `AMRVAC_DIR` set to the intended source checkout, the chosen GNU Fortran/OpenMPI toolchain on `PATH`, and `ARCH=debug-repro`. Existing generated binaries must be rebuilt after changing compiler flags.

For example, for the Coleman LTE case in a fresh checkout:

```sh
export ARCH=debug-repro
export NUM_PROCS=2
export OMP_NUM_THREADS=1
make -C "$AMRVAC_DIR/lib" -j 4 1d1 ARCH=debug-repro
make -C "$AMRVAC_DIR/tests/hd/Coleman_1D" -j 1 -f test.make hdr coleman4_LTE_HHe.log ARCH=debug-repro
```

For the full test suite the usual `make` in `tests/` can use `ARCH=debug-repro`; full-suite coverage with this profile remains untested. The scope and provenance above should accompany any proposed upstream change.
