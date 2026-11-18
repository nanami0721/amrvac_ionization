# Memory-bounded convergence diagnostics

Keep source definitions here, and build/run/analyze only in an external workspace.
Use an external AMRVAC source copy as `AMRVAC_DIR`. The spherical module is 3D3V;
the axial Cartesian module is 1D3V (`Cartesian_1.75D`). In 1D the framework does
not allocate the GLM psi field; the axial solution has Bx=0.

`run.py` defaults to **one MPI rank**. Direct invocation automatically enters
`memory_guard.py` with a **1024 MiB process-tree threshold**. Existing output
directories are rejected. For example, with an externally built executable:

```sh
python run.py /external/task/run /external/task/amrvac --grid spherical --scale 1
```

The guard records aggregate RSS and macOS physical footprint every 0.2 seconds.
It terminates the task process group and verified descendants on a threshold
breach or monitoring failure. It is a sampling watchdog, **not a kernel hard
memory quota**; transient overshoot is possible. Guard metadata for a direct run
is saved beside the run in `<run>-memory`. It requires macOS process-inspection
permissions and fails closed when it cannot monitor the task.

Use the same guard for post-processing:

```sh
python memory_guard.py --limit-mib 768 --out /external/task/analysis-memory -- python analyze.py SOURCE /external/task/run
```

Do not increase mesh size from physical-memory capacity alone. Measure a complete
smaller run first, keep substantial headroom, and do not automatically raise the
threshold after a termination. A killed run without `run.json` is incomplete and
must not enter convergence fits. Keep the guard summary and partial logs.

The 2026-11-16 recovery tests measured about 531 MiB for spherical s=2 with one
rank, versus 1968 MiB with four ranks. The four-rank half-step run crossed 2 GiB
and was terminated. A one-rank queue including half-step and wider-profile
controls plus analysis completed at about 543 MiB. These are local measurements,
not general AMRVAC memory guarantees. The former unmonitored s=4/6 queue is disabled.

The baseline remains TVDLF/MC, three-stage stepping, width 0.18 and t=0.6.
`--width`, `--courant`, `--no-divbfix`, `--extended`, and explicit `--ranks` support
separate diagnostic runs. `analyze.py` reports the original moving-sphere L1/L2,
its interior-shell subset, and fixed radial regions. Do not equate the 1D error
normalization with the spherical volume integral. Wider profiles diagnose
resolution dependence; they do not replace validation of the original profile.
