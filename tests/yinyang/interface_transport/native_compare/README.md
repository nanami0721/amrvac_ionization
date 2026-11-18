# Native divergence: spherical versus Yin-Yang

The spherical test reuses the original native `get_divb` diagnostic, including
state save/ghost refresh/restore. It uses the same exact magnetic transport,
radial shell, MC/TVDLF, GLM alpha=0.5 and radial divergence closure as the
Yin-Yang `verified-fix` tests. Only the angular domain/ghost treatment differs.
The spherical domain places the old Yin longitude seam in its interior.

Build `mod_usr.t` as 3D3V in an external case with an external source copy as
`AMRVAC_DIR`. Run `run.py EXTERNAL_RUN EXTERNAL_EXECUTABLE --scale 1`, and repeat
with scale 2. Direct invocation automatically uses the existing convergence
watchdog: one rank by default, 1024 MiB process-tree monitoring, terminate on
excess or monitor failure. The watchdog is sampled, not a kernel hard quota.
Do not automatically escalate mesh size or raise the threshold.

`--debug` enables a native dump on every step for repeated-call checking;
`--no-native` disables diagnostics. Normal dumps are at t=0,0.3,0.6. Snapshots
are also saved at 0.15 and 0.45. Existing run directories are rejected.

Analyze under the same watchdog, with a 768 MiB threshold:

```sh
python ../convergence/memory_guard.py --limit-mib 768 --out EXTERNAL_GUARD -- python analyze.py SOURCE SPHERICAL_ROOT YINYANG_RUN_ROOT
```

The spherical root contains `s1/`, `s2/`; the Yin-Yang run root contains the
existing `c1-p1-d1-s1/`, `c1-p1-d1-s2/` results and native dumps. The analyzer
compares the forward magnetic case only. It reads one patch at a time and
checks dump coverage, time, native volume and an independent implementation
of the exact-field FV divergence. It does not reconstruct numerical ghosts
in Python or substitute a different divergence operator.

All compared regions lie within the same moving physical sphere of radius
0.54, intersected with r=1..2. Regions are: the full observation sphere; first
radial physical layers; fixed inner r<1.1 and outer r>1.9 bands; the sphere
excluding the first layers; and the fixed interior shell 1.2<r<1.8, subdivided
by the union of 10-degree bands around the original Yin-Yang edges. Apply the
same physical seam mask to the spherical run, where there is no actual seam.
These diagnostic regions overlap; they are not an additive partition.

Yin-Yang integrals retain 16x16 overlap quadrature; spherical volumes are counted
once. Sampling/overlap weights differ, so small numerical differences are not
an isolated interpolation-error measurement. Do not compare the two different
full domains, or treat this local peak as the earlier global peak.

Report numerical divergence, analytic-field discrete divergence, and their
difference separately. Analytic ghosts are used for the reference, so the
difference includes numerical boundary/interpolation effects. First-layer
divergence is constrained by the radial closure and is not an independent
accuracy certificate. A negative reference/defect correlation can cause a
small raw divergence on coarse grids; it does not prove GLM is the cause.
