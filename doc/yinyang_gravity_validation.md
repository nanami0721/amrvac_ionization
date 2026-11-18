# Yin–Yang central gravity and equilibrium validation

Validation date: 2026-11-17. All planned gates passed on the local development branch.

Source: `/Users/nanami/workspace/2608_uawsom/yinyang-recovery-20261116/source`; base commit `6f5e6f10e5da5719dfc160da4591ae369450b603`, including uncommitted changes. No commit or push. All builds and data are external.

## Evidence and scope

- Standard O3 build and complete debug libraries for 1D1V, 1D2V, 1D3V, 2D2V, 2D3V, 3D3V passed. GCC garbage-collection parameters reduced compilation overhead without changing numerical options.
- Three repeated debug cases passed with empty stderr. Central gravity was checked with nonzero velocity for coordinate rotation, radial momentum and energy work.
- Uniform-state 2/4-rank comparison, forward transport against the previous baseline, diagnostics on/off, equilibrium 2/4-rank comparison, and restart at 1 tau through 3 tau versus continuous evolution: all maximum snapshot differences were zero. Reverse transport also completed with empty stderr.
- All three separated equilibria on both grids passed at 1 and 10 sound-crossing times; radial refinement passed at 1 crossing.

## Ten-crossing errors in the shared interior sector

| Grid | Case | Mach RMS | Density L1 (%) | Pressure L1 (%) | Magnetic L1 (%) |
|---|---|---:|---:|---:|---:|
| spherical | 0 | 0.000232391 | 0.4866 | 0.0433 | zero field |
| yinyang | 0 | 0.00011309 | 0.5222 | 0.0680 | zero field |
| spherical | 1 | 0.000478255 | 0.0240 | 0.0370 | 1.5742 |
| yinyang | 1 | 0.00037553 | 0.0239 | 0.0357 | 1.6936 |
| spherical | 2 | 0.000536471 | 0.4151 | 0.0489 | 1.8157 |
| yinyang | 2 | 0.000327048 | 0.4565 | 0.0446 | 1.8631 |

Cases: 0 = nonmagnetic hydrostatic gas; 1 = potential dipole in uniform gas; 2 = hydrostatic gas plus potential dipole. The common region is r=1.2..1.8, global theta=60..120 degrees, phi=100..190 degrees. Criteria: Mach RMS <0.01, Mach maximum <0.05, density/pressure/B relative L1 <2%; B=0 uses absolute magnetic RMS <1e-12. These are engineering acceptance criteria, not a convergence-order claim.

## Radial refinement and seam comparison

At one crossing, doubling radial cells 16 to 32 (angles unchanged) reduced Mach RMS 0.000834189 to 0.000180061, density L1 0.0846% to 0.0116%, and magnetic L1 0.5351% to 0.1295%. This supports a substantial radial-discretization contribution; two radial resolutions do not establish multidimensional asymptotic order.

- spherical: magnetic L1 in the common seam bands 1.8987%, off-seam 1.7270%. The spherical bands use the same geometric mask; they are not spherical grid interfaces.
- yinyang: magnetic L1 in the common seam bands 1.9377%, off-seam 1.7892%. The spherical bands use the same geometric mask; they are not spherical grid interfaces.

## Budgets and interpretation

Maximum normalized discrete mass/energy closure residual across the analyzed runs: 1.5063e-13. This checks accounting against all-direction weighted numerical fluxes and measured source increments. It does **not** prove strict global conservation of the overset grid. Radial flux alone excludes angular boundary fluxes and overlap-weight effects; the spherical control covers a different physical domain, so global totals are not directly comparable. Gravity work and other source energy are separately recorded.

The validated scope is fixed-grid 3D total-energy ideal MHD with GLM, central gravity, and analytic equilibrium boundaries. It does not validate open solar-wind boundaries, a driven corona, CT, B0 splitting, AMR, or performance scaling. This is a stable equilibrium milestone, not a complete solar-corona model.

## Memory and reproducibility

- Physical RAM: 64 GiB. System accounting: anonymous + wired + physical compressor pages, with file cache and swap separately logged. Start threshold 65%, stop threshold 68%; the 0.2 s sampled watchdog is not a kernel hard limit.
- User-authorized compilation cap: 4 GiB; all simulation/analysis caps: 1 GiB. Successful complete build peak RSS 1685.45 MiB; largest simulation peak RSS 551.06 MiB (four ranks); radial-refined run 454.48 MiB. Successful full build/test system peak about 52%; sampled swap usage remained zero.
- Initial 1 GiB compilation attempts were safely terminated; VACPP header-length and obsolete gamma-interface errors were corrected before the successful build. Startup, synthetic system-runtime termination, and small task-runtime termination paths were exercised.
- Exact parameters and executable hashes: per-run parameters.json/run.json and amrvac.par. Compiler options, source diff, hashes, build provenance, matrix.json, suite.json, native diagnostic dumps, analysis.json and guard samples are retained in the external task directory.

External task directory: `/Users/nanami/workspace/2608_uawsom/yinyang-gravity-20261117`.
