# Memory-bounded low-beta equilibrium scan

Use the shared hydrostatic definitions: case 3 is a uniform Cartesian field with
uniform rho=1,p=0.5, no gravity; case 2 is a potential dipole plus hydrostatic gas.
Both fields tilt 20 degrees from the global polar axis. Magnetic units have mu0=1.
`beta_reference=2*p_base/B_pole**2` sets B amplitude while keeping gas unchanged.
For case 3 this is the initial beta everywhere; for case 2 it is the analytic
surface-pole beta, not uniform volume beta. The sampled minimum is recorded.
The default beta=4 exactly retains the previous B_pole=0.5 configuration.

Per patch mesh 16x20x52, blocks 8x4x4, TVDLF/woodward, SSPRK3, GLM alpha=0.5,
CFL 0.3. Compare Yin-Yang (2 ranks) with spherical (1 rank) in the same common
sector, plus global and boundary/seam diagnostic regions. Analytic radial and
spherical angular ghosts; overset angular exchange for Yin-Yang. No CT, Boris,
pressure repair, physical heating/cooling or conduction. This tests equilibrium
accuracy under fixed analytic boundaries, not open-corona behavior.

Scan beta=1,0.1,0.01,0.004, short duration 0.1 sound crossing with outputs every
0.02 crossing. Run all short cases first. Complete positive finite runs may
continue to lower beta even if engineering errors fail. A numerical run failure
blocks smaller beta for that grid/case; monitor failure stops the entire queue.
Extend matched spherical/Yin-Yang pairs to 1 sound crossing (output every 0.1)
only if both short runs meet existing criteria at every saved common-region time:
Mach RMS<0.01, max<0.05, rho/p/B relative L1<0.02. Positivity is checked globally.
An error-threshold failure is distinct from a solver crash. No tolerance changes.

Mach uses reference sound speed sqrt(gamma*0.5); Alfven Mach uses the local
analytic initial |B|/sqrt(rho). Report both plus absolute velocity. Native raw
divB and analytic discrete-reference defects are reported, also normalized by
analytic |B| with shell thickness L=1. Normalization by changing B is not a
resolution-convergence claim. All native dumps must cover each final snapshot.

Build in an external source copy; only test definitions change. Reuse the
previous six-configuration library evidence only after core-source hash equality.
Compile default/debug; repeat debug callbacks for both cases and grids at beta
0.004, check stderr=0, and compare default beta=4 against the previous equilibrium
snapshot. Run `suite.py EXTERNAL_WORKDIR` after external build preparation.

All simulation and analysis processes run under the shared 1 GiB tree watchdog.
Compilation alone may use 4 GiB. System start/stop thresholds are 65%/68% of
anonymous+wired+physical-compressor memory; cache and swap are logged separately.
No concurrent builds or runs. The sampled watchdog is not an OS hard quota.

After short beta=0.01 failures on both grids, `refinement.py` performs one-direction
doubling (radial, theta, phi separately) for the uniform-field case on each grid
at 0.1 crossing. Each run doubles cells, not an 8x full-grid refinement. These
directional sensitivity checks do not establish full-grid convergence order.

Completed on 2026-11-17: all runs remained positive/finite, but short beta=0.01
and 0.004 failed engineering accuracy on both grids. At beta=0.1 the dipole
case also failed the 1-crossing Mach target. See doc/yinyang_low_beta_validation.md.
`refinement.py --time-step-control EXTERNAL_WORKDIR` runs the original grid at
CFL=0.15 for a time-step isolation control. No thresholds are relaxed.
