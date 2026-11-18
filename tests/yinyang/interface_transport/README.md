# Forward/reverse interface transport accuracy

This suite validates the current fixed-grid Yin–Yang ideal MHD+GLM implementation.
The suite includes test diagnostics and radial boundary controls; the core ghost
update orders Yin–Yang exchange before radial divB closure when enabled. All builds,
simulations, logs, analysis output and plots belong outside the source repository.

## Two exact solutions

Define a right-handed frame e1 (equatorial tangent), e2 (global z), e3 (equatorial
normal) at longitude phi0. The reference centre is 1.5 e3. Direction d = +1 or -1,
velocity v=d e1, moving centre xc=1.5 e3+d(t-0.3)e1, width a=0.18, t in [0,0.6].

1. **Entropy pulse:** rho=1+0.5 exp(-|x-xc|^2/a^2), p=0.1, and uniform
   B=0.2e1+0.1e2+0.3e3. Reverse motion changes velocity and starting centre; it
   does not rotate or reverse the magnetic field.
2. **Force-free rotational magnetic pulse:** rho=1, p=0.1,
   alpha=(pi/2) exp(-[(x-xc).e1]^2/a^2),
   B=0.3(cos(alpha)e2+sin(alpha)e3). The strength is constant and B.e1=0.
   Therefore div B=0, curl B=-alpha' B, and J x B=0. Uniform translation is an
   exact finite-amplitude ideal MHD solution, including momentum and total energy.
   No passive-field or small-amplitude approximation is used. This is a planar
   rotational structure extending transversely, not a compact 3D flux rope/loop.

For seam runs phi0=146.25 degrees (Yin's upper longitude edge). At t=0.3 the
centre crosses that edge; d=+1 leaves Yin and d=-1 enters Yin from Yang. The
control has phi0=0, so the centre's neighbourhood lies well inside Yin. The
magnetic pulse is an extended plane: other parts still intersect patch edges
in the control. Only the local observation region is an interior control.

Radial ghosts follow the exact solution at the current boundary time. Angular
ghosts use the actual overset exchange. The shell is r=1..2. Both patches use
16s x 20s x 52s cells for s=1,2,3, with unchanged angular bounds
[33.75,146.25] degrees in theta and [-146.25,146.25] in phi. Blocks are 8x4x4.
The seam directions/cases each have three resolutions; interior controls use s=2.
TVDLF/minmod, three-step, CFL=0.3, GLM, FI ideal-gas EOS; 4 MPI ranks. Values are
sampled at cell centres both for initialization and for analytical comparison.

## Diagnostics and acceptance questions

- Compare at t=0,0.3,0.6. Initial errors should be at floating-point roundoff.
- Require finite fields and positive density/pressure, empty debug runtime stderr.
- Report whether error falls with refinement in **both** transport directions;
  estimate observed orders, without assuming second order for minmod-limited peaks.
- For entropy use |rho-rho_exact| normalized by exact density excess. For magnetic
  transport use |B-B_exact| normalized by |B_exact-0.3e2|. Report global L1/L2 and
  local L1 in a sphere of radius 3a around xc. These different tests' normalized
  errors must not be interpreted as the same physical error.
- Report the contribution to global L1 from fixed 10-degree angular bands on BOTH
  sides of BOTH patch boundary curves (union in physical space). Selecting only
  interior edge rows would bias inward/outward comparisons. Cell-centre membership
  approximates this band; overlap volume weights use 16x16 angular quadrature.
- Report Cartesian B error, |B| error, velocity/density departures and extrema.

The divergence diagnostic includes native angular edge cells. Reconstruct one
angular ghost layer by Cartesian-vector interpolation from the other patch,
then apply centred spherical derivatives. Omit only the first and last radial
layers. Apply the SAME operation to the sampled analytic field at each time and
report (1) numerical divB RMS, (2) analytic-sample divB RMS, and (3) RMS of their
pointwise difference. The latter isolates the evolution defect relative to this
diagnostic baseline; it does not erase the baseline error. This is not the
solver's own discrete divergence operator and not a CT proof. Also report
seam-band RMS/max. No h factor is used to manufacture convergence of the metric.

## Reproduction

Copy mod_usr.t to an external case directory; build with AMRVAC_DIR pointing to
an external source copy containing the deduplicated Yin–Yang implementation.
`run_suite.py RUN_ROOT EXECUTABLE --scale 1` runs four seam cases. Repeat with
`--scale 2` and `--scale 3`; use `--scale 2 --control` for four interior controls.
Use a separate RUN_ROOT with a debug executable and `--debug` for three-step
checks. Existing case directories are rejected to avoid mixed-run output.

`analyze.py SOURCE RUN_ROOT` writes results.json in RUN_ROOT's parent. It uses
NumPy and SciPy. Output generation must also remain outside the repository.
Do not infer performance from these runs: diagnostic work can overlap simulation,
and this suite is an accuracy study rather than a controlled timing benchmark.

## Diagnosis controls (2026-11-16)

`run_suite.py` additionally accepts `--limiter woodward`, `--direction 1`,
`--glm-alpha 0.1`, and `--extended-shell`. The defaults retain the original
minmod, both directions, alpha=0.5 and r=1..2 configuration. Use a **separate
external run root for each variant**. Extended shell uses r=0.75..2.25 and 24s
radial cells, preserving the original radial spacing. This control moves radial
boundaries; compare the same physical region rather than full-domain averages.

`investigate_divergence.py SOURCE RUN_ROOT EXTERNAL_OUTPUT.json --scales 1 2 3`
examines the forward magnetic seam case at t=0 and 0.6. It reports patch-edge,
seam-band, fixed interior-shell and original-shell metrics, peak locations,
and exact-field ghost interpolation errors. Replacing interpolated analytic
angular ghosts by exact analytic ghosts separates ghost interpolation from
interior derivative truncation in this **postprocessing operator**. It does not
replace solver evolution or measure solver-native `get_divb`. Ghost-field RMS
is an unweighted boundary sample RMS; volume divergence RMS uses overlap weights.
The `original_shell_seam` mask excludes one original radial layer at each end,
and is the same physical mask for the moved-boundary comparison.

See `doc/yinyang_accuracy_diagnosis.md` for the controlled results and limits.

The diagnosis runner also accepts `--case 0` or `--case 1` to isolate an exact
solution, and `--courant 0.15` for the half-time-step control. The default is
both cases and CFL=0.3. After the November 16 reboot, the original temporary
artifacts were lost; reconstructed source and new diagnostic runs are preserved
in `/Users/nanami/workspace/2608_uawsom/yinyang-recovery-20261116/`.
The new diagnostic suite is forward-only; it does not replace the full
forward/reverse acceptance suite.

## Native divergence and radial boundaries

`--radial-divbfix` enables AMRVAC's existing discrete radial normal-field closure
on the two radial special boundaries (`boundary_divbfix=2*.true.,4*.false.`).
The Yin–Yang ghost update now exchanges angular ghosts **before** this closure,
and exchanges again afterwards to propagate corrected radial/angular corners.
Initialization takes a world-wide logical OR of the need for this extra pass,
so all ranks execute the same collectives. With no active radial special-boundary
closure there is still only one overset pass.

The alternative `--radial-bc 1` extrapolates the outgoing GLM *perturbation*
relative to the analytic solution and sets the incoming perturbation to zero.
In outward-normal coordinates, Wout=psi+ch*(Bn-Bn_exact). The ghost values use
psi=Wout/2 and Bn=Bn_exact+Wout/(2ch); magnetic energy is adjusted to preserve
thermal pressure. At initialization ch=0 and the exact state is retained.
This is a reference-solution test boundary, not a general nonreflecting MHD or
solar boundary. It must not be combined with `--radial-divbfix`. The tested
linear-extrapolation prototype is not retained as a supported option.

Native diagnostics are enabled by default in this test, and disabled with
`--no-native`. At t=0,0.3,0.6 (every step for debug checks), the diagnostic saves
and restores the full owned state around a ghost refresh, then writes the actual
`get_divb` values for **all physical cells including radial edge cells**. A
reusable state-copy workspace costs additional memory; disable diagnostics for
performance runs. Fixed block ownership is required. Files must be in a fresh
run directory; existing diagnostic files are not overwritten.

Each `native_pP_rRRRR_iIIIIIII.bin` uses native-endian int32 header
`[version=1, patch, world_rank, iteration, number_of_cells]`, float64 time/ch,
then ten float64 values per cell: local r/theta/phi, numerical divergence,
analytic-sample divergence, geometric volume, psi, and the three local magnetic
component errors. The analytic reference uses exact ghosts (no interpolation).
Thus its difference from numerical divergence includes interpolation and boundary
closure effects as well as evolution errors. It is not the earlier postprocessing
operator's baseline. `analyze_native.py SOURCE RUN_ROOT` checks cell coverage,
volume, time and an independently assembled finite-volume reference, then reports
weighted global, radial-boundary and interface metrics. Only inner/outer radial
bands get incoming/outgoing GLM norms; these are not reflection coefficients.

`probe_interpolation.py SOURCE RUN_ROOT OUTPUT.json` compares four-point
bilinear and nine-point biquadratic angular interpolation of the exact magnetic
field, in Cartesian components, at three times and both directions. It is a
**static probe**, not a modified solver, conservation test or positivity test.
See `doc/yinyang_glm_boundary_validation.md` for the accepted results and costs.
