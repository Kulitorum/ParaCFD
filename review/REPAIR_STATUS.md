Latest recovery results: [Final attempt — pressure traction, small-cell increments and independent checks](FINAL_ATTEMPT.md). The real-flight lift remains unresolved.

# Simulation repair status — 6 September 2026

**Latest follow-up:** [FLIGHT_RECOVERY.md](FLIGHT_RECOVERY.md) records additional velocity-reconstruction and gradient-stencil repairs, the final tests, and the unresolved flight-frame requirement. The numerical history below describes the preceding repair stage.

The supplied PlanB STEP and `configs/planb_parakite.json` are unchanged. No additional PlanB rotation, AoA sweep, force multiplier, prescribed lift, or glide-angle correction has been applied. The pre-repair working tree is preserved under `build-paraglider-ui/repair-baseline-20260906/`.

## Implemented numerical repairs

- **R1, vanishing-timestep consistency:** removed the finite embedded-aperture velocity remapping that ran independently of timestep. The transport update now approaches identity as dt approaches zero.
- **R2, pressure consistency:** enabled the complete nonorthogonal gradient in the production pressure operator and flux correction, including AMR interfaces. The nonsymmetric corrected operator uses flexible GMRES with a multigrid preconditioner. Removed the unadvertised pressure-tolerance floor and allowed tighter pressure tolerances in the GUI.
- **R3, one geometric envelope:** cut-cell volumes, apertures, material classification and pressure patches now use the same triangle envelope. Full-envelope classification survives cropped display/flow geometry. Replaced fitted CAD-face planes with exact finite-facet partitions, moved local clipping arithmetic into cell coordinates, stabilized cap construction/volume integration, and prevented clipping intersections from extrapolating beyond their edges. Assign pressure patches by geometric overlap instead of nearest facet centroid. The core rejects nonclosing cell surfaces before allocating the timestep solver.
- **R4, pressure reaction:** the conservative cell path applies the same finite-volume pressure tractions used for reported surface loads. Removed the pressure-dependent geometry-closure compensation. Retained redistribution preserves the integrated impulse. Boundary tractions are included in cropped-domain budgets.
- **R5, solid masking:** classify every connected regular material region instead of assuming all atlas/crop boundaries lead to exterior fluid. Classify and exclude active bricks lying wholly inside the solid, including levels with no surface atlas.
- **R6, shared selection and checks:** sweeps use the selected momentum solver. Added direct production GPU pressure/transport/reaction tests, cropped-solid and interior-brick tests; corrected the airfoil test's sealed-fluid and pressure-tolerance expectations. Conservative cell momentum is now the default across the core, GUI and command-line probes; staggered MAC remains available explicitly for diagnostic comparison. Result status does not imply physical validation.
- **R9, trim preservation:** restore the original placement after trial completion, cancellation and failure; saving during a trial writes the baseline placement.
- **Additional composed-update defect:** the conservative predictor previously restored the last pressure's face/cell difference and then solved for an absolute new pressure. That mixed two pressure-update formulations. The nonincremental predictor now interpolates its pressure-free provisional momentum and applies the newly solved pressure once.

## Verified evidence

A fresh Release build, including the GUI, passed all eight CTests, including 64 evolved conservative-momentum steps and an analytic AMR circulation-contour test. The full, unchanged PlanB production-contract test also passed:

| Check | Result |
|---|---:|
| Maximum cell area-vector closure defect | 1.1305e-10 m² (7.235e-9 h²) |
| Sealed interior fluid DOFs | 0 |
| GPU constant-pressure derivative error | 0 |
| GPU linear-pressure derivative error, worst axis | 2.981e-7 |
| Pressure-force / fluid-impulse relative mismatch | 4.2721e-8 |
| Velocity change at dt=1e-12 s, after eight evolved steps | 1.1642e-10 m/s |

These checks establish specific discrete identities. They do not establish aerodynamic accuracy, long-time stability, moment/energy conservation, or a flight-performance prediction.

## End-to-end evidence and remaining failures

Before the final nonincremental-predictor repair, conservative NACA2412 at +4° produced mean CL 0.65890 at 64 cells/chord and 0.69749 at 128 cells/chord, compared with the recorded XFOIL inviscid reference 0.7414. Both runs completed 0.8 s but failed the existing acceptance criteria (including force-window drift). These are convergence evidence, not accepted validation. A longer MAC airfoil run lost lift and hit its wall-time limit. The first conservative full-PlanB run became unstable near t=0.54 s and was stopped; its large subsequent forces are invalid. The old fine-airfoil run also ended with a CUDA illegal-access error and is not accepted. The subsequent predictor and diagnostic retests are recorded below.

The final nonincremental predictor completed **1.00053 s / 1,031 steps** on unchanged PlanB with bounded velocity (prestep maximum 28.40 m/s). Total force at the end was **[44.1723, 0.1137, -36.3466] N**. Thus the numerical blow-up was removed over this regression interval, but this horizontal-flow case still does not reproduce the requested positive-lift flight state. It is not a converged flight-load estimate.

With the corrected interpolated circulation diagnostic, **NACA2412 +4°, inviscid, 128 cells/chord, t=1.2 s passed the existing acceptance gates without relaxing them**: mean CL **0.71062** versus recorded reference **0.7414** (4.15% lower); current/previous-window drift **0.281%**; RMS CL **0.000312**; circulation-based CL **0.64465** versus final pressure CL **0.71106** (difference 0.06641, below the pre-existing 0.07 limit); pressure residual **7.65e-6**. The new AMR contour test checks both known nonzero curl and zero curl at off-grid contour positions. This one accepted inviscid comparison is not the complete airfoil/finite-wing validation matrix.

A 256-cells/chord run of the repaired predictor also completed 1.2 s with mean CL **0.71591** (3.44% below the reference), bounded maximum 13.95 m/s and 0.326% force-window drift. It still used the old nearest-cell circulation diagnostic and is recorded as a failed overall test, not retroactively relabeled as passed. Its mean lift differs from the 128-cell result by about 0.74%. The earlier CUDA illegal-access error did not recur in this 9,004-step run.

A separate XFOIL check of the exported PlanB section was attempted at zero added rotation. Its panel results were not convergent between 240 and 480 nodes and had large spurious pressure drag; they are rejected and are **not** evidence for a PlanB lift prediction. The downloaded reference executable is confined to the ignored build directory. A NACA2412 sanity check with that executable produced CL 0.7383 at +4°.

The coarse PlanB case has about 18 cells along its midspan chord. It is a regression case, not an adequate independent mesh-convergence study of a wing with this trailing edge and surface detail.

## Flight-frame question, with the geometry frozen

The current solver holds the wing fixed and imposes horizontal airflow. It does not evolve wing/pilot mass under gravity or develop a sink velocity. The approximately 10% sink condition therefore is not present in the current boundary conditions. The user has been asked whether the exported neutral orientation is relative to the horizon or already relative to the incoming air/flight path, and for measured all-up mass and trim airspeed. No answer has been assumed.

`physics_contract_probe --section-planb` exports the unchanged midspan section to `review/planb-neutral-section.csv`. The accompanying plot shows a strongly reflexed mean line; a rough thin-airfoil integration at horizontal airflow gives a small negative section-lift estimate. Because the section is thick and this is a finite 3D wing, that estimate is only a sign/convention diagnostic, not validation or evidence to change the supplied trim. [MIT thin-airfoil notes](https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/d721171c42af48a056aaccec784a6d10_f03_0304.pdf) explain the mean-line dependence. [NASA's relative-wind definition](https://www.grc.nasa.gov/www/k-12/Aero2000/studweb/2-3-2i2.html) distinguishes the flight path from horizon-relative attitude.

## Remaining project gates

1. Extend the passing one-second PlanB stability regression and the accepted NACA comparison to longer runs and the full pressure/transport coupling matrix. Diagnose the fine-case CUDA failure with a reproducible sanitizer run if it recurs.
2. Finish the uniform/AMR circulation, momentum, energy and moment budgets; establish affine-preserving small-cell stabilization. The retained MAC variable-viscosity stress and explicit diffusion timestep still need repair/validation.
3. Complete independent airfoil reference comparisons, three meshes, two timesteps, domain/tessellation sensitivity and a 3D finite-wing benchmark. Do not promote a settled force trace to physical validation.
4. Confirm the flight-frame metadata while preserving neutral line trim; then model and compare the measured glide state, including line/pilot/harness drag and weight.
5. Persist solver/wall options and complete run manifests, GUI/CLI equivalence and packaging checks, as specified in `COMPLETION_PLAN.md`.

## Reproduction

```powershell
$env:CUDA_PATH = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.1'
$env:CudaToolkitDir = $env:CUDA_PATH
cmake --build build-paraglider-ui --config Release -- /m
ctest --test-dir build-paraglider-ui -C Release --output-on-failure
$env:PARACFD_EB_THREADS = '16'
.\build-paraglider-ui\Release\physics_contract_probe.exe --planb
```

Raw diagnostic histories are retained in `review/repair-*.log` (ignored build/run artifacts). Do not interpret initialization pressure impulses or failed/time-limited runs as flight loads.

The final build/source hashes are in `review/repair-fingerprints.json`. `review/repair-only.patch` compares the repaired code with the preserved pre-repair working tree, excluding the user's earlier edits. The final GUI was built; interactive GUI lifecycle tests were not automated.
