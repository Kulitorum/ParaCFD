# Final recovery attempt — 6 September 2026

**The flown PlanB lift is still unresolved.** This attempt fixes two additional numerical defects and improves the flow display. It does not establish that the simulated case can support a person. The supplied neutral-trim placement, relative-flow direction, STEP and original case JSON remain unchanged. There is no lift multiplier, prescribed lift or added PlanB rotation.

## Additional repairs

### Pressure was applied at the wrong locations

The corrected pressure solver already reproduced a linear pressure field on its flux connections. The cell momentum update did not: it used the cell pressure at wall patches and interpolated internal face pressure without the displacement to the actual aperture centroid. These are different locations on cut cells. The new production test measured **0.815970** error for an acceleration that should be exactly 1, despite the pressure solve's own derivative test passing.

Internal pressure traction now includes the displacement from the interpolated cell-centre location to the face centroid. Wall pressure is reconstructed at the wall centroid with the checked same-fluid pressure gradient. The GPU fluid impulse, CPU reference and reported body traction use this same quadrature. Each internal face retains a shared pressure and equal/opposite impulse; surface force is the reaction to the fluid traction.

On the inclined-box fixture, the worst new linear-pressure cell error is **0.000525** after the repair. The integrated linear-pressure force agrees with the displaced-volume identity to **3.55e-8 N**; pressure/body reaction mismatch is below **8.5e-9 relative**. The separate full-PlanB operator test also passes: worst linear-pressure cell error **0.001714**, pressure/body reaction mismatch **3.24e-8 relative**, and affine face-velocity error **6.44e-7 m/s**. These are operator tests, not flight validation.

### Small-cell stabilization erased the existing velocity field

The previous conditioning step averaged the complete velocity state within a group after advection, pressure, diffusion and wall updates. It therefore erased resolved gradients even when the applied update was zero. A fixture with **270 groups / 540 member cells** changed an affine velocity field by **0.0148519 m/s** in all four zero-force/tiny-timestep checks.

Stabilization now redistributes only the integrated momentum increment. Each cell retains its previous velocity and receives the group's conservative increment per unit volume. The implicit wall and pressure updates retain a pre-update state for this operation. The same tests now preserve the field exactly at stored precision. This removes one unphysical damping mechanism; it does not validate the turbulence model.

### Display and diagnostic corrections

The GUI's uniform flow snapshot previously always sampled at the coarsest AMR spacing. It now selects the finest spacing that fits a **2,097,152-cell** display budget, increasing by factors of two when necessary. Arrows, tracers, derived volume fields and recordings consequently retain more of the resolved field. Live AMR slices already sampled the fine solution. No synthetic turbulence was added.

The command-line force summary now reports pressure-only coefficients correctly when wall shear is disabled, instead of selecting unavailable total-force coefficients. Earlier section logs containing `Cl=nan` still contain valid dimensional pressure forces; their coefficients below are explicitly calculated from those forces and the slab reference area.

## Independent check of the unchanged input

`section_panel_check.py` implements a separate constant-source/vortex panel calculation with a trailing-edge Kutta condition. It shares no flow solver or force integration code with ParaCFD. The formulation follows the [Hess–Smith panel equations](https://web.itu.edu.tr/~atares/courses/CA/4.2_HessSmith.html).

Its symmetric NACA0012 zero-angle result is zero; the ±4° cases have equal/opposite lift. NACA2412 at +4° gives **CL 0.74076**, versus the existing sharp-profile reference **0.7414**. Pressure and circulation lift agree and spurious pressure drag is small.

The exact exported PlanB centre-section polyline, translated/scaled but **not rotated**, gives:

| Panels | CL in horizontal flow | Spurious pressure CD |
|---:|---:|---:|
| 200 | -0.104425 | 0.004883 |
| 400 | -0.107124 | 0.002520 |
| 800 | -0.108369 | 0.001473 |
| 1600 | -0.108351 | 0.000791 |

The finest two lift values differ by 0.000018. The finest pressure/circulation lift difference is 0.000281. Two XFOIL direct-panel checks also give approximately -0.1068 and -0.1053. A higher-count XFOIL attempt becomes inconsistent and is rejected. None of these calculations represents viscosity, separation, the whole arched wing, or a complete aircraft.

The span audit is significant: centre and quarter-span sections predict negative lift, while sections at 10% and 90% of projected span predict positive lift. A centre-section result therefore cannot prove the whole-wing lift. A full CAD-to-AVL approximation was attempted using [MIT's AVL](https://web.mit.edu/drela/Public/web/avl/). It produced inconsistent near-field/wake forces, negative induced drag and spurious asymmetry near the folded tips, and is explicitly **rejected**. It is not a validated independent 3D comparison. The raw results and rejection reasons are retained.

## Verification and remaining accuracy gap

All **nine CTests pass**, including the new conditioned-momentum and linear-pressure tests. The final NACA2412 +4° inviscid benchmark passes its unchanged gate: mean **CL 0.71223**, reference **0.7414**, approximately **3.9% low**; force-window drift **0.292%**. Circulation CL is **0.64611**, so the remaining pressure/circulation discrepancy still needs resolution.

The previous unmerged h=0.015625 m PlanB thin-section run diverged at step 78; halving the timestep also diverged. After the pressure-location repair, a 200-step regression at the same fine resolution remained bounded at approximately 15 m/s. This short check preceded the increment-redistribution repair; sustained final-method evidence is the separate merged run below.

Final-method thin-section results at t≈1.2 s, U=10 m/s, rho=1.225 kg/m³, near-inviscid nu=1e-12, Cs=0 and no wall shear:

| Finest h | Minimum volume fraction | Slab width | Fz | CL from q × chord × slab width |
|---:|---:|---:|---:|---:|
| 0.03125 m | 0.005 | 0.125 m | -0.667828 N | -0.03913 |
| 0.015625 m | 0.25 | 0.125 m | -1.332724 N | -0.07809 |

The finer run moves toward the independent section result but remains materially different. The changed merge threshold, slight span variation in the cropped CAD slab, finite domain, discretization error and incomplete approach to steady flow must be separated before calling this a mesh-convergence study. The simulated duration is only about 5.4 chord travel times. The inviscid panel result is a steady, unbounded-domain section reference.

For the complete original PlanB case, the h=0.125 m run with corrected pressure quadrature completed 1.60072 s with total force **[45.0562, 0.1843, -39.1270] N**. It has no conditioning groups, so the subsequent change to group increments does not change that path. The final h=0.0625 m run completed 1.60022 s with force **[29.7965, 0.0262, -21.2061] N**, bounded maximum cell velocity 36.55 m/s and pressure residual 9.09e-6. Raw results are in `final-planb-whole-checked-h0625.log` and `final-results.json`. These endpoint runs are regressions, not a validated lift/drag estimate with a numerical uncertainty bound.

## Build for testing

The updated GUI is:

`C:\CODE\ParaCFD\build-paraglider-ui\Release\paracfd-gui-fixed.exe`

A separately deployed copy, including Qt/OCCT dependencies, is in `build-paraglider-ui/FinalTest/`. Windows prevented replacement of the originally requested `Release/paracfd-gui.exe` because that older GUI remained running. No existing GUI process was terminated. If the requested executable is subsequently replaced, the final fingerprint file records that explicitly.

Load **`configs/planb_parakite_review.json`** to reproduce the finer numerical settings: its only differences from the original JSON are three AMR levels and pressure tolerance 1e-5. Geometry, placement, speed, viscosity, SGS coefficient, wall/default momentum mode and neutral trim are preserved. This is a review case, not a physically validated flight preset.

The separate executable links and its deployed application starts, recognizes the RTX 4090 and loads the previous geometry. The offscreen startup invocation reports zero simulation steps and returns its existing scripted-run failure code; it did not request a CFD run. **Interactive rendering, Run/Pause and recording have not been validated by that check.** CLI and operator validation do not substitute for those GUI checks.

## What is required to finish

1. **Establish the exported flight frame.** The product currently fixes the wing in horizontal +X relative airflow with slip outer side/top/bottom boundaries. It does not simulate descent. Neutral line settings and direction of relative airflow are separate quantities. For a CAD pose referenced to gravity, 10% descent in still air gives upward relative flow with inclination atan(0.1)=5.71°. If the CAD is already referenced to incoming air, adding that inclination would change the supplied AoA. The user has been asked; no offset has been assumed. [NASA defines AoA relative to the airflow.](https://www.grc.nasa.gov/www/k-12/Aero2000/studweb/2-3-2i2.html)
2. **Match the actual flight evidence.** Obtain measured neutral-trim airspeed, all-up mass, air density and the loaded shape corresponding to this export. Implement the confirmed flight-frame convention and consistent far-field conditions while preserving neutral rigging. Do not select an angle or force multiplier to obtain a desired weight.
3. **Close the numerical acceptance gaps.** Finish a matched section study with consistent merge settings, longer averaging, smaller timesteps, three grids and larger-domain checks; resolve the pressure/circulation discrepancy. Add an accepted independent finite-wing comparison before accepting PlanB's whole-wing force. The full AVL approximation in this attempt does not satisfy that requirement.
4. **Complete the flight and product checks.** Validate whole-aircraft force and moment balance, including line/pilot/harness drag where applicable, against the supplied measurements. Verify the actual GUI workflow and save a complete run manifest. The broader acceptance gates remain in `COMPLETION_PLAN.md`.

There is now stronger evidence for both numerical defects and a mismatch between the currently simulated condition and the known flight. The fact that the real wing carries people remains the validation target. It is not evidence that the exported fixed pose in horizontal flow must generate their weight.

## Reproduction and provenance

`final-recovery.patch` records this attempt's production/probe changes relative to the prior repaired source snapshot in `build-paraglider-ui/final-attempt-baseline/`; it is not a patch against a clean historical checkout. Pre-existing uncommitted work is preserved. `final-fingerprints.json` records source, input and executable SHA-256 values; `final-results.json` records the case summaries, checks and unresolved status.

From the repository root, with the configured CUDA/OCCT build available:

```powershell
ctest --test-dir build-paraglider-ui -C Release --output-on-failure
.\build-paraglider-ui\Release\physics_contract_probe.exe --conditioning
.\build-paraglider-ui\Release\physics_contract_probe.exe --planb
.\build-paraglider-ui\Release\naca_validation_probe.exe --profile 2412 --single-angle 4 --max-levels 3 --span-layers 4 --physical-time 1.2 --conservative-cell-momentum --max-wall-seconds 900 --progress-steps 1000
.\build-paraglider-ui\Release\paraglider_case_probe.exe --config configs/planb_parakite_review.json --physical-time 1.6 --sample-every 500
python review/section_panel_check.py
```

The independent AVL script requires the separately downloaded MIT executable under the ignored reference-tools directory. Its full-CAD approximation currently exits nonzero by design because its consistency checks fail. Do not promote those rejected results to a flight reference.
