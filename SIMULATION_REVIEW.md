Latest recovery results: [Final attempt — pressure traction, small-cell increments and independent checks](review/FINAL_ATTEMPT.md). The real-flight lift remains unresolved.

# ParaCFD simulation review — 6 September 2026

The current solver has demonstrable numerical defects large enough to invalidate its lift predictions. The rigid, closed-envelope modelling choice is a reasonable starting point for external wing aerodynamics. The immediate problem is the consistency of geometry, momentum transport, pressure projection, and reported surface force. Increasing resolution or waiting for the force graph to flatten cannot, by itself, establish that these operators solve the intended equations.

This review covers the current working tree at `f1f6552` plus the pre-existing local edits, not just the committed code. I rebuilt the diagnostic executables, ran the four CTests, ran two short PlanB MAC cases and a short collocated case, attempted the current NACA benchmark with a declared time limit, and added an independent diagnostic executable. No simulation implementation or numerical setting was changed. The runnable diagnostics are in [review/](review/README.md); the recovery plan is [COMPLETION_PLAN.md](COMPLETION_PLAN.md).

**User-confirmed acceptance case:** the supplied wing is already at neutral line trim and must work with its existing AoA, with approximately 10% sink. Preserve that orientation and relative-flow angle. The earlier +12° runs were temporary command-line diagnostics, performed before this clarification; they were never saved and are not a proposed remedy or acceptance case. “0°” in the run table means zero added rotation, not a claim that the supplied airfoil has zero aerodynamic incidence. R1–R4 were measured on the unchanged PlanB placement.

## What the current code actually runs

The production chain is STEP validation and SI conversion → placed triangle mesh and OCCT solid → static Cartesian AMR → local cut-cell geometry and agglomeration → composite pressure topology → GPU momentum transport, diffusion, optional wall shear, and projection → triangle pressure/shear integration → GUI convergence and flight estimates.

There are two materially different timesteppers:

| Entry point | Momentum path | Wall treatment |
| --- | --- | --- |
| Ordinary GUI run, checkbox initially off | Staggered MAC, bounded semi-Lagrangian MacCormack away from the surface, separate near-surface/aperture updates | Spalding shear enabled |
| GUI AoA sweep and fixed-attitude glide search | Automatically selects collocated control-volume momentum | Spalding shear enabled |
| Case probe, default | MAC | Enabled unless explicitly disabled |
| NACA probe, default | MAC | Slip; molecular/SGS settings still exist |

See `external_aero_core.cpp:291`, `paraglider_window.cpp:601`, `paraglider_window.cpp:821`, and `naca_validation_probe.cpp:766`. A successful benchmark of one row does not validate the others. The current pressure/Krylov scalar is FP64, with FP32 flow fields and pressure shadow; there is now a recursive coarse hierarchy. Older descriptions of an exclusively FP32 pressure solve or a two-level-only preconditioner are stale.

## Findings, ranked by impact

### R1 — The MAC embedded-boundary update contains a finite filter at zero elapsed time

**Confirmed by executing the current GPU operator. Highest-priority explanation for loss of lift/circulation on the MAC path.**

`DeviceCompositeAmrProjection::transport_embedded_apertures` gathers advected aperture and regular-carrier velocities into component means, then replaces aperture velocities with the average of the two endpoint means. It also blends regular carrier velocities toward these means. These operations do not contain `dt`, viscosity, or wall friction:

- `amr_pressure.cu:390`: `reconstruct_compatible_eb_kernel` writes `0.5 * (ua + ub)`.
- `amr_pressure.cu:394`: `reconstruct_compatible_carrier_kernel` adds `0.5 * (mean - current)`.
- `amr_pressure.cu:1978`: this reconstruction is called on every MAC timestep.

Starting from an actual projected PlanB state, resetting that same state before each trial, and calling only this transport operation with `nu=0`, `Cs=0`, and no wall sink gives:

| State | Timestep | Maximum regular velocity change | Maximum aperture velocity change |
| --- | ---: | ---: | ---: |
| Immediately after initialization | `1e-6 s` | `1.24763 m/s` | `62.96098 m/s` |
| Immediately after initialization | `1e-12 s` | `1.24763 m/s` | `62.96092 m/s` |
| After 100 evolved steps, `t=0.2610 s` | `1e-6 s` | `0.47195 m/s` | `11.47089 m/s` |
| After 100 evolved steps, `t=0.2610 s` | `1e-12 s` | `0.47195 m/s` | `11.47097 m/s` |

This is not just the startup impulse being removed. An evolved, projected state still receives a finite change as elapsed time approaches zero. The comment calling the removed component “pressure-null circulation” is not proof that it is an unphysical mode. Divergence-free circulation is fundamental to lifting flow. A filter can conserve some definition of total momentum and still destroy the spatial velocity distribution required for lift.

**Repair direction:** establish one compatible velocity/momentum representation and a transport update that approaches identity as `dt→0` on admissible states. If stabilization is necessary, derive its consistency and conservation properties. Do not merely remove this averaging and assume the previously suppressed instability is solved. Add an evolved-state zero-time test, affine-state preservation, and a vortex/circulation transport test before comparing lift.

### R2 — The production pressure gradient fails even for linear pressure

**Confirmed from actual PlanB topology and the production coefficient formula. Affects both timesteppers.**

`amr_pressure.cpp:69` forms the compact derivative from pressure differences divided by the normal component of the centroid separation, and computes a tangential correction. `amr_pressure.cu:2135` applies that correction only when `PARACFD_PRESSURE_FULL_NONORTHOGONAL` is present. Normal runs omit it. Regular corrections are also skipped by that branch.

For a Cartesian face normal to X, two cell centroids can differ in Y and Z. If `p=y`, the exact X-normal pressure derivative is zero, but the two-point formula produces `Δy/Δx`. In the PlanB hierarchy:

- All 16,384 coarse/fine tiles produce a derivative error of magnitude **1/3** for the unit Y and Z pressure-gradient tests.
- Cut-cell faces also fail: maximum derivative errors for `p=x`, `p=y`, and `p=z` are **2.94, 3.79, and 3.37**, respectively; area-weighted RMS errors are **0.0725, 0.1057, and 0.0945**.

These are dimensionless errors relative to a unit manufactured gradient, not measured lift percentages. Pairing a symmetric matrix with the same incomplete flux correction can reduce divergence and converge PCG while still giving the wrong pressure derivative. [OpenFOAM's numerical-scheme documentation](https://www.openfoam.com/documentation/user-guide/6-solving/6.2-numerical-schemes) describes the corresponding non-orthogonal correction issue.

**Repair direction:** make the pressure operator and its face corrections pass constant/linear manufactured fields on cut cells and AMR transitions. Choose a consistent corrected discretization or an appropriate multi-point construction, then match the Krylov solver to its actual symmetry. Keep the current SPD operator as a possible preconditioner. Do not restore the expensive nested diagnostic solver as an unmeasured production fix.

### R3 — Fluid control volumes do not geometrically close, despite complete wing surface coverage

**Confirmed by the current topology's own closure calculation.**

For each control volume, the vector sum of its oriented face areas should vanish. `composite_cell_pressure_closure_cpu`, at `amr_advection.cu:922`, reports for the current PlanB solid:

- **4,119** active rows with area-vector defect above `1e-10 m²`;
- maximum defect **0.0270633 m²**;
- sum of defect magnitudes **18.9284 m²**.

A finest Cartesian face is only `0.125² = 0.015625 m²`. The maximum defect exceeds that area. This is not insignificant floating-point noise. The aggregate absolute sum is a diagnostic, not the area of a single missing hole.

The geometry builder constructs local fitted-plane/BSP fragments, while output patches retain the triangle load surface, and adjacent openings are formed from common polygon intersections (`solid_embedded_boundary.cpp:760`, `:1045`, `:1130`). Those representations must agree per cell, including merges. Merely sharing an aperture between two neighbours guarantees equal/opposite graph flux, not closure against each cell's remaining boundary.

The whole closed wing **does** pass a useful check: uniform 1 Pa exterior pressure integrates to essentially zero resultant, about `4e-15 N` per component. Thus a globally closed surface and correct overall force sign do not establish local finite-volume consistency.

**Repair direction:** derive volume, centroid, open faces, and wall area vectors from the same geometric partition. Audit closure and first moments before and after agglomeration. Reject material defects with source-face/cell diagnostics. A pressure-dependent compensation term is not a substitute for closed control volumes. The old handover's closure warning remains relevant.

### R4 — Collocated pressure feedback is not the reaction to the reported CAD pressure force

**Confirmed in a manufactured pressure-impulse test. Directly relevant to GUI AoA/glide searches.**

`ExternalAeroCore` calls `apply_projected_pressure_gradient` on the collocated path. This computes area-weighted averages of face velocity corrections, then applies them to cell velocities (`amr_advection.cu:373`, `:385`, `:1372`). The reported force instead integrates pressure times CAD patch normals (`aero_loads.cpp:41`). Advective flux conservation does not establish equality between those two pressure operations.

With an artificial smooth pressure field supported strictly inside the outer domain, and zero initial cell momentum, the present implementation gives:

| Quantity | X force | Z force |
| --- | ---: | ---: |
| CAD pressure integration | `−42.0564 N` | `189.1171 N` |
| Negative fluid pressure impulse / dt | `−11.0098 N` | `102.3333 N` |
| Difference | `31.0466 N` | `−86.7838 N` |

No outer-boundary pressure traction is present in that fixture. These numbers are an algebraic operator check, not an aerodynamic prediction. They demonstrate that the force returned to a future structural solver need not equal the fluid's reaction.

The existing transient momentum diagnostic also flags a large discrepancy in the short +12° collocated run: about **2,317 N in X and 53 N in Z** at step 100. That diagnostic samples over a finite time interval and omits a complete viscous boundary/time-quadrature treatment, so use it as supporting evidence; the isolated pressure-impulse test is the cleaner finding.

**Repair direction:** derive cell momentum pressure force, face mass-flux correction, and CAD reaction as a compatible set. Include the actual boundary terms and both pressure/viscous reaction in conservation tests. Geometric closure must be repaired alongside this. Do not promote a solver to “load solver” based only on its conservative advection.

### R5 — Cropped solid cases retain fluid cells inside the solid

**Confirmed against OCCT classification; relevant to the NACA validation and potentially half-wing cuts.**

`solid_embedded_boundary.cpp:1022` seeds every regular cell on the atlas boundary as exterior fluid. A cropped-Y airfoil has solid cross-section on that boundary. Seeding it floods the airfoil interior as fluid. Separately, `amr_eb.cpp:97` skips levels without surface-intersecting active bricks, while `amr_pressure.cpp:254` initially marks every active brick cell as fluid; solid classification also needs to cover non-surface bricks.

The zero-incidence NACA slab diagnostic finds **162,560 active sealed states**, representing `0.00968933 m³`. All **32 sampled sealed-state centroids** classify as **inside** the authoritative OCCT solid. The +4° NACA run independently reports 162,528 sealed DOFs. The ordinary full PlanB case has zero sealed states.

The NACA test still explicitly requires `gauges == 1` (`naca_validation_probe.cpp:711` and `:939`). That requirement belongs to the old enclosed-fluid/fabric model; it accidentally rewards this retained interior. An ordinary solid airfoil in an outlet-connected external domain should not require a fluid pressure gauge inside its material.

The current sampled NACA result reports zero lift contribution from perturbing sealed pressure, so this review does **not** claim the interior alone caused its low external lift. It is nevertheless a geometry-model and validation defect, adds unnecessary unknowns, and complicates sampling and boundary behaviour.

**Repair direction:** classify boundary seeds before flooding whenever a crop/symmetry plane intersects the solid; classify whole interior bricks at every active level. Assert that ordinary fluid sample points are exterior and that cropped/full topology agrees where it should. Replace the inherited sealed-interior acceptance rule with the intended solid-domain invariant.

### R6 — Validation and GUI routing do not certify the same model

**Confirmed from source and freshly run tests.**

All four CTests pass in **0.38 s**, but:

- `parity` validates the retained uniform FP64 kernels; its own header says production is different.
- `paraglider_gpu_probe.cpp` creates an empty 64³ EB grid, benchmarks `DeviceEbPressureOperator`, and accepts a positive last output. It never executes the current composite wing timestep.
- `aero_convergence` checks decision/averaging helpers.
- `visualization_fields` checks display sampling/field logic.

The live CTest suite therefore cannot detect R1–R5. Many elaborate manufactured and physical gates claimed in the historical documents no longer have corresponding sources/targets in this checkout. Old executables left in build folders are not reproducible evidence.

The NACA configuration uses projection tolerance `5e-4`, whereas `mechanics_ok` requires residual `≤1.1e-5`. A run can satisfy its solver configuration and fail the benchmark's projection criterion. The NACA probe also defaults to MAC/slip while the GUI automatically forces collocated/shear for sweeps. The current coefficient envelope is `max(0.02, 5% |target|)`, not the historical 15% envelope.

The diagnostic NACA +4° run reached only `t=0.1440 s` of the requested `0.8 s` before its explicit 180-second limit. At the latest force sample, pressure `CL=0.46368`; final circulation-derived `CL=0.37869`; the stored XFOIL target is `0.7414`. It correctly exited **FAIL / CONTROLLED TIMEOUT**. This is **not** a completed coefficient comparison or evidence of a converged 37% error. Its final residual `4.95e-4` and normalized local flux defect `0.00366` also miss the declared mechanics thresholds.

**Repair direction:** make fresh tests execute the exact geometry and solver selected by the product. Preserve benchmark geometry, executable/source/config fingerprints, boundary conditions, wall/SGS model, complete histories, and stopping reason. Align the requested solve tolerance with the declared gate before running; demonstrate force insensitivity rather than simply tightening every number.

### R7 — Near-wall and AMR transport has additional consistency/stability risks

**Confirmed implementation properties; their individual contribution to PlanB lift remains unmeasured.**

The ordinary MAC advection path is point-sampled semi-Lagrangian transport with local bounded correction (`amr_advection.cu:658`, `:662`, `:1715`), not a global conservative momentum-flux update. Near the surface it changes to `side_safe_muscl`, with incomplete-stencil fallbacks. The separate conservative-uniform methods and historical momentum-interface helpers are not what `ExternalAeroCore::step` calls. Equal mass flux across a pressure AMR interface is not momentum refluxing.

The MAC Smagorinsky diffusion kernel uses local viscosity times a velocity Laplacian, with molecular-only viscosity in its protected band (`amr_advection.cu:725` onward). That is not generally the divergence of a variable eddy-viscosity stress. It executes one explicit diffusion step without the diffusion-rate subcycling used by the collocated branch. The MAC timestep uses maxima of velocity components and a separate aperture rate; it does not bound the full multidimensional transport/diffusion operator. `side_safe_muscl` rescales increments when the summed directional Courant number exceeds one, altering the intended update instead of reducing the actual timestep.

`reconstruct_compatible_carrier_kernel` also reads a destination value and then atomically adds a correction while multiple carrier records may share the same destination. The constructor explicitly permits multiplicity two (`amr_pressure.cu:1879`). Atomic addition alone does not make the preceding read part of an atomic transaction. This is a potential schedule-dependent update; repeat trials showed only small aperture differences, so the magnitude of this particular issue is not established.

The collocated branch uses piecewise-constant group redistribution of retained states, including after pressure and wall operations. That can flatten supported gradients. In its wall path, forces are finalized after redistribution (`amr_advection.cu:1355`), so the claimed fluid/shear reaction should be rechecked against the actual composed update, not only an isolated wall kernel.

**Repair direction:** specify the discrete transport and stress operators, make physical face ownership explicit, and test uniform/affine/vortex transport, energy behaviour, timestep convergence, and full-domain momentum balance. Use a single-writer gather or frozen-input delta buffer for shared carriers. Small-cell stabilization must have an accuracy argument; [weighted state redistribution research](https://arxiv.org/abs/2112.12360) is a relevant starting point.

### R8 — A statistically flat force is being asked to stand in for validated aerodynamics

**Confirmed limitation of the current stopping policy and configurations.**

`paraglider_sim_worker.cpp:959` calls the state conservative if volume-weighted RMS divergence is below an absolute `1e-3 s⁻¹`. It then allows adjacent mean-force windows to trigger a converged label after one domain flow-through. Three successive checks reuse heavily overlapping windows; they are not three independent statistical samples. Moment convergence, local mass defects, momentum/traction balance, circulation preservation, mesh error, and domain error are absent from that decision.

The two-level PlanB spacing is `0.125 m`; the GUI sweep's minimum is only 32 measured cells/chord (`paraglider_window.cpp:793`). That is a resource/resolution floor, not evidence that suction peaks, trailing-edge flow, tip vortices, or wall stress are resolved. The wake refinement box is centred on the wing bbox with a fixed radius; check that it contains the full spanwise wake and tips as resolution increases (`amr_grid.cpp:116`).

The Y/Z outer boundaries set normal velocity to zero (`amr_fields.cu:113` and `:120`): these are slip walls/symmetry conditions, not open infinity. They can represent an approximate distant boundary, but their effect on downwash and induced drag requires domain expansion tests. Domain rounding to whole bricks changes actual clearances.

**Repair direction:** distinguish numerical solve convergence, time averaging, discretization convergence, and physical validation in both stored results and UI. Use dimensionalized error budgets, multiple non-overlapping force/moment windows where appropriate, and independent grid/time/domain studies. [NASA's grid-convergence tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html) explains why an iterative residual alone does not measure discretization error.

### R9 — The glide-search trial overwrites the active trim placement

**Confirmed by source tracing; not reproduced by driving the GUI during this review. Directly relevant to preserving the user's neutral trim.**

`paraglider_window.cpp:786` captures the current viewer placement as the search baseline. At `:850`, each trial writes a rotated baseline into `config_.placement` and the viewer. Neither normal completion (`:908` onward) nor cancellation (`:962` onward) restores that baseline. The next search captures the retained trial as its new baseline; an ordinary rebuild uses it (`:1382`), and Save serializes the current placement (`:1078`). Thus trials can accumulate rotations across searches or become the saved case, despite the message that the rigged trim itself has not been redefined.

The helper's coordinate transformation is reasonable for its declared problem of fixed attitude relative to gravity and varying descent direction. That problem changes the wing's angle relative to the flow. It does not implement the user's requirement to preserve the supplied aerodynamic AoA. No such search should be used to compensate for the missing neutral-trim lift.

**Repair direction:** store the immutable source/trim placement separately from trial and display transforms. Restore it on completion, cancellation, and failure, and make save/rebuild use the intended case. Add a regression covering repeated searches and save/reload. For this project acceptance case, retain the existing AoA and derive glide quantities from verified loads in a documented wind/gravity frame.

## Physical assumptions that need a defined scope

A rigid inflated outer surface can produce lift without modelling internal ribs or air pressure. It cannot predict inflation, collapse, aeroelastic twist, or the difference between the exported equilibrium shape and the actual loaded flight shape. Closing vents also changes the inlet/lip geometry; quantify that approximation after the basic external solver works. Reintroducing internal fabric now would add complexity without fixing the demonstrated operator defects.

The setup has +X relative flow, +Z lift, and positive pitch about +Y raises the upstream leading edge. The basic pressure sign and SI units look consistent. `compute_pressure_loads` uses `−p A n` for an outward solid normal. Reference area affects coefficients, not force; its zero default explains unavailable CL/CD, not missing newtons. No automatic “make lift equal weight” term belongs in the CFD.

The user identifies the existing placement/AoA as the neutral-line trim case. Preserve it throughout recovery. The ordinary solver specifies a scalar freestream speed along +X (`paraglider_config.h:13`, `external_aero_core.cpp:229` onward); it has no separately evolved descent or gravity attitude. Wind-aligned coordinates can represent a descending glider without prescribing a vertical freestream component. The relation of that wind frame to gravity must be documented when interpreting flight performance; a missing vertical inlet velocity is not, by itself, a lift defect.

If “10% sink” means vertical sink speed divided by horizontal airspeed is approximately 0.10, the intended steady glide has `gamma=atan(0.10)≈5.71°` and whole-aircraft `L/D≈10`. Its force balance is `L=W cos(gamma)≈0.995W` and `D=W sin(gamma)≈0.0995W`. The sink does not mean 10% less lift, and this calculation is not an instruction to add 5.71° to the supplied AoA. Compare the unchanged case at measured mass, airspeed, density and loaded shape, including line, pilot and harness drag for whole-aircraft L/D. The precise meaning and uncertainty of the approximate sink observation should accompany the final flight comparison.

The existing fixed-attitude helper correctly discloses constant-coefficient speed scaling and absent pitch-moment balance, but its trial-angle workflow is unsuitable for preserving a fixed aerodynamic AoA and has the state-retention defect in R9. A full trim or stability solver is a separate capability; it is not needed to establish that the present neutral-trim aerodynamic case produces credible lift.

At `10 m/s`, the model is a low-speed incompressible problem. A slip-wall benchmark should deliberately validate attached pressure lift, including circulation and trailing-edge behaviour. A Navier–Stokes solver does not normally need a manually imposed lift or circulation value, but its boundary treatment and resolution must establish the appropriate trailing-edge flow. XFOIL uses a panel formulation with a Kutta condition; matching an isolated number from it is not equivalent to validating viscous separation. See the [XFOIL documentation](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt).

Spalding wall shear plus Smagorinsky SGS is not, by itself, a validated high-Reynolds-number airfoil transition/stall model. Choose and validate a wall-resolved, wall-modelled LES, or RANS strategy with its required grid and boundary conditions. Separate attached lift, profile drag, induced drag, and stall claims. The [Turbulence Modeling Resource NACA0012 case](https://tmbwg.github.io/turbmodels/naca0012_val.html) provides a defined viscous validation problem; use its stated conditions rather than mixing datasets.

## Reproduced runs and limits

| Run | Physical endpoint | Result | Interpretation |
| --- | ---: | --- | --- |
| All four current Release CTests | — | 4/4 pass | Existing narrow checks pass; production physics remains untested by them |
| PlanB MAC, unchanged neutral trim, 10 m/s, h=0.125 m | `1.60271 s`, 580 steps | Total `[D,L]=[79.924, −7.479] N` | Reproduces weak lift at the required placement; instantaneous endpoint, not converged |
| Earlier PlanB MAC diagnostic, +12°, same grid and speed | `1.60217 s`, 536 steps | Total `[D,L]=[144.362, 192.705] N` | Historical diagnostic only; not the user's trim case or an accepted fix |
| Earlier PlanB collocated diagnostic, +12°, same grid and speed | `0.084996 s`, 100 steps | Total `[D,L]=[71.576, 440.478] N`, preceding cell maximum `59.1 m/s` | Historical short diagnostic; different elapsed time from MAC, not the acceptance case |
| NACA2412 MAC/slip, +4°, 256 cells/chord | `0.144045 s`, 683 steps | Controlled timeout; flow gate fails | Requested 0.8 s benchmark was not completed |
| Independent operator audit | Initial and evolved PlanB states | R1–R4 reproduced | Direct numerical evidence, independent of a target flight weight |
| NACA solid-mask audit, 0° slab | Preprocessing only | 32/32 sampled sealed centroids inside solid | Confirms R5 |

Full source audit focus: STEP/placement/BVH and closed-solid EB; AMR hierarchy, fields, sampling and exchange; composite projection and both momentum paths; wall/shear and force integration; convergence/sweep/glide helpers and GUI routing; current probe sources and CMake test wiring. Retained uniform kernels and visualization code were checked for relevance and test scope, not re-proven line by line. This is not an exhaustive proof that every defect is listed. No full current FP64 wing polar, completed high-resolution NACA suite, independent external CFD comparison, or matched real-flight validation was performed.

## Recommended direction

Retain the closed-envelope scope, importer, viewer, and useful GPU infrastructure. Put new UI/FSI features behind completion of the numerical contracts. First make geometry close, pressure derivatives reproduce linear fields, transport approach identity at zero time, and CAD loads equal fluid reaction. Validate those changes through an identical small-airfoil test path, then return to the unchanged neutral-trim PlanB case and its measured glide performance.

Keep an established CFD implementation as an independent reference and a possible product backend. If a coherent custom discretization cannot satisfy these small gates, integrating that backend is a practical completion route. Maintaining two differently defective solvers and choosing whichever produces more lift is not a validation strategy.
