Latest recovery results: [Final attempt — pressure traction, small-cell increments and independent checks](review/FINAL_ATTEMPT.md). The real-flight lift remains unresolved.

# ParaCFD completion plan — 6 September 2026

This plan supersedes the historical implementation checklist for future work. Its evidence is [SIMULATION_REVIEW.md](SIMULATION_REVIEW.md). “Finished” first means a reproducible, validated **rigid-wing external aerodynamic solver**. Real-flight trim and structural coupling are subsequent, separately accepted capabilities.

**Fixed requirement:** the user's PlanB wing is already at neutral line trim. Keep its supplied AoA and orientation unchanged. Recover credible lift and approximately 10% sink at that trim; do not search for an angle that makes a defective solver produce lift. If sink means vertical/horizontal airspeed ≈0.10, the comparison target is approximately 10:1 whole-aircraft glide. The benchmark airfoils below have their own prescribed incidences and do not redefine PlanB's trim.

## 1. Establish one reproducible case and one declared solver contract

Preserve the existing working tree and record the exact STEP, config, source/executable hashes, placement, chord/span/area conventions, and solver/wall options. Add a schema-versioned run manifest and machine-readable force/moment history with physical time. Persist execution options currently absent from JSON, especially momentum and wall modes. Route GUI single runs, sweeps, glide trials, and command-line validation through the same case/execution description.

Lock the neutral-trim source placement separately from temporary display/trial transforms. Fix the existing glide-search state retention (R9): completion, cancellation, failure, rebuild and save/reload must preserve the intended baseline. Do not run a PlanB angle search as part of this recovery. Record wind-frame versus gravity-frame conventions without introducing an AoA offset.

Make the current solver variants explicit research options until their gates pass. Remove automatic promotion of the collocated branch to an accepted “load solver”; preserve users' ability to compare it diagnostically. Update CLAUDE/HANDOVER/RESEARCH to one current architecture description, with historical claims clearly archived. Consolidate build targets so stale binaries are not mistaken for available tests.

**Exit evidence:** the GUI and CLI produce the same geometry fingerprint, topology counts, solver choice, and force history for the same small case. The manifest records actual padded domain bounds. The current failing diagnostics reproduce from a fresh build.

## 2. Repair geometry and enforce control-volume identities

Fix solid classification at crop/symmetry boundaries and for whole interior bricks. Replace the old required sealed-airfoil gauge with an exterior-solid invariant. Build fluid volumes, centroids, apertures, wall patches, and shared-face partitions from one consistent representation. Audit area-vector closure and first moments before and after merging; preserve failure provenance. Include the minimum retained volumes and agglomeration extent in case records.

Fixtures: axis-aligned and inclined boxes, curved cylinder, sharp airfoil trailing edge, full and cropped solids, half-wing symmetry, translated/rotated copies, and an object wholly enclosing a coarse brick. Use OCCT selectively as an independent material/volume check, avoiding per-cell production booleans.

**Exit evidence:** zero active ordinary fluid states inside the certified solid on these cases; every shared aperture has one consistent area/centroid and opposite orientation; no finite area-vector closure defects. Constant pressure gives zero resultant on a closed body; linear pressure gives the analytically expected integrated force. Grid-phase changes do not produce topology-dependent force jumps beyond the declared spatial error.

## 3. Make pressure, momentum, stabilization, and loads compatible

Write the discrete equations and ownership rules before changing kernels. Select one result-producing path. The pressure correction must be the derivative represented by the solved operator, and its momentum reaction must match the reported CAD traction. Repair the missing tangential/non-orthogonal terms and AMR interface derivatives. Use an SPD solver only for an SPD operator; keep the recursive preconditioner where appropriate.

Replace the timestep-independent MAC remapping with a derived consistent update, or repair the collocated path as the sole result solver. Do not combine these changes with arbitrary force calibration. Preserve useful conservative flux machinery, but require the composed pressure/transport/wall update to pass. Use single-writer output or frozen-input deltas for shared carriers. Give small-cell redistribution an affine-preservation/accuracy contract and include it in every reaction test.

Implement a consistent variable-viscosity stress and a timestep bound for the actual multidimensional transport/diffusion scheme. If explicit terms need subcycling, apply and report it on the selected production path. Preserve gradients/circulation under grid transfer and account for momentum across AMR interfaces.

**Exit evidence:**

- The new review zero-time test approaches identity on initial and evolved projected states as dt decreases, within storage roundoff.
- Constant and linear pressure derivatives are exact to the relevant FP32/FP64 tolerance, including every AMR orientation and cut-cell geometry.
- Uniform/affine velocity and a transported vortex pass on uniform and AMR grids; measure amplitude, circulation, energy, and momentum error.
- Pressure and wall-fluid impulse match CAD force and moment. A complete control-volume budget includes unsteady storage, all boundary fluxes/tractions, and correct temporal quadrature.
- Halving timestep produces the expected convergence behaviour; it does not strengthen an undocumented filter.
- CPU reference and GPU **production** operators agree on the same fixtures. Tests assert finite outputs and deliberately fail for disabled work/incorrect signs.

**Architecture decision:** review these contracts before investing in another large polar. If neither retained path can meet them with a clear discrete formulation, retain ParaCFD's geometry/UI and use an established finite-volume solver as the result backend. The same benchmark ladder below applies to either route.

## 4. Validate attached airfoil lift through the exact product path

Start with a small reproducible section domain and three angles on NACA0012 (negative, zero, positive). Add a cambered NACA2412 at zero and positive angle. Keep the geometry convention, wall/SGS settings, Reynolds number, and reference source fixed before running. For an inviscid pressure gate, check the attached trailing-edge flow and circulation; for viscous gates, use a defined laminar or turbulent benchmark with matched wall/transition conditions.

Run identical cases through an independent reference solver. Establish both surface pressure distribution and integrated CL; do not fit a lift multiplier. Use 2D circulation as a section diagnostic, not a substitute for 3D whole-wing force. Correct the current NACA gauge and residual-contract defects. Keep a small mandatory CTest plus a longer scheduled validation matrix with archived data.

**Proposed initial acceptance targets, to freeze before examining results:** `|CL−CLref| ≤ max(0.02, 5% |CLref|)` for the declared attached reference; correct symmetry/sign and camber zero-lift shift; no unexplained pressure/circulation discrepancy. Show at least three consistently refined meshes and two smaller timesteps, then increase domain clearances. Seek finest-pair lift changes below 2–3%, timestep effect below 1%, and domain effect below 1%; report estimated uncertainty where asymptotic convergence is not demonstrated. These are project targets, not universal CFD standards or evidence that a particular cells/chord count is sufficient.

Run for a measured transient-removal interval and sufficient complete averaging windows. A wall-time cap yields an incomplete result and fails validation. Set pressure/flux tolerances from measured force sensitivity; the current residual alone is not an acceptance certificate.

## 5. Validate the 3D rigid wing at unchanged neutral trim

Use a simple finite-span wing before the complex PlanB surface. Compare full-span/half-span force and moments, span loading, downwash, and induced drag against a matched independent computation. Prove the symmetry implementation through evolved flow, not just initialization. Extend wake refinement to tips and the full relevant wake.

Then run PlanB at the supplied neutral-trim AoA using the frozen loaded outer envelope. Compare full-wing and half-wing only where symmetry is physically appropriate. Repeat the grid/time/domain studies at that same trim; compute force/moment means with confidence/variability information. Establish tessellation and small-cell-merge sensitivity separately from fluid-grid sensitivity. Add proper open far-field conditions if distant slip walls are too expensive or influential. A future PlanB polar is separate work and is not a prerequisite or substitute for recovering this known flying trim case.

**Exit evidence:** versioned, independently checked PlanB loads at unchanged neutral trim, with CL/CD/Cm conventions, uncertainty, topology and pressure diagnostics, and reproducible resource/time estimates. No design-quality label is granted solely by the auto-pause decision. Quantitative glide requires validated drag; stall remains outside the initial claim.

## 6. Match flight evidence and complete the user workflow

Obtain measured all-up mass, airspeed and density/altitude, glide/sink data, trim/brake/riser setup, and the loaded geometry corresponding to the export. Preserve the user-confirmed trim/AoA and document its coordinate convention. Check `L=W cos(gamma)` and `D=W sin(gamma)` for the complete aircraft, including pilot, harness and line drag. For vertical/horizontal airspeed ≈0.10, compare against `gamma≈5.71°`, `L/D≈10`, and lift ≈99.5% of weight. Treat that as a flight-performance check, not a prescribed CFD force or extra pitch rotation. Rerun CFD near any scaled equilibrium speed at the same AoA, because Reynolds number and separation can change.

Provide a fixed-AoA flight estimate from the verified force direction and magnitude, with gravity-frame bookkeeping separate from the immutable aerodynamic case. The current fixed-attitude angle-search helper addresses a different problem and should not run implicitly. A full trim capability needs aerodynamic moments, mass/CG and riser/line/harness forces, plus their equilibria; flight stability is a further task. Reopen vent or aeroelastic modelling only if the defined external-envelope approximation cannot meet the measured objective.

Finish packaging, fresh-machine build/deployment checks, result export, cancellation/failure recovery, and clear result status: running, time-averaged, numerically verified, physically validated, or incomplete. Display assumptions and reference conventions with the results.

**Version 1 completion gate:** a saved case can be rerun from a clean build without changing trim, GUI and CLI agree, the selected solver passes the manufactured and airfoil/finite-wing gates, the unchanged PlanB trim case has a documented error budget and independent comparison, and the claimed flight prediction matches defined measurements within that combined uncertainty. FSI, inflation, collapse, and unrestricted stall prediction are not implicitly included.

## Repair progress

The latest [flight recovery follow-up](review/FLIGHT_RECOVERY.md) fixes affine velocity transfer at cut/AMR faces and insufficient velocity-gradient support. Full PlanB contracts and the inviscid airfoil benchmark pass, but the known flying/gliding state remains unvalidated. The export's horizon-versus-relative-airflow convention must be resolved before adding descent kinematics; otherwise the user-protected AoA could be changed inadvertently.

The current implementation and measured results are recorded in [review/REPAIR_STATUS.md](review/REPAIR_STATUS.md). Geometry closure, sealed-solid masking, nonorthogonal pressure corrections, zero-time transport, conservative pressure reaction, pressure-predictor consistency and trim restoration have been repaired. The conservative branch is now the default; MAC remains available for diagnostic comparison. The full PlanB operator checks and a NACA2412 inviscid benchmark pass, while the fixed horizontal-flow PlanB case still does not establish the requested flying/gliding state. Confirm the flight-frame metadata before introducing any sink-motion boundary condition; do not rotate the supplied trim to obtain lift.

The next physical mesh ladder should retain exactly the same placement while refining h from 0.125 to 0.0625 and 0.03125 m, with smaller bricks and a narrow surface collar to control memory. This is a proposed convergence study, not a claim that its finest mesh is sufficient.

## Immediate work queue

1. Preserve/run `review/physics_audit_probe.cpp` and turn its failures into focused production tests.
2. Repair solid-mask seeding and local geometric closure (R3/R5).
3. Repair pressure linear consistency and pressure/CAD reaction (R2/R4).
4. Replace the zero-time MAC filter or complete the selected collocated formulation; test the composed update (R1/R7).
5. Unify GUI/probe solver selection, preserve the immutable trim placement, and restore matching end-to-end tests (R6/R9).
6. Run the attached-airfoil → finite-wing → unchanged-trim PlanB validation ladder, then the matched approximately 10% sink flight case (R8 and physical scope).

Feature work, remeshing experiments, and performance tuning should follow the correctness gate they serve. Optimizing an inconsistent operator is not on the critical path to a trustworthy result.
