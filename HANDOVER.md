# ParaCFD handover

## 2026-08-21 reset checkpoint -- read this first

This section supersedes the historical geometry and pressure-solver descriptions below.
Do not restart the old open-fabric/contact-topology work merely because it remains in the
history of this file.

### Current physical and geometry decision

ParaCFD is presently solving the external aerodynamics of an assumed rigid, inflated wing.
For that purpose the ribs, crossports, diagonals, and other internal fabric are irrelevant.
The CFD obstacle is the clean closed outer envelope exported by LeParaglider as
`Test-Data/PlanBParakite-Solid.step`. The vent openings are closed in this export, so OCCT can
provide a proper solid with an unambiguous inside, outside, and outward orientation. Do not
reintroduce the earlier zero-thickness internal-fabric pressure graph unless the user explicitly
changes the physical problem.

The current production geometry path is reported as `closed-solid-triangle-bsp`:

- OCCT loads and validates the closed solid (`69` faces, `4.40065444 m3`, no orientation reversal);
- the solid is tessellated once, placed in the CFD frame, and indexed by the triangle BVH;
- cell/triangle candidates use AABB queries and local polygon/BSP clipping rather than per-cell
  OCCT booleans or a 3-D point-sample flood fill;
- adjacent-cell triangle overlap pieces are coalesced into one conservative common region per
  fragment pair in `src/core/geometry/solid_embedded_boundary.cpp`;
- `src/core/fluid/amr_eb.cpp` now requires `options.closed_solid` and calls
  `build_closed_solid_embedded_boundary`;
- much of the old contact-conforming/open-fabric geometry implementation and its probes has been
  deliberately removed from the working tree.

Do not add isotropic remeshing yet. OCCT tessellation remains deflection-driven; the proposed
`huxingyi/isotropicremesher` integration and a nominal 50 mm target edge length were explicitly
deferred until the closed-solid CFD path works correctly.

The current two-level PlanB build reports approximately `168,938` triangles, `36` active bricks,
`1,313,933` pressure DOFs, `10,471` EB edges, `215,656` surface patches, and `464.71 MiB` persistent
GPU storage. STEP load is about `0.86--0.95 s`; all remaining preprocessing is about `2.2--2.4 s`.
Surface representation is essentially complete (`0.99999997` area coverage, no collapsed
two-sided patches). Shared-control-volume closure is not mathematically clean yet: the diagnostic
still reports `4,119` corrected rows, `sum|dA|=18.9290`, and `max|dA|=0.0270633`. Do not hide that
number; revisit common-face consistency if force accuracy points back to it.

### Current pressure solver decision

The former production pressure path was structurally unusable. Its displayed iteration count was
an outer flexible GCRO/GMRES count, and every outer basis vector launched a `27--50` iteration PCG
solve; every inner iteration included `16` coarse Jacobi sweeps and several synchronizing GPU
reductions. A nominal `330`-iteration step therefore meant thousands of full-grid passes and took
about `14--18 s`. Tuning that nested solver is not the current direction.

Production now uses one symmetric conservative orthogonal/two-point pressure operator and the
matching two-point regular/EB flux correction in `src/core/fluid/amr_pressure.cu`. Deferred
non-orthogonal flux kernels are skipped so the solved matrix and corrected divergence remain
consistent. The old full non-orthogonal GCRO path is retained only for diagnosis by setting:

```
PARACFD_PRESSURE_FULL_NONORTHOGONAL=1
```

Do not enable that flag for an ordinary GUI run: it restores the multi-second nested solver.

The direct FP32 solve has a real precision floor on this highly cut graph. Initialization stalled
near `2.9e-5` even after `5,000` PCG iterations, and an evolved RHS could stall around
`2.4e-4`; requesting `1e-5` is therefore dishonest in the current FP32 formulation. Explicit
defect correction is implemented, but it cannot remove residual components below the stored
FP32 operator/pressure resolution. `configs/planb_parakite.json` is consequently set to the
measured usable tolerance `5e-4`. Do not silently tighten it back to `1e-5` without changing the
numerics or precision and re-benchmarking the complete wing.

### Latest verified run

Command:

```
build-paraglider-ui/Release/paraglider_case_probe.exe \
  --config configs/planb_parakite.json --steps 100 --sample-every 10 \
  --staggered-momentum --projection-tolerance 5e-4
```

The run completed successfully in `12.2 s`, including STEP load and preprocessing. The first
evolved projection is a difficult initialization transient and took about `1.34 s` (`3,000`
accumulated PCG/refinement iterations). It did not recur. Normal early steps took `48--78 ms`
total with `28--59 ms` in projection; by step 100 the result was:

- physical time `0.261087 s` (this is still an early transient, not a converged aerodynamic run);
- regular/special velocity maxima `11.745 / 14.845 m/s`;
- total force `[109.749, 0.397, 17.645] N` and pressure force
  `[104.792, 0.355, 17.609] N`;
- volume-weighted RMS divergence `4.514e-5 1/s`, net integrated flux error
  `-1.564e-2 m3/s`;
- final pressure solve `65` iterations, residual `4.435e-4`;
- final step `49.27 ms`, of which projection was `29.41 ms`.

The initial compact velocity spike was about `70 m/s` and decayed rather than growing. Step-zero
pressure/load is an initialization impulse and must not be interpreted aerodynamically. This run
demonstrates boundedness and useful performance only; it does not yet prove lift convergence,
grid convergence, or validated non-orthogonal spatial accuracy.

All four configured Release tests pass: `parity`, `aero_convergence`, `paraglider_gpu`, and
`visualization_fields`. The rebuilt GUI is
`build-paraglider-ui/Release/paracfd-gui.exe`.

Build command used on this machine:

```
$env:CUDA_PATH='C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.1'
$env:CudaToolkitDir=$env:CUDA_PATH
cmake --build build-paraglider-ui --config Release --target paracfd-gui -- /m /p:TrackFileAccess=false
```

### Exact restart point

1. Run the rebuilt GUI with `configs/planb_parakite.json` and confirm its live timings match the
   probe rather than the obsolete `18 s/step` path.
2. Run long enough to cover physically meaningful flow-through time before judging whether the
   wing produces the expected force. One hundred steps covered only `0.261 s`.
3. If accuracy is wrong, inspect the remaining shared-control-volume closure corrections first.
   Do not return to internal ribs, vent connectivity, point-sampling boxes, or OCCT per-cell
   booleans; those solve the discarded physical model.
4. Solver improvements should target a genuine multigrid/direct SPD implementation and/or a
   precision strategy. Do not optimize or restore nested GCRO -> PCG.
5. Keep the conservative orthogonal flux/operator pairing intact. Any later deferred
   non-orthogonal correction needs a small fixed number of outer defect/Picard updates with a
   measured cost and stability gate, not hundreds of nested Krylov solves.

The working tree contains extensive intentional geometry deletions/changes plus the untracked
`Test-Data/PlanBParakite-Solid.step`. Preserve them; do not reset or clean the tree.

## Current state

ParaCFD is a purpose-built static-geometry paraglider solver. The former building/channel product path has been removed. The repository builds in production FP32 and validation FP64 modes. Its invariant is zero-thickness, two-sided fabric: never extrude it, flood-fill it, close openings, or replace it with a binary solid mask.

STEP import now has an explicit fast/checked split. Normal GUI loading uses OCCT's ordinary
per-face tessellation, marks contact connectivity as not checked, and permits only a prominently
labelled approximate CFD run. The preview BVH may reconcile differently tessellated pieces of the
same known CAD edge; it never welds distinct CAD entities by distance. `Perform detailed geometry
check (slower)` runs the exact BRep contact audit, atomizes partial contacts/junctions, creates a
contact-conforming triangle mesh with shared topology IDs, and then invokes the same production
finite-triangle EB builder. A detailed failure is fail-closed and diagnostic. The per-cell OCCT
General Fuse decomposer is an offline development oracle only and is not a second GUI solver.

The validated/default momentum solver is the face-centred staggered MAC path.
The collocated conservative-control-volume path described in the historical section
below is explicitly experimental because it becomes unstable on high-incidence
PlanB cases. Exact BVH same-side preprocessing is retained. Replacing the actual
O(N²) coarse/fine-link deduplication with a hash reduced PlanB preprocessing from
roughly 75--85 seconds to about 9 seconds without weakening geometry tests.

### Archived collocated imported-wing experiment (not production)

This experiment previously made `ExternalAeroCore` default to the conservative control-volume state and labelled it as such in the GUI. That promotion was reversed after high-incidence imported-wing failures. The face-centred MAC path is now the default and the checkbox opts into this experiment. The experiment's numerical construction and manufactured evidence are retained below to support future diagnosis; they are not current product claims.

The new timestep is conservative control-volume advection -> same-side molecular/Smagorinsky diffusion -> two-sided Spalding wall shear -> pressure-compatible provisional flux reconstruction -> external BC -> composite projection -> projection-consistent pressure feedback. Wall shear freezes each patch's nonlinear Spalding coefficient at the old velocity, assembles the symmetric tangential drag tensor on each affected control volume, and solves its three-component backward-Euler update. Shared wall control volumes therefore cannot overshoot or reverse velocity when many CAD patches contribute, while the per-patch force evaluated at the updated velocity remains exactly opposite the fluid impulse. The reconstruction interpolates the advanced cell state, removes the cell-averaged correction from the preceding projection, and restores that projection's exact face-normal correction. This preserves a projected state exactly while discarding unrelated face-only aperture modes that cannot be represented by one collocated vector. The zero-time projection is retained because it creates a divergence-free impermeable initial velocity, but its pressure is an initialization impulse rather than an evolved aerodynamic load: the GUI now withholds Cp and loads until the first complete transport/LES/wall/projection step, and the case probe labels its step-zero output explicitly. CFL uses the exact maximum `sum(outward volumetric flux)/volume` over every active control volume. No absolute speed cap, runaway guard, or timestep band-aid is used.

The current 47,997-triangle PlanB fixture exposed the explicit wall update by alternating its leading-edge pressure field every few steps and eventually exciting the old full-wing collapse. Disabling wall shear removed the mode, while disabling LES did not. With implicit wall coupling, full-wing production runs remain bounded through 1,500 steps at both -6 and -3 degrees AoA. At -6 degrees the regular/compact maxima are `12.34/12.28 m/s`, delta-Cp range is `1.64`, max/RMS divergence is `1.47e-4 / 3.83e-6 s^-1`, and total/pressure/viscous `[Fx,Fz]` is `[104.61,-85.46] / [97.43,-85.72] / [7.18,0.26] N`. At -3 degrees the corresponding maxima are `12.16/12.24 m/s`, delta-Cp range `1.43`, and total vertical force `-56.61 N`. The normalized area-weighted frame-to-frame delta-Cp RMS decays to `0.0040` and `0.0021`. Wall matrices are pooled only for the roughly 23,000 unique wall-adjacent CVs, not all 4.48 million pressure slots, so persistent FP32 storage is 654.09 MiB in the -6-degree case. These forces are transient and grid-dependent, not validated aerodynamic results.

### Analytic airfoil validation

`naca_step_generator` creates any valid four-digit NACA section as a finite-span OpenCascade STEP wing and can emit the exact matching XFOIL coordinate file. The 0012/2412/4415/4112 suite returns through STEP import, placement, tessellation, BVH, AMR, and two-sided EB preprocessing. Each 0.125 m free-slip span slab requires exact represented surface area, one disconnected sealed-interior pressure gauge, no inactive patch mapping, and a converged composite projection. The accepted production coefficient gate uses a physical 0.15 m finest collar, 256 cells/chord, 0.8 s at 10 m/s, a final 0.2 s mean, slip fabric, and official XFOIL 6.99 inviscid lift. Its declared FP32 envelope is `max(0.03, 15% |CL_XFOIL|)`, plus independent topology, RMS-divergence, integrated-flux, unsteadiness, and pressure/circulation gates. This validates attached pressure lift; it does not claim transition, separated-flow, or drag-polar accuracy. Natural- and forced-transition XFOIL columns and an exploratory smooth-wall result remain visible diagnostics rather than a target selected after a run. The exact 4112 +4-degree cusp uses a topology-valid half-cell phase for the transient; an h/512 topology-only preprocessing regression independently resolves phase zero to one sealed interior component.

The completed +4-degree production suite passes 4/4. ParaCFD/XFOIL-inviscid lift coefficients and errors are: NACA0012 `0.50388/0.4828 (+4.37%)`, NACA2412 `0.70650/0.7414 (-4.71%)`, NACA4415 `0.93052/1.0273 (-9.42%)`, and NACA4112 `0.83188/0.9192 (-9.50%)`. Independent circulation coefficients are `0.49312/0.69069/0.92066/0.82247`. Normalized maximum/local flux defects are at most `1.92e-4` of a finest freestream face flux; normalized volume-weighted RMS divergence is at most `3.0e-8`. Persistent storage is 420--456 MiB and a production step is 68.6--112.6 ms on the RTX 4090. Do not reinterpret this as separated-flow or drag validation.

Those figures describe the earlier qualitative triangle-EB validation path. They are
not a current strict-CAD acceptance result. The new fail-closed local-arrangement path
currently rejects the tracked NACA 2412 case (four OCCT General Fuse warnings and eight
shared-face ownership mismatches), so `naca_geometry` is intentionally red. A separate
closed triangulated-cylinder regression now retains every manufactured microscopic
opening and eliminates all explicit roundoff apertures. Preserving analytic-plane
source-triangle lineage reduced wrong pressure-owner edges from 2,348 to 1,292, but
fitted analytic planes still disagree with exact piecewise traces in adjacent curved
cells. The correct next step is a piecewise trace-to-fragment face-region adapter;
neither an aperture-area cutoff nor a whole-face plus/minus mask is acceptable, because
the latter closes the tested `3.4911e-16 m^2` confirmed-free vent. Until that graph audit
passes, the GUI's default path remains explicitly labelled qualitative preview and no
NACA or paraglider load should be called CAD-certified.

Implemented:

- OpenCascade STEP tessellation in metres with triangle-to-face provenance and optional UVs;
- deterministic OpenCascade four-digit NACA STEP generation and a checked-in NACA 2412 geometry/topology regression;
- affine-safe placement, static double-precision triangle BVH, face-only bbox/domain sizing (STEP wires do not affect the domain), and a shared FP32-canonical zero-origin frame that gives GUI/probes identical clipping coordinates;
- static balanced 2:1 brick AMR, pooled FP32 SoA fields, same-level halos, hash lookup, restriction/prolongation, and conservative four-tile coarse/fine interfaces;
- a sparse cross-brick EB atlas, so brick boundaries are not walls and covered coarse cells are excluded by finest-owner selection;
- smooth-sheet fitting plus a configurable finest-cell fluid-connectivity fallback supporting arbitrary fragment counts without solid/parity classification;
- split face apertures, including multiple disconnected openings on a partially covered aligned Cartesian face;
- conservative same-fluid small-fragment merging, a configurable small-aperture reporting threshold that never changes topology, and explicitly retained/counted pressure-static pockets;
- a composite matrix-free pressure operator containing implicit regular faces plus compact EB and coarse/fine connections;
- CPU divergence/gradient reference operations and a persistent CUDA FP32/FP64 composite projection across all AMR levels;
- pressure gauges for every active fluid component disconnected from the outlet and a two-level additive Galerkin PCG preconditioner that retains compact nonlocal aggregate edges;
- the production bounded face-centred MacCormack/RK2 path with a static exact same-side fabric graph, plus compact same-side EB aperture transport and the composite pressure projection;
- a separately callable one-level finite-volume MAC advection foundation plus CPU/GPU pairwise SoA component transport; regular face volume comes from half of each active adjacent pressure cell, normal tiles replace the covered-side share exactly, and tangential links use clipped staggered-dual overlaps; one unified graph owns 7,023 canonical faces and 7,840 transfers, including 501 multi-axis refinement-corner states, conserving each component to `1.47e-14`; its persistent GPU driver retains 6,144 exact pressure-tile weights, consumes GPU-resident projected mass flux, and scatters 384 coarse aliases with `5.70e-8 / 4.08e-8` FP32 error (zero FP64); it rejects alias ownership collisions and clears 7,920 compact-owned structured directions; a one-level boundary gate now uses prescribed +X inflow momentum, convective X-max export, zero Y/Z normal flux, and exact half-width boundary MAC volumes in both FP32 and FP64; compact EB and multilevel boundary coupling are not wired into the production conservative path yet;
- conservative 2:1 velocity-state synchronization: transported fine MAC faces feed compact flux tiles, and projected tiles scatter back to fine faces plus an aperture-weighted coarse face;
- compact same-side EB aperture transport using structured carrier states whose staggered volume is the sum of the adjacent composite fluid half-volumes (not a blanket `h^3`), bounded minmod-limited MUSCL updates on complete same-axis chains, a first-order fallback at endings/junctions, molecular diffusion, and component-collocated weighted least-squares reconstruction of all nine velocity derivatives for Smagorinsky strain without a full-domain sparse velocity graph; FP32 preprocessing discards modes beyond condition number `1e4`, and the combined update obeys its local same-side stencil extrema without imposing an absolute velocity cap;
- the collocated vector-momentum replacement remains an explicit research experiment: it stores vector momentum on each pressure control volume and has extensive manufactured conservative-transport, diffusion, wall, reconstruction, and parity gates, but high-incidence imported-wing cell/face pressure coupling is not stable enough for production;
- experimental collocated cell-to-flux reconstruction and reverse pressure impulse retain deterministic manufactured coverage: implicit regular faces plus compact regular/fragment, 2:1, and EB connections interpolate at the real face position between control-volume centroids; CPU/GPU reconstruction error is `1.49e-8` FP32 / roundoff FP64, and the manufactured AMR pressure-impulse gate conserves to `3.74e-9` FP32 / `5.12e-15` FP64. These results do not promote that path over the production face-centred solver;
- a compact two-sided fabric-wall work list sourced from every surface patch, including face-aligned sheets with no compact aperture; production solves Spalding's smooth-wall law at the actual wall distance and arbitrary patch orientation, applies the requested shear impulse through exact component dual-volume weights, and retains the equal-and-opposite CAD patch traction; low/high-Re, dissipation, winding, and reaction conservation tests pass in FP32/FP64 without connecting the two fluid sides; the older molecular no-slip sink remains validation-only;
- adaptive one-global-step CFL control from a persistent GPU maximum reduction over regular AMR velocity fields; compact EB states retain their own bounded local update;
- `ExternalAeroCore`, which owns persistent device fields and runs advection -> LES/diffusion -> external BC -> compact EB transport -> coarse/fine synchronization -> composite projection without bulk per-step transfers;
- composite surface-patch mappings to global p+/p- pressure states plus separate pressure, smooth-wall viscous, and total triangle/whole-wing load publication;
- deterministic dynamic normal/parallel/inclined plate gates, developed opened/closed cavity flux validation, and an equal-finest AMR-versus-uniform inclined-plate comparison;
- two-sided pressure/Cp/pressure-force accumulation with winding-invariant force;
- BVH tracer collision and triangle-native delta-Cp render storage;
- a dedicated GUI worker that advances `ExternalAeroCore`, publishes pressure timing/residuals and separated pressure/skin-friction/total forces, and colors STEP triangles by visible-side Cp, Cp+, Cp-, or live delta-Cp;
- a purpose-built compact grid panel whose domain/AMR/solver/reference controls feed the actual `ParagliderConfig`, with bbox span/chord inference, an explicit LE/TE polarity flip, separate `Fit Wing`/`Fit Domain` camera framing, restored arrow/tracer tuning controls, wheel-safe editors, a finest-active-brick scalar slice sized per axis from finest cell spacing, and live CFL, regular/EB peak velocity, divergence, flux-error, and persistent-memory diagnostics;
- a selectable cropped-Y collapse diagnostic that clips open fabric at a span station and requested physical width, snaps upward to whole finest cells, runs the real uniform solver, switches the viewer to a Y-normal 2-D presentation, and overlays the cell/timestep count. EB atlas halos are clamped at physical boundaries so stabilization cannot merge into a phantom halo cell.
- a physical half-wing mode that clips the oriented +Y half without a generated cap, omits original fabric wholly coplanar with the mirror plane while retaining crossing skins, anchors Y-min at the centre as a free-slip mirror plane, and reconstructs whole-wing integrated forces/moments/coefficients while retaining honest half-wing triangle Cp. It is restricted to symmetric geometry, zero sideslip, and cases where suppressing antisymmetric wake modes is acceptable.
- persistent Qt recent-file state: the File menu retains ten valid STEP/config paths, ordinary startup restores the most recent wing/config placement without building CFD, and missing paths are pruned.

GUI auto-pause computes its force/field settle score as soon as ten sampled fields support two comparison halves, but cannot pause until 0.5 streamwise-domain flow-through and 90% of the rolling observation window have elapsed. It then requires three consecutive scores below one plus acceptable RMS divergence. The overlay shows the score during warm-up. On either automatic or manual pause the solver readout reports wall time since `Build CFD Grid + Run`, including preprocessing and initialization. This is deliberately based on physical convective time rather than a resolution-dependent CFD step count.

The run-state overlay now also carries step, physical CFD time, and live/frozen Build-to-run wall time. Ordinary runs and AoA-sweep cases share one bounded auto-pause policy: conservation-aware field settling, convergence of time-weighted mean force between adjacent windows, or a configurable maximum physical flow-through count. The last condition is explicitly labelled `MAX FLOW — NOT CONVERGED`; it retains a labelled full or partial final-window mean when possible and must never be described as convergence. The single-run load panel uses the same mean-selection rules and explicitly labels converged, bounded, partial, or instantaneous endpoint loads. Resuming Play resets the observation epoch, Step cannot queue while Play is active, and evolved projection nonconvergence is a hard solver error that also cancels an active sweep. A GUI AoA sweep uses the current CAD placement as 0 degrees and rebuilds the actual AMR/EB case at every pitch. Rows contain mean whole-wing D/L/L-over-D plus CD/CL when an explicit reference area is present, times, and exit reason. It is cancelable and retains the final case for inspection. Scripted `--aoa-sweep min max step` and `--exit-after-auto-pause` exist for controller runs. The completed current-fixture half-wing production sweep at `[-6,-3,0,3,6,9,12]` degrees gives whole-wing mean `[D,L]` in newtons of `[153.99,-327.40]`, `[128.66,-172.89]`, `[128.35,-14.03]`, `[139.73,134.33]`, `[150.97,294.46]`, `[167.16,441.24]`, and `[192.22,553.87]`. Every case exits as `SETTLED` or `MEAN CONVERGED` after 1.05--1.32 flow-throughs; none uses the maximum-flow ceiling. The monotonic sign change exercises the production solver but is not a grid-converged paraglider polar. Coefficients remain withheld because the config reference area is zero.

The intended future structural integration is loose quasi-steady coupling, not live moving-boundary FSI: XPBD settles geometry -> ParaCFD converges a static CFD state -> per-triangle two-sided pressures/tractions return to XPBD -> repeat. Avoiding repeated STEP serialization, warm-starting fields on nearby geometries, and safely reusing far-field AMR bricks are future performance milestones; EB topology intersected by moved fabric must still be rebuilt rather than reused speculatively.

Reusable FP64 uniform-MAC kernels remain only as CPU/GPU validation references. They are not a second application or result path.

## Current acceptance geometry

`Test-Data/PlanBParakite.step` is exercised through `configs/planb_parakite.json`. Its leading edge points toward source -Y, so the confirmed config applies -90 degrees about Z to map the leading edge to upstream -X, against ParaCFD's fixed +X freestream velocity.

At 2 mm tessellation, three AMR levels, 62.5 mm finest spacing, and `complex_subdivisions=4`, the case currently reports approximately:

- 49,673 triangles / 176 CAD faces;
- 120 active bricks / 3,932,160 active cells;
- 32,338 owned EB fragments, 114,042 face apertures, and 155,324 surface patches;
- zero unresolved cells and 348 pressure-static isolated pockets;
- small positive-area aperture pieces are retained and counted below the configurable `1e-4 h^2` reporting threshold; the former 2,796-piece / 0.000283775 m^2 deletion is no longer permitted because it was not topology- or subdivision-invariant;
- 4,478,434 composite pressure slots with 81,920 coarse/fine and 102,227 EB connections;
- 124.14 MiB pooled FP32 field estimate;
- an initial +X freestream projection converging in about 189 PCG iterations to the tightened 1e-5 global relative residual in roughly 130-170 ms on the RTX 4090;
- 348.68 MiB estimated persistent GPU storage for fields, projection, and conservative 2:1 velocity synchronization (CUDA context/driver allocations excluded).
- 135,707 / 12,165,120 active MAC faces in the PlanB static fabric-protection band;
- about 25-30 ms bounded MacCormack advection, 8-10 ms LES/diffusion, 0.1-0.6 ms compact EB transport, and typically 40-140 ms projection as the warm solve evolves (individual timings vary);
- about 509.3 MiB total persistent estimate for fields, projection, compact EB transport with directed chain links plus least-squares matrices/RHS, fabric-wall work lists, two advection states/masks, and locator.

`paraglider_case_probe` records long imported-wing histories in the exact frame used by the GUI. On the current user-owned 47,997-triangle PlanB working fixture, a full-wing +12-degree production face-MAC run remains bounded through 1,200 steps. Pressure-only `[D,L]` is `[184.05,537.08] N`; regular/compact peaks are `14.31/18.57 m/s`; max/RMS-volume divergence is `1.55e-3 / 5.25e-6 s^-1`; and net flux error is `-3.20e-6 m^3/s`. The last step is 98.0 ms including 54.3 ms projection, with 80 PCG iterations to `9.996e-6` and 518.75 MiB persistent storage. This is a stability/conservation gate, not a converged aerodynamic result. `Test-Data/PlanBParakite.step` contains user working-tree edits and is intentionally not included in solver commits.

The current half-wing PlanB path retains 24,108 clipped triangles after omitting the original centre-plane fabric, uses 60 bricks and about 2.240 million pressure states, builds the core in about 3.1 s, and holds about 259.3 MiB persistent GPU storage. The full/half manufactured normal-plate force gate is currently red at 0.82% versus its 0.5% limit. Initialization agrees to `7.7e-8` relative, while all 2,048 mapped CPU pressure rows and a manufactured RHS agree exactly; investigate the evolved GPU projection/state path rather than clipping or load reconstruction. Symmetry still suppresses antisymmetric wake modes and is not automatically valid for every geometry/flow state; do not claim that a half-domain proves agreement with an unconstrained full-wing transient.

The case probe accepts `--max-levels N` for resolution studies. A four-level case builds about 12.94 million pressure states and uses about 1.41 GiB persistent storage. Without EB LES its compact aft-junction aperture reached 102 m/s while carrying conservative net fragment flux; graph LES reduced that historical peak. The affected nodes lie 10-14 mm from two CAD faces and one is a 39-connection merged same-side fragment, identifying a complex trailing-edge/seam junction rather than a cross-fabric leak. Older fixed-step samples reached unequal physical times and must not be used as convergence evidence.

Use `paraglider_case_probe --physical-time T --max-levels N` for new comparisons. It stops after crossing the requested physical time, emits a machine-readable `[paraglider-case-summary]`, and treats an explicitly supplied `--steps` as a safety ceiling. `--thin-y-fraction 0.5 --thin-y-width W` reproduces the GUI's cropped centre-span diagnostic; `W` is metres and defaults to 0.125. The two-cell PlanB diagnostic remained stable through 8,000 UI steps. Its numerical probe ran 1,000 steps to 3.0 s with mean streamwise speed 10.006 m/s, 0.045% reversed volume, regular/compact peaks 12.52/15.55 m/s, maximum divergence `5.0e-4 s^-1`, and about 25.1 MiB persistent storage. At `t~=0.100 s` on the current working fixture, 2/3/4 levels (finest 0.125/0.0625/0.03125 m) give pressure-only `[Fx,Fz]` of `[180.2,-8.17]`, `[200.6,-5.25]`, and `[203.3,-24.77]` N. The 3-to-4-level streamwise change is only 1.4%, but the vertical component is not converged and this early time is still transient. Do not use the apparent Fx agreement as a force-validation claim.

## Next engineering work

1. Replace the impulsive one-pass freestream initialization with a coupled/converged conservative initialization and validate the early force history.
2. Retain the exact per-face BVH same-side tests; the hash deduplication has already made production preprocessing interactive. Profile further only if new geometry grows this cost materially.
3. Extend the two-level Galerkin preconditioner into a recursive V-cycle, tighten local conservation gates, and add aperture-aware EB reconstruction when fabric reaches a 2:1 interface.
4. Extend the opened-cavity flux test to internal pressure equilibration and resolved inlet/crossport cases.
5. Replace the finest-brick-aware but throttled host slice snapshot with direct CUDA/OpenGL field sampling and add pressure-force vectors.

## Important files

- Geometry: `src/core/geometry/tri_mesh.h`, `step_import.*`, `triangle_bvh.*`, `embedded_boundary.*`
- AMR fields/exchange: `src/core/fluid/amr_grid.*`, `amr_fields.*`, `amr_exchange.*`
- AMR EB/pressure: `src/core/fluid/amr_eb.*`, `amr_pressure.*`, `eb_pressure.*`
- Loads/config: `src/core/aero_loads.*`, `src/core/paraglider_config.*`
- GUI: `src/gui/paraglider_window.*`, `paraglider_sim_worker.*`, `slice_viewer.*`
- Tests/probes: `tools/paraglider_geometry_probe.cpp`, `paraglider_gpu_probe.cpp`, `paraglider_flow_probe.cpp`, `paraglider_probe.cpp`

## Validation policy

Do not bless changed output merely because an operator changed. Validate redesigned geometry and operators against analytical or CPU-double references. Production fields are FP32; geometry, residual accumulation where useful, and force totals remain FP64. Never hide unresolved topology or cross-fabric connections by loosening tolerances.
