# ParaCFD handover

## Current state

ParaCFD is a purpose-built static-geometry paraglider solver. The former building/channel product path has been removed. The repository builds in production FP32 and validation FP64 modes. Its invariant is zero-thickness, two-sided fabric: never extrude it, flood-fill it, close openings, or replace it with a binary solid mask.

Implemented:

- OpenCascade STEP tessellation in metres with triangle-to-face provenance and optional UVs;
- affine-safe placement, static double-precision triangle BVH, face-only bbox/domain sizing (STEP wires do not affect the domain), and a shared FP32-canonical zero-origin frame that gives GUI/probes identical clipping coordinates;
- static balanced 2:1 brick AMR, pooled FP32 SoA fields, same-level halos, hash lookup, restriction/prolongation, and conservative four-tile coarse/fine interfaces;
- a sparse cross-brick EB atlas, so brick boundaries are not walls and covered coarse cells are excluded by finest-owner selection;
- smooth-sheet fitting plus a configurable finest-cell fluid-connectivity fallback supporting arbitrary fragment counts without solid/parity classification;
- split face apertures, including multiple disconnected openings on a partially covered aligned Cartesian face;
- conservative same-fluid small-fragment merging, a configurable/reportable sub-grid aperture cutoff, and explicitly retained/counted pressure-static pockets;
- a composite matrix-free pressure operator containing implicit regular faces plus compact EB and coarse/fine connections;
- CPU divergence/gradient reference operations and a persistent CUDA FP32/FP64 composite projection across all AMR levels;
- pressure gauges for every active fluid component disconnected from the outlet and a two-level additive Galerkin PCG preconditioner that retains compact nonlocal aggregate edges;
- bounded MacCormack/RK2 finest-brick GPU advection, minmod-limited linear transport on complete stencils of a static one-byte same-side six-link fabric graph (first-order fallback at incomplete stencils), and explicit molecular/Smagorinsky diffusion with continuous same-level cross-brick stencils and consistent old/forward states across levels;
- a separately callable one-level finite-volume MAC advection foundation plus CPU/GPU pairwise SoA component transport; regular face volume comes from half of each active adjacent pressure cell, normal tiles replace the covered-side share exactly, and tangential links use clipped staggered-dual overlaps; one unified graph owns 7,023 canonical faces and 7,840 transfers, including 501 multi-axis refinement-corner states, conserving each component to `1.47e-14`; its persistent GPU driver retains 6,144 exact pressure-tile weights, consumes GPU-resident projected mass flux, and scatters 384 coarse aliases with `5.70e-8 / 4.08e-8` FP32 error (zero FP64); it rejects alias ownership collisions and clears 7,920 compact-owned structured directions; a one-level boundary gate now uses prescribed +X inflow momentum, convective X-max export, zero Y/Z normal flux, and exact half-width boundary MAC volumes in both FP32 and FP64; compact EB and multilevel boundary coupling are not wired into the production conservative path yet;
- conservative 2:1 velocity-state synchronization: transported fine MAC faces feed compact flux tiles, and projected tiles scatter back to fine faces plus an aperture-weighted coarse face;
- compact same-side EB aperture transport using structured carrier states whose staggered volume is the sum of the adjacent composite fluid half-volumes (not a blanket `h^3`), bounded minmod-limited MUSCL updates on complete same-axis chains, a first-order fallback at endings/junctions, molecular diffusion, and component-collocated weighted least-squares reconstruction of all nine velocity derivatives for Smagorinsky strain without a full-domain sparse velocity graph; FP32 preprocessing discards modes beyond condition number `1e4`, and the combined update obeys its local same-side stencil extrema without imposing an absolute velocity cap;
- the earlier compact carrier-ring transport remains a negative/reference experiment: it conserves its isolated subset but has no free-stream-preserving mortar to the surrounding staggered dual volumes. Its replacement stores vector momentum directly on each actual pressure control volume (regular cell or split fragment), with projected MAC/aperture velocities supplying mass flux. Regular same-level faces stay implicit in per-level brick CUDA kernels; only 2:1 tiles plus the 110 perimeter/96 EB connections in the membrane fixture are compact. The one-level freestream, conservation, and analytical external-export gates pass. The three-level fixture advances 1,536 coarse/fine tiles, preserves freestream exactly, and closes perturbed-state conservation/parity to `1.42e-7 / 1.49e-8` FP32 and `4.44e-15 / 0` FP64. Physical-domain flags—not missing same-level neighbours—own external momentum boundaries. Conservative molecular diffusion uses the same topology and non-orthogonal pressure conductance; EB/multilevel energy decreases while momentum drift stays `7.45e-8 / 2.39e-9` FP32 and `2.84e-14 / 1.11e-14` FP64, and a face-aligned membrane cannot mix the two uniform sides. It is not called by `ExternalAeroCore` until reconstruction and LES are complete;
- cell-to-flux reconstruction for the conservative state: implicit regular faces plus compact regular/fragment, 2:1, and EB connections interpolate at the real face position between control-volume centroids; CPU/GPU error is `1.49e-8` FP32 / roundoff FP64, a face-aligned sheet remains exactly impermeable with opposite side states, and composite projection reduces the reconstructed multilevel divergence from `1.4` to `4.93e-5 s^-1` FP32 (`3.83e-11` FP64); the reverse post-projection pressure-impulse update remains deliberately unimplemented rather than smoothing the state through a face-average round trip;
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
- persistent Qt recent-file state: the File menu retains ten valid STEP/config paths, ordinary startup restores the most recent wing/config placement without building CFD, and missing paths are pruned.

Reusable FP64 uniform-MAC kernels remain only as CPU/GPU validation references. They are not a second application or result path.

## Current acceptance geometry

`Test-Data/PlanBParakite.step` is exercised through `configs/planb_parakite.json`. Its leading edge points toward source -Y, so the confirmed config applies -90 degrees about Z to map the leading edge to upstream -X, against ParaCFD's fixed +X freestream velocity.

At 2 mm tessellation, three AMR levels, 62.5 mm finest spacing, and `complex_subdivisions=4`, the case currently reports approximately:

- 49,673 triangles / 176 CAD faces;
- 120 active bricks / 3,932,160 active cells;
- 32,338 owned EB fragments, 114,042 face apertures, and 155,324 surface patches;
- zero unresolved cells and 348 pressure-static isolated pockets;
- 2,796 numerical aperture slivers discarded during the confirmed-upstream canonical-frame level-atlas construction, totalling 0.000283775 m^2 of cumulative atlas area at the configurable `1e-4 h^2` cutoff;
- 4,478,434 composite pressure slots with 81,920 coarse/fine and 102,227 EB connections;
- 124.14 MiB pooled FP32 field estimate;
- an initial +X freestream projection converging in about 189 PCG iterations to the tightened 1e-5 global relative residual in roughly 130-170 ms on the RTX 4090;
- 348.68 MiB estimated persistent GPU storage for fields, projection, and conservative 2:1 velocity synchronization (CUDA context/driver allocations excluded).
- 135,707 / 12,165,120 active MAC faces in the PlanB static fabric-protection band;
- about 25-30 ms bounded MacCormack advection, 8-10 ms LES/diffusion, 0.1-0.6 ms compact EB transport, and typically 40-140 ms projection as the warm solve evolves (individual timings vary);
- about 509.3 MiB total persistent estimate for fields, projection, compact EB transport with directed chain links plus least-squares matrices/RHS, fabric-wall work lists, two advection states/masks, and locator.

`paraglider_case_probe` records long imported-wing histories in the exact frame used by the GUI. On the current user-owned 47,997-triangle PlanB working fixture, with the confirmed leading edge facing upstream, default 0.25 same-side fragment merge, `1e-4 h^2` aperture cutoff, bounded fabric-band/MUSCL compact transport, component-collocated least-squares graph strain, actual composite MAC carrier volumes, adaptive CFL, and 1e-5 projection tolerance, 78,299 fragment/component stencils are well-conditioned/full rank and 23,105 wall-side control volumes are active. The production smooth-wall 1,500-step run reaches 3.831115 s. Regular/compact peaks are 12.35/15.57 m/s; pressure, viscous, and total forces are `[173.09, -4.73, -34.56]`, `[7.16, -0.024, 0.033]`, and `[180.25, -4.76, -34.53]` N. Max/RMS-volume divergence is `9.16e-4 / 5.18e-6 s^-1`, net flux error is `-6.31e-5 m^3/s`, and persistent storage is 521.11 MiB. The last step is 60.3 ms, including 17.7 ms projection and 0.34 ms compact/wall transport. The run remains stable beyond the former roughly 1,000-step collapse window without an absolute velocity guard. This is a stability/conservation gate, not a converged aerodynamic result. `Test-Data/PlanBParakite.step` contains user working-tree edits and is intentionally not included in solver commits.

The case probe accepts `--max-levels N` for resolution studies. A four-level case builds about 12.94 million pressure states and uses about 1.41 GiB persistent storage. Without EB LES its compact aft-junction aperture reached 102 m/s while carrying conservative net fragment flux; graph LES reduced that historical peak. The affected nodes lie 10-14 mm from two CAD faces and one is a 39-connection merged same-side fragment, identifying a complex trailing-edge/seam junction rather than a cross-fabric leak. Older fixed-step samples reached unequal physical times and must not be used as convergence evidence.

Use `paraglider_case_probe --physical-time T --max-levels N` for new comparisons. It stops after crossing the requested physical time, emits a machine-readable `[paraglider-case-summary]`, and treats an explicitly supplied `--steps` as a safety ceiling. `--thin-y-fraction 0.5 --thin-y-width W` reproduces the GUI's cropped centre-span diagnostic; `W` is metres and defaults to 0.125. The two-cell PlanB diagnostic remained stable through 8,000 UI steps. Its numerical probe ran 1,000 steps to 3.0 s with mean streamwise speed 10.006 m/s, 0.045% reversed volume, regular/compact peaks 12.52/15.55 m/s, maximum divergence `5.0e-4 s^-1`, and about 25.1 MiB persistent storage. At `t~=0.100 s` on the current working fixture, 2/3/4 levels (finest 0.125/0.0625/0.03125 m) give pressure-only `[Fx,Fz]` of `[180.2,-8.17]`, `[200.6,-5.25]`, and `[203.3,-24.77]` N. The 3-to-4-level streamwise change is only 1.4%, but the vertical component is not converged and this early time is still transient. Do not use the apparent Fx agreement as a force-validation claim.

## Next engineering work

1. Reconstruct the new multilevel cell/fragment-centred momentum state from corrected face fluxes and add same-side LES; switch `ExternalAeroCore` only after the composed constant/momentum/boundary gates pass, then rerun the equal-time convergence gate.
2. Run equal-time PlanB grid/domain studies with the production smooth-wall flux, validate total force convergence, and calibrate/replace the smooth-wall assumption if real fabric roughness matters.
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
