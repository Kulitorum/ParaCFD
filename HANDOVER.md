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
- a staged internal compact-momentum operator reconstructs fragment vectors from aperture/carrier half masses and applies equal-and-opposite embedded-aperture donor transfers (`1.19e-7` FP32, zero FP64 global physical-carrier momentum error); it is not called by `ExternalAeroCore`; a deterministic regular/compact ownership builder emits canonical MAC address, pressure endpoints, area, actual carrier volume, and endpoint flags for every ordinary incident face; the persistent projection bridge uploads those records and combines actual pooled MAC values with independent aperture states on the GPU—110 records/56 true perimeter faces close a parallel-flow fragment fixture exactly in FP32/FP64, while a one-face perturbation produces its analytical `0.5 m^3/s` imbalance; a separate staged driver applies both aperture and perimeter momentum transfers over unique physical carriers (`1.80e-6` FP32, zero FP64 conservation error), preserves parallel compact aperture states, and refreshes same-level brick aliases exactly; the surrounding conservative structured update does not consume the regular-side flux register yet;
- a separate compact two-sided fabric-wall work list sourced from every surface patch, including face-aligned sheets with no compact aperture; it applies an implicit molecular no-slip flux using actual patch area, volume, and normal centroid distance, then updates incident regular carriers/compact apertures without connecting the two sides;
- adaptive one-global-step CFL control from a persistent GPU maximum reduction over regular AMR velocity fields; compact EB states retain their own bounded local update;
- `ExternalAeroCore`, which owns persistent device fields and runs advection -> LES/diffusion -> external BC -> compact EB transport -> coarse/fine synchronization -> composite projection without bulk per-step transfers;
- composite surface-patch mappings to global p+/p- pressure states and pressure-only triangle/whole-wing load publication;
- deterministic dynamic normal/parallel/inclined plate gates, developed opened/closed cavity flux validation, and an equal-finest AMR-versus-uniform inclined-plate comparison;
- two-sided pressure/Cp/pressure-force accumulation with winding-invariant force;
- BVH tracer collision and triangle-native delta-Cp render storage;
- a dedicated GUI worker that advances `ExternalAeroCore`, publishes pressure timing/residuals and pressure-only forces, and colors STEP triangles by visible-side Cp, Cp+, Cp-, or live delta-Cp;
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

`paraglider_case_probe` records long imported-wing histories in the exact frame used by the GUI. On the current user-owned 47,997-triangle PlanB working fixture, with the confirmed leading edge facing upstream, default 0.25 same-side fragment merge, `1e-4 h^2` aperture cutoff, bounded fabric-band/MUSCL compact transport, component-collocated least-squares graph strain, molecular no-slip wall flux, actual composite MAC carrier volumes, adaptive CFL, and 1e-5 projection tolerance, 78,299 fragment/component stencils are well-conditioned/full rank and 23,105 wall-side control volumes are active. A 1,500-step run reaches 3.342 s. Regular/compact peaks are 12.33/17.17 m/s and the still-transient pressure-only force is `[168.25, -4.04, -31.83]` N. Max/RMS-volume divergence are `8.51e-4 / 5.06e-6 s^-1`, and net integrated flux error is `1.22e-5 m^3/s`. The last measured step is 64.4 ms, including 22.3 ms projection and 0.41 ms compact/wall transport; individual pressure iteration counts and timings vary. Persistent storage is 509.32 MiB. The run remains stable beyond the former roughly 1,000-step collapse window without an absolute velocity guard. This is a stability/conservation result, not a converged aerodynamic result. `Test-Data/PlanBParakite.step` contains user working-tree edits and is intentionally not included in solver commits.

The case probe accepts `--max-levels N` for resolution studies. A four-level case builds about 12.94 million pressure states and uses about 1.41 GiB persistent storage. Without EB LES its compact aft-junction aperture reached 102 m/s while carrying conservative net fragment flux; graph LES reduced that historical peak. The affected nodes lie 10-14 mm from two CAD faces and one is a 39-connection merged same-side fragment, identifying a complex trailing-edge/seam junction rather than a cross-fabric leak. Older fixed-step samples reached unequal physical times and must not be used as convergence evidence.

Use `paraglider_case_probe --physical-time T --max-levels N` for new comparisons. It stops after crossing the requested physical time, emits a machine-readable `[paraglider-case-summary]`, and treats an explicitly supplied `--steps` as a safety ceiling. `--thin-y-fraction 0.5 --thin-y-width W` reproduces the GUI's cropped centre-span diagnostic; `W` is metres and defaults to 0.125. The two-cell PlanB diagnostic remained stable through 8,000 UI steps. Its numerical probe ran 1,000 steps to 3.0 s with mean streamwise speed 10.006 m/s, 0.045% reversed volume, regular/compact peaks 12.52/15.55 m/s, maximum divergence `5.0e-4 s^-1`, and about 25.1 MiB persistent storage. At `t~=0.100 s` on the current working fixture, 2/3/4 levels (finest 0.125/0.0625/0.03125 m) give pressure-only `[Fx,Fz]` of `[180.2,-8.17]`, `[200.6,-5.25]`, and `[203.3,-24.77]` N. The 3-to-4-level streamwise change is only 1.4%, but the vertical component is not converged and this early time is still transient. Do not use the apparent Fx agreement as a force-validation claim.

## Next engineering work

1. Add exact multi-axis corner dual volumes and normal-tile GPU gather/scatter, then couple compact EB states into the now-masked pairwise/structured update; switch `ExternalAeroCore` only after composite constant/momentum gates pass, then rerun the equal-time convergence gate.
2. Add a high-Re smooth-fabric wall function and accumulate its viscous traction separately from pressure loads.
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
