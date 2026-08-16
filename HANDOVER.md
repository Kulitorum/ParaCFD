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
- bounded MacCormack/RK2 finest-brick GPU advection, a static one-byte same-side six-link graph in the fabric band, and explicit molecular/Smagorinsky diffusion with continuous same-level cross-brick stencils and consistent old/forward states across levels;
- conservative 2:1 velocity-state synchronization: transported fine MAC faces feed compact flux tiles, and projected tiles scatter back to fine faces plus an aperture-weighted coarse face;
- compact same-side EB aperture transport using structured carrier states, first-order upwinding, and molecular diffusion without a full-domain sparse velocity graph;
- adaptive one-global-step CFL control from a persistent GPU maximum reduction over regular AMR velocity fields; compact EB states retain their own bounded local update;
- `ExternalAeroCore`, which owns persistent device fields and runs advection -> LES/diffusion -> external BC -> compact EB transport -> coarse/fine synchronization -> composite projection without bulk per-step transfers;
- composite surface-patch mappings to global p+/p- pressure states and pressure-only triangle/whole-wing load publication;
- deterministic dynamic normal/parallel/inclined plate gates, developed opened/closed cavity flux validation, and an equal-finest AMR-versus-uniform inclined-plate comparison;
- two-sided pressure/Cp/pressure-force accumulation with winding-invariant force;
- BVH tracer collision and triangle-native delta-Cp render storage;
- a dedicated GUI worker that advances `ExternalAeroCore`, publishes pressure timing/residuals and pressure-only forces, and colors STEP triangles by visible-side Cp, Cp+, Cp-, or live delta-Cp;
- a purpose-built grid panel whose domain/AMR/solver/reference controls feed the actual `ParagliderConfig`, with bbox span/chord inference, an explicit LE/TE polarity flip, separate `Fit Wing`/`Fit Domain` camera framing, and live CFL, regular/EB peak velocity, divergence, flux-error, and persistent-memory diagnostics.

Reusable FP64 uniform-MAC kernels remain only as CPU/GPU validation references. They are not a second application or result path.

## Current acceptance geometry

`Test-Data/PlanBParakite.step` is exercised through `configs/planb_parakite.json`. That exporter uses -Y as forward, so the config applies +90 degrees about Z to map the model to ParaCFD's fixed +X freestream.

At 2 mm tessellation, three AMR levels, 62.5 mm finest spacing, and `complex_subdivisions=4`, the case currently reports approximately:

- 49,673 triangles / 176 CAD faces;
- 120 active bricks / 3,932,160 active cells;
- 32,338 owned EB fragments, 114,042 face apertures, and 155,324 surface patches;
- zero unresolved cells and 348 pressure-static isolated pockets;
- 3,007 numerical aperture slivers discarded during canonical-frame level-atlas construction, totalling 0.000312367 m^2 of cumulative atlas area at the configurable `1e-4 h^2` cutoff;
- 4,478,434 composite pressure slots with 81,920 coarse/fine and 102,227 EB connections;
- 124.14 MiB pooled FP32 field estimate;
- an initial +X freestream projection converging in about 189 PCG iterations to the tightened 1e-5 global relative residual in roughly 130-170 ms on the RTX 4090;
- 348.68 MiB estimated persistent GPU storage for fields, projection, and conservative 2:1 velocity synchronization (CUDA context/driver allocations excluded).
- 135,577 / 12,165,120 active MAC faces in the PlanB static fabric-protection band;
- about 25-30 ms bounded MacCormack advection, 8-10 ms LES/diffusion, 0.1-0.6 ms compact EB transport, and typically 40-140 ms projection as the warm solve evolves (individual timings vary);
- 493.55 MiB total persistent estimate for fields, projection, compact EB transport, two advection states/masks, and locator.

`paraglider_case_probe` records long imported-wing histories in the exact frame used by the GUI. With the default 0.25 same-side fragment merge, `1e-4 h^2` aperture cutoff, precomputed fabric-band link graph, adaptive CFL, and 1e-5 projection tolerance, a 500-step PlanB run reached 0.622 s. The regular peak stayed near 12.1 m/s instead of the previous 56.1 m/s leading-edge runaway. A compact EB aperture state remained near 46.0 m/s. Pressure-only force was still evolving at `[132.0, -0.57, 38.1]` N; max/RMS-volume divergence were `6.10e-4 / 3.90e-6 s^-1`, and net integrated flux error was `-6.02e-5 m^3/s`. The last measured step was 70.1 ms, including 35.6 ms projection and 0.15 ms compact EB transport. Persistent storage is 493.55 MiB. This is a stability/conservation result, not a converged aerodynamic result.

The case probe accepts `--max-levels N` for resolution studies. The current coarse two-level case uses 148.85 MiB and reaches a compact peak of about 29.0 m/s after 300 steps at 1.362 s. The default three-level case uses 493.55 MiB and retains the 46.0 m/s compact peak above. A four-level case builds 12.938 million pressure states with 1.427 GiB persistent storage; after 500 steps at only 0.202 s its regular peak is 13.9 m/s but its compact EB peak has grown to 78.8 m/s. The last four-level step takes 203 ms, including 95 ms projection. Because these samples are at unequal physical times, their forces are not a convergence comparison; they demonstrate that compact EB transport is resolution-sensitive and remains the primary numerical blocker.

## Next engineering work

1. Replace first-order fabric-band and compact-graph updates with higher-order same-side transport and consistent EB Smagorinsky treatment; diagnose and bound the resolution-growing compact PlanB aperture peak with a local acceptance gate.
2. Make general cross-level interpolation consistent with the conservative normal 2:1 flux state and run systematic grid/domain/orientation force-convergence studies.
3. Extend the two-level Galerkin preconditioner into a recursive V-cycle, tighten local conservation gates, and add aperture-aware EB reconstruction when fabric reaches a 2:1 interface.
4. Extend the opened-cavity flux test to internal pressure equilibration and resolved inlet/crossport cases.
5. Replace throttled host visualization snapshots with direct AMR-aware CUDA/OpenGL field sampling and add pressure-force vectors.

## Important files

- Geometry: `src/core/geometry/tri_mesh.h`, `step_import.*`, `triangle_bvh.*`, `embedded_boundary.*`
- AMR fields/exchange: `src/core/fluid/amr_grid.*`, `amr_fields.*`, `amr_exchange.*`
- AMR EB/pressure: `src/core/fluid/amr_eb.*`, `amr_pressure.*`, `eb_pressure.*`
- Loads/config: `src/core/aero_loads.*`, `src/core/paraglider_config.*`
- GUI: `src/gui/paraglider_window.*`, `paraglider_sim_worker.*`, `slice_viewer.*`
- Tests/probes: `tools/paraglider_geometry_probe.cpp`, `paraglider_gpu_probe.cpp`, `paraglider_flow_probe.cpp`, `paraglider_probe.cpp`

## Validation policy

Do not bless changed output merely because an operator changed. Validate redesigned geometry and operators against analytical or CPU-double references. Production fields are FP32; geometry, residual accumulation where useful, and force totals remain FP64. Never hide unresolved topology or cross-fabric connections by loosening tolerances.
