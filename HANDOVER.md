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
- conservative 2:1 velocity-state synchronization: transported fine MAC faces feed compact flux tiles, and projected tiles scatter back to fine faces plus an aperture-weighted coarse face;
- compact same-side EB aperture transport using structured carrier states, bounded minmod-limited MUSCL updates on complete same-axis chains, a first-order fallback at endings/junctions, molecular diffusion, and a graph-normal Smagorinsky estimate without a full-domain sparse velocity graph; the combined update obeys its local same-side stencil extrema without imposing an absolute velocity cap;
- adaptive one-global-step CFL control from a persistent GPU maximum reduction over regular AMR velocity fields; compact EB states retain their own bounded local update;
- `ExternalAeroCore`, which owns persistent device fields and runs advection -> LES/diffusion -> external BC -> compact EB transport -> coarse/fine synchronization -> composite projection without bulk per-step transfers;
- composite surface-patch mappings to global p+/p- pressure states and pressure-only triangle/whole-wing load publication;
- deterministic dynamic normal/parallel/inclined plate gates, developed opened/closed cavity flux validation, and an equal-finest AMR-versus-uniform inclined-plate comparison;
- two-sided pressure/Cp/pressure-force accumulation with winding-invariant force;
- BVH tracer collision and triangle-native delta-Cp render storage;
- a dedicated GUI worker that advances `ExternalAeroCore`, publishes pressure timing/residuals and pressure-only forces, and colors STEP triangles by visible-side Cp, Cp+, Cp-, or live delta-Cp;
- a purpose-built compact grid panel whose domain/AMR/solver/reference controls feed the actual `ParagliderConfig`, with bbox span/chord inference, an explicit LE/TE polarity flip, separate `Fit Wing`/`Fit Domain` camera framing, restored arrow/tracer tuning controls, wheel-safe editors, a finest-active-brick scalar slice sized per axis from finest cell spacing, and live CFL, regular/EB peak velocity, divergence, flux-error, and persistent-memory diagnostics.

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
- about 494.4 MiB total persistent estimate for fields, projection, compact EB transport and its two directed chain-link arrays, two advection states/masks, and locator.

`paraglider_case_probe` records long imported-wing histories in the exact frame used by the GUI. As the pre-MUSCL long-run baseline, with the confirmed leading edge facing upstream, default 0.25 same-side fragment merge, `1e-4 h^2` aperture cutoff, limited fabric-band link transport, graph-normal EB Smagorinsky treatment, adaptive CFL, and 1e-5 projection tolerance, a 500-step PlanB run reaches 0.531 s. The regular peak is 11.99 m/s, while the compact peak is still 46.8 m/s. Pressure-only force is still evolving at `[145.9, 0.10, 10.0]` N; max/RMS-volume divergence are `4.88e-4 / 3.93e-6 s^-1`, and net integrated flux error is `-3.15e-5 m^3/s`. The background mean is `[10.009, -1.3e-5, -2.5e-4]` m/s with zero reverse-flow volume. The last measured step is 73.0 ms, including 37.2 ms projection and 0.13 ms compact EB transport. Current directed compact-chain metadata brings persistent storage to about 494.4 MiB. This is a stability/conservation result, not a converged aerodynamic result; the compact peak remains an explicit blocker rather than being masked by a velocity guard. A new equal-duration long run is still required before quantifying the MUSCL force/stability change.

The case probe accepts `--max-levels N` for resolution studies. A four-level case builds 12.938 million pressure states with 1.427 GiB persistent storage. Without EB LES its compact aft-junction aperture reaches 102 m/s by 0.457 s while carrying conservative net fragment flux; with graph LES it settles near 49 m/s and reaches 0.294 s in 500 steps. The affected nodes lie 10-14 mm from two CAD faces and one is a 39-connection merged same-side fragment, identifying a complex trailing-edge/seam junction rather than a cross-fabric leak. At comparable physical time, three/four-level pressure forces remain roughly 12-14% apart. This demonstrates improved compact stability but incomplete force/grid convergence.

## Next engineering work

1. Extend compact MUSCL chains to multidimensional/junction-aware irregular reconstruction and extend graph-normal EB LES to a full irregular strain tensor/wall treatment.
2. Replace the now linearly consistent cross-level semi-Lagrangian momentum update with a conservative/refluxed finite-volume transport and run systematic grid/domain/orientation force-convergence studies.
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
