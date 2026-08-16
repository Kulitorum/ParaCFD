# ParaCFD handover

## Current state

ParaCFD is being migrated from the building/channel solver to a purpose-built static-geometry paraglider solver. The repository builds in production FP32 and validation FP64 modes. Its invariant is zero-thickness, two-sided fabric: never extrude it, flood-fill it, close openings, or replace it with a binary solid mask.

Implemented:

- OpenCascade STEP tessellation in metres with triangle-to-face provenance and optional UVs;
- affine-safe placement, static double-precision triangle BVH, and face-only bbox/domain sizing (STEP wires do not affect the domain);
- static balanced 2:1 brick AMR, pooled FP32 SoA fields, same-level halos, hash lookup, restriction/prolongation, and conservative four-tile coarse/fine interfaces;
- a sparse cross-brick EB atlas, so brick boundaries are not walls and covered coarse cells are excluded by finest-owner selection;
- smooth-sheet fitting plus a configurable finest-cell fluid-connectivity fallback supporting arbitrary fragment counts without solid/parity classification;
- split face apertures, including multiple disconnected openings on a partially covered aligned Cartesian face;
- conservative same-fluid small-fragment merging and explicitly retained/counted pressure-static pockets;
- a composite matrix-free pressure operator containing implicit regular faces plus compact EB and coarse/fine connections;
- CPU divergence/gradient reference operations and a persistent CUDA FP32/FP64 composite projection across all AMR levels;
- pressure gauges for every active fluid component disconnected from the outlet and a two-level additive Galerkin PCG preconditioner that retains compact nonlocal aggregate edges;
- first-order finest-brick GPU backtrace advection, a static no-cross-fabric protection band, and explicit molecular/Smagorinsky diffusion with continuous same-level cross-brick stencils and a common old-time state across levels;
- `ExternalAeroCore`, which owns persistent device fields and runs advection -> LES/diffusion -> external BC -> composite projection without bulk per-step transfers;
- composite surface-patch mappings to global p+/p- pressure states and pressure-only triangle/whole-wing load publication;
- deterministic dynamic normal/parallel/inclined plate gates, opened/closed cavity pressure connectivity, and an equal-finest AMR-versus-uniform inclined-plate comparison;
- two-sided pressure/Cp/pressure-force accumulation with winding-invariant force;
- BVH tracer collision and triangle-native delta-Cp render storage;
- a dedicated GUI worker that advances `ExternalAeroCore`, publishes pressure timing/residuals and pressure-only forces, and colors the STEP triangles by live delta-Cp; loading a paraglider releases the unrelated legacy channel GPU core.

The old solver remains buildable only as a reference. It is not a paraglider result path.

## Current acceptance geometry

`Test-Data/PlanBParakite.step` is exercised through `configs/planb_parakite.json`. That exporter uses -Y as forward, so the config applies +90 degrees about Z to map the model to ParaCFD's fixed +X freestream.

At 2 mm tessellation, three AMR levels, 62.5 mm finest spacing, and `complex_subdivisions=4`, the case currently reports approximately:

- 49,673 triangles / 176 CAD faces;
- 120 active bricks / 3,932,160 active cells;
- 32,338 owned EB fragments, 117,040 face apertures, and 155,324 surface patches;
- zero unresolved cells and 294 pressure-static isolated pockets;
- 4,484,367 composite pressure slots with 81,920 coarse/fine and 112,563 EB connections;
- 124.14 MiB pooled FP32 field estimate;
- an initial +X freestream projection converging in 168 PCG iterations to a 9.37e-5 global relative residual in about 138 ms on the RTX 4090;
- 347.02 MiB estimated persistent GPU storage for fields plus projection (CUDA context/driver allocations excluded).
- 135,577 / 12,165,120 active MAC faces in the PlanB static fabric-protection band;
- about 5.7 ms advection, 7.9 ms LES/diffusion, and 168 ms projection (192 iterations) for the first post-initialization step after enabling cross-brick stencils;
- 425.74 MiB total persistent estimate for fields, projection, advection scratch/masks, and locator.

## Next engineering work

1. Replace the static first-order fabric protection fallback with a higher-order side-aware reconstruction and make the current interpolatory velocity treatment conservative at 2:1 interfaces.
2. Advect/diffuse compact EB and coarse/fine aperture velocity states and tighten the existing normal/parallel/inclined plate and AMR-versus-uniform dynamic validations.
3. Extend the two-level Galerkin preconditioner into a recursive V-cycle and add aperture-aware EB reconstruction when fabric reaches a 2:1 interface.
4. Extend the opened-cavity topology test into a developed internal mass-flow test once compact aperture velocities receive full advection/diffusion updates.
5. Complete the Qt controls/slices/scene workflow and only then remove building/channel/ground/seabed/porous code.

## Important files

- Geometry: `src/core/geometry/tri_mesh.h`, `step_import.*`, `triangle_bvh.*`, `embedded_boundary.*`
- AMR fields/exchange: `src/core/fluid/amr_grid.*`, `amr_fields.*`, `amr_exchange.*`
- AMR EB/pressure: `src/core/fluid/amr_eb.*`, `amr_pressure.*`, `eb_pressure.*`
- Loads/config: `src/core/aero_loads.*`, `src/core/paraglider_config.*`
- Tests/probes: `tools/paraglider_geometry_probe.cpp`, `paraglider_gpu_probe.cpp`, `paraglider_probe.cpp`

## Validation policy

Do not bless changed output merely because an operator changed. Validate redesigned geometry and operators against analytical or CPU-double references. Production fields are FP32; geometry, residual accumulation where useful, and force totals remain FP64. Never hide unresolved topology or cross-fabric connections by loosening tolerances.
