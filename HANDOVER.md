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
- two-sided pressure/Cp/pressure-force accumulation with winding-invariant force;
- BVH tracer collision and triangle-native delta-Cp render storage;
- a GUI AMR/EB preview that validates pressure-topology construction and refuses to run the unrelated legacy channel timestep.

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

## Next engineering work

1. Implement the new external-aero timestep around the persistent projection: AMR-aware advection and LES with side-safe backtraces, then boundary conditions, projection, and statistics.
2. Retain composite surface-patch-to-pressure-DOF mappings and publish p+/p-/delta-Cp and pressure-only forces from a converged flow.
3. Extend the two-level Galerkin preconditioner into a recursive V-cycle and add aperture-aware EB reconstruction when fabric reaches a 2:1 interface.
4. Add dynamic normal/parallel/inclined plate and opened-cavity tests, followed by AMR-versus-uniform force validation.
5. Complete the Qt workflow and only then remove building/channel/ground/seabed/porous code.

## Important files

- Geometry: `src/core/geometry/tri_mesh.h`, `step_import.*`, `triangle_bvh.*`, `embedded_boundary.*`
- AMR fields/exchange: `src/core/fluid/amr_grid.*`, `amr_fields.*`, `amr_exchange.*`
- AMR EB/pressure: `src/core/fluid/amr_eb.*`, `amr_pressure.*`, `eb_pressure.*`
- Loads/config: `src/core/aero_loads.*`, `src/core/paraglider_config.*`
- Tests/probes: `tools/paraglider_geometry_probe.cpp`, `paraglider_gpu_probe.cpp`, `paraglider_probe.cpp`

## Validation policy

Do not bless changed output merely because an operator changed. Validate redesigned geometry and operators against analytical or CPU-double references. Production fields are FP32; geometry, residual accumulation where useful, and force totals remain FP64. Never hide unresolved topology or cross-fabric connections by loosening tolerances.
