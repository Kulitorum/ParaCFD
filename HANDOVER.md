# ParaCFD handover

## Current state

ParaCFD is being migrated from the building/channel solver to a purpose-built static-geometry paraglider solver. The repository builds, and the new foundation has deterministic CPU/GPU tests. The physical invariant is zero-thickness, two-sided fabric: never extrude it, flood-fill it, close its openings, or replace it with a binary solid mask.

Implemented:

- STEP triangle source-face provenance and optional UV retention;
- affine-safe placed normals;
- robust static double-precision triangle BVH;
- paraglider configuration and automatic far-field bounding box;
- 2:1 static block-AMR hierarchy, CPU/GPU brick hashing, level-pooled FP32 SoA fields, and CPU/CUDA same-level halos;
- ratio-two restriction, prolongation, balancing, and aperture-aware flux matching primitives;
- one-level zero-thickness EB cells with split fluid fragments, split Cartesian-face apertures, coincident blocked faces without tiny fragments, surface patches, explicit unresolved-complexity reporting, and conservative same-side small-fragment merging (including merge chains);
- one-level finite-volume EB pressure operator, CPU PCG/reference path, and CUDA FP32 matrix-free operator/CG path;
- two-sided pressure/Cp/pressure-force accumulation with winding-invariant force and optional reference coefficients;
- BVH segment collision for particles and tracers, plus triangle-native delta-Cp rendering support;
- GPU +X inflow, X-max pressure outlet, and symmetric Y/Z free-slip boundary kernels with no ground;
- geometry, operator, AMR, configuration, and GPU timing probes.

The old solver remains buildable as a reference. It must not be confused with the new aerodynamic path.

## Next engineering work

1. Assemble EB topology across adjacent bricks, including fragment/aperture identity at brick interfaces and physical boundary records.
2. Implement composite AMR pressure coupling. Coarse flux must equal the sum of fine aperture fluxes, and the preconditioner must operate across levels rather than solving levels independently.
3. Port velocity storage, boundary conditions, MacCormack advection, and Smagorinsky LES to the AMR hierarchy. Near fabric, BVH/EB segment-side tests must prevent a backtrace from sampling the opposite pressure side.
4. Add the complete external-aero timestep and manufactured GPU flow cases: normal/parallel/inclined plates, opened cavity, AMR conservation, and AMR-versus-uniform comparison.
5. Replace the legacy Qt workflow with STEP -> Build CFD Grid -> Start Simulation; then add brick/EB overlays, two-sided Cp colouring, pressure-only force readouts, residuals, timings, and memory statistics.
6. Only after the new path is operational, remove building, channel, ground/seabed, porous, old voxel-load, centerline, and building-config code and simplify CMake.

## Important files

- Geometry: `src/core/geometry/tri_mesh.h`, `step_import.cpp`, `triangle_bvh.*`, `embedded_boundary.*`
- AMR: `src/core/fluid/amr_grid.*`, `amr_fields.*`, `amr_exchange.*`
- Pressure: `src/core/fluid/eb_pressure.*`
- Loads/config: `src/core/aero_loads.*`, `src/core/paraglider_config.*`
- Tests/probes: `tools/paraglider_geometry_probe.cpp`, `paraglider_gpu_probe.cpp`, `paraglider_probe.cpp`
- New default: `configs/paraglider.json`

## Validation policy

Do not bless old output merely because the operator changed. Keep mathematically unchanged legacy parity tests. Validate redesigned geometry/operators against analytical or CPU-double references. Production fields are FP32, while geometry, residual accumulation where useful, force totals, and reference calculations may be FP64.

Do not hide unresolved geometry at the finest level. Emit the affected brick/cell and triangle information and reject the CFD grid until it is resolved or the configured refinement limit is intentionally increased.
