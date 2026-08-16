# ParaCFD

ParaCFD is a C++20/CUDA/Qt6 aerodynamic CFD application for stabilized paraglider geometry imported from STEP. It treats canopy skins, ribs, diagonals, and panels as zero-thickness, two-sided impermeable fabric; fluid exists on both sides, intentional CAD openings remain open, and no solid-volume filling or artificial extrusion is used.

## Physical model

The intended workflow is one-way coupling:

```
XPBD equilibrium -> STEP export -> ParaCFD preprocessing -> static CFD run
```

OpenCascade remains isolated in the geometry layer. STEP faces are tessellated in SI metres, with source-face provenance retained on every triangle. Triangle winding defines that triangle's local minus-to-plus normal only. Reversing winding swaps the pressure sides and normal together, leaving the physical pressure force invariant:

```
F = (p_minus - p_plus) A n
```

Missing triangles are openings. There is no separate interior-air material: ram-air cells, crossports, and the external flow contain the same fluid.

## Geometry preprocessing

- `TriMesh` stores positions, indices, source STEP face IDs, and optional per-corner UV coordinates.
- `TriangleBvh` is a static double-precision CPU BVH supporting cell AABB queries, segment intersections, nearest-surface queries, and distance queries.
- Runtime geometry is placed before BVH construction. Normals use the affine inverse transpose.
- CAD preprocessing is static. CUDA timesteps do not traverse CAD, rebuild topology, or intersect triangles.

The OpenCascade linear deflection is exposed in the paraglider configuration. The current default is 2 mm, which is still substantially finer than the initial 62.5 mm CFD surface cells and avoids unnecessary million-triangle previews. It should eventually be derived from finest CFD spacing.

## Grid and fields

The new grid path is static block-structured Cartesian AMR with refinement ratio 2:1. A brick has uniform cells and compact metadata: level, integer coordinate, origin, spacing, six same-level neighbours, parent, children, and boundary/embedded-boundary flags. Per-level integer-coordinate hash tables locate bricks without tree traversal or binary searches.

Fields use pooled, contiguous Structure-of-Arrays storage per AMR level. Production `Real` is `float`; `PARACFD_VALIDATION_FP64=ON` changes the new numerical path to `double` for validation. Geometry preprocessing and global validation calculations may still use double precision. A GPU-resident per-level coordinate hash locates the finest brick and local cell directly.

The implemented AMR path includes same-level halo exchange, ratio-two restriction/prolongation, aperture-aware flux matching, hierarchy balancing, and one matrix-free composite pressure operator spanning all active levels. Coarse/fine faces are compact four-tile connections; ordinary same-level faces remain implicit structured stencils.

## Zero-thickness embedded boundaries

Regular cells retain implicit Cartesian topology. Only cells intersected by fabric receive compact irregular records:

- fluid fragments with independent pressure state, volume, centroid, and parent cell;
- fragment connections through open fluid;
- Cartesian-face apertures with area, centroid, and the fluid fragments they join;
- fabric surface patches with source triangle/face IDs, area, centroid, normal, and plus/minus fragments.

A smooth membrane cut produces two fluid fragments and no connection through fabric. Fabric exactly coincident with a Cartesian face becomes a blocked cut face between two full-volume regular cells; partial coverage is reconstructed into separate connected apertures. Complex finest-level cells use a configurable deterministic local fluid-connectivity graph whose edges are blocked by exact BVH intersections. This supports more than two fragments at ribs, sheet endings, and junctions without inside/outside filling or cross-fabric merging; any topology still unresolved is a hard preprocessing failure.

Small fragments below `min_volume_fraction` are conservatively merged through suitable open apertures to a larger same-side control volume. Fabric is never a legal merge connection. An isolated sealed pocket retains its pressure state and volume but is excluded from the projection until component-wise gauge/nullspace handling is implemented; the preprocessing report counts these states explicitly.

## Pressure projection and aerodynamic loads

The matrix-free finite-volume pressure operator uses actual control-volume volumes and coefficients proportional to open area divided by centre distance. Regular/regular faces retain a structured fast path; EB apertures, blocked coincident faces, and 2:1 coarse/fine tiles are compact special cases. A sparse cross-brick atlas prevents brick boundaries from becoming physical walls and augments one composite pressure graph spanning every level. CPU reference and FP32/FP64 CUDA operator implementations share the topology. A GPU-resident Jacobi-PCG solve and manufactured divergence/gradient correction tests cover both EB and coarse/fine fluxes. GPU boundary kernels prescribe +X freestream, use an X-max pressure reference/zero-gradient velocity, and apply the same free-slip far-field condition at both Y and both Z boundaries; there is no ground branch.

Surface pressure results retain `p_plus`, `p_minus`, `delta_p`, the corresponding Cp values, and pressure force per source triangle. Whole-wing force and moment use double-precision accumulation. With freestream along +X, drag, side, and lift axes are +X, +Y, and +Z. Coefficients are reported only when the user supplies a positive reference area; the present aerodynamic load is explicitly pressure-only because fabric skin friction is not implemented.

## Visualization collision

The Qt/OpenGL viewer still displays the STEP surface. Flow particles and tracers can use the placed triangle BVH for segment-crossing tests and cannot pass through fabric. This does not use an inside-solid or parity classification. The mesh renderer accepts per-source-triangle delta-Cp through an OpenGL shader-storage buffer indexed by primitive ID, so shared CAD vertices do not smear results across fabric panels. The aerodynamic timestep has not yet been wired to publish that buffer automatically.

## Build and test

The RTX 4090 development configuration is:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build
ctest --test-dir build --output-on-failure
```

The repository may already contain a Visual Studio build directory; preserve and use its established generator rather than reconfiguring it in place.

Important probes are:

- `paraglider_geometry_probe`: deterministic geometry, EB topology, pressure/operator, AMR exchange, FP32 parity, and load tests;
- `paraglider_gpu_probe`: CUDA pressure-operator timing and allocation estimates;
- `paraglider_probe`: STEP/config preprocessing diagnostics;
- `parity_probe`: retained CPU/GPU parity coverage for useful legacy kernels.

`configs/paraglider.json` is the new configuration reference.
`configs/planb_parakite.json` is the checked-in PlanB acceptance case; its +90-degree Z placement maps that exporter’s forward direction (-Y) to ParaCFD freestream (+X).

## Current limitations

The repository is in an incremental migration state and must not yet be described as a trustworthy end-to-end paraglider CFD product:

- EB that reaches a 2:1 interface is currently rejected instead of receiving aperture-aware cross-level fragment reconstruction. The tested PlanB hierarchy keeps its EB atlas away from these interfaces.
- The composite solver currently uses Jacobi-PCG, not a geometric multigrid preconditioner; production-size convergence and component-wise pressure gauges still need work.
- The existing MacCormack and Smagorinsky kernels have not yet been ported to AMR-aware, side-aware sampling on the new fields.
- The new external-aerodynamic timestep is not wired end to end through advection, LES, boundary conditions, projection, and statistics.
- The Qt viewer can inspect AMR bricks and owned EB cells and reports whether the composite pressure topology is ready. It deliberately refuses to run the legacy channel timestep for a loaded paraglider, but its controls remain substantially inherited from the building/channel product.
- Skin-friction/wall-model force is absent; reported new-path force is pressure-only.
- Required full-flow validations such as the opening-cavity case and AMR-versus-uniform force comparison remain outstanding.

The legacy channel/building code is retained only as a numerical and visualization reference while these missing replacements are completed. Do not extend it as the new architecture, and do not delete it until the replacement path supplies equivalent working functionality and tests.
