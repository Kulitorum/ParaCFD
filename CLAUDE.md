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

The implemented AMR path includes same-level halo exchange, ratio-two restriction/prolongation, aperture-aware flux matching, hierarchy balancing, and one matrix-free composite pressure operator spanning all active levels. Coarse/fine faces are compact four-tile connections; ordinary same-level faces remain implicit structured stencils. After transport, every compact tile gathers its collocated fine MAC velocity; after projection, corrected tiles scatter to the fine faces and their aperture-weighted mean is written to the coarse face.

`ExternalAeroCore` owns the static hierarchy/topology and persistent device fields. Its current timestep is bounded MacCormack AMR advection with RK2 characteristics, explicit molecular/Smagorinsky diffusion, external boundary conditions, compact EB-aperture transport, conservative coarse/fine synchronization, then the composite projection. A GPU max reduction chooses one global timestep from the finest spacing and current regular-flow velocity; only its scalar result returns to the host. Compact EB transport bounds its own local Courant/diffusion update, so a negligible geometric sliver cannot throttle the whole domain. Backtraces locate the finest active brick through the GPU hash. Trilinear and diffusion stencils cross same-level brick boundaries without clamping, and every level reads common old-time and forward states before correction/commit. Static preprocessing marks face centres within 2.5 local cells of fabric; those faces retain their local value rather than tracing across a sheet. This guarantees side safety but remains deliberately dissipative inside the protection band.

## Zero-thickness embedded boundaries

Regular cells retain implicit Cartesian topology. Only cells intersected by fabric receive compact irregular records:

- fluid fragments with independent pressure state, volume, centroid, and parent cell;
- fragment connections through open fluid;
- Cartesian-face apertures with area, centroid, and the fluid fragments they join;
- fabric surface patches with source triangle/face IDs, area, centroid, normal, and plus/minus fragments.

A smooth membrane cut produces two fluid fragments and no connection through fabric. Fabric exactly coincident with a Cartesian face becomes a blocked cut face between two full-volume regular cells; partial coverage is reconstructed into separate connected apertures. Complex finest-level cells use a configurable deterministic local fluid-connectivity graph whose edges are blocked by exact BVH intersections. This supports more than two fragments at ribs, sheet endings, and junctions without inside/outside filling or cross-fabric merging; any topology still unresolved is a hard preprocessing failure.

Independent EB aperture-normal velocities are advanced on a compact same-side graph. Structured MAC faces adjacent to each fragment supply the carrier state, and each aperture uses a conservative first-order upwind update plus molecular diffusion before projection. The graph contains only actual fluid connections, so it cannot transport through fabric and it does not turn the regular domain into CSR. This first implementation deliberately omits Smagorinsky viscosity and higher-order reconstruction on the compact graph.

Small fragments below `min_volume_fraction` (default 0.25 of a Cartesian cell) are conservatively merged through suitable open apertures to a larger same-side control volume. Fabric is never a legal merge connection. Apertures below `min_aperture_area_fraction * h^2` (default `1e-4`) are finite-resolution clipping slivers: preprocessing removes the matching velocity and pressure connection together and reports their count and cumulative atlas area. This threshold does not close a resolved CAD opening, and can be set to zero for geometry studies. An isolated sealed pocket retains its pressure state and volume but is excluded from the projection until component-wise gauge/nullspace handling is implemented; the preprocessing report counts these states explicitly.

## Pressure projection and aerodynamic loads

The matrix-free finite-volume pressure operator uses actual control-volume volumes and coefficients proportional to open area divided by centre distance. Regular/regular faces retain a structured fast path; EB apertures, blocked coincident faces, and 2:1 coarse/fine tiles are compact special cases. A sparse cross-brick atlas prevents brick boundaries from becoming physical walls and augments one composite pressure graph spanning every level. CPU reference and FP32/FP64 CUDA operator implementations share the topology. A persistent GPU projection performs divergence assembly, RHS formation, one composite PCG solve, pressure scatter, and conservative regular/special-flux correction without bulk field transfers. The PCG preconditioner combines the fine diagonal with an additive level-zero Galerkin aggregate solve; compact nonlocal aggregate edges retain conservative fragment-merge topology. Each active fluid component that cannot reach the pressure outlet receives a deterministic gauge. GPU boundary kernels prescribe +X freestream, use an X-max pressure reference/zero-gradient velocity, and apply the same free-slip far-field condition at both Y and both Z boundaries; there is no ground branch.

Every owned composite surface patch retains its global plus/minus pressure DOFs and CAD triangle/face provenance. Surface pressure results retain `p_plus`, `p_minus`, `delta_p`, the corresponding Cp values, and pressure force per source triangle. Whole-wing force and moment use double-precision accumulation. With freestream along +X, drag, side, and lift axes are +X, +Y, and +Z. Coefficients are reported only when the user supplies a positive reference area; the present aerodynamic load is explicitly pressure-only because fabric skin friction is not implemented.

## Visualization collision

The Qt/OpenGL viewer still displays the STEP surface. Flow particles and tracers can use the placed triangle BVH for segment-crossing tests and cannot pass through fabric. This does not use an inside-solid or parity classification. The mesh renderer accepts per-source-triangle delta-Cp through an OpenGL shader-storage buffer indexed by primitive ID, so shared CAD vertices do not smear results across fabric panels. A dedicated worker now advances `ExternalAeroCore`, throttles surface/conservation downloads, publishes delta-Cp automatically, and reports instantaneous pressure-only loads, CFL, regular/compact-EB velocity maxima, projection convergence, divergence/flux diagnostics, and persistent GPU memory.

The purpose-built paraglider window controls density/viscosity, face-bbox domain margins, base spacing/level count, surface/wing/wake refinement, CFL, Smagorinsky coefficient, projection tolerance/iterations, and optional aerodynamic reference values. `Build CFD Grid` reconstructs the static hierarchy and resets the flow. The internal freestream axis remains +X. On direct STEP load, the longer horizontal face-only bbox axis is treated as span and the shorter as chord; ParaCFD rotates either known `-X`-forward or `-Y`-forward exporter layout into +X flow. A prominent 180-degree LE/TE flip remains because bbox extents cannot determine chord polarity. Config files store the confirmed placement matrix and do not re-infer it.

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
- `paraglider_flow_probe`: dynamic normal/parallel/inclined plates, projection divergence, opened-cavity connectivity, and AMR-versus-uniform force comparison;
- `paraglider_probe`: STEP/config preprocessing diagnostics;
- `paraglider_case_probe`: imported-wing multi-step force, CFL, conservation, and performance history (manual long-running diagnostic);
- `parity_probe`: CPU/GPU parity and immutable snapshots for the reusable FP64 uniform-MAC reference kernels.

`configs/paraglider.json` is the new configuration reference.
`configs/planb_parakite.json` is the checked-in PlanB acceptance case; its +90-degree Z placement maps that exporter’s forward direction (-Y) to ParaCFD freestream (+X).

## Current limitations

ParaCFD is now a paraglider-only application, but it must not yet be described as producing trustworthy end-to-end paraglider aerodynamic results:

- EB that reaches a 2:1 interface is currently rejected instead of receiving aperture-aware cross-level fragment reconstruction. The tested PlanB hierarchy keeps its EB atlas away from these interfaces.
- The composite solver has a two-level additive geometric/Galerkin preconditioner, not yet a complete recursive multigrid hierarchy. The checked-in PlanB initial projection converges at the configured global residual tolerance, but local maximum divergence remains sensitive to the very small irregular control volumes and needs a stricter local/conservation acceptance criterion.
- The new external-aero timestep is wired end to end with bounded MacCormack/RK2 advection. Fabric-near faces still use a conservative static protection fallback. Same-level brick crossings and normal 2:1 interface fluxes are synchronized conservatively; general cross-level interpolation remains nonconservative.
- Compact EB aperture velocities use a first-order same-side graph transport with molecular diffusion. Conservative 0.25-volume merging and the reportable `1e-4 h^2` aperture cutoff eliminate the previously observed compact-state runaway, but higher-order reconstruction and consistent Smagorinsky treatment are still required for force convergence.
- A 500-step PlanB diagnostic keeps the compact peak below the structured-field peak at the end, but the regular peak still grows to about 56 m/s for a 10 m/s inlet and the pressure force is still evolving. Grid/force convergence, fabric-near transport accuracy, and orientation/AoA confirmation remain necessary before the result is aerodynamically trustworthy.
- The Qt viewer displays live velocity/pressure slices, arrows, fabric-stopped tracers, AMR/EB boxes, and triangle-native visible-side Cp, Cp+, Cp-, or delta-Cp. Visualization currently uses throttled host snapshots and OpenGL uploads rather than direct AMR CUDA/OpenGL interop.
- Skin-friction/wall-model force is absent; reported new-path force is pressure-only.
- The opened-cavity gate measures developed bidirectional flow through a deliberately missing fabric face while the closed control has zero represented opening flux. Longer internal-pressure and force-convergence studies on resolved paraglider geometry remain outstanding.

The former building, binary-solid voxel, ground/channel, seabed, porous, scene, and building GUI paths have been removed. Only mathematically reusable uniform-MAC FP64 kernels remain as validation references; they are not linked into the production external-aero timestep.
