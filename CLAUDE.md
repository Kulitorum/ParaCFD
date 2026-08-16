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

`ExternalAeroCore` owns the static hierarchy/topology and persistent device fields. Its current timestep is bounded MacCormack AMR advection with RK2 characteristics, explicit molecular/Smagorinsky diffusion, external boundary conditions, compact EB-aperture transport, conservative coarse/fine synchronization, then the composite projection. A GPU max reduction over physical faces of active bricks chooses one global timestep from the finest spacing and current regular-flow velocity; ghost storage and covered coarse bricks are excluded, and only the scalar result returns to the host. Compact EB transport bounds its own local Courant/diffusion update, so a negligible geometric sliver cannot throttle the whole domain. Backtraces locate the finest active brick through the GPU hash. Trilinear and diffusion stencils cross same-level brick boundaries without clamping. At 2:1 transitions, staggered face values and cell-centred eddy viscosity use one-sided linear prolongation; the reverse MacCormack correction retains the bounded forward RK2 value if its trace changes lattice. Every level reads common old-time and forward states before correction/commit. Static preprocessing marks face centres within 2.5 local cells of fabric and stores six one-byte same-side link bits. Marked faces use bounded minmod-limited linear transport only when the complete upwind/downwind stencil is connected on the same fluid side, with a first-order donor fallback at sheet edges and unavailable level/boundary links; diffusion is likewise restricted to open links. Far faces retain MacCormack/LES. This prevents opposite-side stencil data without per-step triangle traversal while reducing fabric-band diffusion in smooth tangential flow.

## Zero-thickness embedded boundaries

Regular cells retain implicit Cartesian topology. Only cells intersected by fabric receive compact irregular records:

- fluid fragments with independent pressure state, volume, centroid, and parent cell;
- fragment connections through open fluid;
- Cartesian-face apertures with area, centroid, and the fluid fragments they join;
- fabric surface patches with source triangle/face IDs, area, centroid, normal, and plus/minus fragments.

A smooth membrane cut produces two fluid fragments and no connection through fabric. Fabric exactly coincident with a Cartesian face becomes a blocked cut face between two full-volume regular cells; partial coverage is reconstructed into separate connected apertures. Complex finest-level cells use a configurable deterministic local fluid-connectivity graph whose edges are blocked by exact BVH intersections. This supports more than two fragments at ribs, sheet endings, and junctions without inside/outside filling or cross-fabric merging; any topology still unresolved is a hard preprocessing failure.

Independent EB aperture-normal velocities are advanced on a compact same-side graph. Structured MAC faces adjacent to each fragment supply the carrier state. Static preprocessing links nearest same-axis apertures into directed chains; complete upstream/upstream2/downstream chains use a minmod-limited MUSCL update, while sheet endings, junctions, and incomplete chains retain the bounded first-order donor fallback. An area-weighted graph reconstruction estimates all nine velocity-gradient components from connected fragment states; its symmetric tensor supplies `|S| = sqrt(2 Sij Sij)` to the compact-region Smagorinsky viscosity. Molecular diffusion and that eddy viscosity are applied before projection. The combined explicit update is limited to the extrema of its actual same-side stencil, so transport cannot invent a new local maximum or minimum; there is no arbitrary absolute velocity cap, and the subsequent pressure projection remains free to create physically required extrema. The graph contains only actual fluid connections, so it cannot transport through fabric and it does not turn the regular domain into CSR. This axis-column reconstruction still needs a multidimensional least-squares treatment at irregular junctions and an arbitrary-orientation fabric wall model.

Small fragments below `min_volume_fraction` (default 0.25 of a Cartesian cell) are conservatively merged through suitable open apertures to a larger same-side control volume. Fabric is never a legal merge connection. A 0.05 experiment retained substantially stiffer PlanB compact states, so 0.25 remains the stability-oriented default while resolution studies quantify its geometric influence. Apertures below `min_aperture_area_fraction * h^2` (default `1e-4`) are finite-resolution clipping slivers: preprocessing removes the matching velocity and pressure connection together and reports their count and cumulative atlas area. This threshold does not close a resolved CAD opening, and can be set to zero for geometry studies. Isolated sealed fluid components retain their pressure states and volumes; the composite solve assigns one explicit pressure gauge to every component not connected to the pressure-reference outlet.

## Pressure projection and aerodynamic loads

The matrix-free finite-volume pressure operator uses actual control-volume volumes and coefficients proportional to open area divided by centre distance. Regular/regular faces retain a structured fast path; EB apertures, blocked coincident faces, and 2:1 coarse/fine tiles are compact special cases. A sparse cross-brick atlas prevents brick boundaries from becoming physical walls and augments one composite pressure graph spanning every level. CPU reference and FP32/FP64 CUDA operator implementations share the topology. A persistent GPU projection performs divergence assembly, RHS formation, one composite PCG solve, pressure scatter, and conservative regular/special-flux correction without bulk field transfers. The PCG preconditioner combines the fine diagonal with an additive level-zero Galerkin aggregate solve; compact nonlocal aggregate edges retain conservative fragment-merge topology. Each active fluid component that cannot reach the pressure outlet receives a deterministic gauge. GPU boundary kernels prescribe +X freestream, use an X-max pressure reference/zero-gradient velocity, and apply the same free-slip far-field condition at both Y and both Z boundaries; there is no ground branch.

Every owned composite surface patch retains its global plus/minus pressure DOFs and CAD triangle/face provenance. Surface pressure results retain `p_plus`, `p_minus`, `delta_p`, the corresponding Cp values, and pressure force per source triangle. Whole-wing force and moment use double-precision accumulation. With freestream along +X, drag, side, and lift axes are +X, +Y, and +Z. Coefficients are reported only when the user supplies a positive reference area; the present aerodynamic load is explicitly pressure-only because fabric skin friction is not implemented.

## Visualization collision

The Qt/OpenGL viewer still displays the STEP surface. Flow particles and tracers can use the placed triangle BVH for segment-crossing tests and cannot pass through fabric. This does not use an inside-solid or parity classification. The mesh renderer accepts per-source-triangle delta-Cp through an OpenGL shader-storage buffer indexed by primitive ID, so shared CAD vertices do not smear results across fabric panels. A dedicated worker now advances `ExternalAeroCore`, throttles surface/conservation downloads, publishes delta-Cp automatically, and reports instantaneous pressure-only loads, CFL, regular/compact-EB velocity maxima, projection convergence, divergence/flux diagnostics, and persistent GPU memory. The scalar plane samples the finest active AMR brick at every display vertex from the throttled per-level host snapshot; its two axis-dependent dimensions are derived from the domain extent and finest active cell spacing rather than a fixed square resolution, with a bounded visualization-only allocation cap. Arrows and tracers intentionally retain a cheaper coarse 3-D display field. The default camera fits the placed wing; `Fit Wing` and `Fit Domain` explicitly separate close aerodynamic inspection from inspection of the longer downstream wake volume. Live controls expose arrow mode/density/animation speed/uniform size and tracer mode/density/length/width/boring-filter settings; the calmer default arrow animation is visual only. Two-dimensional arrows draw as a depth-independent plane annotation, while 3-D arrows retain ordinary occlusion. Editors ignore ordinary wheel input so scrolling the compact settings dock cannot silently change simulation or visualization values.

The purpose-built paraglider window controls density/viscosity, face-bbox domain margins, base spacing/level count, surface/wing/wake refinement, CFL, Smagorinsky coefficient, projection tolerance/iterations, and optional aerodynamic reference values. `Build CFD Grid + Run` reconstructs the static hierarchy, initializes pressure, and immediately starts continuous CFD; Play/Pause and Step control the resulting worker. The internal freestream velocity remains +X, so the leading edge must face upstream toward -X. On direct STEP load, the longer horizontal face-only bbox axis is treated as span and the shorter as chord; the current exporter convention assumes the negative chord direction points toward the leading edge and maps it to upstream -X. A prominent 180-degree LE/TE flip remains because bbox extents cannot determine chord polarity. Config files store the confirmed placement matrix and do not re-infer it. GUI and command-line preprocessing share one brick-aligned zero-origin framing helper; placement is canonicalized to the FP32 display/production-mesh frame before clipping so both paths build identical EB topology.

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
`configs/planb_parakite.json` is the checked-in PlanB acceptance case; its -90-degree Z placement maps the confirmed leading-edge direction (-Y) to ParaCFD upstream (-X), against the +X freestream velocity.

## Current limitations

ParaCFD is now a paraglider-only application, but it must not yet be described as producing trustworthy end-to-end paraglider aerodynamic results:

- EB that reaches a 2:1 interface is currently rejected instead of receiving aperture-aware cross-level fragment reconstruction. The tested PlanB hierarchy keeps its EB atlas away from these interfaces.
- The composite solver has a two-level additive geometric/Galerkin preconditioner, not yet a complete recursive multigrid hierarchy. The checked-in PlanB initial projection converges at the configured global residual tolerance, but local maximum divergence remains sensitive to the very small irregular control volumes and needs a stricter local/conservation acceptance criterion.
- The new external-aero timestep is wired end to end with bounded MacCormack/RK2 advection. Fabric-near faces use a precomputed same-side link graph with minmod-limited linear transport where the full connected stencil exists, a first-order fallback at incomplete stencils, and link-restricted molecular diffusion. Same-level brick crossings and normal 2:1 interface fluxes are synchronized conservatively. Cross-level face/cell prolongation is now linearly exact to the FP32 geometry-coordinate floor, but semi-Lagrangian momentum advection is still not a globally conservative refluxed finite-volume update.
- Compact EB aperture velocities use bounded minmod-limited MUSCL transport on complete same-axis graph chains, with a first-order donor fallback where a chain ends or branches, plus molecular diffusion and a reconstructed nine-component graph-strain Smagorinsky estimate. Conservative 0.25-volume merging and the reportable `1e-4 h^2` aperture cutoff remove numerical slivers. Tests cover a smooth quadratic characteristic, uniform-state preservation, local bounds, LES dissipation, tangential shear, rigid rotation, and normal extension. This is not yet a least-squares junction reconstruction or a fabric wall model, and force convergence is still required before trusting loads.
- On the current 47,997-triangle PlanB working fixture with the confirmed leading edge facing upstream, the three-level 500-step diagnostic reaches 0.527 s with 12.76/45.45 m/s regular/compact peaks and pressure-only force `[174.66, 0.02, -5.27]` N. The immediately preceding bounded-MUSCL build using only the graph-normal derivative reached a 55.05 m/s compact peak after the same step count, so the reconstructed strain lowers that peak by about 17%. Volume-weighted RMS divergence is `4.18e-6 s^-1`, net integrated flux error is `3.89e-6 m^3/s`, and persistent storage is 496.98 MiB. The compact peak and loads remain unresolved/grid-dependent and are explicit convergence blockers, not values hidden by a velocity guard. Full grid/domain/force convergence and AoA confirmation remain necessary.
- The Qt viewer displays live velocity/pressure slices, arrows, fabric-stopped tracers, AMR/EB boxes, and triangle-native visible-side Cp, Cp+, Cp-, or delta-Cp. The scalar slice is finest-brick aware and its per-axis sample counts follow finest cell spacing, but visualization still uses throttled host snapshots and OpenGL uploads rather than direct CUDA/OpenGL interop; arrows and tracers use a coarse uniform 3-D display copy.
- Skin-friction/wall-model force is absent; reported new-path force is pressure-only.
- The opened-cavity gate measures developed bidirectional flow through a deliberately missing fabric face while the closed control has zero represented opening flux. Longer internal-pressure and force-convergence studies on resolved paraglider geometry remain outstanding.

The former building, binary-solid voxel, ground/channel, seabed, porous, scene, and building GUI paths have been removed. Only mathematically reusable uniform-MAC FP64 kernels remain as validation references; they are not linked into the production external-aero timestep.
