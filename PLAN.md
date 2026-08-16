# Paraglider CFD implementation plan

This file records actual migration status, not a claim that unfinished architecture exists.

## A. Geometry foundation — implemented

- Preserve OpenCascade STEP import and metre conversion.
- Retain source STEP face ID and optional UV data per tessellated triangle.
- Build a static triangle BVH with AABB, segment, nearest-point, and distance queries.
- Add paraglider configuration and external-aero data structures.
- Keep the old solver usable only as a reference during migration.

## B. Uniform brick grid — implemented foundation

- Uniform one-level block hierarchy behaves as one Cartesian grid.
- Compact brick metadata and CPU/GPU per-level coordinate hash.
- Pooled SoA host/device `Real` fields and same-level halo exchange.

The complete external-aero velocity timestep does not yet consume these fields.

## C. One-level zero-thickness EB — implemented foundation

- Exact triangle/cell clipping for a resolved planar membrane.
- Independent plus/minus fluid fragments and fabric patches.
- Split Cartesian-face aperture records and coincident blocked-face handling.
- Explicit complex-cell detection for refinement without cross-fabric merging.
- Conservative same-side small-fragment merge.

Next: connect topology across bricks and improve deterministic subdivision of cells containing multiple coplanar or meeting surface patches.

## D. One-level EB projection — implemented foundation

- Matrix-free finite-volume operator with actual volume, aperture area, and distance.
- Fast implicit regular faces plus compact irregular aperture work.
- CPU-double reference/PCG and CUDA `Real` operator/CG.
- Manufactured independent-side pressure/load and projection tests.

Next: add robust nullspace treatment for every disconnected all-Neumann fluid component and integrate velocity-face correction into the full timestep.

## E. Static AMR — exchange primitives implemented, composite solve pending

- 2:1 hierarchy, balancing, geometry/wake refinement, hash lookup.
- Halo, restriction, prolongation, and aperture-aware reflux helpers.

Pending:

- cross-brick and cross-level EB topology;
- conservative composite divergence and pressure gradient;
- multilevel matrix-free projection/preconditioner;
- AMR-versus-uniform flow validation.

## F. FP32 production — storage/operator foundation implemented

- New production scalar is FP32 with an optional FP64 validation build.
- Geometry preprocessing and force/reference calculations retain FP64 where valuable.
- CUDA pressure operator has measured FP32/FP64 parity and timing probes.

Pending: port the entire timestep and quantify full-case convergence, memory, and performance.

## G. Aerodynamics — manufactured load path implemented

- Per-triangle plus/minus pressure and Cp, delta-Cp, force, whole-wing force/moment.
- Winding-invariant pressure force.
- CL/CD/CS only with explicit positive reference area.
- Pressure drag is labelled pressure-only.

Pending: obtain the pressure states from a complete converged paraglider flow and add surface field transfer/rendering.

## H. GUI migration — pending

The STEP viewer and CUDA/OpenGL infrastructure remain useful. Particle/tracer fabric collision now uses the triangle BVH. The primary controls, scene schema, AMR/EB debug rendering, two-sided surface colouring, and aerodynamic dashboard still require migration.

## I. Obsolete-code removal — intentionally pending

Remove the building, centerline, roof/wall, solid-voxel-load, channel/ground/seabed/porous, building probe/config, and corresponding UI/doc paths only after the new timestep and UI replacements pass their validation cases. Until then they are reference code, not current architecture.

## Definition of trustworthy version 1

- No pressure/velocity connection crosses fabric.
- Openings remain fluid-connected.
- Composite AMR flux is conservative through EB and coarse/fine interfaces.
- GPU projection reaches the declared residual tolerance with a documented norm.
- FP32 agrees with CPU/FP64 reference cases at a justified tolerance.
- Normal, parallel, and inclined plates; opened cavity; AMR interface; and AMR-versus-uniform cases pass.
- End-to-end timestep contains no per-step CPU geometry work or bulk field transfers.
- GUI and logs clearly distinguish pressure drag from total drag and expose unresolved geometry.
