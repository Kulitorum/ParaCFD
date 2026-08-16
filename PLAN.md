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

The external-aero core now consumes these pooled fields directly.

## C. One-level zero-thickness EB — implemented foundation

- Exact triangle/cell clipping for a resolved planar membrane.
- Independent plus/minus fluid fragments and fabric patches.
- Split Cartesian-face aperture records and coincident blocked-face handling.
- Explicit complex-cell detection for refinement without cross-fabric merging.
- Conservative same-side small-fragment merge.

The cross-brick atlas and deterministic multi-fragment fallback are implemented. Partial Cartesian-face coverage retains separate connected apertures. Remaining: exact cross-level EB topology where fabric reaches a 2:1 interface.

## D. One-level EB projection — implemented foundation

- Matrix-free finite-volume operator with actual volume, aperture area, and distance.
- Fast implicit regular faces plus compact irregular aperture work.
- CPU-double reference/PCG and CUDA `Real` operator/CG.
- Manufactured independent-side pressure/load and projection tests.

The composite CPU reference and persistent GPU projection now cover EB divergence, RHS, pressure, regular/special-flux correction, and post-projection divergence. Disconnected active fluid components receive deterministic gauges. Remaining: integrate the projection into the external-aero timestep.

## E. Static AMR — composite pressure foundation implemented

- 2:1 hierarchy, balancing, geometry/wake refinement, hash lookup.
- Finest refinement is reserved for fabric; the configured volumetric near wake stops one level coarser.
- Halo, restriction, prolongation, and aperture-aware reflux helpers.
- Sparse cross-brick EB atlas with finest-owner selection.
- One coupled matrix-free pressure graph spanning every level.
- Compact conservative 2:1 face tiles and EB aperture connections.
- CUDA operator/composite PCG plus manufactured coupled projection tests.
- Persistent GPU divergence/RHS/correction bridge over pooled fields and compact special-flux arrays.
- Additive two-level Galerkin preconditioner with compact nonlocal aggregate connections.
- Deterministic gauges for active pressure components disconnected from the outlet.

Pending:

- aperture-aware cross-level EB topology when fabric reaches a 2:1 interface;
- recursive multigrid/V-cycle acceleration beyond the present two-level preconditioner;
- AMR-versus-uniform flow validation.

## F. FP32 production — storage/operator foundation implemented

- New production scalar is FP32 with an optional FP64 validation build.
- Geometry preprocessing and force/reference calculations retain FP64 where valuable.
- CUDA pressure operator has measured FP32/FP64 parity and timing probes.

The first-order AMR advection/Smagorinsky/external-BC/composite-projection path is GPU resident. Pending: higher-order fabric-side-safe reconstruction, consistent brick/coarse-fine velocity interpolation, and full-case convergence validation.

## G. Aerodynamics — manufactured load path implemented

- Per-triangle plus/minus pressure and Cp, delta-Cp, force, whole-wing force/moment.
- Winding-invariant pressure force.
- CL/CD/CS only with explicit positive reference area.
- Pressure drag is labelled pressure-only.

Composite patches now retain global plus/minus pressure DOFs and publish triangle/whole-wing pressure-only loads. Pending: establish a converged validated paraglider flow and add automatic surface field transfer/rendering.

## H. GUI migration — preview foundation implemented

The STEP viewer and CUDA/OpenGL infrastructure remain useful. Particle/tracer fabric collision uses the triangle BVH. The viewer displays static AMR/owned-EB debug boxes and reports the composite pressure readiness/counts. It refuses to run the legacy timestep for a paraglider. Primary controls, scene schema, per-side surface colouring, and the aerodynamic dashboard still require migration.

## I. Obsolete-code removal — intentionally pending

Remove the building, centerline, roof/wall, solid-voxel-load, channel/ground/seabed/porous, building probe/config, and corresponding UI/doc paths only after the new timestep and UI replacements pass their validation cases. Until then they are reference code, not current architecture.

## Definition of trustworthy version 1

- No pressure/velocity connection crosses fabric.
- Openings remain fluid-connected.
- Composite AMR flux is conservative through EB and coarse/fine interfaces.
- GPU projection reaches the declared residual tolerance with a documented norm.
- FP32 agrees with CPU/FP64 reference cases at a justified tolerance.
- Normal, parallel, and inclined dynamic plates pass; the current equal-finest AMR inclined-plate lift differs from uniform-fine by about 15% (25% initial gate).
- A missing cavity inlet reconnects the internal pressure graph to external air; developed internal mass-flow validation remains pending.
- End-to-end timestep contains no per-step CPU geometry work or bulk field transfers.
- GUI and logs clearly distinguish pressure drag from total drag and expose unresolved geometry.
