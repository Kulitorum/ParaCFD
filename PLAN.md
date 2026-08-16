# Paraglider CFD implementation plan

This file records actual migration status, not a claim that unfinished architecture exists.

## A. Geometry foundation — implemented

- Preserve OpenCascade STEP import and metre conversion.
- Retain source STEP face ID and optional UV data per tessellated triangle.
- Build a static triangle BVH with AABB, segment, nearest-point, and distance queries.
- Add paraglider configuration and external-aero data structures.
- Retain only mathematically reusable FP64 MAC kernels as immutable CPU/GPU validation references.

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

The composite CPU reference and persistent GPU projection now cover EB divergence, RHS, pressure, regular/special-flux correction, and post-projection divergence. Disconnected active fluid components receive deterministic gauges. `ExternalAeroCore` integrates this projection into a persistent GPU timestep.

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
- longer-run force/conservation convergence beyond the current short-case comparison.

## F. FP32 production — storage/operator foundation implemented

- New production scalar is FP32 with an optional FP64 validation build.
- Geometry preprocessing and force/reference calculations retain FP64 where valuable.
- CUDA pressure operator has measured FP32/FP64 parity and timing probes.

The bounded MacCormack/RK2 AMR advection, Smagorinsky, external-BC, and composite-projection path is GPU resident. Same-level interpolation, Smagorinsky gradients, eddy-viscosity access, and diffusion cross brick boundaries without stencil clamping; all levels build consistent forward states before correction and commit. Normal 2:1 flux tiles gather from transported fine MAC faces, then scatter projected values to the fine faces and their aperture-weighted coarse mean. Compact EB aperture velocities advance on their own same-side graph using structured carrier states, first-order upwinding, molecular diffusion, and a graph-normal Smagorinsky estimate. One global timestep is selected from a GPU maximum reduction over physical active-brick faces; ghost storage and covered bricks are excluded, and compact graph updates enforce their local bound. Conservative same-side 0.25-volume merging plus a configurable/reportable `1e-4 h^2` numerical aperture cutoff remove compact numerical slivers. The static fabric band stores a compact six-link same-side graph and uses bounded first-order donor transport plus link-restricted molecular diffusion. PlanB diagnostics show EB LES reduces the four-level aft-junction peak from 102 to about 49 m/s without changing conservation, but equal-time pressure loads remain roughly 12-14% grid-dependent. Pending: higher-order same-side reconstruction in both graph paths, full-tensor EB LES/wall treatment, a local compact-flux acceptance gate, consistent general cross-level interpolation, and full grid/domain/force convergence validation.

## G. Aerodynamics — manufactured load path implemented

- Per-triangle plus/minus pressure and Cp, delta-Cp, force, whole-wing force/moment.
- Winding-invariant pressure force.
- CL/CD/CS only with explicit positive reference area.
- Pressure drag is labelled pressure-only.

Composite patches now retain global plus/minus pressure DOFs and publish triangle/whole-wing pressure-only loads. The GUI worker automatically transfers throttled per-triangle Cp+, Cp-, delta-Cp, and pressure-force summaries. Visible-side and explicit plus/minus surface rendering are implemented. Pending: establish a converged validated paraglider flow.

## H. GUI migration — implemented foundation

The executable now starts a paraglider-only Qt window. Particle/tracer fabric collision uses the triangle BVH. The viewer displays live velocity/pressure slices, arrows, tracers, static AMR/owned-EB debug boxes, and two-sided triangle pressure colours. Direct STEP load infers horizontal span/chord from the face-only bbox, rotates the assumed negative chord direction into the fixed +X freestream, then sizes and centres the domain around the placed wing; LE/TE polarity remains an explicit user confirmation. A shared FP32-canonical zero-origin framing helper makes GUI and probe EB clipping coordinates identical. The default close view and explicit `Fit Wing`/`Fit Domain` controls keep visual framing separate from the physically asymmetric wake domain. Start/Play/Step advance only `ExternalAeroCore`. The panel reports pressure convergence/timing, pressure-only loads, regular/EB velocity, conservation, and persistent-memory diagnostics. Pending: direct AMR CUDA/OpenGL sampling and pressure-force-vector display.

## I. Obsolete-code removal — implemented

The building, centerline, roof/wall, solid voxelizer/load, channel/ground/seabed/porous, old scene/worker/window, building probe/config, old installer, and misleading legacy research paths have been deleted after the replacement release and tests passed. Git history remains the reference if an old numerical idea must be recovered.

## Definition of trustworthy version 1

- No pressure/velocity connection crosses fabric.
- Openings remain fluid-connected.
- Composite AMR flux is conservative through EB and coarse/fine interfaces.
- GPU projection reaches the declared residual tolerance with a documented norm.
- FP32 agrees with CPU/FP64 reference cases at a justified tolerance.
- Normal, parallel, and inclined dynamic plates pass; the current short equal-finest AMR inclined-plate lift differs from uniform-fine by about 0.4%.
- A missing cavity inlet admits developed bidirectional flow while its closed control has zero represented opening flux.
- End-to-end timestep contains no per-step CPU geometry work or bulk field transfers.
- GUI and logs clearly distinguish pressure drag from total drag and expose unresolved geometry.
