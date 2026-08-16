# Numerical basis for the paraglider CFD path

## Geometry and topology

The CAD model is a collection of oriented surface patches, not a watertight solid. A triangle normal defines a local minus-to-plus direction. Each resolved fabric cut partitions a Cartesian control volume into fluid fragments. Pressure belongs to a fragment, not necessarily to a Cartesian cell, because distinct fluid regions within one cell must retain independent states.

An opening is represented by the absence of fabric. No parity test or inside/outside classification is valid for this model. Surface preprocessing therefore asks only which triangles intersect a cell/face and how they partition local fluid connectivity.

Regular cells use implicit Cartesian topology. Irregular cells carry compact fragment, aperture, blocked-face, and fabric-patch arrays. A sheet coincident with a Cartesian face suppresses that implicit face connection while retaining two ordinary full-volume cells; partial coverage produces one record per connected fluid opening. At the finest level, a deterministic local microcell graph may recover multi-fragment topology at sheet endings and junctions. All microcells remain fluid and only exact triangle crossings remove graph edges, so this is neither solid voxelization nor parity filling. Reaching the finest level without a valid topology is a grid-generation error, never permission to connect or merge regions.

## Finite-volume pressure operator

For two connected fluid control volumes `i` and `j`, an open aperture of area `A` and representative centre distance `d` contributes the symmetric conductance

```
g_ij = A / d
```

to the integrated pressure operator. The volume-normalized Laplacian is

```
(L p)_i = (1 / V_i) sum_j g_ij (p_i - p_j)
```

The divergence and velocity correction must use the same oriented aperture fluxes. This compatibility is necessary for conservation. Fabric patches do not create graph edges; their prescribed normal mass flux is zero for a stationary impermeable sheet.

At a 2:1 AMR interface, the geometric coarse aperture flux is replaced by the sum of the corresponding fine-aperture fluxes. The composite operator must enforce this during every matrix-vector product and correction, not by projecting each level independently and reconciling fields afterwards.

The CUDA solve is matrix-free PCG. Its current additive preconditioner applies the fine operator diagonal and restricts residuals geometrically into level-zero cell aggregates. The aggregate operator is Galerkin: every regular, coarse/fine, and EB conductance contributes to its diagonal and edge coefficients. Conservative fragment merges can produce aggregate edges that are not immediate Cartesian neighbours, so those edges remain in a compact nonlocal list rather than being discarded. Components that cannot reach the X-max pressure outlet receive one deterministic positive gauge coefficient; this removes each component's constant nullspace without opening a path through fabric.

## Small control volumes

Tiny cut fragments create severe explicit stability limits and poorly conditioned pressure equations. The initial stabilization joins adjacent tiny fragments only through actual apertures, then merges the group into an appropriate larger fluid neighbour. A fabric surface is not an open aperture, so the operation cannot cross sides. Merge chains resolve to one root, and volumes are accumulated conservatively into that representative. A group with no external fluid aperture remains an explicitly counted pressure-static pocket rather than being connected through fabric.

## AMR layout

Refinement is static for a run and always ratio two. Cells within a brick are uniform. Bricks are located using integer coordinates in per-level hash tables. Field allocations are contiguous per level and per component; normal CUDA kernels use structured brick indexing, while only compact EB work lists use irregular indirection.

The initial global timestep is set by the finest active cell spacing. There is no time subcycling in version 1.

## Advection near fabric

Ordinary MacCormack/semi-Lagrangian sampling is invalid when a backtrace crosses fabric because a geometrically nearby point may be in a disconnected pressure/velocity region. The GPU implementation precomputes a conservative protection mask for every MAC face centre within 2.5 local cells of fabric. Protected values do not backtrace; all other faces use bounded MacCormack advection with RK2 characteristics and the finest-brick hash. With the current global timestep, an unprotected trace cannot reach the sheet. This prevents cross-fabric sampling without per-step triangle traversal, but the protected band is a first-order fallback: the final scheme still needs a higher-order same-side reconstruction or trace clipping. Trilinear, Smagorinsky-gradient, eddy-viscosity, and diffusion stencils fetch across same-level brick boundaries instead of clamping, and every level's forward state is complete before the MacCormack correction reads it. Coarse/fine face reconstruction is still interpolatory rather than conservative and remains a force-convergence blocker.

## Surface forces

For a patch with normal `n` from minus to plus,

```
delta_p = p_minus - p_plus
F_patch = delta_p A n
```

Reversing triangle winding changes `n` to `-n` and swaps plus/minus, leaving `F_patch` unchanged. Source-triangle values are area-weighted sums of their embedded patches. The pressure coefficient is based on the configured freestream dynamic pressure

```
q_inf = 0.5 rho U_inf^2
Cp = (p - p_inf) / q_inf
```

Whole-wing force/moment accumulation uses double precision. Reference coefficients are undefined unless a reference area is explicitly supplied. Until tangential wall shear is implemented, the force result is pressure force and pressure drag only.

## Precision policy

Production velocity, pressure, turbulent viscosity, and temporary GPU fields use FP32 for an Ada/RTX 4090 target. Geometry clipping/BVH construction uses CPU FP64. Validation compares the GPU FP32 operator with analytical or CPU FP64 results; tolerances must follow measured conditioning and convergence rather than being loosened to accept changed output.

## Current evidence and remaining validation

Manufactured tests cover flat/inclined/curved membranes, winding reversal, a deliberate gap, a four-region T-junction, a terminating sheet, a thin trailing edge, disconnected apertures on a partially covered aligned face, same-side small-fragment conservation, disconnected-component gauges, cross-brick EB projection, coupled 2:1 projection, AMR constant-state advection/LES, fabric-protected backtraces, and FP32/FP64 operator parity. Dynamic gates now give downstream pressure force for a normal plate, zero pressure drag for a parallel plate, and positive lift/drag signs at +15 degrees. At equal 0.125 m finest spacing the initial AMR inclined-plate lift differs from uniform-fine by about 15%, while the measured small-case step drops from about 18.6 ms to 10.6 ms. FP32 and FP64 dynamic plate forces agree within about 3e-5 relative; maximum post-projection divergence is about 1e-6 in the FP32 cases and 1e-12 in FP64. Removing the upstream face of a closed fabric box reconnects its internal pressure component to external air, although a developed internal mass-flow measurement still requires full transport of compact aperture velocities.

The checked-in PlanB wing builds 4.48 million composite pressure states with zero unresolved owned topology cells; an initialized freestream projection reaches about a 9.4e-5 global relative residual in 168 iterations on the RTX 4090. With bounded MacCormack enabled, the first following timestep costs about 24.9 ms advection, 7.8 ms LES/diffusion, and 126 ms projection (projection timing varies). Persistent storage is about 488.72 MiB. The transport cost and extra forward-state buffer are explicit optimization targets. Its local maximum divergence and transient pressure-only force are not yet acceptance-quality, so these are architecture/performance evidence rather than aerodynamic validation.
