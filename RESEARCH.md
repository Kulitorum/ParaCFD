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

Ordinary MacCormack/semi-Lagrangian sampling is invalid when a backtrace crosses fabric because a geometrically nearby point may be in a disconnected pressure/velocity region. The initial GPU implementation precomputes a conservative protection mask for every MAC face centre within 2.5 local cells of fabric. Protected values do not backtrace; all other traces use the finest-brick hash. With the current global timestep, an unprotected trace cannot reach the sheet. This prevents cross-fabric sampling without per-step triangle traversal, but it is a first-order fallback: the final scheme still needs a higher-order same-side reconstruction or trace clipping. Current local interpolation and diffusion also clamp at brick edges and therefore require a conservative same-level/coarse-fine reconstruction before force convergence studies.

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

Manufactured tests cover flat/inclined/curved membranes, winding reversal, a deliberate gap, a four-region T-junction, a terminating sheet, a thin trailing edge, disconnected apertures on a partially covered aligned face, same-side small-fragment conservation, disconnected-component gauges, cross-brick EB projection, coupled 2:1 projection, AMR constant-state advection/LES, fabric-protected backtraces, and FP32/FP64 operator parity. The checked-in PlanB wing builds 4.48 million composite pressure states with zero unresolved owned topology cells; an initialized freestream projection reaches about a 9.4e-5 global relative residual in 168 iterations on the RTX 4090. The first following timestep costs about 1.0 ms advection, 1.3 ms LES/diffusion, and 150 ms projection. Its local maximum divergence and transient pressure-only force are not yet acceptance-quality, so these are architecture/performance evidence rather than aerodynamic validation. Normal/parallel/inclined dynamic plate cases, an opened internal cavity, and AMR-versus-uniform force/convergence studies remain required.
