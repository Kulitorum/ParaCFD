# Numerical basis for the paraglider CFD path

## Geometry and topology

The CAD model is a collection of oriented surface patches, not a watertight solid. A triangle normal defines a local minus-to-plus direction. Each resolved fabric cut partitions a Cartesian control volume into fluid fragments. Pressure belongs to a fragment, not necessarily to a Cartesian cell, because distinct fluid regions within one cell must retain independent states.

An opening is represented by the absence of fabric. No parity test or inside/outside classification is valid for this model. Surface preprocessing therefore asks only which triangles intersect a cell/face and how they partition local fluid connectivity.

Regular cells use implicit Cartesian topology. Irregular cells carry compact fragment, aperture, blocked-face, and fabric-patch arrays. A sheet coincident with a Cartesian face suppresses that implicit face connection while retaining two ordinary full-volume cells. If multiple non-coplanar sheets or unresolved sheet endings make local topology ambiguous, refinement is mandatory. Reaching the finest level without resolving the topology is a grid-generation error, never permission to connect or merge regions.

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

## Small control volumes

Tiny cut fragments create severe explicit stability limits and poorly conditioned pressure equations. The initial stabilization merges a fragment below a configurable volume fraction with a larger fluid neighbour through the largest appropriate open aperture. A fabric surface is not an open aperture, so the operation cannot cross sides. Strictly increasing merge chains are resolved to one root, and volumes/centroids/source fluxes are accumulated conservatively into that representative.

## AMR layout

Refinement is static for a run and always ratio two. Cells within a brick are uniform. Bricks are located using integer coordinates in per-level hash tables. Field allocations are contiguous per level and per component; normal CUDA kernels use structured brick indexing, while only compact EB work lists use irregular indirection.

The initial global timestep is set by the finest active cell spacing. There is no time subcycling in version 1.

## Advection near fabric

Ordinary MacCormack/semi-Lagrangian sampling is invalid when a backtrace crosses fabric because a geometrically nearby point may be in a disconnected pressure/velocity region. The AMR port must test the trace segment against static fabric/EB topology and either clip/project the trace on its current side or apply the wall state. A point-distance test alone is insufficient.

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

Manufactured tests cover a flat and inclined membrane, winding reversal, a deliberate gap, complex T/trailing-edge detection, same-side small-fragment conservation, one-level projection, ratio-two exchange conservation, and FP32 operator parity. They do not yet establish end-to-end paraglider accuracy. Normal/parallel/inclined dynamic plate cases, an opened internal cavity, composite-AMR flow conservation, and AMR-versus-uniform force/convergence studies are still required.
