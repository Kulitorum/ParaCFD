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

The global timestep is set by the configured CFL, finest active spacing, and a GPU maximum reduction over current regular AMR velocity. There is no time subcycling in version 1. Compact EB graph transport separately bounds its local Courant and explicit molecular-diffusion contribution; microscopic apertures therefore remain visible diagnostics without forcing an irrelevant domain-wide timestep.

## Advection near fabric

Ordinary MacCormack/semi-Lagrangian sampling is invalid when a backtrace or interpolation stencil crosses fabric because a geometrically nearby point may be in a disconnected pressure/velocity region. Static preprocessing marks every MAC face centre within 2.5 local cells of fabric and tests its six one-cell Cartesian segments with the triangle BVH. Bits 0..5 of one byte store the open same-side links and bit 6 marks the fabric band. Marked faces use a bounded minmod-limited linear reconstruction when their upstream, second-upstream, and downstream segments are all open on the same side; missing sheet-edge, boundary, or level-transition links revert only that local contribution to the first-order donor update. Diffusion is also restricted to open links. All other faces use bounded MacCormack advection with RK2 characteristics, the finest-brick hash, and Smagorinsky diffusion. The global displacement is bounded below one cell, so the unmarked far-field stencil cannot reach the sheet. This prevents cross-fabric donor data without per-step triangle traversal, but the band still needs an arbitrary-orientation wall/LES treatment. Same-level stencils cross brick boundaries instead of clamping, and every level's forward state is complete before the MacCormack correction reads it. At a 2:1 interface, compact normal-flux tiles gather from transported fine MAC faces and scatter projected values back to fine faces plus an aperture-weighted coarse mean. Irregular aperture velocities advance on a separate compact fluid-connectivity graph: adjacent implicit MAC faces and aperture states form the local fragment velocity, followed by first-order upwind transport, molecular diffusion, and a graph-normal Smagorinsky estimate. Neither graph contains a fabric crossing or imposes irregular indirection on far-field cells. General cross-level interpolation, higher-order compact-graph reconstruction, and full-tensor irregular LES remain force-convergence blockers.

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

Manufactured tests cover flat/inclined/curved membranes, winding reversal, a deliberate gap, a four-region T-junction, a terminating sheet, a thin trailing edge, disconnected apertures on a partially covered aligned face, conservative removal of a numerical sub-grid aperture from both pressure and velocity graphs, same-side small-fragment conservation, disconnected-component gauges, cross-brick EB projection, coupled 2:1 projection, AMR constant-state advection/LES, opposite-side rejection and tangential transport on the precomputed fabric link graph, physical-face-only CFL reduction, independent compact-aperture transport, and FP32/FP64 operator parity. Dynamic gates give downstream pressure force for a normal plate, zero pressure drag for a parallel plate, and positive lift/drag signs at +15 degrees. In the current short equal-finest 0.125 m case, AMR inclined-plate lift differs from uniform-fine by about 0.4%. FP32 and FP64 dynamic plate forces agree within about 3e-5 relative; maximum post-projection divergence is about 1e-6 in the FP32 cases and 1e-12 in FP64. After 12 steps, a fabric box with its upstream sheet removed has about 0.0382 m^3/s represented positive inlet flux and 0.0763 m^3/s absolute bidirectional exchange; the otherwise identical closed box remains exactly zero. This validates transport through a real opening, not pressure equilibration or a converged ram-air inlet flow.

The checked-in PlanB wing builds 4.478 million composite pressure states with zero unresolved owned topology cells. GUI and probes canonicalize placement to one FP32-representable, base-brick-aligned frame before clipping; this removes topology changes caused by sub-micrometre display/CLI placement roundoff. At the tightened 1e-5 projection tolerance, the initialized freestream solve takes about 180 iterations and roughly 127 ms on the RTX 4090. With conservative same-side merging below 0.25 cell volume and removal/reporting below `1e-4 h^2`, 2,796 numerical slivers totalling 0.000283775 m^2 of cumulative level-atlas area are removed; a resolved opening remains ordinary fluid connectivity. A confirmed-upstream-orientation 500-step adaptive-CFL diagnostic with limited fabric-band transport and graph-normal EB LES reaches 0.531 s. The regular peak is 11.99 m/s and the compact peak is 46.8 m/s. Pressure-only force is still evolving at `[145.9, 0.10, 10.0]` N, max/volume-weighted-RMS divergence is `4.88e-4 / 3.93e-6 s^-1`, net integrated flux error is `-3.15e-5 m^3/s`, and the last step takes 73.0 ms including 37.2 ms projection and 0.13 ms compact EB transport. Persistent storage remains about 493.57 MiB because the same-side mask is one byte per face. The regular and background field no longer exhibits the former structured leading-edge runaway, but the compact peak is still too high for trusted aerodynamic loads and remains a documented discretization/convergence blocker. Reducing the merge threshold to 0.05 left a 65.7 m/s compact peak after 500 steps and advanced only 0.288 s in the earlier orientation experiment, so 0.25 remains the default pending a less geometry-diffusive stabilization method.

Resolution diagnostics can override only the level count through `paraglider_case_probe --max-levels N`. At two levels (0.125 m finest), the case has 1.316 million pressure states and uses 148.85 MiB. At four levels (0.03125 m finest), it has 12.938 million states and uses 1.427 GiB. A CAD-provenance/flux trace shows the previous four-level 102 m/s maximum connected two conservative fluid fragments 10-14 mm from different aft CAD faces; the high-speed aperture carried about 0.0125 m^3/s and both fragment net fluxes were about 1e-9 m^3/s. One merged same-side fragment had 39 incident apertures. This identifies a complex trailing-edge/seam passage, not fabric penetration or global divergence.

Compact transport now adds `nu_t = (Cs d)^2 |du_n/dn|` from same-side graph endpoint states. It preserves uniform flow and reduces manufactured perturbation energy. With `Cs=0.1`, the default three-level case reaches 0.677 s in 500 steps with regular/compact peaks of 13.1/39.2 m/s and force `[131.9, -0.52, 39.4]` N. The four-level case reaches 0.294 s with peaks 14.7/48.7 m/s and force `[157.4, 0.02, 40.3]` N; before graph LES it reached only 0.202 s with a 78.8 m/s compact peak and later approached 102 m/s. The samples are still at unequal physical times, and comparable-time forces remain roughly 12-14% apart. Thus graph LES fixes the compact growth/stiffness but does not establish force convergence. Higher-order same-side transport, full-tensor EB LES/wall treatment, and an explicit local velocity/flux acceptance criterion remain mandatory.
