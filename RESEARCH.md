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

The global timestep is set by the configured CFL, finest active spacing, and a GPU maximum reduction over current regular AMR velocity. There is no time subcycling in version 1. Compact EB graph transport separately bounds its local Courant and explicit molecular-diffusion contribution. Precomputed nearest same-axis graph links provide upstream, second-upstream, and downstream aperture states for a minmod-limited MUSCL update when the complete directed chain exists; incomplete or junction chains use the donor fallback. The combined explicit update is also restricted to its actual same-side stencil extrema because independently limited advection and diffusion coefficients need not sum to a convex update. This is a monotone discretization rule, not an absolute speed cap; pressure projection can still create a new compact velocity extremum, which remains visible to diagnostics and the following timestep selection.

## Advection near fabric

Ordinary MacCormack/semi-Lagrangian sampling is invalid when a backtrace or interpolation stencil crosses fabric because a geometrically nearby point may be in a disconnected pressure/velocity region. Static preprocessing marks every MAC face centre within 2.5 local cells of fabric and tests its six one-cell Cartesian segments with the triangle BVH. Bits 0..5 of one byte store the open same-side links and bit 6 marks the fabric band. Marked faces use a bounded minmod-limited linear reconstruction when their upstream, second-upstream, and downstream segments are all open on the same side; missing sheet-edge, boundary, or level-transition links revert only that local contribution to the first-order donor update. Diffusion is also restricted to open links. All other faces use bounded MacCormack advection with RK2 characteristics, the finest-brick hash, and Smagorinsky diffusion. The global displacement is bounded below one cell, so the unmarked far-field stencil cannot reach the sheet. This prevents cross-fabric donor data without per-step triangle traversal, but the band still needs an arbitrary-orientation wall/LES treatment. Same-level stencils cross brick boundaries instead of clamping, and every level's forward state is complete before the MacCormack correction reads it. A tangential staggered lattice ends half a cell short of a coarse brick boundary; using the nearest coarse value there produced a 0.0266 m/s defect in an exactly linear characteristic. Cross-level velocity and eddy-viscosity samples now use a one-sided linear extension, and reverse MacCormack correction is disabled locally when its trace changes lattice. The same 1,728-interface manufactured case now has about `4.8e-8` m/s FP32/FP64-coordinate error after both advection and linear diffusion. At a 2:1 interface, compact normal-flux tiles gather from transported fine MAC faces and scatter projected values back to fine faces plus an aperture-weighted coarse mean. Irregular aperture velocities advance on a separate compact fluid-connectivity graph: adjacent implicit MAC faces and aperture states form the local fragment velocity, complete same-axis aperture chains use bounded MUSCL transport, and chain endings/junctions use the donor fallback. Connected edge differences reconstruct all nine velocity-gradient components at graph nodes; the symmetric tensor magnitude supplies compact Smagorinsky dissipation. Neither graph contains a fabric crossing or imposes irregular indirection on far-field cells. Momentum advection remains semi-Lagrangian rather than globally conservative/refluxed; multidimensional least-squares junction reconstruction and an arbitrary-orientation wall treatment remain force-convergence blockers.

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

The current 47,997-triangle PlanB working fixture builds 4.478 million composite pressure states with zero unresolved owned topology cells. GUI and probes canonicalize placement to one FP32-representable, base-brick-aligned frame before clipping; this removes topology changes caused by sub-micrometre display/CLI placement roundoff. At the tightened 1e-5 projection tolerance, the initialized freestream solve takes about 150 iterations and roughly 100 ms on the RTX 4090. With conservative same-side merging below 0.25 cell volume and removal/reporting below `1e-4 h^2`, 2,918 numerical slivers totalling 0.000299041 m^2 of cumulative level-atlas area are removed; a resolved opening remains ordinary fluid connectivity. The reconstructed-strain 500-step run with confirmed upstream orientation and bounded fabric-band/compact transport reaches 0.527 s. Regular/compact peaks are 12.76/45.45 m/s and pressure-only force is still evolving at `[174.66, 0.02, -5.27]` N. Max/volume-weighted-RMS divergence is `1.20e-3 / 4.18e-6 s^-1`, net integrated flux error is `3.89e-6 m^3/s`, and the last measured step takes 55.8 ms including 14.2 ms projection and 0.28 ms compact EB transport. Persistent storage is 496.98 MiB including directed chain links and node gradient accumulators; the same-side structured mask remains one byte per face. The preceding graph-normal MUSCL build reached 55.05 m/s after 500 steps, making the reconstructed-strain peak about 17% lower. The remaining compact peak and grid-dependent loads are still too high for trustworthy aerodynamics. The PlanB STEP has user-owned working-tree edits and is not part of the solver milestone commit.

An equal-time gate now removes adaptive-step count as a convergence confounder. At `t~=0.100 s`, 2/3/4 levels with 0.125/0.0625/0.03125 m finest spacing require 41/89/197 steps and give pressure-only `[Fx,Fz]` of `[180.2,-8.17]`, `[200.6,-5.25]`, and `[203.3,-24.77]` N. Persistent estimates are 149.50/496.98/1442.93 MiB. Streamwise force changes 11.3% from levels 2 to 3 and 1.4% from 3 to 4, while vertical force is strongly non-monotone. The flow is also still in its early transient. This is a reproducible convergence failure baseline for conservative momentum transport and longer-time studies, not evidence that Fx is converged.

The first conservative-momentum kernel advances one uniform MAC level with first-order upwind finite-volume fluxes. Each lattice link computes the same oriented momentum flux from either endpoint, including across a same-level brick boundary; a static fabric link removes the pair from both endpoints. FP32 and FP64 manufactured tests verify exact constant-state preservation, compact-support discrete momentum conservation, duplicate normal-face agreement at a brick boundary, and zero transport onto the opposite side of a membrane. This is an incremental operator test, not the production solver: the kernel has not yet incorporated compact fragment/aperture momentum or the staggered geometry of a 2:1 flux register, and `ExternalAeroCore` continues to use bounded semi-Lagrangian transport.

Resolution diagnostics can override the level count through `paraglider_case_probe --max-levels N` and stop at a common time with `--physical-time T`. At two levels (0.125 m finest), the case has 1.316 million pressure states and uses about 149.5 MiB. At four levels (0.03125 m finest), it has 12.937 million states and uses 1.409 GiB. A CAD-provenance/flux trace shows the previous four-level 102 m/s maximum connected two conservative fluid fragments 10-14 mm from different aft CAD faces; the high-speed aperture carried about 0.0125 m^3/s and both fragment net fluxes were about 1e-9 m^3/s. One merged same-side fragment had 39 incident apertures. This identifies a complex trailing-edge/seam passage, not fabric penetration or global divergence.

Compact transport now reconstructs a row-major `du_i/dx_j` tensor from area-weighted same-side graph endpoint differences and adds `nu_t = (Cs d)^2 sqrt(2 Sij Sij)`. The normal-rate estimate remains a conservative fallback where the graph lacks derivative support. Manufactured checks cover pure tangential shear, rigid rotation, normal extension, uniform flow, perturbation-energy dissipation, and same-side stencil bounds. A separate quadratic-characteristic test exercises the MUSCL formula, and an eight-cell membrane fixture proves that complete compact chains exist and preserve a uniform three-component state on GPU. The current three-level comparison reduces the 500-step compact peak from 55.05 to 45.45 m/s relative to graph-normal LES, but this does not establish force convergence. Multidimensional least-squares junction reconstruction, conservative/refluxed momentum transport, and a fabric wall model remain mandatory.
