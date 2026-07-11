# Stable Fluids in 3D: GPU Pressure Solvers and Accuracy Fixes

## Overview

Stam's "Stable Fluids" (SIGGRAPH 1999) solves incompressible Navier-Stokes with operator splitting: (1) add body forces, (2) advect velocity semi-Lagrangianly, (3) diffuse implicitly, (4) project to divergence-free via a pressure Poisson solve. Every step is unconditionally stable, so the time step is limited only by accuracy, not stability. The price is heavy first-order numerical diffusion. "Real-Time Fluid Dynamics for Games" (Stam, GDC 2003) is the same algorithm with a Gauss-Seidel linear solver and collocated grid; it is a prototype, not a production spec.

For our engineering use (bed shear stress around 3D-printed shapes, 10 x 10 x 5 m domain, cell size h = 0.025-0.05 m, i.e. 200x200x100 = 4M to 400x400x200 = 32M cells, currents 0.5-2.5 m/s, rho = 1025 kg/m3, nu = 1.05e-6 m2/s), three fixes are mandatory relative to the 1999/2003 papers: a staggered MAC grid, MacCormack (BFECC-style) advection with a limiter, and a multigrid-preconditioned CG pressure solver. Vorticity confinement is a graphics device we should NOT use for quantitative work. Smagorinsky LES replaces the implicit diffusion step (molecular viscosity is negligible at this resolution).

## Key equations & results

### 1. Operator splitting and semi-Lagrangian advection (Stam 1999)

Per step: u* = advect(u^n, dt); u** = u* + dt*f; then project u** -> u^{n+1}.
Semi-Lagrangian advection of any quantity q (velocity component, sediment concentration):

q^{n+1}(x) = q^n(x_back),  x_back = x - dt * u(x)   (trilinear interpolation at x_back)

First order in space and time; unconditionally stable because interpolation never creates new extrema. Use a 2nd-order Runge-Kutta backtrace instead of the single Euler step: x_mid = x - 0.5*dt*u(x); x_back = x - dt*u(x_mid). Backtraces that land inside a solid voxel are clipped to the solid face (Fedkiw, Stam & Jensen 2001, Fig. 2). Practical CFL (= u*dt/h): accuracy degrades above CFL ~ 1-5 even though stability is unconditional. At u = 2.5 m/s, h = 0.05 m: dt = 0.02 s at CFL 1, dt = 0.1 s at CFL 5.

Effective numerical viscosity of 1st-order semi-Lagrangian/upwind advection is O(0.5 * |u| * h * (1 - CFL)) ~ 0.03-0.06 m2/s for u = 2.5 m/s, h = 0.05 m — 4-5 orders of magnitude above seawater molecular viscosity (1.05e-6 m2/s) and 10-100x a typical Smagorinsky eddy viscosity. This is THE accuracy problem of this solver family and why step 2 below is required.

### 2. MacCormack advection with clamping (Selle, Fedkiw, Kim, Liu & Rossignac, J. Sci. Comput. 35:350-371, 2008)

Let A = first-order semi-Lagrangian advection forward in time, A_R = the same operator run with reversed velocity (backward in time):

phi_hat^{n+1} = A(phi^n)                    (forward step)
phi_hat^{n}   = A_R(phi_hat^{n+1})          (backward step)
e             = (phi_hat^{n} - phi^n) / 2   (error estimate)
phi^{n+1}     = phi_hat^{n+1} - e  =  phi_hat^{n+1} + (phi^n - phi_hat^n)/2

Two semi-Lagrangian sweeps total (BFECC needs three: it corrects phi^n and re-advects, phi^{n+1} = A(phi^n - e)). Both are 2nd-order accurate in space AND time; MacCormack is equivalent to RK2 for ODEs. Measured convergence order ~2.0 on Zalesak's disc at CFL 0.75 and 1.75 (Selle et al., Tables/Figs 4-6).

Stability limiter (required — "without any limiter the simulation becomes unstable"): clamp phi^{n+1} to [min, max] of the 8 corner values used in the trilinear interpolation of the FIRST advection step; or, preferred by Selle et al., revert to the plain first-order result phi_hat^{n+1} wherever the corrected value exceeds those bounds. Also revert to first-order semi-Lagrangian at any cell whose backtrace characteristic could sample a solid or domain boundary (i.e. within CFL*h of a boundary), because forward/backward rays mixing fluid and wall data give garbage error estimates. Advect each velocity component separately and limit separately.

### 3. Helmholtz projection / pressure Poisson equation (Stam 1999; Fedkiw et al. 2001)

With constant density rho (= 1025 kg/m3), solve:

laplacian(p) = (rho / dt) * div(u*)         [p in Pa]
u^{n+1} = u* - (dt / rho) * grad(p)

Boundary conditions: pure Neumann dp/dn = 0 on solid faces and inlet; Dirichlet p = 0 on the outlet. Discrete 7-point Laplacian on the MAC grid: (sum of fluid-neighbor p) - (number of non-solid neighbors)*p_ijk, divided by h^2; a solid neighbor is simply dropped from the stencil (this IS the Neumann BC). With Neumann on ALL boundaries the system is singular and requires the compatibility condition integral(div u*) dV = 0 — enforce global mass balance by scaling outlet face velocities so outlet flux = inlet flux before every projection. Fixing p = 0 at the outlet removes the null space entirely.

### 4. MAC staggered grid vs collocated (Harlow & Welch 1965; Fedkiw et al. 2001, Appendix A)

Pressure, sediment concentration, and eddy viscosity live at cell centers; u on x-faces, v on y-faces, w on z-faces. Divergence and pressure gradient are then central differences over ONE cell width h, so adjacent pressures are directly coupled and the discrete Laplacian has no spurious null mode. On a collocated grid (Stam 1999/2003) both operators difference over 2h, decoupling the odd/even lattices — a checkerboard pressure mode is invisible to the projection and pollutes the velocity. Fedkiw, Stam & Jensen 2001 explicitly switched Stam's solver to MAC, reporting "improved results ... with less artificial dissipation." Use MAC. (Collocated can be salvaged with Rhie-Chow momentum interpolation, but on a uniform voxel grid MAC is simpler and strictly better.)

### 5. Vorticity confinement (Fedkiw, Stam & Jensen, SIGGRAPH 2001) — know it, don't use it

omega = curl(u);  N = grad(|omega|) / |grad(|omega|)|;  f_conf = epsilon * h * (N x omega)

epsilon > 0 (dimensionless O(0.1-0.5) in graphics practice; the paper gives no number) adds back rotational energy killed by numerical diffusion; the h factor makes the term vanish under refinement. It injects energy at unphysical locations/rates and is not a turbulence model. For our quantitative sediment work: default epsilon = 0 (OFF). MacCormack advection + Smagorinsky LES is the physically defensible combination. If ever used for visualization runs, keep epsilon <= 0.1.

### 6. GPU pressure Poisson solvers — convergence at 4-30M cells

- Jacobi (Harris, GPU Gems ch. 38, 2004): p_new = (p_W + p_E + p_S + p_N + p_D + p_U + alpha*b) * rBeta with alpha = -h^2, rBeta = 1/6 in 3D, b = div(u*)*rho/dt. Harris recommends 40-80 iterations ("not below 20") — but that is a VISUAL criterion. Jacobi's smoothest-mode convergence factor is ~ 1 - pi^2/(2 N^2) for N cells across; for N = 400 a 1e-4 residual reduction needs O(10^5) iterations. 40-80 Jacobi sweeps leave large-scale divergence in the flow around a 1-m obstacle. Use Jacobi only as a bring-up/debug solver.
- Red-black Gauss-Seidel / SOR: two half-sweeps (red cells, then black) are fully data-parallel on GPU. With optimal over-relaxation omega_SOR = 2 / (1 + sin(pi/N)) (N ~ 400 -> omega ~ 1.984), iteration count drops to O(N): roughly 300-600 sweeps for 1e-4 on our grids. Usable, still slow; RBGS with omega = 1 is mainly valuable as a multigrid smoother.
- Geometric multigrid V-cycle: factor-2 coarsening, damped Jacobi smoother with omega = 2/3 (stable for Poisson, fully parallel), 2-3 pre/post smoothing sweeps, direct or heavily-smoothed solve on the coarsest ~4^3 level. Plain MG on irregular voxelized domains can oscillate or diverge; McAdams, Sifakis & Teran 2010 needed 30-40 extra Gauss-Seidel sweeps (15-20 before + 15-20 after the interior sweep) in a band at least 3 cells wide around boundaries just to stabilize it.
- MGPCG (recommended; McAdams et al. 2010): one V-cycle as preconditioner inside conjugate gradient. Reported: residual drops one order of magnitude every ~2 PCG iterations, convergence rate independent of grid size and boundary irregularity; 768^2x1152 voxels in <16 GB; 128^3 solved to 1e4 residual reduction in <0.75 s on a 16-core CPU (2010) — an RTX 4090 does far better. Budget: 8-14 MGPCG iterations to ||r||/||b|| <= 1e-4 at 4-30M cells, tens of ms per solve. Incomplete-Cholesky PCG (the CPU standard, ~20 iterations for visual results in Fedkiw et al. 2001) parallelizes poorly on GPU (sparse triangular solves) — skip it.
- Tolerance guidance: ||r||/||b|| <= 1e-4 relative residual per step is adequate for morphodynamics; tighten to 1e-6 for validation runs. Warm-start CG with the previous step's pressure (halves iterations for quasi-steady tidal flow).

### 7. Solid boundaries in the projection (voxel masks; Batty, Bertails & Bridson, SIGGRAPH 2007)

Binary mask ("voxelized") approach: mark cells overlapping the STL as SOLID; set solid-face normal velocities to the solid velocity (0); drop solid neighbors from the Poisson stencil (Neumann). Free-slip (constrain only u.n = 0, leave tangential velocity untouched) is correct for us: at h = 2.5-5 cm the viscous sublayer (<1 mm) is unresolvable, so voxel no-slip (ghost tangential = -interior, as in Harris 2004) just adds fake drag and staircase dissipation. Impose bed friction through a log-law wall model on the first cell above the bed, not through no-slip.
Variational/masked projection upgrade: Batty et al. 2007 recast projection as kinetic-energy minimization; each face in the Poisson stencil and divergence is weighted by its non-solid volume (later shown slightly better with face-AREA fractions, Ng et al.). System stays symmetric positive (semi-)definite, gives sub-voxel-accurate smooth forces on curved 3D-printed shapes, and removes staircase artifacts at ~zero extra solve cost. Recommended as milestone 2; binary mask for milestone 1.

### 8. Inflow / outflow open boundaries

Inlet (upstream face): Dirichlet velocity. Use a logarithmic profile u(z) = (u_star/0.40) * ln(z/z0) fitted to the target depth-averaged current (0.5-2.5 m/s), z0 = k_s/30 with Nikuradse roughness k_s ~ 2.5*d50. Pressure: Neumann.
Outlet (downstream face): convective/radiation condition (Orlanski, J. Comput. Phys. 21:251, 1976):

du/dt + U_c * du/dn = 0,   U_c = bulk outflow speed
discrete: u_b^{n+1} = u_b^n - (U_c*dt/h) * (u_b^n - u_{b-1}^n)

applied to all velocity components on the outlet plane, followed by the global flux-balance rescale (Sec. 3) and p = 0 Dirichlet in the projection. Lateral walls and free surface (rigid lid): free-slip. For reversing tidal flow, either swap inlet/outlet roles by phase, or make both ends Orlanski + prescribed far-field velocity.

### 9. Smagorinsky LES in this solver class

nu_t = (Cs * Delta)^2 * |S|,  |S| = sqrt(2 * S_ij * S_ij),  S_ij = 0.5*(du_i/dx_j + du_j/dx_i),  Delta = (dx*dy*dz)^(1/3) = h

Cs = 0.17 from Lilly's isotropic-turbulence analysis; reduce to ~0.10 in shear flows (Deardorff); practical range 0.05-0.2; environmental-flow codes commonly use 0.10-0.12. Recommend Cs = 0.10-0.12, Delta = h. Compute |S| at cell centers from MAC face velocities; total viscosity nu_eff = 1.05e-6 + nu_t. Magnitude check: near-bed |S| ~ 50 1/s, h = 0.05 m, Cs = 0.12 -> nu_t = (0.006)^2*50 = 1.8e-3 m2/s. Explicit diffusion stability dt <= h^2/(6*nu_max) = 0.0025/(6*1.8e-3) ~ 0.23 s >> advective dt (0.02-0.1 s), so treat diffusion EXPLICITLY (skip Stam's implicit diffusion solve entirely). Clamp nu_t (e.g. <= h^2/(6*dt)) as a safety net. Note the interaction: semi-Lagrangian numerical diffusion acts as an uncontrolled extra eddy viscosity; with plain 1st-order advection it swamps nu_t, which is another reason MacCormack is mandatory before LES is meaningful.

## Practical guidance for our simulator

Per-step pipeline (velocity): MacCormack-advect u,v,w (revert to 1st-order near solids/boundaries, clamp) -> add forces (buoyancy from sediment-laden density if modeled) -> explicit Smagorinsky diffusion -> outlet Orlanski update + flux balance -> MGPCG projection to 1e-4 -> subtract (dt/rho)*grad(p). Scalars (suspended sediment): MacCormack advection with clamping (clamping also guarantees positivity), settling and diffusion handled in the sediment module.

Concrete defaults: MAC grid, h = 0.05 m (4M cells) for iteration, 0.025 m (32M) for hero runs; dt from CFL <= 2 (dt = 0.04 s at 2.5 m/s, h = 0.05 m); Cs = 0.12; vorticity confinement OFF; free-slip voxel obstacles (variational weights later); MGPCG with damped-Jacobi (omega = 2/3) V-cycle preconditioner, 2-3 smoothing sweeps, extra boundary-band Gauss-Seidel sweeps.

Honest accuracy assessment: this solver family was built for plausibility, not certified accuracy. Published quantitative uses exist — "Fast Fluid Dynamics" for room airflow (Zuo & Chen 2009: ~50x faster than CFD, "acceptable accuracy except localized disparities" vs experiments) and urban pedestrian-wind studies — but nobody validates bed shear stress from a stable-fluids solver directly. Known weaknesses: under-predicted turbulence intensity and recirculation lengths, energy loss at high CFL, staircase geometry errors, and 1st-order pressure BC at solids. Consequences for us: (a) derive bed shear from a log-law fit to the first-cell velocity, never from the raw near-wall gradient; (b) treat results as COMPARATIVE (shape A vs shape B scour ranking, equilibrium morphology trends), not absolute; (c) validate against the canonical monopile scour benchmarks (Sumer & Fredsoe) and tune Cs/wall roughness there before trusting new shapes.

## Sources

- Stam, "Stable Fluids", SIGGRAPH 1999. https://www.josstam.com/publications (PDF via ACM DOI 10.1145/311535.311548)
- Stam, "Real-Time Fluid Dynamics for Games", GDC 2003. https://www.josstam.com/publications
- Selle, Fedkiw, Kim, Liu, Rossignac, "An Unconditionally Stable MacCormack Method", J. Sci. Comput. 35:350-371, 2008. https://faculty.cc.gatech.edu/~jarek/papers/maccormack.pdf
- Fedkiw, Stam, Jensen, "Visual Simulation of Smoke", SIGGRAPH 2001. https://web.stanford.edu/class/cs237d/smoke.pdf
- Harris, "Fast Fluid Dynamics Simulation on the GPU", GPU Gems ch. 38, 2004. https://developer.nvidia.com/gpugems/gpugems/part-vi-beyond-triangles/chapter-38-fast-fluid-dynamics-simulation-gpu
- McAdams, Sifakis, Teran, "A Parallel Multigrid Poisson Solver for Fluids Simulation on Large Grids", SCA 2010, pp. 65-73. https://math.ucdavis.edu/~jteran/papers/MST10.pdf
- Batty, Bertails, Bridson, "A Fast Variational Framework for Accurate Solid-Fluid Coupling", SIGGRAPH 2007. https://www.cs.ubc.ca/labs/imager/tr/2007/Batty_VariationalFluids/
- Orlanski, "A Simple Boundary Condition for Unbounded Hyperbolic Flows", J. Comput. Phys. 21:251-269, 1976. (summary: https://www.flow3d.com/resources/cfd-101/physical-phenomena/outflow-boundary-conditions/)
- Smagorinsky model overview (Cs values, Lilly 0.17, Deardorff 0.1, filter width). https://www.sciencedirect.com/topics/engineering/smagorinsky-model
- Zuo & Chen, "Validation of Fast Fluid Dynamics for Room Airflow" / FFD improvements. https://www.researchgate.net/publication/330505570_Validation_of_Fast_Fluid_Dynamics_for_Room_Airflow
- Bridson, "Fluid Simulation for Computer Graphics" course notes (MAC grid, CFL, projection discretization). https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf

## Verification

Adversarial fact-check performed 2026-07-07 against independently fetched primary sources (paper PDFs text-extracted and searched; quotes below are verbatim from those extractions). Verdict per equation: CONFIRMED / CORRECTED / UNVERIFIABLE.

### 1. Semi-Lagrangian advection (RK2 backtrace) — CONFIRMED

- Functional form q^{n+1}(x) = q^n(x - dt*u) with trilinear interpolation, first order in space and time, unconditionally stable: confirmed in Stam, "Stable Fluids" (PDF fetched from https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/ns.pdf): "the maximum value of the new field is never larger than the largest value of the previous field", and Selle et al. 2008: the straight-line trace with trilinear interpolation "is first order accurate in space and time".
- RK2 backtrace is genuinely in Stam 1999, not a later add-on: "We use both a simple second order Runge-Kutta (RK2) method for the particle trace [14] and an adaptive particle tracer". The midpoint form x_mid = x - 0.5*dt*u(x), x_back = x - dt*u(x_mid) matches Bridson's recommended "Modified Euler, a second order Runge-Kutta (RK2) method" (fluids_notes.pdf, ch. 3).
- Clipping backtraces at solids: Fedkiw, Stam & Jensen 2001 (smoke.pdf, Fig. 2 caption): "Semi-Lagrangian paths that end up in a boundary voxel are clipped against the boundaries' face."
- CFL arithmetic checks: u = 2.5 m/s, h = 0.05 m gives dt = 0.02 s (CFL 1) to 0.1 s (CFL 5). The CFL <= ~5 accuracy heuristic is graphics practice (Foster & Fedkiw 2001; Bridson notes suggest limiting the trace to ~5 cells), not a theorem — fine as stated.

### 2. MacCormack advection with clamp — CONFIRMED

Verified against the full text of Selle, Fedkiw, Kim, Liu & Rossignac (PDF via https://faculty.cc.gatech.edu/~jarek/papers/maccormack.pdf):
- Update equations match exactly: error e = (phi_hat^n - phi^n)/2 and "phi^{n+1}_i = phi_hat^{n+1}_i - e = (phi^n_i + phi_hat^{n+1}_i - lambda*Delta+ phi_hat^{n+1})/2"; i.e. phi^{n+1} = phi_hat^{n+1} + 0.5*(phi^n - phi_hat^n). BFECC needs "three semi-Lagrangian advections"; MacCormack "requiring only twice the effort of the first order accurate scheme". 2nd order in space AND time confirmed; "the MacCormack method is identical to second order accurate Runge-Kutta for ordinary differential equations".
- Limiter confirmed: clamp the final result to min/max of the corner values "used in computing the first advection step"; the alternative "reverts to a first order accurate method when the higher order accurate method would overshoot"; the authors "prefer the reversion approach". Instability without a limiter confirmed: "this error correction step can lead to new extrema and possible instability."
- Reversion near solids confirmed: "automatically revert to the first order accurate scheme when the upwind and downwind building blocks pull data from non-commensurate regions (e.g. if one gets information from the fluid and the other from a solid wall boundary)". (The "within CFL*h of solids" radius is our operationalization of that rule.)
- Zalesak disc: tested at CFL = 0.75 and 1.75; measured MacCormack orders 1.59-2.08, settling at ~2.0-2.1 on the finest grids (paper Figs. 4-6 tables) — "~2.0" is fair, though intermediate grids dip to ~1.6.

### 3. Pressure projection (MAC grid) — CONFIRMED

- Bridson notes Eq. (4.1): "u^{n+1} = u - dt*(1/rho)*grad p"; substituting into the MAC divergence yields laplacian(p) = (rho/dt)*div(u*). Same convention (constant rho = 1025 kg/m3) is dimensionally consistent with p in Pa.
- 7-point stencil, diagonal = number of non-solid neighbors: Bridson §4.3: "the coefficient in front of p_ij is equal to the number of non-solid grid cell neighbours (this is the same in three dimensions)". Independent confirmation in McAdams et al. 2010: their discretization "is derived from the standard 7-point Poisson finite difference stencil using the zero-Neumann boundary condition (p_i'j'k' - p_ijk)/h = 0 to eliminate pressure values at Neumann cells" with "coefficients -6 and 1" in the interior — dropping the solid neighbor IS the dp/dn = 0 BC.
- All-Neumann singularity/compatibility: Bridson §4.3 discusses the "compatibility condition" (net flux through the boundary must vanish) for the singular semi-definite system; rescaling outlet flux to match inlet before the solve is the standard global-mass-balance fix for incompressible outflow (see e.g. FLOW-3D CFD-101 outflow notes; Ferziger & Peric). Pinning p = 0 at the outlet (Dirichlet) removes the null space — standard.

### 4. Jacobi Poisson iteration (3D) — CONFIRMED (two clarifications)

- Harris GPU Gems ch. 38 (fetched from developer.nvidia.com) verified: x^{k+1} = (xL + xR + xB + xT + alpha*b)*rBeta with alpha = -(dx)^2, rBeta = 1/4 for the pressure solve, and exact quotes "we typically use 40 to 80 Jacobi iterations" and "It is not a good idea to go below 20 iterations, because the error is noticeable."
- Clarification (a): Harris ch. 38 is explicitly 2D; the 6-neighbor /6 form in our note is the correct 3D generalization (rBeta = 1/6), cf. GPU Gems 3 ch. 30 for the 3D version.
- Clarification (b): Harris's b is just div(w) (rho/dt absorbed into p); with b = (rho/dt)*div(u*) as written here, p comes out in Pa — consistent with Sec. 3.
- Smooth-mode factor: undamped Jacobi eigenvalue for 1D/tensor Poisson is cos(pi*h) ~ 1 - pi^2/(2N^2) (Strang, MIT 18.086 multigrid notes §6.3, math.mit.edu/classes/18.086/2006/am63.pdf; Demmel CS267 Lecture 25). Refined count: 1e-4 reduction at N = 400 needs ln(1e4)/(pi^2/(2*400^2)) ~ 3e5 sweeps — same order as the "~1e5" claim but ~3x larger. Conclusion unchanged: debug/bring-up solver only.

### 5. MGPCG convergence — CORRECTED

Verified against the full paper text (math.ucdavis.edu/~jteran/papers/MST10.pdf):
- CORRECTION (authors): the paper is by McAdams, **Sifakis**, Teran (SCA 2010, pp. 65-73), not "Sitthi-amorn" (dblp: conf/sca/McAdamsST10). Fixed throughout this note.
- CORRECTION (smoother spec): the paper's smoothing stage is ONE damped Jacobi sweep (omega = 2/3: "Our smoother of choice is the damped Jacobi method with parameter omega = 2/3") over the whole domain per V-cycle stroke — not 2-3 pre/post sweeps. The "30-40 extra Gauss-Seidel sweeps" ("at least 30-40 iterations of boundary smoothing (15-20 iterations before the interior sweep, and 15-20 after) on a boundary band at least 3 cells wide") are what the V-cycle needed to be stable **as a standalone solver**. When the V-cycle is used as the CG preconditioner (the recommended MGPCG), they use only **2 Gauss-Seidel boundary sweeps at the finest level** ("using 2 Gauss-Seidel iterations on the boundary band at the finest level struck the best balance"), doubling per coarser level (2^{l+1} sweeps at level l), on a band 1-3 cells wide. Implement THAT, not 30-40 sweeps.
- Coarsening: standard 8-to-1 (factor 2 per axis) — confirmed.
- Convergence claims confirmed verbatim: "typically reduces the residual by one order of magnitude every 2 iterations"; "our convergence rates are independent of grid size". Paper Table 1 (smoke-past-sphere): 9-13 iterations to ||r|| = 1e-4 across 64^3-512^3 — our 8-14 iteration budget at 4-30M cells is consistent. "residual reduction by a factor of 10^4 at 128^3 resolution in less than 0.75 seconds" (16-core SMP, 2010) and "768^2 x 1152 voxels with a memory footprint less than 16GB" confirmed.
- Warm start: NOT from the paper — the V-cycle preconditioner mandates a zero initial guess internally; warm-starting the outer CG with the previous pressure is standard practice and remains recommended, but cite it as general practice, not McAdams et al.

### 6. Smagorinsky eddy viscosity — CONFIRMED

- Form nu_t = (Cs*Delta)^2*|S|, |S| = sqrt(2*S_ij*S_ij) is the standard Smagorinsky-Lilly closure (any LES reference, e.g. Wikipedia LES / COSMO-LES documentation). Delta = (dx*dy*dz)^{1/3} is Deardorff's proposal; = h for cubic voxels.
- Constants: Lilly's inertial-subrange analysis gives Cs ~ 0.17 (0.17-0.21 depending on Kolmogorov constant); with mean shear this over-damps, and Deardorff's channel-flow LES used Cs = 0.094 ~ 0.1 (sources: COSMO-LES Smagorinsky-Lilly note, cosmo-model.org; Springer TCFD "Smagorinsky constant in LES..."; ScienceDirect Smagorinsky-model overview). Practical range 0.05-0.2 and our recommendation Cs = 0.10-0.12 are consistent.
- Arithmetic verified: (0.12*0.05)^2*50 = (0.006)^2*50 = 1.8e-3 m2/s; explicit 3D diffusion limit dt <= h^2/(6*nu_t) = 0.0025/(6*1.8e-3) = 0.23 s >> advective dt. Explicit treatment is safe as claimed.

### 7. Orlanski convective outflow — CORRECTED (attribution detail)

- Reference confirmed: Orlanski, "A Simple Boundary Condition for Unbounded Hyperbolic Flows", J. Comput. Phys. 21:251-269, 1976 (ScienceDirect 0021-9991(76)90023-1). The radiation form d(phi)/dt + C*d(phi)/dn = 0 (Sommerfeld condition) confirmed.
- CORRECTION: in Orlanski's actual scheme C is NOT a fixed bulk speed — it is a local phase speed computed each step from the interior solution, C = -(dphi/dt)/(dphi/dx) evaluated one point inside at the previous time levels (leapfrog), then clamped to 0 <= C <= h/dt. Using a constant U_c = bulk outflow velocity with first-order upwind time stepping (u_b^{n+1} = u_b^n - (U_c*dt/h)*(u_b^n - u_{b-1}^n)) is a widely used SIMPLIFICATION (the "convective outflow" BC of Ferziger & Peric; also standard in ocean/LES codes, cf. WikiROMS radiation BCs). Keep the simple version (it is robust and adequate with the downstream flux rescale + p = 0), but require U_c*dt/h <= 1 for stability of the explicit upwind update, and cite it as "convective/Orlanski-type", not Orlanski's exact scheme.

### 8. Vorticity confinement — CONFIRMED

Verified against Fedkiw, Stam & Jensen, "Visual Simulation of Smoke" (web.stanford.edu/class/cs237d/smoke.pdf): omega = curl(u) is Eq. (9), N = grad|omega|/|grad|omega|| Eq. (10), and f_conf = epsilon*h*(N x omega) is Eq. (11), with "epsilon > 0 is used to control the amount of small scale detail added back into the flow field" and the h factor so that "as the mesh is refined the physically correct solution is still obtained". Confirmed the paper gives NO numeric epsilon (0.1-0.5 is folklore graphics practice). Keeping epsilon = 0 for quantitative sediment runs is our (sound) engineering decision, not a paper claim.

### Verification sources (independently fetched)

- Stam, Stable Fluids PDF: https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/ns.pdf
- Selle et al., MacCormack PDF: https://faculty.cc.gatech.edu/~jarek/papers/maccormack.pdf
- Fedkiw, Stam, Jensen, Smoke PDF: https://web.stanford.edu/class/cs237d/smoke.pdf
- Harris, GPU Gems ch. 38: https://developer.nvidia.com/gpugems/gpugems/part-vi-beyond-triangles/chapter-38-fast-fluid-dynamics-simulation-gpu
- McAdams, Sifakis, Teran PDF: https://www.math.ucdavis.edu/~jteran/papers/MST10.pdf ; dblp entry: https://dblp.uni-trier.de/rec/conf/sca/McAdamsST10.html
- Bridson, fluids notes PDF: https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf
- Orlanski 1976 (ScienceDirect record): https://www.sciencedirect.com/science/article/abs/pii/0021999176900231 ; WikiROMS radiation BCs: https://www.myroms.org/wiki/Boundary_Conditions
- Smagorinsky constants: https://www.cosmo-model.org/content/model/documentation/newsLetters/newsLetter12/2-langhans.pdf ; https://link.springer.com/article/10.1007/s00162-007-0064-z ; https://en.wikipedia.org/wiki/Large_eddy_simulation
- Jacobi/multigrid convergence factors: https://math.mit.edu/classes/18.086/2006/am63.pdf ; http://people.eecs.berkeley.edu/~demmel/cs267/lecture25/lecture25.html
