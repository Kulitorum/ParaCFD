# Alternative Solvers: Lattice Boltzmann & Others for Scour CFD

## Overview

Our plan-A fluid core is a Jos-Stam semi-Lagrangian "stable fluids" solver + Smagorinsky LES. Its known weaknesses are (1) strong numerical diffusion of the first-order semi-Lagrangian advection, and (2) near-bed shear computed by finite-differencing a velocity field on a staircase voxel bed — exactly the quantity (Shields stress) our whole sediment model keys on. The lattice Boltzmann method (LBM) is the strongest alternative: it is explicit, local, embarrassingly parallel (bandwidth-bound, ideal for an RTX 4090), and delivers the full viscous stress tensor *locally per cell* from the non-equilibrium populations — no gradient differencing at walls. This note gives the implementation-critical LBM equations, GPU performance/memory numbers, license status of existing codes, and a concrete recommendation.

Reference frame: domain 10 x 10 x 5 m, dx = 2.5 cm gives 400 x 400 x 200 = 32e6 cells; u = 0.5–2.5 m/s; nu = 1.05e-6 m2/s; cell Reynolds number Re_cell = u*dx/nu ≈ 12,000–60,000, i.e. wall-modeled LES territory for ANY solver we pick.

## Key equations & results

### 1. LBM–BGK update and equilibrium (Krüger et al., *The Lattice Boltzmann Method*, Springer 2017)

f_i(x + c_i*dt, t + dt) = f_i(x,t) - (dt/tau) * [f_i(x,t) - f_i_eq(x,t)],  i = 0..Q-1

f_i_eq = w_i * rho * [ 1 + (c_i.u)/cs^2 + (c_i.u)^2/(2*cs^4) - (u.u)/(2*cs^2) ],  cs^2 = (1/3)*(dx/dt)^2 (lattice sound speed; cs = 1/sqrt(3) in lattice units).

D3Q19 weights: w0 = 1/3; w(6 axis dirs) = 1/18; w(12 diagonal dirs) = 1/36. D3Q27 weights: 8/27, 2/27 (6 axis), 1/54 (12 edge), 1/216 (8 corner). Macroscopic: rho = sum_i f_i, rho*u = sum_i c_i*f_i.

D3Q19 vs D3Q27: D3Q19 has anisotropy defects in third-order moments that produce spurious secondary flows at high Re (visible in square-duct/jet flows); D3Q27 restores rotational invariance at 27/19 = 1.42x the memory and ~1.4x the bandwidth cost. The cumulant operator (below) is formulated on D3Q27. For a first LES build, D3Q19 + TRT is acceptable; for production-quality high-Re accuracy, D3Q27 cumulant is the reference.

### 2. Viscosity–relaxation relation, stability and Mach limits (Krüger et al. 2017)

nu = cs^2 * (tau/dt - 1/2) * dx^2/dt   →  in lattice units (dx = dt = 1):  nu_lat = (1/3)*(tau - 0.5).

Constraints: tau > 0.5 strictly; plain BGK becomes unstable for tau below ≈ 0.505–0.51 in under-resolved turbulence. Velocity must satisfy u_lat = u_phys*dt/dx ≤ 0.1–0.15 to keep compressibility error O(Ma^2) ≤ ~1% (hard stability limit ≈ 0.3–0.4 for BGK; cumulant survives to ≈ 0.4+).

Unit conversion for us: choose u_lat = 0.1 at u_phys = 2.5 m/s, dx = 2.5 cm → dt = u_lat*dx/u_phys = 1.0e-3 s. Then nu_lat = nu_phys*dt/dx^2 = 1.05e-6*1e-3/6.25e-4 = 1.68e-6 → tau_0 = 0.500005. That is hopelessly below the BGK stability floor — molecular viscosity is unresolvable at our grid, and the run is stabilized entirely by the LES eddy viscosity plus a robust collision operator (TRT/MRT/cumulant). This is standard practice (FluidX3D, waLBerla, VirtualFluids all do exactly this), but it means: **never run plain BGK without LES at our Re**.

### 3. Collision operators: BGK vs TRT vs MRT vs cumulant

- **BGK (SRT)**: 1 relaxation rate; cheapest; unstable at tau→0.5; bounce-back wall location drifts with viscosity.
- **TRT** (Ginzburg/d'Humières): two rates tau+, tau-; with the "magic parameter" **Lambda = (tau+ - 1/2)(tau- - 1/2) = 3/16** the half-way bounce-back wall sits *exactly* mid-link independent of viscosity (Krüger et al. 2017, ch. 10). Nearly BGK cost, much better wall behavior — best cost/benefit for our voxel bed.
- **MRT** (d'Humières et al., Phil. Trans. R. Soc. A 360:437, 2002): relaxes 19/27 moments individually; wider stability, but tunable free rates, mild Galilean-invariance violation and hyper-viscosity artifacts.
- **Cumulant** (Geier et al., Comput. Math. Appl. 70:507–547, 2015, D3Q27): collides statistically independent cumulants; Galilean invariant, no hyper-viscosity, stable at extremely high Re even with tau_0 - 0.5 ~ 1e-6; demonstrated well-resolved flows up to Re ~ 1e6. This is the operator used by VirtualFluids (TU Braunschweig) for environmental/civil flows. Cost ≈ 1.5–2x BGK arithmetic, but the method stays bandwidth-bound on GPU, so wall-clock penalty is mostly the 27/19 memory-traffic factor.

### 4. Smagorinsky LES–LBM, local closed form (Hou et al., arXiv comp-gas/9401004, 1994; as implemented in FluidX3D)

Compute the non-equilibrium momentum flux locally: Pi_ab_neq = sum_i c_ia*c_ib*(f_i - f_i_eq), Q = Pi_neq : Pi_neq (double contraction). Then in lattice units:

tau_eff = 1/2 * ( tau_0 + sqrt( tau_0^2 + 18*sqrt(2)*Cs^2 * sqrt(Q) / rho ) ),  Cs = 0.10–0.17 (Smagorinsky constant; FluidX3D hard-codes the prefactor 16*sqrt(2)/(3*pi^2) = 0.7639, equivalent to Cs ≈ 0.173).

This is algebraic and per-cell — no velocity gradients needed, unlike grid LES. Use Cs ≈ 0.1 in sheared near-bed regions (classical channel-flow value; 0.17 is the isotropic-turbulence value and over-damps shear flows). Van Driest damping is impossible without wall distance; at 2.5 cm cells we are wall-modeled anyway (apply the log-law bed shear model from the sediment notes at the first cell, same as for stable fluids).

### 5. Local strain rate and wall shear stress (Krüger et al. 2017; Chai & Zhao local-scheme papers)

S_ab = - Pi_ab_neq / (2 * rho * cs^2 * tau_eff)  (lattice units: S_ab = -3*Pi_ab_neq/(2*rho*tau_eff)),  second-order accurate, fully local.

Viscous stress: sigma_ab = 2*rho*nu_eff*S_ab. Bed shear stress on a voxel face with unit normal n: traction t_a = sigma_ab*n_b; wall shear tau_w = |t - (t.n)n|. Convert to physical units via tau_phys = tau_lat * rho_phys/rho_lat * (dx/dt)^2. **This locality is LBM's single biggest advantage for scour**: Shields stress theta = tau_w/((rho_s - rho)*g*d50) comes straight from per-cell data with no one-sided finite differences on the staircase bed.

### 6. Boundary conditions on voxel geometry

- **Half-way bounce-back**: f_ibar(x_f, t+dt) = f_i_postcollision(x_f, t); wall is implicitly located halfway along the cut link. Second-order accurate for grid-aligned walls, degrades to ~first order on staircase (inclined/curved) surfaces; staircase corners shed spurious small vortices at high Re. With TRT Lambda = 3/16 the mid-link location is exact and viscosity-independent — always use this, since our bed and STL obstacles are voxelized.
- **Bouzidi interpolated (curved) bounce-back** (Bouzidi, Firdaouss, Lallemand, Phys. Fluids 13:3452, 2001): with q = distance(fluid node → wall)/distance(fluid node → solid node) along link i, post-collision values f*:
  - q < 1/2:  f_ibar(x_f, t+dt) = 2q * f_i*(x_f) + (1 - 2q) * f_i*(x_f - c_i*dt)
  - q ≥ 1/2:  f_ibar(x_f, t+dt) = (1/(2q)) * f_i*(x_f) + ((2q-1)/(2q)) * f_ibar*(x_f)
  Restores second order on curved geometry; requires storing q per cut link (from the STL). Worth it for the printed-shape surface; overkill for the evolving sand bed (which is genuinely voxel-resolved and re-voxelized every Exner step).
- Inflow/outflow: equilibrium or non-equilibrium-extrapolation velocity/pressure BCs; sponge/ramped-viscosity zone of ~10–20 cells at the outlet to kill reflections.

### 7. GPU performance and memory (FluidX3D benchmarks, ProjectPhysX GitHub; Wittmann et al. arXiv:1112.0850)

LBM is memory-bandwidth-bound: MLUPS ≈ efficiency * BW / bytes_per_cell_update, efficiency 80–90% achievable. RTX 4090 (1008 GB/s), D3Q19 with in-place ("Esoteric-Pull"/AA-pattern) streaming:

- Memory: **93 B/cell** (FP32 storage) or **55 B/cell** (FP32 compute / FP16 storage), vs 344 B/cell naive two-lattice + macros.
- Measured FluidX3D throughput on RTX 4090, D3Q19: **5,624 MLUPS (FP32/FP32), 11,091 MLUPS (FP32/FP16S), 11,496 MLUPS (FP32/FP16C)**. A competent custom CUDA kernel should reach 4,000–5,500 MLUPS FP32; D3Q27 cumulant ≈ 0.6–0.7x of that (~3,000–3,800 MLUPS).
- Our 32e6-cell domain: 3.0 GB (D3Q19 FP32) / 1.8 GB (FP16S) / ~4.1 GB (D3Q27 FP32) — trivially fits in 24 GB; even dx = 1.25 cm (256e6 cells) fits at 14 GB in FP16S.
- Wall-clock: dt = 1 ms → 300 s of flow = 3.0e5 steps ≈ 9.6e12 cell-updates ≈ **~30 min at 5 GLUPS (FP32), ~15 min FP16S** per shape evaluation at dx = 2.5 cm. Comparable to (and possibly faster than) our stable-fluids solver, which takes ~20x fewer steps (dt ≈ 20 ms at CFL 2) but pays a multigrid pressure solve every step.

### 8. Published LBM sediment/scour work (feasibility evidence)

- *Sediment transport in turbulent flows with the lattice Boltzmann method*, Computers & Fluids (2018): D3Q27 entropic multi-relaxation (KBC) LBM + grid refinement + Grad off-lattice BCs; suspended + bedload + bed evolution around arbitrarily shaped objects in unidirectional turbulent flow — essentially our exact use case, demonstrating LBM shear fields drive credible morphodynamics.
- *DNS of sediment transport in open channel flow with LBM*, MDPI Fluids 6(6):217 (2021): particle-resolved LBM DNS (research-grade, not our scale).
- DBM–LBM coupled local-scour model around double triangular prisms, JMSE 14(10):941 (2D hydro-morphodynamic).
- Multiple LBM–DEM pile-scour studies (grain-resolved; 1000x too expensive for 10 m domains but validates near-bed LBM shear).

## Practical guidance for our simulator

**SPH**: Lagrangian particles; needs ~10x the DOF of a grid for the same near-bed resolution, noisy pressure (WCSPH) and kernel-deficient boundaries make wall shear the *worst*-computed quantity; turbulence modeling immature. Only compelling for violent free-surface/FSI. **Reject.**

**FLIP/PIC**: hybrid particle–grid built for free-surface liquid animation; velocity gradients inherit particle-transfer noise, so bed shear is worse than the underlying grid solver alone; our domain is fully submerged (no free surface), so FLIP adds cost with zero benefit. **Reject.**

**FVM references (OpenFOAM/sedFoam, REEF3D, Delft3D)**: use as validation cross-checks, not as our GPU core (CPU-bound, hours-per-case).

**Licenses — hard constraints for COBOD (commercial):**
- **FluidX3D: NOT usable.** Custom "free for non-commercial use" license; any COBOD R&D use is commercial use. Do not copy code. Its *published techniques* (Esoteric-Pull streaming, FP16 storage, benchmark methodology) are in open-access papers and freely reimplementable. **Highest license risk — flag to anyone tempted to fork it.**
- **Palabos: AGPLv3** — copyleft triggers even on network-served use; avoid linking.
- **OpenLB: GPLv2; waLBerla: GPLv3** — internal use is legal (copyleft triggers on distribution), but linking either into a simulator we ever ship to partners forces open-sourcing. Treat as reference/validation tools only.
- **STLBM (Latt et al., UNIGE, gitlab.com/UniGeHPFS/stlbm): MIT** — clean-room quality reference for collision models (incl. cumulant) we may freely port to CUDA. **Safest starting point.**
- Writing our own CUDA D3Q19 kernel set (stream-collide + bounce-back + LES) is ~1–2k lines; this is the zero-risk path.

**Recommendation:**
1. **Stable fluids stays plan-A for the shape-ranking milestone**, with two mandatory mitigations: (a) MacCormack/BFECC advection (second order, cuts the numerical diffusion that would otherwise swamp the Smagorinsky term), and (b) bed shear from a log-law wall model using the first-cell velocity (z0 = ks/30, ks = 2.5*d50) — never from raw one-sided gradients on the voxel staircase. For *relative ranking* of shapes under identical numerics, consistent bias largely cancels; validate against the Sumer & Fredsøe live-bed cylinder benchmark (equilibrium S/D ≈ 1.3).
2. **LBM is plan-B, promoted to plan-A for the fluid core if** stable-fluids shear maps prove too diffusive/noisy to discriminate candidate shapes (test: does it reproduce the horseshoe-vortex amplification factor ~3–4x ambient tau_b at a cylinder base?). Architecture requirement now: hide the fluid core behind an interface `FluidCore::step() -> {u(x), tau_b(bed cells)}` so the sediment/Exner module is solver-agnostic.
3. If/when building LBM: **D3Q19 + TRT (Lambda = 3/16) + Smagorinsky (Cs = 0.1), FP32 compute / FP16 storage, in-place streaming, half-way bounce-back on bed, Bouzidi on the printed shape**; upgrade path to D3Q27 cumulant (port from MIT STLBM) if stability or secondary-flow accuracy demands it. Expected cost: ~15–30 min per 300 s flow evaluation at dx = 2.5 cm on the RTX 4090 — same order as stable fluids, with strictly better near-bed shear locality.

## Sources

- FluidX3D — ProjectPhysX (benchmarks, memory model, license): https://github.com/ProjectPhysX/FluidX3D
- Geier, Schönherr, Pasquali, Krafczyk, "The cumulant lattice Boltzmann equation in three dimensions: Theory and validation", Comput. Math. Appl. 70 (2015): https://sciencedirect.com/science/article/pii/S0898122115002126
- Hou, Sterling, Chen, Doolen, "A Lattice Boltzmann Subgrid Model for High Reynolds Number Flows" (1994): https://arxiv.org/pdf/comp-gas/9401004
- Bouzidi, Firdaouss, Lallemand — interpolated bounce-back, via review "LBM for isothermal micro-gaseous flow": https://arxiv.org/pdf/1508.01562
- "Measurements of wall shear stress with the lattice Boltzmann method and staircase approximation of boundaries", Computers & Fluids (2010): https://www.sciencedirect.com/science/article/abs/pii/S0045793010001283
- "Wall orientation and shear stress in the lattice Boltzmann model": https://arxiv.org/pdf/1203.3078
- "Sediment transport in turbulent flows with the lattice Boltzmann method", Computers & Fluids (2018): https://www.sciencedirect.com/science/article/abs/pii/S0045793018302111
- "DNS of Sediment Transport in Turbulent Open Channel Flow Using LBM", Fluids 6(6):217 (2021): https://www.mdpi.com/2311-5521/6/6/217
- "Local scour around double triangular prisms, DBM–LBM coupled model", JMSE 14(10):941: https://doi.org/10.3390/jmse14100941
- Latt et al., "Cross-platform programming model for many-core lattice Boltzmann simulations" (STLBM, MIT license): https://arxiv.org/abs/2010.11751
- Palabos (AGPLv3): https://palabos.unige.ch/ ; VirtualFluids (cumulant LBM): https://www.sciencedirect.com/science/article/pii/S0010465525003121
- Wittmann et al., "Performance engineering for LBM on GPGPUs": https://arxiv.org/pdf/1112.0850
- Krüger, Kusumaatmaja, Kuzmin, Shardt, Silva, Viggen, *The Lattice Boltzmann Method: Principles and Practice*, Springer 2017 (textbook; TRT magic parameter, D3Q19/D3Q27 moments, local stress).

## Verification

Independent adversarial fact-check (2026-07-07). Every equation was re-derived or checked against primary sources fetched independently of this note's citations. Verdicts: CONFIRMED / CORRECTED / UNVERIFIABLE.

### 1. LBM–BGK update + equilibrium, D3Q19/D3Q27 weights — CONFIRMED

- Update rule f_i(x + c_i dt, t+dt) = f_i − (dt/τ)(f_i − f_i^eq) and the second-order equilibrium f_i^eq = w_i ρ [1 + (c_i·u)/cs² + (c_i·u)²/(2cs⁴) − u·u/(2cs²)] match the standard form (e.g. the equivalent 3/c², 9/(2c⁴), 3/(2c²) coefficients with cs² = c²/3 = (1/3)(dx/dt)²) — [Wikipedia: Lattice Boltzmann methods](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods); Hou et al. 1994 eq. (13) (2D analogue with identical coefficients).
- D3Q19 weights 1/3 (rest), 1/18 (6 axis), 1/36 (12 diagonal): confirmed ([Wikipedia](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods)). D3Q27 weights 8/27, 2/27 (6 axis), 1/54 (12 edge), 1/216 (8 corner): confirmed ([imechanica D3Q27 note](https://imechanica.org/sites/default/files/D3Q27.pdf); [ResearchGate D3Q27 figure](https://www.researchgate.net/figure/The-lattice-cell-of-D3Q27-Lattice-Boltzmann-method-The-weight-factors-of-lattice-wk-are_fig1_329953746)). Weights sum to 1 in both sets (sanity check: 1/3 + 6/18 + 12/36 = 1; 8/27 + 12/27 + 12/54 + 8/216 = 1).

### 2. Viscosity–relaxation relation — CORRECTED (dimensional form); all derived numbers CONFIRMED

- The compact physical-units formula as originally extracted, nu = cs²·(τ/dt − 1/2)·dx²/dt, is **dimensionally inconsistent** if cs² = (1/3)(dx/dt)² as defined in §1 (it double-counts (dx/dt)², yielding m⁴/s³). Correct physical-units statement: **ν = cs²(τ − dt/2)** with cs² = (1/3)(dx/dt)², equivalently ν = (1/3)(dx²/dt)(τ/dt − 1/2). Source: standard Chapman–Enskog result, ν = cs²(τ − Δt/2) ([arXiv:1211.0205](https://arxiv.org/pdf/1211.0205); [Wikipedia](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) ν = (τ − 0.5)cs²δt); lattice-units form ν_lat = (1/3)(τ − 0.5) independently confirmed by Hou et al. 1994 eq. (16), ν = (2τ−1)/6.
- τ > 0.5 strictly; BGK instability approaching τ → 0.5 in under-resolved turbulence: confirmed (Hou et al. 1994 report blow-up at Re ≳ 10⁴ on 256² without LES).
- COBOD unit conversion re-computed and confirmed: dt = u_lat·dx/u_phys = 0.1×0.025/2.5 = **1.0e-3 s**; ν_lat = ν_phys·dt/dx² = 1.05e-6×1e-3/6.25e-4 = **1.68e-6**; τ₀ = 0.5 + 3ν_lat = **0.5000050**. u_lat ≤ 0.1–0.15 (Ma² compressibility error) is standard practice guidance (Krüger et al. ch. 7; FluidX3D docs) — reasonable, not a hard constant.

### 3. Smagorinsky LES–LBM effective relaxation time — CONFIRMED (formula); CORRECTED (prefactor digits)

- Re-derived from the primary source [Hou et al. 1994, arXiv:comp-gas/9401004](https://arxiv.org/abs/comp-gas/9401004): eq. (14) Π⁽¹⁾_ab = −(2nτ/3)S_ab, eq. (25) ν_total = ν₀ + CΔ²|S| with |S| = √(2S_ab S_ab), τ_total = 3ν_total + 1/2, Q = Π_ab Π_ab. Solving the resulting quadratic with Δ = 1, C = Cs² gives exactly **τ_eff = ½(τ₀ + √(τ₀² + 18√2·Cs²·√Q/ρ))** — the stated form is correct. (Beware: the arXiv preprint's own printed eq. (26) and its solved form contain typos — a missing √2 from |S| = √(2S:S) and ν₀ written where τ₀ belongs; the form above is the standard corrected one, identical to what FluidX3D documents.)
- FluidX3D documentation ([DOCUMENTATION.md](https://github.com/ProjectPhysX/FluidX3D/blob/master/DOCUMENTATION.md)) states τ = ½(τ₀ + √(τ₀² + (16√2)/(3π²)·√Q/ρ)) — same functional form. Matching prefactors: 18√2·Cs² = 16√2/(3π²) ⟺ Cs² = 8/(27π²) ⟺ **Cs = 0.1733**. Numeric correction: **16√2/(3π²) = 0.76421**, not 0.7639 as stated (4th-digit error). Cs = 0.10–0.17 usage range is standard (Hou et al. tested C = Cs² = 0.0025–0.04, i.e. Cs = 0.05–0.2, and needed Cs² ≥ 0.073–0.084 at Re = 10⁵–10⁶ in the cavity — constant is geometry-dependent).

### 4. Local strain rate and bed shear stress from Π_neq — CONFIRMED

- S_ab = −Π_ab^neq/(2ρcs²τ_eff): confirmed directly by Hou et al. 1994 eq. (14) (Π⁽¹⁾ = −2ρcs²τ·S with cs² = 1/3, i.e. −(2nτ/3)S); lattice form −3Π_neq/(2ρτ_eff) follows with cs² = 1/3. The locality claim (no finite differences) is the explicit point of Hou et al. §3(c) and their concluding remarks.
- σ_ab = 2ρν_eff S_ab is the deviatoric viscous stress (equivalently σ_ab = −(1 − 1/(2τ))Π_ab^neq, consistent with ν = cs²(τ − 1/2)); traction t_a = σ_ab n_b and tangential projection τ_w = |t − (t·n)n| are standard continuum mechanics. Stress conversion τ_phys = τ_lat·(ρ_phys/ρ_lat)·(dx/dt)² is dimensionally exact (stress ~ ρu²). Shields θ = τ_w/((ρ_s−ρ)g·d50) is the standard Shields parameter definition.
- Caveat kept from the check: use the **total** τ_eff (and ν_eff) consistently in both S_ab and σ_ab, as written.

### 5. TRT magic parameter Λ = 3/16 — CONFIRMED (with scope note)

- Λ = (τ⁺ − 1/2)(τ⁻ − 1/2) = 3/16 makes the half-way bounce-back wall location viscosity-independent, and for Poiseuille-type (parabolic) channel flow places the no-slip wall **exactly** midway between solid and fluid nodes (Ginzburg et al., "magic" collision numbers; independently confirmed via [Khirevich et al., J. Comput. Phys. 2015](https://www.sciencedirect.com/science/article/abs/pii/S0021999114007207) and [arXiv:2009.04604](https://arxiv.org/pdf/2009.04604)). Scope note: exactness is proven for parabolic profiles/straight walls; in general geometry Λ = 3/16 gives viscosity-independent (not pointwise-exact) wall placement — still the right default for our voxel bed.

### 6. Bouzidi interpolated bounce-back — CONFIRMED

- Independent statement of the linear BFL scheme ([arXiv:2203.01316](https://arxiv.org/pdf/2203.01316), citing Bouzidi, Firdaouss, Lallemand, Phys. Fluids **13**(11):3452, 2001 — journal reference confirmed via [AIP](https://pubs.aip.org/aip/pof/article-pdf/13/11/3452/12713899/3452_1_online.pdf)): q < 1/2: f_ī(x_f, t+dt) = 2q·f_i*(x_f) + (1−2q)·f_i*(x_ff); q ≥ 1/2: f_ī(x_f, t+dt) = (1/(2q))·f_i*(x_f) + ((2q−1)/(2q))·f_ī*(x_f), with q = |x_f − x_w|/|x_f − x_b| ∈ [0,1] and f* post-collision. Identical to the note (x_ff = x_f − c_i dt is the fluid neighbor away from the wall). Second-order accuracy on curved boundaries confirmed by the same sources; plain bounce-back is first-order at non-grid-aligned boundaries (also stated in Hou et al. 1994, citing Hou et al. [28]).

### 7. GPU performance and memory model — CORRECTED (bytes-per-update vs storage conflation); all measured numbers CONFIRMED

- Confirmed from [FluidX3D README](https://github.com/ProjectPhysX/FluidX3D): RTX 4090 (1008 GB/s theoretical BW), D3Q19: **5,624 MLUPS (FP32/FP32), 11,091 MLUPS (FP32/FP16S), 11,496 MLUPS (FP32/FP16C)**; storage **93 B/cell (FP32)**, **55 B/cell (FP16 compressed)**, vs ~344 B/cell traditional FP64 two-lattice.
- **Correction**: in MLUPS ≈ 0.8–0.9 × BW / bytes_per_cell_update, the denominator is the memory **traffic per cell per step**, which FluidX3D states is **153 B (FP32/FP32)** and **77 B (FP32/FP16)** with in-place (Esoteric-Pull) streaming — *not* the 93/55 B storage footprint. Check: 1008/153 = 6,588 MLUPS ideal → 5,624 measured = 85% efficiency; 1008/0.077 = 13,090 → 11,091 = 85%. Using 93 B would predict ~9,200 MLUPS at 85%, contradicting the measurement.
- Wall-clock arithmetic re-computed and confirmed: 32e6 cells × 3.0e5 steps = 9.6e12 updates → 9.6e12/5.624e9 = 1,707 s ≈ **28 min FP32**; 9.6e12/1.1091e10 = 866 s ≈ **14.4 min FP16S**. Memory: 32e6 × 93 B = 2.98 GB; × 55 B = 1.76 GB; 256e6 × 55 B = 14.1 GB — all fit a 24 GB RTX 4090 as stated.

### 8. Lattice unit scaling recipe — CONFIRMED (with one units clarification)

- dt = u_lat·dx/u_phys_max, ν_lat = ν_phys·dt/dx², τ₀ = 0.5 + 3ν_lat (τ_total = 3ν_total + 1/2 confirmed verbatim in Hou et al. 1994), and acceleration scaling g_lat = g_phys·dt²/dx are all dimensionally exact; g_lat = 9.81×(1e-3)²/0.025 = **3.924e-4** for dx = 2.5 cm, dt = 1 ms — arithmetic confirmed. See also [dbarker.uk LBM parameterization](https://dbarker.uk/posts/lbm-parameterization/) for the same recipe.
- Clarification: the (dx/dt)²·ρ_phys/ρ_lat factor converts **stresses/pressures** (force per area, ~ρu²). Body-force *density* (N/m³) converts by (ρ_phys/ρ_lat)·(dx/dt²)·(1/1) = ρ·dx/dt², and accelerations by dx/dt² (inverse dt²/dx when going physical → lattice, as used for g_lat). The recipe as applied in this note (gravity handled as an acceleration) is correct.

**Net result**: no functional-form errors survive except the two flagged — (a) the compact viscosity formula in §2's extracted form was dimensionally garbled (the in-file text is saved by its lattice-units reading, but use ν = cs²(τ − dt/2) in the spec), and (b) in the GPU model use 153/77 B *traffic* per update, not 93/55 B *storage*, in the MLUPS estimate. One numeric constant fixed: FluidX3D LES prefactor 16√2/(3π²) = 0.76421 (Cs ≈ 0.1733), not 0.7639.
