# Validation Benchmarks & Experimental Datasets

## Overview

This note defines the validation ladder for the C++/CUDA morphodynamic simulator: (1) pure-flow analytic/benchmark tests (cavity, Poiseuille, log-law channel), (2) bluff-body tests (cylinder Cd/St, backward-facing step), (3) sediment-only tests (settling column, Rouse profile), (4) full coupled pile-scour benchmarks (Roulund et al. 2005; Sumer & Fredsoe; Melville & Coleman). Every test lists reference numbers and the accuracy that published RANS/LES scour models actually achieve, so acceptance tolerances are realistic, not aspirational. Key headline: a uniform Cartesian grid of D/10 cells per pile diameter reproduced the Roulund equilibrium scour depth to <4% (REEF3D::CFD), and state-of-the-art models still show up to 30% error on downstream scour and factor 2-4 errors on scour time scale — our tolerances should reflect that.

## Key equations & results

### 1. Lid-driven cavity (Ghia, Ghia & Shin 1982, J. Comp. Phys. 48:387-411)

Square cavity, lid velocity U=1, Re = U*L/nu. Reference u-velocity along the vertical centerline (x=0.5) and v-velocity along the horizontal centerline (y=0.5), 129x129 multigrid solution:

u(y) at x=0.5:

| y      | Re=100    | Re=1000   |
|--------|-----------|-----------|
| 1.0000 |  1.00000  |  1.00000  |
| 0.9766 |  0.84123  |  0.65928  |
| 0.9688 |  0.78871  |  0.57492  |
| 0.9609 |  0.73722  |  0.51117  |
| 0.9531 |  0.68717  |  0.46604  |
| 0.8516 |  0.23151  |  0.33304  |
| 0.7344 |  0.00332  |  0.18719  |
| 0.6172 | -0.13641  |  0.05702  |
| 0.5000 | -0.20581  | -0.06080  |
| 0.4531 | -0.21090  | -0.10648  |
| 0.2813 | -0.15662  | -0.27805  |
| 0.1719 | -0.10150  | -0.38289  |
| 0.1016 | -0.06434  | -0.29730  |
| 0.0703 | -0.04775  | -0.22220  |
| 0.0625 | -0.04192  | -0.20196  |
| 0.0547 | -0.03717  | -0.18109  |
| 0.0000 |  0.00000  |  0.00000  |

v(x) at y=0.5:

| x      | Re=100    | Re=1000   |
|--------|-----------|-----------|
| 1.0000 |  0.00000  |  0.00000  |
| 0.9688 | -0.05906  | -0.21388  |
| 0.9609 | -0.07391  | -0.27669  |
| 0.9531 | -0.08864  | -0.33714  |
| 0.9453 | -0.10313  | -0.39188  |
| 0.9063 | -0.16914  | -0.51550  |
| 0.8594 | -0.22445  | -0.42665  |
| 0.8047 | -0.24533  | -0.31966  |
| 0.5000 |  0.05454  |  0.02526  |
| 0.2344 |  0.17527  |  0.32235  |
| 0.2266 |  0.17507  |  0.33075  |
| 0.1563 |  0.16077  |  0.37095  |
| 0.0938 |  0.12317  |  0.32627  |
| 0.0781 |  0.10890  |  0.30353  |
| 0.0703 |  0.10091  |  0.29012  |
| 0.0625 |  0.09233  |  0.27485  |
| 0.0000 |  0.00000  |  0.00000  |

Key scalar checks: u_min on centerline = -0.21090 at y=0.4531 (Re=100) and -0.38289 at y=0.1719 (Re=1000). (Original table has known misprints in the Re=400/3200 columns; Re=100 and Re=1000 columns are the trusted ones.)

### 2. Laminar Poiseuille / open-channel (analytic)

Plane Poiseuille (gap H, body force g_x): u(y) = (g_x/(2*nu)) * y*(H-y); u_max = g_x*H^2/(8*nu) = 1.5*u_avg.
Laminar free-surface channel (depth h, free-slip top): u(z) = (g_x/nu)*(h*z - z^2/2); u_max = g_x*h^2/(2*nu) at surface; u_avg = g_x*h^2/(3*nu). Target: L2 error < 1-2% (grid convergence 2nd order; semi-Lagrangian advection contributes no error here since steady and unidirectional).

### 3. Turbulent open-channel log law (Nezu & Nakagawa 1993; Soulsby 1997)

Rough bed: u(z)/u_* = (1/kappa) * ln(z/z0), kappa = 0.40 (Soulsby) to 0.412 (Nezu & Nakagawa), z0 = k_s/30, k_s = 2.5*d50 (Nikuradse), valid for z/h < 0.2 and k_s+ = k_s*u_*/nu > 70 (rough regime). Smooth: u+ = (1/0.41)*ln(y+) + 5.0-5.29, y+ > 30.
Depth-averaged: V/u_* = (1/kappa)*(ln(h/z0) - 1). Worked check at our scale: h = 0.4 m, d50 = 0.26 mm -> z0 = 2.17e-5 m; u_* = 0.02 m/s gives V = 0.44 m/s (matches Roulund flume: V = 0.46 m/s). Acceptance: recovered kappa within 5%, u_* from log-fit within 10%.

### 4. Flow past a circular cylinder — Cd and St vs Re

Consensus experimental/numerical reference values (Sumer & Fredsoe 2006 "Hydrodynamics Around Cylindrical Structures"; Williamson 1996; Norberg 2003; Tritton 1959; 2D numerical benchmarks):

| Re    | Cd (mean)       | St            | Regime |
|-------|-----------------|---------------|--------|
| 40    | 1.50-1.60       | none (steady) | twin attached vortices, bubble length L/D ~ 2.2 |
| 100   | 1.33-1.40       | 0.164-0.168 (exp 0.164) | laminar 2D shedding |
| 200   | 1.30-1.40 (2D)  | 0.196-0.20 (2D); exp 0.18-0.19 | 3D mode-A begins at Re~190 |
| 1000  | 0.98-1.10       | ~0.21         | subcritical, turbulent wake |
| 1e4   | 1.10-1.20       | 0.19-0.20     | subcritical |
| 1e5   | ~1.20           | ~0.19         | subcritical (drag crisis at Re = 2e5-5e5, Cd -> 0.2-0.3) |

St = f*D/U (f = shedding frequency), Cd = 2*F_x/(rho*U^2*D*L). Full-scale monopile flows are post-critical (Re > 1e6, Cd ~ 0.6-0.8); a wall-modelled LES on a voxel grid will NOT capture the drag crisis — validate in the ranges Re=100-200 (laminar, tight tolerance: Cd within 5-10%, St within 3-5%) and Re=1e4 (LES, Cd within 15%, St within 10%). Set target Re via artificial nu so the cylinder spans >= 20-40 voxels.

### 5. Backward-facing step (Armaly et al. 1983, JFM 127:473-496; Gartling 1990, IJNMF 11:953-967)

Re = U_avg*D_h/nu, D_h = 2*(inlet channel height); step height S; expansion ratio ~2 (Armaly: 1.94). Primary reattachment length x1/S:
- Re=100: x1/S ~ 3 (Armaly exp).
- Re=200: x1/S = 5.0 (Armaly exp); 2D numerics reproduce 4.96.
- Re=800: Gartling 2D benchmark: x1/S = 12.20; upper-wall separation x4/S = 9.70; upper-wall reattachment x5/S = 20.96 (often quoted as 6.10/4.85/10.48 in units of full channel height H = 2S). Caution: Armaly's EXPERIMENT gives x1/S ~ 14 at Re=800 because the flow is 3D for Re > ~400; compare 2D runs only for Re <= 400, and to Gartling's 2D solution at Re=800.

### 6. Settling column (analytic; Soulsby 1997 "Dynamics of Marine Sands", eq. 102)

w_s = (nu/d) * [ sqrt(10.36^2 + 1.049*D*^3) - 10.36 ],  D* = [g*(s-1)/nu^2]^(1/3) * d, valid all d, natural (irregular) grains.
With rho_s=2650, rho=1025 kg/m3 (s=2.585), nu=1.05e-6 m2/s, g=9.81 m/s2 -> D* = 24160*d(m):

| d50     | D*   | w_s          | theta_cr (Soulsby-Whitehouse) |
|---------|------|--------------|-------------------------------|
| 0.1 mm  | 2.4  | 0.0073 m/s   | 0.080 |
| 0.26 mm | 6.3  | 0.036 m/s    | 0.042 |
| 1.0 mm  | 24.2 | 0.117 m/s    | 0.031 |
| 10 mm   | 242  | 0.40 m/s     | 0.056 |

theta_cr = 0.30/(1 + 1.2*D*) + 0.055*[1 - exp(-0.020*D*)] (Soulsby & Whitehouse 1997). Test: still water column, verify a concentration front advects down at exactly w_s (checks the settling term and mass conservation to <1%).

### 7. Equilibrium Rouse profile (Rouse 1937; Soulsby 1997)

C(z)/C_a = [ ((h-z)/z) * (a/(h-a)) ]^R,  R = w_s/(beta*kappa*u_*), kappa=0.40, beta~1 (beta = sediment/momentum diffusivity ratio, 1-2), valid a < z < h. Derives from eddy diffusivity eps_s = beta*kappa*u_**z*(1-z/h); the simulator's LES diffusivity + settling must reproduce it. Concrete test: h=0.4 m, u_*=0.02 m/s, d50=0.1 mm (w_s=0.0073 m/s) -> R=0.91; fix C=C_a at reference height a=0.05*h and compare steady profile. Suspension is active when u_*/w_s > ~1 (0.1 mm sand: u_* > 0.008 m/s — always, in our 0.5-2.5 m/s currents; 10 mm gravel: never, bedload-only). Acceptance: log-slope of C(z) within 10-15% of R.

### 8. Pile scour benchmarks

Roulund, Sumer, Fredsoe & Michelsen (2005), JFM 534:351-401 — THE primary benchmark:
- Scour experiment: circular pile D = 0.10 m, 3.85 m from upstream end of a 5.65 m x 3.60 m sand pit (0.15 m sand layer), flow depth h = 0.40 m, mean velocity V = 0.46 m/s (ReD = 4.4e4), d50 = 0.26 mm, angle of repose phi = 32 deg, V/Vcr = 1.25 -> live-bed. Measured equilibrium scour at upstream stagnation point = 12.48 cm -> S/D = 1.25, reached within ~2 h.
- Rigid-bed flow experiment/computation: k-omega SST, ~8e5 cells; parameter study delta/D = 0.02-100, ReD = 100-2e6.
- Their morphodynamic run (bedload + sand-slide only, no suspended load): equilibrium scour agrees well upstream; up to 30% discrepancy downstream. Sand slide triggered at bed slope > phi, relaxed to phi-2 deg.

Sumer, Christiansen & Fredsoe (1992) / Sumer & Fredsoe (2002) "The Mechanics of Scour in the Marine Environment":
- Steady-current live-bed equilibrium: S/D = 1.3, standard deviation sigma/D = 0.7 (design envelope S/D = 1.3 + 0.7 = 2.0, DNV-RP-C205 uses 1.3).
- Time scale: T = T* * D^2 / sqrt(g*(s-1)*d50^3), with T* = (1/2000)*(delta/D)*theta^(-2.2) for steady current (delta = boundary-layer thickness ~ h in a flume). Worked example (Roulund case, theta ~ 0.12, delta/D = 4): T* = 0.21 -> T ~ 130 s (scour development curve S_t = S*(1 - exp(-t/T))).

Melville & Coleman (2000) "Bridge Scour" / Melville & Chiew (1999):
- Clear-water design envelope: d_se = 2.4*D for narrow piers (D/h < 0.7); 2.0*sqrt(h*D) for 0.7 < D/h < 5; 4.5*h for D/h > 5. Max scour at V/Vc ~ 1.0.
- Time to equilibrium (clear-water): t_e(days) = 48.26*(D/V)*(V/Vc - 0.4) for h/D > 6; t_e = 30.89*(D/V)*(V/Vc - 0.4)*(h/D)^0.25 for h/D <= 6; valid V/Vc > 0.4.

Baykal et al. (2015, Phil. Trans. R. Soc. A 373:20140104) — k-omega, DTU model, current case = Sumer et al. (2013) Test 1: D = 40 mm, h = 0.08 m (h/D = 2), V = 0.413 m/s, U_f = 0.019 m/s, d50 = 0.17 mm, theta = 0.13, ReD = 1.7e4. Computed equilibrium S/D = 0.91 upstream / 0.52 downstream, in line with Sumer & Fredsoe data for h/D = 2. Switching OFF suspended load reduced upstream scour ~50% at theta = 0.13. Bed shear in the horseshoe-vortex region underpredicted ~30% vs Hjorth (1975) data.

## Grid resolution from published scour CFD (directly applicable)

- REEF3D::CFD, uniform Cartesian grid on the Roulund case (D = 0.1 m): 2.0 cm (= D/5) unstable; 1.5 cm (D/6.7) underpredicts; 1.0 cm (D/10) -> equilibrium 12.95 cm vs measured 12.48 cm (<4% error); 0.75 cm (D/13.3) -> 12.25 cm (2%) at 4.6x cost. 5 hydrodynamic updates per bed update sufficed (fully decoupled underpredicts early scour rate only). => For a voxel solver, D/10 per obstacle diameter is the working minimum, D/13+ preferred.
- TUDflow3D WALE-LES (de Wit, ISSMGE 2024): dx=dy=dz = D/40 near pile (D = 0.15 m -> 3.75 mm), domain 18D x 21D, ~14e6 cells; coarse D/25 grid deviates RMS 8% of final scour depth; morphological acceleration factor 100 accurate (factors 50/200/400 differ ~4%); 5 h lab scour in 5 days -> 0.7 days at factor 800.
- DTU body-fitted k-omega (Roulund 2005; Baykal 2015/2017): 96 cells around pile perimeter, min radial cell 0.02*D at pile, min vertical cell 0.01*D at bed, 14-16 cells over depth h = 2D, total ~1e5 cells (morphology) / 8e5 (rigid-bed flow).
- Structures to resolve: horseshoe vortex extends ~0.5*D upstream and ~0.3*D vertically in currents (Sumer et al. 1997; Roulund 2005) -> need >= 5-10 cells across each => near-bed dz <= 0.03-0.05*D.

## Practical guidance for our simulator

Test order and acceptance tolerances (informed by what published models achieve):
1. Poiseuille/laminar channel: error < 2%. 2. Cavity Re=100 at 64^3-128^3: RMS centerline error < 5% (semi-Lagrangian diffusion will fail Re=1000 unless BFECC/MacCormack advection is added — use Re=1000 as the diagnostic for advection-scheme upgrades, not as a gate). 3. Log-law channel with wall model at 2.5-5 cm cells: u_* within 10%, kappa within 5%. 4. Cylinder: Re=100 Cd 1.33-1.40, St 0.164-0.168 (gate); Re=1e4 LES Cd within 15%, St within 10%. 5. BFS Re<=400 vs Armaly. 6. Settling column: front speed = w_s within 1%. 7. Rouse: slope within 15%. 8. Roulund scour: with 1 cm voxels (D/10) target upstream equilibrium S/D = 1.25 within +/-15%; downstream within +/-30%; time scale within factor 2 (published models are consistently too fast — Baykal et al. 2017 report simulated time scales "significantly smaller" than experiment); scour-hole planform slope ~ phi = 32 deg. Include suspended load: bedload-only underestimates upstream scour by ~50% at theta >= 0.13. Bed shear amplification under the horseshoe vortex: expect to be ~30% low even with good models — do not tune sediment constants to compensate; validate shear on flat bed first (log-law) and accept alpha = tau/tau_infinity peak ~ 3-4 flanking the pile (11 at pile side in potential-flow limit is wrong; measured max amplification ~4-5, Hjorth 1975; Roulund 2005). At our full domain (10x10x5 m, 2.5-5 cm voxels), a 1 m-wide printed shape spans 20-40 cells — inside the validated D/10-D/40 window; features < 0.5 m (10 cells at 5 cm) will be under-resolved: refine to 2.5 cm or treat results as qualitative.

## Sources

- Ghia, Ghia & Shin (1982) "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method", J. Comp. Phys. 48:387-411. Data tables: https://gist.github.com/ivan-pi/3e9326d18a366ffe6a8e5bfda6353219 and https://gist.github.com/ivan-pi/caa6c6737d36a9140fbcf2ea59c78b3c
- Roulund, Sumer, Fredsoe & Michelsen (2005) "Numerical and experimental investigation of flow and scour around a circular pile", JFM 534:351-401. https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/numerical-and-experimental-investigation-of-flow-and-scour-around-a-circular-pile/17FB29E1C6E9351A41E39BAAE1B62C63
- Baykal, Sumer, Fuhrman, Jacobsen & Fredsoe (2015) "Numerical investigation of flow and scour around a vertical circular cylinder", Phil. Trans. R. Soc. A 373:20140104. https://pmc.ncbi.nlm.nih.gov/articles/PMC4275922/
- Baykal et al. (2017) "Numerical simulation of scour and backfilling processes around a circular pile in waves", Coastal Eng. 122:87-107. https://users.metu.edu.tr/cbaykal/manuscript_baykal_etal_CENG_2017.pdf
- REEF3D::CFD monopile scour grid-convergence study (ISSMGE ISFOG proc.): https://www.issmge.org/uploads/publications/108/137/Numerical_Modeling_of_Scour_Development_Around_a_Monopile_Foundation_Using_Computational_Fluid_Dynamics.pdf
- de Wit (ISSMGE proc.) "3D CFD LES process-based scour simulations with morphological acceleration": https://www.issmge.org/uploads/publications/108/123/3D_CFD_LES_process-based_scour_simulations_with_morphological_acceleration.pdf
- Sumer, Christiansen & Fredsoe (1992) "Time scale of scour around a vertical pile", Proc. 2nd ISOPE, 3:308-315; Sumer & Fredsoe (2002) "The Mechanics of Scour in the Marine Environment", World Scientific.
- Soulsby (1997) "Dynamics of Marine Sands", Thomas Telford (w_s eq. 102, theta_cr, log-law).
- Melville & Coleman (2000) "Bridge Scour", Water Resources Publications; Melville & Chiew (1999) "Time scale for local scour at bridge piers", J. Hydraul. Eng. 125(1):59-65.
- Armaly, Durst, Pereira & Schonung (1983) "Experimental and theoretical investigation of backward-facing step flow", JFM 127:473-496; Gartling (1990) IJNMF 11:953-967. Validation summary: https://doc.notus-cfd.org/dc/d98/group__doc__validation__backward__facing__step.html
- Sumer & Fredsoe (2006) "Hydrodynamics Around Cylindrical Structures", World Scientific (Cd/St vs Re); Williamson (1996) Ann. Rev. Fluid Mech. 28:477-539; Norberg (2003) J. Fluids Struct. 17:57-96.

## Verification

Adversarial fact-check (2026-07-07) of the eight key equations. Every functional form, constant, exponent and
derived number was re-derived numerically and checked against sources located independently of this note's citations.
Verdict per equation: CONFIRMED / CORRECTED / UNVERIFIABLE.

1. **Soulsby settling velocity — CONFIRMED.** Independent source (USACE HEC-RAS 2D Sediment technical reference,
   Particle Settling Velocity page) quotes Soulsby (1997) exactly as w_s = (nu/d)*[(10.36^2 + 1.049*Dstar^3)^(1/2) - 10.36],
   Dstar = d*[g(s-1)/nu^2]^(1/3). Recomputed with nu=1.05e-6 m2/s, s=2.585, g=9.81: Dstar = 24160.5*d(m);
   w_s = 0.00726 / 0.0356 / 0.1173 / 0.4028 m/s for d = 0.1 / 0.26 / 1 / 10 mm — all four table values match to the
   quoted rounding. Constant 1.049 is the natural (irregular) grain fit, as stated.
   Source: https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.5/model-description/fall-velocity-and-settling/particle-settling-velocity

2. **Soulsby-Whitehouse critical Shields — CONFIRMED.** Independent source (HEC-RAS Critical Thresholds page) quotes
   theta_cr = 0.30/(1 + 1.2*Dstar) + 0.055*[1 - exp(-0.02*Dstar)] — identical constants (0.30, 1.2, 0.055, 0.020).
   Recomputed: theta_cr = 0.0795 / 0.0416 / 0.0311 / 0.0556 for Dstar = 2.42 / 6.28 / 24.2 / 242 — matches the table
   (0.080 / 0.042 / 0.031 / 0.056). It is a fit to the Shields curve proposed for all Dstar (limit 0.30/(1+1.2 Dstar)
   -> 0.30 as Dstar -> 0).
   Source: https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.6/model-description/critical-thresholds-for-transport-and-erosion

3. **Rough-bed log law — CONFIRMED.** z0 = k_s/30 for rough beds (COHERENS sediment-model documentation eq. 7.6:
   z0 = 0.11*nu/u_* + k_b/30 ~ k_b/30 in the fully rough limit); k_s = 2.5*d50 is the standard Nikuradse grain-roughness
   convention used by Soulsby (1997)/Nielsen (1992); kappa = 0.40 (Soulsby) vs 0.41-0.412 (Nezu & Nakagawa) is the
   accepted spread. The depth-averaged form V/u_* = (1/kappa)*(ln(h/z0) - 1) is confirmed by the equivalent 2-D drag
   coefficient in COHERENS eq. 7.2: Cdb = [kappa/(ln(H/z0) - 1)]^2. Worked value recomputed: h=0.4 m, d50=0.26 mm ->
   z0 = 2.17e-5 m; u_*=0.02 m/s -> V = 0.441 m/s (note says 0.44). Validity limits (log layer z/h < 0.2; fully rough
   k_s+ = u_*k_s/nu > 70) are the standard Nezu & Nakagawa / Nikuradse criteria.
   Source: https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf

4. **Rouse profile — CONFIRMED.** COHERENS eq. (7.128) gives exactly C_eq/C_a = [((H-z)/z)*(a/(H-a))]^R with
   R = w_s/(beta*kappa*u_*), beta = inverse Prandtl-Schmidt number (sediment/momentum diffusivity ratio) — beta in the
   DENOMINATOR, as written here (van Rijn: beta = 1 + 2*(w_s/u_*)^2, capped at 2, hence the quoted range 1-2).
   Test case recomputed: w_s(0.1 mm) = 0.00726 m/s, u_* = 0.02 m/s, kappa = 0.40, beta = 1 -> R = 0.907 ~ 0.91. Match.
   Sources: https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf ; https://en.wikipedia.org/wiki/Rouse_number

5. **Sumer-Fredsoe scour depth & time scale — CONFIRMED.** Independent quote of S/D: "Sumer et al. (1992) stated that
   the mean value and standard deviation for the normalized equilibrium scour depth (S/D) for a vertical cylindrical
   pile in steady current are 1.3 and 0.7, respectively" (Mostafa & Agamy 2011, IJEST 3(11):8160-8178); DNV-OS-J101
   likewise adopts S/D = 1.3. Time scale independently quoted from Sumer et al. (2002) in an Aalborg University study
   (projekter.aau.dk): T* = (1/2000)*(delta/D)*theta^(-2.2), T = T* * D^2/(g(s-1)d50^3)^(1/2),
   s(t)/D = (s_e/D)*(1 - exp(-t/T)) — identical constants and exponent. Worked example recomputed: theta = 0.12,
   delta/D = 4 -> T* = 0.212, T = 128 s (~130 s as stated).
   Sources: https://www.idc-online.com/technical_references/pdfs/chemical_engineering/SCOUR%20AROUND%20SINGLE%20PILE%20AND%20PILE%20GROUPS%20SUBJECTED%20TO%20WAVES%20AND%20CURRENTS.pdf ;
   https://projekter.aau.dk/projekter/files/14765331/Scour_in_a_marine_env._cha._by_currents_and_waves_combined.pdf

6. **Melville clear-water envelope & equilibrium time — CONFIRMED.** Melville's own ISSMGE paper ("Local Scour Depths
   at Bridge Foundations: The New Zealand Methodology") Table 1 gives K_yb = 2.4b for b/y < 0.7; 2*(yb)^(1/2) for
   0.7 < b/y < 5; 4.5y for b/y > 5, and t_e(days) = 48.26*(b/V)*(V/Vc - 0.4) for y/b > 6;
   t_e(days) = 30.89*(b/V)*(V/Vc - 0.4)*(y/b)^0.25 for y/b <= 6, both valid V/Vc > 0.4 — identical to this note.
   Internal consistency check: 30.89*6^0.25 = 48.35 ~ 48.26 (the two branches join at y/b = 6). CAVEAT worth keeping in
   the spec: these t_e equations are CLEAR-WATER only (V/Vc <= 1); under live-bed conditions Melville sets K_t = 1
   (equilibrium reached quickly) — do not extrapolate t_e to V/Vc > 1.
   Source: https://www.issmge.org/uploads/publications/108/110/Local_Scour_Depths_at_Bridge_Foundations_-_New_Zealand_Methodology.pdf

7. **Cylinder Cd/St reference values — CONFIRMED.** Independent benchmark compilations agree: Re=40 steady twin-vortex
   flow (shedding onset Re ~ 47), Cd ~ 1.5-1.6; Re=100: mean Cd 1.32-1.39 across published solutions, St = 0.164-0.165
   (Williamson's experimental 0.164); Re=200: 2D numerics St ~ 0.196-0.20 while experiments give ~0.18-0.19 because
   mode-A three-dimensionality begins near Re ~ 190 (Williamson); Re=1000: Cd ~ 1.0, St ~ 0.21; subcritical
   Re = 1e4-1e5: Cd ~ 1.1-1.2 with St ~ 0.19-0.20 (St = 0.19 found for Re < 4.5e4, Norberg); drag crisis at
   Re ~ 2e5-5e5 with Cd dropping to ~0.2-0.3 minimum. Formula conventions St = fD/U and Cd = 2Fx/(rho U^2 D L)
   (i.e., Fx/(0.5 rho U^2 A), A = D*L) are standard.
   Sources: https://www.researchgate.net/figure/The-Strouhal-number-for-the-flow-past-a-circular-cylinder-at-Re-100_tbl1_309884026 ;
   https://web.stanford.edu/group/tfsa/TF_reports/TF-062_Beaudon.pdf

8. **Backward-facing step — CONFIRMED.** Independent source (arXiv:2507.16509, FVM/MHD BFS benchmark) confirms:
   Armaly Re uses hydraulic diameter D = 2h (twice inlet height) with mean inlet velocity; ER = 1.942; 2D primary
   reattachment x1/S = 3.00 at Re=100 and 4.95 at Re=200 (notus-cfd validation page: Armaly experimental value 5 at
   Re=200); Gartling (1990) Re=800 benchmark x1/S = 12.20 with upper-wall bubble length (x_rs - x_s)/S = 11.26 —
   consistent with the quoted separation 9.70 and reattachment 20.96 (20.96 - 9.70 = 11.26); Gresho et al. reproduce
   the same values. Three-dimensionality of the experiment above Re ~ 400 confirmed by both sources ("from Re=400 the
   flow is 3D"), which is why Armaly's measured x1/S (~14 at Re=800) exceeds all 2D solutions — compare 2D runs to
   Gartling, not to the experiment, at Re=800, exactly as this note instructs.
   Sources: https://arxiv.org/pdf/2507.16509 ; https://doc.notus-cfd.org/dc/d98/group__doc__validation__backward__facing__step.html

No corrections required: all 8 equation sets survived verification with constants, exponents, units, validity ranges
and worked numbers intact.
