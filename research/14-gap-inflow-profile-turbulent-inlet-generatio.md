# Inflow Profile, Synthetic Turbulent Inlet, and Top/Lateral Boundary Conditions

## Overview

The solver notes (08) fix the outlet (Orlanski + global flux rescale + p = 0) and obstacle BCs but leave the inlet, top and lateral faces unspecified. This note closes that gap for the 10 x 10 x 5 m domain (dx = 2.5-5 cm, U_d = 0.5-2.5 m/s, d50 = 0.1-10 mm). Decisions in one line each:

1. **Inlet mean profile**: rough log law u(z) = (u*/0.40) ln(z/z0), z0 = d50/12, with u* mapped from the depth-averaged target U_d over the 5-m domain height; **delta := h_dom = 5 m** (real tidal BLs are 20-45 m thick, so the whole domain sits inside the log layer; use delta = 5 m in the T* = (1/2000)(delta/D) theta^-2.2 timescale of notes 05/11).
2. **Inlet fluctuations**: classic **Synthetic Eddy Method (Jarrin et al. 2006)** with Nezu-Nakagawa (1993) Reynolds-stress profiles, sigma = h/5 = 1.0 m, N ~ 150 eddies (~200 lines of CUDA). Divergence errors are removed for free by the existing MGPCG projection, so the divergence-free SEM variant (Poletto 2013) is optional. For the Roulund benchmark gate and hero runs, replace SEM by a **precursor periodic channel** driven by g_x = u*^2/h (the LES analogue of Baykal et al. 2015/2017 and exactly what Kirkil & Constantinescu did).
3. **Top**: rigid **free-slip lid** (w = 0, d(u,v)/dz = 0, Neumann p and scalars) — the standard in every published scour CFD (Roulund 2005; Liang & Cheng 2005; Baykal 2015/2017). Fr = 0.36 at 2.5 m/s is the marginal upper edge; fine for ranking, flag hero runs above 2.0 m/s.
4. **Lateral**: free-slip walls (Baykal 2017 verbatim), frontal blockage <= 6%, else correct with continuity/Maskell.
5. **Tidal reversal**: swap inlet/outlet faces through a cosine ramp of ~60 s (compressed slack water); do NOT body-force the obstacle domain (wake recycling). Morphology stays live through reversals (transport ~0 near slack anyway); freeze Exner only during initial spin-up.

## Key equations & results

### 1. Inlet mean-velocity profile (log law) and the U_d -> u* mapping
u(z) = (u*/kappa) * ln(z/z0), kappa = 0.40, z0 = d50/12 (= k_s/30 with k_s = 2.5 d50), valid z0 < z <= h_dom = 5 m. [Soulsby 1997, Dynamics of Marine Sands; confirmed in HR Wallingford TR137]

Depth-averaging the log law over [z0, h_dom] gives the mapping from the target depth-averaged speed U_d:
**u* = kappa * U_d / (ln(h_dom/z0) - 1)**, h_dom = 5 m. [Soulsby 1997, eq. 37 form; same identity used in note 11]

Worked values (u*/U_d for h_dom = 5 m, z0 = d50/12):

| d50 (mm) | z0 (m)   | u*/U_d | u* at U_d = 1 m/s | tau_0 = rho u*^2 (Pa, rho = 1025) |
|----------|----------|--------|-------------------|------------------------------------|
| 0.1      | 8.33e-6  | 0.0325 | 0.033 m/s         | 1.09                               |
| 0.2      | 1.67e-5  | 0.0345 | 0.034 m/s         | 1.22                               |
| 1.0      | 8.33e-5  | 0.0400 | 0.040 m/s         | 1.64                               |
| 10       | 8.33e-4  | 0.0519 | 0.052 m/s         | 2.77                               |

Scale linearly in U_d (tau_0 quadratically). Note this grain-roughness u* is lower than the field C_D = 0.0025 value (u* = 0.05 U_bar, note 12) because that default assumes a *rippled* bed (z0 = 6 mm). Use z0 = d50/12 at the inlet: the LES + Exner resolves bedform drag itself; only substitute z0 = 6 mm if ripples are declared subgrid for a given run (then u*/U_d = 0.070, tau_0 = 5.0 Pa at 1 m/s). The profile cap at the lid is u(5 m) = 1.086-1.15 * U_d. Implement u(z) at cell centers; first cell center z = 0.0125-0.025 m >> z0, so no singularity handling needed.

### 2. Inlet Reynolds-stress / TI profiles (targets for the turbulence generator)
Nezu & Nakagawa (1993), open-channel universal profiles (z measured from bed, h = h_dom = 5 m):
- u'_rms/u* = 2.30 exp(-z/h) (streamwise)
- v'_rms/u* = 1.63 exp(-z/h) (spanwise)
- w'_rms/u* = 1.27 exp(-z/h) (vertical)
- k/u*^2 = 4.78 exp(-2 z/h); -u'w' = u*^2 (1 - z/h); other off-diagonals 0.

Resulting ambient streamwise TI (= u'_rms/U_d) for 0.2 mm sand: 7.2% at z = 0.5 m, 4.8% at mid-depth — matching note 12's "ambient 5-10%". To hit the **amplified near-pile target TI = 10-20%** (note 12), multiply the whole R_ij tensor by gamma^2 with gamma = TI_target * U_d / (2.081 u*) (2.081 = 2.30 e^-0.1, anchoring at z = 0.1 h). For 0.2 mm sand and TI = 15%: gamma = 2.1. **Cap gamma at 2.5**: without matching mean shear the excess decays within ~1-2 h of fetch.

### 3. Synthetic Eddy Method (Jarrin et al. 2006) — the inlet fluctuation generator
u_i(x,t) = U_i(z) + (1/sqrt(N)) * SUM_{k=1..N} a_ij * eps_j^k * f_sigma(x - x^k(t))

- f_sigma(x - x^k) = sqrt(V_B / sigma^3) * f((x-x^k)/sigma) * f((y-y^k)/sigma) * f((z-z^k)/sigma), tent function **f(r) = sqrt(3/2) * (1 - |r|) for |r| < 1, else 0**.
- a_ij = Cholesky factor of R_ij: a11 = sqrt(R11); a21 = R21/a11; a22 = sqrt(R22 - a21^2); a31 = R31/a11; a32 = (R32 - a21 a31)/a22; a33 = sqrt(R33 - a31^2 - a32^2). (With our R: a21 = a32 = 0, a31 = -u*^2(1-z/h)/a11.)
- eps_j^k = independent random signs (+1/-1); eddy centers x^k uniform in box B = inlet plane extruded +/- sigma in all directions; convected x^k += U_c dt with **U_c = U_d** (NASA used U_inf/2; anything O(U) works); an eddy leaving B is regenerated at the upstream face with new random (y, z, eps).
- **N ~ V_B/sigma^3** (equivalently N = C * A_inlet/A_eddy with C = 1, THV review eq. 58). For sigma = 1 m: V_B = (2)(12)(6) = 144 m^3 -> **N = 150**.
- Length scale: **sigma = h/5 = 1.0 m** (NASA TM-2018-219966: L = delta/5 reproduced measured Cf = 2.1e-3 within 3%; delta/10 and delta/20 attenuated stresses and biased Cf by 15%). 20-40 cells per sigma at dx = 2.5-5 cm — well resolved. Optionally taper sigma(z) = max(0.41 z, 0.2 m) below z = 2.4 m to mimic near-bed structure shrinking (Nezu-Nakagawa L_x ~ h (z/h)^1/2).
- Artifact to know: SEM injects a spurious spectral peak at f_SEM = U_c/sigma (= 1-2.5 Hz here); it vanishes downstream [NASA TM-2018-219966].
- SEM is **not divergence-free**; our MGPCG projection removes the divergent part each step (this is exactly what DFSEM otherwise fixes; DFSEM's benefit is mainly reduced inlet pressure noise [Poletto et al. 2013]). Keep classic SEM.

### 4. Adjustment (development) length — the reason a precursor exists
Both SEM and digital filtering need **~6 delta of fetch** before mean profile, Reynolds stresses and wall shear are self-consistent (measured via dCf/dx; independent of sigma) [NASA TM-2018-219966]. With delta = 5 m that is 30 m >> our 10 m box: SEM-fed runs carry an O(15%) wall-shear uncertainty at the structure (5-6 m fetch). Kirkil & Constantinescu instead fed instantaneous planes from a **precalculated periodic-channel LES** (Re = 18,000) — fully developed turbulence with essentially zero adjustment length; they showed a steady-profile-only inlet fundamentally changes HSV bimodality. Baykal et al. (2015, 2017) is the RANS analogue: "flow driven by a constant horizontal pressure gradient based on the desired undisturbed friction velocity, U_f = [-(h/rho)(dp/dx)]^0.5", i.e. body force
**g_x = u*^2 / h_dom** (= 2.4e-4 m/s^2 for u* = 0.0345 m/s, h = 5 m),
in a pile-free domain; the developed u (and k, omega) profiles become the inlet of the pile run.

**Precursor recipe (ours)**: same solver, no obstacle, periodic in x and y, free-slip lid, rough-wall bed, g_x = u*^2/h plus a mass-flux PI controller (g_x^(n+1) = g_x^n + 0.1 (u*^2/h) * (U_d - U_bulk^n)/U_d) to lock U_bulk = U_d. Spin-up ~10 eddy turnovers h/u* (~1100 s at U_d = 1 m/s) starting from log profile + Nezu-scaled white noise; then record the inlet plane at 5-10 Hz for ~300 s and replay (linear interpolation in time, loop the library). One library per (U_d, d50) pair, reusable across all shape iterations.

### 5. Top boundary: rigid free-slip lid
w = 0, du/dz = dv/dz = 0, dp/dn = dc/dn = 0. Verbatim standard: "the top boundary ... is treated as a frictionless slip wall, with ... zero gradient" and zero wall-normal velocity [Baykal et al. 2017 §2.2; same in Roulund 2005, Liang & Cheng 2005]. Validity: **Fr = U_d/sqrt(g h_dom)** = 0.14 (1 m/s), 0.29 (2 m/s), **0.36 (2.5 m/s)**; neglected surface deformation ~ U_d^2/(2g) = 0.05 / 0.20 / 0.32 m = 1-6.4% of depth. Roulund's flume benchmark ran at Fr = 0.23 under a lid. Accept the lid for the whole range; treat U_d > 2 m/s results as ranking-only (bed shear near the structure locally biased high by a few % because the lid suppresses the drawdown).

### 6. Lateral boundaries and blockage
Free-slip walls, Neumann for scalars [Baykal et al. 2017 §2.2]. Keep the structure's lateral clearance >= 3.5 structure widths. Blockage ratio **BR = A_frontal/(W h) = A_frontal/50 m^2 <= 0.06** (A_frontal <= 3 m^2, e.g. 3 m wide x 1 m tall unit). If exceeded: mean-flow continuity correction U_eff = U_d/(1 - BR), or Maskell for drag: C_D,corrected = C_D,measured/(1 + theta C_D,measured BR), **theta = 0.96** for high-aspect bluff bodies [Maskell; Cambridge Aeronautical J. extension]. Periodic lateral BCs only in the precursor (they would alias the structure's wake laterally in the obstacle domain).

### 7. Tidal reversal
Do **not** body-force the obstacle domain (streamwise periodicity would recycle the wake into the inlet). Swap faces instead:
U_d(t) = U_max * s(t), s(t) ramping +1 -> -1 with a cosine over **T_ramp = 60 s** (>= 3-6 flow-through times L_x/U_d = 10-20 s at 0.5-1 m/s) — a compressed slack water. During the ramp: inlet profile and SEM amplitude scale with |U_d(t)| (u*(t) = 0.40|U_d(t)|/(ln(h/z0)-1)); when s crosses 0, the x- face becomes the SEM/precursor inlet and the x+ face becomes Orlanski. Clip the Orlanski phase speed C to [0, dx/dt] so it degenerates gracefully at slack; the existing global flux rescale keeps the projection compatible throughout. Morphology may stay live (theta < theta_cr for |U_d| < ~0.4 m/s on 0.2 mm sand, so slack transport is negligible); freeze Exner only during initial warm-up — Baykal et al. (2017) froze the bed for a 10-forcing-period warm-up (their oscillatory inlet used u = U_m sin(2 pi t/T_w) with k_m = 0.0005 U_m^2, i.e. TI ~ 1.8% — too low for us, reinforcing the SEM choice).

## Practical guidance for our simulator

- **Default (shape-ranking runs)**: log-law inlet (eq. 1 table) + SEM (eq. 3) with sigma = 1.0 m, N = 150, R_ij from eq. 2 with gamma chosen for TI 10-15%. CUDA cost: one kernel over the 200x100-400x200 inlet plane, each thread summing <= N compact-support eddies — <0.1 ms/step; eddy state (x^k, eps^k: ~150 x 6 floats) updated host-side. Place the printed unit at x >= 5 m from the inlet.
- **Benchmark gate (Roulund) and hero runs**: precursor periodic channel (eq. 4) — removes the 6-delta adjustment bias that would otherwise sit inside the ~15-20% acceptance band on tau_0 and scour depth. For Roulund's flume, set h_dom = delta = 0.4 m, V = 0.46 m/s, d50 = 0.26 mm, and verify recovered u* = 0.02 m/s (note 11 gate) before enabling morphology.
- **Consistency checks to code**: (i) empty-domain SEM run must hold the log profile and Nezu stresses within 20% at x = 5 m; (ii) flat-bed tau_0 from the wall model must match rho u*^2 from the table within 20% (links to note 12's validation target); (iii) with the lid on, verify depth-mean flux conservation after each Orlanski + rescale step to 1e-6.
- **delta for the scour timescale**: use delta = 5 m (domain height) in T* = (1/2000)(delta/D) theta^-2.2; for benchmark flumes use their depth.
- TI is defined as u'_rms/U_d throughout, matching notes 05/12.

## Sources

- Soulsby, R.L. (1997) *Dynamics of Marine Sands*; log law, z0 = d50/12, depth-average identity — via Soulsby & Clarke, HR Wallingford TR137: https://eprints.hrwallingford.com/558/1/TR137.pdf
- Jarrin, N., Benhamadouche, S., Laurence, D., Prosser, R. (2006) "A synthetic-eddy-method for generating inflow conditions for large-eddy simulations", Int. J. Heat Fluid Flow 27:585-593: https://www.sciencedirect.com/science/article/abs/pii/S0142727X06000282
- Poletto, R., Craft, T., Revell, A. (2013) "A New Divergence Free Synthetic Eddy Method...", Flow Turbul. Combust. 91:519-539: https://www.researchgate.net/publication/256089177
- Smirnov, A., Shi, S., Celik, I. (2001) "Random Flow Generation Technique for LES and Particle-Dynamics Modeling", ASME J. Fluids Eng. 123(2):359-371: https://asmedigitalcollection.asme.org/fluidsengineering/article-abstract/123/2/359/463161
- NASA TM-2018-219966, "Evaluation of Inflow Turbulence Methods in LES" (SEM vs digital filtering; L = delta/5; 6-delta adjustment length): https://ntrs.nasa.gov/api/citations/20180006518/downloads/20180006518.pdf
- Baykal, C., Sumer, B.M., Fuhrman, D.R., Jacobsen, N.G., Fredsoe, J. (2017) "Numerical simulation of scour and backfilling processes around a circular pile in waves", Coastal Eng. (author manuscript, BC section quoted): https://users.metu.edu.tr/cbaykal/manuscript_baykal_etal_CENG_2017.pdf
- Baykal, C. et al. (2015) "Numerical investigation of flow and scour around a vertical circular cylinder", Phil. Trans. R. Soc. A 373:20140104.
- Kirkil, G., Constantinescu, G. (2010) "Flow and turbulence structure around an in-stream rectangular cylinder with scour hole", Water Resour. Res. 46:W11549 (precursor periodic-channel inflow, Re = 18,000): https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010WR009336
- Nezu, I., Nakagawa, H. (1993) *Turbulence in Open-Channel Flows*, IAHR Monograph, Balkema (u'/u* = 2.30 e^-z/h etc.): https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/turbulence-in-openchannel-flows-by-i-nezu-and-h-nakagawa-balkema-1993-286-pp-85-or-55/7BBAA5824CF1848C906D02E14633E869
- "Triple Hill's Vortex Synthetic Eddy Method" (SEM review; N = C A_in/A_eddy, C = 1), arXiv:2111.13757: https://arxiv.org/pdf/2111.13757
- Maskell blockage theory (theta = 0.96) — "Extensions to Maskell's theory for blockage effects on bluff bodies in a closed wind tunnel", Aeronaut. J.: https://www.cambridge.org/core/journals/aeronautical-journal/article/abs/extensions-to-maskells-theory-for-blockage-effects-on-bluff-bodies-in-a-closed-wind-tunnel/10AB99FE436F4318F00251BD609C178D
- Orlanski, I. (1976) "A Simple Boundary Condition for Unbounded Hyperbolic Flows", J. Comput. Phys. 21:251-269 (outlet, per note 08): https://www.sciencedirect.com/science/article/abs/pii/0021999176900231
