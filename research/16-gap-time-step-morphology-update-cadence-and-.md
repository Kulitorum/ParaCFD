# Time-Step, Morphology-Update Cadence, and Run-Length/Ranking Protocol

## Overview

This note closes the numerical time-management gap: it fixes the fluid time step from *accuracy* (semi-Lagrangian MacCormack is stable at any CFL, so stability gives no dt), sets the Exner-update cadence and its shear-averaging window, gives a MORFAC schedule with validity checks, and reconciles scour equilibrium timescales with the 15-30 min/candidate RTX 4090 budget. The central finding: by Sumer & Fredsoe's own T* relation, equilibrium time T across our design matrix spans 66 s (2.5 m/s, 0.2 mm) to 50 h (0.5 m/s, 5 mm) — a factor 2,700. A fixed 300 s evaluation captures ~99% of equilibrium in the fastest case and <1% in the slowest; raw S(300 s) rankings would be timescale artifacts. The fix is a per-scenario run length t_run >= 0.5*T_pred plus exponential extrapolation to S_eq (error ~10-20% when fitting half the record, per Welzel et al. 2019), with MORFAC <= 10 needed only for the U = 0.5 m/s scenarios. Everything closes inside the budget on the 4 M-cell grid (h = 5 cm); the 32 M-cell grid (h = 2.5 cm) is for hero/validation runs only.

## Key equations & results

**E1. Fluid time step from advection accuracy (Selle, Fedkiw, Kim, Liu & Rossignac 2008, J. Sci. Comput. 35:350-371).**
`dt = CFL_adv * h / max|u|`, with `CFL_adv <= 2.0` for ranking runs and `<= 1.0` for hero/validation runs. Evidence: the semi-Lagrangian MacCormack scheme is measured 2nd-order in space and time at BOTH CFL = 0.75 (order 1.96-2.01, Table 2) and CFL = 1.75 (order 2.05-2.21, Table 3); the analytic per-step truncation error is `E(alpha)*dt^3*u_ttt/6` with `E_mac(alpha) = alpha^2/2 - 1` for alpha <= 1 (piecewise-quadratic continuation for alpha in (1,2], (2,3]; error magnitude oscillates but stays bounded to CFL = 3 in their plots). The accuracy killer at high CFL is *not* the truncation error but the mandatory reversion to 1st-order semi-Lagrangian wherever the forward/backward characteristics touch solid cells — at large CFL the backtraced rays cross the voxelized bed/obstacle more often, degrading exactly the near-bed region that sets tau_b. Hence cap CFL at 1-2, not 4-5. Use max|u| over the domain (horseshoe-vortex/contraction amplification makes local speed ~1.5-1.7x U_inf), recomputed every step.
Concrete dt (h = 5 cm, max|u| = 1.7*U): U = 0.5/1.0/1.5/2.5 m/s -> dt = 88/44/29/18 ms at CFL 1.5; at h = 2.5 cm halve these. This brackets the expected 10-40 ms at 2.5 m/s. Explicit Smagorinsky diffusion limit dt <= h^2/(6*nu_t) ~ 0.1-0.25 s and settling limit ws*dt <= h are never binding.

**E2. Cost per fluid step (memory-traffic model, RTX 4090, 1008 GB/s peak, ~75-85% achievable).**
`t_step ~ B_cell * N_cells / BW_eff`. Traffic budget B_cell: MacCormack advection of u,v,w,c (2 semi-Lagrangian sweeps each) ~400-600 B; Smagorinsky + forces ~120 B; warm-started MGPCG projection to 1e-4 (6-12 V-cycle iterations) ~1.4-2.2 kB; sediment/Exner amortized ~50 B. Total ~2-3 kB/cell/step. -> **~12-18 ms/step at 4 M cells (h = 5 cm); ~100-150 ms/step at 32 M cells (h = 2.5 cm)**. Wall-clock closure: 300 s of flow at U = 2.5, h = 2.5 cm, CFL 2 (dt = 20 ms) = 15,000 steps x ~0.12 s = ~30 min — matching the 15-30 min/300 s figure quoted for the LBM alternative (research/09). At h = 5 cm the same 300 s costs 2-4 min; **a 20-min budget buys ~80,000 steps**, i.e. 7,060/3,530/2,350/1,410 fluid-seconds at U = 0.5/1.0/1.5/2.5 m/s (CFL 1.5). Ranking must therefore run on the 4 M grid.

**E3. Bed-update cadence and shear averaging (practice in Delft3D/FLOW-3D/REEF3D; research/10).**
Update Exner every N flow steps using bed shear averaged over the window `t_avg = N*dt ~ 1-2 large-eddy turnovers D/U` (D = structure size 1-3 m -> t_avg = 1-3 s -> N ~ 10-100, consistent with the earlier notes). Implement as an exponential moving average, `tau_avg <- (1 - dt/t_avg)*tau_avg + (dt/t_avg)*tau_b`, so no history storage. Cadence N is a *compute-saving and de-noising* knob only — it does not accelerate morphological time (that is MORFAC's job). Adaptive N from the bed-change limiter: choose N (clip to [5,100]) such that predicted `max|dzb| = M * N * dt * max_rate <= 0.05*dx`; with S_eq ~ 1.3*D = 2.6 m and T = 66 s (U = 2.5, 0.2 mm) the initial rate S_eq/T ~ 0.04 m/s forces N ~ 4-10, while clear-water cases allow N = 100.

**E4. Exner with MORFAC and validity limits (Reyns, Dastgheib, Ranasinghe, Luijendijk, Walstra & Roelvink, ICCE 2014; Ranasinghe et al. 2011, Coastal Eng. 58:806-811).**
`(1 - p) * dzb/dt = M * [ -div(qb) - E + D ]`, p = 0.36. Constraints: (a) per-update bed change `|M*dzb| <= 0.05-0.1*dx` AND `<= 0.1*h_local`; (b) morphological celerity CFL `c_bed * M * N * dt / dx <= 0.9`, c_bed ~ 3*qb/((1-p)*h_dep); (c) critical MORFAC drops by MORE THAN AN ORDER OF MAGNITUDE in reversing/tidal flow vs unidirectional (Reyns et al.: <200 vs ~4000 at Fr = 0.11, coastal grid dx = 7.5-30 m; critical MORFAC decreases exponentially with Fr, scales quasi-linearly with dx, dt secondary). At our dx = 0.025-0.05 m the safe envelope is far smaller: **M <= 10 steady current, M <= 5 once tidal reversal/waves are added** (Delft3D local-scour guidance). Verification: for one candidate per scenario run M = 1 vs M = target and require Brier Skill Score >= 0.95 (Reyns et al. used BSS 0.99 cut-off on idealized cases; van Rijn et al. 2003 call 0.5 "sufficient" for complex cases).

**E5. Scour timescale (Sumer & Fredsoe 2002, The Mechanics of Scour in the Marine Environment):**
`S(t) = S_eq * (1 - exp(-t/T))`; `T = T* * D^2 / sqrt(g*(s-1)*d50^3)`; steady current `T* = (1/2000) * (delta/D) * theta^(-2.2)` (delta = boundary-layer thickness ~ water depth 5 m). Alternatives: Whitehouse 1998 `T* = 0.014 * theta^(-1.29)`; **Larsen & Fuhrman 2023 (Coastal Eng. 185:104356) argue on transport-rate grounds T* ~ theta^(-3/2)** — flatter than theta^-2.2, so Sumer-Fredsoe *overestimates* T at low theta; treat the E6 table's clear-water T values as upper bounds (conservative for run-length planning).

**E6. Worked run-length table** (delta = 5 m, D = 2 m reference, s = 2.585, g = 9.81 m/s^2, nu = 1.05e-6 m^2/s; u_f = 0.4*U/(ln(h/z0)-1), z0 = d50/12, h = 5 m; theta = u_f^2/(g*(s-1)*d50); theta_cr from Soulsby-Whitehouse with D* = 4.8/8.5/24.2/120.8):

| d50 [mm] | U [m/s] | theta | regime (theta/theta_cr) | T [s] | t(0.5T) [s] | fluid s in 20 min | M needed |
|---|---|---|---|---|---|---|---|
| 0.20 | 0.5 | 0.095 | live-bed (1.9) | 78,800 | 39,400 | 7,060 | 6 |
| 0.20 | 1.0 | 0.382 | live-bed (7.8) | 3,730 | 1,870 | 3,530 | 1 |
| 0.20 | 1.5 | 0.859 | live-bed (17) | 627 | 313 | 2,350 | 1 (reaches 3.7T) |
| 0.20 | 2.5 | 2.385 | sheet flow (48) | 66 | 33 | 1,410 | 1 (reaches >5T) |
| 0.35 | 0.5 | 0.060 | live-bed (1.7) | 93,800 | 46,900 | 7,060 | 7 |
| 0.35 | 1.0 | 0.241 | live-bed (6.8) | 4,440 | 2,220 | 3,530 | 1 |
| 0.35 | 1.5 | 0.542 | live-bed (15) | 746 | 373 | 2,350 | 1 |
| 0.35 | 2.5 | 1.504 | live-bed (42) | 79 | 40 | 1,410 | 1 |
| 1.0 | 0.5 | 0.026 | clear-water (0.83) | 126,100 | 63,100 | 7,060 | 9 |
| 1.0 | 1.0 | 0.103 | live-bed (3.3) | 5,970 | 2,990 | 3,530 | 1 |
| 1.0 | 1.5 | 0.231 | live-bed (7.4) | 1,000 | 500 | 2,350 | 1 |
| 1.0 | 2.5 | 0.643 | live-bed (21) | 106 | 53 | 1,410 | 1 |
| 5.0 | 0.5 | 0.007 | no motion (0.14) | — | — | — | skip (short check) |
| 5.0 | 1.0 | 0.029 | clear-water (0.56) | 8,520 | 4,260 | 3,530 | 2 |
| 5.0 | 1.5 | 0.066 | live-bed (1.3) | 1,430 | 715 | 2,350 | 1 |
| 5.0 | 2.5 | 0.183 | live-bed (3.5) | 151 | 76 | 1,410 | 1 |

T scales as D^2 * (delta/D) = D*delta: for D = 1 m halve T, for D = 3 m multiply by 1.5. The gap statement is confirmed: T ~ 500-750 s in energetic live-bed, 22-50 h in the 0.5 m/s clear-water/marginal cases.

**E7. Equilibrium extrapolation reliability (Welzel, Schendel, Hildebrandt & Schlurmann 2019, Coastal Eng. 152:103515, verified full text).**
They fitted `S(t) = a * (1 - 1/(1 + b*t))` (hyperbolic; chosen over Sheppard et al. 2004's form for robustness on noisy records) by least squares, then re-fitted using ONLY THE FIRST 50% of each record: **average |S_end - S_fit,end50| = 18.6% over all tests, 8.7% after excluding 3 tests contaminated by global scour**; their records ran to >= 90% of equilibrium (i.e. record length ~2.3T, so the half-record fit window is ~1.15T). Sheppard, Ozalp et al. equilibrium-rate criterion (lab): quasi-equilibrium when growth < 5% of D per 24 h; translate for the simulator as `(dS/dt) * T_fit / S_eq_fit < 0.05` over the last quarter of the run.

**E8. Spin-up and MORFAC ramp.** Domain flow-through time L/U = 10/U = 4-20 s. Protocol: (i) freeze morphology, run 3-5 flow-throughs (12-100 s fluid time; <2 min wall) until KE and mean tau_b plateau; (ii) enable Exner with M ramped linearly 1 -> M_target over 1-2 flow-throughs (avoids shocking the bed with accelerated transients, per Deltares practice and Reyns et al.'s spin-up caution); (iii) start the S(t)/V(t) recording clock at the end of the ramp.

## Practical guidance for our simulator

**Fluid dt:** auto-step `dt = 1.5 * h / max|u|` (ranking) with hard bounds [5 ms, 250 ms]; drop the prefactor to 1.0 (and h to 2.5 cm) for hero runs. Do not exploit unconditional stability beyond CFL 2 — near-boundary 1st-order reversion and LES temporal resolution (dt/(h/u') ~ 0.3 at CFL 2 already) both degrade, and dt savings are irrelevant because the budget closes at CFL 1.5.

**Cadence + MORFAC are orthogonal knobs.** N (10-100, adaptive from the 0.05*dx limiter, EMA window 1-3 s) removes turbulence noise from tau_b and cuts sediment-module cost; M multiplies morphological time. Combined morphological step = M*N*dt; enforce `M*N*dt*max_rate <= 0.05*dx` first by shrinking N, then, if N hits 5, by shrinking M (adaptive MORFAC). Log the fraction of updates where the limiter clips — if >1%, results are rate-distorted; rerun with lower M.

**Per-scenario recipe (D = 2 m; grid 4 M cells, h = 5 cm; ~15 ms/step):** run the table's t(0.5T) of morphological time at the listed M (M = 6-9 only at U = 0.5 m/s; M = 1-2 elsewhere), i.e. **every cell of the design matrix fits in 15-25 min wall-clock**. Fast live-bed cases (U >= 1.5) reach 2-5T — report S_eq directly from the plateau. For U = 0.5 clear-water cases, if Larsen-Fuhrman theta^-3/2 scaling is adopted, T shrinks ~3-5x and M = 2-3 suffices — run one M = 1 vs M = 6 BSS check before trusting the accelerated results.

**Ranking protocol:** (1) record S(t) (max scour depth in a 1D-annulus around the shape) and trapped volume V(t) (sum of f_pack inside/around the shape) in *morphological* time; (2) fit BOTH `S_eq*(1-exp(-t/T))` and Welzel's `a*(1-1/(1+b*t))` by nonlinear least squares; if the two S_eq estimates disagree by >15%, extend the run (adaptive stopping); (3) rank on fitted S_eq (lower = better protection) and V_eq (higher = better self-ballasting), carrying the fit 95% CIs — two shapes are distinguishable only if CIs do not overlap; (4) never compare raw S at fixed 300 s across scenarios, and always compare candidates at identical t_run/T within a scenario. Expected extrapolation error: ~10-20% at t_run = 0.5T (Welzel), <10% at t_run >= 1T; do not fit records shorter than 0.3T (fit becomes degenerate: S_eq and T trade off along a ridge — check the fit covariance).

**Do not use a fixed 300 s.** 300 s is 4.5T at (2.5 m/s, 0.2 mm) but 0.004T at (0.5 m/s, 1 mm). The run-length table replaces it.

## Sources

- Selle, Fedkiw, Kim, Liu & Rossignac (2008), "An Unconditionally Stable MacCormack Method", J. Sci. Comput. 35:350-371. https://faculty.cc.gatech.edu/~jarek/papers/maccormack.pdf (full text verified: Tables 2-3 convergence orders at CFL 0.75/1.75; E_mac(alpha) truncation analysis; boundary reversion to 1st order)
- Reyns, Dastgheib, Ranasinghe, Luijendijk, Walstra & Roelvink (2014), "Morphodynamic upscaling with the MORFAC approach in tidal conditions: the critical MORFAC", Coastal Engineering Proceedings (ICCE 2014). https://icce-ojs-tamu.tdl.org/icce/index.php/icce/article/download/7269/pdf_702/0 (full text verified: critical MORFAC <200 tidal vs ~4000 unidirectional; BSS 0.99 criterion; dx dominance)
- Ranasinghe, Swinkels, Luijendijk, Roelvink, Bosboom, Stive & Walstra (2011), "Morphodynamic upscaling with the MORFAC approach: Dependencies and sensitivities", Coastal Eng. 58:806-811. https://www.sciencedirect.com/science/article/abs/pii/S0378383911000457
- Sumer & Fredsoe (2002), The Mechanics of Scour in the Marine Environment, World Scientific. (T* = (1/2000)(delta/D)theta^-2.2; S(t) exponential; constants as independently verified in research/05.)
- Larsen & Fuhrman (2023), "Re-parameterization of equilibrium scour depths and time scales for monopiles", Coastal Eng. 185:104356. https://doi.org/10.1016/j.coastaleng.2023.104356 (open access; T* ~ theta^(-3/2) scaling argument)
- Welzel, Schendel, Hildebrandt & Schlurmann (2019), "Scour development around a jacket structure in combined waves and current conditions compared to monopile foundations", Coastal Eng. 152:103515. https://doi.org/10.15488/9846 (full text verified via Hannover repository: Eq. (4) hyperbolic fit; 18.6%/8.7% half-record extrapolation error; S_end = mean of last 25%)
- Sheppard, Odeh & Glasser (2004), "Large scale clear-water local pier scour experiments", J. Hydraul. Eng. 130(10):957-963. (equilibrium-rate criterion; alternative fit form referenced by Welzel)
- Delft3D-FLOW User Manual, Deltares, ch. 11 (MORFAC guidance, bed-change-vs-depth limit). https://content.oss.deltares.nl/delft3d4/Delft3D-FLOW_User_Manual.pdf
- Soulsby (1997), Dynamics of Marine Sands, Thomas Telford (Soulsby-Whitehouse theta_cr(D*); log-law depth-averaged friction used for the table).
- FluidX3D benchmark data for RTX 4090 bandwidth context (research/09). https://github.com/ProjectPhysX/FluidX3D
