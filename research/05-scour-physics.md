# Scour around monopiles: physics & empirical results

## Overview

Local scour at a vertical circular pile in an erodible bed is driven by three flow features that amplify the undisturbed bed shear stress tau_0 (Sumer & Fredsoe 2002; Roulund et al. 2005):

1. **Downflow + horseshoe vortex (HSV)** — the adverse pressure gradient upstream of the pile drives a downflow on the front face; the incoming boundary layer separates and rolls up into a horseshoe vortex wrapping around the pile base. This is the *dominant* scour agent in steady current.
2. **Contraction of streamlines** at the pile flanks — highest bed-shear amplification in current: alpha = tau/tau_0 up to O(10) (measured up to 11 by Hjorth 1975) near phi ~ 45 deg from the upstream stagnation line. Under the HSV in front, alpha ~ 4-5.
3. **Lee-wake vortex shedding** — shed vortices sweep sediment into their cores and carry it downstream. In *waves* this is the dominant mechanism (HSV barely forms in oscillatory flow at low KC); max amplification there is O(4) at the side edges (Sumer, Fredsoe et al., JFM 1997).

Roulund et al. (2005, JFM 534:351-401, k-omega RANS + experiment, pile D = 0.1 m) showed the HSV size and bed-shear amplification **increase with boundary-layer-thickness ratio delta/D**, increase with ReD, and decrease with bed roughness; for very small delta/D (O(0.1)) or low ReD the HSV may not form at all. Scour starts at the HSV in front, the hole then wraps around; the equilibrium upstream slope equals the sediment friction angle (~32 deg for sand), the downstream slope is much flatter (~16-20 deg). Baykal et al. (2015, Phil Trans R Soc A 373:20140104) further showed with the same benchmark (ReD = 1.7e4, D = 4 cm, d50 = 0.17 mm, theta = 0.13): disabling **suspended load** cuts the upstream equilibrium scour depth by ~50% (S/D 0.91 -> 0.46), and disabling **vortex shedding** turns downstream scour into accretion — so a morphodynamic solver needs bedload + suspended load + resolved/parameterized wake unsteadiness to get both sides of the hole right.

**Regimes.** With U = depth-averaged (or delta-averaged) velocity and U_cr the threshold velocity for sediment motion (theta_cr ~ 0.05):
- U/U_cr < 0.5: no scour (even amplified local shear stays subcritical).
- 0.5 < U/U_cr < 1 (theta < theta_cr): **clear-water scour** — no ambient transport; scour grows slowly, no sediment supply into the hole; deepest scour occurs near U/U_cr = 1.
- U/U_cr > 1 (theta > theta_cr): **live-bed scour** — general transport everywhere; equilibrium is a dynamic balance (hole is continuously refilled and re-eroded); equilibrium depth nearly independent of theta, S/D ~ 1.3.

## Key equations & results

All units SI. D = pile diameter [m], S = equilibrium scour depth below ambient bed [m], d50 = median grain size [m], s = rho_s/rho (= 2650/1025 = 2.59 for quartz in seawater), g = 9.81 m/s^2, theta = Shields parameter = Uf^2 / (g (s-1) d50), theta_cr ~ 0.05, Uf = bed friction velocity [m/s].

**E1. Equilibrium scour depth, steady current, live-bed (Sumer, Fredsoe & Christiansen 1992; Sumer & Fredsoe 2002):**
- S/D = 1.3 (mean), standard deviation sigma(S/D) = 0.7.
- Valid: slender pile (D/h < ~0.7, i.e. h/D > ~1.4 so depth does not limit scour), live-bed (theta > theta_cr), delta/D not too small (HSV fully developed; Roulund et al. 2005 show scour depth drops as delta/D -> O(0.1)). Weak ReD dependence in the fully turbulent HSV regime.
- Design practice (e.g. DNV-OS-J101 / DNVGL-ST-0126) uses S/D = 1.3 for monopiles when no site data exist; a conservative envelope is mean + 1 sigma = 2.0.
- Melville & Coleman (2000) envelope for slender piers, uniform sand: S <= 2.4 D (clear-water peak at U/U_cr = 1 can exceed the live-bed mean).

**E2. Flow-intensity ramp, clear-water regime (Breusers, Nicollet & Shen 1977; as recited in Whitehouse 1998, Sumer & Fredsoe 2002):**
- S/D = 2.0 * tanh(h/D) * f(U/U_cr), with f = 0 for U/U_cr <= 0.5; f = (2 U/U_cr - 1) for 0.5 < U/U_cr < 1; f = 1 for U/U_cr >= 1.
- h = water depth. Use this to blend from zero scour to the full live-bed value; Melville & Coleman's flow-intensity factor K_I = U/U_cr (capped at 1) is an equivalent linear alternative.

**E3. Scour in waves (Sumer, Fredsoe & Christiansen 1992, J. Waterway Port Coastal Ocean Eng 118(1):15-31):**
- S/D = 1.3 * [1 - exp(-0.03 * (KC - 6))], for KC >= 6; S = 0 for KC < 6 (no vortex shedding, no HSV -> no scour at a circular pile).
- KC = Um * Tw / D = 2*pi*A/D, Um = near-bed orbital velocity amplitude [m/s], Tw = wave period [s], A = orbital excursion amplitude [m].
- Valid: live-bed (theta > theta_cr), regular waves, 6 <= KC <= ~100 (tests 0.07 <= theta <= 0.19, 7 <= KC <= 34); asymptotes to the current value 1.3 as KC -> infinity. Examples: KC = 10 -> S/D = 0.145; KC = 24 -> S/D = 0.54.

**E4. Combined waves + current (Sumer & Fredsoe 2001, J. Waterway Port Coastal Ocean Eng 127(5):403-411):**
- S/D = (S_c/D) * [1 - exp(-A * (KC - B))], for KC >= B, with S_c/D = 1.3 (current-only value),
- A = 0.03 + 0.75 * Ucw^2.6, B = 6 * exp(-4.7 * Ucw), Ucw = Uc / (Uc + Um),
- Uc = current velocity at z = D/2 above bed. Valid live-bed, KC up to ~30 in the underlying tests. Even a weak current on waves increases scour sharply; for Ucw >= ~0.7 the result is current-dominated (S/D -> 1.3). Example: KC = 6 waves alone give S/D ~ 0.08, but waves against current can still hold S/D ~ 0.2-0.5 depending on Ucw.

**E5. Time scale of scour development (Sumer & Fredsoe 2002):**
- Scour evolves as S(t) = S_eq * (1 - exp(-t/T)) (exponential approach; T = time scale).
- Nondimensional: T = T* * D^2 / sqrt(g * (s-1) * d50^3).
- Steady current: T* = (1/2000) * (delta/D) * theta^(-2.2), delta = boundary-layer/flow depth at the pile.
- Waves: T* = 1e-6 * (KC/theta)^3, valid 0.07 <= theta <= 0.19, 7 <= KC <= 34.
- Alternative (Whitehouse 1998, current): T* = 0.014 * theta^(-1.29).
- The strong negative theta exponent means live-bed scour at storm/tidal-peak conditions is fast (minutes-hours at lab/model scale, hours-days at full scale), clear-water scour is slow (days-weeks).

**E6. Shields parameter and regime test:**
- theta = Uf^2 / (g * (s-1) * d50); theta_cr ~ 0.05 (flat bed; use Soulsby's D*-dependent curve for the 0.1-10 mm range: theta_cr = 0.055 at d50 = 0.2 mm, ~0.03 at 1 mm, ~0.05 at 10 mm).
- Clear-water: theta < theta_cr; live-bed: theta > theta_cr. Local scour still occurs in clear-water because amplification alpha = 4-11 makes alpha*theta exceed theta_cr locally once theta > ~theta_cr/alpha.

**E7. Edge scour at scour protections (Petersen, Sumer, Fredsoe, Raaijmakers & Schouten 2015, Coastal Eng 106:42-72; Petersen, Sumer & Fredsoe, ICSE 2023):**
- Steady current (governing case): a pair of symmetric counter-rotating streamwise vortices forms in the near-bed wake of pile + protection and digs a significant scour hole immediately *downstream* of the berm; equilibrium edge-scour depth and length scale with pile diameter Dp and the berm aspect ratio Ar = hb/Wb (berm height / berm width).
- Waves: edge scour appears on the offshore and onshore sides, driven by wave-boundary-layer streaming and the roughness jump at the stone/sand junction; S/Dp = f(KC, Ar, Wb/Dp, theta), increasing with KC (tests to KC = 17, hb/Dp = 0.2-0.5) and approaching the current-alone values; sand is partly *deposited inside the porous stone layer*. Weak theta dependence in live-bed.
- Field (Egmond aan Zee OWF): edge scour developed over years, equilibrium after ~7-8 years.

## Practical guidance for our simulator

**Which formulas to implement as validation targets (milestone 1, steady current):**
1. **Rigid-bed shear amplification first.** Before enabling morphology, run the voxelized cylinder on a fixed bed and check alpha = tau/tau_0: expect alpha ~ 10 (7-11) at the 45-deg flanks and alpha ~ 4-5 under the HSV upstream. A semi-Lagrangian stable-fluids solver is numerically diffusive and will tend to smear the HSV — Roulund/Baykal needed ~100+ cells around the pile perimeter and RANS closure to get alpha within ~30%. With 2.5-5 cm voxels, a D = 1 m pile gives 60-125 cells per perimeter: adequate at the fine end; expect HSV (and hence upstream scour) under-prediction at the coarse end. Report alpha maps as a standard diagnostic.
2. **Equilibrium depth.** Live-bed runs should converge to S/D = 1.3 +/- 0.7; clear-water runs should follow the E2 ramp (zero below U/U_cr = 0.5, max near U/U_cr = 1, capped near 2.4 D). This brackets acceptable simulator output.
3. **Time scale.** Use E5 to size (and sanity-check) morphodynamic acceleration. Example at our scale: D = 1 m, delta = 5 m, U = 1.5 m/s over d50 = 0.5 mm sand (Uf ~ 0.06 m/s -> theta ~ 0.46, strongly live-bed): T* = (1/2000)*5*0.46^-2.2 = 0.014, T = 0.014 * 1 / sqrt(9.81*1.585*1.25e-10) ~ 310 s — equilibrium in ~15-25 min of physical time: directly simulable without acceleration. Clear-water gravel case: d50 = 10 mm, U = 1.5 m/s -> theta ~ 0.023 < theta_cr: ambient bed immobile, but alpha ~ 4-10 at the pile mobilizes it locally — the solver must therefore apply the *local, amplified* tau in the bedload/Exner step, never the ambient theta.
4. **Regime map for our parameter box** (h = 5 m, Soulsby threshold-current estimates): U_cr ~ 0.39 m/s (d50 = 0.2 mm), ~0.51 m/s (1 mm), ~1.57 m/s (10 mm). So tidal currents 0.5-2.5 m/s are live-bed for sand, and clear-water-to-marginal for the 10 mm gravel end — both regimes must work.
5. **Physics the solver must retain:** suspended load contributes ~50% of upstream equilibrium depth at theta ~ 0.13 (Baykal 2015) and more at our higher theta — do not tune scour with bedload alone. Downstream deposition/scour balance depends on wake unsteadiness; Smagorinsky LES should shed vortices if the grid resolves D by >~20 cells. Include a sand-slide/avalanche limiter at the angle of repose (32 deg) — equilibrium upstream slope equals it.
6. **For COBOD shapes (the actual product):** the relevant empirical anchor is E7, not bare-pile scour. Expect (a) a downstream edge-scour trench from the counter-rotating wake vortex pair, scaling with structure height/width ratio — flat, wide skirts reduce it; (b) deposition *inside* porous/rough geometry via streaming and shelter effects — exactly the self-ballasting behavior we want to maximize; validate qualitatively against Petersen et al. (2015) zone maps. Waves (KC-driven) can wait for the later milestone, but keep KC = Um*Tw/D and formulas E3/E4 in the spec now so the API carries Um, Tw from day one.

## Sources

- Sumer, B.M., Fredsoe, J. & Christiansen, N. (1992). Scour around vertical pile in waves. J. Waterway, Port, Coastal, Ocean Eng. 118(1):15-31. (S/D = 1.3[1-exp(-0.03(KC-6))]; mean 1.3, sigma 0.7 for current.)
- Sumer, B.M. & Fredsoe, J. (2001). Scour around pile in combined waves and current. J. Waterway, Port, Coastal, Ocean Eng. 127(5):403-411. (A, B, Ucw formula.)
- Sumer, B.M. & Fredsoe, J. (2002). The Mechanics of Scour in the Marine Environment. World Scientific. https://www.researchgate.net/publication/270892821_The_Mechanics_of_Scour_in_the_Marine_Environment
- Roulund, A., Sumer, B.M., Fredsoe, J. & Michelsen, J. (2005). Numerical and experimental investigation of flow and scour around a circular pile. J. Fluid Mech. 534:351-401. https://www.researchgate.net/publication/231897347
- Baykal, C., Sumer, B.M., Fuhrman, D.R., Jacobsen, N.G. & Fredsoe, J. (2015). Numerical investigation of flow and scour around a vertical circular cylinder. Phil. Trans. R. Soc. A 373:20140104. https://pmc.ncbi.nlm.nih.gov/articles/PMC4275922/
- Whitehouse, R. (1998). Scour at Marine Structures. Thomas Telford. (T* = 0.014 theta^-1.29.)
- Melville, B.W. & Coleman, S.E. (2000). Bridge Scour. Water Resources Publications. (K-factor method, 2.4 D envelope.)
- Petersen, T.U., Sumer, B.M., Fredsoe, J., Raaijmakers, T.C. & Schouten, J.-J. (2015). Edge scour at scour protections around piles in the marine environment. Coastal Eng. 106:42-72. https://www.sciencedirect.com/science/article/abs/pii/S0378383915001374
- Petersen, T.U., Sumer, B.M. & Fredsoe, J. (2023). Edge scour at scour protections around monopiles in waves. ICSE-11. https://www.issmge.org/uploads/publications/108/137/Edge_scour_at_scour_protections_around_monopiles_in_waves.pdf
- Myrhaug, D. & Ong, M.C. Time scales for scour below pipelines and around vertical piles (HR Wallingford eprint; recites Sumer & Fredsoe T* constants and validity ranges). https://eprints.hrwallingford.com/1105/1/PA_6_2-Myrhaug-D.pdf
- Mostafa, Y.E. & Agamy, A.F. (2011). Scour around single pile and pile groups subjected to waves and currents. IJEST 3(11). (Recites Sumer formulas with constants.) https://www.idc-online.com/technical_references/pdfs/chemical_engineering/SCOUR%20AROUND%20SINGLE%20PILE%20AND%20PILE%20GROUPS%20SUBJECTED%20TO%20WAVES%20AND%20CURRENTS.pdf
- Hartvig, P.A. et al. (Aalborg Univ. project report). Scour in a marine environment characterized by currents and waves. (Independent recitation of T* = (1/2000)(delta/D)theta^-2.2 and tidal-current scour results.) https://projekter.aau.dk/projekter/files/14765331/Scour_in_a_marine_env._cha._by_currents_and_waves_combined.pdf

## Verification

Adversarial fact-check (2026-07-07) of the eight key equations against sources located independently of this note's own citation list. Verification sources actually inspected (full text extracted and read, not abstracts, unless noted):

- [V1] Mostafa & Agamy (2011), "Scour around single pile and pile groups subjected to waves and currents", IJEST 3(11) — full-text PDF, equations (1), (5) and worked examples. https://www.idc-online.com/technical_references/pdfs/chemical_engineering/SCOUR%20AROUND%20SINGLE%20PILE%20AND%20PILE%20GROUPS%20SUBJECTED%20TO%20WAVES%20AND%20CURRENTS.pdf
- [V2] Rudolph et al., "Scour around offshore structures — analysis of field measurements" (ICSE/BAW TC213 proceedings) — full-text PDF, eqs. (3)-(4) and the N7 monopile design worked example. https://izw.baw.de/publikationen/tc213/0/CAS_21.pdf
- [V3] Hartvig et al. (Aalborg) project report — full-text PDF, eqs. [7]-[9], [17], [18], [26]-[34] (independent recitation chain: Offshore Center Denmark survey 2006).
- [V4] Myrhaug & Ong, "Random wave-induced time scales for scour below pipelines and around vertical piles" — full-text PDF, eqs. (3)-(8) with (p, r, s) coefficient table and test-range statement.
- [V5] Melville, B.W. (2008), "The Physics of Local Scour at Bridge Piers" (ICSE keynote) — full-text PDF. https://izw.baw.de/publikationen/tc213/0/k_2.pdf
- [V6] Tom, Draper & Yao (2021), "Estimating seabed shear stress amplification around circular cylinders" (ISFOG) — full-text PDF reviewing Hjorth 1975, Whitehouse 1998, Roulund 2005, Tavouktsoglou 2015 amplification data. https://www.issmge.org/uploads/publications/108/109/Estimating_seabed_shear_stress_amplification_around_circular_cylinders_-_an_observational_method_based_on_laboratory_experiments.pdf
- [V7] Baykal, Sumer, Fuhrman, Jacobsen & Fredsoe (2017), "Numerical simulation of scour and backfilling processes around a circular pile in waves", Coastal Eng. 122 — accepted-manuscript PDF (DTU Orbit). https://backend.orbit.dtu.dk/ws/files/129005572/manuscript_baykal_etal_CENG_submitted_r2_afterProof.pdf
- [V8] Larsen & Fuhrman (2023), "Re-parameterization of equilibrium scour depths and time scales for monopiles", Coastal Eng. 185:104356 — abstract/metadata only (full text paywalled). https://orbit.dtu.dk/en/publications/re-parameterization-of-equilibrium-scour-depths-and-time-scales-f/
- [V9] Escarameia & May (1999), HR Wallingford Report SR521, "Scour around structures in tidal flows" — full-text PDF. https://eprints.hrwallingford.com/436/2/SR521-Scour-structures-tidal-flows-HRWallingford.pdf
- [V10] Sheppard, Melville & Demir (2014), "Evaluation of Existing Equations for Local Scour at Bridge Piers", J. Hydraul. Eng. 140(1) — tabulation of the Breusers et al. (1977) equation (accessed via search excerpts only; paywalled).
- [V11] Baykal et al. (2015), Phil. Trans. R. Soc. A 373:20140104 — open-access full text at PMC. https://pmc.ncbi.nlm.nih.gov/articles/PMC4275922/

**E1 — Equilibrium scour, steady current, live-bed: CONFIRMED.** [V1] recites verbatim: "Sumer et al. (1992) stated that the mean value and standard deviation for the normalized equilibrium scour depth (S/D) for a vertical cylindrical pile in steady current are 1.3 and 0.7, respectively." [V3] eq. [17]: S/D = 1.3, sigma = 0.7, "valid for live bed". [V2] applies exactly the design logic stated here in a field case: "S_c/D = 1.3 + sigma = 2.0 and S_c/D = 1.3 + 2 sigma = 2.7" (and finds even mean+1sigma UNDER-predicted the observed 6.3 m scour at the N7 monopile after 5 years, S_max,obs/D = 1.05 on D = 6 m — treat 2.0 as an envelope for lab-like conditions, not an upper bound for multi-year field exposure). The 2.4·D cap: [V5], in Melville's own words, "the maximum depth of scour at a relatively narrow circular bridge pier is about 2.4 times the pier diameter, irrespective of the flow depth" (narrow-pier class D/h < 0.7, i.e. h/D > ~1.4 — consistent with the slenderness clause given here).

**E2 — Clear-water flow-intensity ramp (Breusers, Nicollet & Shen 1977): CONFIRMED, with one documented caveat on the 2.0 coefficient.** The functional form S/D = C·tanh(h/D)·f(U/U_cr) with piecewise f (0 below U/U_cr = 0.5; 2·U/U_cr − 1 between 0.5 and 1; 1 above 1) is confirmed: [V3] eqs. [7]-[9] recite s/D = f·k·tanh(h/D) with exactly that piecewise f, and [V9] independently confirms the onset threshold: "local scour begins when the undisturbed upstream flow velocity is equal to about half the value corresponding to the threshold of movement" (Breusers et al. 1977 analysis). Coefficient caveat: C = 2.0 is Breusers et al.'s *design* recommendation (as tabulated in [V10]: d_se/b = f1(U/Uc)·[2.0·tanh(h/b)]·Ks·Ktheta, with Ks = Ktheta = 1 for a circular pile); best-estimate recitations use C = 1.5 ([V2] eq. (3): S_max/D = 1.5·tanh(h/D)) or C = 1.3 ([V3]). Keep C = 2.0 as the conservative validation envelope, but expect simulator equilibria nearer 1.3-1.5·tanh(h/D)·f.

**E3 — Wave scour vs KC: CONFIRMED.** [V1] eq. (5) recites S/D = 1.3·[1 − exp(−0.03·(KC − 6))] for KC ≥ 6 verbatim, states "if KC ≤ 6, scour does not occur" (vortex-shedding cutoff KC = 6 for a circular pile), and its worked example reproduces the arithmetic here (KC = 10 -> S/D = 0.147). [V3] eq. [18] identical. KC = Um·Tw/D and KC = 2*pi*A/D: [V4] eqs. (4), (8). Test ranges 0.07 <= theta <= 0.19 and 7 <= KC <= 34: [V4] verbatim ("The tests upon which this formula is based were performed for 0.07 <= theta <= 0.19 and 7 <= KC <= 34").

**E4 — Combined waves + current: CONFIRMED (validity range sharpened).** A = 0.03 + 0.75·Ucw^2.6 and B = 6·exp(−4.7·Ucw) confirmed digit-for-digit in [V2] eq. (4) and [V3] eqs. [27]-[28]. Ucw = Uc/(Uc + Um) confirmed in both ([V3] distinguishes Uw for regular vs Um for irregular waves). Uc measured at z = D/2 above the bed: confirmed by multiple independent recitations of Sumer & Fredsoe (2001) (e.g. China Ocean Eng. 2023 full-scale study; ScienceDirect-indexed wave-current pile-scour papers). Validity: the underlying tests span approximately 4 <= KC <= 26 ([V3] recites "for KC >= 4"; other recitations state 4 < KC < 26) — so "KC up to ~30" as written is slightly generous; use 4-26 as the tested envelope. Current-dominated for Ucw >= ~0.7 checks out numerically (Ucw = 0.7, KC = 10: A = 0.327, B = 0.22, S/D = 1.25 ~ 1.3).

**E5 — Time development S(t) = S_eq·(1 − exp(−t/T)): CONFIRMED.** [V3] eq. [34] verbatim; also the same exponential form in [V2] and [V4] context.

**E6 — Time scale, steady current: CONFIRMED.** T* = (1/2000)·(delta/D)·theta^(−2.2): [V3] recites it twice (section 2.2.1 and eq. [33]) and notes their own tidal-current experiments re-fit gives T* = 0.0022·(delta/D)·theta^(−2.43) at delta/D = 3, "very similar" — i.e. the 1/2000 and −2.2 constants are the canonical Sumer & Fredsoe (2002) values. Nondimensionalization T = T*·D^2/sqrt(g·(s−1)·d50^3): [V3] eq. [31], [V4] (same T* definition for piles and pipelines). Whitehouse (1998) alternative T* = 0.014·theta^(−1.29): [V3] eq. [32] verbatim, attributed to Whitehouse (1998). s = 2650/1025 = 2.585 ~ 2.59: arithmetic checks. CAUTION for the simulator: [V8] (Larsen & Fuhrman 2023) argue on sediment-transport-scaling grounds that T* should go as theta^(−3/2), and that steeper lab-fit exponents like −2.2 "may become unreliable when extrapolated to field scales" — treat the theta exponent as uncertain between −1.3 and −2.2 when validating full-scale time scales, and do not tune the morphodynamic model to the −2.2 curve alone.

**E7 — Time scale, waves: CONFIRMED.** [V4] eq. (3): T* = p·KC^r·theta^(−s) with (p, r, s) = (1e−6, 3, 3) for piles — exactly T* = 1e−6·(KC/theta)^3 — with the same T-to-T* nondimensionalization, "valid for live-bed scour (theta > theta_cr)", tests at 0.07 <= theta <= 0.19 and 7 <= KC <= 34 (verbatim).

**E8 — Shields parameter and bed-shear amplification: CONFIRMED, with Reynolds-number caveat on the peak alpha.** theta = Uf^2/(g·(s−1)·d50) with theta_cr ~ 0.05 ("critical value of motion at a flat bed, i.e. theta_cr ~ 0.05"): [V4] eq. (5) verbatim. Amplification geometry and magnitudes ([V6], reviewing the primary experiments): Hjorth (1975) — maximum at ~45 deg from the flow-attack direction, instantaneous localized alpha ~ 7 (a zone of alpha ~ 3 extends ~0.5·D from the surface); recitations of Hjorth+Baker give the range 7-11 ([V7]: "amplified by a factor of 7-11"); Roulund et al. (2005) — maximum at 45-70 deg, peak alpha ~ 10 at ReD = 2e2, DECREASING with increasing ReD (~6.2 at higher ReD per Tavouktsoglou et al. 2015); Whitehouse (1998) review: typical 4-6 for circular/rectangular structures. Under the horseshoe vortex upstream: alpha ~ 4-5 (Roulund front-of-pile hot-film data as reproduced in [V11] Fig. 2-3, model within ~30%). Waves: side-edge amplification O(4) at KC ~ 10-20 (Sumer et al. 1997 data reproduced in [V7]; lee-wake alpha up to ~2 at KC = 10, vs ~6 max lee-side in steady current). Implication for the rigid-bed validation target in item 1 above: at field-scale ReD expect peak alpha ~ 5-8 rather than 10-11; the 10-11 ceiling stems from low-ReD/instantaneous measurements. The local-mobilization criterion alpha·theta > theta_cr is the standard clear-water scour argument (consistent with [V6]'s alpha_max = tau_local,max/tau_ff framework).

Cross-cutting result: no numeric constant in E1-E8 was found to be garbled; two validity statements were tightened (E4: tested KC range 4-26, not ~30; E2: 2.0 = design coefficient, 1.3-1.5 = best estimate) and two extrapolation cautions added (E6 theta-exponent at field scale; E8 alpha vs ReD).
