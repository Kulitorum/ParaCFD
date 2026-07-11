# Extracting Bed Shear Stress from a CFD Velocity Field

## Overview

Every sediment module downstream of the fluid solver (Shields criterion, bedload, pickup, Exner) is driven by the bed shear stress vector tau_b [N/m2], or equivalently the friction velocity u* = sqrt(tau_b/rho) [m/s]. On a 2.5-5 cm voxel grid we cannot resolve the viscous sublayer (~1 mm) or the grain-scale roughness layer, so tau_b must come from a **wall function**: sample the tangential velocity at a known height above the bed and invert the logarithmic law of the wall. This note gives the exact log-law forms for smooth/transitional/rough regimes, the roughness closure (Nikuradse ks, roughness length z0), the sampling/iteration recipe, the simpler quadratic drag alternative, and the known stair-step artifacts of voxelized beds with mitigations. Constants: kappa = 0.40 (von Karman; OpenFOAM/Fluent use 0.41 — pick 0.40 and stay consistent), rho = 1025 kg/m3, nu = 1.05e-6 m2/s.

## Key equations & results

### 1. Log law over a rough bed (fully rough form)

u(z)/u* = (1/kappa) * ln((z - z_b)/ks) + 8.5        [Nikuradse 1933; Ting & Ganti 2022, Eq. 1]

- u(z) = time-averaged tangential velocity at height z above the theoretical bed level z_b [m/s]
- kappa = 0.40 (dimensionless), B = 8.5 for close-packed uniform sand (fully rough)
- ks = equivalent Nikuradse sand-grain roughness [m]
- Equivalent compact form: **u(z)/u* = (1/kappa) * ln(z/z0)** with **z0 = ks/30**, because exp(-kappa*B) = exp(-0.40*8.5) = e^-3.4 = 1/29.96 ≈ 1/30. z0 is the height where the extrapolated log profile reaches zero velocity.
- Valid fitting range from measured-profile practice: 0.2*ks <= (z - z_b) <= (0.2-0.3)*h, h = water depth (Ting & Ganti 2022, citing Sumer & Fuhrman 2020). For a single-point wall function be stricter: z above the roughness crests, z >= ks (Fluent guidance: first cell centroid above ks) and z <= 0.1-0.2 h.
- Displacement height: theoretical bed z_b sits a fraction of a grain below the grain tops (Kamphuis: 0.7*d90 above the glued-grain base); at our resolution set z_b = the Exner bed surface and ignore the sub-grain offset.

### 2. Roughness regimes (grain Reynolds number ks+)

ks+ = u* * ks / nu:
- **Hydraulically smooth**: ks+ < 5 — roughness buried in viscous sublayer
- **Transitional**: 5 <= ks+ <= 70
- **Fully rough**: ks+ > 70 — tau_b independent of viscosity
(Classical Nikuradse/Schlichting limits; CFD wall functions in OpenFOAM/Fluent use the Cebeci–Bradshaw variant with limits 2.25 and 90, see Eq. 7.)

### 3. Roughness closure for sand/gravel beds

Flat sand bed, grain roughness only:
**ks = 2.5 * d50**, i.e. **z0 = d50/12**  [Soulsby 1997 as used in Soulsby & Clarke 2005, Appendix A: "zo = d50/12"; FLOW-3D scour model uses ks = cs*d50, cs = 2.5 recommended]. Van Rijn prefers kg ≈ 3*d90 (Houwman & van Rijn 1999) — for well-sorted sand these nearly coincide (d90 ≈ 2*d50 gives 6*d50 vs 2.5*d50; sources genuinely disagree within a factor ~2, and tau_b is only logarithmically sensitive).

Total roughness is the **sum of components** (Houwman & van Rijn 1999):
kp = kg + kf + kt (grain + form + sediment-transport roughness), z0 = kp/30.

Rippled bed form roughness:
**kf = a * eta^2 / lambda**, eta = ripple height [m], lambda = ripple length [m],
a = 8 (Nielsen 1992), a = 20 (van Rijn 1993), a = 27.7 (Grant & Madsen 1982) [Houwman & van Rijn 1999].
Example: eta = 2 cm, lambda = 14 cm, a = 8 gives kf = 8*0.0004/0.14 = 0.023 m, z0 = 0.76 mm — ~20x the grain z0. Ripples are sub-grid at 2.5-5 cm cells, so this enters only through z0. Van Rijn recommends a minimum total kp = 0.01 m for field sheet-flow conditions.

### 4. Transitional smooth/rough roughness length (all regimes, one formula)

**z0 = ks/30 * [1 - exp(-u* * ks / (27*nu))] + nu / (9 * u*)**   [Christoffersen & Jonsson 1985; used in REEF3D-type scour models and Fuhrman-group RANS models; Soulsby 1997 Eq. 23-24 gives the simpler sum z0 = ks/30 + nu/(9*u*)]

- Limits: ks+ -> infinity gives z0 = ks/30 (rough); ks+ -> 0 gives z0 = nu/(9*u*) (smooth, equivalent to u+ = (1/kappa) ln(9*z+) i.e. B_smooth ≈ 5.5).
- Because z0 now depends on u*, the wall function becomes implicit (see Eq. 6).

### 5. Wall-function extraction of u* (THE core formula for our code)

Given tangential velocity magnitude U_p at height z_p above the local bed:
**u* = kappa * U_p / ln(z_p / z0)**,   **tau_b = rho * u*^2**, direction = local tangential unit vector u_t/|u_t|.

- Fully rough with known z0 = ks/30: explicit, no iteration.
- Sampling height validity: z_p >= max(ks, 30*nu/u*) and z_p <= 0.1-0.2*h. For our grid, z_p = first cell center = Delta z/2 = 1.25-2.5 cm; y+ = u* z_p/nu ≈ 300-2000 — above the classic RANS 30 < y+ < 300 comfort band but standard in coastal/scour models where the log layer fills the lower 10-20% of the depth (5 m depth -> log layer to ~0.5-1 m).
- **Gravel warning**: d50 = 10 mm gives ks = 25 mm > z_p = 12.5 mm. The first cell center is INSIDE the roughness layer; sample the 2nd or 3rd cell up so z_p >= ~2*ks (>= 5 cm).

### 6. Iteration when z0 depends on u* (transitional/smooth)

Fixed-point loop (converges in 2-5 iterations to 1e-4 relative):
1. u*(0) = kappa*U_p / ln(30*z_p/ks)  (fully rough guess)
2. z0(n) = ks/30*[1-exp(-u*(n)*ks/(27*nu))] + nu/(9*u*(n))
3. u*(n+1) = kappa*U_p / ln(z_p/z0(n)); repeat.
Guard: clamp ln argument to >= e (i.e. z_p/z0 >= 2.718) so u* <= kappa*U_p.

### 7. CFD-code rough wall function for reference (OpenFOAM nutkRoughWallFunction, Cebeci–Bradshaw)

u+ = (1/kappa)*ln(E' * y+), kappa = 0.41, smooth E = 9.8 (Fluent: 9.793), E' = E/fn(ks+):
- ks+ < 2.25: fn = 1 (smooth)
- 2.25 <= ks+ < 90: fn = [(ks+ - 2.25)/87.75 + Cs*ks+]^(sin(0.4258*(ln ks+ - 0.811)))
- ks+ >= 90: fn = 1 + Cs*ks+, roughness constant Cs = 0.5 default
[OpenFOAM source, nutkRoughWallFunctionFvPatchScalarField.C]. With Cs = 0.5 the fully rough limit gives u+ ≈ (1/0.41) ln(19.6*z/ks), i.e. B ≈ 7.3 — slightly lower than Nikuradse's 8.5; another example of source disagreement at the 10-15% level in tau_b. OpenFOAM also rate-limits nut updates to [0.5x, 2x] per iteration "to avoid oscillations" — worth copying.

### 8. Quadratic friction law (simpler alternative / diagnostic)

**tau_b = rho * Cd * U_ref^2** with, for U_ref sampled at height z_ref in the log layer:
**Cd(z_ref) = [kappa / ln(z_ref/z0)]^2** (identical to Eq. 5, just precomputed).
For depth-averaged velocity U_bar in depth h:
**Cd = [0.40 / (ln(h/z0) - 1)]^2**  [Soulsby 1997 Eq. 37; Soulsby & Clarke 2005 Eq. 5/A6]
Smooth-turbulent depth-averaged alternative: CDs = 0.0001615 * exp(6 * Rec^-0.08), Rec = U_bar*h/nu, laminar below Rec = 2000 (tau = 3*rho*nu*U_bar/h) [Soulsby & Clarke 2005, Eqs. A1, A4, A11]. Regime selection: compute rough and smooth tau, take the larger.

## Practical guidance for our simulator

**Recommended scheme.** Per bed-surface column (or bed-adjacent face), each time step:
1. Reconstruct local bed normal n from the Exner height field / smoothed voxel bed (NOT from raw voxel faces).
2. Probe the velocity at a **fixed physical distance along n**: z_p = max(1.5*Delta z, 2*ks) above the bed surface, trilinear-interpolated. Fixed-height probing is the single most effective anti-stair-step measure.
3. Tangential component: u_t = u - (u.n)n, U_p = |u_t|.
4. u* from Eq. 5 with transitional z0 (Eq. 4, loop of Eq. 6). Use ks = 2.5*d50 (+ form roughness kf if/when a ripple parameterization is on).
5. tau_b = rho*u*^2 * u_t/U_p. Feed the **grain-stress fraction** to the sediment module: when z0 includes form roughness, recompute a skin-friction u*' with z0_grain = d50/12 for Shields/bedload (form drag does not move grains); with grain-only z0 they are identical.

**Concrete numbers** (z_p = 1.25 cm, U_p = 1 m/s, flat bed, grain roughness only):

| d50 | ks = 2.5 d50 | z0 = d50/12 | Cd(z_p) | u* (m/s) | tau_b (Pa) | ks+ (regime) |
|-----|------|------|---------|------|------|------|
| 0.1 mm | 0.25 mm | 8.3e-6 m | 3.0e-3 | 0.055 | 3.1 | ~13 transitional |
| 0.5 mm | 1.25 mm | 4.2e-5 m | 4.9e-3 | 0.070 | 5.0 | ~83 rough |
| 1 mm | 2.5 mm | 8.3e-5 m | 6.4e-3 | 0.080 | 6.6 | ~190 rough |
| 10 mm | 25 mm | 8.3e-4 m | invalid at 1.25 cm — probe at z_p >= 5 cm; Cd(5 cm) = 9.5e-3 | 0.098 | 9.8 | ~2300 rough |

So for d50 >= ~0.5 mm at our current speeds the bed is fully rough and Eq. 5 is explicit; below ~0.3 mm and/or at slack tide (u* < 0.03 m/s) the transitional correction matters (10-20% in tau_b) — implement Eq. 4+6 once, it covers everything.

**LES note.** Applying Eq. 5 to the instantaneous filtered velocity is the standard equilibrium wall model (Schumann-type). It yields an instantaneous tau_b vector whose direction fluctuates — good, this drives realistic scour. Time-average or exponentially smooth tau_b over ~1-5 s before Exner if the bed update goes noisy.

**Stair-step / voxel artifacts and mitigations** (Ji et al. 2021; Zhu et al. 2025):
- *Cause*: on Cartesian stair-step boundaries the wall distance of the first fluid cell is irregular (varies with surface position/orientation vs the grid), and skin friction is extremely sensitive to it -> spurious tau_b oscillations with grid-pitch wavelength; on a slope, a periodic error pattern of wavelength ~Delta z/tan(theta). Untreated, this checkerboards the Exner update and prints the grid into the bed.
- *Mitigations used in the literature*:
  1. Sample at fixed normal distance with trilinear/inverse-distance-weighted interpolation from surrounding fluid cells (not the single stair-step-adjacent cell).
  2. **Interpolate u\* rather than velocity** when populating ghost/forcing cells — u* varies far less near the wall than u (Ji et al. 2021, JCP 439:110240).
  3. Evaluate near-wall gradients by weighted least squares over the fluid-cell neighborhood; correct the wall-normal gradient of tangential velocity with the log law instead of one-sided differences (Ji et al. 2021).
  4. Two-layer wall functions tuned for moving sediment beds give smooth tau_b at >80% coarser grids (Zhu et al. 2025).
  5. Post-smooth tau_b tangentially along the bed (3x3 filter) before Exner, and pair with a sand-slide/avalanche limiter — standard in morphodynamic codes.
- Keep the immersed/voxel obstacle boundaries (the printed concrete shapes) on the same treatment: compute their local wall shear with the same fixed-normal-distance log-law probe, with ks of concrete ≈ 1-3 mm (form-cast) unless we decide printed-layer ridges dominate (then ks ≈ layer height).

**Pitfalls.** (a) Never finite-difference du/dz across the first cell to get tau_b — at y+ ~ 1000 that underestimates tau_b severely. (b) Do not let z_p fall below ks after bed accretion — recompute z_p against the moving bed every step. (c) kappa/B/Cs choices shift tau_b by 10-15% between codes; calibrate once against the Soulsby depth-averaged Cd (Eq. 8) in an unobstructed periodic channel and freeze the constants.

## Sources

- Ting, F.C.K. & Ganti, S.M. (2022). "Finding the Bed Shear Stress on a Rough Bed Using the Log Law." J. Waterway, Port, Coastal, Ocean Eng. 148(4). https://par.nsf.gov/servlets/purl/10335702
- Soulsby, R.L. & Clarke, S. (2005). "Bed Shear-stresses Under Combined Waves and Currents on Smooth and Rough Beds." HR Wallingford TR137 (implements Soulsby 1997, "Dynamics of Marine Sands", Eqs. 36-37, 62). https://eprints.hrwallingford.com/558/1/TR137.pdf
- Houwman, K.T. & van Rijn, L.C. (1999). "Flow resistance in the coastal zone." Coastal Engineering 38:261-273. https://www.leovanrijn-sediment.com/papers/P1-1999.pdf
- OpenFOAM nutkRoughWallFunction source (Cebeci-Bradshaw roughness function). https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/MomentumTransportModels/momentumTransportModels/derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkRoughWallFunction/nutkRoughWallFunctionFvPatchScalarField.C
- Ji, C. et al. (2021). "An improved immersed boundary method for turbulent flow simulations on Cartesian grids." J. Comput. Phys. https://www.sciencedirect.com/science/article/abs/pii/S0021999121001352
- Zhu, H. et al. (2025). "Improved immersed boundary/wall modeling method for RANS solver coupled with wall functions: application to Cartesian grid systems." Eng. Appl. Comput. Fluid Mech. https://www.tandfonline.com/doi/full/10.1080/19942060.2025.2486657
- FLOW-3D Sedimentation & Scour Model report (ks = cs*d50, cs = 2.5). https://flow3d.co.kr/wp-content/uploads/FSR-03-14_sedimentation-scour-model.pdf
- Christoffersen, J.B. & Jonsson, I.G. (1985). "Bed friction and dissipation in a combined current and wave motion." Ocean Engineering 12(5):387-423 (transitional z0 formula; as adopted in REEF3D-type scour models).

## Verification

Adversarial fact-check (2026-07-07). Each key equation was re-verified against primary sources fetched independently of this note: the full text of Soulsby & Clarke TR137 (PDF, HR Wallingford eprints), the full text of Houwman & van Rijn 1999 (PDF, leovanrijn-sediment.com), the full text of Ting & Ganti 2022 (PDF, NSF Public Access Repository), the OpenFOAM-dev source tree on GitHub, the USACE HEC-RAS 2D Sediment Transport Technical Reference, the Ansys Fluent 12.0 User's Guide (Sec. 7.3.14), and Flack & Schultz (2014), "Roughness effects on wall-bounded turbulent flows," Phys. Fluids 26, 101305 (USNA copy). All numeric constants, exponents, and worked-table values below were recomputed by hand.

1. **Rough-bed log law (Eq. 1) — CONFIRMED.** Ting & Ganti (2022) Eq. (1) full text confirms u/u* = (1/kappa) ln((y-y0)/ks) + B with B = 8.5 and kappa = 0.40 for fully rough flow, y-intercept = ln(ks) - 8.5*kappa, and the fitting range verbatim: "0.2ks <= y - y0 <= (0.2-0.3)h". Kamphuis y0 = 0.7*d90 confirmed in the same paper. Compact-form equivalence checked arithmetically: exp(0.40*8.5) = e^3.4 = 29.96 ~= 30, so z0 = ks/30. Fully rough threshold ks+ > 70: Flack & Schultz (2014), attributing the onset to Schlichting/Nikuradse. Single-point rule "z above ks": Fluent UG 7.3.14 ("distance from the wall to the centroid of the wall-adjacent cell should be greater than Ks").
2. **Wall-function u* extraction (Eq. 5) — CONFIRMED.** Exact algebraic inversion of the confirmed log law; tau_b = rho*u*^2 is the definition of friction velocity. Sampling constraints consistent with Fluent centroid-above-Ks guidance and with the Ting & Ganti upper fitting bound (0.2-0.3h; the stricter 0.1-0.2h used here is conservative and safe). Table values recomputed: all four rows of Cd, u*, tau_b, ks+ reproduce to the stated precision (e.g. d50 = 0.5 mm: ln(0.0125/4.17e-5) = 5.70, u* = 0.4/5.70 = 0.070 m/s, tau = 1025*0.070^2 = 5.0 Pa, ks+ = 0.070*1.25e-3/1.05e-6 = 83).
3. **Grain roughness closure ks = 2.5*d50, z0 = d50/12 (Eq. 3) — CONFIRMED.** Extracted verbatim from TR137 full text: "zo = d50/12" (App. A, roughness conversion step) and "for a flat, hydrodynamically rough bed of sediment zo = 2.5d/30" (main text). Van Rijn alternative kg ~= 3*D90 confirmed verbatim in Houwman & van Rijn (1999) Sec. 4. The four z0 values check: d50/12 = 8.3e-6, 4.17e-5, 8.3e-5, 8.3e-4 m.
4. **Transitional z0 (Eq. 4) — CONFIRMED.** USACE HEC-RAS 2D Sediment Transport Technical Reference ("Bottom Roughness") gives, citing Christoffersen & Jonsson (1985), exactly: z0 = (ks/30)[1 - exp(-u*ks/(27*nu))] + nu/(9u*). Independently, Roelvink (UN-IHE) lecture notes give the same form with regime bands Re* = 5 and 70. The smooth limit z0 = nu/(9u*) appears verbatim in TR137 Eq. (36) ("z0 = nu/(9u*e)"); the smooth-limit equivalence to u+ = (1/kappa) ln(9 y+) => B_smooth = ln(9)/0.4 = 5.49 checks against the classical smooth-wall value ~5.5. *Caveat*: the side-claim that Soulsby (1997) Eqs. 23-24 give the simpler sum z0 = ks/30 + nu/(9u*) could not be independently checked against the book (no accessible copy); both limiting forms are confirmed, only the exact equation numbers remain unverified.
5. **Ripple form roughness (Eq. 3b) — CONFIRMED.** Houwman & van Rijn (1999) full text: kf = a*eta^2/lambda ("kf = a h (h/l)") with "a in the range from 8 (Nielsen) to 27.7 (Grant and Madsen, 1982)", and explicitly "a = 8 (Nielsen, 1992), a = 20 (Van Rijn, 1993), a = 27.7 (Grant and Madsen, 1982; Li et al., 1996)". Total roughness as sum kp = kg + kf + kt confirmed verbatim (Sec. 4.1). Minimum roughness confirmed: "Van Rijn (1993) recommends a minimum physical roughness value equal to ks = 0.01 m for field" (sheet-flow) conditions. Minor attribution note: the paper credits a = 8 to both Nielsen 1981 and Nielsen 1992 in different passages. Worked example checks: 8*0.02^2/0.14 = 0.0229 m, /30 = 0.76 mm.
6. **Quadratic drag law (Eq. 8) — CONFIRMED.** TR137 full text, Eq. (A6): CDr = [0.40/(ln(h/z0) - 1)]^2; Eq. (39)/(A4): CDs = 0.0001615 exp[6 Rec^-0.08] with Rec = U*h/nu (notation list confirms Rec = Uh/nu); laminar branch: "If Rec <= 2000, then laminar flow" with tau_m = tau_max = 3*rho*nu*U/h (Eq. A9). Cd(z_ref) = [kappa/ln(z_ref/z0)]^2 is the direct log-law identity.
7. **Cebeci-Bradshaw / OpenFOAM roughness function (Eq. 7) — CONFIRMED, one attribution CORRECTED.** Fetched nutkRoughWallFunctionFvPatchScalarField.C (OpenFOAM-dev master): thresholds 2.25 and 90; transitional E' = E/[(KsPlus - 2.25)/87.75 + Cs*KsPlus]^sin(0.4258*(ln KsPlus - 0.811)); fully rough E' = E/(1 + Cs*KsPlus); nut rate-limited to [0.5x, 2x] per iteration — all exactly as stated. Defaults kappa = 0.41, E = 9.8 confirmed in the base class nutWallFunctionFvPatchScalarField.C (lookupOrDefault). **Correction**: in OpenFOAM, Cs (and Ks) are *required* dictionary entries with no code default — the header's usage example shows Cs = 0.5. Cs = 0.5 is the *Ansys Fluent* default ("chosen to reproduce Nikuradse's resistance data for tightly-packed uniform sand-grain roughness", Fluent UG 7.3.14, which also confirms E = 9.793 and the same Cebeci-Bradshaw regimes). The derived fully-rough equivalent B ~= ln(9.8/0.5)/0.41 = 7.26 ~= 7.3 checks.
8. **ks+ regime thresholds (Eq. 2) — CONFIRMED.** Flack & Schultz (2014): "Nikuradse (1933) observed that sand grain roughness was hydraulically smooth for ks+ < 5, fully rough for ks+ > 70, and transitional between those extremes"; Schlichting places fully-rough onset at ks+ = 70. Caveat also noted there: for non-sand-grain roughness types the transitional band varies widely (lower bound 1.4-15, upper 18-70), so the 5/70 limits are specific to uniform-sand ks — appropriate for our quartz-sand beds. CFD-code variant limits 2.25/90 confirmed directly from the OpenFOAM source (item 7).

Independent sources used: https://eprints.hrwallingford.com/558/1/TR137.pdf ; https://www.leovanrijn-sediment.com/papers/P1-1999.pdf ; https://par.nsf.gov/servlets/purl/10335702 ; https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-dev/master/src/MomentumTransportModels/momentumTransportModels/derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkRoughWallFunction/nutkRoughWallFunctionFvPatchScalarField.C (+ .H, + ../nutWallFunction/nutWallFunctionFvPatchScalarField.C) ; https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/latest/model-description/bedform-geometry-and-hydraulic-roughness/bottom-roughness ; https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node250.htm (Fluent 12.0 UG 7.3.14) ; https://www.usna.edu/NAOE/_files/documents/Faculty/schultz/Flack,%20Schultz%20-%20Roughness%20Effects%20on%20Wall,%202014.pdf ; https://ocw.un-ihe.org/pluginfile.php/117578/mod_resource/content/0/Coastal_Sediment_Transport_part1-3.pdf
