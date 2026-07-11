# Settling Velocity & Suspended-Load Transport

## Overview

This note specifies the sediment-settling and suspended-load physics for the COBOD scour simulator: single-grain settling velocity w_s (Stokes, Ferguson & Church 2004, Soulsby 1997), a tabulated w_s for our quartz-in-seawater range (rho_s = 2650 kg/m3, rho_w = 1025 kg/m3, so s = 2.585, s-1 = R = 1.585), hindered settling (Richardson-Zaki), the suspension criterion, the Rouse number/profile (used for validation), the 3D advection-diffusion equation with settling and turbulent Schmidt number, and the erosion (van Rijn pickup) / deposition (w_s * c_b) bed boundary condition that couples the concentration field to the Exner bed update. All formulas use SI units; concentrations c are VOLUMETRIC (m3/m3) unless noted. Convert mass concentration [kg/m3] = rho_s * c.

Key dimensionless group used throughout (van Rijn / Soulsby):

    D* = d * [ (s-1) * g / nu^2 ]^(1/3)      (dimensionless grain size)

For our conditions (g = 9.81, s-1 = 1.585, nu = 1.05e-6 m2/s): D* = 24160 * d, i.e. D* = 2.42 (d = 0.1 mm) ... 241.6 (d = 10 mm).

## Key equations & results

### 1. Stokes law (Ferguson & Church 2004, Eq. 1)

    w_s = R * g * d^2 / (C1 * nu),   C1 = 18 (theoretical, spheres)

Valid ONLY for particle Reynolds number Re = w_s*d/nu < ~1, i.e. quartz d < ~0.1 mm in water. At d = 0.2 mm it already overpredicts by ~35% (Re ≈ 6). Do not use alone for our 0.1-10 mm range.

### 2. Ferguson & Church (2004) universal formula (their Eq. 4) — primary source read

    w_s = R * g * d^2 / ( C1 * nu + sqrt(0.75 * C2 * R * g * d^3) )

- R = (rho_s - rho_w)/rho_w = submerged specific gravity (1.585 for quartz in seawater; 1.65 in fresh water)
- Constants (C1, C2): smooth spheres (18, 0.4); **typical natural sand, sieve diameter (18, 1.0)** — recommended; natural sand, nominal diameter (20, 1.1); very angular grains (24, 1.2).
- Asymptotes: Stokes law for small d, constant drag coefficient C_D = C2 for Re > 1e3. Valid for the entire size range including the transitional band 0.1 < d < 4 mm.
- Accuracy: rms error 4-8% vs Raudkivi data, ~16% vs Hallermeier's 115-point multi-density dataset, 6% vs their own Fraser River sand experiments (d = 0.068-4.36 mm).

### 3. Soulsby (1997) formula ("Dynamics of Marine Sands", Eq. SC 102; reproduced in HEC-RAS 2D sediment docs)

    w_s = (nu / d) * [ sqrt( 10.36^2 + 1.049 * D*^3 ) - 10.36 ]
        (10.36^2 = 107.33)

Calibrated on natural (irregular) marine sand grains; good for 0.063 <= d <= ~2 mm and remains sane to 10 mm. With hindered settling built in (Soulsby 1997):

    w_s(C) = (nu / d) * [ sqrt( 10.36^2 + 1.049 * (1-C)^4.7 * D*^3 ) - 10.36 ]

### 4. Settling velocity table — quartz in seawater (rho_s = 2650, rho_w = 1025, g = 9.81)

nu = 1.05e-6 m2/s (~15 C, per project spec):

| d (mm) | D* | Stokes (m/s) | F&C C1=18,C2=1.0 (m/s) | Soulsby (m/s) | Re = w*d/nu |
|--------|------|--------------|------------------------|----------------|-------------|
| 0.1 | 2.42 | 0.0082 (Re=0.8, marginal) | 0.0070 | **0.0073** | 0.7 |
| 0.2 | 4.83 | 0.0329 (invalid) | 0.0218 | **0.0245** | 4.7 |
| 0.5 | 12.1 | invalid | 0.0681 | **0.0711** | 34 |
| 1.0 | 24.2 | invalid | 0.1226 | **0.1173** | 112 |
| 2.0 | 48.3 | invalid | 0.1918 | **0.1753** | 334 |
| 5.0 | 120.8 | invalid | 0.317 | **0.283** | 1350 |
| 10.0 | 241.6 | invalid | 0.453 | **0.403** | 3840 |

At nu = 1.36e-6 (cold ~5 C seawater) the fine end drops: w_s(0.1 mm) = 0.0057, w_s(0.2 mm) = 0.0202 m/s; d >= 1 mm changes < 3%. F&C and Soulsby agree within ~10% across the whole range.

### 5. Hindered settling — Richardson & Zaki (1954)

    w_s,m = w_s,0 * (1 - C)^n

- C = local volumetric concentration, valid for C = 0 to ~0.3.
- n depends on particle Reynolds number: n = 4.65 (Re < 0.2), decreasing to n ≈ 2.4 (Re > 500). Van Rijn, Bisschop & van Rhee (2019, primary source read) use n = 4 to 5 for sand; HEC-RAS docs: n = 3.75-4.45 for 0.05-0.5 mm sand, ~4.0 typical. Soulsby's variant uses the (1-C)^4.7 factor inside the D*^3 term (Sec. 3).
- Negligible below C ≈ 0.001 (2.65 kg/m3): (1-0.001)^4.7 = 0.995. At C = 0.05, factor = 0.79; at C = 0.1, factor = 0.61.

### 6. Suspension criterion (when does sand go into suspension?)

van Rijn (1984, Part II, J. Hydraul. Eng. 110(11)):

    u*_cr,susp / w_s = 4 / D*      for 1 < D* <= 10
    u*_cr,susp / w_s = 0.4         for D* > 10

Bagnold's stricter criterion (full suspension): u*/w_s > 1. Van Rijn (1993/2012; verified in van Rijn's own "Simple general formulae" note) also gives a Shields-style curve:

    theta_cr,susp = 0.3/(1 + D*) + 0.1*[1 - exp(-0.05*D*)]
    u*_cr,susp = sqrt( theta_cr,susp * (s-1) * g * d50 )

For our sediments (nu = 1.05e-6): u*_cr,susp = 0.012 m/s (0.1 mm), 0.015 (0.2 mm), 0.023 (0.5 mm), 0.036 (1 mm), 0.055 (2 mm), 0.126 m/s (10 mm). With u* ≈ 0.05*U for a sandy seabed, our 0.5-2.5 m/s currents give u* ≈ 0.025-0.125 m/s: 0.1-0.2 mm sand is suspended at nearly all tidal stages; 1-2 mm sand suspends only above U ≈ 0.9-1.4 m/s; 10 mm gravel virtually never suspends (bedload only).

### 7. Rouse number and Rouse profile (equilibrium check case; van Rijn 1993/2012)

Steady balance c*w_s + eps_s * dc/dz = 0 with parabolic diffusivity eps_s = kappa*u**z*(1 - z/h) gives:

    c(z)/c_a = [ ((h - z)/z) * (a/(h - a)) ]^Z ,   Z = w_s / (beta * kappa * u*)

- kappa = 0.4 (von Karman), a = reference level, c_a = reference concentration, h = depth.
- beta = ratio of sediment to fluid mixing = inverse turbulent Schmidt number; van Rijn (1984-II): beta = 1 + 2*(w_s/u*)^2, valid 0.1 < w_s/u* < 1, clipped to 1 <= beta <= 1.5 (Delft3D convention).
- Interpretation (van Rijn): Z = 5 suspension confined below z = 0.1h; Z = 2 up to mid-depth; Z = 1 reaches surface; Z = 0.1 nearly uniform over depth. Use this analytic profile as a unit test of the 3D solver in equilibrium open-channel flow.

### 8. 3D advection-diffusion for suspended concentration (standard; cf. COHERENS ch.7, Delft3D)

    dc/dt + d(u_i c)/dx_i - d(w_s c)/dz = d/dx_i [ (nu_t/sigma_s + nu_mol) dc/dx_i ]

- Settling appears as extra downward advection velocity w_s added to the fluid vertical velocity (z positive up).
- sigma_s = turbulent Schmidt number. Measured values scatter 0.2-1.0 (Gualtieri et al.; APS Phys. Rev. Fluids 7, 014307 notes the controversy); **sigma_s = 0.7 is the common engineering default**; equivalently eps_s = beta*nu_t with beta from Sec. 7 (beta = 1/sigma_s in [1, 1.5]).
- In our Smagorinsky LES: nu_t = (C_s*Delta)^2*|S|, diffusivity = nu_t/sigma_s per cell. Free-surface BC: total vertical flux zero (w_s*c + eps_s*dc/dz = 0).

### 9. Near-bed reference concentration — van Rijn (1984, Part II)

    c_a = 0.015 * (d50 / a) * T^1.5 / D*^0.3      [volumetric, -]

    T = (tau'_b - tau_cr) / tau_cr   (transport-stage parameter, grain-related skin shear stress tau'_b)
    tau_cr = rho_w * (s-1) * g * d50 * theta_cr   (Shields; theta_cr ≈ 0.03-0.06, use Soulsby: theta_cr = 0.3/(1+1.2*D*) + 0.055*[1-exp(-0.02*D*)])

- Reference height a = max(0.5*bedform_height, k_s), with floor a >= 0.01*h. Cap c_a <= 0.05 (vol.) — above this the dilute-suspension assumptions fail.
- E from Sec. 10 and c_a are alternative bed forcings: Dirichlet (impose c = c_a at z = a) vs flux (impose pickup E). For a voxel LES solver the flux form is far more robust.

### 10. Van Rijn pickup (erosion) function — primary sources read

Original (van Rijn 1984, "Sediment pick-up functions", J. Hydraul. Eng. 110(10):1494-1502):

    E = 0.00033 * rho_s * [ (s-1) * g * d50 ]^0.5 * D*^0.3 * T^1.5      [kg/m2/s]

Calibrated for U = 0.5-1.5 m/s, d50 = 100-1500 um. High-velocity extension (van Rijn, Bisschop & van Rhee 2019, J. Hydraul. Eng. 145(1)), recalibrated to 6 m/s on 50-560 um sand:

    E = alpha * rho_s * [ (s-1)*g*d50 ]^0.5 * D*^0.3 * f_D * T^1.5
    alpha = 0.00033 (+-30%);  f_D = 1/theta'  for theta' > 1, else f_D = 1
    (theta' = grain-related Shields parameter; damping represents turbulence collapse + dilatancy)

Accuracy: within factor ~2; underpredicts ~50% for 125-150 um at very high velocity. Valid for clean sand, porosity 0.4 +- 0.03.

### 11. Deposition flux and the bed exchange boundary condition

    D_flux = w_s * c_b        [m/s * (-) = m3/m3 * m/s; multiply by rho_s for kg/m2/s]

c_b = near-bed concentration evaluated at the reference level / first cell above the bed (use hindered w_s(c_b) when c_b > 0.001). FUNWAVE/Cao (1999) optionally multiply by (1 - gamma*c)^2 — a refinement, not needed initially. The bed exchange is imposed as the bottom boundary flux of the advection-diffusion equation:

    net upward flux at bed:  F_bed = E/rho_s - w_s * c_b      [m/s, volumetric]

i.e. add E as a source and w_s*c_b as a sink in each bed-adjacent fluid cell (per horizontal area). Couple to Exner with bed porosity p ≈ 0.4:

    (1 - p) * dz_b/dt = w_s*c_b - E/rho_s - div(q_b)     [q_b = volumetric bedload flux, m2/s]

## Practical guidance for our simulator

1. **Use Soulsby (1997) as the single w_s formula** (one sqrt, no branches, validated for marine sand, matches F&C within 10%). Precompute w_s per grain class at startup from d50, nu, R = 1.585; store as a constant. Keep Ferguson & Church (C1 = 18, C2 = 1.0) in the test suite as a cross-check; both must reproduce the Sec. 4 table to < 1%.
2. **Concrete constants**: nu = 1.05e-6 m2/s (spec baseline; recompute if T = 5 C -> 1.36e-6, only fine sand is sensitive), kappa = 0.4, sigma_s = 0.7, p = 0.4, theta_cr from Soulsby's D* fit, d90 ≈ 2*d50, k_s = 2.5*d50 (grain) or 3*d90.
3. **Split by grain size**: d50 <= 0.5 mm -> full suspended + bedload; d50 >= 2 mm -> bedload/Exner only (w_s = 0.18-0.40 m/s, Rouse Z > 5 for u* < 0.14 m/s; suspension negligible). This lets the suspended-load kernel be skipped for gravel runs — big GPU saving.
4. **Hindered settling**: implement w_s,eff = w_s0*(1 - c)^4.7 clamped at c <= 0.35; it only activates near the bed during vigorous scour (c > 0.001 ≈ 2.65 kg/m3) and costs one powf.
5. **Bed BC**: use the flux form (Sec. 11) with van Rijn pickup INCLUDING the 2019 damping factor f_D = 1/theta' for theta' > 1 — our 2.5 m/s tidal case can reach theta' > 1 for 0.1-0.2 mm sand, where the undamped 1984 function overpredicts erosion. Evaluate tau'_b from the LES wall shear (grain roughness only), not total stress.
6. **Reference level vs voxel size**: with 2.5-5 cm cells and h ≈ 5 m, a = 0.01h = 5 cm ≈ 1 cell — the first cell center doubles naturally as the reference level; apply E and w_s*c_b there. Cap c_a-equivalent concentrations at 0.05 vol.
7. **Validation**: (i) drop test — particle-settling column must reproduce Sec. 4 to <1%; (ii) equilibrium channel — steady 1 m/s flow, d50 = 0.2 mm (Z ≈ w_s/(0.4*u*) ≈ 1.2 at u* = 0.05 m/s) must relax to the Rouse profile of Sec. 7; (iii) mass balance — Exner + suspended + bedload inventories must close to machine precision.
8. **Magnitude sanity**: at U = 1 m/s, h = 5 m, d50 = 0.25 mm, expect depth-integrated suspended transport ≈ 0.4-0.6 kg/s/m (van Rijn TR2004 / simplified formula q_s = 0.012*rho_s*U*d50*Me^2.4*D*^-0.6); total transport scales as U^3-U^5 — small velocity errors triple transport errors.

## Sources

- Ferguson, R.I. & Church, M. (2004), "A Simple Universal Equation for Grain Settling Velocity", J. Sedimentary Research 74(6), 933-937. https://geoweb.uwyo.edu/geol5330/FergusonChurch_GrainSettling_JSR04.pdf (primary, read in full)
- van Rijn, L.C., Bisschop, R. & van Rhee, C. (2019), "Modified Sediment Pick-up Function", J. Hydraulic Engineering 145(1). https://www.leovanrijn-sediment.com/papers/Pickup2019.pdf (primary, read in full; restates van Rijn 1984 pickup)
- van Rijn, L.C., "Simple General Formulae for Sand Transport in Rivers, Estuaries and Coastal Waters". https://www.leovanrijn-sediment.com/papers/Formulaesandtransport.pdf (primary-author note, read in full: Rouse profile, suspension criteria, critical velocities)
- van Rijn, L.C. (1984), "Sediment Transport, Part II: Suspended Load Transport", J. Hydraulic Engineering 110(11), 1613-1641. https://ascelibrary.org/doi/10.1061/%28ASCE%290733-9429%281984%29110%3A11%281613%29
- Soulsby, R. (1997), "Dynamics of Marine Sands", Thomas Telford, London — settling formula & hindered variant as documented in HEC-RAS 2D Sediment Transport docs: https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.5/model-description/fall-velocity-and-settling/particle-settling-velocity and https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.2/model-description/fall-velocity-and-settling/hindered-settling
- FUNWAVE-TVD documentation, "Suspended Sediment Transport Equation (Non-cohesive)" — van Rijn pickup + Cao deposition as implemented in a working model. https://fengyanshi.github.io/build/html/sed_equation.html
- Richardson, J.F. & Zaki, W.N. (1954), "Sedimentation and Fluidisation: Part I", Trans. Inst. Chem. Eng. 32, 35-53.
- Gualtieri et al. / APS: "Controversial turbulent Schmidt number value in particle-laden flows", Phys. Rev. Fluids 7, 014307. https://link.aps.org/accepted/10.1103/PhysRevFluids.7.014307
- COHERENS ocean model documentation, Chapter 7 "Sediment transport model". https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf

## Verification

Independent adversarial fact-check (2026-07-07). Method: every formula was re-derived numerically (all table
values recomputed from scratch) and every constant was checked against primary sources fetched and read
independently of the citations above: van Rijn (1984) J. Hydraul. Eng. 110(11):1613-1641 scanned original
(https://www.leovanrijn-sediment.com/papers/P2-1984b.pdf), van Rijn/Bisschop/van Rhee (2019) author PDF
(https://www.leovanrijn-sediment.com/papers/Pickup2019.pdf), van Rijn "Simple general formulae" note
(https://www.leovanrijn-sediment.com/papers/Formulaesandtransport.pdf), van Rijn (2007) Unified View II
(https://www.leovanrijn-sediment.com/papers/P2-2007b.pdf), Ferguson & Church (2004) JSR 74(6):933-937 full PDF,
COHERENS ch.7 (Eqs. 7.41, 7.116, 7.120-7.125, 7.128, 7.135-7.136, 7.142), van Rijn & Walstra "Modelling of
sand transport in Delft3D" (WL|Delft Z3624, Eqs. 2.2.32, 2.2.41, 2.2.68), HEC-RAS 2D sediment docs,
Miedema Slurry Transport ch. 4.6 (LibreTexts), and APS Phys. Rev. Fluids 7, 014307 abstract.

1. **Soulsby (1997) settling velocity — CONFIRMED.** w_s = (nu/d)[sqrt(10.36^2 + 1.049 D*^3) - 10.36] with
   10.36^2 = 107.33 verified verbatim in HEC-RAS 2D docs and COHERENS Eq. 7.41 (both give constants 10.36 and
   1.049; D* = d[(s-1)g/nu^2]^(1/3)). Caveat: Soulsby's fit is calibrated on natural marine sands ~0.063-2 mm;
   use at 10 mm is extrapolation on the constant-drag asymptote (Sec. 3 of this note already states this — the
   formula stays within ~11% of F&C there, per the recomputed Sec. 4 table). Hindered variant
   (1.049 -> 1.049(1-C)^4.7): functional form and the 1.049 constant confirmed via secondary sources (HEC-RAS
   hindered-settling page names Soulsby 1997 as the only formula with built-in hindered settling; USACE DTIC
   AD1002918 gives X2 = 1.049 with hindered exponent in the 2.5-5 range); the exact exponent 4.7 matches
   Soulsby's recommended Richardson-Zaki n for sand but could not be reproduced verbatim from an accessible
   copy of the book — treat 4.7 as verified-by-consensus, not primary-verified.
2. **Ferguson & Church (2004) — CONFIRMED (primary PDF read in full).** Eq. 4: w = R g d^2/(C1 nu +
   sqrt(0.75 C2 R g d^3)). Constants exactly as stated: (18, 0.4) smooth spheres; (18, 1.0) intermediate
   "grains of varied shape", best with sieve diameter (rmse 16% on Hallermeier's 115-point set); (20, 1.1)
   Dietrich-equivalent natural grains with nominal diameter; (24, 1.2) angular extreme. Fraser River
   experiments: d = 0.068-4.36 mm, rms error 6%. Raudkivi-data rms errors 7-9% for natural-parameter
   combinations. The note's "rms error 6-16%" is correct.
3. **Richardson-Zaki (1954) — CONFIRMED, one boundary caveat.** w_s,m = w_s,0(1-C)^n valid C = 0-0.3
   (Miedema/LibreTexts, reproducing the original R&Z table): n = 4.65 for Re_p < 0.2; n = 4.35 Re^-0.03
   (0.2-1); n = 4.45 Re^-0.1 (1-200); n = 2.39 for large Re. Caveat: the constant-2.39 regime starts at
   Re_p > 200 (original table; 400 in the revised continuous fit; ">500" appears in later restatements) —
   the value ~2.4 is right, the note's Re > 500 boundary is the loosest of the quoted variants. Irrelevant
   for implementation since we fix n = 4.7 (Soulsby) / n = 4-5 (van Rijn et al. 2019, primary-verified:
   "n = exponent (range of 4 to 5) according to Richardson and Zaki (1954)").
4. **Suspension criterion — CONFIRMED (primary scan read).** van Rijn 1984-II Eqs. 8-9 verbatim:
   u*,crs/w_s = 4/D* for 1 < D* <= 10 and u*,crs/w_s = 0.4 for D* > 10 (Bagnold upper limit u*/w_s = 1,
   Eq. 5). Shields-style curve verbatim in van Rijn's Simple-general-formulae note Eq. 3.2:
   theta_cr,suspension = 0.3/(1+D*) + 0.1[1-exp(-0.05 D*)] (and Eq. 3.1 initiation-of-motion
   theta_cr = 0.3/(1+1.2 D*) + 0.055[1-exp(-0.02 D*)] as used in Sec. 9). All Sec. 6 u*_cr,susp values
   recomputed and match to the printed precision.
5. **Rouse profile & number — CONFIRMED (primary scan read).** Z = w_s/(beta kappa u*) is 1984-II Eq. 3;
   beta = 1 + 2(w_s/u*)^2 for 0.1 < w_s/u* < 1 is 1984-II Eq. 22 verbatim (note: in the paper beta itself
   reaches ~3 at w_s/u* = 1 — the clip to [1, 1.5] is a model convention, verified verbatim in Delft3D
   Z3624 Eq. 2.2.32 ("limited to the range 1 < beta < 1.5") and COHERENS Eq. 7.136 ("limited to values
   between 1 and 1.5")). Profile form and Z = 5/2/1/0.1 interpretation verbatim in van Rijn's note
   Eq. 3.14. kappa = 0.4 confirmed.
6. **3D advection-diffusion — equation CONFIRMED, one citation CORRECTED.** Equation, settling as extra
   downward advection, eps_s = nu_t/sigma_s (= beta nu_t), and zero-total-flux surface BC
   (eps_s dc/dz + w_s c = 0) verified verbatim against COHERENS Eqs. 7.116, 7.135, 7.142. sigma_s = 0.7 as
   engineering default stands. CORRECTION: the "measured range 0.2-1.0" should not be attributed to
   Phys. Rev. Fluids 7, 014307 — that paper argues the opposite: direct flux measurements + two-phase
   simulations give sigma_s ~ 3-4 (> 1) for sand in boundary layers, and attribute the classical < 1 values
   to Rouse-profile fits using quiescent-water w_s. The 0.1-1.3 (typically 0.2-1.0) range comes from the
   earlier literature reviewed by Gualtieri et al. (2017, Fluids 2(2):17). Practical impact: sigma_s = 0.7
   remains the defensible default (consistent with van Rijn beta in [1,1.5] near the bed, i.e. sigma_s in
   [0.67,1]), but treat sigma_s as an uncertainty knob spanning ~0.5-1.5 rather than a measured constant,
   and do not cite PRF 7:014307 as support for values < 1.
7. **Reference concentration & pickup — CONFIRMED (primary sources read).** 1984-II Eq. 38 verbatim:
   c_a = 0.015 (D50/a) T^1.5/D*^0.3, volumetric ("solids volume per unit fluid volume, or kg/m3 after
   multiplying by rho_s"); Eq. 37 verbatim: a = 0.5*Delta or a = k_s with a_min = 0.01d. Cap: van Rijn
   (2007-II) Eq. 5 gives "c_a <= 0.05" explicitly (Delft3D/COHERENS additionally cap a at 0.20h/0.10h).
   Pickup: 2019 paper Eq. 1 verbatim E = 0.00033 rho_s[(s-1)g d50]^0.5 D*^0.3 T^1.5 kg/m2/s (calibrated
   0.5-1.5 m/s, d50 = 100-1500 um) and Eq. 5 verbatim E = alpha rho_s[(s-1)g d50]^0.5 D*^0.3 f_D T^1.5,
   alpha = 0.00033 (+-30%), f_D = 1/theta' for theta' > 1 (else 1), recalibrated on 50-560 um sand at
   1.5-6 m/s, "only valid for clean... sand... porosity in the range of 0.4+-0.03", "within a factor of 2",
   underprediction ~50% for 150 um (turbidity-current hindcast) and factor ~2 for 50-125 um at > 4 m/s.
8. **Bed exchange & Exner — CONFIRMED.** Deposition = w_s c_b and net bed flux E/rho_s - w_s c_b verified
   verbatim in COHERENS Eqs. 7.120-7.122 (D_n = w_s,n(a) c_n(a); E, D in m/s for volumetric c). The 2019
   paper's Eqs. 3-4b (read) give the exact mass balance E = (1 - n_i - c_nb) v_e rho_s + w_s,m c_nb rho_s,
   i.e. (1-p) dz_b/dt = w_s c_b - E/rho_s with the dilute approximation c_nb << 1-p (our c cap 0.05 vs
   1-p = 0.6 makes the neglected term <= 8% of the erosion-velocity coefficient — acceptable; keep the
   (1 - p - c_b) form in mind if c_b is allowed to grow). Hindered w_s in the deposition flux matches the
   2019 usage (w_s,m = w_s,0(1-c_nb)^n). Porosity p = 0.4 consistent with 2019's 0.4+-0.03.

Numerical spot-checks (all recomputed independently): D* coefficient 24162/m (note: 24160 ✓); the entire
Sec. 4 w_s table reproduced to the last printed digit for Stokes, F&C (18,1.0) and Soulsby; cold-water values
0.0057/0.0202 m/s ✓ and <3% change for d >= 1 mm ✓; hindered factors 0.995/0.79/0.61 at C = 0.001/0.05/0.1 ✓;
Sec. 6 u*_cr,susp 0.012-0.126 m/s ✓; Sec. 7 example Z = 1.22 ≈ 1.2 ✓; 10.36^2 = 107.33 ✓.
