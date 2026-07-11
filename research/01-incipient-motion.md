# Incipient Motion & Critical Bed Shear Stress

## Overview

Incipient motion is the threshold at which the fluid shear stress on the bed, tau_b, first mobilizes
grains. Everything downstream in the sediment pipeline (bedload flux, pickup/erosion of suspended load,
Exner bed update, and the self-healing "sand settles inside the shape" behaviour COBOD is designing for)
switches on/off at this threshold, so it must be computed per-cell, per-timestep, and corrected for the
local bed slope. The industry-standard framework is the Shields (1936) curve, evaluated in explicit form
via the Soulsby–Whitehouse (1997) fit as a function of the dimensionless grain size D*. This is exactly
what Delft3D, COHERENS, HEC-RAS 2D, and most scour CFD (e.g., Roulund et al. 2005) use. All formulas
below use SI units. Project fluid constants: rho = 1025 kg/m3 (seawater), rho_s = 2650 kg/m3 (quartz),
s = rho_s/rho = 2.585, g = 9.81 m/s2, nu = 1.05e-6 m2/s (project spec; note that seawater at 10 C /
35 ppt actually has nu ~ 1.36e-6 m2/s per Soulsby (1997) — see sensitivity table).

## Key equations & results

### 1. Shields parameter (Shields 1936; Soulsby 1997)

    theta = tau_b / ((rho_s - rho) * g * d)          [dimensionless]

- tau_b = bed shear stress [Pa] = rho * u*^2, with u* = friction velocity [m/s].
- d = grain diameter [m] (use d50 for a single-fraction bed).
- Motion occurs when theta > theta_cr. Source: Shields (1936); Soulsby, "Dynamics of Marine Sands" (1997).

### 2. Dimensionless grain size D* (Soulsby 1997; van Rijn 1984)

    D* = d * [ g * (s - 1) / nu^2 ]^(1/3)            with s = rho_s / rho

For our constants (nu = 1.05e-6): [g(s-1)/nu^2]^(1/3) = 24162 per metre, so D* = 24162 * d.
For seawater at 10 C (nu = 1.36e-6): factor = 20335 per metre. D* removes u* from the classic Shields
diagram (which used grain Reynolds number Re* = u* d / nu) so no iteration is needed.
Source: Soulsby (1997); COHERENS sediment manual eq. 7.29; Bosboom & Stive, Coastal Dynamics 6.3.2.

### 3. Soulsby–Whitehouse (1997) explicit critical Shields parameter — RECOMMENDED

    theta_cr = 0.30 / (1 + 1.2 * D*) + 0.055 * (1 - exp(-0.020 * D*))

Valid for all D* (fit to the full Shields data set plus fine-sediment data; intended range roughly
0.1 < D* < 10000, i.e., silt to cobbles). Asymptotes: theta_cr -> 0.30 as D* -> 0 (raised above Shields'
original curve to fit fine-sand/silt data); theta_cr -> 0.055 as D* -> infinity. Minimum theta_cr ~ 0.030
near D* ~ 15-25 (d50 ~ 0.6-1.0 mm in seawater).
Source: Soulsby & Whitehouse (1997), "Threshold of sediment motion in coastal environments", Pacific
Coasts and Ports '97; reproduced in COHERENS manual eq. 7.32 and HEC-RAS 2D sediment documentation.

DISAGREEMENT NOTE: van Rijn's design note (leovanrijn-sediment.com, 2018) writes the first term as
0.3/(1 + D*) (no 1.2 factor), which raises theta_cr by ~10-14 % for d50 <= 0.2 mm and is negligible
above 1 mm. Brownlie (1981), theta_cr = 0.22*D*^-0.9 + 0.06*10^(-7.7*D*^-0.9), gives +25 % at 0.1 mm,
+4-11 % at 0.2-2 mm, -2 to -4 % at 5-10 mm relative to Soulsby-Whitehouse. Use Soulsby-Whitehouse as
the single implementation; the spread (~10-25 % at the fine end) is the genuine physical uncertainty band.

### 4. Critical shear stress and critical friction velocity

    tau_cr = theta_cr * (rho_s - rho) * g * d        [Pa]
    u*_cr  = sqrt(tau_cr / rho)                      [m/s]

For rho_s - rho = 1625 kg/m3: tau_cr = 15942 * theta_cr * d [Pa, d in m].

Tabulated values, quartz in seawater, nu = 1.05e-6 m2/s (project baseline):

| d50 [mm] | D*    | theta_cr [-] | tau_cr [Pa] | u*_cr [m/s] | U_cr depth-avg, h=5 m [m/s] |
|----------|-------|--------------|-------------|-------------|------------------------------|
| 0.1      | 2.42  | 0.0795       | 0.127       | 0.0111      | 0.37                         |
| 0.2      | 4.83  | 0.0492       | 0.157       | 0.0124      | 0.37                         |
| 0.5      | 12.1  | 0.0312       | 0.248       | 0.0156      | 0.41                         |
| 1.0      | 24.2  | 0.0311       | 0.495       | 0.0220      | 0.52                         |
| 2.0      | 48.3  | 0.0392       | 1.249       | 0.0349      | 0.75                         |
| 5.0      | 120.8 | 0.0521       | 4.156       | 0.0637      | 1.20                         |
| 10.0     | 241.6 | 0.0556       | 8.862       | 0.0930      | 1.58                         |

Sensitivity to viscosity — same table with nu = 1.36e-6 m2/s (real 10 C seawater):
theta_cr = 0.0894 / 0.0553 / 0.0328 / 0.0302 / 0.0366 / 0.0502 / 0.0553 and
tau_cr = 0.143 / 0.176 / 0.262 / 0.481 / 1.168 / 4.005 / 8.813 Pa for the same seven sizes.
Differences are <= 13 % (largest at 0.1 mm); coarse grains are viscosity-insensitive.
The last column is Soulsby (1997)'s threshold depth-averaged current speed,
U_cr = 7 * (h/d)^(1/7) * sqrt(g*(s-1)*d*theta_cr), h = water depth [m] — use it only as a sanity check
of the simulator (with 0.5-2.5 m/s tidal currents, ALL of our d50 range is mobile at the upper end;
5-10 mm gravel is immobile below ~1.2-1.6 m/s).

### 5. Hjulström–Sundborg curve (context only — do NOT implement)

Empirical curve of depth-averaged velocity (defined at ~1 m flow depth) vs grain size with three regimes:
erosion (above upper curve), transport of already-moving grains (between curves), deposition (below lower
curve = settling velocity). Approximate erosion thresholds: minimum ~0.15-0.25 m/s at d ~ 0.1-0.5 mm;
~0.5 m/s at 1 mm; ~1.5-2 m/s at 10 mm; increases again below 0.06 mm because silt/clay are cohesive
(e.g., ~1 m/s for consolidated clay). Deposition thresholds ~ settling velocity: ~0.008 m/s (0.1 mm),
~0.1 m/s (1 mm). Limitations: no depth dependence, mixes velocity with shear stress, historical value
only — it is consistent with but superseded by the Shields/Soulsby framework above (compare the U_cr
column). The erosion-vs-deposition hysteresis it displays IS physically real and emerges automatically
in our model because pickup requires tau_b > tau_cr while settling continues at any tau_b.
Source: Hjulström (1935), Sundborg (1956); en.wikipedia.org/wiki/Hjulström_curve.

### 6. Bed-slope corrections to theta_cr — REQUIRED for scour holes (slopes reach the repose angle)

General 3D formula, Soulsby (1997) eq. (80a), equivalently Roulund et al. (2005) / Fredsøe & Deigaard (1992):

    theta_cr(beta, psi) = theta_cr,flat * [ cos(psi)*sin(beta)
                          + sqrt( cos(beta)^2 * tan(phi_i)^2 - sin(psi)^2 * sin(beta)^2 ) ] / tan(phi_i)

- beta = local bed slope angle relative to horizontal [rad]; psi = angle between the near-bed flow vector
  and the UP-slope (steepest-ascent) direction; phi_i = internal friction / repose angle.
- phi_i = 30-35 deg for sand (use 32 deg, tan(phi_i) = mu_s = 0.63 as in Roulund et al. 2005),
  35-40 deg for angular gravel (van Rijn 2018 note: 30-40 deg).
- Special cases (all recovered exactly from the general formula):
  - Flow directly upslope (psi = 0):    factor = sin(phi_i + beta) / sin(phi_i)   [Damgaard, Whitehouse & Soulsby 1997, J. Hydraulic Eng. 123]
  - Flow directly downslope (psi = 180 deg): factor = sin(phi_i - beta) / sin(phi_i)
  - Pure transverse slope (psi = 90 deg):    factor = cos(beta) * sqrt(1 - tan(beta)^2 / tan(phi_i)^2)   [Lane 1955; van Rijn 1993]
- Numeric examples (phi_i = 32 deg): beta = 10 deg -> up 1.26x, down 0.71x, transverse 0.95x;
  beta = 20 deg -> up 1.49x, down 0.39x, transverse 0.76x; beta = 30 deg -> up 1.67x, down 0.07x,
  transverse 0.33x. Factor -> 0 as beta -> phi_i (grains roll without flow): clamp theta_cr >= 0 and add
  a sand-slide/avalanche routine that relaxes any bed face steeper than phi_i (Roulund et al. 2005 trigger
  avalanching at slope > 32 deg, stop at ~31 deg).
- Multiplicative composition used by van Rijn (1993/2018) when treating slopes as separate longitudinal
  (angle beta2, must be < phi_i) and transverse (angle beta1) components:
  theta_cr = K_long * K_trans * theta_cr,flat with the two special-case factors above. The general
  Soulsby formula is preferred on a voxel grid because each cell has one slope vector, not two angles.
Source: Soulsby (1997) "Dynamics of Marine Sands" sec. 6.3; Roulund et al. (2005) J. Fluid Mech. 534;
van Rijn note "Design of bed protections" (2018); Dey (2003) gives an alternative empirical fit with the
same limiting behaviour — differences are second-order compared to the choice of phi_i.

### 7. Hiding–exposure correction for sand–gravel mixtures (multi-fraction beds only)

For fraction i with diameter d_i in a mixture with median d50, multiply theta_cr(d_i) by xi_i:

Egiazaroff (1965), as adapted by Ashida & Michiue (1972) — COHERENS manual eq. 7.37:

    xi_i = [ log10(19) / ( log10(19) + log10(d_i/d50) ) ]^2      for d_i/d50 >= 0.38889
    xi_i = 0.8429 * (d50 / d_i)                                  for d_i/d50 <  0.38889

xi_i > 1 for fines hidden between coarse grains (up to ~2x), xi_i < 1 for exposed coarse grains
(floor ~0.5-0.6). Alternative: Wu et al. (2000), xi_i = (p_e,i/p_h,i)^-0.6 with hidden/exposed
probabilities p_h,i = sum_m(p_m * d_m/(d_i+d_m)), p_e,i = 1 - p_h,i (p_m = mass fraction of class m).
Field data (Ferret et al. 2018, Marine Geology) show mixture effects can raise sand tau_cr by up to
~75 % and lower gravel tau_cr by up to ~64 % vs uniform beds — a leading-order effect if we ever
simulate bimodal beds, but safely omitted for single-fraction runs.
Source: Egiazaroff (1965); Ashida & Michiue (1972); Wu, Wang & Jia (2000); COHERENS manual eqs. 7.34-7.37.

## Practical guidance for our simulator

1. Implement Soulsby–Whitehouse (eq. 3) as THE theta_cr formula — one branch-free line per cell,
   matches Delft3D/COHERENS/HEC-RAS practice, valid across our whole d50 = 0.1-10 mm range
   (D* = 2.4-242 at nu = 1.05e-6). Precompute D*, theta_cr, tau_cr once per sediment fraction at
   startup (they depend only on constants), store tau_cr,flat as a scalar; only the slope correction
   varies per cell/timestep.
2. Constants to hard-code: g = 9.81, rho = 1025, rho_s = 2650, s = 2.585, rho_s - rho = 1625 kg/m3.
   Make nu a config parameter: 1.05e-6 (project baseline, ~20 C freshwater-like) vs 1.36e-6 (real 10 C
   seawater). It changes tau_cr by <= 13 % at the fine end — expose it, don't bury it.
3. Reference check values for unit tests (nu = 1.05e-6): d50 = 0.5 mm -> D* = 12.08,
   theta_cr = 0.0312, tau_cr = 0.248 Pa, u*_cr = 0.0156 m/s; d50 = 10 mm -> D* = 241.6,
   theta_cr = 0.0556, tau_cr = 8.86 Pa, u*_cr = 0.0930 m/s. Curve minimum theta_cr ~ 0.030 at ~0.8 mm.
4. Apply the slope correction (eq. 6, general Soulsby form) per bed cell each morphodynamic step:
   compute bed normal from the voxel heightfield/level-set gradient, get beta and psi from the near-bed
   velocity vector, use phi_i = 32 deg (sand) with tan(phi_i) = 0.6249. Clamp the corrected theta_cr to
   [0, 2*theta_cr,flat]. Pair it with an avalanching pass (trigger > 32 deg, relax to 31 deg) — without
   it, scour-hole side slopes go unphysically steep and the slope factor goes negative.
5. The threshold defines onset only; excess stress (theta - theta_cr) drives the bedload and pickup
   formulas (separate research notes). Use the SAME theta_cr in both so erosion and transport switch
   together; the deposition flux (settling) is NOT thresholded — that asymmetry reproduces the
   Hjulström transport band and is what lets sand accumulate inside the printed shapes even while the
   free bed nearby erodes.
6. Expect any turbulence-resolving LES to intermittently exceed tau_cr: the Shields curve already
   represents "frequent movement" under turbulent bursts (Deltares 1972 via van Rijn: at the curve,
   ~100 % of locations show frequent particle movement; r = 0.4*theta_cr corresponds to occasional
   movement). Drive erosion with the time-filtered tau_b from the wall model, not instantaneous LES
   fluctuations, or the bed will erode too fast. If a "no damage" criterion is ever needed for the
   printed-armour stability itself, van Rijn's damage-level reduction r = 0.4-0.8 applies.
7. Skip hiding-exposure (eq. 7) for milestone 1 (single d50 per run); implement Egiazaroff/Ashida-Michiue
   as a per-fraction multiplier if/when multi-fraction beds are added.

## Sources

- Soulsby, R.L. (1997). Dynamics of Marine Sands: A Manual for Practical Applications. Thomas Telford, London. https://eprints.hrwallingford.com/412/
- Soulsby, R.L. & Whitehouse, R.J.S. (1997). Threshold of sediment motion in coastal environments. Proc. Pacific Coasts and Ports '97, Christchurch, 149-154.
- COHERENS documentation, Chapter 7: Sediment transport model (eqs. 7.27-7.37, 7.81-7.85). https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf
- HEC-RAS 2D Sediment Transport: Critical Thresholds for Transport and Erosion. https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.6/model-description/critical-thresholds-for-transport-and-erosion
- van Rijn, L.C. (2018). Design of bed protections; stability and movement of cobbles, boulders and rocks. https://www.leovanrijn-sediment.com/papers/Bedprotections2017.pdf
- Damgaard, J.S., Whitehouse, R.J.S. & Soulsby, R.L. (1997). Bed-load sediment transport on steep longitudinal slopes. J. Hydraulic Engineering 123(12), 1130-1138.
- Roulund, A., Sumer, B.M., Fredsøe, J. & Michelsen, J. (2005). Numerical and experimental investigation of flow and scour around a circular pile. J. Fluid Mechanics 534, 351-401.
- Roulund, A., Sutherland, J., Todd, D. & Sterner, J. (2016). Parametric equations for Shields parameter and wave orbital velocity in combined current and irregular waves. ICSE 2016. https://eprints.hrwallingford.com/1134/1/Roulund_ICSE2016.pdf
- Maldonado & Borthwick (2016). Quasi-two-layer morphodynamic model: bed slope-induced morphological diffusion. R. Soc. Open Science. https://arxiv.org/pdf/1607.05820
- Bosboom, J. & Stive, M.J.F. Coastal Dynamics, sec. 6.3.2 Shields curve. https://geo.libretexts.org/Bookshelves/Oceanography/Coastal_Dynamics_(Bosboom_and_Stive)/06:_Sediment_transport/6.03:_Initiation_of_motion/6.3.2:_Shields_curve
- Hjulström curve. https://en.wikipedia.org/wiki/Hjulström_curve
- Ferret, Y. et al. (2018). The hiding-exposure effect revisited. Marine Geology. https://www.sciencedirect.com/science/article/pii/S0025322718302226

## Verification

Adversarial fact-check (2026-07-07). Every equation was re-derived/recomputed numerically (double
precision) and checked against sources located independently of the citations above: the full text of
Soulsby (1997) *Dynamics of Marine Sands* (pdfcoffee.com full-text copy), the TELEMAC/SISYPHE v6.3
reference manual (EDF, opentelemac.org), the Delft3D morphology-kernel Fortran source
(svn.oss.deltares.nl, `comphidexp.f90`), the COHERENS sediment-transport chapter (extracted PDF text),
Miedema *Slurry Transport* sec. 5.1 (eng.libretexts.org), HEC-RAS 2D sediment docs, Coastal Wiki, the
Inch Cape Offshore EIA Annex 10A.4 bed-shear-stress methodology (marine.gov.scot), and Maldonado &
Borthwick (R. Soc. Open Sci., ar5iv copy of arXiv:1607.05820).

1. **Shields parameter — CONFIRMED.** theta = tau_b/((rho_s-rho) g d), tau_b = rho u*^2, motion for
   theta > theta_cr: identical definition in Soulsby (1997) sec. 6.4, COHERENS ch. 7, Miedema 5.1, and
   ISOPE-2013 rock-berm paper (cantab.net) eq. (1). Pedantic constant check: (rho_s-rho)*g =
   1625 x 9.81 = 15941.25 exactly, so "15942" is over-rounded by 0.75 (5e-5 relative); use 15941.25 in
   code. All tabulated values in this note were computed from the exact product and are unaffected.
2. **Dimensionless grain size D* — CONFIRMED.** D* = d[g(s-1)/nu^2]^(1/3): Soulsby (1997) eq. (75)
   (quoted verbatim in the full text and cited by eq. number in the Inch Cape annex and Henry et al.,
   ICCE 2012), COHERENS eq. (7.29), Coastal Wiki. Recomputed factors: 24162.4 m^-1 (nu = 1.05e-6),
   20334.8 m^-1 (nu = 1.36e-6); the note's 24162/20335 are correct roundings. s = 2650/1025 = 2.58537.
   Soulsby's own worked examples use nu = 1.36e-6 m2/s and s = 2.58 for 10 C / 35 ppt seawater,
   confirming the sensitivity-table premise.
3. **Soulsby-Whitehouse theta_cr — CONFIRMED.** theta_cr = 0.30/(1+1.2 D*) + 0.055[1-exp(-0.020 D*)]
   verified verbatim, all four constants, in five independent renderings: Soulsby (1997) eq. (77) full
   text; Inch Cape EIA annex ("Equation 77 in Soulsby (1997)"); COHERENS eq. (7.32); HEC-RAS 2D
   sediment docs eq. (3); Miedema/LibreTexts 5.1. Book states validity "any non-cohesive sediment...
   D* > 0.1, valid in any units". Asymptotes 0.30 (D*->0) and 0.055 (D*->inf) follow algebraically;
   recomputed curve minimum theta_cr = 0.0299 at D* = 17.1 (note's "~0.030 at D* ~ 15-25" is right).
4. **tau_cr / u*_cr values — CONFIRMED.** All 7 (theta_cr, tau_cr, u*_cr) triples at nu = 1.05e-6 and
   all 7 sensitivity values at nu = 1.36e-6 reproduced to the last printed digit by direct recomputation
   (e.g. d50 = 0.5 mm: D* = 12.08, theta_cr = 0.03116, tau_cr = 0.2484 Pa, u*_cr = 0.01557 m/s;
   d50 = 10 mm: theta_cr = 0.05559, tau_cr = 8.862 Pa, u*_cr = 0.0930 m/s).
5. **General bed-slope correction — CONFIRMED.** The identical formula, including psi defined as the
   angle of the current to the UP-slope direction, is implemented as SLOPEFF=2 in SISYPHE v6.3 manual
   eq. (27): theta_c/theta_c0 = [cos psi sin chi + (cos^2 chi tan^2 phi_s - sin^2 psi sin^2 chi)^0.5]
   / tan phi_s. Soulsby (1997) full text sec. 6.4 carries it as eq. (80a) with special cases (80b-d)
   and the avalanching remark for slopes > phi_i. Equivalence to Roulund et al. (2005) verified
   algebraically: their form cos(beta)sqrt(1 - sin^2(alpha)tan^2(beta)/mu_s^2) - cos(alpha)sin(beta)/mu_s
   with alpha measured from the DOWN-slope direction maps exactly onto Soulsby's with psi = 180 - alpha,
   mu_s = tan(phi_i); tan(32 deg) = 0.62487. One caveat: "factor -> 0 as beta -> phi_i" holds for flow
   with a downslope component (90 <= psi <= 180 deg, incl. pure transverse) — exactly the scour-hole
   case — but NOT for upslope flow, where the factor tends to sin(2 phi_i)/sin(phi_i) = 2 cos(phi_i)
   ~ 1.70; the avalanching pass handles that regime regardless.
6. **Slope special cases — CONFIRMED.** Up/down/transverse factors recovered exactly from the general
   eq. (80a) by direct substitution (recomputed at phi_i = 32 deg, beta = 20 deg: 1.487 / 0.392 / 0.764,
   matching 1.49/0.39/0.76; beta = 10 deg: 1.263/0.707/0.945; beta = 30 deg: 1.666/0.066/0.331).
   Independent quotes: Soulsby (1997) eqs. (80b) sin(phi_i+beta)/sin(phi_i) and (80c)
   sin(phi_i-beta)/sin(phi_i); Maldonado & Borthwick eq. (15) sin(phi+beta)/sin(phi). The transverse
   case is the classical Lane (1955) side-slope factor and follows exactly from eq. (80a) at psi = 90.
7. **Hiding-exposure (Egiazaroff / Ashida-Michiue; Wu et al.) — CONFIRMED.** Delft3D morphology kernel
   (comphidexp.f90) implements exactly: dd = d_i/d_m; if dd < 0.38889 then xi = 0.8429/dd else
   xi = [log10(19)/(log10(19)+log10(dd))]^2, with log10(19) = 1.27875360; COHERENS eq. (7.37) carries
   the same constants. Continuity at the cutoff verified numerically: both branches give xi = 2.1675 at
   d_i/d50 = 0.38889 (the historical Ashida-Michiue pair 0.85/0.40 is discontinuous by ~1%; 0.8429 &
   0.38889 is the continuity-exact variant used by Delft3D/COHERENS — implement these). Beware: the
   COHERENS manual prints the branch condition with the ratio inverted (d50/dn < 0.38889); the Delft3D
   source and this note's form (condition on d_i/d50) are the correct, continuous ones. Wu et al. (2000):
   xi = (p_e/p_h)^m with m = -0.6, p_h,i = sum_m p_m d_m/(d_i+d_m), p_e,i = 1 - p_h,i confirmed in both
   COHERENS eqs. (7.34-7.36) ("which Wu et al. (2000) determined as m = -0.6") and the Delft3D source.
8. **Threshold depth-averaged current U_cr — CONFIRMED.** Soulsby (1997) sec. 6.2, eq. SC(72a,b): the
   full text confirms the "Soulsby formula for threshold current speed" is eq. (77) combined with the
   1/7th-power friction law (eq. 34), "valid for any non-cohesive sediment and water conditions for
   which D* > 0.1, and valid in any units", with eq. (72b) = the S-W theta_cr expression under the
   square root. OCR of the equation body drops the leading factor, so the exact form was validated
   numerically instead: U_cr = 7(h/d50)^(1/7) sqrt(g(s-1) d50 theta_cr) reproduces (i) the note's whole
   U_cr column (0.37/0.37/0.41/0.52/0.75/1.20/1.58 m/s at h = 5 m, nu = 1.05e-6) and (ii) the book's own
   Example 6.1 (d50 = 200 um, h = 5 m, 10 C/35 ppt seawater -> U_cr = 0.39 m/s; this recomputation with
   nu = 1.36e-6 gives 0.39 m/s) and is consistent with van Rijn (1984) eq. (71a) (~0.38 m/s, same case).
   Use as a sanity check only, as the note says (it assumes an unrippled flat bed and 1/7 profile).

Verification sources (independent of the citations above):
- Soulsby (1997) full text: https://pdfcoffee.com/dynamics-of-marine-sands-a-manual-for-practical-applications-compress-pdf-free.html
- SISYPHE v6.3 reference manual (eq. 27, SLOPEFF=2): https://www.opentelemac.org/downloads/MANUALS/SISYPHE/sisyphe_user_manual_en_v6p3.pdf
- Delft3D morphology kernel source (Egiazaroff/Ashida-Michiue/Wu): https://svn.oss.deltares.nl/repos/delft3d/trunk/src/utils_gpl/morphology/packages/morphology_kernel/src/comphidexp.f90
- Miedema, Slurry Transport 5.1 (S-W + Brownlie forms): https://eng.libretexts.org/Bookshelves/Civil_Engineering/Slurry_Transport_(Miedema)/05:_Initiation_of_Motion_and_Sediment_Transport/5.01:_Initiation_of_Motion_of_Particles
- Inch Cape Offshore Wind EIA, Annex 10A.4 (quotes Soulsby eqs. 30/37/57/69/70/77 by number): https://marine.gov.scot/datafiles/lot/inch_cape/Environmental%20Statement/Volume%202B%20-%20Appendices/Annex%2010A.4%20-%20Bed%20Sheer%20Stress%20Analysis%20Methodology.pdf
- Thusyanthan, ISOPE-2013, Stability of Rock Berm (Shields definition + S-W): https://www.cantab.net/users/Dr_Thusyanthan/41_2013_Dr_Thusyanthan_Stability%20of%20Rock%20Berm%20under%20Wave%20and%20Current%20Loading_ISOPE2013_TPC_0621.pdf
- Maldonado & Borthwick, R. Soc. Open Sci. (slope factor eq. 15): https://ar5iv.labs.arxiv.org/html/1607.05820
- Coastal Wiki (D* definition): https://www.coastalwiki.org/wiki/Sediment_transport_formulas_for_the_coastal_environment
- Henry et al., ICCE 2012 (cites Soulsby 1997 "Eq. 75" for D*): https://icce-ojs-tamu.tdl.org/icce/index.php/icce/article/download/6510/pdf_513/
