# Bedload Transport & Bed Evolution (Exner Equation)

## Overview

This note specifies the bedload-transport and bed-update layer of the morphodynamic simulator:
(1) dimensionless bedload rate Phi as a function of Shields excess, with the three candidate
formulas (Meyer-Peter & Muller 1948 incl. Wong & Parker 2006 correction; van Rijn 1984;
Engelund-Fredsoe 1976) and their validity ranges over our d50 = 0.1-10 mm span; (2) the Exner
equation with bed porosity; (3) gravity/slope corrections to bedload magnitude and direction;
(4) the avalanching (sand-slide) algorithm that every scour code needs; (5) morphological
acceleration (MORFAC) practice and pitfalls.

Conventions used throughout: qb = volumetric bedload rate per unit width [m2/s], solids volume
only (no pores); tau_b = bed shear stress [N/m2] (skin friction, i.e. grain-related part only);
rho = 1025 kg/m3; rho_s = 2650 kg/m3; s = rho_s/rho = 2.585; g = 9.81 m/s2; nu = 1.05e-6 m2/s
(10 C seawater); d = d50 [m]. Useful derived constants for our water/sand combination:
(s-1)g = 15.55 m/s2; (rho_s - rho)g = 15941 N/m3; D* = 24160 * d (d in metres).

## Key equations & results

### 1. Dimensionless groups (Apsley bedload sheet; van Rijn 1984; Soulsby 1997)

- Dimensionless bedload rate (Einstein parameter):
  Phi = qb / sqrt( (s-1) * g * d50^3 )        [dimensionless]
- Shields parameter:
  theta = tau_b / ( rho * (s-1) * g * d50 ) = u_star^2 / ( (s-1) * g * d50 )
- Dimensionless grain size:
  D* = d50 * ( (s-1) * g / nu^2 )^(1/3)
  For our fluid: D* = 2.42 (d50 = 0.1 mm), 7.2 (0.3 mm), 24.2 (1 mm), 72.5 (3 mm), 242 (10 mm).

### 2. Critical Shields parameter, flat bed (Soulsby & Whitehouse 1997; in Soulsby "Dynamics of Marine Sands", also HEC-RAS 2D manual)

theta_cr = 0.30 / (1 + 1.2*D*) + 0.055 * (1 - exp(-0.020*D*))

Valid for the full silt-to-gravel range. For our sediments (seawater, 10 C):

| d50 [mm] | D* | theta_cr | tau_cr [N/m2] | u*_cr [m/s] |
|---|---|---|---|---|
| 0.1 | 2.4 | 0.080 | 0.13 | 0.011 |
| 0.3 | 7.2 | 0.038 | 0.18 | 0.013 |
| 1.0 | 24 | 0.031 | 0.49 | 0.022 |
| 3.0 | 72 | 0.046 | 2.18 | 0.046 |
| 10.0 | 242 | 0.056 | 8.86 | 0.093 |

(tau_cr = theta_cr * 15941 * d50; u*_cr = sqrt(tau_cr/rho).)

### 3. Meyer-Peter & Muller (1948) and Wong & Parker (2006) correction

Original MPM (Apsley sheet; TELEMAC/Gaia docs):
  Phi = 8 * (theta' - theta_cr)^1.5 ,  theta_cr = 0.047, for theta' > theta_cr else 0.
theta' is the skin-friction Shields number. Database: uniform sediment 0.4-29 mm,
s = 1.25-4.2, mostly theta < ~0.25.

Wong & Parker (2006, J. Hydraulic Eng. 132(11)) reanalysed the MPM plane-bed data and showed
the original form-drag correction was unnecessary; their corrected fits (use either):
  Phi = 3.97 * (theta - 0.0495)^1.5        (exponent fixed at 3/2)
  Phi = 4.93 * (theta - 0.0470)^1.60       (best-fit exponent)
Applies to plane-bed (no dune) transport of uniform sand/gravel with total = skin stress.
The 3.97 version is the recommended default in TELEMAC/Gaia and modern morphodynamic codes.
Note: original MPM with coefficient 8 over-predicts ~2x on plane beds.

### 4. van Rijn (1984, J. Hydraulic Eng. 110(10), his eq. 22)

  Phi = 0.053 / D*^0.3 * T^2.1 ,  T = (theta' - theta_cr) / theta_cr  (transport stage)
Equivalent form: qb = 0.053 * sqrt((s-1)g) * d50^1.5 * T^2.1 / D*^0.3.
theta' uses grain roughness only (Chezy C' = 18*log10(12h/(3*d90))). Calibrated for
0.2 mm <= d50 <= 2 mm (Gaia manual states 0.6-2.0 mm for its implementation); for T >= 3
van Rijn recommends the exponent drop to Phi = 0.100/D*^0.3 * T^1.5. van Rijn's own critical
Shields curve is a 5-segment piecewise function of D* (0.24/D* for D*<=4; 0.14*D*^-0.64 for
4<D*<=10; 0.04*D*^-0.1 for 10<D*<=20; 0.013*D*^0.29 for 20<D*<=150; 0.056 for D*>150) —
nearly identical to Soulsby-Whitehouse; pick ONE curve and use it consistently.

### 5. Engelund & Fredsoe (1976) — the formula used in the reference pile-scour model (Roulund et al. 2005, JFM 534)

Deterministic closed form (Parker CISM notes eq. 3.14):
  Phi = 18.74 * (theta - theta_cr) * ( sqrt(theta) - 0.7*sqrt(theta_cr) ), theta_cr = 0.05.
Full probabilistic version (Fredsoe & Deigaard 1992; used by Roulund et al. 2005):
  Phi = 5 * p * ( sqrt(theta) - 0.7*sqrt(theta_cr) ),
  p = [ 1 + ( (pi/6)*mu_d / (theta - theta_cr) )^4 ]^(-1/4),  mu_d = 0.51 (dynamic friction),
p = fraction of surface grains moving. Advantage: p and theta_cr can be modified locally for
bed slope (this is exactly what Roulund et al. did for pile scour), and it saturates more
realistically at high theta. Derived/validated for sand (~0.2-1 mm); comparable to
Ashida-Michiue Phi = 17*(theta-0.05)*(sqrt(theta)-sqrt(0.05)) validated 0.3-7 mm.

Cross-check example (d50 = 0.3 mm, U = 1.5 m/s, Cd = 0.0025 -> tau_b = 5.77 N/m2,
theta = 1.21): Wong-Parker Phi = 3.97*(1.16)^1.5 = 4.96 -> qb = 4.96 * 2.05e-5 =
1.0e-4 m2/s. At this theta suspended load dominates for fine sand — bedload formulas are
extrapolated beyond their data (theta ~< 0.25-0.5); acceptable for bedload share but the
suspended module must carry the rest.

### 6. Exner equation for bed evolution (sedExnerFoam, GMD 19:2299, 2026; standard)

  (1 - lambda) * dz_b/dt + d(qb_x)/dx + d(qb_y)/dy = D - E

z_b = bed elevation [m]; lambda = bed porosity, use lambda = 0.4 for natural sand
(typical range 0.35-0.45; 1-lambda = 0.6 solids fraction); D, E = suspended-load deposition
and entrainment fluxes [m/s, volume of solids per bed area per time] coupling to the
advection-diffusion module. qb is solids volume flux; if a formula returns bulk volume,
divide by (1-lambda) consistently. Numerics (sedExnerFoam findings): pure central
differencing of div(qb) is oscillatory; use first/second-order upwinding of qb along the
transport direction, explicit Euler or 2nd-order Adams-Bashforth in time, and limit the bed
change per step (keep |dz_b| per morphological step < ~5-10% of the vertical cell size
0.025-0.05 m, i.e. < ~2-5 mm/step).

### 7. Slope effects on bedload magnitude and direction

(a) Threshold modification (Soulsby 1997, "Dynamics of Marine Sands" eq. 80a):
  theta_cr(slope) = theta_cr * [ cos(psi)*sin(beta) + sqrt( cos(beta)^2 * tan(phi_i)^2
                     - sin(psi)^2 * sin(beta)^2 ) ] / tan(phi_i)
beta = local bed slope angle, psi = angle between flow direction and up-slope direction,
phi_i = internal friction angle; Soulsby recommends phi_i = 32 deg for sand (Gaia's default
is 40 deg — too high for loose sand; use 30-32 deg). Downslope flow lowers theta_cr, upslope
raises it; theta_cr -> 0 as beta -> phi_i. This is the physically-complete correction and is
what Roulund et al. (2005) apply (their equivalent Fredsoe & Deigaard 1992 form).

(b) Magnitude correction, longitudinal (Koch & Flokstra 1980; TELEMAC/Gaia default):
  qb_corrected = qb * (1 - beta_KF * dz_b/ds),  beta_KF = 1.3,
dz_b/ds = bed slope along transport direction (positive uphill). Acts like a diffusion term
that smooths bed features.

(c) Direction deviation, transverse (Talmon et al. 1995; TELEMAC/Gaia):
  tan(alpha_new) = tan(alpha_tau) - (1/f(theta)) * dz_b/dn,  f(theta) = beta_2 * sqrt(theta),
beta_2 = 0.85 (Gaia default; calibration studies report optima up to ~1.6). Simpler
equivalent used by Delft3D (Ikeda 1982 / van Rijn 1993):
  qb_n = qb_s * alpha_bn * sqrt(theta_cr/theta) * dz_b/dn, alpha_bn = 1.5 (Delft3D default).
Smaller f(theta) or larger beta_2/alpha_bn = more downslope deflection.

Recommendation: implement (a) + (c); (b) is optional if (a) is in (do not double-count).

### 8. Avalanching / sand-slide (Roulund et al. 2005, JFM; sedExnerFoam 2026)

Submerged angle of repose for sand: phi_s ~ 30-32 deg (tan = 0.58-0.62). Roulund et al.:
activate sand-slide where local slope exceeds phi_s = 32 deg, deactivate once relaxed to
30 deg (hysteresis avoids limit-cycling). Standard grid algorithm (used in Delft3D-derived
and most scour codes), GPU-friendly Jacobi form on the 2D height field:
  for each cell i and each of its 8 neighbours j (distance Delta_ij = dx or dx*sqrt(2)):
    excess = (z_i - z_j) - Delta_ij * tan(phi_repose)
    if excess > 0: transfer dV = 0.5 * k * excess * A_cell from i to j (k ~ 0.5-1.0 relaxation)
  iterate sweeps until max slope <= tan(30 deg); typically < 10-30 iterations/step.
Mass is exactly conserved (antisymmetric transfers). Alternative without iteration
(sedExnerFoam, after Duran Vinent et al. 2019): add a continuous avalanche bedload flux
oriented down the steepest slope,
  |q_av| = q_av0 * [ tanh(tan(beta)) - tanh(tan(phi_r)) ] / [ 1 - tanh(tan(phi_r)) ]
for beta > phi_r, magnitude constant q_av0 chosen large enough to keep slopes near phi_r
(order 10-100x the max physical qb). Either works; the iterative slide is simpler to verify.

### 9. Morphological acceleration factor (MORFAC)

Delft3D practice (D-Morphology manual): multiply the bed-level change from each hydrodynamic
time step by MORFAC = f_mor, so t_morph = f_mor * t_hydro. Sediment fluxes and Exner are
otherwise unchanged. Critical values are regime-dependent (Ranasinghe et al. 2011, Coastal
Eng. 58; ICCE "critical MORFAC" study): unidirectional/steady flow tolerates f_mor up to
O(1000) (reported ~4000 in idealized cases), reversing tidal currents only O(10-200), and
practical coastal work uses 5-20 for tidal, ~1-10 for storms; fluvial studies report ~10%
volume errors even at low MORFAC. Errors arise when the bed changes appreciably within one
hydrodynamic event (tidal phase, vortex-shedding cycle), breaking flow-bed synchronization —
gradients steepen exactly where our simulator cares (scour holes, shoals). Pitfalls: (i) do
not apply MORFAC while the hydrodynamics are still spinning up; (ii) keep
f_mor * |dz_b/step| < ~1% of local depth AND < ~10% of a voxel; (iii) run the sand-slide
AFTER the MORFAC-scaled bed update every step; (iv) for oscillatory/tidal forcing, an
alternative is input-averaging ("representative tide") or Roelvink's mormerge.

## Practical guidance for our simulator

1. Formula choice. Default: Wong & Parker corrected MPM, Phi = 3.97*(theta'-0.0495)^1.5 —
   single power law, no D* dependence, robust from sand to fine gravel (MPM database 0.4-29 mm
   covers our 0.3-10 mm; use with skin-friction theta'). Option 2 (recommended toggle):
   Engelund-Fredsoe with the p-function — it is the formulation validated in the canonical
   pile-scour simulation (Roulund et al. 2005) and integrates cleanly with slope-modified
   theta_cr. Use van Rijn (1984) only for d50 = 0.1-0.3 mm runs where suspended load dominates
   and consistency with van Rijn's reference concentration (suspended module) is wanted.
   Implement all three behind one function pointer; they are each < 10 flops per cell.

2. Threshold: Soulsby-Whitehouse theta_cr(D*) with the table above; multiply by the Soulsby
   slope factor (eq. 7a) with phi_i = 32 deg. Clamp theta_cr(slope) >= 0.25*theta_cr_flat
   for numerical safety on near-repose slopes (avalanching handles the rest).

3. Exner on the 2D height field z_b(x,y) under the voxel domain: face-centered qb with
   upwinding along transport direction, lambda = 0.4, explicit update; convert dz_b to
   voxel fill fractions (partial-cell/height-field hybrid) rather than flipping whole voxels,
   otherwise 2.5-5 cm steps create artificial slope shocks. Magnitude check: qb ~ 1e-4 m2/s
   (0.3 mm sand at theta ~ 1.2) with dx = 0.05 m gives dz_b/dt up to ~ qb/dx/(1-lambda) ~
   3e-3 m/s near the pile: one 5 cm voxel in ~15 s of morphological time — the bed-change
   limiter (< 2-5 mm/step) then sets the morphological time step, not the CFL of the flow.

4. Direction: start qb aligned with local tau_b (from the LES wall stress), add Talmon
   transverse deviation with beta_2 = 0.85. This is what makes scour holes get their correct
   inverted-cone shape.

5. Sand slide: iterative 8-neighbour relaxation, trigger 32 deg, relax to 30 deg, k = 0.5,
   run every morphological step after the Exner update. This is mandatory — without it the
   upstream scour-hole face steepens unphysically and the hole never reaches equilibrium shape.

6. MORFAC: safe default 1 for design iterations with steady current (runs are short anyway);
   allow 5-10 for long tidal-cycle runs, never > 50 with reversing currents; enforce the
   1%-of-depth / 10%-of-voxel bed-change guard per step and log violations.

7. Unit test set: (i) flat bed, uniform theta -> zero dz_b/dt; (ii) MPM vs WP vs E-F curves at
   theta = 0.06/0.1/0.3/1.0 against tabulated Phi; (iii) 1D migrating sediment hump (analytic
   Exner celerity c = (1/(1-lambda)) * dPhi/dtheta * dtheta/dh linearization); (iv) conical
   sand pile relaxing to exactly 30-32 deg; (v) Roulund et al. (2005) pile-scour benchmark
   (0.1 m pile, d50 = 0.26 mm, u_f/u_fc ~ 1.0, S/D ~ 0.9-1.0 equilibrium).

Disagreements to be aware of: MPM's classic coefficient 8 vs Wong-Parker's 3.97 (factor ~2 —
use 3.97); Gaia's default friction angle 40 deg vs Soulsby's 32 deg (use 32); Talmon beta_2
0.85 default vs ~1.6 from calibration studies (expose as tunable); van Rijn's piecewise
theta_cr vs Soulsby-Whitehouse (differences < 15%, pick Soulsby-Whitehouse).

## Sources

- Apsley, D., "Bed-Load Sediment Transport Formulae" (U. Manchester hydraulics notes) — https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.pdf
- Wong, M. & Parker, G. (2006), "Reanalysis and Correction of Bed-Load Relation of Meyer-Peter and Muller Using Their Own Database", J. Hydraul. Eng. 132(11) — https://ascelibrary.org/doi/10.1061/(ASCE)0733-9429(2006)132:11(1159)
- van Rijn, L.C. (1984), "Sediment Transport, Part I: Bed Load Transport", J. Hydraul. Eng. 110(10):1431-1456 — https://ascelibrary.org/doi/10.1061/(ASCE)0733-9429(1984)110:10(1431)
- Parker, G., "CISM lecture notes, Ch. 3: Relations for bedload and suspended load" (Engelund-Fredsoe, Ashida-Michiue, slope effects) — http://hydrolab.illinois.edu/people/parkerg/_private/CourseNotes/CISMnot3.pdf
- Gaia/TELEMAC bedload documentation (MPM & Wong-Parker coefficients, Koch-Flokstra beta = 1.3, Talmon beta_2 = 0.85, Soulsby slope correction) — https://hydro-informatics.com/gaia-bedload
- Chauchat, J. et al. (2026), "sedExnerFoam 2412: a 3D Exner-based sediment transport and morphodynamics model", Geosci. Model Dev. 19:2299 (Exner discretization, avalanche flux, stability) — https://gmd.copernicus.org/articles/19/2299/2026/
- Roulund, A., Sumer, B.M., Fredsoe, J. & Michelsen, J. (2005), "Numerical and experimental investigation of flow and scour around a circular pile", J. Fluid Mech. 534:351-401 (E-F bedload, sand slide 32/30 deg) — https://www.researchgate.net/publication/231897347
- Ranasinghe, R. et al. (2011), "Morphodynamic upscaling with the MORFAC approach: Dependencies and sensitivities", Coastal Engineering 58(8) — https://www.sciencedirect.com/science/article/abs/pii/S0378383911000457
- "Morphodynamic upscaling with the MORFAC approach in tidal conditions: the critical MORFAC" (ICCE) — https://icce-ojs-tamu.tdl.org/icce/index.php/icce/article/download/7269/pdf_702/0
- HEC-RAS 2D Sediment Transport manual, "Critical Thresholds for Transport and Erosion" (Soulsby-Whitehouse curve) — https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.6/model-description/critical-thresholds-for-transport-and-erosion
- Deltares, Delft3D FM D-Morphology User Manual (MORFAC implementation, slope factors) — https://content.oss.deltares.nl/delft3dfm2d3d/D-Morphology_User_Manual.pdf
- Soulsby, R. (1997), "Dynamics of Marine Sands", Thomas Telford (theta_cr curve, slope-adjusted threshold eq. 80a, phi_i = 32 deg)

## Verification

Adversarial fact-check performed 2026-07-07. Every equation was re-derived numerically (Python,
double precision) and checked against sources located independently of the ones cited above
(original Wong & Parker paper text, COHERENS sediment-transport manual, Delft3D-FLOW user manual,
TELEMAC/Sisyphe 6.3 reference manual, Baykal et al. 2015 full text, HEC-RAS manuals, Parker CISM
notes, sedExnerFoam GMD paper). Verdict per equation:

1. **Dimensionless groups (Phi, theta, D*) — CONFIRMED.**
   Definitions match Apsley's bedload sheet and COHERENS Ch. 7 (eqs for d* and theta). Recomputed
   constants: (s-1)g = 15.552 m/s2 (note: 15.55 OK); (rho_s-rho)g = 15941.25 N/m3 (OK);
   D* coefficient ((s-1)g/nu^2)^(1/3) = 24162 m^-1 (note's 24160 is 0.01% off — fine). Table
   values D* = 2.42/7.25/24.2/72.5/241.6 for 0.1/0.3/1/3/10 mm all reproduce.

2. **Soulsby-Whitehouse critical Shields — CONFIRMED.**
   theta_cr = 0.30/(1+1.2 D*) + 0.055[1-exp(-0.020 D*)] verified character-for-character against
   the HEC-RAS 2D Sediment manual ("Critical Thresholds for Transport and Erosion") and COHERENS
   Ch. 7 eq. (7.32). Recomputed table: theta_cr = 0.0795/0.0384/0.0311/0.0455/0.0556 and
   tau_cr = 0.127/0.183/0.495/2.176/8.862 N/m2, u*_cr = 0.011/0.013/0.022/0.046/0.093 m/s —
   all note values correct to quoted precision.

3. **MPM / Wong-Parker — CONFIRMED (against the original paper text).**
   Wong & Parker (2006) eq. (22): q* = 4.93(tau*_b - 0.0470)^1.60 and eq. (24):
   q* = 3.97(tau*_b - 0.0495)^1.50 — exactly as stated (PDF of the JHE paper, EPFL mirror).
   Original MPM Phi = 8(theta-0.047)^1.5 confirmed via HEC-RAS 1D/2D manuals. Database check
   from the WP paper: ETH-up uniform 5.21 & 28.65 mm; ETH-nup Dm 0.38-5.21 mm with R = 0.25-3.22
   (s = 1.25-4.22) — the note's "0.4-29 mm, s = 1.25-4.2" is accurate. WP report the corrected
   fits are LOWER than classic MPM by a factor 2.0-2.5 (note said ~2x — OK), and most plane-bed
   runs have excess Shields stress 0.02-0.10 (so theta' up to ~0.15-0.25; the note's
   "extrapolate with caution above ~0.25" is consistent).

4. **van Rijn (1984) bedload — CONFIRMED.**
   Phi = 0.053 T^2.1 / D*^0.3 confirmed in COHERENS eq. (7.60) (validity 200-2000 um) and the
   Gaia/TELEMAC bedload doc. The high-transport branch and skin-friction closure confirmed in the
   Delft3D-FLOW User Manual eq. (11.247): Sb = 0.1 sqrt(g' D50^3) D*^-0.3 T^1.5 for T >= 3.0,
   with grain-related Chezy Cg,90 = 18 log10(12h/(3 D90)) (eq. 11.252) — identical to the note.
   van Rijn's 5-segment theta_cr(D*) curve (0.24/D*; 0.14 D*^-0.64; 0.04 D*^-0.1; 0.013 D*^0.29;
   0.056) verified against Apsley's sheet.

5. **Engelund-Fredsoe (1976) — CONFIRMED, with one usage caution.**
   Full probabilistic form verified against COHERENS eqs. (7.57)-(7.59): Phi = 5p(sqrt(theta) -
   0.7 sqrt(theta_cr)), p = [1 + (pi*f_d/6 / (theta-theta_cr))^4]^(-1/4), f_d (= mu_d) = 0.51;
   validated on 190/270/930 um sands. mu_d = 0.51 and theta_c0 = 0.05 also confirmed in Baykal
   et al. (2015, Phil. Trans. R. Soc. A 373:20140104, full text PMC4275922). Closed form
   Phi = 18.74(theta-theta_cr)(sqrt(theta)-0.7 sqrt(theta_cr)) confirmed as Parker CISM notes
   eq. (3.14) (fetched PDF), with theta_cr = 0.05; algebra check: 5/((pi/6)*0.51) = 18.72 ~ 18.74.
   Ashida-Michiue coefficient 17 and 0.3-7 mm validity also confirmed (same source).
   CAUTION (numeric check): the 18.74 closed form is a near-threshold linearization (p << 1).
   It matches the full p-form to <1% at theta <= 0.1 but overshoots by ~15% at theta = 0.3 and
   ~3.6x at theta = 1.0 (15.0 vs 4.21). At our design Shields numbers (theta ~ 0.3-1.2) use the
   full p-form, never the closed form — the note's recommendation already points this way.

6. **Exner equation — CONFIRMED.**
   (1-lambda_s) dz_b/dt + div_H(qb) = D - E verified verbatim against sedExnerFoam paper
   (Renaud, Bonamy, Bertrand & Chauchat, GMD 19:2299-2331, 2026 — citation checks out).
   Porosity 0.4 confirmed as the value used by Baykal et al. (2015) (n = 0.4). The per-step bed
   change limit "< 2-5 mm (5-10% of voxel)" is project implementation guidance, not a literature
   constant (arithmetically consistent: 5-10% of 2.5-5 cm = 1.25-5 mm); sedExnerFoam confirms
   upwind (1st/2nd order) discretization of div(qb) and Euler/Adams-Bashforth time stepping.

7. **Slope corrections — CONFIRMED (two citation footnotes).**
   Threshold formula verified against Sisyphe 6.3 Reference Manual eq. (27):
   theta_c/theta_c0 = [cos psi sin chi + (cos^2 chi tan^2 phi_s - sin^2 psi sin^2 chi)^0.5]
   / tan phi_s — identical to the note's Soulsby form (chi = beta). Sisyphe/Gaia default
   friction angle is 40 deg (as the note warns); 32 deg is the Roulund/DTU value for sand.
   Talmon deviation verified against Sisyphe eqs. (28)+(30): tan alpha = tan delta -
   (1/(beta_2 sqrt(theta))) dZ/dn with beta_2 = 0.85 default (higher ~1.5-1.6 for rivers).
   Koch-Flokstra magnitude factor verified against Sisyphe eq. (26) and COHERENS eq. (7.83):
   Qb = Qb0 (1 - beta dZ/ds), beta = 1.3. Footnotes: (i) Sisyphe and COHERENS cite the
   Koch & Flokstra paper as 1981 (XIX IAHR Congress, New Delhi), not 1980 — both years circulate
   (report 1980, congress 1981); equation unaffected. (ii) Cross-check example in Section 5:
   intermediate value should read Phi = 3.97*(1.156)^1.5 = 4.93 (not 4.96); final qb = 1.0e-4
   m2/s unchanged.

8. **Sand-slide avalanching — CONFIRMED (angles/hysteresis); transfer constants are
   implementation choices.**
   The 32/30 deg hysteresis is verified: Baykal et al. (2015, same DTU group, following Roulund
   et al. 2005) state the sand-slide is "activated at positions where the local bed angle exceeds
   the angle of repose phi_s = 32 deg and de-activated once the local bed angle is reduced to
   30.0 deg". tan 32 = 0.62, tan 30 = 0.58 as stated. The continuous avalanche-flux alternative
   |q_av| = q_av0 [tanh(tan beta) - tanh(tan phi_r)]/(1 - tanh(tan phi_r)) verified verbatim
   against the sedExnerFoam GMD paper (after Duran Vinent et al. 2019). The 8-neighbour Jacobi
   transfer dV = 0.5 k [(z_i - z_j) - Delta_ij tan(phi)] A_cell with k ~ 0.5 and <10-30 sweeps is
   the standard mass-conserving grid implementation used across scour codes, but k and the sweep
   count are tuning parameters, not published constants — treat as implementation guidance.

Verification sources (independent of the note's own citations where possible):
- Wong & Parker (2006) full text (EPFL mirror) — https://documents.epfl.ch/users/b/bl/blanckae/www/HydrauliqueFluviale/Papers/2006_JHE_Wong_Parker.pdf
- COHERENS documentation, Ch. 7 "Sediment transport model" (S-W curve eq. 7.32; E-F eqs. 7.57-7.59 with f_d = 0.51; van Rijn eq. 7.60; Koch-Flokstra eqs. 7.83-7.84) — https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf
- Delft3D-FLOW User Manual, Sec. 11.5.8 Van Rijn (1984) (eqs. 11.247-11.252: both T branches, Cg,90 = 18 log10(12h/3D90)) — https://content.oss.deltares.nl/delft3d/Delft3D-FLOW_User_Manual.pdf
- Sisyphe 6.3 Reference Manual (eqs. 26-30: Koch-Flokstra beta = 1.3, Soulsby slope formula, Talmon beta_2 = 0.85; MPM validity 0.4-29 mm; van Rijn 0.20-2.0 mm) — https://www.opentelemac.org/downloads/MANUALS/SISYPHE/sisyphe_user_manual_en_v6p3.pdf
- Baykal, Sumer, Fuhrman, Jacobsen & Fredsoe (2015), Phil. Trans. R. Soc. A 373:20140104, full text (mu_d = 0.51, theta_c0 = 0.05, n = 0.4, sand slide 32/30 deg) — https://pmc.ncbi.nlm.nih.gov/articles/PMC4275922/
- HEC-RAS 2D Sediment manual: Critical Thresholds (S-W curve) and MPM/Wong-Parker pages — https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.6/model-description/critical-thresholds-for-transport-and-erosion
- Parker, G., CISM notes Ch. 3 (fetched PDF; eq. 3.13 Ashida-Michiue coeff. 17, eq. 3.14 E-F closed form 18.74) — http://hydrolab.illinois.edu/people/parkerg/_private/CourseNotes/CISMnot3.pdf
- Renaud, Bonamy, Bertrand & Chauchat (2026), sedExnerFoam 2412, GMD 19:2299-2331 (Exner form, avalanche flux, upwind schemes) — https://gmd.copernicus.org/articles/19/2299/2026/
- Gaia/TELEMAC bedload doc (MPM f_mpm = 8 / 0.047; Wong-Parker 3.97/0.0495; van Rijn 0.053; Koch-Flokstra 1.3; Talmon 0.85) — https://hydro-informatics.com/gaia-bedload
- Apsley bedload sheet (fetched PDF; Phi/theta/D* definitions, MPM, van Rijn eq. 22 and 5-segment theta_cr curve) — https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.pdf
