# Coupling Hydrodynamics + Sediment + Moving Bed in 3D Morphodynamic Codes

## Overview

Every production 3D scour code uses the same **loosely-coupled (operator-split) loop per morphological step**:

1. **Flow step(s)**: solve momentum + continuity (RANS k-omega in REEF3D and the DTU/OpenFOAM scour solvers; hydrostatic sigma-layer in Delft3D-FLOW; VOF+RANS in FLOW-3D) on the *current* bed geometry.
2. **Bed shear stress** on bed-interface cells via a log-law wall function with Nikuradse roughness ks = 2.5 d50 (FLOW-3D default cs = 2.5).
3. **Slope-corrected critical Shields parameter** (Soulsby 1997 in FLOW-3D, Kovacs & Parker 1994 in REEF3D, Roulund et al. 2005 in the DTU solvers).
4. **Sediment fluxes**: bedload vector qb (van Rijn 1984 / Meyer-Peter-Muller / Engelund-Fredsoe), entrainment E (reference concentration or pickup/lifting velocity) and deposition D = ws * c_b for the suspended phase; advance the suspended advection-diffusion-settling equation.
5. **Exner bed update**: (1-n) dzb/dt = -div(qb) + D - E.
6. **Sand-slide/avalanche sweeps** wherever local slope exceeds the angle of repose.
7. **Geometry update** (mesh motion, level-set, or volume-fraction/mask update) and return to 1.

How the moving bed modifies the fluid domain distinguishes the code families:

- **Body-fitted moving mesh (ALE)**: Roulund et al. (2005), the DTU OpenFOAM chain (Jacobsen; Baykal et al. 2015), sedExnerFoam. Exner is solved on the bed boundary patch (finite-area mesh in sedExnerFoam) and the volume mesh deforms. Accurate wall treatment, but mesh tangling limits slopes and it cannot deposit sand *inside* a complex structure.
- **Level-set on a Cartesian grid**: REEF3D represents both the free surface and the scoured bed as level-set functions; solids use the ghost-cell immersed boundary method (Berthelsen & Faltinsen 2008), extrapolating flow into solid cells. Uniform Cartesian cells (dx = 2 cm for a D = 20 cm pile case).
- **Volume/area fractions (voxel fill)**: FLOW-3D (Wei, Brethour, Grünzner & Burnham 2014) describes the packed bed by FAVOR fractional areas/volumes. **This is the direct published precedent for our approach**: each cell carries a packed-sediment volume fraction; cells whose sediment fraction reaches the *critical packing fraction* (default 0.64, commonly 0.55-0.7) become solid obstacle cells; entrainment ("lifting") converts them back to fluid. If no packed bed exists somewhere initially, a packed-bed component is *created automatically where suspended sediment settles* - i.e., deposition inside/around an arbitrary structure is native to the method. Lattice-Boltzmann partially-saturated-cell methods (Noble & Torczynski 1998) are an analogous precedent for fractional solid cells.
- **Two-phase continuum**: SedFoam-2.0 (Chauchat et al. 2017) solves Eulerian mass/momentum for fluid and particle phases; the "bed" is simply where particle volume fraction alpha reaches max packing (~0.6-0.635), with kinetic-theory or mu(I) granular rheology. No Exner, no reference concentration, no avalanche model needed - but it is far too expensive for 10 m domains (used for sheet flow, apron scour, 2D sluice-gate scour).

**Update frequency**: Baykal et al. run morphology and flow with the *same* variable time step ("morphological and hydrodynamic times are equivalent"), which cost ~10 CPU-days on 8 cores per 1 minute of physical scour. Delft3D instead multiplies the bed flux by a **morphological acceleration factor MORFAC** (fMORFAC, 1-100 for tidal problems) every flow step, with a spin-up delay MorStt (e.g., 720 time units) before any bed update, and a positivity guard: if accretion in one step exceeds the water depth the cell is set dry. REEF3D uses adaptive CFL-based time stepping and decouples the sediment time step from the flow step.

## Key equations & results

All plain-text math; theta = Shields parameter, s = rho_s/rho_w (= 2650/1025 = 2.585 for us), D* = dimensionless grain size, n or p = bed porosity (~0.36 at packing fraction 0.64).

1. **Exner bed evolution** (REEF3D, Ahmad et al.; identical in sedExnerFoam):
   `(1 - n) * dzb/dt = -dqb_x/dx - dqb_y/dy - E + D`
   qb in m2/s, E = entrainment rate (m/s), D = deposition rate (m/s). Valid for any bedload/suspended split; n ~ 0.36.

2. **Dimensionless grain size** (Soulsby 1997; used by all codes):
   `D* = d50 * ((s-1)*g / nu^2)^(1/3)`
   With s = 2.585, g = 9.81 m/s2, nu = 1.05e-6 m2/s: D* = 2.42e4 * d50[m] -> D* = 2.4 (d50 = 0.1 mm), 4.8 (0.2 mm), 24 (1 mm), 242 (10 mm).

3. **Critical Shields parameter, Soulsby-Whitehouse (1997)** (FLOW-3D Eq. 7; valid all D*):
   `theta_cr = 0.30/(1 + 1.2*D*) + 0.055*(1 - exp(-0.02*D*))`
   For us: theta_cr = 0.079 (0.1 mm), 0.049 (0.2 mm), 0.031 (1 mm), 0.056 (10 mm). Slope correction (Soulsby 1997, FLOW-3D Eq. 9): `theta_cr_slope = theta_cr * (cos(psi)*sin(beta) + sqrt(cos^2(beta)*tan^2(phi) - sin^2(psi)*sin^2(beta))) / tan(phi)` with beta = bed slope angle, phi = angle of repose (default 32 deg), psi = angle between flow and up-slope direction.

4. **Van Rijn (1984) bedload** (REEF3D, Ahmad et al. Eq. 14; valid d50 = 0.2-2 mm):
   `qb = 0.053 * sqrt((s-1)*g) * d50^1.5 * T^2.1 / D*^0.3`, with transport-stage parameter `T = (tau_b - tau_cr)/tau_cr`
   qb in m2/s. For d50 > 2 mm switch to **Meyer-Peter-Muller (1948)** (FLOW-3D Eqs. 12-13; valid 0.4-29 mm): `qb = B * (theta - theta_cr_slope)^1.5 * cb * sqrt((s-1)*g*d^3)`, B = 8.0 default (5.0-5.7 low transport, up to 13.0 sheet flow); cb = fraction of that species in the bed. Bedload layer thickness (van Rijn 1984, FLOW-3D Eq. 16): `h_b = 0.3 * d * D*^0.7 * T^0.5`.

5. **Van Rijn reference concentration** (Delft3D-FLOW manual Eq. 11.184, after van Rijn et al. 2000):
   `c_a = min(c_a_max, 0.015 * (d50/a) * Ta^1.5 / D*^0.3)` [volumetric, m3 sediment / m3 water]
   at reference height `a = min( max(AksFac*ks, r/2, 0.01*h), 0.20*h )` (Eq. 11.48; ripple height r = 0.025 m constant, AksFac ~ 1, ks = current roughness). Suspended sediment is exchanged with the bed only through source/sink terms in the lowest cell fully above a (the "kmx layer"), assuming a Rouse profile between a and the kmx cell centre; below a everything is bedload. c_a_max ~ 0.05 by volume.

6. **Entrainment (lifting) velocity, Winterwerp et al. (1992)** (FLOW-3D Eq. 10 - the natural E for a voxel-fill model):
   `u_lift = alpha * D*^0.3 * (theta - theta_cr_slope)^1.5 * sqrt((s-1)*g*d)`, alpha = 0.018 default
   Applied along the outward bed normal; E = u_lift * c_pack converts packed fraction to suspended concentration.

7. **Settling velocity, Soulsby (1997)** (FLOW-3D Eq. 11; valid all D*):
   `ws = (nu/d) * ( sqrt(10.36^2 + 1.049*D*^3) - 10.36 )`
   For us: ws = 0.0073 m/s (0.1 mm), 0.024 m/s (0.2 mm), 0.118 m/s (1 mm), 0.40 m/s (10 mm). Deposition D = ws * c_bottom.

8. **Sand-slide regularization** (Baykal et al. 2015; REEF3D uses Kovacs & Parker 1994): where local bed angle > phi_s = 32 deg, redistribute bed volume to downslope neighbours until the angle falls below a *deactivation* angle of 30.0 deg (hysteresis avoids flip-flopping); iterate sweeps to convergence, conserving volume exactly. Baykal additionally low-pass filters (least-squares smoothing) the bed near the pile to kill grid-scale ripples.

**Validation benchmarks and accuracy actually achieved:**
- Roulund et al. (2005), JFM 534: k-omega + Engelund-Fredsoe bedload + Exner + sand slide, D = 0.1 m pile, d50 = 0.26 mm - the canonical current-scour benchmark; equilibrium S/D and time development reproduced to ~10-20%.
- Baykal et al. (2015): D = 4 cm pile, d50 = 0.17 mm, theta = 0.13, h/D = 2, U = 0.413 m/s -> computed S/D = 0.91, in good agreement with Sumer & Fredsoe (2002) data for that depth ratio. **Omitting suspended load halved upstream scour (S/D 0.91 -> ~0.46)** - suspended load is mandatory for fine sand.
- REEF3D (Ahmad et al.): Link (2006) current scour, D = 0.20 m, d50 = 0.97 mm, U = 0.30 m/s, dx = 2 cm uniform: simulated S = 14 cm vs ~14 cm measured (vs 18 cm from Sumer's interpolated curve, i.e. ~20% low); wave scour KC = 4: S = 1.20 cm vs 1.20 cm measured.
- FLOW-3D (Wei et al. 2014): Chatterjee (1994) submerged-jet scour, d50 = 0.76 mm, theta_cr = 0.05, alpha = 0.018, B = 8.0: scour-hole depth slightly under-, dune height slightly over-predicted; Gladstone (1998) bidisperse settling tank reproduced well with all-default constants.

## Practical guidance for our simulator

**Adopt the FLOW-3D architecture, not ALE**: per-cell packed-sediment volume fraction f_pack in [0,1] on the voxel grid; a cell is solid for the fluid solver when f_pack >= f_crit = 0.64 (porosity 0.36). This natively handles sand settling *inside* our printed shapes (self-ballasting) and is the published precedent to cite (Flow Science Report 03-14). Combine with the STL obstacle mask: obstacle fraction + packed fraction + fluid fraction = 1 per cell.

**Per-step sequence** (steps 2-7 on the *morphological* step): flow substeps -> time-average bed shear over the morphology interval -> theta, theta_cr with slope correction -> qb (van Rijn 1984 for d50 <= 2 mm, MPM B = 8.0 above), E (Winterwerp, alpha = 0.018), D = ws*c -> suspended AD step -> Exner/volume-fraction update -> avalanche sweeps (32/30 deg) -> rebuild solid mask + recompute bed normals. Averaging bed shear over 10-100 flow steps (0.5-5 s) before each bed update both filters turbulent noise and gives a 10-100x morphology speedup; Delft3D's MORFAC (keep <= 10 for local scour, and never let one-step accretion exceed the local water depth) is the alternative knob.

**Stability**: enforce a bed-update CFL: per-step bed change |dzb| <= 0.05-0.1 * dx (and <= 10% of local water depth, the criterion from coupled breach modeling); CFL <= 0.9 on the bed-wave celerity c_bed ~ 3*qb/((1-n)*h). The avalanche sweep is *required* to regularize the Exner equation - without it slopes steepen past repose and grid-scale sawtooth instabilities grow (Baykal; REEF3D both report this).

**Mass bookkeeping** (three reservoirs, audit every N steps): bed solid volume = sum(f_pack * dV_cell) * (1 - internal porosity convention - pick ONE: store f_pack as *solid* volume fraction, max 0.64); suspended volume = sum(c * dV_fluid); bedload has no storage (thin-layer) - its divergence goes straight into f_pack. Exchanges: E (bed->suspended), D (suspended->bed), div(qb) (bed<->bed), plus logged boundary in/out fluxes. Conservative (flux-form) discretization of qb and c keeps global drift at round-off.

**Concrete defaults for our range** (rho_s = 2650, rho_w = 1025 kg/m3, nu = 1.05e-6 m2/s): s = 2.585; f_crit = 0.64; phi = 32 deg (deactivate 30 deg); ks = 2.5*d50; reference height a = first fluid cell centre above the bed (0.5*dx = 1.25-2.5 cm satisfies van Rijn's a >= max(0.5 ks, 0.01h) for h ~ 3-5 m); c_a capped at 0.05; theta_cr, ws, D* tables in item 2/3/7 above. Suspended load matters strongly for d50 <= 0.3 mm (Baykal's factor-2) and is negligible for our gravel end (ws = 0.4 m/s >> u* ~ 0.05-0.1 m/s) - keep both paths but expect bedload-only behaviour above ~2 mm. Target validation: Roulund/Baykal pile case (S/D ~ 0.9-1.0 within +-20%) and Link (2006) at D = 0.20 m, d50 = 0.97 mm - the closest analogue to our geometry scale and grain size.

## Sources

- Delft3D-FLOW User Manual (Deltares, ch. 11 Sediment transport and morphology): https://content.oss.deltares.nl/delft3d4/Delft3D-FLOW_User_Manual.pdf
- Wei, Brethour, Grünzner & Burnham (2014), "The Sedimentation Scour Model in FLOW-3D", Flow Science Report 03-14: https://flow3d.co.kr/wp-content/uploads/FSR-03-14_sedimentation-scour-model.pdf
- Chauchat, Cheng, Nagel, Bonamy & Hsu (2017), "SedFoam-2.0: a 3-D two-phase flow numerical model for sediment transport", Geosci. Model Dev. 10, 4367-4392: https://gmd.copernicus.org/articles/10/4367/2017/
- SedFoam/sedExnerFoam (OpenFOAM single-phase + Exner ALE solver), GitHub: https://github.com/SedFoam/sedExnerFoam
- Baykal, Sumer, Fuhrman, Jacobsen & Fredsoe (2015), "Numerical investigation of flow and scour around a vertical circular cylinder", Phil. Trans. R. Soc. A 373: https://pmc.ncbi.nlm.nih.gov/articles/PMC4275922/
- Ahmad, Bihs, Kamath & Arntsen, "3D numerical modelling of pile scour with free surface profile under waves and current using the level set method in model REEF3D": https://eprints.hrwallingford.com/1014/1/PA_3_20-Ahamd-N.pdf
- Roulund, Sumer, Fredsoe & Michelsen (2005), "Numerical and experimental investigation of flow and scour around a circular pile", J. Fluid Mech. 534, 351-401: https://orbit.dtu.dk/en/publications/numerical-and-experimental-investigation-of-flow-and-scour-around-2/
- Noble & Torczynski (1998), "A Lattice-Boltzmann Method for Partially Saturated Computational Cells", Int. J. Mod. Phys. C 9: https://www.worldscientific.com/doi/10.1142/S0129183198001084
- Soulsby (1997), Dynamics of Marine Sands, Thomas Telford (settling velocity, threshold, slope correction - constants as reproduced in FSR-03-14).

## Verification

Adversarial fact-check (2026-07-07) of the eight key equations against sources retrieved independently of this note's citation chain: USACE HEC-RAS 2D Sediment Transport Technical Reference; Apsley's bed-load formula reference sheet (Univ. of Manchester); the full Delft3D-FLOW User Manual PDF (760 pp., equations located and extracted directly) plus the Delft3D morphology Fortran source on Deltares SVN; van Rijn & Walstra (2003), WL|Delft Hydraulics report Z3624; COHERENS model documentation ch. 7; Flow Science Report FSR-03-14 (retrieved and parsed in full); Wong & Parker (2006) MPM reanalysis; sedExnerFoam GMD 19, 2299 (2026); an OpenFOAM scour-solver paper (arXiv:2012.03051) for the Roulund sand-slide details; and a Saint-Venant-Exner analysis paper (arXiv:2102.00852) for bed-wave celerity. All derived numbers (D*, theta_cr, ws tables) were independently recomputed with s = 2.5854, g = 9.81, nu = 1.05e-6.

1. **Exner bed evolution - CONFIRMED.** sedExnerFoam (GMD 19, 2299, 2026) states exactly `(1 - lambda_s) dzb/dt + div_h(qb) = D - E`, algebraically identical to the note's form including the -E/+D signs; the (1-porosity) factor is the standard Exner form throughout the Saint-Venant-Exner literature. n = 0.36 is consistent with FLOW-3D's default critical packing fraction 0.64 (FSR-03-14, applications sections: "The critical packing fraction is 0.64").

2. **Critical Shields (Soulsby-Whitehouse 1997) - CONFIRMED.** `theta_cr = 0.30/(1 + 1.2*D*) + 0.055*(1 - exp(-0.020*D*))` verified character-for-character in three independent places: HEC-RAS 2D Sediment Technical Reference ("Critical Thresholds for Transport and Erosion"), Apsley's Manchester reference sheet, and FSR-03-14 Eq. 7. D* definition confirmed. Recomputation reproduces the note's table: D* = 2.42/4.83/24.2/241.6 and theta_cr = 0.0795/0.0492/0.0311/0.0556 for d50 = 0.1/0.2/1/10 mm. The slope-correction formula matches FSR-03-14 Eq. 9 verbatim (phi default 32 deg).

3. **Van Rijn (1984) bedload - CORRECTED (incomplete as stated).** The stated formula `qb = 0.053*sqrt((s-1)*g)*d50^1.5*T^2.1/D*^0.3` and the 0.2-2 mm validity range are confirmed (Apsley sheet: Phi = (0.053/D*^0.3)*(tau*/tau*_cr - 1)^2.1, "his equation 22"; COHERENS ch. 7 Eq. 7.60: "particles with diameter in the range of 200 to 2000 um"). **However it is only the T < 3 branch.** The full van Rijn (1984a) relation is piecewise (Delft3D-FLOW User Manual Eq. 11.247): `qb = 0.053*sqrt((s-1)*g*d50^3)*T^2.1/D*^0.3 for T < 3.0; qb = 0.100*sqrt((s-1)*g*d50^3)*T^1.5/D*^0.3 for T >= 3.0`. Our flow range (0.5-2.5 m/s over 0.1-0.2 mm sand) will routinely exceed T = 3, so implement both branches. Caveat: van Rijn defines T with the grain-related (skin-friction) shear stress tau'_b, not total bed shear; wall-resolved CFD codes use the resolved wall shear directly, which is the appropriate analogue.

4. **Meyer-Peter-Muller bedload (FLOW-3D form) - CONFIRMED.** FSR-03-14 Eqs. 12-13 verbatim: `Phi_n = B_n*(theta_n - theta_cr,n)^1.5*c_b,n`, with B "generally 5.0 to 5.7 for low transport, around 8.0 for intermediate transport, and up to 13.0 for very high transport (for example, sand in sheet flow under waves and currents)", default 8.0; c_b,n = volume fraction of species n in the bed (added by Flow Science, not in original MPM; original is q* = 8*(theta - 0.047)^1.5 per Apsley sheet). Bedload-layer thickness confirmed as FSR-03-14 Eq. 16: `h_n = 0.3*d*D*^0.7*(theta/theta_cr - 1)^0.5` (= van Rijn 1984a saltation height). Validity range confirmed: MPM 1948 experiments span d = 0.4-30 mm (Wong & Parker 2006 reanalysis: 0.38-28.65 mm), so "0.4-29 mm" is correct.

5. **Van Rijn reference concentration (Delft3D) - CORRECTED (one constant wrong).** Formula and reference height confirmed exactly: Delft3D-FLOW User Manual Eq. 11.184 `c_a = 0.015*(D50/a)*Ta^1.5/D*^0.3` (volumetric, m3/m3, then multiplied by rho_s), capped by CaMax; Eq. 11.48 `a = min[max(AksFac*ks, dr/2, 0.01h), 0.20h]` with "wave-induced ripple height, set to a constant value of 0.025 m" - equation numbers, form and ripple height all check out (also van Rijn & Walstra 2003 Eqs. 2.2.33/2.2.41; COHERENS Eq. 7.125). **But c_a_max ~ 0.05 is wrong as a Delft3D default: the source code default is CaMax = 0.65** (morphology_data_module.f90, `camax = 0.65_fp`), i.e. essentially uncapped up to max packing. A 0.05 cap is a legitimate conservative modeling choice for us, but do not attribute it to Delft3D. Two further pedantic points: (a) manual Eq. 11.184 misprints the cap as `max(ca,max, ...)`; the code applies `min` (ca,max is defined as the maximum). (b) Ta in Delft3D (Eq. 11.186) uses efficiency factors: Ta = (mu_c*tau_b,cw + mu_w*tau_b,w - tau_cr)/tau_cr, i.e. skin-friction stress, not raw total stress.

6. **Entrainment lifting velocity (Winterwerp 1992) - CONFIRMED.** FSR-03-14 Eq. 10 verbatim: `u_lift,n = n_b * alpha_n * d*_n^0.3 * (theta_n - theta_cr,n)^1.5 * sqrt(g*d_n*(s_n-1))`, alpha default 0.018, n_b = outward bed normal; slope-corrected theta'_cr is used when the sloping-bed option is active. Third-party corroboration of the 0.018/8.0 defaults in published FLOW-3D application studies (e.g. offshore-monopile scour, Ocean Eng. 2024). Note: `E = u_lift * c_pack` is our (sound) volumetric bookkeeping interpretation for the voxel model - FSR-03-14 does not print that exact expression; it converts packed bed to suspended sediment at u_lift.

7. **Settling velocity (Soulsby 1997) - CONFIRMED.** HEC-RAS 2D Sediment Technical Reference ("Particle Settling Velocity") gives `ws = (nu/d)*[(10.36^2 + 1.049*D*^3)^0.5 - 10.36]` - both constants confirmed; identical in FSR-03-14 Eq. 11. Recomputation for our fluids: ws = 0.0073/0.0245/0.117/0.403 m/s for 0.1/0.2/1/10 mm - matches the note's table to rounding.

8. **Bed-update stability + avalanche - CONFIRMED with two precisions.** (a) Avalanche 32/30 hysteresis confirmed: Roulund et al. (2005)-type sand-slide "moves the sediment until the bed slope deceeds the angle of repose by two degrees" (arXiv:2012.03051 describing the Roulund model), i.e. trigger at phi = 32 deg, relax to 30 deg - exactly the note's numbers. (b) Bed celerity: from the quasi-linear Saint-Venant-Exner system with a power-law transport qb ~ u^m, the bed eigenvalue is c_bed ~ m*qb/((1-p)*h*(1-Fr^2)) (arXiv:2102.00852, psi = m*qb/q), so `c_bed ~ 3*qb/((1-n)*h)` is the standard m = 3, Fr << 1 approximation - valid for our tidal flows (Fr ~ 0.1-0.4), degrades near Fr = 1. (c) Precision: Delft3D's per-step bed-change guard default is DzMax = 0.05 (5 % of local water depth, user-adjustable; manual keyword table + `dzmax = 0.05_fp` in source) - use 0.05h, not 0.10h, if citing Delft3D. (d) The per-step |dzb| <= 0.05-0.1*dx cap and CFL <= 0.9 are reasonable engineering heuristics but we found no primary literature source stating them in that form - treat as tunable defaults, not published constraints.

**Verification sources (all retrieved independently):**
- USACE HEC-RAS 2D Sediment Transport Technical Reference: Critical Thresholds (Soulsby-Whitehouse) and Particle Settling Velocity (Soulsby): https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.6/model-description/critical-thresholds-for-transport-and-erosion and .../6.5/model-description/fall-velocity-and-settling/particle-settling-velocity
- Apsley, D., Bed-Load Sediment Transport Formulae reference sheet, Univ. of Manchester: https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.pdf
- Delft3D-FLOW User Manual (Deltares, 760 pp.), Eqs. 11.48, 11.49-11.52, 11.184-11.194, 11.247, mor-file keyword table (DzMax, CaMax): https://content.oss.deltares.nl/delft3d4/Delft3D-FLOW_User_Manual.pdf
- Delft3D morphology source, morphology_data_module.f90 (camax = 0.65_fp, dzmax = 0.05_fp): https://svn.oss.deltares.nl/repos/delft3d/tags/delft3dfm/67136/src/utils_gpl/morphology/packages/morphology_data/src/morphology_data_module.f90
- Van Rijn & Walstra (2003), Modelling of Sand Transport in DELFT3D, WL|Delft Hydraulics Z3624, Eqs. 2.2.33, 2.2.41: https://open.rijkswaterstaat.nl/publish/pages/81574/modelling_of_sand_transport_in_delftd.pdf
- COHERENS documentation ch. 7 (van Rijn 1984 Eqs. 7.60, 7.125 with 200-2000 um validity): https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf
- Wei, Brethour, Grunzner & Burnham (2014), FSR-03-14 (full text parsed; Eqs. 1-2, 7-17): https://flow3d.co.kr/wp-content/uploads/FSR-03-14_sedimentation-scour-model.pdf
- Wong & Parker (2006), Reanalysis and correction of the MPM bed-load relation: http://hydrolab.illinois.edu/people/parkerg/_private/Preprints/MPMWongParkerver.12-05.pdf
- sedExnerFoam 2412 (GMD 19, 2299-xxxx, 2026), Exner equation and avalanche flux: https://gmd.copernicus.org/articles/19/2299/2026/
- OpenFOAM scour solver with Bingham sand-slide (Roulund slide model described: relax to repose minus 2 deg): https://arxiv.org/pdf/2012.03051
- Splitting scheme for coupled Saint-Venant-Exner (bed-perturbation celerity psi = m*qb/q): https://arxiv.org/pdf/2102.00852
