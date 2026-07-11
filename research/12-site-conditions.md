# Offshore Wind Site Metocean & Seabed Conditions (North Sea / Baltic)

## Overview

This note fixes the *environmental parameter space* the morphodynamic simulator must cover: water depths, tidal current magnitudes and vertical profiles, wave climate (context only), seabed grain sizes, seawater properties at 5-15 C, monopile and scour-protection dimensions — and the geometric consequence that a 10 x 10 x 5 m voxel domain **cannot** contain a modern monopile plus its protection, so the simulator must target well-chosen sub-problems with amplified inflow.

Headline numbers: North Sea OWF water depths 5-45 m (new farms mostly 20-45 m); depth-averaged peak spring tidal currents 0.5-1.8 m/s (field-verified range across nine UK farms: 0.54-1.77 m/s); seabed d50 from 0.05 mm silty sand to ~20 mm gravel, dominated by 0.2-0.6 mm sand; monopiles for 10-15 MW turbines are 8-11.5 m diameter with rock scour protection 40-60 m across; near-pile flow amplification factor 1.5-2.0 on velocity (up to ~4 on bed shear stress).

## Key equations & results

### 1. Water depths at wind farms
- UK Round 1/2 farms (field survey of 460 monopiles at 9 OWFs): **5-35 m** total range; e.g. Robin Rigg 1-14 m, Greater Gabbard 21-35 m. [Source: WES 10, 2189 (2025)]
- Current-generation North Sea farms (Hornsea, Dogger Bank, Borssele): ~15-45 m; XXL monopiles are marketed for 40-120 m depth. [Guide to an Offshore Wind Farm; ScienceDirect scour review 2024]
- Baltic (Arkona basin etc.): ~25-45 m, essentially **microtidal** (tide of a few cm to 20 cm; currents wind-driven and typically < 0.5 m/s). [Power Technology / Baltic siting studies]
- Simulator note: our 5 m tall domain is a *near-bed slab*, not the full water column — see Section "Practical guidance".

### 2. Depth-averaged tidal currents (the primary forcing)
- Nine UK OWFs, 99th-percentile depth-averaged current U99: **0.54 m/s (Gunfleet Sands) to 1.77 m/s (Robin Rigg)**. [WES 10, 2189 (2025), Table 1]
- Direct measurement near a 4.7 m monopile in Liverpool Bay (12.5-21 m depth): depth-averaged **spring peak 0.8 m/s flood / 0.6 m/s ebb; neap +-0.35 m/s**. [Ocean Science 21, 81 (2025)]
- Southern North Sea typical peak spring: 0.6-1.2 m/s; adding storm surge, design extreme depth-averaged currents reach ~1.5-2.0 m/s at energetic sites.
- Conclusion: the simulator's 0.5-2.5 m/s inflow range covers ambient conditions (0.5-1.8 m/s) *and* near-pile amplified flow (see eq. 7); 2.0-2.5 m/s should be interpreted as amplified/extreme, not typical ambient.

### 3. Log-law velocity profile (use for inflow BC)
U(z) = (u_star / kappa) * ln(z / z0), kappa = 0.40
with z0 = k_s / 30 and Nikuradse grain roughness k_s = 2.5 * d50, i.e. **z0 = d50 / 12** for a flat immobile bed. Valid from a few cm above the bed up to ~0.2-0.3 of the boundary-layer/water depth. [Soulsby, Dynamics of Marine Sands, 1997]
Soulsby's field-calibrated z0 by bed type (Table 7): silt/sand 0.05 mm; **unrippled sand 0.4 mm**; sand/gravel 0.3 mm; **gravel 3 mm**; **rippled sand 6 mm** (rippled beds are the common state for fine/medium sand under tides). Units: mm.

### 4. 1/7-power-law profile (cheap alternative for inflow BC)
U(z) = U_bar * (z / (0.32 * h))^(1/7) for 0 < z <= 0.5 h; U(z) = 1.07 * U_bar for 0.5 h < z <= h,
where U_bar = depth-averaged speed, h = water depth. Matches the log profile within a few percent for typical sandy North Sea beds. [Soulsby 1997, eq. 26-27]
Example (h = 30 m, U_bar = 1.0 m/s): U(0.5 m) = 0.66 m/s, U(1 m) = 0.72 m/s, U(2.5 m) = 0.83 m/s, U(5 m) = 0.91 m/s — i.e. across our 5 m domain height the inflow rises only ~0.66 -> 0.91 m/s. (For shallow sites, h = 12 m: 0.75 / 0.83 / 0.94 / 1.04 m/s at the same heights.) A mild shear profile, not a full boundary layer, is the correct inflow.

### 5. Depth-averaged drag coefficient / bed shear stress (for calibration)
tau_0 = rho * C_D * U_bar^2, with C_D = [0.40 / (ln(h / z0) - 1)]^2; Soulsby's default **C_D = 0.0025** (canonical shelf-sea value). Note: unrippled-sand z0 = 0.4 mm with h = 10-40 m gives C_D = 0.0014-0.0019; **rippled-sand z0 = 6 mm gives C_D = 0.0026-0.0039** — the 0.0025 default effectively assumes a rippled sandy bed in ~30-40 m of water, not a flat unrippled one. [Soulsby 1997]
Example: U_bar = 1.0 m/s, rho = 1027 kg/m3 -> tau_0 ~ 2.6 Pa, u_star ~ 0.050 m/s. The LES + wall model should reproduce this to within ~20% on a flat sand bed — use as a validation target.

### 6. Seawater properties at 5-15 C (S = 35 g/kg, North Sea)
- 5 C: rho = 1027.7 kg/m3, nu = 1.56e-6 m2/s
- 10 C: rho = 1026.9 kg/m3, nu = 1.35e-6 m2/s
- 15 C: rho = 1025.9 kg/m3, nu = 1.19e-6 m2/s
[ITTC 7.5-02-01-03 salt-water (3.5%) tables; MIT Seawater Property Tables 2017]
**Important correction to the project brief:** nu = 1.05e-6 m2/s corresponds to ~19-20 C water; at the stated 5-15 C use **1.2-1.6e-6 m2/s** (recommend default nu = 1.35e-6 m2/s, rho = 1027 kg/m3 at 10 C). This matters for settling velocity and grain Reynolds numbers of the 0.1-0.3 mm fractions (w_s changes ~10-20%).
- Baltic brackish water (S ~ 7-8 g/kg): rho ~ 1005-1006 kg/m3 at 10 C, nu ~ 1.31e-6 m2/s; raises submerged specific gravity s-1 from 1.58 to 1.64 — a ~4% effect on Shields numbers.

### 7. Flow amplification around a monopile (inflow scaling for near-pile sub-domains)
Potential-flow speed along the lateral (90 deg) axis of a cylinder radius R:
u(r) / U_inf = 1 + R^2 / r^2 -> **2.0 at the pile wall (r = R), 1.5 at r = 1.41 R, 1.25 at r = 2 R**.
Measured/computed near-bed values are slightly lower (~1.5-1.7) because of the boundary layer. Bed shear stress amplification alpha = tau / tau_inf (scales ~ velocity squared): **O(4) at the side edges under the horseshoe vortex in waves; locally up to O(10) at ~45 deg from the upstream stagnation line in steady current** (Hjorth-type measurements). A front-face downflow of up to ~0.4 * U_inf feeds the horseshoe vortex. [Sumer, Christiansen & Fredsoe, J. Fluid Mech. 332 (1997); Sumer & Fredsoe, The Mechanics of Scour in the Marine Environment; Roulund et al. 2005]
Practical inflow rule: to represent the protection annulus beside a pile, multiply the ambient current by **1.5-2.0** and raise inflow turbulence intensity to ~10-20% (ambient marine BL is ~5-10%).

### 8. Monopile & scour-protection dimensions (why they don't fit in 10 m)
- 15 MW turbine: monopile ~**10 m diameter** (90 m long, 1850 t) at 40 m depth; up to **11.5 m** (120 m, 2700 t) at 60 m; XXL class spans 8-11 m. 10 MW class: ~7.5-9 m. Legacy farms: 4-6 m (Horns Rev 1: 4.0 m monopile / 4.2 m transition piece; Egmond aan Zee: 4.6 m; Burbo Bank: 4.7 m). [Guide to an Offshore Wind Farm B.2.1]
- Rock protection, Horns Rev 1 as-built: two layers to **radius 9.5 m** from pile centre (19 m footprint on a 4.0 m pile / 4.2 m transition piece ~ 4.5-4.75 D); **filter layer d50 = 0.10 m, 0.5 m thick; armour d50 = 0.40 m, 1.0 m thick**; traditional designs extend to ~**6 D** outer diameter. [JMSE 7(12):440, 2019; Nielsen et al., Horns Rev 1 sinking study]
- BOEM estimate for a 12 MW monopile: protection footprint 0.34 ha = 3400 m2 -> equivalent diameter ~66 m. Rule of thumb for current generation: **footprint 40-60 m diameter, thickness 1.5-2.2 m**.

### Wave climate (context for later milestone)
Nine UK OWFs: 99th-percentile Hs = 1.5-2.7 m. Southern North Sea: annual mean Hs ~ 1-1.5 m; most frequent winter sea state 1.5-2 m; 50-year Hs ~ 6.5-9.5 m (southern) rising to nearly double in the northern North Sea; wind-sea Tp typically 4-8 s, swell 8-14 s. Baltic (FINO2) markedly milder than North Sea (FINO1/3). [WES 10, 2189 (2025); Fraunhofer Windmonitor; Weisse & Gunther 2007]

### Seabed sediments at wind farms
- Nine UK OWFs: d50 from **0.052 mm (cohesive/silty) to ~19.9 mm (medium gravel)**; per-farm examples: Robin Rigg fine-medium sand 0.182-0.268 mm; Lynn & Inner Dowsing coarse sand 0.695-1.951 mm; Teesside fine/silty sand; Humber Gateway sandy gravel with boulders. [WES 10, 2189 (2025)]
- Liverpool Bay (measured): rippled sand, d50 = 0.25 mm, total roughness height k_b = 0.122 m including ripples. [Ocean Science 21, 81 (2025)]
- Irish Sea ORE sites (36 grabs): d50 from 0.25 mm to ~13.6 mm; 44% gravelly sand, 36% sandy gravel. [Frontiers Mar. Sci. 10:1156486 (2023)]
- Defaults: quartz rho_s = 2650 kg/m3; porosity ~0.4. **Recommended canonical test sands: d50 = 0.20 mm (fine sand, rippled), 0.35 mm (medium sand — most common), 1.0 mm (coarse sand), 5 mm (fine gravel).** The project's 0.1-10 mm range is exactly right.

## Practical guidance for our simulator

**1. Parameter matrix.** Ambient depth-averaged current 0.5 / 0.8 / 1.2 / 1.8 m/s (neap, mean, spring, extreme); near-pile amplified cases 1.5-2.5 m/s. Water: rho = 1027 kg/m3, nu = 1.35e-6 m2/s (10 C default; sweep 1.19-1.56e-6 for 15-5 C). Sand: d50 = 0.2 / 0.35 / 1.0 / 5.0 mm, rho_s = 2650 kg/m3.

**2. Inflow BC.** Impose the log-law (eq. 3) over the 5 m domain height with u_star chosen so the *local* speed matches the target (remember eq. 4: the 5 m slab only sees ~0.66-0.91x the depth-average for h = 30 m, or ~0.75-1.04x for h = 12 m). Cheaper: the 1/7 law referenced to the true depth h (a domain parameter, 10-45 m), clipped at z = 5 m. z0 from Soulsby's table (0.4 mm unrippled sand default; 6 mm if ripples are unresolved subgrid).

**3. Domain reality check.** A 15 MW monopile (D ~ 10 m) alone fills the whole 10 m domain width; its protection (40-60 m) needs a >= 100-150 m domain and >= 10-100x more cells. Also the pile wake sheds vortices at f = St * U / D with St ~ 0.2 -> period ~50 s and wavelength >> 10 m: the full pile problem is out of scope by construction. Valid sub-problems:
   - **(a) Single printed unit or small cluster on open seabed** (unit footprint 1-3 m, height 0.5-1.5 m): 3-4 m upstream fetch, blockage < 5% of the 10 x 5 m cross-section. Primary design-iteration case.
   - **(b) Periodic strip/carpet of units**: spanwise-periodic lateral BCs, unit pitch 2-5 m; measures sheltering, inter-unit scour, and sand-trapping efficiency of a repeated tile.
   - **(c) Sector of the protection annulus near the pile wall**: represent the pile as a vertical (slightly curved, R = 4-5.75 m) wall along one domain face; inflow = ambient current x 1.5-2.0 (eq. 7), directed tangentially, turbulence intensity 10-20%; optionally superimpose a downflow component up to 0.4 * U_inf near the wall face to mimic the horseshoe-vortex feed. This is the case that tests whether printed units survive where rock protections see their worst loads (edge scour at the side edges, alpha ~ 4).

**4. What the near-pile amplification means for loading.** Amplification 2.0 on velocity = factor ~4 on bed shear = far above the Shields threshold for all sands: at 2 m/s local flow over 0.35 mm sand, tau_0 ~ 10 Pa vs tau_cr ~ 0.2 Pa — live-bed conditions, mobility ~50. Printed shapes must work in live-bed transport, not just clear-water; validation cases should include both regimes (ambient 0.5 m/s over 1 mm sand is near-threshold; use for clear-water tests).

**5. Baltic variant.** Same geometry, but currents <= 0.5 m/s, no tidal reversal (use uni-directional or slowly-varying wind-driven flow), rho = 1005 kg/m3. This is the low-energy end of the matrix and mainly tests self-ballasting by wave-driven transport later.

## Sources

- Scour variability across offshore wind farms (460 monopiles, 9 UK OWFs), Wind Energy Science 10, 2189 (2025) — https://wes.copernicus.org/articles/10/2189/2025/
- Enhanced bed shear stress and mixing in the tidal wake of an offshore wind turbine monopile, Ocean Science 21, 81 (2025) — https://os.copernicus.org/articles/21/81/2025/
- Soulsby, R.L., "Dynamics of Marine Sands", Thomas Telford, 1997 — https://eprints.hrwallingford.com/412/
- Sumer, B.M., Christiansen, N., Fredsoe, J., "The horseshoe vortex and vortex shedding around a vertical wall-mounted cylinder exposed to waves", J. Fluid Mech. 332 (1997) — https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/horseshoe-vortex-and-vortex-shedding-around-a-vertical-wallmounted-cylinder-exposed-to-waves/7E3F407EEDA8EDF8DD339721DC4C861C
- Guide to an Offshore Wind Farm, B.2.1 Monopile (BVG Associates / The Crown Estate) — https://guidetoanoffshorewindfarm.com/guide/b-balance-of-plant/b-2-turbine-foundation/b-2-1-monopile/
- BOEM, "Offshore Wind Turbine Foundations" white paper — https://www.boem.gov/sites/default/files/documents/environment/Wind-Turbine-Foundations-White%20Paper-Final-White-Paper.pdf
- "Riprap Scour Protection for Monopiles in Offshore Wind Farms", J. Mar. Sci. Eng. 7(12):440 (2019) — https://www.mdpi.com/2077-1312/7/12/440
- ITTC Recommended Procedures 7.5-02-01-03, Fresh Water and Seawater Properties — https://ittc.info/media/4048/75-02-01-03.pdf
- MIT Seawater Property Tables (2017, r2b) — https://web.mit.edu/seawater/2017_MIT_Seawater_Property_Tables_r2b_2023c.pdf
- Characterizing seabed sediments at contrasting offshore renewable energy sites, Front. Mar. Sci. 10:1156486 (2023) — https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2023.1156486/full
- Fraunhofer IEE Windmonitor — Wave heights and accessibility (FINO1/2/3) — https://windmonitor.iee.fraunhofer.de/windmonitor_en/4_Offshore/3_externe_Bedingungen/3_Wellen/
- Weisse & Gunther, "Wave climate and long-term changes for the Southern North Sea, hindcast 1958-2002", Ocean Dynamics 57 (2007) — https://link.springer.com/article/10.1007/s10236-006-0094-x

## Verification

Adversarial fact-check (2026-07-07) of the eight key equations, against sources located independently of this note. Body text above has been patched where corrections were found (sections 4, 5, 6, 8 and Practical guidance item 2).

### 1. Log-law inflow profile — CONFIRMED
- U(z) = (u_star/kappa) ln(z/z0), kappa = 0.40: standard; kappa = 0.40 confirmed in Soulsby & Clarke, HR Wallingford TR137 (2005), https://eprints.hrwallingford.com/558/1/TR137.pdf, which also states verbatim "Convert d50 to bed roughness length (rough flow) z0 = d50/12" — i.e. z0 = k_s/30 with k_s = 2.5 d50. COHERENS sediment-model documentation (RBINS), eq. 7.6, independently gives z0 = k_b/30 (https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf).
- Field z0 table (silt/sand 0.05 mm, unrippled sand 0.4 mm, rippled sand 6 mm, sand/gravel 0.3 mm, gravel 3 mm): verified indirectly but tightly — the companion C_100 column of Soulsby's table (mud 0.0022, mud/sand 0.0030, silt/sand 0.0016, unrippled sand 0.0026, rippled sand 0.0061, sand/shell & sand/gravel 0.0024, gravel 0.0047), reproduced in "Guidance on Setup, Calibration, and Validation of Hydrodynamic, Wave, and Sediment Models for Shelf Seas and Estuaries" (J. Waterway Port Coast. Ocean Eng., via ResearchGate 321944391), maps 1:1 onto exactly these z0 values through C_100 = [0.40/ln(1 m/z0)]^2 (e.g. z0 = 6 mm -> 0.00611; z0 = 3 mm -> 0.00474; z0 = 0.4 mm -> 0.00261; z0 = 0.05 mm -> 0.00163). Validity to ~0.2-0.3 h is the standard log-layer extent.

### 2. 1/7-power-law profile — CORRECTED (example values only; formula is right)
- Functional form confirmed verbatim in Whitehouse, Lam, Richardson & Keel, OMAE2010-20999 (https://pure.manchester.ac.uk/ws/portalfiles/portal/34052203/FULL_TEXT.PDF), citing Soulsby (1997): U(z) = U_bar (z/(0.32 h))^(1/7) for 0 < z < 0.5h; U(z) = 1.07 U_bar for 0.5h < z < h. Self-consistency checks pass: (0.5/0.32)^(1/7) = 1.066 ~ 1.07, and the profile depth-averages to 1.001 U_bar.
- The worked example was WRONG: for h = 30 m, U_bar = 1 m/s the formula gives U(0.5 m) = 0.66, U(1 m) = 0.72, U(2.5 m) = 0.83, U(5 m) = 0.91 m/s. The quoted 0.75 / 0.83 / 0.94 / 1.04 m/s are the values for h = 12 m, not 30 m. Body text fixed. Downstream consequence: a 5 m near-bed slab at a 30 m site sees ~0.66-0.91x the depth-averaged current (not 0.75-1.05x); guidance item 2 fixed accordingly.

### 3. Depth-averaged drag / bed shear stress — CORRECTED (parenthetical CD range; formula and default are right)
- tau_0 = rho C_D U_bar^2 with C_D = [0.40/(ln(h/z0) - 1)]^2 confirmed twice independently: COHERENS eq. 7.2 gives Cdb = [kappa/(ln(H/z0) - 1)]^2 for the depth-averaged (2-D) case; Inch Cape ES Annex 10A.4 (marine.gov.scot) quotes Soulsby (1997) Eq. 37 as CD = (0.4/(1 + ln(z0/d)))^2, algebraically identical.
- Default C_D = 0.0025 is the canonical shelf-sea value (used e.g. in Ocean Science 21, 81 (2025), and traceable to Dronkers 1964 / Soulsby 1997).
- The claim "h = 10-40 m, z0 = 0.4 mm gives 0.0022-0.0028" is WRONG: direct evaluation gives C_D = 0.0019 (h=10 m) to 0.0014 (h=40 m). Rippled-sand z0 = 6 mm gives 0.0039 (h=10 m) to 0.0026 (h=40 m), which is where the 0.0025 default actually lives. Body text fixed.
- Worked numbers check out: U_bar = 1 m/s, rho = 1027, C_D = 0.0025 -> tau_0 = 2.57 Pa, u_star = sqrt(tau_0/rho) = 0.050 m/s.

### 4. Potential-flow amplification beside a cylinder — CONFIRMED
- Classical solution: u_theta = -U_inf (1 + R^2/r^2) sin(theta); along the lateral axis (theta = 90 deg) the speed is U_inf (1 + R^2/r^2) -> 2.00 at r = R, 1.50 at r = sqrt(2) R = 1.41 R, 1.25 at r = 2R. Verified analytically (any potential-flow text, e.g. Batchelor 1967; Wikipedia "Potential flow around a circular cylinder"). The 2x wall value is the standard textbook maximum. Near-bed measured values 1.5-1.7 are plausible (boundary-layer reduction) but not independently pinned to a specific figure — treat as indicative.

### 5. Bed shear stress amplification around a pile — CONFIRMED
- Baykal, Sumer, Fuhrman, Jacobsen & Fredsoe, Coastal Engineering (2017) (DTU Orbit, https://backend.orbit.dtu.dk/ws/files/129005572/manuscript_baykal_etal_CENG_submitted_r2_afterProof.pdf) state verbatim: "while the maximum amplification in the bed shear stress around the pile remains at O(3-4) in the wave cases, it can be as much as O(10) in the case of the steady current", validated against Sumer et al. (1997) wave data and Hjorth (1975, Fig. 6.18b) / Roulund et al. (2005, Fig. 16) steady-current data. Other literature quotes 7-11 for the steady-current horseshoe-vortex maximum, located ~45 deg off the upstream stagnation line (waves: maximum at the side edges) — matching the note.
- Caveat on the downflow: "up to 0.4 U_inf" is consistent with flat-bed measurements/simulations (Roulund et al. 2005) but literature values range up to ~0.8-0.9 U_inf inside a developed scour hole; keep 0.4 U_inf for the initial flat-bed configuration only.

### 6. Seawater properties — CONFIRMED
- ITTC 7.5-02-01-03 Table 3 (standard seawater, S = 35.16 g/kg; PDF fetched and parsed): 5 C: rho = 1027.72, nu = 1.576e-6; 10 C: rho = 1027.00, nu = 1.360e-6; 15 C: rho = 1026.02, nu = 1.189e-6; 20 C: nu = 1.0508e-6 m2/s. The note's values (1027.7/1.56e-6, 1026.9/1.35e-6, 1025.9/1.19e-6) match within ~1%, the small density offsets being exactly the S = 35.00 vs 35.16 g/kg difference (UNESCO EOS-80 at S = 35: 1027.68 / 1026.95 / 1025.97 kg/m3 — computed directly, agrees).
- The headline correction stands and is now exact: nu = 1.05e-6 m2/s is 20 C seawater (ITTC: 1.0508e-6); at 5-15 C use 1.19-1.58e-6.
- Baltic: EOS-80 gives rho(S=7, 10 C) = 1005.2 kg/m3 — confirms ~1005.

### 7. Vortex-shedding frequency — CONFIRMED (with a Reynolds-number caveat)
- f = St U/D with St ~ 0.2 is the standard circular-cylinder value across the subcritical regime; D = 10 m, U = 1 m/s -> f = 0.02 Hz, period 50 s, advective wavelength ~50 m >> 10 m domain. Conclusion (pile-scale wake cannot be resolved in the 10 m domain; replace with amplified mean inflow + 10-20% turbulence intensity) is robust.
- Caveat: at the actual pile Reynolds number Re = UD/nu ~ 7e6 (supercritical/transcritical), smooth-cylinder St rises to ~0.25-0.3 (Schewe 1983; Roshko 1961), i.e. period 33-40 s; marine-growth-roughened piles revert toward St ~ 0.2. Either way the order-of-magnitude argument is unchanged.

### 8. Scour-protection as-built dimensions — CORRECTED (minor: pile diameter)
- Nielsen et al., "Sinking of scour protections at Horns Rev 1 Offshore Wind Farm", Proc. ICCE 2014 (https://icce-ojs-tamu.tdl.org/icce/index.php/icce/article/download/7700/pdf_865), confirms: filter layer stones 0.03/0.10/0.20 m (min/median/max, i.e. d50 = 0.10 m), 0.5 m thick; cover (armour) stones 0.35/0.40/0.55 m (d50 = 0.40 m), two stone layers; protection design diameter "around 6 times the pile diameter (outer limits)"; water depths 6.5-13 m; currents ~0.5 m/s (up to 0.8 m/s in storms). The 9.5 m as-built radius matches JMSE 7(12):440 (2019).
- CORRECTION: the Horns Rev 1 monopile is 4.0 m diameter; 4.2 m is the transition piece ("monopiles with a diameter of 4.0 m (transition piece 4.2 m)"). 19 m footprint = 4.5-4.75 D. Body text fixed.
- BOEM foundations white paper: scour protection for a 12 MW monopile disturbs 0.34 ha (0.85 acres); equivalent diameter sqrt(4*3400/pi) = 65.8 m ~ 66 m — arithmetic checks.

### Verification sources (independent of the note's own citations)
- Soulsby & Clarke (2005), TR137, HR Wallingford — https://eprints.hrwallingford.com/558/1/TR137.pdf
- COHERENS model documentation, Ch. 7 Sediment transport (RBINS) — https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf
- Inch Cape Offshore ES, Annex 10A.4 Bed Shear Stress Analysis — https://marine.gov.scot/datafiles/lot/inch_cape/Environmental%20Statement/Volume%202B%20-%20Appendices/Annex%2010A.4%20-%20Bed%20Sheer%20Stress%20Analysis%20Methodology.pdf
- Whitehouse, Lam, Richardson & Keel (2010), OMAE2010-20999 — https://pure.manchester.ac.uk/ws/portalfiles/portal/34052203/FULL_TEXT.PDF
- Larson et al., Guidance on Setup, Calibration, and Validation of Hydrodynamic, Wave, and Sediment Models — https://www.researchgate.net/publication/321944391
- ITTC 7.5-02-01-03 Fresh Water and Seawater Properties — https://ittc.info/media/4048/75-02-01-03.pdf (Table 3, parsed directly)
- UNESCO EOS-80 equation of state (direct computation for S = 35 and S = 7)
- Baykal, Sumer, Fuhrman, Jacobsen & Fredsoe (2017), Coastal Engineering — https://backend.orbit.dtu.dk/ws/files/129005572/manuscript_baykal_etal_CENG_submitted_r2_afterProof.pdf
- Nielsen et al. (2014), Sinking of scour protections at Horns Rev 1, Proc. ICCE — https://icce-ojs-tamu.tdl.org/icce/index.php/icce/article/download/7700/pdf_865
- BOEM, Offshore Wind Turbine Foundations white paper (0.34 ha / 0.85 acre figure) — https://www.boem.gov/sites/default/files/documents/environment/Wind-Turbine-Foundations-White%20Paper-Final-White-Paper.pdf
