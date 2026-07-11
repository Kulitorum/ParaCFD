# Scour Protection State of the Art & Failure Modes

## Overview

Virtually all installed offshore-wind scour protection is rock armour: a granular filter layer (~0.5 m) plus 2 stone layers of quarry-run armour, laid in a circle of 5-6 pile diameters (Dp) around the monopile (Horns Rev 1: outer diameter = 6 Dp; Egmond aan Zee/OWEZ: ~25 m around a 4.6 m pile ≈ 5.4 Dp). Two design philosophies exist: **static** (no armour motion allowed under design load; threshold-of-motion sizing with a bed-shear amplification factor α = 2-4 for the pile) and **dynamic** (limited stone displacement accepted, quantified by the 3D damage number S3D ≤ 1; gives 1.2-1.8x smaller stones — De Vos et al. 2012). The governing standard is DNV-RP-0618 "Rock scour protection for monopiles" (2022-09); DNV-RP-C212 treats scour geotechnically (scour = full loss of lateral/axial soil resistance down to scour depth).

Documented failure modes (all are test metrics our simulator must reproduce):
1. **Shear failure** — armour entrained where pile-amplified shear stress exceeds the stone threshold (amplification up to ~4 in steady current near a cylinder; design values α = 2-4).
2. **Winnowing** — base sand sucked out through armour pores; countered by geometrically closed filters or geotextile.
3. **Sinking** — Horns Rev 1: protections sank 0.5-1.5 m (locally more) between 2002 and 2005 despite intact armour; wave-driven suction of sand through the whole protection (currents there, ≤0.6-0.7 m/s, were too weak). Lab: sinking in waves ≈ 1.3-2.4 x cover-stone size; thicker/wider protections sink MORE (more pore volume to fill before sand infill from the perimeter re-seals the layer) (Nielsen et al., ICCE 2014).
4. **Edge scour** — secondary scour at the armour perimeter. In current: counter-rotating wake vortices scour a downstream hole whose depth/length scale with Dp and the protection thickness-to-width ratio (Petersen et al. 2015, Coast. Eng. 106). In waves: roughness-change streaming removes sand on/offshore of the berm and deposits it INSIDE the porous stone layer; edge scour grows with KC and berm aspect ratio Ar = hb/Wb, and wave values stay below current values for KC < 17. OWEZ field surveys: armour lowered 0.2-0.8 m adjacent to the pile within one year; stones within 1 Dp displaced to a ring at 1-1.5 Dp.

Alternatives to rock (Deltares HaSPro handbook, 2023): articulated concrete block mattresses; frond mats (>1000 buoyant fronds/m2, ~1 m high, mimic seaweed — reduce local velocity, trap sand into a self-formed berm, and suppress edge scour without geotextile); geotextile sand containers (cheap, but sand lost once fabric abrades/UV-degrades); gabions/rock bags. Eco-engineered units are the direct precedent for COBOD: ARC Marine Reef Cubes (interlocking hollow concrete cubes, 3D-printed surface textures) — world-first full-scale eco scour protection at Rampion OWF (~75,000 units, 2025). Lab evidence on hollow perforated units (Du et al., arXiv 2503.13860): cube units with 3x3 side holes (19.6% side porosity, open bottom) and hemisphere units (31.3% surface porosity) around a monopile reduced scour 37-100% depending on velocity (0.20-0.30 m/s, d50 = 0.235 mm, Shields 0.027-0.062); combined layouts gave 70-100% reduction; distinct sand deposition zones formed behind the arrays; hemispheres "sank into the bed and filled the scour pit" (self-healing), whereas taller cubes suffered edge scour deeper than the unprotected pile scour and were displaced by up to half their width. Geometric lessons: open bottoms + internal voids let units self-ballast and self-seal; low outer profile (small hb/Wb) minimizes edge scour; porous perimeters convert edge erosion into deposition inside the structure.

## Key equations & results

**1) De Vos et al. (2012) dynamic damage formula** (Coast. Eng. 60:286-298, Eq. 9; irregular waves + current, 1/50 scale):

S3D / N^b0 = a0 * Um^3 * Tm-1,0^2 / ( sqrt(g*d) * (s-1)^(3/2) * Dn50^2 )
           + a1 * ( a2 + a3 * (Uc/ws)^2 * (Uc + a4*Um)^2 * sqrt(d) / ( g * Dn50^(3/2) ) )

with b0 = 0.243, a0 = 0.00076, a2 = -0.022, a3 = 0.0079;
a1 = 0 if Uc/sqrt(g*Dn50) < 0.92 AND waves follow current, else a1 = 1;
a4 = 1 (waves following current) or a4 = Ur/6.4 (opposing), Ursell Ur = H*L^2/d^3.
N = number of waves, Um = sqrt(2)*sigma_U from the orbital-velocity spectrum, Tm-1,0 = m-1/m0 spectral period (Tp = 1.107*Tm-1,0 for JONSWAP gamma = 3.3), Uc = depth-averaged current [m/s], d = water depth [m], s = rho_s/rho_w, Dn50 = nominal stone size [m].

**2) Damage number definition & criterion** (De Vos et al. 2012, Eqs. 3-4): S3D,sub = Ve / (Dn50 * pi*Dp^2/4) per subarea; S3D = max(S3D,sub). Ve = eroded volume. S3D = 1 means one armour-layer height (1 Dn50) lost over the subarea. Accept S3D ≤ 1 for an armour layer ≥ 2.5 Dn50 thick; no dynamically stable profile develops for larger damage.

**3) Auxiliary relations** (De Vos et al. 2012): Dn50 = 0.84 * D50; stone fall velocity ws = 1.1 * sqrt( (s-1) * g * D50 ) for D50 ≥ 1 mm. Worked example (d = 20 m, Hm0 = 6.5 m, Tp = 11.2 s, Uc = 1.5 m/s): static D50 = 0.496 m vs dynamic D50 = 0.27-0.46 m for S3D = 1-0.2.

**4) Static sizing via threshold of motion** (Soulsby 1997, "Dynamics of Marine Sands"): depth-averaged threshold current U_cr = 7 * (h/D50)^(1/7) * sqrt( g*(s-1)*D50*theta_cr ), theta_cr ≈ 0.055-0.056 for stones (D* large); Soulsby-Whitehouse general form theta_cr = 0.30/(1+1.2*D*) + 0.055*(1-exp(-0.020*D*)), D* = (g*(s-1)/nu^2)^(1/3) * D50. Static design: size D50 so alpha * theta_ambient ≤ theta_cr with pile amplification alpha = 2-4.

**5) Izbash criterion** (Izbash 1935, via CIRIA C683 Rock Manual / USACE CEM): near-bed critical velocity u_b,cr = E * sqrt( 2*g*((rho_s-rho_w)/rho_w)*D ), E = 0.86 (stone exposed on top of layer), E = 1.2 (stone embedded among others). Useful in our solver because the LES resolves local near-bed velocity directly (no alpha needed).

**6) OPTI-PILE stability parameter** (den Boon et al. 2004, in De Vos et al. 2012): Stab = theta_max/theta_cr with theta_cr = 0.056; Stab < 0.4155 = no movement; 0.4155-0.46 = movement without failure; > 0.46 = failure. (Whitehouse et al. 2006 found the limits need recalibration — treat as indicative.)

**7) Filter rules** (Terzaghi, per CIRIA Rock Manual / DNV-RP-0618): geometrically closed if D15,filter / d85,base ≤ 4-5 (retention) and D15,filter / d15,base ≥ 4-5 (permeability); typical OWF filter D50 ≈ 0.10 m (Horns Rev: 0.03/0.10/0.20 m for D15/D50/D85), armour D50 ≈ 0.40 m (0.35-0.55 m), filter 0.5 m thick, armour 2 layers (~0.8-1.1 m).

**8) Edge scour & sinking scalings**: edge scour S/Dp = f(KC, Ar = hb/Wb, Wb/Dp, theta) — increases with KC and Ar; current-generated edge scour (wake vortex pair) sets the upper bound (Petersen et al. 2015; ISSMGE wave study, Wb/Dp = 4, hb/Dp = 0.2-0.5, KC < 17). Sinking in waves = 1.3-2.4 x cover-stone size (≈0.5-1.0 m at Horns Rev scale; up to 1.3 m in current alone) (Nielsen et al. 2014).

## Practical guidance for our simulator

- **Stone/element stability check**: since our LES resolves the amplified local flow, implement Shields-based mobilization per bed/armour voxel using theta = tau_b/((rho_s-rho_w)*g*D) vs theta_cr(D*) from Soulsby-Whitehouse (Eq. 4), and cross-check with Izbash (Eq. 5) on the resolved near-bed velocity. No empirical alpha factor.
- **Concrete numbers** (rho_s = 2650, rho_w = 1025, s = 2.585, nu = 1.05e-6): at Uc = 2.0 m/s in h = 5 m, Soulsby threshold gives stable D50 ≈ 2 cm for the undisturbed bed; with pile amplification alpha = 3-4 on shear stress the classical static requirement rises to D50 ≈ 6-9 cm. Our d50 = 0.1-10 mm seabed is therefore fully mobile at 1-2.5 m/s — deposition INTO the printed shape must come from wake sheltering and streaming, exactly as observed for porous berms and frond mats.
- **Damage metric**: replicate S3D from the Exner bed: Ve over pile-centred subareas of area pi*Dp^2/4, normalized by a characteristic element size; report max. Target S3D ≤ 1 per design storm/tidal event.
- **Failure modes as simulator outputs**: (a) edge scour depth/extent at the shape perimeter (needs wake vortex pair — check LES resolves counter-rotating streamwise vortices at 2.5-5 cm voxels); (b) winnowing/suction — sub-voxel sand transport through unit openings will need a parameterized porous-flux model (Darcy-Forchheimer inside voxelized voids) since 2.5 cm cells cannot resolve pore flow through 30-200 mm filter stones; (c) sinking via undermining (Exner erosion under obstacle footprints must be allowed, letting units settle).
- **Shape-design targets from the evidence**: open bottom + internal cavities (self-sealing sinking, H-AR behaviour); side perforations of 20-30% porosity that admit sand-laden flow but drop internal velocity below the suspension/settling threshold (ws of 0.1-10 mm quartz sand: ~0.008-0.35 m/s — internal velocities must fall below ~0.5-1 x ws-driven Rouse criterion for net deposition); low aspect ratio hb/Wb (< ~0.2) at the perimeter to limit edge scour; interlocking outer geometry (Reef Cube precedent) so units are not displaced by edge scour (cubes moved up to half their width when tall and free-standing).
- **Validation cases**: (1) De Vos worked example (Eq. 1 vs simulated damage); (2) Petersen edge scour in current, Wb/Dp = 4; (3) Du et al. perforated-unit tests (0.20-0.30 m/s, d50 = 0.235 mm) — directly analogous to our printed shapes and fully within our solver's velocity range.

## Sources

- De Vos, De Rouck, Troch, Frigaard (2012), "Empirical design of scour protection around monopile foundations. Part 2: Dynamic approach", Coastal Engineering 60:286-298 — https://www.vliz.be/imisdocs/publications/234400.pdf
- Nielsen, Hansen et al. (2014), "Sinking of scour protections at Horns Rev 1 Offshore Wind Farm", Proc. ICCE 34 — https://icce-ojs-tamu.tdl.org/icce/article/view/7700
- Petersen, Sumer, Fredsoe, Raaijmakers, Schouten (2015), "Edge scour at scour protections around piles in the marine environment — laboratory and field investigation", Coastal Engineering 106:42-72 — https://www.sciencedirect.com/science/article/abs/pii/S0378383915001374
- Petersen et al., "Edge scour at scour protections around monopiles in waves" (ISSMGE) — https://www.issmge.org/uploads/publications/108/137/Edge_scour_at_scour_protections_around_monopiles_in_waves.pdf
- Du et al. (2025), "An experimental study of using artificial reefs as scour protection around an offshore wind monopile", arXiv:2503.13860 — https://arxiv.org/html/2503.13860v1
- Esteban et al. (2019), "Riprap Scour Protection for Monopiles in Offshore Wind Farms", J. Mar. Sci. Eng. 7(12):440 — https://www.mdpi.com/2077-1312/7/12/440
- "A Review and Design Principle of Fixed-Bottom Foundation Scour Protection Schemes for Offshore Wind Energy", J. Mar. Sci. Eng. 12(4):660 (2024) — https://www.mdpi.com/2077-1312/12/4/660
- DNV-RP-0618 "Rock scour protection for monopiles" (2022) — https://www.dnv.com/energy/standards-guidelines/dnv-rp-0618-rock-scour-protection-for-monopiles/
- DNV-RP-C212 "Offshore soil mechanics and geotechnical engineering" — https://www.dnv.com/energy/standards-guidelines/dnv-rp-c212-offshore-soil-mechanics-and-geotechnical-engineering/
- Deltares (2023), "Handbook of Scour and Cable Protection Methods" (JIP HaSPro) — https://publications.deltares.nl/Deltares250.pdf
- ARC Marine RESP / Reef Cubes at Rampion OWF — https://arcmarine.co.uk/resp-reef-enhancing-scour-protection-for-offshore-wind/ ; https://ocean-energyresources.com/2025/10/16/worlds-first-full-scale-use-of-eco-engineered-scour-protection/
- Subsea Protection Systems, anti-scour frond mattresses — https://www.subseaprotectionsystems.co.uk/anti-scour-frond-mats

## Verification

Independent fact-check (2026-07-07) of the eight key equations. Each equation was checked against the full text of the primary paper (extracted from PDF) and/or at least one independent secondary source; every constant and exponent below was read directly from the cited source, not from this note.

**1) De Vos et al. (2012) dynamic damage formula — CONFIRMED.**
Verified against the full text of De Vos et al. 2012, Coastal Eng. 60:286-298 (VLIZ open PDF), Eqs. 9-14 and 21: b0 = 0.243, a0 = 0.00076, a2 = -0.022, a3 = 0.0079 ("respectively equal to 0.243, 0.00076, −0.022 and 0.0079"); a1 = 0 for Uc/sqrt(g·Dn50) < 0.92 AND waves following current, a1 = 1 for Uc/sqrt(g·Dn50) >= 0.92 OR waves opposing current (Eq. 10 — the note's "else" phrasing is logically equivalent); a4 = 1 (following) or Ur/6.4 (opposing) (Eq. 11); Ur = L^2·H/d^3 (Eq. 13, L computed with Hm0 and Tm-1,0); Um = sqrt(2)·sigma_U (Eq. 14); Tp = 1.107·Tm-1,0 for JONSWAP gamma = 3.3 (Eq. 21). Functional form (Um^3·Tm-1,0^2 / (sqrt(g·d)·(s-1)^(3/2)·Dn50^2) wave term; (Uc/ws)^2·(Uc+a4·Um)^2·sqrt(d)/(g·Dn50^(3/2)) current term) matches exactly. Structure independently re-confirmed via Esteban et al. 2019, J. Mar. Sci. Eng. 7(12):440, Table 2. Worked example also confirmed: static D50 = 0.496 m vs dynamic D50 = 0.27-0.46 m (S3D = 1-0.2), reduction factor 1.2-1.8 (Tables 5-6).

**2) 3D damage number & acceptance criterion — CONFIRMED (with one nuance).**
De Vos et al. 2012 Eq. 3: S3D,sub = Ve / (Dn50 · pi·D^2/4) with D = pile diameter; Eq. 4: S3D = max(S3D,sub); "when S3D,sub = 1 ... the height of the scour protection has decreased over this sub area over a distance equal to Dn50" — all verbatim. Nuance: the paper measures the failure transition (damage level 3→4) at S3D = 1.12 and advises S3D <= 1; the tested armour layer was exactly 2.5·Dn50 thick, and the paper states a dynamically stable profile is NOT possible at that thickness (Section 3.3). So read the criterion as "S3D <= 1, calibrated for a 2.5·Dn50-thick layer" — thicker layers may tolerate larger damage numbers (cf. Chambel et al. 2022, ICSE).

**3) Stone auxiliary relations — CONFIRMED.**
Dn50/D50 = 0.84 (De Vos 2012 Eq. 12; independently in Esteban et al. 2019, which adds Wn50/W50 = 0.593). Fall velocity ws = 1.1·[(s-1)·g·D50]^0.5 for D50 >= 1000 μm (De Vos 2012 Eq. 22); this is van Rijn's (1984) fall-velocity expression for coarse grains (1.1 outside the square root).

**4) Soulsby threshold current + Soulsby-Whitehouse theta_cr — CONFIRMED.**
theta_cr = 0.30/(1 + 1.2·D*) + 0.055·[1 − exp(−0.020·D*)], D* = [g(s−1)/nu^2]^(1/3)·d50, quoted verbatim as Soulsby (1997) Eq. 77 in three independent sources: Intertek/Inch Cape Bed Shear Stress Annex 10A.4 (marine.gov.scot), HEC-RAS 2D sediment docs (hec.usace.army.mil), van Rijn "Simple general formulae for sand transport" (leovanrijn-sediment.com). Asymptote for large D* is 0.055 (from the formula); 0.056 is the classical Shields value used by OPTI-PILE — the note's "0.055-0.056" range is correct. The threshold-current formula U_cr = 7·(h/d)^(1/7)·[g·(s−1)·d·theta_cr]^(1/2) is Soulsby (1997) Eq. 72 ("valid for any non-cohesive sediment and water conditions for which D* > 0.1, and valid in any units" — book text); numerically verified against the book's own Example 6.1 (d50 = 200 μm, h = 5 m, 10 °C, 35 ppt, nu = 1.36e-6, D* = 4.06): the formula as written gives U_cr = 0.390 m/s, exactly the book's answer (0.39 m/s); no other leading constant or exponent reproduces this. Pile amplification alpha = 2, 3, 4 on bed shear stress is exactly what De Vos et al. 2012 (Table 5) vary in their static-design comparison — CONFIRMED.

**5) Izbash critical near-bed velocity — CONFIRMED (with attribution/rounding caveats).**
CIRIA C683 Rock Manual (Ch. 5, Eqs. 5.120-5.121) gives "empirically-derived formulae for exposed and embedded stones on a sill" as u_b^2/(2·g·Δ·D50) = 0.7 (exposed) and 1.4 (embedded), i.e. E = sqrt(0.7) ≈ 0.84 and sqrt(1.4) ≈ 1.18 in the form u_b,cr = E·sqrt(2·g·Δ·D50). The note's E = 0.86 / 1.2 are the standard USACE roundings (EM 1110-2-1601 / NRCS TS14C / HEC-RAS riprap calculator: "Isbash constant C = 0.86 for highly turbulent flow, C = 1.20 for low-turbulence flow" — same numbers, turbulence framing). Direction is physically consistent (exposed stone moves at lower velocity). Caveats to carry into the spec: (a) C683 attributes the equations to Izbash & Khaldre (1970), not Izbash 1935; (b) u_b is the velocity NEAR the stone, not depth-averaged; (c) C683 states validity h/D = 5-10. The single-coefficient form u_cr = 1.7·sqrt(Δ·g·d) (Wikipedia/HandWiki) equals the embedded case (1.2·sqrt(2) = 1.70).

**6) OPTI-PILE stability parameter — CONFIRMED.**
De Vos et al. 2012, Eqs. 7-8 and Section 3.2.1 (secondary source for den Boon et al. 2004, which is not openly available): Stab = theta_max/theta_cr with theta_cr = 0.056 and theta_max = tau_max/(rho_w·g·(s−1)·D50); "the limit Stab = 0.4155 is defined as the transition between no movement and movement without failure and the limit Stab = 0.46 as the transition between movement without failure and failure" — verbatim, including the four-decimal 0.4155. Also verbatim: "In Whitehouse et al. (2006) it was noted that for another test series, the limits for the parameter Stab should be adjusted", and De Vos's own data show Stab "fails to correctly predict the observed damage levels" (underestimates for Uc > 0.2 m/s model scale, opposing currents, high-density stones) — treat as indicative only, as the note says.

**7) Terzaghi filter rules — CONFIRMED.**
Retention: D15,filter/d85,base <= 4-5; permeability: D15,filter/d15,base >= 4-5. CIRIA C683 Eq. 5.272 gives the geometrically tight interface criterion D15f/D85b < 5 (valid for uniformly graded materials, D60/D10 < 10); classical Terzaghi & Peck use <= 4; multiple geotechnical references (USBR DSO-04-06; filter-design literature) state both rules with the 4-to-5 band exactly as written. Note C683 also adds internal-stability rules (D10/D5, D20/D10, D30/D15, D40/D20 all < 3; D60/D10 < 10).

**8) Edge scour & sinking scalings — CONFIRMED.**
Edge scour: Petersen et al. ISSMGE wave paper (open PDF, issmge.org) Eq. 10: S/Dp = f(KC, Ar, Wb/Dp, theta) with Ar = hb/Wb — verbatim, including: tested KC < 17, hb/Dp ~ 0.2-0.5, Wb/Dp = 4 (constant); "with increasing KC number the scour depth increases approaching the values experienced for edge scour in current"; scour increases with Ar; Shields-parameter effect weak in live-bed regime. Current-driven edge scour as upper bound — confirmed. Sinking: Nielsen, Sumer & Petersen, ICCE 2014 (open PDF): "the sinking in waves alone will be around 1.3 to 2.4 times the size of the cover stones — corresponding to around 0.5 to 1.0 m sinking. This number will increase up to 1.3 m in the — unlikely — case of current alone" — verbatim match. Also confirmed from the same paper: Horns Rev 1 protections sank up to 1.5 m (2002-2005), currents "almost never exceeded 0.6 to 0.7 m/s", and thicker protections can sink MORE (longer pore-filling time) — all consistent with the note's failure-modes section.

Verification sources (independent of the note's own citation trail where possible):
- De Vos et al. 2012 full text (PDF via VLIZ) — primary source, all constants read from Eqs. 3-22 and Tables 5-6.
- Esteban et al. 2019, JMSE 7(12):440 (open access) — independent reproduction of De Vos formula and Dn50 = 0.84·D50.
- Soulsby 1997 "Dynamics of Marine Sands" book text (Section 6.2, Eqs. 71-72, 77; Example 6.1) via pdfcoffee full-text copy — threshold current speed numerically validated.
- Intertek Inch Cape Annex 10A.4 (marine.gov.scot) and HEC-RAS 2D sediment-transport docs (hec.usace.army.mil) — Soulsby-Whitehouse theta_cr, D*.
- van Rijn, "Simple general formulae for sand transport" (leovanrijn-sediment.com) — theta_cr fit and fall-velocity/threshold context.
- CIRIA C683 Rock Manual Ch. 5 (kennisbank-waterbouw.tudelft.nl) — Izbash Eqs. 5.120-5.121, filter Eq. 5.272.
- USACE/NRCS riprap references (HEC-RAS riprap calculator docs; NRCS TS14C) — Isbash C = 0.86/1.20.
- Petersen et al., "Edge scour at scour protections around monopiles in waves" (issmge.org, open PDF) — Eq. 10 and parameter ranges.
- Nielsen, Sumer & Petersen 2014, ICCE 34 (icce-ojs-tamu.tdl.org, open PDF) — sinking magnitudes and Horns Rev observations.
