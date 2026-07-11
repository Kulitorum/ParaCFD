# Flow Through Permeable Structures & Sediment Trapping Physics

## Overview

The COBOD scour-protection concept relies on a printed shape that is **porous enough to admit
sediment-laden flow, calm enough inside to let sand deposit, and sheltered enough that the deposit
is not re-eroded at storm-peak flow**. This note collects the quantitative physics for that balance:
(1) settling-basin / sand-trap theory (Camp–Hazen trapping efficiency), (2) flow through porous
screens (velocity reduction vs porosity, with a fully closed-form model), (3) field/lab evidence
from permeable groynes, brushwood/bamboo dams and artificial reefs, and (4) the residence-time vs
settling-time criterion that governs whether a grain entering the structure actually stays.

The governing idea is a **competition of time scales**: a grain deposits if its *settling time*
t_s = H/w_s (H = height of the water column inside the shape it must fall through, w_s = settling
velocity) is shorter than its *residence time* t_r = L/u_in (L = internal flow path length, u_in =
mean interior velocity). Trapping requires t_r >= t_s, i.e. **u_in / L <= w_s / H**, equivalently
**u_in <= w_s · (L/H)**. A porous shell lowers u_in below the ambient current U∞, which is exactly
what makes deposition possible; the shell must also raise the local bed shear-stress deficit enough
that the deposit survives the next flood.

## Key equations & results

### 1. Ideal settling-basin (Hazen 1904 / Camp 1946) trapping efficiency
Treat the interior as a horizontal-flow settling basin with plan area A [m^2] and throughflow
discharge Q [m^3/s]. Define the **overflow / surface-loading rate** (a critical settling velocity):

    v_c = Q / A            [m/s]   (also called surface loading rate)

- Particles with settling velocity **w_s >= v_c are 100% trapped**.
- Particles with **w_s < v_c** are partially trapped with (ideal, quiescent) efficiency
      eta = w_s / v_c = w_s · A / Q      [-],  range 0–1.
- Source: Hazen (1904); Camp (1946). Key result: efficiency depends on **plan area A, not depth** —
  a wide, shallow interior traps as well as a deep one for the same A and Q.

**Camp (1945/46) turbulence correction** (accounts for re-suspension by turbulent eddies, more
realistic for an open structure in a current) — Dobbins–Camp form:

    eta = 1 − exp( − w_s / v_c )   (first-order / well-mixed limit)

or the fuller Camp expression eta = f(w_s/v_c , w_s/u_*) using shear velocity u_*. Use the exponential
form as a conservative estimate: it gives eta ≈ 0.63 when w_s = v_c (vs 1.0 ideal), and eta ≈ 0.39
at w_s = 0.5 v_c. Valid for dilute suspensions (no hindered settling), steady flow.

### 2. Settling velocity of our sediment (Soulsby 1997)
Needed for both v_c comparisons and t_s. Soulsby "Dynamics of Marine Sands" single-grain law:

    w_s = (nu / d) · [ sqrt(10.36^2 + 1.049 · D*^3) − 10.36 ]      [m/s]
    D*  = [ g (s − 1) / nu^2 ]^(1/3) · d       (dimensionless grain size)

with g = 9.81 m/s^2, s = rho_s/rho = 2650/1025 = 2.585, nu = 1.05e-6 m^2/s, d = grain diameter [m].
Computed for our range (seawater, 10 C):

| d50 (mm) | D* | w_s (mm/s) |
|---|---|---|
| 0.10 | 2.42 | 7.3 |
| 0.25 | 6.04 | 33.8 |
| 0.50 | 12.1 | 71 |
| 1.0 | 24.2 | 117 |
| 2.0 | 48.3 | 175 |
| 5.0 | 121 | 283 |
| 10 | 242 | 403 |

Valid 0.1–10 mm (covers 1 <= D* <= ~2e3). Below ~0.1 mm use Stokes w_s = g(s−1)d^2/(18 nu).

### 3. Flow through a porous screen — velocity reduction vs porosity
For a thin permeable barrier normal to the flow, define **porosity beta = A_open / A_gross** (0 = solid,
1 = open). The pressure-loss (resistance) coefficient of the screen (Taylor & Davies 1944):

    k = Δp / (½ rho u^2) = 1/beta^2 − 1        (u = velocity in the pores)

Steiros & Hultmark (2018, JFM 853 R3) give a **closed, self-contained** model that yields both the
drag and, crucially for us, the throughflow (bleed) velocity ratio u/U∞ and the near-wake velocity.
Solve simultaneously (u* = u/U∞):

    C_D = (4/3)·(1 − u*)(2 + u*)/(2 − u*)                     (energy/momentum balance)
    C_D = u*^2 (1/beta^2 − 1) − (4/3)(1 − u*)^3/(2 − u*)^2    (screen resistance)

The far-near-wake mean velocity behind the screen is **U_wake = E · U∞** with E = u*/(2 − u*).
Solving numerically gives our design table:

| beta (porosity) | k = 1/beta^2−1 | bleed u/U∞ | wake U_wake/U∞ | C_D |
|---|---|---|---|---|
| 0.20 | 24.0 | 0.25 | **0.14** | 1.29 |
| 0.23 | 17.9 | 0.28 | **0.17** | 1.27 |
| 0.30 | 10.1 | 0.37 | **0.22** | 1.22 |
| 0.40 | 5.25 | 0.48 | **0.32** | 1.13 |
| 0.50 | 3.00 | 0.59 | **0.42** | 1.00 |
| 0.60 | 1.78 | 0.70 | **0.53** | 0.84 |
| 0.70 | 1.04 | 0.79 | **0.65** | 0.64 |

Validity: model matches data for beta >= 0.23 (shedding suppressed). Solid-plate limit C_D = 4/3 =
1.33 (measured 1.44 ± 0.4). Below **beta ≈ 0.23** a von-Kármán vortex street forms (Castro 1971):
strong wake unsteadiness and re-suspension. **Design implication: a wake interior velocity of
~15–25% of ambient is achievable with beta ≈ 0.2–0.3**, which is the classic optimum band.

### 4. Optimal porosity for shelter (porous-fence / windbreak literature)
Wind-tunnel and field consensus: **maximum leeward velocity reduction and longest sheltered fetch
at porosity beta ≈ 0.2–0.4** (optimum often quoted 0.2–0.35). Above ~0.4 too much flow bleeds
through (weak shelter); below ~0.2 flow separates and a strong recirculation/vortex street forms
that re-suspends the deposit and shortens the sheltered zone. Sheltered length scales with barrier
height h: measurable velocity reduction persists ~10–20 h downstream, >50% reduction out to ~5–8 h
for beta ≈ 0.2–0.3. Sources: Raine & Stevenson (1977); Dong et al. (2007); Cornelis & Gabriels
(2005).

### 5. Emergent-array (vegetation / brushwood) bulk drag and interior velocity
For flow through an array of vertical elements (diameter D, number density n per m^2, solid volume
fraction phi = n·(pi/4)D^2), momentum balance in steady current gives an interior velocity set by
the array drag. Bulk friction (Nepf 1999):

    interior force balance: g·S = (1/2) C_D a U_in^2 / (1 − phi)
    a = n·D  (frontal area per unit volume, [1/m]);  porosity of layer = 1 − phi

Element drag C_D ≈ 1.0 for arrays at rod Reynolds Re_D = U_in·D/nu > ~1000 (friction losses
negligible above Re_D ≈ 1000, Hoerner 1952). Denser arrays (higher a) → lower U_in and stronger
deposition, but also more reflection. Nepf's experiments span phi = 0.01–0.35.

### 6. Wave transmission through permeable (brushwood) dams — Dalrymple/Winterwerp
For the later oscillatory-forcing milestone, wave-height transmission through an emergent porous
dam of thickness L_tot with cylinder array (Dalrymple et al. 1984, as applied by Winterwerp et al.
2020):

    K_t = H_t / H_i = 1 / [ 1 + (4/9)·C_D·D·N·(H_i/L_tot term)·f(kh) ]   (shallow water)

Field result: flexible brushwood bundles gave **up to 80% wave-height reduction (K_t as low as 0.2)**,
on average ~35% reduction (K_t ≈ 0.65). Reflection coefficient K_r ≈ 0.2 at porosity beta ≈ 0.54
(Demak field); K_r = 0.2–0.4 predicted at beta = 0.4 — i.e. porous structures are **low-reflection**,
so they do not scour their own toe the way a solid wall does.

## Practical guidance for our simulator

**Target interior porosity beta ≈ 0.25–0.35 on the flow-facing walls.** This is the sweet spot from
three independent literatures (porous-plate theory §3, wind-fence shelter §4, permeable-dam field
data §5–6): it drops interior mean velocity to ~15–25% of ambient (Steiros–Hultmark E ≈ 0.14–0.22)
while staying above the beta ≈ 0.23 vortex-shedding threshold that would re-suspend the catch.
Fully solid walls (beta→0) maximise wake calm but shed strong vortices, reflect the current (scouring
the toe), and starve the interior of the sediment-laden inflow you want to capture.

**Implement the porous shell as a momentum sink**, not a hard mask. In the Stam solver add a
Darcy–Forchheimer body force in voxels tagged as porous-wall:
`f = −(mu/K)·u − rho·C_F·|u|·u`, or equivalently a thin-screen pressure jump Δp = ½ρ k u^2 with
k = 1/beta^2 − 1 (§3). For beta = 0.3, k = 10.1; for beta = 0.2, k = 24. This reproduces both the
throughflow bleed and the ~80% wake velocity deficit automatically. Calibrate C_D ≈ 1.0 per element
(Re_D > 1000 in our 0.5–2.5 m/s currents with cm-scale printed ribs, so friction losses negligible).

**Trapping criterion to evaluate each candidate shape** (post-process the flow field):
1. Compute interior mean velocity u_in and the vertical fall height H a grain sees inside.
2. Trapping if **u_in <= w_s·(L/H)** (residence >= settling). Using §2 numbers: fine sand d50 =
   0.25 mm has w_s = 34 mm/s and settling time t_s = H/w_s ≈ 15 s for H = 0.5 m; the interior must
   keep u_in low enough that a parcel dwells >15 s. Coarse sand d50 = 1 mm (w_s = 117 mm/s, t_s ≈
   4 s) traps far more easily; d50 = 0.1 mm (w_s = 7 mm/s, t_s ≈ 70 s) is the hard case and needs
   the lowest interior velocities (largest, calmest interior).
3. Cross-check with Camp efficiency: eta = 1 − exp(−w_s·A/Q). For the interior "basin" set Q =
   u_in·A_inlet and A = interior plan area. Aim eta > 0.8, i.e. w_s·A/Q > 1.6, i.e. surface loading
   v_c = Q/A < w_s/1.6. For d50 = 0.25 mm (w_s = 0.034 m/s) that means v_c < 0.021 m/s — the interior
   must present a large horizontal footprint relative to inflow.

**Anti-erosion (self-healing) check at storm flow.** The deposit survives if bed shear stress inside
stays below the critical Shields stress. With interior velocity ~0.15–0.25·U∞ (from §3), bed shear
(∝ velocity^2) is reduced to **~2–6% of the ambient bed shear**. So even at U∞ = 2.5 m/s the interior
"feels" like ~0.4–0.6 m/s — below motion threshold for coarse sand. Verify per shape that interior
tau_b < tau_cr (see the Shields/threshold note) across the full 0.5–2.5 m/s current range.

**Concrete starting geometry for first simulations:** porous walls beta = 0.30 (k = 10), wall
element/rib diameter D = 3–5 cm, interior clear height H ≈ 0.3–0.5 m, interior plan footprint sized
so Q/A < 0.02 m/s at design current. Expect interior velocity ≈ 0.2·U∞, ~80% wave-height damping,
reflection K_r ≈ 0.2–0.3 (low toe-scour risk). Bias the design toward capturing the **coarser end**
(d50 >= 0.25 mm) first — its 5–15 s settling time is comfortably shorter than achievable residence
times; sub-0.1 mm silt needs either much larger/calmer interiors or flocculation and should be
treated as a stretch goal.

## Sources

- Steiros, K. & Hultmark, M. (2018). *Drag on flat plates of arbitrary porosity.* J. Fluid Mech. 853, R3. https://fluids.princeton.edu/pubs/Steiros_2018.pdf — closed-form C_D(beta), bleed and wake velocity, k = 1/beta^2 − 1, shedding threshold beta ≈ 0.23.
- Taylor, G.I. & Davies, R.M. (1944). *The aerodynamics of porous sheets.* ARC R&M 2237 — origin of k = 1/beta^2 − 1.
- Hazen, A. (1904) & Camp, T.R. (1946). Settling-basin overflow-rate theory; Camp turbulence correction. Summary: https://www.particles.org.uk/particle_support/chapter5.htm ; USACE HEC chen-sediment-trap https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/erosion-and-sediment-transport-under-construction/reservoir-sediment-methods/chen-sediment-trap
- Soulsby, R. (1997). *Dynamics of Marine Sands.* Thomas Telford — settling-velocity law, D* definition.
- Winterwerp, J.C. et al. (2020). *Managing erosion of mangrove-mud coasts with permeable dams – lessons learned.* Ecological Engineering. https://www.ecoshape.org/app/uploads/sites/2/2017/08/Managing-erosion-of-mangrove-mud-coasts-with-permeable-dams-%E2%80%93-lessons-learned_Winterwerp-et-al_Ecological-Engineering-22-Sept-20.pdf — dam geometry (poles 0.12–0.15 m dia, 0.6 m spacing, 0.4 m brushwood), K_t up to 0.2, K_r ≈ 0.2 at beta = 0.54, sedimentation 0.15–0.5 m/season, Dalrymple wave-damping formula.
- Uijttewaal, W.S.J. (2005). *Effects of groyne layout on the flow in groyne fields.* J. Hydraulic Eng. 131(9), 782–791. https://www.researchgate.net/publication/245297095 — permeable pile groynes suppress recirculation, reduce mixing-layer turbulence vs solid groynes.
- Nepf, H.M. (1999). *Drag, turbulence, and diffusion in flow through emergent vegetation.* Water Resources Research. https://www.researchgate.net/publication/253447138 — array bulk drag, C_D, phi = 0.01–0.35.
- Cornelis, W.M. & Gabriels, D. (2005); Dong, Z. et al. (2007). Porous-fence optimum porosity 0.2–0.35. https://www.sciencedirect.com/science/article/abs/pii/S0167610522003464
- Zhang, Y. et al. (2025). *Artificial reefs as scour protection around an offshore wind monopile.* arXiv:2503.13860. https://arxiv.org/html/2503.13860v1 — cubic reef beta = 19.6%, hemisphere beta = 31%, wake velocity reduced 50–80%, scour reduced 48–100%, deposition downstream.

## Verification

Independent adversarial fact-check (2026-07-07). Each equation re-derived / re-computed from primary
sources located independently. Primary sources actually opened: the Steiros & Hultmark (2018) JFM PDF
(full text extracted), HEC-RAS sediment documentation for the Soulsby law, the Nazaroff &
Alvarez-Cohen settling-tank lecture notes (Dartmouth ENGS 37) for Camp/Hazen theory, and an
independent Katul-group arXiv paper (2403.12232) restating the Nepf emergent-array balance.

**1. Trapping criterion (residence vs settling) — CONFIRMED.**
t_r = L/u_in, t_s = H/w_s, trap if t_r ≥ t_s ⇒ u_in ≤ w_s·(L/H) is dimensionally consistent and is the
standard residence-time / overflow idealization. It is a worst-case bound (assumes a grain that enters
at the top of the interior must fall the full height H). It is the same physics as the Hazen overflow
criterion in eq 2 (set L/H·A_face geometry). Attribution to Camp (1946) settling-basin theory is
appropriate.

**2. Hazen–Camp trapping efficiency — CONFIRMED.**
Verified against Nazaroff & Alvarez-Cohen settling-tank theory (Dartmouth ENGS 37 notes,
https://cushman.host.dartmouth.edu/courses/engs37/Settling.pdf): overflow rate v_c = Q/A defined on
plan area A = W·L (footprint), NOT cross-section; particles with w_s ≥ v_c are 100% collected; ideal
(quiescent) η = w_s/v_c = w_s·A/Q; efficiency depends on plan area not depth. The turbulent form
η = 1 − exp(−w_s/v_c) is exactly the plug-flow-reactor (transverse-mixing) result and gives
1 − e^(−1) = 0.632 at w_s = v_c and 1 − e^(−0.5) = 0.393 at w_s = 0.5 v_c — both numbers confirmed.
Caveat (already noted in the text): this exponential is the Dobbins/PFR simplification; Camp's full
turbulence expression additionally depends on w_s/u_*.

**3. Soulsby settling velocity — CONFIRMED (formula & constants); CAVEAT on the viscosity/temperature label.**
Formula w_s = (ν/d)·[√(10.36² + 1.049·D*³) − 10.36] with D* = [g(s−1)/ν²]^(1/3)·d verified verbatim
against the HEC-RAS Soulsby (1997) documentation
(https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.5/model-description/fall-velocity-and-settling/particle-settling-velocity):
ω = (ν/d)[(10.36² + 1.049 d*³)^0.5 − 10.36], d* = d(Rg)^(1/3)ν^(−2/3) (identical). Recomputed table with
g = 9.81, s = 2.585, ν = 1.05e-6 reproduces the note exactly: D* = 2.42, 6.04, 12.08, 24.16, 48.32,
120.8, 241.6 and w_s = 7.3, 33.8, 71.1, 117.3, 175.3, 283.4, 402.8 mm/s for d = 0.1–10 mm.
CAVEAT: ν = 1.05e-6 m²/s is NOT seawater at 10 °C. Computed from the Sharqawy/ITTC seawater
correlation, ν(S=35) = 1.576e-6 (5 °C), 1.360e-6 (10 °C), 1.189e-6 (15 °C), 1.051e-6 (20 °C). So
1.05e-6 corresponds to ≈20 °C water; at the project's stated 5–15 °C the correct ν is 1.19–1.58e-6.
Using the true colder-water ν lowers w_s for the fine grains (viscosity-dominated, D* ≲ 10) by roughly
10–20 % (e.g. d50 = 0.1 mm drops from ~7.3 to ~5.9 mm/s at ν = 1.36e-6). Either relabel the table as
"~20 °C" or recompute at the design temperature; the fine-grain (hard) cases are the ones affected.

**4. Porous-screen resistance coefficient — CORRECTED (reference velocity).**
The relation k = Δp/(½ρu²) = 1/β² − 1 and the values (β = 0.2 → 24, 0.3 → 10.11, 0.5 → 3.0) are
correct — confirmed as the full-span limit of Steiros & Hultmark eq (2.15): "k = 1/β² − 1, where
k = C_D U∞²/u²", citing Taylor & Davies (1944). BUT the note mislabels the reference velocity: it says
"u = pore velocity" / "u = velocity in the pores". It is NOT the pore velocity — it is the superficial /
approach (face) velocity u = Q/A_gross. The interstitial pore velocity is u_pore = u_face/β. Energy
derivation (Taylor–Davies, as used by Steiros): Δp = ½ρu_pore² − ½ρu_face² = ½ρu_face²(1/β² − 1), so
k = 1/β² − 1 is referenced to the FACE velocity. If you instead plug the pore velocity into
Δp = ½ρ k u², you overpredict Δp by a factor 1/β² (≈11× at β = 0.3). Fix: in eqs 4 and 7, u is the
superficial (face/approach) velocity, not the pore velocity.

**5. Steiros–Hultmark porous-plate model — CONFIRMED.**
All extracted verbatim from the JFM 853 R3 PDF: C_D = (4/3)(1−u*)(2+u*)/(2−u*) [eq 2.12b];
C_D = u*²(1/β² − 1) − (4/3)(1−u*)³/(2−u*)² [eq 2.15]; E = u*/(2−u*) with U_wake = E·U∞ [eq 2.5,
"U_w = EU∞"]; validity β ≥ 0.23 (Castro 1971 shedding threshold — confirmed in text); solid-plate
C_D = 4/3 = 1.33 predicted vs 1.44 ± 0.4 measured (confirmed verbatim, §3.2).

**6. Wake velocity reduction vs porosity — CONFIRMED.**
Independently solved the coupled system (2.12b)=(2.15) by bisection for E = u*/(2−u*):
β = 0.20 → 0.141, 0.23 → 0.165, 0.30 → 0.224, 0.40 → 0.316, 0.50 → 0.420, 0.60 → 0.533. Matches the
note's 0.14 / 0.17 / 0.22 / 0.32 / 0.42 / 0.53 to rounding. The associated C_D column (1.29, 1.27, 1.22,
1.13, 1.00, 0.84) also reproduced. Bed-shear inference "~2–6 % of ambient" = (0.15–0.25)² is consistent.

**7. Porous-wall momentum sink (implementation) — CONFIRMED.**
Darcy–Forchheimer f = −(μ/K)u − ρ·C_F·|u|·u is the standard form. Thin-screen jump Δp = ½ρ k u² with
k = 1/β² − 1 is correct SUBJECT TO the face-velocity caveat in item 4 (u = superficial velocity). The
critical Reynolds number Re_c ≈ 1000 for cylindrical-element screens/fabrics above which friction losses
are negligible is confirmed (Steiros & Hultmark citing Hoerner 1952, textile Res. J. 22(4)). C_D ≈ 1.0
per rod element at Re_D > 1000 is a reasonable engineering value (isolated-cylinder C_D ≈ 1.0–1.2 over
Re ≈ 10³–2×10⁵).

**8. Emergent-array interior velocity balance — CORRECTED (density range only; formula CONFIRMED).**
The balance g·S = ½·C_D·a·U_in²/(1−φ) with a = n·D and φ = n·(π/4)·D² is confirmed. Verified against an
independent restatement of the Nepf drag law (Buono/Katul, arXiv:2403.12232, eq 7):
S_f = C_d·mD/(1−φ_v)·U²/(2g) with φ_v = mπD²/4 and frontal density = mD — algebraically identical
(multiply by g, a = mD). CORRECTION: "Nepf range φ = 0.01–0.35" is overstated. Nepf's own
characterization of natural aquatic vegetation gives solid volume fraction φ ≈ 0.001–0.01 for marsh
grass (frontal density a ≈ 0.01–0.07 cm⁻¹), and the dimensionless density she spans is ad ≈ 0.001–0.1,
i.e. φ = (π/4)·ad ≈ 0.001–0.08; pneumatophore field/lab canopies sit at φ ≈ 0.005–0.04. So the realistic
Nepf/vegetation range is φ ≈ 0.001–0.05 (up to ~0.1 for dense lab arrays), NOT up to 0.35. Note this is
below the printed-shell porosities of interest here anyway (β = 0.2–0.4 open ⇒ φ_solid = 0.6–0.8), so
the emergent-array formula is being applied well outside its calibrated density range — use it only as a
sparse-rib bulk-drag estimate, and prefer the Steiros–Hultmark screen model (§3) for the solid shell.

Sources used for verification (opened independently):
- Steiros & Hultmark (2018), *Drag on flat plates of arbitrary porosity*, J. Fluid Mech. 853 R3 — full-text PDF, eqs 2.5, 2.12b, 2.15, §3.2. https://fluids.princeton.edu/pubs/Steiros_2018.pdf
- HEC-RAS 2D Sediment docs — Soulsby (1997) settling-velocity formula & d* definition. https://www.hec.usace.army.mil/confluence/rasdocs/d2sd/ras2dsedtr/6.5/model-description/fall-velocity-and-settling/particle-settling-velocity
- Nazaroff & Alvarez-Cohen settling-tank theory (Dartmouth ENGS 37) — overflow rate, ideal & 1−exp(−w/v_c) efficiencies. https://cushman.host.dartmouth.edu/courses/engs37/Settling.pdf
- Buono, Katul et al. (2024), arXiv:2403.12232 — independent restatement of the Nepf emergent-array drag balance (eq 7, φ_v = mπD²/4). https://arxiv.org/pdf/2403.12232
- Sharqawy/ITTC seawater property correlations — kinematic viscosity vs temperature (used to check the ν label). https://ittc.info/media/2017/75-02-01-03.pdf
