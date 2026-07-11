# Deposition/Trapping Validation Benchmarks — Closing the Gap in the Validation Ladder

## Overview

The existing validation ladder (lid-driven cavity, channel log-law, cylinder Cd/St, backward-facing
step, Roulund pile scour) exercises hydrodynamics and **erosion** only. The product metric — how much
sand settles inside/around a printed shape and how fast — runs through the **deposition pathway**:
settling flux, the `D = w_s * c_b` bed boundary condition, pickup–deposition asymmetry, and
settling/erosion **lag** (a concentration profile needs time `~h/w_s` to adapt; models without lag
mis-predict trench migration by metres, see Benchmark A). This note assembles three quantitative
deposition-dominated benchmarks with full SI parameter tables and pass tolerances, plus cheap
closed-form pre-gates. **These must sit between the Roulund scour gate and any shape-ranking
milestone** — a solver can pass every erosion gate and still be off by 2–5x on trapped volume.

Key literature evidence for realistic tolerances: SISYPHE-2D reproduced trench infill to <10% of
trench depth only after adding a settling-lag term (without it, migration error was 4.5–5 m over 10 h
= ~50x the infill error); van Rijn's own SEDPIT engineering model, calibrated (factor 0.5 on upstream
transport), still over-predicted a field trench infill by ~42% (49.6 vs 35 m3/m). So "deposited volume
within ±30%, timescale within factor 2" is a demanding-but-achievable gate; ±10% is state of the art
for tuned 2DV models.

## Key equations & results

All plain-text math. rho_s = 2650 kg/m3, s = rho_s/rho_w, nu ~ 1e-6 m2/s. u*,1 = bed shear velocity
inside the deepened area, h1 = flow depth inside it, h0/v0 upstream depth/velocity, ws = settling velocity.

### 1. Trapping efficiency of a deepened area, cross flow (van Rijn 1987, from 300 SUTRENCH runs)
    e_s = 1 - exp( -A_vr * L * d / h1^2 )
    A_vr = 0.25 * (ws/u*,1) * (1 + 2*ws/u*,1)
e_s = (b0*qs,0 - b1*qs,1,min)/(b0*qs,0); L = effective settling length = B/sin(alpha_1) (B = width
between mid-slope points), d = h1 - h0 = deepening. Derived for v0 = 1 m/s, h0 = 5 m, alpha_0 = 15–90 deg,
d = 2–10 m, B = 50–500 m, ws = 0.0021–0.036 m/s, ks = 0.2 m. Stated accuracy: e_s error ~25% for a 20%
error in approach velocity (Hs/h in 0–0.3, ks/h in 0.02–0.06). e_s -> 0 as d -> 0.
[van Rijn, leovanrijn-sediment.com "Basics of channel deposition/siltation", Sec. 4.2]

### 2. Suspended-load relaxation over a deepened area (Eysink & Vermaas 1983)
    qs,x = (b0/b1)*qs,0 - [ (b0/b1)*qs,0 - qs,1 ] * [ 1 - exp( -A_ev * x / h1 ) ]
    A_ev = 0.015 * (2*ws/u*,1) * [1 + (2*ws/u*,1)] * [1 + 4.1*(ks/h1)^0.25 ]
qs,1 = equilibrium (saturated) transport in the channel; x = distance into the channel. This is the
closed-form answer to "how fast does suspended load relax toward equilibrium when the flow decelerates"
— exactly the physics our advection-diffusion + settling solver must reproduce. [same source]

### 3. Deposition rate per unit channel length immediately after dredging (van Rijn 1987)
    S = ( e_s * qs,0 + e_b * qb,0 ) * sin(alpha_0),   e_b ~ 1.0 (bed load is always trapped)
qs,0, qb,0 in kg/(s·m) per unit channel length; alpha_0 = angle between flow and channel axis. [same source]

### 4. Volume-of-cut rule (Trawle & Herbich 1980; van Rijn 2006/2012) — order-of-magnitude sanity check
    V_deposit_per_year = alpha * V_cut,   alpha = 0.12 +/- 0.06 (US offshore sandy channels, d50 0.1–0.3 mm)
alpha ~ 0.04 for Euro-Maas channel (tidal 0.5–0.7 m/s). A printed shape's annual infill should land in
the same alpha decade when its interior volume is treated as "volume of cut". [same source]

### 5. Morphodynamic skill score — the pass/fail metric for bed profiles
    BSS = 1 - mean( (zb_model - zb_measured)^2 ) / mean( (zb_initial - zb_measured)^2 )
Classification (van Rijn et al. 2003 / Sutherland et al. 2004, Coastal Eng. 51): BSS > 0.8 excellent,
0.6–0.8 good, 0.3–0.6 reasonable/fair, 0–0.3 poor, < 0 bad (worse than "no change" prediction).

### 6. Deposition-pathway unit test (analytic, quiescent/well-mixed)
    dc/dt = -ws*c/h  =>  c(t) = c0 * exp( -ws * t / h )
Well-mixed suspension, depth h, bed BC deposition flux D = ws*c_b with pickup E = 0 (tau_b < tau_cr).
Pass: L2 error < 2% over 3 e-folding times, global sand mass (suspended + bed) conserved to < 0.1%.
Companion: Camp–Hazen mixed-basin trap efficiency eta = 1 - exp(-ws*L/(h*U)) (see research/07).

### 7. Backfilled-equilibrium target for waves (Sumer et al. 2013, JWPCOE 139(1):9-23)
    S_eq/D = 1.3 * [ 1 - exp( -0.03 * (KC - 6) ) ],   KC = Um * Tw / D
Backfilling ends at the SAME equilibrium depth as ordinary scour at the final KC (their central
experimental finding). Backfill timescale >> scour timescale for KC_f < O(10).

### 8. Normalized morphological time for scour/backfill (Baykal et al., ICCE 2014, Eq. 30)
    t* = t * sqrt( g * (s-1) * d50^3 ) / D^2
Used to compare backfill time series across scales.

## The benchmarks

### Benchmark A (primary gate): van Rijn (1986) trench siltation flume, Delft Hydraulics
Widely used to validate Delft3D (Lesser et al. 2004), SISYPHE (Knaapen & Kelly 2011), ROMS/COAWST.
Parameter table (SI), from Knaapen & Kelly's replication of van Rijn (1986), JWPCOE 112(5):

| Quantity | Value |
|---|---|
| Flume | 17 m long x 0.3 m wide x 0.5 m deep |
| Water depth h | 0.255 m |
| Depth-avg current U | 0.18 m/s (following the waves) |
| Regular waves | H = 0.08 m, T = 1.5 s |
| Sand | d50 = 0.10 mm, d90 = 0.13 mm, rho_s = 2650 kg/m3 |
| ws (suspended) | 0.007 m/s (Knaapen prints "0.07 m/s" — typo; Soulsby (1997) formula gives 0.008 m/s for 0.1 mm) |
| Upstream sand feed | 0.0167 kg/(s·m) (maintains equilibrium bed upstream) |
| Trench | depth 0.125 m, side slopes 1:12 (some replications quote 1:10) |
| Measured | bed profiles at t = 0 and t = 10 h; velocity (ADV) + concentration (siphon) profiles at 5 stations at t ~ 0 |

Reported model skill: SISYPHE with settling-lag matched trench position and both slopes; centre infill
under-predicted 6.5 mm (< 10% of the 0.125 m trench depth). Without lag: centre migrated 4.5–5 m too
far and infill badly over/under-predicted — the deposition pathway (not hydrodynamics) is what this
case discriminates.

**Pass criteria (gate):** (i) deposited volume in trench at t = 10 h within ±30% of measured;
(ii) centreline infill depth within ±20% of measured; (iii) trench-centre migration distance within
±25%; (iv) BSS >= 0.3 to pass, >= 0.6 target. Also compare the 5 initial concentration profiles:
depth-integrated suspended transport within factor 1.5 at every station.

**Grid/runtime feasibility:** needs dz ~ 10–12.5 mm (>= 20 cells over 0.255 m) — finer than the 2.5 cm
production voxel. Run 10 x 0.3 x 0.35 m subdomain at 10 mm voxels = ~1.1 M cells (trivial on the 4090);
10 h physical time with morphological acceleration factor 5–10. Because the case includes small waves,
schedule the full case with the oscillatory-forcing milestone; before that, run the current-only
variant of the same geometry and gate on Eq. 1/2 (trapping efficiency within ±30% of e_s formula).

### Benchmark B (product-analog gate): Du et al. (2025), perforated artificial-reef units around a monopile
This is the "perforated hollow unit" study (arXiv:2503.13860) — 3x3 arrays of perforated cubes (C-AR,
surface porosity 19.63%, 3x3 holes per face) and hemispheres (H-AR, porosity 31.25%, 10 holes) around
a pile; directly analogous to our printed shapes. Live-bed test table (SI):

| Quantity | Value |
|---|---|
| Sediment pit | 2.4 m long x 0.5 m wide, 0.20 m deep |
| Water depth | 0.10 m |
| Velocities | U = 0.20 / 0.25 / 0.30 m/s -> theta = 0.027 / 0.043 / 0.062 (theta_cr = 0.043): clear-water, threshold, live-bed |
| Sand | d50 = 0.235 mm, sigma_g = 1.6, s = 2.65 |
| Pile | D = 0.05 m; units 0.05 m cubes/hemispheres, hole dia. 8.33 mm (cube = D/6) and 12.5 mm (hemisphere = D/4) |
| Layout | 3x3 array around pile, 0.01 m gaps (variants: 0.025 m gaps; +2 upstream rows) |
| Duration | 1 h to equilibrium (2 h runs showed no further change) |

Measured equilibrium scour depths (upstream / downstream / lateral, cm):
control 3.0/1.5/3.0 (0.20 m/s), 5.1/3.2/4.6 (0.25), 6.1/4.3/5.7 (0.30);
C-AR 0/0/1.8, 3.2/0/3.2, 3.2/1.9/4.1; H-AR 1.9/0/2.3, 2.9/0/2.9, 3.9/1.4/2.7.
Scour reduction 28–100% (the "37–100%" band cited in research/06); wake velocity reduced 50–80%; gap
jets 1.5–1.6 U0; hole jets ~1.3 U0; C-AR units displaced up to 2.5 cm at 0.30 m/s; sand deposition
mound behind arrays (H-ARs partially self-buried into their own deposit — exactly our self-ballasting
mechanism).

**Pass criteria (gate):** (i) control equilibrium scour depths at all 3 velocities within ±30%;
(ii) protected-case upstream scour within ±30% or ±0.7 cm (measurement floor), and the measured
C-AR-vs-H-AR ranking preserved at each velocity; (iii) net deposition (positive bed change) predicted
downstream of the array with peak deposit height within factor 2; (iv) time to 90% of equilibrium
scour within factor 2 of measured (~1 h). **Grid:** 4–5 mm voxels needed to resolve 8.3 mm holes
(2.4 x 0.5 x 0.25 m at 4 mm = ~4.7 M cells — inside budget); do NOT attempt at 2.5 cm production
resolution.

### Benchmark C (field-scale desk check): Scheveningen trial pipeline trench, North Sea, 1964
From van Rijn's channel-siltation note (SEDPIT example): trench 700 m long, bottom width ~10 m, side
slopes ~1:7, depth ~2 m below seabed, local depth 8 m below MSL; sand d50 = 0.2 mm (d90 = 0.3 mm);
tidal currents perpendicular to trench, peak flood 0.6 m/s / ebb 0.5 m/s; representative Hs = 1 m,
Tp = 7 s. **Measured infill: 12 m3/m after 48 days; 35 m3/m after 173 days** (overall trapping
efficiency ~11%). SEDPIT with calibration factor 0.5 computed 49.6 m3/m at 173 d (+42%). Too large for
a 3D run (domain 10 m); use it to validate our implementation of Eqs. 1–3 as an engineering
cross-check layer (predict 35 m3/m within factor 2) — the same layer we will use to sanity-check
full-scale extrapolations of printed-shape infill.

### Benchmark D (wave milestone): backfilling of a current-generated scour hole (Sumer et al. 2013; Baykal et al. 2014)
D = 0.04 m pile, d50 = 0.17 mm, h/D = 2. Initial condition: steady-current scour, V = 0.41 m/s,
u* = 0.019 m/s, theta = 0.13 (theta_cr = 0.05, live-bed), equilibrium S/D = 0.91. Then waves only:
KC = 10 (Um = 0.225 m/s, Tw = 1.79 s, theta = 0.15) -> S/D decays 0.91 -> 0.50 (50 wave periods) ->
0.35 (100 T) -> 0.25 (190 T) -> equilibrium ~ Eq. 7 value over ~800 T (Sumer Test 32). KC = 20
(Um = 0.20 m/s, Tw = 4 s, theta = 0.09) -> equilibrium S/D ~ 0.45 in ~300 T (Test 27). **Pass:**
equilibrium backfilled S/D within ±25% of Eq. 7; time to 50% backfill within factor 2. This is the
only benchmark that tests deposition INTO an existing scour hole — the self-healing claim.

## Practical guidance for our simulator

1. **Milestone placement.** Insert a new gate "M-DEP" immediately after the Roulund pile-scour gate
   and strictly before any shape-ranking milestone: M-DEP.1 = analytic unit tests (Eq. 6 settling
   column, Camp–Hazen basin, Rouse-profile deposition/pickup balance; mass conservation < 0.1%);
   M-DEP.2 = current-only trench + Eq. 1/2 trapping-efficiency check (±30%); M-DEP.3 = Du et al.
   live-bed array case (the product analog). Benchmark A (full, with waves) and Benchmark D slot into
   the oscillatory-forcing milestone. Shape ranking is unlocked only when M-DEP.2 and M-DEP.3 pass.
2. **Settling/erosion lag is mandatory physics.** Our 3D advection–diffusion + settling solver has lag
   built in *if* deposition is computed as ws*c_b from the local resolved near-bed concentration and
   pickup from local tau_b — never from an equilibrium-transport shortcut. Benchmark A quantifies the
   penalty for getting this wrong (metres of spurious migration).
3. **Tolerances**: deposited/trapped volume ±30%; local deposit/infill depth ±20–30%; timescales
   factor 2; BSS >= 0.3 (pass) / 0.6 (target). Tighter than ±20% on volume is not supported by the
   literature (even calibrated reference models sit at ~±40% in the field).
4. **Concrete numbers for our range** (Soulsby ws at 10 C, nu = 1.35e-6 m2/s is ~10% lower than the
   20 C values here): d50 = 0.1 mm -> ws ~ 0.007–0.008 m/s; 0.235 mm -> ~0.026 m/s; with U = 0.5–2.5 m/s
   and u* ~ 0.02–0.10 m/s, ws/u* spans 0.07–1.3 — i.e. our sands straddle the suspension threshold, so
   both Eq. 1 (suspension trapping) and e_b = 1 bedload trapping matter.
5. **Benchmark grids are special cases**: 10 mm voxels (Benchmark A), 4–5 mm (Benchmark B). Keep voxel
   size a config parameter; both fit in < 5 M cells.
6. **Report trapped volume like the trench literature**: V_trap(t) in m3 per unit width (or per shape),
   plus e_s, so simulator output is directly comparable to Eqs. 1–4 and to future flume tests of
   printed shapes.

## Sources

- van Rijn, L.C. (1986). "Sedimentation of dredged channels by currents and waves." J. Waterway, Port, Coastal and Ocean Eng. 112(5):541-559. (primary flume experiment)
- Knaapen, M.A.F. & Kelly, D.M. (2011). "Modelling sediment transport with hysteresis effects." XVIII Telemac & Mascaret User Club. https://eprints.hrwallingford.com/812/1/HRPP549.pdf (full Benchmark-A parameters + skill numbers)
- van Rijn, L.C. (2025). "Basics of channel deposition/siltation." https://www.leovanrijn-sediment.com/papers/Channelsedimentation2013.pdf (Eqs. 1-4, Scheveningen field case, SEDPIT)
- Du, S. et al. (2025). "An experimental study of using artificial reefs as scour protection around an offshore wind monopile." arXiv:2503.13860. https://arxiv.org/html/2503.13860v1 (Benchmark B)
- Baykal, C., Sumer, B.M., Fuhrman, D.R., Jacobsen, N.G., Fredsoe, J. (2014). "Numerical modeling of backfilling process around monopiles." Proc. ICCE 34. https://icce-ojs-tamu.tdl.org/icce/article/view/7204 (Benchmark D numbers)
- Sumer, B.M., Petersen, T.U., Locatelli, L., Fredsoe, J., Musumeci, R.E., Foti, E. (2013). "Backfilling of a scour hole around a pile in waves and current." JWPCOE 139(1):9-23. https://ascelibrary.org/doi/10.1061/(ASCE)WW.1943-5460.0000161
- Lesser, G.R., Roelvink, J.A., van Kester, J.A.T.M., Stelling, G.S. (2004). "Development and validation of a three-dimensional morphological model." Coastal Eng. 51:883-915. https://www.sciencedirect.com/science/article/abs/pii/S0378383904000870 (Delft3D validation on the same trench cases)
- Sutherland, J., Peet, A.H., Soulsby, R.L. (2004). "Evaluating the performance of morphological models." Coastal Eng. 51:917-939. (BSS classification)
- van Rijn, L.C. (1986). "Mathematical modeling of suspended sediment in nonuniform flows." J. Hydraulic Eng. 112(6):433-455. https://www.leovanrijn-sediment.com/papers/P2-1986b.pdf (scanned; concentration-profile verification data)
