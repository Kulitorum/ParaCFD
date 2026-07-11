# RESEARCH.md — Physics Specification for the Scour-Protection Simulator

This is the single authoritative physics spec for the simulator. Every equation below was
extracted from the literature by a research pass and then **independently re-verified against
primary sources** (constants, exponents, validity ranges). Detailed derivations, tables, and
citations live in `research/01-…16-*.md`; each section links to its note. When implementing,
use the constants exactly as written here — do not "improve" them from memory.

**Product context.** COBOD is developing 3D-printed concrete shapes placed on the seabed
around offshore wind foundations. The shapes must decelerate sediment-laden tidal flow so
sand settles **inside/around** them (self-ballasting, self-healing scour protection). The
simulator ranks candidate shapes by how much sand they trap and keep under increasing flow.

---

## 1. Physical setting & canonical parameters

| Quantity | Value | Notes |
|---|---|---|
| Water density ρ | **1027 kg/m³** | Seawater S=35, 10 °C (ITTC). Config parameter. |
| Kinematic viscosity ν | **1.36e-6 m²/s** | 10 °C seawater. ⚠ 1.05e-6 (used in some note tables) is ~20 °C water; ν must be a config parameter. Range 5–15 °C: 1.56–1.19e-6. |
| Sediment density ρs | 2650 kg/m³ | Quartz. s = ρs/ρ = **2.58**; (s−1)g = 15.50 m/s² |
| Grain sizes d50 | **0.2 / 0.35 / 1.0 / 5.0 mm** | Canonical design matrix (North Sea sands 0.2–0.6 mm dominate). Full supported range 0.1–10 mm. |
| Bed porosity p | **0.36** | Critical packing fraction c_pack = 0.64 (FLOW-3D; range 0.55–0.70) |
| Angle of repose φ | **32°** (sand), 35–40° (gravel) | tan 32° = 0.6249 |
| Depth-avg current U | 0.5–1.8 m/s ambient; **2.0–2.5 m/s = pile-amplified/extreme only** | 99th-pct at 9 UK farms: 0.54–1.77 m/s (research/12) |
| Domain | 10 × 10 × 5 m, uniform voxels h = 5 cm (ranking) / 2.5 cm (hero) | See §9: the domain is a *sub-problem*, not a full monopile |

**Scale reality check (research/12):** modern monopiles are 8–11.5 m diameter with 40–66 m
protection footprints — a 10×10×5 m domain cannot contain one. Valid sub-problems:
(a) a single printed unit or small cluster (footprint 1–3 m, frontal area ≤ 3 m² for ≤ 6% blockage);
(b) a spanwise-periodic strip of units;
(c) a protection-annulus sector near the pile wall, driven with **1.5–2.0× amplified inflow
and 10–20% turbulence intensity** (potential-flow amplification u/U∞ = 1 + R²/r²; bed-shear
amplification O(4), locally up to O(10) at 45° in steady current).

---

## 2. The per-timestep loop (operator-split, as in all production scour codes)

```
fluid step(s)  →  bed shear τ_b (log-law wall function, EMA-averaged)
              →  slope-corrected critical Shields θ_cr
              →  bedload q_b + pickup E + deposition D
              →  bed update (packed-fraction field, Exner-equivalent, × MORFAC)
              →  avalanche sweep (32°/30° hysteresis)
              →  refresh solid mask / interface normals
```
Morphology is updated every N = 5–100 flow steps using time-averaged shear (§7). Reference:
research/10 (Delft3D, sedExnerFoam, REEF3D, FLOW-3D all use this loop).

---

## 3. Fluid solver (research/08, /16)

Stam-style operator splitting, upgraded to engineering grade. All four upgrades are
**mandatory** — plain 1999 stable-fluids has effective numerical viscosity ~0.5·|u|·h ≈
0.03–0.06 m²/s, four orders above seawater and 10–100× any eddy viscosity.

1. **MAC staggered grid** (p at centers, u/v/w on faces). Eliminates checkerboard pressure.
2. **MacCormack semi-Lagrangian advection** (Selle et al. 2008): two SL sweeps,
   `φ_new = φ̂_{n+1} + ½(φ_n − φ̂_n)`; clamp to the min/max of the 8 trilinear corners of the
   first sweep (or revert on overshoot); **always revert to 1st order within CFL·h of solids**.
   2nd-order verified at CFL 0.75 and 1.75. Backtrace with RK2; clip rays at solid faces.
   CFL_adv ≤ **2.0 ranking / 1.0 hero**, dt = CFL·h/max|u| with max|u| ≈ 1.7·U∞
   (→ dt ≈ 18–88 ms at h = 5 cm).
3. **Pressure Poisson via MGPCG** (McAdams, Sifakis & Teran 2010): CG preconditioned with one
   geometric V-cycle (factor-2 coarsening, damped Jacobi ω = 2/3 interior smoother, 2 extra
   Gauss–Seidel boundary-band sweeps at the finest level). Converges one order of magnitude
   per ~2 iterations, grid-size independent; budget 8–14 iterations to ‖r‖/‖b‖ ≤ 1e-4.
   Jacobi (Harris GPU Gems, 40–80 sweeps) is a *visual* criterion — bring-up/debug only.
   For validation-gate runs tighten to ‖r‖/‖b‖ ≤ 1e-6 (iteration count unconstrained);
   post-projection max|div u| scales with the solve tolerance (~tol·U/h), so gates must pair
   the divergence bound with the tolerance used (research/08).
   7-point stencil; solid neighbors drop out (∂p/∂n = 0); pin p = 0 at the outlet and rescale
   outlet flux to match inlet flux before every solve (all-Neumann compatibility).
4. **Smagorinsky LES**: ν_t = (Cs·Δ)²·|S|, |S| = √(2 S_ij S_ij), Δ = h, **Cs = 0.10–0.12**.
   Treat diffusion explicitly (dt limit h²/6ν_t ≈ 0.23 s ≫ advective dt); drop Stam's
   implicit diffusion solve entirely.

**Obstacles:** free-slip on voxelized solids (constrain u·n = 0 only) — the ~1 mm viscous
sublayer is unresolvable at h = 2.5–5 cm, so voxel no-slip adds fake drag. Bed friction enters
only through the wall model (§4). Optional later upgrade: variational/cut-cell projection
(Batty 2007) for sub-voxel STL fidelity.

**Vorticity confinement: OFF** (ε = 0). It is a graphics energy injector and corrupts
quantitative transport.

**Free surface:** rigid free-slip lid (w = 0, ∂u/∂z = ∂v/∂z = 0, Neumann p). Valid: Fr ≤ 0.36
and neglected drawdown ≤ 6.4% of depth at U ≤ 2.5 m/s; flag U > 2 m/s runs as ranking-only.

**Outlet:** Orlanski-*type* convective BC `∂u/∂t + U_c ∂u/∂n = 0` with constant bulk U_c
(explicit upwind, requires U_c·dt/h ≤ 1), then global flux rescale; p = 0 pin in projection.

---

## 4. Bed shear stress from the velocity field (research/04)

**Never finite-difference du/dz across the first cell** (y+ ~ 300–2000 here → severe
underestimate). Use the log-law wall function:

```
u* = κ·U_p / ln(z_p/z0),   κ = 0.40,   τ_b = ρ·u*²,  direction = local tangential unit vector
```
- Probe the **tangential** velocity U_p by trilinear interpolation at fixed normal distance
  `z_p = max(1.5·h, 2·ks)` above the reconstructed bed surface (never the raw stair-step).
- Roughness: `ks = 2.5·d50`, `z0 = ks/30 = d50/12` (flat sand). Optional ripple form-roughness
  kf = a·η²/λ, a = 8–27.7 — if added, feed **grain-only** (skin-friction) stress θ′ computed
  with z0_grain = d50/12 to all transport formulas, not the total-roughness stress.
- Regimes by ks+ = u*·ks/ν: rough > 70 (explicit formula OK — true for d50 ≥ 0.5 mm at ≥ 1 m/s);
  transitional 5–70 (d50 = 0.1–0.2 mm!): use Christoffersen–Jonsson
  `z0 = (ks/30)(1 − exp(−u*·ks/27ν)) + ν/(9u*)`, fixed-point iterate 2–5× from the rough guess.
- Gravel trap: d50 = 10 mm → ks = 25 mm > first-cell height; probe at z_p ≥ 2ks ≈ 5 cm.
- Stair-step mitigation: interpolate u* (not u), smooth τ_b tangentially (3×3) before Exner,
  rate-limit wall updates to [0.5×, 2×] per step (OpenFOAM practice).
- **Time-filter τ_b with an EMA over 1–3 s (≈ N·dt)** before erosion/transport: the Shields
  threshold is calibrated for mean stress, not instantaneous LES fluctuations (research/01).
- Inside cavities (clearance H_c by raycast along the bed normal) — four-branch rule
  (research/13): H_c ≥ 4·z_p0 → standard probe at z_p0 = max(1.5h, 2ks); 6h ≤ H_c < 4·z_p0 →
  probe at H_c/4; 3–6h → mid-gap log-law if Re_gap > 2800 else laminar τ = 6μU/H_c; < 3h →
  freeze erosion and use the sub-grid screen model (§8). Resolved perforations need ≥ 8–10
  cells across.

Calibration anchor (self-consistent with z0 = d50/12, h = 5 m — research/14 table, /16 E6):
flat periodic channel, depth-averaged Cd = [κ/(ln(h/z0) − 1)]²; d50 = 0.2 mm → Cd ≈ 0.0012,
u* ≈ 0.0345 m/s, τ0 ≈ 1.22 Pa at U = 1 m/s; d50 = 1.0 mm → u* ≈ 0.040 m/s. Across
d50 = 0.1–10 mm: u*/U = 0.033–0.052, τ0 = 1.1–2.8 Pa. ⚠ The often-quoted shelf-sea default
Cd ≈ 0.0025 (τ0 ≈ 2.6 Pa, u* ≈ 0.050 m/s) assumes a **rippled** bed (z0 ≈ 6 mm) — do NOT use
it as the flat-sand target or you build a ~2× bed-shear bias into every transport rate.

---

## 5. Threshold of motion (research/01)

```
θ  = τ_b / ((ρs − ρ)·g·d)                       Shields parameter
D* = d·(g(s−1)/ν²)^(1/3)                        (= 20 313·d[m] at ρ = 1027, ν = 1.36e-6)
θ_cr = 0.30/(1 + 1.2·D*) + 0.055·(1 − exp(−0.020·D*))     Soulsby–Whitehouse (1997)
τ_cr = θ_cr·(ρs − ρ)·g·d                        ((ρs−ρ)·g = 15 922 N/m³ at ρ = 1027)
```
⚠ Derived constants (D* factor, (ρs−ρ)g, w_s tables, θ_cr tables) must be **computed from the
config (ρ, ν) at runtime**, never hard-coded. Reference values in research notes 01–03/16 were
tabulated at ρ = 1025, ν = 1.05e-6 (20 °C): there D* = 24 162·d, (ρs−ρ)g = 15 941 N/m³. Unit
tests must state their (ρ, ν) tuple.
Valid for all D* ≥ 0.1 (silt→cobbles). Minimum θ_cr ≈ 0.030 at D* ≈ 17 (d ≈ 0.8 mm).
Physical uncertainty at the fine end is 10–25% between published fits.

**Slope correction (mandatory — scour holes reach the repose angle).** Soulsby (1997) eq. 80a,
general 3D form, per interface cell:

```
θ_cr(β,ψ) = θ_cr,flat · [cos ψ·sin β + √(cos²β·tan²φ − sin²ψ·sin²β)] / tan φ
```
β = bed slope angle, ψ = angle of near-bed flow to the up-slope direction, φ = 32°.
Factor → 0 as β → φ for downslope flow (clamp ≥ 0.1·θ_cr,flat for overhang normals) and
→ 2cos φ ≈ 1.70 for pure upslope flow. Sanity: at β = 20°: 1.49 (upslope) / 0.39 (downslope)
/ 0.76 (transverse).

**Asymmetry rule (the product's core mechanism):** θ_cr gates **erosion and bedload only —
deposition is never thresholded**. This reproduces the Hjulström erosion/deposition hysteresis
and is exactly why sand that settles inside a sheltered shape stays there.

---

## 6. Sediment transport

### 6.1 Settling velocity (research/02)
Soulsby (1997) — single formula for the whole range, matches Ferguson–Church within 10%:
```
w_s = (ν/d)·[√(10.36² + 1.049·D*³) − 10.36]        (10.36² = 107.33)
```
At ν = 1.05e-6 (20 °C tables): w_s = 7.3 / 24.5 / 117 / 283 mm/s for d = 0.1 / 0.2 / 1 / 5 mm.
At 10 °C fine-sand w_s drops ~10–20% — compute from config ν, don't hard-code the table.
Hindered settling: multiply by (1 − c)^4.7 above c ≈ 0.001; clamp c ≤ 0.35.

### 6.2 Suspended load (research/02)
One volumetric concentration field c per active grain size:
```
∂c/∂t + ∇·(u c) − ∂(w_s c)/∂z = ∇·((ν_t/σ_s)∇c)         σ_s = 0.7 (uncertainty knob 0.5–1.5)
```
Settling = extra downward advection (fold into the MacCormack advection velocity for c).
Suspension criterion (van Rijn 1984): u*_cr,susp/w_s = 4/D* for 1 < D* ≤ 10; = 0.4 for D* > 10.
Consequence: 0.1–0.2 mm sand is suspended through most of the tidal range; **d ≥ 2 mm gravel
never suspends → bedload-only path (skip its c field, saves GPU work)**.

**Bed exchange — flux BC in the bed-adjacent cell, never a Dirichlet c_a:**
```
Pickup   E = 0.00033·ρs·[(s−1)g·d50]^0.5·D*^0.3·f_D·T^1.5   [kg/m²/s]
         T = (τ′_b − τ_cr)/τ_cr;   f_D = 1 for θ′ ≤ 1, f_D = 1/θ′ for θ′ > 1  (van Rijn 2019
         damping — REQUIRED at our 1.5–2.5 m/s cases or erosion is overpredicted)
Deposit  D = w_s·c_b   [m/s, volumetric — c is volumetric; ×ρs only if you need kg/m²/s]
         (c_b = concentration in first cell above bed; hindered w_s if c_b > 0.001)
```
Accuracy of E is factor ~2 (α = 0.00033 ± 30%). Reference-concentration cap c_a ≤ 0.05 vol.
At h = 2.5–5 cm the first cell center ≈ van Rijn reference height a = 0.01·h_water — consistent.

### 6.3 Bedload (research/03)
Dimensionless transport: Φ = q_b/√((s−1)g·d50³).
- **Default: Wong & Parker (2006) corrected MPM** — `Φ = 3.97·(θ′ − 0.0495)^1.5`
  (alt. best fit 4.93·(θ′ − 0.047)^1.60). Valid 0.4–29 mm plane bed. ⚠ The classic MPM
  coefficient 8 double-counts form drag — over-predicts ~2×.
- **Toggle: Engelund–Fredsøe (1976)** `Φ = 5p·(√θ − 0.7·√θ_cr)`,
  p = [1 + ((π/6)·0.51/(θ − θ_cr))⁴]^(−1/4) — this is what the Roulund et al. (2005) benchmark
  used; enable it for the M7 validation gate. (Do not use the 18.74 closed form except near
  threshold — it overshoots ~3.6× at θ = 1.)
- van Rijn 1984 (0.2–2 mm): Φ = 0.053·T^2.1/D*^0.3 for T < 3, **Φ = 0.100·T^1.5/D*^0.3 for
  T ≥ 3** (both branches required — our flows exceed T = 3 over fine sand).

Direction: local tangential flow, deflected by transverse slope (Talmon et al. 1995):
`tan α = tan α_τ − (1/(0.85·√θ))·∂z_b/∂n` (β₂ = 0.85, tunable to 1.6).
Bedload formulas are calibrated to θ ≈ 0.25–0.5; at θ > 1 suspended load must carry the
transport (it does, if §6.2 is implemented).

---

## 7. Bed representation & morphodynamics (research/13, /10, /03, /16)

**Decision: a single 3D packed-grain-fraction field `f_pack` on the voxel grid** — not a
heightfield. Published precedent: FLOW-3D's FAVOR sediment model (Flow Science Report 03-14),
which auto-creates packed bed wherever suspension settles, including inside structures —
exactly our self-ballasting mechanism. A column sum of this model reduces *exactly* to the
Exner equation, so all heightfield literature remains applicable.

- Normalized fill `F = f_pack/(0.64·(1 − φ_printed))`. Cell states: PRINTED / PACKED (F ≥ 1,
  hydraulically solid) / INTERFACE (0 < F < 1) / FLUID.
- Interface normal & area from the box-filtered fill: `n_b = −∇F̃/|∇F̃|`, `A_b = |∇F̃|·V_cell`.
  All bed formulas evaluate on INTERFACE cells.
- Deposition: `df_pack/dt = +(A_b/V)·w_s,hindered·c`. Erosion (Winterwerp 1992 / FLOW-3D):
  `u_lift = 0.018·D*^0.3·(θ − θ′_cr)^1.5·√((s−1)g·d)` along n_b.
- Bedload runs as upwind donor transfers of f_pack between INTERFACE cells along the flow
  tangent (MPM flux, layer thickness h_b = 0.3·d·D*^0.7·T^0.5).
- Packed cells are solid: Ergun permeability 3e-11–7.6e-8 m² makes seepage negligible.
  Partially packed cells get an implicit Ergun/Darcy–Forchheimer momentum sink.
- **Avalanching (mandatory for stability):** trigger where slope > 32°, relax to 30°
  (hysteresis, Roulund 2005). Generalize to overhangs via layered columns (multi-layer
  heightfield per (x,y)): pairwise mass-conserving transfers
  dV = 0.5·k·A·(Δz − L·tan 30°), k = 0.5, Gauss–Seidel sweeps to convergence, blocked by
  PRINTED cells.

**Exner form (for tests & reasoning):** `(1 − p)·∂z_b/∂t = −∇·q_b − E/ρs + D`, p = 0.36.
⚠ **Unit convention (dimensional trap):** E is in kg/m²/s (van Rijn pickup, §6.2) and takes
the /ρs; D = w_s·c_b is **already a volumetric flux [m/s]** (c is volumetric) and must NOT be
divided by ρs — doing so undercounts deposition by 2650× and silently kills the product's
core mechanism. q_b is solids volume flux [m²/s]. (research/02 §11, /16 E4.)
Discretize ∇·q_b **upwind along transport direction** (central is unstable — sedExnerFoam).

**Stability limits:** per bed update |Δz_b| ≤ 0.05–0.1·h and ≤ 5% of local water depth
(Delft3D default DzMax = 0.05); morphological CFL ≤ 0.9 on bed celerity c_bed ≈ 3q_b/((1−p)·h_dep).

**Time acceleration (research/16):** bed-update cadence N (flow steps per morphology step,
5–100, EMA shear window 1–3 s) is *de-noising only*; **MORFAC M is the only time accelerator**:
`(1−p)·∂z_b/∂t = M·[…]`, M ≤ 10 steady current, M ≤ 5 tidal/reversing (critical MORFAC is >10×
lower in reversing flow). Enforce M·N·dt·max|dz_b rate| ≤ 0.05·h adaptively. Ramp M from 1 over
1–2 flow-throughs after a 3–5 flow-through frozen-bed spin-up. Verify one M=1 vs M=target run
per scenario (Brier Skill Score ≥ 0.95).

**Run-length protocol (research/16):** the scour/fill timescale
`T = T*·D²/√(g(s−1)d50³)`, T* = (1/2000)(δ/D)·θ^(−2.2) (δ = 5 m domain height) spans 66 s
(U = 2.5, d = 0.2 mm) to ~50 h (U = 0.5 clear-water). The θ exponent is uncertain between
−1.3 (Whitehouse) and −2.2; Larsen & Fuhrman (2023) support θ^(−3/2) — treat clear-water T
values as **upper bounds** for run-length/MORFAC planning and never calibrate the model to
the −2.2 curve alone (research/16 E5, /05 E6). **Fixed-duration runs produce ranking
artifacts.** Run t ≥ 0.5·T of morphological time, then fit both S(t) = S_eq(1 − e^(−t/T)) and
the Welzel hyperbolic a·(1 − 1/(1+bt)); expect 10–20% S_eq error at 0.5T; never fit < 0.3T;
extend if the two fits disagree > 15%. Rank shapes on fitted equilibria (S_eq, V_trapped,eq)
with CIs, at matched t/T within a scenario.

---

## 8. Boundary conditions & forcing (research/14)

- **Inlet:** rough log-law `u(z) = (u*/κ)·ln(z/z0)` with `u* = κ·U_d/(ln(h_dom/z0) − 1)`
  (u*/U_d = 0.033–0.052 over our d50 range; the whole 5 m slab is log-layer since real tidal
  boundary layers are 20–45 m thick). Take δ = 5 m in the T* formula.
- **Inlet turbulence (ranking runs):** Jarrin (2006) Synthetic Eddy Method — σ = 1.0 m (= h/5),
  N = 150 eddies, Reynolds stresses from Nezu–Nakagawa profiles (u′_rms/u* = 2.30·e^(−z/h),
  v′: 1.63, w′: 1.27, −u′w′ = u*²(1 − z/h)), amplitude scaling γ ≤ 2.5 to reach 10–20% TI for
  near-pile scenarios. SEM divergence is cleaned by the existing projection. SEM carries ~15%
  wall-shear bias for ~6δ of fetch — acceptable for ranking, **not for validation gates**.
- **Validation/hero runs:** precursor periodic channel with the same solver (periodic x/y,
  body force g_x = u*²/h, PI mass-flux controller, ~1100 s spin-up), record inlet planes.
- **Top:** rigid free-slip lid. **Sides:** free-slip; blockage ≤ 6% or apply U_eff = U_d/(1−BR).
- **Tidal reversal:** swap inlet/outlet faces through a 60 s cosine ramp at slack (≥ 3–6
  flow-throughs); scale u*(t) with |U_d(t)|; never body-force an obstacle domain.
- **Sub-grid porosity (walls thinner than ~3 cells or perforations < 8 cells across):**
  thin-screen pressure jump `Δp = ½·ρ·k·u_face²` with **k = 1/β² − 1, u = face (superficial)
  velocity — NOT pore velocity** (using pore velocity overpredicts Δp by 1/β², ~11× at β = 0.3),
  or equivalent Darcy–Forchheimer sink in porous voxels.

---

## 9. What makes a shape trap sand (design physics, research/07, /06, /05)

- **Trapping criterion:** a grain deposits when residence time exceeds settling time:
  `u_in ≤ w_s·(L/H)` (L = interior path length, H = fall height). Camp–Hazen efficiency
  `η = 1 − exp(−w_s·A/Q)`: η > 0.8 needs w_s·A/Q > 1.6 — the interior must present a large
  plan footprint relative to throughflow.
- **Porosity optimum β ≈ 0.25–0.35 (open-area ratio):** three independent literatures converge
  here. Interior velocity drops to ~15–25% of ambient (bed shear to 2–6%, since τ ∝ u²) while
  staying above the β ≈ 0.23 threshold below which a von-Kármán vortex street forms and
  re-suspends the catch. Solid walls reflect flow and scour their own toe; permeable ones are
  low-reflection (K_r ≈ 0.2–0.4). Steiros–Hultmark (2018) gives the closed-form
  porosity→wake-velocity map for design pre-screening.
- **Openings matter:** field sediment-trapping structures (brushwood dams, permeable groynes)
  always include openings so sediment-laden flow can *enter*; sedimentation 0.15–0.5 m/season
  observed (Demak).
- **Empirical support:** perforated hollow units around a model pile (19.6–31.3% porosity,
  open bottoms) cut scour 28–100% and deposited sand while self-burying (Du et al. 2025;
  measured table in research/15 Benchmark B — the 37–100% band in research/06 is superseded);
  waves deposit sand *inside* porous rock berms (Petersen 2015). Failure modes to measure:
  edge scour (wake vortex pair, governing in current), winnowing, sinking (Horns Rev 1:
  0.5–1.5 m), displacement of tall unanchored units.
- Deposition KPI sanity number: c = 0.002 vol of 0.2 mm sand fills a still cavity at
  ~0.27 m/h of bed rise.
- Targeting: bias designs to capture d50 ≥ 0.25 mm (t_s ≈ 15 s over 0.5 m); 0.1 mm is hard
  (t_s ≈ 70 s), silt is a stretch goal.

---

## 10. Validation ladder & acceptance gates (research/11, /15)

| # | Test | Reference | Gate |
|---|---|---|---|
| V1 | Lid-driven cavity | Ghia et al. 1982 tables (in research/11) | Re=100: RMS centerline u,v error < 5%. Re=1000 is a **diagnostic** for the advection scheme (log it, don't hard-fail — semi-Lagrangian is marginal there; research/11) |
| V2 | Poiseuille + open-channel log-law | analytic (grain roughness z0 = d50/12 — see §4 anchor) | u* within 10%; recovered κ within 5% |
| V3 | Cylinder Cd/St | Cd 1.33–1.40 & St 0.164–0.168 @ Re100; Cd 1.1–1.2, St ≈ 0.19–0.20 @ Re1e4 | Re=100: Cd within 5–10%, St within 3–5%; Re=1e4 (LES): Cd within 15%, St within 10% |
| V4 | Backward-facing step | x1/S = 5.0 @ Re200 (Armaly) | reattachment within ~20% |
| V5 | Settling column (still water) | c(t) = c0·e^(−w_s·t/h) | L2 < 2%, mass error < 0.1% |
| V6 | Equilibrium Rouse profile | R = w_s/(β·κ·u*) analytic | profile shape within tolerance band |
| V7 | **Pile scour (Roulund et al. 2005)**: D = 0.1 m, h = 0.4 m, V = 0.46 m/s, d50 = 0.26 mm, live-bed V/Vcr = 1.25 | measured S/D = 1.25, ~2 h to equilibrium | S/D upstream ±15%, downstream ±30%, timescale factor 2. Grid: ≥ D/10 cells (REEF3D: D/10 → <4% err; D/5 unstable). Use Engelund–Fredsøe bedload + suspended load (disabling suspended halves scour). |
| V8 | **Deposition gate "M-DEP"** (current-only, pre-waves): Du et al. 2025 perforated units (d50 = 0.235 mm, U = 0.2–0.3 m/s, measured scour-reduction tables in research/15) + current-only trench variant checked against van Rijn's trapping efficiency e_s = 1 − exp(−A_vr·L·d/h1²) | trapped volume ±30%, depth ±20–30%, BSS ≥ 0.3 (target 0.6); e_s within ±30% | **must pass before any shape ranking**. ⚠ The full van Rijn 1986 trench experiment includes waves (H = 0.08 m, T = 1.5 s) + sand feed 0.0167 kg/s/m — it moves to the wave milestone with V9 (0.18 m/s current alone is below threshold for 0.1 mm sand) |
| V9 | Wave backfill (later, with waves) | Sumer 2013: S/D 0.91 → 0.25 in 190 periods at KC = 10 | qualitative backfill + equilibrium at final-KC value |

Benchmarks V7/V8 need special fine grids (4–10 mm voxels — lab scale); production runs use
2.5–5 cm. Anchor empirical expectations: live-bed cylinder scour S/D = 1.3 (σ = 0.7);
clear-water ramp S/D = C·tanh(h/D)·f(U/U_cr) with onset at U/U_cr = 0.5 — C = 2.0 is the
conservative *design envelope*; expect simulated equilibria near C = 1.3–1.5 (research/05
verification); wave scour S/D = 1.3·(1 − e^(−0.03(KC−6))).

---

## 11. Known limitations (be honest in every report)

1. Semi-Lagrangian solvers under-resolve the horseshoe vortex; published models underpredict
   HSV bed shear by ~30% even with body-fitted RANS. Expect to **calibrate** (via pickup α,
   Cs, or an HSV shear multiplier) against V7, then freeze constants.
2. Realistic absolute accuracy after calibration: ±15–30% on depths/volumes, factor 2 on
   timescales. The tool's authority is **comparative ranking** of shapes under identical
   conditions — treat absolute numbers as indicative until flume-validated.
3. Pickup function: factor-2 accuracy; θ_cr fine-sand uncertainty 10–25%; σ_s scatter 0.5–1.5.
4. The 10 m domain cannot host pile-wake shedding (period ~50 s, wavelength ≫ domain);
   near-pile scenarios use amplified mean inflow + elevated TI instead (§1, §8).
5. Voxel stair-stepping biases local τ_b; mitigations in §4 reduce but don't erase it.
6. Plan-B: if calibrated shear maps cannot reproduce the required amplification pattern, the
   fluid core swaps to LBM (D3Q19/27, TRT Λ = 3/16 or cumulant + Smagorinsky; ~15–30 min per
   300 s at h = 2.5 cm on the RTX 4090; wall shear is *local & 2nd-order* in LBM — research/09).
   Keep the fluid behind a `FluidCore` interface from day one.
   ⚠ FluidX3D is license-blocked for commercial use — reimplement techniques, never copy code.

## 12. Note index

| research/ | Topic |
|---|---|
| 01 | Incipient motion, Shields, slope corrections |
| 02 | Settling velocity, suspended load, pickup/deposition |
| 03 | Bedload, Exner, avalanching, MORFAC |
| 04 | Bed shear from CFD (wall functions, roughness) |
| 05 | Scour physics around piles (HSV, S/D, timescales, KC) |
| 06 | Scour protection SOTA & failure modes |
| 07 | Permeable structures & sediment trapping |
| 08 | Stable fluids 3D, MacCormack, MGPCG, LES |
| 09 | LBM alternative, GPU perf, licenses |
| 10 | Morphodynamic coupling in production codes |
| 11 | Validation benchmarks (hard numbers) |
| 12 | Site conditions (metocean, seabed, monopile scale) |
| 13 | 3D packed-fraction bed model (gap-fill) |
| 14 | Inflow/SEM/BCs (gap-fill) |
| 15 | Deposition validation benchmarks (gap-fill) |
| 16 | Timestep/MORFAC/run-length protocol (gap-fill) |
