# RESEARCH.md — Physics & Numerics Specification for WindCFD

This is the single authoritative physics/numerics spec for **WindCFD**, a GPU (CUDA) 3D
incompressible-flow LES solver (`namespace windcfd::core`). Every model below is grounded in the
solver that actually exists in `src/core/fluid/` and `src/core/geometry/`; forward-looking
quantities the solver must still be extended to compute are called out in a clearly-labelled
section (§7). Where a value is quoted, use it exactly as written — do not "improve" constants from
memory. Fluid-relevant derivations live in the `research/` notes (index in §10).

**Product context.** COBOD 3D-prints concrete buildings. Because a printer lays down curved and
straight walls at the same cost, a building's corners can be **sharp or rounded** at will. WindCFD
quantifies how **wind** flows around such buildings in an **atmospheric boundary layer (ABL)** and
what **wind loads** it imposes on the façade and roof, so that **rounded vs sharp corner** designs
can be compared on drag, roof uplift, and peak cladding suctions. This is bluff-body building
aerodynamics in air — not water, not sediment.

> **⚠ Known configuration gap (fix before any wind run).** The code's `Config` defaults still carry
> the **inherited seawater** fluid properties from the fork's origin: `config.h`,
> `stam_fluid_core.h`, `channel_periodic.h`, and `bedshear.h` all default `rho = 1027 kg/m³` and
> `nu = 1.36e-6 m²/s`. For wind these **must** be set to **air** (§1). This is a defaults bug to fix,
> not a description of the intended physics. Likewise the wall model still derives its roughness
> length from a `d50` "grain size" (`z0 = d50/12`); for wind, `z0` must be supplied directly as an
> **aerodynamic terrain roughness length** (§5.2).

---

## 1. Problem & objective

**Objective.** For a 3D-printed building immersed in a turbulent ABL, compute the surface pressure
field and integrate it into the engineering wind loads (§7), then rank **sharp-cornered vs
rounded-cornered** variants of the same building under identical inflow. Rounding corners is expected
to delay/soften separation, narrow the wake, and cut both mean drag and the peak roof-edge suctions
that govern cladding design — WindCFD is the tool that measures that difference.

**Fluid = air.** All physics is for air at roughly sea-level, ~15 °C:

| Quantity | Value | Notes |
|---|---|---|
| Air density ρ | **1.225 kg/m³** | 15 °C, 101.3 kPa. **Config parameter** (currently mis-defaulted to seawater — see callout). |
| Kinematic viscosity ν | **1.5e-5 m²/s** | μ ≈ 1.8e-5 Pa·s; ν = μ/ρ. Weak T-dependence; keep it a config parameter. |
| Reference wind speed U_ref | **10–40 m/s** | Mean/gust design wind at reference height (building/roof height, or 10 m). |
| Building scale B, H | **B ≈ 5–15 m, H ≈ 3–12 m** | Low- to mid-rise printed structures; B = characteristic width. |
| Reynolds number Re = U_ref·B/ν | **~10⁶–10⁷** | Fully turbulent bluff-body regime; loads are largely Re-independent for **sharp** corners but Re-sensitive for **rounded** ones (curved-surface separation). |
| Dynamic pressure q = ½ρU_ref² | **≈ 245 Pa @ 20 m/s**, ≈ 980 Pa @ 40 m/s | Normalises all pressure/force coefficients (§7). |

**Length/velocity scales.** The relevant scales are the building width/height (B, H), the reference
speed U_ref, the ground friction velocity u\* and aerodynamic roughness z0 of the approach terrain
(§5), and the ABL turbulence intensity I_u = σ_u/U (≈ 10–20 % near the ground, decaying with height).

**Computational domain (bluff-body best practice, COST 732 / AIJ).** Size the domain around the
building, not to a fixed box: ≥ **5H** clear upstream of the windward face, ≥ **10–15H** downstream to
let the wake develop, ≥ **5H** lateral and top clearance, and total **blockage < 3 %** (else correct
U_ref by the blockage ratio). Uniform voxels edge `h`; resolve the building with enough cells across
B and H (≥ 10–20) and, for **rounded** corners, fine enough that the voxel staircase does not
masquerade as a sharp corner (§6, §9). The grid, domain, and voxel size are all `Config` fields
(`config.h`).

---

## 2. Governing equations

Incompressible Navier–Stokes for a Newtonian fluid (air, constant ρ, low-Mach ⇒ incompressible):

```
∇·u = 0
∂u/∂t + (u·∇)u = −(1/ρ)∇p + ν ∇²u + f
```

`u = (u,v,w)` [m/s], `p` static pressure [Pa], `f` body force [m/s²] (buoyancy is neglected — neutral
ABL). At Re ~ 10⁶–10⁷ the flow is fully turbulent and the equations are solved in **Large-Eddy
Simulation (LES)** form: apply a spatial box filter of width Δ = h (the grid). The filtered momentum
equation is

```
∂ū/∂t + (ū·∇)ū = −(1/ρ)∇p̄ + ∇·[(ν + ν_t)(∇ū + ∇ūᵀ)] + f
∇·ū = 0
```

The unresolved sub-grid stress is modelled by an eddy viscosity ν_t (§3), folded into the diffusion
term. Overbars are dropped hereafter; all fields are the resolved (filtered) fields.

---

## 3. Turbulence closure (LES)

**Sub-grid model — Smagorinsky** (`turbulence.cu`, `smagorinsky_nut_gpu`):

```
ν_t = (C_s·Δ)² · |S|,   |S| = √(2 S_ij S_ij),   S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i),   Δ = h
```

Strain rates are formed at cell centres from the staggered faces. Default **C_s = 0.10–0.12** for
wall-bounded bluff-body flow; laminar validation gates (cavity, low-Re cylinder) run **C_s = 0**.
ν_t is added to ν in the explicit diffusion of the MAC velocity (`mac_diffuse_gpu`). Diffusion is
treated **explicitly** — the stability limit dt ≤ h²/(6(ν+ν_t)) is far looser than the advective
limit at these scales, so Stam's implicit diffusion solve is dropped entirely.

**Wall-adjacent treatment.** At Re ~ 10⁶–10⁷ with h ~ 0.1–0.5 m the first off-surface cell sits at
y⁺ of order 10³–10⁴ — the viscous sublayer is unresolvable. Therefore:
- **Ground:** an equilibrium **log-law wall function** supplies the surface shear (§5.2), never a
  finite-differenced ∂u/∂z across the first cell.
- **Building surfaces (production):** free-slip on the voxelised solid (constrain only u·n = 0). At
  voxel resolution a no-slip surface would add spurious numerical drag from an unresolved boundary
  layer; for bluff bodies with fixed (corner) separation the loads are pressure-dominated and only
  weakly sensitive to the surface-friction detail. A resolved **no-slip** surface option exists
  (`SolidBC = SOLID_NOSLIP`) and is used only where the boundary layer is resolved (e.g. the laminar
  cylinder validation, which needs it to shed).
- A future wall-function-on-the-building upgrade (and/or cut-cell surfaces, §9) is the path to
  faithful skin friction and rounded-corner separation.

---

## 4. Numerical method (Stam stable-fluids, engineering grade)

Operator-split incompressible solve on a staggered grid. Plain 1999 stable-fluids is unusable here
(its implicit-advection numerical viscosity dwarfs both ν and ν_t); the four upgrades below are
**mandatory**. Reference: `research/08`.

**4.1 Staggered MAC grid** (`mac_grid.h`). Pressure at cell centres; u, v, w on the x-, y-, z-faces.
Eliminates checkerboard pressure decoupling. All indexing, BC-aware ghost fetch, and BC-aware
trilinear sampling are pure `__host__ __device__` helpers so the CUDA kernels and their CPU
reference twins share identical arithmetic (GPU-vs-CPU parity is a hard test rule).

**4.2 MacCormack semi-Lagrangian advection** (`advect.cu`, Selle et al. 2008). Two semi-Lagrangian
sweeps combined as `φ_new = φ̂_{n+1} + ½(φ_n − φ̂_n)`, then **clamped** to the min/max of the eight
trilinear corners of the forward interpolation (overshoot control). Reverts to first-order
semi-Lagrangian within a `band` of cells of any wall/solid (`near_wall`). RK2 backtrace; rays are
clipped to the domain (`clamp_to_domain`). Second-order accurate away from boundaries, unconditionally
stable in the advection step.

**4.3 Pressure projection via MGPCG** (`mgpcg.h`, `project.cu`; McAdams, Sifakis & Teran 2010).
Conjugate gradient preconditioned by **one geometric multigrid V-cycle**: factor-2 coarsening,
damped-Jacobi (ω = 2/3) interior smoother with `pre_post_jacobi = 2` sweeps, extra boundary-band
Gauss–Seidel sweeps (`gs_band = 2`), and a coarsest-level Jacobi solve. Grid-size-independent
convergence; production tolerance ‖r‖/‖b‖ ≤ **1e-4**, tightened to **1e-6** for validation gates.
The 7-point Laplacian drops solid/box neighbours (∂p/∂n = 0). The closed-box case (cavity) is
all-Neumann/singular — handled by a zero-mean RHS and pressure-gauge removal; the open channel pins
`p = 0` at the outlet and rescales outlet flux to inlet flux before each solve (compatibility). A
plain damped-Jacobi path is kept for bring-up/debug only.

**4.4 Time stepping / CFL.** Adaptive `dt = CFL·h / max|u|` with max|u| the peak resolved speed
(around a bluff body, local corner acceleration gives max|u| ≈ 1.5–2·U_ref). Advective CFL ≤ **1**
for LES accuracy (the scheme is stable to ~2). Explicit-diffusion stability dt ≤ h²/(6(ν+ν_t)) is
checked but rarely binding. A fixed dt mode exists for uniform-rate sampling (Strouhal FFT).

**Per-step loop:** advect (MacCormack) → add body force / diffuse with ν+ν_t (explicit) → apply
boundary conditions & any porous sink → **project** (MGPCG) to enforce ∇·u = 0 → (optional) sample
surface pressure/shear for load statistics (§7). Vorticity confinement is **off** (a graphics energy
injector that corrupts quantitative transport).

---

## 5. Boundary conditions & forcing

The closed-box core (`stam_fluid_core.*`, cavity) uses six-face wall BCs (no-slip / free-slip /
moving-lid). The building-aerodynamics core is the **open channel** (`channel_core.*`,
`channel_bc.h`): inlet on xmin, outlet on xmax, free-slip lateral and top, ground wall model on the
floor, and a voxel solid mask for the building.

**5.1 Inflow — ABL mean + synthetic turbulence.**
- **Mean profile** (`channel_inlet_u`, `INLET_LOGLAW`): a rough-wall log law
  `u(z) = (u\*/κ)·ln(z/z0)`, κ = 0.40, referenced to the ground datum so u → 0 at the surface. u\* is
  chosen to flux-match a target reference speed over the domain height (`loglaw_ustar_for_U`). A
  uniform top-hat inlet (`INLET_UNIFORM`) is available for canonical benchmarks (cylinder). *Target
  (§7): add a selectable **power-law** mean profile `u/u_ref = (z/z_ref)^α` (α ≈ 0.11–0.33 by terrain
  category) and terrain-category presets, so the inflow matches a design-code ABL.*
- **Synthetic turbulence — SEM** (`sem_inlet.*`, Jarrin 2006). N compact tent eddies convect through
  a box straddling the inlet; their superposition is coloured by the Cholesky factor of a
  Reynolds-stress tensor to produce divergence-carrying inflow fluctuations that the MGPCG projection
  cleans each step. Amplitude scaling γ (≤ 2.5) drives the inlet turbulence intensity up to the
  10–20 % ABL target. *The stresses currently follow **Nezu–Nakagawa open-channel** profiles
  (`nezu_*`); for wind these should be retuned to an **ABL** turbulence-intensity/length-scale profile
  (Eurocode/ESDU or measured spectra) — see §9.*
- **Precursor inflow** (`precursor.*`): record inlet planes from a precursor periodic channel and
  replay them (interpolated, looped) at the building-domain inlet, removing the SEM adjustment-length
  bias for high-fidelity runs.
- **Periodic driver** (`channel_periodic.*`): a flat-ground channel periodic in x/y, driven by a body
  force g_x = u\*²/h with a **PI mass-flux controller** locking the bulk speed; interior closure is a
  Prandtl mixing length l = κz√(1−z/h) giving a depth-filling neutral-ABL log law. This is both the
  precursor generator and the §8 log-law validation case.

**5.2 Ground wall function — the ABL surface model** (`bedshear.h`, `research/04`). The ground/terrain
enters only through an equilibrium log-law wall model, not through a resolved boundary layer:

```
u\* = κ·U_p / ln(z_p/z0),   κ = 0.40,   τ_w = ρ·u\*²   (tangential to the surface)
```

- Probe the tangential speed U_p by interpolation at a **fixed normal distance** z_p = max(1.5h, 2k_s)
  above the ground (never the raw stair-stepped first cell).
- **Roughness** is the terrain aerodynamic roughness: z0 (open sea ~2e-4 m, open terrain ~0.03 m,
  suburban ~0.1–0.3 m, urban ~0.3–1 m). *The code currently derives z0 = d50/12 from a legacy
  "grain size"; for wind, z0 must be supplied directly (§9).* Roughness Reynolds regimes are handled
  by k_s⁺ = u\*·k_s/ν: fully-rough explicit branch above 70, transitional (5–70) by a
  Christoffersen–Jonsson z0 fixed-point iteration.
- Anti-stair-step mitigations: interpolate u\* (not u), tangential smoothing of τ_w, per-step
  update rate-limiting, and an EMA time-filter of τ_w — so the wall stress reflects the mean, not
  instantaneous LES fluctuations. A near-surface **gap/clearance** classifier (log-law vs laminar
  Couette for thin gaps under overhangs/eaves) is present in the wall model for flow squeezed beneath
  building features.

**5.3 Outflow** (xmax). Orlanski-type convective condition `∂u/∂t + U_c·∂u/∂n = 0` (constant bulk
U_c, explicit upwind, requires U_c·dt/h ≤ 1) followed by a **global flux rescale** to match the inlet
flux; pressure Dirichlet p = 0 pin in the projection. Place ≥ 10–15H downstream so the wake is not
clipped.

**5.4 Top & lateral faces.** Rigid **free-slip** lid and free-slip side walls (zero normal velocity,
zero tangential shear). Keep them ≥ 5H from the building and blockage < 3 %, or correct U_ref by the
blockage ratio.

**5.5 Solid obstacle (the building).** A voxelised solid-cell mask (§6). On solid faces the normal
velocity is zeroed (u·n = 0); the tangential condition is **free-slip** in production, **no-slip**
where a resolved boundary layer is wanted. The mask can be hot-swapped without resetting the flow
(`update_solid`). Perforated/thin features below the ≥ ~8-cell resolution floor are represented by a
**sub-grid porous momentum sink** (`channel_porous.h`): a thin-screen pressure jump
Δp = ½·ρ·k·u_n², k = 1/β² − 1 (β = open-area ratio), applied as a linearized-implicit Forchheimer
drag on the face-normal velocity before projection.

---

## 6. Geometry pipeline (STEP → mesh → voxel mask)

1. **STEP import** (`step_import.*`). OpenCascade reads the CAD building (STEP, native millimetres)
   and triangulates it to a Qt-free `TriMesh` in **metres** (BRepMesh linear deflection ~0.1 mm
   default). `load_step_mesh` returns the merged display/voxelization mesh; `load_step_solids` returns
   one mesh per solid (for future per-piece handling). OCC is isolated in the `windcfd_geometry`
   static lib so the solver stays OCC-free.
2. **Placement** (`model_placement.h`). A single shared affine transform `world(v) = M·v + t`
   (rotation·scale + translation) positions the mesh; `place_model_on_bed` centres it in x/y and rests
   its lowest vertex on the ground. The viewer and the voxelizer use the **same** transform so the
   solid mask lands exactly where the building is drawn.
3. **Ray-parity voxelizer** (`voxelize.*`). Watertight inside test by an axis-aligned **+z ray**,
   even-odd (parity) count of triangle crossings, robust to edge grazing via inclusive barycentric
   tests. Each cell is probed with a K×K×K sub-point grid (default K = 3) and marked solid iff its
   inside-fraction > 0.5 (**supersampled majority fill**); the per-cell fraction is retained as a
   cut-cell volume fraction for later use. A **thin-wall safeguard** keeps sealed any wall thinner
   than half a cell (over-thicken beats leak). Output is the `ChannelBC` cell mask (1 = solid,
   0 = fluid, `g.pidx`); `voxelize_mesh_instances` ORs several placements into one rigid structure.

**Rounded-vs-sharp caveat.** The mask is a binary voxel staircase, so a rounded corner at coarse h is
indistinguishable from a chamfered/sharp one. Faithful rounded-corner separation needs fine h (many
cells across the fillet radius) and/or the retained cut-cell fractions promoted to real **cut-cell
surface physics** (§9).

---

## 7. TARGET wind-load quantities — **to be implemented**

The current solver computes the velocity/pressure fields, the SEM/precursor inflow, the ground wall
stress, and the voxel building mask. It does **not yet** integrate wind loads. This section specifies
the quantities WindCFD must be extended to compute; they are the deliverable of the project.

**7.1 Surface pressure coefficient.**
```
Cp = (p − p∞) / (½ρU_ref²)
```
Sampled on the building envelope (façade + roof) from the pressure field at fluid cells adjacent to
**exposed solid faces** — the exposed faces are already identifiable from the voxel mask. p∞ is the
reference static pressure and U_ref the reference speed (both at a stated reference height). Because
p is defined only to a gauge, use a consistent p∞ (e.g. inlet/freestream reference).

**7.2 Integrated force & moment coefficients.** Integrate surface traction (pressure dominates for
bluff bodies; add wall shear when a surface wall model exists) over the envelope:
```
F = ∮ (−p n + τ_w) dA,     C_D = F_x/(½ρU_ref²A),  C_L = F_z/(½ρU_ref²A),  C_S = F_y/(½ρU_ref²A)
M = ∮ r × (−p n) dA,       C_M = M/(½ρU_ref²A·L_ref)
```
A = reference area (frontal area for drag; plan area for roof uplift); L_ref = a reference length
(H for the base overturning moment). Report drag, lift/uplift, side force, and the base overturning
and torsional moments.

**7.3 Roof uplift.** Net upward roof load = ∮_roof (p∞ − p) dA (external suction), reported as an
uplift force and coefficient. Internal pressure (a GCpi-type assumption for sealed vs dominant-opening
enclosures) is a **parameter to add**, since WindCFD resolves the external flow only.

**7.4 Mean, fluctuating, and peak pressures.** After spin-up, accumulate per-face time statistics:
mean C̄p, RMS C′p, and **peak** (min/max over the record) Cp. Peak **suctions** (large-magnitude
negative Cp) at windward roof edges/corners and along leading façade edges govern cladding/fixing
design and are the headline output. Track true min/max (and/or a peak factor g ≈ 3–4).

**7.5 Rounded-vs-sharp comparison.** For matched inflow, compare between the sharp and rounded
variants: Cp fields, C_D, across-wind C_L RMS and its shedding Strouhal, peak roof-edge suctions,
**separation and reattachment** locations, the roof **corner/conical-vortex** structure, and the wake
width and base pressure. The expected result — rounding delays separation, narrows the wake, and cuts
both mean drag and peak suctions — is precisely what these diagnostics must quantify.

**7.6 ABL design inflow.** Add a selectable design-wind ABL inflow (log-law **and** power-law mean +
matched turbulence-intensity/length-scale profiles by terrain category, §5.1), so loads correspond to
a code-defined wind climate rather than a generic channel profile.

---

## 8. Validation benchmarks

**Current (fluid-core) gates** — these must pass; references in `research/11`.

| # | Test | Reference / target | Gate |
|---|---|---|---|
| V1 | Lid-driven cavity (`cavity.*`) | Ghia, Ghia & Shin 1982 centreline tables | Re = 100: RMS centreline u,v error < 5 %. Re = 1000 is a **diagnostic** for the advection scheme (log it; semi-Lagrangian is marginal there). |
| V2 | Open-channel / neutral-ABL log law (`channel_periodic.*`) | analytic log law (κ, z0) | recover u\* within 10 %; recover κ within 5 % from a log-law fit of the mean profile. |
| V3 | Flow past a cylinder (`cylinder.*`) — bluff-body drag & shedding | Re = 100: Cd ≈ 1.33–1.40, St ≈ 0.164–0.168. Re = 1e4 (LES): Cd ≈ 1.1–1.2, St ≈ 0.19–0.20 | Re = 100: Cd within 5–10 %, St within 3–5 %. Re = 1e4: Cd within 15 %, St within 10 %. Uses a **resolved no-slip** surface so the Kármán street forms. |

**Future wind-engineering validation targets** (add as §7 lands — canonical bluff-body-in-ABL cases
with published surface-pressure and load data):
- **Surface-mounted cube** (Silsoe 6 m cube; Castro & Robins 1977 wind-tunnel cube) — the archetypal
  sharp-cornered bluff body in an ABL, with extensive Cp, separation/reattachment, and roof-suction
  data. The primary rounded-vs-sharp-relevant benchmark.
- **CAARC standard tall building** — the standard rectangular tall-building model with a wind-tunnel
  database of mean/RMS surface pressures and along-/across-wind force and base-moment coefficients.
- **TPU aerodynamic database** — wind-tunnel Cp databases for low- and high-rise buildings, for broad
  envelope-pressure validation.
- **Corner-modification studies** (e.g. chamfered/rounded tall-building corners; Tamura & Miyagi 1999
  and related) — direct literature anchors for the expected rounded-vs-sharp load reduction.

Follow **COST 732** and the **AIJ** CFD-in-wind-engineering guidelines (Tominaga et al. 2008;
Franke et al. 2007) for domain sizing, grid convergence, inflow specification, and reporting.

---

## 9. Known gaps & limitations (state honestly in every report)

1. **Seawater config defaults** (callout at top): ρ, ν, and the roughness input are still fork
   legacies. Set ρ = 1.225, ν = 1.5e-5, and supply z0 as a terrain roughness before trusting any
   number.
2. **No load post-processing yet** (§7): Cp, force/moment coefficients, roof uplift, and
   mean/RMS/peak pressure statistics are unimplemented — the current core produces the fields these
   integrals will consume.
3. **Inlet turbulence is open-channel, not ABL**: the SEM Reynolds stresses use Nezu–Nakagawa
   profiles; retune to an ABL intensity/length-scale spectrum for design-representative inflow.
4. **Building surfaces are free-slip at production resolution**: no wall function on the building
   itself, so skin friction and rounded-corner separation are not yet faithfully captured.
5. **Voxel stair-stepping**: binary voxels cannot distinguish a rounded corner from a sharp one at
   coarse h. Resolving the rounding benefit needs fine grids or promoting the retained cut-cell
   fractions to true cut-cell surface physics (Batty et al. 2007-style variational projection).
6. **Comparative authority**: like all coarse-LES bluff-body tools, absolute Cp/loads carry
   uncertainty; the tool's strength is the **relative** ranking of sharp vs rounded (and other
   design) variants under identical inflow. Freeze constants (C_s, wall z0, inflow) across a
   comparison set.
7. **Swappable core**: the fluid lives behind the `FluidCore` interface (`fluid_core.h`); an LBM
   alternative (`research/09`) remains a drop-in option if the projection-based core proves limiting.
   ⚠ FluidX3D is license-blocked for commercial use — reimplement techniques, never copy code.

---

## 10. Note index

Fluid-relevant research notes retained from the fork:

| research/ | Topic |
|---|---|
| 04 | Wall shear from CFD (log-law wall functions, roughness regimes) — reframed as the ABL ground model |
| 08 | Stable fluids 3D: MAC grid, MacCormack advection, MGPCG projection, LES |
| 09 | Alternative solver (LBM) & GPU/licensing notes |
| 11 | Validation benchmarks (Ghia cavity, cylinder Cd/St, log-law channel — hard numbers) |
| 12 | Site/environmental conditions (approach-flow context) |
| 14 | Inflow: ABL profile, SEM turbulent inlet, precursor, boundary conditions |

**Key citations.** Stam 1999 (stable fluids); Selle et al. 2008 (MacCormack advection);
McAdams, Sifakis & Teran 2010 (MGPCG); Smagorinsky 1963 (LES); Jarrin 2006 (SEM);
Orlanski 1976 (convective outflow); Ghia, Ghia & Shin 1982 (cavity); Batty et al. 2007 (cut-cell,
future). Wind engineering: Richards & Hoxey 1993 (ABL inflow / Silsoe); Tominaga et al. 2008 (AIJ
CFD guidelines); Franke et al. 2007 (COST 732 best practice); the CAARC and TPU building databases;
Eurocode EN 1991-1-4 / ASCE 7 (terrain categories, wind loads).
