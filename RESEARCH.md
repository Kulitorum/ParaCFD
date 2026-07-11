# RESEARCH.md — Physics & Numerics Specification for WindCFD

This is the single authoritative physics/numerics spec for **WindCFD**, a GPU (CUDA) 3D
incompressible-flow LES solver (`namespace windcfd::core`). Every model below is grounded in the
solver that actually exists in `src/core/fluid/`, `src/core/geometry/` and `src/core/windloads.*`.
The wind-load extraction and the printed-building geometry pipeline are now **implemented** (§6, §7);
the remaining physics work (time-averaging the loads, an ABL design inflow) is called out in a
clearly-labelled section (§7.6, §9). Where a value is quoted, use it exactly as written — do not
"improve" constants from memory. Fluid-relevant derivations live in the `research/` notes (index in
§10).

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
> length from a `d50` "grain size" (`z0 = d50/12`); for wind, z0 must be supplied directly as an
> **aerodynamic terrain roughness length** (§5.2). (The wind-load module `windloads.*` already
> defaults to air: ρ = 1.225 kg/m³.)

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
boundary conditions & any porous sink → **project** (MGPCG) to enforce ∇·u = 0 → (optional) integrate
surface pressure into wind loads (§7). Vorticity confinement is **off** (a graphics energy injector
that corrupts quantitative transport).

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
  uniform top-hat inlet (`INLET_UNIFORM`) is available for canonical benchmarks (cylinder) and is
  what the current building configs use. *Target (§7.6): a selectable **power-law** mean profile
  `u/u_ref = (z/z_ref)^α` (α ≈ 0.11–0.33 by terrain category) and terrain-category presets, so the
  inflow matches a design-code ABL.*
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

## 6. Geometry pipeline (printed centerline → footprint → voxel mask)

WindCFD does **not** need a watertight CAD solid — the flow sees only the solid/fluid cell mask. So
for a 3D-printed building the geometry pipeline works directly from the **wall centerline** the printer
follows, thickening it into a solid mask by a distance test. This is inherently **leak-proof** (no
thin-wall gaps a ray-parity voxelizer could tunnel through) and makes the **corner style a first-class
parameter** — exactly the rounded-vs-sharp knob the project exists to study.

1. **STEP import** (`step_import.*`). OpenCascade reads the printed **centerline** STEP (vertical wall
   surface *ribbons* — the swept toolpath, not a closed solid), native millimetres, and triangulates
   it to a Qt-free `TriMesh` in **metres** (BRepMesh linear deflection ~0.1 mm default). OCC is
   isolated in the `windcfd_geometry` static lib so the solver and `windloads.*` stay OCC-free.
2. **Placement** (`model_placement.h`). A single shared affine transform positions the mesh; the
   viewer and the geometry stage use the **same** transform so the solid mask lands exactly where the
   building is drawn. `center_footprint` provides the simple "drop it in the middle of the wind
   tunnel" placement (footprint bbox centre → domain centre in x/y, base at `base_z`).
3. **Horizontal section → 2D footprint** (`building.cpp`, `mesh_horizontal_section`). Intersect the
   centerline ribbons with a horizontal plane at mid-height (a tiny z-nudge avoids exact
   vertex/seam hits) to get section segments, then greedily endpoint-chain them into closed
   **centerline loops** — a `Footprint` of 2D polylines in the xy-plane. This collapses the vertical
   ribbons to the plan-view wall path.
4. **Direct voxelization** (`voxelize_building`). For each cell centre within the wall height band
   (`base_z ≤ z ≤ base_z + wall_height`), the cell is marked **solid iff its distance to the footprint
   loops ≤ half the wall thickness** — i.e. the wall is the ± half-thickness band swept along the
   centerline. Corner treatment:
   - **Rounded corners** (`corner_radius > 0`): the centerline loops are **pre-filleted** by a
     procedural 2D fillet (`round_loop`, tangent circular arc of radius `corner_radius − half`), so the
     ± half-thickness band yields an **outer corner radius ≈ corner_radius**. Configurable radius; the
     bare distance band alone already rounds outer corners by half the wall thickness.
   - **True sharp corners** (`corner_radius ≈ 0`): start from the rounded distance band, then **square
     off each genuine convex corner** by adding an outward **miter wedge** triangle (edge-normal
     bisector, with a miter limit that bevels very acute spikes). This is local and robust on
     non-convex outlines — no global offset polygon to self-intersect — and leaves gently-curved walls
     rounded via a turn-angle threshold.
5. **Flat roof slab** (`voxelize_building`, roof band `top < z ≤ top + roof_thickness`). A **solid**
   flat roof = the **convex hull of the footprint**, dilated outward by the **overhang** (a cell is
   roof-solid if it is inside the hull, or within `half + roof_overhang` of it). Using the convex hull
   guarantees a filled, non-hollow roof regardless of how the section split the footprint into loops,
   and gives the eaves/overhang for free.

Output is the `ChannelBC` cell mask (size `g.p_count()`, indexed `g.pidx`, 1 = solid / 0 = fluid),
fed straight to the solver's obstacle (§5.5). A general **watertight ray-parity voxelizer**
(`voxelize.*`, +z ray even-odd parity, K³-supersampled majority fill, retained cut-cell fractions)
remains available for arbitrary closed STEP solids; the centerline pipeline above is the production
path for printed buildings.

**Rounded-vs-sharp caveat.** The mask is still a binary voxel staircase, so a rounded corner at coarse
h is indistinguishable from a chamfered/sharp one. Faithful rounded-corner separation needs fine h
(many cells across the fillet radius) and/or promoting retained cut-cell fractions to real **cut-cell
surface physics** (§9).

---

## 7. Wind-load extraction — **implemented** (`src/core/windloads.{h,cpp}`)

The solver integrates the engineering wind loads directly from the pressure field over the voxelized
building. The building is the solid-cell mask of §6; its surface is the set of **exposed voxel faces**
— faces between a solid cell and a fluid neighbour (the 6-neighbour test in `compute_wind_loads`).
The whole module is host-only and OCC-free (in `libwindcfd`); wind is assumed along **+x** (the inlet).

**7.1 Surface pressure (exposed-face integration).** For each exposed face the **surface pressure is
taken as the adjacent fluid cell's pressure** (cell-centred, physical Pa). The **reference pressure
`p_ref`** is the mean pressure over the fluid cells of the **upstream inlet slab** (the first ~2 cell
layers at xmin — the freestream). The gauge pressure on the face is `p_face − p_ref`.

**7.2 Pressure force & moment.** Only the **pressure** load is integrated — skin friction is small
for a bluff body with fixed corner separation and is omitted (consistent with the production free-slip
building surface, §3). Summing over exposed faces of area `h²` with outward unit normal `n̂`
(solid → fluid):

```
F = − Σ_faces (p_face − p_ref) · n̂ · h²          (pressure pushes inward)
M = Σ_faces r × dF ,   r = face-centre − building-centroid
```

**7.3 Coefficients.** With dynamic pressure `q = ½ρU_ref²` (ρ, U_ref from `WindLoadParams`, default
air ρ = 1.225):

```
Cp = (p − p_ref)/q            (per exposed face; a per-solid-cell mean Cp is also exported for colouring)
Cd = F_x /(q·A_frontal)       drag  (+x)
Cs = F_y /(q·A_frontal)       side  (+y)
Cl = F_z /(q·A_plan)          lift / roof uplift (+z)
CMx,CMy,CMz = M /(q·A_ref·L_ref)   moments about the building centroid
```

Reference areas are the **projected** solid extents: `A_frontal` = frontal (y–z) projection,
`A_plan` = plan (x–y) projection, each `= (projected cell count)·h²`; `L_ref` = building height
(solid z-extent). Drag/side normalise on the frontal area, lift/uplift on the plan area.

**7.4 Convention check (why Cp needs no extra factor).** The projection updates velocity as
`u −= (dt/ρ) ∇p` (`mac_ops.h::subtract_gradient`, §4.3), so the solved `p` is a **physical static
pressure in Pa**, not a kinematic `p/ρ`. Therefore `Cp = (p − p_ref)/q` with `q = ½ρU_ref²` is
dimensionally correct as written — no extra ρ factor. This is the single easiest place to introduce a
factor-of-ρ error, so it is called out explicitly.

**7.5 Validation of the load path.** On a roughly symmetric test "house" the integration reproduces
the physically expected signatures: **windward stagnation Cp ≈ 1.0** (Bernoulli — the front face
recovers the full dynamic pressure), **Cd ≈ 1.2** (blunt bluff-body drag), and a **side force ≈ 0**
(symmetry). These are the sanity gates on the exposed-face integrator itself.

**7.6 Not yet done — makes the current live loads provisional.**
- **Time-averaging (the key remaining physics step).** `compute_wind_loads` returns an
  **instantaneous** snapshot; the GUI currently reads it off an **unsettled** flow. Trustworthy
  rounded-vs-sharp comparison needs **mean / RMS / peak** loads and Cp accumulated over the
  statistically-steady window after spin-up (peak roof-edge/corner **suctions** — large-magnitude
  negative Cp — govern cladding and are the headline output; track true min/max and/or a peak factor
  g ≈ 3–4). Until this lands, treat all reported coefficients as provisional single-frame values.
- **ABL design inflow (§5.1).** The building configs run a **uniform** inlet; add a selectable design
  ABL profile (log-law **and** power-law mean + matched turbulence-intensity/length-scale by terrain
  category) so loads correspond to a code-defined wind climate rather than a uniform stream.
- **Internal pressure** for roof uplift (a GCpi-type sealed vs dominant-opening assumption) is a
  parameter still to add — WindCFD resolves only the external flow.

---

## 8. Validation benchmarks

**Current (fluid-core) gates** — these must pass; references in `research/11`.

| # | Test | Reference / target | Gate |
|---|---|---|---|
| V1 | Lid-driven cavity (`cavity.*`) | Ghia, Ghia & Shin 1982 centreline tables | Re = 100: RMS centreline u,v error < 5 %. Re = 1000 is a **diagnostic** for the advection scheme (log it; semi-Lagrangian is marginal there). |
| V2 | Open-channel / neutral-ABL log law (`channel_periodic.*`) | analytic log law (κ, z0) | recover u\* within 10 %; recover κ within 5 % from a log-law fit of the mean profile. |
| V3 | Flow past a cylinder (`cylinder.*`) — bluff-body drag & shedding | Re = 100: Cd ≈ 1.33–1.40, St ≈ 0.164–0.168. Re = 1e4 (LES): Cd ≈ 1.1–1.2, St ≈ 0.19–0.20 | Re = 100: Cd within 5–10 %, St within 3–5 %. Re = 1e4: Cd within 15 %, St within 10 %. Uses a **resolved no-slip** surface so the Kármán street forms. |

**Wind-load sanity gates** (the §7 integrator): windward stagnation **Cp ≈ 1.0**, blunt-body
**Cd ≈ 1.2**, and **side force ≈ 0** on a symmetric test house (§7.5).

**Wind-engineering validation targets** (fold in as the §7.6 time-averaging + ABL inflow land —
canonical bluff-body-in-ABL cases with published surface-pressure and load data):
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
   number. (The `windloads.*` module already defaults ρ to air.)
2. **Loads are instantaneous, not time-averaged** (§7.6): `compute_wind_loads` reports a single-frame
   snapshot read off an unsettled flow. Mean/RMS/peak accumulation over the statistically-steady
   window is the key remaining physics step — until it lands, the **rounded-vs-sharp comparison is
   not yet trustworthy**.
3. **Inlet is uniform / open-channel, not an ABL** (§5.1, §7.6): the building configs run a uniform
   inlet, and the SEM Reynolds stresses (when used) follow Nezu–Nakagawa open-channel profiles.
   Add a design-wind ABL mean (log-law + power-law) and retune the turbulence to an ABL
   intensity/length-scale spectrum for design-representative inflow.
4. **Building surfaces are free-slip at production resolution**: no wall function on the building
   itself, so skin friction is omitted (§7.2) and rounded-corner separation is not yet faithfully
   captured.
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
