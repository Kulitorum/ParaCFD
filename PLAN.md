# PLAN.md — Implementation Plan: WindCFD (wind loads on 3D-printed buildings)

GPU (CUDA) 3D incompressible-flow **LES** tool that quantifies how **wind** flows around
3D-printed concrete buildings and the wind **loads** on their facade and roof, comparing
**rounded** vs **sharp** building corners. The fluid is **air** in an atmospheric boundary
layer (ABL) — bluff-body aerodynamics, not a hydraulic problem.

WindCFD was forked from a sediment-transport simulator; **all** sediment / seabed / erosion
code has been removed. What remains is a validated, GPU-resident incompressible fluid core
plus a CAD-to-voxel geometry pipeline and a Qt/GL viewer — the foundation this plan builds a
wind-load workflow on top of. Physics constants come from the referenced notes / benchmark
papers, never from memory.

The plan is written to be executed **one milestone per coding session**. Each forward
milestone (W1…W8) is self-contained, states what to build and why, and ends with an
acceptance criterion that must pass before the next one starts.

---

## 1. Goal & scope

**Question:** for a 3D-printed concrete building in a design wind, how large are the wind
loads on the facade and roof, and how do **rounded** corners change those loads versus
**sharp** corners?

**Approach:** drive an ABL inflow (mean shear + resolved turbulence) past the building
(imported/generated as a voxelized solid), run wall-modelled LES, sample the surface
pressure on the envelope, integrate it into pressure coefficients and force/moment
coefficients, and compare two geometries that differ **only** in corner radius under matched
inflow and grid.

**In scope:** single isolated building (optionally a small cluster), incompressible air,
mean + fluctuating + peak surface pressures, drag/lift/side-force, roof uplift, base moments,
and the rounded-vs-sharp delta. **Out of scope (for now):** compressibility, thermal
stratification / buoyancy, rain/driven-particle effects, full urban context, structural
response (this produces the aerodynamic loads a structural check would consume).

**Fluid:** air — ρ ≈ 1.225 kg/m³, ν ≈ 1.5e-5 m²/s. ⚠ The current `Config` defaults still
carry the inherited **seawater** values (ρ = 1027, ν = 1.36e-6); switching them is the first
forward task (W1).

---

## 2. Current state — what already works

Everything below is **implemented, builds, and has a passing gate or smoke run**. The forward
work in §4 is genuinely future — do not assume any W-milestone feature exists yet.

**Fluid core (complete).**
- MAC staggered grid, MacCormack semi-Lagrangian advection (clamp + near-solid reversion),
  explicit Smagorinsky LES eddy viscosity, MGPCG pressure projection (Jacobi fallback),
  no-slip / free-slip walls, divergence check.
  *Gate that survives:* lid-driven cavity 128³ vs **Ghia (1982)** centreline profiles < 5%
  at Re=100; MGPCG convergence + `max|div u|` tolerances (`gate_M1`).
- Open-channel flow: log-law inlet, Orlanski-type convective outlet + flux balance + pressure
  pin, rigid lid, voxel **obstacle mask**, free-slip **and** no-slip solids.
  *Gate that survives:* circular-cylinder **Cd / Strouhal** at Re=100 & 200 (Kármán street,
  Cd ≈ 1.33–1.4, St ≈ 0.165) with < 0.1% global mass imbalance (`gate_M2`).
- Porous momentum sink (thin-screen Δp model), available on obstacle faces.
- Log-law **wall function** for the ground / ABL surface (rough + transitional-roughness
  branch, EMA time-filtering, τ smoothing, stair-step mitigations) — file `fluid/bedshear.*`,
  functionally the ground-shear/wall model.
- Periodic-channel driver: body force + **PI mass-flux controller** for a homogeneous
  turbulent channel; mixing-length closure for the mean profile.
- **SEM synthetic turbulent inlet** (Jarrin 2006, Nezu–Nakagawa stresses, Cholesky
  colouring) + a **precursor** inlet-plane record/replay.
  *Gate that survives:* flat-wall periodic channel recovers analytic u\*/κ within 10% and
  holds target turbulence intensity at mid-domain without divergence growth (`gate_M3`).

**Geometry pipeline (complete).**
- OpenCascade **STEP import** → triangulated `TriMesh` (mm→m scaling, outward-normal winding
  fix), isolated in the `windcfd_geometry` static lib (the only target that links OCC).
- Watertight **ray-parity voxelizer** (CAD mesh → solid cell mask), with an
  enclosed-volume-vs-mesh-volume sanity check and affine model placement.
- Procedural reference-shape masks (`geometry/shape_masks.h`): plate, ring, open box — a
  starting point for procedurally generated building blocks.

**IO & config (complete).** Hand-rolled VTK **ImageData (.vti)** writer/reader (appended
binary); **JSON** config loader (vendored nlohmann) with a strict allowed-key whitelist;
`namespace windcfd::core`, units SI.

**GUI `windcfd-gui` (complete).** Qt 6 + GL 4.3 slice viewer: threaded solver (GL strictly on
the main thread, render decoupled from stepping → ~60 fps), live field slices, CUDA-GL interop
colour kernel, animated flow **particles / tracers**, **video capture**, clip plane, **STEP
model load** + a full **placement gizmo** (move/rotate/scale then voxelize where placed),
scene save/restore, live inlet-speed / profile controls, and **flow reversal**. Doubles as the
**headless driver** via `--offscreen` + a family of scripted smoke flags (`--load-step`,
`--voxelize`, `--set-u`, `--apply-h`, `--autoclose-ms`, …).

**Build & run (verified).** Configures and builds green with **VS 2022 (MSVC 14.44) + CUDA
13.1 + Qt 6.11.1 + OpenCascade 8.0**, CUDA archs 75;86;89 (Ada = 89 dev box, RTX 4090).
Fast unit suite + gates pass; the headless offscreen smoke run passes and the exe is
double-clickable with the Qt + OCC runtime deployed next to it. There is **no separate CLI
target** yet — the GUI is the run driver.

---

## 3. Architecture & directory layout (current, accurate)

Four CMake targets; the physics core has no Qt and no OCC dependency.

```
WindCFD/
├── CMakeLists.txt            project(WindCFD LANGUAGES CXX CUDA); Qt auto-detect glob
├── PLAN.md  RESEARCH.md  CLAUDE.md  HANDOVER.md  research/
├── configs/                 v1_cavity_re100 · v2_channel_loglaw · v3_cylinder ·
│                            g1_viewer(_full/_highres) · m0_smoke   (+ wind configs = forward work)
├── src/
│   ├── 3rdparty/nlohmann/json.hpp
│   ├── core/                → libwindcfd (static, C++/CUDA, NO Qt/OCC)
│   │   ├── config.{h,cpp}         JSON run config; ⚠ defaults still seawater (W1 target)
│   │   ├── cuda_probe.{h,cu}      device query
│   │   ├── vti_writer.cpp / vti_reader.cpp / vti_synthetic.h   VTK ImageData IO
│   │   ├── fluid/                 the incompressible LES core (FluidCore interface)
│   │   │   ├── mac_grid.h, mac_ops.h, fluid_core.h
│   │   │   ├── advect.cu          MacCormack SL + clamp + solid reversion
│   │   │   ├── turbulence.cu      Smagorinsky ν_t
│   │   │   ├── project.cu, mgpcg.cu   pressure projection (MGPCG + Jacobi fallback)
│   │   │   ├── stam_fluid_core.cu, cavity.{cpp,h}   closed-box core + Ghia gate
│   │   │   ├── channel_*.{cu,cpp,h}  open channel: inlet/Orlanski/lid/mask, pressure,
│   │   │   │                         porous sink, periodic driver, SEM driver, empty-channel
│   │   │   ├── channel_bc.h        inlet profiles (uniform / log-law), BC sampling
│   │   │   ├── cylinder.cu         Cd(control-volume) + St(FFT) — the M2 gate probe
│   │   │   ├── bedshear.{cu,h}, bedshear_ops.h   log-law GROUND/ABL wall model
│   │   │   ├── sem_inlet.{cu,h}, inlet_fluct.h, plane_ops.{cu,h}, periodic_ops.h
│   │   │   └── precursor.{cpp,h}   inlet-plane record/replay
│   │   └── geometry/          → voxelizer in libwindcfd; STEP import split out (below)
│   │       ├── tri_mesh.h, model_placement.h, shape_masks.h
│   │       ├── voxelize.{cpp,h}    watertight ray-parity mesh → solid mask
│   │       └── step_import.{cpp,h} → compiled into windcfd_geometry (OCC-only .cpp)
│   └── gui/                  → windcfd-gui (Qt6+GL) + windcfd_gui_cuda (Qt-free interop lib)
│       ├── main.cpp, main_stub.cpp, main_window.{cpp,h}
│       ├── slice_viewer.{cpp,h}   QOpenGLWidget, GL 4.3, camera, gizmo, clip plane
│       ├── sim_worker.{cpp,h}     steps the FluidCore on a QThread (no GL)
│       ├── sim_setup.{cpp,h}      builds the sim from JSON (looser parser than core Config)
│       ├── scene_io.{cpp,h}       .scn save/restore
│       ├── flow_particles/flow_tracers.{cpp,h}, video_recorder.{cpp,h}, video_settings_dialog.cpp
│       ├── slice_field.{cu,h}     sampler + colourmap (+ CPU reference, parity-tested)
│       ├── slice_gl.{cu,h}        cudaGraphicsGLRegisterBuffer interop (only GL-header TU)
│       └── camera.h, colormap.h, gl_thread_check.h
└── tests/                   GoogleTest + ctest labels (unit + gate_M1/M2/M3, geometry, slice)
```

**Targets:** `libwindcfd` (core), `windcfd_geometry` (OCC STEP import, isolated),
`windcfd_gui_cuda` (Qt-free CUDA-GL interop), `windcfd-gui` (Qt6 viewer / headless driver).

**Memory / performance (RTX 4090, air).** Grid sizes are unchanged by the fluid switch:
~4M cells at h = 5 cm (~0.3 GB) is the working resolution; 32M cells at h = 2.5 cm (~2.5 GB)
is the hero resolution — both trivial on 24 GB, everything stays resident. Per-step cost is a
few tens of ms at 4M cells; long averaging windows for peak-pressure statistics are the real
wall-clock driver (see §5 risks).

**Toolchain (verified on this machine, 2026-07).** CUDA **13.1**
(`…/CUDA/v13.1`, arch 75;86;89); MSVC **14.44** (VS2022 Professional, C++20); CMake 4.1.1;
Qt **6.11.1** (`C:/Qt/6.11.1/msvc2022_64`); OpenCascade **8.0** (`C:/OpenCASCADE-8.0/build2`).
Configure/build/test:
```
cmake -B build -S . -G "Visual Studio 17 2022" -A x64        # or -G Ninja from a VS2022 dev prompt
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure        # gates: -L gate_M<k>; fast: -L unit
```
⚠ **Host-compiler trap:** this box also has VS 18 (MSVC 14.51), which CUDA 13.1 **rejects**.
The VS2022 generator picks 14.44 automatically; with Ninja, build from the **x64 Native Tools
prompt for VS 2022** (or pass `-DCMAKE_CUDA_HOST_COMPILER=…/14.44…/cl.exe`). Never
`-allow-unsupported-compiler`.

---

## 4. Forward milestones

Implement in order (W1→W8); a GUI/reporting sub-task may interleave once its prerequisite
data exists. A milestone is DONE only when its acceptance check passes as a test/smoke run and
the commit lands. Never weaken an acceptance check to pass it — if one looks wrong, flag it.
Keep the `FluidCore` interface clean: load post-processing consumes only the fields the core
already exposes (`{u, v, w, p, ν_t, τ_wall}`).

### W1 — Switch physics defaults to AIR + a wind config (S)
**What:** change `Config` defaults to air (ρ ≈ 1.225, ν ≈ 1.5e-5), add a reference wind speed
/ reference height concept, reinterpret the inherited roughness key (currently `d50`, mapped
to z₀ = d50/12) as an **ABL aerodynamic roughness length z₀** for a terrain category, and add
a first `configs/wind_*.json`. Extend the strict allowed-key whitelist in `config.cpp` (and
the looser `sim_setup` parser) for any new keys.
**Why:** the entire tool must default to air; every downstream coefficient (Cp, Cd, …)
normalises by ρ and U, so wrong fluid properties silently corrupt all loads.
**Acceptance:** `Config` round-trips the air defaults; the surviving `gate_M1/M2/M3` fluid
gates still pass with the air ν (artificial-ν validation cases are unaffected); an
empty-domain run with the wind config is steady and matches its inlet profile at mid-domain
within 5%; the headless smoke run passes on the wind config.

### W2 — ABL inflow profile (M)
**What:** assemble a proper atmospheric boundary-layer inlet from the existing pieces — a
log-law (or power-law) **mean** profile set by z₀ + reference speed/height, plus **SEM**
turbulence at a target intensity/length-scale, optionally seeded from a **precursor** for a
converged spectrum. Expose it as an `inlet_profile: "abl"` mode driven from config.
**Why:** design wind loads depend strongly on the incoming shear and turbulence intensity;
a top-hat inlet under-loads and mislocates separation. This is the reference "design wind".
**Acceptance:** on an empty fetch, the mid-domain mean profile matches the target
log-/power-law within ~5% and the turbulence intensity holds at its target (e.g. 10–20%)
across ≥ several flow-throughs without divergence growth; the profile + TI are documented in
the wind config and reproduced from a fixed seed.

### W3 — Building geometry: rounded vs sharp corner variants (M)
**What:** produce the 3D-printed building as a voxelized solid obstacle in **two variants
that differ only in corner radius** — sharp (r = 0) and rounded (r > 0). Two supported paths:
(i) import a building **STEP** per variant via the existing OCC pipeline + voxelizer; and/or
(ii) procedurally generate a prismatic building block with a `corner_radius` parameter
(extending `shape_masks.h`) so a radius sweep needs no CAD round-trip. Both variants must
share footprint, height, and blockage.
**Why:** the rounded-vs-sharp comparison is the whole point; it is only fair if the two
geometries are identical except at the corners and voxelize on the same grid.
**Acceptance:** both variants voxelize **watertight** (enclosed-volume vs mesh-volume within
tolerance), report matched footprint/height/blockage, and the LES diverts stably around each
with a plausible wake; a corner-radius change is visible in the voxel mask at the working
resolution.

### W4 — Surface pressure sampling → Cp (M)
**What:** identify the building's wetted envelope cells (facade + roof faces adjacent to
fluid), sample the pressure there each step, time-average and accumulate fluctuation
statistics, and convert to a **pressure coefficient** Cp = (p − p_ref)/(½ ρ U_ref²) using a
defined reference pressure and reference velocity. Export per-face Cp maps (VTI/CSV).
**Why:** every wind load derives from the surface pressure field; Cp is the standard,
grid-independent way to report and validate it.
**Acceptance:** on a validation bluff body (e.g. a surface-mounted cube), the mean Cp on the
windward/side/leeward/roof faces reproduces published values within tolerance (windward
stagnation Cp ≈ +0.7…+0.8, separated faces negative); mean **and** RMS Cp fields export
cleanly.

### W5 — Integrated wind loads (M)
**What:** integrate the surface pressure (plus, optionally, wall shear from the ground/wall
model) over the envelope to obtain **drag / lift / side-force coefficients**, **roof uplift**,
and **base overturning moments**; report mean, RMS/fluctuating, and **peak** values over a
statistically converged window. Export force/moment time series to CSV.
**Why:** these coefficients + roof uplift are the deliverable the building's structural design
consumes; peak (not just mean) pressures govern cladding and roof fixings.
**Acceptance:** on the W4 validation building, integrated force coefficients match the
benchmark within tolerance; roof uplift and base moments are reported with mean **and** peak;
force histories export to CSV and a re-run is bit-reproducible.

### W6 — Rounded-vs-sharp comparison study (M)
**What:** run the two W3 variants under **matched** ABL inflow, grid, and run length, then
quantify the load differences — Cp distributions on facade + roof, integrated force/uplift
coefficients, and peak pressures — as a paired delta (rounded − sharp).
**Why:** this answers the project's question and produces the headline result (does rounding
the printed corners reduce drag / roof uplift / peak cladding pressure, and by how much).
**Acceptance:** a paired run identical except for corner radius yields signed, magnitude-bearing
load-difference metrics; both runs are individually reproducible; the report states the effect
of rounding on drag and roof uplift with an uncertainty estimate (see W7/§5 resolution caveat).

### W7 — Validation against a wind-engineering benchmark (M)
**What:** reproduce a standard, well-documented bluff-body ABL benchmark for mean **and**
fluctuating Cp and integrated loads — e.g. the **CAARC** standard tall building, the **Silsoe /
surface-mounted cube**, or an **AIJ / TPU** aerodynamic-database case — matching that case's
specified inflow with the W2 ABL inlet.
**Why:** credibility of the rounded-vs-sharp numbers rests on the pipeline reproducing a case
with published experimental data; ideally validated **before or alongside** W6.
**Acceptance:** mean Cp and force coefficients fall within the benchmark's reported
experimental scatter at the specified grid, and the fluctuating/peak Cp trend is captured
qualitatively; the comparison + grid are documented in `configs/` and a report.

### W8 — Results, reporting & visualization (M)
**What:** a batch runner for the variant/inflow matrix plus an **HTML/CSV comparison report**
(Cp maps, load tables with mean + peak, rounded-vs-sharp deltas); visualization of the
pressure field + surface Cp + streamlines in the GUI, and an **offscreen** render path that
produces report figures with no interactive session.
**Why:** turns runs into a communicable deliverable and makes the comparison repeatable.
**Acceptance:** a single command reproduces the rounded-vs-sharp comparison report end-to-end
with figures; the offscreen path emits the images headlessly; re-running a config reproduces
the report numbers.

**Sequencing.** W1 → W2 → W3 unblock everything. W4 → W5 build the load pipeline on the
validation building; W7 validates it; W6 applies it to the real comparison; W8 packages it.
GUI/reporting increments (pressure-field colouring, Cp overlay, streamline export) can land
opportunistically once W4/W5 produce data.

---

## 5. Risks & open questions

1. **Rounded-corner physics is Reynolds-sensitive.** Sharp-edged bluff-body flow has *fixed*
   separation (Re-independent) — LES at model Re transfers to full scale. But a **rounded**
   corner's separation point *moves* with Re, and full-scale building Re is unreachable; LES
   at reduced Re can misplace separation on the rounded variant. Treat rounded-vs-sharp deltas
   as **model-scale** results, sanity-check against benchmark data (W7), and flag the caveat in
   any report — do not over-claim a full-scale drag reduction.
2. **Voxel resolution vs corner radius.** A rounded corner of radius r is only meaningful when
   r spans **several voxels**; at h = 5 cm a small radius is re-sharpened by the staircase, so
   the comparison collapses. Needs a resolution study and likely the hero grid (h = 2.5 cm) for
   the final W6 numbers; state r/h for every result.
3. **Staircase over-prediction of drag.** The voxel staircase over-predicts Cd (measured
   ~+14% on the cylinder gate; cut-cell/variational projection is deferred). Bias partly
   cancels in the *rounded − sharp difference* only if both run at **matched resolution** —
   enforce that, and report absolute coefficients with the staircase caveat.
4. **Headless / offscreen GL for report figures.** The current headless path uses Qt's
   `qoffscreen` platform, which has **no GL context**, so the interactive drag/render feel and
   any screenshot capture are unverified headlessly. W8's figure export may need an EGL/pbuffer
   offscreen GL context or a CPU raster fallback — confirm the capture path early rather than
   assuming a windowed screenshot works headlessly.
5. **Validation-data availability & inflow match.** A benchmark is only usable if its **inflow
   (z₀, profile, TI, length scale) is fully specified** so the W2 ABL inlet can reproduce it.
   Pick the case (CAARC / Silsoe cube / AIJ / TPU) on that basis; some datasets report loads
   but under-specify the approach flow.
6. **Peak/fluctuating pressures need long, converged windows.** Mean Cp settles quickly; peak
   cladding pressures need many flow-throughs of statistics, which drives wall-clock at hero
   resolution. Budget run length per the averaging requirement, and report the window used +
   convergence of the peak estimate.
7. **Config schema & inherited naming.** The core `Config` uses a **strict key whitelist** and
   inherited names (`d50` for roughness); adding ABL roughness-length, reference-height, and
   corner-radius keys means extending both `config.cpp` and the GUI's looser `sim_setup`
   parser. Decide the ABL roughness parameterization (terrain-category z₀ vs the legacy
   z₀ = d50/12 mapping) in W1 and keep the two parsers in sync.
8. **Wall model on a vertical facade.** The log-law wall function was built for the horizontal
   ground; applying/validating it on the building's vertical walls and roof (where separation,
   not an attached log layer, dominates) needs checking — near separated regions the free-slip
   / resolved-BL treatment may be more appropriate than a log-law probe.

## 6. How implementation sessions should work

- Read `CLAUDE.md`, the target W-milestone, and the relevant `research/` notes / benchmark
  paper first. Physics constants come from those sources verbatim; if a needed constant is
  unspecified, **stop and flag it in the commit** rather than inventing one.
- Every new kernel gets a CPU reference and a GPU-vs-CPU parity test (≈1e-5 rel. for float32
  reductions); keep the `FluidCore` interface clean (load code only *reads* `{u,v,w,p,ν_t,τ}`).
- Gates/acceptance checks are ctest targets or scripted headless smoke runs; a milestone commit
  must show its check output. Determinism: fixed seeds for SEM eddies; document any
  nondeterminism.
- Commit style: Conventional Commits (`feat:`, `fix:`, `test:`, …). Keep the fluid gates
  (`gate_M1/M2/M3`) green — they are the foundation every wind result stands on.
