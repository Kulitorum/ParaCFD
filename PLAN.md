# PLAN.md — Implementation Plan: WindCFD (wind loads on 3D-printed buildings)

GPU (CUDA) 3D incompressible-flow **LES** tool that quantifies how **wind** flows around
3D-printed concrete buildings and the wind **loads** on their facade and roof, comparing
**rounded** vs **sharp** building corners. The fluid is **air** in an atmospheric boundary
layer (ABL) — bluff-body aerodynamics, not a hydraulic problem.

WindCFD was forked from a sediment-transport simulator; **all** sediment / seabed / erosion
code has been removed. What remains is a validated, GPU-resident incompressible fluid core
plus a CAD-to-voxel geometry pipeline and a Qt/GL viewer. On that foundation the wind-load
workflow is now **substantially built**: air defaults are in, a centerline→building solid
pipeline (rounded **and** sharp corner variants) exists, surface Cp + integrated force/moment
loads are extracted, and the interactive GUI drives the whole load-a-house → place → Build →
Start-Simulation flow with live load readouts and Cp colouring. What remains is turning the
**instantaneous** loads into **converged, time-averaged** statistics, adding a real ABL inflow,
validating against a benchmark, and running the rounded-vs-sharp study. Physics constants come
from the referenced notes / benchmark papers, never from memory.

The plan is written to be executed **one milestone per coding session**. Each remaining
milestone (R1…R4) is self-contained, states what to build and why, and ends with an acceptance
criterion that must pass before the next one starts.

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

**Fluid:** air — ρ ≈ 1.225 kg/m³, ν ≈ 1.5e-5 m²/s. ✅ The core `Config` defaults are now air
(ρ = 1.225, ν = 1.5e-5); the old seawater defaults are gone. One caveat remains: the
`configs/building.json` GUI scenario still runs with a **moderate ν = 1e-3** (numerically
gentle at the working grid), not the real air value — retargeting it toward the true building
Re is folded into R2.

---

## 2. Current state — what already works

Everything below is **implemented and builds**. §4 splits it explicitly into DONE milestones
(with the code that realizes each) and the REMAINING milestones (R1…R4) — the remaining work is
genuinely future, do not assume an R-milestone feature exists yet.

**Fluid core (complete, inherited + validated).**
- MAC staggered grid, MacCormack semi-Lagrangian advection (clamp + near-solid reversion),
  explicit Smagorinsky LES eddy viscosity, MGPCG pressure projection (Jacobi fallback),
  no-slip / free-slip walls, divergence check. Validated against **Ghia (1982)** lid-driven
  cavity centreline profiles < 5% at Re=100.
- Open-channel flow: log-law inlet, Orlanski-type convective outlet + flux balance + pressure
  pin, rigid lid, voxel **obstacle mask**, free-slip **and** no-slip solids. Validated against
  circular-cylinder **Cd / Strouhal** (Kármán street) at Re=100 & 200.
- Porous momentum sink (thin-screen Δp model); log-law **wall function** for the ground/ABL
  surface (`fluid/bedshear.*`); periodic-channel driver with **PI mass-flux controller**;
  **SEM synthetic turbulent inlet** (Jarrin 2006) + a **precursor** inlet-plane record/replay.

**Geometry pipeline (complete).**
- OpenCascade **STEP import** → triangulated `TriMesh` (mm→m scaling, outward-normal winding
  fix), isolated in the `windcfd_geometry` static lib (the only target that links OCC).
- Watertight **ray-parity voxelizer** (CAD mesh → solid cell mask) with an
  enclosed-volume-vs-mesh-volume sanity check and affine model placement.
- **Building solid from a 3D-print centerline** (`geometry/building.{h,cpp}`) — see W-DONE B
  below. Procedural reference-shape masks (`geometry/shape_masks.h`) remain as building blocks.

**Wind loads (complete, instantaneous).** `core/windloads.{h,cpp}` integrates the surface
pressure over the building's exposed voxel faces into force/moment coefficients + a per-cell Cp
field — see W-DONE C below.

**IO & config (complete).** Hand-rolled VTK **ImageData (.vti)** writer/reader; **JSON** config
loader (vendored nlohmann) with a strict allowed-key whitelist; **air** defaults;
`namespace windcfd::core`, units SI.

**GUI `windcfd-gui` (complete, incl. the Build workflow).** Qt 6 + GL 4.3 slice viewer:
threaded solver (GL strictly on the main thread, render decoupled from stepping → ~60 fps), live
field slices, CUDA-GL interop colour kernel, animated flow **particles / tracers**, **video
capture**, clip plane, **STEP model load** + a full **placement gizmo** (move/rotate/scale then
voxelize where placed), scene save/restore, live inlet-speed / profile controls, and **flow
reversal**. The wind-load **Build workflow** is wired end-to-end — see W-DONE D below. Doubles as
the **headless driver** via `--offscreen` + scripted smoke flags (`--load-step`, `--set-u`,
`--apply-h`, `--autoclose-ms`, …). `tools/building_probe.cpp` is a headless dev tool that loads a
centerline STEP and voxelizes the building for geometry checks with no GL.

**Build & run (verified).** Configures and builds green with **VS 2022 (MSVC 14.44) + CUDA
13.1 + Qt 6.11.1 + OpenCascade 8.0**, CUDA archs 75;86;89 (Ada = 89 dev box, RTX 4090). The
headless offscreen smoke run passes and the exe is double-clickable with the Qt + OCC runtime
deployed next to it. There is **no separate CLI target** — the GUI (and `building_probe`) are the
run drivers.

---

## 3. Architecture & directory layout (current, accurate)

Four library/exe targets plus a dev tool; the physics core has no Qt and no OCC dependency.

```
WindCFD/
├── CMakeLists.txt            project(WindCFD LANGUAGES CXX CUDA); Qt auto-detect glob
├── PLAN.md  RESEARCH.md  CLAUDE.md  HANDOVER.md  research/
├── configs/                 v1_cavity_re100 · v2_channel_loglaw · v3_cylinder ·
│                            g1_viewer(_full/_highres) · m0_smoke · building  (wind-around-house)
├── src/
│   ├── 3rdparty/nlohmann/json.hpp
│   ├── core/                → libwindcfd (static, C++/CUDA, NO Qt/OCC)
│   │   ├── config.{h,cpp}         JSON run config; ✅ defaults now AIR (rho=1.225, nu=1.5e-5)
│   │   ├── windloads.{h,cpp}      Cp + force/moment integration over the voxelized building
│   │   ├── cuda_probe.{h,cu}      device query
│   │   ├── vti_writer.cpp / vti_reader.cpp / vti_synthetic.h   VTK ImageData IO
│   │   ├── fluid/                 the incompressible LES core (FluidCore interface)
│   │   │   ├── mac_grid.h, mac_ops.h, fluid_core.h
│   │   │   ├── advect.cu          MacCormack SL + clamp + solid reversion
│   │   │   ├── turbulence.cu      Smagorinsky ν_t
│   │   │   ├── project.cu, mgpcg.cu   pressure projection (MGPCG + Jacobi fallback)
│   │   │   ├── stam_fluid_core.cu, cavity.{cpp,h}   closed-box core + Ghia validation
│   │   │   ├── channel_*.{cu,cpp,h}  open channel: inlet/Orlanski/lid/mask, pressure,
│   │   │   │                         porous sink, periodic driver, SEM driver, empty-channel
│   │   │   ├── channel_bc.h        inlet profiles (uniform / log-law), BC sampling
│   │   │   ├── cylinder.cu         Cd(control-volume) + St(FFT) — cylinder validation probe
│   │   │   ├── bedshear.{cu,h}, bedshear_ops.h   log-law GROUND/ABL wall model
│   │   │   ├── sem_inlet.{cu,h}, inlet_fluct.h, plane_ops.{cu,h}, periodic_ops.h
│   │   │   └── precursor.{cpp,h}   inlet-plane record/replay
│   │   └── geometry/          → voxelizer + building in libwindcfd; STEP import split out
│   │       ├── tri_mesh.h, model_placement.h, shape_masks.h
│   │       ├── voxelize.{cpp,h}    watertight ray-parity mesh → solid mask
│   │       ├── building.{cpp,h}    centerline section → thickened wall + hull roof solid mask
│   │       └── step_import.{cpp,h} → compiled into windcfd_geometry (OCC-only .cpp)
│   └── gui/                  → windcfd-gui (Qt6+GL) + windcfd_gui_cuda (Qt-free interop lib)
│       ├── main.cpp, main_stub.cpp, main_window.{cpp,h}   incl. the Build workflow + gizmo undo/redo
│       ├── slice_viewer.{cpp,h}   QOpenGLWidget, GL 4.3, camera, gizmo, clip plane, Cp colouring
│       ├── sim_worker.{cpp,h}     steps the FluidCore on a QThread (no GL); publishes WindLoads + Cp
│       ├── sim_setup.{cpp,h}      builds the sim from JSON (looser parser than core Config)
│       ├── scene_io.{cpp,h}       .scn save/restore
│       ├── flow_particles/flow_tracers.{cpp,h}, video_recorder.{cpp,h}, video_settings_dialog.cpp
│       ├── slice_field.{cu,h}     sampler + colourmap (+ CPU reference, parity-tested)
│       ├── slice_gl.{cu,h}        cudaGraphicsGLRegisterBuffer interop (only GL-header TU)
│       └── camera.h, colormap.h, gl_thread_check.h
└── tools/building_probe.cpp  headless: centerline STEP → footprint → building solid mask (geometry check)
```

**Targets:** `libwindcfd` (core, incl. `windloads` + `building`), `windcfd_geometry` (OCC STEP
import, isolated), `windcfd_gui_cuda` (Qt-free CUDA-GL interop), `windcfd-gui` (Qt6 viewer /
headless driver), `building_probe` (geometry dev tool).

**Memory / performance (RTX 4090, air).** ~4M cells at h = 5 cm (~0.3 GB) is a comfortable
working resolution; the `building` config runs ~1.15M cells at h = 25 cm; 32M cells at h = 2.5 cm
(~2.5 GB) is the hero resolution — all trivial on 24 GB, everything stays resident. Per-step cost
is a few tens of ms at 4M cells; long averaging windows for converged/peak load statistics (R1)
are the real wall-clock driver (see §5 risks).

**Toolchain (verified on this machine, 2026-07).** CUDA **13.1** (arch 75;86;89); MSVC **14.44**
(VS2022 Professional, C++20); CMake 4.1.1; Qt **6.11.1** (`C:/Qt/6.11.1/msvc2022_64`);
OpenCascade **8.0** (`C:/OpenCASCADE-8.0/build2`).
```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build --target windcfd-gui
```
⚠ **Host-compiler trap:** this box also has VS 18 (MSVC 14.51), which CUDA 13.1 **rejects**.
Build from the **x64 Native Tools prompt for VS 2022** (14.44). Never `-allow-unsupported-compiler`.

---

## 4. Milestones

### 4a. COMPLETED — the wind-load foundation

These are implemented and land in the git history; each is verifiable against the code cited.

**W-DONE A — Air physics defaults + empty start.** ✅
`Config` (`src/core/config.h`) now defaults to **air**: ρ = 1.225 kg/m³, ν = 1.5e-5 m²/s, U = 1
m/s; the inherited seawater values are gone. The default solid-cylinder obstacle was deleted, so
a fresh domain starts **empty** (flow develops with no body until a building is Built). The
`WindLoadParams` default (`ρ = 1.225`) matches. *Remaining nuance:* `configs/building.json` still
uses a numerically-gentle `nu = 1e-3` rather than true air ν — retargeting that toward the real
building Re is part of R2.

**W-DONE B — Building geometry pipeline: rounded vs sharp corners.** ✅
`src/core/geometry/building.{h,cpp}` turns a **3D-print centerline** surface (loaded as a
`TriMesh` by the STEP importer) into a solid-cell obstacle **without** meshing a watertight
solid, chosen deliberately because the flow needs only the mask and a distance test is leak-proof:
- `mesh_horizontal_section()` slices the centerline at mid-height into a 2D **footprint** (closed
  loops, greedy segment chaining).
- `voxelize_building()` marks a cell solid when it lies **within ± half the wall thickness** of
  the footprint over the wall height (a distance-field thickening of the wall ribbons), and adds a
  **SOLID flat roof slab** over the **convex hull** of the footprint, dilated by the overhang
  (hull guarantees a filled roof regardless of how the section split into loops).
- **Rounded corners:** a configurable `corner_radius` pre-fillets the centerline with a
  tangent-arc 2D fillet (`round_loop`), so the outer corner radius ≈ radius.
- **True sharp corners:** at `corner_radius = 0` the rounded band is **squared off** at each
  genuine convex corner by adding an outward **miter-wedge** triangle (local + robust on
  non-convex outlines, with a miter limit to bevel spikes). So the **rounded-vs-sharp variants
  both exist** and differ only at the corners.
`tools/building_probe.cpp` exercises this headlessly and prints stats.

**W-DONE C — Surface Cp + integrated wind loads.** ✅
`src/core/windloads.{h,cpp}` (`compute_wind_loads`) walks the solid mask's **exposed faces**
(solid cell adjacent to fluid), takes the surface pressure as the neighbouring fluid cell's
pressure, and integrates F = −Σ (p − p_ref) n̂ h² to get **Fx/Fy/Fz**, moments about the
centroid, and the coefficients **Cd / Cs / Cl** (+ **CMx/CMy/CMz**) against q = ½ρU_ref² and the
projected frontal/plan areas. `p_ref` is the mean pressure over the upstream inlet slab
(freestream). It also emits a **per-cell mean Cp** field for visualisation. Sanity-validated:
windward stagnation **Cp ≈ 1.0** and bluff-body **Cd ≈ 1.2**. *Limitation → R1:* these are
**instantaneous** loads on an unsettled flow; skin friction is not included.

**W-DONE D — Interactive GUI Build workflow.** ✅
`src/gui/main_window.{cpp,h}` + `sim_worker.{cpp,h}` + `slice_viewer.{cpp,h}` implement the full
prepare-the-site → run loop:
- **Load a centerline STEP** (stored as `centerline_mesh_`), **gizmo-place** it (translate /
  rotate / scale) with a **30-level undo/redo** placement history (Ctrl+Z / Ctrl+Shift+Z), a live
  placement readout, and a reset.
- **Set the domain** (Lx/Ly/Lz + voxel h) with an explicit **Apply** — Build/Apply re-voxelize
  the placed house into the **current** domain (no surprise auto-resize).
- **Build** thickens + injects the building solid; a **Start Simulation** master gate holds the
  worker frozen (steps nothing) until pressed, so the site can be prepared with the flow static.
- **Live load readouts** (Cd/Cl/Cs + Cp range, refreshed on the repaint tick from the worker's
  published `WindLoads`) and **Cp colouring** of the building surface + a **legend**.
- Default scenario `configs/building.json` (a ~40×30×15 m wind tunnel at h = 0.25 m, air).

### 4b. REMAINING — in priority order

A milestone is DONE only when its acceptance check passes as a scripted headless smoke run /
probe and the commit lands. Never weaken an acceptance check to pass it — if one looks wrong,
flag it. Keep the `FluidCore` interface clean: load post-processing consumes only the fields the
core already exposes (`{u, v, w, p, ν_t, τ_wall}`).

**R1 — Convergence + time-averaging of the loads (+ radius sweep). [TOP PRIORITY]**
**What:** accumulate the loads over the run instead of reporting the latest step — detect the
statistically-steady window (e.g. monitor Cd/running-mean drift), then report **mean, RMS, and
peak** of Cp and of the force/moment coefficients over that window, and export the force/moment
**time series + summary to CSV**. Add a **batch runner** that sweeps `corner_radius` (and writes
one CSV row per radius) so a rounded-vs-sharp curve can be produced without hand-driving the GUI.
**Why:** the current live loads are **instantaneous on an unsettled flow**, so the headline
rounded-vs-sharp comparison is **not yet trustworthy**. Converged mean + peak statistics are the
deliverable a structural check consumes (peak, not mean, governs cladding + roof fixings).
**Acceptance:** on the building config, the running-mean Cd settles to within a stated tolerance
over the averaging window and is insensitive to the window start; mean/RMS/peak Cp + force
coefficients export to CSV and a re-run with a fixed seed reproduces them; the batch sweep emits a
radius→loads CSV in one command.

**R2 — ABL wind inflow profile (+ real building Re).**
**What:** assemble a proper atmospheric boundary-layer inlet from the existing pieces — a log-law
(or power-law) **mean** profile set by z₀ + reference speed/height, plus **SEM** turbulence at a
target intensity/length-scale (optionally precursor-seeded). Expose it as `inlet_profile: "abl"`
from config, wired into the building scenario (currently **uniform** top-hat inflow). Also
**retarget `configs/building.json`'s ν** from the gentle 1e-3 toward the real air value / building
Re (with the grid resolution that keeps it stable).
**Why:** design wind loads depend strongly on incoming shear and turbulence intensity; a top-hat
inlet under-loads and mislocates separation. This is the reference "design wind".
**Acceptance:** on an empty fetch, the mid-domain mean profile matches the target log-/power-law
within ~5% and the turbulence intensity holds at target across ≥ several flow-throughs without
divergence growth, from a fixed seed; the building run is stable at the retargeted ν.

**R3 — Validation against a wind-engineering benchmark (+ results reporting).**
**What:** reproduce a standard, well-documented bluff-body ABL benchmark for mean **and**
fluctuating Cp and integrated loads — e.g. the **Silsoe / surface-mounted cube**, the **CAARC**
standard tall building, or an **AIJ / TPU** aerodynamic-database case — matching that case's
specified inflow with the R2 ABL inlet, and produce an **HTML/CSV report** (Cp maps + load tables
with mean + peak).
**Why:** credibility of the rounded-vs-sharp numbers rests on the pipeline reproducing a case with
published experimental data; ideally validated **before or alongside** R4.
**Acceptance:** mean Cp and force coefficients fall within the benchmark's reported experimental
scatter at the specified grid, and the fluctuating/peak Cp trend is captured qualitatively; the
comparison + grid are documented in `configs/` and the report reproduces from a config.

**R4 — Rounded-vs-sharp comparison study.**
**What:** run the two W-DONE B variants (sharp r = 0 and rounded r > 0) under **matched** ABL
inflow, grid, and run length, then quantify the load differences — Cp distributions on facade +
roof, integrated force/uplift coefficients, and peak pressures — as a paired delta (rounded −
sharp). Now **unblocked** by the geometry + loads pipeline; **gated on R1** (convergence +
averaging) and best done after R3 (validation).
**Why:** this answers the project's question and produces the headline result (does rounding the
printed corners reduce drag / roof uplift / peak cladding pressure, and by how much).
**Acceptance:** a paired run identical except for corner radius yields signed, magnitude-bearing
load-difference metrics on **converged** statistics; both runs are individually reproducible; the
report states the effect of rounding on drag and roof uplift with an uncertainty estimate and the
r/h caveat (see §5).

**Sequencing.** R1 first — nothing else is trustworthy without converged loads. R2 supplies the
design wind; R3 validates the pipeline; R4 applies it to the real comparison. GUI/reporting
increments (mean/peak overlays, CSV export from the viewer, offscreen figures) can land
opportunistically once R1 produces averaged data.

---

## 5. Risks & open questions

1. **Rounded-corner physics is Reynolds-sensitive.** Sharp-edged bluff-body flow has *fixed*
   separation (Re-independent) — LES at model Re transfers to full scale. But a **rounded**
   corner's separation point *moves* with Re, and full-scale building Re is unreachable; LES at
   reduced Re can misplace separation on the rounded variant. Treat rounded-vs-sharp deltas as
   **model-scale** results, sanity-check against benchmark data (R3), and flag the caveat in any
   report — do not over-claim a full-scale drag reduction. The `building` config's gentle ν = 1e-3
   makes this acute until R2 retargets it.
2. **Voxel resolution vs corner radius.** A rounded corner of radius r is only meaningful when r
   spans **several voxels**; at h = 25 cm (building config) or 5 cm, a small radius is re-sharpened
   by the staircase and the comparison collapses. The distance-field thickening + arc fillet help,
   but the sharp/rounded distinction still needs the corner radius to resolve — likely the hero
   grid (h = 2.5 cm) for the final R4 numbers; state r/h for every result.
3. **Staircase over-prediction of drag.** The voxel staircase over-predicts Cd (measured ~+14% on
   the cylinder validation; cut-cell projection is deferred). Bias partly cancels in the
   *rounded − sharp difference* only if both run at **matched resolution** — enforce that, and
   report absolute coefficients with the staircase caveat.
4. **Loads are instantaneous until R1.** The current Cp/Cd/Cl are a single unsettled step — mean
   Cp settles quickly but peak/fluctuating pressures need long, converged windows, which drives
   wall-clock at hero resolution. Budget run length per the averaging requirement, and report the
   window used + convergence of the peak estimate. This is exactly what R1 fixes.
5. **Headless / offscreen GL for report figures.** The headless path uses Qt's `qoffscreen`
   platform, which has **no GL context**, so screenshot capture is unverified headlessly. R3/R4
   figure export may need an EGL/pbuffer offscreen GL context or a CPU raster fallback — confirm
   the capture path early rather than assuming a windowed screenshot works headlessly. (Cp maps can
   also be exported as VTI/CSV and rendered externally.)
6. **Validation-data availability & inflow match.** A benchmark is only usable if its **inflow
   (z₀, profile, TI, length scale) is fully specified** so the R2 ABL inlet can reproduce it. Pick
   the case (Silsoe cube / CAARC / AIJ / TPU) on that basis; some datasets report loads but
   under-specify the approach flow.
7. **Config schema & inherited naming.** The core `Config` uses a **strict key whitelist** and
   inherited names (`d50` for roughness); adding ABL roughness-length / reference-height / inlet-mode
   keys (R2) means extending both `config.cpp` and the GUI's looser `sim_setup` parser — keep the
   two in sync. Decide the ABL roughness parameterization (terrain-category z₀ vs the legacy
   z₀ = d50/12 mapping) in R2.
8. **Wall model on a vertical facade.** The log-law wall function was built for the horizontal
   ground; applying/validating it on the building's vertical walls and roof (where separation, not
   an attached log layer, dominates) needs checking — near separated regions the free-slip /
   resolved-BL treatment may be more appropriate than a log-law probe. `windloads` currently
   integrates **pressure only**, which sidesteps a wall-shear model for the loads themselves.

## 6. How implementation sessions should work

- Read `CLAUDE.md`, the target R-milestone, and the relevant `research/` notes / benchmark paper
  first. Physics constants come from those sources verbatim; if a needed constant is unspecified,
  **stop and flag it in the commit** rather than inventing one.
- Every new kernel gets a CPU reference and a GPU-vs-CPU parity test (≈1e-5 rel. for float32
  reductions); keep the `FluidCore` interface clean (load code only *reads* `{u,v,w,p,ν_t,τ}`).
- Acceptance checks are scripted headless smoke runs / probes (`building_probe`, `--offscreen`
  smoke flags) or a small validation harness; a milestone commit must show its check output.
  Determinism: fixed seeds for SEM eddies; document any nondeterminism.
- Commit style: Conventional Commits (`feat:`, `fix:`, `test:`, …). Keep the inherited fluid
  validations (cavity / cylinder) green — they are the foundation every wind result stands on.
