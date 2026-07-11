# PLAN.md — Implementation Plan: 3D Morphodynamic Scour-Protection Simulator

C++/CUDA simulator that ranks 3D-printed concrete shapes by their ability to make flowing
sand settle inside them and keep it there. Physics is fully specified in **RESEARCH.md**
(read its relevant section before implementing any milestone; constants come from there,
never from memory). Detailed literature notes live in `research/`.

This plan is written to be executed **one milestone per AI coding session** by cheaper
models. Each milestone is self-contained, lists exactly what to build, and ends with an
acceptance gate that must pass before the next milestone starts.

---

## 1. Target machine & toolchain (verified on this machine, 2026-07)

| Component | Version / Path |
|---|---|
| GPU | NVIDIA RTX 4090, 24 GB, driver 595.79 (CUDA driver API 13.2) |
| CUDA toolkit | **13.1** — `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.1`, `nvcc` on PATH. `CMAKE_CUDA_ARCHITECTURES=89` |
| CPU / RAM | i9-13900K, 128 GB |
| Compiler | MSVC 14.44 (VS2022 Professional), C++20 |
| CMake | 4.1.1 (`C:\Program Files\CMake\bin`). Ninja bundled with VS2022 |
| Qt | **6.11.1** — `C:/Qt/6.11.1/msvc2022_64` (also 6.10.1). Modules: Core, Gui, Widgets, OpenGLWidgets |
| OpenCascade | **8.0** — `C:/OpenCASCADE-8.0/build2` (inc/, win64/vc14/lib, win64/vc14/bin); 7.9 fallback at `C:/OpenCASCADE-7.9.0-vc14-64/opencascade-7.9.0` |
| TBB | 2021.13 — `C:/OpenCASCADE-7.9.0-vc14-64/3rdparty-vc14-64/tbb-2021.13.0-x64` (only if needed) |

Configure/build (mirrors cobod-slicer conventions):
```
cmake -B build -S . -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```
⚠ **Host-compiler trap:** this machine ALSO has VS 18 Community (MSVC 14.51), which CUDA 13.1
**rejects** ("unsupported Microsoft Visual Studio version" — only 2019–2022 supported), and
`vswhere -latest` resolves to it. Build from the **x64 Native Tools Command Prompt for VS 2022
(Professional)**, or pass
`-DCMAKE_CUDA_HOST_COMPILER="C:/Program Files/Microsoft Visual Studio/2022/Professional/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64/cl.exe"`
(and matching `CMAKE_CXX_COMPILER`). Never use `-allow-unsupported-compiler`. Ninja is not on
PATH — use the VS dev prompt or the VS-bundled ninja
(`…/2022/Professional/Common7/IDE/CommonExtensions/Microsoft/CMake/Ninja/ninja.exe`).

CMake: `project(WindCFD LANGUAGES CXX CUDA)`. Copy cobod-slicer's Qt auto-detect glob
(`C:/Qt/6.*/msvc2022_64`, newest with valid Qt6Config) rather than hard-coding a Qt version.
**Wire `enable_testing()` + ctest from M0** (the slicer team regrets not doing this).
Copy `.clang-format` from `C:\CODE\cobod-slicer` (team style: Allman + IndentBraces, tabs,
ColumnLimit 0). Conventional Commits (`feat:`, `fix:`, `test:`, …).

---

## 2. What we simulate

A 10×10×5 m patch of seabed with one printed unit or a small cluster (frontal area ≤ 3 m²
→ blockage ≤ 6%), uniform voxels h = 5 cm ranking / 2.5 cm hero. Modern monopiles (8–11.5 m
dia.) don't fit in this domain — near-pile conditions are represented by **amplified inflow
(1.5–2.0×) + 10–20% turbulence intensity** (RESEARCH.md §1, §8). Scenario matrix: U = 0.5 /
1.0 / 1.5 / 2.5 m/s × d50 = 0.2 / 0.35 / 1.0 / 5.0 mm (skip 5 mm @ 0.5 m/s — no motion).

Primary KPIs per shape: fitted equilibrium **trapped sand volume** inside the shape's control
volume, **retention** under a subsequent higher-flow stage, **edge scour** depth at the
perimeter, sinking/undermining volume beneath the unit.

---

## 3. Architecture

Three CMake targets; the core has no Qt dependency:

```
scour/
├── CMakeLists.txt
├── RESEARCH.md  PLAN.md  CLAUDE.md  research/
├── src/
│   ├── core/                 → libwindcfd (static lib, C++/CUDA, no Qt)
│   │   ├── config.h/.cpp         run config (JSON, vendored nlohmann single header)
│   │   ├── grid.h                MAC grid layout, indexing, SoA field arrays
│   │   ├── fields.cu/.h          allocation, host↔device, ping-pong buffers
│   │   ├── fluid/                FluidCore INTERFACE + StamFluidCore impl
│   │   │   ├── fluid_core.h          step(dt) → {u,v,w, ν_t, τ_b}; swappable (LBM plan-B)
│   │   │   ├── advect.cu             MacCormack semi-Lagrangian + clamp + solid reversion
│   │   │   ├── project.cu            MGPCG (Jacobi fallback for bring-up)
│   │   │   ├── turbulence.cu         Smagorinsky ν_t
│   │   │   ├── bedshear.cu           log-law wall function, EMA filter, cavity rules
│   │   │   │                         (produces τ_b for the Stam core; an LBM core supplies
│   │   │   │                          its own local wall shear natively)
│   │   │   └── boundary.cu           inlet log-law+SEM, Orlanski outlet, lid, mask BCs
│   │   ├── sediment/
│   │   │   ├── suspended.cu          c advection (w_s folded in), ν_t/σ_s diffusion, exchange
│   │   │   ├── bedload.cu            Wong-Parker / Engelund-Fredsøe on interface cells
│   │   │   ├── bedstate.cu           f_pack field, states, interface normals (RESEARCH §7)
│   │   │   └── avalanche.cu          layered-column 32/30° sweeps
│   │   ├── geometry/
│   │   │   ├── step_import.cpp       OpenCascade STEP→mesh (pattern below)
│   │   │   ├── stl_io.cpp            RWStl read/write
│   │   │   └── voxelize.cpp          watertight ray-parity voxelizer → PRINTED mask
│   │   ├── io/
│   │   │   ├── vti_writer.cpp        hand-rolled VTK ImageData (appended binary)
│   │   │   ├── metrics.cpp           KPI extraction, CSV timeseries
│   │   │   └── checkpoint.cpp        full-state save/restore
│   │   └── sim.h/.cpp                the operator-split loop (RESEARCH §2), MORFAC control
│   ├── cli/                  → scour.exe: run config(s), write VTI/CSV/report
│   └── gui/                  → windcfd-gui.exe: Qt6 live viewer (GL on main thread only)
├── tests/                    → GoogleTest; unit + validation gates as ctest labels
├── tools/                    → py scripts to plot CSV/compare gates (matplotlib, dev-only)
└── configs/                  → scenario matrix, benchmark cases (V1–V8)
```

**Memory budget** (fields: u,v,w ×2 ping-pong, p, div, c, f_pack, ν_t, τ_b, mask ≈ 15–20
floats/cell): 4M cells (h=5 cm) ≈ 0.3 GB; 32M cells (h=2.5 cm) ≈ 2.5 GB — trivial on 24 GB.
Keep everything resident; use 128 GB host RAM for async VTI/checkpoint buffering.

**STEP import** (copy the proven cobod-slicer pattern — see the exploration report facts):
`STEPControl_Reader` → `TransferRoots()` → `OneShape()` → `BRepMesh_IncrementalMesh(shape,
0.1)` → per-face `BRep_Tool::Triangulation`, applying `TopLoc_Location` transforms, 1-based→
0-based indices, and a winding fix (all three correctness-critical for watertight
voxelization; reference impl at cobod-slicer `src/app/widgets/render/mesh.cpp:494-522`,
`renderer.cpp:3223-3250`). ⚠ **Winding:** the slicer's `mesh.cpp:517` swaps indices when the
face is **NOT** reversed — its meshes are globally inverted vs the OCC convention. For our
voxelizer, swap when `face.Orientation() == TopAbs_REVERSED` to get outward normals; do NOT
copy the slicer's branch verbatim, and use **|signed mesh volume|** in the watertightness
check. **OCC outputs millimeters — multiply by 0.001.** Minimal OCC toolkits: `TKernel TKMath TKG2d TKG3d
TKGeomBase TKGeomAlgo TKBRep TKTopAlgo TKMesh TKXSBase TKDESTEP` (+`TKDESTL` for STL,
+`TKShHealing` for sewing dirty shells — lift `gentle_repair.h` from the slicer).

**GUI** (follow cobod-slicer's proven approach, not OCC's viewer): `QOpenGLWidget` +
`QOpenGLExtraFunctions`, **GL 4.3 core** requested via `QSurfaceFormat` *before*
`QApplication` (slicer uses 3.3; we need 4.3 for SSBOs/compute + clean CUDA-GL interop).
Reuse slicer patterns: `Camera` class (orbit/pan/zoom, ray picking), chunked `BufferArena` +
GLsync fence recycling, geometry-shader/instancing expansion for particles. CUDA writes
directly into registered GL buffers (`cudaGraphicsGLRegisterBuffer`) — zero copy. Sim runs
on a worker thread; GL strictly on the main thread; hand-off via
`QMetaObject::invokeMethod(..., Qt::QueuedConnection)` (slicer `mainwindow.cpp:1691-1703`
pattern). Reference files to read: slicer `src/main.cpp:474-483`, `opengl.h/.cpp`,
`camera.h`, `shaders/mesh.geom`.

---

## 4. Milestones

Rules: implement in order (GUI track G1/G2 can interleave). A milestone is DONE only when its
gate passes as a ctest and `git commit` lands. Print gate numbers in the test output.
Sizes: S ≈ half a session, M ≈ one session, L ≈ two sessions.
**Manual-gate exception:** G1/G2 interactivity gates (fps, GL-thread assert) are manual —
record measurements + screenshots in the milestone commit; the GL-main-thread check runs as
a Qt debug assertion in windcfd-gui, not ctest. M0's automated check is a `tools/` Python
reader validating the .vti header + payload; opening in ParaView is a one-time manual
confirmation. Each milestone commits the config(s) its gate consumes
(e.g. M1 → `configs/v1_cavity_re100.json`).

### M0 — Scaffold (S)
CMake (3 targets, CUDA arch 89, Qt auto-detect glob), GoogleTest + ctest wiring, config
loader, VTI writer + a synthetic-field smoke test, `.clang-format`, `.gitignore`, git init.
**Gate:** `ctest` green; a generated `.vti` opens in ParaView showing the synthetic field.

### M1 — Fluid core in a closed box (L)
MAC grid, MacCormack advection (with clamp + 1st-order reversion near solids), explicit
Smagorinsky, MGPCG projection (Jacobi fallback switch), no-slip/free-slip walls. Verbatim
specs: RESEARCH §3. Include a divergence-check kernel.
**Gate (V1):** lid-driven cavity, 128³ grid, **Cs = 0 (LES off)**, molecular ν = U·L/Re
treated explicitly, no wall model: RMS centerline error of u(y) and v(x) vs Ghia tables
(research/11; normalize by lid speed) **< 5% at Re=100**; Re=1000 is a diagnostic — log its
error (expect < 8% only at 128³+), do not hard-fail. Pressure: at production tolerance
‖r‖/‖b‖ ≤ 1e-4, MGPCG converges in ≤ 14 iterations and max|div u| < 1e-4·U/h; in validation
mode (tolerance 1e-6, iterations unconstrained) max|div u| < 1e-6·U/h. Closed box is
all-Neumann/singular: subtract mean(b) and pin p = 0 at one cell.

### M2 — Open channel + obstacles (M)
Inlet log-law profile, Orlanski-type outlet + flux balance + p pin, rigid lid, voxel obstacle
mask (hand-built cylinder first), free-slip solids. RESEARCH §3, §8.
**Gate (V3):** cylinder benchmark setup: D = 32 voxels, artificial ν = U·D/Re, domain
≥ 24D × 12D (thin spanwise, periodic or free-slip); measure St by FFT of cross-stream
velocity probed at (5D, 0) over ≥ 20 shedding cycles after 10 flow-throughs. Re=100:
St within 0.164–0.168 (±5%), Cd 1.33–1.40 (±10%); Re≈200 sheds vortices. Global mass
imbalance < 0.1%; a 10×10×5 m empty-channel run is steady and matches the inlet profile
within 5% at mid-domain.

### M3 — Bed shear + SEM inlet turbulence (M)
Log-law wall function with EMA filtering, transitional-roughness branch, stair-step
mitigations, cavity clearance rules (RESEARCH §4). Jarrin SEM at the inlet (RESEARCH §8).
**Also build here (needed by this gate and by V6/V7 later):** periodic x/y boundary option +
body-force channel driver g_x = u*²/h with a PI mass-flux controller
(g_x ← g_x + 0.1·(u*²/h)·(U_d − U_bulk)/U_d); precursor inlet-plane recording/replay
(5–10 Hz, ~300 s library) here or as the first task of M6.
**Gate (V2):** flat-bed periodic channel, h = 5 m, U = 1 m/s: extracted u* within 10% of the
analytic κ·U/(ln(h/z0) − 1) computed from the run's own h and z0 = d50/12 — i.e.
u* = 0.0345 m/s (τ_b ≈ 1.22 Pa) at d50 = 0.2 mm; repeat at d50 = 1.0 mm (u* = 0.040 m/s);
recovered κ from profile fit within 5%. ⚠ Do NOT target the rippled-bed Cd ≈ 0.0025 / u* =
0.050 (see RESEARCH §4). τ_b field smooth (no grid-pitch banding after mitigation); SEM run
maintains 5–10% ambient TI at mid-domain without divergence growth.

### G1 — GUI v1 (M, anytime after M1)
Qt window, GL 4.3 context, slicer-style camera, slice-plane rendering of any field via
CUDA-GL interop, play/pause/step, load config. No editing features.
**Gate:** 60 fps slice view of a running M2 simulation at h = 5 cm; no GL calls off the main
thread (assert via Qt debug).

### M4 — Suspended sediment (M)
Concentration field(s), advection with settling folded in, ν_t/σ_s diffusion, van Rijn
pickup with 2019 damping, w_s·c_b deposition as flux BC, hindered settling. RESEARCH §6.1–6.2.
**Gate (V5, V6):** settling column L2 < 2%, mass error < 0.1%. Rouse: periodic channel
h = 0.4 m, u* = 0.02 m/s, d50 = 0.1 mm (R = 0.91): fitted Rouse exponent of the steady
profile within 15% of R, with c Dirichlet-pinned to c_a at a = 0.05h — the Dirichlet pin is
permitted in this validation test ONLY; production bed exchange stays flux-based
(RESEARCH §6.2).

### M5 — Bed state + morphodynamics (L)
f_pack field, cell states, interface normals/areas, Winterwerp erosion, bedload donor
transfers (Wong-Parker default, Engelund-Fredsøe toggle), slope-corrected θ_cr, layered-column
avalanching, Exner-equivalence bookkeeping, MORFAC + cadence + limiters. RESEARCH §5–§7.
**Gate:** (a-i) still-water deposition: uniform c0 over a flat bed — column sum of f_pack
rises at w_s·c0 (bed elevation at w_s·c0/0.64), match < 0.5%, total mass error < 0.1%;
(a-ii) prescribed bedload q_b(x) = q0·sin(2πx/L) with E = D = 0: bed change matches
(1−p)·∂z_b/∂t = −∂q_b/∂x to < 0.5% after 100 updates (upwind ∇·q_b per RESEARCH §7);
(b) sand-pile relaxation settles to 30–32° everywhere; (c) three-reservoir mass
(bed + suspended + boundary fluxes) conserved to < 0.1% over 10⁴ steps; (d) no motion when
θ < θ_cr everywhere; onset within 10% of Soulsby–Whitehouse threshold.

### M6 — Pile-scour validation (L) — the big calibration gate
Roulund et al. (2005) benchmark at lab scale (Engelund–Fredsøe + suspended load ON).
**Setup:** voxel h = 1 cm (D/10); domain 3.0 × 1.6 × 0.4 m (pile D = 0.1 m centered 1.0 m
from inlet, ≥ 7D lateral clearance) ≈ 19M cells; V = 0.46 m/s, d50 = 0.26 mm, live-bed
V/Vcr = 1.25. Inflow from a precursor periodic channel at h_dom = 0.4 m; **verify recovered
u* = 0.020 m/s ± 10% on a frozen bed before enabling morphology** (research/14). Run ≥ 400 s
morphological time at M = 1 (T ≈ 130 s, research/11); fit S(t) = S_eq(1 − e^(−t/T)); the
gate compares fitted S_eq and T. Calibrate at most: pickup α (±30% band), Cs (0.10–0.12),
optional HSV shear multiplier — then FREEZE constants in `configs/calibrated.json`.
**Gate (V7):** upstream S/D = 1.25 ± 15%; downstream ± 30%; timescale within factor 2;
sensitivity run at 1.5× resolution changes S/D < 10%.

### M7 — Deposition validation "M-DEP" (M) — must pass before any ranking
Du et al. 2025 perforated-unit case (4–5 mm voxels) + current-only trench variant vs van
Rijn's trapping efficiency e_s (RESEARCH §10 V8; the full wave-stirred 1986 trench moves to
M10). **Build the C-AR/H-AR units procedurally as voxel masks** (no STEP import needed —
that's M8). **Move the thin-screen porous model here from M8:** Δp = ½·ρ·k·u_face²,
k = 1/β² − 1 (RESEARCH §8), applied to the perforated faces — at 4 mm voxels the 8.33 mm
holes are ~2 cells across, below the ≥ 8-cell resolution floor (research/13), so model them
sub-grid; do NOT attempt to resolve the hole jets.
**Gate (V8):** trapped volume ± 30%; depth ± 20–30%; BSS ≥ 0.3 (report if ≥ 0.6); measured
C-AR/H-AR ranking reproduced; e_s within ± 30%; qualitative deposition behind/inside the
porous array.

### M8 — Shape-testing harness (L, two sessions)
**Session 1 — geometry pipeline:** STEP/STL import + voxelization with unit-scale and
enclosed-volume gates (|signed mesh volume| vs voxel volume, fail > 2%); reference shapes
(solid ring, flat plate of equal footprint) generated procedurally as STL by
`tools/gen_ref_shapes.py` (no CAD input needed).
**Session 2 — campaign machinery:** scenario-matrix batch runner (run-length protocol of
RESEARCH §7: t ≥ 0.5T, dual-fit S_eq & V_eq, CIs, matched t/T comparisons); derive T,
θ/θ_cr regime, and the MORFAC schedule **at runtime from config (ρ, ν)** using RESEARCH §7
formulas — research/16 E6 (computed at ρ = 1025, ν = 1.05e-6) is a reference example, not a
lookup table; matrix configs generated by `tools/gen_matrix.py`. KPI extraction (trapped
volume in shape CV, retention stage, edge scour, undermining), HTML/CSV comparison report.
`checkpoint.cpp` (full-state save/restore; unit test: save→restore→step bit-identical to
uninterrupted stepping). **Retention protocol:** from the fitted-equilibrium state
(checkpoint), ramp U to the next-higher matrix velocity over 60 s (cosine), run 1·T at the
new U, report retained fraction V_trapped(end)/V_trapped(eq).
**Gate:** end-to-end: the two reference shapes rank correctly with non-overlapping CIs; a
full 4M-cell scenario completes in ≤ 25 min wall-clock; re-running a config is
bit-reproducible.

### G2 — GUI v2 (M)
Volume rendering of c (raymarch slices), bed surface extraction (marching cubes on F = 1
or GPU point splat), particle tracers (instanced billboards), KPI live plots, side-by-side
run comparison.
**Gate:** interactive (≥ 30 fps) on a live 4M-cell run; screenshots exportable for reports.

### G3 — Preparation phase: drop & settle rigid bodies (Jolt) (L) — DONE 2026-07-10
**Status: COMPLETE.** G3.1 (settle core, `scour_physics_tests` gate), G3.2 (windowed drop demo,
MH sign-off — "perfect"), G3.3a (`resettleOnBed` unit-gated), G3.3b (live sink coupling: headless
verified — re-settle fires on bed change, pile sinks, flow re-masks in place, 1600+ steps stable).
Jolt v5.5.0 (prebuilt `/MD`, `find_package(Jolt CONFIG)`); 100/100 unit; gate_G1 PASS. Windowed:
`windcfd-gui --config configs/g1_viewer_full.json --drop-step <XStone.stp> --drop <count>,<seed>[,smin,smax,sand]`.

A pre-simulation "experiment setup" stage: scatter N scaled copies of a printed protection unit
(e.g. Holcim **XStone**, loaded from STEP) above the bed and let them fall + settle under gravity
with body–body + body–bed contact, then voxelize the settled pile as the rigid structure and start
the morphodynamic run. Rigid-body settling is **fully DECOUPLED from the CUDA fluid** — it runs on
the CPU (**Jolt Physics**, MIT) and emits one affine `ModelPlacement` per block, consumed unchanged
by the existing voxelizer + `build_seabed_from_structure` path. Holes are preserved: the settling
**collision proxy** is a convex compound, but the flow still voxelizes the TRUE perforated mesh.
**Isolation:** Jolt is linked only by a new `scour_physics` static lib (mirrors how `windcfd_geometry`
isolates OpenCascade); libwindcfd + every physics gate stay Jolt-free. Bodies persist for the whole
run so the sink-coupling increment is additive.
- **G3.1 — settle core.** `scour_physics` (Jolt) + `load_step_solids()` multi-solid STEP loader
  (each convex piece → a Jolt `ConvexHullShape`; bootstrap = whole-mesh hull). Scatter → settle on a
  rigid bed → `std::vector<ModelPlacement>`. Deterministic (fixed seed).
- **G3.2 — windowed self-driving demo.** A non-interactive-but-VISIBLE scenario that scatters,
  animates the live settle (repeat-draw the mesh at each body transform), commits the union voxel
  mask, adds sand, and runs the flow — MH watches + narrates. CLI `--drop <shape>,<count>,<seed>`.
- **G3.3 — two-way sink coupling.** Persistent Jolt world; when the bed erodes ≥ ~½ cell the ground
  heightfield (from z_b) is rebuilt, blocks re-settle into the scour hole, moved blocks are
  re-voxelized, and the flow re-masks in place (`update_solid`) with neighbour-averaged seeding of
  newly-revealed cells + a dynamic structure-footprint pin.

**Gate (G3.1, automated):** drop ≥ 8 blocks, settle to sleep within the step budget; every body
rests with min-z ≥ bed and no pair interpenetrates beyond a small tolerance; re-running with the
same seed is bit-reproducible. **Gate (G3.2/G3.3, manual + stability):** the committed pile
voxelizes watertight, the flow diverts and stays stable, sand mass conserves to < 1e-5 over the run;
the windowed demo renders the drop + scour without artifacts (MH sign-off).

### M9 — Tidal reversal + ranking campaign (M)
Face-swap reversal with cosine ramp (RESEARCH §8), M ≤ 5 in reversing scenarios; run the
first real shape-candidate campaign; write the comparison report template.
**Gate:** reversal run conserves total sand mass (bed + suspended ± boundary fluxes) to
< 0.1% and shows no monotonic kinetic-energy growth across 4 slack transitions; MORFAC
limiter clip-fraction < 1% (research/16); first ranked report of ≥ 3 candidate shapes
delivered.

### M10 — Waves (L, later)
Oscillatory inlet forcing (KC-parameterized), wave-current combination, backfill benchmark,
plus the **full van Rijn 1986 trench** (waves H = 0.08 m, T = 1.5 s over 0.18 m/s current,
sand feed 0.0167 kg/s/m — deferred from M7 because its infill is wave-stirred).
**Gate (V9):** Sumer 2013 backfill direction & final-KC equilibrium reproduced qualitatively;
S/D vs KC trend matches 1.3·(1 − e^(−0.03(KC−6))) within ±30% at 2–3 KC points; trench
infill volume ±30%, BSS ≥ 0.3.

### M11 — Performance & hero runs (M, optional)
Profile (Nsight), kernel fusion, 2.5 cm hero-grid runs of the top-3 shapes, longer-horizon
storm-sequence scenarios.
**Gate:** ≤ 150 ms/step at 32M cells; hero-vs-ranking grid ordering unchanged for top-3.

---

## 5. Performance budget (RTX 4090, from research/16 traffic model)

| Grid | Cells | ms/step | 300 s flow | 20 min wall-clock buys |
|---|---|---|---|---|
| h = 5 cm | 4.0M | ~12–18 | 2–4 min | 1 400–7 000 fluid-seconds (U-dependent dt 18–88 ms) |
| h = 2.5 cm | 32M | ~100–150 | ~1.5–2 h at hero settings (CFL 1.0, max\|u\| = 1.7·U∞ → dt ≈ 6–30 ms); ~30 min only at CFL 2 / dt = 20 ms (lower bound) | hero runs only |

Every scenario-matrix cell closes in 15–25 min at h = 5 cm using the MORFAC schedule
(M = 6–9 only for U = 0.5 m/s; M = 1–2 otherwise). If a run exceeds budget, check the
adaptive limiter M·N·dt·max|dz_b rate| ≤ 0.05·h before touching anything else.

## 6. Risks & mitigations

1. **HSV under-resolution → under-predicted upstream scour.** Expected; handled by M6
   calibration + honesty rules (RESEARCH §11). If calibration cannot reach the V7 gate,
   promote the LBM plan-B core behind `FluidCore` (research/09; do NOT copy FluidX3D code —
   non-commercial license).
2. **Ranking artifacts from fixed run lengths.** Prevented by the fitted-equilibrium protocol;
   never compare shapes at different t/T.
3. **Stair-step τ_b noise → noisy erosion.** Mitigations in RESEARCH §4; verify in M3 gate.
4. **STEP units/watertightness.** mm→m factor; gentle-repair + sewing before voxelization;
   voxelizer must report enclosed volume vs mesh volume (fail > 2% mismatch).
5. **Two-branch formulas silently truncated** (van Rijn T ≥ 3, transitional roughness,
   slope-correction clamps). Each has a dedicated unit test in its milestone.

## 7. How implementation sessions should work

- Read `CLAUDE.md`, this plan's milestone, and the referenced RESEARCH.md sections first.
- Constants come from RESEARCH.md verbatim; if code needs a constant not specified there,
  stop and flag it in the commit message rather than inventing one.
- Every kernel gets a CPU reference implementation and a GPU-vs-CPU unit test (tolerance
  1e-5 rel. for float32 reductions).
- Gates are ctest targets labeled `gate_M<k>`; a milestone PR/commit must show its gate
  output. Never weaken a gate to pass it — if a gate seems wrong, flag it.
- Determinism: fixed seeds (SEM eddies), no atomics in accumulation paths where avoidable;
  document any nondeterminism.
- Keep `FluidCore` interface pure: sediment code may only consume `{u, v, w, ν_t, τ_b}`.
  For the Stam core, τ_b comes from the log-law probe (`fluid/bedshear.cu`); an LBM core
  supplies its own local wall shear natively.
