# CLAUDE.md — WindCFD

GPU (CUDA) 3D incompressible-flow **LES CFD** tool (C++/CUDA + Qt6) for **wind around
3D-printed concrete buildings**. Owner: MH (COBOD). Namespace `windcfd::core` (core) /
`windcfd::gui` (GUI).

## What WindCFD is for

Quantify how **wind** flows around a 3D-printed concrete building and the resulting wind
**loads on the facade and roof** — a bluff body in an atmospheric boundary layer (ABL). The
motivating question is a 3D-printing design choice: **rounded vs sharp building corners**, and
how that changes the flow field and surface loading.

The fluid is **air/wind** (use ρ ≈ 1.225 kg/m³, ν ≈ 1.5e-5 m²/s in prose/analysis). ⚠ The
code's `Config` defaults do NOT yet reflect this — see **Known gaps** below.

A CAD building (STEP) is voxelized into a solid obstacle; the LES solver develops the flow
around it; the Qt/OpenGL viewer shows live velocity/pressure slices, flow arrows and the model.
Surface pressure-coefficient (Cp) and integrated force/moment load extraction are **not built
yet** (Known gaps).

## Architecture / module map

Layering is deliberate so the physics core stays free of GL/Qt/OpenCascade. Three link
"islands": `libwindcfd` (pure CUDA+C++), `windcfd_geometry` (the only OCC target),
`windcfd_gui_cuda` (the only CUDA-GL-header target); `windcfd-gui` ties them to Qt.

### `libwindcfd` — the solver core (no Qt, no OCC, no GL)
- **Fluid solver — Stam stable-fluids on a MAC grid** (`src/core/fluid/`):
  `mac_grid.h`/`mac_ops.h` (face-staggered grid), `advect.cu` (MacCormack advection),
  `turbulence.cu` (Smagorinsky LES eddy viscosity), `project.cu` + `mgpcg.cu`
  (multigrid-preconditioned CG pressure projection), `stam_fluid_core.{h,cu}` (the core loop),
  `cavity.cpp` (lid-driven cavity driver). Velocity/pressure fields are **double**.
- **`fluid_core.h`** — the swappable `FluidCore` interface (`step()` → dt, `snapshot()`,
  `grid()`). One implementation today (`StamFluidCore` via `ChannelFluidCore`); an alternative
  core must stay drop-in.
- **Open-channel / obstacle flow** (`src/core/fluid/channel_*`): `channel_bc.h` (inlet /
  Orlanski convective outlet / lid / mask BCs, `flow_sign` direction, `solid_mode`
  no-slip/free-slip), `channel_ops.cu` (masked advection + diffusion + Poisson pieces),
  `channel_pressure.{h,cu}` (masked + Dirichlet-outlet MGPCG), `channel_mask.cpp` (cylinder
  mask builder), `channel_core.{h,cu}` (the `ChannelFluidCore`), `cylinder.cu` (drag CV + St
  FFT bluff-body diagnostics), `channel_porous.{h,cu}` (sub-grid porous momentum sink),
  `channel_empty.cpp`.
- **Near-wall / ground (ABL) model** — `bedshear.{h,cu}` + `bedshear_ops.h`: a **log-law wall
  function** (rough + transitional u\*, EMA time-filter, τ smoothing, rate limiting) used as the
  ground/ABL wall model for the LES. Kept as fluid infrastructure.
- **Turbulent inflow**: `channel_periodic.{h,cu}` + `periodic_ops.h` (periodic channel driver:
  body force + **PI mass-flux controller**, mixing-length ν_t), `sem_inlet.{h,cu}` (Jarrin 2006
  **synthetic eddy method** turbulent inlet, Cholesky-coloured Reynolds stresses),
  `channel_sem.cpp` (SEM open-channel driver), `precursor.{h,cpp}` (inlet-plane record/replay),
  `plane_ops.{h,cu}` + `inlet_fluct.h` (the SEM plane hook).
- **Geometry** (`src/core/geometry/`, OCC-free): `tri_mesh.h` (plain `TriMesh`),
  `voxelize.{h,cpp}` (**watertight ray-parity voxelizer** with K³-supersampled majority fill —
  turns a `TriMesh` into a solid obstacle mask; retains per-cell inside-fraction; brute-force,
  no BVH, runs once at load), `model_placement.h` (the one shared display+voxelize transform),
  `shape_masks.h` (procedural masks).
- **IO / config** (`src/core/`): `config.{h,cpp}` (JSON config via vendored nlohmann,
  `src/3rdparty/nlohmann/json.hpp`), `vti_writer.cpp`/`vti_reader.cpp` (VTI fields),
  `cuda_probe.cu`.

### `windcfd_geometry` — STEP import (the ONLY OpenCascade target)
`src/core/geometry/step_import.{h,cpp}`. `step_import.h` is OCC-free (returns a `TriMesh`); all
OCC includes live in the `.cpp`. Mirrors how `windcfd_gui_cuda` isolates GL — so `libwindcfd`
and all tests stay OCC-free. STEP traps handled in the `.cpp`: (a) OCC emits **millimetres** →
×0.001 to metres; (b) reversed-face **winding** fix for outward normals; (c) area-weighted
per-vertex normals.

### `windcfd_gui_cuda` — CUDA-GL interop (the ONLY CUDA target that sees GL headers)
`src/gui/slice_field.cu` (pure sampler + colourmap kernel with a `__host__ __device__`
evaluator → exact GPU-vs-CPU parity, no GL needed) + `slice_gl.cu` (the only GL-header TU:
`cudaGraphicsGLRegisterBuffer`, zero-copy writes into the registered colour VBO, on its own
non-blocking stream). Kept Qt-free so nvcc never sees Qt host flags.

### `windcfd-gui` — Qt6 + OpenGL 4.3 slice viewer
`src/gui/`: `main.cpp` (sets the GL 4.3 `QSurfaceFormat` before `QApplication`; CLI parsing),
`camera.h` (orbit/pan/zoom), `slice_viewer.{h,cpp}` (`QOpenGLWidget` — **all GL on the main
thread**, guarded by `WINDCFD_ASSERT_GL_THREAD()`; slice plane, model mesh, voxel staircase
overlay, flow arrows, clip plane, placement gizmo, legend, axis triad), `sim_worker.{h,cpp}`
(steps the `ChannelFluidCore` on a QThread, **no GL**, publishes device snapshots under a
display mutex; live setters for inlet speed/profile, flow direction, tidal reversal),
`sim_setup.{h,cpp}` (builds the sim from JSON; parses directly so viewer-only keys like U/Re are
allowed; `GridOverride` for the Apply-resize path), `scene_io.{h,cpp}` (Qt-free `.scn`
save/restore — self-contained binary container of recipe + mesh + live fields),
`flow_particles.{h,cpp}` / `flow_tracers.{h,cpp}` (Qt-free CPU tracer field drawn as an
instanced GL arrow glyph), `video_recorder.{h,cpp}` + `video_settings_dialog.{h,cpp}`
(fixed-cadence frame capture), `colormap.h` (the one shared colour ramp), `gl_thread_check.h`,
`main_stub.cpp` (`-DWINDCFD_ENABLE_QT=OFF` fallback for Qt-less hosts).

GUI features: STEP model load auto-voxelizes into the live obstacle (File → Open STEP… or
`--load-step`); left-side settings dock (live sim controls that don't reset, plus a Domain &
resolution Apply that rebuilds at t=0); flow arrows (2D/3D, density/speed/width); clip plane to
see inside hollow structures; model-placement gizmo (translate/rotate/scale, then voxelize where
placed); **flow reversal** (`flow_sign` ±1 / `set_flow_direction`, plus a live tidal-reversal
driver); scene checkpoints; recent files; video recording.

## Build

CMake + Ninja with the **VS 2022 Professional toolset (MSVC 14.44)**. ⚠ Do **not** use the
VS 2026 toolset — CUDA 13.1 rejects it. Put the VS 2022 `vcvars64` environment **and the
VS-bundled Ninja** on PATH first (the VS 2022 x64 dev prompt does both), then:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build --target windcfd-gui
```

- **Targets**: `libwindcfd` (core), `windcfd_geometry` (OCC STEP import), `windcfd_gui_cuda`
  (CUDA-GL interop), `windcfd-gui` (the app).
- **Dependencies**: CUDA 13.1, Qt 6.11.1 (`C:/Qt/6.11.1/msvc2022_64`, auto-detected via glob),
  OpenCascade 8.0 (`C:/OpenCASCADE-8.0/build2`).
- `CMAKE_CUDA_ARCHITECTURES=89` builds only for this Ada RTX 4090 (fast). The CMake default
  (`75;86;89`) is a multi-arch shipping build (~triples device-code compile time).
- **libwindcfd carries `${CUDAToolkit_INCLUDE_DIRS}` PUBLIC** so MSVC-compiled `.cpp` TUs find
  `cuda_runtime.h` (nvcc adds it implicitly for `.cu`, MSVC does not). Don't drop this.
- **DLL deployment**: a CMake POST_BUILD copies the Qt runtime (`windeployqt`) **and** the
  OpenCascade DLL closure next to `windcfd-gui.exe`, so it is double-clickable with no Qt/OCC on
  PATH. A missing OCC DLL (e.g. `freetype`) = `0xC0000135` at startup with no output.
- ⚠ MSVC 14.44 `CL.exe` occasionally crashes mid-compile on an AUTOMOC TU (`-1073741819` /
  0xC0000005). Transient — just re-run the build.

## Run

```
build/windcfd-gui.exe --config configs/g1_viewer.json                                    # windowed
build/windcfd-gui.exe --config configs/g1_viewer.json --offscreen --autoclose-ms 5000    # headless
```
- Windowed runs the interactive viewer. `configs/g1_viewer_full.json` is the full 10×10×5 m @
  h=5 cm grid (~4M cells); `g1_viewer_full_highres.json` is finer.
- **Headless `--offscreen` is solver-only**: the Qt offscreen platform has no GL context here,
  so it runs the CUDA solver but renders **0 frames** (`--autoclose-ms N` exits 0 after N ms).
  Useful for CI / solver debugging, not for visual checks.
- Other CLI flags (see `main.cpp`): `--load-step <path>`, `--model-place dx,dy,dz[,rz[,scale]]`,
  `--set-u <U>`, `--tidal Umax,plateau,ramp`, `--clip <x|y|z|camera>,<frac>`, `--tracers`,
  `--save-scene`/`--load-scene`, `--record`.

## Key conventions & invariants (preserve these)

- **Every CUDA kernel has a CPU reference + a GPU-vs-CPU parity test** at rel. max-norm **1e-5**
  (fluid kernels, MGPCG, channel kernels, the slice colour kernel — the last is exact).
- **Layering**: `libwindcfd` stays Qt-free, OCC-free and GL-free. OpenCascade lives only in
  `windcfd_geometry`; GL headers only in `windcfd_gui_cuda`'s `slice_gl.cu`. Keep new OCC/GL
  code inside those islands.
- **GL thread rule**: all OpenGL runs on the **main thread**; the solver runs on a worker
  QThread with no GL; hand-off is via queued Qt signals. Every viewer GL entry point asserts
  `WINDCFD_ASSERT_GL_THREAD()`.
- **Rendering is decoupled from stepping** (this is what holds ~60 fps): `step()` runs lock-free
  on the worker, which then copies {u,v,w,p} into device snapshot buffers under
  `display_mutex()`; the main-thread slice kernel reads only those snapshots, on the interop's
  own stream. Never lock across `step()` or route the render kernel onto the default stream, or
  fps collapses to the step rate.
- ⚠ **GL depth-state gotcha**: `QPainter::beginNativePainting()` resets GL to default state
  (depth test off, mask undefined). `paintGL` must re-establish depth test + `LESS` +
  `glDepthMask(TRUE)` every frame right before `glClear`, or the z-buffer silently dies.
- **Voxelizer must verify** enclosed volume vs mesh volume (< 2% mismatch). STEP geometry
  arrives in millimetres → ×0.001.
- **Physics constants are not invented from memory.** Values in the code cite an adversarially
  verified spec (`RESEARCH.md`, inherited — see history note). If a needed constant is missing,
  flag it; don't guess.
- **Style**: `.clang-format` = Allman + IndentBraces, tabs, ColumnLimit 0. Conventional Commits
  (`feat:`/`fix:`/`test:`/`docs:`/`refactor:`). One class per file, snake_case filenames, SI
  units everywhere (document a field's units at its declaration).

## Validation configs (`configs/`)

- **`v1_cavity_re100.json`** — lid-driven cavity, Re=100, vs Ghia (1982): verifies the core
  solver (MAC + MacCormack + MGPCG projection).
- **`v2_channel_loglaw.json`** — flat periodic channel: log-law wall function recovers u\*/κ,
  plus the Jarrin SEM synthetic-turbulence inlet. Exercises the ABL wall model + inflow.
- **`v3_cylinder.json`** / **`v3_tune32.json`** — flow past a cylinder, Cd/St (Kármán shedding)
  bluff-body validation. `v3_tune32` is the reduced fast-iteration variant.
- **`g1_viewer.json`** / **`g1_viewer_full.json`** / **`g1_viewer_full_highres.json`** — GUI
  viewer scenarios (reduced / full 4M-cell / high-res).
- **`m0_smoke.json`** — minimal smoke config.

## Known gaps / next

- **Air physics defaults not applied yet.** `Config` (`src/core/config.h`) still ships the
  inherited seawater values: `rho = 1027`, `nu = 1.36e-6`, plus a grain-size `d50` used only as
  a wall-roughness length. These need retargeting to air (ρ≈1.225, ν≈1.5e-5) and an ABL
  roughness for building-wind runs. Treat this as a documented follow-up.
- **No wind-load extraction yet.** Surface **Cp**, integrated **force/moment coefficients** on
  the facade+roof, and the **rounded-vs-sharp corner** comparison are not implemented.
- **No ABL wind inflow profile** wired into the building configs yet (the log-law wall model and
  SEM inlet exist as infrastructure).

## History

WindCFD is forked from a former marine/sediment simulator; all of that sediment, seabed and bed
-evolution simulation has been removed (`src/core/` now holds only `fluid` + `geometry`, and the
build no longer compiles any sediment sources). Companion docs — `RESEARCH.md` (physics &
numerics spec), `PLAN.md` (roadmap, milestones W1–W8), `HANDOVER.md` (current status) — were all
retargeted to WindCFD alongside this file.
