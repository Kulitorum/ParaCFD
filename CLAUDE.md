# CLAUDE.md — WindCFD

GPU (CUDA) 3D incompressible-flow **LES CFD** tool (C++/CUDA + Qt6) for **wind around
3D-printed concrete buildings**. Owner: MH (COBOD). Namespace `windcfd::core` (core) /
`windcfd::gui` (GUI).

## What WindCFD is for

Quantify how **wind** flows around a 3D-printed concrete building and the resulting wind
**loads on the facade and roof** — a bluff body in an atmospheric boundary layer (ABL). The
motivating question is a 3D-printing design choice: **rounded vs sharp building corners**, and
how that changes the flow field and surface loading.

The fluid is **air/wind** (ρ ≈ 1.225 kg/m³, ν ≈ 1.5e-5 m²/s). `Config` defaults now reflect
this (air; `src/core/config.h`).

A 3D-printing **centerline** STEP is loaded, placed with a gizmo, and voxelized into a solid
building (thickened walls + flat roof); the LES solver develops the flow around it; the
Qt/OpenGL viewer shows live velocity/pressure slices, flow arrows, the model, and the building
**Cp-coloured** by surface pressure. Surface pressure-coefficient (Cp) and integrated
force/moment **wind-load** extraction (Cd/Cl/Cs + moments) are **built** — see the
`windloads` module. Remaining caveats are in **Known gaps** (loads are still instantaneous on
an unsettled flow; inflow is uniform, not ABL).

## Architecture / module map

Layering is deliberate so the physics core stays free of GL/Qt/OpenCascade. Three link
"islands": `libwindcfd` (pure CUDA+C++), `windcfd_geometry` (the only OCC target),
`windcfd_gui_cuda` (the only CUDA-GL-header target); `windcfd-gui` ties them to Qt.

### `libwindcfd` — the solver core (no Qt, no OCC, no GL)
- **Graded structured grid** (`src/core/fluid/grid_metrics.{h,cpp}`): `MacGrid` (`mac_grid.h`) carries
  nullable per-axis metric arrays (cell widths + centres + cumulative face coords) with accessors that
  fall back to the exact uniform `h` when null (uniform = byte-identical). `FineCoreSpec` +
  `GridMetrics::generate` build a UNIFORM `h_fine` core inside a geometrically graded coarse far field;
  `GridMetrics` owns host+device copies and hands out `host_view()` (CPU/geometry) vs `device_view()`
  (GPU kernels). Shared world↔index map (`locate_frac`/`grid_fx/cx/…`). See the graded-grid note in
  **Known gaps**.
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
  `channel_pressure.{h,cu}` (masked + Dirichlet-outlet MGPCG), `channel_core.{h,cu}` (the
  `ChannelFluidCore`), `channel_porous.{h,cu}` (sub-grid porous momentum sink),
  `channel_empty.{h,cpp}` (the empty-channel obstacle-free start). The obstacle mask is now
  supplied by the **building** pipeline (below) or an arbitrary voxelized `TriMesh` — the old
  procedural default **cylinder obstacle was removed** (`channel_mask.cpp` cylinder-mask builder
  and `cylinder.cu` Cd/St bluff-body diagnostics are deleted); a fresh sim starts as an **empty
  channel** until a building/model is injected.
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
  turns a watertight `TriMesh` into a solid obstacle mask; retains per-cell inside-fraction;
  brute-force, no BVH, runs once at load), `model_placement.h` (the one shared display+voxelize
  transform — `ModelPlacement` affine `world(v)=M·v+t`; `placed_mesh()` applies a gizmo
  placement to every vertex of a `TriMesh` so section/voxelize land under the drawn model),
  `shape_masks.h` (procedural masks, header-only).
- **Building pipeline** (`src/core/geometry/building.{h,cpp}`) — turns a **3D-printing
  CENTERLINE** STEP (vertical wall ribbons, loaded as a `TriMesh`) into a solid building mask
  **without** a watertight solid, via a distance-field voxelization:
  - `mesh_horizontal_section(mesh, z0)` → intersects the mesh with the horizontal plane `z=z0`
    (`NaN` ⇒ robust mid-height), chains section segments into closed loops → a 2D `Footprint`
    (closed polylines in domain xy, metres).
  - `center_footprint(fp, Lx, Ly)` → translate a footprint's xy-bbox centre to the domain
    centre (used by `building_probe`; the GUI places via `placed_mesh` instead, so no centring).
  - `voxelize_building(fp, prm, g)` → the mask (1=solid/0=fluid, `ChannelBC` cell field,
    `g.pidx`). A cell is **solid** when it lies within **± half `wall_thickness`** of the
    (optionally corner-rounded) footprint over the wall height `[base_z, base_z+wall_height]` —
    the half-width is **floored to one cell's circumradius `0.5·h·√2`** so a wall thinner than
    ~1 cell still voxelizes **watertight** (a printed 50–80 mm COBOD wall stays a solid band even
    at `wall_thickness ≈ h`; walls already ≥ ~1.4 cells keep their exact width, `max()` is a
    no-op) — **plus** a **SOLID convex-hull flat roof slab** over `[top, top+roof_thickness]` dilated by
    `roof_overhang` (hull guarantees a filled roof regardless of how the section split into
    loops). **Rounded** corners: a 2D fillet pre-rounds the centerline (`round_loop`) so the
    outer radius ≈ `corner_radius`. **TRUE sharp** corners (`corner_radius`≈0): start from the
    ±half distance band (which is naturally round) and square off each **convex** corner by
    adding an outward **miter WEDGE** triangle (local + robust on non-convex outlines,
    miter-limited to bevel spikes); curved walls stay rounded via an angle threshold.
  - `BuildingParams`: `wall_thickness` (default **0.08 m** — a COBOD-printed wall), `wall_height`,
    `corner_radius`, `roof_overhang`, `roof_thickness`, `base_z` (all metres).
  - **`tools/building_probe`** (dev CLI, `windcfd_geometry`-linked) — load a centerline STEP →
    section → voxelize headlessly and print stats, to verify the pipeline before the GUI.
- **Wind loads** (`src/core/windloads.{h,cpp}`) — integrate the **pressure** load on the
  voxelized building. `compute_wind_loads(p, solid, g, prm, out_cell_cp)` walks the exposed
  voxel faces (a face between a solid cell and an adjacent fluid cell), takes the surface
  pressure as the first fluid cell's `p` (cell-centred, **physical Pa**), and sums
  `F = -Σ (p_face − p_ref)·n̂·h²` (n̂ = outward solid→fluid normal). `p_ref` = **mean pressure
  over the upstream inlet slab** (freestream). Coefficients use `q = ½ρU²`:
  `Cp=(p−p_ref)/q`, `Cd=Fx/(qA_frontal)`, `Cs=Fy/(qA_frontal)`, `Cl=Fz/(qA_plan)`, plus moments
  `Mx/My/Mz` + `CM*`, `A_frontal`/`A_plan`, `L_ref`=height, `cp_min/max`, and an optional
  **per-cell Cp** field for visualisation. Pressure-only (no skin friction); wind assumed +x.
  Host-only, OCC-free.
  - **`LoadAverager` + `WindLoadStats`** (same file) — turn the **instantaneous** per-step
    `WindLoads` into **converged, time-averaged** statistics. A turbulent (LES) flow never settles
    instantaneously — only its statistics do — so a single snapshot is one random draw. `add()`
    folds each sample (Welford, numerically stable) into a running **mean + RMS(fluctuation) +
    peak(min/max)** of Cd/Cl/Cs plus a running **time-average of the per-cell surface Cp**;
    `result()` reports those + the averaged-Cp range and a **convergence-drift hint** (how far the
    running mean still moves between the window's first half and the whole — small ⇒ converged).
    The window is **user-controlled** (there is deliberately **no automatic convergence
    detection** — the user watches the flow spin up and judges by eye). Host-only.
- **IO / config** (`src/core/`): `config.{h,cpp}` (JSON config via vendored nlohmann,
  `src/3rdparty/nlohmann/json.hpp`; **air defaults** ρ=1.225, ν=1.5e-5, `d50` retained only as a
  wall-roughness length), `vti_writer.cpp`/`vti_reader.cpp` (VTI fields), `cuda_probe.cu`.

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
overlay (Cp-coloured), flow arrows, clip plane, placement gizmo, legend, axis triad),
`main_window.{h,cpp}` (the dock UI: building workflow, gizmo + placement undo/redo, Start gate,
wind-load readout), `sim_worker.{h,cpp}` (steps the `ChannelFluidCore` on a QThread, **no GL**,
publishes device snapshots under a display mutex; integrates + publishes `WindLoads`; live
setters for inlet speed/profile, flow direction, tidal reversal; a **run gate** so it can be
held paused), `sim_setup.{h,cpp}` (builds the sim from JSON; parses directly so viewer-only keys
like U/Re are allowed; `GridOverride` for the Apply-resize path), `scene_io.{h,cpp}` (Qt-free
`.scn` save/restore — self-contained binary container of recipe + the **source STEP** (embedded
as the source of truth; the display mesh is **regenerated** from it on load via
`load_step_mesh_from_memory`, not stored — legacy mesh-blob scenes still load) + live fields; a
`mesh_is_centerline` flag routes a restored model back to the centerline pipeline so **Build works
after a scene load** instead of "no centerline loaded"),
`flow_particles.{h,cpp}` / `flow_tracers.{h,cpp}` (Qt-free CPU tracer field drawn as an
instanced GL arrow glyph), `video_recorder.{h,cpp}` + `video_settings_dialog.{h,cpp}`
(fixed-cadence frame capture), `colormap.h` (the one shared colour ramp), `gl_thread_check.h`,
`main_stub.cpp` (`-DWINDCFD_ENABLE_QT=OFF` fallback for Qt-less hosts).

GUI features:
- **Building workflow** (the primary path): File → **Open centerline STEP…** or
  `--load-centerline` loads a 3D-printing centerline. Dial the wall/roof params in the
  **"Building (centerline → solid)"** dock group, then **Build** → voxelizes the **placed**
  centerline into the **CURRENT** domain (section → footprint → `voxelize_building`) and injects
  it as the obstacle (worker rebuild off-thread). **No auto-resize** — the user sets the domain
  size + voxel size and presses **Apply** (which re-voxelizes the placed house on the new grid).
- **Model placement gizmo**: the loaded centerline (or a plain STEP model) is gizmo-placeable
  (translate/rotate/scale); a **placement readout** shows the live pose; **30-level undo/redo**
  (**Ctrl+Z** / Ctrl+Shift+Z / Ctrl+Y, or dock buttons) makes accidental gizmo edits reversible;
  **Reset placement** returns to centre-on-bed. Gizmo edits only MOVE the model — press
  Build/Apply to re-voxelize the new pose.
- **"Start Simulation" gate**: the worker starts **held (paused)** so the user can prepare the
  site (load, place, Build, set domain) with the flow frozen; pressing **Start** releases it.
  Headless/scripted runs **auto-start** after setup. Display/repaint and Build/voxelize work
  while held.
- **Live wind loads**: a dock readout shows live **Cd / Cl / Cs** and the **Cp** range,
  integrated worker-side from the pressure field over the building surface.
- **Converged time-averaged loads**: those live loads are **instantaneous on a turbulent flow**
  (one random sample). A **"Start averaging"** button begins accumulating converged **mean / RMS
  / peak** Cd/Cl/Cs + a **time-averaged surface Cp** over a user-controlled window; a second
  readout shows the state, sample count, elapsed sim-time in **flow-through times** (Lₓ/U), the
  converged coefficients and a **convergence-drift** hint. **"Start in N h"** defers the start by
  **wall-clock** time for **overnight runs** (start the sim, schedule +6 h, wake to a settled
  multi-hour average). **Stop** freezes the result. No automatic convergence detection — the user
  watches the instantaneous Cd stop trending, then starts the window (see `LoadAverager`).
- **Cp colouring**: the building's voxel staircase is tinted by surface **Cp** (a diverging
  colour ramp) with a legend + a "Colour building by Cp" toggle; once averaging is active the
  tint + legend switch to the **time-averaged** Cp (the meaningful field for comparing suctions).
- Plain-STEP load auto-voxelizes into the live obstacle (File → Open STEP… / `--load-step`);
  left-side settings dock (live sim controls that don't reset; **Input speed U** capped at
  **40 m/s**; a Domain & resolution **Apply** that rebuilds at t=0); flow arrows (2D/3D,
  density/speed/width); clip plane to see inside hollow structures; **flow reversal**
  (`flow_sign` ±1 / `set_flow_direction`, plus a live tidal-reversal driver); scene checkpoints;
  recent files; video recording.

## Build

CMake + Ninja with the **VS 2022 Professional toolset (MSVC 14.44)**. ⚠ Do **not** use the
VS 2026 toolset — CUDA 13.1 rejects it. Put the VS 2022 `vcvars64` environment **and the
VS-bundled Ninja** on PATH first (the VS 2022 x64 dev prompt does both), then:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build --target windcfd-gui
```

- **Targets**: `libwindcfd` (core), `windcfd_geometry` (OCC STEP import), `windcfd_gui_cuda`
  (CUDA-GL interop), `windcfd-gui` (the app), `building_probe` (dev CLI: centerline → mask),
  `parity_probe` (the GPU-vs-CPU + golden parity oracle; `ctest -R parity` / `./build/parity_probe.exe`).
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
build/windcfd-gui.exe --config configs/building.json                                     # windowed building sim
build/windcfd-gui.exe --load-centerline house.stp                                        # loads configs/building.json
build/windcfd-gui.exe --config configs/g1_viewer.json --offscreen --autoclose-ms 5000    # headless (solver only)
```
- **Interactive building workflow**: load a centerline (`--load-centerline` / File ▸ Open
  centerline STEP…) → **gizmo place** it (translate/rotate/scale, undoable) → set the **domain
  size + voxel size** (Apply) → **Build** (voxelize the placed house into the current domain) →
  **Start Simulation** → read the live **Cd/Cl/Cs + Cp** and the Cp-coloured building.
- **`configs/building.json`** — the default air wind-tunnel (~40×30×15 m @ h=0.25 m, ~1.15M
  cells, ρ=1.225). Auto-selected when `--load-centerline` is given with no `--config`.
- Windowed runs the interactive viewer. `configs/g1_viewer*` are **empty-channel** demos (no
  built-in obstacle now that the cylinder was removed); `g1_viewer_full.json` is the full
  10×10×5 m @ h=5 cm grid (~4M cells), `g1_viewer_full_highres.json` finer.
- **Headless `--offscreen` is solver-only**: the Qt offscreen platform has no GL context here,
  so it runs the CUDA solver but renders **0 frames** (`--autoclose-ms N` exits 0 after N ms);
  the worker auto-starts (no Start button). Useful for CI / solver debugging, not visual checks.
- Other CLI flags (see `main.cpp`): `--load-centerline <path>`, `--load-step <path>`,
  `--model-place dx,dy,dz[,rz[,scale]]`, `--set-u <U>`, `--tidal Umax,plateau,ramp`,
  `--clip <x|y|z|camera>,<frac>`, `--tracers`, `--save-scene`/`--load-scene`, `--record`,
  `--average-now` / `--average-in <s>` (headless smoke of load-averaging: start/schedule the
  converged window; the exit-log prints the mean/RMS/peak + flow-through count).

## Key conventions & invariants (preserve these)

- **Every CUDA kernel keeps a matching CPU reference twin** (`*_cpu`, e.g. `ch_poisson_apply_cpu`)
  as the design contract, exercised by the **live parity oracle** `tools/parity_probe.cu` — a
  `ctest` target (`add_test(NAME parity)`, no gtest dependency) that for every MAC kernel checks
  (1) GPU-vs-CPU-twin parity ≤ **1e-5** AND (2) GPU-vs-blessed-golden ≤1e-5 (the golden gate is the
  permanent "grading collapses to uniform" regression). Blessed baseline in `tests/golden/*.f64`;
  re-bless with `parity_probe --bless`. It also covers the graded path directly (GPU-vs-CPU on a
  generated graded grid, MMS Laplacian order, graded MGPCG convergence, the windloads bbox guard,
  voxelizer world→index, and scene round-trip). Run: `ctest --test-dir build -R parity` or
  `./build/parity_probe.exe`. Build it via the VS dev env (see Build). Preserve the `_cpu` twins and
  keep this green across any solver change.
- **Layering**: `libwindcfd` stays Qt-free, OCC-free and GL-free. OpenCascade lives only in
  `windcfd_geometry`; GL headers only in `windcfd_gui_cuda`'s `slice_gl.cu`. Keep new OCC/GL
  code inside those islands. (`building.*` and `windloads.*` are host-only + OCC-free.)
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
- **Voxelizer must verify** enclosed volume vs mesh volume (< 2% mismatch) for watertight solids
  (the building pipeline uses a distance test instead — inherently leak-proof). STEP geometry
  arrives in millimetres → ×0.001.
- **Model placement is single-sourced**: `model_placement.h` is the ONE display+voxelize
  transform; the viewer draws and `placed_mesh()`/the voxelizer section at the same pose.
- **Physics constants are not invented from memory.** Values in the code cite an adversarially
  verified spec (`RESEARCH.md`, inherited — see history note). If a needed constant is missing,
  flag it; don't guess.
- **Style**: `.clang-format` = Allman + IndentBraces, tabs, ColumnLimit 0. Conventional Commits
  (`feat:`/`fix:`/`test:`/`docs:`/`refactor:`). One class per file, snake_case filenames, SI
  units everywhere (document a field's units at its declaration).

## Validation / demo configs (`configs/`)

- **`v1_cavity_re100.json`** — lid-driven cavity, Re=100, vs Ghia (1982): verifies the core
  solver (MAC + MacCormack + MGPCG projection).
- **`v2_channel_loglaw.json`** — flat periodic channel: log-law wall function recovers u\*/κ,
  plus the Jarrin SEM synthetic-turbulence inlet. Exercises the ABL wall model + inflow.
- **`building.json`** — the air wind-tunnel default for the centerline Build workflow
  (~40×30×15 m @ h=0.25 m, ρ=1.225); auto-selected by `--load-centerline`.
- **`g1_viewer.json`** / **`g1_viewer_full.json`** / **`g1_viewer_full_highres.json`** — GUI
  **empty-channel** demos (reduced / full 4M-cell / high-res), no built-in obstacle.
- **`m0_smoke.json`** — minimal smoke config.
- The `v3_cylinder*` bluff-body validation configs were **removed** with the cylinder obstacle.

## Known gaps / next

- **DONE**: air physics defaults (`Config` now ships ρ=1.225, ν=1.5e-5); the **building
  pipeline** (centerline → thickened walls + roof solid, rounded vs sharp corners); the
  **wind-load extraction** (surface **Cp**, integrated **Cd/Cl/Cs** + moments), live in the GUI.
- **DONE — convergence + TIME-AVERAGED loads** (`LoadAverager`/`WindLoadStats`, GUI "Start
  averaging" / "Start in N h" / Stop): the live loads are still instantaneous, but the user can
  now accumulate **converged mean / RMS / peak** Cd/Cl/Cs + a **time-averaged Cp** over a window
  they control (with an overnight wall-clock deferred start), and the building recolours to the
  averaged Cp. **Remaining judgement call**: convergence is **not auto-detected** — the user must
  run long enough (watch the drift hint / the mean plateau, ~10+ flow-through times) before the
  **rounded-vs-sharp corner** comparison is trustworthy. Verify each config's mean is converged
  (small drift, tight RMS) before comparing two runs.
- **No ABL wind inflow profile** — the inlet is currently **uniform**; the log-law wall model
  and SEM inlet exist as infrastructure but aren't wired into the building configs.
- **Building viscosity is a moderate `nu = 1e-3` default** (`configs/building.json`), not the
  true air Re — a numerically-forgiving value, not a physically-resolved one.
- **DONE (mostly) — resolution decouple via the graded structured grid.** The uniform grid welded
  feature resolution to domain size (`h=0.25 m` ⇒ a 0.30 m corner was ~1.2 cells, so rounded ≈ sharp
  and the effect the tool exists to measure was invisible). Fixed with an axis-separable **graded
  structured grid**: a UNIFORM `h_fine` core around the building inside a geometrically graded coarse
  far field (`src/core/fluid/grid_metrics.{h,cpp}` — `FineCoreSpec` + `GridMetrics::generate`).
  `MacGrid` now carries nullable per-axis metric arrays (`dx/dy/dz`, `xc/yc/zc`, `xf/yf/zf`) with
  accessors that fall back to the exact uniform closed form when null, so a uniform grid is
  byte-identical. Every operator, the voxelizers, windloads, the slice sampler, and the flow tracers
  are metric-aware; `adaptive_dt` keys off `h_min()`. Opt-in via a config `"fine_core": {enabled,
  x0..z1, h_fine, growth}` object (parsed in `sim_setup.cpp`; `configs/g1_viewer_graded.json`);
  `.scn` persists the spec + regenerates. All gated by the `parity_probe` oracle (above).
  **windloads asserts the building bbox ⊆ the uniform fine core** (its `h_fine²`-face math is exact
  only there) — enlarge the `fine_core` box if it throws. **Remaining**: GUI dock controls for the
  fine core (task 3.4 — config path works today); the voxel-staircase overlay in `slice_viewer.cpp`
  still draws at uniform spacing (display polish); and the physical **rounded-vs-sharp validation
  runs** (below) are user-driven GPU jobs. Full proposal/design/tasks in
  `openspec/changes/graded-structured-grid/`.
- **Corner-resolution protocol (design D7, `grid_metrics.h` helpers).** The rounded-vs-sharp signal
  rides on the discretization error, so `h_fine` is not free: **size it ≈ r/10** (≥10 cells across
  the corner radius, e.g. 30 mm for a 300 mm radius — `recommended_h_fine()`); a run at **`h_fine ≥
  r/4` is below the resolution floor** (rounded voxelizes ≈ sharp — `below_resolution_floor()`
  warns) and is NOT comparison-grade. Before trusting a delta, run the **grid-convergence gate**:
  rerun the rounded case at a finer `h_fine` (e.g. 30 mm → 20 mm) and confirm Cd + peak-suction Cp
  are stable within the time-averaging RMS band. A/B pairs (rounded vs sharp) MUST share the domain,
  `fine_core` spec, and inflow — differ only in the corner geometry. Judge convergence off the
  **time-averaged** loads (`LoadAverager`), not an instantaneous snapshot.
- **DEAD-END (tried + reverted): open far-field BC.** An opt-in Dirichlet-p=0 "open" mode on the
  top/side faces (to relieve blockage in a small domain) gave textbook Cp at ≤~5% blockage but
  **backflow-diverges at higher blockage** — the regime where it's needed. Not committed; the
  graded grid (big cheap domain) is the path instead. A robust version would need a
  convective/backflow-stabilised outflow, not a raw pressure outlet.
- **Blockage caveat (in the meantime).** With closed free-slip top/side walls, keep the domain
  large enough that frontal **blockage < ~5%** or Cp inflates badly (a small domain gave windward
  Cp ~10 vs the physical ~+1) — the flow nozzles around the building. This is exactly what the
  graded grid lets you afford cheaply.

## History

WindCFD is forked from a former marine/sediment simulator; all of that sediment, seabed and bed
-evolution simulation has been removed (`src/core/` now holds only `fluid` + `geometry`, and the
build no longer compiles any sediment sources). Companion docs — `RESEARCH.md` (physics &
numerics spec), `PLAN.md` (roadmap, milestones W1–W8), `HANDOVER.md` (current status) — were all
retargeted to WindCFD alongside this file.
