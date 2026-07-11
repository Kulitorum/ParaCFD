# HANDOVER.md — WindCFD current status (2026-07-11)

This is a fresh-session pick-up doc. It reflects the code as it stands now: the building
pipeline and wind-load extraction that earlier handovers listed as "planned" are
**implemented and working**. Skim `CLAUDE.md` for the module map; this file is the
current-state snapshot and the next-steps list.

## 1. What WindCFD is + the goal
WindCFD is a GPU (CUDA) 3D incompressible-flow **LES** CFD tool (`C:/CODE/WindCFD`,
C++/CUDA + Qt6, namespace `windcfd::core` / `windcfd::gui`). The fluid is **air/wind**
(ρ ≈ 1.225 kg/m³, ν ≈ 1.5e-5 m²/s). A CAD building is voxelized into a solid obstacle,
the LES solver develops the flow around it, and the Qt/OpenGL viewer shows live
velocity/pressure slices, flow arrows and the model.

The motivating question is a 3D-printing design choice: **rounded vs sharp building
corners**, and how that changes the flow field and the **wind loads on the facade and
roof** — a bluff body in an atmospheric boundary layer (ABL).

## 2. What is implemented now (all DONE, verified against the code)

- **Air fluid defaults.** `Config` (`src/core/config.h`) ships **air**: `rho = 1.225`,
  `nu = 1.5e-5`, `U = 1.0`. (The old inherited seawater values are gone.) The `d50` field
  remains only as a log-law wall-roughness length.
- **No default obstacle.** The old procedural solid-cylinder obstacle was **deleted** — a
  run now starts as an **empty channel**; obstacles come only from a loaded building
  centerline or STEP model.
- **Building pipeline** (`src/core/geometry/building.{h,cpp}`, OCC-free, in `libwindcfd`):
  a 3D-printing **CENTERLINE** STEP is loaded via `windcfd_geometry`'s `step_import` into a
  `TriMesh` of the vertical wall surface ribbons, then:
  - horizontally sectioned (`mesh_horizontal_section`) into **2D footprint loops**;
  - voxelized **directly** (`voxelize_building`) — no watertight solid: a cell is solid when
    it lies within ± half `wall_thickness` of the (optionally corner-rounded) centerline over
    the wall height, which is leak-proof by construction;
  - plus a **solid flat roof** = convex hull of the footprint dilated by the roof overhang,
    over `roof_thickness`.
  - Corners are configurable: a **ROUNDED** outer radius is obtained by pre-filleting the
    centerline (2D fillet); `corner_radius = 0` gives **TRUE SHARP** corners by adding outward
    **miter-wedge** triangles at each convex corner (with a miter limit to bevel spikes).
  - `model_placement.h::placed_mesh()` applies the gizmo transform (translate/rotate/scale)
    to the mesh **before** section+voxelize, so the building is built where it is placed.
  - Dev tool `tools/building_probe.cpp` (CMake target `building_probe`) exercises
    centerline-STEP → building-solid headlessly.
- **Wind loads** (`src/core/windloads.{h,cpp}`, host-only, OCC-free, in `libwindcfd`):
  `compute_wind_loads()` integrates the CFD **pressure** over the building's exposed voxel
  faces (solid↔fluid) → forces `Fx/Fy/Fz` and moments about the centroid, and coefficients
  **Cd** (drag +x), **Cs** (side +y), **Cl** (uplift +z) via projected frontal/plan
  reference areas. Surface pressure at a face is the adjacent fluid cell's pressure (Pa);
  `p_ref` is the **upstream inlet-slab mean** (freestream). Per-cell surface **Cp** is
  returned for visualisation, `Cp = (p − p_ref)/(½ρU²)`. Validated headless on a test house:
  windward stagnation **Cp ≈ 1.0**, **Cd ≈ 1.2**, **Cs ≈ 0**.
- **GUI** (`windcfd-gui`):
  - File ▸ **"Open centerline STEP…"** and the `--load-centerline <path>` CLI flag.
  - A **"Building (centerline → solid)"** dock group: wall thickness, wall height, corner
    radius (0 = sharp), roof overhang, roof thickness, and a **Build** button.
  - The centerline model is **gizmo-placeable** (translate/rotate/scale) with **30-level
    UNDO/REDO** (Ctrl+Z / Ctrl+Shift+Z or Ctrl+Y, plus Undo/Redo buttons), a **live
    placement readout**, and a **Reset placement** button.
  - **Build** voxelizes the **placed** model into the **CURRENT** domain — **no auto-resize**;
    the user sets domain size + voxel size via the Domain/resolution controls. **Apply**
    (domain/resolution) rebuilds the grid at t=0 and re-voxelizes the placed house.
  - A **"Start Simulation"** button **gates stepping**: the worker starts **paused** so the
    site can be prepared (load, place, Build, size the domain) with the flow frozen; display
    and Build/voxelize still work. Headless runs auto-start when `--autoclose-ms` is given.
  - Live **Cd / Cl / Cs** and **Cp-range** readouts.
  - The building voxel overlay is **coloured by Cp** (diverging ramp) with a legend and a
    "Colour building by Cp" toggle.
  - **Input speed U is capped at 40 m/s**.
- **Configs** (`configs/`): **`building.json`** — a ~40×30×15 m AIR wind-tunnel default at
  h=0.25 m (~1.15M cells), **auto-selected** when `--load-centerline` is passed without
  `--config`. The **cylinder validation configs were removed**; `g1_viewer*` are now
  **empty-channel** viewer demos. `v1_cavity_re100.json` and `v2_channel_loglaw.json` remain
  as core/wall-model validation configs.

## 3. The interactive workflow
1. Launch with `configs/building.json` (or just pass `--load-centerline`, which selects it).
2. **File ▸ Open centerline STEP…** — loads the 3D-print centerline; it is auto-Built once.
3. **Gizmo-place / rotate / scale** the model to site it in the tunnel (Undo/Redo as needed).
4. **Set the domain size + voxel size** (Domain/resolution dock) and **Apply** if you changed
   them — Build does **not** resize for you.
5. Set the wall/roof/corner params in the **Building** group and press **Build** to voxelize
   the placed model into the current domain (swap corner radius 0 ↔ >0 to compare sharp vs
   rounded; Build re-voxelizes in place).
6. Press **Start Simulation** to begin advancing the flow.
7. Read **Cd / Cl / Cs** and the **Cp range**; the building overlay tints by surface Cp.

## 4. Build & run
Toolchain: **VS 2022 Professional** (MSVC 14.44) — **NOT VS 2026** (CUDA 13.1 rejects the
2026 toolset) — CUDA 13.1, Qt 6.11.1 (`C:/Qt/6.11.1/msvc2022_64`, glob-detected),
OpenCascade 8.0 (`C:/OpenCASCADE-8.0/build2`). Put the VS 2022 `vcvars64` environment **and
the VS-bundled Ninja** on PATH first (the VS 2022 x64 dev prompt does both), then:
```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build --target windcfd-gui
```
`-DCMAKE_CUDA_ARCHITECTURES=89` is the fast single-arch (Ada RTX 4090) dev build; the CMake
default `75;86;89` is the multi-arch shipping build (~3× device-code compile time). A
POST_BUILD step deploys the Qt (`windeployqt`) and OpenCascade DLL closures next to the exe,
so it is double-clickable. If `CL.exe` crashes mid-AUTOMOC (`-1073741819`), just re-run — it
is transient.

Targets: `libwindcfd` (core, incl. `windloads.cpp` + `building.cpp`), `windcfd_geometry`
(the only OCC target — STEP import), `windcfd_gui_cuda` (the only CUDA-GL target),
`windcfd-gui` (the Qt app), and `building_probe` (dev tool).

Run (windowed):
```
build/windcfd-gui.exe --config configs/building.json
build/windcfd-gui.exe --load-centerline <path-to-centerline.step>   # auto-picks configs/building.json
```
Run (headless smoke — solver steps, exits clean):
```
build/windcfd-gui.exe --config configs/building.json --offscreen --autoclose-ms 5000
```
`--offscreen` has **no GL context**, so it **renders 0 frames** but **does step the CUDA
solver** (and auto-starts the run), exiting 0 after N ms. Use it for CI / solver debugging,
not for visual checks. Other CLI flags (see `src/gui/main.cpp`): `--load-step`,
`--model-place dx,dy,dz[,rz[,scale]]`, `--set-u`, `--tidal`, `--clip`, `--tracers`,
`--save-scene`/`--load-scene`, `--record`.

## 5. Immediate next steps (priority order)
1. **Convergence + time-averaging (highest priority).** The live Cd/Cl/Cs and Cp are
   **instantaneous** values read off an **unsettled** flow, so the rounded-vs-sharp
   comparison is **not yet trustworthy**. Add flow-steadiness detection, then **mean / RMS /
   peak** load accumulation, and a batch/CSV **corner-radius sweep** (headless) so a
   sharp-vs-rounded comparison can be produced repeatably.
2. **ABL wind inflow profile.** Inflow is currently **uniform**; wire in a log-law / power-law
   ABL velocity profile plus inlet turbulence, reusing the existing log-law wall model and the
   Jarrin SEM synthetic-eddy inlet infrastructure.
3. **Building Reynolds number.** `configs/building.json` uses a **moderate `nu = 1.0e-3`**
   (numerical stability aid), not true air ν — revisit once the flow settles and the ABL inlet
   is in, so the building runs at a physical Re.
