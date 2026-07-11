# HANDOVER.md — WindCFD current status (2026-07-11)

## 1. What WindCFD is + the goal
WindCFD is a GPU (CUDA) 3D incompressible-flow LES CFD tool (`C:/CODE/WindCFD`,
namespace `windcfd::core`). The fluid is **air/wind**. The goal is to quantify wind
flow around **3D-printed concrete buildings** and the wind **loads on the facade and
roof**, comparing **rounded vs sharp** corners — bluff-body aerodynamics in an
atmospheric boundary layer (ABL).

⚠ Known gap: the code `Config` still ships **inherited seawater defaults**
(`rho = 1027`, `nu = 1.36e-6` in `src/core/config.h`). These must be switched to air
before any wind result is meaningful — see next steps.

## 2. What was just done
WindCFD was **forked from a scour-protection (sediment/seabed) simulator** and stripped
down to a fluid-only wind tool. Fresh git history, four commits (oldest → newest):

1. `810c122` — baseline: the scour codebase copied in (pre-refactor).
2. `e962030` — remove scour: delete all sediment/morphodynamics/seabed code, the CLI,
   all tests/gates, trim `Config`, collapse `CMakeLists.txt`.
3. `2657600` — de-scour the GUI to a **fluid-only slice viewer**; wire CUDA includes.
4. `1ec2366` — rename ScourProtection → WindCFD (namespace, targets, macros, installer,
   files).

**Removed** (so nothing below misleads): all sediment/morphodynamics/seabed/scour
physics, the command-line tool, all unit tests and ctest gates, and the Jolt
drop/settle feature; docs were trimmed. Everything **builds and runs**.

## 3. Current architecture & targets
Fluid core and infrastructure that remain and build:
- **Fluid solver** — MAC grid, MacCormack advection, Smagorinsky LES, MGPCG pressure
  projection (`src/core/fluid/`).
- **Open channel + voxelized obstacles**, sub-grid **porous momentum sink**, **log-law
  wall function** (ground/ABL wall model), **periodic channel + PI body-force control**,
  **SEM turbulent inlet** + precursor.
- **Geometry** — OpenCascade **STEP import** → mesh, **ray-parity voxelizer** to a solid
  obstacle mask (`src/core/geometry/`).
- **IO/config** — VTI read/write, JSON config loader (`src/core/`).
- **windcfd-gui** — Qt6 + GL 4.3 slice viewer (`src/gui/`) with the flow slice, arrows,
  camera, voxel overlay, STEP loading, and a placement gizmo.

CMake targets: `libwindcfd` (Qt-free static core + CUDA), `windcfd_geometry` (OpenCascade
STEP import, isolated), `windcfd_gui_cuda` (CUDA slice kernel + GL interop), `windcfd-gui`
(the Qt viewer executable).

## 4. How to build & run
Toolchain: **VS 2022 Professional** (MSVC 14.44) — **not VS 2026** (CUDA 13.1 rejects it)
— CUDA 13.1 (`CMAKE_CUDA_ARCHITECTURES` default `75;86;89`), Qt 6.11.1, OpenCascade 8.0.

Primary build (VS2022 generator picks the 14.44 host compiler automatically):
```
cmake -B build -S . -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```
Ninja fallback (must run inside the VS2022 Professional x64 dev environment so nvcc finds
cl 14.44):
```
cmd /c "\"C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat\" && cmake -B build -S . -G Ninja && cmake --build build"
```
Run (windowed):
```
build/Release/windcfd-gui.exe --config configs/g1_viewer.json
```
Run (headless smoke — solver steps, exits clean):
```
build/Release/windcfd-gui.exe --config configs/g1_viewer.json --offscreen --autoclose-ms 5000
```
Note: `--offscreen` has no GL context, so it renders 0 frames, but the CUDA solver does
step and the process exits 0.

## 5. Immediate next steps (the wind-load work)
1. **Switch physics defaults to air** (`rho`/`nu` in `Config`; audit any other inherited
   seawater constants).
2. **ABL wind inflow** — atmospheric-boundary-layer inlet profile (reuse the log-law wall
   model + SEM inlet).
3. **Building geometry** — import/generate **rounded vs sharp** cornered buildings and
   voxelize them.
4. **Surface pressure + loads** — compute facade/roof surface **Cp** and integrate the
   **wind loads** on facade + roof.
5. **Comparison study** — the rounded-vs-sharp load comparison.
