## Why

WindCFD's motivating question — how **rounded vs sharp building corners** change the wind loads —
is currently **undetectable at the grid resolution the tool runs**. The uniform MAC grid welds
feature resolution to domain size (`feature_h ≡ domain / N`), so a domain large enough for free
flow forces cells too coarse to see a corner. Concretely: `configs/building.json` runs at
`voxel_h = 0.25 m`, and the target corner radius is `~0.30 m` → only **~1.2 cells across the
radius**. A "rounded" corner voxelizes to essentially the same staircase as a sharp one; the
effect we exist to measure is invisible. A uniform grid fine enough to fix this (0.03 m over the
whole 40×30×15 m tunnel ≈ 6.7×10¹¹ cells) is impossible.

We need to **decouple feature resolution from domain size**. Because the refinement region is
**known and static** (the building is gizmo-placed before the sim runs) and the **wake is an
artifact we do not need to resolve**, an axis-separable **graded (stretched) structured grid** —
fine core around the building, coarse far field — is the right tool. It keeps the GPU-perfect
Cartesian data structure, the zero-copy slice sampler, and the MGPCG topology, and avoids the
octree/AMR rewrite (which solves a *dynamic/unknown-refinement* problem we do not have).

## What Changes

- **BREAKING (internal)**: `MacGrid`'s single scalar `double h` becomes three per-axis metric
  arrays `dx[i]`, `dy[j]`, `dz[k]` plus cumulative face/centre coordinates. **Uniform is the
  exact special case** where every spacing equals `h` — the migration preserves this so existing
  results reproduce bit-for-bit.
- **All fluid operators become metric-aware**: divergence/gradient use local cell widths;
  MacCormack advection inverts a non-uniform world→index map at departure points; Smagorinsky
  filter width becomes per-cell `(dx·dy·dz)^(1/3)`; the pressure Poisson operator becomes a
  **variable-coefficient** finite-volume Laplacian (`count/h²` → `Σ_f w_f`, `w_f = A_f/(d_f·V)`),
  still SPD so CG/MGPCG stay valid.
- **Grid generation from a fine-core spec**: config/GUI declare a fine-core box + target
  `h_fine` + growth ratio (≤ 1.15) + domain bounds; a generator emits the metric arrays. The
  fine box may auto-track the placed building's bounding box.
- **Grid consumers updated**: the building voxelizer uses **local** cell size (retiring the
  `0.5·h·√2` wall-thickness floor for cells inside the fine core); the CUDA-GL slice sampler
  gains the world→index map; `.scn` scene IO persists (or regenerates) the grid metrics.
- **Physics-validity protocol added**: an `h_fine`-sizing rule tied to corner radius (≥ ~10
  cells/radius; hard floor at `h ≥ r/4`) and a **grid-convergence gate** — rerun the rounded case
  at a finer `h_fine` and confirm loads are grid-independent **before** a sharp-vs-rounded
  comparison is trusted.
- **Tests**: every kernel keeps its 1e-5 GPU-vs-CPU parity test in a metric-aware form; add a
  "grading collapses to uniform" regression and a manufactured-solution (MMS) test on the graded
  Laplacian.

## Capabilities

### New Capabilities
- `graded-grid`: the axis-separable graded structured grid — per-axis metric representation,
  invertible world↔index mapping, generation from a fine-core spec, `.scn` persistence, and the
  uniform-as-special-case invariant.
- `graded-fluid-operators`: all MAC solver kernels (divergence/gradient, MacCormack advection,
  diffusion + Smagorinsky, variable-coefficient pressure projection/MGPCG, BCs) operate correctly
  on variable spacing, with parity + MMS + collapse-to-uniform test coverage.
- `corner-resolution-workflow`: the `h_fine` sizing rule, the resolution floor, and the
  grid-convergence gate that make a rounded-vs-sharp corner comparison physically trustworthy.

### Modified Capabilities
<!-- None: openspec/specs/ is currently empty; there are no existing requirement specs to modify. -->

## Impact

- **⚠ Prerequisite (no existing test oracle)**: CLAUDE.md advertises a "1e-5 GPU-vs-CPU parity
  test per kernel" invariant, but **no test target exists** (CMake builds only `building_probe` +
  `windcfd-gui`; the CPU twins are defined but never run). A parity harness **or** a golden-master
  snapshot must be built as **Phase 0** before the operator rewrite, or "zero behaviour change" is
  unfalsifiable. See `design.md` D4.
- **Core solver** (`src/core/fluid/`): `mac_grid.h` (representation), `advect.cu`,
  `turbulence.cu`, `project.cu`/`channel_pressure.cu`/`mgpcg.cu`, `channel_ops.cu`,
  `channel_core.cu` (**incl. `adaptive_dt` → `h_min`**), `mac_ops.h`, `bedshear.*` (wall-normal
  distance; deferred — not wired into building configs) — ~300 `h` call sites.
- **Wind loads** (`src/core/windloads.cpp`) — **the Cd/Cl/Cp deliverable**: hardcodes `area=h·h`,
  `A_frontal/A_plan=count·area`, `L_ref`, and moment arms in scalar `h`. Feed `h_fine` + **assert
  building bbox ⊆ uniform core** (loads are wrong the moment a surface cell lands in the graded
  transition).
- **Geometry** (`src/core/geometry/`): `voxelize.*`, `building.cpp` (local-cell wall floor — now
  committed, `d648cb6`), `model_placement.h` (world↔index consumers).
- **GUI** (`src/gui/`): `slice_field.cu`/`slice_gl.cu` (sampler world→index), **`flow_particles.cpp`
  + `flow_tracers.cpp`** (CPU tracer `floor(x/h)` locate + `nx·h` domain extent → world→index map),
  `sim_setup.*` (grid generation + Apply path), `scene_io.*` (persist the fine-core recipe, matching
  `7a8f346`'s STEP-as-source-of-truth), `main_window.*` (fine-core controls).
- **Config** (`src/core/config.*`, `configs/building.json`): fine-core box + `h_fine` + growth
  keys.
- **Numerics**: formal accuracy drops from 2nd → ~1st order on stretched cells (accepted at
  growth ≤ 1.15); timestep is set by the smallest cell (more steps per flow-through, but only for
  the refined region).
- **No new external dependencies.** Layering invariants (Qt/OCC/GL islands) unchanged.
