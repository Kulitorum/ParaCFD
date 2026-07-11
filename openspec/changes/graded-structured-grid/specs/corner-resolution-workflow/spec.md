## ADDED Requirements

### Requirement: Fine-core resolution sizing rule
The fine core SHALL be sized so the target corner radius is resolved by at least ~10 cells across
the radius. For a target radius `r`, the recommended starting `h_fine` SHALL be `≈ r/10` (e.g.
`h_fine = 30 mm` for `r = 300 mm`). The building voxelizer SHALL use the **local** cell size when
determining the wall band, so that inside the fine core the sub-cell wall-thickness floor
(`0.5·h·√2`) does not fire and the printed wall (~80 mm) resolves with multiple cells.

#### Scenario: Recommended h_fine for a 300 mm radius
- **WHEN** the target corner radius is 300 mm
- **THEN** the sizing rule yields a starting `h_fine` of ~30 mm (~10 cells across the radius)

#### Scenario: Local cell size drives the wall band in the fine core
- **WHEN** the building is voxelized with the fine core covering the walls
- **THEN** the wall band uses the local (fine) cell size and the sub-cell floor does not alter the
  band for cells inside the core

### Requirement: Resolution floor for a meaningful corner comparison
A rounded-vs-sharp corner comparison MUST NOT be treated as valid when the fine-core cell size is
coarser than one quarter of the corner radius (`h_fine ≥ r/4`), because the voxelized rounded
corner is then indistinguishable from the sharp corner and any measured difference is a
discretization artifact. The workflow SHALL surface this floor to the user.

#### Scenario: Reject an under-resolved comparison
- **WHEN** a rounded-vs-sharp comparison is attempted with `h_fine ≥ r/4`
- **THEN** the workflow flags the run as below the resolution floor and not comparison-grade

#### Scenario: Above the floor is permitted
- **WHEN** `h_fine` is comfortably below `r/4` (e.g. `r/10`)
- **THEN** the run is eligible to proceed to the grid-convergence gate

### Requirement: Grid-convergence gate before trusting a comparison
Before a sharp-vs-rounded load comparison is reported as trustworthy, the rounded case SHALL be
run at two fine-core resolutions (e.g. `h_fine = 30 mm` and a finer `20 mm`) and the integrated
loads (Cd and peak-suction Cp at minimum) compared. The comparison SHALL be considered
grid-converged only when those loads are stable between the two resolutions within the
time-averaging RMS band.

#### Scenario: Converged — comparison earned
- **WHEN** the rounded case's Cd and peak-suction Cp shift between `h_fine` levels by less than the
  averaging RMS band
- **THEN** the grid is declared converged and the sharp-vs-rounded comparison may be reported

#### Scenario: Not converged — refine further
- **WHEN** the loads shift by more than the averaging RMS band between the two resolutions
- **THEN** the result is declared not grid-converged and a finer `h_fine` is required before
  comparison

### Requirement: Load integration is valid on a graded grid
Because the wind-load integration (`windloads`) assumes square exposed faces of area `h²`, the
building's exposed surface MUST lie entirely inside the **uniform** fine core, where every face is
`h_fine²`. The system SHALL feed `h_fine` (the fine-core spacing) to the load integration and
SHALL assert the building's bounding box ⊆ the uniform core; a surface cell in the graded
transition (where per-axis face areas differ and `A_frontal = count·area` breaks) MUST be treated
as an error, not silently integrated.

#### Scenario: Building inside the uniform core integrates correctly
- **WHEN** the building bbox lies within the uniform fine core and loads are integrated with
  `h_fine`
- **THEN** `A_frontal`, `A_plan`, `L_ref`, forces, and moments are computed on `h_fine²` faces

#### Scenario: Surface cell in the transition is rejected
- **WHEN** any exposed building surface cell falls in the graded transition region
- **THEN** the load integration flags an error rather than reporting wrong coefficients

### Requirement: Comparison runs share grid and conditions
A rounded and a sharp run compared as an A/B pair SHALL use the same domain, the same fine-core
specification, and the same inflow/physics settings, differing only in the corner geometry, so
the measured load delta is attributable to the corner shape rather than grid or boundary
differences.

#### Scenario: Matched A/B pair
- **WHEN** a rounded run and a sharp run are compared
- **THEN** their domain bounds, fine-core spec, and inflow/physics settings are identical and only
  the corner radius differs
