## ADDED Requirements

### Requirement: Axis-separable graded grid representation
The `MacGrid` SHALL represent cell spacing as three per-axis metric arrays — `dx[i]`, `dy[j]`,
`dz[k]` (cell widths) plus cumulative face coordinates `xf[]`, `yf[]`, `zf[]` from which cell
centres are derived — such that spacing along each axis depends only on that axis's index. The
grid SHALL remain logically Cartesian: the `pidx`/`uidx`/`vidx`/`widx` layout and index topology
are unchanged. The representation MUST expose both the face-to-centre and centre-to-centre
distance families (they differ on a graded grid and are not equal to `dx[i]`).

#### Scenario: Metrics carry cell geometry
- **WHEN** the grid is queried for cell `i`'s width, its centre coordinate, and the centre-to-
  centre distance to cell `i+1`
- **THEN** each value is returned from the metric arrays with no assumption of a single scalar `h`

#### Scenario: Index topology preserved
- **WHEN** a kernel addresses cells/faces via `pidx`/`uidx`/`vidx`/`widx`
- **THEN** the addressing is identical to the uniform-grid layout (only spacing values differ)

### Requirement: Uniform grid is the exact special case
A grid whose per-axis spacings are all equal to a constant `h` MUST reproduce the previous
uniform-grid behaviour. Constructing a `MacGrid` with uniform metrics SHALL yield operator
results identical (within the existing 1e-5 parity tolerance) to the pre-change uniform solver.

#### Scenario: Collapse-to-uniform regression
- **WHEN** the solver runs with every `dx`/`dy`/`dz` set equal to `h`
- **THEN** velocity, pressure, and derived fields match the uniform-grid reference at rel.
  max-norm ≤ 1e-5

### Requirement: Invertible world↔index mapping
The grid SHALL provide a routine mapping a physical world coordinate to a fractional cell index
and back, using the cumulative face coordinates. The mapping MUST be correct on non-uniform
spacing and is the single shared routine used by advection departure-point lookup and the slice
sampler.

#### Scenario: Round-trip mapping
- **WHEN** a world point inside the domain is mapped to a fractional index and back to world
- **THEN** the recovered coordinate equals the original within floating-point tolerance

#### Scenario: Fine/coarse point lands in the right cell
- **WHEN** a point in the fine core and a point in the coarse far field are each mapped to indices
- **THEN** each index identifies the cell whose face bounds enclose the point

### Requirement: Grid generation from a fine-core specification
The system SHALL generate the metric arrays from a fine-core specification: a fine-core box
(`x0..x1, y0..y1, z0..z1`), a target `h_fine`, a growth ratio, and the domain bounds. Inside the
core, spacing SHALL be uniform at `h_fine`; outside, spacing SHALL grow geometrically toward the
boundaries at no more than the specified growth ratio. The generator SHALL run at Apply/Build
time (rebuild at t=0), consistent with the existing domain-resize flow.

#### Scenario: Fine core is uniform at h_fine
- **WHEN** a grid is generated with a fine-core box and target `h_fine`
- **THEN** every cell whose centre lies inside the box has width `h_fine` on each axis

#### Scenario: Growth ratio is bounded
- **WHEN** the grid transitions from the fine core to the coarse far field
- **THEN** the ratio of adjacent cell widths never exceeds the specified growth ratio

#### Scenario: Apply regenerates metrics
- **WHEN** the user changes the fine-core spec and presses Apply
- **THEN** the metric arrays are regenerated and the sim rebuilds at t=0 on the new grid

### Requirement: Grid consumers use the world↔index mapping
Every display and geometry consumer that locates a physical point in the grid — the CUDA-GL slice
sampler, the CPU flow tracers/particles, and near-building membership tests — SHALL locate points
via the grid's world↔index mapping rather than `floor(x/h)`, and SHALL derive domain extent from
the cumulative face coordinates (`xf[nx]`, `yf[ny]`, `zf[nz]`) rather than `nx·h`.

#### Scenario: Tracer locates in the correct cell on a graded grid
- **WHEN** a flow tracer at a world position in the coarse far field is located
- **THEN** it maps to the enclosing coarse cell (not a `floor(x/h)` cell that ignores grading)

#### Scenario: Domain extent from face coordinates
- **WHEN** any consumer needs the domain length along an axis
- **THEN** it reads the terminal cumulative face coordinate, not cell-count × a scalar `h`

### Requirement: Scene persistence of graded grids
A saved scene (`.scn`) SHALL persist enough to reconstruct the graded grid on load. Legacy
uniform scenes (no grading recipe) MUST load correctly as the all-equal special case.

#### Scenario: Round-trip a graded scene
- **WHEN** a graded-grid scene is saved and reloaded
- **THEN** the reconstructed grid produces the same cell widths, centres, and world↔index mapping

#### Scenario: Legacy uniform scene still loads
- **WHEN** a scene saved before grading (uniform `h` only) is loaded
- **THEN** the grid is reconstructed as a uniform grid with all spacings equal to that `h`
