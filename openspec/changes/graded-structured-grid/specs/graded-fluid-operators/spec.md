## ADDED Requirements

### Requirement: Metric-aware differential operators
The divergence, gradient, and diffusion operators SHALL use local per-axis spacing from the grid
metrics rather than a single scalar `h`. Divergence at a cell SHALL divide face velocities by the
local cell width; the pressure gradient updating a face velocity SHALL divide the pressure
difference by the centre-to-centre distance (not the cell width). Each operator MUST retain a CPU
reference implementation and pass GPU-vs-CPU parity at rel. max-norm ≤ 1e-5.

#### Scenario: Divergence uses local cell width
- **WHEN** divergence is evaluated in a graded region where adjacent cells differ in width
- **THEN** the result matches the finite-volume divergence computed from the local metrics

#### Scenario: Operator parity holds on a graded grid
- **WHEN** any metric-aware operator runs on GPU and on its CPU reference over a graded grid
- **THEN** the two results agree at rel. max-norm ≤ 1e-5

### Requirement: MacCormack advection on non-uniform spacing
Semi-Lagrangian MacCormack advection SHALL locate departure points using the grid's invertible
world→index mapping and interpolate with weights derived from local cell widths. Accuracy MUST be
preserved for growth ratios within the supported bound.

#### Scenario: Departure-point lookup respects grading
- **WHEN** a departure point falls in a coarse cell adjacent to fine cells
- **THEN** the interpolation stencil and weights are those of the enclosing coarse cell, not a
  uniform-`h` neighbourhood

#### Scenario: Advection parity on a graded grid
- **WHEN** advection runs on GPU and on its CPU reference over a graded grid
- **THEN** the two results agree at rel. max-norm ≤ 1e-5

### Requirement: Per-cell Smagorinsky filter width
The Smagorinsky LES eddy-viscosity model SHALL compute its filter width per cell as
`(dx·dy·dz)^(1/3)` from the local metrics, replacing the constant `h`.

#### Scenario: Filter width varies with local cell size
- **WHEN** eddy viscosity is evaluated in the fine core and in the coarse far field
- **THEN** the filter width used is the local `(dx·dy·dz)^(1/3)` of each cell

### Requirement: Variable-coefficient pressure projection
The pressure Poisson solve SHALL use a variable-coefficient finite-volume Laplacian:
`Ap[t] = Σ_f w_f · (p[t] − p[nb_f])` with `w_f = A_f / (d_f · V_cell)` and `diag = Σ_f w_f`, which
reduces exactly to `count/h²` on a uniform grid. The operator MUST remain symmetric positive
definite so CG/MGPCG stay valid. Multigrid restriction and prolongation SHALL be volume-weighted
so the coarse operator remains a good approximation; coarsening remains by-2 in index space.

#### Scenario: Projection produces a divergence-free field on a graded grid
- **WHEN** projection is applied to a divergent velocity field on a graded grid
- **THEN** the resulting field's finite-volume divergence is reduced below `proj_tol`

#### Scenario: Poisson operator collapses to the uniform stencil
- **WHEN** the variable-coefficient operator is evaluated with all spacings equal to `h`
- **THEN** it equals `(count·p − nbsum)/h²` at rel. max-norm ≤ 1e-5

#### Scenario: Multigrid convergence does not regress catastrophically
- **WHEN** the graded MGPCG solves to `proj_tol` with volume-weighted transfers
- **THEN** the iteration count stays within a bounded factor of the uniform-grid baseline

### Requirement: Boundary and wall models use metric distances
Inlet/outlet/lid/mask boundary conditions and the log-law wall model (`bedshear`) SHALL use the
actual metric distances (e.g. the wall-normal spacing of the first off-wall cell) rather than a
constant `h`.

#### Scenario: Wall model reads the true first-cell distance
- **WHEN** the log-law wall function samples the first off-wall cell on a graded grid
- **THEN** it uses that cell's actual wall-normal distance from the metrics

### Requirement: Stable timestep on a graded grid
The adaptive timestep SHALL be limited by the **minimum** cell dimension across the grid, not a
single scalar `h`. Both the advective (`cfl·h/|u|`) and diffusive (`h²/(6ν)`) limits MUST use
`h_min` precomputed from the metric arrays, so the fine core does not go unstable.

#### Scenario: Timestep respects the finest cell
- **WHEN** `adaptive_dt` is evaluated on a graded grid whose smallest cell is in the fine core
- **THEN** the returned dt satisfies the CFL and diffusive limits for that smallest cell

### Requirement: Behaviour-preserving migration oracle
Before the scalar-`h` operators are converted to metric-aware form, the change SHALL provide a
falsifiable oracle: either a parity harness driving each CPU-twin reference against its GPU kernel
at rel. max-norm ≤ 1e-5, or a golden-master snapshot of `{u,v,w,p}` from a canonical run diffed
bit-for-bit. The uniform-initialised solver MUST pass this oracle after the operator rewrite with
no behaviour change.

#### Scenario: Oracle exists and passes on the pre-change solver
- **WHEN** the oracle runs against the current uniform solver
- **THEN** it passes, establishing the baseline

#### Scenario: Oracle still passes after the metric-aware rewrite (uniform metrics)
- **WHEN** the oracle runs after operators are converted, with all spacings equal to `h`
- **THEN** it still passes, proving the rewrite changed no behaviour

### Requirement: Manufactured-solution accuracy check
The change SHALL include a method-of-manufactured-solutions (MMS) test for the graded Laplacian
that quantifies the observed order of accuracy under grading and confirms it degrades gracefully
(≈2nd order at gentle growth, no worse than 1st order formally).

#### Scenario: MMS confirms expected order under grading
- **WHEN** the MMS test runs on a sequence of graded grids with growth ratio ≤ the supported bound
- **THEN** the measured error and order of accuracy are within the documented expectation
