## Context

ParaCFD is a GPU LES solver for wind around 3D-printed concrete buildings; its purpose is to
quantify how **rounded vs sharp corners** change facade/roof loads. The solver runs a Stam
stable-fluids scheme on a **uniform MAC grid** (`MacGrid`, `src/core/fluid/mac_grid.h:45`) with a
single scalar `double h` threaded through ~300 call sites across 38 files. The pressure Poisson
operator is a **constant-coefficient** 7-point stencil — `channel_ops.cu:353`:

```
Ap[t] = (count · p[t] − nbsum) / (g.h · g.h)      // count = # fluid neighbours
diag  =  count / (g.h · g.h)
```

**Current state / the wall we hit.** `configs/building.json`: 40×30×15 m domain, `voxel_h =
0.25 m` (160×120×60 ≈ 1.15M cells). Target corner radius `r ≈ 0.30 m` ⇒ `r/h ≈ 1.2` cells across
the radius. The rounded corner is indistinguishable from sharp after voxelization — the effect we
measure is below the grid's resolution. A uniform grid fine enough to fix this everywhere is
~10¹¹ cells; impossible. The 4090 holds ~90–100M double-precision cells at ~80% VRAM.

**Key constraints that pick the architecture:**
- The refinement region is **known and static** — the building is gizmo-placed *before* the sim.
- The **wake is an artifact** the user explicitly does not need resolved.
- Kernels are GPU-perfect today (coalesced, branchless); we want to keep that.
- Layering islands (Qt/OCC/GL) must stay intact.

**⚠ Corrected assumption (verified against the tree, code review 2026-07-11).** CLAUDE.md
advertises "every CUDA kernel has a CPU reference + GPU-vs-CPU parity test at 1e-5" as a live
invariant. **It is not live.** The CPU twins exist (`ch_poisson_apply_cpu`, `poisson_apply_cpu`,
`mac_advect_cpu`, …) but there is **no test target** — CMake defines only `building_probe` and
`paracfd-gui`; no `add_test`, no `enable_testing`, no gtest, nothing calls the twins. The parity
suite this design originally leaned on as "the oracle" **does not run**. Building an oracle is
therefore **Phase 0** (see D4 + Migration Plan), not an assumption.

## Goals / Non-Goals

**Goals:**
- Decouple feature resolution from domain size: a fine core (~10 cells across a 300 mm radius)
  around the building with a coarse, graded far field, within the 4090 VRAM budget.
- Preserve the Cartesian data structure, coalesced kernels, zero-copy CUDA-GL slice sampler, and
  MGPCG topology.
- Make the uniform grid a provable special case so migration is verifiable against current
  results (no "everything changed, nothing verifies" moment).
- Provide a defensible `h_fine` sizing rule and a grid-convergence gate so a sharp-vs-rounded
  comparison is physics, not discretization noise.

**Non-Goals:**
- No octree / cell-by-cell adaptive mesh refinement.
- No **dynamic** refinement (no solution-feature-following, no moving criteria) — refinement is
  static, set at Apply/Build time.
- No general **curvilinear / body-fitted** mesh (no full Jacobian, no cross-axis metrics).
- No wake resolution target; the far field is deliberately coarse.
- No change to the inflow model or viscosity (`nu=1e-3`, uniform inlet) — separate honesty gaps,
  out of scope here.
- Not (yet) multi-block **nested refinement boxes** — reserved as the escalation path if a single
  graded core cannot span a future scale ratio.

## Decisions

### D1 — Axis-separable (rectilinear tensor-product) grading, not octree, not curvilinear
The 3-D grid is the outer product of three independent 1-D meshes: `dx` depends only on `i`, `dy`
only on `j`, `dz` only on `k`. No tree, no neighbour search, no Jacobian, no cross terms.

- **Why over octree/AMR**: octree earns its cost only when refinement is unknown a priori, moves,
  or needs many levels — none apply here (static, gizmo-placed, ~2–3 effective levels). Octree is
  GPU-hostile (pointer chasing, warp divergence, poor coalescing), needs coarse-fine flux
  matching + hanging-node interpolation + a from-scratch octree multigrid, and would force a
  rewrite of the zero-copy slice sampler (it samples a regular grid). Cost ≫ benefit.
- **Why over curvilinear/body-fitted**: a full metric Jacobian per cell is a different, much
  larger project; separable grading keeps every finite-difference stencil diagonal-friendly.
- **Why over static nested boxes (block-structured)**: nested boxes are the more powerful
  middle option but add coarse↔fine interface interpolation and cross-level projection. Grading
  needs neither and covers a single, centrally-placed building. Nested boxes stay on the shelf as
  the escalation path (see Non-Goals / Open Questions).
- **Trade-off**: grading refines **slabs** (a fine band per axis), whose intersection is a fine
  **box** — you cannot refine "just the corners." For one building near domain centre this is
  fine; the interior solid cells inside the box are masked but still allocated.

### D2 — `MacGrid` representation: metric arrays + cumulative coordinates
Add per-axis `dx[]`, `dy[]`, `dz[]` (cell widths) and cumulative face-coordinate arrays
`xf[]`, `yf[]`, `zf[]` (plus derived cell centres `xc[]` = ½(xf[i]+xf[i+1])). Keep `pidx/uidx/
vidx/widx` and the memory layout **unchanged**. The distinction that carries the subtle bugs:
on a graded grid the **face-to-centre** distance (used for `div`) and the **centre-to-centre**
distance (used for the pressure gradient that updates a face velocity) differ and are *not*
`dx[i]` — they are averages of neighbouring widths. Encode both explicitly.

- **Sub-decision**: **tabulated** spacing (store the arrays), not a closed-form analytic stretch.
  Tabulated supports an arbitrary "uniform fine core + geometric transition + coarse far field"
  directly and inverts via a cumulative lookup; analytic stretching constrains the fine region to
  the function's shape. Cost is 6 small 1-D arrays in device memory — negligible.

### D3 — Variable-coefficient finite-volume pressure projection
Generalise the Poisson operator to per-face weights:

```
Ap[t] = Σ_f  w_f · (p[t] − p[nb_f]),   w_f = A_f / (d_f · V_cell),   diag = Σ_f w_f
```

`count/h²` is exactly the `A_f=h², d_f=h, V=h³` special case. The operator stays **SPD** (keep
the volume-scaled FV form), so CG and the MG smoothers remain valid. Only ~6 kernels change
(residual, Jacobi, R-B Gauss-Seidel, `apply_A`, divergence RHS, pressure-gradient correction) —
they all share the one stencil pattern.

**Global, not local.** Unlike advection/diffusion, the projection is a **global** solve: its MG
hierarchy coarsens across the fine core, the transition ramps, and the far field in a single
sweep. So the Poisson kernels must be **fully metric-aware everywhere** — the "uniform core ⇒
treat as scalar h" reframe (D0) does **not** apply to the solve, even though the coefficients
happen to be uniform inside the core.

- **Care-item A — volume-weighted transfers**: multigrid **restriction/prolongation** should be
  **volume-weighted** so the coarse operator stays a good approximation; naive injection can
  degrade the MG convergence *rate* (not correctness). Coarsening is still by-2 in index space
  (grid stays logically Cartesian).
- **Care-item B — per-level metrics by face-coordinate decimation (verified concern)**: each
  coarser level needs its own spacings, and the correct derivation is to **decimate the cumulative
  face coordinates** — `xf_coarse[i] = xf_fine[2i]`, hence `dx_coarse[i] = dx_fine[2i] +
  dx_fine[2i+1]` — stored per level alongside the coarsened solid mask. It is **not** an average
  of `dx`; averaging silently degrades the coarse operator.
- **Care-item C — interior-transition smoothing (verified concern)**: `ch_gs_band`
  (`channel_ops.cu`) sweeps only the *domain-edge* band (`i<band || i>=nx-band || …`). On a graded
  grid the sharpest coefficient jump is at the **interior transition plane**, which that band
  smoother never touches — the V-cycle leans on Jacobi + coarse correction there. Likely fine;
  instrument iteration counts vs the uniform baseline, and if the solve stalls, extend the R-B GS
  to a band around the transition planes, not just the domain edges.

### D4 — Migration: build an oracle FIRST, then uniform-as-special-case
The single most important risk-reduction move — but the original two-phase plan assumed an oracle
that **does not exist** (see Context ⚠). Corrected to three phases:

- **Phase 0 (build the oracle — do this before touching any operator)**: without it, "rewrite ~300
  sites with zero behaviour change, proven by the suite" is **unfalsifiable** — the exact
  "everything changed, nothing verifies" moment we are trying to avoid. Two acceptable forms:
  - **(a) Reconstruct a parity harness** — a gtest (or minimal `add_test`) target that drives the
    existing CPU twins (`ch_poisson_apply_cpu`, `mac_advect_cpu`, `smagorinsky_nut_cpu`, …) vs
    their GPU kernels at rel. max-norm 1e-5. Restores the invariant CLAUDE.md claims. Preferred.
  - **(b) Golden-master snapshot** — if (a) is too much up front: checksum/snapshot the
    `{u,v,w,p}` fields of a canonical short run (e.g. `configs/g1_viewer.json`, N steps) **before**
    any change, then diff bit-for-bit after Phase A. Cheap, coarse, but a real falsifiable gate.
  - Either way, **pull task 5.1 ("grading collapses to uniform") forward** — it is the natural
    permanent form of this oracle and should exist before Phase A lands.
- **Phase A (mechanical, zero behaviour change)**: add the metric arrays initialised uniform;
  rewrite every operator to read metrics instead of scalar `h`. The Phase-0 oracle MUST still pass
  unchanged. This lands the ~300-call-site churn with a falsifiable green gate as proof.
- **Phase B (the feature)**: generate non-uniform metrics from the fine-core spec; add the MMS
  Laplacian test and the permanent "grading collapses to uniform" regression; wire GUI/config/scene
  and the voxelizer/sampler/windloads/tracer consumers.

**Timestep (`adaptive_dt`) — a scalar-h site that bites in the fine core.** `channel_core.cu:169`
computes `dt_adv = cfl·g.h/maxu` and `dt_diff = g.h²/(6ν)`. On a graded grid dt MUST key off the
**minimum** cell dimension (precompute `h_min` across the metric arrays), or the step is too large
and the fine core goes unstable. (With `nu=1e-3` the advective limit binds, not the diffusive one
— but use `h_min` for both.)

### D5 — Grid generation from a fine-core spec (config + GUI)
Config/GUI declare: fine-core box (`x0..x1, y0..y1, z0..z1`), target `h_fine`, growth ratio
(≤ 1.15), domain bounds. Generator emits uniform spacing inside the core and geometric growth
outward to the boundaries. Fits the existing "set domain + voxel size → Apply → rebuild at t=0"
model; Apply regenerates the metrics. The fine box may auto-track the placed building bbox +
margin (see Open Questions). Transition span from 30 mm → 250 mm at ratio 1.15 is only ~15 cells
per side.

### D6 — Consumers use local cell size / world→index map
The building voxelizer tests membership using per-cell centres and floors the wall band by the
**local** cell size, so inside the fine core the `0.5·h·√2` wall floor (now **committed** to
master — `building.cpp`, commit `d648cb6`, uses scalar `g.h`) effectively never fires and the
80 mm printed wall resolves for real (~2.7 cells at 30 mm). The slice sampler and advection
departure-point lookup share one invertible world→index routine (cumulative-coordinate binary
search or O(1) table).

**Missing consumers the first Impact pass omitted (verified against the tree):**
- **`src/core/windloads.cpp` — the Cd/Cl/Cp deliverable, the whole point of the tool.** It
  hardcodes `area = g.h·g.h` (`:15`), `A_frontal/A_plan = count·area` (`:61–62`), `L_ref =
  (kmax−kmin+1)·g.h` (`:63`), and moment arms `(i+0.5)·h` (`:44–97`). If left on a scalar `g.h`
  (which won't exist), **every load on a graded grid is wrong.** Mitigation that exploits D0: the
  building lives entirely inside the **uniform** fine core, so all exposed faces are `h_fine²`.
  windloads can stay a scalar computation **iff** (a) it is fed `h_fine` (not a stale global `h`),
  and (b) a guard **asserts the building bbox ⊆ the uniform core** so no surface cell lands in the
  graded transition (where per-axis face areas differ and `A_frontal = count·area` breaks). Add
  the guard; it is a load-bearing correctness invariant, not a nicety.
- **`src/gui/flow_particles.cpp` + `flow_tracers.cpp` — display-only, but currently silent-wrong.**
  CPU tracer advection locates particles with `floor(x/g.h)` and samples with `h=g.h`
  (`flow_particles.cpp:36,50`; `flow_tracers.cpp:24,34`); domain extent is `Lx=g.nx·g.h`
  (`:92,127`). On a graded grid these mis-locate → wrong-looking arrows. Route them through the
  world→index map and take domain extent from the cumulative face coordinate `xf[nx]`, not `nx·h`.

### D0 — Reframe: metric complexity is confined to the far field + transition ramps
Because the fine core is **uniform**, the **explicit / local** operators and geometry consumers
that act purely within the building's core — the voxelizer wall floor, near-building membership
tests, windloads, tracer sampling — can treat `h_fine` as a scalar (with the D6 bbox⊆core guard).
True per-axis metric-awareness is only forced where a stencil **crosses the core↔far-field
transition** or sweeps the coarse far field. This **bounds the audit surface** — but note the
explicit exception in D3: the **pressure projection is global** and must be metric-aware
everywhere regardless. State this scoping in the tasks so reviewers know what must be audited vs
what provably reduces to the uniform case.

### D7 — `h_fine` sizing rule + grid-convergence gate
`h_fine` is not a free knob: the rounded-vs-sharp signal sits on top of the discretization error.

- **Rule**: target ≥ ~10 cells across the corner radius ⇒ start at **`h_fine = 30 mm`** for
  `r = 300 mm`. Hard floor: at `h ≥ r/4` (≥ 75 mm) the voxelized rounded corner ≈ sharp and any
  measured difference is an artifact — reject such runs for comparison.
- **Gate**: before trusting a sharp-vs-rounded delta, rerun the rounded case at a finer `h_fine`
  (e.g. 20 mm, ~15 cells/radius) and confirm Cd / peak-suction Cp are grid-independent (within the
  averaging RMS band). Only then is the comparison earned.

## Risks / Trade-offs

- **[Formal accuracy drops 2nd → ~1st order on stretched cells]** → cap growth ratio ≤ 1.15 (near
  2nd-order in practice; universally accepted for bluff-body LES). Add MMS test to quantify.
- **[MG convergence rate degrades on graded metrics]** → volume-weighted restriction/prolongation
  (D3); monitor iteration counts vs the uniform baseline; keep `proj_tol` honest.
- **[Timestep set by the smallest cell]** → CFL + `h²/ν` diffusive limit key off the fine core, so
  more steps per flow-through. Accepted: the cost is paid only for the refined region, not
  everywhere. Report effective dt so the user sees it.
- **[~300 call sites churned with NO existing oracle]** → the parity suite CLAUDE.md advertises is
  not live (Context ⚠). Mitigated by **D4 Phase 0**: build the oracle (parity harness or
  golden-master) *before* Phase A, or the "zero behaviour change" claim is unfalsifiable.
- **[windloads silently wrong on a graded grid — the deliverable]** → feed `h_fine` + assert
  building bbox ⊆ uniform core (D6). This is a correctness invariant; without the guard the reported
  Cd/Cl/Cp are wrong the moment any surface cell lands in the transition ramp.
- **[MG coarse-level metrics derived wrong]** → decimate cumulative face coordinates per level, not
  average `dx` (D3 care-item B).
- **[Interior transition plane under-smoothed]** → `ch_gs_band` only sweeps domain edges; instrument
  MG iterations, extend R-B GS to the transition band if it stalls (D3 care-item C).
- **[Grading refines slabs, not spots (D1 trade-off)]** → fine box around one central building is
  acceptable; escalate to nested boxes only if a future case needs disjoint/off-centre refinement.
- **[Scene compatibility]** → `.scn` must persist metrics or the fine-core recipe; legacy uniform
  scenes load as the all-equal special case (back-compatible by construction).
- **[Subtle face-to-centre vs centre-to-centre distance bugs]** → encode both distance families
  explicitly in `MacGrid` (D2); cover with the collapse-to-uniform + MMS tests.

## Migration Plan

0. **Phase 0 (oracle)**: build a parity harness (CPU twins vs GPU at 1e-5) or a golden-master
   snapshot of a canonical run. Nothing else starts until this can fail. (D4)
1. **Phase A**: land metric arrays (init uniform) + metric-aware operators; oracle passes unchanged
   (rollback = revert; no config/GUI surface touched yet).
2. **Phase B**: grid generator + config/GUI fine-core controls; voxelizer/sampler/scene/windloads/
   tracer consumers; `adaptive_dt` → `h_min`; MMS + collapse-to-uniform tests; a graded
   `configs/building.json` variant.
3. **Validation**: run rounded corner at 30 mm and 20 mm; confirm grid-convergence gate (D7);
   document the converged sharp-vs-rounded comparison protocol.
- **Rollback**: Phase A is behaviour-neutral; Phase B is gated behind the fine-core config keys —
  a uniform config path remains available throughout.
- **Coordination**: rebase on master first — commits `d648cb6` (wall floor, scalar `g.h`) and
  `7a8f346` (STEP embedded in `.scn`) have landed. The Open-Question-3 lean (persist the *recipe*,
  regenerate metrics) fits `7a8f346`'s STEP-as-source-of-truth pattern exactly.

## Open Questions

- Should the fine box **auto-track** the gizmo-placed building bbox (+ margin), or be an explicit
  user-drawn box? (Lean: auto-track with an adjustable margin, overridable.)
- What margin around the building does the fine core need so separation/near-wake is captured
  before grading coarsens — 1 radius? 1 building height? (Set empirically during validation.)
- Does `.scn` store the **metric arrays** (exact reproduction) or the **fine-core recipe**
  (regenerate)? (Lean: recipe, matching the STEP-as-source-of-truth pattern; arrays are derived.)
- MG transfer: is volume-weighted restriction enough, or do strongly-graded cases need a Galerkin
  coarse operator? (Decide from measured iteration counts.)
- ~~Near-wall spacing vs wall model~~ **(resolved — defer)**: `bedshear` / the log-law wall model /
  the SEM inlet are **infrastructure not wired into the building configs today** (the building
  workflow uses a uniform inlet and free-slip faces). So do **not** spend Phase-A effort making the
  wall model graded-correct for a path the building workflow doesn't exercise — keep it
  correct-when-uniform and defer graded near-wall spacing until the ABL inflow is actually wired in
  (a separate change). Revisit `h_fine` vs y⁺ then.
