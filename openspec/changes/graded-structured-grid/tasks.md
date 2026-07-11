# Tasks — graded-structured-grid

Ordering follows the migration (design D4): **Phase 0 builds an oracle that does not yet exist**;
Phase A lands the ~300-call-site churn with **zero behaviour change**, proven by that oracle;
Phase B turns on grading and its consumers. Do not start Phase A until Phase 0 can fail; do not
start Phase B until Phase A's oracle is green.

## 0. Phase 0 — Build the oracle (BLOCKS everything; the parity suite is NOT live today)

- [ ] 0.1 Confirm the gap (done in review): no `add_test`/`enable_testing`/gtest target exists;
      CPU twins (`ch_poisson_apply_cpu`, `mac_advect_cpu`, `smagorinsky_nut_cpu`, …) are defined
      but never run
- [ ] 0.2 Choose the oracle form: (a) reconstruct a gtest/`add_test` parity harness driving CPU
      twins vs GPU kernels at rel. max-norm 1e-5 (preferred — restores the CLAUDE.md invariant),
      or (b) golden-master snapshot of `{u,v,w,p}` from a canonical short run
      (`configs/g1_viewer.json`, N steps) for bit-for-bit diffing
- [ ] 0.3 Implement the chosen oracle against the CURRENT uniform solver and confirm it passes
      (establishes the baseline before any edit)
- [ ] 0.4 Pull task 5.1 forward: wire the "grading collapses to uniform" check as the permanent
      form of this oracle

## 1. Phase A — Representation (uniform-initialised, no behaviour change)

- [ ] 1.1 Add per-axis metric arrays to `MacGrid` (`dx[]`, `dy[]`, `dz[]`, `xf/yf/zf`, derived
      centres) and both distance families (face-to-centre, centre-to-centre); keep `pidx/uidx/
      vidx/widx` and layout unchanged
- [ ] 1.2 Add a uniform-metrics constructor path that initialises all spacings to `h` (the
      special case)
- [ ] 1.3 Add the invertible world↔index mapping routine (cumulative-coordinate lookup) with a
      unit round-trip test

## 2. Phase A — Metric-aware operators (parity suite is the oracle)

- [ ] 2.1 Rewrite divergence/gradient (`mac_ops.h`, `channel_ops.cu`) to use local widths /
      centre-to-centre distances
- [ ] 2.2 Rewrite MacCormack advection (`advect.cu`) to use the world→index mapping and
      local-width interpolation weights
- [ ] 2.3 Rewrite diffusion + Smagorinsky (`turbulence.cu`) for per-cell filter width
      `(dx·dy·dz)^(1/3)`
- [ ] 2.4 Rewrite the pressure Poisson operator (`channel_ops.cu`, `project.cu`,
      `channel_pressure.cu`) to the variable-coefficient FV form `Σ_f w_f·(p−p_nb)`; keep it SPD
- [ ] 2.5 Make MGPCG restriction/prolongation volume-weighted (`mgpcg.cu`); coarsen by-2 in index;
      derive **per-level metrics by decimating cumulative face coords** (`xf_coarse[i]=xf_fine[2i]`
      ⇒ `dx_coarse[i]=dx_fine[2i]+dx_fine[2i+1]`), NOT by averaging `dx` (design D3 care-item B)
- [ ] 2.6 Update `adaptive_dt` (`channel_core.cu:169`) to key off `h_min` (min cell dim across the
      metric arrays) for both `dt_adv` and `dt_diff`
- [ ] 2.7 Keep BCs + `bedshear` wall model **correct-when-uniform only** — defer graded near-wall
      work (not wired into building configs; design Open Question 5 resolved: defer)
- [ ] 2.8 Update every operator's CPU reference to the metric-aware form
- [ ] 2.9 Run the Phase-0 oracle with uniform metrics — MUST pass unchanged (Phase A gate)

## 3. Phase B — Grid generation & config

- [ ] 3.1 Implement the fine-core → metric-array generator (uniform core + geometric growth,
      ratio-bounded)
- [ ] 3.2 Add fine-core config keys (box, `h_fine`, growth ratio) to `config.*` and a graded
      `configs/building.json` variant
- [ ] 3.3 Wire the generator into `sim_setup.*` / the Apply path (rebuild at t=0)
- [ ] 3.4 Add GUI fine-core controls in `main_window.*`; optionally auto-track the placed building
      bbox + margin (see design Open Questions)

## 4. Phase B — Consumers

- [ ] 4.1 Update the building voxelizer (`voxelize.*`, `building.cpp`) to use local cell size for
      the wall band; the committed `0.5·h·√2` floor (`d648cb6`) becomes a no-op inside the fine core
- [ ] 4.2 Update the CUDA-GL slice sampler (`slice_field.cu`, `slice_gl.cu`) to sample via the
      world→index mapping
- [ ] 4.3 **windloads (`windloads.cpp`) — the deliverable**: feed `h_fine` (not a stale global `h`)
      and **add a guard asserting the building bbox ⊆ the uniform fine core** so all surface faces
      are `h_fine²` and `A_frontal/A_plan/L_ref`/moment arms stay valid; fail loudly if a surface
      cell lands in the graded transition
- [ ] 4.4 Route `flow_particles.cpp` + `flow_tracers.cpp` locate/sample through the world→index map;
      take domain extent from `xf[nx]`/`yf[ny]`/`zf[nz]`, not `nx·h`
- [ ] 4.5 Persist/reconstruct the grid in `scene_io.*` (fine-core **recipe**, per design lean +
      `7a8f346` STEP-as-source-of-truth pattern); ensure legacy uniform scenes load as the special
      case

## 5. Phase B — Verification tests

- [ ] 5.1 Add the permanent "grading collapses to uniform" regression (all `dx=h` ⇒ matches
      uniform reference ≤ 1e-5)
- [ ] 5.2 Add the MMS Laplacian order-of-accuracy test on a sequence of graded grids
- [ ] 5.3 Add metric-aware parity tests for advection, diffusion, projection on a genuinely graded
      grid
- [ ] 5.4 Add a scene round-trip test (graded save/reload; legacy uniform load)

## 6. Corner-resolution workflow

- [ ] 6.1 Implement the `h_fine` sizing helper (~`r/10`) and surface the `h ≥ r/4` resolution
      floor to the user
- [ ] 6.2 Document + support the matched-A/B and grid-convergence-gate protocol (rounded at 30 mm
      vs 20 mm; compare Cd + peak-suction Cp within the averaging RMS band)

## 7. Validation & docs

- [ ] 7.1 Run the rounded corner at `h_fine = 30 mm` and `20 mm`; confirm the grid-convergence
      gate; record the converged loads
- [ ] 7.2 Produce one converged sharp-vs-rounded comparison on the graded grid as the acceptance
      demonstration
- [ ] 7.3 Update `CLAUDE.md` / `RESEARCH.md` / `PLAN.md` to reflect the graded grid and the
      resolution/convergence protocol
