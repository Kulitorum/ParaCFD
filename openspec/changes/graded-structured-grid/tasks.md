# Tasks — graded-structured-grid

Ordering follows the migration (design D4): **Phase 0 builds an oracle that does not yet exist**;
Phase A lands the ~300-call-site churn with **zero behaviour change**, proven by that oracle;
Phase B turns on grading and its consumers. Do not start Phase A until Phase 0 can fail; do not
start Phase B until Phase A's oracle is green.

## 0. Phase 0 — Build the oracle (BLOCKS everything; the parity suite is NOT live today)

- [x] 0.1 Confirm the gap (done in review): no `add_test`/`enable_testing`/gtest target exists;
      CPU twins (`ch_poisson_apply_cpu`, `mac_advect_cpu`, `smagorinsky_nut_cpu`, …) are defined
      but never run
      → CONFIRMED: the only match repo-wide was a stale comment in `CMakeLists.txt`; no test target.
- [x] 0.2 Choose the oracle form: (a) reconstruct a gtest/`add_test` parity harness driving CPU
      twins vs GPU kernels at rel. max-norm 1e-5 (preferred — restores the CLAUDE.md invariant),
      or (b) golden-master snapshot of `{u,v,w,p}` from a canonical short run
      (`configs/g1_viewer.json`, N steps) for bit-for-bit diffing
      → Chose (a) **and** (b) together: `tools/parity_probe.cu` (minimal `add_test`, no gtest dep)
      does BOTH per kernel — GPU-vs-CPU parity ≤1e-5 (a) + GPU-vs-blessed-golden ≤1e-5 (b).
- [x] 0.3 Implement the chosen oracle against the CURRENT uniform solver and confirm it passes
      (establishes the baseline before any edit)
      → 32/32 checks green (parity at machine-eps 1e-15→0); `ctest -R parity` passes. Goldens
      blessed in `tests/golden/*.f64` against the unmodified scalar-`h` solver.
- [x] 0.4 Pull task 5.1 forward: wire the "grading collapses to uniform" check as the permanent
      form of this oracle
      → The golden gate IS that regression: goldens are the uniform reference; after the metric
      rewrite (all `dx=h`) the same kernels must still match them ≤1e-5. See 5.1.

## 1. Phase A — Representation (uniform-initialised, no behaviour change)

- [x] 1.1 Add per-axis metric arrays to `MacGrid` (`dx[]`, `dy[]`, `dz[]`, `xf/yf/zf`, derived
      centres) and both distance families (face-to-centre, centre-to-centre); keep `pidx/uidx/
      vidx/widx` and layout unchanged
      → `mac_grid.h`: nullable raw-pointer arrays `dxa/dya/dza/xca/yca/zca/xfa/yfa/zfa` (host copy
      for CPU twins, device copy for GPU) + accessors `dx/dy/dz`, `xc/yc/zc`, `xf/yf/zf`, `dxc/dyc/
      dzc` (centre-to-centre), `gapx/gapy/gapz` (OOB-safe), `Lx/Ly/Lz`, `h_min()`. Layout unchanged.
- [x] 1.2 Add a uniform-metrics constructor path that initialises all spacings to `h` (the
      special case)
      → Every accessor falls back to the exact uniform closed form when its pointer is null, so a
      null-metric `MacGrid` IS the uniform special case, byte-identical to the pre-graded grid.
- [x] 1.3 Add the invertible world↔index mapping routine (cumulative-coordinate lookup) with a
      unit round-trip test
      → `locate_frac()` (binary search over cumulative coords) + `grid_fx/fy/fz` (faces) & `grid_cx/
      cy/cz` (centres) with exact `x/h` uniform fast-path. Round-trip test in `parity_probe`:
      err 5.6e-17 on a geometric-growth grid.

## 2. Phase A — Metric-aware operators (parity suite is the oracle)

- [x] 2.1 Rewrite divergence/gradient (`mac_ops.h`, `channel_ops.cu`) to use local widths /
      centre-to-centre distances
      → `div_cell` (both families): each face flux / local `dx(i)`/`dy(j)`/`dz(k)`. Gradient
      (`k_subgrad`/`subtract_gradient_cpu`, both families): / centre-to-centre `dxc/dyc/dzc`;
      Dirichlet-outlet ghost / edge width `dx(nx-1)`/`dx(0)`.
- [x] 2.2 Rewrite MacCormack advection (`advect.cu`) to use the world→index mapping and
      local-width interpolation weights
      → `node_pos` (advect.cu + channel_ops.cu) → `xf/yf/zf`, `xc/yc/zc`; the `trilerp_*` samplers
      (`mac_grid.h` + `channel_bc.h`) route through `grid_fx/cx/…`; `clamp_to_domain` → `Lx/Ly/Lz`;
      channel `cell_of_point_solid` → `floor(grid_fx…)`. All share the invertible map.
- [x] 2.3 Rewrite diffusion + Smagorinsky (`turbulence.cu`) for per-cell filter width
      `(dx·dy·dz)^(1/3)`
      → Filter width `Cs·cbrt(dx·dy·dz)` (both families); diagonal strains / local width, cross
      strains / 2-cell centre span; diffusion Laplacian = Σ of non-uniform `d2_axis` per axis
      (normal axis face-spaced, tangential axes centre-spaced via `gap*`). channel too.
- [x] 2.4 Rewrite the pressure Poisson operator (`channel_ops.cu`, `project.cu`,
      `channel_pressure.cu`) to the variable-coefficient FV form `Σ_f w_f·(p−p_nb)`; keep it SPD
      → `fv_stencil` (project.cu, all-Neumann) + `fv_stencil_ch` (channel, masked + Dirichlet
      outlet): `w_f = 1/(d_centre·width_normal)`, `diag=Σw_f`, `wnb=Σw_f·p_nb`, `Ap=diag·p−wnb`
      (SPD). apply/residual/jacobi/GS + CPU twins converted; collapses to `(count·p−nbsum)/h²`.
      NOTE: `channel_pressure.cu`/`mgpcg.cu` per-level operators call these kernels, so they are
      metric-aware on the finest level; **per-level coarse metrics are task 2.5 (below)**.
- [x] 2.5 Make MGPCG restriction/prolongation volume-weighted (`mgpcg.cu`); coarsen by-2 in index;
      derive **per-level metrics by decimating cumulative face coords** (`xf_coarse[i]=xf_fine[2i]`
      ⇒ `dx_coarse[i]=dx_fine[2i]+dx_fine[2i+1]`), NOT by averaging `dx` (design D3 care-item B)
      → `ChannelMgpcg::build_level_metrics` (channel_pressure.cu): when the finest carries metric
      arrays, each coarse level gets its own device metrics by **decimating** the finer level's
      cumulative face coords (`k_decimate_faces` + `k_faces_to_metrics`), exactly per care-item B.
      Uniform finest (null metrics) skips it (per-level scalar `h=2^l·h` already correct). VERIFIED:
      graded 22³ 2-level solve converges in **13 iters** to relres 4.8e-6 — a strong preconditioner,
      so index-based transfers suffice (volume-weighting is an available rate refinement, not needed).
      NOTE: the M1 cavity `mgpcg.cu` is left uniform (the cavity validation never uses a graded grid).
- [x] 2.6 Update `adaptive_dt` (`channel_core.cu:169`) to key off `h_min` (min cell dim across the
      metric arrays) for both `dt_adv` and `dt_diff`
      → `MacGrid::h_min()` (=`hmin` if set, else `h`); `adaptive_dt` keys both limits off it.
- [x] 2.7 Keep BCs + `bedshear` wall model **correct-when-uniform only** — defer graded near-wall
      work (not wired into building configs; design Open Question 5 resolved: defer)
      → Left `ch_orlanski` (`Uc·dt/h`), `ch_uplane_flux` (`·h²`), `bedshear`, and loglaw-inlet infra
      on scalar `g.h` — correct on uniform, not wired into the building workflow. Documented defer.
- [x] 2.8 Update every operator's CPU reference to the metric-aware form
      → Automatic: GPU kernels + CPU twins share the same `__host__ __device__` helpers (`div_cell`,
      `fv_stencil*`, `smag_nut`, `diffuse_node`, `d2_axis`, the samplers). Twins verified at
      machine-eps parity + golden.
- [x] 2.9 Run the Phase-0 oracle with uniform metrics — MUST pass unchanged (Phase A gate)
      → **GREEN: 42/42.** Null-metric path: 32/32 golden Δ ≤ 2.6e-16. Uniform-metric-ARRAY path
      (device+host arrays, the metric-reading code): 9/9 golden Δ ≤ 7e-15. Full `windcfd-gui` + `ctest`
      build/pass. Behaviour change: none.

## 3. Phase B — Grid generation & config

- [x] 3.1 Implement the fine-core → metric-array generator (uniform core + geometric growth,
      ratio-bounded)
      → `grid_metrics.{h,cpp}`: `graded_axis_faces()` (pure host, testable) + `GridMetrics` owner
      (host arrays for CPU/geometry consumers + device upload for kernels; `host_view()`/`device_view()`;
      `FineCoreSpec`, `generate()`, `uniform()`, `h_min()`). Verified: 26-cell axis, uniform core at
      h_fine, max ratio 1.1500 ≤ 1.15.
- [x] 3.2 Add fine-core config keys (box, `h_fine`, growth ratio) to `config.*` and a graded
      `configs/building.json` variant
      → Keys parsed in `gui/sim_setup.cpp` (`Config` in config.* is dead code at runtime). A nested
      `"fine_core": { enabled, x0..z1, h_fine, growth }` object → `FineCoreSpec`. Opt-in: absent/
      disabled/`h_fine>=voxel_h` ⇒ uniform (byte-identical). Added `configs/g1_viewer_graded.json`.
- [x] 3.3 Wire the generator into `sim_setup.*` / the Apply path (rebuild at t=0)
      → `build_sim`: when graded, `GridMetrics::generate(fc)` into a `std::shared_ptr` in the
      (copyable) `SimRecipe`; `recipe.grid = metrics->device_view()` (GPU core), `info` dims/extent
      from `host_view()`. Ownership: shared_ptr in the recipe keeps the metric arrays alive across
      every recipe copy (recipe_, worker factory, scene_io) → outlives the core's device-view grid.
      Worker gains `metrics_` + `setMetrics()` + `hostGrid()`; host consumers (wind loads, published
      mask/flow grids, flow-through time via `Lx()`) use `hostGrid()`, GPU path uses the device view.
      main_window: `setMetrics(recipe_.metrics)` at spawn; voxelizer uses `host_view()`. Apply re-reads
      the config's `fine_core` → regenerates the graded grid at t=0. VERIFIED headless: graded run
      42×26×8 @ h_fine=0.05 steps 5176 steps, no crash; uniform run unregressed.
      NOTE: GUI dock fine-core controls + Apply-override fine-core fields are task 3.4 (deferred); the
      grid-dims readout still shows uniform-derived dims on a graded grid (cosmetic, 3.4).
- [ ] 3.4 Add GUI fine-core controls in `main_window.*`; optionally auto-track the placed building
      bbox + margin (see design Open Questions)

## 4. Phase B — Consumers

- [x] 4.1 Update the building voxelizer (`voxelize.*`, `building.cpp`) to use local cell size for
      the wall band; the committed `0.5·h·√2` floor (`d648cb6`) becomes a no-op inside the fine core
      → `building.cpp`: world→cell via the grid map — bbox range `floor(grid_fx/fy/fz)`, cell centres
      `g.xc/yc/zc(i)` (were `(i+0.5)·h`). The `0.5·g.h·√2` wall floor stays correct: `g.h==h_fine`
      is the local cell size in the (uniform) core, and at h_fine=0.03 the 0.08 m wall dominates the
      floor (no-op), resolving ~2.7 cells. `voxelize.cpp` (mesh ray-parity voxelizer): the model must
      sit in the uniform core where `xf(i)=i·h+Δ` (per-axis constant Δ), so shifting the placed
      geometry by −Δ makes the uniform sub-lattice's `floor(x/h)` equal the graded core-cell index —
      reuses all the parity code; uniform grid ⇒ `xfa` null ⇒ Δ=0 (byte-identical). VERIFIED
      (`parity_probe`): a square footprint at world (1,1) voxelizes to solid cells centred on (1,1),
      all in the uniform core — the old `floor(x/h)` would mis-place it.
- [ ] 4.2 Update the CUDA-GL slice sampler (`slice_field.cu`, `slice_gl.cu`) to sample via the
      world→index mapping
- [x] 4.3 **windloads (`windloads.cpp`) — the deliverable**: feed `h_fine` (not a stale global `h`)
      and **add a guard asserting the building bbox ⊆ the uniform fine core** so all surface faces
      are `h_fine²` and `A_frontal/A_plan/L_ref`/moment arms stay valid; fail loudly if a surface
      cell lands in the graded transition
      → `g.h` already carries `h_fine`; under the guard the existing scalar math is provably exact
      (proof: within a uniform sub-block the (i+0.5)h arm proxy differs from the true `xc(i)` by a
      per-axis constant that cancels in `r = face − centroid`). Added the guard: track building
      imin..kmax, assert `dx(i)=dy(j)=dz(k)=g.h` over the bbox, else `throw std::runtime_error`. Caller
      (`sim_worker::maybe_compute_loads`) now passes `hostGrid()` (host metric arrays). VERIFIED in
      `parity_probe`: building in the core integrates; building reaching the graded transition is
      rejected. NOTE: a graded run producing loads over a real building also needs task 4.1 (the
      voxelizer's world→index cell mapping); the guard + host-grid caller are complete.
- [ ] 4.4 Route `flow_particles.cpp` + `flow_tracers.cpp` locate/sample through the world→index map;
      take domain extent from `xf[nx]`/`yf[ny]`/`zf[nz]`, not `nx·h`
- [ ] 4.5 Persist/reconstruct the grid in `scene_io.*` (fine-core **recipe**, per design lean +
      `7a8f346` STEP-as-source-of-truth pattern); ensure legacy uniform scenes load as the special
      case

## 5. Phase B — Verification tests

- [x] 5.1 Add the permanent "grading collapses to uniform" regression (all `dx=h` ⇒ matches
      uniform reference ≤ 1e-5)
      → The `parity_probe` golden gate IS this regression (32 kernels null-path + 9 uniform-metric-
      ARRAY checks, all ≤7e-15 vs the scalar-`h` goldens). `GridMetrics::uniform()` also builds
      explicit all-`h` metrics for it.
- [x] 5.2 Add the MMS Laplacian order-of-accuracy test on a sequence of graded grids
      → `mms_laplacian_order` in `parity_probe`: manufactured sin·sin·sin, FV operator vs analytic
      k²·p on the interior; error 1.36e-2→5.14e-3 as h_fine halves ⇒ **order 1.40** (≳1st, as designed).
- [x] 5.3 Add metric-aware parity tests for advection, diffusion, projection on a genuinely graded
      grid
      → 15 kernels (both families) GPU-vs-CPU on a generated graded 23³ grid, all ≤3e-15. Proves the
      metric-reading kernels are correct on non-uniform spacing (not just collapsing to uniform).
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
