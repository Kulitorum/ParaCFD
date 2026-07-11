# CLAUDE.md — ScourProtection

3D morphodynamic simulator (C++/CUDA + Qt6) that ranks 3D-printed concrete scour-protection
shapes by how much sand they trap from tidal flow. Owner: MH (COBOD).

## Read first
- **PLAN.md** — architecture + milestones. Work on exactly one milestone per session, in order.
- **RESEARCH.md** — the physics spec. Every equation/constant comes from there **verbatim**;
  the notes behind it are in `research/*.md`. Never substitute a constant from memory — the
  spec's values were adversarially verified against primary sources.

## Build & test
Verified working on this machine (M0, 2026-07-07) — the VS2022 generator picks the
14.44 host compiler automatically (CUDA MSBuild integration is installed), avoiding the
VS18 trap. Run from any shell (no dev prompt needed):
```
cmake -B build -S . -G "Visual Studio 17 2022" -A x64      # multi-arch default 75;86;89 (Turing/Ampere/Ada)
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure   # gates: ctest -C Release -L gate_M<k>
# Faster single-arch DEV build (RTX 4090 only): add -DCMAKE_CUDA_ARCHITECTURES=89 (a fresh build dir, or
# it stays cached). Shipping/installer builds use the multi-arch default so they run beyond the 4090.
```
Verified working (M1, 2026-07-07). Fast unit suite (kernels + MGPCG + M0), seconds:
```
ctest --test-dir build -C Release -L unit --output-on-failure
```
M1 gate (`gate_M1`, lid-driven cavity 128³ vs Ghia 1982) runs the RTX 4090 for minutes:
```
ctest --test-dir build -C Release -L gate_M1 --output-on-failure
# or directly: build/Release/scour_m1_gate.exe configs/v1_cavity_re100.json
```
M2 gate (`gate_M2`, cylinder Cd/St @ Re=100 & 200 + open-channel BCs, research/11 §4)
runs the RTX 4090 ~25 min (768×384×4 ≈ 1.2M cells, ~37 ms/step, ~42k steps total):
```
ctest --test-dir build -C Release -L gate_M2 --output-on-failure
# or directly: build/Release/scour_m2_gate.exe configs/v3_cylinder.json
# fast iteration (Re=100 only, reduced): build/Release/scour_m2_gate.exe configs/v3_tune32.json
```
- M2 open channel lives in `src/core/fluid/channel_*.{h,cu,cpp}` (channel_bc.h inlet/
  Orlanski/mask sampling, channel_ops.cu masked MacCormack+diffusion+Poisson pieces,
  channel_pressure.cu mask+Dirichlet-outlet MGPCG, channel_mask.cpp cylinder builder,
  channel_core.cu the FluidCore, cylinder.cu drag(CV)+St(FFT), channel_empty.cpp).
- ⚠ The V3 cylinder uses **no-slip** on the obstacle (`SOLID_NOSLIP`): the resolved BL at
  D=32/Re=100 is what generates the Kármán street and Cd≈1.33. PLAN M2's "free-slip
  solids" (the RESEARCH §3 production default for unresolvable sublayers) is also
  implemented (`SOLID_FREESLIP`) and unit-tested, but free-slip produces NO vortex
  shedding and only pressure drag (~1.0) — it cannot meet the St/Cd gate. Both modes are
  a `ChannelBC::solid_mode` toggle. Voxel staircase over-predicts Cd (1.52 at D=32 vs true
  1.33-1.37); the ±10% gate band (→1.54) accommodates this — the cut-cell/variational
  projection upgrade (RESEARCH §3, Batty 2007) is deferred.
- M2 GPU-vs-CPU kernel parity + ChannelMgpcg convergence: `tests/test_channel_kernels.cu`.
M3 gate (`gate_M3`, V2 log-law: flat-bed periodic channel u*/κ + Jarrin SEM inlet TI)
runs the RTX 4090 ~2 min (two 1-D periodic cases ~40 s each + a 100×100×50 SEM run):
```
ctest --test-dir build -C Release -L gate_M3 --output-on-failure
# or directly: build/Release/scour_m3_gate.exe configs/v2_channel_loglaw.json
```
- M3 lives in `src/core/fluid/`: **bedshear.{h,cu}** + bedshear_ops.h (log-law wall
  function: rough + Christoffersen–Jonsson transitional u*, EMA time-filter, 3×3 τ_b
  smoothing, [0.5×,2×] rate-limit, cavity four-branch rule, bottom momentum sink);
  **channel_periodic.{h,cu}** + periodic_ops.h (V2 driver: body force g_x=u*²/h + PI
  mass-flux controller, mixing-length ν_t=l²|dU/dz|, l=κz√(1−z/Lz)); **sem_inlet.{h,cu}**
  (Jarrin 2006 SEM, Nezu–Nakagawa R_ij, Cholesky colouring); **precursor.{h,cpp}** (5–10 Hz
  inlet-plane record/replay, looped linear interp); **channel_sem.cpp** (SEM open-channel TI
  driver); inlet_fluct.h + plane_ops.{h,cu} (the ChannelFluidCore SEM hook — guarded, so
  M2 stays byte-identical). GPU-vs-CPU parity: `tests/test_bedshear.cu`, `tests/test_m3.cu`.
- ⚠ M3 CLOSURE DEVIATION (flagged): the V2 periodic channel uses the standard open-channel
  Prandtl mixing length l=κz√(1−z/Lz), NOT RESEARCH §3.4 Smagorinsky (Δ=h). Rationale: a
  grid-tied Smagorinsky length (0.11·5cm≈5.5 mm) cannot carry a 5-m log layer; the mixing
  length makes the body-forced steady state an exact log law (research/04). Only κ=0.40
  (spec) is used — no new constant. Smagorinsky remains the default for the obstacle/SEM
  runs. The periodic-channel flow is x/y-homogeneous ⇒ 1-D, so advection and projection are
  exact no-ops there (verified: max|∇·u|~1e-13); the full 3-D turbulent precursor is M6.
- ⚠ For d50=0.2 mm the bed is transitional (ks+≈12), so the wall model's CJ branch yields
  u*≈0.0333 vs the rough-z0 analytic 0.0345 (3.4% low, inside the 10% gate) — this is the
  physically-correct regime, not a miss. d50=1.0 mm is rough (ks+≈74) and matches to 0.4%.
M4 gate (`gate_M4`, V5 settling column + V6 Rouse profile) is FAST (tiny 1-D columns, <1 s
total on the RTX 4090) — kept out of the fast unit suite only as a labelled gate. Verified
PASS (M4, 2026-07-08): settling L2=0.03%, w_s exact, sand-mass (susp+bed) 0.00000%; Rouse
R_fit=0.9071=R_target (0.0% err, R²=1.0).
```
ctest --test-dir build -C Release -L gate_M4 --output-on-failure
# or directly: build/Release/scour_m4_gate.exe configs/m4_sediment.json
```
- M4 suspended sediment lives in `src/core/sediment/`: **sed_physics.h** (pure __host__
  __device__ inlines — Soulsby w_s, D*, Soulsby–Whitehouse θ_cr, τ_cr, van Rijn 2019 pickup
  E with the f_D=1/θ′ damping, hindered (1−c)^4.7, deposition D=w_s·c_b; every derived constant
  computed from config (ρ,ρs,ν) at runtime, g=9.81); **suspended.{h,cu}** (kernels + CPU twins:
  effective-w_s, MacCormack-SL scalar advection with settling folded into vz, CONSERVATIVE
  van-Leer flux advection, ν_t/σ_s diffusion, flux-form bed exchange); **sed_validation.{h,cu}**
  (the V5/V6 gate drivers). GPU-vs-CPU parity + sed_physics tables: `tests/test_suspended.cu`.
- ⚠ M4 SCHEME DEVIATION (flagged): the concentration transport on the sand-budget path uses a
  CONSERVATIVE flux-form (TVD van-Leer, `suspended_advect_cons`), NOT the RESEARCH §6.2 "settling
  folded into the MacCormack advection velocity" (`suspended_advect`, also implemented + parity-
  tested). Reason: MacCormack semi-Lagrangian is 2nd-order accurate but NOT mass-conservative
  (~0.14% drift floor, resolution-independent — measured), which is unacceptable for a sand-MASS
  trapping tool. The conservative form closes the suspended+deposited budget to machine precision
  and, being flux-form, correctly depletes the surface so Rouse develops from any state (SL has a
  spurious uniform fixed point). No physics constant changed. Deposition stays a flux BC and is
  NEVER thresholded (θ_cr gates erosion/pickup only) — RESEARCH §5 asymmetry, the product core.
- ⚠ RESEARCH §10 lists the settling column as c(t)=c0·e^(−w_s t/h); that is the WELL-MIXED lumped
  model. A resolved diffusion-free still-water column is PLUG FLOW (research/11 §6: "front advects
  down at exactly w_s; mass conservation"), so the gate checks the L2 of the numeric vs the exact
  translated slab (settling velocity + advection accuracy) and total-sand conservation — both
  reported; the exponential is printed as a diagnostic only.
- Rouse gate uses β=1 (σ_s=1 ⇒ ε_s=κ·u*·z(1−z/h)) to match the classical analytic (R=0.91 at the
  benchmark ν=1.05e-6/ρ=1025); production suspended diffusion uses σ_s=0.7 (RESEARCH §6.2). The
  Dirichlet pin of c=c_a at a=0.05h is used in THIS validation only (PLAN M4).
M5 gate (`gate_M5`, bed state + morphodynamics: a-i deposition column, a-ii prescribed bedload
vs 1D Exner, b sand-pile avalanche, c three-reservoir mass conservation over 1e4 steps,
d Soulsby–Whitehouse threshold). The (c) 1e4-step run takes ~2–3 min on the RTX 4090; the rest
are seconds. Verified PASS (M5, 2026-07-08): a-i rate err 0.0%/mass 0.0%; a-ii upwind Exner L2
0.31%; b pile 45°→30.41° conv, mass 5e-15; c mass 0.0%; d onset 0.0493 vs θ_cr 0.0492 (0.25%).
Launch in background + poll the log — do NOT run it as one blocking foreground call:
```
ctest --test-dir build -C Release -L gate_M5 --output-on-failure
# or directly: build/Release/scour_m5_gate.exe configs/m5_morpho.json
```
- M5 lives in `src/core/sediment/`: **sed_physics.h** (M5 inlines added — c_pack=0.64, normalized
  fill F, slope-corrected θ_cr Soulsby eq.80a with the 0.1 overhang clamp, Winterwerp u_lift,
  Wong–Parker & Engelund–Fredsøe Φ, q_b, bedload layer thickness); **bedstate.{h,cu}** (kernels +
  CPU twins: f_pack reconstruction from per-column grain thickness G, 3×3×3 box-filtered interface
  geometry n_b=−∇F̃/|∇F̃| + A_b, per-column β/up-slope, bedload flux, upwind-along-transport-direction
  divergence, the Exner column update with burial + MORFAC + |Δz_b|≤0.05h limiter); **avalanche.{h,cu}**
  (8-neighbour mass-conserving sand-slide, 32° trigger / 30° relax); **morpho_validation.{h,cu}**
  (the a–d gate drivers). GPU-vs-CPU parity + physics tables: `tests/test_bedstate.cu`.
- Bed model is the OPEN-BED single-interface Exner-equivalent reduction of the f_pack field
  (RESEARCH §7 eq.6): the morphodynamic state per column is grain thickness G=Σf_pack·h,
  z_b=G/c_pack, and f_pack is the recompacted view. Overhang/layered-column avalanching is the
  general case (M6/M8); the gate exercises the single-layer reduction. The box-filtered interface
  normal/area IS computed and used (the general overhang-capable form), verified flat→β=0 and
  20°-ramp→β≈20°.
- ⚠ M5 divergence choice (flagged, not a deviation): bedload ∇·q_b is upwinded ALONG THE TRANSPORT
  (flow/τ) DIRECTION per research/03 ("upwind of q_b along the transport direction"), NOT a
  sign-of-flux Godunov split. For the aligned load q=q_b·t̂ they coincide; the distinction only
  matters for the a-ii test (fixed +x current, q0·sin magnitude), where directional upwind is the
  clean backward difference (0.31% L2) — a Godunov split hits O(1) errors at the sin sign-changes
  (√(1/nx)≈4%). Single-valued per face ⇒ conservative (verified in c: mass 0.0% over 1e4 steps).
- ⚠ M5 a-i note: the deposition-rate target w_s·c0 (PLAN M5) is the DILUTE limit; the exact rate
  onto a rising bed is w_s·c0/(1−c0/c_pack). The gate uses a dilute, sub-cell bed (c0=0.004) so
  w_s·c0 is accurate to <0.5% (rate err measured 0.0%). Multi-cell bed growth + suspended-sand
  BURIAL (covered cells → bed grain) are exercised + verified mass-conserving by the (c) gate.
  Deposition is NEVER thresholded; θ_cr gates erosion + bedload only (RESEARCH §5 asymmetry).
- MorphoParams carries the toggles the harness needs: `bedload_formula` (0=Wong–Parker default,
  1=Engelund–Fredsøe — the Roulund 2005 form for M6/M7), `morfac` M, `dz_limit_frac` (0.05·h),
  and per-process on/off. Sediment consumes only {u,v,w,ν_t,τ_b}: τ_b (grain-skin, per-column
  vector) is an INPUT here (synthesised in the gate; wired from `fluid/bedshear.cu` in M6).
G1 GUI (`scour-gui`, Qt 6.11.1 + GL 4.3 core slice viewer). Verified working
(G1, 2026-07-08). Qt is ON by default now (`-DSCOUR_ENABLE_QT=ON`); reconfigure once so
CMake picks it up, then build the target. Running the exe needs the Qt bin dir + the
`platforms` plugin on PATH (the VS build does not deploy them):
```
cmake -B build -S . -G "Visual Studio 17 2022" -A x64 -DSCOUR_ENABLE_QT=ON   # arch default 75;86;89
cmake --build build --config Release --target scour-gui
# run (PowerShell): put Qt on PATH so the DLLs + qwindows platform plugin resolve
$env:Path = "C:\Qt\6.11.1\msvc2022_64\bin;" + $env:Path
$env:QT_QPA_PLATFORM_PLUGIN_PATH = "C:\Qt\6.11.1\msvc2022_64\plugins\platforms"
build/Release/scour-gui.exe configs/g1_viewer.json            # interactive
build/Release/scour-gui.exe configs/g1_viewer_full.json       # full 10x10x5 m @ h=5cm (4M cells)
```
gate_G1 (automated half): launches offscreen-capable, runs configs/g1_viewer.json for
10 s, exits 0 on scripted auto-close. ctest inherits your shell env, so set the Qt PATH
first:
```
ctest --test-dir build -C Release -L gate_G1 --output-on-failure
# directly: build/Release/scour-gui.exe configs/g1_viewer.json --autoclose-ms 10000
```
- G1 viz add-ons (2026-07-08, GUI-only — all on the main thread, decoupled render/step untouched;
  verified 59–68 fps on the full 4M-cell fluid viewer AND the seabed scenario, gate_G1 + unit suite
  green). Controls live in a **left-side settings dock** (`main_window.cpp buildControlDock`): a
  scrollable column of QGroupBox sections — **Simulation (live — no reset)** (Play/Pause + Step, Input
  speed U, inlet profile, sediment BC, Add sand), **Domain & resolution (Apply resets)** (editable
  domain + voxel/cell size h + the cell-count readout + Apply — the ONLY control that resets to t=0), and
  **Visualization** (field, slice plane, plane pos, the layer toggles + arrow controls + **Fast sim**
  throttle). The idea mirrors cobod-slicer's grouped left panel (NOT its styling); the old crowded top
  toolbar + View menu are gone. File menu has **Recent Files**:
  - **Fast sim — graphics-update throttle (2026-07-09, GUI-only)**: a **Visualization** combo (Off (live)
    / 0.5 s / 1 s / 3 s / Only when paused) that trades render smoothness for sim throughput when you are
    away from the screen. The ONLY "graphics" work on the sim's critical path is the worker's per-step
    `publish_display()` (a device snapshot + `cudaStreamSynchronize`) + the arrow host D2H; the throttle
    skips both between wall-clock windows. `SimWorker::setDisplayInterval(s)` gates publish_display/host-flow
    to ≤ once per s seconds, but ALWAYS publishes on pause (idle branch) and on a manual Step (req>0) so a
    paused/stepped view is never stale. `MainWindow::setDisplayThrottle` also widens the ~60 Hz repaint timer
    to match (clamped ≤ 0.5 s so a paused/interactive frame still appears promptly). Off (default, 0) ⇒
    publish every step ⇒ **byte-identical to G1** (gate_G1 unaffected). Re-applied across a rebuild/restore
    (`spawnWorker`). Headless smoke `--display-interval <s>`; measured **+7% steps/wall-s on the 4M-cell V000
    seabed offscreen** (worker-side only — a windowed run adds the paintGL/auto-range/arrow-contention saving
    on top). `main_window.cpp (fast_sim_box_, setDisplayThrottle, repaint_timer_)`; `sim_worker.cpp` run loop.
  - **Live sim-fps readout (2026-07-09)**: the lower-left status label appends the measured **simulation
    throughput** ("… dt = 9.55 ms   11.2 sim fps") — steps/second computed on the WORKER across each 8-step
    stats window (paused wall-time excluded via an idle re-baseline), lightly EMA-smoothed, passed as a 4th
    arg on `SimWorker::stats`. It is the SIM rate (rises under Fast sim), distinct from the **render fps** on
    the status bar's right — the two side by side confirm the throttle speedup at a glance. `sim_worker.cpp`
    run loop; `main_window.cpp` stats handler.
  - **Editable domain + resolution → Apply** (2026-07-08, GUI-only): the Simulation group has
    Lx/Ly/Lz (`QDoubleSpinBox`, m) + a single **voxel / cell size h** box (h IS the voxel size — one
    uniform grid: fluid cell = model/bed voxelization spacing; there is NO separate voxel control), a
    live **"→ nx×ny×nz = N cells (~MB)"** readout (recomputed from the SAME `grid_dims_for()` the
    builders use, so it is truthful; guards Apply off above `kMaxCells`=40 M to avoid an OOM), and an
    **Apply** button that rebuilds + **resets the whole sim at t=0** on the new grid. Apply reuses the
    EXISTING setup path (`sim_setup` `build_sim`/`build_seabed_sim`) with a new `GridOverride` that
    replaces ONLY the domain + h — every physics/scenario constant (U, Re, obstacle, d50, MORFAC, α,
    Cs, ν, structure STEP …) is preserved. Three rebuild paths (`main_window.cpp applyGrid`): (a) a
    seabed scenario FILE re-voxelizes its structure + re-inits the sand from config (`build_seabed_sim`);
    (b) a **runtime-CONVERTED seabed** (a fluid viewer turned into a sediment run via "Add sand",
    tracked by `added_sand_m_>0`) re-voxelizes the **loaded model (or the config obstacle) as the rigid
    STRUCTURE + re-applies the stored sand depth** at the new grid (`build_seabed_from_structure`) — so a
    resize keeps the model+sand instead of reverting to the config's default cylinder (the old
    add-sand-provenance gap); (c) a plain fluid viewer keeps its config obstacle and re-voxelizes its
    loaded STEP model as the obstacle. The rebuild uses the SAME full-teardown-then-respawn path as
    "Load seabed scenario…" (`teardownWorkerForReload` stops+joins+deletes the worker, THEN the new
    core/engine is built while NO worker thread is live and with no GL — race-free — and a fresh
    `SimWorker` re-allocates every grid-sized device snapshot; the viewer's `setInfo`/`setWorker`/
    `setBedSurface` re-size the slice, bed surface, voxel overlay + arrow flow-field for the new dims,
    so no stale-size buffer is ever read). `recipe.source_config`/`recipe.is_scenario` carry the
    provenance so Apply knows which builder + config to re-run. Coarser (bigger h) = fewer cells =
    faster to iterate: observed ≈59 steps/s at h=0.10 (0.5 M cells) vs ≈7 at h=0.05 (4 M) on the V000
    seabed; Apply-at-same-h matches the direct-build rate (zero rebuild overhead). Headless smoke:
    `--apply-h <h>` (+ optional `--apply-domain Lx Ly Lz`) drives the exact Apply path after startup
    (mirrors `--load-step`/`--voxelize`). `main_window.cpp applyGrid/applyGridValues/syncGridControls/
    updateGridReadout`; `sim_setup.h GridOverride/grid_dims_for`.
  - **Input speed U (live current)** (2026-07-08, GUI-only): a `QDoubleSpinBox` (m/s) in the Simulation
    group that sets the inlet current **LIVE — no reset**, so you can test scour/wake at different
    currents interactively (the developing flow evolves toward the new U). It drives
    `SimWorker::setInletSpeed` → applied on the WORKER thread at the top of its loop (atomic pending
    value + dirty flag, same pattern as `requestRebuild`) via `ChannelFluidCore::set_inlet_speed`, which
    updates `bc_.U_inlet` + `bc_.Uc` (both read every step: the Dirichlet inlet + Orlanski convective
    speed) — the SAME BC the M2 gate drives, so no gate impact (the gates never call it). The knob is
    ALSO folded into `GridOverride.U` so an **Apply** rebuild keeps the chosen current (build_sim/
    build_seabed_sim override the config U when `ov.U>0`); it is re-applied after an obstacle re-inject
    so a `--voxelize` toggle doesn't revert it (`inlet_override_` latches). `SliceViewer::setReferenceU`
    tracks the fixed colour range + arrow speed fallback when auto-range is off. ⚠ nu is unchanged (the
    fluid viewer's artificial ν=U₀·D/Re is NOT rescaled), so raising U raises the effective Re — the
    point of "different currents", but a big jump can destabilise the coarse cylinder grid (seabed
    nu_fluid is U-independent → fine). Headless smoke: `--set-u <U>` sets it live after startup
    (`main_window.cpp setInputSpeed`; verified: doubling U halved the adaptive dt — the flow genuinely
    sped up — sim stable, exit 0). `main_window.cpp buildControlDock (u_spin_)`.
  - **MORFAC — live bed speed-up (2026-07-09, GUI-only)**: a `QDoubleSpinBox` ("MORFAC (bed speed-up)", 1–50×)
    in the Simulation group that retunes the morphological acceleration **LIVE** so the bed catches up faster
    without a reset — the answer to "make the scour develop sooner". It drives `SimWorker::setMorfac` →
    applied on the WORKER thread (atomic value + dirty flag, same pattern as the other live knobs) via
    `SeabedMorpho::set_morfac` (which sets both `MorphoParams.morfac` + the suspended `sp_.morfac`). The apply
    is guarded `if (morpho_ && morfac_dirty_.exchange(false))` so a value dialed in on a plain fluid viewer
    stays **pending** and lands the moment a bed attaches via "Add sand". The control **reflects the loaded
    scenario's MORFAC** on worker spawn (`scen.params.morfac`, signals blocked so it shows without re-pushing).
    A no-op without a bed (like Sediment BC); MORFAC multiplies morphological time ONLY, so the flow clock `t`
    is unchanged while the bed ages M×t (`bedstate.cu` `M*dt`). ⚠ Physically valid to ~10 steady / ~5 reversing
    (RESEARCH §7); above that the flow↔bed decouple and the per-step |Δz_b|≤0.05h limiter clips, so higher is
    not proportionally faster (the `clip_fraction` diagnostic tracks this). No gate impact (the M-gates set
    morfac in their own configs; the setter is GUI-only). Headless smoke `--morfac <M>` (verified on the
    add-sand seabed conversion: at step 600 z_b range 0.062 m at M=12 vs 0.023 m at M=3 — the bed genuinely
    evolves faster — budget err ~1e-6, stable, exit 0). `main_window.cpp (morfac_spin_, setMorfacValue)`;
    `sim_worker.{h,cpp} (setMorfac)`.
  - **Inlet profile: uniform ↔ boundary-layer (log-law) (2026-07-08)**: a Simulation-group combo (live) that
    switches the inlet between top-hat and a log-law BL profile. **Fixes the top-hat leading-edge scour**:
    with a uniform inlet the near-bed velocity at column i=0 is the full free-stream U (no BL), so the wall
    model probes ≈U just above the sand → maximal τ_b exactly at the sediment's upstream edge → it washes
    away far too fast (i=0 is pinned to the Dirichlet inlet so it never develops a BL; downstream a BL forms
    and shear drops). The log-law inlet `u(z)=(u*/κ)ln((z−bed_datum)/z0)` (referenced to the bed top via the
    new `ChannelBC::bed_datum`, z0=d50/12) makes near-bed u→0 at the sand and is flux-matched to U over the
    fluid depth (`loglaw_ustar_for_U`), so the leading-edge shear matches downstream. **Seabed runs default to
    log-law** (config key `inlet_profile: "loglaw"|"uniform"`; "Add sand" too); fluid viewers stay uniform (M2
    cylinder unchanged — the uniform branch of `channel_inlet_u` is untouched, `bed_datum` defaults 0 so the M3
    log-law is byte-identical, gate verified). Live via `ChannelFluidCore::set_inlet_profile` (recomputes u*;
    `set_inlet_speed` also re-derives u* in log-law mode) → `SimWorker::setInletProfile` (worker-thread apply +
    re-apply after rebuild) → the combo, synced from `recipe.bc.inlet_mode` in `syncGridControls`. Verified on
    V000: near-bed inlet u 0.50→0.369 m/s (τ_b −45%), stable, 60 fps. `sim_setup.cpp assemble_seabed`;
    `channel_bc.h channel_inlet_u/loglaw_ustar_for_U`; `main_window.cpp (inlet_profile_box_)`.
  - **Add sand (2026-07-08, GUI-only; flow-PRESERVING as of the live-conversion update)**: a "Sand depth"
    (m) box + **Add sand** button in the Simulation group that fills the lower N m of every non-rigid column
    with sand (rigid structure/obstacle voxels stay solid) and continues as an erodible-seabed sediment run —
    so even a plain fluid viewer with a voxelized STEP obstacle becomes a live scour sim. ⚠ It **PRESERVES the
    developed flow — NO t=0 reset**: you spin the flow up to steady state, THEN add sand and the wake keeps
    going. It reads the current RIGID structure as the **kind-2 cells of the live flow mask**
    (`SimWorker::copyMask`, excludes existing sand kind-1; falls back to `recipe.base_solid` if no mask yet),
    then `sim_setup.cpp build_seabed_from_structure(..., build_core=false)` assembles JUST the morpho engine +
    initial sand+structure mask + recipe on the SAME grid (via the shared `assemble_seabed`, which skips the
    flow-core construction when `build_core=false`) — **no new core, no teardown**. `main_window.cpp addSand`
    hands the engine + mask to `SimWorker::requestSeabedConversion`, which on the worker thread installs the
    mask via **`ChannelFluidCore::update_solid`** (untouched cells keep {u,v,w,p}; now-solid sand cells are
    zeroed — the SAME in-place re-mask the moving bed uses), enables the bed-inlet mask, switches the inlet to
    the log-law BL, attaches the engine, and starts morphology `spinup` (200) steps later so the re-masked
    flow settles. Sim time/step count are NOT reset. A near-surface "proximity" velocity ramp is deferred (the
    update_solid face-zeroing + re-projection handles the interface; add the ramp if a shock appears). Sediment
    params carry from the current scenario (`SeabedScenario::params` → `last_sp_`) or default to the nominal
    demo values. PRE-CALIBRATION / qualitative. Headless smoke `--add-sand <m>` (verified: g1_viewer_full +
    `--add-sand 0.3` → live conversion, morphology from step 201, at step 400 z_b∈[0.281,0.303] (scour
    signature), **sand mass err 5.2e-07**, 528 steps, exit 0). `main_window.cpp addSand`; `sim_worker.cpp
    requestSeabedConversion/apply_seabed_conversion`; `sim_setup.cpp build_seabed_from_structure/assemble_seabed`.
  - **Recent Files** submenu (last 8 STEP paths **and saved scenes**, MRU-first, persisted via
    `QSettings("COBOD","ScourProtection")` key `recentStepFiles`; missing files pruned on open;
    **Clear Recent**). A recent entry dispatches by extension: `.scn` → full scene restore, else a STEP
    load. Updated on every STEP open + on a manual **Save Scene** (numbered `.scn` checkpoints are NOT
    added). `main_window.cpp addRecentFile/rebuildRecentMenu`.
  - **Scene save / restore (checkpoints) (2026-07-08, `src/gui/scene_io.{h,cpp}`)**: File menu **Save
    Scene As…** (Ctrl+S), **Load Scene…**, and **Auto-save checkpoints** (Off / 500 / 1000 / 2000 / 5000
    steps). A **`.scn`** is ONE self-contained, self-describing binary container (magic `SCOURSCN` + version
    + nlohmann-JSON metadata + TLV field blobs) holding EVERYTHING to reconstruct AND continue a run: the
    setup (`SimRecipe` grid/BC/params/config-obstacle/provenance), the loaded STEP model (embedded triangle
    mesh, display-only), and all LIVE field DATA at the saved step — velocities+pressure (`u/v/w/p`), the
    flow solid mask, and (seabed) the bed grain `G`, suspended `c`, EMA-filtered near-bed shear `emax/emay`
    + rigid `structure`. Big arrays are **float32** (user-chosen; ≈ half size, qualitatively lossless);
    G/emax/emay f64 (tiny), masks u8. One manual **Save Scene As** writes `<name>.scn` (→ Recent Files); an
    auto-save writes `<name>.<step>.scn` every N steps (a restart point; not Recent) — both full, independently
    loadable. **Load Scene** offers a restart-point picker (`list_checkpoints` scans the `<name>.*.scn`
    siblings) and, on pick, tears the sim down and rebuilds it through the SAME build path
    (`make_core` + `SeabedMorpho`), then INJECTS the saved fields — so the flow + bed resume EXACTLY where
    they were (not t=0), morphology continuing immediately (`spinup 0`). **What's saved is the STATEFUL
    device state only** (per-step scratch — `uB/vB/wB`, ν_t, near-solid, all the seabed geometry temporaries
    — is recomputed on the next step; the wall-model EMA `emax/emay` IS saved as it carries `t_avg=2 s`
    memory). New plumbing: `ChannelFluidCore::copy_state_host/load_state_host` + `bc()`/`bed_inlet_mask_on()`;
    `SeabedMorpho::save_state/load_state` + `params()` (H2D also resets `G_last=G`); the worker gathers the
    state **on its own thread (race-free)** and emits `checkpointReady(std::shared_ptr<CheckpointState>, tag)`
    (queued → the GUI writes on the main thread, no big copy crosses the connection), via
    `SimWorker::requestCheckpoint`/`setAutosaveInterval`/`primeCounters`; `MainWindow::saveSceneNow/
    loadSceneFromPath/restoreScene/onCheckpointReady`. Headless smoke: `--save-scene <path>` / `--load-scene
    <path>` / `--autosave <n>` (verified: fluid save@417→restore continues to 945; seabed convert→save→restore
    resumes the scoured bed z_b∈[0.21,0.35] with mass err ~1e-6; auto-save series `rt_seabed.{1800,2100,…}.scn`;
    numbered-checkpoint load resumes at exactly its step). Qt-free `scene_io` (std streams + nlohmann) ⇒ no Qt
    in the physics libs. **`scene_io.cpp` is a `scour-gui` source** (CMake).
  - **Show slice / Show model / Show solid voxels / Show sediment voxels / Show bed** toggles
    (independent; default ON) — `SliceViewer::setShowSlice/setShowModel/setShowSolidVoxels/
    setShowSedimentVoxels/setShowBed` gate the coloured slice plane, the STEP-mesh, the two voxel
    overlays and the seabed surface in `paintGL` (hiding the slice skips its CUDA fill; hiding the bed
    skips its per-paint fill + range work).
  - **Clip plane — see inside hollow structures (2026-07-09, GUI-only)**: a **Clip plane (see inside)**
    dock group (Enable / Orientation combo / Position slider / Flip) that hides SOLIDS on the camera side
    of a movable plane so the interior of a hollow structure is exposed. It cuts the STEP mesh, the voxel-
    solid staircase (both sediment AND rigid kinds) and the seabed surface, but **NEVER the flow slice or
    the arrows** — so the flow field inside the revealed cavity stays on screen ("Solids only", MH's choice).
    Two orientation modes: axis-aligned **X/Y/Z** (the slider shifts the plane along that axis, auto-oriented
    to hide the camera-side half via `sign(eye_a − p)`) and **Face camera** (normal = the view direction, the
    slider pushes it into the scene along the line of sight, swept ±½·diagonal about the camera target); a
    **Flip** toggle swaps the hidden half. Implemented with **`gl_ClipDistance[0]`**: `kVert`/`kMeshVert`
    write `dot(pos, uClipPlane)` every frame but it only cuts while `GL_CLIP_DISTANCE0` is glEnabled — enabled
    per SOLID draw (mesh, voxels, bed) and left OFF for the slice/arrows/box/axes, and always disabled before
    the QPainter overlay. A translucent cyan quad + outline (`drawClipPlaneViz`, blended, no depth write)
    shows where the cut lies. `SliceViewer::setClipEnabled/setClipMode/setClipFraction/setClipFlip` +
    `computeClipPlane` (world-space `(n,d)`, kept half `dot(x,n)+d≥0`); `Camera::forward()/targetPoint()`.
    No physics touched (gate_G1 + 93/93 unit PASS). Headless smoke `--clip <x|y|z|camera>,<frac>[,flip]`
    (verified windowed on the 4M-cell full viewer: all three modes render, ~45–52 fps, exit 0).
  - **Model-placement gizmo — move/rotate/scale a loaded STEP model, then voxelize where placed (2026-07-09,
    GUI-only; UNIFIED manipulator 2026-07-09)**: a **Model placement (gizmo)** dock group (an **Enable
    manipulator** checkbox + **Reset placement** + a live pos/rot/scale readout) with a **single 3D manipulator
    that draws ALL handles at once** — 3 translate arrows, 3 rotate rings and 3 scale cubes (X red / Y green /
    Z blue) plus a grey uniform-scale centre cube — and **whichever handle you grab picks the operation
    dynamically** (no mode switching; all 9 + uniform functions live simultaneously). Handles are laid out at
    distinct radii (arrow tip 0.70·g < ring 1.00·g < scale cube 1.25·g, the cube OUTSIDE the ring) so
    nearest-in-screen picking disambiguates cleanly; the scale cubes are drawn 1.5× larger with a 1.5× pick
    radius (easy to grab) and, being compact point handles, take pick priority over the shaft/rings. Overall
    gizmo size `g = 0.16/3·distance` (1/3-scale; tune `gizmoSize()`/`kScalePos`/`kCubeHalf` in `slice_viewer.cpp`).
    Left-drag a handle to move/rotate/scale; a left-drag away from any handle still orbits the camera (the
    press pick-tests handles first and only grabs on a hit). The active/hovered handle is a `(op,axis)` pair
    (`op` 0 translate / 1 rotate / 2 scale / 3 uniform), highlighted yellow. The model's world transform is the **single shared placement** the DISPLAY and
    the VOXELIZER use, so the next **"Apply" (Domain & resolution) re-voxelizes the obstacle exactly where the
    model was placed** — MH's ask. The placement PERSISTS across the Apply rebuild (captured before teardown,
    re-applied after `setMesh`) and across scene save/restore. The core change is `ModelPlacement` upgraded
    from translation-only to a full **affine** `world(v)=M·v+t` (M = row-major 3×3 rotation·scale, default
    identity ⇒ `place_model_on_bed` + the translation-only path are byte-identical; `voxelize.cpp` applies the
    affine and place-corrects the mesh-volume gate by `|det M|` so the vol-err stays meaningful under scale).
    The gizmo is enabled whenever a **gizmo-editable model** exists — a loaded STEP in a fluid viewer AND a
    **seabed structure** (a real scenario-FILE structure, a fluid viewer CONVERTED via "Add sand", or a
    restored seabed scene). The seabed structure is seeded into the gizmo (`setMesh` + `setModelPlacement`,
    `mesh_override_` OFF ⇒ `hasModelPlacement()` true ⇒ editable) rather than pinned by a fixed override, so
    it can be **repositioned before running**; Apply then re-voxelizes it into the sand at the gizmo placement
    (place-before-run, full restart at t=0). The placement survives Apply for ALL paths: it is captured from
    the VIEWER (not `model_mesh_`, which a scenario structure lacks) and threaded into `build_seabed_sim`
    (new `structure_place` override), `build_seabed_from_structure` (converted), and `setModelXform`/re-voxelize
    (fluid); the **Enable manipulator** state (`gizmo_on_`) persists in the viewer across the rebuild (no
    combo/mode to restore). `updateGizmoUi` gates the group on `hasModelPlacement()` alone (NOT on `seabed_`). ⚠ On a domain-CHANGING Apply the placement is world-fixed (not re-centred to the new domain) —
    consistent with "voxelize where placed"; use Reset placement to re-centre. A resolution-only Apply (same
    domain) reproduces the exact placement, so it is byte-identical to the old centre-on-bed path when unmoved.
    Implementation: `SliceViewer` owns the TRS (`gz_pivot_/gz_t_/gz_rot_(QQuaternion)/gz_scale_` about the
    model bbox centre) + `modelPlacement()/setModelXform()/setModelPlacement()` (the last decomposes an affine
    M=R·S back to rot/scale via column norms); picking is screen-space distance to the projected handle
    segments/rings; Move projects the mouse delta onto the screen-axis, Rotate intersects a world pick-ray with
    the rotation plane for a sign-correct angle, Scale is exp(screen-axis delta) / radial ratio. All GL is
    main-thread flat-line draws reusing `prog_` (depth test OFF so the manipulator stays grabbable), decoupled
    render/step untouched. `SliceViewer::setGizmoMode/resetModelPlacement/modelPlacement/modelXform/mouseRay/
    pickGizmo/dragGizmo/drawGizmo`; `MainWindow::updateGizmoUi/nudgeModelPlacement` + the dock group; scene_io
    persists the affine `m[9]` (absent ⇒ identity, back-compatible). Headless smoke `--model-place
    dx,dy,dz[,rz[,scale]]` (verified: V000 placed t=(6.5,5.8,1.3) rz=25° scale=1.2 → mesh vol 6.707→11.59 m³
    = ×1.728=det, re-voxelizes to 90110 cells, flow diverts, stable; placement survives `--apply-h 0.10`
    (11230 cells, same 2.25% of domain) AND scene save→restore→Apply). gate_G1 + voxel_flow + 93/93 unit PASS.
    ⚠ The GL rendering + drag FEEL are unverified headlessly (offscreen has no GL) — needs a windowed try.
  - **Two-colour voxel overlay (2026-07-08)**: the staircase overlay is SPLIT by cell kind so the rigid
    structure/obstacle reads apart from the erodible sand. The worker publishes a KIND mask (0 fluid,
    1 sediment, 2 rigid solid) via `SimWorker::publish_solid_kind` (sand = flow-solid AND NOT the engine's
    `SeabedMorpho::structure_host()`; a fluid viewer with no bed ⇒ every solid is rigid). `uploadVoxelOverlay`
    extracts boundary faces **PER KIND**: a cell of a given kind emits a face wherever its across-neighbour is
    a DIFFERENT kind (fluid OR the other solid kind), so each kind is its own CLOSED surface — a rigid voxel
    buried in sand keeps its full skin, and hiding one kind reveals the other's complete surface at the sand/
    structure seam (NOT a hollow shell — the earlier union test suppressed the seam faces, so a solid buried in
    sand vanished when sediment was hidden; fixed). It emits SEDIMENT faces first then SOLID faces into one VBO,
    split at `vox_sed_vertex_count_`; `paintGL` draws `[0,split)` **orange** (sediment) and `[split,end)`
    **dark-blue** (solid), each gated by its own toggle. `publish_mask` carries the kind-valued mask; the
    auto-range solid-exclusion uses the separate `disp_solid()`.
  - **Legend** (feature 3): a QPainter colorbar drawn AFTER raw GL (`begin/endNativePainting` fence,
    `drawLegendWith`). Reflects the ACTUAL slice colormap + `[vmin_,vmax_]` range with title+units
    (`|u| [m/s]`, `p [Pa]`, …); a 2nd **`Δz_b [m]` bar** (the bed diverging ramp) in seabed mode.
    Title + tick-label text is **WHITE** (reads over the dark UI / the bar).
  - **Auto-scaling bed range (2026-07-08)**: the bed diverging map + its `Δz_b` legend share ONE half-range
    `bed_half_range_` about z0. With auto-range ON (the same "Auto range" checkbox) it tracks the live peak
    `|z_b − z0|` **expand-fast / shrink-slow** (snap up so the deepest point stays saturated; ~0.8 s EMA
    contraction) with a **1 mm floor**, so sub-mm early scour/deposition already fills the ramp instead of
    vanishing against the old fixed ±25 cm band; OFF ⇒ the fixed `kBedFixedRange`. The legend labels the
    **relative** deviation `Δz_b = z_b − z0` (±half-range) because absolute `z_b ≈ z0` rounds away at the
    legend's 3 sig-figs. `updateBedSurface` computes it (host, from the per-column snapshot) so bar + surface
    always agree. `slice_viewer.cpp updateBedSurface/drawLegendWith`.
  - **Auto colour-range** (default ON, "Auto range" checkbox): the slice colormap, the legend range/ticks
    AND the arrow speed scale all follow a LIVE `[vmin,vmax]` of the displayed field so an accelerating /
    blowing-up region stays visible instead of saturating to red. Computed by a **GPU block-reduction** over
    the worker's device snapshot, FLUID cells only (`SimWorker::disp_solid()` excludes bed/structure/obstacle;
    non-finite cells skipped so a NaN can't poison the scale) — `slice_reduce_gpu`/`slice_reduce_cpu`
    (`slice_field.cu`, GPU-vs-CPU parity `tests/test_slice_field.cu`) run on the interop stream via
    `slice_gl_reduce` and copy back ONLY 3 scalars `[min,max,speed_max]` (NEVER the whole field → 60 fps +
    decoupled step/render intact), throttled to every `kRangeEvery`=5 paints. `applyAutoRange` folds them in
    **expand-fast / shrink-slow** (EMA, ~2 s contraction) so a blow-up lights up the next update but a settled
    field doesn't flicker; magnitude fields pin `vmin=0`, signed fields track both data ends. Off ⇒ the fixed
    per-field `updateRange()` bands. A throttled `[G1] auto-range …` stderr line logs the live range.
  - **Flow arrows** (features 4+5): variomigo-style animated particle-arrows. `flow_particles.{h,cpp}`
    is a Qt-free CPU tracer field (seed→advect→respawn); the viewer draws them as a **GL 4.3 INSTANCED**
    unit-arrow glyph (`glDrawArraysInstanced`, per-instance {pos,dir,speed,alpha}, billboarded to the
    camera, length∝speed). Colour: BY SPEED through the SAME 5-stop ramp as the slice when the slice is
    HIDDEN (matches the legend, reads on the dark bg), but **solid WHITE when the slice is VISIBLE**
    (`uColorMode` uniform) so the arrows contrast against the coloured field instead of blending in. The
    speed scale is the live `auto_speed_max_` (auto-range) so arrow colour reads on the |u| legend.
    The ramp is the ONE source of truth `gui/colormap.h` (`scour_colormap`, shared by `slice_field.cu`,
    the arrow fragment shader's GLSL mirror, and the legend). Master **Show arrows** toggle, **2D/3D**
    mode combo (3D = full {u,v,w} through the volume; 2D = in-plane velocity on the current slice), a
    live **Density** slider (100…8000), an **Arrow speed** slider [0,1] that scales `ArrowView.gain`
    (`gain = 1.5·mult`; 1 = default lively, 0 = frozen) so a fast current (e.g. 8 m/s) doesn't whip the
    arrows about — VISUAL ONLY (`SliceViewer::setArrowSpeedMult`, physics untouched) — and an **Arrow width**
    slider (thickness) that scales ONLY the glyph's lateral extent (the `uWidth` uniform multiplies the
    across-shaft term in `kArrowVert`, leaving the length term untouched) so lower = thinner shaft+head while
    the arrow LENGTH is kept — a declutter knob for a dense field (`setArrowWidthMult`, default 0.5 = half the
    glyph width). The tracers advect
    on the MAIN thread from a host {u,v,w,solid}
    snapshot the worker D2H-publishes at ~12 Hz ONLY when arrows are on (`SimWorker::setWantHostFlow` /
    `withFlowField`, double-buffered under `flow_mtx_`); advection runs under that lock (O(1) worker
    contention), so no big per-frame copy hits the render thread and the 60 fps hold is preserved.
  - ⚠ **DEPTH-STATE GOTCHA**: `QPainter::beginNativePainting()` resets GL to DEFAULT state (depth test
    OFF, depth mask undefined, clear colour black, blending on). So `paintGL` MUST re-establish the
    scene depth state every frame right before `glClear` (`glEnable(GL_DEPTH_TEST)` + `glDepthFunc(LESS)`
    + `glDepthMask(TRUE)` — the mask matters or the depth CLEAR is masked out). If you drop this the
    z-buffer silently dies: voxels paint on top and polygons draw in submission order.
- ⚠ MSVC (14.44) `CL.exe` occasionally crashes mid-compile on the AUTOMOC TU with exit
  `-1073741819` (0xC0000005) — it is transient; just re-run `cmake --build`. Not a code bug.
- GUI code is `src/gui/`: main.cpp (GL 4.3 `QSurfaceFormat` BEFORE `QApplication`,
  slicer main.cpp:474 pattern), camera.h (self-contained orbit/pan/zoom, QMatrix4x4),
  slice_viewer.{h,cpp} (`QOpenGLWidget`+`QOpenGLFunctions_4_3_Core`, ALL GL on the main
  thread — every entry point guards `SCOUR_ASSERT_GL_THREAD()`, a Qt **debug** assertion),
  sim_worker.{h,cpp} (steps the M2 `ChannelFluidCore` on a QThread, NO GL), sim_setup.cpp
  (builds the sim from JSON — parses directly, not via the strict core Config loader, so
  viewer-only keys U/Re/cylinder are allowed), main_window.cpp (queued-signal hand-off,
  cobod-slicer mainwindow.cpp:1691 pattern). main_stub.cpp is the `-DSCOUR_ENABLE_QT=OFF`
  fallback (Qt-less hosts).
- CUDA-GL interop is split into a Qt-free static lib `scour_gui_cuda` so nvcc never sees
  Qt host flags: slice_field.cu (pure sampler+colourmap kernel, one `__host__ __device__`
  evaluator ⇒ the GPU-vs-CPU parity test `tests/test_slice_field.cu` is exact, no GL/Qt
  needed) + slice_gl.cu (the ONLY GL-header TU: `cudaGraphicsGLRegisterBuffer`, zero-copy
  writes into the registered colour VBO, on its OWN non-blocking CUDA stream).
- ⚠ **Rendering is decoupled from stepping** (this is what makes it 60 fps): `step()` runs
  on the worker with NO lock; after each step the worker copies {u,v,w,p} into device
  snapshot buffers under `display_mutex()`; the main-thread slice kernel reads only those
  snapshots, on the interop's own stream. Measured: **~60 fps hold** on both the bundled
  reduced config (200×100×10, sim ~115 steps/s) and the full 10×10×5 m @ h=5 cm (4M cells,
  sim ~20 steps/s ≈ 50 ms/step) — the viewer stays at the 60 Hz repaint cap regardless of
  step cost. If you ever route the render kernel onto the default stream or lock across
  `step()`, fps collapses to the step rate (~14 fps) — keep them decoupled.
- M1 fluid core lives in `src/core/fluid/` (mac_grid.h, advect.cu, turbulence.cu,
  project.cu, mgpcg.cu, stam_fluid_core.cu, cavity.cpp). Fields are **double** at M1
  (2M cells at 128³ is tiny; double makes both the <5% Ghia match and the 1e-6
  divergence gate robust). Switch velocity to float for the 32M-cell hero grids later.
- GPU-vs-CPU kernel parity tests: `tests/test_fluid_kernels.cu` (rel. max-norm 1e-5);
  MGPCG convergence/grid-independence: `tests/test_mgpcg.cu`.
STEP model loading in scour-gui (OpenCascade; increment 1 = LOAD + DISPLAY ONLY, no
voxelization yet). File menu → "Open STEP…" (and "Close model"), or the CLI flag
`--load-step <path>` (loads at startup; smoke-testable headlessly with `--offscreen`).
Loads + triangulates the STEP into a lit GL mesh drawn in the domain (centred on x/y,
sat on the bed z=0) over the live sim. Verified (2026-07-08): cube.step → 12 triangles,
bbox 0.01 m cube; the V000 unit (`Experiments/…/ScourProtection V000_001.stp`) → 70
triangles, bbox 7.2×7.2×2.0 m; both at ~60 fps, exit 0, double-clickable with NO OCC/Qt on
PATH.
```
build/Release/scour-gui.exe --load-step "C:/CODE/cobod-slicer/tests/inputs/cube.step" --offscreen --autoclose-ms 8000
```
- **`scour_geometry`** static lib (`src/core/geometry/step_import.{h,cpp}`) is the ONLY target
  that links OpenCascade — it mirrors how `scour_gui_cuda` isolates the GL/CUDA interop, so
  `libscour` and every physics gate test stay OCC-FREE. `step_import.h` is OCC-free (a plain
  `TriMesh{positions,normals,indices,bbox}`); all OCC includes live in the .cpp. scour-gui
  links it; the physics core does NOT depend on it.
- OCC linkage (OpenCascade 8.0, `C:/OpenCASCADE-8.0/build2`): inc = `${OCC}/inc`, libs =
  `${OCC}/win64/vc14/lib`, DLLs = `${OCC}/win64/vc14/bin`. Minimal STEP+mesh toolkits linked:
  `TKernel TKMath TKG2d TKG3d TKGeomBase TKGeomAlgo TKBRep TKTopAlgo TKMesh TKXSBase TKDESTEP`
  (PRIVATE on the static lib, but CMake still propagates them for the final link). Compile with
  `OCCT_NO_DEBUG OCCT_NO_DEPRECATED` (+ the 8.0 `NCollectionAliases` deprecated-typedef dir).
- ⚠ **DLL deployment (double-clickable):** `windeployqt` handles ONLY Qt, so a CMake POST_BUILD
  (`scour_deploy_occ_dlls()`) copies the OCC runtime **DLL closure** next to the exe: the 23-DLL
  transitive set of the link toolkits (dumpbin `/dependents` — `TKDESTEP.dll` pulls in
  `TKDE/TKXCAF/TKLCAF/TKShHealing → TKV3d/TKService`) **plus 6 3rdparty DLLs**
  (`freetype brotlicommon brotlidec bz2 libpng16 zlib1` — `freetype.dll` is a LOAD-TIME dep of
  `TKService.dll`, so a missing one = `0xC0000135` at startup with no output). A qoffscreen
  platform plugin is also deployed into `platforms/` for the headless `--offscreen` path. Verify
  clean-machine startup by launching with PATH = System32 only (no OCC/Qt/JDK).
- ⚠ STEP → mesh **traps** (in `step_import.cpp`): (a) OCC emits **millimetres** → every node ×0.001
  to METRES; (b) **winding** — swap two triangle indices when `face.Orientation()==TopAbs_REVERSED`
  for OUTWARD normals (this is the OPPOSITE branch from cobod-slicer, whose meshes are globally
  inverted); (c) per-vertex normals = area-weighted average (sum of un-normalised face crosses).
- The STEP loader has its OWN OCC-linked test exe **`scour_geometry_tests`** (`tests/test_step_import.cpp`,
  label `unit`, DISCOVERY_MODE PRE_TEST so discovery runs after the DLL deploy) loading the committed
  `tests/inputs/cube.step` — the physics `scour_tests` stays OCC-free. `74/74 unit` PASS incl. it.
STEP → fluid VOXELIZATION (increment 2 — the STEP model becomes a live solid obstacle the flow
diverts around). Verified (2026-07-08): cube unit gate 0.000% vol err; V000 in the fitting 10×10×5 m
domain → **watertight, 53668 solid cells, vol err 0.021% (<2%), 0 thin cells**; live flow diverts
(max|u| 1.78 m/s > U, mean|u| in solid = 0) at 59 fps on 4M cells.
```
build/Release/scour_voxel_flow_gate.exe               # headless flow-diversion gate (PASS)
ctest --test-dir build -C Release -L voxel_flow       # same, as ctest
# GUI: File→Open STEP… — the loaded model IS the obstacle (auto-voxelized on load; no toggle); or headless:
build/Release/scour-gui.exe --config configs/g1_viewer_full.json --load-step <model.stp> --autoclose-ms 9000
```
- **Voxelizer** (`src/core/geometry/voxelize.{h,cpp}`, in **libscour**, OCC-FREE — needs only the
  mesh): `voxelize_mesh(TriMesh, MacGrid, ModelPlacement, …)` → the `ChannelBC` solid mask
  (`std::vector<unsigned char>`, 1=solid, `g.pidx`, same format as `build_cylinder_mask`). Watertight
  **ray-parity** (axis-aligned +z ray, even-odd triangle-crossing count) as the per-point inside test,
  under a **K³ supersampled MAJORITY fill** (K=3 default; a cell is solid iff inside-fraction > 0.5).
  The per-cell inside-fraction is **retained** (`out_fraction`) — the cut-cell volume fraction for
  later; NO cut-cell physics yet. A **thin-wall safeguard** keeps any isolated sub-cell partial cell
  (0<frac≤0.5 with no solid face-neighbour) solid so a < h barrier can't leak (over-thicken beats
  leak), logging the count. **BRUTE FORCE by design — no BVH/k-d tree** (input is watertight low-poly,
  ≈70 tris; voxelize runs ONCE at load, obstacle is rigid); only accel is a per-triangle xy-AABB cull.
  ⚠ TIE-BREAK: an asymmetric sub-cell nudge (ox≠oy) keeps sample points off axis-aligned edges/
  diagonals so a coplanar-triangle shared edge is owned by exactly one triangle (else parity flips).
- **`TriMesh` moved to OCC-free `core/geometry/tri_mesh.h`** (libscour include path) so the voxelizer
  consumes it with NO scour_geometry/OCC dependency; `step_import.h` includes it. **`ModelPlacement`
  (`core/geometry/model_placement.h`, `place_model_on_bed`)** is the ONE shared display+voxelize
  transform (centre x/y, sit min-z on the bed) — the mask lands exactly where the mesh is drawn.
- **Live injection (cross-thread):** loading a STEP (File→Open STEP… / CLI `--load-step`, `--noslip`
  optional) **auto-voxelizes** the model on the sim grid and injects it as the obstacle — the loaded
  model IS the obstacle (there is no "Model as obstacle" toggle / Model menu; `--voxelize` is a kept
  no-op). The injected mask is the model voxels ALONE — it **REPLACES the config obstacle (the default
  cylinder)**, which returns via "Close model" (`setModelAsObstacle(false)` restores `recipe_.base_solid`).
  In a seabed scenario the bed+structure own the mask, so loading a STEP there only displays the mesh.
  Injection **rebuilds the ChannelFluidCore ON THE WORKER THREAD** via `SimWorker::requestRebuild` (mask queued
  under a mutex + atomic flag; the worker swaps the core at the top of its loop, `make_core` from an
  immutable `SimRecipe` snapshot — the core is NEVER constructed from the main/GL thread). Default
  surface **`SOLID_FREESLIP`** (RESEARCH §3 — the printed structure's BL is unresolvable at these
  voxels). The decoupled render/step + `display_mutex` snapshot pattern is untouched (59 fps held).
- **Voxel overlay (LIVE flow mask):** the "Show voxels" overlay draws the EXPOSED voxel-surface faces
  (staircase) of the **current flow solid mask** — structure ∪ sand bed ∪ config obstacle ∪ ANY
  runtime-marked solid — so every cell the fluid treats as solid is visible (a mis-placed solid voxel
  shows up). The worker snapshots the host mask + bumps `SimWorker::maskGeneration()` ONLY when the mask
  changes (load, obstacle rebuild, each bed `update_solid`) — no per-frame mask D2H; the viewer polls the
  generation each paint and re-extracts the exposed faces (`SliceViewer::setVoxelOverlay`/`copyMask`) only
  then, all GL on the main thread. Drawn in orange over the smooth mesh. (The old STEP-only overlay push in
  `setModelAsObstacle` was dropped — the merged mask the rebuild publishes now supersedes it.)
- **XYZ axis triad (`SliceViewer::drawAxes`/`drawAxesLabels`, "Show axes" toggle, default ON):** three
  coloured axes from the **WORLD ORIGIN (0,0,0)** = the x=0/y=0/z=0 domain corner (X=red, Y=green, Z=blue),
  drawn depth-tested in the scene with **metre tick marks + tick/letter labels** (nice ~2 m spacing for a
  10 m domain) so the user can eyeball world coordinates. NB `place_model_on_bed` centres the model in x/y
  at z=0, so the origin is the domain CORNER, not under the model. Plus a small **camera-aligned orientation
  gizmo** (bottom-left, rotation-only orthographic transform, always on top) so the axis directions read from
  any view even when the origin is off-screen. Letters + tick numbers go through the same QPainter-after-GL
  path as the legend (`projectPoint` world→screen); lines via the flat slice shader (`uFlat`, per-axis colour).
- Tests: `tests/test_voxelize.cpp` (in `scour_tests`, OCC-free procedural-cube volume <2% + inside/
  outside + ×8 scaling gates) and `scour_voxel_flow_gate` (label `voxel_flow`: box → small channel,
  step 400×, assert mean|u| in solid ≈0 and peak fluid speed > U). ⚠ A model LARGER than the domain is
  clipped (near-total blockage → the flow can blow up); this is a domain mismatch, not a voxelizer bug
  — size the domain to the model (V000 needs the full 10×10×5, not the reduced g1_viewer 10×5×0.5).
ERODIBLE-SEABED demo scenario (increment 3 — the FIRST full live morphodynamic loop coupled to the
live flow: fluid → bed-shear wall model τ_b → suspended sediment → bed Exner update → avalanche, on a
1 m sand bed with the voxelized V000 structure seated in it, plus a height-coloured bed surface so you
can WATCH a scour hole develop). ⚠ **PRE-CALIBRATION, QUALITATIVE — NOT a measurement, NOT calibrated**
(α/Cs/HSV calibration is the SUPERVISED M6 gate). Constants NOMINAL; the config + a startup log line say
so. Verified (2026-07-08): flat-bed sanity **mass err 3.5e-16, stable**; procedural-box + V000 both
**stable, mass-conserving (V000: 1e-7 over 1656 steps), scour signature develops** (z_b lowers where τ_b
high, rises where sheltered). V000 offscreen: 53668 structure cells seated at z=1.0 m, 4M cells @ ~15
steps/s, z_b min 1.000→0.956 / max →1.010, max_ustar bounded (no blow-up).
```
# GUI (windowed; offscreen has no GL so it steps but doesn't render):
build/Release/scour-gui.exe --scenario configs/seabed_v000.json   # or pass the file positionally
# File menu → "Load seabed scenario…" loads one at runtime (rebuilds core+engine on the worker).
ctest --test-dir build -C Release -L seabed_flat --output-on-failure   # headless flat+box gate
build/Release/scour_seabed_gate.exe                                    # same, directly
```
- **Morpho loop lives in the WORKER** (`src/gui/sim_worker.cpp`): after `core.step()` (fluid), if the
  bed is active and past `spinup_steps`, `SeabedMorpho::step(u,v,w,ν_t,dt,&new_solid)` runs the whole
  loop and, when the bed has moved > 0.5·h since the last mask build, returns the updated solid mask →
  `core.update_solid()` re-masks the flow IN PLACE (flow preserved, unlike a factory rebuild). Decoupled
  render/step + `display_mutex` snapshot pattern intact; the worker also publishes a per-column z_b
  snapshot for the bed viz. NO GL off the main thread.
- **Engine** (`src/core/sediment/seabed_engine.{h,cpp}`, Qt-free libscour) drives the gate-verified
  kernels UNCHANGED (`suspended.*`, `bedstate.*`, `avalanche.*`, `bedshear_*` smooth/EMA) + THREE thin
  bridge kernels (`src/core/sediment/seabed_morpho.{h,cu}`, GPU-vs-CPU parity `tests/test_seabed.cu`):
  (1) **seabed_bedshear** — τ_b probed at z_p ABOVE the CURRENT bed top z_b=G/c_pack (the gate wall model
  assumes the bed at k=0); (2) **seabed_confine_c** — lift any suspended c below the bed / inside the
  structure back into the first fluid cell (mass-conserving; the gate suspended kernels have no interior
  bed); (3) **seabed_bed_post** — clamp G∈[0,Gmax] + pin the rigid structure-footprint columns. Bed
  coupling to c is the gate `morpho_exner` (deposition+erosion exchange at k_bed). Suspended transport =
  effective-w_s + `suspended_advect_cons` (mass-exact) + ν_t/σ_s diffusion (sub-stepped for stability).
- **Cell kinds reconciled**: STRUCTURE (rigid — its footprint columns are frozen so they neither erode
  nor accrete; buried sand is also naturally shear-protected) OR SAND (erodible f_pack from G, solid to
  the flow where z<z_b) OR FLUID. The flow solid mask = sand ∪ structure, rebuilt on the worker as z_b
  moves. Structure is voxelized on the bed (`place_model_on_bed` + `tz += sand_depth`); the viewer draws
  it at that same translate (`SliceViewer::setMeshTranslate`).
- **Bed viz**: `SliceViewer::setBedSurface` renders a per-column height field z_b(x,y) coloured by
  elevation (blue = scoured, warm = deposited), refreshed each paint from the worker's z_b snapshot;
  reuses the slice shader (per-vertex colour), ALL GL main-thread. Kept alongside the slice + structure
  mesh. The plain `g1_viewer*.json` configs still run FLUID-ONLY (bed off ⇒ G1 unchanged, `gate_G1` PASS).
- ⚠ **STABILITY FIX (root cause found + fixed):** a solid bed that reaches the INLET plane made the
  Dirichlet inlet drive the buried sand cells at i=0 (unremovable inflow) → a corner anomaly that
  coupled to any obstacle wake and blew the flow up (worse with a deeper reservoir). Fix:
  `ChannelFluidCore::set_bed_inlet_mask(true)` zeros the i=0 inlet u-faces into solid cells (guarded, OFF
  by default ⇒ M2/M3 byte-identical — their obstacles never touch the inlet). The FLOW uses an artificial
  viscosity `nu_fluid` (≈U·L/Re, Re~150 — the g1_viewer_full recipe) for coarse-grid stability; the
  SEDIMENT uses molecular ν (`nu_sed`). A very bluff obstacle still over-constricts the coarse M2 core
  (that regime is M6's calibrated turbulent flow); V000 (few % blockage in 10×10) is stable. MORFAC ≤10.
- Increment-4 = **M6** (SUPERVISED pile-scour calibration of THIS exact live pipeline): the τ_b is already
  wired live from the wall model; M6 swaps the uniform inlet for the precursor turbulent channel, freezes
  the bed for the u*=0.020 pre-check, then calibrates α (±30%) / Cs (0.10–0.12) / optional HSV against
  Roulund S/D — never an agent's call. Re-run M1–M5 + gate_G1 after any constant is frozen.
Ninja fallback (must be inside the VS2022 Professional x64 env so nvcc finds cl 14.44):
```
cmd /c "\"C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat\" && cmake -B build -S . -G Ninja && cmake --build build"
```
- Toolchain: MSVC 14.44 (VS2022 **Professional**), CUDA 13.1 (`CMAKE_CUDA_ARCHITECTURES` default `75;86;89`),
  CMake 4.1, Ninja (VS-bundled — not on PATH; use the VS2022 x64 dev prompt).
- ⚠ VS 18 Community (MSVC 14.51) is also installed and CUDA 13.1 **rejects it**. Build from
  the VS 2022 Professional dev prompt or pin `CMAKE_CUDA_HOST_COMPILER` to the 14.44 cl.exe
  (full path in PLAN §1). Never pass `-allow-unsupported-compiler`.
- Qt 6.11.1 auto-detected via glob `C:/Qt/6.*/msvc2022_64` (don't hard-code versions).
- OpenCascade 8.0 at `C:/OpenCASCADE-8.0/build2` (STEP import only; minimal TK list in PLAN §3).

## M6 + M7 + M8 — HARNESSES BUILT 2026-07-09 (calibration SUPERVISED — see M6_REPORT.md)

Status 2026-07-09: the **M6/M7/M8 harnesses are built + committed** (`gate_M6` Roulund pile scour,
`gate_M7` Du deposition, `gate_M8` shape ranking = PASS), plus the sub-grid porous momentum sink
(`channel_porous`), campaign machinery (`core/campaign/run_protocol.h`), and the HSV shear-multiplier
knob. **92/92 unit.** Build + run (each long run: background + poll):
```
build/Release/scour_m6_gate.exe configs/m6_roulund.json [--quick] [--hsv X] [--noslip] [--erosion E] [--morfac M] [--csv out.csv]
build/Release/scour_m7_gate.exe configs/m7_du.json      [--quick] [--U u] [--only control|car|har] [--csv out.csv]
build/Release/scour_m8_gate.exe configs/m8_ranking.json [--quick] [--shapes plate,ring,cup] [--csv out.csv]
ctest --test-dir build -C Release -L gate_M8 --output-on-failure   # gate_M8 + gate_M7 PASS; gate_M6 FAIL (absolute, uncalibrated)
```
- **M6 is calibration-BLOCKED, not code-blocked (READ `M6_REPORT.md`).** Harness verified: frozen-bed
  u\* pre-check PASS (0.0196 vs 0.020), stable on molecular ν, mass-conserving, correct qualitative
  scour — but uncalibrated **S/D ≈ 0.13 ≈ 10× under Roulund's 1.25** (the RESEARCH §11 #1 SL-solver
  HSV under-resolution). Sweeps: 5× erosion → +2 %; no-slip pile → +40 %; HSV multiplier saturates
  (×4→0.27, ×8→0.33, ×12→0.35); **D/20 grid refinement → NO gain (0.185→0.189 at 23 M cells) ⇒ the SL
  advection scheme, not resolution, is the limit.** The sanctioned ±30 % sediment knobs CANNOT close
  a 10× gap. **gate_M7 (Du) PASSES its COMPARATIVE gate** (porous units reduce scour 6–10 % + trap
  sand, stable, mass-exact) with the SAME ~9× absolute deficit — the comparative behaviour is sound.
  **M6_REPORT §5: accept comparative-ranking-only (recommended) or LBM plan-B; nothing frozen;
  `configs/calibrated.json` = PROVISIONAL.** ⚠ CLOSED-mode Exner erosion uses Winterwerp
  `MorphoParams.erosion_coeff` (default 0.018, exposed), NOT the van Rijn α (α drives the M4 suspended
  path + open-sea boundary) — M6_REPORT §6.
- ⚠ Grid-count doc error: M6 at D/10 is **2.9 M cells**, not "≈19 M" (see below / M6_REPORT §2).
- **SUPERVISED**: an agent may build the harness + RUN sims, but the **calibration decisions (α, Cs,
  erosion_coeff, HSV multiplier, grid) are the orchestrator's + MH's — never an agent's**. After any
  freeze, **re-run M1–M5** (all additions are byte-identical by default). M7 must pass before any
  shape ranking is trustworthy; gate_M8 already ranks correctly on the COMPARATIVE metric.

### Original M6/M7 spec (below) — the physics targets remain the authority

### M6 — Pile-scour calibration (`gate_M6` / V7) — the big calibration gate
Roulund et al. (2005) live-bed pile scour. **Engelund–Fredsøe bedload (`MorphoParams.bedload_formula=1`)
+ suspended load BOTH ON** (disabling suspended HALVES the scour — RESEARCH §10 V7).
- **Setup (PLAN M6 / RESEARCH V7):** voxel h = 1 cm (= D/10; REEF3D: D/10 → <4% err, D/5 unstable);
  domain **3.0 × 1.6 × 0.4 m**; pile **D = 0.1 m** centred 1.0 m from inlet, ≥ 7D lateral clearance;
  ≈ **19M cells**. **V = 0.46 m/s, d50 = 0.26 mm, live-bed V/Vcr = 1.25**, h_dom = 0.4 m.
- **Inflow = precursor periodic channel** (this is the real 3-D turbulence M3's 1-D mixing-length
  leg deferred; reuse `precursor.{h,cpp}` record/replay built in M3). **MANDATORY pre-check before
  ANY morphology: freeze the bed, spin up 3–5 flow-throughs, verify recovered u\* = 0.020 m/s ± 10%**
  (research/14). Do NOT enable the Exner update until this passes.
- **Run protocol (RESEARCH §7 / §10):** MORFAC **M is the ONLY time accelerator** (M ≤ 10 steady
  current); the wall-model EMA is de-noising only. Run **≥ 400 s morphological time** (T ≈ 130 s,
  research/11); fit BOTH S(t) = S_eq(1 − e^(−t/T)) AND the Welzel hyperbolic a·(1 − 1/(1+bt));
  never fit < 0.3T, expect 10–20% S_eq error at 0.5T, extend if the two fits disagree > 15%. Verify
  one M=1 vs M=target run agree before trusting accelerated runs.
- **Calibrate at most THREE knobs, within these bands ONLY:** van Rijn pickup **α = 0.00033 ± 30%**
  (RESEARCH §6.2); Smagorinsky **Cs = 0.10–0.12** (RESEARCH §3.4); optional **HSV (horseshoe-vortex)
  shear multiplier**. SL solvers under-resolve the HSV and underpredict its bed shear ~30%
  (RESEARCH §11), so expect to calibrate, then **FREEZE the chosen constants in
  `configs/calibrated.json`**.
- **Gate (V7):** upstream **S/D = 1.25 ± 15%**; downstream ± 30%; fitted timescale within **factor 2**
  of measured; a **sensitivity run at 1.5× resolution changes S/D < 10%**.
- E–F near-threshold p-factor (the Roulund form, already in `sed_physics.h`):
  p = [1 + ((π/6)·0.51/(θ − θ_cr))⁴]^(−1/4). τ_b (grain-skin, per-column vector) is now wired LIVE
  from `fluid/bedshear.cu` (it was synthesised in the M5 gate). Notes: research/05 (HSV/S/D/
  timescales), /11 (benchmarks + T), /14 (inflow/BCs), /16 (MORFAC/run-length).

### M7 — Deposition validation "M-DEP" (`gate_M7` / V8) — MUST PASS BEFORE ANY SHAPE RANKING
Du et al. 2025 perforated-unit case + a current-only trench variant vs van Rijn trapping efficiency
(the wave-stirred van Rijn 1986 trench moves to M10 — 0.18 m/s current alone is sub-threshold for
0.1 mm sand).
- **Setup (PLAN M7 / RESEARCH V8):** 4–5 mm voxels (lab scale); **d50 = 0.235 mm, U = 0.2–0.3 m/s**;
  measured scour-reduction tables in research/15. **Build the C-AR / H-AR units procedurally as voxel
  masks** (no STEP import — that is M8).
- **Thin-screen porous model (moved here from M8):** Δp = ½·ρ·k·u_face², **k = 1/β² − 1**
  (RESEARCH §8), applied to the perforated faces. At 4 mm voxels the 8.33 mm holes are ~2 cells
  across — **below the ≥ 8-cell resolution floor (research/13), so model them SUB-GRID; do NOT
  attempt to resolve the hole jets.**
- **Gate (V8):** trapped volume ± 30%; deposition depth ± 20–30%; **BSS ≥ 0.3** (report if ≥ 0.6;
  target 0.6); measured **C-AR / H-AR ranking reproduced**; van Rijn trapping efficiency
  **e_s = 1 − exp(−A_vr·L·d/h1²) within ± 30%**; qualitative deposition behind/inside the porous
  array. Notes: research/07 (permeable structures/trapping), /15 (Du 2025 tables), /13 (res floor).

## M9 — TIDAL REVERSAL + ranking campaign (`gate_M9`) — BUILT + PASS 2026-07-09

Face-swap tidal reversal (RESEARCH §8) on the live morphodynamic pipeline + the first reversing
shape-ranking campaign. **gate_M9 PASS; 93/93 unit (M2/M3 byte-identical — reversal is `flow_sign=+1`
by default).** Long run: background + poll.
```
build/Release/scour_m9_gate.exe configs/m9_tidal.json [--quick] [--shapes plate,ring,cup] [--morfac M] [--half-cycles N] [--noreverse] [--U u] [--csv out.csv]
ctest --test-dir build -C Release -L gate_M9 --output-on-failure   # ~5–8 min at 0.30 M cells (bg+poll)
python tools/gen_report.py out.csv -o m9_report.html               # ranked HTML report (reversal-stability strip)
```
- **Face-swap reversal is a fluid-core BC change** (`ChannelBC::flow_sign` ±1, `ChannelFluidCore::
  set_flow_direction`): +1 = M2/M3 inlet@xmin / Orlanski-outlet@xmax (**default ⇒ byte-identical**,
  gate-verified); −1 SWAPS the faces (Dirichlet inlet drives −x on xmax, outlet on xmin). `U_inlet`/
  `ustar` stay MAGNITUDES; the sign carries direction. Generalised across the inlet (`k_inlet`), the
  Orlanski outlet (`k_orlanski`, outflow-only clamp mirrors sign), the inlet/outlet flux-rescale, the
  **pressure Dirichlet-p=0 pin** (`ChannelMgpcg::set_outlet_dir` → the `dir_xmax` int is now a signed
  face code: >0 xmax, <0 xmin, threaded through `nb_stencil`/`k_subgrad` + CPU twins), and the
  bed-inlet mask. The driver flips `flow_sign` at slack (|U|→0) so the swap is smooth. `channel_bc.h`,
  `channel_ops.cu`, `channel_pressure.{h,cu}`, `channel_core.{h,cu}`.
- **The SEDIMENT side needs NO change**: the bed-shear wall model sets τ direction from the signed
  velocity (`taux = τ·uc/Up`), so bedload/pickup/deposition follow the reversed flow automatically; the
  suspended conservative advection seals all six faces + the open-sea boundary is per-face upwind
  (already bidirectional).
- **Tidal schedule** (gate): start at +U_max; each half-cycle holds a plateau then cosine-ramps to the
  opposite sign over T_ramp (the compressed slack, ≥ 3–6 flow-throughs, RESEARCH §8), slack at the ramp
  midpoint. 4 half-cycles ⇒ **4 slack transitions**. **MORFAC ≤ 5** (RESEARCH §7 reversing; config M=4).
- **Gate criteria (all PASS):** (1) total sand (bed+suspended) conserved **< 0.1 %** — with the CLOSED
  suspended boundary the conservative advection seals the domain, so it is conserved to **machine
  precision (~1e-14)** through the reversals; the budget CREDITS the structure-footprint pin discard
  (sand cannot accrete under a rigid slab — a physical sink, `SeabedMorpho::structure_discard()`, like
  ±boundary fluxes) so the residual is the true transport test. (2) **no monotonic KE growth** across
  the 4 slacks — measured on the MATCHED-PHASE KE (sampled at each plateau peak |U|=U_max, where the
  bulk-flow KE baseline is invariant): ×1.00–1.01, flat. ⚠ Do NOT gate on the slack-KE (|U|→0): it is
  dominated by the still-developing scour recirculation (the plate's rises ×1.7 over 2 cycles) and
  would false-fail a stable flow — a bed-development transient, not energy pumping (maxu stays 0.96,
  matched-phase KE flat). (3) MORFAC limiter **clip-fraction < 1 %** (`SeabedMorpho::clip_fraction()`
  wires the `morpho_exner` `clip_count`): 0.000 % at M=4.
- **First reversing ranked report (≥3 shapes):** plate / ring / cup on the identical reversing tide,
  ranked by matched-time trapped volume. **ring (2.70e-3 m³) ≫ plate (7.3e-4) > cup (3.8e-4)** — under
  a reversing tide the SYMMETRIC enclosure (ring) shelters its interior from BOTH flow directions and
  keeps the most sand; the square cup is under-resolved on the compact patch. Comparative/qualitative
  (RESEARCH §11; absolute magnitudes inherit the gate_M6 HSV deficit). The gate asserts the enclosing
  PRINCIPLE beats the plate + a distinguishable pair (not a specific enclosure — that would over-fit).
- Uses a CLOSED sediment boundary (finite LOCAL sand budget): the ranking is "which shape gathers/keeps
  sand as the tide flips", which sidesteps the OPEN-mode storm-ACCUMULATION artifact (HANDOVER: higher
  U ⇒ more open-sea supply ⇒ retention ≫ 100 %). `--noreverse` gives a steady control.
- `SeabedMorpho` gained two pure DIAGNOSTICS (no physics feedback ⇒ existing gates unaffected):
  `clip_fraction()`/`reset_clip_stats()` (MORFAC hit-rate) and `structure_discard()` (sand removed by
  the bed clamp+pin, for the open-budget audit). `tests/gate_m9_main.cpp`, `configs/m9_tidal.json`.
- **Interactive GUI driver (scour-gui, GUI-only — gate_G1 PASS, no physics change):** a **Tidal reversal
  (live)** section in the Simulation dock (enable + Peak U_max / Plateau / Slack-ramp spins). When on, the
  worker drives `U_d(t)` each step via `core.set_inlet_speed(|U_d|)` + `core.set_flow_direction(sign)`
  (+ `engine.set_inlet_speed`), face-swapping at slack — so you can watch a scour hole reorganize as the
  tide flips. It takes over "Input speed U" while active (greyed) and restores steady +x flow on disable;
  the status bar shows the live phase ("tide +0.62 m/s (flood/ebb/slack)"). Worker-thread apply (same
  atomic-dirty pattern as the other live controls). Headless smoke `--tidal Umax,plateau,ramp` (verified
  stable on the cylinder viewer AND a converted seabed run). `sim_worker.{h,cpp}` (`setTidalReversal`/
  `tide_u_signed`), `main_window.{h,cpp}` (`syncTidal`, tidal group), `main.cpp`.

## Hard rules
- A milestone is DONE only when its `gate_M<k>` ctest passes. Never weaken a gate to pass it;
  if a gate looks wrong, say so and stop.
- Every CUDA kernel has a CPU reference implementation and a GPU-vs-CPU test.
- Physics constants: RESEARCH.md only. Missing constant → flag it, don't invent it.
- Sediment code consumes only `{u, v, w, ν_t, τ_b}` via the `FluidCore` interface (LBM
  plan-B must stay drop-in; the Stam core's τ_b comes from `fluid/bedshear.cu`).
- Threshold θ_cr gates erosion/bedload only — **deposition is never thresholded** (this
  asymmetry is the product's core mechanism; breaking it silently invalidates all results).
- OCC STEP geometry arrives in millimeters → multiply by 0.001. Voxelizer must verify
  enclosed volume vs mesh volume (< 2% mismatch).
- GUI: all OpenGL on the main thread; sim on worker thread; hand-off via queued signals
  (pattern: cobod-slicer `mainwindow.cpp:1691`).
- Don't copy code from FluidX3D (non-commercial license). GPL sources: reference, don't paste.

## Style & conventions
- `.clang-format` copied from cobod-slicer (Allman + IndentBraces, tabs, ColumnLimit 0).
- Conventional Commits (`feat:`, `fix:`, `test:`, `docs:`, `refactor:`).
- One class per file; snake_case filenames; namespaces `scour::core`, `scour::gui`.
- Units: SI everywhere in code (m, s, kg, Pa); document any field's units at declaration.

## Reference projects on this machine
- `C:\CODE\cobod-slicer` — proven Qt6+OCC patterns: STEP→mesh extraction
  (`src/app/widgets/render/mesh.cpp:494`), GL context setup (`src/main.cpp:474`), camera,
  BufferArena/fence recycling (`opengl.h/.cpp`), threading hand-off (`mainwindow.cpp:1691`).
