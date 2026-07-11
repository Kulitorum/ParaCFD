# HANDOVER.md — Orchestrator handover (refreshed 2026-07-09, GUI model-placement-gizmo session)

> **2026-07-09 session summary (LATEST — feature tip `ebe2d2c`, all pushed, tree clean):** a **3D
> transform gizmo for placing loaded models/structures** + an **arrow-width declutter slider** (GUI-only;
> NO physics change — `gate_G1` + **93/93 unit** + voxel_flow PASS, all M-gates byte-identical). MH asked
> for a "proper 3D-modelling gizmo" to move/rotate/scale a loaded STEP model, then voxelize it *where he
> placed it* on the next Apply. Six commits (`fcf4b10`→`ebe2d2c`):
> - **Unified transform gizmo** (`fcf4b10` first, then generalised in `01b6fb1`): an **Enable manipulator**
>   dock checkbox draws ALL handles at once — 3 translate arrows, 3 rotate rings, 3 scale cubes (X red /
>   Y green / Z blue) + a grey uniform-scale centre — and **whichever handle you grab picks the operation
>   dynamically** (no mode switching). Picking = screen-distance to the projected handle, tiered so the
>   compact scale/uniform cubes win over the shaft/rings; a left-drag off any handle still orbits. Translate
>   projects the mouse delta onto the screen-axis, Rotate intersects a world pick-ray with the rotation plane
>   (sign-correct), Scale is exp(screen-delta)/radial. All GL is main-thread flat lines reusing `prog_`
>   (depth-test OFF, always on top). `ebe2d2c` shrank the whole gizmo to **1/3 size** and made the **scale
>   handles 1.5× larger + moved OUTSIDE the ring** (kScalePos 1.25·g) for an easy pick — per MH.
> - **The core change is `ModelPlacement` upgraded to a full AFFINE** `world(v)=M·v+t` (M = row-major 3×3
>   rotation·scale; default identity ⇒ `place_model_on_bed` + every existing path byte-identical, verified by
>   the voxelize unit gates). `voxelize.cpp` applies the affine + place-corrects the mesh-volume gate by
>   `|det M|`. `SliceViewer` owns the TRS (quaternion rot about the bbox centre) and exposes it as the ONE
>   `ModelPlacement` the display and voxelizer share ⇒ the mask lands exactly under the drawn mesh.
> - **Voxelize-where-placed + survives Apply for ALL paths** (`fcf4b10`, `71fe7fc`, `848d648`): the fluid
>   obstacle re-voxelizes at the gizmo placement on Apply (placement captured before teardown, re-applied
>   after `setMesh`); a **seabed structure is now gizmo-EDITABLE too** (seeded via `setMesh`+`setModelPlacement`,
>   NOT a fixed override) and Apply re-voxelizes it into the sand where placed — `build_seabed_sim` gained a
>   `structure_place` override, threaded from the viewer (captured independent of `model_mesh_`). `71fe7fc`
>   un-ghosted the group when sand is in the scene (`updateGizmoUi` gates on `hasModelPlacement()`, not
>   `seabed_`); `848d648` fixed "gizmo dead + transform reset after Apply" (the override + config re-voxelize).
>   Verified headless: fluid place (scale 1.2 → 11230 cells) + seabed structure (+0.5 z → re-voxelized at
>   z=1.5, was resetting to 1.0) + converted, all survive Apply; scene save/restore round-trips the affine
>   `m[9]` (scene_io extended, back-compatible). ⚠ Domain-CHANGING Apply keeps world-fixed placement (Reset
>   re-centres); resolution-only Apply unmoved is byte-identical to the old centre-on-bed path.
> - **Arrow width (thickness) slider** (`28cd467`): a Visualization-dock slider that scales ONLY the arrow
>   glyph's lateral extent (new `uWidth` uniform on the across-shaft term; length untouched) so lower =
>   thinner/less cluttered — MH's "too dense" fix; default 0.5 (half width). `SliceViewer::setArrowWidthMult`.
> - ⚠ **OPEN — the interactive GL draw + drag FEEL are UNVERIFIED** (offscreen has no GL, so the gizmo never
>   actually painted or took a real mouse drag in my smokes; the placement→voxelize *plumbing* is verified).
>   Needs a windowed pass; the handle radii / pick tolerances / drag sensitivities are single-constant tweaks
>   in `slice_viewer.cpp` (`gizmoSize`, `kTransLen`/`kRingRad`/`kScalePos`/`kCubeHalf`, tier-1/2 pick px).
>   `model_placement.h`, `voxelize.cpp`, `slice_viewer.{h,cpp}`, `main_window.{h,cpp}`, `sim_setup.{h,cpp}`,
>   `scene_io.cpp`, `main.cpp`, `CLAUDE.md`. Nothing else open from this session.

> **2026-07-09 session summary (feature tip `e3002dc`, all pushed, tree clean):** a **GUI
> bed-surface diagnostic** (GUI + one write-only kernel output; NO physics change — `gate_G1` +
> **93/93 unit** + `seabed_flat` PASS, all M-gates byte-identical). Grew out of MH asking why the
> **Concentration** slice legend maxes near ~0.9 m³/m³: that's a near-bed one-cell artifact (settling
> pile-up + the `confine_c` lift stacking a column into one cell; > c_pack=0.64 ⇒ not physical), which
> the auto-range latches onto. MH wanted an honest "where is sand arriving vs leaving" view instead.
> - **Bed colour: Elevation ↔ Exchange rate** (`SliceViewer::setBedColorMode`, Visualization-dock combo):
>   in rate mode the sand surface keeps its z_b geometry but is recoloured by the **net suspended flux
>   delivery − pickup** — blue = pickup/erosion, red = delivery/deposition, tan = balance. Same diverging
>   ramp + VBO/shader path as the elevation colouring; a mode-switched **"bed rate [mm/hr]"** legend with
>   its own auto-scaled ±half-range. MH's read: it's a **leading indicator** (rate leads state) — a preview
>   of where the bed is about to move before Δz_b accumulates.
> - **Captured from the REAL Exner update, not re-derived:** `morpho_exner` gained optional write-only
>   per-column `dep_out`/`ero_out` (applied deposition/erosion grain, post-MORFAC/post-limiter) — null-safe,
>   default null ⇒ byte-identical (M4/M5 `morpho_validation` pass null; `SeabedMorpho` writes never feed
>   physics). `SeabedMorpho::copy_exchange_host` returns the **physical (MORFAC-divided) rate** [m/s];
>   the worker publishes it per step (`disp_exch_`/`copyBedExchange`). ⚠ SCOPE: suspended coupling ONLY
>   (van Rijn dep − Winterwerp ero) — **bedload excluded**, so it complements Elevation mode (total bed
>   change, bedload included). Extended the `morpho_exner` GPU-vs-CPU parity test to cover the new outputs.
>   Headless smoke `--bed-colour <elevation|rate>` (offscreen add-sand conversion: 264 morpho steps, exit 0).
>   `bedstate.{h,cu}`, `seabed_engine.{h,cpp}`, `sim_worker.{h,cpp}`, `slice_viewer.{h,cpp}`,
>   `main_window.cpp`, `main.cpp`, `tests/test_bedstate.cu`. Nothing open from this session.

> **2026-07-09 session summary (tip `34f7699`; superseded as repo tip by `e3002dc` above):** two **GUI
> viewer-inspection controls** (GUI-only; NO physics — `gate_G1` + **93/93 unit** PASS, all M-gates
> untouched/byte-identical). One combined commit `34f7699` (the two features interleave in `main.cpp`/
> `main_window.cpp`/`CLAUDE.md`, and this env has no interactive `git add -p` to split hunks cleanly).
> - **Clip plane — "see inside hollow structures"** (`SliceViewer::setClipEnabled/setClipMode/
>   setClipFraction/setClipFlip`, `computeClipPlane`, `drawClipPlaneViz`; new `Camera::forward()/
>   targetPoint()`): a **Clip plane (see inside)** dock group (Enable / Orientation / Position / Flip)
>   that hides SOLIDS on the camera side — STEP mesh + voxel solids (both kinds) + seabed surface —
>   while the **flow slice + arrows stay visible** (MH chose "solids only"). Modes: axis-aligned
>   **X/Y/Z** (auto-oriented to hide the camera side via `sign(eye_a−p)`) + **Face camera** (normal =
>   view direction, slider pushes it in view depth); Flip swaps the hidden half. Implemented with
>   **`gl_ClipDistance[0]`**: `kVert`/`kMeshVert` write `dot(pos,uClipPlane)`, cut ONLY while
>   `GL_CLIP_DISTANCE0` is glEnabled per SOLID draw and **always disabled before the QPainter overlay**
>   (the depth-state gotcha class). Translucent cyan quad + outline shows the cut. Headless
>   `--clip <x|y|z|camera>,<frac>[,flip]` (windowed 4M-cell: all three modes render, ~45–52 fps, exit 0).
> - **Live MORFAC (bed speed-up)** (`SimWorker::setMorfac` → `SeabedMorpho::set_morfac`;
>   `MainWindow::setMorfacValue`, `morfac_spin_`): a Simulation-group `QDoubleSpinBox` (1–50×) that
>   retunes morphological acceleration LIVE so scour catches up faster (no reset). Worker-thread apply
>   (atomic + dirty), guarded `if (morpho_ && …)` so a value set on a fluid viewer stays pending and
>   lands when "Add sand" attaches a bed; reflects the loaded scenario's MORFAC on spawn (signals
>   blocked); no-op without a bed. MORFAC scales **bed time only** — the flow clock `t` is unchanged.
>   ⚠ Physically ≤10 steady / ≤5 reversing (RESEARCH §7); above that flow↔bed decouple + the per-step
>   |Δz_b|≤0.05h limiter clips (so higher ≠ proportionally faster). Headless `--morfac <M>` (add-sand
>   run: step-600 z_b range **0.062 m at M=12 vs 0.023 m at M=3** — the bed genuinely evolves faster —
>   budget err ~1e-6, stable, exit 0).
> - **Seed of the MORFAC ask (worth remembering):** MH asked whether the status `t = 46.8 s` means 46 s
>   of real progress. Answer: `t` is FLOW (hydrodynamic) time = Σdt; with morphology active the **bed has
>   aged MORFAC×t** (`bedstate.cu` `M*dt`; default seabed M=5 ⇒ ~234 s of bed). `sim fps` is compute
>   throughput, not physical time. Nothing open from this session.

> **2026-07-09 session summary (tip `3fd2ba3`; superseded as repo tip by `34f7699` above):** completed
> **M9 — tidal reversal + ranking campaign (`gate_M9` PASS)** and added its **interactive GUI driver**.
> - **Face-swap reversal is a fluid-core BC change** (`ChannelBC::flow_sign` ±1 + `ChannelFluidCore::
>   set_flow_direction`): +1 = the M2/M3 inlet@xmin / Orlanski-outlet@xmax (**default ⇒ byte-identical**);
>   −1 swaps the faces + re-points the pressure Dirichlet-p=0 pin (`ChannelMgpcg::set_outlet_dir`; the
>   `dir_xmax` int is now a signed face code). `U_inlet`/`ustar` stay magnitudes; the driver flips the
>   sign at slack (|U|→0). **The sediment side needed NO change** (bed shear/bedload follow the signed
>   velocity; suspended advection seals all faces; open-sea boundary is per-face upwind).
> - **gate_M9** (`b193507`, `tests/gate_m9_main.cpp`, `configs/m9_tidal.json`): a reversing cosine-slack
>   tide (MORFAC 4) over 4 slacks. All criteria PASS at 0.30 M cells: total sand conserved **~1e-14**
>   (CLOSED boundary ⇒ machine precision; the budget credits the structure-pin discard via new
>   `SeabedMorpho::structure_discard()`), **no monotonic KE growth** (matched-phase |U|=U_max KE ×1.00 —
>   NOT slack-KE, which is a bed-development transient), MORFAC limiter **clip-fraction 0 %** (new
>   `SeabedMorpho::clip_fraction()`). Ranking under reversal: **ring 2.70e-3 ≫ plate 7.3e-4 > cup 3.8e-4
>   m³** — the symmetric ring shelters its interior from BOTH tidal directions. Comparative (RESEARCH §11).
> - **Interactive GUI driver** (`3fd2ba3`, GUI-only): a **Tidal reversal (live)** dock section (enable +
>   Peak U_max / Plateau / Slack-ramp); the worker drives `U_d(t)` each step via the core/engine live
>   setters, greys "Input speed U" while active, shows the phase in the status bar. Headless `--tidal
>   Umax,plateau,ramp`. gate_G1 PASS.
> - **Verified byte-identical:** gate_M2 (St=0.1620, Cd=1.5184) + gate_M3 both re-gated PASS against a
>   fresh relink; **93/93 unit**; all targets build. The two new `SeabedMorpho` accessors are pure
>   diagnostics (no physics feedback ⇒ M4/M5/M6/M7/M8/seabed unaffected). Also committed a one-line
>   parallel GUI fix (`a061b1c`, box wireframe rebuild on resize — that authoring session was gone).
> - **The M6 calibration decision (below) is STILL the one open SUPERVISED item.** Everything in the
>   overnight M6/M7/M8 block below still stands. ⚠ The stale bits below ("uncommitted suspended work",
>   "tip c6377dc", "92/92") are SUPERSEDED — the open-sea/recycle sediment boundary, M8 porous-cup +
>   retention protocol + `gen_report.py`, and the Fast-sim throttle are all long committed (see `git log`).

> **2026-07-09 overnight session summary (read `M6_REPORT.md` first):** built + committed the
> **M6 pile-scour harness** (`gate_M6`), the **M7 Du-array deposition gate** (`gate_M7`), the **M8
> shape-ranking gate** (`gate_M8`, PASS), the **sub-grid porous momentum sink** (`channel_porous`),
> the **campaign machinery** (`run_protocol.h`), and the **HSV shear-multiplier** calibration knob.
> **93/93 unit tests pass. gate_M8 PASS; gate_M7 PASS (comparative); gate_M6 FAIL (absolute —
> supervised).** KEY RESULT: the M6 harness is correct (u\* pre-check PASS, stable on molecular ν,
> mass-conserving, right qualitative scour) but uncalibrated **S/D ≈ 0.13 is ~10× below Roulund's
> 1.25** — a structural HSV under-resolution (RESEARCH §11 #1), NOT closable by the sanctioned ±30 %
> sediment knobs (5× erosion → +2 %), by no-slip (+40 %), by the HSV multiplier (saturates ≈0.35),
> **and NOT by grid refinement to D/20 (0.185→0.189 at 23 M cells — the SL advection scheme, not the
> grid, is the limit).** M7/Du independently confirms the same ~9× absolute deficit yet PASSES its
> comparative gate (porous units reduce scour + trap sand). **SUPERVISED decision for MH — see
> `M6_REPORT.md` §5: accept comparative-ranking-only (recommended, ships now) or LBM plan-B for
> absolute fidelity. Nothing frozen; `configs/calibrated.json` = PROVISIONAL.** Everything below the
> "2026-07-08" line predates this session and remains valid.
>
> **PRELIMINARY PRODUCT SIGNAL (uncalibrated, COMPARATIVE — the tool's valid regime).** Ranking a
> flat plate / solid ring / solid cup / **porous cup (β≈0.30, the M7 screen model)** by matched-time
> trapped sand volume: the **porous cup traps ~3× the solid cup and ~5.5× the plate, with LESS edge
> scour than the solid cup** (0.034 vs 0.044 m). This is exactly the RESEARCH §9 design physics
> (permeable enclosures admit + settle sediment while decelerating flow; solid walls reflect flow and
> scour their toe) — i.e. the tool reproduces the COBOD self-ballasting hypothesis QUALITATIVELY
> (porous ≫ solid). ⚠ **BUT a porosity sweep (β = 0.10→0.65) shows trapping rising MONOTONICALLY with
> openness (1.6→6.1e-3) and edge scour falling monotonically — the model captures the SOLID-side
> penalty (reflection/toe-scour) but NOT the RESEARCH §9 too-OPEN limit (β<0.23 vortex street
> re-suspends the catch), which needs resolved wake shedding + a retention analysis it cannot do at
> this coarse LES resolution. So the tool currently OVER-favours openness and does not locate the
> β≈0.25–0.35 optimum** — flag before using it to pick a porosity. Fixing this needs the M8 retention
> protocol (rank on retained-under-storm, not just trapped) and/or finer/LBM shedding. Absolute
> numbers are indicative only (gate_M6 HSV deficit). Reproduce: `scour_m8_gate configs/m8_ranking.json
> --shapes plate,ring,cup,pcup --csv out.csv` (render `python tools/gen_report.py out.csv`);
> porosity sweep: `--shapes pcup --beta-por <β>`.
>
> ⚠ The **retention protocol** (`--retention`, PLAN M8 session-2) IS implemented and runs stably
> (incl. a live cosine U-ramp to a storm stage — verified stable + mass-conserving at 2.0 m/s), but
> with the **open-sea boundary it measures storm-ACCUMULATION, not catch-washout**: higher U supplies
> vastly more suspended load (the equilibrium `c_a` scales with u*), so both cup and porous-cup gain
> 38–43× under a 2.5× storm rather than losing their catch (retention ≫ 100 %). Testing true retention
> (the openness penalty) needs a **finite/CLOSED sediment supply** during the storm, or tagging the
> original catch — a follow-up. So the openness optimum is still not resolvable yet; do not read the
> `--retention` output as washout resistance until that setup is added.

## Your role
Lead **orchestrator** of WindCFD (`C:/CODE/WindCFD`). You do not write most code
yourself: you launch/verify **Opus 4.8** implementation agents (Agent tool for one feature,
Workflow for multi-step), verify **every** change yourself (rebuild + `ctest` + read the diff —
an agent's "passed" is a claim, not a fact), and apply judgment where agents must not (physics
calibration, M6/M7). Owner: **MH** (mh@cobod.com, COBOD) — an experienced C++/Qt engineer; talk
to him as a peer.

## Authority documents (read in order, before acting)
1. **CLAUDE.md** — hard rules + the verified build/test commands + a dense, per-feature "Build &
   test" log kept current by every agent. Read it fully; it is the most current technical map.
2. **PLAN.md** — milestones M0–M11 with ctest acceptance gates.
3. **RESEARCH.md** — physics spec (every constant verbatim); the 16 notes behind it are `research/`.
4. **Memory** (auto-loaded): `scour-project`, `user-profile`, `cutcell-sensitivity-study`,
   `mh-forms-visual-hypotheses`.

## First actions
1. Assess ground truth: `git -C C:/CODE/WindCFD log --oneline -20` (feature tip is **`ebe2d2c`**,
   the model-placement gizmo + arrow-width slider — see the latest summary; the bed-exchange-rate, clip-plane
   /MORFAC, and M9 blocks precede it); `git status` (**clean**); build + run the fast suite:
   `cmake --build build --config Release` then
   `ctest --test-dir build -C Release -L unit --output-on-failure` (should be **93/93** now).
2. The working tree is **clean** — everything is committed and pushed. No in-flight parallel work to
   coordinate (the old uncommitted-suspended-work thread is long resolved: the open-sea/recycle
   sediment boundary is committed).
3. There is **no blocking open thread**. The remaining PLAN milestones are **M10 (waves)**, **M11
   (perf/hero)**, **G2 (GUI v2)**, and **M8 session-2 leftovers** (batch matrix runner + checkpoint
   bit-repro test); the **M6 calibration** is the standing SUPERVISED decision for MH (thread #3). Also
   open: render-perf diagnosis (#2) and the cut-cell sensitivity study (#4).

## State at handover (2026-07-08 — historical; see the latest summary at top for current state)
Everything below is **committed and pushed** to `origin/master`
(github.com/COBOD-International/WindCFD, default branch **master**), tip **`c6377dc`**:
- **M0–M5 physics + G1 GUI** — all gate-verified (unchanged): MAC/MacCormack/MGPCG/Smagorinsky
  fluid; open channel + cylinder Cd/St; log-law bed-shear wall model + Jarrin SEM; suspended
  sediment (van Rijn pickup, w_s·c_b deposition, Rouse); bed morphodynamics (f_pack,
  Winterwerp/bedload, avalanche, Exner); Qt6 + GL 4.3 slice viewer.
- **Geometry pipeline in the GUI** (unchanged): STEP loading (OpenCascade → mesh); STEP → fluid
  **voxelization** (majority ray-parity fill, thin-wall safeguard); live obstacle via `--voxelize`
  / "Model as obstacle". Real unit: `Experiments/V000 Code tests/WindCFD V000_001.stp`.
- **Erodible-seabed live scenario** (`--scenario configs/seabed_v000.json`): the FULL live loop
  fluid→τ_b→suspended→bed→avalanche around the voxelized V000. **Still PRE-CALIBRATION / qualitative**
  (artificial low-Re `nu_fluid`, uncalibrated α/Cs, MORFAC=5) — NOT a measurement; real numbers = M6.

### Added in the `c55f045`→`95df6c4` batch (all committed + pushed)
All GUI/viz except the last, which is a physics-BC fix:
- **Input speed U** (Simulation group): sets the inlet current **LIVE, no reset**
  (`ChannelFluidCore::set_inlet_speed` → `SimWorker::setInletSpeed`); also folded into `GridOverride.U`
  so Apply keeps it. Headless smoke `--set-u <U>`.
- **Auto-scaling Δz_b legend**: the bed diverging map + its legend share one half-range that tracks
  the live peak `|z_b−z0|` (expand-fast / shrink-slow, 1 mm floor) so sub-mm early scour is visible;
  the legend now labels the **relative** Δz_b (absolute z_b≈z0 rounded away at 3 sig-figs).
- **Show bed** toggle (gates the seabed surface draw+fill+legend).
- **Two-colour voxel overlay**: "Show voxels" split into **Show solid voxels** (rigid structure/
  obstacle, dark-blue) + **Show sediment voxels** (sand, orange). The worker publishes a **KIND mask**
  (0 fluid / 1 sediment / 2 rigid) via `SimWorker::publish_solid_kind`.
- **Add sand** (Simulation group): fill the lower N m of every non-rigid column with sand (solids
  stay solid) and reset as a sediment run — converts a fluid viewer with a voxelized obstacle into a
  live scour run. `sim_setup.cpp build_seabed_from_structure` + the shared `assemble_seabed`.
  Headless smoke `--add-sand <m>`.
- **Log-law boundary-layer inlet (PHYSICS FIX)**: the top-hat inlet drove full-U near the bed at the
  first column → maximal τ_b → the sediment leading edge over-eroded. Fixed with a log-law BL profile
  referenced to the bed top (new `ChannelBC::bed_datum`, z0=d50/12, u* flux-matched to U). **Seabed runs
  now DEFAULT to log-law**; live "Inlet profile" toggle (Uniform ↔ Boundary layer). See Traps.
- **Arrow speed** slider [0,1] (Visualization): scales only the visual tracer gain (0 frozen, 1 = the
  default lively `1.5×`); physics untouched — calms the arrows on a fast current (MH's 8 m/s case).
- Verified: **gate_M2 PASS** (cylinder Cd/St — `2245 s`, cleared the outlet-fix re-gate), **gate_M3
  PASS** (log-law byte-identical), **gate_G1 PASS**, **80/80 unit**, seabed_flat PASS. V000 seabed
  runs at ~60 fps; log-law inlet cut the near-bed inlet velocity 0.50→0.369 m/s (τ_b −45%).
- Toolchain unchanged: MSVC 14.44 (**VS2022 Professional**), CUDA 13.1 arch 89, Qt 6.11.1, OCC 8.0.

### Added since (all committed + pushed; commits `d62fe6a` → `c6377dc`, tip **`c6377dc`**)
- **Load STEP = auto-obstacle; per-kind voxel overlay; flow-preserving Add sand** (`d62fe6a`):
  loading a STEP now auto-voxelizes + injects it as the obstacle (no toggle); the overlay is split by
  kind (rigid structure/obstacle dark-blue, sand orange); "Add sand" preserves the developed flow via
  in-place `update_solid` (no t=0 reset).
- **Editable domain + resolution "Apply" + seabed readout/resize fixes** (`bf3a7d8`): Simulation group
  gets editable Lx/Ly/Lz + voxel size h with a live cell-count readout and an **Apply** that rebuilds+
  resets at the new grid preserving every physics/scenario constant; a **converted** seabed (Add sand on
  a fluid viewer) now re-voxelizes the loaded model as the rigid STRUCTURE + re-applies the sand at the
  new grid — **this closes the old add-sand provenance gap (open thread #5)**. The seabed step/time status
  readout no longer freezes.
- **Scene save / load + auto-save checkpoints (`.scn`)** (`c6377dc`, THIS session): File menu **Save
  Scene As** (Ctrl+S) / **Load Scene** / **Auto-save** (Off/500/1000/2000/5000 steps). One self-contained
  `.scn` holds the whole setup + embedded STEP mesh + ALL live field data at the step; loading **restores
  and CONTINUES from that step** (flow keeps developing; a scoured bed comes back intact — not t=0). Only
  the STATEFUL device state is saved (velocities/pressure/live solid mask; seabed grain `G`, suspended `c`,
  EMA-filtered near-bed shear, structure); per-step scratch is recomputed on the next step. Big arrays are
  **float32** (≈half size, qualitatively lossless). A manual `<name>.scn` → **Recent Files** (immediate +
  prune-safe); an auto-save writes `<name>.<step>.scn` (a restart point, not Recent; Load offers a picker
  over the siblings). The worker gathers state **on its own thread (race-free)** → queued `checkpointReady`
  → main-thread write (no 60 fps stall). New `src/gui/scene_io.{h,cpp}` (Qt-free); core/engine
  `copy_state_host`/`load_state_host`/`save_state`/`load_state`; `SimWorker::requestCheckpoint/
  setAutosaveInterval/primeCounters`; `MainWindow::saveSceneNow/loadSceneFromPath/restoreScene`. CLI smoke
  `--save-scene`/`--load-scene`/`--autosave`. **`*.scn` + `Scenes/` are gitignored** (a full 4M-cell scene
  is ~90 MB). Verified: fluid save@417→restore→945; seabed convert→save→restore resumes the scoured bed
  (mass err ~1e-6); auto-save series; numbered-checkpoint load resumes at its exact step; **80/80 unit +
  gate_G1 PASS**.

## Open threads (where we stopped)
1. **✅ CLOSED — outlet-fix re-gate.** `gate_M2` (Cd/St) and `gate_M3` both re-run and PASS this
   session, confirming commit `2f5e6bb` (outlet clamp) + all the new work leave M2/M3 intact.
   Nothing to do here.
2. **Render performance (OPEN, mid-diagnosis, unchanged).** The viewport goes **choppy/jerky
   intermittently** — even with just cylinder + a model (NOT seabed), even hiding all geometry, even
   on camera orbit — then self-resolves to 60 fps. MH ruled OUT the voxel overlay. Working diagnosis:
   the render is **GPU-starved by the 4M-cell MGPCG solve** (one physical GPU; decoupling the render/
   sim *threads* does NOT decouple the GPU; MGPCG iters vary with flow → intermittent). Prime added-
   load suspect: the **flow arrows** (worker D2H-copies the whole {u,v,w} field ~12 Hz + main-thread
   advection). NEXT: MH to test **arrows-OFF** and **coarse grid (h=0.10 via Apply)**; if confirmed,
   move the arrow field-copy/advection off the render hot path and/or give the render stream priority.
   Coarse-grid is the practical workaround. (The new **Arrow speed** slider only changes visual pacing,
   not the D2H load, so it does NOT address this.)
3. **M6/M7/M8 harnesses BUILT + committed (2026-07-09); M6 calibration is the open SUPERVISED
   decision.** `gate_M6` (Roulund pile scour), `gate_M7` (Du perforated-unit deposition, sub-grid
   porous model), `gate_M8` (shape ranking — **PASS**), campaign machinery (`run_protocol.h`), and
   the HSV shear-multiplier knob are all in. **Read `M6_REPORT.md`.** The blocker is physical, not
   code: uncalibrated S/D ≈ 0.13 (~10× under 1.25) from the SL solver's HSV under-resolution
   (RESEARCH §11 #1). The sweeps prove the sanctioned sediment knobs can't fix it (5× erosion → +2 %;
   no-slip pile → +40 %; HSV multiplier moves it sub-linearly — even ×8–12 falls short). **MH's
   decision (M6_REPORT §5):** (1) accept comparative-ranking-only (the tool's authority; gate_M8
   already ranks correctly), (2) refine grid to D/20 (`--h 0.005`), (3) freeze an HSV multiplier as a
   documented crutch, or (4) LBM plan-B. **Nothing frozen; `configs/calibrated.json` = PROVISIONAL.**
   Re-run M1–M5 before freezing anything (all additions are byte-identical by default; 93/93 unit
   incl. M5 parity). Remaining M7 work: a full (non-quick) Du run for the C-AR/H-AR ranking (the
   harness is mechanically verified — porous k correct, stable, mass-conserving, deposits sand).
   Architecture note (M6_REPORT §6): the CLOSED-mode Exner erosion uses Winterwerp `erosion_coeff`
   (default 0.018, now exposed), NOT the van Rijn α (α drives the M4 suspended path + open-sea
   boundary) — decide how "α ±30 %" maps to this pipeline.
   Remaining PLAN milestones: **M10 (waves)**, **M11 (perf/hero)**, **G2 (GUI v2)**, and M8 session-2
   leftovers (batch matrix runner + a checkpoint save→restore→step bit-repro test; the HTML report +
   retention protocol are DONE). ✅ **M9 (tidal reversal + campaign) is DONE** (`gate_M9` PASS, tip
   `3fd2ba3`) — see the latest summary at the top.
4. **Cut-cell sensitivity study** — deferred (memory `cutcell-sensitivity-study`): after M6/M7,
   rank shapes at h vs h/2; build the Batty cut-cell projection only if rankings flip.
5. **✅ CLOSED — add-sand provenance on grid Apply** (`bf3a7d8`). Apply on a converted seabed now
   re-voxelizes the loaded model as the rigid structure + re-applies the stored sand depth at the new
   grid (tracked by `added_sand_m_`), instead of reverting to the sand-less config. Nothing to do.
6. **✅ CLOSED — the parallel uncommitted suspended-sediment work is committed.** It became the
   **open-sea / recycle sediment boundary** (`SED_BC_OPEN`/`SED_BC_RECYCLE`, `suspended_boundary_flux`,
   `equilibrium_cb`/Rouse ghost) shipped in `f8e91f6` + `8cc5e3e`, plus the GUI "Sediment BC" combo. Tree
   is clean. Nothing to coordinate.
7. **M9 GUI tidal driver — possible follow-ups (nice-to-have, not blocking).** Peak U_max is an explicit
   spin (does not auto-seed from the loaded config's U — a deliberate predictability call; revisit if MH
   wants auto-seed). MORFAC ≤ 5 for reversing flow is a tooltip note, not enforced. A period+amplitude
   control layout (vs plateau+ramp) was offered as an alternative.
8. **Model-placement gizmo — interactive GL/drag unverified (OPEN, this session).** The whole
   placement→voxelize→Apply pipeline is verified headless, but the gizmo's on-screen DRAW + mouse DRAG feel
   were never exercised (offscreen has no GL). MH needs a windowed pass: confirm all handles render at the
   1/3 size, the scale cubes are grabbable, and each drag direction/sensitivity feels right. Every knob is a
   single constant in `slice_viewer.cpp` — `gizmoSize()` (overall size), `kTransLen`/`kRingRad`/`kScalePos`/
   `kCubeHalf` (handle layout), the tier-1/2 pick radii (px), and the drag rates (translate 1:1, rotate
   ray-plane, scale `exp(along·0.012)`). If a rotate/scale direction reads inverted it's a sign flip there.

## How to work here (MH's norms — also in memory)
- **Delegate to Opus agents, verify everything yourself** (build + ctest + git-diff scope). Efficient
  verify: skip re-running prior gates a git-diff proves are unchanged. **Long gates exceed the 10-min
  Bash timeout → background + poll**, never a blocking call. ⚠ Running windowed GUI smokes WHILE a gate
  runs shares the one GPU and stretches the gate's wall-clock (gate_M2 took 37 min this session vs the
  ~25 min nominal for exactly this reason) — schedule them apart.
- **Push cadence:** MH wants verified work pushed; it's the **company repo** (outward-facing). He
  sometimes drives background agents directly — verify + push whatever lands. Commit or push only when
  MH asks (he does ask, per-batch).
- **Debugging: build MH the observability/viz tools and let HIM analyze the progression / form the
  hypothesis — don't race ahead theorizing** (memory `mh-forms-visual-hypotheses`). For physics
  **calibration** (M6/M7) you apply your own judgment *with* him. The leading-edge-scour fix this
  session followed this: MH observed the artifact, I traced the BC and confirmed + fixed it.
- Git identity is "Claude Agent"; **Conventional Commits**; commits end with a
  `Co-Authored-By:` trailer naming the running model (this session: `Claude Opus 4.8` per the
  system prompt — history also shows earlier `Claude Fable 5`; use whichever model is running).

## Traps (each cost a real debugging session)
- **VS18/14.51 is installed and CUDA 13.1 rejects it** — use VS2022 Professional 14.44.
- **`max_ustar` (bed friction) is NOT the flow-blow-up gauge** — the wall model rate-limits/clamps
  it, so it stays ~0.03 even when the free-stream flow blows up. Log **`max|u|` of the velocity
  field** to detect a flow instability.
- **Top-hat inlet over-shears the sediment LEADING EDGE (fixed this session).** A uniform inlet
  imposes full U at the bed at i=0 (no boundary layer), so the wall model probes ≈U there → maximal
  τ_b at the sediment's upstream edge → it scours away too fast (i=0 is pinned to the Dirichlet inlet
  so it never develops a BL). Fix = the **log-law BL inlet** (default for seabed; profile referenced to
  the bed top via `ChannelBC::bed_datum`). ⚠ Keep `bed_datum` defaulting **0** so M2 (uniform, early-
  returns) and M3 (loglaw, z−0) stay byte-identical — the gate depends on it.
- **The published voxel mask is now a KIND mask (0 fluid / 1 sediment / 2 rigid)**, not binary —
  `SimWorker::publish_solid_kind`. Consumers using `!=0` still work; the two-colour overlay + Add-sand
  read the codes. The auto-range solid-exclusion (`disp_solid()`) is a SEPARATE binary snapshot.
- **Outlet BCs**: a convective/Orlanski outlet must clamp outflow-only and must NOT use a
  multiplicative flux-rescale that can amplify reversed inflow (fix `2f5e6bb`, gate-confirmed).
- **Render vs sim share the GPU** — decoupling the *threads* ≠ decoupling the GPU; a heavy 4M-cell
  solve can starve the render even if dt is fine (see open thread #2).
- Deposition is never thresholded (θ_cr gates erosion/bedload only). Exner: E takes /ρs, D=w_s·c_b
  does NOT. OCC STEP is millimetres → ×0.001. All GL on the main thread; sim on the worker.
- A model larger than the domain clips → near-total blockage → the flow can blow up; size the
  domain to the model.
- **Saved scenes (`.scn`) are large** (a full 4M-cell grid ≈ 90 MB even at float32) — `*.scn` +
  `Scenes/` are **gitignored**; never commit one. A `.scn` restores through the normal build path then
  injects fields, so it only saves the STATEFUL device state; if you add a new stateful buffer to the
  core/engine, extend `scene_io` + `save_state`/`load_state` or a restore will silently start it at zero.
- **Transient AUTOMOC build crash**: `CL.exe exited with code -1073741819` (0xC0000005) on the moc TU is
  a known MSVC 14.44 flake — just re-run `cmake --build`. Not a code bug.
