# M6 — Pile-scour calibration report (Roulund V7)

**Status: harness COMPLETE and verified; calibration is a SUPERVISED decision for MH (this report
presents the data + options, and freezes NOTHING).** Generated 2026-07-09 by the overnight agent
session. All runs: `configs/m6_roulund.json` via `build/Release/scour_m6_gate.exe` (`gate_M6`).

---

## 1. TL;DR

- The M6 harness works end-to-end and is honest: the **frozen-bed u\* pre-check PASSES** (recovered
  ambient u\* = 0.0196 m/s vs the 0.020 target, 2.4 % error), the flow is **stable on molecular ν**
  (no artificial viscosity needed — `maxu ≈ 1.7·U`, MGPCG 5–17 iters), sand mass is conserved to
  ~1e-7, and the **scour signature is qualitatively correct** (upstream deeper than downstream, hole
  forms around the pile, upstream slope steepens toward repose).
- **The uncalibrated equilibrium is S/D ≈ 0.13 — about 10× below Roulund's 1.25.** This is the
  RESEARCH §11 #1 limitation (semi-Lagrangian solvers under-resolve the horseshoe vortex and
  under-predict its bed shear), here far more severe than the ~30 % published RANS shortfall because
  the SL advection is strongly diffusive at D/10.
- **The sanctioned ±30 % sediment knobs cannot close a 10× gap, and the sweeps prove it directly:**
  a **5× increase in the erosion coefficient changed S/D by <3 %** (0.131 → 0.134). The equilibrium
  scour depth is set by the near-pile **flow shear**, not the sediment erosion rate. So α / erosion
  are the wrong levers; the lever is **τ_b** (the Shields stress in the scour ring).
- Two τ_b levers were exercised: the **no-slip pile** BL (0.131 → 0.185, +40 %) and the sanctioned
  **HSV shear multiplier** (results in §4). Neither is a model *default* — both are MH's calibration
  choice.

**Decision required from MH (§5):** accept comparative-ranking-only (the tool's stated authority,
RESEARCH §11 #2), or promote the LBM plan-B (RESEARCH §11 #6) if absolute S/D is required. **Grid
refinement to D/20 was tested and does NOT help** (S/D 0.185→0.189) — the limit is the SL advection
scheme, not resolution. An HSV multiplier can be frozen as an explicit crutch but saturates below
1.25. **No constant has been frozen.**

**Cross-benchmark confirmation (M7 / Du 2025):** the same harness on the Du perforated-unit case
gives control scour 0.7 cm vs measured ~6 cm — the *same* ~9× under-prediction — yet **gate_M7 PASSES
its comparative gate**: both porous units reduce scour vs the bare pile (C-AR 6 %, H-AR 10 %) and
trap sand, stable and mass-conserving (1e-16). The absolute-magnitude deficit is systematic and
consistent; the comparative behaviour (the tool's authority) is sound.

---

## 2. Setup (as built)

| Quantity | Value | Note |
|---|---|---|
| Domain | 3.0 × 1.6 × 0.60 m | 0.40 m water (= Roulund h_dom) **+ 0.20 m sand reservoir** for the scour hole |
| Grid | 300 × 160 × 60 = **2.88 M cells**, h = 0.01 m (= D/10) | ⚠ PLAN/CLAUDE say "≈19 M cells" — that is a ~7× error; the physical D/10 grid is 2.9 M. Flagged, docs to fix. |
| Pile | D = 0.1 m (10 cells) at (1.0 m, 0.8 m) | 10D from inlet, 7.5D lateral clearance |
| Flow | V = 0.46, d50 = 0.26 mm, ρ = 1027, ρs = 2650, ν = 1.36e-6 | live-bed V/Vcr = 1.25; θ_ambient ≈ 0.108 |
| Inflow | flux-matched **log-law BL** inlet (u\* = 0.0209), molecular ν + Smagorinsky Cs = 0.11 | matches Roulund's steady-inlet RANS setup; SEM/precursor available but deferred |
| Bedload | Engelund–Fredsøe + suspended load ON (RESEARCH V7) | |
| Pile surface | free-slip (RESEARCH §3 default) unless `--noslip` | |

**Protocol** (research/14, /16): freeze bed → spin up 5 flow-throughs → verify u\* → ramp MORFAC 1→M
over 2 flow-throughs → run ≥0.5·T of morphological time → fit S(t) = S_eq(1−e^(−t/T)) and Welzel
a(1−1/(1+bt)) → report S/D. T_pred = 130 s (research/11).

### Inflow choice (documented deviation)
RESEARCH V7 prefers a **3-D turbulent precursor**; this harness uses the **flux-matched log-law mean
inlet** (which reproduces the ambient u\* = 0.020 exactly, as the pre-check confirms) — the same
steady-inlet condition Roulund's own k-ω RANS used. The SEM generator (`sem_inlet`, M3-gated) and
`precursor.*` record/replay are wired and available; adding ambient TI is a follow-up. Since the
scour deficit is a **mean-HSV** problem (§3–4), inlet turbulence is not the bottleneck.

---

## 3. Results — the calibration sweeps

MORFAC and run-length differ between the M=1 baseline and the M=4 sweep rows; equilibrium S/D is
MORFAC-independent (M only sets how fast equilibrium is reached), so the rows are comparable. All
rows PASS the u\* pre-check, are stable, and conserve mass.

| Run | pile BC | erosion_coeff | S/D (exp fit) | S/D (Welzel) | note |
|---|---|---|---|---|---|
| **baseline** (M=1) | free-slip | 0.018 | **0.131** | ~0.18 | reference; T_fit = 27.7 s |
| erosion 5× (M=4) | free-slip | 0.090 | 0.134 | 0.186 | **+2 % — erosion-rate insensitive** |
| erosion 3× (M=4) | free-slip | 0.054 | 0.132 | 0.184 | +1 % |
| no-slip (M=4) | **no-slip** | 0.018 | **0.185** | 0.275 | +40 % — the pile BL/HSV matters |
| **D/20** (h=0.005, 23 M, M=6) | no-slip | 0.018 | **0.189** | 0.279 | **grid ×2 finer ⇒ NO gain over D/10 no-slip (0.185)** |
| HSV ×4 (M=4) | free-slip | 0.018 | 0.265 | 0.323 | τ_b ring multiplier — see §4 |
| HSV ×8 (M=4) | free-slip | 0.018 | 0.332 | — | saturating |
| HSV ×12 (M=4) | free-slip | 0.018 | 0.354 | — | saturating (≪ 1.25) |

**The headline finding:** scaling the erosion rate 3–5× barely moves equilibrium S/D. Equilibrium is
where the *local τ_b (Shields)* drops to threshold as the hole deepens — a **flow** balance, not a
sediment-rate one. Confirmed independently by the no-slip result (a pure-flow change, +40 %).

> The Welzel fits sit ~40 % above the exponential fits because the 0.5T runs are short (fat-tail
> extrapolation). A converged calibration run should use ≥1T; reported here for the sweep only.

The T-band gate also fails (T_fit ≈ 26–35 s vs the 65–260 s band): a shallow hole reaches its (low)
τ_b balance quickly, so the fitted timescale is short. This is a *symptom* of the weak HSV, coupled
to the S/D deficit — not an independent problem.

---

## 4. HSV shear-multiplier demo

The sanctioned optional knob (RESEARCH §11 #1, `--hsv X`): τ_b is amplified by ×X in the near-pile
scour ring [R, 2.5R] (`SeabedMorpho::set_shear_multiplier`, GPU-vs-CPU parity-tested). Because the
deficit is a τ_b deficit, this is the lever that *does* move equilibrium S/D. The bracket sweep
quantifies how much amplification reaches the 1.25 target:

| HSV ×X | S/D (exp fit) | vs baseline |
|---|---|---|
| 1 (none) | 0.131 | — |
| 4 | 0.265 | ×2.0 |
| 8 | 0.332 | ×2.5 |
| 12 | 0.354 | ×2.7 |

**The multiplier works but SATURATES far below the target.** Tripling the amplification (×4→×12)
raises S/D only ~34 % (0.265 → 0.354), and the curve is flattening well under 1.25. Reason: a
per-column τ_b boost in a
fixed [R,2.5R] ring deepens the hole until the **repose-angle avalanche** and the finite amplified
zone re-balance — a per-column shear multiplier cannot manufacture the missing **3-D flow structure**
(the deep HSV downflow that the SL solver smears). Extrapolating, reaching 1.25 would need a
physically-absurd multiplier (≫10×), i.e. it is effectively unreachable this way.

**Conclusion:** the scour deficit is a **flow-model** limitation (RESEARCH §11 #1/#6), not a tunable
sediment constant. The honest fidelity paths are **grid refinement** (resolve the HSV — §5 option 2)
or the **LBM plan-B** (§5 option 4). The HSV multiplier stays available as an explicit, documented
calibration crutch if MH accepts comparative-fidelity-with-a-crutch, but it cannot honestly reproduce
the absolute S/D = 1.25 on its own.

---

## 5. Options / decision for MH (calibration is yours, not the agent's)

1. **Accept comparative-ranking-only.** The tool's stated authority is *relative* shape ranking
   under identical conditions (RESEARCH §11 #2); gate_M8 already demonstrates a correct, stable,
   distinguishable ranking. Absolute S/D is flagged indicative. Lowest risk; ships now.
2. **Refine the grid to D/20 (h = 0.005, ~23 M cells). ⚠ TESTED — it does NOT help.** S/D went
   0.185 (D/10 no-slip) → **0.189 (D/20 no-slip)**: no meaningful gain at 8× the cells. The deficit
   is the **semi-Lagrangian advection scheme** (its numerical diffusion smears the HSV at any
   resolution), not the grid. So option 2 is largely ruled out — do not expect resolution alone to
   close the gap. This makes options 1 and 4 the realistic paths.
3. **Freeze an HSV multiplier** at the value §4 shows reaches 1.25 (a calibration crutch, honestly
   documented). Then re-run M1–M5 (they are byte-identical by construction — see §6).
4. **Promote the LBM plan-B** (RESEARCH §11 #6): wall shear is local & 2nd-order in LBM. Largest
   effort; the right move if 2–3 fail.

I recommend **(1) for shipping the ranking tool now** (gate_M8 + gate_M7 already demonstrate correct,
stable, distinguishable *comparative* behaviour), with **(4) the LBM plan-B** as the fidelity path if
absolute S/D accuracy is later required. **(2) is ruled out** (D/20 tested, no gain). **(3)** stays
available only as an explicitly-labelled crutch. **I have frozen nothing.**

`configs/calibrated.json` is written as **PROVISIONAL / UNCALIBRATED** (nominal constants + this
report referenced) so M8's comparative ranking has a config to consume without implying calibration.

---

## 6. Architecture findings (flagged for the record)

- **α does not drive CLOSED-mode scour.** The seabed Exner erosion uses **Winterwerp u_lift**
  (coeff 0.018), not the van Rijn α-pickup (α drives the M4 *suspended* path + the open-sea
  boundary). RESEARCH names α the headline knob, but in this architecture the exposed erosion knob is
  `MorphoParams.erosion_coeff` (default 0.018, byte-identical). Decide whether to (a) keep Winterwerp,
  or (b) route van Rijn pickup into the Exner (a physics change requiring M5 re-gating). Moot for the
  scour magnitude given §3, but relevant to how "α ±30 %" is interpreted.
- **Grid count doc error:** M6 is 2.9 M cells at D/10, not "≈19 M" (PLAN §4 / CLAUDE.md).
- **No earlier gate was traded:** the erosion_coeff and HSV additions default to byte-identical
  behaviour; the full unit suite is **92/92** including the M5 morpho parity tests. The M1–M5 gates
  should still be re-run before any constant is frozen (HANDOVER §5), but no committed change alters
  their inputs.

---

## 7. Reproduce

```
# baseline (M=1, ~20-40 min):
build/Release/scour_m6_gate.exe configs/m6_roulund.json --tag baseline --csv baseline.csv
# sweeps:
build/Release/scour_m6_gate.exe configs/m6_roulund.json --morfac 4 --trun-frac 0.5 --erosion 0.09   # erosion 5x
build/Release/scour_m6_gate.exe configs/m6_roulund.json --morfac 4 --trun-frac 0.5 --noslip          # no-slip pile
build/Release/scour_m6_gate.exe configs/m6_roulund.json --morfac 4 --trun-frac 0.5 --hsv 8           # HSV x8
# fidelity experiment (D/20):
build/Release/scour_m6_gate.exe configs/m6_roulund.json --h 0.005 --morfac 4 --trun-frac 0.5
```
