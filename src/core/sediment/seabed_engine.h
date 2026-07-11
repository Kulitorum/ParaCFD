// seabed_engine.h — the live erodible-seabed morphodynamic loop, assembled for the first time on
// TOP of the live fluid core. Given the post-step fluid fields {u,v,w,ν_t} it advances one
// morphological cadence: near-bed τ_b (wall model) → suspended sediment (settling advection +
// ν_t/σ_s diffusion) → bed Exner update (deposition + Winterwerp erosion + bedload divergence,
// ×MORFAC, |Δz_b|≤0.05h) → avalanche (32/30°) → confine c to the fluid → clamp/pin the rigid
// structure. When the bed has moved > threshold since the last flow-mask build it returns the
// updated solid mask (sand cells with z<z_b OR structure) so the caller re-masks the flow core.
//
// PRE-CALIBRATION DEMO: constants are NOMINAL (α, Cs, MORFAC from config, uncalibrated). This is a
// qualitative scour visualisation, NOT a measurement — α/Cs/HSV calibration is the SUPERVISED M6
// milestone. All gate-verified kernels (suspended.*, bedstate.*, avalanche.*, bedshear.*) are
// reused UNCHANGED; only the three seabed_morpho.* helpers bridge the elevated/moving bed.
//
// Qt-free + GL-free (pure libscour): it drives the *_gpu launchers, owns its device buffers, and
// consumes only the FluidCore's {u,v,w,ν_t} device pointers (PLAN §3). Not thread-safe by itself —
// the GUI runs it entirely on the sim worker thread (never the main/GL thread).
#pragma once

#include "core/fluid/bedshear.h"     // WallParams
#include "core/fluid/mac_grid.h"
#include "core/sediment/bedstate.h"  // MorphoParams
#include "core/sediment/avalanche.h" // AvalancheParams

#include <cstdint>
#include <vector>

namespace scour::core
{
	// Suspended-sediment domain boundary on the streamwise (x) faces (RESEARCH §6; CLAUDE.md).
	// CLOSED   — the M4/M5 zero-flux box (mass-exact; the DEFAULT so the gates stay byte-identical).
	// OPEN     — open-sea equilibrium: the inlet carries the Rouse equilibrium load (this model's own
	//            pickup=deposition fixed point via equilibrium_cb) and the outlet lets it leave freely;
	//            tidal reversal handled by the per-face upwind (suspended_boundary_flux).
	// RECYCLE  — flux-matched recirculating flume: the outlet load is fed back to the inlet, scaled so
	//            injected mass exactly equals extracted mass (mass-conserving cross-check / lab analogue).
	enum SedBoundary : int { SED_BC_CLOSED = 0, SED_BC_OPEN = 1, SED_BC_RECYCLE = 2 };

	// Config-level knobs for the seabed scenario (all NOMINAL / uncalibrated for the demo).
	struct SeabedParams
	{
		double sand_depth = 1.0; // m — initial flat sand reservoir thickness (z_b0)
		double d50 = 0.2e-3;     // m
		double rho = 1027.0;     // kg/m³
		double rho_s = 2650.0;   // kg/m³
		double nu = 1.36e-6;     // m²/s
		double Cs = 0.11;        // Smagorinsky (used by the fluid core; recorded here for logs)
		double alpha = 0.00033;  // van Rijn pickup coefficient (NOMINAL ±30%); drives the M4 suspended
		                         // path + the open-sea equilibrium boundary — NOT the bed-Exner erosion
		double erosion_coeff = 0.018; // Winterwerp u_lift coeff — the de-facto bed-Exner erosion-rate
		                              // calibration knob in CLOSED mode (RESEARCH §7 default 0.018)
		double morfac = 5.0;     // MORFAC M (≤10 steady current)
		double sigma_s = 0.7;    // suspended Schmidt number (RESEARCH §6.2 production)
		int bedload_formula = 1; // 0 = Wong–Parker, 1 = Engelund–Fredsøe (Roulund form)
		int diffusion_on = 1;    // ν_t/σ_s suspended diffusion (needs the core's ν_t)
		double remask_frac = 0.5; // rebuild the flow mask when max|Δz_b| > remask_frac·h
		int sed_bc = SED_BC_CLOSED; // suspended x-face boundary (0 closed / 1 open-sea / 2 recycle)
		double U_inlet = 1.0;       // m/s inlet current — sets the open-sea equilibrium u* (log-law inversion)
	};

	class SeabedMorpho
	{
	public:
		// `g` is the full flow grid; `structure` is the rigid printed-unit solid mask (size
		// p_count, 1=structure), or empty for a bare seabed. Allocates + initialises the bed
		// reservoir, suspended field (clear water), and the structure-footprint pin/mask.
		SeabedMorpho(MacGrid g, SeabedParams sp, const std::vector<unsigned char>& structure);
		~SeabedMorpho();
		SeabedMorpho(const SeabedMorpho&) = delete;
		SeabedMorpho& operator=(const SeabedMorpho&) = delete;

		// The initial flow solid mask = sand reservoir (z < sand_depth) OR structure. The flow core
		// must be constructed with this so the bed + structure are solid to the flow from step 0.
		std::vector<unsigned char> initial_flow_solid() const;

		// Advance one morphological cadence from the live fluid fields (device pointers). `nut` may
		// be null (⇒ no suspended diffusion). Returns true and fills `new_solid` when the bed moved
		// enough to require a flow-mask rebuild. dt is the fluid step's dt (MORFAC is applied inside).
		bool step(const double* u, const double* v, const double* w, const double* nut,
			double dt, std::vector<unsigned char>* new_solid);

		// --- Checkpoint save / restore (scene_io) --------------------------------------------
		// The sediment params this engine runs with (for the scene header + a rebuild after restore).
		SeabedParams params() const { return sp_; }
		// Copy the STATEFUL bed state to host: grain thickness G (per-column, ncol), suspended
		// concentration c (per-cell, np), the EMA-filtered near-bed shear ema_x/ema_y (per-column,
		// carries the wall-model time-filter memory) and the morpho step count. Everything else is
		// per-step scratch, recomputed on the next step().
		void save_state(std::vector<double>& G, std::vector<double>& c,
			std::vector<double>& emax, std::vector<double>& emay, long long& nstep) const;
		// Restore the stateful bed state (H2D). Sizes must match (ncol for G/emax/emay, np for c),
		// else that array is left at its constructed value. Resets the remask baseline G_last = G and
		// the morpho step count. The structure mask is set at construction (pass it to the ctor).
		void load_state(const std::vector<double>& G, const std::vector<double>& c,
			const std::vector<double>& emax, const std::vector<double>& emay, long long nstep);

		// --- live switches (call on the worker thread, between step()s) -----------------------
		// Freeze morphology: when set, step() computes ONLY the near-bed grain-skin τ_b / u* (so the
		// wall model can be probed) and returns false WITHOUT advancing suspended/Exner/avalanche or
		// rebuilding the flow mask (nstep_ stays put). Used by the M6 pre-check: freeze the bed, spin
		// up the flow, and verify the recovered ambient u* before enabling the Exner update
		// (RESEARCH V7 / research/14). Default OFF ⇒ existing gates byte-identical.
		void set_morphology_frozen(bool f) { morph_frozen_ = f; }
		bool morphology_frozen() const { return morph_frozen_; }
		// Live MORFAC (RESEARCH §7): ramped 1→target over 1–2 flow-throughs after the frozen spin-up
		// (research/16 E8). M multiplies morphological time only; the wall-model EMA/cadence do not.
		void set_morfac(double m) { mp_.morfac = m; sp_.morfac = m; }
		double morfac() const { return mp_.morfac; }
		// M6 HSV shear-multiplier field (per-column, size nx·ny): scales the near-pile τ_b to restore
		// the SL-under-resolved horseshoe-vortex amplification (RESEARCH §11 #1; seabed_apply_shear_mult).
		// SUPERVISED calibration knob — set by the gate/orchestrator, never a model default (null ⇒ 1).
		void set_shear_multiplier(const std::vector<double>& mult);
		// Select the suspended x-face boundary (SedBoundary). Rebuilds the open-sea equilibrium ghost.
		void set_boundary_mode(int mode);
		// New inlet current: re-derives the open-sea equilibrium inflow (u*, Rouse profile) for U.
		void set_inlet_speed(double U);
		// PLAN G3.3 (drop/settle SINK coupling): swap the rigid structure mask mid-run — the dropped
		// blocks re-settled into the developing scour hole and were re-voxelized. Updates the device
		// structure + the per-column frozen footprint (columns a block now occupies are pinned; ones it
		// left resume eroding) and RETURNS the new flow solid mask (sand ∪ new structure) so the caller
		// re-masks the flow core in place (update_solid). A mismatched-size mask is ignored (returns
		// empty). Additive — the gates never call it, so they remain byte-identical.
		std::vector<unsigned char> set_structure(const std::vector<unsigned char>& structure);
		// Net suspended sand [m³] imported across the open x-boundaries since construction (Σ inflow −
		// outflow). ~0 in CLOSED/RECYCLE; the open-sea supply/deficit in OPEN mode (open-budget audit).
		double boundary_sand_in() const { return boundary_sand_in_; }

		// MORFAC bed-change limiter (RESEARCH §7 / research/16): the Exner update clamps |Δz_b| ≤
		// dz_limit_frac·h per morphology step. This is the running fraction of column-updates that HIT
		// that clamp since the last reset — the diagnostic the tidal gate (M9) checks stays < 1% (a high
		// clip-fraction ⇒ MORFAC too large for the flow). Counts every column each morpho step (inactive
		// columns never clip), matching how production codes report the DzMax hit-rate. Pure diagnostic,
		// no physics feedback ⇒ existing gates are unaffected.
		double clip_fraction() const { return clip_colsteps_ > 0 ? clip_sum_ / (double)clip_colsteps_ : 0.0; }
		void reset_clip_stats() { clip_sum_ = 0.0; clip_colsteps_ = 0; structure_discard_ = 0.0; }

		// Cumulative sand [m³] REMOVED from the domain by the bed post-process (seabed_bed_post): the
		// [0,Gmax] clamp + the rigid-structure-footprint pin. A rigid slab's footprint cannot accrete
		// sand in this single-interface bed, so deposition landing there is discarded — a physical sink,
		// not a transport leak. Positive ⇒ net removed (add it back to close the budget: V_domain +
		// structure_discard == V0 + boundary_sand_in). Lets the tidal mass gate credit structure pinning
		// exactly (analogous to ± boundary fluxes) so the residual is the true test of reversal/transport
		// conservation. Reset to 0 at construction (and by reset_clip_stats for a clean post-spin-up start).
		double structure_discard() const { return structure_discard_; }

		// --- snapshots / diagnostics (host reads; call on the worker thread) ------------------
		void copy_zb_host(std::vector<float>& zb) const; // per-column bed elevation z_b [m], size nx·ny
		void copy_ustar_host(std::vector<float>& us) const; // per-column friction velocity u* [m/s], size nx·ny (last step's wall model; for the M6 u* pre-check)
		// Per-column NET suspended bed-exchange rate [m/s] from the last step: (delivery − pickup) =
		// (dep − ero)/(c_pack·M·dt), the PHYSICAL (MORFAC-divided, real-time) bed-elevation velocity.
		// +ve ⇒ net deposition (sand arriving), −ve ⇒ net pickup (erosion). Size nx·ny; all-zero before
		// the first morphological step or while the bed is frozen. GUI diagnostic (bed-surface colouring).
		void copy_exchange_host(std::vector<float>& rate) const;
		double bed_volume() const;   // Σ G·h²  [m³ solids]
		double susp_volume() const;  // Σ c·h³  [m³ solids-equivalent]
		double max_ustar() const { return last_max_ustar_; }
		double max_tau() const { return last_max_tau_; }
		double zb_min() const { return last_zb_min_; }
		double zb_max() const { return last_zb_max_; }
		long long morpho_steps() const { return nstep_; }
		MacGrid grid() const { return g_; }
		double sand_depth() const { return sp_.sand_depth; }
		const double* c_device() const { return c_; } // device suspended concentration field (GUI viz snapshot source)
		// The rigid structure solid mask (host, p_count, 1=structure) — empty for a bare seabed. Used
		// by the GUI to colour the voxel overlay by kind: structure cells are rigid solid, the rest of
		// the flow solid mask is erodible sediment. Static (structure never moves).
		const std::vector<unsigned char>& structure_host() const { return structure_host_; }

	private:
		std::vector<unsigned char> build_flow_solid(const std::vector<double>& G_host) const;
		void build_equilibrium_ghost();            // OPEN: fill cin_x* with the Rouse equilibrium profile
		void build_recycle_ghost(const double* u); // RECYCLE: flux-matched outlet→inlet feedback (host)
		void apply_boundary_flux(const double* u, double dt); // open the x-faces per sed_bc + accumulate budget
		void accumulate_bedload_boundary(double dt); // add net inlet−outlet open-x bedload sand to the budget

		MacGrid g_;
		SeabedParams sp_;
		double cpack_ = 0.0, ws0_ = 0.0, Gmax_ = 0.0, remask_thr_ = 0.0;
		int np_ = 0, ncol_ = 0;
		bool has_structure_ = false;
		std::vector<unsigned char> structure_host_; // retained host copy (empty ⇒ bare seabed)
		long long nstep_ = 0;
		bool morph_frozen_ = false; // M6 pre-check: compute τ_b only, don't advance the bed
		bool first_shear_ = true;   // seed the wall-model EMA on the first bedshear call (frozen or not)

		WallParams wp_;
		MorphoParams mp_;
		AvalancheParams ap_;

		// device buffers
		double *G_ = nullptr, *G_fixed_ = nullptr, *G_last_ = nullptr, *Gb_ = nullptr, *dzscr_ = nullptr;
		double *c_ = nullptr, *c2_ = nullptr, *cscr_ = nullptr, *ws_ = nullptr, *Dc_ = nullptr;
		double *taux_ = nullptr, *tauy_ = nullptr, *ustar_ = nullptr, *emax_ = nullptr, *emay_ = nullptr, *tscr_ = nullptr;
		double *fpack_ = nullptr, *Ft_ = nullptr, *nbx_ = nullptr, *nby_ = nullptr, *nbz_ = nullptr, *Ab_ = nullptr;
		double *beta_ = nullptr, *upx_ = nullptr, *upy_ = nullptr, *qx_ = nullptr, *qy_ = nullptr, *divq_ = nullptr;
		double* clip_ = nullptr; // per-column MORFAC-limiter clip flag (1 where |Δz_b| was clamped)
		double* dep_ = nullptr; // per-column applied deposition grain [m] last step (GUI exchange-rate viz)
		double* ero_ = nullptr; // per-column applied erosion grain [m] last step (GUI exchange-rate viz)
		double last_dt_ = 0.0;  // last step's flow dt [s] — for the physical exchange-rate normalisation
		double clip_sum_ = 0.0; long long clip_colsteps_ = 0; // running clip-fraction accumulators
		double structure_discard_ = 0.0; // cumulative sand [m³] removed by the bed clamp+structure pin
		unsigned char *structure_ = nullptr, *frozen_ = nullptr;
		double *shear_mult_ = nullptr; // M6 HSV shear-multiplier field (null ⇒ ×1 everywhere)

		// open-sea / recycle boundary (suspended x-faces): ghost concentration on the yz-plane (ny·nz,
		// indexed k·ny+j) + the per-boundary-cell net-flux scratch (2·ny·nz) and its running budget.
		int ncolf_ = 0; // ny·nz
		double *cin_xmin_ = nullptr, *cin_xmax_ = nullptr, *bflux_ = nullptr;
		double boundary_sand_in_ = 0.0;

		double last_max_ustar_ = 0.0, last_max_tau_ = 0.0, last_zb_min_ = 0.0, last_zb_max_ = 0.0;
	};
}
