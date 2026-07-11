// channel_core.h — the M2 open-channel fluid core: MAC grid + masked MacCormack
// advection + explicit Smagorinsky + masked/Dirichlet MGPCG projection, with a rough
// log-law/uniform inlet, Orlanski-type outlet + global flux rescale, free-slip lid &
// laterals, and a voxel obstacle mask (free-slip or no-slip surface). Implements the
// FluidCore interface (PLAN §3). Per-step pipeline follows RESEARCH §3 / research/08 §.
#pragma once

#include "core/fluid/channel_bc.h"
#include "core/fluid/channel_pressure.h"
#include "core/fluid/fluid_core.h"
#include "core/fluid/inlet_fluct.h"

#include <memory>
#include <vector>

namespace windcfd::core
{
	struct ChannelParams
	{
		double rho = 1.0;    // kg/m^3 (nondimensional cylinder uses 1)
		double nu = 0.01;    // m^2/s molecular (or artificial ν = U·D/Re) — explicit
		double Cs = 0.0;     // Smagorinsky constant (0 = LES off)
		double cfl = 1.0;    // advective CFL target (adaptive dt)
		double safety = 0.9; // dt safety factor
		double proj_tol = 1e-4;
		int proj_max_iter = 60;
		double fixed_dt = 0.0; // >0 => fixed dt (uniform sampling for the St FFT)
		int advect_band = 1;   // nearsolid dilation band (1st-order reversion radius)
	};

	class ChannelFluidCore final : public FluidCore
	{
	public:
		// `solid` is the finest-level cell mask (host); may be all-zero (empty channel).
		ChannelFluidCore(MacGrid grid, ChannelBC bc, ChannelParams pr, const std::vector<unsigned char>& solid);
		~ChannelFluidCore() override;
		ChannelFluidCore(const ChannelFluidCore&) = delete;
		ChannelFluidCore& operator=(const ChannelFluidCore&) = delete;

		double step() override;
		FluidSnapshot snapshot() const override;
		const MacGrid& grid() const override { return g_; }

		// Initialise the whole fluid interior to a uniform u (v=w=0) plus the inlet;
		// optional transverse blip on v to break shedding symmetry.
		void init_uniform(double u0, double v_blip = 0.0);
		// Initialise u to the inlet profile everywhere (steady open-channel start).
		void init_inlet_profile();

		// In-place update of the interior solid mask, PRESERVING the current flow (u,v,w,p): recompute
		// the near-solid band, upload the mask, rebuild the pressure-solver level masks, and zero the
		// now-solid faces. Used by the moving erodible bed (the seabed scenario) to re-mask the flow as
		// z_b changes WITHOUT resetting the developing wake (unlike a factory rebuild). No physics/gate
		// impact: step() is unchanged, and the M2 gate never calls this.
		void update_solid(const std::vector<unsigned char>& solid);

		// Optional inlet-turbulence generator (M3 SEM / precursor). Null ⇒ M2 behaviour
		// unchanged (byte-identical). RESEARCH §8, research/14.
		void set_inlet_fluct(InletFluct* f) { fluct_ = f; }

		// Sub-grid porous/thin-screen resistance field (per-cell k = 1/β²−1; M7 perforated units,
		// RESEARCH §8). Sizes to p_count; k=0 everywhere ⇒ no drag (never allocated ⇒ byte-identical).
		// Applied as an implicit Forchheimer sink after diffusion, before projection (channel_porous.h).
		void set_porous(const std::vector<double>& k_cell);
		bool has_porous() const { return porous_k_ != nullptr; }

		// Open-channel-with-erodible-bed option: when a solid bed reaches the inlet plane (the seabed
		// scenario), the Dirichlet inlet must NOT drive the buried (solid) cells at i=0 — otherwise it
		// injects unremovable inflow into the bed, seeding a corner anomaly that destabilises the wake.
		// Enabling this zeros the i=0 inlet u-faces whose cell is solid. Default OFF ⇒ M2/M3 unchanged
		// (their obstacles never touch the inlet plane, so it is a no-op there anyway).
		void set_bed_inlet_mask(bool on) { bed_inlet_mask_ = on; }

		// Live-update the inlet + convective-outlet speed WITHOUT resetting the developing flow — the
		// running wake evolves toward the new current. The GUI "Input speed" control uses this to test
		// different currents interactively. In log-law mode u* is re-derived so the profile stays
		// flux-matched to U. The M2/M3 gates never call it, so no gate impact. RESEARCH §8.1.
		void set_inlet_speed(double U)
		{
			bc_.U_inlet = U; bc_.Uc = U;
			if (bc_.inlet_mode == INLET_LOGLAW) recompute_loglaw_ustar();
		}
		double inlet_speed() const { return bc_.U_inlet; }

		// Tidal reversal (M9, RESEARCH §8): flip the streamwise flow direction WITHOUT resetting the
		// flow. sign≥0 ⇒ inlet on xmin / Orlanski outlet on xmax (M2/M3 default); sign<0 ⇒ the faces
		// swap (inlet drives −x on xmax, outlet on xmin). U_inlet/ustar stay MAGNITUDES; the driver sets
		// the speed with set_inlet_speed(|U|) and flips this at slack (|U|→0) so the swap is smooth. Also
		// re-points the pressure solver's Dirichlet p=0 pin to the new outlet face. Never called by the
		// M2/M3 gates ⇒ byte-identical there. Pair with a live re-hold of the inlet on the next step().
		void set_flow_direction(int sign)
		{
			bc_.flow_sign = sign < 0 ? -1 : 1;
			solver_->set_outlet_dir(bc_.flow_sign);
		}
		int flow_direction() const { return bc_.flow_sign; }

		// Live-switch the inlet between UNIFORM (top-hat) and a LOG-LAW boundary-layer profile referenced
		// to the bed top (z0 = grain roughness, bed_datum = inlet bed elevation). The BL profile makes the
		// near-bed inlet velocity → 0 at the sand surface, which removes the top-hat's spurious full-U
		// leading-edge bed shear (over-erosion of the sediment's upstream edge). u* is flux-matched to the
		// current U over the fluid depth above the bed. RESEARCH §8.1. GUI "Inlet profile" toggle.
		void set_inlet_profile(bool loglaw, double z0, double bed_datum)
		{
			if (loglaw)
			{
				bc_.inlet_mode = INLET_LOGLAW;
				if (z0 > 0.0) bc_.z0 = z0;
				bc_.bed_datum = bed_datum;
				if (bc_.kappa <= 0.0) bc_.kappa = 0.40;
				recompute_loglaw_ustar();
			}
			else
				bc_.inlet_mode = INLET_UNIFORM;
		}
		bool inlet_is_loglaw() const { return bc_.inlet_mode == INLET_LOGLAW; }

		int last_solve_iters() const { return last_iters_; }
		double last_solve_relres() const { return last_relres_; }
		double last_dt() const { return last_dt_; }
		double inlet_flux() const { return q_in_; }
		double outlet_flux() const { return q_out_; }

		// --- Checkpoint save / restore (scene_io) --------------------------------------------
		// The live boundary descriptor (carries the interactive inlet speed/profile overrides) and
		// whether the bed-inlet mask is active — both part of the restorable state.
		ChannelBC bc() const { return bc_; }
		bool bed_inlet_mask_on() const { return bed_inlet_mask_; }
		// Copy the stateful velocity + pressure fields to host (face-staggered doubles, MAC counts).
		// Scratch (uB/vB/wB, ν_t, near-solid) is recomputed each step and NOT part of a checkpoint.
		void copy_state_host(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p) const;
		// Overwrite the velocity + pressure fields from host arrays (H2D), preserving the solid mask
		// the core was constructed with. Sizes must match the MAC counts, else the field is skipped.
		// Used by a scene restore right after make_core() (replaces the init_uniform state).
		void load_state_host(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& p);

		// Device accessors (const) for probes / drag integration in the driver.
		const double* u_dev() const { return uA_; }
		const double* v_dev() const { return vA_; }
		const double* w_dev() const { return wA_; }
		const double* p_dev() const { return pres_; }
		const double* nut_dev() const { return nut_; } // Smagorinsky ν_t (valid after step() when Cs>0)
		const unsigned char* solid_dev() const { return solid_; }

		// Single-value host reads (small copies).
		double get_u(int i, int j, int k) const;
		double get_v(int i, int j, int k) const;

	private:
		double adaptive_dt();
		void apply_bcs();
		void mask_inlet_solid(); // zero i=0 inlet u-faces into solid cells (bed_inlet_mask_ only)
		bool bed_inlet_mask_ = false;
		// Re-derive the log-law inlet u* so the profile stays flux-matched to U_inlet over the fluid depth
		// H = domain height − bed_datum. Called on a speed change or a profile switch (loglaw only).
		void recompute_loglaw_ustar()
		{
			double H = g_.nz * g_.h - bc_.bed_datum;
			bc_.ustar = loglaw_ustar_for_U(bc_.U_inlet, H, bc_.z0, bc_.kappa);
		}

		MacGrid g_;
		ChannelBC bc_;
		ChannelParams pr_;
		std::unique_ptr<ChannelMgpcg> solver_;
		InletFluct* fluct_ = nullptr; // M3 inlet turbulence (null ⇒ M2 behaviour)

		double *uA_, *vA_, *wA_;
		double *uB_, *vB_, *wB_;
		double *nut_, *scratch_, *planescr_;
		double *rhs_, *pres_, *divscr_;
		unsigned char *solid_, *nearsolid_;
		double* porous_k_ = nullptr; // sub-grid porous resistance k per cell (null ⇒ no screen)
		int last_iters_ = 0;
		double last_relres_ = 0.0, last_dt_ = 0.0, q_in_ = 0.0, q_out_ = 0.0;
	};
}
