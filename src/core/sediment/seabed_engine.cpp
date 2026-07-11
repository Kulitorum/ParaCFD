// seabed_engine.cpp — see seabed_engine.h. Host driver: assembles the gate-verified kernels into
// the live morphodynamic loop. Pure libscour (drives the *_gpu launchers), no Qt/GL.
#include "core/sediment/seabed_engine.h"

#include "core/fluid/bedshear_ops.h" // bedshear_smooth/filter (τ de-noising, reused)
#include "core/fluid/channel_bc.h"   // loglaw_ustar_for_U (open-sea equilibrium u*)
#include "core/fluid/mac_ops.h"      // reductions, scale_add, axpy
#include "core/sediment/seabed_morpho.h"
#include "core/sediment/sed_physics.h"
#include "core/sediment/suspended.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace scour::core
{
	namespace
	{
		double* dalloc0(size_t n) { double* p = nullptr; cudaMalloc(&p, sizeof(double) * n); cudaMemset(p, 0, sizeof(double) * n); return p; }
		double* dfill(size_t n, double v) { double* p = dalloc0(n); std::vector<double> h(n, v); cudaMemcpy(p, h.data(), sizeof(double) * n, cudaMemcpyHostToDevice); return p; }
		unsigned char* ualloc(const std::vector<unsigned char>& h) { unsigned char* p = nullptr; cudaMalloc(&p, h.size()); cudaMemcpy(p, h.data(), h.size(), cudaMemcpyHostToDevice); return p; }
		constexpr double DEG = 3.14159265358979323846 / 180.0;
	}

	SeabedMorpho::SeabedMorpho(MacGrid g, SeabedParams sp, const std::vector<unsigned char>& structure)
		: g_(g), sp_(sp)
	{
		cpack_ = SED_CPACK;
		ws0_ = settling_ws(sp_.d50, sp_.rho, sp_.rho_s, sp_.nu);
		np_ = g_.p_count();
		ncol_ = g_.nx * g_.ny;
		Gmax_ = std::max(0.0, (g_.nz * g_.h - 2.0 * g_.h)) * cpack_;
		remask_thr_ = sp_.remask_frac * g_.h;
		has_structure_ = false;
		if ((int)structure.size() == np_)
			for (unsigned char s : structure) if (s) { has_structure_ = true; break; }
		if (has_structure_) structure_host_ = structure; // retain for the GUI overlay kind-colouring

		// Wall model (grain-skin τ_b). AUTO regime picks rough/transitional by ks+ (d50=0.2mm ⇒
		// transitional, CLAUDE.md M3 note). smooth + EMA de-noise τ (research/04 anti-stair-step).
		wp_.kappa = 0.40; wp_.d50 = sp_.d50; wp_.nu = sp_.nu; wp_.rho = sp_.rho;
		wp_.regime = WALL_AUTO; wp_.cj_iters = 5; wp_.t_avg = 2.0; wp_.smooth = 1;

		mp_.d50 = sp_.d50; mp_.rho = sp_.rho; mp_.rho_s = sp_.rho_s; mp_.nu = sp_.nu; mp_.ws0 = ws0_;
		mp_.cpack = cpack_; mp_.morfac = sp_.morfac; mp_.dz_limit_frac = 0.05; mp_.hindered = 1;
		mp_.erosion_coeff = sp_.erosion_coeff; // Winterwerp coeff (M6 erosion-rate knob; default 0.018)
		mp_.bedload_formula = sp_.bedload_formula;
		mp_.erosion_on = 1; mp_.deposition_on = 1; mp_.bedload_on = 1;

		ap_.tan_trigger = cpack_ * std::tan(32.0 * DEG); // avalanche acts on G = z_b·c_pack
		ap_.tan_repose = cpack_ * std::tan(30.0 * DEG);
		ap_.relax = 0.1; ap_.periodic = 0;

		double G0 = sp_.sand_depth * cpack_;
		G_ = dfill(ncol_, G0); G_fixed_ = dfill(ncol_, G0); G_last_ = dfill(ncol_, G0);
		Gb_ = dalloc0(ncol_); dzscr_ = dalloc0(ncol_);
		c_ = dalloc0(np_); c2_ = dalloc0(np_); cscr_ = dalloc0(np_); ws_ = dfill(np_, ws0_); Dc_ = dalloc0(np_);
		taux_ = dalloc0(ncol_); tauy_ = dalloc0(ncol_); ustar_ = dalloc0(ncol_);
		emax_ = dalloc0(ncol_); emay_ = dalloc0(ncol_); tscr_ = dalloc0(ncol_);
		fpack_ = dalloc0(np_); Ft_ = dalloc0(np_); nbx_ = dalloc0(np_); nby_ = dalloc0(np_); nbz_ = dalloc0(np_); Ab_ = dalloc0(np_);
		beta_ = dalloc0(ncol_); upx_ = dfill(ncol_, 1.0); upy_ = dalloc0(ncol_);
		qx_ = dalloc0(ncol_); qy_ = dalloc0(ncol_); divq_ = dalloc0(ncol_); clip_ = dalloc0(ncol_);
		dep_ = dalloc0(ncol_); ero_ = dalloc0(ncol_); // GUI exchange-rate diagnostic (applied dep/pickup)
		ncolf_ = g_.ny * g_.nz; // yz boundary-plane cell count (open-sea / recycle x-face ghosts)
		cin_xmin_ = dalloc0(ncolf_); cin_xmax_ = dalloc0(ncolf_); bflux_ = dalloc0(2 * ncolf_);

		// Structure solid mask + its per-column footprint (frozen columns: rigid, never erode/accrete).
		std::vector<unsigned char> frozen(ncol_, 0);
		if (has_structure_)
		{
			structure_ = ualloc(structure);
			for (int k = 0; k < g_.nz; ++k) for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i)
				if (structure[g_.pidx(i, j, k)]) frozen[j * g_.nx + i] = 1;
		}
		frozen_ = ualloc(frozen);
		last_zb_min_ = last_zb_max_ = sp_.sand_depth;
		build_equilibrium_ghost(); // seed the open-sea inflow profile (used only when sed_bc == OPEN)
	}

	SeabedMorpho::~SeabedMorpho()
	{
		for (double* p : { G_, G_fixed_, G_last_, Gb_, dzscr_, c_, c2_, cscr_, ws_, Dc_,
			taux_, tauy_, ustar_, emax_, emay_, tscr_, fpack_, Ft_, nbx_, nby_, nbz_, Ab_,
			beta_, upx_, upy_, qx_, qy_, divq_, clip_, dep_, ero_, cin_xmin_, cin_xmax_, bflux_ })
			if (p) cudaFree(p);
		if (structure_) cudaFree(structure_);
		if (frozen_) cudaFree(frozen_);
		if (shear_mult_) cudaFree(shear_mult_);
	}

	void SeabedMorpho::set_shear_multiplier(const std::vector<double>& mult)
	{
		if ((int)mult.size() != ncol_) return; // ignore a mismatched field
		if (!shear_mult_) shear_mult_ = dalloc0(ncol_);
		cudaMemcpy(shear_mult_, mult.data(), sizeof(double) * ncol_, cudaMemcpyHostToDevice);
	}

	std::vector<unsigned char> SeabedMorpho::set_structure(const std::vector<unsigned char>& structure)
	{
		if ((int)structure.size() != np_) return {}; // ignore a mismatched mask (keep the current one)

		bool any = false;
		for (unsigned char s : structure) if (s) { any = true; break; }
		has_structure_ = any;
		structure_host_ = any ? structure : std::vector<unsigned char>{};

		// (Re)upload the device structure mask (allocate on first use).
		if (!structure_) structure_ = ualloc(structure);
		else cudaMemcpy(structure_, structure.data(), (size_t)np_, cudaMemcpyHostToDevice);

		// Recompute the per-column frozen footprint (rigid columns never erode/accrete; a column a block
		// has left is no longer frozen, so its sand resumes morphology).
		std::vector<unsigned char> frozen((size_t)ncol_, 0);
		if (any)
			for (int k = 0; k < g_.nz; ++k) for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i)
				if (structure[g_.pidx(i, j, k)]) frozen[j * g_.nx + i] = 1;
		cudaMemcpy(frozen_, frozen.data(), (size_t)ncol_, cudaMemcpyHostToDevice);

		// New flow solid mask (sand ∪ new structure) for the caller to re-mask the flow (update_solid).
		std::vector<double> G_host((size_t)ncol_);
		cudaMemcpy(G_host.data(), G_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
		return build_flow_solid(G_host);
	}

	std::vector<unsigned char> SeabedMorpho::build_flow_solid(const std::vector<double>& G_host) const
	{
		std::vector<unsigned char> solid((size_t)np_, 0);
		for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i)
		{
			double zb = G_host[j * g_.nx + i] / cpack_;
			int ktop = (int)std::floor(zb / g_.h + 0.5); // cells with centre below z_b ⇒ sand
			if (ktop < 0) ktop = 0; if (ktop > g_.nz) ktop = g_.nz;
			for (int k = 0; k < ktop; ++k) solid[g_.pidx(i, j, k)] = 1;
		}
		if (has_structure_)
		{
			std::vector<unsigned char> str((size_t)np_);
			cudaMemcpy(str.data(), structure_, (size_t)np_, cudaMemcpyDeviceToHost);
			for (int n = 0; n < np_; ++n) if (str[n]) solid[n] = 1;
		}
		return solid;
	}

	std::vector<unsigned char> SeabedMorpho::initial_flow_solid() const
	{
		std::vector<double> G_host((size_t)ncol_, sp_.sand_depth * cpack_);
		return build_flow_solid(G_host);
	}

	bool SeabedMorpho::step(const double* u, const double* v, const double* w, const double* nut,
		double dt, std::vector<unsigned char>* new_solid)
	{
		// (0) near-bed grain-skin τ_b from the LIVE flow, de-noised (3×3 smooth + rate-limit/EMA).
		seabed_bedshear_gpu(u, v, G_, taux_, tauy_, ustar_, g_, wp_, cpack_);
		if (shear_mult_) seabed_apply_shear_mult_gpu(taux_, tauy_, ustar_, shear_mult_, g_); // M6 HSV amplification (before smooth/EMA)
		if (wp_.smooth > 0) bedshear_smooth_gpu(taux_, tauy_, tscr_, g_, wp_.smooth);
		bedshear_filter_gpu(taux_, tauy_, emax_, emay_, g_, wp_, dt, first_shear_ ? 1 : 0);
		first_shear_ = false;
		last_max_ustar_ = reduce_max_abs_gpu(ustar_, ncol_);
		last_max_tau_ = sp_.rho * last_max_ustar_ * last_max_ustar_;

		last_dt_ = dt; // for the physical exchange-rate normalisation (copy_exchange_host)

		// M6 pre-check: bed frozen ⇒ τ_b/u* are refreshed (above) but the bed does not move and no
		// flow-mask rebuild is triggered. Suspended/Exner/avalanche are skipped; nstep_ (the
		// morphological-time clock) does not advance.
		if (morph_frozen_)
		{
			cudaMemset(dep_, 0, sizeof(double) * ncol_); // no exchange while frozen ⇒ zero rate
			cudaMemset(ero_, 0, sizeof(double) * ncol_);
			return false;
		}

		// (1) suspended sediment: hindered settling, conservative settling advection, ν_t/σ_s diffusion.
		suspended_effective_ws_gpu(c_, ws_, g_, ws0_, mp_.hindered);
		suspended_advect_cons_gpu(c_, u, v, w, ws_, c2_, cscr_, g_, dt);
		std::swap(c_, c2_);
		apply_boundary_flux(u, dt); // open-sea inflow / free outflow (or recycle); a no-op when CLOSED
		if (sp_.diffusion_on && nut)
		{
			scale_add_gpu(Dc_, 1.0 / sp_.sigma_s, nut, 0.0, np_); // Dc = ν_t/σ_s
			double Dcmax = reduce_max_abs_gpu(Dc_, np_);
			int nsub = 1;
			if (Dcmax > 0.0) nsub = std::max(1, (int)std::ceil(6.0 * dt * Dcmax / (g_.h * g_.h)));
			nsub = std::min(nsub, 32); // safety cap; if it wants more, the flow dt is already tiny
			ScalarBC sbc; // all-Neumann (zero-flux) — a closed suspended domain (mass conserved)
			double dts = dt / nsub;
			for (int s = 0; s < nsub; ++s) { suspended_diffuse_gpu(c_, c2_, Dc_, g_, sbc, dts); std::swap(c_, c2_); }
		}

		// (2) bed geometry + bedload flux (from the de-noised grain-skin τ = ema).
		fpack_from_G_gpu(G_, fpack_, g_, cpack_);
		bed_boxfilter_gpu(fpack_, nullptr, Ft_, g_, cpack_);
		bed_interface_geom_gpu(Ft_, nbx_, nby_, nbz_, Ab_, g_);
		bed_column_geom_gpu(fpack_, nbx_, nby_, nbz_, beta_, upx_, upy_, g_, cpack_);
		bedload_flux_gpu(emax_, emay_, beta_, upx_, upy_, qx_, qy_, g_, mp_);
		int bl_xopen = (sp_.sed_bc != SED_BC_CLOSED && mp_.bedload_on) ? 1 : 0; // open the streamwise bedload faces
		bedload_div_gpu(qx_, qy_, emax_, emay_, divq_, g_, /*periodic*/ 0, bl_xopen);
		if (bl_xopen) accumulate_bedload_boundary(dt); // net inlet−outlet bedload sand → open budget

		// (3) Exner column update (deposition + Winterwerp erosion + bedload div, ×MORFAC, limiter);
		//     exchanges the eroded/deposited grain with c at k_bed (mass-conserving). Then (4) avalanche.
		cudaMemset(clip_, 0, sizeof(double) * ncol_); // morpho_exner only SETS clip flags; clear first
		morpho_exner_gpu(G_, c_, emax_, emay_, beta_, upx_, upy_, divq_, clip_, g_, mp_, dt, dep_, ero_);
		clip_sum_ += reduce_sum_gpu(clip_, ncol_); clip_colsteps_ += ncol_; // MORFAC-limiter hit-rate diagnostic (M9)
		avalanche_sweep_gpu(G_, Gb_, g_, ap_); std::swap(G_, Gb_);

		// (5) keep suspended sand in the fluid column; clamp G ∈ [0,Gmax] and pin the rigid footprint.
		//     The clamp+pin can discard deposition landing on a rigid footprint (a physical sink); track
		//     the removed grain [m³] so an open-budget audit can credit it exactly (structure_discard()).
		seabed_confine_c_gpu(c_, G_, structure_, g_, cpack_);
		double bed_pre = reduce_sum_gpu(G_, ncol_);
		seabed_bed_post_gpu(G_, G_fixed_, frozen_, g_, Gmax_);
		structure_discard_ += (bed_pre - reduce_sum_gpu(G_, ncol_)) * g_.h * g_.h;
		++nstep_;

		// (6) flow-mask rebuild when the bed has moved > threshold since the last build.
		cudaMemcpy(dzscr_, G_, sizeof(double) * ncol_, cudaMemcpyDeviceToDevice);
		axpy_gpu(-1.0, G_last_, dzscr_, ncol_); // dzscr = G − G_last
		double dzmax = reduce_max_abs_gpu(dzscr_, ncol_) / cpack_;
		if (dzmax > remask_thr_ && new_solid)
		{
			std::vector<double> G_host((size_t)ncol_);
			cudaMemcpy(G_host.data(), G_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
			*new_solid = build_flow_solid(G_host);
			cudaMemcpy(G_last_, G_, sizeof(double) * ncol_, cudaMemcpyDeviceToDevice);
			return true;
		}
		return false;
	}

	void SeabedMorpho::copy_zb_host(std::vector<float>& zb) const
	{
		std::vector<double> G_host((size_t)ncol_);
		cudaMemcpy(G_host.data(), G_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
		zb.resize((size_t)ncol_);
		double mn = 1e30, mx = -1e30;
		for (int n = 0; n < ncol_; ++n) { double z = G_host[n] / cpack_; zb[n] = (float)z; mn = std::min(mn, z); mx = std::max(mx, z); }
		const_cast<SeabedMorpho*>(this)->last_zb_min_ = mn;
		const_cast<SeabedMorpho*>(this)->last_zb_max_ = mx;
	}

	void SeabedMorpho::copy_ustar_host(std::vector<float>& us) const
	{
		std::vector<double> h((size_t)ncol_);
		cudaMemcpy(h.data(), ustar_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
		us.resize((size_t)ncol_);
		for (int n = 0; n < ncol_; ++n) us[n] = (float)h[n];
	}

	void SeabedMorpho::copy_exchange_host(std::vector<float>& rate) const
	{
		rate.resize((size_t)ncol_);
		// Physical (real-time) bed-elevation velocity from the applied suspended exchange this step:
		// (dep − ero grain)/(c_pack·M·dt). Dividing out MORFAC makes it MORFAC-invariant — the pattern +
		// sign (the diagnostic) are identical either way, so we report the honest physical rate. Zero
		// before the first step (last_dt_==0) or with no morfac (guarded), so the bed shows a flat 0.
		double denom = cpack_ * std::max(mp_.morfac, 1e-12) * last_dt_;
		if (denom <= 0.0) { std::fill(rate.begin(), rate.end(), 0.0f); return; }
		std::vector<double> dep((size_t)ncol_), ero((size_t)ncol_);
		cudaMemcpy(dep.data(), dep_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
		cudaMemcpy(ero.data(), ero_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
		double inv = 1.0 / denom;
		for (int n = 0; n < ncol_; ++n) rate[n] = (float)((dep[n] - ero[n]) * inv);
	}

	double SeabedMorpho::bed_volume() const { return reduce_sum_gpu(G_, ncol_) * g_.h * g_.h; }
	double SeabedMorpho::susp_volume() const { return reduce_sum_gpu(c_, np_) * g_.h * g_.h * g_.h; }

	void SeabedMorpho::save_state(std::vector<double>& G, std::vector<double>& c,
		std::vector<double>& emax, std::vector<double>& emay, long long& nstep) const
	{
		G.resize((size_t)ncol_); c.resize((size_t)np_); emax.resize((size_t)ncol_); emay.resize((size_t)ncol_);
		cudaMemcpy(G.data(), G_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
		cudaMemcpy(c.data(), c_, sizeof(double) * np_, cudaMemcpyDeviceToHost);
		cudaMemcpy(emax.data(), emax_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
		cudaMemcpy(emay.data(), emay_, sizeof(double) * ncol_, cudaMemcpyDeviceToHost);
		nstep = nstep_;
	}

	void SeabedMorpho::load_state(const std::vector<double>& G, const std::vector<double>& c,
		const std::vector<double>& emax, const std::vector<double>& emay, long long nstep)
	{
		if ((int)G.size() == ncol_) cudaMemcpy(G_, G.data(), sizeof(double) * ncol_, cudaMemcpyHostToDevice);
		if ((int)c.size() == np_) cudaMemcpy(c_, c.data(), sizeof(double) * np_, cudaMemcpyHostToDevice);
		if ((int)emax.size() == ncol_) cudaMemcpy(emax_, emax.data(), sizeof(double) * ncol_, cudaMemcpyHostToDevice);
		if ((int)emay.size() == ncol_) cudaMemcpy(emay_, emay.data(), sizeof(double) * ncol_, cudaMemcpyHostToDevice);
		cudaMemcpy(G_last_, G_, sizeof(double) * ncol_, cudaMemcpyDeviceToDevice); // remask baseline = restored bed
		nstep_ = nstep;
		cudaDeviceSynchronize();
	}

	void SeabedMorpho::set_boundary_mode(int mode)
	{
		sp_.sed_bc = mode;
		if (mode == SED_BC_OPEN) build_equilibrium_ghost(); // ensure the ghost matches the current U
	}

	void SeabedMorpho::set_inlet_speed(double U)
	{
		sp_.U_inlet = U;
		build_equilibrium_ghost(); // re-derive u*, c_a and the Rouse profile for the new current
	}

	// OPEN: the far-field inlet carries the EQUILIBRIUM suspended load. Near-bed c_a is this model's own
	// pickup=deposition fixed point (equilibrium_cb), decaying upward as the Rouse profile. u* is the
	// log-law friction velocity flux-matched to U over the fluid depth (identical to the seabed log-law
	// inlet BL). The far field is x/y-homogeneous ⇒ a function of k only, broadcast over j into both
	// x-face ghosts (tidal reversal turns the outlet into an inlet using the same equilibrium load).
	void SeabedMorpho::build_equilibrium_ghost()
	{
		if (!cin_xmin_ || !cin_xmax_) return;
		double H = g_.nz * g_.h - sp_.sand_depth; // fluid depth above the flat far-field bed top
		if (H < g_.h) H = g_.h;
		double z0 = (sp_.d50 > 0.0) ? sp_.d50 / 12.0 : 1e-4;
		double kappa = 0.40;
		double ustar = loglaw_ustar_for_U(sp_.U_inlet, H, z0, kappa);
		double tau_eq = sp_.rho * ustar * ustar; // far-field grain-skin bed shear ρ·u*²
		double c_a = equilibrium_cb(tau_eq, sp_.d50, sp_.rho, sp_.rho_s, sp_.nu, sp_.alpha, ws0_);
		double R = rouse_number(ws0_, sp_.sigma_s, kappa, ustar);
		double aref = 0.5 * g_.h; // reference height above the bed = first fluid cell centre
		std::vector<double> ghost((size_t)ncolf_, 0.0);
		for (int k = 0; k < g_.nz; ++k)
		{
			double zeta = (k + 0.5) * g_.h - sp_.sand_depth; // height above the bed top
			double val = zeta <= 0.0 ? 0.0 : rouse_profile(zeta, aref, H, c_a, R);
			for (int j = 0; j < g_.ny; ++j) ghost[(size_t)k * g_.ny + j] = val; // t = k·ny + j
		}
		cudaMemcpy(cin_xmin_, ghost.data(), sizeof(double) * ncolf_, cudaMemcpyHostToDevice);
		cudaMemcpy(cin_xmax_, ghost.data(), sizeof(double) * ncolf_, cudaMemcpyHostToDevice);
	}

	// RECYCLE: feed the outlet-plane load back to the inlet, scaled by ONE factor so the injected inflow
	// mass exactly equals the extracted outflow mass — an exactly-conserving recirculating flume. Done
	// host-side (recycle is the non-default cross-check): the x-planes are tiny (ny·nz) strided copies.
	// The scale is derived from the SAME c/u the kernel will read, so Σ net_flux ≈ 0 to machine precision.
	void SeabedMorpho::build_recycle_ghost(const double* u)
	{
		if (!cin_xmin_ || !cin_xmax_) return;
		int nf = ncolf_;
		size_t db = sizeof(double);
		std::vector<double> cin(nf), cout(nf), uin(nf), uout(nf);
		// c planes (pidx stride = nx): inlet i=0, outlet i=nx−1. u planes (uidx stride = nx+1): i=0, i=nx.
		cudaMemcpy2D(cin.data(), db, c_, g_.nx * db, db, nf, cudaMemcpyDeviceToHost);
		cudaMemcpy2D(cout.data(), db, c_ + (g_.nx - 1), g_.nx * db, db, nf, cudaMemcpyDeviceToHost);
		cudaMemcpy2D(uin.data(), db, u, (g_.nx + 1) * db, db, nf, cudaMemcpyDeviceToHost);
		cudaMemcpy2D(uout.data(), db, u + g_.nx, (g_.nx + 1) * db, db, nf, cudaMemcpyDeviceToHost);
		double Mout = 0.0, Min_raw = 0.0;
		for (int t = 0; t < nf; ++t)
		{
			double ui = uin[t], uo = uout[t];
			if (ui >= 0.0) Min_raw += ui * cout[t]; else Mout += (-ui) * cin[t]; // xmin: inflow recycles outlet
			if (uo >= 0.0) Mout += uo * cout[t]; else Min_raw += (-uo) * cin[t]; // xmax: inflow recycles inlet
		}
		double scale = Min_raw > 1e-300 ? Mout / Min_raw : 0.0;
		std::vector<double> gmin(nf), gmax(nf);
		for (int t = 0; t < nf; ++t) { gmin[t] = scale * cout[t]; gmax[t] = scale * cin[t]; }
		cudaMemcpy(cin_xmin_, gmin.data(), db * nf, cudaMemcpyHostToDevice);
		cudaMemcpy(cin_xmax_, gmax.data(), db * nf, cudaMemcpyHostToDevice);
	}

	void SeabedMorpho::apply_boundary_flux(const double* u, double dt)
	{
		if (sp_.sed_bc == SED_BC_CLOSED) return; // zero-flux box (M4/M5, seabed_flat) — byte-identical
		if (sp_.sed_bc == SED_BC_RECYCLE) build_recycle_ghost(u);
		// OPEN: cin_xmin_/cin_xmax_ already hold the equilibrium Rouse profile (ctor / set_inlet_speed).
		suspended_boundary_flux_gpu(c_, u, cin_xmin_, cin_xmax_, bflux_, g_, dt);
		boundary_sand_in_ += reduce_sum_gpu(bflux_, 2 * ncolf_);
	}

	// Net bedload sand [m³] crossing the open x-boundary this step = M·dt·h·Σ_j (qx_inlet − qx_outlet):
	// the zero-gradient faces feed q at the inlet and export q at the outlet, so the domain gains the
	// difference. Extracts the two strided x-planes of qx_ (col = j·nx + i, stride nx) to host. Limiter-
	// ACCURATE (not machine-exact): the Exner |Δz_b| clamp can trim the applied change at a boundary
	// column (rare on a developed bed), so the open budget closes to that accuracy on the bedload path.
	void SeabedMorpho::accumulate_bedload_boundary(double dt)
	{
		size_t db = sizeof(double);
		std::vector<double> qin(g_.ny), qout(g_.ny);
		cudaMemcpy2D(qin.data(), db, qx_, g_.nx * db, db, g_.ny, cudaMemcpyDeviceToHost);
		cudaMemcpy2D(qout.data(), db, qx_ + (g_.nx - 1), g_.nx * db, db, g_.ny, cudaMemcpyDeviceToHost);
		double si = 0.0, so = 0.0;
		for (int j = 0; j < g_.ny; ++j) { si += qin[j]; so += qout[j]; }
		boundary_sand_in_ += mp_.morfac * dt * g_.h * (si - so);
	}
}
