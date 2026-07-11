// channel_core.cu — see channel_core.h. RESEARCH §3, §8, research/08, research/14.
#include "core/fluid/channel_core.h"
#include "core/fluid/channel_ops.h"
#include "core/fluid/channel_porous.h" // sub-grid porous momentum sink (M7)
#include "core/fluid/mac_ops.h" // reductions

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace windcfd::core
{
	namespace
	{
		double* dalloc0(int n) { double* p = nullptr; cudaMalloc(&p, sizeof(double) * (size_t)n); cudaMemset(p, 0, sizeof(double) * (size_t)n); return p; }
		unsigned char* ualloc0(int n) { unsigned char* p = nullptr; cudaMalloc(&p, (size_t)n); cudaMemset(p, 0, (size_t)n); return p; }
		int maxface(MacGrid g) { return std::max(g.u_count(), std::max(g.v_count(), g.w_count())); }

		struct Bbox { bool any = false; int i0 = 0, i1 = 0, j0 = 0, j1 = 0; };
		Bbox solid_bbox(MacGrid g, const std::vector<unsigned char>& s)
		{
			Bbox b; b.i0 = g.nx; b.j0 = g.ny;
			for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
				if (s[g.pidx(i, j, k)]) { b.any = true; b.i0 = std::min(b.i0, i); b.i1 = std::max(b.i1, i); b.j0 = std::min(b.j0, j); b.j1 = std::max(b.j1, j); }
			return b;
		}
	}

	static Bbox g_bbox; // per last-constructed core (single-core usage in the M2 gate)

	namespace
	{
		// Zero the inlet (i=0) u-faces whose cell is solid — the erodible bed reaches the inlet plane
		// in the seabed scenario, and the Dirichlet inlet must not drive those buried cells.
		__global__ void k_zero_inlet_solid(double* u, const unsigned char* solid, MacGrid g, int njk, int icell, int iplane)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= njk) return;
			int j = t % g.ny, k = t / g.ny;
			if (solid[g.pidx(icell, j, k)]) u[g.uidx(iplane, j, k)] = 0.0;
		}
	}

	void ChannelFluidCore::mask_inlet_solid()
	{
		int njk = g_.ny * g_.nz;
		int iplane = bc_.flow_sign < 0 ? g_.nx : 0;      // reversed tide: inlet face is xmax
		int icell = bc_.flow_sign < 0 ? g_.nx - 1 : 0;   // …and its bed cell is the last x-cell
		k_zero_inlet_solid<<<(njk + 255) / 256, 256>>>(uA_, solid_, g_, njk, icell, iplane);
	}

	ChannelFluidCore::ChannelFluidCore(MacGrid grid, ChannelBC bc, ChannelParams pr, const std::vector<unsigned char>& solid)
		: g_(grid), bc_(bc), pr_(pr)
	{
		solver_ = std::make_unique<ChannelMgpcg>(g_);
		uA_ = dalloc0(g_.u_count()); vA_ = dalloc0(g_.v_count()); wA_ = dalloc0(g_.w_count());
		uB_ = dalloc0(g_.u_count()); vB_ = dalloc0(g_.v_count()); wB_ = dalloc0(g_.w_count());
		nut_ = dalloc0(g_.p_count()); scratch_ = dalloc0(maxface(g_)); planescr_ = dalloc0(g_.ny * g_.nz);
		rhs_ = dalloc0(g_.p_count()); pres_ = dalloc0(g_.p_count()); divscr_ = dalloc0(g_.p_count());
		solid_ = ualloc0(g_.p_count()); nearsolid_ = ualloc0(g_.p_count());

		std::vector<unsigned char> s = solid;
		if ((int)s.size() != g_.p_count()) s.assign((size_t)g_.p_count(), 0);
		std::vector<unsigned char> ns;
		// nearsolid dilation band = advect_band.
		{
			ns.assign((size_t)g_.p_count(), 0);
			int band = pr_.advect_band;
			for (int k = 0; k < g_.nz; ++k) for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i)
			{
				bool near = false;
				for (int dk = -band; dk <= band && !near; ++dk) for (int dj = -band; dj <= band && !near; ++dj) for (int di = -band; di <= band && !near; ++di)
				{
					int ii = i + di, jj = j + dj, kk = k + dk;
					if (ii < 0 || ii >= g_.nx || jj < 0 || jj >= g_.ny || kk < 0 || kk >= g_.nz) continue;
					if (s[g_.pidx(ii, jj, kk)]) near = true;
				}
				ns[g_.pidx(i, j, k)] = near ? 1 : 0;
			}
		}
		cudaMemcpy(solid_, s.data(), (size_t)g_.p_count(), cudaMemcpyHostToDevice);
		cudaMemcpy(nearsolid_, ns.data(), (size_t)g_.p_count(), cudaMemcpyHostToDevice);
		solver_->build_masks(solid_);
		g_bbox = solid_bbox(g_, s);
	}

	ChannelFluidCore::~ChannelFluidCore()
	{
		for (double* p : {uA_, vA_, wA_, uB_, vB_, wB_, nut_, scratch_, planescr_, rhs_, pres_, divscr_}) cudaFree(p);
		cudaFree(solid_); cudaFree(nearsolid_);
		if (porous_k_) cudaFree(porous_k_);
	}

	void ChannelFluidCore::set_porous(const std::vector<double>& k_cell)
	{
		if ((int)k_cell.size() != g_.p_count()) return; // ignore a mismatched field
		if (!porous_k_) porous_k_ = dalloc0(g_.p_count());
		cudaMemcpy(porous_k_, k_cell.data(), sizeof(double) * g_.p_count(), cudaMemcpyHostToDevice);
	}

	void ChannelFluidCore::init_uniform(double u0, double v_blip)
	{
		std::vector<double> hu((size_t)g_.u_count(), u0), hv((size_t)g_.v_count(), 0.0), hw((size_t)g_.w_count(), 0.0);
		if (v_blip != 0.0 && g_bbox.any)
		{
			// One-sided transverse kick just downstream of the obstacle to break the
			// up/down symmetry and seed vortex shedding (deterministic).
			int i0 = g_bbox.i1, i1 = std::min(g_.nx - 1, g_bbox.i1 + 3 * (g_bbox.i1 - g_bbox.i0 + 1));
			int j0 = std::max(0, g_bbox.j0 - (g_bbox.j1 - g_bbox.j0)), j1 = std::min(g_.ny, g_bbox.j1 + (g_bbox.j1 - g_bbox.j0));
			for (int k = 0; k < g_.nz; ++k) for (int j = j0; j <= j1 && j <= g_.ny; ++j) for (int i = i0; i <= i1; ++i)
				hv[g_.vidx(i, j, k)] = v_blip;
		}
		cudaMemcpy(uA_, hu.data(), sizeof(double) * g_.u_count(), cudaMemcpyHostToDevice);
		cudaMemcpy(vA_, hv.data(), sizeof(double) * g_.v_count(), cudaMemcpyHostToDevice);
		cudaMemcpy(wA_, hw.data(), sizeof(double) * g_.w_count(), cudaMemcpyHostToDevice);
		ch_apply_inlet_gpu(uA_, g_, bc_);
		ch_apply_solid_bc_gpu(uA_, vA_, wA_, solid_, g_);
		if (bed_inlet_mask_) mask_inlet_solid();
		cudaMemset(pres_, 0, sizeof(double) * g_.p_count());
	}

	void ChannelFluidCore::init_inlet_profile()
	{
		std::vector<double> hu((size_t)g_.u_count(), 0.0), hv((size_t)g_.v_count(), 0.0), hw((size_t)g_.w_count(), 0.0);
		for (int k = 0; k < g_.nz; ++k) for (int j = 0; j < g_.ny; ++j) for (int i = 0; i <= g_.nx; ++i)
			hu[g_.uidx(i, j, k)] = channel_inlet_u(bc_, g_, k);
		cudaMemcpy(uA_, hu.data(), sizeof(double) * g_.u_count(), cudaMemcpyHostToDevice);
		cudaMemcpy(vA_, hv.data(), sizeof(double) * g_.v_count(), cudaMemcpyHostToDevice);
		cudaMemcpy(wA_, hw.data(), sizeof(double) * g_.w_count(), cudaMemcpyHostToDevice);
		ch_apply_inlet_gpu(uA_, g_, bc_);
		ch_apply_solid_bc_gpu(uA_, vA_, wA_, solid_, g_);
		if (bed_inlet_mask_) mask_inlet_solid();
		cudaMemset(pres_, 0, sizeof(double) * g_.p_count());
	}

	void ChannelFluidCore::update_solid(const std::vector<unsigned char>& solid)
	{
		std::vector<unsigned char> s = solid;
		if ((int)s.size() != g_.p_count()) return; // ignore a mismatched mask (keep the current one)
		// Recompute the near-solid dilation band (same rule as the constructor).
		std::vector<unsigned char> ns((size_t)g_.p_count(), 0);
		int band = pr_.advect_band;
		for (int k = 0; k < g_.nz; ++k) for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i)
		{
			bool near = false;
			for (int dk = -band; dk <= band && !near; ++dk) for (int dj = -band; dj <= band && !near; ++dj) for (int di = -band; di <= band && !near; ++di)
			{
				int ii = i + di, jj = j + dj, kk = k + dk;
				if (ii < 0 || ii >= g_.nx || jj < 0 || jj >= g_.ny || kk < 0 || kk >= g_.nz) continue;
				if (s[g_.pidx(ii, jj, kk)]) near = true;
			}
			ns[g_.pidx(i, j, k)] = near ? 1 : 0;
		}
		cudaMemcpy(solid_, s.data(), (size_t)g_.p_count(), cudaMemcpyHostToDevice);
		cudaMemcpy(nearsolid_, ns.data(), (size_t)g_.p_count(), cudaMemcpyHostToDevice);
		solver_->build_masks(solid_);
		g_bbox = solid_bbox(g_, s);
		ch_apply_solid_bc_gpu(uA_, vA_, wA_, solid_, g_); // zero velocities inside the new solid cells
		cudaDeviceSynchronize();
	}

	double ChannelFluidCore::adaptive_dt()
	{
		double maxu = std::max({reduce_max_abs_gpu(uA_, g_.u_count()), reduce_max_abs_gpu(vA_, g_.v_count()),
			reduce_max_abs_gpu(wA_, g_.w_count()), std::fabs(bc_.U_inlet)});
		double numax = pr_.nu;
		if (pr_.Cs > 0.0) { ch_smagorinsky_gpu(uA_, vA_, wA_, nut_, solid_, g_, bc_, pr_.Cs); numax = pr_.nu + reduce_max_abs_gpu(nut_, g_.p_count()); }
		double dt_adv = (maxu > 1e-12) ? pr_.cfl * g_.h / maxu : 1e30;
		double dt_diff = (numax > 0.0) ? g_.h * g_.h / (6.0 * numax) : 1e30;
		return pr_.safety * std::min(dt_adv, dt_diff);
	}

	double ChannelFluidCore::step()
	{
		double dt = (pr_.fixed_dt > 0.0) ? pr_.fixed_dt : adaptive_dt();
		last_dt_ = dt;

		const double* nutptr = nullptr;
		if (pr_.Cs > 0.0) { ch_smagorinsky_gpu(uA_, vA_, wA_, nut_, solid_, g_, bc_, pr_.Cs); nutptr = nut_; }

		ch_advect_gpu(uA_, vA_, wA_, uB_, vB_, wB_, scratch_, solid_, nearsolid_, g_, bc_, dt);
		ch_apply_solid_bc_gpu(uB_, vB_, wB_, solid_, g_);
		ch_diffuse_gpu(uB_, vB_, wB_, uA_, vA_, wA_, nutptr, solid_, g_, bc_, dt, pr_.nu);

		// M7 sub-grid porous/thin-screen momentum sink (perforated units). No-op if never set.
		if (porous_k_) ch_porous_drag_gpu(uA_, vA_, wA_, porous_k_, g_, dt);

		// M3 inlet turbulence: advance the generator once per step (null ⇒ no-op).
		if (fluct_) fluct_->advance(dt);

		// pre-projection velocity BCs
		ch_apply_inlet_gpu(uA_, g_, bc_);
		if (fluct_) { fluct_->add_u(uA_, g_); fluct_->set_vw(vA_, wA_, g_); } // SEM fluctuations
		ch_apply_solid_bc_gpu(uA_, vA_, wA_, solid_, g_);
		if (bed_inlet_mask_) mask_inlet_solid();
		ch_orlanski_gpu(uA_, g_, bc_, dt);
		if (bc_.flow_sign < 0)
		{
			// Reversed tide: inlet flux on xmax, outlet on xmin (both negative). Rescale the outlet plane
			// to match the inlet flux; the outflow guard is q_out < 0 (mirror of the forward q_out > 0).
			q_in_ = ch_uplane_flux_gpu(uA_, planescr_, g_, g_.nx);
			q_out_ = ch_uplane_flux_gpu(uA_, planescr_, g_, 0);
			if (q_out_ < -1e-9) { double s = q_in_ / q_out_; s = s < 0.5 ? 0.5 : (s > 2.0 ? 2.0 : s); ch_scale_uplane_gpu(uA_, g_, 0, s); }
		}
		else
		{
			q_in_ = ch_uplane_flux_gpu(uA_, planescr_, g_, 0);
			q_out_ = ch_uplane_flux_gpu(uA_, planescr_, g_, g_.nx);
			if (q_out_ > 1e-9) { double s = q_in_ / q_out_; s = s < 0.5 ? 0.5 : (s > 2.0 ? 2.0 : s); ch_scale_uplane_gpu(uA_, g_, g_.nx, s); } // CLAMP the outlet mass-rescale to [0.5,2]: a shrinking q_out (reversed inflow) would otherwise blow up q_in/q_out and amplify the reversal. No-op for M2 (q_out~q_in).
		}

		// projection (Dirichlet p=0 at the outlet — xmax by default, xmin under a reversed tide)
		ch_poisson_rhs_gpu(uA_, vA_, wA_, rhs_, solid_, g_, pr_.rho, dt);
		SolveResult r = solver_->solve(pres_, rhs_, pr_.proj_tol, pr_.proj_max_iter, /*warm_start=*/true);
		ch_subtract_gradient_gpu(uA_, vA_, wA_, pres_, solid_, g_, bc_, pr_.rho, dt, /*dir_xmax=*/bc_.flow_sign < 0 ? -1 : 1);

		// post-projection BCs
		ch_apply_inlet_gpu(uA_, g_, bc_);
		if (fluct_) fluct_->add_u(uA_, g_); // re-hold the turbulent Dirichlet inlet u
		ch_apply_solid_bc_gpu(uA_, vA_, wA_, solid_, g_);
		if (bed_inlet_mask_) mask_inlet_solid();

		last_iters_ = r.iters; last_relres_ = r.relres;
		q_out_ = ch_uplane_flux_gpu(uA_, planescr_, g_, bc_.flow_sign < 0 ? 0 : g_.nx);
		cudaDeviceSynchronize();
		return dt;
	}

	FluidSnapshot ChannelFluidCore::snapshot() const
	{
		FluidSnapshot s; s.grid = g_;
		s.u.resize(g_.u_count()); s.v.resize(g_.v_count()); s.w.resize(g_.w_count()); s.nut.resize(g_.p_count());
		cudaMemcpy(s.u.data(), uA_, sizeof(double) * g_.u_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(s.v.data(), vA_, sizeof(double) * g_.v_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(s.w.data(), wA_, sizeof(double) * g_.w_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(s.nut.data(), nut_, sizeof(double) * g_.p_count(), cudaMemcpyDeviceToHost);
		return s;
	}

	double ChannelFluidCore::get_u(int i, int j, int k) const
	{ double v = 0.0; cudaMemcpy(&v, uA_ + g_.uidx(i, j, k), sizeof(double), cudaMemcpyDeviceToHost); return v; }
	double ChannelFluidCore::get_v(int i, int j, int k) const
	{ double v = 0.0; cudaMemcpy(&v, vA_ + g_.vidx(i, j, k), sizeof(double), cudaMemcpyDeviceToHost); return v; }

	void ChannelFluidCore::copy_state_host(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p) const
	{
		u.resize((size_t)g_.u_count()); v.resize((size_t)g_.v_count()); w.resize((size_t)g_.w_count()); p.resize((size_t)g_.p_count());
		cudaMemcpy(u.data(), uA_, sizeof(double) * g_.u_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(v.data(), vA_, sizeof(double) * g_.v_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(w.data(), wA_, sizeof(double) * g_.w_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(p.data(), pres_, sizeof(double) * g_.p_count(), cudaMemcpyDeviceToHost);
	}

	void ChannelFluidCore::load_state_host(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& p)
	{
		if ((int)u.size() == g_.u_count()) cudaMemcpy(uA_, u.data(), sizeof(double) * g_.u_count(), cudaMemcpyHostToDevice);
		if ((int)v.size() == g_.v_count()) cudaMemcpy(vA_, v.data(), sizeof(double) * g_.v_count(), cudaMemcpyHostToDevice);
		if ((int)w.size() == g_.w_count()) cudaMemcpy(wA_, w.data(), sizeof(double) * g_.w_count(), cudaMemcpyHostToDevice);
		if ((int)p.size() == g_.p_count()) cudaMemcpy(pres_, p.data(), sizeof(double) * g_.p_count(), cudaMemcpyHostToDevice);
		// Re-hold the solid/inlet BCs so the restored faces are consistent with the constructed mask.
		ch_apply_inlet_gpu(uA_, g_, bc_);
		ch_apply_solid_bc_gpu(uA_, vA_, wA_, solid_, g_);
		if (bed_inlet_mask_) mask_inlet_solid();
		cudaDeviceSynchronize();
	}
}
