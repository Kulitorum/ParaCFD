// stam_fluid_core.cu — see stam_fluid_core.h. RESEARCH §3.
#include "core/fluid/stam_fluid_core.h"
#include "core/fluid/mac_ops.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace windcfd::core
{
	namespace
	{
		double* dalloc0(int n)
		{
			double* p = nullptr;
			cudaMalloc(&p, sizeof(double) * (size_t)n);
			cudaMemset(p, 0, sizeof(double) * (size_t)n);
			return p;
		}
		int maxface(MacGrid g) { return std::max(g.u_count(), std::max(g.v_count(), g.w_count())); }
	}

	StamFluidCore::StamFluidCore(MacGrid grid, BC bc, StamParams params)
		: g_(grid), bc_(bc), pr_(params)
	{
		solver_ = std::make_unique<MgpcgSolver>(g_);
		uA_ = dalloc0(g_.u_count()); vA_ = dalloc0(g_.v_count()); wA_ = dalloc0(g_.w_count());
		uB_ = dalloc0(g_.u_count()); vB_ = dalloc0(g_.v_count()); wB_ = dalloc0(g_.w_count());
		pu_ = dalloc0(g_.u_count()); pv_ = dalloc0(g_.v_count()); pw_ = dalloc0(g_.w_count());
		um_ = dalloc0(g_.u_count()); vm_ = dalloc0(g_.v_count()); wm_ = dalloc0(g_.w_count());
		nut_ = dalloc0(g_.p_count());
		scratch_ = dalloc0(maxface(g_));
		rhs_ = dalloc0(g_.p_count()); pres_ = dalloc0(g_.p_count()); divscr_ = dalloc0(g_.p_count());
	}

	StamFluidCore::~StamFluidCore()
	{
		for (double* p : {uA_, vA_, wA_, uB_, vB_, wB_, pu_, pv_, pw_, um_, vm_, wm_, nut_, scratch_, rhs_, pres_, divscr_})
			cudaFree(p);
	}

	double StamFluidCore::step()
	{
		// --- adaptive dt (advective CFL + explicit-diffusion stability) ----------
		double maxu = std::max({reduce_max_abs_gpu(uA_, g_.u_count()),
			reduce_max_abs_gpu(vA_, g_.v_count()), reduce_max_abs_gpu(wA_, g_.w_count())});
		double u_eff = std::max(maxu, bc_.lid_u);
		double numax = pr_.nu;
		if (pr_.Cs > 0.0)
		{
			smagorinsky_nut_gpu(uA_, vA_, wA_, nut_, g_, bc_, pr_.Cs);
			numax = pr_.nu + reduce_max_abs_gpu(nut_, g_.p_count());
		}
		double dt_adv = (u_eff > 1e-12) ? pr_.cfl * g_.h / u_eff : 1e30;
		double dt_diff = (numax > 0.0) ? g_.h * g_.h / (6.0 * numax) : 1e30;
		double dt = pr_.safety * std::min(dt_adv, dt_diff);
		last_dt_ = dt;

		// --- advect A -> B (MacCormack, 1st-order reversion within 1 cell of walls) -
		mac_advect_gpu(uA_, vA_, wA_, uB_, vB_, wB_, scratch_, g_, bc_, dt, /*band=*/1);

		// --- explicit diffusion B -> A ------------------------------------------
		mac_diffuse_gpu(uB_, vB_, wB_, uA_, vA_, wA_, (pr_.Cs > 0.0 ? nut_ : nullptr), g_, bc_, dt, pr_.nu);

		// --- projection ---------------------------------------------------------
		poisson_rhs_gpu(uA_, vA_, wA_, rhs_, g_, pr_.rho, dt);
		SolveResult r = solver_->solve(pres_, rhs_, pr_.proj_tol, pr_.proj_max_iter, /*warm_start=*/true, pr_.use_mg);
		subtract_gradient_gpu(uA_, vA_, wA_, pres_, g_, pr_.rho, dt);
		last_iters_ = r.iters; last_relres_ = r.relres;
		cudaDeviceSynchronize();
		return dt;
	}

	double StamFluidCore::max_velocity_change_and_remark()
	{
		double result = 1e30;
		if (marked_)
		{
			// scratch_ big enough for the largest component; reuse per component.
			cudaMemcpy(scratch_, uA_, sizeof(double) * g_.u_count(), cudaMemcpyDeviceToDevice);
			axpy_gpu(-1.0, um_, scratch_, g_.u_count());
			double du = reduce_max_abs_gpu(scratch_, g_.u_count());
			cudaMemcpy(scratch_, vA_, sizeof(double) * g_.v_count(), cudaMemcpyDeviceToDevice);
			axpy_gpu(-1.0, vm_, scratch_, g_.v_count());
			double dv = reduce_max_abs_gpu(scratch_, g_.v_count());
			cudaMemcpy(scratch_, wA_, sizeof(double) * g_.w_count(), cudaMemcpyDeviceToDevice);
			axpy_gpu(-1.0, wm_, scratch_, g_.w_count());
			double dw = reduce_max_abs_gpu(scratch_, g_.w_count());
			result = std::max({du, dv, dw});
		}
		cudaMemcpy(um_, uA_, sizeof(double) * g_.u_count(), cudaMemcpyDeviceToDevice);
		cudaMemcpy(vm_, vA_, sizeof(double) * g_.v_count(), cudaMemcpyDeviceToDevice);
		cudaMemcpy(wm_, wA_, sizeof(double) * g_.w_count(), cudaMemcpyDeviceToDevice);
		marked_ = true;
		return result;
	}

	ProjectionProbe StamFluidCore::projection_probe(double tol, int max_iter, bool use_mg)
	{
		// Build a realistic divergent field u* from the current state (advect+diffuse)
		// in scratch buffers, then project it cold-start. Live state (A) is untouched.
		double maxu = std::max({reduce_max_abs_gpu(uA_, g_.u_count()),
			reduce_max_abs_gpu(vA_, g_.v_count()), reduce_max_abs_gpu(wA_, g_.w_count())});
		double u_eff = std::max(maxu, bc_.lid_u);
		double numax = pr_.nu;
		if (pr_.Cs > 0.0) { smagorinsky_nut_gpu(uA_, vA_, wA_, nut_, g_, bc_, pr_.Cs); numax = pr_.nu + reduce_max_abs_gpu(nut_, g_.p_count()); }
		double dt_adv = (u_eff > 1e-12) ? pr_.cfl * g_.h / u_eff : 1e30;
		double dt_diff = (numax > 0.0) ? g_.h * g_.h / (6.0 * numax) : 1e30;
		double dt = pr_.safety * std::min(dt_adv, dt_diff);

		mac_advect_gpu(uA_, vA_, wA_, uB_, vB_, wB_, scratch_, g_, bc_, dt, 1);
		mac_diffuse_gpu(uB_, vB_, wB_, pu_, pv_, pw_, (pr_.Cs > 0.0 ? nut_ : nullptr), g_, bc_, dt, pr_.nu);

		ProjectionProbe out;
		out.maxdiv_before = max_abs_divergence_gpu(pu_, pv_, pw_, divscr_, g_);
		poisson_rhs_gpu(pu_, pv_, pw_, rhs_, g_, pr_.rho, dt);
		out.solve = solver_->solve(pres_, rhs_, tol, max_iter, /*warm_start=*/false, use_mg);
		subtract_gradient_gpu(pu_, pv_, pw_, pres_, g_, pr_.rho, dt);
		out.maxdiv_after = max_abs_divergence_gpu(pu_, pv_, pw_, divscr_, g_);
		cudaDeviceSynchronize();
		return out;
	}

	FluidSnapshot StamFluidCore::snapshot() const
	{
		FluidSnapshot s;
		s.grid = g_;
		s.u.resize(g_.u_count()); s.v.resize(g_.v_count()); s.w.resize(g_.w_count()); s.nut.resize(g_.p_count());
		cudaMemcpy(s.u.data(), uA_, sizeof(double) * g_.u_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(s.v.data(), vA_, sizeof(double) * g_.v_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(s.w.data(), wA_, sizeof(double) * g_.w_count(), cudaMemcpyDeviceToHost);
		cudaMemcpy(s.nut.data(), nut_, sizeof(double) * g_.p_count(), cudaMemcpyDeviceToHost);
		return s;
	}
}
