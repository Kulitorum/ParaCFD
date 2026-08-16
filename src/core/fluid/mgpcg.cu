// mgpcg.cu — MGPCG solver implementation. See mgpcg.h. RESEARCH §3.3.
#include "core/fluid/mgpcg.h"
#include "core/fluid/mac_ops.h"

#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>

namespace paracfd::core
{
	namespace
	{
		double* dalloc(int n)
		{
			double* p = nullptr;
			cudaMalloc(&p, sizeof(double) * (size_t)n);
			cudaMemset(p, 0, sizeof(double) * (size_t)n);
			return p;
		}
	}

	MgpcgSolver::MgpcgSolver(MacGrid finest)
	{
		// Build the factor-2 hierarchy: coarsen while every dim is even and stays >= 2.
		MacGrid g = finest;
		grids_.push_back(g);
		while (g.nx % 2 == 0 && g.ny % 2 == 0 && g.nz % 2 == 0 && g.nx / 2 >= 2 && g.ny / 2 >= 2 && g.nz / 2 >= 2)
		{
			MacGrid c;
			c.nx = g.nx / 2; c.ny = g.ny / 2; c.nz = g.nz / 2; c.h = g.h * 2.0;
			grids_.push_back(c);
			g = c;
		}
		for (auto& lg : grids_)
		{
			int n = lg.p_count();
			Lp_.push_back(dalloc(n)); Lrhs_.push_back(dalloc(n));
			Lres_.push_back(dalloc(n)); Ltmp_.push_back(dalloc(n));
		}
		n0_ = grids_[0].p_count();
		b0_ = dalloc(n0_); r_ = dalloc(n0_); z_ = dalloc(n0_); s_ = dalloc(n0_); As_ = dalloc(n0_);
	}

	MgpcgSolver::~MgpcgSolver()
	{
		for (auto p : Lp_) cudaFree(p);
		for (auto p : Lrhs_) cudaFree(p);
		for (auto p : Lres_) cudaFree(p);
		for (auto p : Ltmp_) cudaFree(p);
		cudaFree(b0_); cudaFree(r_); cudaFree(z_); cudaFree(s_); cudaFree(As_);
	}

	void MgpcgSolver::apply_finest(const double* p, double* Ap) { poisson_apply_gpu(p, Ap, grids_[0]); }

	void MgpcgSolver::vcycle(int level)
	{
		MacGrid g = grids_[level];
		int n = g.p_count();
		double* p = Lp_[level];
		double* rhs = Lrhs_[level];
		double* res = Lres_[level];
		double* tmp = Ltmp_[level];

		if (level == (int)grids_.size() - 1)
		{
			// Coarsest: heavy damped-Jacobi smoothing; keep the gauge mean-zero.
			jacobi_smooth_gpu(p, rhs, tmp, g, omega, coarse_sweeps);
			double mean = reduce_sum_gpu(p, n) / n;
			add_const_gpu(p, -mean, n);
			return;
		}

		int gsweeps = 1 << (level + 1); // 2 at finest, doubling per level (RESEARCH §3.3)
		// pre-smooth (jacobi then forward GS band)
		jacobi_smooth_gpu(p, rhs, tmp, g, omega, pre_post_jacobi);
		gs_band_gpu(p, rhs, g, gs_band, gsweeps, true);

		// residual -> restrict to coarser rhs, zero coarse p
		poisson_residual_gpu(p, rhs, res, g);
		MacGrid gc = grids_[level + 1];
		restrict_gpu(res, Lrhs_[level + 1], g, gc);
		cudaMemset(Lp_[level + 1], 0, sizeof(double) * gc.p_count());

		vcycle(level + 1);

		prolong_add_gpu(Lp_[level + 1], p, gc, g);

		// post-smooth (backward GS band then jacobi) — symmetric to the pre-smooth
		gs_band_gpu(p, rhs, g, gs_band, gsweeps, false);
		jacobi_smooth_gpu(p, rhs, tmp, g, omega, pre_post_jacobi);
	}

	void MgpcgSolver::precondition(const double* r, double* z)
	{
		cudaMemcpy(Lrhs_[0], r, sizeof(double) * n0_, cudaMemcpyDeviceToDevice);
		cudaMemset(Lp_[0], 0, sizeof(double) * n0_);
		vcycle(0);
		double mean = reduce_sum_gpu(Lp_[0], n0_) / n0_;
		add_const_gpu(Lp_[0], -mean, n0_);
		cudaMemcpy(z, Lp_[0], sizeof(double) * n0_, cudaMemcpyDeviceToDevice);
	}

	SolveResult MgpcgSolver::solve(double* x, const double* b, double tol, int max_iter, bool warm_start, bool use_mg)
	{
		MacGrid g = grids_[0];
		int n = n0_;
		SolveResult out;

		// b0 = b - mean(b)  (all-Neumann compatibility)
		cudaMemcpy(b0_, b, sizeof(double) * n, cudaMemcpyDeviceToDevice);
		double bmean = reduce_sum_gpu(b0_, n) / n;
		add_const_gpu(b0_, -bmean, n);
		double bnorm = std::sqrt(dot_gpu(b0_, b0_, n));
		if (!(bnorm > 0.0))
		{
			cudaMemset(x, 0, sizeof(double) * n);
			out.converged = true; out.relres = 0.0; return out;
		}

		if (!warm_start) cudaMemset(x, 0, sizeof(double) * n);

		if (!use_mg)
		{
			// Plain damped-Jacobi fallback (bring-up/debug). Iterate to tol; each sweep
			// counts as an iteration. RESEARCH §3.3 (Jacobi = debug-only criterion).
			double relres = 1.0; int it = 0;
			for (; it < max_iter; ++it)
			{
				jacobi_smooth_gpu(x, b0_, z_, g, omega, 1);
				double mean = reduce_sum_gpu(x, n) / n; add_const_gpu(x, -mean, n);
				poisson_residual_gpu(x, b0_, r_, g);
				relres = std::sqrt(dot_gpu(r_, r_, n)) / bnorm;
				if (relres <= tol) { ++it; break; }
			}
			out.iters = it; out.relres = relres; out.converged = relres <= tol; return out;
		}

		// r = b0 - A x
		if (warm_start) { apply_finest(x, As_); cudaMemcpy(r_, b0_, sizeof(double) * n, cudaMemcpyDeviceToDevice); axpy_gpu(-1.0, As_, r_, n); }
		else cudaMemcpy(r_, b0_, sizeof(double) * n, cudaMemcpyDeviceToDevice);

		double relres = std::sqrt(dot_gpu(r_, r_, n)) / bnorm;
		if (relres <= tol) { out.iters = 0; out.relres = relres; out.converged = true; double m = reduce_sum_gpu(x, n) / n; add_const_gpu(x, -m, n); return out; }

		precondition(r_, z_);
		cudaMemcpy(s_, z_, sizeof(double) * n, cudaMemcpyDeviceToDevice);
		double rz = dot_gpu(r_, z_, n);

		int it = 0;
		for (; it < max_iter; )
		{
			++it;
			apply_finest(s_, As_);
			double sAs = dot_gpu(s_, As_, n);
			if (!(sAs > 0.0)) break; // safeguard (should not trigger for SPD Poisson)
			double alpha = rz / sAs;
			axpy_gpu(alpha, s_, x, n);      // x += alpha s
			axpy_gpu(-alpha, As_, r_, n);   // r -= alpha As
			relres = std::sqrt(dot_gpu(r_, r_, n)) / bnorm;
			if (relres <= tol) break;
			precondition(r_, z_);
			double rznew = dot_gpu(r_, z_, n);
			double beta = rznew / rz;
			scale_add_gpu(s_, 1.0, z_, beta, n); // s = z + beta s
			rz = rznew;
		}
		double m = reduce_sum_gpu(x, n) / n; add_const_gpu(x, -m, n); // pin gauge
		out.iters = it; out.relres = relres; out.converged = relres <= tol;
		return out;
	}
}
