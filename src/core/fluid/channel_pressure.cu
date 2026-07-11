// channel_pressure.cu — see channel_pressure.h. RESEARCH §3.3.
#include "core/fluid/channel_pressure.h"
#include "core/fluid/channel_ops.h"
#include "core/fluid/mac_ops.h" // restrict/prolong transfers + reductions

#include <cuda_runtime.h>

#include <cmath>

namespace windcfd::core
{
	namespace
	{
		double* dalloc(int n) { double* p = nullptr; cudaMalloc(&p, sizeof(double) * (size_t)n); cudaMemset(p, 0, sizeof(double) * (size_t)n); return p; }
		unsigned char* ualloc(int n) { unsigned char* p = nullptr; cudaMalloc(&p, (size_t)n); cudaMemset(p, 0, (size_t)n); return p; }
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }

		__global__ void k_coarsen_mask(const unsigned char* fine, unsigned char* coarse, MacGrid gf, MacGrid gc, int nc)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= nc) return;
			int ci = t % gc.nx, cj = (t / gc.nx) % gc.ny, ck = t / (gc.nx * gc.ny);
			int all = 1;
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c)
			{
				int fi = 2 * ci + a, fj = 2 * cj + b, fk = 2 * ck + c;
				unsigned char s = (fi < gf.nx && fj < gf.ny && fk < gf.nz) ? fine[gf.pidx(fi, fj, fk)] : 0;
				if (!s) all = 0;
			}
			coarse[t] = (unsigned char)all;
		}
		__global__ void k_zero_solids(double* p, const unsigned char* solid, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			if (solid[t]) p[t] = 0.0;
		}
	}

	ChannelMgpcg::ChannelMgpcg(MacGrid finest)
	{
		// Isotropic factor-2 coarsening. A thin spanwise axis (e.g. nz=4) is allowed to
		// coarsen to nz=1 (a valid 2-D coarse Poisson operator with no z-coupling), which
		// deepens the hierarchy and keeps the preconditioner strong for quasi-2-D runs.
		MacGrid g = finest; grids_.push_back(g);
		while (g.nx % 2 == 0 && g.ny % 2 == 0 && g.nz % 2 == 0 && g.nx / 2 >= 1 && g.ny / 2 >= 1 && g.nz / 2 >= 1)
		{
			MacGrid c; c.nx = g.nx / 2; c.ny = g.ny / 2; c.nz = g.nz / 2; c.h = g.h * 2.0;
			grids_.push_back(c); g = c;
		}
		for (auto& lg : grids_)
		{
			int n = lg.p_count();
			Lp_.push_back(dalloc(n)); Lrhs_.push_back(dalloc(n)); Ltmp_.push_back(dalloc(n));
			Lsolid_.push_back(ualloc(n));
		}
		n0_ = grids_[0].p_count();
		r_ = dalloc(n0_); z_ = dalloc(n0_); s_ = dalloc(n0_); As_ = dalloc(n0_);
	}
	ChannelMgpcg::~ChannelMgpcg()
	{
		for (auto p : Lp_) cudaFree(p);
		for (auto p : Lrhs_) cudaFree(p);
		for (auto p : Ltmp_) cudaFree(p);
		for (auto p : Lsolid_) cudaFree(p);
		cudaFree(r_); cudaFree(z_); cudaFree(s_); cudaFree(As_);
	}

	void ChannelMgpcg::build_masks(const unsigned char* finest_solid)
	{
		cudaMemcpy(Lsolid_[0], finest_solid, (size_t)n0_, cudaMemcpyDeviceToDevice);
		for (size_t l = 1; l < grids_.size(); ++l)
		{
			int nc = grids_[l].p_count();
			k_coarsen_mask<<<gsz(nc), 256>>>(Lsolid_[l - 1], Lsolid_[l], grids_[l - 1], grids_[l], nc);
		}
	}

	void ChannelMgpcg::apply_finest(const double* p, double* Ap) { ch_poisson_apply_gpu(p, Ap, Lsolid_[0], grids_[0], dir_xmax_); }

	void ChannelMgpcg::vcycle(int level)
	{
		MacGrid g = grids_[level]; int n = g.p_count();
		double* p = Lp_[level]; double* rhs = Lrhs_[level]; double* tmp = Ltmp_[level];
		const unsigned char* solid = Lsolid_[level];
		if (level == (int)grids_.size() - 1)
		{
			ch_jacobi_gpu(p, rhs, tmp, solid, g, dir_xmax_, omega, coarse_sweeps);
			return;
		}
		int gsweeps = 1 << (level + 1); if (gsweeps > 8) gsweeps = 8; // cap doubling on deep (thin-domain) hierarchies
		ch_jacobi_gpu(p, rhs, tmp, solid, g, dir_xmax_, omega, pre_post_jacobi);
		ch_gs_band_gpu(p, rhs, solid, g, dir_xmax_, gs_band, gsweeps, true);
		ch_poisson_residual_gpu(p, rhs, tmp, solid, g, dir_xmax_); // tmp = residual
		MacGrid gc = grids_[level + 1];
		restrict_gpu(tmp, Lrhs_[level + 1], g, gc);
		cudaMemset(Lp_[level + 1], 0, sizeof(double) * gc.p_count());
		vcycle(level + 1);
		prolong_add_gpu(Lp_[level + 1], p, gc, g);
		ch_gs_band_gpu(p, rhs, solid, g, dir_xmax_, gs_band, gsweeps, false);
		ch_jacobi_gpu(p, rhs, tmp, solid, g, dir_xmax_, omega, pre_post_jacobi);
	}
	void ChannelMgpcg::precondition(const double* r, double* z)
	{
		cudaMemcpy(Lrhs_[0], r, sizeof(double) * n0_, cudaMemcpyDeviceToDevice);
		cudaMemset(Lp_[0], 0, sizeof(double) * n0_);
		vcycle(0);
		k_zero_solids<<<gsz(n0_), 256>>>(Lp_[0], Lsolid_[0], grids_[0], n0_);
		cudaMemcpy(z, Lp_[0], sizeof(double) * n0_, cudaMemcpyDeviceToDevice);
	}

	SolveResult ChannelMgpcg::solve(double* x, const double* b, double tol, int max_iter, bool warm_start)
	{
		int n = n0_; SolveResult out;
		double bnorm = std::sqrt(dot_gpu(b, b, n));
		if (!(bnorm > 0.0)) { cudaMemset(x, 0, sizeof(double) * n); out.converged = true; return out; }
		if (!warm_start) cudaMemset(x, 0, sizeof(double) * n);

		// r = b - A x
		if (warm_start) { apply_finest(x, As_); cudaMemcpy(r_, b, sizeof(double) * n, cudaMemcpyDeviceToDevice); axpy_gpu(-1.0, As_, r_, n); }
		else cudaMemcpy(r_, b, sizeof(double) * n, cudaMemcpyDeviceToDevice);

		double relres = std::sqrt(dot_gpu(r_, r_, n)) / bnorm;
		if (relres <= tol) { out.iters = 0; out.relres = relres; out.converged = true; return out; }

		precondition(r_, z_);
		cudaMemcpy(s_, z_, sizeof(double) * n, cudaMemcpyDeviceToDevice);
		double rz = dot_gpu(r_, z_, n);
		int it = 0;
		for (; it < max_iter; )
		{
			++it;
			apply_finest(s_, As_);
			double sAs = dot_gpu(s_, As_, n);
			if (!(sAs > 0.0)) break;
			double alpha = rz / sAs;
			axpy_gpu(alpha, s_, x, n);
			axpy_gpu(-alpha, As_, r_, n);
			relres = std::sqrt(dot_gpu(r_, r_, n)) / bnorm;
			if (relres <= tol) break;
			precondition(r_, z_);
			double rznew = dot_gpu(r_, z_, n);
			double beta = rznew / rz;
			scale_add_gpu(s_, 1.0, z_, beta, n);
			rz = rznew;
		}
		out.iters = it; out.relres = relres; out.converged = relres <= tol;
		return out;
	}
}
