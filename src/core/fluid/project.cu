// project.cu — MAC pressure projection. Divergence, the all-Neumann 7-point Poisson
// operator (solid/box neighbours dropped => dp/dn=0), damped-Jacobi + boundary-band
// Gauss-Seidel smoothers, factor-2 restriction/prolongation, and an MGPCG solver
// (one geometric V-cycle preconditioner inside CG) with a plain-Jacobi fallback.
// RESEARCH §3.3 (MGPCG spec: V-cycle preconditioner, damped Jacobi omega=2/3, factor-2
// coarsening, boundary-band GS sweeps), verified research/08 item 5. Closed box is
// singular: rhs made zero-mean, pressure gauge pinned by mean-removal.
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h"

#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/inner_product.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>

#include <cmath>
#include <cstdio>
#include <vector>

namespace scour::core
{
	namespace
	{
		SCOUR_HD inline int nnb_count(MacGrid g, int i, int j, int k)
		{
			return (i > 0) + (i < g.nx - 1) + (j > 0) + (j < g.ny - 1) + (k > 0) + (k < g.nz - 1);
		}
		SCOUR_HD inline double nb_sum(const double* p, MacGrid g, int i, int j, int k)
		{
			double s = 0.0;
			if (i > 0) s += p[g.pidx(i - 1, j, k)];
			if (i < g.nx - 1) s += p[g.pidx(i + 1, j, k)];
			if (j > 0) s += p[g.pidx(i, j - 1, k)];
			if (j < g.ny - 1) s += p[g.pidx(i, j + 1, k)];
			if (k > 0) s += p[g.pidx(i, j, k - 1)];
			if (k < g.nz - 1) s += p[g.pidx(i, j, k + 1)];
			return s;
		}

		SCOUR_HD inline double div_cell(const double* u, const double* v, const double* w, MacGrid g, int i, int j, int k)
		{
			return (u[g.uidx(i + 1, j, k)] - u[g.uidx(i, j, k)]
				+ v[g.vidx(i, j + 1, k)] - v[g.vidx(i, j, k)]
				+ w[g.widx(i, j, k + 1)] - w[g.widx(i, j, k)]) / g.h;
		}

		__global__ void k_rhs(const double* u, const double* v, const double* w, double* rhs, MacGrid g, double rho, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			rhs[t] = -(rho / dt) * div_cell(u, v, w, g, i, j, k);
		}
		__global__ void k_apply(const double* p, double* Ap, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			Ap[t] = (nnb_count(g, i, j, k) * p[t] - nb_sum(p, g, i, j, k)) / (g.h * g.h);
		}
		__global__ void k_residual(const double* p, const double* rhs, double* r, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			r[t] = rhs[t] - (nnb_count(g, i, j, k) * p[t] - nb_sum(p, g, i, j, k)) / (g.h * g.h);
		}
		__global__ void k_jacobi(const double* p, const double* rhs, double* pout, MacGrid g, double omega, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			int cnt = nnb_count(g, i, j, k);
			if (cnt == 0) { pout[t] = p[t]; return; }
			double diag = cnt / (g.h * g.h);
			double Ap = (cnt * p[t] - nb_sum(p, g, i, j, k)) / (g.h * g.h);
			pout[t] = p[t] + omega * (rhs[t] - Ap) / diag;
		}
		__global__ void k_gs_color(double* p, const double* rhs, MacGrid g, int band, int color, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			if (((i + j + k) & 1) != color) return;
			bool inband = i < band || i >= g.nx - band || j < band || j >= g.ny - band || k < band || k >= g.nz - band;
			if (!inband) return;
			int cnt = nnb_count(g, i, j, k); if (cnt == 0) return;
			p[t] = (rhs[t] * g.h * g.h + nb_sum(p, g, i, j, k)) / cnt;
		}

		// factor-2 transfer: fine cell fi maps to coarse index space cpos = fi/2 - 0.25.
		SCOUR_HD inline void axis_stencil(int fi, int nc, int& C0, int& C1, double& w0, double& w1)
		{
			double cpos = (fi + 0.5) * 0.5 - 0.5;
			int Cb = (int)floor(cpos); double t = cpos - Cb;
			C0 = Cb; C1 = Cb + 1; w0 = 1.0 - t; w1 = t;
			if (C0 < 0) C0 = 0; else if (C0 > nc - 1) C0 = nc - 1;
			if (C1 < 0) C1 = 0; else if (C1 > nc - 1) C1 = nc - 1;
		}
		__global__ void k_restrict(const double* fine, double* coarse, MacGrid gf, MacGrid gc, int nf)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= nf) return;
			int fi = t % gf.nx, fj = (t / gf.nx) % gf.ny, fk = t / (gf.nx * gf.ny);
			int Cx0, Cx1, Cy0, Cy1, Cz0, Cz1; double wx0, wx1, wy0, wy1, wz0, wz1;
			axis_stencil(fi, gc.nx, Cx0, Cx1, wx0, wx1);
			axis_stencil(fj, gc.ny, Cy0, Cy1, wy0, wy1);
			axis_stencil(fk, gc.nz, Cz0, Cz1, wz0, wz1);
			double f = fine[t] * 0.125; // R = (1/8) P^T
			int cx[2] = {Cx0, Cx1}, cy[2] = {Cy0, Cy1}, cz[2] = {Cz0, Cz1};
			double wx[2] = {wx0, wx1}, wy[2] = {wy0, wy1}, wz[2] = {wz0, wz1};
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c)
				atomicAdd(&coarse[gc.pidx(cx[a], cy[b], cz[c])], f * wx[a] * wy[b] * wz[c]);
		}
		__global__ void k_prolong_add(const double* coarse, double* fine, MacGrid gc, MacGrid gf, int nf)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= nf) return;
			int fi = t % gf.nx, fj = (t / gf.nx) % gf.ny, fk = t / (gf.nx * gf.ny);
			int Cx0, Cx1, Cy0, Cy1, Cz0, Cz1; double wx0, wx1, wy0, wy1, wz0, wz1;
			axis_stencil(fi, gc.nx, Cx0, Cx1, wx0, wx1);
			axis_stencil(fj, gc.ny, Cy0, Cy1, wy0, wy1);
			axis_stencil(fk, gc.nz, Cz0, Cz1, wz0, wz1);
			int cx[2] = {Cx0, Cx1}, cy[2] = {Cy0, Cy1}, cz[2] = {Cz0, Cz1};
			double wx[2] = {wx0, wx1}, wy[2] = {wy0, wy1}, wz[2] = {wz0, wz1};
			double v = 0.0;
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c)
				v += wx[a] * wy[b] * wz[c] * coarse[gc.pidx(cx[a], cy[b], cz[c])];
			fine[t] += v;
		}
		__global__ void k_subgrad(double* u, double* v, double* w, const double* p, MacGrid g, double coef)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x;
			int nu = (g.nx + 1) * g.ny * g.nz;
			int nv = g.nx * (g.ny + 1) * g.nz;
			int nw = g.nx * g.ny * (g.nz + 1);
			if (t < nu)
			{
				int i = t % (g.nx + 1), j = (t / (g.nx + 1)) % g.ny, k = t / ((g.nx + 1) * g.ny);
				if (i >= 1 && i <= g.nx - 1) u[t] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i - 1, j, k)]) / g.h;
			}
			else if (t < nu + nv)
			{
				int s = t - nu; int i = s % g.nx, j = (s / g.nx) % (g.ny + 1), k = s / (g.nx * (g.ny + 1));
				if (j >= 1 && j <= g.ny - 1) v[s] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i, j - 1, k)]) / g.h;
			}
			else if (t < nu + nv + nw)
			{
				int s = t - nu - nv; int i = s % g.nx, j = (s / g.nx) % g.ny, k = s / (g.nx * g.ny);
				if (k >= 1 && k <= g.nz - 1) w[s] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i, j, k - 1)]) / g.h;
			}
		}
		__global__ void k_divabs(const double* u, const double* v, const double* w, double* out, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			out[t] = fabs(div_cell(u, v, w, g, i, j, k));
		}
		__global__ void k_axpy(double alpha, const double* x, double* y, int n)
		{ int t = blockIdx.x * blockDim.x + threadIdx.x; if (t < n) y[t] += alpha * x[t]; }
		__global__ void k_scale_add(double* y, double alpha, const double* x, double beta, int n)
		{ int t = blockIdx.x * blockDim.x + threadIdx.x; if (t < n) y[t] = alpha * x[t] + beta * y[t]; }
		__global__ void k_add_const(double* x, double c, int n)
		{ int t = blockIdx.x * blockDim.x + threadIdx.x; if (t < n) x[t] += c; }

		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }

		struct AbsOp { __host__ __device__ double operator()(double v) const { return fabs(v); } };
	}

	// ---- public launchers -------------------------------------------------------
	void poisson_rhs_gpu(const double* u, const double* v, const double* w, double* rhs, MacGrid g, double rho, double dt)
	{ int n = g.p_count(); k_rhs<<<gsz(n), 256>>>(u, v, w, rhs, g, rho, dt, n); }
	void poisson_apply_gpu(const double* p, double* Ap, MacGrid g)
	{ int n = g.p_count(); k_apply<<<gsz(n), 256>>>(p, Ap, g, n); }
	void poisson_residual_gpu(const double* p, const double* rhs, double* r, MacGrid g)
	{ int n = g.p_count(); k_residual<<<gsz(n), 256>>>(p, rhs, r, g, n); }
	void jacobi_smooth_gpu(double* p, const double* rhs, double* scratch, MacGrid g, double omega, int sweeps)
	{
		int n = g.p_count(); double* src = p; double* dst = scratch;
		for (int s = 0; s < sweeps; ++s) { k_jacobi<<<gsz(n), 256>>>(src, rhs, dst, g, omega, n); double* tmp = src; src = dst; dst = tmp; }
		if (src != p) cudaMemcpy(p, src, sizeof(double) * n, cudaMemcpyDeviceToDevice);
	}
	void gs_band_gpu(double* p, const double* rhs, MacGrid g, int band, int sweeps, bool forward)
	{
		int n = g.p_count();
		for (int s = 0; s < sweeps; ++s)
		{
			int c0 = forward ? 0 : 1, c1 = forward ? 1 : 0;
			k_gs_color<<<gsz(n), 256>>>(p, rhs, g, band, c0, n);
			k_gs_color<<<gsz(n), 256>>>(p, rhs, g, band, c1, n);
		}
	}
	void restrict_gpu(const double* fine, double* coarse, MacGrid gf, MacGrid gc)
	{ int nf = gf.p_count(); cudaMemset(coarse, 0, sizeof(double) * gc.p_count()); k_restrict<<<gsz(nf), 256>>>(fine, coarse, gf, gc, nf); }
	void prolong_add_gpu(const double* coarse, double* fine, MacGrid gc, MacGrid gf)
	{ int nf = gf.p_count(); k_prolong_add<<<gsz(nf), 256>>>(coarse, fine, gc, gf, nf); }
	void subtract_gradient_gpu(double* u, double* v, double* w, const double* p, MacGrid g, double rho, double dt)
	{
		int nu = g.u_count(), nv = g.v_count(), nw = g.w_count();
		k_subgrad<<<gsz(nu + nv + nw), 256>>>(u, v, w, p, g, dt / rho);
	}
	double max_abs_divergence_gpu(const double* u, const double* v, const double* w, double* divscratch, MacGrid g)
	{
		int n = g.p_count(); k_divabs<<<gsz(n), 256>>>(u, v, w, divscratch, g, n);
		return reduce_max_abs_gpu(divscratch, n);
	}

	double reduce_sum_gpu(const double* x, int n)
	{ return thrust::reduce(thrust::device, thrust::device_pointer_cast(x), thrust::device_pointer_cast(x) + n, 0.0, thrust::plus<double>()); }
	double reduce_max_abs_gpu(const double* x, int n)
	{
		auto p = thrust::device_pointer_cast(x);
		return thrust::transform_reduce(thrust::device, p, p + n, AbsOp(), 0.0, thrust::maximum<double>());
	}
	double dot_gpu(const double* a, const double* b, int n)
	{ return thrust::inner_product(thrust::device, thrust::device_pointer_cast(a), thrust::device_pointer_cast(a) + n, thrust::device_pointer_cast(b), 0.0); }
	void axpy_gpu(double alpha, const double* x, double* y, int n) { k_axpy<<<gsz(n), 256>>>(alpha, x, y, n); }
	void scale_add_gpu(double* y, double alpha, const double* x, double beta, int n) { k_scale_add<<<gsz(n), 256>>>(y, alpha, x, beta, n); }
	void add_const_gpu(double* x, double c, int n) { k_add_const<<<gsz(n), 256>>>(x, c, n); }

	// ---- CPU references ---------------------------------------------------------
	void poisson_rhs_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, std::vector<double>& rhs, MacGrid g, double rho, double dt)
	{
		rhs.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			rhs[g.pidx(i, j, k)] = -(rho / dt) * div_cell(u.data(), v.data(), w.data(), g, i, j, k);
	}
	void poisson_apply_cpu(const std::vector<double>& p, std::vector<double>& Ap, MacGrid g)
	{
		Ap.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			Ap[g.pidx(i, j, k)] = (nnb_count(g, i, j, k) * p[g.pidx(i, j, k)] - nb_sum(p.data(), g, i, j, k)) / (g.h * g.h);
	}
	void poisson_residual_cpu(const std::vector<double>& p, const std::vector<double>& rhs, std::vector<double>& r, MacGrid g)
	{
		r.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			r[g.pidx(i, j, k)] = rhs[g.pidx(i, j, k)] - (nnb_count(g, i, j, k) * p[g.pidx(i, j, k)] - nb_sum(p.data(), g, i, j, k)) / (g.h * g.h);
	}
	void jacobi_smooth_cpu(std::vector<double>& p, const std::vector<double>& rhs, MacGrid g, double omega, int sweeps)
	{
		std::vector<double> pout(p.size());
		for (int s = 0; s < sweeps; ++s)
		{
			for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			{
				int c = g.pidx(i, j, k); int cnt = nnb_count(g, i, j, k);
				if (cnt == 0) { pout[c] = p[c]; continue; }
				double diag = cnt / (g.h * g.h);
				double Ap = (cnt * p[c] - nb_sum(p.data(), g, i, j, k)) / (g.h * g.h);
				pout[c] = p[c] + omega * (rhs[c] - Ap) / diag;
			}
			p = pout;
		}
	}
	void gs_band_cpu(std::vector<double>& p, const std::vector<double>& rhs, MacGrid g, int band, int sweeps, bool forward)
	{
		for (int s = 0; s < sweeps; ++s)
			for (int pass = 0; pass < 2; ++pass)
			{
				int color = forward ? pass : (1 - pass);
				for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
				{
					if (((i + j + k) & 1) != color) continue;
					bool inband = i < band || i >= g.nx - band || j < band || j >= g.ny - band || k < band || k >= g.nz - band;
					if (!inband) continue;
					int cnt = nnb_count(g, i, j, k); if (cnt == 0) continue;
					p[g.pidx(i, j, k)] = (rhs[g.pidx(i, j, k)] * g.h * g.h + nb_sum(p.data(), g, i, j, k)) / cnt;
				}
			}
	}
	void restrict_cpu(const std::vector<double>& fine, std::vector<double>& coarse, MacGrid gf, MacGrid gc)
	{
		coarse.assign(gc.p_count(), 0.0);
		for (int fk = 0; fk < gf.nz; ++fk) for (int fj = 0; fj < gf.ny; ++fj) for (int fi = 0; fi < gf.nx; ++fi)
		{
			int Cx0, Cx1, Cy0, Cy1, Cz0, Cz1; double wx0, wx1, wy0, wy1, wz0, wz1;
			axis_stencil(fi, gc.nx, Cx0, Cx1, wx0, wx1);
			axis_stencil(fj, gc.ny, Cy0, Cy1, wy0, wy1);
			axis_stencil(fk, gc.nz, Cz0, Cz1, wz0, wz1);
			double f = fine[gf.pidx(fi, fj, fk)] * 0.125;
			int cx[2] = {Cx0, Cx1}, cy[2] = {Cy0, Cy1}, cz[2] = {Cz0, Cz1};
			double wx[2] = {wx0, wx1}, wy[2] = {wy0, wy1}, wz[2] = {wz0, wz1};
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c)
				coarse[gc.pidx(cx[a], cy[b], cz[c])] += f * wx[a] * wy[b] * wz[c];
		}
	}
	void prolong_add_cpu(const std::vector<double>& coarse, std::vector<double>& fine, MacGrid gc, MacGrid gf)
	{
		for (int fk = 0; fk < gf.nz; ++fk) for (int fj = 0; fj < gf.ny; ++fj) for (int fi = 0; fi < gf.nx; ++fi)
		{
			int Cx0, Cx1, Cy0, Cy1, Cz0, Cz1; double wx0, wx1, wy0, wy1, wz0, wz1;
			axis_stencil(fi, gc.nx, Cx0, Cx1, wx0, wx1);
			axis_stencil(fj, gc.ny, Cy0, Cy1, wy0, wy1);
			axis_stencil(fk, gc.nz, Cz0, Cz1, wz0, wz1);
			int cx[2] = {Cx0, Cx1}, cy[2] = {Cy0, Cy1}, cz[2] = {Cz0, Cz1};
			double wx[2] = {wx0, wx1}, wy[2] = {wy0, wy1}, wz[2] = {wz0, wz1};
			double val = 0.0;
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c)
				val += wx[a] * wy[b] * wz[c] * coarse[gc.pidx(cx[a], cy[b], cz[c])];
			fine[gf.pidx(fi, fj, fk)] += val;
		}
	}
	void subtract_gradient_cpu(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, const std::vector<double>& p, MacGrid g, double rho, double dt)
	{
		double coef = dt / rho;
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 1; i <= g.nx - 1; ++i)
			u[g.uidx(i, j, k)] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i - 1, j, k)]) / g.h;
		for (int k = 0; k < g.nz; ++k) for (int j = 1; j <= g.ny - 1; ++j) for (int i = 0; i < g.nx; ++i)
			v[g.vidx(i, j, k)] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i, j - 1, k)]) / g.h;
		for (int k = 1; k <= g.nz - 1; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			w[g.widx(i, j, k)] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i, j, k - 1)]) / g.h;
	}
	double max_abs_divergence_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, MacGrid g)
	{
		double m = 0.0;
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			m = std::max(m, std::fabs(div_cell(u.data(), v.data(), w.data(), g, i, j, k)));
		return m;
	}
}
