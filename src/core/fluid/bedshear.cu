// bedshear.cu — see bedshear_ops.h / bedshear.h. Log-law wall-function kernels + CPU
// twins. RESEARCH §4, research/04. Every kernel shares the __host__ __device__ physics
// inlines in bedshear.h so the GPU-vs-CPU parity tests certify only launch/index/memory.
#include "core/fluid/bedshear_ops.h"

#include <cuda_runtime.h>

#include <cmath>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }
		PARACFD_HD inline int pl(MacGrid g, int i, int j) { return j * g.nx + i; }
		PARACFD_HD inline int wrap(int i, int n) { int r = i % n; return r < 0 ? r + n : r; }

		// ---- raw τ_b from the field (probe → u* → τ) --------------------------
		PARACFD_HD inline void bedshear_cell(const double* u, const double* v, MacGrid g, WallParams wp,
			int i, int j, double& taux, double& tauy, double& us)
		{
			double ks = wall_ks(wp.d50);
			double zp = wall_zp(g.h, ks);
			double uc, vc;
			bed_sample_uv(u, v, g, i, j, zp, uc, vc);
			double Up = sqrt(uc * uc + vc * vc);
			us = (Up > 1e-12) ? wall_ustar(Up, zp, wp.d50, wp.nu, wp.kappa, wp.regime, wp.cj_iters) : 0.0;
			double tau = wp.rho * us * us;
			if (Up > 1e-12) { taux = tau * uc / Up; tauy = tau * vc / Up; }
			else { taux = 0.0; tauy = 0.0; }
		}

		__global__ void k_compute(const double* u, const double* v, double* taux, double* tauy, double* ustar,
			MacGrid g, WallParams wp, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = t / g.nx;
			double tx, ty, us; bedshear_cell(u, v, g, wp, i, j, tx, ty, us);
			taux[t] = tx; tauy[t] = ty; ustar[t] = us;
		}

		// ---- 3×3 tangential box smooth (one pass, Neumann edges) --------------
		PARACFD_HD inline double smooth_cell(const double* f, MacGrid g, int i, int j)
		{
			double s = 0.0; int c = 0;
			for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di)
			{
				int ii = i + di, jj = j + dj;
				if (ii < 0 || ii >= g.nx || jj < 0 || jj >= g.ny) continue;
				s += f[pl(g, ii, jj)]; ++c;
			}
			return c ? s / c : f[pl(g, i, j)];
		}
		__global__ void k_smooth(const double* in, double* out, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = t / g.nx;
			out[t] = smooth_cell(in, g, i, j);
		}

		// ---- rate-limit + EMA ------------------------------------------------
		PARACFD_HD inline void filter_cell(double tx, double ty, double& ex, double& ey, WallParams wp, double dt, int first)
		{
			if (first) { ex = tx; ey = ty; return; }
			double mnew = sqrt(tx * tx + ty * ty);
			double mold = sqrt(ex * ex + ey * ey);
			// per-step magnitude rate-limit to [rate_lo, rate_hi]×previous (research/04 §7).
			if (mold > 1e-30 && mnew > 1e-30)
			{
				double r = mnew / mold, rc = r;
				if (rc < wp.rate_lo) rc = wp.rate_lo; else if (rc > wp.rate_hi) rc = wp.rate_hi;
				double s = rc / r; tx *= s; ty *= s;
			}
			double a = dt / wp.t_avg; if (a > 1.0) a = 1.0; if (a < 0.0) a = 0.0;
			ex = (1.0 - a) * ex + a * tx;
			ey = (1.0 - a) * ey + a * ty;
		}
		__global__ void k_filter(const double* taux, const double* tauy, double* emax, double* emay,
			WallParams wp, double dt, int first, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			double ex = emax[t], ey = emay[t];
			filter_cell(taux[t], tauy[t], ex, ey, wp, dt, first);
			emax[t] = ex; emay[t] = ey;
		}

		// ---- momentum sink on the k=0 faces ----------------------------------
		__global__ void k_sink_u(double* u, const double* emax, MacGrid g, double coef, int njx)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= njx) return;
			int i = t % (g.nx + 1), j = t / (g.nx + 1);
			int ia = wrap(i - 1, g.nx), ib = wrap(i, g.nx);
			double tf = 0.5 * (emax[pl(g, ia, j)] + emax[pl(g, ib, j)]);
			u[g.uidx(i, j, 0)] -= coef * tf;
		}
		__global__ void k_sink_v(double* v, const double* emay, MacGrid g, double coef, int nvj)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= nvj) return;
			int i = t % g.nx, j = t / g.nx; // j in [0,ny]
			int ja = wrap(j - 1, g.ny), jb = wrap(j, g.ny);
			double tf = 0.5 * (emay[pl(g, i, ja)] + emay[pl(g, i, jb)]);
			v[g.vidx(i, j, 0)] -= coef * tf;
		}
	}

	// ================= GPU launchers =====================================
	void bedshear_compute_gpu(const double* u, const double* v, double* taux, double* tauy, double* ustar,
		MacGrid g, WallParams wp)
	{
		int n = g.nx * g.ny;
		k_compute<<<gsz(n), 256>>>(u, v, taux, tauy, ustar, g, wp, n);
	}

	void bedshear_smooth_gpu(double* fx, double* fy, double* scratch, MacGrid g, int passes)
	{
		int n = g.nx * g.ny;
		for (int p = 0; p < passes; ++p)
		{
			k_smooth<<<gsz(n), 256>>>(fx, scratch, g, n);
			cudaMemcpy(fx, scratch, sizeof(double) * n, cudaMemcpyDeviceToDevice);
			k_smooth<<<gsz(n), 256>>>(fy, scratch, g, n);
			cudaMemcpy(fy, scratch, sizeof(double) * n, cudaMemcpyDeviceToDevice);
		}
	}

	void bedshear_filter_gpu(const double* taux, const double* tauy, double* emax, double* emay,
		MacGrid g, WallParams wp, double dt, int first)
	{
		int n = g.nx * g.ny;
		k_filter<<<gsz(n), 256>>>(taux, tauy, emax, emay, wp, dt, first, n);
	}

	void bedshear_apply_sink_gpu(double* u, double* v, const double* emax, const double* emay,
		MacGrid g, double rho, double dt)
	{
		double coef = dt / (rho * g.h);
		int njx = (g.nx + 1) * g.ny, nvj = g.nx * (g.ny + 1);
		k_sink_u<<<gsz(njx), 256>>>(u, emax, g, coef, njx);
		k_sink_v<<<gsz(nvj), 256>>>(v, emay, g, coef, nvj);
	}

	// ================= CPU reference twins ===============================
	void bedshear_compute_cpu(const std::vector<double>& u, const std::vector<double>& v,
		std::vector<double>& taux, std::vector<double>& tauy, std::vector<double>& ustar,
		MacGrid g, WallParams wp)
	{
		int n = g.nx * g.ny; taux.assign(n, 0.0); tauy.assign(n, 0.0); ustar.assign(n, 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			double tx, ty, us; bedshear_cell(u.data(), v.data(), g, wp, i, j, tx, ty, us);
			int p = pl(g, i, j); taux[p] = tx; tauy[p] = ty; ustar[p] = us;
		}
	}

	void bedshear_smooth_cpu(std::vector<double>& fx, std::vector<double>& fy, MacGrid g, int passes)
	{
		int n = g.nx * g.ny; std::vector<double> tmp(n);
		for (int p = 0; p < passes; ++p)
		{
			for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) tmp[pl(g, i, j)] = smooth_cell(fx.data(), g, i, j);
			fx = tmp;
			for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) tmp[pl(g, i, j)] = smooth_cell(fy.data(), g, i, j);
			fy = tmp;
		}
	}

	void bedshear_filter_cpu(const std::vector<double>& taux, const std::vector<double>& tauy,
		std::vector<double>& emax, std::vector<double>& emay, MacGrid g, WallParams wp, double dt, int first)
	{
		int n = g.nx * g.ny;
		for (int t = 0; t < n; ++t)
		{
			double ex = emax[t], ey = emay[t];
			filter_cell(taux[t], tauy[t], ex, ey, wp, dt, first);
			emax[t] = ex; emay[t] = ey;
		}
	}

	void bedshear_apply_sink_cpu(std::vector<double>& u, std::vector<double>& v,
		const std::vector<double>& emax, const std::vector<double>& emay, MacGrid g, double rho, double dt)
	{
		double coef = dt / (rho * g.h);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i)
		{
			int ia = wrap(i - 1, g.nx), ib = wrap(i, g.nx);
			double tf = 0.5 * (emax[pl(g, ia, j)] + emax[pl(g, ib, j)]);
			u[g.uidx(i, j, 0)] -= coef * tf;
		}
		for (int j = 0; j <= g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			int ja = wrap(j - 1, g.ny), jb = wrap(j, g.ny);
			double tf = 0.5 * (emay[pl(g, i, ja)] + emay[pl(g, i, jb)]);
			v[g.vidx(i, j, 0)] -= coef * tf;
		}
	}
}
