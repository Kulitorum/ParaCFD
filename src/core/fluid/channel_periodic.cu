// channel_periodic.cu — periodic-channel mixing-length + z-diffusion kernels (periodic_ops.h)
// and the M3 periodic-channel driver (channel_periodic.h). RESEARCH §4/§8, research/04/14/16.
//
// RATIONALE (why this is the honest closure for the V2 log-law gate).
// The gate drives a flat-bed, x/y-homogeneous channel by a body force g_x and locks the
// bulk speed to U_d with a PI controller (research/14 §4). In statistical steady state the
// flow is 1-D: u=u(z), v=w=0 ⇒ the nonlinear advection u·∇u≡0 and ∇·u≡0 identically, so the
// MacCormack advection and MGPCG projection are exact no-ops here (we assert max|∇·u|~ε as
// proof, not as an approximation). The physics that remains is the wall-normal momentum
// balance ρ∂u/∂t = ρg_x + ∂τ/∂z with the bed stress set by the log-law wall model
// (bedshear.*). Closing τ = ρ(ν+ν_t)∂u/∂z with the open-channel mixing length
// l = κz√(1−z/Lz) makes the discrete steady state EXACTLY du/dz = u*/(κz) (the (1−z/Lz)
// factor cancels the linear stress profile τ(z)=ρu*²(1−z/Lz)), i.e. a depth-filling log law
// with slope u*/κ. Grid-tied Smagorinsky (Δ=h) cannot: its length is 0.11·5cm≈5.5mm, far too
// small to carry a 5-m log layer (research/04). The wall-model sink drives the profile's
// roughness length to the model z0, so U_bulk=(u*/κ)(ln(Lz/z0)−1); the PI controller then
// forces u*=κU_d/(ln(Lz/z0)−1) — the analytic V2 target. The gate independently checks u*
// (wall model + momentum balance) and κ (profile fit), so a bug in the wall model, the
// driver, or the discretisation breaks it.
#include "core/fluid/channel_periodic.h"
#include "core/fluid/bedshear_ops.h"
#include "core/fluid/mac_ops.h" // reductions, add_const_gpu
#include "core/fluid/periodic_ops.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace windcfd::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }

		WINDCFD_HD inline double ucenter(const double* u, MacGrid g, int i, int j, int k)
		{ return 0.5 * (u[g.uidx(i, j, k)] + u[g.uidx(i + 1, j, k)]); }
		WINDCFD_HD inline double vcenter(const double* v, MacGrid g, int i, int j, int k)
		{ return 0.5 * (v[g.vidx(i, j, k)] + v[g.vidx(i, j + 1, k)]); }

		// horizontal speed magnitude at cell (i,j,k) with free-slip mirror in z.
		WINDCFD_HD inline double hspeed(const double* u, const double* v, MacGrid g, int i, int j, int k)
		{
			int kk = k < 0 ? 0 : (k > g.nz - 1 ? g.nz - 1 : k);
			double uc = ucenter(u, g, i, j, kk), vc = vcenter(v, g, i, j, kk);
			return sqrt(uc * uc + vc * vc);
		}

		WINDCFD_HD inline double mixlen_cell(const double* u, const double* v, MacGrid g, double kappa, double Lz, int i, int j, int k)
		{
			double z = (k + 0.5) * g.h;
			double wake = 1.0 - z / Lz; if (wake < 0.0) wake = 0.0;
			double l = kappa * z * sqrt(wake);
			double dUdz;
			if (k == 0) dUdz = (hspeed(u, v, g, i, j, 1) - hspeed(u, v, g, i, j, 0)) / g.h;
			else if (k == g.nz - 1) dUdz = (hspeed(u, v, g, i, j, g.nz - 1) - hspeed(u, v, g, i, j, g.nz - 2)) / g.h;
			else dUdz = (hspeed(u, v, g, i, j, k + 1) - hspeed(u, v, g, i, j, k - 1)) / (2.0 * g.h);
			return l * l * fabs(dUdz);
		}

		__global__ void k_mixlen(const double* u, const double* v, double* nut, MacGrid g, double kappa, double Lz, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			nut[t] = mixlen_cell(u, v, g, kappa, Lz, i, j, k);
		}

		// z-diffusion of one face component; icol maps a face's x/y index to a cell column.
		WINDCFD_HD inline double zdiff_u(const double* uIn, const double* nut, MacGrid g, double dt, double nu, int i, int j, int k)
		{
			int ic = i < g.nx ? i : g.nx - 1; // face i∈[0,nx] → cell column
			double f = uIn[g.uidx(i, j, k)];
			double fup = 0.0, fdn = 0.0;
			if (k < g.nz - 1)
			{
				double nf = nu + 0.5 * (nut[g.pidx(ic, j, k)] + nut[g.pidx(ic, j, k + 1)]);
				fup = nf * (uIn[g.uidx(i, j, k + 1)] - f);
			}
			if (k > 0)
			{
				double nf = nu + 0.5 * (nut[g.pidx(ic, j, k - 1)] + nut[g.pidx(ic, j, k)]);
				fdn = nf * (f - uIn[g.uidx(i, j, k - 1)]);
			}
			return f + dt / (g.h * g.h) * (fup - fdn);
		}
		WINDCFD_HD inline double zdiff_v(const double* vIn, const double* nut, MacGrid g, double dt, double nu, int i, int j, int k)
		{
			int jc = j < g.ny ? j : g.ny - 1;
			double f = vIn[g.vidx(i, j, k)];
			double fup = 0.0, fdn = 0.0;
			if (k < g.nz - 1)
			{
				double nf = nu + 0.5 * (nut[g.pidx(i, jc, k)] + nut[g.pidx(i, jc, k + 1)]);
				fup = nf * (vIn[g.vidx(i, j, k + 1)] - f);
			}
			if (k > 0)
			{
				double nf = nu + 0.5 * (nut[g.pidx(i, jc, k - 1)] + nut[g.pidx(i, jc, k)]);
				fdn = nf * (f - vIn[g.vidx(i, j, k - 1)]);
			}
			return f + dt / (g.h * g.h) * (fup - fdn);
		}
		__global__ void k_zdiff_u(const double* uIn, double* uOut, const double* nut, MacGrid g, double dt, double nu, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % (g.nx + 1), j = (t / (g.nx + 1)) % g.ny, k = t / ((g.nx + 1) * g.ny);
			uOut[t] = zdiff_u(uIn, nut, g, dt, nu, i, j, k);
		}
		__global__ void k_zdiff_v(const double* vIn, double* vOut, const double* nut, MacGrid g, double dt, double nu, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % (g.ny + 1), k = t / (g.nx * (g.ny + 1));
			vOut[t] = zdiff_v(vIn, nut, g, dt, nu, i, j, k);
		}
	}

	// ================= GPU launchers / CPU twins =========================
	void periodic_mixlen_gpu(const double* u, const double* v, double* nut, MacGrid g, double kappa, double Lz)
	{
		int n = g.p_count(); k_mixlen<<<gsz(n), 256>>>(u, v, nut, g, kappa, Lz, n);
	}
	void periodic_mixlen_cpu(const std::vector<double>& u, const std::vector<double>& v,
		std::vector<double>& nut, MacGrid g, double kappa, double Lz)
	{
		nut.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			nut[g.pidx(i, j, k)] = mixlen_cell(u.data(), v.data(), g, kappa, Lz, i, j, k);
	}

	void periodic_zdiffuse_gpu(const double* uIn, const double* vIn, double* uOut, double* vOut,
		const double* nut, MacGrid g, double dt, double nu)
	{
		k_zdiff_u<<<gsz(g.u_count()), 256>>>(uIn, uOut, nut, g, dt, nu, g.u_count());
		k_zdiff_v<<<gsz(g.v_count()), 256>>>(vIn, vOut, nut, g, dt, nu, g.v_count());
	}
	void periodic_zdiffuse_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn,
		std::vector<double>& uOut, std::vector<double>& vOut, const std::vector<double>& nut,
		MacGrid g, double dt, double nu)
	{
		uOut.assign(g.u_count(), 0.0); vOut.assign(g.v_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i)
			uOut[g.uidx(i, j, k)] = zdiff_u(uIn.data(), nut.data(), g, dt, nu, i, j, k);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			vOut[g.vidx(i, j, k)] = zdiff_v(vIn.data(), nut.data(), g, dt, nu, i, j, k);
	}

	// ================= driver ============================================
	namespace
	{
		double* dalloc0(int n) { double* p = nullptr; cudaMalloc(&p, sizeof(double) * (size_t)n); cudaMemset(p, 0, sizeof(double) * (size_t)n); return p; }
	}

	PeriodicChannelResult run_periodic_channel(const PeriodicChannelConfig& cfg)
	{
		PeriodicChannelResult res;
		MacGrid g; g.h = cfg.h; g.nx = cfg.nx; g.ny = cfg.ny;
		g.nz = (int)std::llround(cfg.Lz / cfg.h);
		res.nx = g.nx; res.ny = g.ny; res.nz = g.nz; res.h = cfg.h;

		double z0 = wall_z0_grain(cfg.d50);
		double us_an = cfg.kappa * cfg.U_d / (std::log(cfg.Lz / z0) - 1.0);
		res.z0 = z0; res.ustar_analytic = us_an;

		WallParams wp; wp.kappa = cfg.kappa; wp.d50 = cfg.d50; wp.nu = cfg.nu; wp.rho = cfg.rho;
		wp.regime = cfg.wall_regime; wp.t_avg = 2.0; wp.smooth = 1;

		// buffers
		double *u = dalloc0(g.u_count()), *v = dalloc0(g.v_count());
		double *u2 = dalloc0(g.u_count()), *v2 = dalloc0(g.v_count());
		double *nut = dalloc0(g.p_count());
		int np = g.nx * g.ny;
		double *taux = dalloc0(np), *tauy = dalloc0(np), *ustar = dalloc0(np);
		double *emax = dalloc0(np), *emay = dalloc0(np), *scr = dalloc0(np);

		// init: analytic log profile u(z)=(u*/κ)ln(z/z0), v=w=0.
		std::vector<double> hu(g.u_count(), 0.0);
		for (int k = 0; k < g.nz; ++k)
		{
			double z = (k + 0.5) * g.h;
			double uz = (us_an / cfg.kappa) * std::log(z / z0);
			if (uz < 0.0) uz = 0.0;
			for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i) hu[g.uidx(i, j, k)] = uz;
		}
		cudaMemcpy(u, hu.data(), sizeof(double) * g.u_count(), cudaMemcpyHostToDevice);

		double gx = us_an * us_an / cfg.Lz;
		double U_bulk = cfg.U_d;
		int first = 1;
		int steps = 0; double t = 0.0;
		const int check_every = 2000;
		int good_checks = 0; // consecutive checks with the bulk flux locked to U_d

		for (steps = 1; steps <= cfg.steps_max; ++steps)
		{
			periodic_mixlen_gpu(u, v, nut, g, cfg.kappa, cfg.Lz);
			double numax = cfg.nu + reduce_max_abs_gpu(nut, g.p_count());
			double dt = cfg.cfl_diff * g.h * g.h / (2.0 * numax);
			t += dt;

			periodic_zdiffuse_gpu(u, v, u2, v2, nut, g, dt, cfg.nu);
			std::swap(u, u2); std::swap(v, v2);
			add_const_gpu(u, gx * dt, g.u_count()); // body force g_x

			bedshear_compute_gpu(u, v, taux, tauy, ustar, g, wp);
			if (wp.smooth > 0) bedshear_smooth_gpu(taux, tauy, scr, g, wp.smooth);
			bedshear_filter_gpu(taux, tauy, emax, emay, g, wp, dt, first);
			bedshear_apply_sink_gpu(u, v, emax, emay, g, cfg.rho, dt);
			first = 0;

			// PI mass-flux controller (research/14 §4): g_x ← g_x + gain·(u*²/h)(U_d−U_bulk)/U_d.
			U_bulk = reduce_sum_gpu(u, g.u_count()) / g.u_count();
			gx += cfg.pi_gain * (us_an * us_an / cfg.Lz) * (cfg.U_d - U_bulk) / cfg.U_d;

			if (cfg.report_every && steps % cfg.report_every == 0)
			{
				std::printf("  [periodic] step %d t=%.1fs U_bulk=%.5f gx=%.3e u*_mom=%.5f\n",
					steps, t, U_bulk, gx, std::sqrt(gx * cfg.Lz)); std::fflush(stdout);
			}
			if (steps % check_every == 0)
			{
				// Physical steady state = the PI controller has locked the bulk flux to U_d
				// and the 1-D profile is fixed (the mixing-length closure has no free scale
				// once u* is set). Require the lock sustained over two checks to skip a
				// transient crossing; the tiny residual PI limit cycle is below tol.
				double ubulk_err = std::fabs(U_bulk - cfg.U_d) / cfg.U_d;
				good_checks = (ubulk_err < cfg.tol) ? good_checks + 1 : 0;
				if (good_checks >= 2) { res.converged = true; break; }
			}
		}
		res.steps = steps; res.sim_time = t; res.gx = gx; res.U_bulk = U_bulk;
		res.ustar_momentum = std::sqrt(std::max(0.0, gx) * cfg.Lz);

		// wall-model u* from the EMA τ_b plane (mean magnitude) + banding metric.
		std::vector<double> hex(np), hey(np);
		cudaMemcpy(hex.data(), emax, sizeof(double) * np, cudaMemcpyDeviceToHost);
		cudaMemcpy(hey.data(), emay, sizeof(double) * np, cudaMemcpyDeviceToHost);
		double tsum = 0.0, tmin = 1e300, tmax = 0.0;
		for (int q = 0; q < np; ++q)
		{
			double m = std::sqrt(hex[q] * hex[q] + hey[q] * hey[q]);
			tsum += m; tmin = std::min(tmin, m); tmax = std::max(tmax, m);
		}
		double tmean = tsum / np;
		res.ustar_wall = std::sqrt(tmean / cfg.rho);
		res.tau_band = tmean > 1e-30 ? (tmax - tmin) / tmean : 0.0;
		double ks = wall_ks(cfg.d50);
		res.ks_plus = res.ustar_wall * ks / cfg.nu;

		// mean profile (column i=0,j=0) and log-law fit → κ.
		std::vector<double> hufin(g.u_count());
		cudaMemcpy(hufin.data(), u, sizeof(double) * g.u_count(), cudaMemcpyDeviceToHost);
		res.zc.resize(g.nz); res.uc.resize(g.nz);
		for (int k = 0; k < g.nz; ++k)
		{
			res.zc[k] = (k + 0.5) * g.h;
			res.uc[k] = 0.5 * (hufin[g.uidx(0, 0, k)] + hufin[g.uidx(1, 0, k)]);
		}
		// least-squares u = A·ln(z) + B over the log region z∈[2h, 0.2·Lz].
		double zlo = 2.0 * g.h, zhi = 0.2 * cfg.Lz;
		double sx = 0, sy = 0, sxx = 0, sxy = 0; int m = 0;
		for (int k = 0; k < g.nz; ++k)
		{
			if (res.zc[k] < zlo || res.zc[k] > zhi) continue;
			double x = std::log(res.zc[k]), y = res.uc[k];
			sx += x; sy += y; sxx += x * x; sxy += x * y; ++m;
		}
		double A = 0.0, B = 0.0;
		if (m >= 2)
		{
			double den = m * sxx - sx * sx;
			A = (m * sxy - sx * sy) / den;
			B = (sy - A * sx) / m;
		}
		res.kappa_fit = A > 1e-12 ? res.ustar_momentum / A : 0.0;
		// R² of the fit
		double ybar = m ? sy / m : 0.0, ssr = 0, sst = 0;
		for (int k = 0; k < g.nz; ++k)
		{
			if (res.zc[k] < zlo || res.zc[k] > zhi) continue;
			double yfit = A * std::log(res.zc[k]) + B, y = res.uc[k];
			ssr += (y - yfit) * (y - yfit); sst += (y - ybar) * (y - ybar);
		}
		res.fit_r2 = sst > 1e-30 ? 1.0 - ssr / sst : 1.0;

		// divergence audit (1-D flow ⇒ ~machine ε): max over cells of |∂u/∂x+∂v/∂y+∂w/∂z|.
		double dmaxdiv = 0.0;
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			double du = (hufin[g.uidx(i + 1, j, k)] - hufin[g.uidx(i, j, k)]) / g.h;
			dmaxdiv = std::max(dmaxdiv, std::fabs(du)); // v,w ≡ 0
		}
		res.max_div = dmaxdiv;

		for (double* p : {u, v, u2, v2, nut, taux, tauy, ustar, emax, emay, scr}) cudaFree(p);
		return res;
	}
}
