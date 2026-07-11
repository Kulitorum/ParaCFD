// sed_validation.cu — M4 gate drivers: still-water settling column (V5) and equilibrium Rouse
// profile (V6). Both drive the shipped suspended-sediment kernels on the GPU. RESEARCH §6.1–6.2,
// research/02 §7, research/11 §6–7. Analytic references and the honest reading of RESEARCH §10's
// exponential are documented at each driver.
#include "core/sediment/sed_validation.h"
#include "core/sediment/suspended.h"
#include "core/sediment/sed_physics.h"
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h" // reduce_sum_gpu, reduce_max_abs_gpu

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace scour::core
{
	namespace
	{
		double* dalloc0(size_t n) { double* p = nullptr; cudaMalloc(&p, sizeof(double) * n); cudaMemset(p, 0, sizeof(double) * n); return p; }

		// least-squares slope/intercept of y = A·x + B, plus R².
		void lsq(const std::vector<double>& x, const std::vector<double>& y, double& A, double& B, double& r2)
		{
			int m = (int)x.size();
			double sx = 0, sy = 0, sxx = 0, sxy = 0;
			for (int q = 0; q < m; ++q) { sx += x[q]; sy += y[q]; sxx += x[q] * x[q]; sxy += x[q] * y[q]; }
			double den = m * sxx - sx * sx;
			A = den != 0 ? (m * sxy - sx * sy) / den : 0.0;
			B = m ? (sy - A * sx) / m : 0.0;
			double ybar = m ? sy / m : 0.0, ssr = 0, sst = 0;
			for (int q = 0; q < m; ++q) { double yf = A * x[q] + B; ssr += (y[q] - yf) * (y[q] - yf); sst += (y[q] - ybar) * (y[q] - ybar); }
			r2 = sst > 1e-300 ? 1.0 - ssr / sst : 1.0;
		}
	}

	// ================= V5: still-water settling column ================================
	SettlingColumnResult run_settling_column(const SettlingColumnConfig& cfg)
	{
		SettlingColumnResult res;
		MacGrid g; g.h = cfg.h; g.nx = cfg.nx; g.ny = cfg.ny; g.nz = (int)std::llround(cfg.Lz / cfg.h);
		res.nx = g.nx; res.ny = g.ny; res.nz = g.nz;

		double ws0 = settling_ws(cfg.d50, cfg.rho, cfg.rho_s, cfg.nu);
		res.ws = ws0;
		double dt = cfg.cfl * g.h / ws0; res.dt = dt;
		double Vcell = g.h * g.h * g.h;

		int np = g.p_count();
		// initial Gaussian slab, homogeneous in x/y.
		std::vector<double> hc(np, 0.0);
		for (int k = 0; k < g.nz; ++k)
		{
			double z = (k + 0.5) * g.h;
			double val = cfg.c0 * std::exp(-((z - cfg.z_center) * (z - cfg.z_center)) / (2.0 * cfg.sigma * cfg.sigma));
			for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) hc[g.pidx(i, j, k)] = val;
		}
		double* c = dalloc0(np); cudaMemcpy(c, hc.data(), sizeof(double) * np, cudaMemcpyHostToDevice);
		double* c2 = dalloc0(np);
		double* scratch = dalloc0(np);
		double* ws = dalloc0(np); { std::vector<double> hws(np, ws0); cudaMemcpy(ws, hws.data(), sizeof(double) * np, cudaMemcpyHostToDevice); }
		double* uz = dalloc0(g.u_count()), *vz = dalloc0(g.v_count()), *wz = dalloc0(g.w_count()); // zero fluid velocity
		int ncol = g.nx * g.ny;
		double* bedflux = dalloc0(ncol);

		SedBedParams bp; bp.d50 = cfg.d50; bp.rho = cfg.rho; bp.rho_s = cfg.rho_s; bp.nu = cfg.nu;
		bp.ws0 = ws0; bp.hindered = 0; // constant w_s ⇒ exact-translate reference (hindered unit-tested separately)

		double M0 = reduce_sum_gpu(c, np) * Vcell;
		double M_bed = 0.0; // cumulative sand mass delivered to the bed [volume]
		double mass_err_max = 0.0;

		double t = 0.0;
		bool checked = false;
		double t_check = cfg.travel_check / ws0;
		double t_end = 1.4 * cfg.z_center / ws0; // long enough for the slab to reach and deposit at the bed
		int steps = 0;
		while (t < t_end)
		{
			// (1) settle: conservative flux-form advection with vz = w_fluid − w_s = −w_s. Zero
			// advective flux through the bed/surface ⇒ sand piles conservatively at the bed cell,
			// where the deposition sink drains it — no clamp, mass conserved to machine precision.
			suspended_advect_cons_gpu(c, uz, vz, wz, ws, c2, scratch, g, dt);
			std::swap(c, c2);
			// (2) bed exchange (deposition; still water ⇒ pickup E=0). Accumulate the applied flux.
			suspended_bed_exchange_gpu(c, nullptr, bedflux, g, bp, dt);
			double db = reduce_sum_gpu(bedflux, ncol);       // Σ Δ(c·h) per column [m]
			M_bed += -db * g.h * g.h;                        // deposited (Δsuspended<0 ⇒ M_bed>0)
			t += dt; ++steps;

			double M_susp = reduce_sum_gpu(c, np) * Vcell;
			double err = std::fabs(M_susp + M_bed - M0) / M0; // total sand (suspended + deposited)
			if (err > mass_err_max) mass_err_max = err;

			if (!checked && t >= t_check)
			{
				checked = true; res.t_check = t;
				std::vector<double> cc(np); cudaMemcpy(cc.data(), c, sizeof(double) * np, cudaMemcpyDeviceToHost);
				double z_shift = cfg.z_center - ws0 * t;
				double num = 0, den = 0, cmsum = 0, wsum = 0;
				res.zc.resize(g.nz); res.c_num.resize(g.nz); res.c_exact.resize(g.nz);
				for (int k = 0; k < g.nz; ++k)
				{
					double z = (k + 0.5) * g.h;
					double cn = cc[g.pidx(0, 0, k)];
					double ce = cfg.c0 * std::exp(-((z - z_shift) * (z - z_shift)) / (2.0 * cfg.sigma * cfg.sigma));
					res.zc[k] = z; res.c_num[k] = cn; res.c_exact[k] = ce;
					num += (cn - ce) * (cn - ce); den += ce * ce;
					cmsum += cn * z; wsum += cn;
				}
				res.l2_translate = den > 0 ? std::sqrt(num / den) : 0.0;
				double z_cm = wsum > 0 ? cmsum / wsum : cfg.z_center;
				res.v_centroid = (cfg.z_center - z_cm) / t;
				res.v_centroid_err = std::fabs(res.v_centroid - ws0) / ws0;
			}
		}
		double M_susp_end = reduce_sum_gpu(c, np) * Vcell;
		res.mass_err = mass_err_max;
		res.susp_frac_end = M_susp_end / M0;

		for (double* p : { c, c2, scratch, ws, uz, vz, wz, bedflux }) cudaFree(p);
		return res;
	}

	// ================= V6: equilibrium Rouse profile ==================================
	RouseChannelResult run_rouse_channel(const RouseChannelConfig& cfg)
	{
		RouseChannelResult res;
		double a = cfg.a_frac * cfg.h_dom;
		double H = cfg.h_dom - a;
		int nz = cfg.nz;
		double dz = H / nz;
		MacGrid g; g.h = dz; g.nx = cfg.nx; g.ny = cfg.ny; g.nz = nz;
		res.nz = nz;

		double ws0 = settling_ws(cfg.d50, cfg.rho, cfg.rho_s, cfg.nu);
		res.ws = ws0;
		res.R_target = ws0 / (cfg.beta * cfg.kappa * cfg.u_star);

		int np = g.p_count();
		int ncol = g.nx * g.ny;

		// diffusivity D_c(z) = ε_s = β·κ·u*·z·(1−z/h_dom), z = physical height above the real bed.
		std::vector<double> hDc(np, 0.0), hws(np, ws0);
		for (int k = 0; k < g.nz; ++k)
		{
			double z = a + (k + 0.5) * dz;
			double eps = cfg.beta * cfg.kappa * cfg.u_star * z * (1.0 - z / cfg.h_dom);
			if (eps < 0) eps = 0;
			for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) hDc[g.pidx(i, j, k)] = eps;
		}
		double maxDc = 0; for (double d : hDc) maxDc = std::max(maxDc, d);
		double dt_diff = 0.5 * dz * dz / std::max(maxDc, 1e-30);
		double dt_adv = 0.4 * dz / ws0;
		double dt = std::min(dt_diff, dt_adv); res.dt = dt;

		double* c = dalloc0(np);
		double* c2 = dalloc0(np);
		double* scratch = dalloc0(np);
		double* Dc = dalloc0(np); cudaMemcpy(Dc, hDc.data(), sizeof(double) * np, cudaMemcpyHostToDevice);
		double* ws = dalloc0(np); cudaMemcpy(ws, hws.data(), sizeof(double) * np, cudaMemcpyHostToDevice);
		double* uz = dalloc0(g.u_count()), *vz = dalloc0(g.v_count()), *wz = dalloc0(g.w_count()); // zero fluid velocity

		// Seed a NON-uniform profile (Rouse shape with a deliberately WRONG exponent) so the scheme
		// must relax to the true equilibrium — a uniform field is a spurious fixed point of the
		// settling advection (its flux divergence vanishes), so seeding uniform would never develop
		// the profile. The recovered R must come from the settling/diffusion balance, not the seed.
		{
			std::vector<double> hc(np, cfg.c_a);
			double R_seed = 0.5 * res.R_target; // wrong on purpose
			for (int k = 0; k < g.nz; ++k)
			{
				double z = a + (k + 0.5) * dz;
				double val = cfg.c_a * std::pow(((cfg.h_dom - z) / z) * (a / (cfg.h_dom - a)), R_seed);
				for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) hc[g.pidx(i, j, k)] = val;
			}
			cudaMemcpy(c, hc.data(), sizeof(double) * np, cudaMemcpyHostToDevice);
		}
		std::vector<double> cA_plane(ncol, cfg.c_a); // k=0 plane is contiguous [0,ncol)

		// Conservative advection carries ZERO settling flux across the surface (no sand above) so the
		// top depletes and the gradient develops (a flux-form scheme, unlike SL, has no spurious
		// uniform fixed point). Diffusion: zero-flux (Neumann) surface, Dirichlet c_a at the bed z=a.
		ScalarBC sbc_dif; sbc_dif.zmin = SC_DIRICHLET; sbc_dif.v_zmin = cfg.c_a; // zmax Neumann

		std::vector<double> prev(np, 0.0), cur(np, 0.0);
		double t = 0.0; int steps = 0; bool converged = false;
		const int check_every = 500;
		for (steps = 1; steps <= cfg.steps_max; ++steps)
		{
			suspended_advect_cons_gpu(c, uz, vz, wz, ws, c2, scratch, g, dt);
			std::swap(c, c2);
			suspended_diffuse_gpu(c, c2, Dc, g, sbc_dif, dt);
			std::swap(c, c2);
			cudaMemcpy(c, cA_plane.data(), sizeof(double) * ncol, cudaMemcpyHostToDevice); // pin bed cell (z=a)
			t += dt;

			if (steps % check_every == 0)
			{
				cudaMemcpy(cur.data(), c, sizeof(double) * np, cudaMemcpyDeviceToHost);
				double dmax = 0, cmax = 0;
				for (int q = 0; q < np; ++q) { dmax = std::max(dmax, std::fabs(cur[q] - prev[q])); cmax = std::max(cmax, std::fabs(cur[q])); }
				prev = cur;
				if (cmax > 0 && dmax / cmax < cfg.tol && steps > check_every) { converged = true; break; }
			}
		}
		res.steps = steps; res.sim_time = t; res.converged = converged;

		std::vector<double> hc(np); cudaMemcpy(hc.data(), c, sizeof(double) * np, cudaMemcpyDeviceToHost);
		res.zc.resize(g.nz); res.c_num.resize(g.nz); res.c_rouse.resize(g.nz);
		std::vector<double> fx, fy;
		double zlo = a + 2.0 * dz, zhi = 0.5 * cfg.h_dom; // fit window: clear of the pinned bed and the surface
		for (int k = 0; k < g.nz; ++k)
		{
			double z = a + (k + 0.5) * dz;
			double cn = hc[g.pidx(0, 0, k)];
			double rouse = cfg.c_a * std::pow(((cfg.h_dom - z) / z) * (a / (cfg.h_dom - a)), res.R_target);
			res.zc[k] = z; res.c_num[k] = cn; res.c_rouse[k] = rouse;
			if (z >= zlo && z <= zhi && cn > 1e-12)
			{
				fx.push_back(std::log((cfg.h_dom - z) / z));
				fy.push_back(std::log(cn));
			}
		}
		double A = 0, B = 0, r2 = 0;
		if (fx.size() >= 2) lsq(fx, fy, A, B, r2);
		res.R_fit = A; res.fit_r2 = r2;
		res.R_err = res.R_target > 0 ? std::fabs(res.R_fit - res.R_target) / res.R_target : 0.0;

		for (double* p : { c, c2, scratch, Dc, ws, uz, vz, wz }) cudaFree(p);
		return res;
	}
}
