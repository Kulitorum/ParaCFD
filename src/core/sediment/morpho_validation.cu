// morpho_validation.cu — M5 gate drivers (a-i deposition, a-ii prescribed bedload, b sand-pile,
// c three-reservoir mass conservation, d threshold of motion). Each drives the shipped kernels
// (bedstate.*, avalanche.*, suspended.*) on the GPU. RESEARCH §5, §7; research/13, /03, /01, /16.
#include "core/sediment/morpho_validation.h"
#include "core/sediment/bedstate.h"
#include "core/sediment/avalanche.h"
#include "core/sediment/suspended.h"
#include "core/sediment/sed_physics.h"
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h" // reduce_sum_gpu

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace scour::core
{
	namespace
	{
		double* dalloc0(size_t n) { double* p = nullptr; cudaMalloc(&p, sizeof(double) * n); cudaMemset(p, 0, sizeof(double) * n); return p; }
		double* dfill(size_t n, double val) { double* p = dalloc0(n); std::vector<double> h(n, val); cudaMemcpy(p, h.data(), sizeof(double) * n, cudaMemcpyHostToDevice); return p; }
		double* dupload(const std::vector<double>& h) { double* p = dalloc0(h.size()); cudaMemcpy(p, h.data(), sizeof(double) * h.size(), cudaMemcpyHostToDevice); return p; }
		const double DEG = 3.14159265358979323846 / 180.0;
	}

	// ================= (a-i) still-water deposition column ============================
	DepositionColumnResult run_deposition_column(const DepositionColumnConfig& cfg)
	{
		DepositionColumnResult res;
		MacGrid g; g.h = cfg.h; g.nx = cfg.nx; g.ny = cfg.ny; g.nz = (int)std::llround(cfg.Lz / cfg.h);
		res.nx = g.nx; res.ny = g.ny; res.nz = g.nz;
		double cpack = SED_CPACK, Vcell = g.h * g.h * g.h;
		double ws0 = settling_ws(cfg.d50, cfg.rho, cfg.rho_s, cfg.nu); res.ws = ws0;
		double dt = cfg.cfl * g.h / ws0; res.dt = dt;
		int np = g.p_count(), ncol = g.nx * g.ny;

		double* c = dfill(np, cfg.c0);
		double* c2 = dalloc0(np), *scratch = dalloc0(np);
		double* ws = dfill(np, ws0);
		double* uz = dalloc0(g.u_count()), *vz = dalloc0(g.v_count()), *wz = dalloc0(g.w_count());
		double* G = dalloc0(ncol);
		double* taub = dalloc0(ncol);           // still water: no shear
		double* beta = dalloc0(ncol), *upx = dfill(ncol, 1.0), *upy = dalloc0(ncol);
		double* divz = dalloc0(ncol);

		MorphoParams p; p.d50 = cfg.d50; p.rho = cfg.rho; p.rho_s = cfg.rho_s; p.nu = cfg.nu; p.ws0 = ws0;
		p.cpack = cpack; p.morfac = 1.0; p.hindered = 0; p.deposition_on = 1; p.erosion_on = 0; p.bedload_on = 0;

		double V0 = reduce_sum_gpu(c, np) * Vcell;
		double t1 = 0.35 * cfg.Lz / ws0, t2 = 0.75 * cfg.Lz / ws0, t_end = 1.3 * cfg.Lz / ws0;
		double G_t1 = 0, G_t2 = 0, ta1 = 0, ta2 = 0; bool got1 = false, got2 = false;
		double mass_err = 0.0, t = 0.0;
		while (t < t_end)
		{
			suspended_advect_cons_gpu(c, uz, vz, wz, ws, c2, scratch, g, dt); std::swap(c, c2);
			morpho_exner_gpu(G, c, taub, taub, beta, upx, upy, divz, nullptr, g, p, dt);
			t += dt;
			double Vb = reduce_sum_gpu(G, ncol) * g.h * g.h;
			double Vs = reduce_sum_gpu(c, np) * Vcell;
			double e = std::fabs(Vb + Vs - V0) / V0; if (e > mass_err) mass_err = e;
			if (!got1 && t >= t1) { got1 = true; ta1 = t; res.t_check = t; G_t1 = reduce_sum_gpu(G, ncol) / ncol; }
			if (!got2 && t >= t2) { got2 = true; ta2 = t; G_t2 = reduce_sum_gpu(G, ncol) / ncol; }
		}
		double rate = (G_t2 - G_t1) / (ta2 - ta1);        // measured dG/dt at the actual sample times
		res.G_num = G_t2; res.G_expected = ws0 * cfg.c0 * (ta2 - ta1) + G_t1; // for reporting
		res.rate_err = std::fabs(rate - ws0 * cfg.c0) / (ws0 * cfg.c0);
		res.zb_rate = ws0 * cfg.c0 / cpack;
		res.mass_err = mass_err;
		for (double* q : { c, c2, scratch, ws, uz, vz, wz, G, taub, beta, upx, upy, divz }) cudaFree(q);
		return res;
	}

	// ================= (a-ii) prescribed sinusoidal bedload vs 1D Exner ===============
	BedloadExnerResult run_bedload_exner(const BedloadExnerConfig& cfg)
	{
		BedloadExnerResult res;
		MacGrid g; g.nx = cfg.nx; g.ny = 1; g.nz = 8; g.h = cfg.L / cfg.nx;
		double cpack = 1.0 - cfg.porosity;
		int ncol = g.nx;
		double k = 2.0 * 3.14159265358979323846 / cfg.L;
		double dt = cfg.dt;
		if (dt <= 0.0) dt = 0.2 * g.h * cpack * cfg.L / (cfg.updates * cfg.q0 * 2.0 * 3.14159265358979323846); // Δz_b,peak≈0.2h
		res.nx = g.nx; res.updates = cfg.updates; res.dt = dt;

		std::vector<double> hqx(ncol), hqy(ncol, 0.0);
		for (int i = 0; i < g.nx; ++i) { double x = (i + 0.5) * g.h; hqx[i] = cfg.q0 * std::sin(k * x); }
		double* qx = dupload(hqx), *qy = dupload(hqy);
		double* dirx = dfill(ncol, 1.0), *diry = dalloc0(ncol); // uniform +x transport (1D Exner)
		double* divq = dalloc0(ncol);
		bedload_div_gpu(qx, qy, dirx, diry, divq, g, /*periodic*/ 1);

		double* G = dfill(ncol, 4.0 * g.h * cpack); // z_b0 = 4h (headroom)
		double* taub = dalloc0(ncol);
		double* beta = dalloc0(ncol), *upx = dfill(ncol, 1.0), *upy = dalloc0(ncol);
		std::vector<double> G0(ncol); cudaMemcpy(G0.data(), G, sizeof(double) * ncol, cudaMemcpyDeviceToHost);

		MorphoParams p; p.cpack = cpack; p.morfac = 1.0; p.deposition_on = 0; p.erosion_on = 0; p.bedload_on = 1;
		for (int s = 0; s < cfg.updates; ++s)
			morpho_exner_gpu(G, nullptr, taub, taub, beta, upx, upy, divq, nullptr, g, p, dt);

		std::vector<double> Ge(ncol); cudaMemcpy(Ge.data(), G, sizeof(double) * ncol, cudaMemcpyDeviceToHost);
		double num = 0, den = 0, dzmax = 0;
		for (int i = 0; i < g.nx; ++i)
		{
			double x = (i + 0.5) * g.h;
			double dz_num = (Ge[i] - G0[i]) / cpack;
			double dz_exact = -(cfg.updates * dt / cpack) * cfg.q0 * k * std::cos(k * x);
			num += (dz_num - dz_exact) * (dz_num - dz_exact); den += dz_exact * dz_exact;
			dzmax = std::max(dzmax, std::fabs(dz_num));
		}
		res.l2_err = den > 0 ? std::sqrt(num / den) : 0.0; res.dz_max = dzmax;
		for (double* q : { qx, qy, dirx, diry, divq, G, taub, beta, upx, upy }) cudaFree(q);
		return res;
	}

	// ================= (b) sand-pile avalanche relaxation ============================
	SandpileResult run_sandpile(const SandpileConfig& cfg)
	{
		SandpileResult res;
		MacGrid g; g.nx = cfg.nx; g.ny = cfg.ny; g.nz = 1; g.h = cfg.h;
		int ncol = g.nx * g.ny;
		double cx = (g.nx - 1) * 0.5, cy = (g.ny - 1) * 0.5;
		double slope0 = std::tan(cfg.init_slope_deg * DEG);
		std::vector<double> z(ncol, 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			double r = std::sqrt((i - cx) * (i - cx) + (j - cy) * (j - cy)) * g.h;
			double zz = cfg.peak - slope0 * r; if (zz < 0) zz = 0;
			z[j * g.nx + i] = zz;
		}
		res.init_max_slope_deg = std::atan(avalanche_max_tanslope(z, g, 0)) / DEG;
		double m0 = 0; for (double v : z) m0 += v;

		AvalancheParams ap; ap.tan_trigger = std::tan(cfg.phi_trigger_deg * DEG);
		ap.tan_repose = std::tan(cfg.phi_repose_deg * DEG); ap.relax = cfg.relax; ap.periodic = 0;

		double* za = dupload(z), *zb = dalloc0(ncol);
		std::vector<double> cur = z;
		int sweeps = 0; bool converged = false;
		// hysteresis activation: only avalanche if the pile is steeper than the trigger angle.
		if (avalanche_max_tanslope(z, g, 0) > ap.tan_trigger)
		{
			for (sweeps = 1; sweeps <= cfg.max_sweeps; ++sweeps)
			{
				avalanche_sweep_gpu(za, zb, g, ap); std::swap(za, zb);
				if (sweeps % 200 == 0 || sweeps == cfg.max_sweeps)
				{
					cudaMemcpy(cur.data(), za, sizeof(double) * ncol, cudaMemcpyDeviceToHost);
					if (avalanche_max_tanslope(cur, g, 0) <= ap.tan_repose * 1.02) { converged = true; break; }
				}
			}
		}
		cudaMemcpy(cur.data(), za, sizeof(double) * ncol, cudaMemcpyDeviceToHost);
		res.sweeps = sweeps; res.converged = converged;
		res.final_max_slope_deg = std::atan(avalanche_max_tanslope(cur, g, 0)) / DEG;
		double m1 = 0; for (double v : cur) m1 += v;
		res.mass_err = m0 > 0 ? std::fabs(m1 - m0) / m0 : std::fabs(m1 - m0);
		cudaFree(za); cudaFree(zb);
		return res;
	}

	// ================= (c) three-reservoir mass conservation =========================
	MassConservationResult run_mass_conservation(const MassConservationConfig& cfg)
	{
		MassConservationResult res;
		MacGrid g; g.nx = cfg.nx; g.ny = cfg.ny; g.nz = cfg.nz; g.h = cfg.h;
		double cpack = SED_CPACK, Vcell = g.h * g.h * g.h;
		int np = g.p_count(), ncol = g.nx * g.ny;
		double ws0 = settling_ws(cfg.d50, cfg.rho, cfg.rho_s, cfg.nu);
		double dt = cfg.cfl * g.h / ws0; res.dt = dt; res.steps = cfg.steps;

		double tcr = tau_cr(cfg.d50, cfg.rho, cfg.rho_s, cfg.nu);
		double tau_peak = cfg.tau_peak > 0 ? cfg.tau_peak : 2.0 * tcr; // live-bed

		// steady grain-skin shear in +x, varying across x (some columns below threshold → deposit).
		std::vector<double> htx(ncol), hty(ncol, 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			htx[j * g.nx + i] = tau_peak * (0.2 + 0.8 * (double)i / (g.nx - 1));
		double* taubx = dupload(htx), *tauby = dupload(hty);

		// initial thin bed (z_b≈3h ⇒ k_bed=3) + suspended c0 in the fluid above.
		double* G = dfill(ncol, 3.0 * g.h * cpack);
		std::vector<double> hc(np, 0.0);
		for (int k = 4; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) hc[g.pidx(i, j, k)] = cfg.c0;
		double* c = dupload(hc), *c2 = dalloc0(np), *scratch = dalloc0(np);
		double* ws = dfill(np, ws0);
		double* uz = dalloc0(g.u_count()), *vz = dalloc0(g.v_count()), *wz = dalloc0(g.w_count());
		double* fpack = dalloc0(np), *Ft = dalloc0(np);
		double* nbx = dalloc0(np), *nby = dalloc0(np), *nbz = dalloc0(np), *Ab = dalloc0(np);
		double* beta = dalloc0(ncol), *upx = dfill(ncol, 1.0), *upy = dalloc0(ncol);
		double* qx = dalloc0(ncol), *qy = dalloc0(ncol), *divq = dalloc0(ncol);

		MorphoParams p; p.d50 = cfg.d50; p.rho = cfg.rho; p.rho_s = cfg.rho_s; p.nu = cfg.nu; p.ws0 = ws0;
		p.cpack = cpack; p.morfac = 1.0; p.hindered = 1; p.deposition_on = 1; p.erosion_on = 1; p.bedload_on = 1;

		AvalancheParams ap; ap.tan_repose = cpack * std::tan(30.0 * DEG); ap.tan_trigger = cpack * std::tan(32.0 * DEG);
		ap.relax = 0.1; ap.periodic = 1; // avalanche directly on G (threshold scaled by c_pack)
		double* Gb = dalloc0(ncol);

		double V0 = reduce_sum_gpu(G, ncol) * g.h * g.h + reduce_sum_gpu(c, np) * Vcell; res.V0 = V0;
		double mass_err = 0.0;
		for (int s = 0; s < cfg.steps; ++s)
		{
			suspended_advect_cons_gpu(c, uz, vz, wz, ws, c2, scratch, g, dt); std::swap(c, c2);
			fpack_from_G_gpu(G, fpack, g, cpack);
			bed_boxfilter_gpu(fpack, nullptr, Ft, g, cpack);
			bed_interface_geom_gpu(Ft, nbx, nby, nbz, Ab, g);
			bed_column_geom_gpu(fpack, nbx, nby, nbz, beta, upx, upy, g, cpack);
			bedload_flux_gpu(taubx, tauby, beta, upx, upy, qx, qy, g, p);
			bedload_div_gpu(qx, qy, taubx, tauby, divq, g, /*periodic*/ 1); // upwind along the flow (τ) direction
			morpho_exner_gpu(G, c, taubx, tauby, beta, upx, upy, divq, nullptr, g, p, dt);
			avalanche_sweep_gpu(G, Gb, g, ap); std::swap(G, Gb);
			if (s % 200 == 0 || s == cfg.steps - 1)
			{
				double Vb = reduce_sum_gpu(G, ncol) * g.h * g.h, Vs = reduce_sum_gpu(c, np) * Vcell;
				double e = std::fabs(Vb + Vs - V0) / V0; if (e > mass_err) mass_err = e;
			}
		}
		double Vb = reduce_sum_gpu(G, ncol) * g.h * g.h, Vs = reduce_sum_gpu(c, np) * Vcell;
		res.mass_err = mass_err; res.vbed_frac_end = Vb / V0; res.vsusp_frac_end = Vs / V0;
		for (double* q : { taubx, tauby, G, Gb, c, c2, scratch, ws, uz, vz, wz, fpack, Ft, nbx, nby, nbz, Ab, beta, upx, upy, qx, qy, divq }) cudaFree(q);
		return res;
	}

	// ================= (d) threshold of motion ======================================
	ThresholdResult run_threshold(const ThresholdConfig& cfg)
	{
		ThresholdResult res;
		double Ds = grain_Dstar(cfg.d50, cfg.rho, cfg.rho_s, cfg.nu);
		double tcr_flat = theta_cr_sw(Ds); res.theta_cr_flat = tcr_flat;
		double gam = sed_submerged_gamma(cfg.rho, cfg.rho_s);
		double ws0 = settling_ws(cfg.d50, cfg.rho, cfg.rho_s, cfg.nu);

		MacGrid g; g.nx = 1; g.ny = 1; g.nz = 4; g.h = cfg.h;
		int ncol = 1;
		double cpack = SED_CPACK;
		double* G = dalloc0(ncol);
		double* taubx = dalloc0(ncol), *tauby = dalloc0(ncol);
		double* beta = dalloc0(ncol), *upx = dfill(ncol, 1.0), *upy = dalloc0(ncol);
		double* divz = dalloc0(ncol);

		MorphoParams p; p.d50 = cfg.d50; p.rho = cfg.rho; p.rho_s = cfg.rho_s; p.nu = cfg.nu; p.ws0 = ws0;
		p.cpack = cpack; p.morfac = 1.0; p.deposition_on = 0; p.erosion_on = 1; p.bedload_on = 1;
		double dt = 0.5 * g.h / std::max(ws0, 1e-6);

		double G0 = 2.0 * g.h * cpack;
		double onset = 0.0; bool no_motion_below = true;
		for (int s = 0; s < cfg.nsweep; ++s)
		{
			double frac = cfg.theta_lo_frac + (cfg.theta_hi_frac - cfg.theta_lo_frac) * s / (cfg.nsweep - 1);
			double theta = frac * tcr_flat;
			double tau = theta * gam * cfg.d50;
			cudaMemcpy(G, &G0, sizeof(double), cudaMemcpyHostToDevice);
			double txv = tau; cudaMemcpy(taubx, &txv, sizeof(double), cudaMemcpyHostToDevice);
			for (int it = 0; it < cfg.erode_steps; ++it)
				morpho_exner_gpu(G, nullptr, taubx, tauby, beta, upx, upy, divz, nullptr, g, p, dt);
			double Gend = 0; cudaMemcpy(&Gend, G, sizeof(double), cudaMemcpyDeviceToHost);
			double dG = std::fabs(Gend - G0);
			bool moved = dG > 1e-14 * G0;
			if (moved && onset == 0.0) onset = theta;
			if (theta <= tcr_flat && moved) no_motion_below = false;
		}
		res.theta_onset = onset;
		res.onset_err = tcr_flat > 0 ? std::fabs(onset - tcr_flat) / tcr_flat : 0.0;
		res.no_motion_below = no_motion_below;

		// bedload gate: q_b = 0 at 0.9·θ_cr, > 0 at 1.1·θ_cr (β=0 ⇒ slope factor = 1).
		auto qmag_at = [&](double theta) -> double {
			double tau = theta * gam * cfg.d50;
			double txv = tau, zero = 0.0;
			cudaMemcpy(taubx, &txv, sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(tauby, &zero, sizeof(double), cudaMemcpyHostToDevice);
			double* qx = dalloc0(ncol), *qy = dalloc0(ncol);
			bedload_flux_gpu(taubx, tauby, beta, upx, upy, qx, qy, g, p);
			double hqx = 0, hqy = 0; cudaMemcpy(&hqx, qx, sizeof(double), cudaMemcpyDeviceToHost); cudaMemcpy(&hqy, qy, sizeof(double), cudaMemcpyDeviceToHost);
			cudaFree(qx); cudaFree(qy);
			return std::sqrt(hqx * hqx + hqy * hqy);
		};
		double q_lo = qmag_at(0.9 * tcr_flat), q_hi = qmag_at(1.1 * tcr_flat);
		res.bedload_gate_ok = (q_lo == 0.0) && (q_hi > 0.0);

		for (double* q : { G, taubx, tauby, beta, upx, upy, divz }) cudaFree(q);
		return res;
	}
}
