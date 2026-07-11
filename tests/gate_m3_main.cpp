// gate_m3_main.cpp — M3 acceptance gate (V2): bed-shear log-law wall function + SEM inlet.
// Consumes configs/v2_channel_loglaw.json. Prints "gate_M3" lines and returns 0 on PASS.
//
// PASS requires (PLAN M3 / RESEARCH §4 anchor, §8):
//   * flat-bed periodic channel, h_dom=5 m, U=1 m/s, per d50 case: extracted u* within 10%
//     of κ·U/(ln(h_dom/z0)−1), z0=d50/12  (u*≈0.0345 @0.2 mm, 0.040 @1.0 mm);
//   * recovered κ within 5% from the log-law profile fit;
//   * τ_b smooth over the bed plane (no grid-pitch banding);
//   * SEM run holds 5–10% ambient TI at mid-domain with no divergence growth.
#include "core/fluid/channel_periodic.h"
#include "core/fluid/channel_sem.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace scour::core;

namespace
{
	double jd(const nlohmann::json& j, const char* k, double d) { return j.contains(k) ? j.at(k).get<double>() : d; }
	int ji(const nlohmann::json& j, const char* k, int d) { return j.contains(k) ? j.at(k).get<int>() : d; }
	nlohmann::json jsub(const nlohmann::json& j, const char* k) { return j.contains(k) ? j.at(k) : nlohmann::json::object(); }
}

int main(int argc, char** argv)
{
	nlohmann::json j = nlohmann::json::object();
	if (argc >= 2)
	{
		try { std::ifstream in(argv[1]); std::ostringstream ss; ss << in.rdbuf(); j = nlohmann::json::parse(ss.str()); }
		catch (const std::exception& e) { std::fprintf(stderr, "gate_M3: bad config %s: %s\n", argv[1], e.what()); }
	}

	nlohmann::json jp = jsub(j, "periodic");
	nlohmann::json jg = jsub(j, "gates");
	double ustar_tol = jd(jg, "ustar_tol", 0.10);
	double kappa_tol = jd(jg, "kappa_tol", 0.05);
	double tau_band_max = jd(jg, "tau_band_max", 0.05);
	double ti_lo = jd(jg, "ti_lo", 0.05), ti_hi = jd(jg, "ti_hi", 0.10);
	double div_growth = jd(jg, "div_growth", 1.5), div_max_norm = jd(jg, "div_max_norm", 1e-2);

	std::vector<double> d50s;
	if (j.contains("cases")) for (auto& c : j.at("cases")) d50s.push_back(jd(c, "d50", 0.2e-3));
	if (d50s.empty()) { d50s.push_back(0.2e-3); d50s.push_back(1.0e-3); }

	std::printf("gate_M3: bed-shear log-law wall function + Jarrin SEM (V2). RESEARCH §4/§8.\n");
	std::fflush(stdout);

	bool pass_ustar = true, pass_kappa = true, pass_band = true;
	for (double d50 : d50s)
	{
		PeriodicChannelConfig pc;
		pc.h = jd(jp, "h", 0.05); pc.Lz = jd(jp, "Lz", 5.0);
		pc.nx = ji(jp, "nx", 8); pc.ny = ji(jp, "ny", 8);
		pc.U_d = jd(jp, "U_d", 1.0); pc.nu = jd(jp, "nu", 1.36e-6);
		pc.rho = jd(jp, "rho", 1027.0); pc.kappa = jd(jp, "kappa", 0.40);
		pc.d50 = d50; pc.wall_regime = ji(jp, "wall_regime", 0);
		pc.steps_max = ji(jp, "steps_max", 400000);
		PeriodicChannelResult r = run_periodic_channel(pc);

		double eu_wall = std::fabs(r.ustar_wall - r.ustar_analytic) / r.ustar_analytic;
		double eu_mom = std::fabs(r.ustar_momentum - r.ustar_analytic) / r.ustar_analytic;
		double ek = std::fabs(r.kappa_fit - pc.kappa) / pc.kappa;
		bool pu = (eu_wall < ustar_tol), pk = (ek < kappa_tol), pb = (r.tau_band < tau_band_max);
		pass_ustar = pass_ustar && pu; pass_kappa = pass_kappa && pk; pass_band = pass_band && pb;

		std::printf("gate_M3: [periodic d50=%.3fmm] nz=%d z0=%.3e ks+=%.1f conv=%d steps=%d t=%.0fs\n",
			d50 * 1e3, r.nz, r.z0, r.ks_plus, (int)r.converged, r.steps, r.sim_time);
		std::printf("gate_M3:   u*_analytic=%.5f  u*_wall=%.5f (err %.1f%%)  u*_mom=%.5f (err %.1f%%)  [tol %.0f%%] %s\n",
			r.ustar_analytic, r.ustar_wall, eu_wall * 100, r.ustar_momentum, eu_mom * 100, ustar_tol * 100, pu ? "PASS" : "FAIL");
		std::printf("gate_M3:   kappa_fit=%.4f (err %.1f%%, R2=%.4f)  [tol %.0f%%] %s   U_bulk=%.5f gx=%.3e\n",
			r.kappa_fit, ek * 100, r.fit_r2, kappa_tol * 100, pk ? "PASS" : "FAIL", r.U_bulk, r.gx);
		std::printf("gate_M3:   tau_band=%.2e (flat-bed uniformity)  [< %.2f] %s   max|div|=%.2e\n",
			r.tau_band, tau_band_max, pb ? "PASS" : "FAIL", r.max_div);
		std::fflush(stdout);
	}

	// ---- SEM open channel: TI + divergence ----
	nlohmann::json js = jsub(j, "sem");
	ChannelSemConfig sc;
	sc.Lx = jd(js, "Lx", 10.0); sc.Ly = jd(js, "Ly", 10.0); sc.Lz = jd(js, "Lz", 5.0);
	sc.h = jd(js, "h", 0.1); sc.U = jd(js, "U", 1.0); sc.d50 = jd(js, "d50", 0.2e-3);
	sc.nu = jd(js, "nu", 1.36e-6); sc.Cs = jd(js, "Cs", 0.11); sc.gamma = jd(js, "gamma", 1.0);
	sc.N = ji(js, "N", 150); sc.sigma = jd(js, "sigma", 1.0);
	sc.spinup_flowthroughs = jd(js, "spinup_flowthroughs", 3.0);
	sc.record_flowthroughs = jd(js, "record_flowthroughs", 3.0);
	sc.sample_every = ji(js, "sample_every", 8);
	sc.z_ref_frac = jd(js, "z_ref_frac", 0.1);
	ChannelSemResult sr = run_channel_sem(sc);

	double div_norm = sr.div_max / (sc.U / sc.h);
	bool pass_ti = (sr.TI_ref >= ti_lo && sr.TI_ref <= ti_hi);
	bool pass_div = (div_norm < div_max_norm) && (sr.div_first < 1e-30 || sr.div_second <= div_growth * sr.div_first + 1e-12);

	std::printf("gate_M3: [SEM] %dx%dx%d h=%.3f u*=%.4f gamma=%.2f steps=%d samples=%d\n",
		sr.nx, sr.ny, sr.nz, sr.h, sr.ustar, sc.gamma, sr.steps, sr.samples);
	std::printf("gate_M3:   TI@z_ref=%.4f  TI_lowerband=%.4f  [in %.2f-%.2f] %s\n",
		sr.TI_ref, sr.TI_band_mean, ti_lo, ti_hi, pass_ti ? "PASS" : "FAIL");
	std::printf("gate_M3:   max|div|=%.2e (norm %.2e < %.1e)  half-means %.2e->%.2e (<%.1fx)  %s\n",
		sr.div_max, div_norm, div_max_norm, sr.div_first, sr.div_second, div_growth, pass_div ? "PASS" : "FAIL");
	std::fflush(stdout);

	bool pass = pass_ustar && pass_kappa && pass_band && pass_ti && pass_div;
	std::printf("gate_M3: [%s] u*   [%s] kappa   [%s] tau_smooth   [%s] SEM_TI   [%s] no_div_growth\n",
		pass_ustar ? "PASS" : "FAIL", pass_kappa ? "PASS" : "FAIL", pass_band ? "PASS" : "FAIL",
		pass_ti ? "PASS" : "FAIL", pass_div ? "PASS" : "FAIL");
	std::printf("gate_M3: RESULT %s\n", pass ? "PASS" : "FAIL");
	return pass ? 0 : 1;
}
