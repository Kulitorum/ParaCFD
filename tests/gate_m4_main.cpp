// gate_m4_main.cpp — M4 acceptance gate (V5 settling column + V6 Rouse profile). Consumes
// configs/m4_sediment.json. Prints "gate_M4" lines and returns 0 on PASS. RESEARCH §6.1–6.2,
// research/02 §7, research/11 §6–7.
//
// PASS requires (PLAN M4):
//   * settling column: L2(numeric − exact translate) < 2% AND sand-mass error < 0.1%
//     (still-water column is plug flow: the slab settles at exactly w_s and deposits via the
//      flux BC D=w_s·c_b — research/11 §6. RESEARCH §10's c0·e^(−w_s t/h) is the well-mixed
//      lumped idealisation and is reported as a diagnostic only.);
//   * Rouse: fitted Rouse exponent within 15% of R = w_s/(β·κ·u*)  (=0.91 at the benchmark ν).
#include "core/sediment/sed_validation.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

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
		catch (const std::exception& e) { std::fprintf(stderr, "gate_M4: bad config %s: %s\n", argv[1], e.what()); }
	}
	nlohmann::json js = jsub(j, "settling"), jr = jsub(j, "rouse"), jg = jsub(j, "gates");

	double l2_max = jd(jg, "l2_max", 0.02);
	double mass_err_max = jd(jg, "mass_err_max", 0.001);
	double vcen_tol = jd(jg, "v_centroid_tol", 0.02);
	double R_err_max = jd(jg, "R_err_max", 0.15);
	double r2_min = jd(jg, "r2_min", 0.98);

	std::printf("gate_M4: suspended sediment — settling column (V5) + Rouse profile (V6). RESEARCH §6.1-6.2.\n");
	std::fflush(stdout);

	// ---- V5: settling column ----
	SettlingColumnConfig sc;
	sc.Lz = jd(js, "Lz", 1.0); sc.h = jd(js, "h", 0.005); sc.nx = ji(js, "nx", 4); sc.ny = ji(js, "ny", 4);
	sc.d50 = jd(js, "d50", 0.2e-3); sc.rho = jd(js, "rho", 1027.0); sc.rho_s = jd(js, "rho_s", 2650.0);
	sc.nu = jd(js, "nu", 1.36e-6); sc.c0 = jd(js, "c0", 0.01); sc.z_center = jd(js, "z_center", 0.7);
	sc.sigma = jd(js, "sigma", 0.06); sc.cfl = jd(js, "cfl", 0.7); sc.travel_check = jd(js, "travel_check", 0.30);
	SettlingColumnResult sr = run_settling_column(sc);

	// diagnostic: the well-mixed lumped model c(t)=c0·e^(−w_s t/h) vs the (correct) plug-flow linear
	// depth-average 1−w_s t/h, at the checkpoint time — to make the RESEARCH §10 distinction explicit.
	double x = sr.ws * sr.t_check / sc.Lz;
	std::printf("gate_M4: [V5 settling] %dx%dx%d w_s=%.5f m/s dt=%.4f s  t_check=%.2f s\n",
		sr.nx, sr.ny, sr.nz, sr.ws, sr.dt, sr.t_check);
	std::printf("gate_M4:   L2(num vs exact translate)=%.4f%%  [< %.1f%%] %s\n",
		sr.l2_translate * 100, l2_max * 100, sr.l2_translate < l2_max ? "PASS" : "FAIL");
	std::printf("gate_M4:   centroid speed=%.5f m/s (err %.3f%% vs w_s) [< %.1f%%] %s\n",
		sr.v_centroid, sr.v_centroid_err * 100, vcen_tol * 100, sr.v_centroid_err < vcen_tol ? "PASS" : "FAIL");
	std::printf("gate_M4:   sand-mass error (suspended+deposited thru bed flux BC)=%.5f%%  [< %.2f%%] %s   susp_frac_end=%.3f\n",
		sr.mass_err * 100, mass_err_max * 100, sr.mass_err < mass_err_max ? "PASS" : "FAIL", sr.susp_frac_end);
	std::printf("gate_M4:   [diag] depth-avg plug-flow(1−w_s t/h)=%.3f vs well-mixed e^(−w_s t/h)=%.3f (RESEARCH §10 is the latter)\n",
		1.0 - x, std::exp(-x));
	std::fflush(stdout);

	bool pass_l2 = sr.l2_translate < l2_max;
	bool pass_vcen = sr.v_centroid_err < vcen_tol;
	bool pass_mass = sr.mass_err < mass_err_max; // total sand (suspended + deposited) conservation

	// ---- V6: Rouse profile ----
	RouseChannelConfig rc;
	rc.h_dom = jd(jr, "h_dom", 0.4); rc.u_star = jd(jr, "u_star", 0.02); rc.d50 = jd(jr, "d50", 0.1e-3);
	rc.rho = jd(jr, "rho", 1025.0); rc.rho_s = jd(jr, "rho_s", 2650.0); rc.nu = jd(jr, "nu", 1.05e-6);
	rc.kappa = jd(jr, "kappa", 0.40); rc.beta = jd(jr, "beta", 1.0); rc.a_frac = jd(jr, "a_frac", 0.05);
	rc.c_a = jd(jr, "c_a", 1.0e-3); rc.nz = ji(jr, "nz", 80); rc.nx = ji(jr, "nx", 4); rc.ny = ji(jr, "ny", 4);
	rc.tol = jd(jr, "tol", 1e-6); rc.steps_max = ji(jr, "steps_max", 400000);
	RouseChannelResult rr = run_rouse_channel(rc);

	std::printf("gate_M4: [V6 Rouse] nz=%d w_s=%.5f m/s dt=%.4f s steps=%d t=%.0fs conv=%d\n",
		rr.nz, rr.ws, rr.dt, rr.steps, rr.sim_time, (int)rr.converged);
	std::printf("gate_M4:   R_target=%.4f (=w_s/(β·κ·u*))  R_fit=%.4f (R2=%.4f)  err=%.2f%%  [< %.0f%%] %s\n",
		rr.R_target, rr.R_fit, rr.fit_r2, rr.R_err * 100, R_err_max * 100, rr.R_err < R_err_max ? "PASS" : "FAIL");
	std::fflush(stdout);

	bool pass_R = rr.R_err < R_err_max;
	bool pass_r2 = rr.fit_r2 > r2_min;

	bool pass = pass_l2 && pass_vcen && pass_mass && pass_R && pass_r2;
	std::printf("gate_M4: [%s] L2  [%s] centroid=w_s  [%s] mass<0.1%%  [%s] Rouse_R±15%%  [%s] fit_R2\n",
		pass_l2 ? "PASS" : "FAIL", pass_vcen ? "PASS" : "FAIL", pass_mass ? "PASS" : "FAIL",
		pass_R ? "PASS" : "FAIL", pass_r2 ? "PASS" : "FAIL");
	std::printf("gate_M4: RESULT %s\n", pass ? "PASS" : "FAIL");
	return pass ? 0 : 1;
}
