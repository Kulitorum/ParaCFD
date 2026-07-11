// gate_m1_main.cpp — the M1 acceptance gate (V1): 128^3 lid-driven cavity vs the
// Ghia et al. (1982) tables (RESEARCH §10 V1, research/11). Consumes
// configs/v1_cavity_re100.json. Prints "gate_M1" lines and returns 0 on PASS.
//
// PASS requires (PLAN M1):
//   * Re=100 RMS centreline error of u(y) and v(x) < 5% (normalised by lid speed)
//   * MGPCG reaches ||r||/||b|| <= 1e-4 in <= 14 iterations, and max|div u| < 1e-4 U/h
//   * validation-tolerance projection gives max|div u| < 1e-6 U/h
//   * Re=1000 is logged as a diagnostic (never hard-fails).
#include "core/fluid/cavity.h"

#include <nlohmann/json.hpp>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using namespace scour::core;

namespace
{
	double jget(const nlohmann::json& j, const char* k, double d) { return j.contains(k) ? j.at(k).get<double>() : d; }
	int jgeti(const nlohmann::json& j, const char* k, int d) { return j.contains(k) ? j.at(k).get<int>() : d; }
}

int main(int argc, char** argv)
{
	// Config (optional path in argv[1]).
	int N = 128, max_steps = 60000, check_interval = 1000, diag_max_steps = 30000;
	double Re = 100.0, Re_diag = 1000.0, rms_tol = 0.05, steady_tol = 1e-5;
	int mgpcg_max_iters = 14;
	double div_tol_prod = 1e-4, div_tol_valid = 1e-6;
	if (argc >= 2)
	{
		try
		{
			std::ifstream in(argv[1]); std::ostringstream ss; ss << in.rdbuf();
			auto j = nlohmann::json::parse(ss.str());
			N = jgeti(j, "N", N); Re = jget(j, "Re", Re); Re_diag = jget(j, "Re_diag", Re_diag);
			rms_tol = jget(j, "rms_tol", rms_tol); max_steps = jgeti(j, "max_steps", max_steps);
			check_interval = jgeti(j, "check_interval", check_interval); steady_tol = jget(j, "steady_tol", steady_tol);
			diag_max_steps = jgeti(j, "diag_max_steps", diag_max_steps);
			mgpcg_max_iters = jgeti(j, "mgpcg_max_iters", mgpcg_max_iters);
			div_tol_prod = jget(j, "div_tol_prod", div_tol_prod); div_tol_valid = jget(j, "div_tol_valid", div_tol_valid);
		}
		catch (const std::exception& e) { std::fprintf(stderr, "gate_M1: bad config %s: %s\n", argv[1], e.what()); }
	}
	const double U_over_h = 1.0 * N; // U=1, h=1/N
	const double div_bound_prod = div_tol_prod * U_over_h;
	const double div_bound_valid = div_tol_valid * U_over_h;

	std::printf("gate_M1: lid-driven cavity %d^3, Re=%.0f (Cs=0, nu=U*L/Re)\n", N, Re);
	std::fflush(stdout);

	CavityMetrics m = run_cavity(N, Re, max_steps, check_interval, steady_tol, /*probe=*/true);

	std::printf("gate_M1: Re=%.0f steps=%d dt=%.3e final_change=%.2e\n", Re, m.steps, m.dt_last, m.final_change);
	std::printf("gate_M1: RMS u(y) error = %.4f (%.2f%%)  [< %.2f required]\n", m.rms_u, 100 * m.rms_u, rms_tol);
	std::printf("gate_M1: RMS v(x) error = %.4f (%.2f%%)  [< %.2f required]\n", m.rms_v, 100 * m.rms_v, rms_tol);
	std::printf("gate_M1: centreline u_min = %.5f (Ghia Re=100: -0.21090)\n", m.umin_center);
	std::printf("gate_M1: run MGPCG max iters (warm) = %d, max relres = %.2e\n", m.run_max_iters, m.run_max_relres);

	const ProjectionProbe& pp = m.probe_prod;
	const ProjectionProbe& pv = m.probe_valid;
	std::printf("gate_M1: MGPCG production tol=1e-4 -> iters=%d relres=%.2e ; max|div| before=%.3e after=%.3e (bound %.3e)\n",
		pp.solve.iters, pp.solve.relres, pp.maxdiv_before, pp.maxdiv_after, div_bound_prod);
	std::printf("gate_M1: MGPCG validation      -> iters=%d relres=%.2e ; max|div| after=%.3e (bound %.3e)\n",
		pv.solve.iters, pv.solve.relres, pv.maxdiv_after, div_bound_valid);
	std::printf("gate_M1: Jacobi fallback sweeps=%d relres=%.2e (debug path OK)\n", m.probe_jacobi.solve.iters, m.probe_jacobi.solve.relres);

	bool pass_rms = (m.rms_u < rms_tol) && (m.rms_v < rms_tol);
	bool pass_iters = (pp.solve.iters <= mgpcg_max_iters) && pp.solve.converged && (m.run_max_iters <= mgpcg_max_iters);
	bool pass_div_prod = pp.maxdiv_after < div_bound_prod;
	bool pass_div_valid = pv.maxdiv_after < div_bound_valid;

	std::printf("gate_M1: [%s] RMS   [%s] MGPCG<=%d iters   [%s] div<%gU/h prod   [%s] div<%gU/h valid\n",
		pass_rms ? "PASS" : "FAIL", pass_iters ? "PASS" : "FAIL", mgpcg_max_iters,
		pass_div_prod ? "PASS" : "FAIL", div_tol_prod, pass_div_valid ? "PASS" : "FAIL", div_tol_valid);

	// Re=1000 diagnostic (log only, never hard-fail — RESEARCH §11).
	CavityMetrics d = run_cavity(N, Re_diag, diag_max_steps, check_interval, steady_tol, /*probe=*/false);
	std::printf("gate_M1: [DIAG] Re=%.0f steps=%d RMS u=%.2f%% v=%.2f%% u_min=%.5f (Ghia: -0.38289) — diagnostic, not a hard-fail\n",
		Re_diag, d.steps, 100 * d.rms_u, 100 * d.rms_v, d.umin_center);

	bool pass = pass_rms && pass_iters && pass_div_prod && pass_div_valid;
	std::printf("gate_M1: RESULT %s\n", pass ? "PASS" : "FAIL");
	return pass ? 0 : 1;
}
