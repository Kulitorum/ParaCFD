// gate_m5_main.cpp — M5 acceptance gate (bed state + morphodynamics). Consumes configs/
// m5_morpho.json. Prints "gate_M5" lines and returns 0 on PASS. RESEARCH §5, §7; PLAN M5.
//
// PASS requires (PLAN M5):
//   (a-i)  still-water deposition: Σf_pack·h rises at w_s·c0 (z_b at w_s·c0/0.64) < 0.5%,
//          three-reservoir sand-mass error < 0.1%;
//   (a-ii) prescribed bedload q_b(x)=q0·sin(2πx/L), E=D=0: (1−p)∂z_b/∂t=−∂q_b/∂x (upwind ∇·q_b)
//          matched < 0.5% after 100 updates;
//   (b)    sand-pile relaxes to the 30–32° repose band everywhere (mass conserved);
//   (c)    three-reservoir mass conserved < 0.1% over 1e4 steps;
//   (d)    no motion below Soulsby–Whitehouse θ_cr; onset within 10%.
#include "core/sediment/morpho_validation.h"

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
	const char* PF(bool b) { return b ? "PASS" : "FAIL"; }
}

int main(int argc, char** argv)
{
	nlohmann::json j = nlohmann::json::object();
	if (argc >= 2)
	{
		try { std::ifstream in(argv[1]); std::ostringstream ss; ss << in.rdbuf(); j = nlohmann::json::parse(ss.str()); }
		catch (const std::exception& e) { std::fprintf(stderr, "gate_M5: bad config %s: %s\n", argv[1], e.what()); }
	}
	nlohmann::json ja = jsub(j, "deposition"), jb = jsub(j, "bedload_exner"), jp = jsub(j, "sandpile"),
		jc = jsub(j, "mass"), jt = jsub(j, "threshold"), jg = jsub(j, "gates");

	double rate_max = jd(jg, "rate_err_max", 0.005);
	double mass_max = jd(jg, "mass_err_max", 0.001);
	double l2_max = jd(jg, "l2_err_max", 0.005);
	double slope_lo = jd(jg, "slope_lo_deg", 29.5), slope_hi = jd(jg, "slope_hi_deg", 32.0);
	double onset_max = jd(jg, "onset_err_max", 0.10);

	std::printf("gate_M5: bed state + morphodynamics (RESEARCH §5, §7). f_pack + states + interface geom +\n");
	std::printf("gate_M5: Winterwerp erosion + Wong-Parker/E-F bedload + slope θ_cr + avalanche + Exner + MORFAC.\n");
	std::fflush(stdout);

	// ---- (a-i) still-water deposition ----
	DepositionColumnConfig da;
	da.Lz = jd(ja, "Lz", 1.0); da.h = jd(ja, "h", 0.01); da.nx = ji(ja, "nx", 4); da.ny = ji(ja, "ny", 4);
	da.d50 = jd(ja, "d50", 0.2e-3); da.rho = jd(ja, "rho", 1027.0); da.rho_s = jd(ja, "rho_s", 2650.0);
	da.nu = jd(ja, "nu", 1.36e-6); da.c0 = jd(ja, "c0", 0.004); da.cfl = jd(ja, "cfl", 0.7);
	DepositionColumnResult ra = run_deposition_column(da);
	bool pa_rate = ra.rate_err < rate_max, pa_mass = ra.mass_err < mass_max;
	std::printf("gate_M5: [a-i deposition] %dx%dx%d w_s=%.5f m/s dt=%.4f s\n", ra.nx, ra.ny, ra.nz, ra.ws, ra.dt);
	std::printf("gate_M5:   dG/dt vs w_s*c0: err=%.4f%% [< %.1f%%] %s   z_b rise=%.3e m/s (=w_s*c0/0.64)\n",
		ra.rate_err * 100, rate_max * 100, PF(pa_rate), ra.zb_rate);
	std::printf("gate_M5:   three-reservoir mass err=%.5f%% [< %.2f%%] %s\n", ra.mass_err * 100, mass_max * 100, PF(pa_mass));
	std::fflush(stdout);

	// ---- (a-ii) prescribed bedload ----
	BedloadExnerConfig db;
	db.L = jd(jb, "L", 10.0); db.nx = ji(jb, "nx", 1024); db.q0 = jd(jb, "q0", 1.0e-6);
	db.porosity = jd(jb, "porosity", 0.36); db.updates = ji(jb, "updates", 100); db.dt = jd(jb, "dt", 0.0);
	BedloadExnerResult rb = run_bedload_exner(db);
	bool pb = rb.l2_err < l2_max;
	std::printf("gate_M5: [a-ii bedload Exner] nx=%d updates=%d dt=%.4f s dz_max=%.3e m\n", rb.nx, rb.updates, rb.dt, rb.dz_max);
	std::printf("gate_M5:   upwind (1-p)dz_b/dt=-dq/dx L2 err=%.4f%% [< %.1f%%] %s\n", rb.l2_err * 100, l2_max * 100, PF(pb));
	std::fflush(stdout);

	// ---- (b) sand-pile ----
	SandpileConfig dp;
	dp.nx = ji(jp, "nx", 61); dp.ny = ji(jp, "ny", 61); dp.h = jd(jp, "h", 0.05); dp.peak = jd(jp, "peak", 2.0);
	dp.init_slope_deg = jd(jp, "init_slope_deg", 45.0); dp.phi_trigger_deg = jd(jp, "phi_trigger_deg", 32.0);
	dp.phi_repose_deg = jd(jp, "phi_repose_deg", 30.0); dp.relax = jd(jp, "relax", 0.1); dp.max_sweeps = ji(jp, "max_sweeps", 40000);
	SandpileResult rp = run_sandpile(dp);
	bool pp = rp.converged && rp.final_max_slope_deg >= slope_lo && rp.final_max_slope_deg <= slope_hi && rp.mass_err < mass_max;
	std::printf("gate_M5: [b sand-pile] sweeps=%d init_slope=%.1f deg -> final=%.2f deg [%.1f-%.1f] conv=%d mass_err=%.2e\n",
		rp.sweeps, rp.init_max_slope_deg, rp.final_max_slope_deg, slope_lo, slope_hi, (int)rp.converged, rp.mass_err);
	std::printf("gate_M5:   relaxes to repose band + mass-conserving: %s\n", PF(pp));
	std::fflush(stdout);

	// ---- (c) three-reservoir mass conservation ----
	MassConservationConfig dc;
	dc.nx = ji(jc, "nx", 16); dc.ny = ji(jc, "ny", 16); dc.nz = ji(jc, "nz", 20); dc.h = jd(jc, "h", 0.05);
	dc.d50 = jd(jc, "d50", 0.35e-3); dc.rho = jd(jc, "rho", 1027.0); dc.rho_s = jd(jc, "rho_s", 2650.0);
	dc.nu = jd(jc, "nu", 1.36e-6); dc.c0 = jd(jc, "c0", 0.002); dc.tau_peak = jd(jc, "tau_peak", 0.0);
	dc.steps = ji(jc, "steps", 10000); dc.cfl = jd(jc, "cfl", 0.5);
	MassConservationResult rc = run_mass_conservation(dc);
	bool pc = rc.mass_err < mass_max;
	std::printf("gate_M5: [c mass conservation] %d steps dt=%.4f s  V_bed_end=%.3f V_susp_end=%.3f (of V0)\n",
		rc.steps, rc.dt, rc.vbed_frac_end, rc.vsusp_frac_end);
	std::printf("gate_M5:   V_bed+V_susp conserved: err=%.5f%% [< %.2f%%] %s\n", rc.mass_err * 100, mass_max * 100, PF(pc));
	std::fflush(stdout);

	// ---- (d) threshold of motion ----
	ThresholdConfig dt_;
	dt_.d50 = jd(jt, "d50", 0.2e-3); dt_.rho = jd(jt, "rho", 1025.0); dt_.rho_s = jd(jt, "rho_s", 2650.0);
	dt_.nu = jd(jt, "nu", 1.05e-6); dt_.nsweep = ji(jt, "nsweep", 200); dt_.h = jd(jt, "h", 0.05);
	dt_.erode_steps = ji(jt, "erode_steps", 20);
	ThresholdResult rt = run_threshold(dt_);
	bool pd = rt.no_motion_below && rt.onset_err < onset_max && rt.bedload_gate_ok;
	std::printf("gate_M5: [d threshold] theta_cr(S-W)=%.4f  onset=%.4f  err=%.2f%% [< %.0f%%]\n",
		rt.theta_cr_flat, rt.theta_onset, rt.onset_err * 100, onset_max * 100);
	std::printf("gate_M5:   no motion below theta_cr=%s  onset within 10%%=%s  bedload gate=%s\n",
		PF(rt.no_motion_below), PF(rt.onset_err < onset_max), PF(rt.bedload_gate_ok));
	std::fflush(stdout);

	bool pass = pa_rate && pa_mass && pb && pp && pc && pd;
	std::printf("gate_M5: [%s] a-i.rate [%s] a-i.mass [%s] a-ii.exner [%s] b.pile [%s] c.mass [%s] d.threshold\n",
		PF(pa_rate), PF(pa_mass), PF(pb), PF(pp), PF(pc), PF(pd));
	std::printf("gate_M5: RESULT %s\n", pass ? "PASS" : "FAIL");
	return pass ? 0 : 1;
}
