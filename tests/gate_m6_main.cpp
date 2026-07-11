// gate_m6_main.cpp — M6 pile-scour calibration gate (Roulund et al. 2005 / RESEARCH V7).
//
// Live morphodynamic loop (fluid → τ_b wall model → suspended + bedload → Exner → avalanche)
// around a voxelized vertical pile in an erodible sand bed, at the Roulund lab benchmark:
//   D = 0.1 m, water depth h_dom = 0.4 m, V = 0.46 m/s, d50 = 0.26 mm, live-bed V/Vcr = 1.25,
//   voxel h = 1 cm (= D/10 — REEF3D: D/10 → <4% err; D/5 unstable, research/11).
// Domain 3.0 × 1.6 × (h_dom + sand_depth) m; pile centred 1.0 m (10D) from the inlet with
// ≥ 7D lateral clearance. (NB the PLAN/CLAUDE "≈19M cells" figure is an error: at h=0.01 the
// grid is 300×160×~60 ≈ 2.9M cells — flagged in the M6 report.)
//
// Protocol (RESEARCH §7 / research/14, /16):
//   Phase 1  frozen bed, spin up ≥ 5 flow-throughs, verify recovered ambient u* = 0.020 ± 10%
//            (research/14 gate) BEFORE enabling morphology.
//   Phase 2  ramp MORFAC 1→target over ~2 flow-throughs, then run ≥ 0.5·T of morphological time,
//            recording S(t) = upstream/max/downstream scour depth in the pile annulus.
//   Fit      BOTH S(t)=S_eq(1−e^(−t/T)) and Welzel a(1−1/(1+bt)); report S_eq/D and T.
// Gate (V7): upstream S/D = 1.25 ± 15%; timescale within factor 2 of the measured ~130 s; the
//            u* pre-check within 10%; flow stable; sand mass (bed+susp) conserved (CLOSED box).
//
// ⚠ SUPERVISED CALIBRATION (CLAUDE.md / HANDOVER): this harness RUNS the sim and REPORTS S/D vs
// the V7 band. The calibration knobs — van Rijn α (±30%), Smagorinsky Cs (0.10–0.12), the
// Winterwerp erosion_coeff (the de-facto bed-Exner erosion-rate knob in CLOSED mode; α drives the
// M4 suspended path, not this Exner), optional HSV shear multiplier — are config/CLI-settable and
// SWEPT, but the choice of which to FREEZE is the orchestrator's + MH's, NEVER this program's. It
// prints the numbers; it does not edit configs/calibrated.json.
//
// Links only libscour (Qt-free, GL-free). CLI: gate_m6 <config.json> [--tag T] [--csv path]
//   [--alpha A] [--cs C] [--erosion E] [--morfac M] [--h HH] [--nu-fluid NF] [--hsv X]
//   [--trun-frac F] [--spinup-ft S] [--quick].
#include "core/fluid/channel_core.h"
#include "core/fluid/channel_mask.h"
#include "core/fluid/mac_ops.h"
#include "core/sediment/seabed_engine.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace scour::core;

namespace
{
	double jd(const nlohmann::json& j, const char* k, double d) { return j.contains(k) ? j.at(k).get<double>() : d; }
	int ji(const nlohmann::json& j, const char* k, int d) { return j.contains(k) ? j.at(k).get<int>() : d; }
	bool jb(const nlohmann::json& j, const char* k, bool d) { return j.contains(k) ? j.at(k).get<bool>() : d; }
	nlohmann::json jsub(const nlohmann::json& j, const char* k) { return j.contains(k) ? j.at(k) : nlohmann::json::object(); }
	const char* PF(bool b) { return b ? "PASS" : "FAIL"; }

	// CLI double/flag override lookup (--key value / --flag).
	bool has_flag(int argc, char** argv, const char* f)
	{ for (int i = 1; i < argc; ++i) if (!std::strcmp(argv[i], f)) return true; return false; }
	double opt_d(int argc, char** argv, const char* f, double d)
	{ for (int i = 1; i < argc - 1; ++i) if (!std::strcmp(argv[i], f)) return std::atof(argv[i + 1]); return d; }
	std::string opt_s(int argc, char** argv, const char* f, const std::string& d)
	{ for (int i = 1; i < argc - 1; ++i) if (!std::strcmp(argv[i], f)) return argv[i + 1]; return d; }

	// One S(t) sample (all in morphological time / metres).
	struct Sample { double t = 0, S_up = 0, S_max = 0, S_down = 0, zb_min = 0, maxu = 0, ustar_amb = 0; };

	// ---- Nonlinear fits by 1-D grid search over the timescale, linear amplitude per candidate ----
	// S(t) = Seq·(1 − e^(−t/T)); for fixed T, Seq = Σ S_i b_i / Σ b_i², b_i = 1−e^(−t_i/T).
	struct FitExp { double Seq = 0, T = 0, sse = 1e300; };
	FitExp fit_exp(const std::vector<Sample>& s, double (*get)(const Sample&))
	{
		FitExp best;
		if (s.size() < 4) return best;
		double tmax = s.back().t;
		for (int it = 0; it < 400; ++it)
		{
			double T = 0.02 * tmax * std::pow(2000.0, it / 399.0); // 0.02·tmax … 40·tmax, log-spaced
			double num = 0, den = 0;
			for (const auto& p : s) { double b = 1.0 - std::exp(-p.t / T); num += get(p) * b; den += b * b; }
			if (den <= 0) continue;
			double Seq = num / den, sse = 0;
			for (const auto& p : s) { double f = Seq * (1.0 - std::exp(-p.t / T)); double e = get(p) - f; sse += e * e; }
			if (sse < best.sse) { best.sse = sse; best.Seq = Seq; best.T = T; }
		}
		return best;
	}
	// Welzel a·(1 − 1/(1+b·t)); for fixed b, a = Σ S_i g_i / Σ g_i², g_i = 1−1/(1+b t_i). Report S_eq=a.
	struct FitWel { double a = 0, b = 0, sse = 1e300; };
	FitWel fit_welzel(const std::vector<Sample>& s, double (*get)(const Sample&))
	{
		FitWel best;
		if (s.size() < 4) return best;
		double tmax = s.back().t;
		for (int it = 0; it < 400; ++it)
		{
			double b = (0.02 / tmax) * std::pow(2000.0, it / 399.0); // 1/(50·tmax) … 40/tmax
			double num = 0, den = 0;
			for (const auto& p : s) { double g = 1.0 - 1.0 / (1.0 + b * p.t); num += get(p) * g; den += g * g; }
			if (den <= 0) continue;
			double a = num / den, sse = 0;
			for (const auto& p : s) { double f = a * (1.0 - 1.0 / (1.0 + b * p.t)); double e = get(p) - f; sse += e * e; }
			if (sse < best.sse) { best.sse = sse; best.a = a; best.b = b; }
		}
		return best;
	}
	double g_up(const Sample& s) { return s.S_up; }
	double g_max(const Sample& s) { return s.S_max; }
}

int main(int argc, char** argv)
{
	// ---- config ----
	nlohmann::json j = nlohmann::json::object();
	std::string cfg_path;
	if (argc >= 2 && argv[1][0] != '-')
	{
		cfg_path = argv[1];
		try { std::ifstream in(cfg_path); std::ostringstream ss; ss << in.rdbuf(); j = nlohmann::json::parse(ss.str()); }
		catch (const std::exception& e) { std::fprintf(stderr, "gate_M6: bad config %s: %s\n", cfg_path.c_str(), e.what()); }
	}
	nlohmann::json jg = jsub(j, "gates");

	const bool quick = has_flag(argc, argv, "--quick");
	std::string tag = opt_s(argc, argv, "--tag", quick ? "quick" : "roulund");
	std::string csv = opt_s(argc, argv, "--csv", "");

	// Geometry (m). Domain z = water depth + sand reservoir (the scour hole eats into the reservoir).
	double Lx = jd(j, "domain_x", 3.0), Ly = jd(j, "domain_y", 1.6);
	double water = jd(j, "water_depth", 0.4), sand = jd(j, "sand_depth", 0.2);
	double h = opt_d(argc, argv, "--h", jd(j, "voxel_h", 0.01));
	double D = jd(j, "pile_D", 0.1), pxc = jd(j, "pile_x", 1.0), pyc = jd(j, "pile_y", 0.5 * Ly);
	if (quick) { h = std::max(h, 0.02); }             // coarser + smaller for a fast pipeline smoke
	double Lz = water + sand;

	// Physics.
	double U = jd(j, "U", 0.46), rho = jd(j, "rho", 1027.0), rho_s = jd(j, "rho_s", 2650.0);
	double d50 = jd(j, "d50", 0.26e-3), nu = jd(j, "nu", 1.36e-6);
	double Cs = opt_d(argc, argv, "--cs", jd(j, "Cs", 0.11));
	double alpha = opt_d(argc, argv, "--alpha", jd(j, "alpha", 0.00033));
	double erosion = opt_d(argc, argv, "--erosion", jd(j, "erosion_coeff", 0.018));
	double morfac_target = opt_d(argc, argv, "--morfac", jd(j, "morfac", 1.0));
	int bedload_formula = ji(j, "bedload_formula", 1); // Engelund–Fredsøe (Roulund form) per V7
	bool susp_diff = jb(j, "suspended_diffusion", true);
	double sigma_s = jd(j, "sigma_s", 0.7);
	double nu_fluid = opt_d(argc, argv, "--nu-fluid", jd(j, "nu_fluid", 0.0)); // 0 ⇒ molecular nu
	if (nu_fluid <= 0.0) nu_fluid = nu;
	int solid_mode = ji(j, "solid_mode", SOLID_FREESLIP); // RESEARCH §3: free-slip on voxel solids
	if (has_flag(argc, argv, "--noslip")) solid_mode = SOLID_NOSLIP; // sweep the pile/bed BL (HSV strength)
	if (has_flag(argc, argv, "--freeslip")) solid_mode = SOLID_FREESLIP;
	double hsv = opt_d(argc, argv, "--hsv", jd(j, "hsv_mult", 1.0)); // reserved (unused if 1.0)
	double cfl = jd(j, "cfl", 0.7), proj_tol = jd(j, "proj_tol", 1e-4);
	int proj_max_iter = ji(j, "proj_max_iter", 80);

	// Run protocol.
	double spinup_ft = opt_d(argc, argv, "--spinup-ft", jd(j, "spinup_flowthroughs", 5.0));
	double ramp_ft = jd(j, "ramp_flowthroughs", 2.0);
	double T_pred = jd(j, "T_pred", 130.0);
	double trun_frac = opt_d(argc, argv, "--trun-frac", jd(j, "trun_frac_T", 0.6));
	double sample_dt = jd(j, "sample_dt_s", 2.0);
	if (quick) { spinup_ft = 2.0; ramp_ft = 1.0; trun_frac = 0.03; sample_dt = 0.5; }

	// Gate tolerances.
	double us_target = jd(jg, "ustar_target", 0.020), us_tol = jd(jg, "ustar_tol", 0.10);
	double SD_target = jd(jg, "SD_target", 1.25), SD_tol = jd(jg, "SD_tol", 0.15);
	double T_factor = jd(jg, "T_factor", 2.0);
	double mass_tol = jd(jg, "mass_err_max", 0.05);

	// ---- grid ----
	MacGrid g;
	g.nx = std::max(8, (int)std::llround(Lx / h));
	g.ny = std::max(8, (int)std::llround(Ly / h));
	g.nz = std::max(8, (int)std::llround(Lz / h));
	g.h = h;
	long long ncell = (long long)g.nx * g.ny * g.nz;

	std::printf("========================================================================\n");
	std::printf("gate_M6 (Roulund V7) tag=%s  SUPERVISED CALIBRATION — reports S/D, does not freeze constants\n", tag.c_str());
	std::printf("gate_M6: grid %dx%dx%d h=%.4f m  (%.2fM cells)  domain %.2fx%.2fx%.2f m (water %.2f + sand %.2f)\n",
		g.nx, g.ny, g.nz, g.h, ncell / 1e6, Lx, Ly, Lz, water, sand);
	std::printf("gate_M6: pile D=%.3f m (%.0f cells) at (%.2f,%.2f); U=%.3f d50=%.3f mm | Cs=%.3f alpha=%.5f erosion=%.4f morfac=%.1f\n",
		D, D / h, pxc, pyc, U, d50 * 1e3, Cs, alpha, erosion, morfac_target);
	std::printf("gate_M6: nu_fluid=%.2e (%s) solid=%s bedload=%s susp_diff=%d hsv=%.2f | spinup=%.1fft ramp=%.1fft trun=%.2fT (T_pred=%.0fs)\n",
		nu_fluid, nu_fluid == nu ? "molecular" : "artificial", solid_mode == SOLID_FREESLIP ? "free-slip" : "no-slip",
		bedload_formula == 1 ? "E-F" : "W-P", (int)susp_diff, hsv, spinup_ft, ramp_ft, trun_frac, T_pred);
	std::fflush(stdout);

	// ---- pile mask (full-z cylinder = the rigid structure/foundation) ----
	std::vector<unsigned char> pile;
	int npile = build_cylinder_mask(g, pxc, pyc, 0.5 * D, pile);

	// ---- morpho engine + fluid core ----
	SeabedParams sp;
	sp.sand_depth = sand; sp.d50 = d50; sp.rho = rho; sp.rho_s = rho_s; sp.nu = nu;
	sp.Cs = Cs; sp.alpha = alpha; sp.erosion_coeff = erosion; sp.morfac = 1.0; // start at M=1 (ramp up)
	sp.sigma_s = sigma_s; sp.bedload_formula = bedload_formula; sp.diffusion_on = susp_diff ? 1 : 0;
	sp.sed_bc = SED_BC_CLOSED; // Roulund flume = closed box (mass-exact sand budget)
	SeabedMorpho engine(g, sp, pile);

	ChannelBC bc;
	bc.inlet_mode = INLET_LOGLAW;           // developed BL inflow (Roulund used a log-law inlet)
	bc.z0 = d50 / 12.0; bc.kappa = 0.40; bc.bed_datum = sand; // BL referenced to the flat bed top
	bc.ustar = loglaw_ustar_for_U(U, water, bc.z0, bc.kappa);
	bc.U_inlet = U; bc.Uc = U;
	bc.solid_mode = solid_mode;
	bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;

	ChannelParams pr;
	pr.rho = rho; pr.nu = nu_fluid; pr.Cs = Cs; pr.cfl = cfl; pr.safety = 0.9;
	pr.proj_tol = proj_tol; pr.proj_max_iter = proj_max_iter; pr.advect_band = 1;

	ChannelFluidCore core(g, bc, pr, engine.initial_flow_solid());
	core.set_bed_inlet_mask(true);
	core.init_inlet_profile(); // start from the log-law inlet profile everywhere (steady open channel)

	std::printf("gate_M6: pile cells=%d  inlet u*=%.5f m/s (flux-matched log-law, target ambient %.3f)\n",
		npile, bc.ustar, us_target);
	std::fflush(stdout);

	// Column geometry for the scour-depth annulus (host, one-time).
	const int ncol = g.nx * g.ny;
	std::vector<unsigned char> frozen(ncol, 0);
	std::vector<double> rdist(ncol, 0), colx(ncol, 0);
	for (int jj = 0; jj < g.ny; ++jj)
		for (int ii = 0; ii < g.nx; ++ii)
		{
			int c = jj * g.nx + ii;
			double x = (ii + 0.5) * g.h, y = (jj + 0.5) * g.h;
			colx[c] = x;
			rdist[c] = std::sqrt((x - pxc) * (x - pxc) + (y - pyc) * (y - pyc));
			// frozen footprint = any column with a pile (structure) cell
			for (int k = 0; k < g.nz; ++k) if (pile[g.pidx(ii, jj, k)]) { frozen[c] = 1; break; }
		}
	const double Rp = 0.5 * D;

	// M6 HSV shear multiplier (supervised knob): amplify τ_b in the near-pile scour ring [R,2.5R] to
	// restore the SL-under-resolved horseshoe-vortex bed shear (RESEARCH §11 #1; the erosion-rate knob
	// does NOT move equilibrium S/D — the flow shear does, see the M6 report). hsv=1 ⇒ no-op.
	if (hsv != 1.0)
	{
		std::vector<double> mult((size_t)ncol, 1.0);
		int nh = 0;
		for (int c = 0; c < ncol; ++c)
			if (!frozen[c] && rdist[c] > Rp && rdist[c] < 2.5 * Rp) { mult[c] = hsv; ++nh; }
		engine.set_shear_multiplier(mult);
		std::printf("gate_M6: HSV shear multiplier ×%.2f applied to %d near-pile columns [R,2.5R]\n", hsv, nh);
		std::fflush(stdout);
	}

	auto measure = [&](Sample& s)
	{
		std::vector<float> zb; engine.copy_zb_host(zb);
		std::vector<float> us; engine.copy_ustar_host(us);
		double zb0 = sand;
		double up = zb0, mx = zb0, dn = zb0, gmin = zb0; // track min z_b (⇒ max scour)
		for (int c = 0; c < ncol; ++c)
		{
			if (frozen[c]) continue;
			double z = zb[c];
			if (rdist[c] > Rp && rdist[c] < 3.0 * Rp) // pile annulus
			{
				mx = std::min(mx, z);
				if (colx[c] <= pxc) up = std::min(up, z); else dn = std::min(dn, z);
			}
			gmin = std::min(gmin, z);
		}
		s.S_up = zb0 - up; s.S_max = zb0 - mx; s.S_down = zb0 - dn; s.zb_min = gmin;
		// ambient u*: upstream band (x∈[0.03,0.15]·plausible), away from lateral walls & pile wake
		double sum = 0; int n = 0;
		int i_lo = std::max(2, (int)std::llround(0.03 / g.h)), i_hi = std::max(i_lo + 1, (int)std::llround(0.15 / g.h));
		int j_lo = g.ny / 4, j_hi = g.ny - g.ny / 4;
		for (int ii = i_lo; ii <= i_hi && ii < g.nx; ++ii)
			for (int jj = j_lo; jj < j_hi; ++jj)
			{
				int c = jj * g.nx + ii; if (frozen[c] || rdist[c] < 3.0 * Rp) continue;
				sum += us[c]; ++n;
			}
		s.ustar_amb = n ? sum / n : 0.0;
	};

	// ---- Phase 1: frozen-bed spin-up + u* pre-check ----
	engine.set_morphology_frozen(true);
	double spinup_time = spinup_ft * (Lx / U);
	double t_flow = 0.0; long long step = 0;
	double maxu_seen = 0.0; bool stable = true;
	std::vector<unsigned char> new_solid;
	std::printf("gate_M6: --- Phase 1: frozen-bed spin-up %.1f s flow (%.1f flow-throughs) ---\n", spinup_time, spinup_ft);
	std::fflush(stdout);
	while (t_flow < spinup_time)
	{
		double dt = core.step();
		engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid); // frozen: τ_b only
		t_flow += dt; ++step;
		if (step % 200 == 0 || t_flow >= spinup_time)
		{
			double mu = reduce_max_abs_gpu(core.u_dev(), g.u_count());
			maxu_seen = std::max(maxu_seen, mu);
			if (!std::isfinite(mu) || mu > 20.0 * U) stable = false;
			Sample s; measure(s);
			std::printf("gate_M6:   t_flow=%6.1fs step=%lld maxu=%.3f ambient u*=%.5f (iters=%d)\n",
				t_flow, step, mu, s.ustar_amb, core.last_solve_iters());
			std::fflush(stdout);
			if (!stable) break;
		}
	}
	Sample pre; measure(pre);
	double us_err = us_target > 0 ? std::fabs(pre.ustar_amb - us_target) / us_target : 1.0;
	bool pass_ustar = stable && (us_err <= us_tol);
	std::printf("gate_M6: [u* pre-check] recovered ambient u*=%.5f vs target %.5f  err=%.1f%% [< %.0f%%] %s\n",
		pre.ustar_amb, us_target, us_err * 100, us_tol * 100, PF(pass_ustar));
	std::fflush(stdout);

	// ---- Phase 2: MORFAC ramp + morphology, recording S(t) ----
	engine.set_morphology_frozen(false);
	double V0 = engine.bed_volume() + engine.susp_volume();
	double ramp_time = ramp_ft * (Lx / U);
	double t_morph = 0.0, t_morph_target = trun_frac * T_pred;
	double next_sample = 0.0, mass_err = 0.0;
	std::vector<Sample> series;
	std::printf("gate_M6: --- Phase 2: ramp M 1->%.1f over %.1fs, then run %.1fs morphological (%.2fT) ---\n",
		morfac_target, ramp_time, t_morph_target, trun_frac);
	std::fflush(stdout);

	double t_since_unfreeze = 0.0;
	long long morpho_updates = 0, remasks = 0;
	while (t_morph < t_morph_target && stable)
	{
		double dt = core.step();
		// MORFAC ramp (linear 1→target over ramp_time of flow), then hold.
		double M = morfac_target;
		if (t_since_unfreeze < ramp_time)
			M = 1.0 + (morfac_target - 1.0) * (t_since_unfreeze / ramp_time);
		engine.set_morfac(M);
		bool re = engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid);
		if (re) { core.update_solid(new_solid); ++remasks; }
		++morpho_updates;
		t_since_unfreeze += dt;
		t_morph += M * dt; // morphological time = Σ M·dt
		++step;

		if (t_morph >= next_sample)
		{
			Sample s; s.t = t_morph; measure(s);
			double mu = reduce_max_abs_gpu(core.u_dev(), g.u_count());
			s.maxu = mu; maxu_seen = std::max(maxu_seen, mu);
			double V = engine.bed_volume() + engine.susp_volume();
			mass_err = std::max(mass_err, V0 > 0 ? std::fabs(V - V0) / V0 : 0.0);
			if (!std::isfinite(mu) || mu > 20.0 * U) stable = false;
			series.push_back(s);
			std::printf("gate_M6:   t_morph=%7.1fs (M=%.2f) S_up/D=%.3f S_max/D=%.3f S_down/D=%.3f zb_min=%.4f maxu=%.3f mass_err=%.2e\n",
				s.t, M, s.S_up / D, s.S_max / D, s.S_down / D, s.zb_min, mu, mass_err);
			std::fflush(stdout);
			next_sample += sample_dt;
		}
	}

	// ---- fits ----
	FitExp fe_up = fit_exp(series, g_up), fe_mx = fit_exp(series, g_max);
	FitWel fw_up = fit_welzel(series, g_up);
	double SD_exp = fe_up.Seq / D, SD_wel = fw_up.a / D;
	double SD_last = series.empty() ? 0.0 : series.back().S_up / D;
	bool fit_ok = series.size() >= 4 && series.back().t > 0.3 * fe_up.T; // need ≥0.3T for a non-degenerate fit
	double fit_disagree = (fe_up.Seq > 0) ? std::fabs(fe_up.Seq - fw_up.a) / fe_up.Seq : 1.0;

	std::printf("gate_M6: --- fits (upstream S) ---\n");
	std::printf("gate_M6:   exp:    S_eq/D=%.3f  T=%.1f s\n", SD_exp, fe_up.T);
	std::printf("gate_M6:   Welzel: S_eq/D=%.3f  b=%.4f\n", SD_wel, fw_up.b);
	std::printf("gate_M6:   last-sample S_up/D=%.3f  S_max/D(exp)=%.3f  fits disagree %.0f%% (extend if >15%%)\n",
		SD_last, fe_mx.Seq / D, fit_disagree * 100);
	std::fflush(stdout);

	// ---- CSV ----
	if (!csv.empty())
	{
		std::ofstream out(csv);
		out << "# gate_M6 " << tag << " Cs=" << Cs << " alpha=" << alpha << " erosion=" << erosion
			<< " morfac=" << morfac_target << " h=" << h << " nu_fluid=" << nu_fluid << "\n";
		out << "# u*_ambient=" << pre.ustar_amb << " S_eq/D(exp)=" << SD_exp << " T=" << fe_up.T
			<< " S_eq/D(welzel)=" << SD_wel << "\n";
		out << "t_morph_s,S_up,S_max,S_down,zb_min,maxu,ustar_amb\n";
		for (const auto& s : series)
			out << s.t << "," << s.S_up << "," << s.S_max << "," << s.S_down << "," << s.zb_min
				<< "," << s.maxu << "," << s.ustar_amb << "\n";
		std::printf("gate_M6: wrote S(t) CSV -> %s\n", csv.c_str());
	}

	// ---- gate (V7) ----
	bool pass_SD = std::fabs(SD_exp - SD_target) <= SD_tol * SD_target && fit_ok;
	bool pass_T = fe_up.T >= T_pred / T_factor && fe_up.T <= T_pred * T_factor && fit_ok;
	bool pass_mass = mass_err < mass_tol;
	bool pass = pass_ustar && stable && pass_SD && pass_T && pass_mass;

	std::printf("========================================================================\n");
	std::printf("gate_M6: u*_precheck=%s  stable=%s(maxu=%.3f)  mass=%s(%.2e<%.0e)\n",
		PF(pass_ustar), PF(stable), maxu_seen, PF(pass_mass), mass_err, mass_tol);
	std::printf("gate_M6: S/D(exp)=%.3f [%.3f±%.0f%%]=%s  T=%.1fs [%.0f..%.0f]=%s  (fit_ok=%d, samples=%zu)\n",
		SD_exp, SD_target, SD_tol * 100, PF(pass_SD), fe_up.T, T_pred / T_factor, T_pred * T_factor, PF(pass_T),
		(int)fit_ok, series.size());
	std::printf("gate_M6: RESULT %s  %s\n", pass ? "PASS" : "FAIL",
		pass ? "" : "(uncalibrated / needs MH calibration — see M6 report; NOT a harness failure per se)");
	std::fflush(stdout);
	return pass ? 0 : 1;
}
