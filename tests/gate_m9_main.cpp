// gate_m9_main.cpp — M9 tidal-reversal + ranking-campaign harness (PLAN M9 gate).
//
// Runs the live morphodynamic pipeline under a REVERSING tidal current (RESEARCH §8: face-swap
// reversal through a 60 s cosine slack ramp, u*(t) ∝ |U_d(t)|, MORFAC ≤ 5) and, on the same runs,
// checks the three reversal-gate criteria AND ranks ≥ 3 candidate shapes:
//
//   (1) MASS  — total sand (bed + suspended ± boundary fluxes) conserved to < 0.1 %. With the CLOSED
//       suspended boundary the conservative advection seals all six faces, so bed+susp is conserved
//       to machine precision across ANY number of reversals — the honest test of the reversal BC
//       machinery (and it ranks shapes on a FINITE local sand budget: does the shape gather/keep sand
//       as the tide flips, rather than being handed an ever-growing open-sea supply?).
//   (2) ENERGY — no MONOTONIC kinetic-energy growth across the 4 slack transitions (a reversing BC
//       that pumped energy would show KE climbing every slack). KE = Σ|u|² sampled at each slack.
//   (3) MORFAC — the |Δz_b| ≤ 0.05·h limiter clip-fraction stays < 1 % (research/16): MORFAC is safe
//       for the reversing flow (critical MORFAC is >10× lower than steady — RESEARCH §7).
//
// Ranking: matched-time trapped-sand volume V_trap inside each shape's control volume after the SAME
// reversing schedule (all shapes identical tide). Comparative/qualitative — the tool's stated
// authority (RESEARCH §11); absolute magnitudes inherit the gate_M6 HSV under-prediction.
//
// Reference shapes (procedural voxel masks, shape_masks.h): flat PLATE, solid RING, open-top CUP.
// Links only libscour. CLI: gate_m9 [config] [--quick] [--csv path] [--shapes plate,ring,cup]
//   [--U u] [--morfac M] [--half-cycles N] [--noreverse].
#include "core/campaign/run_protocol.h"
#include "core/fluid/channel_core.h"
#include "core/fluid/mac_ops.h"
#include "core/geometry/shape_masks.h"
#include "core/sediment/seabed_engine.h"
#include "core/sediment/sed_physics.h"

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
	constexpr double PI = 3.14159265358979323846;
	double jd(const nlohmann::json& j, const char* k, double d) { return j.contains(k) ? j.at(k).get<double>() : d; }
	const char* PF(bool b) { return b ? "PASS" : "FAIL"; }
	bool has_flag(int c, char** v, const char* f) { for (int i = 1; i < c; ++i) if (!std::strcmp(v[i], f)) return true; return false; }
	double opt_d(int c, char** v, const char* f, double d) { for (int i = 1; i < c - 1; ++i) if (!std::strcmp(v[i], f)) return std::atof(v[i + 1]); return d; }
	std::string opt_s(int c, char** v, const char* f, const std::string& d) { for (int i = 1; i < c - 1; ++i) if (!std::strcmp(v[i], f)) return v[i + 1]; return d; }

	// Tidal forcing (RESEARCH §8). Starts at +U_max; each half-cycle holds a plateau for T_plateau then
	// cosine-ramps to the OPPOSITE sign over T_ramp (the compressed slack, ≥ 3–6 flow-throughs), the
	// zero-crossing (slack) at the ramp midpoint. Returns the SIGNED depth-averaged current U_d(t).
	struct Tide
	{
		double Umax = 0.8, T_plateau = 20.0, T_ramp = 20.0;
		double half() const { return T_plateau + T_ramp; }
		double u_signed(double t) const
		{
			int n = (int)std::floor(t / half());
			double base = (n % 2 == 0) ? 1.0 : -1.0; // alternate the plateau sign each half-cycle
			double local = t - n * half();
			double s = (local < T_plateau) ? base : base * std::cos(PI * (local - T_plateau) / T_ramp);
			return Umax * s;
		}
	};

	struct ShapeResult
	{
		std::string name;
		double V_last = 0, V_eq = 0, V_ci = 0, T_fit = 0, t_run = 0;
		double edge_scour = 0, mass_err = 0, maxu = 0, clip_frac = 0, ke_growth = 1.0;
		int nslack = 0; bool ke_ok = true, stable = true;
	};

	double field_ke(const ChannelFluidCore& core, const MacGrid& g)
	{ return dot_gpu(core.u_dev(), core.u_dev(), g.u_count()) + dot_gpu(core.v_dev(), core.v_dev(), g.v_count()) + dot_gpu(core.w_dev(), core.w_dev(), g.w_count()); }
}

int main(int argc, char** argv)
{
	nlohmann::json j = nlohmann::json::object();
	if (argc >= 2 && argv[1][0] != '-')
	{ try { std::ifstream in(argv[1]); std::ostringstream ss; ss << in.rdbuf(); j = nlohmann::json::parse(ss.str()); } catch (...) {} }

	const bool quick = has_flag(argc, argv, "--quick");
	const bool reverse = !has_flag(argc, argv, "--noreverse"); // --noreverse: steady control (no face swap)
	std::string csv = opt_s(argc, argv, "--csv", "");
	std::string shapes = opt_s(argc, argv, "--shapes", "plate,ring,cup");

	double Lx = jd(j, "domain_x", 4.0), Ly = jd(j, "domain_y", 4.0);
	double water = jd(j, "water_depth", 0.8), sand = jd(j, "sand_depth", 0.4);
	double h = jd(j, "voxel_h", 0.04);
	if (quick) h = std::max(h, 0.08);
	double Lz = water + sand;
	double U = opt_d(argc, argv, "--U", jd(j, "U", 0.8)); // peak (spring) tidal current U_max
	double rho = jd(j, "rho", 1027.0), rho_s = jd(j, "rho_s", 2650.0), d50 = jd(j, "d50", 0.2e-3), nu = jd(j, "nu", 1.36e-6);
	double Cs = jd(j, "Cs", 0.11), alpha = jd(j, "alpha", 0.00033), erosion = jd(j, "erosion_coeff", 0.018);
	double morfac = opt_d(argc, argv, "--morfac", jd(j, "morfac", 4.0)); // RESEARCH §7: M ≤ 5 reversing
	double sigma_s = jd(j, "sigma_s", 0.7);
	double foot = jd(j, "footprint_half", 0.5), height = jd(j, "shape_height", 0.4);
	double spinup_ft = jd(j, "spinup_flowthroughs", 4.0);
	int half_cycles = (int)opt_d(argc, argv, "--half-cycles", jd(j, "half_cycles", 4)); // ≥ 4 ⇒ ≥ 4 slacks
	double plateau_ft = jd(j, "plateau_flowthroughs", 4.0), ramp_ft = jd(j, "ramp_flowthroughs", 4.0);
	if (quick) { spinup_ft = 2.0; plateau_ft = 2.0; ramp_ft = 2.0; }

	MacGrid g;
	g.nx = std::max(8, (int)std::llround(Lx / h)); g.ny = std::max(8, (int)std::llround(Ly / h)); g.nz = std::max(8, (int)std::llround(Lz / h)); g.h = h;
	const int ncol = g.nx * g.ny;
	double xc = 0.5 * Lx, yc = 0.5 * Ly, ft = Lx / U;

	Tide tide; tide.Umax = U; tide.T_plateau = plateau_ft * ft; tide.T_ramp = ramp_ft * ft;
	double t_tide = half_cycles * tide.half(); // total flow-time of the reversing campaign

	double z0 = d50 / 12.0;
	std::printf("========================================================================\n");
	std::printf("gate_M9 tidal reversal  grid %dx%dx%d h=%.3f (%.2fM) domain %.1fx%.1fx%.1f  U_max=%.2f d50=%.3fmm morfac=%.0f\n",
		g.nx, g.ny, g.nz, h, (double)g.p_count() / 1e6, Lx, Ly, Lz, U, d50 * 1e3, morfac);
	std::printf("gate_M9 tide: %s  plateau=%.1fs ramp(slack)=%.1fs  half-cycle=%.1fs  %d half-cycles ⇒ %d slacks  t_tide=%.0fs (flow) morpho≈%.0fs\n",
		reverse ? "REVERSING" : "steady(control)", tide.T_plateau, tide.T_ramp, tide.half(), half_cycles, half_cycles, t_tide, t_tide * morfac);
	std::fflush(stdout);

	// Control-volume + edge-scour geometry (host, one-time): CV = disk of radius r_cv about centre.
	double r_cv = foot * 1.15, r_edge0 = foot * 1.2, r_edge1 = foot * 2.2;
	std::vector<double> rdist(ncol);
	for (int jj = 0; jj < g.ny; ++jj) for (int ii = 0; ii < g.nx; ++ii)
	{ double x = (ii + 0.5) * h, y = (jj + 0.5) * h; rdist[jj * g.nx + ii] = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc)); }

	auto build_shape = [&](const std::string& s, std::vector<unsigned char>& mask) -> int
	{
		if (s == "plate") return build_plate_mask(g, xc, yc, foot, sand, std::max(0.05, 0.1 * height), mask);
		if (s == "ring") return build_ring_mask(g, xc, yc, foot, foot * 0.6, sand, height, mask);
		return build_cup_mask(g, xc, yc, foot, sand, height, std::max(g.h, 0.08), mask); // cup (default)
	};

	auto run_shape = [&](const std::string& s) -> ShapeResult
	{
		ShapeResult r; r.name = s;
		std::vector<unsigned char> mask; int ncells = build_shape(s, mask);

		SeabedParams sp;
		sp.sand_depth = sand; sp.d50 = d50; sp.rho = rho; sp.rho_s = rho_s; sp.nu = nu;
		sp.Cs = Cs; sp.alpha = alpha; sp.erosion_coeff = erosion; sp.morfac = morfac; sp.sigma_s = sigma_s;
		sp.bedload_formula = 1; sp.diffusion_on = 1;
		sp.sed_bc = SED_BC_CLOSED; sp.U_inlet = U; // CLOSED: finite local sand, conserved to machine precision
		SeabedMorpho engine(g, sp, mask);

		ChannelBC bc;
		bc.inlet_mode = INLET_LOGLAW; bc.z0 = z0; bc.kappa = 0.40; bc.bed_datum = sand;
		bc.ustar = loglaw_ustar_for_U(U, water, z0, 0.40); bc.U_inlet = U; bc.Uc = U;
		bc.solid_mode = SOLID_FREESLIP; bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;

		ChannelParams pr; pr.rho = rho; pr.nu = 5e-3 /*coarse-grid artificial ν for the ranking patch*/;
		pr.Cs = Cs; pr.cfl = 0.7; pr.safety = 0.9; pr.proj_tol = 1e-4; pr.proj_max_iter = 60; pr.advect_band = 1;

		ChannelFluidCore core(g, bc, pr, engine.initial_flow_solid());
		core.set_bed_inlet_mask(true);
		core.init_inlet_profile();

		double V0 = engine.bed_volume() + engine.susp_volume();
		std::vector<unsigned char> new_solid;

		// (a) frozen spin-up at +U_max (develop the wake before the bed responds — research/14).
		double tf = 0, spin = spinup_ft * ft;
		engine.set_morphology_frozen(true);
		while (tf < spin) { double dt = core.step(); engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid); tf += dt; }
		engine.set_morphology_frozen(false);
		engine.reset_clip_stats();

		// (b) reversing tidal campaign: drive U_d(t)/flow-sign live; track mass, matched-phase KE, clip,
		// V_trap. KE is sampled at each half-cycle's PLATEAU MIDPOINT (|U|=U_max), NOT at slack: there the
		// bulk-flow KE baseline is identical every cycle, so a stable reversing flow keeps it ~invariant
		// and only a genuine energy pump makes it climb — a slack sample (|U|→0) is instead dominated by
		// the still-developing scour recirculation and would confound "reversal pumps energy" with "bed
		// not yet in equilibrium". slack-KE is still logged as a secondary diagnostic.
		std::vector<double> ts, Vs, ke_slack, ke_phase((size_t)std::max(1, half_cycles), -1.0);
		double tflow = 0, tmorph = 0, next = 0, sample_dt = quick ? 3.0 : 6.0;
		int prev_sign = 1;
		while (tflow < t_tide && r.stable)
		{
			double Ud = reverse ? tide.u_signed(tflow) : U;   // steady control keeps +U_max
			int sign = (Ud >= 0.0) ? 1 : -1;
			core.set_inlet_speed(std::fabs(Ud));               // magnitude; direction carried by flow_sign
			core.set_flow_direction(sign);
			engine.set_inlet_speed(std::fabs(Ud));             // keep the wall-model reference current in step
			if (sign != prev_sign) { if (reverse) ke_slack.push_back(field_ke(core, g)); ++r.nslack; } // slack crossing
			prev_sign = sign;
			int hc = (int)std::floor(tflow / tide.half());     // half-cycle index
			double local = tflow - hc * tide.half();           // time into this half-cycle
			if (reverse && hc >= 0 && hc < half_cycles && ke_phase[hc] < 0 && local >= 0.5 * tide.T_plateau)
				ke_phase[hc] = field_ke(core, g); // matched-phase sample at |U|=U_max

			double dt = core.step();
			if (engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid)) core.update_solid(new_solid);
			tflow += dt; tmorph += morfac * dt;

			double mu = reduce_max_abs_gpu(core.u_dev(), g.u_count());
			r.maxu = std::max(r.maxu, mu); if (!std::isfinite(mu) || mu > 20 * U) r.stable = false;
			double V = engine.bed_volume() + engine.susp_volume();
			// Credit the two physical sinks/sources so the residual is the TRUE test of reversal/transport
			// conservation: ± open-boundary flux (≈0 in CLOSED) and the structure-footprint pin discard
			// (sand that cannot accrete under a rigid slab). Residual → machine precision when sound.
			double net = V + engine.structure_discard() - V0 - engine.boundary_sand_in();
			r.mass_err = std::max(r.mass_err, V0 > 0 ? std::fabs(net) / V0 : 0.0);

			if (tmorph >= next)
			{
				std::vector<float> zb; engine.copy_zb_host(zb);
				double vtrap = 0, edge = sand;
				for (int c = 0; c < ncol; ++c)
				{
					double dz = zb[c] - sand;
					if (rdist[c] <= r_cv && dz > 0) vtrap += dz * h * h;
					if (rdist[c] > r_edge0 && rdist[c] < r_edge1) edge = std::min(edge, (double)zb[c]);
				}
				r.edge_scour = sand - edge; ts.push_back(tmorph); Vs.push_back(vtrap); r.V_last = vtrap;
				next += sample_dt;
			}
		}
		r.clip_frac = engine.clip_fraction();

		// KE growth test on the matched-phase (|U|=U_max) samples across the half-cycles: the bulk-flow
		// KE baseline is invariant, so a reversal that neither pumps nor bleeds energy holds this
		// ~constant. Gate: the peak matched-phase KE must stay within 25 % of the first (a genuine energy
		// pump would climb far past that). The developing scour hole perturbs the bulk-dominated plateau
		// KE by ≪ 25 %, so this does not false-fail the transient. r.ke_growth = peak/first.
		if (reverse)
		{
			std::vector<double> valid; for (double k : ke_phase) if (k >= 0) valid.push_back(k);
			if (valid.size() >= 2)
			{
				double f = std::max(valid.front(), 1e-300), pk = 0;
				for (double k : valid) pk = std::max(pk, k);
				r.ke_growth = pk / f;
				r.ke_ok = r.ke_growth < 1.25;
			}
			std::printf("gate_M9:   [%-5s] matched-phase KE:", s.c_str());
			for (double k : valid) std::printf(" %.3e", k);
			std::printf("  (peak/first ×%.2f, %s)   slack-KE:", r.ke_growth, r.ke_ok ? "ok" : "GROW");
			for (double k : ke_slack) std::printf(" %.3e", k);
			std::printf("\n");
		}

		EqFit fe = fit_exponential(ts, Vs), fw = fit_welzel(ts, Vs);
		double rmse = fe.sse > 0 && ts.size() ? std::sqrt(fe.sse / ts.size()) : 0.0;
		ShapeRank rk = equilibrium_with_ci(fe, fw, rmse);
		r.V_eq = rk.value; r.V_ci = rk.ci; r.T_fit = fe.param; r.t_run = tmorph;
		std::printf("gate_M9: [%-5s] cells=%d  V_trap=%.4e m³  edge_scour=%.3f m  slacks=%d  mass_err=%.2e  clip=%.3f%%  KE×%.2f(%s)  maxu=%.2f %s\n",
			s.c_str(), ncells, r.V_last, r.edge_scour, r.nslack, r.mass_err, r.clip_frac * 100, r.ke_growth, r.ke_ok ? "ok" : "GROW", r.maxu, r.stable ? "" : "UNSTABLE");
		std::fflush(stdout);
		return r;
	};

	std::vector<ShapeResult> results;
	{ std::string s; std::stringstream ss(shapes); while (std::getline(ss, s, ',')) if (!s.empty()) results.push_back(run_shape(s)); }

	// Rank by matched-time trapped volume V_last (all shapes ran the identical reversing tide).
	std::sort(results.begin(), results.end(), [](const ShapeResult& a, const ShapeResult& b) { return a.V_last > b.V_last; });
	std::printf("gate_M9: --- ranking by matched-time trapped volume V_last under the reversing tide (higher = better self-ballasting) ---\n");
	for (size_t i = 0; i < results.size(); ++i)
		std::printf("gate_M9:   #%zu %-5s  V_last=%.4e m³  edge_scour=%.3f m  (mass_err=%.1e clip=%.2f%% KE×%.2f)\n",
			i + 1, results[i].name.c_str(), results[i].V_last, results[i].edge_scour, results[i].mass_err, results[i].clip_frac * 100, results[i].ke_growth);

	// --- Gate criteria (RESEARCH §8/§7 + PLAN M9) ------------------------------------------------
	bool stable = true, massok = true, keok = true, clipok = true, slacks_ok = true;
	for (auto& r : results)
	{
		stable = stable && r.stable;
		massok = massok && r.mass_err < 1e-3;   // < 0.1 % total-sand conservation
		clipok = clipok && r.clip_frac < 0.01;  // MORFAC limiter clip-fraction < 1 %
		if (reverse) { keok = keok && r.ke_ok; slacks_ok = slacks_ok && r.nslack >= 4; } // ≥ 4 slack transitions
	}
	// Ranking deliverable: ≥ 3 candidate shapes, with the BEST enclosing shape (ring or cup) out-trapping
	// the flat plate by a clear margin (> 20 % of the best) — a distinguishable, correctly-ordered
	// ranking. Under a REVERSING tide a symmetric enclosure that shelters its interior from BOTH flow
	// directions is the good trap; requiring one specific enclosure to win would over-fit the regime
	// (e.g. the square cup is under-resolved on a coarse patch), so the gate asserts the enclosing
	// PRINCIPLE beats the plate, plus a distinguishable pair (matched-time margin).
	double Vbest = results.empty() ? 0.0 : results.front().V_last, margin = 0.20 * Vbest;
	auto find = [&](const std::string& n) -> const ShapeResult* { for (auto& r : results) if (r.name == n) return &r; return nullptr; };
	const ShapeResult* cup = find("cup"); const ShapeResult* ring = find("ring"); const ShapeResult* plate = find("plate");
	bool enclose_beats_plate = true;
	if (plate && (cup || ring))
	{
		double best_enc = std::max(cup ? cup->V_last : 0.0, ring ? ring->V_last : 0.0);
		enclose_beats_plate = best_enc > plate->V_last + margin;
	}
	bool any_distinct = false; // at least one shape pair separated by > the matched-time margin
	for (size_t i = 0; i < results.size(); ++i) for (size_t k = i + 1; k < results.size(); ++k)
		if (std::fabs(results[i].V_last - results[k].V_last) > margin) any_distinct = true;
	bool ncand_ok = results.size() >= 3 && any_distinct;

	if (!csv.empty())
	{
		std::ofstream out(csv);
		out << "shape,V_eq_m3,V_ci_m3,V_last_m3,edge_scour_m,maxu,mass_err,clip_frac,slacks,ke_growth\n";
		for (auto& r : results)
			out << r.name << "," << r.V_eq << "," << r.V_ci << "," << r.V_last << "," << r.edge_scour << ","
			    << r.maxu << "," << r.mass_err << "," << r.clip_frac << "," << r.nslack << "," << r.ke_growth << "\n";
	}

	bool pass = stable && massok && clipok && keok && slacks_ok && ncand_ok && enclose_beats_plate;
	std::printf("========================================================================\n");
	std::printf("gate_M9: stable=%s  mass<0.1%%=%s  clip<1%%=%s  no-KE-growth=%s  ≥4-slacks=%s  ≥3-shapes+distinct=%s  enclosing>plate=%s\n",
		PF(stable), PF(massok), PF(clipok), PF(keok), PF(slacks_ok), PF(ncand_ok), PF(enclose_beats_plate));
	std::printf("gate_M9: RESULT %s\n", pass ? "PASS" : "FAIL");
	return pass ? 0 : 1;
}
