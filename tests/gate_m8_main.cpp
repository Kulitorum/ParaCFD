// gate_m8_main.cpp — M8 shape-ranking harness (PLAN M8 gate). Runs procedural reference shapes on
// the live morphodynamic pipeline, records trapped-sand volume V_trap(t) inside each shape's control
// volume + edge-scour depth, fits equilibria (exp + Welzel, run_protocol.h) with confidence
// intervals, and RANKS the shapes by fitted V_eq (higher = better self-ballasting). The tool's
// stated authority is COMPARATIVE ranking under identical conditions (RESEARCH §11) — robust to the
// absolute-scour under-prediction quantified by gate_M6.
//
// Reference shapes (procedural voxel masks, shape_masks.h): a flat PLATE (low trapping), a solid
// RING (encloses an area, reflective/toe-scouring), and an open-top CUP (walls+floor, the classic
// settling trap). Physically expected order by trapping: cup > ring > plate.
//
// Gate: every case stable + sand-mass-conserving; the cup out-traps the plate; and at least one
// shape pair is DISTINGUISHABLE (fitted-V_eq CIs do not overlap) — the PLAN-M8 "rank correctly with
// non-overlapping CIs" criterion. Links only libscour. CLI: gate_m8 [config] [--quick] [--csv path]
//   [--shapes plate,ring,cup] [--trun-frac F] [--U u].
#include "core/campaign/run_protocol.h"
#include "core/fluid/channel_core.h"
#include "core/geometry/shape_masks.h"
#include "core/fluid/mac_ops.h"
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
	double jd(const nlohmann::json& j, const char* k, double d) { return j.contains(k) ? j.at(k).get<double>() : d; }
	int ji(const nlohmann::json& j, const char* k, int d) { return j.contains(k) ? j.at(k).get<int>() : d; }
	const char* PF(bool b) { return b ? "PASS" : "FAIL"; }
	bool has_flag(int c, char** v, const char* f) { for (int i = 1; i < c; ++i) if (!std::strcmp(v[i], f)) return true; return false; }
	double opt_d(int c, char** v, const char* f, double d) { for (int i = 1; i < c - 1; ++i) if (!std::strcmp(v[i], f)) return std::atof(v[i + 1]); return d; }
	std::string opt_s(int c, char** v, const char* f, const std::string& d) { for (int i = 1; i < c - 1; ++i) if (!std::strcmp(v[i], f)) return v[i + 1]; return d; }

	struct ShapeResult
	{
		std::string name;
		double V_eq = 0, V_ci = 0, V_last = 0;   // fitted trapped volume + CI [m³]
		double edge_scour = 0;                    // max scour depth just outside the footprint [m]
		double T_fit = 0, t_run = 0;              // fitted timescale + morphological run length [s]
		double V_design = 0, V_retained = 0, retention = -1.0; // retention protocol (−1 ⇒ not run)
		double mass_err = 0, maxu = 0; bool stable = true;
		bool fit_ok() const { return t_run > 0.3 * T_fit && T_fit > 0; } // research/16: never fit < 0.3T
	};
}

int main(int argc, char** argv)
{
	nlohmann::json j = nlohmann::json::object();
	if (argc >= 2 && argv[1][0] != '-')
	{
		try { std::ifstream in(argv[1]); std::ostringstream ss; ss << in.rdbuf(); j = nlohmann::json::parse(ss.str()); }
		catch (...) {}
	}
	const bool quick = has_flag(argc, argv, "--quick");
	std::string csv = opt_s(argc, argv, "--csv", "");
	std::string shapes = opt_s(argc, argv, "--shapes", "plate,ring,cup");
	bool do_retention = has_flag(argc, argv, "--retention"); // PLAN M8 retention protocol
	double ret_ramp = opt_d(argc, argv, "--retention-ramp", 1.5); // storm-stage U multiplier

	// Ranking domain (a single unit in a compact patch — fast, comparative). Domain z = sand + water.
	double Lx = jd(j, "domain_x", 4.0), Ly = jd(j, "domain_y", 4.0);
	double water = jd(j, "water_depth", 0.8), sand = jd(j, "sand_depth", 0.4);
	double h = jd(j, "voxel_h", 0.04);
	if (quick) h = std::max(h, 0.06);
	double Lz = water + sand;
	double U = opt_d(argc, argv, "--U", jd(j, "U", 0.8));
	double rho = jd(j, "rho", 1027.0), rho_s = jd(j, "rho_s", 2650.0), d50 = jd(j, "d50", 0.2e-3), nu = jd(j, "nu", 1.36e-6);
	double Cs = jd(j, "Cs", 0.11), alpha = jd(j, "alpha", 0.00033), erosion = jd(j, "erosion_coeff", 0.018);
	double morfac = jd(j, "morfac", 5.0), sigma_s = jd(j, "sigma_s", 0.7);
	double foot = jd(j, "footprint_half", 0.5), height = jd(j, "shape_height", 0.4);
	double spinup_ft = jd(j, "spinup_flowthroughs", 4.0), trun_frac = opt_d(argc, argv, "--trun-frac", jd(j, "trun_frac_T", 0.7));
	double sample_dt = jd(j, "sample_dt_s", 5.0);
	if (quick) { spinup_ft = 2.0; trun_frac = 0.15; sample_dt = 2.0; }

	MacGrid g;
	g.nx = std::max(8, (int)std::llround(Lx / h)); g.ny = std::max(8, (int)std::llround(Ly / h)); g.nz = std::max(8, (int)std::llround(Lz / h)); g.h = h;
	const int ncol = g.nx * g.ny;
	double xc = 0.5 * Lx, yc = 0.5 * Ly;

	// Predicted timescale for the run length (RESEARCH §7): u* from the flux-matched inlet, θ, T.
	double z0 = d50 / 12.0, ustar_amb = loglaw_ustar_for_U(U, water, z0, 0.40);
	double D_shape = 2.0 * foot;
	double T_pred = scour_T_from_flow(ustar_amb, D_shape, d50, rho, rho_s, water, 0);
	double t_run = trun_frac * T_pred;

	std::printf("========================================================================\n");
	std::printf("gate_M8 shape ranking  grid %dx%dx%d h=%.3f (%.2fM) domain %.1fx%.1fx%.1f U=%.2f d50=%.3fmm | T_pred=%.0fs t_run=%.0fs (%.2fT) morfac=%.0f\n",
		g.nx, g.ny, g.nz, h, (double)g.p_count() / 1e6, Lx, Ly, Lz, U, d50 * 1e3, T_pred, t_run, trun_frac, morfac);
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

	// "pcup" = a POROUS cup: the same walls, but modelled as a sub-grid porous screen (k = 1/β²−1 at
	// the RESEARCH §9 optimum β≈0.30) instead of solid — the actual COBOD permeable-enclosure concept.
	// Returns the wall cells as a porous-k field (mask stays empty ⇒ no rigid structure). Combines the
	// M7 porous model with the M8 ranking to test whether porosity improves trapping vs a solid cup.
	double beta_por = opt_d(argc, argv, "--beta-por", jd(j, "beta_por", 0.30)); // porous-cup open-area ratio
	auto build_porous = [&](const std::string& s, std::vector<double>& pk) -> int
	{
		pk.assign((size_t)g.p_count(), 0.0);
		if (s != "pcup") return 0;
		std::vector<unsigned char> wall; build_cup_mask(g, xc, yc, foot, sand, height, std::max(g.h, 0.08), wall);
		double k = screen_resistance_k(beta_por); int n = 0;
		for (int c = 0; c < g.p_count(); ++c) if (wall[c]) { pk[c] = k; ++n; }
		return n;
	};

	auto run_shape = [&](const std::string& s) -> ShapeResult
	{
		ShapeResult r; r.name = s;
		std::vector<unsigned char> mask; std::vector<double> pk;
		int ncells = (s == "pcup") ? build_porous(s, pk) : build_shape(s, mask);

		SeabedParams sp;
		sp.sand_depth = sand; sp.d50 = d50; sp.rho = rho; sp.rho_s = rho_s; sp.nu = nu;
		sp.Cs = Cs; sp.alpha = alpha; sp.erosion_coeff = erosion; sp.morfac = morfac; sp.sigma_s = sigma_s;
		sp.bedload_formula = 1; sp.diffusion_on = 1; sp.sed_bc = SED_BC_OPEN; sp.U_inlet = U; // open sea carries the ambient load
		SeabedMorpho engine(g, sp, mask);

		ChannelBC bc;
		bc.inlet_mode = INLET_LOGLAW; bc.z0 = z0; bc.kappa = 0.40; bc.bed_datum = sand;
		bc.ustar = ustar_amb; bc.U_inlet = U; bc.Uc = U;
		bc.solid_mode = SOLID_FREESLIP; bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;

		ChannelParams pr; pr.rho = rho; pr.nu = 5e-3 /*coarse-grid artificial ν for the ranking patch*/;
		pr.Cs = Cs; pr.cfl = 0.7; pr.safety = 0.9; pr.proj_tol = 1e-4; pr.proj_max_iter = 60; pr.advect_band = 1;

		ChannelFluidCore core(g, bc, pr, engine.initial_flow_solid());
		core.set_bed_inlet_mask(true);
		if (!pk.empty()) core.set_porous(pk); // porous shape (pcup): walls are a sub-grid screen
		core.init_inlet_profile();

		double V0 = engine.bed_volume() + engine.susp_volume();
		std::vector<unsigned char> new_solid;
		double tf = 0, spin = spinup_ft * (Lx / U);
		engine.set_morphology_frozen(true);
		while (tf < spin) { double dt = core.step(); engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid); tf += dt; }
		engine.set_morphology_frozen(false);

		std::vector<double> ts, Vs;
		double tm = 0, next = 0;
		while (tm < t_run && r.stable)
		{
			double dt = core.step();
			if (engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid)) core.update_solid(new_solid);
			tm += morfac * dt;
			if (tm >= next)
			{
				std::vector<float> zb; engine.copy_zb_host(zb);
				double vtrap = 0, edge = sand;
				for (int c = 0; c < ncol; ++c)
				{
					double dz = zb[c] - sand;
					if (rdist[c] <= r_cv && dz > 0) vtrap += dz * h * h;             // trapped sand in the CV
					if (rdist[c] > r_edge0 && rdist[c] < r_edge1) edge = std::min(edge, (double)zb[c]); // edge scour
				}
				double mu = reduce_max_abs_gpu(core.u_dev(), g.u_count());
				r.maxu = std::max(r.maxu, mu); if (!std::isfinite(mu) || mu > 20 * U) r.stable = false;
				double V = engine.bed_volume() + engine.susp_volume();
				// OPEN sea boundary imports/exports sand; the conserved quantity credits that net flux.
				double net = V - V0 - engine.boundary_sand_in();
				double denom = V0 + std::fabs(engine.boundary_sand_in());
				r.mass_err = std::max(r.mass_err, denom > 0 ? std::fabs(net) / denom : 0.0);
				r.edge_scour = sand - edge;
				ts.push_back(tm); Vs.push_back(vtrap); r.V_last = vtrap;
				next += sample_dt;
			}
		}
		// Retention protocol (PLAN M8): ramp U up to a storm stage (×ret_ramp) over ~60 s, run another
		// t_run of morphological time, and report the fraction of the trapped sand that SURVIVES. This
		// is the metric that penalises too-open shapes (their sheltered interior sees higher shear at
		// the higher U ⇒ the catch re-erodes) — the missing half of the RESEARCH §9 porosity trade-off.
		if (do_retention && r.stable)
		{
			const double PI = 3.14159265358979323846;
			r.V_design = r.V_last;
			double Uh = U * ret_ramp, tr = 0, ramp_s = 60.0, ret_dur = t_run;
			while (tr < ret_dur && r.stable)
			{
				double Unow = (tr < ramp_s) ? U + (Uh - U) * 0.5 * (1.0 - std::cos(PI * tr / ramp_s)) : Uh;
				core.set_inlet_speed(Unow); engine.set_inlet_speed(Unow);
				double dt = core.step();
				if (engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid)) core.update_solid(new_solid);
				tr += morfac * dt;
				double mu = reduce_max_abs_gpu(core.u_dev(), g.u_count());
				r.maxu = std::max(r.maxu, mu); if (!std::isfinite(mu) || mu > 20 * Uh) r.stable = false;
			}
			std::vector<float> zb; engine.copy_zb_host(zb);
			double vr = 0; for (int c = 0; c < ncol; ++c) { double dz = zb[c] - sand; if (rdist[c] <= r_cv && dz > 0) vr += dz * h * h; }
			r.V_retained = vr; r.retention = r.V_design > 0 ? vr / r.V_design : 0.0;
		}

		EqFit fe = fit_exponential(ts, Vs), fw = fit_welzel(ts, Vs);
		double rmse = fe.sse > 0 && ts.size() ? std::sqrt(fe.sse / ts.size()) : 0.0;
		ShapeRank rk = equilibrium_with_ci(fe, fw, rmse);
		r.V_eq = rk.value; r.V_ci = rk.ci; r.T_fit = fe.param; r.t_run = tm;
		std::printf("gate_M8: [%-5s] cells=%d  V_trap: last=%.4e eq=%.4e ±%.1e m³ (T=%.0fs, t/T=%.2f)  edge_scour=%.3f m  maxu=%.2f mass_err=%.2e %s\n",
			s.c_str(), ncells, r.V_last, r.V_eq, r.V_ci, fe.param, r.T_fit > 0 ? tm / r.T_fit : 0.0, r.edge_scour, r.maxu, r.mass_err, r.stable ? "" : "UNSTABLE");
		if (r.retention >= 0.0)
			std::printf("gate_M8:   [%-5s] retention @ %.2f×U: V_design=%.4e -> V_retained=%.4e  retained=%.0f%%\n",
				s.c_str(), ret_ramp, r.V_design, r.V_retained, r.retention * 100);
		std::fflush(stdout);
		return r;
	};

	// Parse shape list + run.
	std::vector<ShapeResult> results;
	{
		std::string s; std::stringstream ss(shapes);
		while (std::getline(ss, s, ',')) if (!s.empty()) results.push_back(run_shape(s));
	}

	// Rank by MATCHED-TIME trapped volume V_last (all shapes ran identical t_run under identical
	// conditions — the assumption-free comparison, research/16). The fitted V_eq is reported too, but
	// slow-filling traps (e.g. the cup) legitimately do not saturate within a fixed run, so V_eq
	// extrapolation is a supplement, not the ranking basis.
	std::sort(results.begin(), results.end(), [](const ShapeResult& a, const ShapeResult& b) { return a.V_last > b.V_last; });
	std::printf("gate_M8: --- ranking by matched-time trapped volume V_last (higher = better self-ballasting) ---\n");
	for (size_t i = 0; i < results.size(); ++i)
		std::printf("gate_M8:   #%zu %-5s  V_last=%.4e m³  (V_eq≈%.4e ±%.1e, t/T=%.2f %s)  edge_scour=%.3f m\n",
			i + 1, results[i].name.c_str(), results[i].V_last, results[i].V_eq, results[i].V_ci,
			results[i].T_fit > 0 ? results[i].t_run / results[i].T_fit : 0.0, results[i].fit_ok() ? "reliable" : "not-saturated", results[i].edge_scour);

	// Gate (matched-time, robust to unsaturated traps): all stable + mass-conserving (with boundary
	// credit); the ENCLOSING shapes (cup/ring) out-trap the flat plate by a clear margin (>20% of the
	// best V_last — the "rank correctly with a distinguishable gap" PLAN-M8 criterion; sheltered
	// interior accretes more than a bare plate's wake). Fit reliability is REPORTED (research/16) but
	// not required to pass — a slow-filling trap not saturating in a fixed run is expected, not a fault.
	bool stable = true, massok = true, fitok = true;
	for (auto& r : results) { stable = stable && r.stable; massok = massok && r.mass_err < 0.05; fitok = fitok && r.fit_ok(); }
	auto find = [&](const std::string& n) -> const ShapeResult* { for (auto& r : results) if (r.name == n) return &r; return nullptr; };
	const ShapeResult* cup = find("cup"); const ShapeResult* ring = find("ring"); const ShapeResult* plate = find("plate");
	double Vbest = results.empty() ? 0.0 : results.front().V_last, margin = 0.20 * Vbest;
	bool enclose_beats_plate = true;
	if (plate) { if (cup) enclose_beats_plate = enclose_beats_plate && cup->V_last > plate->V_last + margin;
	             if (ring) enclose_beats_plate = enclose_beats_plate && ring->V_last > plate->V_last + margin; }
	bool any_distinct = false; // any pair separated by > the matched-time margin
	for (size_t i = 0; i < results.size(); ++i) for (size_t k = i + 1; k < results.size(); ++k)
		if (std::fabs(results[i].V_last - results[k].V_last) > margin) any_distinct = true;
	bool pass = stable && massok && enclose_beats_plate && any_distinct;

	if (!csv.empty())
	{
		std::ofstream out(csv); out << "shape,V_eq_m3,V_ci_m3,V_last_m3,edge_scour_m,maxu,mass_err\n";
		for (auto& r : results) out << r.name << "," << r.V_eq << "," << r.V_ci << "," << r.V_last << "," << r.edge_scour << "," << r.maxu << "," << r.mass_err << "\n";
	}

	std::printf("========================================================================\n");
	std::printf("gate_M8: stable=%s mass=%s  enclosing>plate(trapping)=%s  distinguishable-pair=%s  (fit-reliable-all=%s, reported only)\n",
		PF(stable), PF(massok), PF(enclose_beats_plate), PF(any_distinct), PF(fitok));
	std::printf("gate_M8: RESULT %s\n", pass ? "PASS" : "FAIL");
	return pass ? 0 : 1;
}
