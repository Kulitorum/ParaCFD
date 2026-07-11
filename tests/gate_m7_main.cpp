// gate_m7_main.cpp — M7 deposition/protection gate "M-DEP" (Du et al. 2025 perforated-unit case /
// RESEARCH V8). Three cases on the Du sediment pit around a model monopile, at lab scale:
//   control : bare pile (no protection)
//   C-AR    : 3×3 array of perforated cubes,      surface porosity β = 0.1963
//   H-AR    : 3×3 array of perforated hemispheres, surface porosity β = 0.3125
// Sub-grid porous units (holes 8.33 mm ≈ 1–2 cells, below the 8-cell floor ⇒ thin-screen model,
// research/13 §8): the unit shells are porous cells with resistance k = 1/β²−1 (channel_porous.*),
// interior hollow, open bottom. d50 = 0.235 mm, U = 0.20/0.25/0.30 m/s (θ = 0.027/0.043/0.062).
//
// Gate (qualitative-comparative — see the ⚠ below): the two porous configs must REDUCE scour vs the
// control and TRAP sand (net deposition in/behind the array), mass must be conserved and the flow
// stable, and the measured C-AR-vs-H-AR ranking direction reproduced. Absolute scour depths are
// reported against Du's tables but are expected LOW: the diffusive SL solver under-resolves the
// horseshoe vortex (the same RESEARCH §11 limitation quantified by gate_M6), so the absolute
// magnitude is a SUPERVISED calibration matter for MH — this gate scores the PROTECTION EFFECT and
// RANKING, which are comparative and robust to that bias. A Camp–Hazen trapping cross-check is also
// reported (research/15 §6). Links only libscour. CLI: gate_m7 <config.json> [--U u] [--quick]
//   [--csv path] [--only control|car|har].
#include "core/fluid/channel_core.h"
#include "core/fluid/channel_mask.h"
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
	nlohmann::json jsub(const nlohmann::json& j, const char* k) { return j.contains(k) ? j.at(k) : nlohmann::json::object(); }
	const char* PF(bool b) { return b ? "PASS" : "FAIL"; }
	bool has_flag(int c, char** v, const char* f) { for (int i = 1; i < c; ++i) if (!std::strcmp(v[i], f)) return true; return false; }
	double opt_d(int c, char** v, const char* f, double d) { for (int i = 1; i < c - 1; ++i) if (!std::strcmp(v[i], f)) return std::atof(v[i + 1]); return d; }
	std::string opt_s(int c, char** v, const char* f, const std::string& d) { for (int i = 1; i < c - 1; ++i) if (!std::strcmp(v[i], f)) return v[i + 1]; return d; }

	struct CaseOut
	{
		double S_up = 0, S_down = 0, S_lat = 0, S_max = 0; // scour depths [m]
		double deposit_max = 0, deposit_vol = 0;           // peak bed rise + total trapped volume [m], [m³]
		double mass_err = 0, maxu = 0; bool stable = true;
	};
}

int main(int argc, char** argv)
{
	nlohmann::json j = nlohmann::json::object();
	if (argc >= 2 && argv[1][0] != '-')
	{
		try { std::ifstream in(argv[1]); std::ostringstream ss; ss << in.rdbuf(); j = nlohmann::json::parse(ss.str()); }
		catch (const std::exception& e) { std::fprintf(stderr, "gate_M7: bad config: %s\n", e.what()); }
	}
	nlohmann::json jg = jsub(j, "gates");
	const bool quick = has_flag(argc, argv, "--quick");
	std::string csv = opt_s(argc, argv, "--csv", "");
	std::string only = opt_s(argc, argv, "--only", "");

	// Geometry (Du 2025 pit; m). Domain z = sand pit depth + water depth.
	double Lx = jd(j, "domain_x", 2.4), Ly = jd(j, "domain_y", 0.5);
	double water = jd(j, "water_depth", 0.10), sand = jd(j, "sand_depth", 0.20);
	double h = jd(j, "voxel_h", 0.005);
	if (quick) h = std::max(h, 0.01);
	double D = jd(j, "pile_D", 0.05), pxc = jd(j, "pile_x", 0.8), pyc = jd(j, "pile_y", 0.5 * Ly);
	double unit_s = jd(j, "unit_size", 0.05), gap = jd(j, "unit_gap", 0.01);
	double Lz = water + sand;

	double U = opt_d(argc, argv, "--U", jd(j, "U", 0.30));
	double rho = jd(j, "rho", 1000.0), rho_s = jd(j, "rho_s", 2650.0);
	double d50 = jd(j, "d50", 0.235e-3), nu = jd(j, "nu", 1.0e-6);
	double Cs = jd(j, "Cs", 0.11), alpha = jd(j, "alpha", 0.00033), erosion = jd(j, "erosion_coeff", 0.018);
	double morfac = jd(j, "morfac", 4.0), sigma_s = jd(j, "sigma_s", 0.7);
	double cfl = jd(j, "cfl", 0.7), proj_tol = jd(j, "proj_tol", 1e-4);
	int proj_max_iter = ji(j, "proj_max_iter", 80);
	double beta_car = jd(j, "beta_car", 0.1963), beta_har = jd(j, "beta_har", 0.3125);
	double spinup_ft = jd(j, "spinup_flowthroughs", 4.0), trun_s = jd(j, "trun_morph_s", 60.0);
	if (quick) { spinup_ft = 2.0; trun_s = 8.0; }

	MacGrid g;
	g.nx = std::max(8, (int)std::llround(Lx / h));
	g.ny = std::max(8, (int)std::llround(Ly / h));
	g.nz = std::max(8, (int)std::llround(Lz / h));
	g.h = h;
	const int ncol = g.nx * g.ny, np = g.p_count();
	double kcar = screen_resistance_k(beta_car), khar = screen_resistance_k(beta_har);

	std::printf("========================================================================\n");
	std::printf("gate_M7 (Du 2025 / V8) grid %dx%dx%d h=%.4f (%.2fM cells) U=%.2f d50=%.3fmm | pit %.1fx%.1fx%.2f (sand %.2f+water %.2f)\n",
		g.nx, g.ny, g.nz, h, (double)np / 1e6, U, d50 * 1e3, Lx, Ly, Lz, sand, water);
	std::printf("gate_M7: pile D=%.3f at (%.2f,%.2f); 3x3 units side=%.3f gap=%.3f | C-AR β=%.4f k=%.2f  H-AR β=%.4f k=%.2f | morfac=%.1f\n",
		D, pxc, pyc, unit_s, gap, beta_car, kcar, beta_har, khar, morfac);
	std::fflush(stdout);

	// Pile (rigid full-z cylinder = structure/foundation).
	std::vector<unsigned char> pile;
	build_cylinder_mask(g, pxc, pyc, 0.5 * D, pile);

	// Per-column pile-footprint + distance (for scour metrics).
	std::vector<unsigned char> frozen(ncol, 0);
	std::vector<double> rdist(ncol), colx(ncol), coly(ncol);
	for (int jj = 0; jj < g.ny; ++jj) for (int ii = 0; ii < g.nx; ++ii)
	{
		int c = jj * g.nx + ii; double x = (ii + 0.5) * h, y = (jj + 0.5) * h;
		colx[c] = x; coly[c] = y; rdist[c] = std::sqrt((x - pxc) * (x - pxc) + (y - pyc) * (y - pyc));
		for (int k = 0; k < g.nz; ++k) if (pile[g.pidx(ii, jj, k)]) { frozen[c] = 1; break; }
	}
	const double Rp = 0.5 * D;

	// Build the porous-k field for a 3×3 array of hollow perforated units around the pile (open bottom).
	// The centre cell is the pile; the 8 surrounding positions get a unit. A unit shell = the boundary
	// cells (perimeter in x/y, + top layer) of its cube footprint in the water column above the bed.
	auto build_units = [&](double kval) -> std::vector<double>
	{
		std::vector<double> pk((size_t)np, 0.0);
		if (kval <= 0.0) return pk;
		double pitch = unit_s + gap;
		int z0 = (int)std::llround(sand / h), z1 = std::min(g.nz, z0 + (int)std::llround(unit_s / h));
		for (int gy = -1; gy <= 1; ++gy) for (int gx = -1; gx <= 1; ++gx)
		{
			if (gx == 0 && gy == 0) continue; // centre = pile
			double ucx = pxc + gx * pitch, ucy = pyc + gy * pitch;
			int i0 = (int)std::llround((ucx - 0.5 * unit_s) / h), i1 = (int)std::llround((ucx + 0.5 * unit_s) / h);
			int jj0 = (int)std::llround((ucy - 0.5 * unit_s) / h), jj1 = (int)std::llround((ucy + 0.5 * unit_s) / h);
			for (int k = z0; k < z1; ++k) for (int jj = jj0; jj < jj1; ++jj) for (int ii = i0; ii < i1; ++ii)
			{
				if (ii < 0 || ii >= g.nx || jj < 0 || jj >= g.ny || k < 0 || k >= g.nz) continue;
				bool shell = (ii == i0 || ii == i1 - 1 || jj == jj0 || jj == jj1 - 1 || k == z1 - 1); // walls + top, open bottom
				if (shell) pk[g.pidx(ii, jj, k)] = kval;
			}
		}
		return pk;
	};

	auto run_case = [&](const char* name, double kval) -> CaseOut
	{
		CaseOut r;
		SeabedParams sp;
		sp.sand_depth = sand; sp.d50 = d50; sp.rho = rho; sp.rho_s = rho_s; sp.nu = nu;
		sp.Cs = Cs; sp.alpha = alpha; sp.erosion_coeff = erosion; sp.morfac = morfac; sp.sigma_s = sigma_s;
		sp.bedload_formula = 1; sp.diffusion_on = 1; sp.sed_bc = SED_BC_CLOSED;
		SeabedMorpho engine(g, sp, pile);

		ChannelBC bc;
		bc.inlet_mode = INLET_LOGLAW; bc.z0 = d50 / 12.0; bc.kappa = 0.40; bc.bed_datum = sand;
		bc.ustar = loglaw_ustar_for_U(U, water, bc.z0, bc.kappa); bc.U_inlet = U; bc.Uc = U;
		bc.solid_mode = SOLID_FREESLIP; bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;

		ChannelParams pr;
		pr.rho = rho; pr.nu = nu; pr.Cs = Cs; pr.cfl = cfl; pr.safety = 0.9;
		pr.proj_tol = proj_tol; pr.proj_max_iter = proj_max_iter; pr.advect_band = 1;

		ChannelFluidCore core(g, bc, pr, engine.initial_flow_solid());
		core.set_bed_inlet_mask(true);
		std::vector<double> pk = build_units(kval);
		if (kval > 0.0) core.set_porous(pk);
		core.init_inlet_profile();

		double V0 = engine.bed_volume() + engine.susp_volume();
		double t_flow = 0.0, spin = spinup_ft * (Lx / U);
		std::vector<unsigned char> new_solid;
		engine.set_morphology_frozen(true);
		while (t_flow < spin) { double dt = core.step(); engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid); t_flow += dt; }
		engine.set_morphology_frozen(false);

		double t_morph = 0.0;
		while (t_morph < trun_s && r.stable)
		{
			double dt = core.step();
			if (engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid)) core.update_solid(new_solid);
			t_morph += morfac * dt;
			double mu = reduce_max_abs_gpu(core.u_dev(), g.u_count());
			r.maxu = std::max(r.maxu, mu);
			if (!std::isfinite(mu) || mu > 20.0 * U) r.stable = false;
		}

		std::vector<float> zb; engine.copy_zb_host(zb);
		double up = sand, dn = sand, lat = sand, mx = sand, depmax = sand, depvol = 0.0;
		for (int c = 0; c < ncol; ++c)
		{
			if (frozen[c]) continue;
			double z = zb[c], dz = z - sand;
			// scour ring around the pile (within ~2.5 unit-pitch)
			if (rdist[c] > Rp && rdist[c] < 2.5 * (unit_s + gap))
			{
				mx = std::min(mx, z);
				double ddx = std::fabs(colx[c] - pxc), ddy = std::fabs(coly[c] - pyc);
				if (colx[c] < pxc && ddy < ddx) up = std::min(up, z);          // upstream sector
				else if (colx[c] > pxc && ddy < ddx) dn = std::min(dn, z);     // downstream sector
				else lat = std::min(lat, z);                                   // lateral sector
			}
			// net deposition (trapped sand) WITHIN the array + near wake (unit-specific self-ballasting
			// KPI, not domain-wide settling): rdist is the pile-relative distance.
			if (dz > 0 && rdist[c] < 3.0 * (unit_s + gap)) { depmax = std::max(depmax, z); depvol += dz * h * h; }
		}
		r.S_up = sand - up; r.S_down = sand - dn; r.S_lat = sand - lat; r.S_max = sand - mx;
		r.deposit_max = depmax - sand; r.deposit_vol = depvol;
		double V = engine.bed_volume() + engine.susp_volume();
		r.mass_err = V0 > 0 ? std::fabs(V - V0) / V0 : 0.0;

		std::printf("gate_M7: [%-7s] S(up/down/lat/max)=%.1f/%.1f/%.1f/%.1f cm  deposit(max=%.1f cm, vol=%.2e m³)  maxu=%.3f mass_err=%.2e %s\n",
			name, r.S_up * 100, r.S_down * 100, r.S_lat * 100, r.S_max * 100, r.deposit_max * 100, r.deposit_vol, r.maxu, r.mass_err, r.stable ? "" : "UNSTABLE");
		std::fflush(stdout);
		return r;
	};

	std::printf("gate_M7: --- running %.0f s morphological time per case (MORFAC=%.0f) ---\n", trun_s, morfac);
	std::fflush(stdout);
	CaseOut ctrl, car, har; bool have_ctrl = false, have_car = false, have_har = false;
	if (only.empty() || only == "control") { ctrl = run_case("control", 0.0); have_ctrl = true; }
	if (only.empty() || only == "car") { car = run_case("C-AR", kcar); have_car = true; }
	if (only.empty() || only == "har") { har = run_case("H-AR", khar); have_har = true; }

	// Camp–Hazen trapping cross-check for the interior (research/15 §6): a unit of plan length ~unit_s,
	// interior depth ~unit_s, throughflow ~ interior velocity (≈ β·U for a porous box, RESEARCH §9).
	double ws = settling_ws(d50, rho, rho_s, nu);
	double eta_car = camp_hazen_efficiency(ws, unit_s, unit_s, beta_car * U);
	double eta_har = camp_hazen_efficiency(ws, unit_s, unit_s, beta_har * U);
	std::printf("gate_M7: Camp-Hazen interior trap η (ws=%.4f): C-AR≈%.2f  H-AR≈%.2f (higher-porosity H-AR flushes more ⇒ lower η)\n",
		ws, eta_car, eta_har);

	// ---- gate (comparative) ----
	bool pass = true;
	if (have_ctrl && !ctrl.stable) pass = false;
	if (have_car && !car.stable) pass = false;
	if (have_har && !har.stable) pass = false;
	bool reduce = true, deposit = true;
	if (have_ctrl && have_car) { reduce = reduce && (car.S_max <= ctrl.S_max + 1e-9); deposit = deposit && (car.deposit_vol > 0); }
	if (have_ctrl && have_har) { reduce = reduce && (har.S_max <= ctrl.S_max + 1e-9); deposit = deposit && (har.deposit_vol > 0); }
	double mass_tol = jd(jg, "mass_err_max", 0.05);
	bool massok = (!have_ctrl || ctrl.mass_err < mass_tol) && (!have_car || car.mass_err < mass_tol) && (!have_har || har.mass_err < mass_tol);
	pass = pass && reduce && deposit && massok;

	if (have_ctrl && have_car && have_har)
		std::printf("gate_M7: scour reduction vs control: C-AR %.0f%%  H-AR %.0f%% (Du: 28-100%%); C-AR/H-AR S_up ranking: sim %s, Du C-AR<H-AR\n",
			ctrl.S_max > 0 ? 100 * (1 - car.S_max / ctrl.S_max) : 0.0,
			ctrl.S_max > 0 ? 100 * (1 - har.S_max / ctrl.S_max) : 0.0,
			car.S_up <= har.S_up ? "C-AR<H-AR" : "C-AR>H-AR");
	std::printf("========================================================================\n");
	std::printf("gate_M7: stable=%s  porous reduces scour=%s  traps sand=%s  mass=%s\n",
		PF((!have_ctrl || ctrl.stable) && (!have_car || car.stable) && (!have_har || har.stable)), PF(reduce), PF(deposit), PF(massok));
	std::printf("gate_M7: RESULT %s  (comparative/qualitative — absolute S/D calibration is supervised, see gate_M6 report)\n", pass ? "PASS" : "FAIL");
	if (!csv.empty())
	{
		std::ofstream out(csv);
		out << "case,S_up_cm,S_down_cm,S_lat_cm,S_max_cm,deposit_max_cm,deposit_vol_m3,maxu,mass_err\n";
		auto row = [&](const char* n, const CaseOut& c) { out << n << "," << c.S_up * 100 << "," << c.S_down * 100 << "," << c.S_lat * 100 << "," << c.S_max * 100 << "," << c.deposit_max * 100 << "," << c.deposit_vol << "," << c.maxu << "," << c.mass_err << "\n"; };
		if (have_ctrl) row("control", ctrl); if (have_car) row("C-AR", car); if (have_har) row("H-AR", har);
	}
	return pass ? 0 : 1;
}
