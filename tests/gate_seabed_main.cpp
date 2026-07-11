// gate_seabed_main.cpp — headless end-to-end verification of the live erodible-seabed
// morphodynamic loop (fluid → τ_b wall model → suspended sediment → bed Exner → avalanche),
// with NO Qt/GL. Two cases on a 0.5 m sand bed under a current:
//   A) flat bed, no structure — must stay STABLE (no blow-up) and MASS-CONSERVING (bed+suspended).
//   B) a procedural box structure seated on the bed — must stay stable/conserving AND develop a
//      scour signature (the bed LOWERS somewhere: min z_b drops below z_b0).
// This is the PRE-CALIBRATION demo pipeline exercised headlessly (the GUI run adds the visuals).
// Constants are NOMINAL/uncalibrated — qualitative only. Links only libscour.
#include "core/fluid/channel_core.h"
#include "core/fluid/mac_ops.h" // reduce_max_abs_gpu
#include "core/sediment/seabed_engine.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace scour::core;

namespace
{
	struct CaseResult
	{
		bool stable = true;
		double mass_err = 0.0, maxu = 0.0;
		double zb0 = 0.0, zb_min = 0.0, zb_max = 0.0;
	};

	CaseResult run_case(const char* name, bool with_structure)
	{
		CaseResult r;
		MacGrid g; g.nx = 60; g.ny = 40; g.nz = 40; g.h = 0.05; // 96k cells; domain 3.0×2.0×2.0 m
		const double U = 0.5, sand_depth = 0.5;

		SeabedParams sp;
		sp.sand_depth = sand_depth; sp.d50 = 0.2e-3; sp.rho = 1027.0; sp.rho_s = 2650.0; sp.nu = 1.36e-6;
		sp.Cs = 0.11; sp.alpha = 0.00033; sp.morfac = 5.0; sp.bedload_formula = 1; sp.diffusion_on = 1;

		// A modest obstacle seated on the bed (a few % of the span, like a printed unit in the wide
		// production domain). A very bluff obstacle over-constricts the coarse M2 core — that regime
		// is M6's calibrated turbulent flow, not this qualitative demo.
		std::vector<unsigned char> structure;
		if (with_structure)
		{
			structure.assign((size_t)g.p_count(), 0);
			int i0 = 26, i1 = 30, j0 = 18, j1 = 22;                         // 4×4 cells in plan
			int k0 = (int)std::llround(sand_depth / g.h), k1 = k0 + 5;      // 5 cells tall, on the bed
			for (int k = k0; k < k1 && k < g.nz; ++k) for (int j = j0; j < j1; ++j) for (int i = i0; i < i1; ++i)
				structure[g.pidx(i, j, k)] = 1;
		}

		SeabedMorpho engine(g, sp, structure);

		ChannelBC bc;
		bc.inlet_mode = INLET_UNIFORM; bc.U_inlet = U; bc.Uc = U;
		bc.solid_mode = SOLID_NOSLIP;
		bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;

		ChannelParams pr;
		// Artificial fluid viscosity (ν = U·L/Re, Re≈150 on the obstacle scale) keeps the coarse-grid
		// LES flow stable — the g1_viewer_full recipe. FLOW closure only; the SEDIMENT uses molecular ν
		// (sp.nu) for its grain physics. Qualitative demo, not a calibrated Re.
		pr.rho = sp.rho; pr.nu = 5e-3; pr.Cs = sp.Cs; pr.cfl = 0.7; pr.safety = 0.9;
		pr.proj_tol = 1e-4; pr.proj_max_iter = 60; pr.advect_band = 1;

		ChannelFluidCore core(g, bc, pr, engine.initial_flow_solid());
		core.set_bed_inlet_mask(true); // the erodible bed reaches the inlet plane (no buried inflow)
		core.init_uniform(U, 0.0);

		const int spinup = 120, nsteps = 500;
		double V0 = engine.bed_volume() + engine.susp_volume();
		std::vector<unsigned char> new_solid;
		int remasks = 0;

		for (int s = 0; s < nsteps; ++s)
		{
			double dt = core.step();
			if (s >= spinup && engine.step(core.u_dev(), core.v_dev(), core.w_dev(), core.nut_dev(), dt, &new_solid))
			{
				core.update_solid(new_solid);
				++remasks;
			}
			if (s % 50 == 0 || s == nsteps - 1)
			{
				double maxu = reduce_max_abs_gpu(core.u_dev(), g.u_count());
				r.maxu = std::max(r.maxu, maxu);
				double V = engine.bed_volume() + engine.susp_volume();
				r.mass_err = std::max(r.mass_err, V0 > 0 ? std::fabs(V - V0) / V0 : 0.0);
				if (!std::isfinite(maxu) || maxu > 20.0 * U) r.stable = false;
			}
		}

		std::vector<float> zb; engine.copy_zb_host(zb);
		r.zb0 = sand_depth; r.zb_min = engine.zb_min(); r.zb_max = engine.zb_max();

		std::printf("[seabed:%s] steps=%d spinup=%d remasks=%d morfac=%.1f | maxu=%.3f m/s | "
			"mass_err=%.3e | z_b: min=%.4f max=%.4f (z_b0=%.3f) | max_ustar=%.4f\n",
			name, nsteps, spinup, remasks, sp.morfac, r.maxu, r.mass_err, r.zb_min, r.zb_max, r.zb0, engine.max_ustar());
		return r;
	}
}

int main()
{
	std::printf("PRE-CALIBRATION DEMO — constants nominal, results qualitative\n");

	CaseResult a = run_case("flat", false);
	CaseResult b = run_case("box", true);

	// Flat bed: stable, mass-conserving, and the bed stays close to flat (only mild features).
	bool passA = a.stable && (a.mass_err < 5e-3) && (std::fabs(a.zb_max - a.zb0) < 0.15) && (std::fabs(a.zb0 - a.zb_min) < 0.15);
	// Box: stable, conserving (looser — the rigid footprint pin trades a little conservation), and a
	// scour signature (the bed lowers near the obstacle).
	bool passB = b.stable && (b.mass_err < 5e-2) && (b.zb_min < b.zb0 - 0.002);

	std::printf("seabed-flat: %s (stable=%d mass_err=%.2e<5e-3 flatness ok)\n",
		passA ? "PASS" : "FAIL", (int)a.stable, a.mass_err);
	std::printf("seabed-box:  %s (stable=%d mass_err=%.2e<5e-2 scour: z_b_min=%.4f < %.4f)\n",
		passB ? "PASS" : "FAIL", (int)b.stable, b.mass_err, b.zb_min, b.zb0 - 0.002);

	bool pass = passA && passB;
	std::printf("gate_seabed: %s\n", pass ? "PASS" : "FAIL");
	return pass ? 0 : 1;
}
