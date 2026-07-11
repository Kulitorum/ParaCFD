// channel_periodic.h — M3 periodic-channel driver (RESEARCH §4 anchor, §8; research/14
// §4, research/16). A flat-bed open-channel strip, periodic in x/y, driven by a constant
// body force g_x = u*²/h_dom with a PI mass-flux controller locking U_bulk = U_d
// (research/14 §4: Baykal-type precursor). Bed friction is supplied by the log-law wall
// model (bedshear.*) applied as an equilibrium momentum sink; the interior turbulence
// closure is a Prandtl mixing length l = κ·z·√(1−z/h_dom) — the standard open-channel
// eddy viscosity that yields a depth-filling log law (grid-tied Smagorinsky Δ=h cannot,
// research/04). See channel_periodic.cpp for the full rationale and the exact steady-state
// momentum balance the gate checks.
//
// Gate V2 (PLAN M3): extract u* within 10% of κ·U_d/(ln(h_dom/z0)−1), z0 = d50/12; recover
// κ within 5% from a log-law fit of the mean profile; τ_b smooth (flat bed → uniform).
#pragma once

#include "core/fluid/bedshear.h"

#include <vector>

namespace windcfd::core
{
	struct PeriodicChannelConfig
	{
		double h = 0.05;      // voxel [m]
		double Lz = 5.0;      // domain height h_dom [m]
		int nx = 8, ny = 8;   // transverse extent (homogeneous — kept small)
		double U_d = 1.0;     // target depth-averaged speed [m/s]
		double d50 = 0.2e-3;  // grain size [m] → z0 = d50/12, ks = 2.5·d50
		double nu = 1.36e-6;  // molecular viscosity [m²/s] (RESEARCH §1)
		double rho = 1027.0;  // density [kg/m³]
		double kappa = 0.40;  // von Kármán (RESEARCH §4)
		int wall_regime = WALL_AUTO;
		double pi_gain = 0.1; // PI mass-flux controller gain (research/14 §4)
		double cfl_diff = 0.6;// explicit z-diffusion safety
		double tol = 1e-4;    // steady-state tol: |U_bulk−U_d|/U_d lock threshold
		int steps_max = 400000;
		int report_every = 0; // 0 = silent; else stdout progress cadence
	};

	struct PeriodicChannelResult
	{
		int nx = 0, ny = 0, nz = 0, steps = 0;
		double h = 0, z0 = 0, ks_plus = 0, sim_time = 0;
		double ustar_analytic = 0; // κ·U_d/(ln(Lz/z0)−1), z0 = d50/12  (rough reference)
		double ustar_wall = 0;     // √(mean|τ_b|/ρ) from the EMA wall model
		double ustar_momentum = 0; // √(g_x·Lz)  (streamwise momentum balance)
		double kappa_fit = 0;      // recovered from the log-law profile fit
		double U_bulk = 0, gx = 0;
		double tau_band = 0;       // (max−min)/mean of |τ_b| over the bed plane (banding)
		double max_div = 0;        // max|∇·u| (1D flow ⇒ ~machine ε; proves projection moot)
		double fit_r2 = 0;         // R² of the log-law profile fit
		std::vector<double> zc, uc; // mean profile (cell-centred)
		bool converged = false;
	};

	PeriodicChannelResult run_periodic_channel(const PeriodicChannelConfig& cfg);
}
