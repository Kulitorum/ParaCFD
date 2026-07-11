// sed_validation.h — M4 acceptance-gate drivers (V5 settling column, V6 Rouse profile).
// Both run the real suspended-sediment kernels (suspended.*) on the GPU so the gate exercises
// the shipped code path, not a bespoke reimplementation. RESEARCH §6.1–6.2, research/02 §7,
// research/11 §6–7. See sed_validation.cu for the exact analytic references each gate checks.
#pragma once

#include <vector>

namespace scour::core
{
	// ---- V5: still-water settling column ------------------------------------------------------
	// A smooth (Gaussian) sediment slab released in quiescent water settles at w_s. The physically
	// correct behaviour of a RESOLVED still-water column is plug flow: the slab TRANSLATES downward
	// at exactly w_s with mass conserved, then deposits on the bed via the flux BC D = w_s·c_b
	// (research/11 §6: "a concentration front advects down at exactly w_s; mass conservation").
	// [NB: RESEARCH §10's c(t)=c0·e^(−w_s t/h) is the *well-mixed lumped* idealisation; it does not
	//  describe a resolved diffusion-free column — see sed_validation.cu. Both are reported.]
	struct SettlingColumnConfig
	{
		double Lz = 1.0;       // column height [m]
		double h = 0.005;      // voxel edge [m] (nz = Lz/h)
		int nx = 4, ny = 4;    // homogeneous transverse strip
		double d50 = 0.2e-3;   // grain size [m]
		double rho = 1027.0, rho_s = 2650.0, nu = 1.36e-6; // RESEARCH §1 (config)
		double c0 = 0.01;      // slab peak volumetric concentration
		double z_center = 0.7; // initial slab center height [m]
		double sigma = 0.06;   // initial slab Gaussian half-width [m]
		double cfl = 0.7;      // settling advective Courant w_s·dt/h
		double travel_check = 0.30; // airborne descent [m] at which the L2/translate check is taken
	};
	struct SettlingColumnResult
	{
		int nx = 0, ny = 0, nz = 0;
		double ws = 0;               // Soulsby w_s(d50) [m/s]
		double dt = 0;
		double t_check = 0;          // time of the airborne L2 checkpoint [s]
		double l2_translate = 0;     // ‖c − c_exact_translate‖₂ / ‖c_exact‖₂ at the checkpoint  (gate: <2%)
		double v_centroid = 0;       // fitted centroid descent speed [m/s]
		double v_centroid_err = 0;   // |v_centroid − w_s| / w_s                                  (settling correctness)
		double mass_err = 0;         // max |M_susp + M_bed − M0| / M0 over the whole run (settling + deposition; gate <0.1%)
		double susp_frac_end = 0;    // suspended fraction remaining at end (→ ~0 as all deposits)
		std::vector<double> zc, c_num, c_exact; // checkpoint profile (column i=0)
	};
	SettlingColumnResult run_settling_column(const SettlingColumnConfig& cfg);

	// ---- V6: equilibrium Rouse profile in a periodic channel ----------------------------------
	// Column of height (h_dom − a) with the eddy diffusivity ε_s = β·κ·u*·z·(1−z/h_dom) imposed as
	// D_c = ν_t/σ_s (β=1 ⇒ σ_s=1 to match the classical Rouse analytic; production default σ_s=0.7).
	// Settling folded into advection + this diffusion relaxes to  c/c_a = [((h−z)/z)(a/(h−a))]^R,
	// R = w_s/(β·κ·u*)  (research/11 §7: h=0.4, u*=0.02, d50=0.1 mm ⇒ R=0.91). c is Dirichlet-pinned
	// to c_a at a=0.05h (permitted in THIS validation only; production bed exchange stays flux-based).
	struct RouseChannelConfig
	{
		double h_dom = 0.4;    // flow depth [m]
		double u_star = 0.02;  // friction velocity [m/s]
		double d50 = 0.1e-3;   // grain size [m]
		double rho = 1025.0, rho_s = 2650.0, nu = 1.05e-6; // 20 °C reference (research/11 §7 ⇒ R=0.91)
		double kappa = 0.40;   // von Kármán (RESEARCH §4)
		double beta = 1.0;     // sediment/momentum diffusivity ratio (Rouse analytic uses β=1)
		double a_frac = 0.05;  // reference height a = a_frac·h_dom
		double c_a = 1.0e-3;   // reference (pinned) concentration; small ⇒ hindered settling negligible
		int nz = 80;           // vertical cells over [a, h_dom]
		int nx = 4, ny = 4;
		double tol = 1e-6;     // steady-state relative-change threshold
		int steps_max = 400000;
	};
	struct RouseChannelResult
	{
		int nz = 0, steps = 0;
		double ws = 0;               // Soulsby w_s(d50) [m/s]
		double R_target = 0;         // w_s/(β·κ·u*)  (computed from config; ≈0.91 at the benchmark ν)
		double R_fit = 0;            // slope of ln c vs ln((h−z)/z) over the fit window
		double R_err = 0;            // |R_fit − R_target| / R_target                              (gate: <15%)
		double fit_r2 = 0;
		double dt = 0, sim_time = 0;
		bool converged = false;
		std::vector<double> zc, c_num, c_rouse; // steady profile + analytic (scaled to the same c_a)
	};
	RouseChannelResult run_rouse_channel(const RouseChannelConfig& cfg);
}
