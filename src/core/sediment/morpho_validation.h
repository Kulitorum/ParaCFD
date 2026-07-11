// morpho_validation.h — M5 acceptance-gate drivers (RESEARCH §5, §7; PLAN M5 gate a–d). Each
// driver runs the SHIPPED morphodynamics kernels (bedstate.*, avalanche.*, suspended.*) on the GPU
// so the gate exercises the real code path, not a bespoke reimplementation.
//
//   (a-i)  still-water deposition: column grain rises at w_s·c0 (bed elevation w_s·c0/0.64) < 0.5%,
//          three-reservoir sand-mass error < 0.1%.
//   (a-ii) prescribed bedload q_b(x)=q0·sin(2πx/L), E=D=0: (1−p)∂z_b/∂t = −∂q_b/∂x (upwind ∇·q_b)
//          matched < 0.5% after 100 updates.
//   (b)    sand-pile relaxes to the 30–32° repose band everywhere (mass conserved).
//   (c)    three-reservoir mass (bed + suspended + boundary) conserved < 0.1% over 1e4 steps.
//   (d)    no motion below Soulsby–Whitehouse θ_cr; onset within 10%.
#pragma once

#include <vector>

namespace scour::core
{
	// ---- (a-i) still-water deposition column ----------------------------------------------------
	struct DepositionColumnConfig
	{
		double Lz = 1.0, h = 0.01;
		int nx = 4, ny = 4;
		double d50 = 0.2e-3, rho = 1027.0, rho_s = 2650.0, nu = 1.36e-6;
		double c0 = 0.004;    // uniform initial volumetric concentration (dilute ⇒ D≈w_s·c0 < 0.5%)
		double cfl = 0.7;     // settling Courant w_s·dt/h
	};
	struct DepositionColumnResult
	{
		int nx = 0, ny = 0, nz = 0;
		double ws = 0, dt = 0, t_check = 0;
		double G_num = 0;        // deposited grain thickness Σf_pack·h at t_check [m]
		double G_expected = 0;   // w_s·c0·t_check [m]
		double rate_err = 0;     // |G_num − G_expected|/G_expected           (gate: < 0.5%)
		double zb_rate = 0;      // bed-elevation rise rate w_s·c0/0.64 [m/s]
		double mass_err = 0;     // max |V_bed+V_susp − V0|/V0 over the run    (gate: < 0.1%)
	};
	DepositionColumnResult run_deposition_column(const DepositionColumnConfig& cfg);

	// ---- (a-ii) prescribed sinusoidal bedload vs 1D Exner --------------------------------------
	struct BedloadExnerConfig
	{
		double L = 10.0;      // domain length (one full sine period) [m]
		int nx = 1024;        // columns (upwind O(h) ⇒ nx≳600 for < 0.5% L2)
		double q0 = 1.0e-6;   // bedload amplitude [m²/s]
		double porosity = 0.36;
		int updates = 100;    // morphology updates
		double dt = 0.0;      // per-update dt [s]; 0 ⇒ auto (keeps Δz_b/step ≪ limiter)
	};
	struct BedloadExnerResult
	{
		int nx = 0, updates = 0;
		double dt = 0, l2_err = 0; // ‖Δz_b,num − Δz_b,exact‖₂/‖Δz_b,exact‖₂   (gate: < 0.5%)
		double dz_max = 0;         // max |Δz_b| over the field [m]
	};
	BedloadExnerResult run_bedload_exner(const BedloadExnerConfig& cfg);

	// ---- (b) sand-pile avalanche relaxation ---------------------------------------------------
	struct SandpileConfig
	{
		int nx = 61, ny = 61;
		double h = 0.05;
		double peak = 2.0;          // initial cone peak elevation [m]
		double init_slope_deg = 45; // initial cone flank slope
		double phi_trigger_deg = 32, phi_repose_deg = 30;
		double relax = 0.1;
		int max_sweeps = 40000;
	};
	struct SandpileResult
	{
		int sweeps = 0;
		double init_max_slope_deg = 0, final_max_slope_deg = 0;
		double mass_err = 0;   // |ΣG_end − ΣG_0|/ΣG_0  (avalanche is mass-conserving)
		bool converged = false;
	};
	SandpileResult run_sandpile(const SandpileConfig& cfg);

	// ---- (c) three-reservoir mass conservation over 1e4 steps ---------------------------------
	struct MassConservationConfig
	{
		int nx = 16, ny = 16, nz = 20;
		double h = 0.05;
		double d50 = 0.35e-3, rho = 1027.0, rho_s = 2650.0, nu = 1.36e-6;
		double c0 = 0.002;     // initial suspended concentration in the fluid column
		double tau_peak = 0.0; // peak grain-skin shear [Pa]; 0 ⇒ auto (≈2·τ_cr, live-bed)
		int steps = 10000;
		double cfl = 0.5;
	};
	struct MassConservationResult
	{
		int steps = 0;
		double dt = 0, V0 = 0;
		double mass_err = 0;   // max |V_bed+V_susp − V0|/V0 over the run  (gate: < 0.1%)
		double vbed_frac_end = 0, vsusp_frac_end = 0;
	};
	MassConservationResult run_mass_conservation(const MassConservationConfig& cfg);

	// ---- (d) threshold of motion (Soulsby–Whitehouse) -----------------------------------------
	struct ThresholdConfig
	{
		double d50 = 0.2e-3, rho = 1025.0, rho_s = 2650.0, nu = 1.05e-6; // 20 °C ⇒ θ_cr≈0.0492
		double theta_lo_frac = 0.5, theta_hi_frac = 1.5; // sweep window in units of θ_cr,flat
		int nsweep = 200;
		double h = 0.05;
		int erode_steps = 20;
	};
	struct ThresholdResult
	{
		double theta_cr_flat = 0;
		double theta_onset = 0;   // smallest θ producing bed motion
		double onset_err = 0;     // |θ_onset − θ_cr,flat|/θ_cr,flat        (gate: < 10%)
		bool no_motion_below = false; // ΔG==0 for every θ ≤ θ_cr,flat      (gate)
		bool bedload_gate_ok = false; // q_b=0 at 0.9·θ_cr, q_b>0 at 1.1·θ_cr
	};
	ThresholdResult run_threshold(const ThresholdConfig& cfg);
}
