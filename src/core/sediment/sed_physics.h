// sed_physics.h — pure __host__ __device__ sediment-physics inlines shared by the CUDA
// kernels and their CPU reference twins (CLAUDE.md hard rule: identical arithmetic on both
// sides, GPU-vs-CPU parity tested to rel. max-norm 1e-5). Every formula here is RESEARCH.md
// §5/§6 verbatim; every DERIVED constant (D*, (ρs−ρ)g, θ_cr, w_s, τ_cr) is COMPUTED from the
// config (ρ, ρs, ν, d) at call time — never hard-coded (RESEARCH §5 ⚠, research/02 §Verify).
//
// Units: SI (m, s, kg, Pa). Concentrations c are VOLUMETRIC (m³/m³); mass conc = ρs·c.
#pragma once

#ifdef __CUDACC__
#define SCOUR_HD __host__ __device__
#else
#define SCOUR_HD
#endif

#include <cmath>

namespace scour::core
{
	// Standard gravity [m/s²]. Not a "physics constant to invent": it is g, and RESEARCH's
	// verified derived numbers pin it — (ρs−ρ)g = 15922 N/m³ (ρ=1027) and (s−1)g = 15.50 m/s²
	// both reproduce with g = 9.81 (RESEARCH §1/§5). Kept as a named constant, not a magic literal.
	constexpr double SED_G = 9.81;

	// π to double precision (CUDA device code has no guaranteed M_PI). Used by Engelund–Fredsøe.
	constexpr double SED_PI = 3.14159265358979323846;

	// Critical packing fraction c_pack = 0.64 (RESEARCH §1/§7, FLOW-3D default = random close
	// packing). Bed porosity p = 1 − c_pack = 0.36. The Exner solids-fraction factor (1−p) = 0.64.
	constexpr double SED_CPACK = 0.64;

	// Submerged specific gravity s = ρs/ρ (RESEARCH §1; s = 2.58 at ρ=1027, ρs=2650).
	SCOUR_HD inline double sed_s(double rho, double rho_s) { return rho_s / rho; }

	// Submerged unit weight (ρs−ρ)·g [N/m³]  (RESEARCH §5: 15922 at ρ=1027).
	SCOUR_HD inline double sed_submerged_gamma(double rho, double rho_s) { return (rho_s - rho) * SED_G; }

	// (s−1)·g [m/s²]  (RESEARCH §1: 15.50 at ρ=1027) = submerged_gamma/ρ.
	SCOUR_HD inline double sed_s1g(double rho, double rho_s) { return (rho_s - rho) * SED_G / rho; }

	// Dimensionless grain size D* = d·[(s−1)g/ν²]^(1/3)  (RESEARCH §5; = 20313·d at ρ=1027, ν=1.36e-6).
	SCOUR_HD inline double grain_Dstar(double d, double rho, double rho_s, double nu)
	{
		double s1g = sed_s1g(rho, rho_s);
		return d * cbrt(s1g / (nu * nu));
	}

	// Soulsby (1997) settling velocity (RESEARCH §6.1): w_s = (ν/d)[√(10.36²+1.049 D*³) − 10.36].
	// 10.36² = 107.3296. Single formula for the whole range; matches Ferguson–Church within 10%.
	SCOUR_HD inline double settling_ws_dstar(double d, double nu, double Dstar)
	{
		double Ds3 = Dstar * Dstar * Dstar;
		return (nu / d) * (sqrt(107.3296 + 1.049 * Ds3) - 10.36);
	}
	SCOUR_HD inline double settling_ws(double d, double rho, double rho_s, double nu)
	{
		return settling_ws_dstar(d, nu, grain_Dstar(d, rho, rho_s, nu));
	}

	// Hindered-settling factor (RESEARCH §6.1): multiply w_s by (1−c)^4.7 above c≈0.001; clamp c≤0.35.
	// Below 0.001 the factor is >0.995 → return 1 exactly (RESEARCH §6.1 / research/02 §5).
	SCOUR_HD inline double hindered_factor(double c)
	{
		double cc = c;
		if (cc < 0.0) cc = 0.0;
		if (cc > 0.35) cc = 0.35; // clamp (dilute-suspension validity)
		if (cc <= 0.001) return 1.0;
		return pow(1.0 - cc, 4.7);
	}

	// Soulsby–Whitehouse (1997) critical Shields (RESEARCH §5):
	// θ_cr = 0.30/(1+1.2 D*) + 0.055·(1 − exp(−0.020 D*)). Valid all D* ≥ 0.1.
	SCOUR_HD inline double theta_cr_sw(double Dstar)
	{
		return 0.30 / (1.0 + 1.2 * Dstar) + 0.055 * (1.0 - exp(-0.020 * Dstar));
	}

	// Critical bed shear τ_cr = θ_cr·(ρs−ρ)g·d [Pa]  (RESEARCH §5).
	SCOUR_HD inline double tau_cr(double d, double rho, double rho_s, double nu)
	{
		double Ds = grain_Dstar(d, rho, rho_s, nu);
		return theta_cr_sw(Ds) * sed_submerged_gamma(rho, rho_s) * d;
	}

	// Shields parameter θ = τ/((ρs−ρ)g·d)  (RESEARCH §5).
	SCOUR_HD inline double shields(double tau, double d, double rho, double rho_s)
	{
		return tau / (sed_submerged_gamma(rho, rho_s) * d);
	}

	// van Rijn (2019) high-velocity damping (RESEARCH §6.2): f_D = 1 for θ′≤1, else 1/θ′.
	SCOUR_HD inline double f_damp(double theta_prime)
	{
		return theta_prime > 1.0 ? 1.0 / theta_prime : 1.0;
	}

	// van Rijn (2019) pickup E [kg/m²/s]  (RESEARCH §6.2):
	//   E = α·ρs·[(s−1)g·d50]^0.5·D*^0.3·f_D·T^1.5,  T = (τ′_b − τ_cr)/τ_cr,  α = 0.00033 (±30%).
	// θ_cr GATES erosion: T ≤ 0 ⇒ E = 0 (RESEARCH §5 asymmetry rule — deposition is NOT thresholded).
	// τ_prime is the GRAIN-related (skin-friction) bed shear stress [Pa] (RESEARCH §6.2, research/02 §5).
	SCOUR_HD inline double vanrijn_pickup(double tau_prime, double d, double rho, double rho_s, double nu, double alpha)
	{
		double Ds = grain_Dstar(d, rho, rho_s, nu);
		double gam = sed_submerged_gamma(rho, rho_s);
		double tcr = theta_cr_sw(Ds) * gam * d;
		double T = (tau_prime - tcr) / tcr;
		if (T <= 0.0) return 0.0; // clear-water: no pickup (θ_cr gate)
		double theta_p = tau_prime / (gam * d);       // grain Shields θ′
		double fD = f_damp(theta_p);                   // van Rijn 2019 damping
		double s1g_d = sed_s1g(rho, rho_s) * d;        // (s−1)g·d
		return alpha * rho_s * sqrt(s1g_d) * pow(Ds, 0.3) * fD * pow(T, 1.5);
	}

	// Deposition flux D = w_s·c_b [m/s, volumetric]  (RESEARCH §6.2 — NEVER divided by ρs;
	// c is volumetric, so D is already a volumetric flux). NEVER thresholded (RESEARCH §5).
	SCOUR_HD inline double deposition_flux(double ws, double c_b) { return ws * c_b; }

	// Net upward volumetric bed flux F_bed = E/ρs − w_s·c_b [m/s]  (RESEARCH §6.2, research/02 §11).
	// Applied to the bed-adjacent cell as a source dc/dt += F_bed/h (flux per cell height).
	SCOUR_HD inline double bed_net_flux(double E_kg, double rho_s, double ws, double c_b)
	{
		return E_kg / rho_s - ws * c_b;
	}

	// ==================== open-sea equilibrium inflow (RESEARCH §6; no new constant) ============
	// The far-field seabed upstream of the structure carries an EQUILIBRIUM suspended load: a flat
	// bed neither erodes nor deposits when this model's OWN van Rijn pickup balances its OWN settling
	// deposition, E(τ_eq)/ρs = w_s·c_b. So the near-bed equilibrium VOLUMETRIC concentration is
	//   c_b,eq = E(τ_eq)/(ρs·w_s).
	// This is NOT a new empirical constant (e.g. the van Rijn c_a 0.015 reference) — it is the fixed
	// point of the pickup+deposition pair already in this file. τ_eq is the far-field grain-skin bed
	// shear ρ·u*_eq². Clear-water (E=0, sub-threshold current) ⇒ 0: a current that can't mobilise sand
	// carries no suspended load, the correct degenerate case. ws = clear-water w_s(d50) (dilute limit).
	SCOUR_HD inline double equilibrium_cb(double tau_eq, double d, double rho, double rho_s, double nu, double alpha, double ws)
	{
		if (ws <= 0.0) return 0.0;
		double E = vanrijn_pickup(tau_eq, d, rho, rho_s, nu, alpha); // kg/m²/s, θ_cr-gated
		return E / (rho_s * ws);
	}

	// Rouse number R = w_s·σ_s/(κ·u*) (RESEARCH §6.2). ε_s = ν_t/σ_s with a parabolic ν_t = κ u* z(1−z/h)
	// ⇒ the steady settling/diffusion balance is the classical Rouse with exponent w_s/(β κ u*), β = 1/σ_s.
	// The V6 gate uses β=1 (σ_s=1) to match the classical analytic; production uses σ_s=0.7 (CLAUDE.md M4).
	SCOUR_HD inline double rouse_number(double ws, double sigma_s, double kappa, double ustar)
	{
		if (ustar <= 0.0) return 0.0;
		return ws * sigma_s / (kappa * ustar);
	}

	// Rouse (1937) equilibrium suspended-concentration profile (RESEARCH §6, research/02):
	//   c(z) = c_a·[ (a/z)·((h−z)/(h−a)) ]^R,   a = reference (near-bed) height, h = flow depth.
	// The SAME formula the V6 Rouse gate checks (sed_validation.cu). z clamped into (a,h) so the
	// endpoints stay finite: z≤a ⇒ c_a (near-bed anchor), z≥h ⇒ 0 (vanishes at the free surface).
	SCOUR_HD inline double rouse_profile(double z, double a, double h, double c_a, double R)
	{
		if (z <= a) return c_a;
		if (z >= h) return 0.0;
		double frac = (a / z) * ((h - z) / (h - a));
		if (frac <= 0.0) return 0.0;
		return c_a * pow(frac, R);
	}

	// ==================== M7 (M-DEP): trapping-efficiency engineering cross-checks =============
	// These closed forms are the deposition-validation "engineering layer" (research/15 §1,§6): the
	// simulator's resolved trapped volume is gated against them (±30%). No new empirical constant —
	// w_s and u* come from the config-derived physics above.

	// van Rijn (1987) trapping efficiency of a deepened area in cross-flow (research/15 Eq.1, from
	// 300 SUTRENCH runs): e_s = 1 − exp(−A_vr·L·d/h1²), A_vr = 0.25·(ws/u*1)·(1 + 2·ws/u*1). L =
	// effective settling length, d = deepening (h1−h0), h1 = flow depth in the deepened area, u*1 =
	// bed shear velocity there. e_s → 0 as d → 0. Stated accuracy ~25% (research/15 §1).
	SCOUR_HD inline double vanrijn_trap_efficiency(double ws, double ustar1, double L, double d, double h1)
	{
		if (ustar1 <= 0.0 || h1 <= 0.0 || d <= 0.0) return 0.0;
		double r = ws / ustar1;
		double A_vr = 0.25 * r * (1.0 + 2.0 * r);
		double e = 1.0 - exp(-A_vr * L * d / (h1 * h1));
		if (e < 0.0) e = 0.0; if (e > 1.0) e = 1.0;
		return e;
	}

	// Camp–Hazen settling-basin trap efficiency (research/07, research/15 §6 companion):
	// η = 1 − exp(−w_s·A/Q) = 1 − exp(−w_s·L/(h·U)) for a rectangular basin of plan length L, depth h,
	// throughflow U (A/Q = L/(h·U)). η > 0.8 ⇒ w_s·L/(h·U) > 1.6 (RESEARCH §9 design criterion).
	SCOUR_HD inline double camp_hazen_efficiency(double ws, double L, double h, double U)
	{
		if (h <= 0.0 || U <= 0.0) return 0.0;
		double e = 1.0 - exp(-ws * L / (h * U));
		if (e < 0.0) e = 0.0; if (e > 1.0) e = 1.0;
		return e;
	}

	// Thin-screen (sub-grid perforation) resistance coefficient k = 1/β² − 1 (Steiros–Hultmark 2018 /
	// Taylor–Davies; RESEARCH §8, research/13 §8), β = open-area ratio. The pressure jump across a
	// perforated face is Δp = ½·ρ·k·u_face² (u = superficial/face velocity, NOT pore velocity — using
	// pore velocity overpredicts by 1/β²). β→1 ⇒ k→0 (fully open); β→0 ⇒ k→∞ (solid).
	SCOUR_HD inline double screen_resistance_k(double beta_open)
	{
		if (beta_open <= 1e-6) return 1e12;      // ~solid
		if (beta_open >= 1.0) return 0.0;        // fully open
		return 1.0 / (beta_open * beta_open) - 1.0;
	}

	// ==================== M5: bed state + morphodynamics (RESEARCH §5, §7) =====================

	// Normalized fill F = f_pack/(c_pack·(1−phi_obs)) ∈ [0,1] (RESEARCH §7, research/13 §1).
	// Cell state: PRINTED (phi_obs≈1) / PACKED (F≥1) / INTERFACE (0<F<1) / FLUID (F=0).
	SCOUR_HD inline double normalized_fill(double f_pack, double phi_obs)
	{
		double denom = SED_CPACK * (1.0 - phi_obs);
		return denom > 1e-300 ? f_pack / denom : 1.0;
	}

	// Slope-corrected critical Shields (Soulsby 1997 eq. 80a, general 3D form; RESEARCH §5,
	// research/01 §6). β = bed slope angle (n_b to vertical), ψ = angle of near-bed flow to the
	// up-slope direction, φ = repose angle (rad). θ_cr(β,ψ) = θ_cr,flat·[cosψ·sinβ +
	// √(cos²β·tan²φ − sin²ψ·sin²β)] / tanφ. The √-argument is clamped ≥0 (downslope/overhang) and
	// the factor floored at 0.1 (overhang normals β>φ, RESEARCH §5 clamp — avalanching handles
	// those). Sanity (φ=32°, β=20°): 1.49 up / 0.39 down / 0.76 transverse.
	SCOUR_HD inline double theta_cr_slope(double theta_cr_flat, double beta, double psi, double phi_rep)
	{
		double tanphi = tan(phi_rep);
		double cb = cos(beta), sb = sin(beta), cp = cos(psi), sp = sin(psi);
		double arg = cb * cb * tanphi * tanphi - sp * sp * sb * sb;
		if (arg < 0.0) arg = 0.0;
		double factor = (cp * sb + sqrt(arg)) / tanphi;
		if (factor < 0.1) factor = 0.1; // overhang clamp (RESEARCH §5)
		return theta_cr_flat * factor;
	}

	// Winterwerp (1992) / FLOW-3D lifting velocity for erosion of a packed cell along n_b
	// (RESEARCH §7, research/13 §4): u_lift = coeff·D*^0.3·(θ − θ′_cr)^1.5·√((s−1)g·d) [m/s].
	// θ_cr GATES erosion: θ ≤ θ′_cr ⇒ u_lift = 0. The bed recedes at u_lift; grain volume flux to
	// suspension = u_lift·c_pack [m/s]. θ′_cr is the slope-corrected critical Shields.
	// `coeff` defaults to the RESEARCH §7 / Winterwerp value 0.018; it is exposed as a calibration
	// knob (MorphoParams.erosion_coeff) for the M6 supervised gate — the de-facto erosion-rate knob in
	// the bed-Exner path (the van Rijn α-pickup of §6.2 drives the M4 suspended path, not this Exner).
	SCOUR_HD inline double winterwerp_ulift(double theta, double theta_cr, double d, double rho, double rho_s, double nu, double coeff = 0.018)
	{
		if (theta <= theta_cr) return 0.0;
		double Ds = grain_Dstar(d, rho, rho_s, nu);
		double s1g_d = sed_s1g(rho, rho_s) * d;
		return coeff * pow(Ds, 0.3) * pow(theta - theta_cr, 1.5) * sqrt(s1g_d);
	}

	// Bedload dimensionless rate Φ = q_b/√((s−1)g·d³) (RESEARCH §6.3, research/03).
	// Default: Wong & Parker (2006) corrected MPM, Φ = 3.97·(θ′ − 0.0495)^1.5 (θ′>0.0495 else 0).
	// The 0.0495 is Wong–Parker's own reference (≈ θ_cr,SW of 0.2 mm sand); the on/off gate is the
	// slope-corrected Soulsby θ_cr applied by the caller (RESEARCH §5 threshold gates bedload).
	SCOUR_HD inline double bedload_phi_wp(double theta_prime)
	{
		if (theta_prime <= 0.0495) return 0.0;
		return 3.97 * pow(theta_prime - 0.0495, 1.5);
	}

	// Toggle: Engelund–Fredsøe (1976) probabilistic form (RESEARCH §6.3, research/03 §5 — the
	// formulation Roulund et al. 2005 used): Φ = 5p·(√θ − 0.7·√θ_cr),
	// p = [1 + ((π/6)·0.51/(θ − θ_cr))^4]^(−1/4). θ ≤ θ_cr (or √θ − 0.7√θ_cr ≤ 0) ⇒ Φ = 0.
	SCOUR_HD inline double bedload_phi_ef(double theta, double theta_cr)
	{
		if (theta <= theta_cr) return 0.0;
		double sq = sqrt(theta) - 0.7 * sqrt(theta_cr);
		if (sq <= 0.0) return 0.0;
		double base = (SED_PI / 6.0) * 0.51 / (theta - theta_cr);
		double b4 = base * base; b4 *= b4;
		double p = pow(1.0 + b4, -0.25);
		return 5.0 * p * sq;
	}

	// Dimensional bedload rate q_b = Φ·√((s−1)g·d³) [m²/s] (solids volume flux per width, RESEARCH §6.3).
	SCOUR_HD inline double bedload_qb(double phi, double d, double rho, double rho_s)
	{
		double s1g = sed_s1g(rho, rho_s);
		return phi * sqrt(s1g * d * d * d);
	}

	// Bedload layer thickness h_b = 0.3·d·D*^0.7·T^0.5, T = (θ−θ_cr)/θ_cr (RESEARCH §7, research/13 §5).
	SCOUR_HD inline double bedload_layer_thickness(double d, double Dstar, double T)
	{
		if (T <= 0.0) return 0.0;
		return 0.3 * d * pow(Dstar, 0.7) * sqrt(T);
	}
}
