// run_protocol.h — M8 campaign machinery core (RESEARCH §7, research/16): the run-length /
// timescale / MORFAC protocol and the equilibrium-fit + shape-ranking maths, as pure host inlines
// (no CUDA, no Qt) so they are unit-testable and reusable by the ranking gate and the batch runner.
//
// The scour/fill timescale spans a factor ~2700 across the design matrix, so fixed-duration runs
// produce ranking artifacts (research/16 E6). The protocol: derive T per scenario, run ≥ 0.5·T of
// morphological time, fit BOTH S(t)=S_eq(1−e^(−t/T)) and Welzel a(1−1/(1+bt)), extend if they
// disagree > 15%, and rank shapes on fitted equilibria at matched t/T with confidence intervals.
//
// Constants come from RESEARCH.md verbatim; every derived number is computed from the config
// (ρ, ρs, ν) via sed_physics.h — never hard-coded.
#pragma once

#include "core/sediment/sed_physics.h"

#include <cmath>
#include <vector>

namespace scour::core
{
	// ---- Scour timescale (Sumer & Fredsøe 2002; RESEARCH §7, research/05 E5/E6, /16 E5) ----------
	// Nondimensional T* for steady current: T* = (1/2000)·(δ/D)·θ^(−2.2), δ = boundary-layer/flow
	// depth. Alternatives kept for the uncertainty band (research/16 E5): Whitehouse θ^(−1.29),
	// Larsen–Fuhrman θ^(−3/2). `exponent` selects −2.2 (0, default), −1.29 (1), −1.5 (2).
	inline double scour_Tstar(double theta, double delta_over_D, int exponent = 0)
	{
		if (theta <= 0.0) return 0.0;
		if (exponent == 1) return 0.014 * std::pow(theta, -1.29);      // Whitehouse 1998
		if (exponent == 2) return (1.0 / 2000.0) * delta_over_D * std::pow(theta, -1.5); // Larsen–Fuhrman
		return (1.0 / 2000.0) * delta_over_D * std::pow(theta, -2.2);  // Sumer–Fredsøe (canonical)
	}
	// Dimensional T = T*·D²/√(g(s−1)·d50³) [s] (RESEARCH §7; research/05 E5).
	inline double scour_T(double Tstar, double D, double d50, double rho, double rho_s)
	{
		double s1g = sed_s1g(rho, rho_s);
		double denom = std::sqrt(s1g * d50 * d50 * d50);
		return denom > 0 ? Tstar * D * D / denom : 0.0;
	}
	// Convenience: full T from flow quantities (u* → θ → T*, δ = flow depth).
	inline double scour_T_from_flow(double ustar, double D, double d50, double rho, double rho_s, double delta, int exponent = 0)
	{
		double tau = rho * ustar * ustar;
		double theta = shields(tau, d50, rho, rho_s);
		double Ts = scour_Tstar(theta, D > 0 ? delta / D : 0.0, exponent);
		return scour_T(Ts, D, d50, rho, rho_s);
	}

	// ---- MORFAC / bed-change limiter (RESEARCH §7, research/16 E3/E4) ----------------------------
	// Largest MORFAC that keeps the per-morphology-step bed change within the limiter:
	// M·N·dt·max_rate ≤ frac·h (frac = 0.05 default). Returns M clamped to [1, Mcap] (Mcap = 10
	// steady current / 5 reversing). max_rate = peak |dz_b/dt| [m/s] (unaccelerated).
	inline double morfac_limit(double max_rate, double N, double dt, double h, double frac = 0.05, double Mcap = 10.0)
	{
		if (max_rate <= 0.0 || N * dt <= 0.0) return Mcap;
		double M = frac * h / (max_rate * N * dt);
		if (M < 1.0) M = 1.0; if (M > Mcap) M = Mcap;
		return M;
	}

	// ---- Equilibrium fits (research/16 E5/E7) ----------------------------------------------------
	struct EqFit { double S_eq = 0, param = 0, sse = 1e300, r2 = 0; };

	// S(t) = S_eq·(1 − e^(−t/T)); grid-search T, linear S_eq per T. `param` = T.
	inline EqFit fit_exponential(const std::vector<double>& t, const std::vector<double>& S)
	{
		EqFit best; size_t n = t.size(); if (n < 4 || S.size() != n) return best;
		double tmax = t.back(), Smean = 0; for (double s : S) Smean += s; Smean /= n;
		double sstot = 0; for (double s : S) sstot += (s - Smean) * (s - Smean);
		for (int it = 0; it < 500; ++it)
		{
			double T = 0.02 * tmax * std::pow(2000.0, it / 499.0);
			double num = 0, den = 0;
			for (size_t i = 0; i < n; ++i) { double b = 1.0 - std::exp(-t[i] / T); num += S[i] * b; den += b * b; }
			if (den <= 0) continue;
			double Seq = num / den, sse = 0;
			for (size_t i = 0; i < n; ++i) { double f = Seq * (1.0 - std::exp(-t[i] / T)); sse += (S[i] - f) * (S[i] - f); }
			if (sse < best.sse) { best.sse = sse; best.S_eq = Seq; best.param = T; best.r2 = sstot > 0 ? 1.0 - sse / sstot : 0.0; }
		}
		return best;
	}
	// Welzel a·(1 − 1/(1+b·t)); grid-search b, linear a. `param` = b, `S_eq` = a.
	inline EqFit fit_welzel(const std::vector<double>& t, const std::vector<double>& S)
	{
		EqFit best; size_t n = t.size(); if (n < 4 || S.size() != n) return best;
		double tmax = t.back(), Smean = 0; for (double s : S) Smean += s; Smean /= n;
		double sstot = 0; for (double s : S) sstot += (s - Smean) * (s - Smean);
		for (int it = 0; it < 500; ++it)
		{
			double b = (0.02 / tmax) * std::pow(2000.0, it / 499.0);
			double num = 0, den = 0;
			for (size_t i = 0; i < n; ++i) { double gg = 1.0 - 1.0 / (1.0 + b * t[i]); num += S[i] * gg; den += gg * gg; }
			if (den <= 0) continue;
			double a = num / den, sse = 0;
			for (size_t i = 0; i < n; ++i) { double f = a * (1.0 - 1.0 / (1.0 + b * t[i])); sse += (S[i] - f) * (S[i] - f); }
			if (sse < best.sse) { best.sse = sse; best.S_eq = a; best.param = b; best.r2 = sstot > 0 ? 1.0 - sse / sstot : 0.0; }
		}
		return best;
	}

	// Fit-disagreement fraction (research/16: extend the run if > 0.15).
	inline double fit_disagreement(const EqFit& a, const EqFit& b)
	{
		double m = std::fabs(a.S_eq) > 1e-30 ? a.S_eq : 1e-30;
		return std::fabs(a.S_eq - b.S_eq) / std::fabs(m);
	}

	// ---- Ranking with confidence intervals (research/16 ranking protocol) ------------------------
	// Two shapes are distinguishable only if their fitted-equilibrium CIs do not overlap. A crude but
	// honest CI from the two independent fits: half-width = max(|exp−welzel|, rmse-based floor). This
	// captures the research/16 "10–20% at 0.5T" extrapolation spread without overstating precision.
	struct ShapeRank { double value = 0, ci = 0; };
	inline ShapeRank equilibrium_with_ci(const EqFit& expfit, const EqFit& welfit, double sample_rmse)
	{
		ShapeRank r;
		r.value = 0.5 * (expfit.S_eq + welfit.S_eq);
		double spread = std::fabs(expfit.S_eq - welfit.S_eq);
		r.ci = std::fmax(spread, sample_rmse); // 1-σ-ish half width
		return r;
	}
	// Distinguishable iff the CIs (k-σ) do not overlap.
	inline bool distinguishable(const ShapeRank& a, const ShapeRank& b, double ksigma = 1.0)
	{
		return std::fabs(a.value - b.value) > ksigma * (a.ci + b.ci);
	}

	// ---- Morphodynamic skill + retention (research/15 §5, PLAN M7/M8/M9 gates) -------------------
	// Brier Skill Score of a modelled bed vs measured, relative to the "no-change" (initial) baseline:
	//   BSS = 1 − mean((z_model − z_meas)²) / mean((z_initial − z_meas)²).
	// Classification (Sutherland 2004): >0.8 excellent, 0.6–0.8 good, 0.3–0.6 reasonable, 0–0.3 poor,
	// <0 bad (worse than predicting no change). Arrays must be equal length.
	inline double brier_skill_score(const std::vector<double>& model, const std::vector<double>& measured, const std::vector<double>& initial)
	{
		size_t n = measured.size();
		if (n == 0 || model.size() != n || initial.size() != n) return -1e30;
		double num = 0, den = 0;
		for (size_t i = 0; i < n; ++i) { double e = model[i] - measured[i]; num += e * e; double d = initial[i] - measured[i]; den += d * d; }
		if (den <= 0) return num <= 0 ? 1.0 : -1e30; // no initial change to beat
		return 1.0 - num / den;
	}

	// Retention fraction (PLAN M8 retention protocol): trapped volume after a higher-flow stage vs at
	// the fitted equilibrium of the design stage. 1 ⇒ fully retained; <1 ⇒ partial washout.
	inline double retention_fraction(double V_after, double V_eq) { return V_eq > 0 ? V_after / V_eq : 0.0; }
}
