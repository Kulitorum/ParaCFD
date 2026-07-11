// test_campaign.cpp — M8 campaign machinery unit tests: run-length/timescale (RESEARCH §7),
// MORFAC limiter, equilibrium fits (exp + Welzel) and CI-based shape distinguishability. Host-only
// (run_protocol.h is pure inlines); part of the fast `unit` suite.
#include "core/campaign/run_protocol.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using namespace scour::core;

TEST(Campaign, ScourTimescaleRoulund)
{
	// research/16 E6 / research/11 §8 worked example: θ≈0.12, δ/D=4 ⇒ T*≈0.21, and for the Roulund
	// grains (D=0.1, d50=0.26 mm, ρ=1027) T ≈ 130 s.
	double Ts = scour_Tstar(0.12, 4.0, 0);
	EXPECT_NEAR(Ts, 0.212, 0.01);
	double T = scour_T(Ts, 0.1, 0.26e-3, 1027.0, 2650.0);
	EXPECT_GT(T, 110.0); EXPECT_LT(T, 150.0); // ≈ 128 s
	// Alternatives are flatter at low θ (Whitehouse/Larsen-Fuhrman overestimate less time at low θ).
	EXPECT_LT(scour_Tstar(0.03, 4.0, 2), scour_Tstar(0.03, 4.0, 0)); // θ^-1.5 < θ^-2.2 at θ<1
}

TEST(Campaign, MorfacLimiter)
{
	// A fast bed rate forces M→1; a tiny rate allows the cap.
	EXPECT_NEAR(morfac_limit(1.0, 10, 0.01, 0.05, 0.05, 10.0), 1.0, 1e-9);   // huge rate ⇒ clamp low
	EXPECT_NEAR(morfac_limit(1e-9, 10, 0.01, 0.05, 0.05, 10.0), 10.0, 1e-9); // ~no motion ⇒ cap
	// Intermediate: frac·h/(rate·N·dt) = 0.05·0.05/(1e-4·10·0.01) = 0.0025/1e-5 = 250 → clamped to cap 10.
	EXPECT_NEAR(morfac_limit(1e-4, 10, 0.01, 0.05, 0.05, 10.0), 10.0, 1e-9);
}

TEST(Campaign, EquilibriumFitsRecoverSyntheticCurve)
{
	// Each fitter recovers its OWN functional form to high accuracy.
	// Exponential data S(t) = 2.0·(1 − e^(−t/50)) sampled to 100 s (2T).
	std::vector<double> t, Se;
	for (int i = 1; i <= 50; ++i) { double tt = 2.0 * i; t.push_back(tt); Se.push_back(2.0 * (1.0 - std::exp(-tt / 50.0))); }
	EqFit fe = fit_exponential(t, Se);
	EXPECT_NEAR(fe.S_eq, 2.0, 0.05);
	EXPECT_NEAR(fe.param, 50.0, 5.0); // recovered T
	EXPECT_GT(fe.r2, 0.999);
	// Welzel/hyperbolic data S(t) = 1.5·(1 − 1/(1+0.05·t)) — fit_welzel recovers a=1.5, b=0.05.
	std::vector<double> Sw;
	for (double tt : t) Sw.push_back(1.5 * (1.0 - 1.0 / (1.0 + 0.05 * tt)));
	EqFit fw = fit_welzel(t, Sw);
	EXPECT_NEAR(fw.S_eq, 1.5, 0.05);
	EXPECT_NEAR(fw.param, 0.05, 0.01);
	EXPECT_GT(fw.r2, 0.999);
	// Cross-form disagreement is REAL (research/16): the hyperbolic has a fatter tail, so fitting it
	// to exponential data extrapolates a higher asymptote — the protocol extends the run when >15%.
	EqFit fw_on_exp = fit_welzel(t, Se);
	EXPECT_GT(fit_disagreement(fe, fw_on_exp), 0.0);
}

TEST(Campaign, RankingDistinguishability)
{
	// Two shapes with clearly separated equilibria + small CIs ⇒ distinguishable; overlapping ⇒ not.
	EqFit a_e{1.0, 50, 0, 0.99}, a_w{1.05, 0.02, 0, 0.99};
	EqFit b_e{2.0, 60, 0, 0.99}, b_w{2.05, 0.02, 0, 0.99};
	ShapeRank ra = equilibrium_with_ci(a_e, a_w, 0.02);
	ShapeRank rb = equilibrium_with_ci(b_e, b_w, 0.02);
	EXPECT_TRUE(distinguishable(ra, rb, 1.0)); // ~1.0 vs ~2.0, CIs ~0.05 ⇒ clearly separable
	// Nearly-equal shapes with wide CIs ⇒ NOT distinguishable.
	EqFit c_e{1.0, 50, 0, 0.9}, c_w{1.4, 0.02, 0, 0.9}; // 40% fit disagreement ⇒ big CI
	ShapeRank rc = equilibrium_with_ci(c_e, c_w, 0.05);
	EXPECT_FALSE(distinguishable(ra, rc, 1.0));
}

TEST(Campaign, BrierSkillScoreAndRetention)
{
	// Perfect model ⇒ BSS = 1; no-change model (= initial) ⇒ BSS = 0; anti-skill ⇒ < 0.
	std::vector<double> meas{1.0, 2.0, 3.0}, init{0.0, 0.0, 0.0};
	EXPECT_NEAR(brier_skill_score(meas, meas, init), 1.0, 1e-12);
	EXPECT_NEAR(brier_skill_score(init, meas, init), 0.0, 1e-12);
	std::vector<double> worse{-1.0, -2.0, -3.0}; // reversed sign ⇒ worse than no-change
	EXPECT_LT(brier_skill_score(worse, meas, init), 0.0);
	// A "reasonable" model (halves the error) ⇒ BSS = 0.75.
	std::vector<double> half{0.5, 1.0, 1.5};
	EXPECT_NEAR(brier_skill_score(half, meas, init), 0.75, 1e-12);
	// Retention.
	EXPECT_NEAR(retention_fraction(0.8, 1.0), 0.8, 1e-12);
	EXPECT_DOUBLE_EQ(retention_fraction(0.5, 0.0), 0.0);
}
