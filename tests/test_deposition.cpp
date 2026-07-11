// test_deposition.cpp — M7 (M-DEP) engineering cross-check layer: value/limit checks for the
// van Rijn (1987) trapping efficiency, Camp–Hazen basin efficiency, and the thin-screen resistance
// k = 1/β²−1 (RESEARCH §8, research/15 §1/§6, research/13 §8). Host-only (sed_physics.h inlines
// compile as plain C++); part of the fast `unit` suite.
#include "core/sediment/sed_physics.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace scour::core;

TEST(Deposition, VanRijnTrapEfficiencyLimitsAndValue)
{
	// d → 0 ⇒ e_s → 0 (no deepening, nothing trapped).
	EXPECT_DOUBLE_EQ(vanrijn_trap_efficiency(0.02, 0.02, 1.0, 0.0, 0.4), 0.0);
	// Hand value: ws/u*1 = 1 ⇒ A_vr = 0.25·1·(1+2) = 0.75; e_s = 1−exp(−0.75·1·0.1/0.16) = 0.3739…
	double e = vanrijn_trap_efficiency(0.02, 0.02, 1.0, 0.1, 0.4);
	EXPECT_NEAR(e, 1.0 - std::exp(-0.75 * 0.1 / 0.16), 1e-12);
	EXPECT_NEAR(e, 0.37422, 1e-4); // 1 − exp(−0.46875)
	// Monotone increasing in ws/u*1 (heavier / slower flow traps more) and in d.
	EXPECT_GT(vanrijn_trap_efficiency(0.04, 0.02, 1.0, 0.1, 0.4), e);
	EXPECT_GT(vanrijn_trap_efficiency(0.02, 0.02, 1.0, 0.2, 0.4), e);
	// Bounded in [0,1].
	double big = vanrijn_trap_efficiency(0.1, 0.01, 100.0, 5.0, 0.4);
	EXPECT_LE(big, 1.0); EXPECT_GE(big, 0.0);
}

TEST(Deposition, CampHazenBasinEfficiency)
{
	// η = 1 − exp(−ws·L/(h·U)). The RESEARCH §9 design criterion: η > 0.8 needs ws·L/(h·U) > 1.6.
	double eta = camp_hazen_efficiency(0.02, 1.6, 0.4, 0.05); // ws·L/(hU) = 0.02·1.6/(0.4·0.05) = 1.6
	EXPECT_NEAR(eta, 1.0 - std::exp(-1.6), 1e-12);
	EXPECT_GT(eta, 0.79); EXPECT_LT(eta, 0.81); // ≈ 0.798
	// Faster throughflow ⇒ lower efficiency; zero settling ⇒ zero.
	EXPECT_LT(camp_hazen_efficiency(0.02, 1.6, 0.4, 0.5), eta);
	EXPECT_DOUBLE_EQ(camp_hazen_efficiency(0.0, 1.6, 0.4, 0.05), 0.0);
}

TEST(Deposition, ScreenResistanceK)
{
	EXPECT_DOUBLE_EQ(screen_resistance_k(1.0), 0.0);               // fully open ⇒ no jump
	EXPECT_NEAR(screen_resistance_k(0.30), 1.0 / 0.09 - 1.0, 1e-9); // β=0.30 ⇒ k≈10.11
	EXPECT_NEAR(screen_resistance_k(0.25), 15.0, 1e-9);            // β=0.25 ⇒ k=15
	// Optimum porosity band β≈0.25–0.35 (RESEARCH §9): k in ~[7.2, 15].
	EXPECT_GT(screen_resistance_k(0.35), 7.0);
	EXPECT_LT(screen_resistance_k(0.35), 8.0);
	EXPECT_GT(screen_resistance_k(0.05), 100.0); // near-solid ⇒ large resistance
}
