// test_suspended.cu — M4 suspended-sediment: GPU-vs-CPU kernel parity (rel. max-norm 1e-5,
// CLAUDE.md hard rule) for effective-settling, scalar advection, diffusion, and bed exchange;
// plus sed_physics.h value checks against the RESEARCH §5/§6 / research/02 §4 tables.
#include "core/sediment/suspended.h"
#include "core/sediment/sed_physics.h"
#include "gpu_testutil.h"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

using namespace scour::core;

namespace
{
	MacGrid grid() { MacGrid g; g.nx = 6; g.ny = 8; g.nz = 10; g.h = 0.5; return g; }
	std::vector<double> uni(int n, unsigned seed, double lo, double hi)
	{
		std::mt19937 rng(seed); std::uniform_real_distribution<double> d(lo, hi);
		std::vector<double> v(n); for (auto& x : v) x = d(rng); return v;
	}
	BC vbc0() { BC b; b.xmin = b.xmax = b.ymin = b.ymax = b.zmin = b.zmax = WALL_FREESLIP; b.lid_u = 0.0; return b; }
}

// ---- (0) effective settling w_s,eff = w_s0·(1−c)^4.7 -------------------------------------
TEST(Suspended, EffectiveWsParity)
{
	MacGrid g = grid();
	auto hc = uni(g.p_count(), 11, 0.0, 0.34); // span the hindered regime up to the 0.35 clamp
	std::vector<double> ref; suspended_effective_ws_cpu(hc, ref, g, 0.03, 1);
	DevVec<double> c(hc), ws((size_t)g.p_count());
	suspended_effective_ws_gpu(c.p, ws.p, g, 0.03, 1);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(ws.download(), ref), 1e-5);
}

// ---- (1) MacCormack scalar advection with settling folded in ----------------------------
TEST(Suspended, AdvectParityNeumann)
{
	MacGrid g = grid();
	auto hc = uni(g.p_count(), 1, 0.0, 1.0);
	auto hu = uni(g.u_count(), 2, -0.1, 0.1);
	auto hv = uni(g.v_count(), 3, -0.1, 0.1);
	auto hw = uni(g.w_count(), 4, -0.1, 0.1);
	auto hws = uni(g.p_count(), 5, 0.01, 0.05);
	BC vbc = vbc0(); ScalarBC sbc; double dt = 0.05;
	for (int band : {0, 1})
	{
		std::vector<double> ref; suspended_advect_cpu(hc, hu, hv, hw, hws, ref, g, vbc, sbc, dt, band);
		DevVec<double> c(hc), u(hu), v(hv), w(hw), ws(hws), out((size_t)g.p_count()), scr((size_t)g.p_count());
		suspended_advect_gpu(c.p, u.p, v.p, w.p, ws.p, out.p, scr.p, g, vbc, sbc, dt, band);
		cudaDeviceSynchronize();
		EXPECT_LT(rel_maxnorm(out.download(), ref), 1e-5) << "band=" << band;
	}
}

TEST(Suspended, AdvectConservativeParity)
{
	MacGrid g = grid();
	auto hc = uni(g.p_count(), 41, 0.0, 1.0);
	auto hu = uni(g.u_count(), 42, -0.2, 0.2);
	auto hv = uni(g.v_count(), 43, -0.2, 0.2);
	auto hw = uni(g.w_count(), 44, -0.2, 0.2);
	auto hws = uni(g.p_count(), 45, 0.01, 0.05);
	double dt = 0.05;
	std::vector<double> ref; suspended_advect_cons_cpu(hc, hu, hv, hw, hws, ref, g, dt);
	DevVec<double> c(hc), u(hu), v(hv), w(hw), ws(hws), out((size_t)g.p_count()), scr((size_t)g.p_count());
	suspended_advect_cons_gpu(c.p, u.p, v.p, w.p, ws.p, out.p, scr.p, g, dt);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(out.download(), ref), 1e-5);
}

TEST(Suspended, AdvectConservativeConservesMass)
{
	// zero flux through all domain faces ⇒ Σc is invariant to machine precision.
	MacGrid g = grid();
	auto hc = uni(g.p_count(), 51, 0.1, 1.0);
	auto hu = uni(g.u_count(), 52, -0.2, 0.2);
	auto hv = uni(g.v_count(), 53, -0.2, 0.2);
	auto hw = uni(g.w_count(), 54, -0.2, 0.2);
	std::vector<double> hws(g.p_count(), 0.03);
	double dt = 0.05, m0 = 0; for (double x : hc) m0 += x;
	DevVec<double> c(hc), u(hu), v(hv), w(hw), ws(hws), out((size_t)g.p_count()), scr((size_t)g.p_count());
	std::vector<double> cur = hc;
	for (int s = 0; s < 20; ++s)
	{
		DevVec<double> cc(cur);
		suspended_advect_cons_gpu(cc.p, u.p, v.p, w.p, ws.p, out.p, scr.p, g, dt);
		cudaDeviceSynchronize(); cur = out.download();
	}
	double m1 = 0; for (double x : cur) m1 += x;
	EXPECT_LT(std::fabs(m1 - m0) / m0, 1e-12);
}

TEST(Suspended, AdvectParityDirichletTop)
{
	MacGrid g = grid();
	auto hc = uni(g.p_count(), 6, 0.0, 1.0);
	auto hu = uni(g.u_count(), 7, -0.05, 0.05);
	auto hv = uni(g.v_count(), 8, -0.05, 0.05);
	auto hw = uni(g.w_count(), 9, -0.05, 0.05);
	auto hws = uni(g.p_count(), 10, 0.02, 0.04);
	BC vbc = vbc0();
	ScalarBC sbc; sbc.zmax = SC_DIRICHLET; sbc.v_zmax = 0.0; sbc.zmin = SC_DIRICHLET; sbc.v_zmin = 0.123;
	double dt = 0.05;
	std::vector<double> ref; suspended_advect_cpu(hc, hu, hv, hw, hws, ref, g, vbc, sbc, dt, 0);
	DevVec<double> c(hc), u(hu), v(hv), w(hw), ws(hws), out((size_t)g.p_count()), scr((size_t)g.p_count());
	suspended_advect_gpu(c.p, u.p, v.p, w.p, ws.p, out.p, scr.p, g, vbc, sbc, dt, 0);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(out.download(), ref), 1e-5);
}

// ---- (2) variable-coefficient diffusion --------------------------------------------------
TEST(Suspended, DiffuseParity)
{
	MacGrid g = grid();
	auto hc = uni(g.p_count(), 21, 0.0, 1.0);
	auto hDc = uni(g.p_count(), 22, 0.0, 0.02);
	ScalarBC sbc; sbc.zmin = SC_DIRICHLET; sbc.v_zmin = 0.5; double dt = 0.5;
	std::vector<double> ref; suspended_diffuse_cpu(hc, ref, hDc, g, sbc, dt);
	DevVec<double> c(hc), Dc(hDc), out((size_t)g.p_count());
	suspended_diffuse_gpu(c.p, out.p, Dc.p, g, sbc, dt);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(out.download(), ref), 1e-5);
}

// ---- (3) flux-form bed exchange (pickup + deposition) -----------------------------------
TEST(Suspended, BedExchangeParity)
{
	MacGrid g = grid();
	auto hc = uni(g.p_count(), 31, 0.0, 0.1);
	auto htau = uni(g.p_count(), 32, 0.0, 2.0); // grain skin shear [Pa]: some below, some above τ_cr
	SedBedParams p; p.d50 = 0.2e-3; p.rho = 1027; p.rho_s = 2650; p.nu = 1.36e-6; p.alpha = 0.00033;
	p.ws0 = settling_ws(p.d50, p.rho, p.rho_s, p.nu); p.hindered = 1;
	double dt = 0.1;
	std::vector<double> cc = hc, bf; suspended_bed_exchange_cpu(cc, htau, bf, g, p, dt);
	DevVec<double> c(hc), tau(htau), bflux((size_t)(g.nx * g.ny));
	suspended_bed_exchange_gpu(c.p, tau.p, bflux.p, g, p, dt);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(c.download(), cc), 1e-5);
	EXPECT_LT(rel_maxnorm(bflux.download(), bf), 1e-5);
}

// ---- (4) open-x boundary advective flux (equilibrium inflow / outflow) -------------------
TEST(Suspended, BoundaryFluxParity)
{
	MacGrid g = grid();
	int nyz = g.ny * g.nz;
	double dt = 0.05;
	// seed 61: forward flow (inflow xmin / outflow xmax); seed 62: reversing flow (both branches).
	for (int seed : {61, 62})
	{
		double ulo = seed == 61 ? 0.0 : -0.2;
		auto hc = uni(g.p_count(), seed, 0.0, 1.0);
		auto hu = uni(g.u_count(), seed + 100, ulo, 0.2);
		auto gmin = uni(nyz, seed + 200, 0.0, 0.05);
		auto gmax = uni(nyz, seed + 300, 0.0, 0.05);
		std::vector<double> cc = hc, nf_ref;
		suspended_boundary_flux_cpu(cc, hu, gmin, gmax, &nf_ref, g, dt);
		DevVec<double> c(hc), u(hu), dmin(gmin), dmax(gmax), nf((size_t)(2 * nyz));
		suspended_boundary_flux_gpu(c.p, u.p, dmin.p, dmax.p, nf.p, g, dt);
		cudaDeviceSynchronize();
		EXPECT_LT(rel_maxnorm(c.download(), cc), 1e-5) << "seed=" << seed;
		EXPECT_LT(rel_maxnorm(nf.download(), nf_ref), 1e-5) << "seed=" << seed;
		// budget self-consistency: Σ net_flux [m³] == h³·Σ(Δc) over the whole field (exact accounting).
		double dsum = 0; for (int t = 0; t < g.p_count(); ++t) dsum += (cc[t] - hc[t]);
		double nfsum = 0; for (double x : nf_ref) nfsum += x;
		double vol = g.h * g.h * g.h;
		EXPECT_NEAR(nfsum, dsum * vol, 1e-11 * (std::fabs(dsum * vol) + 1e-30)) << "seed=" << seed;
	}
}

// ---- sed_physics.h value checks (RESEARCH §5/§6, research/02 §4) -------------------------
TEST(SedPhysics, SoulsbySettlingTable)
{
	// research/02 §4 table at 20 °C water (ρ=1025, ν=1.05e-6): w_s = 7.3/24.5/117.3/283 mm/s.
	double rho = 1025, rho_s = 2650, nu = 1.05e-6;
	EXPECT_NEAR(settling_ws(0.1e-3, rho, rho_s, nu), 0.00727, 0.00727 * 0.02);
	EXPECT_NEAR(settling_ws(0.2e-3, rho, rho_s, nu), 0.0245, 0.0245 * 0.02);
	EXPECT_NEAR(settling_ws(1.0e-3, rho, rho_s, nu), 0.1173, 0.1173 * 0.02);
	EXPECT_NEAR(settling_ws(5.0e-3, rho, rho_s, nu), 0.283, 0.283 * 0.03);
}

TEST(SedPhysics, DstarAndThetaCr)
{
	double rho = 1025, rho_s = 2650, nu = 1.05e-6;
	// D* = 24162·d(m) at these conditions (RESEARCH §5).
	EXPECT_NEAR(grain_Dstar(1e-3, rho, rho_s, nu), 24.162, 0.05);
	// Soulsby–Whitehouse table (research/11): θ_cr = 0.080/0.042/0.031/0.056 at d = 0.1/0.26/1/10 mm.
	EXPECT_NEAR(theta_cr_sw(grain_Dstar(0.1e-3, rho, rho_s, nu)), 0.080, 0.002);
	EXPECT_NEAR(theta_cr_sw(grain_Dstar(0.26e-3, rho, rho_s, nu)), 0.042, 0.002);
	EXPECT_NEAR(theta_cr_sw(grain_Dstar(1.0e-3, rho, rho_s, nu)), 0.031, 0.002);
	EXPECT_NEAR(theta_cr_sw(grain_Dstar(10.0e-3, rho, rho_s, nu)), 0.056, 0.002);
}

TEST(SedPhysics, HinderedSettling)
{
	// (1−c)^4.7: 0.995 / 0.79 / 0.61 at c = 0.001 / 0.05 / 0.1 (research/02 §5).
	EXPECT_NEAR(hindered_factor(0.001), 1.0, 1e-9);   // ≤0.001 ⇒ exactly 1 (RESEARCH §6.1)
	EXPECT_NEAR(hindered_factor(0.05), 0.79, 0.01);
	EXPECT_NEAR(hindered_factor(0.1), 0.61, 0.01);
	EXPECT_NEAR(hindered_factor(0.40), std::pow(1.0 - 0.35, 4.7), 1e-9); // clamp c≤0.35
}

TEST(SedPhysics, VanRijnPickupGateAndDamping)
{
	double rho = 1027, rho_s = 2650, nu = 1.36e-6, d = 0.2e-3, alpha = 0.00033;
	double tcr = tau_cr(d, rho, rho_s, nu);
	// θ_cr gates erosion: below τ_cr ⇒ E = 0 (RESEARCH §5 asymmetry rule).
	EXPECT_EQ(vanrijn_pickup(0.5 * tcr, d, rho, rho_s, nu, alpha), 0.0);
	EXPECT_EQ(vanrijn_pickup(tcr, d, rho, rho_s, nu, alpha), 0.0);
	// Above τ_cr ⇒ E > 0 and increases with stress.
	double E1 = vanrijn_pickup(2.0 * tcr, d, rho, rho_s, nu, alpha);
	double E2 = vanrijn_pickup(4.0 * tcr, d, rho, rho_s, nu, alpha);
	EXPECT_GT(E1, 0.0);
	EXPECT_GT(E2, E1);
	// f_D damping: f_D = 1 for θ′≤1, 1/θ′ for θ′>1 (van Rijn 2019).
	EXPECT_DOUBLE_EQ(f_damp(0.5), 1.0);
	EXPECT_DOUBLE_EQ(f_damp(1.0), 1.0);
	EXPECT_DOUBLE_EQ(f_damp(2.0), 0.5);
}
