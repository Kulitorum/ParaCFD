// test_bedshear.cu — GPU-vs-CPU parity for the M3 bed-shear wall-function kernels and the
// periodic-channel closure kernels (rel. max-norm 1e-5, CLAUDE.md), plus direct checks of
// the log-law inversion, the Christoffersen–Jonsson fixed point, the cavity four-branch
// rule, and the anti-stair-step smoothing mitigation. RESEARCH §4, research/04.
#include "core/fluid/bedshear_ops.h"
#include "core/fluid/periodic_ops.h"
#include "gpu_testutil.h"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

using namespace scour::core;

namespace
{
	MacGrid grid() { MacGrid g; g.nx = 8; g.ny = 10; g.nz = 12; g.h = 0.05; return g; }
	std::vector<double> randf(int n, unsigned seed, double amp = 1.0)
	{
		std::mt19937 rng(seed); std::uniform_real_distribution<double> d(-amp, amp);
		std::vector<double> v(n); for (auto& x : v) x = d(rng); return v;
	}
	const double TOL = 1e-5;
	WallParams wp_default()
	{
		WallParams wp; wp.kappa = 0.40; wp.d50 = 0.2e-3; wp.nu = 1.36e-6; wp.rho = 1027.0;
		wp.regime = WALL_AUTO; wp.cj_iters = 5; wp.t_avg = 2.0; wp.smooth = 1; return wp;
	}
}

TEST(BedShear, RoughInversionAnalytic)
{
	// u* = κ·U_p/ln(z_p/z0). Hand value.
	double kappa = 0.40, Up = 0.6, zp = 0.075, z0 = 0.2e-3 / 12.0;
	double us = wall_ustar_rough(Up, zp, z0, kappa);
	double expect = kappa * Up / std::log(zp / z0);
	EXPECT_NEAR(us, expect, 1e-12);
	// Forced-rough regime path agrees.
	double us2 = wall_ustar(Up, zp, 0.2e-3, 1.36e-6, kappa, WALL_ROUGH, 5);
	EXPECT_NEAR(us2, expect, 1e-12);
}

TEST(BedShear, TransitionalFixedPointConsistent)
{
	// d50 = 0.2 mm at U_p = 0.6 m/s, zp = 0.075 m is transitional (ks+ ≈ 12): the returned
	// u* must satisfy the Christoffersen–Jonsson fixed point u* = κU_p/ln(zp/z0_cj(u*)).
	double kappa = 0.40, Up = 0.6, zp = 0.075, d50 = 0.2e-3, nu = 1.36e-6, ks = 2.5 * d50;
	double us = wall_ustar(Up, zp, d50, nu, kappa, WALL_TRANSITIONAL, 12);
	double ksplus = us * ks / nu;
	EXPECT_GT(ksplus, 5.0); EXPECT_LT(ksplus, 70.0); // genuinely transitional
	double z0 = wall_z0_cj(us, ks, nu);
	double residual = std::fabs(us - kappa * Up / std::log(zp / z0)) / us;
	EXPECT_LT(residual, 1e-4);
	// Transitional z0 is smaller than the rough z0 here ⇒ slightly lower u* (research/04).
	double us_rough = wall_ustar_rough(Up, zp, ks / 30.0, kappa);
	EXPECT_LT(us, us_rough);
}

TEST(BedShear, CavityFourBranch)
{
	double h = 0.05, ks = 2.5 * 0.2e-3, zp0 = wall_zp(h, ks), nu = 1.36e-6;
	EXPECT_EQ(cavity_branch(10.0 * zp0, h, zp0, 0.5, nu), CAV_STANDARD);
	// QUARTER band 6h ≤ Hc < 4·zp0 is only non-empty when zp0 > 1.5h (coarse grid / large
	// grains, so z_p is set by 2ks): use zp0=0.075 with a finer h=0.01 → band [0.06,0.30).
	EXPECT_EQ(cavity_branch(0.15, 0.01, 0.075, 0.5, nu), CAV_QUARTER);
	// 3h ≤ Hc < 6h with high vs low gap Reynolds number.
	double Hc = 4.0 * h;
	EXPECT_EQ(cavity_branch(Hc, h, zp0, 1.0, nu), CAV_GAP_TURB);       // Re_gap = 1·0.2/nu ≫ 2800
	EXPECT_EQ(cavity_branch(Hc, h, zp0, 1e-3, nu), CAV_GAP_LAMINAR);   // Re_gap = 0.2 ≪ 2800
	EXPECT_EQ(cavity_branch(2.0 * h, h, zp0, 0.5, nu), CAV_FROZEN);    // Hc < 3h
	// Gap closures positive; laminar τ = 6μU/Hc.
	double tl = cavity_tau_gap_laminar(0.1, Hc, 1027.0, nu);
	EXPECT_NEAR(tl, 6.0 * 1027.0 * nu * 0.1 / Hc, 1e-15);
	EXPECT_GT(cavity_tau_gap_turb(1.0, Hc, ks / 30.0, 1027.0, 0.40), 0.0);
}

TEST(BedShear, ComputeParity)
{
	MacGrid g = grid(); WallParams wp = wp_default();
	auto u = randf(g.u_count(), 11, 1.0), v = randf(g.v_count(), 12, 0.3);
	int np = g.nx * g.ny;
	std::vector<double> rtx, rty, rus;
	bedshear_compute_cpu(u, v, rtx, rty, rus, g, wp);
	DevVec<double> du(u), dv(v), dtx(np), dty(np), dus(np);
	bedshear_compute_gpu(du.p, dv.p, dtx.p, dty.p, dus.p, g, wp);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dtx.download(), rtx), TOL);
	EXPECT_LT(rel_maxnorm(dty.download(), rty), TOL);
	EXPECT_LT(rel_maxnorm(dus.download(), rus), TOL);
}

TEST(BedShear, SmoothFilterSinkParity)
{
	MacGrid g = grid(); WallParams wp = wp_default();
	int np = g.nx * g.ny;
	auto tx = randf(np, 21, 2.0), ty = randf(np, 22, 1.0);
	// smooth
	auto ctx = tx, cty = ty; bedshear_smooth_cpu(ctx, cty, g, 2);
	DevVec<double> gtx(tx), gty(ty), scr(np);
	bedshear_smooth_gpu(gtx.p, gty.p, scr.p, g, 2); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(gtx.download(), ctx), TOL);
	EXPECT_LT(rel_maxnorm(gty.download(), cty), TOL);
	// filter (EMA + rate-limit), starting from a previous EMA state
	auto ex = randf(np, 23, 1.0), ey = randf(np, 24, 1.0);
	auto cex = ex, cey = ey; bedshear_filter_cpu(tx, ty, cex, cey, g, wp, 0.02, 0);
	DevVec<double> gex(ex), gey(ey), dtx2(tx), dty2(ty);
	bedshear_filter_gpu(dtx2.p, dty2.p, gex.p, gey.p, g, wp, 0.02, 0); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(gex.download(), cex), TOL);
	EXPECT_LT(rel_maxnorm(gey.download(), cey), TOL);
	// sink
	auto u = randf(g.u_count(), 25), vv = randf(g.v_count(), 26);
	auto cu = u, cv = vv; bedshear_apply_sink_cpu(cu, cv, cex, cey, g, 1027.0, 0.02);
	DevVec<double> gu(u), gv(vv), sex(cex), sey(cey);
	bedshear_apply_sink_gpu(gu.p, gv.p, sex.p, sey.p, g, 1027.0, 0.02); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(gu.download(), cu), TOL);
	EXPECT_LT(rel_maxnorm(gv.download(), cv), TOL);
}

TEST(BedShear, SmoothingReducesBanding)
{
	// Synthetic grid-pitch banding (alternating columns) is damped by the 3×3 smoother
	// (anti-stair-step mitigation, research/04 §mitigations).
	MacGrid g = grid(); int np = g.nx * g.ny;
	std::vector<double> tx(np), ty(np, 0.0);
	for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) tx[j * g.nx + i] = 1.0 + ((i % 2) ? 0.5 : -0.5);
	auto band = [&](const std::vector<double>& f) {
		double mn = 1e9, mx = -1e9, s = 0; for (double v : f) { mn = std::min(mn, v); mx = std::max(mx, v); s += v; }
		return (mx - mn) / (s / f.size()); };
	double before = band(tx);
	bedshear_smooth_cpu(tx, ty, g, 2);
	double after = band(tx);
	EXPECT_LT(after, before * 0.6);
}

TEST(PeriodicOps, MixlenParity)
{
	MacGrid g = grid();
	auto u = randf(g.u_count(), 31, 1.0), v = randf(g.v_count(), 32, 0.2);
	std::vector<double> rn; periodic_mixlen_cpu(u, v, rn, g, 0.40, g.nz * g.h);
	DevVec<double> du(u), dv(v), dn(g.p_count());
	periodic_mixlen_gpu(du.p, dv.p, dn.p, g, 0.40, g.nz * g.h); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dn.download(), rn), TOL);
}

TEST(PeriodicOps, ZDiffuseParity)
{
	MacGrid g = grid();
	auto u = randf(g.u_count(), 41, 1.0), v = randf(g.v_count(), 42, 0.5);
	auto nut = randf(g.p_count(), 43, 0.01);
	for (auto& x : nut) x = std::fabs(x); // ν_t ≥ 0
	std::vector<double> ru, rv; periodic_zdiffuse_cpu(u, v, ru, rv, nut, g, 0.01, 1.36e-6);
	DevVec<double> du(u), dv(v), dn(nut), ou(g.u_count()), ov(g.v_count());
	periodic_zdiffuse_gpu(du.p, dv.p, ou.p, ov.p, dn.p, g, 0.01, 1.36e-6); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(ou.download(), ru), TOL);
	EXPECT_LT(rel_maxnorm(ov.download(), rv), TOL);
}
