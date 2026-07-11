// test_fluid_kernels.cu — GPU-vs-CPU parity for every M1 MAC-solver kernel
// (rel. max-norm tol 1e-5, CLAUDE.md hard rule). Fields are double, so the kernels
// and their serial CPU twins should agree to ~1e-12; the test guards launch/index/
// memory correctness on a small non-cubic grid.
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h"
#include "gpu_testutil.h"

#include <gtest/gtest.h>

#include <random>
#include <vector>

using namespace scour::core;

namespace
{
	MacGrid small_grid()
	{
		MacGrid g; g.nx = 6; g.ny = 8; g.nz = 10; g.h = 0.05; return g;
	}
	BC cavity_bc() { BC bc; bc.lid_u = 1.3; return bc; } // exercises the lid ghost

	std::vector<double> randf(int n, unsigned seed, double amp = 1.0)
	{
		std::mt19937 rng(seed);
		std::uniform_real_distribution<double> d(-amp, amp);
		std::vector<double> v(n);
		for (auto& x : v) x = d(rng);
		return v;
	}
	const double TOL = 1e-5;
}

TEST(FluidKernels, MacCormackAdvection)
{
	MacGrid g = small_grid(); BC bc = cavity_bc();
	auto u = randf(g.u_count(), 1), v = randf(g.v_count(), 2), w = randf(g.w_count(), 3);
	double dt = 0.02;

	std::vector<double> ru, rv, rw;
	mac_advect_cpu(u, v, w, ru, rv, rw, g, bc, dt, 1);

	DevVec<double> du(u), dv(v), dw(w);
	DevVec<double> ou(g.u_count()), ov(g.v_count()), ow(g.w_count());
	DevVec<double> scratch((size_t)std::max({g.u_count(), g.v_count(), g.w_count()}));
	mac_advect_gpu(du.p, dv.p, dw.p, ou.p, ov.p, ow.p, scratch.p, g, bc, dt, 1);
	cudaDeviceSynchronize();

	EXPECT_LT(rel_maxnorm(ou.download(), ru), TOL);
	EXPECT_LT(rel_maxnorm(ov.download(), rv), TOL);
	EXPECT_LT(rel_maxnorm(ow.download(), rw), TOL);
}

TEST(FluidKernels, SmagorinskyNut)
{
	MacGrid g = small_grid(); BC bc = cavity_bc();
	auto u = randf(g.u_count(), 4), v = randf(g.v_count(), 5), w = randf(g.w_count(), 6);
	std::vector<double> rn;
	smagorinsky_nut_cpu(u, v, w, rn, g, bc, 0.12);
	DevVec<double> du(u), dv(v), dw(w), dn(g.p_count());
	smagorinsky_nut_gpu(du.p, dv.p, dw.p, dn.p, g, bc, 0.12);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dn.download(), rn), TOL);
}

TEST(FluidKernels, ExplicitDiffusion)
{
	MacGrid g = small_grid(); BC bc = cavity_bc();
	auto u = randf(g.u_count(), 7), v = randf(g.v_count(), 8), w = randf(g.w_count(), 9);
	auto nut = randf(g.p_count(), 10, 0.5);
	for (auto& x : nut) x = std::fabs(x); // nu_t >= 0
	double dt = 1e-3, nu = 1e-2;
	std::vector<double> ru, rv, rw;
	mac_diffuse_cpu(u, v, w, ru, rv, rw, &nut, g, bc, dt, nu);
	DevVec<double> du(u), dv(v), dw(w), dn(nut), ou(g.u_count()), ov(g.v_count()), ow(g.w_count());
	mac_diffuse_gpu(du.p, dv.p, dw.p, ou.p, ov.p, ow.p, dn.p, g, bc, dt, nu);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(ou.download(), ru), TOL);
	EXPECT_LT(rel_maxnorm(ov.download(), rv), TOL);
	EXPECT_LT(rel_maxnorm(ow.download(), rw), TOL);
}

TEST(FluidKernels, PoissonRhsAndApply)
{
	MacGrid g = small_grid();
	auto u = randf(g.u_count(), 11), v = randf(g.v_count(), 12), w = randf(g.w_count(), 13);
	auto p = randf(g.p_count(), 14);
	double rho = 1027.0, dt = 5e-3;
	std::vector<double> r_rhs, r_ap;
	poisson_rhs_cpu(u, v, w, r_rhs, g, rho, dt);
	poisson_apply_cpu(p, r_ap, g);
	DevVec<double> du(u), dv(v), dw(w), dp(p), drhs(g.p_count()), dap(g.p_count());
	poisson_rhs_gpu(du.p, dv.p, dw.p, drhs.p, g, rho, dt);
	poisson_apply_gpu(dp.p, dap.p, g);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(drhs.download(), r_rhs), TOL);
	EXPECT_LT(rel_maxnorm(dap.download(), r_ap), TOL);
}

TEST(FluidKernels, JacobiAndGsBand)
{
	MacGrid g = small_grid();
	auto p0 = randf(g.p_count(), 15);
	auto rhs = randf(g.p_count(), 16);
	// Jacobi (3 sweeps)
	{
		auto pc = p0;
		jacobi_smooth_cpu(pc, rhs, g, 2.0 / 3.0, 3);
		DevVec<double> dp(p0), drhs(rhs), dscr(g.p_count());
		jacobi_smooth_gpu(dp.p, drhs.p, dscr.p, g, 2.0 / 3.0, 3);
		cudaDeviceSynchronize();
		EXPECT_LT(rel_maxnorm(dp.download(), pc), TOL);
	}
	// Gauss-Seidel band (2 sweeps, band=2, forward)
	{
		auto pc = p0;
		gs_band_cpu(pc, rhs, g, 2, 2, true);
		DevVec<double> dp(p0), drhs(rhs);
		gs_band_gpu(dp.p, drhs.p, g, 2, 2, true);
		cudaDeviceSynchronize();
		// GS colour order differs (parallel red-black vs serial), but a red-black GS
		// sweep is order-independent within a colour -> identical result.
		EXPECT_LT(rel_maxnorm(dp.download(), pc), TOL);
	}
}

TEST(FluidKernels, RestrictProlong)
{
	MacGrid gf; gf.nx = 8; gf.ny = 8; gf.nz = 12; gf.h = 0.05;
	MacGrid gc; gc.nx = 4; gc.ny = 4; gc.nz = 6; gc.h = 0.10;
	auto fine = randf(gf.p_count(), 17);
	auto coarse0 = randf(gc.p_count(), 18);
	// restrict
	{
		std::vector<double> rc;
		restrict_cpu(fine, rc, gf, gc);
		DevVec<double> df(fine), dc(gc.p_count());
		restrict_gpu(df.p, dc.p, gf, gc);
		cudaDeviceSynchronize();
		EXPECT_LT(rel_maxnorm(dc.download(), rc), TOL);
	}
	// prolong-add
	{
		auto fc = fine;
		prolong_add_cpu(coarse0, fc, gc, gf);
		DevVec<double> dco(coarse0), dfi(fine);
		prolong_add_gpu(dco.p, dfi.p, gc, gf);
		cudaDeviceSynchronize();
		EXPECT_LT(rel_maxnorm(dfi.download(), fc), TOL);
	}
}

TEST(FluidKernels, SubtractGradientAndDivergence)
{
	MacGrid g = small_grid();
	auto u = randf(g.u_count(), 19), v = randf(g.v_count(), 20), w = randf(g.w_count(), 21);
	auto p = randf(g.p_count(), 22);
	double rho = 1027.0, dt = 5e-3;
	// subtract gradient
	auto uc = u, vc = v, wc = w;
	subtract_gradient_cpu(uc, vc, wc, p, g, rho, dt);
	DevVec<double> du(u), dv(v), dw(w), dp(p);
	subtract_gradient_gpu(du.p, dv.p, dw.p, dp.p, g, rho, dt);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(du.download(), uc), TOL);
	EXPECT_LT(rel_maxnorm(dv.download(), vc), TOL);
	EXPECT_LT(rel_maxnorm(dw.download(), wc), TOL);

	// max|div| (scalar)
	double cpu = max_abs_divergence_cpu(uc, vc, wc, g);
	DevVec<double> scr(g.p_count());
	double gpu = max_abs_divergence_gpu(du.p, dv.p, dw.p, scr.p, g);
	cudaDeviceSynchronize();
	EXPECT_NEAR(gpu, cpu, 1e-5 * std::max(1.0, cpu));
}
