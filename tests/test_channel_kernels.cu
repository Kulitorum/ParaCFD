// test_channel_kernels.cu — GPU-vs-CPU parity for the M2 open-channel masked kernels
// (rel. max-norm tol 1e-5, CLAUDE.md hard rule) + a ChannelMgpcg convergence check.
#include "core/fluid/channel_ops.h"
#include "core/fluid/channel_pressure.h"
#include "core/fluid/mac_ops.h"
#include "gpu_testutil.h"

#include <gtest/gtest.h>

#include <random>
#include <vector>

using namespace scour::core;

namespace
{
	MacGrid small_grid() { MacGrid g; g.nx = 8; g.ny = 10; g.nz = 6; g.h = 0.05; return g; }
	ChannelBC chan_bc()
	{
		ChannelBC bc; bc.inlet_mode = INLET_UNIFORM; bc.U_inlet = 1.3; bc.Uc = 1.3;
		bc.solid_mode = SOLID_NOSLIP; return bc;
	}
	std::vector<double> randf(int n, unsigned seed, double amp = 1.0)
	{
		std::mt19937 rng(seed); std::uniform_real_distribution<double> d(-amp, amp);
		std::vector<double> v(n); for (auto& x : v) x = d(rng); return v;
	}
	// a solid block in the interior + its 1-cell dilation
	void make_masks(MacGrid g, std::vector<unsigned char>& solid, std::vector<unsigned char>& near)
	{
		solid.assign(g.p_count(), 0);
		for (int k = 1; k < g.nz - 1; ++k) for (int j = 3; j < 6; ++j) for (int i = 3; i < 6; ++i) solid[g.pidx(i, j, k)] = 1;
		near.assign(g.p_count(), 0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			bool nn = false;
			for (int dk = -1; dk <= 1 && !nn; ++dk) for (int dj = -1; dj <= 1 && !nn; ++dj) for (int di = -1; di <= 1 && !nn; ++di)
			{
				int ii = i + di, jj = j + dj, kk = k + dk;
				if (ii < 0 || ii >= g.nx || jj < 0 || jj >= g.ny || kk < 0 || kk >= g.nz) continue;
				if (solid[g.pidx(ii, jj, kk)]) nn = true;
			}
			near[g.pidx(i, j, k)] = nn ? 1 : 0;
		}
	}
	const double TOL = 1e-5;
}

TEST(ChannelKernels, Advection)
{
	MacGrid g = small_grid(); ChannelBC bc = chan_bc();
	std::vector<unsigned char> solid, near; make_masks(g, solid, near);
	auto u = randf(g.u_count(), 1), v = randf(g.v_count(), 2), w = randf(g.w_count(), 3);
	double dt = 0.02;
	std::vector<double> ru, rv, rw;
	ch_advect_cpu(u, v, w, ru, rv, rw, solid, near, g, bc, dt);
	DevVec<double> du(u), dv(v), dw(w), ou(g.u_count()), ov(g.v_count()), ow(g.w_count());
	DevVec<double> scratch((size_t)std::max({g.u_count(), g.v_count(), g.w_count()}));
	DevVec<unsigned char> ds(solid), dn(near);
	ch_advect_gpu(du.p, dv.p, dw.p, ou.p, ov.p, ow.p, scratch.p, ds.p, dn.p, g, bc, dt);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(ou.download(), ru), TOL);
	EXPECT_LT(rel_maxnorm(ov.download(), rv), TOL);
	EXPECT_LT(rel_maxnorm(ow.download(), rw), TOL);
}

TEST(ChannelKernels, Smagorinsky)
{
	MacGrid g = small_grid(); ChannelBC bc = chan_bc();
	std::vector<unsigned char> solid, near; make_masks(g, solid, near);
	auto u = randf(g.u_count(), 4), v = randf(g.v_count(), 5), w = randf(g.w_count(), 6);
	std::vector<double> rn;
	ch_smagorinsky_cpu(u, v, w, rn, solid, g, bc, 0.12);
	DevVec<double> du(u), dv(v), dw(w), dn(g.p_count()); DevVec<unsigned char> ds(solid);
	ch_smagorinsky_gpu(du.p, dv.p, dw.p, dn.p, ds.p, g, bc, 0.12);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dn.download(), rn), TOL);
}

TEST(ChannelKernels, Diffusion)
{
	MacGrid g = small_grid(); ChannelBC bc = chan_bc();
	std::vector<unsigned char> solid, near; make_masks(g, solid, near);
	auto u = randf(g.u_count(), 7), v = randf(g.v_count(), 8), w = randf(g.w_count(), 9);
	auto nut = randf(g.p_count(), 10, 0.5); for (auto& x : nut) x = std::fabs(x);
	double dt = 1e-3, nu = 1e-2;
	std::vector<double> ru, rv, rw;
	ch_diffuse_cpu(u, v, w, ru, rv, rw, &nut, solid, g, bc, dt, nu);
	DevVec<double> du(u), dv(v), dw(w), dnut(nut), ou(g.u_count()), ov(g.v_count()), ow(g.w_count());
	DevVec<unsigned char> ds(solid);
	ch_diffuse_gpu(du.p, dv.p, dw.p, ou.p, ov.p, ow.p, dnut.p, ds.p, g, bc, dt, nu);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(ou.download(), ru), TOL);
	EXPECT_LT(rel_maxnorm(ov.download(), rv), TOL);
	EXPECT_LT(rel_maxnorm(ow.download(), rw), TOL);
}

TEST(ChannelKernels, PoissonRhsApplyJacobi)
{
	MacGrid g = small_grid(); ChannelBC bc = chan_bc();
	std::vector<unsigned char> solid, near; make_masks(g, solid, near);
	auto u = randf(g.u_count(), 11), v = randf(g.v_count(), 12), w = randf(g.w_count(), 13);
	auto p = randf(g.p_count(), 14);
	for (int c = 0; c < g.p_count(); ++c) if (solid[c]) p[c] = 0.0;
	double rho = 1027.0, dt = 5e-3; int dir = 1;
	std::vector<double> r_rhs, r_ap;
	ch_poisson_rhs_cpu(u, v, w, r_rhs, solid, g, rho, dt);
	ch_poisson_apply_cpu(p, r_ap, solid, g, dir);
	DevVec<double> du(u), dv(v), dw(w), dp(p), drhs(g.p_count()), dap(g.p_count());
	DevVec<unsigned char> ds(solid);
	ch_poisson_rhs_gpu(du.p, dv.p, dw.p, drhs.p, ds.p, g, rho, dt);
	ch_poisson_apply_gpu(dp.p, dap.p, ds.p, g, dir);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(drhs.download(), r_rhs), TOL);
	EXPECT_LT(rel_maxnorm(dap.download(), r_ap), TOL);

	auto pc = p; ch_jacobi_cpu(pc, r_rhs, solid, g, dir, 2.0 / 3.0, 3);
	DevVec<double> dp2(p), drhs2(r_rhs), dscr(g.p_count());
	ch_jacobi_gpu(dp2.p, drhs2.p, dscr.p, ds.p, g, dir, 2.0 / 3.0, 3);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dp2.download(), pc), TOL);
}

TEST(ChannelKernels, SubtractGradientAndDiv)
{
	MacGrid g = small_grid(); ChannelBC bc = chan_bc();
	std::vector<unsigned char> solid, near; make_masks(g, solid, near);
	auto u = randf(g.u_count(), 19), v = randf(g.v_count(), 20), w = randf(g.w_count(), 21);
	auto p = randf(g.p_count(), 22);
	double rho = 1027.0, dt = 5e-3; int dir = 1;
	auto uc = u, vc = v, wc = w;
	ch_subtract_gradient_cpu(uc, vc, wc, p, solid, g, bc, rho, dt, dir);
	DevVec<double> du(u), dv(v), dw(w), dp(p); DevVec<unsigned char> ds(solid);
	ch_subtract_gradient_gpu(du.p, dv.p, dw.p, dp.p, ds.p, g, bc, rho, dt, dir);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(du.download(), uc), TOL);
	EXPECT_LT(rel_maxnorm(dv.download(), vc), TOL);
	EXPECT_LT(rel_maxnorm(dw.download(), wc), TOL);

	double cpu = ch_max_div_cpu(uc, vc, wc, solid, g);
	DevVec<double> scr(g.p_count());
	double gpu = ch_max_div_gpu(du.p, dv.p, dw.p, scr.p, ds.p, g);
	cudaDeviceSynchronize();
	EXPECT_NEAR(gpu, cpu, 1e-5 * std::max(1.0, cpu));
}

TEST(ChannelMgpcgSolver, ConvergesWithMaskAndDirichlet)
{
	MacGrid g; g.nx = 32; g.ny = 32; g.nz = 16; g.h = 0.05;
	std::vector<unsigned char> solid(g.p_count(), 0);
	// a solid cylinder-ish block
	for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
	{
		double x = i - 12.0, y = j - 16.0; if (x * x + y * y < 25.0) solid[g.pidx(i, j, k)] = 1;
	}
	auto b = randf(g.p_count(), 99);
	for (int c = 0; c < g.p_count(); ++c) if (solid[c]) b[c] = 0.0;
	DevVec<double> db(b), dx(g.p_count()), dr(g.p_count());
	DevVec<unsigned char> ds(solid);
	ChannelMgpcg solver(g);
	solver.build_masks(ds.p);
	SolveResult r = solver.solve(dx.p, db.p, 1e-8, 200, /*warm=*/false);
	// verify residual r = b - A x is small relative to ||b||
	ch_poisson_residual_gpu(dx.p, db.p, dr.p, ds.p, g, 1);
	cudaDeviceSynchronize();
	double rn = std::sqrt(dot_gpu(dr.p, dr.p, g.p_count()));
	double bn = std::sqrt(dot_gpu(db.p, db.p, g.p_count()));
	EXPECT_TRUE(r.converged);
	EXPECT_LT(rn / bn, 1e-6);
	EXPECT_LT(r.iters, 60);
}
