// test_mgpcg.cu — MGPCG solver: correctness on a manufactured Poisson problem,
// grid-independent iteration count (the RESEARCH §3.3 headline: ~1 order per 2 iters),
// and the plain-Jacobi fallback path. Singular all-Neumann box.
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h"
#include "core/fluid/mgpcg.h"
#include "gpu_testutil.h"

#include <gtest/gtest.h>

#include <random>
#include <vector>

using namespace scour::core;

namespace
{
	std::vector<double> zero_mean(std::vector<double> v)
	{
		double s = 0.0; for (double x : v) s += x; s /= v.size();
		for (double& x : v) x -= s; return v;
	}
	std::vector<double> randf(int n, unsigned seed)
	{
		std::mt19937 rng(seed); std::uniform_real_distribution<double> d(-1.0, 1.0);
		std::vector<double> v(n); for (auto& x : v) x = d(rng); return v;
	}
}

// A x = b with b = A * p_exact must recover p_exact (both mean-zero).
TEST(Mgpcg, RecoversManufacturedSolution)
{
	MacGrid g; g.nx = 32; g.ny = 32; g.nz = 32; g.h = 1.0 / 32;
	int n = g.p_count();
	auto pe = zero_mean(randf(n, 101));
	DevVec<double> dpe(pe), db(n), dx(n);
	poisson_apply_gpu(dpe.p, db.p, g); // b = A p_exact (exactly mean-zero)
	cudaDeviceSynchronize();

	MgpcgSolver solver(g);
	SolveResult r = solver.solve(dx.p, db.p, 1e-9, 60, /*warm_start=*/false, /*use_mg=*/true);
	EXPECT_TRUE(r.converged) << "relres=" << r.relres << " iters=" << r.iters;
	EXPECT_LE(r.relres, 1e-9);
	EXPECT_LE(r.iters, 25);
	EXPECT_LT(rel_maxnorm(dx.download(), pe), 1e-6);
}

// Gate criterion: MGPCG reaches ||r||/||b|| <= 1e-4 in <= 14 iterations at 128^3.
TEST(Mgpcg, GateIterationBudget128)
{
	MacGrid g; g.nx = 128; g.ny = 128; g.nz = 128; g.h = 1.0 / 128;
	int n = g.p_count();
	auto pe = zero_mean(randf(n, 202));
	DevVec<double> dpe(pe), db(n), dx(n);
	poisson_apply_gpu(dpe.p, db.p, g);
	cudaDeviceSynchronize();
	MgpcgSolver solver(g);
	SolveResult r = solver.solve(dx.p, db.p, 1e-4, 30, false, true);
	std::printf("[mgpcg 128^3] iters=%d relres=%.3e\n", r.iters, r.relres);
	EXPECT_TRUE(r.converged);
	EXPECT_LE(r.iters, 14);
}

// Grid independence: iteration count to 1e-6 must not blow up with N.
TEST(Mgpcg, GridIndependentIterations)
{
	int prev = 0;
	for (int N : {32, 64})
	{
		MacGrid g; g.nx = N; g.ny = N; g.nz = N; g.h = 1.0 / N;
		int n = g.p_count();
		auto pe = zero_mean(randf(n, 303 + N));
		DevVec<double> dpe(pe), db(n), dx(n);
		poisson_apply_gpu(dpe.p, db.p, g);
		cudaDeviceSynchronize();
		MgpcgSolver solver(g);
		SolveResult r = solver.solve(dx.p, db.p, 1e-6, 40, false, true);
		std::printf("[mgpcg %d^3] iters=%d relres=%.3e\n", N, r.iters, r.relres);
		EXPECT_TRUE(r.converged);
		EXPECT_LE(r.iters, 20);
		if (prev) EXPECT_LE(r.iters, prev + 4); // near grid-independent
		prev = r.iters;
	}
}

// The bring-up/debug Jacobi fallback must at least reduce the residual.
TEST(Mgpcg, JacobiFallbackReducesResidual)
{
	MacGrid g; g.nx = 16; g.ny = 16; g.nz = 16; g.h = 1.0 / 16;
	int n = g.p_count();
	auto pe = zero_mean(randf(n, 404));
	DevVec<double> dpe(pe), db(n), dx(n);
	poisson_apply_gpu(dpe.p, db.p, g);
	cudaDeviceSynchronize();
	MgpcgSolver solver(g);
	SolveResult r = solver.solve(dx.p, db.p, 1e-3, 4000, false, /*use_mg=*/false);
	std::printf("[jacobi 16^3] sweeps=%d relres=%.3e\n", r.iters, r.relres);
	EXPECT_LT(r.relres, 1.0); // converging (Jacobi is slow by design)
}
