// test_seabed.cu — GPU-vs-CPU parity (rel. max-norm 1e-5, CLAUDE.md hard rule) for the three
// erodible-seabed bridge kernels in core/sediment/seabed_morpho.*: the elevated/moving-bed grain-
// skin τ_b probe, the suspended-field fluid confinement, and the bed clamp/pin. The physics is
// shared with the CPU twin through the bedshear.h inlines, so these certify launch/index/memory.
#include "core/sediment/seabed_morpho.h"
#include "core/sediment/sed_physics.h" // SED_CPACK
#include "gpu_testutil.h"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

using namespace scour::core;

namespace
{
	MacGrid grid() { MacGrid g; g.nx = 6; g.ny = 5; g.nz = 24; g.h = 0.05; return g; }
	std::vector<double> uni(int n, unsigned seed, double lo, double hi)
	{
		std::mt19937 rng(seed); std::uniform_real_distribution<double> d(lo, hi);
		std::vector<double> v(n); for (auto& x : v) x = d(rng); return v;
	}
}

TEST(Seabed, BedshearParity)
{
	MacGrid g = grid();
	double cpack = SED_CPACK;
	auto u = uni(g.u_count(), 11, -0.6, 0.9);
	auto v = uni(g.v_count(), 12, -0.4, 0.4);
	// per-column bed elevations spread through the reservoir (grain thickness G = z_b·cpack).
	auto zb = uni(g.nx * g.ny, 13, 0.3, 0.8);
	std::vector<double> G(g.nx * g.ny);
	for (size_t i = 0; i < G.size(); ++i) G[i] = zb[i] * cpack;

	WallParams wp; wp.d50 = 0.2e-3; wp.nu = 1.36e-6; wp.rho = 1027.0; wp.regime = WALL_AUTO;

	std::vector<double> tx_c, ty_c, us_c;
	seabed_bedshear_cpu(u, v, G, tx_c, ty_c, us_c, g, wp, cpack);

	DevVec<double> du(u), dv(v), dG(G), dtx(g.nx * g.ny), dty(g.nx * g.ny), dus(g.nx * g.ny);
	seabed_bedshear_gpu(du.p, dv.p, dG.p, dtx.p, dty.p, dus.p, g, wp, cpack);
	cudaDeviceSynchronize();

	EXPECT_LT(rel_maxnorm(dtx.download(), tx_c), 1e-5);
	EXPECT_LT(rel_maxnorm(dty.download(), ty_c), 1e-5);
	EXPECT_LT(rel_maxnorm(dus.download(), us_c), 1e-5);
	// sanity: some column produced a non-zero u* (flow above the bed drives shear).
	double mx = 0; for (double x : us_c) mx = std::max(mx, std::fabs(x));
	EXPECT_GT(mx, 0.0);
}

TEST(Seabed, ConfineConservesAndClearsBelowBed)
{
	MacGrid g = grid();
	double cpack = SED_CPACK;
	int np = g.p_count(), ncol = g.nx * g.ny;
	auto c = uni(np, 21, 0.0, 0.01);
	auto zb = uni(ncol, 22, 0.3, 0.7);
	std::vector<double> G(ncol);
	for (int i = 0; i < ncol; ++i) G[i] = zb[i] * cpack;
	// a small structure block sitting on the bed in one column region.
	std::vector<unsigned char> structure(np, 0);
	for (int k = 16; k < 19; ++k) structure[g.pidx(2, 2, k)] = 1;

	double total0 = 0; for (double x : c) total0 += x;

	std::vector<double> c_cpu = c;
	seabed_confine_c_cpu(c_cpu, G, &structure, g, cpack);

	DevVec<double> dc(c), dG(G); DevVec<unsigned char> ds(structure);
	seabed_confine_c_gpu(dc.p, dG.p, ds.p, g, cpack);
	cudaDeviceSynchronize();
	auto c_gpu = dc.download();

	EXPECT_LT(rel_maxnorm(c_gpu, c_cpu), 1e-5);

	// mass conserved by the lift, and no c left below the first fluid cell.
	double total1 = 0; for (double x : c_cpu) total1 += x;
	EXPECT_NEAR(total1, total0, 1e-12 * (total0 + 1e-30));
	for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
	{
		int kbed = (int)std::floor((G[j * g.nx + i] / cpack) / g.h);
		for (int k = 0; k < kbed; ++k) EXPECT_DOUBLE_EQ(c_cpu[g.pidx(i, j, k)], 0.0);
	}
}

TEST(Seabed, BedPostClampAndPin)
{
	MacGrid g = grid();
	int ncol = g.nx * g.ny;
	auto G = uni(ncol, 31, -0.1, 2.0); // deliberately out of range to exercise the clamp
	auto Gfix = uni(ncol, 32, 0.2, 0.5);
	std::vector<unsigned char> frozen(ncol, 0);
	for (int i = 0; i < ncol; i += 3) frozen[i] = 1;
	double Gmax = (g.nz * g.h - 2.0 * g.h) * SED_CPACK;

	std::vector<double> G_cpu = G;
	seabed_bed_post_cpu(G_cpu, Gfix, &frozen, g, Gmax);

	DevVec<double> dG(G), dGf(Gfix); DevVec<unsigned char> df(frozen);
	seabed_bed_post_gpu(dG.p, dGf.p, df.p, g, Gmax);
	cudaDeviceSynchronize();

	EXPECT_LT(rel_maxnorm(dG.download(), G_cpu), 1e-5);
	for (int i = 0; i < ncol; ++i)
	{
		if (frozen[i]) EXPECT_DOUBLE_EQ(G_cpu[i], Gfix[i]);
		else { EXPECT_GE(G_cpu[i], 0.0); EXPECT_LE(G_cpu[i], Gmax); }
	}
}

TEST(Seabed, ShearMultiplierParity)
{
	// M6 HSV shear multiplier: τ → m·τ, u* → √m·u* (per column). GPU-vs-CPU parity + the scaling law.
	MacGrid g = grid();
	int ncol = g.nx * g.ny;
	auto taux = uni(ncol, 41, -1.0, 1.0), tauy = uni(ncol, 42, -1.0, 1.0), us = uni(ncol, 43, 0.0, 0.1);
	auto mult = uni(ncol, 44, 0.5, 4.0); // per-column amplification (incl. >1)
	std::vector<double> tx = taux, ty = tauy, u = us;
	seabed_apply_shear_mult_cpu(tx, ty, u, mult, g);
	DevVec<double> dtx(taux), dty(tauy), du(us), dm(mult);
	seabed_apply_shear_mult_gpu(dtx.p, dty.p, du.p, dm.p, g);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dtx.download(), tx), 1e-5);
	EXPECT_LT(rel_maxnorm(dty.download(), ty), 1e-5);
	EXPECT_LT(rel_maxnorm(du.download(), u), 1e-5);
	EXPECT_NEAR(tx[3], taux[3] * mult[3], 1e-12);       // τ scaled by m
	EXPECT_NEAR(u[3], us[3] * std::sqrt(mult[3]), 1e-12); // u* scaled by √m
}
