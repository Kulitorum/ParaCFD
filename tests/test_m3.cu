// test_m3.cu — SEM plane GPU-vs-CPU parity + Cholesky sanity (RESEARCH §8, research/14 §3)
// and precursor record/replay round-trip (research/14 §4).
#include "core/fluid/precursor.h"
#include "core/fluid/sem_inlet.h"
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
	std::vector<double> signs(int n, unsigned seed)
	{
		std::mt19937 rng(seed); std::uniform_int_distribution<int> d(0, 1);
		std::vector<double> v(n); for (auto& x : v) x = d(rng) ? 1.0 : -1.0; return v;
	}
}

TEST(Sem, CholeskyReproducesStresses)
{
	double us = 0.0345, h = 5.0, gamma = 1.4, z = 0.7;
	double a11, a22, a31, a33; sem_cholesky(us, z, h, gamma, a11, a22, a31, a33);
	double R11 = nezu_urms(us, z, h, gamma); R11 *= R11;
	double R22 = nezu_vrms(us, z, h, gamma); R22 *= R22;
	double R33 = nezu_wrms(us, z, h, gamma); R33 *= R33;
	double R13 = nezu_uw(us, z, h, gamma);
	EXPECT_NEAR(a11 * a11, R11, 1e-12);
	EXPECT_NEAR(a22 * a22, R22, 1e-12);
	EXPECT_NEAR(a11 * a31, R13, 1e-12);            // <u'w'>
	EXPECT_NEAR(a31 * a31 + a33 * a33, R33, 1e-12); // <w'w'>
	EXPECT_LT(R13, 0.0);                            // −u'w' = u*²(1−z/h) > 0
}

TEST(Sem, PlaneParity)
{
	MacGrid g = grid();
	SemParams p; p.sigma = 1.0; p.N = 40; p.ustar = 0.0345; p.h_dom = g.nz * g.h;
	p.U_d = 1.0; p.gamma = 1.5; p.Ly = g.ny * g.h; p.Lz = g.nz * g.h;
	// eddy state inside the box straddling the inlet plane.
	auto ex = uni(p.N, 1, -p.sigma, p.sigma);
	auto ey = uni(p.N, 2, -p.sigma, p.Ly + p.sigma);
	auto ez = uni(p.N, 3, -p.sigma, p.Lz + p.sigma);
	auto e1 = signs(p.N, 4), e2 = signs(p.N, 5), e3 = signs(p.N, 6);
	std::vector<double> ru, rv, rw;
	sem_planes_cpu(ex, ey, ez, e1, e2, e3, ru, rv, rw, g, p);
	DevVec<double> dex(ex), dey(ey), dez(ez), d1(e1), d2(e2), d3(e3);
	DevVec<double> du((size_t)g.ny * g.nz), dv((size_t)(g.ny + 1) * g.nz), dw((size_t)g.ny * (g.nz + 1));
	sem_planes_gpu(dex.p, dey.p, dez.p, d1.p, d2.p, d3.p, du.p, dv.p, dw.p, g, p);
	cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(du.download(), ru), 1e-5);
	EXPECT_LT(rel_maxnorm(dv.download(), rv), 1e-5);
	EXPECT_LT(rel_maxnorm(dw.download(), rw), 1e-5);
	// SEM injects non-zero fluctuations.
	double s = 0; for (double x : ru) s += std::fabs(x); EXPECT_GT(s, 0.0);
}

TEST(Precursor, RecordAndInterpolate)
{
	MacGrid g; g.nx = 4; g.ny = 3; g.nz = 4; g.h = 0.5;
	PrecursorLibrary lib(g, 10.0); // 10 Hz → record_dt = 0.1 s
	size_t nu = (size_t)g.ny * g.nz, nv = (size_t)(g.ny + 1) * g.nz, nw = (size_t)g.ny * (g.nz + 1);
	const int M = 5;
	for (int m = 0; m < M; ++m)
	{
		std::vector<double> up(nu, (double)m), vp(nv, 2.0 * m), wp(nw, -1.0 * m);
		bool kept = lib.maybe_record(m * 0.1, up, vp, wp);
		EXPECT_TRUE(kept);
	}
	// a duplicate-time record before the next slot is rejected.
	std::vector<double> up(nu, 99.0), vp(nv, 99.0), wp(nw, 99.0);
	EXPECT_FALSE(lib.maybe_record(0.35, up, vp, wp));
	EXPECT_EQ(lib.size(), M);

	std::vector<double> su, sv, sw;
	lib.sample(0.0, su, sv, sw); EXPECT_NEAR(su[0], 0.0, 1e-12);
	lib.sample(0.15, su, sv, sw); EXPECT_NEAR(su[0], 1.5, 1e-12);  // between frame 1 and 2
	lib.sample(0.1, su, sv, sw); EXPECT_NEAR(sv[0], 2.0, 1e-12);   // exact frame 1 → 2·1
	// seamless loop: period = M·0.1 = 0.5 s ⇒ t=0.5 wraps to frame 0.
	lib.sample(0.5, su, sv, sw); EXPECT_NEAR(su[0], 0.0, 1e-12);
	// t=0.45 interpolates frame 4 → frame 0 (4·(0.5)+0·0.5 = 2.0).
	lib.sample(0.45, su, sv, sw); EXPECT_NEAR(su[0], 2.0, 1e-12);
}

TEST(Precursor, ReplayUploadsInterpolatedPlane)
{
	MacGrid g; g.nx = 4; g.ny = 3; g.nz = 4; g.h = 0.5;
	PrecursorLibrary lib(g, 10.0);
	size_t nu = (size_t)g.ny * g.nz, nv = (size_t)(g.ny + 1) * g.nz, nw = (size_t)g.ny * (g.nz + 1);
	for (int m = 0; m < 4; ++m)
		lib.maybe_record(m * 0.1, std::vector<double>(nu, (double)m), std::vector<double>(nv, 0.0), std::vector<double>(nw, 0.0));
	PrecursorReplay rep(lib);
	rep.advance(0.15); // clock = 0.15 → interp frame1/frame2 → 1.5
	EXPECT_NEAR(rep.u_plane()[0], 1.5, 1e-12);
}
