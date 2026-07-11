// test_porous.cu — GPU-vs-CPU parity for the sub-grid porous momentum sink (channel_porous.*),
// plus the physical limits of the drag (β=1 open ⇒ no-op; larger k ⇒ more deceleration; the
// implicit form never reverses or overshoots). CLAUDE.md hard rule: every kernel has a CPU twin
// and a GPU-vs-CPU test (rel. max-norm 1e-5 for doubles here → tighter).
#include "core/fluid/channel_porous.h"
#include "core/sediment/sed_physics.h" // screen_resistance_k

#include <gtest/gtest.h>

#include <cuda_runtime.h>

#include <cmath>
#include <random>
#include <vector>

using namespace scour::core;

namespace
{
	double* up(const std::vector<double>& h)
	{ double* d = nullptr; cudaMalloc(&d, sizeof(double) * h.size()); cudaMemcpy(d, h.data(), sizeof(double) * h.size(), cudaMemcpyHostToDevice); return d; }
	void down(std::vector<double>& h, const double* d)
	{ cudaMemcpy(h.data(), d, sizeof(double) * h.size(), cudaMemcpyDeviceToHost); }
}

TEST(Porous, DragLimitsAndMonotonicity)
{
	// β = 1 (k=0) ⇒ velocity unchanged.
	EXPECT_DOUBLE_EQ(porous_face_drag(0.5, 0.0, 0.01, 0.02), 0.5);
	// Positive k decelerates but keeps sign; larger k ⇒ smaller magnitude.
	double u0 = 0.8;
	double k1 = screen_resistance_k(0.35), k2 = screen_resistance_k(0.25); // k2 > k1
	double a = porous_face_drag(u0, k1, 0.01, 0.004), b = porous_face_drag(u0, k2, 0.01, 0.004);
	EXPECT_GT(a, 0.0); EXPECT_LT(a, u0);
	EXPECT_LT(b, a); // stiffer screen ⇒ slower
	// Symmetric in sign.
	EXPECT_NEAR(porous_face_drag(-u0, k1, 0.01, 0.004), -a, 1e-15);
	// Implicit form is unconditionally stable: never overshoots zero even for a huge k·dt.
	double huge = porous_face_drag(1.0, 1e6, 1.0, 0.001);
	EXPECT_GE(huge, 0.0); EXPECT_LT(huge, 1.0);
}

TEST(Porous, GpuMatchesCpu)
{
	MacGrid g; g.nx = 12; g.ny = 10; g.nz = 8; g.h = 0.004;
	int nu = g.u_count(), nv = g.v_count(), nw = g.w_count(), np = g.p_count();
	std::mt19937 rng(1234);
	std::uniform_real_distribution<double> uni(-1.0, 1.0), pos(0.0, 1.0);

	std::vector<double> u(nu), v(nv), w(nw), pk(np, 0.0);
	for (auto& x : u) x = uni(rng);
	for (auto& x : v) x = uni(rng);
	for (auto& x : w) x = uni(rng);
	// A porous block in the middle third with a realistic β≈0.3 resistance; rest open.
	for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		if (i >= 4 && i < 8 && j >= 3 && j < 7) pk[g.pidx(i, j, k)] = screen_resistance_k(0.30);

	std::vector<double> uc = u, vc = v, wc = w;
	ch_porous_drag_cpu(uc, vc, wc, pk, g, 0.02);

	double *du = up(u), *dv = up(v), *dw = up(w), *dpk = up(pk);
	ch_porous_drag_gpu(du, dv, dw, dpk, g, 0.02);
	cudaDeviceSynchronize();
	std::vector<double> ug(nu), vg(nv), wg(nw);
	down(ug, du); down(vg, dv); down(wg, dw);
	cudaFree(du); cudaFree(dv); cudaFree(dw); cudaFree(dpk);

	double maxerr = 0.0;
	for (int i = 0; i < nu; ++i) maxerr = std::max(maxerr, std::fabs(ug[i] - uc[i]));
	for (int i = 0; i < nv; ++i) maxerr = std::max(maxerr, std::fabs(vg[i] - vc[i]));
	for (int i = 0; i < nw; ++i) maxerr = std::max(maxerr, std::fabs(wg[i] - wc[i]));
	EXPECT_LT(maxerr, 1e-12);

	// Sanity: velocities inside the porous block are reduced vs the input (drag acted).
	int fc = g.uidx(6, 5, 4);
	EXPECT_LT(std::fabs(ug[fc]), std::fabs(u[fc]) + 1e-15);
}
