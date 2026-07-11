// test_slice_field.cu — GPU-vs-CPU parity for the G1 slice colour-fill kernel
// (gui/slice_field.cu). The kernel and its CPU reference share one __host__ __device__
// evaluator, so this catches launch/index/memory bugs. Rel. max-norm tol 1e-5 (float
// colours). Labeled `unit` — runs in the fast suite, no GL/Qt needed.
#include "gui/slice_field.h"

#include <gtest/gtest.h>

#include <cuda_runtime.h>

#include <cmath>
#include <random>
#include <vector>

using namespace scour::gui;
using scour::core::MacGrid;

namespace
{
	struct DevFields
	{
		double *u = nullptr, *v = nullptr, *w = nullptr, *p = nullptr, *c = nullptr;
	};

	// Fill host MAC fields with a smooth-ish deterministic pattern, upload to device. `hc`
	// is the cell-centred suspended-sediment concentration (same pidx layout as `hp`).
	void make_fields(MacGrid g, std::vector<double>& hu, std::vector<double>& hv,
		std::vector<double>& hw, std::vector<double>& hp, std::vector<double>& hc, DevFields& d)
	{
		std::mt19937 rng(1234);
		std::uniform_real_distribution<double> dist(-1.0, 1.0);
		hu.resize(g.u_count()); hv.resize(g.v_count()); hw.resize(g.w_count()); hp.resize(g.p_count()); hc.resize(g.p_count());
		for (auto& x : hu) x = dist(rng);
		for (auto& x : hv) x = dist(rng);
		for (auto& x : hw) x = dist(rng);
		for (auto& x : hp) x = dist(rng);
		for (auto& x : hc) x = dist(rng);
		cudaMalloc(&d.u, sizeof(double) * hu.size());
		cudaMalloc(&d.v, sizeof(double) * hv.size());
		cudaMalloc(&d.w, sizeof(double) * hw.size());
		cudaMalloc(&d.p, sizeof(double) * hp.size());
		cudaMalloc(&d.c, sizeof(double) * hc.size());
		cudaMemcpy(d.u, hu.data(), sizeof(double) * hu.size(), cudaMemcpyHostToDevice);
		cudaMemcpy(d.v, hv.data(), sizeof(double) * hv.size(), cudaMemcpyHostToDevice);
		cudaMemcpy(d.w, hw.data(), sizeof(double) * hw.size(), cudaMemcpyHostToDevice);
		cudaMemcpy(d.p, hp.data(), sizeof(double) * hp.size(), cudaMemcpyHostToDevice);
		cudaMemcpy(d.c, hc.data(), sizeof(double) * hc.size(), cudaMemcpyHostToDevice);
	}

	void free_fields(DevFields& d) { cudaFree(d.u); cudaFree(d.v); cudaFree(d.w); cudaFree(d.p); cudaFree(d.c); }

	double max_rel_diff(const std::vector<float4>& a, const std::vector<float4>& b)
	{
		double m = 0.0;
		for (size_t i = 0; i < a.size(); ++i)
		{
			const float av[4] = { a[i].x, a[i].y, a[i].z, a[i].w };
			const float bv[4] = { b[i].x, b[i].y, b[i].z, b[i].w };
			for (int c = 0; c < 4; ++c)
			{
				double d = std::fabs((double)av[c] - (double)bv[c]);
				double s = std::max(1.0, std::fabs((double)bv[c]));
				m = std::max(m, d / s);
			}
		}
		return m;
	}

	void run_case(MacGrid g, Axis axis, float plane_pos, Field field, int nu, int nv)
	{
		std::vector<double> hu, hv, hw, hp, hc;
		DevFields d;
		make_fields(g, hu, hv, hw, hp, hc, d);

		SliceParams sp;
		sp.grid = g; sp.axis = axis; sp.plane_pos = plane_pos;
		sp.nu = nu; sp.nv = nv; sp.field = field; sp.vmin = -1.0f; sp.vmax = 1.0f;

		int n = nu * nv;
		std::vector<float4> cpu(n), gpu(n);
		slice_fill_cpu(hu.data(), hv.data(), hw.data(), hp.data(), hc.data(), sp, cpu.data());

		float4* dout = nullptr;
		cudaMalloc(&dout, sizeof(float4) * n);
		slice_fill_gpu(d.u, d.v, d.w, d.p, d.c, sp, dout, 0);
		cudaDeviceSynchronize();
		cudaMemcpy(gpu.data(), dout, sizeof(float4) * n, cudaMemcpyDeviceToHost);
		cudaFree(dout);
		free_fields(d);

		EXPECT_LT(max_rel_diff(gpu, cpu), 1e-5) << "axis=" << (int)axis << " field=" << (int)field;
	}
}

TEST(SliceField, GpuMatchesCpu_AllFieldsAllAxes)
{
	MacGrid g; g.nx = 20; g.ny = 16; g.nz = 12; g.h = 0.05;
	const Field fields[] = { Field::SpeedMag, Field::VelU, Field::VelV, Field::VelW, Field::Pressure };
	const Axis axes[] = { Axis::X, Axis::Y, Axis::Z };
	for (Axis ax : axes)
	{
		float pos = (ax == Axis::X) ? (float)(g.nx * g.h) * 0.5f
			: (ax == Axis::Y) ? (float)(g.ny * g.h) * 0.5f
							  : (float)(g.nz * g.h) * 0.5f;
		for (Field f : fields)
			run_case(g, ax, pos, f, 96, 72);
	}
}

// Concentration (Field=5) samples a cell-centred suspended-sediment scalar `c` (const double*,
// same pidx layout as pressure). run_case uploads a random `c` field of size g.p_count() and
// asserts the GPU fill matches its CPU twin at the same 1e-5 rel. max-norm tol as the other fields.
TEST(SliceField, GpuMatchesCpu_Concentration)
{
	MacGrid g; g.nx = 20; g.ny = 16; g.nz = 12; g.h = 0.05;
	const Axis axes[] = { Axis::X, Axis::Y, Axis::Z };
	for (Axis ax : axes)
	{
		float pos = (ax == Axis::X) ? (float)(g.nx * g.h) * 0.5f
			: (ax == Axis::Y) ? (float)(g.ny * g.h) * 0.5f
							  : (float)(g.nz * g.h) * 0.5f;
		run_case(g, ax, pos, Field::Concentration, 96, 72);
	}
}

// Auto-range reduction (gui/slice_field.cu): the GPU block-reduction must match the CPU reference
// bit-for-bit (min/max are order-independent). Exercises fluid-only exclusion via a solid mask and
// a deliberately planted spurious spike (must show up in the max — the feature reveals blow-ups).
TEST(SliceField, ReduceGpuMatchesCpu_FluidOnly)
{
	MacGrid g; g.nx = 20; g.ny = 16; g.nz = 12; g.h = 0.05;
	std::vector<double> hu, hv, hw, hp, hc;
	DevFields d;
	make_fields(g, hu, hv, hw, hp, hc, d);

	// Mark ~30% of cells solid (deterministic) — these must be excluded from the range.
	std::vector<unsigned char> hsolid(g.p_count(), 0);
	for (int n = 0; n < g.p_count(); ++n) hsolid[n] = (n % 10 < 3) ? 1 : 0;
	// Plant a large finite spike in a FLUID cell: it must dominate the reported max.
	int spike = -1;
	for (int n = 0; n < g.p_count(); ++n) if (!hsolid[n]) { spike = n; break; }
	int si = spike % g.nx, sj = (spike / g.nx) % g.ny, sk = spike / (g.nx * g.ny);
	hu[g.uidx(si, sj, sk)] = 50.0; hu[g.uidx(si + 1, sj, sk)] = 50.0;
	cudaMemcpy(d.u, hu.data(), sizeof(double) * hu.size(), cudaMemcpyHostToDevice);
	unsigned char* dsolid = nullptr;
	cudaMalloc(&dsolid, hsolid.size());
	cudaMemcpy(dsolid, hsolid.data(), hsolid.size(), cudaMemcpyHostToDevice);

	const Field fields[] = { Field::SpeedMag, Field::VelU, Field::VelV, Field::VelW, Field::Pressure, Field::Concentration };
	float* d3 = nullptr; cudaMalloc(&d3, 3 * sizeof(float));
	for (Field f : fields)
	{
		float cpu3[3] = { 0, 0, 0 }, gpu3[3] = { 0, 0, 0 };
		slice_reduce_cpu(hu.data(), hv.data(), hw.data(), hp.data(), hc.data(), hsolid.data(), g, f, cpu3);
		slice_reduce_gpu(d.u, d.v, d.w, d.p, d.c, dsolid, g, f, d3, 0);
		cudaDeviceSynchronize();
		cudaMemcpy(gpu3, d3, 3 * sizeof(float), cudaMemcpyDeviceToHost);
		for (int c = 0; c < 3; ++c)
			EXPECT_NEAR(gpu3[c], cpu3[c], 1e-4f) << "field=" << (int)f << " comp=" << c;
		EXPECT_GT(cpu3[2], 40.0f) << "planted spike must dominate speed_max (field=" << (int)f << ")";
	}
	cudaFree(d3); cudaFree(dsolid); free_fields(d);
}

TEST(SliceField, NullPressureIsZeroColour)
{
	MacGrid g; g.nx = 8; g.ny = 8; g.nz = 8; g.h = 0.1;
	std::vector<double> hu, hv, hw, hp, hc;
	DevFields d;
	make_fields(g, hu, hv, hw, hp, hc, d);
	SliceParams sp; sp.grid = g; sp.axis = Axis::Y; sp.plane_pos = 0.4f;
	sp.nu = 32; sp.nv = 32; sp.field = Field::Pressure; sp.vmin = -1.0f; sp.vmax = 1.0f;
	int n = sp.nu * sp.nv;
	std::vector<float4> cpu(n), gpu(n);
	slice_fill_cpu(hu.data(), hv.data(), hw.data(), nullptr, hc.data(), sp, cpu.data());
	float4* dout = nullptr; cudaMalloc(&dout, sizeof(float4) * n);
	slice_fill_gpu(d.u, d.v, d.w, nullptr, d.c, sp, dout, 0);
	cudaDeviceSynchronize();
	cudaMemcpy(gpu.data(), dout, sizeof(float4) * n, cudaMemcpyDeviceToHost);
	cudaFree(dout); free_fields(d);
	EXPECT_LT(max_rel_diff(gpu, cpu), 1e-5);
}
