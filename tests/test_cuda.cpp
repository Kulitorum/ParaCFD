// test_cuda.cpp — CUDA toolchain probe: GPU SAXPY vs CPU reference (rel tol 1e-5).
#include "core/cuda_probe.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

using namespace scour::core;

TEST(Cuda, DevicePresent)
{
	std::string err;
	const std::string name = cuda_device_name(&err);
	// The M0 gate runs on the RTX 4090; a missing device is a real failure.
	ASSERT_FALSE(name.empty()) << "CUDA device query failed: " << err;
	std::printf("[cuda] device 0: %s\n", name.c_str());
}

TEST(Cuda, SaxpyMatchesCpuReference)
{
	const int n = 1 << 20; // 1,048,576 elements
	std::vector<float> x(n), y(n);
	for (int i = 0; i < n; ++i)
	{
		x[i] = std::sin(0.001f * i) * 3.0f;
		y[i] = std::cos(0.0007f * i) * 2.0f - 1.0f;
	}
	const float a = 2.75f;

	std::vector<float> ref, gpu;
	saxpy_cpu(a, x, y, ref);

	std::string err;
	ASSERT_TRUE(saxpy_gpu(a, x, y, gpu, &err)) << err;
	ASSERT_EQ(gpu.size(), ref.size());

	// GPU contracts a*x+y into a single-rounding FMA while the CPU rounds twice;
	// at points where a*x ~ -y this cancellation makes pointwise relative error
	// ill-defined. Compare in the max-norm normalized by the field scale (the same
	// normalized-error convention the physics gates use), rel. tol 1e-5.
	double max_abs = 0.0;
	double scale = 0.0;
	for (int i = 0; i < n; ++i)
	{
		max_abs = std::max(max_abs, std::fabs(static_cast<double>(gpu[i]) - static_cast<double>(ref[i])));
		scale = std::max(scale, std::fabs(static_cast<double>(ref[i])));
	}
	ASSERT_GT(scale, 0.0);
	const double rel = max_abs / scale;
	EXPECT_LT(rel, 1e-5) << "max-norm GPU-vs-CPU error normalized by field scale (max_abs=" << max_abs << ", scale=" << scale << ")";
}
