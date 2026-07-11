// gpu_testutil.h — tiny RAII device-buffer + comparison helpers for the GPU-vs-CPU
// kernel tests (CLAUDE.md: every CUDA kernel has a CPU reference + a GPU-vs-CPU test).
#pragma once

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <vector>

template <class T>
struct DevVec
{
	T* p = nullptr;
	size_t n = 0;
	DevVec() = default;
	explicit DevVec(size_t n_) { alloc(n_); }
	explicit DevVec(const std::vector<T>& h) { alloc(h.size()); upload(h); }
	void alloc(size_t n_) { n = n_; cudaMalloc(&p, n * sizeof(T)); cudaMemset(p, 0, n * sizeof(T)); }
	void upload(const std::vector<T>& h) { cudaMemcpy(p, h.data(), n * sizeof(T), cudaMemcpyHostToDevice); }
	std::vector<T> download() const { std::vector<T> h(n); cudaMemcpy(h.data(), p, n * sizeof(T), cudaMemcpyDeviceToHost); return h; }
	~DevVec() { if (p) cudaFree(p); }
	DevVec(const DevVec&) = delete;
	DevVec& operator=(const DevVec&) = delete;
};

// Relative error in the max-norm, normalised by the field scale (the convention the
// physics gates use). Returns 0 for two all-zero fields.
inline double rel_maxnorm(const std::vector<double>& a, const std::vector<double>& b)
{
	double maxabs = 0.0, scale = 0.0;
	size_t n = std::min(a.size(), b.size());
	for (size_t i = 0; i < n; ++i)
	{
		maxabs = std::max(maxabs, std::fabs(a[i] - b[i]));
		scale = std::max(scale, std::fabs(b[i]));
	}
	return scale > 0.0 ? maxabs / scale : maxabs;
}
