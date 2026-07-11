// cuda_probe.cu — trivial SAXPY kernel + CPU reference (M0 CUDA toolchain probe).
#include "core/cuda_probe.h"

#include <cuda_runtime.h>

#include <string>
#include <vector>

namespace windcfd::core
{
	void saxpy_cpu(float a, const std::vector<float>& x, const std::vector<float>& y, std::vector<float>& out)
	{
		const std::size_t n = x.size();
		out.resize(n);
		for (std::size_t i = 0; i < n; ++i)
			out[i] = a * x[i] + y[i];
	}

	namespace
	{
		__global__ void saxpy_kernel(float a, const float* x, const float* y, float* out, int n)
		{
			int i = blockIdx.x * blockDim.x + threadIdx.x;
			if (i < n)
				out[i] = a * x[i] + y[i];
		}

		bool cuda_ok(cudaError_t e, const char* what, std::string* err)
		{
			if (e == cudaSuccess)
				return true;
			if (err)
				*err = std::string(what) + ": " + cudaGetErrorString(e);
			return false;
		}
	}

	bool saxpy_gpu(float a, const std::vector<float>& x, const std::vector<float>& y, std::vector<float>& out, std::string* err)
	{
		const int n = static_cast<int>(x.size());
		out.assign(x.size(), 0.0f);
		if (n == 0)
			return true;

		int dev_count = 0;
		if (!cuda_ok(cudaGetDeviceCount(&dev_count), "cudaGetDeviceCount", err))
			return false;
		if (dev_count < 1)
		{
			if (err)
				*err = "no CUDA device found";
			return false;
		}

		float* dx = nullptr;
		float* dy = nullptr;
		float* dout = nullptr;
		const std::size_t bytes = static_cast<std::size_t>(n) * sizeof(float);
		bool ok = true;

		ok = ok && cuda_ok(cudaMalloc(&dx, bytes), "cudaMalloc dx", err);
		ok = ok && cuda_ok(cudaMalloc(&dy, bytes), "cudaMalloc dy", err);
		ok = ok && cuda_ok(cudaMalloc(&dout, bytes), "cudaMalloc dout", err);
		ok = ok && cuda_ok(cudaMemcpy(dx, x.data(), bytes, cudaMemcpyHostToDevice), "cudaMemcpy dx", err);
		ok = ok && cuda_ok(cudaMemcpy(dy, y.data(), bytes, cudaMemcpyHostToDevice), "cudaMemcpy dy", err);

		if (ok)
		{
			const int block = 256;
			const int grid = (n + block - 1) / block;
			saxpy_kernel<<<grid, block>>>(a, dx, dy, dout, n);
			ok = cuda_ok(cudaGetLastError(), "kernel launch", err);
			ok = ok && cuda_ok(cudaDeviceSynchronize(), "cudaDeviceSynchronize", err);
		}

		ok = ok && cuda_ok(cudaMemcpy(out.data(), dout, bytes, cudaMemcpyDeviceToHost), "cudaMemcpy out", err);

		cudaFree(dx);
		cudaFree(dy);
		cudaFree(dout);
		return ok;
	}

	std::string cuda_device_name(std::string* err)
	{
		int dev_count = 0;
		if (!cuda_ok(cudaGetDeviceCount(&dev_count), "cudaGetDeviceCount", err) || dev_count < 1)
		{
			if (err && err->empty())
				*err = "no CUDA device found";
			return {};
		}
		cudaDeviceProp prop{};
		if (!cuda_ok(cudaGetDeviceProperties(&prop, 0), "cudaGetDeviceProperties", err))
			return {};
		return std::string(prop.name);
	}

	int cuda_probe_usable(std::string* name, int* cc_major, int* cc_minor, std::string* err)
	{
		int dev_count = 0;
		if (!cuda_ok(cudaGetDeviceCount(&dev_count), "cudaGetDeviceCount", err) || dev_count < 1)
		{
			if (err && err->empty())
				*err = "no CUDA device found";
			return 1; // no device, or the driver is too old (cudaGetDeviceCount surfaces both)
		}
		cudaDeviceProp prop{};
		if (!cuda_ok(cudaGetDeviceProperties(&prop, 0), "cudaGetDeviceProperties", err))
			return 1;
		if (name)
			*name = prop.name;
		if (cc_major)
			*cc_major = prop.major;
		if (cc_minor)
			*cc_minor = prop.minor;

		// Actually EXECUTE one of this build's kernels. The driver + device are fine at this point (the
		// calls above passed), so the only realistic failure of a 3-element SAXPY is that this GPU's
		// compute capability is not among the compiled archs and no PTX is embedded to JIT — i.e.
		// cudaErrorNoKernelImageForDevice. That is exactly the "GPU too old to run our kernels" case.
		const std::vector<float> x{ 1.0f, 2.0f, 3.0f }, y{ 1.0f, 1.0f, 1.0f };
		std::vector<float> out;
		std::string kerr;
		if (!saxpy_gpu(2.0f, x, y, out, &kerr))
		{
			if (err)
				*err = kerr;
			return 2; // a device is present, but it cannot run this build's kernels (arch too old)
		}
		return 0; // usable
	}
}
