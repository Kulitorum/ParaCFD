// cuda_probe.h — minimal CUDA toolchain probe for the M0 scaffold.
//
// A trivial SAXPY kernel with a CPU reference, so the M0 build actually compiles
// and runs device code on the RTX 4090 and the GPU-vs-CPU test harness pattern
// (mandated for every kernel, CLAUDE.md) exists from day one.
#pragma once

#include <string>
#include <vector>

namespace paracfd::core
{
	// out[i] = a*x[i] + y[i], CPU reference.
	void saxpy_cpu(float a, const std::vector<float>& x, const std::vector<float>& y, std::vector<float>& out);

	// out[i] = a*x[i] + y[i] on the GPU. Returns false and fills *err on any CUDA
	// error (no device, alloc/launch/copy failure). out is resized to x.size().
	bool saxpy_gpu(float a, const std::vector<float>& x, const std::vector<float>& y, std::vector<float>& out, std::string* err = nullptr);

	// Name of CUDA device 0, or empty string if none / error (message in *err).
	std::string cuda_device_name(std::string* err = nullptr);

	// Startup GPU usability probe. Fills `name` + the compute capability (cc_major.cc_minor) when a
	// device is present, and actually LAUNCHES a trivial kernel to confirm this build's compiled
	// architectures can execute on it. Returns:
	//   0 = usable;
	//   1 = no CUDA device / driver too old (reason in *err);
	//   2 = a device is present but too old to run this build's kernels — its compute capability is not
	//       among the compiled archs, so the launch returned "no kernel image" (name/cc_* are still set).
	// This catches the "device exists but our sm_XX kernels can't run on it" case that cuda_device_name
	// misses, so the GUI can show a clear "GPU not supported" message instead of dying at the first kernel.
	int cuda_probe_usable(std::string* name, int* cc_major, int* cc_minor, std::string* err = nullptr);
}
