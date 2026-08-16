#include "core/fluid/amr_eb.h"

#include <cuda_runtime.h>

#include <limits>
#include <stdexcept>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		void check(cudaError_t error, const char* operation)
		{
			if (error != cudaSuccess) throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(error));
		}
		__global__ void gather_kernel(const Real* level, Real* eb, const std::uint64_t* field, int n)
		{
			const int q = blockIdx.x * blockDim.x + threadIdx.x; if (q < n && field[q] != ~std::uint64_t(0)) eb[q] = level[field[q]];
		}
		__global__ void scatter_kernel(const Real* eb, Real* level, const std::uint64_t* field, int n)
		{
			const int q = blockIdx.x * blockDim.x + threadIdx.x; if (q < n && field[q] != ~std::uint64_t(0)) level[field[q]] = eb[q];
		}
	}

	DeviceOneLevelPressureMap::DeviceOneLevelPressureMap(const OneLevelEmbeddedBoundary& composite, const EbPressureSystem& system)
	{
		if (system.eb != &composite.topology) throw std::invalid_argument("pressure system/topology mismatch");
		storage_size_ = system.storage_size; std::vector<std::uint64_t> map(storage_size_, ~std::uint64_t(0));
		for (int cell = 0; cell < composite.topology.grid.cell_count(); ++cell)
		{
			const int dof = system.cell_dof[cell]; if (dof >= 0) map[dof] = static_cast<std::uint64_t>(composite.cell_field_index[cell]);
		}
		if (storage_size_ > 0)
		{
			check(cudaMalloc(&field_index_, map.size() * sizeof(std::uint64_t)), "cudaMalloc one-level pressure map");
			check(cudaMemcpy(field_index_, map.data(), map.size() * sizeof(std::uint64_t), cudaMemcpyHostToDevice), "upload one-level pressure map");
			bytes_ = map.size() * sizeof(std::uint64_t);
		}
	}

	DeviceOneLevelPressureMap::~DeviceOneLevelPressureMap() { if (field_index_) cudaFree(field_index_); }
	void DeviceOneLevelPressureMap::gather(const Real* level, Real* eb) const
	{
		if (!storage_size_) return; gather_kernel<<<(storage_size_ + 255) / 256, 256>>>(level, eb, field_index_, storage_size_); check(cudaDeviceSynchronize(), "gather one-level pressure");
	}
	void DeviceOneLevelPressureMap::scatter(const Real* eb, Real* level) const
	{
		if (!storage_size_) return; scatter_kernel<<<(storage_size_ + 255) / 256, 256>>>(eb, level, field_index_, storage_size_); check(cudaDeviceSynchronize(), "scatter one-level pressure");
	}
}
