#pragma once

#include "core/fluid/amr_grid.h"
#include "core/fluid/real.h"

#include <cstddef>
#include <cstdint>
#include <vector>

#ifdef __CUDACC__
#define PARACFD_AMR_HD __host__ __device__
#else
#define PARACFD_AMR_HD
#endif

namespace paracfd::core
{
	struct BrickFieldLayout
	{
		int brick_size = 0, ghost = 0;
		int cell_n = 0, u_nx = 0, v_ny = 0, w_nz = 0;
		std::size_t cell_stride = 0, u_stride = 0, v_stride = 0, w_stride = 0;

		static BrickFieldLayout make(int brick_size, int ghost);
		PARACFD_AMR_HD std::size_t cell_index(int brick, int i, int j, int k) const { i+=ghost;j+=ghost;k+=ghost;return static_cast<std::size_t>(brick)*cell_stride+(static_cast<std::size_t>(k)*cell_n+j)*cell_n+i; }
		PARACFD_AMR_HD std::size_t u_index(int brick, int i, int j, int k) const { i+=ghost;j+=ghost;k+=ghost;return static_cast<std::size_t>(brick)*u_stride+(static_cast<std::size_t>(k)*cell_n+j)*u_nx+i; }
		PARACFD_AMR_HD std::size_t v_index(int brick, int i, int j, int k) const { i+=ghost;j+=ghost;k+=ghost;return static_cast<std::size_t>(brick)*v_stride+(static_cast<std::size_t>(k)*v_ny+j)*cell_n+i; }
		PARACFD_AMR_HD std::size_t w_index(int brick, int i, int j, int k) const { i+=ghost;j+=ghost;k+=ghost;return static_cast<std::size_t>(brick)*w_stride+(static_cast<std::size_t>(k)*cell_n+j)*cell_n+i; }
	};

	struct AmrHostLevelFields
	{
		BrickFieldLayout layout;
		std::vector<Real> u, v, w, p, nut, temp;
	};

	// Non-owning view of one persistent device pool. It is intentionally plain data so
	// composite CUDA launchers can upload a compact array of level views without gaining
	// access to DeviceAmrFields ownership or introducing per-brick allocations.
	struct DeviceAmrFieldLevelView
	{
		BrickFieldLayout layout;
		int brick_count = 0;
		Real *u = nullptr, *v = nullptr, *w = nullptr, *p = nullptr, *nut = nullptr, *temp = nullptr;
		const int* neighbors = nullptr;
		const std::uint32_t* flags = nullptr;
	};

	class AmrHostFields
	{
	public:
		explicit AmrHostFields(const AmrHierarchy& hierarchy);
		std::vector<AmrHostLevelFields>& levels() { return levels_; }
		const std::vector<AmrHostLevelFields>& levels() const { return levels_; }
		std::size_t bytes() const;
		void exchange_same_level_pressure_halos(const AmrHierarchy& hierarchy);

	private:
		std::vector<AmrHostLevelFields> levels_;
	};

	// Persistent per-level GPU allocations. Each field is one allocation per level, laid out
	// [brick][contiguous brick field], rather than thousands of per-brick allocations.
	class DeviceAmrFields
	{
	public:
		explicit DeviceAmrFields(const AmrHierarchy& hierarchy);
		~DeviceAmrFields();
		DeviceAmrFields(const DeviceAmrFields&) = delete;
		DeviceAmrFields& operator=(const DeviceAmrFields&) = delete;

		void upload(const AmrHostFields& host);
		void download(AmrHostFields& host) const;
		void download_pressure(AmrHostFields& host) const;
		void fill(Real value);
		void initialize_freestream(Real speed);
		// External aerodynamic BC: prescribed +X inflow, zero-pressure X-max,
		// and symmetric free-slip far field on both Y and both Z boundaries.
		// There is deliberately no ground-plane branch.
		void apply_external_aero_boundaries(Real freestream_speed);
		void exchange_same_level_pressure_halos();
		// Physical faces of active bricks only: ghost storage and covered coarse bricks
		// cannot spuriously throttle the one-global-step CFL reduction.
		double max_abs_velocity() const; // one scalar D2H reduction
		std::size_t bytes() const { return bytes_; }
		Real* pressure(int level) { return levels_[level].p; }
		const Real* pressure(int level) const { return levels_[level].p; }
		int level_count() const { return static_cast<int>(levels_.size()); }
		DeviceAmrFieldLevelView level_view(int level)
		{
			auto& d=levels_[level];return {d.layout,d.brick_count,d.u,d.v,d.w,d.p,d.nut,d.temp,d.neighbors,d.flags};
		}

	private:
		struct Level
		{
			BrickFieldLayout layout;
			int brick_count = 0;
			Real *u = nullptr, *v = nullptr, *w = nullptr, *p = nullptr, *nut = nullptr, *temp = nullptr;
			int* neighbors = nullptr; // [brick][6]
			std::uint32_t* flags = nullptr; // [brick]
		};
		std::vector<Level> levels_;
		Real* max_abs_scratch_ = nullptr;
		std::size_t bytes_ = 0;
	};

	// Shared positive max-absolute reduction for persistent Real device arrays. `scratch`
	// is one persistent device Real owned by the caller; only the scalar result is copied.
	double max_abs_device_values(const Real* values, std::size_t count, Real* scratch);

	// Compact device-side finest-brick locator. Each AMR level owns an open-addressed
	// integer-coordinate table; a point probes levels finest-to-coarsest. Normal sampling
	// kernels therefore do not traverse an octree or binary-search graded coordinates.
	struct GpuAmrPoint { float x = 0, y = 0, z = 0; };
	struct GpuBrickLocation
	{
		int level = -1, brick = -1;
		int i = 0, j = 0, k = 0;
	};
	struct GpuBrickLookupEntry { int x = 0, y = 0, z = 0, brick = -1; };
	struct GpuBrickRecord
	{
		float origin_x = 0, origin_y = 0, origin_z = 0, h = 0;
		std::uint32_t flags = 0;
	};
	struct GpuAmrLevelView
	{
		const GpuBrickLookupEntry* lookup = nullptr;
		const GpuBrickRecord* bricks = nullptr;
		std::uint32_t lookup_mask = 0;
		int lookup_size = 0, brick_count = 0;
		float h = 0;
	};
	struct GpuAmrHierarchyView
	{
		const GpuAmrLevelView* levels = nullptr;
		int level_count = 0, brick_size = 0;
		GpuAmrPoint domain_lo{}, domain_hi{};
	};

	class DeviceAmrLocator
	{
	public:
		explicit DeviceAmrLocator(const AmrHierarchy& hierarchy);
		~DeviceAmrLocator();
		DeviceAmrLocator(const DeviceAmrLocator&) = delete;
		DeviceAmrLocator& operator=(const DeviceAmrLocator&) = delete;

		// `points` and `locations` are device arrays. This batched entry point is also a
		// validation helper; CFD kernels can consume view() and use the same locator inline.
		void locate_points(const GpuAmrPoint* points, GpuBrickLocation* locations, int count) const;
		GpuAmrHierarchyView view() const { return view_; }
		std::size_t bytes() const { return bytes_; }

	private:
		struct Allocation
		{
			GpuBrickLookupEntry* lookup = nullptr;
			GpuBrickRecord* bricks = nullptr;
		};
		std::vector<Allocation> allocations_;
		GpuAmrLevelView* device_levels_ = nullptr;
		GpuAmrHierarchyView view_{};
		std::size_t bytes_ = 0;
	};
}
