#pragma once

#include "core/aero_loads.h"
#include "core/fluid/amr_advection.h"
#include "core/fluid/amr_eb.h"
#include "core/fluid/amr_pressure.h"
#include "core/paraglider_config.h"

#include <cstddef>
#include <memory>

namespace paracfd::core
{
	struct ExternalAeroStepStats
	{
		double dt = 0.0;
		double physical_time = 0.0;
		double advection_ms = 0.0;
		double turbulence_ms = 0.0;
		double projection_ms = 0.0;
		double gpu_step_ms = 0.0;
		AmrGpuSolveResult pressure;
		bool first_order_fabric_protection = true;
		bool les_applied = false;
	};

	// GPU-native static-geometry paraglider flow core. CAD/BVH/EB work happens once in
	// the constructor. initialize() and step() retain all fields and pressure topology on
	// the device; only scalar PCG reductions return to the host during a timestep.
	class ExternalAeroCore
	{
	public:
		ExternalAeroCore(const TriMesh& placed_wing, const TriangleBvh& bvh,
			const ParagliderConfig& config);
		~ExternalAeroCore();
		ExternalAeroCore(const ExternalAeroCore&) = delete;
		ExternalAeroCore& operator=(const ExternalAeroCore&) = delete;

		ExternalAeroStepStats initialize();
		ExternalAeroStepStats step();
		AerodynamicLoads pressure_loads(double pressure_reference = 0.0) const;
		void download_fields(AmrHostFields& host) const;
		double max_abs_divergence() const; // throttled validation/statistics download

		const AmrHierarchy& hierarchy() const { return hierarchy_; }
		const AmrEmbeddedBoundaryAtlas& embedded_boundary() const { return embedded_boundary_; }
		const CompositeAmrPressureSystem& pressure_system() const { return pressure_system_; }
		double physical_time() const { return physical_time_; }
		bool initialized() const { return initialized_; }
		std::size_t gpu_bytes() const;
		std::size_t protected_face_count() const;
		std::size_t active_face_count() const;

	private:
		ExternalAeroStepStats project(bool warm_start);

		ParagliderConfig config_;
		std::size_t source_triangle_count_ = 0;
		AmrHierarchy hierarchy_;
		AmrEmbeddedBoundaryAtlas embedded_boundary_;
		CompositeAmrPressureSystem pressure_system_;
		std::unique_ptr<DeviceAmrFields> fields_;
		std::unique_ptr<DeviceAmrAdvection> advection_;
		std::unique_ptr<DeviceCompositeAmrProjection> projection_;
		double physical_time_ = 0.0;
		bool initialized_ = false;
	};
}
