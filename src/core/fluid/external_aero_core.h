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
		double embedded_transport_ms = 0.0;
		double projection_ms = 0.0;
		double gpu_step_ms = 0.0;
		double max_abs_velocity = 0.0;
		double max_abs_regular_velocity = 0.0;
		double max_abs_special_velocity = 0.0;
		double max_embedded_cfl_rate = 0.0; // 1/s, compact EB graph transport
		double cfl_velocity = 0.0;
		double effective_cfl = 0.0;
		AmrGpuSolveResult pressure;
		bool side_safe_fabric_transport = true;
		bool les_applied = false;
		bool embedded_transport_applied = false;
	};
	struct ExternalAeroConservationStats
	{
		double max_abs_divergence = 0.0;            // 1/s
		double volume_weighted_rms_divergence = 0.0; // 1/s
		double worst_control_volume = 0.0;          // m^3 at max_abs_divergence
		double max_integrated_flux_error = 0.0;     // m^3/s for one control volume
		double absolute_integrated_flux_error = 0.0; // sum |divergence * volume|, m^3/s
		double net_integrated_flux_error = 0.0;     // sum divergence * volume, m^3/s
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
		void download_special_fluxes(CompositeAmrFluxes& host) const;
		void download_pressure(std::vector<double>& host) const; // throttled validation/debug download
		ExternalAeroConservationStats conservation_stats() const; // throttled validation/statistics download
		double max_abs_divergence() const; // throttled validation/statistics download

		const AmrHierarchy& hierarchy() const { return hierarchy_; }
		const AmrEmbeddedBoundaryAtlas& embedded_boundary() const { return embedded_boundary_; }
		const CompositeAmrPressureSystem& pressure_system() const { return pressure_system_; }
		const ParagliderConfig& config() const { return config_; }
		double physical_time() const { return physical_time_; }
		bool initialized() const { return initialized_; }
		std::size_t gpu_bytes() const;
		std::size_t protected_face_count() const;
		std::size_t active_face_count() const;
		int embedded_high_order_stencil_count() const;
		int embedded_least_squares_full_rank_count() const;

	private:
		ExternalAeroStepStats project(bool warm_start, double dt);

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
