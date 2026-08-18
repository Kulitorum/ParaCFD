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
		double max_diffusion_rate = 0.0; // 1/s, explicit finite-volume graph spectral bound
		double cfl_velocity = 0.0;
		double effective_cfl = 0.0;
		int diffusion_substeps = 1;
		AmrGpuSolveResult pressure;
		bool side_safe_fabric_transport = true;
		bool les_applied = false;
		bool embedded_transport_applied = false;
		bool smooth_fabric_wall_applied = false;
		bool conservative_cell_momentum = false;
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
	struct ExternalAeroExecutionOptions
	{
		// The face-centred MAC solver is the validated production path. The
		// collocated control-volume transport remains an explicit experiment.
		bool conservative_cell_momentum = false;
		bool smooth_fabric_wall = true;
		bool pressure_impulse = true;
		// Test harness only: reproduce pressure behavior on the historical sampled
		// topology. Any front end enabling this must label the result qualitative.
		bool allow_unsafe_same_fragment_patches = false;
		// Keep the conservative two-point A/d pressure flux when the legacy preview
		// topology cannot support its deferred non-orthogonal WLS correction. Strict
		// CAD-certified runs leave this false and reject the topology instead.
		bool use_qualitative_first_order_orthogonal_pressure = false;
		// Static geometry callback supplied by front ends that link
		// paracfd_geometry. The solver library itself remains OpenCascade-free.
		ExactCellDecomposer exact_cell_decomposer = nullptr;
	};
	// GPU-native static-geometry paraglider flow core. CAD/BVH/EB work happens once in
	// the constructor. initialize() and step() retain all fields and pressure topology on
	// the device; only scalar PCG reductions return to the host during a timestep.
	class ExternalAeroCore
	{
	public:
		ExternalAeroCore(const TriMesh& placed_wing, const TriangleBvh& bvh,
			const ParagliderConfig& config, ExternalAeroExecutionOptions options = {});
		~ExternalAeroCore();
		ExternalAeroCore(const ExternalAeroCore&) = delete;
		ExternalAeroCore& operator=(const ExternalAeroCore&) = delete;

		ExternalAeroStepStats initialize();
		ExternalAeroStepStats step();
		AerodynamicLoads pressure_loads(double pressure_reference = 0.0) const;
		AerodynamicLoads aerodynamic_loads(double pressure_reference = 0.0) const;
		void download_fields(AmrHostFields& host) const;
		void download_special_fluxes(CompositeAmrFluxes& host) const;
		void download_pressure(std::vector<double>& host) const; // throttled validation/debug download
		bool download_cell_momentum_state(CompositeCellMomentumState& host) const; // diagnostic only
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
		int fabric_wall_node_count() const;
		int clamped_cell_flux_interpolation_count() const;
		int pressure_closure_correction_count() const;
		double max_pressure_closure_acceleration() const;
		bool uses_conservative_cell_momentum() const { return use_conservative_cell_momentum_; }
		bool uses_half_wing_symmetry() const { return config_.domain.half_wing_symmetry; }

	private:
		ExternalAeroStepStats project(bool warm_start, double dt, bool sync_coarse_fine = true);
		AerodynamicLoads unscaled_pressure_loads(double pressure_reference) const;

		ParagliderConfig config_;
		std::size_t source_triangle_count_ = 0;
		AmrHierarchy hierarchy_;
		AmrEmbeddedBoundaryAtlas embedded_boundary_;
		CompositeAmrPressureSystem pressure_system_;
		std::unique_ptr<DeviceAmrFields> fields_;
		std::unique_ptr<DeviceAmrAdvection> advection_;
		std::unique_ptr<DeviceCompositeAmrProjection> projection_;
		std::unique_ptr<DeviceCompositeCellMomentumTransport> cell_momentum_;
		double physical_time_ = 0.0;
		double symmetry_plane_y_ = 0.0;
		bool initialized_ = false;
		bool smooth_fabric_wall_applied_ = false;
		bool use_conservative_cell_momentum_ = false;
		bool enable_smooth_fabric_wall_ = true;
		bool enable_pressure_impulse_ = true;
	};
}
