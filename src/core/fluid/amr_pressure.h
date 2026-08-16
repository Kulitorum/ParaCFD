// amr_pressure.h — matrix-free composite regular-region pressure operator for static AMR.
//
// Ordinary cells retain implicit Cartesian connectivity within and between same-level bricks.
// Only 2:1 coarse/fine interfaces are stored explicitly, so this is not a full-domain CSR graph.
#pragma once

#include "core/fluid/amr_fields.h"
#include "core/fluid/real.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace paracfd::core
{
	struct DeviceCompositeAmrLevelView;
	struct AmrEmbeddedBoundaryAtlas;
	struct CoarseFinePressureConnection
	{
		int coarse_dof = -1;
		int fine_dof = -1;
		double open_area = 0.0;
		double centre_distance = 0.0;
		std::int8_t axis = 0;
		std::int8_t direction = 1; // +1: first DOF is lower-axis; -1: second DOF is lower-axis
		Vec3d face_centroid{};
	};
	struct CompositePressureGauge { int dof = -1; double coefficient = 0.0; };
	struct CompositeSurfacePressurePatch
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		double area = 0.0;
		Vec3d centroid{};
		Vec3d normal{}; // local triangle normal, minus -> plus
		int plus_dof = -1;
		int minus_dof = -1;
	};

	struct CompositeAmrPressureSystem
	{
		const AmrHierarchy* hierarchy = nullptr;
		int brick_size = 0;
		int storage_size = 0;
		bool pressure_outlet_xmax = true;
		std::vector<int> level_offset;
		std::vector<unsigned char> active;
		std::vector<std::uint8_t> cut_face_mask; // +X/+Y/+Z bits on regular Cartesian DOFs
		std::vector<double> volume;
		// Active pressure DOF -> containing level-0 Cartesian cell DOF. This compact
		// geometric aggregation map feeds the matrix-free two-level preconditioner.
		std::vector<int> preconditioner_aggregate;
		std::vector<CoarseFinePressureConnection> coarse_fine;
		std::vector<CoarseFinePressureConnection> embedded;
		// Static CAD provenance and two-sided pressure mapping. These records are
		// consumed only when publishing loads/visualization, not by timestep kernels.
		std::vector<CompositeSurfacePressurePatch> surface_patches;
		// One compact pressure reference for each active fluid component that cannot
		// reach the X-max Dirichlet outlet. Empty for intentionally pure-Neumann tests.
		std::vector<CompositePressureGauge> gauges;

		int dof(int level, int brick, int i, int j, int k) const;
		void apply_cpu(const std::vector<double>& pressure, std::vector<double>& output) const;
		void diagonal_cpu(std::vector<double>& diagonal) const;
	};

	CompositeAmrPressureSystem build_composite_amr_pressure_system(const AmrHierarchy& hierarchy,
		bool pressure_outlet_xmax = true);
	CompositeAmrPressureSystem build_composite_amr_pressure_system(const AmrHierarchy& hierarchy,
		const AmrEmbeddedBoundaryAtlas& embedded_boundary, bool pressure_outlet_xmax = true);

	struct CompositeAmrFluxes
	{
		std::vector<double> coarse_fine_velocity; // one +axis velocity per fine aperture tile
		std::vector<double> embedded_velocity;    // one +axis velocity per EB aperture
	};

	CompositeAmrFluxes make_zero_composite_fluxes(const CompositeAmrPressureSystem& system);
	void composite_amr_divergence_cpu(const CompositeAmrPressureSystem& system, const AmrHostFields& fields,
		const CompositeAmrFluxes& special_flux, std::vector<double>& divergence);
	void composite_amr_projection_rhs_cpu(const CompositeAmrPressureSystem& system, const AmrHostFields& fields,
		const CompositeAmrFluxes& special_flux, double rho, double dt, std::vector<double>& rhs);
	void composite_amr_correct_fluxes_cpu(const CompositeAmrPressureSystem& system, const std::vector<double>& pressure,
		double rho, double dt, AmrHostFields& fields, CompositeAmrFluxes& special_flux);

	class DeviceCompositeAmrPressureOperator
	{
	public:
		explicit DeviceCompositeAmrPressureOperator(const CompositeAmrPressureSystem& system);
		~DeviceCompositeAmrPressureOperator();
		DeviceCompositeAmrPressureOperator(const DeviceCompositeAmrPressureOperator&) = delete;
		DeviceCompositeAmrPressureOperator& operator=(const DeviceCompositeAmrPressureOperator&) = delete;
		void apply(const Real* pressure, Real* output) const;
		int storage_size() const { return storage_size_; }
		std::size_t bytes() const { return bytes_; }

	private:
		struct Allocation { int* neighbors = nullptr; std::uint32_t* flags = nullptr; };
		DeviceCompositeAmrLevelView* levels_ = nullptr;
		std::vector<Allocation> allocations_;
		std::vector<int> brick_counts_;
		int *coarse_dof_ = nullptr, *fine_dof_ = nullptr;
		Real* coefficient_ = nullptr;
		int* gauge_dof_ = nullptr;
		Real* gauge_coefficient_ = nullptr;
		unsigned char *active_ = nullptr;
		std::uint8_t *cut_face_mask_ = nullptr;
		int level_count_ = 0, connection_count_ = 0, gauge_count_ = 0, storage_size_ = 0, brick_size_ = 0;
		bool outlet_ = true;
		std::size_t bytes_ = 0;
	};

	struct AmrGpuSolveResult { int iterations = 0; double relative_residual = 0.0; bool converged = false; };

	// Composite Jacobi-PCG: every level and every coarse/fine connection participates
	// in one solve. A geometric multilevel preconditioner can replace Jacobi without
	// changing this conservative composite operator.
	class DeviceCompositeAmrPressureSolver
	{
	public:
		explicit DeviceCompositeAmrPressureSolver(const CompositeAmrPressureSystem& system);
		~DeviceCompositeAmrPressureSolver();
		DeviceCompositeAmrPressureSolver(const DeviceCompositeAmrPressureSolver&) = delete;
		DeviceCompositeAmrPressureSolver& operator=(const DeviceCompositeAmrPressureSolver&) = delete;
		AmrGpuSolveResult solve(Real* pressure, const Real* rhs, double tolerance, int max_iterations, bool warm_start = false);
		std::size_t bytes() const { return bytes_; }

	private:
		DeviceCompositeAmrPressureOperator op_;
		int n_ = 0;
		Real *r_ = nullptr, *z_ = nullptr, *direction_ = nullptr, *Ad_ = nullptr, *diagonal_ = nullptr;
		unsigned char* active_ = nullptr;
		int *aggregate_ = nullptr, *base_neighbors_ = nullptr, *base_extra_a_ = nullptr, *base_extra_b_ = nullptr;
		Real *base_positive_ = nullptr, *base_diagonal_ = nullptr, *base_extra_coefficient_ = nullptr, *base_rhs_ = nullptr, *base_x_ = nullptr, *base_tmp_ = nullptr;
		int base_offset_ = 0, base_bricks_ = 0, base_cells_ = 0, base_extra_count_ = 0, brick_size_ = 0;
		std::size_t bytes_ = 0;
		void apply_preconditioner(const Real* residual, Real* output);
	};

	struct DeviceCompositeAmrFluxLevelView;

	// Persistent GPU bridge from pooled brick MAC velocities to the composite pressure
	// graph. Regular faces stay implicit and structured; only coarse/fine and EB aperture
	// velocities occupy the compact `special_velocity_` array. A project() call performs
	// divergence -> RHS -> coupled pressure solve -> flux correction without bulk host
	// transfers or per-step CPU geometry work. PCG scalar reductions are host-orchestrated.
	class DeviceCompositeAmrProjection
	{
	public:
		DeviceCompositeAmrProjection(const CompositeAmrPressureSystem& system, DeviceAmrFields& fields);
		~DeviceCompositeAmrProjection();
		DeviceCompositeAmrProjection(const DeviceCompositeAmrProjection&) = delete;
		DeviceCompositeAmrProjection& operator=(const DeviceCompositeAmrProjection&) = delete;

		void clear_special_fluxes();
		void initialize_special_freestream(Real speed);
		// Refresh each 2:1 tile velocity from its collocated fine MAC face after
		// transport. EB aperture velocities remain independent compact states.
		void sync_coarse_fine_from_fields();
		// First-order same-side graph transport for split EB aperture states. Only
		// DOFs incident to EB apertures are represented; there is no full-domain CSR.
		void transport_embedded_apertures(Real dt, Real molecular_nu);
		void upload_special_fluxes(const CompositeAmrFluxes& host);
		void download_special_fluxes(CompositeAmrFluxes& host) const;
		void compute_divergence();
		void build_projection_rhs(Real rho, Real dt);
		void correct_fluxes(Real rho, Real dt);
		AmrGpuSolveResult project(Real rho, Real dt, double tolerance, int max_iterations, bool warm_start = false);
		void download_divergence(std::vector<Real>& host) const;
		void download_pressure(std::vector<Real>& host) const;
		Real* pressure() { return pressure_; }
		const Real* pressure() const { return pressure_; }
		const Real* rhs() const { return rhs_; }
		std::size_t bytes() const { return bytes_; }

	private:
		DeviceCompositeAmrPressureSolver solver_;
		DeviceAmrFields* fields_ = nullptr;
		DeviceCompositeAmrFluxLevelView* levels_ = nullptr;
		unsigned char* active_ = nullptr;
		std::uint8_t* cut_face_mask_ = nullptr;
		Real *volume_ = nullptr, *integrated_ = nullptr, *divergence_ = nullptr, *rhs_ = nullptr, *pressure_ = nullptr;
		int *first_dof_ = nullptr, *second_dof_ = nullptr;
		std::int8_t *direction_ = nullptr, *axis_ = nullptr;
		Real *open_area_ = nullptr, *centre_distance_ = nullptr, *special_velocity_ = nullptr;
		int *cf_fine_level_ = nullptr, *cf_group_ = nullptr, *cf_group_level_ = nullptr;
		std::uint64_t *cf_fine_index_ = nullptr, *cf_group_index_ = nullptr;
		std::int8_t* cf_group_axis_ = nullptr;
		Real *cf_group_area_ = nullptr, *cf_group_sum_ = nullptr;
		int *eb_node_a_ = nullptr, *eb_node_b_ = nullptr;
		int *eb_carrier_node_ = nullptr, *eb_carrier_level_ = nullptr;
		std::uint64_t* eb_carrier_index_ = nullptr;
		std::int8_t* eb_carrier_axis_ = nullptr;
		Real *eb_carrier_area_ = nullptr, *eb_node_axis_sum_ = nullptr,
			*eb_node_axis_weight_ = nullptr, *eb_transport_length_ = nullptr,
			*eb_transport_scratch_ = nullptr;
		std::vector<int> brick_counts_;
		int storage_size_ = 0, brick_size_ = 0, level_count_ = 0;
		int coarse_fine_count_ = 0, coarse_fine_group_count_ = 0, special_count_ = 0;
		int embedded_count_ = 0, embedded_node_count_ = 0, embedded_carrier_count_ = 0;
		bool outlet_ = true;
		std::size_t bytes_ = 0;
	};
}
