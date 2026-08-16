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
		std::vector<CoarseFinePressureConnection> coarse_fine;
		std::vector<CoarseFinePressureConnection> embedded;

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
		unsigned char *active_ = nullptr;
		std::uint8_t *cut_face_mask_ = nullptr;
		int level_count_ = 0, connection_count_ = 0, storage_size_ = 0, brick_size_ = 0;
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
		std::size_t bytes_ = 0;
	};
}
