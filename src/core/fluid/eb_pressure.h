// eb_pressure.h — matrix-free finite-volume pressure operator for zero-thickness EB topology.
//
// Regular-to-regular faces are never stored: CPU/GPU apply them as a Cartesian stencil. Only
// split apertures use compact records. Rows are integrated flux balances (coefficient A_open/d),
// while reported divergence is divided by the actual control-volume volume.
#pragma once

#include "core/fluid/real.h"
#include "core/geometry/embedded_boundary.h"

#include <vector>

namespace paracfd::core
{
	struct EbPressureSystem
	{
		const EmbeddedBoundary* eb = nullptr;
		int storage_size = 0; // cell_count + fragment_count; split-cell slots are inactive
		std::vector<int> cell_dof;
		std::vector<int> fragment_dof; // merged fragments alias their conservative target
		std::vector<double> volume;     // aggregated actual fluid volume per active DOF
		std::vector<Vec3d> centroid;
		std::vector<unsigned char> active;
		bool pressure_outlet_xmax = true;

		int dof(FragmentRef ref) const;
		void apply_cpu(const std::vector<double>& p, std::vector<double>& Ap) const;
		void diagonal_cpu(std::vector<double>& diagonal) const;
		double volume_weighted_mean(const std::vector<double>& x) const;
	};

	EbPressureSystem build_eb_pressure_system(const EmbeddedBoundary& eb, bool pressure_outlet_xmax = true);

	struct EbFaceFluxes
	{
		// Normal velocity [m/s], positive along +axis. Faces touching split cells are ignored here
		// and replaced by `aperture_velocity` entries.
		std::vector<double> x, y, z;
		std::vector<double> aperture_velocity;
	};

	EbFaceFluxes make_zero_fluxes(const EmbeddedBoundary& eb);
	void eb_divergence_cpu(const EbPressureSystem& system, const EbFaceFluxes& flux, std::vector<double>& divergence);
	void eb_projection_rhs_cpu(const EbPressureSystem& system, const EbFaceFluxes& flux, double rho, double dt, std::vector<double>& rhs);
	void eb_correct_fluxes_cpu(const EbPressureSystem& system, const std::vector<double>& pressure, double rho, double dt, EbFaceFluxes& flux);

	struct EbCpuSolveResult { int iterations=0; double relative_residual=0; bool converged=false; };
	EbCpuSolveResult solve_eb_pressure_cpu(const EbPressureSystem& system, const std::vector<double>& rhs,
		std::vector<double>& pressure, double tolerance=1e-10, int max_iterations=1000);

	// GPU-resident FP32 production operator; reductions/validation remain double on the host.
	class DeviceEbPressureOperator
	{
	public:
		explicit DeviceEbPressureOperator(const EbPressureSystem& system);
		~DeviceEbPressureOperator();
		DeviceEbPressureOperator(const DeviceEbPressureOperator&)=delete;
		DeviceEbPressureOperator& operator=(const DeviceEbPressureOperator&)=delete;
		void apply(const Real* pressure, Real* output) const;
		int storage_size() const { return storage_size_; }
		std::size_t bytes() const { return bytes_; }

	private:
		UniformEbGrid grid_;
		int storage_size_=0, aperture_count_=0;
		int *cell_dof_=nullptr,*fragment_dof_=nullptr;
		std::uint8_t* cut_face_mask_=nullptr;
		FragmentRef *aperture_a_=nullptr,*aperture_b_=nullptr;
		Real *aperture_area_=nullptr,*aperture_distance_=nullptr;
		std::size_t bytes_=0;
		bool outlet_=true;
	};

	struct EbGpuSolveResult
	{
		int iterations = 0;
		double relative_residual = 0.0;
		bool converged = false;
	};

	// GPU-resident Jacobi-preconditioned CG for the one-level composite EB operator. Field vectors
	// stay on device; only scalar convergence reductions return to the host, matching the retained
	// MGPCG control pattern. Static disconnected all-Neumann components must be gauge-pinned during
	// preprocessing; the initial external-aero path uses the X-max pressure outlet.
	class DeviceEbPressureSolver
	{
	public:
		explicit DeviceEbPressureSolver(const EbPressureSystem& system);
		~DeviceEbPressureSolver();
		DeviceEbPressureSolver(const DeviceEbPressureSolver&) = delete;
		DeviceEbPressureSolver& operator=(const DeviceEbPressureSolver&) = delete;
		EbGpuSolveResult solve(Real* pressure, const Real* rhs, double tolerance, int max_iterations, bool warm_start = false);
		std::size_t bytes() const { return bytes_; }

	private:
		DeviceEbPressureOperator op_;
		int n_ = 0;
		Real *r_ = nullptr, *z_ = nullptr, *direction_ = nullptr, *Ad_ = nullptr, *diagonal_ = nullptr;
		unsigned char* active_ = nullptr;
		std::size_t bytes_ = 0;
	};
}
