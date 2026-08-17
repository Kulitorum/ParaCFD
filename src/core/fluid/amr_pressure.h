// amr_pressure.h — matrix-free composite regular-region pressure operator for static AMR.
//
// Ordinary cells retain implicit Cartesian connectivity within and between same-level bricks.
// Only 2:1 coarse/fine interfaces are stored explicitly, so this is not a full-domain CSR graph.
#pragma once

#include "core/fluid/amr_fields.h"
#include "core/fluid/real.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace paracfd::core
{
	namespace detail
	{
		PARACFD_AMR_HD inline Real compact_min(Real a, Real b) { return a < b ? a : b; }
		PARACFD_AMR_HD inline Real compact_max(Real a, Real b) { return a > b ? a : b; }
		PARACFD_AMR_HD inline Real compact_abs(Real a) { return a < Real(0) ? -a : a; }
		PARACFD_AMR_HD inline Real compact_minmod(Real a, Real b)
		{
			return a * b <= Real(0) ? Real(0) : (compact_abs(a) < compact_abs(b) ? a : b);
		}
		// Squared Smagorinsky strain magnitude for a row-major velocity gradient:
		// gradient[velocity_component * 3 + derivative_axis]. This form avoids a
		// host/device sqrt in the shared manufactured tests; timestep kernels take
		// the square root after reconstructing the compact EB graph gradient.
		PARACFD_AMR_HD inline Real compact_strain_magnitude_squared(const Real* gradient)
		{
			const Real sxx = gradient[0], syy = gradient[4], szz = gradient[8];
			const Real sxy = Real(0.5) * (gradient[1] + gradient[3]);
			const Real sxz = Real(0.5) * (gradient[2] + gradient[6]);
			const Real syz = Real(0.5) * (gradient[5] + gradient[7]);
			return Real(2) * (sxx*sxx + syy*syy + szz*szz
				+ Real(2) * (sxy*sxy + sxz*sxz + syz*syz));
		}
		// Smooth-wall Spalding law expressed as friction velocity. Solving
		// Re_y = u+ y+(u+) gives the viscous limit continuously and avoids an
		// arbitrary laminar/log-layer switch. Geometry supplies actual wall distance.
		PARACFD_AMR_HD inline Real smooth_wall_friction_velocity(Real tangential_speed,
			Real wall_distance, Real molecular_nu)
		{
			if(!(tangential_speed>Real(0))||!(wall_distance>Real(0))||!(molecular_nu>Real(0)))return Real(0);
			const Real reynolds=tangential_speed*wall_distance/molecular_nu,kappa=Real(0.41),coefficient=Real(0.1185998568);Real lower=Real(0),upper=Real(80);
			for(int iteration=0;iteration<36;++iteration){const Real uplus=Real(0.5)*(lower+upper),x=kappa*uplus,x2=x*x,x3=x2*x;const Real remainder=x<Real(0.5)?x2*x2*(Real(1.0/24.0)+x*Real(1.0/120.0)+x2*Real(1.0/720.0)):exp(x)-Real(1)-x-Real(0.5)*x2-x3*Real(1.0/6.0);const Real yplus=uplus+coefficient*remainder;if(uplus*yplus<reynolds)lower=uplus;else upper=uplus;}
			return tangential_speed/compact_max(Real(0.5)*(lower+upper),Real(1e-12));
		}
		// One compact-aperture transport update. `upstream*` and `downstream` are ordered
		// along the sign of `current`; complete_chain selects MUSCL or the endpoint fallback.
		PARACFD_AMR_HD inline Real bounded_compact_transport_update(Real current,
			Real lower_node, Real upper_node, Real upstream, Real upstream2, Real downstream,
			Real courant, Real diffusion, bool complete_chain)
		{
			Real stencil_min = compact_min(current, compact_min(lower_node, upper_node));
			Real stencil_max = compact_max(current, compact_max(lower_node, upper_node));
			Real advected;
			if (complete_chain)
			{
				const Real delta = current - upstream;
				const Real slope = compact_minmod(delta, downstream - current);
				const Real upstream_slope = compact_minmod(upstream - upstream2, delta);
				advected = current - courant * delta
					- Real(0.5) * courant * (Real(1) - courant) * (slope - upstream_slope);
				stencil_min = compact_min(stencil_min, compact_min(upstream, compact_min(upstream2, downstream)));
				stencil_max = compact_max(stencil_max, compact_max(upstream, compact_max(upstream2, downstream)));
			}
			else
			{
				const Real node_upstream = current >= Real(0) ? lower_node : upper_node;
				advected = current + courant * (node_upstream - current);
			}
			const Real candidate = advected + diffusion * (lower_node + upper_node - Real(2) * current);
			return compact_max(stencil_min, compact_min(stencil_max, candidate));
		}
		PARACFD_AMR_HD inline Real bounded_compact_tangential_update(Real normal_update,
			Real current, Real lower_node, Real upper_node, const Real* advector,
			const Real* gradient, int component, Real dt)
		{
			Real tangential_rate = Real(0);
			for (int derivative_axis=0; derivative_axis<3; ++derivative_axis)
				if (derivative_axis != component)
					tangential_rate += advector[derivative_axis] * gradient[component*3+derivative_axis];
			const Real candidate = normal_update - dt*tangential_rate;
			const Real stencil_min = compact_min(normal_update, compact_min(current, compact_min(lower_node,upper_node)));
			const Real stencil_max = compact_max(normal_update, compact_max(current, compact_max(lower_node,upper_node)));
			return compact_max(stencil_min,compact_min(stencil_max,candidate));
		}
		PARACFD_AMR_HD inline Real compact_diffusive_momentum_transfer(
			Real mass_a, Real mass_b, Real velocity_a, Real velocity_b, Real fraction)
		{
			return fraction*(mass_a*mass_b/(mass_a+mass_b))*(velocity_b-velocity_a);
		}
	}

	struct DeviceCompositeAmrLevelView;
	struct AmrEmbeddedBoundaryAtlas;
	struct CoarseFinePressureConnection
	{
		int coarse_dof = -1;
		int fine_dof = -1;
		double open_area = 0.0;
		double centre_distance = 0.0;
		double normal_distance = 0.0; // centroid separation projected onto the aperture normal
		std::int8_t axis = 0;
		std::int8_t direction = 1; // +1: first DOF is lower-axis; -1: second DOF is lower-axis
		Vec3d face_centroid{};
	};
	inline double pressure_gradient_factor(const CoarseFinePressureConnection& connection)
	{
		return connection.normal_distance /
			(connection.centre_distance * connection.centre_distance);
	}
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
	struct SmoothFabricWallPatchLoad
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		Vec3d centroid{};
		Vec3d force{}; // force exerted by the fluid on this fabric patch [N]
	};
	// Host-only lookup retained for diagnostics and slice rendering. The timestep never
	// traverses these vectors: regular CUDA kernels and compact EB work lists remain
	// unchanged. Entries align with AmrEmbeddedBoundaryAtlas::levels.
	struct CompositeEbPressureSamplingMap
	{
		int level = -1;
		std::vector<int> cell_dof;     // atlas cell -> regular composite DOF, or -1 for split cells
		std::vector<int> fragment_dof; // atlas fragment -> resolved same-side composite DOF
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
		std::vector<Vec3d> centroid; // merged fluid-control-volume centroid
		// Active pressure DOF -> containing level-0 Cartesian cell DOF. This compact
		// geometric aggregation map feeds the matrix-free two-level preconditioner.
		std::vector<int> preconditioner_aggregate;
		std::vector<CoarseFinePressureConnection> coarse_fine;
		std::vector<CoarseFinePressureConnection> embedded;
		// Static CAD provenance and two-sided pressure mapping. These records are
		// consumed only when publishing loads/visualization, not by timestep kernels.
		std::vector<CompositeSurfacePressurePatch> surface_patches;
		std::vector<CompositeEbPressureSamplingMap> eb_sampling_maps;
		// One compact pressure reference for each active fluid component that cannot
		// reach the X-max Dirichlet outlet. Empty for intentionally pure-Neumann tests.
		std::vector<CompositePressureGauge> gauges;

		int dof(int level, int brick, int i, int j, int k) const;
		// Select the actual fluid-control-volume pressure at a world point. In a
		// fabric-cut cell this returns the local fragment DOF instead of the inactive
		// parent Cartesian slot, so opposite sides remain visually independent.
		int pressure_dof_at_point(const AmrEmbeddedBoundaryAtlas& embedded_boundary, Vec3d point) const;
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
	double composite_mac_carrier_volume(const CompositeAmrPressureSystem& system, int dof_a, int dof_b);
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
		// Bounded same-side graph transport for split EB aperture states. Complete
		// same-axis chains use minmod-limited MUSCL; endings and junctions use the donor
		// fallback. Only DOFs incident to EB apertures are represented, so there is no
		// full-domain CSR. A component-collocated multidimensional least-squares graph
		// reconstruction supplies all nine irregular-region velocity gradients. The
		// combined update obeys local same-side stencil bounds.
		void transport_embedded_apertures(Real dt, Real molecular_nu, Real smagorinsky_cs,
			bool apply_molecular_fabric_wall = true);
		// Production arbitrary-orientation smooth-wall model. It uses Spalding's law
		// at the actual patch-to-fragment distance and retains the opposite patch load.
		void apply_smooth_fabric_wall_model(Real dt, Real molecular_nu);
		// Throttled results/diagnostics for the wall path. Patch loads are the
		// equal-and-opposite reaction of the modelled fluid shear; momentum excludes rho.
		void download_smooth_fabric_wall_loads(Real rho,
			std::vector<SmoothFabricWallPatchLoad>& host) const;
		std::array<double, 3> fabric_wall_physical_momentum() const;
		// Validation path for the internal conservative compact momentum contract. It gathers
		// fragment-centred vector momentum from actual aperture/carrier dual masses,
		// applies one equal-and-opposite first-order donor transfer per EB aperture,
		// and scatters the increments back to the existing staggered states. Ordinary
		// regular/compact perimeter fluxes are deliberately not included yet, so this
		// is not a standalone production transport step.
		void transport_embedded_internal_momentum(Real dt);
		// Validation-only D2H reduction of the compact fragment-side momentum map.
		std::array<double, 3> embedded_transport_momentum() const;
		// Validation-only GPU reduction of the integrated mass-flux imbalance at
		// compact EB pressure nodes. This combines independent aperture states with
		// the exact ordinary MAC faces on the regular/compact perimeter.
		double max_embedded_transport_mass_imbalance() const;
		int embedded_regular_connection_count() const { return embedded_regular_connection_count_; }
		void upload_special_fluxes(const CompositeAmrFluxes& host);
		void download_special_fluxes(CompositeAmrFluxes& host) const;
		double max_abs_special_velocity() const;
		// Maximum |u_aperture| / transport_length over compact EB states [1/s].
		// This is their explicit graph-transport Courant rate; tiny-area apertures do
		// not spuriously constrain dt merely because their point velocity is large.
		double max_embedded_cfl_rate() const;
		int embedded_high_order_stencil_count() const { return embedded_high_order_stencil_count_; }
		int embedded_least_squares_full_rank_count() const { return embedded_least_squares_full_rank_count_; }
		int fabric_wall_node_count() const { return fabric_wall_node_count_; }
		// Validation-only D2H diagnostic. Production timesteps never call this.
		void download_embedded_node_gradients(std::vector<Real>& gradients,
			std::vector<std::uint8_t>* component_rank = nullptr) const;
		void compute_divergence();
		void build_projection_rhs(Real rho, Real dt);
		void correct_fluxes(Real rho, Real dt);
		AmrGpuSolveResult project(Real rho, Real dt, double tolerance, int max_iterations, bool warm_start = false);
		void download_divergence(std::vector<Real>& host) const;
		void download_pressure(std::vector<Real>& host) const;
		Real* pressure() { return pressure_; }
		const Real* pressure() const { return pressure_; }
		const Real* rhs() const { return rhs_; }
		// First segment of the compact special-flux array: one +axis velocity for
		// every fine-owned 2:1 pressure aperture. Exposed read-only so conservative
		// momentum transport can use exactly the projected mass flux on the GPU.
		const Real* coarse_fine_velocity_device() const { return special_velocity_; }
		int coarse_fine_velocity_count() const { return coarse_fine_count_; }
		Real* embedded_velocity_device() { return special_velocity_ + coarse_fine_count_; }
		const Real* embedded_velocity_device() const { return special_velocity_ + coarse_fine_count_; }
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
		Real *open_area_ = nullptr, *centre_distance_ = nullptr, *pressure_gradient_factor_ = nullptr,
			*special_velocity_ = nullptr, *max_abs_scratch_ = nullptr;
		int *cf_fine_level_ = nullptr, *cf_group_ = nullptr, *cf_group_level_ = nullptr;
		std::uint64_t *cf_fine_index_ = nullptr, *cf_group_index_ = nullptr;
		std::int8_t* cf_group_axis_ = nullptr;
		Real *cf_group_area_ = nullptr, *cf_group_sum_ = nullptr;
		int *eb_node_a_ = nullptr, *eb_node_b_ = nullptr;
		int *eb_negative_edge_ = nullptr, *eb_positive_edge_ = nullptr;
		int *eb_carrier_node_ = nullptr, *eb_carrier_level_ = nullptr;
		std::uint64_t* eb_carrier_index_ = nullptr;
		std::int8_t* eb_carrier_axis_ = nullptr;
		int *eb_regular_lower_node_ = nullptr, *eb_regular_upper_node_ = nullptr,
			*eb_regular_level_ = nullptr;
		std::uint64_t* eb_regular_index_ = nullptr;
		std::int8_t* eb_regular_axis_ = nullptr;
		Real *eb_carrier_area_ = nullptr, *eb_carrier_mass_ = nullptr,
			*eb_carrier_unmapped_mass_ = nullptr, *eb_node_axis_sum_ = nullptr,
			*eb_node_axis_weight_ = nullptr, *eb_node_gradient_sum_ = nullptr,
			*eb_node_gradient_inverse_ = nullptr, *eb_edge_ls_displacement_ = nullptr,
			*eb_edge_ls_weight_ = nullptr, *eb_transport_length_ = nullptr,
			*eb_transport_scratch_ = nullptr, *eb_diffusion_rate_ = nullptr,
			*eb_node_diffusion_sum_ = nullptr, *eb_node_neighbor_count_ = nullptr,
			*eb_diffusion_scratch_ = nullptr, *eb_regular_area_ = nullptr,
			*eb_node_flux_balance_ = nullptr;
		int *wall_edge_node_a_ = nullptr, *wall_edge_node_b_ = nullptr;
		int* wall_patch_node_ = nullptr;
		int *wall_carrier_node_ = nullptr, *wall_carrier_level_ = nullptr;
		std::uint64_t* wall_carrier_index_ = nullptr;
		std::int8_t* wall_carrier_axis_ = nullptr;
		int* wall_unique_carrier_level_ = nullptr;
		std::uint64_t* wall_unique_carrier_index_ = nullptr;
		std::int8_t* wall_unique_carrier_axis_ = nullptr;
		Real *wall_carrier_mass_ = nullptr, *wall_node_rate_ = nullptr,
			*wall_node_axis_sum_ = nullptr, *wall_node_axis_weight_ = nullptr,
			*wall_node_velocity_delta_ = nullptr, *wall_patch_normal_ = nullptr,
			*wall_patch_area_ = nullptr,
			*wall_patch_distance_ = nullptr, *wall_patch_force_per_density_ = nullptr,
			*wall_unique_carrier_mass_ = nullptr;
		double* wall_momentum_scratch_ = nullptr;
		std::vector<std::uint32_t> wall_patch_source_triangle_id_, wall_patch_source_face_id_;
		std::vector<Vec3d> wall_patch_centroid_;
		std::vector<int> brick_counts_;
		std::vector<std::uint8_t> embedded_gradient_rank_;
		int storage_size_ = 0, brick_size_ = 0, level_count_ = 0;
		int coarse_fine_count_ = 0, coarse_fine_group_count_ = 0, special_count_ = 0;
		int embedded_count_ = 0, embedded_node_count_ = 0, embedded_carrier_count_ = 0,
			embedded_regular_connection_count_ = 0;
		int embedded_high_order_stencil_count_ = 0;
		int embedded_least_squares_full_rank_count_ = 0;
		int fabric_wall_node_count_ = 0, fabric_wall_carrier_count_ = 0,
			fabric_wall_patch_count_ = 0, fabric_wall_unique_carrier_count_ = 0;
		bool outlet_ = true;
		std::size_t bytes_ = 0;
	};
}
