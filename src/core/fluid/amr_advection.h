#pragma once

#include "core/fluid/amr_fields.h"
#include "core/geometry/triangle_bvh.h"

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace paracfd::core
{
	struct CompositeAmrPressureSystem;
	class DeviceAmrMacFaceMap;

	struct AmrMacFaceAddress
	{
		int level = -1;
		int brick = -1;
		int component = -1;
		int i = -1, j = -1, k = -1;
	};
	// Same-level bricks duplicate their shared component-normal MAC face. Compact
	// topology owns that physical face once, on the negative-coordinate brick.
	AmrMacFaceAddress canonical_amr_mac_face_address(const AmrHierarchy& hierarchy,
		AmrMacFaceAddress address);
	// Volume owned by the regular portion of a staggered face control volume: one
	// half of each active same-level pressure cell adjacent along the component axis.
	// A normal 2:1 interface tile adds/replaces the missing covered-side contribution
	// separately; physical-boundary faces intentionally own one half cell.
	double regular_mac_dual_volume(const AmrHierarchy& hierarchy,
		const AmrMacFaceAddress& address);

	// One normal-velocity tile at a 2:1 interface. The fine interface face owns the
	// transported state; the coincident coarse face is an area-mean alias. Transport
	// links from the coarse/fine interior faces meet at this tile, whose dual volume
	// spans half a coarse cell and half a fine cell. This keeps the four fine fluxes
	// explicit instead of evolving an independent coarse interface value.
	struct NormalMomentumInterfaceTile
	{
		AmrMacFaceAddress coarse_interface;
		AmrMacFaceAddress fine_interface;
		AmrMacFaceAddress coarse_interior;
		AmrMacFaceAddress fine_interior;
		double area = 0.0;
		double dual_volume = 0.0;
		std::int8_t direction = 1;
	};

	std::vector<NormalMomentumInterfaceTile> build_normal_momentum_interface_tiles(
		const CompositeAmrPressureSystem& system);

	struct TangentialMomentumInterfaceConnection
	{
		struct FluxSource
		{
			int coarse_fine_connection = -1;
			double overlap_area = 0.0;
		};
		AmrMacFaceAddress coarse;
		AmrMacFaceAddress fine;
		double open_area = 0.0;
		double centre_distance = 0.0;
		std::int8_t interface_axis = 0;
		std::int8_t component = 0;
		std::int8_t direction = 1;
		std::vector<FluxSource> flux_sources;
	};

	// Intersect the staggered component-dual rectangles on both sides of every 2:1
	// interface. This is not assumed to be four links: face-centred tangential lattices
	// meet on grid lines, so exact overlap clipping and same-level face canonicalization
	// are required to avoid double counting at child/brick boundaries.
	std::vector<TangentialMomentumInterfaceConnection> build_tangential_momentum_interface_connections(
		const CompositeAmrPressureSystem& system);

	struct AmrMacFaceDirectionalLink
	{
		AmrMacFaceAddress face;
		std::int8_t transport_axis = 0;
		std::int8_t direction = 1;
	};
	std::vector<AmrMacFaceDirectionalLink> build_momentum_interface_replaced_links(
		const CompositeAmrPressureSystem& system);

	struct CompositeAmrMomentumInterfaceConnection
	{
		int lower_node = -1;
		int upper_node = -1;
		double open_area = 0.0;
		std::int8_t transport_axis = 0;
		std::int8_t component = 0;
	};

	struct CompositeAmrMomentumCoarseAlias
	{
		AmrMacFaceAddress coarse_face;
		int fine_owned_node = -1;
		double area = 0.0;
	};

	struct CompositeAmrMomentumMassFluxSource
	{
		int momentum_connection = -1;
		int coarse_fine_connection = -1;
		double overlap_area = 0.0;
	};

	// One compact ownership graph for every 2:1 momentum interface. Normal fine
	// tiles and tangential coarse/fine states share the same canonical node table,
	// including faces at intersections of two or three refinement boundaries.
	// Connections are oriented lower -> upper along transport_axis. Coarse normal
	// faces are aliases only and are reconstructed from fine-owned tile states.
	struct CompositeAmrMomentumInterfaceTopology
	{
		std::vector<AmrMacFaceAddress> nodes;
		std::vector<double> dual_volume;
		std::vector<std::uint8_t> interface_axis_mask;
		std::vector<CompositeAmrMomentumInterfaceConnection> connections;
		std::vector<CompositeAmrMomentumMassFluxSource> mass_flux_sources;
		std::vector<CompositeAmrMomentumCoarseAlias> coarse_aliases;
	};

	CompositeAmrMomentumInterfaceTopology build_composite_amr_momentum_interface_topology(
		const CompositeAmrPressureSystem& system);
	std::vector<double> composite_amr_momentum_connection_velocities_cpu(
		const CompositeAmrMomentumInterfaceTopology& topology,
		const std::vector<double>& node_velocity,
		const std::vector<double>& coarse_fine_velocity);

	// Persistent GPU indirection only for compact interface/EB work. Ordinary regular
	// cells remain on structured brick kernels. Addresses must be canonical and unique,
	// making scatter race-free; the field allocations themselves stay owned by
	// DeviceAmrFields.
	class DeviceAmrMacFaceMap
	{
	public:
		DeviceAmrMacFaceMap(DeviceAmrFields& fields,
			const std::vector<AmrMacFaceAddress>& addresses);
		~DeviceAmrMacFaceMap();
		DeviceAmrMacFaceMap(const DeviceAmrMacFaceMap&) = delete;
		DeviceAmrMacFaceMap& operator=(const DeviceAmrMacFaceMap&) = delete;

		void gather(Real* compact) const;
		void scatter(const Real* compact) const;
		int size() const { return size_; }
		std::size_t bytes() const { return bytes_; }

	private:
		DeviceAmrFieldLevelView* levels_ = nullptr;
		int* level_ = nullptr;
		std::uint64_t* index_ = nullptr;
		std::int8_t* component_ = nullptr;
		int size_ = 0, level_count_ = 0;
		std::size_t bytes_ = 0;
	};

	// Geometry-independent conservative transport contract used by the composite
	// momentum path. A node is one velocity-component dual control volume (a regular
	// MAC face, a 2:1 tile, or a compact EB aperture state). Every connection is
	// oriented a -> b and contributes one shared upwind momentum flux to both nodes.
	// Fabric is impermeable by topology: no connection is emitted through a patch.
	struct PairwiseMomentumConnection
	{
		int a = -1;
		int b = -1;
		double open_area = 0.0;
		double normal_velocity = 0.0; // signed a -> b
	};

	void conservative_pairwise_momentum_cpu(const std::vector<double>& dual_volume,
		const std::vector<PairwiseMomentumConnection>& connections, double dt,
		std::vector<double>& velocity_x, std::vector<double>& velocity_y,
		std::vector<double>& velocity_z);
	void conservative_pairwise_scalar_cpu(const std::vector<double>& dual_volume,
		const std::vector<PairwiseMomentumConnection>& connections, double dt,
		std::vector<double>& velocity);

	// Ordinary same-level pressure-face connection incident to at least one fluid DOF
	// that also participates in the compact EB aperture graph. These records define the
	// missing regular/compact perimeter ownership without introducing a full-domain graph.
	struct CompositeEbMomentumRegularConnection
	{
		int lower_dof = -1;
		int upper_dof = -1;
		AmrMacFaceAddress face;
		double open_area = 0.0;
		double carrier_volume = 0.0;
		bool lower_embedded_node = false;
		bool upper_embedded_node = false;
	};
	std::vector<CompositeEbMomentumRegularConnection>
		build_composite_eb_momentum_regular_connections(const CompositeAmrPressureSystem& system);

	// Cell/fragment-centred conservative momentum state. Pressure projection still
	// owns staggered regular-face and compact aperture mass fluxes; those fluxes
	// transport this vector state between the actual fluid control volumes. This
	// avoids inventing a staggered dual-volume mortar at the regular/fragment edge.
	// The first implementation is deliberately one uniform AMR level: ordinary
	// connections remain an implicit brick kernel and only EB aperture/perimeter
	// connections are materialized.
	struct CompositeCellMomentumState
	{
		std::vector<Real> x, y, z;
	};

	void conservative_composite_cell_momentum_cpu(
		const CompositeAmrPressureSystem& system, const AmrHostFields& fields,
		const std::vector<double>& embedded_velocity, double dt,
		CompositeCellMomentumState& state, bool external_aero = false,
		double freestream_speed = 0.0);

	class DeviceCompositeCellMomentumTransport
	{
	public:
		DeviceCompositeCellMomentumTransport(const CompositeAmrPressureSystem& system,
			DeviceAmrFields& fields);
		~DeviceCompositeCellMomentumTransport();
		DeviceCompositeCellMomentumTransport(const DeviceCompositeCellMomentumTransport&) = delete;
		DeviceCompositeCellMomentumTransport& operator=(const DeviceCompositeCellMomentumTransport&) = delete;

		void upload_state(const CompositeCellMomentumState& state);
		void download_state(CompositeCellMomentumState& state) const;
		void step(const Real* embedded_velocity, Real dt, bool external_aero = false,
			Real freestream_speed = Real(0));
		std::array<double, 3> momentum() const;
		int regular_compact_connection_count() const { return regular_count_; }
		int embedded_connection_count() const { return embedded_count_; }
		std::size_t bytes() const;

	private:
		const CompositeAmrPressureSystem* system_ = nullptr;
		DeviceAmrFields* fields_ = nullptr;
		std::unique_ptr<DeviceAmrMacFaceMap> regular_face_map_;
		Real *x_ = nullptr, *y_ = nullptr, *z_ = nullptr;
		Real *delta_x_ = nullptr, *delta_y_ = nullptr, *delta_z_ = nullptr;
		Real *volume_ = nullptr, *regular_velocity_ = nullptr;
		unsigned char *active_ = nullptr, *cut_face_mask_ = nullptr,
			*compact_plus_mask_ = nullptr;
		int *embedded_a_ = nullptr, *embedded_b_ = nullptr;
		Real* embedded_area_ = nullptr;
		int *regular_a_ = nullptr, *regular_b_ = nullptr;
		Real* regular_area_ = nullptr;
		std::vector<double> volume_host_;
		std::vector<unsigned char> active_host_;
		int storage_size_ = 0, base_cell_count_ = 0, level_offset_ = 0;
		int embedded_count_ = 0, regular_count_ = 0;
		std::size_t bytes_ = 0;
	};

	// GPU twin of conservative_pairwise_momentum_cpu. Static volumes and topology are
	// uploaded once; only the per-step signed connection velocities are supplied by the
	// caller. Integrated momentum increments are accumulated pairwise in persistent SoA
	// scratch, so the two endpoints receive exactly opposite transfers.
	class DevicePairwiseMomentumTransport
	{
	public:
		DevicePairwiseMomentumTransport(const std::vector<double>& dual_volume,
			const std::vector<PairwiseMomentumConnection>& connections);
		~DevicePairwiseMomentumTransport();
		DevicePairwiseMomentumTransport(const DevicePairwiseMomentumTransport&) = delete;
		DevicePairwiseMomentumTransport& operator=(const DevicePairwiseMomentumTransport&) = delete;

		void step(Real* velocity_x, Real* velocity_y, Real* velocity_z,
			const Real* connection_normal_velocity, Real dt);
		void step_scalar(Real* velocity, const Real* connection_normal_velocity, Real dt);
		int node_count() const { return node_count_; }
		int connection_count() const { return connection_count_; }
		std::size_t bytes() const { return bytes_; }

	private:
		int *a_ = nullptr, *b_ = nullptr;
		Real *volume_ = nullptr, *area_ = nullptr;
		Real *delta_x_ = nullptr, *delta_y_ = nullptr, *delta_z_ = nullptr;
		int node_count_ = 0, connection_count_ = 0;
		std::size_t bytes_ = 0;
	};

	// Staged GPU finite-volume update for the compact EB graph plus its exact
	// ordinary-face perimeter. Fragment-centred vectors are reconstructed from
	// unique physical MAC carriers and embedded aperture dual masses. Every
	// aperture/perimeter transfer is applied with equal and opposite endpoint
	// momentum increments, then scattered back to the shared staggered states.
	// The surrounding structured fluxes are not part of this class yet, so it is
	// a validation/coupling layer rather than the production advection driver.
	class DeviceCompositeEbMomentumTransport
	{
	public:
		DeviceCompositeEbMomentumTransport(const CompositeAmrPressureSystem& system,
			DeviceAmrFields& fields);
		~DeviceCompositeEbMomentumTransport();
		DeviceCompositeEbMomentumTransport(const DeviceCompositeEbMomentumTransport&) = delete;
		DeviceCompositeEbMomentumTransport& operator=(const DeviceCompositeEbMomentumTransport&) = delete;

		void step(Real* embedded_velocity, Real dt);
		// Validation-only D2H total over each unique physical MAC carrier and EB aperture.
		std::array<double, 3> momentum(const Real* embedded_velocity) const;
		int node_count() const { return node_count_; }
		int embedded_connection_count() const { return embedded_count_; }
		int regular_connection_count() const { return regular_count_; }
		int carrier_count() const { return carrier_count_; }
		int same_level_alias_count() const { return alias_count_; }
		std::size_t bytes() const;

	private:
		std::unique_ptr<DeviceAmrMacFaceMap> carrier_map_;
		std::unique_ptr<DeviceAmrMacFaceMap> alias_map_;
		Real *carrier_state_ = nullptr, *carrier_delta_ = nullptr;
		int* alias_source_ = nullptr;
		Real* alias_value_ = nullptr;
		int *incidence_node_ = nullptr, *incidence_carrier_ = nullptr;
		std::int8_t* incidence_component_ = nullptr;
		Real* incidence_half_mass_ = nullptr;
		int *embedded_a_ = nullptr, *embedded_b_ = nullptr;
		std::int8_t* embedded_axis_ = nullptr;
		Real *embedded_area_ = nullptr, *embedded_mass_ = nullptr;
		int *regular_a_ = nullptr, *regular_b_ = nullptr, *regular_carrier_ = nullptr;
		Real* regular_area_ = nullptr;
		Real *node_sum_ = nullptr, *node_weight_ = nullptr, *node_delta_ = nullptr;
		std::vector<double> carrier_mass_host_, embedded_mass_host_;
		std::vector<std::int8_t> carrier_component_host_, embedded_axis_host_;
		int node_count_ = 0, incidence_count_ = 0, carrier_count_ = 0, alias_count_ = 0;
		int embedded_count_ = 0, regular_count_ = 0;
		std::size_t bytes_ = 0;
	};

	// Persistent GPU driver for the compact 2:1 interface graph. It gathers each
	// canonical pooled MAC state once, advances one pairwise scalar momentum flux per
	// connection, scatters the owned states, then rebuilds coarse normal-face aliases
	// as aperture-area means. The caller supplies one +transport-axis velocity for
	// every connection; deriving those velocities from composite mass flux is the next
	// coupling layer.
	class DeviceCompositeAmrMomentumInterfaceTransport
	{
	public:
		DeviceCompositeAmrMomentumInterfaceTransport(const CompositeAmrPressureSystem& system,
			DeviceAmrFields& fields);
		~DeviceCompositeAmrMomentumInterfaceTransport();
		DeviceCompositeAmrMomentumInterfaceTransport(const DeviceCompositeAmrMomentumInterfaceTransport&) = delete;
		DeviceCompositeAmrMomentumInterfaceTransport& operator=(const DeviceCompositeAmrMomentumInterfaceTransport&) = delete;

		void step(const Real* connection_normal_velocity, Real dt);
		// Derive normal-component advectors from the endpoint MAC states and
		// tangential advectors from the conservative coarse/fine pressure flux tiles.
		void step_from_composite_flux(const Real* coarse_fine_velocity, Real dt);
		int node_count() const { return node_count_; }
		int connection_count() const { return connection_count_; }
		int coarse_alias_group_count() const { return alias_group_count_; }
		int mass_flux_source_count() const { return mass_flux_source_count_; }
		std::size_t bytes() const;

	private:
		std::unique_ptr<DeviceAmrMacFaceMap> node_map_;
		std::unique_ptr<DeviceAmrMacFaceMap> alias_map_;
		std::unique_ptr<DevicePairwiseMomentumTransport> transport_;
		Real *node_state_ = nullptr, *connection_normal_velocity_ = nullptr;
		Real *alias_sum_ = nullptr, *alias_value_ = nullptr;
		int *alias_node_ = nullptr, *alias_group_ = nullptr;
		Real *alias_area_ = nullptr, *alias_group_area_ = nullptr;
		int *connection_lower_ = nullptr, *connection_upper_ = nullptr;
		std::int8_t *connection_axis_ = nullptr, *connection_component_ = nullptr;
		Real* connection_area_ = nullptr;
		int *mass_flux_connection_ = nullptr, *mass_flux_source_ = nullptr;
		Real* mass_flux_area_ = nullptr;
		int node_count_ = 0, connection_count_ = 0;
		int alias_count_ = 0, alias_group_count_ = 0;
		int mass_flux_source_count_ = 0;
		std::size_t bytes_ = 0;
	};

	// Production AMR advection layer. Backtraces locate the finest active brick
	// through the integer-coordinate GPU hash. Static preprocessing marks faces within
	// `protection_cells * h` of fabric and stores their six Cartesian same-side links.
	// Those faces use bounded, minmod-limited same-side transport (with a first-order
	// fallback wherever the required links are unavailable) plus link-restricted molecular
	// diffusion; far-field faces retain bounded RK2/MacCormack transport and LES. This
	// permits tangential transport near a sheet without per-step triangle traversal or
	// opposite-side stencil sampling.
	// Trilinear and diffusion stencils resolve same-level brick crossings through the
	// GPU coordinate hash. At 2:1 transitions, staggered face and cell values use a
	// one-sided linear prolongation rather than brick-local clamping; the reverse
	// MacCormack correction falls back to the bounded forward RK2 value when its trace
	// changes lattice. This is linearly consistent, while the semi-Lagrangian momentum
	// update itself is not yet a globally conservative finite-volume/refluxed scheme.
	class DeviceAmrAdvection
	{
	public:
		DeviceAmrAdvection(const AmrHierarchy& hierarchy, const TriangleBvh& fabric,
			double protection_cells = 2.5,
			const CompositeAmrPressureSystem* composite_system = nullptr);
		~DeviceAmrAdvection();
		DeviceAmrAdvection(const DeviceAmrAdvection&) = delete;
		DeviceAmrAdvection& operator=(const DeviceAmrAdvection&) = delete;

		void advect(DeviceAmrFields& fields, Real dt);
		// Incremental conservative-momentum foundation. This finite-volume MAC update
		// uses one shared upwind flux per same-level lattice link and is currently
		// restricted to one uniform level. Fabric-band links can suppress a flux, but
		// compact EB aperture coupling and 2:1 flux registers are not yet included. It
		// remains separate from production advect() until both pieces are complete.
		void advect_conservative_uniform(DeviceAmrFields& fields, Real dt);
		// One-level validation path with the production external-aero momentum BC:
		// prescribed +X inflow donor state, convective X-max outflow, and zero normal
		// mass flux at the Y/Z free-slip boundaries. Boundary MAC dual volumes and
		// edge/corner face areas use their exact half-width geometry.
		void advect_conservative_uniform_external(DeviceAmrFields& fields, Real dt,
			Real freestream_speed);
		void diffuse_smagorinsky(DeviceAmrFields& fields, Real molecular_nu, Real cs, Real dt);
		std::size_t protected_face_count() const { return protected_faces_; }
		std::size_t active_face_count() const { return active_faces_; }
		std::size_t replaced_interface_link_count() const { return replaced_interface_links_; }
		std::size_t bytes() const { return bytes_ + locator_.bytes(); }

	private:
		struct Level
		{
			BrickFieldLayout layout;
			int brick_count = 0;
			Real *u = nullptr, *v = nullptr, *w = nullptr;
			Real *forward_u = nullptr, *forward_v = nullptr, *forward_w = nullptr;
			// Bits 0..5 are -/+ xyz links; bit 6 marks the near-fabric band.
			unsigned char *links_u = nullptr, *links_v = nullptr, *links_w = nullptr;
		};
		DeviceAmrLocator locator_;
		DeviceAmrFieldLevelView *device_views_ = nullptr, *device_forward_views_ = nullptr;
		std::vector<Level> levels_;
		std::size_t protected_faces_ = 0, active_faces_ = 0, replaced_interface_links_ = 0, bytes_ = 0;
	};
}
