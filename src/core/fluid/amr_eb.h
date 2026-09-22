// amr_eb.h — bridge between block-pooled AMR fields and zero-thickness EB topology.
//
// The first implementation deliberately targets one complete uniform level. Geometry is
// preprocessed in level-global cell coordinates, so AMR brick boundaries are never mistaken
// for physical boundaries and a membrane lying on a brick interface remains impermeable.
#pragma once

#include "core/fluid/amr_fields.h"
#include "core/fluid/eb_pressure.h"

#include <cstddef>
#include <vector>

namespace paracfd::core
{
	struct OneLevelEmbeddedBoundary
	{
		EmbeddedBoundary topology;
		BrickFieldLayout field_layout;
		// Global Cartesian cell -> location in level.p[brick][cell-with-halos].
		std::vector<std::size_t> cell_field_index;

		bool ready_for_flow() const { return topology.ready_for_flow(); }
	};

	OneLevelEmbeddedBoundary build_one_level_embedded_boundary(const AmrHierarchy& hierarchy,
		const TriMesh& mesh, const TriangleBvh& bvh, const EmbeddedBoundaryBuildOptions& options = {});

	struct AmrEbLevelAtlas
	{
		int level = -1;
		Int3 minimum_brick_coord{};
		Int3 maximum_brick_coord{};
		EmbeddedBoundary topology;
		std::vector<unsigned char> owned_cell; // active cells whose finest owner is this level
		std::size_t quality_agglomerations = 0; // same-fluid merges restoring usable aperture geometry
		std::size_t static_subcell_components = 0; // isolated one-root states with no positive-area aperture
		std::size_t static_subcell_aggregates = 0;
		double static_subcell_volume = 0;
		std::size_t face_state_retained_small_roots = 0;
		double face_state_retained_small_volume = 0;
		double minimum_face_state_retained_volume_fraction = 1;
		double maximum_aggregate_span_cells = 0; // constituent-centroid span / h after all merges
	};

	struct AmrEmbeddedBoundaryAtlas
	{
		std::vector<AmrEbLevelAtlas> levels;
		// Entire active bricks inside material need a mask even when their level
		// has no surface atlas (e.g. a coarse brick enclosed by a large solid).
		std::vector<std::vector<unsigned char>> solid_bricks;
		std::size_t unresolved_count() const;
		std::size_t owned_unresolved_count() const;
		// Covered cells are retained as a one-cell topology halo. A rejection there
		// invalidates the face/fragment contract seen by its owned neighbour even though
		// it owns no pressure DOF on this level. Ownership is a storage policy, not a
		// licence to ignore a geometry failure.
		bool ready_for_flow() const { return unresolved_count() == 0; }
	};

	// Build one sparse rectangular topology atlas around the active EB bricks of each
	// level, plus a one-cell halo. Geometry never sees an individual brick boundary as
	// a physical domain boundary; cells covered by a finer owner are retained only as
	// topology halo and excluded through `owned_cell`.
	AmrEmbeddedBoundaryAtlas build_amr_embedded_boundary_atlas(const AmrHierarchy& hierarchy,
		const TriMesh& mesh, const TriangleBvh& bvh, const EmbeddedBoundaryBuildOptions& options = {});

	// Host reference bridge. Irregular pressure entries remain in `pressure`; regular entries
	// are gathered/scattered from the level pool without changing the compact EB DOF numbering.
	void gather_one_level_pressure(const OneLevelEmbeddedBoundary& composite, const EbPressureSystem& system,
		const AmrHostLevelFields& level, std::vector<Real>& pressure);
	void scatter_one_level_pressure(const OneLevelEmbeddedBoundary& composite, const EbPressureSystem& system,
		const std::vector<Real>& pressure, AmrHostLevelFields& level);

	// Device bridge performs the same operation entirely on the GPU. This is an incremental
	// compatibility layer; no D2H/H2D field transfer is introduced into a timestep.
	class DeviceOneLevelPressureMap
	{
	public:
		DeviceOneLevelPressureMap(const OneLevelEmbeddedBoundary& composite, const EbPressureSystem& system);
		~DeviceOneLevelPressureMap();
		DeviceOneLevelPressureMap(const DeviceOneLevelPressureMap&) = delete;
		DeviceOneLevelPressureMap& operator=(const DeviceOneLevelPressureMap&) = delete;
		void gather(const Real* level_pressure, Real* eb_pressure) const;
		void scatter(const Real* eb_pressure, Real* level_pressure) const;
		std::size_t bytes() const { return bytes_; }

	private:
		std::uint64_t* field_index_ = nullptr;
		int storage_size_ = 0;
		std::size_t bytes_ = 0;
	};
}
