// amr_grid.h — static block-structured Cartesian AMR metadata.
//
// Refinement is always 2:1. Cells are uniform inside a brick; no cell/octree pointers are used.
// Point location descends compact per-level integer-coordinate hash tables and returns the finest
// active brick. CUDA kernels receive the same flat lookup representation (see amr_fields.cu).
#pragma once

#include "core/geometry/triangle_bvh.h"
#include "core/paraglider_config.h"

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace paracfd::core
{
	struct Int3
	{
		int x = 0, y = 0, z = 0;
		bool operator==(const Int3& o) const { return x == o.x && y == o.y && z == o.z; }
	};

	enum BrickFlags : std::uint32_t
	{
		BRICK_XMIN = 1u << 0, BRICK_XMAX = 1u << 1, BRICK_YMIN = 1u << 2,
		BRICK_YMAX = 1u << 3, BRICK_ZMIN = 1u << 4, BRICK_ZMAX = 1u << 5,
		BRICK_EMBEDDED_BOUNDARY = 1u << 6, BRICK_COVERED = 1u << 7
	};

	struct BrickMetadata
	{
		int level = 0;
		Int3 coord{};
		Vec3d origin{};
		double h = 0.0;
		std::array<int, 6> same_level_neighbor{{-1, -1, -1, -1, -1, -1}};
		int parent = -1;
		std::array<int, 8> children{{-1, -1, -1, -1, -1, -1, -1, -1}};
		std::uint32_t flags = 0;

		bool active() const { return (flags & BRICK_COVERED) == 0; }
		bool embedded_boundary() const { return (flags & BRICK_EMBEDDED_BOUNDARY) != 0; }
	};

	struct BrickLookupEntry
	{
		Int3 coord{};
		int brick_id = -1; // -1 is an empty hash slot
	};

	struct AmrLevel
	{
		int level = 0;
		double h = 0.0;
		std::vector<BrickMetadata> bricks;
		std::vector<BrickLookupEntry> lookup; // power-of-two open-addressed table
		std::uint32_t lookup_mask = 0;
		std::size_t active_bricks = 0;
	};

	struct BrickLocation
	{
		int level = -1;
		int brick = -1;
		Int3 cell{}; // interior cell index in [0,brick_size)
		bool found() const { return brick >= 0; }
	};

	class AmrHierarchy
	{
	public:
		static AmrHierarchy uniform(const Aabb3d& requested_domain, double cell_size, int brick_size = 32, int ghost_cells = 1);
		static AmrHierarchy build_static(const Aabb3d& requested_domain, const TriMesh& wing, const TriangleBvh& bvh,
			const AmrConfig& config, bool anchor_y_min = false);

		const Aabb3d& domain() const { return domain_; }
		int brick_size() const { return brick_size_; }
		int ghost_cells() const { return ghost_cells_; }
		double finest_cell_size() const { return levels_.empty() ? 0.0 : levels_.back().h; }
		const std::vector<AmrLevel>& levels() const { return levels_; }
		std::vector<AmrLevel>& levels() { return levels_; }
		BrickLocation locate_finest(Vec3d point) const;
		int find_brick(int level, Int3 coord) const;
		bool is_two_to_one_balanced() const;
		std::size_t active_brick_count() const;
		std::size_t active_cell_count() const;

	private:
		static std::uint64_t coord_hash(Int3 c);
		void initialize_base(const Aabb3d& requested_domain, double h, bool anchor_y_min = false);
		void rebuild_level_tables_and_metadata(const TriangleBvh* bvh);
		void refine_brick(int level, int brick_id);
		void enforce_balance();
		bool point_in_domain(Vec3d p) const;

		Aabb3d domain_;
		int brick_size_ = 32;
		int ghost_cells_ = 1;
		int max_levels_ = 1;
		std::vector<AmrLevel> levels_;
	};
}
