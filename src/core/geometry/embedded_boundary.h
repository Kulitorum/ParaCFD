// embedded_boundary.h — zero-thickness, two-sided fabric topology on a uniform Cartesian level.
//
// Ordinary cells/faces remain implicit. Only cells divided by fabric and faces touching them are
// represented explicitly. A regular fragment reference is encoded from its Cartesian cell; an
// irregular reference addresses the compact fragment array. There is deliberately no "solid side".
#pragma once

#include "core/geometry/triangle_bvh.h"

#include <cstdint>
#include <string>
#include <vector>

#ifdef __CUDACC__
#define PARACFD_EB_HD __host__ __device__
#else
#define PARACFD_EB_HD
#endif

namespace paracfd::core
{
	using FragmentRef = std::int32_t;
	constexpr FragmentRef invalid_fragment = 0;
	inline FragmentRef regular_fragment(int cell) { return -static_cast<FragmentRef>(cell + 1); }
	inline FragmentRef irregular_fragment(int fragment) { return static_cast<FragmentRef>(fragment + 1); }
	inline bool fragment_is_regular(FragmentRef r) { return r < 0; }
	inline int regular_fragment_cell(FragmentRef r) { return -static_cast<int>(r) - 1; }
	inline int irregular_fragment_index(FragmentRef r) { return static_cast<int>(r) - 1; }

	struct UniformEbGrid
	{
		Vec3d origin{};
		int nx = 0, ny = 0, nz = 0;
		double h = 0.0;
		PARACFD_EB_HD int cell_count() const { return nx * ny * nz; }
		PARACFD_EB_HD int cell_index(int i, int j, int k) const { return (k * ny + j) * nx + i; }
		std::array<int, 3> cell_coord(int cell) const { return {cell % nx, (cell / nx) % ny, cell / (nx * ny)}; }
		Aabb3d cell_box(int i, int j, int k) const;
		Vec3d cell_centroid(int cell) const;
	};

	enum class EbCellState : std::uint8_t { regular = 0, split = 1, unresolved = 2 };

	struct EbCellTopology
	{
		EbCellState state = EbCellState::regular;
		int first_fragment = -1;
		std::uint16_t fragment_count = 0;
		Vec3d plane_normal{}; // topology plane, minus -> plus
		double plane_offset = 0.0; // dot(n,x)-offset
		int sampled_voxel_offset = -1; // fallback topology, resolution^3 FragmentRefs
		std::uint8_t sampled_resolution = 0;
	};

	struct FluidFragment
	{
		int parent_cell = -1;
		double volume = 0.0;
		Vec3d centroid{};
		std::int8_t side = 0; // -1 / +1 relative to cell topology plane
		FragmentRef merge_target = invalid_fragment; // self when unmerged; never crosses fabric
		int pressure_dof = -1;
		int connection_offset = 0, connection_count = 0;
		bool pressure_static = false; // isolated sealed pocket: retained state, excluded from projection
	};

	struct FragmentConnection
	{
		FragmentRef fragment_a = invalid_fragment;
		FragmentRef fragment_b = invalid_fragment;
		double open_area = 0.0;
		Vec3d face_centroid{};
		double centre_distance = 0.0;
		std::int8_t axis = 0; // 0=x,1=y,2=z
	};

	struct FaceAperture
	{
		int parent_face_cell = -1; // lower-index cell, or boundary owner
		std::int8_t axis = 0;
		double area = 0.0;
		Vec3d centroid{};
		FragmentRef fragment_a = invalid_fragment;
		FragmentRef fragment_b = invalid_fragment; // invalid at a physical boundary
	};

	struct SurfacePatch
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		double area = 0.0;
		Vec3d centroid{};
		Vec3d normal{}; // triangle winding: minus -> plus
		FragmentRef plus_fragment = invalid_fragment;
		FragmentRef minus_fragment = invalid_fragment;
	};

	struct UnresolvedEbCell
	{
		int parent_cell = -1;
		std::vector<std::uint32_t> source_triangles;
		std::string reason;
	};

	struct EmbeddedBoundary
	{
		UniformEbGrid grid;
		std::vector<EbCellTopology> cells; // compact flags/offset per Cartesian cell
		// Bits 0/1/2 suppress the implicit +X/+Y/+Z regular face of this cell.
		// This handles fabric exactly coincident with a Cartesian face without inventing
		// zero-volume fragments. Any remaining open pieces use explicit apertures.
		std::vector<std::uint8_t> cut_face_mask;
		std::vector<int> irregular_cells;
		std::vector<FluidFragment> fragments;
		std::vector<FragmentConnection> connections; // only faces touching irregular topology
		std::vector<FaceAperture> apertures;
		std::vector<SurfacePatch> patches;
		std::vector<FragmentRef> sampled_voxel_fragments;
		std::vector<UnresolvedEbCell> unresolved;
		double min_volume_fraction = 0.05;

		FragmentRef fragment_for_side(int cell, int side) const;
		Vec3d fragment_centroid(FragmentRef ref) const;
		double fragment_volume(FragmentRef ref) const;
		bool ready_for_flow() const { return unresolved.empty(); }
	};

	struct EmbeddedBoundaryBuildOptions
	{
		// A CAD canopy is curved, so a CFD cell commonly contains several non-coplanar
		// tessellation triangles belonging to one smooth sheet. Fit those triangles to
		// one local topology plane only while both bounds below hold. Ribs, converging
		// skins, and other distinct sheets exceed these bounds and remain unresolved.
		double smooth_sheet_angle_tolerance = 0.26; // radians (about 15 degrees)
		double smooth_sheet_distance_fraction = 0.15; // max vertex distance from fitted plane / h
		double coplanar_angle_tolerance = 1e-6; // exact Cartesian-face classification only
		double coplanar_distance_tolerance = 1e-8; // metres, exact-face classification only
		double surface_coverage_tolerance = 0.05; // projected coverage of fitted plane/cell section
		// Optional finest-level fallback for cells with junctions, terminating sheets,
		// or multiple surfaces. It reconstructs fluid connectivity (never solid fill)
		// on an N^3 local graph whose edges are blocked by exact BVH intersections.
		int complex_subdivisions = 0;
		double min_volume_fraction = 0.05;
	};

	EmbeddedBoundary build_embedded_boundary(const TriMesh& mesh, const TriangleBvh& bvh,
		const UniformEbGrid& grid, const EmbeddedBoundaryBuildOptions& options = {});
}
