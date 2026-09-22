// embedded_boundary.h — closed-solid topology on a uniform Cartesian level.
//
// Ordinary fluid cells/faces remain implicit. Only cells cut by the solid and faces touching them are
// represented explicitly. A regular fragment reference is encoded from its Cartesian cell; an
// irregular reference addresses the compact exterior-fluid fragment array.
#pragma once

#include "core/geometry/closed_solid.h"

#include <array>
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

	enum class EbCellState : std::uint8_t
	{
		regular = 0,
		split = 1,
		unresolved = 2,
		solid = 3 // inactive aerodynamic body: no fluid fragment or pressure DOF
	};

	struct EbCellTopology
	{
		EbCellState state = EbCellState::regular;
		int first_fragment = -1;
		std::uint16_t fragment_count = 0;
		int exact_locator_index = -1; // convex cut-cell fragment locator, host only
	};

	struct EbOrientedBoundaryTriangle
	{
		Vec3d a{}, b{}, c{}; // outward orientation of one closed fluid fragment
	};

	struct EbConvexHalfspace
	{
		Vec3d outward_normal{};
		double offset=0.0; // inside when dot(normal,point) <= offset
	};

	struct EbConvexRegion
	{
		std::vector<EbConvexHalfspace> halfspaces;
	};

	struct ExactEbFragmentBoundary
	{
		Vec3d interior_witness{};
		std::vector<EbOrientedBoundaryTriangle> triangles;
		// Fast polygonal cut cells are unions of convex BSP atoms. Their native
		// halfspaces give a more accurate point locator than rebuilding a winding
		// shell from independently tiled atom faces.
		std::vector<EbConvexRegion> convex_regions;
	};

	struct ExactEbCellLocator
	{
		std::vector<ExactEbFragmentBoundary> fragments; // local fragment order
		double ambiguity_tolerance = 0.0;
		double winding_tolerance = 1.0e-6;
	};

	struct FluidFragment
	{
		int parent_cell = -1;
		double volume = 0.0;
		Vec3d centroid{};
		FragmentRef merge_target = invalid_fragment; // self when unmerged; never crosses the solid wall
		int pressure_dof = -1;
		int connection_offset = 0, connection_count = 0;
		bool pressure_static = false; // isolated sealed pocket: retained state, excluded from projection
		// A topology-safe sub-threshold root that could not be agglomerated without
		// reversing another aperture. The production face-centred solver retains its
		// independent pressure DOF and advances only conservative aperture/dual fluxes;
		// no explicit update is divided by this small volume.
		bool face_state_retained = false;
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

	struct BoundaryFaceAperture
	{
		FragmentRef fragment = invalid_fragment;
		std::int8_t axis = 0, direction = 1; // outward axis normal
		double area = 0;
		Vec3d centroid{};
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
		// True when source_triangles is the complete local display-mesh input used while
		// diagnosing this cell, rather than a producer-certified causal subset.  UI
		// diagnostics must not paint this broad context as confirmed bad geometry.
		bool source_triangles_are_candidates = false;
	};

	struct EmbeddedBoundary
	{
		UniformEbGrid grid;
		std::vector<EbCellTopology> cells; // compact flags/offset per Cartesian cell
		// Bits 0/1/2 suppress the implicit +X/+Y/+Z regular face of this cell.
		// This handles the solid wall exactly coincident with a Cartesian face without inventing
		// zero-volume fragments. Any remaining open pieces use explicit apertures.
		std::vector<std::uint8_t> cut_face_mask;
		std::vector<int> irregular_cells;
		std::vector<FluidFragment> fragments;
		std::vector<FragmentConnection> connections; // only faces touching irregular topology
		std::vector<FaceAperture> apertures;
		std::vector<BoundaryFaceAperture> boundary_apertures;
		std::vector<SurfacePatch> patches;
		std::vector<ExactEbCellLocator> exact_locators;
		std::vector<UnresolvedEbCell> unresolved;
		double min_volume_fraction = 0.25;
		double min_aperture_area_fraction = 1e-4;
		// Apertures below the configured reporting scale remain physical flux paths.
		// These counters make their numerical cost visible without changing topology.
		std::size_t retained_subgrid_apertures = 0;
		double retained_subgrid_aperture_area = 0.0;
		// Every atom-pair interface is one polygon shared by both sides and remains
		// conservative as later planes partition it. These counters expose roundoff in
		// that common partition; it cannot redefine the solid or disable a CFD cell.
		std::size_t reconciled_facet_pairs = 0;
		double reconciled_facet_area_residual = 0.0;
		double maximum_reconciled_facet_area_residual = 0.0;
		// Grazing plane cuts whose discarded child is below the cell-volume
		// resolution are not committed. The parent remains a closed convex atom.
		std::size_t suppressed_subresolution_atom_splits = 0;
		double suppressed_subresolution_atom_volume = 0.0;
		double maximum_suppressed_subresolution_atom_volume = 0.0;

		FragmentRef fragment_containing_point(int cell, Vec3d point,
			double tolerance = 0.0) const;
		Vec3d fragment_centroid(FragmentRef ref) const;
		double fragment_volume(FragmentRef ref) const;
		bool ready_for_flow() const { return unresolved.empty(); }
	};

	struct EmbeddedBoundaryBuildOptions
	{
		double min_volume_fraction = 0.25;
		// Agglomeration is preferred, but a real exterior-fluid connector around the solid
		// can touch several independent pressure sides and have no legal merge target.
		// The production face-centred MAC path may retain such a raw root only when all
		// of its apertures remain positive, owned, and signed-bracketed. Collocated
		// volume-state transport must leave this disabled until it has redistribution.
		bool retain_signed_bracketed_small_roots_for_face_state = false;
		// Reporting scale for small positive-area apertures, as a fraction of h^2.
		// It is not a topology threshold: a physical opening remains connected even
		// when one or every subdivision piece lies below this value.  Coalescing is
		// legal only with an independent topology/provenance certificate and exact
		// preservation of area and first moment; the current implementation keeps
		// every piece.
		double min_aperture_area_fraction = 1e-4;
		// Required certificate that STEP import produced one valid closed solid.
		ClosedSolidGeometryPtr closed_solid;
	};

	// Fast rigid-body path. Exact triangle/box clipping supplies the wall surface;
	// locally grouped CAD-face planes form convex fluid fragments without CAD Booleans.
	EmbeddedBoundary build_closed_solid_embedded_boundary(const TriMesh& display_mesh,
		const TriangleBvh& broad_phase, const ClosedSolidGeometry& solid,
		const UniformEbGrid& grid, const EmbeddedBoundaryBuildOptions& options = {});
}
