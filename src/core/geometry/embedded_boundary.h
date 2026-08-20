// embedded_boundary.h — zero-thickness, two-sided fabric topology on a uniform Cartesian level.
//
// Ordinary cells/faces remain implicit. Only cells divided by fabric and faces touching them are
// represented explicitly. A regular fragment reference is encoded from its Cartesian cell; an
// irregular reference addresses the compact fragment array. There is deliberately no "solid side".
#pragma once

#include "core/geometry/exact_cell_decomposition.h"
#include "core/geometry/local_surface_arrangement.h"

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
		int arrangement_index = -1; // exact finite-triangle arrangement for complex cells
		int exact_locator_index = -1; // OCC-produced closed fragment boundaries, host only
		std::uint32_t source_face_id = ~std::uint32_t{0}; // common CAD face for an analytic local sheet
		// Range in EmbeddedBoundary::analytic_plane_support_triangles.  The analytic
		// plane fit deliberately remains smooth, while this exact finite-triangle
		// lineage lets a shared-face transaction reuse its once-derived physical trace.
		int plane_support_offset = -1;
		std::uint32_t plane_support_count = 0;
	};

	struct EbOrientedBoundaryTriangle
	{
		Vec3d a{}, b{}, c{}; // outward orientation of one closed fluid fragment
	};

	struct ExactEbFragmentBoundary
	{
		Vec3d interior_witness{};
		std::vector<EbOrientedBoundaryTriangle> triangles;
	};

	struct ExactEbCellLocator
	{
		std::vector<ExactEbFragmentBoundary> fragments; // local fragment order
		double ambiguity_tolerance = 0.0;
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
		int surface_side_offset = 0, surface_side_count = 0; // host preprocessing provenance
		// A topology-safe sub-threshold root that could not be agglomerated without
		// reversing another aperture. The production face-centred solver retains its
		// independent pressure DOF and advances only conservative aperture/dual fluxes;
		// no explicit update is divided by this small volume.
		bool face_state_retained = false;
	};

	struct FragmentSurfaceSide
	{
		std::uint32_t source_face_id = 0;
		std::uint8_t side_mask = 0; // bit 0: minus, bit 1: plus; both is valid around a sheet termination
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
		// Host-preprocessing ownership for a Cartesian-aligned patch staged before the
		// neighbouring cell topology is known. Exact shared-face assembly replaces only
		// records with this matching owner; production patches keep -1.
		int provisional_face_cell = -1;
		std::int8_t provisional_face_axis = -1;
	};

	struct UnresolvedEbCell
	{
		int parent_cell = -1;
		std::vector<std::uint32_t> source_triangles;
		std::string reason;
		// True when source_triangles is the complete local fabric input used while
		// diagnosing this cell, rather than a producer-certified causal subset.  UI
		// diagnostics must not paint this broad context as confirmed bad geometry.
		bool source_triangles_are_candidates = false;
	};

	// Exact duplicate CAD triangles on one Cartesian face represent one physical
	// zero-thickness barrier.  The canonical transaction retains one deterministic
	// load carrier and records every coincident source here so that the discarded
	// provenance is observable instead of silently double-counting wall area.
	struct CoincidentSurfaceProvenance
	{
		int parent_face_cell = -1;
		std::int8_t axis = -1;
		std::uint32_t canonical_source_triangle_id = 0;
		std::vector<std::uint32_t> coincident_source_triangle_ids;
		std::vector<std::uint32_t> coincident_source_face_ids;
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
		std::vector<FragmentSurfaceSide> fragment_surface_sides;
		std::vector<FragmentConnection> connections; // only faces touching irregular topology
		std::vector<FaceAperture> apertures;
		std::vector<SurfacePatch> patches;
		std::vector<FragmentRef> sampled_voxel_fragments;
		std::vector<LocalSurfaceArrangement> arrangements;
		std::vector<ExactEbCellLocator> exact_locators;
		std::vector<std::uint32_t> analytic_plane_support_triangles;
		std::vector<UnresolvedEbCell> unresolved;
		std::vector<CoincidentSurfaceProvenance> coincident_surface_provenance;
		double min_volume_fraction = 0.25;
		double min_aperture_area_fraction = 1e-4;
		// Apertures below the configured reporting scale remain physical flux paths.
		// These counters make their numerical cost visible without changing topology.
		std::size_t retained_subgrid_apertures = 0;
		double retained_subgrid_aperture_area = 0.0;
		// Reserved for a future topology-certified coalescing policy.  Measure alone
		// must never increment these counters or remove an aperture.
		std::size_t discarded_subgrid_apertures = 0;
		double discarded_subgrid_aperture_area = 0.0;
		// Positive common-refinement carriers removed only by the auditable triple
		// certificate: conflicting local sides on a closed fan-2 CAD face plus a
		// fabric-crossing finite-volume leg. Their area remains in the explicit
		// accepted+rejected face-partition conservation invariant.
		std::size_t rejected_cross_fabric_apertures = 0;
		double rejected_cross_fabric_aperture_area = 0.0;
		std::size_t recovered_subgrid_surface_sides = 0;
		double maximum_surface_side_recovery_distance = 0.0;
		// Qualitative preview only: a rasterized sliver can have no represented
		// pressure sample on one side. Its impermeable topology is retained, but the
		// unmappable patch is omitted from pressure-load/wall integration and counted.
		std::size_t diagnostic_unmapped_surface_patches = 0;
		double diagnostic_unmapped_surface_area = 0.0;
		// A sampled cell may use one pressure fragment on both sides only when its
		// fluid component reaches a geometrically confirmed free fabric edge inside
		// that cell.  Rejected counters make seam/junction raster leaks auditable.
		std::size_t accepted_free_edge_same_fragment_patches = 0;
		double accepted_free_edge_same_fragment_area = 0.0;
		std::size_t rejected_same_fragment_patches = 0;
		double rejected_same_fragment_area = 0.0;
		std::size_t ambiguous_edge_collapse_cells = 0;
		std::size_t diagnostic_unverified_same_fragment_patches = 0;
		double diagnostic_unverified_same_fragment_area = 0.0;
		std::size_t exact_decomposition_cells = 0;
		double exact_decomposition_milliseconds = 0.0;

		FragmentRef fragment_for_side(int cell, int side) const;
		FragmentRef fragment_containing_point(int cell, Vec3d point,
			double tolerance = 0.0) const;
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
		double min_volume_fraction = 0.25;
		// Agglomeration is preferred, but a real fluid connector around a fabric edge
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
		// Pressure-solver diagnostics sometimes need the historical, leaky topology
		// to reproduce a numerical failure independently of preprocessing. This is
		// intentionally absent from product configuration and defaults fail-closed.
		bool allow_unverified_same_fragment_patches_for_diagnostics = false;
		// Optional OCC-free callback seam for the OCCT General Fuse development
		// oracle.  Production preprocessing uses the deterministic finite-triangle
		// LocalSurfaceArrangement for complex cells; merely supplying this callback
		// must not divert those cells through General Fuse.  Tests and offline
		// diagnostics that intentionally exercise the oracle must also enable the
		// explicit flag below.  The function implementation lives in
		// paracfd_geometry; libparacfd stores and consumes only value-type results,
		// preserving the one-way CAD -> static-topology dependency.
		ExactCellDecomposer exact_cell_decomposer = nullptr;
		bool use_exact_cell_decomposer_as_development_oracle = false;
	};

	EmbeddedBoundary build_embedded_boundary(const TriMesh& mesh, const TriangleBvh& bvh,
		const UniformEbGrid& grid, const EmbeddedBoundaryBuildOptions& options = {});
}
