// local_surface_arrangement.h -- exact local topology of finite zero-thickness triangles.
//
// This is CPU-only preprocessing for the small minority of Cartesian cells containing
// seams, junctions, terminating sheets, or several fabric panels.  It constructs the
// arrangement of the triangles' physical support planes inside one cell, then takes a
// finite-polygon common refinement on every support facet.  It joins convex atoms wherever
// that facet has positive-area open fluid.  It therefore
// represents the complement of the finite triangle union; it never performs solid fill
// and it has no two-fragment assumption.
#pragma once

#include "core/geometry/triangle_bvh.h"

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace paracfd::core
{
	struct LocalSurfaceTriangle
	{
		Vec3d a{}, b{}, c{};
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0; // provenance only; never used for connectivity
		std::array<FabricEdgeKind, 3> edge_kind{
			FabricEdgeKind::unknown, FabricEdgeKind::unknown, FabricEdgeKind::unknown};
		// Optional caller-provided seam group.  Equal non-sentinel values identify
		// half-edges belonging to the same certified attachment.  Geometry is still
		// checked; this provenance never joins separated fluid on its own.
		std::array<std::uint64_t, 3> edge_attachment_id{
			std::numeric_limits<std::uint64_t>::max(),
			std::numeric_limits<std::uint64_t>::max(),
			std::numeric_limits<std::uint64_t>::max()};
	};

	struct LocalSurfaceArrangementOptions
	{
		// Coincident planes closer than contact_tolerance are one geometric plane.  Two
		// nearly coincident planes inside ambiguity_tolerance are rejected rather than
		// silently sewn together.  Zero selects a scale-aware FP64 default; callers that
		// received FP32 CAD coordinates should pass the BVH contact/clearance tolerances.
		double contact_tolerance = 0.0;
		double ambiguity_tolerance = 0.0;
		double angular_contact_tolerance = 0.0;
		double angular_ambiguity_tolerance = 0.0;
		// Cartesian-face alignment is a numerical clipping property, not a CAD
		// attachment/healing tolerance.  Keeping it separate prevents a near-face
		// triangle accepted by a generous CAD contact tolerance from being emitted as
		// face-coplanar fabric without a matching physical barrier.  Zero selects a
		// scale-aware FP64 default.
		double cartesian_face_tolerance = 0.0;
		double angular_cartesian_face_tolerance = 0.0;
		double minimum_area_fraction = 1e-13;
		double minimum_volume_fraction = 1e-15;
		int maximum_planes = 64;
		int maximum_atoms = 4096;
	};

	struct LocalArrangementPlane
	{
		Vec3d normal{}; // deterministic canonical orientation
		double offset = 0.0; // dot(normal,x)-offset
		// Original indices in the build input, including any earlier entries skipped
		// as degenerate/outside.  These are diagnostics/provenance, not compressed
		// PreparedTriangle indices.
		std::vector<int> support_triangles;
		// Stable mesh triangle provenance corresponding to `support_triangles`.
		// Unlike the build-input indices above, these IDs can be joined to the
		// once-derived physical trace on a neighbouring Cartesian cell face.
		std::vector<std::uint32_t> support_source_triangles;
		// Retained for diagnostics/file compatibility.  Finite triangle edges are handled
		// by the 2-D common refinement on the support facet; extending every edge to an
		// infinite 3-D partition plane causes an unrelated O(N^3) atom arrangement.
		bool has_edge_limiter_role = false;
	};

	struct LocalArrangementFace
	{
		std::vector<Vec3d> polygon;
		// >=0: local arrangement plane.  -1..-6: cell -X,+X,-Y,+Y,-Z,+Z.
		int plane_id = -1;
		// Combinatorial identity of the original split cap.  Both sides receive the
		// same token and retain it through later clipping, so adjacency does not have
		// to be inferred from independently rounded sign vectors.
		std::uint64_t lineage = 0;
	};

	struct LocalArrangementAtom
	{
		double volume = 0.0;
		Vec3d centroid{};
		// -1/+1 for a represented split; 0 when a grazing sub-volume-tolerance
		// split was deliberately retained as one atom and therefore straddles the plane.
		std::vector<std::int8_t> plane_side;
		std::vector<LocalArrangementFace> faces;
		int fragment = -1;
	};

	struct LocalArrangementFragment
	{
		double volume = 0.0;
		Vec3d centroid{};
		std::vector<int> atoms;
	};

	struct LocalArrangementSurfacePatch
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		double area = 0.0;
		Vec3d centroid{};
		Vec3d normal{}; // input winding, minus -> plus
		std::vector<Vec3d> polygon;
		int plus_fragment = -1;
		int minus_fragment = -1;
	};

	// An exact, connected piece of open fluid on one Cartesian cell face.  Neighbouring
	// cells form FaceApertures by intersecting pieces on their common face.  Separate
	// polygons are deliberately retained even when they have the same fragment pair.
	struct LocalArrangementBoundaryAperture
	{
		std::int8_t axis = 0;
		bool upper = false;
		double area = 0.0;
		Vec3d centroid{};
		std::vector<Vec3d> polygon;
		int fragment = -1;
		double contact_tolerance = 0.0;
	};

	// Fabric exactly coincident with a Cartesian face has only one local fluid side.  The
	// neighbouring cell supplies the other side when the per-cell arrangements are joined.
	struct LocalArrangementBoundarySurfacePatch
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		std::int8_t axis = 0;
		bool upper = false;
		double area = 0.0;
		Vec3d centroid{};
		Vec3d normal{};
		std::vector<Vec3d> polygon;
		int interior_fragment = -1;
		bool interior_is_plus = false;
		double contact_tolerance = 0.0;
	};

	// An edge whose attachment/free status is geometrically ambiguous but which lies
	// wholly on a Cartesian cell face.  It cannot alter connectivity inside this cell,
	// so resolution is deferred to the two-cell shared-face transaction.  Interior
	// unknown edges are rejected during local arrangement construction.
	struct LocalArrangementBoundaryEdgeHazard
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		std::int8_t axis = 0;
		bool upper = false;
		Vec3d a{}, b{};
		double contact_tolerance = 0.0;
		double ambiguity_tolerance = 0.0;
	};

	struct LocalArrangementBoundaryOverlap
	{
		std::int8_t axis = 0;
		double area = 0.0;
		Vec3d centroid{};
		std::vector<Vec3d> polygon;
		int fragment_a = -1;
		int fragment_b = -1;
	};

	// One common subdivision of fabric lying exactly on a face shared by two
	// Cartesian cells.  It is emitted once, even when the two cells partition that
	// face differently.  `owner_is_a` selects the lower-coordinate cell and is
	// therefore invariant to the order in which the two inputs are passed.
	struct LocalArrangementBoundarySurfaceOverlap
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		std::int8_t axis = 0;
		double area = 0.0;
		Vec3d centroid{};
		Vec3d normal{}; // follows side A's winding; minus -> plus
		std::vector<Vec3d> polygon;
		int plus_fragment = -1;
		int minus_fragment = -1;
		bool owner_is_a = false;
	};

	struct LocalArrangementBoundarySurfaceAssembly
	{
		bool valid = false;
		std::string error;
		double area = 0.0;
		std::vector<LocalArrangementBoundarySurfaceOverlap> patches;
	};

	// A physical finite-triangle trace on a Cartesian face.  A transverse sheet
	// contributes a segment; it subdivides adjacent-fluid ownership but has zero
	// blocked area on the face.  Provenance is retained for diagnostics only.
	struct LocalArrangementSharedFaceTrace
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		Vec3d a{}, b{};
		Vec3d normal{}; // source triangle winding, retained for diagnostics
		// Build-local reciprocal half-edge identity when this trace is exactly an
		// attached mesh edge.  Equal non-sentinel values certify one physical line;
		// geometric proximity alone never grants that identity.
		std::uint64_t attachment_id = std::numeric_limits<std::uint64_t>::max();
	};

	// Positive-area fabric coincident with the Cartesian face.  The polygon is a
	// physical no-flux barrier.  Unlike a transverse trace, it removes its covered
	// area from the open Cartesian-face aperture set.
	struct LocalArrangementSharedFaceBarrier
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		Vec3d normal{}; // input winding, minus -> plus
		std::vector<Vec3d> polygon;
	};

	// One polygon in a complete adjacent-fluid ownership cover of the shared face.
	// Callers form this cover from boundary apertures (fragment) and boundary
	// surface patches (interior_fragment).  The negative-axis cell is side A and
	// the positive-axis cell is side B.  Fragment labels are opaque: connectivity
	// is decided by the canonical geometry, never by CAD face/source identity.
	constexpr int no_shared_face_fragment = std::numeric_limits<int>::min();
	struct LocalArrangementSharedFaceOwner
	{
		std::vector<Vec3d> polygon;
		int fragment = no_shared_face_fragment;
	};

	// Structural description of the fluid regions adjacent to one Cartesian face.
	// Unlike LocalArrangementSharedFaceOwner this is not a polygon reconstructed from
	// an already clipped cell boundary.  It is the same implicit partition used to
	// construct the cell's pressure fragments, so independently clipped polygons never
	// become the authority for cross-cell connectivity.
	struct LocalArrangementSharedFacePartitionPlane
	{
		Vec3d normal{};
		double offset = 0.0; // dot(normal, world_point) - offset
		// Original source triangles supporting this local plane.  The shared-face
		// transaction uses this lineage to reuse an exact physical face trace when
		// independently normalized plane equations differ by roundoff.
		std::vector<std::uint32_t> support_triangles;
	};

	struct LocalArrangementSharedFaceRegion
	{
		// One entry per plane, in the same order as `planes`.  -1/+1 selects a
		// half-space.  Zero is the explicit wildcard used by a local arrangement when
		// an arithmetic-volume grazing split was deliberately retained unsplit.
		std::vector<std::int8_t> plane_side;
		int fragment = no_shared_face_fragment;
	};

	struct LocalArrangementSharedFaceSide
	{
		std::vector<LocalArrangementSharedFacePartitionPlane> planes;
		std::vector<LocalArrangementSharedFaceRegion> regions;
		// Direction from the common face into this cell along `axis`.  For the
		// lower-coordinate (A) cell this is -1; for the upper-coordinate (B) cell
		// this is +1.  It selects the fluid half-space of a support plane exactly
		// coplanar with the Cartesian face without an epsilon-offset point sample.
		std::int8_t inward_axis_sign = 0;
	};

	struct LocalArrangementSharedFaceAperture
	{
		double area = 0.0;
		Vec3d centroid{};
		std::vector<Vec3d> polygon;
		int fragment_a = no_shared_face_fragment;
		int fragment_b = no_shared_face_fragment;
	};

	struct LocalArrangementSharedFaceBlockedSurface
	{
		std::uint32_t source_triangle_id = 0;
		std::uint32_t source_face_id = 0;
		double area = 0.0;
		Vec3d centroid{};
		Vec3d normal{}; // selected barrier winding, minus -> plus
		std::vector<Vec3d> polygon;
		int plus_fragment = no_shared_face_fragment;
		int minus_fragment = no_shared_face_fragment;
	};

	struct LocalArrangementSharedFaceOptions
	{
		// A finite supporting line may over-subdivide the face beyond the source
		// segment.  This changes neither topology nor measure, but the cap keeps a
		// malformed or pathologically dense input transactional and bounded.
		int maximum_tiles = 16384;
	};

	struct LocalArrangementSharedFaceAssembly
	{
		bool valid = false;
		std::string error;
		std::vector<LocalArrangementSharedFaceAperture> apertures;
		std::vector<LocalArrangementSharedFaceBlockedSurface> blocked_surfaces;
		double face_area = 0.0;
		double open_area = 0.0;
		double blocked_area = 0.0;
		double area_conservation_error = 0.0;
		double first_moment_conservation_error = 0.0;
		double area_conservation_tolerance = 0.0;
		double first_moment_conservation_tolerance = 0.0;
	};

	struct LocalSurfaceArrangement
	{
		Aabb3d cell;
		bool valid = false;
		std::string error;
		std::vector<LocalArrangementPlane> planes;
		std::vector<LocalArrangementAtom> atoms;
		std::vector<LocalArrangementFragment> fragments;
		std::vector<LocalArrangementSurfacePatch> surface_patches;
		std::vector<LocalArrangementBoundaryAperture> boundary_apertures;
		std::vector<LocalArrangementBoundarySurfacePatch> boundary_surface_patches;
		std::vector<LocalArrangementBoundaryEdgeHazard> boundary_edge_hazards;
		double contact_tolerance = 0.0;
		double ambiguity_tolerance = 0.0;
		double cartesian_face_tolerance = 0.0;
		double angular_cartesian_face_tolerance = 0.0;
		double represented_surface_area = 0.0;
		double clipped_input_surface_area = 0.0;
		double volume_conservation_error = 0.0;
		double first_moment_conservation_error = 0.0;
		double volume_conservation_tolerance = 0.0;
		double first_moment_conservation_tolerance = 0.0;
	};

	LocalSurfaceArrangement build_local_surface_arrangement(
		const Aabb3d& cell, const std::vector<LocalSurfaceTriangle>& triangles,
		const LocalSurfaceArrangementOptions& options = {});

	// Host-side classifier for probes/visualisation and Cartesian-face staging.  A point
	// on an actual finite fabric patch always returns -1.  A point only on an artificial
	// partition plane is accepted when every incident atom belongs to the same fragment;
	// callers selecting a fabric side should offset along its normal.
	int locate_local_arrangement_fragment(const LocalSurfaceArrangement& arrangement,
		Vec3d point, double tolerance = 0.0);

	// Intersect two convex coplanar Cartesian-face pieces.  This is the staging primitive
	// used to create exact cross-cell apertures without an N x N face raster.
	bool intersect_local_boundary_apertures(const LocalArrangementBoundaryAperture& a,
		const LocalArrangementBoundaryAperture& b, LocalArrangementBoundaryOverlap& overlap,
		double tolerance = 0.0);

	// Assemble two lists belonging to opposite sides of one common Cartesian face.
	// Every positive-area piece must have a provenance-matched counterpart and the
	// total area must agree on both sides; otherwise the result is invalid.
	LocalArrangementBoundarySurfaceAssembly assemble_local_boundary_surface_patches(
		const std::vector<LocalArrangementBoundarySurfacePatch>& a,
		const std::vector<LocalArrangementBoundarySurfacePatch>& b,
		double tolerance = 0.0);

	// Construct one canonical subdivision for a face shared by two Cartesian cells.
	// `face_square` lies on `axis`; side A is the negative-axis cell and side B the
	// positive-axis cell.  Physical traces/barriers are supplied once, while each
	// side supplies a complete polygonal adjacent-fluid ownership cover.  Every
	// positive canonical tile is retained.  The routine fails transactionally when
	// a tile lacks exactly one owner on either side, has overlapping barriers, or
	// violates area/first-moment conservation.
	LocalArrangementSharedFaceAssembly assemble_canonical_local_shared_face(
		const std::array<Vec3d, 4>& face_square, std::int8_t axis,
		const std::vector<LocalArrangementSharedFaceTrace>& traces,
		const std::vector<LocalArrangementSharedFaceBarrier>& barriers,
		const std::vector<LocalArrangementSharedFaceOwner>& owners_a,
		const std::vector<LocalArrangementSharedFaceOwner>& owners_b,
		const LocalArrangementSharedFaceOptions& options = {});

	// Production shared-face transaction.  All non-coplanar support planes from both
	// cells are intersected with the common face to make one canonical 2-D tiling.
	// Every positive-area tile is labelled by its complete half-space sign vector and
	// must match exactly one region on each side.  Fabric barriers alone remove open
	// area; a transverse support plane partitions ownership but never blocks flux.
	//
	// The polygon-cover overload above remains useful for diagnostics and malformed
	// cover tests.  Production topology should use this structural overload because it
	// does not depend on lossy, independently clipped boundary-owner polygons.
	LocalArrangementSharedFaceAssembly assemble_canonical_local_shared_face(
		const std::array<Vec3d, 4>& face_square, std::int8_t axis,
		const std::vector<LocalArrangementSharedFaceBarrier>& barriers,
		const LocalArrangementSharedFaceSide& side_a,
		const LocalArrangementSharedFaceSide& side_b,
		const LocalArrangementSharedFaceOptions& options = {});

	// Production overload with the once-derived physical transverse traces.  A
	// certified attached trace supplies the canonical representative line for every
	// incident support plane, preventing arithmetic slivers at seams and junctions.
	LocalArrangementSharedFaceAssembly assemble_canonical_local_shared_face(
		const std::array<Vec3d, 4>& face_square, std::int8_t axis,
		const std::vector<LocalArrangementSharedFaceTrace>& traces,
		const std::vector<LocalArrangementSharedFaceBarrier>& barriers,
		const LocalArrangementSharedFaceSide& side_a,
		const LocalArrangementSharedFaceSide& side_b,
		const LocalArrangementSharedFaceOptions& options = {});
}
