// triangle_bvh.h — static, dependency-free acceleration over a placed triangle mesh.
//
// Construction and queries use double precision because this is preprocessing and collision
// geometry, not a CFD timestep field. The BVH owns a compact copy of its triangles so callers may
// replace/free the source TriMesh after construction. Runtime CFD uploads only derived float EB data.
#pragma once

#include "core/geometry/tri_mesh.h"

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace paracfd::core
{
	struct Vec3d
	{
		double x = 0.0, y = 0.0, z = 0.0;
		double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
		double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
	};

	inline Vec3d operator+(Vec3d a, Vec3d b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
	inline Vec3d operator-(Vec3d a, Vec3d b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
	inline Vec3d operator*(Vec3d a, double s) { return {a.x * s, a.y * s, a.z * s}; }
	inline Vec3d operator/(Vec3d a, double s) { return {a.x / s, a.y / s, a.z / s}; }
	inline double dot(Vec3d a, Vec3d b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
	inline Vec3d cross(Vec3d a, Vec3d b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
	inline double length2(Vec3d a) { return dot(a, a); }
	Vec3d normalized(Vec3d a);

	struct Aabb3d
	{
		Vec3d lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
		Vec3d hi{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};

		void expand(Vec3d p);
		void expand(const Aabb3d& b);
		bool valid() const;
		bool overlaps(const Aabb3d& b) const;
		double distance2(Vec3d p) const;
	};

	struct BvhTriangle
	{
		Vec3d a, b, c;
		std::uint32_t triangle_id = 0;
		std::uint32_t source_face_id = 0;
	};

	struct SegmentHit
	{
		bool hit = false;
		double t = 0.0; // p(t)=a+t(b-a), [0,1]
		double u = 0.0, v = 0.0; // barycentric weights for b,c; weight(a)=1-u-v
		Vec3d position{};
		Vec3d geometric_normal{}; // follows triangle winding
		std::uint32_t triangle_id = 0;
		std::uint32_t source_face_id = 0;
	};

	struct NearestSurfacePoint
	{
		bool found = false;
		double distance = std::numeric_limits<double>::infinity();
		Vec3d point{};
		Vec3d geometric_normal{};
		std::uint32_t triangle_id = 0;
		std::uint32_t source_face_id = 0;
	};

	// Static topology of one tessellation edge.  `attached` includes ordinary
	// within-face tessellation edges as well as CAD seams and non-manifold fabric
	// junctions.  Only `confirmed_free` may be used to justify fluid connectivity
	// around a zero-thickness sheet termination.  Geometry that falls inside the
	// source-coordinate uncertainty band is deliberately `unknown`, never silently
	// treated as an opening or a closed seam.
	enum class FabricEdgeKind : std::uint8_t
	{
		attached = 0,
		confirmed_free = 1,
		unknown = 2
	};

	// Geometric role of one triangle half-edge.  Openness and role are deliberately
	// separate: both an ordinary tessellation diagonal and a sewn CAD seam are
	// `attached`, but only a junction requires general local topology.  `none` is the
	// role of a geometrically confirmed free edge.  `unknown` is fail-closed evidence,
	// never an opening.
	enum class FabricEdgeRole : std::uint8_t
	{
		none = 0,
		tessellation_interior = 1,
		manifold_seam = 2,
		junction = 3,
		unknown = 4
	};

	// Diagnostic reason retained for fail-closed edge certificates. This is not used
	// to relax topology; it makes imported-CAD failures reproducible without rerunning
	// the classifier under a debugger.
	enum class FabricEdgeUnknownReason : std::uint8_t
	{
		none = 0,
		invalid_triangle = 1,
		source_provenance = 2,
		indexed_fan_mismatch = 3,
		proximity_ambiguity = 4,
		mixed_open_attachment = 5,
		cad_fan_mismatch = 6,
		no_geometric_evidence = 7,
		nonreciprocal_attachment = 8
	};

	struct FabricEdgeCertificate
	{
		static constexpr std::uint64_t no_attachment =
			std::numeric_limits<std::uint64_t>::max();

		FabricEdgeKind kind = FabricEdgeKind::unknown;
		FabricEdgeRole role = FabricEdgeRole::unknown;
		FabricEdgeUnknownReason unknown_reason = FabricEdgeUnknownReason::invalid_triangle;
		// Build-local equality token.  Every reciprocal half-edge in one certified
		// attachment has the same value; it is not a persistent CAD identifier.
		std::uint64_t attachment_id = no_attachment;
		// Stable source-CAD contact curve, when OpenCascade certified that a nominally
		// free BRep edge lies on another trimmed face. Unlike attachment_id this does
		// not assert polygon-half-edge equality and is safe for curved/nonconforming
		// independent tessellations.
		std::uint64_t cad_contact_id = TriMesh::kNoCadContactId;
		// Number of incident half-sheet sectors certified along the edge.  A free
		// boundary has one, a manifold continuation/seam has two, and a junction has
		// three or more.  Zero means that incidence itself is unresolved.
		std::uint16_t incident_fan_degree = 0;
		// Tolerance actually needed to certify this relation.  The no-argument BVH
		// tolerance accessors remain the coordinate-precision baseline.
		double contact_tolerance = 0.0;
		double clearance_tolerance = 0.0;
	};

	class TriangleBvh
	{
	public:
		TriangleBvh() = default;
		explicit TriangleBvh(const TriMesh& mesh) { build(mesh); }

		void build(const TriMesh& mesh, std::uint32_t leaf_size = 8);
		bool empty() const { return triangles_.empty(); }
		std::size_t triangle_count() const { return triangles_.size(); }
		const BvhTriangle& triangle(std::uint32_t original_id) const { return triangles_by_id_[original_id]; }
		FabricEdgeCertificate edge_certificate(std::uint32_t original_triangle_id,
			int local_edge) const;
		FabricEdgeKind edge_kind(std::uint32_t original_triangle_id, int local_edge) const;
		FabricEdgeRole edge_role(std::uint32_t original_triangle_id, int local_edge) const
		{
			return edge_certificate(original_triangle_id, local_edge).role;
		}
		std::uint64_t edge_attachment_id(std::uint32_t original_triangle_id, int local_edge) const
		{
			return edge_certificate(original_triangle_id, local_edge).attachment_id;
		}
		std::uint16_t edge_incident_fan_degree(std::uint32_t original_triangle_id,
			int local_edge) const
		{
			return edge_certificate(original_triangle_id, local_edge).incident_fan_degree;
		}
		double edge_contact_tolerance(std::uint32_t original_triangle_id, int local_edge) const
		{
			return edge_certificate(original_triangle_id, local_edge).contact_tolerance;
		}
		double edge_clearance_tolerance(std::uint32_t original_triangle_id, int local_edge) const
		{
			return edge_certificate(original_triangle_id, local_edge).clearance_tolerance;
		}
		// Baseline thresholds for geometry without per-edge CAD provenance. Contact is
		// double-precision arithmetic coincidence; float-only meshes retain a separate,
		// wider FP32 ambiguity/clearance band. Certified CAD tolerance is per edge above.
		double edge_contact_tolerance() const { return edge_contact_tolerance_; }
		double edge_clearance_tolerance() const { return edge_clearance_tolerance_; }

		// Exact triangle/AABB overlap, not merely BVH-node candidates. Results are original mesh
		// triangle IDs, sorted for deterministic preprocessing/tests.
		void query_aabb(const Aabb3d& box, std::vector<std::uint32_t>& out_triangle_ids) const;
		std::vector<std::uint32_t> query_aabb(const Aabb3d& box) const;

		// Nearest two-sided intersection on the closed segment. `t_min` is useful for ignoring a
		// launch point already on fabric; it is expressed as segment fraction, not metres.
		SegmentHit intersect_segment(Vec3d a, Vec3d b, double t_min = 1e-10, double t_max = 1.0) const;
		// Closed-segment proximity query used only while constructing static topology.  Unlike
		// Moller-Trumbore crossing, this also catches a segment that is coplanar with, ends on,
		// or follows a triangle edge.  The caller supplies a geometric roundoff tolerance; it
		// is not a finite-thickness surface model.
		bool segment_touches_surface(Vec3d a, Vec3d b, double tolerance) const;
		NearestSurfacePoint nearest(Vec3d p, double max_distance = std::numeric_limits<double>::infinity()) const;
		NearestSurfacePoint nearest_on_face(Vec3d p, std::uint32_t source_face_id,
			double max_distance = std::numeric_limits<double>::infinity()) const;
		double distance(Vec3d p, double max_distance = std::numeric_limits<double>::infinity()) const;

		static bool triangle_intersects_aabb(const BvhTriangle& tri, const Aabb3d& box);

	private:
		struct Primitive
		{
			BvhTriangle tri;
			Aabb3d bounds;
			Vec3d centroid;
		};
		struct Node
		{
			Aabb3d bounds;
			std::uint32_t first = 0, count = 0;
			std::uint32_t left = 0, right = 0;
			bool leaf() const { return count != 0; }
		};

		std::uint32_t build_node(std::uint32_t first, std::uint32_t count, std::uint32_t leaf_size);
		std::vector<Primitive> primitives_; // reordered into leaf ranges
		std::vector<Node> nodes_;
		std::vector<BvhTriangle> triangles_; // reordered copy (statistics/empty)
		std::vector<BvhTriangle> triangles_by_id_; // direct original-id lookup
		std::vector<std::array<FabricEdgeCertificate, 3>> edge_certificate_by_triangle_;
		double edge_contact_tolerance_ = 0.0;
		double edge_clearance_tolerance_ = 0.0;
	};
}
