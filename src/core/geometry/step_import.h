// step_import.h — OpenCascade STEP-file import to a Qt-free triangle mesh in SI METRES.
//
// This header is deliberately free of any OpenCascade (and Qt) include so the GUI and the
// CFD preprocessing and the viewer can consume the mesh without inheriting OCC's include path or host
// flags. All OCC usage is isolated in step_import.cpp, compiled into the `paracfd_geometry`
// static lib (see CMakeLists.txt); libparacfd and the physics gates stay OCC-free.
//
// ⚠ Units: OCC emits geometry in MILLIMETRES; load_step_mesh() scales every coordinate by
// 0.001 so the returned TriMesh is in METRES, ready for the SI simulation domain and
// zero-thickness embedded-boundary preprocessing. The float positions remain the
// display/GPU representation; the importer also retains the transformed OCCT coordinates
// in TriMesh::positions_fp64 for robust static geometry preprocessing.
//
// Only triangulated TopoDS_Face entities contribute geometry. Standalone STEP edges/wires
// (for example suspension lines) are intentionally ignored, so they cannot enlarge the
// aerodynamic bbox or become impermeable CFD surfaces.
#pragma once

#include "core/geometry/geometry_quality.h"
#include "core/geometry/tri_mesh.h"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace paracfd::core
{
	// STEP representation-item names attached to one transferred BRep face.  A face can
	// legitimately inherit more than one name (for example an ADVANCED_FACE name and the
	// enclosing SHELL_BASED_SURFACE_MODEL name), so this is deliberately a set-like list
	// rather than one lossy string.  The importer sorts and de-duplicates it.
	struct StepSourceFaceMetadata
	{
		std::uint32_t source_face_id = 0;
		std::vector<std::string> representation_names;
	};

	// Exact pre-tessellation CAD contact graph. Reciprocal source edges with the same full
	// OCCT extent are canonicalized into one physical curve before this OCC-free hand-off.
	// Each curve owns one common 3D sample chain, containing the union of its source edges'
	// tessellation samples, and every participating face records that same chain in its UV
	// coordinates. Proximity between independently generated triangles is deliberately not
	// a contact criterion.
	enum class StepContactUseKind : std::uint8_t
	{
		trim_boundary = 0,
		face_interior = 1
	};

	struct StepContactFaceUse
	{
		std::uint32_t source_face_id = 0;
		StepContactUseKind kind = StepContactUseKind::trim_boundary;
		// A boundary contributes one half-sheet sector; a curve in the interior of a
		// face contributes two. Their sum is the required discrete fan degree.
		std::uint8_t sector_count = 1;
		std::vector<std::array<double, 2>> sample_uv;
	};

	struct StepContactCurve
	{
		// Deterministic physical-curve ID. This is the smallest deterministic TopoDS edge
		// traversal ID among the exact reciprocal source edges below.
		std::uint64_t id = 0;
		// Canonical source edge retained for compatibility and concise diagnostics.
		std::uint32_t source_edge_id = TriMesh::kNoCadEdgeId;
		// All exact full-extent source-edge copies represented by this physical curve,
		// sorted by their deterministic TopoDS traversal IDs.
		std::vector<std::uint32_t> source_edge_ids;
		double tolerance_m = 0.0;
		std::uint32_t fan_degree = 0;
		// Parameters on source_edge_id's canonical OCCT curve, in the same deterministic
		// order as sample_positions_m. The order may be descending when that is required
		// to make the 3D endpoint orientation deterministic.
		std::vector<double> source_parameters;
		std::vector<std::array<double, 3>> sample_positions_m;
		// OCC-free discrete identities shared by every face-local use of each sample.
		// Exact endpoint/T-node junctions on different physical curves deliberately reuse
		// one ID; otherwise IDs identify unique graph nodes. They are stable for a
		// deterministic import and tessellation configuration.
		std::vector<std::uint32_t> sample_topology_ids;
		std::vector<StepContactFaceUse> uses;
	};

	struct StepCadContactGraph
	{
		std::vector<StepContactCurve> curves;
		// A face use whose canonical 3D contact samples could not be represented by
		// one unique, trim-valid UV chain.  The source geometry and the other contact
		// records remain useful for viewing and diagnostics, but a conforming face
		// tessellation must not consume this graph until every failure is resolved.
		struct UvProjectionFailure
		{
			std::uint64_t curve_id = TriMesh::kNoCadContactId;
			std::uint32_t source_face_id = 0;
			StepContactUseKind kind = StepContactUseKind::trim_boundary;
			std::string detail;
		};
		std::vector<UvProjectionFailure> unresolved_uv_projections;
		// Exact positive-length overlaps for which the current graph cannot yet form one
		// conforming atomic chain (for example one full CAD edge against two subedges).
		// These are never proximity guesses. A non-empty list makes the graph unsuitable
		// for constrained tessellation until the intervals are atomized or the CAD is
		// changed; callers must not silently treat the overlapping curves independently.
		struct PartialOverlap
		{
			std::uint32_t source_edge_a = TriMesh::kNoCadEdgeId;
			std::uint32_t source_edge_b = TriMesh::kNoCadEdgeId;
			std::uint32_t owner_face_a = 0;
			std::uint32_t owner_face_b = 0;
			bool edge_a_fully_covered = false;
			bool edge_b_fully_covered = false;
		};
		std::vector<PartialOverlap> unresolved_partial_overlaps;
		// Trim uses for which no full-extent oriented face-edge pcurve was
		// available. These were projected only from unique ON candidates with a
		// native-tolerance 3D round trip; the count remains visible for audits.
		std::uint32_t trim_surface_fallback_count = 0;

		bool conforming_ready() const
		{
			return unresolved_partial_overlaps.empty()
				&& unresolved_uv_projections.empty();
		}
	};

	struct StepGeometry
	{
		TriMesh mesh;
		GeometryQualityReport quality;
		StepCadContactGraph contacts;
		// Indexed in the same deterministic TopExp face traversal used by
		// TriMesh::source_face_ids.  No OpenCascade type crosses this boundary.
		std::vector<StepSourceFaceMetadata> source_faces;
	};

	// Full STEP preprocessing result.  Geometry warnings are deliberately separate from the
	// import error string: representable source-quality findings are retained in `mesh`, reported
	// as warning_run, and do not turn a successful STEP read into an error. Clients can apply the
	// shared reversible exclusion policy without modifying the STEP file.
	StepGeometry load_step_geometry(const std::string& path, double deflection_mm = 2.0,
		std::string* error = nullptr);

	// Read + triangulate a STEP file into a TriMesh (metres, face-local normals). Triangles
	// retain source-face IDs and per-half-edge CAD boundary provenance without exposing OCC
	// types. A face-boundary segment whose OCCT polygon mapping is unavailable or ambiguous is
	// explicitly marked `unknown_boundary`; it is never collapsed into the same sentinel as an
	// ordinary internal tessellation edge/opening. Face/edge IDs are stable across different tessellation tolerances of the same
	// transferred topology, but are not persistent STEP entity identifiers across CAD rewrites.
	//   deflection_mm — BRepMesh linear tolerance in MILLIMETRES (OCC's native unit; the
	//                   paraglider default is 2 mm). Smaller = finer mesh.
	// On any failure (unreadable file, empty/failed transfer, no triangulable faces) the
	// returned TriMesh is empty() and, if `error` is non-null, it holds a human-readable reason.
	TriMesh load_step_mesh(const std::string& path, double deflection_mm = 2.0, std::string* error = nullptr);

}
