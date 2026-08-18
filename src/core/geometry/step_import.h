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

#include <cstdint>
#include <string>
#include <string_view>
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

	struct StepGeometry
	{
		TriMesh mesh;
		GeometryQualityReport quality;
		// Indexed in the same deterministic TopExp face traversal used by
		// TriMesh::source_face_ids.  No OpenCascade type crosses this boundary.
		std::vector<StepSourceFaceMetadata> source_faces;
	};

	// Narrow producer-role classifier used by the import policy.  Only a normalized name
	// beginning "Mini-rib <integer>" (optionally followed by a space-delimited description)
	// matches; generic ribs and arbitrary strings containing "rib" do not. Exposed so the
	// exact policy has a cheap deterministic unit test.
	bool is_unsupported_mini_rib_representation_name(std::string_view name);

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
