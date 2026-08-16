// step_import.h — OpenCascade STEP-file import to a Qt-free triangle mesh in SI METRES.
//
// This header is deliberately free of any OpenCascade (and Qt) include so the GUI and the
// CFD preprocessing and the viewer can consume the mesh without inheriting OCC's include path or host
// flags. All OCC usage is isolated in step_import.cpp, compiled into the `paracfd_geometry`
// static lib (see CMakeLists.txt); libparacfd and the physics gates stay OCC-free.
//
// ⚠ Units: OCC emits geometry in MILLIMETRES; load_step_mesh() scales every coordinate by
// 0.001 so the returned TriMesh is in METRES, ready for the SI simulation domain and
// zero-thickness embedded-boundary preprocessing.
//
// Only triangulated TopoDS_Face entities contribute geometry. Standalone STEP edges/wires
// (for example suspension lines) are intentionally ignored, so they cannot enlarge the
// aerodynamic bbox or become impermeable CFD surfaces.
#pragma once

#include "core/geometry/tri_mesh.h"

#include <string>

namespace paracfd::core
{
	// Read + triangulate a STEP file into a TriMesh (metres, outward normals).
	//   deflection_mm — BRepMesh linear tolerance in MILLIMETRES (OCC's native unit; the
	//                   paraglider default is 2 mm). Smaller = finer mesh.
	// On any failure (unreadable file, empty/failed transfer, no triangulable faces) the
	// returned TriMesh is empty() and, if `error` is non-null, it holds a human-readable reason.
	TriMesh load_step_mesh(const std::string& path, double deflection_mm = 2.0, std::string* error = nullptr);

}
