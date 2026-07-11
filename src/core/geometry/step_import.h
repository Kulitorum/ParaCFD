// step_import.h — OpenCascade STEP-file import to a Qt-free triangle mesh in SI METRES.
//
// This header is deliberately free of any OpenCascade (and Qt) include so the GUI and the
// future voxelizer can consume the mesh without inheriting OCC's include path or host
// flags. All OCC usage is isolated in step_import.cpp, compiled into the `scour_geometry`
// static lib (see CMakeLists.txt); libscour and the physics gates stay OCC-free.
//
// ⚠ Units: OCC emits geometry in MILLIMETRES; load_step_mesh() scales every coordinate by
// 0.001 so the returned TriMesh is in METRES, ready for the SI simulation domain and (later)
// the watertight voxelizer.
//
// The TriMesh struct itself moved to the OCC-free core/geometry/tri_mesh.h so libscour's
// voxelizer can consume it without inheriting scour_geometry's OpenCascade dependency.
#pragma once

#include "core/geometry/tri_mesh.h"

#include <string>
#include <vector>

namespace scour::core
{
	// Read + triangulate a STEP file into a TriMesh (metres, outward normals).
	//   deflection_mm — BRepMesh linear tolerance in MILLIMETRES (OCC's native unit; the
	//                   PLAN §3 default is ~0.1 mm). Smaller = finer mesh.
	// On any failure (unreadable file, empty/failed transfer, no triangulable faces) the
	// returned TriMesh is empty() and, if `error` is non-null, it holds a human-readable reason.
	TriMesh load_step_mesh(const std::string& path, double deflection_mm = 0.1, std::string* error = nullptr);

	// Read + triangulate a STEP file into ONE mesh PER SOLID (metres, outward normals) — the
	// convex-piece decomposition for the drop/settle preparation phase (PLAN G3). Each returned
	// mesh is one convex collision piece: its vertices become a Jolt ConvexHullShape, while the
	// merged load_step_mesh() result stays the display + voxelization geometry (so the animal
	// holes are preserved for the flow — only the settling COLLISION uses the convex compound).
	// A STEP with no solids (shell/face-only) falls back to a single whole-shape mesh. On any
	// failure the returned vector is empty() and, if `error` is non-null, holds the reason.
	std::vector<TriMesh> load_step_solids(const std::string& path, double deflection_mm = 0.1, std::string* error = nullptr);
}
