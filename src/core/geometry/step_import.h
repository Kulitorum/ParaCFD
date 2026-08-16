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
// The TriMesh struct itself moved to the OCC-free core/geometry/tri_mesh.h so libparacfd's
// voxelizer can consume it without inheriting paracfd_geometry's OpenCascade dependency.
#pragma once

#include "core/geometry/tri_mesh.h"

#include <string>
#include <vector>

namespace paracfd::core
{
	// Read + triangulate a STEP file into a TriMesh (metres, outward normals).
	//   deflection_mm — BRepMesh linear tolerance in MILLIMETRES (OCC's native unit; the
	//                   paraglider default is 2 mm). Smaller = finer mesh.
	// On any failure (unreadable file, empty/failed transfer, no triangulable faces) the
	// returned TriMesh is empty() and, if `error` is non-null, it holds a human-readable reason.
	TriMesh load_step_mesh(const std::string& path, double deflection_mm = 2.0, std::string* error = nullptr);

	// Same as load_step_mesh but sourced from the raw STEP file BYTES held in memory (the exact
	// contents of the .stp), not a path. A saved scene embeds the STEP (the source of truth) and
	// reconstructs the triangulation from it on load — the mesh is a derived artifact, so we store
	// the STEP and regenerate, never the other way round. Empty() on failure (as load_step_mesh).
	TriMesh load_step_mesh_from_memory(const std::vector<unsigned char>& step_bytes, double deflection_mm = 2.0, std::string* error = nullptr);

	// Read + triangulate a STEP file into ONE mesh PER SOLID (metres, outward normals) — the
	// convex-piece decomposition for the drop/settle preparation phase (PLAN G3). Each returned
	// mesh is one convex collision piece: its vertices become a Jolt ConvexHullShape, while the
	// merged load_step_mesh() result stays the display + voxelization geometry (so the animal
	// holes are preserved for the flow — only the settling COLLISION uses the convex compound).
	// A STEP with no solids (shell/face-only) falls back to a single whole-shape mesh. On any
	// failure the returned vector is empty() and, if `error` is non-null, holds the reason.
	std::vector<TriMesh> load_step_solids(const std::string& path, double deflection_mm = 0.1, std::string* error = nullptr);
}
