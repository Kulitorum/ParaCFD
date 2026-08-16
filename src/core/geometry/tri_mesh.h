// tri_mesh.h — a plain, dependency-free triangle mesh in SI METRES.
//
// This struct is the hand-off contract between the OpenCascade-backed STEP importer and
// the OCC-free paraglider geometry/CFD core. It intentionally contains no OpenCascade or Qt
// types, so BVH/AMR/EB preprocessing can stay isolated from the CAD dependency.
//
// ⚠ Units: METRES. The STEP importer scales OCC's native millimetres by 0.001 before
// filling this; every consumer (display, BVH, AMR, and EB) works in metres.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace paracfd::core
{
	// A lit triangle mesh in SI METRES with face-oriented per-vertex normals. No global
	// watertight/outward shell orientation is required. Flat SoA
	// arrays: positions/normals are 3 floats per vertex (x,y,z); indices are triangle vertex
	// indices, 3 per triangle (0-based). bbox_min/max are the axis-aligned bounds in metres.
	struct TriMesh
	{
		std::vector<float> positions;       // 3 * vertex_count, metres
		std::vector<float> normals;         // 3 * vertex_count, unit, area-weighted per vertex
		// Optional face-local CAD parameters. STEP import duplicates vertices per source face,
		// therefore one (u,v) pair per mesh vertex retains the BRep tessellation provenance
		// without exposing OpenCascade types to the CFD core. Missing UV data is stored as NaN.
		std::vector<float> vertex_uv;       // 2 * vertex_count, source-face parameter coordinates
		std::vector<std::uint32_t> indices; // 3 * triangle_count
		// Stable within one imported STEP shape: the zero-based TopoDS face traversal index that
		// produced each triangle. This is the key used to accumulate CFD patches back to CAD faces.
		std::vector<std::uint32_t> source_face_ids; // triangle_count

		std::array<float, 3> bbox_min{ { 0.0f, 0.0f, 0.0f } };
		std::array<float, 3> bbox_max{ { 0.0f, 0.0f, 0.0f } };

		std::size_t vertex_count() const { return positions.size() / 3; }
		std::size_t triangle_count() const { return indices.size() / 3; }
		bool empty() const { return indices.empty(); }
		bool has_uv() const { return vertex_uv.size() == 2 * vertex_count(); }
		bool has_face_provenance() const { return source_face_ids.size() == triangle_count(); }

		std::array<float, 3> bbox_size() const
		{
			return { { bbox_max[0] - bbox_min[0], bbox_max[1] - bbox_min[1], bbox_max[2] - bbox_min[2] } };
		}
	};
}
