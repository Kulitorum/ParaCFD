// tri_mesh.h — a plain, dependency-free triangle mesh in SI METRES.
//
// This struct is the hand-off contract between the (OpenCascade-backed) STEP importer in
// the `windcfd_geometry` lib and the voxelizer in `libwindcfd`. It lives here — in libwindcfd's
// include path, free of any OpenCascade AND Qt include — so BOTH libs can consume it
// without a libwindcfd→windcfd_geometry (or →OCC) dependency. step_import.h includes it.
//
// ⚠ Units: METRES. The STEP importer scales OCC's native millimetres by 0.001 before
// filling this; every consumer (display + voxelizer) works in metres.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace windcfd::core
{
	// A lit triangle mesh in SI METRES with outward-facing per-vertex normals. Flat SoA
	// arrays: positions/normals are 3 floats per vertex (x,y,z); indices are triangle vertex
	// indices, 3 per triangle (0-based). bbox_min/max are the axis-aligned bounds in metres.
	struct TriMesh
	{
		std::vector<float> positions;       // 3 * vertex_count, metres
		std::vector<float> normals;         // 3 * vertex_count, unit, area-weighted per vertex
		std::vector<std::uint32_t> indices; // 3 * triangle_count

		std::array<float, 3> bbox_min{ { 0.0f, 0.0f, 0.0f } };
		std::array<float, 3> bbox_max{ { 0.0f, 0.0f, 0.0f } };

		std::size_t vertex_count() const { return positions.size() / 3; }
		std::size_t triangle_count() const { return indices.size() / 3; }
		bool empty() const { return indices.empty(); }

		std::array<float, 3> bbox_size() const
		{
			return { { bbox_max[0] - bbox_min[0], bbox_max[1] - bbox_min[1], bbox_max[2] - bbox_min[2] } };
		}
	};
}
