// mesh_selection.h — deterministic triangle subsets of an OCC-free TriMesh.
#pragma once

#include "core/geometry/tri_mesh.h"

#include <cstdint>
#include <span>
#include <vector>

namespace paracfd::core
{
	// `source_triangle_ids[q]` is the triangle in the input mesh that produced triangle q
	// in `mesh`.  Keeping this mapping explicit lets callers scatter CFD loads back to an
	// unfiltered display/source mesh when that becomes useful.
	struct TriMeshSelection
	{
		TriMesh mesh;
		std::vector<std::uint32_t> source_triangle_ids;
	};

	// Return a compact mesh containing exactly the selected input triangles, in input order.
	// Vertices unused by the selection are removed and the bbox is recomputed from referenced
	// vertices, so excluded suspension/artifact geometry cannot affect domain placement.
	TriMeshSelection select_triangles(const TriMesh& source,
		std::span<const std::uint8_t> keep_triangle);

	// Convenience wrapper for issue reports, whose triangle IDs are naturally an exclusion list.
	TriMeshSelection exclude_triangles(const TriMesh& source,
		std::span<const std::uint32_t> excluded_triangle_ids);
}
