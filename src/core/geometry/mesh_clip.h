// Geometry-only debug helpers. Clipping keeps the original open-sheet semantics:
// it trims triangles at the slab boundaries and never creates closing/cap faces.
#pragma once

#include "core/geometry/tri_mesh.h"

namespace paracfd::core
{
	// Retain the portion of `mesh` inside lower <= coordinate[axis] <= upper.
	// Output triangles preserve winding and source-face provenance. The returned bbox
	// spans the requested slab along `axis`, allowing a deliberately thin CFD domain
	// even when no retained vertex happens to lie exactly on one slab boundary.
	TriMesh clip_mesh_to_axis_slab(const TriMesh& mesh, int axis, double lower, double upper);
}
