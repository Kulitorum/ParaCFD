// Geometry-only debug helpers. Clipping keeps the original open-sheet semantics:
// it trims triangles at the slab boundaries and never creates closing/cap faces.
#pragma once

#include "core/geometry/tri_mesh.h"

namespace paracfd::core
{
	// Retain the portion of `mesh` inside lower <= coordinate[axis] <= upper.
	// Output triangles preserve winding, exact CPU positions, and source-face provenance.
	// When the input has discrete topology IDs, retained vertices keep them and a cut
	// point shared by copies of the same topological edge receives one deterministic new
	// ID. Distinct topological edges are never welded merely because they coincide.
	// A retained subsegment of an original CAD half-edge keeps its CAD edge ID/count;
	// clip-plane edges and triangulation fan diagonals are explicitly marked unavailable.
	// The returned bbox
	// spans the requested slab along `axis`, allowing a deliberately thin CFD domain
	// even when no retained vertex happens to lie exactly on one slab boundary.
	// A mirror-domain crop may discard an original surface lying wholly in either
	// symmetry plane: the physical mirror condition already supplies its no-normal-flow
	// boundary, and the omitted opposite-side pressure volume is outside the domain.
	// Crossing sheets are still clipped at the plane and remain open; no cap is created.
	TriMesh clip_mesh_to_axis_slab(const TriMesh& mesh, int axis, double lower, double upper,
		bool discard_lower_coplanar = false, bool discard_upper_coplanar = false);
}
