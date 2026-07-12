// voxelize.h — watertight ray-parity voxelizer: a loaded STEP TriMesh → a solid-cell mask
// on the MAC grid, so the live flow diverts around the printed obstacle (PLAN §3).
//
// OCC-free (lives in libwindcfd): it needs only the mesh geometry, never OpenCascade. The mesh
// is placed with the SHARED core/geometry/model_placement.h transform, so the mask lands
// exactly where the viewer draws the model.
//
// Algorithm (RESEARCH/PLAN §3 "watertight ray-parity voxelizer", + owner's supersampled
// majority-fill steer):
//   * Per-sub-point inside test = axis-aligned +z ray, even-odd (parity) count of triangle
//     crossings. Robust to edge grazing via inclusive barycentric (a missed edge is caught by
//     the neighbouring triangle, keeping the crossing count's parity correct).
//   * SUPERSAMPLED MAJORITY fill: each candidate cell is probed with a KxKxK grid of
//     sub-points; the cell is solid iff its inside-fraction > 0.5. The per-cell inside-fraction
//     is RETAINED (out_fraction) — it is the cut-cell volume fraction for later use; no cut-cell
//     physics is implemented here (deferred).
//   * THIN-WALL safeguard: the >0.5 rule would delete any wall thinner than ~half a cell,
//     turning a flow barrier into a leaky one. After the majority fill, every partial cell
//     (0<fraction≤0.5) with NO solid face-neighbour — i.e. an isolated sub-cell sheet with no
//     solid backing — is kept solid so the barrier stays sealed (over-thicken beats leak).
//     The count is logged.
// Units: SI METRES; mask is the ChannelBC format (cell field, 1=solid 0=fluid, g.pidx).
#pragma once

#include "core/fluid/mac_grid.h"
#include "core/geometry/model_placement.h"
#include "core/geometry/tri_mesh.h"

#include <vector>

namespace windcfd::core
{
	// Voxelize `mesh` (placed by `place`) onto grid `g`. Returns the 1=solid/0=fluid cell mask
	// (size g.p_count(), indexed g.pidx). Optional outputs:
	//   out_mesh_volume  — |signed mesh volume| via the divergence theorem (m³).
	//   out_voxel_volume — solid_count · h³ (m³), the volume of the binary mask the sim uses.
	//   out_fraction     — per-cell inside-fraction in [0,1] (size g.p_count()); the retained
	//                      cut-cell volume fraction. Boundary cells carry their partial value.
	//   out_thin_cells   — number of thin (< h) cells kept solid by the leak safeguard.
	//   supersample (K)  — sub-points per axis for the majority test (default 3 ⇒ K³ = 27).
	// Logs a one-line summary (tris, solid cells, % of domain, mesh vs voxel volume + err) and,
	// if any, the thin-wall safeguard count.
	std::vector<unsigned char> voxelize_mesh(const TriMesh& mesh, MacGrid g, const ModelPlacement& place,
		double* out_mesh_volume = nullptr, double* out_voxel_volume = nullptr,
		std::vector<float>* out_fraction = nullptr, int* out_thin_cells = nullptr, int supersample = 3);

	// |signed mesh volume| by the divergence theorem (signed-tet sum over triangles). Placement
	// -invariant (translation does not change volume). Metres³. Exposed for the volume gate/test.
	double mesh_signed_volume_abs(const TriMesh& mesh);

	// Voxelize the SAME mesh at each of `placements` and OR the results into one solid-cell mask —
	// the settled drop/settle pile (PLAN G3) as a single rigid structure. Size g.p_count(),
	// 1=solid/0=fluid, indexed g.pidx. `out_solid_count` (optional) receives the solid-cell tally.
	// Each placement runs the full voxelize_mesh (per-instance thin-wall reseal + log); at a block–
	// block contact the two masks may over-thicken by a cell (acceptable for a rigid structure).
	std::vector<unsigned char> voxelize_mesh_instances(const TriMesh& mesh, MacGrid g,
		const std::vector<ModelPlacement>& placements, int* out_solid_count = nullptr, int supersample = 3);

	// Solidify ENCLOSED voids: mark every fluid cell (solid==0) that is NOT reachable from the open
	// domain boundary as solid. Turns the hollow interior of a voxelized building (walls + roof, with
	// trapped fluid inside) into a filled block — removing the full-cost enclosed-cavity pressure solve
	// and the spurious inner-wall load faces. A 6-connected flood from the boundary defines "outside";
	// any fluid the flood cannot reach is a sealed pocket. The GROUND face (k==0) is NOT a flood seed:
	// a building sits ON the ground, so its interior floor touches k==0 and seeding there would leak the
	// flood into the interior. Exterior ground cells are still reached via the fluid volume above them,
	// so excluding k==0 is safe (it only ever fills genuinely sealed cavities; an OPEN structure — a gap,
	// a tunnel, no roof — stays connected to the boundary and is left fluid). In-place on `solid`;
	// returns the number of cells newly filled. Host-only, OCC-free; `g` supplies dims + g.pidx only.
	int seal_enclosed_voids(std::vector<unsigned char>& solid, const MacGrid& g);
}
