// building.h — turn a 3D-printing CENTERLINE surface (vertical wall ribbons, loaded as a
// TriMesh by the STEP importer) into a solid-cell obstacle mask for the CFD, WITHOUT building
// a watertight solid. We horizontally section the centerline to a 2D footprint, then mark a
// cell solid when it lies within half the wall thickness of that footprint — optionally after
// rounding the footprint corners to a configured radius — over the wall height, plus a flat
// overhanging roof slab on top.
//
// Rationale (owner's steer): the only thing the flow needs is the solid/fluid mask, so we
// thicken the centerline directly (a distance test) instead of extruding+meshing+booleaning a
// real solid just to run inside/outside. A distance test is also inherently leak-proof (no
// thin-wall gaps). "distance <= half thickness" naturally rounds the OUTER corners by half the
// wall thickness; a larger configured radius is obtained by pre-filleting the centerline.
//
// OCC-free (libparacfd): consumes only the TriMesh + grid. The centerline STEP itself is read
// by paracfd_geometry's load_step_mesh; corner rounding here is a procedural 2D fillet.
// Units: SI METRES. Mask is the ChannelBC format (cell field, 1=solid 0=fluid, g.pidx).
#pragma once

#include "core/fluid/mac_grid.h"
#include "core/geometry/tri_mesh.h"

#include <array>
#include <vector>

namespace paracfd::core
{
	// A 2D building footprint: closed polyline loops in the xy-plane (metres), each a wall
	// centerline path. Produced by horizontally sectioning the centerline surface.
	using Loop2D = std::vector<std::array<double, 2>>;
	struct Footprint
	{
		std::vector<Loop2D> loops;
		std::array<double, 2> bbox_min{ { 0.0, 0.0 } };
		std::array<double, 2> bbox_max{ { 0.0, 0.0 } };

		bool empty() const { return loops.empty(); }
		std::array<double, 2> bbox_size() const { return { { bbox_max[0] - bbox_min[0], bbox_max[1] - bbox_min[1] } }; }
		std::size_t point_count() const;
	};

	// Parameters for building a wall+roof mask from a footprint. Distances in metres.
	struct BuildingParams
	{
		double wall_thickness = 0.08; // full wall thickness (band = +/- half of this about the centerline); COBOD-printed wall
		double wall_height = 3.00;    // wall top above base_z
		double corner_radius = 0.00;  // OUTER corner radius; 0 => "sharp" (currently min radius = half wall)
		double roof_overhang = 0.50;  // roof extends this far beyond the OUTER wall face (0 => flush)
		double roof_thickness = 0.20; // flat roof slab thickness (0 => no roof)
		double base_z = 0.00;         // ground level: domain z of the wall base
		bool   roof = true;           // false => NO roof slab: voxelize only the wall/surface band (e.g. a
		                              // wing/airfoil profile — an open extruded surface, not a capped building)
	};

	// Intersect `mesh` with the horizontal plane z = z0 (mesh coords, metres) and chain the
	// section segments into closed loops -> a 2D footprint (in mesh xy). Pass z0 = NaN to use
	// the mesh mid-height (robust to top/bottom caps). Logs a one-line summary.
	Footprint mesh_horizontal_section(const TriMesh& mesh, double z0);

	// Return a copy of `fp` translated so its xy bbox centre sits at the domain centre
	// (Lx/2, Ly/2). Base stays at z via BuildingParams.base_z. The simplest "place the building
	// in the middle of the wind tunnel" transform; a full gizmo placement can come later.
	Footprint center_footprint(const Footprint& fp, double Lx, double Ly);

	// Voxelize a building (walls thickened from the footprint + a flat overhanging roof) onto
	// grid `g`. The footprint must already be in DOMAIN xy coordinates (see center_footprint).
	// Returns the 1=solid/0=fluid cell mask (size g.p_count(), indexed g.pidx). Optional:
	//   out_solid_count — number of solid cells (walls + roof).
	// Logs a one-line summary (loops, wall cells, roof cells, % of domain).
	std::vector<unsigned char> voxelize_building(const Footprint& fp, const BuildingParams& prm,
		MacGrid g, int* out_solid_count = nullptr);
}
