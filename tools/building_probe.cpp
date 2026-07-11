// building_probe.cpp — dev tool: load a 3D-printing centerline STEP, section it to a 2D
// footprint, and voxelize a thickened wall + flat roof building. Prints stats so the geometry
// pipeline can be verified headlessly before it is wired into the GUI.
//
// usage: building_probe <centerline.stp> [wall_thickness] [wall_height] [corner_radius]
//                       [roof_overhang] [roof_thickness] [voxel_h]   (all metres)
#include "core/fluid/mac_grid.h"
#include "core/geometry/building.h"
#include "core/geometry/step_import.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace windcfd::core;

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::printf("usage: building_probe <centerline.stp> [wall_thickness] [wall_height] "
					"[corner_radius] [roof_overhang] [roof_thickness] [voxel_h]  (metres)\n");
		return 2;
	}
	const std::string path = argv[1];
	BuildingParams prm;
	if (argc > 2) prm.wall_thickness = std::atof(argv[2]);
	if (argc > 3) prm.wall_height = std::atof(argv[3]);
	if (argc > 4) prm.corner_radius = std::atof(argv[4]);
	if (argc > 5) prm.roof_overhang = std::atof(argv[5]);
	if (argc > 6) prm.roof_thickness = std::atof(argv[6]);
	double h = (argc > 7) ? std::atof(argv[7]) : 0.10;

	std::string err;
	TriMesh mesh = load_step_mesh(path, 0.5, &err);
	if (mesh.empty())
	{
		std::printf("[probe] load failed: %s\n", err.c_str());
		return 1;
	}
	std::printf("[probe] loaded %zu tris, bbox=[%.3f %.3f %.3f]..[%.3f %.3f %.3f] m\n",
		mesh.triangle_count(), (double)mesh.bbox_min[0], (double)mesh.bbox_min[1], (double)mesh.bbox_min[2],
		(double)mesh.bbox_max[0], (double)mesh.bbox_max[1], (double)mesh.bbox_max[2]);

	Footprint fp = mesh_horizontal_section(mesh, std::nan(""));
	if (fp.empty())
	{
		std::printf("[probe] empty footprint (no closed loops from the section)\n");
		return 1;
	}

	// Domain sized around the building: footprint + generous margins, base at z=0.
	const auto sz = fp.bbox_size();
	const double margin = std::max(sz[0], sz[1]) * 1.5 + 2.0;
	const double Lx = sz[0] + 2.0 * margin, Ly = sz[1] + 2.0 * margin;
	const double Lz = prm.base_z + prm.wall_height + prm.roof_thickness + 3.0;
	MacGrid g;
	g.h = h;
	g.nx = (int)std::ceil(Lx / h);
	g.ny = (int)std::ceil(Ly / h);
	g.nz = (int)std::ceil(Lz / h);
	std::printf("[probe] domain %.1f x %.1f x %.1f m @ h=%.3f -> %d x %d x %d = %d cells\n",
		Lx, Ly, Lz, h, g.nx, g.ny, g.nz, g.p_count());

	Footprint placed = center_footprint(fp, g.nx * h, g.ny * h);
	int solid = 0;
	std::vector<unsigned char> mask = voxelize_building(placed, prm, g, &solid);
	std::printf("[probe] DONE: %d solid cells (%.2f%% of domain)\n",
		solid, 100.0 * solid / std::max(1, g.p_count()));
	return solid > 0 ? 0 : 3;
}
