// test_voxelize.cpp — OCC-free unit test for the ray-parity voxelizer (scour::core::
// voxelize_mesh, in libscour). Builds a unit-cube TriMesh directly in code (no OpenCascade),
// so the physics unit suite covers the volume + inside/outside gates without linking OCC.
//
// Gates (PLAN §3 / CLAUDE.md HARD RULE): |mesh_vol − voxel_vol| / mesh_vol < 2% on a grid fine
// enough that the axis-aligned cube tiles exactly; a cell at the cube centre is solid and a far
// corner is fluid; the solid-cell count scales ~×8 when h halves.
#include "core/fluid/mac_grid.h"
#include "core/geometry/model_placement.h"
#include "core/geometry/voxelize.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace
{
	using scour::core::MacGrid;
	using scour::core::ModelPlacement;
	using scour::core::TriMesh;

	// An axis-aligned cube [0,L]^3 as 12 triangles. Winding is irrelevant to ray-parity and to
	// |signed volume|, so any consistent triangulation works.
	TriMesh make_cube(float L)
	{
		TriMesh m;
		const float v[8][3] = {
			{ 0, 0, 0 }, { L, 0, 0 }, { L, L, 0 }, { 0, L, 0 },
			{ 0, 0, L }, { L, 0, L }, { L, L, L }, { 0, L, L }
		};
		for (auto& p : v) { m.positions.push_back(p[0]); m.positions.push_back(p[1]); m.positions.push_back(p[2]); }
		const unsigned q[6][4] = {
			{ 0, 3, 2, 1 }, // -z (bottom)
			{ 4, 5, 6, 7 }, // +z (top)
			{ 0, 1, 5, 4 }, // -y
			{ 2, 3, 7, 6 }, // +y
			{ 1, 2, 6, 5 }, // +x
			{ 0, 4, 7, 3 }  // -x
		};
		for (auto& f : q)
		{
			m.indices.push_back(f[0]); m.indices.push_back(f[1]); m.indices.push_back(f[2]);
			m.indices.push_back(f[0]); m.indices.push_back(f[2]); m.indices.push_back(f[3]);
		}
		m.bbox_min = { { 0, 0, 0 } };
		m.bbox_max = { { L, L, L } };
		return m;
	}
}

TEST(Voxelize, CubeVolumeGateAndInsideOutside)
{
	const float L = 0.01f; // 10 mm cube, mesh volume = 1e-6 m^3
	TriMesh cube = make_cube(L);

	// Zero placement so the cube sits at the origin, tiling cells [0,10) at h=1 mm.
	ModelPlacement place; // {0,0,0}
	MacGrid g; g.h = 0.001; g.nx = 16; g.ny = 16; g.nz = 16;

	double mesh_vol = 0.0, voxel_vol = 0.0;
	int thin = 0;
	std::vector<float> frac;
	std::vector<unsigned char> mask = scour::core::voxelize_mesh(cube, g, place, &mesh_vol, &voxel_vol, &frac, &thin);

	ASSERT_EQ((int)mask.size(), g.p_count());
	ASSERT_EQ((int)frac.size(), g.p_count());

	// Mesh volume by the divergence theorem = L^3.
	EXPECT_NEAR(mesh_vol, (double)L * L * L, 1e-12);

	// <2% volume gate (axis-aligned cube tiles exactly → expect ~0%).
	const double err = std::abs(mesh_vol - voxel_vol) / mesh_vol;
	EXPECT_LT(err, 0.02) << "mesh_vol=" << mesh_vol << " voxel_vol=" << voxel_vol;

	// A cell at the cube centre (~0.005 m) is solid; a far corner is fluid.
	EXPECT_TRUE(mask[g.pidx(5, 5, 5)] != 0) << "cube-centre cell should be solid";
	EXPECT_TRUE(mask[g.pidx(15, 15, 15)] == 0) << "far-corner cell should be fluid";

	// The retained inside-fraction is 1.0 for a fully-interior cell.
	EXPECT_NEAR(frac[g.pidx(5, 5, 5)], 1.0f, 1e-6f);

	int solid = 0;
	for (unsigned char c : mask) solid += (c != 0);
	EXPECT_EQ(solid, 1000) << "10x10x10 cells should tile the 10 mm cube at h=1 mm";
}

TEST(Voxelize, SolidCountScalesByEightWhenResolutionDoubles)
{
	const float L = 0.01f;
	TriMesh cube = make_cube(L);
	ModelPlacement place;

	auto count_solid = [&](double h, int n) -> int
	{
		MacGrid g; g.h = h; g.nx = n; g.ny = n; g.nz = n;
		std::vector<unsigned char> mask = scour::core::voxelize_mesh(cube, g, place);
		int solid = 0;
		for (unsigned char c : mask) solid += (c != 0);
		return solid;
	};

	const int coarse = count_solid(0.001, 16);   // 10 cells across → 1000
	const int fine = count_solid(0.0005, 32);    // 20 cells across → 8000
	EXPECT_GT(coarse, 0);
	EXPECT_NEAR((double)fine / (double)coarse, 8.0, 0.05);
}

// voxelize_mesh_instances — the union voxelizer that turns a settled drop/settle pile (many
// placements of one mesh) into a single rigid-structure mask (PLAN G3.2).
TEST(Voxelize, InstancesUnionDisjointOverlappingAndEmpty)
{
	const float L = 0.01f; // 10 mm cube
	TriMesh cube = make_cube(L);
	MacGrid g; g.h = 0.001; g.nx = 32; g.ny = 32; g.nz = 32;

	// Two cubes placed apart (1 mm gap) → disjoint → union = 2 × 1000 cells.
	ModelPlacement a;                // origin: cells x ∈ [0,10)
	ModelPlacement b; b.tx = 0.011;  // +11 mm: cells x ∈ [11,21)
	int count = -1;
	std::vector<unsigned char> um = scour::core::voxelize_mesh_instances(cube, g, { a, b }, &count);
	ASSERT_EQ((int)um.size(), g.p_count());
	EXPECT_EQ(count, 2000);
	int solid = 0;
	for (unsigned char c : um) solid += (c != 0);
	EXPECT_EQ(solid, count) << "out_solid_count must match the mask tally";
	EXPECT_TRUE(um[g.pidx(5, 5, 5)] != 0) << "first cube solid";
	EXPECT_TRUE(um[g.pidx(15, 5, 5)] != 0) << "second cube solid";
	EXPECT_TRUE(um[g.pidx(10, 5, 5)] == 0) << "gap between the two cubes is fluid";

	// Two IDENTICAL placements → OR is idempotent → union = one cube = 1000 cells.
	int count2 = -1;
	scour::core::voxelize_mesh_instances(cube, g, { a, a }, &count2);
	EXPECT_EQ(count2, 1000);

	// No placements → empty (all-fluid) mask of the right size.
	int count3 = -1;
	std::vector<unsigned char> em = scour::core::voxelize_mesh_instances(cube, g, {}, &count3);
	EXPECT_EQ(count3, 0);
	EXPECT_EQ((int)em.size(), g.p_count());
}
