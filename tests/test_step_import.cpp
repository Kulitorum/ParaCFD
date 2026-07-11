// test_step_import.cpp — unit test for the OpenCascade STEP → TriMesh loader
// (scour::core::load_step_mesh, in the scour_geometry lib). Loads the bundled
// tests/inputs/cube.step — a 10 mm cube (0,0,0)→(10,10,10) in the STEP file — and asserts
// the mesh is non-empty and its bounding box is a 0.01 m cube (the mm→m scaling trap).
#include "core/geometry/step_import.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

namespace
{
	// Path to the committed asset. STEP_CUBE_PATH is injected by CMake (absolute); fall back
	// to a relative path so the exe is also runnable straight from the build dir if needed.
#ifndef STEP_CUBE_PATH
#define STEP_CUBE_PATH "tests/inputs/cube.step"
#endif

	std::string cube_path() { return std::string(STEP_CUBE_PATH); }
}

TEST(StepImport, CubeLoadsAsMetreScaleMesh)
{
	std::string err;
	scour::core::TriMesh mesh = scour::core::load_step_mesh(cube_path(), 0.1, &err);

	ASSERT_FALSE(mesh.empty()) << "load_step_mesh failed: " << err << " (path=" << cube_path() << ")";
	EXPECT_GT(mesh.triangle_count(), 0u);
	EXPECT_EQ(mesh.indices.size() % 3, 0u);
	EXPECT_EQ(mesh.positions.size(), mesh.normals.size());

	// A closed box triangulates to at least 12 triangles (2 per face * 6 faces).
	EXPECT_GE(mesh.triangle_count(), 12u);

	// The STEP cube spans 0..10 mm on each axis → 0..0.01 m after the mm→m scaling.
	const auto size = mesh.bbox_size();
	const float expected = 0.010f; // metres
	const float tol = 1e-4f;
	EXPECT_NEAR(size[0], expected, tol) << "x extent (m)";
	EXPECT_NEAR(size[1], expected, tol) << "y extent (m)";
	EXPECT_NEAR(size[2], expected, tol) << "z extent (m)";

	EXPECT_NEAR(mesh.bbox_min[0], 0.0f, tol);
	EXPECT_NEAR(mesh.bbox_min[1], 0.0f, tol);
	EXPECT_NEAR(mesh.bbox_min[2], 0.0f, tol);
	EXPECT_NEAR(mesh.bbox_max[0], expected, tol);
	EXPECT_NEAR(mesh.bbox_max[1], expected, tol);
	EXPECT_NEAR(mesh.bbox_max[2], expected, tol);

	// Every per-vertex normal must be unit length (area-weighted then normalised).
	for (std::size_t v = 0; v < mesh.vertex_count(); ++v)
	{
		const float nx = mesh.normals[3 * v + 0];
		const float ny = mesh.normals[3 * v + 1];
		const float nz = mesh.normals[3 * v + 2];
		EXPECT_NEAR(std::sqrt(nx * nx + ny * ny + nz * nz), 1.0f, 1e-3f);
	}
}

TEST(StepImport, MissingFileReturnsEmptyWithError)
{
	std::string err;
	scour::core::TriMesh mesh = scour::core::load_step_mesh("does_not_exist_12345.step", 0.1, &err);
	EXPECT_TRUE(mesh.empty());
	EXPECT_FALSE(err.empty());
}

// load_step_solids() — the per-solid convex-piece decomposition for the G3 drop/settle phase.
// The cube is ONE closed solid, so it must come back as exactly one piece whose geometry matches
// the whole-shape load_step_mesh() result (same triangle count, same 0.01 m box).
TEST(StepImport, CubeLoadsAsSingleSolid)
{
	std::string err;
	std::vector<scour::core::TriMesh> solids = scour::core::load_step_solids(cube_path(), 0.1, &err);

	ASSERT_FALSE(solids.empty()) << "load_step_solids failed: " << err << " (path=" << cube_path() << ")";
	EXPECT_EQ(solids.size(), 1u) << "cube is a single solid ⇒ one convex piece";

	const scour::core::TriMesh& s = solids.front();
	EXPECT_GE(s.triangle_count(), 12u);
	EXPECT_EQ(s.positions.size(), s.normals.size());

	scour::core::TriMesh whole = scour::core::load_step_mesh(cube_path(), 0.1, &err);
	ASSERT_FALSE(whole.empty());
	EXPECT_EQ(s.triangle_count(), whole.triangle_count()) << "single-solid == whole-shape geometry";

	const auto size = s.bbox_size();
	const float tol = 1e-4f;
	EXPECT_NEAR(size[0], 0.010f, tol);
	EXPECT_NEAR(size[1], 0.010f, tol);
	EXPECT_NEAR(size[2], 0.010f, tol);
}

TEST(StepImport, SolidsMissingFileReturnsEmptyWithError)
{
	std::string err;
	std::vector<scour::core::TriMesh> solids = scour::core::load_step_solids("does_not_exist_12345.step", 0.1, &err);
	EXPECT_TRUE(solids.empty());
	EXPECT_FALSE(err.empty());
}
