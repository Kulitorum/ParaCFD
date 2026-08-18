// Regression probe for fail-closed EB construction transactions.
//
// This deliberately keeps the tests outside embedded_boundary.cpp:
//  1. an exact decomposer rejects only a covered/unowned fine-level halo cell;
//     the composite atlas must still reject the build;
//  2. one exact shared face commits, then the next shared face rejects; no topology
//     record belonging to either rejected cell may survive the transaction.

#include "core/fluid/amr_eb.h"
#include "core/fluid/amr_grid.h"
#include "core/geometry/exact_cell_decomposition.h"
#include "core/geometry/triangle_bvh.h"
#include "core/paraglider_config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

using namespace paracfd::core;

namespace
{
	int failures = 0;

	void check(bool condition, const char* description)
	{
		std::printf("[eb-transaction] %-68s %s\n", description,
			condition ? "PASS" : "FAIL");
		if (!condition) ++failures;
	}

	TriMesh quad(Vec3d a, Vec3d b, Vec3d c, Vec3d d)
	{
		TriMesh mesh;
		mesh.bbox_min = {1e30f, 1e30f, 1e30f};
		mesh.bbox_max = {-1e30f, -1e30f, -1e30f};
		for (Vec3d point : {a, b, c, d})
		{
			mesh.positions.insert(mesh.positions.end(), {static_cast<float>(point.x),
				static_cast<float>(point.y), static_cast<float>(point.z)});
			mesh.positions_fp64.insert(mesh.positions_fp64.end(), {point.x, point.y, point.z});
			for (int axis = 0; axis < 3; ++axis)
			{
				mesh.bbox_min[axis] = std::min(mesh.bbox_min[axis], static_cast<float>(point[axis]));
				mesh.bbox_max[axis] = std::max(mesh.bbox_max[axis], static_cast<float>(point[axis]));
			}
		}
		mesh.indices = {0, 1, 2, 0, 2, 3};
		mesh.source_face_ids = {0, 0};
		return mesh;
	}

	ExactCellDecomposition synthetic_open_cell(const ExactCellInput& input,
		bool omit_lower_x_aperture);

	ExactCellDecomposition reject_exact_halo_cell(const ExactCellInput&)
	{
		ExactCellDecomposition rejected;
		rejected.errors.push_back("intentional exact-decomposer rejection in unowned topology halo");
		return rejected;
	}

	ExactCellPoint point(double x, double y, double z)
	{
		return {x, y, z};
	}

	void append_oriented_quad(std::vector<ExactCellOrientedTriangle>& triangles,
		ExactCellPoint a, ExactCellPoint b, ExactCellPoint c, ExactCellPoint d)
	{
		triangles.push_back({{a, b, c}});
		triangles.push_back({{a, c, d}});
	}

	ExactCellPlanarRegion face_region(ExactCellPoint a, ExactCellPoint b,
		ExactCellPoint c, ExactCellPoint d, double h)
	{
		ExactCellPlanarRegion region;
		region.outer_loop = {a, b, c, d};
		region.convex_pieces = {{{a, b, c}}, {{a, c, d}}};
		region.convex_piece_area_sum = h * h;
		region.convex_piece_first_moment = {
			0.25 * h * h * (a[0] + b[0] + c[0] + d[0]),
			0.25 * h * h * (a[1] + b[1] + c[1] + d[1]),
			0.25 * h * h * (a[2] + b[2] + c[2] + d[2])};
		return region;
	}

	ExactCellDecomposition synthetic_open_cell(const ExactCellInput& input,
		bool omit_lower_x_aperture)
	{
		const double x0 = input.cell_min[0], y0 = input.cell_min[1], z0 = input.cell_min[2];
		const double x1 = input.cell_max[0], y1 = input.cell_max[1], z1 = input.cell_max[2];
		const double h = x1 - x0;
		const ExactCellPoint p000 = point(x0, y0, z0), p001 = point(x0, y0, z1);
		const ExactCellPoint p010 = point(x0, y1, z0), p011 = point(x0, y1, z1);
		const ExactCellPoint p100 = point(x1, y0, z0), p101 = point(x1, y0, z1);
		const ExactCellPoint p110 = point(x1, y1, z0), p111 = point(x1, y1, z1);

		ExactCellDecomposition result;
		ExactCellFragment fragment;
		fragment.id = 0;
		fragment.volume = h * h * h;
		fragment.centroid = fragment.interior_witness = point(
			0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1));
		// CCW as seen from outside the box.
		append_oriented_quad(fragment.boundary_triangles, p000, p001, p011, p010); // -X
		append_oriented_quad(fragment.boundary_triangles, p100, p110, p111, p101); // +X
		append_oriented_quad(fragment.boundary_triangles, p000, p100, p101, p001); // -Y
		append_oriented_quad(fragment.boundary_triangles, p010, p011, p111, p110); // +Y
		append_oriented_quad(fragment.boundary_triangles, p000, p010, p110, p100); // -Z
		append_oriented_quad(fragment.boundary_triangles, p001, p101, p111, p011); // +Z
		result.fragments.push_back(std::move(fragment));

		auto aperture = [&](int axis, bool upper, ExactCellPoint a, ExactCellPoint b,
			ExactCellPoint c, ExactCellPoint d)
		{
			if (axis == 0 && !upper && omit_lower_x_aperture) return;
			ExactCellBoxAperture value;
			value.axis = static_cast<std::int8_t>(axis);
			value.upper = upper;
			value.area = h * h;
			value.centroid = point(0.25 * (a[0] + b[0] + c[0] + d[0]),
				0.25 * (a[1] + b[1] + c[1] + d[1]),
				0.25 * (a[2] + b[2] + c[2] + d[2]));
			value.region = face_region(a, b, c, d, h);
			value.fragment = 0;
			result.box_apertures.push_back(std::move(value));
		};
		aperture(0, false, p000, p001, p011, p010);
		aperture(0, true,  p100, p110, p111, p101);
		aperture(1, false, p000, p100, p101, p001);
		aperture(1, true,  p010, p011, p111, p110);
		aperture(2, false, p000, p010, p110, p100);
		aperture(2, true,  p001, p101, p111, p011);
		result.expected_cell_volume = result.fragment_volume_sum = h * h * h;
		return result;
	}

	ExactCellDecomposition reject_third_cell_lower_face(const ExactCellInput& input)
	{
		return synthetic_open_cell(input, input.cell_min[0] >= 2.0 - 1e-12);
	}

	int referenced_cell(const EmbeddedBoundary& eb, FragmentRef ref)
	{
		if (ref == invalid_fragment) return -1;
		if (fragment_is_regular(ref)) return regular_fragment_cell(ref);
		const int fragment = irregular_fragment_index(ref);
		return fragment >= 0 && fragment < static_cast<int>(eb.fragments.size())
			? eb.fragments[fragment].parent_cell : -1;
	}

	bool record_references_unresolved(const EmbeddedBoundary& eb, FragmentRef ref)
	{
		const int cell = referenced_cell(eb, ref);
		return cell >= 0 && cell < static_cast<int>(eb.cells.size())
			&& eb.cells[cell].state == EbCellState::unresolved;
	}

	bool rejected_cells_are_transactionally_clean(const EmbeddedBoundary& eb)
	{
		for (int cell : eb.irregular_cells)
			if (cell >= 0 && cell < static_cast<int>(eb.cells.size())
				&& eb.cells[cell].state == EbCellState::unresolved) return false;
		for (const FluidFragment& fragment : eb.fragments)
			if (fragment.parent_cell >= 0 && fragment.parent_cell < static_cast<int>(eb.cells.size())
				&& eb.cells[fragment.parent_cell].state == EbCellState::unresolved) return false;
		for (const FragmentConnection& connection : eb.connections)
			if (record_references_unresolved(eb, connection.fragment_a)
				|| record_references_unresolved(eb, connection.fragment_b)) return false;
		for (const FaceAperture& aperture : eb.apertures)
			if (record_references_unresolved(eb, aperture.fragment_a)
				|| record_references_unresolved(eb, aperture.fragment_b)) return false;
		for (const SurfacePatch& patch : eb.patches)
			if (record_references_unresolved(eb, patch.plus_fragment)
				|| record_references_unresolved(eb, patch.minus_fragment)) return false;
		for (int cell = 0; cell < static_cast<int>(eb.cells.size()); ++cell)
			if (eb.cells[cell].state == EbCellState::unresolved
				&& (eb.cells[cell].fragment_count != 0 || eb.cells[cell].arrangement_index >= 0
					|| eb.cells[cell].exact_locator_index >= 0)) return false;
		return true;
	}

	void halo_rejection_test()
	{
		// Build one exact-rejected topology cell, then place that cell in an atlas as
		// a covered topology halo. This isolates the atlas contract from the separate
		// coarse/fine aperture policy: ownership must never erase a geometry failure.
		const TriMesh topology_mesh = quad({0.2, 0.2, 0.5}, {0.8, 0.2, 0.5},
			{0.8, 0.8, 0.5}, {0.2, 0.8, 0.5});
		const TriangleBvh topology_bvh(topology_mesh);
		EmbeddedBoundaryBuildOptions options;
		options.exact_cell_decomposer = &reject_exact_halo_cell;
		options.complex_subdivisions = 8;
		options.min_volume_fraction = 1e-12;
		AmrEbLevelAtlas level;
		level.level = 1;
		level.topology = build_embedded_boundary(topology_mesh, topology_bvh,
			{{0, 0, 0}, 1, 1, 1, 1.0}, options);
		level.owned_cell.assign(1, 0);
		AmrEmbeddedBoundaryAtlas atlas;
		atlas.levels.push_back(std::move(level));

		std::size_t owned_unresolved = 0, halo_unresolved = 0;
		for (const AmrEbLevelAtlas& level : atlas.levels)
			for (const UnresolvedEbCell& problem : level.topology.unresolved)
			{
				const bool owned = problem.parent_cell >= 0
					&& problem.parent_cell < static_cast<int>(level.owned_cell.size())
					&& level.owned_cell[problem.parent_cell] != 0;
				owned ? ++owned_unresolved : ++halo_unresolved;
			}
		check(owned_unresolved == 0 && halo_unresolved > 0,
			"fixture rejects exact topology only in an unowned fine-level halo");
		check(!atlas.ready_for_flow(),
			"an unresolved topology halo makes the composite atlas non-flow-ready");
	}

	void shared_face_rollback_test()
	{
		const TriMesh mesh = quad({0, 0, 0.5}, {3, 0, 0.5},
			{3, 1, 0.5}, {0, 1, 0.5});
		const TriangleBvh bvh(mesh);
		EmbeddedBoundaryBuildOptions options;
		options.exact_cell_decomposer = &reject_third_cell_lower_face;
		options.min_volume_fraction = 1e-12;
		const EmbeddedBoundary eb = build_embedded_boundary(mesh, bvh,
			{{0, 0, 0}, 3, 1, 1, 1.0}, options);

		check(!eb.ready_for_flow() && eb.cells[1].state == EbCellState::unresolved
			&& eb.cells[2].state == EbCellState::unresolved,
			"second shared-face rejection invalidates both participating cells");
		check(rejected_cells_are_transactionally_clean(eb),
			"rejected cells retain no fragments, patches, apertures, or connections");
	}
}

int main()
{
	halo_rejection_test();
	shared_face_rollback_test();
	std::printf("[eb-transaction] failures=%d\n", failures);
	return failures == 0 ? 0 : 1;
}
