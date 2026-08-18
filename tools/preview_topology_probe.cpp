#include "core/geometry/embedded_boundary.h"
#include "core/geometry/triangle_bvh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <queue>
#include <vector>

using namespace paracfd::core;

namespace
{
	int failures = 0;

	void check(bool condition, const char* description)
	{
		std::printf("[preview-topology] %-68s %s\n", description,
			condition ? "PASS" : "FAIL");
		if (!condition)
			++failures;
	}

	bool near(double a, double b, double tolerance = 1e-12)
	{
		return std::abs(a - b) <= tolerance *
			std::max({ 1.0, std::abs(a), std::abs(b) });
	}

	TriMesh mesh_from_quads(const std::vector<std::array<Vec3d, 4>>& quads)
	{
		TriMesh mesh;
		mesh.bbox_min = { 1e30f, 1e30f, 1e30f };
		mesh.bbox_max = { -1e30f, -1e30f, -1e30f };
		for (std::size_t face = 0; face < quads.size(); ++face)
		{
			const std::uint32_t base = static_cast<std::uint32_t>(mesh.vertex_count());
			for (const Vec3d point : quads[face])
			{
				mesh.positions.insert(mesh.positions.end(), {
					static_cast<float>(point.x), static_cast<float>(point.y),
					static_cast<float>(point.z) });
				mesh.positions_fp64.insert(mesh.positions_fp64.end(), {
					point.x, point.y, point.z });
				for (int axis = 0; axis < 3; ++axis)
				{
					mesh.bbox_min[axis] = std::min(mesh.bbox_min[axis],
						static_cast<float>(point[axis]));
					mesh.bbox_max[axis] = std::max(mesh.bbox_max[axis],
						static_cast<float>(point[axis]));
				}
			}
			mesh.indices.insert(mesh.indices.end(), {
				base, base + 1, base + 2, base, base + 2, base + 3 });
			mesh.source_face_ids.push_back(static_cast<std::uint32_t>(face));
			mesh.source_face_ids.push_back(static_cast<std::uint32_t>(face));
		}
		return mesh;
	}

	EmbeddedBoundaryBuildOptions preview_options()
	{
		EmbeddedBoundaryBuildOptions options;
		// These are the current interactive GUI preview settings. Keep this probe
		// explicit: it freezes the legacy behavior while that path is quarantined.
		options.complex_subdivisions = 4;
		options.min_volume_fraction = 0.005;
		options.min_aperture_area_fraction = 1e-4;
		options.allow_unverified_same_fragment_patches_for_diagnostics = true;
		options.exact_cell_decomposer = nullptr;
		return options;
	}

	EmbeddedBoundary build_preview(const TriMesh& mesh,
		UniformEbGrid grid = { { 0, 0, 0 }, 1, 1, 1, 1.0 })
	{
		const TriangleBvh bvh(mesh);
		return build_embedded_boundary(mesh, bvh, grid, preview_options());
	}

	bool connected(const EmbeddedBoundary& boundary, FragmentRef begin, FragmentRef target)
	{
		if (begin == invalid_fragment || target == invalid_fragment)
			return false;
		std::queue<FragmentRef> pending;
		std::vector<FragmentRef> visited;
		pending.push(begin);
		while (!pending.empty())
		{
			const FragmentRef current = pending.front();
			pending.pop();
			if (current == target)
				return true;
			if (std::find(visited.begin(), visited.end(), current) != visited.end())
				continue;
			visited.push_back(current);
			for (const FragmentConnection& connection : boundary.connections)
			{
				if (connection.fragment_a == current)
					pending.push(connection.fragment_b);
				else if (connection.fragment_b == current)
					pending.push(connection.fragment_a);
			}
		}
		return false;
	}

	double patch_area(const EmbeddedBoundary& boundary)
	{
		double area = 0.0;
		for (const SurfacePatch& patch : boundary.patches)
			area += patch.area;
		return area;
	}

	void report(const char* name, const EmbeddedBoundary& boundary)
	{
		std::printf("[preview-topology] %s: ready=%d fragments=%zu apertures=%zu "
			"connections=%zu patches=%zu arrangements=%zu sampled=%zu "
			"accepted-same=%zu diagnostic-same=%zu unmapped=%zu rejected-same=%zu\n",
			name, boundary.ready_for_flow() ? 1 : 0, boundary.fragments.size(),
			boundary.apertures.size(), boundary.connections.size(), boundary.patches.size(),
			boundary.arrangements.size(), boundary.sampled_voxel_fragments.size(),
			boundary.accepted_free_edge_same_fragment_patches,
			boundary.diagnostic_unverified_same_fragment_patches,
			boundary.diagnostic_unmapped_surface_patches,
			boundary.rejected_same_fragment_patches);
	}

	void opening_case()
	{
		// Two pieces of one membrane leave a resolved 0.5 h opening across the middle.
		// Its width deliberately spans two of the preview's four sampling rows.
		const TriMesh mesh = mesh_from_quads({
			{{ { 0.5, 0, 0 }, { 0.5, 0.25, 0 }, { 0.5, 0.25, 1 }, { 0.5, 0, 1 } }},
			{{ { 0.5, 0.75, 0 }, { 0.5, 1, 0 }, { 0.5, 1, 1 }, { 0.5, 0.75, 1 } }} });
		const EmbeddedBoundary boundary = build_preview(mesh,
			UniformEbGrid{ { 0, 0, 0 }, 2, 1, 1, 1.0 });
		report("opening", boundary);

		const FragmentRef left = boundary.fragment_containing_point(0, { 0.25, 0.5, 0.5 });
		const FragmentRef right = boundary.fragment_containing_point(1, { 1.5, 0.5, 0.5 });
		check(boundary.ready_for_flow(), "opening topology is flow-ready");
		check(boundary.fragments.size() == 1 && boundary.patches.size() == 4,
			"opening retains one fluid fragment and four triangle patches");
		check(boundary.apertures.size() == 1 && boundary.connections.size() == 1 &&
			near(boundary.apertures.front().area, 1.0),
			"opening cell exports one full-area aperture to its regular neighbour");
		check(near(patch_area(boundary), 0.5), "opening retains exactly 0.5 h^2 of fabric");
		check(left != invalid_fragment && connected(boundary, left, right),
			"fluid connects from one side to the other only through the gap");
		check(boundary.accepted_free_edge_same_fragment_patches == boundary.patches.size(),
			"confirmed free edges account for every same-fragment opening patch");
	}

	void t_junction_case()
	{
		const TriMesh mesh = mesh_from_quads({
			{{ { 0, 0, 0.5 }, { 1, 0, 0.5 }, { 1, 1, 0.5 }, { 0, 1, 0.5 } }},
			{{ { 0.5, 0, 0 }, { 0.5, 1, 0 }, { 0.5, 1, 0.5 }, { 0.5, 0, 0.5 } }} });
		const EmbeddedBoundary boundary = build_preview(mesh);
		report("T-junction", boundary);

		const FragmentRef below_left = boundary.fragment_containing_point(0, { 0.25, 0.25, 0.25 });
		const FragmentRef below_right = boundary.fragment_containing_point(0, { 0.75, 0.25, 0.25 });
		const FragmentRef above_left = boundary.fragment_containing_point(0, { 0.25, 0.25, 0.75 });
		const FragmentRef above_right = boundary.fragment_containing_point(0, { 0.75, 0.25, 0.75 });
		check(boundary.ready_for_flow(), "T-junction topology is flow-ready");
		check(boundary.fragments.size() == 3 && boundary.patches.size() == 4,
			"T-junction retains three fluid sectors and four triangle patches");
		check(boundary.apertures.empty() && boundary.connections.empty() &&
			boundary.sampled_voxel_fragments.size() == 64,
			"single-cell T-junction has zero inter-cell apertures and a 4^3 preview map");
		check(near(patch_area(boundary), 1.5), "T-junction fabric area is conservative");
		check(below_left != below_right && !connected(boundary, below_left, below_right),
			"lower sectors have no connection through the vertical fabric");
		check(above_left == above_right && below_left != above_left && below_right != above_left,
			"upper sector is connected while both lower fabric sides stay independent");
	}

	void closed_wedge_case()
	{
		const TriMesh mesh = mesh_from_quads({
			{{ { 0, 0, 0.75 }, { 0.5, 0, 0.5 }, { 0.5, 1, 0.5 }, { 0, 1, 0.75 } }},
			{{ { 0, 0, 0.25 }, { 0, 1, 0.25 }, { 0.5, 1, 0.5 }, { 0.5, 0, 0.5 } }} });
		const EmbeddedBoundary boundary = build_preview(mesh);
		report("closed wedge", boundary);

		// Query represented preview-sample locations instead of the cusp ambiguity band.
		const FragmentRef interior = boundary.fragment_containing_point(0, { 0.125, 0.5, 0.375 });
		const FragmentRef exterior = boundary.fragment_containing_point(0, { 0.875, 0.5, 0.375 });
		bool distinct_patch_sides = true;
		for (const SurfacePatch& patch : boundary.patches)
			distinct_patch_sides = distinct_patch_sides &&
				patch.plus_fragment != invalid_fragment &&
				patch.minus_fragment != invalid_fragment &&
				patch.plus_fragment != patch.minus_fragment;
		check(boundary.ready_for_flow(), "closed wedge topology is flow-ready");
		check(boundary.fragments.size() == 2 && boundary.patches.size() == 4,
			"closed wedge retains two fluid fragments and four triangle patches");
		check(boundary.apertures.empty() && boundary.connections.empty() &&
			boundary.sampled_voxel_fragments.size() == 64,
			"single-cell wedge has zero inter-cell apertures and a 4^3 preview map");
		check(interior != invalid_fragment && exterior != invalid_fragment &&
			interior != exterior && !connected(boundary, interior, exterior),
			"closed seam has no interior/exterior cross-fabric connection");
		check(!distinct_patch_sides &&
			boundary.diagnostic_unverified_same_fragment_patches == 1,
			"under-resolved wedge cusp is exposed as one preview-only diagnostic patch");
	}

	void diagnostic_same_fragment_case()
	{
		// The nested, unrelated face historically authorizes a sampled side collapse
		// only in the explicitly unsafe qualitative-preview mode. This diagnostic is
		// intentionally frozen so the future strict replacement can delete it visibly.
		const TriMesh mesh = mesh_from_quads({
			{{ { 0.5, 0, 0 }, { 0.5, 0.5, 0 }, { 0.5, 0.5, 1 }, { 0.5, 0, 1 } }},
			{{ { 0.5, 0.1, 0.2 }, { 0.5, 0.4, 0.2 }, { 0.5, 0.4, 0.8 }, { 0.5, 0.1, 0.8 } }} });
		const EmbeddedBoundary boundary = build_preview(mesh);
		report("diagnostic same-fragment", boundary);

		check(boundary.ready_for_flow(), "diagnostic preview explicitly accepts the legacy collapse");
		check(boundary.fragments.size() == 1 && boundary.patches.size() == 4,
			"diagnostic case retains one sampled fluid fragment and four patches");
		check(boundary.apertures.empty() && boundary.connections.empty() &&
			boundary.sampled_voxel_fragments.size() == 64,
			"diagnostic case keeps the concrete zero-aperture 4^3 preview layout");
		check(boundary.accepted_free_edge_same_fragment_patches == 2 &&
			boundary.diagnostic_unverified_same_fragment_patches == 2 &&
			boundary.diagnostic_unverified_same_fragment_area > 0.0,
			"unsafe same-fragment acceptance remains observable in diagnostic counters");
		check(boundary.rejected_same_fragment_patches == 0 &&
			boundary.diagnostic_unmapped_surface_patches == 0,
			"accepted diagnostic collapse is neither rejected nor silently unmapped");
	}
}

int main()
{
	opening_case();
	t_junction_case();
	closed_wedge_case();
	diagnostic_same_fragment_case();
	if (failures != 0)
		std::fprintf(stderr, "[preview-topology] %d regression check(s) failed\n", failures);
	else
		std::printf("[preview-topology] all qualitative-preview topology checks passed\n");
	return failures == 0 ? 0 : 1;
}
