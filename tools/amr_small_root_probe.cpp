// Deterministic regression for a signed-bracketed tiny connector that cannot be
// agglomerated into either neighbour without reversing the remaining aperture.
//
// The exact-cell callback is an established test seam: the trigger mesh only selects
// the two Cartesian cells, while the callback supplies conservative value-type cell
// decompositions.  The resulting fragments and shared-face common refinement then
// travel through the production EB transaction and AMR stabilization code.

#include "core/fluid/amr_advection.h"
#include "core/fluid/amr_eb.h"
#include "core/fluid/amr_grid.h"
#include "core/fluid/amr_pressure.h"
#include "core/geometry/exact_cell_decomposition.h"
#include "core/geometry/triangle_bvh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace paracfd::core;

namespace
{
	int failures = 0;

	void check(bool condition, const char* description)
	{
		std::printf("[amr-small-root] %-68s %s\n", description,
			condition ? "PASS" : "FAIL");
		if (!condition) ++failures;
	}

	bool near(double a, double b, double tolerance = 1e-12)
	{
		return std::abs(a - b) <= tolerance *
			std::max({1.0, std::abs(a), std::abs(b)});
	}

	ExactCellPoint point(double x, double y, double z)
	{
		return {x, y, z};
	}

	Vec3d vector(ExactCellPoint p)
	{
		return {p[0], p[1], p[2]};
	}

	ExactCellPlanarRegion region(ExactCellPoint a, ExactCellPoint b,
		ExactCellPoint c, ExactCellPoint d)
	{
		ExactCellPlanarRegion out;
		out.outer_loop = {a, b, c, d};
		out.convex_pieces = {{{a, b, c}}, {{a, c, d}}};
		for (const auto& triangle : out.convex_pieces)
		{
			const Vec3d av = vector(triangle[0]);
			const Vec3d bv = vector(triangle[1]);
			const Vec3d cv = vector(triangle[2]);
			const double area = 0.5 * std::sqrt(length2(cross(bv - av, cv - av)));
			const Vec3d centroid = (av + bv + cv) / 3.0;
			out.convex_piece_area_sum += area;
			out.convex_piece_first_moment[0] += area * centroid.x;
			out.convex_piece_first_moment[1] += area * centroid.y;
			out.convex_piece_first_moment[2] += area * centroid.z;
		}
		return out;
	}

	void append_oriented_quad(std::vector<ExactCellOrientedTriangle>& triangles,
		ExactCellPoint a, ExactCellPoint b, ExactCellPoint c, ExactCellPoint d)
	{
		triangles.push_back({{a, b, c}});
		triangles.push_back({{a, c, d}});
	}

	ExactCellFragment box_fragment(int id, ExactCellPoint lo, ExactCellPoint hi)
	{
		const ExactCellPoint p000 = point(lo[0], lo[1], lo[2]);
		const ExactCellPoint p001 = point(lo[0], lo[1], hi[2]);
		const ExactCellPoint p010 = point(lo[0], hi[1], lo[2]);
		const ExactCellPoint p011 = point(lo[0], hi[1], hi[2]);
		const ExactCellPoint p100 = point(hi[0], lo[1], lo[2]);
		const ExactCellPoint p101 = point(hi[0], lo[1], hi[2]);
		const ExactCellPoint p110 = point(hi[0], hi[1], lo[2]);
		const ExactCellPoint p111 = point(hi[0], hi[1], hi[2]);

		ExactCellFragment out;
		out.id = id;
		out.volume = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
		out.centroid = out.interior_witness = point(
			0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]),
			0.5 * (lo[2] + hi[2]));
		append_oriented_quad(out.boundary_triangles, p000, p001, p011, p010); // -X
		append_oriented_quad(out.boundary_triangles, p100, p110, p111, p101); // +X
		append_oriented_quad(out.boundary_triangles, p000, p100, p101, p001); // -Y
		append_oriented_quad(out.boundary_triangles, p010, p011, p111, p110); // +Y
		append_oriented_quad(out.boundary_triangles, p000, p010, p110, p100); // -Z
		append_oriented_quad(out.boundary_triangles, p001, p101, p111, p011); // +Z
		return out;
	}

	void append_boundary_apertures(ExactCellDecomposition& out, int fragment,
		ExactCellPoint fragment_lo, ExactCellPoint fragment_hi,
		ExactCellPoint cell_lo, ExactCellPoint cell_hi)
	{
		const double tolerance = 1e-12;
		auto append = [&](int axis, bool upper, ExactCellPoint a,
			ExactCellPoint b, ExactCellPoint c, ExactCellPoint d)
		{
			ExactCellBoxAperture aperture;
			aperture.axis = static_cast<std::int8_t>(axis);
			aperture.upper = upper;
			aperture.region = region(a, b, c, d);
			aperture.area = aperture.region.convex_piece_area_sum;
			for (int q = 0; q < 3; ++q)
				aperture.centroid[q] = aperture.region.convex_piece_first_moment[q]
					/ aperture.area;
			aperture.fragment = fragment;
			out.box_apertures.push_back(std::move(aperture));
		};

		const double x0 = fragment_lo[0], y0 = fragment_lo[1], z0 = fragment_lo[2];
		const double x1 = fragment_hi[0], y1 = fragment_hi[1], z1 = fragment_hi[2];
		if (std::abs(x0 - cell_lo[0]) <= tolerance)
			append(0, false, point(x0,y0,z0), point(x0,y0,z1),
				point(x0,y1,z1), point(x0,y1,z0));
		if (std::abs(x1 - cell_hi[0]) <= tolerance)
			append(0, true, point(x1,y0,z0), point(x1,y1,z0),
				point(x1,y1,z1), point(x1,y0,z1));
		if (std::abs(y0 - cell_lo[1]) <= tolerance)
			append(1, false, point(x0,y0,z0), point(x1,y0,z0),
				point(x1,y0,z1), point(x0,y0,z1));
		if (std::abs(y1 - cell_hi[1]) <= tolerance)
			append(1, true, point(x0,y1,z0), point(x0,y1,z1),
				point(x1,y1,z1), point(x1,y1,z0));
		if (std::abs(z0 - cell_lo[2]) <= tolerance)
			append(2, false, point(x0,y0,z0), point(x0,y1,z0),
				point(x1,y1,z0), point(x1,y0,z0));
		if (std::abs(z1 - cell_hi[2]) <= tolerance)
			append(2, true, point(x0,y0,z1), point(x1,y0,z1),
				point(x1,y1,z1), point(x0,y1,z1));
	}

	void append_patch(ExactCellDecomposition& out, std::uint64_t source_triangle,
		std::uint64_t source_face, ExactCellPoint a, ExactCellPoint b,
		ExactCellPoint c, ExactCellPoint d, ExactCellPoint normal,
		int plus_fragment, int minus_fragment)
	{
		ExactCellSurfacePatch patch;
		patch.source_triangle_id = source_triangle;
		patch.source_face_id = source_face;
		patch.region = region(a, b, c, d);
		patch.area = patch.region.convex_piece_area_sum;
		for (int q = 0; q < 3; ++q)
			patch.centroid[q] = patch.region.convex_piece_first_moment[q] / patch.area;
		patch.normal = normal;
		patch.plus_fragment = plus_fragment;
		patch.minus_fragment = minus_fragment;
		out.surface_patches.push_back(std::move(patch));
	}

	ExactCellDecomposition connector_decomposition(const ExactCellInput& input,
		double tiny_width)
	{
		const ExactCellPoint cell_lo = input.cell_min;
		const ExactCellPoint cell_hi = input.cell_max;
		ExactCellDecomposition out;
		out.expected_cell_volume = 1.0;
		out.fragment_volume_sum = 1.0;
		out.input_clipped_surface_area_sum = 1.0;
		out.history_surface_area_sum = 1.0;

		if (cell_lo[0] < 0.5)
		{
			// A and B partition the lower cell across Z. Both independently expose
			// half of the X=1 face and are opposite pressure sides of real fabric.
			const ExactCellPoint a_lo = cell_lo;
			const ExactCellPoint a_hi = point(cell_hi[0], cell_hi[1], 0.5);
			const ExactCellPoint b_lo = point(cell_lo[0], cell_lo[1], 0.5);
			const ExactCellPoint b_hi = cell_hi;
			out.fragments.push_back(box_fragment(0, a_lo, a_hi));
			out.fragments.push_back(box_fragment(1, b_lo, b_hi));
			append_boundary_apertures(out, 0, a_lo, a_hi, cell_lo, cell_hi);
			append_boundary_apertures(out, 1, b_lo, b_hi, cell_lo, cell_hi);
			append_patch(out, 0, 0, point(cell_lo[0],cell_lo[1],0.5),
				point(cell_hi[0],cell_lo[1],0.5), point(cell_hi[0],cell_hi[1],0.5),
				point(cell_lo[0],cell_hi[1],0.5), point(0,0,1), 1, 0);
		}
		else
		{
			// T is a 0.001 h^3 slab. Its full X=1 aperture is split by the
			// neighbouring A/B partition, giving the two-edge aperture fan. U is
			// separated from T by fabric, so it is never a legal donor.
			const ExactCellPoint t_lo = cell_lo;
			const ExactCellPoint t_hi = point(cell_lo[0] + tiny_width,
				cell_hi[1], cell_hi[2]);
			const ExactCellPoint u_lo = point(cell_lo[0] + tiny_width,
				cell_lo[1], cell_lo[2]);
			const ExactCellPoint u_hi = cell_hi;
			out.fragments.push_back(box_fragment(0, t_lo, t_hi));
			out.fragments.push_back(box_fragment(1, u_lo, u_hi));
			append_boundary_apertures(out, 0, t_lo, t_hi, cell_lo, cell_hi);
			append_boundary_apertures(out, 1, u_lo, u_hi, cell_lo, cell_hi);
			const double x = cell_lo[0] + tiny_width;
			append_patch(out, 2, 1, point(x,cell_lo[1],cell_lo[2]),
				point(x,cell_hi[1],cell_lo[2]), point(x,cell_hi[1],cell_hi[2]),
				point(x,cell_lo[1],cell_hi[2]), point(1,0,0), 1, 0);
		}
		return out;
	}

	ExactCellDecomposition tiny_connector_decomposer(const ExactCellInput& input)
	{
		return connector_decomposition(input, 1e-3);
	}

	ExactCellDecomposition unbracketed_connector_decomposer(const ExactCellInput& input)
	{
		// The raw upper centroid is only 0.5e-12 h from the face, inside the
		// production 1e-10 h signed-bracketing exclusion. Retention must remain
		// fail-closed even though all records are otherwise finite and positive.
		return connector_decomposition(input, 1e-12);
	}

	TriMesh trigger_mesh()
	{
		// The synthetic exact decomposer owns the manufactured topology. These two
		// interior quads merely make both cells part of the sparse EB atlas.
		TriMesh mesh;
		mesh.bbox_min = {1e30f, 1e30f, 1e30f};
		mesh.bbox_max = {-1e30f, -1e30f, -1e30f};
		const std::array<std::array<Vec3d,4>,2> quads{{
			{{{0.2,0.2,0.2},{0.4,0.2,0.2},{0.4,0.4,0.2},{0.2,0.4,0.2}}},
			{{{1.2,0.2,0.2},{1.4,0.2,0.2},{1.4,0.4,0.2},{1.2,0.4,0.2}}}
		}};
		for (std::size_t face = 0; face < quads.size(); ++face)
		{
			const std::uint32_t base = static_cast<std::uint32_t>(mesh.vertex_count());
			for (Vec3d p : quads[face])
			{
				mesh.positions.insert(mesh.positions.end(), {static_cast<float>(p.x),
					static_cast<float>(p.y), static_cast<float>(p.z)});
				mesh.positions_fp64.insert(mesh.positions_fp64.end(), {p.x,p.y,p.z});
				for (int axis = 0; axis < 3; ++axis)
				{
					mesh.bbox_min[axis] = std::min(mesh.bbox_min[axis],
						static_cast<float>(p[axis]));
					mesh.bbox_max[axis] = std::max(mesh.bbox_max[axis],
						static_cast<float>(p[axis]));
				}
			}
			mesh.indices.insert(mesh.indices.end(), {base,base+1,base+2,base,base+2,base+3});
			mesh.source_face_ids.push_back(static_cast<std::uint32_t>(face));
			mesh.source_face_ids.push_back(static_cast<std::uint32_t>(face));
		}
		return mesh;
	}

	AmrEmbeddedBoundaryAtlas build_fixture(bool retain,
		ExactCellDecomposer decomposer = &tiny_connector_decomposer)
	{
		const TriMesh mesh = trigger_mesh();
		const TriangleBvh bvh(mesh);
		const AmrHierarchy hierarchy = AmrHierarchy::uniform(
			{{0,0,0},{2,1,1}}, 1.0, 1, 1);
		EmbeddedBoundaryBuildOptions options;
		options.exact_cell_decomposer = decomposer;
		options.min_volume_fraction = 0.05;
		options.retain_signed_bracketed_small_roots_for_face_state = retain;
		return build_amr_embedded_boundary_atlas(hierarchy, mesh, bvh, options);
	}

	void regression()
	{
		const AmrEmbeddedBoundaryAtlas strict = build_fixture(false);
		check(!strict.ready_for_flow() && strict.owned_unresolved_count() == 1,
			"without face-state retention the tiny fan fails closed");
		const AmrEmbeddedBoundaryAtlas raw_unbracketed = build_fixture(true,
			&unbracketed_connector_decomposer);
		check(!raw_unbracketed.ready_for_flow() &&
			raw_unbracketed.owned_unresolved_count() == 1,
			"retention never accepts a raw aperture inside bracketing tolerance");

		const AmrEmbeddedBoundaryAtlas retained = build_fixture(true);
		check(retained.ready_for_flow() && retained.levels.size() == 1,
			"signed-bracketed face-state policy makes the atlas flow-ready");
		if (retained.levels.empty()) return;
		const AmrEbLevelAtlas& level = retained.levels.front();
		const EmbeddedBoundary& eb = level.topology;
		check(level.face_state_retained_small_roots == 1 &&
			near(level.face_state_retained_small_volume, 1e-3) &&
			near(level.minimum_face_state_retained_volume_fraction, 1e-3),
			"exactly one 0.001 h^3 pressure root is retained and reported");

		const AmrHierarchy hierarchy = AmrHierarchy::uniform(
			{{0,0,0},{2,1,1}}, 1.0, 1, 1);
		CompositeAmrPressureBuildOptions pressure_options;
		pressure_options.unsupported_nonorthogonal_correction =
			UnsupportedNonorthogonalCorrectionPolicy::qualitative_preview_first_order_orthogonal;
		const CompositeAmrPressureSystem pressure_system =
			build_composite_amr_pressure_system(hierarchy, retained, pressure_options);
		check(pressure_system.face_state_retained_small_root_count == 1,
			"pressure system propagates the mapped retained-root count");

		AmrEmbeddedBoundaryAtlas unrepresentable_atlas = retained;
		for (FluidFragment& fragment : unrepresentable_atlas.levels.front().topology.fragments)
			if (fragment.face_state_retained)
				fragment.volume = 0.5 * static_cast<double>(std::numeric_limits<Real>::min());
		bool unrepresentable_atlas_rejected = false;
		try
		{
			(void)build_composite_amr_pressure_system(
				hierarchy, unrepresentable_atlas, pressure_options);
		}
		catch (const std::runtime_error& error)
		{
			unrepresentable_atlas_rejected =
				std::string(error.what()).find("positive normal production Real") !=
				std::string::npos;
		}
		check(unrepresentable_atlas_rejected,
			"pressure construction rejects a retained root below normal Real range");

		AmrEmbeddedBoundaryAtlas stale_metadata = retained;
		stale_metadata.levels.front().face_state_retained_small_roots = 0;
		bool stale_metadata_rejected = false;
		try
		{
			(void)build_composite_amr_pressure_system(
				hierarchy, stale_metadata, pressure_options);
		}
		catch (const std::runtime_error& error)
		{
			stale_metadata_rejected =
				std::string(error.what()).find("does not match mapped pressure roots") !=
				std::string::npos;
		}
		check(stale_metadata_rejected,
			"pressure construction rejects stale retained-root metadata");

		bool cpu_collocated_rejected = false;
		try
		{
			CompositeCellMomentumState state;
			diffuse_composite_cell_momentum_cpu(pressure_system, 1e-5, 1e-3, state);
		}
		catch (const std::invalid_argument& error)
		{
			cpu_collocated_rejected =
				std::string(error.what()).find("face-state-retained small pressure roots") !=
				std::string::npos;
		}
		check(cpu_collocated_rejected,
			"CPU collocated operators reject retained small roots");

		DeviceAmrFields device_fields(hierarchy);
		CompositeAmrPressureSystem unrepresentable_system = pressure_system;
		int smallest_active = -1;
		for (int q = 0; q < unrepresentable_system.storage_size; ++q)
			if (unrepresentable_system.active[q] &&
				(smallest_active < 0 || unrepresentable_system.volume[q] <
					unrepresentable_system.volume[smallest_active]))
				smallest_active = q;
		if (smallest_active >= 0)
			unrepresentable_system.volume[smallest_active] =
				0.5 * static_cast<double>(std::numeric_limits<Real>::min());
		bool projection_pack_rejected = false;
		try
		{
			DeviceCompositeAmrProjection projection(unrepresentable_system, device_fields);
		}
		catch (const std::invalid_argument& error)
		{
			projection_pack_rejected =
				std::string(error.what()).find("positive normal production Real") !=
				std::string::npos;
		}
		check(smallest_active >= 0 && projection_pack_rejected,
			"GPU projection rejects an active volume below normal Real range");

		bool collocated_rejected = false;
		try
		{
			DeviceCompositeCellMomentumTransport transport(pressure_system, device_fields);
		}
		catch (const std::invalid_argument& error)
		{
			collocated_rejected =
				std::string(error.what()).find("face-state-retained small pressure roots") !=
				std::string::npos;
		}
		check(collocated_rejected,
			"experimental collocated transport rejects retained small roots");

		int retained_fragment = -1;
		int retained_count = 0;
		double fragment_volume_sum = 0.0;
		Vec3d fragment_first_moment{};
		for (int q = 0; q < static_cast<int>(eb.fragments.size()); ++q)
		{
			const FluidFragment& fragment = eb.fragments[q];
			fragment_volume_sum += fragment.volume;
			fragment_first_moment = fragment_first_moment + fragment.centroid * fragment.volume;
			if (fragment.face_state_retained)
			{
				retained_fragment = q;
				++retained_count;
			}
		}
		check(retained_count == 1 && retained_fragment >= 0 &&
			eb.fragments[retained_fragment].merge_target == irregular_fragment(retained_fragment),
			"tiny connector remains an independent pressure root");
		check(near(fragment_volume_sum, 2.0) &&
			near(fragment_first_moment.x, 2.0) &&
			near(fragment_first_moment.y, 1.0) &&
			near(fragment_first_moment.z, 1.0),
			"fixture volume and first moment remain conservative");

		int incident = 0;
		bool bracketed = true;
		double fan_area = 0.0;
		std::vector<FragmentRef> neighbours;
		const FragmentRef retained_ref = irregular_fragment(retained_fragment);
		for (const FaceAperture& aperture : eb.apertures)
		{
			if (aperture.fragment_a != retained_ref && aperture.fragment_b != retained_ref)
				continue;
			++incident;
			fan_area += aperture.area;
			const FragmentRef other = aperture.fragment_a == retained_ref
				? aperture.fragment_b : aperture.fragment_a;
			neighbours.push_back(other);
			const Vec3d lower = eb.fragment_centroid(aperture.fragment_a);
			const Vec3d upper = eb.fragment_centroid(aperture.fragment_b);
			bracketed = bracketed &&
				aperture.centroid[aperture.axis] > lower[aperture.axis] &&
				aperture.centroid[aperture.axis] < upper[aperture.axis];
		}
		std::sort(neighbours.begin(), neighbours.end());
		neighbours.erase(std::unique(neighbours.begin(), neighbours.end()), neighbours.end());
		std::printf("[amr-small-root] retained fan records=%d neighbours=%zu total-apertures=%zu\n",
			incident, neighbours.size(), eb.apertures.size());
		bool donor_short_circuit = false;
		if (neighbours.size() == 2)
			for (const FaceAperture& aperture : eb.apertures)
				donor_short_circuit = donor_short_circuit ||
					((aperture.fragment_a == neighbours[0] && aperture.fragment_b == neighbours[1]) ||
					 (aperture.fragment_a == neighbours[1] && aperture.fragment_b == neighbours[0]));
		check(incident >= 2 && neighbours.size() == 2 && bracketed &&
			near(fan_area, 1.0) && !donor_short_circuit,
			"both positive-area fan apertures survive with signed bracketing");
		bool every_pairwise_merge_reverses_survivor = neighbours.size() == 2;
		if (neighbours.size() == 2)
			for (FragmentRef donor : neighbours)
			{
				const double retained_volume = eb.fragment_volume(retained_ref);
				const double donor_volume = eb.fragment_volume(donor);
				const double merged_x = (retained_volume * eb.fragment_centroid(retained_ref).x +
					donor_volume * eb.fragment_centroid(donor).x) /
					(retained_volume + donor_volume);
				every_pairwise_merge_reverses_survivor =
					every_pairwise_merge_reverses_survivor && merged_x < 1.0;
			}
		check(every_pairwise_merge_reverses_survivor,
			"merging into either donor reverses the other aperture distance");

		bool independent_patch_sides = true;
		double patch_area = 0.0;
		std::vector<std::uint32_t> patch_faces;
		for (const SurfacePatch& patch : eb.patches)
		{
			independent_patch_sides = independent_patch_sides &&
				patch.plus_fragment != patch.minus_fragment;
			patch_area += patch.area;
			patch_faces.push_back(patch.source_face_id);
		}
		std::sort(patch_faces.begin(), patch_faces.end());
		patch_faces.erase(std::unique(patch_faces.begin(), patch_faces.end()), patch_faces.end());
		std::printf("[amr-small-root] patch records=%zu distinct-sides=%d\n",
			eb.patches.size(), independent_patch_sides ? 1 : 0);
		check(patch_faces.size() == 2 && near(patch_area, 2.0) && independent_patch_sides,
			"both fabric sheets retain independent plus/minus pressure sides");
	}
}

int main()
{
	regression();
	std::printf("[amr-small-root] failures=%d\n", failures);
	return failures == 0 ? 0 : 1;
}
