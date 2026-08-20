#include "core/geometry/uv_constraint_mesh.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using namespace paracfd::core;

namespace
{
	void require(bool condition, const std::string& message)
	{
		if (!condition) throw std::runtime_error(message);
	}

	void require_close(double actual, double expected, double tolerance, const char* message)
	{
		if (std::abs(actual - expected) > tolerance)
			throw std::runtime_error(std::string(message) + ": actual=" + std::to_string(actual)
				+ " expected=" + std::to_string(expected));
	}

	UvConstraintMesh square_mesh()
	{
		return UvConstraintMesh({
			{{0.0, 0.0}}, {{1.0, 0.0}}, {{1.0, 1.0}}, {{0.0, 1.0}}},
			{{{0, 1, 2}}, {{0, 2, 3}}});
	}

	UvConstraintMesh one_to_two_mesh()
	{
		// The bottom trim boundary is already represented as two edges, while
		// square_mesh() represents the same geometric boundary with one edge.
		return UvConstraintMesh({
			{{0.0, 0.0}}, {{0.5, 0.0}}, {{1.0, 0.0}}, {{1.0, 1.0}}, {{0.0, 1.0}}},
			{{{4, 0, 1}}, {{4, 1, 2}}, {{4, 2, 3}}});
	}

	UvConstraintMesh annulus_mesh()
	{
		return UvConstraintMesh({
			{{0.0, 0.0}}, {{1.0, 0.0}}, {{1.0, 1.0}}, {{0.0, 1.0}},
			{{0.35, 0.35}}, {{0.65, 0.35}}, {{0.65, 0.65}}, {{0.35, 0.65}}},
			{
				{{0, 1, 5}}, {{0, 5, 4}},
				{{1, 2, 6}}, {{1, 6, 5}},
				{{2, 3, 7}}, {{2, 7, 6}},
				{{3, 0, 4}}, {{3, 4, 7}}
			});
	}

	UvConstraintMesh strip_mesh(const std::vector<double>& x_coordinates,
		double coordinate_tolerance)
	{
		std::vector<UvVertex> vertices;
		std::vector<UvTriangle> triangles;
		vertices.reserve(2 * x_coordinates.size());
		for (double x : x_coordinates)
		{
			vertices.push_back({{x, 0.0}});
			vertices.push_back({{x, 1.0}});
		}
		for (std::size_t column = 0; column + 1 < x_coordinates.size(); ++column)
		{
			const std::uint32_t bottom_left = static_cast<std::uint32_t>(2 * column);
			const std::uint32_t top_left = bottom_left + 1;
			const std::uint32_t bottom_right = bottom_left + 2;
			const std::uint32_t top_right = bottom_left + 3;
			triangles.push_back({{bottom_left, bottom_right, top_right}});
			triangles.push_back({{bottom_left, top_right, top_left}});
		}
		return UvConstraintMesh(std::move(vertices), std::move(triangles),
			{coordinate_tolerance});
	}

	std::vector<std::uint32_t> constrained_neighbours(const UvConstraintMesh& mesh,
		std::uint32_t vertex)
	{
		std::vector<std::uint32_t> result;
		for (const UvEdgeView& edge : mesh.edges())
			if (!edge.constraint_ids.empty()
				&& (edge.vertices[0] == vertex || edge.vertices[1] == vertex))
				result.push_back(edge.vertices[0] == vertex ? edge.vertices[1] : edge.vertices[0]);
		std::sort(result.begin(), result.end());
		return result;
	}

	void require_chain(const UvConstraintMesh& mesh, std::uint32_t topology_a,
		std::uint32_t topology_b, std::uint32_t topology_c, std::uint64_t constraint)
	{
		const auto a = mesh.vertex_for_topology_id(topology_a);
		const auto b = mesh.vertex_for_topology_id(topology_b);
		const auto c = mesh.vertex_for_topology_id(topology_c);
		require(a && b && c, "canonical topology chain is missing a sample");
		require(mesh.edge_has_constraint(*a, *b, constraint)
			&& mesh.edge_has_constraint(*b, *c, constraint),
			"canonical samples are not joined by the requested constraint edges");
	}

	bool same_mesh(const UvConstraintMesh& a, const UvConstraintMesh& b)
	{
		if (a.vertices().size() != b.vertices().size()
			|| a.triangles().size() != b.triangles().size()) return false;
		for (std::size_t i = 0; i < a.vertices().size(); ++i)
		{
			const UvVertex& lhs = a.vertices()[i];
			const UvVertex& rhs = b.vertices()[i];
			if (std::bit_cast<std::uint64_t>(lhs.uv.u) != std::bit_cast<std::uint64_t>(rhs.uv.u)
				|| std::bit_cast<std::uint64_t>(lhs.uv.v) != std::bit_cast<std::uint64_t>(rhs.uv.v)
				|| lhs.topology_id != rhs.topology_id) return false;
		}
		for (std::size_t i = 0; i < a.triangles().size(); ++i)
			if (a.triangles()[i].vertices != b.triangles()[i].vertices) return false;
		const auto a_edges = a.edges();
		const auto b_edges = b.edges();
		if (a_edges.size() != b_edges.size()) return false;
		for (std::size_t i = 0; i < a_edges.size(); ++i)
			if (a_edges[i].vertices != b_edges[i].vertices
				|| a_edges[i].boundary != b_edges[i].boundary
				|| a_edges[i].constraint_ids != b_edges[i].constraint_ids) return false;
		return true;
	}

	void t_junction_test()
	{
		UvConstraintMesh mesh = square_mesh();
		std::string error;
		const std::vector<UvConstraintSample> horizontal{
			{{0.1, 0.5}, 100}, {{0.9, 0.5}, 101}};
		require(mesh.insert_polyline(horizontal, 10, &error), "horizontal insertion: " + error);
		const std::vector<UvConstraintSample> vertical{
			{{0.5, 0.9}, 103}, {{0.5, 0.5}, 102}};
		require(mesh.insert_polyline(vertical, 20, &error), "T stem insertion: " + error);
		require(mesh.validate(&error), "T-junction validation: " + error);
		require_close(mesh.area(), 1.0, 1.0e-14, "T-junction changed face area");
		require(mesh.boundary_loop_count() == 1, "T-junction changed outer boundary topology");

		const auto left = mesh.vertex_for_topology_id(100);
		const auto right = mesh.vertex_for_topology_id(101);
		const auto junction = mesh.vertex_for_topology_id(102);
		const auto stem = mesh.vertex_for_topology_id(103);
		require(left && right && junction && stem, "T-junction lost canonical topology IDs");
		require(mesh.edge_has_constraint(*left, *junction, 10)
			&& mesh.edge_has_constraint(*junction, *right, 10),
			"splitting the horizontal constraint did not preserve its tag");
		require(mesh.edge_has_constraint(*stem, *junction, 20),
			"T stem is not a constrained edge");
		require(constrained_neighbours(mesh, *junction).size() == 3,
			"interior T vertex does not have constrained degree three");
	}

	void one_to_two_sampling_test()
	{
		std::string error;
		UvConstraintMesh missing_union = one_to_two_mesh();
		const UvConstraintMesh before_failure = missing_union;
		const std::vector<UvConstraintSample> insufficient{
			{{0.0, 0.0}, 200}, {{1.0, 0.0}, 202}};
		require(!missing_union.insert_polyline(insufficient, 30, &error),
			"one-to-two mismatch was silently accepted without the fine-side sample");
		require(error.find("include it in the canonical sample chain") != std::string::npos,
			"one-to-two mismatch did not produce the sample-union diagnostic");
		require(same_mesh(missing_union, before_failure),
			"failed one-to-two insertion was not transactional");

		const std::vector<UvConstraintSample> shared_chain{
			{{0.0, 0.0}, 200}, {{0.5, 0.0}, 201}, {{1.0, 0.0}, 202}};
		UvConstraintMesh coarse = square_mesh();
		UvConstraintMesh fine = one_to_two_mesh();
		require(coarse.insert_polyline(shared_chain, 30, &error), "coarse chain insertion: " + error);
		require(fine.insert_polyline(shared_chain, 30, &error), "fine chain insertion: " + error);
		require_chain(coarse, 200, 201, 202, 30);
		require_chain(fine, 200, 201, 202, 30);
		require(coarse.vertices().size() == 5 && fine.vertices().size() == 5,
			"one-to-two normalization created a triangulation-dependent chain vertex");
		for (const UvConstraintMesh* mesh : {&coarse, &fine})
		{
			const auto middle = mesh->vertex_for_topology_id(201);
			require(middle.has_value(), "one-to-two midpoint lacks its canonical topology ID");
			unsigned boundary_constraint_edges = 0;
			for (const UvEdgeView& edge : mesh->edges())
				if (edge.boundary && std::find(edge.constraint_ids.begin(), edge.constraint_ids.end(), 30)
					!= edge.constraint_ids.end()) ++boundary_constraint_edges;
			require(boundary_constraint_edges == 2,
				"one-to-two insertion did not retain two trim-boundary pieces");
		}

		std::uint32_t allocator = 1000;
		require(coarse.assign_unclaimed_topology_ids(allocator, &error),
			"coarse local topology assignment: " + error);
		require(fine.assign_unclaimed_topology_ids(allocator, &error),
			"fine local topology assignment: " + error);
		for (const UvConstraintMesh* mesh : {&coarse, &fine})
			for (const UvVertex& vertex : mesh->vertices())
				require(vertex.topology_id != UvVertex::no_topology_id,
					"topology export retained an unassigned vertex");
	}

	void hole_and_boundary_test()
	{
		UvConstraintMesh mesh = annulus_mesh();
		std::string error;
		require(mesh.validate(&error), "initial annulus validation: " + error);
		require(mesh.boundary_loop_count() == 2, "annulus did not expose outer and inner loops");
		require(!mesh.contains({0.5, 0.5}), "hole was incorrectly classified as triangulated face");
		const double original_area = mesh.area();
		require_close(original_area, 0.91, 1.0e-14, "annulus area is incorrect");

		const std::vector<UvConstraintSample> bridge{
			{{0.5, 0.0}, 300}, {{0.5, 0.35}, 301}};
		require(mesh.insert_polyline(bridge, 40, &error), "annulus boundary insertion: " + error);
		require(mesh.validate(&error), "annulus post-insertion validation: " + error);
		require(mesh.boundary_loop_count() == 2,
			"splitting the outer/hole boundaries changed the number of trim loops");
		require_close(mesh.area(), original_area, 1.0e-14,
			"constraint insertion filled or clipped part of the hole");
		require(!mesh.contains({0.5, 0.5}), "constraint insertion filled the hole");
		const auto outer = mesh.vertex_for_topology_id(300);
		const auto inner = mesh.vertex_for_topology_id(301);
		require(outer && inner && mesh.edge_has_constraint(*outer, *inner, 40),
			"outer-to-hole boundary constraint was not inserted as one exact edge");

		const UvConstraintMesh before_rejection = mesh;
		const std::vector<UvConstraintSample> into_hole{
			{{0.15, 0.5}, 302}, {{0.5, 0.5}, 303}};
		require(!mesh.insert_polyline(into_hole, 41, &error),
			"constraint endpoint inside the hole was accepted");
		require(error.find("inside one of its holes") != std::string::npos,
			"hole rejection did not identify the trimmed-domain failure");
		require(same_mesh(mesh, before_rejection), "hole rejection was not transactional");
	}

	UvConstraintMesh deterministic_fixture()
	{
		UvConstraintMesh mesh = square_mesh();
		std::string error;
		const std::vector<UvConstraintSample> polyline{
			{{0.08, 0.28}, 400}, {{0.35, 0.62}, 401}, {{0.88, 0.73}, 402}};
		require(mesh.insert_polyline(polyline, 50, &error),
			"determinism fixture insertion: " + error);
		return mesh;
	}

	void deterministic_replay_test()
	{
		const UvConstraintMesh first = deterministic_fixture();
		const UvConstraintMesh second = deterministic_fixture();
		require(same_mesh(first, second),
			"identical UV input and constraints did not produce bitwise-identical topology");
	}

	void exact_existing_vertex_provenance_test()
	{
		const auto make_mesh = []
		{
			return UvConstraintMesh({
				{{0.0, 0.0}}, {{1.0, 0.0}}, {{1.0, 1.0}}, {{0.0, 1.0}}},
				{{{0, 1, 2}}, {{0, 2, 3}}}, {1.1e-9});
		};
		// This coordinate is farther than the UV vertex tolerance from vertex zero,
		// yet lies within that tolerance of both incident bottom/diagonal edges.  A
		// coordinate-only insertion is therefore intentionally ambiguous.
		const UvPoint p{2.0e-9, 1.0e-9};
		std::string error;
		UvConstraintMesh ambiguous = make_mesh();
		const std::vector<UvConstraintSample> coordinate_only{
			{p, 800}, {{1.0, 0.0}, 801}};
		require(!ambiguous.insert_polyline(coordinate_only, 55, &error)
			&& error.find("multiple distinct mesh edges") != std::string::npos,
			"manufactured near-vertex sample did not exercise ambiguous UV classification");

		UvConstraintMesh claimed = make_mesh();
		const std::vector<UvConstraintSample> exact_claim{
			{p, 800, 0}, {{1.0, 0.0}, 801, 1}};
		require(claimed.insert_polyline(exact_claim, 55, &error),
			"exact existing-vertex provenance was not honored: " + error);
		const auto first = claimed.vertex_for_topology_id(800);
		const auto second = claimed.vertex_for_topology_id(801);
		require(first && second && *first == 0 && *second == 1
			&& claimed.edge_has_constraint(0, 1, 55),
			"exact provenance did not claim/tag the named original vertices");

		UvConstraintMesh boundary = make_mesh();
		const std::vector<UvConstraintSample> exact_boundary{
			{p, 802, UvConstraintSample::no_existing_vertex, true},
			{{1.0, 0.0}, 803, 1, true}};
		require(boundary.insert_polyline(exact_boundary, 56, &error),
			"exact trim-boundary provenance did not disambiguate boundary/interior edges: "
			+ error);
		const auto boundary_split = boundary.vertex_for_topology_id(802);
		const auto boundary_end = boundary.vertex_for_topology_id(803);
		require(boundary_split && boundary_end
			&& boundary.edge_has_constraint(*boundary_split, *boundary_end, 56),
			"trim-boundary provenance did not split/tag the unique boundary edge");
	}

	void parameter_tolerance_test()
	{
		std::string error;
		// These crossings are close as dimensionless parameters but separated by 0.025
		// native UV units, well beyond the 0.001 coordinate tolerance. Comparing the raw
		// parameters to a coordinate tolerance used to reject this valid long segment.
		UvConstraintMesh long_strip = strip_mesh({0.0, 400.0, 400.1, 1000.0}, 1.0e-3);
		require(long_strip.validate(&error), "long-strip validation: " + error);
		const std::vector<UvConstraintSample> resolved_crossings{
			{{0.0, 0.25}, 450}, {{1000.0, 0.25}, 451}};
		require(long_strip.insert_polyline(resolved_crossings, 55, &error),
			"dimensionless crossing separation: " + error);
		const auto long_start = long_strip.vertex_for_topology_id(450);
		const auto long_finish = long_strip.vertex_for_topology_id(451);
		require(long_start && long_finish
			&& long_strip.edge_has_constraint(*long_start, *long_finish, 55),
			"resolved long-segment crossings did not produce the requested edge");

		// Conversely, these intersections are less than one coordinate tolerance apart
		// along a short segment. They must be treated as one ambiguous parameter event.
		UvConstraintMesh short_strip = strip_mesh(
			{0.0, 0.0004, 0.0004005, 0.001}, 1.0e-6);
		require(short_strip.validate(&error), "short-strip validation: " + error);
		const UvConstraintMesh before_rejection = short_strip;
		const std::vector<UvConstraintSample> ambiguous_crossings{
			{{0.0, 0.25}, 452}, {{0.001, 0.25}, 453}};
		require(!short_strip.insert_polyline(ambiguous_crossings, 56, &error),
			"coordinate-indistinguishable crossing parameters were accepted");
		require(error.find("same parameter") != std::string::npos,
			"ambiguous crossing rejection did not identify its parameter collision");
		require(same_mesh(short_strip, before_rejection),
			"ambiguous crossing rejection was not transactional");

		// Orientation is area-valued, so its threshold must scale with the
		// actual edge length rather than an arbitrary one-UV-unit floor.
		UvConstraintMesh tiny_chart({{{0.0,0.0}},{{1.0e-3,0.0}},{{0.0,1.0e-3}}},
			{{{0,1,2}}},UvConstraintOptions{1.0e-6});
		require(tiny_chart.validate(&error),"small-chart validation: "+error);
		const std::vector<UvConstraintSample> tiny_segment{
			{{2.0e-4,2.0e-4},454},{{4.0e-4,2.0e-4},455}};
		require(tiny_chart.insert_polyline(tiny_segment,57,&error),
			"scale-neutral small-chart point location: "+error);
		require(tiny_chart.vertex_for_topology_id(454).has_value()
			&&tiny_chart.vertex_for_topology_id(455).has_value(),
			"small-chart constraint samples were not retained");
	}

	void geometric_validation_test()
	{
		std::string error;
		UvConstraintMesh crossing_edges({
			{{0.0, 0.0}}, {{2.0, 0.0}}, {{1.0, 2.0}},
			{{0.0, 1.0}}, {{1.0, -1.0}}, {{2.0, 1.0}}},
			{{{0, 1, 2}}, {{3, 4, 5}}});
		require(!crossing_edges.validate(&error),
			"geometrically crossing triangle edges were accepted");
		require(error.find("geometric mesh edges") != std::string::npos,
			"crossing-edge rejection did not identify the geometric edge invariant");
		const UvConstraintMesh before_crossing_rejection = crossing_edges;
		const std::vector<UvConstraintSample> segment{
			{{0.8, 0.2}, 520}, {{1.2, 0.2}, 521}};
		require(!crossing_edges.insert_polyline(segment, 62, &error),
			"constraint insertion accepted a mesh with crossing edges");
		require(same_mesh(crossing_edges, before_crossing_rejection),
			"crossing-edge input rejection was not transactional");

		// The inner triangle has no edge intersection with the outer triangle, so this
		// specifically exercises the positive-area overlap gate rather than edge crossing.
		UvConstraintMesh contained_triangle({
			{{0.0, 0.0}}, {{3.0, 0.0}}, {{0.0, 3.0}},
			{{0.5, 0.5}}, {{1.0, 0.5}}, {{0.5, 1.0}}},
			{{{0, 1, 2}}, {{3, 4, 5}}});
		require(!contained_triangle.validate(&error),
			"positive-area contained triangle overlap was accepted");
		require(error.find("positive-area overlap") != std::string::npos,
			"contained-triangle rejection did not identify positive-area overlap");
		const UvConstraintMesh before_overlap_rejection = contained_triangle;
		require(!contained_triangle.insert_polyline(segment, 63, &error),
			"constraint insertion accepted positive-area overlapping triangles");
		require(same_mesh(contained_triangle, before_overlap_rejection),
			"overlapping-triangle input rejection was not transactional");
	}

	void malformed_input_test()
	{
		std::string error;
		UvConstraintMesh clockwise({
			{{0.0, 0.0}}, {{1.0, 0.0}}, {{1.0, 1.0}}, {{0.0, 1.0}}},
			{{{0, 2, 1}}, {{0, 3, 2}}});
		const UvConstraintMesh before_rejection = clockwise;
		const std::vector<UvConstraintSample> segment{
			{{0.2, 0.25}, 500}, {{0.8, 0.25}, 501}};
		require(!clockwise.insert_polyline(segment, 60, &error),
			"clockwise input triangulation was accepted");
		require(error.find("not strictly CCW") != std::string::npos,
			"clockwise rejection did not identify the winding invariant");
		require(same_mesh(clockwise, before_rejection),
			"malformed-input rejection changed the source mesh");

		UvConstraintMesh non_manifold({
			{{0.0, 0.0}}, {{1.0, 0.0}}, {{0.5, 1.0}}, {{0.5, -1.0}}, {{0.25, 0.5}}},
			{{{0, 1, 2}}, {{1, 0, 3}}, {{0, 1, 4}}});
		require(!non_manifold.validate(&error), "three-sided mesh edge was accepted");
		require(error.find("non-manifold") != std::string::npos,
			"three-sided edge rejection did not identify non-manifold incidence");

		UvConstraintMesh same_direction({
			{{0.0, 0.0}}, {{1.0, 0.0}}, {{0.5, 1.0}}, {{0.25, 0.5}}},
			{{{0, 1, 2}}, {{0, 1, 3}}});
		require(!same_direction.validate(&error),
			"same-direction incidence on a nominal interior edge was accepted");
		require(error.find("inconsistent adjacent triangle winding") != std::string::npos,
			"same-direction edge rejection did not identify the orientation invariant");

		UvConstraintMesh exhausted = square_mesh();
		const UvConstraintMesh before_exhaustion = exhausted;
		std::uint32_t allocator = UvVertex::no_topology_id - 1;
		const std::uint32_t allocator_before = allocator;
		require(!exhausted.assign_unclaimed_topology_ids(allocator, &error),
			"topology allocator exhaustion was accepted");
		require(error.find("allocator exhausted") != std::string::npos,
			"allocator exhaustion did not produce its diagnostic");
		require(allocator == allocator_before && same_mesh(exhausted, before_exhaustion),
			"failed topology allocation partially mutated the mesh or caller cursor");

		UvConstraintMesh last_available_id({
			{{0.0, 0.0}, 700}, {{1.0, 0.0}, 701},
			{{1.0, 1.0}, 702}, {{0.0, 1.0}}},
			{{{0, 1, 2}}, {{0, 2, 3}}});
		std::uint32_t last_allocator = UvVertex::no_topology_id - 1;
		require(last_available_id.assign_unclaimed_topology_ids(last_allocator, &error),
			"last available topology ID was not assignable: " + error);
		require(last_allocator == UvVertex::no_topology_id
			&& last_available_id.vertices()[3].topology_id == UvVertex::no_topology_id - 1,
			"topology allocator did not stop exactly at the reserved sentinel");

		UvConstraintMesh missing_sample_id = square_mesh();
		const UvConstraintMesh before_missing_id = missing_sample_id;
		const std::vector<UvConstraintSample> invalid_id_segment{
			{{0.2, 0.2}, UvVertex::no_topology_id}, {{0.8, 0.2}, 703}};
		require(!missing_sample_id.insert_polyline(invalid_id_segment, 64, &error),
			"reserved topology sentinel was accepted as a sample ID");
		require(error.find("no canonical topology ID") != std::string::npos
			&& same_mesh(missing_sample_id, before_missing_id),
			"reserved sample-ID rejection was not diagnostic and transactional");
	}
}

int main()
{
	try
	{
		t_junction_test();
		std::puts("[uv-constraint] interior T-junction: PASS");
		one_to_two_sampling_test();
		std::puts("[uv-constraint] one-to-two canonical sampling: PASS");
		hole_and_boundary_test();
		std::puts("[uv-constraint] holes and trim boundaries: PASS");
		deterministic_replay_test();
		std::puts("[uv-constraint] deterministic replay: PASS");
		exact_existing_vertex_provenance_test();
		std::puts("[uv-constraint] exact existing-vertex provenance: PASS");
		parameter_tolerance_test();
		std::puts("[uv-constraint] scale-aware parameter tolerances: PASS");
		geometric_validation_test();
		std::puts("[uv-constraint] geometric planar validation: PASS");
		malformed_input_test();
		std::puts("[uv-constraint] malformed input rejection: PASS");
		std::puts("UV constraint mesh: PASS");
		return 0;
	}
	catch (const std::exception& error)
	{
		std::fprintf(stderr, "UV constraint mesh: FAIL: %s\n", error.what());
		return 1;
	}
}
