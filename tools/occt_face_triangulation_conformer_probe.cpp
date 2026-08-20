#include "core/geometry/occt_face_triangulation_conformer.h"

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Surface.hxx>
#include <Poly_Triangle.hxx>
#include <Precision.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <gp.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
	using namespace paracfd::core;

	void require(bool condition, const std::string& message)
	{
		if (!condition) throw std::runtime_error(message);
	}

	void require_close(double actual, double expected, double tolerance,
		const std::string& message)
	{
		if (std::abs(actual - expected) > tolerance)
			throw std::runtime_error(message + ": actual=" + std::to_string(actual)
				+ " expected=" + std::to_string(expected));
	}

	TopoDS_Face rectangle_face(double u0, double u1, double v0, double v1)
	{
		BRepBuilderAPI_MakeFace maker(gp_Pln(gp::XOY()), u0, u1, v0, v1);
		require(maker.IsDone(), "could not manufacture rectangular plane face");
		return maker.Face();
	}

	TopoDS_Wire rectangle_wire(double x0, double y0, double x1, double y1,
		bool clockwise)
	{
		std::vector<gp_Pnt> points{{x0, y0, 0.0}, {x1, y0, 0.0},
			{x1, y1, 0.0}, {x0, y1, 0.0}};
		if (clockwise) std::reverse(points.begin(), points.end());
		BRepBuilderAPI_MakePolygon polygon;
		for (const gp_Pnt& point : points) polygon.Add(point);
		polygon.Close();
		require(polygon.IsDone(), "could not manufacture rectangular wire");
		return polygon.Wire();
	}

	TopoDS_Face face_with_square_hole()
	{
		BRepBuilderAPI_MakeFace maker(gp_Pln(gp::XOY()),
			rectangle_wire(0.0, 0.0, 10.0, 10.0, false), true);
		require(maker.IsDone(), "could not manufacture hole fixture outer face");
		maker.Add(rectangle_wire(4.0, 4.0, 6.0, 6.0, true));
		require(maker.IsDone(), "could not add hole fixture inner wire");
		return maker.Face();
	}

	occ::handle<Poly_Triangulation> make_triangulation(const TopoDS_Face& face,
		const std::vector<UvPoint>& uv,
		const std::vector<std::array<int, 3>>& triangles)
	{
		TopLoc_Location surface_location;
		const occ::handle<Geom_Surface> surface = BRep_Tool::Surface(face,
			surface_location);
		require(!surface.IsNull(), "manufactured face has no surface");
		const occ::handle<Poly_Triangulation> mesh = new Poly_Triangulation(
			static_cast<int>(uv.size()), static_cast<int>(triangles.size()), true, false);
		for (std::size_t index = 0; index < uv.size(); ++index)
		{
			gp_Pnt world = surface->Value(uv[index].u, uv[index].v);
			world.Transform(surface_location.Transformation());
			mesh->SetNode(static_cast<int>(index + 1), world);
			mesh->SetUVNode(static_cast<int>(index + 1),
				gp_Pnt2d(uv[index].u, uv[index].v));
		}
		for (std::size_t index = 0; index < triangles.size(); ++index)
			mesh->SetTriangle(static_cast<int>(index + 1), Poly_Triangle(
				triangles[index][0], triangles[index][1], triangles[index][2]));
		return mesh;
	}

	occ::handle<Poly_Triangulation> square_mesh(const TopoDS_Face& face,
		double u0 = 0.0, double u1 = 10.0, double v0 = 0.0, double v1 = 10.0)
	{
		// Deliberately mix one CW and one CCW input triangle. The conformer must
		// normalize them before UvConstraintMesh sees the source topology.
		return make_triangulation(face,
			{{u0, v0}, {u1, v0}, {u1, v1}, {u0, v1}},
			{{1, 3, 2}, {1, 3, 4}});
	}

	std::array<double, 3> plane_position(double u, double v)
	{
		return {u, v, 0.0};
	}

	OcctFaceConstraintSample plane_sample(double u, double v, std::uint32_t id)
	{
		return {{u, v}, id, plane_position(u, v), 1.0e-10};
	}

	std::uint32_t vertex_with_topology(const OcctConformedFaceMesh& mesh,
		std::uint32_t topology_id)
	{
		for (std::size_t index = 0; index < mesh.topology_ids.size(); ++index)
			if (mesh.topology_ids[index] == topology_id)
				return static_cast<std::uint32_t>(index);
		throw std::runtime_error("missing topology ID " + std::to_string(topology_id));
	}

	std::vector<std::uint32_t> vertices_with_topology(
		const OcctConformedFaceMesh& mesh, std::uint32_t topology_id)
	{
		std::vector<std::uint32_t> vertices;
		for (std::size_t index = 0; index < mesh.topology_ids.size(); ++index)
			if (mesh.topology_ids[index] == topology_id)
				vertices.push_back(static_cast<std::uint32_t>(index));
		return vertices;
	}

	std::uint32_t vertex_at_uv(const OcctConformedFaceMesh& mesh, double u, double v)
	{
		for (std::size_t index = 0; index < mesh.uv.size(); ++index)
			if (mesh.uv[index].u == u && mesh.uv[index].v == v)
				return static_cast<std::uint32_t>(index);
		throw std::runtime_error("missing exact output UV vertex");
	}

	bool edge_has_tag(const OcctConformedFaceMesh& mesh, std::uint32_t a,
		std::uint32_t b, std::uint64_t tag, bool* boundary = nullptr)
	{
		if (a > b) std::swap(a, b);
		for (const OcctConformedFaceEdge& edge : mesh.edges)
		{
			auto vertices = edge.vertices;
			if (vertices[0] > vertices[1]) std::swap(vertices[0], vertices[1]);
			if (vertices != std::array<std::uint32_t, 2>{a, b}) continue;
			if (boundary) *boundary = edge.boundary;
			return std::binary_search(edge.constraint_ids.begin(),
				edge.constraint_ids.end(), tag);
		}
		return false;
	}

	long double uv_area2(const OcctConformedFaceMesh& mesh,
		const std::array<std::uint32_t, 3>& triangle)
	{
		const UvPoint& a = mesh.uv[triangle[0]];
		const UvPoint& b = mesh.uv[triangle[1]];
		const UvPoint& c = mesh.uv[triangle[2]];
		return (static_cast<long double>(b.u) - a.u)
			* (static_cast<long double>(c.v) - a.v)
			- (static_cast<long double>(b.v) - a.v)
			* (static_cast<long double>(c.u) - a.u);
	}

	double uv_area(const OcctConformedFaceMesh& mesh)
	{
		long double area2 = 0.0;
		for (const auto& triangle : mesh.triangles) area2 += uv_area2(mesh, triangle);
		return static_cast<double>(0.5L * area2);
	}

	bool same_output(const OcctConformedFaceMesh& a,
		const OcctConformedFaceMesh& b);

	std::size_t boundary_loop_count(const OcctConformedFaceMesh& mesh)
	{
		std::map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
		for (const OcctConformedFaceEdge& edge : mesh.edges)
			if (edge.boundary)
			{
				adjacency[edge.vertices[0]].push_back(edge.vertices[1]);
				adjacency[edge.vertices[1]].push_back(edge.vertices[0]);
			}
		std::set<std::uint32_t> visited;
		std::size_t loops = 0;
		for (const auto& [start, neighbours] : adjacency)
		{
			require(neighbours.size() == 2, "output boundary is not a closed loop");
			if (visited.contains(start)) continue;
			++loops;
			std::vector<std::uint32_t> stack{start};
			while (!stack.empty())
			{
				const std::uint32_t vertex = stack.back();
				stack.pop_back();
				if (!visited.insert(vertex).second) continue;
				for (std::uint32_t neighbour : adjacency[vertex])
					if (!visited.contains(neighbour)) stack.push_back(neighbour);
			}
		}
		return loops;
	}

	bool point_inside_triangle(const UvPoint& p, const UvPoint& a,
		const UvPoint& b, const UvPoint& c)
	{
		const auto orient = [](const UvPoint& x, const UvPoint& y, const UvPoint& z)
		{
			return (y.u - x.u) * (z.v - x.v) - (y.v - x.v) * (z.u - x.u);
		};
		return orient(a, b, p) >= 0.0 && orient(b, c, p) >= 0.0
			&& orient(c, a, p) >= 0.0;
	}

	bool contains(const OcctConformedFaceMesh& mesh, const UvPoint& point)
	{
		for (const auto& triangle : mesh.triangles)
			if (point_inside_triangle(point, mesh.uv[triangle[0]],
				mesh.uv[triangle[1]], mesh.uv[triangle[2]])) return true;
		return false;
	}

	void planar_t_contact_test()
	{
		const TopoDS_Face face = rectangle_face(0.0, 10.0, 0.0, 10.0);
		const auto input = square_mesh(face);
		const std::vector<OcctFaceConstraintPolyline> constraints{
			{10, {plane_sample(1.0, 5.0, 100), plane_sample(9.0, 5.0, 101)}},
			{20, {plane_sample(5.0, 9.0, 103), plane_sample(5.0, 5.0, 102)}}};
		OcctConformedFaceMesh output;
		std::string error;
		require(conform_occt_face_triangulation(face, input, {}, constraints,
			output, error), "planar T conforming failed: " + error);
		const std::uint32_t left = vertex_with_topology(output, 100);
		const std::uint32_t right = vertex_with_topology(output, 101);
		const std::uint32_t junction = vertex_with_topology(output, 102);
		const std::uint32_t stem = vertex_with_topology(output, 103);
		require(edge_has_tag(output, left, junction, 10)
			&& edge_has_tag(output, junction, right, 10)
			&& edge_has_tag(output, stem, junction, 20),
			"planar T lost a constraint edge/tag");
		unsigned constrained_degree = 0;
		for (const OcctConformedFaceEdge& edge : output.edges)
			if ((edge.vertices[0] == junction || edge.vertices[1] == junction)
				&& !edge.constraint_ids.empty()) ++constrained_degree;
		require(constrained_degree == 3, "planar T junction does not have degree three");
		require_close(uv_area(output), 100.0, 1.0e-12,
			"planar T changed trimmed area");
		for (const auto& triangle : output.triangles)
			require(uv_area2(output, triangle) > 0.0L,
				"planar T output is not consistently CCW in UV");
		const std::size_t unclaimed = static_cast<std::size_t>(std::count(
			output.topology_ids.begin(), output.topology_ids.end(),
			UvVertex::no_topology_id));
		require(unclaimed == 4, "ordinary source vertices were globally assigned topology IDs");
	}

	void boundary_split_test()
	{
		const TopoDS_Face face = rectangle_face(0.0, 10.0, 0.0, 10.0);
		const std::vector<OcctFaceConstraintPolyline> constraints{
			{30, {plane_sample(0.0, 0.0, 200), plane_sample(10.0, 0.0, 202)}},
			{31, {plane_sample(5.0, 0.0, 201), plane_sample(5.0, 4.0, 203)}}};
		OcctConformedFaceMesh output;
		std::string error;
		require(conform_occt_face_triangulation(face, square_mesh(face), {},
			constraints, output, error), "boundary split failed: " + error);
		const std::uint32_t a = vertex_with_topology(output, 200);
		const std::uint32_t split = vertex_with_topology(output, 201);
		const std::uint32_t b = vertex_with_topology(output, 202);
		const std::uint32_t interior = vertex_with_topology(output, 203);
		bool a_boundary = false, b_boundary = false;
		require(edge_has_tag(output, a, split, 30, &a_boundary)
			&& edge_has_tag(output, split, b, 30, &b_boundary)
			&& a_boundary && b_boundary,
			"splitting a tagged boundary edge did not preserve its tag/boundary state");
		require(edge_has_tag(output, split, interior, 31),
			"boundary contact stem lost its tag");
	}

	void exact_source_node_provenance_test()
	{
		const TopoDS_Face face = rectangle_face(0.0, 10.0, 0.0, 10.0);
		// The first sample's face-chart coordinate is deliberately offset from source
		// node zero.  Its independent 3-D certificate proves that it represents that
		// original node; UV insertion must claim the node by identity rather than ask
		// the near-junction coordinate classifier to choose among incident edges.
		OcctFaceConstraintSample first{{2.0e-7, 0.9e-7}, 220,
			plane_position(0.0, 0.0), 3.0e-7, 0};
		OcctFaceConstraintSample last = plane_sample(10.0, 0.0, 221);
		last.source_node_index = 1;
		const std::vector<OcctFaceConstraintPolyline> constraints{
			{32, {first, last}}};
		OcctConformedFaceMesh output;
		std::string error;
		const bool conformed = conform_occt_face_triangulation(face,
			square_mesh(face), {}, constraints, output, error);
		require(conformed,
			"exact original-node conforming failed: " + error);
		const std::uint32_t a = vertex_with_topology(output, 220);
		const std::uint32_t b = vertex_with_topology(output, 221);
		require(a == 0 && b == 1 && edge_has_tag(output, a, b, 32),
			"source-node provenance did not claim/tag the exact original edge");
		require(output.uv[a] == UvPoint{0.0, 0.0}
			&& output.world_positions[a] == plane_position(0.0, 0.0),
			"source-node provenance moved the authoritative original vertex");
	}

	void hole_preservation_test()
	{
		const TopoDS_Face face = face_with_square_hole();
		const std::vector<UvPoint> uv{{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0},
			{0.0, 10.0}, {4.0, 4.0}, {6.0, 4.0}, {6.0, 6.0}, {4.0, 6.0}};
		const std::vector<std::array<int, 3>> triangles{
			{1, 2, 6}, {1, 6, 5}, {2, 3, 7}, {2, 7, 6},
			{3, 4, 8}, {3, 8, 7}, {4, 1, 5}, {4, 5, 8}};
		const std::vector<OcctFaceConstraintPolyline> constraints{
			{40, {plane_sample(5.0, 0.0, 300), plane_sample(5.0, 4.0, 301)}}};
		OcctConformedFaceMesh output;
		std::string error;
		require(conform_occt_face_triangulation(face,
			make_triangulation(face, uv, triangles), {}, constraints, output, error),
			"hole conforming failed: " + error);
		require(boundary_loop_count(output) == 2,
			"conforming changed outer+hole boundary loop count");
		require_close(uv_area(output), 96.0, 1.0e-12,
			"conforming filled or clipped the hole");
		require(!contains(output, {5.0, 5.0}),
			"conforming filled the triangulation's hole");
	}

	void reversed_face_test()
	{
		const TopoDS_Face forward = rectangle_face(0.0, 10.0, 0.0, 10.0);
		TopoDS_Face reversed = forward;
		reversed.Reverse();
		const auto input = square_mesh(forward);
		OcctConformedFaceMesh a, b;
		std::string error;
		require(conform_occt_face_triangulation(forward, input, {}, {}, a, error),
			"forward face conforming failed: " + error);
		require(conform_occt_face_triangulation(reversed, input, {}, {}, b, error),
			"reversed face conforming failed: " + error);
		require(!a.face_reversed && b.face_reversed,
			"face orientation parity was not retained explicitly");
		require(a.uv == b.uv && a.world_positions == b.world_positions
			&& a.triangles == b.triangles && a.topology_ids == b.topology_ids,
			"reversing a face changed its normalized UV triangulation");
		for (const auto& triangle : b.triangles)
			require(uv_area2(b, triangle) > 0.0L,
				"reversed face output is not normalized CCW in UV");
	}

	void curved_round_trip_test()
	{
		const occ::handle<Geom_CylindricalSurface> cylinder =
			new Geom_CylindricalSurface(gp_Ax3(gp::Origin(), gp::DZ()), 2.0);
		BRepBuilderAPI_MakeFace maker(cylinder, 0.0, 1.5, 0.0, 5.0,
			Precision::Confusion());
		require(maker.IsDone(), "could not manufacture cylindrical face");
		const TopoDS_Face face = maker.Face();
		const auto input = square_mesh(face, 0.0, 1.5, 0.0, 5.0);
		const auto exact_position = [&](double u, double v)
		{
			const gp_Pnt p = cylinder->Value(u, v);
			return std::array<double, 3>{p.X(), p.Y(), p.Z()};
		};
		const std::vector<OcctFaceConstraintPolyline> constraints{{50, {
			{{0.2, 2.0}, 400, exact_position(0.2, 2.0), 1.0e-10},
			{{1.2, 2.0}, 401, exact_position(1.2, 2.0), 1.0e-10}}}};
		OcctConformedFaceMesh output;
		std::string error;
		require(conform_occt_face_triangulation(face, input, {}, constraints,
			output, error), "curved conforming failed: " + error);
		for (std::uint32_t id : {400u, 401u})
		{
			const auto& p = output.world_positions[vertex_with_topology(output, id)];
			require_close(std::hypot(p[0], p[1]), 2.0, 2.0e-14,
				"curved topology vertex is not on the exact cylinder");
		}
		for (std::size_t vertex = 0; vertex < output.uv.size(); ++vertex)
		{
			const gp_Pnt expected = cylinder->Value(output.uv[vertex].u,
				output.uv[vertex].v);
			const gp_Pnt actual(output.world_positions[vertex][0],
				output.world_positions[vertex][1], output.world_positions[vertex][2]);
			require(expected.Distance(actual) <= 1.0e-10,
				"curved output vertex did not round-trip through the exact surface");
		}
	}

	void periodic_cylinder_seam_test()
	{
		const double period = 2.0 * std::acos(-1.0);
		const occ::handle<Geom_CylindricalSurface> cylinder =
			new Geom_CylindricalSurface(gp_Ax3(gp::Origin(), gp::DZ()), 2.0);
		BRepBuilderAPI_MakeFace maker(cylinder, 0.0, period, 0.0, 5.0,
			Precision::Confusion());
		require(maker.IsDone(), "could not manufacture periodic cylinder face");
		const TopoDS_Face face = maker.Face();
		const auto sample = [&](double chart_u, double v, std::uint32_t id)
		{
			const gp_Pnt canonical = cylinder->Value(0.0, v);
			return OcctFaceConstraintSample{{chart_u, v}, id,
				{canonical.X(), canonical.Y(), canonical.Z()}, 1.0e-10};
		};
		const std::vector<OcctFaceConstraintPolyline> constraints{
			{90, {sample(0.0, 0.0, 800), sample(0.0, 2.5, 801),
				sample(0.0, 5.0, 802)}, 9000},
			{90, {sample(period, 0.0, 800), sample(period, 2.5, 801),
				sample(period, 5.0, 802)}, 9001}};
		OcctConformedFaceMesh output;
		std::string error;
		require(conform_occt_face_triangulation(face,
			square_mesh(face, 0.0, period, 0.0, 5.0), {}, constraints,
			output, error), "periodic cylinder seam conforming failed: " + error);

		for (std::uint32_t id : {800u, 801u, 802u})
		{
			const std::vector<std::uint32_t> aliases =
				vertices_with_topology(output, id);
			require(aliases.size() == 2,
				"periodic seam did not emit two occurrences of physical topology ID "
				+ std::to_string(id));
			const UvPoint& a_uv = output.uv[aliases[0]];
			const UvPoint& b_uv = output.uv[aliases[1]];
			require((a_uv.u == 0.0 && b_uv.u == period)
				|| (a_uv.u == period && b_uv.u == 0.0),
				"periodic topology aliases lost their distinct UV occurrences");
			for (unsigned axis = 0; axis < 3; ++axis)
				require(std::bit_cast<std::uint64_t>(
					output.world_positions[aliases[0]][axis])
					== std::bit_cast<std::uint64_t>(
						output.world_positions[aliases[1]][axis]),
					"periodic topology aliases are not bit-identical in world space");
		}
		for (double u : {0.0, period})
		{
			const std::uint32_t a = vertex_at_uv(output, u, 0.0);
			const std::uint32_t split = vertex_at_uv(output, u, 2.5);
			const std::uint32_t b = vertex_at_uv(output, u, 5.0);
			bool first_boundary = false, second_boundary = false;
			require(edge_has_tag(output, a, split, 90, &first_boundary)
				&& edge_has_tag(output, split, b, 90, &second_boundary)
				&& first_boundary && second_boundary,
				"periodic seam branch lost its boundary constraint chain");
		}
		require_close(uv_area(output), period * 5.0, 1.0e-12,
			"periodic seam insertion changed chart area");
		for (const auto& triangle : output.triangles)
			require(uv_area2(output, triangle) > 0.0L,
				"periodic seam output is not consistently CCW in UV");
	}

	void malformed_periodic_alias_test()
	{
		const double period = 2.0 * std::acos(-1.0);
		const occ::handle<Geom_CylindricalSurface> cylinder =
			new Geom_CylindricalSurface(gp_Ax3(gp::Origin(), gp::DZ()), 2.0);
		BRepBuilderAPI_MakeFace maker(cylinder, 0.0, period, 0.0, 5.0,
			Precision::Confusion());
		require(maker.IsDone(), "could not manufacture malformed-alias fixture");
		const TopoDS_Face face = maker.Face();
		const auto input = square_mesh(face, 0.0, period, 0.0, 5.0);
		const auto sample = [&](double chart_u, double v, std::uint32_t id,
			double tolerance = 1.0e-10)
		{
			const gp_Pnt canonical = cylinder->Value(0.0, v);
			return OcctFaceConstraintSample{{chart_u, v}, id,
				{canonical.X(), canonical.Y(), canonical.Z()}, tolerance};
		};
		OcctConformedFaceMesh output;
		output.uv = {{91.0, 92.0}};
		output.world_positions = {{{93.0, 94.0, 95.0}}};
		output.triangles = {{{7, 8, 9}}};
		output.topology_ids = {77};
		output.edges = {{{7, 8}, true, {1234}}};
		output.face_reversed = true;
		output.source_node_count = 987;
		const OcctConformedFaceMesh sentinel = output;
		std::string error;
		const auto expect_rejection = [&](const auto& constraints,
			const std::string& diagnostic, const std::string& description)
		{
			require(!conform_occt_face_triangulation(face, input, {}, constraints,
				output, error) && error.find(diagnostic) != std::string::npos,
				description + " was not rejected diagnostically: " + error);
			require(same_output(output, sentinel),
				description + " changed caller output");
		};

		const std::vector<OcctFaceConstraintPolyline> implicit_alias{
			{91, {sample(0.0, 1.0, 810), sample(0.0, 4.0, 811)}},
			{91, {sample(period, 1.0, 810), sample(period, 4.0, 811)}}};
		expect_rejection(implicit_alias, "explicit chart branch",
			"implicit periodic alias");

		const std::vector<OcctFaceConstraintPolyline> reused_branch{
			{91, {sample(0.0, 1.0, 810), sample(0.0, 4.0, 811)}, 77},
			{91, {sample(period, 1.0, 810), sample(period, 4.0, 811)}, 77}};
		expect_rejection(reused_branch, "reuses one chart branch",
			"reused periodic branch identity");

		const std::vector<OcctFaceConstraintPolyline> missing_branch{
			{91, {sample(0.0, 1.0, 810), sample(0.0, 4.0, 811)}, 77},
			{91, {sample(period, 1.0, 810), sample(period, 4.0, 811)}}};
		expect_rejection(missing_branch, "explicit chart branch",
			"partially named periodic alias");

		const gp_Pnt shared = cylinder->Value(0.0, 1.0);
		const OcctFaceConstraintSample offset_alias{{period, 1.01}, 820,
			{shared.X(), shared.Y(), shared.Z()}, 0.02};
		const std::vector<OcctFaceConstraintPolyline> nonperiodic_offset{
			{92, {sample(0.0, 1.0, 820), sample(0.0, 4.0, 821)}, 78},
			{92, {offset_alias, sample(period, 4.0, 822)}, 79}};
		expect_rejection(nonperiodic_offset, "not periodic aliases",
			"explicit nonperiodic chart offset");
	}

	void canonical_cross_face_test()
	{
		const TopoDS_Face left_face = rectangle_face(0.0, 10.0, 0.0, 10.0);
		const TopoDS_Face right_face = rectangle_face(10.0, 20.0, 0.0, 10.0);
		const std::vector<OcctFaceConstraintPolyline> shared{{60, {
			plane_sample(10.0, 0.0, 500), plane_sample(10.0, 5.0, 501),
			plane_sample(10.0, 10.0, 502)}}};
		OcctConformedFaceMesh left, right;
		std::string error;
		require(conform_occt_face_triangulation(left_face, square_mesh(left_face), {},
			shared, left, error), "left canonical face failed: " + error);
		require(conform_occt_face_triangulation(right_face,
			square_mesh(right_face, 10.0, 20.0, 0.0, 10.0), {}, shared,
			right, error), "right canonical face failed: " + error);
		for (std::uint32_t id : {500u, 501u, 502u})
		{
			const auto& a = left.world_positions[vertex_with_topology(left, id)];
			const auto& b = right.world_positions[vertex_with_topology(right, id)];
			for (unsigned axis = 0; axis < 3; ++axis)
				require(std::bit_cast<std::uint64_t>(a[axis])
					== std::bit_cast<std::uint64_t>(b[axis]),
					"canonical cross-face position is not bit-identical");
		}
	}

	void tolerance_isolation_test()
	{
		const TopoDS_Face face = rectangle_face(0.0, 10.0, 0.0, 10.0);
		const std::vector<OcctFaceConstraintPolyline> constraints{
			{70, {
				{{1.0, 2.0}, 700, {1.0, 2.0, 1.0}, 1.0},
				{{9.0, 2.0}, 701, {9.0, 2.0, 1.0}, 1.0}}},
			{71, {
				plane_sample(1.0, 7.0, 702),
				{{9.0, 7.0}, 703, {9.0, 7.0, 0.5}, 0.0}}}};
		OcctConformedFaceMesh output;
		std::string error;
		require(!conform_occt_face_triangulation(face, square_mesh(face), {},
			constraints, output, error),
			"one loose contact tolerance validated an unrelated strict contact");
		require(error.find("polyline 1 sample 1") != std::string::npos
			&& error.find("misses") != std::string::npos,
			"isolated-tolerance failure did not identify the strict sample: " + error);
	}

	TopoDS_Face anisotropic_face(occ::handle<Geom_BezierSurface>& surface)
	{
		TColgp_Array2OfPnt poles(1, 2, 1, 2);
		poles.SetValue(1, 1, gp_Pnt(0.0, 0.0, 0.0));
		poles.SetValue(2, 1, gp_Pnt(1000.0, 0.0, 0.0));
		poles.SetValue(1, 2, gp_Pnt(0.0, 1.0, 0.0));
		poles.SetValue(2, 2, gp_Pnt(1000.0, 1.0, 0.0));
		surface = new Geom_BezierSurface(poles);
		BRepBuilderAPI_MakeFace maker(surface, 0.0, 1.0, 0.0, 1.0,
			Precision::Confusion());
		require(maker.IsDone(), "could not manufacture anisotropic Bezier face");
		return maker.Face();
	}

	TopoDS_Face reparameterized_plane(double u0, double u1, double v0, double v1,
		occ::handle<Geom_BSplineSurface>& surface)
	{
		TColgp_Array2OfPnt poles(1, 2, 1, 2);
		poles.SetValue(1, 1, gp_Pnt(0.0, 0.0, 0.0));
		poles.SetValue(2, 1, gp_Pnt(10.0, 0.0, 0.0));
		poles.SetValue(1, 2, gp_Pnt(0.0, 10.0, 0.0));
		poles.SetValue(2, 2, gp_Pnt(10.0, 10.0, 0.0));
		TColStd_Array1OfReal u_knots(1, 2), v_knots(1, 2);
		u_knots.SetValue(1, u0);
		u_knots.SetValue(2, u1);
		v_knots.SetValue(1, v0);
		v_knots.SetValue(2, v1);
		TColStd_Array1OfInteger u_multiplicities(1, 2), v_multiplicities(1, 2);
		u_multiplicities.Init(2);
		v_multiplicities.Init(2);
		surface = new Geom_BSplineSurface(poles, u_knots, v_knots,
			u_multiplicities, v_multiplicities, 1, 1);
		BRepBuilderAPI_MakeFace maker(surface, u0, u1, v0, v1,
			Precision::Confusion());
		require(maker.IsDone(), "could not manufacture reparameterized plane face");
		return maker.Face();
	}

	void conform_reparameterized_face(double u0, double u1, double v0, double v1,
		std::uint32_t first_id, std::uint64_t constraint_id,
		const std::string& description)
	{
		occ::handle<Geom_BSplineSurface> surface;
		const TopoDS_Face face = reparameterized_plane(u0, u1, v0, v1, surface);
		const double u_mid = u0 + 0.5 * (u1 - u0);
		const auto sample = [&](double u, std::uint32_t id)
		{
			const gp_Pnt p = surface->Value(u, v0);
			return OcctFaceConstraintSample{{u, v0}, id,
				{p.X(), p.Y(), p.Z()}, 0.0};
		};
		const std::vector<OcctFaceConstraintPolyline> constraints{{constraint_id, {
			sample(u0, first_id), sample(u_mid, first_id + 1),
			sample(u1, first_id + 2)}}};
		OcctConformedFaceMesh output;
		std::string error;
		require(conform_occt_face_triangulation(face,
			square_mesh(face, u0, u1, v0, v1), {}, constraints, output, error),
			description + " conforming failed: " + error);
		const std::uint32_t a = vertex_with_topology(output, first_id);
		const std::uint32_t split = vertex_with_topology(output, first_id + 1);
		const std::uint32_t b = vertex_with_topology(output, first_id + 2);
		bool first_boundary = false, second_boundary = false;
		require(edge_has_tag(output, a, split, constraint_id, &first_boundary)
			&& edge_has_tag(output, split, b, constraint_id, &second_boundary)
			&& first_boundary && second_boundary,
			description + " collapsed or lost its reparameterized boundary chain");
		for (const auto& triangle : output.triangles)
			require(uv_area2(output, triangle) > 0.0L,
				description + " output is not consistently CCW in UV");
	}

	void reparameterized_uv_tolerance_test()
	{
		const double tiny_u0 = 1.0;
		const double tiny_u1 = tiny_u0 + 1.0e-12;
		const double tiny_v0 = -2.0;
		const double tiny_v1 = tiny_v0 + 2.0e-12;
		conform_reparameterized_face(tiny_u0, tiny_u1, tiny_v0, tiny_v1,
			730, 73, "tiny-span UV chart");

		const double offset_u0 = 1.0e8;
		const double offset_u1 = offset_u0 + 2.0e-5;
		const double offset_v0 = -1.0e8;
		const double offset_v1 = offset_v0 + 4.0e-5;
		conform_reparameterized_face(offset_u0, offset_u1, offset_v0, offset_v1,
			740, 74, "large-offset/small-span UV chart");

		const double unresolved_u0 = 1.0e12;
		const double unresolved_u1 = unresolved_u0 + 1.0e-3;
		occ::handle<Geom_BSplineSurface> unresolved_surface;
		const TopoDS_Face unresolved_face = reparameterized_plane(unresolved_u0,
			unresolved_u1, 0.0, 1.0, unresolved_surface);
		OcctConformedFaceMesh output;
		output.uv = {{91.0, 92.0}};
		output.world_positions = {{{93.0, 94.0, 95.0}}};
		output.triangles = {{{7, 8, 9}}};
		output.topology_ids = {77};
		output.edges = {{{7, 8}, true, {1234}}};
		output.face_reversed = true;
		output.source_node_count = 987;
		const OcctConformedFaceMesh sentinel = output;
		std::string error;
		require(!conform_occt_face_triangulation(unresolved_face,
			square_mesh(unresolved_face, unresolved_u0, unresolved_u1, 0.0, 1.0),
			{}, {}, output, error) && error.find("numerically unresolved")
			!= std::string::npos,
			"materially under-resolved UV chart was not rejected: " + error);
		require(same_output(output, sentinel),
			"under-resolved UV chart changed caller output");
	}

	void anisotropic_uv_tolerance_test()
	{
		occ::handle<Geom_BezierSurface> surface;
		const TopoDS_Face face = anisotropic_face(surface);
		constexpr double split_u = 5.0e-8;
		const std::vector<UvPoint> uv{{0.0, 0.0}, {split_u, 0.0}, {1.0, 0.0},
			{1.0, 1.0}, {0.0, 1.0}};
		const std::vector<std::array<int, 3>> triangles{
			{5, 1, 2}, {5, 2, 3}, {5, 3, 4}};
		const auto sample = [&](double u, std::uint32_t id)
		{
			const gp_Pnt p = surface->Value(u, 0.0);
			return OcctFaceConstraintSample{{u, 0.0}, id,
				{p.X(), p.Y(), p.Z()}, 0.0};
		};
		const std::vector<OcctFaceConstraintPolyline> constraints{{72, {
			sample(0.0, 710), sample(split_u, 711), sample(1.0, 712)}}};
		OcctConformedFaceMesh output;
		std::string error;
		require(conform_occt_face_triangulation(face,
			make_triangulation(face, uv, triangles), {}, constraints, output, error),
			"anisotropic chart collapsed sensitive-axis vertices: " + error);
		const std::uint32_t a = vertex_with_topology(output, 710);
		const std::uint32_t split = vertex_with_topology(output, 711);
		const std::uint32_t b = vertex_with_topology(output, 712);
		require(edge_has_tag(output, a, split, 72)
			&& edge_has_tag(output, split, b, 72),
			"anisotropic boundary constraint lost its split/tag topology");
	}

	void sentinel_boundary_tag_test()
	{
		const TopoDS_Face face = rectangle_face(0.0, 10.0, 0.0, 10.0);
		const double nan = std::numeric_limits<double>::quiet_NaN();
		const auto unclaimed = [&](double u, double v)
		{
			return OcctFaceConstraintSample{{u, v}, UvVertex::no_topology_id,
				{nan, nan, nan}, nan};
		};
		const std::vector<OcctFaceConstraintPolyline> constraints{
			{80, {unclaimed(0.0, 0.0), unclaimed(10.0, 0.0)}},
			{81, {unclaimed(5.0, 0.0), unclaimed(5.0, 4.0)}}};
		OcctConformedFaceMesh output;
		std::string error;
		require(conform_occt_face_triangulation(face, square_mesh(face), {},
			constraints, output, error),
			"sentinel boundary tagging failed: " + error);
		require(std::all_of(output.topology_ids.begin(), output.topology_ids.end(),
			[](std::uint32_t id) { return id == UvVertex::no_topology_id; }),
			"temporary face-local IDs leaked into conformer output");
		const std::uint32_t a = vertex_at_uv(output, 0.0, 0.0);
		const std::uint32_t split = vertex_at_uv(output, 5.0, 0.0);
		const std::uint32_t b = vertex_at_uv(output, 10.0, 0.0);
		const std::uint32_t stem = vertex_at_uv(output, 5.0, 4.0);
		require(edge_has_tag(output, a, split, 80)
			&& edge_has_tag(output, split, b, 80)
			&& edge_has_tag(output, split, stem, 81),
			"unclaimed boundary split did not retain its constraint tags");
	}

	bool same_output(const OcctConformedFaceMesh& a,
		const OcctConformedFaceMesh& b)
	{
		if (a.uv.size() != b.uv.size()) return false;
		for (std::size_t i = 0; i < a.uv.size(); ++i)
			if (!(a.uv[i] == b.uv[i])) return false;
		if (a.world_positions != b.world_positions || a.triangles != b.triangles
			|| a.topology_ids != b.topology_ids || a.face_reversed != b.face_reversed
			|| a.source_node_count != b.source_node_count || a.edges.size() != b.edges.size())
			return false;
		for (std::size_t i = 0; i < a.edges.size(); ++i)
			if (a.edges[i].vertices != b.edges[i].vertices
				|| a.edges[i].boundary != b.edges[i].boundary
				|| a.edges[i].constraint_ids != b.edges[i].constraint_ids) return false;
		return true;
	}

	void transactional_rejection_test()
	{
		const TopoDS_Face face = rectangle_face(0.0, 10.0, 0.0, 10.0);
		OcctConformedFaceMesh output;
		output.uv = {{91.0, 92.0}};
		output.world_positions = {{{93.0, 94.0, 95.0}}};
		output.triangles = {{{7, 8, 9}}};
		output.topology_ids = {77};
		output.edges = {{{7, 8}, true, {1234}}};
		output.face_reversed = true;
		output.source_node_count = 987;
		const OcctConformedFaceMesh sentinel = output;
		std::string error;

		const std::vector<OcctFaceConstraintPolyline> off_surface{{70, {
			{{2.0, 2.0}, 600, {2.0, 2.0, 1.0}, 0.0},
			plane_sample(8.0, 2.0, 601)}}};
		require(!conform_occt_face_triangulation(face, square_mesh(face), {},
			off_surface, output, error) && error.find("misses") != std::string::npos,
			"off-surface canonical point was not rejected diagnostically");
		require(same_output(output, sentinel),
			"off-surface failure changed caller output");

		const std::vector<OcctFaceConstraintPolyline> conflict{{71, {
			plane_sample(2.0, 2.0, 610), plane_sample(8.0, 2.0, 610)}}};
		require(!conform_occt_face_triangulation(face, square_mesh(face), {},
			conflict, output, error) && error.find("conflicting") != std::string::npos,
			"conflicting canonical topology position was not rejected");
		require(same_output(output, sentinel),
			"canonical conflict changed caller output");

		const occ::handle<Poly_Triangulation> missing_uv =
			new Poly_Triangulation(4, 2, false, false);
		require(!conform_occt_face_triangulation(face, missing_uv, {}, {}, output, error)
			&& error.find("no UV nodes") != std::string::npos,
			"missing UV data was not rejected");
		require(same_output(output, sentinel), "missing UV failure changed caller output");

		const TopoDS_Face overlap_face = rectangle_face(0.0, 3.0, 0.0, 3.0);
		const auto overlap = make_triangulation(overlap_face,
			{{0.0, 0.0}, {3.0, 0.0}, {0.0, 3.0},
				{0.5, 0.5}, {1.0, 0.5}, {0.5, 1.0}},
			{{1, 2, 3}, {4, 5, 6}});
		require(!conform_occt_face_triangulation(overlap_face, overlap, {}, {},
			output, error) && error.find("overlap") != std::string::npos,
			"positive-area overlapping UV triangles were not rejected");
		require(same_output(output, sentinel), "overlap failure changed caller output");

		const TopoDS_Face non_manifold_face = rectangle_face(0.0, 3.0, -2.0, 2.0);
		const auto non_manifold = make_triangulation(non_manifold_face,
			{{0.0, 0.0}, {3.0, 0.0}, {1.0, 1.0}, {1.0, -1.0}, {2.0, 1.0}},
			{{1, 2, 3}, {2, 1, 4}, {1, 2, 5}});
		require(!conform_occt_face_triangulation(non_manifold_face, non_manifold,
			{}, {}, output, error) && error.find("non-manifold") != std::string::npos,
			"three-triangle UV edge incidence was not rejected as non-manifold");
		require(same_output(output, sentinel),
			"non-manifold failure changed caller output");
	}
}

int main()
{
	try
	{
		planar_t_contact_test();
		std::puts("[occt-face-conformer] planar interior T contact: PASS");
		boundary_split_test();
		std::puts("[occt-face-conformer] boundary split/tag inheritance: PASS");
		exact_source_node_provenance_test();
		std::puts("[occt-face-conformer] exact original-node provenance: PASS");
		hole_preservation_test();
		std::puts("[occt-face-conformer] trimmed hole preservation: PASS");
		reversed_face_test();
		std::puts("[occt-face-conformer] reversed face orientation: PASS");
		curved_round_trip_test();
		std::puts("[occt-face-conformer] curved exact-surface round trip: PASS");
		periodic_cylinder_seam_test();
		std::puts("[occt-face-conformer] periodic cylinder seam alias: PASS");
		canonical_cross_face_test();
		std::puts("[occt-face-conformer] cross-face canonical identity: PASS");
		tolerance_isolation_test();
		std::puts("[occt-face-conformer] per-contact tolerance isolation: PASS");
		anisotropic_uv_tolerance_test();
		std::puts("[occt-face-conformer] anisotropic UV tolerance: PASS");
		reparameterized_uv_tolerance_test();
		std::puts("[occt-face-conformer] reparameterized UV tolerance/resolution: PASS");
		sentinel_boundary_tag_test();
		std::puts("[occt-face-conformer] unclaimed boundary tagging: PASS");
		transactional_rejection_test();
		std::puts("[occt-face-conformer] transactional rejection: PASS");
		malformed_periodic_alias_test();
		std::puts("[occt-face-conformer] malformed periodic aliases: PASS");
		std::puts("OCCT face triangulation conformer: PASS");
		return 0;
	}
	catch (const std::exception& exception)
	{
		std::fprintf(stderr, "OCCT face triangulation conformer: FAIL: %s\n",
			exception.what());
		return 1;
	}
}
