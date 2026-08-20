#include "core/geometry/occt_conforming_mesh_builder.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Precision.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp.hxx>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
	using namespace paracfd::core;

	void require(bool condition, const std::string& message)
	{
		if (!condition) throw std::runtime_error(message);
	}

	TopoDS_Compound compound(const std::vector<TopoDS_Face>& faces)
	{
		BRep_Builder builder;
		TopoDS_Compound result;
		builder.MakeCompound(result);
		for (const auto& face : faces) builder.Add(result, face);
		return result;
	}

	OcctContactTopology make_topology(const TopoDS_Shape& shape,
		const std::vector<TopoDS_Face>& faces)
	{
		BRepMesh_IncrementalMesh mesher(shape, 0.01, false, 0.20, true);
		require(mesher.IsDone(), "OCCT did not mesh the manufactured shape");
		OcctContactTopology result;
		std::string error;
		const bool built = build_occt_contact_topology(shape, faces, {}, result, error);
		require(built,
			"contact topology fixture failed: " + error);
		return result;
	}

	std::vector<TopoDS_Face> faces_of(const TopoDS_Shape& shape)
	{
		std::vector<TopoDS_Face> result;
		for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
			result.push_back(TopoDS::Face(face.Current()));
		return result;
	}

	bool share_edge(const TopoDS_Face& a, const TopoDS_Face& b)
	{
		for (TopExp_Explorer ea(a, TopAbs_EDGE); ea.More(); ea.Next())
			for (TopExp_Explorer eb(b, TopAbs_EDGE); eb.More(); eb.Next())
				if (ea.Current().IsSame(eb.Current())) return true;
		return false;
	}

	std::vector<TopoDS_Face> adjacent_box_faces(const TopoDS_Shape& box)
	{
		const auto faces = faces_of(box);
		for (std::size_t a = 0; a < faces.size(); ++a)
			for (std::size_t b = a + 1; b < faces.size(); ++b)
				if (share_edge(faces[a], faces[b])) return {faces[a], faces[b]};
		throw std::runtime_error("box fixture has no adjacent face pair");
	}

	std::size_t boundary_loop_count(const OcctConformingFaceMesh& face)
	{
		std::map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
		for (const auto& edge : face.edge_provenance)
			if (edge.boundary)
			{
				adjacency[edge.vertices[0]].push_back(edge.vertices[1]);
				adjacency[edge.vertices[1]].push_back(edge.vertices[0]);
			}
		std::set<std::uint32_t> visited;
		std::size_t loops = 0;
		for (const auto& [start, neighbours] : adjacency)
		{
			require(neighbours.size() == 2,
				"manufactured face boundary is not a two-regular loop graph");
			if (visited.contains(start)) continue;
			++loops;
			std::vector<std::uint32_t> stack{start};
			while (!stack.empty())
			{
				const std::uint32_t vertex = stack.back();
				stack.pop_back();
				if (!visited.insert(vertex).second) continue;
				for (const auto next : adjacency[vertex]) stack.push_back(next);
			}
		}
		return loops;
	}

	bool bits_equal(const std::array<double, 3>& a,
		const std::array<double, 3>& b)
	{
		for (std::size_t axis = 0; axis < 3; ++axis)
			if (std::bit_cast<std::uint64_t>(a[axis])
				!= std::bit_cast<std::uint64_t>(b[axis])) return false;
		return true;
	}

	const OcctConformingFaceMesh& output_face(const OcctConformingMesh& mesh,
		std::uint32_t face_id)
	{
		const auto found = std::find_if(mesh.faces.begin(), mesh.faces.end(),
			[&](const auto& face) { return face.source_face_id == face_id; });
		require(found != mesh.faces.end(), "missing conformed source face");
		return *found;
	}

	void tag_namespace_zero_collision()
	{
		constexpr std::uint64_t cad = occt_cad_edge_constraint_id(0);
		constexpr auto atom = occt_contact_atom_constraint_id(0);
		static_assert(atom.has_value());
		static_assert(cad != *atom);
		static_assert(occt_cad_edge_from_constraint(cad).value() == 0);
		static_assert(occt_contact_atom_from_constraint(*atom).value() == 0);
		static_assert(!occt_contact_atom_from_constraint(cad));
		static_assert(!occt_cad_edge_from_constraint(*atom));
		require(!occt_contact_atom_constraint_id(occt_constraint_payload_mask + 1),
			"out-of-range atom tag payload was accepted");
	}

	void ordinary_shared_boundary_and_fan_one()
	{
		const TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 8.0, 6.0).Shape();
		const auto faces = adjacent_box_faces(box);
		const auto topology = make_topology(box, faces);
		const auto shared = std::find_if(topology.atomization.atoms.begin(),
			topology.atomization.atoms.end(), [](const auto& atom)
			{
				return atom.fan_degree == 2 && atom.face_uses.size() == 2;
			});
		require(shared != topology.atomization.atoms.end(),
			"adjacent box faces did not produce a shared boundary atom");
		const auto open = std::find_if(topology.atomization.atoms.begin(),
			topology.atomization.atoms.end(), [](const auto& atom)
			{
				return atom.is_open_boundary();
			});
		require(open != topology.atomization.atoms.end(),
			"adjacent-face fixture did not retain an explicit fan-one boundary atom");

		OcctConformingMesh output;
		std::string error;
		const bool built = build_occt_conforming_mesh(topology, output, error);
		require(built,
			"ordinary shared-boundary build failed: " + error);
		std::vector<std::array<double, 3>> shared_positions;
		bool saw_fan_one_tag = false;
		for (const auto& face : output.faces)
		{
			for (std::size_t vertex = 0; vertex < face.mesh.topology_ids.size(); ++vertex)
				for (const auto& sample : shared->samples)
					if (face.mesh.topology_ids[vertex] == sample.topology_id)
						shared_positions.push_back(face.mesh.world_positions[vertex]);
			for (const auto& edge : face.edge_provenance)
				if (edge.contact_atom_id == open->id)
				{
					require(edge.source_cad_edge_id.has_value() && edge.boundary,
						"fan-one boundary lost its CAD/open-boundary lineage");
					saw_fan_one_tag = true;
				}
		}
		require(shared_positions.size() >= 4,
			"shared boundary topology IDs did not occur on both faces");
		for (const auto& sample : shared->samples)
		{
			std::vector<std::array<double, 3>> copies;
			for (const auto& face : output.faces)
				for (std::size_t vertex = 0; vertex < face.mesh.topology_ids.size(); ++vertex)
					if (face.mesh.topology_ids[vertex] == sample.topology_id)
						copies.push_back(face.mesh.world_positions[vertex]);
			require(copies.size() == 2 && bits_equal(copies[0], copies[1]),
				"ordinary shared boundary is not bit-identical across faces");
		}
		require(saw_fan_one_tag, "fan-one atom tag was not published");
	}

	TopoDS_Face rectangle_xy(double x0, double x1, double y0, double y1,
		double z = 0.0)
	{
		BRepBuilderAPI_MakePolygon polygon;
		polygon.Add(gp_Pnt(x0, y0, z));
		polygon.Add(gp_Pnt(x1, y0, z));
		polygon.Add(gp_Pnt(x1, y1, z));
		polygon.Add(gp_Pnt(x0, y1, z));
		polygon.Close();
		return BRepBuilderAPI_MakeFace(polygon.Wire(), true).Face();
	}

	TopoDS_Face rectangle_xz(double x0, double x1, double y,
		double z0, double z1)
	{
		BRepBuilderAPI_MakePolygon polygon;
		polygon.Add(gp_Pnt(x0, y, z0));
		polygon.Add(gp_Pnt(x1, y, z0));
		polygon.Add(gp_Pnt(x1, y, z1));
		polygon.Add(gp_Pnt(x0, y, z1));
		polygon.Close();
		return BRepBuilderAPI_MakeFace(polygon.Wire(), true).Face();
	}

	std::uint32_t bottom_rib_edge(const OcctContactTopology& topology,
		double expected_z = 0.0)
	{
		for (const auto& edge : topology.edges)
		{
			if (edge.owner_face_ids != std::vector<std::uint32_t>{1}) continue;
			BRepAdaptor_Curve curve(edge.edge);
			const gp_Pnt a = curve.Value(curve.FirstParameter());
			const gp_Pnt b = curve.Value(curve.LastParameter());
			if (std::abs(a.Y() - 5.0) < 1.0e-9
				&& std::abs(b.Y() - 5.0) < 1.0e-9
				&& std::abs(a.Z() - expected_z) < 1.0e-9
				&& std::abs(b.Z() - expected_z) < 1.0e-9)
				return edge.source_edge_id;
		}
		throw std::runtime_error("could not identify manufactured rib attachment edge");
	}

	void interior_t_fan_three()
	{
		const TopoDS_Face skin = rectangle_xy(0.0, 10.0, 0.0, 10.0);
		const TopoDS_Face rib = rectangle_xz(2.0, 8.0, 5.0, 0.0, 3.0);
		const TopoDS_Compound shape = compound({skin, rib});
		auto topology = make_topology(shape, {skin, rib});
		const std::uint32_t rib_edge = bottom_rib_edge(topology);
		auto atom = std::find_if(topology.atomization.atoms.begin(),
			topology.atomization.atoms.end(), [&](const auto& candidate)
			{
				return std::any_of(candidate.source_spans.begin(),
					candidate.source_spans.end(), [&](const auto& span)
						{ return span.source_edge_id == rib_edge; });
			});
		require(atom != topology.atomization.atoms.end() && atom->fan_degree == 1,
			"manufactured rib boundary atom is missing");
		OcctContactAtomFaceUse skin_use;
		skin_use.source_face_id = 0;
		skin_use.location = OcctFaceIntervalLocation::interior;
		skin_use.sector_count = 2;
		atom->face_uses.push_back(std::move(skin_use));
		atom->fan_degree = 3;

		OcctConformingMesh output;
		std::string error;
		const bool built = build_occt_conforming_mesh(topology, output, error);
		require(built,
			"interior T fan-three build failed: " + error);
		bool skin_interior = false, rib_boundary = false;
		for (const auto& edge : output_face(output, 0).edge_provenance)
			if (edge.contact_atom_id == atom->id && !edge.boundary
				&& !edge.source_cad_edge_id) skin_interior = true;
		for (const auto& edge : output_face(output, 1).edge_provenance)
			if (edge.contact_atom_id == atom->id && edge.boundary
				&& edge.source_cad_edge_id == rib_edge) rib_boundary = true;
		require(skin_interior && rib_boundary,
			"fan-three atom did not retain distinct interior and boundary lineage");
	}

	void exact_operation_tolerance_is_per_face_use()
	{
		constexpr double offset = 5.0e-5;
		constexpr double certified_bound = 6.0e-5;
		const TopoDS_Face skin = rectangle_xy(0.0, 10.0, 0.0, 10.0);
		const TopoDS_Face rib = rectangle_xz(2.0, 8.0, 5.0, offset, 3.0);
		const TopoDS_Compound shape = compound({skin, rib});
		auto topology = make_topology(shape, {skin, rib});
		const std::uint32_t rib_edge = bottom_rib_edge(topology, offset);
		auto atom = std::find_if(topology.atomization.atoms.begin(),
			topology.atomization.atoms.end(), [&](const auto& candidate)
			{
				return std::any_of(candidate.source_spans.begin(),
					candidate.source_spans.end(), [&](const auto& span)
						{ return span.source_edge_id == rib_edge; });
			});
		require(atom != topology.atomization.atoms.end() && atom->fan_degree == 1,
			"offset rib boundary atom is missing");
		OcctContactAtomFaceUse skin_use;
		skin_use.source_face_id = 0;
		skin_use.location = OcctFaceIntervalLocation::interior;
		skin_use.sector_count = 2;
		atom->face_uses.push_back(std::move(skin_use));
		atom->fan_degree = 3;

		OcctConformingMesh output;
		std::string error;
		require(!build_occt_conforming_mesh(topology, output, error)
			&& error.find("interior atom projection failed") != std::string::npos,
			"an underbound exact face use was not rejected: " + error);
		require(output.faces.empty(),
			"underbound face-use rejection published a partial conforming mesh");

		auto& certified_use = atom->face_uses.back();
		certified_use.exact_operation_tolerance = certified_bound;
		require(build_occt_conforming_mesh(topology, output, error),
			"the exact per-face operation certificate was not consumed: " + error);
		bool skin_interior = false;
		for (const auto& edge : output_face(output, 0).edge_provenance)
			if (edge.contact_atom_id == atom->id && !edge.boundary)
				skin_interior = true;
		require(skin_interior,
			"operation-certified interior contact was not inserted into its target face");
	}

	void periodic_seam_and_closed_endpoint()
	{
		const double period = 2.0 * std::acos(-1.0);
		const occ::handle<Geom_CylindricalSurface> surface =
			new Geom_CylindricalSurface(gp_Ax3(gp::Origin(), gp::DZ()), 2.0);
		BRepBuilderAPI_MakeFace maker(surface, 0.0, period, 0.0, 5.0,
			Precision::Confusion());
		require(maker.IsDone(), "could not manufacture periodic cylinder face");
		const TopoDS_Face face = maker.Face();
		const auto topology = make_topology(face, {face});
		const auto seam = std::find_if(topology.atomization.atoms.begin(),
			topology.atomization.atoms.end(), [](const auto& atom)
			{
				return atom.fan_degree == 2 && atom.face_uses.size() == 1
					&& atom.face_uses[0].boundary_occurrences.size() == 2;
			});
		require(seam != topology.atomization.atoms.end(),
			"periodic face did not publish its two seam occurrences");
		const auto closed = std::find_if(topology.atomization.atoms.begin(),
			topology.atomization.atoms.end(), [](const auto& atom)
			{
				return atom.samples.size() >= 3
					&& atom.samples.front().topology_id
						== atom.samples.back().topology_id;
			});
		require(closed != topology.atomization.atoms.end(),
			"periodic fixture did not retain a closed endpoint chain");
		const std::size_t closed_atom_index = static_cast<std::size_t>(
			std::distance(topology.atomization.atoms.begin(), closed));
		OcctConformingMesh output;
		std::string error;
		const bool built = build_occt_conforming_mesh(topology, output, error);
		require(built,
			"periodic seam/closed endpoint build failed: " + error);
		const auto& mesh = output.faces.front().mesh;
		bool retained_wrapped_boundary_node = false;
		for (const auto& edge : output.faces.front().edge_provenance)
			if (edge.source_cad_edge_id == 1u
				&& (edge.vertices[0] == 4u || edge.vertices[1] == 4u))
				retained_wrapped_boundary_node = true;
		require(mesh.source_node_count > 4 && retained_wrapped_boundary_node,
			"monotonic periodic CAD chain did not retain original boundary node 4");
		for (const auto& sample : seam->samples)
		{
			std::size_t aliases = 0;
			for (const auto id : mesh.topology_ids)
				if (id == sample.topology_id) ++aliases;
			require(aliases == 2,
				"periodic seam sample did not retain two explicit UV aliases");
		}
		std::size_t closed_aliases = 0;
		for (const auto id : mesh.topology_ids)
			if (id == closed->samples.front().topology_id) ++closed_aliases;
		require(closed_aliases >= 2,
			"closed periodic boundary endpoint was collapsed in the UV chart");

		// Manufacture the exact case which previously aliased provisional branch
		// identities: two contact atoms cover consecutive pieces of one periodic CAD
		// boundary occurrence.  The untouched cylinder happens to publish the closed
		// occurrence as one atom, so it cannot exercise that identity collision by
		// itself.
		auto split_topology = topology;
		auto& first_atom = split_topology.atomization.atoms[closed_atom_index];
		const OcctContactAtom original_atom = first_atom;
		const std::size_t split_sample = original_atom.samples.size() / 2;
		require(split_sample > 0 && split_sample + 1 < original_atom.samples.size(),
			"closed periodic atom has no exact interior split sample");
		OcctContactAtom second_atom = original_atom;
		second_atom.id = split_topology.atomization.atoms.size();
		first_atom.samples.resize(split_sample + 1);
		second_atom.samples.erase(second_atom.samples.begin(),
			second_atom.samples.begin() + static_cast<std::ptrdiff_t>(split_sample));
		for (std::size_t span_index = 0;
			span_index < original_atom.source_spans.size(); ++span_index)
		{
			const auto& original_span = original_atom.source_spans[span_index];
			const auto source = std::find_if(
				original_atom.samples[split_sample].source_locations.begin(),
				original_atom.samples[split_sample].source_locations.end(),
				[&](const auto& location)
				{
					return location.source_edge_id == original_span.source_edge_id;
				});
			require(source
					!= original_atom.samples[split_sample].source_locations.end(),
				"periodic split sample lacks exact source-edge provenance");
			first_atom.source_spans[span_index].parameter_end = source->parameter;
			second_atom.source_spans[span_index].parameter_begin = source->parameter;
		}
		for (std::size_t face_use_index = 0;
			face_use_index < original_atom.face_uses.size(); ++face_use_index)
		{
			const auto& original_use = original_atom.face_uses[face_use_index];
			for (std::size_t occurrence_index = 0;
				occurrence_index < original_use.boundary_occurrences.size();
				++occurrence_index)
			{
				const auto& original_occurrence =
					original_use.boundary_occurrences[occurrence_index];
				auto& first_occurrence = first_atom.face_uses[face_use_index]
					.boundary_occurrences[occurrence_index];
				auto& second_occurrence = second_atom.face_uses[face_use_index]
					.boundary_occurrences[occurrence_index];
				const double split_target =
					original_occurrence.sample_target_parameters[split_sample];
				const auto mapping_source = std::find_if(
					original_atom.samples[split_sample].source_locations.begin(),
					original_atom.samples[split_sample].source_locations.end(),
					[&](const auto& location)
					{
						return location.source_edge_id
							== original_occurrence.mapping_source_edge_id;
					});
				require(mapping_source
						!= original_atom.samples[split_sample].source_locations.end(),
					"periodic split sample lacks exact mapping-edge provenance");
				first_occurrence.source_parameter_end = mapping_source->parameter;
				first_occurrence.target_parameter_end = split_target;
				first_occurrence.sample_target_parameters.resize(split_sample + 1);
				second_occurrence.source_parameter_begin = mapping_source->parameter;
				second_occurrence.target_parameter_begin = split_target;
				second_occurrence.sample_target_parameters.erase(
					second_occurrence.sample_target_parameters.begin(),
					second_occurrence.sample_target_parameters.begin()
						+ static_cast<std::ptrdiff_t>(split_sample));
			}
		}
		std::vector<OcctContactTopologyMeshNodeClaim> shared_split_claims;
		for (auto& claim : split_topology.mesh_node_claims)
		{
			if (claim.contact_atom_id != original_atom.id) continue;
			if (claim.atom_sample_index > split_sample)
			{
				claim.contact_atom_id = second_atom.id;
				claim.atom_sample_index -= static_cast<std::uint32_t>(split_sample);
			}
			else if (claim.atom_sample_index == split_sample)
			{
				auto duplicate = claim;
				duplicate.contact_atom_id = second_atom.id;
				duplicate.atom_sample_index = 0;
				shared_split_claims.push_back(std::move(duplicate));
			}
		}
		split_topology.mesh_node_claims.insert(split_topology.mesh_node_claims.end(),
			shared_split_claims.begin(), shared_split_claims.end());
		split_topology.atomization.atoms.push_back(std::move(second_atom));

		OcctConformingMesh split_output;
		require(build_occt_conforming_mesh(split_topology, split_output, error),
			"two atoms on one periodic CAD occurrence cross-attached chart branches: "
				+ error);
	}

	TopoDS_Face face_with_hole()
	{
		BRepBuilderAPI_MakePolygon outer_polygon;
		outer_polygon.Add(gp_Pnt(0, 0, 0));
		outer_polygon.Add(gp_Pnt(10, 0, 0));
		outer_polygon.Add(gp_Pnt(10, 10, 0));
		outer_polygon.Add(gp_Pnt(0, 10, 0));
		outer_polygon.Close();
		BRepBuilderAPI_MakePolygon inner_polygon;
		inner_polygon.Add(gp_Pnt(3, 3, 0));
		inner_polygon.Add(gp_Pnt(3, 7, 0));
		inner_polygon.Add(gp_Pnt(7, 7, 0));
		inner_polygon.Add(gp_Pnt(7, 3, 0));
		inner_polygon.Close();
		BRepBuilderAPI_MakeFace maker(outer_polygon.Wire(), true);
		maker.Add(inner_polygon.Wire());
		require(maker.IsDone(), "could not manufacture planar face with hole");
		return maker.Face();
	}

	void hole_preservation()
	{
		const TopoDS_Face face = face_with_hole();
		const auto topology = make_topology(face, {face});
		OcctConformingMesh output;
		std::string error;
		const bool built = build_occt_conforming_mesh(topology, output, error);
		require(built,
			"hole-preserving build failed: " + error);
		require(output.faces.size() == 1 && boundary_loop_count(output.faces.front()) == 2,
			"conforming mesh did not preserve both outer and hole boundary loops");
	}

	OcctConformingMesh sentinel_output()
	{
		OcctConformingMesh result;
		result.shared_topology_node_count = 4321;
		OcctConformingFaceMesh face;
		face.source_face_id = 987;
		face.mesh.uv = {{91, 92}};
		face.mesh.world_positions = {{{93, 94, 95}}};
		face.mesh.triangles = {{{7, 8, 9}}};
		face.mesh.topology_ids = {77};
		face.mesh.edges = {{{7, 8}, true, {1234}}};
		face.edge_provenance = {{{7, 8}, true, 4u, 5u}};
		result.faces.push_back(std::move(face));
		return result;
	}

	void require_sentinel(const OcctConformingMesh& output)
	{
		require(output.shared_topology_node_count == 4321 && output.faces.size() == 1
			&& output.faces[0].source_face_id == 987
			&& output.faces[0].mesh.uv == std::vector<UvPoint>{{91, 92}}
			&& output.faces[0].edge_provenance.size() == 1
			&& output.faces[0].edge_provenance[0].source_cad_edge_id == 4u
			&& output.faces[0].edge_provenance[0].contact_atom_id == 5u,
			"failed build changed the caller's output transaction");
	}

	void rollback_tests()
	{
		const TopoDS_Face face = rectangle_xy(0, 10, 0, 10);
		auto topology = make_topology(face, {face});
		const std::uint32_t damaged = topology.atomization.atoms.front()
			.samples.front().topology_id;
		for (auto& atom : topology.atomization.atoms)
			for (auto& sample : atom.samples)
				if (sample.topology_id == damaged)
					sample.canonical_world_position[2] += 1.0;
		OcctConformingMesh output = sentinel_output();
		std::string error;
		require(!build_occt_conforming_mesh(topology, output, error)
			&& error.find("conforming failed") != std::string::npos,
			"per-face exact-surface rejection did not fail diagnostically: " + error);
		require_sentinel(output);

		const TopoDS_Face first = rectangle_xy(0, 10, 0, 10);
		const TopoDS_Face second = rectangle_xy(20, 30, 0, 10);
		const TopoDS_Compound two = compound({first, second});
		auto two_topology = make_topology(two, {first, second});
		BRepTools::Clean(two_topology.faces[1].face);
		output = sentinel_output();
		require(!build_occt_conforming_mesh(two_topology, output, error)
			&& error.find("source face 1 has no OCCT triangulation") != std::string::npos,
			"global second-face rejection did not identify missing mesh: " + error);
		require_sentinel(output);
	}
}

int main()
{
	try
	{
		tag_namespace_zero_collision();
		std::cout << "[PASS] disjoint CAD0/contact-atom0 tags\n";
		ordinary_shared_boundary_and_fan_one();
		std::cout << "[PASS] ordinary shared boundary and fan-one lineage\n";
		interior_t_fan_three();
		std::cout << "[PASS] interior T fan-three\n";
		exact_operation_tolerance_is_per_face_use();
		std::cout << "[PASS] exact operation tolerance is isolated per face use\n";
		periodic_seam_and_closed_endpoint();
		std::cout << "[PASS] periodic seam aliases and closed endpoint\n";
		hole_preservation();
		std::cout << "[PASS] hole preservation\n";
		rollback_tests();
		std::cout << "[PASS] per-face and global rollback\n";
		std::cout << "All OCCT conforming-mesh builder checks passed.\n";
		return EXIT_SUCCESS;
	}
	catch (const std::exception& exception)
	{
		std::cerr << "[FAIL] " << exception.what() << '\n';
		return EXIT_FAILURE;
	}
}
