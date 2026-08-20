#include "core/geometry/occt_conforming_import_assembly.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRep_Builder.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Precision.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Ax3.hxx>
#include <gp.hxx>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
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
		throw std::runtime_error("manufactured box has no adjacent faces");
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
		OcctContactTopology topology;
		std::string error;
		require(build_occt_contact_topology(shape, faces, {}, topology, error),
			"contact topology failed: " + error);
		return topology;
	}

	OcctConformingImportAssembly assemble(const OcctContactTopology& topology)
	{
		OcctConformingMesh conformed;
		std::string error;
		const bool conformed_ok = build_occt_conforming_mesh(topology, conformed, error);
		require(conformed_ok,
			"conforming mesh failed: " + error);
		OcctConformingImportAssembly output;
		const bool assembled = assemble_occt_conforming_import(topology, conformed,
			output, error);
		require(assembled,
			"import assembly failed: " + error);
		return output;
	}

	bool fp64_bits_equal(const TriMesh& mesh, std::size_t a, std::size_t b)
	{
		for (std::size_t axis = 0; axis < 3; ++axis)
			if (std::bit_cast<std::uint64_t>(mesh.positions_fp64[3 * a + axis])
				!= std::bit_cast<std::uint64_t>(mesh.positions_fp64[3 * b + axis]))
				return false;
		return true;
	}

	void ordinary_and_shared_ids()
	{
		const TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 8.0, 6.0).Shape();
		const auto faces = adjacent_box_faces(box);
		const auto topology = make_topology(box, faces);
		const auto output = assemble(topology);
		const TriMesh& mesh = output.mesh;
		require(output.contacts.audit_status == StepCadContactAuditStatus::not_performed,
			"assembly falsely claimed that its input contact discovery was exhaustive");
		require(!output.contacts.conforming_ready(),
			"assembly without an exhaustive-discovery certificate was marked CFD-ready");
		require(!mesh.empty() && mesh.has_fp64_positions() && mesh.has_uv()
			&& mesh.has_topology_vertex_ids() && mesh.has_face_provenance()
			&& mesh.has_cad_edge_provenance()
			&& mesh.has_cad_edge_contact_provenance()
			&& mesh.has_cad_edge_atom_provenance(),
			"assembled TriMesh hand-off arrays are incomplete");

		std::map<std::uint32_t, std::vector<std::size_t>> copies;
		std::set<std::uint32_t> ordinary;
		for (std::size_t vertex = 0; vertex < mesh.vertex_count(); ++vertex)
		{
			const auto id = mesh.topology_vertex_ids[vertex];
			if (id < topology.atomization.topology_node_count)
				copies[id].push_back(vertex);
			else require(ordinary.insert(id).second,
				"ordinary face-local topology ID was reused");
		}
		for (const auto& [unused_id, vertices] : copies)
		{
			(void)unused_id;
			for (std::size_t copy = 1; copy < vertices.size(); ++copy)
				require(fp64_bits_equal(mesh, vertices.front(), vertices[copy]),
					"shared topology position lost bit identity after SI conversion");
		}
		std::set<std::uint64_t> public_ids;
		for (const auto& curve : output.contacts.curves)
		{
			require(curve.fan_degree >= 2 && curve.id != TriMesh::kNoCadContactId,
				"fan-one atom leaked into public contact graph");
			public_ids.insert(curve.id);
		}
		bool saw_fan_one_atom = false;
		bool saw_public_atom = false;
		for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
			for (unsigned side = 0; side < 3; ++side)
			{
				const auto id = mesh.cad_edge_contact_id(triangle, side);
				if (id == TriMesh::kNoCadContactId)
					require(mesh.cad_edge_certified_fan_degree(triangle, side) == 0,
						"non-contact half-edge carries a fan degree");
				else require(public_ids.contains(id)
					&& mesh.cad_edge_certified_fan_degree(triangle, side) >= 2,
					"mesh contact label is absent from public graph");
				const auto atom_id = mesh.cad_edge_atom_id(triangle, side);
				if (atom_id == TriMesh::kNoCadEdgeAtomId) continue;
				require(atom_id < topology.atomization.atoms.size()
					&& topology.atomization.atoms[static_cast<std::size_t>(atom_id)].id
						== atom_id,
					"mesh exact-atom sidecar references an unknown atom");
				const auto& atom = topology.atomization.atoms[
					static_cast<std::size_t>(atom_id)];
				if (atom.fan_degree == 1)
				{
					saw_fan_one_atom = true;
					require(id == TriMesh::kNoCadContactId
						&& mesh.cad_edge_certified_fan_degree(triangle, side) == 0,
						"fan-one atom was incorrectly published as a CAD contact");
				}
				else
				{
					saw_public_atom = true;
					require(id == atom_id
						&& mesh.cad_edge_certified_fan_degree(triangle, side)
							== atom.fan_degree,
						"fan-two-or-more atom changed its public contact certificate");
				}
			}
		require(saw_fan_one_atom && saw_public_atom,
			"manufactured assembly did not exercise both fan-one and public atoms");
	}

	void periodic_branches()
	{
		const double period = 2.0 * std::acos(-1.0);
		const occ::handle<Geom_CylindricalSurface> surface =
			new Geom_CylindricalSurface(gp_Ax3(gp::Origin(), gp::DZ()), 2.0);
		BRepBuilderAPI_MakeFace maker(surface, 0.0, period, 0.0, 5.0,
			Precision::Confusion());
		require(maker.IsDone(), "could not manufacture periodic cylinder face");
		const TopoDS_Face face = maker.Face();
		auto topology = make_topology(face, {face});
		OcctConformingMesh conformed;
		std::string error;
		require(build_occt_conforming_mesh(topology, conformed, error),
			"periodic conforming mesh failed: " + error);
		OcctConformingImportAssembly output;
		require(assemble_occt_conforming_import(topology, conformed, output, error),
			"periodic assembly failed: " + error);
		auto seam = std::find_if(topology.atomization.atoms.begin(),
			topology.atomization.atoms.end(), [](const auto& atom)
			{
				return atom.fan_degree == 2 && atom.face_uses.size() == 1
					&& atom.face_uses.front().boundary_occurrences.size() == 2;
			});
		require(seam != topology.atomization.atoms.end(),
			"periodic fixture produced no exact two-branch seam atom");
		std::set<std::uint32_t> occurrence_ids;
		for (const auto& branch : output.contact_uv_branches)
			if (branch.contact_atom_id == seam->id)
			{
				require(branch.sample_uv.size() == seam->samples.size(),
					"periodic UV branch lost samples");
				occurrence_ids.insert(branch.boundary_occurrence_id);
			}
		require(occurrence_ids.size() == 2,
			"periodic seam did not preserve two oriented UV branches");
		bool periodic_contact = false;
		for (std::size_t triangle = 0; triangle < output.mesh.triangle_count(); ++triangle)
			for (unsigned side = 0; side < 3; ++side)
				if (output.mesh.cad_edge_contact_id(triangle, side) == seam->id)
				{
					require(output.mesh.cad_edge_atom_id(triangle, side) == seam->id,
						"periodic public contact lost its exact atom identity");
					require(output.mesh.cad_edge_is_periodic_seam(triangle, side),
						"periodic seam contact lost its CAD seam flag");
					periodic_contact = true;
				}
		require(periodic_contact, "periodic seam has no published mesh contact segments");

		// Removing the exact oriented occurrence sidecar forces the interior-contact
		// path to consult the per-face topology-ID index.  A periodic seam has two
		// distinct UV aliases for each shared topology ID, so it must remain rejected
		// rather than silently selecting one chart branch.
		seam->face_uses.front().boundary_occurrences.clear();
		OcctConformingImportAssembly rejected;
		require(!assemble_occt_conforming_import(topology, conformed, rejected, error)
			&& error.find("multiple UV-chart branches") != std::string::npos,
			"indexed UV lookup accepted an unresolved periodic chart alias: " + error);
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

	TopoDS_Face rectangle_xy()
	{
		return rectangle_xy(0.0, 10.0, 0.0, 8.0);
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

	std::uint32_t bottom_rib_edge(const OcctContactTopology& topology)
	{
		for (const auto& edge : topology.edges)
		{
			if (edge.owner_face_ids != std::vector<std::uint32_t>{1}) continue;
			BRepAdaptor_Curve curve(edge.edge);
			const gp_Pnt a = curve.Value(curve.FirstParameter());
			const gp_Pnt b = curve.Value(curve.LastParameter());
			if (std::abs(a.Y() - 5.0) < 1.0e-9
				&& std::abs(b.Y() - 5.0) < 1.0e-9
				&& std::abs(a.Z()) < 1.0e-9 && std::abs(b.Z()) < 1.0e-9)
				return edge.source_edge_id;
		}
		throw std::runtime_error("could not identify manufactured rib attachment edge");
	}

	void reordered_dense_indexes_and_interior_uv()
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
		const std::uint64_t atom_id = atom->id;

		OcctConformingMesh conformed;
		std::string error;
		require(build_occt_conforming_mesh(topology, conformed, error),
			"interior T fan-three conforming failed: " + error);
		std::reverse(topology.faces.begin(), topology.faces.end());
		std::reverse(topology.edges.begin(), topology.edges.end());
		std::reverse(conformed.faces.begin(), conformed.faces.end());

		OcctConformingImportAssembly output;
		require(assemble_occt_conforming_import(topology, conformed, output, error),
			"indexed reordered assembly failed: " + error);
		const auto curve = std::find_if(output.contacts.curves.begin(),
			output.contacts.curves.end(), [&](const auto& candidate)
				{ return candidate.id == atom_id; });
		require(curve != output.contacts.curves.end() && curve->fan_degree == 3,
			"reordered lookup lost the fan-three public contact");
		const auto interior = std::find_if(curve->uses.begin(), curve->uses.end(),
			[](const auto& use)
			{
				return use.source_face_id == 0
					&& use.kind == StepContactUseKind::face_interior;
			});
		require(interior != curve->uses.end()
			&& interior->sample_uv.size()
				== topology.atomization.atoms[static_cast<std::size_t>(atom_id)].samples.size()
			&& std::all_of(interior->sample_uv.begin(), interior->sample_uv.end(),
				[](const auto& uv)
					{ return std::isfinite(uv[0]) && std::isfinite(uv[1]); }),
			"indexed face-local UV chain is incomplete");
		bool saw_public_atom = false;
		for (std::size_t triangle = 0; triangle < output.mesh.triangle_count(); ++triangle)
			for (unsigned side = 0; side < 3; ++side)
				if (output.mesh.cad_edge_atom_id(triangle, side) == atom_id)
				{
					require(output.mesh.cad_edge_contact_id(triangle, side) == atom_id
						&& output.mesh.cad_edge_certified_fan_degree(triangle, side) == 3,
						"indexed assembly changed atom/public fan sidecars");
					saw_public_atom = true;
				}
		require(saw_public_atom, "fan-three atom has no assembled mesh segments");
	}

	void reversed_winding_and_transaction()
	{
		TopoDS_Face face = rectangle_xy();
		face.Reverse();
		BRep_Builder builder;
		TopoDS_Compound shape;
		builder.MakeCompound(shape);
		builder.Add(shape, face);
		auto topology = make_topology(shape, {face});
		OcctConformingMesh conformed;
		std::string error;
		require(build_occt_conforming_mesh(topology, conformed, error),
			"reversed face conforming failed: " + error);
		OcctConformingImportAssembly output;
		require(assemble_occt_conforming_import(topology, conformed, output, error),
			"reversed face assembly failed: " + error);
		for (std::size_t vertex = 0; vertex < output.mesh.vertex_count(); ++vertex)
			require(output.mesh.normals[3 * vertex + 2] < -0.99f,
				"reversed face winding did not produce the oriented negative normal");

		OcctConformingImportAssembly sentinel;
		sentinel.mesh.positions = {91.0f, 92.0f, 93.0f};
		sentinel.contacts.trim_surface_fallback_count = 77;
		sentinel.contact_uv_branches.push_back({123, 4, 5, 0,
			StepContactUseKind::trim_boundary, {{{6, 7}}}});
		conformed.faces.front().mesh.world_positions.front()[0] += 1.0;
		require(!assemble_occt_conforming_import(topology, conformed, sentinel, error)
			&& !error.empty(), "damaged canonical shared position was accepted");
		require(sentinel.mesh.positions == std::vector<float>{91.0f, 92.0f, 93.0f}
			&& sentinel.contacts.trim_surface_fallback_count == 77
			&& sentinel.contact_uv_branches.size() == 1
			&& sentinel.contact_uv_branches.front().contact_atom_id == 123,
			"failed assembly changed the caller transaction");
	}
}

int main()
{
	try
	{
		ordinary_and_shared_ids();
		std::cout << "[PASS] SI TriMesh, ordinary IDs, and public contacts\n";
		periodic_branches();
		std::cout << "[PASS] periodic UV branches and seam provenance\n";
		reordered_dense_indexes_and_interior_uv();
		std::cout << "[PASS] reordered dense indexes and interior UV lookup\n";
		reversed_winding_and_transaction();
		std::cout << "[PASS] reversed winding and transactional rejection\n";
		std::cout << "All conforming import assembly checks passed.\n";
		return EXIT_SUCCESS;
	}
	catch (const std::exception& exception)
	{
		std::cerr << "[FAIL] " << exception.what() << '\n';
		return EXIT_FAILURE;
	}
}
