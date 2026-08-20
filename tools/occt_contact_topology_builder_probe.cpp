// Manufactured checks for deterministic OCCT topology/contact assembly.
#include "core/geometry/occt_contact_topology_builder.h"
#include "core/geometry/occt_trimmed_edge_face.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BezierCurve.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Iterator.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <gp.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace
{
	using namespace paracfd::core;
	using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

	[[noreturn]] void fail(const std::string& message)
	{
		std::fprintf(stderr, "[occt-contact-topology] FAIL: %s\n", message.c_str());
		std::exit(1);
	}

	void require(bool condition, const std::string& message)
	{
		if (!condition) fail(message);
	}

	TopoDS_Face rectangle(double x0, double x1, double y0, double y1)
	{
		BRepBuilderAPI_MakeFace maker(gp_Pln(gp::XOY()), x0, x1, y0, y1);
		require(maker.IsDone(), "could not manufacture rectangle face");
		return maker.Face();
	}

	TopoDS_Face triangle(gp_Pnt a, gp_Pnt b, gp_Pnt c)
	{
		BRepBuilderAPI_MakePolygon polygon;
		polygon.Add(a); polygon.Add(b); polygon.Add(c); polygon.Close();
		require(polygon.IsDone(), "could not manufacture triangle wire");
		BRepBuilderAPI_MakeFace maker(polygon.Wire(), true);
		require(maker.IsDone(), "could not manufacture triangle face");
		return maker.Face();
	}

	TopoDS_Face quadrilateral(gp_Pnt a, gp_Pnt b, gp_Pnt c, gp_Pnt d)
	{
		BRepBuilderAPI_MakePolygon polygon;
		polygon.Add(a); polygon.Add(b); polygon.Add(c); polygon.Add(d); polygon.Close();
		require(polygon.IsDone(), "could not manufacture quadrilateral wire");
		BRepBuilderAPI_MakeFace maker(polygon.Wire(), true);
		require(maker.IsDone(), "could not manufacture quadrilateral face");
		return maker.Face();
	}

	TopoDS_Edge tolerant_linear_edge(const TopoDS_Vertex& first,
		const TopoDS_Vertex& last, const gp_Pnt& curve_first,
		const gp_Pnt& curve_last)
	{
		TColgp_Array1OfPnt poles(1, 2);
		poles.SetValue(1, curve_first);
		poles.SetValue(2, curve_last);
		const Handle(Geom_BezierCurve) curve = new Geom_BezierCurve(poles);
		BRepBuilderAPI_MakeEdge maker(curve, first, last, 0.0, 1.0);
		require(maker.IsDone(), "could not manufacture tolerant CAD edge");
		return maker.Edge();
	}

	TopoDS_Shape compound(const std::vector<TopoDS_Face>& faces)
	{
		TopoDS_Compound shape;
		BRep_Builder builder;
		builder.MakeCompound(shape);
		for (const auto& face : faces) builder.Add(shape, face);
		return shape;
	}

	void mesh(const TopoDS_Shape& shape, double deflection = 0.05)
	{
		BRepMesh_IncrementalMesh mesher(shape, deflection, false, 0.1, false);
		mesher.Perform();
		require(mesher.IsDone(), "manufactured shape did not mesh");
	}

	std::array<gp_Pnt, 2> endpoints(const TopoDS_Edge& edge)
	{
		BRepAdaptor_Curve curve(edge);
		return {{curve.Value(curve.FirstParameter()), curve.Value(curve.LastParameter())}};
	}

	bool same_point(const gp_Pnt& a, const gp_Pnt& b, double tolerance = 1.0e-8)
	{
		return a.SquareDistance(b) <= tolerance * tolerance;
	}

	std::uint32_t edge_between(const TopoDS_Shape& shape, gp_Pnt a, gp_Pnt b)
	{
		ShapeMap edges;
		TopExp::MapShapes(shape, TopAbs_EDGE, edges);
		for (Standard_Integer index = 1; index <= edges.Extent(); ++index)
		{
			const auto edge = TopoDS::Edge(edges(index));
			if (BRep_Tool::Degenerated(edge)) continue;
			const auto ends = endpoints(edge);
			if ((same_point(ends[0], a) && same_point(ends[1], b))
				|| (same_point(ends[0], b) && same_point(ends[1], a)))
				return static_cast<std::uint32_t>(index - 1);
		}
		fail("could not find manufactured edge by exact endpoints");
	}

	bool face_owns_edge(const TopoDS_Face& face, const TopoDS_Edge& edge)
	{
		for (TopExp_Explorer occurrence(face, TopAbs_EDGE); occurrence.More();
			occurrence.Next())
			if (occurrence.Current().IsSame(edge)) return true;
		return false;
	}

	std::vector<OcctContactTopologyEdgeFaceInterval> certify_pairs(
		const TopoDS_Shape& shape, const std::vector<TopoDS_Face>& faces,
		const std::vector<std::pair<std::uint32_t, std::uint32_t>>& pairs)
	{
		ShapeMap edges;
		TopExp::MapShapes(shape, TopAbs_EDGE, edges);
		std::vector<OcctContactTopologyEdgeFaceInterval> result;
		for (const auto [edge_id, face_id] : pairs)
		{
			require(edge_id < static_cast<std::uint32_t>(edges.Extent())
				&& face_id < faces.size(), "certification pair is out of range");
			const TopoDS_Edge edge = TopoDS::Edge(edges(static_cast<int>(edge_id + 1)));
			require(!BRep_Tool::Degenerated(edge), "cannot certify a degenerate test edge");
			require(!face_owns_edge(faces[face_id], edge),
				"test certification pair incorrectly names an owner face");
			const auto common = exact_trimmed_edge_face_common(edge, faces[face_id]);
			if (!common.valid())
			{
				std::string detail;
				for (const auto& error : common.errors) detail += " " + error;
				fail("exact manufactured edge/face certification failed:" + detail);
			}
			require(!common.intervals.empty(), "certified test pair has no positive interval");
			for (const auto& interval : common.intervals)
				result.push_back({edge_id, face_id,
					common.source_parameter_first, common.source_parameter_last, interval});
		}
		return result;
	}

	OcctContactTopology build(const TopoDS_Shape& shape,
		const std::vector<TopoDS_Face>& faces,
		const std::vector<OcctContactTopologyEdgeFaceInterval>& intervals)
	{
		OcctContactTopology result;
		std::string error;
		require(build_occt_contact_topology(shape, faces, intervals, result, error), error);
		return result;
	}

	const OcctContactAtom* atom_with_edge(const OcctContactTopology& topology,
		std::uint32_t edge_id, double midpoint_x)
	{
		for (const auto& atom : topology.atomization.atoms)
		{
			const auto span = std::find_if(atom.source_spans.begin(), atom.source_spans.end(),
				[&](const auto& candidate) { return candidate.source_edge_id == edge_id; });
			if (span == atom.source_spans.end()) continue;
			const auto source = std::find_if(topology.atomizer_source_edges.begin(),
				topology.atomizer_source_edges.end(), [&](const auto& candidate)
				{ return candidate.source_edge_id == edge_id; });
			require(source != topology.atomizer_source_edges.end(), "atom lost source edge");
			BRepAdaptor_Curve curve(source->edge);
			const gp_Pnt middle = curve.Value(0.5
				* (span->parameter_begin + span->parameter_end));
			if (std::abs(middle.X() - midpoint_x) < 1.0e-7) return &atom;
		}
		return nullptr;
	}

	void install_target_midpoint_triangulation(const TopoDS_Face& face)
	{
		// Planar rectangle [0,2]x[-1,1], with one face-interior node exactly on
		// the manufactured contact line y=0.  It is deliberately not an edge-
		// polygon node of the independent source face.
		const Handle(Poly_Triangulation) triangulation =
			new Poly_Triangulation(5, 4, true, false);
		const std::array<gp_Pnt, 5> nodes{{
			{0.0, -1.0, 0.0}, {2.0, -1.0, 0.0}, {2.0, 1.0, 0.0},
			{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}}};
		const std::array<gp_Pnt2d, 5> uv{{
			{0.0, -1.0}, {2.0, -1.0}, {2.0, 1.0},
			{0.0, 1.0}, {1.0, 0.0}}};
		for (Standard_Integer node = 1; node <= 5; ++node)
		{
			triangulation->SetNode(node, nodes[static_cast<std::size_t>(node - 1)]);
			triangulation->SetUVNode(node, uv[static_cast<std::size_t>(node - 1)]);
		}
		triangulation->SetTriangle(1, Poly_Triangle(1, 2, 5));
		triangulation->SetTriangle(2, Poly_Triangle(2, 3, 5));
		triangulation->SetTriangle(3, Poly_Triangle(3, 4, 5));
		triangulation->SetTriangle(4, Poly_Triangle(4, 1, 5));
		triangulation->Deflection(0.05);
		BRep_Builder builder;
		builder.UpdateFace(face, triangulation, true);
	}

	std::string fingerprint(const OcctContactTopology& topology)
	{
		std::ostringstream out;
		out << std::setprecision(std::numeric_limits<double>::max_digits10);
		out << topology.edges.size() << ' ' << topology.faces.size() << ' '
			<< topology.reciprocal_intervals.size() << ' '
			<< topology.exact_junctions.size() << ' '
			<< topology.atomization.topology_node_count << ' '
			<< topology.mesh_node_claims.size() << '\n';
		for (const auto& face : topology.faces)
		{
			out << 'f' << face.source_face_id << ':';
			for (const auto& occurrence : face.boundary_occurrences)
				out << occurrence.occurrence_id << ',' << occurrence.source_edge_id << ','
					<< occurrence.wire_id << ',' << static_cast<int>(occurrence.orientation) << ';';
			out << '\n';
		}
		for (const auto& reciprocal : topology.reciprocal_intervals)
			out << 'r' << reciprocal.source_edge_a << ',' << reciprocal.source_edge_b << ','
				<< reciprocal.parameter_a_begin << ',' << reciprocal.parameter_a_end << ','
				<< reciprocal.parameter_b_begin << ',' << reciprocal.parameter_b_end << '\n';
		for (const auto& junction : topology.exact_junctions)
		{
			out << 'j' << junction.junction_id << ':';
			for (const auto& incidence : junction.incidences)
				out << incidence.source_edge_id << ',' << incidence.parameter << ';';
			out << junction.canonical_world_position[0] << ','
				<< junction.canonical_world_position[1] << ','
				<< junction.canonical_world_position[2] << '\n';
		}
		for (const auto& atom : topology.atomization.atoms)
		{
			out << 'a' << atom.id << ',' << atom.fan_degree << ':';
			for (const auto& span : atom.source_spans)
				out << span.source_edge_id << ',' << span.parameter_begin << ','
					<< span.parameter_end << ';';
			out << '|';
			for (const auto& sample : atom.samples)
				out << sample.canonical_world_position[0] << ','
					<< sample.canonical_world_position[1] << ','
					<< sample.canonical_world_position[2] << ',' << sample.topology_id << ';';
			out << '\n';
		}
		for (const auto& claim : topology.mesh_node_claims)
			out << 'm' << claim.source_face_id << ',' << claim.source_node_index << ','
				<< claim.boundary_occurrence_id << ',' << claim.contact_atom_id << ','
				<< claim.atom_sample_index << ',' << claim.topology_id << '\n';
		return out.str();
	}

	void fan_one_opening()
	{
		std::vector<TopoDS_Face> faces{rectangle(0.0, 2.0, 0.0, 1.0)};
		const TopoDS_Shape shape = compound(faces);
		mesh(shape);
		const auto topology = build(shape, faces, {});
		require(topology.edges.size() == 4, "rectangle global edge index is incomplete");
		ShapeMap indexed_edges;
		TopExp::MapShapes(shape, TopAbs_EDGE, indexed_edges);
		for (std::size_t edge = 0; edge < topology.edges.size(); ++edge)
			require(topology.edges[edge].source_edge_id == edge
				&& topology.edges[edge].edge.IsSame(indexed_edges(static_cast<int>(edge + 1))),
				"published source-edge ID differs from TopExp::MapShapes traversal");
		for (std::size_t occurrence = 0;
			occurrence < topology.faces[0].boundary_occurrences.size(); ++occurrence)
			require(topology.faces[0].boundary_occurrences[occurrence].occurrence_id
					== occurrence,
				"face occurrence ID is not its exact oriented lookup index");
		require(topology.atomization.atoms.size() == 4,
			"rectangle should produce four unsplit boundary atoms");
		for (const auto& atom : topology.atomization.atoms)
			require(atom.is_open_boundary(), "fan-one opening was incorrectly healed");
	}

	struct SplitFixture
	{
		TopoDS_Shape shape;
		std::vector<TopoDS_Face> faces;
		std::uint32_t long_edge = 0;
		std::uint32_t short_a = 0;
		std::uint32_t short_b = 0;
		std::vector<OcctContactTopologyEdgeFaceInterval> intervals;
	};

	SplitFixture split_fixture()
	{
		SplitFixture fixture;
		fixture.faces = {rectangle(0.0, 2.0, 0.0, 1.0),
			triangle({0.0, 0.0, 0.0}, {0.4, -1.0, 0.0}, {1.0, 0.0, 0.0}),
			triangle({1.0, 0.0, 0.0}, {1.6, -1.0, 0.0}, {2.0, 0.0, 0.0})};
		fixture.shape = compound(fixture.faces);
		mesh(fixture.shape);
		fixture.long_edge = edge_between(fixture.shape, {0, 0, 0}, {2, 0, 0});
		fixture.short_a = edge_between(fixture.shape, {0, 0, 0}, {1, 0, 0});
		fixture.short_b = edge_between(fixture.shape, {1, 0, 0}, {2, 0, 0});
		fixture.intervals = certify_pairs(fixture.shape, fixture.faces,
			{{fixture.long_edge, 1u}, {fixture.long_edge, 2u},
				{fixture.short_a, 0u}, {fixture.short_b, 0u}});
		return fixture;
	}

	void long_against_two_short_and_asymmetric_samples()
	{
		const SplitFixture fixture = split_fixture();
		const auto topology = build(fixture.shape, fixture.faces, fixture.intervals);
		const OcctContactAtom* left = atom_with_edge(topology, fixture.long_edge, 0.5);
		const OcctContactAtom* right = atom_with_edge(topology, fixture.long_edge, 1.5);
		require(left && right && left != right,
			"long edge was not atomized against its two exact short contacts");
		require(left->fan_degree == 2 && right->fan_degree == 2,
			"split contact atoms did not retain two fluid-sheet sectors");
		auto has_source = [](const OcctContactAtom& atom, std::uint32_t edge)
		{
			return std::any_of(atom.source_spans.begin(), atom.source_spans.end(),
				[&](const auto& span) { return span.source_edge_id == edge; });
		};
		require(has_source(*left, fixture.short_a) && has_source(*right, fixture.short_b),
			"reciprocal short-edge provenance was not retained");
		const auto certified_use = std::find_if(fixture.intervals.begin(),
			fixture.intervals.end(), [&](const auto& interval)
			{
				return interval.source_edge_id == fixture.long_edge
					&& interval.target_face_id == 1u;
			});
		require(certified_use != fixture.intervals.end()
			&& certified_use->interval.exact_operation_tolerance > 0.0,
			"manufactured Common did not retain its result-edge tolerance");
		const auto target_face_use = std::find_if(left->face_uses.begin(),
			left->face_uses.end(), [](const auto& use)
				{ return use.source_face_id == 1u; });
		require(target_face_use != left->face_uses.end()
			&& target_face_use->exact_operation_tolerance
				>= certified_use->interval.exact_operation_tolerance,
			"topology builder did not propagate the interval's local operation bound "
			"to its target face use");
		require(left->samples.size() >= 2,
			"contact atom lost its exact endpoint samples");
		for (const auto& sample : left->samples)
		{
			bool has_long = false, has_short = false;
			for (const auto& source : sample.source_locations)
			{
				has_long = has_long || source.source_edge_id == fixture.long_edge;
				has_short = has_short || source.source_edge_id == fixture.short_a;
			}
			require(has_long && has_short,
				"union-chain sample was not recovered on every reciprocal edge copy");
		}

		const auto short_occurrence = std::find_if(
			topology.faces[1].boundary_occurrences.begin(),
			topology.faces[1].boundary_occurrences.end(), [&](const auto& occurrence)
			{ return occurrence.source_edge_id == fixture.short_a; });
		require(short_occurrence != topology.faces[1].boundary_occurrences.end(),
			"target face lost its exact short-edge boundary occurrence");
		const auto boundary_claim = std::find_if(topology.mesh_node_claims.begin(),
			topology.mesh_node_claims.end(), [&](const auto& claim)
			{
				return claim.source_face_id == 1u
					&& claim.boundary_occurrence_id == short_occurrence->occurrence_id
					&& claim.contact_atom_id == left->id
					&& claim.atom_sample_index < left->samples.size()
					&& left->samples[claim.atom_sample_index].topology_id
						== claim.topology_id;
			});
		require(boundary_claim != topology.mesh_node_claims.end(),
			"target boundary tessellation node did not retain exact occurrence/node provenance");
	}

	void target_face_mid_contact_sample_is_global()
	{
		// The source face is vertical, so only its lower boundary lies in the
		// horizontal target face.  The target's center tessellation vertex therefore
		// lies on an already-certified *interior* edge/face contact, without creating
		// a second CAD edge or a proximity-discovered relationship.
		std::vector<TopoDS_Face> faces{
			quadrilateral({0.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
				{2.0, 0.0, 1.0}, {0.0, 0.0, 1.0}),
			rectangle(0.0, 2.0, -1.0, 1.0)};
		const TopoDS_Shape shape = compound(faces);
		mesh(shape);
		const std::uint32_t source_edge = edge_between(shape,
			{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0});
		install_target_midpoint_triangulation(faces[1]);
		auto intervals = certify_pairs(shape, faces, {{source_edge, 1u}});
		// A deliberately loose operation certificate is local to the target face.
		// It may select a target node, but must never become atom-wide canonical
		// uncertainty inherited by the independent owner face.
		for (auto& interval : intervals)
			interval.interval.exact_operation_tolerance = 0.25;
		const auto topology = build(shape, faces, intervals);

		const auto source = std::find_if(topology.atomizer_source_edges.begin(),
			topology.atomizer_source_edges.end(), [&](const auto& candidate)
			{ return candidate.source_edge_id == source_edge; });
		require(source != topology.atomizer_source_edges.end(),
			"mid-contact fixture lost its source edge");
		BRepAdaptor_Curve curve(source->edge);
		bool retained_midpoint = false;
		for (const auto& sample : source->samples)
			retained_midpoint = retained_midpoint
				|| same_point(curve.Value(sample.parameter), {1.0, 0.0, 0.0});
		require(retained_midpoint,
			"target-face interior tessellation node was not added to the global source chain");

		const OcctContactAtom* atom = atom_with_edge(topology, source_edge, 1.0);
		require(atom && atom->fan_degree == 3 && atom->face_uses.size() == 2,
			"interior contact did not retain one owner sector plus two target sectors");
		const auto midpoint = std::find_if(atom->samples.begin(), atom->samples.end(),
			[](const auto& sample)
			{
				return same_point({sample.canonical_world_position[0],
					sample.canonical_world_position[1],
					sample.canonical_world_position[2]}, {1.0, 0.0, 0.0});
			});
		require(midpoint != atom->samples.end()
			&& midpoint->topology_id != std::numeric_limits<std::uint32_t>::max(),
			"target-derived midpoint was not atomized as one global topology sample");
		const std::uint32_t midpoint_index = static_cast<std::uint32_t>(
			std::distance(atom->samples.begin(), midpoint));
		const auto midpoint_claim = std::find_if(topology.mesh_node_claims.begin(),
			topology.mesh_node_claims.end(), [&](const auto& claim)
			{
				return claim.source_face_id == 1u && claim.source_node_index == 4u
					&& claim.boundary_occurrence_id
						== OcctContactBoundaryOccurrence::no_occurrence_id
					&& claim.contact_atom_id == atom->id
					&& claim.atom_sample_index == midpoint_index
					&& claim.topology_id == midpoint->topology_id;
			});
		require(midpoint_claim != topology.mesh_node_claims.end(),
			"target interior tessellation midpoint lost its exact original node provenance");
		require(midpoint->canonical_position_tolerance < 0.01,
			"target-face operation tolerance leaked into the atom-wide sample bound");
		const auto target_use = std::find_if(atom->face_uses.begin(), atom->face_uses.end(),
			[](const auto& use) { return use.source_face_id == 1u; });
		require(target_use != atom->face_uses.end()
			&& target_use->exact_operation_tolerance >= 0.25,
			"target-face operation tolerance was not retained on its isolated face use");
		for (const auto& use : atom->face_uses)
			for (const auto& occurrence : use.boundary_occurrences)
				require(occurrence.sample_target_parameters.size() == atom->samples.size(),
					"one fan use did not inherit the global contact segmentation");
	}

	void t_junction()
	{
		const TopoDS_Edge skin_contact = BRepBuilderAPI_MakeEdge(
			gp_Lin({100.0, 0.0, 0.0}, {1.0, 0.0, 0.0}), -100.0, -98.0).Edge();
		const TopoDS_Edge rib_contact = BRepBuilderAPI_MakeEdge(
			gp_Lin({1.0, 1000.0, 0.0}, {0.0, -1.0, 0.0}), 1000.0, 1001.0).Edge();
		TopoDS_Vertex skin_first, skin_last, rib_first, rib_last;
		TopExp::Vertices(skin_contact, skin_first, skin_last, true);
		TopExp::Vertices(rib_contact, rib_first, rib_last, true);
		const TopoDS_Vertex upper_right = BRepBuilderAPI_MakeVertex({2.0, 1.0, 0.0});
		const TopoDS_Vertex upper_left = BRepBuilderAPI_MakeVertex({0.0, 1.0, 0.0});
		BRepBuilderAPI_MakeWire skin_wire;
		skin_wire.Add(skin_contact);
		skin_wire.Add(BRepBuilderAPI_MakeEdge(skin_last, upper_right).Edge());
		skin_wire.Add(BRepBuilderAPI_MakeEdge(upper_right, upper_left).Edge());
		skin_wire.Add(BRepBuilderAPI_MakeEdge(upper_left, skin_first).Edge());
		const TopoDS_Vertex rib_side = BRepBuilderAPI_MakeVertex({1.2, -1.0, 0.0});
		BRepBuilderAPI_MakeWire rib_wire;
		rib_wire.Add(rib_contact);
		rib_wire.Add(BRepBuilderAPI_MakeEdge(rib_last, rib_side).Edge());
		rib_wire.Add(BRepBuilderAPI_MakeEdge(rib_side, rib_first).Edge());
		require(skin_wire.IsDone() && rib_wire.IsDone(),
			"could not manufacture nonmatching-parameter T wires");
		BRepBuilderAPI_MakeFace skin_face(skin_wire.Wire(), true);
		BRepBuilderAPI_MakeFace rib_face(rib_wire.Wire(), true);
		require(skin_face.IsDone() && rib_face.IsDone(),
			"could not manufacture nonmatching-parameter T faces");
		std::vector<TopoDS_Face> faces{skin_face.Face(), rib_face.Face()};
		const TopoDS_Shape shape = compound(faces);
		mesh(shape);
		const auto topology = build(shape, faces, {});
		const std::uint32_t long_edge = edge_between(shape, {0, 0, 0}, {2, 0, 0});
		const std::uint32_t rib_edge = edge_between(shape, {1, 0, 0}, {1, -1, 0});
		std::set<std::uint32_t> topology_ids;
		std::set<std::uint32_t> source_edges;
		std::optional<double> long_parameter, rib_parameter;
		for (const auto& atom : topology.atomization.atoms)
			for (const auto& sample : atom.samples)
				if (std::abs(sample.canonical_world_position[0] - 1.0) < 1.0e-8
					&& std::abs(sample.canonical_world_position[1]) < 1.0e-8)
				{
					topology_ids.insert(sample.topology_id);
					for (const auto& source : sample.source_locations)
					{
						source_edges.insert(source.source_edge_id);
						if (source.source_edge_id == long_edge)
							long_parameter = source.parameter;
						if (source.source_edge_id == rib_edge)
							rib_parameter = source.parameter;
					}
				}
		require(topology_ids.size() == 1 && source_edges.contains(long_edge)
			&& source_edges.size() >= 3 && long_parameter && rib_parameter
			&& std::abs(*long_parameter + 99.0) <= 1.0e-10
			&& std::abs(*rib_parameter - 1000.0) <= 1.0e-10,
			"zero-fuzzy T junction did not create one exact multi-edge topology node");
	}

	void shared_topology_vertex_suppresses_analytic_duplicates()
	{
		TopoDS_Vertex origin = BRepBuilderAPI_MakeVertex({0.0, 0.0, 0.0});
		BRep_Builder tolerance_builder;
		tolerance_builder.UpdateVertex(origin, 5.0e-5);
		auto vertex = [](double x, double y, double z)
		{
			return BRepBuilderAPI_MakeVertex({x, y, z}).Vertex();
		};
		constexpr double offset = 1.0e-5;
		const TopoDS_Vertex upper_end = vertex(1.0, offset - 1.0, 0.0);
		const TopoDS_Vertex upper_far = vertex(-1.0, -1.0, 0.0);
		const TopoDS_Vertex lower_end = vertex(1.0, 1.0 - offset, 0.0);
		const TopoDS_Vertex lower_far = vertex(-1.0, 1.0, 0.0);
		const TopoDS_Edge upper_contact = tolerant_linear_edge(origin, upper_end,
			{0.0, offset, 0.0}, {1.0, offset - 1.0, 0.0});
		const TopoDS_Edge lower_contact = tolerant_linear_edge(origin, lower_end,
			{0.0, -offset, 0.0}, {1.0, 1.0 - offset, 0.0});
		BRepBuilderAPI_MakeWire upper_wire;
		upper_wire.Add(upper_contact);
		upper_wire.Add(BRepBuilderAPI_MakeEdge(upper_end, upper_far).Edge());
		upper_wire.Add(BRepBuilderAPI_MakeEdge(upper_far, origin).Edge());
		BRepBuilderAPI_MakeWire lower_wire;
		lower_wire.Add(lower_contact);
		lower_wire.Add(BRepBuilderAPI_MakeEdge(lower_end, lower_far).Edge());
		lower_wire.Add(BRepBuilderAPI_MakeEdge(lower_far, origin).Edge());
		require(upper_wire.IsDone() && lower_wire.IsDone(),
			"could not manufacture shared-edge triangle wires");
		BRepBuilderAPI_MakeFace upper_face(upper_wire.Wire(), true);
		BRepBuilderAPI_MakeFace lower_face(lower_wire.Wire(), true);
		require(upper_face.IsDone() && lower_face.IsDone(),
			"could not manufacture shared-edge triangle faces");
		std::vector<TopoDS_Face> faces{upper_face.Face(), lower_face.Face()};
		const TopoDS_Shape shape = compound(faces);
		mesh(shape);
		const auto topology = build(shape, faces, {});

		std::size_t origin_junctions = 0;
		std::size_t origin_incidence_count = 0;
		for (const auto& junction : topology.exact_junctions)
			if (same_point(gp_Pnt(junction.canonical_world_position[0],
				junction.canonical_world_position[1],
				junction.canonical_world_position[2]), {0.0, 0.0, 0.0}))
			{
				++origin_junctions;
				origin_incidence_count = std::max(origin_incidence_count,
					junction.incidences.size());
			}
		require(origin_junctions == 1 && origin_incidence_count == 4,
			"analytic edge intersections duplicated one exact multi-edge topology vertex");
	}

	void seam_and_closed_edge()
	{
		const TopoDS_Shape shape = BRepPrimAPI_MakeCylinder(2.0, 5.0).Shape();
		mesh(shape, 0.15);
		TopoDS_Face lateral;
		for (TopExp_Explorer face(shape, TopAbs_FACE); face.More() && lateral.IsNull();
			face.Next())
		{
			std::vector<TopoDS_Edge> occurrences;
			for (TopExp_Explorer wires(face.Current(), TopAbs_WIRE); wires.More(); wires.Next())
				for (TopoDS_Iterator child(wires.Current(), true, true); child.More(); child.Next())
					if (child.Value().ShapeType() == TopAbs_EDGE)
						occurrences.push_back(TopoDS::Edge(child.Value()));
			for (std::size_t a = 0; a < occurrences.size(); ++a)
				for (std::size_t b = a + 1; b < occurrences.size(); ++b)
					if (occurrences[a].IsSame(occurrences[b]))
						lateral = TopoDS::Face(face.Current());
		}
		require(!lateral.IsNull(), "could not find cylindrical periodic face");
		std::vector<TopoDS_Face> faces{lateral};
		const auto topology = build(shape, faces, {});
		std::uint32_t seam_id = std::numeric_limits<std::uint32_t>::max();
		std::map<std::uint32_t, int> occurrence_counts;
		for (const auto& occurrence : topology.faces[0].boundary_occurrences)
			++occurrence_counts[occurrence.source_edge_id];
		for (const auto& [edge_id, occurrences] : occurrence_counts)
			if (occurrences == 2) seam_id = edge_id;
		require(seam_id != std::numeric_limits<std::uint32_t>::max(),
			"oriented occurrence table lost periodic seam branches");
		TopLoc_Location triangulation_location;
		const Handle(Poly_Triangulation) triangulation =
			BRep_Tool::Triangulation(faces[0], triangulation_location);
		require(!triangulation.IsNull(), "periodic face lost its source triangulation");
		std::set<std::uint32_t> boundary_polygon_nodes;
		for (const auto& occurrence : topology.faces[0].boundary_occurrences)
		{
			TopoDS_Edge reversed = occurrence.edge;
			reversed.Reverse();
			const std::array<Handle(Poly_PolygonOnTriangulation), 2> polygons{{
				BRep_Tool::PolygonOnTriangulation(occurrence.edge,
					triangulation, triangulation_location),
				BRep_Tool::PolygonOnTriangulation(reversed,
					triangulation, triangulation_location)}};
			require(!polygons[0].IsNull() || !polygons[1].IsNull(),
				"periodic boundary occurrence has no source edge polygon");
			for (const auto& polygon : polygons)
				if (!polygon.IsNull())
					for (Standard_Integer index = 1; index <= polygon->NbNodes(); ++index)
					{
						const Standard_Integer node = polygon->Node(index);
						require(node > 0 && node <= triangulation->NbNodes(),
							"periodic edge polygon references an invalid source node");
						boundary_polygon_nodes.insert(static_cast<std::uint32_t>(node - 1));
					}
		}
		for (const std::uint32_t node : boundary_polygon_nodes)
			require(std::any_of(topology.mesh_node_claims.begin(),
				topology.mesh_node_claims.end(), [&](const auto& claim)
				{
					return claim.source_face_id == 0u && claim.source_node_index == node
						&& claim.boundary_occurrence_id
							!= OcctContactBoundaryOccurrence::no_occurrence_id;
				}), "periodic boundary polygon node lost its explicit original-node claim");
		const auto seam_atom = std::find_if(topology.atomization.atoms.begin(),
			topology.atomization.atoms.end(), [&](const auto& atom)
			{
				return std::any_of(atom.source_spans.begin(), atom.source_spans.end(),
					[&](const auto& span) { return span.source_edge_id == seam_id; });
			});
		require(seam_atom != topology.atomization.atoms.end()
			&& seam_atom->fan_degree == 2
			&& seam_atom->face_uses.size() == 1
			&& seam_atom->face_uses[0].boundary_occurrences.size() == 2,
			"periodic seam did not retain two oriented branches/two sectors");
		for (const auto& occurrence : seam_atom->face_uses[0].boundary_occurrences)
			require(std::any_of(topology.mesh_node_claims.begin(),
				topology.mesh_node_claims.end(), [&](const auto& claim)
				{
					return claim.source_face_id == 0u
						&& claim.boundary_occurrence_id == occurrence.occurrence_id
						&& claim.contact_atom_id == seam_atom->id;
				}), "periodic seam occurrence lost its independent mesh-node claims");

		bool found_closed_identity = false;
		for (const auto& atom : topology.atomization.atoms)
			if (atom.fan_degree == 1 && atom.samples.size() > 2
				&& atom.samples.front().topology_id == atom.samples.back().topology_id)
			{
				const auto first_claim = std::find_if(topology.mesh_node_claims.begin(),
					topology.mesh_node_claims.end(), [&](const auto& claim)
					{
						return claim.contact_atom_id == atom.id
							&& claim.atom_sample_index == 0u;
					});
				const auto last_claim = std::find_if(topology.mesh_node_claims.begin(),
					topology.mesh_node_claims.end(), [&](const auto& claim)
					{
						return claim.contact_atom_id == atom.id
							&& claim.atom_sample_index == atom.samples.size() - 1;
					});
				found_closed_identity = first_claim != topology.mesh_node_claims.end()
					&& last_claim != topology.mesh_node_claims.end()
					&& first_claim->source_node_index != last_claim->source_node_index;
			}
		require(found_closed_identity,
			"closed circular edge first/last samples lost distinct original UV branches");
	}

	void deterministic_reorder_and_transaction()
	{
		const SplitFixture fixture = split_fixture();
		const auto forward = build(fixture.shape, fixture.faces, fixture.intervals);
		auto reversed_intervals = fixture.intervals;
		std::reverse(reversed_intervals.begin(), reversed_intervals.end());
		const auto reversed = build(fixture.shape, fixture.faces, reversed_intervals);
		require(fingerprint(forward) == fingerprint(reversed),
			"reordering exact interval records changed topology/atomization output");

		OcctContactTopology preserved;
		preserved.edges.push_back({12345u, {}, {678u}});
		const std::string before = fingerprint(preserved);
		auto malformed = fixture.intervals;
		require(!malformed.empty(), "split fixture produced no exact intervals");
		malformed.front().interval.location = OcctFaceIntervalLocation::boundary;
		malformed.front().interval.target_sector_count = 1;
		malformed.front().interval.boundary_occurrence_count = 1;
		malformed.front().interval.target_boundary_occurrence_ids[0] =
			std::numeric_limits<std::uint32_t>::max();
		std::string error;
		require(!build_occt_contact_topology(fixture.shape, fixture.faces, malformed,
			preserved, error) && !error.empty(),
			"malformed occurrence mapping did not fail closed");
		require(fingerprint(preserved) == before,
			"failed topology build mutated its caller-owned output transaction");
	}
}

int main()
{
	fan_one_opening();
	std::puts("[occt-contact-topology] fan-one/open boundary: PASS");
	long_against_two_short_and_asymmetric_samples();
	std::puts("[occt-contact-topology] long-vs-two-short + asymmetric samples: PASS");
	target_face_mid_contact_sample_is_global();
	std::puts("[occt-contact-topology] target-face midpoint global sample: PASS");
	t_junction();
	std::puts("[occt-contact-topology] zero-fuzzy T junction: PASS");
	shared_topology_vertex_suppresses_analytic_duplicates();
	std::puts("[occt-contact-topology] shared topology vertex duplicate suppression: PASS");
	seam_and_closed_edge();
	std::puts("[occt-contact-topology] periodic seam + closed endpoint identity: PASS");
	deterministic_reorder_and_transaction();
	std::puts("[occt-contact-topology] deterministic reorder + transaction: PASS");
	return 0;
}
