// Deterministic manufactured checks for exact CAD contact interval atomization.
#include "core/geometry/occt_contact_atomizer.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRep_Builder.hxx>
#include <GeomConvert.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Circle.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Pnt.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace
{
	using namespace paracfd::core;

	[[noreturn]] void fail(const std::string& message)
	{
		std::fprintf(stderr, "[occt-contact-atomizer] FAIL: %s\n", message.c_str());
		std::exit(1);
	}

	void require(bool condition, const std::string& message)
	{
		if (!condition) fail(message);
	}

	bool near(double a, double b, double tolerance = 1.0e-9)
	{
		return std::abs(a - b) <= tolerance;
	}

	std::array<double, 2> range_of(const TopoDS_Edge& edge)
	{
		BRepAdaptor_Curve curve(edge);
		return {{curve.FirstParameter(), curve.LastParameter()}};
	}

	TopoDS_Edge line(gp_Pnt a, gp_Pnt b)
	{
		return BRepBuilderAPI_MakeEdge(a, b).Edge();
	}

	TopoDS_Edge nonlinear_line()
	{
		// A geometrically straight 0..10 segment whose quadratic Bezier parameter is
		// deliberately nonlinear: x(t) = 2t + 8t^2.
		TColgp_Array1OfPnt poles(1, 3);
		poles.SetValue(1, gp_Pnt(0.0, 0.0, 0.0));
		poles.SetValue(2, gp_Pnt(1.0, 0.0, 0.0));
		poles.SetValue(3, gp_Pnt(10.0, 0.0, 0.0));
		const occ::handle<Geom_BezierCurve> curve = new Geom_BezierCurve(poles);
		return BRepBuilderAPI_MakeEdge(curve).Edge();
	}

	TopoDS_Edge reparameterized_line(double first, double last)
	{
		TColgp_Array1OfPnt poles(1, 2);
		poles.SetValue(1, gp_Pnt(0.0, 0.0, 0.0));
		poles.SetValue(2, gp_Pnt(10.0, 0.0, 0.0));
		const occ::handle<Geom_BezierCurve> bezier = new Geom_BezierCurve(poles);
		const occ::handle<Geom_BSplineCurve> bspline =
			GeomConvert::CurveToBSplineCurve(bezier);
		TColStd_Array1OfReal knots(1, bspline->NbKnots());
		bspline->Knots(knots);
		const double old_first = knots(knots.Lower());
		const double old_last = knots(knots.Upper());
		for (int knot = knots.Lower(); knot <= knots.Upper(); ++knot)
			knots.SetValue(knot, first + (knots(knot) - old_first)
				/ (old_last - old_first) * (last - first));
		bspline->SetKnots(knots);
		return BRepBuilderAPI_MakeEdge(bspline, first, last).Edge();
	}

	TopoDS_Edge closed_circle()
	{
		const occ::handle<Geom_Circle> circle = new Geom_Circle(
			gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 2.0);
		return BRepBuilderAPI_MakeEdge(circle, 0.0, 2.0 * std::numbers::pi).Edge();
	}

	OcctContactSourceEdge source(std::uint32_t edge_id, const TopoDS_Edge& edge,
		std::uint32_t owner_face)
	{
		OcctContactSourceEdge output;
		output.source_edge_id = edge_id;
		output.edge = edge;
		output.owner_face_uses.push_back({owner_face,
			OcctFaceIntervalLocation::boundary, 1, {}});
		return output;
	}

	void add_owner_occurrence(OcctContactSourceEdge& edge, std::uint32_t face_id,
		std::uint32_t occurrence_id)
	{
		auto use = std::find_if(edge.owner_face_uses.begin(), edge.owner_face_uses.end(),
			[&](const auto& candidate) { return candidate.source_face_id == face_id; });
		if (use == edge.owner_face_uses.end())
		{
			edge.owner_face_uses.push_back({face_id,
				OcctFaceIntervalLocation::boundary, 1, {}});
			use = std::prev(edge.owner_face_uses.end());
		}
		const auto range = range_of(edge.edge);
		OcctContactBoundaryOccurrence occurrence;
		occurrence.occurrence_id = occurrence_id;
		occurrence.target_source_edge_id = edge.source_edge_id;
		occurrence.orientation = static_cast<std::int8_t>(edge.edge.Orientation());
		occurrence.target_edge = edge.edge;
		occurrence.source_parameter_begin = range[0];
		occurrence.source_parameter_end = range[1];
		occurrence.target_parameter_begin = range[0];
		occurrence.target_parameter_end = range[1];
		use->boundary_occurrences.push_back(std::move(occurrence));
	}

	OcctContactEdgeFaceInterval interior_interval(const OcctContactSourceEdge& source_edge,
		std::uint32_t target_face, double normalized_begin, double normalized_end)
	{
		const auto source_range = range_of(source_edge.edge);
		OcctContactEdgeFaceInterval output;
		output.source_edge_id = source_edge.source_edge_id;
		output.target_face_id = target_face;
		output.source_parameter_first = source_range[0];
		output.source_parameter_last = source_range[1];
		output.interval.begin = normalized_begin;
		output.interval.end = normalized_end;
		output.interval.location = OcctFaceIntervalLocation::interior;
		output.interval.target_sector_count = 2;
		return output;
	}

	OcctContactEdgeFaceInterval boundary_interval(const OcctContactSourceEdge& source_edge,
		std::uint32_t target_face, double normalized_begin, double normalized_end,
		std::uint32_t occurrence_id, const TopoDS_Edge& target_edge,
		double source_mapping_begin, double source_mapping_end,
		std::uint32_t target_source_edge_id)
	{
		const auto source_range = range_of(source_edge.edge);
		const auto target_range = range_of(target_edge);
		OcctContactEdgeFaceInterval output;
		output.source_edge_id = source_edge.source_edge_id;
		output.target_face_id = target_face;
		output.source_parameter_first = source_range[0];
		output.source_parameter_last = source_range[1];
		output.interval.begin = normalized_begin;
		output.interval.end = normalized_end;
		output.interval.location = OcctFaceIntervalLocation::boundary;
		output.interval.target_sector_count = 1;
		output.interval.boundary_occurrence_count = 1;
		output.interval.target_boundary_occurrence_ids[0] = occurrence_id;
		output.interval.target_boundary_orientations[0] =
			static_cast<std::int8_t>(target_edge.Orientation());
		output.interval.target_parameter_begin[0] = target_range[0];
		output.interval.target_parameter_end[0] = target_range[1];
		output.interval.target_mapping_source_begin[0] = source_mapping_begin;
		output.interval.target_mapping_source_end[0] = source_mapping_end;
		output.target_boundary_occurrences[0] = target_edge;
		output.target_boundary_source_edge_ids[0] = target_source_edge_id;
		return output;
	}

	OcctContactReciprocalInterval reciprocal(const OcctContactSourceEdge& a,
		double a_begin, double a_end, const OcctContactSourceEdge& b,
		double b_begin, double b_end)
	{
		return {a.source_edge_id, b.source_edge_id, a_begin, a_end, b_begin, b_end, 0.0};
	}

	OcctContactAtomization run(const std::vector<OcctContactSourceEdge>& edges,
		const std::vector<OcctContactEdgeFaceInterval>& intervals = {},
		const std::vector<OcctContactReciprocalInterval>& overlaps = {},
		const std::vector<OcctContactExactJunction>& junctions = {})
	{
		OcctContactAtomization output;
		std::string error;
		const OcctContactAtomizerInput input{edges, intervals, overlaps, junctions};
		require(atomize_occt_contacts(input, output, error), error);
		return output;
	}

	const OcctContactAtom* atom_with_span(const OcctContactAtomization& result,
		std::uint32_t edge, double low, double high)
	{
		for (const auto& atom : result.atoms)
			for (const auto& span : atom.source_spans)
				if (span.source_edge_id == edge
					&& near(std::min(span.parameter_begin, span.parameter_end), low)
					&& near(std::max(span.parameter_begin, span.parameter_end), high)) return &atom;
		return nullptr;
	}

	std::optional<std::uint32_t> topology_at(const OcctContactAtomization& result,
		std::array<double, 3> point)
	{
		std::optional<std::uint32_t> found;
		for (const auto& atom : result.atoms)
			for (const auto& sample : atom.samples)
				if (near(sample.canonical_world_position[0], point[0])
					&& near(sample.canonical_world_position[1], point[1])
					&& near(sample.canonical_world_position[2], point[2]))
				{
					if (found) require(*found == sample.topology_id,
						"one exact point received inconsistent topology IDs");
					found = sample.topology_id;
				}
		return found;
	}

	bool same_atomization(const OcctContactAtomization& a,
		const OcctContactAtomization& b)
	{
		if (a.topology_node_count != b.topology_node_count
			|| a.atoms.size() != b.atoms.size()) return false;
		for (std::size_t atom = 0; atom < a.atoms.size(); ++atom)
		{
			const auto& lhs = a.atoms[atom];
			const auto& rhs = b.atoms[atom];
			if (lhs.id != rhs.id || lhs.fan_degree != rhs.fan_degree
				|| lhs.source_spans.size() != rhs.source_spans.size()
				|| lhs.face_uses.size() != rhs.face_uses.size()
				|| lhs.samples.size() != rhs.samples.size()) return false;
			for (std::size_t span = 0; span < lhs.source_spans.size(); ++span)
				if (lhs.source_spans[span].source_edge_id != rhs.source_spans[span].source_edge_id
					|| lhs.source_spans[span].parameter_begin != rhs.source_spans[span].parameter_begin
					|| lhs.source_spans[span].parameter_end != rhs.source_spans[span].parameter_end)
					return false;
			for (std::size_t use = 0; use < lhs.face_uses.size(); ++use)
			{
				if (lhs.face_uses[use].source_face_id != rhs.face_uses[use].source_face_id
					|| lhs.face_uses[use].location != rhs.face_uses[use].location
					|| lhs.face_uses[use].sector_count != rhs.face_uses[use].sector_count
					|| lhs.face_uses[use].exact_operation_tolerance
						!= rhs.face_uses[use].exact_operation_tolerance
					|| lhs.face_uses[use].boundary_occurrences.size()
						!= rhs.face_uses[use].boundary_occurrences.size()) return false;
				for (std::size_t occurrence = 0;
					occurrence < lhs.face_uses[use].boundary_occurrences.size(); ++occurrence)
				{
					const auto& left = lhs.face_uses[use].boundary_occurrences[occurrence];
					const auto& right = rhs.face_uses[use].boundary_occurrences[occurrence];
					if (left.source_edge_id != right.source_edge_id
						|| left.mapping_source_edge_id != right.mapping_source_edge_id
						|| left.occurrence_id != right.occurrence_id
						|| left.orientation != right.orientation
						|| left.source_parameter_begin != right.source_parameter_begin
						|| left.source_parameter_end != right.source_parameter_end
						|| left.target_parameter_begin != right.target_parameter_begin
						|| left.target_parameter_end != right.target_parameter_end
						|| left.sample_target_parameters != right.sample_target_parameters)
						return false;
				}
			}
			for (std::size_t sample = 0; sample < lhs.samples.size(); ++sample)
				if (lhs.samples[sample].canonical_world_position
						!= rhs.samples[sample].canonical_world_position
					|| lhs.samples[sample].canonical_position_tolerance
						!= rhs.samples[sample].canonical_position_tolerance
					|| lhs.samples[sample].topology_id != rhs.samples[sample].topology_id
					|| lhs.samples[sample].source_locations.size()
						!= rhs.samples[sample].source_locations.size()
					|| !std::equal(lhs.samples[sample].source_locations.begin(),
						lhs.samples[sample].source_locations.end(),
						rhs.samples[sample].source_locations.begin(), [](const auto& left,
							const auto& right)
							{
								return left.source_edge_id == right.source_edge_id
									&& left.parameter == right.parameter;
							})) return false;
		}
		return true;
	}

	void partial_long_edge()
	{
		const TopoDS_Edge target = line({2, 0, 0}, {8, 0, 0});
		std::vector<OcctContactSourceEdge> edges{
			source(10, line({0, 0, 0}, {10, 0, 0}), 0),
			source(11, target, 1)};
		add_owner_occurrence(edges[1], 1, 7);
		const auto range = range_of(edges[0].edge);
		const auto target_range = range_of(edges[1].edge);
		// An inherited tessellation node within native resolution of the exact split
		// must not move that split's canonical coordinate.
		edges[0].samples.push_back({range[0] + 0.2 * (range[1] - range[0])
			- 5.0e-8, 1.0e-7});
		std::vector<OcctContactEdgeFaceInterval> intervals{
			boundary_interval(edges[0], 1, 0.2, 0.8, 7, target,
				range[0] + 0.2 * (range[1] - range[0]),
				range[0] + 0.8 * (range[1] - range[0]), 11)};
		const auto result = run(edges, intervals,
			{reciprocal(edges[0], range[0] + 0.2 * (range[1] - range[0]),
				range[0] + 0.8 * (range[1] - range[0]), edges[1],
				target_range[0], target_range[1])});
		require(result.atoms.size() == 3, "partial long edge was not split into three atoms");
		require(result.atoms[0].fan_degree == 1 && result.atoms[0].is_open_boundary(),
			"partial-contact leading gap was not retained as fan one");
		require(result.atoms[1].fan_degree == 2 && !result.atoms[1].is_open_boundary(),
			"partial-contact covered interval has the wrong fan");
		require(result.atoms[2].fan_degree == 1 && result.atoms[2].is_open_boundary(),
			"partial-contact trailing gap was not retained as fan one");
		require(topology_at(result, {{2, 0, 0}}).has_value(),
			"nearby tessellation sample displaced the exact interval endpoint");
	}

	void one_long_two_subedges()
	{
		std::vector<OcctContactSourceEdge> edges{
			source(20, line({0, 0, 0}, {10, 0, 0}), 0),
			source(21, line({0, 0, 0}, {4, 0, 0}), 1),
			source(22, line({4, 0, 0}, {10, 0, 0}), 1)};
		const auto a = range_of(edges[0].edge);
		const auto b = range_of(edges[1].edge);
		const auto c = range_of(edges[2].edge);
		std::vector<OcctContactReciprocalInterval> overlaps{
			reciprocal(edges[0], a[0], a[0] + 0.4 * (a[1] - a[0]),
				edges[1], b[0], b[1]),
			reciprocal(edges[0], a[0] + 0.4 * (a[1] - a[0]), a[1],
				edges[2], c[0], c[1])};
		const auto result = run(edges, {}, overlaps);
		require(result.atoms.size() == 2, "one-long/two-subedge contact did not canonicalize to two atoms");
		for (const auto& atom : result.atoms)
			require(atom.source_spans.size() == 2 && atom.fan_degree == 2,
				"one-long/two-subedge atom lost a reciprocal span or face sector");
		require(topology_at(result, {{4, 0, 0}}).has_value(),
			"one-long/two-subedge handoff lost its shared topology node");
		std::reverse(edges.begin(), edges.end());
		std::reverse(overlaps.begin(), overlaps.end());
		const auto reordered = run(edges, {}, overlaps);
		require(same_atomization(result, reordered),
			"atomization depends on source/reciprocal input ordering");
	}

	void t_junction_fan_three()
	{
		std::vector<OcctContactSourceEdge> edges{
			source(30, line({0, 0, 0}, {10, 0, 0}), 0),
			source(31, line({5, 0, 0}, {5, 4, 0}), 1)};
		std::vector<OcctContactEdgeFaceInterval> intervals{
			interior_interval(edges[1], 0, 0.0, 1.0)};
		const auto skin_range = range_of(edges[0].edge);
		const auto rib_range = range_of(edges[1].edge);
		OcctContactExactJunction junction;
		junction.junction_id = 900;
		// The exact operation's canonical junction may lie inside, rather than at
		// the centre of, the incident edges' native tolerance tubes. Publishing the
		// replacement coordinate must retain the full accepted group bound.
		junction.canonical_world_position = {{5, 2.0e-5, 0}};
		junction.geometric_tolerance = 5.0e-5;
		junction.incidences = {{30, 0.5 * (skin_range[0] + skin_range[1])},
			{31, rib_range[0]}};
		const auto result = run(edges, intervals, {}, {junction});
		const OcctContactAtom* rib = atom_with_span(result, 31, rib_range[0], rib_range[1]);
		require(rib && rib->fan_degree == 3 && rib->face_uses.size() == 2,
			"interior T-junction did not form a three-sector atom");
		const auto topology = topology_at(result, {{5, 2.0e-5, 0}});
		require(topology.has_value(), "T-junction did not publish a canonical topology node");
		std::size_t copies = 0;
		for (const auto& atom : result.atoms)
			for (const auto& sample : atom.samples)
				if (sample.topology_id == *topology)
				{
					++copies;
					require(sample.canonical_position_tolerance >= 9.0e-5,
						"displaced exact-junction coordinate under-reported its "
						"canonical-to-member certificate bound");
				}
		require(copies == 3, "T-junction node was not shared by both skin atoms and the rib atom");
	}

	void shared_plus_target_fan_four()
	{
		std::vector<OcctContactSourceEdge> edges{
			source(40, line({0, 0, 0}, {10, 0, 0}), 0),
			source(41, line({10, 0, 0}, {0, 0, 0}), 1)};
		BRep_Builder tolerance_builder;
		tolerance_builder.UpdateEdge(edges[1].edge, 1.0e-5);
		const auto a = range_of(edges[0].edge);
		const auto b = range_of(edges[1].edge);
		OcctContactBoundaryOccurrence owner_occurrence;
		owner_occurrence.occurrence_id = 300;
		owner_occurrence.target_source_edge_id = 41;
		owner_occurrence.orientation = static_cast<std::int8_t>(edges[1].edge.Orientation());
		owner_occurrence.target_edge = edges[1].edge;
		owner_occurrence.source_parameter_begin = b[0];
		owner_occurrence.source_parameter_end = b[1];
		owner_occurrence.target_parameter_begin = b[0];
		owner_occurrence.target_parameter_end = b[1];
		edges[1].owner_face_uses[0].boundary_occurrences.push_back(owner_occurrence);
		auto target = boundary_interval(edges[0], 1, 0.0, 1.0, 300,
			edges[1].edge, a[0], a[1], 41);
		// The exact target projection is allowed to differ from the owner edge's
		// endpoint parameter within that edge's native CAD tolerance.  This mirrors
		// imported BReps whose reciprocal copies use opposite parameterizations.
		target.interval.target_parameter_begin[0] = b[1] - 1.0e-6;
		target.interval.target_parameter_end[0] = b[0];
		const auto result = run(edges, {interior_interval(edges[0], 2, 0.0, 1.0),
			target}, {reciprocal(edges[0], a[0], a[1], edges[1], b[1], b[0])});
		require(result.atoms.size() == 1 && result.atoms[0].fan_degree == 4,
			"shared edge plus interior target did not form a fan-four atom");
		require(result.atoms[0].face_uses.size() == 3,
			"fan-four atom double-counted or lost a participating face");
		const auto face1 = std::find_if(result.atoms[0].face_uses.begin(),
			result.atoms[0].face_uses.end(), [](const auto& use)
				{ return use.source_face_id == 1; });
		require(face1 != result.atoms[0].face_uses.end()
			&& face1->boundary_occurrences.size() == 1,
			"reciprocal owner and target certificate emitted duplicate chart branches");
		const auto& occurrence = face1->boundary_occurrences.front();
		require(occurrence.source_edge_id == 41
			&& occurrence.mapping_source_edge_id == 41
			&& near(occurrence.target_parameter_begin, b[1])
			&& near(occurrence.target_parameter_end, b[0]),
			"duplicate occurrence did not canonicalize to its direct target-edge branch");
	}

	void deliberate_gap()
	{
		std::vector<OcctContactSourceEdge> edges{
			source(50, line({0, 0, 0}, {10, 0, 0}), 0),
			source(51, line({0, 0, 0}, {4, 0, 0}), 1),
			source(52, line({6, 0, 0}, {10, 0, 0}), 1)};
		add_owner_occurrence(edges[1], 1, 51);
		add_owner_occurrence(edges[2], 1, 52);
		const auto range = range_of(edges[0].edge);
		const auto left_range = range_of(edges[1].edge);
		const auto right_range = range_of(edges[2].edge);
		const double p4 = range[0] + 0.4 * (range[1] - range[0]);
		const double p6 = range[0] + 0.6 * (range[1] - range[0]);
		std::vector<OcctContactEdgeFaceInterval> intervals{
			boundary_interval(edges[0], 1, 0.0, 0.4, 51,
				edges[1].edge, range[0], p4, 51),
			boundary_interval(edges[0], 1, 0.6, 1.0, 52,
				edges[2].edge, p6, range[1], 52)};
		const auto result = run(edges, intervals,
			{reciprocal(edges[0], range[0], p4, edges[1],
				left_range[0], left_range[1]),
			 reciprocal(edges[0], p6, range[1], edges[2],
				right_range[0], right_range[1])});
		require(result.atoms.size() == 3 && result.atoms[0].fan_degree == 2
			&& result.atoms[1].fan_degree == 1 && result.atoms[2].fan_degree == 2,
			"deliberate uncovered interval was healed or lost");
	}

	void boundary_occurrence_handoff()
	{
		std::vector<OcctContactSourceEdge> edges{
			source(60, line({0, 0, 0}, {10, 0, 0}), 0),
			source(61, line({0, 0, 0}, {5, 0, 0}), 1),
			source(62, line({5, 0, 0}, {10, 0, 0}), 1)};
		add_owner_occurrence(edges[1], 1, 100);
		add_owner_occurrence(edges[2], 1, 101);
		const auto range = range_of(edges[0].edge);
		const auto left_range = range_of(edges[1].edge);
		const auto right_range = range_of(edges[2].edge);
		const double middle = 0.5 * (range[0] + range[1]);
		std::vector<OcctContactEdgeFaceInterval> intervals{
			boundary_interval(edges[0], 1, 0.0, 0.5, 100,
				edges[1].edge, range[0], middle, 61),
			boundary_interval(edges[0], 1, 0.5, 1.0, 101,
				edges[2].edge, middle, range[1], 62)};
		const auto result = run(edges, intervals,
			{reciprocal(edges[0], range[0], middle, edges[1],
				left_range[0], left_range[1]),
			 reciprocal(edges[0], middle, range[1], edges[2],
				right_range[0], right_range[1])});
		require(result.atoms.size() == 2,
			"boundary occurrence handoff was incorrectly coalesced");
		for (std::size_t atom = 0; atom < 2; ++atom)
		{
			require(result.atoms[atom].fan_degree == 2
				&& result.atoms[atom].face_uses.size() == 2,
				"boundary handoff atom has the wrong face fan");
			const auto& target = result.atoms[atom].face_uses[1];
			require(target.boundary_occurrences.size() == 1
				&& target.boundary_occurrences[0].occurrence_id == 100 + atom,
				"boundary handoff lost its oriented occurrence identity");
		}
	}

	void reversed_source_orientation()
	{
		const TopoDS_Edge forward = line({0, 0, 0}, {10, 0, 0});
		TopoDS_Edge reversed = forward;
		reversed.Reverse();
		std::vector<OcctContactSourceEdge> forward_edges{source(70, forward, 0)};
		std::vector<OcctContactSourceEdge> reversed_edges{source(70, reversed, 0)};
		const auto a = run(forward_edges, {interior_interval(forward_edges[0], 1, 0.2, 0.8)});
		const auto b = run(reversed_edges, {interior_interval(reversed_edges[0], 1, 0.2, 0.8)});
		require(a.atoms.size() == b.atoms.size()
			&& a.topology_node_count == b.topology_node_count,
			"reversing a TopoDS source orientation changed atom counts");
		for (std::size_t atom = 0; atom < a.atoms.size(); ++atom)
		{
			require(a.atoms[atom].fan_degree == b.atoms[atom].fan_degree
				&& a.atoms[atom].source_spans.size() == b.atoms[atom].source_spans.size()
				&& a.atoms[atom].samples.size() == b.atoms[atom].samples.size(),
				"reversing a TopoDS source orientation changed atom structure");
			for (std::size_t sample = 0; sample < a.atoms[atom].samples.size(); ++sample)
				require(a.atoms[atom].samples[sample].canonical_world_position
					== b.atoms[atom].samples[sample].canonical_world_position
					&& a.atoms[atom].samples[sample].topology_id
						== b.atoms[atom].samples[sample].topology_id,
					"reversing source orientation changed canonical samples");
		}
	}

	void nonlinear_parameter_recovery()
	{
		std::vector<OcctContactSourceEdge> edges{
			source(80, line({0, 0, 0}, {10, 0, 0}), 0),
			source(81, nonlinear_line(), 1)};
		const auto linear_range = range_of(edges[0].edge);
		const auto nonlinear_range = range_of(edges[1].edge);
		add_owner_occurrence(edges[1], 3, 201);
		edges[0].samples.push_back({linear_range[0]
			+ 0.2 * (linear_range[1] - linear_range[0]), 0.0});
		std::vector<OcctContactEdgeFaceInterval> intervals{
			interior_interval(edges[0], 2, 0.0, 0.5),
			boundary_interval(edges[0], 3, 0.0, 1.0, 201, edges[1].edge,
				linear_range[0], linear_range[1], 81)};
		std::vector<OcctContactReciprocalInterval> overlaps{
			reciprocal(edges[0], linear_range[0], linear_range[1], edges[1],
				nonlinear_range[0], nonlinear_range[1])};
		const auto result = run(edges, intervals, overlaps);
		require(result.atoms.size() == 2,
			"nonlinear reciprocal parameterization did not receive the propagated split");
		const auto nonlinear_left = std::find_if(result.atoms[0].source_spans.begin(),
			result.atoms[0].source_spans.end(), [](const auto& span)
				{ return span.source_edge_id == 81; });
		require(nonlinear_left != result.atoms[0].source_spans.end(),
			"nonlinear reciprocal atom lost its source span");
		const double nonlinear_split = std::max(nonlinear_left->parameter_begin,
			nonlinear_left->parameter_end);
		require(nonlinear_split > 0.6 && nonlinear_split < 0.75
			&& !near(nonlinear_split, 0.5, 1.0e-3),
			"reciprocal split used endpoint-linear parameter interpolation");

		const auto face3 = std::find_if(result.atoms[0].face_uses.begin(),
			result.atoms[0].face_uses.end(), [](const auto& use)
				{ return use.source_face_id == 3; });
		require(face3 != result.atoms[0].face_uses.end()
			&& face3->boundary_occurrences.size() == 1,
			"nonlinear target boundary occurrence was not retained");
		const double target_split = face3->boundary_occurrences[0].target_parameter_end;
		require(target_split > 0.6 && target_split < 0.75
			&& !near(target_split, 0.5, 1.0e-3),
			"target boundary atom endpoint used linear parameter interpolation");
		const auto& sample_parameters =
			face3->boundary_occurrences[0].sample_target_parameters;
		require(sample_parameters.size() == result.atoms[0].samples.size()
			&& sample_parameters.size() == 3,
			"target occurrence did not receive parameters for the complete union chain");
		require(sample_parameters[1] > 0.3 && sample_parameters[1] < 0.5
			&& !near(sample_parameters[1], 0.2, 1.0e-3),
			"target occurrence interior sample used endpoint-linear interpolation");
	}

	void transitive_split_and_union_samples()
	{
		// IDs deliberately make B<->C sort before A<->B.  The x=3 split and A-only
		// x=1.7 tessellation sample must still traverse the complete A<->B<->C chain.
		std::vector<OcctContactSourceEdge> edges{
			source(102, line({0, 0, 0}, {10, 0, 0}), 0),
			source(101, nonlinear_line(), 1),
			source(100, line({0, 0, 0}, {10, 0, 0}), 2)};
		const auto a = range_of(edges[0].edge);
		const auto b = range_of(edges[1].edge);
		const auto c = range_of(edges[2].edge);
		edges[0].samples.push_back({a[0] + 0.17 * (a[1] - a[0]), 0.0});
		std::vector<OcctContactReciprocalInterval> overlaps{
			reciprocal(edges[1], b[0], b[1], edges[2], c[0], c[1]),
			reciprocal(edges[0], a[0], a[1], edges[1], b[0], b[1])};
		const auto result = run(edges,
			{interior_interval(edges[0], 3, 0.0, 0.3)}, overlaps);
		require(result.atoms.size() == 2,
			"transitive reciprocal chain did not propagate the terminal split to every copy");
		for (const auto& atom : result.atoms)
			require(atom.source_spans.size() == 3,
				"transitive physical atom did not canonicalize all three source copies");
		const auto sample = std::find_if(result.atoms[0].samples.begin(),
			result.atoms[0].samples.end(), [](const auto& value)
				{ return near(value.canonical_world_position[0], 1.7); });
		require(sample != result.atoms[0].samples.end()
			&& sample->source_locations.size() == 3,
			"asymmetric tessellation sample was not recovered on every reciprocal member");
		std::vector<std::uint32_t> sampled_edges;
		for (const auto& location : sample->source_locations)
			sampled_edges.push_back(location.source_edge_id);
		require(sampled_edges == std::vector<std::uint32_t>({100, 101, 102}),
			"union-chain sample source locations are incomplete or nondeterministic");

		std::reverse(edges.begin(), edges.end());
		std::reverse(overlaps.begin(), overlaps.end());
		const auto reordered = run(edges,
			{interior_interval(edges[2], 3, 0.0, 0.3)}, overlaps);
		require(same_atomization(result, reordered),
			"transitive split/sample result depends on input ordering");
	}

	void closed_periodic_endpoint_junction()
	{
		std::vector<OcctContactSourceEdge> edges{
			source(105, closed_circle(), 0)};
		const auto range = range_of(edges[0].edge);
		edges[0].samples.push_back({range[0] + 0.25 * (range[1] - range[0]), 0.0});
		edges[0].owner_face_uses[0].location = OcctFaceIntervalLocation::interior;
		edges[0].owner_face_uses[0].sector_count = 2;

		OcctContactBoundaryOccurrence forward;
		forward.occurrence_id = 400;
		forward.target_source_edge_id = 105;
		forward.orientation = static_cast<std::int8_t>(TopAbs_FORWARD);
		forward.target_edge = edges[0].edge;
		forward.source_parameter_begin = range[0];
		forward.source_parameter_end = range[1];
		forward.target_parameter_begin = range[0];
		forward.target_parameter_end = range[1];
		OcctContactBoundaryOccurrence reverse = forward;
		reverse.occurrence_id = 401;
		reverse.orientation = static_cast<std::int8_t>(TopAbs_REVERSED);
		reverse.target_edge.Reverse();
		reverse.target_parameter_begin = range[1];
		reverse.target_parameter_end = range[0];
		edges[0].owner_face_uses[0].boundary_occurrences = {forward, reverse};

		BRepAdaptor_Curve curve(edges[0].edge);
		const gp_Pnt endpoint = curve.Value(range[0]);
		OcctContactExactJunction junction;
		junction.junction_id = 910;
		junction.canonical_world_position = {{endpoint.X(), endpoint.Y(), endpoint.Z()}};
		junction.incidences = {{105, range[0]}, {105, range[1]}};
		const auto result = run(edges, {}, {}, {junction});
		require(result.atoms.size() == 1 && result.atoms[0].fan_degree == 2
			&& result.atoms[0].samples.size() == 3,
			"closed periodic edge was not retained as one two-sector loop atom");
		require(result.atoms[0].samples.front().topology_id
				== result.atoms[0].samples.back().topology_id,
			"closed edge's distinct first/last parameters did not share one topology node");
		require(result.atoms[0].face_uses.size() == 1
			&& result.atoms[0].face_uses[0].boundary_occurrences.size() == 2,
			"periodic seam lost one oriented pcurve branch");
		for (const auto& occurrence : result.atoms[0].face_uses[0].boundary_occurrences)
			require(occurrence.sample_target_parameters.size() == 3,
				"periodic pcurve branch lacks an aligned parameter for each loop endpoint");
	}

	void parameter_scale_resolution()
	{
		std::vector<OcctContactSourceEdge> tiny{
			source(110, reparameterized_line(0.0, 1.0e-12), 0)};
		const auto tiny_result = run(tiny, {interior_interval(tiny[0], 1, 0.25, 0.75)});
		require(tiny_result.atoms.size() == 3,
			"fixed absolute parameter tolerance erased a tiny finite source span");

		const double offset = 1.0e12;
		std::vector<OcctContactSourceEdge> shifted{
			source(111, reparameterized_line(offset, offset + 16.0), 0)};
		const auto shifted_result = run(shifted,
			{interior_interval(shifted[0], 1, 0.25, 0.75)});
		require(shifted_result.atoms.size() == 3,
			"absolute parameter magnitude erased a resolved large-offset source span");
	}

	void local_operation_bound_and_canonical_member_bound()
	{
		std::vector<OcctContactSourceEdge> single{
			source(120, line({0, 0, 0}, {10, 0, 0}), 0)};
		auto certified = interior_interval(single[0], 1, 0.25, 0.75);
		certified.interval.exact_operation_tolerance = 7.5e-6;
		const auto split = run(single, {certified});
		require(split.atoms.size() == 3,
			"operation-bound interval did not retain its two exact split points");
		for (const auto& atom : split.atoms)
		{
			const auto owner = std::find_if(atom.face_uses.begin(), atom.face_uses.end(),
				[](const auto& use) { return use.source_face_id == 0; });
			require(owner != atom.face_uses.end()
				&& owner->exact_operation_tolerance == 0.0,
				"target operation bound leaked into the source owner face use");
			const auto target = std::find_if(atom.face_uses.begin(), atom.face_uses.end(),
				[](const auto& use) { return use.source_face_id == 1; });
			const bool covered = atom.fan_degree == 3;
			require((target != atom.face_uses.end()) == covered,
				"target face use does not match its certified subspan");
			if (target != atom.face_uses.end())
			{
				require(target->exact_operation_tolerance == 7.5e-6,
					"exact operation bound was not isolated on its target face use");
			}
		}

		std::vector<OcctContactSourceEdge> copies{
			source(121, line({0, 0, 0}, {10, 0, 0}), 0),
			source(122, line({0, 0, 0}, {10, 0, 0}), 1)};
		BRep_Builder tolerance_builder;
		tolerance_builder.UpdateEdge(copies[0].edge, 2.0e-5);
		tolerance_builder.UpdateEdge(copies[1].edge, 3.0e-5);
		const auto a = range_of(copies[0].edge);
		const auto b = range_of(copies[1].edge);
		copies[1].samples.push_back({b[0] + 0.3 * (b[1] - b[0]), 0.0});
		const auto joined = run(copies, {},
			{reciprocal(copies[0], a[0], a[1], copies[1], b[0], b[1])});
		const auto sample = std::find_if(joined.atoms.front().samples.begin(),
			joined.atoms.front().samples.end(), [](const auto& value)
				{ return near(value.canonical_world_position[0], 3.0); });
		require(sample != joined.atoms.front().samples.end(),
			"reciprocal native-bound sample was not retained");
		require(sample->canonical_position_tolerance >= 5.0e-5,
			"canonical sample published only one member's native tolerance instead of "
			"the proven canonical-to-member bound");
	}

	void transactional_rejection()
	{
		std::vector<OcctContactSourceEdge> edges{
			source(90, line({0, 0, 0}, {10, 0, 0}), 0),
			source(91, line({0, 1, 0}, {10, 1, 0}), 1)};
		const auto a = range_of(edges[0].edge);
		const auto b = range_of(edges[1].edge);
		OcctContactAtomization output;
		output.topology_node_count = 123;
		output.atoms.push_back({999});
		std::string error;
		auto invalid_operation = interior_interval(edges[0], 2, 0.2, 0.8);
		invalid_operation.interval.exact_operation_tolerance =
			std::numeric_limits<double>::quiet_NaN();
		const std::vector<OcctContactEdgeFaceInterval> invalid_intervals{invalid_operation};
		const OcctContactAtomizerInput invalid_operation_input{
			edges, invalid_intervals, {}, {}};
		require(!atomize_occt_contacts(invalid_operation_input, output, error)
			&& !error.empty(),
			"non-finite exact-operation certificate did not fail closed");
		require(output.topology_node_count == 123 && output.atoms.size() == 1
			&& output.atoms[0].id == 999,
			"invalid exact-operation certificate modified caller-owned output");
		error.clear();
		const std::vector<OcctContactReciprocalInterval> overlaps{
			reciprocal(edges[0], a[0], a[1], edges[1], b[0], b[1])};
		const OcctContactAtomizerInput input{edges, {}, overlaps, {}};
		require(!atomize_occt_contacts(input, output, error) && !error.empty(),
			"invalid reciprocal certificate did not fail closed");
		require(output.topology_node_count == 123 && output.atoms.size() == 1
			&& output.atoms[0].id == 999,
			"failed atomization modified the caller's existing output");
	}
}

int main()
{
	partial_long_edge();
	std::puts("[occt-contact-atomizer] partial long edge + fan-one gaps: PASS");
	one_long_two_subedges();
	std::puts("[occt-contact-atomizer] one long / two reciprocal subedges: PASS");
	t_junction_fan_three();
	std::puts("[occt-contact-atomizer] exact T-junction fan three: PASS");
	shared_plus_target_fan_four();
	std::puts("[occt-contact-atomizer] shared edge + interior target fan four: PASS");
	deliberate_gap();
	std::puts("[occt-contact-atomizer] deliberate gap remains fan one: PASS");
	boundary_occurrence_handoff();
	std::puts("[occt-contact-atomizer] boundary occurrence handoff: PASS");
	reversed_source_orientation();
	std::puts("[occt-contact-atomizer] reversed source orientation: PASS");
	nonlinear_parameter_recovery();
	std::puts("[occt-contact-atomizer] nonlinear exact parameter recovery: PASS");
	transitive_split_and_union_samples();
	std::puts("[occt-contact-atomizer] transitive splits + asymmetric union samples: PASS");
	closed_periodic_endpoint_junction();
	std::puts("[occt-contact-atomizer] closed periodic endpoint junction + branches: PASS");
	parameter_scale_resolution();
	std::puts("[occt-contact-atomizer] tiny/large-offset parameter resolution: PASS");
	local_operation_bound_and_canonical_member_bound();
	std::puts("[occt-contact-atomizer] local operation + canonical-member bounds: PASS");
	transactional_rejection();
	std::puts("[occt-contact-atomizer] transactional rejection: PASS");
	std::puts("[occt-contact-atomizer] all manufactured checks passed");
	return 0;
}
