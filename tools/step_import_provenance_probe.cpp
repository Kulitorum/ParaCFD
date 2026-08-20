#include "core/geometry/model_placement.h"
#include "core/geometry/mesh_clip.h"
#include "core/geometry/step_import.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <map>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

using namespace paracfd::core;

namespace
{
	int failures = 0;

	void check(bool condition, const char* label)
	{
		std::printf("%s  %s\n", condition ? "PASS" : "FAIL", label);
		if (!condition) ++failures;
	}

	std::uint64_t edge_key(std::uint32_t a, std::uint32_t b)
	{
		if (a > b) std::swap(a, b);
		return (static_cast<std::uint64_t>(a) << 32) | b;
	}

	bool positions_are_bit_identical(const std::array<double, 3>& a,
		const std::array<double, 3>& b)
	{
		for (std::size_t axis = 0; axis < 3; ++axis)
			if (std::bit_cast<std::uint64_t>(a[axis])
				!= std::bit_cast<std::uint64_t>(b[axis])) return false;
		return true;
	}

	struct TopologySummary
	{
		bool structurally_valid = true;
		std::size_t boundary_segments = 0;
		std::size_t labelled_boundary_segments = 0;
		std::size_t unknown_boundary_segments = 0;
		std::size_t labelled_interior_segments = 0;
		std::set<std::pair<std::uint32_t, std::uint32_t>> face_edge_uses;
		std::map<std::uint32_t, std::uint32_t> edge_incidence;
		std::map<std::uint32_t, double> edge_tolerance_m;
		std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> periodic_seam_segments;
	};

	TopologySummary summarize(const TriMesh& mesh)
	{
		TopologySummary result;
		if (!mesh.has_face_provenance() || !mesh.has_cad_edge_provenance())
		{
			result.structurally_valid = false;
			return result;
		}

		using EdgeCounts = std::unordered_map<std::uint64_t, unsigned>;
		std::unordered_map<std::uint32_t, EdgeCounts> counts_by_face;
		for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
		{
			const std::uint32_t face = mesh.source_face_ids[triangle];
			for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
			{
				const std::uint32_t a = mesh.indices[3 * triangle + half_edge];
				const std::uint32_t b = mesh.indices[3 * triangle + (half_edge + 1) % 3];
				if (a >= mesh.vertex_count() || b >= mesh.vertex_count() || a == b)
				{
					result.structurally_valid = false;
					continue;
				}
				++counts_by_face[face][edge_key(a, b)];
			}
		}

		for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
		{
			const std::uint32_t face = mesh.source_face_ids[triangle];
			for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
			{
				const std::uint32_t a = mesh.indices[3 * triangle + half_edge];
				const std::uint32_t b = mesh.indices[3 * triangle + (half_edge + 1) % 3];
				const unsigned occurrences = counts_by_face[face][edge_key(a, b)];
				const std::uint32_t edge = mesh.cad_edge_id(triangle, half_edge);
				const CadEdgeProvenanceState state =
					mesh.cad_edge_provenance_state(triangle, half_edge);
				if (state != CadEdgeProvenanceState::none
					&& state != CadEdgeProvenanceState::known
					&& state != CadEdgeProvenanceState::unknown_boundary)
					result.structurally_valid = false;
				const std::uint32_t incidence = mesh.cad_edge_incident_face_count(triangle, half_edge);
				const double tolerance_m = mesh.cad_edge_tolerance(triangle, half_edge);
				const bool periodic_seam = mesh.cad_edge_is_periodic_seam(triangle, half_edge);
				if (occurrences == 1)
				{
					++result.boundary_segments;
					if (state == CadEdgeProvenanceState::known) ++result.labelled_boundary_segments;
					else if (state == CadEdgeProvenanceState::unknown_boundary)
						++result.unknown_boundary_segments;
					else result.structurally_valid = false;
				}
				else if (edge != TriMesh::kNoCadEdgeId)
				{
					++result.labelled_interior_segments;
				}
				if (state != CadEdgeProvenanceState::known)
				{
					if (edge != TriMesh::kNoCadEdgeId || incidence != 0 || tolerance_m != 0.0 || periodic_seam)
						result.structurally_valid = false;
					continue;
				}
				if (edge == TriMesh::kNoCadEdgeId) result.structurally_valid = false;
				if (incidence == 0 || !std::isfinite(tolerance_m) || tolerance_m < 0.0)
					result.structurally_valid = false;
				result.face_edge_uses.emplace(face, edge);
				auto [found, inserted] = result.edge_incidence.emplace(edge, incidence);
				if (!inserted && found->second != incidence) result.structurally_valid = false;
				auto [tolerance, tolerance_inserted] = result.edge_tolerance_m.emplace(edge, tolerance_m);
				if (!tolerance_inserted && tolerance->second != tolerance_m) result.structurally_valid = false;
				if (periodic_seam)
				{
					if (incidence < 2) result.structurally_valid = false;
					++result.periodic_seam_segments[{face, edge}];
				}
			}
		}
		for (const auto& [face_edge, segments] : result.periodic_seam_segments)
			if (segments < 2) result.structurally_valid = false;
		return result;
	}

	bool contact_graph_is_structurally_valid(const StepCadContactGraph& graph)
	{
		auto span_key=[](const StepCadContactGraph::UnresolvedEdgeFaceSpan& span)
		{
			return std::tuple{span.source_edge_id,span.target_face_id,
				span.source_parameter_begin,span.source_parameter_end,
				span.source_parameter_tolerance,span.endpoint_begin_m,span.endpoint_end_m,
				static_cast<unsigned>(span.kind),static_cast<unsigned>(span.target_sector_count),
				static_cast<unsigned>(span.target_boundary_occurrence_count),
				span.target_boundary_occurrence_ids,span.target_boundary_orientations,
				span.target_parameter_begin,span.target_parameter_end,
				static_cast<unsigned>(span.issue),span.reason};
		};
		if(!std::is_sorted(graph.unresolved_edge_face_spans.begin(),
			graph.unresolved_edge_face_spans.end(),[&](const auto& a,const auto& b)
			{return span_key(a)<span_key(b);}))return false;
		for(std::size_t span_id=0;span_id<graph.unresolved_edge_face_spans.size();++span_id)
		{
			const auto& span=graph.unresolved_edge_face_spans[span_id];
			if(span.source_edge_id==TriMesh::kNoCadEdgeId||span.reason.empty()
				||!std::isfinite(span.source_parameter_begin)
				||!std::isfinite(span.source_parameter_end)
				||!(span.source_parameter_end>span.source_parameter_begin)
				||!std::isfinite(span.source_parameter_tolerance)
				||span.source_parameter_tolerance<0.0)return false;
			for(const double coordinate:span.endpoint_begin_m)if(!std::isfinite(coordinate))return false;
			for(const double coordinate:span.endpoint_end_m)if(!std::isfinite(coordinate))return false;
			if(span.target_boundary_occurrence_count>2)return false;
			for(std::size_t occurrence=0;occurrence<2;++occurrence)
			{
				const bool active=occurrence<span.target_boundary_occurrence_count;
				if(active)
				{
					if(span.target_boundary_occurrence_ids[occurrence]
						==std::numeric_limits<std::uint32_t>::max()
						||!std::isfinite(span.target_parameter_begin[occurrence])
						||!std::isfinite(span.target_parameter_end[occurrence]))return false;
					const auto orientation=span.target_boundary_orientations[occurrence];
					if(orientation!=0&&orientation!=1)return false;
				}
				else if(span.target_boundary_occurrence_ids[occurrence]
					!=std::numeric_limits<std::uint32_t>::max())return false;
			}
			if(span.kind==StepEdgeFaceSpanKind::unclassified)
			{
				if(span.target_sector_count!=0
					||span.target_boundary_occurrence_count!=0
					||span.issue!=StepEdgeFaceSpanIssue::certification_failure)return false;
			}
			else if(span.target_sector_count<1||span.target_sector_count>2
				||(span.kind==StepEdgeFaceSpanKind::trim_boundary
					&&span.target_boundary_occurrence_count!=1)
				||(span.kind==StepEdgeFaceSpanKind::face_interior
					&&span.target_boundary_occurrence_count!=0
					&&span.target_boundary_occurrence_count!=2))return false;
			if(span_id&&span_key(graph.unresolved_edge_face_spans[span_id-1])==span_key(span))
				return false;
		}
		using FailureKey = std::tuple<std::uint64_t, std::uint32_t, StepContactUseKind>;
		std::set<FailureKey> unresolved_uses;
		for (const StepCadContactGraph::UvProjectionFailure& failure :
			graph.unresolved_uv_projections)
		{
			if (failure.curve_id == TriMesh::kNoCadContactId || failure.detail.empty()
				|| !unresolved_uses.emplace(failure.curve_id, failure.source_face_id,
					failure.kind).second) return false;
		}
		std::set<FailureKey> matched_failures;
		struct TopologyPoint
		{
			std::array<double, 3> position{};
			double tolerance_m = 0.0;
		};
		std::map<std::uint32_t, TopologyPoint> topology_points;
		std::optional<std::uint64_t> previous_atom_id;
		for (const StepContactCurve& curve : graph.curves)
		{
			if ((previous_atom_id && curve.id <= *previous_atom_id)
				|| curve.id == TriMesh::kNoCadContactId
				|| curve.source_edge_id == TriMesh::kNoCadEdgeId
				|| curve.fan_degree < 2 || curve.sample_positions_m.size() < 2
				|| curve.source_parameters.size() != curve.sample_positions_m.size()
				|| curve.sample_topology_ids.size() != curve.sample_positions_m.size()
				|| curve.source_edge_ids.empty()
				|| !std::is_sorted(curve.source_edge_ids.begin(), curve.source_edge_ids.end())
				|| std::adjacent_find(curve.source_edge_ids.begin(), curve.source_edge_ids.end())
					!= curve.source_edge_ids.end()
				|| curve.source_edge_id != curve.source_edge_ids.front()
				|| curve.uses.empty() || !std::isfinite(curve.tolerance_m)
				|| curve.tolerance_m < 0.0) return false;
			previous_atom_id = curve.id;
			for (std::uint32_t source_edge : curve.source_edge_ids)
				if (source_edge == TriMesh::kNoCadEdgeId) return false;
			bool parameters_increasing = true, parameters_decreasing = true;
			std::map<std::uint32_t, std::size_t> curve_topology_ids;
			for (std::size_t sample = 0; sample < curve.sample_positions_m.size(); ++sample)
			{
				if (!std::isfinite(curve.source_parameters[sample])) return false;
				if (sample > 0)
				{
					parameters_increasing = parameters_increasing
						&& curve.source_parameters[sample] > curve.source_parameters[sample - 1];
					parameters_decreasing = parameters_decreasing
						&& curve.source_parameters[sample] < curve.source_parameters[sample - 1];
				}
				const auto& point = curve.sample_positions_m[sample];
				if (!std::isfinite(point[0]) || !std::isfinite(point[1])
					|| !std::isfinite(point[2])) return false;
				auto [first_sample, first_on_curve] = curve_topology_ids.emplace(
					curve.sample_topology_ids[sample], sample);
				if (!first_on_curve)
				{
					// A closed atom deliberately repeats its first topology node at the
					// end. No other within-atom topology-node repetition is valid.
					const bool intentional_closure = first_sample->second == 0
						&& sample + 1 == curve.sample_positions_m.size()
						&& positions_are_bit_identical(point,
							curve.sample_positions_m.front());
					if (!intentional_closure) return false;
				}
				auto [existing, inserted] = topology_points.emplace(curve.sample_topology_ids[sample],
					TopologyPoint{ point, curve.tolerance_m });
				if (!inserted)
				{
					// A shared topology node must already be canonicalized exactly. Downstream
					// BVH/EB code deliberately refuses to turn native CAD tolerance into an
					// accidental geometry weld.
					if (!positions_are_bit_identical(point, existing->second.position))
						return false;
				}
			}
			if (!parameters_increasing && !parameters_decreasing) return false;
			std::uint32_t sectors = 0;
			std::set<std::uint32_t> faces;
			for (const StepContactFaceUse& use : curve.uses)
			{
				if (!faces.insert(use.source_face_id).second || use.sector_count < 1
					|| use.sector_count > 2) return false;
				const FailureKey use_key{curve.id, use.source_face_id, use.kind};
				const bool unresolved = unresolved_uses.contains(use_key);
				if (unresolved)
				{
					if (!use.sample_uv.empty()) return false;
					matched_failures.insert(use_key);
				}
				else if (use.sample_uv.size() != curve.sample_positions_m.size())
					return false;
				sectors += use.sector_count;
				for (const auto& uv : use.sample_uv)
					if (!std::isfinite(uv[0]) || !std::isfinite(uv[1])) return false;
			}
			if (sectors != curve.fan_degree) return false;
		}
		std::uint32_t expected_topology_id = 0;
		for (const auto& [topology_id, point] : topology_points)
		{
			(void)point;
			if (topology_id != expected_topology_id++) return false;
		}
		return matched_failures == unresolved_uses;
	}

	bool mesh_contact_labels_match_graph(const TriMesh& mesh, const StepCadContactGraph& graph)
	{
		if (!mesh.has_cad_edge_contact_provenance() || !mesh.has_topology_vertex_ids()
			|| !mesh.has_fp64_positions()) return false;

		using ContactSegment = std::tuple<std::uint64_t, std::uint32_t, std::uint32_t>;
		std::map<std::uint64_t, std::uint32_t> fan_by_curve;
		std::map<ContactSegment, std::uint32_t> expected_segment_incidence;
		std::map<ContactSegment, std::uint32_t> actual_segment_incidence;
		std::map<std::uint32_t, std::array<double, 3>> graph_topology_positions;
		for (const StepContactCurve& curve : graph.curves)
		{
			if (!fan_by_curve.emplace(curve.id, curve.fan_degree).second) return false;
			for (std::size_t sample = 0; sample < curve.sample_positions_m.size(); ++sample)
			{
				auto [position, inserted] = graph_topology_positions.emplace(
					curve.sample_topology_ids[sample], curve.sample_positions_m[sample]);
				if (!inserted && !positions_are_bit_identical(position->second,
					curve.sample_positions_m[sample])) return false;
			}
			for (std::size_t segment = 1; segment < curve.sample_topology_ids.size(); ++segment)
			{
				std::uint32_t a = curve.sample_topology_ids[segment - 1];
				std::uint32_t b = curve.sample_topology_ids[segment];
				if (a == b) return false;
				if (a > b) std::swap(a, b);
				if (!expected_segment_incidence.emplace(ContactSegment{ curve.id, a, b },
					curve.fan_degree).second) return false;
			}
		}

		std::map<std::uint32_t, std::array<double, 3>> mesh_topology_positions;
		for (std::size_t vertex = 0; vertex < mesh.vertex_count(); ++vertex)
		{
			const std::uint32_t topology = mesh.topology_vertex_id(vertex);
			const std::array<double, 3> point = mesh.vertex_position_double(vertex);
			auto [position, inserted] = mesh_topology_positions.emplace(topology, point);
			if (!inserted && !positions_are_bit_identical(position->second, point)) return false;
		}
		for (const auto& [topology, point] : graph_topology_positions)
		{
			const auto found = mesh_topology_positions.find(topology);
			if (found == mesh_topology_positions.end()
				|| !positions_are_bit_identical(found->second, point)) return false;
		}

		std::set<std::uint64_t> seen_contacts;
		for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
			for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
			{
				const std::uint64_t contact = mesh.cad_edge_contact_id(triangle, half_edge);
				const std::uint32_t fan = mesh.cad_edge_certified_fan_degree(triangle, half_edge);
				if (contact == TriMesh::kNoCadContactId)
				{
					if (fan != 0) return false;
					continue;
				}
				const auto found = fan_by_curve.find(contact);
				if (found == fan_by_curve.end() || found->second != fan) return false;
				const std::uint32_t vertex_a = mesh.indices[3 * triangle + half_edge];
				const std::uint32_t vertex_b = mesh.indices[3 * triangle + (half_edge + 1) % 3];
				if (vertex_a >= mesh.vertex_count() || vertex_b >= mesh.vertex_count()) return false;
				std::uint32_t topology_a = mesh.topology_vertex_id(vertex_a);
				std::uint32_t topology_b = mesh.topology_vertex_id(vertex_b);
				if (topology_a == topology_b) return false;
				if (topology_a > topology_b) std::swap(topology_a, topology_b);
				const ContactSegment segment{ contact, topology_a, topology_b };
				if (!expected_segment_incidence.contains(segment)) return false;
				++actual_segment_incidence[segment];
				seen_contacts.insert(contact);
			}
		return actual_segment_incidence == expected_segment_incidence
			&& seen_contacts.size() == fan_by_curve.size()
			&& std::all_of(fan_by_curve.begin(), fan_by_curve.end(),
				[&](const auto& entry) { return seen_contacts.contains(entry.first); });
	}

	struct ContactGraphCounts
	{
		std::size_t source_edges = 0;
		std::size_t samples = 0;
		std::size_t shared_topology_nodes = 0;
	};

	ContactGraphCounts contact_graph_counts(const StepCadContactGraph& graph)
	{
		ContactGraphCounts result;
		std::map<std::uint32_t, std::set<std::uint64_t>> topology_atom_uses;
		for (const StepContactCurve& curve : graph.curves)
		{
			result.source_edges += curve.source_edge_ids.size();
			result.samples += curve.sample_positions_m.size();
			for (std::uint32_t topology : curve.sample_topology_ids)
				topology_atom_uses[topology].insert(curve.id);
		}
		result.shared_topology_nodes = static_cast<std::size_t>(std::count_if(
			topology_atom_uses.begin(), topology_atom_uses.end(),
			[](const auto& entry) { return entry.second.size() > 1; }));
		return result;
	}

	StepContactFaceUse manufactured_contact_use(std::uint32_t face,
		std::size_t sample_count)
	{
		StepContactFaceUse use;
		use.source_face_id = face;
		use.kind = StepContactUseKind::trim_boundary;
		use.sector_count = 1;
		use.sample_uv.resize(sample_count);
		for (std::size_t sample = 0; sample < sample_count; ++sample)
			use.sample_uv[sample] = { static_cast<double>(sample), static_cast<double>(face) };
		return use;
	}

	StepContactCurve manufactured_contact_atom(std::uint64_t atom_id,
		std::uint32_t source_edge_id, std::vector<double> parameters,
		std::vector<std::array<double, 3>> positions,
		std::vector<std::uint32_t> topology_ids)
	{
		StepContactCurve atom;
		atom.id = atom_id;
		atom.source_edge_id = source_edge_id;
		atom.source_edge_ids = { source_edge_id };
		atom.tolerance_m = 1e-9;
		atom.fan_degree = 2;
		atom.source_parameters = std::move(parameters);
		atom.sample_positions_m = std::move(positions);
		atom.sample_topology_ids = std::move(topology_ids);
		atom.uses = { manufactured_contact_use(0, atom.sample_positions_m.size()),
			manufactured_contact_use(1, atom.sample_positions_m.size()) };
		return atom;
	}

	struct ExpectedEdgeFaceSpan
	{
		std::uint32_t edge=0;
		std::uint32_t face=0;
		double begin=0.0;
		double end=0.0;
		StepEdgeFaceSpanKind kind=StepEdgeFaceSpanKind::unclassified;
		std::uint8_t sectors=0;
		StepEdgeFaceSpanIssue issue=StepEdgeFaceSpanIssue::certification_failure;
	};

	std::optional<StepEdgeFaceSpanKind> parse_span_kind(const std::string& value)
	{
		if(value=="boundary")return StepEdgeFaceSpanKind::trim_boundary;
		if(value=="interior")return StepEdgeFaceSpanKind::face_interior;
		if(value=="unclassified")return StepEdgeFaceSpanKind::unclassified;
		return std::nullopt;
	}

	std::optional<StepEdgeFaceSpanIssue> parse_span_issue(const std::string& value)
	{
		if(value=="partial")return StepEdgeFaceSpanIssue::partial_coverage;
		if(value=="mixed")return StepEdgeFaceSpanIssue::mixed_sector_count;
		if(value=="handoff")return StepEdgeFaceSpanIssue::boundary_occurrence_handoff;
		if(value=="failure")return StepEdgeFaceSpanIssue::certification_failure;
		return std::nullopt;
	}

	bool expected_span_matches(const ExpectedEdgeFaceSpan& expected,
		const StepCadContactGraph::UnresolvedEdgeFaceSpan& actual)
	{
		const double arithmetic=256.0*std::numeric_limits<double>::epsilon()
			*std::max({1.0,std::abs(expected.begin),std::abs(expected.end),
				std::abs(actual.source_parameter_begin),std::abs(actual.source_parameter_end)});
		const double tolerance=std::max(actual.source_parameter_tolerance,arithmetic);
		return actual.source_edge_id==expected.edge&&actual.target_face_id==expected.face
			&&actual.kind==expected.kind&&actual.target_sector_count==expected.sectors
			&&actual.issue==expected.issue
			&&std::abs(actual.source_parameter_begin-expected.begin)<=tolerance
			&&std::abs(actual.source_parameter_end-expected.end)<=tolerance;
	}

	auto unresolved_span_signature(const StepCadContactGraph::UnresolvedEdgeFaceSpan& span)
	{
		return std::tuple{span.source_edge_id,span.target_face_id,
			span.source_parameter_begin,span.source_parameter_end,
			span.source_parameter_tolerance,span.endpoint_begin_m,span.endpoint_end_m,
			static_cast<unsigned>(span.kind),static_cast<unsigned>(span.target_sector_count),
			static_cast<unsigned>(span.target_boundary_occurrence_count),
			span.target_boundary_occurrence_ids,span.target_boundary_orientations,
			span.target_parameter_begin,span.target_parameter_end,
			static_cast<unsigned>(span.issue),span.reason};
	}

	bool exact_contact_topology_equal(const StepCadContactGraph& a,
		const StepCadContactGraph& b,std::string* mismatch=nullptr,
		bool compare_uv_failures=true)
	{
		auto fail=[&](std::string detail)
		{
			if(mismatch)*mismatch=std::move(detail);
			return false;
		};
		if(a.audit_status!=b.audit_status||a.audit_failure_reason!=b.audit_failure_reason)
			return fail("contact audit status/reason mismatch");
		if(a.curves.size()!=b.curves.size()
			||a.unresolved_edge_face_spans.size()!=b.unresolved_edge_face_spans.size()
			||a.unresolved_partial_overlaps.size()!=b.unresolved_partial_overlaps.size()
			||(compare_uv_failures
				&&a.unresolved_uv_projections.size()!=b.unresolved_uv_projections.size()))
		{
			std::ostringstream detail;
			detail<<"count mismatch: curves "<<a.curves.size()<<"/"<<b.curves.size()
				<<", spans "<<a.unresolved_edge_face_spans.size()<<"/"
				<<b.unresolved_edge_face_spans.size()<<", overlaps "
				<<a.unresolved_partial_overlaps.size()<<"/"
				<<b.unresolved_partial_overlaps.size()<<", UV failures "
				<<a.unresolved_uv_projections.size()<<"/"
				<<b.unresolved_uv_projections.size();
			return fail(detail.str());
		}
		for(std::size_t curve_id=0;curve_id<a.curves.size();++curve_id)
		{
			const auto& x=a.curves[curve_id];const auto& y=b.curves[curve_id];
			if(x.id!=y.id||x.source_edge_id!=y.source_edge_id
				||x.source_edge_ids!=y.source_edge_ids||x.fan_degree!=y.fan_degree
				||x.uses.size()!=y.uses.size())
			{
				std::ostringstream detail;detail<<"curve "<<curve_id<<" identity/fan mismatch";
				return fail(detail.str());
			}
			for(std::size_t use=0;use<x.uses.size();++use)
				if(x.uses[use].source_face_id!=y.uses[use].source_face_id
					||x.uses[use].kind!=y.uses[use].kind
					||x.uses[use].sector_count!=y.uses[use].sector_count)
				{
					std::ostringstream detail;detail<<"curve "<<curve_id<<" use "<<use
						<<" mismatch";return fail(detail.str());
				}
		}
		for(std::size_t span=0;span<a.unresolved_edge_face_spans.size();++span)
			if(unresolved_span_signature(a.unresolved_edge_face_spans[span])
				!=unresolved_span_signature(b.unresolved_edge_face_spans[span]))
			{
				const auto& x=a.unresolved_edge_face_spans[span];
				const auto& y=b.unresolved_edge_face_spans[span];
				std::ostringstream detail;detail<<std::setprecision(17)
					<<"span "<<span<<" mismatch: edge/face "<<x.source_edge_id<<"/"
					<<x.target_face_id<<" vs "<<y.source_edge_id<<"/"<<y.target_face_id
					<<", range ["<<x.source_parameter_begin<<", "<<x.source_parameter_end
					<<"] vs ["<<y.source_parameter_begin<<", "<<y.source_parameter_end
					<<"], tolerance "<<x.source_parameter_tolerance<<" vs "
					<<y.source_parameter_tolerance<<", issue "
					<<static_cast<unsigned>(x.issue)<<" vs "<<static_cast<unsigned>(y.issue)
					<<", reason '"<<x.reason<<"' vs '"<<y.reason<<"'";
				return fail(detail.str());
			}
		for(std::size_t overlap=0;overlap<a.unresolved_partial_overlaps.size();++overlap)
		{
			const auto& x=a.unresolved_partial_overlaps[overlap];
			const auto& y=b.unresolved_partial_overlaps[overlap];
			if(std::tie(x.source_edge_a,x.source_edge_b,x.owner_face_a,x.owner_face_b,
				x.edge_a_fully_covered,x.edge_b_fully_covered)
				!=std::tie(y.source_edge_a,y.source_edge_b,y.owner_face_a,y.owner_face_b,
					y.edge_a_fully_covered,y.edge_b_fully_covered))
			{
				std::ostringstream detail;detail<<"partial overlap "<<overlap<<" mismatch";
				return fail(detail.str());
			}
		}
		if(compare_uv_failures)for(std::size_t failure=0;
			failure<a.unresolved_uv_projections.size();++failure)
		{
			const auto& x=a.unresolved_uv_projections[failure];
			const auto& y=b.unresolved_uv_projections[failure];
			if(std::tie(x.curve_id,x.source_face_id,x.kind,x.detail)
				!=std::tie(y.curve_id,y.source_face_id,y.kind,y.detail))
			{
				std::ostringstream detail;detail<<"UV failure "<<failure<<" mismatch";
				return fail(detail.str());
			}
		}
		return true;
	}

	void set_contact_thread_environment(const char* value)
	{
#if defined(_WIN32)
		_putenv_s("PARACFD_STEP_CONTACT_THREADS",value?value:"");
#else
		if(value)setenv("PARACFD_STEP_CONTACT_THREADS",value,1);
		else unsetenv("PARACFD_STEP_CONTACT_THREADS");
#endif
	}
}

int main(int argc, char** argv)
{
	const std::string path = argc > 1 ? argv[1] : "Test-Data/NACA2412_C1_S2p03.step";
	bool quality_only = false;
	bool compare_contact_workers=false;
	std::optional<double> comparison_deflection_mm;
	std::optional<std::size_t> expected_components, expected_warnings, expected_contacts,
		expected_source_edges, expected_shared_nodes_min, expected_partial_overlaps,
		expected_uv_failures,expected_edge_face_spans;
	std::vector<ExpectedEdgeFaceSpan> expected_spans;
	for (int argument = 2; argument < argc; ++argument)
	{
		const std::string option = argv[argument];
		if (option == "--quality-only") quality_only = true;
		else if(option=="--compare-contact-workers")compare_contact_workers=true;
		else if(option=="--compare-deflection-mm"&&argument+1<argc)
			comparison_deflection_mm=std::stod(argv[++argument]);
		else if (option == "--expect-components" && argument + 1 < argc)
			expected_components = static_cast<std::size_t>(std::stoull(argv[++argument]));
		else if (option == "--expect-warnings" && argument + 1 < argc)
			expected_warnings = static_cast<std::size_t>(std::stoull(argv[++argument]));
		else if (option == "--expect-contacts" && argument + 1 < argc)
			expected_contacts = static_cast<std::size_t>(std::stoull(argv[++argument]));
		else if (option == "--expect-source-edges" && argument + 1 < argc)
			expected_source_edges = static_cast<std::size_t>(std::stoull(argv[++argument]));
		else if (option == "--expect-shared-nodes-min" && argument + 1 < argc)
			expected_shared_nodes_min = static_cast<std::size_t>(std::stoull(argv[++argument]));
		else if (option == "--expect-partial-overlaps" && argument + 1 < argc)
			expected_partial_overlaps = static_cast<std::size_t>(std::stoull(argv[++argument]));
		else if (option == "--expect-uv-failures" && argument + 1 < argc)
			expected_uv_failures = static_cast<std::size_t>(std::stoull(argv[++argument]));
		else if(option=="--expect-edge-face-spans"&&argument+1<argc)
			expected_edge_face_spans=static_cast<std::size_t>(std::stoull(argv[++argument]));
		else if(option=="--expect-edge-face-span"&&argument+7<argc)
		{
			ExpectedEdgeFaceSpan span;
			span.edge=static_cast<std::uint32_t>(std::stoul(argv[++argument]));
			span.face=static_cast<std::uint32_t>(std::stoul(argv[++argument]));
			span.begin=std::stod(argv[++argument]);span.end=std::stod(argv[++argument]);
			const auto kind=parse_span_kind(argv[++argument]);
			span.sectors=static_cast<std::uint8_t>(std::stoul(argv[++argument]));
			const auto issue=parse_span_issue(argv[++argument]);
			if(!kind||!issue)
			{
				std::fprintf(stderr,"invalid --expect-edge-face-span kind or issue\n");
				return 2;
			}
			span.kind=*kind;span.issue=*issue;expected_spans.push_back(span);
		}
		else
		{
			std::fprintf(stderr, "usage: step_import_provenance_probe [file.step] "
				"[--quality-only] [--expect-components N] [--expect-warnings N] "
				"[--expect-contacts N] [--expect-source-edges N] "
				"[--expect-shared-nodes-min N] [--expect-partial-overlaps N] "
				"[--expect-uv-failures N] [--expect-edge-face-spans N] "
				"[--compare-deflection-mm MM] "
				"[--compare-contact-workers] "
				"[--expect-edge-face-span EDGE FACE BEGIN END KIND SECTORS ISSUE]\n");
			return 2;
		}
	}
	StepCadContactGraph default_contact_graph;
	check(!default_contact_graph.conforming_ready(),
		"a default/unbuilt CAD contact graph is never certified ready");
	GeometryQualityReport default_quality;
	check(default_quality.flow_eligibility()==GeometryFlowEligibility::ready_with_warnings,
		"an unperformed connectivity audit is never reported as clean");
	StepCadContactGraph completed_contact_graph;
	completed_contact_graph.audit_status=StepCadContactAuditStatus::complete;
	check(completed_contact_graph.conforming_ready(),
		"a completed clean CAD contact audit is conforming-ready");
	completed_contact_graph.audit_failure_reason="stale manufactured failure";
	check(!completed_contact_graph.conforming_ready(),
		"a retained contact-audit failure reason prevents strict readiness");
	completed_contact_graph.audit_failure_reason.clear();
	StepCadContactGraph uncertain_bridge;
	uncertain_bridge.audit_status=StepCadContactAuditStatus::complete;
	StepCadContactGraph::UnresolvedEdgeFaceSpan bridge_span;
	bridge_span.source_edge_id=3;bridge_span.target_face_id=7;
	bridge_span.source_parameter_begin=0;bridge_span.source_parameter_end=1;
	bridge_span.endpoint_begin_m={0,0,0};bridge_span.endpoint_end_m={1,0,0};
	bridge_span.issue=StepEdgeFaceSpanIssue::certification_failure;
	bridge_span.reason="manufactured contact bridge could not be certified";
	uncertain_bridge.unresolved_edge_face_spans.push_back(bridge_span);
	GeometryQualityReport uncertain_components;
	uncertain_components.connectivity_status=GeometryConnectivityStatus::unknown;
	GeometryIssue apparent_island;
	apparent_island.kind=GeometryIssueKind::disconnected_fabric_component;
	apparent_island.source_triangle_ids={4,5};
	uncertain_components.issues.push_back(apparent_island);
	check(!uncertain_bridge.conforming_ready()
		&&default_excluded_triangle_ids(uncertain_components).empty(),
		"a failed certificate bridging apparent components blocks readiness and automatic exclusion");
	StepCadContactGraph split_source_contract;
	split_source_contract.curves.push_back(manufactured_contact_atom(2, 7,
		{ 0.0, 0.5 }, { { { 0.0, 0.0, 0.0 } }, { { 0.5, 0.0, 0.0 } } }, { 0, 1 }));
	split_source_contract.curves.push_back(manufactured_contact_atom(5, 7,
		{ 0.5, 1.0 }, { { { 0.5, 0.0, 0.0 } }, { { 1.0, 0.0, 0.0 } } }, { 1, 2 }));
	check(contact_graph_is_structurally_valid(split_source_contract),
		"stable public atom IDs may have fan-1 gaps and one source may span several atoms");
	StepCadContactGraph closed_atom_contract;
	closed_atom_contract.curves.push_back(manufactured_contact_atom(0, 8,
		{ 0.0, 0.3, 0.7, 1.0 },
		{ { { 0.0, 0.0, 0.0 } }, { { 1.0, 0.0, 0.0 } },
			{ { 0.0, 1.0, 0.0 } }, { { 0.0, 0.0, 0.0 } } },
		{ 0, 1, 2, 0 }));
	closed_atom_contract.curves.front().uses.resize(1);
	closed_atom_contract.curves.front().uses.front().kind = StepContactUseKind::face_interior;
	closed_atom_contract.curves.front().uses.front().sector_count = 2;
	check(contact_graph_is_structurally_valid(closed_atom_contract),
		"a closed/seam atom may repeat its first topology node and contribute two sectors");
	std::optional<std::string> saved_contact_threads;
	if(const char* configured=std::getenv("PARACFD_STEP_CONTACT_THREADS"))
		saved_contact_threads=std::string(configured);
	if(compare_contact_workers)set_contact_thread_environment("8");
	std::string preview_error;
	const StepGeometry preview_geometry=load_step_geometry_preview(path,2.0,&preview_error);
	check(!preview_geometry.mesh.empty(),"interactive STEP preview imports without an exact audit");
	check(preview_geometry.contacts.audit_status==StepCadContactAuditStatus::not_performed
		&&!preview_geometry.contacts.conforming_ready(),
		"interactive STEP preview is never presented as exact-contact certified");
	check(preview_geometry.quality.connectivity_status==GeometryConnectivityStatus::unknown
		&&preview_geometry.quality.flow_eligibility()==GeometryFlowEligibility::ready_with_warnings,
		"interactive STEP preview remains explicitly approximate and warning-eligible");
	check(default_excluded_triangle_ids(preview_geometry.quality).empty(),
		"interactive STEP preview never invents automatic exclusions before the detailed check");
	std::string error;
	const StepGeometry coarse_geometry = load_step_geometry(path, 2.0, &error);
	const TriMesh& coarse = coarse_geometry.mesh;
	check(!coarse.empty(), "tracked STEP fixture imports");
	if (coarse.empty())
	{
		if(compare_contact_workers)set_contact_thread_environment(
			saved_contact_threads?saved_contact_threads->c_str():nullptr);
		std::fprintf(stderr, "%s\n", error.c_str());
		return 1;
	}
	check(coarse_geometry.source_faces.size() == 1 + *std::max_element(
		coarse.source_face_ids.begin(), coarse.source_face_ids.end()),
		"OCC-free source-face metadata parallels deterministic face IDs");
	check(contact_graph_is_structurally_valid(coarse_geometry.contacts),
		"exact CAD contact graph has common samples, per-face UV chains, and correct fans");
	check(mesh_contact_labels_match_graph(coarse, coarse_geometry.contacts),
		"conforming mesh contact segments match atom IDs, fan incidence, and topology positions");
	const ContactGraphCounts contact_counts = contact_graph_counts(coarse_geometry.contacts);
	if (expected_components) check(coarse_geometry.quality.exact_connected_component_count
		== *expected_components, "fixture has the expected exact CAD component count");
	if (expected_warnings) check(coarse_geometry.quality.warning_count() == *expected_warnings,
		"fixture has the expected geometry warning count");
	if (expected_contacts) check(coarse_geometry.contacts.curves.size() == *expected_contacts,
		"fixture has the expected deterministic contact-curve count");
	if (expected_source_edges) check(contact_counts.source_edges == *expected_source_edges,
		"canonical curves retain every expected source-edge provenance record");
	if (expected_shared_nodes_min) check(contact_counts.shared_topology_nodes
		>= *expected_shared_nodes_min,
		"exact curve junctions share topology IDs across canonical physical curves");
	if (expected_partial_overlaps)
		check(coarse_geometry.contacts.unresolved_partial_overlaps.size()
			== *expected_partial_overlaps,
			"fixture has the expected exact partial-contact overlap count");
	if (expected_uv_failures)
		check(coarse_geometry.contacts.unresolved_uv_projections.size()
			== *expected_uv_failures,
			"fixture has the expected unresolved contact UV count");
	if(expected_edge_face_spans)
		check(coarse_geometry.contacts.unresolved_edge_face_spans.size()
			==*expected_edge_face_spans,"fixture has the expected unresolved edge/face span count");
	for(const ExpectedEdgeFaceSpan& expected:expected_spans)
		check(std::any_of(coarse_geometry.contacts.unresolved_edge_face_spans.begin(),
			coarse_geometry.contacts.unresolved_edge_face_spans.end(),[&](const auto& actual)
			{return expected_span_matches(expected,actual);}),
			"fixture retains the expected exact unresolved edge/face parameter span");
	if(compare_contact_workers)
	{
		set_contact_thread_environment("1");
		std::string serial_error;
		const StepGeometry serial=load_step_geometry(path,2.0,&serial_error);
		set_contact_thread_environment(saved_contact_threads
			?saved_contact_threads->c_str():nullptr);
		check(!serial.mesh.empty(),"single-worker STEP contact audit succeeds");
		check(!serial.mesh.empty()
			&&exact_contact_topology_equal(coarse_geometry.contacts,serial.contacts)
			&&coarse_geometry.quality.exact_connected_component_count
				==serial.quality.exact_connected_component_count,
			"single- and multi-worker exact contact audits are deterministic");
	}
	if(comparison_deflection_mm)
	{
		std::string comparison_error;
		const StepGeometry comparison=load_step_geometry(path,*comparison_deflection_mm,
			&comparison_error);
		check(!comparison.mesh.empty(),"comparison-deflection STEP import succeeds");
		std::string topology_mismatch;
		const bool topology_equal=!comparison.mesh.empty()
			&&exact_contact_topology_equal(coarse_geometry.contacts,comparison.contacts,
				&topology_mismatch,false);
		check(topology_equal
			&&coarse_geometry.quality.exact_connected_component_count
				==comparison.quality.exact_connected_component_count,
			"exact contact certificates and unresolved spans are tessellation-deflection invariant");
		if(!topology_equal)std::fprintf(stderr,"deflection comparison detail: %s\n",
			topology_mismatch.c_str());
	}
	if (!expected_partial_overlaps && !expected_uv_failures&&!expected_edge_face_spans)
		check(coarse_geometry.contacts.conforming_ready(),
			"contact graph has no unresolved edge/face spans, partial overlaps, or UV projections");
	if (quality_only)
	{
		std::printf("STEP quality: audit=%u ready=%s connectivity=%u components=%zu contacts=%zu source-edges=%zu samples=%zu "
			"shared-nodes=%zu edge-face-spans=%zu partial-overlaps=%zu UV-failures=%zu trim-UV-fallbacks=%u "
			"largest-area=%.12g m^2 warnings=%zu blockers=%zu\n",
			static_cast<unsigned>(coarse_geometry.contacts.audit_status),
			coarse_geometry.contacts.conforming_ready()?"yes":"no",
			static_cast<unsigned>(coarse_geometry.quality.connectivity_status),
			coarse_geometry.quality.exact_connected_component_count,
			coarse_geometry.contacts.curves.size(),
			contact_counts.source_edges, contact_counts.samples,
			contact_counts.shared_topology_nodes,
			coarse_geometry.contacts.unresolved_edge_face_spans.size(),
			coarse_geometry.contacts.unresolved_partial_overlaps.size(),
			coarse_geometry.contacts.unresolved_uv_projections.size(),
			coarse_geometry.contacts.trim_surface_fallback_count,
			coarse_geometry.quality.largest_component_area_m2,
			coarse_geometry.quality.warning_count(), coarse_geometry.quality.blocker_count());
		if(!coarse_geometry.contacts.audit_failure_reason.empty())
			std::printf("  contact audit failure: %s\n",
				coarse_geometry.contacts.audit_failure_reason.c_str());
		for(const auto& span:coarse_geometry.contacts.unresolved_edge_face_spans)
			std::printf("  unresolved edge/face span edge %u face %u raw [%.17g, %.17g] "
				"tol %.3g kind=%u sectors=%u occurrences=%u ids=[%u,%u] issue=%u "
				"p0=[%.9g %.9g %.9g] "
				"p1=[%.9g %.9g %.9g]: %s\n",span.source_edge_id,span.target_face_id,
				span.source_parameter_begin,span.source_parameter_end,
				span.source_parameter_tolerance,static_cast<unsigned>(span.kind),
				static_cast<unsigned>(span.target_sector_count),
				static_cast<unsigned>(span.target_boundary_occurrence_count),
				span.target_boundary_occurrence_ids[0],span.target_boundary_occurrence_ids[1],
				static_cast<unsigned>(span.issue),
				span.endpoint_begin_m[0],span.endpoint_begin_m[1],span.endpoint_begin_m[2],
				span.endpoint_end_m[0],span.endpoint_end_m[1],span.endpoint_end_m[2],
				span.reason.c_str());
		for (const StepCadContactGraph::PartialOverlap& overlap :
			coarse_geometry.contacts.unresolved_partial_overlaps)
			std::printf("  unresolved partial contact edges %u / %u, owner faces %u / %u "
				"(full coverage %s / %s)\n",
				overlap.source_edge_a, overlap.source_edge_b,
				overlap.owner_face_a, overlap.owner_face_b,
				overlap.edge_a_fully_covered ? "yes" : "no",
				overlap.edge_b_fully_covered ? "yes" : "no");
		for (const StepCadContactGraph::UvProjectionFailure& failure :
			coarse_geometry.contacts.unresolved_uv_projections)
			std::printf("  unresolved UV contact curve %llu, face %u, kind=%s: %s\n",
				static_cast<unsigned long long>(failure.curve_id), failure.source_face_id,
				failure.kind == StepContactUseKind::trim_boundary ? "boundary" : "interior",
				failure.detail.c_str());
		for (const GeometryIssue& issue : coarse_geometry.quality.issues)
		{
			std::printf("  issue=%llu kind=%u faces=%zu triangles=%zu area=%.12g m^2 boundary-lines=%zu "
				"boundary-length=%.12g m\n", static_cast<unsigned long long>(issue.id),
				static_cast<unsigned>(issue.kind), issue.source_face_ids.size(),
				issue.source_triangle_ids.size(), issue.surface_area_m2,
				issue.boundary_polylines.size(), issue.boundary_length_m);
			if (!issue.source_representation_names.empty())
			{
				std::printf("    representation names:");
				for (const std::string& name : issue.source_representation_names)
					std::printf(" [%s]", name.c_str());
				std::printf("\n");
			}
			std::printf("    faces:");
			for (std::uint32_t face : issue.source_face_ids) std::printf(" %u", face);
			std::printf("\n    triangles:");
			for (std::uint32_t triangle : issue.source_triangle_ids) std::printf(" %u", triangle);
			std::printf("\n");
		}
		return failures == 0 ? 0 : 1;
	}
	check(coarse_geometry.quality.exact_connected_component_count == 1
		&& coarse_geometry.quality.issues.empty()
		&& coarse_geometry.quality.flow_eligibility() == GeometryFlowEligibility::clean,
		"closed STEP fixture is one exact CAD-connected fabric component");
	GeometryQualityReport warning_contract;
	GeometryIssue warning;
	warning.kind = GeometryIssueKind::disconnected_fabric_component;
	warning.severity = GeometryIssueSeverity::warning_run;
	warning_contract.issues.push_back(warning);
	check(warning_contract.flow_eligibility() == GeometryFlowEligibility::ready_with_warnings
		&& warning_contract.warning_count() == 1 && warning_contract.blocker_count() == 0,
		"representable geometry warnings retain flow eligibility");

	check(coarse.has_fp64_positions() && coarse.positions_fp64.size() == coarse.positions.size(),
		"STEP import retains one FP64 coordinate beside every FP32 coordinate");
	bool rounded_copy = true;
	for (std::size_t coordinate = 0; coordinate < coarse.positions.size(); ++coordinate)
		rounded_copy = rounded_copy
			&& coarse.positions[coordinate] == static_cast<float>(coarse.positions_fp64[coordinate]);
	check(rounded_copy, "FP32 display positions are direct rounded copies of OCCT FP64 positions");

	const TopologySummary coarse_topology = summarize(coarse);
	check(coarse_topology.structurally_valid, "CAD edge IDs and incident-face counts are structurally consistent");
	check(coarse_topology.boundary_segments > 0
		&& coarse_topology.labelled_boundary_segments == coarse_topology.boundary_segments,
		"every face-local tessellation boundary segment maps to a CAD edge");
	check(coarse_topology.labelled_interior_segments == 0,
		"internal tessellation half-edges retain the unavailable sentinel");
	check(std::all_of(coarse_topology.periodic_seam_segments.begin(),
		coarse_topology.periodic_seam_segments.end(),
		[](const auto& entry) { return entry.second >= 2; }),
		"both periodic-seam polygon chains carry one CAD edge ID with incidence at least two");
	check(std::any_of(coarse_topology.edge_incidence.begin(), coarse_topology.edge_incidence.end(),
		[](const auto& entry) { return entry.second == 2; }),
		"closed fixture exposes CAD edges incident on two source faces");

	// Face-local mesh copies remain separate storage, but conforming contact samples carry
	// one exact topology ID (and, above, a bit-identical FP64 position) across all uses.
	std::map<std::uint32_t, std::set<std::uint32_t>> topology_faces;
	for (std::size_t triangle = 0; triangle < coarse.triangle_count(); ++triangle)
		for (unsigned corner = 0; corner < 3; ++corner)
		{
			const std::uint32_t vertex = coarse.indices[3 * triangle + corner];
			topology_faces[coarse.topology_vertex_id(vertex)].insert(
				coarse.source_face_ids[triangle]);
		}
	check(std::any_of(topology_faces.begin(), topology_faces.end(),
		[](const auto& entry) { return entry.second.size() > 1; }),
		"face-local vertex copies share exact topology IDs across contact uses");

	const TriMesh fine = load_step_mesh(path, 0.5, &error);
	const TopologySummary fine_topology = summarize(fine);
	check(!fine.empty() && fine_topology.structurally_valid,
		"same STEP topology imports at a second tessellation tolerance");
	check(fine_topology.face_edge_uses == coarse_topology.face_edge_uses
		&& fine_topology.edge_incidence == coarse_topology.edge_incidence
		&& fine_topology.edge_tolerance_m == coarse_topology.edge_tolerance_m
		&& fine_topology.periodic_seam_segments.size() == coarse_topology.periodic_seam_segments.size()
		&& std::equal(fine_topology.periodic_seam_segments.begin(),
			fine_topology.periodic_seam_segments.end(), coarse_topology.periodic_seam_segments.begin(),
			[](const auto& a, const auto& b) { return a.first == b.first; }),
		"source face/edge IDs and incidence remain stable across re-tessellation");

	TriMesh programmatic;
	programmatic.positions = { 1.25f, -2.5f, 3.75f };
	const auto promoted = programmatic.vertex_position_double(0);
	check(!programmatic.has_fp64_positions() && promoted[0] == 1.25 && promoted[1] == -2.5
		&& promoted[2] == 3.75,
		"programmatic mesh accessor promotes FP32 when the exact sidecar is absent");

	TriMesh exact_programmatic = programmatic;
	exact_programmatic.positions_fp64 = { 1.25 + 1e-11, -2.5 - 2e-11, 3.75 + 3e-11 };
	ModelPlacement placement;
	placement.tx = 10.0;
	placement.ty = -20.0;
	placement.tz = 30.0;
	const TriMesh placed = placed_mesh(exact_programmatic, placement);
	check(placed.has_fp64_positions()
		&& std::abs(placed.positions_fp64[0] - (11.25 + 1e-11)) < 1e-14
		&& std::abs(placed.positions_fp64[1] - (-22.5 - 2e-11)) < 1e-14
		&& std::abs(placed.positions_fp64[2] - (33.75 + 3e-11)) < 1e-14,
		"model placement transforms retained CPU geometry without an FP32 round trip");

	TriMesh tolerance_mesh;
	tolerance_mesh.positions = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	tolerance_mesh.positions_fp64 = { 0, 0, 1e-11, 1, 0, 1e-11, 0, 1, 1e-11 };
	tolerance_mesh.indices = { 0, 1, 2 };
	tolerance_mesh.triangle_cad_edge_provenance_states = {
		static_cast<std::uint8_t>(CadEdgeProvenanceState::known),
		static_cast<std::uint8_t>(CadEdgeProvenanceState::none),
		static_cast<std::uint8_t>(CadEdgeProvenanceState::known) };
	tolerance_mesh.triangle_cad_edge_ids = { 7, TriMesh::kNoCadEdgeId, 8 };
	tolerance_mesh.triangle_cad_edge_incident_face_counts = { 2, 0, 1 };
	tolerance_mesh.triangle_cad_edge_tolerances = { 1e-6, 0.0, 2e-6 };
	tolerance_mesh.triangle_cad_edge_is_periodic_seam = { 0, 0, 0 };
	ModelPlacement anisotropic;
	anisotropic.m[0] = 2.0;
	anisotropic.m[4] = 3.0;
	anisotropic.m[8] = 0.5;
	const TriMesh scaled = placed_mesh(tolerance_mesh, anisotropic);
	check(scaled.cad_edge_provenance_state(0, 0) == CadEdgeProvenanceState::known
		&& scaled.cad_edge_provenance_state(0, 1) == CadEdgeProvenanceState::none
		&& scaled.cad_edge_tolerance(0, 0) >= 3e-6
		&& scaled.cad_edge_tolerance(0, 0) < 3e-6 * (1.0 + 1e-12)
		&& scaled.cad_edge_tolerance(0, 1) == 0.0
		&& scaled.cad_edge_tolerance(0, 2) >= 6e-6,
		"affine placement scales CAD tolerances by its largest singular value");

	const TriMesh clipped = clip_mesh_to_axis_slab(tolerance_mesh, 0, 0.25, 0.75);
	bool retained_source_edge = false, introduced_sentinel = false;
	for (std::size_t triangle = 0; triangle < clipped.triangle_count(); ++triangle)
		for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
		{
			const std::uint32_t id = clipped.cad_edge_id(triangle, half_edge);
			if (id == 7)
				retained_source_edge = retained_source_edge
					|| (clipped.cad_edge_incident_face_count(triangle, half_edge) == 2
						&& clipped.cad_edge_tolerance(triangle, half_edge) == 1e-6);
			else if (id == TriMesh::kNoCadEdgeId)
				introduced_sentinel = introduced_sentinel
					|| (clipped.cad_edge_incident_face_count(triangle, half_edge) == 0
						&& clipped.cad_edge_tolerance(triangle, half_edge) == 0.0);
		}
	check(clipped.has_fp64_positions() && clipped.has_cad_edge_provenance()
		&& retained_source_edge && introduced_sentinel,
		"slab clipping retains exact/source-edge metadata and marks synthetic edges unavailable");

	TriMesh unknown_boundary_mesh = tolerance_mesh;
	unknown_boundary_mesh.triangle_cad_edge_provenance_states[0] =
		static_cast<std::uint8_t>(CadEdgeProvenanceState::unknown_boundary);
	unknown_boundary_mesh.triangle_cad_edge_ids[0] = TriMesh::kNoCadEdgeId;
	unknown_boundary_mesh.triangle_cad_edge_incident_face_counts[0] = 0;
	unknown_boundary_mesh.triangle_cad_edge_tolerances[0] = 0.0;
	const TriMesh clipped_unknown = clip_mesh_to_axis_slab(unknown_boundary_mesh, 0, 0.25, 0.75);
	bool retained_unknown = false;
	for (std::size_t triangle = 0; triangle < clipped_unknown.triangle_count(); ++triangle)
		for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
			retained_unknown = retained_unknown || clipped_unknown.cad_edge_provenance_state(
				triangle, half_edge) == CadEdgeProvenanceState::unknown_boundary;
	check(retained_unknown,
		"slab clipping preserves an inherited unresolved CAD boundary distinctly from synthetic edges");

	std::printf("STEP provenance: vertices=%zu triangles=%zu CAD edges=%zu face-edge uses=%zu "
		"UV-failures=%zu trim-UV-fallbacks=%u\n",
		coarse.vertex_count(), coarse.triangle_count(), coarse_topology.edge_incidence.size(),
		coarse_topology.face_edge_uses.size(),
		coarse_geometry.contacts.unresolved_uv_projections.size(),
		coarse_geometry.contacts.trim_surface_fallback_count);
	return failures == 0 ? 0 : 1;
}
