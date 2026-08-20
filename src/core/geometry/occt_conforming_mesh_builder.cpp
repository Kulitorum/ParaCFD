#include "core/geometry/occt_conforming_mesh_builder.h"

#include "core/geometry/step_face_uv_projection.h"

#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

namespace paracfd::core
{
	namespace
	{
		using FaceId = std::uint32_t;
		using AtomId = std::uint64_t;
		using FaceAtom = std::pair<FaceId, AtomId>;
		using FaceCad = std::pair<FaceId, std::uint32_t>;
		using FaceOccurrence = std::pair<FaceId, std::uint32_t>;

		struct FaceStage
		{
			const OcctContactTopologyFace* source = nullptr;
			std::vector<OcctFaceConstraintPolyline> cad_constraints;
			std::vector<OcctFaceConstraintPolyline> atom_constraints;
			std::uint64_t next_boundary_branch_id = 0;
			std::set<AtomId> expected_atoms;
			std::set<std::pair<std::uint32_t, AtomId>> expected_boundary_pairs;
			std::map<std::tuple<std::uint32_t, std::uint64_t, std::uint32_t>,
				std::uint32_t>
				source_node_by_occurrence_topology;
		};

		struct PreparedChain
		{
			std::vector<OcctFaceConstraintSample> samples;
		};

		bool finite(double value)
		{
			return std::isfinite(value);
		}

		bool finite(const std::array<double, 3>& value)
		{
			return finite(value[0]) && finite(value[1]) && finite(value[2]);
		}

		bool valid_orientation(std::int8_t orientation)
		{
			return orientation == static_cast<std::int8_t>(TopAbs_FORWARD)
				|| orientation == static_cast<std::int8_t>(TopAbs_REVERSED);
		}

		double parameter_roundoff(double first, double last, double parameter)
		{
			const double scale = std::max({std::abs(first), std::abs(last),
				std::abs(parameter), std::abs(last - first),
				std::numeric_limits<double>::min()});
			return 64.0 * std::numeric_limits<double>::epsilon() * scale;
		}

		bool same_position_bits(const std::array<double, 3>& a,
			const std::array<double, 3>& b)
		{
			for (std::size_t axis = 0; axis < 3; ++axis)
				if (std::bit_cast<std::uint64_t>(a[axis])
					!= std::bit_cast<std::uint64_t>(b[axis])) return false;
			return true;
		}

		bool same_uv_bits(const UvPoint& a, const UvPoint& b)
		{
			return std::bit_cast<std::uint64_t>(a.u)
					== std::bit_cast<std::uint64_t>(b.u)
				&& std::bit_cast<std::uint64_t>(a.v)
					== std::bit_cast<std::uint64_t>(b.v);
		}

		bool split_chart_aliases(std::vector<OcctFaceConstraintSample> samples,
			std::vector<PreparedChain>& output, std::string& error)
		{
			output.clear();
			if (samples.size() < 2)
			{
				error = "boundary constraint has fewer than two samples";
				return false;
			}
			std::size_t begin = 0;
			while (begin + 1 < samples.size())
			{
				std::map<std::uint32_t, UvPoint> claimed;
				std::size_t end = begin;
				for (; end < samples.size(); ++end)
				{
					const auto& sample = samples[end];
					if (sample.topology_id == UvVertex::no_topology_id) continue;
					const auto found = claimed.find(sample.topology_id);
					if (found != claimed.end() && !same_uv_bits(found->second, sample.uv))
						break;
					claimed.emplace(sample.topology_id, sample.uv);
				}
				if (end == begin + 1)
				{
					error = "periodic closed boundary constraint needs an exact interior sample"
						" between aliased endpoints";
					return false;
				}
				const std::size_t inclusive_end = end == samples.size()
					? samples.size() - 1 : end - 1;
				PreparedChain chain;
				chain.samples.assign(samples.begin() + static_cast<std::ptrdiff_t>(begin),
					samples.begin() + static_cast<std::ptrdiff_t>(inclusive_end + 1));
				output.push_back(std::move(chain));
				if (end == samples.size()) break;
				begin = inclusive_end;
			}
			return true;
		}

		bool boundary_uv(const FaceStage& stage,
			const OcctContactTopologyBoundaryOccurrence& occurrence,
			const OcctContactAtomBoundaryUse& use,
			const OcctContactAtom& atom, double exact_operation_tolerance,
			std::vector<OcctFaceConstraintSample>& samples, std::string& error)
		{
			const OcctContactTopologyFace& face = *stage.source;
			samples.clear();
			if (use.sample_target_parameters.size() != atom.samples.size())
			{
				error = "boundary occurrence parameter count does not match its atom samples";
				return false;
			}
			double first = 0.0, last = 0.0;
			const Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(
				occurrence.edge, face.face, first, last);
			if (pcurve.IsNull() || !finite(first) || !finite(last) || !(last > first))
			{
				error = "boundary occurrence has no finite exact face pcurve";
				return false;
			}
			samples.reserve(atom.samples.size());
			std::vector<std::size_t> order(atom.samples.size());
			std::iota(order.begin(), order.end(), std::size_t{0});
			std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
			{
				return std::tie(use.sample_target_parameters[a], a)
					< std::tie(use.sample_target_parameters[b], b);
			});
			// Atom sample order is physical/canonical and may wrap a closed edge as
			// 0, last-delta, ..., first+delta, last.  A face-chart polyline must instead
			// follow this occurrence's native pcurve parameter monotonically; otherwise
			// the wrap is emitted as a long segment across the periodic chart and skips
			// the original boundary nodes which are already present later in the atom.
			for (std::size_t ordered_index = 0;
				ordered_index < order.size(); ++ordered_index)
			{
				const std::size_t index = order[ordered_index];
				const double parameter = use.sample_target_parameters[index];
				const double roundoff = parameter_roundoff(first, last, parameter);
				if (!finite(parameter) || parameter < first - roundoff
					|| parameter > last + roundoff)
				{
					std::ostringstream message;
					message << "boundary occurrence sample " << index
						<< " parameter is outside its exact pcurve range";
					error = message.str();
					return false;
				}
				const gp_Pnt2d uv = pcurve->Value(parameter);
				if (!finite(uv.X()) || !finite(uv.Y()))
				{
					error = "boundary occurrence pcurve produced non-finite UV coordinates";
					return false;
				}
				const auto& source = atom.samples[index];
				const auto claim = stage.source_node_by_occurrence_topology.find(
					{occurrence.occurrence_id, atom.id, static_cast<std::uint32_t>(index)});
				// The Boolean-operation certificate belongs to this atom/face use only.
				// Carry it with this constraint sample so it cannot widen source nodes or
				// another face incident on the same canonical topology vertex.
				samples.push_back({{uv.X(), uv.Y()}, source.topology_id,
					source.canonical_world_position,
					source.canonical_position_tolerance + exact_operation_tolerance,
					claim == stage.source_node_by_occurrence_topology.end()
						? OcctFaceConstraintSample::no_source_node_index : claim->second,
					occurrence.occurrence_id, atom.id,
					static_cast<std::uint32_t>(index)});
			}
			return true;
		}

		bool interior_uv(const FaceStage& stage,
			const OcctContactAtom& atom, double exact_operation_tolerance,
			std::vector<OcctFaceConstraintSample>& samples, std::string& error)
		{
			const OcctContactTopologyFace& face = *stage.source;
			std::vector<gp_Pnt> points;
			std::vector<double> tolerances;
			points.reserve(atom.samples.size());
			tolerances.reserve(atom.samples.size());
			for (const auto& sample : atom.samples)
			{
				points.emplace_back(sample.canonical_world_position[0],
					sample.canonical_world_position[1],
					sample.canonical_world_position[2]);
				// This is an acceptance certificate for an already-exact contact span,
				// not a projection/search or proximity tolerance. Keep it isolated to the
				// current face-use samples.
				tolerances.push_back(sample.canonical_position_tolerance
					+ exact_operation_tolerance);
			}
			std::vector<std::array<double, 2>> uv;
			if (!detail::project_trim_valid_face_uv(points, face.face, {}, false,
				tolerances, uv, error))
			{
				error = "interior atom projection failed: " + error;
				return false;
			}
			if (uv.size() != atom.samples.size())
			{
				error = "interior atom projection returned the wrong sample count";
				return false;
			}
			samples.clear();
			samples.reserve(atom.samples.size());
			for (std::size_t index = 0; index < atom.samples.size(); ++index)
			{
				const auto claim = stage.source_node_by_occurrence_topology.find(
					{OcctContactBoundaryOccurrence::no_occurrence_id, atom.id,
						static_cast<std::uint32_t>(index)});
				samples.push_back({{uv[index][0], uv[index][1]},
					atom.samples[index].topology_id,
					atom.samples[index].canonical_world_position,
					atom.samples[index].canonical_position_tolerance
						+ exact_operation_tolerance,
					claim == stage.source_node_by_occurrence_topology.end()
						? OcctFaceConstraintSample::no_source_node_index : claim->second,
					OcctContactBoundaryOccurrence::no_occurrence_id, atom.id,
					static_cast<std::uint32_t>(index)});
			}
			return true;
		}

		bool append_boundary_constraints(FaceStage& stage,
			const OcctContactTopologyBoundaryOccurrence& occurrence,
			const OcctContactAtomBoundaryUse& use, const OcctContactAtom& atom,
			double exact_operation_tolerance, std::string& error)
		{
			std::vector<OcctFaceConstraintSample> samples;
			if (!boundary_uv(stage, occurrence, use, atom,
				exact_operation_tolerance, samples, error))
				return false;
			std::vector<PreparedChain> pieces;
			if (!split_chart_aliases(std::move(samples), pieces, error)) return false;
			const auto atom_tag = occt_contact_atom_constraint_id(atom.id);
			if (!atom_tag)
			{
				error = "contact atom ID exceeds the constraint-tag payload";
				return false;
			}
			const std::uint64_t cad_tag =
				occt_cad_edge_constraint_id(occurrence.source_edge_id);
			for (const PreparedChain& piece : pieces)
			{
				if (stage.next_boundary_branch_id
					== OcctFaceConstraintPolyline::no_chart_branch_id)
				{
					error = "boundary constraint chart-branch identity space is exhausted";
					return false;
				}
				// This provisional identity exists only to pair the CAD-edge and
				// contact-atom copies of this exact chain piece.  Occurrence IDs repeat
				// across all atoms spanning one CAD edge, and piece indices restart for
				// every atom; deriving the identity from either would cross-attach
				// unrelated source mesh nodes before physical chart normalization.
				const std::uint64_t branch_id = stage.next_boundary_branch_id++;
				stage.cad_constraints.push_back({cad_tag, piece.samples, branch_id});
				stage.atom_constraints.push_back({*atom_tag, piece.samples, branch_id});
			}
			stage.expected_atoms.insert(atom.id);
			stage.expected_boundary_pairs.emplace(occurrence.source_edge_id, atom.id);
			return true;
		}

		bool append_interior_constraint(FaceStage& stage,
			const OcctContactAtom& atom, double exact_operation_tolerance,
			std::string& error)
		{
			std::vector<OcctFaceConstraintSample> samples;
			if (!interior_uv(stage, atom, exact_operation_tolerance,
				samples, error)) return false;
			const auto tag = occt_contact_atom_constraint_id(atom.id);
			if (!tag)
			{
				error = "contact atom ID exceeds the constraint-tag payload";
				return false;
			}
			stage.atom_constraints.push_back({*tag, std::move(samples),
				OcctFaceConstraintPolyline::no_chart_branch_id});
			stage.expected_atoms.insert(atom.id);
			return true;
		}

		bool normalize_chart_branches(
			std::vector<OcctFaceConstraintPolyline>& constraints, std::string& error)
		{
			// A periodic physical vertex can occur at more than one UV location. The
			// conformer requires every polyline meeting one such UV occurrence to name
			// the same branch, and the other UV occurrence to name a different branch.
			// First join the deliberate CAD/atom duplicate of each boundary chain, then
			// join chains which meet at one aliased topology vertex. A one-chart-location
			// junction does not join branches: doing so would reconnect the two halves of
			// a closed periodic boundary through their ordinary midpoint.
			std::vector<std::size_t> parent(constraints.size());
			for (std::size_t index = 0; index < parent.size(); ++index) parent[index] = index;
			auto root = [&](std::size_t value)
			{
				std::size_t current = value;
				while (parent[current] != current) current = parent[current];
				while (parent[value] != value)
				{
					const std::size_t next = parent[value];
					parent[value] = current;
					value = next;
				}
				return current;
			};
			auto join = [&](std::size_t a, std::size_t b)
			{
				a = root(a);
				b = root(b);
				if (a != b) parent[std::max(a, b)] = std::min(a, b);
			};

			std::map<std::uint64_t, std::size_t> first_provisional;
			for (std::size_t polyline = 0; polyline < constraints.size(); ++polyline)
			{
				const auto branch = constraints[polyline].chart_branch_id;
				if (branch == OcctFaceConstraintPolyline::no_chart_branch_id) continue;
				const auto [found, inserted] = first_provisional.emplace(branch, polyline);
				if (!inserted) join(found->second, polyline);
			}

			struct PhysicalOccurrence
			{
				std::size_t polyline = 0;
				UvPoint uv;
			};
			std::map<std::uint32_t, std::vector<PhysicalOccurrence>> physical;
			for (std::size_t polyline = 0; polyline < constraints.size(); ++polyline)
			{
				if (constraints[polyline].chart_branch_id
					== OcctFaceConstraintPolyline::no_chart_branch_id) continue;
				for (const auto& sample : constraints[polyline].samples)
					if (sample.topology_id != UvVertex::no_topology_id)
						physical[sample.topology_id].push_back({polyline, sample.uv});
			}
			for (const auto& [topology_id, occurrences] : physical)
			{
				(void)topology_id;
				std::vector<std::vector<std::size_t>> clusters;
				for (std::size_t occurrence = 0; occurrence < occurrences.size(); ++occurrence)
				{
					auto cluster = std::find_if(clusters.begin(), clusters.end(),
						[&](const auto& members)
						{
							return same_uv_bits(occurrences[members.front()].uv,
								occurrences[occurrence].uv);
						});
					if (cluster == clusters.end()) clusters.push_back({occurrence});
					else cluster->push_back(occurrence);
				}
				if (clusters.size() <= 1) continue;
				for (const auto& cluster : clusters)
					for (std::size_t member = 1; member < cluster.size(); ++member)
						join(occurrences[cluster.front()].polyline,
							occurrences[cluster[member]].polyline);
				for (std::size_t a = 0; a < clusters.size(); ++a)
					for (std::size_t b = a + 1; b < clusters.size(); ++b)
						if (root(occurrences[clusters[a].front()].polyline)
							== root(occurrences[clusters[b].front()].polyline))
						{
							error = "one boundary polyline spans incompatible periodic chart aliases";
							return false;
						}
			}

			std::map<std::size_t, std::uint64_t> branch_by_root;
			for (std::size_t polyline = 0; polyline < constraints.size(); ++polyline)
			{
				const auto provisional = constraints[polyline].chart_branch_id;
				if (provisional == OcctFaceConstraintPolyline::no_chart_branch_id) continue;
				const std::size_t component = root(polyline);
				auto [stored, inserted] = branch_by_root.emplace(component, provisional);
				if (!inserted) stored->second = std::min(stored->second, provisional);
			}
			std::set<std::uint64_t> assigned;
			for (const auto& [unused, branch] : branch_by_root)
			{
				(void)unused;
				if (branch == OcctFaceConstraintPolyline::no_chart_branch_id
					|| !assigned.insert(branch).second)
				{
					error = "periodic chart branch identity is not deterministic and unique";
					return false;
				}
			}
			for (std::size_t polyline = 0; polyline < constraints.size(); ++polyline)
				if (constraints[polyline].chart_branch_id
					!= OcctFaceConstraintPolyline::no_chart_branch_id)
					constraints[polyline].chart_branch_id = branch_by_root[root(polyline)];
			return true;
		}

		bool decode_edges(const OcctContactTopology& topology, const FaceStage& stage,
			OcctConformingFaceMesh& face_output, std::string& error)
		{
			std::set<std::uint32_t> known_edges;
			for (const auto& edge : topology.edges) known_edges.insert(edge.source_edge_id);
			std::set<AtomId> observed_atoms;
			std::set<std::pair<std::uint32_t, AtomId>> observed_boundary_pairs;
			face_output.edge_provenance.clear();
			face_output.edge_provenance.reserve(face_output.mesh.edges.size());
			for (const auto& edge : face_output.mesh.edges)
			{
				OcctConformedMeshEdgeProvenance decoded;
				decoded.vertices = edge.vertices;
				decoded.boundary = edge.boundary;
				for (const std::uint64_t tag : edge.constraint_ids)
				{
					if (const auto cad = occt_cad_edge_from_constraint(tag))
					{
						if (decoded.source_cad_edge_id && *decoded.source_cad_edge_id != *cad)
						{
							error = "one conformed edge carries more than one CAD-edge tag";
							return false;
						}
						if (!known_edges.contains(*cad))
						{
							error = "conformed edge carries an unknown CAD-edge tag";
							return false;
						}
						decoded.source_cad_edge_id = *cad;
					}
					else if (const auto atom = occt_contact_atom_from_constraint(tag))
					{
						if (*atom >= topology.atomization.atoms.size()
							|| topology.atomization.atoms[static_cast<std::size_t>(*atom)].id
								!= *atom)
						{
							error = "conformed edge carries an unknown contact-atom tag";
							return false;
						}
						if (decoded.contact_atom_id && *decoded.contact_atom_id != *atom)
						{
							error = "one conformed edge carries more than one contact-atom tag";
							return false;
						}
						decoded.contact_atom_id = *atom;
					}
					else
					{
						error = "conformed edge carries an unknown constraint-tag namespace";
						return false;
					}
				}
				if (edge.boundary && !decoded.source_cad_edge_id)
				{
					error = "a final face boundary edge has no exact CAD-edge tag";
					return false;
				}
				if (decoded.source_cad_edge_id && !edge.boundary)
				{
					error = "an interior conformed edge incorrectly carries a CAD-boundary tag";
					return false;
				}
				if (edge.boundary && !decoded.contact_atom_id)
				{
					error = "a final face boundary edge has no exact contact-atom tag";
					return false;
				}
				if (decoded.contact_atom_id)
				{
					observed_atoms.insert(*decoded.contact_atom_id);
					if (edge.boundary)
					{
						if (!decoded.source_cad_edge_id)
						{
							error = "a boundary contact atom has no CAD-edge lineage";
							return false;
						}
						observed_boundary_pairs.emplace(*decoded.source_cad_edge_id,
							*decoded.contact_atom_id);
					}
				}
				face_output.edge_provenance.push_back(std::move(decoded));
			}
			if (observed_atoms != stage.expected_atoms)
			{
				error = "conformed face contact-atom tag incidence differs from its exact fan";
				return false;
			}
			if (observed_boundary_pairs != stage.expected_boundary_pairs)
			{
				error = "conformed boundary CAD-edge/contact-atom incidence differs from"
					" its exact contact topology";
				return false;
			}
			return true;
		}
	}

	bool build_occt_conforming_mesh(const OcctContactTopology& topology,
		OcctConformingMesh& output, std::string& error,
		OcctFaceConformerOptions options)
	{
		error.clear();
		try
		{
			if (!finite(options.geometric_tolerance)
				|| options.geometric_tolerance < 0.0)
			{
				error = "conforming-mesh geometric tolerance is invalid";
				return false;
			}
			std::map<FaceId, FaceStage> stages;
			std::set<std::uint32_t> edge_ids;
			std::set<FaceOccurrence> expected_occurrences;
			for (const auto& edge : topology.edges)
			{
				if (edge.edge.IsNull() || !edge_ids.insert(edge.source_edge_id).second)
				{
					error = "contact topology has a null or duplicate source CAD edge";
					return false;
				}
			}
			for (const auto& face : topology.faces)
			{
				if (face.face.IsNull() || stages.contains(face.source_face_id))
				{
					error = "contact topology has a null or duplicate source face";
					return false;
				}
				FaceStage stage;
				stage.source = &face;
				stages.emplace(face.source_face_id, std::move(stage));
				for (std::size_t occurrence_index = 0;
					occurrence_index < face.boundary_occurrences.size(); ++occurrence_index)
				{
					const auto& occurrence = face.boundary_occurrences[occurrence_index];
					if (occurrence.occurrence_id != occurrence_index
						|| occurrence.edge.IsNull()
						|| !edge_ids.contains(occurrence.source_edge_id)
						|| !valid_orientation(occurrence.orientation))
					{
						error = "source face has a malformed boundary occurrence table";
						return false;
					}
					if (!BRep_Tool::Degenerated(occurrence.edge))
						expected_occurrences.emplace(face.source_face_id,
							occurrence.occurrence_id);
				}
			}
			std::map<std::pair<FaceId, std::uint32_t>, std::uint32_t>
				topology_by_face_node;
			for (const auto& claim : topology.mesh_node_claims)
			{
				auto stage = stages.find(claim.source_face_id);
				if (stage == stages.end()
					|| claim.source_node_index
						== OcctFaceConstraintSample::no_source_node_index
					|| claim.topology_id >= topology.atomization.topology_node_count
					|| claim.contact_atom_id >= topology.atomization.atoms.size()
					|| topology.atomization.atoms[
						static_cast<std::size_t>(claim.contact_atom_id)].id
						!= claim.contact_atom_id
					|| claim.atom_sample_index >= topology.atomization.atoms[
						static_cast<std::size_t>(claim.contact_atom_id)].samples.size()
					|| topology.atomization.atoms[
						static_cast<std::size_t>(claim.contact_atom_id)].samples[
							claim.atom_sample_index].topology_id != claim.topology_id
					|| (claim.boundary_occurrence_id
						!= OcctContactBoundaryOccurrence::no_occurrence_id
						&& (claim.boundary_occurrence_id
							>= stage->second.source->boundary_occurrences.size()
							|| stage->second.source->boundary_occurrences[
								claim.boundary_occurrence_id].occurrence_id
								!= claim.boundary_occurrence_id)))
				{
					error = "contact topology has an invalid original face-mesh node claim";
					return false;
				}
				const auto node_key = std::make_pair(claim.source_face_id,
					claim.source_node_index);
				auto [node, node_inserted] = topology_by_face_node.emplace(node_key,
					claim.topology_id);
				if (!node_inserted && node->second != claim.topology_id)
				{
					error = "one original face-mesh node is claimed by conflicting topology IDs";
					return false;
				}
				const auto claim_key = std::make_tuple(claim.boundary_occurrence_id,
					claim.contact_atom_id, claim.atom_sample_index);
				auto [stored, inserted] =
					stage->second.source_node_by_occurrence_topology.emplace(claim_key,
						claim.source_node_index);
				if (!inserted && stored->second != claim.source_node_index)
				{
					error = "one contact sample/face occurrence names multiple original mesh nodes";
					return false;
				}
			}

			if (topology.atomization.atoms.size() > occt_constraint_payload_mask)
			{
				error = "contact atom count exceeds the constraint-tag payload";
				return false;
			}
			std::vector<std::optional<std::array<double, 3>>> topology_positions(
				topology.atomization.topology_node_count);
			std::set<FaceOccurrence> observed_occurrences;
			for (std::size_t atom_index = 0;
				atom_index < topology.atomization.atoms.size(); ++atom_index)
			{
				const auto& atom = topology.atomization.atoms[atom_index];
				if (atom.id != atom_index || !occt_contact_atom_constraint_id(atom.id))
				{
					error = "contact atoms are not dense or exceed the tag payload";
					return false;
				}
				if (atom.samples.size() < 2 || atom.source_spans.empty()
					|| atom.face_uses.empty())
				{
					error = "contact atom has no finite chain, source span, or face use";
					return false;
				}
				std::set<std::uint32_t> atom_source_edges;
				for (const auto& span : atom.source_spans)
				{
					if (!edge_ids.contains(span.source_edge_id)
						|| !finite(span.parameter_begin) || !finite(span.parameter_end)
						|| span.parameter_begin == span.parameter_end
						|| !atom_source_edges.insert(span.source_edge_id).second)
					{
						error = "contact atom has a malformed or duplicate source span";
						return false;
					}
				}
				for (const auto& sample : atom.samples)
				{
					if (!finite(sample.canonical_world_position)
						|| !finite(sample.canonical_position_tolerance)
						|| sample.canonical_position_tolerance < 0.0
						|| sample.topology_id >= topology_positions.size())
					{
						error = "contact atom has an invalid canonical topology sample";
						return false;
					}
					auto& position = topology_positions[sample.topology_id];
					if (!position) position = sample.canonical_world_position;
					else if (!same_position_bits(*position, sample.canonical_world_position))
					{
						error = "one contact topology ID has conflicting canonical position bits";
						return false;
					}
				}
				std::uint32_t fan = 0;
				std::set<FaceId> atom_faces;
				for (const auto& face_use : atom.face_uses)
				{
					auto stage = stages.find(face_use.source_face_id);
					if (stage == stages.end()
						|| !atom_faces.insert(face_use.source_face_id).second
						|| !finite(face_use.exact_operation_tolerance)
						|| face_use.exact_operation_tolerance < 0.0
						|| (face_use.sector_count != 1 && face_use.sector_count != 2)
						|| (face_use.location != OcctFaceIntervalLocation::boundary
							&& face_use.location != OcctFaceIntervalLocation::interior))
					{
						error = "contact atom has a malformed or duplicate face use";
						return false;
					}
					fan += face_use.sector_count;
					std::set<std::uint32_t> use_occurrences;
					for (const auto& use : face_use.boundary_occurrences)
					{
						if (!atom_source_edges.contains(use.source_edge_id)
							|| !valid_orientation(use.orientation)
							|| !use_occurrences.insert(use.occurrence_id).second
							|| use.occurrence_id
								>= stage->second.source->boundary_occurrences.size())
						{
							error = "contact atom has a malformed boundary occurrence use";
							return false;
						}
						const auto& occurrence =
							stage->second.source->boundary_occurrences[use.occurrence_id];
						if (occurrence.occurrence_id != use.occurrence_id
							|| occurrence.orientation != use.orientation)
						{
							error = "atom boundary use does not match its exact face occurrence";
							return false;
						}
						if (!append_boundary_constraints(stage->second, occurrence, use,
							atom, face_use.exact_operation_tolerance, error))
						{
							error = "face " + std::to_string(face_use.source_face_id)
								+ ", atom " + std::to_string(atom.id) + ": " + error;
							return false;
						}
						observed_occurrences.emplace(face_use.source_face_id,
							use.occurrence_id);
					}
					if (face_use.boundary_occurrences.empty())
					{
						if (face_use.location != OcctFaceIntervalLocation::interior)
						{
							error = "a boundary face use has no exact boundary occurrence";
							return false;
						}
						if (!append_interior_constraint(stage->second, atom,
							face_use.exact_operation_tolerance, error))
						{
							error = "face " + std::to_string(face_use.source_face_id)
								+ ", atom " + std::to_string(atom.id) + ": " + error;
							return false;
						}
					}
				}
				if (fan != atom.fan_degree || fan == 0)
				{
					error = "contact atom fan degree disagrees with its distinct face sectors";
					return false;
				}
			}
			if (std::any_of(topology_positions.begin(), topology_positions.end(),
				[](const auto& value) { return !value.has_value(); }))
			{
				error = "contact topology IDs are not dense";
				return false;
			}
			if (observed_occurrences != expected_occurrences)
			{
				error = "contact atoms do not cover every exact source-face boundary occurrence";
				return false;
			}

			OcctConformingMesh candidate;
			candidate.shared_topology_node_count =
				topology.atomization.topology_node_count;
			candidate.faces.reserve(stages.size());
			for (auto& [face_id, stage] : stages)
			{
				TopLoc_Location location;
				const Handle(Poly_Triangulation) triangulation =
					BRep_Tool::Triangulation(stage.source->face, location);
				if (triangulation.IsNull())
				{
					error = "source face " + std::to_string(face_id)
						+ " has no OCCT triangulation";
					return false;
				}
				std::vector<OcctFaceConstraintPolyline> constraints;
				constraints.reserve(stage.cad_constraints.size()
					+ stage.atom_constraints.size());
				constraints.insert(constraints.end(), stage.cad_constraints.begin(),
					stage.cad_constraints.end());
				constraints.insert(constraints.end(), stage.atom_constraints.begin(),
					stage.atom_constraints.end());
				if (!normalize_chart_branches(constraints, error))
				{
					error = "source face " + std::to_string(face_id)
						+ " chart branches failed: " + error;
					return false;
				}
				OcctConformingFaceMesh face_output;
				face_output.source_face_id = face_id;
				std::string face_error;
				if (!conform_occt_face_triangulation(stage.source->face, triangulation,
					location, constraints, face_output.mesh, face_error, options))
				{
					error = "source face " + std::to_string(face_id)
						+ " conforming failed: " + face_error;
					return false;
				}
				if (!decode_edges(topology, stage, face_output, error))
				{
					error = "source face " + std::to_string(face_id)
						+ " provenance failed: " + error;
					return false;
				}
				candidate.faces.push_back(std::move(face_output));
			}
			output = std::move(candidate);
			return true;
		}
		catch (const Standard_Failure& failure)
		{
			const char* message = failure.what();
			error = std::string("OpenCascade conforming-mesh failure: ")
				+ (message ? message : "unknown Standard_Failure");
			return false;
		}
		catch (const std::exception& exception)
		{
			error = std::string("conforming-mesh failure: ") + exception.what();
			return false;
		}
	}
}
