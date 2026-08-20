#include "core/geometry/occt_conforming_import_assembly.h"

#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Pnt2d.hxx>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <limits>
#include <iomanip>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

namespace paracfd::core
{
	namespace
	{
		constexpr double mm_to_m = 0.001;
		constexpr std::uint32_t no_u32 = std::numeric_limits<std::uint32_t>::max();

		using LocalEdge = std::array<std::uint32_t, 2>;
		using ContactSegment = std::tuple<std::uint64_t, std::uint32_t, std::uint32_t>;
		using FaceCad = std::pair<std::uint32_t, std::uint32_t>;

		struct IndexedUv
		{
			UvPoint uv{};
			bool all_finite = true;
			bool ambiguous = false;
		};

		struct IndexedConformingFace
		{
			const OcctConformingFaceMesh* face = nullptr;
			std::map<std::uint32_t, IndexedUv> shared_uv;
		};

		struct AssemblyLookup
		{
			std::vector<const OcctContactTopologyFace*> topology_faces;
			std::vector<IndexedConformingFace> conformed_faces;
			std::vector<const OcctContactTopologyEdge*> topology_edges;
		};

		LocalEdge edge_key(std::uint32_t a, std::uint32_t b)
		{
			if (b < a) std::swap(a, b);
			return {{a, b}};
		}

		ContactSegment contact_segment_key(std::uint64_t atom,
			std::uint32_t a, std::uint32_t b)
		{
			if (b < a) std::swap(a, b);
			return {atom, a, b};
		}

		bool finite(double value)
		{
			return std::isfinite(value);
		}

		bool finite(const std::array<double, 3>& value)
		{
			return finite(value[0]) && finite(value[1]) && finite(value[2]);
		}

		bool finite(const UvPoint& value)
		{
			return finite(value.u) && finite(value.v);
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

		bool valid_orientation(std::int8_t orientation)
		{
			return orientation == static_cast<std::int8_t>(TopAbs_FORWARD)
				|| orientation == static_cast<std::int8_t>(TopAbs_REVERSED);
		}

		double edge_tolerance_m(const OcctContactTopologyEdge& edge)
		{
			const double scale = std::abs(
				edge.edge.Location().Transformation().ScaleFactor());
			return std::max(0.0, BRep_Tool::Tolerance(edge.edge)) * scale * mm_to_m;
		}

		bool boundary_branch_uv(const OcctContactTopologyFace& face,
			const OcctContactTopologyBoundaryOccurrence& occurrence,
			const OcctContactAtomBoundaryUse& use, const OcctContactAtom& atom,
			std::vector<std::array<double, 2>>& output, std::string& error)
		{
			if (use.sample_target_parameters.size() != atom.samples.size())
			{
				error = "periodic/contact boundary occurrence has the wrong parameter count";
				return false;
			}
			double first = 0.0, last = 0.0;
			const Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(
				occurrence.edge, face.face, first, last);
			if (pcurve.IsNull() || !finite(first) || !finite(last) || !(last > first))
			{
				error = "exact boundary occurrence has no finite face pcurve";
				return false;
			}
			output.clear();
			output.reserve(atom.samples.size());
			for (std::size_t sample = 0; sample < atom.samples.size(); ++sample)
			{
				const double parameter = use.sample_target_parameters[sample];
				const double scale = std::max({std::abs(first), std::abs(last),
					std::abs(parameter), std::abs(last - first),
					std::numeric_limits<double>::min()});
				const double roundoff = 64.0 * std::numeric_limits<double>::epsilon()
					* scale;
				if (!finite(parameter) || parameter < first - roundoff
					|| parameter > last + roundoff)
				{
					error = "exact boundary sample parameter escaped its pcurve range";
					return false;
				}
				const gp_Pnt2d uv = pcurve->Value(parameter);
				if (!finite(uv.X()) || !finite(uv.Y()))
				{
					error = "exact boundary pcurve produced a non-finite UV sample";
					return false;
				}
				output.push_back({{uv.X(), uv.Y()}});
			}
			return true;
		}

		bool face_mesh_uv_chain(const IndexedConformingFace& indexed_face,
			const OcctContactAtom& atom,
			std::vector<std::array<double, 2>>& output, std::string& error)
		{
			output.clear();
			output.reserve(atom.samples.size());
			for (const auto& sample : atom.samples)
			{
				const auto found = indexed_face.shared_uv.find(sample.topology_id);
				if (found == indexed_face.shared_uv.end())
				{
					error = "contact topology sample is absent from its conformed face";
					return false;
				}
				if (!found->second.all_finite)
				{
					error = "contact topology sample has no finite conformed-face UV";
					return false;
				}
				if (found->second.ambiguous)
				{
					error = "interior contact has unresolved multiple UV-chart branches";
					return false;
				}
				output.push_back({{found->second.uv.u, found->second.uv.v}});
			}
			return true;
		}

		bool make_uv_branches(const AssemblyLookup& lookup,
			const OcctContactAtom& atom,
			const OcctContactAtomFaceUse& use,
			std::vector<OcctConformingContactUvBranch>& branches,
			std::vector<std::array<double, 2>>& logical_uv, std::string& error)
		{
			if (use.source_face_id >= lookup.topology_faces.size()
				|| use.source_face_id >= lookup.conformed_faces.size())
			{
				error = "contact face use is absent from topology or conformed mesh";
				return false;
			}
			const auto* face = lookup.topology_faces[use.source_face_id];
			const auto& indexed_mesh = lookup.conformed_faces[use.source_face_id];
			if (!face || !indexed_mesh.face)
			{
				error = "contact face use is absent from topology or conformed mesh";
				return false;
			}
			const StepContactUseKind kind = use.location
				== OcctFaceIntervalLocation::boundary
				? StepContactUseKind::trim_boundary : StepContactUseKind::face_interior;
			logical_uv.clear();
			if (use.boundary_occurrences.empty())
			{
				OcctConformingContactUvBranch branch;
				branch.contact_atom_id = atom.id;
				branch.source_face_id = use.source_face_id;
				branch.kind = kind;
				if (!face_mesh_uv_chain(indexed_mesh, atom,
					branch.sample_uv, error)) return false;
				logical_uv = branch.sample_uv;
				branches.push_back(std::move(branch));
				return true;
			}

			std::vector<const OcctContactAtomBoundaryUse*> ordered;
			ordered.reserve(use.boundary_occurrences.size());
			for (const auto& occurrence_use : use.boundary_occurrences)
				ordered.push_back(&occurrence_use);
			std::sort(ordered.begin(), ordered.end(), [](const auto* a, const auto* b)
			{
				return std::tie(a->occurrence_id, a->source_edge_id, a->orientation)
					< std::tie(b->occurrence_id, b->source_edge_id, b->orientation);
			});
			for (const auto* occurrence_use : ordered)
			{
				if (!occurrence_use
					|| occurrence_use->occurrence_id >= face->boundary_occurrences.size())
				{
					error = "contact face use references an invalid boundary occurrence";
					return false;
				}
				const auto& occurrence =
					face->boundary_occurrences[occurrence_use->occurrence_id];
				if (occurrence.occurrence_id != occurrence_use->occurrence_id
					|| occurrence.source_edge_id != occurrence_use->source_edge_id
					|| occurrence.orientation != occurrence_use->orientation
					|| !valid_orientation(occurrence.orientation))
				{
					error = "contact branch disagrees with its exact oriented occurrence";
					return false;
				}
				OcctConformingContactUvBranch branch;
				branch.contact_atom_id = atom.id;
				branch.source_face_id = use.source_face_id;
				branch.boundary_occurrence_id = occurrence.occurrence_id;
				branch.boundary_orientation = occurrence.orientation;
				branch.kind = kind;
				if (!boundary_branch_uv(*face, occurrence, *occurrence_use, atom,
					branch.sample_uv, error)) return false;
				if (logical_uv.empty()) logical_uv = branch.sample_uv;
				branches.push_back(std::move(branch));
			}
			if (logical_uv.empty())
			{
				error = "contact face use produced no exact UV branch";
				return false;
			}
			return true;
		}

		bool validate_inputs(const OcctContactTopology& topology,
			const OcctConformingMesh& conformed,
			std::vector<std::array<double, 3>>& canonical_positions,
			AssemblyLookup& lookup,
			std::string& error)
		{
			if (topology.faces.empty() || conformed.faces.empty()
				|| conformed.faces.size() != topology.faces.size())
			{
				error = "conforming import requires one conformed mesh for every source face";
				return false;
			}
			if (conformed.shared_topology_node_count
				!= topology.atomization.topology_node_count)
			{
				error = "conformed mesh and contact atomization disagree on shared-node count";
				return false;
			}
			canonical_positions.assign(conformed.shared_topology_node_count,
				{{0.0, 0.0, 0.0}});
			std::vector<bool> known(canonical_positions.size(), false);
			for (std::size_t atom_index = 0;
				atom_index < topology.atomization.atoms.size(); ++atom_index)
			{
				const auto& atom = topology.atomization.atoms[atom_index];
				if (atom.id != atom_index || atom.samples.size() < 2
					|| atom.face_uses.empty() || atom.fan_degree == 0)
				{
					error = "contact atom table is not dense or contains an empty atom";
					return false;
				}
				std::uint32_t fan = 0;
				std::set<std::uint32_t> faces;
				for (const auto& use : atom.face_uses)
				{
					if (!faces.insert(use.source_face_id).second
						|| (use.sector_count != 1 && use.sector_count != 2))
					{
						error = "contact atom has a duplicate/malformed face use";
						return false;
					}
					fan += use.sector_count;
				}
				if (fan != atom.fan_degree)
				{
					error = "contact atom fan degree disagrees with its face sectors";
					return false;
				}
				for (const auto& sample : atom.samples)
				{
					if (sample.topology_id >= canonical_positions.size()
						|| !finite(sample.canonical_world_position))
					{
						error = "contact atom contains an invalid canonical sample";
						return false;
					}
					auto& position = canonical_positions[sample.topology_id];
					if (!known[sample.topology_id])
					{
						position = sample.canonical_world_position;
						known[sample.topology_id] = true;
					}
					else if (!same_position_bits(position,
						sample.canonical_world_position))
					{
						error = "one shared topology ID has conflicting canonical position bits";
						return false;
					}
				}
			}
			if (std::any_of(known.begin(), known.end(), [](bool value) { return !value; }))
			{
				error = "shared topology IDs are not dense";
				return false;
			}

			lookup.topology_faces.assign(topology.faces.size(), nullptr);
			for (const auto& face : topology.faces)
			{
				if (face.face.IsNull() || face.source_face_id >= topology.faces.size()
					|| lookup.topology_faces[face.source_face_id])
				{
					error = "contact topology contains a null or duplicate source face";
					return false;
				}
				lookup.topology_faces[face.source_face_id] = &face;
			}
			if (std::any_of(lookup.topology_faces.begin(), lookup.topology_faces.end(),
				[](const auto* face) { return !face; }))
			{
				error = "contact topology source face IDs are not dense";
				return false;
			}

			lookup.conformed_faces.assign(topology.faces.size(), {});
			for (const auto& face : conformed.faces)
			{
				if (face.source_face_id >= lookup.topology_faces.size()
					|| !lookup.topology_faces[face.source_face_id]
					|| lookup.conformed_faces[face.source_face_id].face)
				{
					error = "conformed mesh contains an unknown or duplicate source face";
					return false;
				}
				const std::size_t vertices = face.mesh.world_positions.size();
				if (vertices == 0 || face.mesh.uv.size() != vertices
					|| face.mesh.topology_ids.size() != vertices
					|| face.edge_provenance.size() != face.mesh.edges.size())
				{
					error = "conformed face has inconsistent vertex/edge arrays";
					return false;
				}
				auto& indexed = lookup.conformed_faces[face.source_face_id];
				indexed.face = &face;
				for (std::size_t vertex = 0; vertex < vertices; ++vertex)
				{
					const std::uint32_t topology_id = face.mesh.topology_ids[vertex];
					if (topology_id == UvVertex::no_topology_id) continue;
					const UvPoint& uv = face.mesh.uv[vertex];
					auto [found, inserted] = indexed.shared_uv.try_emplace(
						topology_id, IndexedUv{uv, finite(uv), false});
					if (inserted) continue;
					if (!finite(uv)) found->second.all_finite = false;
					else if (found->second.all_finite
						&& !same_uv_bits(found->second.uv, uv))
						found->second.ambiguous = true;
				}
			}
			if (std::any_of(lookup.conformed_faces.begin(),
				lookup.conformed_faces.end(),
				[](const auto& face) { return !face.face; }))
			{
				error = "conforming import requires one conformed mesh for every source face";
				return false;
			}

			lookup.topology_edges.assign(topology.edges.size(), nullptr);
			for (const auto& edge : topology.edges)
			{
				if (edge.edge.IsNull() || edge.source_edge_id >= topology.edges.size()
					|| lookup.topology_edges[edge.source_edge_id])
				{
					error = "contact topology contains a null or duplicate source CAD edge";
					return false;
				}
				lookup.topology_edges[edge.source_edge_id] = &edge;
			}
			if (std::any_of(lookup.topology_edges.begin(), lookup.topology_edges.end(),
				[](const auto* edge) { return !edge; }))
			{
				error = "contact topology source CAD edge IDs are not dense";
				return false;
			}
			return true;
		}

		bool assemble_graph(const OcctContactTopology& topology,
			const AssemblyLookup& lookup, StepCadContactGraph& graph,
			std::vector<OcctConformingContactUvBranch>& branches,
			std::string& error)
		{
			for (const auto& atom : topology.atomization.atoms)
			{
				StepContactCurve curve;
				curve.id = atom.id;
				curve.fan_degree = atom.fan_degree;
				for (const auto& span : atom.source_spans)
					curve.source_edge_ids.push_back(span.source_edge_id);
				std::sort(curve.source_edge_ids.begin(), curve.source_edge_ids.end());
				curve.source_edge_ids.erase(std::unique(curve.source_edge_ids.begin(),
					curve.source_edge_ids.end()), curve.source_edge_ids.end());
				if (curve.source_edge_ids.empty())
				{
					error = "contact atom has no source CAD edge span";
					return false;
				}
				curve.source_edge_id = curve.source_edge_ids.front();
				const OcctContactAtomSourceSpan* canonical_span = nullptr;
				for (const auto& span : atom.source_spans)
					if (span.source_edge_id == curve.source_edge_id)
					{
						if (canonical_span)
						{
							error = "contact atom repeats its canonical source-edge span";
							return false;
						}
						canonical_span = &span;
					}
				if (!canonical_span || !finite(canonical_span->parameter_begin)
					|| !finite(canonical_span->parameter_end)
					|| canonical_span->parameter_begin == canonical_span->parameter_end)
				{
					error = "contact atom has no finite canonical source-edge span";
					return false;
				}
				for (const auto source_edge_id : curve.source_edge_ids)
				{
					const auto* edge = source_edge_id < lookup.topology_edges.size()
						? lookup.topology_edges[source_edge_id] : nullptr;
					if (!edge || edge->edge.IsNull())
					{
						error = "contact atom references an unknown source CAD edge";
						return false;
					}
					curve.tolerance_m = std::max(curve.tolerance_m,
						edge_tolerance_m(*edge));
				}
				for (std::size_t sample_index = 0;
					sample_index < atom.samples.size(); ++sample_index)
				{
					const auto& sample = atom.samples[sample_index];
					curve.sample_positions_m.push_back({{
						sample.canonical_world_position[0] * mm_to_m,
						sample.canonical_world_position[1] * mm_to_m,
						sample.canonical_world_position[2] * mm_to_m}});
					curve.sample_topology_ids.push_back(sample.topology_id);
					curve.tolerance_m = std::max(curve.tolerance_m,
						sample.canonical_position_tolerance * mm_to_m);
					std::vector<double> parameters;
					for (const auto& source : sample.source_locations)
						if (source.source_edge_id == curve.source_edge_id)
							parameters.push_back(source.parameter);
					std::sort(parameters.begin(), parameters.end());
					std::vector<double> unique_parameters;
					for (const double value : parameters)
					{
						const double scale = std::max({std::abs(value),
							std::abs(canonical_span->parameter_begin),
							std::abs(canonical_span->parameter_end),
							std::abs(canonical_span->parameter_end
								- canonical_span->parameter_begin),
							std::numeric_limits<double>::min()});
						const double representable_roundoff = 64.0
							* std::numeric_limits<double>::epsilon() * scale;
						if (unique_parameters.empty()
							|| std::abs(value - unique_parameters.back())
								> representable_roundoff)
							unique_parameters.push_back(value);
					}
					parameters = std::move(unique_parameters);
					if (parameters.empty()
						|| std::any_of(parameters.begin(), parameters.end(),
							[](double value) { return !finite(value); }))
					{
						error = "contact sample has no canonical source-edge parameter";
						return false;
					}
					double parameter = parameters.front();
					if (parameters.size() > 1)
					{
						// A closed CAD edge legitimately maps its one topological seam
						// vertex to both finite-range endpoints.  The atom chain keeps
						// two endpoint samples with one shared topology ID; select the
						// endpoint consistent with the atom's physical direction.  A
						// multiply-mapped interior sample remains ambiguous and fails.
						if (sample_index != 0 && sample_index + 1 != atom.samples.size())
						{
							std::ostringstream message;
							message << "atom " << atom.id << " interior sample "
								<< sample_index
								<< " has conflicting canonical-edge parameters";
							message << std::setprecision(17);
							for (const double value : parameters) message << ' ' << value;
							error = message.str();
							return false;
						}
						const double endpoint = sample_index == 0
							? canonical_span->parameter_begin
							: canonical_span->parameter_end;
						const auto closest = std::min_element(parameters.begin(), parameters.end(),
							[&](double a, double b)
							{
								return std::abs(a - endpoint) < std::abs(b - endpoint);
							});
						parameter = *closest;
						if (parameter != endpoint)
						{
							error = "closed-edge contact endpoint does not retain its exact span parameter";
							return false;
						}
					}
					curve.source_parameters.push_back(parameter);
				}
				std::uint32_t logical_fan = 0;
				for (const auto& face_use : atom.face_uses)
				{
					StepContactFaceUse use;
					use.source_face_id = face_use.source_face_id;
					use.kind = face_use.location == OcctFaceIntervalLocation::boundary
						? StepContactUseKind::trim_boundary : StepContactUseKind::face_interior;
					use.sector_count = face_use.sector_count;
					if (!make_uv_branches(lookup, atom, face_use,
						branches, use.sample_uv, error)) return false;
					curve.uses.push_back(std::move(use));
					logical_fan += face_use.sector_count;
				}
				if (logical_fan != atom.fan_degree)
				{
					error = "assembled contact face uses disagree with atom fan degree";
					return false;
				}
				// A fan-one atom is an intentional fabric boundary/opening.  It remains
				// represented by CAD provenance and the UV audit sidecar, but is not a
				// public multi-sheet contact curve and must not reach EB as one.
				if (atom.fan_degree >= 2) graph.curves.push_back(std::move(curve));
			}
			return true;
		}

		bool assemble_mesh(const OcctContactTopology& topology,
			const OcctConformingMesh& conformed,
			const AssemblyLookup& lookup,
			const std::vector<std::array<double, 3>>& canonical_positions,
			TriMesh& mesh, std::string& error)
		{
			std::map<std::uint32_t, std::uint32_t> occurrence_count_by_edge;
			std::map<FaceCad, std::uint32_t> occurrence_count_by_face_edge;
			for (const auto& face : topology.faces)
				for (const auto& occurrence : face.boundary_occurrences)
				{
					if (occurrence.source_edge_id >= lookup.topology_edges.size()
						|| !lookup.topology_edges[occurrence.source_edge_id])
					{
						error = "source face occurrence references an unknown CAD edge";
						return false;
					}
					++occurrence_count_by_edge[occurrence.source_edge_id];
					++occurrence_count_by_face_edge[{face.source_face_id,
						occurrence.source_edge_id}];
				}

			std::map<ContactSegment, std::uint64_t> expected_contact_incidence;
			std::map<ContactSegment, std::uint64_t> observed_contact_incidence;
			for (const auto& atom : topology.atomization.atoms)
				for (std::size_t sample = 1; sample < atom.samples.size(); ++sample)
				{
					const auto key = contact_segment_key(atom.id,
						atom.samples[sample - 1].topology_id,
						atom.samples[sample].topology_id);
					expected_contact_incidence[key] += atom.fan_degree;
				}

			std::uint64_t next_topology_id = conformed.shared_topology_node_count;
			std::map<std::uint32_t, std::array<double, 3>> published_shared_positions;
			for (const auto& face : conformed.faces)
			{
				const std::size_t vertex_count = face.mesh.world_positions.size();
				if (mesh.vertex_count() > std::numeric_limits<std::uint32_t>::max()
					- vertex_count)
				{
					error = "conforming import exceeds the 32-bit mesh-index space";
					return false;
				}
				const std::uint32_t base = static_cast<std::uint32_t>(mesh.vertex_count());
				for (std::size_t vertex = 0; vertex < vertex_count; ++vertex)
				{
					const auto& world = face.mesh.world_positions[vertex];
					const auto& uv = face.mesh.uv[vertex];
					if (!finite(world) || !finite(uv))
					{
						error = "conformed face contains a non-finite vertex or UV";
						return false;
					}
					std::uint32_t topology_id = face.mesh.topology_ids[vertex];
					if (topology_id == UvVertex::no_topology_id)
					{
						if (next_topology_id >= UvVertex::no_topology_id)
						{
							error = "conforming import exhausted ordinary topology IDs";
							return false;
						}
						topology_id = static_cast<std::uint32_t>(next_topology_id++);
					}
					else
					{
						if (topology_id >= canonical_positions.size()
							|| !same_position_bits(world, canonical_positions[topology_id]))
						{
							error = "conformed shared vertex differs from its canonical position bits";
							return false;
						}
						auto [found, inserted] = published_shared_positions.emplace(
							topology_id, world);
						if (!inserted && !same_position_bits(found->second, world))
						{
							error = "one shared topology ID has non-identical conformed positions";
							return false;
						}
					}
					mesh.positions_fp64.insert(mesh.positions_fp64.end(), {
						world[0] * mm_to_m, world[1] * mm_to_m, world[2] * mm_to_m});
					mesh.positions.insert(mesh.positions.end(), {
						static_cast<float>(world[0] * mm_to_m),
						static_cast<float>(world[1] * mm_to_m),
						static_cast<float>(world[2] * mm_to_m)});
					mesh.vertex_uv.insert(mesh.vertex_uv.end(), {
						static_cast<float>(uv.u), static_cast<float>(uv.v)});
					mesh.topology_vertex_ids.push_back(topology_id);
				}

				std::map<LocalEdge, const OcctConformedMeshEdgeProvenance*> provenance;
				for (std::size_t edge_index = 0;
					edge_index < face.edge_provenance.size(); ++edge_index)
				{
					const auto& edge = face.edge_provenance[edge_index];
					if (edge.vertices != face.mesh.edges[edge_index].vertices
						|| edge.vertices[0] >= vertex_count
						|| edge.vertices[1] >= vertex_count
						|| edge.vertices[0] == edge.vertices[1]
						|| !provenance.emplace(edge_key(edge.vertices[0], edge.vertices[1]),
							&edge).second)
					{
						error = "conformed face has malformed or duplicate edge provenance";
						return false;
					}
					if (edge.boundary != face.mesh.edges[edge_index].boundary)
					{
						error = "decoded edge boundary state differs from conformed edge table";
						return false;
					}
				}
				std::map<LocalEdge, std::uint32_t> local_incidence;
				for (const auto& triangle : face.mesh.triangles)
				{
					for (const auto vertex : triangle)
						if (vertex >= vertex_count)
						{
							error = "conformed triangle references an invalid local vertex";
							return false;
						}
					for (std::size_t side = 0; side < 3; ++side)
						++local_incidence[edge_key(triangle[side],
							triangle[(side + 1) % 3])];
				}
				if (local_incidence.size() != provenance.size())
				{
					error = "conformed edge table is not exactly the triangle edge set";
					return false;
				}
				for (const auto& [key, edge] : provenance)
				{
					const auto incidence = local_incidence.find(key);
					if (incidence == local_incidence.end()
						|| incidence->second != (edge->boundary ? 1u : 2u)
						|| (edge->boundary && !edge->source_cad_edge_id)
						|| (!edge->boundary && edge->source_cad_edge_id))
					{
						error = "conformed edge incidence or CAD-boundary lineage is inconsistent";
						return false;
					}
				}

				for (const auto& source_triangle : face.mesh.triangles)
				{
					std::array<std::uint32_t, 3> triangle = source_triangle;
					if (face.mesh.face_reversed) std::swap(triangle[1], triangle[2]);
					for (const auto vertex : triangle) mesh.indices.push_back(base + vertex);
					mesh.source_face_ids.push_back(face.source_face_id);
					for (std::size_t side = 0; side < 3; ++side)
					{
						const std::uint32_t a = triangle[side];
						const std::uint32_t b = triangle[(side + 1) % 3];
						const auto found = provenance.find(edge_key(a, b));
						if (found == provenance.end())
						{
							error = "oriented triangle half-edge is absent from provenance table";
							return false;
						}
						const auto& edge = *found->second;
						const bool known = edge.boundary
							&& edge.source_cad_edge_id.has_value();
						mesh.triangle_cad_edge_provenance_states.push_back(
							static_cast<std::uint8_t>(known
								? CadEdgeProvenanceState::known
								: CadEdgeProvenanceState::none));
						if (known)
						{
							const auto* cad = *edge.source_cad_edge_id
								< lookup.topology_edges.size()
								? lookup.topology_edges[*edge.source_cad_edge_id]
								: nullptr;
							if (!cad
								|| occurrence_count_by_edge[*edge.source_cad_edge_id] == 0)
							{
								error = "boundary half-edge references an unowned CAD edge";
								return false;
							}
							mesh.triangle_cad_edge_ids.push_back(*edge.source_cad_edge_id);
							mesh.triangle_cad_edge_incident_face_counts.push_back(
								occurrence_count_by_edge[*edge.source_cad_edge_id]);
							mesh.triangle_cad_edge_tolerances.push_back(
								edge_tolerance_m(*cad));
							mesh.triangle_cad_edge_is_periodic_seam.push_back(
								occurrence_count_by_face_edge[{face.source_face_id,
									*edge.source_cad_edge_id}] == 2 ? 1u : 0u);
						}
						else
						{
							mesh.triangle_cad_edge_ids.push_back(TriMesh::kNoCadEdgeId);
							mesh.triangle_cad_edge_incident_face_counts.push_back(0);
							mesh.triangle_cad_edge_tolerances.push_back(0.0);
							mesh.triangle_cad_edge_is_periodic_seam.push_back(0);
						}

						std::uint64_t public_contact = TriMesh::kNoCadContactId;
						std::uint32_t public_fan = 0;
						std::uint64_t exact_atom = TriMesh::kNoCadEdgeAtomId;
						if (edge.contact_atom_id)
						{
							if (*edge.contact_atom_id >= topology.atomization.atoms.size()
								|| topology.atomization.atoms[
									static_cast<std::size_t>(*edge.contact_atom_id)].id
									!= *edge.contact_atom_id)
							{
								error = "half-edge references an unknown contact atom";
								return false;
							}
							const auto& atom = topology.atomization.atoms[
								static_cast<std::size_t>(*edge.contact_atom_id)];
							exact_atom = atom.id;
							const std::uint32_t topology_a =
								mesh.topology_vertex_ids[base + a];
							const std::uint32_t topology_b =
								mesh.topology_vertex_ids[base + b];
							const auto key = contact_segment_key(atom.id,
								topology_a, topology_b);
							if (!expected_contact_incidence.contains(key))
							{
								error = "contact-tagged mesh edge is not an exact atom segment";
								return false;
							}
							++observed_contact_incidence[key];
							if (atom.fan_degree >= 2)
							{
								public_contact = atom.id;
								public_fan = atom.fan_degree;
							}
						}
						mesh.triangle_cad_edge_contact_ids.push_back(public_contact);
						mesh.triangle_cad_edge_certified_fan_degrees.push_back(public_fan);
						mesh.triangle_cad_edge_atom_ids.push_back(exact_atom);
					}
				}
			}
			if (observed_contact_incidence != expected_contact_incidence)
			{
				error = "conformed contact half-edge incidence disagrees with exact atom fan";
				return false;
			}
			if (published_shared_positions.size() != canonical_positions.size())
			{
				error = "one or more canonical shared topology nodes were not published";
				return false;
			}
			return true;
		}

		bool finish_mesh(TriMesh& mesh, std::string& error)
		{
			if (mesh.indices.empty() || mesh.positions_fp64.size() != 3 * mesh.vertex_count()
				|| mesh.vertex_uv.size() != 2 * mesh.vertex_count()
				|| mesh.topology_vertex_ids.size() != mesh.vertex_count()
				|| mesh.source_face_ids.size() != mesh.triangle_count()
				|| !mesh.has_cad_edge_provenance()
				|| !mesh.has_cad_edge_contact_provenance()
				|| !mesh.has_cad_edge_atom_provenance())
			{
				error = "assembled conforming mesh has inconsistent hand-off arrays";
				return false;
			}
			std::vector<double> accumulated(mesh.positions_fp64.size(), 0.0);
			for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
			{
				const std::uint32_t i0 = mesh.indices[3 * triangle];
				const std::uint32_t i1 = mesh.indices[3 * triangle + 1];
				const std::uint32_t i2 = mesh.indices[3 * triangle + 2];
				const double* p0 = &mesh.positions_fp64[3 * i0];
				const double* p1 = &mesh.positions_fp64[3 * i1];
				const double* p2 = &mesh.positions_fp64[3 * i2];
				const double e10 = p1[0] - p0[0], e11 = p1[1] - p0[1],
					e12 = p1[2] - p0[2];
				const double e20 = p2[0] - p0[0], e21 = p2[1] - p0[1],
					e22 = p2[2] - p0[2];
				const std::array<double, 3> normal{{
					e11 * e22 - e12 * e21,
					e12 * e20 - e10 * e22,
					e10 * e21 - e11 * e20}};
				const double magnitude = std::sqrt(normal[0] * normal[0]
					+ normal[1] * normal[1] + normal[2] * normal[2]);
				if (!finite(magnitude) || magnitude <= 0.0)
				{
					error = "conformed mesh contains a zero-area or non-finite triangle";
					return false;
				}
				for (const auto vertex : {i0, i1, i2})
					for (std::size_t axis = 0; axis < 3; ++axis)
						accumulated[3 * vertex + axis] += normal[axis];
			}
			mesh.normals.resize(mesh.positions.size());
			for (std::size_t vertex = 0; vertex < mesh.vertex_count(); ++vertex)
			{
				const double x = accumulated[3 * vertex];
				const double y = accumulated[3 * vertex + 1];
				const double z = accumulated[3 * vertex + 2];
				const double length = std::sqrt(x * x + y * y + z * z);
				if (!finite(length) || length <= 0.0)
				{
					error = "conformed mesh contains an unreferenced/degenerate vertex";
					return false;
				}
				mesh.normals[3 * vertex] = static_cast<float>(x / length);
				mesh.normals[3 * vertex + 1] = static_cast<float>(y / length);
				mesh.normals[3 * vertex + 2] = static_cast<float>(z / length);
			}
			std::array<double, 3> low{{
				std::numeric_limits<double>::max(),
				std::numeric_limits<double>::max(),
				std::numeric_limits<double>::max()}};
			std::array<double, 3> high{{
				-std::numeric_limits<double>::max(),
				-std::numeric_limits<double>::max(),
				-std::numeric_limits<double>::max()}};
			for (std::size_t vertex = 0; vertex < mesh.vertex_count(); ++vertex)
				for (std::size_t axis = 0; axis < 3; ++axis)
				{
					const double value = mesh.positions_fp64[3 * vertex + axis];
					low[axis] = std::min(low[axis], value);
					high[axis] = std::max(high[axis], value);
				}
			mesh.bbox_min = {{static_cast<float>(low[0]), static_cast<float>(low[1]),
				static_cast<float>(low[2])}};
			mesh.bbox_max = {{static_cast<float>(high[0]), static_cast<float>(high[1]),
				static_cast<float>(high[2])}};
			return true;
		}
	}

	bool assemble_occt_conforming_import(const OcctContactTopology& topology,
		const OcctConformingMesh& conformed,
		OcctConformingImportAssembly& output, std::string& error)
	{
		error.clear();
		try
		{
			std::vector<std::array<double, 3>> canonical_positions;
			AssemblyLookup lookup;
			if (!validate_inputs(topology, conformed, canonical_positions,
				lookup, error))
				return false;
			OcctConformingImportAssembly candidate;
			if (!assemble_graph(topology, lookup, candidate.contacts,
				candidate.contact_uv_branches, error)) return false;
			if (!assemble_mesh(topology, conformed, lookup, canonical_positions,
				candidate.mesh, error)) return false;
			if (!finish_mesh(candidate.mesh, error)) return false;
			output = std::move(candidate);
			return true;
		}
		catch (const Standard_Failure& failure)
		{
			const char* message = failure.what();
			error = std::string("OpenCascade conforming import assembly failure: ")
				+ (message ? message : "unknown Standard_Failure");
			return false;
		}
		catch (const std::exception& exception)
		{
			error = std::string("conforming import assembly failure: ")
				+ exception.what();
			return false;
		}
	}
}
