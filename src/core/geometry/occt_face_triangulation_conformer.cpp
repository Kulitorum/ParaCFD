#include "core/geometry/occt_face_triangulation_conformer.h"

#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <Poly_Triangle.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <utility>

namespace paracfd::core
{
	namespace
	{
		struct SurfaceChart
		{
			occ::handle<Geom_Surface> surface;
			TopLoc_Location location;
			double position_tolerance = 0.0;
			double uv_tolerance = 0.0;
			bool u_periodic = false;
			bool v_periodic = false;
			double u_period = 0.0;
			double v_period = 0.0;
		};

		struct CanonicalRecord
		{
			std::array<double, 3> position{};
			double tolerance = 0.0;
		};

		struct PreparedConstraints
		{
			std::vector<std::vector<UvConstraintSample>> samples;
			// Internal UvConstraintMesh IDs are unique per chart vertex. Entries here
			// translate temporary IDs either to an explicitly aliased physical topology
			// ID or back to the unclaimed sentinel on export.
			std::map<std::uint32_t, std::uint32_t> output_topology_by_internal;
		};

		bool finite(double value)
		{
			return std::isfinite(value);
		}

		bool finite(const UvPoint& point)
		{
			return finite(point.u) && finite(point.v);
		}

		bool finite(const std::array<double, 3>& point)
		{
			return finite(point[0]) && finite(point[1]) && finite(point[2]);
		}

		std::array<double, 3> coordinates(const gp_Pnt& point)
		{
			return {point.X(), point.Y(), point.Z()};
		}

		gp_Pnt point(const std::array<double, 3>& coordinates)
		{
			return {coordinates[0], coordinates[1], coordinates[2]};
		}

		bool same_position(const std::array<double, 3>& a,
			const std::array<double, 3>& b)
		{
			return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
		}

		bool finite_transform(const gp_Trsf& transform)
		{
			for (int row = 1; row <= 3; ++row)
				for (int column = 1; column <= 4; ++column)
					if (!finite(transform.Value(row, column))) return false;
			const double scale = std::abs(transform.ScaleFactor());
			return finite(scale) && scale > std::numeric_limits<double>::min();
		}

		double numerical_position_tolerance(const gp_Pnt& a, const gp_Pnt& b)
		{
			const double scale = std::max({std::abs(a.X()), std::abs(a.Y()), std::abs(a.Z()),
				std::abs(b.X()), std::abs(b.Y()), std::abs(b.Z()), a.Distance(b)});
			return 512.0 * std::numeric_limits<double>::epsilon()
				* std::max(scale, std::numeric_limits<double>::min());
		}

		double representable_spacing(double value)
		{
			const double up = std::nextafter(value,
				std::numeric_limits<double>::infinity());
			const double down = std::nextafter(value,
				-std::numeric_limits<double>::infinity());
			double spacing = 0.0;
			if (finite(up)) spacing = std::max(spacing, std::abs(up - value));
			if (finite(down)) spacing = std::max(spacing, std::abs(value - down));
			return spacing;
		}

		gp_Pnt evaluate_world(const SurfaceChart& chart, const UvPoint& uv)
		{
			gp_Pnt evaluated = chart.surface->Value(uv.u, uv.v);
			evaluated.Transform(chart.location.Transformation());
			return evaluated;
		}

		long double oriented_area2(const UvPoint& a, const UvPoint& b,
			const UvPoint& c)
		{
			return (static_cast<long double>(b.u) - a.u)
				* (static_cast<long double>(c.v) - a.v)
				- (static_cast<long double>(b.v) - a.v)
				* (static_cast<long double>(c.u) - a.u);
		}

		bool initialize_chart(const TopoDS_Face& face,
			const OcctFaceConformerOptions& options,
			SurfaceChart& chart, std::string& error)
		{
			chart.surface = BRep_Tool::Surface(face, chart.location);
			if (chart.surface.IsNull())
			{
				error = "face has no exact geometric surface";
				return false;
			}
			if (!finite_transform(chart.location.Transformation()))
			{
				error = "face has an invalid surface location transform";
				return false;
			}

			const double location_scale = std::abs(
				chart.location.Transformation().ScaleFactor());
			// The face and supplied caller-wide tolerance bound independent error
			// sources, so they add. A loose per-contact certificate is deliberately not
			// folded into this base: it must never validate another contact or source node.
			chart.position_tolerance = std::max(Precision::Confusion(),
				std::max(0.0, BRep_Tool::Tolerance(face)) * location_scale)
				+ options.geometric_tolerance;

			double u_min = 0.0, u_max = 0.0, v_min = 0.0, v_max = 0.0;
			BRepTools::UVBounds(face, u_min, u_max, v_min, v_max);
			if (!finite(u_min) || !finite(u_max) || !finite(v_min) || !finite(v_max)
				|| !(u_max > u_min) || !(v_max > v_min))
			{
				error = "face has non-finite, reversed, or zero-area UV bounds";
				return false;
			}
			const double u_span = u_max - u_min;
			const double v_span = v_max - v_min;
			// Surface parameters are dimensionless and arbitrarily scaled/translated.
			// Precision::PConfusion is therefore not a valid lower bound here. Bound the
			// arithmetic uncertainty by a small number of representable endpoint steps
			// plus span-relative roundoff; unlike epsilon*abs(parameter), this does not
			// pretend every operation loses hundreds of bits merely because a tiny chart
			// happens to live near parameter 1 (or another large offset).
			const double endpoint_spacing = std::max({
				representable_spacing(u_min), representable_spacing(u_max),
				representable_spacing(v_min), representable_spacing(v_max)});
			const double span_roundoff = 32.0
				* std::numeric_limits<double>::epsilon() * std::max(u_span, v_span);
			const double arithmetic_uv_tolerance =
				8.0 * endpoint_spacing + span_roundoff;
			BRepAdaptor_Surface adaptor(face, true);
			const double u_resolution = adaptor.UResolution(chart.position_tolerance);
			const double v_resolution = adaptor.VResolution(chart.position_tolerance);
			// UvConstraintMesh presently has one isotropic coordinate tolerance. It
			// must be safe in both chart axes, so the more sensitive (smaller positive
			// resolution) axis is authoritative. Taking max here would let an insensitive
			// axis collapse distinct vertices in the sensitive one.
			double conservative_resolution = std::numeric_limits<double>::infinity();
			if (finite(u_resolution) && u_resolution > 0.0)
				conservative_resolution = std::min(conservative_resolution, u_resolution);
			if (finite(v_resolution) && v_resolution > 0.0)
				conservative_resolution = std::min(conservative_resolution, v_resolution);
			if (!finite(conservative_resolution)) conservative_resolution = 0.0;
			// Arithmetic uncertainty and the native surface-resolution certificate are
			// independent error sources, so they add rather than one hiding the other.
			chart.uv_tolerance = arithmetic_uv_tolerance + conservative_resolution;
			if (!finite(chart.position_tolerance) || !(chart.position_tolerance > 0.0)
				|| !finite(chart.uv_tolerance) || !(chart.uv_tolerance > 0.0))
			{
				error = "face tolerance could not be mapped to a finite UV tolerance";
				return false;
			}
			constexpr double maximum_tolerance_fraction = 1.0 / 64.0;
			if (!(chart.uv_tolerance < maximum_tolerance_fraction * u_span)
				|| !(chart.uv_tolerance < maximum_tolerance_fraction * v_span))
			{
				std::ostringstream message;
				message << "face UV chart is numerically unresolved: coordinate tolerance "
					<< chart.uv_tolerance << " consumes a material fraction of chart spans "
					<< u_span << " x " << v_span;
				error = message.str();
				return false;
			}
			chart.u_periodic = chart.surface->IsUPeriodic();
			chart.v_periodic = chart.surface->IsVPeriodic();
			if (chart.u_periodic)
			{
				chart.u_period = chart.surface->UPeriod();
				if (!finite(chart.u_period) || !(chart.u_period > chart.uv_tolerance))
				{
					error = "face reports an invalid U period";
					return false;
				}
			}
			if (chart.v_periodic)
			{
				chart.v_period = chart.surface->VPeriod();
				if (!finite(chart.v_period) || !(chart.v_period > chart.uv_tolerance))
				{
					error = "face reports an invalid V period";
					return false;
				}
			}
			return true;
		}

		bool periodic_alias(const SurfaceChart& chart, const UvPoint& a,
			const UvPoint& b)
		{
			bool shifted = false;
			auto axis_matches = [&](double first, double second, bool periodic,
				double period)
			{
				const double delta = second - first;
				if (!periodic) return std::abs(delta) <= chart.uv_tolerance;
				const double turns = std::round(delta / period);
				if (!finite(turns)
					|| std::abs(delta - turns * period) > chart.uv_tolerance) return false;
				shifted = shifted || std::abs(turns) >= 1.0;
				return true;
			};
			return axis_matches(a.u, b.u, chart.u_periodic, chart.u_period)
				&& axis_matches(a.v, b.v, chart.v_periodic, chart.v_period)
				&& shifted;
		}

		bool preflight_constraints(const SurfaceChart& chart,
			const occ::handle<Poly_Triangulation>& triangulation,
			const TopLoc_Location& triangulation_location,
			std::span<const OcctFaceConstraintPolyline> constraints,
			std::map<std::uint32_t, CanonicalRecord>& canonical,
			PreparedConstraints& prepared,
			std::string& error)
		{
			canonical.clear();
			prepared = {};
			prepared.samples.resize(constraints.size());
			struct Occurrence
			{
				std::size_t polyline = 0;
				std::size_t sample = 0;
				UvPoint uv;
				std::uint64_t branch =
					OcctFaceConstraintPolyline::no_chart_branch_id;
				std::uint32_t source_node_index =
					OcctFaceConstraintSample::no_source_node_index;
				std::uint32_t source_boundary_occurrence_id =
					std::numeric_limits<std::uint32_t>::max();
				std::uint64_t source_contact_atom_id =
					std::numeric_limits<std::uint64_t>::max();
				std::uint32_t source_atom_sample_index =
					std::numeric_limits<std::uint32_t>::max();
			};
			std::map<std::uint32_t, std::vector<Occurrence>> physical_occurrences;
			std::map<std::uint32_t, std::uint32_t> physical_id_by_source_node;
			std::set<std::uint32_t> used_internal_ids;
			if (triangulation.IsNull())
			{
				error = "face has no current Poly_Triangulation";
				return false;
			}

			// Validate every geometric/canonical certificate first. This pass is kept
			// independent of UV occurrence allocation so caller order cannot change which
			// periodic branch receives a temporary planar identity.
			for (std::size_t polyline_index = 0;
				polyline_index < constraints.size(); ++polyline_index)
			{
				const OcctFaceConstraintPolyline& polyline = constraints[polyline_index];
				if (polyline.samples.size() < 2)
				{
					error = "constraint polyline " + std::to_string(polyline_index)
						+ " has fewer than two samples";
					return false;
				}
				prepared.samples[polyline_index].resize(polyline.samples.size());
				for (std::size_t sample_index = 0;
					sample_index < polyline.samples.size(); ++sample_index)
				{
					const OcctFaceConstraintSample& sample = polyline.samples[sample_index];
					if (!finite(sample.uv))
					{
						error = "constraint polyline " + std::to_string(polyline_index)
							+ " sample " + std::to_string(sample_index)
							+ " has non-finite UV coordinates";
						return false;
					}
					const gp_Pnt evaluated = evaluate_world(chart, sample.uv);
					if (!finite(coordinates(evaluated)))
					{
						error = "constraint polyline " + std::to_string(polyline_index)
							+ " sample " + std::to_string(sample_index)
							+ " has a non-finite exact-surface evaluation";
						return false;
					}
					if (sample.topology_id == UvVertex::no_topology_id) continue;
					if (!finite(sample.canonical_world_position)
						|| !finite(sample.canonical_position_tolerance)
						|| sample.canonical_position_tolerance < 0.0)
					{
						error = "constraint polyline " + std::to_string(polyline_index)
							+ " sample " + std::to_string(sample_index)
							+ " has an invalid canonical position or tolerance";
						return false;
					}
					const auto found = canonical.find(sample.topology_id);
					if (found == canonical.end())
						canonical.emplace(sample.topology_id, CanonicalRecord{
							sample.canonical_world_position,
							sample.canonical_position_tolerance});
					else
					{
						if (!same_position(found->second.position,
							sample.canonical_world_position))
						{
							error = "topology ID " + std::to_string(sample.topology_id)
								+ " has conflicting canonical world positions";
							return false;
						}
						found->second.tolerance = std::max(found->second.tolerance,
							sample.canonical_position_tolerance);
					}
					const gp_Pnt canonical_point = point(sample.canonical_world_position);
					const double distance = evaluated.Distance(canonical_point);
					const double tolerance = chart.position_tolerance
						+ sample.canonical_position_tolerance
						+ numerical_position_tolerance(evaluated, canonical_point);
					if (!finite(distance) || distance > tolerance)
					{
						std::ostringstream message;
						message << "constraint polyline " << polyline_index << " sample "
							<< sample_index << " canonical position misses the exact located"
							<< " face surface by " << distance << " model units (limit "
							<< tolerance << ")";
						error = message.str();
						return false;
					}
					if (sample.source_node_index
						!= OcctFaceConstraintSample::no_source_node_index)
					{
						if (sample.source_node_index
							>= static_cast<std::uint32_t>(triangulation->NbNodes()))
						{
							error = "constraint sample names an out-of-range exact source mesh node";
							return false;
						}
						gp_Pnt source_node = triangulation->Node(
							static_cast<Standard_Integer>(sample.source_node_index + 1));
						source_node.Transform(triangulation_location.Transformation());
						const gp_Pnt2d source_node_parameter = triangulation->UVNode(
							static_cast<Standard_Integer>(sample.source_node_index + 1));
						const double source_uv_distance = std::hypot(
							source_node_parameter.X() - sample.uv.u,
							source_node_parameter.Y() - sample.uv.v);
						const UvPoint source_node_uv{source_node_parameter.X(),
							source_node_parameter.Y()};
						// Exact source-node provenance is intentionally stronger than a
						// near-UV lookup: contact construction can carry a slightly displaced
						// parameter while its independent 3-D certificate still identifies the
						// authoritative source node.  The exception is a genuine periodic chart
						// alias.  Reusing the node from U=0 at U=period would erase the second
						// seam occurrence instead of creating its required face-local alias.
						const bool distinct_periodic_occurrence =
							finite(source_node_uv) && finite(source_uv_distance)
							&& source_uv_distance > chart.uv_tolerance
							&& periodic_alias(chart, sample.uv, source_node_uv);
						if (!finite(source_node_parameter.X())
							|| !finite(source_node_parameter.Y())
							|| !finite(source_uv_distance)
							|| distinct_periodic_occurrence)
						{
							std::ostringstream message;
							message.precision(17);
							message << "constraint sample's exact source mesh node belongs to a"
								" different UV chart occurrence: sample [" << sample.uv.u << ','
								<< sample.uv.v << "] node [" << source_node_parameter.X() << ','
								<< source_node_parameter.Y() << "] distance "
								<< source_uv_distance << " (limit " << chart.uv_tolerance << ')';
							message << " provenance [occurrence "
								<< sample.source_boundary_occurrence_id << ", atom "
								<< sample.source_contact_atom_id << ", atom sample "
								<< sample.source_atom_sample_index << ", source node "
								<< sample.source_node_index << ']';
							error = message.str();
							return false;
						}
						const double node_distance = source_node.Distance(canonical_point);
						const double node_tolerance = chart.position_tolerance
							+ sample.canonical_position_tolerance
							+ numerical_position_tolerance(source_node, canonical_point);
						if (!finite(coordinates(source_node)) || !finite(node_distance)
							|| node_distance > node_tolerance)
						{
							std::ostringstream message;
							message << "constraint sample's exact source mesh node misses its canonical"
								" position by " << node_distance << " model units (limit "
								<< node_tolerance << ')';
							error = message.str();
							return false;
						}
						auto [claimed, inserted] = physical_id_by_source_node.emplace(
							sample.source_node_index, sample.topology_id);
						if (!inserted && claimed->second != sample.topology_id)
						{
							error = "one exact source mesh node is claimed by distinct physical topology IDs";
							return false;
						}
					}
					used_internal_ids.insert(sample.topology_id);
					physical_occurrences[sample.topology_id].push_back({polyline_index,
						sample_index, sample.uv, polyline.chart_branch_id,
						sample.source_node_index, sample.source_boundary_occurrence_id,
						sample.source_contact_atom_id,
						sample.source_atom_sample_index});
				}
			}

			std::uint64_t next_temporary =
				static_cast<std::uint64_t>(UvVertex::no_topology_id) - 1u;
			auto allocate_temporary = [&](std::uint32_t output_topology)
				-> std::optional<std::uint32_t>
			{
				while (next_temporary > 0u
					&& used_internal_ids.contains(static_cast<std::uint32_t>(next_temporary)))
					--next_temporary;
				if (used_internal_ids.contains(static_cast<std::uint32_t>(next_temporary)))
					return std::nullopt;
				const std::uint32_t internal = static_cast<std::uint32_t>(next_temporary);
				used_internal_ids.insert(internal);
				prepared.output_topology_by_internal.emplace(internal, output_topology);
				if (next_temporary > 0u) --next_temporary;
				return internal;
			};
			struct AssignedUv
			{
				UvPoint uv;
				std::uint32_t internal_id = UvVertex::no_topology_id;
			};
			std::vector<AssignedUv> assigned_uv;

			// Allocate one planar identity per distinct UV occurrence of a physical
			// topology vertex. Multiple occurrences are accepted only as explicitly named,
			// distinct branches of a periodic chart.
			for (auto& [physical_id, occurrences] : physical_occurrences)
			{
				struct Cluster
				{
					UvPoint uv;
					std::vector<std::size_t> occurrence_indices;
					std::set<std::uint64_t> branches;
					std::uint32_t internal_id = UvVertex::no_topology_id;
					std::optional<std::uint32_t> source_node_index;
					std::optional<std::size_t> source_node_polyline;
					std::optional<std::size_t> source_node_sample;
					std::optional<Occurrence> source_node_occurrence;
				};
				std::vector<Cluster> clusters;
				for (std::size_t occurrence_index = 0;
					occurrence_index < occurrences.size(); ++occurrence_index)
				{
					const Occurrence& occurrence = occurrences[occurrence_index];
					std::optional<std::size_t> matching_cluster;
					for (std::size_t cluster = 0; cluster < clusters.size(); ++cluster)
						if (std::hypot(occurrence.uv.u - clusters[cluster].uv.u,
							occurrence.uv.v - clusters[cluster].uv.v) <= chart.uv_tolerance)
						{
							if (matching_cluster)
							{
								error = "topology ID " + std::to_string(physical_id)
									+ " lies within UV tolerance of multiple chart occurrences";
								return false;
							}
							matching_cluster = cluster;
						}
					if (!matching_cluster)
					{
						clusters.push_back({occurrence.uv, {}, {}});
						matching_cluster = clusters.size() - 1;
					}
					Cluster& cluster = clusters[*matching_cluster];
					if (occurrence.source_node_index
						!= OcctFaceConstraintSample::no_source_node_index)
					{
						if (cluster.source_node_index
							&& *cluster.source_node_index != occurrence.source_node_index)
						{
							const Occurrence& existing = *cluster.source_node_occurrence;
							const gp_Pnt2d existing_node_uv = triangulation->UVNode(
								static_cast<Standard_Integer>(*cluster.source_node_index + 1));
							const gp_Pnt2d conflicting_node_uv = triangulation->UVNode(
								static_cast<Standard_Integer>(occurrence.source_node_index + 1));
							std::ostringstream message;
							message.precision(17);
							message << "one physical topology sample/chart branch names multiple exact"
								<< " source mesh nodes: topology " << physical_id
								<< " branch " << occurrence.branch << " uv ["
								<< occurrence.uv.u << ',' << occurrence.uv.v << "] existing node "
								<< *cluster.source_node_index << " from polyline "
								<< cluster.source_node_polyline.value_or(
									std::numeric_limits<std::size_t>::max()) << " sample "
								<< cluster.source_node_sample.value_or(
									std::numeric_limits<std::size_t>::max()) << " constraint "
								<< constraints[cluster.source_node_polyline.value_or(0)].constraint_id
								<< " provenance [occurrence "
								<< existing.source_boundary_occurrence_id << ", atom "
								<< existing.source_contact_atom_id << ", atom sample "
								<< existing.source_atom_sample_index << ", node uv "
								<< existing_node_uv.X() << ',' << existing_node_uv.Y() << ']'
								<< "; conflicting node " << occurrence.source_node_index
								<< " from polyline " << occurrence.polyline << " sample "
								<< occurrence.sample << " constraint "
								<< constraints[occurrence.polyline].constraint_id
								<< " provenance [occurrence "
								<< occurrence.source_boundary_occurrence_id << ", atom "
								<< occurrence.source_contact_atom_id << ", atom sample "
								<< occurrence.source_atom_sample_index << ", node uv "
								<< conflicting_node_uv.X() << ',' << conflicting_node_uv.Y() << ']';
							error = message.str();
							return false;
						}
						cluster.source_node_index = occurrence.source_node_index;
						cluster.source_node_polyline = occurrence.polyline;
						cluster.source_node_sample = occurrence.sample;
						cluster.source_node_occurrence = occurrence;
					}
					cluster.occurrence_indices.push_back(occurrence_index);
					cluster.branches.insert(occurrence.branch);
				}
				if (clusters.size() > 1)
				{
					std::set<std::uint64_t> unique_branches;
					for (const Cluster& cluster : clusters)
					{
						if (cluster.branches.size() != 1
							|| cluster.branches.contains(
								OcctFaceConstraintPolyline::no_chart_branch_id))
						{
							error = "topology ID " + std::to_string(physical_id)
								+ " occurs at multiple UV locations without one explicit"
								+ " chart branch identity per occurrence";
							return false;
						}
						if (!unique_branches.insert(*cluster.branches.begin()).second)
						{
							error = "topology ID " + std::to_string(physical_id)
								+ " reuses one chart branch identity at multiple UV locations";
							return false;
						}
					}
					for (std::size_t a = 0; a < clusters.size(); ++a)
						for (std::size_t b = a + 1; b < clusters.size(); ++b)
							if (!periodic_alias(chart, clusters[a].uv, clusters[b].uv))
							{
								error = "topology ID " + std::to_string(physical_id)
									+ " has explicitly branched UV locations which are not"
									+ " periodic aliases";
								return false;
							}
				}
				for (std::size_t cluster_index = 0;
					cluster_index < clusters.size(); ++cluster_index)
				{
					Cluster& cluster = clusters[cluster_index];
					if (cluster_index == 0) cluster.internal_id = physical_id;
					else
					{
						const std::optional<std::uint32_t> allocated =
							allocate_temporary(physical_id);
						if (!allocated)
						{
							error = "temporary periodic-branch topology ID space is exhausted";
							return false;
						}
						cluster.internal_id = *allocated;
					}
					for (const AssignedUv& assigned : assigned_uv)
						if (std::hypot(cluster.uv.u - assigned.uv.u,
							cluster.uv.v - assigned.uv.v) <= chart.uv_tolerance
							&& assigned.internal_id != cluster.internal_id)
						{
							error = "distinct physical topology IDs occupy one UV vertex";
							return false;
						}
					assigned_uv.push_back({cluster.uv, cluster.internal_id});
					for (std::size_t occurrence_index : cluster.occurrence_indices)
					{
						const Occurrence& occurrence = occurrences[occurrence_index];
						prepared.samples[occurrence.polyline][occurrence.sample] =
							{occurrence.uv, cluster.internal_id,
								cluster.source_node_index.value_or(
									UvConstraintSample::no_existing_vertex),
								constraints[occurrence.polyline].chart_branch_id
									!= OcctFaceConstraintPolyline::no_chart_branch_id};
					}
				}
			}

			// Sentinel samples still need temporary identities while UvConstraintMesh
			// inserts/splits them. Reuse an authoritative chart vertex when one already
			// exists at that UV, otherwise translate the temporary ID back to sentinel.
			for (std::size_t polyline_index = 0;
				polyline_index < constraints.size(); ++polyline_index)
				for (std::size_t sample_index = 0;
					sample_index < constraints[polyline_index].samples.size(); ++sample_index)
				{
					const OcctFaceConstraintSample& sample =
						constraints[polyline_index].samples[sample_index];
					if (sample.topology_id != UvVertex::no_topology_id) continue;
					if (sample.source_node_index
						!= OcctFaceConstraintSample::no_source_node_index)
					{
						error = "an unclaimed constraint sample cannot name an exact source mesh node";
						return false;
					}
					std::vector<std::uint32_t> matching_ids;
					for (const AssignedUv& assigned : assigned_uv)
						if (std::hypot(sample.uv.u - assigned.uv.u,
							sample.uv.v - assigned.uv.v) <= chart.uv_tolerance
							&& std::find(matching_ids.begin(), matching_ids.end(),
								assigned.internal_id) == matching_ids.end())
							matching_ids.push_back(assigned.internal_id);
					if (matching_ids.size() > 1)
					{
						error = "unclaimed constraint sample lies within UV tolerance of"
							" multiple chart vertices";
						return false;
					}
					std::uint32_t internal = UvVertex::no_topology_id;
					if (!matching_ids.empty()) internal = matching_ids.front();
					else
					{
						const std::optional<std::uint32_t> allocated =
							allocate_temporary(UvVertex::no_topology_id);
						if (!allocated)
						{
							error = "temporary face-local topology ID space is exhausted";
							return false;
						}
						internal = *allocated;
						assigned_uv.push_back({sample.uv, internal});
					}
					prepared.samples[polyline_index][sample_index] = {sample.uv, internal,
						UvConstraintSample::no_existing_vertex,
						constraints[polyline_index].chart_branch_id
							!= OcctFaceConstraintPolyline::no_chart_branch_id};
				}
			return true;
		}

		bool import_uv_mesh(const SurfaceChart& chart,
			const occ::handle<Poly_Triangulation>& triangulation,
			const TopLoc_Location& triangulation_location,
			std::optional<UvConstraintMesh>& mesh, std::size_t& source_node_count,
			std::string& error)
		{
			if (triangulation.IsNull())
			{
				error = "face has no current Poly_Triangulation";
				return false;
			}
			if (!finite_transform(triangulation_location.Transformation()))
			{
				error = "triangulation has an invalid location transform";
				return false;
			}
			if (!triangulation->HasUVNodes())
			{
				error = "Poly_Triangulation has no UV nodes";
				return false;
			}
			if (triangulation->NbNodes() < 3 || triangulation->NbTriangles() < 1)
			{
				error = "Poly_Triangulation has fewer than three nodes or no triangles";
				return false;
			}
			if (static_cast<std::uint64_t>(triangulation->NbNodes())
				>= static_cast<std::uint64_t>(UvVertex::no_topology_id))
			{
				error = "Poly_Triangulation exceeds uint32 vertex-index capacity";
				return false;
			}

			std::vector<UvVertex> vertices;
			vertices.reserve(static_cast<std::size_t>(triangulation->NbNodes()));
			for (int node = 1; node <= triangulation->NbNodes(); ++node)
			{
				const gp_Pnt2d parameter = triangulation->UVNode(node);
				const UvPoint uv{parameter.X(), parameter.Y()};
				if (!finite(uv))
				{
					error = "Poly_Triangulation UV node " + std::to_string(node)
						+ " is non-finite";
					return false;
				}
				gp_Pnt tessellated = triangulation->Node(node);
				tessellated.Transform(triangulation_location.Transformation());
				const gp_Pnt evaluated = evaluate_world(chart, uv);
				const double distance = tessellated.Distance(evaluated);
				const double tolerance = chart.position_tolerance
					+ numerical_position_tolerance(tessellated, evaluated);
				if (!finite(distance) || distance > tolerance)
				{
					std::ostringstream message;
					message << "Poly_Triangulation node " << node
						<< " does not round-trip through its face UV and locations: residual "
						<< distance << " model units, limit " << tolerance;
					error = message.str();
					return false;
				}
				vertices.push_back({uv, UvVertex::no_topology_id});
			}

			std::vector<UvTriangle> triangles;
			triangles.reserve(static_cast<std::size_t>(triangulation->NbTriangles()));
			for (int triangle_index = 1;
				triangle_index <= triangulation->NbTriangles(); ++triangle_index)
			{
				int a = 0, b = 0, c = 0;
				triangulation->Triangle(triangle_index).Get(a, b, c);
				if (a < 1 || b < 1 || c < 1 || a > triangulation->NbNodes()
					|| b > triangulation->NbNodes() || c > triangulation->NbNodes()
					|| a == b || b == c || c == a)
				{
					error = "Poly_Triangulation triangle " + std::to_string(triangle_index)
						+ " has invalid or repeated node indices";
					return false;
				}
				UvTriangle triangle{{static_cast<std::uint32_t>(a - 1),
					static_cast<std::uint32_t>(b - 1),
					static_cast<std::uint32_t>(c - 1)}};
				const long double area2 = oriented_area2(vertices[triangle.vertices[0]].uv,
					vertices[triangle.vertices[1]].uv,
					vertices[triangle.vertices[2]].uv);
				if (area2 == 0.0L)
				{
					error = "Poly_Triangulation triangle " + std::to_string(triangle_index)
						+ " is degenerate in UV";
					return false;
				}
				if (area2 < 0.0L)
					std::swap(triangle.vertices[1], triangle.vertices[2]);
				triangles.push_back(triangle);
			}

			UvConstraintMesh candidate(std::move(vertices), std::move(triangles),
				UvConstraintOptions{chart.uv_tolerance});
			std::string validation_error;
			if (!candidate.validate(&validation_error))
			{
				error = "invalid Poly_Triangulation UV topology: " + validation_error;
				return false;
			}
			mesh.emplace(std::move(candidate));
			source_node_count = static_cast<std::size_t>(triangulation->NbNodes());
			return true;
		}

		bool export_mesh(const TopoDS_Face& face, const SurfaceChart& chart,
			const UvConstraintMesh& mesh,
			const std::map<std::uint32_t, CanonicalRecord>& canonical,
			const std::map<std::uint32_t, std::uint32_t>&
				output_topology_by_internal,
			std::size_t source_node_count, OcctConformedFaceMesh& output,
			std::string& error)
		{
			OcctConformedFaceMesh candidate;
			candidate.face_reversed = face.Orientation() == TopAbs_REVERSED;
			candidate.source_node_count = source_node_count;
			candidate.uv.reserve(mesh.vertices().size());
			candidate.world_positions.reserve(mesh.vertices().size());
			candidate.topology_ids.reserve(mesh.vertices().size());
			for (std::size_t vertex_index = 0;
				vertex_index < mesh.vertices().size(); ++vertex_index)
			{
				const UvVertex& vertex = mesh.vertices()[vertex_index];
				const gp_Pnt evaluated = evaluate_world(chart, vertex.uv);
				if (!finite(coordinates(evaluated)))
				{
					error = "exact surface evaluation is non-finite at output vertex "
						+ std::to_string(vertex_index);
					return false;
				}
				std::array<double, 3> world = coordinates(evaluated);
				std::uint32_t output_topology_id = vertex.topology_id;
				const auto translated =
					output_topology_by_internal.find(vertex.topology_id);
				if (translated != output_topology_by_internal.end())
					output_topology_id = translated->second;
				if (output_topology_id != UvVertex::no_topology_id)
				{
					const auto found = canonical.find(output_topology_id);
					if (found == canonical.end())
					{
						error = "output vertex claims unknown topology ID "
							+ std::to_string(output_topology_id);
						return false;
					}
					const gp_Pnt canonical_point = point(found->second.position);
					const double distance = evaluated.Distance(canonical_point);
					const double tolerance = chart.position_tolerance
						+ found->second.tolerance
						+ numerical_position_tolerance(evaluated, canonical_point);
					if (!finite(distance) || distance > tolerance)
					{
						std::ostringstream message;
						message << "output topology vertex " << output_topology_id
							<< " failed its exact-surface round trip: residual " << distance
							<< " model units, limit " << tolerance;
						error = message.str();
						return false;
					}
					world = found->second.position;
				}
				candidate.uv.push_back(vertex.uv);
				candidate.world_positions.push_back(world);
				candidate.topology_ids.push_back(output_topology_id);
			}
			candidate.triangles.reserve(mesh.triangles().size());
			for (const UvTriangle& triangle : mesh.triangles())
				candidate.triangles.push_back(triangle.vertices);
			const std::vector<UvEdgeView> edges = mesh.edges();
			candidate.edges.reserve(edges.size());
			for (const UvEdgeView& edge : edges)
				candidate.edges.push_back({edge.vertices, edge.boundary,
					edge.constraint_ids});
			output = std::move(candidate);
			return true;
		}
	}

	bool conform_occt_face_triangulation(const TopoDS_Face& face,
		const occ::handle<Poly_Triangulation>& triangulation,
		const TopLoc_Location& triangulation_location,
		std::span<const OcctFaceConstraintPolyline> constraints,
		OcctConformedFaceMesh& output, std::string& error,
		OcctFaceConformerOptions options)
	{
		error.clear();
		if (face.IsNull())
		{
			error = "cannot conform a null face";
			return false;
		}
		if (face.Orientation() != TopAbs_FORWARD
			&& face.Orientation() != TopAbs_REVERSED)
		{
			error = "face orientation must be FORWARD or REVERSED";
			return false;
		}
		if (!finite(options.geometric_tolerance)
			|| options.geometric_tolerance < 0.0)
		{
			error = "face conformer geometric tolerance is invalid";
			return false;
		}

		try
		{
			SurfaceChart chart;
			if (!initialize_chart(face, options, chart, error)) return false;
			std::map<std::uint32_t, CanonicalRecord> canonical;
			PreparedConstraints prepared;
			if (!preflight_constraints(chart, triangulation, triangulation_location,
				constraints, canonical, prepared, error))
				return false;

			std::optional<UvConstraintMesh> mesh;
			std::size_t source_node_count = 0;
			if (!import_uv_mesh(chart, triangulation, triangulation_location,
				mesh, source_node_count, error)) return false;
			for (std::size_t polyline_index = 0;
				polyline_index < constraints.size(); ++polyline_index)
			{
				const OcctFaceConstraintPolyline& polyline = constraints[polyline_index];
				std::string insertion_error;
				if (!mesh->insert_polyline(prepared.samples[polyline_index],
					polyline.constraint_id,
					&insertion_error))
				{
					std::ostringstream report;
					report << "constraint polyline " << polyline_index << " (ID "
						<< polyline.constraint_id << ") was rejected: " << insertion_error
						<< "; exact sample provenance";
					for (std::size_t sample_index = 0;
						sample_index < polyline.samples.size(); ++sample_index)
					{
						const auto& sample = polyline.samples[sample_index];
						report << " [" << sample_index << ": topology " << sample.topology_id
							<< ", node " << sample.source_node_index << ", occurrence "
							<< sample.source_boundary_occurrence_id << ", atom "
							<< sample.source_contact_atom_id << ", atom-sample "
							<< sample.source_atom_sample_index << ']';
					}
					error = report.str();
					return false;
				}
			}
			std::string validation_error;
			if (!mesh->validate(&validation_error))
			{
				error = "conformed UV mesh is invalid: " + validation_error;
				return false;
			}

			OcctConformedFaceMesh candidate;
			if (!export_mesh(face, chart, *mesh, canonical,
				prepared.output_topology_by_internal, source_node_count,
				candidate, error)) return false;
			output = std::move(candidate);
			return true;
		}
		catch (const Standard_Failure& failure)
		{
			const char* message = failure.GetMessageString();
			error = std::string("OpenCascade face-conforming failure: ")
				+ (message ? message : "unknown Standard_Failure");
			return false;
		}
		catch (const std::exception& exception)
		{
			error = std::string("face-conforming failure: ") + exception.what();
			return false;
		}
	}
}
