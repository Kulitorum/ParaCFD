// occt_contact_atomizer.cpp -- interval-only CAD contact atomization.
#include "core/geometry/occt_contact_atomizer.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Curve.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

namespace paracfd::core
{
	namespace
	{
		constexpr std::uint32_t no_id = std::numeric_limits<std::uint32_t>::max();
		constexpr std::uint64_t no_id64 = std::numeric_limits<std::uint64_t>::max();

		bool finite(double value) { return std::isfinite(value); }

		bool finite(const std::array<double, 3>& point)
		{
			return finite(point[0]) && finite(point[1]) && finite(point[2]);
		}

		gp_Pnt as_point(const std::array<double, 3>& point)
		{
			return {point[0], point[1], point[2]};
		}

		std::array<double, 3> as_array(const gp_Pnt& point)
		{
			return {{point.X(), point.Y(), point.Z()}};
		}

		bool lexicographically_less(const gp_Pnt& a, const gp_Pnt& b)
		{
			if (a.X() != b.X()) return a.X() < b.X();
			if (a.Y() != b.Y()) return a.Y() < b.Y();
			return a.Z() < b.Z();
		}

		std::string edge_label(std::uint32_t edge)
		{
			return "source edge " + std::to_string(edge);
		}

		double parameter_roundoff(double first, double last)
		{
			const auto ulp = [](double value)
			{
				return std::abs(std::nextafter(value,
					std::numeric_limits<double>::infinity()) - value);
			};
			const double span = std::abs(last - first);
			return std::max({8.0 * ulp(first), 8.0 * ulp(last),
				128.0 * std::numeric_limits<double>::epsilon() * span});
		}

		bool resolved_parameter_tolerance(const BRepAdaptor_Curve& curve,
			double first, double last, double linear_tolerance,
			double& tolerance, std::string& error)
		{
			const double span = std::abs(last - first);
			double resolution = 0.0;
			try { resolution = std::abs(curve.Resolution(linear_tolerance)); }
			catch (...) { resolution = 0.0; }
			if (!finite(resolution)) resolution = 0.0;
			tolerance = std::max(parameter_roundoff(first, last), resolution);
			if (finite(tolerance) && tolerance >= 0.0 && span > 0.0
				&& tolerance < 1.0e-3 * span) return true;
			std::ostringstream report;
			report.precision(17);
			report << "parameter tolerance " << tolerance
				<< " is unresolved relative to finite curve span " << span;
			error = report.str();
			return false;
		}

		bool valid_boundary_orientation(std::int8_t orientation)
		{
			return orientation == static_cast<std::int8_t>(TopAbs_FORWARD)
				|| orientation == static_cast<std::int8_t>(TopAbs_REVERSED);
		}

		struct DisjointSet
		{
			explicit DisjointSet(std::size_t count) : parent(count), rank(count, 0)
			{
				std::iota(parent.begin(), parent.end(), std::size_t{0});
			}

			std::size_t find(std::size_t value)
			{
				if (parent[value] != value) parent[value] = find(parent[value]);
				return parent[value];
			}

			void merge(std::size_t a, std::size_t b)
			{
				a = find(a);
				b = find(b);
				if (a == b) return;
				if (rank[a] < rank[b]) std::swap(a, b);
				parent[b] = a;
				if (rank[a] == rank[b]) ++rank[a];
			}

			std::vector<std::size_t> parent;
			std::vector<std::uint8_t> rank;
		};

		struct EdgeData
		{
			const OcctContactSourceEdge* source = nullptr;
			std::uint32_t id = no_id;
			TopoDS_Edge edge;
			occ::handle<Geom_Curve> curve;
			TopLoc_Location location;
			double first = 0.0;
			double last = 0.0;
			double parameter_tolerance = 0.0;
			double geometric_tolerance = 0.0;
			std::vector<double> cuts;
		};

		gp_Pnt edge_value(const EdgeData& edge, double parameter)
		{
			gp_Pnt point = edge.curve->Value(parameter);
			point.Transform(edge.location.Transformation());
			return point;
		}

		bool inside_parameter(const EdgeData& edge, double parameter)
		{
			return parameter >= edge.first - edge.parameter_tolerance
				&& parameter <= edge.last + edge.parameter_tolerance;
		}

		double clamped_parameter(const EdgeData& edge, double parameter)
		{
			if (std::abs(parameter - edge.first) <= edge.parameter_tolerance) return edge.first;
			if (std::abs(parameter - edge.last) <= edge.parameter_tolerance) return edge.last;
			return std::clamp(parameter, edge.first, edge.last);
		}

		void normalize_cuts(EdgeData& edge)
		{
			for (double& cut : edge.cuts) cut = clamped_parameter(edge, cut);
			std::sort(edge.cuts.begin(), edge.cuts.end());
			std::vector<double> normalized;
			for (double cut : edge.cuts)
			{
				if (normalized.empty() || cut - normalized.back() > edge.parameter_tolerance)
					normalized.push_back(cut);
				else if (cut == edge.last)
					normalized.back() = edge.last;
			}
			if (normalized.empty() || normalized.front() != edge.first)
				normalized.insert(normalized.begin(), edge.first);
			if (normalized.back() != edge.last) normalized.push_back(edge.last);
			edge.cuts = std::move(normalized);
		}

		struct ExactIntervalData
		{
			const OcctContactEdgeFaceInterval* source = nullptr;
			std::size_t edge = 0;
			double begin = 0.0;
			double end = 0.0;
		};

		struct OverlapData
		{
			const OcctContactReciprocalInterval* source = nullptr;
			std::size_t edge_a = 0;
			std::size_t edge_b = 0;
		};

		struct JunctionData
		{
			const OcctContactExactJunction* source = nullptr;
			std::vector<std::pair<std::size_t, double>> incidences;
		};

		struct Piece
		{
			std::size_t edge = 0;
			double begin = 0.0;
			double end = 0.0;
		};

		struct AtomSampleBuild
		{
			double canonical_parameter = 0.0;
			gp_Pnt position;
			double tolerance = 0.0;
			double certificate_tolerance = 0.0;
			std::vector<OcctContactAtomSampleSource> origins;
			std::vector<std::uint64_t> junction_ids;
			// Exact junctions outrank atom breakpoints, which outrank inherited
			// tessellation samples when choosing a bit-identical shared coordinate.
			std::uint8_t canonical_priority = 0;
			OcctContactAtomSampleSource canonical_origin;
		};

		struct AtomBuild
		{
			OcctContactAtom atom;
			std::vector<AtomSampleBuild> sample_metadata;
			// Exact oriented target occurrence retained only during OCCT preprocessing.
			std::map<std::pair<std::uint32_t, std::uint32_t>, TopoDS_Edge> boundary_edges;
		};

		bool parameter_close(const EdgeData& edge, double a, double b)
		{
			return std::abs(a - b) <= edge.parameter_tolerance;
		}

		bool range_contains(double begin, double end, double parameter, double tolerance)
		{
			if (end < begin) std::swap(begin, end);
			return parameter >= begin - tolerance && parameter <= end + tolerance;
		}

		bool range_contains_piece(double begin, double end, const Piece& piece, double tolerance)
		{
			if (end < begin) std::swap(begin, end);
			return piece.begin >= begin - tolerance && piece.end <= end + tolerance;
		}

		bool project_parameter(const EdgeData& edge, const gp_Pnt& world,
			double range_begin, double range_end, double linear_tolerance,
			double& parameter, gp_Pnt& projected, std::string& error)
		{
			if (range_end < range_begin) std::swap(range_begin, range_end);
			gp_Pnt local = world;
			local.Transform(edge.location.Transformation().Inverted());
			GeomAPI_ProjectPointOnCurve projection(local, edge.curve, range_begin, range_end);
			std::vector<std::pair<double, gp_Pnt>> accepted;
			auto append_candidate = [&](double candidate)
			{
				candidate = std::clamp(candidate, range_begin, range_end);
				const gp_Pnt point = edge_value(edge, candidate);
				if (point.SquareDistance(world) <= linear_tolerance * linear_tolerance)
					accepted.emplace_back(candidate, point);
			};
			append_candidate(range_begin);
			append_candidate(range_end);
			for (int candidate = 1; candidate <= projection.NbPoints(); ++candidate)
				append_candidate(projection.Parameter(candidate));
			std::sort(accepted.begin(), accepted.end(), [](const auto& a, const auto& b)
				{ return a.first < b.first; });
			accepted.erase(std::unique(accepted.begin(), accepted.end(), [&](const auto& a,
				const auto& b) { return std::abs(a.first - b.first) <= edge.parameter_tolerance; }),
				accepted.end());
			if (accepted.empty())
			{
				error = "exact bounded target curve has no parameter within the native CAD bound";
				return false;
			}
			if (accepted.size() != 1)
			{
				error = "exact bounded target curve has multiple parameters within the native CAD bound";
				return false;
			}
			parameter = accepted.front().first;
			projected = accepted.front().second;
			return true;
		}

		bool overlap_map_parameter(const EdgeData& from, const EdgeData& to,
			double from_begin, double from_end, double to_begin, double to_end,
			double parameter, double tolerance, double& mapped, std::string& error)
		{
			if (parameter_close(from, parameter, from_begin))
			{
				mapped = to_begin;
				return true;
			}
			if (parameter_close(from, parameter, from_end))
			{
				mapped = to_end;
				return true;
			}
			if (!range_contains(from_begin, from_end, parameter, from.parameter_tolerance))
			{
				error = "internal overlap mapping requested outside its exact source interval";
				return false;
			}
			const gp_Pnt point = edge_value(from, parameter);
			gp_Pnt projected;
			if (!project_parameter(to, point, to_begin, to_end, tolerance, mapped,
				projected, error)) return false;
			if (point.SquareDistance(projected) > tolerance * tolerance)
			{
				error = "caller-supplied reciprocal interval does not map within its native CAD bound";
				return false;
			}
			return true;
		}

		double native_pair_tolerance(const EdgeData& a, const EdgeData& b,
			double additional)
		{
			return a.geometric_tolerance + b.geometric_tolerance + additional
				+ Precision::Confusion();
		}

		bool validate_overlap_endpoints(const EdgeData& a, const EdgeData& b,
			const OcctContactReciprocalInterval& overlap, std::string& error)
		{
			const double tolerance = native_pair_tolerance(a, b, overlap.geometric_tolerance);
			const gp_Pnt a0 = edge_value(a, overlap.parameter_a_begin);
			const gp_Pnt a1 = edge_value(a, overlap.parameter_a_end);
			const gp_Pnt b0 = edge_value(b, overlap.parameter_b_begin);
			const gp_Pnt b1 = edge_value(b, overlap.parameter_b_end);
			if (a0.SquareDistance(b0) > tolerance * tolerance
				|| a1.SquareDistance(b1) > tolerance * tolerance)
			{
				error = "reciprocal interval endpoint correspondence is outside the native CAD bound";
				return false;
			}
			return true;
		}

		std::optional<std::size_t> matching_piece(const std::vector<std::size_t>& pieces,
			const std::vector<Piece>& all, const EdgeData& edge, double begin, double end)
		{
			if (end < begin) std::swap(begin, end);
			for (std::size_t piece_id : pieces)
			{
				const Piece& piece = all[piece_id];
				if (parameter_close(edge, piece.begin, begin)
					&& parameter_close(edge, piece.end, end)) return piece_id;
			}
			return std::nullopt;
		}

		bool occurrence_parameter(const EdgeData& source_edge, double source_parameter,
			const OcctContactBoundaryOccurrence& occurrence, double& target_parameter,
			std::string& error)
		{
			if (occurrence.target_edge.IsNull())
			{
				error = "boundary occurrence is missing its exact oriented TopoDS edge";
				return false;
			}
			TopLoc_Location target_location;
			Standard_Real target_first = 0.0, target_last = 0.0;
			const auto target_curve = BRep_Tool::Curve(occurrence.target_edge,
				target_location, target_first, target_last);
			if (target_curve.IsNull() || !finite(target_first) || !finite(target_last)
				|| !(target_last > target_first))
			{
				error = "boundary occurrence has no finite nondegenerate 3D curve";
				return false;
			}
			const double mapped_first = std::min(occurrence.target_parameter_begin,
				occurrence.target_parameter_end);
			const double mapped_last = std::max(occurrence.target_parameter_begin,
				occurrence.target_parameter_end);
			const double target_tolerance = std::max(0.0,
				BRep_Tool::Tolerance(occurrence.target_edge))
				* std::abs(target_location.Transformation().ScaleFactor());
			const double tolerance = source_edge.geometric_tolerance + target_tolerance
				+ Precision::Confusion();
			BRepAdaptor_Curve target_adaptor(occurrence.target_edge);
			double target_parameter_tolerance = 0.0;
			if (!resolved_parameter_tolerance(target_adaptor, target_first, target_last,
				target_tolerance + Precision::Confusion(), target_parameter_tolerance,
				error))
			{
				error = "boundary occurrence " + error;
				return false;
			}
			if (mapped_first < target_first - target_parameter_tolerance
				|| mapped_last > target_last + target_parameter_tolerance
				|| !(mapped_last > mapped_first)
				|| target_parameter_tolerance >= 1.0e-3 * (mapped_last - mapped_first))
			{
				error = "boundary occurrence target parameter range is unresolved or outside its exact edge";
				return false;
			}
			const gp_Pnt source_world = edge_value(source_edge, source_parameter);
			auto target_value = [&](double parameter)
			{
				gp_Pnt point = target_curve->Value(parameter);
				point.Transform(target_location.Transformation());
				return point;
			};
			if (parameter_close(source_edge, source_parameter,
				occurrence.source_parameter_begin))
				target_parameter = occurrence.target_parameter_begin;
			else if (parameter_close(source_edge, source_parameter,
				occurrence.source_parameter_end))
				target_parameter = occurrence.target_parameter_end;
			else
			{
				EdgeData target;
				target.edge = occurrence.target_edge;
				target.curve = target_curve;
				target.location = target_location;
				target.first = target_first;
				target.last = target_last;
				target.parameter_tolerance = target_parameter_tolerance;
				target.geometric_tolerance = target_tolerance;
				gp_Pnt projected;
				if (!project_parameter(target, source_world, mapped_first, mapped_last,
					tolerance, target_parameter, projected, error))
				{
					error = "could not recover an atom endpoint on its exact boundary occurrence: "
						+ error;
					return false;
				}
			}
			const gp_Pnt target_world = target_value(target_parameter);
			if (source_world.SquareDistance(target_world) > tolerance * tolerance)
			{
				error = "atom endpoint is outside its exact boundary occurrence's native CAD bound";
				return false;
			}
			return true;
		}

		bool append_boundary_slice(const EdgeData& source_edge,
			const OcctContactBoundaryOccurrence& occurrence,
			double atom_begin, double atom_end,
			std::vector<OcctContactAtomBoundaryUse>& output, std::string& error)
		{
			if (occurrence.occurrence_id == OcctContactBoundaryOccurrence::no_occurrence_id)
				return true;
			if (!finite(occurrence.source_parameter_begin)
				|| !finite(occurrence.source_parameter_end)
				|| !finite(occurrence.target_parameter_begin)
				|| !finite(occurrence.target_parameter_end)
				|| std::abs(occurrence.source_parameter_end
					- occurrence.source_parameter_begin) <= source_edge.parameter_tolerance)
			{
				error = "boundary occurrence has an invalid source/target parameter mapping";
				return false;
			}
			if (!range_contains(occurrence.source_parameter_begin,
				occurrence.source_parameter_end, atom_begin, source_edge.parameter_tolerance)
				|| !range_contains(occurrence.source_parameter_begin,
					occurrence.source_parameter_end, atom_end, source_edge.parameter_tolerance))
			{
				error = "boundary occurrence parameter mapping does not cover its contact atom";
				return false;
			}
			double target_begin = 0.0, target_end = 0.0;
			if (!occurrence_parameter(source_edge, atom_begin, occurrence,
				target_begin, error)
				|| !occurrence_parameter(source_edge, atom_end, occurrence,
					target_end, error)) return false;
			OcctContactAtomBoundaryUse use;
			use.source_edge_id = occurrence.target_source_edge_id;
			use.mapping_source_edge_id = source_edge.id;
			use.occurrence_id = occurrence.occurrence_id;
			use.orientation = occurrence.orientation;
			use.source_parameter_begin = atom_begin;
			use.source_parameter_end = atom_end;
			use.target_parameter_begin = target_begin;
			use.target_parameter_end = target_end;
			output.push_back(use);
			return true;
		}

		std::vector<OcctContactBoundaryOccurrence> exact_occurrences(
			const ExactIntervalData& interval)
		{
			std::vector<OcctContactBoundaryOccurrence> output;
			const auto& source = interval.source->interval;
			for (std::size_t occurrence = 0;
				occurrence < source.boundary_occurrence_count; ++occurrence)
			{
				OcctContactBoundaryOccurrence item;
				item.occurrence_id = source.target_boundary_occurrence_ids[occurrence];
				item.target_source_edge_id =
					interval.source->target_boundary_source_edge_ids[occurrence];
				item.orientation = source.target_boundary_orientations[occurrence];
				item.target_edge = interval.source->target_boundary_occurrences[occurrence];
				// target_parameter_begin/end are certified for this emitted clipped
				// interval.  Pair them with this interval's raw source endpoints, not
				// with the enclosing provenance mapping range.
				item.source_parameter_begin = interval.begin;
				item.source_parameter_end = interval.end;
				item.target_parameter_begin = source.target_parameter_begin[occurrence];
				item.target_parameter_end = source.target_parameter_end[occurrence];
				output.push_back(item);
			}
			return output;
		}

		bool boundary_use_less(const OcctContactAtomBoundaryUse& a,
			const OcctContactAtomBoundaryUse& b)
		{
			const auto key = [](const OcctContactAtomBoundaryUse& use)
			{
				// Prefer the target occurrence's own edge mapping over an indirect
				// reciprocal certificate when canonicalizing duplicate provenance.
				return std::tuple(use.occurrence_id,
					use.mapping_source_edge_id != use.source_edge_id,
					use.mapping_source_edge_id, use.orientation, use.source_edge_id,
					use.source_parameter_begin, use.source_parameter_end,
					use.target_parameter_begin, use.target_parameter_end);
			};
			return key(a) < key(b);
		}

		bool source_location_less(const OcctContactAtomSampleSource& a,
			const OcctContactAtomSampleSource& b)
		{
			return std::tie(a.source_edge_id, a.parameter)
				< std::tie(b.source_edge_id, b.parameter);
		}

		bool source_location_equal(const OcctContactAtomSampleSource& a,
			const OcctContactAtomSampleSource& b)
		{
			return a.source_edge_id == b.source_edge_id && a.parameter == b.parameter;
		}

		std::optional<std::size_t> edge_index(const std::map<std::uint32_t, std::size_t>& ids,
			std::uint32_t id)
		{
			const auto found = ids.find(id);
			return found == ids.end() ? std::nullopt
				: std::optional<std::size_t>(found->second);
		}
	}

	bool atomize_occt_contacts(const OcctContactAtomizerInput& input,
		OcctContactAtomization& output, std::string& error)
	{
		error.clear();
		OcctContactAtomization result;
		try
		{
			if (input.source_edges.empty())
			{
				error = "contact atomizer requires at least one source edge";
				return false;
			}

			std::vector<const OcctContactSourceEdge*> ordered_sources;
			ordered_sources.reserve(input.source_edges.size());
			for (const auto& source : input.source_edges) ordered_sources.push_back(&source);
			std::sort(ordered_sources.begin(), ordered_sources.end(), [](const auto* a, const auto* b)
				{ return a->source_edge_id < b->source_edge_id; });

			std::vector<EdgeData> edges;
			std::map<std::uint32_t, std::size_t> edge_by_id;
			for (const OcctContactSourceEdge* source : ordered_sources)
			{
				if (!source || source->source_edge_id == no_id || source->edge.IsNull())
				{
					error = "contact atomizer received a null or unnumbered source edge";
					return false;
				}
				if (!edge_by_id.emplace(source->source_edge_id, edges.size()).second)
				{
					error = "duplicate " + edge_label(source->source_edge_id);
					return false;
				}
				if (source->owner_face_uses.empty())
				{
					error = edge_label(source->source_edge_id) + " has no owner face use";
					return false;
				}
				for (const auto& use : source->owner_face_uses)
					if (use.sector_count < 1 || use.sector_count > 2
						|| (use.location == OcctFaceIntervalLocation::boundary
							&& use.sector_count != 1)
						|| (use.location == OcctFaceIntervalLocation::interior
							&& use.sector_count != 2))
					{
						error = edge_label(source->source_edge_id)
							+ " has inconsistent owner face location/sector metadata";
						return false;
					}

				EdgeData edge;
				edge.source = source;
				edge.id = source->source_edge_id;
				edge.edge = source->edge;
				Standard_Real first = 0.0, last = 0.0;
				edge.curve = BRep_Tool::Curve(edge.edge, edge.location, first, last);
				if (edge.curve.IsNull() || !finite(first) || !finite(last) || !(last > first))
				{
					error = edge_label(edge.id) + " has no finite nondegenerate 3D curve range";
					return false;
				}
				edge.first = first;
				edge.last = last;
				const double scale = std::abs(edge.location.Transformation().ScaleFactor());
				edge.geometric_tolerance = std::max(0.0, BRep_Tool::Tolerance(edge.edge)) * scale;
				BRepAdaptor_Curve adaptor(edge.edge);
				if (!resolved_parameter_tolerance(adaptor, first, last,
					edge.geometric_tolerance + Precision::Confusion(),
					edge.parameter_tolerance, error))
				{
					error = edge_label(edge.id) + " " + error;
					return false;
				}
				for (const auto& use : source->owner_face_uses)
					for (const auto& occurrence : use.boundary_occurrences)
						if (occurrence.occurrence_id
								== OcctContactBoundaryOccurrence::no_occurrence_id
							|| occurrence.target_source_edge_id == no_id
							|| !valid_boundary_orientation(occurrence.orientation)
							|| occurrence.target_edge.IsNull()
							|| !finite(occurrence.source_parameter_begin)
							|| !finite(occurrence.source_parameter_end)
							|| !finite(occurrence.target_parameter_begin)
							|| !finite(occurrence.target_parameter_end)
							|| !range_contains(occurrence.source_parameter_begin,
								occurrence.source_parameter_end, first,
								edge.parameter_tolerance)
							|| !range_contains(occurrence.source_parameter_begin,
								occurrence.source_parameter_end, last,
								edge.parameter_tolerance))
						{
							error = edge_label(edge.id)
								+ " has invalid owner boundary-occurrence metadata";
							return false;
						}
				edge.cuts = {first, last};
				for (const auto& sample : source->samples)
					if (!finite(sample.parameter) || !finite(sample.tolerance)
						|| sample.tolerance < 0.0 || !inside_parameter(edge, sample.parameter))
					{
						error = edge_label(edge.id) + " has an invalid retained source sample";
						return false;
					}
				edges.push_back(std::move(edge));
			}
			for (const EdgeData& edge : edges)
				for (const auto& use : edge.source->owner_face_uses)
					for (const auto& occurrence : use.boundary_occurrences)
					{
						const auto target = edge_index(edge_by_id,
							occurrence.target_source_edge_id);
						if (!target
							|| !edges[*target].edge.IsSame(occurrence.target_edge))
						{
							error = edge_label(edge.id)
								+ " owner occurrence does not match its exact target source edge";
							return false;
						}
					}

			std::vector<ExactIntervalData> intervals;
			intervals.reserve(input.edge_face_intervals.size());
			for (const auto& source : input.edge_face_intervals)
			{
				const auto index = edge_index(edge_by_id, source.source_edge_id);
				if (!index)
				{
					error = "edge/face interval references an unknown source edge";
					return false;
				}
				EdgeData& edge = edges[*index];
				if (!finite(source.source_parameter_first)
					|| !finite(source.source_parameter_last)
					|| !(source.source_parameter_last > source.source_parameter_first)
					|| !parameter_close(edge, source.source_parameter_first, edge.first)
					|| !parameter_close(edge, source.source_parameter_last, edge.last))
				{
					error = edge_label(edge.id)
						+ " edge/face interval uses a different finite source range";
					return false;
				}
				const auto& interval = source.interval;
				if (!finite(interval.begin) || !finite(interval.end)
					|| interval.begin < 0.0 || interval.end > 1.0
					|| !(interval.end > interval.begin)
					|| interval.target_sector_count == 0
					|| !finite(interval.exact_operation_tolerance)
					|| interval.exact_operation_tolerance < 0.0)
				{
					error = edge_label(edge.id) + " has an invalid exact edge/face interval";
					return false;
				}
				if ((interval.location == OcctFaceIntervalLocation::boundary
						&& interval.target_sector_count != 1)
					|| (interval.location == OcctFaceIntervalLocation::interior
						&& interval.target_sector_count != 2)
					|| interval.boundary_occurrence_count > 2)
				{
					error = edge_label(edge.id)
						+ " exact interval has inconsistent location/sector metadata";
					return false;
				}
				const double raw_begin = edge.first + (edge.last - edge.first) * interval.begin;
				const double raw_end = edge.first + (edge.last - edge.first) * interval.end;
				if (!(raw_end > raw_begin + edge.parameter_tolerance))
				{
					error = edge_label(edge.id) + " exact interval is below native resolution";
					return false;
				}
				for (std::size_t occurrence = 0;
					occurrence < interval.boundary_occurrence_count; ++occurrence)
				{
					const double map_begin = interval.target_mapping_source_begin[occurrence];
					const double map_end = interval.target_mapping_source_end[occurrence];
					if (interval.target_boundary_occurrence_ids[occurrence]
							== OcctContactBoundaryOccurrence::no_occurrence_id
						|| source.target_boundary_source_edge_ids[occurrence] == no_id
						|| !edge_by_id.contains(
							source.target_boundary_source_edge_ids[occurrence])
						|| !valid_boundary_orientation(
							interval.target_boundary_orientations[occurrence])
						|| source.target_boundary_occurrences[occurrence].IsNull()
						|| !finite(map_begin) || !finite(map_end)
						|| !finite(interval.target_parameter_begin[occurrence])
						|| !finite(interval.target_parameter_end[occurrence])
						|| std::abs(map_end - map_begin) <= edge.parameter_tolerance
						|| !range_contains(map_begin, map_end, raw_begin,
							edge.parameter_tolerance)
						|| !range_contains(map_begin, map_end, raw_end,
							edge.parameter_tolerance))
					{
						error = edge_label(edge.id)
							+ " exact boundary occurrence has an invalid parameter map";
						return false;
					}
					const auto target_index = edge_index(edge_by_id,
						source.target_boundary_source_edge_ids[occurrence]);
					if (!target_index || !edges[*target_index].edge.IsSame(
						source.target_boundary_occurrences[occurrence]))
					{
						error = edge_label(edge.id)
							+ " exact boundary occurrence does not match its target source edge";
						return false;
					}
				}
				edge.cuts.push_back(raw_begin);
				edge.cuts.push_back(raw_end);
				intervals.push_back({&source, *index, raw_begin, raw_end});
			}

			std::sort(intervals.begin(), intervals.end(), [&](const auto& a, const auto& b)
			{
				return std::tie(edges[a.edge].id, a.begin, a.end, a.source->target_face_id,
					a.source->interval.target_sector_count,
					a.source->interval.target_boundary_occurrence_ids)
					< std::tie(edges[b.edge].id, b.begin, b.end, b.source->target_face_id,
						b.source->interval.target_sector_count,
						b.source->interval.target_boundary_occurrence_ids);
			});

			std::vector<OverlapData> overlaps;
			overlaps.reserve(input.reciprocal_intervals.size());
			for (const auto& source : input.reciprocal_intervals)
			{
				const auto a = edge_index(edge_by_id, source.source_edge_a);
				const auto b = edge_index(edge_by_id, source.source_edge_b);
				if (!a || !b || *a == *b)
				{
					error = "reciprocal interval references an unknown or identical edge pair";
					return false;
				}
				if (!finite(source.parameter_a_begin) || !finite(source.parameter_a_end)
					|| !finite(source.parameter_b_begin) || !finite(source.parameter_b_end)
					|| !finite(source.geometric_tolerance) || source.geometric_tolerance < 0.0
					|| !inside_parameter(edges[*a], source.parameter_a_begin)
					|| !inside_parameter(edges[*a], source.parameter_a_end)
					|| !inside_parameter(edges[*b], source.parameter_b_begin)
					|| !inside_parameter(edges[*b], source.parameter_b_end)
					|| std::abs(source.parameter_a_end - source.parameter_a_begin)
						<= edges[*a].parameter_tolerance
					|| std::abs(source.parameter_b_end - source.parameter_b_begin)
						<= edges[*b].parameter_tolerance
					|| edges[*a].parameter_tolerance >= 1.0e-3
						* std::abs(source.parameter_a_end - source.parameter_a_begin)
					|| edges[*b].parameter_tolerance >= 1.0e-3
						* std::abs(source.parameter_b_end - source.parameter_b_begin))
				{
					error = "reciprocal interval has an invalid parameter range or tolerance";
					return false;
				}
				if (!validate_overlap_endpoints(edges[*a], edges[*b], source, error)) return false;
				edges[*a].cuts.push_back(source.parameter_a_begin);
				edges[*a].cuts.push_back(source.parameter_a_end);
				edges[*b].cuts.push_back(source.parameter_b_begin);
				edges[*b].cuts.push_back(source.parameter_b_end);
				overlaps.push_back({&source, *a, *b});
			}

			std::sort(overlaps.begin(), overlaps.end(), [&](const auto& lhs, const auto& rhs)
			{
				const auto key = [&](const OverlapData& overlap)
				{
					const auto& source = *overlap.source;
					return std::tuple(std::min(source.source_edge_a, source.source_edge_b),
						std::max(source.source_edge_a, source.source_edge_b),
						std::min(source.parameter_a_begin, source.parameter_a_end),
						std::min(source.parameter_b_begin, source.parameter_b_end));
				};
				return key(lhs) < key(rhs);
			});

			std::vector<JunctionData> junctions;
			std::set<std::uint64_t> junction_ids;
			for (const auto& source : input.exact_junctions)
			{
				if (source.junction_id == no_id64
					|| !junction_ids.insert(source.junction_id).second
					|| source.incidences.size() < 2 || !finite(source.canonical_world_position)
					|| !finite(source.geometric_tolerance) || source.geometric_tolerance < 0.0)
				{
					error = "exact junction has an invalid/duplicate ID, point, tolerance or incidence count";
					return false;
				}
				JunctionData junction;
				junction.source = &source;
				for (const auto& incidence : source.incidences)
				{
					const auto index = edge_index(edge_by_id, incidence.source_edge_id);
					if (!index || !finite(incidence.parameter)
						|| !inside_parameter(edges[*index], incidence.parameter))
					{
						error = "exact junction has an invalid source-edge incidence";
						return false;
					}
					const EdgeData& edge = edges[*index];
					const bool duplicate = std::any_of(junction.incidences.begin(),
						junction.incidences.end(), [&](const auto& existing)
						{
							return existing.first == *index && parameter_close(edge,
								existing.second, incidence.parameter);
						});
					if (duplicate)
					{
						error = "exact junction repeats one resolved source-edge parameter";
						return false;
					}
					const double tolerance = edge.geometric_tolerance
						+ source.geometric_tolerance + Precision::Confusion();
					if (edge_value(edge, incidence.parameter).SquareDistance(
							as_point(source.canonical_world_position)) > tolerance * tolerance)
					{
						error = "exact junction point is outside a source edge's native CAD bound";
						return false;
					}
					edges[*index].cuts.push_back(incidence.parameter);
					junction.incidences.emplace_back(*index, incidence.parameter);
				}
				std::sort(junction.incidences.begin(), junction.incidences.end(),
					[&](const auto& a, const auto& b)
					{
						return std::tie(edges[a.first].id, a.second)
							< std::tie(edges[b.first].id, b.second);
					});
				junctions.push_back(std::move(junction));
			}
			std::sort(junctions.begin(), junctions.end(), [](const auto& a, const auto& b)
				{ return a.source->junction_id < b.source->junction_id; });

			for (EdgeData& edge : edges) normalize_cuts(edge);

			// Every split on one certified reciprocal span is also a split on its copy.
			// Iterate to a fixed point because a split can traverse a chain of three or
			// more differently parameterized copies. Mapping occurs only inside supplied
			// overlaps; this is not contact discovery.
			bool split_fixed_point = overlaps.empty();
			for (std::size_t pass = 0; pass <= edges.size() && !split_fixed_point; ++pass)
			{
				std::size_t cut_count_before = 0;
				for (const EdgeData& edge : edges) cut_count_before += edge.cuts.size();
				for (const OverlapData& overlap : overlaps)
				{
					const auto& source = *overlap.source;
					EdgeData& a = edges[overlap.edge_a];
					EdgeData& b = edges[overlap.edge_b];
					const double tolerance = native_pair_tolerance(a, b,
						source.geometric_tolerance);
					const std::vector<double> cuts_a = a.cuts;
					const std::vector<double> cuts_b = b.cuts;
					for (double cut : cuts_a)
						if (range_contains(source.parameter_a_begin, source.parameter_a_end,
							cut, a.parameter_tolerance))
						{
							double mapped = 0.0;
							if (!overlap_map_parameter(a, b, source.parameter_a_begin,
								source.parameter_a_end, source.parameter_b_begin,
								source.parameter_b_end, cut, tolerance, mapped, error)) return false;
							b.cuts.push_back(mapped);
						}
					for (double cut : cuts_b)
						if (range_contains(source.parameter_b_begin, source.parameter_b_end,
							cut, b.parameter_tolerance))
						{
							double mapped = 0.0;
							if (!overlap_map_parameter(b, a, source.parameter_b_begin,
								source.parameter_b_end, source.parameter_a_begin,
								source.parameter_a_end, cut, tolerance, mapped, error)) return false;
							a.cuts.push_back(mapped);
						}
				}
				for (EdgeData& edge : edges) normalize_cuts(edge);
				std::size_t cut_count_after = 0;
				for (const EdgeData& edge : edges) cut_count_after += edge.cuts.size();
				split_fixed_point = cut_count_after == cut_count_before;
			}
			if (!split_fixed_point)
			{
				error = "reciprocal contact split propagation did not reach a deterministic fixed point";
				return false;
			}

			std::vector<Piece> pieces;
			std::vector<std::vector<std::size_t>> pieces_by_edge(edges.size());
			for (std::size_t edge_id = 0; edge_id < edges.size(); ++edge_id)
				for (std::size_t cut = 1; cut < edges[edge_id].cuts.size(); ++cut)
				{
					const double begin = edges[edge_id].cuts[cut - 1];
					const double end = edges[edge_id].cuts[cut];
					if (!(end > begin + edges[edge_id].parameter_tolerance)) continue;
					if (edges[edge_id].parameter_tolerance >= 1.0e-3 * (end - begin))
					{
						std::ostringstream report;
						report.precision(17);
						report << edge_label(edges[edge_id].id)
							<< " atom parameter span is unresolved at native CAD tolerance"
							<< ": begin=" << begin << ", end=" << end
							<< ", span=" << (end - begin)
							<< ", parameter_tolerance="
							<< edges[edge_id].parameter_tolerance
							<< ", edge_range=[" << edges[edge_id].first << ','
							<< edges[edge_id].last << ']';
						error = report.str();
						return false;
					}
					pieces_by_edge[edge_id].push_back(pieces.size());
					pieces.push_back({edge_id, begin, end});
				}
			if (pieces.empty())
			{
				error = "contact atomizer produced no positive-length source interval";
				return false;
			}

			DisjointSet physical(pieces.size());
			for (const OverlapData& overlap : overlaps)
			{
				const auto& source = *overlap.source;
				const EdgeData& a = edges[overlap.edge_a];
				const EdgeData& b = edges[overlap.edge_b];
				const double tolerance = native_pair_tolerance(a, b, source.geometric_tolerance);
				std::set<std::size_t> paired_b;
				for (std::size_t piece_a_id : pieces_by_edge[overlap.edge_a])
				{
					const Piece& piece_a = pieces[piece_a_id];
					if (!range_contains_piece(source.parameter_a_begin, source.parameter_a_end,
						piece_a, a.parameter_tolerance)) continue;
					double b0 = 0.0, b1 = 0.0;
					if (!overlap_map_parameter(a, b, source.parameter_a_begin,
						source.parameter_a_end, source.parameter_b_begin,
						source.parameter_b_end, piece_a.begin, tolerance, b0, error)
						|| !overlap_map_parameter(a, b, source.parameter_a_begin,
							source.parameter_a_end, source.parameter_b_begin,
							source.parameter_b_end, piece_a.end, tolerance, b1, error)) return false;
					const auto piece_b = matching_piece(pieces_by_edge[overlap.edge_b], pieces,
						b, b0, b1);
					if (!piece_b)
					{
						error = "reciprocal contact spans were not split into matching physical atoms";
						return false;
					}
					physical.merge(piece_a_id, *piece_b);
					paired_b.insert(*piece_b);
				}
				for (std::size_t piece_b_id : pieces_by_edge[overlap.edge_b])
					if (range_contains_piece(source.parameter_b_begin, source.parameter_b_end,
						pieces[piece_b_id], b.parameter_tolerance)
						&& !paired_b.contains(piece_b_id))
					{
						error = "reciprocal contact atom mapping was not symmetric";
						return false;
					}
			}

			std::map<std::size_t, std::vector<std::size_t>> groups_by_root;
			for (std::size_t piece = 0; piece < pieces.size(); ++piece)
				groups_by_root[physical.find(piece)].push_back(piece);
			std::vector<std::vector<std::size_t>> groups;
			for (auto& [unused, members] : groups_by_root)
			{
				(void)unused;
				std::sort(members.begin(), members.end(), [&](std::size_t a, std::size_t b)
				{
					return std::tie(edges[pieces[a].edge].id, pieces[a].begin, pieces[a].end)
						< std::tie(edges[pieces[b].edge].id, pieces[b].begin, pieces[b].end);
				});
				for (std::size_t i = 1; i < members.size(); ++i)
					if (pieces[members[i - 1]].edge == pieces[members[i]].edge)
					{
						error = "transitive reciprocal intervals mapped two atoms of one source edge together";
						return false;
					}
				groups.push_back(std::move(members));
			}
			std::sort(groups.begin(), groups.end(), [&](const auto& a, const auto& b)
			{
				const Piece& pa = pieces[a.front()];
				const Piece& pb = pieces[b.front()];
				return std::tie(edges[pa.edge].id, pa.begin, pa.end)
					< std::tie(edges[pb.edge].id, pb.begin, pb.end);
			});

			std::vector<AtomBuild> builds;
			builds.reserve(groups.size());
			for (const auto& members : groups)
			{
				const Piece& canonical_piece = pieces[members.front()];
				const EdgeData& canonical_edge = edges[canonical_piece.edge];
				const gp_Pnt canonical_low = edge_value(canonical_edge, canonical_piece.begin);
				const gp_Pnt canonical_high = edge_value(canonical_edge, canonical_piece.end);
				const bool canonical_increasing = !lexicographically_less(canonical_high,
					canonical_low);
				const gp_Pnt physical_begin = canonical_increasing ? canonical_low : canonical_high;
				const gp_Pnt physical_end = canonical_increasing ? canonical_high : canonical_low;

				AtomBuild build;
				std::map<std::size_t, OcctContactAtomSourceSpan> span_by_edge;
				for (std::size_t piece_id : members)
				{
					const Piece& piece = pieces[piece_id];
					const EdgeData& edge = edges[piece.edge];
					const gp_Pnt low = edge_value(edge, piece.begin);
					const gp_Pnt high = edge_value(edge, piece.end);
					const double forward_distance = low.Distance(physical_begin)
						+ high.Distance(physical_end);
					const double reverse_distance = high.Distance(physical_begin)
						+ low.Distance(physical_end);
					const double tolerance = native_pair_tolerance(canonical_edge, edge, 0.0);
					if (std::min(forward_distance, reverse_distance) > 2.0 * tolerance)
					{
						error = "reciprocal physical atom endpoints do not share one canonical direction";
						return false;
					}
					const bool increasing = forward_distance <= reverse_distance;
					span_by_edge[piece.edge] = {edge.id,
						increasing ? piece.begin : piece.end,
						increasing ? piece.end : piece.begin};
				}
				for (const auto& [unused, span] : span_by_edge)
				{
					(void)unused;
					build.atom.source_spans.push_back(span);
				}
				std::sort(build.atom.source_spans.begin(), build.atom.source_spans.end(),
					[](const auto& a, const auto& b)
					{
						return std::tie(a.source_edge_id, a.parameter_begin, a.parameter_end)
							< std::tie(b.source_edge_id, b.parameter_begin, b.parameter_end);
					});

				struct FaceAccumulator
				{
					std::uint8_t sectors = 0;
					OcctFaceIntervalLocation location = OcctFaceIntervalLocation::boundary;
					std::vector<OcctContactAtomBoundaryUse> occurrences;
					double exact_operation_tolerance = 0.0;
				};
				std::map<std::uint32_t, FaceAccumulator> face_uses;
				auto append_face_use = [&](std::size_t edge_index_value, double atom_begin,
					double atom_end, std::uint32_t face_id, OcctFaceIntervalLocation location,
					std::uint8_t sectors,
					const std::vector<OcctContactBoundaryOccurrence>& occurrences,
					double exact_operation_tolerance) -> bool
				{
					if (!finite(exact_operation_tolerance)
						|| exact_operation_tolerance < 0.0)
					{
						error = "face use has an invalid exact-operation tolerance";
						return false;
					}
					FaceAccumulator& accumulated = face_uses[face_id];
					if (sectors > accumulated.sectors
						|| (sectors == accumulated.sectors
							&& location == OcctFaceIntervalLocation::interior))
						accumulated.location = location;
					accumulated.sectors = std::max(accumulated.sectors, sectors);
					accumulated.exact_operation_tolerance = std::max(
						accumulated.exact_operation_tolerance,
						exact_operation_tolerance);
					for (const auto& occurrence : occurrences)
					{
						const auto key = std::make_pair(face_id, occurrence.occurrence_id);
						auto [stored, inserted] = build.boundary_edges.emplace(key,
							occurrence.target_edge);
						if (!inserted && !stored->second.IsEqual(occurrence.target_edge))
						{
							error = "one face occurrence ID resolved to different oriented target edges";
							return false;
						}
						if (!append_boundary_slice(edges[edge_index_value], occurrence,
							atom_begin, atom_end,
							accumulated.occurrences, error)) return false;
					}
					return true;
				};

				for (std::size_t piece_id : members)
				{
					const Piece& piece = pieces[piece_id];
					const EdgeData& edge = edges[piece.edge];
					const auto span = span_by_edge[piece.edge];
					for (const auto& use : edge.source->owner_face_uses)
						if (!append_face_use(piece.edge, span.parameter_begin, span.parameter_end,
							use.source_face_id, use.location, use.sector_count,
							use.boundary_occurrences, 0.0)) return false;
					const double midpoint = 0.5 * (piece.begin + piece.end);
					for (const ExactIntervalData& interval : intervals)
						if (interval.edge == piece.edge
							&& range_contains(interval.begin, interval.end, midpoint,
								edge.parameter_tolerance))
						{
							const auto& source = *interval.source;
							if (!append_face_use(piece.edge, span.parameter_begin,
								span.parameter_end, source.target_face_id,
								source.interval.location, source.interval.target_sector_count,
								exact_occurrences(interval),
								source.interval.exact_operation_tolerance)) return false;
						}
				}
				for (auto& [face, accumulated] : face_uses)
				{
					std::sort(accumulated.occurrences.begin(), accumulated.occurrences.end(),
						boundary_use_less);
					std::vector<OcctContactAtomBoundaryUse> unique_occurrences;
					for (const auto& occurrence : accumulated.occurrences)
					{
						if (!unique_occurrences.empty()
							&& unique_occurrences.back().occurrence_id == occurrence.occurrence_id)
						{
							const auto& previous = unique_occurrences.back();
							const auto previous_mapping = edge_index(edge_by_id,
								previous.mapping_source_edge_id);
							const auto current_mapping = edge_index(edge_by_id,
								occurrence.mapping_source_edge_id);
							const auto target_index = edge_index(edge_by_id,
								previous.source_edge_id);
							const auto target_occurrence = build.boundary_edges.find(
								{face, occurrence.occurrence_id});
							bool consistent = previous.source_edge_id == occurrence.source_edge_id
								&& previous.orientation == occurrence.orientation
								&& previous_mapping && current_mapping && target_index
								&& target_occurrence != build.boundary_edges.end()
								&& edges[*target_index].edge.IsSame(target_occurrence->second);
							if (consistent)
							{
								const EdgeData& previous_edge = edges[*previous_mapping];
								const EdgeData& current_edge = edges[*current_mapping];
								const EdgeData& target_edge = edges[*target_index];
								consistent = inside_parameter(previous_edge,
										previous.source_parameter_begin)
									&& inside_parameter(previous_edge,
										previous.source_parameter_end)
									&& inside_parameter(current_edge,
										occurrence.source_parameter_begin)
									&& inside_parameter(current_edge,
										occurrence.source_parameter_end)
									&& inside_parameter(target_edge,
										previous.target_parameter_begin)
									&& inside_parameter(target_edge,
										previous.target_parameter_end)
									&& inside_parameter(target_edge,
										occurrence.target_parameter_begin)
									&& inside_parameter(target_edge,
										occurrence.target_parameter_end);
								if (consistent)
								{
									const gp_Pnt previous_source_begin = edge_value(previous_edge,
										previous.source_parameter_begin);
									const gp_Pnt previous_source_end = edge_value(previous_edge,
										previous.source_parameter_end);
									const gp_Pnt current_source_begin = edge_value(current_edge,
										occurrence.source_parameter_begin);
									const gp_Pnt current_source_end = edge_value(current_edge,
										occurrence.source_parameter_end);
									const gp_Pnt previous_target_begin = edge_value(target_edge,
										previous.target_parameter_begin);
									const gp_Pnt previous_target_end = edge_value(target_edge,
										previous.target_parameter_end);
									const gp_Pnt current_target_begin = edge_value(target_edge,
										occurrence.target_parameter_begin);
									const gp_Pnt current_target_end = edge_value(target_edge,
										occurrence.target_parameter_end);
									const double source_bound = native_pair_tolerance(
										previous_edge, current_edge, 0.0);
									const double previous_target_bound = native_pair_tolerance(
										previous_edge, target_edge, 0.0);
									const double current_target_bound = native_pair_tolerance(
										current_edge, target_edge, 0.0);
									const double target_bound = native_pair_tolerance(
										target_edge, target_edge, 0.0);
									const double previous_target_delta =
										previous.target_parameter_end
										- previous.target_parameter_begin;
									const double current_target_delta =
										occurrence.target_parameter_end
										- occurrence.target_parameter_begin;
									consistent = std::signbit(previous_target_delta)
										== std::signbit(current_target_delta)
										&& std::abs(previous.target_parameter_begin
											- occurrence.target_parameter_begin)
											<= target_edge.parameter_tolerance
										&& std::abs(previous.target_parameter_end
											- occurrence.target_parameter_end)
											<= target_edge.parameter_tolerance
										&& previous_source_begin.SquareDistance(current_source_begin)
											<= source_bound * source_bound
										&& previous_source_end.SquareDistance(current_source_end)
											<= source_bound * source_bound
										&& previous_source_begin.SquareDistance(previous_target_begin)
											<= previous_target_bound * previous_target_bound
										&& previous_source_end.SquareDistance(previous_target_end)
											<= previous_target_bound * previous_target_bound
										&& current_source_begin.SquareDistance(current_target_begin)
											<= current_target_bound * current_target_bound
										&& current_source_end.SquareDistance(current_target_end)
											<= current_target_bound * current_target_bound
										&& previous_target_begin.SquareDistance(current_target_begin)
											<= target_bound * target_bound
										&& previous_target_end.SquareDistance(current_target_end)
											<= target_bound * target_bound;
								}
							}
							if (!consistent)
							{
								std::ostringstream report;
								report.precision(17);
								report << "duplicate face occurrence produced inconsistent atom parameter branches"
									<< ": face=" << face << ", occurrence="
									<< occurrence.occurrence_id;
								auto append_use = [&](const char* label,
									const OcctContactAtomBoundaryUse& use)
								{
									report << "; " << label << "={source_edge="
										<< use.source_edge_id << ", mapping_source_edge="
										<< use.mapping_source_edge_id << ", orientation="
										<< static_cast<int>(use.orientation)
										<< ", source=[" << use.source_parameter_begin
										<< ',' << use.source_parameter_end << "]"
										<< ", target=[" << use.target_parameter_begin
										<< ',' << use.target_parameter_end << "]";
									const auto source_index = edge_index(edge_by_id,
										use.mapping_source_edge_id);
									if (source_index)
									{
										const gp_Pnt begin_point = edge_value(edges[*source_index],
											use.source_parameter_begin);
										const gp_Pnt end_point = edge_value(edges[*source_index],
											use.source_parameter_end);
										report << ", source_xyz=[(" << begin_point.X() << ','
											<< begin_point.Y() << ',' << begin_point.Z()
											<< "),(" << end_point.X() << ',' << end_point.Y()
											<< ',' << end_point.Z() << ")]";
									}
									report << '}';
								};
								append_use("previous", previous);
								append_use("current", occurrence);
								error = report.str();
								return false;
							}
							// Sorted source_edge_id makes the retained provenance deterministic.
							continue;
						}
						unique_occurrences.push_back(occurrence);
					}
					accumulated.occurrences = std::move(unique_occurrences);
					OcctContactAtomFaceUse face_use;
					face_use.source_face_id = face;
					face_use.location = accumulated.location;
					face_use.sector_count = accumulated.sectors;
					face_use.boundary_occurrences = std::move(accumulated.occurrences);
					face_use.exact_operation_tolerance =
						accumulated.exact_operation_tolerance;
					build.atom.face_uses.push_back(std::move(face_use));
					build.atom.fan_degree += accumulated.sectors;
				}
				if (build.atom.fan_degree == 0)
				{
					error = "physical contact atom has no incident face sector";
					return false;
				}

				auto append_sample = [&](std::size_t source_edge_index, double source_parameter,
					double source_tolerance, std::uint8_t canonical_priority,
					std::optional<std::uint64_t> junction_id) -> bool
				{
					const EdgeData& source_edge = edges[source_edge_index];
					const gp_Pnt source_point = edge_value(source_edge, source_parameter);
					const double tolerance = source_edge.geometric_tolerance
						+ canonical_edge.geometric_tolerance + source_tolerance
						+ Precision::Confusion();
					double canonical_parameter = 0.0;
					gp_Pnt canonical_point;
					if (source_edge_index == canonical_piece.edge)
					{
						canonical_parameter = source_parameter;
						canonical_point = source_point;
					}
					else if (source_point.SquareDistance(canonical_low) <= tolerance * tolerance)
					{
						canonical_parameter = canonical_piece.begin;
						canonical_point = canonical_low;
					}
					else if (source_point.SquareDistance(canonical_high) <= tolerance * tolerance)
					{
						canonical_parameter = canonical_piece.end;
						canonical_point = canonical_high;
					}
					else if (!project_parameter(canonical_edge, source_point,
						canonical_piece.begin, canonical_piece.end, tolerance,
						canonical_parameter, canonical_point, error))
					{
						error = "could not recover a reciprocal source sample on its canonical atom: "
							+ error;
						return false;
					}
					if (source_point.SquareDistance(canonical_point) > tolerance * tolerance)
					{
						error = "reciprocal source sample is outside the canonical atom's native CAD bound";
						return false;
					}
					AtomSampleBuild sample;
					sample.canonical_parameter = canonical_parameter;
					sample.position = canonical_point;
					sample.tolerance = tolerance;
					sample.certificate_tolerance = std::max(source_edge.geometric_tolerance,
						source_tolerance);
					sample.canonical_priority = canonical_priority;
					sample.canonical_origin = {source_edge.id,
						clamped_parameter(source_edge, source_parameter)};
					sample.origins.push_back(sample.canonical_origin);
					if (junction_id) sample.junction_ids.push_back(*junction_id);
					build.sample_metadata.push_back(std::move(sample));
					return true;
				};

				for (std::size_t piece_id : members)
				{
					const Piece& piece = pieces[piece_id];
					const EdgeData& edge = edges[piece.edge];
					if (!append_sample(piece.edge, piece.begin, 0.0, 2, std::nullopt)
						|| !append_sample(piece.edge, piece.end, 0.0, 2, std::nullopt)) return false;
					for (const auto& sample : edge.source->samples)
						if (range_contains(piece.begin, piece.end, sample.parameter,
							edge.parameter_tolerance)
							&& !append_sample(piece.edge, sample.parameter, sample.tolerance,
								0, std::nullopt)) return false;
					for (const JunctionData& junction : junctions)
						for (const auto& incidence : junction.incidences)
							if (incidence.first == piece.edge
								&& range_contains(piece.begin, piece.end, incidence.second,
									edge.parameter_tolerance)
								&& !append_sample(piece.edge, incidence.second,
									junction.source->geometric_tolerance, 3,
									junction.source->junction_id)) return false;
				}

				std::sort(build.sample_metadata.begin(), build.sample_metadata.end(),
					[](const auto& a, const auto& b)
					{
						if (a.canonical_parameter != b.canonical_parameter)
							return a.canonical_parameter < b.canonical_parameter;
						return source_location_less(a.origins.front(), b.origins.front());
					});
				std::vector<AtomSampleBuild> unique_samples;
				for (AtomSampleBuild sample : build.sample_metadata)
				{
					if (!unique_samples.empty()
						&& std::abs(unique_samples.back().canonical_parameter
							- sample.canonical_parameter) <= canonical_edge.parameter_tolerance
						&& unique_samples.back().position.SquareDistance(sample.position)
							<= std::max(unique_samples.back().tolerance, sample.tolerance)
							* std::max(unique_samples.back().tolerance, sample.tolerance))
					{
						auto& previous = unique_samples.back();
						previous.tolerance = std::max(previous.tolerance, sample.tolerance);
						previous.certificate_tolerance = std::max(
							previous.certificate_tolerance, sample.certificate_tolerance);
						if (sample.canonical_priority > previous.canonical_priority
							|| (sample.canonical_priority == previous.canonical_priority
								&& source_location_less(sample.canonical_origin,
									previous.canonical_origin)))
						{
							previous.canonical_priority = sample.canonical_priority;
							previous.canonical_origin = sample.canonical_origin;
							previous.canonical_parameter = sample.canonical_parameter;
							previous.position = sample.position;
						}
						previous.origins.insert(previous.origins.end(), sample.origins.begin(),
							sample.origins.end());
						previous.junction_ids.insert(previous.junction_ids.end(),
							sample.junction_ids.begin(), sample.junction_ids.end());
						continue;
					}
					unique_samples.push_back(std::move(sample));
				}
				for (auto& sample : unique_samples)
				{
					// The common chain is consumed by every reciprocal face use.  Recover
					// an exact bounded parameter on every member edge for every union
					// sample, including samples inherited from only one tessellation.
					if (!sample.junction_ids.empty())
					{
						const std::uint64_t selected = *std::min_element(sample.junction_ids.begin(),
							sample.junction_ids.end());
						const auto found = std::find_if(junctions.begin(), junctions.end(),
							[&](const JunctionData& junction)
								{ return junction.source->junction_id == selected; });
						if (found == junctions.end())
						{
							error = "contact sample references an unknown exact junction";
							return false;
						}
						sample.position = as_point(found->source->canonical_world_position);
						sample.tolerance = std::max(sample.tolerance,
							found->source->geometric_tolerance + Precision::Confusion());
						sample.certificate_tolerance = std::max(sample.certificate_tolerance,
							found->source->geometric_tolerance);
					}
					const double canonical_tolerance = sample.tolerance;
					double aggregate_tolerance = canonical_tolerance;
					for (std::size_t piece_id : members)
					{
						const Piece& piece = pieces[piece_id];
						const EdgeData& edge = edges[piece.edge];
						const double tolerance = canonical_tolerance
							+ edge.geometric_tolerance + Precision::Confusion();
						double parameter = 0.0;
						gp_Pnt recovered;
						const gp_Pnt low = edge_value(edge, piece.begin);
						const gp_Pnt high = edge_value(edge, piece.end);
						if (sample.position.SquareDistance(low) <= tolerance * tolerance)
						{
							parameter = piece.begin;
							recovered = low;
						}
						else if (sample.position.SquareDistance(high) <= tolerance * tolerance)
						{
							parameter = piece.end;
							recovered = high;
						}
						else if (!project_parameter(edge, sample.position, piece.begin,
							piece.end, tolerance, parameter, recovered, error))
						{
							error = "could not recover a union-chain sample on every reciprocal member: "
								+ error;
							return false;
						}
						if (sample.position.SquareDistance(recovered) > tolerance * tolerance)
						{
							error = "union-chain sample is outside a reciprocal member's native CAD bound";
							return false;
						}
						sample.origins.push_back({edge.id, clamped_parameter(edge, parameter)});
						aggregate_tolerance = std::max(aggregate_tolerance, tolerance);
						sample.certificate_tolerance = std::max(sample.certificate_tolerance,
							tolerance);
					}
					sample.tolerance = aggregate_tolerance;
					std::sort(sample.origins.begin(), sample.origins.end(), source_location_less);
					sample.origins.erase(std::unique(sample.origins.begin(), sample.origins.end(),
						source_location_equal), sample.origins.end());
					std::sort(sample.junction_ids.begin(), sample.junction_ids.end());
					sample.junction_ids.erase(std::unique(sample.junction_ids.begin(),
						sample.junction_ids.end()), sample.junction_ids.end());
				}
				if (!canonical_increasing) std::reverse(unique_samples.begin(), unique_samples.end());
				if (unique_samples.size() < 2)
				{
					error = "physical contact atom has fewer than two distinct exact samples";
					return false;
				}
				build.sample_metadata = std::move(unique_samples);
				for (const auto& sample : build.sample_metadata)
					build.atom.samples.push_back({as_array(sample.position),
						sample.certificate_tolerance, no_id,
						sample.origins});
				builds.push_back(std::move(build));
			}

			// Equal source locations and explicit exact junction IDs are the only ways
			// topology nodes are shared.  Coincident but unrelated coordinates stay apart.
			std::vector<std::pair<std::size_t, std::size_t>> flat_location;
			for (std::size_t atom = 0; atom < builds.size(); ++atom)
				for (std::size_t sample = 0; sample < builds[atom].atom.samples.size(); ++sample)
					flat_location.emplace_back(atom, sample);
			DisjointSet topology(flat_location.size());
			struct OriginOccurrence
			{
				std::uint32_t edge = no_id;
				double parameter = 0.0;
				std::size_t flat = 0;
			};
			std::vector<OriginOccurrence> origins;
			std::map<std::uint64_t, std::vector<std::size_t>> flats_by_junction;
			for (std::size_t flat = 0; flat < flat_location.size(); ++flat)
			{
				const auto [atom, sample] = flat_location[flat];
				for (const auto& origin : builds[atom].sample_metadata[sample].origins)
					origins.push_back({origin.source_edge_id, origin.parameter, flat});
				for (std::uint64_t junction : builds[atom].sample_metadata[sample].junction_ids)
					flats_by_junction[junction].push_back(flat);
			}
			std::sort(origins.begin(), origins.end(), [](const auto& a, const auto& b)
			{
				return std::tie(a.edge, a.parameter, a.flat)
					< std::tie(b.edge, b.parameter, b.flat);
			});
			for (std::size_t begin = 0; begin < origins.size();)
			{
				std::size_t end = begin + 1;
				const auto index = edge_index(edge_by_id, origins[begin].edge);
				while (end < origins.size() && origins[end].edge == origins[begin].edge
					&& index && std::abs(origins[end].parameter - origins[begin].parameter)
						<= edges[*index].parameter_tolerance) ++end;
				for (std::size_t item = begin + 1; item < end; ++item)
					topology.merge(origins[begin].flat, origins[item].flat);
				begin = end;
			}
			for (const auto& [unused, flats] : flats_by_junction)
			{
				(void)unused;
				for (std::size_t item = 1; item < flats.size(); ++item)
					topology.merge(flats.front(), flats[item]);
			}

			std::map<std::uint64_t, const OcctContactExactJunction*> junction_by_id;
			for (const JunctionData& junction : junctions)
				junction_by_id.emplace(junction.source->junction_id, junction.source);
			struct TopologyGroup
			{
				std::size_t root = 0;
				std::vector<std::size_t> flats;
				std::pair<std::uint32_t, double> key{no_id, 0.0};
				std::array<double, 3> position{};
				double canonical_position_tolerance = 0.0;
			};
			std::map<std::size_t, TopologyGroup> topology_groups;
			for (std::size_t flat = 0; flat < flat_location.size(); ++flat)
			{
				const std::size_t root = topology.find(flat);
				auto& group = topology_groups[root];
				group.root = root;
				group.flats.push_back(flat);
			}
			std::vector<TopologyGroup*> ordered_topology;
			for (auto& [unused, group] : topology_groups)
			{
				(void)unused;
				bool has_key = false;
				bool has_preferred_origin = false;
				std::uint8_t preferred_priority = 0;
				OcctContactAtomSampleSource preferred_origin;
				const OcctContactExactJunction* selected_junction = nullptr;
				for (std::size_t flat : group.flats)
				{
					const auto [atom, sample] = flat_location[flat];
					const auto& metadata = builds[atom].sample_metadata[sample];
					if (!has_preferred_origin
						|| metadata.canonical_priority > preferred_priority
						|| (metadata.canonical_priority == preferred_priority
							&& source_location_less(metadata.canonical_origin,
								preferred_origin)))
					{
						has_preferred_origin = true;
						preferred_priority = metadata.canonical_priority;
						preferred_origin = metadata.canonical_origin;
					}
					for (const auto& origin : metadata.origins)
						if (!has_key || std::tie(origin.source_edge_id, origin.parameter)
							< std::tie(group.key.first, group.key.second))
						{
							group.key = {origin.source_edge_id, origin.parameter};
							has_key = true;
						}
					for (std::uint64_t junction_id : metadata.junction_ids)
					{
						const auto found = junction_by_id.find(junction_id);
						if (found != junction_by_id.end()
							&& (!selected_junction
								|| junction_id < selected_junction->junction_id))
							selected_junction = found->second;
					}
				}
				if (!has_key)
				{
					error = "contact topology node has no exact source parameter identity";
					return false;
				}
				if (selected_junction)
					group.position = selected_junction->canonical_world_position;
				else if (has_preferred_origin)
				{
					const auto index = edge_index(edge_by_id, preferred_origin.source_edge_id);
					group.position = as_array(edge_value(edges[*index], preferred_origin.parameter));
				}
				else
				{
					error = "contact topology node has no canonical exact source location";
					return false;
				}
				for (std::size_t flat : group.flats)
				{
					const auto [atom, sample] = flat_location[flat];
					const auto& metadata = builds[atom].sample_metadata[sample];
					const double tolerance = metadata.tolerance
						+ (selected_junction ? selected_junction->geometric_tolerance : 0.0)
						+ Precision::Confusion();
					if (metadata.position.SquareDistance(as_point(group.position))
						> tolerance * tolerance)
					{
						error = "joined topology samples disagree beyond their native CAD bounds";
						return false;
					}
					// The topology node publishes one bit-identical coordinate to every
					// participating atom. Retain the largest per-member bound used to
					// validate that replacement; publishing only the individual source
					// tolerance would understate a displaced exact-junction/canonical
					// selection.
					group.canonical_position_tolerance = std::max(
						group.canonical_position_tolerance, tolerance);
				}
				ordered_topology.push_back(&group);
			}
			std::sort(ordered_topology.begin(), ordered_topology.end(), [](const auto* a,
				const auto* b) { return a->key < b->key; });
			if (ordered_topology.size() >= no_id)
			{
				error = "contact atomizer exhausted the 32-bit topology ID space";
				return false;
			}
			for (std::size_t id = 0; id < ordered_topology.size(); ++id)
				for (std::size_t flat : ordered_topology[id]->flats)
				{
					const auto [atom, sample] = flat_location[flat];
					builds[atom].atom.samples[sample].topology_id =
						static_cast<std::uint32_t>(id);
					builds[atom].atom.samples[sample].canonical_world_position =
						ordered_topology[id]->position;
					builds[atom].atom.samples[sample].canonical_position_tolerance = std::max(
						builds[atom].atom.samples[sample].canonical_position_tolerance,
						ordered_topology[id]->canonical_position_tolerance);
				}

			// Publish an exact target-occurrence parameter for every common-chain
			// sample. Endpoint pairs only anchor the two atom endpoints; every interior
			// parameter is recovered uniquely on the bounded oriented occurrence.
			for (AtomBuild& build : builds)
				for (OcctContactAtomFaceUse& face_use : build.atom.face_uses)
					for (OcctContactAtomBoundaryUse& occurrence : face_use.boundary_occurrences)
					{
						const auto stored = build.boundary_edges.find({face_use.source_face_id,
							occurrence.occurrence_id});
						if (stored == build.boundary_edges.end() || stored->second.IsNull())
						{
							error = "atom boundary use lost its exact oriented target occurrence";
							return false;
						}
						EdgeData target;
						target.edge = stored->second;
						Standard_Real target_first = 0.0, target_last = 0.0;
						target.curve = BRep_Tool::Curve(target.edge, target.location,
							target_first, target_last);
						if (target.curve.IsNull() || !finite(target_first) || !finite(target_last)
							|| !(target_last > target_first))
						{
							error = "atom boundary occurrence has no finite target curve";
							return false;
						}
						target.first = target_first;
						target.last = target_last;
						target.geometric_tolerance = std::max(0.0,
							BRep_Tool::Tolerance(target.edge))
							* std::abs(target.location.Transformation().ScaleFactor());
						BRepAdaptor_Curve target_adaptor(target.edge);
						if (!resolved_parameter_tolerance(target_adaptor, target_first,
							target_last, target.geometric_tolerance + Precision::Confusion(),
							target.parameter_tolerance, error))
						{
							error = "atom boundary occurrence " + error;
							return false;
						}
						const double atom_target_first = std::min(occurrence.target_parameter_begin,
							occurrence.target_parameter_end);
						const double atom_target_last = std::max(occurrence.target_parameter_begin,
							occurrence.target_parameter_end);
						if (!(atom_target_last > atom_target_first)
							|| atom_target_first < target.first - target.parameter_tolerance
							|| atom_target_last > target.last + target.parameter_tolerance
							|| target.parameter_tolerance >= 1.0e-3
								* (atom_target_last - atom_target_first))
						{
							error = "atom boundary occurrence target span is unresolved";
							return false;
						}
						occurrence.sample_target_parameters.clear();
						occurrence.sample_target_parameters.reserve(build.atom.samples.size());
						for (std::size_t sample_id = 0;
							sample_id < build.atom.samples.size(); ++sample_id)
						{
							const OcctContactAtomSample& sample = build.atom.samples[sample_id];
							const gp_Pnt world = as_point(sample.canonical_world_position);
							const double linear_tolerance = sample.canonical_position_tolerance
								+ target.geometric_tolerance + Precision::Confusion();
							double parameter = 0.0;
							gp_Pnt recovered;
							if (sample_id == 0)
							{
								parameter = occurrence.target_parameter_begin;
								recovered = edge_value(target, parameter);
							}
							else if (sample_id + 1 == build.atom.samples.size())
							{
								parameter = occurrence.target_parameter_end;
								recovered = edge_value(target, parameter);
							}
							else if (!project_parameter(target, world, atom_target_first,
								atom_target_last, linear_tolerance, parameter, recovered, error))
							{
								error = "cannot recover a union sample on its exact target occurrence: "
									+ error;
								return false;
							}
							if (world.SquareDistance(recovered)
								> linear_tolerance * linear_tolerance)
							{
								error = "target-occurrence sample disagrees with its canonical world point";
								return false;
							}
							occurrence.sample_target_parameters.push_back(parameter);
						}
					}

			for (std::size_t atom = 0; atom < builds.size(); ++atom)
			{
				builds[atom].atom.id = atom;
				result.atoms.push_back(std::move(builds[atom].atom));
			}
			result.topology_node_count = static_cast<std::uint32_t>(ordered_topology.size());
		}
		catch (const Standard_Failure& failure)
		{
			error = std::string("OpenCascade contact atomization exception: ")
				+ failure.GetMessageString();
			return false;
		}
		catch (const std::exception& exception)
		{
			error = std::string("contact atomization exception: ") + exception.what();
			return false;
		}
		catch (...)
		{
			error = "unknown contact atomization exception";
			return false;
		}

		output = std::move(result);
		return true;
	}
}
