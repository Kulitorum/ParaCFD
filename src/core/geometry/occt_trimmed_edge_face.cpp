// Exact trimmed edge/face coverage and face-sector classification.

#include "core/geometry/occt_trimmed_edge_face.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <IntTools_CommonPrt.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <IntTools_Range.hxx>
#include <NCollection_List.hxx>
#include <NCollection_Sequence.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		using Range = std::pair<double, double>;

		double parameter_tolerance(double first, double last)
		{
			return std::max(Precision::PConfusion(),
				128.0 * std::numeric_limits<double>::epsilon()
				* std::max({ 1.0, std::abs(first), std::abs(last), std::abs(last - first) }));
		}

		bool normalize_range(Range& range, double first, double last, double tolerance,
			std::string& error)
		{
			if (range.second < range.first) std::swap(range.first, range.second);
			if (!std::isfinite(range.first) || !std::isfinite(range.second))
			{
				error = "OpenCascade returned a non-finite edge parameter range";
				return false;
			}
			if (range.first < first - tolerance || range.second > last + tolerance)
			{
				std::ostringstream report;
				report.precision(17);
				report << "exact edge mapping [" << range.first << ',' << range.second
					<< "] escaped finite range [" << first << ',' << last
					<< "] beyond parameter tolerance " << tolerance;
				error = report.str();
				return false;
			}
			range.first = std::clamp(range.first, first, last);
			range.second = std::clamp(range.second, first, last);
			return true;
		}

		double edge_world_tolerance(const TopoDS_Edge& edge)
		{
			return std::max(0.0, BRep_Tool::Tolerance(edge))
				* std::abs(edge.Location().Transformation().ScaleFactor());
		}

		double mapping_parameter_tolerance(const BRepAdaptor_Curve& curve,
			double first, double last, double linear_tolerance)
		{
			const double arithmetic = parameter_tolerance(first, last);
			const double resolved = std::abs(curve.Resolution(linear_tolerance));
			return std::isfinite(resolved) ? std::max(arithmetic, resolved) : arithmetic;
		}

		std::vector<Range> merged_ranges(std::vector<Range> ranges, double tolerance)
		{
			std::sort(ranges.begin(), ranges.end());
			std::vector<Range> result;
			for (const Range& range : ranges)
			{
				if (!(range.second > range.first + tolerance)) continue;
				if (result.empty() || range.first > result.back().second + tolerance)
					result.push_back(range);
				else
					result.back().second = std::max(result.back().second, range.second);
			}
			return result;
		}

		bool ranges_have_positive_overlap(std::vector<Range> ranges, double tolerance)
		{
			for (Range& range : ranges)
				if (range.second < range.first) std::swap(range.first, range.second);
			std::sort(ranges.begin(), ranges.end());
			for (std::size_t i = 1; i < ranges.size(); ++i)
				if (ranges[i].first < ranges[i - 1].second - tolerance) return true;
			return false;
		}

		bool ranges_cover(const std::vector<Range>& ranges, double first, double last,
			double tolerance)
		{
			double covered = first;
			for (const Range& range : ranges)
			{
				if (range.second < first || range.first > last) continue;
				if (range.first > covered + tolerance) return false;
				covered = std::max(covered, std::min(last, range.second));
			}
			return covered >= last - tolerance;
		}

		struct EdgeMapping
		{
			std::vector<Range> source_ranges;
			bool candidate_fully_covered = false;
			double source_parameter_tolerance = 0.0;
			std::string error;
		};

		EdgeMapping exact_edge_mapping(const TopoDS_Edge& source,
			const TopoDS_Edge& candidate, bool require_candidate_full)
		{
			EdgeMapping result;
			try
			{
				BRepAdaptor_Curve source_curve(source), candidate_curve(candidate);
				const double source_first = source_curve.FirstParameter();
				const double source_last = source_curve.LastParameter();
				const double candidate_first = candidate_curve.FirstParameter();
				const double candidate_last = candidate_curve.LastParameter();
				if (!std::isfinite(source_first) || !std::isfinite(source_last)
					|| !(source_last > source_first) || !std::isfinite(candidate_first)
					|| !std::isfinite(candidate_last) || !(candidate_last > candidate_first))
				{
					result.error = "cannot map a null or non-finite edge parameter range";
					return result;
				}

				// Put the candidate first.  For the Boolean fragments this is the trimmed
				// subset, which also avoids OCCT's analytic line/line path extending the
				// second range past its trimmed endpoint when the longer carrier is first.
				// Shape identities below remain authoritative because IntTools may still
				// reorder unlike curve types internally.
				IntTools_EdgeEdge intersection(candidate, source);
				intersection.SetFuzzyValue(0.0);
				intersection.Perform();
				if (!intersection.IsDone())
				{
					result.error = "zero-fuzzy IntTools_EdgeEdge did not complete";
					return result;
				}

				std::vector<Range> candidate_ranges;
				for (NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
					intersection.CommonParts()); common.More(); common.Next())
				{
					if (common.Value().Type() != TopAbs_EDGE) continue;
					const bool source_is_first = common.Value().Edge1().IsSame(source)
						&& common.Value().Edge2().IsSame(candidate);
					const bool source_is_second = common.Value().Edge1().IsSame(candidate)
						&& common.Value().Edge2().IsSame(source);
					if (!source_is_first && !source_is_second)
					{
						result.error = "IntTools common part does not identify the supplied edges";
						return result;
					}
					double first_begin = 0.0, first_end = 0.0;
					common.Value().Range1(first_begin, first_end);
					if (source_is_first)
						result.source_ranges.emplace_back(first_begin, first_end);
					else
						candidate_ranges.emplace_back(first_begin, first_end);
					std::size_t range2_count = 0;
					for (NCollection_Sequence<IntTools_Range>::Iterator range(
						common.Value().Ranges2()); range.More(); range.Next())
					{
						if (source_is_first)
							candidate_ranges.emplace_back(range.Value().First(),
								range.Value().Last());
						else
							result.source_ranges.emplace_back(range.Value().First(),
								range.Value().Last());
						++range2_count;
					}
					if (range2_count != 1)
					{
						result.error = "exact edge overlap has a non-unique parameter mapping";
						return result;
					}
				}

				const double linear_tolerance = std::max(Precision::Confusion(),
					edge_world_tolerance(source) + edge_world_tolerance(candidate));
				const double source_tolerance = mapping_parameter_tolerance(source_curve,
					source_first, source_last, linear_tolerance);
				result.source_parameter_tolerance = source_tolerance;
				for (Range& range : result.source_ranges)
					if (!normalize_range(range, source_first, source_last, source_tolerance,
						result.error))
					{
						result.error = "source-side " + result.error;
						return result;
					}
				if (ranges_have_positive_overlap(result.source_ranges, source_tolerance))
				{
					result.error = "exact edge overlap maps one candidate to overlapping source intervals";
					return result;
				}
				result.source_ranges = merged_ranges(std::move(result.source_ranges),
					source_tolerance);

				const double candidate_tolerance = mapping_parameter_tolerance(candidate_curve,
					candidate_first, candidate_last, linear_tolerance);
				for (Range& range : candidate_ranges)
					if (!normalize_range(range, candidate_first, candidate_last,
						candidate_tolerance, result.error))
					{
						result.error = "candidate-side " + result.error;
						return result;
					}
				candidate_ranges = merged_ranges(std::move(candidate_ranges),
					candidate_tolerance);
				result.candidate_fully_covered = ranges_cover(candidate_ranges,
					candidate_first, candidate_last, candidate_tolerance);
				if (require_candidate_full && !result.candidate_fully_covered)
					result.error = "Boolean common edge did not map completely back to the source edge";
			}
			catch (const Standard_Failure& failure)
			{
				result.error = std::string("OpenCascade edge-mapping exception: ")
					+ failure.GetMessageString();
			}
			catch (const std::exception& exception)
			{
				result.error = std::string("edge-mapping exception: ") + exception.what();
			}
			return result;
		}

		std::vector<TopoDS_Edge> boundary_occurrences(const TopoDS_Face& face)
		{
			std::vector<TopoDS_Edge> result;
			for (TopExp_Explorer wires(face, TopAbs_WIRE); wires.More(); wires.Next())
				// Preserve the wire occurrence's orientation and location on every child.
				// Disabling cumulative traversal silently maps a located face's boundary
				// edges in the unlocated model frame and loses FORWARD/REVERSED seam uses.
				for (TopoDS_Iterator child(wires.Current(), true, true); child.More();
					child.Next())
					if (child.Value().ShapeType() == TopAbs_EDGE)
						result.push_back(TopoDS::Edge(child.Value()));
			return result;
		}

		bool contains_parameter(const std::vector<Range>& ranges, double parameter,
			double tolerance)
		{
			return std::any_of(ranges.begin(), ranges.end(), [&](const Range& range)
			{
				return parameter >= range.first - tolerance
					&& parameter <= range.second + tolerance;
			});
		}

		struct BoundaryCoverage
		{
			TopoDS_Edge occurrence;
			std::vector<Range> ranges;
		};
	}

	OcctTrimmedEdgeFaceCommonResult exact_trimmed_edge_face_common(
		const TopoDS_Edge& source_edge, const TopoDS_Face& target_face)
	{
		OcctTrimmedEdgeFaceCommonResult result;
		auto fail = [&](std::string error)
		{
			result.intervals.clear();
			result.coverage = OcctTrimmedEdgeCoverage::none;
			result.errors.push_back(std::move(error));
		};

		try
		{
			if (source_edge.IsNull() || target_face.IsNull())
			{
				fail("source edge and target face must both be non-null");
				return result;
			}
			if (!BRepCheck_Analyzer(source_edge, true).IsValid()
				|| !BRepCheck_Analyzer(target_face, true).IsValid())
			{
				fail("source edge or target face is not a valid BRep");
				return result;
			}
			BRepAdaptor_Curve source_curve(source_edge);
			const double first = source_curve.FirstParameter();
			const double last = source_curve.LastParameter();
			if (!std::isfinite(first) || !std::isfinite(last) || !(last > first))
			{
				fail("source edge does not have one finite positive parameter range");
				return result;
			}
			result.source_parameter_first = first;
			result.source_parameter_last = last;
			const double tolerance = parameter_tolerance(first, last);
			double span_tolerance = tolerance;

			NCollection_List<TopoDS_Shape> arguments, tools;
			arguments.Append(source_edge);
			tools.Append(target_face);
			BRepAlgoAPI_Common common;
			common.SetArguments(arguments);
			common.SetTools(tools);
			common.SetFuzzyValue(0.0);
			common.SetNonDestructive(true);
			common.SetRunParallel(false);
			common.SetUseOBB(true);
			common.SetToFillHistory(false);
			common.Build();
			if (!common.IsDone() || common.HasErrors() || common.HasWarnings())
			{
				std::ostringstream report;
				if (!common.IsDone()) report << "BRepAlgoAPI_Common did not complete; ";
				if (common.HasErrors()) common.DumpErrors(report);
				if (common.HasWarnings()) common.DumpWarnings(report);
				fail("exact trimmed edge/face Common failed closed: " + report.str());
				return result;
			}
			if (common.Shape().IsNull())
			{
				fail("exact trimmed edge/face Common returned a null result shape");
				return result;
			}

			std::vector<Range> covered_ranges;
			TopTools_MapOfShape seen_common_edges;
			for (TopExp_Explorer edges(common.Shape(), TopAbs_EDGE); edges.More(); edges.Next())
			{
				const TopoDS_Edge edge = TopoDS::Edge(edges.Current());
				if (!seen_common_edges.Add(edge) || BRep_Tool::Degenerated(edge)) continue;
				const EdgeMapping mapping = exact_edge_mapping(source_edge, edge, true);
				if (!mapping.error.empty())
				{
					fail("cannot map a Boolean common edge: " + mapping.error);
					return result;
				}
				if (mapping.source_ranges.size() != 1)
				{
					fail("a Boolean common edge maps to multiple disjoint source intervals");
					return result;
				}
				span_tolerance = std::max(span_tolerance,
					mapping.source_parameter_tolerance);
				covered_ranges.push_back(mapping.source_ranges.front());
			}
			if (covered_ranges.empty())
			{
				result.coverage = OcctTrimmedEdgeCoverage::none;
				result.normalized_parameter_tolerance = tolerance / (last - first);
				return result;
			}

			std::vector<BoundaryCoverage> boundaries;
			for (const TopoDS_Edge& occurrence : boundary_occurrences(target_face))
			{
				if (BRep_Tool::Degenerated(occurrence)) continue;
				const EdgeMapping mapping = exact_edge_mapping(source_edge, occurrence, false);
				if (!mapping.error.empty())
				{
					fail("cannot classify a target-face boundary occurrence: " + mapping.error);
					return result;
				}
				if (!mapping.source_ranges.empty())
				{
					span_tolerance = std::max(span_tolerance,
						mapping.source_parameter_tolerance);
					boundaries.push_back({ occurrence, mapping.source_ranges });
				}
			}
			result.normalized_parameter_tolerance = span_tolerance / (last - first);
			if (ranges_have_positive_overlap(covered_ranges, span_tolerance))
			{
				fail("distinct Boolean common edges overlap on the source parameter range");
				return result;
			}
			covered_ranges = merged_ranges(std::move(covered_ranges), span_tolerance);

			std::vector<double> breakpoints;
			for (const Range& range : covered_ranges)
			{
				breakpoints.push_back(range.first);
				breakpoints.push_back(range.second);
				for (const BoundaryCoverage& boundary : boundaries)
					for (const Range& boundary_range : boundary.ranges)
					{
						const double begin = std::max(range.first, boundary_range.first);
						const double end = std::min(range.second, boundary_range.second);
						if (end > begin + span_tolerance)
						{
							breakpoints.push_back(begin);
							breakpoints.push_back(end);
						}
					}
			}
			std::sort(breakpoints.begin(), breakpoints.end());
			breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end(),
				[&](double a, double b) { return std::abs(a - b) <= span_tolerance; }),
				breakpoints.end());

			const double inverse_range = 1.0 / (last - first);
			for (std::size_t point = 1; point < breakpoints.size(); ++point)
			{
				const double begin = breakpoints[point - 1], end = breakpoints[point];
				if (!(end > begin + span_tolerance)) continue;
				const double midpoint = 0.5 * (begin + end);
				if (!contains_parameter(covered_ranges, midpoint, span_tolerance)) continue;

				std::vector<const BoundaryCoverage*> occurrences;
				for (const BoundaryCoverage& boundary : boundaries)
					if (contains_parameter(boundary.ranges, midpoint, span_tolerance))
						occurrences.push_back(&boundary);

				OcctFaceIntervalLocation location = OcctFaceIntervalLocation::interior;
				std::uint8_t sectors = 2;
				if (occurrences.size() == 1)
				{
					if (!BRep_Tool::IsClosed(occurrences.front()->occurrence, target_face))
					{
						location = OcctFaceIntervalLocation::boundary;
						sectors = 1;
					}
				}
				else if (occurrences.size() == 2)
				{
					const TopoDS_Edge& a = occurrences[0]->occurrence;
					const TopoDS_Edge& b = occurrences[1]->occurrence;
					if (!a.IsSame(b) || !BRep_Tool::IsClosed(a, target_face)
						|| a.Orientation() == b.Orientation())
					{
						fail("two unrelated target-face boundary occurrences cover one source interval");
						return result;
					}
				}
				else if (occurrences.size() > 2)
				{
					fail("more than two target-face boundary occurrences cover one source interval");
					return result;
				}

				OcctTrimmedEdgeFaceInterval interval;
				interval.begin = std::clamp((begin - first) * inverse_range, 0.0, 1.0);
				interval.end = std::clamp((end - first) * inverse_range, 0.0, 1.0);
				interval.location = location;
				interval.target_sector_count = sectors;
				interval.boundary_occurrence_count = static_cast<std::uint8_t>(
					occurrences.size());
				// Boundary occurrences carry distinct face pcurves.  Even when two
				// collinear trim edges meet and have identical one-sector classification,
				// their common vertex is a required pcurve handoff and must survive in the
				// interval graph.  Interior atoms have no occurrence identity and may merge.
				if (interval.boundary_occurrence_count == 0
					&& !result.intervals.empty()
					&& result.intervals.back().boundary_occurrence_count == 0
					&& result.intervals.back().location == interval.location
					&& result.intervals.back().target_sector_count == interval.target_sector_count
					&& result.intervals.back().boundary_occurrence_count
						== interval.boundary_occurrence_count
					&& std::abs(result.intervals.back().end - interval.begin)
						<= result.normalized_parameter_tolerance)
					result.intervals.back().end = interval.end;
				else
					result.intervals.push_back(interval);
			}

			if (result.intervals.empty())
			{
				fail("positive-length Common coverage produced no classified source interval");
				return result;
			}
			result.coverage = ranges_cover(covered_ranges, first, last, span_tolerance)
				? OcctTrimmedEdgeCoverage::full : OcctTrimmedEdgeCoverage::partial;
		}
		catch (const Standard_Failure& failure)
		{
			fail(std::string("OpenCascade trimmed edge/face exception: ")
				+ failure.GetMessageString());
		}
		catch (const std::exception& exception)
		{
			fail(std::string("trimmed edge/face exception: ") + exception.what());
		}
		return result;
	}
}
