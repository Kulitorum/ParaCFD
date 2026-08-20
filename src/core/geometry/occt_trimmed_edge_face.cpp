// Exact trimmed edge/face coverage and face-sector classification.

#include "core/geometry/occt_trimmed_edge_face.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BOPAlgo_Alerts.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <GeomConvert.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IntTools_CommonPrt.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <IntTools_EdgeFace.hxx>
#include <IntTools_Range.hxx>
#include <Message_Alert.hxx>
#include <Message_Report.hxx>
#include <NCollection_List.hxx>
#include <NCollection_Sequence.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <TopLoc_Location.hxx>

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
		constexpr double kMaximumResolvedParameterFraction=1.0e-3;

		double parameter_tolerance(double first, double last)
		{
			const double span=std::abs(last-first);
			auto ulp=[](double value)
			{
				return std::abs(std::nextafter(value,
					std::numeric_limits<double>::infinity())-value);
			};
			return std::max({8.0*ulp(first),8.0*ulp(last),
				128.0*std::numeric_limits<double>::epsilon()*span});
		}

		bool parameter_tolerance_is_resolved(double tolerance,double first,double last,
			const char* label,std::string& error)
		{
			const double span=std::abs(last-first);
			if(std::isfinite(tolerance)&&tolerance>=0.0&&span>0.0
				&&tolerance<kMaximumResolvedParameterFraction*span)return true;
			std::ostringstream report;
			report.precision(17);
			report<<label<<" parameter tolerance "<<tolerance
				<<" is not resolved relative to finite span "<<span;
			error=report.str();
			return false;
		}

		TopoDS_Edge normalized_parameter_copy(const TopoDS_Edge& source,double first,
			double last,std::string& error)
		{
			TopLoc_Location location;Standard_Real curve_first=0.0,curve_last=0.0;
			const Handle(Geom_Curve) curve=BRep_Tool::Curve(source,location,curve_first,curve_last);
			if(curve.IsNull())
			{
				error="source edge has no 3D curve to reparameterize";
				return {};
			}
			Handle(Geom_BSplineCurve) bspline=
				Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
			if(!bspline.IsNull()
				&&(first>bspline->FirstParameter()||last<bspline->LastParameter()))
				bspline->Segment(first,last);
			if(bspline.IsNull())
			{
				const Handle(Geom_TrimmedCurve) trimmed=new Geom_TrimmedCurve(curve,first,last);
				bspline=GeomConvert::CurveToBSplineCurve(trimmed);
			}
			if(bspline.IsNull()||bspline->NbKnots()<2)
			{
				error="source edge could not be converted exactly to a finite B-spline";
				return {};
			}
			TColStd_Array1OfReal knots(1,bspline->NbKnots());bspline->Knots(knots);
			const double knot_first=knots(knots.Lower());
			const double knot_last=knots(knots.Upper());
			if(!std::isfinite(knot_first)||!std::isfinite(knot_last)||!(knot_last>knot_first))
			{
				error="converted source edge has no finite positive knot range";
				return {};
			}
			for(Standard_Integer knot=knots.Lower();knot<=knots.Upper();++knot)
				knots.SetValue(knot,(knots(knot)-knot_first)/(knot_last-knot_first));
			bspline->SetKnots(knots);
			BRepBuilderAPI_MakeEdge builder(bspline,0.0,1.0);
			if(!builder.IsDone())
			{
				error="could not build normalized finite source edge";
				return {};
			}
			TopoDS_Edge result=builder.Edge();
			if(!location.IsIdentity())result=TopoDS::Edge(result.Moved(location));
			BRep_Builder topology_builder;
			topology_builder.UpdateEdge(result,BRep_Tool::Tolerance(source));
			result.Orientation(source.Orientation());
			return result;
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

		bool unique_edge_parameter(const gp_Pnt& world_point,const TopoDS_Edge& edge,
			double first,double last,double parameter_tolerance_value,
			double linear_tolerance,double& parameter,std::string& error)
		{
			TopLoc_Location location;Standard_Real curve_first=0.0,curve_last=0.0;
			const Handle(Geom_Curve) curve=BRep_Tool::Curve(edge,location,curve_first,curve_last);
			if(curve.IsNull()||!std::isfinite(first)||!std::isfinite(last)||!(last>first))
			{
				error="bounded edge has no finite positive 3D curve range";
				return false;
			}
			gp_Pnt local_point=world_point;
			local_point.Transform(location.Transformation().Inverted());
			GeomAPI_ProjectPointOnCurve projection(local_point,curve,first,last);
			std::vector<double> parameters;
			for(Standard_Integer candidate=1;candidate<=projection.NbPoints();++candidate)
			{
				const double value=projection.Parameter(candidate);
				if(!std::isfinite(value)||value<first-parameter_tolerance_value
					||value>last+parameter_tolerance_value)continue;
				const double bounded=std::clamp(value,first,last);
				gp_Pnt point=curve->Value(bounded);point.Transform(location.Transformation());
				if(point.Distance(world_point)<=linear_tolerance)parameters.push_back(bounded);
			}
			std::sort(parameters.begin(),parameters.end());
			parameters.erase(std::unique(parameters.begin(),parameters.end(),
				[&](double a,double b){return std::abs(a-b)<=parameter_tolerance_value;}),
				parameters.end());
			if(parameters.size()!=1)
			{
				std::ostringstream report;report<<"point has "<<parameters.size()
					<<" native-bounded edge parameters (expected exactly one)";
				error=report.str();
				return false;
			}
			parameter=parameters.front();return true;
		}

		bool same_parameter_endpoint(double a,double b)
		{
			if(a==b)return true;
			if(!std::isfinite(a)||!std::isfinite(b))return false;
			// Endpoint values independently evaluated by OCCT may differ by a handful
			// of representable doubles.  This fixed ULP comparison is arithmetic identity,
			// not a CAD/world tolerance: it cannot grow with an unrelated edge or face.
			double cursor=a;
			for(unsigned step=0;step<8;++step)
			{
				cursor=std::nextafter(cursor,b);
				if(cursor==b)return true;
			}
			return false;
		}

		bool positive_parameter_span(const Range& range)
		{
			return range.second>range.first
				&&!same_parameter_endpoint(range.first,range.second);
		}

		std::vector<Range> exact_range_union(std::vector<Range> ranges)
		{
			std::sort(ranges.begin(), ranges.end());
			std::vector<Range> result;
			for (const Range& range : ranges)
			{
				if (!positive_parameter_span(range)) continue;
				if (result.empty() || (range.first > result.back().second
					&& !same_parameter_endpoint(range.first,result.back().second)))
					result.push_back(range);
				else
					result.back().second = std::max(result.back().second, range.second);
			}
			return result;
		}

		bool ranges_have_positive_overlap(std::vector<Range> ranges)
		{
			for (Range& range : ranges)
				if (range.second < range.first) std::swap(range.first, range.second);
			std::sort(ranges.begin(), ranges.end());
			for (std::size_t i = 1; i < ranges.size(); ++i)
				if (ranges[i].first < ranges[i - 1].second
					&&!same_parameter_endpoint(ranges[i].first,ranges[i-1].second)) return true;
			return false;
		}

		bool exact_ranges_cover(const std::vector<Range>& ranges,double first,double last)
		{
			double covered = first;
			for (const Range& range : ranges)
			{
				if (range.second < first || range.first > last) continue;
				if (range.first > covered
					&&!same_parameter_endpoint(range.first,covered)) return false;
				covered = std::max(covered, std::min(last, range.second));
			}
			return covered>=last||same_parameter_endpoint(covered,last);
		}

		bool has_only_unable_to_orient_warnings(const BRepAlgoAPI_Common& common)
		{
			if(!common.HasWarnings())return true;
			const Handle(Message_Report)& report=common.GetReport();
			if(report.IsNull())return false;
			const auto& warnings=report->GetAlerts(Message_Warning);
			if(warnings.IsEmpty())return false;
			for(NCollection_List<Handle(Message_Alert)>::Iterator warning(warnings);
				warning.More();warning.Next())
			{
				const Handle(Message_Alert)& alert=warning.Value();
				if(alert.IsNull()
					||!alert->IsKind(STANDARD_TYPE(BOPAlgo_AlertUnableToOrientTheShape)))
					return false;
			}
			return true;
		}

		bool edge_face_fallback_ranges(const TopoDS_Edge& source,const TopoDS_Face& face,
			double first,double last,double tolerance,std::vector<Range>& ranges,
			std::string& error)
		{
			IntTools_EdgeFace intersection;
			intersection.SetEdge(source);intersection.SetFace(face);
			intersection.SetRange(first,last);
			// IntTools_EdgeFace clamps this to Precision::Confusion and combines it with
			// the native edge/face tolerances. It is used only to audit an otherwise empty
			// Boolean result carrying the one explicitly accepted orientation warning.
			intersection.SetFuzzyValue(0.0);
			intersection.UseQuickCoincidenceCheck(false);
			intersection.Perform();
			if(!intersection.IsDone())
			{
				std::ostringstream report;report<<"IntTools_EdgeFace fallback failed with status "
					<<intersection.ErrorStatus();error=report.str();return false;
			}
			for(NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
				intersection.CommonParts());common.More();common.Next())
			{
				if(common.Value().Type()!=TopAbs_EDGE)continue;
				double begin=0.0,end=0.0;common.Value().Range1(begin,end);
				Range range{begin,end};
				if(!normalize_range(range,first,last,tolerance,error))return false;
				if(positive_parameter_span(range))ranges.push_back(range);
			}
			return true;
		}

		struct EdgeMapping
		{
			struct MappedRange
			{
				Range source;
				double candidate_begin=0.0;
				double candidate_end=0.0;
			};
			std::vector<Range> source_ranges;
			std::vector<MappedRange> mapped_ranges;
			bool candidate_fully_covered = false;
			double source_parameter_tolerance = 0.0;
			double candidate_parameter_tolerance = 0.0;
			double linear_tolerance = 0.0;
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

				const double linear_tolerance = std::max(Precision::Confusion(),
					edge_world_tolerance(source) + edge_world_tolerance(candidate));
				result.linear_tolerance=linear_tolerance;
				const double source_tolerance = mapping_parameter_tolerance(source_curve,
					source_first, source_last, linear_tolerance);
				if(!parameter_tolerance_is_resolved(source_tolerance,source_first,source_last,
					"source",result.error))return result;
				result.source_parameter_tolerance = source_tolerance;
				const double candidate_tolerance = mapping_parameter_tolerance(candidate_curve,
					candidate_first, candidate_last, linear_tolerance);
				if(!parameter_tolerance_is_resolved(candidate_tolerance,candidate_first,
					candidate_last,"candidate",result.error))return result;
				result.candidate_parameter_tolerance=candidate_tolerance;

				// IntTools can reorder its inputs, and its analytic line/line path can return
				// a carrier-curve range extending beyond one trimmed endpoint.  Range1/2
				// identity is therefore not accepted as the mapping.  Treat every reported
				// range as a source-range hypothesis, then recover both target endpoints by
				// unique projection onto the finite candidate edge at native tolerance.
				IntTools_EdgeEdge intersection(candidate, source);
				intersection.SetFuzzyValue(0.0);
				intersection.Perform();
				if (!intersection.IsDone())
				{
					result.error = "zero-fuzzy IntTools_EdgeEdge did not complete";
					return result;
				}

				std::vector<EdgeMapping::MappedRange> mapped_hypotheses;
				bool saw_edge_common=false;
				for (NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
					intersection.CommonParts()); common.More(); common.Next())
				{
					if (common.Value().Type() != TopAbs_EDGE) continue;
					saw_edge_common=true;
					const bool source_is_first=common.Value().Edge1().IsSame(source)
						&&common.Value().Edge2().IsSame(candidate);
					const bool source_is_second=common.Value().Edge1().IsSame(candidate)
						&&common.Value().Edge2().IsSame(source);
					double first_begin = 0.0, first_end = 0.0;
					common.Value().Range1(first_begin, first_end);
					std::vector<Range> second_ranges;
					for (NCollection_Sequence<IntTools_Range>::Iterator range(
						common.Value().Ranges2()); range.More(); range.Next())
						second_ranges.emplace_back(range.Value().First(),range.Value().Last());
					auto append_hypothesis=[&](Range source_range)
					{
						std::string ignored_error;
						if(!normalize_range(source_range,source_first,source_last,
							source_tolerance,ignored_error)
							||!positive_parameter_span(source_range))return false;
						EdgeMapping::MappedRange mapped;mapped.source=source_range;
						if(!unique_edge_parameter(source_curve.Value(source_range.first),candidate,
							candidate_first,candidate_last,candidate_tolerance,linear_tolerance,
							mapped.candidate_begin,ignored_error)
							||!unique_edge_parameter(source_curve.Value(source_range.second),candidate,
								candidate_first,candidate_last,candidate_tolerance,linear_tolerance,
								mapped.candidate_end,ignored_error))return false;
						mapped_hypotheses.push_back(mapped);
						return true;
					};
					auto append_reported_mapping=[&](Range source_range,Range candidate_range)
					{
						std::string ignored_error;
						if(!normalize_range(source_range,source_first,source_last,
							source_tolerance,ignored_error)
							||!normalize_range(candidate_range,candidate_first,candidate_last,
								candidate_tolerance,ignored_error)
							||!positive_parameter_span(source_range)
							||!positive_parameter_span(candidate_range))
							return false;
						const gp_Pnt source_begin=source_curve.Value(source_range.first);
						const gp_Pnt source_end=source_curve.Value(source_range.second);
						const gp_Pnt candidate_begin=candidate_curve.Value(candidate_range.first);
						const gp_Pnt candidate_end=candidate_curve.Value(candidate_range.second);
						EdgeMapping::MappedRange mapped;mapped.source=source_range;
						if(source_begin.SquareDistance(candidate_begin)
							+source_end.SquareDistance(candidate_end)
							<=source_begin.SquareDistance(candidate_end)
								+source_end.SquareDistance(candidate_begin))
						{
							mapped.candidate_begin=candidate_range.first;
							mapped.candidate_end=candidate_range.second;
						}
						else
						{
							mapped.candidate_begin=candidate_range.second;
							mapped.candidate_end=candidate_range.first;
						}
						if(source_begin.Distance(candidate_curve.Value(mapped.candidate_begin))
							>linear_tolerance
							||source_end.Distance(candidate_curve.Value(mapped.candidate_end))
								>linear_tolerance)return false;
						mapped_hypotheses.push_back(mapped);return true;
					};
					bool mapped=false;
					if(source_is_first&&second_ranges.size()==1)
						mapped=append_reported_mapping({first_begin,first_end},
							second_ranges.front());
					else if(source_is_second&&second_ranges.size()==1)
						mapped=append_reported_mapping(second_ranges.front(),
							{first_begin,first_end});
					if(!mapped)mapped=append_hypothesis({first_begin,first_end});
					if(!mapped)for(const Range& second_range:second_ranges)
						if((mapped=append_hypothesis(second_range)))break;
					// OCCT's analytic line/line common can associate Range1 with the
					// wrong supplied edge when the source is wholly inside a longer
					// candidate. The finite full-source fallback is accepted only when
					// both of its endpoints uniquely round-trip onto that candidate.
					if(!mapped)append_hypothesis({source_first,source_last});
				}
				std::sort(mapped_hypotheses.begin(),mapped_hypotheses.end(),
					[](const auto& a,const auto& b)
					{return std::tie(a.source,a.candidate_begin,a.candidate_end)
						<std::tie(b.source,b.candidate_begin,b.candidate_end);});
				mapped_hypotheses.erase(std::unique(mapped_hypotheses.begin(),
					mapped_hypotheses.end(),[](const auto& a,const auto& b)
					{return same_parameter_endpoint(a.source.first,b.source.first)
						&&same_parameter_endpoint(a.source.second,b.source.second)
						&&same_parameter_endpoint(a.candidate_begin,b.candidate_begin)
						&&same_parameter_endpoint(a.candidate_end,b.candidate_end);}),
					mapped_hypotheses.end());
				if(saw_edge_common&&mapped_hypotheses.empty())
				{
					result.error="exact edge overlap has no unique native-bounded parameter mapping";
					return result;
				}
				std::vector<Range> candidate_ranges;
				for(const auto& mapped:mapped_hypotheses)
				{
					result.source_ranges.push_back(mapped.source);
					candidate_ranges.emplace_back(std::min(mapped.candidate_begin,
						mapped.candidate_end),std::max(mapped.candidate_begin,mapped.candidate_end));
					result.mapped_ranges.push_back(mapped);
				}
				std::sort(result.mapped_ranges.begin(),result.mapped_ranges.end(),
					[](const auto& a,const auto& b){return a.source<b.source;});
				if (ranges_have_positive_overlap(result.source_ranges))
				{
					result.error = "exact edge overlap maps one candidate to overlapping source intervals";
					return result;
				}
				result.source_ranges=exact_range_union(std::move(result.source_ranges));
				candidate_ranges=exact_range_union(std::move(candidate_ranges));
				result.candidate_fully_covered=exact_ranges_cover(candidate_ranges,
					candidate_first,candidate_last);
				if (require_candidate_full && !result.candidate_fully_covered)
				{
					std::ostringstream report;report.precision(17);
					report<<"Boolean common edge did not map completely back to its finite range ["
						<<candidate_first<<','<<candidate_last<<"] from";
					for(const Range& range:candidate_ranges)
						report<<" ["<<range.first<<','<<range.second<<']';
					result.error=report.str();
				}
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

		bool contains_parameter(const std::vector<Range>& ranges,double parameter)
		{
			return std::any_of(ranges.begin(), ranges.end(), [&](const Range& range)
			{
				return (parameter>=range.first
					||same_parameter_endpoint(parameter,range.first))
					&&(parameter<=range.second
						||same_parameter_endpoint(parameter,range.second));
			});
		}

		struct BoundaryCoverage
		{
			TopoDS_Edge occurrence;
			std::uint32_t occurrence_id=0;
			std::vector<Range> ranges;
			std::vector<EdgeMapping::MappedRange> mapped_ranges;
			double source_parameter_tolerance=0.0;
			double target_parameter_tolerance=0.0;
			double linear_tolerance=0.0;
		};

		struct CommonCoverageFragment
		{
			Range source;
			// Native tolerance carried by this particular BRepAlgo Common result
			// edge, in the located model/world frame.
			double exact_operation_tolerance = 0.0;
		};

		bool unique_target_parameter_at_source(const BRepAdaptor_Curve& source_curve,
			double source_parameter,const TopoDS_Edge& target_occurrence,
			const EdgeMapping::MappedRange& mapping,double source_parameter_tolerance,
			double target_parameter_tolerance,
			double linear_tolerance,double& target_parameter,std::string& error)
		{
			if(std::abs(source_parameter-mapping.source.first)<=source_parameter_tolerance)
			{
				target_parameter=mapping.candidate_begin;return true;
			}
			if(std::abs(source_parameter-mapping.source.second)<=source_parameter_tolerance)
			{
				target_parameter=mapping.candidate_end;return true;
			}
			TopLoc_Location location;Standard_Real target_first=0.0,target_last=0.0;
			const Handle(Geom_Curve) target_curve=BRep_Tool::Curve(target_occurrence,
				location,target_first,target_last);
			if(target_curve.IsNull())
			{
				error="target boundary occurrence has no 3D curve";
				return false;
			}
			const double mapped_first=std::min(mapping.candidate_begin,mapping.candidate_end);
			const double mapped_last=std::max(mapping.candidate_begin,mapping.candidate_end);
			if(!std::isfinite(mapped_first)||!std::isfinite(mapped_last)
				||!(mapped_last>mapped_first))
			{
				error="target boundary occurrence has no finite positive mapped range";
				return false;
			}
			gp_Pnt source_point=source_curve.Value(source_parameter);
			gp_Pnt local_source=source_point;
			local_source.Transform(location.Transformation().Inverted());
			GeomAPI_ProjectPointOnCurve projection(local_source,target_curve,
				mapped_first,mapped_last);
			std::vector<double> parameters;
			for(Standard_Integer candidate=1;candidate<=projection.NbPoints();++candidate)
			{
				const double parameter=projection.Parameter(candidate);
				if(!std::isfinite(parameter)||parameter<mapped_first-target_parameter_tolerance
					||parameter>mapped_last+target_parameter_tolerance)continue;
				gp_Pnt point=target_curve->Value(std::clamp(parameter,mapped_first,mapped_last));
				point.Transform(location.Transformation());
				if(point.Distance(source_point)>linear_tolerance)continue;
				parameters.push_back(std::clamp(parameter,mapped_first,mapped_last));
			}
			std::sort(parameters.begin(),parameters.end());
			parameters.erase(std::unique(parameters.begin(),parameters.end(),
				[&](double a,double b){return std::abs(a-b)<=target_parameter_tolerance;}),
				parameters.end());
			if(parameters.size()!=1)
			{
				std::ostringstream report;
				report<<"emitted source endpoint has "<<parameters.size()
					<<" native-bounded target parameters (expected exactly one)";
				error=report.str();
				return false;
			}
			target_parameter=parameters.front();
			return true;
		}
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
			if (source_edge.IsNull())
			{
				fail("source edge must be non-null");
				return result;
			}
			if (!BRepCheck_Analyzer(source_edge, true).IsValid())
			{
				fail("source edge is not a valid BRep");
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
			std::string tolerance_error;
			if(!parameter_tolerance_is_resolved(tolerance,first,last,"source",
				tolerance_error))
			{
				fail(std::move(tolerance_error));
				return result;
			}
			TopoDS_Edge audit_source=source_edge;
			const double source_span=last-first;
			if(source_span<1.0e-8||tolerance>1.0e-10*source_span)
			{
				std::string normalization_error;
				audit_source=normalized_parameter_copy(source_edge,first,last,
					normalization_error);
				if(audit_source.IsNull())
				{
					fail("cannot normalize a poorly conditioned source parameterization: "
						+normalization_error);
					return result;
				}
			}
			BRepAdaptor_Curve audit_source_curve(audit_source);
			const double audit_first=audit_source_curve.FirstParameter();
			const double audit_last=audit_source_curve.LastParameter();
			const double audit_span=audit_last-audit_first;
			const double audit_mapping_tolerance=parameter_tolerance(audit_first,audit_last);
			if(!parameter_tolerance_is_resolved(audit_mapping_tolerance,audit_first,audit_last,
				"audit source",tolerance_error))
			{
				fail(std::move(tolerance_error));
				return result;
			}
			if(target_face.IsNull()||!BRepCheck_Analyzer(target_face,true).IsValid())
			{
				fail("target face is null or is not a valid BRep");
				return result;
			}

			NCollection_List<TopoDS_Shape> arguments, tools;
			arguments.Append(audit_source);
			tools.Append(target_face);
			BRepAlgoAPI_Common common;
			common.SetArguments(arguments);
			common.SetTools(tools);
			common.SetFuzzyValue(0.0);
			common.SetNonDestructive(true);
			common.SetRunParallel(false);
			// Do not let a pre-existing visualization triangulation participate in this
			// exact CAD certificate.  OCCT's OBB shortcut can be tessellation-backed and,
			// on a meshed STEP face with multiple disjoint common spans, was observed to
			// discard the short span while the identical unmeshed BRep retained it.
			common.SetUseOBB(false);
			common.SetToFillHistory(false);
			common.Build();
			const bool only_orientation_warnings=has_only_unable_to_orient_warnings(common);
			if (!common.IsDone() || common.HasErrors() || !only_orientation_warnings)
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
			std::vector<CommonCoverageFragment> common_fragments;
			TopTools_MapOfShape seen_common_edges;
			for (TopExp_Explorer edges(common.Shape(), TopAbs_EDGE); edges.More(); edges.Next())
			{
				const TopoDS_Edge edge = TopoDS::Edge(edges.Current());
				if (!seen_common_edges.Add(edge) || BRep_Tool::Degenerated(edge)) continue;
				const EdgeMapping mapping = exact_edge_mapping(audit_source, edge, true);
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
				const double result_edge_tolerance = BRep_Tool::Tolerance(edge)
					* std::abs(edge.Location().Transformation().ScaleFactor());
				if (!std::isfinite(result_edge_tolerance)
					|| result_edge_tolerance < 0.0)
				{
					fail("a Boolean common result edge has an invalid native tolerance");
					return result;
				}
				covered_ranges.push_back(mapping.source_ranges.front());
				common_fragments.push_back({mapping.source_ranges.front(),
					result_edge_tolerance});
			}
			std::vector<BoundaryCoverage> boundaries;
			const std::vector<TopoDS_Edge> target_occurrences=boundary_occurrences(target_face);
			for (std::size_t occurrence_id=0;occurrence_id<target_occurrences.size();++occurrence_id)
			{
				const TopoDS_Edge& occurrence=target_occurrences[occurrence_id];
				if (BRep_Tool::Degenerated(occurrence)) continue;
				const EdgeMapping mapping = exact_edge_mapping(audit_source, occurrence, false);
				if (!mapping.error.empty())
				{
					fail("cannot classify a target-face boundary occurrence: " + mapping.error);
					return result;
				}
				if (!mapping.source_ranges.empty())
				{
					boundaries.push_back({ occurrence,static_cast<std::uint32_t>(occurrence_id),
						mapping.source_ranges,mapping.mapped_ranges,
						mapping.source_parameter_tolerance,
						mapping.candidate_parameter_tolerance,mapping.linear_tolerance });
				}
			}
			if(common.HasWarnings())
			{
				// UnableToOrientTheShape can accompany either an empty result or a
				// partially emitted set of otherwise valid one-dimensional fragments.
				// Independently audit every warning-bearing Common, not just an empty
				// one, so a retained result edge cannot hide another omitted span.
				std::string fallback_error;
				std::vector<Range> fallback_ranges;
				if(!edge_face_fallback_ranges(audit_source,target_face,audit_first,
					audit_last,audit_mapping_tolerance,fallback_ranges,fallback_error))
				{
					fail("exact trimmed edge/face Common emitted an orientation warning; "
						"independent EdgeFace audit failed: "+fallback_error);
					return result;
				}
				fallback_ranges=exact_range_union(std::move(fallback_ranges));
				for (const Range& common_range : covered_ranges)
					if (!exact_ranges_cover(fallback_ranges,common_range.first,
						common_range.second))
					{
						fail("warning-bearing Boolean Common contains a span not confirmed "
							"by the independent native EdgeFace audit");
						return result;
					}
				// EdgeFace returns source ranges rather than result edges.  Their
				// certificate is therefore exactly the already-carried native source
				// edge and target-face bounds; there is no separate operation-result
				// tolerance to invent or spread to other face uses.
				for (const Range& range : fallback_ranges)
				{
					covered_ranges.push_back(range);
					common_fragments.push_back({range, 0.0});
				}
			}
			result.normalized_parameter_tolerance=audit_mapping_tolerance/audit_span;
			if(covered_ranges.empty())
			{
				result.coverage = OcctTrimmedEdgeCoverage::none;
				return result;
			}
			// OCCT can emit duplicate 1-D fragments, together with
			// BOPAlgo_AlertUnableToOrientTheShape, when a source edge lies on a trim
			// boundary.  Orientation is irrelevant to interval coverage.  Every fragment
			// above has already mapped uniquely and completely to the finite source edge,
			// so unioning overlapping ranges is an exact set operation rather than a
			// tolerance-based repair. Any other warning still failed closed above.
			covered_ranges=exact_range_union(std::move(covered_ranges));

			std::vector<double> breakpoints;
			for (const Range& range : covered_ranges)
			{
				breakpoints.push_back(range.first);
				breakpoints.push_back(range.second);
				// Preserve every operation-result fragment boundary.  This prevents
				// a loose result edge on one fragment from donating its bound to an
				// adjacent subspan certified by a different result edge.
				for (const CommonCoverageFragment& fragment : common_fragments)
				{
					const double fragment_begin = std::max(range.first,
						fragment.source.first);
					const double fragment_end = std::min(range.second,
						fragment.source.second);
					if (positive_parameter_span({fragment_begin,fragment_end}))
					{
						breakpoints.push_back(fragment_begin);
						breakpoints.push_back(fragment_end);
					}
				}
				for (const BoundaryCoverage& boundary : boundaries)
					for (const Range& boundary_range : boundary.ranges)
					{
						const double begin = std::max(range.first, boundary_range.first);
						const double end = std::min(range.second, boundary_range.second);
						if (positive_parameter_span({begin,end}))
						{
							breakpoints.push_back(begin);
							breakpoints.push_back(end);
						}
					}
			}
			std::sort(breakpoints.begin(), breakpoints.end());
			breakpoints.erase(std::unique(breakpoints.begin(),breakpoints.end(),
				[](double a,double b){return same_parameter_endpoint(a,b);}),
				breakpoints.end());

			const double inverse_range=1.0/audit_span;
			for (std::size_t point = 1; point < breakpoints.size(); ++point)
			{
				const double begin = breakpoints[point - 1], end = breakpoints[point];
				if (!positive_parameter_span({begin,end})) continue;
				const double midpoint = 0.5 * (begin + end);
				if (!contains_parameter(covered_ranges,midpoint)) continue;

				std::vector<const BoundaryCoverage*> occurrences;
				for (const BoundaryCoverage& boundary : boundaries)
					if (contains_parameter(boundary.ranges,midpoint))
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
				interval.begin=std::clamp((begin-audit_first)*inverse_range,0.0,1.0);
				interval.end=std::clamp((end-audit_first)*inverse_range,0.0,1.0);
				interval.location = location;
				interval.target_sector_count = sectors;
				interval.boundary_occurrence_count = static_cast<std::uint8_t>(
					occurrences.size());
				bool has_operation_certificate = false;
				for (const CommonCoverageFragment& fragment : common_fragments)
				{
					// Aggregate only fragments with positive parameter overlap. A
					// shared endpoint contributes no uncertainty to either neighbour.
					// Native result-edge bounds stay local after exact interval topology
					// has already been established; they never expand that topology.
					if (!(std::min(end, fragment.source.second)
						> std::max(begin, fragment.source.first))) continue;
					interval.exact_operation_tolerance = std::max(
						interval.exact_operation_tolerance,
						fragment.exact_operation_tolerance);
					has_operation_certificate = true;
				}
				if (!has_operation_certificate)
				{
					fail("covered source interval lost its exact operation certificate");
					return result;
				}
				for(std::size_t occurrence_id=0;occurrence_id<occurrences.size();++occurrence_id)
				{
					const BoundaryCoverage& occurrence=*occurrences[occurrence_id];
					std::vector<const EdgeMapping::MappedRange*> mappings;
					for(const auto& mapping:occurrence.mapped_ranges)
						if(contains_parameter({mapping.source},midpoint))
							mappings.push_back(&mapping);
					if(mappings.size()!=1)
					{
						fail("target boundary occurrence does not have one exact parameter mapping for an interval");
						return result;
					}
					const auto& mapping=*mappings.front();
					interval.target_boundary_occurrence_ids[occurrence_id]=occurrence.occurrence_id;
					interval.target_boundary_orientations[occurrence_id]=static_cast<std::int8_t>(
						occurrence.occurrence.Orientation());
					std::string endpoint_error;
					if(!unique_target_parameter_at_source(audit_source_curve,begin,
						occurrence.occurrence,mapping,occurrence.source_parameter_tolerance,
						occurrence.target_parameter_tolerance,
						occurrence.linear_tolerance,
						interval.target_parameter_begin[occurrence_id],endpoint_error)
						||!unique_target_parameter_at_source(audit_source_curve,end,
							occurrence.occurrence,mapping,occurrence.source_parameter_tolerance,
							occurrence.target_parameter_tolerance,
							occurrence.linear_tolerance,
							interval.target_parameter_end[occurrence_id],endpoint_error))
					{
						fail("cannot recover target parameters at emitted interval endpoints: "
							+endpoint_error);
						return result;
					}
					interval.target_mapping_source_begin[occurrence_id]=first
						+(mapping.source.first-audit_first)*inverse_range*source_span;
					interval.target_mapping_source_end[occurrence_id]=first
						+(mapping.source.second-audit_first)*inverse_range*source_span;
				}
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
					&& result.intervals.back().exact_operation_tolerance
						== interval.exact_operation_tolerance
					&&same_parameter_endpoint(result.intervals.back().end,interval.begin))
					result.intervals.back().end = interval.end;
				else
					result.intervals.push_back(interval);
			}

			if (result.intervals.empty())
			{
				fail("positive-length Common coverage produced no classified source interval");
				return result;
			}
			result.coverage=exact_ranges_cover(covered_ranges,audit_first,audit_last)
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
