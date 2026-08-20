#include "core/geometry/occt_trimmed_edge_face.h"
#include "core/geometry/step_import.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBndLib.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Line.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <NCollection_IndexedMap.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace
{
	using paracfd::core::OcctFaceIntervalLocation;
	using paracfd::core::OcctTrimmedEdgeCoverage;
	using paracfd::core::OcctTrimmedEdgeFaceCommonResult;

	bool near(double a, double b, double tolerance = 2.0e-10)
	{
		return std::abs(a - b) <= tolerance;
	}

	void require(bool condition, const std::string& message)
	{
		if (condition) return;
		std::fprintf(stderr, "occt_trimmed_edge_face_probe: %s\n", message.c_str());
		std::exit(1);
	}

	TopoDS_Wire rectangle_wire(double x0, double y0, double x1, double y1,
		double z, bool clockwise)
	{
		std::vector<gp_Pnt> points{
			{ x0, y0, z }, { x1, y0, z }, { x1, y1, z }, { x0, y1, z }
		};
		if (clockwise) std::reverse(points.begin(), points.end());
		BRepBuilderAPI_MakePolygon polygon;
		for (const gp_Pnt& point : points) polygon.Add(point);
		polygon.Close();
		require(polygon.IsDone(), "could not manufacture a rectangle wire");
		return polygon.Wire();
	}

	TopoDS_Face rectangle_face(double scale, bool hole)
	{
		const gp_Pln plane(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
		BRepBuilderAPI_MakeFace face(plane,
			rectangle_wire(0.0, 0.0, 10.0 * scale, 10.0 * scale, 0.0, false), true);
		require(face.IsDone(), "could not manufacture the outer face");
		if (hole)
			face.Add(rectangle_wire(4.0 * scale, 4.0 * scale, 6.0 * scale,
				6.0 * scale, 0.0, true));
		require(face.IsDone(), "could not manufacture the face hole");
		return face.Face();
	}

	TopoDS_Face split_boundary_rectangle_face()
	{
		BRepBuilderAPI_MakePolygon polygon;
		for (const gp_Pnt& point : std::vector<gp_Pnt>{
			{ 0.0, 0.0, 0.0 }, { 5.0, 0.0, 0.0 }, { 10.0, 0.0, 0.0 },
			{ 10.0, 10.0, 0.0 }, { 0.0, 10.0, 0.0 } })
			polygon.Add(point);
		polygon.Close();
		require(polygon.IsDone(), "could not manufacture a split-boundary wire");
		BRepBuilderAPI_MakeFace face(gp_Pln(gp_Pnt(0.0, 0.0, 0.0),
			gp_Dir(0.0, 0.0, 1.0)), polygon.Wire(), true);
		require(face.IsDone(), "could not manufacture a split-boundary face");
		return face.Face();
	}

	TopoDS_Edge line_edge(double x0, double y0, double x1, double y1,
		double native_tolerance = 0.0)
	{
		TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(gp_Pnt(x0, y0, 0.0),
			gp_Pnt(x1, y1, 0.0));
		if (native_tolerance > 0.0)
		{
			BRep_Builder builder;
			builder.UpdateEdge(edge, native_tolerance);
		}
		return edge;
	}

	void require_valid(const OcctTrimmedEdgeFaceCommonResult& result,
		const std::string& label)
	{
		if (result.valid()) return;
		std::string message = label + " failed:";
		for (const std::string& error : result.errors) message += " " + error;
		require(false, message);
	}

	void check_interval(const OcctTrimmedEdgeFaceCommonResult& result, std::size_t index,
		double begin, double end, OcctFaceIntervalLocation location, std::uint8_t sectors,
		std::uint8_t occurrences)
	{
		require(index < result.intervals.size(), "missing expected interval");
		const auto& interval = result.intervals[index];
		const double endpoint_tolerance=std::max(2.0e-10,
			4.0*result.normalized_parameter_tolerance);
		if(!near(interval.begin,begin,endpoint_tolerance)
			||!near(interval.end,end,endpoint_tolerance))
		{
			char detail[256];std::snprintf(detail,sizeof(detail),
				"normalized interval endpoints differ: got [%.17g,%.17g], expected [%.17g,%.17g]",
				interval.begin,interval.end,begin,end);
			require(false,detail);
		}
		require(interval.location == location && interval.target_sector_count == sectors
			&& interval.boundary_occurrence_count == occurrences,
			"interval face-sector classification is incorrect");
		require(std::isfinite(interval.exact_operation_tolerance)
			&& interval.exact_operation_tolerance >= 0.0,
			"interval has no finite non-negative exact-operation certificate");
		const std::uint32_t unused=std::numeric_limits<std::uint32_t>::max();
		for(std::size_t occurrence=0;occurrence<2;++occurrence)
		{
			if(occurrence<occurrences)
			{
				require(interval.target_boundary_occurrence_ids[occurrence]!=unused,
					"active boundary interval did not retain its target occurrence ID");
				require(interval.target_boundary_orientations[occurrence]
					==static_cast<std::int8_t>(TopAbs_FORWARD)
					||interval.target_boundary_orientations[occurrence]
					==static_cast<std::int8_t>(TopAbs_REVERSED),
					"active boundary occurrence has no usable orientation");
				require(std::isfinite(interval.target_parameter_begin[occurrence])
					&&std::isfinite(interval.target_parameter_end[occurrence])
					&&std::isfinite(interval.target_mapping_source_begin[occurrence])
					&&std::isfinite(interval.target_mapping_source_end[occurrence]),
					"active boundary occurrence did not retain finite parameter mapping");
				const double raw_begin=result.source_parameter_first
					+begin*(result.source_parameter_last-result.source_parameter_first);
				const double raw_end=result.source_parameter_first
					+end*(result.source_parameter_last-result.source_parameter_first);
				const double tolerance=result.normalized_parameter_tolerance
					*(result.source_parameter_last-result.source_parameter_first);
				require(interval.target_mapping_source_begin[occurrence]<=raw_begin+tolerance
					&&interval.target_mapping_source_end[occurrence]>=raw_end-tolerance,
					"target occurrence mapping does not cover the emitted source interval");
			}
			else require(interval.target_boundary_occurrence_ids[occurrence]==unused,
				"inactive boundary occurrence slot must retain its sentinel identity");
		}
	}

	void full_interior_test()
	{
		const auto result = paracfd::core::exact_trimmed_edge_face_common(
			line_edge(1.0, 5.0, 9.0, 5.0), rectangle_face(1.0, false));
		require_valid(result, "full interior");
		require(result.coverage == OcctTrimmedEdgeCoverage::full,
			"interior edge should have full coverage");
		require(result.intervals.size() == 1, "full interior should be one interval");
		check_interval(result, 0, 0.0, 1.0, OcctFaceIntervalLocation::interior, 2, 0);
		require(result.intervals[0].exact_operation_tolerance > 0.0,
			"BRepAlgo Common result-edge tolerance was not retained on its covered subspan");
	}

	void partial_gap_test()
	{
		const auto result = paracfd::core::exact_trimmed_edge_face_common(
			line_edge(1.0, 5.0, 9.0, 5.0), rectangle_face(1.0, true));
		require_valid(result, "partial with hole");
		require(result.coverage == OcctTrimmedEdgeCoverage::partial,
			"a face hole must leave partial source-edge coverage");
		require(result.intervals.size() == 2, "face hole should create two intervals");
		check_interval(result, 0, 0.0, 3.0 / 8.0,
			OcctFaceIntervalLocation::interior, 2, 0);
		check_interval(result, 1, 5.0 / 8.0, 1.0,
			OcctFaceIntervalLocation::interior, 2, 0);
	}

	void meshed_partial_gap_test()
	{
		TopoDS_Face face=rectangle_face(1.0,true);
		BRepMesh_IncrementalMesh mesh(face,0.25);mesh.Perform();
		const auto result=paracfd::core::exact_trimmed_edge_face_common(
			line_edge(1.0,5.0,9.0,5.0),face);
		require_valid(result,"meshed partial with hole");
		require(result.coverage==OcctTrimmedEdgeCoverage::partial
			&&result.intervals.size()==2,
			"visualization meshing must not remove a disjoint exact Common span");
		check_interval(result,0,0.0,3.0/8.0,
			OcctFaceIntervalLocation::interior,2,0);
		check_interval(result,1,5.0/8.0,1.0,
			OcctFaceIntervalLocation::interior,2,0);
	}

	TopoDS_Face unrelated_loose_boundary_with_gap_face(double gap)
	{
		const gp_Pnt a(0.0,0.0,0.0),b(2.0,0.0,0.0),c(2.0,-1.0,0.0);
		const gp_Pnt d(10.0,-1.0,0.0),e(10.0,1.0,0.0),f(0.0,1.0,0.0);
		const TopoDS_Edge loose_boundary=line_edge(a.X(),a.Y(),b.X(),b.Y(),1.0e-3);
		BRepBuilderAPI_MakeWire outer;
		outer.Add(loose_boundary);
		outer.Add(BRepBuilderAPI_MakeEdge(b,c));
		outer.Add(BRepBuilderAPI_MakeEdge(c,d));
		outer.Add(BRepBuilderAPI_MakeEdge(d,e));
		outer.Add(BRepBuilderAPI_MakeEdge(e,f));
		outer.Add(BRepBuilderAPI_MakeEdge(f,a));
		require(outer.IsDone(),"could not manufacture loose-boundary outer wire");
		BRepBuilderAPI_MakeFace face(gp_Pln(gp_Pnt(0.0,0.0,0.0),
			gp_Dir(0.0,0.0,1.0)),outer.Wire(),true);
		require(face.IsDone(),"could not manufacture loose-boundary face");
		face.Add(rectangle_wire(5.0-0.5*gap,-0.5,5.0+0.5*gap,0.5,0.0,true));
		require(face.IsDone(),"could not add narrow exact hole to loose-boundary face");
		return face.Face();
	}

	void unrelated_loose_occurrence_preserves_gap_test()
	{
		constexpr double gap=4.0e-4;
		const auto result=paracfd::core::exact_trimmed_edge_face_common(
			line_edge(-1.0,0.0,11.0,0.0),unrelated_loose_boundary_with_gap_face(gap));
		require_valid(result,"unrelated loose boundary with narrow exact hole");
		require(result.coverage==OcctTrimmedEdgeCoverage::partial,
			"a loose target occurrence must not promote a positive hole gap to full coverage");
		require(result.intervals.size()==3,
			"narrow hole gap and exact boundary handoff must both remain represented");
		check_interval(result,0,1.0/12.0,3.0/12.0,
			OcctFaceIntervalLocation::boundary,1,1);
		check_interval(result,1,3.0/12.0,(6.0-0.5*gap)/12.0,
			OcctFaceIntervalLocation::interior,2,0);
		check_interval(result,2,(6.0+0.5*gap)/12.0,11.0/12.0,
			OcctFaceIntervalLocation::interior,2,0);
		// The final interval above reaches x=10; this explicit endpoint assertion makes
		// the sub-tolerance positive gap the focus rather than the outside source tails.
		require(result.intervals[1].end<result.intervals[2].begin,
			"positive hole gap was merged by an unrelated occurrence tolerance");
	}

	OcctTrimmedEdgeFaceCommonResult boundary_result(double scale, bool reversed,
		double native_tolerance)
	{
		TopoDS_Edge source = line_edge(-2.0 * scale, 0.0, 12.0 * scale, 0.0,
			native_tolerance);
		if (reversed) source.Reverse();
		return paracfd::core::exact_trimmed_edge_face_common(source,
			rectangle_face(scale, false));
	}

	void boundary_orientation_and_scale_test()
	{
		const auto baseline = boundary_result(1.0, false, 1.0e-6);
		require_valid(baseline, "boundary baseline");
		require(baseline.coverage == OcctTrimmedEdgeCoverage::partial,
			"long source edge should be only partially covered by the face boundary");
		require(baseline.intervals.size() == 1, "boundary overlap should be one interval");
		check_interval(baseline, 0, 1.0 / 7.0, 6.0 / 7.0,
			OcctFaceIntervalLocation::boundary, 1, 1);

		const auto reversed = boundary_result(1.0, true, 1.0e-6);
		require_valid(reversed, "reversed boundary");
		check_interval(reversed, 0, 1.0 / 7.0, 6.0 / 7.0,
			OcctFaceIntervalLocation::boundary, 1, 1);

		const auto scaled = boundary_result(1000.0, false, 1.0e-3);
		require_valid(scaled, "scaled boundary with native tolerance");
		check_interval(scaled, 0, 1.0 / 7.0, 6.0 / 7.0,
			OcctFaceIntervalLocation::boundary, 1, 1);
	}

	void source_shorter_than_boundary_test()
	{
		const auto result=paracfd::core::exact_trimmed_edge_face_common(
			line_edge(2.0,0.0,8.0,0.0),rectangle_face(1.0,false));
		require_valid(result,"source shorter than target boundary occurrence");
		require(result.coverage==OcctTrimmedEdgeCoverage::full
			&&result.intervals.size()==1,
			"source clipped inside a longer target boundary must be fully covered");
		check_interval(result,0,0.0,1.0,OcctFaceIntervalLocation::boundary,1,1);
	}

	void located_face_boundary_test()
	{
		gp_Trsf transform;
		transform.SetTranslation(gp_Vec(37.0, -12.0, 4.0));
		const TopLoc_Location location(transform);
		const TopoDS_Face face = TopoDS::Face(rectangle_face(1.0, false).Moved(location));
		const TopoDS_Edge source = TopoDS::Edge(
			line_edge(-2.0, 0.0, 12.0, 0.0).Moved(location));
		const auto result = paracfd::core::exact_trimmed_edge_face_common(source, face);
		require_valid(result, "located face boundary");
		require(result.coverage == OcctTrimmedEdgeCoverage::partial
			&& result.intervals.size() == 1,
			"located face boundary should retain one partial interval");
		check_interval(result, 0, 1.0 / 7.0, 6.0 / 7.0,
			OcctFaceIntervalLocation::boundary, 1, 1);
	}

	void split_collinear_boundary_test()
	{
		const auto result = paracfd::core::exact_trimmed_edge_face_common(
			line_edge(-2.0, 0.0, 12.0, 0.0), split_boundary_rectangle_face());
		require_valid(result, "split collinear boundary");
		require(result.coverage == OcctTrimmedEdgeCoverage::partial,
			"split boundary should partially cover the longer source edge");
		require(result.intervals.size() == 2,
			"adjacent trim occurrences must retain their shared pcurve breakpoint");
		check_interval(result, 0, 1.0 / 7.0, 0.5,
			OcctFaceIntervalLocation::boundary, 1, 1);
		check_interval(result, 1, 0.5, 6.0 / 7.0,
			OcctFaceIntervalLocation::boundary, 1, 1);
		require(result.intervals[0].target_boundary_occurrence_ids[0]
			!=result.intervals[1].target_boundary_occurrence_ids[0],
			"split trim intervals must retain distinct target occurrence identities");
	}

	void exact_shared_vertex_adjacency_test()
	{
		const auto result=paracfd::core::exact_trimmed_edge_face_common(
			line_edge(0.0,0.0,10.0,0.0),split_boundary_rectangle_face());
		require_valid(result,"exact shared-vertex boundary adjacency");
		require(result.coverage==OcctTrimmedEdgeCoverage::full,
			"two exactly adjacent trim occurrences must cover the full source edge");
		require(result.intervals.size()==2,
			"shared topology must retain the exact pcurve handoff without opening a gap");
		check_interval(result,0,0.0,0.5,
			OcctFaceIntervalLocation::boundary,1,1);
		check_interval(result,1,0.5,1.0,
			OcctFaceIntervalLocation::boundary,1,1);
		require(result.intervals[0].end==result.intervals[1].begin,
			"authoritative shared-vertex adjacency did not retain one exact endpoint");
	}

	TopoDS_Edge parameterized_line_edge(double first,double last)
	{
		TColgp_Array1OfPnt poles(1,2);poles.SetValue(1,gp_Pnt(-2.0,5.0,0.0));
		poles.SetValue(2,gp_Pnt(12.0,5.0,0.0));
		TColStd_Array1OfReal knots(1,2);knots.SetValue(1,first);knots.SetValue(2,last);
		TColStd_Array1OfInteger multiplicities(1,2);
		multiplicities.SetValue(1,2);multiplicities.SetValue(2,2);
		const Handle(Geom_BSplineCurve) curve=new Geom_BSplineCurve(poles,knots,
			multiplicities,1);
		return BRepBuilderAPI_MakeEdge(curve,first,last);
	}

	Handle(Geom_BSplineCurve) nonlinear_line_curve(double x0,double control_x,double x1,
		double first,double last)
	{
		TColgp_Array1OfPnt poles(1,3);
		poles.SetValue(1,gp_Pnt(x0,0.0,0.0));
		poles.SetValue(2,gp_Pnt(control_x,0.0,0.0));
		poles.SetValue(3,gp_Pnt(x1,0.0,0.0));
		TColStd_Array1OfReal knots(1,2);knots.SetValue(1,first);knots.SetValue(2,last);
		TColStd_Array1OfInteger multiplicities(1,2);
		multiplicities.SetValue(1,3);multiplicities.SetValue(2,3);
		return new Geom_BSplineCurve(poles,knots,multiplicities,2);
	}

	struct NonlinearBoundaryFixture
	{
		TopoDS_Edge source;
		TopoDS_Face target;
	};

	NonlinearBoundaryFixture clipped_nonlinear_boundary_fixture()
	{
		const Handle(Geom_BSplineCurve) curve=nonlinear_line_curve(0.0,0.5,10.0,
			10.0,20.0);
		auto parameter_at_x=[](double x)
		{return 10.0+10.0*(-1.0+std::sqrt(1.0+36.0*x))/18.0;};
		BRepBuilderAPI_MakeWire wire;
		wire.Add(BRepBuilderAPI_MakeEdge(curve,10.0,20.0));
		wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(10.0,0.0,0.0),gp_Pnt(10.0,10.0,0.0)));
		wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(10.0,10.0,0.0),gp_Pnt(0.0,10.0,0.0)));
		wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0.0,10.0,0.0),gp_Pnt(0.0,0.0,0.0)));
		require(wire.IsDone(),"could not manufacture nonlinear boundary wire");
		BRepBuilderAPI_MakeFace face(gp_Pln(gp_Pnt(0.0,0.0,0.0),
			gp_Dir(0.0,0.0,1.0)),wire.Wire(),true);
		require(face.IsDone(),"could not manufacture nonlinear boundary face");
		return {BRepBuilderAPI_MakeEdge(curve,parameter_at_x(2.0),parameter_at_x(8.0)),
			face.Face()};
	}

	std::vector<TopoDS_Edge> face_boundary_occurrences(const TopoDS_Face& face)
	{
		std::vector<TopoDS_Edge> result;
		for(TopExp_Explorer wires(face,TopAbs_WIRE);wires.More();wires.Next())
			for(TopoDS_Iterator child(wires.Current(),true,true);child.More();child.Next())
				if(child.Value().ShapeType()==TopAbs_EDGE)
					result.push_back(TopoDS::Edge(child.Value()));
		return result;
	}

	void clipped_nonlinear_target_mapping_test()
	{
		const NonlinearBoundaryFixture fixture=clipped_nonlinear_boundary_fixture();
		const TopoDS_Edge& source=fixture.source;
		const TopoDS_Face& target=fixture.target;
		const auto result=paracfd::core::exact_trimmed_edge_face_common(source,target);
		require_valid(result,"clipped nonlinear reciprocal parameterization");
		require(result.coverage==OcctTrimmedEdgeCoverage::full
			&&result.intervals.size()==1,
			"clipped nonlinear source must remain one full target-boundary atom");
		check_interval(result,0,0.0,1.0,OcctFaceIntervalLocation::boundary,1,1);
		const auto occurrences=face_boundary_occurrences(target);
		BRepAdaptor_Curve source_curve(source);
		for(const auto& interval:result.intervals)
		{
			require(interval.boundary_occurrence_count==1,
				"nonlinear split-boundary atom must identify one target occurrence");
			const std::uint32_t occurrence_id=interval.target_boundary_occurrence_ids[0];
			require(occurrence_id<occurrences.size(),
				"nonlinear split-boundary target occurrence ID is out of range");
			BRepAdaptor_Curve target_curve(occurrences[occurrence_id]);
			const double raw_begin=result.source_parameter_first+interval.begin
				*(result.source_parameter_last-result.source_parameter_first);
			const double raw_end=result.source_parameter_first+interval.end
				*(result.source_parameter_last-result.source_parameter_first);
			require(source_curve.Value(raw_begin).Distance(
				target_curve.Value(interval.target_parameter_begin[0]))<=1.0e-6
				&&source_curve.Value(raw_end).Distance(
					target_curve.Value(interval.target_parameter_end[0]))<=1.0e-6,
				"target parameters do not correspond to the emitted clipped nonlinear endpoints");
		}
	}

	void parameterization_test()
	{
		for(const auto& [first,last,label]:std::vector<std::tuple<double,double,std::string>>{
			{0.0,1.0,"unit"},{18.0,32.0,"shifted"},
			{1.0e-12,2.0e-12,"tiny"},{1.0e9,1.0e9+0.1,"large-offset"}})
		{
			std::fprintf(stderr,"[trimmed probe]   parameter range %s\n",label.c_str());
			const auto result=paracfd::core::exact_trimmed_edge_face_common(
				parameterized_line_edge(first,last),rectangle_face(1.0,false));
			require_valid(result,label+" curve parameterization");
			require(result.coverage==OcctTrimmedEdgeCoverage::partial,
				label+" parameter source should be partially covered");
			check_interval(result,0,1.0/7.0,6.0/7.0,
				OcctFaceIntervalLocation::interior,2,0);
		}
		const double unresolved_first=1.0e12;
		const double unresolved_last=unresolved_first+0.5;
		std::fprintf(stderr,"[trimmed probe]   unresolved parameter range\n");
		const auto unresolved=paracfd::core::exact_trimmed_edge_face_common(
			parameterized_line_edge(unresolved_first,unresolved_last),rectangle_face(1.0,false));
		require(!unresolved.valid()&&unresolved.intervals.empty(),
			"numerically unresolved parameter span must fail closed");
	}

	void periodic_seam_test()
	{
		const TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(2.0, 5.0).Shape();
		TopoDS_Face lateral;
		for (TopExp_Explorer faces(cylinder, TopAbs_FACE); faces.More(); faces.Next())
		{
			const TopoDS_Face candidate = TopoDS::Face(faces.Current());
			if (BRepAdaptor_Surface(candidate).GetType() == GeomAbs_Cylinder)
			{
				lateral = candidate;
				break;
			}
		}
		require(!lateral.IsNull(), "could not find manufactured cylinder lateral face");
		TopoDS_Edge seam;
		for (TopExp_Explorer edges(lateral, TopAbs_EDGE); edges.More(); edges.Next())
		{
			const TopoDS_Edge candidate = TopoDS::Edge(edges.Current());
			if (!BRep_Tool::Degenerated(candidate) && BRep_Tool::IsClosed(candidate, lateral))
			{
				seam = candidate;
				break;
			}
		}
		require(!seam.IsNull(), "could not find manufactured periodic seam edge");
		const auto result = paracfd::core::exact_trimmed_edge_face_common(seam, lateral);
		require_valid(result, "periodic seam");
		require(result.coverage == OcctTrimmedEdgeCoverage::full
			&& result.intervals.size() == 1, "periodic seam should be one full interval");
		check_interval(result, 0, 0.0, 1.0, OcctFaceIntervalLocation::interior, 2, 2);
		require(result.intervals[0].target_boundary_occurrence_ids[0]
			!=result.intervals[0].target_boundary_occurrence_ids[1]
			&&result.intervals[0].target_boundary_orientations[0]
				!=result.intervals[0].target_boundary_orientations[1],
			"periodic seam must retain two distinct oppositely oriented pcurve occurrences");
	}

	void empty_and_failure_test()
	{
		const auto empty = paracfd::core::exact_trimmed_edge_face_common(
			line_edge(-2.0, 12.0, 12.0, 12.0), rectangle_face(1.0, false));
		require_valid(empty, "disjoint edge");
		require(empty.coverage == OcctTrimmedEdgeCoverage::none
			&& empty.intervals.empty(), "disjoint edge should be valid empty coverage");

		const auto invalid = paracfd::core::exact_trimmed_edge_face_common(
			TopoDS_Edge{}, rectangle_face(1.0, false));
		require(!invalid.valid() && invalid.intervals.empty(),
			"null input must fail closed without partial intervals");

		TopoDS_Face invalid_face;BRep_Builder builder;builder.MakeFace(invalid_face);
		const auto invalid_target=paracfd::core::exact_trimmed_edge_face_common(
			line_edge(-2.0,5.0,12.0,5.0),invalid_face);
		require(!invalid_target.valid()&&invalid_target.intervals.empty()
			&&std::isfinite(invalid_target.source_parameter_first)
			&&std::isfinite(invalid_target.source_parameter_last)
			&&invalid_target.source_parameter_last>invalid_target.source_parameter_first,
			"invalid target must fail closed while retaining the finite source range");
	}

	TopoDS_Face planar_face_from_edges(const std::vector<TopoDS_Edge>& edges,
		const char* label)
	{
		BRepBuilderAPI_MakeWire wire;
		for(const TopoDS_Edge& edge:edges)wire.Add(edge);
		require(wire.IsDone(),std::string("could not manufacture ")+label+" wire");
		BRepBuilderAPI_MakeFace face(wire.Wire(),true);
		require(face.IsDone(),std::string("could not manufacture ")+label+" face");
		return face.Face();
	}

	void manufactured_fan4_import_test()
	{
		const gp_Pnt a(2.0,0.0,0.0),b(8.0,0.0,0.0);
		const TopoDS_Edge shared=BRepBuilderAPI_MakeEdge(a,b);
		const TopoDS_Face upper=planar_face_from_edges({shared,
			BRepBuilderAPI_MakeEdge(b,gp_Pnt(8.0,0.0,3.0)),
			BRepBuilderAPI_MakeEdge(gp_Pnt(8.0,0.0,3.0),gp_Pnt(2.0,0.0,3.0)),
			BRepBuilderAPI_MakeEdge(gp_Pnt(2.0,0.0,3.0),a)},"upper owner");
		TopoDS_Edge reversed_shared=shared;reversed_shared.Reverse();
		const TopoDS_Face lower=planar_face_from_edges({reversed_shared,
			BRepBuilderAPI_MakeEdge(a,gp_Pnt(2.0,0.0,-3.0)),
			BRepBuilderAPI_MakeEdge(gp_Pnt(2.0,0.0,-3.0),gp_Pnt(8.0,0.0,-3.0)),
			BRepBuilderAPI_MakeEdge(gp_Pnt(8.0,0.0,-3.0),b)},"lower owner");
		BRepBuilderAPI_MakeFace target_builder(gp_Pln(gp_Pnt(0.0,0.0,0.0),
			gp_Dir(0.0,0.0,1.0)),rectangle_wire(0.0,-5.0,10.0,5.0,0.0,false),true);
		require(target_builder.IsDone(),"could not manufacture fan-4 target face");
		const TopoDS_Face target=target_builder.Face();
		BRep_Builder builder;TopoDS_Shell owner_shell;builder.MakeShell(owner_shell);
		builder.Add(owner_shell,upper);builder.Add(owner_shell,lower);
		TopoDS_Compound compound;builder.MakeCompound(compound);
		builder.Add(compound,owner_shell);builder.Add(compound,target);

		const auto nonce=std::chrono::high_resolution_clock::now().time_since_epoch().count();
		const std::filesystem::path path=std::filesystem::temp_directory_path()
			/("paracfd_fan4_"+std::to_string(nonce)+".step");
		STEPControl_Writer writer;
		require(writer.Transfer(compound,STEPControl_AsIs)==IFSelect_RetDone,
			"could not transfer manufactured fan-4 shape to STEP");
		require(writer.Write(path.string().c_str())==IFSelect_RetDone,
			"could not write manufactured fan-4 STEP");
		std::string error;
		const paracfd::core::StepGeometry geometry=
			paracfd::core::load_step_geometry(path.string(),0.25,&error);
		std::error_code remove_error;std::filesystem::remove(path,remove_error);
		require(!geometry.mesh.empty(),"manufactured fan-4 STEP import failed: "+error);
		const bool found=std::any_of(geometry.contacts.curves.begin(),
			geometry.contacts.curves.end(),[](const paracfd::core::StepContactCurve& curve)
		{
			unsigned sectors=0;for(const auto& use:curve.uses)sectors+=use.sector_count;
			return curve.fan_degree==4&&sectors==4&&curve.uses.size()==3;
		});
		if(!found)
		{
			std::fprintf(stderr,"[trimmed probe] fan4 graph curves=%zu unresolved=%zu\n",
				geometry.contacts.curves.size(),geometry.contacts.unresolved_edge_face_spans.size());
			for(const auto& curve:geometry.contacts.curves)
			{
				unsigned sectors=0;for(const auto& use:curve.uses)sectors+=use.sector_count;
				std::fprintf(stderr,"  curve %llu fan=%u uses=%zu sectors=%u sources=%zu\n",
					static_cast<unsigned long long>(curve.id),curve.fan_degree,
					curve.uses.size(),sectors,curve.source_edge_ids.size());
			}
			for(const auto& span:geometry.contacts.unresolved_edge_face_spans)
				std::fprintf(stderr,"  unresolved edge=%u face=%u kind=%u [%.17g,%.17g]: %s\n",
					span.source_edge_id,span.target_face_id,static_cast<unsigned>(span.kind),
					span.source_parameter_begin,span.source_parameter_end,span.reason.c_str());
		}
		require(found,"shared two-owner edge plus a target-face interior was not certified as fan 4");
	}

	struct StepPairAudit
	{
		OcctTrimmedEdgeFaceCommonResult result;
		std::array<double,6> edge_bounds{};
	};

	StepPairAudit audit_step_pair(const char* path,std::uint32_t edge_id,
		std::uint32_t face_id,bool mesh_first)
	{
		STEPControl_Reader reader;
		require(reader.ReadFile(path)==IFSelect_RetDone,"could not read requested STEP file");
		require(reader.TransferRoots()>0,"could not transfer requested STEP roots");
		const TopoDS_Shape shape=reader.OneShape();
		if(mesh_first)
		{
			BRepMesh_IncrementalMesh mesher(shape,2.0);mesher.Perform();
		}
		NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher> edges;
		TopExp::MapShapes(shape,TopAbs_EDGE,edges);
		std::vector<TopoDS_Face> faces;
		for(TopExp_Explorer face(shape,TopAbs_FACE);face.More();face.Next())
			faces.push_back(TopoDS::Face(face.Current()));
		require(edge_id<static_cast<std::uint32_t>(edges.Extent())&&face_id<faces.size(),
			"requested zero-based edge/face ID is out of range");
		const TopoDS_Edge edge=TopoDS::Edge(edges(static_cast<Standard_Integer>(edge_id+1)));
		Bnd_Box edge_box;BRepBndLib::AddOptimal(edge,edge_box,false,true);
		StepPairAudit audit;
		edge_box.Get(audit.edge_bounds[0],audit.edge_bounds[1],audit.edge_bounds[2],
			audit.edge_bounds[3],audit.edge_bounds[4],audit.edge_bounds[5]);
		audit.result=paracfd::core::exact_trimmed_edge_face_common(edge,faces[face_id]);
		return audit;
	}

	int inspect_step_pair(const char* path,std::uint32_t edge_id,std::uint32_t face_id,bool mesh_first)
	{
		const StepPairAudit audit=audit_step_pair(path,edge_id,face_id,mesh_first);
		std::printf("edge bbox=[%.17g %.17g %.17g]-[%.17g %.17g %.17g]\n",
			audit.edge_bounds[0],audit.edge_bounds[1],audit.edge_bounds[2],
			audit.edge_bounds[3],audit.edge_bounds[4],audit.edge_bounds[5]);
		const auto& result=audit.result;
		std::printf("edge=%u face=%u valid=%u coverage=%u range=[%.17g,%.17g] tol=%.17g intervals=%zu\n",
			edge_id,face_id,result.valid()?1u:0u,static_cast<unsigned>(result.coverage),
			result.source_parameter_first,result.source_parameter_last,
			result.normalized_parameter_tolerance,result.intervals.size());
		for(const auto& interval:result.intervals)
			std::printf("  normalized [%.17g,%.17g] kind=%u sectors=%u occurrences=%u\n",
				interval.begin,interval.end,static_cast<unsigned>(interval.location),
				static_cast<unsigned>(interval.target_sector_count),
				static_cast<unsigned>(interval.boundary_occurrence_count));
		for(const std::string& error:result.errors)std::printf("  error: %s\n",error.c_str());
		return result.valid()?0:1;
	}

	int compare_meshed_step_pair(const char* path,std::uint32_t edge_id,std::uint32_t face_id)
	{
		const StepPairAudit unmeshed=audit_step_pair(path,edge_id,face_id,false);
		const StepPairAudit meshed=audit_step_pair(path,edge_id,face_id,true);
		require(unmeshed.edge_bounds==meshed.edge_bounds,
			"visualization meshing changed the exact source-edge bounds");
		require_valid(unmeshed.result,"unmeshed STEP pair");
		require_valid(meshed.result,"meshed STEP pair");
		require(unmeshed.result.coverage==meshed.result.coverage
			&&unmeshed.result.intervals.size()==meshed.result.intervals.size(),
			"visualization meshing changed exact STEP interval coverage");
		const double tolerance=std::max({unmeshed.result.normalized_parameter_tolerance,
			meshed.result.normalized_parameter_tolerance,
			256.0*std::numeric_limits<double>::epsilon()});
		for(std::size_t interval=0;interval<unmeshed.result.intervals.size();++interval)
		{
			const auto& a=unmeshed.result.intervals[interval];
			const auto& b=meshed.result.intervals[interval];
			require(std::abs(a.begin-b.begin)<=tolerance
				&&std::abs(a.end-b.end)<=tolerance&&a.location==b.location
				&&a.target_sector_count==b.target_sector_count
				&&a.boundary_occurrence_count==b.boundary_occurrence_count,
				"visualization meshing changed an exact STEP interval or sector");
		}
		std::printf("occt_trimmed_edge_face_probe: meshed/unmeshed STEP pair PASS "
			"(edge=%u face=%u intervals=%zu)\n",edge_id,face_id,
			meshed.result.intervals.size());
		return 0;
	}
}

int main(int argc,char** argv)
{
	if((argc==5||argc==6)&&std::string(argv[1])=="--inspect-step")
		return inspect_step_pair(argv[2],static_cast<std::uint32_t>(std::stoul(argv[3])),
			static_cast<std::uint32_t>(std::stoul(argv[4])),argc==6&&std::string(argv[5])=="--mesh-first");
	if(argc==5&&std::string(argv[1])=="--compare-meshed-step")
		return compare_meshed_step_pair(argv[2],
			static_cast<std::uint32_t>(std::stoul(argv[3])),
			static_cast<std::uint32_t>(std::stoul(argv[4])));
	if(argc!=1)
	{
		std::fprintf(stderr,"usage: occt_trimmed_edge_face_probe "
			"[--inspect-step FILE EDGE FACE [--mesh-first] | "
			"--compare-meshed-step FILE EDGE FACE]\n");
		return 2;
	}
	auto run=[](const char* label,auto test)
	{
		std::fprintf(stderr,"[trimmed probe] %s\n",label);test();
	};
	run("full interior",full_interior_test);
	run("partial gap",partial_gap_test);
	run("meshed partial gap",meshed_partial_gap_test);
	run("unrelated loose occurrence preserves gap",
		unrelated_loose_occurrence_preserves_gap_test);
	run("boundary orientation and scale",boundary_orientation_and_scale_test);
	run("source shorter than boundary",source_shorter_than_boundary_test);
	run("located face boundary",located_face_boundary_test);
	run("split collinear boundary",split_collinear_boundary_test);
	run("exact shared-vertex adjacency",exact_shared_vertex_adjacency_test);
	run("parameterization",parameterization_test);
	run("clipped nonlinear mapping",clipped_nonlinear_target_mapping_test);
	run("periodic seam",periodic_seam_test);
	run("empty and failure",empty_and_failure_test);
	run("manufactured fan 4 import",manufactured_fan4_import_test);
	std::printf("occt_trimmed_edge_face_probe: PASS\n");
	return 0;
}
