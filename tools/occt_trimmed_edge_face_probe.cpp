#include "core/geometry/occt_trimmed_edge_face.h"

#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Geom_Line.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
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
		require(near(interval.begin, begin) && near(interval.end, end),
			"normalized interval endpoints differ from the exact manufactured result");
		require(interval.location == location && interval.target_sector_count == sectors
			&& interval.boundary_occurrence_count == occurrences,
			"interval face-sector classification is incorrect");
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
	}

	void parameterization_test()
	{
		const Handle(Geom_Line) line = new Geom_Line(gp_Pnt(-20.0, 5.0, 0.0),
			gp_Dir(1.0, 0.0, 0.0));
		TopoDS_Edge source = BRepBuilderAPI_MakeEdge(line, 18.0, 32.0);
		const auto result = paracfd::core::exact_trimmed_edge_face_common(source,
			rectangle_face(1.0, false));
		require_valid(result, "shifted curve parameterization");
		require(result.coverage == OcctTrimmedEdgeCoverage::partial,
			"shifted parameter source should be partially covered");
		check_interval(result, 0, 1.0 / 7.0, 6.0 / 7.0,
			OcctFaceIntervalLocation::interior, 2, 0);
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
	}
}

int main()
{
	full_interior_test();
	partial_gap_test();
	boundary_orientation_and_scale_test();
	located_face_boundary_test();
	split_collinear_boundary_test();
	parameterization_test();
	periodic_seam_test();
	empty_and_failure_test();
	std::printf("occt_trimmed_edge_face_probe: PASS\n");
	return 0;
}
