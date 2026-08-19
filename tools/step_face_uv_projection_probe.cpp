// Manufactured periodic/seamed OCCT regressions for the production STEP
// contact-chain UV projector. No external STEP fixture is required.

#include "core/geometry/step_face_uv_projection.h"

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Precision.hxx>
#include <TopAbs_State.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pnt2d.hxx>
#include <gp.hxx>

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

namespace
{
	constexpr double kRadius = 10.0;
	constexpr double kHeight = 20.0;
	constexpr double kLinearTolerance = 1.0e-6;

	TopoDS_Face cylindrical_face(double u_min, double u_max)
	{
		const Handle(Geom_CylindricalSurface) surface =
			new Geom_CylindricalSurface(gp_Ax3(gp::Origin(), gp::DZ()), kRadius);
		BRepBuilderAPI_MakeFace make(surface, u_min, u_max, 0.0, kHeight,
			Precision::Confusion());
		return make.IsDone() ? make.Face() : TopoDS_Face{};
	}

	gp_Pnt world_value(const TopoDS_Face& face, double u, double v)
	{
		TopLoc_Location location;
		const Handle(Geom_Surface) surface = BRep_Tool::Surface(face, location);
		gp_Pnt point = surface->Value(u, v);
		point.Transform(location.Transformation());
		return point;
	}

	std::vector<gp_Pnt> points_at(const TopoDS_Face& face,
		const std::vector<double>& u_values, double v)
	{
		std::vector<gp_Pnt> result;
		for (double u : u_values) result.push_back(world_value(face, u, v));
		return result;
	}

	TopoDS_Edge top_boundary_occurrence(const TopoDS_Face& face)
	{
		for (TopExp_Explorer edge(face, TopAbs_EDGE); edge.More(); edge.Next())
		{
			const TopoDS_Edge occurrence = TopoDS::Edge(edge.Current());
			double first = 0.0, last = 0.0;
			const Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(
				occurrence, face, first, last);
			if (pcurve.IsNull() || !std::isfinite(first) || !std::isfinite(last)) continue;
			const gp_Pnt2d midpoint = pcurve->Value(0.5 * (first + last));
			if (std::abs(midpoint.Y() - kHeight) <= kLinearTolerance)
				return occurrence;
		}
		return {};
	}

	bool validate_success(const char* label, const TopoDS_Face& face,
		const std::vector<gp_Pnt>& points,
		const std::vector<std::array<double, 2>>& uv,
		bool require_on, bool require_periodic_continuity)
	{
		bool ok = uv.size() == points.size();
		for (std::size_t sample = 0; ok && sample < uv.size(); ++sample)
		{
			const double u = uv[sample][0], v = uv[sample][1];
			if (!std::isfinite(u) || !std::isfinite(v)) { ok = false; break; }
			const TopAbs_State state = BRepClass_FaceClassifier(face,
				gp_Pnt2d(u, v), 1.0e-7, true).State();
			if (require_on ? state != TopAbs_ON : state != TopAbs_IN && state != TopAbs_ON)
			{
				ok = false;
				break;
			}
			if (world_value(face, u, v).Distance(points[sample]) > kLinearTolerance)
			{
				ok = false;
				break;
			}
			if (require_periodic_continuity && sample > 0
				&& std::abs(u - uv[sample - 1][0]) >= std::numbers::pi)
			{
				ok = false;
				break;
			}
		}
		std::printf("[%s] %s; samples=%zu\n", label, ok ? "PASS" : "FAIL", uv.size());
		return ok;
	}

	bool shifted_periodic_success_regression()
	{
		const double u_min = 1.5 * std::numbers::pi;
		const double u_max = 2.5 * std::numbers::pi;
		const TopoDS_Face face = cylindrical_face(u_min, u_max);
		if (face.IsNull())
		{
			std::printf("[shifted periodic interior] FAIL: face construction\n");
			return false;
		}
		const std::vector<double> parameters = {1.8 * std::numbers::pi,
			1.95 * std::numbers::pi, 2.05 * std::numbers::pi,
			2.2 * std::numbers::pi};
		const std::vector<gp_Pnt> interior_points = points_at(face, parameters, 0.5 * kHeight);
		std::vector<std::array<double, 2>> interior_uv, replay_uv;
		std::string error, replay_error;
		const bool projected = paracfd::core::detail::project_trim_valid_face_uv(
			interior_points, face, {}, false, kLinearTolerance, interior_uv, error);
		const bool replayed = paracfd::core::detail::project_trim_valid_face_uv(
			interior_points, face, {}, false, kLinearTolerance, replay_uv, replay_error);
		bool ok = projected && replayed && error.empty() && replay_error.empty()
			&& interior_uv == replay_uv
			&& validate_success("shifted periodic interior", face, interior_points,
				interior_uv, false, true);
		if (!projected || !replayed)
			std::printf("  diagnostic: %s%s%s\n", error.c_str(),
				error.empty() || replay_error.empty() ? "" : " | ", replay_error.c_str());

		const TopoDS_Edge boundary = top_boundary_occurrence(face);
		const std::vector<gp_Pnt> boundary_points = points_at(face, parameters, kHeight);
		std::vector<std::array<double, 2>> boundary_uv;
		error.clear();
		const bool boundary_projected = !boundary.IsNull()
			&& paracfd::core::detail::project_trim_valid_face_uv(boundary_points, face,
				{boundary}, true, kLinearTolerance, boundary_uv, error);
		ok = boundary_projected
			&& validate_success("exact boundary pcurve", face, boundary_points,
				boundary_uv, true, true) && ok;
		if (!boundary_projected)
			std::printf("[exact boundary pcurve] FAIL: %s\n", boundary.IsNull()
				? "top edge occurrence not found" : error.c_str());
		return ok;
	}

	bool full_period_failure_regressions()
	{
		const TopoDS_Face face = cylindrical_face(0.0, 2.0 * std::numbers::pi);
		if (face.IsNull()) return false;
		const std::vector<double> crossing_parameters = {1.8 * std::numbers::pi,
			1.95 * std::numbers::pi, 2.05 * std::numbers::pi,
			2.2 * std::numbers::pi};
		std::vector<std::array<double, 2>> uv;
		std::string error;
		const bool crossing_projected = paracfd::core::detail::project_trim_valid_face_uv(
			points_at(face, crossing_parameters, 0.5 * kHeight), face, {}, false,
			kLinearTolerance, uv, error);
		const bool discontinuity = !crossing_projected && uv.empty()
			&& error.find("continuous branch") != std::string::npos;
		std::printf("[full-period seam crossing] %s; diagnostic=%s\n",
			discontinuity ? "PASS" : "FAIL", error.c_str());

		const std::vector<gp_Pnt> seam_points = {
			world_value(face, 0.0, 4.0), world_value(face, 0.0, 10.0),
			world_value(face, 0.0, 16.0)};
		uv.clear();
		error.clear();
		const bool seam_projected = paracfd::core::detail::project_trim_valid_face_uv(
			seam_points, face, {}, false, kLinearTolerance, uv, error);
		const bool ambiguity = !seam_projected && uv.empty()
			&& error.find("multiple") != std::string::npos;
		std::printf("[full-period seam ambiguity] %s; diagnostic=%s\n",
			ambiguity ? "PASS" : "FAIL", error.c_str());
		return discontinuity && ambiguity;
	}
}

int main()
{
	const bool success = shifted_periodic_success_regression();
	const bool failures = full_period_failure_regressions();
	std::printf("step face UV projection probe: %s\n",
		success && failures ? "PASS" : "FAIL");
	return success && failures ? 0 : 1;
}
