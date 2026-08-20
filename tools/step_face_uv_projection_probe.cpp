// Manufactured periodic/seamed OCCT regressions for the production STEP
// contact-chain UV projector. No external STEP fixture is required.

#include "core/geometry/step_face_uv_projection.h"

#include <BRep_Builder.hxx>
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
#include <gp_Vec.hxx>
#include <gp.hxx>

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <utility>
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
		const std::vector<double> boundary_parameters = {u_min,
			1.8 * std::numbers::pi, 1.95 * std::numbers::pi,
			2.05 * std::numbers::pi, 2.2 * std::numbers::pi, u_max};
		const std::vector<gp_Pnt> boundary_points = points_at(face,
			boundary_parameters, kHeight);
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

	bool combined_native_boundary_tolerance_regression()
	{
		const double u_min = 1.5 * std::numbers::pi;
		const double u_max = 2.5 * std::numbers::pi;
		const TopoDS_Face face = cylindrical_face(u_min, u_max);
		TopoDS_Edge boundary = top_boundary_occurrence(face);
		if (face.IsNull() || boundary.IsNull())
		{
			std::printf("[combined native boundary tolerance] FAIL: construction\n");
			return false;
		}
		constexpr double source_native_tolerance = 3.0e-6;
		constexpr double target_native_tolerance = 2.0e-6;
		BRep_Builder builder;
		builder.UpdateEdge(boundary, target_native_tolerance);
		const double actual_target_tolerance = BRep_Tool::Tolerance(boundary)
			* std::abs(boundary.Location().Transformation().ScaleFactor());
		const double actual_face_tolerance = BRep_Tool::Tolerance(face)
			* std::abs(face.Location().Transformation().ScaleFactor());
		const double combined_tolerance = source_native_tolerance
			+ actual_target_tolerance + actual_face_tolerance;
		const double accepted_offset = 4.0e-6;
		if (!(accepted_offset > source_native_tolerance)
			|| !(accepted_offset < combined_tolerance))
		{
			std::printf("[combined native boundary tolerance] FAIL: invalid test bounds\n");
			return false;
		}

		const std::vector<double> parameters = {u_min, 1.8 * std::numbers::pi,
			2.2 * std::numbers::pi, u_max};
		std::vector<gp_Pnt> points = points_at(face, parameters, kHeight);
		for (gp_Pnt& point : points) point.Translate(gp_Vec(0.0, 0.0, accepted_offset));
		std::vector<std::array<double, 2>> uv, replay_uv;
		std::string error, replay_error;
		bool used_exact = false, replay_used_exact = false;
		const bool projected = paracfd::core::detail::project_trim_valid_face_uv(
			points, face, {boundary}, true, source_native_tolerance, uv, error,
			&used_exact);
		const bool replayed = paracfd::core::detail::project_trim_valid_face_uv(
			points, face, {boundary}, true, source_native_tolerance, replay_uv,
			replay_error, &replay_used_exact);
		bool accepted = projected && replayed && used_exact && replay_used_exact
			&& error.empty() && replay_error.empty() && uv == replay_uv
			&& uv.size() == points.size();
		const std::string acceptance_error = error.empty() ? replay_error : error;
		for (std::size_t sample = 0; accepted && sample < uv.size(); ++sample)
		{
			const double residual = world_value(face, uv[sample][0], uv[sample][1])
				.Distance(points[sample]);
			accepted = residual > source_native_tolerance
				&& residual <= combined_tolerance;
		}

		std::vector<gp_Pnt> rejected_points = points_at(face, parameters, kHeight);
		for (gp_Pnt& point : rejected_points)
			point.Translate(gp_Vec(0.0, 0.0, 2.0 * combined_tolerance));
		uv.clear();
		error.clear();
		const bool rejected = !paracfd::core::detail::project_trim_valid_face_uv(
			rejected_points, face, {boundary}, true, source_native_tolerance, uv, error)
			&& uv.empty()
			&& error.find("curve-residual-min=") != std::string::npos
			&& error.find("trim-state IN/ON/OUT/UNKNOWN=") != std::string::npos
			&& error.find("surface-residual-min=") != std::string::npos
			&& error.find("curve-surface-residual-min=") != std::string::npos
			&& error.find("sample-limit=") != std::string::npos;
		const bool ok = accepted && rejected;
		std::printf("[combined native boundary tolerance] %s; source=%g target=%g "
			"combined=%g accepted-offset=%g\n", ok ? "PASS" : "FAIL",
			source_native_tolerance, actual_target_tolerance, combined_tolerance,
			accepted_offset);
		if (!accepted)
			std::printf("  acceptance diagnostic: %s\n", acceptance_error.c_str());
		if (!rejected) std::printf("  rejection diagnostic: %s\n", error.c_str());
		return ok;
	}

	bool per_sample_tolerance_isolation_regression()
	{
		const double u_min = 1.5 * std::numbers::pi;
		const double u_max = 2.5 * std::numbers::pi;
		const TopoDS_Face face = cylindrical_face(u_min, u_max);
		const TopoDS_Edge boundary = top_boundary_occurrence(face);
		if (face.IsNull() || boundary.IsNull())
		{
			std::printf("[per-sample tolerance isolation] FAIL: construction\n");
			return false;
		}
		constexpr double loose_tolerance = 1.0e-4;
		constexpr double tight_tolerance = 1.0e-8;
		constexpr double offset = 4.0e-5;
		const std::vector<double> parameters = {1.8 * std::numbers::pi,
			2.2 * std::numbers::pi};

		auto path_is_isolated = [&](const char* label, std::vector<gp_Pnt> points,
			bool boundary_use, const std::vector<TopoDS_Edge>& boundary_edges)
		{
			std::vector<std::array<double, 2>> uv, scalar_uv;
			std::string error, scalar_error;
			bool scalar_used_exact = false;
			const bool all_loose =
				paracfd::core::detail::project_trim_valid_face_uv(points, face,
					boundary_edges, boundary_use, loose_tolerance, scalar_uv,
					scalar_error, &scalar_used_exact);
			const bool loose_then_tight =
				paracfd::core::detail::project_trim_valid_face_uv(points, face,
					boundary_edges, boundary_use,
					std::vector<double>{loose_tolerance, tight_tolerance}, uv, error);
			const bool rejected_second = !loose_then_tight && uv.empty()
				&& error.find("sample 1") != std::string::npos;
			const std::string second_error = error;

			uv.clear();
			error.clear();
			const bool tight_then_loose =
				paracfd::core::detail::project_trim_valid_face_uv(points, face,
					boundary_edges, boundary_use,
					std::vector<double>{tight_tolerance, loose_tolerance}, uv, error);
			const bool rejected_first = !tight_then_loose && uv.empty()
				&& error.find("sample 0") != std::string::npos;
			const bool exact_path = !boundary_use || scalar_used_exact;
			const bool ok = all_loose && scalar_error.empty()
				&& scalar_uv.size() == points.size() && exact_path
				&& rejected_second && rejected_first;
			std::printf("[%s per-sample tolerance isolation] %s\n", label,
				ok ? "PASS" : "FAIL");
			if (!all_loose)
				std::printf("  all-loose diagnostic: %s\n", scalar_error.c_str());
			if (!rejected_second)
				std::printf("  loose/tight diagnostic: %s\n", second_error.c_str());
			if (!rejected_first)
				std::printf("  tight/loose diagnostic: %s\n", error.c_str());
			return ok;
		};

		std::vector<gp_Pnt> interior_points = points_at(face, parameters,
			0.5 * kHeight);
		for (gp_Pnt& point : interior_points)
		{
			gp_Vec radial(point.X(), point.Y(), 0.0);
			radial.Normalize();
			point.Translate(radial.Multiplied(offset));
		}
		std::vector<gp_Pnt> boundary_points = points_at(face, parameters, kHeight);
		for (gp_Pnt& point : boundary_points)
			point.Translate(gp_Vec(0.0, 0.0, offset));

		return path_is_isolated("interior", std::move(interior_points), false, {})
			&& path_is_isolated("exact boundary", std::move(boundary_points), true,
				{boundary});
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
	const bool native_tolerance = combined_native_boundary_tolerance_regression();
	const bool sample_isolation = per_sample_tolerance_isolation_regression();
	const bool failures = full_period_failure_regressions();
	std::printf("step face UV projection probe: %s\n",
		success && native_tolerance && sample_isolation && failures ? "PASS" : "FAIL");
	return success && native_tolerance && sample_isolation && failures ? 0 : 1;
}
