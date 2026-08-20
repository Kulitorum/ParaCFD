#include "core/geometry/step_face_uv_projection.h"

#include <BRepAdaptor_Surface.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Extrema_ExtPS.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <TopAbs_State.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <utility>

namespace paracfd::core::detail
{
	namespace
	{
		struct SurfaceChart
		{
			Handle(Geom_Surface) surface;
			TopLoc_Location location;
			gp_Trsf inverse_location;
			double location_scale = 1.0;
			double u_min = 0.0;
			double u_max = 0.0;
			double v_min = 0.0;
			double v_max = 0.0;
			double numeric_u_tolerance = 0.0;
			double numeric_v_tolerance = 0.0;
			double numeric_uv_tolerance = 0.0;
			double face_native_tolerance = 0.0;
			bool u_periodic = false;
			bool v_periodic = false;
			double u_period = 0.0;
			double v_period = 0.0;
		};

		struct SampleProjectionTolerance
		{
			double geometric_tolerance = 0.0;
			double u_tolerance = 0.0;
			double v_tolerance = 0.0;
			double uv_tolerance = 0.0;
		};

		struct UvCandidate
		{
			double u = 0.0;
			double v = 0.0;
			double distance = 0.0;
		};

		using CandidateRow = std::vector<UvCandidate>;
		using CandidateTable = std::vector<CandidateRow>;
		struct CandidateAudit
		{
			std::size_t extrema = 0;
			std::size_t aliases = 0;
			std::size_t trim_valid = 0;
			std::size_t residual_valid = 0;
			std::size_t regular = 0;
			double minimum_residual = std::numeric_limits<double>::infinity();
		};

		bool finite_bounds(double lo, double hi)
		{
			return std::isfinite(lo) && std::isfinite(hi) && hi >= lo;
		}

		double representable_spacing(double value)
		{
			const double up = std::nextafter(value,
				std::numeric_limits<double>::infinity());
			const double down = std::nextafter(value,
				-std::numeric_limits<double>::infinity());
			double spacing = 0.0;
			if (std::isfinite(up)) spacing = std::max(spacing, std::abs(up - value));
			if (std::isfinite(down)) spacing = std::max(spacing, std::abs(value - down));
			return spacing;
		}

		bool initialize_chart(const TopoDS_Face& face, SurfaceChart& chart,
			std::string& error)
		{
			chart.surface = BRep_Tool::Surface(face, chart.location);
			if (chart.surface.IsNull())
			{
				error = "face has no geometric surface";
				return false;
			}
			BRepTools::UVBounds(face, chart.u_min, chart.u_max, chart.v_min, chart.v_max);
			if (!finite_bounds(chart.u_min, chart.u_max)
				|| !finite_bounds(chart.v_min, chart.v_max))
			{
				error = "face has non-finite or reversed UV bounds";
				return false;
			}
			chart.location_scale = std::abs(
				chart.location.Transformation().ScaleFactor());
			if (!std::isfinite(chart.location_scale)
				|| !(chart.location_scale > std::numeric_limits<double>::min()))
			{
				error = "face has an invalid location scale";
				return false;
			}
			chart.face_native_tolerance = BRep_Tool::Tolerance(face)
				* chart.location_scale;
			if (!std::isfinite(chart.face_native_tolerance)
				|| chart.face_native_tolerance < 0.0)
			{
				error = "face has an invalid native tolerance";
				return false;
			}
			const double u_span = chart.u_max - chart.u_min;
			const double v_span = chart.v_max - chart.v_min;
			// Surface parameters are arbitrarily scaled and translated. A fixed
			// Precision::PConfusion() floor is dimensionally invalid here and can consume
			// an entire small chart. Bound only representable endpoint spacing and
			// span-relative arithmetic, matching the face conformer's policy.
			chart.numeric_u_tolerance = 8.0 * std::max(
				representable_spacing(chart.u_min), representable_spacing(chart.u_max))
				+ 32.0 * std::numeric_limits<double>::epsilon() * u_span;
			chart.numeric_v_tolerance = 8.0 * std::max(
				representable_spacing(chart.v_min), representable_spacing(chart.v_max))
				+ 32.0 * std::numeric_limits<double>::epsilon() * v_span;
			chart.numeric_uv_tolerance = std::max(chart.numeric_u_tolerance,
				chart.numeric_v_tolerance);
			if (!(chart.numeric_u_tolerance > 0.0)
				|| !(chart.numeric_v_tolerance > 0.0))
			{
				error = "face UV chart has no representable coordinate resolution";
				return false;
			}
			chart.inverse_location = chart.location.Transformation().Inverted();
			chart.u_periodic = chart.surface->IsUPeriodic();
			chart.v_periodic = chart.surface->IsVPeriodic();
			if (chart.u_periodic)
			{
				chart.u_period = chart.surface->UPeriod();
				if (!std::isfinite(chart.u_period)
					|| !(chart.u_period > chart.numeric_uv_tolerance))
				{
					error = "face reports an invalid U period";
					return false;
				}
			}
			if (chart.v_periodic)
			{
				chart.v_period = chart.surface->VPeriod();
				if (!std::isfinite(chart.v_period)
					|| !(chart.v_period > chart.numeric_uv_tolerance))
				{
					error = "face reports an invalid V period";
					return false;
				}
			}
			return true;
		}

		bool initialize_sample_tolerances(const TopoDS_Face& face,
			const SurfaceChart& chart,
			const std::vector<double>& source_geometric_tolerances_mm,
			double additional_target_native_tolerance,
			std::vector<SampleProjectionTolerance>& tolerances, std::string& error)
		{
			tolerances.clear();
			tolerances.reserve(source_geometric_tolerances_mm.size());
			if (!std::isfinite(additional_target_native_tolerance)
				|| additional_target_native_tolerance < 0.0)
			{
				error = "target has an invalid additional native tolerance";
				return false;
			}
			BRepAdaptor_Surface adaptor(face, true);
			for (std::size_t sample = 0;
				sample < source_geometric_tolerances_mm.size(); ++sample)
			{
				const double source_tolerance = source_geometric_tolerances_mm[sample];
				if (!std::isfinite(source_tolerance) || source_tolerance < 0.0)
				{
					error = "sample " + std::to_string(sample)
						+ " has an invalid geometric tolerance";
					return false;
				}
				SampleProjectionTolerance tolerance;
				const double certified_tolerance = source_tolerance
					+ additional_target_native_tolerance
					+ chart.face_native_tolerance;
				if (!std::isfinite(certified_tolerance))
				{
					error = "sample " + std::to_string(sample)
						+ " geometric tolerance overflows";
					return false;
				}
				tolerance.geometric_tolerance = std::max(Precision::Confusion(),
					certified_tolerance);
				const double local_linear_tolerance = tolerance.geometric_tolerance
					/ chart.location_scale;
				const double resolved_u = adaptor.UResolution(local_linear_tolerance);
				const double resolved_v = adaptor.VResolution(local_linear_tolerance);
				tolerance.u_tolerance = chart.numeric_u_tolerance;
				tolerance.v_tolerance = chart.numeric_v_tolerance;
				if (std::isfinite(resolved_u) && resolved_u > 0.0)
					tolerance.u_tolerance = std::max(tolerance.u_tolerance, resolved_u);
				if (std::isfinite(resolved_v) && resolved_v > 0.0)
					tolerance.v_tolerance = std::max(tolerance.v_tolerance, resolved_v);
				tolerance.uv_tolerance = std::max(tolerance.u_tolerance,
					tolerance.v_tolerance);
				constexpr double maximum_tolerance_fraction = 1.0 / 64.0;
				if (!(tolerance.u_tolerance < maximum_tolerance_fraction
						* (chart.u_max - chart.u_min))
					|| !(tolerance.v_tolerance < maximum_tolerance_fraction
						* (chart.v_max - chart.v_min)))
				{
					error = "sample " + std::to_string(sample)
						+ " tolerance consumes a material fraction of the UV chart";
					return false;
				}
				if ((chart.u_periodic && !(chart.u_period > tolerance.uv_tolerance))
					|| (chart.v_periodic && !(chart.v_period > tolerance.uv_tolerance)))
				{
					error = "sample " + std::to_string(sample)
						+ " tolerance is too large for the periodic UV chart";
					return false;
				}
				tolerances.push_back(tolerance);
			}
			return true;
		}

		TopAbs_State trim_state(const TopoDS_Face& face, double u, double v, double tolerance)
		{
			const BRepClass_FaceClassifier classifier(face, gp_Pnt2d(u, v), tolerance, true);
			return classifier.State();
		}

		bool re_evaluates_to(const SurfaceChart& chart, const gp_Pnt& world_point,
			double u, double v, double geometric_tolerance, double& distance)
		{
			if (!std::isfinite(u) || !std::isfinite(v)) return false;
			gp_Pnt evaluated = chart.surface->Value(u, v);
			evaluated.Transform(chart.location.Transformation());
			distance = evaluated.Distance(world_point);
			return std::isfinite(distance) && distance <= geometric_tolerance;
		}

		bool regular_parameterization(const SurfaceChart& chart, double u, double v)
		{
			gp_Pnt point;
			gp_Vec derivative_u, derivative_v;
			chart.surface->D1(u, v, point, derivative_u, derivative_v);
			derivative_u.Transform(chart.location.Transformation());
			derivative_v.Transform(chart.location.Transformation());
			const double magnitude_u = derivative_u.Magnitude();
			const double magnitude_v = derivative_v.Magnitude();
			const double product = magnitude_u * magnitude_v;
			if (!std::isfinite(magnitude_u) || !std::isfinite(magnitude_v)
				|| !std::isfinite(product) || !(product > std::numeric_limits<double>::min()))
				return false;
			const double cross_magnitude = derivative_u.Crossed(derivative_v).Magnitude();
			return std::isfinite(cross_magnitude)
				&& cross_magnitude > 4096.0 * std::numeric_limits<double>::epsilon()
					* (magnitude_u * magnitude_u + magnitude_v * magnitude_v);
		}

		bool same_candidate(const UvCandidate& a, const UvCandidate& b,
			double u_tolerance, double v_tolerance)
		{
			return std::abs(a.u - b.u) <= u_tolerance
				&& std::abs(a.v - b.v) <= v_tolerance;
		}

		void insert_candidate(CandidateRow& row, const UvCandidate& candidate,
			double u_tolerance, double v_tolerance)
		{
			for (UvCandidate& existing : row)
				if (same_candidate(existing, candidate, u_tolerance, v_tolerance))
				{
					existing.distance = std::min(existing.distance, candidate.distance);
					return;
				}
			row.push_back(candidate);
		}

		void sort_candidates(CandidateRow& row)
		{
			std::sort(row.begin(), row.end(), [](const UvCandidate& a, const UvCandidate& b)
			{
				if (a.u != b.u) return a.u < b.u;
				if (a.v != b.v) return a.v < b.v;
				return a.distance < b.distance;
			});
		}

		bool periodic_aliases(double value, double lo, double hi, bool periodic,
			double period, double tolerance, std::vector<double>& aliases,
			std::string& error)
		{
			aliases.clear();
			if (!std::isfinite(value)) return true;
			if (!periodic)
			{
				if (value >= lo - tolerance && value <= hi + tolerance)
					aliases.push_back(value);
				return true;
			}
			const double first_shift_real = std::ceil((lo - tolerance - value) / period);
			const double last_shift_real = std::floor((hi + tolerance - value) / period);
			if (!std::isfinite(first_shift_real) || !std::isfinite(last_shift_real))
			{
				error = "could not enumerate periodic UV aliases";
				return false;
			}
			if (last_shift_real < first_shift_real) return true;
			if (last_shift_real - first_shift_real > 64.0)
			{
				error = "periodic UV bounds contain more than 65 equivalent aliases";
				return false;
			}
			const long long first_shift = static_cast<long long>(first_shift_real);
			const long long last_shift = static_cast<long long>(last_shift_real);
			for (long long shift = first_shift; shift <= last_shift; ++shift)
				aliases.push_back(value + static_cast<double>(shift) * period);
			return true;
		}

		bool append_surface_aliases(const TopoDS_Face& face, const SurfaceChart& chart,
			const SampleProjectionTolerance& tolerance,
			const gp_Pnt& world_point, double base_u, double base_v,
			bool require_on, CandidateRow& row, bool& rejected_singular,
			CandidateAudit& audit, std::string& error)
		{
			std::vector<double> u_values, v_values;
			if (!periodic_aliases(base_u, chart.u_min, chart.u_max, chart.u_periodic,
				chart.u_period, tolerance.uv_tolerance, u_values, error)
				|| !periodic_aliases(base_v, chart.v_min, chart.v_max, chart.v_periodic,
					chart.v_period, tolerance.uv_tolerance, v_values, error)) return false;
			for (double u : u_values)
				for (double v : v_values)
				{
					++audit.aliases;
					const TopAbs_State state = trim_state(face, u, v,
						tolerance.uv_tolerance);
					if (require_on ? state != TopAbs_ON
						: state != TopAbs_IN && state != TopAbs_ON) continue;
					++audit.trim_valid;
					double distance = 0.0;
					const bool residual_valid = re_evaluates_to(chart, world_point, u, v,
						tolerance.geometric_tolerance, distance);
					if (std::isfinite(distance)) audit.minimum_residual = std::min(
						audit.minimum_residual, distance);
					if (!residual_valid) continue;
					++audit.residual_valid;
					if (!regular_parameterization(chart, u, v))
					{
						rejected_singular = true;
						continue;
					}
					++audit.regular;
					insert_candidate(row, {u, v, distance}, tolerance.u_tolerance,
						tolerance.v_tolerance);
				}
			return true;
		}

		bool collect_surface_candidates(const std::vector<gp_Pnt>& world_points,
			const TopoDS_Face& face, const SurfaceChart& chart,
			const std::vector<SampleProjectionTolerance>& tolerances, bool require_on,
			CandidateTable& candidates, std::string& error)
		{
			candidates.clear();
			candidates.resize(world_points.size());
			BRepAdaptor_Surface adaptor(face, true);
			for (std::size_t sample = 0; sample < world_points.size(); ++sample)
			{
				const SampleProjectionTolerance& tolerance = tolerances[sample];
				Extrema_ExtPS projection;
				projection.Initialize(adaptor, chart.u_min, chart.u_max, chart.v_min,
					chart.v_max, tolerance.u_tolerance, tolerance.v_tolerance);
				projection.SetFlag(Extrema_ExtFlag_MIN);
				projection.SetAlgo(Extrema_ExtAlgo_Grad);
				projection.Perform(world_points[sample]);
				bool rejected_singular = false;
				CandidateAudit audit;
				for (int candidate = 1; projection.IsDone()
					&& candidate <= projection.NbExt(); ++candidate)
				{
					++audit.extrema;
					double u = 0.0, v = 0.0;
					projection.Point(candidate).Parameter(u, v);
					if (!append_surface_aliases(face, chart, tolerance,
						world_points[sample], u, v,
						require_on, candidates[sample], rejected_singular, audit, error)) return false;
				}
				sort_candidates(candidates[sample]);
				if (candidates[sample].empty())
				{
					std::ostringstream message;
					message << "sample " << sample << (rejected_singular
						? " lies only on a rank-deficient surface chart"
						: " has no trim-valid IN/ON surface projection within tolerance")
						<< " (extrema=" << audit.extrema << ", aliases=" << audit.aliases
						<< ", trim-valid=" << audit.trim_valid << ", residual-valid="
						<< audit.residual_valid << ", regular=" << audit.regular
						<< ", min-residual=" << audit.minimum_residual
						<< " mm, tolerance=" << tolerance.geometric_tolerance << " mm)";
					error = message.str();
					return false;
				}
			}
			return true;
		}

		bool transition_allowed(const SurfaceChart& chart,
			const SampleProjectionTolerance& from_tolerance,
			const SampleProjectionTolerance& to_tolerance,
			const UvCandidate& from, const UvCandidate& to)
		{
			const double u_tolerance = std::min(from_tolerance.u_tolerance,
				to_tolerance.u_tolerance);
			const double v_tolerance = std::min(from_tolerance.v_tolerance,
				to_tolerance.v_tolerance);
			if (chart.u_periodic
				&& std::abs(to.u - from.u) >= 0.5 * chart.u_period - u_tolerance)
				return false;
			if (chart.v_periodic
				&& std::abs(to.v - from.v) >= 0.5 * chart.v_period - v_tolerance)
				return false;
			return true;
		}

		bool select_consistent_branch(const CandidateTable& candidates,
			const SurfaceChart& chart,
			const std::vector<SampleProjectionTolerance>& tolerances,
			std::vector<std::array<double, 2>>& result,
			std::string& error)
		{
			result.clear();
			if (candidates.empty()) return true;
			std::vector<unsigned char> previous_paths(candidates.front().size(), 1u);
			std::vector<std::vector<std::size_t>> predecessors(candidates.size());
			predecessors.front().assign(candidates.front().size(), 0);
			for (std::size_t sample = 1; sample < candidates.size(); ++sample)
			{
				std::vector<unsigned char> next_paths(candidates[sample].size(), 0u);
				predecessors[sample].assign(candidates[sample].size(), 0);
				for (std::size_t to = 0; to < candidates[sample].size(); ++to)
					for (std::size_t from = 0; from < candidates[sample - 1].size(); ++from)
					{
						if (previous_paths[from] == 0u
							|| !transition_allowed(chart, tolerances[sample - 1],
							tolerances[sample], candidates[sample - 1][from],
								candidates[sample][to])) continue;
						if (next_paths[to] == 0u)
							predecessors[sample][to] = from;
						next_paths[to] = static_cast<unsigned char>(std::min<unsigned int>(2u,
							static_cast<unsigned int>(next_paths[to]) + previous_paths[from]));
					}
				previous_paths = std::move(next_paths);
			}
			std::size_t best = 0;
			unsigned int full_paths = 0u;
			for (std::size_t candidate = 0; candidate < previous_paths.size(); ++candidate)
			{
				if (previous_paths[candidate] == 0u) continue;
				if (full_paths == 0u) best = candidate;
				full_paths = std::min<unsigned int>(2u, full_paths + previous_paths[candidate]);
			}
			if (full_paths == 0u)
			{
				error = (chart.u_periodic || chart.v_periodic)
					? "periodic UV candidates have no continuous branch"
					: "UV candidates have no continuous branch";
				return false;
			}
			if (full_paths > 1u)
			{
				error = (chart.u_periodic || chart.v_periodic)
					? "periodic UV chart has multiple equally valid continuous branches"
					: "UV chart has multiple equally valid continuous branches";
				return false;
			}
			std::vector<std::size_t> selected(candidates.size(), 0);
			selected.back() = best;
			for (std::size_t sample = candidates.size() - 1; sample > 0; --sample)
				selected[sample - 1] = predecessors[sample][selected[sample]];
			result.reserve(candidates.size());
			for (std::size_t sample = 0; sample < candidates.size(); ++sample)
				result.push_back({candidates[sample][selected[sample]].u,
					candidates[sample][selected[sample]].v});
			return true;
		}

		enum class PcurveStatus
		{
			unavailable,
			success,
			failure
		};

		struct BoundaryPcurveAudit
		{
			std::size_t projected_parameters = 0;
			std::size_t in_range_parameters = 0;
			std::size_t trim_in = 0;
			std::size_t trim_on = 0;
			std::size_t trim_out = 0;
			std::size_t trim_unknown = 0;
			std::size_t accepted = 0;
			double minimum_curve_residual = std::numeric_limits<double>::infinity();
			double minimum_surface_residual = std::numeric_limits<double>::infinity();
			double minimum_curve_surface_residual =
				std::numeric_limits<double>::infinity();
		};

		bool is_oriented_face_edge_occurrence(const TopoDS_Face& face,
			const TopoDS_Edge& edge)
		{
			for (TopExp_Explorer occurrence(face, TopAbs_EDGE); occurrence.More();
				occurrence.Next())
				if (TopoDS::Edge(occurrence.Current()).IsEqual(edge)) return true;
			return false;
		}

		double native_edge_tolerance(const TopoDS_Edge& edge)
		{
			const double scale = std::abs(edge.Location().Transformation().ScaleFactor());
			if (!std::isfinite(scale) || !(scale > std::numeric_limits<double>::min()))
				return std::numeric_limits<double>::infinity();
			return std::max(0.0, BRep_Tool::Tolerance(edge)) * scale;
		}

		void count_trim_state(BoundaryPcurveAudit& audit, TopAbs_State state)
		{
			switch (state)
			{
			case TopAbs_IN: ++audit.trim_in; break;
			case TopAbs_ON: ++audit.trim_on; break;
			case TopAbs_OUT: ++audit.trim_out; break;
			default: ++audit.trim_unknown; break;
			}
		}

		std::string boundary_audit_message(std::size_t sample,
			const BoundaryPcurveAudit& audit, double sample_tolerance,
			double curve_surface_tolerance)
		{
			std::ostringstream message;
			message << "exact boundary pcurve rejected sample " << sample
				<< " (parameters=" << audit.projected_parameters
				<< ", in-range=" << audit.in_range_parameters
				<< ", curve-residual-min=" << audit.minimum_curve_residual << " mm"
				<< ", trim-state IN/ON/OUT/UNKNOWN=" << audit.trim_in << "/"
				<< audit.trim_on << "/" << audit.trim_out << "/"
				<< audit.trim_unknown
				<< ", surface-residual-min=" << audit.minimum_surface_residual << " mm"
				<< ", curve-surface-residual-min="
				<< audit.minimum_curve_surface_residual << " mm"
				<< ", sample-limit=" << sample_tolerance << " mm"
				<< ", curve-surface-limit=" << curve_surface_tolerance << " mm)";
			return message.str();
		}

		PcurveStatus project_exact_pcurve(const std::vector<gp_Pnt>& world_points,
			const TopoDS_Face& face, const TopoDS_Edge& edge, const SurfaceChart& chart,
			const std::vector<double>& source_geometric_tolerances_mm,
			std::vector<std::array<double, 2>>& result, std::string& error)
		{
			// CurveOnSurface is authoritative only for the exact oriented occurrence
			// supplied by the face wire. IsSame() would lose orientation and could select
			// the wrong pcurve branch of a seam, so require IsEqual() membership first.
			if (!is_oriented_face_edge_occurrence(face, edge))
			{
				error = "exact boundary edge is not an oriented occurrence of the face";
				return PcurveStatus::failure;
			}
			double pcurve_first = 0.0, pcurve_last = 0.0;
			const Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(
				edge, face, pcurve_first, pcurve_last);
			if (pcurve.IsNull()) return PcurveStatus::unavailable;
			TopLoc_Location curve_location;
			double curve_first = 0.0, curve_last = 0.0;
			const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, curve_location,
				curve_first, curve_last);
			if (curve.IsNull() || !finite_bounds(curve_first, curve_last)
				|| !finite_bounds(pcurve_first, pcurve_last))
			{
				error = "exact boundary edge has an invalid 3D curve or pcurve range";
				return PcurveStatus::failure;
			}
			if (!BRep_Tool::SameParameter(edge) || !BRep_Tool::SameRange(edge))
			{
				error = "exact boundary edge is not SameParameter/SameRange";
				return PcurveStatus::failure;
			}
			const double range_tolerance = std::max(Precision::PConfusion(),
				128.0 * std::numeric_limits<double>::epsilon()
					* std::max({1.0, std::abs(curve_first), std::abs(curve_last),
						std::abs(pcurve_first), std::abs(pcurve_last)}));
			if (std::abs(curve_first - pcurve_first) > range_tolerance
				|| std::abs(curve_last - pcurve_last) > range_tolerance)
			{
				error = "exact boundary edge 3D/pcurve parameter ranges differ";
				return PcurveStatus::failure;
			}
			const double first = std::max(curve_first, pcurve_first);
			const double last = std::min(curve_last, pcurve_last);
			if (!(last >= first))
			{
				error = "exact boundary edge 3D/pcurve parameter ranges do not overlap";
				return PcurveStatus::failure;
			}
			const double parameter_tolerance = std::max(Precision::PConfusion(),
				128.0 * std::numeric_limits<double>::epsilon()
					* std::max({1.0, std::abs(first), std::abs(last)}));
			const double edge_tolerance = native_edge_tolerance(edge);
			if (!std::isfinite(edge_tolerance))
			{
				error = "exact boundary edge has an invalid native tolerance";
				return PcurveStatus::failure;
			}
			std::vector<SampleProjectionTolerance> tolerances;
			if (!initialize_sample_tolerances(face, chart,
				source_geometric_tolerances_mm, edge_tolerance, tolerances, error))
				return PcurveStatus::failure;
			// The public chain was evaluated on the certified source curve. A distinct
			// target trim occurrence may differ from it by the sum of their native CAD
			// bounds. The caller-supplied bound is the source-chain side of that exact
			// certificate; the target edge and face bounds are read here, where their
			// actual OCCT occurrences are available.
			const double curve_surface_tolerance = std::max(Precision::Confusion(),
				edge_tolerance + chart.face_native_tolerance);
			CandidateTable candidates(world_points.size());
			for (std::size_t sample = 0; sample < world_points.size(); ++sample)
			{
				const SampleProjectionTolerance& tolerance = tolerances[sample];
				const double sample_tolerance = tolerance.geometric_tolerance;
				BoundaryPcurveAudit audit;
				gp_Pnt local_point = world_points[sample];
				local_point.Transform(curve_location.Transformation().Inverted());
				GeomAPI_ProjectPointOnCurve projection(local_point, curve, first, last);
				std::vector<double> parameters;
				for (int candidate = 1; candidate <= projection.NbPoints(); ++candidate)
					parameters.push_back(projection.Parameter(candidate));
				gp_Pnt first_world = curve->Value(first);
				first_world.Transform(curve_location.Transformation());
				if (first_world.Distance(world_points[sample]) <= sample_tolerance)
					parameters.push_back(first);
				gp_Pnt last_world = curve->Value(last);
				last_world.Transform(curve_location.Transformation());
				if (last_world.Distance(world_points[sample]) <= sample_tolerance)
					parameters.push_back(last);
				std::sort(parameters.begin(), parameters.end());
				parameters.erase(std::unique(parameters.begin(), parameters.end(),
					[&](double a, double b) { return std::abs(a - b) <= parameter_tolerance; }),
					parameters.end());
				audit.projected_parameters = parameters.size();
				for (double parameter : parameters)
				{
					if (!std::isfinite(parameter) || parameter < first - parameter_tolerance
						|| parameter > last + parameter_tolerance) continue;
					++audit.in_range_parameters;
					gp_Pnt curve_world = curve->Value(parameter);
					curve_world.Transform(curve_location.Transformation());
					const double curve_distance = curve_world.Distance(world_points[sample]);
					if (std::isfinite(curve_distance))
						audit.minimum_curve_residual = std::min(
							audit.minimum_curve_residual, curve_distance);
					const gp_Pnt2d uv = pcurve->Value(parameter);
					if (!std::isfinite(uv.X()) || !std::isfinite(uv.Y())) continue;
					count_trim_state(audit,
						trim_state(face, uv.X(), uv.Y(), tolerance.uv_tolerance));
					gp_Pnt surface_world = chart.surface->Value(uv.X(), uv.Y());
					surface_world.Transform(chart.location.Transformation());
					const double surface_distance =
						surface_world.Distance(world_points[sample]);
					const double curve_surface_distance =
						surface_world.Distance(curve_world);
					if (std::isfinite(surface_distance))
						audit.minimum_surface_residual = std::min(
							audit.minimum_surface_residual, surface_distance);
					if (std::isfinite(curve_surface_distance))
						audit.minimum_curve_surface_residual = std::min(
							audit.minimum_curve_surface_residual, curve_surface_distance);
					if (!std::isfinite(curve_distance) || !std::isfinite(surface_distance)
						|| !std::isfinite(curve_surface_distance)
						|| curve_distance > sample_tolerance
						|| surface_distance > sample_tolerance
						|| curve_surface_distance > curve_surface_tolerance) continue;
					// The pcurve is itself the face wire's oriented boundary definition.
					// Requiring the generic classifier to rediscover TopAbs_ON is redundant
					// and rejects valid endpoint/tolerance-near values. Membership, range,
					// SameParameter/SameRange and both independent 3D checks above are the
					// stronger proof; the classifier state remains in failure diagnostics.
					++audit.accepted;
					insert_candidate(candidates[sample], {uv.X(), uv.Y(),
						std::max({curve_distance, surface_distance,
							curve_surface_distance})}, tolerance.u_tolerance,
							tolerance.v_tolerance);
				}
				sort_candidates(candidates[sample]);
				if (candidates[sample].empty())
				{
					error = boundary_audit_message(sample, audit, sample_tolerance,
						curve_surface_tolerance);
					return PcurveStatus::failure;
				}
			}
			if (!select_consistent_branch(candidates, chart, tolerances, result, error))
			{
				error = "exact boundary pcurve: " + error;
				return PcurveStatus::failure;
			}
			return PcurveStatus::success;
		}

		bool same_chart(const std::vector<std::array<double, 2>>& a,
			const std::vector<std::array<double, 2>>& b,
			const std::vector<SampleProjectionTolerance>& tolerances)
		{
			if (a.size() != b.size() || a.size() != tolerances.size()) return false;
			for (std::size_t sample = 0; sample < a.size(); ++sample)
				if (std::abs(a[sample][0] - b[sample][0])
						> tolerances[sample].u_tolerance
					|| std::abs(a[sample][1] - b[sample][1])
						> tolerances[sample].v_tolerance) return false;
			return true;
		}
	}

	bool project_trim_valid_face_uv(const std::vector<gp_Pnt>& world_points,
		const TopoDS_Face& face, const std::vector<TopoDS_Edge>& exact_boundary_edges,
		bool boundary_use,
		const std::vector<double>& sample_geometric_tolerances_mm,
		std::vector<std::array<double, 2>>& sample_uv, std::string& error,
		bool* used_exact_boundary_pcurve)
	{
		sample_uv.clear();
		error.clear();
		if (used_exact_boundary_pcurve) *used_exact_boundary_pcurve = false;
		if (face.IsNull())
		{
			error = "cannot project onto a null face";
			return false;
		}
		if (sample_geometric_tolerances_mm.size() != world_points.size())
		{
			error = "sample geometric tolerance count does not match point count";
			return false;
		}
		SurfaceChart chart;
		if (!initialize_chart(face, chart, error)) return false;
		std::vector<SampleProjectionTolerance> tolerances;
		if (!initialize_sample_tolerances(face, chart,
			sample_geometric_tolerances_mm, 0.0, tolerances, error)) return false;
		if (world_points.empty()) return true;

		if (boundary_use)
		{
			bool found_pcurve = false;
			std::vector<std::array<double, 2>> exact_result;
			for (std::size_t edge = 0; edge < exact_boundary_edges.size(); ++edge)
			{
				std::vector<std::array<double, 2>> candidate_result;
				std::string candidate_error;
				const PcurveStatus status = project_exact_pcurve(world_points, face,
					exact_boundary_edges[edge], chart,
					sample_geometric_tolerances_mm, candidate_result, candidate_error);
				if (status == PcurveStatus::unavailable) continue;
				found_pcurve = true;
				if (status == PcurveStatus::failure)
				{
					error = "exact boundary edge " + std::to_string(edge)
						+ " failed: " + candidate_error;
					return false;
				}
				if (exact_result.empty()) exact_result = std::move(candidate_result);
				else if (!same_chart(exact_result, candidate_result, tolerances))
				{
					error = "exact boundary pcurves select different UV branches";
					return false;
				}
			}
			if (found_pcurve)
			{
				sample_uv = std::move(exact_result);
				if (used_exact_boundary_pcurve) *used_exact_boundary_pcurve = true;
				return true;
			}
		}

		CandidateTable candidates;
		if (!collect_surface_candidates(world_points, face, chart, tolerances,
			boundary_use,
			candidates, error))
			return false;
		return select_consistent_branch(candidates, chart, tolerances, sample_uv, error);
	}

	bool project_trim_valid_face_uv(const std::vector<gp_Pnt>& world_points,
		const TopoDS_Face& face, const std::vector<TopoDS_Edge>& exact_boundary_edges,
		bool boundary_use, double geometric_tolerance_mm,
		std::vector<std::array<double, 2>>& sample_uv, std::string& error,
		bool* used_exact_boundary_pcurve)
	{
		const double tolerance = std::max(0.0, geometric_tolerance_mm);
		return project_trim_valid_face_uv(world_points, face, exact_boundary_edges,
			boundary_use, std::vector<double>(world_points.size(), tolerance), sample_uv,
			error, used_exact_boundary_pcurve);
	}
}
