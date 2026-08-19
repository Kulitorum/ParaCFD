// Exact transactional normalization of positive-area coplanar source-face overlaps.
//
// This OCCT-facing component is deliberately separate from STEP import.  It implements a
// fail-closed transaction suitable for zero-thickness fabric whose source model contains
// coincident panel coverage (for example two cell-local copies of a shared rib):
//
//   1. detect planar overlap from exact OCCT face/face common parts, never tessellation
//      proximity, while exposing generic curved candidates as an incomplete preflight;
//   2. form connected coplanar overlap clusters;
//   3. build each cluster's explicit Boolean membership cells as non-overlapping faces;
//   4. retain every original source-face contributor and its orientation parity; and
//   5. prove that the atomic coverage neither adds nor removes fabric and is counted once.
//
// There is no sewing, fuzzy tolerance, gap healing, same-domain unification, tessellation,
// or STEP-path mutation here.  A failed invariant rejects the complete transaction.

#include "core/geometry/occt_surface_normalization.h"

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GeomConvert_BSplineSurfaceToBezierSurface.hxx>
#include <GeomConvert_SurfToAnaSurf.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Surface.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace paracfd::core
{
	namespace
	{
	using Clock = std::chrono::steady_clock;
	using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

	std::vector<TopoDS_Face> unique_faces(const TopoDS_Shape& shape)
	{
		ShapeMap faces;
		TopExp::MapShapes(shape, TopAbs_FACE, faces);
		std::vector<TopoDS_Face> result;
		result.reserve(static_cast<std::size_t>(faces.Extent()));
		for (Standard_Integer face = 1; face <= faces.Extent(); ++face)
			result.push_back(TopoDS::Face(faces(face)));
		return result;
	}

	TopoDS_Shape face_compound(const std::vector<TopoDS_Face>& faces)
	{
		TopoDS_Compound result;
		BRep_Builder builder;
		builder.MakeCompound(result);
		for (const TopoDS_Face& face : faces) builder.Add(result, face);
		return result;
	}

	struct Vec3
	{
		double x = 0.0, y = 0.0, z = 0.0;
	};

	Vec3 operator+(Vec3 a, Vec3 b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
	Vec3 operator-(Vec3 a, Vec3 b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
	Vec3 operator*(Vec3 a, double scale) { return { a.x * scale, a.y * scale, a.z * scale }; }
	double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
	double length(Vec3 vector) { return std::sqrt(std::max(0.0, dot(vector, vector))); }
	Vec3 normalized(Vec3 vector)
	{
		const double magnitude = length(vector);
		return magnitude > 0.0 ? vector * (1.0 / magnitude) : Vec3{};
	}
	Vec3 as_vec(const gp_Dir& direction)
	{
		return { direction.X(), direction.Y(), direction.Z() };
	}
	Vec3 as_vec(const gp_Pnt& point) { return { point.X(), point.Y(), point.Z() }; }

	Vec3 canonical_direction(Vec3 direction)
	{
		direction = normalized(direction);
		const double component[3] = { direction.x, direction.y, direction.z };
		for (double value : component)
			if (std::abs(value) > 64.0 * std::numeric_limits<double>::epsilon())
				return value < 0.0 ? direction * -1.0 : direction;
		return direction;
	}

	struct SurfaceMeasure
	{
		double area = 0.0; // mm^2
		Vec3 centroid{};   // mm
	};

	SurfaceMeasure surface_measure(const TopoDS_Shape& shape)
	{
		SurfaceMeasure result;
		const std::vector<TopoDS_Face> faces = unique_faces(shape);
		long double area = 0.0L, x = 0.0L, y = 0.0L, z = 0.0L;
		for (const TopoDS_Face& face : faces)
		{
			GProp_GProps properties;
			// The no-Eps overload uses a non-adaptive integration path. On Rib 1's
			// multiply-trimmed planar face that path overstates area by 164.613 mm^2,
			// falsely suggesting a source-only sliver which exact Common/Cut cannot find.
			// Conservation must use the adaptive overload with an explicit relative error.
			BRepGProp::SurfaceProperties(face, properties, 1.0e-12, false);
			const double face_area = std::abs(properties.Mass());
			if (!(face_area > 0.0) || !std::isfinite(face_area)) continue;
			const gp_Pnt centre = properties.CentreOfMass();
			area += static_cast<long double>(face_area);
			x += static_cast<long double>(face_area) * centre.X();
			y += static_cast<long double>(face_area) * centre.Y();
			z += static_cast<long double>(face_area) * centre.Z();
		}
		result.area = static_cast<double>(area);
		if (area > 0.0L)
			result.centroid = { static_cast<double>(x / area), static_cast<double>(y / area),
				static_cast<double>(z / area) };
		return result;
	}

	double boundary_length(const TopoDS_Face& face)
	{
		GProp_GProps properties;
		BRepGProp::LinearProperties(face, properties);
		return std::abs(properties.Mass());
	}

	double maximum_tolerance(const TopoDS_Shape& shape)
	{
		double result = 0.0;
		ShapeMap faces, edges, vertices;
		TopExp::MapShapes(shape, TopAbs_FACE, faces);
		TopExp::MapShapes(shape, TopAbs_EDGE, edges);
		TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
		for (Standard_Integer index = 1; index <= faces.Extent(); ++index)
			result = std::max(result, BRep_Tool::Tolerance(TopoDS::Face(faces(index))));
		for (Standard_Integer index = 1; index <= edges.Extent(); ++index)
			result = std::max(result, BRep_Tool::Tolerance(TopoDS::Edge(edges(index))));
		for (Standard_Integer index = 1; index <= vertices.Extent(); ++index)
			result = std::max(result, BRep_Tool::Tolerance(TopoDS::Vertex(vertices(index))));
		return result;
	}

	double numeric_area_tolerance(double area, double perimeter, double native_tolerance,
		double length_scale, std::size_t operation_count = 1)
	{
		// All terms retain their physical dimensions.  A one-model-unit floor makes
		// this decision depend on the STEP author's unit scale and can hide an entire
		// small duplicate face.  Native BRep tolerance remains authoritative; the
		// floating-point term is relative to the actual local extent.
		const double scale = std::max(length_scale,
			std::sqrt(std::max(0.0, area)));
		if (!(scale > 0.0) || !std::isfinite(scale)) return 0.0;
		const double tolerance = std::max(std::max(0.0, native_tolerance),
			64.0 * std::numeric_limits<double>::epsilon() * scale);
		const double boundary_uncertainty = 64.0 * tolerance * std::max(perimeter, scale);
		const double roundoff = 4096.0 * std::numeric_limits<double>::epsilon()
			* std::max(area, scale * scale);
		return static_cast<double>(std::max<std::size_t>(1, operation_count))
			* std::max(boundary_uncertainty, roundoff);
	}

	struct Interval
	{
		double lower = std::numeric_limits<double>::infinity();
		double upper = -std::numeric_limits<double>::infinity();
	};

	Interval interval_value(double value)
	{
		return { value, value };
	}

	Interval outward_interval(double lower, double upper)
	{
		return { std::nextafter(lower, -std::numeric_limits<double>::infinity()),
			std::nextafter(upper, std::numeric_limits<double>::infinity()) };
	}

	Interval interval_add(Interval a, Interval b)
	{
		return outward_interval(a.lower + b.lower, a.upper + b.upper);
	}

	Interval interval_subtract(Interval a, Interval b)
	{
		return outward_interval(a.lower - b.upper, a.upper - b.lower);
	}

	Interval interval_multiply(Interval a, Interval b)
	{
		const double products[4] = { a.lower * b.lower, a.lower * b.upper,
			a.upper * b.lower, a.upper * b.upper };
		return outward_interval(*std::min_element(std::begin(products), std::end(products)),
			*std::max_element(std::begin(products), std::end(products)));
	}

	Interval interval_hull(const std::vector<Interval>& values)
	{
		Interval result;
		for (Interval value : values)
		{
			result.lower = std::min(result.lower, value.lower);
			result.upper = std::max(result.upper, value.upper);
		}
		return result;
	}

	std::array<Interval, 3> interval_cross(const std::array<Interval, 3>& a,
		const std::array<Interval, 3>& b)
	{
		return {
			interval_subtract(interval_multiply(a[1], b[2]), interval_multiply(a[2], b[1])),
			interval_subtract(interval_multiply(a[2], b[0]), interval_multiply(a[0], b[2])),
			interval_subtract(interval_multiply(a[0], b[1]), interval_multiply(a[1], b[0]))
		};
	}

	struct BezierPatchEnvelope
	{
		std::array<double, 3> lower{};
		std::array<double, 3> upper{};
		std::vector<gp_Pnt> poles;
		std::array<Interval, 3> normal;
		bool normal_certified = false;
	};

	struct BezierSupportEnvelope
	{
		bool certified = false;
		std::vector<BezierPatchEnvelope> patches;
	};

	bool append_bezier_envelope(const occ::handle<Geom_BezierSurface>& patch,
		double native_tolerance, std::vector<BezierPatchEnvelope>& destination)
	{
		if (patch.IsNull()) return false;
		BezierPatchEnvelope envelope;
		envelope.lower.fill(std::numeric_limits<double>::max());
		envelope.upper.fill(-std::numeric_limits<double>::max());
		double coordinate_scale = 1.0;
		for (int u = 1; u <= patch->NbUPoles(); ++u)
			for (int v = 1; v <= patch->NbVPoles(); ++v)
			{
				const double weight = patch->Weight(u, v);
				const gp_Pnt pole = patch->Pole(u, v);
				if (!(weight > gp::Resolution()) || !std::isfinite(weight)
					|| !std::isfinite(pole.X()) || !std::isfinite(pole.Y())
					|| !std::isfinite(pole.Z())) return false;
				envelope.poles.push_back(pole);
				const double coordinate[3] = { pole.X(), pole.Y(), pole.Z() };
				for (int axis = 0; axis < 3; ++axis)
				{
					envelope.lower[axis] = std::min(envelope.lower[axis], coordinate[axis]);
					envelope.upper[axis] = std::max(envelope.upper[axis], coordinate[axis]);
					coordinate_scale = std::max(coordinate_scale, std::abs(coordinate[axis]));
				}
			}
		const double padding = std::max({ native_tolerance, Precision::Confusion(),
			128.0 * std::numeric_limits<double>::epsilon() * coordinate_scale });
		for (int axis = 0; axis < 3; ++axis)
		{
			envelope.lower[axis] -= padding;
			envelope.upper[axis] += padding;
		}

		// Bound the two rational derivatives in homogeneous coordinates.  Bezier
		// basis functions form a non-negative partition of unity, so every surface
		// value and derivative lies in the interval hull of its (derivative) control
		// values.  The common positive W^4 denominator does not affect whether two
		// normals can be parallel.
		const int u_count = patch->NbUPoles();
		const int v_count = patch->NbVPoles();
		const gp_Pnt reference = patch->Pole(1, 1);
		std::vector<std::array<Interval, 4>> homogeneous(
			static_cast<std::size_t>(u_count * v_count));
		auto control = [&](int u, int v) -> std::array<Interval, 4>&
		{
			return homogeneous[static_cast<std::size_t>((u - 1) * v_count + (v - 1))];
		};
		for (int u = 1; u <= u_count; ++u)
			for (int v = 1; v <= v_count; ++v)
			{
				const gp_Pnt pole = patch->Pole(u, v);
				const Interval weight = interval_value(patch->Weight(u, v));
				control(u, v) = {
					interval_multiply(weight, interval_subtract(interval_value(pole.X()),
						interval_value(reference.X()))),
					interval_multiply(weight, interval_subtract(interval_value(pole.Y()),
						interval_value(reference.Y()))),
					interval_multiply(weight, interval_subtract(interval_value(pole.Z()),
						interval_value(reference.Z()))),
					weight
				};
			}
		std::array<Interval, 4> value_hull;
		std::array<Interval, 4> derivative_u;
		std::array<Interval, 4> derivative_v;
		for (int component = 0; component < 4; ++component)
		{
			std::vector<Interval> values, u_values, v_values;
			values.reserve(homogeneous.size());
			u_values.reserve(static_cast<std::size_t>(std::max(0, u_count - 1) * v_count));
			v_values.reserve(static_cast<std::size_t>(u_count * std::max(0, v_count - 1)));
			for (int u = 1; u <= u_count; ++u)
				for (int v = 1; v <= v_count; ++v)
				{
					values.push_back(control(u, v)[component]);
					if (u < u_count)
						u_values.push_back(interval_multiply(interval_value(patch->UDegree()),
							interval_subtract(control(u + 1, v)[component],
								control(u, v)[component])));
					if (v < v_count)
						v_values.push_back(interval_multiply(interval_value(patch->VDegree()),
							interval_subtract(control(u, v + 1)[component],
								control(u, v)[component])));
				}
			value_hull[component] = interval_hull(values);
			derivative_u[component] = interval_hull(u_values);
			derivative_v[component] = interval_hull(v_values);
		}
		if (value_hull[3].lower > 0.0 && u_count > 1 && v_count > 1)
		{
			std::array<Interval, 3> euclidean_u, euclidean_v;
			for (int component = 0; component < 3; ++component)
			{
				euclidean_u[component] = interval_subtract(
					interval_multiply(derivative_u[component], value_hull[3]),
					interval_multiply(value_hull[component], derivative_u[3]));
				euclidean_v[component] = interval_subtract(
					interval_multiply(derivative_v[component], value_hull[3]),
					interval_multiply(value_hull[component], derivative_v[3]));
			}
			envelope.normal = interval_cross(euclidean_u, euclidean_v);
			envelope.normal_certified = std::all_of(envelope.normal.begin(),
				envelope.normal.end(), [](Interval value)
				{ return std::isfinite(value.lower) && std::isfinite(value.upper); });
		}
		destination.push_back(std::move(envelope));
		return true;
	}

	bool subdivide_bezier_envelopes(const occ::handle<Geom_BezierSurface>& patch,
		int depth, double native_tolerance, std::vector<BezierPatchEnvelope>& destination)
	{
		if (patch.IsNull()) return false;
		if (depth <= 0) return append_bezier_envelope(patch, native_tolerance, destination);
		double u_first = 0.0, u_last = 0.0, v_first = 0.0, v_last = 0.0;
		patch->Bounds(u_first, u_last, v_first, v_last);
		const double u_middle = 0.5 * (u_first + u_last);
		const double v_middle = 0.5 * (v_first + v_last);
		const double u_bounds[3] = { u_first, u_middle, u_last };
		const double v_bounds[3] = { v_first, v_middle, v_last };
		for (int u = 0; u < 2; ++u)
			for (int v = 0; v < 2; ++v)
			{
				occ::handle<Geom_BezierSurface> child =
					occ::down_cast<Geom_BezierSurface>(patch->Copy());
				if (child.IsNull()) return false;
				child->Segment(u_bounds[u], u_bounds[u + 1], v_bounds[v], v_bounds[v + 1]);
				if (!subdivide_bezier_envelopes(child, depth - 1, native_tolerance,
					destination)) return false;
			}
		return true;
	}

	BezierSupportEnvelope make_bezier_support_envelope(const BRepAdaptor_Surface& surface,
		double native_tolerance)
	{
		BezierSupportEnvelope result;
		try
		{
			constexpr int subdivision_depth = 1;
			if (surface.GetType() == GeomAbs_BSplineSurface)
			{
				const occ::handle<Geom_BSplineSurface> spline = surface.BSpline();
				if (spline.IsNull() || surface.IsUPeriodic() || surface.IsVPeriodic()) return result;
				GeomConvert_BSplineSurfaceToBezierSurface converter(spline,
					surface.FirstUParameter(), surface.LastUParameter(),
					surface.FirstVParameter(), surface.LastVParameter(),
					Precision::PConfusion());
				for (int u = 1; u <= converter.NbUPatches(); ++u)
					for (int v = 1; v <= converter.NbVPatches(); ++v)
						if (!subdivide_bezier_envelopes(converter.Patch(u, v),
							subdivision_depth, native_tolerance, result.patches))
							return {};
			}
			else if (surface.GetType() == GeomAbs_BezierSurface)
			{
				occ::handle<Geom_BezierSurface> patch =
					occ::down_cast<Geom_BezierSurface>(surface.Bezier()->Copy());
				if (patch.IsNull()) return result;
				patch->Segment(surface.FirstUParameter(), surface.LastUParameter(),
					surface.FirstVParameter(), surface.LastVParameter());
				if (!subdivide_bezier_envelopes(patch, subdivision_depth,
					native_tolerance, result.patches)) return {};
			}
			else return result;
			result.certified = !result.patches.empty();
		}
		catch (const Standard_Failure&)
		{
			return {};
		}
		catch (const std::exception&)
		{
			return {};
		}
		return result;
	}

	struct PlanarFace
	{
		bool planar = false;
		gp_Pln plane{};
		Vec3 surface_normal{};
		Vec3 oriented_normal{};
		SurfaceMeasure measure{};
		double perimeter = 0.0;
		double tolerance = 0.0;
		Bnd_Box bounds{};
		BezierSupportEnvelope support_envelope;
	};

	PlanarFace inspect_face(const TopoDS_Face& face)
	{
		PlanarFace result;
		result.measure = surface_measure(face);
		result.perimeter = boundary_length(face);
		result.tolerance = maximum_tolerance(face);
		BRepBndLib::AddOptimal(face, result.bounds, false, true);
		BRepAdaptor_Surface surface(face, true);
		result.support_envelope = make_bezier_support_envelope(surface, result.tolerance);
		if (surface.GetType() != GeomAbs_Plane) return result;
		result.plane = surface.Plane();
		result.planar = true;
		result.surface_normal = normalized(as_vec(result.plane.Axis().Direction()));
		result.oriented_normal = result.surface_normal;
		if (face.Orientation() == TopAbs_REVERSED)
			result.oriented_normal = result.oriented_normal * -1.0;
		else if (face.Orientation() != TopAbs_FORWARD)
			result.planar = false;
		return result;
	}

	TopoDS_Wire copied_wire(const TopoDS_Wire& wire, bool reverse = false)
	{
		BRepBuilderAPI_Copy copy(wire, true, false);
		if (!copy.IsDone() || copy.Shape().IsNull() || copy.Shape().ShapeType() != TopAbs_WIRE)
			throw std::runtime_error("failed to copy planar trim wire");
		TopoDS_Wire result = TopoDS::Wire(copy.Shape());
		if (reverse) result.Reverse();
		return result;
	}

	TopoDS_Face rehost_planar_face(const TopoDS_Face& source, const gp_Pln& plane,
		bool reverse_trim_wires)
	{
		const TopoDS_Wire outer = BRepTools::OuterWire(source);
		if (outer.IsNull()) throw std::runtime_error("planar source face has no outer wire");
		BRepBuilderAPI_MakeFace make(plane, copied_wire(outer, reverse_trim_wires), true);
		if (!make.IsDone()) throw std::runtime_error("failed to rebuild outer wire on canonical plane");
		for (TopExp_Explorer wire(source, TopAbs_WIRE); wire.More(); wire.Next())
			if (!wire.Current().IsSame(outer)) make.Add(copied_wire(TopoDS::Wire(wire.Current()),
				reverse_trim_wires));
		if (!make.IsDone()) throw std::runtime_error("failed to add inner wire on canonical plane");
		return make.Face();
	}

	bool same_plane(const PlanarFace& a, const PlanarFace& b)
	{
		if (!a.planar || !b.planar) return false;
		const Vec3 cross{ a.surface_normal.y * b.surface_normal.z
			- a.surface_normal.z * b.surface_normal.y,
			a.surface_normal.z * b.surface_normal.x - a.surface_normal.x * b.surface_normal.z,
			a.surface_normal.x * b.surface_normal.y - a.surface_normal.y * b.surface_normal.x };
		if (length(cross) > std::max(Precision::Angular(),
			64.0 * std::numeric_limits<double>::epsilon())) return false;
		const double tolerance = std::max({ Precision::Confusion(), a.tolerance, b.tolerance });
		return a.plane.Distance(b.plane.Location()) <= tolerance
			&& b.plane.Distance(a.plane.Location()) <= tolerance;
	}

	struct BooleanCommon
	{
		bool completed = false;
		bool warnings = false;
		SurfaceMeasure measure{};
		std::string report;
	};

	struct BooleanParts
	{
		bool completed = false;
		bool warnings = false;
		TopoDS_Shape shape;
		std::vector<TopoDS_Face> faces;
		SurfaceMeasure measure{};
		std::string report;
	};

	template<class Operation>
	BooleanParts exact_boolean_parts(const TopoDS_Face& argument, const TopoDS_Face& tool)
	{
		BooleanParts result;
		try
		{
			NCollection_List<TopoDS_Shape> arguments, tools;
			arguments.Append(argument);
			tools.Append(tool);
			Operation operation;
			operation.SetArguments(arguments);
			operation.SetTools(tools);
			operation.SetFuzzyValue(0.0);
			operation.SetNonDestructive(true);
			operation.SetRunParallel(false);
			operation.SetUseOBB(true);
			operation.SetToFillHistory(true);
			operation.Build();
			// Boolean warnings mean the returned topology is not a proof.  This helper
			// is a normalization correctness gate, so warning-degraded results fail the
			// transaction exactly like errors instead of being used to infer "no overlap".
			result.completed = operation.IsDone() && !operation.HasErrors()
				&& !operation.HasWarnings();
			result.warnings = operation.HasWarnings();
			if (operation.HasErrors() || operation.HasWarnings())
			{
				std::ostringstream report;
				if (operation.HasErrors()) operation.DumpErrors(report);
				if (operation.HasWarnings()) operation.DumpWarnings(report);
				result.report = report.str();
			}
			if (result.completed && !operation.Shape().IsNull())
			{
				result.shape = operation.Shape();
				result.faces = unique_faces(result.shape);
				result.measure = surface_measure(result.shape);
			}
		}
		catch (const Standard_Failure& failure)
		{
			result.report = std::string("OpenCascade exception: ") + failure.GetMessageString();
		}
		catch (const std::exception& exception)
		{
			result.report = std::string("exception: ") + exception.what();
		}
		return result;
	}

	BooleanCommon exact_common(const TopoDS_Face& a, const TopoDS_Face& b)
	{
		const BooleanParts parts = exact_boolean_parts<BRepAlgoAPI_Common>(a, b);
		BooleanCommon result;
		result.completed = parts.completed;
		result.warnings = parts.warnings;
		result.measure = parts.measure;
		result.report = parts.report;
		return result;
	}

	struct DisjointSet
	{
		explicit DisjointSet(std::size_t size) : parent(size), rank(size, 0)
		{
			std::iota(parent.begin(), parent.end(), std::size_t{ 0 });
		}
		std::size_t find(std::size_t value)
		{
			std::size_t root = value;
			while (parent[root] != root) root = parent[root];
			while (parent[value] != value)
			{
				const std::size_t next = parent[value];
				parent[value] = root;
				value = next;
			}
			return root;
		}
		void join(std::size_t a, std::size_t b)
		{
			a = find(a); b = find(b); if (a == b) return;
			if (rank[a] < rank[b]) std::swap(a, b);
			parent[b] = a;
			if (rank[a] == rank[b]) ++rank[a];
		}
		std::vector<std::size_t> parent;
		std::vector<unsigned char> rank;
	};

	struct OverlapPair
	{
		std::size_t a = 0, b = 0;
		SurfaceMeasure common{};
		double area_tolerance = 0.0;
		std::string warnings;
	};

	struct UnsupportedOverlap
	{
		std::size_t a = 0, b = 0;
		SurfaceMeasure common{};
		double area_tolerance = 0.0;
	};

	struct SupportScreen
	{
		bool possible_same_domain = true;
		bool exact_confirmation_required = false;
		std::string warning;
	};

	bool envelopes_intersect(const BezierPatchEnvelope& a, const BezierPatchEnvelope& b)
	{
		for (int axis = 0; axis < 3; ++axis)
			if (a.upper[axis] < b.lower[axis] || b.upper[axis] < a.lower[axis]) return false;
		return true;
	}

	bool normals_may_be_parallel(const BezierPatchEnvelope& a,
		const BezierPatchEnvelope& b)
	{
		if (!a.normal_certified || !b.normal_certified) return true;
		const std::array<Interval, 3> cross = interval_cross(a.normal, b.normal);
		return std::all_of(cross.begin(), cross.end(), [](Interval component)
		{
			return component.lower <= 0.0 && component.upper >= 0.0;
		});
	}

	bool envelope_intersects_bounds(const BezierPatchEnvelope& envelope,
		const Bnd_Box& bounds)
	{
		if (bounds.IsVoid()) return true; // indeterminate is retained
		double x_min = 0.0, y_min = 0.0, z_min = 0.0;
		double x_max = 0.0, y_max = 0.0, z_max = 0.0;
		bounds.Get(x_min, y_min, z_min, x_max, y_max, z_max);
		const double lower[3] = { x_min, y_min, z_min };
		const double upper[3] = { x_max, y_max, z_max };
		for (int axis = 0; axis < 3; ++axis)
			if (envelope.upper[axis] < lower[axis] || upper[axis] < envelope.lower[axis])
				return false;
		return true;
	}

	bool patch_can_lie_on_plane(const BezierPatchEnvelope& patch, const gp_Pln& plane,
		double tolerance)
	{
		// A positive-weight rational Bezier patch lies in a plane on an open set only
		// if its full analytic patch lies in that plane.  Bernstein independence then
		// requires every weighted control pole to lie in the same plane slab.
		return std::all_of(patch.poles.begin(), patch.poles.end(),
			[&](const gp_Pnt& pole) { return plane.Distance(pole) <= tolerance; });
	}

	// Cheap, conservative support-surface screen for every non-coplanar AABB candidate
	// involving a non-planar face.  Negative results come only from exact Bezier convex-
	// hull/normal separation or the all-control-poles planar-patch certificate.  Identical
	// handles and equal elementary supports request exact confirmation; generic indeterminate
	// representations remain explicitly unverified by this planar-only helper.
	SupportScreen plausible_same_surface(const TopoDS_Face& a, const TopoDS_Face& b,
		const PlanarFace& info_a, const PlanarFace& info_b, double tolerance)
	{
		SupportScreen result;
		try
		{
			TopLoc_Location location_a, location_b;
			const occ::handle<Geom_Surface>& local_a = BRep_Tool::Surface(a, location_a);
			const occ::handle<Geom_Surface>& local_b = BRep_Tool::Surface(b, location_b);
			if (!local_a.IsNull() && !local_b.IsNull()
				&& local_a == local_b && location_a.IsEqual(location_b))
			{
				result.exact_confirmation_required = true;
				return result;
			}

			// GeomConvert's predicate is an analytic quadric/quadric equality test, not
			// a point sample.  Use it only as a positive certificate; exact Common still
			// decides whether the two trimmed domains overlap with positive area.
			const occ::handle<Geom_Surface> world_a = BRep_Tool::Surface(a);
			const occ::handle<Geom_Surface> world_b = BRep_Tool::Surface(b);
			if (!world_a.IsNull() && !world_b.IsNull()
				&& GeomConvert_SurfToAnaSurf::IsSame(world_a, world_b, tolerance))
			{
				result.exact_confirmation_required = true;
				return result;
			}

			if (info_a.support_envelope.certified && info_b.support_envelope.certified)
			{
				result.possible_same_domain = std::any_of(
					info_a.support_envelope.patches.begin(),
					info_a.support_envelope.patches.end(), [&](const BezierPatchEnvelope& patch_a)
					{
						return std::any_of(info_b.support_envelope.patches.begin(),
							info_b.support_envelope.patches.end(),
							[&](const BezierPatchEnvelope& patch_b)
							{
								return envelopes_intersect(patch_a, patch_b)
									&& normals_may_be_parallel(patch_a, patch_b);
							});
					});
				return result;
			}

			if (info_a.planar != info_b.planar)
			{
				const PlanarFace& planar = info_a.planar ? info_a : info_b;
				const PlanarFace& curved = info_a.planar ? info_b : info_a;
				if (curved.support_envelope.certified)
				{
					result.possible_same_domain = std::any_of(
						curved.support_envelope.patches.begin(),
						curved.support_envelope.patches.end(),
						[&](const BezierPatchEnvelope& patch)
						{
							return envelope_intersects_bounds(patch, planar.bounds)
								&& patch_can_lie_on_plane(patch, planar.plane, tolerance);
						});
					return result;
				}
			}

		}
		catch (const Standard_Failure& failure)
		{
			result.warning = std::string("support-envelope test raised OpenCascade "
				"exception; retaining exact-Common candidate: ") + failure.GetMessageString();
		}
		catch (const std::exception& exception)
		{
			result.warning = std::string("support-envelope test raised exception; "
				"retaining exact-Common candidate: ") + exception.what();
		}
		return result;
	}

	struct Detection
	{
		std::vector<OverlapPair> pairs;
		std::vector<UnsupportedOverlap> unsupported;
		std::vector<std::vector<std::size_t>> clusters;
		std::vector<std::string> warnings;
		std::vector<std::string> errors;
		std::size_t broad_phase_candidates = 0;
		std::size_t coplanar_candidates = 0;
		std::size_t nonplanar_screened_pairs = 0;
		std::size_t nonplanar_support_rejected_pairs = 0;
		std::size_t nonplanar_unverified_pairs = 0;
		std::size_t nonplanar_exact_common_pairs = 0;
		double support_screen_milliseconds = 0.0;
	};

	Detection detect_clusters(const std::vector<TopoDS_Face>& faces,
		const std::vector<PlanarFace>& info)
	{
		Detection result;
		DisjointSet components(faces.size());
		for (std::size_t a = 0; a < faces.size(); ++a)
			for (std::size_t b = a + 1; b < faces.size(); ++b)
			{
				if (info[a].bounds.IsOut(info[b].bounds)) continue;
				++result.broad_phase_candidates;
				const bool coplanar = same_plane(info[a], info[b]);
				bool nonplanar_candidate = false;
				if (coplanar) ++result.coplanar_candidates;
				else if (!info[a].planar || !info[b].planar)
				{
					nonplanar_candidate = true;
					++result.nonplanar_screened_pairs;
					const double support_tolerance = std::max(Precision::Confusion(),
						info[a].tolerance + info[b].tolerance);
					const auto screen_begin = Clock::now();
					const SupportScreen screen = plausible_same_surface(faces[a], faces[b],
						info[a], info[b], support_tolerance);
					result.support_screen_milliseconds += std::chrono::duration<double,
						std::milli>(Clock::now() - screen_begin).count();
					if (!screen.warning.empty())
					{
						std::ostringstream warning;
						warning << "support screen for source faces " << a << " / " << b
							<< ": " << screen.warning;
						result.warnings.push_back(warning.str());
					}
					if (!screen.possible_same_domain)
					{
						++result.nonplanar_support_rejected_pairs;
						continue;
					}
					if (!screen.exact_confirmation_required)
					{
						++result.nonplanar_unverified_pairs;
						continue;
					}
					++result.nonplanar_exact_common_pairs;
				}
				else
					continue; // distinct planar supports cannot have a positive-area common region

				const BooleanCommon common = exact_common(faces[a], faces[b]);
				if (!common.completed)
				{
					std::ostringstream error;
					error << "exact common failed while "
						<< (coplanar ? "detecting coplanar overlap" :
							"screening unsupported non-planar overlap")
						<< " for source faces " << a << " / " << b;
					if (!common.report.empty()) error << ": " << common.report;
					result.errors.push_back(error.str());
					continue;
				}
				const double length_scale = std::sqrt(std::max(info[a].measure.area,
					info[b].measure.area));
				const double area_tolerance = numeric_area_tolerance(
					std::max(info[a].measure.area, info[b].measure.area),
					info[a].perimeter + info[b].perimeter,
					std::max(info[a].tolerance, info[b].tolerance), length_scale);
				if (!(common.measure.area > area_tolerance)) continue;
				if (nonplanar_candidate)
				{
					result.unsupported.push_back({ a, b, common.measure, area_tolerance });
					std::ostringstream error;
					error << std::setprecision(12)
						<< "unsupported positive-area overlap involving non-planar source faces "
						<< a << " / " << b << ": common area " << common.measure.area
						<< " exceeds decision tolerance " << area_tolerance;
					result.errors.push_back(error.str());
				}
				else
				{
					result.pairs.push_back({ a, b, common.measure, area_tolerance,
						common.warnings ? common.report : std::string{} });
					components.join(a, b);
				}
				if (common.warnings)
				{
					std::ostringstream warning;
					warning << "exact Common warning for source faces " << a << " / " << b;
					if (!common.report.empty()) warning << ": " << common.report;
					result.warnings.push_back(warning.str());
				}
			}
		if (result.nonplanar_unverified_pairs > 0)
		{
			std::ostringstream warning;
			warning << "generic non-planar overlap preflight is incomplete for "
				<< result.nonplanar_unverified_pairs
				<< " candidate pair(s); planar normalization proceeded without claiming "
					"those pairs are overlap-free";
			result.warnings.push_back(warning.str());
		}

		std::map<std::size_t, std::vector<std::size_t>> by_root;
		for (const OverlapPair& pair : result.pairs)
		{
			by_root[components.find(pair.a)].push_back(pair.a);
			by_root[components.find(pair.b)].push_back(pair.b);
		}
		for (auto& [root, cluster] : by_root)
		{
			(void)root;
			std::sort(cluster.begin(), cluster.end());
			cluster.erase(std::unique(cluster.begin(), cluster.end()), cluster.end());
			if (cluster.size() > 1) result.clusters.push_back(std::move(cluster));
		}
		std::sort(result.clusters.begin(), result.clusters.end());
		return result;
	}

	struct Accumulation
	{
		long double area = 0.0L;
		long double moment_x = 0.0L, moment_y = 0.0L, moment_z = 0.0L;
	};

	void append(Accumulation& sum, const SurfaceMeasure& measure, Vec3 reference)
	{
		const long double area = measure.area;
		sum.area += area;
		sum.moment_x += area * static_cast<long double>(measure.centroid.x - reference.x);
		sum.moment_y += area * static_cast<long double>(measure.centroid.y - reference.y);
		sum.moment_z += area * static_cast<long double>(measure.centroid.z - reference.z);
	}

	long double moment_difference(const Accumulation& a, const Accumulation& b)
	{
		const long double x = a.moment_x - b.moment_x;
		const long double y = a.moment_y - b.moment_y;
		const long double z = a.moment_z - b.moment_z;
		return std::sqrt(x * x + y * y + z * z);
	}

	Vec3 box_centre(const std::vector<std::size_t>& source_ids,
		const std::vector<PlanarFace>& info, double& length_scale)
	{
		double lo[3] = { std::numeric_limits<double>::max(),
			std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
		double hi[3] = { -std::numeric_limits<double>::max(),
			-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() };
		for (std::size_t source : source_ids)
		{
			if (info[source].bounds.IsVoid()) continue;
			double face_lo[3]{}, face_hi[3]{};
			info[source].bounds.Get(face_lo[0], face_lo[1], face_lo[2],
				face_hi[0], face_hi[1], face_hi[2]);
			for (int axis = 0; axis < 3; ++axis)
			{
				lo[axis] = std::min(lo[axis], face_lo[axis]);
				hi[axis] = std::max(hi[axis], face_hi[axis]);
			}
		}
		const Vec3 centre{ 0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]),
			0.5 * (lo[2] + hi[2]) };
		length_scale = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0])
			+ (hi[1] - lo[1]) * (hi[1] - lo[1])
			+ (hi[2] - lo[2]) * (hi[2] - lo[2]));
		return centre;
	}

	struct Contributor
	{
		std::size_t source_face = 0;
		int orientation_parity = 1;
		bool constructed_membership = false;
	};

	struct Atom
	{
		TopoDS_Face face;
		SurfaceMeasure measure{};
		double perimeter = 0.0;
		double tolerance = 0.0;
		std::vector<Contributor> contributors;
	};

	struct ClusterResult
	{
		bool valid = false;
		std::vector<std::size_t> source_faces;
		std::vector<Atom> atoms;
		std::vector<std::string> failures;
		std::string warnings;
		double input_area_sum = 0.0;
		double physical_union_area = 0.0;
		double input_max_tolerance = 0.0;
		double output_max_tolerance = 0.0;
		double output_tolerance_limit = 0.0;
		double maximum_rehost_area_error = 0.0;
		double maximum_rehost_moment_error = 0.0;
		double maximum_source_area_error = 0.0;
		double maximum_source_moment_error = 0.0;
		double maximum_atom_overlap_area = 0.0;
		double milliseconds = 0.0;
	};

	ClusterResult normalize_cluster(const std::vector<std::size_t>& source_ids,
		const std::vector<TopoDS_Face>& all_faces, const std::vector<PlanarFace>& info)
	{
		ClusterResult result;
		result.source_faces = source_ids;
		if (source_ids.size() < 2)
		{
			result.failures.push_back("overlap cluster contains fewer than two source faces");
			return result;
		}

		double length_scale = 0.0;
		const Vec3 reference = box_centre(source_ids, info, length_scale);
		const Vec3 physical_normal = canonical_direction(info[source_ids.front()].surface_normal);
		for (std::size_t source : source_ids)
		{
			result.input_area_sum += info[source].measure.area;
			result.input_max_tolerance = std::max(result.input_max_tolerance,
				info[source].tolerance);
			if (!same_plane(info[source_ids.front()], info[source]))
				result.failures.push_back("transitive overlap cluster is not on one native-tolerance plane");
			if (std::abs(dot(info[source].oriented_normal, physical_normal))
				< 1.0 - 256.0 * std::numeric_limits<double>::epsilon())
				result.failures.push_back("source face has no stable orientation parity to cluster plane");
		}
		if (!result.failures.empty()) return result;
		std::vector<TopoDS_Face> canonical_sources;
		canonical_sources.reserve(source_ids.size());
		for (std::size_t local_source = 0; local_source < source_ids.size(); ++local_source)
		{
			const std::size_t source_id = source_ids[local_source];
			const bool reverse_trim_wires = dot(info[source_id].oriented_normal,
				info[source_ids.front()].surface_normal) < 0.0;
			TopoDS_Face rebuilt = rehost_planar_face(all_faces[source_id],
				info[source_ids.front()].plane, reverse_trim_wires);
			if (!BRepCheck_Analyzer(rebuilt, true).IsValid())
				result.failures.push_back("source face " + std::to_string(source_id)
					+ " is invalid after canonical-plane rehosting");
			const SurfaceMeasure rebuilt_measure = surface_measure(rebuilt);
			const double area_tolerance = numeric_area_tolerance(info[source_id].measure.area,
				info[source_id].perimeter, std::max(info[source_id].tolerance,
					maximum_tolerance(rebuilt)), length_scale, 2);
			Accumulation expected, actual;
			append(expected, info[source_id].measure, reference);
			append(actual, rebuilt_measure, reference);
			const double area_error = static_cast<double>(std::abs(expected.area - actual.area));
			const double moment_error = static_cast<double>(moment_difference(expected, actual));
			result.maximum_rehost_area_error = std::max(result.maximum_rehost_area_error,
				area_error);
			result.maximum_rehost_moment_error = std::max(result.maximum_rehost_moment_error,
				moment_error);
			if (area_error > area_tolerance
				|| moment_error > 4.0 * area_tolerance
					* std::max(1.0, length_scale))
			{
				std::ostringstream failure;
				failure << std::setprecision(12)
					<< "canonical-plane rehosting changed coverage of source face " << source_id
					<< ": area expected/actual " << static_cast<double>(expected.area) << " / "
					<< static_cast<double>(actual.area) << " mm^2, moment error "
					<< moment_error;
				result.failures.push_back(failure.str());
			}
			canonical_sources.push_back(std::move(rebuilt));
		}
		if (!result.failures.empty()) return result;

		// Build the Boolean OR as explicit membership cells. General Fuse alone is not
		// sufficient here: on Two-Cells-fixed OCCT legitimately returns the dominant Rib 0
		// image while omitting Rib 1's narrow but positive-area unique strip. Binary Common
		// and Cut operations make every membership cell explicit, and the conservation gates
		// below reject any Boolean operation that still loses such a strip.
		struct DraftAtom
		{
			TopoDS_Face face;
			std::set<std::size_t> local_sources;
		};
		std::vector<DraftAtom> drafts{ { canonical_sources.front(), { 0 } } };
		auto append_warning = [&](const std::string& warning)
		{
			if (warning.empty()) return;
			if (!result.warnings.empty()) result.warnings += '\n';
			result.warnings += warning;
		};
		auto operation_area_tolerance = [&](const TopoDS_Face& a, const TopoDS_Face& b)
		{
			const PlanarFace a_info = inspect_face(a), b_info = inspect_face(b);
			return numeric_area_tolerance(std::max(a_info.measure.area, b_info.measure.area),
				a_info.perimeter + b_info.perimeter,
				std::max(a_info.tolerance, b_info.tolerance), length_scale);
		};
		auto append_parts = [&](std::vector<DraftAtom>& destination,
			const BooleanParts& parts, const std::set<std::size_t>& contributors,
			const char* operation)
		{
			if (!parts.completed)
			{
				result.failures.push_back(std::string(operation) + " failed: " + parts.report);
				return;
			}
			append_warning(parts.report);
			for (const TopoDS_Face& face : parts.faces)
			{
				const SurfaceMeasure measure = surface_measure(face);
				if (!(measure.area > 0.0) || !std::isfinite(measure.area)) continue;
				destination.push_back({ face, contributors });
			}
		};

		const auto begin = Clock::now();
		for (std::size_t local_source = 1; local_source < source_ids.size(); ++local_source)
		{
			const TopoDS_Face& source = canonical_sources[local_source];
			const std::vector<DraftAtom> previous = std::move(drafts);
			std::vector<DraftAtom> next;

			// Split every previously accepted disjoint atom into old-only and shared cells.
			for (const DraftAtom& old : previous)
			{
				const BooleanParts common = exact_boolean_parts<BRepAlgoAPI_Common>(old.face, source);
				if (!common.completed)
				{
					result.failures.push_back("explicit membership Common failed: " + common.report);
					continue;
				}
				append_warning(common.report);
				if (!(common.measure.area > operation_area_tolerance(old.face, source)))
				{
					next.push_back(old);
					continue;
				}
				const BooleanParts old_only = exact_boolean_parts<BRepAlgoAPI_Cut>(old.face, source);
				append_parts(next, old_only, old.local_sources, "old-only Cut");
				std::set<std::size_t> shared_sources = old.local_sources;
				shared_sources.insert(local_source);
				append_parts(next, common, shared_sources, "shared Common");
			}

			// Subtract the complete old union from the new source. Previous atoms are already
			// pairwise disjoint, so sequential exact cuts retain precisely the new-only cells.
			std::vector<TopoDS_Face> remainder{ source };
			for (const DraftAtom& old : previous)
			{
				std::vector<TopoDS_Face> cut_remainder;
				for (const TopoDS_Face& piece : remainder)
				{
					const BooleanParts common = exact_boolean_parts<BRepAlgoAPI_Common>(piece, old.face);
					if (!common.completed)
					{
						result.failures.push_back("new-only membership Common failed: "
							+ common.report);
						continue;
					}
					append_warning(common.report);
					if (!(common.measure.area > operation_area_tolerance(piece, old.face)))
					{
						cut_remainder.push_back(piece);
						continue;
					}
					const BooleanParts difference = exact_boolean_parts<BRepAlgoAPI_Cut>(piece,
						old.face);
					if (!difference.completed)
					{
						result.failures.push_back("new-only Cut failed: " + difference.report);
						continue;
					}
					append_warning(difference.report);
					for (const TopoDS_Face& face : difference.faces)
						if (surface_measure(face).area > 0.0) cut_remainder.push_back(face);
				}
				remainder = std::move(cut_remainder);
				if (remainder.empty()) break;
			}
			for (const TopoDS_Face& face : remainder)
				next.push_back({ face, { local_source } });
			drafts = std::move(next);
		}
		result.milliseconds = std::chrono::duration<double, std::milli>(
			Clock::now() - begin).count();
		if (drafts.empty())
		{
			result.failures.push_back("explicit all-cells union produced no atomic faces");
			return result;
		}

		std::vector<TopoDS_Face> draft_faces;
		draft_faces.reserve(drafts.size());
		for (const DraftAtom& draft : drafts) draft_faces.push_back(draft.face);
		if (!BRepCheck_Analyzer(face_compound(draft_faces), true).IsValid())
			result.failures.push_back("explicit all-cells union is not a valid face compound");

		result.atoms.reserve(drafts.size());
		for (std::size_t atom_index = 0; atom_index < drafts.size(); ++atom_index)
		{
			TopoDS_Face atom_face = drafts[atom_index].face;
			PlanarFace atom_info = inspect_face(atom_face);
			if (!atom_info.planar || !same_plane(info[source_ids.front()], atom_info))
			{
				result.failures.push_back("result atom " + std::to_string(atom_index)
					+ " is not on the cluster plane");
				continue;
			}
			if (!BRepCheck_Analyzer(atom_face, true).IsValid())
				result.failures.push_back("result atom " + std::to_string(atom_index)
					+ " is not a valid face");
			if (dot(atom_info.oriented_normal, physical_normal) < 0.0)
			{
				atom_face.Reverse();
				atom_info.oriented_normal = atom_info.oriented_normal * -1.0;
			}

			Atom atom;
			atom.face = atom_face;
			atom.measure = atom_info.measure;
			atom.perimeter = atom_info.perimeter;
			atom.tolerance = atom_info.tolerance;
			result.output_max_tolerance = std::max(result.output_max_tolerance, atom.tolerance);

			std::set<std::size_t> geometric_sources;
			for (std::size_t local_source = 0; local_source < source_ids.size(); ++local_source)
			{
				const std::size_t source_id = source_ids[local_source];
				const BooleanCommon common = exact_common(atom.face, all_faces[source_id]);
				if (!common.completed)
				{
					std::string failure="atom/source coverage Common failed for atom "
						+std::to_string(atom_index)+" and source "+std::to_string(source_id);
					if(!common.report.empty())failure+=": "+common.report;
					result.failures.push_back(std::move(failure));
					continue;
				}
				const double area_tolerance = numeric_area_tolerance(atom.measure.area,
					atom.perimeter + info[source_id].perimeter,
					std::max(atom.tolerance, info[source_id].tolerance), length_scale, 2);
				Accumulation atom_sum, common_sum;
				append(atom_sum, atom.measure, reference);
				append(common_sum, common.measure, reference);
				const double moment_tolerance = 4.0 * area_tolerance
					* std::max(1.0, length_scale);
				if (std::abs(common.measure.area - atom.measure.area) <= area_tolerance
					&& moment_difference(common_sum, atom_sum) <= moment_tolerance)
					geometric_sources.insert(local_source);
			}

			const std::set<std::size_t>& constructed = drafts[atom_index].local_sources;
			if (geometric_sources.empty())
				result.failures.push_back("result atom " + std::to_string(atom_index)
					+ " adds coverage not fully owned by any source face");
			if (constructed != geometric_sources)
			{
				std::ostringstream mismatch;
				mismatch << "constructed/geometric provenance mismatch for atom " << atom_index
					<< " (constructed";
				for (std::size_t source : constructed) mismatch << ' ' << source_ids[source];
				mismatch << "; geometric";
				for (std::size_t source : geometric_sources) mismatch << ' ' << source_ids[source];
				mismatch << ')';
				result.failures.push_back(mismatch.str());
			}
			for (std::size_t local_source : geometric_sources)
			{
				const std::size_t source_id = source_ids[local_source];
				const int parity = dot(info[source_id].oriented_normal, physical_normal) >= 0.0
					? 1 : -1;
				atom.contributors.push_back({ source_id, parity,
					constructed.contains(local_source) });
			}
			result.physical_union_area += atom.measure.area;
			result.atoms.push_back(std::move(atom));
		}

		// Boolean intersection vertices are recomputed on the same exact plane. Permit only
		// a coordinate-roundoff-sized increase; this remains orders of magnitude below any
		// geometric healing tolerance and catches the global split's observed inflation.
		const double tolerance_growth_allowance = 256.0
			* std::numeric_limits<double>::epsilon() * std::max(1.0, length_scale);
		const double tolerance_limit = std::max(result.input_max_tolerance,
			Precision::Confusion()) + tolerance_growth_allowance;
		result.output_tolerance_limit = tolerance_limit;
		if (result.output_max_tolerance > tolerance_limit)
		{
			std::ostringstream failure;
			failure << std::setprecision(12) << "output tolerance "
				<< result.output_max_tolerance << " mm exceeds native bound "
				<< tolerance_limit << " mm (input max " << result.input_max_tolerance << " mm)";
			result.failures.push_back(failure.str());
		}

		// A source must be exactly reconstructible from the non-overlapping atoms that name it.
		for (std::size_t source_id : source_ids)
		{
			Accumulation expected, reconstructed;
			append(expected, info[source_id].measure, reference);
			std::size_t contributing_atoms = 0;
			for (const Atom& atom : result.atoms)
				if (std::any_of(atom.contributors.begin(), atom.contributors.end(),
					[source_id](const Contributor& contributor)
					{ return contributor.source_face == source_id; }))
				{
					append(reconstructed, atom.measure, reference);
					++contributing_atoms;
				}
			const double area_tolerance = numeric_area_tolerance(info[source_id].measure.area,
				info[source_id].perimeter, std::max(result.input_max_tolerance,
					result.output_max_tolerance), length_scale, 1 + contributing_atoms);
			const double moment_tolerance = 4.0 * area_tolerance * std::max(1.0, length_scale);
			const double area_error = static_cast<double>(std::abs(
				reconstructed.area - expected.area));
			const double moment_error = static_cast<double>(moment_difference(
				reconstructed, expected));
			result.maximum_source_area_error = std::max(result.maximum_source_area_error,
				area_error);
			result.maximum_source_moment_error = std::max(result.maximum_source_moment_error,
				moment_error);
			if (area_error > area_tolerance || moment_error > moment_tolerance)
			{
				std::ostringstream failure;
				failure << std::setprecision(12) << "source face " << source_id
					<< " is not reconstructed by contributor atoms: area expected/actual "
					<< static_cast<double>(expected.area) << " / "
					<< static_cast<double>(reconstructed.area) << " mm^2, moment error "
					<< moment_error;
				result.failures.push_back(failure.str());
			}
		}

		// The explicit cell arrangement must contain each physical region once, not once per source.
		for (std::size_t a = 0; a < result.atoms.size(); ++a)
			for (std::size_t b = a + 1; b < result.atoms.size(); ++b)
			{
				const BooleanCommon common = exact_common(result.atoms[a].face, result.atoms[b].face);
				if (!common.completed)
				{
					std::string failure="atom/atom disjointness Common failed for atoms "
						+std::to_string(a)+" / "+std::to_string(b);
					if(!common.report.empty())failure+=": "+common.report;
					result.failures.push_back(std::move(failure));
					continue;
				}
				const double area_tolerance = numeric_area_tolerance(std::max(
					result.atoms[a].measure.area, result.atoms[b].measure.area),
					result.atoms[a].perimeter + result.atoms[b].perimeter,
					std::max(result.atoms[a].tolerance, result.atoms[b].tolerance),
					length_scale);
				result.maximum_atom_overlap_area = std::max(result.maximum_atom_overlap_area,
					common.measure.area);
				if (common.measure.area > area_tolerance)
				{
					std::ostringstream failure;
					failure << std::setprecision(12) << "result atoms " << a << " / " << b
						<< " retain positive-area overlap " << common.measure.area << " mm^2";
					result.failures.push_back(failure.str());
				}
			}

		std::stable_sort(result.atoms.begin(), result.atoms.end(), [](const Atom& a, const Atom& b)
		{
			auto ids = [](const Atom& atom)
			{
				std::vector<std::size_t> result;
				for (const Contributor& contributor : atom.contributors)
					result.push_back(contributor.source_face);
				return result;
			};
			const std::vector<std::size_t> a_ids = ids(a), b_ids = ids(b);
			if (a_ids != b_ids) return a_ids < b_ids;
			return std::tie(a.measure.centroid.x, a.measure.centroid.y, a.measure.centroid.z,
				a.measure.area, a.perimeter, a.tolerance)
				< std::tie(b.measure.centroid.x, b.measure.centroid.y,
					b.measure.centroid.z, b.measure.area, b.perimeter, b.tolerance);
		});

		result.valid = result.failures.empty();
		return result;
	}

	} // namespace

	OcctSurfaceNormalizationResult normalize_occt_surface_faces(
		const std::vector<TopoDS_Face>& source_faces)
	{
		OcctSurfaceNormalizationResult result;
		result.stats.source_face_count = source_faces.size();
		try
		{
			std::vector<PlanarFace> info;
			info.reserve(source_faces.size());
			for (std::size_t source = 0; source < source_faces.size(); ++source)
			{
				if (source_faces[source].IsNull())
				{
					result.errors.push_back("source face " + std::to_string(source)
						+ " is null");
					return result;
				}
				if(!BRepCheck_Analyzer(source_faces[source],true).IsValid())
				{
					result.errors.push_back("source face "+std::to_string(source)
						+" is not a valid BRep face");
					return result;
				}
				info.push_back(inspect_face(source_faces[source]));
				if(!(info.back().measure.area>0.0)
					||!std::isfinite(info.back().measure.area))
				{
					result.errors.push_back("source face "+std::to_string(source)
						+" has non-positive or non-finite area");
					return result;
				}
				result.stats.planar_source_face_count += info.back().planar;
				result.stats.input_area_sum += info.back().measure.area;
			}

			const auto detection_begin = Clock::now();
			const Detection detection = detect_clusters(source_faces, info);
			result.stats.detection_milliseconds = std::chrono::duration<double, std::milli>(
				Clock::now() - detection_begin).count();
			result.stats.broad_phase_candidate_pair_count = detection.broad_phase_candidates;
			result.stats.coplanar_candidate_pair_count = detection.coplanar_candidates;
			result.stats.nonplanar_screened_pair_count = detection.nonplanar_screened_pairs;
			result.stats.nonplanar_support_rejected_pair_count =
				detection.nonplanar_support_rejected_pairs;
			result.stats.nonplanar_unverified_pair_count =
				detection.nonplanar_unverified_pairs;
			result.stats.nonplanar_exact_common_pair_count =
				detection.nonplanar_exact_common_pairs;
			result.stats.support_screen_milliseconds = detection.support_screen_milliseconds;
			result.stats.positive_area_overlap_pair_count = detection.pairs.size();
			result.stats.overlap_cluster_count = detection.clusters.size();
			result.warnings = detection.warnings;
			result.errors = detection.errors;

			result.overlap_pairs.reserve(detection.pairs.size());
			for (const OverlapPair& pair : detection.pairs)
			{
				OcctSurfaceOverlapPairAudit audit;
				audit.source_a = pair.a;
				audit.source_b = pair.b;
				audit.common_area = pair.common.area;
				audit.common_centroid = { pair.common.centroid.x, pair.common.centroid.y,
					pair.common.centroid.z };
				audit.positive_area_tolerance = pair.area_tolerance;
				if (!pair.warnings.empty()) audit.warnings.push_back(pair.warnings);
				result.overlap_pairs.push_back(std::move(audit));
			}
			result.unsupported_overlaps.reserve(detection.unsupported.size());
			for (const UnsupportedOverlap& overlap : detection.unsupported)
			{
				result.unsupported_overlaps.push_back({ overlap.a, overlap.b,
					info[overlap.a].planar, info[overlap.b].planar, overlap.common.area,
					{ overlap.common.centroid.x, overlap.common.centroid.y,
						overlap.common.centroid.z }, overlap.area_tolerance });
			}

			std::vector<ClusterResult> normalized_clusters;
			normalized_clusters.reserve(detection.clusters.size());
			result.clusters.reserve(detection.clusters.size());
			for (std::size_t cluster_id = 0; cluster_id < detection.clusters.size(); ++cluster_id)
			{
				ClusterResult normalized;
				try
				{
					normalized = normalize_cluster(detection.clusters[cluster_id],
						source_faces, info);
				}
				catch (const Standard_Failure& failure)
				{
					normalized.source_faces = detection.clusters[cluster_id];
					normalized.failures.push_back(std::string("OpenCascade exception: ")
						+ failure.GetMessageString());
				}
				catch (const std::exception& exception)
				{
					normalized.source_faces = detection.clusters[cluster_id];
					normalized.failures.push_back(std::string("exception: ") + exception.what());
				}

				OcctSurfaceNormalizationClusterAudit audit;
				audit.source_face_ids = normalized.source_faces;
				audit.atomic_face_count = normalized.atoms.size();
				for (const Atom& atom : normalized.atoms)
				{
					if (atom.contributors.size() > 1) ++audit.shared_source_atom_count;
					else ++audit.single_source_atom_count;
				}
				audit.input_area_sum = normalized.input_area_sum;
				audit.once_covered_area = normalized.physical_union_area;
				audit.duplicate_area_removed = normalized.input_area_sum
					- normalized.physical_union_area;
				audit.input_maximum_tolerance = normalized.input_max_tolerance;
				audit.output_maximum_tolerance = normalized.output_max_tolerance;
				audit.output_tolerance_limit = normalized.output_tolerance_limit;
				audit.maximum_rehost_area_error = normalized.maximum_rehost_area_error;
				audit.maximum_rehost_first_moment_error =
					normalized.maximum_rehost_moment_error;
				audit.maximum_source_area_error = normalized.maximum_source_area_error;
				audit.maximum_source_first_moment_error =
					normalized.maximum_source_moment_error;
				audit.maximum_atom_overlap_area = normalized.maximum_atom_overlap_area;
				audit.arrangement_milliseconds = normalized.milliseconds;
				if (!normalized.warnings.empty()) audit.warnings.push_back(normalized.warnings);
				audit.errors = normalized.failures;
				if (!audit.valid())
					for (const std::string& error : audit.errors)
						result.errors.push_back("overlap cluster " + std::to_string(cluster_id)
							+ ": " + error);
				result.stats.arrangement_milliseconds += normalized.milliseconds;
				result.clusters.push_back(std::move(audit));
				normalized_clusters.push_back(std::move(normalized));
			}

			std::vector<std::size_t> source_cluster(source_faces.size(), source_faces.size());
			for (std::size_t cluster = 0; cluster < detection.clusters.size(); ++cluster)
				for (std::size_t source : detection.clusters[cluster])
					source_cluster[source] = cluster;
			result.stats.untouched_face_count = static_cast<std::size_t>(std::count(
				source_cluster.begin(), source_cluster.end(), source_faces.size()));

			// All output is staged until every detection and cluster gate has passed.  This
			// prevents callers from accidentally consuming a partially normalized model.
			if (!result.errors.empty()) return result;
			for (std::size_t source = 0; source < source_faces.size(); ++source)
			{
				const std::size_t cluster = source_cluster[source];
				if (cluster == source_faces.size())
				{
					const PlanarFace& source_info = info[source];
					result.normalized_faces.push_back({ source_faces[source], { { source, 1 } },
						source_info.measure.area, { source_info.measure.centroid.x,
							source_info.measure.centroid.y, source_info.measure.centroid.z },
						source_info.tolerance });
					continue;
				}
				if (source != detection.clusters[cluster].front()) continue;
				for (const Atom& atom : normalized_clusters[cluster].atoms)
				{
					const std::size_t output_id = result.normalized_faces.size();
					OcctNormalizedSurfaceFace output;
					output.face = atom.face;
					output.area = atom.measure.area;
					output.centroid = { atom.measure.centroid.x, atom.measure.centroid.y,
						atom.measure.centroid.z };
					output.maximum_tolerance = atom.tolerance;
					for (const Contributor& contributor : atom.contributors)
						output.contributors.push_back({ contributor.source_face,
							static_cast<std::int8_t>(contributor.orientation_parity) });
					result.normalized_faces.push_back(std::move(output));
					result.clusters[cluster].normalized_face_indices.push_back(output_id);
				}
			}

			result.stats.normalized_face_count = result.normalized_faces.size();
			for (const OcctNormalizedSurfaceFace& face : result.normalized_faces)
				result.stats.once_covered_area += face.area;
			result.stats.duplicate_area_removed = result.stats.input_area_sum
				- result.stats.once_covered_area;
		}
		catch (const Standard_Failure& failure)
		{
			result.normalized_faces.clear();
			result.stats.normalized_face_count = 0;
			result.errors.push_back(std::string("OpenCascade exception: ")
				+ failure.GetMessageString());
		}
		catch (const std::exception& exception)
		{
			result.normalized_faces.clear();
			result.stats.normalized_face_count = 0;
			result.errors.push_back(std::string("exception: ") + exception.what());
		}
		return result;
	}
}
