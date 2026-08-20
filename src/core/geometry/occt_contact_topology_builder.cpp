// Deterministic construction of exact CAD contact atoms from OCCT topology.

#include "core/geometry/occt_contact_topology_builder.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Geom2dAPI_ProjectPointOnCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Curve.hxx>
#include <IntTools_CommonPrt.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_Sequence.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		using ShapeIndexedMap =
			NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
		constexpr std::uint32_t no_id = std::numeric_limits<std::uint32_t>::max();
		constexpr double maximum_parameter_fraction = 1.0e-3;

		bool finite(double value) { return std::isfinite(value); }

		double location_scale(const TopoDS_Shape& shape)
		{
			return std::abs(shape.Location().Transformation().ScaleFactor());
		}

		double edge_tolerance(const TopoDS_Edge& edge)
		{
			return std::max(0.0, BRep_Tool::Tolerance(edge)) * location_scale(edge);
		}

		double vertex_tolerance(const TopoDS_Vertex& vertex)
		{
			return std::max(0.0, BRep_Tool::Tolerance(vertex)) * location_scale(vertex);
		}

		double arithmetic_parameter_tolerance(double first, double last)
		{
			auto ulp = [](double value)
			{
				return std::abs(std::nextafter(value,
					std::numeric_limits<double>::infinity()) - value);
			};
			const double span = std::abs(last - first);
			return std::max({8.0 * ulp(first), 8.0 * ulp(last),
				128.0 * std::numeric_limits<double>::epsilon() * span});
		}

		struct EdgeCurve
		{
			TopoDS_Edge edge;
			BRepAdaptor_Curve curve;
			double first = 0.0;
			double last = 0.0;
			double parameter_tolerance = 0.0;
			double geometric_tolerance = 0.0;
		};

		bool make_edge_curve(const TopoDS_Edge& edge, EdgeCurve& output,
			std::string& error)
		{
			if (edge.IsNull() || BRep_Tool::Degenerated(edge))
			{
				error = "source edge is null or degenerate";
				return false;
			}
			EdgeCurve result;
			result.edge = edge;
			result.curve.Initialize(edge);
			result.first = result.curve.FirstParameter();
			result.last = result.curve.LastParameter();
			result.geometric_tolerance = edge_tolerance(edge);
			if (!finite(result.first) || !finite(result.last)
				|| !(result.last > result.first))
			{
				error = "edge has no finite positive parameter range";
				return false;
			}
			const double arithmetic = arithmetic_parameter_tolerance(result.first,
				result.last);
			const double resolved = std::abs(result.curve.Resolution(
				result.geometric_tolerance + Precision::Confusion()));
			result.parameter_tolerance = finite(resolved)
				? std::max(arithmetic, resolved) : arithmetic;
			if (!finite(result.parameter_tolerance)
				|| result.parameter_tolerance < 0.0
				|| result.parameter_tolerance >= maximum_parameter_fraction
					* (result.last - result.first))
			{
				error = "edge parameterization is unresolved at its native CAD tolerance";
				return false;
			}
			output = std::move(result);
			return true;
		}

		std::array<double, 3> point_array(const gp_Pnt& point)
		{
			return {{point.X(), point.Y(), point.Z()}};
		}

		gp_Pnt point(const std::array<double, 3>& value)
		{
			return {value[0], value[1], value[2]};
		}

		bool parameter_inside(const EdgeCurve& edge, double parameter)
		{
			return finite(parameter)
				&& parameter >= edge.first - edge.parameter_tolerance
				&& parameter <= edge.last + edge.parameter_tolerance;
		}

		bool close_parameter(const EdgeCurve& edge, double a, double b)
		{
			return std::abs(a - b) <= edge.parameter_tolerance;
		}

		bool range_contains(double a, double b, double value, double tolerance)
		{
			if (b < a) std::swap(a, b);
			return value >= a - tolerance && value <= b + tolerance;
		}

		bool valid_orientation(std::int8_t value)
		{
			return value == static_cast<std::int8_t>(TopAbs_FORWARD)
				|| value == static_cast<std::int8_t>(TopAbs_REVERSED);
		}

		bool unique_parameter(const gp_Pnt& world, const EdgeCurve& edge,
			double linear_tolerance, double& parameter, std::string& error)
		{
			TopLoc_Location location;
			Standard_Real curve_first = 0.0, curve_last = 0.0;
			const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge.edge, location,
				curve_first, curve_last);
			if (curve.IsNull())
			{
				error = "edge has no 3D curve for bounded parameter recovery";
				return false;
			}
			gp_Pnt local = world;
			local.Transform(location.Transformation().Inverted());
			GeomAPI_ProjectPointOnCurve projection(local, curve, edge.first, edge.last);
			std::vector<double> candidates;
			for (Standard_Integer index = 1; index <= projection.NbPoints(); ++index)
			{
				double candidate = projection.Parameter(index);
				if (!parameter_inside(edge, candidate)) continue;
				candidate = std::clamp(candidate, edge.first, edge.last);
				gp_Pnt recovered = curve->Value(candidate);
				recovered.Transform(location.Transformation());
				if (recovered.SquareDistance(world)
					<= linear_tolerance * linear_tolerance)
					candidates.push_back(candidate);
			}
			// A periodic endpoint has two valid curve parameters.  Preserve its exact
			// endpoint branch only when the world point itself selects that endpoint.
			for (double endpoint : {edge.first, edge.last})
			{
				const gp_Pnt endpoint_point = edge.curve.Value(endpoint);
				if (endpoint_point.SquareDistance(world)
					<= linear_tolerance * linear_tolerance)
					candidates.push_back(endpoint);
			}
			std::sort(candidates.begin(), candidates.end());
			candidates.erase(std::unique(candidates.begin(), candidates.end(),
				[&](double a, double b)
				{
					return std::abs(a - b) <= edge.parameter_tolerance;
				}), candidates.end());
			if (candidates.size() != 1)
			{
				std::ostringstream report;
				report << "world sample has " << candidates.size()
					<< " native-bounded edge parameters (expected exactly one)";
				error = report.str();
				return false;
			}
			parameter = candidates.front();
			return true;
		}

		enum class BoundedProjection : std::uint8_t
		{
			none,
			unique,
			ambiguous
		};

		// Recover a parameter only on the already-certified source subspan.  In
		// particular, do not project on the complete support curve and then select a
		// nearby periodic/self-near branch: a target tessellation node is useful only
		// when this bounded span identifies exactly one native-bounded source point.
		BoundedProjection unique_parameter_on_span(const gp_Pnt& world,
			const EdgeCurve& edge, double span_begin, double span_end,
			double linear_tolerance, double& parameter)
		{
			if (!finite(world.X()) || !finite(world.Y()) || !finite(world.Z())
				|| !finite(span_begin) || !finite(span_end)
				|| !finite(linear_tolerance) || linear_tolerance < 0.0)
				return BoundedProjection::none;
			if (span_end < span_begin) std::swap(span_begin, span_end);
			span_begin = std::clamp(span_begin, edge.first, edge.last);
			span_end = std::clamp(span_end, edge.first, edge.last);
			if (!(span_end > span_begin)) return BoundedProjection::none;

			TopLoc_Location location;
			Standard_Real curve_first = 0.0, curve_last = 0.0;
			const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge.edge, location,
				curve_first, curve_last);
			if (curve.IsNull()) return BoundedProjection::none;
			gp_Pnt local = world;
			local.Transform(location.Transformation().Inverted());
			GeomAPI_ProjectPointOnCurve projection(local, curve, span_begin, span_end);
			std::vector<double> candidates;
			for (Standard_Integer index = 1; index <= projection.NbPoints(); ++index)
			{
				double candidate = projection.Parameter(index);
				if (!finite(candidate)
					|| !range_contains(span_begin, span_end, candidate,
						edge.parameter_tolerance))
					continue;
				candidate = std::clamp(candidate, span_begin, span_end);
				gp_Pnt recovered = curve->Value(candidate);
				recovered.Transform(location.Transformation());
				if (recovered.SquareDistance(world)
					<= linear_tolerance * linear_tolerance)
					candidates.push_back(candidate);
			}
			for (double endpoint : {span_begin, span_end})
			{
				const gp_Pnt endpoint_point = edge.curve.Value(endpoint);
				if (endpoint_point.SquareDistance(world)
					<= linear_tolerance * linear_tolerance)
					candidates.push_back(endpoint);
			}
			std::sort(candidates.begin(), candidates.end());
			candidates.erase(std::unique(candidates.begin(), candidates.end(),
				[&](double a, double b)
				{
					return std::abs(a - b) <= edge.parameter_tolerance;
				}), candidates.end());
			if (candidates.empty()) return BoundedProjection::none;
			if (candidates.size() != 1) return BoundedProjection::ambiguous;
			parameter = candidates.front();
			return BoundedProjection::unique;
		}

		bool same_face(const TopoDS_Face& a, const TopoDS_Face& b)
		{
			return !a.IsNull() && !b.IsNull() && a.IsSame(b);
		}

		double polygon_pcurve_error(
			const Handle(Poly_PolygonOnTriangulation)& polygon,
			const Handle(Poly_Triangulation)& triangulation,
			const Handle(Geom2d_Curve)& pcurve, double first, double last)
		{
			if (polygon.IsNull() || triangulation.IsNull() || pcurve.IsNull())
				return std::numeric_limits<double>::infinity();
			double maximum = 0.0;
			for (Standard_Integer node = 1; node <= polygon->NbNodes(); ++node)
			{
				const Standard_Integer triangulation_node = polygon->Node(node);
				if (triangulation_node <= 0
					|| triangulation_node > triangulation->NbNodes())
					return std::numeric_limits<double>::infinity();
				const gp_Pnt2d uv = triangulation->UVNode(triangulation_node);
				double distance = std::numeric_limits<double>::infinity();
				if (polygon->HasParameters())
				{
					const double parameter = polygon->Parameter(node);
					if (!finite(parameter) || parameter < first || parameter > last)
						return std::numeric_limits<double>::infinity();
					distance = uv.Distance(pcurve->Value(parameter));
				}
				else
				{
					Geom2dAPI_ProjectPointOnCurve projection(uv, pcurve, first, last);
					if (projection.NbPoints() < 1)
						return std::numeric_limits<double>::infinity();
					distance = projection.LowerDistance();
				}
				if (!finite(distance)) return std::numeric_limits<double>::infinity();
				maximum = std::max(maximum, distance);
			}
			return maximum;
		}

		bool occurrence_polygon_for_pcurve(const TopoDS_Edge& occurrence,
			const TopoDS_Face& face,
			const Handle(Poly_Triangulation)& triangulation,
			const TopLoc_Location& location,
			Handle(Poly_PolygonOnTriangulation)& selected, std::string& error)
		{
			double first = 0.0, last = 0.0;
			const Handle(Geom2d_Curve) pcurve =
				BRep_Tool::CurveOnSurface(occurrence, face, first, last);
			if (pcurve.IsNull() || !finite(first) || !finite(last) || !(last > first))
			{
				error = "boundary occurrence has no finite exact face pcurve";
				return false;
			}
			const Handle(Poly_PolygonOnTriangulation) direct =
				BRep_Tool::PolygonOnTriangulation(occurrence, triangulation, location);
			TopoDS_Edge reversed = occurrence;
			reversed.Reverse();
			const Handle(Poly_PolygonOnTriangulation) opposite =
				BRep_Tool::PolygonOnTriangulation(reversed, triangulation, location);
			if (direct.IsNull() && opposite.IsNull())
			{
				selected.Nullify();
				return true;
			}
			if (direct.IsNull() || direct == opposite)
			{
				selected = direct.IsNull() ? opposite : direct;
				return true;
			}
			const double direct_error = polygon_pcurve_error(direct, triangulation,
				pcurve, first, last);
			const double opposite_error = polygon_pcurve_error(opposite, triangulation,
				pcurve, first, last);
			if (!finite(direct_error) && !finite(opposite_error))
			{
				error = "neither seam polygon follows the occurrence's exact face pcurve";
				return false;
			}
			if (direct_error == opposite_error)
			{
				error = "seam polygon branch is ambiguous in the exact face chart";
				return false;
			}
			selected = direct_error < opposite_error ? direct : opposite;
			return true;
		}

		double face_uv_identity_tolerance(const TopoDS_Face& face,
			const std::array<double, 2>& uv, double exact_position_tolerance)
		{
			const double scale = std::max({std::abs(uv[0]), std::abs(uv[1]), 1.0});
			const double arithmetic = 128.0
				* std::numeric_limits<double>::epsilon() * scale;
			const double position_tolerance = std::max({Precision::Confusion(),
				std::max(0.0, BRep_Tool::Tolerance(face)) * location_scale(face),
				exact_position_tolerance});
			BRepAdaptor_Surface surface(face, true);
			const double u_resolution = surface.UResolution(position_tolerance);
			const double v_resolution = surface.VResolution(position_tolerance);
			double resolution = std::numeric_limits<double>::infinity();
			if (finite(u_resolution) && u_resolution > 0.0)
				resolution = std::min(resolution, u_resolution);
			if (finite(v_resolution) && v_resolution > 0.0)
				resolution = std::min(resolution, v_resolution);
			if (!finite(resolution)) resolution = 0.0;
			return arithmetic + resolution;
		}

		struct JunctionCandidate
		{
			std::uint64_t stable_key = 0;
			std::uint8_t kind = 0; // topology vertex before analytic edge/edge vertex
			std::array<double, 3> position{};
			double tolerance = 0.0;
			std::vector<OcctContactJunctionIncidence> incidences;
		};

		struct MeshNodeClaimSeed
		{
			std::uint32_t source_face_id = no_id;
			std::uint32_t source_node_index = no_id;
			std::uint32_t boundary_occurrence_id =
				OcctContactBoundaryOccurrence::no_occurrence_id;
			std::uint32_t source_edge_id = no_id;
			double source_parameter = 0.0;
			std::array<double, 2> face_uv{};
			std::array<double, 3> world_position{};
			double position_tolerance = 0.0;
		};

		bool incidence_less(const OcctContactJunctionIncidence& a,
			const OcctContactJunctionIncidence& b)
		{
			return std::tie(a.source_edge_id, a.parameter)
				< std::tie(b.source_edge_id, b.parameter);
		}

		bool incidences_equal_exact(
			const std::vector<OcctContactJunctionIncidence>& a,
			const std::vector<OcctContactJunctionIncidence>& b)
		{
			return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
				[](const auto& left, const auto& right)
				{
					return left.source_edge_id == right.source_edge_id
						&& left.parameter == right.parameter;
				});
		}

		struct EdgeBounds
		{
			std::size_t source_index = 0;
			std::array<double, 3> low{};
			std::array<double, 3> high{};
		};

		bool make_bounds(const TopoDS_Edge& edge, std::size_t source_index,
			EdgeBounds& output, std::string& error)
		{
			Bnd_Box box;
			BRepBndLib::AddOptimal(edge, box, false, true);
			if (box.IsVoid() || box.IsWhole())
			{
				error = "could not build a finite exact edge AABB";
				return false;
			}
			Standard_Real xmin = 0.0, ymin = 0.0, zmin = 0.0;
			Standard_Real xmax = 0.0, ymax = 0.0, zmax = 0.0;
			box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
			if (!finite(xmin) || !finite(ymin) || !finite(zmin)
				|| !finite(xmax) || !finite(ymax) || !finite(zmax))
			{
				error = "exact edge AABB contains a non-finite coordinate";
				return false;
			}
			output = {source_index, {{xmin, ymin, zmin}}, {{xmax, ymax, zmax}}};
			return true;
		}

		bool overlaps_yz(const EdgeBounds& a, const EdgeBounds& b)
		{
			return a.low[1] <= b.high[1] && b.low[1] <= a.high[1]
				&& a.low[2] <= b.high[2] && b.low[2] <= a.high[2];
		}

		using EdgePair = std::pair<std::uint32_t, std::uint32_t>;

		EdgePair ordered_pair(std::uint32_t a, std::uint32_t b)
		{
			return a < b ? EdgePair{a, b} : EdgePair{b, a};
		}

		void canonicalize_reciprocal(OcctContactReciprocalInterval& interval)
		{
			if (interval.source_edge_b < interval.source_edge_a)
			{
				std::swap(interval.source_edge_a, interval.source_edge_b);
				std::swap(interval.parameter_a_begin, interval.parameter_b_begin);
				std::swap(interval.parameter_a_end, interval.parameter_b_end);
			}
			const auto forward = std::make_tuple(interval.parameter_a_begin,
				interval.parameter_a_end, interval.parameter_b_begin,
				interval.parameter_b_end);
			const auto reversed = std::make_tuple(interval.parameter_a_end,
				interval.parameter_a_begin, interval.parameter_b_end,
				interval.parameter_b_begin);
			if (reversed < forward)
			{
				std::swap(interval.parameter_a_begin, interval.parameter_a_end);
				std::swap(interval.parameter_b_begin, interval.parameter_b_end);
			}
		}

		bool reciprocal_less(const OcctContactReciprocalInterval& a,
			const OcctContactReciprocalInterval& b)
		{
			return std::tie(a.source_edge_a, a.source_edge_b,
				a.parameter_a_begin, a.parameter_a_end,
				a.parameter_b_begin, a.parameter_b_end,
				a.geometric_tolerance)
				< std::tie(b.source_edge_a, b.source_edge_b,
					b.parameter_a_begin, b.parameter_a_end,
					b.parameter_b_begin, b.parameter_b_end,
					b.geometric_tolerance);
		}
	}

	bool build_occt_contact_topology(const TopoDS_Shape& shape,
		std::span<const TopoDS_Face> source_faces,
		std::span<const OcctContactTopologyEdgeFaceInterval> exact_intervals,
		OcctContactTopology& output, std::string& error)
	{
		error.clear();
		OcctContactTopology result;
		try
		{
			if (shape.IsNull() || source_faces.empty())
			{
				error = "contact topology requires a non-null shape and at least one source face";
				return false;
			}

			ShapeIndexedMap shape_edges, shape_faces, shape_vertices;
			TopExp::MapShapes(shape, TopAbs_EDGE, shape_edges);
			TopExp::MapShapes(shape, TopAbs_FACE, shape_faces);
			TopExp::MapShapes(shape, TopAbs_VERTEX, shape_vertices);
			if (shape_edges.Extent() <= 0)
			{
				error = "contact topology shape contains no edges";
				return false;
			}
			if (static_cast<std::uint64_t>(shape_edges.Extent()) >= no_id)
			{
				error = "contact topology exhausted the 32-bit source-edge ID space";
				return false;
			}

			// Preserve the caller's deterministic face IDs, but require each one to be
			// a unique exact face of the supplied parent shape.
			result.faces.reserve(source_faces.size());
			for (std::size_t face_id = 0; face_id < source_faces.size(); ++face_id)
			{
				const TopoDS_Face& face = source_faces[face_id];
				if (face.IsNull() || shape_faces.FindIndex(face) <= 0)
				{
					error = "source face is null or is not part of the supplied shape";
					return false;
				}
				for (std::size_t previous = 0; previous < face_id; ++previous)
					if (same_face(source_faces[previous], face))
					{
						error = "source face list repeats one orientation-insensitive TopoDS face";
						return false;
					}
				result.faces.push_back({static_cast<std::uint32_t>(face_id), face, {}});
			}

			result.edges.reserve(static_cast<std::size_t>(shape_edges.Extent()));
			for (Standard_Integer index = 1; index <= shape_edges.Extent(); ++index)
				result.edges.push_back({static_cast<std::uint32_t>(index - 1),
					TopoDS::Edge(shape_edges(index)), {}});

			std::vector<std::map<std::uint32_t,
				std::vector<const OcctContactTopologyBoundaryOccurrence*>>> occurrences_by_face;
			occurrences_by_face.resize(source_faces.size());
			for (std::size_t face_id = 0; face_id < source_faces.size(); ++face_id)
			{
				auto& face_output = result.faces[face_id];
				std::uint32_t wire_id = 0;
				for (TopExp_Explorer wires(source_faces[face_id], TopAbs_WIRE);
					wires.More(); wires.Next(), ++wire_id)
				{
					for (TopoDS_Iterator child(wires.Current(), true, true); child.More();
						child.Next())
					{
						if (child.Value().ShapeType() != TopAbs_EDGE) continue;
						const TopoDS_Edge occurrence = TopoDS::Edge(child.Value());
						const Standard_Integer edge_index = shape_edges.FindIndex(occurrence);
						if (edge_index <= 0)
						{
							error = "oriented face-wire occurrence is absent from the global edge index";
							return false;
						}
						const auto occurrence_id = static_cast<std::uint32_t>(
							face_output.boundary_occurrences.size());
						face_output.boundary_occurrences.push_back({occurrence_id,
							static_cast<std::uint32_t>(edge_index - 1), wire_id,
							static_cast<std::int8_t>(occurrence.Orientation()), occurrence});
					}
				}
				for (const auto& occurrence : face_output.boundary_occurrences)
					occurrences_by_face[face_id][occurrence.source_edge_id]
						.push_back(&occurrence);
				for (const auto& [edge_id, occurrences] : occurrences_by_face[face_id])
				{
					if (occurrences.size() > 2)
					{
						error = "one source edge has more than two oriented occurrences on one face";
						return false;
					}
					result.edges[edge_id].owner_face_ids.push_back(
						static_cast<std::uint32_t>(face_id));
				}
			}

			std::vector<EdgeCurve> source_curves(result.edges.size());
			std::vector<int> atomizer_index_by_edge(result.edges.size(), -1);
			std::vector<MeshNodeClaimSeed> mesh_node_claim_seeds;
			for (const auto& indexed_edge : result.edges)
			{
				const auto edge_id = indexed_edge.source_edge_id;
				if (indexed_edge.owner_face_ids.empty() || BRep_Tool::Degenerated(indexed_edge.edge))
					continue;
				std::string curve_error;
				if (!make_edge_curve(indexed_edge.edge, source_curves[edge_id], curve_error))
				{
					error = "source edge " + std::to_string(edge_id) + ": " + curve_error;
					return false;
				}

				OcctContactSourceEdge source;
				source.source_edge_id = edge_id;
				source.edge = indexed_edge.edge;
				TopoDS_Vertex edge_first_vertex, edge_last_vertex;
				TopExp::Vertices(indexed_edge.edge, edge_first_vertex, edge_last_vertex, true);
				const double first_vertex_tolerance = edge_first_vertex.IsNull()
					? 0.0 : vertex_tolerance(edge_first_vertex);
				const double last_vertex_tolerance = edge_last_vertex.IsNull()
					? 0.0 : vertex_tolerance(edge_last_vertex);
				const bool closed_topology_edge = !edge_first_vertex.IsNull()
					&& !edge_last_vertex.IsNull()
					&& edge_first_vertex.IsSame(edge_last_vertex);
				for (std::uint32_t face_id : indexed_edge.owner_face_ids)
				{
					const auto found = occurrences_by_face[face_id].find(edge_id);
					if (found == occurrences_by_face[face_id].end()
						|| found->second.empty() || found->second.size() > 2)
					{
						error = "owned source edge lost its oriented face occurrence";
						return false;
					}
					const bool seam = found->second.size() == 2;
					if (seam && !BRep_Tool::IsClosed(indexed_edge.edge,
						source_faces[face_id]))
					{
						error = "two same-face edge occurrences are not an OCCT seam";
						return false;
					}
					OcctContactOwnerFaceUse use;
					use.source_face_id = face_id;
					use.location = seam ? OcctFaceIntervalLocation::interior
						: OcctFaceIntervalLocation::boundary;
					use.sector_count = seam ? 2u : 1u;
					for (const auto* occurrence : found->second)
					{
						if (!occurrence || !valid_orientation(occurrence->orientation))
						{
							error = "owned source edge has an invalid oriented occurrence";
							return false;
						}
						use.boundary_occurrences.push_back({occurrence->occurrence_id,
							edge_id, occurrence->orientation, occurrence->edge,
							source_curves[edge_id].first, source_curves[edge_id].last,
							source_curves[edge_id].first, source_curves[edge_id].last});

						TopLoc_Location location;
						const Handle(Poly_Triangulation) triangulation =
							BRep_Tool::Triangulation(source_faces[face_id], location);
						if (triangulation.IsNull()) continue;
						Handle(Poly_PolygonOnTriangulation) polygon;
						if (!occurrence_polygon_for_pcurve(occurrence->edge,
							source_faces[face_id], triangulation, location, polygon, error))
						{
							error = "cannot select the exact occurrence polygon: " + error;
							return false;
						}
						if (polygon.IsNull()) continue;
						const double scale = std::abs(location.Transformation().ScaleFactor());
						// Face tolerance is used only to certify which original mesh node this
						// exact edge parameter represents.  Once recovered, the retained sample is
						// reevaluated on the source curve and carries no face-wide uncertainty.
						const double face_tolerance = std::max(0.0,
							BRep_Tool::Tolerance(source_faces[face_id])) * scale;
						for (Standard_Integer node = 1; node <= polygon->NbNodes(); ++node)
						{
							const Standard_Integer triangulation_node = polygon->Node(node);
							if (triangulation_node <= 0
								|| triangulation_node > triangulation->NbNodes())
							{
								error = "edge polygon references an invalid face-triangulation node";
								return false;
							}
							gp_Pnt world = triangulation->Node(triangulation_node);
							world.Transform(location.Transformation());
							if (!finite(world.X()) || !finite(world.Y()) || !finite(world.Z()))
							{
								error = "edge polygon references a non-finite face-triangulation node";
								return false;
							}
							const double base_linear_tolerance =
								source_curves[edge_id].geometric_tolerance
								+ face_tolerance + Precision::Confusion();
							double parameter = 0.0;
							// OCCT can store the same native parameter for both copies of a
							// closed polygon endpoint.  The ordered polygon identities are the
							// exact chart branches: retain them as first/last rather than
							// collapsing two original UV nodes onto one atom sample. The polygon
							// node array is already the selected occurrence representation; face
							// orientation must not reverse its endpoints a second time.
							if (closed_topology_edge && node == 1)
								parameter = source_curves[edge_id].first;
							else if (closed_topology_edge && node == polygon->NbNodes())
								parameter = source_curves[edge_id].last;
							else if (polygon->HasParameters()) parameter = polygon->Parameter(node);
							else
							{
								const bool first_node = node == 1;
								const bool last_node = node == polygon->NbNodes();
								// An exact edge-polygon endpoint names the corresponding BRep
								// vertex.  Its triangulated position is therefore certified by
								// that vertex's native tolerance even when it is not exactly on
								// the analytic edge.  Do not apply this endpoint allowance to
								// ordinary interior polygon nodes.
								const double endpoint_search_tolerance = base_linear_tolerance
									+ std::max(first_node ? first_vertex_tolerance : 0.0,
										last_node ? last_vertex_tolerance : 0.0);
								const bool at_first = source_curves[edge_id].curve.Value(
									source_curves[edge_id].first).SquareDistance(world)
									<= endpoint_search_tolerance * endpoint_search_tolerance;
								const bool at_last = source_curves[edge_id].curve.Value(
									source_curves[edge_id].last).SquareDistance(world)
									<= endpoint_search_tolerance * endpoint_search_tolerance;
								if (first_node && at_first)
									parameter = source_curves[edge_id].first;
								else if (last_node && at_last)
									parameter = source_curves[edge_id].last;
								else if (first_node && at_last)
									parameter = source_curves[edge_id].last;
								else if (last_node && at_first)
									parameter = source_curves[edge_id].first;
								else
								{
									std::string projection_error;
									if (!unique_parameter(world, source_curves[edge_id],
										base_linear_tolerance,
										parameter, projection_error))
									{
										error = "cannot recover an owner edge-polygon parameter: "
											+ projection_error;
										return false;
									}
								}
							}
							if (!parameter_inside(source_curves[edge_id], parameter))
							{
								error = "owner edge-polygon parameter escaped the finite source edge";
								return false;
							}
							parameter = std::clamp(parameter, source_curves[edge_id].first,
								source_curves[edge_id].last);
							const gp_Pnt canonical = source_curves[edge_id].curve.Value(parameter);
							double linear_tolerance = base_linear_tolerance;
							if (std::abs(parameter - source_curves[edge_id].first)
								<= source_curves[edge_id].parameter_tolerance)
								linear_tolerance += first_vertex_tolerance;
							if (std::abs(parameter - source_curves[edge_id].last)
								<= source_curves[edge_id].parameter_tolerance)
								linear_tolerance += last_vertex_tolerance;
							if (canonical.SquareDistance(world)
								> linear_tolerance * linear_tolerance)
							{
								error = "owner edge-polygon node misses its exact source curve beyond native CAD bounds";
								return false;
							}
							source.samples.push_back({parameter, 0.0});
							mesh_node_claim_seeds.push_back({static_cast<std::uint32_t>(face_id),
								static_cast<std::uint32_t>(triangulation_node - 1),
								occurrence->occurrence_id, edge_id, parameter,
								{triangulation->UVNode(triangulation_node).X(),
									triangulation->UVNode(triangulation_node).Y()},
								point_array(world), linear_tolerance});
						}
					}
					source.owner_face_uses.push_back(std::move(use));
				}
				std::sort(source.owner_face_uses.begin(), source.owner_face_uses.end(),
					[](const auto& a, const auto& b)
					{
						return a.source_face_id < b.source_face_id;
					});
				std::sort(source.samples.begin(), source.samples.end(),
					[](const auto& a, const auto& b)
					{
						return std::tie(a.parameter, a.tolerance)
							< std::tie(b.parameter, b.tolerance);
					});
				std::vector<OcctContactSourceSample> unique_samples;
				for (const auto& sample : source.samples)
				{
					if (!unique_samples.empty() && close_parameter(source_curves[edge_id],
						unique_samples.back().parameter, sample.parameter))
					{
						unique_samples.back().tolerance = std::max(
							unique_samples.back().tolerance, sample.tolerance);
						continue;
					}
					unique_samples.push_back(sample);
				}
				source.samples = std::move(unique_samples);
				atomizer_index_by_edge[edge_id] = static_cast<int>(
					result.atomizer_source_edges.size());
				result.atomizer_source_edges.push_back(std::move(source));
			}
			if (result.atomizer_source_edges.empty())
			{
				error = "contact topology contains no nondegenerate face-owned source edge";
				return false;
			}

			// Resolve exact interval occurrence IDs against the exact oriented table.
			// Input order is intentionally discarded before reciprocal canonicalization.
			std::vector<OcctContactTopologyEdgeFaceInterval> ordered_intervals(
				exact_intervals.begin(), exact_intervals.end());
			std::sort(ordered_intervals.begin(), ordered_intervals.end(),
				[](const auto& a, const auto& b)
				{
					return std::tie(a.source_edge_id, a.target_face_id,
						a.source_parameter_first, a.source_parameter_last,
						a.interval.begin, a.interval.end, a.interval.location,
						a.interval.target_sector_count,
						a.interval.boundary_occurrence_count,
						a.interval.target_boundary_occurrence_ids,
						a.interval.target_boundary_orientations,
						a.interval.target_parameter_begin,
						a.interval.target_parameter_end,
						a.interval.exact_operation_tolerance)
						< std::tie(b.source_edge_id, b.target_face_id,
							b.source_parameter_first, b.source_parameter_last,
							b.interval.begin, b.interval.end, b.interval.location,
							b.interval.target_sector_count,
							b.interval.boundary_occurrence_count,
							b.interval.target_boundary_occurrence_ids,
							b.interval.target_boundary_orientations,
							b.interval.target_parameter_begin,
							b.interval.target_parameter_end,
							b.interval.exact_operation_tolerance);
				});

			for (const auto& exact : ordered_intervals)
			{
				if (exact.source_edge_id >= result.edges.size()
					|| atomizer_index_by_edge[exact.source_edge_id] < 0
					|| exact.target_face_id >= result.faces.size())
				{
					error = "exact interval references an invalid or non-owned source edge/target face";
					return false;
				}
				const EdgeCurve& source = source_curves[exact.source_edge_id];
				if (!finite(exact.source_parameter_first)
					|| !finite(exact.source_parameter_last)
					|| !close_parameter(source, exact.source_parameter_first, source.first)
					|| !close_parameter(source, exact.source_parameter_last, source.last))
				{
					error = "exact interval raw source range differs from its indexed edge";
					return false;
				}
				const auto& interval = exact.interval;
				if (!finite(interval.begin) || !finite(interval.end)
					|| interval.begin < 0.0 || interval.end > 1.0
					|| !(interval.end > interval.begin)
					|| !finite(interval.exact_operation_tolerance)
					|| interval.exact_operation_tolerance < 0.0
					|| (interval.location == OcctFaceIntervalLocation::boundary
						&& (interval.target_sector_count != 1
							|| interval.boundary_occurrence_count != 1))
					|| (interval.location == OcctFaceIntervalLocation::interior
						&& (interval.target_sector_count != 2
							|| (interval.boundary_occurrence_count != 0
								&& interval.boundary_occurrence_count != 2))))
				{
					error = "exact interval has inconsistent normalized range/sector/occurrence metadata";
					return false;
				}
				const double raw_begin = source.first
					+ (source.last - source.first) * interval.begin;
				const double raw_end = source.first
					+ (source.last - source.first) * interval.end;
				if (!(raw_end > raw_begin + source.parameter_tolerance))
				{
					error = "exact interval has no positive native parameter length";
					return false;
				}

				OcctContactEdgeFaceInterval staged;
				staged.source_edge_id = exact.source_edge_id;
				staged.target_face_id = exact.target_face_id;
				staged.source_parameter_first = source.first;
				staged.source_parameter_last = source.last;
				staged.interval = interval;
				for (std::size_t branch = 0;
					branch < interval.boundary_occurrence_count; ++branch)
				{
					const std::uint32_t occurrence_id =
						interval.target_boundary_occurrence_ids[branch];
					const auto& face = result.faces[exact.target_face_id];
					if (occurrence_id >= face.boundary_occurrences.size())
					{
						error = "exact interval references an unknown target-face occurrence ID";
						return false;
					}
					const auto& occurrence = face.boundary_occurrences[occurrence_id];
					if (occurrence.occurrence_id != occurrence_id
						|| occurrence.orientation
							!= interval.target_boundary_orientations[branch]
						|| !valid_orientation(occurrence.orientation)
						|| occurrence.source_edge_id >= source_curves.size()
						|| atomizer_index_by_edge[occurrence.source_edge_id] < 0)
					{
						error = "exact interval target occurrence does not resolve uniquely to an owned edge";
						return false;
					}
					const EdgeCurve& target = source_curves[occurrence.source_edge_id];
					const double target_begin = interval.target_parameter_begin[branch];
					const double target_end = interval.target_parameter_end[branch];
					const double map_begin = interval.target_mapping_source_begin[branch];
					const double map_end = interval.target_mapping_source_end[branch];
					if (!parameter_inside(target, target_begin)
						|| !parameter_inside(target, target_end)
						|| !finite(map_begin) || !finite(map_end)
						|| std::abs(map_end - map_begin) <= source.parameter_tolerance
						|| !range_contains(map_begin, map_end, raw_begin,
							source.parameter_tolerance)
						|| !range_contains(map_begin, map_end, raw_end,
							source.parameter_tolerance))
					{
						error = "exact interval has an ambiguous or out-of-range target mapping";
						return false;
					}
					const gp_Pnt source_begin_point = source.curve.Value(raw_begin);
					const gp_Pnt source_end_point = source.curve.Value(raw_end);
					const gp_Pnt target_begin_point = target.curve.Value(target_begin);
					const gp_Pnt target_end_point = target.curve.Value(target_end);
					const double geometric_tolerance = source.geometric_tolerance
						+ target.geometric_tolerance + Precision::Confusion();
					if (source_begin_point.SquareDistance(target_begin_point)
							> geometric_tolerance * geometric_tolerance
						|| source_end_point.SquareDistance(target_end_point)
							> geometric_tolerance * geometric_tolerance)
					{
						error = "exact interval endpoint mapping violates native CAD bounds";
						return false;
					}
					staged.target_boundary_occurrences[branch] = occurrence.edge;
					staged.target_boundary_source_edge_ids[branch] =
						occurrence.source_edge_id;

					if (occurrence.source_edge_id != exact.source_edge_id)
					{
						OcctContactReciprocalInterval reciprocal;
						reciprocal.source_edge_a = exact.source_edge_id;
						reciprocal.source_edge_b = occurrence.source_edge_id;
						reciprocal.parameter_a_begin = raw_begin;
						reciprocal.parameter_a_end = raw_end;
						reciprocal.parameter_b_begin = target_begin;
						reciprocal.parameter_b_end = target_end;
						reciprocal.geometric_tolerance = 0.0;
						canonicalize_reciprocal(reciprocal);
						result.reciprocal_intervals.push_back(reciprocal);
					}
				}

				bool duplicate = false;
				if (!result.atomizer_edge_face_intervals.empty())
				{
					auto& previous = result.atomizer_edge_face_intervals.back();
					if (previous.source_edge_id == staged.source_edge_id
						&& previous.target_face_id == staged.target_face_id)
					{
						const double previous_begin = source.first
							+ (source.last - source.first) * previous.interval.begin;
						const double previous_end = source.first
							+ (source.last - source.first) * previous.interval.end;
						const bool same_range = std::abs(previous_begin - raw_begin)
								<= source.parameter_tolerance
							&& std::abs(previous_end - raw_end)
								<= source.parameter_tolerance;
						const bool positive_overlap = std::min(previous_end, raw_end)
							> std::max(previous_begin, raw_begin) + source.parameter_tolerance;
						if (same_range)
						{
							const auto& a = previous.interval;
							const auto& b = staged.interval;
							duplicate = a.location == b.location
								&& a.target_sector_count == b.target_sector_count
								&& a.boundary_occurrence_count == b.boundary_occurrence_count
								&& a.target_boundary_occurrence_ids
									== b.target_boundary_occurrence_ids
								&& a.target_boundary_orientations
									== b.target_boundary_orientations;
							if (duplicate)
								for (std::size_t branch = 0;
									branch < b.boundary_occurrence_count; ++branch)
								{
									const auto occurrence_id =
										b.target_boundary_occurrence_ids[branch];
									const auto& occurrence = result.faces[staged.target_face_id]
										.boundary_occurrences[occurrence_id];
									const EdgeCurve& target =
										source_curves[occurrence.source_edge_id];
									duplicate = duplicate
										&& std::abs(a.target_parameter_begin[branch]
											- b.target_parameter_begin[branch])
											<= target.parameter_tolerance
										&& std::abs(a.target_parameter_end[branch]
											- b.target_parameter_end[branch])
											<= target.parameter_tolerance
										&& std::abs(a.target_mapping_source_begin[branch]
											- b.target_mapping_source_begin[branch])
											<= source.parameter_tolerance
										&& std::abs(a.target_mapping_source_end[branch]
											- b.target_mapping_source_end[branch])
											<= source.parameter_tolerance;
								}
							if (!duplicate)
							{
								error = "one exact source/target span has conflicting occurrence metadata";
								return false;
							}
							previous.interval.exact_operation_tolerance = std::max(
								previous.interval.exact_operation_tolerance,
								staged.interval.exact_operation_tolerance);
						}
						else if (positive_overlap)
						{
							error = "exact intervals overlap ambiguously on one source edge/target face";
							return false;
						}
					}
				}
				if (!duplicate) result.atomizer_edge_face_intervals.push_back(std::move(staged));
			}

			// A target face can already have a tessellation vertex in the interior of an
			// exact contact curve even though that vertex is absent from the source
			// edge's own polygon.  The face conformer must not later remove such a vertex
			// or create a hanging node.  Harvest those vertices globally, before
			// atomization, so the atomizer's union chain propagates the split to every
			// reciprocal edge and every incident face use.
			//
			// This is not contact discovery.  Candidate source-edge/target-face pairs and
			// their finite parameter spans come exclusively from the exact certificates
			// staged above.  Tessellated positions only select a unique parameter on that
			// bounded span; the eventual canonical point is reevaluated on the source CAD
			// curve by the atomizer.
			std::vector<std::vector<std::size_t>> intervals_by_target_face(
				result.faces.size());
			for (std::size_t interval_id = 0;
				interval_id < result.atomizer_edge_face_intervals.size(); ++interval_id)
			{
				const auto& staged = result.atomizer_edge_face_intervals[interval_id];
				// Boundary and seam uses already inherit the exact node identities from
				// their target edge polygons.  The full-face scan is needed only when the
				// certified contact lies strictly inside the target face.
				if (staged.interval.location != OcctFaceIntervalLocation::interior
					|| staged.interval.boundary_occurrence_count != 0)
					continue;
				const auto face_id = staged.target_face_id;
				if (face_id >= intervals_by_target_face.size())
				{
					error = "staged exact interval lost its target face";
					return false;
				}
				intervals_by_target_face[face_id].push_back(interval_id);
			}
			for (std::size_t face_id = 0;
				face_id < intervals_by_target_face.size(); ++face_id)
			{
				if (intervals_by_target_face[face_id].empty()) continue;
				TopLoc_Location location;
				const Handle(Poly_Triangulation) triangulation =
					BRep_Tool::Triangulation(source_faces[face_id], location);
				if (triangulation.IsNull()) continue;
				const double scale = std::abs(location.Transformation().ScaleFactor());
				const double target_face_tolerance = std::max(0.0,
					BRep_Tool::Tolerance(source_faces[face_id])) * scale;
				std::vector<std::pair<std::uint32_t, gp_Pnt>> target_nodes;
				target_nodes.reserve(static_cast<std::size_t>(triangulation->NbNodes()));
				for (Standard_Integer node = 1; node <= triangulation->NbNodes(); ++node)
				{
					gp_Pnt world = triangulation->Node(node);
					world.Transform(location.Transformation());
					if (finite(world.X()) && finite(world.Y()) && finite(world.Z()))
						target_nodes.emplace_back(static_cast<std::uint32_t>(node - 1), world);
				}

				for (std::size_t interval_id : intervals_by_target_face[face_id])
				{
					const auto& staged = result.atomizer_edge_face_intervals[interval_id];
					const EdgeCurve& source = source_curves[staged.source_edge_id];
					const double raw_begin = source.first
						+ (source.last - source.first) * staged.interval.begin;
					const double raw_end = source.first
						+ (source.last - source.first) * staged.interval.end;
					const double operation_tolerance =
						staged.interval.exact_operation_tolerance;
					const double linear_tolerance = source.geometric_tolerance
						+ target_face_tolerance + operation_tolerance;
					const int source_index = atomizer_index_by_edge[staged.source_edge_id];
					if (source_index < 0
						|| static_cast<std::size_t>(source_index)
							>= result.atomizer_source_edges.size()
						|| !finite(linear_tolerance) || linear_tolerance < 0.0)
					{
						error = "staged exact interval cannot enrich its source sample chain";
						return false;
					}
					auto& samples = result.atomizer_source_edges[
						static_cast<std::size_t>(source_index)].samples;
					for (const auto& [target_node_index, world] : target_nodes)
					{
						double parameter = 0.0;
						if (unique_parameter_on_span(world, source, raw_begin, raw_end,
							linear_tolerance, parameter) != BoundedProjection::unique)
							continue;
						if (!range_contains(raw_begin, raw_end, parameter,
							source.parameter_tolerance))
							continue;
						// The target-face/operation bound selected this sample for this
						// already-certified face use; it is not uncertainty in the canonical
						// source-curve point.  Keep it isolated on OcctContactAtomFaceUse.
						// Publishing it as a source-sample tolerance would incorrectly widen
						// every unrelated face use of the resulting atom.
						parameter = std::clamp(parameter, raw_begin, raw_end);
						samples.push_back({parameter, 0.0});
						mesh_node_claim_seeds.push_back({static_cast<std::uint32_t>(face_id),
							target_node_index, OcctContactBoundaryOccurrence::no_occurrence_id,
							staged.source_edge_id, parameter, {}, point_array(world),
							linear_tolerance});
					}
				}
			}

			// Deterministically union target-derived samples with all owner-polygon
			// samples.  A source parameter is the identity; tessellated 3-D positions are
			// deliberately not retained.
			for (auto& source : result.atomizer_source_edges)
			{
				const EdgeCurve& curve = source_curves[source.source_edge_id];
				std::sort(source.samples.begin(), source.samples.end(),
					[](const auto& a, const auto& b)
					{
						return std::tie(a.parameter, a.tolerance)
							< std::tie(b.parameter, b.tolerance);
					});
				std::vector<OcctContactSourceSample> unique_samples;
				for (const auto& sample : source.samples)
				{
					if (!unique_samples.empty() && close_parameter(curve,
						unique_samples.back().parameter, sample.parameter))
					{
						unique_samples.back().tolerance = std::max(
							unique_samples.back().tolerance, sample.tolerance);
						continue;
					}
					unique_samples.push_back(sample);
				}
				source.samples = std::move(unique_samples);
			}

			std::sort(result.reciprocal_intervals.begin(),
				result.reciprocal_intervals.end(), reciprocal_less);
			std::vector<OcctContactReciprocalInterval> unique_reciprocals;
			for (const auto& reciprocal : result.reciprocal_intervals)
			{
				bool duplicate = false;
				if (!unique_reciprocals.empty())
				{
					const auto& previous = unique_reciprocals.back();
					if (previous.source_edge_a == reciprocal.source_edge_a
						&& previous.source_edge_b == reciprocal.source_edge_b)
					{
						const auto& a = source_curves[reciprocal.source_edge_a];
						const auto& b = source_curves[reciprocal.source_edge_b];
						duplicate = std::abs(previous.parameter_a_begin
								- reciprocal.parameter_a_begin) <= a.parameter_tolerance
							&& std::abs(previous.parameter_a_end
								- reciprocal.parameter_a_end) <= a.parameter_tolerance
							&& std::abs(previous.parameter_b_begin
								- reciprocal.parameter_b_begin) <= b.parameter_tolerance
							&& std::abs(previous.parameter_b_end
								- reciprocal.parameter_b_end) <= b.parameter_tolerance;
					}
				}
				if (!duplicate) unique_reciprocals.push_back(reciprocal);
			}
			result.reciprocal_intervals = std::move(unique_reciprocals);

			// Exact topological vertices are authoritative even if their incident edge
			// parameterizations use distinct endpoint values (closed edge first/last).
			std::map<std::uint32_t, JunctionCandidate> topology_junctions;
			for (const auto& source : result.atomizer_source_edges)
			{
				const EdgeCurve& edge = source_curves[source.source_edge_id];
				TopoDS_Vertex first_vertex, last_vertex;
				TopExp::Vertices(source.edge, first_vertex, last_vertex, true);
				if (first_vertex.IsNull() || last_vertex.IsNull())
				{
					error = "nondegenerate face-owned edge has no two endpoint occurrences";
					return false;
				}
				const Standard_Integer first_id = shape_vertices.FindIndex(first_vertex);
				const Standard_Integer last_id = shape_vertices.FindIndex(last_vertex);
				if (first_id <= 0 || last_id <= 0)
				{
					error = "source edge endpoint is absent from the global topology-vertex index";
					return false;
				}
				auto append_vertex = [&](const TopoDS_Vertex& vertex,
					std::uint32_t vertex_id, double parameter)
				{
					JunctionCandidate& candidate = topology_junctions[vertex_id];
					candidate.kind = 0;
					candidate.stable_key = vertex_id;
					candidate.position = point_array(BRep_Tool::Pnt(vertex));
					candidate.tolerance = std::max(candidate.tolerance,
						vertex_tolerance(vertex));
					candidate.incidences.push_back({source.source_edge_id, parameter});
				};

				// TopExp reports oriented first/last vertices. Resolve their actual native
				// curve endpoint parameters by exact endpoint distance, including reversal.
				const gp_Pnt curve_first = edge.curve.Value(edge.first);
				const gp_Pnt curve_last = edge.curve.Value(edge.last);
				const gp_Pnt vertex_first = BRep_Tool::Pnt(first_vertex);
				const gp_Pnt vertex_last = BRep_Tool::Pnt(last_vertex);
				if (first_vertex.IsSame(last_vertex))
				{
					append_vertex(first_vertex, static_cast<std::uint32_t>(first_id - 1),
						edge.first);
					append_vertex(last_vertex, static_cast<std::uint32_t>(last_id - 1),
						edge.last);
				}
				else
				{
					const double forward = curve_first.SquareDistance(vertex_first)
						+ curve_last.SquareDistance(vertex_last);
					const double reverse = curve_first.SquareDistance(vertex_last)
						+ curve_last.SquareDistance(vertex_first);
					const double tolerance = edge.geometric_tolerance
						+ vertex_tolerance(first_vertex) + vertex_tolerance(last_vertex)
						+ Precision::Confusion();
					if (std::min(forward, reverse) > 2.0 * tolerance * tolerance)
					{
						error = "topology vertices do not bound their edge within native CAD tolerances";
						return false;
					}
					if (forward <= reverse)
					{
						append_vertex(first_vertex, static_cast<std::uint32_t>(first_id - 1),
							edge.first);
						append_vertex(last_vertex, static_cast<std::uint32_t>(last_id - 1),
							edge.last);
					}
					else
					{
						append_vertex(last_vertex, static_cast<std::uint32_t>(last_id - 1),
							edge.first);
						append_vertex(first_vertex, static_cast<std::uint32_t>(first_id - 1),
							edge.last);
					}
				}
			}

			std::vector<JunctionCandidate> junctions;
			for (auto& [unused, candidate] : topology_junctions)
			{
				(void)unused;
				std::sort(candidate.incidences.begin(), candidate.incidences.end(),
					incidence_less);
				candidate.incidences.erase(std::unique(candidate.incidences.begin(),
					candidate.incidences.end(), [&](const auto& a, const auto& b)
					{
						return a.source_edge_id == b.source_edge_id
							&& close_parameter(source_curves[a.source_edge_id],
								a.parameter, b.parameter);
					}), candidate.incidences.end());
				if (candidate.incidences.size() >= 2) junctions.push_back(std::move(candidate));
			}

			// One original face-triangulation node can be named by several certified
			// contact chains (for example, where two trimmed boundary edges meet).  That
			// exact face-local node identity is stronger than a UV proximity test and must
			// union the chains before topology IDs are published.  It does not merge nearby
			// nodes: every incidence below independently passed its edge/face certificate.
			std::map<std::pair<std::uint32_t, std::uint32_t>,
				std::vector<const MeshNodeClaimSeed*>> seeds_by_face_node;
			for (const auto& seed : mesh_node_claim_seeds)
				seeds_by_face_node[{seed.source_face_id, seed.source_node_index}]
					.push_back(&seed);
			for (const auto& [face_node, seeds] : seeds_by_face_node)
			{
				if (seeds.size() < 2) continue;
				JunctionCandidate candidate;
				candidate.kind = 1;
				candidate.stable_key = (static_cast<std::uint64_t>(face_node.first) << 32)
					| face_node.second;
				candidate.position = seeds.front()->world_position;
				for (const MeshNodeClaimSeed* seed : seeds)
				{
					if (!seed || seed->source_edge_id >= source_curves.size()
						|| point(seed->world_position).SquareDistance(point(candidate.position)) != 0.0)
					{
						error = "one exact face-mesh node has inconsistent stored provenance";
						return false;
					}
					candidate.tolerance = std::max(candidate.tolerance,
						seed->position_tolerance);
					candidate.incidences.push_back({seed->source_edge_id,
						seed->source_parameter});
				}
				std::sort(candidate.incidences.begin(), candidate.incidences.end(),
					incidence_less);
				candidate.incidences.erase(std::unique(candidate.incidences.begin(),
					candidate.incidences.end(), [&](const auto& a, const auto& b)
					{
						return a.source_edge_id == b.source_edge_id
							&& close_parameter(source_curves[a.source_edge_id],
								a.parameter, b.parameter);
					}), candidate.incidences.end());
				if (candidate.incidences.size() < 2) continue;
				const bool already_unioned = std::any_of(junctions.begin(), junctions.end(),
					[&](const JunctionCandidate& existing)
					{
						return std::all_of(candidate.incidences.begin(),
							candidate.incidences.end(), [&](const auto& incidence)
							{
								return std::any_of(existing.incidences.begin(),
									existing.incidences.end(), [&](const auto& known)
									{
										return incidence.source_edge_id == known.source_edge_id
											&& close_parameter(source_curves[incidence.source_edge_id],
												incidence.parameter, known.parameter);
									});
							});
					});
				if (!already_unioned) junctions.push_back(std::move(candidate));
			}

			// Index exact shared TopoDS vertices by every pair of their incident source
			// edges.  IntTools can additionally report the mathematical curve/curve
			// intersection a tiny distance from such a vertex when the BRep curves meet
			// only within that vertex's native tolerance.  That numerical hit represents
			// the already-declared topological vertex; emitting both creates an
			// unresolvable microscopic atom.  The index provides exact topology identity;
			// distance is used below only to verify that an IntTools hit lies inside the
			// identified vertex's own CAD bound, never to discover or heal contact.
			std::map<EdgePair, std::vector<std::size_t>> shared_vertex_junctions;
			for (std::size_t junction = 0; junction < junctions.size(); ++junction)
			{
				std::vector<std::uint32_t> incident_edges;
				for (const auto& incidence : junctions[junction].incidences)
					incident_edges.push_back(incidence.source_edge_id);
				std::sort(incident_edges.begin(), incident_edges.end());
				incident_edges.erase(std::unique(incident_edges.begin(),
					incident_edges.end()), incident_edges.end());
				for (std::size_t a = 0; a < incident_edges.size(); ++a)
					for (std::size_t b = a + 1; b < incident_edges.size(); ++b)
						shared_vertex_junctions[{incident_edges[a], incident_edges[b]}]
							.push_back(junction);
			}

			std::set<EdgePair> reciprocal_pairs;
			for (const auto& reciprocal : result.reciprocal_intervals)
				reciprocal_pairs.insert(ordered_pair(reciprocal.source_edge_a,
					reciprocal.source_edge_b));

			// Sweep-and-prune is only a candidate generator for zero-length exact
			// edge/edge junctions. Positive overlaps are accepted exclusively through
			// the caller-supplied edge/face interval certificates above.
			std::vector<EdgeBounds> bounds;
			bounds.reserve(result.atomizer_source_edges.size());
			for (std::size_t source_index = 0;
				source_index < result.atomizer_source_edges.size(); ++source_index)
			{
				EdgeBounds bound;
				std::string bounds_error;
				if (!make_bounds(result.atomizer_source_edges[source_index].edge,
					source_index, bound, bounds_error))
				{
					error = "cannot index an edge for exact junctions: " + bounds_error;
					return false;
				}
				bounds.push_back(bound);
			}
			std::sort(bounds.begin(), bounds.end(), [&](const auto& a, const auto& b)
			{
				const auto edge_a = result.atomizer_source_edges[a.source_index].source_edge_id;
				const auto edge_b = result.atomizer_source_edges[b.source_index].source_edge_id;
				return std::tie(a.low[0], a.high[0], edge_a)
					< std::tie(b.low[0], b.high[0], edge_b);
			});
			std::vector<std::size_t> active;
			std::uint64_t intersection_key = 0;
			for (std::size_t sorted_index = 0; sorted_index < bounds.size(); ++sorted_index)
			{
				const EdgeBounds& current = bounds[sorted_index];
				active.erase(std::remove_if(active.begin(), active.end(),
					[&](std::size_t previous)
					{
						return bounds[previous].high[0] < current.low[0];
					}), active.end());
				for (std::size_t previous : active)
				{
					const EdgeBounds& other = bounds[previous];
					if (!overlaps_yz(current, other)) continue;
					const auto& source_a = result.atomizer_source_edges[current.source_index];
					const auto& source_b = result.atomizer_source_edges[other.source_index];
					const EdgePair pair = ordered_pair(source_a.source_edge_id,
						source_b.source_edge_id);
					if (reciprocal_pairs.contains(pair)) continue;

					IntTools_EdgeEdge intersection(source_a.edge, source_b.edge);
					intersection.SetFuzzyValue(0.0);
					intersection.Perform();
					if (!intersection.IsDone())
					{
						error = "zero-fuzzy edge/edge junction intersection did not complete";
						return false;
					}
					const EdgeCurve& edge_a = source_curves[source_a.source_edge_id];
					const EdgeCurve& edge_b = source_curves[source_b.source_edge_id];
					bool positive_overlap = false;
					std::vector<std::pair<double, double>> zero_edge_parameters;
					for (NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
						intersection.CommonParts()); common.More(); common.Next())
					{
						if (common.Value().Type() != TopAbs_EDGE) continue;
						const bool first_is_a = common.Value().Edge1().IsSame(source_a.edge)
							&& common.Value().Edge2().IsSame(source_b.edge);
						const bool first_is_b = common.Value().Edge1().IsSame(source_b.edge)
							&& common.Value().Edge2().IsSame(source_a.edge);
						if (!first_is_a && !first_is_b)
						{
							error = "zero-fuzzy edge overlap did not retain its exact input-edge identity";
							return false;
						}
						double first_begin = 0.0, first_end = 0.0;
						common.Value().Range1(first_begin, first_end);
						bool has_second_range = false;
						for (NCollection_Sequence<IntTools_Range>::Iterator second(
							common.Value().Ranges2()); second.More(); second.Next())
						{
							has_second_range = true;
							const double second_begin = second.Value().First();
							const double second_end = second.Value().Last();
							const double a_begin = first_is_a ? first_begin : second_begin;
							const double a_end = first_is_a ? first_end : second_end;
							const double b_begin = first_is_a ? second_begin : first_begin;
							const double b_end = first_is_a ? second_end : first_end;
							if (std::abs(a_end - a_begin) > edge_a.parameter_tolerance
								&& std::abs(b_end - b_begin) > edge_b.parameter_tolerance)
								positive_overlap = true;
							else
								zero_edge_parameters.emplace_back(0.5 * (a_begin + a_end),
									0.5 * (b_begin + b_end));
						}
						const double first_tolerance = first_is_a
							? edge_a.parameter_tolerance : edge_b.parameter_tolerance;
						if (!has_second_range
							&& std::abs(first_end - first_begin) > first_tolerance)
							positive_overlap = true;
					}
					if (positive_overlap)
					{
						error = "positive-length overlap between source edges "
							+ std::to_string(pair.first) + " and "
							+ std::to_string(pair.second)
							+ " has no supplied exact edge/face certificate";
						return false;
					}
					auto append_junction = [&](double parameter_a,
						double parameter_b) -> bool
					{
						if (!parameter_inside(edge_a, parameter_a)
							|| !parameter_inside(edge_b, parameter_b))
						{
							error = "zero-fuzzy edge junction returned an out-of-range parameter";
							return false;
						}
						const gp_Pnt point_a = edge_a.curve.Value(parameter_a);
						const gp_Pnt point_b = edge_b.curve.Value(parameter_b);
						const double tolerance = edge_a.geometric_tolerance
							+ edge_b.geometric_tolerance + Precision::Confusion();
						if (point_a.SquareDistance(point_b) > tolerance * tolerance)
						{
							std::ostringstream report;
							report.precision(17);
							report << "zero-fuzzy edge junction between source edges "
								<< source_a.source_edge_id << " and "
								<< source_b.source_edge_id
								<< " violates native CAD bounds: parameters "
								<< parameter_a << ", " << parameter_b
								<< "; distance=" << point_a.Distance(point_b)
								<< "; bound=" << tolerance;
							error = report.str();
							return false;
						}
						const auto shared_vertices = shared_vertex_junctions.find(pair);
						if (shared_vertices != shared_vertex_junctions.end())
							for (std::size_t junction_id : shared_vertices->second)
							{
								const JunctionCandidate& vertex = junctions[junction_id];
								const gp_Pnt vertex_point = point(vertex.position);
								const double vertex_bound = vertex.tolerance;
								if (vertex_point.SquareDistance(point_a)
										<= vertex_bound * vertex_bound
									&& vertex_point.SquareDistance(point_b)
										<= vertex_bound * vertex_bound)
									return true;
							}
						JunctionCandidate candidate;
						candidate.kind = 2;
						candidate.stable_key = intersection_key++;
						candidate.position = point_array(source_a.source_edge_id
							< source_b.source_edge_id ? point_a : point_b);
						candidate.tolerance = tolerance;
						candidate.incidences = {{source_a.source_edge_id, parameter_a},
							{source_b.source_edge_id, parameter_b}};
						std::sort(candidate.incidences.begin(), candidate.incidences.end(),
							incidence_less);
						junctions.push_back(std::move(candidate));
						return true;
					};
					for (NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
						intersection.CommonParts()); common.More(); common.Next())
					{
						if (common.Value().Type() != TopAbs_VERTEX) continue;
						const bool first_is_a = common.Value().Edge1().IsSame(source_a.edge)
							&& common.Value().Edge2().IsSame(source_b.edge);
						const bool first_is_b = common.Value().Edge1().IsSame(source_b.edge)
							&& common.Value().Edge2().IsSame(source_a.edge);
						if (!first_is_a && !first_is_b)
						{
							error = "zero-fuzzy edge junction did not retain its exact input-edge identity";
							return false;
						}
						const double first_parameter = common.Value().VertexParameter1();
						const double second_parameter = common.Value().VertexParameter2();
						const double parameter_a = first_is_a
							? first_parameter : second_parameter;
						const double parameter_b = first_is_a
							? second_parameter : first_parameter;
						if (!append_junction(parameter_a, parameter_b)) return false;
					}
					for (const auto [parameter_a, parameter_b] : zero_edge_parameters)
						if (!append_junction(parameter_a, parameter_b)) return false;
				}
				active.push_back(sorted_index);
			}

			std::sort(junctions.begin(), junctions.end(), [](const auto& a, const auto& b)
			{
				if (!incidences_equal_exact(a.incidences, b.incidences))
					return std::lexicographical_compare(a.incidences.begin(), a.incidences.end(),
						b.incidences.begin(), b.incidences.end(), incidence_less);
				return std::tie(a.kind, a.stable_key, a.position)
					< std::tie(b.kind, b.stable_key, b.position);
			});
			std::vector<JunctionCandidate> unique_junctions;
			for (auto& candidate : junctions)
			{
				bool duplicate = false;
				if (!unique_junctions.empty()
					&& unique_junctions.back().incidences.size()
						== candidate.incidences.size())
				{
					const auto& previous = unique_junctions.back();
					duplicate = true;
					for (std::size_t incidence = 0;
						incidence < candidate.incidences.size(); ++incidence)
					{
						const auto& a = previous.incidences[incidence];
						const auto& b = candidate.incidences[incidence];
						if (a.source_edge_id != b.source_edge_id
							|| !close_parameter(source_curves[a.source_edge_id],
								a.parameter, b.parameter))
						{
							duplicate = false;
							break;
						}
					}
					if (duplicate)
					{
						const double tolerance = previous.tolerance + candidate.tolerance
							+ Precision::Confusion();
						if (point(previous.position).SquareDistance(point(candidate.position))
							> tolerance * tolerance)
						{
							error = "duplicate exact junctions disagree beyond native CAD bounds";
							return false;
						}
						unique_junctions.back().tolerance = std::max(previous.tolerance,
							candidate.tolerance);
					}
				}
				if (!duplicate) unique_junctions.push_back(std::move(candidate));
			}
			if (unique_junctions.size() >= std::numeric_limits<std::uint64_t>::max())
			{
				error = "contact topology exhausted the exact-junction ID space";
				return false;
			}
			for (std::size_t junction_id = 0; junction_id < unique_junctions.size(); ++junction_id)
			{
				auto& candidate = unique_junctions[junction_id];
				result.exact_junctions.push_back({static_cast<std::uint64_t>(junction_id),
					std::move(candidate.incidences), candidate.position, candidate.tolerance});
			}

			const OcctContactAtomizerInput atomizer_input{
				result.atomizer_source_edges,
				result.atomizer_edge_face_intervals,
				result.reciprocal_intervals,
				result.exact_junctions};
			std::string atomizer_error;
			if (!atomize_occt_contacts(atomizer_input, result.atomization, atomizer_error))
			{
				error = "contact atomization failed: " + atomizer_error;
				return false;
			}

			// Resolve each original face-mesh node through the atomizer's final union
			// topology.  Reciprocal copies and atom-boundary duplicates must all name one
			// topology ID; coordinate coincidence alone is never used here.
			struct ParameterTopology
			{
				double parameter = 0.0;
				std::uint64_t atom_id = std::numeric_limits<std::uint64_t>::max();
				std::uint32_t sample_index = no_id;
				std::uint32_t topology_id = no_id;
			};
			std::vector<std::vector<ParameterTopology>> topology_by_edge(result.edges.size());
			for (const auto& atom : result.atomization.atoms)
				for (std::size_t sample_index = 0; sample_index < atom.samples.size(); ++sample_index)
				{
					const auto& sample = atom.samples[sample_index];
					for (const auto& source : sample.source_locations)
					{
						if (source.source_edge_id >= topology_by_edge.size()
							|| sample.topology_id == no_id || sample_index >= no_id)
						{
							error = "atomized contact sample has invalid source/topology provenance";
							return false;
						}
						topology_by_edge[source.source_edge_id].push_back(
							{source.parameter, atom.id,
								static_cast<std::uint32_t>(sample_index), sample.topology_id});
					}
				}
			for (auto& samples : topology_by_edge)
			{
				std::sort(samples.begin(), samples.end(), [](const auto& a, const auto& b)
				{
					return std::tie(a.parameter, a.atom_id, a.sample_index, a.topology_id)
						< std::tie(b.parameter, b.atom_id, b.sample_index, b.topology_id);
				});
			}

			for (const auto& seed : mesh_node_claim_seeds)
			{
				if (seed.source_face_id >= result.faces.size()
					|| seed.source_edge_id >= topology_by_edge.size()
					|| seed.source_node_index == no_id
					|| !finite(seed.source_parameter))
				{
					error = "original face-mesh node claim has invalid source provenance";
					return false;
				}
				const auto& curve = source_curves[seed.source_edge_id];
				const auto& candidates = topology_by_edge[seed.source_edge_id];
				const double low = seed.source_parameter - curve.parameter_tolerance;
				const double high = seed.source_parameter + curve.parameter_tolerance;
				auto candidate = std::lower_bound(candidates.begin(), candidates.end(), low,
					[](const ParameterTopology& value, double parameter)
					{
						return value.parameter < parameter;
					});
				std::optional<std::uint32_t> topology_id;
				std::vector<ParameterTopology> matches;
				for (; candidate != candidates.end() && candidate->parameter <= high; ++candidate)
				{
					if (!topology_id) topology_id = candidate->topology_id;
					else if (*topology_id != candidate->topology_id)
					{
						error = "one original face-mesh node sample resolves to conflicting contact topology IDs";
						return false;
					}
					matches.push_back(*candidate);
				}
				if (!topology_id)
				{
					error = "original face-mesh node sample was lost during contact atomization";
					return false;
				}
				std::sort(matches.begin(), matches.end(), [](const auto& a, const auto& b)
				{
					return std::tie(a.atom_id, a.sample_index, a.topology_id)
						< std::tie(b.atom_id, b.sample_index, b.topology_id);
				});
				matches.erase(std::unique(matches.begin(), matches.end(),
					[](const auto& a, const auto& b)
					{
						return std::tie(a.atom_id, a.sample_index, a.topology_id)
							== std::tie(b.atom_id, b.sample_index, b.topology_id);
					}), matches.end());

				const bool boundary_seed = seed.boundary_occurrence_id
					!= OcctContactBoundaryOccurrence::no_occurrence_id;
				const OcctContactTopologyFace* seed_face = nullptr;
				const OcctContactTopologyBoundaryOccurrence* seed_occurrence = nullptr;
				Handle(Geom2d_Curve) seed_pcurve;
				double seed_pcurve_first = 0.0, seed_pcurve_last = 0.0;
				double seed_uv_tolerance = 0.0;
				if (boundary_seed)
				{
					if (seed.source_face_id >= result.faces.size()
						|| seed.boundary_occurrence_id
							>= result.faces[seed.source_face_id].boundary_occurrences.size())
					{
						error = "boundary mesh-node seed lost its exact face occurrence";
						return false;
					}
					seed_face = &result.faces[seed.source_face_id];
					seed_occurrence = &seed_face->boundary_occurrences[
						seed.boundary_occurrence_id];
					seed_pcurve = BRep_Tool::CurveOnSurface(seed_occurrence->edge,
						seed_face->face, seed_pcurve_first, seed_pcurve_last);
					seed_uv_tolerance = face_uv_identity_tolerance(seed_face->face,
						seed.face_uv, seed.position_tolerance);
					if (seed_pcurve.IsNull() || !finite(seed_pcurve_first)
						|| !finite(seed_pcurve_last) || !(seed_pcurve_last > seed_pcurve_first)
						|| !finite(seed_uv_tolerance) || !(seed_uv_tolerance > 0.0))
					{
						error = "boundary mesh-node seed has no resolved exact face chart";
						return false;
					}
				}

				bool emitted_claim = false;
				for (std::size_t begin = 0; begin < matches.size();)
				{
					std::size_t end = begin + 1;
					while (end < matches.size()
						&& matches[end].atom_id == matches[begin].atom_id) ++end;
					const ParameterTopology* selected = &matches[begin];
					if (boundary_seed)
					{
						if (matches[begin].atom_id >= result.atomization.atoms.size())
						{
							error = "mesh-node claim candidate references an unknown contact atom";
							return false;
						}
						const auto& atom = result.atomization.atoms[
							static_cast<std::size_t>(matches[begin].atom_id)];
						const auto face_use = std::find_if(atom.face_uses.begin(),
							atom.face_uses.end(), [&](const auto& use)
							{
								return use.source_face_id == seed.source_face_id;
							});
						if (face_use == atom.face_uses.end())
						{
							begin = end;
							continue;
						}
						const auto boundary_use = std::find_if(
							face_use->boundary_occurrences.begin(),
							face_use->boundary_occurrences.end(), [&](const auto& use)
							{
								return use.occurrence_id == seed.boundary_occurrence_id;
							});
						if (boundary_use == face_use->boundary_occurrences.end())
						{
							begin = end;
							continue;
						}
						selected = nullptr;
						std::ostringstream candidate_report;
						candidate_report.precision(17);
						for (std::size_t candidate_index = begin;
							candidate_index < end; ++candidate_index)
						{
							const auto& candidate_match = matches[candidate_index];
							if (candidate_match.sample_index
								>= boundary_use->sample_target_parameters.size())
							{
								error = "mesh-node claim candidate lost its occurrence parameter";
								return false;
							}
							const double target_parameter =
								boundary_use->sample_target_parameters[
									candidate_match.sample_index];
							if (!finite(target_parameter)
								|| target_parameter < seed_pcurve_first
								|| target_parameter > seed_pcurve_last)
							{
								error = "mesh-node claim candidate escaped its occurrence pcurve";
								return false;
							}
							const gp_Pnt2d candidate_uv =
								seed_pcurve->Value(target_parameter);
							const double uv_distance = std::hypot(
								candidate_uv.X() - seed.face_uv[0],
								candidate_uv.Y() - seed.face_uv[1]);
							candidate_report << " sample " << candidate_match.sample_index
								<< " uv [" << candidate_uv.X() << ',' << candidate_uv.Y()
								<< "] distance " << uv_distance;
							if (!finite(uv_distance) || uv_distance > seed_uv_tolerance)
								continue;
							if (selected)
							{
								error = "one exact occurrence UV names multiple samples in one contact atom";
								return false;
							}
							selected = &candidate_match;
						}
						if (!selected)
						{
							std::ostringstream message;
							message.precision(17);
							message << "no contact-atom sample matches original face "
								<< seed.source_face_id << " occurrence "
								<< seed.boundary_occurrence_id << " node "
								<< seed.source_node_index << " UV [" << seed.face_uv[0]
								<< ',' << seed.face_uv[1] << "] in atom "
								<< matches[begin].atom_id << " (limit "
								<< seed_uv_tolerance << "):" << candidate_report.str();
							error = message.str();
							return false;
						}
					}
					else if (end - begin > 1)
					{
						error = "one interior face-mesh node sample resolves ambiguously within one contact atom";
						return false;
					}
					result.mesh_node_claims.push_back({seed.source_face_id,
						seed.source_node_index, seed.boundary_occurrence_id,
						selected->atom_id, selected->sample_index,
						selected->topology_id});
					emitted_claim = true;
					begin = end;
				}
				if (!emitted_claim)
				{
					error = "original face-mesh node has no atom on its exact face occurrence";
					return false;
				}
			}
			std::sort(result.mesh_node_claims.begin(), result.mesh_node_claims.end(),
				[](const auto& a, const auto& b)
				{
					return std::tie(a.source_face_id, a.source_node_index,
						a.boundary_occurrence_id, a.contact_atom_id,
						a.atom_sample_index, a.topology_id)
						< std::tie(b.source_face_id, b.source_node_index,
							b.boundary_occurrence_id, b.contact_atom_id,
							b.atom_sample_index, b.topology_id);
				});
			result.mesh_node_claims.erase(std::unique(result.mesh_node_claims.begin(),
				result.mesh_node_claims.end(), [](const auto& a, const auto& b)
				{
					return std::tie(a.source_face_id, a.source_node_index,
						a.boundary_occurrence_id, a.contact_atom_id,
						a.atom_sample_index, a.topology_id)
						== std::tie(b.source_face_id, b.source_node_index,
							b.boundary_occurrence_id, b.contact_atom_id,
							b.atom_sample_index, b.topology_id);
				}), result.mesh_node_claims.end());

			std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t>
				topology_by_face_node;
			std::map<std::tuple<std::uint32_t, std::uint32_t, std::uint64_t,
				std::uint32_t>,
				std::uint32_t> node_by_face_occurrence_topology;
			for (const auto& claim : result.mesh_node_claims)
			{
				const auto node_key = std::make_pair(claim.source_face_id,
					claim.source_node_index);
				auto [node, node_inserted] = topology_by_face_node.emplace(node_key,
					claim.topology_id);
				if (!node_inserted && node->second != claim.topology_id)
				{
					error = "one original face-mesh node belongs to conflicting contact topology groups: face "
						+ std::to_string(claim.source_face_id) + " node "
						+ std::to_string(claim.source_node_index) + " topology IDs "
						+ std::to_string(node->second) + " and "
						+ std::to_string(claim.topology_id);
					return false;
				}
				const auto occurrence_key = std::make_tuple(claim.source_face_id,
					claim.boundary_occurrence_id, claim.contact_atom_id,
					claim.atom_sample_index);
				auto [occurrence, occurrence_inserted] =
					node_by_face_occurrence_topology.emplace(occurrence_key,
						claim.source_node_index);
				if (!occurrence_inserted && occurrence->second != claim.source_node_index)
				{
					error = "one face/occurrence contact sample names multiple original mesh nodes: face "
						+ std::to_string(claim.source_face_id) + " occurrence "
						+ std::to_string(claim.boundary_occurrence_id) + " atom "
						+ std::to_string(claim.contact_atom_id) + " sample "
						+ std::to_string(claim.atom_sample_index) + " nodes "
						+ std::to_string(occurrence->second) + " and "
						+ std::to_string(claim.source_node_index);
					return false;
				}
			}
		}
		catch (const Standard_Failure& failure)
		{
			error = std::string("OpenCascade contact topology exception: ")
				+ failure.GetMessageString();
			return false;
		}
		catch (const std::exception& exception)
		{
			error = std::string("contact topology exception: ") + exception.what();
			return false;
		}
		catch (...)
		{
			error = "unknown contact topology exception";
			return false;
		}

		output = std::move(result);
		return true;
	}
}
