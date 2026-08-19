#include "core/geometry/exact_cell_decomposition.h"

#include <BOPAlgo_CellsBuilder.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Message_Report.hxx>
#include <NCollection_List.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_Triangle.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace paracfd::core
{
namespace
{
using EdgeKey = std::pair<std::uint32_t, std::uint32_t>;

ExactCellPoint point(const gp_Pnt &p) { return { p.X(), p.Y(), p.Z() }; }

gp_Pnt point(const ExactCellPoint &p) { return { p[0], p[1], p[2] }; }

ExactCellPoint subtract(const ExactCellPoint &a, const ExactCellPoint &b)
{
  return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

ExactCellPoint add(const ExactCellPoint &a, const ExactCellPoint &b)
{
  return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

ExactCellPoint multiply(const ExactCellPoint &a, double scale)
{
  return { a[0] * scale, a[1] * scale, a[2] * scale };
}

ExactCellPoint cross(const ExactCellPoint &a, const ExactCellPoint &b)
{
  return { a[1] * b[2] - a[2] * b[1],
	a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}

double dot(const ExactCellPoint &a, const ExactCellPoint &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const ExactCellPoint &value) { return std::sqrt(dot(value, value)); }

ExactCellPoint normalized(const ExactCellPoint &value)
{
  const double length = norm(value);
  if (!(length > 0.0)) throw std::runtime_error("zero triangle normal");
  return { value[0] / length, value[1] / length, value[2] / length };
}

bool finite(const ExactCellPoint &value)
{
  return std::isfinite(value[0]) && std::isfinite(value[1]) && std::isfinite(value[2]);
}

double polygon_area(const std::vector<ExactCellPoint> &polygon)
{
  if (polygon.size() < 3) return 0.0;
  ExactCellPoint twice_area{};
  for (std::size_t i = 0; i < polygon.size(); ++i)
    twice_area = add(twice_area, cross(polygon[i], polygon[(i + 1) % polygon.size()]));
  return 0.5 * norm(twice_area);
}

std::vector<ExactCellPoint> clip_polygon_axis(const std::vector<ExactCellPoint> &input,
    int axis, double bound, bool keep_upper)
{
  std::vector<ExactCellPoint> output;
  if (input.empty()) return output;
  auto inside = [&](const ExactCellPoint &point)
  { return keep_upper ? point[axis] >= bound : point[axis] <= bound; };
  ExactCellPoint previous = input.back();
  bool previous_inside = inside(previous);
  for (const ExactCellPoint &current : input)
  {
    const bool current_inside = inside(current);
    if (current_inside != previous_inside)
    {
      const double denominator = current[axis] - previous[axis];
      if (denominator != 0.0)
      {
        const double t = (bound - previous[axis]) / denominator;
        output.push_back(add(previous, multiply(subtract(current, previous), t)));
      }
    }
    if (current_inside) output.push_back(current);
    previous = current;
    previous_inside = current_inside;
  }
  return output;
}

double clipped_triangle_area(const ExactCellInput &input, const ExactCellTriangle &triangle)
{
  std::vector<ExactCellPoint> polygon{input.vertices[triangle.vertices[0]],
      input.vertices[triangle.vertices[1]], input.vertices[triangle.vertices[2]]};
  for (int axis = 0; axis < 3 && !polygon.empty(); ++axis)
  {
    polygon = clip_polygon_axis(polygon, axis, input.cell_min[axis], true);
    polygon = clip_polygon_axis(polygon, axis, input.cell_max[axis], false);
  }
  return polygon_area(polygon);
}

bool clipped_segment_has_positive_length(const ExactCellPoint &a,
    const ExactCellPoint &b, const ExactCellInput &input, double length_tolerance)
{
  const ExactCellPoint direction = subtract(b, a);
  double lower = 0.0, upper = 1.0;
  for (int axis = 0; axis < 3; ++axis)
  {
	if (std::abs(direction[axis]) <= std::numeric_limits<double>::min())
	{
	  if (a[axis] < input.cell_min[axis] || a[axis] > input.cell_max[axis])
		return false;
	  continue;
	}
	double first = (input.cell_min[axis] - a[axis]) / direction[axis];
	double last = (input.cell_max[axis] - a[axis]) / direction[axis];
	if (first > last) std::swap(first, last);
	lower = std::max(lower, first);
	upper = std::min(upper, last);
	if (upper < lower) return false;
  }
  return (upper - lower) * norm(direction) > length_tolerance;
}

ExactCellPoint input_normal(const ExactCellInput &input,
    const ExactCellTriangle &triangle)
{
  const ExactCellPoint &a = input.vertices[triangle.vertices[0]];
  const ExactCellPoint &b = input.vertices[triangle.vertices[1]];
  const ExactCellPoint &c = input.vertices[triangle.vertices[2]];
  return normalized(cross(subtract(b, a), subtract(c, a)));
}

ExactCellPoint oriented_face_normal(const TopoDS_Face &face)
{
  BRepAdaptor_Surface surface(face, true);
  const double u = 0.5 * (surface.FirstUParameter() + surface.LastUParameter());
  const double v = 0.5 * (surface.FirstVParameter() + surface.LastVParameter());
  gp_Pnt p;
  gp_Vec du, dv;
  surface.D1(u, v, p, du, dv);
  gp_Vec n = du.Crossed(dv);
  if (face.Orientation() == TopAbs_REVERSED) n.Reverse();
  const double magnitude = n.Magnitude();
  if (!(magnitude > 0.0)) throw std::runtime_error("OCCT face has a zero normal");
  return { n.X() / magnitude, n.Y() / magnitude, n.Z() / magnitude };
}

std::vector<std::vector<ExactCellPoint>> boundary_loops(const TopoDS_Face &face)
{
  std::vector<std::vector<ExactCellPoint>> loops;
  for (TopExp_Explorer wires(face, TopAbs_WIRE); wires.More(); wires.Next())
  {
	const TopoDS_Wire wire = TopoDS::Wire(wires.Current());
	std::vector<ExactCellPoint> loop;
	for (BRepTools_WireExplorer edges(wire, face); edges.More(); edges.Next())
	{
	  const TopoDS_Vertex vertex = edges.CurrentVertex();
	  if (!vertex.IsNull()) loop.push_back(point(BRep_Tool::Pnt(vertex)));
	}
	if (loop.size() >= 3) loops.push_back(std::move(loop));
  }
  return loops;
}

bool point_less(const ExactCellPoint &a, const ExactCellPoint &b)
{
  if (a[0] != b[0]) return a[0] < b[0];
  if (a[1] != b[1]) return a[1] < b[1];
  return a[2] < b[2];
}

void canonicalize_loop(std::vector<ExactCellPoint> &loop)
{
  if (loop.empty()) return;
  auto first = std::min_element(loop.begin(), loop.end(), point_less);
  std::rotate(loop.begin(), first, loop.end());
}

std::vector<std::array<ExactCellPoint, 3>> face_convex_pieces(const TopoDS_Face &face)
{
  TopLoc_Location location;
  const Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);
  if (triangulation.IsNull() || triangulation->NbTriangles() <= 0)
	throw std::runtime_error("planar face has no conservative triangulation");
  const gp_Trsf transform = location.Transformation();
  const bool reversed = face.Orientation() == TopAbs_REVERSED;
  std::vector<std::array<ExactCellPoint, 3>> pieces;
  pieces.reserve(triangulation->NbTriangles());
  for (Standard_Integer index = 1; index <= triangulation->NbTriangles(); ++index)
  {
	Standard_Integer ia, ib, ic;
	triangulation->Triangle(index).Get(ia, ib, ic);
	if (reversed) std::swap(ib, ic);
	const ExactCellPoint a = point(triangulation->Node(ia).Transformed(transform));
	const ExactCellPoint b = point(triangulation->Node(ib).Transformed(transform));
	const ExactCellPoint c = point(triangulation->Node(ic).Transformed(transform));
	const double area = 0.5 * norm(cross(subtract(b, a), subtract(c, a)));
	if (!(area > 0.0) || !std::isfinite(area))
	  throw std::runtime_error("planar face triangulation contains a degenerate triangle");
	pieces.push_back({ a, b, c });
  }
  return pieces;
}

ExactCellPlanarRegion planar_region(const TopoDS_Face &face)
{
  ExactCellPlanarRegion region;
  const TopoDS_Wire outer = BRepTools::OuterWire(face);
  if (outer.IsNull()) throw std::runtime_error("planar face has no outer wire");
  for (TopExp_Explorer wires(face, TopAbs_WIRE); wires.More(); wires.Next())
  {
	const TopoDS_Wire wire = TopoDS::Wire(wires.Current());
	std::vector<ExactCellPoint> loop;
	for (BRepTools_WireExplorer edges(wire, face); edges.More(); edges.Next())
	{
	  const TopoDS_Vertex vertex = edges.CurrentVertex();
	  if (!vertex.IsNull()) loop.push_back(point(BRep_Tool::Pnt(vertex)));
	}
	canonicalize_loop(loop);
	if (wire.IsSame(outer))
	{
	  if (loop.size() < 3) throw std::runtime_error("planar face outer wire has fewer than three vertices");
	  region.outer_loop = std::move(loop);
	}
	else if (loop.size() >= 3) region.hole_loops.push_back(std::move(loop));
	else if (!loop.empty()) region.internal_slit_wires.push_back(std::move(loop));
	else ++region.empty_internal_wire_count;
  }
  if (region.outer_loop.size() < 3) throw std::runtime_error("planar face outer wire was not enumerated");
  std::sort(region.hole_loops.begin(), region.hole_loops.end(),
      [](const auto &a, const auto &b) { return point_less(a.front(), b.front()); });
  std::sort(region.internal_slit_wires.begin(), region.internal_slit_wires.end(),
      [](const auto &a, const auto &b) { return point_less(a.front(), b.front()); });

  region.convex_pieces = face_convex_pieces(face);
  for (const auto &piece : region.convex_pieces)
  {
	const ExactCellPoint &a = piece[0], &b = piece[1], &c = piece[2];
	const double area = 0.5 * norm(cross(subtract(b, a), subtract(c, a)));
	region.convex_piece_area_sum += area;
	region.convex_piece_first_moment = add(region.convex_piece_first_moment,
	    multiply(add(add(a, b), c), area / 3.0));
  }
  return region;
}

double oriented_solid_angle(const ExactCellOrientedTriangle &triangle,
    const ExactCellPoint &query)
{
  const ExactCellPoint a = subtract(triangle.vertices[0], query);
  const ExactCellPoint b = subtract(triangle.vertices[1], query);
  const ExactCellPoint c = subtract(triangle.vertices[2], query);
  const double la = norm(a), lb = norm(b), lc = norm(c);
  const double numerator = dot(a, cross(b, c));
  const double denominator = la * lb * lc + dot(a, b) * lc +
      dot(b, c) * la + dot(c, a) * lb;
  return 2.0 * std::atan2(numerator, denominator);
}

double generalized_winding(const ExactCellFragment &fragment,
    const ExactCellPoint &query)
{
  double solid_angle = 0.0;
  for (const ExactCellOrientedTriangle &triangle : fragment.boundary_triangles)
	solid_angle += oriented_solid_angle(triangle, query);
  return solid_angle / (4.0 * std::acos(-1.0));
}

bool classified_inside(const TopoDS_Solid &solid, const ExactCellPoint &candidate)
{
  BRepClass3d_SolidClassifier classifier(solid, point(candidate), Precision::Confusion());
  return classifier.State() == TopAbs_IN;
}

std::optional<ExactCellPoint> certified_interior_witness(const TopoDS_Solid &solid,
    const ExactCellFragment &fragment, double cell_extent)
{
  if (classified_inside(solid, fragment.centroid)) return fragment.centroid;
  const double minimum_offset = 32.0 * Precision::Confusion();
  for (const ExactCellOrientedTriangle &triangle : fragment.boundary_triangles)
  {
	const ExactCellPoint &a = triangle.vertices[0];
	const ExactCellPoint &b = triangle.vertices[1];
	const ExactCellPoint &c = triangle.vertices[2];
	const ExactCellPoint area_vector = cross(subtract(b, a), subtract(c, a));
	const double twice_area = norm(area_vector);
	if (!(twice_area > 0.0)) continue;
	const ExactCellPoint outward = multiply(area_vector, 1.0 / twice_area);
	const ExactCellPoint centre = multiply(add(add(a, b), c), 1.0 / 3.0);
	const double local_scale = std::min(cell_extent, std::sqrt(0.5 * twice_area));
	for (double fraction : { 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6 })
	{
	  const double offset = std::max(minimum_offset, fraction * local_scale);
	  if (!(offset < 0.25 * local_scale)) continue;
	  const ExactCellPoint candidate = subtract(centre, multiply(outward, offset));
	  if (classified_inside(solid, candidate)) return candidate;
	}
  }
  return std::nullopt;
}

struct SurfaceMeasure
{
  double area = 0.0;
  ExactCellPoint centroid{};
};

SurfaceMeasure surface_measure(const TopoDS_Face &face)
{
  GProp_GProps properties;
  BRepGProp::SurfaceProperties(face, properties);
  return { std::abs(properties.Mass()), point(properties.CentreOfMass()) };
}

void validate_region(const ExactCellPlanarRegion &region, const SurfaceMeasure &measure,
    double area_tolerance, double length_scale, const std::string &label,
    ExactCellDecomposition &result)
{
  const double area_error = std::abs(region.convex_piece_area_sum - measure.area);
  if (area_error > area_tolerance)
	result.errors.push_back(label + " convex-piece area mismatch: residual=" + std::to_string(area_error));
  const ExactCellPoint expected = multiply(measure.centroid, measure.area);
  const double moment_error = norm(subtract(region.convex_piece_first_moment, expected));
  if (moment_error > area_tolerance * std::max(1.0, length_scale))
	result.errors.push_back(label + " convex-piece first-moment mismatch: residual=" + std::to_string(moment_error));
}

struct VolumeMeasure
{
  double volume = 0.0;
  ExactCellPoint centroid{};
};

VolumeMeasure volume_measure(const TopoDS_Solid &solid)
{
  GProp_GProps properties;
  BRepGProp::VolumeProperties(solid, properties);
  return { std::abs(properties.Mass()), point(properties.CentreOfMass()) };
}

std::vector<TopoDS_Face> history_images(BOPAlgo_CellsBuilder &builder,
    const TopoDS_Face &source, const TopTools_IndexedMapOfShape &all_faces)
{
  std::vector<TopoDS_Face> result;
  TopTools_MapOfShape seen;
  const NCollection_List<TopoDS_Shape> &modified = builder.Modified(source);
  for (NCollection_List<TopoDS_Shape>::Iterator it(modified); it.More(); it.Next())
  {
	const TopoDS_Shape &image = it.Value();
	if (image.ShapeType() != TopAbs_FACE || all_faces.FindIndex(image) == 0 || !seen.Add(image)) continue;
	result.push_back(TopoDS::Face(image));
  }
  // OCCT history may contain both modified pieces and an unchanged source image.
  // Never infer that the latter vanished merely because Modified() was non-empty.
  if (all_faces.FindIndex(source) != 0 && seen.Add(source)) result.push_back(source);
  return result;
}

std::vector<TopoDS_Edge> history_edge_images(BOPAlgo_CellsBuilder &builder,
    const TopoDS_Edge &source, const TopTools_IndexedMapOfShape &all_edges)
{
  std::vector<TopoDS_Edge> result;
  TopTools_MapOfShape seen;
  const NCollection_List<TopoDS_Shape> &modified = builder.Modified(source);
  for (NCollection_List<TopoDS_Shape>::Iterator it(modified); it.More(); it.Next())
  {
	const TopoDS_Shape &image = it.Value();
	if (image.ShapeType() != TopAbs_EDGE || all_edges.FindIndex(image) == 0
	    || !seen.Add(image)) continue;
	result.push_back(TopoDS::Edge(image));
  }
  if (all_edges.FindIndex(source) != 0 && seen.Add(source)) result.push_back(source);
  return result;
}

struct BuiltTriangleGeometry
{
  std::vector<TopoDS_Face> faces;
  std::vector<std::array<TopoDS_Edge, 3>> edges;
};

BuiltTriangleGeometry build_triangle_geometry(const ExactCellInput &input)
{
  std::vector<TopoDS_Vertex> vertices;
  vertices.reserve(input.vertices.size());
  for (const ExactCellPoint &vertex : input.vertices)
	vertices.push_back(BRepBuilderAPI_MakeVertex(point(vertex)).Vertex());

  std::map<EdgeKey, TopoDS_Edge> edges;
  BuiltTriangleGeometry result;
  result.faces.reserve(input.triangles.size());
  result.edges.reserve(input.triangles.size());
  for (const ExactCellTriangle &triangle : input.triangles)
  {
	BRepBuilderAPI_MakeWire make_wire;
	std::array<TopoDS_Edge, 3> triangle_edges;
	for (int side = 0; side < 3; ++side)
	{
	  const std::uint32_t a = triangle.vertices[side];
	  const std::uint32_t b = triangle.vertices[(side + 1) % 3];
	  const EdgeKey key = std::minmax(a, b);
	  auto found = edges.find(key);
	  if (found == edges.end())
	  {
		const TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(
		    vertices[key.first], vertices[key.second])
		                             .Edge();
		found = edges.emplace(key, edge).first;
	  }
	  TopoDS_Edge directed = found->second;
	  if (a != key.first) directed.Reverse();
	  triangle_edges[side] = directed;
	  make_wire.Add(directed);
	}
	if (!make_wire.IsDone()) throw std::runtime_error("failed to build triangle wire");
	BRepBuilderAPI_MakeFace make_face(make_wire.Wire(), true);
	if (!make_face.IsDone()) throw std::runtime_error("failed to build triangle face");
	TopoDS_Face face = make_face.Face();
	const ExactCellPoint wanted = input_normal(input, triangle);
	if (dot(oriented_face_normal(face), wanted) < 0.0) face.Reverse();
	result.faces.push_back(face);
	result.edges.push_back(std::move(triangle_edges));
  }
  return result;
}

bool validate_input(const ExactCellInput &input, ExactCellDecomposition &result)
{
  for (int axis = 0; axis < 3; ++axis)
	if (!std::isfinite(input.cell_min[axis]) || !std::isfinite(input.cell_max[axis]) || !(input.cell_max[axis] > input.cell_min[axis]))
	  result.errors.push_back("cell bounds must be finite and strictly increasing");
  for (std::size_t i = 0; i < input.vertices.size(); ++i)
	if (!finite(input.vertices[i]))
	  result.errors.push_back("vertex " + std::to_string(i) + " is not finite");

  std::set<std::uint64_t> triangle_ids;
  const double extent = std::max({ input.cell_max[0] - input.cell_min[0],
      input.cell_max[1] - input.cell_min[1], input.cell_max[2] - input.cell_min[2] });
  const double minimum_cross = std::numeric_limits<double>::epsilon() * std::max(1.0, extent * extent) * 64.0;
  for (std::size_t i = 0; i < input.triangles.size(); ++i)
  {
	const ExactCellTriangle &triangle = input.triangles[i];
	if (!triangle_ids.insert(triangle.triangle_id).second)
	  result.errors.push_back("duplicate source triangle ID " + std::to_string(triangle.triangle_id));
	bool valid_indices = true;
	for (std::uint32_t index : triangle.vertices)
	  valid_indices = valid_indices && index < input.vertices.size();
	if (!valid_indices)
	{
	  result.errors.push_back("triangle " + std::to_string(i) + " has an out-of-range vertex index");
	  continue;
	}
	if (triangle.vertices[0] == triangle.vertices[1] || triangle.vertices[1] == triangle.vertices[2] || triangle.vertices[2] == triangle.vertices[0] || norm(cross(subtract(input.vertices[triangle.vertices[1]], input.vertices[triangle.vertices[0]]), subtract(input.vertices[triangle.vertices[2]], input.vertices[triangle.vertices[0]]))) <= minimum_cross)
	  result.errors.push_back("triangle " + std::to_string(i) + " is degenerate");
	for (int edge = 0; edge < 3; ++edge)
	{
	  const bool has_contact = triangle.cad_contact_ids[edge]
	      != ExactCellTriangle::no_cad_contact_id;
	  const std::uint16_t fan = triangle.certified_fan_degrees[edge];
	  if (has_contact != (fan >= 2))
		result.errors.push_back("triangle " + std::to_string(i) + " edge "
		    + std::to_string(edge)
		    + " has an incomplete CAD contact/fan certificate");
	}
  }
  return result.errors.empty();
}

std::vector<std::pair<int, TopoDS_Face>> incidences(const TopoDS_Face &image,
    const std::vector<TopoDS_Solid> &solids)
{
  std::vector<std::pair<int, TopoDS_Face>> result;
  for (int fragment = 0; fragment < static_cast<int>(solids.size()); ++fragment)
	for (TopExp_Explorer faces(solids[fragment], TopAbs_FACE); faces.More(); faces.Next())
	{
	  const TopoDS_Face occurrence = TopoDS::Face(faces.Current());
	  if (occurrence.IsSame(image)) result.emplace_back(fragment, occurrence);
	}
  return result;
}

bool lies_on_box_face(const TopoDS_Face &face, const ExactCellInput &input,
    int *axis_out = nullptr, bool *upper_out = nullptr)
{
  Bnd_Box bounds;
  BRepBndLib::AddOptimal(face, bounds, false, false);
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const std::array<double, 3> lo{ xmin, ymin, zmin };
  const std::array<double, 3> hi{ xmax, ymax, zmax };
  const double scale = std::max({ 1.0, input.cell_max[0] - input.cell_min[0],
      input.cell_max[1] - input.cell_min[1], input.cell_max[2] - input.cell_min[2] });
  const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() * scale;
  for (int axis = 0; axis < 3; ++axis)
  {
	if (std::abs(lo[axis] - input.cell_min[axis]) <= tolerance && std::abs(hi[axis] - input.cell_min[axis]) <= tolerance)
	{
	  if (axis_out) *axis_out = axis;
	  if (upper_out) *upper_out = false;
	  return true;
	}
	if (std::abs(lo[axis] - input.cell_max[axis]) <= tolerance && std::abs(hi[axis] - input.cell_max[axis]) <= tolerance)
	{
	  if (axis_out) *axis_out = axis;
	  if (upper_out) *upper_out = true;
	  return true;
	}
  }
  return false;
}

std::string triangle_label(const ExactCellTriangle &triangle)
{
  return "source triangle " + std::to_string(triangle.triangle_id);
}

void populate_empty_cell(const ExactCellInput &input, ExactCellDecomposition &result)
{
  BRepPrimAPI_MakeBox make_box(point(input.cell_min), point(input.cell_max));
  const TopoDS_Solid box = make_box.Solid();
  const double cell_extent = std::max({ input.cell_max[0] - input.cell_min[0],
      input.cell_max[1] - input.cell_min[1], input.cell_max[2] - input.cell_min[2] });
  BRepMesh_IncrementalMesh mesher(box,
      std::max(10.0 * Precision::Confusion(), 1.0e-6 * cell_extent), false, 0.1, false);
  if (!mesher.IsDone())
    throw std::runtime_error("could not triangulate an empty exact Cartesian cell");

  ExactCellFragment fragment;
  fragment.id = 0;
  fragment.volume = result.expected_cell_volume;
  for (int axis = 0; axis < 3; ++axis)
    fragment.centroid[axis] = fragment.interior_witness[axis]
        = 0.5 * (input.cell_min[axis] + input.cell_max[axis]);

  for (TopExp_Explorer faces(box, TopAbs_FACE); faces.More(); faces.Next())
  {
    const TopoDS_Face occurrence = TopoDS::Face(faces.Current());
    const auto boundary = face_convex_pieces(occurrence);
    for (const auto &piece : boundary)
      fragment.boundary_triangles.push_back({ piece });

    int axis = -1;
    bool upper = false;
    if (!lies_on_box_face(occurrence, input, &axis, &upper))
      throw std::runtime_error("could not identify an empty cell box face");
    const SurfaceMeasure measure = surface_measure(occurrence);
    ExactCellBoxAperture aperture;
    aperture.axis = static_cast<std::int8_t>(axis);
    aperture.upper = upper;
    aperture.area = measure.area;
    aperture.centroid = measure.centroid;
    aperture.region = planar_region(occurrence);
    aperture.fragment = 0;
    result.box_apertures.push_back(std::move(aperture));
  }
  if (fragment.boundary_triangles.empty())
    throw std::runtime_error("empty exact Cartesian cell has no retained boundary");
  const double winding = generalized_winding(fragment, fragment.interior_witness);
  if (!std::isfinite(winding) || std::abs(std::abs(winding) - 1.0) > 1.0e-7)
    throw std::runtime_error("empty exact Cartesian cell retained boundary is invalid");

  result.fragments.push_back(std::move(fragment));
  result.fragment_volume_sum = result.expected_cell_volume;
  result.requested_fuzzy_tolerance = 0.0;
  result.effective_fuzzy_tolerance = Precision::Confusion();
}
} // namespace

ExactCellDecomposition decompose_exact_cell(const ExactCellInput &input)
{
  ExactCellDecomposition result;
  result.expected_cell_volume = (input.cell_max[0] - input.cell_min[0]) * (input.cell_max[1] - input.cell_min[1]) * (input.cell_max[2] - input.cell_min[2]);
  if (!validate_input(input, result)) return result;

  try
  {
	if (input.triangles.empty())
	{
	  populate_empty_cell(input, result);
	  return result;
	}
	const BuiltTriangleGeometry triangle_geometry = build_triangle_geometry(input);
	const std::vector<TopoDS_Face> &triangle_faces = triangle_geometry.faces;
	BRepPrimAPI_MakeBox make_box(point(input.cell_min), point(input.cell_max));
	const TopoDS_Solid box = make_box.Solid();
	std::vector<TopoDS_Face> box_faces;
	for (TopExp_Explorer faces(box, TopAbs_FACE); faces.More(); faces.Next())
	  box_faces.push_back(TopoDS::Face(faces.Current()));

	NCollection_List<TopoDS_Shape> arguments;
	arguments.Append(box);
	for (const TopoDS_Face &face : triangle_faces) arguments.Append(face);

	BOPAlgo_CellsBuilder builder;
	builder.SetArguments(arguments);
	builder.SetNonDestructive(true);
	builder.SetRunParallel(input.run_parallel);
	builder.SetUseOBB(true);
	builder.SetFuzzyValue(0.0); // Zero requested; OCCT exposes its unavoidable floor below.
	result.requested_fuzzy_tolerance = 0.0;
	result.effective_fuzzy_tolerance = builder.FuzzyValue();
	const auto begin = std::chrono::steady_clock::now();
	builder.Perform();
	const auto end = std::chrono::steady_clock::now();
	result.general_fuse_milliseconds =
	    std::chrono::duration<double, std::milli>(end - begin).count();
	if (builder.HasErrors())
	{
	  std::ostringstream report;
	  builder.DumpErrors(report);
	  result.errors.push_back("OCCT General Fuse failed:\n" + report.str());
	  return result;
	}
	if (builder.HasWarnings())
	{
	  std::ostringstream report;
	  builder.DumpWarnings(report);
	  const std::string message="OCCT General Fuse warning:\n"+report.str();
	  result.warnings.push_back(message);
	  // A warning can denote discarded construction debris, but it can also mean a
	  // membrane split was lost while area and volume still happen to balance.  Until
	  // fragment connectivity is independently cross-checked, warning-only results
	  // are not safe enough to become pressure topology.
	  result.errors.push_back("General Fuse emitted warnings; cell topology is not accepted");
	  return result;
	}

	builder.AddAllToResult();
	const TopoDS_Shape &all = builder.Shape();
	if (all.IsNull())
	{
	  result.errors.push_back("OCCT General Fuse returned a null all-parts result");
	  return result;
	}
	if (!BRepCheck_Analyzer(all, true).IsValid())
	  result.errors.push_back("OCCT General Fuse returned invalid BRep topology");
	const double cell_extent = std::max({ input.cell_max[0] - input.cell_min[0],
	    input.cell_max[1] - input.cell_min[1], input.cell_max[2] - input.cell_min[2] });
	BRepMesh_IncrementalMesh result_mesher(all,
	    std::max(10.0 * Precision::Confusion(), 1.0e-6 * cell_extent), false, 0.1, false);
	if (!result_mesher.IsDone())
	{
	  result.errors.push_back("could not triangulate the exact cell decomposition for OCC-free retention");
	  return result;
	}

	TopTools_IndexedMapOfShape all_faces;
	TopExp::MapShapes(all, TopAbs_FACE, all_faces);
	TopTools_IndexedMapOfShape all_edges;
	TopExp::MapShapes(all, TopAbs_EDGE, all_edges);
	TopTools_IndexedMapOfShape solid_map;
	TopExp::MapShapes(all, TopAbs_SOLID, solid_map);
	std::vector<TopoDS_Solid> solids;
	solids.reserve(solid_map.Extent());
	ExactCellPoint first_moment{};
	for (int index = 1; index <= solid_map.Extent(); ++index)
	{
	  const TopoDS_Solid solid = TopoDS::Solid(solid_map(index));
	  if (!BRepCheck_Analyzer(solid, true).IsValid())
		result.errors.push_back("fragment " + std::to_string(index - 1) + " is not a valid OCCT solid");
	  const VolumeMeasure measure = volume_measure(solid);
	  if (!(measure.volume > 0.0) || !std::isfinite(measure.volume))
		result.errors.push_back("fragment " + std::to_string(index - 1) + " has non-positive or non-finite volume");
	  ExactCellFragment fragment;
	  fragment.id = index - 1;
	  fragment.volume = measure.volume;
	  fragment.centroid = measure.centroid;
	  for (TopExp_Explorer faces(solid, TopAbs_FACE); faces.More(); faces.Next())
	  {
		const TopoDS_Face occurrence = TopoDS::Face(faces.Current());
		if (occurrence.Orientation() == TopAbs_INTERNAL) continue;
		if (occurrence.Orientation() == TopAbs_EXTERNAL)
		{
		  result.errors.push_back("fragment " + std::to_string(index - 1)
		      + " has an EXTERNAL boundary-face occurrence");
		  continue;
		}
		const auto boundary = face_convex_pieces(occurrence);
		for (const auto &piece : boundary)
		  fragment.boundary_triangles.push_back({ piece });
	  }
	  if (fragment.boundary_triangles.empty())
		result.errors.push_back("fragment " + std::to_string(index - 1)
		    + " has no retained closed boundary triangles");
	  else
	  {
		const std::optional<ExactCellPoint> witness =
		    certified_interior_witness(solid, fragment, cell_extent);
		if (!witness)
		  result.errors.push_back("fragment " + std::to_string(index - 1)
		      + " has no BRepClass3d-certified interior witness");
		else
		{
		  fragment.interior_witness = *witness;
		  const double winding = generalized_winding(fragment, *witness);
		  if (!std::isfinite(winding) || std::abs(std::abs(winding) - 1.0) > 1.0e-7)
			result.errors.push_back("fragment " + std::to_string(index - 1)
			    + " retained boundary does not classify its certified interior witness; winding="
			    + std::to_string(winding));
		}
	  }
	  result.fragments.push_back(std::move(fragment));
	  result.fragment_volume_sum += measure.volume;
	  for (int axis = 0; axis < 3; ++axis)
		first_moment[axis] += measure.volume * measure.centroid[axis];
	  solids.push_back(solid);
	}
	if (solids.empty()) result.errors.push_back("General Fuse produced no fluid fragments");
	for (std::size_t expected = 0; expected < result.fragments.size(); ++expected)
	{
	  const ExactCellPoint &witness = result.fragments[expected].interior_witness;
	  int containing = -1;
	  for (std::size_t candidate = 0; candidate < result.fragments.size(); ++candidate)
	  {
		const double magnitude = std::abs(generalized_winding(result.fragments[candidate], witness));
		if (magnitude <= 1.0e-7) continue;
		if (std::abs(magnitude - 1.0) > 1.0e-7 || containing >= 0)
		{
		  containing = -2;
		  break;
		}
		containing = static_cast<int>(candidate);
	  }
	  if (containing != static_cast<int>(expected))
		result.errors.push_back("fragment " + std::to_string(expected)
		    + " retained locator does not uniquely own its certified interior witness");
	}

	result.volume_conservation_residual =
	    result.fragment_volume_sum - result.expected_cell_volume;
	result.relative_volume_conservation_residual = result.expected_cell_volume > 0.0
	    ? result.volume_conservation_residual / result.expected_cell_volume
	    : 0.0;
	ExactCellPoint expected_moment{};
	for (int axis = 0; axis < 3; ++axis)
	  expected_moment[axis] = result.expected_cell_volume * 0.5 * (input.cell_min[axis] + input.cell_max[axis]);
	result.first_moment_conservation_residual = norm(subtract(first_moment, expected_moment));
	result.relative_first_moment_conservation_residual=result.expected_cell_volume>0.0&&cell_extent>0.0
	    ?result.first_moment_conservation_residual/(result.expected_cell_volume*cell_extent):0.0;
	const double relative_error = std::abs(result.relative_volume_conservation_residual);
	if (relative_error > 1.0e-8)
	  result.errors.push_back("fragment volume conservation failed: relative residual=" + std::to_string(result.relative_volume_conservation_residual));
	else if (relative_error > 1.0e-11)
	  result.warnings.push_back("fragment volume conservation is degraded: relative residual=" + std::to_string(result.relative_volume_conservation_residual));
	if(result.relative_first_moment_conservation_residual>1.0e-8)
	  result.errors.push_back("fragment first-moment conservation failed: relative residual="+std::to_string(result.relative_first_moment_conservation_residual));
	else if(result.relative_first_moment_conservation_residual>1.0e-11)
	  result.warnings.push_back("fragment first-moment conservation is degraded: relative residual="+std::to_string(result.relative_first_moment_conservation_residual));

	TopTools_MapOfShape fabric_images;
	TopTools_MapOfShape claimed_fabric_images;
	std::vector<TopoDS_Face> retained_fabric_images;
	std::array<double,6> blocked_box_area{};
	std::array<ExactCellPoint,6> blocked_box_moment{};
	const double area_tolerance=std::max(1024.0*std::numeric_limits<double>::epsilon()*cell_extent*cell_extent,
	    8.0*result.effective_fuzzy_tolerance*cell_extent);
	for (std::size_t triangle_index = 0; triangle_index < triangle_faces.size(); ++triangle_index)
	{
	  const ExactCellTriangle &source = input.triangles[triangle_index];
	  const ExactCellPoint normal = input_normal(input, source);
	  const double expected_area=clipped_triangle_area(input,source);
	  result.input_clipped_surface_area_sum+=expected_area;
	  double retained_area=0.0;
	  const std::vector<TopoDS_Face> images =
	      history_images(builder, triangle_faces[triangle_index], all_faces);
	  for (const TopoDS_Face &image : images)
	  {
		const auto incident = incidences(image, solids);
		if (incident.empty()) continue; // the portion of a tool lying outside the box
		if (!claimed_fabric_images.Add(image))
		{
		  result.errors.push_back(triangle_label(source)
		      + " shares a positive-area General Fuse history image with another input triangle");
		  continue;
		}
		if (fabric_images.Add(image)) retained_fabric_images.push_back(image);
		const SurfaceMeasure measure = surface_measure(image);
		if (!(measure.area > 0.0)||!std::isfinite(measure.area))
		{
		  result.errors.push_back(triangle_label(source)+" produced a non-positive or non-finite history area");
		  continue;
		}
		retained_area+=measure.area;result.history_surface_area_sum+=measure.area;
		int boundary_axis=-1;bool boundary_upper=false;
		if(lies_on_box_face(image,input,&boundary_axis,&boundary_upper))
		{
		  const int side=2*boundary_axis+(boundary_upper?1:0);
		  blocked_box_area[side]+=measure.area;
		  blocked_box_moment[side]=add(blocked_box_moment[side],multiply(measure.centroid,measure.area));
		}
		std::set<int> plus, minus;
		for (const auto &[fragment, occurrence] : incident)
		{
		  if (occurrence.Orientation() == TopAbs_INTERNAL)
		  {
			plus.insert(fragment);
			minus.insert(fragment);
			continue;
		  }
		  if (occurrence.Orientation() == TopAbs_EXTERNAL)
		  {
			result.errors.push_back(triangle_label(source)
			    + " produced an EXTERNAL solid-face occurrence");
			continue;
		  }
		  // Boundary-face orientation is outward from the solid.  An outward
		  // normal aligned with the input plus direction therefore has the
		  // solid (fluid fragment) on the minus side, and vice versa.
		  const double alignment = dot(oriented_face_normal(occurrence), normal);
		  if (alignment > 0.5) minus.insert(fragment);
		  else if (alignment < -0.5) plus.insert(fragment);
		  else result.errors.push_back(triangle_label(source) + " produced a history image with an inconsistent normal");
		}
		if (plus.size() > 1 || minus.size() > 1)
		  result.errors.push_back(triangle_label(source) + " has more than one adjacent fragment on one side");
		const bool on_boundary = lies_on_box_face(image, input);
		if ((!on_boundary && (plus.empty() || minus.empty())) || (on_boundary && plus.empty() && minus.empty()))
		  result.errors.push_back(triangle_label(source) + " has incomplete General Fuse fragment adjacency");
		ExactCellPlanarRegion region = planar_region(image);
		validate_region(region, measure, area_tolerance, cell_extent,
		    triangle_label(source) + " history image", result);
		ExactCellSurfacePatch patch;
		patch.source_triangle_id = source.triangle_id;
		patch.source_face_id = source.source_face_id;
		patch.area = measure.area;
		patch.centroid = measure.centroid;
		patch.normal = normal;
		patch.region = std::move(region);
		patch.plus_fragment = plus.empty() ? -1 : *plus.begin();
		patch.minus_fragment = minus.empty() ? -1 : *minus.begin();
		patch.boundary_axis = static_cast<std::int8_t>(boundary_axis);
		patch.boundary_upper = boundary_upper;
		result.surface_patches.push_back(std::move(patch));
	  }
	  if(std::abs(retained_area-expected_area)>area_tolerance)
	    result.errors.push_back(triangle_label(source)+" clipped-area history mismatch: expected="+
	        std::to_string(expected_area)+", retained="+std::to_string(retained_area));
	}

	// A BRep-certified contact is stronger evidence than polygonal coincidence.
	// General Fuse may accept two independently tessellated sheets that miss by a
	// small amount, but doing so would turn an attached seam into an aerodynamic
	// opening.  Audit the result topology itself: an image of the certified source
	// edge must be shared by at least the declared number of retained fabric-face
	// sectors.  Topological edge identity is used here; the CAD tolerance is never
	// repurposed as a sewing/fuzzy radius.
	// Contact certification predates tessellation and remains meaningful below
	// OCCT's Boolean confusion floor.  Therefore only arithmetic roundoff may
	// classify a clipped segment as zero length; the kernel tolerance must not
	// erase the very certificate this audit is meant to protect.
	const double contact_length_tolerance = 512.0
	    * std::numeric_limits<double>::epsilon() * std::max(1.0, cell_extent);
	for (std::size_t triangle_index = 0; triangle_index < input.triangles.size();
	    ++triangle_index)
	{
	  const ExactCellTriangle &source = input.triangles[triangle_index];
	  for (int local_edge = 0; local_edge < 3; ++local_edge)
	  {
		const std::uint16_t expected_fan = source.certified_fan_degrees[local_edge];
		if (expected_fan < 2 || source.cad_contact_ids[local_edge]
		    == ExactCellTriangle::no_cad_contact_id) continue;
		const ExactCellPoint &a = input.vertices[source.vertices[local_edge]];
		const ExactCellPoint &b = input.vertices[source.vertices[(local_edge + 1) % 3]];
		if (!clipped_segment_has_positive_length(a, b, input,
		    contact_length_tolerance)) continue;

		++result.certified_contact_segments_checked;
		std::size_t minimum_reproduced_fan = std::numeric_limits<std::size_t>::max();
		bool found_retained_segment = false;
		const auto edge_images = history_edge_images(builder,
		    triangle_geometry.edges[triangle_index][local_edge], all_edges);
		for (const TopoDS_Edge &edge_image : edge_images)
		{
		  std::size_t incidence = 0;
		  for (const TopoDS_Face &fabric_image : retained_fabric_images)
		  {
			bool contains = false;
			for (TopExp_Explorer edges(fabric_image, TopAbs_EDGE); edges.More(); edges.Next())
			  if (edges.Current().IsSame(edge_image))
			  {
				contains = true;
				break;
			  }
			if (contains) ++incidence;
		  }
		  // Zero incidence identifies a split edge piece outside the Cartesian
		  // cell. Every positive-length piece retained on an in-cell fabric image
		  // must carry the full certified fan, not merely one lucky sub-piece.
		  if (incidence > 0)
		  {
			found_retained_segment = true;
			minimum_reproduced_fan = std::min(minimum_reproduced_fan, incidence);
		  }
		}
		if (!found_retained_segment || minimum_reproduced_fan < expected_fan)
		{
		  ++result.certified_contact_segments_unreproduced;
		  result.errors.push_back(triangle_label(source) + " edge "
		      + std::to_string(local_edge) + " CAD contact "
		      + std::to_string(source.cad_contact_ids[local_edge])
		      + " was not explicitly reproduced by General Fuse: expected fan="
		      + std::to_string(expected_fan) + ", minimum reproduced fan="
		      + (found_retained_segment ? std::to_string(minimum_reproduced_fan)
		                                  : std::string("none")));
		}
	  }
	}

	std::array<double,6> open_box_area{};
	std::array<ExactCellPoint,6> open_box_moment{};
	for (const TopoDS_Face &source_face : box_faces)
	{
	  int axis = -1;
	  bool upper = false;
	  if (!lies_on_box_face(source_face, input, &axis, &upper))
	  {
		result.errors.push_back("could not identify an original box face");
		continue;
	  }
	  for (const TopoDS_Face &image : history_images(builder, source_face, all_faces))
	  {
		// A positive-area image shared with a fabric tool is a blocked
		// Cartesian-face patch, not an open aperture.
		if (fabric_images.Contains(image)) continue;
		const auto incident = incidences(image, solids);
		std::set<int> adjacent;
		for (const auto &entry : incident) adjacent.insert(entry.first);
		if (adjacent.empty()) continue;
		if (adjacent.size() != 1)
		{
		  result.errors.push_back("box-face image has ambiguous fragment adjacency");
		  continue;
		}
		const SurfaceMeasure measure = surface_measure(image);
		if (!(measure.area > 0.0)||!std::isfinite(measure.area))
		{
		  result.errors.push_back("box-face aperture has non-positive or non-finite area");
		  continue;
		}
		const int side=2*axis+(upper?1:0);open_box_area[side]+=measure.area;
		open_box_moment[side]=add(open_box_moment[side],multiply(measure.centroid,measure.area));
		ExactCellPlanarRegion region = planar_region(image);
		validate_region(region, measure, area_tolerance, cell_extent,
		    "box-face aperture", result);
		ExactCellBoxAperture aperture;
		aperture.axis = static_cast<std::int8_t>(axis);
		aperture.upper = upper;
		aperture.area = measure.area;
		aperture.centroid = measure.centroid;
		aperture.region = std::move(region);
		aperture.fragment = *adjacent.begin();
		result.box_apertures.push_back(std::move(aperture));
	  }
	}
	for(int axis=0;axis<3;++axis)for(int upper=0;upper<2;++upper)
	{
	  const int side=2*axis+upper;double expected_area=1.0;
	  ExactCellPoint expected_centroid{};
	  for(int q=0;q<3;++q)
	  {
	    expected_centroid[q]=q==axis?(upper?input.cell_max[q]:input.cell_min[q])
	        :0.5*(input.cell_min[q]+input.cell_max[q]);
	    if(q!=axis)expected_area*=input.cell_max[q]-input.cell_min[q];
	  }
	  const double represented_area=open_box_area[side]+blocked_box_area[side];
	  const double area_residual=std::abs(represented_area-expected_area);
	  result.maximum_box_face_area_residual=std::max(result.maximum_box_face_area_residual,area_residual);
	  const ExactCellPoint represented_moment=add(open_box_moment[side],blocked_box_moment[side]);
	  const double moment_residual=norm(subtract(represented_moment,multiply(expected_centroid,expected_area)));
	  result.maximum_box_face_first_moment_residual=std::max(result.maximum_box_face_first_moment_residual,moment_residual);
	  if(area_residual>area_tolerance)
	    result.errors.push_back("box face "+std::to_string(side)+" area partition failed: residual="+std::to_string(area_residual));
	  const double moment_tolerance=area_tolerance*std::max(1.0,cell_extent);
	  if(moment_residual>moment_tolerance)
	    result.errors.push_back("box face "+std::to_string(side)+" first-moment partition failed: residual="+std::to_string(moment_residual));
	}

  }
  catch (const Standard_Failure &failure)
  {
	result.errors.push_back(std::string("OCCT exception: ") + (failure.GetMessageString() ? failure.GetMessageString() : "unknown failure"));
  }
  catch (const std::exception &failure)
  {
	result.errors.push_back(std::string("cell decomposition exception: ") + failure.what());
  }
  return result;
}
} // namespace paracfd::core
