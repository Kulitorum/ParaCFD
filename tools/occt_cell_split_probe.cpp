// Standalone acceptance probe for the OCC-isolated exact cell decomposer.
// This translation unit intentionally includes no OpenCascade header.

#include "core/geometry/exact_cell_decomposition.h"
#include "core/geometry/embedded_boundary.h"
#include "core/geometry/triangle_bvh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace
{
using paracfd::core::ExactCellDecomposition;
using paracfd::core::ExactCellInput;
using paracfd::core::ExactCellPoint;
using paracfd::core::ExactCellTriangle;
using paracfd::core::EmbeddedBoundary;
using paracfd::core::EmbeddedBoundaryBuildOptions;
using paracfd::core::FragmentRef;
using paracfd::core::TriMesh;
using paracfd::core::TriangleBvh;
using paracfd::core::UniformEbGrid;
using paracfd::core::Vec3d;

ExactCellPoint add(const ExactCellPoint &a, const ExactCellPoint &b)
{
  return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

ExactCellPoint multiply(const ExactCellPoint &a, double value)
{
  return { a[0] * value, a[1] * value, a[2] * value };
}

struct InputBuilder
{
  ExactCellInput input{ { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 } };
  std::uint64_t next_triangle_id = 1;

  void grid(const ExactCellPoint &origin, const ExactCellPoint &u,
      const ExactCellPoint &v, int nu, int nv, std::uint64_t source_face_id)
  {
	const std::uint32_t first = static_cast<std::uint32_t>(input.vertices.size());
	for (int j = 0; j <= nv; ++j)
	  for (int i = 0; i <= nu; ++i)
		input.vertices.push_back(add(origin, add(multiply(u, static_cast<double>(i) / nu), multiply(v, static_cast<double>(j) / nv))));
	auto vertex = [&](int i, int j)
	{
	  return first + static_cast<std::uint32_t>(j * (nu + 1) + i);
	};
	for (int j = 0; j < nv; ++j)
	  for (int i = 0; i < nu; ++i)
	  {
		const std::uint32_t p00 = vertex(i, j), p10 = vertex(i + 1, j);
		const std::uint32_t p01 = vertex(i, j + 1), p11 = vertex(i + 1, j + 1);
		input.triangles.push_back({ { p00, p10, p11 }, next_triangle_id++, source_face_id });
		input.triangles.push_back({ { p00, p11, p01 }, next_triangle_id++, source_face_id });
	  }
  }
};

struct Expected
{
  std::vector<double> volumes;
  double surface_area = 0.0;
  bool expect_same_fragment_patch = false;
  bool expect_separating_patch = false;
};

bool close(double a, double b, double relative = 2.0e-10)
{
  return std::abs(a - b) <= relative * std::max({ 1.0, std::abs(a), std::abs(b) });
}

bool check_case(const std::string &name, const ExactCellInput &input,
    const Expected &expected)
{
  const ExactCellDecomposition result = paracfd::core::decompose_exact_cell(input);
  bool pass = result.valid();
  std::vector<double> volumes;
  for (const auto &fragment : result.fragments) volumes.push_back(fragment.volume);
  std::sort(volumes.begin(), volumes.end());
  std::vector<double> wanted = expected.volumes;
  std::sort(wanted.begin(), wanted.end());
  pass = pass && volumes.size() == wanted.size();
  for (std::size_t i = 0; i < std::min(volumes.size(), wanted.size()); ++i)
	pass = pass && close(volumes[i], wanted[i]);

  const double patch_area = std::accumulate(result.surface_patches.begin(),
      result.surface_patches.end(), 0.0, [](double sum, const auto &patch)
      { return sum + patch.area; });
  pass = pass && close(patch_area, expected.surface_area);
  bool has_same = false, has_separating = false;
  for (const auto &patch : result.surface_patches)
  {
	pass = pass && patch.plus_fragment >= 0 && patch.minus_fragment >= 0;
	pass = pass && !patch.region.outer_loop.empty() && !patch.region.convex_pieces.empty();
	pass = pass && close(patch.region.convex_piece_area_sum, patch.area);
	has_same = has_same || patch.plus_fragment == patch.minus_fragment;
	has_separating = has_separating || patch.plus_fragment != patch.minus_fragment;
	if (patch.plus_fragment >= 0 && patch.minus_fragment >= 0 && patch.plus_fragment != patch.minus_fragment)
	{
	  const auto &plus = result.fragments[static_cast<std::size_t>(patch.plus_fragment)];
	  const auto &minus = result.fragments[static_cast<std::size_t>(patch.minus_fragment)];
	  double plus_side = 0.0, minus_side = 0.0;
	  for (int axis = 0; axis < 3; ++axis)
	  {
		plus_side += (plus.centroid[axis] - patch.centroid[axis]) * patch.normal[axis];
		minus_side += (minus.centroid[axis] - patch.centroid[axis]) * patch.normal[axis];
	  }
	  pass = pass && plus_side > 0.0 && minus_side < 0.0;
	}
  }
  pass = pass && (!expected.expect_same_fragment_patch || has_same);
  pass = pass && (!expected.expect_separating_patch || has_separating);

  std::array<double, 6> boundary_area{};
  for (const auto &aperture : result.box_apertures)
  {
	const int side = 2 * aperture.axis + (aperture.upper ? 1 : 0);
	pass = pass && side >= 0 && side < 6 && aperture.fragment >= 0;
	pass = pass && !aperture.region.outer_loop.empty() && !aperture.region.convex_pieces.empty();
	pass = pass && close(aperture.region.convex_piece_area_sum, aperture.area);
	if (side >= 0 && side < 6) boundary_area[side] += aperture.area;
  }
  for (double area : boundary_area) pass = pass && close(area, 1.0);
  pass = pass && std::abs(result.relative_volume_conservation_residual) < 1.0e-10;
  pass = pass && result.first_moment_conservation_residual < 1.0e-10;
  pass = pass && result.requested_fuzzy_tolerance == 0.0;

  std::cout << '[' << name << "] triangles=" << input.triangles.size()
            << " fragments=" << result.fragments.size()
            << " patches=" << result.surface_patches.size()
            << " apertures=" << result.box_apertures.size()
            << " surface-area=" << std::setprecision(12) << patch_area
            << " volume-residual=" << result.volume_conservation_residual
            << " requested-fuzzy=" << result.requested_fuzzy_tolerance
            << " effective-fuzzy=" << result.effective_fuzzy_tolerance
	    << " time-ms=" << result.general_fuse_milliseconds
            << " status=" << (pass ? "PASS" : "FAIL") << '\n';
  for (const std::string &warning : result.warnings)
	std::cout << "  warning: " << warning << '\n';
  for (const std::string &error : result.errors)
	std::cout << "  error: " << error << '\n';
  if (!pass)
  {
	std::cout << "  volumes=";
	for (double volume : volumes) std::cout << ' ' << volume;
	std::cout << " expected=";
	for (double volume : wanted) std::cout << ' ' << volume;
	std::cout << " boundary-areas=";
	for (double area : boundary_area) std::cout << ' ' << area;
	std::cout << '\n';
  }
  return pass;
}

ExactCellInput membrane(int nu, int nv)
{
  InputBuilder builder;
  builder.grid({ 0.5, -0.1, -0.1 }, { 0.0, 1.2, 0.0 }, { 0.0, 0.0, 1.2 },
      nu, nv, 10);
  return std::move(builder.input);
}

ExactCellInput finite_patch()
{
  InputBuilder builder;
  builder.grid({ 0.5, 0.3, 0.3 }, { 0.0, 0.4, 0.0 }, { 0.0, 0.0, 0.4 }, 1, 1, 20);
  return std::move(builder.input);
}

ExactCellInput gap()
{
  InputBuilder builder;
  builder.grid({ 0.5, -0.1, -0.1 }, { 0.0, 1.2, 0.0 }, { 0.0, 0.0, 0.55 }, 1, 1, 30);
  builder.grid({ 0.5, -0.1, 0.55 }, { 0.0, 1.2, 0.0 }, { 0.0, 0.0, 0.55 }, 1, 1, 31);
  return std::move(builder.input);
}

ExactCellInput tee()
{
  InputBuilder builder;
  builder.grid({ 0.5, -0.1, -0.1 }, { 0.0, 1.2, 0.0 }, { 0.0, 0.0, 1.2 }, 1, 1, 40);
  builder.grid({ 0.5, 0.5, -0.1 }, { 0.6, 0.0, 0.0 }, { 0.0, 0.0, 1.2 }, 1, 1, 41);
  return std::move(builder.input);
}

ExactCellInput cross_junction()
{
  InputBuilder builder;
  builder.grid({ 0.5, -0.1, -0.1 }, { 0.0, 1.2, 0.0 }, { 0.0, 0.0, 1.2 }, 1, 1, 50);
  builder.grid({ -0.1, 0.5, -0.1 }, { 1.2, 0.0, 0.0 }, { 0.0, 0.0, 1.2 }, 1, 1, 51);
  return std::move(builder.input);
}

ExactCellInput inclined_membrane(int nu = 1, int nv = 1)
{
  InputBuilder builder;
  // x = 0.25 + 0.5 y.  The clipped sheet has area sqrt(1.25), and
  // its two fluid regions each have volume 0.5.
  builder.grid({ 0.25, 0.0, 0.0 }, { 0.5, 1.0, 0.0 }, { 0.0, 0.0, 1.0 },
      nu, nv, 70);
  return std::move(builder.input);
}

ExactCellInput reverse_winding(ExactCellInput input)
{
  for (ExactCellTriangle &triangle : input.triangles)
	std::swap(triangle.vertices[1], triangle.vertices[2]);
  return input;
}

ExactCellInput boundary_fabric(double x, double y0, double y1,
    double z0 = 0.0, double z1 = 1.0, std::uint64_t source_face_id = 80)
{
  InputBuilder builder;
  builder.grid({ x, y0, z0 }, { 0.0, y1 - y0, 0.0 },
      { 0.0, 0.0, z1 - z0 }, 1, 1, source_face_id);
  return std::move(builder.input);
}

ExactCellInput boundary_triangle_with_open_surround()
{
  InputBuilder builder;
  builder.input.vertices = {
      { 0.0, 0.30, 0.30 }, { 0.0, 0.70, 0.30 }, { 0.0, 0.50, 0.70 } };
  builder.input.triangles.push_back({ { 0, 1, 2 }, 1, 90 });
  return std::move(builder.input);
}

ExactCellInput duplicate_triangles()
{
  InputBuilder builder;
  builder.input.vertices = {
      { 0.5, 0.1, 0.1 }, { 0.5, 0.9, 0.1 }, { 0.5, 0.1, 0.9 } };
  builder.input.triangles.push_back({ { 0, 1, 2 }, 1, 100 });
  builder.input.triangles.push_back({ { 0, 1, 2 }, 2, 101 });
  return std::move(builder.input);
}

ExactCellInput overlapping_triangles()
{
  InputBuilder builder;
  builder.input.vertices = {
      { 0.5, 0.1, 0.1 }, { 0.5, 0.9, 0.1 }, { 0.5, 0.1, 0.9 },
      { 0.5, 0.2, 0.2 }, { 0.5, 0.7, 0.2 }, { 0.5, 0.2, 0.7 } };
  builder.input.triangles.push_back({ { 0, 1, 2 }, 1, 110 });
  builder.input.triangles.push_back({ { 3, 4, 5 }, 2, 111 });
  return std::move(builder.input);
}

ExactCellInput near_confusion_sheets()
{
  InputBuilder builder;
  builder.grid({ 0.5, 0.0, 0.0 }, { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 1.0 }, 1, 1, 120);
  // OCCT reports an effective Precision::Confusion() floor of 1e-7 m here.
  // Two physical sheets only 2.5e-8 m apart are therefore ambiguous and must
  // be rejected, never silently collapsed into one membrane.
  builder.grid({ 0.500000025, 0.0, 0.0 }, { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 1.0 }, 1, 1, 121);
  return std::move(builder.input);
}

ExactCellInput neighbour_inclined_sheet(bool upper_cell, int nu, int nv)
{
  InputBuilder builder;
  builder.input.cell_min = { upper_cell ? 1.0 : 0.0, 0.0, 0.0 };
  builder.input.cell_max = { upper_cell ? 2.0 : 1.0, 1.0, 1.0 };
  builder.grid({ upper_cell ? 1.0 : 0.0, 0.0, 0.25 },
      { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.5 }, nu, nv,
      upper_cell ? 131 : 130);
  return std::move(builder.input);
}

void print_result_diagnostics(const ExactCellDecomposition &result)
{
  for (const std::string &warning : result.warnings)
	std::cout << "  warning: " << warning << '\n';
  for (const std::string &error : result.errors)
	std::cout << "  error: " << error << '\n';
}

bool report_custom_case(const std::string &name,
    const ExactCellDecomposition &result, bool pass)
{
  std::cout << '[' << name << "] fragments=" << result.fragments.size()
            << " patches=" << result.surface_patches.size()
            << " apertures=" << result.box_apertures.size()
            << " valid=" << (result.valid() ? "yes" : "no")
            << " status=" << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass) print_result_diagnostics(result);
  return pass;
}

ExactCellPoint patch_pressure_force(const ExactCellDecomposition &result)
{
  ExactCellPoint force{};
  auto pressure = [&](std::int32_t fragment)
  {
	if (fragment < 0 || fragment >= static_cast<std::int32_t>(result.fragments.size()))
	  return 0.0;
	const ExactCellPoint &centroid = result.fragments[static_cast<std::size_t>(fragment)].centroid;
	const double signed_distance = centroid[0] - 0.25 - 0.5 * centroid[1];
	return signed_distance < 0.0 ? 3.0 : 1.0;
  };
  for (const auto &patch : result.surface_patches)
  {
	const double delta_p = pressure(patch.minus_fragment) - pressure(patch.plus_fragment);
	for (int axis = 0; axis < 3; ++axis)
	  force[axis] += delta_p * patch.area * patch.normal[axis];
  }
  return force;
}

bool check_reversed_winding()
{
  const ExactCellDecomposition forward =
      paracfd::core::decompose_exact_cell(inclined_membrane());
  const ExactCellDecomposition reversed =
      paracfd::core::decompose_exact_cell(reverse_winding(inclined_membrane()));
  const ExactCellPoint forward_force = patch_pressure_force(forward);
  const ExactCellPoint reversed_force = patch_pressure_force(reversed);
  bool side_swap = true;
  auto signed_side = [](const ExactCellPoint &point)
  { return point[0] - 0.25 - 0.5 * point[1]; };
  for (const auto &patch : forward.surface_patches)
	side_swap = side_swap && patch.plus_fragment >= 0 && patch.minus_fragment >= 0
	    && signed_side(forward.fragments[static_cast<std::size_t>(patch.plus_fragment)].centroid) > 0.0
	    && signed_side(forward.fragments[static_cast<std::size_t>(patch.minus_fragment)].centroid) < 0.0;
  for (const auto &patch : reversed.surface_patches)
	side_swap = side_swap && patch.plus_fragment >= 0 && patch.minus_fragment >= 0
	    && signed_side(reversed.fragments[static_cast<std::size_t>(patch.plus_fragment)].centroid) < 0.0
	    && signed_side(reversed.fragments[static_cast<std::size_t>(patch.minus_fragment)].centroid) > 0.0;
  double forward_area = 0.0, reversed_area = 0.0;
  for (const auto &patch : forward.surface_patches) forward_area += patch.area;
  for (const auto &patch : reversed.surface_patches) reversed_area += patch.area;
  bool pass = forward.valid() && reversed.valid() && side_swap
      && forward.fragments.size() == 2 && reversed.fragments.size() == 2
      && close(forward_area, std::sqrt(1.25)) && close(reversed_area, forward_area)
      && close(forward_force[0], 2.0) && close(forward_force[1], -1.0)
      && close(forward_force[2], 0.0)
      && close(reversed_force[0], forward_force[0])
      && close(reversed_force[1], forward_force[1])
      && close(reversed_force[2], forward_force[2]);
  std::cout << "[reversed-winding] F=" << forward_force[0] << ','
            << forward_force[1] << ',' << forward_force[2]
            << " reversed-F=" << reversed_force[0] << ',' << reversed_force[1]
            << ',' << reversed_force[2] << " status=" << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass)
  {
	std::cout << "  forward diagnostics:\n";
	print_result_diagnostics(forward);
	std::cout << "  reversed diagnostics:\n";
	print_result_diagnostics(reversed);
  }
  return pass;
}

struct SideMeasure
{
  double area = 0.0;
  ExactCellPoint first_moment{};
  std::vector<std::array<double, 4>> pieces;
};

SideMeasure box_side_measure(const ExactCellDecomposition &result, int axis, bool upper)
{
  SideMeasure measure;
  for (const auto &aperture : result.box_apertures)
	if (aperture.axis == axis && aperture.upper == upper)
	{
	  measure.area += aperture.area;
	  for (int q = 0; q < 3; ++q)
		measure.first_moment[q] += aperture.area * aperture.centroid[q];
	  measure.pieces.push_back({ aperture.area, aperture.centroid[0],
	      aperture.centroid[1], aperture.centroid[2] });
	}
  std::sort(measure.pieces.begin(), measure.pieces.end());
  return measure;
}

bool check_boundary_membrane()
{
  const ExactCellDecomposition result =
      paracfd::core::decompose_exact_cell(boundary_fabric(1.0, 0.0, 1.0));
  bool patch_sides = !result.surface_patches.empty();
  double area = 0.0;
  for (const auto &patch : result.surface_patches)
  {
	area += patch.area;
	patch_sides = patch_sides && patch.boundary_axis == 0 && patch.boundary_upper
	    && patch.minus_fragment >= 0 && patch.plus_fragment < 0;
  }
  const SideMeasure blocked_side = box_side_measure(result, 0, true);
  const SideMeasure opposite_side = box_side_measure(result, 0, false);
  const bool pass = result.valid() && result.fragments.size() == 1 && patch_sides
      && close(area, 1.0) && close(blocked_side.area, 0.0)
      && close(opposite_side.area, 1.0);
  return report_custom_case("boundary-coincident-membrane", result, pass);
}

bool check_planar_region_hole()
{
  const ExactCellDecomposition result = paracfd::core::decompose_exact_cell(
      boundary_triangle_with_open_surround());
  bool found_hole = false;
  for (const auto &aperture : result.box_apertures)
	if (aperture.axis == 0 && !aperture.upper)
	  found_hole = found_hole || !aperture.region.hole_loops.empty();
  const SideMeasure side = box_side_measure(result, 0, false);
  const double triangle_area = 0.08;
  const bool pass = result.valid() && found_hole && result.fragments.size() == 1
      && close(side.area, 1.0 - triangle_area)
      && close(side.first_moment[0], 0.0)
      && close(side.first_moment[1], 0.5 - triangle_area * 0.5)
      && close(side.first_moment[2], 0.5 - triangle_area * (1.3 / 3.0));
  return report_custom_case("planar-region-hole", result, pass);
}

bool check_rejected_case(const std::string &name, const ExactCellInput &input)
{
  const ExactCellDecomposition result = paracfd::core::decompose_exact_cell(input);
  return report_custom_case(name, result, !result.valid() && !result.errors.empty());
}

std::int64_t quantize(double value)
{
  return static_cast<std::int64_t>(std::llround(value * 1.0e10));
}

std::string aggregate_signature(const ExactCellDecomposition &result)
{
  std::vector<std::string> records;
  for (const auto &fragment : result.fragments)
  {
	std::ostringstream record;
	record << "F:" << quantize(fragment.volume);
	for (double coordinate : fragment.centroid) record << ':' << quantize(coordinate);
	records.push_back(record.str());
  }
  auto fragment_signature = [&](std::int32_t fragment)
  {
	if (fragment < 0) return std::string("outside");
	const auto &value = result.fragments[static_cast<std::size_t>(fragment)];
	std::ostringstream record;
	record << quantize(value.volume);
	for (double coordinate : value.centroid) record << ':' << quantize(coordinate);
	return record.str();
  };
  for (const auto &patch : result.surface_patches)
  {
	std::ostringstream record;
	record << "P:" << patch.source_triangle_id << ':' << patch.source_face_id
	       << ':' << quantize(patch.area);
	for (double coordinate : patch.centroid) record << ':' << quantize(coordinate);
	for (double coordinate : patch.normal) record << ':' << quantize(coordinate);
	record << ':' << fragment_signature(patch.plus_fragment)
	       << ':' << fragment_signature(patch.minus_fragment);
	records.push_back(record.str());
  }
  for (const auto &aperture : result.box_apertures)
  {
	std::ostringstream record;
	record << "A:" << static_cast<int>(aperture.axis) << ':' << aperture.upper
	       << ':' << quantize(aperture.area);
	for (double coordinate : aperture.centroid) record << ':' << quantize(coordinate);
	record << ':' << fragment_signature(aperture.fragment);
	records.push_back(record.str());
  }
  std::sort(records.begin(), records.end());
  std::ostringstream signature;
  for (const std::string &record : records) signature << record << '\n';
  return signature.str();
}

bool check_triangle_permutation()
{
  ExactCellInput original = tee();
  ExactCellInput permuted = original;
  std::reverse(permuted.triangles.begin(), permuted.triangles.end());
  const ExactCellDecomposition a = paracfd::core::decompose_exact_cell(original);
  const ExactCellDecomposition b = paracfd::core::decompose_exact_cell(permuted);
  const bool pass = a.valid() && b.valid()
      && aggregate_signature(a) == aggregate_signature(b);
  std::cout << "[triangle-order-permutation] signature-bytes="
            << aggregate_signature(a).size() << " status="
            << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass)
  {
	std::cout << "  original signature:\n" << aggregate_signature(a)
	          << "  permuted signature:\n" << aggregate_signature(b);
	print_result_diagnostics(a);
	print_result_diagnostics(b);
  }
  return pass;
}

bool check_neighbour_partition_match()
{
  const ExactCellDecomposition lower = paracfd::core::decompose_exact_cell(
      neighbour_inclined_sheet(false, 1, 1));
  const ExactCellDecomposition upper = paracfd::core::decompose_exact_cell(
      neighbour_inclined_sheet(true, 3, 2));
  const SideMeasure a = box_side_measure(lower, 0, true);
  const SideMeasure b = box_side_measure(upper, 0, false);
  bool pass = lower.valid() && upper.valid() && close(a.area, 1.0)
      && close(b.area, 1.0) && close(a.area, b.area);
  for (int axis = 0; axis < 3; ++axis)
	pass = pass && close(a.first_moment[axis], b.first_moment[axis])
	    && close(a.first_moment[axis], 0.5 * (axis == 0 ? 2.0 : 1.0));
  pass = pass && a.pieces.size() == b.pieces.size();
  for (std::size_t piece = 0; piece < std::min(a.pieces.size(), b.pieces.size()); ++piece)
	for (int value = 0; value < 4; ++value)
	  pass = pass && close(a.pieces[piece][value], b.pieces[piece][value]);
  std::cout << "[neighbour-mismatched-triangulation] lower-pieces=" << a.pieces.size()
            << " upper-pieces=" << b.pieces.size() << " area=" << a.area
            << '/' << b.area << " status=" << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass)
  {
	std::cout << "  lower diagnostics:\n";
	print_result_diagnostics(lower);
	std::cout << "  upper diagnostics:\n";
	print_result_diagnostics(upper);
  }
  return pass;
}

TriMesh mesh_from_input(const ExactCellInput &input)
{
  TriMesh mesh;
  mesh.positions_fp64.reserve(3 * input.vertices.size());
  mesh.positions.reserve(3 * input.vertices.size());
  for (const ExactCellPoint &point : input.vertices)
  {
	for (double coordinate : point)
	{
	  mesh.positions_fp64.push_back(coordinate);
	  mesh.positions.push_back(static_cast<float>(coordinate));
	}
  }
  for (const ExactCellTriangle &triangle : input.triangles)
  {
	mesh.indices.insert(mesh.indices.end(), triangle.vertices.begin(), triangle.vertices.end());
	mesh.source_face_ids.push_back(static_cast<std::uint32_t>(triangle.source_face_id));
  }
  if (!input.vertices.empty())
  {
	for (int axis = 0; axis < 3; ++axis)
	{
	  double minimum = input.vertices.front()[axis], maximum = minimum;
	  for (const ExactCellPoint &point : input.vertices)
	  {
		minimum = std::min(minimum, point[axis]);
		maximum = std::max(maximum, point[axis]);
	  }
	  mesh.bbox_min[axis] = static_cast<float>(minimum);
	  mesh.bbox_max[axis] = static_cast<float>(maximum);
	}
  }
  return mesh;
}

TriMesh certified_mismatched_tee_mesh()
{
  InputBuilder builder;
  builder.grid({ 0.5, 0.0, 0.0 }, { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 1.0 }, 1, 1, 140);
  constexpr double gap = 5.0e-5;
  builder.grid({ 0.5 + gap, 0.5, 0.0 }, { 0.5 - gap, 0.0, 0.0 },
      { 0.0, 0.0, 1.0 }, 1, 1, 141);
  TriMesh mesh = mesh_from_input(builder.input);
  const std::size_t half_edges = mesh.indices.size();
  mesh.triangle_cad_edge_provenance_states.assign(half_edges,
      static_cast<std::uint8_t>(paracfd::core::CadEdgeProvenanceState::none));
  mesh.triangle_cad_edge_ids.assign(half_edges, TriMesh::kNoCadEdgeId);
  mesh.triangle_cad_edge_incident_face_counts.assign(half_edges, 0);
  mesh.triangle_cad_edge_tolerances.assign(half_edges, 0.0);
  mesh.triangle_cad_edge_is_periodic_seam.assign(half_edges, 0);
  mesh.triangle_cad_edge_contact_ids.assign(half_edges, TriMesh::kNoCadContactId);
  mesh.triangle_cad_edge_certified_fan_degrees.assign(half_edges, 0);

  // The rib's u=0 boundary is triangle 3, local edge 2 (p01 -> p00).
  // CAD certifies that this curve belongs to a skin/rib fan of degree three,
  // even though independent tessellation leaves the polygons 50 microns apart.
  const std::size_t contact = 3 * 3 + 2;
  mesh.triangle_cad_edge_provenance_states[contact] =
      static_cast<std::uint8_t>(paracfd::core::CadEdgeProvenanceState::known);
  mesh.triangle_cad_edge_ids[contact] = 900;
  mesh.triangle_cad_edge_incident_face_counts[contact] = 3;
  mesh.triangle_cad_edge_tolerances[contact] = 2.0 * gap;
  mesh.triangle_cad_edge_contact_ids[contact] = 9000;
  mesh.triangle_cad_edge_certified_fan_degrees[contact] = 3;
  return mesh;
}

void initialize_contact_provenance(TriMesh &mesh)
{
  const std::size_t half_edges = mesh.indices.size();
  mesh.triangle_cad_edge_provenance_states.assign(half_edges,
      static_cast<std::uint8_t>(paracfd::core::CadEdgeProvenanceState::none));
  mesh.triangle_cad_edge_ids.assign(half_edges, TriMesh::kNoCadEdgeId);
  mesh.triangle_cad_edge_incident_face_counts.assign(half_edges, 0);
  mesh.triangle_cad_edge_tolerances.assign(half_edges, 0.0);
  mesh.triangle_cad_edge_is_periodic_seam.assign(half_edges, 0);
  mesh.triangle_cad_edge_contact_ids.assign(half_edges, TriMesh::kNoCadContactId);
  mesh.triangle_cad_edge_certified_fan_degrees.assign(half_edges, 0);
}

void certify_contact(TriMesh &mesh, std::size_t triangle, int edge,
    std::uint32_t cad_edge, std::uint64_t contact, double tolerance,
    std::uint32_t fan)
{
  const std::size_t offset = 3 * triangle + static_cast<std::size_t>(edge);
  mesh.triangle_cad_edge_provenance_states[offset] =
      static_cast<std::uint8_t>(paracfd::core::CadEdgeProvenanceState::known);
  mesh.triangle_cad_edge_ids[offset] = cad_edge;
  mesh.triangle_cad_edge_incident_face_counts[offset] = fan;
  mesh.triangle_cad_edge_tolerances[offset] = tolerance;
  mesh.triangle_cad_edge_contact_ids[offset] = contact;
  mesh.triangle_cad_edge_certified_fan_degrees[offset] = fan;
}

TriMesh certified_conforming_tee_mesh()
{
  TriMesh mesh = mesh_from_input(tee());
  initialize_contact_provenance(mesh);
  // The second grid's u=0 boundary is triangle 3, local edge 2.
  certify_contact(mesh, 3, 2, 901, 9001, 1.0e-10, 3);
  return mesh;
}

TriMesh certified_cell_face_mismatched_tee_mesh()
{
  InputBuilder builder;
  // The skin crosses both cells. The rib begins on their shared x=1 face,
  // but is displaced 50 microns from the skin in y.
  builder.grid({ 0.0, 0.5, 0.0 }, { 2.0, 0.0, 0.0 },
      { 0.0, 0.0, 1.0 }, 2, 1, 142);
  constexpr double gap = 5.0e-5;
  builder.grid({ 1.0, 0.5 + gap, 0.0 }, { 0.5, 0.0, 0.0 },
      { 0.0, 0.0, 1.0 }, 1, 1, 143);
  TriMesh mesh = mesh_from_input(builder.input);
  initialize_contact_provenance(mesh);
  // Four skin triangles precede the two rib triangles. The rib u=0 edge is
  // triangle 5, local edge 2 and lies exactly on the Cartesian cell face.
  certify_contact(mesh, 5, 2, 902, 9002, 2.0 * gap, 3);
  return mesh;
}

EmbeddedBoundary build_exact_eb(const ExactCellInput &input,
    UniformEbGrid grid = { { 0.0, 0.0, 0.0 }, 1, 1, 1, 1.0 })
{
  TriMesh mesh = mesh_from_input(input);
  TriangleBvh bvh(mesh);
  EmbeddedBoundaryBuildOptions options;
  options.exact_cell_decomposer = &paracfd::core::decompose_exact_cell;
  return paracfd::core::build_embedded_boundary(mesh, bvh, grid, options);
}

ExactCellDecomposition disagreeing_shared_face_decomposer(const ExactCellInput &input)
{
  ExactCellDecomposition result = paracfd::core::decompose_exact_cell(input);
  // Deliberately corrupt only the right cell's ownership of its -X apertures.
  // Geometry and aperture measure remain identical, so an adapter that silently
  // seals centroid-incompatible carriers would incorrectly report success.
  if (result.valid() && input.cell_min[0] > 0.5 && result.fragments.size() == 2)
    for (auto &aperture : result.box_apertures)
      if (aperture.axis == 0 && !aperture.upper)
        aperture.fragment = 1 - aperture.fragment;
  return result;
}

bool check_shared_face_fragment_disagreement()
{
  InputBuilder builder;
  builder.grid({ -0.1, 0.0, 0.5 }, { 2.2, 0.0, 0.0 },
      { 0.0, 1.0, 0.0 }, 2, 1, 175);
  TriMesh mesh = mesh_from_input(builder.input);
  TriangleBvh bvh(mesh);
  EmbeddedBoundaryBuildOptions options;
  options.exact_cell_decomposer = &disagreeing_shared_face_decomposer;
  const EmbeddedBoundary boundary = paracfd::core::build_embedded_boundary(mesh, bvh,
      { { 0.0, 0.0, 0.0 }, 2, 1, 1, 1.0 }, options);
  bool named_disagreement = false;
  for (const auto &unresolved : boundary.unresolved)
    named_disagreement = named_disagreement
        || unresolved.reason.find("ownership disagrees") != std::string::npos
        || unresolved.reason.find("partitions disagree") != std::string::npos;
  const bool pass = !boundary.ready_for_flow() && boundary.unresolved.size() == 2
      && named_disagreement && boundary.rejected_cross_fabric_apertures == 0
      && close(boundary.rejected_cross_fabric_aperture_area, 0.0);
  std::cout << "[shared-face-fragment-disagreement] ready="
            << (boundary.ready_for_flow() ? "yes" : "no")
            << " unresolved=" << boundary.unresolved.size() << " status="
            << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass)
    for (const auto &unresolved : boundary.unresolved)
      std::cout << "  unresolved cell " << unresolved.parent_cell << ": "
                << unresolved.reason << '\n';
  return pass;
}

double eb_volume(const EmbeddedBoundary &boundary)
{
  return std::accumulate(boundary.fragments.begin(), boundary.fragments.end(), 0.0,
      [](double sum, const auto &fragment) { return sum + fragment.volume; });
}

double eb_patch_area(const EmbeddedBoundary &boundary)
{
  return std::accumulate(boundary.patches.begin(), boundary.patches.end(), 0.0,
      [](double sum, const auto &patch) { return sum + patch.area; });
}

bool check_shared_face_fabric_adapter(bool partial)
{
  const ExactCellInput input = boundary_fabric(1.0, 0.0, partial ? 0.4 : 1.0,
      0.0, 1.0, partial ? 151 : 150);
  const EmbeddedBoundary boundary = build_exact_eb(input,
      { { 0.0, 0.0, 0.0 }, 2, 1, 1, 1.0 });
  double open_area = 0.0;
  Vec3d open_first_moment{};
  int open_pieces = 0;
  bool valid_refs = true;
  for (const auto &aperture : boundary.apertures)
	if (aperture.parent_face_cell == 0 && aperture.axis == 0)
	{
	  open_area += aperture.area;
	  open_first_moment = open_first_moment + aperture.centroid * aperture.area;
	  ++open_pieces;
	  valid_refs = valid_refs && aperture.fragment_a != paracfd::core::invalid_fragment
	      && aperture.fragment_b != paracfd::core::invalid_fragment;
	}
  bool patch_sides = true;
  for (const auto &patch : boundary.patches)
	patch_sides = patch_sides && patch.plus_fragment != paracfd::core::invalid_fragment
	    && patch.minus_fragment != paracfd::core::invalid_fragment
	    && patch.plus_fragment != patch.minus_fragment;
  const double expected_open = partial ? 0.6 : 0.0;
  const double expected_patch = partial ? 0.4 : 1.0;
  bool pass = boundary.ready_for_flow() && valid_refs && patch_sides
      && close(open_area, expected_open) && close(eb_patch_area(boundary), expected_patch);
  if (partial)
  {
	pass = pass && open_pieces > 0 && close(open_first_moment.x, 0.6)
	    && close(open_first_moment.y, 0.42) && close(open_first_moment.z, 0.30);
  }
  else
  {
	pass = pass && open_pieces == 0 && (boundary.cut_face_mask[0] & 1u) != 0;
  }
  std::cout << '[' << (partial ? "shared-face-partial-fabric" : "shared-face-full-fabric")
            << "] open-pieces=" << open_pieces << " open-area=" << open_area
            << " patch-area=" << eb_patch_area(boundary)
            << " unresolved=" << boundary.unresolved.size() << " status="
            << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass)
	for (const auto &unresolved : boundary.unresolved)
	  std::cout << "  unresolved cell " << unresolved.parent_cell << ": "
	            << unresolved.reason << '\n';
  return pass;
}

bool check_certified_mismatched_tee_adapter()
{
  TriMesh mesh = certified_mismatched_tee_mesh();
  TriangleBvh bvh(mesh);
  EmbeddedBoundaryBuildOptions options;
  options.exact_cell_decomposer = &paracfd::core::decompose_exact_cell;
  const EmbeddedBoundary boundary = paracfd::core::build_embedded_boundary(mesh, bvh,
      { { 0.0, 0.0, 0.0 }, 1, 1, 1, 1.0 }, options);
  bool rib_reconnects = false;
  for (const auto &patch : boundary.patches)
	if (patch.source_face_id == 141)
	  rib_reconnects = rib_reconnects || patch.plus_fragment == patch.minus_fragment;
  // With only a tessellated gap, the exact BRep cannot legitimately invent the
  // missing attachment.  A certified fan-3 contact therefore must either yield
  // all three sectors or fail the cell closed for refinement/CAD repair.
  const bool reconstructed = boundary.ready_for_flow()
      && boundary.fragments.size() >= 3 && !rib_reconnects;
  const bool rejected = !boundary.ready_for_flow() && !boundary.unresolved.empty();
  const bool pass = reconstructed || rejected;
  std::cout << "[eb-certified-mismatched-T] ready="
            << (boundary.ready_for_flow() ? "yes" : "no")
            << " fragments=" << boundary.fragments.size()
            << " rib-reconnects=" << (rib_reconnects ? "yes" : "no")
            << " status=" << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass || rejected)
	for (const auto &unresolved : boundary.unresolved)
	  std::cout << "  unresolved cell " << unresolved.parent_cell << ": "
	            << unresolved.reason << '\n';
  return pass;
}

bool check_certified_conforming_tee_adapter()
{
  TriMesh mesh = certified_conforming_tee_mesh();
  TriangleBvh bvh(mesh);
  EmbeddedBoundaryBuildOptions options;
  options.exact_cell_decomposer = &paracfd::core::decompose_exact_cell;
  const EmbeddedBoundary boundary = paracfd::core::build_embedded_boundary(mesh, bvh,
      { { 0.0, 0.0, 0.0 }, 1, 1, 1, 1.0 }, options);
  const bool pass = boundary.ready_for_flow() && boundary.fragments.size() == 3;
  std::cout << "[eb-certified-conforming-T] ready="
            << (boundary.ready_for_flow() ? "yes" : "no")
            << " fragments=" << boundary.fragments.size() << " status="
            << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass)
    for (const auto &unresolved : boundary.unresolved)
      std::cout << "  unresolved cell " << unresolved.parent_cell << ": "
                << unresolved.reason << '\n';
  return pass;
}

bool check_certified_cell_face_mismatched_tee_adapter()
{
  TriMesh mesh = certified_cell_face_mismatched_tee_mesh();
  TriangleBvh bvh(mesh);
  EmbeddedBoundaryBuildOptions options;
  options.exact_cell_decomposer = &paracfd::core::decompose_exact_cell;
  const EmbeddedBoundary boundary = paracfd::core::build_embedded_boundary(mesh, bvh,
      { { 0.0, 0.0, 0.0 }, 2, 1, 1, 1.0 }, options);
  bool names_contact = false;
  for (const auto &unresolved : boundary.unresolved)
    names_contact = names_contact || unresolved.reason.find("CAD contact 9002")
        != std::string::npos;
  const bool pass = !boundary.ready_for_flow() && names_contact;
  std::cout << "[eb-certified-cell-face-mismatched-T] ready="
            << (boundary.ready_for_flow() ? "yes" : "no")
            << " unresolved=" << boundary.unresolved.size() << " status="
            << (pass ? "PASS" : "FAIL") << '\n';
  if (!pass)
    for (const auto &unresolved : boundary.unresolved)
      std::cout << "  unresolved cell " << unresolved.parent_cell << ": "
                << unresolved.reason << '\n';
  return pass;
}

bool check_embedded_boundary_integration()
{
  bool pass = true;
  const EmbeddedBoundary membrane_eb = build_exact_eb(membrane(1, 1));
  bool membrane_sides = !membrane_eb.patches.empty();
  for (const auto &patch : membrane_eb.patches)
	membrane_sides = membrane_sides && patch.plus_fragment != patch.minus_fragment;
  const FragmentRef membrane_minus = membrane_eb.fragment_containing_point(0, { 0.25, 0.5, 0.5 });
  const FragmentRef membrane_plus = membrane_eb.fragment_containing_point(0, { 0.75, 0.5, 0.5 });
  const bool membrane_pass = membrane_eb.ready_for_flow() && membrane_eb.exact_decomposition_cells == 1
      && membrane_eb.exact_locators.size() == 1 && membrane_eb.fragments.size() == 2
      && close(eb_volume(membrane_eb), 1.0) && close(eb_patch_area(membrane_eb), 1.0)
      && membrane_sides && membrane_minus != paracfd::core::invalid_fragment
      && membrane_plus != paracfd::core::invalid_fragment && membrane_minus != membrane_plus;
  std::cout << "[eb-integration-membrane] fragments=" << membrane_eb.fragments.size()
            << " patches=" << membrane_eb.patches.size() << " status="
            << (membrane_pass ? "PASS" : "FAIL") << '\n';
  pass = membrane_pass && pass;

  const EmbeddedBoundary opening_eb = build_exact_eb(gap());
  const FragmentRef opening_left = opening_eb.fragment_containing_point(0, { 0.25, 0.5, 0.5 });
  const FragmentRef opening_right = opening_eb.fragment_containing_point(0, { 0.75, 0.5, 0.5 });
  const bool opening_pass = opening_eb.ready_for_flow() && opening_eb.fragments.size() == 1
      && close(eb_volume(opening_eb), 1.0) && close(eb_patch_area(opening_eb), 0.9)
      && opening_left != paracfd::core::invalid_fragment && opening_left == opening_right;
  std::cout << "[eb-integration-opening] fragments=" << opening_eb.fragments.size()
            << " area=" << eb_patch_area(opening_eb) << " status="
            << (opening_pass ? "PASS" : "FAIL") << '\n';
  pass = opening_pass && pass;

  const EmbeddedBoundary tee_eb = build_exact_eb(tee());
  const std::array<FragmentRef, 4> tee_samples{
      tee_eb.fragment_containing_point(0, { 0.25, 0.25, 0.5 }),
      tee_eb.fragment_containing_point(0, { 0.75, 0.25, 0.5 }),
      tee_eb.fragment_containing_point(0, { 0.25, 0.75, 0.5 }),
      tee_eb.fragment_containing_point(0, { 0.75, 0.75, 0.5 }) };
  const bool tee_pass = tee_eb.ready_for_flow() && tee_eb.fragments.size() == 3
      && close(eb_volume(tee_eb), 1.0) && close(eb_patch_area(tee_eb), 1.5)
      && tee_samples[0] == tee_samples[2] && tee_samples[0] != tee_samples[1]
      && tee_samples[0] != tee_samples[3] && tee_samples[1] != tee_samples[3];
  std::cout << "[eb-integration-T] fragments=" << tee_eb.fragments.size()
            << " status=" << (tee_pass ? "PASS" : "FAIL") << '\n';
  pass = tee_pass && pass;

  const EmbeddedBoundary cross_eb = build_exact_eb(cross_junction());
  std::array<FragmentRef, 4> cross_samples{
      cross_eb.fragment_containing_point(0, { 0.25, 0.25, 0.5 }),
      cross_eb.fragment_containing_point(0, { 0.75, 0.25, 0.5 }),
      cross_eb.fragment_containing_point(0, { 0.25, 0.75, 0.5 }),
      cross_eb.fragment_containing_point(0, { 0.75, 0.75, 0.5 }) };
  std::sort(cross_samples.begin(), cross_samples.end());
  const bool cross_pass = cross_eb.ready_for_flow() && cross_eb.fragments.size() == 4
      && close(eb_volume(cross_eb), 1.0) && close(eb_patch_area(cross_eb), 2.0)
      && std::adjacent_find(cross_samples.begin(), cross_samples.end()) == cross_samples.end();
  std::cout << "[eb-integration-cross] fragments=" << cross_eb.fragments.size()
            << " status=" << (cross_pass ? "PASS" : "FAIL") << '\n';
  pass = cross_pass && pass;

  InputBuilder shared;
  shared.grid({ -0.1, 0.0, 0.5 }, { 2.2, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, 2, 1, 60);
  const EmbeddedBoundary shared_eb = build_exact_eb(shared.input,
      { { 0.0, 0.0, 0.0 }, 2, 1, 1, 1.0 });
  double shared_area = 0.0;
  int shared_count = 0;
  bool shared_refs = true;
  for (const auto &aperture : shared_eb.apertures)
	if (aperture.parent_face_cell == 0 && aperture.axis == 0)
	{
	  shared_area += aperture.area;
	  ++shared_count;
	  Vec3d lower = aperture.centroid, upper = aperture.centroid;
	  lower.x -= 1.0e-5;
	  upper.x += 1.0e-5;
	  shared_refs = shared_refs
	      && aperture.fragment_a == shared_eb.fragment_containing_point(0, lower)
	      && aperture.fragment_b == shared_eb.fragment_containing_point(1, upper);
	}
  const bool conservation_pass = shared_eb.ready_for_flow() && shared_eb.fragments.size() == 4
      && close(eb_volume(shared_eb), 2.0) && close(eb_patch_area(shared_eb), 2.0)
      && shared_count >= 2 && close(shared_area, 1.0) && shared_refs;
  std::cout << "[eb-integration-shared-face] apertures=" << shared_count
            << " area=" << shared_area << " status="
            << (conservation_pass ? "PASS" : "FAIL") << '\n';
  pass = conservation_pass && pass;
  pass = check_shared_face_fabric_adapter(false) && pass;
  pass = check_shared_face_fabric_adapter(true) && pass;
  pass = check_shared_face_fragment_disagreement() && pass;
  pass = check_certified_conforming_tee_adapter() && pass;
  pass = check_certified_mismatched_tee_adapter() && pass;
  pass = check_certified_cell_face_mismatched_tee_adapter() && pass;
  return pass;
}
} // namespace

int main()
{
  bool pass = true;
  pass = check_case("empty-cell", InputBuilder{}.input, { { 1.0 }, 0.0, false, false }) && pass;
  pass = check_case("membrane", membrane(1, 1), { { 0.5, 0.5 }, 1.0, false, true }) && pass;
  pass = check_case("inclined-membrane", inclined_membrane(),
             { { 0.5, 0.5 }, std::sqrt(1.25), false, true }) && pass;
  pass = check_reversed_winding() && pass;
  pass = check_case("finite-patch", finite_patch(), { { 1.0 }, 0.16, true, false }) && pass;
  pass = check_case("finite-gap", gap(), { { 1.0 }, 0.9, true, false }) && pass;
  pass = check_case("T-junction", tee(), { { 0.25, 0.25, 0.5 }, 1.5, false, true }) && pass;
  pass = check_case("cross-junction", cross_junction(),
             { { 0.25, 0.25, 0.25, 0.25 }, 2.0, false, true }) &&
      pass;
  pass = check_case("membrane-12", membrane(2, 3), { { 0.5, 0.5 }, 1.0, false, true }) && pass;
  pass = check_case("membrane-50", membrane(5, 5), { { 0.5, 0.5 }, 1.0, false, true }) && pass;
  pass = check_boundary_membrane() && pass;
  pass = check_planar_region_hole() && pass;
  pass = check_rejected_case("reject-duplicate-triangles", duplicate_triangles()) && pass;
  pass = check_rejected_case("reject-overlapping-triangles", overlapping_triangles()) && pass;
  pass = check_rejected_case("reject-near-confusion-sheets", near_confusion_sheets()) && pass;
  pass = check_triangle_permutation() && pass;
  pass = check_neighbour_partition_match() && pass;
  pass = check_embedded_boundary_integration() && pass;
  std::cout << "[occt-cell-split-summary] status=" << (pass ? "PASS" : "FAIL") << '\n';
  return pass ? 0 : 1;
}
