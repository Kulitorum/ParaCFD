#include "core/fluid/external_aero_core.h"
#include "core/geometry/mesh_clip.h"
#include "core/geometry/naca_theory.h"
#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"
#include "core/paraglider_config.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace paracfd::core;

namespace
{
constexpr double pi = 3.1415926535897932384626433832795;
constexpr double radians_per_degree = pi / 180.0;
int validation_span_layers = 0;
int validation_topology_refinement_levels = -1;
bool validation_local_amr = true;
bool validation_preprocess_only = false;
bool validation_diagnose_connectivity = false;
double validation_max_wall_seconds = 0;
int validation_max_steps = 0;
int validation_progress_steps = 500;
double validation_surface_distance = .15;
double validation_min_volume_fraction = std::numeric_limits<double>::quiet_NaN();
double validation_smagorinsky_cs = std::numeric_limits<double>::quiet_NaN();
double validation_kinematic_viscosity = std::numeric_limits<double>::quiet_NaN();
double validation_projection_tolerance = std::numeric_limits<double>::quiet_NaN();
double validation_tessellation_deflection_mm = std::numeric_limits<double>::quiet_NaN();
double validation_upstream_margin = std::numeric_limits<double>::quiet_NaN();
double validation_downstream_margin = std::numeric_limits<double>::quiet_NaN();
double validation_vertical_margin = std::numeric_limits<double>::quiet_NaN();
std::string validation_step_override;
std::string validation_experimental_benchmark;

enum class ExperimentalBenchmark
{
  none,
  nasa_tmr_0012,
  naca_tr460_cambered,
  nrel_4415,
  cambered_suite
};
struct CaseResult
{
  double mean_cl = 0, mean_pressure_cl = 0, mean_viscous_cl = 0, rms_cl = 0, final_cl = 0, circulation_cl = 0, divergence = 0, normalized_divergence = 0, rms_divergence = 0, normalized_rms_divergence = 0, residual = 0, coverage = 0, wall_seconds = 0;
  double mean_pressure_cd = 0, mean_pitching_cm = 0;
  double previous_mean_cl = 0, mean_relative_drift = std::numeric_limits<double>::infinity(), mean_window_coverage = 0;
  double max_integrated_flux_error = 0, normalized_integrated_flux_error = 0, net_integrated_flux_error = 0;
  double preprocess_seconds = 0, gpu_mib = 0, final_step_ms = 0, final_projection_ms = 0, max_recovery_cells = 0, physical_time = 0, finest_h = 0;
  double sealed_pressure_span = 0, sealed_pressure_cl_delta = 0;
  int steps = 0, iterations = 0;
  int pressure_dofs = 0;
  std::size_t gauges = 0, recovered_sides = 0, active_bricks = 0, sealed_dofs = 0;
  bool converged = false, timed_out = false;
};

struct TimedCoefficients
{
  double time = 0, cl = 0, pressure_cl = 0, viscous_cl = 0;
  double pressure_cd = 0, pitching_cm = 0;
};

struct CoefficientWindow
{
  bool ready = false;
  double mean = 0, pressure_mean = 0, viscous_mean = 0;
  double pressure_cd_mean = 0, pitching_cm_mean = 0, rms = 0, coverage = 0;
};

CoefficientWindow coefficient_window(const std::deque<TimedCoefficients> &samples, double lo, double hi)
{
  CoefficientWindow out;
  if (!(hi > lo) || samples.size() < 2) return out;
  double integral = 0, pressure_integral = 0, viscous_integral = 0;
  double pressure_cd_integral = 0, pitching_cm_integral = 0, square_integral = 0;
  int segments = 0;
  for (std::size_t q = 0; q + 1 < samples.size(); ++q)
  {
	const TimedCoefficients &a = samples[q], &b = samples[q + 1];
	if (!(b.time > a.time)) continue;
	const double left = std::max(lo, a.time), right = std::min(hi, b.time);
	if (!(right > left)) continue;
	const double ta = (left - a.time) / (b.time - a.time), tb = (right - a.time) / (b.time - a.time), dt = right - left;
	auto lerp = [](double x, double y, double t) { return x + (y - x) * t; };
	const double ca = lerp(a.cl, b.cl, ta), cb = lerp(a.cl, b.cl, tb);
	const double pa = lerp(a.pressure_cl, b.pressure_cl, ta), pb = lerp(a.pressure_cl, b.pressure_cl, tb);
	const double va = lerp(a.viscous_cl, b.viscous_cl, ta), vb = lerp(a.viscous_cl, b.viscous_cl, tb);
	const double da = lerp(a.pressure_cd, b.pressure_cd, ta), db = lerp(a.pressure_cd, b.pressure_cd, tb);
	const double ma = lerp(a.pitching_cm, b.pitching_cm, ta), mb = lerp(a.pitching_cm, b.pitching_cm, tb);
	integral += .5 * dt * (ca + cb);
	pressure_integral += .5 * dt * (pa + pb);
	viscous_integral += .5 * dt * (va + vb);
	pressure_cd_integral += .5 * dt * (da + db);
	pitching_cm_integral += .5 * dt * (ma + mb);
	square_integral += dt * (ca * ca + ca * cb + cb * cb) / 3.0;
	out.coverage += dt;
	++segments;
  }
  if (segments < 3 || out.coverage < .98 * (hi - lo)) return out;
  out.ready = true;
  out.mean = integral / out.coverage;
  out.pressure_mean = pressure_integral / out.coverage;
  out.viscous_mean = viscous_integral / out.coverage;
  out.pressure_cd_mean = pressure_cd_integral / out.coverage;
  out.pitching_cm_mean = pitching_cm_integral / out.coverage;
  out.rms = std::sqrt(std::max(0.0, square_integral / out.coverage - out.mean * out.mean));
  return out;
}
struct LiftReference
{
  double thin = 0, xfoil_inviscid = std::numeric_limits<double>::quiet_NaN(), xfoil_free_transition = std::numeric_limits<double>::quiet_NaN();
  double xfoil_forced_turbulent = std::numeric_limits<double>::quiet_NaN();
  double experimental = std::numeric_limits<double>::quiet_NaN();
  double experimental_pressure_cd = std::numeric_limits<double>::quiet_NaN();
  double experimental_pitching_cm = std::numeric_limits<double>::quiet_NaN();
  double pressure_cd_tolerance = .015, pitching_cm_tolerance = .015;
  double target = 0;
};

double interpolate_experimental_component(const NacaExperimentalLiftPolar &polar,
    double angle_degrees, int component)
{
  if (!polar.points || polar.point_count == 0
      || angle_degrees < polar.points[0].angle_degrees
      || angle_degrees > polar.points[polar.point_count - 1].angle_degrees)
    return std::numeric_limits<double>::quiet_NaN();
  auto value = [component](const NacaExperimentalLiftPoint &point)
  {
    if (component == 1) return point.pressure_drag_coefficient;
    if (component == 2) return point.moment_coefficient;
    return point.cl;
  };
  for (std::size_t q = 0; q < polar.point_count; ++q)
  {
    if (angle_degrees == polar.points[q].angle_degrees) return value(polar.points[q]);
    if (q + 1 < polar.point_count && angle_degrees < polar.points[q + 1].angle_degrees)
    {
      const auto &a = polar.points[q], &b = polar.points[q + 1];
      const double t = (angle_degrees - a.angle_degrees)
          / (b.angle_degrees - a.angle_degrees);
      return value(a) + t * (value(b) - value(a));
    }
  }
  return value(polar.points[polar.point_count - 1]);
}

void diagnose_connectivity(const CompositeAmrPressureSystem &system, const AmrEmbeddedBoundaryAtlas &eb, const TriangleBvh &bvh, const TriMesh &wing)
{
  const Aabb3d geometry_box{ { wing.bbox_min[0], wing.bbox_min[1], wing.bbox_min[2] }, { wing.bbox_max[0], wing.bbox_max[1], wing.bbox_max[2] } };
  Vec3d point = (geometry_box.lo + geometry_box.hi) * .5;
  // The centre of a rotated, cambered bounding box is not guaranteed to lie inside
  // the closed validation airfoil. Find the two skin intersections of its mid-chord
  // vertical line and seed the pressure graph halfway between them. This ray parity
  // is diagnostic-only; production EB construction remains entirely two-sided.
  const double ray_padding = system.hierarchy->finest_cell_size();
  const Vec3d ray_a{ point.x, point.y, geometry_box.lo.z - ray_padding };
  const Vec3d ray_b{ point.x, point.y, geometry_box.hi.z + ray_padding };
  std::vector<double> ray_hits;
  double t_min = 0;
  for (int q = 0; q < 32; ++q)
  {
	const SegmentHit hit = bvh.intersect_segment(ray_a, ray_b, t_min, 1.0);
	if (!hit.hit) break;
	if (ray_hits.empty() || std::abs(hit.position.z - ray_hits.back()) > 1e-9) ray_hits.push_back(hit.position.z);
	t_min = hit.t + 1e-9;
  }
  if (ray_hits.size() >= 2) point.z = .5 * (ray_hits.front() + ray_hits.back());
  std::printf("  cavity diagnostic seed [%.6g %.6g %.6g] from %zu mid-chord skin hits\n", point.x, point.y, point.z, ray_hits.size());
  const int start = system.pressure_dof_at_point(eb, point);
  if (start < 0 || start >= system.storage_size || !system.active[start])
  {
	std::printf("  cavity probe [%.6g %.6g %.6g] has no active pressure DOF\n", point.x, point.y, point.z);
	return;
  }
  std::printf("  cavity pressure DOF %d initial component=%s\n", start, start < static_cast<int>(system.freestream_connected.size()) && system.freestream_connected[start] ? "freestream" : "quiescent");
  struct Edge
  {
	int neighbour = -1;
	char kind = '?';
	int index = -1;
  };
  std::vector<std::vector<Edge>> graph(system.storage_size);
  auto add = [&](int a, int b, char kind, int index)
  {
	if (a < 0 || b < 0 || a == b || !system.active[a] || !system.active[b]) return;
	graph[a].push_back({ b, kind, index });
	graph[b].push_back({ a, kind, index });
  };
  const int bs = system.brick_size;
  for (int level = 0; level < static_cast<int>(system.hierarchy->levels().size()); ++level)
  {
	const AmrLevel &source = system.hierarchy->levels()[level];
	for (int brick = 0; brick < static_cast<int>(source.bricks.size()); ++brick)
	{
	  const BrickMetadata &meta = source.bricks[brick];
	  if (!meta.active()) continue;
	  for (int k = 0; k < bs; ++k)
		for (int j = 0; j < bs; ++j)
		  for (int i = 0; i < bs; ++i)
		  {
			const int a = system.dof(level, brick, i, j, k);
			if (!system.active[a]) continue;
			for (int axis = 0; axis < 3; ++axis)
			{
			  if (system.cut_face_mask[a] & (1u << axis)) continue;
			  int c[3] = { i, j, k }, b = -1;
			  if (++c[axis] < bs) b = system.dof(level, brick, c[0], c[1], c[2]);
			  else
			  {
				const int neighbour = meta.same_level_neighbor[2 * axis + 1];
				if (neighbour >= 0 && source.bricks[neighbour].active())
				{
				  c[axis] = 0;
				  b = system.dof(level, neighbour, c[0], c[1], c[2]);
				}
			  }
			  add(a, b, 'R', axis);
			}
		  }
	}
  }
  for (int q = 0; q < static_cast<int>(system.coarse_fine.size()); ++q) add(system.coarse_fine[q].coarse_dof, system.coarse_fine[q].fine_dof, 'C', q);
  for (int q = 0; q < static_cast<int>(system.embedded.size()); ++q) add(system.embedded[q].coarse_dof, system.embedded[q].fine_dof, 'E', q);

  std::vector<int> component(system.storage_size, -1);
  int component_count = 0;
  for (int seed = 0; seed < system.storage_size; ++seed)
  {
	if (!system.active[seed] || component[seed] >= 0) continue;
	std::deque<int> work{ seed };
	component[seed] = component_count;
	std::size_t dofs = 0;
	double volume = 0;
	Vec3d lo{ std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
	Vec3d hi{ -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
	bool external = false;
	while (!work.empty())
	{
	  const int a = work.front();
	  work.pop_front();
	  ++dofs;
	  volume += system.volume[a];
	  lo.x = std::min(lo.x, system.centroid[a].x); lo.y = std::min(lo.y, system.centroid[a].y); lo.z = std::min(lo.z, system.centroid[a].z);
	  hi.x = std::max(hi.x, system.centroid[a].x); hi.y = std::max(hi.y, system.centroid[a].y); hi.z = std::max(hi.z, system.centroid[a].z);
	  external = external || (a < static_cast<int>(system.freestream_connected.size()) && system.freestream_connected[a]);
	  for (const Edge &edge : graph[a]) if (component[edge.neighbour] < 0) { component[edge.neighbour] = component_count; work.push_back(edge.neighbour); }
	}
	if (!external)
	{
	  bool gauged = false;
	  for (const CompositePressureGauge &gauge : system.gauges) if (gauge.dof >= 0 && component[gauge.dof] == component_count) gauged = true;
	  std::printf("  isolated component %d: dofs=%zu volume=%.9g bbox=[%.6g %.6g %.6g]-[%.6g %.6g %.6g] gauge=%d seed=%d\n", component_count, dofs, volume, lo.x, lo.y, lo.z, hi.x, hi.y, hi.z, gauged ? 1 : 0, seed);
	}
	++component_count;
  }
  std::printf("  pressure graph has %d connected component(s)\n", component_count);

  std::vector<int> predecessor(system.storage_size, -1), predecessor_edge(system.storage_size, -1);
  std::vector<char> predecessor_kind(system.storage_size, '?');
  std::deque<int> queue;
  predecessor[start] = start;
  queue.push_back(start);
  int outlet = -1;
  while (!queue.empty() && outlet < 0)
  {
	const int a = queue.front();
	queue.pop_front();
	if (system.centroid[a].x > system.hierarchy->domain().hi.x - 1.01 * system.hierarchy->levels().front().h)
	{
	  outlet = a;
	  break;
	}
	for (const Edge &edge : graph[a])
	  if (predecessor[edge.neighbour] < 0)
	  {
		predecessor[edge.neighbour] = a;
		predecessor_edge[edge.neighbour] = edge.index;
		predecessor_kind[edge.neighbour] = edge.kind;
		queue.push_back(edge.neighbour);
	  }
  }
  if (outlet < 0)
  {
	std::printf("  cavity probe DOF %d is isolated from X-max\n", start);
	return;
  }
  std::vector<int> path;
  for (int q = outlet;; q = predecessor[q])
  {
	path.push_back(q);
	if (q == start) break;
  }
  std::reverse(path.begin(), path.end());
  std::printf("  cavity probe DOF %d REACHES X-max through %zu graph edges\n", start, path.size() - 1);
  Aabb3d near = geometry_box;
  const double near_padding = 2 * system.hierarchy->finest_cell_size();
  near.lo = near.lo - Vec3d{near_padding, near_padding, near_padding};
  near.hi = near.hi + Vec3d{near_padding, near_padding, near_padding};
  auto contains = [&](Vec3d p) { return p.x >= near.lo.x && p.x <= near.hi.x && p.y >= near.lo.y && p.y <= near.hi.y && p.z >= near.lo.z && p.z <= near.hi.z; };
  int shown = 0, fabric_crossings = 0;
  const double topology_tolerance = 1e-9 * system.hierarchy->finest_cell_size();
  for (std::size_t q = 1; q < path.size(); ++q)
  {
	const Vec3d a = system.centroid[path[q - 1]], b = system.centroid[path[q]];
	const bool touches_fabric = bvh.segment_touches_surface(a, b, topology_tolerance);
	if (touches_fabric) ++fabric_crossings;
	if (!touches_fabric && (!contains(a) || !contains(b))) continue;
	if (!touches_fabric) continue;
	std::printf("    %c[%d] %d -> %d [%.6g %.6g %.6g] -> [%.6g %.6g %.6g]", predecessor_kind[path[q]], predecessor_edge[path[q]], path[q - 1], path[q], a.x, a.y, a.z, b.x, b.y, b.z);
	if (predecessor_kind[path[q]] == 'E')
	{
	  const CoarseFinePressureConnection &edge = system.embedded[predecessor_edge[path[q]]];
	  const SegmentHit hit = bvh.intersect_segment(a, b, 1e-9, 1.0 - 1e-9);
	  std::printf(" face=[%.6g %.6g %.6g] axis=%d A=%.6g direct-hit=%d", edge.face_centroid.x, edge.face_centroid.y, edge.face_centroid.z, edge.axis, edge.open_area, hit.hit ? 1 : 0);
	  for (const AmrEbLevelAtlas &atlas : eb.levels)
	  {
		const EmbeddedBoundary &topology = atlas.topology;
		for (const FaceAperture &aperture : topology.apertures)
		{
		  if (aperture.axis != edge.axis || std::abs(aperture.area - edge.open_area) > 1e-13 || length2(aperture.centroid - edge.face_centroid) > 1e-20) continue;
		  const auto cell_a = topology.grid.cell_coord(aperture.parent_face_cell);
		  auto cell_b = cell_a;
		  ++cell_b[aperture.axis];
		  const int neighbour_cell = topology.grid.cell_index(cell_b[0], cell_b[1], cell_b[2]);
		  std::printf(" atlas=L%d h=%.6g cells=(%d,%d,%d)/(%d,%d,%d) sample=%u/%u fragments=%u/%u cellFaces=%u/%u", atlas.level, topology.grid.h, cell_a[0], cell_a[1], cell_a[2], cell_b[0], cell_b[1], cell_b[2], topology.cells[aperture.parent_face_cell].sampled_resolution, topology.cells[neighbour_cell].sampled_resolution, topology.cells[aperture.parent_face_cell].fragment_count, topology.cells[neighbour_cell].fragment_count, topology.cells[aperture.parent_face_cell].source_face_id, topology.cells[neighbour_cell].source_face_id);
		  auto print_fragment = [&](FragmentRef ref)
		  {
			if (fragment_is_regular(ref))
			{
			  std::printf(" regular(%d)", regular_fragment_cell(ref));
			  return;
			}
			const int index = irregular_fragment_index(ref);
			const FluidFragment &fragment = topology.fragments[index];
			std::printf(" frag(%d cell=%d sides=%d faces=", index, fragment.parent_cell, fragment.surface_side_count);
			for (int side = 0; side < fragment.surface_side_count; ++side)
			{
			  const FragmentSurfaceSide &surface = topology.fragment_surface_sides[fragment.surface_side_offset + side];
			  std::printf("%s%u:%u", side ? "," : "", surface.source_face_id, surface.side_mask);
			}
			std::printf(")");
		  };
		  std::printf(" refs=");
		  print_fragment(aperture.fragment_a);
		  print_fragment(aperture.fragment_b);
		}
	  }
	}
	std::printf("\n");
	if (++shown == 60)
	{
	  std::printf("    ... path output truncated\n");
	  break;
	}
  }
  std::printf("  graph path has %d centroid-leg contact(s) with fabric\n", fabric_crossings);
}

LiftReference lift_reference(const NacaFourDigitDefinition &definition, double angle_degrees,
    bool smooth_wall, const NacaExperimentalLiftPolar *experiment = nullptr,
    const NacaExperimentalLiftSummary *summary = nullptr)
{
  LiftReference reference;
  reference.thin = naca_thin_airfoil_lift_coefficient(definition, angle_degrees * radians_per_degree);
  reference.experimental = experiment
      ? experimental_lift_coefficient(*experiment, angle_degrees)
      : (summary ? summary->lift_curve_slope_per_degree
          * (angle_degrees - summary->zero_lift_angle_degrees)
          : std::numeric_limits<double>::quiet_NaN());
  // The NREL/OSU archive reports surface-pressure-integrated drag and quarter-
  // chord moment, matching ParaCFD's current pressure-only loads. NASA TMR's CD
  // is total drag, so it is deliberately not compared until skin friction is
  // a validated production quantity.
  if (experiment == &nrel_naca4415_clean_re1m())
  {
    reference.experimental_pressure_cd =
        interpolate_experimental_component(*experiment, angle_degrees, 1);
    reference.experimental_pitching_cm =
        interpolate_experimental_component(*experiment, angle_degrees, 2);
  }
  const NacaXfoilLiftPolar *polar = naca_xfoil_lift_polar(definition.code);
  if (!polar)
  {
	reference.target = reference.experimental;
	return reference;
  }
  auto interpolate = [&](int column)
  {
	const auto value = [&](int q)
	{
	  if (column == 1) return polar->points[q].viscous_cl;
	  if (column == 2) return polar->points[q].forced_turbulent_cl;
	  return polar->points[q].inviscid_cl;
	};
	int a = 0, b = 1;
	if (angle_degrees > polar->points[1].angle_degrees)
	{
	  a = 1;
	  b = 2;
	}
	const double va = value(a), vb = value(b);
	if (!std::isfinite(va) || !std::isfinite(vb)) return std::numeric_limits<double>::quiet_NaN();
	const double t = (angle_degrees - polar->points[a].angle_degrees) / (polar->points[b].angle_degrees - polar->points[a].angle_degrees);
	return va + t * (vb - va);
  };
  reference.xfoil_inviscid = interpolate(0);
  reference.xfoil_free_transition = interpolate(1);
  reference.xfoil_forced_turbulent = interpolate(2);
  reference.target = std::isfinite(reference.experimental) ? reference.experimental
      : (smooth_wall ? reference.xfoil_forced_turbulent : reference.xfoil_inviscid);
  return reference;
}

ModelPlacement pitch_placement(ModelPlacement source, double degrees)
{
  const double radians = degrees * radians_per_degree, c = std::cos(radians), s = std::sin(radians), rotation[9] = { c, 0, s, 0, 1, 0, -s, 0, c };
  ModelPlacement out = source;
  for (int row = 0; row < 3; ++row)
	for (int column = 0; column < 3; ++column) out.m[3 * row + column] = rotation[3 * row] * source.m[column] + rotation[3 * row + 1] * source.m[3 + column] + rotation[3 * row + 2] * source.m[6 + column];
  out.tx = out.ty = out.tz = 0;
  return out;
}

double mesh_area(const TriMesh &mesh)
{
  double area = 0;
  for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
  {
	const std::uint32_t ia = mesh.indices[3 * triangle], ib = mesh.indices[3 * triangle + 1], ic = mesh.indices[3 * triangle + 2];
	const Vec3d a{ mesh.positions[3 * ia], mesh.positions[3 * ia + 1], mesh.positions[3 * ia + 2] }, b{ mesh.positions[3 * ib], mesh.positions[3 * ib + 1], mesh.positions[3 * ib + 2] }, c{ mesh.positions[3 * ic], mesh.positions[3 * ic + 1], mesh.positions[3 * ic + 2] };
	area += 0.5 * std::sqrt(length2(cross(b - a, c - a)));
  }
  return area;
}

Vec3d nearest_cell_velocity(const AmrHierarchy &hierarchy, const AmrHostFields &fields, Vec3d point)
{
  const BrickLocation location = hierarchy.locate_finest(point);
  if (!location.found()) return {};
  const AmrHostLevelFields &level = fields.levels()[location.level];
  const int b = location.brick, i = location.cell.x, j = location.cell.y, k = location.cell.z;
  return { 0.5 * (static_cast<double>(level.u[level.layout.u_index(b, i, j, k)]) + static_cast<double>(level.u[level.layout.u_index(b, i + 1, j, k)])), 0.5 * (static_cast<double>(level.v[level.layout.v_index(b, i, j, k)]) + static_cast<double>(level.v[level.layout.v_index(b, i, j + 1, k)])), 0.5 * (static_cast<double>(level.w[level.layout.w_index(b, i, j, k)]) + static_cast<double>(level.w[level.layout.w_index(b, i, j, k + 1)])) };
}

double circulation_cl(const ExternalAeroCore &core, const TriMesh &wing, double speed, double chord)
{
  AmrHostFields fields(core.hierarchy());
  core.download_fields(fields);
  const Aabb3d domain = core.hierarchy().domain();
  const double h = core.hierarchy().finest_cell_size(), pad = .35 * chord, x0 = std::max(domain.lo.x + h, wing.bbox_min[0] - pad), x1 = std::min(domain.hi.x - h, wing.bbox_max[0] + pad), z0 = std::max(domain.lo.z + h, wing.bbox_min[2] - pad), z1 = std::min(domain.hi.z - h, wing.bbox_max[2] + pad), y = .5 * (wing.bbox_min[1] + wing.bbox_max[1]);
  constexpr int samples = 512;
  double circulation = 0;
  for (int q = 0; q < samples; ++q)
  {
	const double t = (q + .5) / samples, x = x0 + t * (x1 - x0), z = z0 + t * (z1 - z0);
	circulation += nearest_cell_velocity(core.hierarchy(), fields, { x, y, z0 }).x * (x1 - x0) / samples;
	circulation += nearest_cell_velocity(core.hierarchy(), fields, { x1, y, z }).z * (z1 - z0) / samples;
	circulation -= nearest_cell_velocity(core.hierarchy(), fields, { x, y, z1 }).x * (x1 - x0) / samples;
	circulation -= nearest_cell_velocity(core.hierarchy(), fields, { x0, y, z }).z * (z1 - z0) / samples;
  }
  return -2 * circulation / (speed * chord);
}

CaseResult run_case(const std::string &config_path, const std::string &step_path,
    double angle, double target_time, double average_window, int max_levels,
    int subdivisions, double phase_fraction, bool smooth_wall,
    bool conservative_momentum, double benchmark_reynolds)
{
  ParagliderConfig config;
  std::string error;
  if (!load_paraglider_config(config_path, config, &error)) throw std::runtime_error(error);
  config.step_path = step_path;
  if (std::isfinite(validation_tessellation_deflection_mm)) config.tessellation_deflection_mm = validation_tessellation_deflection_mm;
  config.amr.complex_subdivisions = subdivisions;
  if (validation_topology_refinement_levels >= 0) config.amr.topology_refinement_levels = validation_topology_refinement_levels;
  config.solver.projection_max_iterations = std::max(config.solver.projection_max_iterations, 3000);
  config.solver.smagorinsky_cs = validation_smagorinsky_cs;
  if (std::isfinite(validation_kinematic_viscosity))
    config.freestream.nu = validation_kinematic_viscosity;
  else if (std::isfinite(benchmark_reynolds) && benchmark_reynolds > 0)
    config.freestream.nu = config.freestream.speed / benchmark_reynolds;
  if (std::isfinite(validation_projection_tolerance)) config.solver.projection_tolerance = validation_projection_tolerance;
  if (std::isfinite(validation_upstream_margin)) config.domain.upstream_margin = validation_upstream_margin;
  if (std::isfinite(validation_downstream_margin)) config.domain.downstream_margin = validation_downstream_margin;
  if (std::isfinite(validation_vertical_margin)) config.domain.vertical_margin = validation_vertical_margin;
  config.placement = pitch_placement(config.placement, angle);
  TriMesh source = load_step_mesh(step_path, config.tessellation_deflection_mm, &error);
  if (source.empty()) throw std::runtime_error(error);
  constexpr double requested_width = .125;
  double h = 0, width = 0;
  int layers = 0;
  if (validation_local_amr)
  {
	config.amr.base_cell_size = .03125;
	config.amr.max_levels = max_levels;
	config.amr.brick_size = 4;
	config.amr.wing_refinement_distance = .05;
	config.amr.surface_refinement_distance = validation_surface_distance;
	config.amr.wake_length = 1.5;
	config.amr.wake_radius = .15;
	if (std::isfinite(validation_min_volume_fraction)) config.amr.min_volume_fraction = validation_min_volume_fraction;
	h = config.amr.base_cell_size / (1 << std::max(0, max_levels - 1));
	width = requested_width;
	layers = 4;
  }
  else
  {
	const int ratio = 1 << std::max(0, max_levels - 1);
	h = config.amr.base_cell_size / ratio;
	layers = validation_span_layers > 0 ? validation_span_layers : std::max(4, static_cast<int>(std::ceil(requested_width / h - 1e-9)));
	width = layers * h;
	config.amr.base_cell_size = h;
	config.amr.max_levels = 1;
	config.amr.brick_size = layers;
  }
  ModelPlacement orientation = config.placement;
  orientation.tx = orientation.ty = orientation.tz = 0;
  const TriMesh oriented = placed_mesh(source, orientation);
  const double centre = .5 * (oriented.bbox_min[1] + oriented.bbox_max[1]);
  TriMesh clipped = clip_mesh_to_axis_slab(oriented, 1, centre - .5 * width, centre + .5 * width);
  if (clipped.empty()) throw std::runtime_error("NACA validation span crop is empty");
  config.domain.lateral_margin = 0;
  config.domain.upstream_margin += phase_fraction * h;
  config.reference.area = width;
  config.reference.length = 1;
  config.placement = frame_wing_for_external_domain(clipped, ModelPlacement{}, config.domain.upstream_margin, 0, config.domain.vertical_margin, config.amr.base_cell_size * config.amr.brick_size);
  TriMesh wing = placed_mesh(clipped, config.placement);
  TriangleBvh bvh(wing);
  const auto start = std::chrono::steady_clock::now();
  const auto elapsed_seconds = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); };
  const double reynolds = config.freestream.nu > 0
      ? config.freestream.speed * config.reference.length / config.freestream.nu
      : std::numeric_limits<double>::infinity();
  std::printf("  building %s alpha=%+.4f deg (%s, Re=%.7g, nu=%.7g m2/s, flow-collar h=%.7g m, width=%.5g m, grid-phase=%.3fh)\n",
      std::filesystem::path(step_path).filename().string().c_str(), angle,
      smooth_wall ? "turbulent smooth wall" : "slip wall", reynolds,
      config.freestream.nu, h, width, phase_fraction);
  ExternalAeroExecutionOptions execution;
  execution.exact_cell_decomposer = &decompose_exact_cell;
  execution.use_exact_cell_decomposer_as_development_oracle = true;
  execution.smooth_fabric_wall = smooth_wall;
  execution.conservative_cell_momentum = conservative_momentum;
  ExternalAeroCore core(wing, bvh, config, execution);
  const auto preprocessing_end = std::chrono::steady_clock::now();
  CaseResult result;
  result.preprocess_seconds = std::chrono::duration<double>(preprocessing_end - start).count();
  result.gpu_mib = core.gpu_bytes() / (1024.0 * 1024.0);
  result.finest_h = core.hierarchy().finest_cell_size();
  result.pressure_dofs = core.pressure_system().storage_size;
  result.active_bricks = core.hierarchy().active_brick_count();
  result.gauges = core.pressure_system().gauges.size();
  std::printf("  preprocessing complete in %.2f s: finest h=%.7g m, %zu bricks, %d pressure DOFs, %.1f MiB GPU estimate\n", result.preprocess_seconds, result.finest_h, result.active_bricks, result.pressure_dofs, result.gpu_mib);
  double represented = 0;
  for (const auto &patch : core.pressure_system().surface_patches) represented += patch.area;
  result.coverage = represented / mesh_area(wing);
  for (const auto &level : core.embedded_boundary().levels)
  {
	result.recovered_sides += level.topology.recovered_subgrid_surface_sides;
	result.max_recovery_cells = std::max(result.max_recovery_cells, level.topology.maximum_surface_side_recovery_distance / level.topology.grid.h);
  }
  if (validation_diagnose_connectivity) diagnose_connectivity(core.pressure_system(), core.embedded_boundary(), bvh, wing);
  if (validation_preprocess_only)
  {
	result.converged = true;
	result.wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
	return result;
  }
  ExternalAeroStepStats stats = core.initialize();
  if (!stats.pressure.converged)
  {
	result.wall_seconds = elapsed_seconds();
	std::printf("  initialization FAILED after %.2f s: PCG %d iterations, residual %.3g\n", result.wall_seconds, stats.pressure.iterations, stats.pressure.relative_residual);
	return result;
  }
  std::printf("  initialization converged: PCG %d iterations, residual %.3g\n", stats.pressure.iterations, stats.pressure.relative_residual);
  std::deque<TimedCoefficients> coefficient_samples;
  auto record_load = [&]()
  {
	const AerodynamicLoads loads = core.aerodynamic_loads();
	const double cl = loads.viscous_loads_valid ? loads.cl : loads.cl_pressure;
	if (!loads.force_coefficients_valid || !std::isfinite(cl)) throw std::runtime_error("NACA validation load coefficient is unavailable");
	const double dynamic_pressure = .5 * config.freestream.rho
	    * config.freestream.speed * config.freestream.speed;
	const double pitching_cm = loads.pressure_moment.y
	    / (dynamic_pressure * config.reference.area * config.reference.length);
	coefficient_samples.push_back({core.physical_time(), cl, loads.cl_pressure,
	    loads.viscous_loads_valid ? loads.cl_viscous : 0.0,
	    loads.cd_pressure, pitching_cm});
	result.final_cl = cl;
	// No fixed sample-count averaging: adaptive CFL makes samples nonuniform in
	// physical time. Retain only the two adjacent reporting windows plus a small
	// interpolation margin.
	while (coefficient_samples.size() > 2 && coefficient_samples[1].time < core.physical_time() - 2.1 * average_window)
	  coefficient_samples.pop_front();
  };
  record_load();
  while (core.physical_time() < target_time)
  {
	stats = core.step();
	++result.steps;
	result.physical_time = core.physical_time();
	if (!stats.pressure.converged)
	{
	  result.wall_seconds = elapsed_seconds();
	  std::printf("  step %d FAILED at t=%.6g s: PCG %d iterations, residual %.3g\n", result.steps, result.physical_time, stats.pressure.iterations, stats.pressure.relative_residual);
	  return result;
	}
	const bool final = core.physical_time() >= target_time;
	if (result.steps % 20 == 0 || final) record_load();
	if (validation_progress_steps > 0 && result.steps % validation_progress_steps == 0)
	{
	  const AerodynamicLoads preview = core.aerodynamic_loads();
	  const double cl = preview.viscous_loads_valid ? preview.cl : preview.cl_pressure;
	  std::printf("  progress step=%d t=%.6g/%.6g s dt=%.3g CL=%+.5f vmax=%.4g regular/EB=%.4g/%.4g PCG=%d/%.3g GPU=%.2f ms wall=%.1f s\n", result.steps, result.physical_time, target_time, stats.dt, cl, stats.max_abs_velocity, stats.max_abs_regular_velocity, stats.max_abs_special_velocity, stats.pressure.iterations, stats.pressure.relative_residual, stats.gpu_step_ms, elapsed_seconds());
	}
	const bool step_limit = validation_max_steps > 0 && result.steps >= validation_max_steps;
	const bool wall_limit = validation_max_wall_seconds > 0 && elapsed_seconds() >= validation_max_wall_seconds;
	if (!final && (step_limit || wall_limit))
	{
	  result.timed_out = true;
	  std::printf("  controlled stop at step=%d t=%.6g s after %.1f wall s (%s limit)\n", result.steps, result.physical_time, elapsed_seconds(), step_limit ? "step" : "wall-time");
	  break;
	}
  }
  const double end_time = core.physical_time();
  const CoefficientWindow current = coefficient_window(coefficient_samples, end_time - average_window, end_time);
  const CoefficientWindow previous = coefficient_window(coefficient_samples, end_time - 2 * average_window, end_time - average_window);
  result.converged = !result.timed_out && end_time >= target_time && current.ready && previous.ready;
  if (current.ready)
  {
	result.mean_cl = current.mean;
	result.mean_pressure_cl = current.pressure_mean;
	result.mean_viscous_cl = current.viscous_mean;
	result.mean_pressure_cd = current.pressure_cd_mean;
	result.mean_pitching_cm = current.pitching_cm_mean;
	result.rms_cl = current.rms;
	result.mean_window_coverage = current.coverage;
	if (previous.ready)
	{
	  result.previous_mean_cl = previous.mean;
	  result.mean_relative_drift = std::abs(current.mean - previous.mean) /
		std::max(.05, std::abs(.5 * (current.mean + previous.mean)));
	}
  }
  else
  {
	result.mean_cl = result.mean_pressure_cl = result.mean_viscous_cl =
	    result.mean_pressure_cd = result.mean_pitching_cm = result.rms_cl =
	    std::numeric_limits<double>::quiet_NaN();
  }
  result.circulation_cl = circulation_cl(core, wing, config.freestream.speed, config.reference.length);
  const ExternalAeroConservationStats conservation = core.conservation_stats();
  result.divergence = conservation.max_abs_divergence;
  result.normalized_divergence = result.divergence * result.finest_h / config.freestream.speed;
  result.rms_divergence = conservation.volume_weighted_rms_divergence;
  result.normalized_rms_divergence = result.rms_divergence * result.finest_h / config.freestream.speed;
  result.max_integrated_flux_error = conservation.max_integrated_flux_error;
  result.normalized_integrated_flux_error = conservation.max_integrated_flux_error /
    (config.freestream.speed * result.finest_h * result.finest_h);
  result.net_integrated_flux_error = conservation.net_integrated_flux_error;
  result.residual = stats.pressure.relative_residual;
  result.iterations = stats.pressure.iterations;
  result.final_step_ms = stats.gpu_step_ms;
  result.final_projection_ms = stats.projection_ms;
  std::vector<double> pressure;
  core.download_pressure(pressure);
  double sealed_min = std::numeric_limits<double>::infinity(), sealed_max = -std::numeric_limits<double>::infinity(), sealed_sum = 0;
  for (int q = 0; q < core.pressure_system().storage_size; ++q)
	if (core.pressure_system().active[q] && !core.pressure_system().freestream_connected[q])
	{
	  ++result.sealed_dofs;
	  sealed_min = std::min(sealed_min, pressure[q]);
	  sealed_max = std::max(sealed_max, pressure[q]);
	  sealed_sum += pressure[q];
	}
  if (result.sealed_dofs)
  {
	result.sealed_pressure_span = sealed_max - sealed_min;
	const double sealed_mean = sealed_sum / result.sealed_dofs;
	std::vector<double> quiet_pressure = pressure;
	for (int q = 0; q < core.pressure_system().storage_size; ++q)
	  if (core.pressure_system().active[q] && !core.pressure_system().freestream_connected[q]) quiet_pressure[q] = sealed_mean;
	const AerodynamicLoads actual = compute_pressure_loads(core.pressure_system(), pressure, wing.triangle_count(), config.freestream, config.reference);
	const AerodynamicLoads quiet = compute_pressure_loads(core.pressure_system(), quiet_pressure, wing.triangle_count(), config.freestream, config.reference);
	result.sealed_pressure_cl_delta = actual.cl_pressure - quiet.cl_pressure;
  }
  result.wall_seconds = elapsed_seconds();
  return result;
}

bool mechanics_ok(const CaseResult &value)
{
  // The pointwise maximum is retained as a sliver-cell diagnostic, while the
  // Volume-weighted divergence is the global projection gate.  Retain a strict
  // local diagnostic too, but normalize its actual flux imbalance by U*h^2.
  // The old absolute 1e-8 m^3/s cutoff changed meaning with resolution and was
  // therefore not a valid convergence criterion.  2.5e-4 permits at most 0.025%
  // of one finest-cell freestream face flux in any pressure control volume.
  return value.converged && !value.timed_out && value.gauges == 1 &&
    std::abs(value.coverage - 1) < 1e-7 && value.normalized_divergence < 2.5e-4 &&
    value.normalized_rms_divergence < 1e-5 && value.normalized_integrated_flux_error < 2.5e-4 &&
    value.residual <= 1.1e-5 && value.rms_cl < 0.025 && value.mean_relative_drift < .01 &&
    std::abs(value.final_cl - value.circulation_cl) < 0.07;
}

double coefficient_tolerance(double target)
{
  // The primary experimental pre-stall lift gate. This combines the stated
  // wind-tunnel uncertainty with a declared coarse Cartesian discretization
  // allowance; resolution/domain convergence remain separate mandatory gates.
  return std::max(0.02, 0.05 * std::abs(target));
}

bool reference_coefficients_ok(const LiftReference &reference,
    const CaseResult &value)
{
  if (!std::isfinite(reference.target)
      || std::abs(value.mean_cl - reference.target) > coefficient_tolerance(reference.target))
    return false;
  if (std::isfinite(reference.experimental_pressure_cd)
      && std::abs(value.mean_pressure_cd - reference.experimental_pressure_cd)
          > reference.pressure_cd_tolerance)
    return false;
  if (std::isfinite(reference.experimental_pitching_cm)
      && std::abs(value.mean_pitching_cm - reference.experimental_pitching_cm)
          > reference.pitching_cm_tolerance)
    return false;
  return true;
}

void print_case(const char *code, double alpha, const LiftReference &reference, const CaseResult &value)
{
  std::printf("NACA%s alpha=%+8.4f  CL mean[P+V]/prev/final/circ=%+9.5f[%+.5f%+.5f]/%+9.5f/%+9.5f/%+9.5f  drift/rms=%.3g/%.3g  EXP=%+8.4f XFOIL I/F/T=%+8.4f/%+8.4f/%+8.4f target=%+8.4f err=%+.5f tol=%.5f thin=%+8.4f  div[max/rms]=%.3g/%.3g (%.3g/%.3g U/h) flux[max/net]=%.3g/%.3g (%.3g U*h^2) PCG=%d/%.3g gauges=%zu coverage=%.9f recovered=%zu/max%.3fh sealed=%zu dP=%.3g dCL=%.4g bricks/DOFs=%zu/%d GPU=%.1fMiB step/proj=%.2f/%.2fms prep=%.2fs steps=%d t=%.4fs wall=%.2fs  %s%s\n", code, alpha, value.mean_cl, value.mean_pressure_cl, value.mean_viscous_cl, value.previous_mean_cl, value.final_cl, value.circulation_cl, value.mean_relative_drift, value.rms_cl, reference.experimental, reference.xfoil_inviscid, reference.xfoil_free_transition, reference.xfoil_forced_turbulent, reference.target, value.mean_cl - reference.target, coefficient_tolerance(reference.target), reference.thin, value.divergence, value.rms_divergence, value.normalized_divergence, value.normalized_rms_divergence, value.max_integrated_flux_error, value.net_integrated_flux_error, value.normalized_integrated_flux_error, value.iterations, value.residual, value.gauges, value.coverage, value.recovered_sides, value.max_recovery_cells, value.sealed_dofs, value.sealed_pressure_span, value.sealed_pressure_cl_delta, value.active_bricks, value.pressure_dofs, value.gpu_mib, value.final_step_ms, value.final_projection_ms, value.preprocess_seconds, value.steps, value.physical_time, value.wall_seconds, mechanics_ok(value) ? "FLOW PASS" : "FLOW FAIL", value.timed_out ? " [CONTROLLED TIMEOUT]" : "");
  if (std::isfinite(reference.experimental_pressure_cd)
      || std::isfinite(reference.experimental_pitching_cm))
    std::printf("         CDp=%+.6f (EXP=%+.6f, tol=%.4f), Cm,c/4=%+.6f (EXP=%+.6f, tol=%.4f)  %s\n",
        value.mean_pressure_cd, reference.experimental_pressure_cd,
        reference.pressure_cd_tolerance, value.mean_pitching_cm,
        reference.experimental_pitching_cm, reference.pitching_cm_tolerance,
        reference_coefficients_ok(reference, value) ? "AERO PASS" : "AERO FAIL");
}
} // namespace

int main(int argc, char **argv)
{
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  std::string config_path = "configs/naca_validation.json", profile_filter;
  // A half-cell origin phase is the deterministic production validation phase:
  // analytical leading/trailing-edge vertices then cannot lie exactly on a
  // Cartesian pressure-graph plane.  --phase-fraction remains available for
  // the explicit phase-sensitivity diagnostic.
  double target_time = .8, average_window = .2, single_angle = std::numeric_limits<double>::quiet_NaN(), single_phase = .5;
  int max_levels = 4, subdivisions = 16;
  bool phase_check = false, smooth_wall = false, conservative_momentum = false;
  for (int i = 1; i < argc; ++i)
  {
	const std::string argument = argv[i];
	if (argument == "--config" && i + 1 < argc) config_path = argv[++i];
	else if (argument == "--profile" && i + 1 < argc) profile_filter = argv[++i];
	else if (argument == "--single-angle" && i + 1 < argc) single_angle = std::atof(argv[++i]);
	else if (argument == "--phase-fraction" && i + 1 < argc) single_phase = std::atof(argv[++i]);
	else if (argument == "--span-layers" && i + 1 < argc) validation_span_layers = std::atoi(argv[++i]);
	else if (argument == "--topology-refinement-levels" && i + 1 < argc) validation_topology_refinement_levels = std::atoi(argv[++i]);
	else if (argument == "--local-amr") validation_local_amr = true;
	else if (argument == "--uniform") validation_local_amr = false;
	else if (argument == "--physical-time" && i + 1 < argc) target_time = std::atof(argv[++i]);
	else if (argument == "--average-window" && i + 1 < argc) average_window = std::atof(argv[++i]);
	else if (argument == "--max-levels" && i + 1 < argc) max_levels = std::atoi(argv[++i]);
	else if (argument == "--complex-subdivisions" && i + 1 < argc) subdivisions = std::atoi(argv[++i]);
	else if (argument == "--max-wall-seconds" && i + 1 < argc) validation_max_wall_seconds = std::atof(argv[++i]);
	else if (argument == "--max-steps" && i + 1 < argc) validation_max_steps = std::atoi(argv[++i]);
	else if (argument == "--progress-steps" && i + 1 < argc) validation_progress_steps = std::atoi(argv[++i]);
	else if (argument == "--surface-distance" && i + 1 < argc) validation_surface_distance = std::atof(argv[++i]);
	else if (argument == "--min-volume-fraction" && i + 1 < argc) validation_min_volume_fraction = std::atof(argv[++i]);
	else if (argument == "--smagorinsky-cs" && i + 1 < argc) validation_smagorinsky_cs = std::atof(argv[++i]);
	else if (argument == "--kinematic-viscosity" && i + 1 < argc) validation_kinematic_viscosity = std::atof(argv[++i]);
	else if (argument == "--projection-tolerance" && i + 1 < argc) validation_projection_tolerance = std::atof(argv[++i]);
	else if (argument == "--tessellation-mm" && i + 1 < argc) validation_tessellation_deflection_mm = std::atof(argv[++i]);
	else if (argument == "--upstream-margin" && i + 1 < argc) validation_upstream_margin = std::atof(argv[++i]);
	else if (argument == "--downstream-margin" && i + 1 < argc) validation_downstream_margin = std::atof(argv[++i]);
	else if (argument == "--vertical-margin" && i + 1 < argc) validation_vertical_margin = std::atof(argv[++i]);
	else if (argument == "--step-file" && i + 1 < argc) validation_step_override = argv[++i];
	else if (argument == "--experimental-benchmark" && i + 1 < argc) validation_experimental_benchmark = argv[++i];
	else if (argument == "--phase-check") phase_check = true;
	else if (argument == "--no-phase-check") phase_check = false;
	else if (argument == "--preprocess-only") validation_preprocess_only = true;
	else if (argument == "--diagnose-connectivity") validation_diagnose_connectivity = true;
	else if (argument == "--smooth-wall") smooth_wall = true;
	else if (argument == "--no-smooth-wall") smooth_wall = false;
	else if (argument == "--conservative-cell-momentum") conservative_momentum = true;
	else if (argument == "--staggered-momentum") conservative_momentum = false;
	else
	{
	  std::fprintf(stderr, "usage: naca_validation_probe [--config file] [--experimental-benchmark cambered|naca-tr460|nrel-4415|nasa-tmr-0012] [--profile 0012|0021|2412|4412|6412|4415|4112] [--step-file exact.step] [--single-angle deg --phase-fraction cells] [--local-amr|--uniform --span-layers N] [--physical-time s] [--average-window s] [--max-levels N] [--topology-refinement-levels N] [--complex-subdivisions N] [--surface-distance m] [--min-volume-fraction 0..1] [--upstream-margin m] [--downstream-margin m] [--vertical-margin m] [--tessellation-mm mm] [--smagorinsky-cs Cs] [--kinematic-viscosity m2/s] [--projection-tolerance value] [--max-wall-seconds s] [--max-steps N] [--progress-steps N] [--phase-check|--no-phase-check] [--preprocess-only] [--diagnose-connectivity] [--smooth-wall] [--conservative-cell-momentum]\n");
	  return 2;
	}
  }
  ExperimentalBenchmark benchmark = ExperimentalBenchmark::none;
  if (!validation_experimental_benchmark.empty())
  {
	if (validation_experimental_benchmark == "nasa-tmr-0012")
	{
	  benchmark = ExperimentalBenchmark::nasa_tmr_0012;
	  profile_filter = "0012";
	}
	else if (validation_experimental_benchmark == "naca-tr460")
	  benchmark = ExperimentalBenchmark::naca_tr460_cambered;
	else if (validation_experimental_benchmark == "nrel-4415")
	{
	  benchmark = ExperimentalBenchmark::nrel_4415;
	  profile_filter = "4415";
	}
	else if (validation_experimental_benchmark == "cambered")
	  benchmark = ExperimentalBenchmark::cambered_suite;
	else
	{
	  std::fprintf(stderr, "unknown experimental benchmark '%s'\n", validation_experimental_benchmark.c_str());
	  return 2;
	}
	smooth_wall = true;
  }
  if (!(target_time > 0 && average_window > 0 && 2 * average_window < target_time) || max_levels < 1 || validation_topology_refinement_levels < -1 || validation_topology_refinement_levels > 2 || max_levels + std::max(0, validation_topology_refinement_levels) > 10 || subdivisions < 2 || validation_surface_distance < 0 || (std::isfinite(validation_min_volume_fraction) && !(validation_min_volume_fraction > 0 && validation_min_volume_fraction < 1)) || (std::isfinite(validation_upstream_margin) && !(validation_upstream_margin > 0)) || (std::isfinite(validation_downstream_margin) && !(validation_downstream_margin > 0)) || (std::isfinite(validation_vertical_margin) && !(validation_vertical_margin > 0)) || (std::isfinite(validation_smagorinsky_cs) && validation_smagorinsky_cs < 0) || (std::isfinite(validation_kinematic_viscosity) && validation_kinematic_viscosity < 0) || (std::isfinite(validation_projection_tolerance) && !(validation_projection_tolerance > 0)) || (std::isfinite(validation_tessellation_deflection_mm) && !(validation_tessellation_deflection_mm > 0)) || validation_max_wall_seconds < 0 || validation_max_steps < 0 || validation_progress_steps < 0 || (validation_span_layers > 0 && validation_span_layers < 4))
  {
	std::fprintf(stderr, "invalid NACA validation controls\n");
	return 2;
  }
  // Slip-wall comparison uses the Euler equations, matching XFOIL's inviscid
  // column.  Molecular viscosity plus Smagorinsky with a slip wall is neither
  // that model nor XFOIL's viscous boundary layer and previously made the
  // accepted coefficient gate compare unlike PDEs. Explicit CLI overrides stay
  // available for sensitivity studies and the exploratory smooth-wall path.
  if (!std::isfinite(validation_smagorinsky_cs)) validation_smagorinsky_cs = smooth_wall ? .10 : 0.0;
  if (!std::isfinite(validation_kinematic_viscosity) && !smooth_wall) validation_kinematic_viscosity = 0.0;
  const std::array<const char *, 7> profiles = {
      "2412", "4412", "6412", "4415", "0012", "0021", "4112" };
  int failures = 0, matched_profiles = 0;
  double total_wall = 0;
  ParagliderConfig reported_config;
  std::string reported_config_error;
  if (!load_paraglider_config(config_path, reported_config, &reported_config_error))
  {
	std::fprintf(stderr, "could not load validation config for report metadata: %s\n",
	    reported_config_error.c_str());
	return 2;
  }
  const int reported_levels = max_levels + std::max(0, validation_topology_refinement_levels);
  // run_case overrides the local-AMR base spacing to 0.03125 m. Uniform
  // validation instead starts from the selected config's base spacing before
  // collapsing the requested level ratio into one uniform level.
  const double reported_base_h = validation_local_amr
      ? .03125 : reported_config.amr.base_cell_size;
  const double reported_h = reported_base_h / (1 << (reported_levels - 1));
  const std::string projection_tolerance_label = std::isfinite(validation_projection_tolerance) ? std::to_string(validation_projection_tolerance) : "config";
  const std::string tessellation_label = std::isfinite(validation_tessellation_deflection_mm) ? std::to_string(validation_tessellation_deflection_mm) : "config";
  const std::string viscosity_label=std::isfinite(validation_kinematic_viscosity)?std::to_string(validation_kinematic_viscosity):"config";
  const std::string merge_label=std::isfinite(validation_min_volume_fraction)?std::to_string(validation_min_volume_fraction):"config";
  const char *reference_name = "NONE (XFOIL diagnostic only)";
  if (benchmark == ExperimentalBenchmark::nasa_tmr_0012)
    reference_name = nasa_tmr_naca0012_ladson_120_grit().id;
  else if (benchmark == ExperimentalBenchmark::naca_tr460_cambered)
    reference_name = "NACA-TR-460 wind-tunnel tables";
  else if (benchmark == ExperimentalBenchmark::nrel_4415)
    reference_name = nrel_naca4415_clean_re1m().id;
  else if (benchmark == ExperimentalBenchmark::cambered_suite)
    reference_name = "TR-460 + NREL/OSU cambered experimental suite";
  std::printf("ParaCFD NACA validation: solver=%s, %s flow-collar h=%.7g m, collar=%.4g m, topology=%d^3, min-volume=%s h^3, t=%.3g s, mean-window=%.3g s, tessellation=%s mm, nu=%s m2/s, Cs=%.4g, projection-tol=%s, primary-reference=%s\n", conservative_momentum ? "collocated-conservative" : "face-centred-MAC", validation_local_amr ? "local-AMR" : "uniform", reported_h, validation_surface_distance, subdivisions, merge_label.c_str(), target_time, average_window, tessellation_label.c_str(), viscosity_label.c_str(), validation_smagorinsky_cs, projection_tolerance_label.c_str(), reference_name);
  for (const char *code : profiles)
  {
	if (!profile_filter.empty() && profile_filter != code) continue;
	const NacaExperimentalLiftPolar *experimental_polar = nullptr;
	const NacaExperimentalLiftSummary *experimental_summary = nullptr;
	if (benchmark == ExperimentalBenchmark::nasa_tmr_0012)
	{
	  if (std::string(code) != "0012") continue;
	  experimental_polar = &nasa_tmr_naca0012_ladson_120_grit();
	}
	else if (benchmark == ExperimentalBenchmark::nrel_4415)
	{
	  if (std::string(code) != "4415") continue;
	  experimental_polar = &nrel_naca4415_clean_re1m();
	}
	else if (benchmark == ExperimentalBenchmark::naca_tr460_cambered
	    || benchmark == ExperimentalBenchmark::cambered_suite)
	{
	  if (std::string(code) == "4415"
	      && benchmark == ExperimentalBenchmark::cambered_suite)
	    experimental_polar = &nrel_naca4415_clean_re1m();
	  else
	    experimental_summary = naca_tr460_lift_summary(code);
	  if (!experimental_polar && !experimental_summary) continue;
	}
	else if (!naca_xfoil_lift_polar(code)) continue;
	++matched_profiles;
	try
	{
	  NacaFourDigitDefinition definition;
	  std::string error;
	  if (!parse_naca_four_digit(code, definition, &error)) throw std::runtime_error(error);
	  const NacaXfoilLiftPolar *polar = naca_xfoil_lift_polar(code);
	  if (!polar && !experimental_summary && !experimental_polar)
	    throw std::runtime_error("missing validation reference");
	  std::string step;
	  if (!validation_step_override.empty()) step = validation_step_override;
	  else if (experimental_polar == &nasa_tmr_naca0012_ladson_120_grit())
	    step = "Test-Data/NASA-TMR-NACA0012-C1-S2p03.step";
	  else if (experimental_polar == &nrel_naca4415_clean_re1m())
	    step = "Test-Data/NREL-NACA4415-C1-S2p03.step";
	  else if (experimental_summary)
	    step = "Test-Data/NACA" + std::string(code) + "-TR460-C1-S2p03.step";
	  else
	    step = "Test-Data/NACA" + std::string(code) + "_C1_S2p03.step";
	  const double benchmark_reynolds = experimental_polar
	      ? experimental_polar->reynolds
	      : (experimental_summary ? experimental_summary->reynolds
	          : std::numeric_limits<double>::quiet_NaN());
	  if (!std::filesystem::exists(step)) throw std::runtime_error("missing " + step);
	  if (std::isfinite(single_angle))
	  {
		const LiftReference reference = lift_reference(definition, single_angle, smooth_wall, experimental_polar, experimental_summary);
		const CaseResult value = run_case(config_path, step, single_angle,
		    target_time, average_window, max_levels, subdivisions, single_phase,
		    smooth_wall, conservative_momentum, benchmark_reynolds);
		total_wall += value.wall_seconds;
		if (validation_preprocess_only)
		{
		  const bool topology_ok = value.gauges == 1 && std::abs(value.coverage - 1) < 1e-7;
		  std::printf("NACA%s alpha=%+8.4f preprocessing: gauges=%zu coverage=%.9f recovered=%zu/max%.3fh bricks/DOFs=%zu/%d GPU=%.1fMiB prep=%.2fs  %s\n", code, single_angle, value.gauges, value.coverage, value.recovered_sides, value.max_recovery_cells, value.active_bricks, value.pressure_dofs, value.gpu_mib, value.preprocess_seconds, topology_ok ? "TOPOLOGY PASS" : "TOPOLOGY FAIL");
		  if (!topology_ok) ++failures;
		  continue;
		}
		print_case(code, single_angle, reference, value);
		if (!mechanics_ok(value) || !reference_coefficients_ok(reference, value)) ++failures;
		continue;
	  }
	  std::array<double,3> angles{};
	  if (experimental_polar == &nrel_naca4415_clean_re1m())
	    angles = { -4.1, 0.0, 4.1 };
	  else if (experimental_polar == &nasa_tmr_naca0012_ladson_120_grit())
	    angles = { -0.01, 2.15, 4.11 };
	  else if (experimental_summary)
	    angles = { experimental_summary->zero_lift_angle_degrees,
	        experimental_summary->zero_lift_angle_degrees + 2.0,
	        experimental_summary->zero_lift_angle_degrees + 4.0 };
	  else
	    angles = { polar->points[0].angle_degrees, polar->points[1].angle_degrees,
	        polar->points[2].angle_degrees };
	  std::array<CaseResult, 3> cases;
	  for (int q = 0; q < 3; ++q)
	  {
		const double alpha = angles[q];
		const LiftReference reference = lift_reference(definition, alpha, smooth_wall, experimental_polar, experimental_summary);
		// Avoid exact CAD vertices landing on Cartesian pressure-graph planes.  This
		// changes only the grid origin, not the geometry or cell size.
		const double phase = .5;
		cases[q] = run_case(config_path, step, alpha, target_time, average_window,
		    max_levels, subdivisions, phase, smooth_wall, conservative_momentum,
		    benchmark_reynolds);
		total_wall += cases[q].wall_seconds;
		print_case(code, alpha, reference, cases[q]);
		if (!mechanics_ok(cases[q]) || !reference_coefficients_ok(reference, cases[q])) ++failures;
	  }
	  const double reference_low = lift_reference(definition, angles[0], smooth_wall, experimental_polar, experimental_summary).target;
	  const double reference_mid = lift_reference(definition, angles[1], smooth_wall, experimental_polar, experimental_summary).target;
	  const double reference_high = lift_reference(definition, angles[2], smooth_wall, experimental_polar, experimental_summary).target;
	  const double angle_span=angles[2]-angles[0];
	  const double slope = (cases[2].mean_cl - cases[0].mean_cl) / angle_span, reference_slope = (reference_high - reference_low) / angle_span, reference_zero = angles[1] - reference_mid / reference_slope, fit_zero = angles[1] - cases[1].mean_cl / slope, slope_error = std::abs(slope / reference_slope - 1), zero_error = std::abs(fit_zero - reference_zero);
	  // alpha0 is extrapolated outside the sampled attached positive-angle range,
	  // especially for NACA4415; one degree is a more honest secondary gate than
	  // treating that extrapolation as a directly simulated zero-lift point.
	  bool fit_ok = std::isfinite(slope) && slope_error <= .15 && zero_error <= 1.1;
	  std::printf("NACA%s fit: dCL/dalpha=%.6f/deg (reference %.6f, err %.2f%%), alpha0=%.4f deg (reference %.4f, delta %.3f)  %s\n", code, slope, reference_slope, 100 * slope_error, fit_zero, reference_zero, fit_zero - reference_zero, fit_ok ? "FIT PASS" : "FIT FAIL");
	  if (!fit_ok) ++failures;
	  if (phase_check)
	  {
		const double alpha = angles[2];
		const LiftReference reference = lift_reference(definition, alpha, smooth_wall, experimental_polar, experimental_summary);
		const double base_phase = .5;
		const CaseResult shifted = run_case(config_path, step, alpha, target_time,
		    average_window, max_levels, subdivisions,
		    base_phase == 0.0 ? .5 : 0.0, smooth_wall, conservative_momentum,
		    benchmark_reynolds);
		total_wall += shifted.wall_seconds;
		print_case(code, alpha, reference, shifted);
		const double spread = std::abs(shifted.mean_cl - cases[2].mean_cl);
		const bool phase_ok = mechanics_ok(shifted)
		    && reference_coefficients_ok(reference, shifted) && spread <= .035;
		std::printf("NACA%s half-cell phase spread: %.6f CL  %s\n", code, spread, phase_ok ? "PHASE PASS" : "PHASE FAIL");
		if (!phase_ok) ++failures;
	  }
	}
	catch (const std::exception &exception)
	{
	  std::fprintf(stderr, "NACA%s ERROR: %s\n", code, exception.what());
	  ++failures;
	}
  }
  if (matched_profiles == 0)
  {
	std::fprintf(stderr, "unknown NACA validation profile '%s'\n", profile_filter.c_str());
	return 2;
  }
  std::printf("NACA validation: %s (%d failures, %.2f s accumulated case wall time)\n", failures ? "FAIL" : "PASS", failures, total_wall);
  return failures ? 1 : 0;
}
