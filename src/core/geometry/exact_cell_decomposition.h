// exact_cell_decomposition.h -- exact CPU decomposition of one Cartesian cell.
//
// The implementation is OpenCascade-backed, but this public contract deliberately uses
// only standard C++ value types.  It belongs to paracfd_geometry; CFD/GUI targets need not
// inherit any OCCT include directories or types.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace paracfd::core
{
using ExactCellPoint = std::array<double, 3>;

// One planar BRep image with an explicit outer boundary, zero or more holes, and
// a conservative triangulation of the represented region.  Consumers must use
// convex_pieces for geometric work; the loops are retained for diagnostics and
// make hole ownership auditable instead of relying on an untagged wire list.
struct ExactCellPlanarRegion
{
  std::vector<ExactCellPoint> outer_loop;
  std::vector<std::vector<ExactCellPoint>> hole_loops;
  // Zero-area closed wires (occasionally a one-vertex degenerate loop) left by a
  // terminating sheet on the region. They are audit data, not holes and never
  // subtract aperture area.
  std::vector<std::vector<ExactCellPoint>> internal_slit_wires;
  std::size_t empty_internal_wire_count = 0;
  std::vector<std::array<ExactCellPoint, 3>> convex_pieces;
  double convex_piece_area_sum = 0.0;
  ExactCellPoint convex_piece_first_moment{};
};

// Outward-oriented boundary triangle of one closed fluid fragment.  These are
// retained as plain values so the OCC-free core can perform generalized-winding
// membership queries after every OpenCascade object has been discarded.
struct ExactCellOrientedTriangle
{
  std::array<ExactCellPoint, 3> vertices{};
};

struct ExactCellTriangle
{
  static constexpr std::uint64_t no_cad_contact_id =
      std::numeric_limits<std::uint64_t>::max();

  std::array<std::uint32_t, 3> vertices{};
  std::uint64_t triangle_id = 0;
  std::uint64_t source_face_id = 0;
  // Exact BRep edge/face contacts are certified before tessellation.  A contact
  // can therefore be real even when independently tessellated faces miss each
  // other geometrically.  The exact cell decomposer must either reproduce the
  // declared number of half-sheet sectors on one common result edge or reject
  // the cell; these values are audit constraints, never fuzzy-sewing radii.
  std::array<std::uint64_t, 3> cad_contact_ids{
      no_cad_contact_id, no_cad_contact_id, no_cad_contact_id };
  std::array<std::uint16_t, 3> certified_fan_degrees{};
};

struct ExactCellInput
{
  ExactCellPoint cell_min{};
  ExactCellPoint cell_max{};
  // Indexed vertices make conformity explicit: triangles in one sheet must use the
  // same vertex IDs on shared edges.  Independent sheets may intersect without
  // sharing vertices; General Fuse imprints those intersections.
  std::vector<ExactCellPoint> vertices;
  std::vector<ExactCellTriangle> triangles;
  bool run_parallel = false;
};

struct ExactCellFragment
{
  std::int32_t id = -1;
  double volume = 0.0;
  ExactCellPoint centroid{};
  // A BRepClass3d-certified point strictly inside this fragment.  A non-convex
  // control volume's true volume centroid need not lie inside, so membership
  // validation and retained-locator reachability use this witness instead.
  ExactCellPoint interior_witness{};
  std::vector<ExactCellOrientedTriangle> boundary_triangles;
};

struct ExactCellSurfacePatch
{
  std::uint64_t source_triangle_id = 0;
  std::uint64_t source_face_id = 0;
  double area = 0.0;
  ExactCellPoint centroid{};
  ExactCellPoint normal{}; // input winding, minus -> plus
  ExactCellPlanarRegion region;
  std::int32_t plus_fragment = -1;
  std::int32_t minus_fragment = -1;
  // A fabric image coincident with a Cartesian cell face has only one local
  // fluid side.  The neighbouring cell supplies the other side during the
  // conservative shared-face transaction.
  std::int8_t boundary_axis = -1;
  bool boundary_upper = false;
};

struct ExactCellBoxAperture
{
  std::int8_t axis = 0; // 0=x, 1=y, 2=z
  bool upper = false;
  double area = 0.0;
  ExactCellPoint centroid{};
  ExactCellPlanarRegion region;
  // The fluid fragment immediately inside this cell face.  The neighbouring cell
  // supplies the other fragment when two cell-face partitions are intersected.
  std::int32_t fragment = -1;
};

struct ExactCellDecomposition
{
  std::vector<ExactCellFragment> fragments;
  std::vector<ExactCellSurfacePatch> surface_patches;
  std::vector<ExactCellBoxAperture> box_apertures;
  std::vector<std::string> warnings;
  std::vector<std::string> errors;

  double expected_cell_volume = 0.0;
  double fragment_volume_sum = 0.0;
  // Signed sum(fragment volumes) - box volume and the corresponding Euclidean
  // first-moment residual.  These are computed from exact BRep properties.
  double volume_conservation_residual = 0.0;
  double relative_volume_conservation_residual = 0.0;
  double first_moment_conservation_residual = 0.0;
  double relative_first_moment_conservation_residual = 0.0;
  // Independent fabric-history audit. Total cell volume can remain exact even when
  // General Fuse loses a zero-thickness tool image, so this is a separate invariant.
  double input_clipped_surface_area_sum = 0.0;
  double history_surface_area_sum = 0.0;
  // Six Cartesian-face partitions are audited as open aperture + boundary fabric.
  double maximum_box_face_area_residual = 0.0;
  double maximum_box_face_first_moment_residual = 0.0;
  // The API requests zero.  OCCT clamps that request to its kernel-wide
  // Precision::Confusion() floor; exposing both values keeps that limitation
  // auditable instead of describing it as literal exact arithmetic.
  double requested_fuzzy_tolerance = 0.0;
  double effective_fuzzy_tolerance = 0.0;
  double general_fuse_milliseconds = 0.0;
  std::size_t certified_contact_segments_checked = 0;
  std::size_t certified_contact_segments_unreproduced = 0;

  bool valid() const { return errors.empty(); }
};

using ExactCellDecomposer = ExactCellDecomposition (*)(const ExactCellInput &);

// Build the complement of finite, conforming, zero-thickness triangle sheets in one
// axis-aligned box.  OCCT General Fuse is always passed fuzzy tolerance zero; the
// result reports OCCT's unavoidable Precision::Confusion() effective floor.
// A terminating sheet can have plus_fragment == minus_fragment: its two local sides
// reconnect around the free edge and are one global fluid component.
ExactCellDecomposition decompose_exact_cell(const ExactCellInput &input);
} // namespace paracfd::core
