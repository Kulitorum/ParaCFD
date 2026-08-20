// occt_face_triangulation_conformer.h -- exact-surface face-local mesh conforming.
//
// This is an internal OpenCascade-facing preprocessing API. It deliberately exposes
// TopoDS and Poly types and therefore belongs in paracfd_geometry, not in the OCCT-free
// CFD core.
#pragma once

#include "core/geometry/uv_constraint_mesh.h"

#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS_Face.hxx>

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace paracfd::core
{
	struct OcctFaceConstraintSample
	{
		static constexpr std::uint32_t no_source_node_index =
			std::numeric_limits<std::uint32_t>::max();

		// Parameters in this face's native surface chart.
		UvPoint uv;
		// Stable ID shared by every face which owns this physical CAD contact sample.
		// no_topology_id is also valid: it requests a face-local split/tag vertex whose
		// output remains deliberately unclaimed. In that case both canonical fields
		// below are ignored.
		std::uint32_t topology_id = UvVertex::no_topology_id;
		// Canonical position after all TopLoc_Location transforms, in OCCT model units.
		// The exact supplied bits, rather than a second surface evaluation, are emitted
		// for a claimed vertex after the surface round-trip certificate passes.
		std::array<double, 3> canonical_world_position{};
		// Additional source-side native CAD certificate tolerance in world/model units.
		// It is added to this face's independent native bound for this sample only;
		// it never loosens another sample, a source tessellation node, or the UV chart.
		double canonical_position_tolerance = 0.0;
		// Zero-based identity in this face's original Poly_Triangulation when exact
		// CAD preprocessing proved that this contact sample came from that node.
		// This is carried into UV insertion as an identity, never a search radius.
		std::uint32_t source_node_index = no_source_node_index;
		// Exact builder provenance retained for diagnostics. These identities never
		// participate in geometric matching or broaden a certificate tolerance.
		std::uint32_t source_boundary_occurrence_id =
			std::numeric_limits<std::uint32_t>::max();
		std::uint64_t source_contact_atom_id =
			std::numeric_limits<std::uint64_t>::max();
		std::uint32_t source_atom_sample_index =
			std::numeric_limits<std::uint32_t>::max();
	};

	struct OcctFaceConstraintPolyline
	{
		static constexpr std::uint64_t no_chart_branch_id =
			std::numeric_limits<std::uint64_t>::max();

		std::uint64_t constraint_id = 0;
		std::vector<OcctFaceConstraintSample> samples;
		// Explicit identity of one occurrence in this face's UV chart. Most chains
		// need no branch ID. A physical curve lying on a periodic chart cut has two
		// distinct UV occurrences; give those polylines distinct non-sentinel IDs so
		// the conformer can use separate planar insertion vertices while emitting one
		// shared physical topology ID/canonical position. Alias permission is never
		// inferred merely from coincident 3D points.
		std::uint64_t chart_branch_id = no_chart_branch_id;
	};

	struct OcctConformedFaceEdge
	{
		std::array<std::uint32_t, 2> vertices{};
		bool boundary = false;
		std::vector<std::uint64_t> constraint_ids;
	};

	struct OcctConformedFaceMesh
	{
		std::vector<UvPoint> uv;
		std::vector<std::array<double, 3>> world_positions;
		// These triangles are consistently CCW in the face's UV chart. Consumers which
		// require the oriented TopoDS face normal swap indices 1/2 when face_reversed.
		std::vector<std::array<std::uint32_t, 3>> triangles;
		// Ordinary tessellation vertices deliberately retain no_topology_id. Only a
		// prescribed canonical contact sample receives a shared topology ID. One
		// physical ID can occur at more than one output vertex only for caller-declared
		// aliases on distinct branches of a periodic face chart; those vertices still
		// carry the identical canonical world-position bits.
		std::vector<std::uint32_t> topology_ids;
		std::vector<OcctConformedFaceEdge> edges;
		bool face_reversed = false;
		std::size_t source_node_count = 0;
	};

	struct OcctFaceConformerOptions
	{
		// Optional caller-wide source-side certificate tolerance in world/model units.
		// It is added to the independent native face bound. Unlike a per-sample bound,
		// this intentionally applies to source-node and contact round trips alike.
		double geometric_tolerance = 0.0;
	};

	// Transactionally insert prescribed contact polylines into one existing OCCT face
	// triangulation. The current UV triangulation supplies the trimmed domain, including
	// every hole loop. It is normalized to a consistently CCW planar triangulation and
	// passed through UvConstraintMesh, which rejects non-manifold incidence, geometric
	// edge crossings and positive-area overlap before any output is published.
	//
	// Every output position is re-evaluated on the exact located OCCT surface. Claimed
	// contact vertices are then replaced by their canonical world positions only after a
	// native-tolerance round trip. The operation leaves `output` byte-for-byte logically
	// unchanged on failure; only `error` is updated.
	bool conform_occt_face_triangulation(
		const TopoDS_Face& face,
		const occ::handle<Poly_Triangulation>& triangulation,
		const TopLoc_Location& triangulation_location,
		std::span<const OcctFaceConstraintPolyline> constraints,
		OcctConformedFaceMesh& output,
		std::string& error,
		OcctFaceConformerOptions options = {});
}
