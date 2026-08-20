// occt_conforming_mesh_builder.h -- transactional exact-CAD mesh assembly.
//
// This OpenCascade-facing layer translates an already certified contact
// atomization into face-local UV constraints, conforms every source face, and
// publishes decoded CAD/contact provenance only after the whole model succeeds.
#pragma once

#include "core/geometry/occt_contact_topology_builder.h"
#include "core/geometry/occt_face_triangulation_conformer.h"

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace paracfd::core
{
	// Constraint tags are deliberately outside the 32-bit topology-ID space.
	// The two high-bit namespaces are disjoint even when both payload IDs are zero.
	inline constexpr std::uint64_t occt_cad_edge_constraint_bit =
		std::uint64_t{1} << 62;
	inline constexpr std::uint64_t occt_contact_atom_constraint_bit =
		std::uint64_t{1} << 63;
	inline constexpr std::uint64_t occt_constraint_namespace_mask =
		occt_cad_edge_constraint_bit | occt_contact_atom_constraint_bit;
	inline constexpr std::uint64_t occt_constraint_payload_mask =
		occt_cad_edge_constraint_bit - 1u;

	constexpr std::uint64_t occt_cad_edge_constraint_id(
		std::uint32_t source_edge_id) noexcept
	{
		return occt_cad_edge_constraint_bit | source_edge_id;
	}

	constexpr std::optional<std::uint64_t> occt_contact_atom_constraint_id(
		std::uint64_t atom_id) noexcept
	{
		if (atom_id > occt_constraint_payload_mask) return std::nullopt;
		return occt_contact_atom_constraint_bit | atom_id;
	}

	constexpr std::optional<std::uint32_t> occt_cad_edge_from_constraint(
		std::uint64_t constraint_id) noexcept
	{
		if ((constraint_id & occt_constraint_namespace_mask)
			!= occt_cad_edge_constraint_bit
			|| (constraint_id & occt_constraint_payload_mask)
			> std::numeric_limits<std::uint32_t>::max()) return std::nullopt;
		return static_cast<std::uint32_t>(
			constraint_id & occt_constraint_payload_mask);
	}

	constexpr std::optional<std::uint64_t> occt_contact_atom_from_constraint(
		std::uint64_t constraint_id) noexcept
	{
		if ((constraint_id & occt_constraint_namespace_mask)
			!= occt_contact_atom_constraint_bit) return std::nullopt;
		return constraint_id & occt_constraint_payload_mask;
	}

	struct OcctConformedMeshEdgeProvenance
	{
		std::array<std::uint32_t, 2> vertices{};
		bool boundary = false;
		std::optional<std::uint32_t> source_cad_edge_id;
		// Dense by contract: when present, this directly indexes the input
		// topology.atomization.atoms vector for fan degree and face incidence.
		std::optional<std::uint64_t> contact_atom_id;
	};

	struct OcctConformingFaceMesh
	{
		std::uint32_t source_face_id = std::numeric_limits<std::uint32_t>::max();
		OcctConformedFaceMesh mesh;
		// Indexed exactly like mesh.edges.
		std::vector<OcctConformedMeshEdgeProvenance> edge_provenance;
	};

	struct OcctConformingMesh
	{
		std::vector<OcctConformingFaceMesh> faces;
		// IDs below this count belong to exact contact samples. Ordinary face-local
		// vertices deliberately remain UvVertex::no_topology_id for the importer to
		// allocate globally after this transaction.
		std::uint32_t shared_topology_node_count = 0;
	};

	// Build every source face from the exact contact atomization. Boundary pcurves
	// are evaluated only at atomizer-certified target parameters; no nearest-point
	// recovery or endpoint interpolation is performed here. On any rejection,
	// `output` remains logically unchanged and only `error` is updated.
	bool build_occt_conforming_mesh(const OcctContactTopology& topology,
		OcctConformingMesh& output, std::string& error,
		OcctFaceConformerOptions options = {});
}
