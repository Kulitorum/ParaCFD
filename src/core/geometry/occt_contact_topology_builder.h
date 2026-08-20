// occt_contact_topology_builder.h -- deterministic OCCT topology -> contact atoms.
//
// This is the one place which translates a transferred, already-meshed BRep into
// the exact source-edge/face-occurrence records consumed by the contact atomizer.
// Contact is established only by shared TopoDS topology or by the exact trimmed
// edge/face interval certificates supplied by the caller.  No distance weld or
// triangle-overlap discovery is performed here.
#pragma once

#include "core/geometry/occt_contact_atomizer.h"

#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace paracfd::core
{
	struct OcctContactTopologyEdgeFaceInterval
	{
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		std::uint32_t target_face_id = std::numeric_limits<std::uint32_t>::max();
		double source_parameter_first = 0.0;
		double source_parameter_last = 0.0;
		OcctTrimmedEdgeFaceInterval interval;
	};

	struct OcctContactTopologyBoundaryOccurrence
	{
		// Zero-based traversal position in this face's complete oriented wire-edge
		// occurrence list.  Therefore faces[f].boundary_occurrences[id] is the exact
		// lookup used by OcctTrimmedEdgeFaceInterval.
		std::uint32_t occurrence_id =
			OcctContactBoundaryOccurrence::no_occurrence_id;
		// Exactly the zero-based orientation-insensitive TopExp::MapShapes edge ID.
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		std::uint32_t wire_id = std::numeric_limits<std::uint32_t>::max();
		std::int8_t orientation = 0;
		TopoDS_Edge edge;
	};

	struct OcctContactTopologyFace
	{
		std::uint32_t source_face_id = std::numeric_limits<std::uint32_t>::max();
		TopoDS_Face face;
		std::vector<OcctContactTopologyBoundaryOccurrence> boundary_occurrences;
	};

	struct OcctContactTopologyEdge
	{
		// Exactly the zero-based orientation-insensitive TopExp::MapShapes edge ID.
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		TopoDS_Edge edge;
		std::vector<std::uint32_t> owner_face_ids;
	};

	struct OcctContactTopologyMeshNodeClaim
	{
		// Face-local, zero-based index in the original Poly_Triangulation.  This is
		// provenance, not a nearest-vertex hint: conforming insertion must claim this
		// exact existing node after its native/exact 3-D certificate has been checked.
		std::uint32_t source_face_id = std::numeric_limits<std::uint32_t>::max();
		std::uint32_t source_node_index = std::numeric_limits<std::uint32_t>::max();
		// Boundary polygon claims name their exact oriented occurrence.  Interior
		// edge/face contacts use no_occurrence_id.
		std::uint32_t boundary_occurrence_id =
			OcctContactBoundaryOccurrence::no_occurrence_id;
		std::uint64_t contact_atom_id = std::numeric_limits<std::uint64_t>::max();
		std::uint32_t atom_sample_index = std::numeric_limits<std::uint32_t>::max();
		std::uint32_t topology_id = std::numeric_limits<std::uint32_t>::max();
	};

	struct OcctContactTopology
	{
		std::vector<OcctContactTopologyFace> faces;
		std::vector<OcctContactTopologyEdge> edges;

		// Retained staging records make the transaction auditable and give the face
		// conformer direct access to the exact oriented occurrence identities used by
		// every atom.  They are also the exact input used to produce `atomization`.
		std::vector<OcctContactSourceEdge> atomizer_source_edges;
		std::vector<OcctContactEdgeFaceInterval> atomizer_edge_face_intervals;
		std::vector<OcctContactReciprocalInterval> reciprocal_intervals;
		std::vector<OcctContactExactJunction> exact_junctions;
		OcctContactAtomization atomization;
		std::vector<OcctContactTopologyMeshNodeClaim> mesh_node_claims;
	};

	// Build and atomize a complete exact CAD contact topology. `source_faces` must
	// use the caller's deterministic source-face IDs and every supplied source edge
	// ID must be the zero-based ID from TopExp::MapShapes(shape, TopAbs_EDGE).
	// Existing target-face tessellation nodes on an already-certified contact span
	// are unioned into its source sample chain before atomization; they never create
	// a contact and their canonical positions are reevaluated on the source CAD edge.
	//
	// The operation is transactional: on any malformed occurrence, non-unique
	// mapping, unresolved exact junction or atomizer failure, `output` is unchanged.
	bool build_occt_contact_topology(const TopoDS_Shape& shape,
		std::span<const TopoDS_Face> source_faces,
		std::span<const OcctContactTopologyEdgeFaceInterval> exact_intervals,
		OcctContactTopology& output, std::string& error);
}
