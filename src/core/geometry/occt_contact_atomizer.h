// occt_contact_atomizer.h -- exact interval atomization for CAD fabric contacts.
//
// This API is intentionally OpenCascade-facing.  It consumes contacts which have
// already been certified by exact CAD operations and never discovers contacts by
// proximity.  Keeping it in paracfd_geometry prevents TopoDS types from leaking into
// the OCCT-free CFD core.
#pragma once

#include "core/geometry/occt_trimmed_edge_face.h"

#include <TopoDS_Edge.hxx>

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace paracfd::core
{
	struct OcctContactBoundaryOccurrence
	{
		static constexpr std::uint32_t no_occurrence_id =
			std::numeric_limits<std::uint32_t>::max();

		std::uint32_t occurrence_id = no_occurrence_id;
		// Orientation-insensitive TopExp::MapShapes ID of this target-face edge.
		// This is intentionally distinct from the source edge whose certified
		// edge/face interval is being mapped onto the occurrence.
		std::uint32_t target_source_edge_id =
			std::numeric_limits<std::uint32_t>::max();
		std::int8_t orientation = 0;
		// The exact oriented edge occurrence on the target face.  Atom endpoint
		// parameters are recovered on this curve under native CAD bounds; they are
		// never inferred by linearly interpolating the endpoint parameter pairs.
		TopoDS_Edge target_edge;
		// Parameter correspondence for this exact oriented face-edge occurrence.
		// The source interval may be increasing or decreasing; target parameters are
		// stored in corresponding order and may therefore also be decreasing.
		double source_parameter_begin = 0.0;
		double source_parameter_end = 0.0;
		double target_parameter_begin = 0.0;
		double target_parameter_end = 0.0;
	};

	struct OcctContactOwnerFaceUse
	{
		std::uint32_t source_face_id = 0;
		OcctFaceIntervalLocation location = OcctFaceIntervalLocation::boundary;
		// Per-face sector count, not a count to add for duplicate records.  The
		// atomizer takes the maximum contribution for a face, preventing an exact
		// target certificate and a reciprocal owner copy from double-counting it.
		std::uint8_t sector_count = 1;
		std::vector<OcctContactBoundaryOccurrence> boundary_occurrences;
	};

	struct OcctContactSourceSample
	{
		double parameter = 0.0;
		// Optional source-side native geometric bound in model/world units.  It is
		// used only when de-duplicating this sample on its exact support curve.
		double tolerance = 0.0;
	};

	struct OcctContactSourceEdge
	{
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		TopoDS_Edge edge;
		std::vector<OcctContactOwnerFaceUse> owner_face_uses;
		// Existing face-tessellation samples to retain in every shared contact chain.
		// Atom endpoints, interval transitions and exact junctions are added even if
		// they are absent here.
		std::vector<OcctContactSourceSample> samples;
	};

	struct OcctContactEdgeFaceInterval
	{
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		std::uint32_t target_face_id = 0;
		// The raw finite source range returned beside `interval` by
		// exact_trimmed_edge_face_common().  interval.begin/end are normalized to it.
		double source_parameter_first = 0.0;
		double source_parameter_last = 0.0;
		OcctTrimmedEdgeFaceInterval interval;
		// Actual oriented target-face occurrences indexed exactly like interval's
		// occurrence metadata.  The importer resolves these deterministic IDs on the
		// already-certified target face; the atomizer performs no occurrence search.
		std::array<TopoDS_Edge, 2> target_boundary_occurrences;
		// Global source-edge IDs of the exact target occurrences above.  A source
		// certificate edge and its target boundary copy commonly use opposite or
		// otherwise unrelated parameterizations.
		std::array<std::uint32_t, 2> target_boundary_source_edge_ids{{
			std::numeric_limits<std::uint32_t>::max(),
			std::numeric_limits<std::uint32_t>::max()}};
	};

	struct OcctContactReciprocalInterval
	{
		std::uint32_t source_edge_a = std::numeric_limits<std::uint32_t>::max();
		std::uint32_t source_edge_b = std::numeric_limits<std::uint32_t>::max();
		// Corresponding exact endpoints.  Each pair may be increasing or decreasing,
		// allowing reciprocal TopoDS occurrences and differently parameterized copies.
		double parameter_a_begin = 0.0;
		double parameter_a_end = 0.0;
		double parameter_b_begin = 0.0;
		double parameter_b_end = 0.0;
		// Additional certificate tolerance in model/world units.  This does not make
		// a contact: it is used only to validate/map a caller-supplied exact overlap.
		double geometric_tolerance = 0.0;
	};

	struct OcctContactJunctionIncidence
	{
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		double parameter = 0.0;
	};

	struct OcctContactExactJunction
	{
		// Stable caller identity.  IDs need only be unique inside one atomization.
		std::uint64_t junction_id = std::numeric_limits<std::uint64_t>::max();
		std::vector<OcctContactJunctionIncidence> incidences;
		// Exact canonical world position in OCCT model units.  It is copied bit-for-
		// bit into every output sample at the junction after native-bound validation.
		std::array<double, 3> canonical_world_position{};
		double geometric_tolerance = 0.0;
	};

	struct OcctContactAtomizerInput
	{
		std::span<const OcctContactSourceEdge> source_edges;
		std::span<const OcctContactEdgeFaceInterval> edge_face_intervals;
		std::span<const OcctContactReciprocalInterval> reciprocal_intervals;
		std::span<const OcctContactExactJunction> exact_junctions;
	};

	struct OcctContactAtomSourceSpan
	{
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		// Raw parameters in the atom's deterministic physical direction.  The values
		// may descend even though the underlying finite edge range is increasing.
		double parameter_begin = 0.0;
		double parameter_end = 0.0;
	};

	struct OcctContactAtomBoundaryUse
	{
		// Exact CAD edge used by this target-face boundary occurrence.
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		// Independent edge whose edge/face certificate mapped onto the occurrence.
		// It may use the opposite or an otherwise unrelated parameterization.
		std::uint32_t mapping_source_edge_id =
			std::numeric_limits<std::uint32_t>::max();
		std::uint32_t occurrence_id = OcctContactBoundaryOccurrence::no_occurrence_id;
		std::int8_t orientation = 0;
		// source_parameter_* are on mapping_source_edge_id; target_parameter_* are
		// on source_edge_id. Both pairs are clipped to this atom and follow its
		// deterministic physical direction. Either parameter pair may descend.
		double source_parameter_begin = 0.0;
		double source_parameter_end = 0.0;
		double target_parameter_begin = 0.0;
		double target_parameter_end = 0.0;
		// Exact target-occurrence parameters aligned one-for-one with the enclosing
		// atom's union sample chain.  They are recovered by unique bounded projection,
		// never endpoint-linear interpolation.  Together with the enclosing face ID
		// and occurrence_id, this is sufficient to select the exact oriented pcurve.
		std::vector<double> sample_target_parameters;
	};

	struct OcctContactAtomFaceUse
	{
		std::uint32_t source_face_id = 0;
		OcctFaceIntervalLocation location = OcctFaceIntervalLocation::boundary;
		std::uint8_t sector_count = 1;
		std::vector<OcctContactAtomBoundaryUse> boundary_occurrences;
		// Additional per-face/use bound from the exact operation that certified
		// this atom against this target face. Owner-topology uses remain zero. This
		// bound is deliberately isolated from every other face use of the atom.
		double exact_operation_tolerance = 0.0;
	};

	struct OcctContactAtomSampleSource
	{
		std::uint32_t source_edge_id = std::numeric_limits<std::uint32_t>::max();
		double parameter = 0.0;
	};

	struct OcctContactAtomSample
	{
		std::array<double, 3> canonical_world_position{};
		// Additive source-side native CAD certificate bound in model/world units.
		// The face conformer adds its independent face/occurrence bound per use.
		double canonical_position_tolerance = 0.0;
		std::uint32_t topology_id = std::numeric_limits<std::uint32_t>::max();
		// Every exact source location represented by this common chain node.  These
		// records allow the importer to map the node into each face chart without a
		// second contact search.
		std::vector<OcctContactAtomSampleSource> source_locations;
	};

	struct OcctContactAtom
	{
		// Dense deterministic ID after sorting physical atoms by source-edge/span key.
		std::uint64_t id = 0;
		std::vector<OcctContactAtomSourceSpan> source_spans;
		std::vector<OcctContactAtomFaceUse> face_uses;
		std::vector<OcctContactAtomSample> samples;
		std::uint32_t fan_degree = 0;

		// A fan-one atom is intentionally retained: it is a real fabric boundary/opening,
		// not a missing contact to heal or bridge.
		bool is_open_boundary() const { return fan_degree == 1; }
	};

	struct OcctContactAtomization
	{
		std::vector<OcctContactAtom> atoms;
		std::uint32_t topology_node_count = 0;
	};

	// Transactionally split and canonicalize already-certified exact CAD contacts.
	// No Boolean, distance/proximity contact search, fuzzy weld, or gap closure occurs
	// here.  Reciprocal overlaps and exact endpoint/T junctions must be supplied by the
	// caller.  On failure `output` is unchanged and `error` describes the rejected input.
	bool atomize_occt_contacts(const OcctContactAtomizerInput& input,
		OcctContactAtomization& output, std::string& error);
}
