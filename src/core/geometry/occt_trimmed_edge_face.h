// occt_trimmed_edge_face.h -- exact source-edge intervals covered by a trimmed face.
//
// This is an internal OpenCascade-facing preprocessing API.  It deliberately exposes
// TopoDS types and therefore belongs in paracfd_geometry, not in the OCCT-free CFD core.
#pragma once

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace paracfd::core
{
	enum class OcctTrimmedEdgeCoverage : std::uint8_t
	{
		none = 0,
		partial = 1,
		full = 2
	};

	enum class OcctFaceIntervalLocation : std::uint8_t
	{
		boundary = 1,
		interior = 2
	};

	struct OcctTrimmedEdgeFaceInterval
	{
		// Parameters are normalized against the source edge's finite underlying-curve
		// range.  They are monotone in that curve parameter and intentionally independent
		// of the TopoDS_Edge FORWARD/REVERSED occurrence orientation.
		double begin = 0.0;
		double end = 0.0;
		OcctFaceIntervalLocation location = OcctFaceIntervalLocation::interior;
		// One sector for a physical trim boundary, two for a face-interior curve.  A
		// genuine periodic seam has two face-side sectors and is classified as interior.
		std::uint8_t target_sector_count = 2;
		// Number of actual wire-edge occurrences covering this interval.  This is zero
		// for an interior curve, one for a physical boundary, and normally two for a seam.
		// Adjacent boundary intervals are intentionally not coalesced across an occurrence
		// transition, even when their counts and sector classifications are identical.
		std::uint8_t boundary_occurrence_count = 0;
		// Zero-based deterministic wire-edge occurrence IDs in the target face's
		// traversal.  Unlike TopoDS_Edge::IsSame, these distinguish the two oriented
		// pcurve branches of a periodic seam. Only the first
		// boundary_occurrence_count entries are active.
		std::array<std::uint32_t, 2> target_boundary_occurrence_ids{{
			std::numeric_limits<std::uint32_t>::max(),
			std::numeric_limits<std::uint32_t>::max() }};
		std::array<std::int8_t, 2> target_boundary_orientations{{ 0, 0 }};
		// Raw parameters on each oriented target edge occurrence corresponding to
		// this source interval. Parameter order follows increasing source parameter;
		// it may therefore be descending on a reversed target occurrence.
		std::array<double, 2> target_parameter_begin{{ 0.0, 0.0 }};
		std::array<double, 2> target_parameter_end{{ 0.0, 0.0 }};
		std::array<double, 2> target_mapping_source_begin{{ 0.0, 0.0 }};
		std::array<double, 2> target_mapping_source_end{{ 0.0, 0.0 }};
		// Additional, local geometric bound certified by the exact CAD operation
		// which produced this covered source subspan.  For BRepAlgoAPI_Common this
		// is the maximum native tolerance of only the result edges overlapping this
		// interval.  It is not a search/fuzzy tolerance and must never be used to
		// create contacts outside this already-certified interval.
		double exact_operation_tolerance = 0.0;
	};

	struct OcctTrimmedEdgeFaceCommonResult
	{
		std::vector<OcctTrimmedEdgeFaceInterval> intervals;
		OcctTrimmedEdgeCoverage coverage = OcctTrimmedEdgeCoverage::none;
		double source_parameter_first = 0.0;
		double source_parameter_last = 0.0;
		double normalized_parameter_tolerance = 0.0;
		std::vector<std::string> errors;

		bool valid() const { return errors.empty(); }
	};

	// Intersect a finite source edge with the *trimmed* target face.  BRepAlgoAPI_Common
	// establishes the actual face coverage first; every resulting edge is then mapped
	// back to the source parameter range by zero-fuzzy IntTools_EdgeEdge.  Boundary versus
	// interior classification comes only from exact overlaps with actual edge occurrences
	// in the target face's wires -- never from a nearest-point or 3-D face classifier.
	//
	// The operation is transactional.  Boolean warnings, Boolean/IntTools failures,
	// incomplete or non-unique mappings, invalid inputs, and ambiguous coincident boundary
	// occurrences all return no intervals and at least one error.
	OcctTrimmedEdgeFaceCommonResult exact_trimmed_edge_face_common(
		const TopoDS_Edge& source_edge, const TopoDS_Face& target_face);
}
