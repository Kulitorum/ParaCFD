// occt_surface_normalization.h -- transactional exact normalization of overlapping faces.
//
// This is an internal OpenCascade-facing API in paracfd_geometry.  Unlike step_import.h,
// it intentionally exposes TopoDS_Face and therefore must only be included by targets that
// already opt into OCCT headers.  STEP import does not call this component yet.
#pragma once

#include <TopoDS_Face.hxx>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace paracfd::core
{
	// One source face whose physical support is represented by a normalized output face.
	// orientation_parity is +1 when its oriented normal agrees with the output carrier and
	// -1 when it opposes it.  It is provenance for side-labelled scalar data; it must never
	// multiply an extensive force or flux assigned once to the physical output face.
	struct OcctSurfaceContributor
	{
		std::size_t source_face_id = 0;
		std::int8_t orientation_parity = 1;
	};

	struct OcctNormalizedSurfaceFace
	{
		TopoDS_Face face;
		std::vector<OcctSurfaceContributor> contributors;
		double area = 0.0; // squared OCCT model units
		std::array<double, 3> centroid{}; // OCCT model units
		double maximum_tolerance = 0.0; // OCCT model units
	};

	struct OcctSurfaceOverlapPairAudit
	{
		std::size_t source_a = 0;
		std::size_t source_b = 0;
		double common_area = 0.0;
		std::array<double, 3> common_centroid{};
		double positive_area_tolerance = 0.0;
		std::vector<std::string> warnings;
	};

	// A positive-area exact Common involving at least one non-planar source face.  The
	// normalizer is deliberately planar-only: these records make unsupported coincident
	// curved coverage explicit and cause the complete transaction to fail closed.
	struct OcctUnsupportedSurfaceOverlap
	{
		std::size_t source_a = 0;
		std::size_t source_b = 0;
		bool source_a_planar = false;
		bool source_b_planar = false;
		double common_area = 0.0;
		std::array<double, 3> common_centroid{};
		double positive_area_tolerance = 0.0;
	};

	struct OcctSurfaceNormalizationClusterAudit
	{
		std::vector<std::size_t> source_face_ids;
		// Filled only for a valid global transaction.  Indices refer to normalized_faces.
		std::vector<std::size_t> normalized_face_indices;
		std::size_t atomic_face_count = 0;
		std::size_t single_source_atom_count = 0;
		std::size_t shared_source_atom_count = 0;
		double input_area_sum = 0.0;
		double once_covered_area = 0.0;
		double duplicate_area_removed = 0.0;
		double input_maximum_tolerance = 0.0;
		double output_maximum_tolerance = 0.0;
		double output_tolerance_limit = 0.0;
		double maximum_rehost_area_error = 0.0;
		double maximum_rehost_first_moment_error = 0.0;
		double maximum_source_area_error = 0.0;
		double maximum_source_first_moment_error = 0.0;
		double maximum_atom_overlap_area = 0.0;
		double arrangement_milliseconds = 0.0;
		std::vector<std::string> warnings;
		std::vector<std::string> errors;

		bool valid() const { return errors.empty(); }
	};

	struct OcctSurfaceNormalizationStats
	{
		std::size_t source_face_count = 0;
		std::size_t planar_source_face_count = 0;
		std::size_t broad_phase_candidate_pair_count = 0;
		std::size_t coplanar_candidate_pair_count = 0;
		// AABB candidates involving at least one non-planar face enter a conservative
		// underlying-support test.  Certified-separated pairs are rejected; only identical
		// handles and equal elementary supports reach Common; generic survivors are unverified.
		std::size_t nonplanar_screened_pair_count = 0;
		std::size_t nonplanar_support_rejected_pair_count = 0;
		// Generic curved supports are intentionally not intersected exactly.  This count
		// exposes candidates for which the planar-only helper makes no overlap claim.
		std::size_t nonplanar_unverified_pair_count = 0;
		std::size_t nonplanar_exact_common_pair_count = 0;
		std::size_t positive_area_overlap_pair_count = 0;
		std::size_t overlap_cluster_count = 0;
		std::size_t untouched_face_count = 0;
		std::size_t normalized_face_count = 0;
		double input_area_sum = 0.0;
		double once_covered_area = 0.0;
		double duplicate_area_removed = 0.0;
		double support_screen_milliseconds = 0.0;
		double detection_milliseconds = 0.0;
		double arrangement_milliseconds = 0.0;
	};

	struct OcctSurfaceNormalizationResult
	{
		// Fail-closed planar transaction: this vector is populated only if every detected
		// planar cluster passes all gates and no cheaply confirmed unsupported non-planar
		// overlap is found.  Generic curved candidates remain observable through the unverified
		// count/warning. Untouched faces are the exact input TopoDS_Face, one-for-one. A cluster
		// replaces its members at the lowest source ID in deterministic geometry/provenance order.
		std::vector<OcctNormalizedSurfaceFace> normalized_faces;
		std::vector<OcctSurfaceOverlapPairAudit> overlap_pairs;
		std::vector<OcctUnsupportedSurfaceOverlap> unsupported_overlaps;
		std::vector<OcctSurfaceNormalizationClusterAudit> clusters;
		OcctSurfaceNormalizationStats stats;
		std::vector<std::string> warnings;
		std::vector<std::string> errors;

		// `valid` means every operation this deliberately planar transaction attempted
		// passed. Importers requiring a complete overlap certificate must additionally
		// require `fully_verified()`: generic curved-support candidates are deliberately
		// bounded and reported rather than subjected to unbounded Boolean work.
		bool valid() const { return errors.empty(); }
		bool fully_verified() const
		{
			return valid() && stats.nonplanar_unverified_pair_count == 0;
		}
	};

	// Normalize exact positive-area overlaps among coplanar faces without sewing, fuzzy
	// tolerance, gap healing, same-domain unification, or tessellation.  Every cluster is
	// accepted only after coverage/provenance, area, first-moment, pairwise-disjointness,
	// BRep validity, coplanarity, and native-tolerance-growth checks pass.  Positive-area
	// coincident overlap on an identical surface handle or equal elementary support is
	// confirmed by exact Common and rejects the transaction.  Other curved candidates are
	// reported as explicitly unverified: this planar-only API makes no completeness claim
	// for generic non-planar coincidence and never spends unbounded time proving its absence.
	OcctSurfaceNormalizationResult normalize_occt_surface_faces(
		const std::vector<TopoDS_Face>& source_faces);
}
