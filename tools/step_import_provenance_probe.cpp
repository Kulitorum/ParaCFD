#include "core/geometry/model_placement.h"
#include "core/geometry/mesh_clip.h"
#include "core/geometry/step_import.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

using namespace paracfd::core;

namespace
{
	int failures = 0;

	void check(bool condition, const char* label)
	{
		std::printf("%s  %s\n", condition ? "PASS" : "FAIL", label);
		if (!condition) ++failures;
	}

	std::uint64_t edge_key(std::uint32_t a, std::uint32_t b)
	{
		if (a > b) std::swap(a, b);
		return (static_cast<std::uint64_t>(a) << 32) | b;
	}

	struct TopologySummary
	{
		bool structurally_valid = true;
		std::size_t boundary_segments = 0;
		std::size_t labelled_boundary_segments = 0;
		std::size_t unknown_boundary_segments = 0;
		std::size_t labelled_interior_segments = 0;
		std::set<std::pair<std::uint32_t, std::uint32_t>> face_edge_uses;
		std::map<std::uint32_t, std::uint32_t> edge_incidence;
		std::map<std::uint32_t, double> edge_tolerance_m;
		std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> periodic_seam_segments;
	};

	TopologySummary summarize(const TriMesh& mesh)
	{
		TopologySummary result;
		if (!mesh.has_face_provenance() || !mesh.has_cad_edge_provenance())
		{
			result.structurally_valid = false;
			return result;
		}

		using EdgeCounts = std::unordered_map<std::uint64_t, unsigned>;
		std::unordered_map<std::uint32_t, EdgeCounts> counts_by_face;
		for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
		{
			const std::uint32_t face = mesh.source_face_ids[triangle];
			for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
			{
				const std::uint32_t a = mesh.indices[3 * triangle + half_edge];
				const std::uint32_t b = mesh.indices[3 * triangle + (half_edge + 1) % 3];
				if (a >= mesh.vertex_count() || b >= mesh.vertex_count() || a == b)
				{
					result.structurally_valid = false;
					continue;
				}
				++counts_by_face[face][edge_key(a, b)];
			}
		}

		for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
		{
			const std::uint32_t face = mesh.source_face_ids[triangle];
			for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
			{
				const std::uint32_t a = mesh.indices[3 * triangle + half_edge];
				const std::uint32_t b = mesh.indices[3 * triangle + (half_edge + 1) % 3];
				const unsigned occurrences = counts_by_face[face][edge_key(a, b)];
				const std::uint32_t edge = mesh.cad_edge_id(triangle, half_edge);
				const CadEdgeProvenanceState state =
					mesh.cad_edge_provenance_state(triangle, half_edge);
				if (state != CadEdgeProvenanceState::none
					&& state != CadEdgeProvenanceState::known
					&& state != CadEdgeProvenanceState::unknown_boundary)
					result.structurally_valid = false;
				const std::uint32_t incidence = mesh.cad_edge_incident_face_count(triangle, half_edge);
				const double tolerance_m = mesh.cad_edge_tolerance(triangle, half_edge);
				const bool periodic_seam = mesh.cad_edge_is_periodic_seam(triangle, half_edge);
				if (occurrences == 1)
				{
					++result.boundary_segments;
					if (state == CadEdgeProvenanceState::known) ++result.labelled_boundary_segments;
					else if (state == CadEdgeProvenanceState::unknown_boundary)
						++result.unknown_boundary_segments;
					else result.structurally_valid = false;
				}
				else if (edge != TriMesh::kNoCadEdgeId)
				{
					++result.labelled_interior_segments;
				}
				if (state != CadEdgeProvenanceState::known)
				{
					if (edge != TriMesh::kNoCadEdgeId || incidence != 0 || tolerance_m != 0.0 || periodic_seam)
						result.structurally_valid = false;
					continue;
				}
				if (edge == TriMesh::kNoCadEdgeId) result.structurally_valid = false;
				if (incidence == 0 || !std::isfinite(tolerance_m) || tolerance_m < 0.0)
					result.structurally_valid = false;
				result.face_edge_uses.emplace(face, edge);
				auto [found, inserted] = result.edge_incidence.emplace(edge, incidence);
				if (!inserted && found->second != incidence) result.structurally_valid = false;
				auto [tolerance, tolerance_inserted] = result.edge_tolerance_m.emplace(edge, tolerance_m);
				if (!tolerance_inserted && tolerance->second != tolerance_m) result.structurally_valid = false;
				if (periodic_seam)
				{
					if (incidence < 2) result.structurally_valid = false;
					++result.periodic_seam_segments[{face, edge}];
				}
			}
		}
		for (const auto& [face_edge, segments] : result.periodic_seam_segments)
			if (segments < 2) result.structurally_valid = false;
		return result;
	}
}

int main(int argc, char** argv)
{
	const std::string path = argc > 1 ? argv[1] : "Test-Data/NACA2412_C1_S2p03.step";
	bool quality_only = false;
	std::optional<std::size_t> expected_mini_rib_faces;
	for (int argument = 2; argument < argc; ++argument)
	{
		const std::string option = argv[argument];
		if (option == "--quality-only") quality_only = true;
		else if (option == "--expect-mini-rib-faces" && argument + 1 < argc)
			expected_mini_rib_faces = static_cast<std::size_t>(std::stoull(argv[++argument]));
		else
		{
			std::fprintf(stderr, "usage: step_import_provenance_probe [file.step] "
				"[--quality-only] [--expect-mini-rib-faces N]\n");
			return 2;
		}
	}
	std::string error;
	const StepGeometry coarse_geometry = load_step_geometry(path, 2.0, &error);
	const TriMesh& coarse = coarse_geometry.mesh;
	check(!coarse.empty(), "tracked STEP fixture imports");
	if (coarse.empty())
	{
		std::fprintf(stderr, "%s\n", error.c_str());
		return 1;
	}
	check(is_unsupported_mini_rib_representation_name("Mini-rib 12")
		&& is_unsupported_mini_rib_representation_name("  MINI-RIB   7  ")
		&& is_unsupported_mini_rib_representation_name("Mini-rib 12 repaired")
		&& !is_unsupported_mini_rib_representation_name("Rib 12")
		&& !is_unsupported_mini_rib_representation_name("Mini-rib")
		&& !is_unsupported_mini_rib_representation_name("Mini-rib 12foo")
		&& !is_unsupported_mini_rib_representation_name("Not a Mini-rib 12"),
		"mini-rib role matching is narrow, normalized, and deterministic");
	check(coarse_geometry.source_faces.size() == 1 + *std::max_element(
		coarse.source_face_ids.begin(), coarse.source_face_ids.end()),
		"OCC-free source-face metadata parallels deterministic face IDs");
	GeometryQualityReport connected_semantic_policy;
	connected_semantic_policy.exact_connected_component_count = 1;
	GeometryIssue connected_named_surface;
	connected_named_surface.kind = GeometryIssueKind::unsupported_named_surface;
	connected_named_surface.source_triangle_ids = {5, 4, 5};
	connected_semantic_policy.issues.push_back(std::move(connected_named_surface));
	check(default_excluded_triangle_ids(connected_semantic_policy)
		== std::vector<std::uint32_t>({4, 5}),
		"semantic exclusion remains active when the named artifact is topologically connected");
	if (expected_mini_rib_faces)
	{
		std::set<std::uint32_t> semantic_faces, semantic_triangles;
		std::set<std::string> semantic_names;
		for (const GeometryIssue& issue : coarse_geometry.quality.issues)
			if (issue.kind == GeometryIssueKind::unsupported_named_surface)
			{
				semantic_faces.insert(issue.source_face_ids.begin(), issue.source_face_ids.end());
				semantic_triangles.insert(issue.source_triangle_ids.begin(),
					issue.source_triangle_ids.end());
				semantic_names.insert(issue.source_representation_names.begin(),
					issue.source_representation_names.end());
			}
		std::size_t metadata_matches = 0;
		for (const StepSourceFaceMetadata& face : coarse_geometry.source_faces)
			metadata_matches += std::any_of(face.representation_names.begin(),
				face.representation_names.end(), is_unsupported_mini_rib_representation_name);
		check(semantic_faces.size() == *expected_mini_rib_faces
			&& metadata_matches == *expected_mini_rib_faces,
			"STEP representation names classify the expected mini-rib source faces");
		check(!semantic_triangles.empty() && !semantic_names.empty()
			&& std::all_of(semantic_names.begin(), semantic_names.end(),
				is_unsupported_mini_rib_representation_name),
			"semantic diagnostic retains triangle and exact matching-name provenance");
		const std::vector<std::uint32_t> default_excluded =
			default_excluded_triangle_ids(coarse_geometry.quality);
		check(std::includes(default_excluded.begin(), default_excluded.end(),
			semantic_triangles.begin(), semantic_triangles.end()),
			"default exclusion policy contains every semantic mini-rib triangle");
		check(default_excluded.size() == semantic_triangles.size(),
			"overlapping semantic and disconnected-component suggestions are de-duplicated");
	}
	if (quality_only)
	{
		std::printf("STEP quality: components=%zu largest-area=%.12g m^2 warnings=%zu blockers=%zu\n",
			coarse_geometry.quality.exact_connected_component_count,
			coarse_geometry.quality.largest_component_area_m2,
			coarse_geometry.quality.warning_count(), coarse_geometry.quality.blocker_count());
		for (const GeometryIssue& issue : coarse_geometry.quality.issues)
		{
			std::printf("  issue=%llu kind=%u faces=%zu triangles=%zu area=%.12g m^2 boundary-lines=%zu "
				"boundary-length=%.12g m\n", static_cast<unsigned long long>(issue.id),
				static_cast<unsigned>(issue.kind), issue.source_face_ids.size(),
				issue.source_triangle_ids.size(), issue.surface_area_m2,
				issue.boundary_polylines.size(), issue.boundary_length_m);
			if (!issue.source_representation_names.empty())
			{
				std::printf("    representation names:");
				for (const std::string& name : issue.source_representation_names)
					std::printf(" [%s]", name.c_str());
				std::printf("\n");
			}
			std::printf("    faces:");
			for (std::uint32_t face : issue.source_face_ids) std::printf(" %u", face);
			std::printf("\n    triangles:");
			for (std::uint32_t triangle : issue.source_triangle_ids) std::printf(" %u", triangle);
			std::printf("\n");
		}
		return failures == 0 ? 0 : 1;
	}
	check(coarse_geometry.quality.exact_connected_component_count == 1
		&& coarse_geometry.quality.issues.empty()
		&& coarse_geometry.quality.flow_eligibility() == GeometryFlowEligibility::clean,
		"closed STEP fixture is one exact CAD-connected fabric component");
	GeometryQualityReport warning_contract;
	GeometryIssue warning;
	warning.kind = GeometryIssueKind::disconnected_fabric_component;
	warning.severity = GeometryIssueSeverity::warning_run;
	warning_contract.issues.push_back(warning);
	check(warning_contract.flow_eligibility() == GeometryFlowEligibility::ready_with_warnings
		&& warning_contract.warning_count() == 1 && warning_contract.blocker_count() == 0,
		"representable geometry warnings retain flow eligibility");

	check(coarse.has_fp64_positions() && coarse.positions_fp64.size() == coarse.positions.size(),
		"STEP import retains one FP64 coordinate beside every FP32 coordinate");
	bool rounded_copy = true;
	for (std::size_t coordinate = 0; coordinate < coarse.positions.size(); ++coordinate)
		rounded_copy = rounded_copy
			&& coarse.positions[coordinate] == static_cast<float>(coarse.positions_fp64[coordinate]);
	check(rounded_copy, "FP32 display positions are direct rounded copies of OCCT FP64 positions");

	const TopologySummary coarse_topology = summarize(coarse);
	check(coarse_topology.structurally_valid, "CAD edge IDs and incident-face counts are structurally consistent");
	check(coarse_topology.boundary_segments > 0
		&& coarse_topology.labelled_boundary_segments == coarse_topology.boundary_segments,
		"every face-local tessellation boundary segment maps to a CAD edge");
	check(coarse_topology.labelled_interior_segments == 0,
		"internal tessellation half-edges retain the unavailable sentinel");
	check(std::all_of(coarse_topology.periodic_seam_segments.begin(),
		coarse_topology.periodic_seam_segments.end(),
		[](const auto& entry) { return entry.second >= 2; }),
		"both periodic-seam polygon chains carry one CAD edge ID with incidence at least two");
	check(std::any_of(coarse_topology.edge_incidence.begin(), coarse_topology.edge_incidence.end(),
		[](const auto& entry) { return entry.second == 2; }),
		"closed fixture exposes CAD edges incident on two source faces");

	// The importer deliberately duplicates tessellation nodes per source face. Confirm that
	// a shared geometric CAD boundary has separate mesh vertices, while its CAD edge ID joins
	// the provenance across those copies.
	using QuantizedPoint = std::tuple<long long, long long, long long>;
	std::map<QuantizedPoint, std::set<std::uint32_t>> coordinate_faces;
	for (std::size_t triangle = 0; triangle < coarse.triangle_count(); ++triangle)
		for (unsigned corner = 0; corner < 3; ++corner)
		{
			const std::uint32_t vertex = coarse.indices[3 * triangle + corner];
			const auto p = coarse.vertex_position_double(vertex);
			coordinate_faces[{std::llround(p[0] * 1e8), std::llround(p[1] * 1e8),
				std::llround(p[2] * 1e8)}].insert(coarse.source_face_ids[triangle]);
		}
	check(std::any_of(coordinate_faces.begin(), coordinate_faces.end(),
		[](const auto& entry) { return entry.second.size() > 1; }),
		"face-local duplicate vertices retain shared CAD topology through edge IDs");

	const TriMesh fine = load_step_mesh(path, 0.5, &error);
	const TopologySummary fine_topology = summarize(fine);
	check(!fine.empty() && fine_topology.structurally_valid,
		"same STEP topology imports at a second tessellation tolerance");
	check(fine_topology.face_edge_uses == coarse_topology.face_edge_uses
		&& fine_topology.edge_incidence == coarse_topology.edge_incidence
		&& fine_topology.edge_tolerance_m == coarse_topology.edge_tolerance_m
		&& fine_topology.periodic_seam_segments.size() == coarse_topology.periodic_seam_segments.size()
		&& std::equal(fine_topology.periodic_seam_segments.begin(),
			fine_topology.periodic_seam_segments.end(), coarse_topology.periodic_seam_segments.begin(),
			[](const auto& a, const auto& b) { return a.first == b.first; }),
		"source face/edge IDs and incidence remain stable across re-tessellation");

	TriMesh programmatic;
	programmatic.positions = { 1.25f, -2.5f, 3.75f };
	const auto promoted = programmatic.vertex_position_double(0);
	check(!programmatic.has_fp64_positions() && promoted[0] == 1.25 && promoted[1] == -2.5
		&& promoted[2] == 3.75,
		"programmatic mesh accessor promotes FP32 when the exact sidecar is absent");

	TriMesh exact_programmatic = programmatic;
	exact_programmatic.positions_fp64 = { 1.25 + 1e-11, -2.5 - 2e-11, 3.75 + 3e-11 };
	ModelPlacement placement;
	placement.tx = 10.0;
	placement.ty = -20.0;
	placement.tz = 30.0;
	const TriMesh placed = placed_mesh(exact_programmatic, placement);
	check(placed.has_fp64_positions()
		&& std::abs(placed.positions_fp64[0] - (11.25 + 1e-11)) < 1e-14
		&& std::abs(placed.positions_fp64[1] - (-22.5 - 2e-11)) < 1e-14
		&& std::abs(placed.positions_fp64[2] - (33.75 + 3e-11)) < 1e-14,
		"model placement transforms retained CPU geometry without an FP32 round trip");

	TriMesh tolerance_mesh;
	tolerance_mesh.positions = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	tolerance_mesh.positions_fp64 = { 0, 0, 1e-11, 1, 0, 1e-11, 0, 1, 1e-11 };
	tolerance_mesh.indices = { 0, 1, 2 };
	tolerance_mesh.triangle_cad_edge_provenance_states = {
		static_cast<std::uint8_t>(CadEdgeProvenanceState::known),
		static_cast<std::uint8_t>(CadEdgeProvenanceState::none),
		static_cast<std::uint8_t>(CadEdgeProvenanceState::known) };
	tolerance_mesh.triangle_cad_edge_ids = { 7, TriMesh::kNoCadEdgeId, 8 };
	tolerance_mesh.triangle_cad_edge_incident_face_counts = { 2, 0, 1 };
	tolerance_mesh.triangle_cad_edge_tolerances = { 1e-6, 0.0, 2e-6 };
	tolerance_mesh.triangle_cad_edge_is_periodic_seam = { 0, 0, 0 };
	ModelPlacement anisotropic;
	anisotropic.m[0] = 2.0;
	anisotropic.m[4] = 3.0;
	anisotropic.m[8] = 0.5;
	const TriMesh scaled = placed_mesh(tolerance_mesh, anisotropic);
	check(scaled.cad_edge_provenance_state(0, 0) == CadEdgeProvenanceState::known
		&& scaled.cad_edge_provenance_state(0, 1) == CadEdgeProvenanceState::none
		&& scaled.cad_edge_tolerance(0, 0) >= 3e-6
		&& scaled.cad_edge_tolerance(0, 0) < 3e-6 * (1.0 + 1e-12)
		&& scaled.cad_edge_tolerance(0, 1) == 0.0
		&& scaled.cad_edge_tolerance(0, 2) >= 6e-6,
		"affine placement scales CAD tolerances by its largest singular value");

	const TriMesh clipped = clip_mesh_to_axis_slab(tolerance_mesh, 0, 0.25, 0.75);
	bool retained_source_edge = false, introduced_sentinel = false;
	for (std::size_t triangle = 0; triangle < clipped.triangle_count(); ++triangle)
		for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
		{
			const std::uint32_t id = clipped.cad_edge_id(triangle, half_edge);
			if (id == 7)
				retained_source_edge = retained_source_edge
					|| (clipped.cad_edge_incident_face_count(triangle, half_edge) == 2
						&& clipped.cad_edge_tolerance(triangle, half_edge) == 1e-6);
			else if (id == TriMesh::kNoCadEdgeId)
				introduced_sentinel = introduced_sentinel
					|| (clipped.cad_edge_incident_face_count(triangle, half_edge) == 0
						&& clipped.cad_edge_tolerance(triangle, half_edge) == 0.0);
		}
	check(clipped.has_fp64_positions() && clipped.has_cad_edge_provenance()
		&& retained_source_edge && introduced_sentinel,
		"slab clipping retains exact/source-edge metadata and marks synthetic edges unavailable");

	TriMesh unknown_boundary_mesh = tolerance_mesh;
	unknown_boundary_mesh.triangle_cad_edge_provenance_states[0] =
		static_cast<std::uint8_t>(CadEdgeProvenanceState::unknown_boundary);
	unknown_boundary_mesh.triangle_cad_edge_ids[0] = TriMesh::kNoCadEdgeId;
	unknown_boundary_mesh.triangle_cad_edge_incident_face_counts[0] = 0;
	unknown_boundary_mesh.triangle_cad_edge_tolerances[0] = 0.0;
	const TriMesh clipped_unknown = clip_mesh_to_axis_slab(unknown_boundary_mesh, 0, 0.25, 0.75);
	bool retained_unknown = false;
	for (std::size_t triangle = 0; triangle < clipped_unknown.triangle_count(); ++triangle)
		for (unsigned half_edge = 0; half_edge < 3; ++half_edge)
			retained_unknown = retained_unknown || clipped_unknown.cad_edge_provenance_state(
				triangle, half_edge) == CadEdgeProvenanceState::unknown_boundary;
	check(retained_unknown,
		"slab clipping preserves an inherited unresolved CAD boundary distinctly from synthetic edges");

	std::printf("STEP provenance: vertices=%zu triangles=%zu CAD edges=%zu face-edge uses=%zu\n",
		coarse.vertex_count(), coarse.triangle_count(), coarse_topology.edge_incidence.size(),
		coarse_topology.face_edge_uses.size());
	return failures == 0 ? 0 : 1;
}
